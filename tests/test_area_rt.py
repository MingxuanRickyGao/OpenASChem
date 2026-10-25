import ast
import importlib.util
import unittest
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pandas as pd

CODE = Path(__file__).resolve().parents[1] / "src" / "openaschem"
spec = importlib.util.spec_from_file_location(
    "area_test_chromatogram_processing", CODE / "chromatogram_processing.py"
)
chrom = importlib.util.module_from_spec(spec)
spec.loader.exec_module(chrom)


class AreaRTTests(unittest.TestCase):
    def test_nonuniform_rt_quadratic_has_analytic_area(self):
        rt = np.array([0, 0.1, 0.3, 0.8, 1.1, 1.5, 2.0])
        signal = np.array([0, *rt[1:6] ** 2, 0])
        actual = chrom.quantify(signal, signal, 1, 5, "area", rt_list=rt)
        self.assertAlmostEqual(actual, (1.5**3 - 0.1**3) / 3)

    def test_same_physical_peak_at_different_sampling_densities(self):
        areas = []
        for count in [5, 9, 17]:
            peak_rt = np.linspace(-1, 1, count)
            rt = np.r_[-2, peak_rt, 2]
            signal = np.r_[0, 1 - peak_rt**2, 0]
            areas.append(chrom.quantify(signal, signal, 1, count, rt_list=rt))
        np.testing.assert_allclose(areas, [4 / 3] * 3, rtol=1e-12)

    def test_area_scales_with_elapsed_time_and_not_time_origin(self):
        signal = np.array([0.0, 0.0, 10.0, 20.0, 10.0, 0.0, 0.0])
        rt = np.arange(7) * 0.1
        area = chrom.quantify(signal, signal, 1, 5, rt_list=rt)
        self.assertAlmostEqual(chrom.quantify(signal, signal, 1, 5, rt_list=rt * 5), area * 5)
        self.assertAlmostEqual(chrom.quantify(signal, signal, 1, 5, rt_list=rt + 100), area)

    def test_baseline_subtraction_preserved(self):
        signal = np.array([7.0, 7.0, 17.0, 27.0, 17.0, 7.0, 7.0])
        corrected = signal - 7
        rt = np.arange(7) * 0.4
        self.assertAlmostEqual(
            chrom.quantify(signal, signal, 1, 5, rt_list=rt),
            chrom.quantify(corrected, corrected, 1, 5, rt_list=rt),
        )

    def test_intensity_mode_keeps_original_positional_api(self):
        smooth = np.array([0.0, 5.0, 10.0, 5.0, 0.0])
        original = np.array([0.0, 1.0, 20.0, 1.0, 0.0])
        self.assertEqual(chrom.quantify(smooth, original, 1, 3, "intensity"), 20)

    def test_area_without_rt_is_rejected(self):
        with self.assertRaisesRegex(ValueError, "rt_list is required"):
            chrom.quantify(np.zeros(5), np.zeros(5), 1, 3)

    def test_invalid_rt_is_rejected(self):
        for rt in [
            [0, 1],
            [[0, 1, 2, 3, 4]],
            [0, 1, 1, 3, 4],
            [0, 2, 1, 3, 4],
            [0, 1, np.nan, 3, 4],
        ]:
            with self.subTest(rt=rt), self.assertRaises(ValueError):
                chrom.quantify(np.zeros(5), np.zeros(5), 1, 3, rt_list=rt)

    def test_auxiliary_peak_function_uses_rt_area(self):
        signal = np.array([0.0, 0.0, 10.0, 20.0, 10.0, 0.0, 0.0])
        rt = np.arange(7) * 0.3
        left, right, area = chrom.find_best_peak_intensity_with_bg_subtraction(signal, rt)
        self.assertAlmostEqual(area, chrom.quantify(signal, signal, left, right, rt_list=rt))

    def test_pipeline_uses_each_original_run_rt(self):
        # Execute the production peak-identification function with identity standard
        # alignment, so this regression isolates the four quantification call sites.
        tree = ast.parse((CODE / "open_identify_utils.py").read_text(encoding="utf-8"))
        function = next(
            n for n in tree.body if isinstance(n, ast.FunctionDef) and n.name == "identify_peaks"
        )
        signal = np.array([0.0, 0.0, 10.0, 20.0, 10.0, 0.0, 0.0])
        identity = pd.DataFrame({name: np.arange(7) for name in ["S1", "S2", "S3"]})
        ns = dict(
            np=np,
            pd=pd,
            find_top_peaks=chrom.find_top_peaks,
            quantify=chrom.quantify,
            get_openaschem_logger=lambda: SimpleNamespace(info=lambda *a: None),
            align_one_compound_standards=lambda *a: identity,
            time_vector_to_index_vector=lambda *a: identity,
            interpolate_index_vector=lambda v: v,
            index_vector_to_score_vector=lambda index, values: pd.DataFrame(values),
        )
        exec(
            compile(
                ast.Module(body=[function], type_ignores=[]),
                str(CODE / "open_identify_utils.py"),
                "exec",
            ),
            ns,
        )
        standards = [
            SimpleNamespace(rt_list=(100 + np.arange(7) * dt).tolist()) for dt in [0.5, 0.6, 0.7]
        ]
        experimental = SimpleNamespace(rt_list=(100 + np.arange(7) * 0.8).tolist())
        standard_xics = [{"compound": signal.copy()} for _ in range(3)]
        std, exp = ns[function.name](
            pd.DataFrame({"id": ["compound"]}),
            {"COMPOUND_ID": "id"},
            standards,
            standard_xics,
            standard_xics,
            {"run": experimental},
            {"run": {"compound": signal}},
            {"run": {"compound": signal}},
            ["run"],
            1,
            0,
            5,
            *(
                [None]
                if any(a.arg == "getAlignedTimesFast_function" for a in function.args.args)
                else []
            ),
            None,
            None,
            1,
            3,
            100,
            0,
            0,
            999,
            10,
            5,
            "area",
        )
        reference = std["compound"][0]
        for i in range(3):
            self.assertAlmostEqual(
                reference[f"q{i + 1}"],
                chrom.quantify(signal, signal, 2, 4, rt_list=standards[i].rt_list),
            )
        self.assertAlmostEqual(
            exp["compound"][0][2]["run"],
            chrom.quantify(signal, signal, 2, 4, rt_list=experimental.rt_list),
        )


if __name__ == "__main__":
    unittest.main(verbosity=2)
