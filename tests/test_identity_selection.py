import copy
import sys
import unittest
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from test_interactive_report import fixture

import openaschem.open_identify_utils as u
from openaschem.chromatogram_processing import extract_xic
from openaschem.interactive_report import build_report_data
from openaschem.peak_identity import normalize_parent_library
from openaschem.peak_identity_evidence import vector_extract
from openaschem.peak_identity_selection import select_candidates


def candidates():
    return pd.DataFrame(
        [
            dict(
                key="CMP_1",
                compound="CMP",
                rt_s=100.0,
                height=100.0,
                relative_height=1.0,
                mass_error_ppm=2.0,
                isotope_loss=0.1,
            ),
            dict(
                key="CMP_12",
                compound="CMP",
                rt_s=130.0,
                height=90.0,
                relative_height=0.9,
                mass_error_ppm=-0.5,
                isotope_loss=0.2,
            ),
        ]
    ).set_index("key", drop=False)


class IdentitySelectionTests(unittest.TestCase):
    def test_mass_supported_alternative_can_beat_highest_peak(self):
        selected, _, _, _ = select_candidates(candidates(), mass_center=-0.5)
        self.assertEqual(selected.index.tolist(), ["CMP_12"])

    def test_selection_never_reads_pool_or_fc(self):
        f = candidates()
        first = select_candidates(f, mass_center=-0.5)[0]
        f["pool"] = ["A", "B"]
        f["fc"] = [100.0, 0.1]
        second = select_candidates(f, mass_center=-0.5)[0]
        self.assertEqual(first.index.tolist(), second.index.tolist())

    def test_competitive_distinct_peaks_are_deferred(self):
        f = candidates()
        f["mass_error_ppm"] = -0.5
        f["isotope_loss"] = 0.1
        selected, ambiguous, _, _ = select_candidates(f, mass_center=-0.5)
        self.assertTrue(selected.empty)
        self.assertEqual(ambiguous.compound.tolist(), ["CMP"])

    def test_nearby_local_maxima_do_not_create_distinct_rt_ambiguity(self):
        f = candidates()
        f["mass_error_ppm"] = -0.5
        f["isotope_loss"] = 0.1
        f["rt_s"] = [100.0, 102.0]
        self.assertEqual(len(select_candidates(f, mass_center=-0.5)[0]), 1)

    def test_too_weak_alternative_cannot_replace_reference(self):
        f = candidates()
        f.loc["CMP_12", "relative_height"] = 0.4
        self.assertEqual(select_candidates(f, mass_center=-0.5)[0].index.tolist(), ["CMP_1"])

    def test_empty_candidate_set_is_reportable(self):
        f = candidates().iloc[:0]
        selected, ambiguous, _, center = select_candidates(f)
        self.assertTrue(selected.empty and ambiguous.empty)
        self.assertIsNone(center)

    def test_invalid_settings_rejected(self):
        for kwargs in [dict(floor=0), dict(floor=np.nan), dict(floor=1.1), dict(min_margin=-1)]:
            with self.assertRaises(ValueError):
                select_candidates(candidates(), **kwargs)

    def test_vector_extraction_matches_production_ppm_and_da(self):
        spectra = SimpleNamespace(
            spectra=[
                (np.array([99.9998, 100.0, 100.0002, 101.0]), np.array([2.0, 4.0, 8.0, 16.0])),
                (np.array([100.0, 101.0]), np.array([3.0, 5.0])),
            ]
        )
        for unit, tol in [("ppm", 3.0), ("Da", 0.0003)]:
            values, _ = vector_extract(spectra, np.array([100.0, 101.0]), unit, tol)
            for j, target in enumerate([100.0, 101.0]):
                np.testing.assert_array_equal(
                    values[j], extract_xic(spectra.spectra, target, unit, tol)
                )

    def test_parent_mass_normalization_uses_formula_or_explicit_neutral_mass(self):
        p = normalize_parent_library(
            pd.DataFrame(
                {"SGC ID for Component": ["P"], "formula": ["C2H6"], "Comp Acc Mass": [999.0]}
            )
        )
        self.assertAlmostEqual(p.MONOISOTOPIC_MASS.iloc[0], 30.04695019, places=6)
        self.assertEqual(p["SGC ID for Compotent"].iloc[0], "P")

    def test_multi_digit_candidate_and_underscored_id_export(self):
        cid = "CMP_with_underscore"
        key = cid + "_12"
        prefix = "STD_pool"
        columns = {
            "COMPOUND_ID": "id",
            "COMPOUND_MZ": "mass",
            "POOL_NAME": "pool",
            "SMILES": "smiles",
            "FORMULA": "formula",
            "POOL_ID": "pool_id",
        }
        output = {
            "COMPOUND_ID": "id",
            "COMPOUND_MZ": "mz",
            "POOL_NAME": "ChemicalID",
            "SMILES": "smiles",
            "FORMULA": "formula",
            "POOL_ID": "pool_id",
            "RT": "rt",
            "TARGET_PROTEIN": "ProtID",
            "ADDUCT_TYPE": "adduct",
        }
        library = pd.DataFrame(
            {
                "id": [cid],
                "mass": [100.0],
                "pool": ["pool"],
                "smiles": ["CC"],
                "formula": ["C2H6"],
                "pool_id": ["P"],
            }
        )
        matrix = pd.DataFrame(
            {
                **{f"{prefix}_{i}": [10.0] for i in [1, 2, 3]},
                **{f"T_pool_{i}": [50.0] for i in [1, 2, 3]},
            },
            index=[key],
        )
        standard = {
            cid: [
                {
                    **{f"lb{i}": float(j) for i in [1, 2, 3]},
                    **{f"rb{i}": float(j + 2) for i in [1, 2, 3]},
                    **{f"q{i}": 10.0 for i in [1, 2, 3]},
                }
                for j in range(12)
            ]
        }
        exported = u.build_quantification_table(
            library, columns, [], output, matrix, standard, "positive", prefix, False
        )
        self.assertEqual(exported["id"].iloc[0], cid)
        self.assertEqual(exported.ReferencePeakID.iloc[0], key)
        self.assertEqual(exported.rt.iloc[0], 12.0)
        fc = pd.DataFrame({"T_pool": [5.00001]}, index=[key])
        result = u.build_hit_table(
            library,
            columns,
            [],
            output,
            fc,
            [0.001],
            matrix,
            standard,
            {},
            "positive",
            5.0,
            prefix,
            False,
        )
        self.assertEqual(result.ReferencePeakID.iloc[0], key)
        self.assertGreater(float(result.fold_change.iloc[0]), 5.0)

    def test_empty_matrix_can_be_exported(self):
        matrix = u.build_quantification_matrix(
            {}, {}, "STD_pool", 5000.0, 13.0, 10.0, experimental_runs=["T_pool_1"]
        )
        self.assertEqual(
            matrix.columns.tolist(), ["STD_pool_1", "STD_pool_2", "STD_pool_3", "T_pool_1"]
        )
        self.assertTrue(matrix.empty)

    def test_report_hit_mapping_uses_reference_ids_after_row_sorting(self):
        args = fixture()
        cid = "CMP_with_underscore"
        other = "OTHER"
        key = cid + "_1"
        otherkey = other + "_1"
        for name in [
            "standard_peak_boundaries_and_quantities",
            "runs_peak_boundaries_and_quantities",
        ]:
            args[name][other] = copy.deepcopy(args[name][cid])
        for name in ["XICs_standard", "orig_XICs_standard"]:
            for item in args[name]:
                item[other] = item[cid].copy()
        for name in ["XICs_runs", "orig_XICs_runs"]:
            for item in args[name].values():
                item[other] = item[cid].copy()
        args["quantification_matrix"] = pd.concat(
            [
                args["quantification_matrix"],
                args["quantification_matrix"].rename(index={key: otherkey}),
            ]
        )
        args["fold_changes"] = pd.DataFrame(
            {"T1_pool": [10.0, 1.0], "T10_pool": [1.0, 20.0]}, index=[key, otherkey]
        )
        args["results_compound"] = pd.DataFrame(
            {
                "ReferencePeakID": [otherkey, key],
                "ProtID": ["T10_pool", "T1_pool"],
                "fold_change": ["20", "10"],
                "p_value": [0.002, 0.001],
            }
        )
        data = build_report_data(**args)
        records = {r["key"]: r for r in data["records"]}
        self.assertEqual(records[key]["targets"], ["T1_pool"])
        self.assertEqual(records[otherkey]["targets"], ["T10_pool"])


if __name__ == "__main__":
    unittest.main()
