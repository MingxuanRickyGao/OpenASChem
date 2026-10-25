import base64
import sys
import tempfile
import unittest
import zlib
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from openaschem.interactive_report import build_report_data, encode_trace, write_interactive_report


def fixture():
    cid = "CMP_with_underscore"
    key = cid + "_1"
    prefix = "STD_pool"
    runs = [f"{sample}_pool_{r}" for sample in ["T1", "T10"] for r in [1, 2, 3]]
    axis = [10.0, 11.0, 12.0, 13.0, 14.0]
    trace = np.array([0, 1, 4, 1, 0.0])
    cols = ["mz", "SGC ID for Component", "formula", "SMILES", "ProtID"]
    output_columns = {
        "COMPOUND_ID": cols[1],
        "COMPOUND_MZ": cols[0],
        "FORMULA": cols[2],
        "SMILES": cols[3],
        "TARGET_PROTEIN": cols[4],
    }
    matrix = pd.DataFrame(
        {name: [1.234567890123] for name in [f"{prefix}_{i}" for i in [1, 2, 3]] + runs},
        index=[key],
    )
    standard = {
        cid: [
            {
                **{f"lb{i}": 11.0 for i in [1, 2, 3]},
                **{f"rb{i}": 13.0 for i in [1, 2, 3]},
                **{f"q{i}": 6.0 for i in [1, 2, 3]},
            }
        ]
    }
    experimental = {
        cid: [
            [{run: 11.0 for run in runs}, {run: 13.0 for run in runs}, {run: 6.0 for run in runs}]
        ]
    }
    return dict(
        pool="pool",
        quantification_method="area",
        standard_prefix=prefix,
        quantification_matrix=matrix,
        fold_changes=pd.DataFrame({"T1_pool": [10.0]}, index=[key]),
        results_compound=pd.DataFrame(
            {"ProtID": ["T1_pool"], "fold_change": ["10.0"], "p_value": [0.001]}
        ),
        results_all=pd.DataFrame(
            {
                "SGC ID for Component": [cid],
                "mz": [123.456789],
                "formula": ["C1H2"],
                "SMILES": ["CC"],
            }
        ),
        output_columns=output_columns,
        spectra_standard=[SimpleNamespace(rt_list=axis) for _ in [1, 2, 3]],
        XICs_standard=[{cid: trace} for _ in [1, 2, 3]],
        orig_XICs_standard=[{cid: trace * 2} for _ in [1, 2, 3]],
        spectra_runs={run: SimpleNamespace(rt_list=axis) for run in runs},
        XICs_runs={run: {cid: trace} for run in runs},
        orig_XICs_runs={run: {cid: trace * 2} for run in runs},
        standard_peak_boundaries_and_quantities=standard,
        runs_peak_boundaries_and_quantities=experimental,
    )


class InteractiveReportTests(unittest.TestCase):
    def test_trace_codec_preserves_float32_values(self):
        values = np.array([0, 1.23456789, 1e8], dtype=float)
        decoded = np.frombuffer(
            zlib.decompress(base64.b64decode(encode_trace(values))), dtype="<f4"
        )
        np.testing.assert_array_equal(decoded, values.astype("<f4"))

    def test_bounds_quantities_and_exact_sample_identities_are_preserved(self):
        args = fixture()
        data = build_report_data(**args)
        r = data["records"][0]
        self.assertEqual(r["compound"], "CMP_with_underscore")
        self.assertEqual(r["targets"], ["T1_pool"])
        self.assertIn("T10_pool", [run["group"] for run in data["runs"]])
        self.assertEqual(r["bounds"], [[11.0, 13.0]] * 9)
        self.assertEqual(r["quantities"][0], 1.234567890123)
        self.assertEqual(r["peak_quantities"][0], 6.0)
        self.assertEqual(data["runs"][0]["rt"], [10.0, 11.0, 12.0, 13.0, 14.0])

    def test_empty_hit_list_keeps_quantified_records(self):
        args = fixture()
        args["fold_changes"] = args["fold_changes"].iloc[:0]
        args["results_compound"] = args["results_compound"].iloc[:0]
        data = build_report_data(**args)
        self.assertEqual(len(data["records"]), 1)
        self.assertFalse(data["records"][0]["hit"])

    def test_empty_quantification_is_reportable(self):
        args = fixture()
        args["quantification_matrix"] = args["quantification_matrix"].iloc[:0]
        args["fold_changes"] = args["fold_changes"].iloc[:0]
        args["results_compound"] = args["results_compound"].iloc[:0]
        data = build_report_data(**args)
        self.assertEqual(data["records"], [])

    def test_invalid_trace_rt_length_is_rejected(self):
        args = fixture()
        args["XICs_standard"][0]["CMP_with_underscore"] = np.zeros(3)
        with self.assertRaisesRegex(ValueError, "length mismatch"):
            build_report_data(**args)

    def test_report_escapes_embedded_metadata_and_embeds_plotly(self):
        args = fixture()
        args["results_all"].loc[0, "SMILES"] = "</script><img src=x onerror=alert(1)>"
        with tempfile.TemporaryDirectory() as out:
            info = write_interactive_report(out=out, **args)
            html = Path(info["path"]).read_text(encoding="utf-8")
            self.assertNotIn("</script><img src=x", html)
            self.assertIn("\\u003c/script>", html)
            self.assertNotIn("/*__DATA__*/", html)
            self.assertNotIn("/*__PLOTLY__*/", html)
            self.assertEqual(info["hits"], 1)


if __name__ == "__main__":
    unittest.main()
