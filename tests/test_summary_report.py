import csv
import json
import sys
import tempfile
import unittest
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from openaschem.summary_report import (
    build_experiment_summary,
    build_pool_summary,
    relative_link,
    write_experiment_summary,
)

COLS = {
    "COMPOUND_ID": "id",
    "TARGET_PROTEIN": "target",
    "COMPOUND_MZ": "mz",
    "RT": "rt",
    "FORMULA": "formula",
}


def entry(pool="P1", compound="CMP", targets=None, fc="10", **kwargs):
    frame = pd.DataFrame(
        [
            dict(
                id=compound,
                target=targets or "T1_" + pool,
                fold_change=fc,
                p_value=0.000123,
                mz=123.456,
                rt=2.0,
                formula="C6H12",
                ReferencePeakID=compound + "_12",
            )
        ]
    )
    return build_pool_summary(
        pool=pool,
        results_compound=frame,
        output_columns=COLS,
        input_compounds=20,
        run_count=9,
        quantified_compounds=15,
        selection_summary=dict(mode="identity", ambiguous_compounds=2, isotope_veto_candidates=3),
        **kwargs,
    )


class SummaryReportTests(unittest.TestCase):
    def test_multiple_targets_keep_corresponding_fc_and_full_names(self):
        e = entry(targets="T1_P1;T10_P1", fc="10;20")
        self.assertEqual(e["hit_compounds"], 1)
        self.assertEqual(e["hit_target_rows"], 2)
        self.assertEqual(
            [(h["target_label"], h["fold_change"]) for h in e["hits"]],
            [("T1", 10.0), ("T10", 20.0)],
        )
        self.assertEqual(e["hits"][1]["target"], "T10_P1")
        self.assertEqual(e["hits"][0]["rt_s"], 120.0)

    def test_duplicate_compound_across_pools_is_not_collapsed(self):
        data = build_experiment_summary([entry("P2"), entry("P1")])
        self.assertEqual(data["counts"]["pool_compound_hits"], 2)
        self.assertEqual(data["counts"]["unique_hit_compounds"], 1)
        self.assertEqual(len(data["hits"]), 2)
        self.assertEqual([r["pool"] for r in data["pools"]], ["P1", "P2"])

    def test_only_exact_pool_suffix_is_removed(self):
        e = entry(pool="P_with_underscores", targets="T1_P_with_underscores;T10_other", fc="3;4")
        self.assertEqual([h["target_label"] for h in e["hits"]], ["T1", "T10_other"])

    def test_mismatched_target_fc_annotations_are_rejected(self):
        with self.assertRaises(ValueError):
            entry(targets="T1_P1;T10_P1", fc="10")

    def test_empty_pool_still_appears_in_summary(self):
        e = build_pool_summary(
            pool="P1",
            results_compound=pd.DataFrame(),
            output_columns=COLS,
            input_compounds=20,
            run_count=9,
            quantified_compounds=0,
        )
        data = build_experiment_summary([e])
        self.assertEqual(data["counts"]["pools"], 1)
        self.assertEqual(data["counts"]["hit_target_rows"], 0)
        self.assertEqual(data["pools"][0]["quantified_compounds"], 0)

    def test_no_pool_summary_is_valid(self):
        data = build_experiment_summary([], analysis_status="in_progress")
        self.assertEqual(data["counts"]["pools"], 0)
        self.assertEqual(data["analysis_status"], "in_progress")

    def test_old_output_files_do_not_enter_current_invocation(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            (root / "OLD_OpenASChem_results.csv").write_text("stale data")
            result = write_experiment_summary(root, [entry()])
            data = json.loads((root / "OpenASChem_summary.json").read_text(encoding="utf8"))
            self.assertEqual(result["pools"], 1)
            self.assertEqual([p["pool"] for p in data["pools"]], ["P1"])

    def test_csv_retains_all_target_rows_and_numeric_values(self):
        with tempfile.TemporaryDirectory() as directory:
            write_experiment_summary(
                directory, [entry(targets="T1_P1;T10_P1", fc="10.123456789;20")]
            )
            with (Path(directory) / "OpenASChem_all_hits.csv").open(
                encoding="utf-8-sig", newline=""
            ) as stream:
                rows = list(csv.DictReader(stream))
            self.assertEqual(len(rows), 2)
            self.assertEqual(float(rows[0]["fold_change"]), 10.123456789)
            self.assertEqual(float(rows[0]["p_value"]), 0.000123)

    def test_reference_deep_link_supports_special_characters_and_multidigit_id(self):
        link = relative_link("P one/report.html", reference="CMP & name_12")
        self.assertEqual(link, "P%20one/report.html#reference=CMP%20%26%20name_12")

    def test_missing_reference_id_links_by_compound(self):
        e = entry()
        self.assertIn("#reference=CMP_12", e["hits"][0]["report_url"])
        self.assertEqual(relative_link("P1.html", compound="CMP"), "P1.html#compound=CMP")

    def test_summary_cannot_escape_output_directory(self):
        for name in ["../outside.html", "/outside.html", "C:\\outside.html", "https://example.com"]:
            with self.assertRaises(ValueError):
                relative_link(name)

    def test_duplicate_pool_is_rejected(self):
        with self.assertRaises(ValueError):
            build_experiment_summary([entry(), entry()])

    def test_embedded_metadata_cannot_close_script(self):
        with tempfile.TemporaryDirectory() as directory:
            write_experiment_summary(
                directory, [entry(compound="</script><script>alert(1)</script>")]
            )
            text = (Path(directory) / "OpenASChem_summary.html").read_text(encoding="utf8")
            self.assertNotIn("</script><script>alert(1)", text)
            self.assertIn("\\u003c/script>", text)


if __name__ == "__main__":
    unittest.main()
