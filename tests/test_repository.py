"""Repository and installed-entry-point integration checks."""

import json
import os
import subprocess
import sys
import tempfile
import unittest
from importlib import resources
from pathlib import Path

from openaschem import __version__

ROOT = Path(__file__).resolve().parents[1]


class RepositoryTests(unittest.TestCase):
    def test_report_assets_are_packaged(self):
        package = resources.files("openaschem")
        for name in ["interactive_report_template.html", "summary_report_template.html"]:
            self.assertIn(
                "<!doctype html>",
                package.joinpath("templates", name).read_text(encoding="utf8").lower(),
            )
        self.assertTrue(package.joinpath("assets", "logo.png").is_file())
        self.assertTrue(package.joinpath("licenses", "Apache-2.0.txt").is_file())

    def test_module_cli_version_and_legacy_command_names(self):
        env = {**os.environ, "PYTHONPATH": str(ROOT / "src"), "PYTHONIOENCODING": "utf8"}
        result = subprocess.run(
            [sys.executable, "-m", "openaschem", "--version"],
            env=env,
            capture_output=True,
            encoding="utf8",
        )
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn(__version__, result.stdout)
        result = subprocess.run(
            [sys.executable, "-m", "openaschem", "--help"],
            env=env,
            capture_output=True,
            encoding="utf8",
        )
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn("openidentify", result.stdout)
        self.assertIn("openfilling", result.stdout)

    def test_synthetic_multi_pool_cli_and_summary(self):
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary)
            result = subprocess.run(
                [
                    sys.executable,
                    str(ROOT / "examples/generate_example.py"),
                    "--out",
                    str(directory / "input"),
                ],
                capture_output=True,
                text=True,
            )
            self.assertEqual(result.returncode, 0, result.stderr)
            env = {**os.environ, "PYTHONPATH": str(ROOT / "src"), "PYTHONIOENCODING": "utf8"}
            result = subprocess.run(
                [
                    sys.executable,
                    "-m",
                    "openaschem",
                    "openidentify",
                    "--file_dir",
                    str(directory / "input/raw"),
                    "--chemical_pool_library",
                    str(directory / "input/library.csv"),
                    "--standard_prefix",
                    "STD",
                    "--out",
                    str(directory / "results"),
                    "--n_threads",
                    "1",
                ],
                env=env,
                capture_output=True,
                encoding="utf8",
                timeout=120,
            )
            self.assertEqual(result.returncode, 0, result.stderr)
            data = json.loads(
                (directory / "results/OpenASChem_summary.json").read_text(encoding="utf8")
            )
            self.assertEqual(data["analysis_status"], "complete")
            self.assertEqual(data["counts"]["pools"], 2)
            self.assertEqual(data["counts"]["pool_compound_hits"], 2)
            self.assertEqual({h["target_label"] for h in data["hits"]}, {"T1"})
            for entry in data["pools"]:
                self.assertTrue((directory / "results" / entry["report_url"]).exists())
                self.assertTrue((directory / "results" / entry["review_url"]).exists())
