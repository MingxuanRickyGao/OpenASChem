"""Validate the GUI-to-CLI handoff without creating a Tk window."""

import importlib.util
import tempfile
import unittest
from pathlib import Path

import click

from openaschem.cli import openaschem

HAS_TK = importlib.util.find_spec("tkinter") is not None
if HAS_TK:
    from openaschem.gui import command_arguments


@unittest.skipUnless(HAS_TK, "Tk is optional for command-line installations")
class GuiArgumentsTests(unittest.TestCase):
    def test_click_validation_keeps_arguments_for_subprocess(self):
        with tempfile.TemporaryDirectory() as directory:
            raw = Path(directory)
            library = raw / "library.csv"
            library.write_text("example", encoding="utf8")
            args = command_arguments(
                openaschem.commands["openidentify"],
                {
                    "file_dir": str(raw),
                    "chemical_pool_library": str(library),
                    "out": str(raw / "results"),
                    "n_threads": "2",
                },
            )
            self.assertIn("--out", args)
            self.assertIn(str(raw / "results"), args)
            self.assertEqual(args[args.index("--n_threads") + 1], "2")

    def test_invalid_settings_fail_before_launch(self):
        with self.assertRaises(click.ClickException):
            command_arguments(
                openaschem.commands["openidentify"], {"reference_relative_height": "2"}
            )
