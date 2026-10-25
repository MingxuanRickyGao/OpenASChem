"""Check GUI settings reach the analysis and advanced controls stay optional."""

import tempfile
import tkinter as tk
import unittest
from pathlib import Path
from unittest.mock import patch

import click
from click.testing import CliRunner

from openaschem.cli import openaschem
from openaschem.gui import ADVANCED_OPTIONS, HELP_TEXT, OpenASChemApp, command_arguments
from openaschem.parameters import validate_rt_limits


class RTOptionTests(unittest.TestCase):
    def test_defaults_and_nondefault_limits_reach_analysis(self):
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary)
            library = directory / "library.csv"
            library.write_text("test", encoding="utf8")
            base = [
                "openidentify",
                "--file_dir",
                temporary,
                "--chemical_pool_library",
                str(library),
                "--out",
                str(directory / "results"),
            ]
            for extra, limits in [
                ([], (50, 300)),
                (["--rt_left_limit", "25", "--rt_right_limit", "350"], (25, 350)),
            ]:
                with patch("openaschem.open_identify.open_identify") as analyse:
                    result = CliRunner().invoke(openaschem, base + extra)
                    self.assertEqual(result.exit_code, 0, result.output)
                    self.assertEqual(analyse.call_args.kwargs["rt_left_limit"], limits[0])
                    self.assertEqual(analyse.call_args.kwargs["rt_right_limit"], limits[1])

    def test_invalid_windows_fail_before_analysis(self):
        for left, right in [(300, 50), (50, 50), (-1, 300), (0, float("inf")), (float("nan"), 300)]:
            with self.subTest(left=left, right=right), self.assertRaises(ValueError):
                validate_rt_limits(left, right)


class GUIOptionTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        try:
            root = tk.Tk()
            root.withdraw()
            root.update_idletasks()
            root.destroy()
        except tk.TclError as error:
            raise unittest.SkipTest(f"Desktop display unavailable: {error}") from error

    def setUp(self):
        original = tk.Tk.__init__

        def hidden_init(window, *args, **kwargs):
            original(window, *args, **kwargs)
            window.withdraw()

        self.patch = patch.object(tk.Tk, "__init__", hidden_init)
        self.patch.start()
        self.app = OpenASChemApp()

    def tearDown(self):
        self.app.update_idletasks()
        self.app.destroy()
        self.patch.stop()

    def test_paths_empty_and_parent_library_absent(self):
        for name in ["file_dir", "chemical_pool_library", "out"]:
            self.assertEqual(self.app.values[name].get(), "")
        self.assertNotIn("identity_parent_library", self.app.values)
        self.assertTrue(ADVANCED_OPTIONS.isdisjoint(self.app.option_labels))

    def test_advanced_values_persist_and_are_forwarded(self):
        self.assertEqual(self.app.values["rt_left_limit"].get(), "50")
        self.assertEqual(self.app.values["rt_right_limit"].get(), "300")
        self.app._show_advanced()
        window = self.app.advanced_window
        window.withdraw()
        for name in ADVANCED_OPTIONS:
            self.assertEqual(self.app.option_labels[name].winfo_toplevel(), window)
        self.app.values["rt_left_limit"].set("25")
        self.app.values["reference_relative_height"].set("0.7")
        self.app._show_advanced()
        window.withdraw()
        self.assertIs(self.app.advanced_window, window)
        self.assertEqual(self.app.values["rt_left_limit"].get(), "25")
        with tempfile.TemporaryDirectory() as temporary:
            library = Path(temporary) / "library.csv"
            library.write_text("test", encoding="utf8")
            values = {key: value.get() for key, value in self.app.values.items()}
            values.update(
                file_dir=temporary,
                chemical_pool_library=str(library),
                out=str(Path(temporary) / "out"),
            )
            args = command_arguments(self.app.command, values)
            with self.app.command.make_context("openidentify", args) as context:
                self.assertEqual(context.params["rt_left_limit"], 25)
                self.assertEqual(context.params["reference_relative_height"], 0.7)
                self.assertIsNone(context.params["identity_parent_library"])
            values["rt_right_limit"] = "20"
            with self.assertRaises(click.BadParameter):
                command_arguments(self.app.command, values)

    def test_every_option_has_working_help(self):
        self.app._show_advanced()
        self.app.advanced_window.withdraw()
        self.assertEqual(set(self.app.values), set(self.app.option_labels))
        for name in self.app.values:
            self.assertTrue(HELP_TEXT[name])
            tooltip = self.app.tooltips[name]
            tooltip._schedule()
            self.assertIsNotNone(tooltip.pending)
            tooltip.hide()
            tooltip.show()
            self.assertIsNotNone(tooltip.window)
            self.assertEqual(tooltip.window.winfo_children()[0].cget("text"), HELP_TEXT[name])
            tooltip.hide()
            self.assertIsNone(tooltip.window)

    def test_file_filling_also_has_empty_paths_and_help(self):
        app = OpenASChemApp("openfilling")
        try:
            self.assertEqual(app.values["file_dir"].get(), "")
            self.assertEqual(app.values["out"].get(), "")
            self.assertEqual(set(app.tooltips), set(app.values))
        finally:
            app.update_idletasks()
            app.destroy()
