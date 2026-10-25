"""Tk interfaces generated from the CLI options, sharing defaults and validation."""

import os
import queue
import subprocess
import sys
import threading
import tkinter as tk
from pathlib import Path
from tkinter import filedialog, messagebox, scrolledtext, ttk

import click

from ._version import __version__
from .cli import openaschem
from .parameters import validate_rt_limits

LABELS = {
    "file_dir": "Raw data directory",
    "chemical_pool_library": "Compound library",
    "out": "Output directory",
    "standard_prefix": "Standard prefix",
    "mz_unit": "m/z unit",
    "mz_tol": "m/z tolerance",
    "rt_tol": "RT tolerance / s",
    "acquisition_mode": "Acquisition mode",
    "quantification_method": "Quantification method",
    "background_intensity": "Background intensity",
    "min_intensity_standard": "Minimum standard intensity",
    "min_fold_change": "Minimum fold change",
    "min_fold_change_single": "Pairwise fold change",
    "p_value_cutoff": "Maximum P value",
    "p_value_cutoff_single": "Pairwise maximum P value",
    "n_threads": "Worker processes",
    "peak_selection": "Reference peak selection",
    "reference_relative_height": "Reference relative height",
    "reference_ambiguity_margin": "Reference ambiguity margin",
    "rt_left_limit": "RT left limit / s",
    "rt_right_limit": "RT right limit / s",
    "n_replicates": "Number of replicates",
}

PATH_OPTIONS = {"file_dir", "chemical_pool_library", "out"}
ADVANCED_OPTIONS = {
    "reference_relative_height",
    "reference_ambiguity_margin",
    "rt_left_limit",
    "rt_right_limit",
}
GUI_EXCLUDED_OPTIONS = {"identity_parent_library"}
HELP_TEXT = {
    "file_dir": "Folder containing centroided mzML or mzXML runs. For identification, include the standard and experimental runs to analyse.",
    "chemical_pool_library": "Compound library (CSV, TSV or XLSX) with compound identifiers, masses and pool assignments. This library is also used automatically for parent-isotope checks.",
    "out": "Folder for analysis results and interactive HTML reports. For file filling, this folder stores the filling log; copied raw files are written to the raw-data folder.",
    "standard_prefix": "Filename prefix identifying the standard runs, e.g. STD100nM. Each pool requires three standard replicates with suffixes _1, _2 and _3.",
    "mz_unit": "Unit used for m/z tolerance: ppm is relative mass error; Da is an absolute mass difference.",
    "mz_tol": "Maximum m/z deviation for extracting compound ion signals, using the selected m/z unit.",
    "rt_tol": "Retention-time tolerance in seconds used when matching peaks between the standard runs.",
    "acquisition_mode": "Ionisation polarity of the MS data. Positive or negative mode determines the expected compound ion masses.",
    "quantification_method": "Intensity uses peak height; area integrates the baseline-corrected signal over its original retention times.",
    "background_intensity": "Lower floor for quantified signals. Missing or weaker signals are replaced by this value; units follow the quantification method.",
    "min_intensity_standard": "Minimum standard peak quantity required in each of the three standard runs. Units follow the quantification method.",
    "min_fold_change": "Minimum enrichment: mean signal in the target runs divided by the mean of all other experimental runs (with a small denominator offset).",
    "min_fold_change_single": "Minimum enrichment required when the target is compared separately with every other experimental sample group.",
    "p_value_cutoff": "Maximum one-sided t-test P value for the target replicates versus all other experimental runs.",
    "p_value_cutoff_single": "Maximum one-sided t-test P value for each target-versus-other-sample comparison.",
    "n_threads": "Number of worker processes used for parallel analysis. Reduce this value if memory is limited.",
    "peak_selection": "identity evaluates candidate reference peaks using standard-run mass, isotope and peak-height evidence. top1 uses the legacy highest-peak method.",
    "reference_relative_height": "Minimum candidate height relative to the strongest accepted standard candidate (0 < value <= 1). Used by identity selection; default 0.5.",
    "reference_ambiguity_margin": "Minimum score separation between distinct competing reference peaks. A smaller separation leaves the compound ambiguous instead of forcing a reference assignment. Used by identity selection; default 0.15.",
    "rt_left_limit": "Lower reference RT boundary, in seconds (default 50). The mean peak start across the three standard runs must be greater than this value. This does not trim the raw XIC.",
    "rt_right_limit": "Upper reference RT boundary, in seconds (default 300). The mean peak end across the three standard runs must be less than this value. This does not trim the raw XIC.",
    "n_replicates": "Expected number of replicates per sample for file filling. Missing replicate files are copied from an existing run; these copies are not independent experimental replicates.",
}


class Tooltip:
    """Show delayed help without interrupting mouse or keyboard interaction."""

    def __init__(self, widget, text):
        self.widget = widget
        self.text = text
        self.pending = None
        self.window = None
        widget.bind("<Enter>", self._schedule, add="+")
        for event in ("<Leave>", "<ButtonPress>", "<Destroy>"):
            widget.bind(event, self.hide, add="+")

    def _schedule(self, event=None):
        self.hide()
        self.pending = self.widget.after(400, self.show)

    def show(self):
        self.pending = None
        if self.window is not None or not self.widget.winfo_exists():
            return
        self.window = tk.Toplevel(self.widget)
        self.window.wm_overrideredirect(True)
        ttk.Label(
            self.window, text=self.text, wraplength=420, padding=8, relief="solid", borderwidth=1
        ).pack()
        self.window.update_idletasks()
        x = min(
            self.widget.winfo_rootx() + 12,
            self.widget.winfo_screenwidth() - self.window.winfo_reqwidth() - 8,
        )
        y = min(
            self.widget.winfo_rooty() + self.widget.winfo_height() + 4,
            self.widget.winfo_screenheight() - self.window.winfo_reqheight() - 8,
        )
        self.window.geometry(f"+{max(0, x)}+{max(0, y)}")

    def hide(self, event=None):
        if self.pending is not None:
            self.widget.after_cancel(self.pending)
            self.pending = None
        if self.window is not None:
            self.window.destroy()
            self.window = None


def command_arguments(command, values):
    """Build and validate CLI arguments without invoking the analysis callback."""
    args = []
    for option in command.params:
        if not isinstance(option, click.Option) or option.hidden:
            continue
        value = values.get(option.name, "").strip()
        if value:
            args.extend([option.opts[0], value])
    with command.make_context(command.name, list(args)) as context:
        if "rt_left_limit" in context.params:
            try:
                validate_rt_limits(
                    context.params["rt_left_limit"], context.params["rt_right_limit"]
                )
            except ValueError as error:
                raise click.BadParameter(str(error), param_hint="RT limits") from error
        return args


class OpenASChemApp(tk.Tk):
    """One form and one main-thread log consumer for either supported command."""

    def __init__(self, command_name="openidentify"):
        super().__init__()
        self.title(f"OpenASChem {__version__}")
        self.command = openaschem.commands[command_name]
        self.values = {}
        self.option_labels = {}
        self.fields = {}
        self.tooltips = {}
        self.advanced_window = None
        self.events = queue.Queue()
        self.running = False
        self.columnconfigure(1, weight=1)
        self._make_form()
        self.after(100, self._poll)

    def _make_form(self):
        from PIL import Image, ImageTk

        with Image.open(Path(__file__).parent / "assets" / "logo.png") as image:
            self.window_icon = ImageTk.PhotoImage(image, master=self)
            self.iconphoto(True, self.window_icon)
            self.logo = ImageTk.PhotoImage(
                image.resize((225, 175), Image.Resampling.LANCZOS), master=self
            )
        ttk.Label(self, image=self.logo).grid(row=0, column=0, columnspan=3, pady=8)
        ttk.Label(self, text="Mingxuan Gao | mingxuan.gao@utoronto.ca").grid(
            row=1, column=0, columnspan=3, pady=4
        )
        options = [
            p
            for p in self.command.params
            if isinstance(p, click.Option) and not p.hidden and p.name not in GUI_EXCLUDED_OPTIONS
        ]
        for option in options:
            default = (
                "" if option.name in PATH_OPTIONS or option.default is None else str(option.default)
            )
            self.values[option.name] = tk.StringVar(self, value=default)
        self.advanced_options = [option for option in options if option.name in ADVANCED_OPTIONS]
        basic_options = [option for option in options if option.name not in ADVANCED_OPTIONS]
        self._make_fields(self, basic_options, start=2)
        next_row = len(basic_options) + 2
        if self.advanced_options:
            self.advanced_button = ttk.Button(
                self, text="Advanced options...", command=self._show_advanced
            )
            self.advanced_button.grid(row=next_row, column=0, columnspan=3, pady=6)
            self.tooltips["advanced"] = Tooltip(
                self.advanced_button,
                "Optional reference-peak settings. Defaults are used even when this window is not opened.",
            )
            next_row += 1
        self.run_button = ttk.Button(self, text="Run", command=self._start)
        self.run_button.grid(row=next_row, column=0, columnspan=3, pady=10)

    def _make_fields(self, container, options, start=0):
        container.columnconfigure(1, weight=1)
        for row, option in enumerate(options, start=start):
            name = option.name
            label = ttk.Label(container, text=LABELS.get(name, name.replace("_", " ")))
            label.grid(row=row, column=0, sticky="e", padx=8, pady=3)
            self.option_labels[name] = label
            self.tooltips[name] = Tooltip(label, HELP_TEXT.get(name, option.help))
            value = self.values[name]
            if isinstance(option.type, click.Choice):
                field = ttk.Combobox(
                    container, textvariable=value, values=option.type.choices, state="readonly"
                )
            else:
                field = ttk.Entry(container, textvariable=value, width=48)
            self.fields[name] = field
            field.grid(row=row, column=1, sticky="ew", padx=4, pady=3)
            if isinstance(option.type, click.Path):
                ttk.Button(container, text="Browse", command=lambda n=name: self._browse(n)).grid(
                    row=row, column=2, padx=6
                )

    def _show_advanced(self):
        if self.advanced_window is None or not self.advanced_window.winfo_exists():
            window = self.advanced_window = tk.Toplevel(self)
            window.title("OpenASChem advanced options")
            window.transient(self)
            window.protocol("WM_DELETE_WINDOW", window.withdraw)
            ttk.Label(
                window,
                text="Defaults are suitable for routine use. Hover over a name for help.",
                wraplength=480,
            ).grid(row=0, column=0, columnspan=3, padx=10, pady=10)
            self._make_fields(window, self.advanced_options, start=1)
            ttk.Button(window, text="Done", command=window.withdraw).grid(
                row=len(self.advanced_options) + 1, column=0, columnspan=3, pady=10
            )
        self.advanced_window.deiconify()
        self.advanced_window.lift()

    def _browse(self, name):
        if name in {"file_dir", "out"}:
            path = filedialog.askdirectory(parent=self)
        else:
            path = filedialog.askopenfilename(
                parent=self,
                filetypes=[("Compound libraries", "*.csv *.tsv *.xlsx"), ("All files", "*.*")],
            )
        if path:
            self.values[name].set(path)

    def _start(self):
        if self.running:
            return
        try:
            args = command_arguments(self.command, {k: v.get() for k, v in self.values.items()})
        except click.ClickException as error:
            messagebox.showerror("Invalid settings", error.format_message(), parent=self)
            return
        window = tk.Toplevel(self)
        window.title("OpenASChem logs")
        self.logs = scrolledtext.ScrolledText(window, width=110, height=32)
        self.logs.pack(fill="both", expand=True)
        # Keep the log widget alive until the job completes.
        window.protocol("WM_DELETE_WINDOW", window.withdraw)
        self.running = True
        self.run_button.configure(state="disabled")
        threading.Thread(target=self._run, args=(args,), daemon=True).start()

    def _run(self, args):
        try:
            kwargs = {"creationflags": subprocess.CREATE_NO_WINDOW} if os.name == "nt" else {}
            with subprocess.Popen(
                [sys.executable, "-m", "openaschem", self.command.name, *args],
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
                text=True,
                encoding="utf8",
                errors="replace",
                env={
                    **os.environ,
                    "PYTHONIOENCODING": "utf8",
                    "PYTHONPATH": os.pathsep.join(
                        filter(
                            None,
                            [
                                str(Path(__file__).resolve().parents[1]),
                                os.environ.get("PYTHONPATH"),
                            ],
                        )
                    ),
                },
                **kwargs,
            ) as process:
                for line in process.stdout:
                    self.events.put(("log", line))
                self.events.put(("done", process.wait()))
        except Exception as error:
            self.events.put(("error", str(error)))

    def _poll(self):
        try:
            while True:
                kind, value = self.events.get_nowait()
                if kind == "log":
                    self.logs.insert(tk.END, value)
                    self.logs.see(tk.END)
                else:
                    self.running = False
                    self.run_button.configure(state="normal")
                    if kind == "error" or value != 0:
                        messagebox.showerror(
                            "Analysis failed",
                            str(value) if kind == "error" else "See logs for details.",
                            parent=self,
                        )
                    else:
                        message = (
                            "Analysis complete. Open OpenASChem_summary.html in the output directory."
                            if self.command.name == "openidentify"
                            else "File filling complete."
                        )
                        messagebox.showinfo("Complete", message, parent=self)
        except queue.Empty:
            pass
        self.after(100, self._poll)


def main():
    OpenASChemApp().mainloop()


def filling_main():
    OpenASChemApp("openfilling").mainloop()
