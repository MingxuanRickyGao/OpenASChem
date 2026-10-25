"""Relocate and verify the bundled runtime without downloading dependencies."""

import argparse
import json
import os
import runpy
import shlex
import shutil
import sys
from pathlib import Path


def finalize(prefix):
    if os.environ.get("OPENASCHEM_INSTALL_TEST") == "1":

        def deny_network(event, arguments):
            if event in {"socket.connect", "socket.getaddrinfo", "urllib.Request"}:
                raise RuntimeError(f"Installer attempted network access: {event}")

        sys.addaudithook(deny_network)
    prefix = Path(prefix).resolve()
    environment = prefix / "environment"
    if Path(sys.prefix).resolve() != environment:
        raise RuntimeError("Run this script with the bundled Python")
    marker = prefix / "openaschem-install.json"
    marker.write_text(
        json.dumps(
            dict(
                application="OpenASChem",
                version="1.3.0",
                status="installing",
                installation_mode="bundled",
            )
        )
    )
    unpack = environment / (
        "Scripts/conda-unpack-script.py" if os.name == "nt" else "bin/conda-unpack"
    )
    print("Preparing the bundled Python environment...", flush=True)
    original_arguments, original_path = sys.argv, list(sys.path)
    try:
        sys.argv = [str(unpack)]
        sys.path.insert(0, str(unpack.parent))
        runpy.run_path(str(unpack), run_name="__main__")
    finally:
        sys.argv, sys.path[:] = original_arguments, original_path
    import tkinter  # noqa: F401
    from importlib.metadata import version

    import lxml.etree  # noqa: F401
    import numba  # noqa: F401
    import numpy  # noqa: F401
    import pandas  # noqa: F401
    import plotly  # noqa: F401
    import scipy  # noqa: F401

    import openaschem
    import openaschem.open_identify  # noqa: F401

    bindir = prefix / "bin"
    bindir.mkdir(exist_ok=True)
    if os.name == "nt":
        (bindir / "openaschem.cmd").write_text(
            '@echo off\r\n"%~dp0..\\environment\\python.exe" -E -m openaschem %*\r\n'
        )
    else:
        python = environment / "bin/python"
        commands = {
            "openaschem": "-m openaschem",
            "openaschem-gui": "-c " + shlex.quote("from openaschem.gui import main; main()"),
            "openaschem-filling-gui": "-c "
            + shlex.quote("from openaschem.gui import filling_main; filling_main()"),
        }
        for name, command in commands.items():
            launcher = bindir / name
            launcher.write_text(f'#!/bin/sh\nexec {shlex.quote(str(python))} -E {command} "$@"\n')
            launcher.chmod(0o755)
        if os.environ.get("OPENASCHEM_NO_SHORTCUTS") != "1":
            icon = prefix / "share/icons/openaschem.png"
            icon.parent.mkdir(parents=True, exist_ok=True)
            shutil.copyfile(Path(openaschem.__file__).parent / "assets/logo.png", icon)
            executable = str(bindir / "openaschem-gui")
            for old, new in [("\\", "\\\\"), ('"', '\\"'), ("`", "\\`"), ("$", "\\$"), ("%", "%%")]:
                executable = executable.replace(old, new)
            entry = (
                "[Desktop Entry]\nType=Application\nName=OpenASChem 1.3.0\n"
                "Comment=Affinity selection mass spectrometry analysis\n"
                f'Exec="{executable}"\nIcon={icon}\nTerminal=false\nCategories=Science;Education;\n'
            )
            applications = (
                Path(os.environ.get("XDG_DATA_HOME", Path.home() / ".local/share")) / "applications"
            )
            applications.mkdir(parents=True, exist_ok=True)
            (applications / "openaschem.desktop").write_text(entry)
            desktop = Path.home() / "Desktop"
            if shutil.which("xdg-user-dir"):
                import subprocess

                result = subprocess.run(["xdg-user-dir", "DESKTOP"], capture_output=True, text=True)
                candidate = Path(result.stdout.strip())
                if result.returncode == 0 and candidate.is_absolute() and candidate != Path.home():
                    desktop = candidate
            if desktop.is_dir() and desktop != Path.home():
                shortcut = desktop / "OpenASChem.desktop"
                shortcut.write_text(entry)
                shortcut.chmod(0o755)
    metadata = dict(
        application="OpenASChem",
        version=version("openaschem"),
        status="complete",
        installation_mode="bundled",
        python=sys.executable,
        dependencies={
            name: version(name)
            for name in ["numpy", "scipy", "numba", "pandas", "pyteomics", "plotly"]
        },
    )
    marker.write_text(json.dumps(metadata, indent=2))
    print(f"OpenASChem installation complete: {prefix}", flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--prefix", required=True)
    options = parser.parse_args()
    if os.name == "nt":
        sys.stdout = open(
            Path(options.prefix) / "installation.log", "a", encoding="utf8", buffering=1
        )
        sys.stderr = sys.stdout
    finalize(options.prefix)
