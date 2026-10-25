# Install OpenASChem 1.3.0

Download an installer from the project's
[GitHub Release](https://github.com/huhehaotecrystal/OpenASChem/releases/tag/v1.3.0).
The installers include Python and dependencies and work without downloading
additional packages. Do not unpack them into a separate Python environment.

## Windows x64

Run `OpenASChem-1.3.0-Windows-x64-Setup.exe` and follow the wizard. Windows 10
or later on x64 is required. A user installation needs no administrator rights.
The installer creates Start Menu entries and offers a desktop shortcut using
the project logo. Keep the desktop-shortcut checkbox selected if desired.

The default location is `%LOCALAPPDATA%\Programs\OpenASChem-1.3.0`.
Uninstall through Windows Settings or the Start Menu uninstall entry. Analysis
outputs should be saved separately from the installation directory.

For command-line use, run `<install directory>\bin\openaschem.cmd --help`.

## Linux x86_64

```shell
bash OpenASChem-1.3.0-Linux-x86_64.sh
# Optional custom directory, including directories containing spaces:
bash OpenASChem-1.3.0-Linux-x86_64.sh --prefix "$PWD/OpenASChem 1.3.0"
```

The default location is `~/.local/share/openaschem/1.3.0`. Launch OpenASChem from
the application menu, or run `<install directory>/bin/openaschem-gui`.
For CLI use, run `<install directory>/bin/openaschem --help`.

A graphical desktop is needed for the GUI. The CLI can run without a display.
The application entry uses the supplied logo; a desktop shortcut is also created
when a desktop directory is available. Some desktops require selecting
"Allow Launching" on a new shortcut. The bundled runtime is tested on Ubuntu
22.04 x86_64; other Linux distributions may need compatible system libraries.

## Verify the download

Compare the installer SHA256 with `SHA256SUMS.txt` from the same release. On
Windows use `Get-FileHash -Algorithm SHA256 <installer>`; on Linux use
`sha256sum -c SHA256SUMS.txt` after downloading the files listed there.

## First launch

Select the raw-data directory, compound library and output directory in the GUI.
See [usage](USAGE.md) for input naming and report interpretation. The first
analysis may take longer while Numba compiles its numerical routines locally.
