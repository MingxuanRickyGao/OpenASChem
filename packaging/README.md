# Building release installers

The release installers bundle Python, Tk and scientific dependencies. End users
do not need Python, Conda or an internet connection to install them. Runtime
archives and binary installers belong in GitHub Releases, not in this repository.

1. Build the current wheel with `python -m build`.
2. Prepare a relocatable Python runtime for each platform with `conda-pack`,
   installing the current `openaschem[gui]` wheel before packing. Include package
   licenses and conda package license inventories. Retain conda-pack's relocation
   metadata. The runtime must contain `openaschem-finalize.py`, copied from
   `finalize.py` in this directory.
3. For Windows, use Inno Setup 6 and a staged runtime directory:

   ```powershell
   ISCC.exe /DRuntimeDir=C:\build\windows-runtime /O=C:\build\release packaging\OpenASChem.iss
   ```

   `src/openaschem/assets/openaschem.ico` is used for the setup executable,
   installed application entry and Windows shortcuts. Keep its multi-resolution
   images when replacing the project logo.

4. For Linux x86_64, wrap the runtime archive:

   ```shell
   python packaging/build_linux.py linux-environment.tar.gz --out release/OpenASChem-1.3.0-Linux-x86_64.sh
   ```

   Installation creates a user application entry using the packaged PNG logo.
   A desktop shortcut is added when a desktop directory exists; desktop systems
   may require users to allow launching a new shortcut.

5. Test installation with network access disabled and without an external Python
   in PATH. Check CLI imports, both GUI forms, logo shortcuts and uninstallation
   on the appropriate platform. Setting `OPENASCHEM_INSTALL_TEST=1` denies Python
   network operations in the finalization step. `OPENASCHEM_NO_SHORTCUTS=1` avoids
   creating Linux shortcuts in automated tests.
6. Record fresh SHA256 hashes and upload the installers, checksum file, wheel
   and source distribution to the `v1.3.0` GitHub Release. Write `SHA256SUMS.txt`
   with LF line endings on every platform (for example, write ASCII bytes with
   `Path.write_bytes`); CRLF makes Linux `sha256sum` treat the carriage return as
   part of each filename. Include only assets actually uploaded to that release,
   and exclude the checksum file itself. Check it on Linux with
   `sha256sum --ignore-missing --check SHA256SUMS.txt`.

The source installation instructions in the main README remain available for
developers. There is no separate dependency-downloading installer.
