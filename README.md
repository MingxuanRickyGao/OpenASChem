# OpenASChem

Compound-centric analysis of high-throughput pull-down affinity selection mass
spectrometry (AS-MS) screens. Version **1.3.0**.

OpenASChem connects pooled-compound standards, cross-run peak quantification,
target enrichment and interactive candidate review. It processes centroided
mzML/mzXML files and CSV/TSV/XLSX compound libraries, supports multiple pools,
and produces offline XIC/barplot reports and an experiment-wide overview.

## Install

For Windows and Linux x86_64, download an installer from
[Releases](https://github.com/huhehaotecrystal/OpenASChem/releases/latest):

- **Windows:** run `OpenASChem-1.3.0-Windows-x64-Setup.exe`. The wizard offers a
  desktop shortcut and creates Start Menu entries with the OpenASChem logo.
- **Linux:** run `bash OpenASChem-1.3.0-Linux-x86_64.sh`, optionally adding
  `--prefix /your/install/directory`. Launch the application shortcut or
  `<install directory>/bin/openaschem-gui`.

Both installers include Python and all required dependencies; installation
requires no downloads or existing Python/Conda environment. A graphical desktop
is required for the GUI. See [installation details](docs/INSTALLATION.md).

### Install from source

Python 3.9–3.12 is supported. From this source checkout:

```shell
python -m pip install ".[gui]"
```

For command-line use without the GUI, run `python -m pip install .`.
Alignment runs in Python/Numba; R and rpy2 are not required.
The GUI needs Tk (provided by many Python installers; on Ubuntu/Debian,
install the matching `python3-tk` system package).

Alternatively, run `conda env create -f openaschem_environment.yml` from this
directory, followed by `conda activate openaschem`.

## Run

```shell
openaschem --version
openaschem openidentify --help
openaschem-gui
```

Existing launchers remain available: `python OpenASChem.py`,
`python openaschem_gui.py` and `python openfilling_gui.py`.
`python -m openaschem` provides the same CLI as `openaschem`.

For a portable two-pool example:

```shell
python examples/generate_example.py --out example_data
openaschem openidentify --file_dir example_data/raw --chemical_pool_library example_data/library.csv --standard_prefix STD --out example_results --n_threads 1
```

This generates artificial Gaussian MS1 data to exercise the software; it is
not a scientific accuracy benchmark. Open `example_results/OpenASChem_summary.html`
in a recent Chrome or Edge browser to filter hits, view cross-run XICs/barplots,
zoom the plots and export PDFs.

See [input conventions and outputs](docs/USAGE.md),
[reference selection](docs/PEAK_IDENTITY.md),
[alignment provenance](docs/ALIGNMENT_BACKEND.md), and
[validation](docs/VALIDATION.md).

## Development

```shell
python -m pip install -e ".[gui,dev]"
ruff check .
ruff format --check .
python -m pytest
python -m build
```

Code lives in `src/openaschem/`; HTML templates and the logo are included in
installed packages. Tests and synthetic fixtures live in `tests/`.
Large experimental data, development caches and generated reports are excluded.
See [CONTRIBUTING.md](CONTRIBUTING.md) and [CHANGELOG.md](CHANGELOG.md).

## License and attribution

OpenASChem-authored code is provided under the [MIT license](LICENSE).
The adapted DreamDIA alignment component retains Apache-2.0 attribution;
see [NOTICE](NOTICE) and [third-party notices](docs/THIRD_PARTY.md).
