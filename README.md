# OpenASChem

Compound-centric analysis of high-throughput pull-down affinity selection mass
spectrometry (AS-MS) screens. Version **1.3.0**.

OpenASChem connects pooled-compound standards, cross-run peak quantification,
target enrichment and interactive candidate review. It processes mzML/mzXML files and CSV/TSV/XLSX compound libraries, supports multiple pools,
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

## Test with experimental mixed-standard data

The experimental mixed-standard dataset is available on Zenodo:
[DOI: 10.5281/zenodo.23851067](https://doi.org/10.5281/zenodo.23851067).
Download **both** `standard_mix_rawdata.zip` and `standard_mix_library.tsv`.
These are experimental LC-MS measurements, separate from the synthetic example above.

### Prepare the inputs

Extract the ZIP archive. For the commands below, place the 12 extracted mzXML
files directly in `standard_mix_test/raw/`, and place the library in
`standard_mix_test/standard_mix_library.tsv`. Keep the original run filenames:

- Standards: `STD_mix_1.mzXML`, `STD_mix_2.mzXML`, `STD_mix_3.mzXML`.
- Experimental groups: `Sample1_mix_1.mzXML` through `Sample1_mix_3.mzXML`,
  with the same three-replicate pattern for `Sample2` and `Sample3`.

The raw-data input must point to the folder containing these files, rather than
an enclosing folder or the ZIP archive. Use a separate output folder.

### Run in the GUI

Launch OpenASChem and set the following fields:

| GUI field | Value |
| --- | --- |
| Raw data directory | The extracted folder containing all 12 mzXML files |
| Compound library | `standard_mix_library.tsv` |
| Output directory | A separate folder, e.g. `standard_mix_results` |
| Standard prefix | `STD` |
| Acquisition mode | `positive` |
| Quantification method | `intensity` |
| Reference peak selection | `identity` |
| Minimum fold change | `5` |
| Maximum P value | `0.01` |
| Pairwise fold change | `1` |
| Pairwise maximum P value | `0.01` |
| Worker processes | `2` |

Keep all remaining settings, including advanced options, at their defaults,
and click **Run**. In particular, change the default standard prefix from
`STD100nM` to `STD`, and the default maximum P value from `0.05` to `0.01`.

### Run from the command line

If the `openaschem` command is available in your environment, run:

```shell
openaschem openidentify --file_dir standard_mix_test/raw --chemical_pool_library standard_mix_test/standard_mix_library.tsv --standard_prefix STD --out standard_mix_results --peak_selection identity --min_fold_change 5 --p_value_cutoff 0.01 --min_fold_change_single 1 --p_value_cutoff_single 0.01 --n_threads 2
```

Run this command from the parent directory of `standard_mix_test`, or replace
the input/output paths with your own paths. Quote paths containing spaces.

### Inspect the results

Open `standard_mix_results/OpenASChem_summary.html` in a recent Chrome or Edge
browser. Follow a compound's report link to inspect its standard and experimental
XICs, selected peak boundaries, and sample-level or per-run barplots. The plots
support RT zoom; use the PDF export control to print or save the plots as PDF.
Keep the summary and linked per-pool HTML reports in the same output folder.

For the previously validated v1.3.0 mixed-standard benchmark, these settings
produced **866 accepted references and 248 hits**, with all hits from Pool A
(`A01`) and none from Pools B or C. The source-pool labels are stored in the
library column `Orginal Pooled Well` and exported as
`Orginal Pooled Well_LIBRARY`; all three source pools share the analysis pool
ID `mix`. Check `mix_OpenASChem_results.csv` for hits and
`mix_OpenASChem_matrix.csv` for quantified compounds. Accepted references include
compounds that do not meet the hit thresholds, so the matrix has more rows than
the hit table.

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
