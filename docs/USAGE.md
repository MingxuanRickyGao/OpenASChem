# Inputs, analysis and reports

## Inputs

Use centroided positive- or negative-mode mzXML/mzML data. Each pool needs three
pooled-compound standard runs and experimental sample groups. Current filenames
follow `GROUP_POOL_REPLICATE.mzXML` (or `.mzML`), e.g. `STD_P1_1.mzXML`,
`T1_P1_1.mzXML` and `T10_P1_1.mzXML`. Use `--standard_prefix STD` here.
Replicate IDs are the last underscore-separated component; pool names occupy
the middle component. Keep pool IDs distinct (avoid one being a substring of
another). Include only raw files for this analysis in the input directory.

The CSV, TSV or XLSX analysis library requires:

| Column | Meaning |
| --- | --- |
| `SGC ID for Component` | Compound identifier |
| `Comp Acc Mass` | Neutral monoisotopic mass |
| `Pooled Well` | Pool identifier matching raw filenames |
| `SMILES` | Structure string |
| `formula` | Molecular formula |
| `SGC ID for Pool` | Pool metadata identifier |

Masses are converted to singly charged [M+H]+ or [M−H]− using the existing mass
convention. Source libraries are not overwritten. Additional columns appear
in exports with a `_LIBRARY` suffix.

## Analysis

```shell
openaschem openidentify --file_dir raw --chemical_pool_library library.csv --standard_prefix STD --out results --n_threads 2
```

Use `openaschem openidentify --help` for defaults. `--quantification_method
intensity` uses apex intensity; `area` integrates over original RT coordinates.
Choose mass, intensity, FC and P-value thresholds for the experiment.

Default `--peak_selection identity` evaluates qualifying standard candidates
using mass error, isotope evidence and relative height, deferring ambiguous
references. `--peak_selection top1` uses the preceding highest-peak reference
selection. It does not restore CV filtering, index-based width limits, R runtime
requirements or automatic static PDFs. See [PEAK_IDENTITY.md](PEAK_IDENTITY.md).

Candidate prominence / own height must be at least 0.10. The standard width
cutoff is 10 seconds between half-prominence boundaries in each original RT
grid. Experimental widths are also limited in elapsed seconds. Standard
alignment does not add general standard-to-experimental-run alignment.
See [ALIGNMENT_BACKEND.md](ALIGNMENT_BACKEND.md).

## Outputs

`OpenASChem_summary.html` combines completed pools from the current invocation,
with counts, target distributions, filters, sorting and filtered CSV export.
Refresh during processing to see newly completed pools. Use separate output
directories for concurrent analyses. Companion files are
`OpenASChem_summary.json`, `OpenASChem_pool_summary.csv` and
`OpenASChem_all_hits.csv`.

Each pool produces result/matrix CSVs and an accepted-reference HTML report.
Identity mode adds an all-candidate review report and reference evidence tables.
`ReferencePeakID` identifies the exact selected standard candidate; summary links
open that reference. Repeated IDs in different pools remain separate entries.
Multiple target annotations expand into corresponding target/FC rows; reported
P is carried from the original result row.

Open reports in recent Chrome or Edge without a web server or network access.
Inspect raw/smoothed XICs, boundaries, overlay/stacked views, RT zoom and target
or per-run bars. PDF export uses the browser print dialog and **Save as PDF**,
with XIC and barplot on separate landscape pages. Keep linked report and summary
files together when moving outputs. Display intensities are compressed float32;
analysis quantities, boundaries and RT coordinates retain their precision.

## Legacy file filling

`openaschem openfilling` and `openaschem-filling-gui` are retained. This utility
**copies an existing raw file** to a missing replicate filename and writes
`file_filling_logs.tsv`. Copies are not independent experimental replicates
and must be tracked when interpreting downstream tests.

## Python imports after cleanup

Modules now live under `openaschem`, e.g.
`from openaschem.open_identify import open_identify`.
Internal function names use snake_case. CLI commands/options and CSV/HTML output
names are preserved. Root launchers support existing source-checkout commands;
installed entry points are preferred.
