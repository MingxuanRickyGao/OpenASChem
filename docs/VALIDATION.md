# Validation of 1.3.0 integration

Validated with Python 3.9 and the installed scientific dependencies on Windows.

## Repository cleanup

- 57 tests passed locally, including the original 52 checks, packaged resources,
  module CLI, a two-pool synthetic raw-data run and GUI argument validation.
- Forty-three numerical functions retain identical normalized Python ASTs after
  accounting for names, docstrings and package-relative imports. HTML templates
  are byte-identical to the preceding release.
- Full raw-data reruns retained mixed-standard 866 references / 248 A calls /
  zero B/C calls, and ChemDivA01 387 references / 2 calls. Candidate identities
  and fold changes match the preceding release. These runs check cleanup
  compatibility and are not new accuracy or controlled runtime benchmarks.
- The wheel was installed in a separate virtual environment, inheriting the
  existing scientific dependencies. Installed entry points, report assets and
  a synthetic multi-pool CLI run passed.
- CI is configured for Windows/Linux and Python 3.9, 3.11 and 3.12. That hosted
  matrix has not been executed as part of this local cleanup.

## Earlier scientific and reporting checks

The multi-pool reporting update was checked separately:

- All 52 unit/regression tests passed (39 existing plus 13 summary checks).
- Actual CLI processing of generated mzXML data passed for both one and two
  pools, including per-pool reports and automatic overview generation.
- Regenerated reports from previously validated eight-pool real-data evidence:
  216 runs, 2,082 accepted references and 15 hits (9 with FC>=5). Raw extraction
  and hit statistics were not rerun for this rendering check.
- Twenty browser checks passed: single/multiple pools, real-data counts,
  pool/target/FC/search filters, empty/invalid filters, actual filtered CSV
  download, pagination/sorting, embedded exact-reference XIC/barplots,
  zoom and PDF preparation, direct links and missing-reference handling.
  No JavaScript exceptions or external network requests occurred.
- The real-data summary is about 22 KB and loads individual pool reports
  on demand. XIC data are not copied into the overview.

The following earlier checks validate the scientific core, which this
reporting update does not change:

- 39 unit/regression checks passed, including R-free CLI import, recorded
  alignment references, original-RT area integration, shoulder/width filtering,
  identity selection, ambiguous references, multi-digit reference IDs, exact
  table-to-report annotation joins, and empty-candidate handling.
- Actual CLI analysis entry point, complete raw-data processing: mixed-standard
  12 runs / 897 compounds -> 866 accepted references, 248 A-pool hits, zero B/C
  calls. Hit identities and FC values match the experimental recommendation.
- Actual CLI analysis entry point, complete raw-data processing: ChemDivA01
  27 runs / 464 compounds -> 387 accepted references and 2 hits; hit identities
  and FC values match the previous real-data test.
- Using the integrated core, replayed saved raw-standard evidence and accepted
  candidate boundaries for mixed-standard and all eight real-data pools. All
  selected reference IDs and scores match the prior experiments. This replay
  checks integration against previously extracted evidence; it is not a second
  full raw-data run of all 216 real runs. The previous eight-pool real-data
  experiment retained 15 FC>2 calls and 9 FC>5 calls without tuning on those data.
- With only the mixed-standard input library as the parent search, selected
  reference IDs match the larger-parent-library recommendation as well.
- GUI hidden startup checked version and new control defaults. Eight browser
  checks cover accepted/review reports, inspectable ambiguous/isotope-vetoed
  candidates, empty filters and real-data plots. An actual two-page PDF export
  includes the ambiguous reference status. No external report requests occurred.

The 248 A / 0 B / 0 C benchmark counts remain exploratory parameter-optimization
results, not independently estimated specificity or FDR. Real-data transfer
preserved the candidate list and does not demonstrate a reduction in real-data
false positives without truth labels. See PEAK_IDENTITY.md for interpretation.

## Bundled installer validation (2026-10-25)

- Windows x64: installed in a path containing spaces without an external Python
  in PATH, with Python network operations denied during setup; all 63 regression
  checks passed. Both GUI windows loaded the project logo; the installed ICO
  matched the packaged logo asset. Silent uninstallation completed successfully.
- Linux x86_64: installed under Ubuntu 22.04 / WSL2 in an isolated network
  namespace, with Python network operations additionally denied; all 63
  regression checks passed. The installed logo matched the source asset and
  the CLI/GUI launchers were present. Desktop shortcut appearance was not
  visually tested across Linux desktop environments.
- Runtime relocation and application imports were tested from the installed
  environments. The installers contain the current GUI and reference-analysis
  code. This packaging check does not add a new scientific accuracy benchmark.

Large validation datasets, reports and trusted analysis caches are maintained
outside this source distribution. Portable alignment fixtures and regression
tests are included; original development datasets are not bundled.
