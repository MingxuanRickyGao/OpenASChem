# Changelog

### 1.3.0 public release

- Publish Windows x64 and Linux x86_64 installers with bundled Python and all dependencies.
- Use standard installer names; retire the separate dependency-downloading installers.
- Use the project logo for Windows setup, installed shortcuts and GUI window icons, and Linux application entries.
- Add installer build scripts, installation documentation and SHA256 download verification.
- Set project and citation links to `huhehaotecrystal/OpenASChem`.

### 1.3.0 GUI update

- Empty input/output path defaults, dedicated advanced settings window, and hover help for every option.
- Move reference relative height and ambiguity margin to advanced settings; expose standard RT limits (50-300 seconds by default) through GUI, CLI and API.
- Use the input compound library automatically for parent-isotope checks in the GUI; keep the optional full-library CLI/API argument.
- Bundle Python and all analysis dependencies in fully offline Windows/Linux installers.

### 1.3.0 installer compatibility

- Constrain Pyteomics to >=4.7,<5: Pyteomics 5.0 fails to import under Python 3.9.

### 1.3.0 repository cleanup

- Installable `src/openaschem` package with `pyproject.toml`, CLI and GUI entry
  points, packaged report templates/logo, source distribution and wheel.
- Unified Ruff formatting/import rules, snake_case internal functions,
  contribution guide, citation metadata and Windows/Linux CI configuration.
- Removed historical R helper, unused static PDF plotting modules/wrappers and
  unused imports. Matplotlib, Seaborn and scikit-learn are no longer required.
- Shared report number conversion and a common GUI form/log runner. GUI options
  derive from the CLI, and background workers communicate through a main-thread
  queue. CLI validation preserves arguments passed to the subprocess.
- Preserved established CLI commands/options, output names and scientific rules.
  Python implementation imports now use the `openaschem` package namespace.
- Added portable synthetic multi-pool example and distribution/GUI handoff tests;
  retained alignment fixtures and third-party provenance/license notices.

### 1.3.0 reporting update

- Automatic offline `OpenASChem_summary.html` for single- and multi-pool
  analyses, updated as each pool completes; current invocation only.
- Pool statistics, target distribution, combined hit filtering/sorting,
  pagination, filtered CSV export and companion CSV/JSON outputs.
- On-demand embedded XIC/barplot viewer and exact-reference URL links;
  interactive zoom and existing PDF export remain available.
- GUI completion message now points to the experiment overview.
- Reporting only: reference selection, enrichment and hit calls are unchanged.

## 1.3.0 — 2026-10-22

- Default standard-only identity scoring over multiple accepted reference peaks;
  configurable relative-height and ambiguity thresholds; legacy top1 mode.
- Formula-based isotope evidence, pool-level mass-error calibration and explicit
  parent M+1 interference checks with an optional larger parent library.
- Candidate decision CSVs, formula/mass review, versioned selection summary and
  an interactive all-candidate report for ambiguous/excluded references.
- CLI and GUI options for reference selection and optional parent libraries.
- Explicit ReferencePeakID in output tables; fixed multi-digit reference parsing
  and hit annotation mapping after row sorting. FC output preserves sufficient
  precision for values just above a threshold.
- Support for an empty accepted-reference set without losing report generation.
- Existing Python-only alignment, actual-RT area integration, shoulder filtering,
  elapsed-second peak-width limits, interactive XIC/barplots and browser PDF
  export are included from the current 1.2.0 development tree.
- Original 1.2.0 directory and its uncommitted changes are retained separately.
