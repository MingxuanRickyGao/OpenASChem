# Python alignment backend

The runtime alignment backend is `src/openaschem/dialignr_numba.py`, adapted from the local
DreamDIA-v5 `DIAlignR.py` implementation by retaining its numerical kernels.
Source snapshot SHA-256: `66f3845c869dec8344918fdf6cbf157eb10e4b5709e10e85fa50adde1ca10be9`.

Supported path: hybrid alignment, mean normalization, `dotProductMasked`
similarity, Savitzky–Golay 11-point / order-4 smoothing, affine gap penalties,
overlap alignment and a soft linear RT constraint. Defaults match the previous
`paramsDIAlignR()`-based standard-run alignment path. The module also retains
DreamDIA's lowess entry point; OpenASChem currently uses the linear entry point.
This is not an implementation of the complete DIAlignR R API.

Standard S2 and S3 are aligned to S1 as before. The RT path, interpolation,
reference peak selection and downstream quantification interfaces remain in
OpenASChem. This change does not add standard-to-experimental-run alignment.

Numba is required in the runtime environment. It JIT-compiles on first use and
caches compiled kernels; a read-only installation can set `NUMBA_CACHE_DIR` to
a writable directory. The GUI uses its own Python executable for the CLI.

Neither R nor rpy2 is imported or initialized by the runtime. Historical R wrappers are not shipped. Legacy
CLI `--r_home` is hidden and ignored so old command lines remain accepted.

Upstream algorithm reference: Gupta et al., DIAlignR, Molecular & Cellular
Proteomics (2019), DOI 10.1074/mcp.TIR118.001132.

## Migration validation (2026-10-21)

Compared against the existing wrapper with DIAlignR 2.12.0, using the unchanged
OpenASChem peak filters and original RT-aware quantification:

- Ten synthetic alignment cases match the R reference paths. These include
  shifts, zero/missing signals, unequal sampling grids and a nonidentity fit;
  the reference fixture is included in `tests/fixtures` and needs no R to run.
- ChemDivA01: all 464 compounds in three standard runs, giving 928 pairwise
  alignments. Individual alignment paths are not bitwise identical to R; some
  differences remain after mapping times to scan indices. Nevertheless, all
  accepted compounds, all six boundaries per accepted standard peak, and all
  three quantities are exactly equal: 383 peaks in intensity mode and 380 in
  area mode. This establishes agreement on this dataset, not universal parity.
- CLI smoke tests used 60 compounds, three standard runs and six experimental
  runs (two targets with three replicates). Both intensity and area workflows
  completed and produced 60-row quantification matrices while `R_HOME` was
  unset and a Python import guard explicitly rejected `rpy2`. The area test
  also checked acceptance of the ignored legacy `--r_home` option.

Run the portable regression tests with `python -m unittest discover -s tests`.
The earlier RT-aware area-integration regression tests were also retained.
Numba removes the R runtime requirement; this migration does not claim a
speed improvement over DIAlignR's compiled C++ implementation.
