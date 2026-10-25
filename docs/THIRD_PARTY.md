# Third-party components

`src/openaschem/dialignr_numba.py` derives from a local DreamDIA-v5 development
snapshot. Its checksum, numerical subset and validation are documented in
[ALIGNMENT_BACKEND.md](ALIGNMENT_BACKEND.md). DreamDIA's project license is
[Apache-2.0](https://github.com/xmuyulab/DreamDIA/blob/main/LICENSE); a copy is
included under `src/openaschem/licenses/`. The OpenASChem MIT license remains in `LICENSE`,
with component attribution in `NOTICE`.

Plotly.py supplies the Plotly.js bundle embedded in generated offline reports.
The installed bundle's license header is preserved. Other scientific packages
are dependencies rather than vendored source; their licenses remain with them.
