# Alignment regression fixture

`alignment_r_2_12.npz` contains ten deterministic synthetic XIC pairs and
their output paths from DIAlignR 2.12.0 using the original OpenASChem
`dream_align_linear` wrapper, defaults from `paramsDIAlignR()`, and the stored
linear fit/adaptive RT parameters. The cases cover zero and missing signals,
positive and negative shifts, intensity scaling, multiple peaks, nonidentity
linear fit, noise, and unequal RT grids. The inputs were generated with
NumPy RNG seed 7143 on 2026-10-21. This fixture contains no experimental data.
