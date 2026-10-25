"""Experimental standard-only peak selector; never uses pool labels or FC."""

import numpy as np


def select_candidates(
    features,
    *,
    floor=0.5,
    alpha=1.0,
    beta=1.0,
    gamma=1.0,
    max_mass=float("inf"),
    min_margin=0.15,
    mass_center=None,
):
    if not np.isfinite(floor) or not 0 < floor <= 1:
        raise ValueError("Reference relative-height floor must be in (0, 1].")
    if not np.isfinite(min_margin) or min_margin < 0:
        raise ValueError("Reference ambiguity margin must be finite and nonnegative.")
    f = features.copy()
    if f.empty:
        for column in ["mass_distance", "eligible", "identity_score", "score_margin", "ambiguous"]:
            f[column] = []
        return f.copy(), f.copy(), f, mass_center
    if mass_center is None:
        mass_center = float(f.loc[f.relative_height >= 0.8, "mass_error_ppm"].median())
    f["mass_distance"] = abs(f.mass_error_ppm - mass_center)
    f["eligible"] = (f.relative_height >= floor) & (f.mass_distance <= max_mass)
    f["identity_score"] = (
        -alpha * f.isotope_loss - beta * f.mass_distance + gamma * np.log(f.relative_height)
    )
    eligible = f[f.eligible].sort_values(
        ["compound", "identity_score", "height"], ascending=[True, False, False]
    )
    best = eligible.drop_duplicates("compound").copy()
    lookup = best.set_index("compound")
    competitors = eligible[
        (abs(eligible.rt_s - eligible.compound.map(lookup.rt_s)) > 5)
        & (eligible.height >= eligible.compound.map(lookup.height) * 0.5)
    ]
    second = competitors.groupby("compound").identity_score.max()
    best["score_margin"] = best.identity_score - best.compound.map(second)
    best["ambiguous"] = best.score_margin.notna() & (best.score_margin < min_margin)
    return best[~best.ambiguous].copy(), best[best.ambiguous].copy(), f, mass_center
