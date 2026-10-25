"""Reusable extraction of the unchanged mixed-standard identity features."""

import math

import numpy as np
import pandas as pd
from pyteomics import mass

OFFSETS = {
    "mono": 0.0,
    "C13": 1.00335483507,
    "C13_2": 2.00670967014,
    "Cl37": 1.99704992,
    "Br81": 1.9979521,
    "S34": 1.99579583,
    "minus_C13": -1.00335483507,
}


def extraction_bounds(targets, mz_unit="ppm", mz_tol=3.0):
    targets = np.asarray(targets, dtype=float)
    if mz_unit == "ppm":
        return targets / (1 + mz_tol / 1e6), targets * (1 + mz_tol / 1e6)
    if mz_unit == "Da":
        return targets - mz_tol / 2, targets + mz_tol / 2
    raise ValueError("m/z unit must be ppm or Da")


def vector_extract(spectra, targets, mz_unit="ppm", mz_tol=3.0):
    lo, hi = extraction_bounds(targets, mz_unit, mz_tol)
    values = np.zeros((len(targets), len(spectra.spectra)))
    weighted = np.zeros_like(values)
    for i, (m, intensities) in enumerate(spectra.spectra):
        left = np.searchsorted(m, lo, side="left")
        right = np.searchsorted(m, hi, side="right")
        sums = np.r_[0.0, np.cumsum(intensities, dtype=np.float64)]
        weights = np.r_[0.0, np.cumsum(m * intensities, dtype=np.float64)]
        values[:, i] = sums[right] - sums[left]
        weighted[:, i] = weights[right] - weights[left]
    return values, weighted


def candidate_features(inputs, peaks, evidence, library):
    d = inputs
    a = peaks
    e = evidence
    lib = library.set_index("SGC ID for Component", drop=False)
    lookup = {cid: i for i, cid in enumerate(e["ids"])}
    records = []
    for cid, candidates in a["standard"].items():
        ix = lookup[cid]
        target = e["mz"][ix]
        comp = mass.Composition(formula=lib.loc[cid, "formula"])
        expected = {
            "C13": comp.get("C", 0) * 0.0107 / 0.9893,
            "C13_2": math.comb(comp.get("C", 0), 2) * (0.0107 / 0.9893) ** 2,
            "Cl37": comp.get("Cl", 0) * 0.2424 / 0.7576,
            "Br81": comp.get("Br", 0) * 0.4931 / 0.5069,
            "S34": comp.get("S", 0) * 0.0425 / 0.9499,
        }
        for j, s in enumerate(candidates):
            rec = dict(
                compound=cid, candidate=j + 1, key=f"{cid}_{j + 1}", formula=lib.loc[cid, "formula"]
            )
            heights = []
            apexes = []
            errors = []
            losses = []
            ratios = {k: [] for k in expected}
            corrs = {k: [] for k in expected}
            predictions = {}
            for i in range(3):
                rt = np.asarray(d["spectra_standard"][i].rt_list)
                mask = (rt >= s[f"lb{i + 1}"]) & (rt <= s[f"rb{i + 1}"])
                trace = e["evidence"][i]["mono"]["smooth"][ix]
                apex = np.flatnonzero(mask)[np.argmax(trace[mask])]
                heights.append(trace[apex])
                apexes.append(rt[apex])
                window = (rt >= s[f"lb{i + 1}"] - 1) & (rt <= s[f"rb{i + 1}"] + 1)
                denominator = e["evidence"][i]["mono"]["raw"][ix][mask].sum()
                centroid = e["evidence"][i]["mono"]["weighted"][ix][mask].sum() / max(
                    denominator, 1.0
                )
                errors.append((centroid - target) / target * 1e6)
                runloss = []
                for channel in expected:
                    lo, hi = extraction_bounds(
                        [target + OFFSETS[channel]], e.get("mz_unit", "ppm"), e.get("mz_tol", 3.0)
                    )
                    predicted = sum(
                        v for k, v in expected.items() if lo[0] <= target + OFFSETS[k] <= hi[0]
                    )
                    predictions[channel] = predicted
                    iso = e["evidence"][i][channel]["raw"][ix]
                    ratio = iso[mask].sum() / max(denominator, 1.0)
                    ratios[channel].append(ratio)
                    sm = e["evidence"][i][channel]["smooth"][ix]
                    corr = (
                        float(np.corrcoef(trace[window], sm[window])[0, 1])
                        if np.std(trace[window]) > 0 and np.std(sm[window]) > 0
                        else 0.0
                    )
                    corrs[channel].append(corr)
                    if channel == "C13" or (
                        predicted >= 0.025 and channel in ["Cl37", "Br81", "S34"]
                    ):
                        runloss.append(abs(np.log((ratio + 0.005) / (predicted + 0.005))))
                losses.append(np.mean(runloss))
            rec.update(
                rt_s=np.mean(apexes),
                height=sum(heights),
                mass_error_ppm=np.median(errors),
                mass_error_abs_ppm=np.median(np.abs(errors)),
                isotope_loss=np.median(losses),
            )
            for i in range(3):
                rec[f"std{i + 1}_height"] = heights[i]
                rec[f"std{i + 1}_mass_error_ppm"] = errors[i]
                rec[f"std{i + 1}_isotope_loss"] = losses[i]
                rec[f"std{i + 1}_C13_corr"] = corrs["C13"][i]
            for channel in expected:
                rec[channel + "_ratio"] = np.median(ratios[channel])
                rec[channel + "_expected"] = predictions[channel]
                rec[channel + "_corr"] = np.median(corrs[channel])
            records.append(rec)
    if not records:
        return pd.DataFrame(
            columns=[
                "key",
                "compound",
                "candidate",
                "formula",
                "rt_s",
                "height",
                "mass_error_ppm",
                "isotope_loss",
                "relative_height",
            ]
        ).set_index("key", drop=False)
    f = pd.DataFrame(records).set_index("key", drop=False)
    f["relative_height"] = f.height / f.groupby("compound").height.transform("max")
    return f


def isotope_veto(features, peaks, inputs, evidence, parent_library, center):
    f = features
    a = peaks
    d = inputs
    e = evidence
    full = parent_library.copy()
    full["nC"] = [mass.Composition(formula=x).get("C", 0) for x in full.formula]
    parent_mz = full.MONOISOTOPIC_MASS.to_numpy() + e.get("ion_shift", 1.007285)
    lookup = {cid: i for i, cid in enumerate(e["ids"])}
    records = []
    for cid, g in f.groupby("compound", sort=False):
        ix = lookup[cid]
        target = e["mz"][ix]
        ppm = (parent_mz + 1.00335483507 - target) / target * 1e6
        lo, hi = extraction_bounds([target], e.get("mz_unit", "ppm"), e.get("mz_tol", 3.0))
        for parent_ix in np.flatnonzero(
            (parent_mz + 1.00335483507 >= lo[0]) & (parent_mz + 1.00335483507 <= hi[0])
        ):
            pr = full.iloc[parent_ix]
            nc = pr.nC
            if nc < 2:
                continue
            for key, row in g.iterrows():
                s = a["standard"][cid][int(row.candidate) - 1]
                correlations = []
                ratios = []
                next_ratios = []
                for i in range(3):
                    rt = np.asarray(d["spectra_standard"][i].rt_list)
                    mask = (rt >= s[f"lb{i + 1}"]) & (rt <= s[f"rb{i + 1}"])
                    window = (rt >= s[f"lb{i + 1}"] - 2) & (rt <= s[f"rb{i + 1}"] + 2)
                    parent = e["evidence"][i]["minus_C13"]["smooth"][ix]
                    mono = e["evidence"][i]["mono"]["smooth"][ix]
                    corr = (
                        float(np.corrcoef(parent[window], mono[window])[0, 1])
                        if np.std(parent[window]) > 0 and np.std(mono[window]) > 0
                        else 0.0
                    )
                    correlations.append(corr)
                    prsum = e["evidence"][i]["minus_C13"]["raw"][ix][mask].sum()
                    mosum = e["evidence"][i]["mono"]["raw"][ix][mask].sum()
                    nextsum = e["evidence"][i]["C13"]["raw"][ix][mask].sum()
                    ratios.append(mosum / max(prsum, 1.0))
                    next_ratios.append(nextsum / max(mosum, 1.0))
                ratio = float(np.median(ratios))
                next_ratio = float(np.median(next_ratios))
                corr = float(np.median(correlations))
                predicted = nc * 0.0107 / 0.9893
                predicted_next = (nc - 1) / 2 * 0.0107 / 0.9893
                iso_residual = abs(row.mass_error_ppm - ppm[parent_ix] - center)
                mono_residual = abs(row.mass_error_ppm - center)
                veto = (
                    corr >= 0.98
                    and abs(np.log((ratio + 0.001) / (predicted + 0.001))) < np.log(1.3)
                    and abs(np.log((next_ratio + 0.001) / (predicted_next + 0.001))) < np.log(1.4)
                    and iso_residual <= 1.0
                    and mono_residual - iso_residual >= 1.0
                )
                records.append(
                    dict(
                        key=key,
                        compound=cid,
                        parent=pr["SGC ID for Compotent"],
                        rt_s=row.rt_s,
                        correlation=corr,
                        observed_ratio=ratio,
                        predicted_ratio=predicted,
                        observed_next_ratio=next_ratio,
                        predicted_next_ratio=predicted_next,
                        isotope_mass_residual=iso_residual,
                        mono_mass_residual=mono_residual,
                        veto=bool(veto),
                    )
                )
    return pd.DataFrame(
        records,
        columns=[
            "key",
            "compound",
            "parent",
            "rt_s",
            "correlation",
            "observed_ratio",
            "predicted_ratio",
            "observed_next_ratio",
            "predicted_next_ratio",
            "isotope_mass_residual",
            "mono_mass_residual",
            "veto",
        ],
    )
