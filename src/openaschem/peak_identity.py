"""Standard-only identity evidence and auditable reference selection for v1.3.

This module never uses experimental enrichment or compound-pool membership
to rank a reference. Scores are heuristics, not identification probabilities.
"""

from pathlib import Path

import numpy as np
import pandas as pd
from pyteomics import mass

from .chromatogram_processing import smooth_array
from .peak_identity_evidence import OFFSETS, candidate_features, isotope_veto, vector_extract
from .peak_identity_selection import select_candidates


def read_parent_library(path):
    path = Path(path)
    if path.suffix.lower() == ".tsv":
        return pd.read_csv(path, sep="\t")
    if path.suffix.lower() == ".csv":
        return pd.read_csv(path)
    if path.suffix.lower() == ".xlsx":
        return pd.read_excel(path)
    raise ValueError("Parent library must be CSV, TSV or XLSX.")


def normalize_parent_library(table):
    table = table.copy()
    if "SGC ID for Compotent" not in table:
        table = table.rename(columns={"SGC ID for Component": "SGC ID for Compotent"})
    required = {"SGC ID for Compotent", "formula"}
    if not required.issubset(table):
        raise ValueError("Parent library requires compound IDs and formula columns.")
    table = table.dropna(subset=["SGC ID for Compotent", "formula"]).drop_duplicates(
        "SGC ID for Compotent"
    )
    calculated = []
    for _, row in table.iterrows():
        try:
            calculated.append(mass.calculate_mass(formula=row["formula"]))
        except Exception as error:
            raise ValueError(
                f"Invalid formula for parent {row['SGC ID for Compotent']}: {row['formula']}"
            ) from error
    if "MONOISOTOPIC_MASS" in table:
        table["MONOISOTOPIC_MASS"] = pd.to_numeric(
            table["MONOISOTOPIC_MASS"], errors="coerce"
        ).fillna(pd.Series(calculated, index=table.index))
    else:
        table["MONOISOTOPIC_MASS"] = calculated
    return table


def extract_standard_evidence(
    spectra_standard, chemical_pool, chemical_pool_columns, mz_unit, mz_tol, acquisition_mode
):
    ids = chemical_pool[chemical_pool_columns["COMPOUND_ID"]].tolist()
    mz = chemical_pool[chemical_pool_columns["COMPOUND_MZ"]].to_numpy(dtype=float)
    evidence = {}
    for i, spectra in enumerate(spectra_standard):
        channels = {}
        for name, offset in OFFSETS.items():
            values, weighted = vector_extract(spectra, mz + offset, mz_unit, mz_tol)
            channels[name] = {"raw": values, "weighted": weighted, "smooth": smooth_array(values)}
        evidence[i] = channels
    return dict(
        ids=ids,
        mz=mz,
        offsets=OFFSETS,
        evidence=evidence,
        mz_unit=mz_unit,
        mz_tol=mz_tol,
        ion_shift=1.007285 if acquisition_mode == "positive" else -1.007285,
    )


def select_reference_peaks(
    chemical_pool,
    chemical_pool_columns,
    spectra_standard,
    standard_peaks,
    evidence,
    parent_library,
    *,
    relative_height=0.5,
    ambiguity_margin=0.15,
):
    normalized = chemical_pool.rename(
        columns={
            chemical_pool_columns["COMPOUND_ID"]: "SGC ID for Component",
            chemical_pool_columns["FORMULA"]: "formula",
        }
    )
    features = candidate_features(
        {"spectra_standard": spectra_standard}, {"standard": standard_peaks}, evidence, normalized
    )
    center = (
        float(features.loc[features.relative_height >= 0.8, "mass_error_ppm"].median())
        if len(features)
        else None
    )
    veto_table = (
        isotope_veto(
            features,
            {"standard": standard_peaks},
            {"spectra_standard": spectra_standard},
            evidence,
            parent_library,
            center,
        )
        if len(features)
        else pd.DataFrame(columns=["key", "compound", "parent", "veto"])
    )
    veto_keys = set(veto_table.loc[veto_table.veto.astype(bool), "key"])
    selected, ambiguous, scored, _ = select_candidates(
        features.loc[~features.index.isin(veto_keys)],
        floor=relative_height,
        min_margin=ambiguity_margin,
        mass_center=center,
    )
    audit = features.copy()
    audit["isotope_veto"] = audit.index.isin(veto_keys)
    audit["identity_score"] = (
        (-audit.isotope_loss - abs(audit.mass_error_ppm - center) + np.log(audit.relative_height))
        if center is not None
        else []
    )
    audit["mass_distance"] = abs(audit.mass_error_ppm - center) if center is not None else []
    audit["selection_status"] = "alternative_not_selected"
    audit.loc[audit.relative_height < relative_height, "selection_status"] = "below_relative_height"
    audit.loc[audit.index.isin(veto_keys), "selection_status"] = "isotope_interference"
    audit.loc[audit.compound.isin(ambiguous.compound) & ~audit.isotope_veto, "selection_status"] = (
        "ambiguous"
    )
    audit.loc[audit.index.isin(selected.index), "selection_status"] = "selected"
    margins = selected["score_margin"].to_dict()
    margins.update(ambiguous["score_margin"].to_dict())
    audit["score_margin"] = audit.index.map(margins)
    decisions = []
    selected_ids = set(selected.compound)
    ambiguous_ids = set(ambiguous.compound)
    for cid in normalized["SGC ID for Component"]:
        rows = audit[audit.compound.eq(cid)]
        chosen = selected[selected.compound.eq(cid)]
        status = (
            "selected"
            if cid in selected_ids
            else "ambiguous"
            if cid in ambiguous_ids
            else "no_eligible_reference"
            if len(rows)
            else "no_standard_peak"
        )
        decisions.append(
            dict(
                compound=cid,
                status=status,
                candidate_count=len(rows),
                selected_reference=chosen.key.iloc[0] if len(chosen) else None,
                selected_rt_s=float(chosen.rt_s.iloc[0]) if len(chosen) else None,
            )
        )
    return dict(
        selected=selected,
        ambiguous=ambiguous,
        candidates=audit,
        decisions=pd.DataFrame(decisions),
        isotope_evidence=veto_table,
        mass_center_ppm=center,
    )


def library_mass_audit(chemical_pool, chemical_pool_columns, acquisition_mode):
    shift = 1.007285 if acquisition_mode == "positive" else -1.007285
    records = []
    for _, row in chemical_pool.iterrows():
        formula = row[chemical_pool_columns["FORMULA"]]
        try:
            calculated = mass.calculate_mass(formula=formula)
        except Exception as error:
            raise ValueError(
                f"Invalid formula for {row[chemical_pool_columns['COMPOUND_ID']]}: {formula}"
            ) from error
        neutral = float(row[chemical_pool_columns["COMPOUND_MZ"]]) - shift
        difference = neutral - calculated
        records.append(
            dict(
                compound=row[chemical_pool_columns["COMPOUND_ID"]],
                formula=formula,
                library_neutral_mass=neutral,
                formula_neutral_mass=calculated,
                difference_da=difference,
                needs_review=abs(difference) > 0.001,
            )
        )
    return pd.DataFrame(records)
