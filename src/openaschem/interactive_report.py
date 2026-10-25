"""Offline, compound-centric XIC and quantification reports.

Only report generation changes here; peak assignments, quantitative matrices
and hit statistics are supplied by the analysis pipeline without recalculation.
"""

from __future__ import annotations

import base64
import json
import time
import zlib
from functools import lru_cache
from pathlib import Path

import numpy as np

from .report_utils import finite_number as _number


def encode_trace(values):
    """Losslessly compress float32 display intensities; quantitation stays float64."""
    array = np.asarray(values, dtype="<f4")
    if array.ndim != 1 or not np.isfinite(array).all():
        raise ValueError("Report traces must be finite one-dimensional arrays.")
    return base64.b64encode(zlib.compress(array.tobytes(), level=1)).decode("ascii")


@lru_cache(maxsize=1)
def _plotly_bundle():
    # Scatter/bar reports need no map resources. Keep the document fully offline.
    import re

    from plotly.offline import get_plotlyjs

    return re.sub(
        r'https://unpkg\.com/maplibre-gl@[^"\s]+/dist/maplibre-gl\.css',
        "data:text/css,",
        get_plotlyjs(),
    ).replace("</script", "<\\/script")


def build_report_data(
    *,
    pool,
    quantification_method,
    standard_prefix,
    quantification_matrix,
    fold_changes,
    results_compound,
    output_columns,
    results_all,
    spectra_standard,
    XICs_standard,
    orig_XICs_standard,
    spectra_runs,
    XICs_runs,
    orig_XICs_runs,
    standard_peak_boundaries_and_quantities,
    runs_peak_boundaries_and_quantities,
    selection_annotations=None,
    selection_summary=None,
    report_kind="accepted",
):
    """Build a deterministic payload, retaining exact sample identities and bounds."""
    id_column = output_columns["COMPOUND_ID"]
    metadata = {}
    for _, row in results_all.iterrows():
        metadata.setdefault(str(row[id_column]), row)
    hits = {}
    if len(results_compound) != len(fold_changes):
        raise ValueError("Hit annotations and fold-change rows must have matching lengths.")
    annotations = (
        results_compound.set_index("ReferencePeakID", drop=False)
        if "ReferencePeakID" in results_compound
        else None
    )
    if annotations is not None and (
        not annotations.index.is_unique or set(annotations.index) != set(fold_changes.index)
    ):
        raise ValueError("Hit reference IDs must uniquely match fold-change rows.")
    for key, (_, legacy_row) in zip(fold_changes.index, results_compound.iterrows()):
        row = annotations.loc[key] if annotations is not None else legacy_row
        target_text = row.get(output_columns["TARGET_PROTEIN"], "")
        hits[str(key)] = dict(
            targets=[t.strip() for t in str(target_text).split(";") if t.strip()],
            fold_change=str(row.get("fold_change", "")),
            p_value=_number(row.get("p_value")),
        )

    standard_names = [f"{standard_prefix}_{i}" for i in (1, 2, 3)]
    experiment_names = sorted(spectra_runs)
    runs = []
    for i, name in enumerate(standard_names):
        runs.append(
            dict(
                id=name,
                group=standard_prefix,
                label=f"Standard {i + 1}",
                standard=True,
                rt=np.asarray(spectra_standard[i].rt_list, dtype=float).tolist(),
            )
        )
    for name in experiment_names:
        # Only strip the replicate suffix. Substring matching can confuse targets.
        group = name.rsplit("_", 1)[0]
        label = group[: -len(pool) - 1] if group.endswith("_" + pool) else group
        runs.append(
            dict(
                id=name,
                group=group,
                label=label,
                standard=False,
                rt=np.asarray(spectra_runs[name].rt_list, dtype=float).tolist(),
            )
        )
    for run in runs:
        rt = np.asarray(run["rt"])
        if not len(rt) or not np.isfinite(rt).all() or np.any(np.diff(rt) <= 0):
            raise ValueError(f"Invalid RT axis in report run: {run['id']}")

    compounds, records = {}, []
    for key in quantification_matrix.index:
        compound, rank_text = str(key).rsplit("_", 1)
        rank = int(rank_text) - 1
        if compound not in compounds:
            traces = []
            for i in range(3):
                traces.append(
                    dict(
                        raw=encode_trace(orig_XICs_standard[i][compound]),
                        smooth=encode_trace(XICs_standard[i][compound]),
                    )
                )
            for name in experiment_names:
                traces.append(
                    dict(
                        raw=encode_trace(orig_XICs_runs[name][compound]),
                        smooth=encode_trace(XICs_runs[name][compound]),
                    )
                )
            for run, trace in zip(runs, traces):
                # Validate lengths against full original RT, without cropping peaks.
                raw_array = (
                    orig_XICs_standard[standard_names.index(run["id"])][compound]
                    if run["standard"]
                    else orig_XICs_runs[run["id"]][compound]
                )
                smooth_array = (
                    XICs_standard[standard_names.index(run["id"])][compound]
                    if run["standard"]
                    else XICs_runs[run["id"]][compound]
                )
                if len(raw_array) != len(run["rt"]) or len(smooth_array) != len(run["rt"]):
                    raise ValueError(f"Trace/RT length mismatch: {compound}, {run['id']}")
            row = metadata.get(compound, {})
            compounds[compound] = dict(
                id=compound,
                mz=_number(row.get(output_columns["COMPOUND_MZ"])),
                formula=str(row.get(output_columns["FORMULA"], "")),
                smiles=str(row.get(output_columns["SMILES"], "")),
                traces=traces,
            )
        standard = standard_peak_boundaries_and_quantities[compound][rank]
        experimental = runs_peak_boundaries_and_quantities[compound][rank]
        bounds, peak_quantities, quantities = [], [], []
        for i, name in enumerate(standard_names, 1):
            bounds.append([_number(standard[f"lb{i}"]), _number(standard[f"rb{i}"])])
            peak_quantities.append(_number(standard[f"q{i}"]))
            quantities.append(_number(quantification_matrix.at[key, name]))
        for name in experiment_names:
            bounds.append([_number(experimental[0][name]), _number(experimental[1][name])])
            peak_quantities.append(_number(experimental[2][name]))
            quantities.append(_number(quantification_matrix.at[key, name]))
        hit = hits.get(str(key))
        records.append(
            dict(
                key=str(key),
                compound=compound,
                rank=rank + 1,
                hit=bool(hit),
                targets=hit["targets"] if hit else [],
                fold_change=hit["fold_change"] if hit else "",
                p_value=hit["p_value"] if hit else None,
                bounds=bounds,
                quantities=quantities,
                peak_quantities=peak_quantities,
            )
        )
        if selection_annotations and str(key) in selection_annotations:
            records[-1]["selection"] = selection_annotations[str(key)]
    records.sort(key=lambda r: (not r["hit"], r["compound"], r["rank"]))
    return dict(
        schema_version=1,
        pool=str(pool),
        quantification_method=quantification_method,
        runs=runs,
        compounds=compounds,
        records=records,
        display_precision="float32 intensities; original RT and quantities preserved",
        encoding="base64-zlib-float32-little-endian",
        selection_summary=selection_summary,
        report_kind=report_kind,
    )


def write_interactive_report(*, out, filename=None, **kwargs):
    """Write one self-contained report per pool; no browser or PDF renderer runs."""
    started = time.perf_counter()
    data = build_report_data(**kwargs)
    template = (Path(__file__).parent / "templates" / "interactive_report_template.html").read_text(
        encoding="utf-8"
    )
    payload = json.dumps(data, ensure_ascii=False, separators=(",", ":"), allow_nan=False)
    payload = (
        payload.replace("<", "\\u003c").replace("\u2028", "\\u2028").replace("\u2029", "\\u2029")
    )
    html = template.replace("/*__PLOTLY__*/", _plotly_bundle()).replace("/*__DATA__*/", payload)
    path = Path(out) / (filename or f"{kwargs['pool']}_OpenASChem_report.html")
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(html, encoding="utf-8")
    return dict(
        path=str(path),
        records=len(data["records"]),
        hits=sum(r["hit"] for r in data["records"]),
        bytes=path.stat().st_size,
        seconds=time.perf_counter() - started,
    )
