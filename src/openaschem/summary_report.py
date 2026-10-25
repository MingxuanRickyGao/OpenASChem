"""Lightweight offline experiment summary; never recalculates hit statistics."""

import csv
import json
from pathlib import Path, PurePosixPath
from urllib.parse import quote

from ._version import __version__
from .report_utils import finite_number as number

SUMMARY_HTML = "OpenASChem_summary.html"
HIT_FIELDS = [
    "pool",
    "compound",
    "reference_peak",
    "target",
    "target_label",
    "fold_change",
    "p_value",
    "mz",
    "rt_s",
    "formula",
]
POOL_FIELDS = [
    "pool",
    "input_compounds",
    "runs",
    "quantified_compounds",
    "hit_compounds",
    "hit_target_rows",
    "ambiguous_compounds",
    "isotope_veto_candidates",
    "mass_center_ppm",
    "selection_mode",
]


def relative_link(filename, *, reference=None, compound=None):
    if filename is None:
        return None
    text = str(filename).replace("\\", "/")
    path = PurePosixPath(text)
    if path.is_absolute() or ".." in path.parts or ":" in text:
        raise ValueError("Summary report links must stay inside the output directory.")
    href = quote(text, safe="/-_.")
    if reference is not None:
        href += "#reference=" + quote(str(reference), safe="")
    elif compound is not None:
        href += "#compound=" + quote(str(compound), safe="")
    return href


def build_pool_summary(
    *,
    pool,
    results_compound,
    output_columns,
    input_compounds,
    run_count,
    quantified_compounds,
    selection_summary=None,
    report_file=None,
    review_file=None,
    output_rt_in_minutes=True,
):
    """Expand reported target/FC pairs while preserving compound and pool identity."""
    pool = str(pool)
    selection = selection_summary or {}
    hits = []
    report_file = report_file or f"{pool}_OpenASChem_report.html"
    for row in results_compound.to_dict("records"):
        compound = str(row[output_columns["COMPOUND_ID"]])
        raw_targets = row[output_columns["TARGET_PROTEIN"]]
        targets = [target.strip() for target in str(raw_targets).split(";")]
        changes = str(row["fold_change"]).split(";")
        if len(targets) != len(changes) or any(not x or x == "nan" for x in targets):
            raise ValueError(f"Target/FC annotations do not match for {pool}/{compound}.")
        ref = row.get("ReferencePeakID")
        ref = str(ref) if ref is not None and str(ref) != "nan" else None
        for target, value in zip(targets, changes):
            fc = number(value)
            if fc is None:
                raise ValueError(f"Invalid reported FC for {pool}/{compound}.")
            label = target[: -len(pool) - 1] if target.endswith("_" + pool) else target
            rt = number(row.get(output_columns.get("RT", "rt")))
            hits.append(
                dict(
                    pool=pool,
                    compound=compound,
                    reference_peak=ref,
                    target=target,
                    target_label=label,
                    fold_change=fc,
                    p_value=number(row.get("p_value")),
                    mz=number(row.get(output_columns.get("COMPOUND_MZ", "mz"))),
                    rt_s=rt * 60 if rt is not None and output_rt_in_minutes else rt,
                    formula=str(row.get(output_columns.get("FORMULA", "formula"), "")),
                    report_url=relative_link(report_file, reference=ref, compound=compound),
                    review_url=relative_link(review_file, reference=ref, compound=compound),
                )
            )
    return dict(
        pool=pool,
        input_compounds=int(input_compounds),
        runs=int(run_count),
        quantified_compounds=int(quantified_compounds),
        hit_compounds=len({hit["compound"] for hit in hits}),
        hit_target_rows=len(hits),
        ambiguous_compounds=int(selection.get("ambiguous_compounds", 0)),
        isotope_veto_candidates=int(selection.get("isotope_veto_candidates", 0)),
        mass_center_ppm=number(selection.get("mass_center_ppm")),
        selection_mode=str(selection.get("mode", "unspecified")),
        report_url=relative_link(report_file),
        review_url=relative_link(review_file),
        hits=hits,
    )


def build_experiment_summary(pools, analysis_status="complete"):
    pools = sorted(pools, key=lambda entry: entry["pool"])
    if len({entry["pool"] for entry in pools}) != len(pools):
        raise ValueError("Each pool must appear once in the experiment summary.")
    hits = [hit for entry in pools for hit in entry["hits"]]
    counts = dict(
        pools=len(pools),
        runs=sum(p["runs"] for p in pools),
        input_compound_entries=sum(p["input_compounds"] for p in pools),
        quantified_compound_entries=sum(p["quantified_compounds"] for p in pools),
        pool_compound_hits=sum(p["hit_compounds"] for p in pools),
        unique_hit_compounds=len({h["compound"] for h in hits}),
        hit_target_rows=len(hits),
        ambiguous_compounds=sum(p["ambiguous_compounds"] for p in pools),
        isotope_veto_candidates=sum(p["isotope_veto_candidates"] for p in pools),
    )
    return dict(
        schema_version=1,
        version=__version__,
        analysis_status=analysis_status,
        counts=counts,
        pools=[{k: v for k, v in entry.items() if k != "hits"} for entry in pools],
        hits=hits,
    )


def write_experiment_summary(out, pools, analysis_status="complete"):
    """Use only explicitly completed pools from this invocation, not old files."""
    directory = Path(out)
    directory.mkdir(parents=True, exist_ok=True)
    data = build_experiment_summary(pools, analysis_status)
    for name, rows, fields in [
        ("OpenASChem_all_hits.csv", data["hits"], HIT_FIELDS),
        ("OpenASChem_pool_summary.csv", data["pools"], POOL_FIELDS),
    ]:
        with (directory / name).open("w", encoding="utf-8-sig", newline="") as stream:
            writer = csv.DictWriter(stream, fieldnames=fields, extrasaction="ignore")
            writer.writeheader()
            writer.writerows(rows)
    (directory / "OpenASChem_summary.json").write_text(
        json.dumps(data, ensure_ascii=False, indent=2, allow_nan=False), encoding="utf8"
    )
    payload = json.dumps(data, ensure_ascii=False, separators=(",", ":"), allow_nan=False)
    payload = (
        payload.replace("<", "\\u003c").replace("\u2028", "\\u2028").replace("\u2029", "\\u2029")
    )
    template = (Path(__file__).parent / "templates" / "summary_report_template.html").read_text(
        encoding="utf8"
    )
    path = directory / SUMMARY_HTML
    path.write_text(template.replace("/*__DATA__*/", payload), encoding="utf8")
    return dict(
        path=str(path),
        pools=len(pools),
        hit_compounds=data["counts"]["pool_compound_hits"],
        hit_target_rows=len(data["hits"]),
        bytes=path.stat().st_size,
    )
