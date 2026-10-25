"""Compound-centric AS-MS analysis workflow."""

import json
import math
import multiprocessing
import os

from ._version import __version__
from .interactive_report import write_interactive_report
from .ms_file_processing import load_rawdata
from .open_identify_utils import (
    build_hit_table,
    build_quantification_matrix,
    build_quantification_table,
    calculate_fold_changes,
    collect_experimental_runs,
    extract_pool_xics,
    identify_peaks,
    load_experimental_runs,
)
from .parameters import DEFAULT_RT_LEFT_LIMIT, DEFAULT_RT_RIGHT_LIMIT, validate_rt_limits
from .peak_identity import (
    extract_standard_evidence,
    library_mass_audit,
    normalize_parent_library,
    read_parent_library,
    select_reference_peaks,
)
from .summary_report import build_pool_summary, write_experiment_summary
from .utils import (
    check_raw_data_files,
    get_openaschem_logger,
    load_chemical_pool,
    split_list,
    welcome,
)


def open_identify(
    file_dir,
    chemical_pool_library,
    standard_prefix,
    out,
    mz_unit,
    mz_tol,
    rt_tol,
    acquisition_mode,
    quantification_method,
    background_intensity,
    min_intensity_standard,
    min_fold_change,
    min_fold_change_single,
    p_value_cutoff,
    p_value_cutoff_single,
    n_threads,
    peak_selection="identity",
    reference_relative_height=0.5,
    reference_ambiguity_margin=0.15,
    identity_parent_library=None,
    rt_left_limit=DEFAULT_RT_LEFT_LIMIT,
    rt_right_limit=DEFAULT_RT_RIGHT_LIMIT,
):
    rt_left_limit, rt_right_limit = validate_rt_limits(rt_left_limit, rt_right_limit)
    if peak_selection not in ("identity", "top1"):
        raise ValueError("peak_selection must be identity or top1.")
    if not math.isfinite(reference_relative_height) or not 0 < reference_relative_height <= 1:
        raise ValueError("reference_relative_height must be in (0, 1].")
    if not math.isfinite(reference_ambiguity_margin) or reference_ambiguity_margin < 0:
        raise ValueError("reference_ambiguity_margin must be finite and nonnegative.")
    if n_threads < 1 or not math.isfinite(mz_tol) or mz_tol <= 0:
        raise ValueError("n_threads and mz_tol must be positive.")
    logger = get_openaschem_logger()
    welcome(logger)

    logger.info("Output files will be at: %s" % out)
    if not os.path.exists(out):
        os.makedirs(out, exist_ok=True)

    # static parameters for OpenASChem
    logger.info("Alignment backend: Python/Numba (DIAlignR-compatible hybrid alignment)")
    slope = 1
    intercept = 0
    adaptive_rt = 5
    merge_order_standard = [("S1", "S2"), ("S1", "S3")]
    key_runs_standard = ["S1", "S1"]
    queue_size = 512
    chemical_pool_columns = {
        "COMPOUND_ID": "SGC ID for Component",
        "COMPOUND_MZ": "Comp Acc Mass",
        "POOL_NAME": "Pooled Well",
        "SMILES": "SMILES",
        "FORMULA": "formula",
        "POOL_ID": "SGC ID for Pool",
    }
    output_columns = {
        "COMPOUND_MZ": "mz",
        "RT": "rt",
        "POOL_NAME": "ChemicalID",
        "TARGET_PROTEIN": "ProtID",
        "SMILES": "SMILES",
        "ADDUCT_TYPE": "Adducts",
        "FORMULA": "formula",
        "COMPOUND_ID": "SGC ID for Component",
        "POOL_ID": "SGC ID for Pool",
    }
    n_compounds_for_output = 100
    max_peak_width = 10
    max_standard_peak_width = 10.0  # seconds, evaluated in each standard run
    peak_refinement_index_range = 5
    delta_rt_cutoff = 13
    n_candidate_peaks = 10000 if peak_selection == "identity" else 1
    epsilon = 1
    output_RT_in_minutes = True
    output_results_suffix = "_OpenASChem_results.csv"
    output_matrix_suffix = "_OpenASChem_matrix.csv"

    logger.info("Load chemical pool library: %s" % chemical_pool_library)
    chemical_pool_total = load_chemical_pool(
        chemical_pool_library, acquisition_mode, chemical_pool_columns
    )
    if peak_selection == "identity":
        parents = normalize_parent_library(
            read_parent_library(identity_parent_library)
            if identity_parent_library
            else chemical_pool_total
        )
    bonus_columns = [
        col
        for col in list(chemical_pool_total.columns)
        if col not in list(chemical_pool_columns.values())
    ]

    logger.info("Load raw data...")
    raw_data_info = check_raw_data_files(file_dir, standard_prefix, logger)
    completed_pools = []
    overview = write_experiment_summary(
        out, completed_pools, analysis_status="in_progress" if raw_data_info else "complete"
    )

    with multiprocessing.Manager() as mgr:
        for pool_idx, pool in enumerate(raw_data_info):
            logger.info(
                "Start processing (%d / %d) pool: %s..." % (pool_idx + 1, len(raw_data_info), pool)
            )

            standard_files, experimental_files, experimental_runs = raw_data_info[pool]
            chemical_pool = chemical_pool_total[
                chemical_pool_total[chemical_pool_columns["POOL_NAME"]] == pool
            ]
            if chemical_pool.empty:
                raise ValueError(f"No library compounds match raw-data pool {pool}.")
            standard_prefix_pool = standard_prefix + "_" + pool

            spectra_standard = [
                load_rawdata(os.path.join(file_dir, file)) for file in standard_files
            ]
            XICs_standard_package = [
                extract_pool_xics(
                    spectra_run, chemical_pool, mz_unit, mz_tol, chemical_pool_columns
                )
                for spectra_run in spectra_standard
            ]
            XICs_standard = [xic[0] for xic in XICs_standard_package]
            orig_XICs_standard = [xic[1] for xic in XICs_standard_package]
            identity_evidence = (
                extract_standard_evidence(
                    spectra_standard,
                    chemical_pool,
                    chemical_pool_columns,
                    mz_unit,
                    mz_tol,
                    acquisition_mode,
                )
                if peak_selection == "identity"
                else None
            )

            spectra_runs, XICs_runs, orig_XICs_runs = mgr.dict(), mgr.dict(), mgr.dict()
            spectra_and_XIC_queue = mgr.JoinableQueue(queue_size)
            loaders = []
            experimental_files_chunks = split_list(experimental_files, n_threads)

            for i, experimental_files_chunk in enumerate(experimental_files_chunks):
                if experimental_files_chunk:
                    p = multiprocessing.Process(
                        target=load_experimental_runs,
                        args=(
                            experimental_files_chunk,
                            spectra_and_XIC_queue,
                            file_dir,
                            chemical_pool,
                            mz_unit,
                            mz_tol,
                            chemical_pool_columns,
                        ),
                    )
                    loaders.append(p)
                    p.daemon = True
                    p.start()

            collector = multiprocessing.Process(
                target=collect_experimental_runs,
                args=(
                    spectra_and_XIC_queue,
                    spectra_runs,
                    XICs_runs,
                    orig_XICs_runs,
                    experimental_runs,
                ),
            )
            collector.daemon = True
            collector.start()

            for p in loaders:
                p.join()
            collector.join()

            spectra_runs, XICs_runs, orig_XICs_runs = (
                dict(spectra_runs),
                dict(XICs_runs),
                dict(orig_XICs_runs),
            )

            logger.info("Identify peaks...")
            standard_peak_boundaries_and_quantities, runs_peak_boundaries_and_quantities = (
                identify_peaks(
                    chemical_pool,
                    chemical_pool_columns,
                    spectra_standard,
                    XICs_standard,
                    orig_XICs_standard,
                    spectra_runs,
                    XICs_runs,
                    orig_XICs_runs,
                    experimental_runs,
                    slope,
                    intercept,
                    adaptive_rt,
                    merge_order_standard,
                    key_runs_standard,
                    n_candidate_peaks,
                    rt_tol,
                    n_compounds_for_output,
                    min_intensity_standard,
                    rt_left_limit,
                    rt_right_limit,
                    max_standard_peak_width,
                    peak_refinement_index_range,
                    quantification_method,
                )
            )

            quantification_matrix = build_quantification_matrix(
                standard_peak_boundaries_and_quantities,
                runs_peak_boundaries_and_quantities,
                standard_prefix_pool,
                background_intensity,
                delta_rt_cutoff,
                max_peak_width,
                experimental_runs=experimental_runs,
            )

            all_candidate_matrix = quantification_matrix
            selection_annotations = {}
            selection_summary = dict(
                version=__version__,
                rt_left_limit=rt_left_limit,
                rt_right_limit=rt_right_limit,
                mode=peak_selection,
                selected_compounds=len(quantification_matrix),
                ambiguous_compounds=0,
                isotope_veto_candidates=0,
                mass_center_ppm=None,
            )
            if peak_selection == "identity":
                selection = select_reference_peaks(
                    chemical_pool,
                    chemical_pool_columns,
                    spectra_standard,
                    standard_peak_boundaries_and_quantities,
                    identity_evidence,
                    parents,
                    relative_height=reference_relative_height,
                    ambiguity_margin=reference_ambiguity_margin,
                )
                quantification_matrix = all_candidate_matrix.loc[selection["selected"].index]
                audit = library_mass_audit(chemical_pool, chemical_pool_columns, acquisition_mode)
                audit.to_csv(os.path.join(out, f"{pool}_OpenASChem_library_check.csv"), index=False)
                if audit.needs_review.any():
                    logger.warning(
                        "%d library entries have a formula/mass discrepancy; see library_check.csv. Source masses were not changed.",
                        int(audit.needs_review.sum()),
                    )
                for name, key in [
                    ("reference_candidates", "candidates"),
                    ("reference_decisions", "decisions"),
                    ("ambiguous", "ambiguous"),
                    ("isotope_evidence", "isotope_evidence"),
                ]:
                    selection[key].to_csv(
                        os.path.join(out, f"{pool}_OpenASChem_{name}.csv"), index=False
                    )
                for _, row in selection["candidates"].iterrows():
                    selection_annotations[row["key"]] = dict(
                        status=row.selection_status,
                        score=float(row.identity_score)
                        if math.isfinite(row.identity_score)
                        else None,
                        mass_error_ppm=float(row.mass_error_ppm),
                        isotope_loss=float(row.isotope_loss),
                    )
                selection_summary.update(
                    selected_compounds=len(quantification_matrix),
                    ambiguous_compounds=len(selection["ambiguous"]),
                    isotope_veto_candidates=int(selection["candidates"].isotope_veto.sum()),
                    mass_center_ppm=selection["mass_center_ppm"],
                    relative_height_floor=reference_relative_height,
                    ambiguity_margin=reference_ambiguity_margin,
                    parent_library=str(identity_parent_library)
                    if identity_parent_library
                    else "input library",
                    score="-isotope_loss -abs(mass_error-center) +ln(relative_height)",
                )
            logger.info(
                "Reference selection: %s; %d accepted; %d ambiguous; %d isotope-vetoed candidate peaks.",
                peak_selection,
                selection_summary["selected_compounds"],
                selection_summary["ambiguous_compounds"],
                selection_summary["isotope_veto_candidates"],
            )
            with open(
                os.path.join(out, f"{pool}_OpenASChem_selection.json"), "w", encoding="utf-8"
            ) as stream:
                json.dump(selection_summary, stream, indent=2, allow_nan=False)

            fold_changes, p_values = calculate_fold_changes(
                quantification_matrix,
                standard_prefix_pool,
                min_fold_change,
                min_fold_change_single,
                p_value_cutoff,
                p_value_cutoff_single,
                background_intensity,
                epsilon,
            )

            logger.info("Writing output...")
            results_compound = build_hit_table(
                chemical_pool,
                chemical_pool_columns,
                bonus_columns,
                output_columns,
                fold_changes,
                p_values,
                quantification_matrix,
                standard_peak_boundaries_and_quantities,
                runs_peak_boundaries_and_quantities,
                acquisition_mode,
                min_fold_change,
                standard_prefix_pool,
                output_RT_in_minutes,
            )
            results_compound.to_csv(
                os.path.join(out, "%s%s" % (pool, output_results_suffix)), index=False
            )

            results_all = build_quantification_table(
                chemical_pool,
                chemical_pool_columns,
                bonus_columns,
                output_columns,
                quantification_matrix,
                standard_peak_boundaries_and_quantities,
                acquisition_mode,
                standard_prefix_pool,
                output_RT_in_minutes,
            )
            results_all.to_csv(
                os.path.join(out, "%s%s" % (pool, output_matrix_suffix)), index=False
            )

            report = write_interactive_report(
                out=out,
                pool=pool,
                quantification_method=quantification_method,
                standard_prefix=standard_prefix_pool,
                quantification_matrix=quantification_matrix,
                fold_changes=fold_changes,
                results_compound=results_compound,
                output_columns=output_columns,
                results_all=results_all,
                spectra_standard=spectra_standard,
                XICs_standard=XICs_standard,
                orig_XICs_standard=orig_XICs_standard,
                spectra_runs=spectra_runs,
                XICs_runs=XICs_runs,
                orig_XICs_runs=orig_XICs_runs,
                standard_peak_boundaries_and_quantities=standard_peak_boundaries_and_quantities,
                runs_peak_boundaries_and_quantities=runs_peak_boundaries_and_quantities,
                selection_annotations=selection_annotations,
                selection_summary=selection_summary,
            )
            logger.info(
                "Interactive report: %s (%d records, %d candidate hits; %.2f s)",
                report["path"],
                report["records"],
                report["hits"],
                report["seconds"],
            )
            if peak_selection == "identity":
                # Audit every candidate, including excluded and ambiguous references.
                # Only final accepted hits carry hit flags in this companion report.
                review_metadata = build_quantification_table(
                    chemical_pool,
                    chemical_pool_columns,
                    bonus_columns,
                    output_columns,
                    all_candidate_matrix,
                    standard_peak_boundaries_and_quantities,
                    acquisition_mode,
                    standard_prefix_pool,
                    output_RT_in_minutes,
                )
                review = write_interactive_report(
                    out=out,
                    filename=f"{pool}_OpenASChem_review.html",
                    pool=pool,
                    quantification_method=quantification_method,
                    standard_prefix=standard_prefix_pool,
                    quantification_matrix=all_candidate_matrix,
                    fold_changes=fold_changes,
                    results_compound=results_compound,
                    output_columns=output_columns,
                    results_all=review_metadata,
                    spectra_standard=spectra_standard,
                    XICs_standard=XICs_standard,
                    orig_XICs_standard=orig_XICs_standard,
                    spectra_runs=spectra_runs,
                    XICs_runs=XICs_runs,
                    orig_XICs_runs=orig_XICs_runs,
                    standard_peak_boundaries_and_quantities=standard_peak_boundaries_and_quantities,
                    runs_peak_boundaries_and_quantities=runs_peak_boundaries_and_quantities,
                    selection_annotations=selection_annotations,
                    selection_summary=selection_summary,
                    report_kind="review",
                )
                logger.info("All-candidate review report: %s", review["path"])

            completed_pools.append(
                build_pool_summary(
                    pool=pool,
                    results_compound=results_compound,
                    output_columns=output_columns,
                    input_compounds=len(chemical_pool),
                    run_count=len(standard_files) + len(experimental_runs),
                    quantified_compounds=len(quantification_matrix),
                    selection_summary=selection_summary,
                    report_file=f"{pool}_OpenASChem_report.html",
                    review_file=f"{pool}_OpenASChem_review.html"
                    if peak_selection == "identity"
                    else None,
                    output_rt_in_minutes=output_RT_in_minutes,
                )
            )
            overview = write_experiment_summary(
                out,
                completed_pools,
                analysis_status="complete"
                if len(completed_pools) == len(raw_data_info)
                else "in_progress",
            )

    logger.info(
        "Experiment summary: %s (%d pools, %d hit compounds, %d hit-target rows)",
        overview["path"],
        overview["pools"],
        overview["hit_compounds"],
        overview["hit_target_rows"],
    )
    logger.info("Done!")
