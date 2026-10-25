"""Alignment, quantification and hit-table helpers."""

import os

import numpy as np
import pandas as pd
from scipy import stats

from .chromatogram_processing import extract_xic, find_top_peaks, quantify, smooth_array
from .ms_file_processing import load_rawdata
from .utils import get_openaschem_logger, jitter_array


def extract_pool_xics(all_spectra, chemical_pool, mz_unit, mz_tol, chemical_pool_columns):
    orig_XICs, XICs = {}, {}
    for compound_id, compound_mz in zip(
        chemical_pool[chemical_pool_columns["COMPOUND_ID"]],
        chemical_pool[chemical_pool_columns["COMPOUND_MZ"]],
    ):
        orig_xic = np.array(extract_xic(all_spectra.spectra, compound_mz, mz_unit, mz_tol))
        xic = smooth_array(np.array([orig_xic]))[0]

        orig_XICs[compound_id] = orig_xic
        XICs[compound_id] = xic

    return XICs, orig_XICs


def load_experimental_runs(
    experimental_files_chunk,
    spectra_and_XIC_queue,
    file_dir,
    chemical_pool,
    mz_unit,
    mz_tol,
    chemical_pool_columns,
):
    for experimental_file in experimental_files_chunk:
        spectra_run = load_rawdata(os.path.join(file_dir, experimental_file))
        XICs_run, orig_XICs_run = extract_pool_xics(
            spectra_run, chemical_pool, mz_unit, mz_tol, chemical_pool_columns
        )
        spectra_and_XIC_queue.put([experimental_file, spectra_run, XICs_run, orig_XICs_run])


def collect_experimental_runs(
    spectra_and_XIC_queue, spectra_runs, XICs_runs, orig_XICs_runs, experimental_runs
):
    logger = get_openaschem_logger()
    for i, _ in enumerate(experimental_runs):
        spectra_and_XIC_run = spectra_and_XIC_queue.get()
        spectra_and_XIC_queue.task_done()
        experimental_file, spectra_run, XICs_run, orig_XICs_run = spectra_and_XIC_run
        if experimental_file.endswith(".mzXML"):
            spectra_runs[experimental_file[:-6]] = spectra_run
            XICs_runs[experimental_file[:-6]] = XICs_run
            orig_XICs_runs[experimental_file[:-6]] = orig_XICs_run
        else:
            spectra_runs[experimental_file[:-5]] = spectra_run
            XICs_runs[experimental_file[:-5]] = XICs_run
            orig_XICs_runs[experimental_file[:-5]] = orig_XICs_run

        logger.info("(%d / %d) files loaded..." % (i + 1, len(experimental_runs)))


def do_align_linear(
    rt_list_ref, ms2_inten_lists_ref, rt_list_exp, ms2_inten_lists_exp, global_fit_arguments
):
    """DIAlignR-compatible hybrid alignment using the bundled Python/Numba backend."""
    from .dialignr_numba import get_aligned_times_linear_numba

    slope, intercept, adaptive_rt = global_fit_arguments
    return get_aligned_times_linear_numba(
        rt_list_ref,
        ms2_inten_lists_ref,
        rt_list_exp,
        ms2_inten_lists_exp,
        slope,
        intercept,
        adaptive_rt,
    )


def merge_aligned_times(aligned_times, merge_order, key_runs):
    time_vectors = aligned_times[merge_order[0]]
    time_vectors = time_vectors.fillna(-1)
    for merge_id, key_run in zip(merge_order[1:], key_runs[1:]):
        time_vectors = pd.merge(time_vectors, aligned_times[merge_id], on=key_run, how="left")
        time_vectors = time_vectors.fillna(-1)
    time_vectors[time_vectors == -1] = pd.NA

    return time_vectors


def align_one_compound_standards(
    rt_list_1,
    rt_list_2,
    rt_list_3,
    XIC_1,
    XIC_2,
    XIC_3,
    slope,
    intercept,
    adaptive_rt,
    merge_order,
    key_runs,
):
    aligned_vector12 = do_align_linear(
        rt_list_1, [XIC_1], rt_list_2, [XIC_2], [slope, intercept, adaptive_rt]
    )
    aligned_vector13 = do_align_linear(
        rt_list_1, [XIC_1], rt_list_3, [XIC_3], [slope, intercept, adaptive_rt]
    )
    aligned_vectors = {
        ("S1", "S2"): pd.DataFrame(aligned_vector12, columns=["S1", "S2"]),
        ("S1", "S3"): pd.DataFrame(aligned_vector13, columns=["S1", "S3"]),
    }

    aligned_res_compound = merge_aligned_times(aligned_vectors, merge_order, key_runs)
    return aligned_res_compound


def aligned_time_to_index(aligned_time, middle_rts, rt_tolerance=1.65):
    middle_rt_middle_points = [
        (middle_rts[i] + middle_rts[i + 1]) / 2 for i in range(len(middle_rts) - 1)
    ]
    middle_rt_middle_points = list(np.round(middle_rt_middle_points, 2))

    middle_rts = list(np.round(middle_rts, 2))

    aligned_index = []
    for time_point in aligned_time:
        if pd.isna(time_point):
            aligned_index.append(pd.NA)
        else:
            try:
                aligned_index.append(middle_rts.index(time_point))
            except ValueError:
                try:
                    aligned_index.append(middle_rt_middle_points.index(time_point))
                except ValueError:
                    bias = np.abs(time_point - np.array(middle_rts))
                    if min(bias) > rt_tolerance:
                        aligned_index.append(pd.NA)
                    else:
                        aligned_index.append(np.argmin(bias))
    return aligned_index


def time_vector_to_index_vector(time_vector, original_time_vectors, rt_tolerance):
    index_vector = {}
    for run_name in list(time_vector.columns):
        aligned_time = list(time_vector[run_name])
        middle_rts = original_time_vectors[run_name]
        index_vector[run_name] = aligned_time_to_index(aligned_time, middle_rts, rt_tolerance / 2)

    return pd.DataFrame(index_vector)


def interpolate_index_vector(index_vector):
    interpolated_index_vector = {}
    for run_name in list(index_vector.columns):
        interpolated_index_vector[run_name] = []

        last_not_na = -1
        na_count = 0

        for time_point in index_vector[run_name]:
            if pd.isna(time_point):
                na_count += 1
            else:
                if last_not_na == -1:
                    if time_point == 0:
                        interpolated_index_vector[run_name].extend([pd.NA] * na_count + [0])
                    else:
                        padding = list(range(time_point))
                        if len(padding) <= na_count:
                            interpolated_index_vector[run_name].extend(
                                [pd.NA] * (na_count - len(padding)) + padding
                            )
                        else:
                            interpolated_index_vector[run_name].extend(
                                padding[time_point - na_count :]
                            )
                        interpolated_index_vector[run_name].append(time_point)

                    na_count = 0

                else:
                    if na_count == 0:
                        interpolated_index_vector[run_name].append(time_point)
                    else:
                        padding = [
                            int(ii) for ii in np.linspace(last_not_na, time_point, na_count + 2)[1:]
                        ]
                        interpolated_index_vector[run_name].extend(padding)
                        na_count = 0

                last_not_na = time_point

        if na_count != 0:
            interpolated_index_vector[run_name].extend([pd.NA] * na_count)

    return pd.DataFrame(interpolated_index_vector).dropna().astype(int)


def index_vector_to_score_vector(index_vector, original_intensity_vectors):
    score_vector = {}
    for run_name in list(index_vector.columns):
        score_vector[run_name] = original_intensity_vectors[run_name][index_vector[run_name].values]
    return pd.DataFrame(score_vector)


def identify_peaks(
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
    rt_tolerance,
    n_compounds_for_output,
    min_intensity_standard,
    rt_left_limit,
    rt_right_limit,
    max_standard_peak_width,
    peak_refinement_index_range,
    quantification_method,
):
    logger = get_openaschem_logger()

    standard_peak_boundaries_and_quantities = {}
    runs_peak_boundaries_and_quantities = {}

    for i, compound in enumerate(chemical_pool[chemical_pool_columns["COMPOUND_ID"]]):
        # 1. Align the standard runs and pick candidate reference peaks
        aligned_res_standard = align_one_compound_standards(
            spectra_standard[0].rt_list,
            spectra_standard[1].rt_list,
            spectra_standard[2].rt_list,
            XICs_standard[0][compound],
            XICs_standard[1][compound],
            XICs_standard[2][compound],
            slope,
            intercept,
            adaptive_rt,
            merge_order_standard,
            key_runs_standard,
        )

        index_vector = interpolate_index_vector(
            time_vector_to_index_vector(
                aligned_res_standard,
                {
                    "S1": spectra_standard[0].rt_list,
                    "S2": spectra_standard[1].rt_list,
                    "S3": spectra_standard[2].rt_list,
                },
                rt_tolerance,
            )
        )
        aligned_intensity_vector = index_vector_to_score_vector(
            index_vector,
            {
                "S1": XICs_standard[0][compound],
                "S2": XICs_standard[1][compound],
                "S3": XICs_standard[2][compound],
            },
        )
        cross_run_intensity_vector = aligned_intensity_vector.sum(axis=1)
        peak_boundaries = find_top_peaks(cross_run_intensity_vector.values, n_candidate_peaks)

        if len(peak_boundaries) == 0:
            if i and (i + 1) % n_compounds_for_output == 0:
                logger.info(
                    "Peak identification of (%d / %d) compounds done!"
                    % (i + 1, chemical_pool.shape[0])
                )
            continue

        # 2. Filter the reference peaks of the standard runs to get the valid ones
        ## (1) None of the peaks should have a width greater than the cutoff
        ## (2) The reference peaks should be within a reasonable RT range
        ## (3) The intensity of the reference peaks should be greater than the cutoff
        standard_peak_boundaries_and_quantities[compound] = []
        for left_boundary, right_boundary in peak_boundaries:
            lb1, rb1 = (
                index_vector["S1"].iloc[left_boundary],
                index_vector["S1"].iloc[right_boundary],
            )
            lb2, rb2 = (
                index_vector["S2"].iloc[left_boundary],
                index_vector["S2"].iloc[right_boundary],
            )
            lb3, rb3 = (
                index_vector["S3"].iloc[left_boundary],
                index_vector["S3"].iloc[right_boundary],
            )

            # Standard peak widths are elapsed seconds, not scan-index counts.
            # Apply the limit independently in all three original RT grids.
            standard_widths = (
                spectra_standard[0].rt_list[rb1] - spectra_standard[0].rt_list[lb1],
                spectra_standard[1].rt_list[rb2] - spectra_standard[1].rt_list[lb2],
                spectra_standard[2].rt_list[rb3] - spectra_standard[2].rt_list[lb3],
            )
            if any(width > max_standard_peak_width for width in standard_widths):
                continue

            q1 = quantify(
                XICs_standard[0][compound],
                orig_XICs_standard[0][compound],
                lb1,
                rb1,
                quantification_method,
                rt_list=spectra_standard[0].rt_list,
            )
            q2 = quantify(
                XICs_standard[1][compound],
                orig_XICs_standard[1][compound],
                lb2,
                rb2,
                quantification_method,
                rt_list=spectra_standard[1].rt_list,
            )
            q3 = quantify(
                XICs_standard[2][compound],
                orig_XICs_standard[2][compound],
                lb3,
                rb3,
                quantification_method,
                rt_list=spectra_standard[2].rt_list,
            )

            lb_mean = (
                spectra_standard[0].rt_list[lb1]
                + spectra_standard[1].rt_list[lb2]
                + spectra_standard[2].rt_list[lb3]
            ) / 3
            rb_mean = (
                spectra_standard[0].rt_list[rb1]
                + spectra_standard[1].rt_list[rb2]
                + spectra_standard[2].rt_list[rb3]
            ) / 3

            if (
                q1 > min_intensity_standard
                and q2 > min_intensity_standard
                and q3 > min_intensity_standard
                and lb_mean > rt_left_limit
                and rb_mean < rt_right_limit
            ):
                standard_peak_boundaries_and_quantities[compound].append(
                    {
                        "lb1": spectra_standard[0].rt_list[lb1],
                        "rb1": spectra_standard[0].rt_list[rb1],
                        "q1": q1,
                        "lb2": spectra_standard[1].rt_list[lb2],
                        "rb2": spectra_standard[1].rt_list[rb2],
                        "q2": q2,
                        "lb3": spectra_standard[2].rt_list[lb3],
                        "rb3": spectra_standard[2].rt_list[rb3],
                        "q3": q3,
                    }
                )

        if len(standard_peak_boundaries_and_quantities[compound]) == 0:
            _ = standard_peak_boundaries_and_quantities.pop(compound, None)
            if i and (i + 1) % n_compounds_for_output == 0:
                logger.info(
                    "Peak identification of (%d / %d) compounds done!"
                    % (i + 1, chemical_pool.shape[0])
                )
            continue

        # 3. Find the corresponding peaks in experimental runs
        runs_peak_boundaries_and_quantities[compound] = []
        for candidate_peak in standard_peak_boundaries_and_quantities[compound]:
            lb1, rb1 = candidate_peak["lb1"], candidate_peak["rb1"]

            lb_runs, rb_runs, q_runs = {}, {}, {}
            for run in experimental_runs:
                lb_index = np.argmin(np.abs(np.array(spectra_runs[run].rt_list) - lb1))
                rb_index = np.argmin(np.abs(np.array(spectra_runs[run].rt_list) - rb1))
                lb_index_extended = max(0, lb_index - peak_refinement_index_range)
                rb_index_extended = min(
                    len(spectra_runs[run].rt_list), rb_index + peak_refinement_index_range
                )
                refined_peak = find_top_peaks(
                    XICs_runs[run][compound][lb_index_extended:rb_index_extended], 1
                )
                if len(refined_peak) != 0:
                    lb_index_refined_pseudo, rb_index_refined_pseudo = refined_peak[0]
                else:
                    lb_index_refined_pseudo, rb_index_refined_pseudo = None, None

                if lb_index_refined_pseudo is not None and rb_index_refined_pseudo is not None:
                    lb_index_refined = list(range(lb_index_extended, rb_index_extended + 1))[
                        lb_index_refined_pseudo
                    ]
                    rb_index_refined = list(range(lb_index_extended, rb_index_extended + 1))[
                        rb_index_refined_pseudo
                    ]
                else:
                    lb_index_refined, rb_index_refined = lb_index, rb_index

                q_run = quantify(
                    XICs_runs[run][compound],
                    orig_XICs_runs[run][compound],
                    lb_index_refined,
                    rb_index_refined,
                    quantification_method,
                    rt_list=spectra_runs[run].rt_list,
                )
                q_runs[run] = q_run

                lb_runs[run] = spectra_runs[run].rt_list[lb_index_refined]
                rb_runs[run] = spectra_runs[run].rt_list[rb_index_refined]

            runs_peak_boundaries_and_quantities[compound].append([lb_runs, rb_runs, q_runs])

        if i and (i + 1) % n_compounds_for_output == 0:
            logger.info(
                "Peak identification of (%d / %d) compounds done!" % (i + 1, chemical_pool.shape[0])
            )

    return standard_peak_boundaries_and_quantities, runs_peak_boundaries_and_quantities


def build_quantification_matrix(
    standard_peak_boundaries_and_quantities,
    runs_peak_boundaries_and_quantities,
    standard_prefix,
    background_intensity,
    delta_rt_cutoff,
    max_peak_width,
    experimental_runs=None,
):
    matrix = {
        "%s_1" % standard_prefix: [],
        "%s_2" % standard_prefix: [],
        "%s_3" % standard_prefix: [],
    }
    if experimental_runs is None:
        experimental_runs = (
            list(next(iter(runs_peak_boundaries_and_quantities.values()))[0][2])
            if runs_peak_boundaries_and_quantities
            else []
        )
    for run in experimental_runs:
        matrix[run] = []
    valid_compounds = []

    for compound in standard_peak_boundaries_and_quantities:
        for j, candidate in enumerate(standard_peak_boundaries_and_quantities[compound]):
            valid_compounds.append(compound + "_%d" % (j + 1))
            matrix["%s_1" % standard_prefix].append(candidate["q1"])
            matrix["%s_2" % standard_prefix].append(candidate["q2"])
            matrix["%s_3" % standard_prefix].append(candidate["q3"])

            for run in experimental_runs:
                standard_apex_rt = (candidate["lb1"] + candidate["rb1"]) / 2
                run_apex_rt = (
                    runs_peak_boundaries_and_quantities[compound][j][0][run]
                    + runs_peak_boundaries_and_quantities[compound][j][1][run]
                ) / 2
                peak_width = (
                    runs_peak_boundaries_and_quantities[compound][j][1][run]
                    - runs_peak_boundaries_and_quantities[compound][j][0][run]
                )

                if (
                    abs(standard_apex_rt - run_apex_rt) < delta_rt_cutoff
                    and peak_width < max_peak_width
                ):
                    matrix[run].append(runs_peak_boundaries_and_quantities[compound][j][2][run])
                else:
                    matrix[run].append(background_intensity)

    matrix = pd.DataFrame(matrix, index=valid_compounds)
    matrix[matrix <= background_intensity] = background_intensity
    return matrix


def calculate_fold_changes(
    quantification_matrix,
    standard_prefix,
    min_fold_change,
    min_fold_change_single,
    p_value_cutoff,
    p_value_cutoff_single,
    background_intensity,
    epsilon,
):
    filtered_quantification_matrix = quantification_matrix.copy()

    # 1. Calculate fold changes
    experimental_runs = [
        run
        for run in list(filtered_quantification_matrix.columns)
        if run not in ["%s_1" % standard_prefix, "%s_2" % standard_prefix, "%s_3" % standard_prefix]
    ]
    experimental_samples = list(set([run[:-2] for run in experimental_runs]))

    fold_changes = []
    for sample in experimental_samples:
        target_runs = [sample + "_1", sample + "_2", sample + "_3"]
        other_runs = [run for run in experimental_runs if run not in target_runs]
        target_quant = filtered_quantification_matrix.loc[:, target_runs].mean(axis=1)
        other_quant = filtered_quantification_matrix.loc[:, other_runs].mean(axis=1)
        fold_change = target_quant / (other_quant + epsilon)
        fold_changes.append(fold_change)
    fold_changes = pd.DataFrame(fold_changes).T
    fold_changes.columns = experimental_samples

    filtered_fold_changes = fold_changes[fold_changes.gt(min_fold_change).any(axis=1)]

    # 2. Filter candidates
    valid_compound_candidates, valid_compound_p_values = [], []
    for compound_candidate in list(filtered_fold_changes.index):
        valid_flag = 1

        fold_change_part = filtered_fold_changes.loc[compound_candidate, :]

        for sample, fold_change_value in zip(list(fold_change_part.index), fold_change_part.values):
            if fold_change_value > min_fold_change:
                # All the replicates should have positive signals.
                target_runs = [sample + "_1", sample + "_2", sample + "_3"]
                target_quantifications = filtered_quantification_matrix.loc[
                    compound_candidate, target_runs
                ]
                if (target_quantifications > background_intensity).sum() < 3:
                    valid_flag = 0
                    break

                # One vs All statistical tests
                other_runs = [run for run in experimental_runs if run not in target_runs]
                other_quantifications = filtered_quantification_matrix.loc[
                    compound_candidate, other_runs
                ]
                t_stat, p_value = stats.ttest_ind(
                    jitter_array(target_quantifications), jitter_array(other_quantifications)
                )
                if t_stat > 0:
                    p_value = p_value / 2
                else:
                    p_value = 1
                if p_value > p_value_cutoff:
                    valid_flag = 0
                    break

                # Each one vs one fold change needs to be greater than a cutoff value, and each one vs one t-test has to be confident.
                other_samples = [i for i in experimental_samples if i != sample]
                for control_sample in other_samples:
                    control_runs = [
                        control_sample + "_1",
                        control_sample + "_2",
                        control_sample + "_3",
                    ]
                    control_quantifications = filtered_quantification_matrix.loc[
                        compound_candidate, control_runs
                    ]
                    control_mean = control_quantifications.mean() + epsilon
                    control_std = control_quantifications.std()

                    fold_change_single_sample = target_quantifications.mean() / control_mean
                    if control_std == 0:
                        t_stat_single_sample, p_value_single_sample = 1, 0
                    else:
                        t_stat_single_sample, p_value_single_sample = stats.ttest_ind(
                            target_quantifications, control_quantifications
                        )
                    if t_stat_single_sample > 0:
                        p_value_single_sample = p_value_single_sample / 2
                    else:
                        p_value_single_sample = 1

                    if (
                        fold_change_single_sample < min_fold_change_single
                        or p_value_single_sample > p_value_cutoff_single
                    ):
                        valid_flag = 0
                        break
                break

        if valid_flag == 1:
            valid_compound_candidates.append(compound_candidate)
            valid_compound_p_values.append(p_value)

    return filtered_fold_changes.loc[valid_compound_candidates, :], valid_compound_p_values


def build_hit_table(
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
    standard_prefix,
    output_RT_in_minutes,
):
    # initiate output
    experimental_runs = [
        name
        for name in quantification_matrix.columns
        if name not in [f"{standard_prefix}_{i}" for i in (1, 2, 3)]
    ]
    output_df = {run: [] for run in experimental_runs}
    other_columns = [
        output_columns["COMPOUND_MZ"],
        output_columns["RT"],
        "%s_1" % standard_prefix,
        "%s_2" % standard_prefix,
        "%s_3" % standard_prefix,
        "fold_change",
        output_columns["POOL_NAME"],
        output_columns["TARGET_PROTEIN"],
        output_columns["SMILES"],
        output_columns["ADDUCT_TYPE"],
        output_columns["FORMULA"],
        output_columns["COMPOUND_ID"],
        output_columns["POOL_ID"],
    ]
    for column in other_columns:
        output_df[column] = []

    for column in bonus_columns:
        output_df[column + "_LIBRARY"] = []

    output_df["p_value"] = p_values

    compounds = list(fold_changes.index)
    for candidate_compound in compounds:
        compound, rank_text = candidate_compound.rsplit("_", 1)
        candidate_id = int(rank_text) - 1

        compound_params = chemical_pool[
            chemical_pool[chemical_pool_columns["COMPOUND_ID"]] == compound
        ]
        output_df[output_columns["COMPOUND_ID"]].append(compound)
        output_df[output_columns["COMPOUND_MZ"]].append(
            compound_params.loc[:, chemical_pool_columns["COMPOUND_MZ"]].iloc[0]
        )
        output_df[output_columns["POOL_NAME"]].append(
            compound_params.loc[:, chemical_pool_columns["POOL_NAME"]].iloc[0]
        )
        output_df[output_columns["SMILES"]].append(
            compound_params.loc[:, chemical_pool_columns["SMILES"]].iloc[0]
        )
        output_df[output_columns["FORMULA"]].append(
            compound_params.loc[:, chemical_pool_columns["FORMULA"]].iloc[0]
        )
        output_df[output_columns["POOL_ID"]].append(
            compound_params.loc[:, chemical_pool_columns["POOL_ID"]].iloc[0]
        )
        if acquisition_mode == "positive":
            output_df[output_columns["ADDUCT_TYPE"]].append("[M+H]+")
        else:
            output_df[output_columns["ADDUCT_TYPE"]].append("[M-H]-")

        for column in bonus_columns:
            output_df[column + "_LIBRARY"].append(compound_params.loc[:, column].iloc[0])

        quantification_record = quantification_matrix.loc[candidate_compound, :]
        for run, quant_value in zip(list(quantification_record.index), quantification_record):
            output_df[run].append(quant_value)

        fold_change_record = fold_changes.loc[candidate_compound, :]
        valid_fold_change = fold_change_record[fold_change_record > min_fold_change]
        target_run = ";".join(list(valid_fold_change.index))
        fold_change_value = ";".join(["%.12g" % i for i in valid_fold_change])
        output_df[output_columns["TARGET_PROTEIN"]].append(target_run)
        output_df["fold_change"].append(fold_change_value)

        standard_data = standard_peak_boundaries_and_quantities[compound][candidate_id]
        output_df[output_columns["RT"]].append((standard_data["lb1"] + standard_data["rb1"]) / 2)

    output_df = pd.DataFrame(output_df)
    output_df["ReferencePeakID"] = compounds
    output_df = output_df.sort_values(by=output_columns["TARGET_PROTEIN"])
    column_order = (
        [output_columns["COMPOUND_MZ"], output_columns["RT"]]
        + experimental_runs
        + [
            "%s_1" % standard_prefix,
            "%s_2" % standard_prefix,
            "%s_3" % standard_prefix,
            "fold_change",
            "p_value",
            output_columns["POOL_NAME"],
            output_columns["TARGET_PROTEIN"],
            output_columns["SMILES"],
            output_columns["ADDUCT_TYPE"],
            output_columns["FORMULA"],
            output_columns["COMPOUND_ID"],
            output_columns["POOL_ID"],
        ]
        + [column + "_LIBRARY" for column in bonus_columns]
        + ["ReferencePeakID"]
    )

    output_df = output_df.loc[:, column_order]

    if output_RT_in_minutes:
        output_df[output_columns["RT"]] /= 60

    return output_df


def build_quantification_table(
    chemical_pool,
    chemical_pool_columns,
    bonus_columns,
    output_columns,
    quantification_matrix,
    standard_peak_boundaries_and_quantities,
    acquisition_mode,
    standard_prefix,
    output_RT_in_minutes,
):
    output_df = quantification_matrix.copy()
    experimental_runs = [
        i
        for i in list(output_df.columns)
        if i not in ["%s_1" % standard_prefix, "%s_2" % standard_prefix, "%s_3" % standard_prefix]
    ]

    all_candidate_compounds = list(quantification_matrix.index)
    output_df[output_columns["COMPOUND_ID"]] = [
        x.rsplit("_", 1)[0] for x in all_candidate_compounds
    ]
    output_df["candidate_id"] = [int(x.rsplit("_", 1)[1]) - 1 for x in all_candidate_compounds]
    output_df["ReferencePeakID"] = all_candidate_compounds
    output_df["standard_total_intensity"] = (
        output_df["%s_1" % standard_prefix]
        + output_df["%s_2" % standard_prefix]
        + output_df["%s_3" % standard_prefix]
    )
    output_df = output_df.sort_values(
        by=[output_columns["COMPOUND_ID"], "standard_total_intensity"], ascending=[True, False]
    )
    output_df = output_df.drop_duplicates(subset=output_columns["COMPOUND_ID"], keep="first")

    static_df = {
        output_columns["COMPOUND_MZ"]: [],
        output_columns["RT"]: [],
        output_columns["POOL_NAME"]: [],
        output_columns["SMILES"]: [],
        output_columns["ADDUCT_TYPE"]: [],
        output_columns["FORMULA"]: [],
        output_columns["POOL_ID"]: [],
    }
    for column in bonus_columns:
        static_df[column + "_LIBRARY"] = []

    for compound, candidate_id in zip(
        output_df[output_columns["COMPOUND_ID"]], output_df["candidate_id"]
    ):
        compound_params = chemical_pool[
            chemical_pool[chemical_pool_columns["COMPOUND_ID"]] == compound
        ]
        static_df[output_columns["COMPOUND_MZ"]].append(
            compound_params.loc[:, chemical_pool_columns["COMPOUND_MZ"]].iloc[0]
        )
        static_df[output_columns["POOL_NAME"]].append(
            compound_params.loc[:, chemical_pool_columns["POOL_NAME"]].iloc[0]
        )
        static_df[output_columns["SMILES"]].append(
            compound_params.loc[:, chemical_pool_columns["SMILES"]].iloc[0]
        )
        static_df[output_columns["FORMULA"]].append(
            compound_params.loc[:, chemical_pool_columns["FORMULA"]].iloc[0]
        )
        static_df[output_columns["POOL_ID"]].append(
            compound_params.loc[:, chemical_pool_columns["POOL_ID"]].iloc[0]
        )
        if acquisition_mode == "positive":
            static_df[output_columns["ADDUCT_TYPE"]].append("[M+H]+")
        else:
            static_df[output_columns["ADDUCT_TYPE"]].append("[M-H]-")

        for column in bonus_columns:
            static_df[column + "_LIBRARY"].append(compound_params.loc[:, column].iloc[0])

        standard_data = standard_peak_boundaries_and_quantities[compound][candidate_id]
        static_df[output_columns["RT"]].append((standard_data["lb1"] + standard_data["rb1"]) / 2)

    static_df = pd.DataFrame(static_df)
    output_df.reset_index(inplace=True, drop=True)
    output_df.drop(["candidate_id", "standard_total_intensity"], axis=1, inplace=True)

    output_df = pd.concat([output_df, static_df], axis=1)

    column_order = (
        [output_columns["COMPOUND_MZ"], output_columns["RT"]]
        + experimental_runs
        + [
            "%s_1" % standard_prefix,
            "%s_2" % standard_prefix,
            "%s_3" % standard_prefix,
            output_columns["POOL_NAME"],
            output_columns["SMILES"],
            output_columns["ADDUCT_TYPE"],
            output_columns["FORMULA"],
            output_columns["COMPOUND_ID"],
            output_columns["POOL_ID"],
        ]
        + [column + "_LIBRARY" for column in bonus_columns]
        + ["ReferencePeakID"]
    )
    output_df = output_df.loc[:, column_order]

    if output_RT_in_minutes:
        output_df[output_columns["RT"]] /= 60

    return output_df
