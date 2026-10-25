"""Targeted XIC extraction, smoothing and peak quantification."""

import bisect
from typing import List, Tuple

import numpy as np
from scipy.integrate import simpson
from scipy.signal import find_peaks, peak_widths


def extract_xic(
    spectra: List[Tuple[np.ndarray, np.ndarray]], mz_to_extract: float, mz_unit: str, mz_tol: float
) -> List[float]:
    """
    Calculate the Extracted Ion Chromatogram (XIC) for a given m/z value.

    Parameters:
    - spectra (List[Tuple[np.ndarray, np.ndarray]]): List of tuples containing m/z arrays and intensity arrays.
    - mz_to_extract (float): The m/z value to extract.
    - mz_unit (str): The unit for m/z tolerance ('Da' or 'ppm').
    - mz_tol (float): The tolerance for m/z extraction.

    Returns:
    - List[float]: The XIC intensity values for each spectrum.

    Example usage:
    >>> spectra = [(np.array([100, 150, 200, 250, 300]), np.array([10, 20, 30, 40, 50]))]
    >>> extract_xic(spectra, 200, "Da", 10)
    [30]
    """
    if mz_to_extract == -1:
        return [0.0 for _ in spectra]

    if mz_unit == "Da":
        extract_width = [mz_to_extract - mz_tol / 2, mz_to_extract + mz_tol / 2]
    elif mz_unit == "ppm":
        extract_width = [mz_to_extract / (1 + mz_tol / 10**6), mz_to_extract * (1 + mz_tol / 10**6)]
    else:
        raise ValueError("Invalid mz_unit. Use 'Da' or 'ppm'.")

    xic = []
    for mz_array, intensity_array in spectra:
        start_idx = bisect.bisect_left(mz_array, extract_width[0])
        end_idx = bisect.bisect_right(mz_array, extract_width[1])
        xic_value = sum(intensity_array[start_idx:end_idx])
        xic.append(xic_value)

    return xic


def smooth_array(a_array):
    if a_array.shape[1] < 1:
        return a_array

    new_arr = np.zeros_like(a_array)
    for i in range(a_array.shape[1]):
        if i == 0:
            new_arr[:, i] = (2 * a_array[:, i] + a_array[:, i + 1]) / 3
        elif i == a_array.shape[1] - 1:
            new_arr[:, i] = (2 * a_array[:, i] + a_array[:, i - 1]) / 3
        else:
            new_arr[:, i] = 0.5 * a_array[:, i] + 0.25 * (a_array[:, i - 1] + a_array[:, i + 1])
    return new_arr


def calc_area(chrom_inten, rt_list):
    rt_list_diff = np.array(rt_list)[1:] - np.array(rt_list)[:-1]
    trapezoids = 0
    for i in range(len(rt_list_diff)):
        trapezoids += (chrom_inten[i] + chrom_inten[i + 1]) * rt_list_diff[i]
    return trapezoids / 2


def find_best_peak_intensity(chromatogram):
    if sum(chromatogram) == 0:
        return None, None
    # Step 1: Detect peaks
    peaks, properties = find_peaks(chromatogram, height=0, prominence=1)

    # Step 2: Identify the peak with the highest intensity
    best_peak_index = np.argmax(properties["peak_heights"])
    best_peak = peaks[best_peak_index]

    # Step 3: Calculate the width of the best peak to find boundaries
    widths, width_heights, left_ips, right_ips = peak_widths(
        chromatogram, [best_peak], rel_height=0.5
    )

    # Step 4: Determine the boundaries of the best peak
    left_boundary = int(left_ips[0])
    right_boundary = int(right_ips[0])

    return left_boundary, right_boundary


def find_candidate_peaks(chromatogram, min_relative_prominence=0.10):
    """Detect local maxima with enough prominence relative to their own height.

    The default excludes shallow ripples on a broad peak without imposing a
    global intensity fraction that would discard isolated weak peaks. Set the
    threshold to zero to reproduce the previous candidate detector.
    """
    if not np.isfinite(min_relative_prominence) or not 0 <= min_relative_prominence <= 1:
        raise ValueError("min_relative_prominence must be between 0 and 1.")
    peaks, properties = find_peaks(chromatogram, height=0, prominence=1)
    relative_prominence = properties["prominences"] / properties["peak_heights"]
    keep = relative_prominence >= min_relative_prominence
    properties = {name: values[keep] for name, values in properties.items()}
    properties["relative_prominences"] = relative_prominence[keep]
    return peaks[keep], properties


def find_top_peaks(chromatogram, num_peaks=3, *, min_relative_prominence=0.10):
    """Rank eligible candidates by height and return half-prominence bounds.

    Relative prominence is filtered before selecting the highest ``num_peaks``.
    The analysis pipeline requests one peak; this helper retains its top-3 API.
    """
    peaks, properties = find_candidate_peaks(chromatogram, min_relative_prominence)
    num_peaks = min(num_peaks, len(peaks))
    if num_peaks == 0:
        return []
    sorted_indices = np.argsort(properties["peak_heights"])[-num_peaks:]
    top_peaks = peaks[sorted_indices]
    _, _, left_ips, right_ips = peak_widths(chromatogram, top_peaks, rel_height=0.5)
    return [(int(left), int(right)) for left, right in zip(left_ips, right_ips)]


def quantify(
    chromatogram, orig_chromatogram, left_boundary, right_boundary, method="area", *, rt_list=None
):
    """Quantify a peak; area uses the original run's RT axis (seconds).

    ``rt_list`` is required for area mode. Peak-height mode retains its existing
    positional API and returns the original, unsmoothed maximum intensity.
    """
    if method == "area":
        if rt_list is None:
            raise ValueError("rt_list is required for area quantification.")
        rt_array = np.asarray(rt_list, dtype=float)
        if rt_array.ndim != 1 or len(rt_array) != len(chromatogram):
            raise ValueError("rt_list must be one-dimensional and match the chromatogram length.")
        if not np.all(np.isfinite(rt_array)) or np.any(np.diff(rt_array) <= 0):
            raise ValueError("rt_list must contain finite, strictly increasing retention times.")
        # Step 5: Estimate the baseline (background) around the peak
        # Taking points before the left boundary and after the right boundary for baseline estimation
        baseline_region = np.concatenate(
            (
                chromatogram[max(0, left_boundary - 5) : left_boundary],
                chromatogram[right_boundary + 1 : min(len(chromatogram), right_boundary + 6)],
            )
        )
        baseline_value = np.median(baseline_region)  # Using median for a robust estimate

        # Step 6: Subtract the baseline from the peak region
        peak_region = chromatogram[left_boundary : right_boundary + 1]
        peak_region_bg_subtracted = peak_region - baseline_value
        peak_region_bg_subtracted[peak_region_bg_subtracted < 0] = 0  # Ensure no negative values

        # Step 7: Calculate the quantity (area under the curve) of the background-subtracted peak
        quantity = simpson(
            peak_region_bg_subtracted, x=rt_array[left_boundary : right_boundary + 1]
        )
    elif method == "intensity":
        quantity = max(orig_chromatogram[left_boundary : right_boundary + 1])
    else:
        raise Exception("Wrong quantification method: %s!" % method)
    return quantity


def find_best_peak_intensity_with_bg_subtraction(chromatogram, rt_list):
    # Step 1: Detect peaks
    peaks, properties = find_peaks(chromatogram, height=0, prominence=1)

    # Step 2: Identify the peak with the highest intensity
    best_peak_index = np.argmax(properties["peak_heights"])
    best_peak = peaks[best_peak_index]

    # Step 3: Calculate the width of the best peak to find boundaries
    widths, width_heights, left_ips, right_ips = peak_widths(
        chromatogram, [best_peak], rel_height=0.5
    )

    # Step 4: Determine the boundaries of the best peak
    left_boundary = int(left_ips[0])
    right_boundary = int(right_ips[0])

    quantity = quantify(
        chromatogram, chromatogram, left_boundary, right_boundary, method="area", rt_list=rt_list
    )

    return left_boundary, right_boundary, quantity
