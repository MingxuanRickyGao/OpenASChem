# SPDX-License-Identifier: Apache-2.0
# Source: DreamDIA-v5/DIAlignR.py (local development snapshot).
# See ALIGNMENT_BACKEND.md for the supported subset and source checksum.
"""
Numba implementation of the DIAlignR getAlignedTimesCpp subset used by DreamDIA.

This module intentionally implements only the current DreamDIA path:
hybrid alignment, mean normalization, dotProductMasked similarity, soft global
RT constraint, affine gap penalties, and overlap alignment.
"""

from __future__ import annotations

import math
from typing import Mapping

import numpy as np

try:
    from numba import njit
except ImportError:  # pragma: no cover - exercised only in environments without numba.
    njit = None


DEFAULT_KERNEL_LEN = 11
DEFAULT_POLY_ORD = 4
DEFAULT_ADAPTIVE_RT = 0.0
DEFAULT_GO_FACTOR = 0.125
DEFAULT_GE_FACTOR = 40.0
DEFAULT_COS_ANGLE_THRESH = 0.3
DEFAULT_DOT_PROD_THRESH = 0.96
DEFAULT_GAP_QUANTILE = 0.5
DEFAULT_SAMPLES_FOR_GRADIENT = 1.0

_SS = 0
_DM = 1
_DA = 2
_DB = 3
_TM = 4
_TA = 5
_TB = 6
_LM = 7
_LA = 8
_LB = 9

_MAT_M = 0
_MAT_A = 1
_MAT_B = 2

_SG_COEFFS = {
    (11, 4): np.array(
        [
            0.041958,
            -0.0699301,
            -0.034965,
            0.034965,
            0.0699301,
            0.041958,
            -0.034965,
            -0.104895,
            -0.0699301,
            0.20979,
            0.916084,
            -0.0699301,
            0.104895,
            0.0699301,
            -0.034965,
            -0.111888,
            -0.104895,
            0.0,
            0.174825,
            0.34965,
            0.412587,
            0.20979,
            -0.034965,
            0.0699301,
            0.02331,
            -0.0536131,
            -0.0815851,
            -0.02331,
            0.11655,
            0.291375,
            0.412587,
            0.34965,
            -0.0699301,
            0.034965,
            -0.034965,
            -0.0536131,
            -0.02331,
            0.04662,
            0.13986,
            0.2331,
            0.296037,
            0.291375,
            0.174825,
            -0.104895,
            0.0699301,
            -0.111888,
            -0.0815851,
            0.04662,
            0.18648,
            0.27972,
            0.296037,
            0.2331,
            0.11655,
            0.0,
            -0.034965,
            0.041958,
            -0.104895,
            -0.02331,
            0.13986,
            0.27972,
            0.333333,
            0.27972,
            0.13986,
            -0.02331,
            -0.104895,
            0.041958,
        ],
        dtype=np.float64,
    ),
    (11, 3): np.array(
        [
            -0.0839161,
            0.0559441,
            0.0909091,
            0.0559441,
            -0.013986,
            -0.0839161,
            -0.118881,
            -0.0839161,
            0.0559441,
            0.335664,
            0.79021,
            0.0559441,
            -0.020979,
            -0.0559441,
            -0.0559441,
            -0.027972,
            0.020979,
            0.0839161,
            0.153846,
            0.223776,
            0.286713,
            0.335664,
            0.0909091,
            -0.0559441,
            -0.102564,
            -0.0745921,
            0.002331,
            0.102564,
            0.200466,
            0.270396,
            0.286713,
            0.223776,
            0.0559441,
            0.0559441,
            -0.0559441,
            -0.0745921,
            -0.0268065,
            0.0606061,
            0.160839,
            0.247086,
            0.292541,
            0.270396,
            0.153846,
            -0.0839161,
            -0.013986,
            -0.027972,
            0.002331,
            0.0606061,
            0.130536,
            0.195804,
            0.240093,
            0.247086,
            0.200466,
            0.0839161,
            -0.118881,
            -0.0839161,
            0.020979,
            0.102564,
            0.160839,
            0.195804,
            0.207459,
            0.195804,
            0.160839,
            0.102564,
            0.020979,
            -0.0839161,
        ],
        dtype=np.float64,
    ),
    (9, 4): np.array(
        [
            0.027195,
            -0.0660451,
            0.003885,
            0.0629371,
            0.034965,
            -0.0582751,
            -0.0971251,
            0.135975,
            0.956488,
            -0.0660451,
            0.149573,
            0.0143745,
            -0.13986,
            -0.128205,
            0.0874126,
            0.398213,
            0.548563,
            0.135975,
            0.003885,
            0.0143745,
            -0.042735,
            -0.0407925,
            0.0699301,
            0.262238,
            0.432012,
            0.398213,
            -0.0971251,
            0.0629371,
            -0.13986,
            -0.0407925,
            0.157343,
            0.314685,
            0.354312,
            0.262238,
            0.0874126,
            -0.0582751,
            0.034965,
            -0.128205,
            0.0699301,
            0.314685,
            0.417249,
            0.314685,
            0.0699301,
            -0.128205,
            0.03496,
        ],
        dtype=np.float64,
    ),
    (9, 3): np.array(
        [
            -0.0707071,
            0.0808081,
            0.0808081,
            0.0,
            -0.0909091,
            -0.121212,
            -0.020202,
            0.282828,
            0.858586,
            0.0808081,
            -0.0707071,
            -0.10101,
            -0.0454545,
            0.0606061,
            0.181818,
            0.282828,
            0.328283,
            0.282828,
            0.0808081,
            -0.10101,
            -0.103175,
            0.00865801,
            0.168831,
            0.311688,
            0.371573,
            0.282828,
            -0.020202,
            0.0,
            -0.0454545,
            0.00865801,
            0.116883,
            0.233766,
            0.313853,
            0.311688,
            0.181818,
            -0.121212,
            -0.0909091,
            0.0606061,
            0.168831,
            0.233766,
            0.255411,
            0.233766,
            0.168831,
            0.0606061,
            -0.0909091,
        ],
        dtype=np.float64,
    ),
    (7, 3): np.array(
        [
            -0.047619,
            0.0952381,
            0.0238095,
            -0.0952381,
            -0.0952381,
            0.190476,
            0.928571,
            0.0952381,
            -0.166667,
            -0.0952381,
            0.142857,
            0.380952,
            0.452381,
            0.190476,
            0.0238095,
            -0.0952381,
            0.047619,
            0.285714,
            0.452381,
            0.380952,
            -0.0952381,
            -0.0952381,
            0.142857,
            0.285714,
            0.333333,
            0.285714,
            0.142857,
            -0.0952381,
        ],
        dtype=np.float64,
    ),
    (7, 2): np.array(
        [
            0.119048,
            -0.0714286,
            -0.142857,
            -0.0952381,
            0.0714286,
            0.357143,
            0.761905,
            -0.0714286,
            0.0,
            0.0714286,
            0.142857,
            0.214286,
            0.285714,
            0.357143,
            -0.142857,
            0.0714286,
            0.214286,
            0.285714,
            0.285714,
            0.214286,
            0.0714286,
            -0.0952381,
            0.142857,
            0.285714,
            0.333333,
            0.285714,
            0.142857,
            -0.0952381,
        ],
        dtype=np.float64,
    ),
}


def _require_numba() -> None:
    if njit is None:
        raise ImportError("numba is required for DIAlignR.py numba alignment.")


def _as_fragment_time_array(values) -> np.ndarray:
    arr = np.asarray(values, dtype=np.float64)
    if arr.ndim != 2:
        raise ValueError("MS2 intensities must be a 2D array/list with shape (n_fragments, n_rt).")
    return np.ascontiguousarray(arr)


def _get_savgol_coeffs(kernel_len: int, poly_ord: int) -> np.ndarray:
    if kernel_len == 0:
        return np.empty(0, dtype=np.float64)
    try:
        return _SG_COEFFS[(kernel_len, poly_ord)]
    except KeyError as exc:
        raise ValueError(
            f"No embedded DIAlignR Savitzky-Golay coefficients for kernel_len={kernel_len}, "
            f"poly_ord={poly_ord}."
        ) from exc


def get_aligned_times_linear_numba(
    rt_ref,
    ms2_ref,
    rt_exp,
    ms2_exp,
    slope: float,
    intercept: float,
    adaptive_rt: float,
    *,
    kernel_len: int = DEFAULT_KERNEL_LEN,
    poly_ord: int = DEFAULT_POLY_ORD,
) -> np.ndarray:
    """Align two XIC groups with a linear global RT constraint."""

    rt_ref_arr = np.asarray(rt_ref, dtype=np.float64)
    bp = rt_ref_arr * float(slope) + float(intercept)
    return get_aligned_times_cpp_numba(
        rt_ref_arr,
        ms2_ref,
        rt_exp,
        ms2_exp,
        bp,
        adaptive_rt,
        kernel_len=kernel_len,
        poly_ord=poly_ord,
    )


def get_aligned_times_lowess_numba(
    rt_ref,
    ms2_ref,
    rt_exp,
    ms2_exp,
    lowess_params,
    adaptive_rt: float,
    *,
    kernel_len: int = DEFAULT_KERNEL_LEN,
    poly_ord: int = DEFAULT_POLY_ORD,
) -> np.ndarray:
    """Align two XIC groups with a lowess/global interpolation constraint."""

    rt_ref_arr = np.asarray(rt_ref, dtype=np.float64)
    rt_exp_arr = np.asarray(rt_exp, dtype=np.float64)
    x, y = _extract_lowess_xy(lowess_params)
    bp = _interp_approxfun_ties_mean(rt_ref_arr, x, y)
    if np.any(~np.isfinite(bp)) or np.any(bp <= 0):
        bp = np.linspace(rt_exp_arr[0], rt_exp_arr[-1], bp.size)
    return get_aligned_times_cpp_numba(
        rt_ref_arr,
        ms2_ref,
        rt_exp_arr,
        ms2_exp,
        bp,
        adaptive_rt,
        kernel_len=kernel_len,
        poly_ord=poly_ord,
    )


def get_aligned_times_cpp_numba(
    rt_ref,
    ms2_ref,
    rt_exp,
    ms2_exp,
    bp,
    adaptive_rt: float,
    *,
    kernel_len: int = DEFAULT_KERNEL_LEN,
    poly_ord: int = DEFAULT_POLY_ORD,
    go_factor: float = DEFAULT_GO_FACTOR,
    ge_factor: float = DEFAULT_GE_FACTOR,
    cos_angle_thresh: float = DEFAULT_COS_ANGLE_THRESH,
    dot_prod_thresh: float = DEFAULT_DOT_PROD_THRESH,
    gap_quantile: float = DEFAULT_GAP_QUANTILE,
    samples_for_gradient: float = DEFAULT_SAMPLES_FOR_GRADIENT,
) -> np.ndarray:
    """Numba equivalent of the DreamDIA-used DIAlignR::getAlignedTimesCpp path.

    Parameters are intentionally limited to the current DreamDIA defaults:
    alignType="hybrid", normalization="mean", simType="dotProductMasked",
    OverlapAlignment=True, hardConstrain=False.
    """

    _require_numba()
    rt_ref_arr = np.ascontiguousarray(np.asarray(rt_ref, dtype=np.float64))
    rt_exp_arr = np.ascontiguousarray(np.asarray(rt_exp, dtype=np.float64))
    ms2_ref_arr = _as_fragment_time_array(ms2_ref)
    ms2_exp_arr = _as_fragment_time_array(ms2_exp)
    bp_arr = np.ascontiguousarray(np.asarray(bp, dtype=np.float64))

    if rt_ref_arr.ndim != 1 or rt_exp_arr.ndim != 1 or bp_arr.ndim != 1:
        raise ValueError("rt_ref, rt_exp, and bp must be 1D arrays.")
    if ms2_ref_arr.shape[1] != rt_ref_arr.size:
        raise ValueError("ms2_ref shape must be (n_fragments, len(rt_ref)).")
    if ms2_exp_arr.shape[1] != rt_exp_arr.size:
        raise ValueError("ms2_exp shape must be (n_fragments, len(rt_exp)).")
    if ms2_ref_arr.shape[0] != ms2_exp_arr.shape[0]:
        raise ValueError("ms2_ref and ms2_exp must contain the same number of fragments.")
    if bp_arr.size != rt_ref_arr.size:
        raise ValueError("bp must have the same length as rt_ref.")

    coeffs = _get_savgol_coeffs(kernel_len, poly_ord)
    return _get_aligned_times_cpp_core(
        rt_ref_arr,
        ms2_ref_arr,
        rt_exp_arr,
        ms2_exp_arr,
        bp_arr,
        float(adaptive_rt),
        int(kernel_len),
        coeffs,
        float(go_factor),
        float(ge_factor),
        float(cos_angle_thresh),
        float(dot_prod_thresh),
        float(gap_quantile),
        float(samples_for_gradient),
    )


def _extract_lowess_xy(lowess_params) -> tuple[np.ndarray, np.ndarray]:
    if isinstance(lowess_params, Mapping):
        x = lowess_params["x"]
        y = lowess_params["y"]
    else:
        x = lowess_params["x"]
        y = lowess_params["y"]
    return np.asarray(x, dtype=np.float64), np.asarray(y, dtype=np.float64)


def _interp_extrapolate(x_new: np.ndarray, x: np.ndarray, y: np.ndarray) -> np.ndarray:
    order = np.argsort(x)
    x = x[order]
    y = y[order]
    out = np.interp(x_new, x, y)
    if x.size < 2:
        return out

    left = x_new < x[0]
    if np.any(left):
        slope = (y[1] - y[0]) / (x[1] - x[0])
        out[left] = y[0] + (x_new[left] - x[0]) * slope

    right = x_new > x[-1]
    if np.any(right):
        slope = (y[-1] - y[-2]) / (x[-1] - x[-2])
        out[right] = y[-1] + (x_new[right] - x[-1]) * slope
    return out


def _interp_approxfun_ties_mean(x_new: np.ndarray, x: np.ndarray, y: np.ndarray) -> np.ndarray:
    order = np.argsort(x)
    x = x[order]
    y = y[order]
    if x.size < 2:
        return np.full(x_new.shape, np.nan, dtype=np.float64)

    unique_x = []
    unique_y = []
    start = 0
    while start < x.size:
        end = start + 1
        while end < x.size and x[end] == x[start]:
            end += 1
        unique_x.append(x[start])
        unique_y.append(float(np.mean(y[start:end])))
        start = end

    x = np.asarray(unique_x, dtype=np.float64)
    y = np.asarray(unique_y, dtype=np.float64)
    if x.size < 2:
        return np.full(x_new.shape, np.nan, dtype=np.float64)

    out = np.interp(x_new, x, y)
    out[(x_new < x[0]) | (x_new > x[-1])] = np.nan
    return out


if njit is not None:

    @njit(cache=True)
    def _quantile_type7(values, p):
        vec = np.sort(values.copy())
        n = vec.size
        m = 1.0 - p
        h = n * p + m
        j = int(math.floor(h))
        g = h - j
        if j <= 0:
            return vec[0]
        if j >= n:
            return vec[n - 1]
        return (1.0 - g) * vec[j - 1] + g * vec[j]

    @njit(cache=True)
    def _smooth_trace(trace, coeffs, frame_size):
        n = trace.size
        output = trace.copy()
        if frame_size == 0 or frame_size > n:
            return output

        mid = frame_size // 2
        current = 0
        for i in range(mid + 1):
            start = current - i
            value = 0.0
            base = (i + 1) * frame_size - 1
            for j in range(frame_size):
                value += trace[start + j] * coeffs[base - j]
            if value < 0.0:
                value = 0.0
            output[current] = value
            current += 1

        steady_end = n - mid
        while current != steady_end:
            start = current - mid
            value = 0.0
            base = mid * frame_size
            for j in range(frame_size):
                value += trace[start + j] * coeffs[base + j]
            if value < 0.0:
                value = 0.0
            output[current] = value
            current += 1

        for i in range(mid - 1, -1, -1):
            start = current - (frame_size - i - 1)
            value = 0.0
            base = i * frame_size
            for j in range(frame_size):
                value += trace[start + j] * coeffs[base + j]
            if value < 0.0:
                value = 0.0
            output[current] = value
            current += 1
        return output

    @njit(cache=True)
    def _smooth_traces(traces, coeffs, frame_size):
        if frame_size == 0:
            return traces.copy()
        out = np.empty_like(traces)
        for frag_idx in range(traces.shape[0]):
            out[frag_idx, :] = _smooth_trace(traces[frag_idx, :], coeffs, frame_size)
        return out

    @njit(cache=True)
    def _mean_normalize(traces):
        total = 0.0
        for i in range(traces.shape[0]):
            for j in range(traces.shape[1]):
                total += traces[i, j]
        mean_value = total / (traces.shape[0] * traces.shape[1])
        out = np.empty_like(traces)
        denom = mean_value + 1e-8
        for i in range(traces.shape[0]):
            for j in range(traces.shape[1]):
                out[i, j] = traces[i, j] / denom
        return out

    @njit(cache=True)
    def _build_dot_product_masked_similarity(ms2_ref, ms2_exp, cos_angle_thresh, dot_prod_thresh):
        ref = _mean_normalize(ms2_ref)
        exp = _mean_normalize(ms2_exp)
        n_frag = ref.shape[0]
        n_ref = ref.shape[1]
        n_exp = exp.shape[1]

        sim = np.zeros((n_ref, n_exp), dtype=np.float64)
        ref_mag = np.zeros(n_ref, dtype=np.float64)
        exp_mag = np.zeros(n_exp, dtype=np.float64)
        for i in range(n_ref):
            value = 0.0
            for frag_idx in range(n_frag):
                value += ref[frag_idx, i] * ref[frag_idx, i]
            ref_mag[i] = math.sqrt(value)
        for j in range(n_exp):
            value = 0.0
            for frag_idx in range(n_frag):
                value += exp[frag_idx, j] * exp[frag_idx, j]
            exp_mag[j] = math.sqrt(value)

        cosine2 = np.zeros((n_ref, n_exp), dtype=np.float64)
        for i in range(n_ref):
            for j in range(n_exp):
                dot = 0.0
                for frag_idx in range(n_frag):
                    dot += ref[frag_idx, i] * exp[frag_idx, j]
                sim[i, j] = dot
                cos_value = dot / ((ref_mag[i] + 1e-6) * (exp_mag[j] + 1e-6))
                if cos_value > 1.0:
                    cos_value = 1.0
                elif cos_value < -1.0:
                    cos_value = -1.0

                abs_cos = abs(cos_value)
                if abs_cos < 1.57:
                    lookup_idx = int(math.floor(abs_cos * 100.0))
                    if lookup_idx < 0:
                        lookup_idx = 0
                    elif lookup_idx > 156:
                        lookup_idx = 156
                    cosine2[i, j] = math.cos(2.0 * math.acos(lookup_idx / 100.0))
                else:
                    cosine2[i, j] = -1.0

        quant = _quantile_type7(sim.ravel(), dot_prod_thresh)
        for i in range(n_ref):
            for j in range(n_exp):
                if sim[i, j] < quant:
                    allow_value = 1.0
                else:
                    allow_value = cosine2[i, j]
                if allow_value <= cos_angle_thresh:
                    sim[i, j] = 0.0
        return sim

    @njit(cache=True)
    def _apply_hybrid_constraint(sim, rt_ref, rt_exp, bp, adaptive_rt, samples_for_gradient):
        n_ref = sim.shape[0]
        n_exp = sim.shape[1]
        if n_ref < 2 or n_exp < 2:
            return sim

        sampling_time_ref = (rt_ref[-1] - rt_ref[0]) / (n_ref - 1)
        sampling_time_exp = (rt_exp[-1] - rt_exp[0]) / (n_exp - 1)
        if sampling_time_ref == 0.0 or sampling_time_exp == 0.0:
            return sim
        no_beef = int(math.ceil(adaptive_rt / sampling_time_ref))
        max_val = np.max(sim)
        constrain_value = -2.0 * max_val / samples_for_gradient

        out = sim.copy()
        for i in range(n_ref):
            for j in range(n_exp):
                dist = abs((bp[i] - rt_exp[j]) / sampling_time_exp)
                rounded = math.floor(dist + 0.5)
                if rounded > no_beef:
                    out[i, j] += constrain_value * (dist - no_beef)
        return out

    @njit(cache=True)
    def _affine_align(sim, go, ge):
        n_ref = sim.shape[0]
        n_exp = sim.shape[1]
        row_size = n_ref + 1
        col_size = n_exp + 1
        neg_inf = -np.inf

        m_mat = np.zeros((row_size, col_size), dtype=np.float64)
        a_mat = np.zeros((row_size, col_size), dtype=np.float64)
        b_mat = np.zeros((row_size, col_size), dtype=np.float64)
        traceback = np.zeros((3, row_size, col_size), dtype=np.int64)

        for i in range(row_size):
            m_mat[i, 0] = neg_inf
            b_mat[i, 0] = neg_inf
            traceback[_MAT_M, i, 0] = _SS
            traceback[_MAT_B, i, 0] = _SS
        for j in range(col_size):
            m_mat[0, j] = neg_inf
            a_mat[0, j] = neg_inf
            traceback[_MAT_M, 0, j] = _SS
            traceback[_MAT_A, 0, j] = _SS
        m_mat[0, 0] = 0.0

        for i in range(1, row_size):
            a_mat[i, 0] = 0.0
            traceback[_MAT_A, i, 0] = _TA
        if row_size > 1:
            traceback[_MAT_A, 1, 0] = _TM
        for j in range(1, col_size):
            b_mat[0, j] = 0.0
            traceback[_MAT_B, 0, j] = _LB
        if col_size > 1:
            traceback[_MAT_B, 0, 1] = _LM

        for i in range(1, row_size):
            for j in range(1, col_size):
                s_value = sim[i - 1, j - 1]
                diago = m_mat[i - 1, j - 1] + s_value
                insert_a = a_mat[i - 1, j - 1] + s_value
                insert_b = b_mat[i - 1, j - 1] + s_value
                if insert_a >= diago and insert_a >= insert_b:
                    traceback[_MAT_M, i, j] = _DA
                    m_mat[i, j] = insert_a
                if insert_b >= diago and insert_b >= insert_a:
                    traceback[_MAT_M, i, j] = _DB
                    m_mat[i, j] = insert_b
                if diago >= insert_a and diago >= insert_b:
                    traceback[_MAT_M, i, j] = _DM
                    m_mat[i, j] = diago

                a_from_m = m_mat[i - 1, j] - go
                a_from_a = a_mat[i - 1, j] - ge
                a_from_b = b_mat[i - 1, j] - go
                if a_from_a >= a_from_m and a_from_a >= a_from_b:
                    traceback[_MAT_A, i, j] = _TA
                    a_mat[i, j] = a_from_a
                if a_from_b >= a_from_m and a_from_b >= a_from_a:
                    traceback[_MAT_A, i, j] = _TB
                    a_mat[i, j] = a_from_b
                if a_from_m >= a_from_a and a_from_m >= a_from_b:
                    traceback[_MAT_A, i, j] = _TM
                    a_mat[i, j] = a_from_m

                b_from_m = m_mat[i, j - 1] - go
                b_from_a = a_mat[i, j - 1] - go
                b_from_b = b_mat[i, j - 1] - ge
                if b_from_a >= b_from_m and b_from_a >= b_from_b:
                    traceback[_MAT_B, i, j] = _LA
                    b_mat[i, j] = b_from_a
                if b_from_b >= b_from_m and b_from_b >= b_from_a:
                    traceback[_MAT_B, i, j] = _LB
                    b_mat[i, j] = b_from_b
                if b_from_m >= b_from_a and b_from_m >= b_from_b:
                    traceback[_MAT_B, i, j] = _LM
                    b_mat[i, j] = b_from_m

        return _traceback_affine_overlap(m_mat, a_mat, b_mat, traceback)

    @njit(cache=True)
    def _traceback_affine_overlap(m_mat, a_mat, b_mat, traceback):
        row_size = m_mat.shape[0]
        col_size = m_mat.shape[1]
        n_ref = row_size - 1
        n_exp = col_size - 1

        row_idx = 0
        col_idx = col_size - 1
        mat_name = _MAT_M
        max_score = -np.inf

        for i in range(row_size):
            value = m_mat[i, col_size - 1]
            if value >= max_score:
                row_idx = i
                col_idx = col_size - 1
                mat_name = _MAT_M
                max_score = value
            elif a_mat[i, col_size - 1] >= max_score:
                row_idx = i
                col_idx = col_size - 1
                mat_name = _MAT_A
                max_score = a_mat[i, col_size - 1]
            elif b_mat[i, col_size - 1] >= max_score:
                row_idx = i
                col_idx = col_size - 1
                mat_name = _MAT_B
                max_score = b_mat[i, col_size - 1]

        for j in range(col_size):
            value = m_mat[row_size - 1, j]
            if value >= max_score:
                row_idx = row_size - 1
                col_idx = j
                mat_name = _MAT_M
                max_score = value
            elif a_mat[row_size - 1, j] >= max_score:
                row_idx = row_size - 1
                col_idx = j
                mat_name = _MAT_A
                max_score = a_mat[row_size - 1, j]
            elif b_mat[row_size - 1, j] >= max_score:
                row_idx = row_size - 1
                col_idx = j
                mat_name = _MAT_B
                max_score = b_mat[row_size - 1, j]

        max_len = n_ref + n_exp + 2
        idx_a = np.zeros(max_len, dtype=np.int64)
        idx_b = np.zeros(max_len, dtype=np.int64)
        count = 0

        if row_idx != n_ref:
            for i in range(n_ref, row_idx, -1):
                idx_a[count] = i
                idx_b[count] = 0
                count += 1
        elif col_idx != n_exp:
            for j in range(n_exp, col_idx, -1):
                idx_a[count] = 0
                idx_b[count] = j
                count += 1

        trace_pointer = traceback[mat_name, row_idx, col_idx]
        while trace_pointer != _SS:
            if trace_pointer == _DM:
                idx_a[count] = row_idx
                idx_b[count] = col_idx
                count += 1
                row_idx -= 1
                col_idx -= 1
                mat_name = _MAT_M
            elif trace_pointer == _DA:
                idx_a[count] = row_idx
                idx_b[count] = col_idx
                count += 1
                row_idx -= 1
                col_idx -= 1
                mat_name = _MAT_A
            elif trace_pointer == _DB:
                idx_a[count] = row_idx
                idx_b[count] = col_idx
                count += 1
                row_idx -= 1
                col_idx -= 1
                mat_name = _MAT_B
            elif trace_pointer == _TM:
                idx_a[count] = row_idx
                idx_b[count] = 0
                count += 1
                row_idx -= 1
                mat_name = _MAT_M
            elif trace_pointer == _TA:
                idx_a[count] = row_idx
                idx_b[count] = 0
                count += 1
                row_idx -= 1
                mat_name = _MAT_A
            elif trace_pointer == _TB:
                idx_a[count] = row_idx
                idx_b[count] = 0
                count += 1
                row_idx -= 1
                mat_name = _MAT_B
            elif trace_pointer == _LM:
                idx_a[count] = 0
                idx_b[count] = col_idx
                count += 1
                col_idx -= 1
                mat_name = _MAT_M
            elif trace_pointer == _LA:
                idx_a[count] = 0
                idx_b[count] = col_idx
                count += 1
                col_idx -= 1
                mat_name = _MAT_A
            elif trace_pointer == _LB:
                idx_a[count] = 0
                idx_b[count] = col_idx
                count += 1
                col_idx -= 1
                mat_name = _MAT_B
            else:
                break
            trace_pointer = traceback[mat_name, row_idx, col_idx]

        out_a = np.empty(count, dtype=np.int64)
        out_b = np.empty(count, dtype=np.int64)
        for i in range(count):
            out_a[i] = idx_a[count - i - 1]
            out_b[i] = idx_b[count - i - 1]
        return out_a, out_b

    @njit(cache=True)
    def _interpolate_negative(values):
        out = values.copy()
        n = out.size
        start = 0
        while start < n and out[start] < 0.0:
            start += 1
        if start >= n:
            return out
        i = start
        while i < n - 1:
            if out[i] >= 0.0 and out[i + 1] < 0.0:
                left_idx = i
                j = i + 1
                while j < n and out[j] < 0.0:
                    j += 1
                if j >= n:
                    break
                base = out[left_idx]
                target = out[j]
                gaps = j - left_idx
                for k in range(left_idx + 1, j):
                    out[k] = base + (k - left_idx) * (target - base) / gaps
                i = j
            else:
                i += 1
        return out

    @njit(cache=True)
    def _round2(value):
        return np.round(value * 100.0) / 100.0

    @njit(cache=True)
    def _expand_indices_to_times(rt_ref, rt_exp, idx_a, idx_b):
        n = idx_a.size
        t_ref = np.empty(n, dtype=np.float64)
        t_exp = np.empty(n, dtype=np.float64)
        for i in range(n):
            if idx_a[i] != 0:
                t_ref[i] = rt_ref[idx_a[i] - 1]
            else:
                t_ref[i] = -1.0
            if idx_b[i] != 0:
                t_exp[i] = rt_exp[idx_b[i] - 1]
            else:
                t_exp[i] = -1.0

        t_ref = _interpolate_negative(t_ref)
        t_exp = _interpolate_negative(t_exp)

        keep = 0
        for i in range(n):
            if idx_a[i] != 0:
                keep += 1

        out = np.empty((keep, 2), dtype=np.float64)
        out_idx = 0
        for i in range(n):
            if idx_a[i] != 0:
                if t_ref[i] < 0.0:
                    out[out_idx, 0] = np.nan
                else:
                    out[out_idx, 0] = _round2(t_ref[i])
                if t_exp[i] < 0.0:
                    out[out_idx, 1] = np.nan
                else:
                    out[out_idx, 1] = _round2(t_exp[i])
                out_idx += 1
        return out

    @njit(cache=True)
    def _get_aligned_times_cpp_core(
        rt_ref,
        ms2_ref,
        rt_exp,
        ms2_exp,
        bp,
        adaptive_rt,
        kernel_len,
        coeffs,
        go_factor,
        ge_factor,
        cos_angle_thresh,
        dot_prod_thresh,
        gap_quantile,
        samples_for_gradient,
    ):
        smoothed_ref = _smooth_traces(ms2_ref, coeffs, kernel_len)
        smoothed_exp = _smooth_traces(ms2_exp, coeffs, kernel_len)
        sim = _build_dot_product_masked_similarity(
            smoothed_ref, smoothed_exp, cos_angle_thresh, dot_prod_thresh
        )
        gap_penalty = max(0.01, _quantile_type7(sim.ravel(), gap_quantile))
        sim = _apply_hybrid_constraint(sim, rt_ref, rt_exp, bp, adaptive_rt, samples_for_gradient)
        idx_a, idx_b = _affine_align(sim, gap_penalty * go_factor, gap_penalty * ge_factor)
        return _expand_indices_to_times(rt_ref, rt_exp, idx_a, idx_b)

else:
    _get_aligned_times_cpp_core = None
