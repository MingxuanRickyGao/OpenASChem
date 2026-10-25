"""Regression checks for shoulder filtering and physical standard widths."""

import sys
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import openaschem.open_identify_utils as u
from openaschem.chromatogram_processing import find_candidate_peaks, find_top_peaks


class PeakSelectionTests(unittest.TestCase):
    def test_shallow_ripple_is_not_an_independent_peak(self):
        signal = np.array([0, 20, 80, 100, 99, 101, 80, 20, 0.0])
        old, _ = find_candidate_peaks(signal, 0)
        new, _ = find_candidate_peaks(signal)
        np.testing.assert_array_equal(old, [3, 5])
        np.testing.assert_array_equal(new, [5])

    def test_isolated_weak_peak_is_preserved(self):
        peaks, _ = find_candidate_peaks([0, 0, 1000, 0, 0, 5, 0, 0])
        np.testing.assert_array_equal(peaks, [2, 5])

    def test_resolved_neighboring_peaks_are_preserved(self):
        peaks, _ = find_candidate_peaks([0, 50, 100, 60, 20, 60, 80, 40, 0])
        np.testing.assert_array_equal(peaks, [2, 6])

    def test_filtering_precedes_top_n_selection(self):
        signal = np.array([0, 20, 100, 99, 101, 20, 0, 0, 5, 0.0])
        bounds = find_top_peaks(signal, 2)
        self.assertEqual(len(bounds), 2)
        maxima = sorted(max(signal[left : right + 1]) for left, right in bounds)
        self.assertEqual(maxima, [5, 101])

    def test_empty_zero_and_monotone_traces_have_no_peaks(self):
        for signal in [[], [0] * 15, list(range(15)), list(range(15, 0, -1))]:
            self.assertEqual(find_top_peaks(signal, 3), [])

    def test_invalid_relative_prominence_is_rejected(self):
        for threshold in [-0.1, 1.1, np.nan]:
            with self.assertRaises(ValueError):
                find_candidate_peaks([0, 2, 0], threshold)

    def identify_on_rt_grids(self, intervals):
        signal = 1e6 * np.exp(-0.5 * ((np.arange(41) - 20) / 5) ** 2)
        spectra = [SimpleNamespace(rt_list=(100 + np.arange(41) * dt).tolist()) for dt in intervals]
        xics = [{"compound": signal} for _ in range(3)]
        identity = pd.DataFrame({name: np.arange(41) for name in ["S1", "S2", "S3"]})
        with (
            patch.object(u, "align_one_compound_standards", return_value=identity),
            patch.object(u, "time_vector_to_index_vector", return_value=identity),
        ):
            standard, _ = u.identify_peaks(
                pd.DataFrame({"id": ["compound"]}),
                {"COMPOUND_ID": "id"},
                spectra,
                xics,
                xics,
                {},
                {},
                {},
                [],
                1,
                0,
                5,
                [("S1", "S2"), ("S1", "S3")],
                ["S1", "S1"],
                1,
                0.45,
                100,
                100000,
                50,
                300,
                10.0,
                5,
                "intensity",
            )
        return standard, find_top_peaks(signal, 1)[0]

    def test_broad_peak_passes_seconds_limit_despite_more_than_ten_indices(self):
        standard, (left, right) = self.identify_on_rt_grids([0.5, 0.6, 0.7])
        self.assertGreater(right - left, 10)
        self.assertIn("compound", standard)
        self.assertLess(standard["compound"][0]["rb3"] - standard["compound"][0]["lb3"], 10)

    def test_width_limit_is_checked_in_every_original_rt_grid(self):
        standard, _ = self.identify_on_rt_grids([0.5, 0.6, 1.0])
        self.assertNotIn("compound", standard)


if __name__ == "__main__":
    unittest.main()
