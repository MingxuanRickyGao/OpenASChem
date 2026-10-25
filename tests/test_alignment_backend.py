"""Run with python -m unittest discover -s tests (no R installation needed)."""

import subprocess
import sys
import unittest
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from openaschem.open_identify_utils import align_one_compound_standards, do_align_linear


class AlignmentBackendTests(unittest.TestCase):
    def test_r_2_12_reference_paths(self):
        # Ten deterministic synthetic pairs: identity, shifts, zero/missing
        # signals, amplitude, multiple peaks, global fit, noise, unequal grids.
        with np.load(Path(__file__).parent / "fixtures/alignment_r_2_12.npz") as z:
            for name in z["names"]:
                with self.subTest(case=name):
                    rt_ref, rt_exp = z[name + "_rt_ref"], z[name + "_rt_exp"]
                    result = do_align_linear(
                        rt_ref,
                        [z[name + "_signal_ref"]],
                        rt_exp,
                        [z[name + "_signal_exp"]],
                        z[name + "_fit"],
                    )
                    np.testing.assert_allclose(result, z[name + "_expected"], rtol=0, atol=1e-10)

    def test_three_standard_interface(self):
        rt = np.arange(81) * 0.5
        signal = np.exp(-0.5 * ((rt - 20) / 1.3) ** 2) * 1000
        result = align_one_compound_standards(
            rt, rt, rt, signal, signal, signal, 1, 0, 5, [("S1", "S2"), ("S1", "S3")], ["S1", "S1"]
        )
        np.testing.assert_allclose(result.to_numpy(dtype=float), np.column_stack([rt] * 3))

    def test_trace_rt_length_mismatch_is_rejected(self):
        with self.assertRaisesRegex(ValueError, "ms2_ref shape"):
            do_align_linear(np.arange(20), [np.zeros(19)], np.arange(20), [np.zeros(20)], [1, 0, 5])

    def test_cli_import_and_help_without_r(self):
        code = """import sys, runpy
class NoR:
    def find_spec(self, fullname, path=None, target=None):
        if fullname == 'rpy2' or fullname.startswith('rpy2.'):
            raise RuntimeError('rpy2 is forbidden')
sys.meta_path.insert(0, NoR())
sys.argv = ['OpenASChem.py', 'openidentify', '--help']
runpy.run_path('OpenASChem.py', run_name='__main__')
"""
        result = subprocess.run(
            [sys.executable, "-B", "-c", code], cwd=ROOT, capture_output=True, text=True
        )
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn("--quantification_method", result.stdout)
        self.assertNotIn("--r_home", result.stdout)


if __name__ == "__main__":
    unittest.main()
