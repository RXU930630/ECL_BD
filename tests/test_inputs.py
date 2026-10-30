"""Independent analytic checks of PSD parsing, moments and configuration."""

import json
import math
from pathlib import Path
import tempfile
import unittest

import numpy as np

from persson_contact.config import Options, load_problem
from persson_contact.psd import PowerSpectrum, read_psd


class InputFixture(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.directory = Path(self.tmp.name)
        self.path = self.directory / "surface.psd"
        self.path.write_text("surface\n1D0\n1\n100\n999\n0.8\n0 0\n1 -2\n2 -4\n")


class SpectrumTests(InputFixture):
    def test_fortran_numbers_and_header_preservation(self):
        psd = read_psd(self.path)
        self.assertEqual(psd.h0_header, 999)
        self.assertAlmostEqual(psd.evaluate(math.sqrt(10)), 0.1)
        np.testing.assert_allclose(psd.evaluate([1, 10, 100]), [1, 0.01, 0.0001])
        self.assertFalse(psd.log10_q.flags.writeable)
        self.assertFalse(psd.log10_c.flags.writeable)

    def test_analytic_power_law_and_logarithmic_primitive(self):
        # C=q^-2: integral q*C dq=ln(q1/q0), integral q^3*C dq=(q1^2-q0^2)/2.
        psd = read_psd(self.path)
        self.assertAlmostEqual(psd.moment(1), math.log(100), places=12)
        self.assertAlmostEqual(psd.moment(1, lower=2, upper=35), math.log(35/2), places=12)
        self.assertAlmostEqual(psd.moment(3), (100**2 - 1)/2, places=8)
        self.assertAlmostEqual(psd.hrms, math.sqrt(2*math.pi*math.log(100)), places=12)
        self.assertAlmostEqual(psd.rms_slope, math.sqrt(math.pi*(100**2-1)), places=10)
        self.assertEqual(psd.moment(3, 2, 2), 0)

    def test_scale_c_does_not_change_original_data(self):
        psd = read_psd(self.path, scale_c=9)
        np.testing.assert_array_equal(psd.log10_c, [0, -2, -4])
        self.assertAlmostEqual(psd.evaluate(10), 0.09)
        self.assertAlmostEqual(psd.hrms / read_psd(self.path).hrms, 3)

    def test_adjacent_bounds_with_identical_rounded_logs(self):
        psd = read_psd(self.path)
        lower = 10.0
        upper = math.nextafter(lower, math.inf)
        self.assertGreater(upper, lower)
        self.assertEqual(math.log(lower), math.log(upper))
        # For C=q^-2, the exact order-3 integral is (upper^2-lower^2)/2.
        exact = (upper-lower)*(upper+lower)/2
        result = psd.moment(3, lower, upper)
        self.assertGreaterEqual(result, 0.0)
        self.assertLessEqual(abs(result-exact), np.finfo(float).eps*psd.moment(3))

    def test_out_of_range_rejected_instead_of_extrapolated(self):
        psd = read_psd(self.path)
        for q in (0, -1, 0.99, 101, float("nan")):
            with self.subTest(q=q), self.assertRaises(ValueError):
                psd.evaluate(q)
        with self.assertRaises(ValueError):
            psd.moment(1, 5, 2)
        with self.assertRaises(ValueError):
            psd.moment(1, 1, 1000)
        self.assertAlmostEqual(psd.evaluate(100*(1+1e-14)), psd.evaluate(100))

    def test_invalid_samples_and_cutoffs_rejected(self):
        for x, y, q0, qr, q1 in [
            ([0, 0], [0, 1], 1, 1, 10),
            ([1, 0], [0, 1], 1, 1, 10),
            ([0, 1], [0, float("nan")], 1, 1, 10),
            ([0, 1], [0, 1], 1, 1, 100),
            ([0, 1], [0, 1], 1, 20, 10),
            ([0, 1], [0, 1], 0, 1, 10),
        ]:
            with self.subTest(x=x, y=y, q0=q0, qr=qr, q1=q1), self.assertRaises(ValueError):
                PowerSpectrum(x, y, q0, qr, q1)
        for scale in (0, -1, float("nan"), True):
            with self.subTest(scale=scale), self.assertRaises(ValueError):
                read_psd(self.path, scale)

    def test_piecewise_moments_use_each_segment(self):
        # C=1 from q=1..10; C=100/q^2 from 10..100.
        psd = PowerSpectrum([0, 1, 2], [0, 0, -2], 1, 1, 100)
        expected = (10**2 - 1)/2 + 100*math.log(10)
        self.assertAlmostEqual(psd.moment(1), expected, places=10)

    def test_numeric_spectrum_retains_rolloff_below_integration_interval(self):
        self.path.write_text("historical spectrum\n10\n1\n100\n1\n1\n1 -2\n2 -4\n")
        original = self.path.read_bytes()
        psd = read_psd(self.path)
        self.assertEqual(psd.qr, 1.0)
        self.assertEqual(psd.q0, 10.0)
        self.assertAlmostEqual(psd.moment(1), math.log(10), places=12)
        self.assertAlmostEqual(psd.moment(3), (100**2-10**2)/2, places=9)
        self.assertEqual(self.path.read_bytes(), original)


class ConfigurationTests(InputFixture):
    def config(self):
        return {"psd": {"path": "surface.psd", "tail_is_extrapolated": True},
                "description": "Demonstration only",
                "material": {"young_modulus_pa": 1e6, "poisson_ratio": 0.49},
                "pressures_pa": [1e3, 1e4, 1e5],
                "options": {"gamma": 0.4, "n_q": 128, "pdf_afact": 0.5}}

    def load(self, raw):
        path = self.directory / "problem.json"
        path.write_text(json.dumps(raw))
        return load_problem(path)

    def test_relative_path_and_single_effective_modulus_correction(self):
        problem, options = self.load(self.config())
        self.assertAlmostEqual(problem.effective_modulus_pa, 1e6/(1-0.49**2))
        self.assertTrue(problem.psd_tail_is_extrapolated)
        self.assertEqual(problem.description, "Demonstration only")
        self.assertEqual(options.n_q, 128)
        raw = self.config()
        raw["material"] = {"effective_modulus_pa": 2e6}
        self.assertEqual(self.load(raw)[0].effective_modulus_pa, 2e6)
        raw["material"] = {"young_modulus_pa": 1e6, "poisson_ratio": 0.5}
        self.assertAlmostEqual(self.load(raw)[0].effective_modulus_pa, 4e6/3)

    def test_ambiguous_material_invalid_pressure_and_units_fail(self):
        for change in ("modulus", "pressure", "units", "gamma", "grid", "afact", "unknown"):
            raw = self.config()
            if change == "modulus": raw["material"]["effective_modulus_pa"] = 2e6
            if change == "pressure": raw["pressures_pa"] = [0, 1]
            if change == "units": raw["psd"]["q_unit"] = "cycles/mm"
            if change == "gamma": del raw["options"]["gamma"]
            if change == "grid": raw["options"]["n_q"] = 10
            if change == "afact": raw["options"]["pdf_afact"] = -1
            if change == "unknown": raw["options"]["nqq"] = 128
            with self.subTest(change=change), self.assertRaises(ValueError):
                self.load(raw)

    def test_research_options_explicit_and_validated(self):
        for args in ((-0.1, 128, 0.5), (1.1, 128, 0.5), (0.4, True, 0.5),
                     (0.4, 128, float("inf"))):
            with self.subTest(args=args), self.assertRaises(ValueError):
                Options(*args)


if __name__ == "__main__":
    unittest.main()
