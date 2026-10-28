"""Pressure generation and an analytical test of the critical u1 definition."""
import math
import unittest

import numpy as np
from scipy.special import erf, exp1

from persson_contact.config import pressure_grid, Options
from persson_contact.core import PerssonModel
from persson_contact.psd import PowerSpectrum


class PressureScanTests(unittest.TestCase):
    def test_defaults_keep_count_bounds_and_reference(self):
        p = pressure_grid({'reference_pa': 1e5})
        self.assertEqual(len(p), 110)
        self.assertEqual(p[0], 10.)
        self.assertEqual(p[-1], 2e5)
        self.assertEqual(np.count_nonzero(p == 1e5), 1)
        self.assertTrue(np.all(np.diff(p) > 0))

    def test_reference_at_upper_bound_and_small_point_count(self):
        for factor in (1., 1.001, 2., 10.):
            p = pressure_grid({'reference_pa': 20., 'upper_factor': factor, 'number_of_points': 3})
            self.assertEqual(len(p), 3)
            self.assertIn(20., p)
            self.assertEqual(p[-1], 20*factor)
            self.assertTrue(np.all(np.diff(p) > 0))

    def test_bad_scan_or_area_options_rejected(self):
        for raw in ({'reference_pa': 0}, {'reference_pa': 1., 'upper_factor': .5},
                    {'reference_pa': 1., 'number_of_points': 2},
                    {'reference_pa': 1., 'number_of_points': True}):
            with self.subTest(raw=raw), self.assertRaises(ValueError):
                pressure_grid(raw)
        for opts in ({'critical_relative_contact_area': 1.}, {'critical_area_tolerance': 0.}):
            with self.assertRaises(ValueError):
                Options(.4, 64, 1., **opts)


class CriticalSeparationTests(unittest.TestCase):
    def model(self):
        # C=c*q^-4, q0=1, q1~infinity, gamma=1 has a known conditional ubar.
        c = .01
        logq1 = 25/math.log(10)
        psd = PowerSpectrum([0., logq1], [math.log10(c), math.log10(c)-4*logq1],
                            1., 1., math.exp(25))
        return PerssonModel(psd, 1e6, 1., n_q=256)

    def test_critical_u1_against_analytical_subspectrum_derivative(self):
        model = self.model()
        r = .1
        p = r*1e6*math.sqrt(math.pi*.01)
        result = model.critical_separation(p)
        self.assertEqual(result.status, 'found')
        self.assertLessEqual(abs(result.area_deviation), .005)
        mu = math.log(result.magnification)
        area = erf(r/math.sqrt(mu))
        area_prime = -r*math.exp(-r*r/mu)/(math.sqrt(math.pi)*mu**1.5)
        k = math.sqrt(math.pi*.01)/result.magnification
        ubar = k*exp1(2*r/area)
        u1 = ubar*(1-area/area_prime) + k*math.exp(-2*r/area)
        np.testing.assert_allclose(result.mean_gap_at_scale_m, ubar, rtol=1e-8)
        np.testing.assert_allclose(result.u1_m, u1, rtol=3e-3)
        self.assertEqual(result.critical_separation_m, result.u1_m)
        self.assertGreater(result.u1_m, result.mean_gap_at_scale_m)

    def test_no_crossing_is_not_replaced_by_nearest_endpoint(self):
        result = self.model().critical_separation(1e9)
        self.assertEqual(result.status, 'threshold_not_reached')
        self.assertIsNone(result.critical_separation_m)


if __name__ == '__main__':
    unittest.main()
