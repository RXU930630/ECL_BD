"""Independent analytical/adaptive checks of the dry-contact physics."""

import unittest
from pathlib import Path

import numpy as np
from scipy.integrate import quad
from scipy.special import erf, exp1

from persson_contact.core import PerssonModel, separation_kernel
from persson_contact.psd import PowerSpectrum, read_psd


class InverseFourthSpectrum:
    """Analytic C(q)=c*q^-4, independent of production PSD integration."""

    q0 = 1.0
    q1 = 100.0
    c = 1e-4

    def evaluate(self, q):
        return self.c * np.asarray(q) ** -4

    def moment(self, order, lower=None, upper=None):
        lo = self.q0 if lower is None else lower
        hi = self.q1 if upper is None else upper
        if order == 3:
            return self.c * np.log(hi / lo)
        power = order - 3
        return self.c * (hi ** power - lo ** power) / power


class TestCore(unittest.TestCase):
    def test_kernel_against_independent_adaptive_integral(self):
        for gamma in (0.0, 0.4, 1.0):
            for z in (1e-5, 0.01, 0.1, 1.0, 4.0):
                with self.subTest(gamma=gamma, z=z):
                    expected = quad(
                        lambda x: (gamma + 3 * (1-gamma) * erf(x)**2)
                        * np.exp(-x*x) / x, z, 15.0, epsabs=1e-14,
                    )[0]
                    np.testing.assert_allclose(separation_kernel(z, gamma), expected, rtol=2e-7)
            self.assertEqual(separation_kernel(np.inf, gamma), 0.0)

    def test_standard_area_analytic_limits_and_modulus_scaling(self):
        spectrum = InverseFourthSpectrum()
        model = PerssonModel(spectrum, 1e6, 0.4, n_q=64)
        pressures = np.array([0, 100, 2000, 1e8])
        expected = erf(pressures / (1e6 * np.sqrt(np.pi*spectrum.c*np.log(100))))
        np.testing.assert_allclose(model.area(pressures), expected, rtol=1e-14)
        self.assertEqual(model.area(100, q=1.0), 1.0)
        self.assertEqual(model.area(0, q=1.0), 0.0)
        other = PerssonModel(spectrum, 2e6, 1.0, n_q=64)
        np.testing.assert_allclose(other.area(2*pressures), model.area(pressures))

    def test_gamma_one_mean_gap_against_analytic_infinite_band_limit(self):
        # For C=c*q^-4, q0=1 and q1->infinity, gamma=1:
        # h=sqrt(pi*c)*E1(2*p/(E*sqrt(pi*c))). Cutoff error here is <1e-10.
        spectrum = InverseFourthSpectrum()
        spectrum.q1 = np.exp(25)
        spectrum.c = 0.01
        p = 1e6*np.sqrt(np.pi*spectrum.c)*0.1
        expected = np.sqrt(np.pi*spectrum.c)*exp1(0.2)
        coarse = PerssonModel(spectrum, 1e6, 1., n_q=32).mean_gap(p)
        fine = PerssonModel(spectrum, 1e6, 1., n_q=128).mean_gap(p)
        np.testing.assert_allclose(fine, expected, rtol=2e-8)
        self.assertLess(abs(fine-expected), abs(coarse-expected))

    def test_mean_gap_monotonic_pressure_and_gamma_effect(self):
        model = PerssonModel(InverseFourthSpectrum(), 1e6, .4, n_q=128)
        pressures = np.geomspace(100, 1e5, 12)
        gaps = model.mean_gap(pressures)
        self.assertTrue(np.all(gaps > 0))
        self.assertTrue(np.all(np.diff(gaps) < 0))
        uncorrected = PerssonModel(model.psd, 1e6, 1., n_q=128)
        self.assertLess(model.mean_gap(100), uncorrected.mean_gap(100))
        np.testing.assert_allclose(model.mean_gap(10000),
            PerssonModel(model.psd, 2e6, .4, n_q=128).mean_gap(20000), rtol=1e-14)

    def test_mean_gap_split_at_psd_knots_against_independent_quad(self):
        spectrum = PowerSpectrum(
            log10_q=np.array([0., .8, 1.01, 1.03, 1.8, 3.]),
            log10_c=np.array([-4., -6., -7.5, -6.4, -10., -14.]),
            q0=1., qr=1., q1=1000.,
        )
        model = PerssonModel(spectrum, 1e6, .4, n_q=64)
        p = 8000.

        def integrand(logq):
            q = np.exp(logq)
            i3 = spectrum.moment(3, lower=1., upper=q)
            return q**3*spectrum.evaluate(q)/np.sqrt(i3)*separation_kernel(
                p/(1e6*np.sqrt(np.pi*i3)), .4)

        bounds = spectrum.log10_q*np.log(10.)
        reference = sum(quad(integrand, a, b, epsabs=1e-13, limit=200)[0]
                        for a, b in zip(bounds[:-1], bounds[1:]))
        np.testing.assert_allclose(model.mean_gap(p), reference, rtol=2e-6)

    def test_pdf_separate_atom_mass_conditioning_and_reported_moment(self):
        model = PerssonModel(InverseFourthSpectrum(), 1e6, .4, n_q=64)
        # CDF bins resolve the sharp near-zero tail where trapz(point PDF) fails.
        gaps = np.r_[0., np.geomspace(1e-40, 1e12, 3000)]
        result = model.gap_distribution(2000., gaps)
        diag = result.diagnostics
        self.assertEqual(result.contact_fraction, model.area(2000.))
        self.assertTrue(np.all(result.continuous_density_per_m >= 0))
        self.assertTrue(np.all(np.diff(result.continuous_cdf) >= -1e-15))
        self.assertEqual(result.continuous_cdf[0], 0)
        self.assertAlmostEqual(result.contact_fraction + diag.resolved_continuous_mass
                               + diag.unresolved_continuous_mass, 1.+diag.quadrature_mass_error, places=14)
        self.assertAlmostEqual(result.continuous_cdf[-1], diag.resolved_continuous_mass, places=12)
        self.assertAlmostEqual(result.positive_cdf[-1], diag.resolved_continuous_mass/(1-result.contact_fraction), places=12)
        np.testing.assert_allclose(result.positive_density_per_m,
                                   result.continuous_density_per_m/(1-result.contact_fraction))
        self.assertAlmostEqual(diag.first_moment_relative_error,
                               diag.analytic_first_moment_m/result.mean_gap_m-1)
        self.assertEqual(diag.approximation, "fortran_method1_reflected_gaussian")
        self.assertGreater(diag.boundary_unresolved_mass, 0)
        self.assertTrue(any("Initial scale area <0.99" in note for note in diag.notes))

    def test_truncated_gap_grid_is_not_renormalized(self):
        model = PerssonModel(InverseFourthSpectrum(), 1e6, .4, n_q=64)
        result = model.gap_distribution(2000., np.linspace(.02, .025, 50))
        self.assertLess(result.diagnostics.integrated_continuous_mass_on_grid, .5)
        self.assertLess(result.diagnostics.sampled_continuous_mass, .5)
        self.assertTrue(any("misses mass" in note for note in result.diagnostics.notes))

    def test_real_psd_reference_means_and_unresolved_probability_accounting(self):
        path = Path(__file__).resolve().parents[1] / 'examples/data/roughness_psd.dat'
        spectrum = read_psd(path)
        model = PerssonModel(spectrum, 1e7, .4, n_q=128)
        # Independent adaptive quadrature split at all input PSD knots; these
        # references are not results from the original Fortran executable.
        reference = [2.0315728006966916e-6, 6.528334657692934e-7, 1.1354135965527952e-8]
        np.testing.assert_allclose(model.mean_gap([1e4, 1e5, 1e6]), reference, rtol=2e-7)
        grid = np.r_[0., np.geomspace(1e-14, 20*spectrum.hrms, 1000)]
        result = model.gap_distribution(1e6, grid)
        diag = result.diagnostics
        self.assertGreater(diag.unresolved_continuous_mass, 0.)
        self.assertAlmostEqual(result.contact_fraction + diag.resolved_continuous_mass
                               + diag.unresolved_continuous_mass, 1.+diag.quadrature_mass_error, places=14)
        self.assertLessEqual(result.positive_cdf[-1], 1.)
        self.assertTrue(np.all(np.diff(result.continuous_cdf) >= -1e-15))
        self.assertAlmostEqual(diag.resolved_continuous_mass + diag.unresolved_continuous_mass,
                               diag.quadrature_continuous_mass, places=14)

    def test_invalid_model_inputs(self):
        for kwargs in ({"effective_modulus_pa": 0}, {"gamma": -1}, {"gamma": 1.1},
                       {"n_q": 2}, {"n_q": 32.5}, {"pdf_afact": -1},
                       {"pdf_log_step": 0}, {"pdf_log_step": True}, {"pdf_min_u1_m": -1},
                       {"pdf_max_relative_contact_area": 1},
                       {"pdf_min_relative_contact_area": .9999}):
            with self.subTest(kwargs=kwargs):
                arguments = dict(effective_modulus_pa=1e6, gamma=.4, n_q=64)
                arguments.update(kwargs)
                with self.assertRaises(ValueError):
                    PerssonModel(InverseFourthSpectrum(), **arguments)

    def test_invalid_load_gap_and_full_contact_cases(self):
        model = PerssonModel(InverseFourthSpectrum(), 1e6, .4, n_q=64)
        with self.assertRaises(ValueError):
            model.mean_gap(0)
        with self.assertRaises(ValueError):
            model.area(-1)
        with self.assertRaises(ValueError):
            model.area(1, q=101)
        with self.assertRaises(ValueError):
            model.gap_distribution(2000, [1e-4, 1e-5])
        with self.assertRaisesRegex(ValueError, "full contact"):
            model.gap_distribution(1e10, [0, 1])


if __name__ == "__main__":
    unittest.main()
