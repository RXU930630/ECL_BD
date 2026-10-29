"""Real Fortran tail regressions: no synthetic reference or fitted parameters."""
import hashlib
import json
from pathlib import Path
import unittest

import numpy as np

from persson_contact.config import Options
from persson_contact.core import PerssonModel, ModelResolutionError
from persson_contact.psd import read_psd


class FortranGapReferenceTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        root = Path(__file__).resolve().parents[1]
        cls.psd_path = root/'examples/data/roughness_psd.dat'
        cls.psd = read_psd(cls.psd_path)
        cls.fixture = json.loads((root/'tests/data/fortran_gap_pdf.json').read_text())

    def test_all_moduli_pressures_and_far_tails(self):
        self.assertEqual(hashlib.sha256(self.psd_path.read_bytes()).hexdigest(),
                         self.fixture['psd_sha256'])
        self.assertEqual(len(self.fixture['cases']), 10)
        for sample in self.fixture['cases']:
            with self.subTest(case=sample['name']):
                case = sample['case']
                model = PerssonModel(self.psd, case['effective_modulus_pa'], .4,
                                     pdf_log_step=case['delzeta'])
                reference = np.array(sample['density_per_m'])
                pdf = model.gap_distribution(case['pressure_pa'], sample['gaps_m'])
                # Includes the very last point, at densities down to 1e-18/m.
                # Remaining kernel-table error is strongly amplified there.
                np.testing.assert_allclose(pdf.continuous_density_per_m, reference,
                                           rtol=.05, atol=0)
                self.assertLess(abs(pdf.diagnostics.analytic_first_moment_m
                                    /sample['first_moment_m']-1), 3e-5)
                diag = pdf.diagnostics
                self.assertAlmostEqual(diag.resolved_continuous_mass+diag.unresolved_continuous_mass,
                                       diag.quadrature_continuous_mass, places=14)
                self.assertAlmostEqual(diag.unresolved_continuous_mass,
                    diag.excluded_high_area_mass+diag.excluded_low_area_mass+
                    diag.excluded_small_u1_mass+diag.boundary_unresolved_mass, places=14)

    def test_tail_disappears_without_clipping_output_gaps(self):
        # Old code gave isolated Gaussian peaks above this reference tail.
        # Extending the requested output domain must preserve natural decay.
        model = PerssonModel(self.psd, 4e6, .4)
        gaps = np.array([5e-6, 8e-6, 1e-5, 2e-5, 5e-5, 1e-4])
        result = model.gap_distribution(1e5, gaps)
        self.assertTrue(np.all(np.diff(result.continuous_density_per_m) <= 0))
        self.assertLess(result.continuous_density_per_m[-1], 1e-250)
        self.assertGreater(result.diagnostics.excluded_high_area_mass, 0)
        # The mass is neither renormalized nor moved into the contact atom.
        self.assertEqual(result.contact_fraction, model.area(1e5))
        self.assertLess(result.continuous_cdf[-1], 1-result.contact_fraction)

    def test_pdf_scale_resolution_is_separate_from_energy_resolution(self):
        gaps = np.geomspace(1e-12, 1e-5, 50)
        coarse = PerssonModel(self.psd, 4e6, .4, n_q=32).gap_distribution(1e5, gaps)
        fine = PerssonModel(self.psd, 4e6, .4, n_q=256).gap_distribution(1e5, gaps)
        np.testing.assert_array_equal(coarse.continuous_density_per_m, fine.continuous_density_per_m)
        for values in ({'pdf_log_step':0}, {'pdf_min_u1_m':0},
                       {'pdf_max_relative_contact_area':1},
                       {'pdf_min_relative_contact_area':.9999}):
            with self.assertRaises(ValueError):
                Options(.4,64,1.,**values)
        with self.assertRaises(ModelResolutionError):
            PerssonModel(self.psd,4e6,.4,pdf_log_step=1e-6).gap_distribution(1e5,gaps)


if __name__ == '__main__':
    unittest.main()
