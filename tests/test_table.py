"""Database checks using actual constitutive calculations, including off-grid queries."""
from pathlib import Path
import unittest

import numpy as np

from persson_contact.config import load_problem
from persson_contact.core import PerssonModel
from persson_contact.table import build_contact_table


class ContactTableTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        config = Path(__file__).resolve().parents[1] / 'examples' / 'dry_contact.json'
        problem, options = load_problem(config)
        cls.model = PerssonModel(problem.psd, problem.effective_modulus_pa,
                                 options.gamma, n_q=128, pdf_afact=options.pdf_afact)
        cls.table = build_contact_table(cls.model, np.geomspace(1e4, 1e6, 48))

    def test_nodes_and_off_grid_queries(self):
        table = self.table
        response = table.evaluate(table.mean_gaps_m)
        np.testing.assert_allclose(response.pressure_pa, table.pressures_pa, rtol=1e-12)
        h = (table.mean_gaps_m[10:-10] + table.mean_gaps_m[11:-9]) / 2
        response = table.evaluate(h)
        np.testing.assert_allclose(self.model.mean_gap(response.pressure_pa), h, rtol=5e-4)
        np.testing.assert_allclose(response.contact_fraction, self.model.area(response.pressure_pa), rtol=1e-14)

    def test_no_silent_extrapolation(self):
        lo, hi = self.table.valid_gap_range_m
        for value in (lo * 0.9, hi * 1.1, float('nan'), -1.0):
            with self.assertRaises(ValueError):
                self.table.evaluate(value)

    def test_invalid_pressure_range(self):
        for p in ([1, 1, 2], [3, 2, 1], [0, 1, 2], [1, 2]):
            with self.assertRaises(ValueError):
                build_contact_table(self.model, p)

    def test_pdf_uses_same_pressure_and_area(self):
        h = float(np.mean(self.table.mean_gaps_m[20:22]))
        hrms = float(self.model.psd.hrms)
        response = self.table.evaluate(h, gaps_m=np.r_[0., np.geomspace(hrms * 1e-6, 20 * hrms, 300)])
        self.assertAlmostEqual(response.contact_fraction, response.gap_distribution.contact_fraction, places=13)
        self.assertLess(abs(response.interpolation_mean_gap_error_m) / h, 5e-4)


if __name__ == '__main__':
    unittest.main()
