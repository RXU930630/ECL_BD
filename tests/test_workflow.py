"""User-facing switches, exclusive material inputs and real output contracts."""
from copy import deepcopy
import json
from pathlib import Path
import tempfile
import unittest
import numpy as np

from persson_contact.material import parse_material
from persson_contact.workflow import run_project


class MaterialInterfaceTests(unittest.TestCase):
    def test_curve_needs_no_fixed_modulus(self):
        with tempfile.TemporaryDirectory() as directory:
            material = parse_material({'model': 'viscoelastic', 'poisson_ratio': .5,
                                       'master_curve_file': 'curve.txt'}, directory)
            self.assertIsNone(material.effective_modulus_pa)
            self.assertEqual(material.master_curve_file, Path(directory) / 'curve.txt')
            with self.assertRaisesRegex(NotImplementedError, '尚未实现'):
                material.elastic_modulus()

    def test_material_sources_are_exclusive(self):
        for key in ('young_modulus_pa', 'effective_modulus_pa'):
            with self.subTest(key=key), self.assertRaisesRegex(ValueError, '不应同时'):
                parse_material({'model': 'viscoelastic', 'poisson_ratio': .5,
                                'master_curve_file': 'curve.txt', key: 1e6})
        with self.assertRaises(ValueError):
            parse_material({'model': 'elastic', 'young_modulus_pa': 1e6,
                            'poisson_ratio': .5, 'master_curve_file': 'curve.txt'})

    def test_invalid_curve_specifications(self):
        for raw in ({'model': 'viscoelastic', 'poisson_ratio': .5},
                    {'model': 'viscoelastic', 'master_curve_file': 'x'},
                    {'model': 'viscoelastic', 'poisson_ratio': .6, 'master_curve_file': 'x'},
                    {'model': 'unknown'}):
            with self.subTest(raw=raw), self.assertRaises(ValueError):
                parse_material(raw)


class WorkflowTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.directory = Path(self.temp.name)
        (self.directory / 'surface.dat').write_text('surface\n1\n1\n100\n1\n1\n0 -8\n1 -10\n2 -12\n')
        self.parameters = {
            'description': 'workflow test',
            'psd': {'path': 'surface.dat'},
            'material': {'model': 'elastic', 'young_modulus_pa': 3e6, 'poisson_ratio': .5},
            'pressures_pa': [100., 300., 1000.],
            'options': {'gamma': .4, 'n_q': 32, 'pdf_afact': 1., 'n_gap': 64},
            'output_directory': 'result', 'leakage': {}, 'mixed_friction': {},
        }

    def test_contact_only_writes_results_and_preserves_parameter_provenance(self):
        summary = run_project(self.parameters, base_directory=self.directory,
                              compute_gap_pdf=False, save_plots=False)
        result = self.directory / 'result'
        self.assertTrue((result / 'contact_table.csv').is_file())
        self.assertFalse(list(result.glob('gap_*.csv')))
        saved = json.loads((result / 'diagnostics.json').read_text())
        self.assertEqual(saved['effective_modulus_pa'], 4e6)
        self.assertEqual(saved['parameters_snapshot'], self.parameters)
        self.assertFalse(summary['compute_gap_pdf'])
        with self.assertRaises(FileExistsError):
            run_project(self.parameters, base_directory=self.directory,
                        compute_gap_pdf=False, save_plots=False)

    def test_disabled_tasks_do_not_require_inputs_or_create_outputs(self):
        self.assertIsNone(run_project({}, base_directory=self.directory, run_contact=False))
        self.assertFalse((self.directory / 'result').exists())

    def test_positive_density_export_keeps_zero_origin_probability_bin(self):
        summary = run_project(self.parameters, base_directory=self.directory,
                              compute_gap_pdf=True, save_plots=False)
        for record in summary['pdf_diagnostics']:
            bins = np.genfromtxt(self.directory / 'result' / record['bin_file'], delimiter=',', names=True)
            density = np.genfromtxt(self.directory / 'result' / record['file'], delimiter=',', names=True)
            self.assertEqual(bins['left_microgap_m'][0], 0.)
            self.assertTrue(np.all(density['microgap_m'] > 0))
            self.assertAlmostEqual(np.sum(bins['nominal_area_probability']), record['exported_bin_mass'], places=12)

    def test_reserved_applications_fail_before_contact_calculation(self):
        for switches in ({'run_mixed_friction': True},
                         {'run_leakage': True, 'run_mixed_friction': True}):
            with self.subTest(switches=switches), self.assertRaises(NotImplementedError):
                run_project(self.parameters, base_directory=self.directory, **switches)
        self.assertFalse((self.directory / 'result').exists())

    def test_curve_selection_does_not_fall_back_to_elastic_or_read_psd(self):
        raw = deepcopy(self.parameters)
        raw['material'] = {'model': 'viscoelastic', 'poisson_ratio': .5,
                           'master_curve_file': 'future_curve.txt'}
        raw['psd']['path'] = 'also_not_present.dat'
        with self.assertRaisesRegex(NotImplementedError, '主曲线'):
            run_project(raw, base_directory=self.directory)
        self.assertFalse((self.directory / 'result').exists())


if __name__ == '__main__':
    unittest.main()
