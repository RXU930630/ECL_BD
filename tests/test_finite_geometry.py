"""Physical surface extent, load conversion and finite edge interpretation."""
from copy import deepcopy
import csv
import json
from pathlib import Path
import tempfile
import unittest

import numpy as np

from parameters import PARAMETERS
from persson_contact.geometry import FiniteCylindricalPad,line_load,transverse_width
from persson_contact.macro_output import prepare_macro,run_macro_project
from persson_contact.shear_fit_output import prepare_fit,_signature

ROOT=Path(__file__).resolve().parents[1]


def parameters():
    p=deepcopy(PARAMETERS)
    p['macro_contact']['geometry']='finite_cylindrical_pad'
    p['macro_contact'].pop('line_load_n_m')
    p['macro_contact']['normal_force_n']=100.
    p['material']['young_modulus_pa']=1.8e6
    return p


class FiniteGeometryTests(unittest.TestCase):
    def test_actual_circular_sag_and_projected_or_arc_extent(self):
        g=FiniteCylindricalPad(.17,1000,2.4e6,.1,.03)
        self.assertEqual(g.domain_half_width_m,.015)
        x=np.array([-.015,0,.015])
        np.testing.assert_allclose(g.profile(x),.17-np.sqrt(.17**2-x*x),atol=1e-17)
        arc=FiniteCylindricalPad(.17,1000,2.4e6,.1,.03,'arc_length')
        self.assertAlmostEqual(arc.domain_half_width_m,.17*np.sin(.03/.34))
        self.assertLess(arc.domain_half_width_m,g.domain_half_width_m)
        with self.assertRaisesRegex(ValueError,'outside'):g.profile([.0151])

    def test_load_width_single_source_and_ignored_numerical_halfwidth(self):
        p=parameters()
        self.assertEqual(transverse_width(p),.1)
        for force in (40,100,160):
            p['macro_contact']['normal_force_n']=force
            self.assertEqual(line_load(p),force/.1)
        p['macro_contact']['domain_half_width_m']=None
        args=prepare_macro(p,ROOT)
        np.testing.assert_allclose(args[3].edges_m[[0,-1]],[-.015,.015])
        p['contact_friction']['transverse_width_m']=.01
        with self.assertRaisesRegex(ValueError,'conflicts'):prepare_macro(p,ROOT)
        p['contact_friction']['transverse_width_m']=None
        p['macro_contact']['line_load_n_m']=1600.
        with self.assertRaisesRegex(ValueError,'not both'):prepare_macro(p,ROOT)

    def test_finite_smooth_full_contact_not_rejected_as_bad_far_field(self):
        p=parameters();p['macro_contact']['contact_model']='smooth'
        p['macro_contact']['normal_force_n']=1000.
        p['macro_contact']['geometry_parameters']['finite_cylindrical_pad']['sliding_length_m']=.01
        p['macro_contact']['solver']['cells']=65
        with tempfile.TemporaryDirectory() as d:
            p['output_directory']=str(Path(d)/'finite')
            result=run_macro_project(p,base_directory=ROOT,save_plots=False)
            self.assertFalse(result['far_field_edge_pressure_check_applied'])
            self.assertGreater(result['edge_pressure_ratio'],.1)
            self.assertLess(result['load_relative_error'],1e-10)
            self.assertLess(result['geometry_relative_residual'],1e-9)
            self.assertIsNone(result['hertz_smooth_pressure_L1_relative_error'])
            self.assertEqual(result['geometry']['surface_interval_m'],[-.005,.005])

    def test_fit_csv_total_forces_use_geometry_not_legacy_width(self):
        p=parameters()
        with tempfile.TemporaryDirectory() as d:
            path=Path(d)/'experiment.csv'
            with path.open('w',newline='') as stream:
                writer=csv.DictWriter(stream,fieldnames=['speed_m_s','friction_coefficient','normal_force_n'])
                writer.writeheader()
                for force in (40,100,160):writer.writerow({'speed_m_s':1e-6,'friction_coefficient':.1,'normal_force_n':force})
            p['shear_strength_fit']['experiment_file']=str(path)
            result=prepare_fit(p,ROOT)
            self.assertEqual(result[4],.1)
            self.assertEqual([e.line_load_n_m for e in result[2]],[400,1000,1600])
            self.assertEqual([q['macro_contact']['line_load_n_m'] for q in result[-1]],[400,1000,1600])
            self.assertTrue(all('normal_force_n' not in q['macro_contact'] for q in result[-1]))
            self.assertIsNone(p['contact_friction']['model_parameters']['constant_shear_strength']['shear_strength_pa'])

    def test_saved_model_signature_uses_actual_active_geometry(self):
        p=parameters();d={'parameters_snapshot':p,'psd_sha256':'same'}
        original=_signature(d)
        p['macro_contact']['domain_half_width_m']=.099
        p['macro_contact']['geometry_parameters']['cylinder_plane']['radius_m']=.8
        self.assertEqual(_signature(d),original)
        p['macro_contact']['geometry_parameters']['finite_cylindrical_pad']['sliding_length_m']=.04
        self.assertNotEqual(_signature(d),original)


if __name__=='__main__':unittest.main()
