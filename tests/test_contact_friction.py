"""Real-area integration, relative-motion direction, dissipation and extensibility."""
from copy import deepcopy
import json
from pathlib import Path
import tempfile
import unittest

import numpy as np

from persson_contact.contact_friction import (ConstantShearStrength,LogGaussianShearStrength,
    CustomShearStrength,ShearState,evaluate_contact_friction,make_shear_strength)
from persson_contact.workflow import run_project

ROOT=Path(__file__).resolve().parents[1]


def sample(law,velocity=.2,**kwargs):
    return evaluate_contact_friction([0,1,3,4],[.25,.5,0],[2,3,0],[.1,.1,.2],velocity,law,**kwargs)


class ContactStrengthTests(unittest.TestCase):
    def test_constant_nonuniform_cell_real_area_and_units(self):
        r=sample(ConstantShearStrength(8),transverse_width_m=.5);s=r.summary
        self.assertAlmostEqual(s['real_contact_area_per_width_m'],1.25)
        self.assertAlmostEqual(s['normal_line_load_n_m'],8)
        self.assertAlmostEqual(s['friction_force_magnitude_per_width_n_m'],10)
        self.assertAlmostEqual(s['kinetic_contact_friction_coefficient'],1.25)
        self.assertAlmostEqual(s['friction_force_magnitude_n'],5)
        self.assertAlmostEqual(s['normal_force_n'],4)
        self.assertAlmostEqual(s['real_contact_area_m2'],.625)
        self.assertAlmostEqual(s['dissipated_power_w'],1)
        np.testing.assert_allclose(r.profile['real_contact_pressure_pa'],[8,6,0])
        np.testing.assert_allclose(r.profile['upper_contact_traction_pa'],[-2,-4,0])

    def test_reversal_action_reaction_and_common_velocity_translation(self):
        law=ConstantShearStrength(8)
        forward=sample(law,.3);backward=sample(law,-.3)
        np.testing.assert_array_equal(forward.profile['upper_contact_traction_pa'],-backward.profile['upper_contact_traction_pa'])
        self.assertEqual(forward.summary['dissipated_power_per_width_w_m'],backward.summary['dissipated_power_per_width_w_m'])
        u0,u1=.1,.4;r=sample(law,u1-u0)
        power=-(r.summary['upper_friction_force_per_width_n_m']*u1+r.summary['lower_friction_force_per_width_n_m']*u0)
        self.assertAlmostEqual(power,r.summary['dissipated_power_per_width_w_m'])
        translated=sample(law,(u1+2)-(u0+2))
        self.assertEqual(translated.summary['friction_force_magnitude_per_width_n_m'],r.summary['friction_force_magnitude_per_width_n_m'])
        self.assertIsNone(sample(law).summary['friction_force_magnitude_n'])

    def test_zero_relative_speed_is_not_zero_static_friction(self):
        r=sample(ConstantShearStrength(8),0);s=r.summary
        self.assertEqual(s['status'],'no_relative_sliding_static_unmodelled')
        self.assertEqual(s['friction_force_magnitude_per_width_n_m'],0)
        self.assertEqual(s['dissipated_power_per_width_w_m'],0)
        self.assertIsNone(s['kinetic_contact_friction_coefficient']);self.assertIsNone(s['static_friction_limit_n_m'])
        np.testing.assert_allclose(r.profile['nominal_sliding_strength_pa'],[2,4,0])
        # Exact zero only; do not suppress a small but nonzero sliding speed.
        self.assertEqual(sample(ConstantShearStrength(8),1e-30).summary['status'],'sliding')
        self.assertEqual(sample(ConstantShearStrength(0),1).summary['kinetic_contact_friction_coefficient'],0)

    def test_log_gaussian_peak_decades_reversal_and_zero_limit(self):
        law=LogGaussianShearStrength(2,10,1e-3,2)
        peak=sample(law,1e-3).profile['real_contact_shear_strength_pa'][0]
        left=sample(law,1e-5).profile['real_contact_shear_strength_pa'][0]
        right=sample(law,1e-1).profile['real_contact_shear_strength_pa'][0]
        self.assertEqual(peak,12);self.assertAlmostEqual(left,2+10*np.exp(-.5));self.assertAlmostEqual(left,right)
        self.assertEqual(sample(law,-1e-3).summary['friction_force_magnitude_per_width_n_m'],sample(law,1e-3).summary['friction_force_magnitude_per_width_n_m'])
        self.assertEqual(sample(law,0).profile['real_contact_shear_strength_pa'][0],2)
        with np.errstate(all='raise'):
            for v in [1e-250,1e250]:self.assertTrue(np.isfinite(sample(law,v).summary['friction_force_magnitude_per_width_n_m']))

    def test_custom_pressure_temperature_interface_and_no_contact_evaluation(self):
        seen=[]
        def callback(state,parameters):
            seen.append(state)
            return parameters['tau0']+.3*state.real_contact_pressure_pa+2*(state.temperature_k-298.15)
        law=CustomShearStrength(callback,{'tau0':2})
        r=sample(law,temperature_k=299.15)
        # Integral a*(4 + .3*p/a) dx = 4*A_real + .3*W.
        self.assertAlmostEqual(r.summary['friction_force_magnitude_per_width_n_m'],4*1.25+.3*8)
        self.assertEqual(len(seen[0].x_m),2)
        self.assertTrue(np.all(seen[0].real_contact_pressure_pa>0))
        r=evaluate_contact_friction([0,1],[0],[0],[1],1,CustomShearStrength(lambda *args:(_ for _ in ()).throw(RuntimeError())))
        self.assertEqual(r.summary['status'],'no_real_contact');self.assertIsNone(r.summary['kinetic_contact_friction_coefficient'])

    def test_invalid_stresses_states_and_callback_shapes(self):
        for tau in [None,True,-1,float('nan'),float('inf')]:
            with self.assertRaises(ValueError):ConstantShearStrength(tau)
        for values in [(0,1,0,1),(0,1,1,0),(0,-1,1,1)]:
            with self.assertRaises(ValueError):LogGaussianShearStrength(*values)
        for callback in [lambda s,p:-1,lambda s,p:float('nan'),lambda s,p:[1,2,3]]:
            with self.assertRaises(ValueError):sample(CustomShearStrength(callback))
        for kwargs in [{'temperature_k':0},{'transverse_width_m':0},{'temperature_k':[298,299]}]:
            with self.assertRaises(ValueError):sample(ConstantShearStrength(1),**kwargs)
        with self.assertRaises(ValueError):evaluate_contact_friction([0,1],[0],[1],[1],1,ConstantShearStrength(1))
        with self.assertRaises(ValueError):evaluate_contact_friction([0,1],[1],[1],[-1],1,ConstantShearStrength(1))

    def test_only_selected_preset_is_read_and_custom_source_is_fresh(self):
        raw={'model':'constant_shear_strength','model_parameters':{'constant_shear_strength':{'shear_strength_pa':8},'custom':None}}
        self.assertEqual(sample(make_shear_strength(raw)).summary['friction_force_magnitude_per_width_n_m'],10)
        with tempfile.TemporaryDirectory() as directory:
            path=Path(directory)/'law.py';path.write_text('def strength(state,parameters):\n    return 4\n')
            raw={'model':'custom','model_parameters':{'custom':{'module_file':'law.py','function':'strength','parameters':{}}}}
            first=make_shear_strength(raw,directory)
            # Same filename/length can share a filesystem timestamp: hash source, no stale pyc.
            path.write_text('def strength(state,parameters):\n    return 8\n')
            second=make_shear_strength(raw,directory)
            self.assertNotEqual(first.metadata['module_sha256'],second.metadata['module_sha256'])
            self.assertEqual(sample(second).summary['friction_force_magnitude_per_width_n_m'],10)
            self.assertEqual(sample(first).summary['friction_force_magnitude_per_width_n_m'],5)


class ContactFrictionWorkflowTests(unittest.TestCase):
    def parameters(self):
        from parameters import PARAMETERS
        p=deepcopy(PARAMETERS);p['psd']['path']=str(ROOT/'examples/data/roughness_psd.dat')
        p['macro_contact']['solver']['cells']=65;p['macro_contact']['contact_table']['points']=601
        p['contact_friction']['model_parameters']['constant_shear_strength']['shear_strength_pa']=1e5
        return p

    def test_auto_normal_dependency_and_unchanged_normal_state(self):
        with tempfile.TemporaryDirectory() as directory:
            p=self.parameters();p['output_directory']='normal'
            run_project(p,base_directory=directory,run_macro_contact=True,save_plots=False)
            p['output_directory']='friction'
            r=run_project(p,base_directory=directory,run_contact=False,run_contact_friction=True,save_plots=False)
            out=Path(directory)/'friction'
            self.assertEqual((out/'macro_geometry.csv').read_bytes(),(Path(directory)/'normal/macro_geometry.csv').read_bytes())
            cases=list(__import__('csv').DictReader((out/'contact_friction_cases.csv').open()))
            self.assertEqual(len(cases),5);self.assertEqual(cases[0]['kinetic_contact_friction_coefficient'],'')
            self.assertEqual(cases[1]['friction_force_magnitude_n'],'')
            self.assertEqual(r['contact_friction']['case_count'],5)
            self.assertFalse((out/'leakage_cases.csv').exists())

    def test_missing_experiment_and_invalid_combinations_fail_before_output(self):
        from parameters import PARAMETERS
        with tempfile.TemporaryDirectory() as directory:
            p=deepcopy(PARAMETERS);p['output_directory']='bad'
            with self.assertRaisesRegex(ValueError,'no experimental value'):
                run_project(p,base_directory=directory,run_contact_friction=True,save_plots=False)
            self.assertFalse((Path(directory)/'bad').exists())
            p=self.parameters();p['output_directory']='bad'
            with self.assertRaises(NotImplementedError):run_project(p,base_directory=directory,run_contact_friction=True,run_leakage=True)
            self.assertFalse((Path(directory)/'bad').exists())

    def test_log_model_and_custom_example_integrate_through_parameters(self):
        with tempfile.TemporaryDirectory() as directory:
            p=self.parameters();p['contact_friction']['model']='gaussian_log_velocity'
            p['contact_friction']['model_parameters']['gaussian_log_velocity']={
                'baseline_pa':1e4,'amplitude_pa':1e5,'peak_velocity_m_s':1e-3,'width_decades':1.5}
            p['output_directory']='gaussian';p['contact_friction']['relative_velocities_m_s']=[1e-3]
            run_project(p,base_directory=directory,run_contact_friction=True,save_plots=False)
            p['contact_friction']['model']='custom';p['contact_friction']['model_parameters']['custom']={
                'module_file':str(ROOT/'examples/contact_shear_model.py'),'function':'shear_strength','parameters':{'shear_strength_pa':1.1e5}}
            p['output_directory']='custom'
            run_project(p,base_directory=directory,run_contact_friction=True,save_plots=False)
            a=json.loads((Path(directory)/'gaussian/contact_friction_diagnostics.json').read_text())
            b=json.loads((Path(directory)/'custom/contact_friction_diagnostics.json').read_text())
            self.assertEqual(a['cases'][0]['friction_force_magnitude_per_width_n_m'],b['cases'][0]['friction_force_magnitude_per_width_n_m'])


if __name__=='__main__':unittest.main()
