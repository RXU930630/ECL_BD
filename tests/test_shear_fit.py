"""Known-strength recovery, load normalization, input units, cached/simulated I/O."""
from copy import deepcopy
import csv
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

import numpy as np

from persson_contact.shear_fit import Experiment,read_experiments,fit_constant_strength
from persson_contact.shear_fit_output import read_forward_run
from persson_contact.workflow import run_project


def write_csv(path,rows):
    with Path(path).open('w',encoding='utf-8-sig',newline='') as stream:
        w=csv.DictWriter(stream,fieldnames=list(rows[0]),lineterminator='\n');w.writeheader();w.writerows(rows)


def make_forward(path,load=100.,speeds=(1e-4,1e-3,1e-2),direction=1,material='same'):
    path=Path(path);path.mkdir(parents=True)
    cases=[]
    for i,speed in enumerate(speeds):
        area=2e-4/(i+1);mu_fluid=.02*(i+1)
        cases.append({'lower_wall_velocity_m_s':direction*speed,'real_contact_area_per_width_m':area,
            'lower_fluid_force_n_m':-direction*mu_fluid*load,
            'contact_line_load_n_m':load*.2,'fluid_line_load_n_m':load*.8,
            'load_relative_error':0.,'geometry_relative_residual':0.,'mass_relative_residual':0.})
    diag={'workflow_status':'complete','synthetic_demo':True,'target_line_load_n_m':load,'cases':cases,
        'psd_sha256':'synthetic','parameters_snapshot':{'material':{'synthetic':material}},
        'limitations':['Synthetic inverse-problem fixture, no physical solve.']}
    write_csv(path/'mixed_cases.csv',cases)
    (path/'mixed_diagnostics.json').write_text(json.dumps(diag))
    return cases


class ConstantFitTests(unittest.TestCase):
    def test_exact_strength_background_and_total_load_normalization(self):
        tau=7e5;area=np.array([2e-4,1e-4,4e-5]);load=np.array([100.,200.,50.])
        fluid=np.array([.01,.05,.1]);background=np.array([.02,.01,.04]);c=area/load
        rows=[Experiment(i+2,10.**(-4+i),w,ff+bg+tau*ci,bg) for i,(w,ff,bg,ci) in enumerate(zip(load,fluid,background,c))]
        summary,records=fit_constant_strength(rows,c,fluid)
        self.assertAlmostEqual(summary['fitted_shear_strength_pa']/tau,1.,places=14)
        np.testing.assert_allclose([r['residual_experiment_minus_model'] for r in records],0.,atol=1e-15)

    def test_reference_invariance_including_sensitivity_selection(self):
        c=np.array([2e-8,3e-11,2e-12]);fluid=np.array([.01,.02,.03]);tau=3e5
        rows=[Experiment(i+2,.001*(i+1),100.,f+tau*ci) for i,(ci,f) in enumerate(zip(c,fluid))]
        for ref in (1e2,1e6,1e10):
            s,r=fit_constant_strength(rows,c,fluid,reference_strength_pa=ref)
            self.assertAlmostEqual(s['fitted_shear_strength_pa']/tau,1.,places=12)
            self.assertEqual([v['used_for_fit'] for v in r],[True,True,False])

    def test_weighted_solution_and_conditional_standard_error(self):
        c=np.array([1e-6,2e-6,3e-6]);sigma=np.array([.01,.02,.01]);fluid=np.array([.02,.03,.01])
        y=np.array([.31,.58,.95])
        rows=[Experiment(i+2,.001*(i+1),100.,f+v,standard_uncertainty_mu=e) for i,(v,f,e) in enumerate(zip(y,fluid,sigma))]
        expected=np.sum(c*y/sigma**2)/np.sum(c*c/sigma**2)
        error=1/np.sqrt(np.sum(c*c/sigma**2))
        s,_=fit_constant_strength(rows,c,fluid,weighting='standard_uncertainty_mu')
        self.assertAlmostEqual(s['fitted_shear_strength_pa']/expected,1.,places=14)
        self.assertAlmostEqual(s['conditional_standard_error_pa']/error,1.,places=14)

    def test_held_out_outlier_does_not_change_fit(self):
        rows=[Experiment(2,.001,100.,.35),Experiment(3,.002,100.,.65),Experiment(4,.003,100.,100.,use_for_fit=False)]
        s,r=fit_constant_strength(rows,[1e-6,2e-6,3e-6],[.05,.05,.05])
        self.assertAlmostEqual(s['fitted_shear_strength_pa'],3e5)
        self.assertFalse(r[-1]['used_for_fit']);self.assertGreater(s['held_out_or_excluded_rmse_mu'],90.)

    def test_nonnegative_boundary_and_no_identifiable_contact(self):
        s,_=fit_constant_strength([Experiment(2,.001,100.,.01)],[1e-6],[.1])
        self.assertEqual(s['fitted_shear_strength_pa'],0.);self.assertTrue(s['nonnegative_boundary_active'])
        self.assertIsNone(s['conditional_standard_error_pa'])
        for c in ([0.],[1e-15]):
            with self.assertRaisesRegex(ValueError,'没有可拟合'):fit_constant_strength([Experiment(2,.001,100.,.1)],c,[.01])

    def test_missing_uncertainty_speed_filter_and_invalid_values(self):
        with self.assertRaisesRegex(ValueError,'standard_uncertainty_mu'):
            fit_constant_strength([Experiment(2,.001,100.,.1)],[1e-6],[.01],weighting='standard_uncertainty_mu')
        rows=[Experiment(2,.001,100.,.4),Experiment(3,.01,100.,5.)]
        s,r=fit_constant_strength(rows,[1e-6,1e-6],[.1,.1],speed_range_m_s=[None,.002])
        self.assertAlmostEqual(s['fitted_shear_strength_pa'],3e5);self.assertFalse(r[1]['used_for_fit'])
        for changes in ({'speed_m_s':0.},{'line_load_n_m':-1.},{'friction_coefficient':float('nan')},{'use_for_fit':'false'}):
            kw={'row_number':2,'speed_m_s':.001,'line_load_n_m':100.,'friction_coefficient':.1,**changes}
            with self.assertRaises(ValueError):Experiment(**kw)


class ExperimentCSVTests(unittest.TestCase):
    def test_force_width_units_bom_and_validation_rows(self):
        with tempfile.TemporaryDirectory() as d:
            path=Path(d)/'e.csv'
            write_csv(path,[{'speed_m_s':.001,'friction_coefficient':.2,'normal_force_n':2.,'use_for_fit':'false'}])
            rows=read_experiments(path,default_line_load_n_m=100.,transverse_width_m=.01)
            self.assertEqual(rows[0].line_load_n_m,200.);self.assertFalse(rows[0].use_for_fit)
            with self.assertRaisesRegex(ValueError,'transverse_width_m'):read_experiments(path,default_line_load_n_m=100.)
            write_csv(path,[{'speed_m_s':.001,'friction_coefficient':.2,'normal_force_n':2.,'line_load_n_m':100.}])
            with self.assertRaisesRegex(ValueError,'不一致'):read_experiments(path,default_line_load_n_m=100.,transverse_width_m=.01)

    def test_empty_duplicate_unknown_and_ragged_csv(self):
        texts=['speed_m_s,friction_coefficient\n','speed_m_s,speed_m_s,friction_coefficient\n1,1,.1\n',
            'speed_m_s,friction_coefficient,unknown\n1,.1,x\n','speed_m_s,friction_coefficient\n1,.1,extra\n',
            'speed_m_s,friction_coefficient\nNaN,.1\n']
        with tempfile.TemporaryDirectory() as d:
            path=Path(d)/'e.csv'
            for text in texts:
                path.write_text(text)
                with self.subTest(text=text),self.assertRaises(ValueError):read_experiments(path,default_line_load_n_m=100.)


class FitWorkflowTests(unittest.TestCase):
    def configure(self,d,loads=(100.,),direction=1):
        from parameters import PARAMETERS
        p=deepcopy(PARAMETERS);p['output_directory']='fit'
        p['contact_friction']['transverse_width_m']=.01
        directories=[];experiments=[]
        for i,load in enumerate(loads):
            directory=Path(d)/f'forward_{i}';cases=make_forward(directory,load,direction=direction);directories.append(str(directory))
            for j,c in enumerate(cases):
                experiments.append({'speed_m_s':abs(c['lower_wall_velocity_m_s']),'line_load_n_m':load,
                    'friction_coefficient':.02*(j+1)+4e5*c['real_contact_area_per_width_m']/load,
                    'use_for_fit':j<2})
        write_csv(Path(d)/'experiment.csv',experiments)
        p['shear_strength_fit'].update(source='saved_forward',experiment_file='experiment.csv',sliding_direction=direction)
        p['shear_strength_fit']['source_parameters']['saved_forward']['directories']=directories
        return p

    def test_saved_multiple_loads_reversal_and_reusable_exported_experiment(self):
        with tempfile.TemporaryDirectory() as d:
            p=self.configure(d,(100.,200.),direction=-1)
            # Saved mode must not read present-day PSD/material or known tau.
            p['psd']['path']='does_not_exist.dat'
            result=run_project(p,base_directory=d,run_contact=False,run_shear_strength_fit=True,save_plots=False)
            self.assertAlmostEqual(result['fit']['fitted_shear_strength_pa'],4e5)
            self.assertEqual(result['fit']['training_points'],4);self.assertTrue(result['synthetic_demo'])
            rows=read_experiments(Path(d)/'fit/experiment_used.csv',default_line_load_n_m=100.,transverse_width_m=.01)
            self.assertEqual(len(rows),6)
            self.assertIsNone(p['contact_friction']['model_parameters']['constant_shear_strength']['shear_strength_pa'])

    def test_stale_cache_different_physics_missing_state_and_nonempty_output(self):
        for mode in ('stale','physics','psd_scale','missing','nonempty'):
            with self.subTest(mode=mode),tempfile.TemporaryDirectory() as d:
                p=self.configure(d,(100.,200.))
                if mode=='stale':
                    path=Path(d)/'forward_0/mixed_cases.csv';path.write_text(path.read_text().replace('-2.0,','-3.0,'))
                if mode=='physics':
                    path=Path(d)/'forward_1/mixed_diagnostics.json';v=json.loads(path.read_text());v['parameters_snapshot']['material']={'different':1};path.write_text(json.dumps(v))
                if mode=='psd_scale':
                    path=Path(d)/'forward_1/mixed_diagnostics.json';v=json.loads(path.read_text());v['parameters_snapshot']['psd']={'scale_c':2.};path.write_text(json.dumps(v))
                if mode=='missing':p['shear_strength_fit']['sliding_direction']=-1
                if mode=='nonempty':
                    out=Path(d)/'fit';out.mkdir();(out/'preserve.txt').write_text('keep')
                with self.assertRaises((ValueError,FileExistsError)):
                    run_project(p,base_directory=d,run_contact=False,run_shear_strength_fit=True,save_plots=False)
                self.assertFalse((Path(d)/'fit/fitted_contact_friction.json').exists())

    def test_simulate_groups_loads_and_overrides_reference_only(self):
        with tempfile.TemporaryDirectory() as d:
            p=self.configure(d,(100.,200.));original=deepcopy(p)
            root=Path(__file__).resolve().parents[1]
            p['psd']['path']=str(root/'examples/data/roughness_psd.dat')
            p['shear_strength_fit']['source']='simulate'
            calls=[]
            def native_fixture(config,*,base_directory,save_plots):
                calls.append(deepcopy(config))
                make_forward(config['output_directory'],config['macro_contact']['line_load_n_m'],
                    config['mixed_friction']['lower_wall_velocities_m_s'])
            with patch('persson_contact.shear_fit_output.run_mixed_project',side_effect=native_fixture):
                result=run_project(p,base_directory=d,run_contact=False,run_shear_strength_fit=True,save_plots=False)
            self.assertEqual(len(calls),2);self.assertAlmostEqual(result['fit']['fitted_shear_strength_pa'],4e5)
            self.assertTrue(all(c['contact_friction']['model_parameters']['constant_shear_strength']['shear_strength_pa']==1e6 for c in calls))
            self.assertEqual(p['contact_friction'],original['contact_friction'])

    def test_failed_forward_preserves_failure_status_no_parameter(self):
        with tempfile.TemporaryDirectory() as d:
            p=self.configure(d)
            p['psd']['path']=str(Path(__file__).resolve().parents[1]/'examples/data/roughness_psd.dat')
            p['shear_strength_fit']['source']='simulate'
            with patch('persson_contact.shear_fit_output.run_mixed_project',side_effect=ValueError('solver failure')):
                with self.assertRaisesRegex(ValueError,'solver failure'):
                    run_project(p,base_directory=d,run_contact=False,run_shear_strength_fit=True,save_plots=False)
            self.assertEqual(json.loads((Path(d)/'fit/fit_diagnostics.json').read_text())['workflow_status'],'forward_failed')
            self.assertFalse((Path(d)/'fit/fitted_contact_friction.json').exists())

    def test_roundoff_duplicate_speeds_and_loads_share_one_state(self):
        with tempfile.TemporaryDirectory() as d:
            p=self.configure(d)
            path=Path(d)/'experiment.csv'
            with path.open(encoding='utf-8-sig') as f:rows=list(csv.DictReader(f))
            extra=dict(rows[0]);extra['speed_m_s']=.0001*(1+1e-14);extra['line_load_n_m']=100.*(1+1e-14)
            rows.append(extra);write_csv(path,rows)
            p['psd']['path']=str(Path(__file__).resolve().parents[1]/'examples/data/roughness_psd.dat')
            p['shear_strength_fit']['source']='simulate';calls=[]
            def forward(config,*,base_directory,save_plots):
                calls.append(config)
                make_forward(config['output_directory'],config['macro_contact']['line_load_n_m'],config['mixed_friction']['lower_wall_velocities_m_s'])
            with patch('persson_contact.shear_fit_output.run_mixed_project',side_effect=forward):
                r=run_project(p,base_directory=d,run_contact=False,run_shear_strength_fit=True,save_plots=False)
            self.assertEqual(len(calls),1);self.assertEqual(len(calls[0]['mixed_friction']['lower_wall_velocities_m_s']),3)
            self.assertAlmostEqual(r['fit']['fitted_shear_strength_pa']/4e5,1.,places=12)


if __name__=='__main__':unittest.main()
