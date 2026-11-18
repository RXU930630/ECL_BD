"""Check scientific diagnostic decisions and preservation of physical geometry."""
from copy import deepcopy
from pathlib import Path
import tempfile
import unittest

from persson_contact.mixed_checks import check_settings,compare_cases,run_checks
from persson_contact.shear_fit import Experiment,fit_constant_strength


def case(speed=1e-4):
    return {'lower_wall_velocity_m_s':speed,'real_contact_area_per_width_m':1e-5,
        'minimum_mean_gap_m':1e-6,'maximum_fluid_pressure_pa':1e4,'fluid_load_fraction':.5,
        'lower_resisting_friction_coefficient':.12,'lower_resisting_solid_coefficient':.1,
        'lower_resisting_fluid_sliding_coefficient':.03,'lower_resisting_fluid_pressure_coefficient':-.01,
        'edge_pressure_ratio':1e-6}


class NumericalAuditTests(unittest.TestCase):
    def setup_audit(self,root,geometry='finite_cylindrical_pad'):
        p={'output_directory':str(root),'macro_contact':{'geometry':geometry,'domain_half_width_m':.015,
            'geometry_parameters':{'finite_cylindrical_pad':{'sliding_length_m':.03}},
            'solver':{'cells':129}},'mixed_friction':{'lower_wall_velocities_m_s':[1e-4]},
            'options':{'n_q':128,'pdf_log_step':.03}}
        summary={'cases':[case()],'target_line_load_n_m':100.,'far_field_edge_pressure_check_applied':geometry=='cylinder_plane'}
        return p,summary

    def test_physical_pad_does_not_grow_and_inputs_are_preserved(self):
        with tempfile.TemporaryDirectory() as directory:
            root=Path(directory);p,summary=self.setup_audit(root);original=deepcopy(p);calls=[]
            def runner(config,**kwargs):
                calls.append(deepcopy(config));return {'cases':[case()]}
            result=run_checks(p,root,summary,{},check_settings({'enabled':True}),runner,root)
            self.assertEqual(p,original)
            self.assertEqual(len(calls),2)
            self.assertTrue(all(c['macro_contact']['geometry_parameters']==p['macro_contact']['geometry_parameters'] for c in calls))
            self.assertTrue(all(c['macro_contact']['domain_half_width_m']==.015 for c in calls))
            self.assertEqual(result['probes'][-1]['status'],'not_applicable')
            self.assertEqual(result['status'],'passed_at_selected_cases')

    def test_small_edge_pressure_does_not_hide_domain_sensitivity(self):
        with tempfile.TemporaryDirectory() as directory:
            root=Path(directory);p,summary=self.setup_audit(root,'cylinder_plane')
            def runner(config,**kwargs):
                row=case()
                if config['macro_contact']['domain_half_width_m']>.015:
                    row['lower_resisting_friction_coefficient']*=1.1
                return {'cases':[row]}
            result=run_checks(p,root,summary,{},check_settings({'enabled':True}),runner,root)
            self.assertEqual(result['status'],'needs_attention')
            domain=[r for r in result['comparisons'] if r['check']=='domain' and r['quantity']=='lower_resisting_friction_coefficient']
            self.assertFalse(domain[0]['within_tolerance'])

    def test_resource_cap_and_failed_refinement_never_count_as_a_pass(self):
        for fail in (False,True):
            with self.subTest(fail=fail),tempfile.TemporaryDirectory() as directory:
                root=Path(directory);p,summary=self.setup_audit(root)
                def runner(config,**kwargs):
                    if fail:raise ValueError('disconnected constitutive branch')
                    return {'cases':[case()]}
                settings=check_settings({'enabled':True,'maximum_cells':129})
                result=run_checks(p,root,summary,{},settings,runner,root)
                self.assertEqual(result['status'],'needs_attention' if fail else 'incomplete')

    def test_absolute_limits_handle_tiny_contributions_and_keep_signed_terms(self):
        a=case();b=case();a['lower_resisting_solid_coefficient']=0.;b['lower_resisting_solid_coefficient']=1e-8
        rows=compare_cases([a],[b],check_settings({}),100.)
        self.assertTrue(all(r['within_tolerance'] for r in rows))
        self.assertTrue(any(r['baseline_value']<0 for r in rows))
        with self.assertRaises(ValueError):check_settings({'enabled':'True'})


class FitConsistencyTests(unittest.TestCase):
    def test_incompatible_load_groups_are_flagged_without_changing_global_fit(self):
        experiments=[];sensitivities=[]
        for load,tau in [(100.,1e5),(200.,5e5)]:
            for index,c in enumerate([1e-6,2e-6]):
                experiments.append(Experiment(len(experiments)+2,10.**(-4+index),load,.01+tau*c))
                sensitivities.append(c)
        summary,records=fit_constant_strength(experiments,sensitivities,[.01]*4)
        self.assertAlmostEqual(summary['fitted_shear_strength_pa'],3e5)
        self.assertEqual(summary['fit_consistency_status'],'needs_attention')
        self.assertAlmostEqual(summary['separate_tau_max_min_ratio'],5.)
        self.assertEqual([round(g['separately_fitted_strength_pa']) for g in summary['per_load_diagnostics']],[100000,500000])
        self.assertAlmostEqual(sum(r['fit_leverage'] for r in records),1.)

    def test_exact_fit_is_not_labelled_physically_validated(self):
        e=[Experiment(2,.001,100.,.21),Experiment(3,.01,100.,.41)]
        summary,_=fit_constant_strength(e,[1e-6,2e-6],[.01,.01])
        self.assertEqual(summary['fit_consistency_status'],'no_diagnostic_threshold_exceeded')
        self.assertIn('not a physically validated',summary['interpretation'])


if __name__=='__main__':unittest.main()
