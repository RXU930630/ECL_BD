"""Shear-factor conventions and independent conservation/analytic limits."""
from copy import deepcopy
import json
from pathlib import Path
import tempfile
import unittest

import numpy as np

from persson_contact.shear import ShearOptions,ShearFlowLaw,shear_length,read_shear_options
from persson_contact.moving_flow import solve_fixed_contact_motion
from persson_contact.workflow import run_project


class ShearFactorTests(unittest.TestCase):
    def test_original_fortran_shear_kernel_reference(self):
        fixture=json.loads((Path(__file__).parent/'data/fortran_shear_flow.json').read_text())
        self.assertEqual(len(fixture['cases']),13)
        for c in fixture['cases']:
            phi=float(shear_length(c['mean_gap_m'],c['rms_height_m'],c['percolation_gap_m']))/c['rms_height_m']
            if c['original_cutoff_only']:
                np.testing.assert_allclose(phi,c['python_shear_flow_factor'],rtol=1e-14)
                self.assertLess(abs(phi/c['fortran_shear_flow_factor']-1),5e-6)
            else:
                np.testing.assert_allclose(phi,c['fortran_shear_flow_factor'],rtol=1e-14)

    def test_smooth_blocked_and_large_gap_limits(self):
        h=np.array([0.,.1,.5,1.,10.])
        np.testing.assert_array_equal(shear_length(h,0.,0.),np.zeros_like(h))
        np.testing.assert_allclose(shear_length(h,1.,.5)[:3],h[:3])
        self.assertTrue(np.all(shear_length(h,1.,.5)<=h))
        large=float(shear_length(1e4,1.,0.))
        np.testing.assert_allclose(large,3/(2e4),rtol=2e-8)
        law=ShearFlowLaw(ShearOptions(rough_surface='lower',upper_wall_velocity_m_s=2.),1.,.5)
        self.assertEqual(float(law.source(.1)),0.) # stationary rough wall, below percolation gap
        moving_rough=ShearFlowLaw(ShearOptions(rough_surface='upper',upper_wall_velocity_m_s=2.),1.,.5)
        self.assertEqual(float(moving_rough.source(.1)),.2)

    def test_wall_ownership_galilean_and_velocity_reversal(self):
        h=np.array([.1,1.,10.]);rms=1.;hc=.5
        lower=ShearFlowLaw(ShearOptions(rough_surface='lower',lower_wall_velocity_m_s=.3,
                                     upper_wall_velocity_m_s=-.2),rms,hc)
        upper=ShearFlowLaw(ShearOptions(rough_surface='upper',lower_wall_velocity_m_s=-.2,
                                     upper_wall_velocity_m_s=.3),rms,hc)
        np.testing.assert_allclose(lower.source(h),upper.source(h),rtol=1e-14)
        c0,c1=lower.coefficients(h);np.testing.assert_allclose(c0+c1,h,rtol=1e-14)
        shifted=ShearFlowLaw(ShearOptions(lower_wall_velocity_m_s=.8,upper_wall_velocity_m_s=.3),rms,hc)
        np.testing.assert_allclose(shifted.source(h)-lower.source(h),.5*h,rtol=1e-14)
        reversed_law=ShearFlowLaw(ShearOptions(lower_wall_velocity_m_s=-.3,upper_wall_velocity_m_s=.2),rms,hc)
        np.testing.assert_allclose(reversed_law.source(h),-lower.source(h),rtol=1e-14)
        common=ShearFlowLaw(ShearOptions(lower_wall_velocity_m_s=.3,upper_wall_velocity_m_s=.3),rms,hc)
        np.testing.assert_allclose(common.source(h),.3*h,rtol=1e-14)

    def test_invalid_configuration_and_inputs(self):
        for changes in ({'rough_surface':'both'},{'upper_wall_velocity_m_s':True},
                        {'lower_wall_velocity_m_s':np.inf},{'pressure_floor_pa':np.nan}):
            with self.assertRaises(ValueError):ShearOptions(**changes)
        with self.assertRaises(NotImplementedError):ShearOptions(model='invented')
        with self.assertRaises(ValueError):read_shear_options({'unknown':1})
        with self.assertRaises(ValueError):shear_length(-1.,1.,0.)
        with self.assertRaises(ValueError):shear_length(1.,-1.,0.)


class MovingConservationTests(unittest.TestCase):
    def test_uniform_gap_pressure_motion_zero_net_and_reversal(self):
        r=lambda x:np.full_like(x,2.)
        s=lambda x:np.full_like(x,.5)
        for pin,pout in ((1.,0.),(0.,0.),(0.,1.),(0.,2.)):
            v=solve_fixed_contact_motion(r,s,0.,1.,pin,pout)
            expected=(pin-pout)/2+.5
            np.testing.assert_allclose(v.flux_per_width,expected,rtol=1e-12,atol=1e-14)
            np.testing.assert_allclose(v.fluid_pressure_pa,pin+(pout-pin)*v.x_m,rtol=1e-12,atol=1e-14)
        with self.assertRaises(ValueError):solve_fixed_contact_motion(lambda x:np.full_like(x,np.inf),s,0,1,0,0)

    def test_resistance_weighted_source_against_exact_integrals(self):
        r=lambda x:1+x;s=lambda x:1-x
        v=solve_fixed_contact_motion(r,s,0.,1.,2.,1.,tolerance=1e-7)
        expected=(1+2/3)/1.5
        np.testing.assert_allclose(v.flux_per_width,expected,rtol=1e-6)
        x=v.x_m;p=2+(1-expected)*x-expected*x*x/2-x**3/3
        np.testing.assert_allclose(v.fluid_pressure_pa,p,atol=3e-7,rtol=0)
        self.assertLess(v.face_flux_relative_spread,1e-9)

    def test_equal_boundary_pressures_allow_internal_extremum(self):
        v=solve_fixed_contact_motion(lambda x:np.ones_like(x),lambda x:1-2*x,
                                    0.,1.,0.,0.,tolerance=1e-7)
        self.assertAlmostEqual(v.flux_per_width,0.,places=14)
        np.testing.assert_allclose(v.fluid_pressure_pa,v.x_m-v.x_m*v.x_m,atol=2e-7,rtol=0)
        self.assertGreater(np.max(v.fluid_pressure_pa),.249)
        with self.assertRaisesRegex(ValueError,'cavitation'):
            solve_fixed_contact_motion(lambda x:np.ones_like(x),lambda x:2*x-1,0,1,0,0)


class ShearWorkflowTests(unittest.TestCase):
    def config(self):
        from parameters import PARAMETERS
        root=Path(__file__).resolve().parents[1]
        cfg=deepcopy(PARAMETERS);cfg['psd']['path']=str(root/'examples/data/roughness_psd.dat')
        cfg.pop('contact_pressure');cfg['pressures_pa']=[1e4,5e4,1e5]
        cfg['leakage'].update(geometry='square_ring',coupling='fixed_contact',
                              inlet_pressures_pa=[0.,1e4],compute_equivalent_gas=False)
        cfg['leakage']['geometry_parameters']['square_ring']['reference_contact_pressure_pa']=5e4
        cfg['output_directory']='moving'
        return cfg

    def test_table_only_and_spatial_combined_flow(self):
        cfg=self.config()
        with tempfile.TemporaryDirectory() as directory:
            result=run_project(cfg,base_directory=directory,run_leakage=True,
                compute_shear_flow=True,compute_gap_pdf=False,save_plots=False)
            self.assertEqual(result['workflow_status'],'complete')
            out=Path(directory)/'moving'
            table=np.genfromtxt(out/'shear_flow_table.csv',delimiter=',',names=True)
            diag=json.loads((out/'leakage_diagnostics.json').read_text())
            self.assertTrue(diag['motion_transport_enabled'])
            self.assertGreater(diag['cases'][0]['flow_rate_m3_s'],0.)
            middle=table[1];width=cfg['leakage']['seal_width_m'];length=.001
            expected=middle['motion_flux_per_width_m2_s']*width
            np.testing.assert_allclose(diag['cases'][0]['flow_rate_m3_s'],expected,rtol=1e-12)
            self.assertIsNone(diag['cases'][0]['pdf_bound_min_flow_m3_s'])
            profile=np.genfromtxt(out/'leakage_profile_001.csv',delimiter=',',names=True)
            np.testing.assert_allclose(profile['combined_flux_per_width_m2_s'],
                profile['motion_flux_per_width_m2_s']+profile['pressure_flux_per_width_m2_s'],rtol=1e-13)
            self.assertGreater(diag['cases'][1]['flow_rate_m3_s'],expected)
            cfg['output_directory']='table_only'
            standalone=run_project(cfg,base_directory=directory,run_contact=False,run_leakage=False,
                compute_shear_flow=True,compute_gap_pdf=False,save_plots=False)
            self.assertIn('shear_flow',standalone)
            self.assertFalse((Path(directory)/'table_only/leakage_cases.csv').exists())

    def test_unsupported_combinations_fail_before_output(self):
        for change in ('slip','gas','unloading','estimates','both_rough'):
            cfg=self.config()
            if change=='slip':cfg['interface']['wall_slip']='navier'
            if change=='gas':cfg['leakage']['fluid']='gas'
            if change=='unloading':cfg['leakage']['coupling']='local_unloading'
            if change=='estimates':cfg['leakage']['compute_equivalent_gas']=True
            if change=='both_rough':cfg['shear_flow']['rough_surface']='both'
            with self.subTest(change=change),tempfile.TemporaryDirectory() as d:
                with self.assertRaises((ValueError,NotImplementedError)):
                    run_project(cfg,base_directory=d,run_leakage=True,compute_shear_flow=True,
                                compute_gap_pdf=False,save_plots=False)
                self.assertFalse((Path(d)/'moving').exists())

    def test_zero_velocity_preserves_pressure_only_result(self):
        cfg=self.config();cfg['shear_flow']['upper_wall_velocity_m_s']=0.
        cfg['leakage']['coupling']='local_unloading'
        with tempfile.TemporaryDirectory() as d:
            rates=[]
            for enabled in (False,True):
                cfg['output_directory']=str(enabled)
                run_project(cfg,base_directory=d,run_leakage=True,compute_shear_flow=enabled,
                            compute_gap_pdf=False,save_plots=False)
                saved=json.loads((Path(d)/str(enabled)/'leakage_diagnostics.json').read_text())
                rates.append([row['flow_rate_m3_s'] for row in saved['cases']])
            np.testing.assert_array_equal(rates[0],rates[1])


if __name__=='__main__':unittest.main()
