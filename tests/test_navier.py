"""Independent planar BVP, slip-aware homogenization and real entry paths."""
from copy import deepcopy
import json
from pathlib import Path
import tempfile
from types import SimpleNamespace
import unittest

import numpy as np
from scipy.integrate import quad

from persson_contact.interface import NewtonianNoSlip,NewtonianNavierSlip,make_interface
from persson_contact.flow import homogenize_bins,pressure_flow
from persson_contact.core import PerssonModel
from persson_contact.psd import read_psd
from persson_contact.leakage import prepare_leakage,make_transport
from persson_contact.transport_table import TransportTable
from persson_contact.workflow import run_project


class LocalNavierTests(unittest.TestCase):
    def test_no_slip_equal_slip_and_single_wall_limits(self):
        eta=.003;h=2e-6;g=-3e5
        base=NewtonianNoSlip(eta)
        for b in (0.,1e-8,h,1e3*h):
            fluid=NewtonianNavierSlip(eta,b,b)
            np.testing.assert_allclose(fluid.pressure_mobility(h),h**3/(12*eta)*(1+6*b/h),rtol=3e-15)
            np.testing.assert_allclose(fluid.flux(h,g),-g*fluid.pressure_mobility(h),rtol=1e-14)
        one=NewtonianNavierSlip(eta,1e8*h,0)
        np.testing.assert_allclose(one.pressure_mobility(h)/base.pressure_mobility(h),4.,rtol=1e-8)
        for fluid in (base,NewtonianNavierSlip(eta,0,0)):
            self.assertEqual(fluid.flux(0,g,1,0),0)
            self.assertTrue(np.isnan(fluid.local_response(0,g,1,0).lower_wall_shear_stress_pa))
        np.testing.assert_allclose(base.flux([h,2*h],g,.1,-.5),
            NewtonianNavierSlip(eta,0,0).flux([h,2*h],g,.1,-.5),rtol=1e-15)

    def test_asymmetric_flow_and_tractions_against_independent_boundary_solution(self):
        eta=.04;h=3e-6
        for b0,b1 in ((0.,0.),(0.,h),(2*h,0.),(.1*h,3*h)):
            for g in (-2e6,0.,3e6):
                u0,u1=.01,-.02
                law=NewtonianNavierSlip(eta,b0,b1)
                v=law.local_response(h,g,u0,u1)
                a=g/(2*eta)
                c0,c1=np.linalg.solve([[1.,-b0],[1.,h+b1]],
                    [u0,u1-a*(h*h+2*b1*h)])
                velocity=lambda z:a*z*z+c1*z+c0
                expected=quad(velocity,0,h,epsabs=1e-30)[0]
                np.testing.assert_allclose(v.flux_per_width_m2_s,expected,rtol=3e-14)
                np.testing.assert_allclose([v.lower_fluid_velocity_m_s,v.upper_fluid_velocity_m_s],
                    [velocity(0),velocity(h)],rtol=1e-13,atol=1e-17)
                np.testing.assert_allclose([v.lower_wall_shear_stress_pa,v.upper_wall_shear_stress_pa],
                    [eta*c1,-eta*(2*a*h+c1)],rtol=1e-13,atol=1e-12)
                np.testing.assert_allclose(v.lower_wall_shear_stress_pa+v.upper_wall_shear_stress_pa,
                    -g*h,rtol=1e-12,atol=1e-12)
                dissipation=eta*quad(lambda z:(2*a*z+c1)**2,0,h)[0]+eta*b0*c1*c1+eta*b1*(2*a*h+c1)**2
                power=-g*expected-v.lower_wall_shear_stress_pa*u0-v.upper_wall_shear_stress_pa*u1
                np.testing.assert_allclose(power,dissipation,rtol=1e-13)

    def test_wall_exchange_and_galilean_shift(self):
        h=1e-6;g=3e6;u0=.003;u1=-.008
        law=NewtonianNavierSlip(.02,h,3*h)
        a=law.local_response(h,g,u0,u1)
        b=NewtonianNavierSlip(.02,3*h,h).local_response(h,g,u1,u0)
        np.testing.assert_allclose(a.flux_per_width_m2_s,b.flux_per_width_m2_s,rtol=1e-14)
        np.testing.assert_allclose(a.lower_wall_shear_stress_pa,b.upper_wall_shear_stress_pa,rtol=1e-14)
        shifted=law.local_response(h,g,u0+.1,u1+.1)
        np.testing.assert_allclose(shifted.flux_per_width_m2_s-a.flux_per_width_m2_s,h*.1,rtol=1e-14)
        np.testing.assert_allclose(shifted.lower_wall_shear_stress_pa,a.lower_wall_shear_stress_pa,rtol=1e-13)

    def test_reserved_constitutive_inputs_and_validation(self):
        law=NewtonianNavierSlip(.02,1e-8,0)
        np.testing.assert_array_equal(law.viscosity([0,1e5],[0,1e8],[200,500]),[.02,.02])
        a=law.local_response(1e-6,2e5,fluid_pressure_pa=0,temperature_k=200)
        b=law.local_response(1e-6,2e5,fluid_pressure_pa=1e8,temperature_k=500)
        self.assertEqual(a.flux_per_width_m2_s,b.flux_per_width_m2_s)
        for value in (-1,None,True,np.nan,np.inf):
            with self.assertRaises(ValueError):NewtonianNavierSlip(.02,value,0)
        for kwargs in ({'shear_rate_s_inv':-1},{'fluid_pressure_pa':-1},{'temperature_k':0}):
            with self.assertRaises(ValueError):law.viscosity(**kwargs)

    def test_independent_selectors_read_only_active_presets(self):
        cfg={'rheology':'newtonian','wall_slip':'no_slip',
             'rheology_parameters':{'newtonian':{'viscosity_pa_s':.02},'future':'unused'},
             'slip_parameters':{'no_slip':{},'navier':{'bad':'unused'}}}
        self.assertIsInstance(make_interface(cfg),NewtonianNoSlip)
        cfg['wall_slip']='navier'
        with self.assertRaises(ValueError):make_interface(cfg)
        cfg['slip_parameters']['navier']={'lower_slip_length_m':1e-8,'upper_slip_length_m':0}
        self.assertTrue(make_interface(cfg).has_slip)
        with self.assertRaises(ValueError):make_interface({**cfg,'viscosity_pa_s':.1})
        self.assertEqual(make_interface({'viscosity_pa_s':.02}).viscosity_pa_s,.02)


class SlipTransportTests(unittest.TestCase):
    def test_uniform_gap_contact_mixture_preserves_absolute_mobility(self):
        u=2e-6;law=NewtonianNavierSlip(.02,u,3*u)
        # One nonzero probability bin has log-midpoint u; reference h differs.
        edges=np.array([0.,u/2,2*u]);cdf=np.array([0.,0.,.8])
        result=homogenize_bins(edges,cdf,.2,u/3,law,.5)
        expected=float(law.geometric_mobility(u))*.6
        np.testing.assert_allclose(result['geometric'],expected,rtol=1e-12)
        np.testing.assert_allclose(result['factor'],expected/(u/3)**3,rtol=1e-12)

    def test_low_pressure_continuation_uses_slip_law(self):
        law=NewtonianNavierSlip(.02,1e-6,2e-6)
        class Model:
            def area(self,p):return np.asarray(p)/1e6
        p=np.array([1e3,1e4,1e5]);h=np.array([3e-6,2e-6,1e-6])
        k=law.geometric_mobility(h)*.5
        table=TransportTable(p,h,k,k,Model(),law)
        low=10.;hlow=table.mean_gap(low)
        expected=k[0]*law.geometric_mobility(hlow)/law.geometric_mobility(h[0])
        np.testing.assert_allclose(table.geometric_conductivity(low),expected,rtol=1e-14)
        self.assertTrue(np.isinf(table.geometric_conductivity(0)))

    def test_critical_junction_uses_selected_wall_law(self):
        law=NewtonianNavierSlip(.02,1e-6,2e-6)
        class Model:
            def pressure_for_relative_contact_area(self,a):return 2e5
            def critical_separation(self,p,a,tol):return SimpleNamespace(status='found',critical_separation_m=1e-6)
        raw={'mode':'spatial','fluid':'liquid','geometry':{'type':'square_ring',
             'length_m':.001,'reference_contact_pressure_pa':1e5},'seal_width_m':.01,
             'inlet_pressure_pa':1e4,'conductance_model':'critical_junction','percolation_policy':'physical_zero'}
        plan=prepare_leakage(raw,'.',4e6)
        rows=[{'pressure_pa':p,'mean_gap_m':h} for p,h in [(1e3,3e-6),(1e4,2e-6),(1e5,1e-6)]]
        table=make_transport(plan,Model(),law,rows,.42)
        np.testing.assert_allclose(table.k,np.full(3,law.geometric_mobility(1e-6)),rtol=1e-14)

    def test_navier_workflow_and_gas_estimate_guard(self):
        from parameters import PARAMETERS
        root=Path(__file__).resolve().parents[1]
        cfg=deepcopy(PARAMETERS);cfg['psd']['path']=str(root/'examples/data/roughness_psd.dat')
        cfg.pop('contact_pressure');cfg['pressures_pa']=[1e4,5e4,1e5]
        cfg['leakage'].update(geometry='square_ring',inlet_pressures_pa=[0.,1e4],coupling='fixed_contact')
        cfg['leakage']['geometry_parameters']['square_ring']['reference_contact_pressure_pa']=5e4
        cfg['interface']['wall_slip']='navier'
        with tempfile.TemporaryDirectory() as directory:
            cfg['output_directory']='rejected'
            with self.assertRaisesRegex(ValueError,'compute_equivalent_gas=False'):
                run_project(cfg,base_directory=directory,run_leakage=True,compute_gap_pdf=False,save_plots=False)
            self.assertFalse((Path(directory)/'rejected').exists())
            cfg['leakage']['compute_equivalent_gas']=False
            cfg['output_directory']='slip'
            result=run_project(cfg,base_directory=directory,run_leakage=True,compute_gap_pdf=False,save_plots=False)
            self.assertEqual(result['workflow_status'],'complete')
            diag=json.loads((Path(directory)/'slip/flow_diagnostics.json').read_text())
            self.assertEqual(diag['interface']['wall_slip'],'navier')
            cases=json.loads((Path(directory)/'slip/leakage_diagnostics.json').read_text())['cases']
            self.assertEqual(cases[0]['flow_rate_m3_s'],0.)
            law=make_interface(cfg['interface'])
            model=PerssonModel(read_psd(cfg['psd']['path']),4e6,.4,
                               pdf_scale_method=cfg['options']['pdf_scale_method'])
            k=pressure_flow(model,5e4,law).geometric_conductivity_m3
            expected=.01/.001*k/(12*.1)*1e4
            np.testing.assert_allclose(cases[1]['flow_rate_m3_s'],expected,rtol=1e-9)
            self.assertGreater(k,pressure_flow(model,5e4,NewtonianNoSlip(.1)).geometric_conductivity_m3)


if __name__=='__main__':unittest.main()
