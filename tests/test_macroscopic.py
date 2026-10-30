"""Analytic elasticity/Hertz limits and independent module-2+6 coupling checks."""
from copy import deepcopy
import json
from pathlib import Path
import tempfile
import unittest

import numpy as np
from scipy.integrate import quad

from persson_contact.macroscopic import HalfSpaceLine,HertzLine,solve_smooth_contact
from persson_contact.dry_contact import RoughContactLaw,solve_rough_contact
from persson_contact.core import PerssonModel
from persson_contact.psd import read_psd
from persson_contact.table import build_contact_table
from persson_contact.workflow import run_project

ROOT=Path(__file__).resolve().parents[1]


def table(model,reference,points=601):
    focus=reference.peak_pressure_pa*1e-3
    p=np.unique(np.r_[np.geomspace(1e-100,focus,points//3),np.geomspace(focus,2*reference.peak_pressure_pa,points-points//3+1)])
    return build_contact_table(model,p)


class ElasticLineTests(unittest.TestCase):
    def test_integrated_log_kernel_against_adaptive_load_integral(self):
        el=HalfSpaceLine(.001,33,4e6)
        p=2e5*(1+.2*np.cos(el.x_m/.001*np.pi));u=el.displacement(p)
        for i in [0,7,16,25,32]:
            x=el.x_m[i]
            expected=0.
            for j,(a,b) in enumerate(zip(el.edges_m[:-1],el.edges_m[1:])):
                expected+=p[j]*quad(lambda s:np.log(el.reference_length_m/abs(x-s)),a,b,
                    points=[x] if a<x<b else None,epsabs=1e-15,epsrel=1e-11)[0]
            expected*=2/(np.pi*4e6)
            np.testing.assert_allclose(u[i],expected,rtol=2e-12)
        np.testing.assert_array_equal(el.compliance_m_pa,el.compliance_m_pa.T)
        self.assertGreater(np.linalg.eigvalsh(el.compliance_m_pa).min(),0.)

    def test_nonlocality_linearity_and_reference_change(self):
        el=HalfSpaceLine(.001,33,4e6);p=np.zeros(33);p[16]=1e5
        u=el.displacement(p)
        self.assertTrue(np.all(u>0));np.testing.assert_allclose(u,u[::-1],rtol=2e-14)
        np.testing.assert_allclose(el.displacement(3*p),3*u,rtol=2e-14)
        other=HalfSpaceLine(.001,33,4e6,.008)
        expected=2*el.line_load(p)/(np.pi*4e6)*np.log(4)
        np.testing.assert_allclose(other.displacement(p)-u,expected,rtol=3e-13)

    def test_bad_inputs_and_inconsistent_material(self):
        for kwargs in [(.001,16,4e6),(.001,33,0),(.001,True,4e6),(None,33,4e6)]:
            with self.assertRaises(ValueError):HalfSpaceLine(*kwargs)
        with self.assertRaises(ValueError):HalfSpaceLine(.001,33,4e6,.001)
        el=HalfSpaceLine(.001,33,4e6)
        with self.assertRaises(ValueError):el.displacement([1.,2.])
        with self.assertRaises(ValueError):solve_smooth_contact(el,HertzLine(.005,100,2e6))
        with self.assertRaises(ValueError):solve_smooth_contact(HalfSpaceLine(1e-5,33,4e6),HertzLine(.005,100,4e6))


class SmoothMacroTests(unittest.TestCase):
    def test_hertz_load_complementarity_and_profile(self):
        ref=HertzLine(.005,100,4e6)
        results=[]
        for n in [129,513]:
            el=HalfSpaceLine(3*ref.half_width_m,n,4e6);s=solve_smooth_contact(el,ref)
            error=np.sum(abs(s.pressure_pa-ref.cell_pressure(el.edges_m)))*el.dx_m/100
            results.append(error)
            self.assertLess(s.load_relative_error,1e-12)
            self.assertGreaterEqual(s.pressure_pa.min(),0)
            self.assertLess(s.equilibrium_relative_residual,1e-9)
            np.testing.assert_allclose(s.pressure_pa,s.pressure_pa[::-1],rtol=1e-10,atol=1e-6)
            active=s.pressure_pa>0
            self.assertLess(np.max(abs(s.mean_gap_m[active])),1e-12)
            self.assertGreaterEqual(s.mean_gap_m.min(),-1e-12)
        self.assertLess(results[-1],results[0]);self.assertLess(results[-1],.0003)
        scale=ref.half_width_m**2/(2*ref.radius_m)
        self.assertLess(np.max(abs(s.mean_gap_m-ref.gap(s.x_m)))/scale,.001)

    def test_smooth_pressure_and_gap_are_gauge_invariant(self):
        ref=HertzLine(.005,100,4e6);L=3*ref.half_width_m
        a=solve_smooth_contact(HalfSpaceLine(L,129,4e6),ref)
        b=solve_smooth_contact(HalfSpaceLine(L,129,4e6,4*L),ref)
        np.testing.assert_allclose(a.pressure_pa,b.pressure_pa,rtol=1e-10,atol=1e-6)
        np.testing.assert_allclose(a.mean_gap_m,b.mean_gap_m,rtol=1e-9,atol=1e-16)
        shift=2*100/(np.pi*4e6)*np.log(2)
        np.testing.assert_allclose(b.rigid_offset_m-a.rigid_offset_m,-shift,rtol=1e-12)


class RoughMacroTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.ref=HertzLine(.005,100,4e6)
        cls.model=PerssonModel(read_psd(ROOT/'examples/data/roughness_psd.dat'),4e6,.4)
        cls.table=table(cls.model,cls.ref,1201)

    def test_rough_contact_uses_one_gap_and_balances_load(self):
        el=HalfSpaceLine(3*self.ref.half_width_m,129,4e6)
        s=solve_rough_contact(el,self.ref,RoughContactLaw(self.table))
        self.assertLess(s.load_relative_error,1e-9);self.assertLess(s.equilibrium_relative_residual,1e-9)
        self.assertGreater(s.mean_gap_m.min(),0)
        np.testing.assert_allclose(s.mean_gap_m,s.geometric_profile_m+s.rigid_offset_m+s.displacement_m,atol=1e-13)
        # Independent evaluation of the energy-derived local law, not a second interpolation.
        gap=self.model.mean_gap(s.pressure_pa)
        self.assertLess(np.max(abs(gap-s.mean_gap_m))/self.model.psd.hrms,1e-4)
        np.testing.assert_allclose(s.relative_contact_area,self.model.area(s.pressure_pa),rtol=1e-14)
        np.testing.assert_allclose(s.pressure_pa,s.pressure_pa[::-1],rtol=1e-7,atol=1e-5)

    def test_rough_gauge_and_modulus_load_scaling(self):
        L=3*self.ref.half_width_m;law=RoughContactLaw(self.table)
        a=solve_rough_contact(HalfSpaceLine(L,129,4e6),self.ref,law)
        b=solve_rough_contact(HalfSpaceLine(L,129,4e6,4*L),self.ref,law)
        np.testing.assert_allclose(a.pressure_pa,b.pressure_pa,rtol=1e-8,atol=1e-6)
        np.testing.assert_allclose(a.mean_gap_m,b.mean_gap_m,rtol=1e-9,atol=1e-15)
        scaled=PerssonModel(self.model.psd,2.4e6,.4)
        t=build_contact_table(scaled,self.table.pressures_pa*.6)
        c=solve_rough_contact(HalfSpaceLine(L,129,2.4e6),HertzLine(.005,60,2.4e6),RoughContactLaw(t))
        np.testing.assert_allclose(a.pressure_pa*.6,c.pressure_pa,rtol=1e-8,atol=1e-6)
        np.testing.assert_allclose(a.mean_gap_m,c.mean_gap_m,rtol=1e-8,atol=1e-15)

    def test_reducing_roughness_approaches_smooth_solution(self):
        el=HalfSpaceLine(3*self.ref.half_width_m,129,4e6);smooth=solve_smooth_contact(el,self.ref)
        errors=[]
        for amplitude in [1.,.5,.25,.1]:
            model=PerssonModel(read_psd(ROOT/'examples/data/roughness_psd.dat',scale_c=amplitude**2),4e6,.4)
            s=solve_rough_contact(el,self.ref,RoughContactLaw(table(model,self.ref,1201)))
            errors.append(np.sum(abs(s.pressure_pa-smooth.pressure_pa))*el.dx_m/100)
        self.assertTrue(np.all(np.diff(errors)<0));self.assertLess(errors[-1],.006)

    def test_small_table_tail_is_bounded_and_high_load_cannot_extrapolate(self):
        law=RoughContactLaw(self.table)
        p,dp=law.pressure_and_derivative([law.h_max_m*2])
        self.assertEqual(p[0],0);self.assertEqual(dp[0],0)
        with self.assertRaises(ValueError):law.pressure_and_derivative([law.h_min_m*.5])
        with self.assertRaisesRegex(ValueError,'did not converge'):
            solve_rough_contact(HalfSpaceLine(3*self.ref.half_width_m,65,4e6),self.ref,law,max_iterations=1)


class MacroWorkflowTests(unittest.TestCase):
    def test_rough_workflow_outputs_and_no_fluid_coupling(self):
        from parameters import PARAMETERS
        p=deepcopy(PARAMETERS);p['psd']['path']=str(ROOT/'examples/data/roughness_psd.dat')
        p['macro_contact']['solver']['cells']=65;p['macro_contact']['contact_table']['points']=601
        with tempfile.TemporaryDirectory() as directory:
            p['output_directory']='macro'
            with self.assertRaisesRegex(NotImplementedError,'完整流体耦合'):
                run_project(p,base_directory=directory,run_macro_contact=True,run_leakage=True)
            self.assertFalse((Path(directory)/'macro').exists())
            r=run_project(p,base_directory=directory,run_macro_contact=True,save_plots=False)
            out=Path(directory)/'macro';self.assertTrue((out/'macro_geometry.csv').is_file())
            self.assertTrue((out/'macro_contact_table.csv').is_file())
            self.assertFalse((out/'leakage_cases.csv').exists())
            self.assertEqual(r['workflow_status'],'complete')
            saved=json.loads((out/'macro_diagnostics.json').read_text())
            self.assertEqual(saved['psd_sha256'],r['psd_sha256'])
            data=np.genfromtxt(out/'macro_geometry.csv',delimiter=',',names=True)
            np.testing.assert_array_equal(data['fluid_pressure_pa'],0.)
            np.testing.assert_allclose(data['display_deformed_upper_surface_m']-data['display_lower_surface_m'],data['mean_gap_m'],atol=1e-16)
            with self.assertRaises(FileExistsError):run_project(p,base_directory=directory,run_macro_contact=True,save_plots=False)

    def test_smooth_standalone_does_not_need_psd_and_invalid_presets_fail(self):
        from parameters import PARAMETERS
        p=deepcopy(PARAMETERS);p.pop('psd');p['macro_contact']['contact_model']='smooth'
        p['macro_contact']['solver']['cells']=65;p['output_directory']='smooth'
        with tempfile.TemporaryDirectory() as directory:
            r=run_project(p,base_directory=directory,run_macro_contact=True,save_plots=False)
            self.assertEqual(r['contact_model'],'smooth');self.assertIsNone(r['psd_sha256'])
            p['macro_contact']['elastic_model']='finite_layer';p['output_directory']='bad'
            with self.assertRaises(NotImplementedError):run_project(p,base_directory=directory,run_macro_contact=True,save_plots=False)
            self.assertFalse((Path(directory)/'bad').exists())

    def test_bad_nested_controls_and_table_bounds_are_explicit(self):
        from parameters import PARAMETERS
        with tempfile.TemporaryDirectory() as directory:
            for block,key,value in [('solver','relative_tolerance',None),('solver','cells',True),
                    ('solver','unexpected',1),('contact_table','minimum_pressure_pa',1.),
                    ('contact_table','maximum_pressure_pa',1e3)]:
                p=deepcopy(PARAMETERS);p['psd']['path']=str(ROOT/'examples/data/roughness_psd.dat')
                p['output_directory']='bad';p['macro_contact'][block][key]=value
                with self.subTest(block=block,key=key),self.assertRaises(ValueError):
                    run_project(p,base_directory=directory,run_macro_contact=True,save_plots=False)
                self.assertFalse((Path(directory)/'bad').exists())


if __name__=='__main__':unittest.main()
