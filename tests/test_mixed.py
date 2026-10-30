"""Independent smooth limits, conservative cavitation and coupled invariants."""
from copy import deepcopy
import json
from pathlib import Path
import tempfile
import unittest

import numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq

from persson_contact.fluid_friction import shear_friction_factors,wall_tractions,scale_inverse_gap
from persson_contact.interface import NewtonianNoSlip
from persson_contact.mixed import film_flux,solve_film
from persson_contact.mixed import solve_mixed
from persson_contact.macroscopic import HalfSpaceLine,HertzLine
from persson_contact.core import PerssonModel
from persson_contact.psd import read_psd
from persson_contact.table import build_contact_table
from persson_contact.dry_contact import RoughContactLaw,solve_rough_contact
from persson_contact.workflow import run_project


class ShearFrictionTests(unittest.TestCase):
    def test_shear_cutoff_is_explicit_and_does_not_reclassify_contact(self):
        class Samples:
            pdf_scale_method='conservative'
            def _conservative_pdf_scales(self,p):
                return np.array([.5e-9,2e-9]),np.array([.2,.3]),np.ones(2),.01,.02,.03
        result=scale_inverse_gap(Samples(),1.,1e-9)
        self.assertAlmostEqual(result['inverse_gap_m_inv'],.3/2e-9)
        self.assertAlmostEqual(result['uncut_inverse_gap_m_inv'],.2/.5e-9+.3/2e-9)
        self.assertAlmostEqual(result['retained_probability'],.5)
        self.assertAlmostEqual(result['excluded_shear_probability'],.2)
        self.assertAlmostEqual(result['unrepresented_noncontact_probability'],.06)

    def test_smooth_wall_exact_local_response_and_force_balance(self):
        h=np.array([1e-6,2e-6,5e-6]);g=np.array([1e8,-2e8,0.]);eta=.1
        f=shear_friction_factors(h,1/h,0.,0.)
        reference=NewtonianNoSlip(eta).local_response(h,g,.02,-.01)
        for rough in ('lower','upper'):
            tr=wall_tractions(h,f,eta,g,.02,-.01,rough_surface=rough)
            np.testing.assert_allclose(tr['lower_total_pa'],reference.lower_wall_shear_stress_pa,rtol=2e-15)
            np.testing.assert_allclose(tr['upper_total_pa'],reference.upper_wall_shear_stress_pa,rtol=2e-15)

    def test_archive_interpolation_kernel(self):
        path=Path(__file__).parent/'data/fortran_shear_friction.json'
        reference=json.loads(path.read_text())
        for row in reference['cases']:
            f=shear_friction_factors(row['h'],row['inverse'],row['rms'],row['hc'])
            np.testing.assert_allclose([f.phi_f,f.phi_fs,f.phi_fp],row['factors'],rtol=3e-14,atol=1e-15)

    def test_rough_force_balance_reversal_and_saturation(self):
        h=np.array([.2,.5,1.,10.]);f=shear_friction_factors(h,2/h,1.,.5)
        for rough in ('lower','upper'):
            tr=wall_tractions(h,f,1.,3.,.1,.4,rough_surface=rough)
            np.testing.assert_allclose(tr['lower_total_pa']+tr['upper_total_pa'],-3*h)
            rev=wall_tractions(h,f,1.,-3.,-.1,-.4,rough_surface=rough)
            np.testing.assert_allclose(rev['lower_total_pa'],-tr['lower_total_pa'])
            partial=wall_tractions(h,f,1.,0.,.1,.4,rough_surface=rough,saturation=.3)
            np.testing.assert_allclose(partial['lower_sliding_pa'],.3*tr['lower_sliding_pa'])
        np.testing.assert_allclose(f.phi_fp[:2],0.)
        np.testing.assert_allclose(f.phi_fs[:2],6.)

    def test_invalid_factor_inputs(self):
        for args in ((0,1,1,1),(1,-1,1,1),(1,1,-1,1),(1,1,1,-1)):
            with self.assertRaises(ValueError):shear_friction_factors(*args)


class ConservativeCavitationTests(unittest.TestCase):
    @staticmethod
    def smooth(eta,u0,u1):
        def transport(h):
            return h**3/(12*eta),h**2/(4*eta),h*(u0+u1)/2,np.full_like(h,(u0+u1)/2)
        return transport

    def test_parallel_poiseuille_couette_exact_and_stationary_zero(self):
        n=51;dx=.001/n;h=np.full(n,2e-6);x=(np.arange(n)+.5)*dx
        for u0,u1,pin,pout in ((.01,0,2e4,1e4),(0,0,0,0),(-.01,0,1e4,2e4)):
            sol=solve_film(h,dx,self.smooth(.1,u0,u1),pressure_scale_pa=1e5,
                inlet_pressure_pa=pin,outlet_pressure_pa=pout)
            np.testing.assert_allclose(sol.pressure_pa,pin+(pout-pin)*x/.001,rtol=1e-12,atol=1e-8)
            expected=h[0]*(u0+u1)/2-h[0]**3/(12*.1)*(pout-pin)/.001
            np.testing.assert_allclose(sol.face_flux_m2_s,expected,rtol=1e-12,atol=1e-23)
            np.testing.assert_array_equal(sol.saturation,np.ones(n))

    def test_cylinder_rigid_jfo_against_independent_integrals(self):
        L=.001;R=.005;h0=1e-6;eta=.1;speed=.01
        gap=lambda x:h0+x*x/(2*R)
        rupture=brentq(lambda r:quad(lambda x:(gap(x)-gap(r))/gap(x)**3,-L,r,epsabs=.1)[0],0,L,xtol=1e-15)
        q_exact=speed/2*gap(rupture);errors=[]
        for n in (65,129,257):
            dx=2*L/n;x=-L+(np.arange(n)+.5)*dx;h=gap(x)
            v=solve_film(h,dx,self.smooth(eta,speed,0),pressure_scale_pa=1e5)
            errors.append(abs(np.mean(v.face_flux_m2_s)/q_exact-1))
            self.assertGreaterEqual(v.pressure_pa.min(),0)
            self.assertTrue(np.all((v.saturation>=0)&(v.saturation<=1)))
            self.assertEqual(float(np.max(v.pressure_pa*(1-v.saturation))),0.)
            self.assertLess(np.ptp(v.face_flux_m2_s)/q_exact,1e-10)
            exact_p=np.array([12*eta*speed/2*quad(lambda t:(gap(t)-gap(rupture))/gap(t)**3,-L,float(z),epsabs=.1)[0] if z<rupture else 0. for z in x])
            if n==257:
                self.assertLess(np.sum(abs(v.pressure_pa-exact_p))/np.sum(exact_p),.01)
                self.assertLess(errors[-1],.003)
                rev=solve_film(h,dx,self.smooth(eta,-speed,0),pressure_scale_pa=1e5)
                np.testing.assert_allclose(rev.pressure_pa,v.pressure_pa[::-1],rtol=1e-10,atol=1e-7)
                np.testing.assert_allclose(rev.saturation,v.saturation[::-1],rtol=1e-10)
                np.testing.assert_allclose(rev.face_flux_m2_s,-v.face_flux_m2_s[::-1],rtol=1e-10)
        self.assertTrue(errors[2]<errors[1]<errors[0])

    def test_flux_jacobian_on_full_and_cavitated_branches(self):
        h=np.array([1.1,1.2,1.5,1.8,2.]);psi=np.array([.3,.2,-.1,-.3,-.5])
        t=self.smooth(.1,.2,0);q,jh,jp,_,_=film_flux(h,psi,.1,t,2.,0.,0.)
        eps=1e-6
        for i in range(len(h)):
            dh=np.zeros_like(h);dh[i]=eps
            num=(film_flux(h+dh,psi,.1,t,2.,0.,0.)[0]-film_flux(h-dh,psi,.1,t,2.,0.,0.)[0])/(2*eps)
            np.testing.assert_allclose(jh[:,i],num,rtol=3e-9,atol=3e-9)
            num=(film_flux(h,psi+dh,.1,t,2.,0.,0.)[0]-film_flux(h,psi-dh,.1,t,2.,0.,0.)[0])/(2*eps)
            np.testing.assert_allclose(jp[:,i],num,rtol=3e-9,atol=3e-9)

    def test_blocked_mobility_is_rejected_not_floored(self):
        def blocked(h):return np.zeros_like(h),np.zeros_like(h),h,np.ones_like(h)
        with self.assertRaisesRegex(ValueError,'blocked/trapped'):
            solve_film(np.ones(10),.1,blocked,pressure_scale_pa=1.)


class MixedConfigurationTests(unittest.TestCase):
    def test_missing_strength_and_unsupported_physics_fail_before_output(self):
        from parameters import PARAMETERS
        root=Path(__file__).resolve().parents[1]
        for mode in ('missing_strength','slip','incompatible_task','bad_boundary','bad_saturation'):
            p=deepcopy(PARAMETERS);p['output_directory']='mixed'
            p['psd']['path']=str(root/'examples/data/roughness_psd.dat')
            if mode!='missing_strength':p['contact_friction']['model_parameters']['constant_shear_strength']['shear_strength_pa']=1e5
            if mode=='slip':p['interface']['wall_slip']='navier'
            if mode=='bad_boundary':p['mixed_friction']['inlet_pressure_pa']=-1
            if mode=='bad_saturation':p['mixed_friction']['cavitation_friction']='delete_shear'
            with self.subTest(mode=mode),tempfile.TemporaryDirectory() as d:
                with self.assertRaises((ValueError,NotImplementedError)):
                    run_project(p,base_directory=d,run_mixed_friction=True,
                        run_leakage=mode=='incompatible_task',save_plots=False)
                self.assertFalse((Path(d)/'mixed').exists())


class CoupledNormalTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        root=Path(__file__).resolve().parents[1]
        cls.model=PerssonModel(read_psd(root/'examples/data/roughness_psd.dat'),4e6,.4,n_q=128,
            pdf_scale_method='conservative')
        pressures=np.unique(np.r_[np.geomspace(1e-100,100,101),np.geomspace(100,2e6,500)])
        cls.law=RoughContactLaw(build_contact_table(cls.model,pressures))
        cls.hertz=HertzLine(.005,100.,4e6)
        cls.elastic=HalfSpaceLine(.0012,65,4e6)
        cls.dry=solve_rough_contact(cls.elastic,cls.hertz,cls.law)

    def test_zero_speed_recovers_dry_contact(self):
        transport=ConservativeCavitationTests.smooth(.1,0.,0.)
        v=solve_mixed(self.elastic,self.hertz,self.law,transport,self.dry)
        np.testing.assert_allclose(v.mean_gap_m,self.dry.mean_gap_m,rtol=1e-12)
        np.testing.assert_allclose(v.contact_pressure_pa,self.dry.pressure_pa,rtol=1e-12)
        np.testing.assert_array_equal(v.fluid_pressure_pa,0.)

    def test_fluid_load_geometry_gap_law_and_reversal(self):
        el=self.elastic;hz=self.hertz
        results=[]
        for speed in (.001,-.001):
            v=solve_mixed(el,hz,self.law,ConservativeCavitationTests.smooth(.1,speed,0.),self.dry)
            total=v.contact_pressure_pa+v.fluid_pressure_pa
            np.testing.assert_allclose(el.line_load(total),hz.line_load_n_m,rtol=1e-8)
            self.assertGreater(el.line_load(v.fluid_pressure_pa),0.)
            self.assertLess(el.line_load(v.contact_pressure_pa),hz.line_load_n_m)
            self.assertLess(v.geometry_relative_residual,1e-8)
            self.assertLess(v.mass_relative_residual,1e-6)
            np.testing.assert_allclose(v.mean_gap_m,hz.profile(el.x_m)+v.rigid_offset_m+el.displacement(total),atol=1e-12)
            active=v.contact_pressure_pa>1e-80
            self.assertLess(np.max(abs(self.model.mean_gap(v.contact_pressure_pa[active])-v.mean_gap_m[active]))/self.model.psd.hrms,1e-4)
            results.append(v)
        np.testing.assert_allclose(results[0].fluid_pressure_pa,results[1].fluid_pressure_pa[::-1],rtol=1e-9,atol=1e-6)
        np.testing.assert_allclose(results[0].mean_gap_m,results[1].mean_gap_m[::-1],rtol=1e-9)

    def test_pressure_reference_length_invariance(self):
        el=HalfSpaceLine(.0012,65,4e6,reference_length_m=.0048)
        dry=solve_rough_contact(el,self.hertz,self.law)
        tr=ConservativeCavitationTests.smooth(.1,.001,0.)
        a=solve_mixed(self.elastic,self.hertz,self.law,tr,self.dry)
        b=solve_mixed(el,self.hertz,self.law,tr,dry)
        np.testing.assert_allclose(a.mean_gap_m,b.mean_gap_m,rtol=1e-8)
        np.testing.assert_allclose(a.fluid_pressure_pa,b.fluid_pressure_pa,rtol=1e-8,atol=1e-6)


if __name__=='__main__':unittest.main()
