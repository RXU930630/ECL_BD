# Reproducible material and friction verification, without experimental results.
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'src'))
import tempfile,unittest,json
from pathlib import Path
import numpy as np
from scipy.special import erf
from scipy.integrate import solve_ivp
from numpy.polynomial.legendre import leggauss
from persson_contact.viscoelastic_material import MasterCurve
from persson_contact.viscoelastic_friction import PerssonViscoelasticFriction
from persson_contact.psd import PowerSpectrum

class AnalyticSLS:
    nu=.5
    def storage_loss(self,omega):
        x=np.asarray(omega)*.01
        return 1e6+9e6*x*x/(1+x*x),9e6*x/(1+x*x)

class ViscoelasticTests(unittest.TestCase):
    def psd(self):return PowerSpectrum(np.array([3.,6.]),np.array([-18.,-30.]),1e3,1e3,1e6)
    def files(self,root):
        f=np.geomspace(1e-10,1e12,2201);ep,loss=AnalyticSLS().storage_loss(2*np.pi*f)
        path=root/'curve.dat';np.savetxt(path,np.c_[f,ep/1e6,loss/1e6],header='SLS',comments='')
        shift=root/'shift.dat';shift.write_text('SLS\n10 10 1\n20 1 1\n30 .1 1\n');return path,shift
    def test_units_shifts_and_sign_convention(self):
        with tempfile.TemporaryDirectory() as d:
            p,s=self.files(Path(d));c=MasterCurve(p,.5,reference_temperature_c=20,temperature_c=10,shift_factors_file=s)
            w=np.array([1.,100.,1e4]);ep,loss=c.storage_loss(w);a,b=AnalyticSLS().storage_loss(w*10)
            np.testing.assert_allclose(ep,a,rtol=2e-5);np.testing.assert_allclose(loss,b,rtol=2e-5)
            np.testing.assert_allclose(c.complex_modulus(-w),np.conj(c.complex_modulus(w)))
            self.assertTrue(np.all(c.complex_modulus(w).imag>0))
    def test_explicit_types_and_frequency_bounds(self):
        with tempfile.TemporaryDirectory() as d:
            p,s=self.files(Path(d))
            y=MasterCurve(p,.5,reference_temperature_c=20,temperature_c=20)
            g=MasterCurve(p,.5,reference_temperature_c=20,temperature_c=20,modulus_type='shear')
            np.testing.assert_allclose(g.storage_loss(10)[0],3*y.storage_loss(10)[0])
            with self.assertRaisesRegex(ValueError,'coverage'):y.storage_loss(1e20)
            with self.assertRaises(ValueError):MasterCurve(p,.5,reference_temperature_c=10,temperature_c=10,shift_factors_file=s)
            with self.assertRaises(ValueError):MasterCurve(p,.5,reference_temperature_c=20,temperature_c=10)
            with self.assertRaises(ValueError):y.storage_loss(0)
    def test_bad_data_not_silently_reordered_or_clipped(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'bad';p.write_text('bad\n2 1 1\n1 1 1\n')
            with self.assertRaises(ValueError):MasterCurve(p,.5,reference_temperature_c=20,temperature_c=20)
            p.write_text('bad\n1 1 -1\n2 1 1\n')
            with self.assertRaises(ValueError):MasterCurve(p,.5,reference_temperature_c=20,temperature_c=20)
    def test_elastic_limit_recovers_erf_area_and_zero_hysteresis(self):
        class Elastic:
            nu=.5
            def storage_loss(self,w):return np.full_like(w,2e6),np.zeros_like(w)
        psd=self.psd();m=PerssonViscoelasticFriction(psd,Elastic(),q_intervals=64,angular_points=16)
        for p in [1e3,1e4,1e5]:
            r=m.evaluate(p,.001,1e6);expected=erf(p/(2e6/(1-.5**2)*np.sqrt(np.pi*psd.moment(3))))
            self.assertAlmostEqual(r['relative_contact_area'],expected,places=13)
            self.assertEqual(r['mu_viscoelastic'],0);self.assertAlmostEqual(r['mu_total'],1e6*expected/p)
    def test_full_contact_quadrant_factor_and_pressure_normalization(self):
        class ConstantKernel:
            nu=.5
            def storage_loss(self,w):return np.full_like(w,2e6),np.full_like(w,3e5)
        psd=self.psd();r=PerssonViscoelasticFriction(psd,ConstantKernel(),q_intervals=64,angular_points=32).evaluate(1e12,.01)
        expected=2*3e5*psd.moment(3)/(1e12*(1-.5**2))
        np.testing.assert_allclose(r['mu_viscoelastic'],expected,rtol=1e-12)
    def test_against_independent_continuous_ode_quadrature(self):
        psd=self.psd();material=AnalyticSLS();p=1e4;speed=.001
        x,w=leggauss(96);theta=(x+1)*np.pi/4;weights=w*np.pi/4;cos=np.cos(theta)
        def derivative(t,state):
            q=np.exp(t);ep,loss=material.storage_loss(q*speed*cos)
            square=4*np.sum(weights*(ep*ep+loss*loss))/(1-.5**2)**2
            damping=4*np.sum(weights*cos*loss)/(1-.5**2)
            P=1. if state[0]<=0 else erf(1/(2*np.sqrt(state[0])))
            return [q**4*psd.evaluate(q)*square/(8*p*p),q**4*psd.evaluate(q)*P*(.5+.5*P*P)*damping/(2*p)]
        ref=solve_ivp(derivative,[np.log(psd.q0),np.log(psd.q1)],[0.,0.],method='DOP853',rtol=2e-10,atol=1e-12)
        r=PerssonViscoelasticFriction(psd,material,q_intervals=2048,angular_points=96).evaluate(p,speed)
        np.testing.assert_allclose(r['mu_viscoelastic'],ref.y[1,-1],rtol=3e-5)
        np.testing.assert_allclose(r['relative_contact_area'],erf(1/(2*np.sqrt(ref.y[0,-1]))),rtol=3e-5)
        self.assertAlmostEqual(r['mu_viscoelastic'],PerssonViscoelasticFriction(psd,material,q_intervals=2048,angular_points=96).evaluate(p,-speed)['mu_viscoelastic'])

    def test_independent_workflow_and_coupling_guards(self):
        from persson_contact.workflow import run_project
        with tempfile.TemporaryDirectory() as d:
            root=Path(d);curve,shift=self.files(root)
            (root/'psd.dat').write_text('synthetic\n1e3\n1e3\n1e6\n0\n1\n3 -18\n6 -30\n')
            p={'material':{'model':'viscoelastic','master_curve_file':str(curve),'poisson_ratio':.5},
                'psd':{'path':'psd.dat'},'macro_contact':{'deformable_body':'lower'},
                'mixed_friction':{'rough_surface':'upper'},
                'contact_friction':{'model':'constant_shear_strength','model_parameters':{'constant_shear_strength':{'shear_strength_pa':1e5}}},
                'viscoelastic':{'reference_temperature_c':20.,'temperature_c':20.,'shift_factors_file':str(shift),
                    'pressures_pa':[1e4],'speeds_m_s':[.001],'q_intervals':128,'angular_points':16,'relative_tolerance':.01},
                'output_directory':'run'}
            result=run_project(p,base_directory=root,run_contact=False,run_viscoelastic_friction=True,save_plots=False)
            self.assertEqual(result['workflow_status'],'complete')
            r=result['cases'][0];self.assertGreater(r['mu_viscoelastic'],0)
            self.assertAlmostEqual(r['mu_total'],r['mu_viscoelastic']+r['mu_adhesive']+r['mu_fluid'])
            self.assertTrue((root/'run/viscoelastic_results.xlsx').exists())
            p['output_directory']='forbidden'
            with self.assertRaises(NotImplementedError):run_project(p,base_directory=root,run_contact=False,run_mixed_friction=True,save_plots=False)
            self.assertFalse((root/'forbidden').exists())
            p['mixed_friction']['rough_surface']='lower'
            with self.assertRaises(NotImplementedError):run_project(p,base_directory=root,run_contact=False,run_viscoelastic_friction=True,save_plots=False)
            self.assertFalse((root/'forbidden').exists())

if __name__=='__main__':unittest.main()
