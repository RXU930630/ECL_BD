"""Independent formulas and grid checks for the conservative numerical path."""
from copy import deepcopy
import json
from pathlib import Path
import tempfile
import unittest

import numpy as np

from persson_contact.config import Options
from persson_contact.core import PerssonModel
from persson_contact.flow import bruggeman_conductivity,pressure_flow
from persson_contact.interface import NewtonianNoSlip
from persson_contact.psd import PowerSpectrum,read_psd
from persson_contact.workflow import run_project

ROOT=Path(__file__).resolve().parents[1]


class ExactMomentTests(unittest.TestCase):
    def test_vector_moments_against_power_law_primitives(self):
        # C=1e-4 q^-4: includes an exactly logarithmic integral.
        psd=PowerSpectrum([0,1,2],[-4,-8,-12],1,1,100)
        lo=np.array([1.,3.,10.,90.]);hi=np.array([2.,40.,20.,100.])
        for order in [1,3,5]:
            k=order-3
            exact=1e-4*np.log(hi/lo) if k==0 else 1e-4*(hi**k-lo**k)/k
            np.testing.assert_allclose(psd.moments(order,lo,hi),exact,rtol=5e-14)
        np.testing.assert_allclose(psd.moments(3,lo[:,None],hi[None,:].max(axis=1)),
                                   1e-4*np.log(100/lo[:,None]),rtol=5e-14)
        self.assertEqual(psd.moments(1,10.,10.),0.)
        with self.assertRaises(ValueError):psd.moments(1,20.,10.)

    def test_small_terminal_intervals_are_not_prefix_differences(self):
        psd=PowerSpectrum([0,1,2],[-4,-8,-12],1,1,100)
        lo=np.array([99.,99.9999,100*(1-1e-10)])
        exact=1e-4*np.log1p((100-lo)/lo)
        np.testing.assert_allclose(psd.moments(3,lo,100.),exact,rtol=1e-5)
        self.assertTrue(np.all(psd.moments(1,lo,100.)>0))


class StableClosureTests(unittest.TestCase):
    def test_small_percolation_margin_against_exact_binary_medium(self):
        for z in [.1,.49,.5-1e-9,.5-1e-12,.5-1e-14]:
            w=1-z
            for scale in [1e-200,1.,1e200]:
                value,_,_=bruggeman_conductivity([scale],[w],z,2.)
                # Exact solution expressed in the supplied probabilities.
                exact=scale*(w-z)
                np.testing.assert_allclose(value,exact,rtol=2e-11)

    def test_near_infinite_phase_threshold(self):
        for infinite in [.1,.49,.5-1e-12,.5-1e-14]:
            w=1-infinite
            value,_,_=bruggeman_conductivity([1.],[w],0,2.,infinite)
            np.testing.assert_allclose(value,1/(w-infinite),rtol=2e-11)


class ConservativePDFTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):cls.psd=read_psd(ROOT/'examples/data/roughness_psd.dat')

    def model(self,**kwargs):
        return PerssonModel(self.psd,4e6,.4,pdf_scale_method='conservative',**kwargs)

    def test_probability_budget_and_first_moment_are_not_fitted(self):
        model=self.model(pdf_log_step=.02)
        grid=np.r_[0.,np.geomspace(1e-14,40*self.psd.hrms,400)]
        pdf=model.gap_distribution(1e5,grid);d=pdf.diagnostics
        self.assertEqual(d.approximation,'conservative_reflected_gaussian')
        self.assertEqual(d.quadrature_mass_error,0.)
        self.assertAlmostEqual(d.resolved_continuous_mass+d.unresolved_continuous_mass,
                               1-model.area(1e5),places=14)
        self.assertTrue(np.all(pdf.continuous_density_per_m>=0))
        self.assertTrue(np.all(np.diff(pdf.continuous_cdf)>=-1e-15))
        self.assertGreater(d.unresolved_continuous_mass,.001)
        self.assertGreater(abs(d.first_moment_relative_error),.001)
        self.assertLess(pdf.continuous_cdf[-1],1-model.area(1e5))
        # One-pressure cache may not accidentally reuse the preceding load.
        other=model.gap_distribution(5e4,grid)
        self.assertFalse(np.array_equal(other.continuous_cdf,pdf.continuous_cdf))
        np.testing.assert_array_equal(model.gap_distribution(1e5,grid).continuous_cdf,
                                      pdf.continuous_cdf)

    def test_scale_refinement_and_modulus_scaling(self):
        a=pressure_flow(self.model(pdf_log_step=.02),1e5,NewtonianNoSlip(.1))
        b=pressure_flow(self.model(pdf_log_step=.01),1e5,NewtonianNoSlip(.1))
        c=pressure_flow(self.model(pdf_log_step=.005),1e5,NewtonianNoSlip(.1))
        self.assertLess(abs(a.geometric_conductivity_m3/c.geometric_conductivity_m3-1),.001)
        self.assertLess(abs(b.geometric_conductivity_m3/c.geometric_conductivity_m3-1),.001)
        scaled=PerssonModel(self.psd,2.4e6,.4,pdf_scale_method='conservative',pdf_log_step=.02)
        r=pressure_flow(scaled,6e4,NewtonianNoSlip(.2))
        np.testing.assert_allclose(r.geometric_conductivity_m3,a.geometric_conductivity_m3,rtol=1e-11)
        np.testing.assert_allclose(r.hydraulic_conductivity_m3_pa_s,
                                   a.hydraulic_conductivity_m3_pa_s/2,rtol=1e-11)

    def test_empty_retained_range_and_blocked_state(self):
        model=self.model(pdf_log_step=.02)
        p=model.pressure_for_relative_contact_area(.9995)
        pdf=model.gap_distribution(p,[0.,self.psd.hrms])
        self.assertEqual(pdf.diagnostics.resolved_continuous_mass,0.)
        self.assertAlmostEqual(pdf.diagnostics.unresolved_continuous_mass,1-model.area(p))
        self.assertEqual(pressure_flow(model,2e5,NewtonianNoSlip(.1)).status,'blocked')

    def test_historical_api_default_and_invalid_choice(self):
        grid=np.array([0.,1e-8,1e-6])
        a=PerssonModel(self.psd,4e6,.4).gap_distribution(1e5,grid)
        b=PerssonModel(self.psd,4e6,.4,pdf_scale_method='fortran').gap_distribution(1e5,grid)
        np.testing.assert_array_equal(a.continuous_cdf,b.continuous_cdf)
        with self.assertRaises(ValueError):self.model(pdf_log_step=1e-6).gap_distribution(1e5,grid)
        with self.assertRaises(ValueError):PerssonModel(self.psd,4e6,.4,pdf_scale_method='unknown')
        with self.assertRaises(ValueError):Options(.4,64,1.,pdf_scale_method='unknown')

    def test_parameter_and_workflow_provenance(self):
        from parameters import PARAMETERS
        params=deepcopy(PARAMETERS)
        params['psd']['path']=str(ROOT/'examples/data/roughness_psd.dat')
        params.pop('contact_pressure');params['pressures_pa']=[1e4,5e4,1e5]
        params['options'].update(pdf_scale_method='conservative',pdf_log_step=.04)
        params['output_directory']='check'
        with tempfile.TemporaryDirectory() as directory:
            run_project(params,base_directory=directory,compute_pressure_flow=True,
                        compute_gap_pdf=False,save_plots=False)
            out=Path(directory)/'check'
            flow=json.loads((out/'flow_diagnostics.json').read_text())
            self.assertEqual(flow['pdf_scale_method'],'conservative')
            self.assertEqual(flow['pdf_log_step'],.04)
            self.assertIn('gap-bin',flow['convergence_scope'])
            contact=json.loads((out/'diagnostics.json').read_text())
            self.assertEqual(contact['options']['pdf_scale_method'],'conservative')


if __name__=='__main__':unittest.main()
