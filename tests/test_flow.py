"""Analytic fluid/EMT limits, missing-mass bounds and independent flow grids."""
from copy import deepcopy
import json
from pathlib import Path
import tempfile
import unittest

import numpy as np

from persson_contact.applications import LeakageParameters,run_leakage
from persson_contact.core import PerssonModel
from persson_contact.flow import (FlowOptions,bruggeman_conductivity,effective_dimension,
                                  homogenize_bins,pressure_flow)
from persson_contact.interface import NewtonianNoSlip,make_interface
from persson_contact.psd import read_psd
from persson_contact.workflow import run_project


class InterfaceTests(unittest.TestCase):
    def test_parallel_plate_poiseuille_couette_and_contact(self):
        fluid=NewtonianNoSlip(.1);gap=2e-6
        self.assertAlmostEqual(fluid.pressure_mobility(gap),gap**3/1.2)
        self.assertEqual(fluid.pressure_mobility(0),0)
        self.assertEqual(fluid.flux(0,1e5,1,0),0)
        self.assertAlmostEqual(fluid.flux(gap,3e5,.4,0),-gap**3*3e5/1.2+gap*.2)
        np.testing.assert_allclose(fluid.flux([gap,2*gap],-3e5),
                                  -fluid.flux([gap,2*gap],3e5),rtol=1e-14)
        for value in (None,0,-1,True,float('nan')):
            with self.assertRaises(ValueError):NewtonianNoSlip(value)
        with self.assertRaises(ValueError):fluid.pressure_mobility(-1)
        for cfg in ({'viscosity_pa_s':1,'wall_slip':'navier'},
                    {'viscosity_pa_s':1,'rheology':'carreau'}):
            with self.assertRaises(NotImplementedError):make_interface(cfg)


class EffectiveMediumTests(unittest.TestCase):
    def test_uniform_open_gap_and_closed_phase_analytical_solution(self):
        for area in (0,.1,.3,.499,.5,.8):
            k=3e-18
            value,residual,_=bruggeman_conductivity([k],[1-area],area,2.)
            expected=k*max(0,1-2*area)
            np.testing.assert_allclose(value,expected,rtol=1e-11,atol=0)
            self.assertLess(abs(residual),1e-11)
        # Equal populations of two positive phases in 2D: geometric mean.
        value,_,_=bruggeman_conductivity([1e-20,1e-16],[.5,.5],0,2.)
        np.testing.assert_allclose(value,1e-18,rtol=1e-11)

    def test_percolation_dimension_and_invalid_mass(self):
        self.assertEqual(effective_dimension(0,.42),2)
        self.assertAlmostEqual(effective_dimension(.42,.42),1/.58)
        self.assertEqual(bruggeman_conductivity([1.],[.58],.42,1/.58)[0],0)
        with self.assertRaises(ValueError):bruggeman_conductivity([1.],[.9],.2,2)
        with self.assertRaises(ValueError):bruggeman_conductivity([-1.],[1.],0,2)
        with self.assertRaises(ValueError):FlowOptions(gap_points=50)
        self.assertTrue(np.isinf(bruggeman_conductivity([1.],[.3],.1,2.,.6)[0]))

    def test_cdf_missing_mass_bounds_do_not_renormalize(self):
        edges=np.array([0,1e-6,2e-6,4e-6]);cdf=np.array([0,.1,.35,.6])
        result=homogenize_bins(edges,cdf,.2,1e-6,NewtonianNoSlip(.1),.5)
        self.assertAlmostEqual(result['missing'],.2)
        self.assertAlmostEqual(result['zero'],.4)
        self.assertGreater(result['upper'],result['geometric'])
        # Actual finite unknown phase must lie between zero and infinity limits.
        k=np.array([.5e-6,np.sqrt(2)*1e-6,np.sqrt(8)*1e-6])**3
        middle,_,_=bruggeman_conductivity(np.r_[k,(3e-6)**3],np.r_[np.diff(cdf),.2],.2,2)
        self.assertLessEqual(result['geometric'],middle)
        self.assertGreaterEqual(result['upper'],middle)
        bad=cdf.copy();bad[-1]=.9
        with self.assertRaises(ValueError):homogenize_bins(edges,bad,.2,1e-6,NewtonianNoSlip(1))


class RealSurfaceFlowTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.root=Path(__file__).resolve().parents[1]
        cls.psd=read_psd(cls.root/'examples/data/roughness_psd.dat')
        cls.model=PerssonModel(cls.psd,4e6,.4)

    def test_flow_viscosity_scaling_convergence_and_blocking(self):
        first=pressure_flow(self.model,1e5,NewtonianNoSlip(.1))
        second=pressure_flow(self.model,1e5,NewtonianNoSlip(.2))
        self.assertEqual(first.pressure_flow_factor,second.pressure_flow_factor)
        self.assertAlmostEqual(first.hydraulic_conductivity_m3_pa_s,
                               2*second.hydraulic_conductivity_m3_pa_s)
        self.assertLess(first.quadrature_relative_difference,.002)
        self.assertGreater(first.upper_pressure_flow_factor,first.pressure_flow_factor)
        from unittest.mock import patch
        with patch.object(self.model,'mean_gap',return_value=0.):
            with self.assertRaisesRegex(ValueError,'underflow'):
                pressure_flow(self.model,1e5,NewtonianNoSlip(.1))
        blocked=pressure_flow(self.model,2e5,NewtonianNoSlip(.1))
        self.assertEqual(blocked.status,'blocked')
        self.assertEqual(blocked.pressure_flow_factor,0)
        self.assertEqual(blocked.upper_pressure_flow_factor,0)
        uncertain=pressure_flow(self.model,1.95e5,NewtonianNoSlip(.1))
        self.assertEqual(uncertain.status,'unresolved_connectivity')
        self.assertEqual(uncertain.geometric_conductivity_m3,0)
        self.assertGreater(uncertain.upper_geometric_conductivity_m3,0)
        cfg=LeakageParameters(1e-3,1e-2,1e4,0)
        q=run_leakage(first,cfg,NewtonianNoSlip(.1))
        expected=10*1e4*first.geometric_conductivity_m3/1.2
        self.assertAlmostEqual(q['flow_rate_m3_s']/expected,1.)
        reverse=run_leakage(first,LeakageParameters(1e-3,1e-2,0,1e4),NewtonianNoSlip(.1))
        self.assertEqual(reverse['flow_rate_m3_s'],-q['flow_rate_m3_s'])
        self.assertLessEqual(reverse['min_flow_rate_m3_s'],reverse['max_flow_rate_m3_s'])
        self.assertEqual(run_leakage(first,LeakageParameters(1,1,1,1),NewtonianNoSlip(.1))['flow_rate_m3_s'],0)
        from dataclasses import replace
        unbounded=replace(first,upper_bound_unbounded=True,upper_geometric_conductivity_m3=None)
        zero=run_leakage(unbounded,LeakageParameters(1,1,1,1),NewtonianNoSlip(.1))
        self.assertEqual(zero['max_flow_rate_m3_s'],0)
        self.assertFalse(zero['bound_unbounded'])

    def test_original_fortran_closure_references(self):
        fixture=json.loads((self.root/'tests/data/fortran_pressure_flow.json').read_text())
        self.assertEqual(len(fixture['cases']),9)
        for case in fixture['cases']:
            with self.subTest(case=case['name']):
                model=PerssonModel(self.psd,case['effective_modulus_pa'],.4,
                                   pdf_log_step=case['delzeta'])
                result=pressure_flow(model,case['pressure_pa'],NewtonianNoSlip(.1))
                np.testing.assert_allclose(result.geometric_conductivity_m3,
                    case['reference_geometric_conductivity_m3'],rtol=1e-4,atol=0)

    def test_nonconverged_flow_marks_partial_output_as_failed(self):
        from parameters import PARAMETERS
        with tempfile.TemporaryDirectory() as directory:
            params=deepcopy(PARAMETERS)
            params['psd']['path']=str(self.root/'examples/data/roughness_psd.dat')
            params.pop('contact_pressure');params['pressures_pa']=[1e4,5e4,1e5]
            params['output_directory']='failed'
            params['pressure_flow']['quadrature_relative_tolerance']=1e-12
            params['pressure_flow']['max_gap_refinements']=0
            with self.assertRaisesRegex(ValueError,'not converged'):
                run_project(params,base_directory=directory,compute_pressure_flow=True,
                            compute_gap_pdf=False,save_plots=False)
            diag=json.loads((Path(directory)/'failed/diagnostics.json').read_text())
            self.assertEqual(diag['workflow_status'],'flow_failed')
            self.assertIn('flow_error',diag)
            self.assertFalse((Path(directory)/'failed/pressure_flow_table.csv').exists())

    def test_leakage_dependencies_and_pdf_export_are_independent(self):
        from parameters import PARAMETERS
        with tempfile.TemporaryDirectory() as directory:
            params=deepcopy(PARAMETERS)
            params['psd']['path']=str(self.root/'examples/data/roughness_psd.dat')
            params.pop('contact_pressure');params['pressures_pa']=[1e4,1e5,2e5]
            params['output_directory']='flow'
            summary=run_project(params,base_directory=directory,run_contact=False,run_leakage=True,
                                compute_gap_pdf=False,compute_pressure_flow=False,save_plots=False)
            out=Path(directory)/'flow'
            self.assertEqual(summary['workflow_status'],'complete')
            self.assertFalse(list(out.glob('gap_pdf_*.csv')))
            self.assertTrue((out/'pressure_flow_table.csv').is_file())
            self.assertTrue((out/'leakage_table.csv').is_file())
            diag=json.loads((out/'flow_diagnostics.json').read_text())
            self.assertEqual(len(diag['results']),3)
            self.assertEqual(diag['results'][-1]['status'],'blocked')
            params['output_directory']='bad';params['leakage']['seal_length_m']=None
            with self.assertRaises(ValueError):run_project(params,base_directory=directory,run_leakage=True)
            self.assertFalse((Path(directory)/'bad').exists())


if __name__=='__main__':unittest.main()
