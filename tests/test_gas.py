"""Gas dimensions, independent analytic limits and actual 500/5000 entry paths."""
from copy import deepcopy
from dataclasses import replace
import csv
import json
from pathlib import Path
import tempfile
import unittest

import numpy as np

from persson_contact.ballistic import ballistic_pressure_flow
from persson_contact.core import PerssonModel
from persson_contact.flow import homogenize_bins,pressure_flow
from persson_contact.gas import GasProperties,BallisticGas,equivalent_gas_estimates
from persson_contact.interface import NewtonianNoSlip
from persson_contact.psd import read_psd
from persson_contact.leakage import prepare_leakage,solve_case,run_spatial_leakage
from persson_contact.transport_table import TransportTable
from persson_contact.workflow import run_project


class AnalyticContact:
    def area(self,p):return np.asarray(p)/1e6


class GasLawTests(unittest.TestCase):
    def test_original_fortran_quadratic_closure_references(self):
        root=Path(__file__).resolve().parents[1]
        fixture=json.loads((root/'tests/data/fortran_ballistic_flow.json').read_text())
        self.assertEqual(len(fixture['cases']),2)
        psd=read_psd(root/'examples/data/roughness_psd.dat')
        for case in fixture['cases']:
            model=PerssonModel(psd,case['effective_modulus_pa'],.4,pdf_log_step=case['log_step'])
            value=ballistic_pressure_flow(model,case['pressure_pa'],BallisticGas(GasProperties()))
            np.testing.assert_allclose(value.geometric_conductivity_m2,case['fortran_K2_m2'],rtol=6e-4)
            a,b=case['refinements']
            self.assertLess(abs(a['fortran_K2_m2']/b['fortran_K2_m2']-1),1e-8)
            self.assertLess(abs(b['same_pdf_nonnegative_u2_root_relative_difference']),4e-6)
        model=PerssonModel(psd,2.4e6,.4)
        with self.assertRaises(ValueError):pressure_flow(model,1e4,BallisticGas(GasProperties()))
        with self.assertRaises(ValueError):ballistic_pressure_flow(model,1e4,NewtonianNoSlip(.001))

    def test_ballistic_parallel_gap_and_contact_mixture(self):
        gas=GasProperties();law=BallisticGas(gas);u=2e-6
        expected=gas.mean_molecular_speed_m_s*u*u/(2*gas.thermal_energy_j)
        self.assertAlmostEqual(law.pressure_mobility(u)/expected,1.)
        self.assertEqual(law.pressure_mobility(0),0.)
        edges=np.array([0.,2*u,4*u])
        for a in (0.,.2,.49,.5,.8):
            result=homogenize_bins(edges,np.array([0.,1-a,1-a]),a,u,law,.5)
            np.testing.assert_allclose(result['geometric'],u*u*max(1-2*a,0),rtol=1e-10,atol=0)
        quadratic=homogenize_bins(edges,np.array([0.,.8,.8]),.2,u,law,.5)['geometric']
        cubic=homogenize_bins(edges,np.array([0.,.8,.8]),.2,u,NewtonianNoSlip(.001),.5)['geometric']
        self.assertGreater(abs(quadratic/cubic**(2/3)-1),.1)

    def test_equivalent_estimates_match_historical_500_columns(self):
        # Actual user output at 414370 Pa: square-unit liquid/diffusive/ballistic columns.
        result=equivalent_gas_estimates(1.1429e-16,414370.,.0013814,.0013814,
            .001,414370.,0.,GasProperties())
        np.testing.assert_allclose(result['estimated_diffusive_gas_reference_flow_m3_s'],1.9290e-14,rtol=3e-5)
        np.testing.assert_allclose(result['estimated_ballistic_gas_reference_flow_m3_s'],2.1673e-13,rtol=5e-5)
        np.testing.assert_allclose(result['equivalent_hydraulic_gap_m'],1.4903e-8,rtol=3e-5)
        reverse=equivalent_gas_estimates(-1.1429e-16,-414370.,.0013814,.0013814,
            .001,0.,414370.,GasProperties())
        for key in ('estimated_diffusive_gas_reference_flow_m3_s','estimated_ballistic_gas_reference_flow_m3_s'):
            self.assertEqual(result[key],-reverse[key])
        zero=equivalent_gas_estimates(0.,0.,1.,1.,.001,0.,0.,GasProperties())
        self.assertEqual(zero['estimated_ballistic_gas_reference_flow_m3_s'],0.)
        self.assertIsNone(zero['equivalent_hydraulic_gap_m'])
        blocked=equivalent_gas_estimates(0.,1.,1.,1.,.001,1.,0.,GasProperties())
        self.assertEqual(blocked['equivalent_hydraulic_gap_m'],0.)
        self.assertEqual(blocked['estimated_diffusive_gas_reference_flow_m3_s'],0.)

    def test_gas_pressure_reference_temperature_and_invalid_inputs(self):
        gas=GasProperties();self.assertAlmostEqual(gas.temperature_k,298.15)
        self.assertEqual(float(gas.absolute_pressure(-1e5)),0.)
        with self.assertRaises(ValueError):gas.absolute_pressure(-100001.)
        for kwargs in ({'temperature_c':-273.15},{'temperature_c':True},{'temperature_c':None},
                       {'reference_pressure_pa':0},{'mean_molecular_speed_m_s':np.nan}):
            with self.assertRaises(ValueError):GasProperties(**kwargs)
        q=gas.reference_volume(1e10)
        self.assertEqual(replace(gas,reference_pressure_pa=2e5).reference_volume(1e10),q/2)


class BallisticSpatialTests(unittest.TestCase):
    def test_original_5000_spatial_reference_with_shared_K2(self):
        from persson_contact.mechanics import ReferenceMechanics
        fixture=json.loads((Path(__file__).parent/'data/fortran_ballistic_spatial.json').read_text())
        shared=fixture['shared_table'];geom=fixture['geometry']
        gas=GasProperties(**fixture['gas_parameters'])
        table=TransportTable(shared['pressures_pa'],shared['mean_gaps_m'],shared['K2_m2'],
            shared['K2_m2'],AnalyticContact(),BallisticGas(gas),interpolation='linear',
            percolation_policy='legacy_spc_plateau',plateau_pressure_pa=shared['plateau_pressure_pa'])
        plan=self.plan();plan.width_m=geom['width_m'];plan.gas=gas
        plan.percolation_policy='legacy_spc_plateau';plan.coupling='local_unloading'
        x=np.array(geom['x_m']);plan.mechanics=ReferenceMechanics(x,np.array(geom['contact_reference_pa']),
            np.zeros_like(x),geom['metadata'],'linear')
        plan.solver.update(nodes=401,tolerance=1e-6,max_nodes=80001)
        for case in fixture['cases']:
            result=solve_case(plan,table,case['inlet_pressure_pa'])
            n=result['flux']*plan.width_m;q=gas.reference_volume(n)
            target=case['refinements'][-1]
            np.testing.assert_allclose(q,target['fortran_reference_volume_m3_s'],rtol=.002)
            np.testing.assert_allclose(n,target['fortran_molecular_flow_per_s'],rtol=.002)
            self.assertLess(result['pressure_drop_error'],1e-8)

    def plan(self,pin=20000.,out=0.,mode=5000):
        return prepare_leakage({'mode':'spatial','fortran_mode':mode,
            'compute_equivalent_gas':mode==500,'gas':{},'geometry':{'type':'square_ring',
                'length_m':.002,'reference_contact_pressure_pa':1e5},
            'seal_width_m':.01,'inlet_pressure_pa':pin,'outlet_pressure_pa':out,
            'coupling':'fixed_contact','percolation_policy':'physical_zero'},'.',4e6)

    def table(self,gas=None,blocked=False):
        p=np.array([1.,1e4,1e5,2e5]);h=1e-6/(1+p/1e5)
        k=np.full(4,0. if blocked else 4e-12)
        return TransportTable(p,h,k,k,AnalyticContact(),BallisticGas(gas or GasProperties()))

    def test_conserved_molecular_flux_and_reference_volume(self):
        table=self.table();gas=table.interface.gas
        for a,b in ((20000.,0.),(0.,20000.),(123.,123.)):
            plan=self.plan(a,b);result=solve_case(plan,table,a)
            expected=.01/.002*4e-12*gas.mean_molecular_speed_m_s/(2*gas.thermal_energy_j)*(a-b)
            np.testing.assert_allclose(result['flux']*.01,expected,rtol=1e-10,atol=0)
        with tempfile.TemporaryDirectory() as directory:
            data=run_spatial_leakage(self.plan(),table,directory)
            row=data['cases'][0]
            expected=.01/.002*4e-12*gas.mean_molecular_speed_m_s/(2*gas.reference_pressure_pa)*20000.
            np.testing.assert_allclose(row['gas_reference_volume_flow_m3_s'],expected,rtol=1e-10)
            self.assertEqual(row['flow_quantity'],'gas_volume_at_reference_pressure')
            self.assertIn('molecular_flow_rate_per_s',row)
            self.assertNotIn('estimated_ballistic_gas_reference_flow_m3_s',row)
            with (Path(directory)/'leakage_transport_table.csv').open() as f:
                names=next(csv.reader(f))
            self.assertIn('raw_geometric_conductivity_m2',names)
            self.assertNotIn('raw_geometric_conductivity_m3',names)

    def test_ballistic_viscosity_independence_and_speed_temperature_scaling(self):
        gas=GasProperties();base=solve_case(self.plan(),self.table(gas),20000.)['flux']
        for alternative in (replace(gas,viscosity_pa_s=1e-3),replace(gas,ambient_pressure_pa=2e5)):
            np.testing.assert_allclose(solve_case(self.plan(),self.table(alternative),20000.)['flux'],base,rtol=1e-12)
        fast=replace(gas,mean_molecular_speed_m_s=2*gas.mean_molecular_speed_m_s)
        np.testing.assert_allclose(solve_case(self.plan(),self.table(fast),20000.)['flux'],2*base,rtol=1e-12)
        warm=replace(gas,temperature_c=2*gas.temperature_k-273.15)
        np.testing.assert_allclose(solve_case(self.plan(),self.table(warm),20000.)['flux'],base/2,rtol=1e-12)
        blocked=solve_case(self.plan(),self.table(blocked=True),20000.)
        self.assertEqual(blocked['flux'],0.)
        self.assertEqual(blocked['status'],'blocked_pressure_not_unique')

    def test_mode_validation_and_liquid_estimate_switch(self):
        raw=self.plan().original
        for change in ({'fortran_mode':True},{'fortran_mode':5001},{'transport':'liquid'},
                       {'gas':{'temperature_c':-300}},{'outlet_pressure_pa':-1.1e5},
                       {'conductance_model':'critical_junction'}):
            with self.assertRaises(ValueError):prepare_leakage({**raw,**change},'.',4e6)
        liquid_plan=self.plan(mode=500)
        p=np.array([1.,1e4,1e5,2e5]);h=1e-6/(1+p/1e5)
        table=TransportTable(p,h,np.full(4,1e-18),np.full(4,1e-18),AnalyticContact(),NewtonianNoSlip(.001))
        with tempfile.TemporaryDirectory() as directory:
            row=run_spatial_leakage(liquid_plan,table,directory)['cases'][0]
            self.assertIn('estimated_ballistic_gas_reference_flow_m3_s',row)
            self.assertNotIn('molecular_flow_rate_per_s',row)
            liquid_plan.compute_equivalent_gas=False
            row=run_spatial_leakage(liquid_plan,table,directory)['cases'][0]
            self.assertNotIn('estimated_ballistic_gas_reference_flow_m3_s',row)

    def test_fluid_selector_and_legacy_conflicts(self):
        raw=dict(self.plan().original)
        raw.pop('fortran_mode')
        for fluid,mode in (('liquid',500),('gas',5000)):
            plan=prepare_leakage({**raw,'fluid':fluid},'.',4e6)
            self.assertEqual(plan.fluid,fluid)
            self.assertEqual(plan.fortran_mode,mode)
            compatible=prepare_leakage({**raw,'fluid':fluid,'fortran_mode':mode},'.',4e6)
            self.assertEqual(compatible.fortran_mode,mode)
        for fluid in (None,True,500,'ballistic','Gas',[],{}):
            with self.subTest(fluid=fluid),self.assertRaisesRegex(ValueError,'fluid'):
                prepare_leakage({**raw,'fluid':fluid},'.',4e6)
        for change in ({'fluid':'gas','fortran_mode':500},
                       {'fluid':'liquid','fortran_mode':5000},
                       {'fluid':'gas','transport':'liquid'}):
            with self.subTest(change=change),self.assertRaisesRegex(ValueError,'conflicts'):
                prepare_leakage({**raw,**change},'.',4e6)


class GasWorkflowTests(unittest.TestCase):
    def test_historical_case_reads_original_inputs_with_new_selector(self):
        import hashlib
        from persson_contact.config import problem_from_mapping
        root=Path(__file__).resolve().parents[1]
        PARAMETERS=json.loads((root/'tests/data/historical_seal/parameters.json').read_text())
        problem,options=problem_from_mapping({k:v for k,v in PARAMETERS.items()
            if k in ('description','psd','material','pressures_pa','options')},root)
        self.assertAlmostEqual(problem.effective_modulus_pa,2.07e6/.75)
        self.assertEqual(options.pdf_log_step,.005)
        self.assertEqual(len(problem.pressures_pa),150)
        plan=prepare_leakage(PARAMETERS['leakage'],root,problem.effective_modulus_pa)
        self.assertAlmostEqual(plan.width_m,2*np.pi*.009525)
        self.assertEqual(plan.fluid,'liquid')
        self.assertEqual(plan.inlet_pressures_pa[0],414370.)
        self.assertEqual(plan.mechanics.metadata['sha256'],
            '59c5f8a6a5d6f24e694fab536888afe515804d378c26269fcad9ac979dcb1e94')
        self.assertEqual(hashlib.sha256((root/PARAMETERS['psd']['path']).read_bytes()).hexdigest(),
            'abcc26e1d55e576148c30e123b232763fcc4d881683638121e06c0af9f3580d4')

    def test_parameters_switch_runs_both_modes_and_preserves_liquid_schema(self):
        from parameters import PARAMETERS
        root=Path(__file__).resolve().parents[1]
        with tempfile.TemporaryDirectory() as directory:
            for mode in (500,5000):
                cfg=deepcopy(PARAMETERS);cfg['psd']['path']=str(root/'examples/data/roughness_psd.dat')
                cfg.pop('contact_pressure');cfg['pressures_pa']=[1e4,5e4,1e5]
                cfg['leakage']['geometry']={'type':'square_ring','length_m':.001,'reference_contact_pressure_pa':5e4}
                cfg['leakage']['fluid']='liquid' if mode==500 else 'gas'
                cfg['leakage']['inlet_pressures_pa']=[0.,1e4]
                cfg['leakage']['coupling']='fixed_contact';cfg['output_directory']=str(mode)
                if mode==5000:cfg.pop('interface') # No liquid viscosity is needed for independent gas.
                result=run_project(cfg,base_directory=directory,run_leakage=True,
                    compute_gap_pdf=False,save_plots=False)
                out=Path(directory)/str(mode);self.assertEqual(result['workflow_status'],'complete')
                diag=json.loads((out/'leakage_diagnostics.json').read_text())
                self.assertEqual(diag['fortran_mode'],mode)
                self.assertEqual(diag['fluid'],cfg['leakage']['fluid'])
                self.assertEqual(diag['cases'][0]['fluid'],cfg['leakage']['fluid'])
                self.assertEqual(diag['cases'][0]['flow_rate_m3_s'],0.)
                if mode==5000:
                    self.assertTrue((out/'ballistic_flow_table.csv').exists())
                    self.assertFalse((out/'pressure_flow_table.csv').exists())
                    self.assertIn('molecular_flow_rate_per_s',diag['cases'][1])
                else:
                    self.assertTrue((out/'pressure_flow_table.csv').exists())
                    self.assertIn('estimated_ballistic_gas_reference_flow_m3_s',diag['cases'][1])


if __name__=='__main__':unittest.main()
