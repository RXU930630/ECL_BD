"""Spatial leakage: physical limits, coordinate preservation and module contracts."""
from copy import deepcopy
import json
from pathlib import Path
import tempfile
import unittest

import numpy as np

from persson_contact.interface import NewtonianNoSlip
from persson_contact.leakage import prepare_leakage,solve_case,run_spatial_leakage
from persson_contact.mechanics import make_mechanics
from persson_contact.profiles import read_pressure_profile
from persson_contact.steady import solve_steady_flow
from persson_contact.transport_table import TransportTable


class AnalyticContact:
    def area(self,p):return np.asarray(p)/1e6


def synthetic_table(policy='physical_zero'):
    p=np.array([1.,1e4,1e5,2e5,4e5])
    h=1e-6/(1+p/1e5)
    return TransportTable(p,h,np.full(5,1e-18),np.full(5,2e-18),AnalyticContact(),NewtonianNoSlip(.1),
        percolation_policy=policy,plateau_pressure_pa=1e5 if policy=='legacy_spc_plateau' else None)


class PressureProfileTests(unittest.TestCase):
    def test_nonuniform_csv_units_and_reference_interface_load(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'p.csv';p.write_text('x_m,contact_pressure_pa,interface_reference_pressure_pa\n0,100,10\n.01,300,20\n.1,200,30\n')
            v=read_pressure_profile(p)
            np.testing.assert_array_equal(v.x_m,[0,.01,.1])
            np.testing.assert_array_equal(v.interface_reference_pressure_pa,[10,20,30])
            mech=make_mechanics({'type':'prescribed_profile','pressure_file':'p.csv'},d,4e6)
            np.testing.assert_allclose(mech.total_reference([0,.01,.1]),[110,320,230],rtol=1e-14)
            self.assertEqual(len(v.metadata['sha256']),64)
            p.write_text('x_m,contact_pressure_pa\n0,100\n.1,300\n.01,200\n')
            with self.assertRaises(ValueError):read_pressure_profile(p)

    def test_original_fortran_format_and_explicit_legacy_resampling(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'profile';p.write_text('FEM example\n.069\n0 .1\n.0001 .2\n.001 .1\n')
            v=read_pressure_profile(p,'fortran')
            self.assertEqual(v.metadata['reference_external_fluid_pressure_pa'],69000.)
            np.testing.assert_allclose(v.contact_pressure_pa,[1e5,2e5,1e5])
            old=read_pressure_profile(p,'fortran','legacy_uniform_zero_padding')
            self.assertEqual(len(old.x_m),5)
            self.assertAlmostEqual(old.x_m[-1],4*.001/3)
            np.testing.assert_array_equal(old.contact_pressure_pa,[0,1e5,2e5,1e5,0])

    def test_hertz_width_and_square_ring(self):
        h=make_mechanics({'type':'oring_hertz','cross_section_diameter_m':.003,'peak_pressure_pa':1e5},'.',4e6)
        self.assertAlmostEqual(h.length_m,.00015)
        np.testing.assert_allclose(h.evaluate([0,h.length_m/2,h.length_m])[0],[0,1e5,0],atol=1e-8)
        s=make_mechanics({'type':'square_ring','length_m':.001,'reference_contact_pressure_pa':1e5},'.',4e6)
        np.testing.assert_array_equal(s.evaluate([0,.0004,.001])[0],[1e5]*3)


class StationaryEquationTests(unittest.TestCase):
    def test_uniform_mobility_analytic_flux_and_reversal(self):
        sigma=1e-18;length=.003
        r=lambda x,p:np.full_like(np.asarray(x),1/sigma)
        for a,b in ((10000.,0.),(0.,10000.),(100.,100.)):
            v=solve_steady_flow(r,0,length,a,b)
            np.testing.assert_allclose(v.flux_per_width_m2_s,sigma*(a-b)/length,rtol=1e-10,atol=0)
            np.testing.assert_allclose(v.fluid_pressure_pa,a+(b-a)*v.x_m/length,rtol=1e-10,atol=1e-9)
            self.assertLess(v.pressure_drop_relative_error,1e-9)

    def test_unloading_dependent_mobility_against_exact_integral(self):
        # sigma(pc)=sigma0*exp(-b*pc), pc=T-pf.
        sigma0=2e-18;b=2e-5;total=1e5;pin=3e4;length=.002
        r=lambda x,p:np.exp(b*(total-p))/sigma0
        v=solve_steady_flow(r,0,length,pin,0,tolerance=1e-7)
        expected=sigma0*np.exp(-b*total)*np.expm1(b*pin)/(b*length)
        np.testing.assert_allclose(v.flux_per_width_m2_s,expected,rtol=2e-7)
        self.assertLess(v.face_flux_relative_spread,1e-5)

    def test_narrow_feature_preserved_by_required_geometry_nodes(self):
        center=.3333;width=1e-5
        resistance=lambda x,p:1+np.maximum(1-abs(x-center)/width,0)/width
        result=solve_steady_flow(resistance,0,1,1,0,
            required_points=[center-width,center,center+width])
        # Unit background integral plus triangular peak area one.
        self.assertAlmostEqual(result.flux_per_width_m2_s,.5,places=9)
        for point in (center-width,center,center+width):self.assertIn(point,result.x_m)

    def test_invalid_or_unresolved_flow_is_not_exported_as_converged(self):
        with self.assertRaises(ValueError):solve_steady_flow(lambda x,p:np.zeros_like(x),0,1,1,0)
        with self.assertRaises(ValueError):solve_steady_flow(lambda x,p:np.full_like(x,np.inf),0,1,1,0)


class SpatialApplicationTests(unittest.TestCase):
    def plan(self,kind='square_ring',coupling='local_unloading'):
        geom={'type':kind,'length_m':.001,'reference_contact_pressure_pa':1e5}
        if kind=='oring_hertz':geom={'type':kind,'length_m':.001,'peak_pressure_pa':1e5}
        return prepare_leakage({'mode':'spatial','geometry':geom,'seal_width_m':.01,
            'inlet_pressures_pa':[0,1e4,2e4],'outlet_pressure_pa':0,'coupling':coupling,
            'percolation_policy':'physical_zero','compute_upper_bound':True},'.',4e6)

    def test_all_three_geometries_and_unloading_state(self):
        table=synthetic_table()
        for kind in ('square_ring','oring_hertz'):
            plan=self.plan(kind)
            plan.outlet_pressure_pa=-20000.
            v=solve_case(plan,table,-10000.)
            np.testing.assert_allclose(v['flux'],1e4/((12*.1)/1e-18*.001),rtol=.001)
            np.testing.assert_allclose(v['contact'],np.maximum(plan.mechanics.total_reference(v['x'])-v['fluid'],0),atol=1e-8)
        with tempfile.TemporaryDirectory() as d:
            path=Path(d)/'p.csv';path.write_text('x_m,contact_pressure_pa\n0,100000\n.0001,150000\n.001,100000\n')
            plan=self.plan();plan.mechanics=make_mechanics({'type':'prescribed_profile','pressure_file':str(path)},d,4e6)
            v=solve_case(plan,table,10000.)
            np.testing.assert_allclose(v['flux'],1e4/(1.2e18*.001),rtol=.001)
            result=run_spatial_leakage(plan,table,d)
            self.assertEqual(len(result['cases']),3)
            self.assertEqual(result['cases'][0]['flow_rate_m3_s'],0)
            self.assertTrue((Path(d)/'leakage_profile_002.csv').is_file())
            for row in result['cases']:
                self.assertLessEqual(row['pdf_bound_min_flow_m3_s'],row['flow_rate_m3_s'])
                self.assertGreaterEqual(row['pdf_bound_max_flow_m3_s'],row['flow_rate_m3_s'])

    def test_blocked_branch_reports_nonunique_pressure(self):
        original=synthetic_table()
        table=TransportTable(original.p,original.h,np.zeros(5),np.zeros(5),AnalyticContact(),NewtonianNoSlip(.1))
        v=solve_case(self.plan(),table,1e4)
        self.assertEqual(v['status'],'blocked_pressure_not_unique')
        self.assertEqual(v['flux'],0)
        self.assertTrue(np.any(np.isnan(v['fluid'])))


    def test_real_fortran_500_kernel_references(self):
        from persson_contact.mechanics import ReferenceMechanics
        fixture=json.loads((Path(__file__).parent/'data/fortran_spatial_leakage.json').read_text())
        table=TransportTable(fixture['pressures_pa'],fixture['mean_gaps_m'],fixture['geometric_conductivity_m3'],
            fixture['geometric_conductivity_m3'],AnalyticContact(),NewtonianNoSlip(.1),interpolation='linear',
            percolation_policy='legacy_spc_plateau',plateau_pressure_pa=fixture['plateau_pressure_pa'])
        self.assertEqual(len(fixture['cases']),6)
        for case in fixture['cases']:
            with self.subTest(geometry=case['geometry'],inlet=case['inlet_pressure_pa']):
                plan=self.plan();plan.solver['tolerance']=1e-6
                plan.mechanics=ReferenceMechanics(np.array(case['x_m']),np.array(case['reference_contact_pressure_pa']),
                    np.array(case['reference_interface_pressure_pa']),case['mechanics_metadata'],'linear')
                v=solve_case(plan,table,case['inlet_pressure_pa'])
                np.testing.assert_allclose(v['flux']*fixture['width_m'],case['fortran_flow_rate_m3_s'],rtol=.002)
                self.assertLess(v['pressure_drop_error'],1e-8)

    def test_spatial_workflow_and_file_preflight(self):
        from parameters import PARAMETERS
        from persson_contact.workflow import run_project
        root=Path(__file__).resolve().parents[1]
        with tempfile.TemporaryDirectory() as d:
            cfg=deepcopy(PARAMETERS);cfg['psd']['path']=str(root/'examples/data/roughness_psd.dat')
            cfg.pop('contact_pressure');cfg['pressures_pa']=[1e4,1e5,2e5]
            cfg['leakage']['geometry']={'type':'square_ring','length_m':.001,'reference_contact_pressure_pa':1e5}
            cfg['leakage']['inlet_pressures_pa']=[0,1e4];cfg['output_directory']='spatial'
            result=run_project(cfg,base_directory=d,run_contact=False,run_leakage=True,
                               compute_gap_pdf=False,save_plots=False)
            out=Path(d)/'spatial'
            self.assertEqual(result['workflow_status'],'complete')
            self.assertEqual(result['leakage']['case_count'],2)
            self.assertGreater(result['leakage_database_extension']['added_nodes'],0)
            used=np.genfromtxt(out/'leakage_transport_table.csv',delimiter=',',names=True,encoding='utf-8')
            self.assertTrue(np.any((used['raw_geometric_conductivity_m3']==0)&(used['used_geometric_conductivity_m3']>0)))
            profile=np.genfromtxt(out/'leakage_profile_001.csv',delimiter=',',names=True,encoding='utf-8')
            self.assertEqual(profile['fluid_pressure_pa'][0],1e4)
            self.assertEqual(profile['fluid_pressure_pa'][-1],0)
            np.testing.assert_allclose(profile['contact_pressure_pa']+profile['fluid_pressure_pa'],1e5)
            cfg['output_directory']='bad';cfg['leakage']['geometry']={'type':'prescribed_profile','pressure_file':'absent.csv'}
            with self.assertRaises(FileNotFoundError):run_project(cfg,base_directory=d,run_leakage=True)
            self.assertFalse((Path(d)/'bad').exists())

    def test_explicit_parameters_and_unsupported_gas(self):
        raw={'mode':'spatial','geometry':{'type':'square_ring','length_m':.001,'reference_contact_pressure_pa':1e5},
             'seal_width_m':.01,'inlet_pressure_pa':1e4}
        v=prepare_leakage(raw,'.',4e6)
        self.assertEqual(v.method,'effective_medium')
        for update in ({'transport':'ballistic_gas'},{'inlet_pressures_pa':[1,2]},
                       {'coupling':'full_fem'},{'spc':1.1}):
            with self.assertRaises((ValueError,NotImplementedError)):prepare_leakage({**raw,**update},'.',4e6)

    def test_preset_selector_does_not_read_inactive_files(self):
        raw={'mode':'spatial','geometry':'square_ring','geometry_parameters':{
            'square_ring':{'length_m':.001,'reference_contact_pressure_pa':1e5},
            'prescribed_profile':{'pressure_file':'not_present.csv'}},
            'seal_width_m':.01,'inlet_pressure_pa':1e4}
        plan=prepare_leakage(raw,'.',4e6)
        self.assertEqual(plan.mechanics.metadata['type'],'square_ring')
        raw['geometry']='prescribed_profile'
        with self.assertRaises(FileNotFoundError):prepare_leakage(raw,'.',4e6)

    def test_low_pressure_and_compatibility_plateau_are_explicit(self):
        table=synthetic_table('legacy_spc_plateau')
        self.assertEqual(table.geometric_conductivity(2e5),table.geometric_conductivity(1e5))
        self.assertGreater(table.mean_gap(.1),table.h[0])
        self.assertEqual(table.resistance(0),0)
        self.assertTrue(np.isnan(table.mean_gap(0)))
        with self.assertRaises(ValueError):table.geometric_conductivity(1e6)


if __name__=='__main__':unittest.main()
