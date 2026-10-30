#!/usr/bin/env python3
"""Validate nonlocal line elasticity, Hertz and coupled Persson dry contact."""
import argparse
from copy import deepcopy
import csv
import hashlib
import json
from pathlib import Path
import subprocess
import sys
import time

import numpy as np
from scipy.integrate import quad

PROJECT=Path(__file__).resolve().parents[1]
sys.path[:0]=[str(PROJECT/'src'),str(PROJECT)]
from persson_contact.core import PerssonModel
from persson_contact.macroscopic import HalfSpaceLine,HertzLine,solve_smooth_contact
from persson_contact.dry_contact import RoughContactLaw,solve_rough_contact
from persson_contact.psd import read_psd
from persson_contact.table import build_contact_table
from persson_contact.workflow import run_project


def write(path,rows):
    with path.open('w',newline='') as out:
        w=csv.DictWriter(out,fieldnames=list(rows[0]),lineterminator='\n');w.writeheader();w.writerows(rows)


def make_table(model,ref,points=1201,high=None):
    focus=ref.peak_pressure_pa*1e-3
    p=np.unique(np.r_[np.geomspace(1e-100,focus,points//3),np.geomspace(focus,high or 2e6,points-points//3+1)])
    return build_contact_table(model,p)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--no-plots',action='store_true');args=parser.parse_args()
    out=args.output;out.mkdir(parents=True,exist_ok=False);start=time.monotonic()
    ref=HertzLine(.005,100.,4e6);a=ref.half_width_m;scale=a*a/(2*ref.radius_m)
    psd_path=PROJECT/'examples/data/roughness_psd.dat';model=PerssonModel(read_psd(psd_path),4e6,.4)
    table=make_table(model,ref);law=RoughContactLaw(table)
    checks=[]
    def check(name,value,limit):
        checks.append({'name':name,'value':float(value),'limit':limit,'passed':bool(value<=limit)})
    elastic=HalfSpaceLine(3*a,65,4e6)
    prescribed=1e5*(1+.2*np.cos(np.pi*elastic.x_m/(3*a)))
    u=elastic.displacement(prescribed);given=[]
    for i in [0,7,16,32,48,57,64]:
        x=elastic.x_m[i];expected=0.
        for j,(l,r) in enumerate(zip(elastic.edges_m[:-1],elastic.edges_m[1:])):
            expected+=prescribed[j]*quad(lambda s:np.log(elastic.reference_length_m/abs(x-s)),l,r,
                points=[x] if l<x<r else None,epsabs=1e-15,epsrel=1e-11)[0]
        expected*=2/(np.pi*4e6)
        given.append({'x_m':x,'displacement_m':u[i],'adaptive_reference_m':expected,'relative_error':abs(u[i]/expected-1)})
    write(out/'given_pressure_displacement.csv',given)
    check('log_kernel_against_independent_adaptive_integral',max(r['relative_error'] for r in given),1e-10)
    smooth=[];rough=[];solutions={}
    for n in [65,129,257,513,1025]:
        el=HalfSpaceLine(3*a,n,4e6);s=solve_smooth_contact(el,ref);r=solve_rough_contact(el,ref,law)
        smooth.append({'cells':n,'pressure_L1_relative_error':np.sum(abs(s.pressure_pa-ref.cell_pressure(el.edges_m)))*el.dx_m/100,
            'peak_pressure_relative_error':abs(max(s.pressure_pa)/ref.peak_pressure_pa-1),
            'gap_max_error_over_hertz_scale':np.max(abs(s.mean_gap_m-ref.gap(s.x_m)))/scale,
            'load_relative_error':s.load_relative_error,'equilibrium_residual':s.equilibrium_relative_residual})
        actual=model.mean_gap(r.pressure_pa)
        rough.append({'cells':n,'peak_pressure_pa':max(r.pressure_pa),'minimum_mean_gap_m':r.minimum_gap_m,
            'load_relative_error':r.load_relative_error,'geometry_relative_residual':r.equilibrium_relative_residual,
            'module2_gap_error_over_hrms':np.max(abs(actual-r.mean_gap_m))/model.psd.hrms,'edge_pressure_ratio':r.edge_pressure_ratio})
        solutions[n]=(el,s,r)
    write(out/'smooth_convergence.csv',smooth);write(out/'rough_convergence.csv',rough)
    check('smooth_finest_pressure_L1',smooth[-1]['pressure_L1_relative_error'],.001)
    check('smooth_finest_gap_profile',smooth[-1]['gap_max_error_over_hertz_scale'],.002)
    check('smooth_load_and_complementarity',max(max(v['load_relative_error'],v['equilibrium_residual']) for v in smooth),1e-8)
    check('rough_load_and_geometry',max(max(v['load_relative_error'],v['geometry_relative_residual']) for v in rough),1e-8)
    check('rough_default_peak_vs_finest',abs(rough[2]['peak_pressure_pa']/rough[-1]['peak_pressure_pa']-1),.001)
    check('rough_default_min_gap_vs_finest',abs(rough[2]['minimum_mean_gap_m']/rough[-1]['minimum_mean_gap_m']-1),.001)
    check('rough_energy_law_vs_tabulation',max(v['module2_gap_error_over_hrms'] for v in rough),1e-4)
    domains=[]
    for fac in [2.,3.,4.]:
        n=int(round(257*fac/3));n+=1-n%2
        el=HalfSpaceLine(fac*a,n,4e6);r=solve_rough_contact(el,ref,law)
        domains.append({'half_width_over_hertz_a':fac,'cells':n,'peak_pressure_pa':max(r.pressure_pa),
            'minimum_mean_gap_m':r.minimum_gap_m,'edge_pressure_ratio':r.edge_pressure_ratio})
    write(out/'domain_convergence.csv',domains)
    check('domain_min_gap_variation',max(r['minimum_mean_gap_m'] for r in domains)/min(r['minimum_mean_gap_m'] for r in domains)-1,.001)
    tables=[];el=solutions[257][0]
    for count in [301,601,1201,2401]:
        t=make_table(model,ref,count);r=solve_rough_contact(el,ref,RoughContactLaw(t))
        tables.append({'table_points':len(t.pressures_pa),'peak_pressure_pa':max(r.pressure_pa),'minimum_mean_gap_m':r.minimum_gap_m,
            'energy_gap_error_over_hrms':np.max(abs(model.mean_gap(r.pressure_pa)-r.mean_gap_m))/model.psd.hrms})
    write(out/'table_convergence.csv',tables)
    check('default_table_peak_vs_finer',abs(tables[2]['peak_pressure_pa']/tables[-1]['peak_pressure_pa']-1),.0001)
    base=solutions[257][2];alternate=solve_rough_contact(HalfSpaceLine(3*a,257,4e6,12*a),ref,law)
    gauge={'reference_ratio':2.,'maximum_pressure_change_over_peak':float(np.max(abs(base.pressure_pa-alternate.pressure_pa))/max(base.pressure_pa)),
        'maximum_gap_change_over_hrms':float(np.max(abs(base.mean_gap_m-alternate.mean_gap_m))/model.psd.hrms),
        'offset_change_m':alternate.rigid_offset_m-base.rigid_offset_m,
        'analytic_offset_change_m':-2*100/(np.pi*4e6)*np.log(2)}
    check('reference_length_pressure_invariance',gauge['maximum_pressure_change_over_peak'],1e-8)
    check('reference_length_gap_invariance',gauge['maximum_gap_change_over_hrms'],1e-8)
    scaled_model=PerssonModel(model.psd,2.4e6,.4)
    scaled_table=build_contact_table(scaled_model,table.pressures_pa*.6)
    scaled=solve_rough_contact(HalfSpaceLine(3*a,257,2.4e6),HertzLine(.005,60.,2.4e6),RoughContactLaw(scaled_table))
    check('load_modulus_scaling_pressure',float(np.max(abs(scaled.pressure_pa-base.pressure_pa*.6))/max(scaled.pressure_pa)),1e-8)
    check('load_modulus_scaling_gap',float(np.max(abs(scaled.mean_gap_m-base.mean_gap_m))/model.psd.hrms),1e-8)
    limits=[];smooth_ref=solutions[257][1]
    for amplitude in [1.,.5,.25,.1]:
        m=PerssonModel(read_psd(psd_path,scale_c=amplitude**2),4e6,.4)
        t=make_table(m,ref,high=2*ref.peak_pressure_pa);r=solve_rough_contact(el,ref,RoughContactLaw(t))
        limits.append({'roughness_amplitude_multiplier':amplitude,'hrms_m':m.psd.hrms,
            'pressure_L1_difference_from_smooth':np.sum(abs(r.pressure_pa-smooth_ref.pressure_pa))*el.dx_m/100,
            'minimum_mean_gap_m':r.minimum_gap_m,'truncated_far_tail_cells':int(np.sum(r.pressure_pa==0)),
            'fixed_gap_tail_load_bound_n_m':float(t.pressures_pa[0]*np.ptp(el.edges_m)),
            'geometry_residual':r.equilibrium_relative_residual,'load_error':r.load_relative_error})
    write(out/'roughness_to_smooth.csv',limits)
    check('roughness_reduction_monotonically_approaches_smooth',0. if np.all(np.diff([r['pressure_L1_difference_from_smooth'] for r in limits])<0) else 1.,0.)
    check('roughness_tenth_amplitude_pressure_difference',limits[-1]['pressure_L1_difference_from_smooth'],.005)
    from parameters import PARAMETERS
    for kind in ['smooth','persson']:
        p=deepcopy(PARAMETERS);p['psd']['path']=str(psd_path);p['macro_contact']['contact_model']=kind
        p['output_directory']=str((out/kind).resolve())
        run_project(p,base_directory=PROJECT,run_macro_contact=True,save_plots=not args.no_plots)
    write(out/'default_rough_profile.csv',[{'x_m':x,'pressure_pa':p,'mean_gap_m':h,'displacement_m':u,'relative_contact_area':a}
        for x,p,h,u,a in zip(base.x_m,base.pressure_pa,base.mean_gap_m,base.displacement_m,base.relative_contact_area)])
    summary={'base_commit':subprocess.check_output(['git','rev-parse','HEAD'],cwd=PROJECT,text=True).strip(),
        'source_hashes':{str(p.relative_to(PROJECT)):hashlib.sha256(p.read_bytes()).hexdigest()
            for p in sorted((PROJECT/'src/persson_contact').glob('*.py'))},
        'psd_sha256':hashlib.sha256(psd_path.read_bytes()).hexdigest(),'checks':checks,
        'passed':all(c['passed'] for c in checks),'gauge_check':gauge,
        'hertz_half_width_m':ref.half_width_m,'hertz_peak_pressure_pa':ref.peak_pressure_pa,
        'elapsed_seconds':time.monotonic()-start,
        'scope':'Independent dry macro elasticity and module-2+6 validation only. No fluid pressure, friction, finite-thickness or whole-ring structural model.'}
    (out/'summary.json').write_text(json.dumps(summary,indent=2,allow_nan=False)+'\n')
    if not args.no_plots:plot(out,smooth,rough,domains,tables,limits)
    print(json.dumps({'passed':summary['passed'],'checks':len(checks),'elapsed_seconds':summary['elapsed_seconds']}))
    return 0 if summary['passed'] else 2


def plot(out,smooth,rough,domains,tables,limits):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig,axes=plt.subplots(2,2,figsize=(11,7.5),layout='constrained')
    axes[0,0].loglog([r['cells'] for r in smooth],[r['pressure_L1_relative_error'] for r in smooth],'o-')
    axes[0,0].set(xlabel='Elastic cells',ylabel='Smooth pressure L1 error vs Hertz')
    axes[0,1].semilogx([r['cells'] for r in rough],[1e6*(r['minimum_mean_gap_m']/rough[-1]['minimum_mean_gap_m']-1) for r in rough],'o-')
    axes[0,1].set(xlabel='Elastic cells',ylabel='Minimum gap change from finest [ppm]')
    axes[1,0].loglog([r['table_points'] for r in tables],[r['energy_gap_error_over_hrms'] for r in tables],'o-')
    axes[1,0].set(xlabel='Persson table points',ylabel='Contact-table gap error / RMS height')
    axes[1,1].loglog([r['roughness_amplitude_multiplier'] for r in limits],[r['pressure_L1_difference_from_smooth'] for r in limits],'o-')
    axes[1,1].set(xlabel='Roughness amplitude multiplier',ylabel='Pressure L1 difference from smooth')
    for ax in axes.flat:ax.grid(alpha=.25)
    fig.suptitle('Module 6 and module 2+6: independent dry-contact validation')
    fig.savefig(out/'convergence.png',dpi=160);plt.close(fig)


if __name__=='__main__':raise SystemExit(main())
