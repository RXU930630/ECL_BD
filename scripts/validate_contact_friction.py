#!/usr/bin/env python3
"""Validate dry real-area sliding friction; demonstration strengths are not experimental data."""
import argparse
from copy import deepcopy
import csv
import hashlib
import json
from pathlib import Path
import subprocess
import sys

import numpy as np

PROJECT=Path(__file__).resolve().parents[1]
sys.path[:0]=[str(PROJECT/'src'),str(PROJECT)]
from persson_contact.contact_friction import ConstantShearStrength,LogGaussianShearStrength,CustomShearStrength,evaluate_contact_friction
from persson_contact.macroscopic import HalfSpaceLine,HertzLine
from persson_contact.dry_contact import RoughContactLaw,solve_rough_contact
from persson_contact.core import PerssonModel
from persson_contact.psd import read_psd
from persson_contact.table import build_contact_table
from persson_contact.workflow import run_project


def write(path,rows):
    with path.open('w',newline='') as stream:
        w=csv.DictWriter(stream,fieldnames=list(rows[0]),lineterminator='\n');w.writeheader();w.writerows(rows)


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--no-plots',action='store_true');args=parser.parse_args()
    out=args.output;out.mkdir(parents=True,exist_ok=False)
    checks=[]
    def check(name,value,limit):checks.append({'name':name,'value':float(value),'limit':limit,'passed':bool(value<=limit)})
    def synthetic(law,v=.2,**kw):
        return evaluate_contact_friction([0,1,3,4],[.25,.5,0],[2,3,0],[.1,.1,.2],v,law,**kw)
    law=ConstantShearStrength(8);a=synthetic(law);b=synthetic(law,-.2);zero=synthetic(law,0)
    check('analytic_nonuniform_real_area_force',abs(a.summary['friction_force_magnitude_per_width_n_m']/10-1),1e-13)
    check('analytic_normalization',abs(a.summary['kinetic_contact_friction_coefficient']/1.25-1),1e-13)
    check('action_reaction',abs(a.summary['upper_friction_force_per_width_n_m']+a.summary['lower_friction_force_per_width_n_m']),1e-13)
    check('reversal_upper_force',abs(a.summary['upper_friction_force_per_width_n_m']+b.summary['upper_friction_force_per_width_n_m']),1e-13)
    check('positive_dissipation',max(0.,-a.summary['dissipated_power_per_width_w_m'],-b.summary['dissipated_power_per_width_w_m']),0.)
    check('zero_speed_sliding_force',zero.summary['friction_force_magnitude_per_width_n_m'],0.)
    check('zero_speed_static_coefficient_not_invented',0. if zero.summary['kinetic_contact_friction_coefficient'] is None and zero.summary['static_friction_limit_n_m'] is None else 1.,0.)
    callback=CustomShearStrength(lambda s,p:2+.3*s.real_contact_pressure_pa+2*(s.temperature_k-298.15))
    c=synthetic(callback,temperature_k=299.15)
    check('pressure_temperature_callback_integral',abs(c.summary['friction_force_magnitude_per_width_n_m']/(4*1.25+.3*8)-1),1e-13)
    gaussian=LogGaussianShearStrength(2,10,1e-3,2)
    check('log_gaussian_peak',abs(synthetic(gaussian,1e-3).profile['real_contact_shear_strength_pa'][0]/12-1),1e-13)
    left=synthetic(gaussian,1e-5).profile['real_contact_shear_strength_pa'][0]
    right=synthetic(gaussian,1e-1).profile['real_contact_shear_strength_pa'][0]
    check('log_gaussian_symmetric_decades',abs(left/right-1),1e-13)
    ref=HertzLine(.005,100.,4e6);psd_path=PROJECT/'examples/data/roughness_psd.dat'
    model=PerssonModel(read_psd(psd_path),4e6,.4)
    focus=ref.peak_pressure_pa*1e-3
    table=build_contact_table(model,np.unique(np.r_[np.geomspace(1e-100,focus,400),np.geomspace(focus,2e6,802)]))
    resolution=[];base=None
    for cells in [65,129,257,513]:
        el=HalfSpaceLine(1.2e-3,cells,4e6);normal=solve_rough_contact(el,ref,RoughContactLaw(table))
        r=evaluate_contact_friction(normal.cell_edges_m,normal.relative_contact_area,normal.pressure_pa,normal.mean_gap_m,.01,ConstantShearStrength(1e5))
        resolution.append({'cells':cells,**r.summary})
        if cells==257:base=normal
    write(out/'spatial_convergence.csv',resolution)
    check('default_coefficient_change_to_finer_mesh',abs(resolution[2]['kinetic_contact_friction_coefficient']/resolution[-1]['kinetic_contact_friction_coefficient']-1),.001)
    widths=[]
    for width in [None,.01,.02]:
        r=evaluate_contact_friction(base.cell_edges_m,base.relative_contact_area,base.pressure_pa,base.mean_gap_m,.01,ConstantShearStrength(1e5),transverse_width_m=width)
        widths.append(r.summary)
    write(out/'width_scaling.csv',widths)
    check('width_doubles_total_force',abs(widths[2]['friction_force_magnitude_n']/(2*widths[1]['friction_force_magnitude_n'])-1),1e-13)
    velocities=np.r_[-.01,0.,np.geomspace(1e-8,1.,49)]
    sweep=[]
    for name,law in [('constant',ConstantShearStrength(1e5)),('log_gaussian',LogGaussianShearStrength(1e4,1e5,1e-3,1.5))]:
        for v in velocities:
            r=evaluate_contact_friction(base.cell_edges_m,base.relative_contact_area,base.pressure_pa,base.mean_gap_m,v,law)
            sweep.append({'model':name,**r.summary})
    write(out/'velocity_sweep.csv',sweep)
    moving=[r for r in sweep if r['model']=='constant' and r['speed_m_s']>0]
    expected=1e5*resolution[2]['real_contact_area_per_width_m']
    check('constant_strength_speed_independence',max(abs(r['friction_force_magnitude_per_width_n_m']/expected-1) for r in moving),1e-13)
    from parameters import PARAMETERS
    params=deepcopy(PARAMETERS);params['psd']['path']=str(psd_path)
    params['contact_friction']['model_parameters']['constant_shear_strength']['shear_strength_pa']=1e5
    params['contact_friction']['transverse_width_m']=.01
    params['output_directory']=str((out/'example_constant').resolve())
    run_project(params,base_directory=PROJECT,run_contact_friction=True,save_plots=not args.no_plots)
    metadata={'base_commit':subprocess.check_output(['git','rev-parse','HEAD'],cwd=PROJECT,text=True).strip(),
        'source_hashes':{str(p.relative_to(PROJECT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted((PROJECT/'src/persson_contact').glob('*.py'))},
        'psd_sha256':hashlib.sha256(psd_path.read_bytes()).hexdigest(),
        'demonstration_only':True,'experimental_shear_strength_supplied':False,
        'demonstration_strength_pa':1e5,'checks':checks,'passed':all(r['passed'] for r in checks),
        'scope':'Dry solid-contact sliding friction on elastic normal state. No fluid friction, static solve or viscoelastic hysteresis.'}
    (out/'summary.json').write_text(json.dumps(metadata,indent=2,allow_nan=False)+'\n')
    print(json.dumps({'passed':metadata['passed'],'checks':len(checks)}))
    return 0 if metadata['passed'] else 2


if __name__=='__main__':raise SystemExit(main())
