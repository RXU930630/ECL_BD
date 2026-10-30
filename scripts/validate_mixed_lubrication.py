#!/usr/bin/env python3
"""Audit a completed mixed run; reuse its immutable raw samples, not a fitted result.

The workflow itself must first run to completion. This audit checks the PSD
hash, rebuilds the normal constitutive law, varies the spatial domain/grid and
subsamples the transport table. It also compares a rigid smooth cavitating
cylinder to independent quadrature, and exports shear-cutoff sensitivity.
"""
import argparse
from copy import deepcopy
import csv
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.integrate import quad
from scipy.interpolate import PchipInterpolator
from scipy.optimize import brentq

from persson_contact.dry_contact import RoughContactLaw,solve_rough_contact
from persson_contact.fluid_friction import wall_tractions,scale_inverse_gap
from persson_contact.macroscopic import HalfSpaceLine
from persson_contact.macro_output import _csv
from persson_contact.mixed import solve_film,solve_mixed_continuation
from persson_contact.mixed_output import prepare_mixed,make_contact_database
from persson_contact.mixed_transport import MixedTransport
from persson_contact.core import PerssonModel


def read_rows(path):
    def value(text):
        if text=='':return None
        if text in ('True','False'):return text=='True'
        return float(text)
    with path.open() as stream:return [{k:value(v) for k,v in r.items()} for r in csv.DictReader(stream)]


def sampled_transport(model,law,interface,rows,rough,critical,cutoff):
    t=MixedTransport.__new__(MixedTransport)
    t.model=model;t.contact_law=law;t.eta=interface.viscosity_pa_s;t.rms=model.psd.hrms
    t.rough_surface=rough;t.minimum_shear_gap_m=cutoff
    t.hc=float(model.mean_gap(model.pressure_for_relative_contact_area(critical)))
    t._set_samples(rows)
    return t


def metrics(v,elastic,hertz,transport,speed,friction):
    h=v.mean_gap_m;pf=v.fluid_pressure_pa;dx=elastic.dx_m
    gradients=np.r_[pf[0]/(dx/2),np.diff(pf)/dx,-pf[-1]/(dx/2)]
    tr=wall_tractions(h,transport.factors(h),transport.eta,(gradients[:-1]+gradients[1:])/2,
        speed,0.,rough_surface=transport.rough_surface,saturation=v.saturation)
    from persson_contact.contact_friction import evaluate_contact_friction
    a=transport.model.area(v.contact_pressure_pa)
    solid=evaluate_contact_friction(elastic.edges_m,a,v.contact_pressure_pa,h,-speed,friction.law)
    fluid=float(np.sum(tr['lower_total_pa'])*dx)
    fs=solid.summary['lower_friction_force_per_width_n_m']
    return {'cells':len(elastic.x_m),'domain_half_width_m':float(elastic.edges_m[-1]),
        'fluid_load_fraction':float(np.sum(pf)*dx/hertz.line_load_n_m),
        'minimum_gap_m':float(min(h)),'flux_per_width_m2_s':float(np.mean(v.face_flux_m2_s)),
        'solid_force_n_m':fs,'fluid_force_n_m':fluid,
        'coefficient':float(-np.sign(speed)*(fs+fluid)/hertz.line_load_n_m),
        'load_relative_error':v.load_relative_error,'geometry_relative_residual':v.geometry_relative_residual,
        'mass_relative_residual':v.mass_relative_residual,'solver_residual':v.residual}


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--run-directory',required=True,type=Path)
    parser.add_argument('--output',required=True,type=Path)
    parser.add_argument('--base-directory',type=Path,default=Path(__file__).resolve().parents[1])
    args=parser.parse_args();base=args.base_directory.resolve();out=args.output.resolve()
    if out.exists():raise FileExistsError('Use a new audit output directory')
    diag=json.loads((args.run_directory/'mixed_diagnostics.json').read_text())
    parameters=diag['parameters_snapshot']
    (raw,speeds,pin,pout,tol,limit,points,rough,interface,macro,solver,hertz,
        elastic,friction,opts,flow)=prepare_mixed(parameters,base)
    if pin!=0 or pout!=0:raise ValueError('This audit is for zero-pressure fully flooded reservoir examples')
    model,table,psd_sha=make_contact_database(parameters,base,macro,solver,hertz,elastic,opts)
    if psd_sha!=diag['psd_sha256']:raise ValueError('PSD changed since source run')
    law=RoughContactLaw(table)
    rows=read_rows(args.run_directory/'mixed_transport_table.csv')
    cutoff=parameters['fluid_friction']['minimum_shear_gap_m']
    t=sampled_transport(model,law,interface,rows,rough,opts.critical_relative_contact_area,cutoff)
    speed=.01
    checks={};grids=[];states={}
    for n in (129,257,513):
        el=HalfSpaceLine(macro['domain_half_width_m'],n,hertz.effective_modulus_pa)
        dry=solve_rough_contact(el,hertz,law)
        v,steps=solve_mixed_continuation(el,hertz,law,t.at_velocity,speed,dry,tolerance=tol,max_iterations=limit)
        states[n]=(el,v);record=metrics(v,el,hertz,t,speed,friction);record['continuation_states']=len(steps)
        grids.append(record);print('grid',record,flush=True)
    coarse,fine=grids[1:]
    grid_change={key:abs(fine[key]/coarse[key]-1) for key in ('coefficient','fluid_load_fraction','minimum_gap_m','flux_per_width_m2_s')}
    checks['spatial_257_to_513_all_reported_changes_below_3percent']=max(grid_change.values())<.03
    checks['coupled_load_geometry_mass_residuals']=all(r['load_relative_error']<1e-7 and r['geometry_relative_residual']<1e-7 and r['mass_relative_residual']<1e-6 for r in grids)
    el,v=states[257]
    other=HalfSpaceLine(1.5*macro['domain_half_width_m'],385,hertz.effective_modulus_pa)
    dry=solve_rough_contact(other,hertz,law)
    larger,_=solve_mixed_continuation(other,hertz,law,t.at_velocity,speed,dry,tolerance=tol,max_iterations=limit)
    domain=[grids[1],metrics(larger,other,hertz,t,speed,friction)]
    domain_change=abs(domain[1]['coefficient']/domain[0]['coefficient']-1)
    checks['domain_enlargement_coefficient_change_below_3percent']=domain_change<.03
    # Coarsen raw samples while preserving the endpoints and explicit blend/critical anchors.
    indices=set(range(0,len(rows),2))|{len(rows)-1}
    for anchor in (t.hc,3*t.rms,4*t.rms):indices.add(int(np.argmin(abs(t.gaps-anchor))))
    lower=sampled_transport(model,law,interface,[rows[i] for i in sorted(indices)],rough,opts.critical_relative_contact_area,cutoff)
    dry=solve_rough_contact(el,hertz,law)
    vc,coarse_steps=solve_mixed_continuation(el,hertz,law,lower.at_velocity,speed,dry,tolerance=tol,max_iterations=limit)
    transport_records=[{'transport_points':len(lower.rows),**metrics(vc,el,hertz,lower,speed,friction),'continuation_states':len(coarse_steps)},
        {'transport_points':len(t.rows),**grids[1]}]
    table_change=abs(transport_records[0]['coefficient']/transport_records[1]['coefficient']-1)
    checks['transport_coarsening_coefficient_change_below_2percent']=table_change<.02
    sensitivity=[]
    for factor,key in ((.5,'phi_f_half_shear_cutoff'),(1.,'phi_f'),(2.,'phi_f_double_shear_cutoff')):
        ts=sampled_transport(model,law,interface,rows,rough,opts.critical_relative_contact_area,cutoff*factor)
        ts._f=PchipInterpolator(ts.gaps,[r[key] for r in rows],extrapolate=False)
        sensitivity.append({'minimum_shear_gap_m':cutoff*factor,**metrics(v,el,hertz,ts,speed,friction)})
    checks['cutoff_changes_friction_without_changing_normal_solution']=len(set(r['fluid_load_fraction'] for r in sensitivity))==1
    checks['zero_boundary_fluid_force_opposes_moving_plane']=all(r['fluid_force_n_m']*speed<=0 for r in sensitivity)
    # Independent rigid-cylinder JFO quadrature, no contact/elasticity/EMT used.
    L=.001;R=.005;h0=1e-6;eta=.1;U=.01
    gap=lambda x:h0+x*x/(2*R)
    rupture=brentq(lambda r:quad(lambda x:(gap(x)-gap(r))/gap(x)**3,-L,r,epsabs=.1)[0],0,L,xtol=1e-15)
    q_exact=U/2*gap(rupture);jfo=[]
    for n in (65,129,257,513):
        dx=2*L/n;x=-L+(np.arange(n)+.5)*dx;h=gap(x)
        def transport(h):return h**3/(12*eta),h*h/(4*eta),U*h/2,np.full_like(h,U/2)
        rigid=solve_film(h,dx,transport,pressure_scale_pa=1e5)
        exact_p=np.array([12*eta*U/2*quad(lambda z:(gap(z)-gap(rupture))/gap(z)**3,-L,float(xx),epsabs=.1)[0] if xx<rupture else 0. for xx in x])
        jfo.append({'cells':n,'flux_per_width_m2_s':float(np.mean(rigid.face_flux_m2_s)),
            'relative_flux_error':abs(float(np.mean(rigid.face_flux_m2_s))/q_exact-1),
            'pressure_L1_relative_error':float(np.sum(abs(rigid.pressure_pa-exact_p))/np.sum(exact_p)),
            'mass_relative_residual':rigid.conservation_relative_residual})
    checks['independent_JFO_finest_flux_and_pressure_errors_below_0p2percent']=max(jfo[-1]['relative_flux_error'],jfo[-1]['pressure_L1_relative_error'])<.002
    refined=PerssonModel(model.psd,model.effective_modulus_pa,model.gamma,n_q=opts.n_q,
        pdf_scale_method='conservative',pdf_log_step=opts.pdf_log_step/2,
        pdf_max_relative_contact_area=opts.pdf_max_relative_contact_area,
        pdf_min_relative_contact_area=opts.pdf_min_relative_contact_area,pdf_min_u1_m=opts.pdf_min_u1_m)
    scale_records=[]
    for p in (1e4,1e5):
        a=scale_inverse_gap(model,p,cutoff);b=scale_inverse_gap(refined,p,cutoff)
        scale_records.append({'pressure_pa':p,'coarse_inverse_gap_m_inv':a['inverse_gap_m_inv'],
            'fine_inverse_gap_m_inv':b['inverse_gap_m_inv'],
            'relative_change':abs(b['inverse_gap_m_inv']/a['inverse_gap_m_inv']-1)})
    checks['scale_step_halving_inverse_moment_change_below_3percent']=max(r['relative_change'] for r in scale_records)<.03
    positive=next(c for c in diag['cases'] if c['lower_wall_velocity_m_s']==.01)
    checks['saved_case_reproduced_at_257_cells']=(
        abs(positive['fluid_load_fraction']-grids[1]['fluid_load_fraction'])<1e-9 and
        abs(positive['lower_resisting_friction_coefficient']-grids[1]['coefficient'])<1e-9)
    negative=next((c for c in diag['cases'] if c['lower_wall_velocity_m_s']==-.01),None)
    if negative is not None:
        checks['signed_velocity_reversal']=abs(positive['lower_total_force_n_m']+negative['lower_total_force_n_m'])/abs(positive['lower_total_force_n_m'])<1e-7
    checks['all_case_cavitation_complementarity']=all(c['complementarity_max_pa']==0 and 0<=c['minimum_saturation']<=1 for c in diag['cases'])
    paths=['src/persson_contact/'+name+'.py' for name in ('mixed','mixed_transport','mixed_output','fluid_friction','macroscopic','dry_contact','core','flow','shear','interface')]
    paths+=['scripts/validate_mixed_lubrication.py','scripts/verify_shear_friction_reference.py','parameters.py','main.py']
    report={'purpose':'Programming verification with explicitly illustrative tau=0.1MPa and shear cutoff=1nm; no experimental calibration.',
        'checks':checks,'all_checks_pass':all(checks.values()),'spatial_relative_changes':grid_change,
        'domain_coefficient_relative_change':domain_change,'transport_coefficient_relative_change':table_change,
        'independent_JFO_reference':{'rupture_m':rupture,'flux_per_width_m2_s':q_exact},
        'case_input':parameters,'source_sha256':{path:hashlib.sha256((base/path).read_bytes()).hexdigest() for path in paths},
        'source_run_diagnostics_sha256':hashlib.sha256((args.run_directory/'mixed_diagnostics.json').read_bytes()).hexdigest(),
        'source_transport_table_sha256':hashlib.sha256((args.run_directory/'mixed_transport_table.csv').read_bytes()).hexdigest()}
    out.mkdir(parents=True)
    _csv(out/'spatial_convergence.csv',grids);_csv(out/'domain_convergence.csv',domain)
    _csv(out/'transport_convergence.csv',transport_records);_csv(out/'shear_cutoff_sensitivity.csv',sensitivity)
    _csv(out/'rigid_JFO_convergence.csv',jfo);_csv(out/'scale_step_convergence.csv',scale_records)
    (out/'summary.json').write_text(json.dumps(report,indent=2,ensure_ascii=False,allow_nan=False)+'\n')
    print(json.dumps(checks,indent=2),flush=True)
    if not all(checks.values()):raise AssertionError('Mixed verification has failed checks; inspect summary.json')


if __name__=='__main__':main()
