#!/usr/bin/env python3
"""Independent convergence and physical-limit audit; no Fortran executable needed."""
import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
import subprocess
import sys
import time

import numpy as np
from scipy.integrate import quad
from scipy.special import erf,ndtr

PROJECT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(PROJECT/'src'))
from persson_contact.core import PerssonModel,separation_kernel
from persson_contact.flow import FlowOptions,bruggeman_conductivity,homogenize_bins,pressure_flow
from persson_contact.interface import NewtonianNoSlip
from persson_contact.moving_flow import solve_fixed_contact_motion
from persson_contact.psd import read_psd
from persson_contact.shear import ShearFlowLaw,ShearOptions,shear_length


def write_csv(path,rows):
    with path.open('w',newline='') as stream:
        writer=csv.DictWriter(stream,fieldnames=list(rows[0]),lineterminator='\n')
        writer.writeheader();writer.writerows(rows)


def energy_reference(psd,pressure,modulus=4e6,gamma=.4):
    """Adaptive q integration, split at all input knots; independent of Gauss rule."""
    def f(t):
        q=psd.q0*np.exp(t*t)
        i3=psd.moment(3,psd.q0,min(q,psd.q1))
        return 2*t*q**3*psd.evaluate(q)/np.sqrt(i3)*separation_kernel(
            pressure/(modulus*np.sqrt(np.pi*i3)),gamma)
    q=10.**psd.log10_q
    q=q[(q>psd.q0)&(q<psd.q1)]
    bounds=np.r_[0.,np.sqrt(np.log(q/psd.q0)),np.sqrt(np.log(psd.q1/psd.q0))]
    pieces=[quad(f,a,b,epsabs=1e-17,epsrel=2e-10,limit=150)
            for a,b in zip(bounds[:-1],bounds[1:])]
    return math.fsum(v for v,e in pieces),math.fsum(e for v,e in pieces)


def physical_limits(psd):
    checks=[]
    def check(name,error,tolerance):
        checks.append({'name':name,'error':float(error),'tolerance':tolerance,'passed':bool(error<=tolerance)})
    # Independent adaptive integral of the energy kernel.
    kernel=[]
    for gamma in [0.,.4,1.]:
        for z in [1e-5,.01,.1,1.,4.]:
            expected=quad(lambda x:(gamma+3*(1-gamma)*erf(x)**2)*np.exp(-x*x)/x,
                          z,15.,epsabs=1e-15,epsrel=2e-11)[0]
            value=separation_kernel(z,gamma)
            kernel.append({'gamma':gamma,'z':z,'relative_error':abs(value/expected-1)})
    check('energy_kernel_vs_adaptive_integral',max(r['relative_error'] for r in kernel),2e-7)
    # Exact two-phase effective-medium solutions, including near-zero/infinite roots.
    root=[]
    for area in [.1,.49,.5-1e-9,.5-1e-12,.5-1e-14]:
        for scale in [1e-200,1.,1e200]:
            value,_,_=bruggeman_conductivity([scale],[1-area],area,2.)
            exact=scale*((1-area)-area)
            root.append({'zero_probability':area,'local_k':scale,'relative_error':abs(value/exact-1)})
    check('binary_closed_medium_near_percolation',max(r['relative_error'] for r in root),2e-11)
    errors=[]
    for infinite in [.1,.49,.5-1e-12,.5-1e-14]:
        value,_,_=bruggeman_conductivity([1.],[1-infinite],0,2.,infinite)
        errors.append(abs(value*((1-infinite)-infinite)-1))
    check('binary_infinite_medium_near_threshold',max(errors),2e-11)
    value,_,_=bruggeman_conductivity([1e-20,1e-16],[.5,.5],0,2.)
    check('two_equal_positive_phases_geometric_mean',abs(value/1e-18-1),2e-11)
    # Prescribed narrow Gaussian gap, no contact: phi=1-3/2*(rms/h)^2+O((rms/h)^4).
    asymptotic=[]
    for ratio in [5.,10.,20.,40.]:
        h=ratio*1e-6;sigma=1e-6
        edges=np.linspace(0.,h+10*sigma,16385)
        cdf=ndtr((edges-h)/sigma)-ndtr((-edges-h)/sigma)
        r=homogenize_bins(edges,cdf,0.,h,NewtonianNoSlip(.1),.5)
        coefficient=(1-r['factor'])*ratio**2
        asymptotic.append({'gap_over_rms':ratio,'phi_p':r['factor'],'coefficient':coefficient})
    check('large_gap_gaussian_second_order_coefficient',abs(asymptotic[-1]['coefficient']/1.5-1),.01)
    # Analytic smooth parallel plates and odd pressure-flow response.
    fluid=NewtonianNoSlip(.1);h=2e-6;gradient=3e7
    exact=-h**3*gradient/(12*.1)
    check('parallel_plate_poiseuille',abs(fluid.flux(h,gradient)/exact-1),1e-13)
    check('pressure_reversal',abs(fluid.flux(h,gradient)+fluid.flux(h,-gradient))/abs(exact),1e-13)
    check('positive_pressure_dissipation',max(0.,fluid.flux(h,gradient)*gradient),0.)
    # Shear amplitude bounds, branch continuity, smooth limit and Galilean shift.
    heights=np.geomspace(1e-12,1e-3,500);rms=psd.hrms;hc=1e-7
    ls=shear_length(heights,rms,hc)
    check('shear_length_between_zero_and_gap',max(0.,float(np.max(-ls)),float(np.max(ls-heights))),0.)
    check('shear_smooth_wall_limit',float(np.max(abs(shear_length(heights,0,hc)))),0.)
    near=shear_length(np.array([hc*(1-1e-9),hc,hc*(1+1e-9)]),rms,hc)
    check('shear_branch_continuity',float(np.max(abs(near/hc-1))),2e-9)
    law=ShearFlowLaw(ShearOptions(lower_wall_velocity_m_s=.02,upper_wall_velocity_m_s=-.01),rms,hc)
    shifted=ShearFlowLaw(ShearOptions(lower_wall_velocity_m_s=.32,upper_wall_velocity_m_s=.29),rms,hc)
    check('shear_common_velocity_translation',float(np.max(abs((shifted.source(heights)-law.source(heights))/(.3*heights)-1))),1e-13)
    # Manufactured nonuniform moving-flow solution: integrate quadratic r*s exactly.
    length=1e-3;r0=1e14;s0=1e-8;beta=.5e-8;pin=2e4;pout=1e4
    exact=(pin-pout+r0*length*(1.5*s0+5*beta/6))/(1.5*r0*length)
    motion=[]
    for nodes in [21,41,81]:
        r=solve_fixed_contact_motion(lambda x:r0*(1+x/length),lambda x:s0+beta*x/length,
            0,length,pin,pout,nodes=nodes,tolerance=1e-7,max_nodes=20001,output_points=201)
        motion.append({'initial_nodes':nodes,'solver_nodes':r.solver_nodes,
            'relative_flux_error':abs(r.flux_per_width/exact-1),'face_flux_spread':r.face_flux_relative_spread})
    check('manufactured_nonuniform_motion_flux',max(r['relative_flux_error'] for r in motion),1e-6)
    check('manufactured_motion_mass_conservation',max(r['face_flux_spread'] for r in motion),1e-10)
    return {'checks':checks,'kernel':kernel,'root':root,'gaussian_asymptotic':asymptotic,'motion':motion}


def convergence_checks(energy,scale,bins):
    """Checks for resolved default K; upper bounds near percolation stay separate."""
    checks=[]
    def add(name,value,tolerance):
        checks.append({'name':name,'relative_difference':float(value),'tolerance':tolerance,
                       'passed':bool(value<=tolerance)})
    add('energy_n256_vs_adaptive',max(r['relative_difference'] for r in energy if r['n_q']==256),1e-7)
    for p in [1e3,1e4,1e5]:
        samples=[r for r in scale if r['method']=='conservative' and r['pressure_pa']==p]
        current=next(r for r in samples if r['pdf_log_step']==.01)
        fine=next(r for r in samples if r['pdf_log_step']==.0025)
        add(f'conservative_scale_K_at_{p:g}_Pa',abs(current['K_m3']/fine['K_m3']-1),.002)
        samples=[r for r in bins if r['pressure_pa']==p]
        current=next(r for r in samples if r['gap_points']==2049)
        fine=next(r for r in samples if r['gap_points']==8193)
        add(f'gap_bin_K_at_{p:g}_Pa',abs(current['K_m3']/fine['K_m3']-1),.002)
    return checks


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--psd',type=Path,default=PROJECT/'examples/data/roughness_psd.dat')
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--no-plots',action='store_true')
    args=parser.parse_args();out=args.output
    out.mkdir(parents=True,exist_ok=False)
    started=time.monotonic();psd=read_psd(args.psd)
    model=PerssonModel(psd,4e6,.4)
    pc=model.pressure_for_relative_contact_area(.42);hc=model.mean_gap(pc)
    pressures=[1e3,1e4,1e5,.98*pc,1.01*pc]
    physical=physical_limits(psd)
    (out/'physical_limits.json').write_text(json.dumps(physical,indent=2)+'\n')
    reference={p:energy_reference(psd,p) for p in [140.,1e4,1e5]}
    energy=[]
    for nq in [32,64,128,256,512]:
        m=PerssonModel(psd,4e6,.4,n_q=nq)
        for p,(h,error) in reference.items():
            value=m.mean_gap(p)
            energy.append({'n_q':nq,'pressure_pa':p,'gap_m':value,'reference_m':h,
                'relative_difference':abs(value/h-1),'quad_estimated_absolute_error_m':error})
    write_csv(out/'energy_convergence.csv',energy)
    scale=[]
    options=FlowOptions(quadrature_relative_tolerance=1e-4,max_gap_refinements=4)
    for method in ['fortran','conservative']:
        for step in [.04,.02,.01,.005,.0025]:
            m=PerssonModel(psd,4e6,.4,pdf_scale_method=method,pdf_log_step=step)
            for p in pressures:
                r=pressure_flow(m,p,NewtonianNoSlip(.1),options)
                diag=m.gap_distribution(p,[0.,40*psd.hrms]).diagnostics
                scale.append({'method':method,'pdf_log_step':step,'pressure_pa':p,
                    'mean_gap_m':r.mean_gap_m,'K_m3':r.geometric_conductivity_m3,
                    'upper_K_m3':r.upper_geometric_conductivity_m3,'phi_p':r.pressure_flow_factor,
                    'phi_s':float(shear_length(r.mean_gap_m,psd.hrms,hc)/psd.hrms),
                    'unrepresented_probability':r.unrepresented_noncontact_probability,
                    'mass_error':diag.quadrature_mass_error,'first_moment_relative_difference':diag.first_moment_relative_error,
                    'gap_bin_relative_difference':r.quadrature_relative_difference,'gap_points':r.gap_points,'status':r.status})
                write_csv(out/'scale_convergence.csv',scale)
                print(method,step,p,r.status,flush=True)
    bins=[]
    m=PerssonModel(psd,4e6,.4,pdf_scale_method='conservative')
    for p in [1e3,1e4,1e5]:
        for count in [129,257,513,1025,2049,4097,8193]:
            r=pressure_flow(m,p,NewtonianNoSlip(.1),FlowOptions(gap_points=count,
                quadrature_relative_tolerance=.99,max_gap_refinements=0))
            bins.append({'pressure_pa':p,'gap_points':count,'K_m3':r.geometric_conductivity_m3})
    write_csv(out/'gap_bin_convergence.csv',bins)
    cutoff=[]
    for p in [1e5,.98*pc]:
        for minimum in [1e-11,1e-12,1e-13]:
            m=PerssonModel(psd,4e6,.4,pdf_scale_method='conservative',pdf_min_u1_m=minimum)
            r=pressure_flow(m,p,NewtonianNoSlip(.1),options)
            cutoff.append({'pressure_pa':p,'pdf_min_u1_m':minimum,'K_m3':r.geometric_conductivity_m3,
                'upper_K_m3':r.upper_geometric_conductivity_m3,
                'unrepresented_probability':r.unrepresented_noncontact_probability,'status':r.status})
    write_csv(out/'cutoff_sensitivity.csv',cutoff)
    change=[]
    for method in ['fortran','conservative']:
        for p in pressures:
            samples=[r for r in scale if r['method']==method and r['pressure_pa']==p]
            fine=samples[-1]
            for r in samples:
                change.append({'method':method,'pressure_pa':p,'step':r['pdf_log_step'],
                    'K_change_to_finest':r['K_m3']/fine['K_m3']-1 if fine['K_m3']>0 else None,
                    'upper_K_change_to_finest':r['upper_K_m3']/fine['upper_K_m3']-1 if fine['upper_K_m3']>0 else None})
    sha=hashlib.sha256(args.psd.read_bytes()).hexdigest()
    numerical=convergence_checks(energy,scale,bins)
    summary={'schema_version':1,'base_commit':subprocess.check_output(['git','rev-parse','HEAD'],cwd=PROJECT,text=True).strip(),
        'source_hashes':{str(p.relative_to(PROJECT)):hashlib.sha256(p.read_bytes()).hexdigest()
             for p in sorted((PROJECT/'src/persson_contact').glob('*.py'))},
        'psd_sha256':sha,'effective_modulus_pa':4e6,'gamma':.4,'viscosity_pa_s':.1,
        'critical_relative_contact_area':.42,'percolation_pressure_pa':pc,
        'percolation_gap_m':hc,'rms_height_m':psd.hrms,
        'physical_limits_passed':all(c['passed'] for c in physical['checks']),
        'well_conditioned_K_convergence_passed':all(c['passed'] for c in numerical),
        'convergence_checks':numerical,
        'physical_checks':physical['checks'],'scale_changes':change,
        'energy_n256_max_relative_difference':max(r['relative_difference'] for r in energy if r['n_q']==256),
        'elapsed_seconds':time.monotonic()-started,
        'limitations':['Finest-grid differences are empirical convergence evidence, not rigorous error bounds.',
            'Scale cutoffs, missing-phase bounds, closure physics and reflected-PDF moment mismatch remain model uncertainty.',
            'Near percolation upper-bound relative changes require separate attention even when default K is zero.',
            'No Fortran matching, gap rescaling, probability renormalization or large-gap splice is imposed.']}
    (out/'summary.json').write_text(json.dumps(summary,indent=2,allow_nan=False)+'\n')
    if not args.no_plots:plot(out,energy,scale,bins)
    print(json.dumps({'physical_limits_passed':summary['physical_limits_passed'],
        'elapsed_seconds':summary['elapsed_seconds']}),flush=True)
    return 0 if summary['physical_limits_passed'] and summary['well_conditioned_K_convergence_passed'] else 2


def plot(out,energy,scale,bins):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig,axs=plt.subplots(2,2,figsize=(11,7.5),constrained_layout=True)
    for p in [140.,1e4,1e5]:
        r=[x for x in energy if x['pressure_pa']==p]
        axs[0,0].loglog([x['n_q'] for x in r],[max(x['relative_difference'],1e-16) for x in r],'o-',label=f'{p:g} Pa')
    axs[0,0].set(xlabel='Energy n_q',ylabel='Difference from adaptive integral')
    for method,style in [('fortran','--'),('conservative','-')]:
        for p,color in [(1e4,'C1'),(1e5,'C2')]:
            r=[x for x in scale if x['method']==method and x['pressure_pa']==p]
            axs[0,1].loglog([x['pdf_log_step'] for x in r[:-1]],
                [abs(x['K_m3']/r[-1]['K_m3']-1) for x in r[:-1]],'o'+style,color=color,label=f'{method}, {p:g} Pa')
    axs[0,1].set(xlabel='PDF log-scale step',ylabel='K change from step 0.0025')
    for p in [1e3,1e4,1e5]:
        r=[x for x in bins if x['pressure_pa']==p]
        axs[1,0].loglog([x['gap_points'] for x in r[:-1]],[abs(x['K_m3']/r[-1]['K_m3']-1) for x in r[:-1]],'o-',label=f'{p:g} Pa')
    axs[1,0].set(xlabel='Gap-bin nodes (conservative scales)',ylabel='K change from 8193 nodes')
    for method,style in [('fortran','--'),('conservative','-')]:
        r=[x for x in scale if x['method']==method and x['pressure_pa']==1e5]
        axs[1,1].semilogx([x['pdf_log_step'] for x in r],[x['phi_p'] for x in r],'o'+style,label=method)
    axs[1,1].set(xlabel='PDF log-scale step',ylabel='Pressure factor at 0.1 MPa contact')
    for ax in axs.flat:ax.grid(alpha=.25);ax.legend(fontsize=8)
    fig.suptitle('Independent numerical convergence | same PSD, E*=4 MPa, gamma=0.4')
    fig.savefig(out/'convergence.png',dpi=160);plt.close(fig)


if __name__=='__main__':raise SystemExit(main())
