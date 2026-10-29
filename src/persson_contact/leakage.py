"""Leakage application: geometry selection, table assembly and output.

Shared conservative pressure/flux equations are in steady.py, reference mechanics in mechanics.py,
profile input in profiles.py and local conductance in module 4.
"""
import csv
from dataclasses import dataclass
import json
from pathlib import Path
import numpy as np

from .interface import positive
from .boundaries import read_pressure_boundaries
from .mechanics import make_mechanics
from .steady import solve_steady_flow
from .transport_table import TransportTable


@dataclass
class LeakagePlan:
    mechanics: object
    width_m: float
    inlet_pressures_pa: np.ndarray
    outlet_pressure_pa: float
    coupling: str
    method: str
    percolation_policy: str
    spc: float
    critical_alpha: float
    interpolation: str
    solver: dict
    compute_upper_bound: bool
    original: dict


def prepare_leakage(raw,base_directory,effective_modulus_pa):
    if not isinstance(raw,dict):raise ValueError('leakage must be a dictionary')
    allowed={'mode','geometry','geometry_parameters','seal_width_m','inlet_pressure_pa','inlet_pressures_pa','outlet_pressure_pa',
             'coupling','conductance_model','percolation_policy','spc','critical_alpha','database_interpolation',
             'solver','compute_upper_bound','transport'}
    if set(raw)-allowed:raise ValueError(f'Unknown leakage parameters: {sorted(set(raw)-allowed)}')
    if raw.get('mode')!='spatial':raise ValueError('Spatial leakage requires mode=spatial')
    if raw.get('transport','liquid')!='liquid':
        raise NotImplementedError('Gas/ballistic leakage (legacy 5000) requires separate validation; no liquid fallback')
    geometry=raw.get('geometry',{})
    if isinstance(geometry,str):
        presets=raw.get('geometry_parameters',{})
        if not isinstance(presets,dict) or geometry not in presets or not isinstance(presets[geometry],dict):
            raise ValueError('Selected geometry needs a matching geometry_parameters dictionary')
        if 'type' in presets[geometry]:raise ValueError('Geometry preset type is supplied by the geometry selector')
        geometry={'type':geometry,**presets[geometry]}
    mech=make_mechanics(geometry,base_directory,effective_modulus_pa)
    boundary=read_pressure_boundaries(raw);inlet=boundary.inlet_pressures_pa;out=boundary.outlet_pressure_pa
    coupling=raw.get('coupling','local_unloading');method=raw.get('conductance_model','effective_medium')
    if coupling not in ('fixed_contact','local_unloading'):raise ValueError('coupling must be fixed_contact or local_unloading')
    if method not in ('effective_medium','critical_junction'):raise ValueError('conductance_model must be effective_medium or critical_junction')
    policy=raw.get('percolation_policy','legacy_spc_plateau')
    if policy not in ('physical_zero','legacy_spc_plateau'):raise ValueError('Invalid percolation_policy')
    spc=positive(raw.get('spc',.8),'spc')
    if spc>=1:raise ValueError('spc must be below 1')
    interpolation=raw.get('database_interpolation','log_pchip')
    if interpolation not in ('log_pchip','linear'):raise ValueError('Invalid database_interpolation')
    if not isinstance(raw.get('solver',{}),dict):raise ValueError('leakage.solver must be a dictionary')
    config=dict(raw.get('solver',{}));defaults={'nodes':201,'output_points':1001,'tolerance':1e-5,'max_nodes':20001}
    if set(config)-set(defaults):raise ValueError('Unknown leakage solver options')
    defaults.update(config)
    for key in ('nodes','output_points','max_nodes'):
        if type(defaults[key]) is not int or defaults[key]<5:raise ValueError(f'{key} must be an integer >=5')
    if defaults['max_nodes']<defaults['nodes'] or defaults['max_nodes']>100001 or defaults['output_points']>100001:
        raise ValueError('Invalid maximum leakage mesh size')
    if not 0<positive(defaults['tolerance'],'tolerance')<.01:raise ValueError('Invalid leakage solver tolerance')
    upper=raw.get('compute_upper_bound',True)
    if type(upper) is not bool:raise ValueError('compute_upper_bound must be True or False')
    return LeakagePlan(mech,positive(raw.get('seal_width_m'),'seal_width_m'),inlet,out,coupling,method,
        policy,spc,positive(raw.get('critical_alpha',1.),'critical_alpha'),interpolation,defaults,upper,raw)


def required_pressures(plan,model,configured,percolation):
    """Retain requested nodes and add needed upper range and plateau node."""
    needed=plan.mechanics.maximum_reference_pressure(plan.coupling!='fixed_contact')
    if plan.coupling!='fixed_contact':needed-=min(min(plan.inlet_pressures_pa),plan.outlet_pressure_pa)
    perc=model.pressure_for_relative_contact_area(percolation)
    plateau=plan.spc*perc if plan.percolation_policy=='legacy_spc_plateau' else None
    # A plateau only needs database coverage up to itself if it lies within the
    # load range; otherwise no material parameter should inflate the pressure scan.
    active_plateau=plateau if plateau is not None and plateau<=max(needed,max(configured)) else None
    maximum=max(needed*1.01,max(configured))
    p=np.asarray(configured,dtype=float)
    if maximum>p[-1]:
        extra=max(3,int(np.ceil(np.log(maximum/p[-1])/np.median(np.diff(np.log(p))))))
        p=np.r_[p,np.geomspace(p[-1],maximum,extra+1)[1:]]
    if active_plateau is not None:
        if active_plateau<p[0]:raise ValueError('Database minimum exceeds spc plateau pressure')
        p=np.unique(np.r_[p,active_plateau])
    return p,active_plateau,perc


def make_transport(plan,model,interface,pressure_rows,percolation,area_tolerance=.005):
    p=np.array([r['pressure_pa'] for r in pressure_rows]);h=np.array([r['mean_gap_m'] for r in pressure_rows])
    perc=model.pressure_for_relative_contact_area(percolation)
    plateau=plan.spc*perc if plan.percolation_policy=='legacy_spc_plateau' else None
    # If the requested scan lies below the cap, it has no effect and is omitted.
    policy=plan.percolation_policy
    if plateau is not None and plateau>max(p):policy='physical_zero';plateau=None
    if plan.method=='effective_medium':
        k=np.array([r['geometric_conductivity_m3'] for r in pressure_rows])
        upper=np.array([np.inf if r['upper_bound_unbounded'] else r['upper_geometric_conductivity_m3'] for r in pressure_rows])
    else:
        values=[model.critical_separation(float(pi),percolation,area_tolerance) for pi in p]
        invalid=[v.status for v in values if v.status not in ('found','threshold_not_reached')]
        if invalid:raise ValueError(f'Critical-junction database has unresolved u1: {invalid[0]}')
        k=np.array([(plan.critical_alpha*v.critical_separation_m)**3 if v.status=='found' else 0 for v in values]);upper=k.copy()
    return TransportTable(p,h,k,upper,model,interface,interpolation=plan.interpolation,
        percolation_policy=policy,plateau_pressure_pa=plateau)


def _csv(path,rows):
    with Path(path).open('w',encoding='utf-8',newline='') as stream:
        writer=csv.DictWriter(stream,fieldnames=list(rows[0]),lineterminator='\n');writer.writeheader();writer.writerows(rows)


def solve_case(plan,table,inlet,*,upper=False):
    mech=plan.mechanics;out=plan.outlet_pressure_pa
    def contact(x,pf):
        return mech.contact_pressure(x,pf,plan.coupling)
    def resistance(x,pf):return table.resistance(contact(x,pf),upper)
    # A barrier blocked even at maximum boundary pressure guarantees no flow.
    xcheck=np.unique(np.r_[mech.integration_points,np.linspace(mech.x_m[0],mech.x_m[-1],4097)])
    best=contact(xcheck,np.full_like(xcheck,max(inlet,out)))
    barrier=np.isinf(table.resistance(best,upper))
    if inlet!=out and np.any(barrier):
        x=np.linspace(mech.x_m[0],mech.x_m[-1],plan.solver['output_points']);pf=np.full_like(x,np.nan)
        first,last=xcheck[barrier][0],xcheck[barrier][-1]
        pf[x<first]=inlet;pf[x>last]=out;pf[0]=inlet;pf[-1]=out
        return {'x':x,'fluid':pf,'contact':contact(x,pf),'flux':0.,'status':'blocked_pressure_not_unique',
                'solver_nodes':0,'equation_residual':0.,'pressure_drop_error':0.,'flux_spread':0.,'iterations':0,'flux_check_pressure_fraction':None,'grid_difference':0.}
    def checked_resistance(x,pf):
        r=resistance(x,pf)
        if np.any(np.isinf(r)):
            raise ValueError('Physical-zero mobility encountered inside coupled region; no uniquely conducting branch resolved. '
                             'Use a structural/opening model or explicitly choose the legacy_spc_plateau approximation.')
        return r
    solution=solve_steady_flow(checked_resistance,mech.x_m[0],mech.x_m[-1],float(inlet),out,required_points=mech.integration_points,**plan.solver)
    pc=contact(solution.x_m,solution.fluid_pressure_pa)
    return {'x':solution.x_m,'fluid':solution.fluid_pressure_pa,'contact':pc,'flux':solution.flux_per_width_m2_s,
            'status':'zero_pressure_drop' if inlet==out else 'converged',
            'solver_nodes':solution.solver_nodes,'equation_residual':solution.maximum_equation_residual,
            'pressure_drop_error':solution.pressure_drop_relative_error,'flux_spread':solution.face_flux_relative_spread,
            'iterations':solution.iterations,'flux_check_pressure_fraction':solution.flux_check_pressure_fraction,'grid_difference':solution.grid_relative_flux_difference}


def run_spatial_leakage(plan,table,output,*,plot=False):
    output=Path(output);rows=[];profiles=[]
    database=[]
    used=table.geometric_conductivity(table.p);used_upper=table.geometric_conductivity(table.p,True)
    for j,p in enumerate(table.p):
        database.append({'contact_pressure_pa':float(p),'mean_gap_m':float(table.h[j]),
            'relative_contact_area':float(table.model.area(p)),
            'raw_geometric_conductivity_m3':float(table.k[j]),
            'used_geometric_conductivity_m3':float(used[j]),
            'used_upper_geometric_conductivity_m3':float(used_upper[j]) if np.isfinite(used_upper[j]) else None,
            'conductance_model':plan.method,'percolation_policy':plan.percolation_policy})
    _csv(output/'leakage_transport_table.csv',database)
    for i,inlet in enumerate(plan.inlet_pressures_pa):
        result=solve_case(plan,table,float(inlet))
        upper=None;upper_error=None
        if plan.compute_upper_bound:
            try:upper=solve_case(plan,table,float(inlet),upper=True)
            except ValueError as exc:upper_error=str(exc)
        x,pf,pc=result['x'],result['fluid'],result['contact'];c0,f0=plan.mechanics.evaluate(x)
        h=table.mean_gap(pc[np.isfinite(pc)])
        gaps=np.full_like(pc,np.nan);gaps[np.isfinite(pc)]=h
        area=np.full_like(pc,np.nan);mask=np.isfinite(pc);area[mask]=table.model.area(pc[mask])
        low=(pc>0)&(pc<table.p[0]);opened=(pc<=0)&np.isfinite(pc)
        k=np.full_like(pc,np.nan);k[mask]=table.geometric_conductivity(pc[mask])
        data=[]
        for j in range(len(x)):
            finite=lambda value:float(value) if np.isfinite(value) else None
            data.append({'x_m':float(x[j]),'reference_contact_pressure_pa':float(c0[j]),
                'reference_interface_pressure_pa':float(f0[j]),'fluid_pressure_pa':finite(pf[j]),
                'contact_pressure_pa':finite(pc[j]),'mean_gap_m':finite(gaps[j]),
                'relative_contact_area':finite(area[j]),'geometric_conductivity_m3':finite(k[j]),
                'open_reservoir_limit':bool(opened[j]),'low_pressure_continuation':bool(low[j])})
        file=f'leakage_profile_{i:03d}.csv';_csv(output/file,data)
        q=result['flux']*plan.width_m;qupper=upper['flux']*plan.width_m if upper is not None else None
        row={'case_index':i,'inlet_pressure_pa':float(inlet),'outlet_pressure_pa':plan.outlet_pressure_pa,
             'flow_rate_m3_s':q,'pdf_bound_min_flow_m3_s':None if qupper is None else min(q,qupper),
             'pdf_bound_max_flow_m3_s':None if qupper is None else max(q,qupper),
             'status':('converged_open_limit' if result['status']=='converged' and np.any(opened) else result['status']),'geometry':plan.mechanics.metadata['type'],
             'coupling':plan.coupling,'conductance_model':plan.method,'percolation_policy':plan.percolation_policy,
             'seal_length_m':plan.mechanics.length_m,'seal_width_m':plan.width_m,
             'solver_nodes':result['solver_nodes'],'maximum_equation_residual':result['equation_residual'],
             'pressure_drop_relative_error':result['pressure_drop_error'],'face_flux_relative_spread':result['flux_spread'],
             'flux_check_pressure_fraction':result['flux_check_pressure_fraction'],
             'grid_relative_flux_difference':result['grid_difference'],
             'maximum_open_excess_pressure_pa':float(np.nanmax(np.maximum(pf-c0-f0,0))),
             'open_limit_points':int(sum(opened)),'low_pressure_continuation_points':int(sum(low)),
             'profile_file':file,'upper_bound_status':upper['status'] if upper is not None else None,
             'upper_bound_error':upper_error}
        rows.append(row);profiles.append(data)
        print(f"泄漏工况 {i+1}/{len(plan.inlet_pressures_pa)}: pin={inlet:g} Pa, Q={q:.6g} m^3/s, {row['status']}",flush=True)
    _csv(output/'leakage_cases.csv',rows)
    summary={'scope':'Incompressible Newtonian liquid; one-dimensional local reference-traction unloading. '
        'No nonlocal FEM compliance, cavitation or ballistic gas.',
        'geometry':plan.mechanics.metadata,'reference_external_load_policy':'Held fixed for all interface-pressure cases; '
        'reference interface pressure is included only when explicitly supplied. Exterior FEM load is not added again.',
        'transport_database':table.metadata,'transport_file':'leakage_transport_table.csv','solver':plan.solver,'parameters':plan.original,
        'open_region_policy':'Reservoir (zero-resistance) limit for pc<=0; no finite gap reported. '
            'A full open-gap deformation model is required for physical aperture and fully open seals.',
        'pdf_bounds_scope':('No PDF-based bound for critical-junction model; equal bounds do not quantify its model error. ' if plan.method=='critical_junction' else '')+'Same local approximation with missing-PDF conductivity bounds; not statistical confidence intervals.',
        'cases':rows}
    (output/'leakage_diagnostics.json').write_text(json.dumps(summary,indent=2,allow_nan=False)+'\n')
    if plot:_plot(output,rows,profiles)
    return summary


def _plot(output,rows,profiles):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig,axes=plt.subplots(1,3,figsize=(15,4),layout='constrained')
    order=np.argsort([r['inlet_pressure_pa'] for r in rows])
    pin=np.array([rows[i]['inlet_pressure_pa']/1e6 for i in order])
    axes[0].plot(pin,[rows[i]['flow_rate_m3_s'] for i in order],'o-',label='Baseline')
    lower=[np.nan if rows[i]['pdf_bound_min_flow_m3_s'] is None else rows[i]['pdf_bound_min_flow_m3_s'] for i in order]
    upper=[np.nan if rows[i]['pdf_bound_max_flow_m3_s'] is None else rows[i]['pdf_bound_max_flow_m3_s'] for i in order]
    axes[0].fill_between(pin,lower,upper,alpha=.2,label='Missing-PDF bounds')
    axes[0].legend(fontsize=8)
    axes[0].set_title('Leakage pressure scan')
    axes[0].set(xlabel='Inlet pressure (MPa)',ylabel='Signed leakage (m^3/s)')
    chosen=len(rows)-1;data=profiles[chosen];x=[v['x_m']*1e3 for v in data]
    for key,label in (('reference_contact_pressure_pa','FEM/reference contact'),('contact_pressure_pa','Contact'),('fluid_pressure_pa','Fluid')):
        axes[1].plot(x,[np.nan if v[key] is None else v[key]/1e6 for v in data],label=label)
    axes[1].set(xlabel='Position (mm)',ylabel='Pressure (MPa)',title=f"State at inlet {rows[chosen]['inlet_pressure_pa']/1e6:g} MPa");axes[1].legend(fontsize=8)
    axes[2].plot(x,[np.nan if v['mean_gap_m'] is None else v['mean_gap_m']*1e6 for v in data])
    axes[2].set(xlabel='Position (mm)',ylabel='Resolved mean gap (um)',title='Open-limit gaps are undefined')
    for ax in axes:ax.grid(alpha=.25)
    fig.savefig(output/'leakage.png',dpi=180);plt.close(fig)
