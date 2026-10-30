"""Configuration, coupled velocity cases, geometry and friction outputs."""
import hashlib
import json
from pathlib import Path

import numpy as np

from .config import Options
from .contact_friction import finite,evaluate_contact_friction
from .contact_friction_output import prepare_contact_friction
from .core import PerssonModel
from .dry_contact import RoughContactLaw,solve_rough_contact
from .flow import FlowOptions
from .fluid_friction import wall_tractions
from .interface import make_interface,positive
from .macro_output import prepare_macro,_keys,_csv
from .mixed import MixedSolution,solve_mixed_continuation
from .mixed_transport import MixedTransport
from .psd import read_psd
from .table import build_contact_table


def prepare_mixed(parameters,base):
    raw=parameters.get('mixed_friction',{})
    _keys(raw,{'lower_wall_velocities_m_s','rough_surface','inlet_pressure_pa','outlet_pressure_pa',
        'cavitation_model','cavitation_friction','transport_points','relative_tolerance','max_iterations'},'mixed_friction')
    if raw.get('cavitation_model','mass_conserving')!='mass_conserving':
        raise NotImplementedError('Choose mass_conserving cavitation')
    if raw.get('cavitation_friction','uniform_wetted_fraction')!='uniform_wetted_fraction':
        raise NotImplementedError('Only the explicit uniform_wetted_fraction cavitation friction approximation is implemented')
    velocities=raw.get('lower_wall_velocities_m_s')
    if not isinstance(velocities,(list,tuple)) or not velocities:raise ValueError('Provide lower_wall_velocities_m_s')
    speeds=tuple(finite(v,'lower_wall_velocity_m_s') for v in velocities)
    pin=finite(raw.get('inlet_pressure_pa',0.),'inlet_pressure_pa',minimum=0.)
    pout=finite(raw.get('outlet_pressure_pa',0.),'outlet_pressure_pa',minimum=0.)
    tol=positive(raw.get('relative_tolerance',1e-8),'mixed relative_tolerance')
    limit=raw.get('max_iterations',100);points=raw.get('transport_points',161)
    if tol>=.001 or type(limit) is not int or limit<1:raise ValueError('Invalid mixed solver settings')
    if type(points) is not int or not 41<=points<=2001:raise ValueError('transport_points must be 41..2001')
    rough=raw.get('rough_surface','upper')
    if rough not in ('upper','lower'):raise ValueError('rough_surface must be upper or lower')
    interface=make_interface(parameters.get('interface',{}))
    if interface.has_slip:raise NotImplementedError('Mixed friction requires no-slip Newtonian liquid')
    macro,solver,hertz,elastic,representation,kind=prepare_macro(parameters,base)
    if kind!='persson' or representation!='upper_elastic_lower_rigid':
        raise NotImplementedError('Mixed workflow currently requires Persson contact and upper_elastic_lower_rigid cylinder/plane geometry')
    friction=prepare_contact_friction(parameters.get('contact_friction',{}),base)
    opts=Options(**parameters.get('options',{}))
    if opts.pdf_scale_method!='conservative':raise ValueError('Mixed friction requires pdf_scale_method=conservative')
    flow_options=FlowOptions(**parameters.get('pressure_flow',{}))
    fluid_friction=parameters.get('fluid_friction',{})
    _keys(fluid_friction,{'model','minimum_shear_gap_m'},'fluid_friction')
    if fluid_friction.get('model','persson_scale_u1')!='persson_scale_u1':
        raise NotImplementedError('Choose fluid_friction.model=persson_scale_u1')
    positive(fluid_friction.get('minimum_shear_gap_m'),'minimum_shear_gap_m')
    return raw,speeds,pin,pout,tol,limit,points,rough,interface,macro,solver,hertz,elastic,friction,opts,flow_options


def make_contact_database(parameters,base,macro,solver,hertz,elastic,opts):
    psd_raw=parameters.get('psd',{})
    _keys(psd_raw,{'path','scale_c','tail_is_extrapolated'},'psd')
    path=Path(psd_raw['path']);path=path if path.is_absolute() else base/path
    psd=read_psd(path,psd_raw.get('scale_c',1.))
    model=PerssonModel(psd,hertz.effective_modulus_pa,opts.gamma,n_q=opts.n_q,
        pdf_afact=opts.pdf_afact,pdf_log_step=opts.pdf_log_step,
        pdf_max_relative_contact_area=opts.pdf_max_relative_contact_area,
        pdf_min_relative_contact_area=opts.pdf_min_relative_contact_area,
        pdf_min_u1_m=opts.pdf_min_u1_m,pdf_scale_method=opts.pdf_scale_method)
    settings={'minimum_pressure_pa':1e-100,'maximum_pressure_pa':2e6,'points':1201,**macro.get('contact_table',{})}
    _keys(settings,{'minimum_pressure_pa','maximum_pressure_pa','points'},'macro contact_table')
    low=positive(settings['minimum_pressure_pa'],'minimum_pressure_pa');high=positive(settings['maximum_pressure_pa'],'maximum_pressure_pa');n=settings['points']
    if type(n) is not int or not 101<=n<=20001 or not low<high or high<=hertz.peak_pressure_pa:
        raise ValueError('Invalid mixed contact-table range or points')
    if low*np.ptp(elastic.edges_m)/hertz.line_load_n_m>solver['relative_tolerance']:
        raise ValueError('Lower contact-table minimum pressure to control omitted far-tail load')
    focus=max(low,min(high,hertz.peak_pressure_pa*1e-3))
    pressures=np.unique(np.r_[np.geomspace(low,focus,n//3),np.geomspace(focus,high,n-n//3+1)])
    return model,build_contact_table(model,pressures),hashlib.sha256(path.read_bytes()).hexdigest()


def run_mixed_project(parameters,*,base_directory,save_plots=True):
    base=Path(base_directory).resolve()
    (raw,speeds,pin,pout,tol,limit,points,rough,interface,macro,solver,hertz,
        elastic,friction,opts,flow_options)=prepare_mixed(parameters,base)
    name=parameters.get('output_directory')
    if not isinstance(name,str) or not name.strip():raise ValueError('Provide output_directory')
    out=Path(name);out=out if out.is_absolute() else base/out
    if out.exists() and (not out.is_dir() or any(out.iterdir())):raise FileExistsError(f'输出目录非空，请使用新的 output_directory: {out}')
    model,table,psd_sha=make_contact_database(parameters,base,macro,solver,hertz,elastic,opts)
    law=RoughContactLaw(table)
    dry=solve_rough_contact(elastic,hertz,law,tolerance=solver['relative_tolerance'],max_iterations=solver['max_iterations'])
    transport=MixedTransport(model,law,interface,flow_options,points=points,
        percolation_relative_contact_area=opts.critical_relative_contact_area,rough_surface=rough,
        minimum_shear_gap_m=parameters['fluid_friction']['minimum_shear_gap_m'])
    cases=[];profiles=[];solutions=[];factors_rows=[]
    dx=elastic.dx_m
    for index,speed in enumerate(speeds):
        continuation=[]
        if speed==0 and pin==0 and pout==0:
            # Exact dry limit; no artificial tiny speed or transport floor.
            v=MixedSolution(dry.mean_gap_m,dry.pressure_pa,np.zeros_like(dry.pressure_pa),
                np.ones_like(dry.pressure_pa),dry.displacement_m,dry.rigid_offset_m,
                np.zeros(len(dry.pressure_pa)+1),0,dry.equilibrium_relative_residual,
                dry.equilibrium_relative_residual,dry.load_relative_error,0.)
        else:
            v,continuation=solve_mixed_continuation(elastic,hertz,law,transport.at_velocity,speed,dry,
                inlet_pressure_pa=pin,outlet_pressure_pa=pout,tolerance=tol,max_iterations=limit)
        total=v.contact_pressure_pa+v.fluid_pressure_pa
        edge=float(max(total[0],total[-1])/max(total))
        # Positive imposed edge pressure need not decay to zero; check only zero-reservoir cases.
        if pin==pout==0 and edge>solver['edge_pressure_ratio_tolerance']:
            raise ValueError('Mixed pressure remains significant at edges; enlarge macro domain')
        a=model.area(v.contact_pressure_pa)
        solid=evaluate_contact_friction(elastic.edges_m,a,v.contact_pressure_pa,v.mean_gap_m,-speed,
            friction.law,temperature_k=friction.temperature_k,transverse_width_m=friction.transverse_width_m)
        ff=transport.factors(v.mean_gap_m)
        pf=v.fluid_pressure_pa
        gradients=np.r_[(pf[0]-pin)/(dx/2),np.diff(pf)/dx,(pout-pf[-1])/(dx/2)]
        grad=(gradients[:-1]+gradients[1:])/2
        tr=wall_tractions(v.mean_gap_m,ff,interface.viscosity_pa_s,grad,speed,0.,
            rough_surface=rough,saturation=v.saturation)
        lower_fluid=float(np.sum(tr['lower_total_pa'])*dx)
        upper_fluid=float(np.sum(tr['upper_total_pa'])*dx)
        lower_solid=solid.summary['lower_friction_force_per_width_n_m']
        upper_solid=solid.summary['upper_friction_force_per_width_n_m']
        slope=np.gradient(v.mean_gap_m,elastic.x_m,edge_order=2)
        upper_macro_pressure=float(-np.sum(pf*slope)*dx)
        wf=float(np.sum(pf)*dx);wc=float(np.sum(v.contact_pressure_pa)*dx)
        q=float(np.mean(v.face_flux_m2_s))
        signed_lower=lower_solid+lower_fluid
        width=friction.transverse_width_m
        case={'case_index':index,'lower_wall_velocity_m_s':speed,'upper_wall_velocity_m_s':0.,
            'relative_velocity_m_s':-speed,'contact_line_load_n_m':wc,'fluid_line_load_n_m':wf,
            'fluid_load_fraction':wf/hertz.line_load_n_m,
            'real_contact_area_per_width_m':float(np.sum(a)*dx),
            'minimum_mean_gap_m':float(min(v.mean_gap_m)),'maximum_fluid_pressure_pa':float(max(pf)),
            'minimum_saturation':float(min(v.saturation)),'cavitated_cells':int(np.sum(v.saturation<1-1e-10)),
            'flux_per_width_m2_s':q,'flow_rate_m3_s':None if width is None else q*width,
            'lower_solid_force_n_m':lower_solid,'lower_fluid_force_n_m':lower_fluid,
            'lower_fluid_sliding_force_n_m':float(np.sum(tr['lower_sliding_pa'])*dx),
            'lower_fluid_pressure_flow_force_n_m':float(np.sum(tr['lower_pressure_flow_pa'])*dx),
            'lower_total_force_n_m':signed_lower,
            'lower_total_force_n':None if width is None else signed_lower*width,
            'lower_resisting_friction_coefficient':None if speed==0 else -np.sign(speed)*signed_lower/hertz.line_load_n_m,
            'upper_solid_force_n_m':upper_solid,'upper_fluid_wall_force_n_m':upper_fluid,
            'upper_macroscopic_fluid_pressure_force_n_m':upper_macro_pressure,
            'load_relative_error':v.load_relative_error,'geometry_relative_residual':v.geometry_relative_residual,
            'mass_relative_residual':v.mass_relative_residual,'solver_residual':v.residual,
            'complementarity_max_pa':float(np.max(pf*(1-v.saturation))),
            'edge_pressure_ratio':edge,'iterations':v.iterations,'continuation_states':len(continuation),
            'profile_file':f'mixed_profile_{index:03d}.csv'}
        cases.append(case);solutions.append(v)
        rows=[]
        for i,x in enumerate(elastic.x_m):
            rows.append({'x_m':x,'cell_left_m':elastic.edges_m[i],'cell_right_m':elastic.edges_m[i+1],
                'initial_upper_surface_m':hertz.profile(elastic.x_m)[i],
                'rigid_positioned_upper_surface_m':hertz.profile(elastic.x_m)[i]+v.rigid_offset_m,
                'relative_elastic_displacement_m':v.displacement_m[i],
                'deformed_upper_surface_m':v.mean_gap_m[i],'lower_surface_m':0.,
                'mean_gap_m':v.mean_gap_m[i],'contact_pressure_pa':v.contact_pressure_pa[i],
                'fluid_pressure_pa':pf[i],'relative_contact_area':a[i],'liquid_saturation':v.saturation[i],
                'left_face_flux_m2_s':v.face_flux_m2_s[i],'right_face_flux_m2_s':v.face_flux_m2_s[i+1],
                'phi_f':ff.phi_f[i],'phi_fs':ff.phi_fs[i],'phi_fp':ff.phi_fp[i],
                'lower_solid_traction_pa':solid.profile['lower_contact_traction_pa'][i],
                **{k:values[i] for k,values in tr.items()}})
        profiles.append(rows)
        print(f'混合润滑 V_plane={speed:g} m/s：流体承载 {wf/hertz.line_load_n_m:.3%}，最小间隙 {min(v.mean_gap_m):.5g} m，空化单元 {case["cavitated_cells"]}',flush=True)
    outputs=['mixed_cases.csv','mixed_transport_table.csv','mixed_contact_table.csv','mixed_diagnostics.json',*[c['profile_file'] for c in cases]]
    summary={'workflow_status':'complete','scope':'Steady 1D isothermal constant-Newtonian no-slip cylinder/plane elastic mixed lubrication; connected liquid-mobility branch.',
        'reference_frame':'Cylinder-following; upper elastic cylinder velocity=0, lower rigid plane velocity specified; signed flux in this frame.',
        'pressure_convention':'Gauge Pa, cavitation pressure=0; fully flooded reservoirs at boundary inflow. No negative-pressure clipping.',
        'normal_equations':'h=h0+x^2/(2R)+C(pc+pf); pc=Persson(h); integral(pc+pf)dx=Wprime.',
        'cavitation':'Conservative finite-volume pressure/saturation complementarity p>=0, 0<=theta<=1, p(1-theta)=0; harmonic mobility and upwind theta.',
        'cavitation_friction':'Explicit approximation: sliding traction multiplied by liquid saturation; not a microscopically validated cavitation shear theory.',
        'reported_body':'Moving lower plane: lower_resisting_friction_coefficient includes solid and averaged fluid wall traction. Its macroscopic slope is zero.',
        'force_scope':'Sum of modelled solid contact and continuum fluid contributions; boundary-film shear below minimum_shear_gap_m is not included.',
        'upper_force_note':'Upper wall traction includes the roughness-scale form drag assigned by the interpolation. Macroscopic fluid pressure projection is exported separately; upper solid normal-pressure slope force is not included in a total upper drag.',
        'transport':transport.metadata,'solid_friction':friction.law.metadata,
        'transverse_width_m':friction.transverse_width_m,
        'target_line_load_n_m':hertz.line_load_n_m,'psd_sha256':psd_sha,
        'far_tail_fixed_gap_load_bound_n_m':float(table.pressures_pa[0]*np.ptp(elastic.edges_m)),
        'limitations':['No viscoelastic hysteresis, heating, non-Newtonian fluid, slip, transient motion or starvation.',
            'No disconnected trapped-pocket solution; physical-zero mobility is never replaced by an spc plateau.',
            'Scale-inverse-gap and large-gap interpolation are constitutive approximations; numerical convergence does not validate them.',
            'Local unknown-PDF conductivity bounds do not bound the nonlinear coupled pressure or friction.',
            'Only mean deformed geometry is reconstructed, not individual roughness phases.'],
        'cases':cases,'parameters_snapshot':parameters,'outputs':outputs}
    out.mkdir(parents=True,exist_ok=True)
    _csv(out/'mixed_cases.csv',cases);_csv(out/'mixed_transport_table.csv',transport.rows)
    _csv(out/'mixed_contact_table.csv',[{'pressure_pa':p,'mean_gap_m':h,'relative_contact_area':a}
        for p,h,a in zip(table.pressures_pa,table.mean_gaps_m,table.contact_fractions)])
    for case,rows in zip(cases,profiles):_csv(out/case['profile_file'],rows)
    if save_plots:
        plot_mixed(out,elastic,hertz,solutions,cases)
        outputs+=['mixed_profiles.png','mixed_friction.png','mixed_geometry.png']
    (out/'mixed_diagnostics.json').write_text(json.dumps(summary,ensure_ascii=False,indent=2,allow_nan=False)+'\n')
    return summary


def plot_mixed(out,elastic,hertz,solutions,cases):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    x=elastic.x_m*1e3
    fig,axes=plt.subplots(3,1,figsize=(9,9),sharex=True,layout='constrained')
    for v,c in zip(solutions,cases):
        label=f'U_plane={c["lower_wall_velocity_m_s"]:g} m/s'
        line=axes[0].plot(x,v.fluid_pressure_pa/1e6,label=label)[0]
        axes[0].plot(x,v.contact_pressure_pa/1e6,'--',color=line.get_color())
        axes[1].plot(x,v.mean_gap_m*1e6,label=label)
        axes[2].plot(x,v.saturation,label=label)
    axes[0].set(ylabel='Pressure [MPa]',title='Solid dashed; fluid solid lines')
    axes[1].set(ylabel='Mean gap / deformed upper [um]');axes[1].axhline(0,color='k',lw=.8)
    axes[2].set(xlabel='x [mm]',ylabel='Liquid saturation',ylim=(-.02,1.02))
    for ax in axes:ax.grid(alpha=.2);ax.legend(fontsize=8)
    fig.savefig(out/'mixed_profiles.png',dpi=170);plt.close(fig)
    fig,axes=plt.subplots(1,2,figsize=(10,4),layout='constrained')
    moving=sorted([c for c in cases if c['lower_wall_velocity_m_s']>0],key=lambda c:c['lower_wall_velocity_m_s'])
    if moving:
        speeds=[c['lower_wall_velocity_m_s'] for c in moving]
        axes[0].semilogx(speeds,[c['lower_resisting_friction_coefficient'] for c in moving],'o-',label='Solid + fluid')
        axes[0].semilogx(speeds,[-c['lower_solid_force_n_m']/hertz.line_load_n_m for c in moving],'o--',label='Solid')
        axes[0].semilogx(speeds,[-c['lower_fluid_force_n_m']/hertz.line_load_n_m for c in moving],'o--',label='Fluid')
        axes[1].semilogx(speeds,[c['fluid_load_fraction'] for c in moving],'o-')
    axes[0].set(xlabel='Lower plane speed [m/s]',ylabel='Resisting coefficient',title='Force on moving lower plane');axes[0].legend()
    axes[1].set(xlabel='Lower plane speed [m/s]',ylabel='Fluid load fraction',ylim=(0,1))
    for ax in axes:ax.grid(alpha=.2)
    fig.savefig(out/'mixed_friction.png',dpi=170);plt.close(fig)
    selected=next((i for i,c in enumerate(cases) if c['lower_wall_velocity_m_s']>0),len(cases)-1)
    v=solutions[selected]
    fig,axes=plt.subplots(1,2,figsize=(11,4),layout='constrained')
    for ax in axes:
        ax.plot(x,hertz.profile(elastic.x_m)*1e6,'--',color='grey',label='Unloaded mean cylinder')
        ax.plot(x,(hertz.profile(elastic.x_m)+v.rigid_offset_m)*1e6,':',label='Rigidly positioned cylinder')
        ax.plot(x,v.mean_gap_m*1e6,label='Deformed mean cylinder')
        ax.axhline(0,color='k',lw=1,label='Rigid plane')
        ax.set(xlabel='x [mm]',ylabel='z [um]');ax.grid(alpha=.2)
    radius=1.5*hertz.half_width_m
    central=abs(elastic.x_m)<=radius
    axes[1].set_xlim(-radius*1e3,radius*1e3)
    axes[1].set_ylim(-.05,max(float(np.max(v.mean_gap_m[central]))*1e6,.1))
    axes[0].set_title('Mean geometry; axes have different scales')
    axes[1].set_title(f'Contact-region detail: U_plane={cases[selected]["lower_wall_velocity_m_s"]:g} m/s')
    axes[0].legend(fontsize=8)
    fig.savefig(out/'mixed_geometry.png',dpi=170);plt.close(fig)
