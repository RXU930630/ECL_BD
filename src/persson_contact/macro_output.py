"""Independent dry module-6 workflow and mean-interface geometry exports."""
import csv
from dataclasses import asdict
import hashlib
import json
from pathlib import Path

import numpy as np

from .config import Options
from .core import PerssonModel
from .macroscopic import HalfSpaceLine,HertzLine,solve_smooth_contact
from .dry_contact import RoughContactLaw,solve_rough_contact
from .material import parse_material
from .interface import positive
from .psd import read_psd
from .table import build_contact_table


def _keys(raw,allowed,name):
    if not isinstance(raw,dict) or set(raw)-set(allowed):
        raise ValueError(f'Unknown or invalid {name} parameters')


def _csv(path,rows):
    with path.open('w',newline='',encoding='utf-8') as stream:
        w=csv.DictWriter(stream,fieldnames=list(rows[0]),lineterminator='\n');w.writeheader();w.writerows(rows)


def prepare_macro(parameters,base):
    raw=parameters.get('macro_contact',{})
    _keys(raw,{'geometry','geometry_parameters','elastic_model','body_representation','contact_model',
        'line_load_n_m','domain_half_width_m','reference_length_m','solver','contact_table'},'macro_contact')
    if raw.get('elastic_model','half_space_plane_strain')!='half_space_plane_strain':
        raise NotImplementedError('Only plane-strain elastic half-space response is implemented')
    if raw.get('geometry')!='cylinder_plane':raise NotImplementedError('Macro contact currently requires cylinder_plane')
    presets=raw.get('geometry_parameters',{})
    if not isinstance(presets,dict):raise ValueError('geometry_parameters must be a dictionary')
    geometry=presets.get('cylinder_plane',{})
    _keys(geometry,{'radius_m'},'cylinder_plane')
    representation=raw.get('body_representation','upper_elastic_lower_rigid')
    if representation not in ('upper_elastic_lower_rigid','relative_only'):
        raise NotImplementedError('Choose upper_elastic_lower_rigid or relative_only; separate two-body displacements are not implemented')
    model=raw.get('contact_model','persson')
    if model not in ('smooth','persson'):raise ValueError('contact_model must be smooth or persson')
    effective=parse_material(parameters.get('material',{}),base).elastic_modulus()
    hertz=HertzLine(geometry.get('radius_m'),raw.get('line_load_n_m'),effective)
    user_solver=raw.get('solver',{})
    _keys(user_solver,{'cells','relative_tolerance','max_iterations','edge_pressure_ratio_tolerance'},'macro solver')
    solver={'cells':257,'relative_tolerance':1e-9,'max_iterations':80,'edge_pressure_ratio_tolerance':1e-4,**user_solver}
    _keys(solver,{'cells','relative_tolerance','max_iterations','edge_pressure_ratio_tolerance'},'macro solver')
    for key in ('relative_tolerance','edge_pressure_ratio_tolerance'):
        solver[key]=positive(solver[key],key)
    if (type(solver['max_iterations']) is not int or solver['max_iterations']<1 or
        not solver['relative_tolerance']<.01 or not solver['edge_pressure_ratio_tolerance']<.01):
        raise ValueError('Invalid macro solver tolerance/iteration limit')
    elastic=HalfSpaceLine(raw.get('domain_half_width_m'),solver['cells'],effective,raw.get('reference_length_m'))
    if elastic.edges_m[-1]<=hertz.half_width_m:
        raise ValueError('Macro domain must extend beyond the Hertz half-width')
    return raw,solver,hertz,elastic,representation,model


def run_macro_project(parameters,*,base_directory,save_plots=True):
    base=Path(base_directory).resolve()
    raw,solver,hertz,elastic,representation,kind=prepare_macro(parameters,base)
    name=parameters.get('output_directory')
    if not isinstance(name,str) or not name.strip():raise ValueError('Provide output_directory')
    out=Path(name);out=out if out.is_absolute() else base/out
    if out.exists() and (not out.is_dir() or any(out.iterdir())):
        raise FileExistsError(f'输出目录非空，请使用新的 output_directory: {out}')
    table=None;psd_sha=None;closure=None
    if kind=='persson':
        psd_raw=parameters.get('psd',{})
        _keys(psd_raw,{'path','scale_c','tail_is_extrapolated'},'psd')
        path=Path(psd_raw['path']);path=path if path.is_absolute() else base/path
        psd=read_psd(path,psd_raw.get('scale_c',1.))
        psd_sha=hashlib.sha256(path.read_bytes()).hexdigest()
        options=Options(**parameters.get('options',{}))
        model=PerssonModel(psd,hertz.effective_modulus_pa,options.gamma,n_q=options.n_q,
            pdf_afact=options.pdf_afact,pdf_log_step=options.pdf_log_step,
            pdf_max_relative_contact_area=options.pdf_max_relative_contact_area,
            pdf_min_relative_contact_area=options.pdf_min_relative_contact_area,
            pdf_min_u1_m=options.pdf_min_u1_m,pdf_scale_method=options.pdf_scale_method)
        user_table=raw.get('contact_table',{})
        _keys(user_table,{'minimum_pressure_pa','maximum_pressure_pa','points'},'macro contact_table')
        settings={'minimum_pressure_pa':1e-100,'maximum_pressure_pa':2e6,'points':1201,**user_table}
        _keys(settings,{'minimum_pressure_pa','maximum_pressure_pa','points'},'macro contact_table')
        n=settings['points'];low=positive(settings['minimum_pressure_pa'],'minimum_pressure_pa');high=positive(settings['maximum_pressure_pa'],'maximum_pressure_pa')
        if type(n) is not int or not 101<=n<=20001 or not low<high:
            raise ValueError('Invalid macro contact-table range or point count')
        if high<=hertz.peak_pressure_pa:raise ValueError('Macro contact-table maximum must exceed the Hertz peak')
        if low*np.ptp(elastic.edges_m)/hertz.line_load_n_m>solver['relative_tolerance']:
            raise ValueError('Macro pressure-table minimum gives too large a far-tail load bound; lower minimum_pressure_pa')
        # Broad low-load tail plus a dense operating range; no extrapolation.
        focus=max(low,min(high,hertz.peak_pressure_pa*1e-3))
        pressures=np.unique(np.r_[np.geomspace(low,focus,n//3),np.geomspace(focus,high,n-n//3+1)])
        table=build_contact_table(model,pressures)
        result=solve_rough_contact(elastic,hertz,RoughContactLaw(table),
            tolerance=solver['relative_tolerance'],max_iterations=solver['max_iterations'])
        exact_gap=np.full_like(result.pressure_pa,np.nan)
        resolved=result.pressure_pa>0
        exact_gap[resolved]=model.mean_gap(result.pressure_pa[resolved])
        closure=exact_gap-result.mean_gap_m
    else:
        result=solve_smooth_contact(elastic,hertz,tolerance=solver['relative_tolerance'],max_iterations=solver['max_iterations'])
    if result.edge_pressure_ratio>solver['edge_pressure_ratio_tolerance']:
        raise ValueError('Macro pressure remains significant at domain edges; enlarge domain_half_width_m')
    if kind=='smooth':
        exact_gap=np.full_like(result.x_m,np.nan)
    rows=[]
    for i,x in enumerate(result.x_m):
        rows.append({'x_m':x,'cell_left_m':result.cell_edges_m[i],'cell_right_m':result.cell_edges_m[i+1],
            'initial_upper_surface_m':result.geometric_profile_m[i],
            'rigid_positioned_upper_surface_m':result.geometric_profile_m[i]+result.rigid_offset_m,
            'relative_elastic_displacement_m':result.displacement_m[i],
            'display_deformed_upper_surface_m':result.mean_gap_m[i],
            'display_lower_surface_m':0.,'mean_gap_m':result.mean_gap_m[i],
            'contact_pressure_pa':result.pressure_pa[i],'fluid_pressure_pa':0.,
            'relative_contact_area':result.relative_contact_area[i],
            'geometry_closure_error_m':result.mean_gap_m[i]-result.geometric_profile_m[i]-result.rigid_offset_m-result.displacement_m[i],
            'persson_mean_gap_m':None if table is None or not np.isfinite(exact_gap[i]) else exact_gap[i],
            'contact_table_gap_error_m':None if closure is None or not np.isfinite(closure[i]) else closure[i]})
    output_names=['macro_geometry.csv','macro_diagnostics.json']
    if table is not None:output_names.append('macro_contact_table.csv')
    if save_plots:output_names+=['macro_geometry.png','macro_geometry_true_scale.png','macro_response.png']
    scale=hertz.half_width_m**2/(2*hertz.radius_m)
    summary={'workflow_status':'complete','scope':'Dry 1D plane-strain cylinder/plane mean-interface contact. No fluid/thermal/friction coupling.',
        'contact_model':kind,'elastic_model':'half_space_plane_strain','body_representation':representation,
        'representation_note':('Upper elastic parabolic cylinder against rigid plane; E* is the one-body plane-strain modulus.'
            if representation=='upper_elastic_lower_rigid' else
            'Only combined relative compliance is known; displayed upper profile is equivalent, not a uniquely determined individual-body shape.'),
        'effective_modulus_pa':hertz.effective_modulus_pa,'radius_m':hertz.radius_m,
        'target_line_load_n_m':hertz.line_load_n_m,'actual_line_load_n_m':result.line_load_n_m,
        'load_relative_error':result.load_relative_error,'geometry_relative_residual':result.equilibrium_relative_residual,
        'pressure_min_pa':float(min(result.pressure_pa)),'pressure_max_pa':float(max(result.pressure_pa)),
        'minimum_mean_gap_m':result.minimum_gap_m,'rigid_offset_m':result.rigid_offset_m,
        'edge_pressure_ratio':result.edge_pressure_ratio,'iterations':result.iterations,'solver':solver,
        'domain_m':[float(elastic.edges_m[0]),float(elastic.edges_m[-1])],
        'reference_length_m':elastic.reference_length_m,
        'displacement_reference':'u=2/(pi E*) integral p(s) log(Lref/|x-s|) ds. Absolute u and h0 depend on Lref; pressure and mean gap do not.',
        'geometry_convention':'Unloaded upper z=g=x^2/(2R); rigidly positioned z=h0+g; deformed upper z=h0+g+u=h. Lower z=0. Positive u opens gap.',
        'microscopic_roughness':'Curves are mean surfaces, not a PSD phase reconstruction. Persson mean gap equals the macro gap; it is not added twice.',
        'hertz_reference':{'half_width_m':hertz.half_width_m,'peak_pressure_pa':hertz.peak_pressure_pa},
        'hertz_smooth_pressure_L1_relative_error':None if kind!='smooth' else float(np.sum(abs(result.pressure_pa-hertz.cell_pressure(elastic.edges_m)))*elastic.dx_m/hertz.line_load_n_m),
        'hertz_smooth_gap_max_error_over_scale':None if kind!='smooth' else float(np.max(abs(result.mean_gap_m-hertz.gap(result.x_m)))/scale),
        'contact_table_range_pa':None if table is None else [float(table.pressures_pa[0]),float(table.pressures_pa[-1])],
        'contact_table_points':None if table is None else len(table.pressures_pa),
        'contact_table_max_gap_error_m':None if closure is None else float(np.nanmax(abs(closure))),
        'contact_table_max_gap_error_over_hrms':None if closure is None else float(np.nanmax(abs(closure))/table.model.psd.hrms),
        'far_tail_zero_pressure_cells':0 if table is None else int(np.sum(result.pressure_pa==0)),
        'far_tail_fixed_gap_load_bound_n_m':None if table is None else float(table.pressures_pa[0]*np.ptp(elastic.edges_m)),
        'far_tail_policy':'Beyond the maximum table gap use zero pressure; omitted positive pressure is at most p_min per cell. The reported bound is at fixed gaps, not a global coupled error bound.',
        'assumptions':['Linear elastic half-space and parabolic local cylinder geometry; not finite thickness or whole-O-ring FEM.',
            'PSD describes residual roughness about the prescribed mean profile; spectral/macro scale separation is a model input.',
            'One-dimensional line load in N/m, not axisymmetric ball contact or total ring load.',
            'Solver residual, table error, spatial/domain convergence and constitutive-model uncertainty are different checks.'],
        'psd_sha256':psd_sha,'parameters_snapshot':parameters,'outputs':output_names}
    out.mkdir(parents=True,exist_ok=True)
    _csv(out/'macro_geometry.csv',rows)
    if table is not None:
        _csv(out/'macro_contact_table.csv',[{'pressure_pa':p,'mean_gap_m':h,'relative_contact_area':a}
            for p,h,a in zip(table.pressures_pa,table.mean_gaps_m,table.contact_fractions)])
    if save_plots:plot_macro(out,result,hertz,representation)
    (out/'macro_diagnostics.json').write_text(json.dumps(summary,ensure_ascii=False,indent=2,allow_nan=False)+'\n')
    print(f'宏观干接触已保存到 {out.resolve()}；载荷误差 {result.load_relative_error:.3g}，最小平均间隙 {result.minimum_gap_m:.6g} m。')
    return summary


def plot_macro(out,result,hertz,representation):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    x=result.x_m*1e3;h=result.mean_gap_m*1e6;g=result.geometric_profile_m*1e6
    fig,ax=plt.subplots(figsize=(10,4.6),layout='constrained')
    ax.plot(x,g,'--',color='grey',label='Unloaded mean profile g')
    ax.plot(x,g+result.rigid_offset_m*1e6,':',color='C1',label='Rigidly positioned g+h0')
    ax.plot(x,h,color='C0',lw=2,label='Deformed mean upper surface' if representation=='upper_elastic_lower_rigid' else 'Equivalent deformed relative profile')
    ax.axhline(0,color='black',lw=2,label='Rigid/reference lower surface')
    ax.fill_between(x,0,h,where=h>=0,alpha=.15,color='C0',label='Mean interfacial gap')
    points=ax.scatter(x,h,c=result.relative_contact_area,s=9,cmap='viridis',vmin=0,vmax=1,zorder=4)
    fig.colorbar(points,ax=ax,label='Relative real contact area')
    zoom=ax.inset_axes([.64,.48,.33,.36])
    zoom.plot(x,h,color='C0');zoom.axhline(0,color='black',lw=1)
    zoom.fill_between(x,0,h,where=h>=0,color='C0',alpha=.15)
    zoom.set_xlim(-1.3*hertz.half_width_m*1e3,1.3*hertz.half_width_m*1e3)
    zoom.set_ylim(-.02,max(5*max(result.minimum_gap_m,0)*1e6,.05*hertz.half_width_m**2/(2*hertz.radius_m)*1e6))
    zoom.set_title('Contact detail: z [um], x [mm]',fontsize=8);zoom.tick_params(labelsize=7);zoom.grid(alpha=.2)
    ax.set(xlabel='x [mm]',ylabel='z [um]',title='Dry macro deformation: different horizontal/vertical scales')
    ax.grid(alpha=.2);ax.legend(fontsize=8);fig.savefig(out/'macro_geometry.png',dpi=170);plt.close(fig)
    fig,ax=plt.subplots(figsize=(10,3.5),layout='constrained')
    ax.plot(x,result.geometric_profile_m*1e3,'--',color='grey',label='Unloaded mean profile')
    ax.plot(x,result.mean_gap_m*1e3,color='C0',label='Deformed mean/equivalent surface')
    ax.axhline(0,color='black',label='Lower reference surface')
    ax.set(xlabel='x [mm]',ylabel='z [mm]',title='True geometric aspect ratio (mean surfaces)')
    ax.set_aspect('equal',adjustable='box');ax.set_ylim(-.05*np.ptp(x),.2*np.ptp(x))
    ax.legend(fontsize=8);ax.grid(alpha=.2);fig.savefig(out/'macro_geometry_true_scale.png',dpi=170);plt.close(fig)
    fig,axes=plt.subplots(3,1,figsize=(10,9),sharex=True,layout='constrained')
    axes[0].plot(x,result.pressure_pa/1e6,label='Solved dry contact')
    axes[0].plot(x,hertz.cell_pressure(result.cell_edges_m)/1e6,'--',label='Smooth Hertz reference')
    axes[0].set_ylabel('Pressure [MPa]');axes[0].legend()
    axes[1].plot(x,result.mean_gap_m*1e6,label='Mean gap')
    axes[1].set_ylabel('Mean gap [um]')
    axes[2].plot(x,result.displacement_m*1e6,label='Relative elastic displacement')
    axes[2].set(xlabel='x [mm]',ylabel='Displacement [um]')
    for ax in axes:ax.grid(alpha=.2)
    fig.savefig(out/'macro_response.png',dpi=170);plt.close(fig)
