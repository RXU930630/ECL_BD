"""Experiment-to-shear-strength task, with fresh or reusable mixed solutions."""
from copy import deepcopy
from dataclasses import asdict
import csv
import hashlib
import json
from pathlib import Path

import numpy as np

from .contact_friction import finite
from .macro_output import _csv,_keys
from .mixed_output import prepare_mixed,run_mixed_project
from .shear_fit import read_experiments,fit_constant_strength
from .geometry import transverse_width,line_load


def _path(name,base):
    if not isinstance(name,str) or not name.strip():raise ValueError('Provide a nonempty file/directory path')
    p=Path(name);return (p if p.is_absolute() else base/p).resolve()


def _signature(diag):
    """Saved physical model is authoritative; current PSD/material is not read."""
    p=diag.get('parameters_snapshot',{});macro=p.get('macro_contact',{});mixed=p.get('mixed_friction',{})
    selected=macro.get('geometry')
    geometry={k:macro.get(k) for k in ('geometry','elastic_model','body_representation','contact_model')}
    geometry['geometry_parameters']=macro.get('geometry_parameters',{}).get(selected,{})
    if selected!='finite_cylindrical_pad':geometry['domain_half_width_m']=macro.get('domain_half_width_m')
    return json.dumps({'material':p.get('material'),'psd_sha256':diag.get('psd_sha256'),
        'psd_settings':{k:v for k,v in p.get('psd',{}).items() if k!='path'},
        'interface':p.get('interface'),'options':p.get('options'),
        'pressure_flow':p.get('pressure_flow'),
        'macro':geometry,
        'mixed':{k:mixed.get(k) for k in ('rough_surface','inlet_pressure_pa','outlet_pressure_pa','cavitation_model','cavitation_friction')},
        'fluid_friction':p.get('fluid_friction')},sort_keys=True,allow_nan=False)


def read_forward_run(directory):
    """Read successful native mixed outputs; check CSV against diagnostics.

    Use nominal total imposed line load for both fluid coefficient and real
    area/load sensitivity. Never use the solid-only normal load or cached
    'contact coefficient', which would have a different denominator.
    """
    directory=Path(directory)
    diagnostic=directory/'mixed_diagnostics.json';table=directory/'mixed_cases.csv'
    diag=json.loads(diagnostic.read_text(encoding='utf-8'))
    if diag.get('workflow_status')!='complete' or not isinstance(diag.get('cases'),list):
        raise ValueError(f'混合润滑结果未完成：{directory}')
    load=finite(diag.get('target_line_load_n_m'),'saved target line load',minimum=0,strict=True)
    with table.open(encoding='utf-8-sig',newline='') as stream:csv_rows=list(csv.DictReader(stream))
    if len(csv_rows)!=len(diag['cases']):raise ValueError('Saved mixed CSV and diagnostics have different case counts')
    rows=[]
    keys=('lower_wall_velocity_m_s','real_contact_area_per_width_m','lower_fluid_force_n_m',
        'contact_line_load_n_m','fluid_line_load_n_m','load_relative_error','geometry_relative_residual','mass_relative_residual')
    for case,raw in zip(diag['cases'],csv_rows):
        values={key:finite(case.get(key),'saved '+key) for key in keys}
        for key,v in values.items():
            if not np.isclose(v,finite(raw.get(key),'saved CSV '+key),rtol=1e-12,atol=1e-30):
                raise ValueError(f'Saved mixed CSV and diagnostics disagree: {key}')
        if (values['real_contact_area_per_width_m']<0 or
            min(values['contact_line_load_n_m'],values['fluid_line_load_n_m'])<0 or
            abs((values['contact_line_load_n_m']+values['fluid_line_load_n_m'])/load-1)>1e-6 or
            any(abs(values[k])>1e-6 for k in ('load_relative_error','geometry_relative_residual','mass_relative_residual'))):
            raise ValueError('Saved forward solution fails area/load/conservation checks; do not fit this result')
        speed=values['lower_wall_velocity_m_s']
        if speed==0:continue
        area=values['real_contact_area_per_width_m']
        rows.append({'lower_wall_velocity_m_s':speed,'line_load_n_m':load,
            'real_contact_area_per_width_m':area,'area_over_total_load_pa_inv':area/load,
            'fluid_coefficient':float(-np.sign(speed)*values['lower_fluid_force_n_m']/load),
            'source_directory':str(directory.resolve())})
    provenance={'directory':str(directory.resolve()),'diagnostics_sha256':hashlib.sha256(diagnostic.read_bytes()).hexdigest(),
        'cases_sha256':hashlib.sha256(table.read_bytes()).hexdigest(),'psd_sha256':diag.get('psd_sha256'),
        'synthetic_demo':bool(diag.get('synthetic_demo',False)),
        'parameters_snapshot':diag.get('parameters_snapshot'),'limitations':diag.get('limitations',[])}
    return rows,provenance,_signature(diag)


def prepare_fit(parameters,base):
    # Both input and the Excel result workbook need this dependency. Check
    # before any expensive forward solve, including when input is legacy CSV.
    import openpyxl
    raw=parameters.get('shear_strength_fit',{})
    _keys(raw,{'model','experiment_file','source','source_parameters','reference_strength_pa',
        'sliding_direction','weighting','speed_range_m_s','minimum_solid_coefficient_at_1mpa'},'shear_strength_fit')
    if raw.get('model','constant_shear_strength')!='constant_shear_strength':raise NotImplementedError('Only constant shear-strength fitting is implemented')
    source=raw.get('source','simulate')
    if source not in ('simulate','saved_forward'):raise ValueError('fit source must be simulate or saved_forward')
    presets=raw.get('source_parameters',{})
    if not isinstance(presets,dict):raise ValueError('source_parameters must be a dictionary')
    active=presets.get(source,{})
    _keys(active,set() if source=='simulate' else {'directories'},'fit source_parameters')
    direction=raw.get('sliding_direction',1)
    if type(direction) is not int or direction not in (-1,1):raise ValueError('sliding_direction must be +1 or -1')
    ref=finite(raw.get('reference_strength_pa',1e6),'reference_strength_pa',minimum=0,strict=True)
    threshold=finite(raw.get('minimum_solid_coefficient_at_1mpa',1e-5),'minimum_solid_coefficient_at_1mpa',minimum=0)
    weighting=raw.get('weighting','equal')
    if weighting not in ('equal','standard_uncertainty_mu'):raise ValueError('Invalid fit weighting')
    speed_range=raw.get('speed_range_m_s',[None,None])
    if not isinstance(speed_range,(list,tuple)) or len(speed_range)!=2:raise ValueError('Provide fit speed_range_m_s=[minimum,maximum]')
    speed_range=[None if v is None else finite(v,'fit speed range',minimum=0,strict=True) for v in speed_range]
    if all(v is not None for v in speed_range) and speed_range[0]>speed_range[1]:raise ValueError('Invalid fit speed range')
    experiment_path=_path(raw.get('experiment_file'),base)
    width=transverse_width(parameters)
    macro=parameters.get('macro_contact',{})
    default_load=line_load(parameters) if macro.get('line_load_n_m') is not None or macro.get('normal_force_n') is not None else None
    experiments=read_experiments(experiment_path,
        default_line_load_n_m=default_load,transverse_width_m=width)
    # Validate membership and measurement uncertainties before expensive work.
    fit_constant_strength(experiments,np.ones(len(experiments))/ref,np.zeros(len(experiments)),
        reference_strength_pa=ref,weighting=weighting,speed_range_m_s=speed_range,
        minimum_solid_coefficient_at_1mpa=0.)
    directories=[];jobs=[]
    if source=='saved_forward':
        names=active.get('directories')
        if not isinstance(names,list) or not names:raise ValueError('saved_forward requires a nonempty directories list')
        directories=[_path(name,base) for name in names]
        if len(set(directories))!=len(directories):raise ValueError('Duplicate saved forward directory')
    else:
        loads=[]
        for value in sorted(set(e.line_load_n_m for e in experiments)):
            if not loads or not np.isclose(value,loads[-1],rtol=1e-9,atol=0.):loads.append(value)
        for load in loads:
            p=deepcopy(parameters)
            p['macro_contact'].pop('normal_force_n',None)
            p['macro_contact']['line_load_n_m']=load
            velocities=[]
            for value in sorted(set(e.speed_m_s for e in experiments if np.isclose(e.line_load_n_m,load,rtol=1e-9,atol=0.))):
                if not velocities or not np.isclose(value,velocities[-1],rtol=1e-9,atol=0.):velocities.append(value)
            p['mixed_friction']['lower_wall_velocities_m_s']=[direction*v for v in velocities]
            p['contact_friction']['model']='constant_shear_strength'
            p['contact_friction']['model_parameters']={'constant_shear_strength':{'shear_strength_pa':ref}}
            prepare_mixed(p,base)
            psd=_path(p['psd']['path'],base)
            if not psd.is_file():raise FileNotFoundError(psd)
            jobs.append(p)
    return raw,source,experiments,experiment_path,width,ref,direction,weighting,threshold,speed_range,directories,jobs


def _match(experiments,responses,direction):
    matched=[]
    for e in experiments:
        hits=[r for r in responses if np.isclose(r['line_load_n_m'],e.line_load_n_m,rtol=1e-9,atol=0.) and
              np.isclose(r['lower_wall_velocity_m_s'],direction*e.speed_m_s,rtol=1e-9,atol=0.)]
        if len(hits)!=1:raise ValueError(f'实验第{e.row_number}行需要唯一正向解：speed={e.speed_m_s:g} m/s、line_load={e.line_load_n_m:g} N/m；找到{len(hits)}个。不自动插值，请按实验速度/载荷重算。')
        matched.append(hits[0])
    return matched


def run_shear_fit_project(parameters,*,base_directory,save_plots=True):
    base=Path(base_directory).resolve()
    raw,source,experiments,experiment_file,width,ref,direction,weighting,threshold,speed_range,directories,jobs=prepare_fit(parameters,base)
    out=_path(parameters.get('output_directory'),base)
    if out.exists() and (not out.is_dir() or any(out.iterdir())):raise FileExistsError(f'输出目录非空，请选择新目录：{out}')
    responses=[];provenance=[];signatures=set()
    def collect(directory):
        rows,origin,signature=read_forward_run(directory)
        responses.extend(rows);provenance.append(origin);signatures.add(signature)
        if len(signatures)>1:raise ValueError('Saved forward runs use different material/PSD/fluid/geometry/closure settings; fit them separately')
    if source=='saved_forward':
        for directory in directories:collect(directory)
        matched=_match(experiments,responses,direction)
    else:
        out.mkdir(parents=True,exist_ok=True)
        try:
            for i,p in enumerate(jobs):
                target=out/'forward'/f'load_{i:03d}';p['output_directory']=str(target)
                print(f'拟合参考计算 {i+1}/{len(jobs)}：Wprime={p["macro_contact"]["line_load_n_m"]:g} N/m，tau_ref={ref:g} Pa',flush=True)
                run_mixed_project(p,base_directory=base,save_plots=save_plots)
                collect(target)
            matched=_match(experiments,responses,direction)
        except Exception as exc:
            (out/'fit_diagnostics.json').write_text(json.dumps({'workflow_status':'forward_failed','error':str(exc),
                'completed_forward_directories':[p['directory'] for p in provenance],
                'note':'No fitted material parameter was produced. Successful forward runs may be reused with source=saved_forward.'},indent=2,ensure_ascii=False)+'\n')
            raise
    try:
        summary,records=fit_constant_strength(experiments,[r['area_over_total_load_pa_inv'] for r in matched],
            [r['fluid_coefficient'] for r in matched],reference_strength_pa=ref,weighting=weighting,
            minimum_solid_coefficient_at_1mpa=threshold,speed_range_m_s=speed_range)
    except Exception as exc:
        if source=='simulate':
            (out/'fit_diagnostics.json').write_text(json.dumps({'workflow_status':'fit_failed','error':str(exc),
                'completed_forward_directories':[p['directory'] for p in provenance],
                'note':'No fitted material parameter was produced; successful forward runs may be reused.'},indent=2,ensure_ascii=False)+'\n')
        raise
    for record,response in zip(records,matched):
        record.update(real_contact_area_per_width_m=response['real_contact_area_per_width_m'],
            real_contact_area_m2=None if width is None else response['real_contact_area_per_width_m']*width,
            total_normal_force_n=None if width is None else record['line_load_n_m']*width)
    snippet={'model':'constant_shear_strength','model_parameters':{'constant_shear_strength':{'shear_strength_pa':summary['fitted_shear_strength_pa']}}}
    outputs=['fit_results.xlsx','fit_cases.csv','fit_summary.json','fit_diagnostics.json','fitted_contact_friction.json','forward_response.csv','experiment_used.csv']
    if save_plots:outputs.append('shear_strength_fit.png')
    diagnostic={'workflow_status':'complete','task':'constant_interface_shear_strength_fit',
        'equation':'mu = mu_fluid + mu_background + tau * (A_real / TOTAL normal load)',
        'source':source,'sliding_direction':direction,'transverse_width_m':width,
        'normal_state_recalculated_during_linear_fit':False,
        'synthetic_demo':any(p['synthetic_demo'] for p in provenance),
        'experiment_file':str(experiment_file),'experiment_sha256':hashlib.sha256(experiment_file.read_bytes()).hexdigest(),
        'source_forward_runs':provenance,'parameters_snapshot':parameters,'fit':summary,'outputs':outputs}
    out.mkdir(parents=True,exist_ok=True)
    _csv(out/'fit_cases.csv',records);_csv(out/'forward_response.csv',responses)
    _csv(out/'experiment_used.csv',[{k:v for k,v in asdict(e).items() if k!='row_number'} for e in experiments])
    from .experiment_excel import write_fit_workbook
    write_fit_workbook(out/'fit_results.xlsx',summary,records,responses,
        [{k:v for k,v in asdict(e).items() if k!='row_number'} for e in experiments])
    for name,data in [('fit_summary.json',summary),('fitted_contact_friction.json',snippet),('fit_diagnostics.json',diagnostic)]:
        (out/name).write_text(json.dumps(data,indent=2,ensure_ascii=False,allow_nan=False)+'\n')
    if save_plots:plot_fit(out,records)
    print(f'拟合完成：tau={summary["fitted_shear_strength_pa"]:.8g} Pa ({summary["fitted_shear_strength_pa"]/1e6:.6g} MPa)，{summary["training_points"]}个参与点。输出：{out}',flush=True)
    if summary['nonnegative_boundary_active']:
        print('提示：最优解落在非负边界0；流体/背景可能高于实验，请检查分项。这不表示真实剪切强度已精确确定为0。',flush=True)
    return diagnostic


def plot_fit(out,records):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig,axes=plt.subplots(1,2,figsize=(11,4.5),layout='constrained')
    loads=sorted(set(r['line_load_n_m'] for r in records))
    for load in loads:
        group=sorted([r for r in records if r['line_load_n_m']==load],key=lambda r:r['speed_m_s'])
        v=[r['speed_m_s'] for r in group];label=f"W'={load:g} N/m"
        line=axes[0].semilogx(v,[r['predicted_total_coefficient'] for r in group],'-',label=label+' fitted total')[0]
        color=line.get_color()
        axes[0].semilogx(v,[r['fitted_solid_coefficient'] for r in group],'--',color=color,label=label+' solid')
        axes[0].semilogx(v,[r['fluid_coefficient']+r['background_friction_coefficient'] for r in group],':',color=color,label=label+' fluid + supplied background')
        for used,marker in ((True,'o'),(False,'s')):
            points=[r for r in group if r['used_for_fit']==used]
            if not points:continue
            axes[0].scatter([r['speed_m_s'] for r in points],[r['friction_coefficient'] for r in points],
                marker=marker,edgecolors=color,facecolors=color if used else 'none',label=label+(' training' if used else 'held out/excluded'))
        axes[1].semilogx(v,[r['residual_experiment_minus_model'] for r in group],'o',color=color,label=label)
    axes[0].set(xlabel='Speed [m/s]',ylabel='Resisting friction coefficient',title='Constant interface shear strength')
    axes[1].set(xlabel='Speed [m/s]',ylabel='Experiment - prediction',title='Residuals; no fitted offset');axes[1].axhline(0,color='k',lw=.8)
    for ax in axes:ax.grid(alpha=.2);ax.legend(fontsize=7)
    fig.savefig(Path(out)/'shear_strength_fit.png',dpi=160);plt.close(fig)
