"""Sequential, isolated mixed-lubrication loads; scalar workflows stay unchanged."""
from copy import deepcopy
import json
from pathlib import Path

from .geometry import line_load, transverse_width
from .interface import positive
from .macro_output import _csv


def batch_field(parameters):
    macro=parameters.get('macro_contact',{})
    fields=[key for key in ('line_load_n_m','normal_force_n') if isinstance(macro.get(key),(list,tuple))]
    if not fields:return None
    if len(fields)!=1 or (macro.get('line_load_n_m') is not None and macro.get('normal_force_n') is not None):
        raise ValueError('批量载荷只能填写 line_load_n_m 或 normal_force_n，不能同时给出')
    key=fields[0];values=macro[key]
    if not values:raise ValueError('载荷列表不能为空')
    checked=[positive(v,key) for v in values]
    if len(set(checked))!=len(checked):raise ValueError('载荷列表不能重复')
    return key


def run_mixed_batch(parameters,*,base_directory,save_plots,runner):
    from .mixed_output import prepare_mixed
    from .mixed_checks import write_mixed_workbook
    base=Path(base_directory).resolve();key=batch_field(parameters)
    name=parameters.get('output_directory')
    if not isinstance(name,str) or not name.strip():raise ValueError('Provide output_directory')
    out=Path(name);out=out if out.is_absolute() else base/out
    if out.exists() and (not out.is_dir() or any(out.iterdir())):
        raise FileExistsError(f'输出目录非空：{out}')
    jobs=[]
    for i,value in enumerate(parameters['macro_contact'][key]):
        p=deepcopy(parameters);p['macro_contact'][key]=positive(value,key)
        p['output_directory']=str(out/f'load_{i:03d}')
        prepare_mixed(p,base)  # validate every load before creating any result
        path=Path(p['psd']['path']);path=path if path.is_absolute() else base/path
        if not path.is_file():raise FileNotFoundError(path)
        jobs.append(p)
    manifest={'workflow_status':'running','scope':'Sequential independent mixed-lubrication loads',
        'parameters_snapshot':parameters,'loads':[]}
    def checkpoint():
        (out/'mixed_diagnostics.json').write_text(json.dumps(manifest,ensure_ascii=False,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    out.mkdir(parents=True);checkpoint();summaries=[];cases=[];rows=[]
    for i,p in enumerate(jobs):
        load=line_load(p);width=transverse_width(p)
        row={'load_index':i,'line_load_n_m':load,'normal_force_n':None if width is None else width*load,
            'directory':f'load_{i:03d}','status':'running'}
        manifest['loads'].append(row);checkpoint()
        print(f'正向载荷 {i+1}/{len(jobs)}：Wprime={load:g} N/m；{row["directory"]}',flush=True)
        try:s=runner(p,base_directory=base,save_plots=save_plots)
        except Exception as exc:
            row.update(status='failed',error=str(exc));manifest['workflow_status']='failed';checkpoint();raise
        summaries.append(s);row.update(status='complete',numerical_validation=s['numerical_validation']['status'])
        for c in s['cases']:
            combined={'load_index':i,'line_load_n_m':load,'normal_force_n':row['normal_force_n'],**c}
            if 'profile_file' in combined:combined['profile_file']=row['directory']+'/'+combined['profile_file']
            cases.append(combined)
        rows.append(dict(row));checkpoint()
    validation=deepcopy(summaries[0]['numerical_validation'])
    statuses=[s['numerical_validation']['status'] for s in summaries]
    validation['status']=next((v for v in ('needs_attention','incomplete','not_checked') if v in statuses),'passed_at_selected_cases')
    validation['scope']='Per-load independent numerical checks; inspect each load directory. Not a physical validation.'
    for name in ('comparisons','probes'):
        validation[name]=[{'load_index':i,**v} for i,s in enumerate(summaries) for v in s['numerical_validation'][name]]
    validation['warnings']=[f'load_{i:03d}: {w}' for i,s in enumerate(summaries) for w in s['numerical_validation']['warnings']]
    manifest.update(workflow_status='complete',cases=cases,numerical_validation=validation,
        solid_friction=summaries[0]['solid_friction'],outputs=['mixed_cases.csv','mixed_loads.csv','mixed_results.xlsx','mixed_diagnostics.json','mixed_numerical_checks.json'])
    _csv(out/'mixed_cases.csv',cases);_csv(out/'mixed_loads.csv',rows)
    (out/'mixed_numerical_checks.json').write_text(json.dumps(validation,ensure_ascii=False,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    if validation['comparisons']:
        _csv(out/'mixed_numerical_checks.csv',validation['comparisons'])
        manifest['outputs'].append('mixed_numerical_checks.csv')
    write_mixed_workbook(out/'mixed_results.xlsx',manifest)
    from openpyxl import load_workbook
    from .experiment_excel import _sheet
    book=load_workbook(out/'mixed_results.xlsx')
    _sheet(book,'载荷列表',list(rows[0]),[[r[k] for k in rows[0]] for r in rows])
    book.save(out/'mixed_results.xlsx');book.close();checkpoint()
    if save_plots:
        import matplotlib.pyplot as plt
        fig,axes=plt.subplots(1,2,figsize=(12,5))
        for i,s in enumerate(summaries):
            selected=[c for c in s['cases'] if c['lower_wall_velocity_m_s']!=0]
            if not selected:continue
            v=[abs(c['lower_wall_velocity_m_s']) for c in selected]
            label=f"W'={rows[i]['line_load_n_m']:g} N/m"
            for field,style in [('lower_resisting_friction_coefficient','-'),('lower_resisting_solid_coefficient','--'),('lower_resisting_fluid_coefficient',':')]:
                axes[0].semilogx(v,[c[field] for c in selected],style,color=f'C{i%10}',label=label+' '+field.replace('lower_resisting_',''))
            axes[1].semilogx(v,[c['fluid_load_fraction'] for c in selected],label=label)
        for ax in axes:ax.set_xlabel('Speed magnitude [m/s]');ax.grid(alpha=.2);ax.legend()
        axes[0].set_ylabel('Resisting friction coefficient');axes[1].set_ylabel('Fluid load fraction')
        fig.tight_layout();fig.savefig(out/'mixed_friction.png',dpi=160);plt.close(fig)
        manifest['outputs'].append('mixed_friction.png');checkpoint()
    return manifest
