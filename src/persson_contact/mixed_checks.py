"""Independent mixed-film refinement checks; never alter the requested solution."""
from copy import deepcopy
import json
from pathlib import Path

import numpy as np

from .macro_output import _csv


def check_settings(raw):
    defaults={'enabled':False,'representative_cases':3,'relative_tolerance':.02,
              'absolute_mu_tolerance':1e-4,'absolute_load_fraction_tolerance':.001,
              'check_domain':True,'maximum_cells':1025}
    if not isinstance(raw,dict) or set(raw)-set(defaults):
        raise ValueError('Unknown mixed_friction.numerical_checks settings')
    settings={**defaults,**raw}
    for key in ('enabled','check_domain'):
        if type(settings[key]) is not bool:raise ValueError(f'numerical_checks.{key} must be True/False')
    for key,low,high in [('representative_cases',1,20),('maximum_cells',17,2049)]:
        value=settings[key]
        if type(value) is not int or not low<=value<=high:raise ValueError(f'Invalid numerical_checks.{key}')
    for key in ('relative_tolerance','absolute_mu_tolerance','absolute_load_fraction_tolerance'):
        value=settings[key]
        if (isinstance(value,bool) or not isinstance(value,(int,float,np.integer,np.floating))
                or not np.isfinite(value) or not 0<float(value)<1):
            raise ValueError(f'Invalid numerical_checks.{key}')
        settings[key]=float(value)
    return settings


def representative_cases(cases,count):
    unique={float(c['lower_wall_velocity_m_s']):c for c in reversed(cases)}
    ordered=sorted(unique.values(),key=lambda c:(abs(c['lower_wall_velocity_m_s']),c['lower_wall_velocity_m_s']))
    indices=np.unique(np.rint(np.linspace(0,len(ordered)-1,min(count,len(ordered)))).astype(int))
    return [ordered[i] for i in indices]


def compare_cases(baseline,probe,settings,load):
    """Mixed absolute/relative limits; tiny areas use sensitivity at 1 MPa."""
    metrics={
        'real_contact_area_per_width_m':load*settings['absolute_mu_tolerance']/1e6,
        'minimum_mean_gap_m':1e-12,
        'maximum_fluid_pressure_pa':1e-6,
        'fluid_load_fraction':settings['absolute_load_fraction_tolerance'],
        **{key:settings['absolute_mu_tolerance'] for key in (
            'lower_resisting_friction_coefficient','lower_resisting_solid_coefficient',
            'lower_resisting_fluid_sliding_coefficient','lower_resisting_fluid_pressure_coefficient')}}
    rows=[]
    for a in baseline:
        matches=[b for b in probe if b['lower_wall_velocity_m_s']==a['lower_wall_velocity_m_s']]
        if len(matches)!=1:raise ValueError('Refinement output has missing or duplicate velocities')
        b=matches[0]
        for key,absolute in metrics.items():
            if a[key] is None and b[key] is None:continue  # zero-speed friction is undefined
            av,bv=float(a[key]),float(b[key]);scale=max(abs(av),abs(bv))
            difference=abs(bv-av);limit=absolute+settings['relative_tolerance']*scale
            rows.append({'lower_wall_velocity_m_s':float(a['lower_wall_velocity_m_s']),
                'quantity':key,'baseline_value':av,'check_value':bv,'absolute_difference':difference,
                'relative_difference':None if scale==0 else difference/scale,
                'allowed_difference':limit,'within_tolerance':bool(difference<=limit)})
    return rows


def run_checks(parameters,base_directory,summary,cache,settings,runner,out):
    report={'status':'not_checked','settings':settings,'comparisons':[],'probes':[],
        'selected_velocities_m_s':[],
        'scope':'One refinement per check, selected velocities only; not a rigorous error bound or physical-model validation.',
        'primary_result_policy':'Requested grid/domain/settings retained. Checks are separate runs; no automatic geometry changes.',
        'warnings':[]}
    if settings['enabled']:
        selected=representative_cases(summary['cases'],settings['representative_cases'])
        speeds=[c['lower_wall_velocity_m_s'] for c in selected]
        report['selected_velocities_m_s']=speeds
        macro=parameters['macro_contact'];cells=macro.get('solver',{}).get('cells',257)
        physical=macro['geometry']=='finite_cylindrical_pad'
        report['boundary_kind']='physical_pad_edges' if physical else 'numerical_far_field_truncation'
        boundary_bad=summary['far_field_edge_pressure_check_applied'] and any(
            c['edge_pressure_ratio']>macro.get('solver',{}).get('edge_pressure_ratio_tolerance',1e-4)
            for c in summary['cases'])
        if boundary_bad:report['warnings'].append('主计算的远场边缘压力偏高；保留结果供诊断，不能视为区间收敛。')
        for kind in ('mesh','tables_and_spectral_quadrature','domain'):
            p=deepcopy(parameters);m=p['macro_contact'];mixed=p['mixed_friction']
            p['options'].setdefault('pdf_log_step',.01)
            mixed['lower_wall_velocities_m_s']=speeds
            mixed['numerical_checks']={'enabled':False}
            target=out/'numerical_checks'/kind;p['output_directory']=str(target.resolve())
            numerical=m.setdefault('solver',{})
            probe_cache=None
            record={'check':kind,'status':'pending','directory':str(target)}
            report['probes'].append(record)
            if kind=='domain' and physical:
                record.update(status='not_applicable',reason='Physical pad edges are fixed; actual geometry is never enlarged.')
                continue
            if kind=='domain' and (not settings['check_domain'] or
                    parameters['mixed_friction'].get('inlet_pressure_pa',0.)!=0. or
                    parameters['mixed_friction'].get('outlet_pressure_pa',0.)!=0.):
                record.update(status='not_checked',reason='Domain check disabled or nonzero reservoir pressures; not a zero-pressure far-field test.')
                continue
            if kind in ('mesh','domain'):
                requested=2*cells-1 if kind=='mesh' else 2*cells+1
                if requested>settings['maximum_cells']:
                    record.update(status='not_checked',reason=f'Requires {requested} cells, exceeding maximum_cells; no false pass at unchanged resolution.')
                    continue
                numerical['cells']=requested
                probe_cache=dict(cache)  # same constitutive model, different spatial grid only
                if kind=='domain':
                    half=m.get('domain_half_width_m')
                    if half is None:
                        half=float(summary['geometry']['surface_interval_m'][1])
                    m['domain_half_width_m']=2*half
                    if m.get('reference_length_m') is not None:
                        m['reference_length_m']=max(m['reference_length_m'],2*m['domain_half_width_m'])
            else:
                table=m.setdefault('contact_table',{})
                table['points']=2*table.get('points',1201)-1
                mixed['transport_points']=2*mixed.get('transport_points',161)-1
                p['options']['n_q']*=2;p['options']['pdf_log_step']/=2
                if table['points']>20001 or mixed['transport_points']>2001:
                    record.update(status='not_checked',reason='Refined tables exceed supported limits; unchanged tables are not a convergence check.')
                    continue
            record['settings']={'cells':numerical.get('cells',257),
                'domain_half_width_m':m.get('domain_half_width_m'),
                'contact_table_points':m.get('contact_table',{}).get('points',1201),
                'transport_points':mixed.get('transport_points',161),
                'n_q':p['options']['n_q'],'pdf_log_step':p['options']['pdf_log_step']}
            print(f'数值核查 {kind}：{len(speeds)} 个速度，主计算结果保持不变。',flush=True)
            try:
                result=runner(p,base_directory=base_directory,save_plots=False,
                    _cache=probe_cache,_defer_edge=True)
                rows=compare_cases(selected,result['cases'],settings,summary['target_line_load_n_m'])
                for row in rows:row['check']=kind
                report['comparisons'].extend(rows)
                record.update(status='compared',within_tolerance=all(row['within_tolerance'] for row in rows))
            except (ValueError,NotImplementedError,RuntimeError) as exc:
                record.update(status='failed',error=str(exc))
                report['warnings'].append(f'{kind} 核查未完成：{exc}')
        completed=[r for r in report['probes'] if r['status']=='compared']
        failed=boundary_bad or any(r['status']=='failed' or (r['status']=='compared' and not r['within_tolerance']) for r in report['probes'])
        incomplete=any(r['status']=='not_checked' for r in report['probes'])
        report['status']='needs_attention' if failed else 'incomplete' if incomplete else 'passed_at_selected_cases'
        if any(not r['within_tolerance'] for r in completed):
            report['warnings'].append('加密或扩大区间后有结果变化超过设定容限；请查看差异表，调整主计算设置后重新核查。')
        if incomplete:report['warnings'].append('部分核查未执行；不能将其视为通过。')
        if physical:report['warnings'].append('有限胶块只验证固定物理边界内的数值精度；这不验证边界条件本身是否符合实验。')
        print(f'数值核查状态：{report["status"]}。详见 mixed_numerical_checks.json / mixed_results.xlsx。',flush=True)
        for message in report['warnings']:print('提示：'+message,flush=True)
    (out/'mixed_numerical_checks.json').write_text(json.dumps(report,ensure_ascii=False,indent=2,allow_nan=False)+'\n', encoding='utf-8')
    if report['comparisons']:_csv(out/'mixed_numerical_checks.csv',report['comparisons'])
    return report


def write_mixed_workbook(path,summary):
    from .experiment_excel import _workbook,_sheet
    book=_workbook();cases=summary['cases'];validation=summary['numerical_validation']
    _sheet(book,'计算概况',('项目','内容'),[
        ('数值核查状态',validation['status']),('解释','非线性方程已求解，不等于物理模型已经验证。'),
        ('主结果','输入网格和物理几何的正向解；核查结果不自动替换主结果。'),
        ('剪切强度模型',json.dumps(summary['solid_friction'],ensure_ascii=False)),
        ('shear flow','混合润滑内部自动计算，COMPUTE_SHEAR_FLOW无需另开。'),
        ('面积单位','real_contact_area_per_width_m为Ar/B；乘轴向宽度B得到m²。'),
        ('核查范围',validation['scope']),
        *[('提示',v) for v in validation['warnings']]])
    _sheet(book,'正向摩擦与面积',list(cases[0]),[[c[k] for k in cases[0]] for c in cases])
    rows=validation['comparisons']
    if rows:_sheet(book,'数值核查差异',list(rows[0]),[[r[k] for k in rows[0]] for r in rows])
    _sheet(book,'核查执行情况',('check','status','details'),[
        (p['check'],p['status'],json.dumps(p,ensure_ascii=False)) for p in validation['probes']])
    for sheet in book:
        for row in sheet:
            for cell in row:
                if isinstance(cell.value,str):cell.data_type='s'
    book.save(path);book.close()
