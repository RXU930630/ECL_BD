"""Dispatch contact, interface laws, pressure flow and fixed-state leakage."""
import hashlib
import json
from pathlib import Path

from .applications import LeakageParameters
from .cli import run_contact as calculate_contact
from .config import problem_from_mapping
from .core import PerssonModel
from .flow import FlowOptions
from .flow_output import run_pressure_flow
from .interface import make_interface


def run_project(parameters, *, base_directory, run_contact=True, run_leakage=False,
                run_mixed_friction=False, compute_gap_pdf=True, save_plots=True,
                compute_pressure_flow=False):
    switches = (run_contact, run_leakage, run_mixed_friction, compute_gap_pdf,
                save_plots,compute_pressure_flow)
    if any(type(value) is not bool for value in switches):
        raise ValueError('所有功能开关必须为 True 或 False。')
    if not any(switches[:3]):
        print('没有启用计算功能；未创建任何结果文件。')
        return None
    base = Path(base_directory).resolve()
    problem_keys = {'description', 'psd', 'material', 'pressures_pa', 'contact_pressure', 'options'}
    allowed = problem_keys | {'output_directory', 'leakage', 'mixed_friction','interface','pressure_flow'}
    if not isinstance(parameters, dict) or set(parameters) - allowed:
        raise ValueError('parameters.py 中存在未知参数分组。')
    for name in ('leakage', 'mixed_friction','interface','pressure_flow'):
        if name in parameters and not isinstance(parameters[name], dict):
            raise ValueError(f'{name} 必须为参数字典。')
    if run_mixed_friction:
        raise NotImplementedError('混合摩擦 RUN_MIXED_FRICTION 尚未实现。')
    flow_enabled=compute_pressure_flow or run_leakage
    interface=flow_options=leakage=None
    # Check application inputs before PSD work or creation of any output.
    if flow_enabled:
        interface=make_interface(parameters.get('interface',{}))
        raw_flow=parameters.get('pressure_flow',{})
        unknown=set(raw_flow)-set(FlowOptions.__dataclass_fields__)
        if unknown:
            raise ValueError(f'Unknown pressure_flow parameters: {sorted(unknown)}')
        flow_options=FlowOptions(**raw_flow)
    if run_leakage:
        raw=dict(parameters.get('leakage',{}))
        legacy=raw.pop('fluid_viscosity_pa_s',None)
        if legacy is not None:
            raise ValueError('请仅在 interface.viscosity_pa_s 中给出黏度。')
        if set(raw)!=set(LeakageParameters.__dataclass_fields__):
            raise ValueError('leakage 需要 seal_length_m、seal_width_m、inlet_pressure_pa、outlet_pressure_pa。')
        leakage=LeakageParameters(**raw)
    problem, options = problem_from_mapping({key: value for key, value in parameters.items()
                                           if key in problem_keys}, base)
    name = parameters.get('output_directory')
    if not isinstance(name, str) or not name.strip():
        raise ValueError('请在 parameters.py 中给出 output_directory。')
    output = Path(name)
    if not output.is_absolute():
        output = base / output
    psd = Path(parameters['psd']['path'])
    if not psd.is_absolute():
        psd = base / psd
    snapshot = json.dumps(parameters, ensure_ascii=False, sort_keys=True, allow_nan=False)
    summary=calculate_contact(problem, options, output, plot=save_plots,
        compute_gap_pdf=compute_gap_pdf, provenance={
            'config_format': 'Python parameters.py','parameters_snapshot': parameters,
            'parameters_sha256': hashlib.sha256(snapshot.encode('utf-8')).hexdigest(),
            'psd_sha256': hashlib.sha256(psd.read_bytes()).hexdigest(),
            'material_input': parameters['material'],
            'functions': {'contact': True,'pressure_flow':flow_enabled,
                          'leakage':run_leakage,'mixed_friction':False},
            'requested_contact_output':run_contact,
            'workflow_status':'contact_complete_flow_pending' if flow_enabled else 'complete',
        })
    if flow_enabled:
        model=PerssonModel(problem.psd,problem.effective_modulus_pa,options.gamma,
            n_q=options.n_q,pdf_afact=options.pdf_afact,pdf_log_step=options.pdf_log_step,
            pdf_max_relative_contact_area=options.pdf_max_relative_contact_area,
            pdf_min_relative_contact_area=options.pdf_min_relative_contact_area,
            pdf_min_u1_m=options.pdf_min_u1_m)
        print('计算压力流动系数；PDF 导出开关不影响内部所需的分布计算。')
        try:
            flow=run_pressure_flow(model,problem.pressures_pa,interface,flow_options,output,
                percolation=options.critical_relative_contact_area,leakage=leakage,plot=save_plots)
        except Exception as exc:
            summary['workflow_status']='flow_failed'
            summary['flow_error']=str(exc)
            (output/'diagnostics.json').write_text(json.dumps(summary,indent=2,allow_nan=False)+'\n')
            raise
        summary['workflow_status']='complete'
        summary['pressure_flow']={'file':'pressure_flow_table.csv','diagnostics_file':'flow_diagnostics.json',
            'status_counts':flow['status_counts'],'leakage_file':'leakage_table.csv' if run_leakage else None}
        (output/'diagnostics.json').write_text(json.dumps(summary,indent=2,allow_nan=False)+'\n')
        print(f'压力流动结果已保存到 {output.resolve()}')
    return summary
