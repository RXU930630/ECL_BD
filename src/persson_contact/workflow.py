"""Dispatch contact, liquid flow, equivalent gas estimates and ballistic leakage."""
import hashlib
import json
from dataclasses import replace
from pathlib import Path

from .applications import LeakageParameters
from .cli import run_contact as calculate_contact
from .config import problem_from_mapping
from .core import PerssonModel
from .flow import FlowOptions
from .flow_output import run_pressure_flow
from .interface import make_interface
from .leakage import prepare_leakage, required_pressures, make_transport, run_spatial_leakage, validate_liquid_estimates
from .gas import BallisticGas
from .shear import read_shear_options,ShearFlowLaw
from .shear_output import run_shear_flow


def run_project(parameters, *, base_directory, run_contact=True, run_leakage=False,
                run_mixed_friction=False, compute_gap_pdf=True, save_plots=True,
                compute_pressure_flow=False,compute_shear_flow=False):
    switches = (run_contact, run_leakage, run_mixed_friction, compute_gap_pdf,
                save_plots,compute_pressure_flow,compute_shear_flow)
    if any(type(value) is not bool for value in switches):
        raise ValueError('所有功能开关必须为 True 或 False。')
    if not any(switches[:3]) and not compute_shear_flow:
        print('没有启用计算功能；未创建任何结果文件。')
        return None
    base = Path(base_directory).resolve()
    problem_keys = {'description', 'psd', 'material', 'pressures_pa', 'contact_pressure', 'options'}
    allowed = problem_keys | {'output_directory', 'leakage', 'mixed_friction','interface','pressure_flow','shear_flow'}
    if not isinstance(parameters, dict) or set(parameters) - allowed:
        raise ValueError('parameters.py 中存在未知参数分组。')
    for name in ('leakage', 'mixed_friction','interface','pressure_flow','shear_flow'):
        if name in parameters and not isinstance(parameters[name], dict):
            raise ValueError(f'{name} 必须为参数字典。')
    if run_mixed_friction:
        raise NotImplementedError('混合摩擦 RUN_MIXED_FRICTION 尚未实现。')
    flow_enabled=compute_pressure_flow or run_leakage or compute_shear_flow
    shear_options=read_shear_options(parameters.get('shear_flow',{})) if compute_shear_flow else None
    interface=flow_options=leakage=None
    # Check application inputs before PSD work or creation of any output.
    if flow_enabled:
        raw_flow=parameters.get('pressure_flow',{})
        unknown=set(raw_flow)-set(FlowOptions.__dataclass_fields__)
        if unknown:
            raise ValueError(f'Unknown pressure_flow parameters: {sorted(unknown)}')
        flow_options=FlowOptions(**raw_flow)
    spatial=run_leakage and parameters.get('leakage',{}).get('mode')=='spatial'
    if run_leakage and not spatial:
        raw=dict(parameters.get('leakage',{}))
        legacy=raw.pop('fluid_viscosity_pa_s',None)
        if legacy is not None:
            raise ValueError('请仅在 interface.viscosity_pa_s 中给出黏度。')
        if set(raw)!=set(LeakageParameters.__dataclass_fields__):
            raise ValueError('leakage 需要 seal_length_m、seal_width_m、inlet_pressure_pa、outlet_pressure_pa。')
        leakage=LeakageParameters(**raw)
    problem, options = problem_from_mapping({key: value for key, value in parameters.items()
                                           if key in problem_keys}, base)
    model=None;plan=None;database_extension=None
    if flow_enabled:
        model=PerssonModel(problem.psd,problem.effective_modulus_pa,options.gamma,
            n_q=options.n_q,pdf_afact=options.pdf_afact,pdf_log_step=options.pdf_log_step,
            pdf_max_relative_contact_area=options.pdf_max_relative_contact_area,
            pdf_min_relative_contact_area=options.pdf_min_relative_contact_area,
            pdf_min_u1_m=options.pdf_min_u1_m)
    if spatial:
        plan=prepare_leakage(parameters['leakage'],base,problem.effective_modulus_pa)
        requested=problem.pressures_pa
        expanded,plateau,percolation_pressure=required_pressures(plan,model,requested,options.critical_relative_contact_area)
        database_extension={'requested_maximum_pa':float(max(requested)), 'actual_maximum_pa':float(max(expanded)),
            'added_nodes':len(expanded)-len(requested),'percolation_pressure_pa':percolation_pressure,'plateau_pressure_pa':plateau}
        problem=replace(problem,pressures_pa=expanded)
    if flow_enabled:
        interface=(BallisticGas(plan.gas) if plan is not None and plan.fortran_mode==5000
                   else make_interface(parameters.get('interface',{})))
        if plan is not None:validate_liquid_estimates(plan,interface)
    shear=None
    if compute_shear_flow:
        if parameters.get('leakage',{}).get('fluid')=='gas' or (plan is not None and plan.fortran_mode==5000):
            raise NotImplementedError('Shear flow currently supports Newtonian liquid, not ballistic gas')
        if getattr(interface,'has_slip',False):
            raise NotImplementedError('The shear interpolation currently requires no-slip walls; Navier pressure flow remains available')
        if run_leakage and not spatial:
            raise ValueError('Shear-driven leakage requires the spatial leakage configuration')
        if spatial and shear_options.moving:
            if plan.coupling!='fixed_contact':
                raise NotImplementedError('Moving-wall spatial flow currently requires coupling=fixed_contact; local unloading is not coupled yet')
            if plan.method!='effective_medium':
                raise NotImplementedError('Moving-wall spatial flow requires effective_medium conductance')
            if plan.compute_equivalent_gas:
                raise ValueError('Set compute_equivalent_gas=False for moving-wall flow; those estimates assume pressure-only flow')
        shear=ShearFlowLaw.from_model(shear_options,model,options.critical_relative_contact_area)
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
                          'leakage':run_leakage,'mixed_friction':False,'shear_flow':compute_shear_flow},
            'requested_contact_output':run_contact,
            'leakage_database_extension':database_extension,
            'workflow_status':'contact_complete_flow_pending' if flow_enabled else 'complete',
        })
    if flow_enabled:
        print('计算压力流动系数；PDF 导出开关不影响内部所需的分布计算。')
        phase='flow'
        try:
            flow=run_pressure_flow(model,problem.pressures_pa,interface,flow_options,output,
                percolation=options.critical_relative_contact_area,leakage=leakage,plot=save_plots)
            if shear is not None:
                phase='shear_flow'
                sheared=run_shear_flow(shear,flow['results'],output,plot=save_plots)
                summary['shear_flow']={'file':sheared['table_file'],'diagnostics_file':'shear_diagnostics.json',
                                       'model':shear.options.model,'rough_surface':shear.options.rough_surface}
            if spatial:
                phase='leakage'
                summary['workflow_status']='flow_complete_leakage_pending'
                (output/'diagnostics.json').write_text(json.dumps(summary,indent=2,allow_nan=False)+'\n')
                table=make_transport(plan,model,interface,flow['results'],options.critical_relative_contact_area,
                                     options.critical_area_tolerance,shear_law=shear)
                leakage_result=run_spatial_leakage(plan,table,output,plot=save_plots)
                summary['leakage']={'file':'leakage_cases.csv','diagnostics_file':'leakage_diagnostics.json',
                    'case_count':len(leakage_result['cases']),'scope':leakage_result['scope']}
        except Exception as exc:
            summary['workflow_status']=phase+'_failed'
            summary[phase+'_error']=str(exc)
            (output/'diagnostics.json').write_text(json.dumps(summary,indent=2,allow_nan=False)+'\n')
            raise
        summary['workflow_status']='complete'
        summary['pressure_flow']={'file':flow['table_file'],'diagnostics_file':'flow_diagnostics.json',
            'status_counts':flow['status_counts'],'leakage_file':('leakage_cases.csv' if spatial else 'leakage_table.csv') if run_leakage else None}
        (output/'diagnostics.json').write_text(json.dumps(summary,indent=2,allow_nan=False)+'\n')
        print(f'压力流动结果已保存到 {output.resolve()}')
    return summary
