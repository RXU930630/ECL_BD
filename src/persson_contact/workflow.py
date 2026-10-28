"""Dispatch user-facing switches without pretending unsupported physics ran."""
import hashlib
import json
from pathlib import Path

from .cli import run_contact as calculate_contact
from .config import problem_from_mapping


def run_project(parameters, *, base_directory, run_contact=True, run_leakage=False,
                run_mixed_friction=False, compute_gap_pdf=True, save_plots=True):
    switches = (run_contact, run_leakage, run_mixed_friction, compute_gap_pdf, save_plots)
    if any(type(value) is not bool for value in switches):
        raise ValueError('所有功能开关必须为 True 或 False。')
    if not any(switches[:3]):
        print('没有启用计算功能；未创建任何结果文件。')
        return None
    base = Path(base_directory).resolve()
    problem_keys = {'description', 'psd', 'material', 'pressures_pa', 'contact_pressure', 'options'}
    allowed = problem_keys | {'output_directory', 'leakage', 'mixed_friction'}
    if not isinstance(parameters, dict) or set(parameters) - allowed:
        raise ValueError('parameters.py 中存在未知参数分组。')
    for name in ('leakage', 'mixed_friction'):
        if name in parameters and not isinstance(parameters[name], dict):
            raise ValueError(f'{name} 必须为参数字典。')
    # Fail before PSD work/output creation when only a reserved application is selected.
    missing = []
    if run_leakage:
        missing.append('泄漏 RUN_LEAKAGE')
    if run_mixed_friction:
        missing.append('混合摩擦 RUN_MIXED_FRICTION')
    if missing:
        raise NotImplementedError('、'.join(missing) + ' 尚未实现；请先仅开启 RUN_CONTACT。')
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
    return calculate_contact(problem, options, output, plot=save_plots,
                             compute_gap_pdf=compute_gap_pdf, provenance={
                                 'config_format': 'Python parameters.py',
                                 'parameters_snapshot': parameters,
                                 'parameters_sha256': hashlib.sha256(snapshot.encode('utf-8')).hexdigest(),
                                 'psd_sha256': hashlib.sha256(psd.read_bytes()).hexdigest(),
                                 'material_input': parameters['material'],
                                 'functions': {'contact': True, 'leakage': False, 'mixed_friction': False},
                             })
