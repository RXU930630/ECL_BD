#!/usr/bin/env python3
"""Finish a recorded failed fit from complete saved forward states, in its original directory."""
import argparse
import ast
from copy import deepcopy
import json
from pathlib import Path
import runpy
import shutil
import sys
import tempfile


def recover_fit(project):
    project=Path(project).resolve()
    parameters=deepcopy(runpy.run_path(str(project/'parameters.py'))['PARAMETERS'])
    out=Path(parameters['output_directory'])
    out=out if out.is_absolute() else project/out
    failure=out/'fit_diagnostics.json'
    status=json.loads(failure.read_text(encoding='utf-8-sig')).get('workflow_status')
    if status not in ('forward_failed','fit_failed'):
        raise ValueError('Recovery requires a recorded failed fit; completed results are never overwritten.')
    directories=[]
    for directory in sorted((out/'forward').glob('load_*')):
        diagnostic=directory/'mixed_diagnostics.json'
        try:diag=json.loads(diagnostic.read_text(encoding='utf-8-sig'))
        except (OSError,ValueError):continue
        if diag.get('workflow_status')=='complete' and (directory/'mixed_cases.csv').is_file():
            directories.append(str(directory.resolve()))
    if not directories:raise ValueError('No complete saved forward states found; recalculate the incomplete cases after updating.')
    save_plots=True
    for node in ast.parse((project/'main.py').read_text(encoding='utf-8-sig')).body:
        if isinstance(node,ast.Assign) and any(isinstance(t,ast.Name) and t.id=='SAVE_PLOTS' for t in node.targets):
            if isinstance(node.value,ast.Constant) and type(node.value.value) is bool:save_plots=node.value.value
    backup=project/'runtime_backup'/'fit_before_recovery'
    if backup.exists():raise FileExistsError('fit_before_recovery already exists; preserve or rename the backup before retrying.')
    sys.path.insert(0,str(project/'src'))
    from persson_contact.workflow import run_project
    parameters['shear_strength_fit'].update(source='saved_forward',source_parameters={
        'saved_forward':{'directories':directories}})
    # Compute and validate all output files before touching the failed run.
    with tempfile.TemporaryDirectory(prefix='.fit-recovery-',dir=out) as temporary:
        target=Path(temporary)/'fit'
        parameters['output_directory']=str(target)
        result=run_project(parameters,base_directory=project,run_contact=False,
                           run_shear_strength_fit=True,save_plots=save_plots)
        result['parameters_snapshot']['output_directory']=str(out)
        (target/'fit_diagnostics.json').write_text(
            json.dumps(result,indent=2,ensure_ascii=False,allow_nan=False)+'\n',encoding='utf-8')
        names=result['outputs'];backup.mkdir(parents=True)
        existing={name for name in names if (out/name).is_file()}
        for name in existing:shutil.copyfile(out/name,backup/name)
        published=[]
        try:
            for name in names:
                (target/name).replace(out/name);published.append(name)
        except Exception:
            for name in published:
                if name in existing:shutil.copyfile(backup/name,out/name)
                else:(out/name).unlink(missing_ok=True)
            raise
    return result


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('project',type=Path,help='Updated previous project directory with your original settings')
    args=parser.parse_args()
    try:result=recover_fit(args.project)
    except (ValueError,OSError) as error:parser.error(str(error))
    print('Fit recovered in the original output directory; settings and forward states preserved.')
    print(f'tau = {result["fit"]["fitted_shear_strength_pa"]:.8g} Pa; inspect fit_results.xlsx and diagnostics.')
    return 0


if __name__=='__main__':raise SystemExit(main())
