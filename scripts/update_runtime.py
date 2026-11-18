#!/usr/bin/env python3
"""Update runtime modules while preserving the target's settings, inputs and results."""
import argparse
import hashlib
import json
from pathlib import Path
import shutil


def update_runtime(source_project,target_project):
    source=Path(source_project).resolve();target=Path(target_project).resolve()
    if source==target:raise ValueError('Choose your previous project directory, not this new package.')
    runtime=target/'src/persson_contact';new=source/'src/persson_contact'
    if not (runtime/'__init__.py').is_file() or not all((target/n).is_file() for n in ('main.py','parameters.py')):
        raise ValueError('Target must be a complete ECL_BD project with main.py, parameters.py and src.')
    files=sorted(new.glob('*.py'))
    if not files or not (new/'__init__.py').is_file():raise ValueError('Source runtime is missing.')
    protected={n:hashlib.sha256((target/n).read_bytes()).hexdigest() for n in ('main.py','parameters.py')}
    backup=target/'runtime_backup'
    if backup.exists():raise FileExistsError('runtime_backup already exists. Preserve or rename it before another update.')
    backup.mkdir()
    previous=sorted(runtime.glob('*.py'))
    for path in previous:shutil.copyfile(path,backup/path.name)
    try:
        for path in files:shutil.copyfile(path,runtime/path.name)
        if any(hashlib.sha256((target/n).read_bytes()).hexdigest()!=digest for n,digest in protected.items()):
            raise RuntimeError('Configuration changed unexpectedly.')
    except Exception:
        old_names={p.name for p in previous}
        for path in files:
            if path.name not in old_names:(runtime/path.name).unlink(missing_ok=True)
        for path in previous:shutil.copyfile(backup/path.name,runtime/path.name)
        raise
    record={'updated_modules':[p.name for p in files],'preserved_configuration_sha256':protected,
            'scope':'Only src/persson_contact/*.py updated; main.py, parameters.py, inputs and results preserved.'}
    (backup/'update_record.json').write_text(json.dumps(record,indent=2,ensure_ascii=False)+'\n',encoding='utf-8')
    return record


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('target',type=Path,help='Previous project directory containing your own settings')
    args=parser.parse_args()
    try:record=update_runtime(Path(__file__).resolve().parents[1],args.target)
    except (ValueError,FileExistsError) as error:parser.error(str(error))
    print(f'Updated {len(record["updated_modules"])} runtime modules. Settings, inputs and results preserved.')
    print('Old runtime modules are in runtime_backup. Run main.py in your previous project directory.')
    return 0


if __name__=='__main__':raise SystemExit(main())
