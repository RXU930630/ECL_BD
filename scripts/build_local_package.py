#!/usr/bin/env python3
"""Create a deterministic runtime ZIP using an explicit allowlist.

No historical results, validation fixtures, Git metadata, environment, caches,
user experiments, source attachments, binaries, or installed dependencies.
"""
import argparse
import hashlib
import json
from pathlib import Path
import re
import subprocess
import tomllib
import zipfile


def package_files(root):
    root=Path(root)
    names=['main.py','parameters.py','pyproject.toml','requirements.txt',
        'docs/local_run.md','docs/shear_strength_fit.md','inputs/experiment_template.csv',
        'examples/data/roughness_psd.dat','examples/data/seal_pressure_profile.csv',
        'examples/local_interface.py','examples/contact_shear_model.py','examples/dry_contact.json',
        'examples/fit_demo.py','scripts/check_installation.py','scripts/build_local_package.py']
    names+=sorted(str(p.relative_to(root)).replace('\\','/') for p in (root/'src/persson_contact').rglob('*.py') if '__pycache__' not in p.parts)
    files={name:(root/name).read_bytes() for name in names}
    readme=root/'README_LOCAL.md'
    files['README.md']=(readme if readme.exists() else root/'README.md').read_bytes()
    return files


def build_package(root,output):
    root=Path(root).resolve();output=Path(output).resolve()
    if output.exists():raise FileExistsError(f'不覆盖已有包：{output}')
    version=tomllib.loads((root/'pyproject.toml').read_text())['project']['version']
    runtime=re.search(r'__version__\s*=\s*"([^"]+)"',(root/'src/persson_contact/__init__.py').read_text())
    if runtime is None or runtime[1]!=version:raise ValueError('Runtime/project versions differ')
    files=package_files(root)
    try:commit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=root,stderr=subprocess.DEVNULL).decode().strip()
    except (OSError,subprocess.CalledProcessError):commit=None
    metadata={'project':'ECL_BD','version':version,'channel':'local-trial-release-candidate',
        'source_base_commit':commit,'status':'Numerically checked candidate; first local experimental trial pending.',
        'scope':'Runtime source and input templates only; no historical calculation results or experimental data.',
        'files':{name:{'sha256':hashlib.sha256(data).hexdigest(),'bytes':len(data)} for name,data in sorted(files.items())}}
    files['PACKAGE_MANIFEST.json']=(json.dumps(metadata,indent=2,ensure_ascii=False)+'\n').encode()
    prefix=f'ECL_BD-{version}/'
    output.parent.mkdir(parents=True,exist_ok=True)
    with zipfile.ZipFile(output,'w',compression=zipfile.ZIP_DEFLATED,compresslevel=9) as archive:
        for name,data in sorted(files.items()):
            info=zipfile.ZipInfo(prefix+name,date_time=(2026,1,1,0,0,0));info.compress_type=zipfile.ZIP_DEFLATED
            info.external_attr=0o100644<<16;info.create_system=3
            archive.writestr(info,data,compress_type=zipfile.ZIP_DEFLATED,compresslevel=9)
    with zipfile.ZipFile(output) as archive:
        if archive.testzip() is not None:raise ValueError('ZIP integrity check failed')
        for name,entry in metadata['files'].items():
            if hashlib.sha256(archive.read(prefix+name)).hexdigest()!=entry['sha256']:raise ValueError('Package hash mismatch')
    checksum=hashlib.sha256(output.read_bytes()).hexdigest()
    output.with_suffix(output.suffix+'.sha256').write_text(checksum+'  '+output.name+'\n')
    return {'path':str(output),'sha256':checksum,'version':version,'files':len(files),'bytes':output.stat().st_size}


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root',type=Path,default=Path(__file__).resolve().parents[1])
    parser.add_argument('--output',required=True,type=Path)
    args=parser.parse_args();print(json.dumps(build_package(args.root,args.output),indent=2))


if __name__=='__main__':main()
