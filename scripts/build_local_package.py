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
        'docs/local_run.md','docs/shear_strength_fit.md','docs/finite_geometry.md',
        'inputs/experiment_template.xlsx','inputs/experiment_template.csv','scripts/create_experiment_workbook.py',
        'examples/data/roughness_psd.dat','examples/data/seal_pressure_profile.csv',
        'examples/local_interface.py','examples/contact_shear_model.py','examples/dry_contact.json',
        'examples/fit_demo.py','scripts/check_installation.py','scripts/build_local_package.py']
    names+=sorted(str(p.relative_to(root)).replace('\\','/') for p in (root/'src/persson_contact').rglob('*.py') if '__pycache__' not in p.parts)
    files={name:(root/name).read_bytes() for name in names}
    readme=root/'README_LOCAL.md'
    files['README.md']=(readme if readme.exists() else root/'README.md').read_bytes()
    return files


def trial_preset(files,preset):
    """Configure the same runtime for the agreed local cylinder experiment."""
    if preset is None:return files
    if preset!='silicone_cylinder':raise ValueError('Unknown local trial preset')
    files=dict(files)
    text=files['parameters.py'].decode()
    changes=[('当前为Hertz示例；功能开关在 main.py。','当前为有限圆柱曲面硅油实验工况；拟合已在 main.py 开启。'),
        ("'description': 'Hertz example.'","'description': 'Finite cylindrical rubber pad, silicone oil; local trial.'"),
        ("'young_modulus_pa': 3.0e6","'young_modulus_pa': 1.8e6"),
        ("'geometry': 'cylinder_plane'","'geometry': 'finite_cylindrical_pad'"),
        ("'line_load_n_m': 100.,","'normal_force_n': 100.,"),
        ("'domain_half_width_m': 1.2e-3","'domain_half_width_m': None"),
        ("'lower_wall_velocities_m_s': [0., 1e-4, 1e-3, 1e-2]","'lower_wall_velocities_m_s': [1e-6, 3e-6, 1e-5, 3e-5, 1e-4, 3e-4, 1e-3, 3e-3, 1e-2]"),
        ("'viscosity_pa_s': 0.1","'viscosity_pa_s': 1.0"),
        ("OUTPUT_DIRECTORY = 'results/leakage'","OUTPUT_DIRECTORY = 'results/silicone_fit_excel_01'")]
    for before,after in changes:
        if text.count(before)!=1:raise ValueError(f'Cannot apply trial configuration: {before}')
        text=text.replace(before,after)
    text=text.replace("'normal_force_n': 100.,         # N/m；可用'normal_force_n':100.替换本行（不能同时给），程序自动除以B",
        "'normal_force_n': 100.,        # 默认总力N；CSV逐行的40/100/160N覆盖此值，自动除以轴向长度B")
    entry=files['main.py'].decode().replace('RUN_LEAKAGE = True','RUN_LEAKAGE = False').replace(
        'RUN_SHEAR_STRENGTH_FIT = False','RUN_SHEAR_STRENGTH_FIT = True')
    compile(text,'parameters.py','exec');compile(entry,'main.py','exec')
    files['parameters.py']=text.encode();files['main.py']=entry.encode()
    files['README.md']=("本包已预填本次硅油实验几何、材料与黏度，并开启独立拟合。请用Excel填写 inputs/experiment.xlsx 后运行；"
        "该文件不含你的实测摩擦数据。默认100N仅用于CSV未给载荷的行。详见 docs/finite_geometry.md。\n\n".encode()+files['README.md'])
    return files


def build_package(root,output,*,preset=None):
    root=Path(root).resolve();output=Path(output).resolve()
    if output.exists():raise FileExistsError(f'不覆盖已有包：{output}')
    version=tomllib.loads((root/'pyproject.toml').read_text())['project']['version']
    runtime=re.search(r'__version__\s*=\s*"([^"]+)"',(root/'src/persson_contact/__init__.py').read_text())
    if runtime is None or runtime[1]!=version:raise ValueError('Runtime/project versions differ')
    files=trial_preset(package_files(root),preset)
    try:commit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=root,stderr=subprocess.DEVNULL).decode().strip()
    except (OSError,subprocess.CalledProcessError):commit=None
    metadata={'project':'ECL_BD','version':version,'channel':'local-trial-release-candidate',
        'source_base_commit':commit,'configuration_preset':preset,
        'status':'Numerically checked candidate; first local experimental trial pending.',
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
    parser.add_argument('--preset',choices=['silicone_cylinder'],default=None)
    args=parser.parse_args();print(json.dumps(build_package(args.root,args.output,preset=args.preset),indent=2))


if __name__=='__main__':main()
