#!/usr/bin/env python3
"""Compile the archive's no-slip shear-friction kernel, not its full mixed solver.

Use the FIRST uploaded September-2026 ZIP, not the later 2025 leakage source.
Only extracts algebraic interpolation lines; no proprietary source package is
republished. The resulting fixture tests the kernel, not inverse-gap closure,
cavitation, spatial solve or experimental friction.
"""
import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import tempfile
import zipfile

import numpy as np

from persson_contact.fluid_friction import shear_friction_factors


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--archive',required=True,type=Path)
    parser.add_argument('--compiler',default='gfortran')
    parser.add_argument('--output',required=True,type=Path)
    args=parser.parse_args()
    archive_sha=hashlib.sha256(args.archive.read_bytes()).hexdigest()
    expected='ae0dfddc8a5b8a8e5058d8cbe74f67fd00102bc76df42f6c734bfd41cea2b968'
    if archive_sha!=expected:raise ValueError('This is not the first uploaded archive')
    with zipfile.ZipFile(args.archive) as z:
        member=next(n for n in z.namelist() if n.endswith('/AI.ContactMechanics.14.09.2026.f90'))
        source=z.read(member)
    lines=source.decode().splitlines()
    # First no-slip interpolation used by the macroscopic mixed code, rather
    # than the later tabulated slip/shear-thinning fit.
    start=next(i for i,s in enumerate(lines) if 'Phifs=3.0*sepu/uii' in s and '!' not in s)
    assert 'Phifp=0.0' in lines[start-2] or 'Phifp=0.0' in lines[start-3]
    kernel='\n'.join(s.strip() for s in lines[start-2:start+2])
    program='''program reference
implicit real*8(a-h,o-z)
integer n,i
read(*,*) n
do i=1,n
read(*,*) sepu, avui1, rms, baruc
forfluidflow=3.0*rms*rms/(1.0+1.0)
uii=1.0/avui1
'''+kernel+'''
Phif=sepu*avui1
write(*,'(3(es25.16e3,1x))') Phif,Phifs,Phifp
end do
end program
'''
    inputs=[]
    for ratio in (.1,.5,.5000001,1.,3.,20.):
        for inverse_factor in (.8,2.):
            h=ratio*1e-6
            inputs.append([h,inverse_factor/h,1e-6,.5e-6])
    with tempfile.TemporaryDirectory() as directory:
        root=Path(directory);src=root/'reference.f90';exe=root/'reference'
        src.write_text(program)
        subprocess.run([args.compiler,'-O2','-ffree-line-length-none',str(src),'-o',str(exe)],check=True)
        text=str(len(inputs))+'\n'+'\n'.join(' '.join(f'{v:.17e}' for v in row) for row in inputs)+'\n'
        result=subprocess.run([str(exe)],input=text,text=True,capture_output=True,check=True)
    actual=np.loadtxt(result.stdout.splitlines());cases=[];max_error=0.
    for row,factors in zip(inputs,actual):
        h,inv,rms,hc=row;f=shear_friction_factors(h,inv,rms,hc)
        python=np.array([f.phi_f,f.phi_fs,f.phi_fp])
        error=float(np.max(abs(python-factors)/np.maximum(abs(factors),1e-300)))
        max_error=max(max_error,error)
        cases.append({'h':h,'inverse':inv,'rms':rms,'hc':hc,'factors':factors.tolist()})
    report={'archive_sha256':archive_sha,'source_sha256':hashlib.sha256(source).hexdigest(),
        'archive_member':member,'extracted_kernel_start_line':start-1,
        'kernel':kernel,'compiled_program_sha256':hashlib.sha256(program.encode()).hexdigest(),
        'scope':'No-slip isotropic interpolation kernel only; phi_f before original clamp/index fit. Not full Fortran mixed-lubrication equivalence.',
        'maximum_relative_difference':max_error,'cases':cases}
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
    print(f'{len(cases)} original-kernel cases; max relative difference {max_error:.3g}')
    if max_error>3e-14:raise AssertionError('Fortran kernel mismatch')


if __name__=='__main__':main()
