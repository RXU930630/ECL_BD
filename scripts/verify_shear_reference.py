#!/usr/bin/env python3
"""Check the actual Sealing3 shear interpolation; no whole-program claim."""
import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import sys

import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
from persson_contact.shear import shear_length

SOURCE_SHA='99d3b2e2c47a663ffe8b65bbdb7415db5eeebb434debd0ef404006600e861854'


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source',required=True)
    parser.add_argument('--compiler',default='gfortran')
    parser.add_argument('--output',required=True)
    parser.add_argument('--fixture')
    args=parser.parse_args()
    data=Path(args.source).read_bytes()
    if hashlib.sha256(data).hexdigest()!=SOURCE_SHA:raise ValueError('Unrecognized original source')
    text=data.decode();routine=text.index('      SUBROUTINE Sealing3()')
    start=text.index('      phis=RMShc/ubar+',routine)
    end=text.index('      ucheck=',start)
    block=text[start:end]
    out=Path(args.output).resolve();out.mkdir(parents=True,exist_ok=False)
    driver="""program shear_reference
implicit real*8(a-h,o-z)
open(10,file='input.dat')
open(20,file='output.dat')
read(10,*) n
ellipsgamma=1d0
do i=1,n
read(10,*) ubar,RMShc,baruc11
diffsep=ubar-baruc11
"""+block+"""write(20,'(ES26.17E3)') phis
enddo
close(10)
close(20)
end program
"""
    (out/'reference.f90').write_text(driver)
    run=subprocess.run([args.compiler,'-O2','-ffree-line-length-none','reference.f90','-o','reference'],cwd=out,capture_output=True,text=True)
    if run.returncode:raise RuntimeError(run.stderr)
    cases=[]
    for rms,hc in ((1.4e-6,2.5e-8),(2e-6,5e-7)):
        for h in (hc/2,hc,hc+1e-9,rms,3*rms,100*rms):
            cases.append({'mean_gap_m':h,'rms_height_m':rms,'percolation_gap_m':hc})
    cases.append({'mean_gap_m':5e-7+5e-11,'rms_height_m':2e-6,'percolation_gap_m':5e-7})
    (out/'input.dat').write_text(str(len(cases))+'\n'+'\n'.join(
        f"{c['mean_gap_m']:.17g} {c['rms_height_m']:.17g} {c['percolation_gap_m']:.17g}" for c in cases)+'\n')
    subprocess.run([str(out/'reference')],cwd=out,check=True,capture_output=True)
    values=np.loadtxt(out/'output.dat')
    for case,value in zip(cases,values):
        h,rms,hc=(case[k] for k in ('mean_gap_m','rms_height_m','percolation_gap_m'))
        python=float(shear_length(h,rms,hc))/rms
        case.update(fortran_shear_flow_factor=float(value),python_shear_flow_factor=python,
                    relative_difference=python/float(value)-1,
                    original_cutoff_only=0<h-hc<1e-10)
    regular=[abs(c['relative_difference']) for c in cases if not c['original_cutoff_only']]
    assert max(regular)<1e-13
    summary={'source_sha256':SOURCE_SHA,'kernel_sha256':hashlib.sha256(block.encode()).hexdigest(),
        'driver_sha256':hashlib.sha256(driver.encode()).hexdigest(),
        'compiler':subprocess.check_output([args.compiler,'--version'],text=True).splitlines()[0],
        'scope':'Unmodified Sealing3 phi_s kernel, isotropic and no slip; supplied mean/RMS/percolation gaps. Not a full mixed-lubrication run.',
        'branch_difference':'Python uses max(h-hc,0); original uses a fixed positive 1e-10 m switch. One explicit near-cutoff case records the difference.',
        'cases':cases}
    (out/'comparison.json').write_text(json.dumps(summary,indent=2)+'\n')
    if args.fixture:Path(args.fixture).write_text(json.dumps(summary,indent=2)+'\n')
    print(json.dumps({'cases':len(cases),'max_regular_relative_difference':max(regular),
                      'near_cutoff_relative_difference':cases[-1]['relative_difference']}))


if __name__=='__main__':main()
