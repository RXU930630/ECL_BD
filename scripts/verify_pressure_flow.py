#!/usr/bin/env python3
"""Compile the actual isotropic Sealing3 P(u) iteration as an isolated harness.

No original contact calculation is changed or rerun. Feed its previously
completed PDF output into the unchanged Fortran loop and the Python root solver,
then separately compare the end-to-end Python PDF + CDF-bin flow integration.
"""
import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import time

import numpy as np

import compare_fortran_reference as ref
from persson_contact.core import PerssonModel
from persson_contact.flow import bruggeman_conductivity,effective_dimension,pressure_flow
from persson_contact.interface import NewtonianNoSlip
from persson_contact.psd import read_psd

SOURCE_SHA='99d3b2e2c47a663ffe8b65bbdb7415db5eeebb434debd0ef404006600e861854'


def harness(source):
    original=Path(source).read_bytes()
    if hashlib.sha256(original).hexdigest()!=SOURCE_SHA:
        raise ValueError('Unrecognized original Fortran version')
    text=original.decode()
    routine=text[text.index('      SUBROUTINE Sealing3()'):]
    begin=routine.index('      do 2277 itt=1,NittSeal')
    end=routine.index(' 2277 continue',begin)+len(' 2277 continue')
    block=routine[begin:end]
    prefix='''program pressure_flow_reference
implicit real*8 (a-h,o-z)
real*8 Probabilityrms(100000)
open(10,file='input.dat')
read(10,*) nurange1,averagesep,relArea,percolation,uend,NittSeal,feedback1,Pumin,Palpha
if(nurange1.gt.100000) stop 2
xnDstart=2d0
xnDstop=1d0/(1d0-percolation)
Aratioperc=relArea/percolation
xnDrun=xnDstop
if(Aratioperc.lt.1d0) xnDrun=xnDstart*(1d0-Aratioperc)+xnDstop*Aratioperc
ellipsgamma=1d0
phip0=1d0
p0=0d0
Deltau=1d0
istop=1
! The extracted loop reads vecSeparationu1(istop); declaration below is scalar-array bridge.
flowsigmaeff1=averagesep**3
flowsigmaeff1y=flowsigmaeff1
AreaContact=1d0
 do jj=1,nurange1
 read(10,*) dummy,Probabilityrms(jj)
 AreaContact=AreaContact-Probabilityrms(jj)
 enddo
close(10)
open(3718,file='iteration.log')
'''
    # Add the original loop's array access, retaining every byte of that loop.
    prefix=prefix.replace('real*8 Probabilityrms(100000)', 'real*8 Probabilityrms(100000),vecSeparationu1(1),vecAria(1)')
    prefix=prefix.replace('flowsigmaeff1=averagesep**3', 'vecSeparationu1(1)=uend\nvecAria(1)=relArea\nflowsigmaeff1=averagesep**3')
    suffix='''
404 format(ES24.16,1X,I8,4(1X,ES24.16))
close(3718)
open(20,file='output.dat')
write(20,'(4ES27.17E3)') flowsigmaeff1,flowsigmaeff1y,AreaContact,sb
close(20)
end program
'''
    return prefix+block+suffix,block


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source',required=True);parser.add_argument('--matrix',required=True)
    parser.add_argument('--psd',required=True);parser.add_argument('--output',required=True)
    parser.add_argument('--compiler',default='gfortran')
    parser.add_argument('--fixture',help='Write compact flow regression references')
    args=parser.parse_args();output=Path(args.output).resolve()
    output.mkdir(parents=True,exist_ok=False)
    source,block=harness(args.source)
    (output/'flow_reference.f90').write_text(source)
    compiled=subprocess.run([args.compiler,'-O3','-ffree-line-length-none',
                             'flow_reference.f90','-o','flow_reference'],cwd=output,
                             capture_output=True,text=True,timeout=60)
    (output/'compile.log').write_text(compiled.stdout+compiled.stderr)
    if compiled.returncode:raise RuntimeError('Fortran flow harness compile failed')
    psd=read_psd(args.psd);records=[];fixture=[];skipped=[]
    for path in sorted(Path(args.matrix).glob('*/manifest.json')):
        manifest=json.loads(path.read_text());case=manifest.get('case',{})
        if manifest.get('status')!='validated' or case.get('mode')!=3:continue
        area_data=ref.read_numeric(path.parent/ref.SUMMARY_NAME,5)[0]
        _,pdf=ref.read_numeric(path.parent/ref.PDF_NAME,2,pressure_header=True,minimum_rows=40000)
        p=case['pressure_pa'];h=area_data[2];a=area_data[4]
        # High-precision scale-area output is unavailable; compute the exact
        # source right-rectangle final area instead of its E11.4 printout.
        from scipy.special import erf
        step=case['delzeta'];n=int(np.log(psd.q1/psd.q0)/step)
        q=psd.q0*np.exp(step*np.arange(1,n+1))
        i3=np.sum(psd.evaluate(q)*q**4)*step
        a=float(erf(p/(case['effective_modulus_pa']*np.sqrt(3.141592654*i3))))
        widths=ref.read_numeric(path.parent/'out.1u1.2u1min',2)
        uend=float(widths[widths[:,0]>=1e-12,0][-1])
        alpha=float(np.log(pdf[-1,0]/pdf[0,0])/(len(pdf)-1))
        weights=pdf[:,0]*pdf[:,1]*alpha
        if 1-sum(weights)>=1-1/effective_dimension(a,.42):
            # A finite-iteration source floor cannot validate a zero root;
            # verify blocked phases separately with exact analytical tests.
            skipped.append({'name':path.parent.name,'reason':'Zero-conductance phase exceeds threshold; exact blocked root checked analytically.'})
            continue
        directory=output/path.parent.name;directory.mkdir()
        header=f'{len(pdf)} {h:.17g} {a:.17g} .42 {uend:.17g} 6000 .9 1e-12 {alpha:.17g}'
        (directory/'input.dat').write_text(header+'\n'+'\n'.join(
            f'{u:.17g} {w:.17g}' for u,w in zip(pdf[:,0],weights))+'\n')
        start=time.monotonic()
        run=subprocess.run([str(output/'flow_reference')],cwd=directory,capture_output=True,text=True,timeout=180)
        if run.returncode:raise RuntimeError('Fortran pressure flow iteration failed')
        reference=np.loadtxt(directory/'output.dat');k=reference[0]
        # Same-PDF closure isolation, with physical nonnegative u^3 law.
        python,_,_=bruggeman_conductivity(pdf[:,0]**3,weights,float(1-sum(weights)),effective_dimension(a,.42))
        model=PerssonModel(psd,case['effective_modulus_pa'],.4,pdf_log_step=step)
        response=pressure_flow(model,p,NewtonianNoSlip(.1))
        record={'name':path.parent.name,**case,'reference_geometric_conductivity_m3':float(k),
            'reference_mean_gap_m':float(h),'reference_relative_contact_area':a,
            'same_pdf_root_relative_difference':float(python/k-1),
            'end_to_end_relative_difference':float(response.geometric_conductivity_m3/k-1),
            'pressure_flow_factor_relative_difference':float(response.pressure_flow_factor/(k/h**3)-1),
            'python_result':response.as_dict(),'reference_phase_mass':float(sum(weights)),
            'run_seconds':time.monotonic()-start,'reference_pdf_sha256':ref.sha256(path.parent/ref.PDF_NAME),
            'input_sha256':ref.sha256(directory/'input.dat'),'output_sha256':ref.sha256(directory/'output.dat')}
        records.append(record);print(record['name'],record['same_pdf_root_relative_difference'],record['end_to_end_relative_difference'],flush=True)
        fixture.append({key:record[key] for key in ('name','effective_modulus_pa','pressure_pa','delzeta',
            'reference_geometric_conductivity_m3','reference_mean_gap_m','reference_pdf_sha256')})
    if len(records)!=9:
        raise ValueError(f'Expected all 9 conducting reference cases, found {len(records)}')
    result={'original_source_sha256':SOURCE_SHA,'loop_sha256':hashlib.sha256(block.encode()).hexdigest(),
            'harness_sha256':ref.sha256(output/'flow_reference.f90'),
            'binary_sha256':ref.sha256(output/'flow_reference'),
            'compiler':subprocess.run([args.compiler,'--version'],capture_output=True,text=True,check=True).stdout.splitlines()[0],
            'skipped_blocked_cases':skipped,
            'scope':'Actual Sealing3 P(u) loop only, isotropic. Original full iflag008 pipeline not rerun. '
            'Numerical controls: 6000 iterations, feedback=.9; no physics equation changes.',
            'cases':records,'case_count':len(records)}
    (output/'comparison.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    if args.fixture:
        Path(args.fixture).write_text(json.dumps({'psd_sha256':ref.sha256(args.psd),
            'original_source_sha256':SOURCE_SHA,'cases':fixture},indent=2)+'\n')


if __name__=='__main__':main()
