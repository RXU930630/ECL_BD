#!/usr/bin/env python3
"""Verify K2 against the untouched Sealing3gas loop on an original dense PDF.

This isolates original Fortran arithmetic. PDF/contact inputs must be produced
independently by the original program; no whole-program gas validation is claimed.
"""
import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import sys

import numpy as np
from scipy.special import erf

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
import compare_fortran_reference as ref
from verify_pressure_flow import harness as liquid_harness,SOURCE_SHA
from persson_contact.ballistic import ballistic_pressure_flow
from persson_contact.core import PerssonModel
from persson_contact.flow import bruggeman_conductivity,effective_dimension
from persson_contact.gas import BallisticGas,GasProperties
from persson_contact.psd import read_psd


def harness(source):
    raw=Path(source).read_bytes()
    if hashlib.sha256(raw).hexdigest()!=SOURCE_SHA:
        raise ValueError('Unrecognized original Fortran source')
    text=raw.decode();start=text.index('      SUBROUTINE Sealing3gas()')
    stop=text.index('      SUBROUTINE Sealing3()',start)
    routine=text[start:stop];begin=routine.index('      do 2277 itt=1,NittSeal')
    end=routine.index(' 2277 continue',begin)+len(' 2277 continue')
    block=routine[begin:end]
    code,liquid_block=liquid_harness(source)
    code=code.replace(liquid_block,block).replace('averagesep**3','averagesep**2')
    return code,block


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source',required=True);parser.add_argument('--pdf-directories',nargs='+',required=True)
    parser.add_argument('--psd',required=True);parser.add_argument('--pressure-pa',type=float,required=True)
    parser.add_argument('--effective-modulus-pa',type=float,required=True)
    parser.add_argument('--log-step',type=float,default=.01)
    parser.add_argument('--output',required=True);parser.add_argument('--compiler',default='gfortran')
    parser.add_argument('--fixture')
    args=parser.parse_args();out=Path(args.output).resolve();out.mkdir(parents=True,exist_ok=False)
    code,block=harness(args.source);(out/'ballistic_reference.f90').write_text(code)
    cmd=[args.compiler,'-O3','-ffree-line-length-none','ballistic_reference.f90','-o','ballistic_reference']
    compilation=subprocess.run(cmd,cwd=out,capture_output=True,text=True,timeout=60)
    (out/'compile.log').write_text(compilation.stdout+compilation.stderr)
    if compilation.returncode:raise RuntimeError('Original ballistic loop compilation failed')
    psd=read_psd(args.psd);p=args.pressure_pa;step=args.log_step
    q=psd.q0*np.exp(step*np.arange(1,int(np.log(psd.q1/psd.q0)/step)+1))
    area=float(erf(p/(args.effective_modulus_pa*np.sqrt(3.141592654*np.sum(psd.evaluate(q)*q**4)*step))))
    model=PerssonModel(psd,args.effective_modulus_pa,.4,pdf_log_step=step)
    independent=ballistic_pressure_flow(model,p,BallisticGas(GasProperties()))
    records=[]
    for index,name in enumerate(args.pdf_directories):
        directory=Path(name)
        _,pdf=ref.read_numeric(directory/ref.PDF_NAME,2,pressure_header=True,minimum_rows=1000)
        summary=ref.read_numeric(directory/ref.SUMMARY_NAME,5)[0];h=summary[2]
        widths=ref.read_numeric(directory/'out.1u1.2u1min',2)
        uend=float(widths[widths[:,0]>=1e-12,0][-1])
        alpha=float(np.log(pdf[-1,0]/pdf[0,0])/(len(pdf)-1))
        weights=pdf[:,0]*pdf[:,1]*alpha
        same_pdf,_,_=bruggeman_conductivity(pdf[:,0]**2,weights,float(1-sum(weights)),effective_dimension(area,.42))
        refinement=[]
        for iterations in (12000,24000):
            work=out/f'pdf{index}_n{len(pdf)}_i{iterations}';work.mkdir()
            header=f'{len(pdf)} {h:.17g} {area:.17g} .42 {uend:.17g} {iterations} .9 1e-12 {alpha:.17g}'
            (work/'input.dat').write_text(header+'\n'+'\n'.join(f'{u:.17g} {w:.17g}' for u,w in zip(pdf[:,0],weights))+'\n')
            run=subprocess.run([str(out/'ballistic_reference')],cwd=work,capture_output=True,text=True,timeout=180)
            if run.returncode:raise RuntimeError('Original ballistic loop failed')
            result=np.loadtxt(work/'output.dat');k=float(result[0])
            refinement.append({'iterations':iterations,'fortran_K2_m2':k,
                'same_pdf_nonnegative_u2_root_relative_difference':same_pdf/k-1,
                'input_sha256':ref.sha256(work/'input.dat'),'output_sha256':ref.sha256(work/'output.dat')})
        record={'pdf_points':len(pdf),'pressure_pa':p,'effective_modulus_pa':args.effective_modulus_pa,
            'log_step':step,'original_pdf_sha256':ref.sha256(directory/ref.PDF_NAME),
            'fortran_K2_m2':refinement[-1]['fortran_K2_m2'],
            'python_independent_K2_m2':independent.geometric_conductivity_m2,
            'independent_relative_difference':independent.geometric_conductivity_m2/refinement[-1]['fortran_K2_m2']-1,
            'refinements':refinement}
        records.append(record);print(json.dumps(record),flush=True)
    summary={'source_sha256':SOURCE_SHA,'loop_sha256':hashlib.sha256(block.encode()).hexdigest(),
        'driver_sha256':ref.sha256(out/'ballistic_reference.f90'),
        'binary_sha256':ref.sha256(out/'ballistic_reference'),
        'compiler':subprocess.run([args.compiler,'--version'],capture_output=True,text=True,check=True).stdout.splitlines()[0],
        'scope':'Untouched Sealing3gas P(u) iteration with independent original PDFs. '
                'Original u^2-u_end^2 subtraction retained; Python local u^2 remains nonnegative.',
        'psd_sha256':ref.sha256(args.psd),'cases':records}
    (out/'comparison.json').write_text(json.dumps(summary,indent=2)+'\n')
    if args.fixture:Path(args.fixture).write_text(json.dumps(summary,indent=2)+'\n')


if __name__=='__main__':main()
