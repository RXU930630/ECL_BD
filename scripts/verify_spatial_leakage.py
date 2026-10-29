#!/usr/bin/env python3
"""Compare spatial solver with actual Fortran-500 pressure-unloading updates.

Feeds identical reference traction and conductivity tables into an isolated
original-source numerical kernel. This tests the spatial solve, not the original
profile reader, contact quadrature or gas-output approximations.
"""
import argparse
from copy import deepcopy
import hashlib
import json
from pathlib import Path
import subprocess
import time
import numpy as np

import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
from parameters import PARAMETERS
from persson_contact.config import problem_from_mapping
from persson_contact.core import PerssonModel
from persson_contact.interface import NewtonianNoSlip
from persson_contact.leakage import prepare_leakage,make_transport,solve_case

SHA='99d3b2e2c47a663ffe8b65bbdb7415db5eeebb434debd0ef404006600e861854'


def harness(path):
    raw=Path(path).read_bytes()
    if hashlib.sha256(raw).hexdigest()!=SHA:raise ValueError('Unrecognized original Fortran source')
    text=raw.decode();start=text.index('      SUBROUTINE SealingLargeFluidPressureSealDeformationReadInpRing()')
    text=text[start:]
    first=text.index('      s=0.0\n      do 2 ix=1,Nx')
    last=text.index('  2   continue',first)+len('  2   continue')
    resistance=text[first:last]
    first=text.index('      do 3 ix=1,Nx',last)
    stop='      if (pcontx(ix).lt.0.0) pcontx(ix)=0.0'
    last=text.index(stop,first)+len(stop)
    unloading=text[first:last]
    formula='      dotQ=pa*xL*(RMShc**3.0)/(12.0*viscosity*delx*xintegral(Nx))'
    if formula not in text:raise ValueError('Original leakage formula missing')
    code='''program reference_500
implicit real*8 (a-h,o-z)
real*8 pc(5000),sigmaeff(5000),sigmaeff1(5000),pcontx(20000),pcontx0(20000),xintegral(20000)
open(10,file='input.dat')
read(10,*) Nx,Npa,Np,xL,p00,RMShc,viscosity,paMax,old
if(Nx.gt.20000.or.Np.gt.5000) stop 2
delx=xL/(Nx+1d0)
delpa=paMax/Npa
do i=1,Np
read(10,*) pc(i),sigmaeff(i)
enddo
do ix=1,Nx
read(10,*) pcontx0(ix)
pcontx(ix)=pcontx0(ix)
enddo
close(10)
do ipa=1,Npa
pa=ipa*delpa
'''+resistance+'\n'+unloading+'\n   3 continue\nenddo\n'+formula+'''
open(20,file='output.dat')
write(20,'(ES26.17E3)') dotQ
do ix=1,Nx
pf=pcontx0(ix)-pcontx(ix)
if(pcontx(ix).le.0d0) pf=pa
write(20,'(3ES26.17E3)') ix*delx,pcontx(ix),pf
enddo
close(20)
end program
subroutine splint(xa,ya,y2,n,x,y)
implicit real*8 (a-h,o-z)
real*8 xa(n),ya(n),y2(n)
klo=1
khi=n
do while(khi-klo.gt.1)
k=(khi+klo)/2
if(xa(k).gt.x) then
khi=k
else
klo=k
endif
enddo
h=xa(khi)-xa(klo)
y=ya(klo)+(x-xa(klo))*(ya(khi)-ya(klo))/h
return
end
'''
    return code,{'resistance_update_sha256':hashlib.sha256(resistance.encode()).hexdigest(),
                 'unloading_update_sha256':hashlib.sha256(unloading.encode()).hexdigest()}


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source',required=True);parser.add_argument('--flow-data',required=True)
    parser.add_argument('--output',required=True);parser.add_argument('--compiler',default='gfortran')
    parser.add_argument('--profile',default='examples/data/seal_pressure_profile.csv')
    parser.add_argument('--fixture')
    parser.add_argument('--profile-format',default='csv',choices=['csv','fortran'])
    args=parser.parse_args();out=Path(args.output).resolve();out.mkdir(parents=True,exist_ok=False)
    code,fragments=harness(args.source);(out/'reference_500.f90').write_text(code)
    cmd=[args.compiler,'-O3','-ffree-line-length-none','reference_500.f90','-o','reference_500']
    compilation=subprocess.run(cmd,cwd=out,capture_output=True,text=True,timeout=60)
    (out/'compile.log').write_text(compilation.stdout+compilation.stderr)
    if compilation.returncode:raise RuntimeError('Fortran 500 kernel compile failed')
    problem,o=problem_from_mapping({k:v for k,v in PARAMETERS.items() if k in ('psd','material','contact_pressure','options','description')},Path.cwd())
    model=PerssonModel(problem.psd,problem.effective_modulus_pa,o.gamma)
    flow=json.loads(Path(args.flow_data).read_text());records=[];fixture_cases=[]
    geometries=[('square',{'type':'square_ring','length_m':1e-3,'reference_contact_pressure_pa':1e5},[1e4,5e4]),
                ('hertz',{'type':'oring_hertz','cross_section_diameter_m':.003,'peak_pressure_pa':1e5},[1e4,5e4]),
                ('profile',{'type':'prescribed_profile','pressure_file':str(Path(args.profile).resolve()),
                  'file_format':args.profile_format,'interpolation':'linear'},[34500.,69000.] if args.profile_format=='fortran' else [1e4,5e4])]
    for name,geom,pressures in geometries:
        raw=deepcopy(PARAMETERS['leakage']);raw['geometry']=geom;raw['database_interpolation']='linear'
        raw['solver']['tolerance']=1e-6
        plan=prepare_leakage(raw,Path.cwd(),model.effective_modulus_pa)
        table=make_transport(plan,model,NewtonianNoSlip(.1),flow['results'],.42)
        for pin in pressures:
            python=solve_case(plan,table,pin)
            refinements=[]
            for nx,npa in ((501,20000),(1001,40000),(2001,80000)):
                work=out/f'{name}_{int(pin)}_n{nx}';work.mkdir()
                length=plan.mechanics.length_m;peak=float(max(plan.mechanics.total_reference(plan.mechanics.x_m)))
                xx=plan.mechanics.x_m[0]+np.arange(1,nx+1)*length/(nx+1)
                reference=plan.mechanics.total_reference(xx)
                k=table.geometric_conductivity(table.p)
                header=f'{nx} {npa} {len(table.p)} {length:.17g} {peak:.17g} {model.psd.hrms:.17g} .1 {pin:.17g} .9'
                body='\n'.join(f'{p/peak:.17g} {value/model.psd.hrms**3:.17g}' for p,value in zip(table.p,k))
                (work/'input.dat').write_text(header+'\n'+body+'\n'+'\n'.join(f'{v:.17g}' for v in reference)+'\n')
                start=time.monotonic();run=subprocess.run([str(out/'reference_500')],cwd=work,capture_output=True,text=True,timeout=240)
                if run.returncode:raise RuntimeError(f'Fortran spatial run failed: {work}')
                lines=(work/'output.dat').read_text().splitlines();square=float(lines[0]);total=square*plan.width_m/length
                data=np.loadtxt(work/'output.dat',skiprows=1)
                px=np.interp(data[:,0]+plan.mechanics.x_m[0],python['x'],python['fluid'])
                refinements.append({'nodes':nx,'pressure_steps':npa,'feedback':.9,'fortran_flow_rate_m3_s':total,
                    'relative_rate_difference':python['flux']*plan.width_m/total-1,
                    'maximum_pressure_difference_pa':float(max(abs(px-data[:,2]))),
                    'run_seconds':time.monotonic()-start,
                    'input_sha256':hashlib.sha256((work/'input.dat').read_bytes()).hexdigest(),
                    'output_sha256':hashlib.sha256((work/'output.dat').read_bytes()).hexdigest()})
                print(name,pin,nx,'relative Q difference',refinements[-1]['relative_rate_difference'],flush=True)
            record={'geometry':name,'inlet_pressure_pa':pin,'python_flow_rate_m3_s':python['flux']*plan.width_m,
                    'refinements':refinements,'geometry_parameters':geom,
                    'python_pressure_drop_error':python['pressure_drop_error'],'python_flux_spread':python['flux_spread']}
            records.append(record)
            fixture_cases.append({'geometry':name,'inlet_pressure_pa':pin,
                'x_m':plan.mechanics.x_m.tolist(),
                'reference_contact_pressure_pa':plan.mechanics.contact_reference_pa.tolist(),
                'reference_interface_pressure_pa':plan.mechanics.interface_reference_pa.tolist(),
                'mechanics_metadata':plan.mechanics.metadata,
                'fortran_flow_rate_m3_s':refinements[-1]['fortran_flow_rate_m3_s']})
    meta={'original_source_sha256':SHA,'kernel_fragments':fragments,
        'harness_sha256':hashlib.sha256(code.encode()).hexdigest(),
        'binary_sha256':hashlib.sha256((out/'reference_500').read_bytes()).hexdigest(),
        'flow_data_sha256':hashlib.sha256(Path(args.flow_data).read_bytes()).hexdigest(),
        'scope':'Original 500 resistance and unloading arithmetic; shared pressure/conductance data. '
          'No whole original-program or gas validation. Fortran continuation and endpoint discretization retained.',
        'cases':records}
    if args.fixture:
        fixture={'source_sha256':SHA,'pressures_pa':table.p.tolist(),'mean_gaps_m':table.h.tolist(),
                 'geometric_conductivity_m3':table.k.tolist(),'width_m':plan.width_m,
                 'plateau_pressure_pa':table.plateau,'cases':fixture_cases}
        Path(args.fixture).write_text(json.dumps(fixture,indent=2,allow_nan=False)+'\n')
    (out/'comparison.json').write_text(json.dumps(meta,indent=2,allow_nan=False)+'\n')


if __name__=='__main__':main()
