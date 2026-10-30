#!/usr/bin/env python3
"""Compare the spatial solver with the original 5000 resistance/unloading loop.

The independent K2 generation is tested separately. Both spatial solvers receive
the same K2 table, reference traction and controlled gas constants here.
"""
import argparse
from copy import deepcopy
import hashlib
import json
from pathlib import Path
import subprocess
import sys
import time

import numpy as np

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
from persson_contact.psd import read_psd
from persson_contact.core import PerssonModel
from persson_contact.gas import BallisticGas
from persson_contact.leakage import prepare_leakage,make_transport,solve_case
from verify_pressure_flow import SOURCE_SHA


def harness(path):
    raw=Path(path).read_bytes()
    if hashlib.sha256(raw).hexdigest()!=SOURCE_SHA:raise ValueError('Unrecognized original source')
    text=raw.decode();start=text.index('      SUBROUTINE SealingLargeGasPressureSealDeformationReadInpRing()')
    stop=text.index('      SUBROUTINE SealingLargeFluidPressureSealDeformationReadInpRing()',start)
    routine=text[start:stop]
    first=routine.index('      s=0.0\n      do 2 ix=1,Nx')
    last=routine.index('  2   continue',first)+len('  2   continue');resistance=routine[first:last]
    first=routine.index('      do 3 ix=1,Nx',last)
    end='      if (pcontx(ix).lt.0.0) pcontx(ix)=0.0'
    last=routine.index(end,first)+len(end);unloading=routine[first:last]
    formula='      dotQ=pa*xL*barv*(RMShc**2.0)/(2.0*ETemp*delx*xintegral(Nx))'
    assert formula in routine
    volume='      xV=dotQ*Etemp/patm';assert volume in routine
    code='''program reference_5000
implicit real*8 (a-h,o-z)
real*8 pc(5000),sigmaeff(5000),sigmaeff1(5000),pcontx(20000),pcontx0(20000),xintegral(20000)
open(10,file='input.dat')
read(10,*) Nx,Npa,Np,xL,p00,RMShc,barv,ETemp,patm,paMax,old
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
'''+resistance+'\n'+unloading+'\n   3 continue\nenddo\n'+formula+'\n'+volume+'''
open(20,file='output.dat')
write(20,'(2ES26.17E3)') dotQ,xV
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
    return code,{'resistance_loop_sha256':hashlib.sha256(resistance.encode()).hexdigest(),
                 'unloading_loop_sha256':hashlib.sha256(unloading.encode()).hexdigest()}


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source',required=True);parser.add_argument('--flow-data',required=True)
    parser.add_argument('--leakage-diagnostics',required=True);parser.add_argument('--psd',required=True)
    parser.add_argument('--effective-modulus-pa',type=float,required=True)
    parser.add_argument('--output',required=True);parser.add_argument('--compiler',default='gfortran')
    parser.add_argument('--fixture')
    parser.add_argument('--finest-nodes',type=int,choices=(2001,4001,8001),default=2001,
                        help='Refine the sensitive high-pressure case beyond the basic three levels')
    args=parser.parse_args();out=Path(args.output).resolve();out.mkdir(parents=True,exist_ok=False)
    code,hashes=harness(args.source);(out/'reference_5000.f90').write_text(code)
    compilation=subprocess.run([args.compiler,'-O3','-ffree-line-length-none','reference_5000.f90','-o','reference_5000'],cwd=out,capture_output=True,text=True,timeout=60)
    (out/'compile.log').write_text(compilation.stdout+compilation.stderr)
    if compilation.returncode:raise RuntimeError('Original 5000 kernel compilation failed')
    raw=json.loads(Path(args.leakage_diagnostics).read_text())['parameters']
    model=PerssonModel(read_psd(args.psd),args.effective_modulus_pa,.4)
    plan=prepare_leakage(raw,Path.cwd(),model.effective_modulus_pa)
    rows=json.loads(Path(args.flow_data).read_text());table=make_transport(plan,model,BallisticGas(plan.gas),rows,.42)
    records=[]
    for pin in (138000.,414370.):
        python=solve_case(plan,table,pin);length=plan.mechanics.length_m;peak=plan.mechanics.maximum_reference_pressure()
        refinement=[]
        levels=[(501,20000),(1001,40000),(2001,80000)]
        if pin==414370.:
            levels.extend((n,(n-1)*40) for n in (4001,8001) if n<=args.finest_nodes)
        for nx,npa in levels:
            work=out/f'p{int(pin)}_n{nx}';work.mkdir()
            x=plan.mechanics.x_m[0]+np.arange(1,nx+1)*length/(nx+1)
            pp=np.unique(np.r_[table.p,table.plateau]) if table.plateau is not None else table.p
            header=f'{nx} {npa} {len(pp)} {length:.17g} {peak:.17g} {model.psd.hrms:.17g} {plan.gas.mean_molecular_speed_m_s:.17g} {plan.gas.thermal_energy_j:.17g} {plan.gas.reference_pressure_pa:.17g} {pin:.17g} .9'
            k=table.geometric_conductivity(pp)
            body='\n'.join(f'{p/peak:.17g} {ki/model.psd.hrms**2:.17g}' for p,ki in zip(pp,k))
            (work/'input.dat').write_text(header+'\n'+body+'\n'+'\n'.join(f'{v:.17g}' for v in plan.mechanics.total_reference(x))+'\n')
            t=time.monotonic();run=subprocess.run([str(out/'reference_5000')],cwd=work,capture_output=True,text=True,timeout=180)
            if run.returncode:raise RuntimeError('Original 5000 kernel failed')
            first=np.loadtxt(work/'output.dat',max_rows=1);data=np.loadtxt(work/'output.dat',skiprows=1)
            q=float(first[1]*plan.width_m/length)
            expected=plan.gas.reference_volume(python['flux']*plan.width_m)
            pf=np.interp(data[:,0]+plan.mechanics.x_m[0],python['x'],python['fluid'])
            refinement.append({'Nx':nx,'Npa':npa,'fortran_reference_volume_m3_s':q,
                'fortran_molecular_flow_per_s':float(first[0]*plan.width_m/length),
                'relative_volume_difference':expected/q-1,
                'maximum_fluid_pressure_error_pa':float(max(abs(pf-data[:,2]))),
                'run_seconds':time.monotonic()-t,
                'input_sha256':hashlib.sha256((work/'input.dat').read_bytes()).hexdigest(),
                'output_sha256':hashlib.sha256((work/'output.dat').read_bytes()).hexdigest()})
            print(pin,nx,json.dumps(refinement[-1]),flush=True)
        records.append({'inlet_pressure_pa':pin,'python_reference_volume_m3_s':expected,'refinements':refinement})
    summary={'source_sha256':SOURCE_SHA,'kernel_hashes':hashes,
        'driver_sha256':hashlib.sha256(code.encode()).hexdigest(),
        'binary_sha256':hashlib.sha256((out/'reference_5000').read_bytes()).hexdigest(),
        'compiler':subprocess.run([args.compiler,'--version'],capture_output=True,text=True,check=True).stdout.splitlines()[0],
        'scope':'Original 5000 resistance, pressure unloading and gas leakage arithmetic; shared K2 data and reference traction. Gas constants supplied identically to isolate arithmetic; not a complete original-program gas run.',
        'cases':records}
    (out/'comparison.json').write_text(json.dumps(summary,indent=2)+'\n')
    if args.fixture:
        summary['shared_table']={'pressures_pa':table.p.tolist(),'mean_gaps_m':table.h.tolist(),
            'K2_m2':table.k.tolist(),'upper_K2_m2':table.ku.tolist(),'plateau_pressure_pa':table.plateau}
        summary['geometry']={'x_m':plan.mechanics.x_m.tolist(),
            'contact_reference_pa':plan.mechanics.contact_reference_pa.tolist(),
            'metadata':plan.mechanics.metadata,'width_m':plan.width_m}
        summary['gas_parameters']=raw['gas']
        Path(args.fixture).write_text(json.dumps(summary,indent=2,allow_nan=False)+'\n')


if __name__=='__main__':main()
