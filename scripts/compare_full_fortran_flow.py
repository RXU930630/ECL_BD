"""Compare a complete full-program Fortran factor table with the current Python model."""
from pathlib import Path
import argparse,shutil,hashlib
import csv,json,sys
import numpy as np
PROJECT=Path(__file__).resolve().parents[1]
sys.path[:0]=[str(PROJECT/'src'),str(PROJECT/'scripts')]
from persson_contact.core import PerssonModel
from persson_contact.psd import read_psd
from persson_contact.interface import NewtonianNoSlip
from persson_contact.flow import pressure_flow
from persson_contact.shear import ShearOptions,ShearFlowLaw,shear_length
from compare_fortran_reference import fortran_float,first_tokens

parser=argparse.ArgumentParser(description=__doc__)
parser.add_argument('--reference-dir',required=True)
parser.add_argument('--output',required=True)
parser.add_argument('--allow-diagnostic-print-failure',action='store_true')
args=parser.parse_args()
root=Path(args.reference_dir).resolve();dest=Path(args.output).resolve()
if dest.exists():raise FileExistsError('Refusing to overwrite comparison output')
name=root.name
m=json.loads((root/'manifest.json').read_text())
assert m['run']['returncode']==0 or (args.allow_diagnostic_print_failure and m['run']['returncode']==2)
assert m['run']['timed_out'] is False
assert m['source']['original_sha256']=='99d3b2e2c47a663ffe8b65bbdb7415db5eeebb434debd0ef404006600e861854'
assert hashlib.sha256((root/m['source']['file']).read_bytes()).hexdigest()==m['source']['compiled_sha256']
for filename,expected in m['used_input_hashes'].items():
 assert hashlib.sha256((root/filename).read_bytes()).hexdigest()==expected,filename
if m['run']['returncode']==2:
 assert 'Expected REAL for item 1 in formatted transfer, got INTEGER' in (root/'run.log').read_text()
f=next(root.glob('out.FLUID.FLOWFACTORS.*'))
lines=[line.strip() for line in f.read_text().splitlines() if line.strip()]
count=int(lines[0]);full=np.array([[fortran_float(v) for v in line.split()] for line in lines[1:]])
assert full.shape==(count,8) and np.all(np.isfinite(full))
def read(path,n):
 out=[]
 for line in path.read_text().splitlines():
  if len(line.split())==n:out.append([fortran_float(v) for v in line.split()])
 return np.asarray(out)
steps=read(next(root.glob('out.1ipressurestep*')),3)
n=len(steps);assert n>0
native=full[-n:];np.testing.assert_allclose(native[:,0],steps[:,2],rtol=6e-6)
footer=[fortran_float(line.strip()) for line in (root/'WORK.1logpressure.2averageSeparation.3PressureFlowFactor.4Area').read_text().splitlines() if len(line.split())==1]
assert len(footer)==2
hc_ft=footer[0]
psd=read_psd(root/'IN_ROADTYPE.q0.qr.q1.h0.H.1logWavevector.2logPowerSpectrum')
mat=first_tokens(root/'IN.sigma.gammabar.xnu.ElasticModul.yieldstress')
num=first_tokens(root/'IN.mathematical')
eta=first_tokens(root/'IN.FLUID.viscosity.and.TrippNumber')[0]
model=PerssonModel(psd,mat[3]/(1-mat[2]**2),num[8],pdf_log_step=num[0],pdf_afact=num[14])
law=ShearFlowLaw.from_model(ShearOptions(),model,.42)
q=psd.q0*np.exp(num[0]*np.arange(1,int(np.log(psd.q1/psd.q0)/num[0])+1))
rms_ft=float(np.sqrt(2*float(np.float32(3.141592654))*num[0]*np.sum(psd.evaluate(q)*q*q)))
raw_path=root/'WORK.SEALING.using.P.1p0.2Normalized.averageSeparation.3Normalized.Sigma.4PressureFlowFactor'
raw=read(raw_path,4);assert len(raw)==n
rows=[]
dest.mkdir(parents=True,exist_ok=False)
for i,(state,step) in enumerate(zip(native,steps)):
 p=float(step[2]);h=float(model.mean_gap(p));hf=float(step[1])*rms_ft
 shear=float(shear_length(h,psd.hrms,law.percolation_gap_m)/psd.hrms)
 shared=float(shear_length(hf,rms_ft,hc_ft)/rms_ft)
 r={'pressure_pa':p,'fortran_gap_m':hf,'python_gap_m':h,
    'fortran_phi_p':state[2],'fortran_raw_phi_p':raw[i,3],
    'fortran_phi_s':state[3],'python_phi_s':shear,
    'shared_fortran_state_python_phi_s':shared,
    'fortran_relative_area':state[5],'python_relative_area':float(model.area(p)),
    'fortran_gap_over_rms':float(step[1]),
    'fortran_K_final_m3':state[2]*hf**3,'fortran_K_raw_m3':raw[i,3]*hf**3,
    'fortran_branch':'asymptotic' if step[1]>=num[1] else 'direct_pdf'}
 try:
  flow=pressure_flow(model,p,NewtonianNoSlip(eta))
  r.update(python_phi_p=flow.pressure_flow_factor,python_K_m3=flow.geometric_conductivity_m3,
           python_pressure_status=flow.status,python_error='')
 except ValueError as e:r.update(python_phi_p=None,python_K_m3=None,python_pressure_status='failed',python_error=str(e))
 rows.append(r)
 print('PYTHON',name,i+1,'/',n,p,r['python_pressure_status'],flush=True)
 with (dest/'comparison.csv').open('w',newline='') as out:
  w=csv.DictWriter(out,fieldnames=list(rows[0]),lineterminator='\n');w.writeheader();w.writerows(rows)
 (dest/'rows.json').write_text(json.dumps(rows,indent=2)+'\n')
def metric(key_py,key_ft,threshold=0.):
 good=[r for r in rows if r.get(key_py) is not None and abs(r[key_ft])>threshold]
 if not good:return None
 rel=np.array([r[key_py]/r[key_ft]-1 for r in good])
 return {'count':len(good),'max_abs_relative':float(np.max(abs(rel))),
         'max_abs_difference':float(max(abs(r[key_py]-r[key_ft]) for r in good)),
         'worst_pressure_pa':good[int(np.argmax(abs(rel)))]['pressure_pa']}
summary={'archive':m['archive'],'source':m['source'],'run':m['run'],'compiler':m['compiler'],
 'controls':{'material':mat,'numerical':num,'fluid_viscosity_pa_s':eta},
 'counts':{'native':n,'extended':count-n,'total':count},
 'fortran_rms_reconstructed_m':rms_ft,'python_rms_m':psd.hrms,
 'fortran_percolation_gap_m':hc_ft,'python_percolation_gap_m':law.percolation_gap_m,
 'metrics':{'mean_gap':metric('python_gap_m','fortran_gap_m'),
            'shear_end_to_end':metric('python_phi_s','fortran_phi_s'),
            'shear_shared_state':metric('shared_fortran_state_python_phi_s','fortran_phi_s'),
            'pressure_end_to_end':metric('python_phi_p','fortran_phi_p',1e-6),
            'pressure_phi_above_001':metric('python_phi_p','fortran_phi_p',.01),
            'pressure_vs_raw':metric('python_phi_p','fortran_raw_phi_p',1e-6)},
 'conducting_K_relative_errors':[{'pressure_pa':r['pressure_pa'],'branch':r['fortran_branch'],
    'vs_final':r['python_K_m3']/r['fortran_K_final_m3']-1,
    'vs_raw':r['python_K_m3']/r['fortran_K_raw_m3']-1} for r in rows
     if r['python_K_m3'] is not None and r['fortran_phi_p']>.01],
 'python_failures':[r for r in rows if r['python_error']],
 'rows':rows}
(dest/'comparison.json').write_text(json.dumps(summary,indent=2)+'\n')
np.savetxt(dest/'fortran_full_flow_table.csv',full,delimiter=',',
 header='pressure_pa,mean_gap_m,phi_p,phi_s,xmult,area_direct,area_boundary,area_u1min',comments='')
print('SUMMARY',json.dumps({k:v for k,v in summary.items() if k not in ('rows','source','archive','run')}),flush=True)

shutil.copyfile(f,dest/'fortran_flow_factors.txt')
