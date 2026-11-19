"""One-command synthetic dry viscoelastic benchmark. No real rubber is calibrated."""
from pathlib import Path
import argparse,sys
root=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(root/'src'))
from persson_contact.workflow import run_project

parser=argparse.ArgumentParser(description=__doc__)
parser.add_argument('--output',type=Path,default=root/'results/run')
parser.add_argument('--temperature',type=float,default=20.)
args=parser.parse_args()
p={
 'material':{'model':'viscoelastic','master_curve_file':'master_curve.dat','poisson_ratio':.5},
 'psd':{'path':'roughness.dat','scale_c':1.,'tail_is_extrapolated':False},
 'macro_contact':{'deformable_body':'lower'},
 'mixed_friction':{'rough_surface':'upper'},
 'contact_friction':{'model':'constant_shear_strength','model_parameters':{'constant_shear_strength':{'shear_strength_pa':1e5}}},
 'viscoelastic':{'reference_temperature_c':20.,'temperature_c':args.temperature,
     'shift_factors_file':'shift_factors.dat','modulus_type':'young','frequency_policy':'error',
     'pressures_pa':[1e4,1e5,1e6], 'speeds_m_s':[1e-6,1e-5,1e-4,1e-3,1e-2,.1],
     'q_intervals':512,'angular_points':64,'kappa':.5,'relative_tolerance':.002},
 'output_directory':str(args.output.resolve())}
result=run_project(p,base_directory=Path(__file__).resolve().parent,run_contact=False,run_viscoelastic_friction=True)
print('Synthetic dry benchmark:',result['numerical_status'])
