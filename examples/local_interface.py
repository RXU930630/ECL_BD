"""Run with python examples/local_interface.py; no PSD or contact table needed."""
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
from persson_contact.interface import NewtonianNoSlip,NewtonianNavierSlip


def main():
    laws={
        'no_slip':NewtonianNoSlip(.1),
        'equal_slip_10nm':NewtonianNavierSlip(.1,1e-8,1e-8),
        'upper_slip_1um':NewtonianNavierSlip(.1,0,1e-6),
    }
    print('case,q_per_width_m2_s,fluid_on_lower_wall_Pa,fluid_on_upper_wall_Pa')
    for name,law in laws.items():
        v=law.local_response(1e-6,-1e8,0.,.1,fluid_pressure_pa=1e5,temperature_k=298.15)
        print(f'{name},{float(v.flux_per_width_m2_s):.12g},'
              f'{float(v.lower_wall_shear_stress_pa):.12g},'
              f'{float(v.upper_wall_shear_stress_pa):.12g}')


if __name__=='__main__':main()
