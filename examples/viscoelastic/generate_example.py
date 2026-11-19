"""Generate an explicitly SYNTHETIC standard-linear-solid, in Persson Hz/MPa format."""
from pathlib import Path
import numpy as np

root=Path(__file__).resolve().parent
frequency=np.geomspace(1e-10,1e12,2201)
x=2*np.pi*frequency*.01
storage=1e6+9e6*x*x/(1+x*x);loss=9e6*x/(1+x*x)
with (root/'master_curve.dat').open('w',encoding='utf-8') as f:
    f.write('SYNTHETIC_SLS_E0_1MPa_Einf_10MPa_tau_0.01s_Tref_20C\n')
    np.savetxt(f,np.c_[frequency,storage/1e6,loss/1e6],fmt='%.15e')
(root/'shift_factors.dat').write_text('SYNTHETIC_SLS_Tref_20C\n10 10 1\n20 1 1\n30 0.1 1\n',encoding='utf-8')
(root/'roughness.dat').write_text('SYNTHETIC_power_law\n1e3\n1e3\n1e6\n0\n1\n3 -18\n6 -30\n',encoding='utf-8')
