#!/usr/bin/env python3
"""Local dependencies, complete runtime imports and one exact smooth-film check."""
from importlib import import_module
from pathlib import Path
import sys

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'src'))


def main():
    if sys.version_info<(3,11):
        print('需要Python 3.11或更新版本，推荐Python 3.12。');return 2
    try:
        import numpy as np
        import scipy
        import matplotlib
        import persson_contact
        from persson_contact.interface import NewtonianNoSlip
        from persson_contact.shear_fit import Experiment,fit_constant_strength
        for path in sorted((ROOT/'src/persson_contact').glob('*.py')):
            if path.stem!='__main__':import_module('persson_contact.'+path.stem)
        response=NewtonianNoSlip(.1).local_response(1e-6,0.,0.,.01)
        np.testing.assert_allclose(response.flux_per_width_m2_s,5e-9,rtol=1e-14)
        result,_=fit_constant_strength([Experiment(2,.01,100.,.25)],np.array([1e-6]),np.array([.05]))
        np.testing.assert_allclose(result['fitted_shear_strength_pa'],2e5,rtol=1e-14)
    except ModuleNotFoundError as exc:
        print(f'缺少依赖：{exc.name}。请执行 "{sys.executable}" -m pip install -r "{ROOT / "requirements.txt"}"');return 2
    print(f'ECL_BD {persson_contact.__version__}; Python {sys.version.split()[0]}')
    print(f'numpy {np.__version__}; scipy {scipy.__version__}; matplotlib {matplotlib.__version__}')
    print('运行模块、光滑流动及剪切强度拟合自检通过。此检查不替代实验验证。')
    return 0


if __name__=='__main__':raise SystemExit(main())
