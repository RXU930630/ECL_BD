"""Module 3: local fluid/interface laws, independent of contact homogenization."""
from dataclasses import dataclass
import math

import numpy as np


def positive(value, name):
    try:
        number=float(value)
    except (TypeError, ValueError) as exc:
        raise ValueError(f'{name} must be positive and finite') from exc
    if isinstance(value, (bool, np.bool_)) or not math.isfinite(number) or number <= 0:
        raise ValueError(f'{name} must be positive and finite')
    return number


@dataclass(frozen=True)
class NewtonianNoSlip:
    """Incompressible isoviscous fluid between two no-slip walls, SI units."""
    viscosity_pa_s: float

    def __post_init__(self):
        object.__setattr__(self, 'viscosity_pa_s', positive(self.viscosity_pa_s, 'viscosity_pa_s'))

    @staticmethod
    def _gaps(gap_m):
        gaps = np.asarray(gap_m, dtype=float)
        if np.any(~np.isfinite(gaps)) or np.any(gaps < 0):
            raise ValueError('gap_m must be finite and nonnegative')
        return gaps

    def pressure_mobility(self, gap_m):
        """Local mobility [m^3/(Pa s)]; zero in true contact."""
        return self._gaps(gap_m)**3 / (12*self.viscosity_pa_s)

    def flux(self, gap_m, pressure_gradient_pa_m, lower_wall_velocity_m_s=0.,
             upper_wall_velocity_m_s=0.):
        """Volume flux per width [m^2/s], Poiseuille + local Couette flow.

        This is a LOCAL wall law; no homogenized shear-flow factor is implied.
        """
        gaps = self._gaps(gap_m)
        gradient, lower, upper = (np.asarray(value, dtype=float) for value in
            (pressure_gradient_pa_m, lower_wall_velocity_m_s, upper_wall_velocity_m_s))
        if any(np.any(~np.isfinite(value)) for value in (gradient, lower, upper)):
            raise ValueError('Pressure gradient and wall velocities must be finite')
        return -self.pressure_mobility(gaps)*gradient + gaps*(lower+upper)/2


def make_interface(parameters):
    if not isinstance(parameters, dict):
        raise ValueError('interface must be a parameter dictionary')
    unknown = set(parameters)-{'rheology', 'wall_slip', 'viscosity_pa_s'}
    if unknown:
        raise ValueError(f'Unknown interface parameters: {sorted(unknown)}')
    if parameters.get('rheology', 'newtonian') != 'newtonian':
        raise NotImplementedError('Only newtonian rheology is implemented; no silent fallback.')
    if parameters.get('wall_slip', 'no_slip') != 'no_slip':
        raise NotImplementedError('Only no_slip walls are implemented; no silent fallback.')
    if 'viscosity_pa_s' not in parameters:
        raise ValueError('interface requires viscosity_pa_s')
    return NewtonianNoSlip(parameters['viscosity_pa_s'])
