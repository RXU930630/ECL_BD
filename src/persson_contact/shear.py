"""Module 4: isotropic, single-rough-surface Newtonian shear transport.

Uses the original Sealing3 Persson interpolation, not an independently derived
PDF shear EMT. Positive phi_s is an amplitude; rough-surface ownership fixes sign.
"""
from dataclasses import asdict,dataclass
import math

import numpy as np


def _finite(value,name):
    if isinstance(value,(bool,np.bool_)):
        raise ValueError(f'{name} must be finite')
    try:value=float(value)
    except (TypeError,ValueError) as exc:raise ValueError(f'{name} must be finite') from exc
    if not math.isfinite(value):raise ValueError(f'{name} must be finite')
    return value


@dataclass(frozen=True)
class ShearOptions:
    model: str = 'persson_interpolation'
    rough_surface: str = 'lower'
    lower_wall_velocity_m_s: float = 0.
    upper_wall_velocity_m_s: float = .01
    pressure_gradient_pa_m: float = 0.
    pressure_floor_pa: float = 0.

    def __post_init__(self):
        if self.model!='persson_interpolation':
            raise NotImplementedError('Only isotropic persson_interpolation shear flow is implemented')
        if self.rough_surface not in ('lower','upper'):
            raise ValueError('rough_surface must be lower or upper (one rough wall, one smooth wall)')
        for name in ('lower_wall_velocity_m_s','upper_wall_velocity_m_s',
                     'pressure_gradient_pa_m','pressure_floor_pa'):
            object.__setattr__(self,name,_finite(getattr(self,name),name))

    @property
    def moving(self):
        return self.lower_wall_velocity_m_s!=0 or self.upper_wall_velocity_m_s!=0


def read_shear_options(raw):
    if not isinstance(raw,dict):raise ValueError('shear_flow must be a dictionary')
    unknown=set(raw)-set(ShearOptions.__dataclass_fields__)
    if unknown:raise ValueError(f'Unknown shear_flow parameters: {sorted(unknown)}')
    return ShearOptions(**raw)


def shear_length(mean_gap_m,rms_height_m,percolation_gap_m):
    """L_s=hrms*phi_s, with gamma=1 (isotropy); bounded in [0,h].

    phi_s=[hrms/h+2 max(h-hc,0)/(3 hrms)]^-1.
    The exact below-hc branch replaces the original fixed 1e-10 m cutoff.
    Smooth-surface limit hrms=0 has L_s=0.
    """
    h=np.asarray(mean_gap_m,dtype=float)
    rms=_finite(rms_height_m,'rms_height_m');hc=_finite(percolation_gap_m,'percolation_gap_m')
    if np.any(~np.isfinite(h)) or np.any(h<0) or rms<0 or hc<0:
        raise ValueError('Mean gap, RMS height and percolation gap must be finite and nonnegative')
    if rms==0:return np.zeros_like(h)
    return h/(1+2*(h/rms)*(np.maximum(h-hc,0)/rms)/3)


@dataclass(frozen=True)
class ShearFlowLaw:
    options: ShearOptions
    rms_height_m: float
    percolation_gap_m: float

    @classmethod
    def from_model(cls,options,model,critical_relative_contact_area):
        pc=model.pressure_for_relative_contact_area(critical_relative_contact_area)
        return cls(options,float(model.psd.hrms),float(model.mean_gap(pc)))

    def coefficients(self,mean_gap_m):
        h=np.asarray(mean_gap_m,dtype=float)
        length=shear_length(h,self.rms_height_m,self.percolation_gap_m)
        signed=length if self.options.rough_surface=='upper' else -length
        return (h-signed)/2,(h+signed)/2

    def source(self,mean_gap_m):
        lower,upper=self.coefficients(mean_gap_m)
        return lower*self.options.lower_wall_velocity_m_s+upper*self.options.upper_wall_velocity_m_s

    @property
    def metadata(self):
        return {**asdict(self.options),'rms_height_m':self.rms_height_m,
                'percolation_gap_m':self.percolation_gap_m,
                'scope':'Single rough wall, isotropic constant-Newtonian no-slip interpolation; no friction or load/film-thickness coupling.',
                'flux_convention':'q=-K/(12 eta)*dp/dx+(U0+U1)*h/2+sign*(U1-U0)*hrms*phi_s/2; sign=-1 for rough lower, +1 for rough upper.',
                'coordinate_convention':'Both wall velocities along +x; lower/upper refer to film-normal z, not inlet/outlet.'}
