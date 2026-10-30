"""Module 3: constant-Newtonian fluid and independent Navier wall boundaries.

Local pressure is absolute [Pa], temperature [K], shear-rate magnitude [1/s].
The viscosity interface reserves these arguments; this release uses constant eta.
"""
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


def nonnegative(value, name):
    if isinstance(value, (bool, np.bool_)):
        raise ValueError(f'{name} must be nonnegative and finite')
    try:
        number=float(value)
    except (TypeError,ValueError) as exc:
        raise ValueError(f'{name} must be nonnegative and finite') from exc
    if not math.isfinite(number) or number<0:
        raise ValueError(f'{name} must be nonnegative and finite')
    return number


@dataclass(frozen=True)
class LocalFlowResponse:
    flux_per_width_m2_s: object
    lower_wall_shear_stress_pa: object
    upper_wall_shear_stress_pa: object
    lower_fluid_velocity_m_s: object
    upper_fluid_velocity_m_s: object


@dataclass(frozen=True)
class NewtonianNoSlip:
    """Incompressible constant-viscosity liquid. Walls at z=0 and z=h."""
    viscosity_pa_s: float
    gap_power = 3

    def __post_init__(self):
        object.__setattr__(self, 'viscosity_pa_s', positive(self.viscosity_pa_s, 'viscosity_pa_s'))

    @property
    def slip_lengths_m(self):
        return 0.,0.

    @property
    def has_slip(self):
        return any(b>0 for b in self.slip_lengths_m)

    @property
    def metadata(self):
        return {'rheology':'newtonian','viscosity_model':'constant',
                'wall_slip':'navier' if isinstance(self,NewtonianNavierSlip) else 'no_slip',
                'viscosity_pa_s':self.viscosity_pa_s,
                'lower_slip_length_m':self.slip_lengths_m[0],
                'upper_slip_length_m':self.slip_lengths_m[1]}

    @property
    def inverse_mobility_prefactor(self):
        return 12*self.viscosity_pa_s

    @staticmethod
    def _gaps(gap_m):
        gaps = np.asarray(gap_m, dtype=float)
        if np.any(~np.isfinite(gaps)) or np.any(gaps < 0):
            raise ValueError('gap_m must be finite and nonnegative')
        return gaps

    def viscosity(self, shear_rate_s_inv=0., fluid_pressure_pa=1e5, temperature_k=298.15):
        """Constitutive hook; p is fluid absolute pressure, NOT contact pressure.

        No pressure/temperature dependence or non-Newtonian response is enabled.
        Future models must also adapt homogenization and the spatial solver.
        """
        rate,p,t=np.broadcast_arrays(np.asarray(shear_rate_s_inv,dtype=float),
            np.asarray(fluid_pressure_pa,dtype=float),np.asarray(temperature_k,dtype=float))
        if (any(np.any(~np.isfinite(v)) for v in (rate,p,t)) or
                np.any(rate<0) or np.any(p<0) or np.any(t<=0)):
            raise ValueError('Require finite shear-rate magnitude >=0, fluid absolute pressure >=0 and temperature >0 K')
        return np.full(rate.shape,self.viscosity_pa_s)

    def geometric_mobility(self, gap_m):
        """12 eta times local pressure mobility [m^3], including wall slip."""
        h=self._gaps(gap_m);b0,b1=self.slip_lengths_m
        if b0==0 and b1==0:
            return h**3
        denominator=h+b0+b1
        return h**3*(1+3*(b0+b1)/denominator)+12*h**2*(b0/denominator)*b1

    def pressure_mobility(self, gap_m):
        """Local mobility [m^3/(Pa s)]; true contact always has zero flux."""
        return self.geometric_mobility(gap_m)/self.inverse_mobility_prefactor

    def local_response(self, gap_m, pressure_gradient_pa_m,
                       lower_wall_velocity_m_s=0., upper_wall_velocity_m_s=0., *,
                       fluid_pressure_pa=1e5, temperature_k=298.15):
        """Exact planar Poiseuille/Couette solution with independent Navier walls.

        u(0)-U0=b0*u'(0); u(h)-U1=-b1*u'(h).
        Tractions are exerted BY FLUID ON each wall in the positive x direction:
        lower=eta*u'(0), upper=-eta*u'(h). At h=0 q=0, velocities/tractions
        are NaN: there is no fluid layer or implied solid-contact friction law.
        This local response does not supply homogenized shear-flow factors.
        """
        h,g,u0,u1,eta=np.broadcast_arrays(self._gaps(gap_m),
            np.asarray(pressure_gradient_pa_m,dtype=float),
            np.asarray(lower_wall_velocity_m_s,dtype=float),
            np.asarray(upper_wall_velocity_m_s,dtype=float),
            self.viscosity(0.,fluid_pressure_pa,temperature_k))
        if any(np.any(~np.isfinite(v)) for v in (g,u0,u1)):
            raise ValueError('Pressure gradient and wall velocities must be finite')
        b0,b1=self.slip_lengths_m
        denominator=h+b0+b1
        safe=np.where(h>0,denominator,1.)
        entrainment=h*((h/2+b1)/safe*u0+(h/2+b0)/safe*u1)
        q=-self.geometric_mobility(h)*g/(12*eta)+entrainment
        # eta*A, where u(z)=U0+b0*A+A*z+g*z^2/(2 eta).
        shear0=(eta*(u1-u0)-g*h*(h/2+b1))/safe
        shear1=-shear0-g*h
        velocity0=u0+b0*shear0/eta
        velocity1=u1+b1*shear1/eta
        defined=lambda v:np.where(h>0,v,np.nan)
        return LocalFlowResponse(np.where(h>0,q,0.),defined(shear0),defined(shear1),
                                 defined(velocity0),defined(velocity1))

    def flux(self, gap_m, pressure_gradient_pa_m, lower_wall_velocity_m_s=0.,
             upper_wall_velocity_m_s=0., *, fluid_pressure_pa=1e5, temperature_k=298.15):
        return self.local_response(gap_m,pressure_gradient_pa_m,lower_wall_velocity_m_s,
            upper_wall_velocity_m_s,fluid_pressure_pa=fluid_pressure_pa,
            temperature_k=temperature_k).flux_per_width_m2_s


@dataclass(frozen=True)
class NewtonianNavierSlip(NewtonianNoSlip):
    lower_slip_length_m: float = 0.
    upper_slip_length_m: float = 0.

    def __post_init__(self):
        super().__post_init__()
        for name in ('lower_slip_length_m','upper_slip_length_m'):
            object.__setattr__(self,name,nonnegative(getattr(self,name),name))

    @property
    def slip_lengths_m(self):
        return self.lower_slip_length_m,self.upper_slip_length_m


def make_interface(parameters):
    """Select one rheology and one wall law; read only their active presets.

    Old flat viscosity_pa_s remains supported, but cannot coexist with an active
    rheology_parameters preset. This avoids ambiguous viscosity sources.
    """
    if not isinstance(parameters, dict):
        raise ValueError('interface must be a parameter dictionary')
    allowed={'rheology','wall_slip','viscosity_pa_s','rheology_parameters','slip_parameters'}
    if set(parameters)-allowed:
        raise ValueError(f'Unknown interface parameters: {sorted(set(parameters)-allowed)}')
    rheology=parameters.get('rheology','newtonian')
    if rheology!='newtonian':
        raise NotImplementedError('Only newtonian constant viscosity is implemented')
    viscosity=parameters.get('viscosity_pa_s')
    if 'rheology_parameters' in parameters:
        presets=parameters['rheology_parameters']
        if not isinstance(presets,dict) or not isinstance(presets.get(rheology),dict):
            raise ValueError('Selected rheology needs a matching rheology_parameters dictionary')
        active=presets[rheology]
        if set(active)!={'viscosity_pa_s'}:
            raise ValueError('newtonian requires only viscosity_pa_s')
        if 'viscosity_pa_s' in parameters:
            raise ValueError('Give viscosity once, in rheology_parameters or the legacy flat field')
        viscosity=active['viscosity_pa_s']
    if viscosity is None:
        raise ValueError('newtonian requires viscosity_pa_s')
    wall=parameters.get('wall_slip','no_slip')
    if wall not in ('no_slip','navier'):
        raise NotImplementedError('wall_slip must be no_slip or navier')
    active={}
    if 'slip_parameters' in parameters:
        presets=parameters['slip_parameters']
        if not isinstance(presets,dict) or not isinstance(presets.get(wall),dict):
            raise ValueError('Selected wall law needs a matching slip_parameters dictionary')
        active=presets[wall]
    if wall=='no_slip':
        if active:raise ValueError('no_slip takes no slip parameters')
        return NewtonianNoSlip(viscosity)
    if set(active)!={'lower_slip_length_m','upper_slip_length_m'}:
        raise ValueError('navier requires lower_slip_length_m and upper_slip_length_m')
    return NewtonianNavierSlip(viscosity,**active)
