"""Isothermal gas properties, local ballistic law and Fortran-500 estimates.

Gas volumes are always labelled with their reference pressure and temperature.
The estimates reuse a liquid solution; they do not solve compressible flow.
"""
from dataclasses import dataclass
import math
import numpy as np

from .interface import NewtonianNoSlip, positive


BOLTZMANN_PA_M3_K = 1.380649e-23


@dataclass(frozen=True)
class GasProperties:
    temperature_c: float = 25.0
    mean_molecular_speed_m_s: float = 471.0
    viscosity_pa_s: float = 1.82e-5
    ambient_pressure_pa: float = 1e5
    reference_pressure_pa: float = 1e5

    def __post_init__(self):
        if isinstance(self.temperature_c, bool):
            raise ValueError('temperature_c must be a finite number above absolute zero')
        try:
            temperature=float(self.temperature_c)
        except (TypeError,ValueError) as exc:
            raise ValueError('temperature_c must be a finite number above absolute zero') from exc
        if not math.isfinite(temperature) or temperature<=-273.15:
            raise ValueError('temperature_c must be above absolute zero')
        object.__setattr__(self,'temperature_c',temperature)
        for name in ('mean_molecular_speed_m_s','viscosity_pa_s','ambient_pressure_pa','reference_pressure_pa'):
            object.__setattr__(self,name,positive(getattr(self,name),name))

    @property
    def temperature_k(self):
        return self.temperature_c+273.15

    @property
    def thermal_energy_j(self):
        return BOLTZMANN_PA_M3_K*self.temperature_k

    def absolute_pressure(self,gauge_pressure_pa):
        value=np.asarray(gauge_pressure_pa,dtype=float)+self.ambient_pressure_pa
        if np.any(~np.isfinite(value)) or np.any(value<0):
            raise ValueError('Gas boundary absolute pressures must be nonnegative')
        return value

    def reference_volume(self,number_rate):
        return number_rate*self.thermal_energy_j/self.reference_pressure_pa

    def mean_free_path(self,absolute_pressure_pa):
        pressure=np.asarray(absolute_pressure_pa,dtype=float)
        if np.any(~np.isfinite(pressure)) or np.any(pressure<0):
            raise ValueError('Absolute pressure must be finite and nonnegative')
        numerator=3*np.pi*self.viscosity_pa_s*self.mean_molecular_speed_m_s/8
        with np.errstate(divide='ignore'):
            return numerator/pressure


def read_gas_properties(parameters):
    if not isinstance(parameters,dict):
        raise ValueError('leakage.gas must be a dictionary')
    unknown=set(parameters)-set(GasProperties.__dataclass_fields__)
    if unknown:
        raise ValueError(f'Unknown gas parameters: {sorted(unknown)}')
    return GasProperties(**parameters)


@dataclass(frozen=True)
class BallisticGas:
    """Original 5000 wall-collision law j_N=-vbar*u^2/(2*kBT)*grad(p)."""
    gas: GasProperties
    gap_power = 2

    @property
    def inverse_mobility_prefactor(self):
        return 2*self.gas.thermal_energy_j/self.gas.mean_molecular_speed_m_s

    def pressure_mobility(self,gap_m):
        return NewtonianNoSlip._gaps(gap_m)**2/self.inverse_mobility_prefactor

    def geometric_mobility(self,gap_m):
        return NewtonianNoSlip._gaps(gap_m)**2


def equivalent_gas_estimates(liquid_flow_m3_s,pressure_drop_pa,length_m,width_m,
                             liquid_viscosity_pa_s,inlet_gauge_pa,outlet_gauge_pa,gas):
    """500 auxiliary estimates at the SAME liquid-derived contact state.

    K3_seal=12*eta*L*Q/(W*DeltaP), u_eff=K3_seal^(1/3).
    Diffusive: Qref=Qliq*(eta_liq/eta_gas)*pmean_abs/pref.
    Ballistic estimate: Qref=Qliq*12*eta_liq*vbar/(2*pref*u_eff).
    """
    q,dp=float(liquid_flow_m3_s),float(pressure_drop_pa)
    if not np.all(np.isfinite([q,dp,inlet_gauge_pa,outlet_gauge_pa])):
        raise ValueError('Equivalent gas estimates require finite liquid rates and pressures')
    length=positive(length_m,'seal_length_m');width=positive(width_m,'seal_width_m')
    eta=positive(liquid_viscosity_pa_s,'liquid_viscosity_pa_s')
    boundary=gas.absolute_pressure([inlet_gauge_pa,outlet_gauge_pa])
    if q!=0 and (dp==0 or q*dp<0):
        raise ValueError('Liquid flow must have the sign of the pressure drop')
    gap=None if dp==0 else (12*eta*length*abs(q)/(width*abs(dp)))**(1/3)
    diffusive=q*(eta/gas.viscosity_pa_s)*float(np.mean(boundary))/gas.reference_pressure_pa
    ballistic=0. if q==0 else q*12*eta*gas.mean_molecular_speed_m_s/(2*gas.reference_pressure_pa*gap)
    return {'equivalent_hydraulic_gap_m':gap,
        'estimated_diffusive_gas_reference_flow_m3_s':diffusive,
        'estimated_ballistic_gas_reference_flow_m3_s':ballistic,
        'estimated_diffusive_gas_reference_flow_mm3_h':diffusive*3.6e12,
        'estimated_ballistic_gas_reference_flow_mm3_h':ballistic*3.6e12,
        'estimated_diffusive_molecular_flow_per_s':diffusive*gas.reference_pressure_pa/gas.thermal_energy_j,
        'estimated_ballistic_molecular_flow_per_s':ballistic*gas.reference_pressure_pa/gas.thermal_energy_j,
        'gas_reference_pressure_pa':gas.reference_pressure_pa,'gas_temperature_k':gas.temperature_k}
