"""Fixed-contact leakage demonstration; mixed-friction coupling remains reserved."""
from dataclasses import dataclass
import math

from .interface import positive


@dataclass(frozen=True)
class LeakageParameters:
    seal_length_m: float
    seal_width_m: float
    inlet_pressure_pa: float
    outlet_pressure_pa: float

    def __post_init__(self):
        for name in ('seal_length_m','seal_width_m'):
            object.__setattr__(self,name,positive(getattr(self,name),name))
        for name in ('inlet_pressure_pa','outlet_pressure_pa'):
            value=getattr(self,name)
            if value is None or isinstance(value,bool) or not math.isfinite(float(value)):
                raise ValueError(f'{name} must be finite')
            object.__setattr__(self,name,float(value))


def run_leakage(contact_flow, parameters, interface):
    """Uniform, fixed contact; incompressible Q=(width/length)*sigma*Delta p.

    Fluid pressures do not alter the supplied dry contact pressure or PDF.
    Gauge pressures may be negative, and reversing Delta p reverses Q.
    """
    if isinstance(parameters,dict):
        parameters=LeakageParameters(**parameters)
    delta=parameters.inlet_pressure_pa-parameters.outlet_pressure_pa
    factor=parameters.seal_width_m/parameters.seal_length_m*delta
    rate=factor*contact_flow.hydraulic_conductivity_m3_pa_s
    upper=(0. if factor==0 else None if contact_flow.upper_bound_unbounded else
           factor*contact_flow.upper_geometric_conductivity_m3/(12*interface.viscosity_pa_s))
    return {'pressure_pa':contact_flow.pressure_pa,'pressure_drop_pa':delta,
            'flow_rate_m3_s':rate,
            'min_flow_rate_m3_s':None if upper is None else min(rate,upper),
            'max_flow_rate_m3_s':None if upper is None else max(rate,upper),
            'bound_unbounded':contact_flow.upper_bound_unbounded and factor!=0,
            'status':contact_flow.status}


def run_mixed_friction(contact_response, parameters):
    raise NotImplementedError('混合摩擦求解器尚未实现，当前只预留接口。')
