"""Module 7: prescribed liquid-pressure boundary cases (common gauge, Pa)."""
from dataclasses import dataclass
import numpy as np


@dataclass(frozen=True)
class PressureBoundaryCases:
    inlet_pressures_pa: np.ndarray
    outlet_pressure_pa: float


def read_pressure_boundaries(parameters):
    if ('inlet_pressure_pa' in parameters)==('inlet_pressures_pa' in parameters):
        raise ValueError('Provide inlet_pressure_pa OR inlet_pressures_pa')
    values=parameters.get('inlet_pressures_pa',parameters.get('inlet_pressure_pa'))
    if any(isinstance(v,(bool,np.bool_)) for v in np.atleast_1d(np.asarray(values,dtype=object)).ravel()):
        raise ValueError('Pressure boundaries must be numbers, not booleans')
    inlet=np.atleast_1d(np.asarray(values,dtype=float))
    if inlet.ndim!=1 or len(inlet)==0 or np.any(~np.isfinite(inlet)) or len(np.unique(inlet))!=len(inlet):
        raise ValueError('inlet_pressures_pa must be a nonempty finite list without duplicates')
    value=parameters.get('outlet_pressure_pa',0.)
    if value is None or isinstance(value,(bool,np.bool_)):raise ValueError('outlet_pressure_pa must be a number')
    out=float(value)
    if not np.isfinite(out):raise ValueError('outlet_pressure_pa must be finite')
    inlet.flags.writeable=False
    return PressureBoundaryCases(inlet,out)
