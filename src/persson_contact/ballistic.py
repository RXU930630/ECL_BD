"""Independent quadratic effective-medium transport for Fortran mode 5000."""
from dataclasses import asdict,dataclass
import math

from .flow import FlowOptions,integrate_mobility
from .gas import BallisticGas


@dataclass(frozen=True)
class BallisticFlowResult:
    pressure_pa: float
    mean_gap_m: float
    relative_contact_area: float
    effective_dimension: float
    geometric_conductivity_m2: float
    molecular_mobility_per_pa_s: float
    pressure_flow_factor: float
    upper_geometric_conductivity_m2: float | None
    upper_pressure_flow_factor: float | None
    upper_bound_unbounded: bool
    retained_pdf_probability: float
    unrepresented_noncontact_probability: float
    zero_conductance_probability: float
    closure_residual: float
    root_iterations: int
    gap_points: int
    quadrature_relative_difference: float
    status: str
    percolation_relative_contact_area: float
    unresolved_flow_policy: str

    def as_dict(self):
        return asdict(self)


def ballistic_pressure_flow(model,pressure_pa,interface,options=None,*,percolation_relative_contact_area=.42):
    if not isinstance(interface,BallisticGas):
        raise ValueError('ballistic_pressure_flow requires BallisticGas')
    options=options or FlowOptions()
    p,a,h,fine,count,difference=integrate_mobility(model,pressure_pa,interface,options,
        percolation_relative_contact_area=percolation_relative_contact_area)
    upper=fine['upper'];unbounded=not math.isfinite(upper)
    return BallisticFlowResult(p,h,a,fine['dimension'],fine['geometric'],
        fine['geometric']/interface.inverse_mobility_prefactor,fine['factor'],
        None if unbounded else upper,None if unbounded else (upper/h**2 if h>0 else 0.),unbounded,
        fine['mass'],fine['missing'],fine['zero'],fine['residual'],fine['iterations'],
        count,difference,('blocked' if upper==0 else 'unresolved_connectivity')
        if fine['geometric']==0 else 'conducting',percolation_relative_contact_area,options.unresolved_flow_policy)
