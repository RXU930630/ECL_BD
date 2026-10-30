"""Module 6: nonlocal plane-strain half-space response and dry line contact.

Piecewise constant tractions on uniform cells; the logarithmic self term is
integrated analytically. Positive displacement opens the relative gap. For the
upper-elastic/lower-rigid representation the upper surface moves in +z.
No fluid coupling, finite thickness, axisymmetry or whole-ring deformation.
"""
from dataclasses import dataclass
import math

import numpy as np
from scipy.linalg import solve,toeplitz


def _positive(value,name):
    if isinstance(value,(bool,np.bool_)):
        raise ValueError(f'{name} must be finite and positive')
    try:result=float(value)
    except (TypeError,ValueError) as exc:raise ValueError(f'{name} must be finite and positive') from exc
    if not np.isfinite(result) or result<=0:raise ValueError(f'{name} must be finite and positive')
    return result


@dataclass(frozen=True)
class HertzLine:
    radius_m: float
    line_load_n_m: float
    effective_modulus_pa: float

    def __post_init__(self):
        for name in ('radius_m','line_load_n_m','effective_modulus_pa'):
            object.__setattr__(self,name,_positive(getattr(self,name),name))

    @property
    def half_width_m(self):
        return math.sqrt(4*self.line_load_n_m*self.radius_m/(math.pi*self.effective_modulus_pa))

    @property
    def peak_pressure_pa(self):
        return 2*self.line_load_n_m/(math.pi*self.half_width_m)

    def profile(self,x_m):return np.asarray(x_m,dtype=float)**2/(2*self.radius_m)

    def cell_pressure(self,edges_m):
        a=self.half_width_m
        z=np.clip(np.asarray(edges_m)/a,-1,1)
        primitive=.5*(z*np.sqrt(np.maximum(0,1-z*z))+np.arcsin(z))
        return self.peak_pressure_pa*a*np.diff(primitive)/np.diff(edges_m)

    def gap(self,x_m):
        y=abs(np.asarray(x_m,dtype=float))/self.half_width_m
        gap=np.zeros_like(y);outside=y>1
        t=y[outside]
        gap[outside]=self.half_width_m**2/(2*self.radius_m)*(t*np.sqrt(t*t-1)-np.arccosh(t))
        return gap


class HalfSpaceLine:
    """Relative normal compliance, u=2/(pi E*) int p(s) log(Lref/|x-s|) ds.

    Lref fixes the arbitrary displacement datum of 2D line loading. Requiring
    Lref >= domain length keeps a convenient positive compliance matrix; a
    load-controlled offset absorbs a different reference choice.
    """
    def __init__(self,half_width_m,cells,effective_modulus_pa,reference_length_m=None):
        half=_positive(half_width_m,'domain half-width')
        if type(cells) is not int or not 17<=cells<=2049:
            raise ValueError('cells must be an integer in [17,2049]')
        self.effective_modulus_pa=_positive(effective_modulus_pa,'effective modulus')
        self.edges_m=np.linspace(-half,half,cells+1)
        self.x_m=(self.edges_m[:-1]+self.edges_m[1:])/2
        self.dx_m=2*half/cells
        self.reference_length_m=2*half if reference_length_m is None else _positive(reference_length_m,'reference length')
        if self.reference_length_m<2*half:
            raise ValueError('reference_length_m must be at least the full domain length')
        r=np.arange(cells,dtype=float)
        def f(t):
            return t*np.log(abs(t))
        integrated=self.dx_m*(1+math.log(self.reference_length_m/self.dx_m)
                              -f(r+.5)+f(r-.5))
        self.compliance_m_pa=2/(math.pi*self.effective_modulus_pa)*toeplitz(integrated)

    def displacement(self,pressure_pa):
        p=np.asarray(pressure_pa,dtype=float)
        if p.shape!=self.x_m.shape or np.any(~np.isfinite(p)):
            raise ValueError('Pressure must be a finite vector on the elastic cells')
        return self.compliance_m_pa@p

    def line_load(self,pressure_pa):return float(self.dx_m*np.sum(pressure_pa))


@dataclass(frozen=True)
class MacroContactSolution:
    x_m: np.ndarray
    cell_edges_m: np.ndarray
    pressure_pa: np.ndarray
    relative_contact_area: np.ndarray
    geometric_profile_m: np.ndarray
    displacement_m: np.ndarray
    mean_gap_m: np.ndarray
    rigid_offset_m: float
    line_load_n_m: float
    iterations: int
    equilibrium_relative_residual: float
    load_relative_error: float
    status: str
    minimum_gap_m: float
    edge_pressure_ratio: float


def solve_smooth_contact(elastic,hertz,*,tolerance=1e-9,max_iterations=2000):
    """Load-constrained convex contact problem; active-set complementarity solve."""
    if tolerance<=0 or not np.isfinite(tolerance) or type(max_iterations) is not int or max_iterations<1:
        raise ValueError('Invalid contact tolerance/iteration limit')
    if not math.isclose(elastic.effective_modulus_pa,hertz.effective_modulus_pa,rel_tol=1e-13):
        raise ValueError('Elastic operator and geometry reference must use the same E*')
    if elastic.edges_m[-1]<=hertz.half_width_m:
        raise ValueError('Domain must extend beyond the Hertz contact half-width')
    x=elastic.x_m;g=hertz.profile(x);a=hertz.half_width_m;p0=hertz.peak_pressure_pa
    scale=a*a/(2*hertz.radius_m)
    c=elastic.compliance_m_pa*p0/scale;shape=g/scale
    target=hertz.line_load_n_m/(elastic.dx_m*p0)
    active=abs(x)<a
    p=np.zeros_like(x);offset=0.
    for iteration in range(1,max_iterations+1):
        indices=np.flatnonzero(active);n=len(indices)
        if n==0:raise ValueError('Smooth contact active set became empty')
        matrix=np.empty((n+1,n+1));matrix[:n,:n]=c[np.ix_(indices,indices)]
        matrix[:n,n]=1.;matrix[n,:n]=1.;matrix[n,n]=0.
        values=solve(matrix,np.r_[-shape[indices],target],assume_a='sym')
        p[:]=0;p[indices]=values[:-1];offset=float(values[-1])
        if np.min(p)<-tolerance:
            active[int(np.argmin(p))]=False;continue
        gap=shape+offset+c@p
        inactive=np.flatnonzero(~active)
        if len(inactive) and np.min(gap[inactive]) < -tolerance:
            active[inactive[int(np.argmin(gap[inactive]))]]=True;continue
        if np.any(p<0):
            active[np.flatnonzero(p<0)]=False;continue
        break
    else:raise ValueError('Smooth contact active set did not converge')
    pressure=p*p0;u=elastic.displacement(pressure);h=g+offset*scale+u
    residual=max(float(np.max(abs(h[active])))/scale,
                 max(0.,-float(np.min(h)))/scale,max(0.,-float(np.min(p))))
    load_error=abs(elastic.line_load(pressure)/hertz.line_load_n_m-1)
    if max(residual,load_error)>10*tolerance:raise ValueError('Smooth contact equilibrium check failed')
    return MacroContactSolution(x.copy(),elastic.edges_m.copy(),pressure,(pressure>0).astype(float),g,u,h,
        offset*scale,elastic.line_load(pressure),iteration,residual,load_error,'smooth_converged',
        float(np.min(h)),float(max(pressure[0],pressure[-1])/np.max(pressure)))
