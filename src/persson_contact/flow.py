"""Module 4: isotropic pressure-flow homogenization, no macroscopic geometry.

Public liquid pressure_flow returns K3 in m^3 (the Fortran flowsigmaeff1).
Its conservative PDF integrator is shared with quadratic ballistic transport.
Hydraulic conductivity is K/(12 eta); phi_p=K/h^3. Shear flow is reserved.
The missing-PDF policy matches Fortran Sealing3 but is explicitly labelled:
zero hydraulic conductivity is not an additional assertion of real contact.
"""
from dataclasses import asdict, dataclass
import math

import numpy as np
from scipy.optimize import brentq
from scipy.special import expit

from .interface import positive


@dataclass(frozen=True)
class FlowOptions:
    closure: str = 'fortran_isotropic_bruggeman'
    gap_points: int = 2049
    quadrature_relative_tolerance: float = .002
    max_gap_refinements: int = 2
    unresolved_flow_policy: str = 'zero_conductance_with_upper_bound'

    def __post_init__(self):
        if self.closure != 'fortran_isotropic_bruggeman':
            raise NotImplementedError('Only fortran_isotropic_bruggeman pressure flow is implemented')
        if self.unresolved_flow_policy != 'zero_conductance_with_upper_bound':
            raise ValueError('Use explicit unresolved_flow_policy=zero_conductance_with_upper_bound')
        for name, minimum in (('gap_points',129), ('max_gap_refinements',0)):
            value = getattr(self,name)
            if type(value) is not int or value < minimum:
                raise ValueError(f'{name} must be an integer >= {minimum}')
        if self.gap_points > 65537 or self.max_gap_refinements > 4:
            raise ValueError('Flow resolution exceeds supported limits')
        if not 0 < positive(self.quadrature_relative_tolerance,'quadrature_relative_tolerance') < 1:
            raise ValueError('quadrature_relative_tolerance must lie in (0,1)')


def effective_dimension(relative_contact_area, percolation_relative_contact_area):
    """Sealing3 xnDrun: interpolate from 2 to 1/(1-percolation) with area."""
    a, pc = float(relative_contact_area), float(percolation_relative_contact_area)
    if not math.isfinite(a) or not 0 <= a <= 1 or not math.isfinite(pc) or not 0 < pc < 1:
        raise ValueError('Invalid relative contact area or percolation threshold')
    return 2 + (1/(1-pc)-2)*min(a/pc,1.)


def bruggeman_conductivity(conductivities, probabilities, zero_probability, dimension,
                          infinite_probability=0.):
    """Solve scalar effective-medium closure with explicit zero/infinite phases.

    E[(k-K)/(k+(dimension-1)K)]=0. The infinite phase is used ONLY to
    bound unknown PDF probability, never as the default physical fluid law.
    Return (K, dimensionless closure residual, root iterations). No floors.
    """
    k, w = np.asarray(conductivities,dtype=float), np.asarray(probabilities,dtype=float)
    if (k.ndim != 1 or k.shape != w.shape or np.any(~np.isfinite(k))
            or np.any(k <= 0) or np.any(~np.isfinite(w)) or np.any(w < 0)):
        raise ValueError('Require positive finite conductivities and nonnegative probabilities')
    z, inf, d = float(zero_probability), float(infinite_probability), float(dimension)
    if not all(math.isfinite(v) for v in (z,inf,d)) or min(z,inf) < 0 or d <= 1:
        raise ValueError('Invalid zero/infinite probability or effective dimension')
    if not math.isclose(float(sum(w))+z+inf,1.,rel_tol=0,abs_tol=2e-12):
        raise ValueError('Effective-medium phase probabilities must sum to one; no renormalization')
    lam = d-1
    if z >= lam/d or math.isclose(z,lam/d,rel_tol=0,abs_tol=8*np.finfo(float).eps):
        return 0., 0., 0
    if inf >= 1/d or math.isclose(inf,1/d,rel_tol=0,abs_tol=8*np.finfo(float).eps):
        return math.inf, 0., 0
    if not np.any(w > 0):
        raise ValueError('No finite conducting phase resolves the effective medium')
    active=w>0; logk=np.log(k[active]); w=w[active]
    loglam=math.log(lam)

    # Avoid subtracting 1 from a saturated sigmoid sum near percolation.
    # The small phase margin, not the O(1) terms, sets a tiny finite root.
    low_margin=math.fsum([*w,inf,-z/lam])
    high_margin=math.fsum([z/lam,*(w/lam),-inf])

    def residual(log_value):
        if low_margin <= high_margin:
            return low_margin-d/lam*float(w @ expit(loglam+log_value-logk))
        return -high_margin+d/lam*float(w @ expit(logk-loglam-log_value))

    # A log-conductivity root remains well conditioned over many gap decades.
    left, right = float(min(logk)-80), float(max(logk)+80)
    root, result = brentq(residual,left,right,xtol=1e-12,full_output=True)
    value = math.exp(root)
    if not result.converged:
        raise ValueError('Bruggeman pressure-flow solve did not converge')
    return value,float(residual(root)),result.iterations


@dataclass(frozen=True)
class PressureFlowResult:
    pressure_pa: float
    mean_gap_m: float
    relative_contact_area: float
    effective_dimension: float
    geometric_conductivity_m3: float
    hydraulic_conductivity_m3_pa_s: float
    pressure_flow_factor: float
    upper_geometric_conductivity_m3: float | None
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


def homogenize_bins(gap_edges_m, cdf, relative_contact_area, mean_gap_m, interface,
                    percolation_relative_contact_area=.42):
    """Conservative CDF bins, log-midpoint local mobility, explicit missing mass.

    Lower/default: assign unrepresented noncontact mass to zero conductance.
    Upper: same mass assigned infinite conductance, true contact stays zero.
    Neither case changes the contact model or renormalizes the resolved PDF.
    """
    edges, cdf = np.asarray(gap_edges_m,dtype=float),np.asarray(cdf,dtype=float)
    if (edges.ndim!=1 or edges.size<3 or edges.shape!=cdf.shape or edges[0]!=0
            or np.any(~np.isfinite(edges)) or np.any(np.diff(edges)<=0)
            or np.any(~np.isfinite(cdf)) or cdf[0]!=0 or np.any(np.diff(cdf)<0)):
        raise ValueError('Require increasing gap edges starting at zero and monotone CDF starting at zero')
    a=float(relative_contact_area);d=effective_dimension(a,percolation_relative_contact_area)
    h=positive(mean_gap_m,'mean_gap_m')
    probabilities=np.diff(cdf);mass=float(sum(probabilities))
    missing=1-a-mass
    if missing < -2e-12:
        raise ValueError('Resolved PDF exceeds noncontact probability; refine contact quadrature')
    if missing < 0:
        raise ValueError('Resolved PDF mass exceeds standard noncontact area; decrease pdf_log_step')
    centers=np.empty(len(edges)-1);centers[0]=edges[1]/2
    centers[1:]=np.sqrt(edges[1:-1]*edges[2:])
    # Normalize the local law; K3 scales with h^3, ballistic K2 with h^2.
    power=interface.gap_power
    reference=interface.geometric_mobility(h)
    conductivity=interface.geometric_mobility(centers)/reference
    zero=1-mass
    value,residual,iterations=bruggeman_conductivity(conductivity,probabilities,zero,d)
    upper,_,_=bruggeman_conductivity(conductivity,probabilities,a,d,missing)
    return {'geometric':value*reference, 'upper':upper*reference,
            'factor':value*reference/h**power, 'mass':mass,'missing':missing,'zero':zero,
            'dimension':d,'residual':residual,'iterations':iterations}


def integrate_mobility(model, pressure_pa, interface, options=None, *, percolation_relative_contact_area=.42):
    """Shared conservative PDF integration for cubic liquid and quadratic gas laws."""
    options = options or FlowOptions()
    p=float(model._pressures(pressure_pa));a=model.area(p);h=model.mean_gap(p)
    if a<1. and h<=0:
        raise ValueError('Mean-gap underflow cannot be interpreted as a blocked flow state')
    if a==1.:
        return p,a,h,{'geometric':0.,'upper':0.,'factor':0.,'mass':0.,'missing':0.,
            'zero':1.,'dimension':effective_dimension(a,percolation_relative_contact_area),
            'residual':0.,'iterations':0},0,0.
    start=max(model.pdf_min_u1_m*.01,np.finfo(float).tiny)
    stop=max(40*model.psd.hrms,10*h,start*1e6)
    count=options.gap_points
    for refinement in range(options.max_gap_refinements+1):
        if count>65537:
            raise ValueError('Flow quadrature exceeded 65537 gap points')
        edges=np.r_[0.,np.geomspace(start,stop,count-1)]
        distribution=model.gap_distribution(p,edges)
        fine=homogenize_bins(edges,distribution.continuous_cdf,a,h,interface,
                            percolation_relative_contact_area)
        take=np.unique(np.r_[np.arange(0,count,2),count-1])
        coarse=homogenize_bins(edges[take],distribution.continuous_cdf[take],a,h,interface,
                              percolation_relative_contact_area)
        difference=abs(fine['geometric']-coarse['geometric'])/fine['geometric'] if fine['geometric'] else 0.
        if math.isfinite(fine['upper']) and fine['upper']>0:
            difference=max(difference,abs(fine['upper']-coarse['upper'])/fine['upper'])
        if difference<=options.quadrature_relative_tolerance:
            break
        count=2*count-1
    else:
        raise ValueError('Pressure-flow quadrature not converged; increase pressure_flow.gap_points')
    return p,a,h,fine,count,difference


def pressure_flow(model, pressure_pa, interface, options=None, *, percolation_relative_contact_area=.42):
    """Liquid K3 table. Ballistic callers must use ballistic_pressure_flow for K2 units."""
    if interface.gap_power!=3:
        raise ValueError('pressure_flow requires cubic liquid mobility; use ballistic_pressure_flow')
    options=options or FlowOptions()
    p,a,h,fine,count,difference=integrate_mobility(model,pressure_pa,interface,options,
        percolation_relative_contact_area=percolation_relative_contact_area)
    upper=fine['upper'];unbounded=not math.isfinite(upper)
    return PressureFlowResult(p,h,a,fine['dimension'],fine['geometric'],
        fine['geometric']/(12*interface.viscosity_pa_s),fine['factor'],
        None if unbounded else upper,None if unbounded else (upper/h**3 if h>0 else 0.),unbounded,
        fine['mass'],fine['missing'],fine['zero'],fine['residual'],fine['iterations'],
        count,difference,('blocked' if upper==0 else 'unresolved_connectivity')
        if fine['geometric']==0 else 'conducting',
        percolation_relative_contact_area,options.unresolved_flow_policy)
