"""Persson shear friction factors, distinct from shear/pressure FLOW factors.

Isotropic, one rough wall, constant-Newtonian liquid, no slip. The interpolation
kernel follows the original September-2026 archive. How <1/u> is obtained is an
explicit closure choice; it is not an integral through a broadened PDF at u=0.
"""
from dataclasses import dataclass

import numpy as np

from .interface import positive
from .shear import _finite


@dataclass(frozen=True)
class FrictionFactors:
    phi_f: object
    phi_fs: object
    phi_fp: object


def shear_friction_factors(mean_gap_m, inverse_gap_m_inv, rms_height_m, percolation_gap_m):
    """Original no-slip interpolation; forfluidflow=3*hrms^2/2 (isotropy).

    phi_f=h I; phi_fp=h(h-hc)/(h^2+3 sigma^2/2) above hc, else 0;
    phi_fs=3 h I/[1+2(h-hc)h^2 I/sigma^2] above hc, else 3 h I.
    I is a supplied nominal-area inverse-gap moment, not 1/mean_gap.
    Smooth limit is explicit: (1,0,1). No original phi_f>=1 clamp or
    pressure-index fit is applied to the supplied inverse moment.
    """
    h,inv=np.broadcast_arrays(np.asarray(mean_gap_m,float),np.asarray(inverse_gap_m_inv,float))
    rms=_finite(rms_height_m,'rms_height_m');hc=_finite(percolation_gap_m,'percolation_gap_m')
    if (np.any(~np.isfinite(h)) or np.any(h<=0) or np.any(~np.isfinite(inv))
            or np.any(inv<0) or rms<0 or hc<0):
        raise ValueError('Require positive gaps and nonnegative finite inverse moment, RMS and hc')
    if rms==0:
        return FrictionFactors(np.ones_like(h),np.zeros_like(h),np.ones_like(h))
    delta=np.maximum(h-hc,0.)
    f=h*inv
    fs=3*f/(1+2*(delta/rms)*(h/rms)*f)
    fp=(delta/h)/(1+1.5*(rms/h)**2)
    return FrictionFactors(f,fs,fp)


def scale_inverse_gap(model,pressure_pa,minimum_shear_gap_m=0.):
    """I=sum(-Delta A/u1), conservative Persson magnification-scale closure.

    This is the archive's iuinverse=1 TYPE of approximation, evaluated using
    the improved conservative scale cells. It is not the broadened gap PDF's
    inverse moment (that PDF has positive density at zero and a divergent
    unregularized inverse moment). Existing area/u1 cutoffs remain explicit.
    Missing noncontact probability is neither contact nor renormalized away.
    """
    if model.pdf_scale_method!='conservative':
        raise ValueError('Shear friction scale_u1 requires options.pdf_scale_method=conservative')
    p=positive(pressure_pa,'pressure_pa')
    cutoff=_finite(minimum_shear_gap_m,'minimum_shear_gap_m')
    if cutoff<0:raise ValueError('minimum_shear_gap_m must be nonnegative')
    centers,weights,_,high,low,small=model._conservative_pdf_scales(p)
    keep=centers>=cutoff
    moment=lambda threshold:float(np.sum(weights[centers>=threshold]/centers[centers>=threshold]))
    return {'inverse_gap_m_inv':moment(cutoff),
        'uncut_inverse_gap_m_inv':moment(0.),
        'half_cutoff_inverse_gap_m_inv':moment(cutoff/2),
        'double_cutoff_inverse_gap_m_inv':moment(2*cutoff),
        'excluded_shear_probability':float(np.sum(weights[~keep])),
        'retained_probability':float(np.sum(weights)),
        'unrepresented_noncontact_probability':float(high+low+small),
        'excluded_small_u1_probability':float(small)}


def wall_tractions(mean_gap_m, factors, viscosity_pa_s, pressure_gradient_pa_m,
                   lower_velocity_m_s, upper_velocity_m_s, *, rough_surface='upper', saturation=1.):
    """Fluid-on-wall +x tractions per NOMINAL area, not internal fluid stresses.

    Couette: t0=eta*(U1-U0)*(phi_f+phi_fs)/h, t1=-t0.
    Pressure: smooth wall coefficient phi_fp; rough wall 2-phi_fp.
    Therefore t0+t1=-h*dp/dx. Rough-wall result includes the microscopic
    pressure/form-drag contribution of this averaged closure. Macroscopic
    surface slope pressure traction is separate. In cavitated cells theta
    scales the sliding term (uniform wetted-striation approximation); it is
    NOT a derived rough-surface cavitation friction theory.
    """
    if rough_surface not in ('lower','upper'):
        raise ValueError('rough_surface must be lower or upper')
    eta=positive(viscosity_pa_s,'viscosity_pa_s')
    h,g,u0,u1,theta,f,fs,fp=np.broadcast_arrays(np.asarray(mean_gap_m,float),
        np.asarray(pressure_gradient_pa_m,float),np.asarray(lower_velocity_m_s,float),
        np.asarray(upper_velocity_m_s,float),np.asarray(saturation,float),
        np.asarray(factors.phi_f,float),np.asarray(factors.phi_fs,float),np.asarray(factors.phi_fp,float))
    if (any(np.any(~np.isfinite(v)) for v in (h,g,u0,u1,theta,f,fs,fp)) or
        np.any(h<=0) or np.any((theta<0)|(theta>1)) or np.any(f<0) or np.any(fs<0)
        or np.any((fp<0)|(fp>1))):
        raise ValueError('Invalid wall-traction state/factors')
    sliding=theta*eta*(u1-u0)*(f+fs)/h
    c0=fp if rough_surface=='upper' else 2-fp
    c1=2-fp if rough_surface=='upper' else fp
    p0=-c0*h*g/2;p1=-c1*h*g/2
    return {'lower_sliding_pa':sliding,'upper_sliding_pa':-sliding,
        'lower_pressure_flow_pa':p0,'upper_pressure_flow_pa':p1,
        'lower_total_pa':sliding+p0,'upper_total_pa':-sliding+p1}
