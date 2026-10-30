"""Module 2+6 dry-contact coupling; the elastic operator stays independent.

The same mean gap satisfies geometry and the Persson constitutive lookup.
This is a dry load-controlled solve, not a fluid/solid mixed-lubrication solver.
"""
import numpy as np
from scipy.linalg import solve
from scipy.optimize import brentq
from scipy.interpolate import PchipInterpolator

from .macroscopic import MacroContactSolution,solve_smooth_contact


class RoughContactLaw:
    """Bounded module-2 lookup; tiny-pressure far tail is zero with a load bound.

    For h>h_max, monotonicity gives 0<p_true<p_min. No positive pressure
    floor or extrapolation is used. Callers must expose p_min * domain length
    as a fixed-gap omitted-load bound (not a global solution error bound).
    """
    def __init__(self,table):
        self.table=table
        self.h_min_m,self.h_max_m=table.valid_gap_range_m
        self._logp=PchipInterpolator(table.mean_gaps_m[::-1],np.log(table.pressures_pa[::-1]),extrapolate=False)
        self._derivative=self._logp.derivative()

    def pressure_and_derivative(self,h_m):
        h=np.asarray(h_m,dtype=float)
        if np.any(~np.isfinite(h)) or np.any(h<self.h_min_m):
            raise ValueError('Macro gap outside Persson table: extend macro_contact.contact_table pressure range')
        inside=h<=self.h_max_m
        p=np.zeros_like(h);derivative=np.zeros_like(h)
        p[inside]=np.exp(self._logp(h[inside]))
        derivative[inside]=p[inside]*self._derivative(h[inside])
        return p,derivative


def solve_rough_contact(elastic,hertz,law,*,tolerance=1e-9,max_iterations=80):
    """Newton solve in mean gaps and one rigid offset, at prescribed line load.

    h=h0+g+C p_c(h), dx sum(p_c)=W'. No fluid pressure or reference-pressure
    unloading is used. The outer-gap variables avoid a pressure floor in tails.
    """
    if not np.isfinite(tolerance) or tolerance<=0 or type(max_iterations) is not int or max_iterations<1:
        raise ValueError('Invalid rough-contact solver controls')
    smooth=solve_smooth_contact(elastic,hertz,tolerance=min(tolerance,1e-10))
    x=elastic.x_m;g=hertz.profile(x);p0=hertz.peak_pressure_pa
    scale=max(hertz.half_width_m**2/(2*hertz.radius_m),law.table.model.psd.hrms)
    c=elastic.compliance_m_pa*p0/scale;shape=g/scale
    lo,hi=law.h_min_m/scale,law.h_max_m/scale
    base=(g+smooth.displacement_m)/scale
    margin=1e-12*max(1.,hi)
    left=lo-float(np.min(base))+margin;right=hi-float(np.min(base))+margin
    mass_scale=elastic.dx_m*p0/hertz.line_load_n_m
    def initial_load(z):
        return float(np.sum(law.pressure_and_derivative((base+z)*scale)[0])*elastic.dx_m/hertz.line_load_n_m-1)
    if right<=left or initial_load(left)<0 or initial_load(right)>0:
        raise ValueError('Persson table cannot bracket macro load; extend its pressure range or check the domain')
    offset=brentq(initial_load,left,right,xtol=1e-13)
    y=base+offset
    def evaluate(h,z):
        pressure,derivative=law.pressure_and_derivative(h*scale)
        p=pressure/p0;dp=derivative*scale/p0
        residual=np.r_[h-shape-z-c@p,mass_scale*np.sum(p)-1]
        return residual,p,dp
    for iteration in range(1,max_iterations+1):
        residual,p,dp=evaluate(y,offset);norm=float(np.max(abs(residual)))
        if norm<=tolerance:break
        n=len(x);jac=np.empty((n+1,n+1))
        jac[:n,:n]=np.eye(n)-c*dp[None,:];jac[:n,n]=-1.
        jac[n,:n]=mass_scale*dp;jac[n,n]=0.
        step=solve(jac,-residual)
        alpha=1.
        downward=step[:-1]<0
        if np.any(downward):alpha=min(alpha,.99*float(np.min((y[downward]-lo)/-step[:-1][downward])))
        for _ in range(40):
            trial=y+alpha*step[:-1];z=offset+alpha*step[-1]
            candidate,_,_=evaluate(trial,z)
            if np.max(abs(candidate)) <= (1-1e-4*alpha)*norm:
                y,offset=trial,float(z);break
            alpha*=.5
        else:raise ValueError('Rough macro contact Newton line search failed; check table and resolution')
    else:raise ValueError('Rough macro contact did not converge')
    pressure=p*p0;u=elastic.displacement(pressure);h=y*scale
    geometric=float(np.max(abs(h-g-offset*scale-u)))/scale
    load_error=abs(elastic.line_load(pressure)/hertz.line_load_n_m-1)
    return MacroContactSolution(x.copy(),elastic.edges_m.copy(),pressure,
        np.asarray(law.table.model.area(pressure)),g,u,h,offset*scale,elastic.line_load(pressure),
        iteration,geometric,load_error,'rough_converged',float(np.min(h)),
        float(max(pressure[0],pressure[-1])/np.max(pressure)))
