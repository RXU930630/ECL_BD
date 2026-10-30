"""Steady 1D load-controlled elastic mixed lubrication with conservative cavitation.

Unknown mean gap, fluid pressure/saturation and rigid offset. Contact pressure
is the bounded Persson p_c(h) law; elasticity responds to p_c+p_f. No local
reference-pressure unloading, negative-pressure clipping or mobility floor.
"""
from dataclasses import dataclass

import numpy as np
from scipy.linalg import solve


def _newton_step(jacobian,residual):
    """Row/column equilibration; physical residual and equations stay unchanged."""
    row=np.max(abs(jacobian),axis=1)
    if np.any(row==0):raise ValueError('Newton matrix has a singular conservation row')
    scaled=jacobian/row[:,None]
    column=np.max(abs(scaled),axis=0)
    if np.any(column==0):raise ValueError('Newton matrix has an undetermined pressure/saturation column')
    try:return solve(scaled/column[None,:],-residual/row)/column
    except np.linalg.LinAlgError as exc:
        raise ValueError('Newton matrix is singular; check transport connectivity and continuation') from exc


@dataclass(frozen=True)
class FilmSolution:
    pressure_pa: np.ndarray
    saturation: np.ndarray
    face_flux_m2_s: np.ndarray
    iterations: int
    conservation_relative_residual: float


def film_flux(h,psi,dx,transport,pressure_scale,pin,pout):
    """Finite-volume flux/Jacobian with harmonic mobility and upwind saturation.

    psi>=0: p=P*psi, theta=1; -1<=psi<0: p=0, theta=1+psi.
    Boundary pressures are prescribed, reservoir saturation=1 at inflow.
    Returns exact analytic derivatives of this discrete residual.
    """
    h=np.asarray(h,float);psi=np.asarray(psi,float);n=len(h)
    m,dm,s,ds=transport(h)
    if np.any(m<=0):
        raise ValueError('Mixed solver currently requires connected positive mobility; a physical-zero blocked/trapped pocket needs a separate pocket model. No mobility floor was added.')
    full=psi>=0;p=pressure_scale*np.maximum(psi,0.);theta=1+np.minimum(psi,0.)
    dp=pressure_scale*full;dt=(~full).astype(float)
    mf=np.r_[m[0],2*m[:-1]*m[1:]/(m[:-1]+m[1:]),m[-1]]
    sf=np.r_[s[0],(s[:-1]+s[1:])/2,s[-1]]
    grad=np.r_[(p[0]-pin)/(dx/2),np.diff(p)/dx,(pout-p[-1])/(dx/2)]
    donor=np.where(sf>=0,np.arange(n+1)-1,np.arange(n+1))
    inside=(donor>=0)&(donor<n)
    th=np.ones(n+1);th[inside]=theta[donor[inside]]
    q=-mf*grad+sf*th
    jh=np.zeros((n+1,n));jp=np.zeros_like(jh)
    jh[0,0]=-dm[0]*grad[0]+ds[0]*th[0]
    jh[-1,-1]=-dm[-1]*grad[-1]+ds[-1]*th[-1]
    i=np.arange(n-1);den=(m[:-1]+m[1:])**2
    jh[i+1,i]=-2*m[1:]**2/den*dm[:-1]*grad[1:-1]+ds[:-1]*th[1:-1]/2
    jh[i+1,i+1]=-2*m[:-1]**2/den*dm[1:]*grad[1:-1]+ds[1:]*th[1:-1]/2
    jp[0,0]=-mf[0]*dp[0]/(dx/2);jp[-1,-1]=mf[-1]*dp[-1]/(dx/2)
    jp[i+1,i]=mf[1:-1]*dp[:-1]/dx;jp[i+1,i+1]=-mf[1:-1]*dp[1:]/dx
    faces=np.flatnonzero(inside)
    jp[faces,donor[inside]]+=sf[inside]*dt[donor[inside]]
    return q,jh,jp,p,theta


def solve_film(h,dx,transport,*,pressure_scale_pa,inlet_pressure_pa=0.,outlet_pressure_pa=0.,
               tolerance=1e-10,max_iterations=100):
    """Rigid prescribed-gap film; useful independently and for analytical checks."""
    if (np.ndim(h)!=1 or len(h)<3 or not np.isfinite(dx) or dx<=0 or
        not np.isfinite(pressure_scale_pa) or pressure_scale_pa<=0 or
        not np.isfinite(tolerance) or tolerance<=0 or type(max_iterations) is not int or max_iterations<1
        or not all(np.isfinite(v) and v>=0 for v in (inlet_pressure_pa,outlet_pressure_pa))):
        raise ValueError('Invalid film geometry, pressure boundaries or solver controls')
    h=np.asarray(h,float);n=len(h);psi=np.linspace(inlet_pressure_pa,outlet_pressure_pa,n)/pressure_scale_pa
    m,_,s,_=transport(h)
    scale=max(float(np.max(abs(s))),float(np.max(m))*pressure_scale_pa/(n*dx),np.finfo(float).tiny)
    for iteration in range(1,max_iterations+1):
        q,_,jp,p,theta=film_flux(h,psi,dx,transport,pressure_scale_pa,inlet_pressure_pa,outlet_pressure_pa)
        residual=np.diff(q)/scale;norm=float(np.max(abs(residual)))
        if norm<tolerance:break
        step=_newton_step(np.diff(jp,axis=0)/scale,residual)
        alpha=1.
        down=step<0
        if np.any(down):alpha=min(alpha,.995*float(np.min((psi[down]+1)/-step[down])))
        for _ in range(50):
            candidate=psi+alpha*step
            qt=film_flux(h,candidate,dx,transport,pressure_scale_pa,inlet_pressure_pa,outlet_pressure_pa)[0]
            if np.max(abs(np.diff(qt)/scale))<=(1-1e-4*alpha)*norm:
                psi=candidate;break
            alpha*=.5
        else:raise ValueError('Conservative film active-set line search failed')
    else:raise ValueError('Conservative film did not converge')
    return FilmSolution(p,theta,q,iteration,float(np.ptp(q)/max(float(np.max(abs(q))),np.finfo(float).tiny)))


@dataclass(frozen=True)
class MixedSolution:
    mean_gap_m: np.ndarray
    contact_pressure_pa: np.ndarray
    fluid_pressure_pa: np.ndarray
    saturation: np.ndarray
    displacement_m: np.ndarray
    rigid_offset_m: float
    face_flux_m2_s: np.ndarray
    iterations: int
    residual: float
    geometry_relative_residual: float
    load_relative_error: float
    mass_relative_residual: float


def solve_mixed(elastic,hertz,contact_law,transport,initial,*,inlet_pressure_pa=0.,outlet_pressure_pa=0.,
                tolerance=1e-8,max_iterations=100,trace=None):
    """Semismooth Newton, load control, p>=0, 0<=theta<=1, p(1-theta)=0.

    Gauge zero is both reservoir reference and cavitation pressure in this
    release; positive reservoir boundaries are allowed. Geometry is stationary
    in the cylinder-following frame; the upper cylinder has zero velocity.
    """
    if (not all(np.isfinite(v) and v>=0 for v in (inlet_pressure_pa,outlet_pressure_pa))
        or not np.isfinite(tolerance) or tolerance<=0 or type(max_iterations) is not int or max_iterations<1):
        raise ValueError('Invalid mixed boundaries or solver controls')
    n=len(elastic.x_m);p0=hertz.peak_pressure_pa
    hs=max(hertz.half_width_m**2/(2*hertz.radius_m),contact_law.table.model.psd.hrms)
    g=hertz.profile(elastic.x_m);c=elastic.compliance_m_pa*p0/hs
    mass_scale=elastic.dx_m*p0/hertz.line_load_n_m
    y=np.asarray(initial.mean_gap_m,float)/hs;offset=initial.rigid_offset_m/hs
    film=solve_film(y*hs,elastic.dx_m,transport,pressure_scale_pa=p0,
        inlet_pressure_pa=inlet_pressure_pa,outlet_pressure_pa=outlet_pressure_pa,tolerance=min(tolerance,1e-10))
    psi=film.pressure_pa/p0+film.saturation-1.
    m,_,s,_=transport(y*hs)
    qs=max(float(np.max(abs(s))),float(np.min(m))*p0/np.ptp(elastic.edges_m),np.finfo(float).tiny)
    # Characteristic Reynolds scale is held fixed during each nonlinear solve.
    def evaluate(v,z,w,jacobian=False):
        pc,dpc=contact_law.pressure_and_derivative(v*hs)
        q,jh,jp,pf,theta=film_flux(v*hs,w,elastic.dx_m,transport,p0,inlet_pressure_pa,outlet_pressure_pa)
        total=(pc+pf)/p0
        r=np.r_[v-g/hs-z-c@total,np.diff(q)/qs,mass_scale*np.sum(total)-1]
        if not jacobian:return r,(pc,pf,theta,q)
        dpc=dpc*hs/p0;dpf=(w>=0).astype(float)
        j=np.zeros((2*n+1,2*n+1))
        j[:n,:n]=np.eye(n)-c*dpc[None,:];j[:n,n:2*n]=-c*dpf[None,:];j[:n,-1]=-1
        j[n:2*n,:n]=np.diff(jh,axis=0)*hs/qs;j[n:2*n,n:2*n]=np.diff(jp,axis=0)/qs
        j[-1,:n]=mass_scale*dpc;j[-1,n:2*n]=mass_scale*dpf
        return r,j
    for iteration in range(1,max_iterations+1):
        residual,j=evaluate(y,offset,psi,True);norm=float(np.max(abs(residual)))
        if trace is not None:trace({'iteration':iteration,'residual':norm,'min_psi':float(min(psi)),
            'min_gap_m':float(min(y)*hs),'fluid_load_fraction':float(np.sum(np.maximum(psi,0))*mass_scale)})
        qcheck=film_flux(y*hs,psi,elastic.dx_m,transport,p0,inlet_pressure_pa,outlet_pressure_pa)[0]
        flux_spread=float(np.ptp(qcheck)/max(float(np.max(abs(qcheck))),np.finfo(float).tiny))
        if norm<=tolerance and flux_spread<=max(100*tolerance,1e-10):break
        step=_newton_step(j,residual);alpha=1.
        down=step[:n]<0
        if np.any(down):alpha=min(alpha,.995*float(np.min((y[down]-contact_law.h_min_m/hs)/-step[:n][down])))
        downp=step[n:2*n]<0
        if np.any(downp):alpha=min(alpha,.995*float(np.min((psi[downp]+1)/-step[n:2*n][downp])))
        for _ in range(50):
            ty=y+alpha*step[:n];tp=psi+alpha*step[n:2*n];tz=offset+alpha*step[-1]
            try:rn,_=evaluate(ty,tz,tp)
            except ValueError:rn=np.array([np.inf])
            if np.max(abs(rn))<=(1-1e-4*alpha)*norm:
                y,psi,offset=ty,tp,float(tz);break
            alpha*=.5
        else:raise ValueError('Mixed lubrication Newton line search failed; check connectivity, pressure table and continuation')
    else:raise ValueError('Mixed lubrication did not converge')
    residual,(pc,pf,theta,q)=evaluate(y,offset,psi)
    u=elastic.displacement(pc+pf)
    return MixedSolution(y*hs,pc,pf,theta,u,offset*hs,q,iteration,float(np.max(abs(residual))),
        float(np.max(abs(y*hs-g-offset*hs-u))/hs),
        abs(float(np.sum(pc+pf)*elastic.dx_m/hertz.line_load_n_m)-1),
        float(np.ptp(q)/max(float(np.max(abs(q))),np.finfo(float).tiny)))


def solve_mixed_continuation(elastic,hertz,contact_law,transport_factory,velocity,initial,*,
                             inlet_pressure_pa=0.,outlet_pressure_pa=0.,tolerance=1e-8,max_iterations=100):
    """Adaptive load-preserving continuation of velocity and boundary pressures.

    Try the requested state directly; on a nonlinear failure insert midpoint
    driving states. Every accepted final state solves the original equations;
    intermediate states are numerical initialization, not altered parameters.
    """
    accepted=[]
    def advance(start,factor0,factor1,depth):
        try:
            result=solve_mixed(elastic,hertz,contact_law,transport_factory(velocity*factor1),start,
                inlet_pressure_pa=inlet_pressure_pa*factor1,outlet_pressure_pa=outlet_pressure_pa*factor1,
                tolerance=tolerance,max_iterations=max_iterations)
            accepted.append({'driving_fraction':factor1,'iterations':result.iterations})
            return result
        except ValueError as exc:
            if depth>=12 or not ('Newton' in str(exc) or 'did not converge' in str(exc)):
                raise
            midpoint=(factor0+factor1)/2
            middle=advance(start,factor0,midpoint,depth+1)
            return advance(middle,midpoint,factor1,depth+1)
    return advance(initial,0.,1.,0),accepted
