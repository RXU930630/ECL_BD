"""Modules 8/9: conservative steady one-dimensional liquid-flow equations.

Unknowns are interior pressures and one constant volume flux. Newton solves
face pressure-drop = flux * face resistance; a second grid checks flux accuracy.
Geometry, contact and conductivity are supplied by the resistance callback.
"""
from dataclasses import dataclass
import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import spsolve


@dataclass
class SteadySolution:
    x_m: np.ndarray
    fluid_pressure_pa: np.ndarray
    flux_per_width_m2_s: float
    solver_nodes: int
    maximum_equation_residual: float
    pressure_drop_relative_error: float
    face_flux_relative_spread: float
    iterations: int
    flux_check_pressure_fraction: float
    grid_relative_flux_difference: float


def solve_steady_flow(resistance,x_start,x_end,inlet_pa,outlet_pa,*,nodes=201,
                      tolerance=1e-5,max_nodes=20001,output_points=1001,required_points=None):
    length=float(x_end-x_start);dp=float(inlet_pa-outlet_pa)
    if not np.all(np.isfinite([x_start,x_end,inlet_pa,outlet_pa,tolerance])) or tolerance<=0 or length<=0 or nodes<5 or output_points<5:raise ValueError('Invalid stationary flow grid')
    required=np.array([]) if required_points is None else np.asarray(required_points,dtype=float)
    if required.ndim!=1 or np.any(~np.isfinite(required)) or np.any(required<x_start) or np.any(required>x_end):
        raise ValueError('Required geometry mesh points must lie inside the domain')
    required=np.unique(np.r_[x_start,required,x_end])
    def sampling_grid(number):
        uniform=np.linspace(x_start,x_end,number)
        nearest=np.clip(np.searchsorted(required,uniform),0,len(required)-1)
        distance=np.minimum(abs(uniform-required[nearest]),abs(uniform-required[np.maximum(nearest-1,0)]))
        uniform=uniform[distance>32*np.finfo(float).eps*max(length,abs(x_start),abs(x_end))]
        return np.sort(np.r_[uniform,required])
    x=sampling_grid(nodes)
    if dp==0:
        xx=sampling_grid(output_points)
        return SteadySolution(xx,np.full_like(xx,inlet_pa),0.,len(x),0.,0.,0.,0,1.,0.)
    previous=None;previous_flux=None;total_iterations=0
    mesh_tolerance=max(10*tolerance,1e-6)
    while len(x)<=max_nodes:
        count=len(x);faces=(x[:-1]+x[1:])/2;dx=np.diff(x)/length
        y=1-(x-x_start)/length if previous is None else np.interp(x,previous[0],previous[1])
        m=count-1
        def raw_r(ymid):
            value=np.broadcast_to(np.asarray(resistance(faces,outlet_pa+dp*np.clip(ymid,0,1)),dtype=float),faces.shape)
            if np.any(~np.isfinite(value)) or np.any(value<0):
                raise ValueError('No finite resolved resistance: blocked branch requires separate treatment')
            return value
        rscale=float(max(raw_r((y[:-1]+y[1:])/2)))
        if rscale<=0:raise ValueError('Seal fully open: a structural/open-gap model is required for finite leakage')
        def evaluate(yn,flux):
            r=raw_r((yn[:-1]+yn[1:])/2)*dx/rscale
            return yn[:-1]-yn[1:]-flux*r,r
        q=(1/sum(raw_r((y[:-1]+y[1:])/2)*dx/rscale) if previous_flux is None
           else previous_flux*length*rscale/dp)
        for iteration in range(80):
            total_iterations+=1
            f,r=evaluate(y,q);norm=float(np.linalg.norm(f))
            if max(abs(f))<=max(tolerance*1e-3,1e-12)/m:break
            mid=(y[:-1]+y[1:])/2;step=1e-6
            dr=(raw_r(np.minimum(mid+step,1))-raw_r(np.maximum(mid-step,0)))
            dr/=np.maximum(np.minimum(mid+step,1)-np.maximum(mid-step,0),1e-12)
            dr*=dx/(2*rscale)
            left=1-q*dr;right=-1-q*dr
            rr=np.r_[np.arange(1,m),np.arange(m-1),np.arange(m)]
            cc=np.r_[np.arange(m-1),np.arange(m-1),np.full(m,m-1)]
            vv=np.r_[left[1:],right[:-1],-r]
            jac=coo_matrix((vv,(rr,cc)),shape=(m,m)).tocsc()
            update=spsolve(jac,-f)
            if np.any(~np.isfinite(update)):raise ValueError('Singular spatial leakage Jacobian')
            accepted=False
            for backtrack in range(35):
                damping=2.**(-backtrack)
                trial=y.copy();trial[1:-1]+=damping*update[:-1];trial_q=q+damping*update[-1]
                if trial_q<=0 or min(trial)<-1e-10 or max(trial)>1+1e-10 or max(np.diff(trial))>1e-10:continue
                trial_f,_=evaluate(trial,trial_q)
                if np.linalg.norm(trial_f)<norm:
                    y,q=trial,trial_q;accepted=True;break
            if not accepted:raise ValueError('Steady leakage Newton line search failed; no converged rate exported')
        else:raise ValueError('Steady leakage Newton iteration limit exceeded')
        flux=float(q*dp/(length*rscale))
        difference=abs(flux-previous_flux)/abs(flux) if previous_flux is not None else np.inf
        if previous_flux is not None and difference<=mesh_tolerance:break
        previous=(x,y);previous_flux=flux;x=np.unique(np.r_[x,faces])
    else:raise ValueError('Steady leakage spatial grid not converged; increase solver.max_nodes')
    f,r=evaluate(y,q)
    balance=abs(float(q*sum(r))-1)
    drops=y[:-1]-y[1:]
    active=(r>0)&(abs(drops)>1e-6)
    local=drops[active]/r[active]
    spread=float(max(abs(local/q-1))) if len(local) else 0.
    coverage=float(sum(abs(drops[active]))/sum(abs(drops)))
    if balance>max(20*tolerance,1e-6) or spread>max(100*tolerance,.002):
        raise ValueError(f'Leakage conservation check failed (drop={balance:g}, flux={spread:g})')
    xx=sampling_grid(output_points)
    pressure=outlet_pa+dp*np.interp(xx,x,np.clip(y,0,1))
    return SteadySolution(xx,pressure,flux,count,float(max(abs(f))),balance,spread,
                          total_iterations,coverage,float(difference))
