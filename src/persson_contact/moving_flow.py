"""Fixed-contact one-dimensional conserved pressure plus motion-driven transport.

For prescribed r(x) and s(x): dp/dx=r(x)[s(x)-q]. A signed constant q
and a nonmonotone pressure are allowed. No clipping or cavitation approximation.
"""
import numpy as np

from .steady import SteadySolution


def solve_fixed_contact_motion(resistance,source,x_start,x_end,inlet_pa,outlet_pa,*,
                               nodes=201,tolerance=1e-5,max_nodes=20001,
                               output_points=1001,required_points=None,pressure_floor_pa=0.):
    a,b,pin,pout=map(float,(x_start,x_end,inlet_pa,outlet_pa))
    if (not np.all(np.isfinite([a,b,pin,pout,tolerance,pressure_floor_pa])) or b<=a or
            nodes<5 or output_points<5 or max_nodes<nodes or tolerance<=0):
        raise ValueError('Invalid fixed-contact moving-flow grid')
    required=np.asarray([] if required_points is None else required_points,dtype=float)
    if required.ndim!=1 or np.any(~np.isfinite(required)) or np.any(required<a) or np.any(required>b):
        raise ValueError('Required motion mesh points must lie inside domain')
    required=np.unique(np.r_[a,required,b]);length=b-a
    def grid(n):
        uniform=np.linspace(a,b,n)
        nearest=np.clip(np.searchsorted(required,uniform),0,len(required)-1)
        distance=np.minimum(abs(uniform-required[nearest]),abs(uniform-required[np.maximum(nearest-1,0)]))
        return np.sort(np.r_[uniform[distance>32*np.finfo(float).eps*max(length,abs(a),abs(b))],required])
    x=grid(nodes);previous=None;iterations=0
    while len(x)<=max_nodes:
        iterations+=1
        face=(x[:-1]+x[1:])/2;dx=np.diff(x)
        r=np.broadcast_to(np.asarray(resistance(face),dtype=float),face.shape)
        s=np.broadcast_to(np.asarray(source(face),dtype=float),face.shape)
        if np.any(~np.isfinite(r)) or np.any(r<=0) or np.any(~np.isfinite(s)):
            raise ValueError('Moving-wall flow requires finite positive mobility and a finite gap in every interior face')
        # Normalize resistance before summation; no arbitrary conductivity floor.
        weight=(r/np.max(r))*dx
        total=float(np.sum(weight));scale_r=float(np.max(r))
        q_pressure=(pin-pout)/(scale_r*total)
        q=q_pressure+float(np.sum(weight*s))/total
        drops=scale_r*weight*(q-s)
        pressure=pin-np.r_[0.,np.cumsum(drops)]
        scale_q=max(abs(q),abs(q_pressure),float(np.max(abs(s))),np.finfo(float).tiny)
        scale_p=max(abs(pin),abs(pout),float(np.max(abs(pressure))),float(np.sum(abs(drops))),1.)
        if previous is not None:
            difference=abs(q-previous[2])/scale_q
            pdiff=float(np.max(abs(pressure-np.interp(x,previous[0],previous[1]))))/scale_p
            if max(difference,pdiff)<=max(10*tolerance,1e-6):break
        previous=(x,pressure,q);x=np.unique(np.r_[x,face])
    else:raise ValueError('Moving-flow spatial grid not converged; increase solver.max_nodes')
    closure=abs(float(pressure[-1])-pout)/scale_p
    face_q=(pressure[:-1]-pressure[1:])/(r*dx)+s
    spread=float(np.max(abs(face_q-q)))/scale_q
    residual=float(np.max(abs((pressure[:-1]-pressure[1:])-drops)))/scale_p
    if max(closure,spread)>max(100*tolerance,1e-6):
        raise ValueError('Moving-flow conservation check failed')
    floor_roundoff=max(128*np.finfo(float).eps*max(scale_p,abs(pressure_floor_pa)),1e-10)
    if np.min(pressure)<pressure_floor_pa-floor_roundoff:
        raise ValueError('Fully flooded moving-flow pressure falls below pressure_floor_pa; cavitation is not implemented. Pressure was not clipped.')
    xx=np.unique(np.r_[grid(output_points),x[np.r_[True,np.diff(np.sign(np.diff(pressure)))!=0,True]]])
    # Preserve endpoints exactly; roundoff is checked above, not used as a floor.
    pressure[-1]=pout
    output=np.interp(xx,x,pressure)
    return SteadySolution(xx,output,float(q),len(x),residual,closure,spread,iterations,1.,float(max(difference,pdiff)))
