"""Module 4: bounded conductivity/contact lookup for spatial solvers."""
import numpy as np
from scipy.interpolate import PchipInterpolator


class TransportTable:
    def __init__(self,pressures,mean_gaps,conductivity,upper_conductivity,model,interface,
                 *,interpolation='log_pchip',percolation_policy='physical_zero',plateau_pressure_pa=None):
        self.p=np.array(pressures,dtype=float,copy=True);self.h=np.array(mean_gaps,dtype=float,copy=True)
        self.k=np.array(conductivity,dtype=float,copy=True);self.ku=np.array(upper_conductivity,dtype=float,copy=True)
        if (len(self.p)<3 or any(v.shape!=self.p.shape for v in (self.h,self.k,self.ku))
                or np.any(np.diff(self.p)<=0) or np.any(self.p<=0)
                or np.any(~np.isfinite(self.p)) or np.any(~np.isfinite(self.h))
                or np.any(self.h<=0) or np.any(~np.isfinite(self.k)) or np.any(self.k<0)
                or np.any(np.isnan(self.ku)) or np.any(self.ku<self.k)):
            raise ValueError('Invalid pressure/contact/transport database')
        if interpolation not in ('log_pchip','linear'):raise ValueError('Unsupported transport interpolation')
        if percolation_policy not in ('physical_zero','legacy_spc_plateau'):
            raise ValueError('percolation_policy must be physical_zero or legacy_spc_plateau')
        for array in (self.p,self.h,self.k,self.ku):array.flags.writeable=False
        self.model,self.interface=model,interface;self.interpolation=interpolation
        self.gap_power=interface.gap_power
        self.policy=percolation_policy;self.plateau=plateau_pressure_pa
        self._lookup_cache={}
        if self.policy=='legacy_spc_plateau':
            if self.plateau is None or not self.p[0]<=self.plateau<=self.p[-1]:
                raise ValueError('Transport database must cover the spc plateau pressure')
            if self._lookup(np.array([self.plateau]),self.k)[0]<=0:
                raise ValueError('spc plateau is already blocked; reduce spc or choose physical_zero')
        self._gap=PchipInterpolator(np.log(self.p),np.log(self.h),extrapolate=False)
        self.low_slope=float((self.h[0]-self.h[1])/np.log(self.p[1]/self.p[0]))
        if self.low_slope<=0:raise ValueError('Low-pressure mean-gap continuation requires decreasing h(p)')
        self.metadata={'pressure_min_pa':float(self.p[0]),'pressure_max_pa':float(self.p[-1]),
            'interpolation':interpolation,'percolation_policy':self.policy,'plateau_pressure_pa':self.plateau,
            'geometric_conductivity_unit':f'm^{self.gap_power}',
            'low_pressure_rule':'Matched logarithmic mean gap; conductance follows the selected local mobility ratio below database minimum.',
            'open_region_rule':'pc<=0: zero resistance reservoir limit, gap undefined; no finite aperture or global deformation calculated.'}

    def _lookup(self,p,values):
        left=np.clip(np.searchsorted(self.p,p,side='right')-1,0,len(self.p)-2)
        unbounded=np.isinf(values[left])|np.isinf(values[left+1])
        if self.interpolation=='linear':
            result=np.interp(p,self.p,np.where(np.isinf(values),0,values))
        else:
            result=np.zeros_like(p,dtype=float)
            key=id(values)
            cached=self._lookup_cache.get(key)
            if cached is None:
                finite=np.isfinite(values)&(values>0);pp=self.p[finite];vv=values[finite]
                fit=PchipInterpolator(np.log(pp),np.log(vv),extrapolate=False) if len(pp)>=2 else None
                cached=pp,vv,fit;self._lookup_cache[key]=cached
            pp,vv,fit=cached
            if fit is not None:
                inside=(p>=pp[0])&(p<=pp[-1])&~unbounded
                result[inside]=np.exp(fit(np.log(p[inside])))
            elif len(pp)==1:result[p==pp[0]]=vv[0]
            a,b=values[left],values[left+1]
            bridge=np.isfinite(a)&np.isfinite(b)&((a==0)!=(b==0))
            result[bridge]=a[bridge]+(b[bridge]-a[bridge])*(p[bridge]-self.p[left[bridge]])/(self.p[left[bridge]+1]-self.p[left[bridge]])
        result[unbounded]=np.inf
        exact_index=np.clip(np.searchsorted(self.p,p),0,len(self.p)-1)
        exact=p==self.p[exact_index]
        result[exact]=values[exact_index[exact]]
        return result

    def mean_gap(self,pressure):
        p=np.asarray(pressure,dtype=float)
        if np.any(~np.isfinite(p)) or np.any(p>self.p[-1]*(1+1e-12)):
            raise ValueError('Contact pressure exceeds transport database; extend the database')
        safe=np.clip(p,self.p[0],self.p[-1]);h=np.exp(self._gap(np.log(safe)))
        low=(p>0)&(p<self.p[0])
        h=np.where(low,self.h[0]+self.low_slope*(np.log(self.p[0])-np.log(np.maximum(p,np.finfo(float).tiny))),h)
        return np.where(p>0,h,np.nan)

    def geometric_conductivity(self,pressure,upper=False):
        p=np.asarray(pressure,dtype=float)
        if np.any(~np.isfinite(p)) or np.any(p>self.p[-1]*(1+1e-12)):
            raise ValueError('Contact pressure exceeds transport database; extend the database')
        query=np.clip(p,self.p[0],self.p[-1])
        if self.policy=='legacy_spc_plateau':query=np.minimum(query,self.plateau)
        values=self.ku if upper else self.k
        k=self._lookup(query.ravel(),values).reshape(p.shape)
        if upper:k=np.maximum(k,self._lookup(query.ravel(),self.k).reshape(p.shape))
        low=(p>0)&(p<self.p[0])
        if np.any(low):
            continued=self.mean_gap(np.maximum(p,np.finfo(float).tiny))
            k=np.where(low,values[0]*self.interface.geometric_mobility(continued)/
                       self.interface.geometric_mobility(self.h[0]),k)
        return np.where(p<=0,np.inf,k)

    def resistance(self,pressure,upper=False):
        k=self.geometric_conductivity(pressure,upper)
        r=np.full_like(k,np.inf)
        np.divide(self.interface.inverse_mobility_prefactor,k,out=r,where=k>0)
        return r
