"""Module 1: pressure-profile readers. Coordinates are never silently rebuilt."""
from dataclasses import dataclass
import hashlib
from pathlib import Path
import numpy as np


@dataclass(frozen=True)
class PressureProfile:
    x_m: np.ndarray
    contact_pressure_pa: np.ndarray
    interface_reference_pressure_pa: np.ndarray
    metadata: dict

    def __post_init__(self):
        x,p,f=(np.asarray(v,dtype=float).copy() for v in
               (self.x_m,self.contact_pressure_pa,self.interface_reference_pressure_pa))
        if (x.ndim!=1 or len(x)<2 or x.shape!=p.shape or x.shape!=f.shape
                or any(np.any(~np.isfinite(v)) for v in (x,p,f))
                or np.any(np.diff(x)<=0) or np.any(p<0) or max(p)<=0):
            raise ValueError('Pressure profile needs >=2 increasing x coordinates and finite nonnegative contact pressures')
        for name,value in (('x_m',x),('contact_pressure_pa',p),('interface_reference_pressure_pa',f)):
            value.flags.writeable=False;object.__setattr__(self,name,value)


def read_pressure_profile(path, format='csv', sampling='preserve'):
    path=Path(path)
    text=path.read_text(encoding='utf-8-sig')
    metadata={'path':str(path.resolve()),'sha256':hashlib.sha256(path.read_bytes()).hexdigest(),
              'format':format,'sampling':sampling,'reference_external_fluid_pressure_pa':None}
    if format=='csv':
        data=np.genfromtxt(path,delimiter=',',names=True,encoding='utf-8-sig')
        names=data.dtype.names or ()
        if set(names)-{'x_m','contact_pressure_pa','interface_reference_pressure_pa'} or not {'x_m','contact_pressure_pa'}<=set(names):
            raise ValueError('CSV columns: x_m,contact_pressure_pa[,interface_reference_pressure_pa]')
        x=np.atleast_1d(data['x_m']);p=np.atleast_1d(data['contact_pressure_pa'])
        f=np.atleast_1d(data['interface_reference_pressure_pa']) if 'interface_reference_pressure_pa' in names else np.zeros_like(p)
    elif format=='fortran':
        lines=[line.strip() for line in text.splitlines() if line.strip()]
        if len(lines)<4:raise ValueError('Incomplete Fortran profile')
        metadata['name']=lines[0]
        metadata['reference_external_fluid_pressure_pa']=float(lines[1].split()[0].replace('D','E'))*1e6
        data=np.array([[float(v.replace('D','E').replace('d','e')) for v in line.split()] for line in lines[2:]])
        if data.ndim!=2 or data.shape[1]!=2:raise ValueError('Fortran profile requires x [m], contact pressure [MPa]')
        x=data[:,0];p=data[:,1]*1e6;f=np.zeros_like(p)
    else:raise ValueError('profile format must be csv or fortran')
    profile=PressureProfile(x,p,f,metadata)
    if sampling=='legacy_uniform_zero_padding':
        n=len(x);step=(x[-1]-x[0])/n
        metadata={**metadata,'original_points':n,'original_length_m':float(x[-1]-x[0]),
                  'warning':'Explicit Fortran compatibility: input positions replaced, zero endpoints added; domain length changes.'}
        return PressureProfile(x[0]+np.arange(n+2)*step,np.r_[0,p,0],np.r_[0,f,0],metadata)
    if sampling!='preserve':raise ValueError('profile_sampling must be preserve or legacy_uniform_zero_padding')
    return profile
