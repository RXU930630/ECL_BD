"""Module 6: reference traction providers and local rough-contact unloading.

External FEM load is metadata. Only the explicitly supplied reference interface
pressure is added to reference contact pressure. No global compliance is assumed.
"""
from dataclasses import dataclass
from pathlib import Path
import numpy as np
from scipy.interpolate import PchipInterpolator,PPoly

from .interface import positive
from .profiles import PressureProfile,read_pressure_profile


@dataclass
class ReferenceMechanics:
    x_m: np.ndarray
    contact_reference_pa: np.ndarray
    interface_reference_pa: np.ndarray
    metadata: dict
    interpolation: str='pchip'

    def __post_init__(self):
        checked=PressureProfile(self.x_m,self.contact_reference_pa,self.interface_reference_pa,self.metadata)
        self.x_m=checked.x_m;self.contact_reference_pa=checked.contact_pressure_pa
        self.interface_reference_pa=checked.interface_reference_pressure_pa
        if self.interpolation not in ('pchip','linear'):raise ValueError('Profile interpolation must be pchip or linear')
        self._contact=PchipInterpolator(self.x_m,self.contact_reference_pa,extrapolate=False)
        self._fluid=PchipInterpolator(self.x_m,self.interface_reference_pa,extrapolate=False)

    @property
    def length_m(self):return float(self.x_m[-1]-self.x_m[0])

    def evaluate(self,x):
        x=np.asarray(x,dtype=float)
        if np.any(x<self.x_m[0]-1e-14*self.length_m) or np.any(x>self.x_m[-1]+1e-14*self.length_m):
            raise ValueError('Mechanics query outside pressure-profile domain')
        x=np.clip(x,self.x_m[0],self.x_m[-1])
        kind=self.metadata['type']
        if kind=='oring_hertz':
            center=(self.x_m[-1]+self.x_m[0])/2
            c=self.metadata['peak_pressure_pa']*np.sqrt(np.maximum(0,1-((x-center)/(self.length_m/2))**2))
            f=np.zeros_like(x)
        elif kind=='square_ring':
            c=np.full_like(x,self.metadata['reference_contact_pressure_pa']);f=np.zeros_like(x)
        elif kind=='gaussian':
            c=self.metadata['peak_pressure_pa']*np.exp(-.5*((x-self.length_m/2)/self.metadata['sigma_m'])**2);f=np.zeros_like(x)
        elif self.interpolation=='linear':
            c=np.interp(x,self.x_m,self.contact_reference_pa);f=np.interp(x,self.x_m,self.interface_reference_pa)
        else:c=self._contact(x);f=self._fluid(x)
        return c,f

    def total_reference(self,x):
        c,f=self.evaluate(x);return c+f

    @property
    def integration_points(self):
        kind=self.metadata['type'];start,end=self.x_m[0],self.x_m[-1]
        if kind=='prescribed_profile':
            points=self.x_m
            if self.interpolation=='pchip':
                poly=PPoly(self._contact.c+self._fluid.c,self.x_m,extrapolate=False)
                roots=poly.derivative().roots(extrapolate=False)
                roots=roots[np.isfinite(roots)&(roots>start)&(roots<end)]
                if len(roots):
                    nearest=np.clip(np.searchsorted(points,roots),0,len(points)-1)
                    distances=np.minimum(abs(roots-points[nearest]),abs(roots-points[np.maximum(nearest-1,0)]))
                    roots=roots[distances>32*np.finfo(float).eps*self.length_m]
                    points=np.unique(np.r_[points,roots])
            return points
        if kind=='gaussian':
            center=(start+end)/2;sigma=self.metadata['sigma_m']
            return np.unique(np.clip(center+sigma*np.array([-3,-2,-1,-.5,0,.5,1,2,3]),start,end))
        return np.array([start,(start+end)/2,end])

    def maximum_reference_pressure(self,include_interface=True):
        c,f=self.evaluate(self.integration_points)
        return float(max(c+f if include_interface else c))

    def contact_pressure(self,x,fluid_pressure_pa,coupling='local_unloading'):
        if coupling=='fixed_contact':return self.evaluate(x)[0]
        if coupling!='local_unloading':raise NotImplementedError('Nonlocal FEM compliance is not implemented')
        return np.maximum(self.total_reference(x)-np.asarray(fluid_pressure_pa),0.)


def make_mechanics(config,base_directory,effective_modulus_pa):
    if not isinstance(config,dict):raise ValueError('leakage.geometry must be a dictionary')
    kind=config.get('type')
    common={'type','interpolation','reference_external_fluid_pressure_pa'}
    fields={'prescribed_profile':{'pressure_file','file_format','profile_sampling'},
            'oring_hertz':{'length_m','peak_pressure_pa','cross_section_diameter_m'},
            'square_ring':{'length_m','reference_contact_pressure_pa'},
            'gaussian':{'length_m','peak_pressure_pa','sigma_m'}}
    if kind not in fields:raise ValueError('geometry.type must be prescribed_profile, oring_hertz, square_ring or gaussian')
    allowed=common|fields[kind]
    if set(config)-allowed:raise ValueError(f'Unknown parameters for {kind}: {sorted(set(config)-allowed)}')
    interpolation=config.get('interpolation','pchip')
    if kind=='prescribed_profile':
        if not isinstance(config.get('pressure_file'),str):raise ValueError('prescribed_profile requires pressure_file')
        path=Path(config['pressure_file'])
        if not path.is_absolute():path=Path(base_directory)/path
        data=read_pressure_profile(path,config.get('file_format','csv'),config.get('profile_sampling','preserve'))
        metadata={**data.metadata,'type':kind}
        if 'reference_external_fluid_pressure_pa' in config:
            metadata['reference_external_fluid_pressure_pa']=float(config['reference_external_fluid_pressure_pa'])
        return ReferenceMechanics(data.x_m,data.contact_pressure_pa,data.interface_reference_pressure_pa,metadata,interpolation)
    if kind not in ('oring_hertz','square_ring','gaussian'):
        raise ValueError('geometry.type must be prescribed_profile, oring_hertz, square_ring or gaussian')
    metadata={'type':kind,'reference_external_fluid_pressure_pa':config.get('reference_external_fluid_pressure_pa')}
    if kind=='oring_hertz':
        peak=positive(config.get('peak_pressure_pa'),'peak_pressure_pa')
        if ('length_m' in config)==('cross_section_diameter_m' in config):
            raise ValueError('Hertz geometry needs length_m OR cross_section_diameter_m')
        length=(positive(config['length_m'],'length_m') if 'length_m' in config else
                2*positive(config['cross_section_diameter_m'],'cross_section_diameter_m')*peak/effective_modulus_pa)
        metadata.update(peak_pressure_pa=peak,cross_section_diameter_m=config.get('cross_section_diameter_m'))
    else:
        length=positive(config.get('length_m'),'length_m')
        if kind=='square_ring':metadata['reference_contact_pressure_pa']=positive(config.get('reference_contact_pressure_pa'),'reference_contact_pressure_pa')
        else:metadata.update(peak_pressure_pa=positive(config.get('peak_pressure_pa'),'peak_pressure_pa'),sigma_m=positive(config.get('sigma_m'),'sigma_m'))
    x=np.linspace(0,length,101)
    if kind=='square_ring':p=np.full_like(x,metadata['reference_contact_pressure_pa'])
    elif kind=='oring_hertz':p=peak*np.sqrt(np.maximum(0,1-(2*x/length-1)**2))
    else:p=metadata['peak_pressure_pa']*np.exp(-.5*((x-length/2)/metadata['sigma_m'])**2)
    return ReferenceMechanics(x,p,np.zeros_like(x),metadata,interpolation)
