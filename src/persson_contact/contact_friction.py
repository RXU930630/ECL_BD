"""Module 5: sliding shear strength on real solid-contact area.

Magnitude laws are independent of direction. The upper traction opposes
V=U_upper-U_lower, the lower traction is opposite. No static-friction solve,
fluid friction, viscoelastic hysteresis, heat solve or tangential deformation.
"""
from dataclasses import asdict,dataclass
import hashlib
import importlib.util
import json
import sys
from pathlib import Path
from typing import Protocol

import numpy as np


def finite(value,name,*,minimum=None,strict=False):
    if isinstance(value,(bool,np.bool_)):
        raise ValueError(f'{name} must be a finite number')
    try:value=float(value)
    except (TypeError,ValueError) as exc:raise ValueError(f'Provide {name}; no experimental value is assumed') from exc
    if not np.isfinite(value) or (minimum is not None and (value<minimum or strict and value==minimum)):
        raise ValueError(f'Invalid {name}')
    return value


@dataclass(frozen=True)
class ShearState:
    """Arrays on real-contact cells only; pressures in Pa, T in K, speed in m/s.

    nominal_contact_pressure_pa is p_c; real_contact_pressure_pa is p_c/a,
    a mean within contact spots, not the pressure at each microscopic point.
    """
    x_m: np.ndarray
    relative_velocity_m_s: np.ndarray
    nominal_contact_pressure_pa: np.ndarray
    real_contact_pressure_pa: np.ndarray
    temperature_k: np.ndarray
    relative_contact_area: np.ndarray
    mean_gap_m: np.ndarray


class ShearStrengthLaw(Protocol):
    def strength_pa(self,state:ShearState): ...
    @property
    def metadata(self): ...


@dataclass(frozen=True)
class ConstantShearStrength:
    shear_strength_pa: float

    def __post_init__(self):
        object.__setattr__(self,'shear_strength_pa',finite(self.shear_strength_pa,'shear_strength_pa',minimum=0))

    def strength_pa(self,state):return self.shear_strength_pa

    @property
    def metadata(self):return {'model':'constant_shear_strength',**asdict(self),'stress_basis':'real_contact_area'}


@dataclass(frozen=True)
class LogGaussianShearStrength:
    baseline_pa: float
    amplitude_pa: float
    peak_velocity_m_s: float
    width_decades: float

    def __post_init__(self):
        for key in ('baseline_pa','amplitude_pa'):
            object.__setattr__(self,key,finite(getattr(self,key),key,minimum=0))
        object.__setattr__(self,'peak_velocity_m_s',finite(self.peak_velocity_m_s,'peak_velocity_m_s',
            minimum=0,strict=True))
        object.__setattr__(self,'width_decades',finite(self.width_decades,'width_decades',minimum=0,strict=True))

    def strength_pa(self,state):
        speed=abs(state.relative_velocity_m_s)
        z=np.full_like(speed,np.inf,dtype=float)
        positive=speed>0
        # Difference of logs avoids an overflowing velocity ratio.
        z[positive]=(np.log10(speed[positive])-np.log10(self.peak_velocity_m_s))/self.width_decades
        with np.errstate(over='ignore',under='ignore'):
            return self.baseline_pa+self.amplitude_pa*np.exp(-.5*z*z)

    @property
    def metadata(self):
        return {'model':'gaussian_log_velocity',**asdict(self),'width_units':'decades',
            'stress_basis':'real_contact_area',
            'zero_velocity_limit':'baseline_pa; not static friction'}


class CustomShearStrength:
    """User callback(state, parameters) returns scalar/vector strength magnitude [Pa]."""
    def __init__(self,callback,parameters=None,*,provenance=None):
        if not callable(callback):raise ValueError('Custom shear law must be callable')
        if parameters is None:parameters={}
        if not isinstance(parameters,dict):raise ValueError('Custom parameters must be a dictionary')
        # A copy also enforces a reproducible parameter snapshot.
        try:self.parameters=json.loads(json.dumps(parameters,allow_nan=False))
        except (TypeError,ValueError) as exc:raise ValueError('Custom parameters must be finite JSON values') from exc
        self.callback=callback;self.provenance=provenance or {}

    def strength_pa(self,state):return self.callback(state,self.parameters)

    @property
    def metadata(self):
        return {'model':'custom','parameters':self.parameters,'stress_basis':'real_contact_area',**self.provenance}


def make_shear_strength(raw,base_directory=Path('.')):
    if not isinstance(raw,dict):raise ValueError('contact_friction must be a dictionary')
    name=raw.get('model','constant_shear_strength')
    presets=raw.get('model_parameters',{})
    if not isinstance(presets,dict):raise ValueError('model_parameters must be a dictionary')
    selected=presets.get(name,{})
    if not isinstance(selected,dict):raise ValueError('Selected shear parameters must be a dictionary')
    fields={'constant_shear_strength':{'shear_strength_pa'},
        'gaussian_log_velocity':{'baseline_pa','amplitude_pa','peak_velocity_m_s','width_decades'},
        'custom':{'module_file','function','parameters'}}
    if name not in fields:raise NotImplementedError(f'Unknown contact shear-strength model: {name}')
    if set(selected)!=fields[name]:raise ValueError(f'{name} requires exactly {sorted(fields[name])}')
    if name=='constant_shear_strength':return ConstantShearStrength(selected['shear_strength_pa'])
    if name=='gaussian_log_velocity':
        return LogGaussianShearStrength(selected['baseline_pa'],selected['amplitude_pa'],selected['peak_velocity_m_s'],selected['width_decades'])
    module_file=selected['module_file'];function=selected['function']
    if not isinstance(module_file,str) or not module_file or not isinstance(function,str) or not function.isidentifier():
        raise ValueError('Provide a custom module_file and a function name')
    path=Path(module_file);path=path if path.is_absolute() else Path(base_directory)/path
    path=path.resolve()
    if path.suffix!='.py' or not path.is_file():raise ValueError('Custom module_file must be an existing Python file')
    source=path.read_bytes();digest=hashlib.sha256(source).hexdigest()
    spec=importlib.util.spec_from_file_location('contact_shear_'+digest,path)
    module=importlib.util.module_from_spec(spec)
    previous=sys.modules.get(spec.name);sys.modules[spec.name]=module
    try:
        # Execute the exact hashed source; avoid timestamp-based stale pyc reuse.
        exec(compile(source,str(path),'exec'),module.__dict__)
    except Exception:
        if previous is None:sys.modules.pop(spec.name,None)
        else:sys.modules[spec.name]=previous
        raise
    callback=getattr(module,function,None)
    return CustomShearStrength(callback,selected['parameters'],provenance={
        'module_file':str(path),'module_sha256':digest,'function':function})


@dataclass
class ContactFrictionResult:
    summary: dict
    profile: dict


def evaluate_contact_friction(cell_edges_m,relative_contact_area,nominal_contact_pressure_pa,mean_gap_m,
                              relative_velocity_m_s,law:ShearStrengthLaw,*,temperature_k=298.15,transverse_width_m=None):
    edges=np.asarray(cell_edges_m,dtype=float)
    if edges.ndim!=1 or len(edges)<2 or np.any(~np.isfinite(edges)) or np.any(np.diff(edges)<=0):
        raise ValueError('Contact friction requires increasing finite cell edges')
    shape=(len(edges)-1,)
    arrays=[]
    for raw,name in [(relative_contact_area,'relative_contact_area'),(nominal_contact_pressure_pa,'contact pressure'),
                     (mean_gap_m,'mean gap'),(temperature_k,'temperature')]:
        try:value=np.broadcast_to(np.asarray(raw,dtype=float),shape).copy()
        except (TypeError,ValueError) as exc:raise ValueError(f'Invalid {name} shape') from exc
        if np.any(~np.isfinite(value)):raise ValueError(f'{name} must be finite')
        arrays.append(value)
    area,pressure,gap,temp=arrays
    if np.any((area<0)|(area>1)) or np.any(pressure<0) or np.any(temp<=0):
        raise ValueError('Require 0<=a<=1, nonnegative pressure and positive Kelvin temperature')
    gap_roundoff=128*np.finfo(float).eps*max(float(np.ptp(edges)),float(np.max(abs(gap))),np.finfo(float).tiny)
    if np.any(gap < -gap_roundoff):raise ValueError('Mean gap is negative beyond numerical roundoff')
    if np.any((area==0)&(pressure>0)) or np.any((area>0)&(pressure==0)):
        raise ValueError('Inconsistent nonadhesive contact pressure and relative contact area')
    velocity=finite(relative_velocity_m_s,'relative_velocity_m_s')
    width=None if transverse_width_m is None else finite(transverse_width_m,'transverse_width_m',minimum=0,strict=True)
    dx=np.diff(edges);x=(edges[:-1]+edges[1:])/2;active=area>0
    real_pressure=np.zeros(shape);np.divide(pressure,area,out=real_pressure,where=active)
    if np.any(~np.isfinite(real_pressure)):raise ValueError('Real-contact mean pressure overflows')
    strength=np.zeros(shape)
    if np.any(active):
        state=ShearState(x[active],np.full(np.sum(active),velocity),pressure[active],real_pressure[active],
            temp[active],area[active],gap[active])
        try:values=np.broadcast_to(np.asarray(law.strength_pa(state),dtype=float),(int(np.sum(active)),))
        except (TypeError,ValueError) as exc:raise ValueError('Shear callback must return scalar or one value per real-contact cell') from exc
        if np.any(~np.isfinite(values)) or np.any(values<0):raise ValueError('Shear strength must be finite and nonnegative')
        strength[active]=values
    nominal_capacity=area*strength
    upper=-np.sign(velocity)*nominal_capacity;lower=-upper
    magnitude=float(dx@abs(upper));upper_force=float(dx@upper);lower_force=-upper_force
    normal=float(dx@pressure);real_area=float(dx@area)
    moving=velocity!=0
    coefficient=magnitude/normal if moving and normal>0 else None
    power=magnitude*abs(velocity)
    def total(value):return None if width is None else value*width
    quantities=[normal,real_area,magnitude,upper_force,lower_force,power]
    if width is not None:quantities += [q*width for q in quantities]
    if not all(np.isfinite(q) for q in quantities):raise ValueError('Contact friction force/power integral overflow')
    status='no_real_contact' if real_area==0 else 'sliding' if moving else 'no_relative_sliding_static_unmodelled'
    summary={'relative_velocity_m_s':velocity,'speed_m_s':abs(velocity),'status':status,
        'normal_line_load_n_m':normal,'real_contact_area_per_width_m':real_area,
        'friction_force_magnitude_per_width_n_m':magnitude,
        'upper_friction_force_per_width_n_m':upper_force,'lower_friction_force_per_width_n_m':lower_force,
        'kinetic_contact_friction_coefficient':coefficient,'dissipated_power_per_width_w_m':power,
        'transverse_width_m':width,'normal_force_n':total(normal),'real_contact_area_m2':total(real_area),
        'friction_force_magnitude_n':total(magnitude),'upper_friction_force_n':total(upper_force),
        'lower_friction_force_n':total(lower_force),'dissipated_power_w':total(power),
        'static_friction_limit_n_m':None}
    profile={'x_m':x,'cell_left_m':edges[:-1],'cell_right_m':edges[1:],
        'relative_contact_area':area,'nominal_contact_pressure_pa':pressure,'real_contact_pressure_pa':real_pressure,
        'mean_gap_m':gap,'temperature_k':temp,'relative_velocity_m_s':np.full(shape,velocity),
        'real_contact_shear_strength_pa':strength,'nominal_sliding_strength_pa':nominal_capacity,
        'upper_contact_traction_pa':upper,'lower_contact_traction_pa':lower,
        'dissipated_power_per_nominal_area_w_m2':abs(upper)*abs(velocity)}
    return ContactFrictionResult(summary,profile)
