"""Measured finite surface geometry, with the existing half-space compliance.

Finite surface extent is NOT a finite-body elasticity model. Thickness is
recorded only; backing, bonding constraints and prestress remain unmodelled.
"""
from dataclasses import dataclass
import math

import numpy as np

from .macroscopic import HertzLine,_positive


@dataclass(frozen=True)
class FiniteCylindricalPad(HertzLine):
    axial_length_m: float
    sliding_length_m: float
    length_definition: str = 'projected'
    rubber_thickness_m: float | None = None
    physical_boundaries = True

    def __post_init__(self):
        super().__post_init__()
        for name in ('axial_length_m','sliding_length_m'):
            object.__setattr__(self,name,_positive(getattr(self,name),name))
        if self.length_definition not in ('projected','arc_length'):
            raise ValueError('length_definition must be projected or arc_length')
        if self.length_definition=='projected' and self.sliding_length_m>=2*self.radius_m:
            raise ValueError('Projected sliding length must be smaller than the cylinder diameter')
        if self.length_definition=='arc_length' and self.sliding_length_m>=math.pi*self.radius_m:
            raise ValueError('Arc length must describe less than a semicircle')
        if self.rubber_thickness_m is not None:
            object.__setattr__(self,'rubber_thickness_m',_positive(self.rubber_thickness_m,'rubber_thickness_m'))

    @property
    def domain_half_width_m(self):
        return self.sliding_length_m/2 if self.length_definition=='projected' else self.radius_m*math.sin(self.sliding_length_m/(2*self.radius_m))

    def profile(self,x_m):
        x=np.asarray(x_m,float)
        if np.any(~np.isfinite(x)) or np.any(abs(x)>self.domain_half_width_m*(1+1e-12)):
            raise ValueError('Coordinate lies outside the actual pad surface; geometry is not extended')
        # Stable circular sag: avoids subtracting nearly equal radii near x=0.
        return x*x/(self.radius_m+np.sqrt((self.radius_m-x)*(self.radius_m+x)))

    @property
    def metadata(self):
        return {'type':'finite_cylindrical_pad','outer_radius_m':self.radius_m,
            'axial_length_m':self.axial_length_m,'sliding_length_m':self.sliding_length_m,
            'length_definition':self.length_definition,'rubber_thickness_m_record_only':self.rubber_thickness_m,
            'projected_sliding_length_m':2*self.domain_half_width_m,
            'surface_interval_m':[-self.domain_half_width_m,self.domain_half_width_m],
            'profile':'circular arc; x^2/(R+sqrt(R^2-x^2))',
            'boundary_interpretation':'physical surface edges, not adjustable far-field truncation',
            'elastic_approximation':'unchanged plane-strain half-space compliance; finite thickness, axial end effects, backing/bonding and prestress ignored'}


def finite_geometry_parameters(parameters):
    raw=parameters.get('macro_contact',{})
    if raw.get('geometry')!='finite_cylindrical_pad':return None
    presets=raw.get('geometry_parameters',{})
    if not isinstance(presets,dict) or not isinstance(presets.get('finite_cylindrical_pad'),dict):
        raise ValueError('finite_cylindrical_pad requires its geometry_parameters block')
    values=presets['finite_cylindrical_pad']
    allowed={'outer_radius_m','axial_length_m','sliding_length_m','length_definition','rubber_thickness_m'}
    if set(values)-allowed:raise ValueError('Unknown finite_cylindrical_pad parameter')
    # Validate shape independently of material/load; no files or solver work.
    return FiniteCylindricalPad(values.get('outer_radius_m'),1.,1.,values.get('axial_length_m'),
        values.get('sliding_length_m'),values.get('length_definition','projected'),values.get('rubber_thickness_m'))


def transverse_width(parameters):
    """Single source: finite geometry axial length; legacy width remains valid."""
    shape=finite_geometry_parameters(parameters)
    legacy=parameters.get('contact_friction',{}).get('transverse_width_m')
    if legacy is not None:legacy=_positive(legacy,'transverse_width_m')
    if shape is None:return legacy
    if legacy is not None and not math.isclose(legacy,shape.axial_length_m,rel_tol=1e-12):
        raise ValueError('transverse_width_m conflicts with geometry.axial_length_m; leave the legacy width None')
    return shape.axial_length_m


def line_load(parameters):
    raw=parameters.get('macro_contact',{})
    force=raw.get('normal_force_n');linear=raw.get('line_load_n_m')
    if force is not None:
        if linear is not None:raise ValueError('Give normal_force_n OR line_load_n_m, not both')
        width=transverse_width(parameters)
        if width is None:raise ValueError('normal_force_n requires geometry axial length or transverse_width_m')
        return _positive(force,'normal_force_n')/width
    return _positive(linear,'line_load_n_m')


def friction_parameters(parameters):
    return {**parameters.get('contact_friction',{}),'transverse_width_m':transverse_width(parameters)}


def geometry_metadata(shape,elastic):
    if getattr(shape,'physical_boundaries',False):return shape.metadata
    return {'type':'cylinder_plane','radius_m':shape.radius_m,'profile':'parabolic x^2/(2R)',
        'surface_interval_m':[float(elastic.edges_m[0]),float(elastic.edges_m[-1])],
        'boundary_interpretation':'numerical far-field truncation',
        'elastic_approximation':'plane-strain half-space'}
