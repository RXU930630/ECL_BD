"""Validated SI input for the elastic, nonadhesive Persson workflow."""

from dataclasses import dataclass
import json
import math
from pathlib import Path

import numpy as np

from .psd import PowerSpectrum, _positive_finite, read_psd
from .material import parse_material


def _integer(value: int, name: str, minimum: int) -> int:
    if isinstance(value, (bool, np.bool_)) or not isinstance(value, (int, np.integer)) or value < minimum:
        raise ValueError(f"{name} must be an integer >= {minimum}")
    return int(value)


@dataclass(frozen=True)
class Options:
    """Explicit numerical and model choices; no inferred research parameters."""

    gamma: float
    n_q: int
    pdf_afact: float
    n_gap: int = 256
    critical_relative_contact_area: float = 0.42
    critical_area_tolerance: float = 0.005
    pdf_log_step: float = 0.01
    pdf_max_relative_contact_area: float = 0.999
    pdf_min_relative_contact_area: float = 0.0001
    pdf_min_u1_m: float = 1e-12
    pdf_scale_method: str = 'fortran'

    def __post_init__(self) -> None:
        if self.pdf_scale_method not in ('fortran', 'conservative'):
            raise ValueError('pdf_scale_method must be fortran or conservative')
        for name in ("gamma", "pdf_afact", "pdf_log_step", "pdf_max_relative_contact_area",
                     "pdf_min_relative_contact_area", "pdf_min_u1_m"):
            value = getattr(self, name)
            if isinstance(value, (bool, np.bool_)) or not math.isfinite(float(value)):
                raise ValueError(f"{name} must be a finite number")
            object.__setattr__(self, name, float(value))
        if not 0 <= self.gamma <= 1:
            raise ValueError("gamma must lie in [0, 1]")
        if self.pdf_afact < 0:
            raise ValueError("pdf_afact must be nonnegative")
        if self.pdf_log_step <= 0 or self.pdf_min_u1_m <= 0:
            raise ValueError("pdf_log_step and pdf_min_u1_m must be positive")
        if not 0 < self.pdf_min_relative_contact_area < self.pdf_max_relative_contact_area < 1:
            raise ValueError("PDF area cutoffs must satisfy 0 < min < max < 1")
        object.__setattr__(self, "n_q", _integer(self.n_q, "n_q", 32))
        object.__setattr__(self, "n_gap", _integer(self.n_gap, "n_gap", 32))
        for name in ('critical_relative_contact_area', 'critical_area_tolerance'):
            value = getattr(self, name)
            if isinstance(value, bool) or not math.isfinite(float(value)):
                raise ValueError(f'{name} must be finite')
            object.__setattr__(self, name, float(value))
        if not 0 < self.critical_relative_contact_area < 1:
            raise ValueError('critical_relative_contact_area must lie in (0,1)')
        if not 0 < self.critical_area_tolerance < min(self.critical_relative_contact_area,
                                                     1-self.critical_relative_contact_area):
            raise ValueError('critical_area_tolerance must be positive and smaller than either area fraction')


@dataclass(frozen=True)
class Problem:
    """Dry-contact inputs; effective modulus already includes both compliances."""

    psd: PowerSpectrum
    effective_modulus_pa: float
    pressures_pa: np.ndarray
    psd_tail_is_extrapolated: bool = False
    description: str = ""
    reference_pressure_pa: float | None = None

    def __post_init__(self) -> None:
        modulus = _positive_finite(self.effective_modulus_pa, "effective_modulus_pa")
        pressures = np.array(self.pressures_pa, dtype=float, copy=True)
        if (pressures.ndim != 1 or pressures.size == 0
                or not np.all(np.isfinite(pressures)) or np.any(pressures <= 0)
                or not np.all(np.diff(pressures) > 0)):
            raise ValueError("pressures_pa must be a nonempty, strictly increasing list of positive finite values")
        if not isinstance(self.psd_tail_is_extrapolated, bool):
            raise ValueError("tail_is_extrapolated must be true or false")
        if not isinstance(self.description, str):
            raise ValueError("description must be a string")
        pressures.setflags(write=False)
        object.__setattr__(self, "effective_modulus_pa", modulus)
        object.__setattr__(self, "pressures_pa", pressures)


def _keys(data: dict, allowed: set[str], context: str) -> None:
    if not isinstance(data, dict):
        raise ValueError(f"{context} must be a JSON object")
    unknown = set(data) - allowed
    if unknown:
        raise ValueError(f"Unknown {context} field(s): {', '.join(sorted(unknown))}")


def pressure_grid(spec):
    """Exactly N logarithmic nodes, with reference and both bounds retained.

    A nearby interior log node is replaced by the reference pressure. Thus the
    grid is approximately logarithmic, with no extra node or interpolation of
    the user's requested reference state. The automatic lower bound is 1e-4 p_ref.
    """
    _keys(spec, {'reference_pa', 'upper_factor', 'number_of_points'}, 'contact_pressure')
    if 'reference_pa' not in spec:
        raise ValueError('contact_pressure requires reference_pa')
    reference = _positive_finite(spec['reference_pa'], 'reference_pa')
    factor = _positive_finite(spec.get('upper_factor', 2.0), 'upper_factor')
    if factor < 1:
        raise ValueError('upper_factor must be >=1 to cover the reference pressure')
    count = _integer(spec.get('number_of_points', 110), 'number_of_points', 3)
    lower, upper = reference * 1e-4, reference * factor
    if not (math.isfinite(upper) and lower > 0):
        raise ValueError('pressure bounds exceed floating-point range')
    values = np.geomspace(lower, upper, count)
    if factor == 1:
        values[-1] = reference
    else:
        index = 1 + int(np.argmin(abs(np.log(values[1:-1]) - math.log(reference))))
        values[index] = reference
    if not np.all(np.diff(values) > 0):
        raise ValueError('pressure grid is not representable at this resolution')
    return values


def load_problem(path: str | Path) -> tuple[Problem, Options]:
    """Load a JSON configuration; resolve PSD paths relative to that JSON file.

    ``material`` contains either effective_modulus_pa OR young_modulus_pa and
    poisson_ratio.  The latter assumes a deformable half-space against a rigid
    counterface and computes E*=E/(1-nu^2) once.  Already-effective moduli are
    never corrected a second time.  Fluid/geometry data are not required.
    """
    path = Path(path)
    try:
        raw = json.loads(path.read_text(encoding="utf-8-sig"))
    except json.JSONDecodeError as exc:
        raise ValueError(f"Invalid problem JSON: {exc}") from exc
    return problem_from_mapping(raw, path.parent)


def problem_from_mapping(raw, base_directory=Path('.')) -> tuple[Problem, Options]:
    """Shared validation for Python parameters and the original JSON interface."""
    _keys(raw, {"psd", "material", "pressures_pa", "contact_pressure", "options", "description"}, "problem")
    for field in ("psd", "material", "options"):
        if field not in raw:
            raise ValueError(f"Missing required problem field: {field}")
    if ('pressures_pa' in raw) == ('contact_pressure' in raw):
        raise ValueError('Provide contact_pressure OR pressures_pa, not both')
    pressures = pressure_grid(raw['contact_pressure']) if 'contact_pressure' in raw else raw['pressures_pa']
    material = parse_material(raw['material'], base_directory)
    effective_modulus = material.elastic_modulus()
    psd_config = raw["psd"]
    _keys(psd_config, {"path", "scale_c", "tail_is_extrapolated", "q_unit", "c_unit"}, "psd")
    if not isinstance(psd_config.get("path"), str) or not psd_config["path"].strip():
        raise ValueError("psd.path must be a nonempty filename")
    if psd_config.get("q_unit", "m^-1") != "m^-1" or psd_config.get("c_unit", "m^4") != "m^4":
        raise ValueError("This PSD format requires q_unit='m^-1' and c_unit='m^4'; convert units explicitly")
    psd_path = Path(psd_config["path"])
    if not psd_path.is_absolute():
        psd_path = Path(base_directory) / psd_path
    psd = read_psd(psd_path, psd_config.get("scale_c", 1.0))
    opts = raw["options"]
    _keys(opts, {"gamma", "n_q", "pdf_afact", "n_gap", "critical_relative_contact_area", "critical_area_tolerance",
                  "pdf_log_step", "pdf_max_relative_contact_area",
                  "pdf_min_relative_contact_area", "pdf_min_u1_m", "pdf_scale_method"}, "options")
    for field in ("gamma", "n_q", "pdf_afact"):
        if field not in opts:
            raise ValueError(f"Missing required options field: {field}")
    options = Options(**opts)
    problem = Problem(psd, effective_modulus, pressures,
                      psd_config.get("tail_is_extrapolated", False), raw.get("description", ""),
                      float(raw['contact_pressure']['reference_pa']) if 'contact_pressure' in raw else None)
    return problem, options


read_problem = load_problem
