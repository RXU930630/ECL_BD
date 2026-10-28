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

    def __post_init__(self) -> None:
        for name in ("gamma", "pdf_afact"):
            value = getattr(self, name)
            if isinstance(value, (bool, np.bool_)) or not math.isfinite(float(value)):
                raise ValueError(f"{name} must be a finite number")
            object.__setattr__(self, name, float(value))
        if not 0 <= self.gamma <= 1:
            raise ValueError("gamma must lie in [0, 1]")
        if self.pdf_afact < 0:
            raise ValueError("pdf_afact must be nonnegative")
        object.__setattr__(self, "n_q", _integer(self.n_q, "n_q", 32))
        object.__setattr__(self, "n_gap", _integer(self.n_gap, "n_gap", 32))


@dataclass(frozen=True)
class Problem:
    """Dry-contact inputs; effective modulus already includes both compliances."""

    psd: PowerSpectrum
    effective_modulus_pa: float
    pressures_pa: np.ndarray
    psd_tail_is_extrapolated: bool = False
    description: str = ""

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
    _keys(raw, {"psd", "material", "pressures_pa", "options", "description"}, "problem")
    for field in ("psd", "material", "pressures_pa", "options"):
        if field not in raw:
            raise ValueError(f"Missing required problem field: {field}")
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
    _keys(opts, {"gamma", "n_q", "pdf_afact", "n_gap"}, "options")
    for field in ("gamma", "n_q", "pdf_afact"):
        if field not in opts:
            raise ValueError(f"Missing required options field: {field}")
    options = Options(**opts)
    problem = Problem(psd, effective_modulus, raw["pressures_pa"],
                      psd_config.get("tail_is_extrapolated", False), raw.get("description", ""))
    return problem, options


read_problem = load_problem
