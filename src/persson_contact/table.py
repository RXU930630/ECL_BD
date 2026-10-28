"""Build and query a bounded, monotone dry-contact database."""
from dataclasses import dataclass
from typing import Any

import numpy as np
from scipy.interpolate import PchipInterpolator


@dataclass(frozen=True)
class ContactResponse:
    mean_gap_m: Any
    pressure_pa: Any
    contact_fraction: Any
    gap_distribution: Any = None
    interpolation_mean_gap_error_m: float | None = None


@dataclass
class ContactTable:
    model: Any
    pressures_pa: np.ndarray
    mean_gaps_m: np.ndarray
    contact_fractions: np.ndarray

    def __post_init__(self):
        p, h, a = (np.asarray(v, dtype=float).copy() for v in
                   (self.pressures_pa, self.mean_gaps_m, self.contact_fractions))
        if p.ndim != 1 or p.size < 3 or p.shape != h.shape or p.shape != a.shape:
            raise ValueError('A contact table requires at least three matching one-dimensional arrays.')
        if not all(np.all(np.isfinite(v)) for v in (p, h, a)):
            raise ValueError('Contact-table values must be finite.')
        if np.any(p <= 0) or np.any(np.diff(p) <= 0):
            raise ValueError('Pressures must be positive and strictly increasing.')
        if np.any(h <= 0) or np.any(np.diff(h) >= 0):
            raise ValueError('Mean gaps must be positive and strictly decreasing; narrow the pressure range or refine the model grid.')
        if np.any((a < 0) | (a > 1)) or np.any(np.diff(a) < -1e-12):
            raise ValueError('Contact fractions must lie in [0,1] and increase with pressure.')
        self.pressures_pa, self.mean_gaps_m, self.contact_fractions = p, h, a
        # Interpolate log pressure to keep its positivity, without extrapolation.
        self._log_pressure = PchipInterpolator(h[::-1], np.log(p[::-1]), extrapolate=False)
        for array in (p, h, a):
            array.flags.writeable = False

    @property
    def valid_gap_range_m(self):
        return float(self.mean_gaps_m[-1]), float(self.mean_gaps_m[0])

    def evaluate(self, mean_gap_m, gaps_m=None):
        h = np.asarray(mean_gap_m, dtype=float)
        low, high = self.valid_gap_range_m
        if not np.all(np.isfinite(h)) or np.any(h < low) or np.any(h > high):
            raise ValueError(f'Mean gap is outside table range [{low:.8g}, {high:.8g}] m; extend the table explicitly.')
        p = np.exp(self._log_pressure(h))
        # Area is cheap to evaluate. Use the same pressure for area and PDF,
        # rather than introducing an independent interpolation of a(h).
        a = self.model.area(p)
        distribution = None
        gap_error = None
        if gaps_m is not None:
            if h.ndim:
                raise ValueError('Query a single mean gap when requesting a microgap distribution.')
            distribution = self.model.gap_distribution(float(p), gaps_m)
            gap_error = float(distribution.mean_gap_m - h)
        scalar = h.ndim == 0
        return ContactResponse(float(h) if scalar else h,
                               float(p) if scalar else p,
                               float(a) if scalar else a,
                               distribution, gap_error)


def build_contact_table(model, pressures_pa):
    p = np.asarray(pressures_pa, dtype=float)
    if p.ndim != 1 or p.size < 3 or not np.all(np.isfinite(p)) or np.any(p <= 0) or np.any(np.diff(p) <= 0):
        raise ValueError('Supply at least three positive, strictly increasing pressures in Pa.')
    return ContactTable(model, p, model.mean_gap(p), model.area(p))


def evaluate_contact(table, mean_gap_m, gaps_m=None):
    return table.evaluate(mean_gap_m, gaps_m=gaps_m)
