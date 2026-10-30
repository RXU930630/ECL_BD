"""Read the supplied Persson PSD format without smoothing or extrapolation.

The columns are log10(q / m^-1) and log10(C / m^4).  C is the
two-dimensional isotropic height PSD, so h_rms^2 = 2*pi*integral(q*C dq).
Each pair of adjacent rows defines a power law (linear in log-log space).
"""

from dataclasses import dataclass
import math
from pathlib import Path

import numpy as np


def _positive_finite(value: float, name: str) -> float:
    if isinstance(value, (bool, np.bool_)):
        raise ValueError(f"{name} must be a finite positive number")
    try:
        result = float(value)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"{name} must be a finite positive number") from exc
    if not math.isfinite(result) or result <= 0:
        raise ValueError(f"{name} must be a finite positive number")
    return result


@dataclass(frozen=True)
class PowerSpectrum:
    """An immutable copy of the original logarithmic PSD samples in SI units.

    Header h0 and H are retained as metadata, not used to replace the measured
    spectrum.  ``scale_c`` scales C only; stored input columns remain unchanged.
    Evaluation is restricted to the declared q0..q1 integration interval.
    """

    log10_q: np.ndarray
    log10_c: np.ndarray
    q0: float
    qr: float
    q1: float
    roadtype: str = ""
    h0_header: float = 0.0
    H_header: float = 0.0
    scale_c: float = 1.0

    def __post_init__(self) -> None:
        x = np.array(self.log10_q, dtype=float, copy=True)
        y = np.array(self.log10_c, dtype=float, copy=True)
        if x.ndim != 1 or y.shape != x.shape or x.size < 2:
            raise ValueError("PSD requires matching one-dimensional columns with at least two rows")
        if not np.all(np.isfinite(x)) or not np.all(np.isfinite(y)):
            raise ValueError("PSD logarithmic columns must be finite")
        if not np.all(np.diff(x) > 0):
            raise ValueError("PSD wavevectors must be strictly increasing, without duplicates")
        q0 = _positive_finite(self.q0, "q0 (m^-1)")
        qr = _positive_finite(self.qr, "qr (m^-1)")
        q1 = _positive_finite(self.q1, "q1 (m^-1)")
        # Numerical spectra use the sampled curve on q0..q1. A nominal
        # roll-off below that interval remains valid header metadata.
        if qr > q1 or q0 >= q1:
            raise ValueError("PSD cutoffs must satisfy 0 < q0 < q1 and 0 < qr <= q1")
        if x[0] > math.log10(q0) + 1e-12 or x[-1] < math.log10(q1) - 1e-12:
            raise ValueError("PSD samples must cover q0..q1; extrapolation is not automatic")
        scale = _positive_finite(self.scale_c, "scale_c")
        with np.errstate(over="ignore", under="ignore", invalid="ignore"):
            q = np.power(10.0, x)
            c = np.power(10.0, y + math.log10(scale))
        if not np.all(np.isfinite(q)) or np.any(q <= 0):
            raise ValueError("PSD wavevectors are not representable as positive SI values")
        if not np.all(np.isfinite(c)) or np.any(c <= 0):
            raise ValueError("PSD amplitudes are not representable as positive SI values")
        for field in ("h0_header", "H_header"):
            try:
                value = float(getattr(self, field))
            except (TypeError, ValueError) as exc:
                raise ValueError(f"{field} must be finite") from exc
            if not math.isfinite(value):
                raise ValueError(f"{field} must be finite")
            object.__setattr__(self, field, value)
        if self.h0_header < 0:
            raise ValueError("h0_header cannot be negative")
        x.setflags(write=False)
        y.setflags(write=False)
        for name, value in (("log10_q", x), ("log10_c", y), ("q0", q0),
                            ("qr", qr), ("q1", q1), ("scale_c", scale)):
            object.__setattr__(self, name, value)

    @property
    def q(self) -> np.ndarray:
        """Input sample wavevectors in m^-1."""
        return np.power(10.0, self.log10_q)

    @property
    def c(self) -> np.ndarray:
        """Scaled input sample PSD values in m^4."""
        return np.power(10.0, self.log10_c + math.log10(self.scale_c))

    def _check_interval(self, values: np.ndarray) -> np.ndarray:
        if not np.all(np.isfinite(values)) or np.any(values <= 0):
            raise ValueError("q must contain finite positive wavevectors in m^-1")
        if np.any(values < self.q0 * (1 - 1e-12)) or np.any(values > self.q1 * (1 + 1e-12)):
            raise ValueError(f"q is outside the PSD interval [{self.q0:g}, {self.q1:g}] m^-1")
        return np.clip(values, self.q0, self.q1)

    def evaluate(self, q: float | np.ndarray) -> float | np.ndarray:
        """Interpolate C(q) in log-log space; reject out-of-range requests."""
        values = self._check_interval(np.asarray(q, dtype=float))
        log_c = np.interp(np.log10(values), self.log10_q, self.log10_c)
        result = np.power(10.0, log_c + math.log10(self.scale_c))
        return float(result) if result.ndim == 0 else result

    def moment(self, order: float, lower: float | None = None,
               upper: float | None = None) -> float:
        """Exactly integrate integral(q**order * C(q) dq), without a 2*pi factor.

        Integrals use the same piecewise power laws as ``evaluate``.  Partial
        segments and the logarithmic primitive at exponent -1 are supported.
        """
        if isinstance(order, (bool, np.bool_)) or not math.isfinite(float(order)):
            raise ValueError("moment order must be finite")
        order = float(order)
        a = self.q0 if lower is None else float(lower)
        b = self.q1 if upper is None else float(upper)
        bounds = self._check_interval(np.array([a, b]))
        a, b = map(float, bounds)
        if b < a:
            raise ValueError("moment upper bound must be at least the lower bound")
        if b == a:
            return 0.0
        la, lb = math.log(a), math.log(b)
        # Distinct adjacent floats may have the same rounded logarithm.
        # Such an interval has zero width on the integration grid.
        if la == lb:
            return 0.0
        log_q = self.log10_q * math.log(10.0)
        log_c = self.log10_c * math.log(10.0) + math.log(self.scale_c)
        knots = np.concatenate(([la], log_q[(log_q > la) & (log_q < lb)], [lb]))
        pieces = []
        for left, right in zip(knots[:-1], knots[1:]):
            # The midpoint selects the segment even when a bound equals a knot.
            i = int(np.searchsorted(log_q, (left + right) / 2, side="right") - 1)
            i = min(max(i, 0), len(log_q) - 2)
            slope = (log_c[i + 1] - log_c[i]) / (log_q[i + 1] - log_q[i])
            log_at_left = log_c[i] + slope * (left - log_q[i])
            width = right - left
            z = (order + 1 + slope) * width
            # log(expm1(z)/z) remains stable for nearly logarithmic primitives.
            if abs(z) < 1e-6:
                log_ratio = math.log1p(z / 2 + z*z / 6 + z*z*z / 24)
            elif z > 50:
                log_ratio = z + math.log1p(-math.exp(-z)) - math.log(z)
            elif z < -50:
                log_ratio = math.log1p(-math.exp(z)) - math.log(-z)
            else:
                log_ratio = math.log(math.expm1(z) / z)
            log_integral = log_at_left + (order + 1) * left + math.log(width) + log_ratio
            try:
                pieces.append(math.exp(log_integral))
            except OverflowError as exc:
                raise ValueError("PSD moment exceeds floating-point range") from exc
        result = math.fsum(pieces)
        if not math.isfinite(result):
            raise ValueError("PSD moment exceeds floating-point range")
        return result

    @property
    def hrms(self) -> float:
        """RMS height in m, computed from C rather than the header h0."""
        return math.sqrt(2 * math.pi * self.moment(1))

    @property
    def rms_slope(self) -> float:
        """Dimensionless RMS surface gradient."""
        return math.sqrt(2 * math.pi * self.moment(3))


PSD = PowerSpectrum


def _number(token: str, line: int) -> float:
    try:
        result = float(token.replace("D", "E").replace("d", "e"))
    except ValueError as exc:
        raise ValueError(f"Invalid PSD number on line {line}: {token!r}") from exc
    if not math.isfinite(result):
        raise ValueError(f"Non-finite PSD number on line {line}")
    return result


def read_psd(path: str | Path, scale_c: float = 1.0) -> PowerSpectrum:
    """Read ROADTYPE, q0, qr, q1, h0, H followed by log10(q), log10(C).

    This reader has a fixed SI/isotropic-2D convention.  It does not infer
    units, smooth data, fit H, or extend the supplied sampling range.
    """
    lines = Path(path).read_text(encoding="utf-8-sig").splitlines()
    if len(lines) < 8 or not lines[0].strip():
        raise ValueError("PSD requires ROADTYPE, five numeric header rows and at least two data rows")
    header = []
    for i in range(1, 6):
        tokens = lines[i].split()
        if not tokens:
            raise ValueError(f"Empty PSD header on line {i + 1}")
        header.append(_number(tokens[0], i + 1))
    rows = []
    for lineno, line in enumerate(lines[6:], start=7):
        content = line.split("#", 1)[0].split("!", 1)[0].strip()
        if not content:
            continue
        tokens = content.split()
        if len(tokens) != 2:
            raise ValueError(f"PSD line {lineno} must contain exactly two logarithmic columns")
        rows.append([_number(token, lineno) for token in tokens])
    if len(rows) < 2:
        raise ValueError("PSD requires at least two numerical data rows")
    data = np.asarray(rows)
    q0, qr, q1, h0, hurst = header
    return PowerSpectrum(data[:, 0], data[:, 1], q0, qr, q1,
                         lines[0].strip(), h0, hurst, scale_c)
