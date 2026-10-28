"""Elastic, nonadhesive Persson contact and an experimental microgap PDF.

The PSD convention is a two-dimensional isotropic C(q) [m**4], q [1/m].
``gamma`` is the elastic-energy correction, *not* a friction coefficient.
The PDF is a reflected-Gaussian approximation inspired by the supplied Fortran
code. Its first moment need not equal the energy-derived mean gap; that error
is reported, never removed by rescaling gaps or renormalizing sampled data.
"""

from dataclasses import asdict, dataclass
from functools import lru_cache
from typing import Protocol

import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.integrate import cumulative_simpson, trapezoid
from scipy.interpolate import PchipInterpolator
from scipy.special import erf, exp1, ndtr


class Spectrum(Protocol):
    q0: float
    q1: float

    def evaluate(self, q): ...

    def moment(self, order, lower=None, upper=None): ...


class ModelResolutionError(ValueError):
    """The selected discretization cannot produce a physical PDF mixture."""


@lru_cache(maxsize=64)
def _gauss_rule(count):
    return leggauss(count)


@lru_cache(maxsize=1)
def _correction_integral():
    """J(z)=integral_z^inf erf(x)**2 exp(-x*x)/x dx, 0<=z<=8.

    Positive reverse integration avoids subtracting nearly equal integrals in
    the tail. At x>=8, replacing erf(x)**2 by 1 has negligible double-precision
    error. This cache is independent of materials and PSDs.
    """
    x = np.linspace(0.0, 8.0, 16385)
    f = np.zeros_like(x)
    f[1:] = erf(x[1:]) ** 2 * np.exp(-x[1:] ** 2) / x[1:]
    reverse = cumulative_simpson(f[::-1], x=-x[::-1], initial=0.0)
    tail = reverse[::-1] + 0.5 * exp1(64.0)
    return PchipInterpolator(x, tail, extrapolate=False)


def separation_kernel(z, gamma):
    r"""F(z)=integral_z^inf [gamma+3(1-gamma)erf²(x)]e^-x² dx/x.

    gamma=1 uses the exact exponential-integral expression. The correction
    term for other gamma values is cached and interpolated. F(0) is infinite
    when gamma>0; F(inf)=0.
    """
    gamma = float(gamma)
    if not np.isfinite(gamma) or not 0 <= gamma <= 1:
        raise ValueError("gamma must be finite and lie in [0, 1]")
    z = np.asarray(z, dtype=float)
    if np.any(np.isnan(z)) or np.any(z < 0):
        raise ValueError("kernel argument must be nonnegative")
    with np.errstate(over="ignore"):
        exact = 0.5 * exp1(z * z)
    if gamma == 1:
        result = exact
    else:
        j = np.array(exact, copy=True)
        inside = z <= 8.0
        j[inside] = _correction_integral()(z[inside])
        result = 3.0 * (1.0 - gamma) * j
        if gamma > 0:
            result = result + gamma * exact
    return float(result) if result.ndim == 0 else result


@dataclass(frozen=True)
class GapDiagnostics:
    expected_continuous_mass: float
    resolved_continuous_mass: float
    sampled_continuous_mass: float
    expected_mean_gap_m: float
    analytic_first_moment_m: float
    sampled_first_moment_m: float
    first_moment_relative_error: float
    sampled_first_moment_relative_error: float
    magnification_intervals: int
    unresolved_continuous_mass: float
    integrated_continuous_mass_on_grid: float
    unresolved_mean_gap_contribution_m: float
    approximation: str = "experimental_reflected_gaussian"
    notes: tuple[str, ...] = ()

    def as_dict(self):
        return asdict(self)


@dataclass(frozen=True)
class GapDistribution:
    pressure_pa: float
    mean_gap_m: float
    contact_fraction: float
    gaps_m: np.ndarray
    continuous_density_per_m: np.ndarray
    positive_density_per_m: np.ndarray
    diagnostics: GapDiagnostics
    continuous_cdf: np.ndarray
    positive_cdf: np.ndarray

    @property
    def relative_contact_area(self):
        """A/A0; contact_fraction remains available as a legacy API alias."""
        return self.contact_fraction


@dataclass(frozen=True)
class CriticalSeparation:
    pressure_pa: float
    target_relative_contact_area: float
    status: str
    magnification: float | None = None
    wavevector_per_m: float | None = None
    relative_contact_area: float | None = None
    area_deviation: float | None = None
    mean_gap_at_scale_m: float | None = None
    u1_m: float | None = None
    critical_separation_m: float | None = None
    u1_refinement_relative_difference: float | None = None


class PerssonModel:
    """Homogeneous elastic half-space against a rigid rough surface.

    ``effective_modulus_pa`` is E/(1-nu**2), or the combined effective modulus
    for two elastic bodies. Standard Persson area is used, with no empirical
    "exact area" fitting. gamma affects the energy/mean-gap relation only.
    """

    def __init__(self, psd: Spectrum, effective_modulus_pa, gamma, n_q=256,
                 pdf_afact=1.0):
        self.psd = psd
        self.effective_modulus_pa = float(effective_modulus_pa)
        self.gamma = float(gamma)
        self.pdf_afact = float(pdf_afact)
        if not np.isfinite(self.effective_modulus_pa) or self.effective_modulus_pa <= 0:
            raise ValueError("effective_modulus_pa must be positive and finite")
        if not np.isfinite(self.gamma) or not 0 <= self.gamma <= 1:
            raise ValueError("gamma must lie in [0, 1]")
        if isinstance(n_q, bool) or int(n_q) != n_q or n_q < 32:
            raise ValueError("n_q must be an integer >=32")
        if not np.isfinite(self.pdf_afact) or self.pdf_afact < 0:
            raise ValueError("pdf_afact must be finite and nonnegative")
        self.n_q = int(n_q)
        self.q0 = float(psd.q0)
        self.q1 = float(psd.q1)
        if not 0 < self.q0 < self.q1 or not np.isfinite(self.q1):
            raise ValueError("PSD must have finite 0<q0<q1")
        self._total_i3 = float(psd.moment(3))
        if not np.isfinite(self._total_i3) or self._total_i3 <= 0:
            raise ValueError("PSD slope moment must be positive and finite")
        self._legendre_nodes, self._legendre_weights = leggauss(self.n_q)
        self._quadrature_cache = {}

    @staticmethod
    def _pressures(pressure_pa, *, allow_zero=False):
        p = np.asarray(pressure_pa, dtype=float)
        if (np.any(~np.isfinite(p)) or np.any(p < 0)
                or (not allow_zero and np.any(p == 0))):
            boundary = "nonnegative" if allow_zero else "positive"
            raise ValueError(f"pressure_pa must be finite and {boundary}")
        return p

    def area(self, pressure_pa, q=None):
        """Relative contact area A/A0; q optionally sets observation cutoff.

        At q=q0 there is no resolved roughness, so area is one for p>0.
        At p=0 this implementation defines area as zero.
        """
        p = self._pressures(pressure_pa, allow_zero=True)
        if q is None:
            integral = np.asarray(self._total_i3)
        else:
            q = np.asarray(q, dtype=float)
            if np.any(~np.isfinite(q)) or np.any((q < self.q0) | (q > self.q1)):
                raise ValueError("q must lie within the PSD's [q0, q1]")
            integral = np.asarray([
                self.psd.moment(3, lower=self.q0, upper=float(qi))
                if qi > self.q0 else 0.0 for qi in q.flat
            ]).reshape(q.shape)
        p, integral = np.broadcast_arrays(p, integral)
        ratio = np.zeros_like(p)
        np.divide(p, self.effective_modulus_pa * np.sqrt(np.pi * integral),
                  out=ratio, where=integral > 0)
        result = np.where(integral > 0, erf(ratio), np.where(p > 0, 1.0, 0.0))
        return float(result) if result.ndim == 0 else result

    def _quadrature(self, lower):
        """Integrate in t=sqrt(log(q/lower)); cancels lower-end singularity."""
        lower = float(lower)
        cached = self._quadrature_cache.get(lower)
        if cached is not None:
            return cached
        stop = np.sqrt(np.log(self.q1 / lower))
        # Split at every input PSD knot: global high-order Gauss integration
        # across log-log slope changes otherwise jitters as ``lower`` moves,
        # contaminating the magnification-volume differences used for u1.
        knots = getattr(self.psd, "log10_q", None)
        if knots is None:
            bounds = np.asarray([0.0, stop])
        else:
            knots = 10.0 ** np.asarray(knots, dtype=float)
            inside = knots[(knots > lower) & (knots < self.q1)]
            bounds = np.r_[0.0, np.sqrt(np.log(inside / lower)), stop]
        nodes, weights = [], []
        for left, right in zip(bounds[:-1], bounds[1:]):
            count = max(8, int(np.ceil(self.n_q * (right - left) / stop)))
            gx, gw = _gauss_rule(count)
            nodes.append(left + 0.5 * (right - left) * (gx + 1.0))
            weights.append(0.5 * (right - left) * gw)
        t = np.concatenate(nodes)
        quadrature_weights = np.concatenate(weights)
        q = lower * np.exp(t * t)
        # Floating-point exponentiation must not move the last node out of band.
        q = np.minimum(q, self.q1)
        i3 = np.asarray([self.psd.moment(3, lower=lower, upper=float(qi)) for qi in q])
        if np.any(i3 <= 0) or np.any(~np.isfinite(i3)):
            raise ModelResolutionError("invalid local PSD moment; check spectral range/resolution")
        coefficient = (quadrature_weights * 2.0 * t
                       * q ** 3 * self.psd.evaluate(q) / np.sqrt(i3))
        scale = np.sqrt(np.pi * i3)
        self._quadrature_cache[lower] = (coefficient, scale)
        return coefficient, scale

    def _mean_gap_above(self, pressure_pa, lower):
        if lower >= self.q1:
            return 0.0
        coefficient, scale = self._quadrature(lower)
        argument = (pressure_pa / self.effective_modulus_pa) / scale
        return float(coefficient @ separation_kernel(argument, self.gamma))

    def mean_gap(self, pressure_pa):
        """Energy-derived mean gap [m] for strictly positive nominal pressure."""
        p = self._pressures(pressure_pa)
        result = np.asarray([self._mean_gap_above(float(pi), self.q0) for pi in p.flat])
        result = result.reshape(p.shape)
        return float(result) if result.ndim == 0 else result

    def critical_separation(self, pressure_pa, target=0.42, area_tolerance=0.005):
        """Approximate u_c=u1 near A(zeta)/A0=target, with alpha fixed to one.

        Only a bracketed threshold is accepted. Bisection stops at the stated
        absolute area tolerance, not an exact root. u1 is evaluated from local
        finite-volume derivatives of A*ubar with respect to A, with a half-step
        comparison to expose derivative sensitivity. This retains the existing
        experimental separation model, not a validated leakage-rate solver.
        """
        p_array = self._pressures(pressure_pa)
        if p_array.ndim != 0:
            raise ValueError('critical_separation requires one scalar pressure')
        p = float(p_array)
        target, tolerance = float(target), float(area_tolerance)
        if not (np.isfinite(target) and 0 < target < 1
                and np.isfinite(tolerance) and 0 < tolerance < min(target, 1-target)):
            raise ValueError('invalid critical relative contact area or absolute tolerance')
        if self.area(p) > target:
            return CriticalSeparation(p, target, 'threshold_not_reached')
        low, high = np.log(self.q0), np.log(self.q1)
        for _ in range(80):
            center = (low + high) / 2
            q = float(np.clip(np.exp(center), self.q0, self.q1))
            a = self.area(p, q)
            if abs(a - target) <= tolerance:
                break
            if a > target:
                low = center
            else:
                high = center
        else:
            return CriticalSeparation(p, target, 'area_resolution_insufficient')
        mean = self._mean_gap_above(p/a, q)
        step = min(np.log(self.q1/self.q0)/(4*self.n_q),
                   (center-np.log(self.q0))/4, (np.log(self.q1)-center)/4)

        def local_u1(width):
            ql, qr = np.exp(center-width), np.exp(center+width)
            al, ar = self.area(p, ql), self.area(p, qr)
            if al <= ar:
                return None
            ul = self._mean_gap_above(p/al, ql)
            ur = self._mean_gap_above(p/ar, qr)
            value = (al*ul-ar*ur)/(al-ar)
            return float(value) if np.isfinite(value) and value > 0 else None

        coarse = local_u1(step) if step > 0 else None
        fine = local_u1(step/2) if step > 0 else None
        difference = None
        for _ in range(8):
            if fine is None or coarse is None:
                break
            difference = abs(fine-coarse)/fine
            if difference <= 1e-3:
                break
            step /= 2
            coarse, fine = fine, local_u1(step/2)
        difference = abs(fine-coarse)/fine if fine is not None and coarse is not None else None
        resolved = fine is not None and difference is not None and difference <= 1e-3
        status = 'found' if resolved else ('u1_unresolved' if fine is None else 'u1_derivative_unconverged')
        return CriticalSeparation(
            p, target, status, q/self.q0, q, a, a-target, mean, fine,
            fine if resolved else None, difference)

    def gap_distribution(self, pressure_pa, gaps_m):
        """Approximate microgap distribution, including a separate zero atom.

        ``continuous_density_per_m`` describes the numerically resolved part
        of noncontact area. Its full-axis mass plus explicitly reported
        ``unresolved_continuous_mass`` equals 1-A/A0. The positive density is
        divided by the full noncontact probability, and therefore its mass is
        1-unresolved_mass/(1-A/A0). Neither curve is renormalized to hide loss.
        Exact Gaussian CDF differences provide interval probabilities; trapz
        of point densities can be unreliable for sharp near-zero components.

        A finite-volume form of u1=ubar+A*ubar'/A' replaces the Fortran's finite
        differences. Each magnification interval has mass Delta A, exactly
        partitioning noncontact area. Reflected-Gaussian broadening follows the
        source's width model; its first-moment discrepancy is exposed below.
        """
        p_array = self._pressures(pressure_pa)
        if p_array.ndim != 0:
            raise ValueError("gap_distribution requires one scalar pressure")
        p = float(p_array)
        gaps = np.asarray(gaps_m, dtype=float)
        if (gaps.ndim != 1 or gaps.size < 2 or np.any(~np.isfinite(gaps))
                or np.any(gaps < 0) or np.any(np.diff(gaps) <= 0)):
            raise ValueError("gaps_m must be a strictly increasing nonnegative 1D grid")

        # Squared log spacing resolves the initial A~1 boundary layer.
        t = np.linspace(0.0, np.sqrt(np.log(self.q1 / self.q0)), self.n_q + 1)
        q = self.q0 * np.exp(t * t)
        q[0], q[-1] = self.q0, self.q1
        a = self.area(p, q)
        contact = float(a[-1])
        expected_mass = 1.0 - contact
        if expected_mass == 0:
            raise ValueError("area rounds to full contact; conditional positive-gap PDF is undefined")
        h = self.mean_gap(p)
        ubar = np.asarray([self._mean_gap_above(p / ai, float(qi))
                           for ai, qi in zip(a, q)])
        volume = a * ubar
        weights = a[:-1] - a[1:]
        selected = weights > 0.0
        if not np.any(selected):
            raise ModelResolutionError("no resolved noncontact-area intervals; increase n_q")
        centers = (volume[:-1][selected] - volume[1:][selected]) / weights[selected]
        weights = weights[selected]
        if np.any(~np.isfinite(centers)):
            raise ModelResolutionError(
                "nonfinite u1 from magnification-volume differences; inspect inputs"
            )
        midpoint_q = np.sqrt(q[:-1][selected] * q[1:][selected])
        # Near q1, conditional mean gaps can underflow before Delta A does.
        # Report their noncontact mass separately; never transfer it to contact.
        # Negative-volume intervals are only tolerated at roundoff-level area.
        negative = centers < 0
        if np.sum(weights[negative]) > 1e-10 * expected_mass:
            raise ModelResolutionError(
                "nonphysical negative u1 carries resolved area; increase n_q or inspect inputs"
            )
        unresolved = centers <= 0
        unresolved_mass = float(np.sum(weights[unresolved]))
        centers, weights, midpoint_q = (values[~unresolved]
                                        for values in (centers, weights, midpoint_q))
        variance_roughness = 2.0 * np.pi * np.asarray([
            self.psd.moment(1, lower=float(qi), upper=self.q1) for qi in midpoint_q
        ])
        if self.pdf_afact == 0:
            widths = np.sqrt(variance_roughness)
        else:
            # Stable even when center**2 would underflow.
            widths = centers / np.sqrt(self.pdf_afact + centers ** 2 / variance_roughness)
        if np.any(widths <= 0) or np.any(~np.isfinite(widths)):
            raise ModelResolutionError("invalid broadening widths; inspect PSD and pdf_afact")
        density = np.zeros_like(gaps)
        cdf = np.zeros_like(gaps)
        # Accumulate one interval at a time to bound memory for large user grids.
        for weight, center, width in zip(weights, centers, widths):
            # Overflowed squared standardized distances are valid zero tails.
            with np.errstate(over="ignore", under="ignore"):
                density += weight * (
                    np.exp(-0.5 * ((gaps - center) / width) ** 2)
                    + np.exp(-0.5 * ((gaps + center) / width) ** 2)
                ) / (np.sqrt(2.0 * np.pi) * width)
                cdf += weight * (ndtr((gaps - center) / width)
                                 - ndtr((-gaps - center) / width))
        if np.any(~np.isfinite(density)):
            raise ModelResolutionError("density is unresolved near zero gap; use a positive grid/CDF bins")
        conditional = density / expected_mass
        ratio = centers / widths
        reflected_means = (centers * erf(ratio / np.sqrt(2.0))
                           + widths * np.sqrt(2.0 / np.pi) * np.exp(-0.5 * ratio ** 2))
        analytic_first = float(weights @ reflected_means)
        sampled_mass = float(trapezoid(density, gaps))
        sampled_first = float(trapezoid(density * gaps, gaps))
        notes = [
            "Experimental reflected-Gaussian approximation; not validated against a complete Fortran run.",
            "Contact probability is a separate Dirac atom at u=0, not a finite-density sample.",
            "Finite-volume magnification discretization differs from the original finite differences.",
            "Broadening can increase the first moment; no density/gap rescaling is applied.",
            "Near-zero density components may be much narrower than the grid: "
            "use analytic CDF differences for bin probabilities, not trapz of density samples.",
        ]
        if unresolved_mass:
            notes.append(
                f"Unresolved noncontact probability {unresolved_mass:.8g} is reported separately "
                "(mean-gap underflow or roundoff); it is neither contact nor renormalized away."
            )
        if abs(sampled_mass - expected_mass) > 0.01 * expected_mass:
            notes.append("The supplied gap grid misses mass or under-resolves the density; extend/refine it.")
        if h <= 0:
            raise ModelResolutionError("mean gap underflow prevents PDF moment diagnostics")
        diagnostics = GapDiagnostics(
            expected_continuous_mass=expected_mass,
            resolved_continuous_mass=float(np.sum(weights)),
            sampled_continuous_mass=sampled_mass,
            expected_mean_gap_m=h,
            analytic_first_moment_m=analytic_first,
            sampled_first_moment_m=sampled_first,
            first_moment_relative_error=(analytic_first - h) / h,
            sampled_first_moment_relative_error=(sampled_first - h) / h,
            magnification_intervals=len(weights),
            unresolved_continuous_mass=unresolved_mass,
            integrated_continuous_mass_on_grid=float(cdf[-1] - cdf[0]),
            unresolved_mean_gap_contribution_m=float(h - weights @ centers),
            notes=tuple(notes),
        )
        return GapDistribution(p, h, contact, gaps.copy(), density, conditional,
                               diagnostics, cdf, cdf / expected_mass)
