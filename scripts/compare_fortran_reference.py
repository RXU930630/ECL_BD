#!/usr/bin/env python3
"""Compare completed original-Fortran dry runs with the Python implementation.

This reads reference outputs; it never starts or modifies a Fortran run.
Material/numerical values identify archived parameters, not a calibrated
experiment. The Python gap PDF remains an experimental approximation.
"""

import argparse
from dataclasses import asdict
import hashlib
import json
import math
from pathlib import Path
import re
import sys
import tempfile
import unittest

import numpy as np
from scipy.integrate import trapezoid

from persson_contact.core import PerssonModel
from persson_contact.psd import read_psd


PSD_NAME = "IN_ROADTYPE.q0.qr.q1.h0.H.1logWavevector.2logPowerSpectrum"
AREA_NAME = "out_No.Adhesion.1pressure.2A.3Aexact"
MEAN_NAME = "RealUnits.1pressure.Pa.2average.interfacial.separation.meter.3A.4Aexact"
SUMMARY_NAME = "out.1pressure.2baruFromP.3baruFromEnergy.4AreaFromP.5AreaFromEnergy"
PDF_NAME = "out.1u.2ProbabilityInterfacialSeparation"
MATERIAL_NAME = "IN.sigma.gammabar.xnu.ElasticModul.yieldstress"
MISSING_E = re.compile(r"^([+-]?(?:\d+(?:\.\d*)?|\.\d+))([+-]\d{3,})$")


def fortran_float(token):
    """Parse D exponents and Fortran Ew.d three-digit exponents missing E."""
    token = token.strip().replace("D", "E").replace("d", "e")
    if "e" not in token.lower():
        match = MISSING_E.fullmatch(token)
        if match:
            token = match[1] + "e" + match[2]
    number = float(token)
    if not math.isfinite(number):
        raise ValueError(f"nonfinite Fortran field: {token!r}")
    return number


def read_numeric(path, columns, *, pressure_header=False, minimum_rows=1):
    rows = []
    header = None
    for lineno, line in enumerate(Path(path).read_text(encoding="utf-8-sig").splitlines(), 1):
        line = line.strip()
        if not line or line.startswith(("#", "!")):
            continue
        fields = line.split()
        try:
            values = [fortran_float(field) for field in fields]
        except ValueError as exc:
            raise ValueError(f"{path}:{lineno}: invalid numeric output ({exc})") from exc
        if pressure_header and header is None and not rows and len(values) == 1:
            header = values[0]
            continue
        if len(values) != columns:
            raise ValueError(f"{path}:{lineno}: expected {columns} columns, got {len(values)}")
        rows.append(values)
    if pressure_header and header is None:
        raise ValueError(f"{path}: missing one-column pressure header")
    if len(rows) < minimum_rows:
        raise ValueError(f"{path}: incomplete output ({len(rows)} rows, expected >={minimum_rows})")
    result = np.asarray(rows, dtype=float)
    return (header, result) if pressure_header else result


def first_tokens(path):
    return [fortran_float(line.split()[0]) for line in Path(path).read_text().splitlines()
            if line.strip() and not line.lstrip().startswith(("#", "!"))]


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def completed_reference(directory, model, args):
    directory = Path(directory).resolve()
    manifests = [directory / name for name in ("manifest.json", "run_manifest.json")
                 if (directory / name).is_file()]
    if len(manifests) != 1:
        raise ValueError(f"{directory}: require exactly one manifest.json or run_manifest.json")
    manifest = json.loads(manifests[0].read_text())
    source = manifest.get('source', {})
    if source.get('modified') and source.get('variant') not in ('real-erf', 'real-erf+dense-pdf-output'):
        raise ValueError(f'{directory}: unsupported modified reference source')
    if source.get('file') and source.get('compiled_sha256'):
        if sha256(directory/source['file']) != source['compiled_sha256']:
            raise ValueError(f'{directory}: reference source hash changed')
    run = manifest.get("run", {})
    if (manifest.get("status") != "validated" or run.get("returncode") != 0
            or run.get("timed_out") is not False
            or manifest.get("validation", {}).get("passed") is not True):
        raise ValueError(f"{directory}: reference is not a successfully validated completed run")
    material = first_tokens(directory / MATERIAL_NAME)
    numeric = first_tokens(directory / "IN.mathematical")
    smoothing = first_tokens(directory / "IN_npolyCq.npolyab.npolyE")
    if len(material) < 5 or len(numeric) < 18:
        raise ValueError(f"{directory}: truncated material/numerical inputs")
    if smoothing[0] != 0:
        raise ValueError(f"{directory}: PSD smoothing is active; cannot compare against unsmoothed Python PSD")
    modulus = material[3] / (1.0 - material[2] ** 2)
    for name, source_value, python_value in (
        ("effective modulus", modulus, args.effective_modulus_pa),
        ("gamma", numeric[8], args.gamma), ("afact", numeric[14], args.pdf_afact),
        ("PSD scale", numeric[17], args.scale_c),
    ):
        if not math.isclose(source_value, python_value, rel_tol=1e-12, abs_tol=0):
            raise ValueError(f"{directory}: source {name}={source_value} differs from Python {python_value}")
    reference_psd = read_psd(directory / PSD_NAME, scale_c=args.scale_c)
    if (reference_psd.q0 != model.q0 or reference_psd.q1 != model.q1
            or not np.array_equal(reference_psd.log10_q, model.psd.log10_q)
            or not np.array_equal(reference_psd.log10_c, model.psd.log10_c)):
        raise ValueError(f"{directory}: PSD/cutoffs differ from --psd")
    nn = int(np.log(model.q1 / model.q0) / numeric[0])
    return {
        "directory": str(directory), "manifest": str(manifests[0]),
        "manifest_sha256": sha256(manifests[0]),
        "run_status": manifest["status"], "nominal_pressure_pa": material[0],
        "source": source,
        "code_variant": source.get('variant', 'legacy_manifest_unspecified'),
        "compiler": manifest.get('compiler'),
        "compile_flags": manifest.get('compile_flags'),
        "binary": manifest.get('binary'),
        "validated_outputs": manifest.get('validation', {}).get('outputs'),
        "delzeta": numeric[0], "source_wavevector_nodes": nn,
        "source_actual_max_q_per_m": float(model.q0 * np.exp(nn * numeric[0])),
        "input_sha256": {name: sha256(directory / name) for name in
                         (PSD_NAME, MATERIAL_NAME, "IN.mathematical", "IN.flagparameters")},
    }


def errors(reference, calculated):
    reference, calculated = np.broadcast_arrays(np.asarray(reference), np.asarray(calculated))
    difference = calculated - reference
    eligible = reference != 0
    with np.errstate(over="ignore", divide="ignore", invalid="ignore"):
        relative = difference[eligible] / reference[eligible]
    overflowed = int(np.sum(~np.isfinite(relative)))

    def rms(values):
        maximum = float(np.max(np.abs(values)))
        return maximum * float(np.sqrt(np.mean((values/maximum)**2))) if maximum else 0.0

    return {
        "count": int(reference.size),
        "max_absolute_error": float(np.max(np.abs(difference))),
        "rms_absolute_error": rms(difference),
        "max_relative_error_nonzero_reference": float(np.max(np.abs(relative))) if relative.size and not overflowed else None,
        "rms_relative_error_nonzero_reference": rms(relative) if relative.size and not overflowed else None,
        "unrepresentable_relative_error_count": overflowed,
        "zero_reference_count": int(np.sum(~eligible)),
    }


def csv(output, filename, names, data):
    np.savetxt(output / filename, np.column_stack(data), delimiter=",", header=",".join(names), comments="")


def compare_area(directory, model, args, output):
    provenance = completed_reference(directory, model, args)
    values = read_numeric(Path(directory) / AREA_NAME, 3, minimum_rows=1000)
    p, reference, fitted = values.T
    if np.any(np.diff(p) <= 0) or np.any((reference < 0) | (reference > 1)):
        raise ValueError("invalid reference pressure/standard-area curve")
    calculated = model.area(p)
    csv(output, "area_comparison.csv", ["pressure_pa", "fortran_standard_area", "python_standard_area",
                                         "python_minus_fortran", "fortran_empirical_fit_not_baseline"],
        [p, reference, calculated, calculated-reference, fitted])
    return {"provenance": provenance, "standard_area": errors(reference, calculated)}, (p, reference, calculated)


def compare_mode2(directory, model, args, output):
    provenance = completed_reference(directory, model, args)
    values = read_numeric(Path(directory) / MEAN_NAME, 4, minimum_rows=108)
    p, reference_h, reference_a, fitted = values.T
    if np.any(np.diff(p) <= 0) or np.any(reference_h < 0):
        raise ValueError("invalid reference pressure/mean-gap curve")
    calculated_h = model.mean_gap(p)
    calculated_a = model.area(p)
    csv(output, "mean_gap_comparison.csv", ["pressure_pa", "fortran_mean_gap_m", "python_mean_gap_m",
        "gap_difference_m", "fortran_standard_area", "python_standard_area", "fortran_empirical_fit_not_baseline"],
        [p, reference_h, calculated_h, calculated_h-reference_h, reference_a, calculated_a, fitted])
    return {"provenance": provenance, "mean_gap": errors(reference_h, calculated_h),
            "standard_area": errors(reference_a, calculated_a)}, (p, reference_h, calculated_h)


def compare_mode3(directory, model, args, output):
    provenance = completed_reference(directory, model, args)
    summary = read_numeric(Path(directory) / SUMMARY_NAME, 5)
    if len(summary) != 1:
        raise ValueError("mode3 comparison requires one pressure block; found multiple summaries")
    pressure, source_first, source_h, source_contact_from_pdf, source_contact = summary[0]
    header, values = read_numeric(Path(directory) / PDF_NAME, 2, pressure_header=True, minimum_rows=400)
    if not math.isclose(header, pressure, rel_tol=5e-5) or not math.isclose(
            pressure, provenance["nominal_pressure_pa"], rel_tol=5e-5):
        raise ValueError("PDF pressure header, summary and archived nominal pressure disagree")
    gaps, source_density = values.T
    if np.any(gaps <= 0) or np.any(np.diff(gaps) <= 0):
        raise ValueError("Fortran PDF gap grid must be positive and strictly increasing")
    result = model.gap_distribution(float(pressure), gaps)
    source_sparse_bin = np.diff(gaps) * (source_density[:-1]+source_density[1:]) / 2.0
    python_bin = np.diff(result.continuous_cdf)
    csv(output, "pdf_density_comparison.csv", ["microgap_m", "fortran_nominal_density_per_m",
        "python_nominal_density_per_m", "python_minus_fortran", "python_resolved_nominal_cdf"],
        [gaps, source_density, result.continuous_density_per_m,
         result.continuous_density_per_m-source_density, result.continuous_cdf])
    csv(output, "pdf_bin_comparison.csv", ["left_microgap_m", "right_microgap_m",
        "fortran_sparse_trapezoid_bin_estimate", "python_analytic_cdf_bin_probability",
        "python_minus_sparse_estimate"], [gaps[:-1], gaps[1:], source_sparse_bin,
            python_bin, python_bin-source_sparse_bin])
    sparse_mass = float(trapezoid(source_density, gaps))
    sparse_first = float(trapezoid(gaps*source_density, gaps))
    source_moment_error = (source_first-source_h)/source_h if source_h else None
    notes = [
        "Fortran PDF output is thinned (first point and every 100th internal point); its trapezoid integrals are estimates.",
        "Internal Fortran summary integrals are reported separately and never replaced by sparse-output quadrature.",
        "Both compared densities carry nominal-area probability, not conditional noncontact normalization.",
        "Python bin masses use analytic Gaussian CDF differences; Fortran bin masses are sparse trapezoid estimates.",
        "The Python PDF uses finite-volume magnification intervals; original Fortran uses finite differences and cutoffs.",
        "No density was normalized or clipped; the Python PDF remains experimental.",
        "Agreement of scalar contact area or mean gap does not establish PDF agreement.",
    ]
    if np.any(source_density < 0):
        notes.append("Negative reference density samples are present and retained; investigate before using this PDF.")
    source_internal_mass = 1.0-source_contact_from_pdf
    return {
        "provenance": provenance, "pressure_pa": float(pressure),
        "mean_gap": errors(source_h, result.mean_gap_m),
        "standard_area": errors(source_contact, result.contact_fraction),
        "source_internal_summary": {
            "mean_gap_from_energy_m": float(source_h), "pdf_first_moment_m": float(source_first),
            "pdf_first_moment_relative_error": float(source_moment_error) if source_moment_error is not None else None,
            "contact_from_pdf": float(source_contact_from_pdf), "standard_contact": float(source_contact),
            "continuous_mass": float(source_internal_mass),
            "mass_partition_error_using_standard_contact": float(source_contact + source_internal_mass - 1),
        },
        "source_sparse_output": {"rows": len(values), "trapezoid_continuous_mass": sparse_mass,
            "trapezoid_first_moment_m": sparse_first,
            "trapezoid_mass_minus_internal_mass": sparse_mass-source_internal_mass,
            "negative_density_samples": int(np.sum(source_density < 0))},
        "python_pdf_diagnostics": asdict(result.diagnostics),
        "nominal_density_errors_on_source_grid": errors(source_density, result.continuous_density_per_m),
        "bin_error_against_sparse_fortran_estimate": errors(source_sparse_bin, python_bin),
        "notes": notes,
    }, (gaps, source_density, result.continuous_density_per_m, source_sparse_bin, python_bin)


def plot_comparisons(output, area, mean, pdf):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(2, 2, figsize=(11, 8), layout="constrained")
    if area is not None:
        p, original, python = area
        axes[0, 0].plot(p, original, label="Fortran standard")
        axes[0, 0].plot(p, python, "--", label="Python standard")
        axes[0, 0].set(xlabel="Pressure (Pa)", ylabel="A/A0", title="Contact area")
    if mean is not None:
        p, original, python = mean
        axes[0, 1].semilogx(p, original*1e6, label="Fortran")
        axes[0, 1].semilogx(p, python*1e6, "--", label="Python")
        axes[0, 1].set(xlabel="Pressure (Pa)", ylabel="Mean gap (micrometres)", title="Mean separation")
    if pdf is not None:
        gaps, original, python, source_bin, python_bin = pdf
        axes[1, 0].semilogx(gaps*1e6, original, label="Fortran sparse samples")
        axes[1, 0].semilogx(gaps*1e6, python, "--", label="Python experimental")
        axes[1, 0].set(xlabel="Microgap (micrometres)", ylabel="Nominal density (1/m)", title="Gap PDF: experimental")
        middle = np.sqrt(gaps[:-1]*gaps[1:])*1e6
        axes[1, 1].semilogx(middle, source_bin, label="Fortran sparse trapezoid")
        axes[1, 1].semilogx(middle, python_bin, "--", label="Python analytic CDF bins")
        axes[1, 1].set(xlabel="Microgap (micrometres)", ylabel="Interval probability", title="Probability bins")
    for ax in axes.flat:
        if ax.lines:
            ax.grid(alpha=.25)
            ax.legend()
        else:
            ax.text(.5, .5, "No completed reference supplied", ha="center", va="center", transform=ax.transAxes)
            ax.set_axis_off()
    fig.suptitle("Archived input parameters; not calibrated experiment. PDF validation remains provisional.")
    fig.savefig(output / "fortran_comparison.png", dpi=180)
    plt.close(fig)


def parser():
    result = argparse.ArgumentParser(description=__doc__)
    result.add_argument("--mode2-dir", type=Path)
    result.add_argument("--mode3-dir", type=Path)
    result.add_argument("--area-dir", type=Path)
    result.add_argument("--psd", type=Path)
    result.add_argument("--output", type=Path, help="New output directory; existing paths are refused")
    result.add_argument("--effective-modulus-pa", type=float, default=4e6)
    result.add_argument("--gamma", type=float, default=.4)
    result.add_argument("--n-q", type=int, default=256)
    result.add_argument("--pdf-afact", type=float, default=1.)
    result.add_argument("--scale-c", type=float, default=1.)
    result.add_argument("--plot", action="store_true")
    result.add_argument("--self-test", action="store_true", help="Run parser/completeness tests only; no physics workload")
    return result


class ParserChecks(unittest.TestCase):
    def test_fortran_exponents(self):
        for text, expected in (("1.23D+04", 12300), ("-.1234-102", -1.234e-103),
                               ("0.5432+101", 5.432e100), ("1.0", 1.)):
            self.assertAlmostEqual(fortran_float(text)/expected, 1.)
        for text in ("***********", "NaN", "Infinity"):
            with self.assertRaises(ValueError):
                fortran_float(text)

    def test_pressure_header_and_incomplete_output(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "pdf"
            path.write_text("1.0D+5\n1.0-100 2.0+100\n2.0-100 3.0+100\n")
            pressure, data = read_numeric(path, 2, pressure_header=True, minimum_rows=2)
            self.assertEqual(pressure, 1e5)
            self.assertEqual(data.shape, (2, 2))
            with self.assertRaises(ValueError):
                read_numeric(path, 2, pressure_header=True, minimum_rows=3)

    def test_zero_reference_metrics_do_not_invent_relative_agreement(self):
        result = errors([0., 1.], [1., 2.])
        self.assertEqual(result["zero_reference_count"], 1)
        self.assertEqual(result["max_absolute_error"], 1.)
        self.assertEqual(result["max_relative_error_nonzero_reference"], 1.)
        self.assertIsNone(errors([0.], [1.])["max_relative_error_nonzero_reference"])
        huge = errors([1e-300], [1e100])
        self.assertEqual(huge["unrepresentable_relative_error_count"], 1)
        self.assertIsNone(huge["max_relative_error_nonzero_reference"])
        self.assertTrue(math.isfinite(errors([1.], [1e300])["rms_absolute_error"]))

    def test_pending_reference_is_rejected_before_reading_outputs(self):
        with tempfile.TemporaryDirectory() as directory:
            (Path(directory)/"run_manifest.json").write_text(json.dumps({"status": "running"}))
            with self.assertRaisesRegex(ValueError, "not a successfully validated completed run"):
                completed_reference(directory, None, None)


def main(argv=None):
    arg_parser = parser()
    args = arg_parser.parse_args(argv)
    if args.self_test:
        suite = unittest.defaultTestLoader.loadTestsFromTestCase(ParserChecks)
        return 0 if unittest.TextTestRunner(verbosity=2).run(suite).wasSuccessful() else 1
    if not args.psd or not args.output or not any((args.area_dir, args.mode2_dir, args.mode3_dir)):
        arg_parser.error("provide --psd, --output and at least one completed reference directory")
    if args.output.exists():
        arg_parser.error("--output must be a new directory")
    try:
        spectrum = read_psd(args.psd, scale_c=args.scale_c)
        model = PerssonModel(spectrum, args.effective_modulus_pa, args.gamma,
                             n_q=args.n_q, pdf_afact=args.pdf_afact)
        # Check every run before creating output or doing expensive PDF work.
        for directory in (args.area_dir, args.mode2_dir, args.mode3_dir):
            if directory:
                completed_reference(directory, model, args)
        args.output.mkdir(parents=True, exist_ok=False)
        summary = {
            "label": "Archived parameters; not a calibrated experiment",
            "scope": "Standard elastic dry area, energy mean separation, experimental source-inspired gap PDF",
            "parameters": {key: getattr(args, key) for key in
                           ("effective_modulus_pa", "gamma", "n_q", "pdf_afact", "scale_c")},
            "psd": {"path": str(args.psd.resolve()), "sha256": sha256(args.psd)},
            "caveats": ["Fortran E11.4 output rounding limits pointwise comparison precision.",
                "Source right-rectangle log-wavevector integration differs from Python exact PSD moments and segmented quadrature.",
                "Source upper wavevector is truncated to an integer number of delzeta intervals.",
                "No arbitrary agreement threshold is applied; scalar and probability errors are reported separately."],
            "comparisons": {},
        }
        area_plot = mean_plot = pdf_plot = None
        if args.area_dir:
            summary["comparisons"]["area"], area_plot = compare_area(args.area_dir, model, args, args.output)
        if args.mode2_dir:
            summary["comparisons"]["mode2"], mean_plot = compare_mode2(args.mode2_dir, model, args, args.output)
        if args.mode3_dir:
            summary["comparisons"]["mode3"], pdf_plot = compare_mode3(args.mode3_dir, model, args, args.output)
        (args.output / "comparison.json").write_text(json.dumps(summary, indent=2, allow_nan=False)+"\n")
        if args.plot:
            plot_comparisons(args.output, area_plot, mean_plot, pdf_plot)
        print(f"Comparison saved to {args.output.resolve()}; inspect scalar and PDF diagnostics separately.")
        return 0
    except (OSError, ValueError, TypeError, KeyError, IndexError) as exc:
        print(f"Comparison failed: {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
