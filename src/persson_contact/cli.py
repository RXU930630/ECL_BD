"""Run a dry-contact case from an explicit JSON configuration."""
import argparse
import csv
from dataclasses import asdict, is_dataclass
import hashlib
import json
from pathlib import Path
import textwrap

import numpy as np

from .config import load_problem
from .core import PerssonModel
from .table import build_contact_table
from . import __version__


def _jsonable(value):
    if is_dataclass(value):
        return _jsonable(asdict(value))
    if isinstance(value, dict):
        return {str(k): _jsonable(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [_jsonable(v) for v in value]
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, np.generic):
        return value.item()
    return value


def run_contact(problem, options, output, *, plot=False, compute_gap_pdf=True, provenance=None):
    """Shared contact calculation/output entry for Python and JSON workflows."""
    output = Path(output)
    model = PerssonModel(problem.psd, problem.effective_modulus_pa,
                         gamma=options.gamma, n_q=options.n_q,
                         pdf_afact=options.pdf_afact)
    table = build_contact_table(model, problem.pressures_pa)
    output.mkdir(parents=True, exist_ok=True)
    if any(output.iterdir()):
        raise FileExistsError(f'输出目录非空，请使用新的 output_directory: {output}')
    np.savetxt(output / 'contact_table.csv',
               np.column_stack((table.pressures_pa, table.mean_gaps_m, table.contact_fractions)),
               delimiter=',', header='nominal_contact_pressure_Pa,mean_gap_m,relative_contact_area', comments='')
    critical = [model.critical_separation(float(p), options.critical_relative_contact_area,
                                         options.critical_area_tolerance) for p in table.pressures_pa]
    with (output / 'critical_separation.csv').open('w', newline='', encoding='utf-8') as stream:
        records = [asdict(value) for value in critical]
        writer = csv.DictWriter(stream, fieldnames=records[0].keys())
        writer.writeheader()
        writer.writerows(records)
    hrms = float(problem.psd.hrms)
    smax = max(20 * hrms, 10 * float(np.max(table.mean_gaps_m)))
    smin = max(hrms * 1e-7, np.finfo(float).tiny)
    gaps = np.r_[0.0, np.geomspace(smin, smax, options.n_gap - 1)]
    pdf_diagnostics = []
    for index, pressure in enumerate(table.pressures_pa if compute_gap_pdf else []):
        if table.contact_fractions[index] == 1.0:
            pdf_diagnostics.append({'pressure_pa': float(pressure), 'status': 'full_contact_to_floating_point_precision; conditional positive-gap density undefined'})
            continue
        # Density is defined/sampled for u>0. Evaluate CDF there too, then add
        # its known zero origin for bin probabilities without sampling a
        # potentially singular density at u=0.
        distribution = model.gap_distribution(float(pressure), gaps[1:])
        continuous_cdf = np.r_[0., distribution.continuous_cdf]
        positive_cdf = np.r_[0., distribution.positive_cdf]
        filename = f'gap_pdf_{index:03d}.csv'
        # The Dirac atom is separate, and arbitrarily narrow positive-gap
        # components make density(0) unsuitable as a tabulated point sample.
        np.savetxt(output / filename,
                   np.column_stack((distribution.gaps_m, distribution.continuous_density_per_m,
                                    distribution.positive_density_per_m,
                                    distribution.continuous_cdf, distribution.positive_cdf)),
                   delimiter=',', header='microgap_m,nominal_area_density_per_m,conditional_positive_density_per_m,resolved_nominal_area_cdf,resolved_conditional_positive_cdf', comments='')
        bin_filename = f'gap_bins_{index:03d}.csv'
        np.savetxt(output / bin_filename,
                   np.column_stack((gaps[:-1], gaps[1:], np.diff(continuous_cdf),
                                    np.diff(positive_cdf),
                                    np.diff(continuous_cdf) / np.diff(gaps),
                                    np.diff(positive_cdf) / np.diff(gaps))),
                   delimiter=',', header='left_microgap_m,right_microgap_m,nominal_area_probability,conditional_positive_probability,bin_average_nominal_density_per_m,bin_average_conditional_density_per_m', comments='')
        pdf_diagnostics.append({'pressure_pa': float(pressure), 'relative_contact_area': float(distribution.contact_fraction),
                                'file': filename, 'bin_file': bin_filename,
                                'exported_bin_mass': float(continuous_cdf[-1]),
                                'diagnostics': _jsonable(distribution.diagnostics)})
    summary = {
        'description': problem.description,
        **(provenance or {}),
        'package_version': __version__,
        'model': 'standard_Persson_area_gamma_energy_reflected_Gaussian_FV',
        'scope': 'nonadhesive, linearly elastic, isotropic 2D PSD; local half-space response',
        'validation_limit': 'Full PDF comparison against Fortran is pending; gap PDF is experimental, with explicit mass and first-moment diagnostics.',
        'probability_representation': 'Use CDF-difference gap_bins probabilities for integration. Pointwise density CSV excludes s=0 and can under-resolve narrow peaks; sampled/trapezoid diagnostics are not probability-normalization certificates.',
        'effective_modulus_pa': float(problem.effective_modulus_pa),
        'psd_hrms_m': hrms, 'psd_rms_slope': float(problem.psd.rms_slope),
        'psd_cutoffs_m_inverse': {'q0': problem.psd.q0, 'qr': problem.psd.qr, 'q1': problem.psd.q1},
        'tail_is_extrapolated': problem.psd_tail_is_extrapolated,
        'options': _jsonable(options),
        'compute_gap_pdf': compute_gap_pdf,
        'valid_mean_gap_range_m': table.valid_gap_range_m,
        'pressure_scan': {'minimum_pa': float(table.pressures_pa[0]),
                          'maximum_pa': float(table.pressures_pa[-1]),
                          'number_of_points': len(table.pressures_pa),
                          'reference_pa': problem.reference_pressure_pa,
                          'reference_index': next((i for i, p in enumerate(table.pressures_pa)
                                                   if p == problem.reference_pressure_pa), None),
                          'pressures_pa': table.pressures_pa.tolist()},
        'area_convention': 'A is real contact area; A/A0 is relative contact area. No dimensional A is calculated without nominal area A0.',
        'critical_separation': {'definition': 'uc=u1 near critical relative contact area; alpha=1',
                                'target_relative_contact_area': options.critical_relative_contact_area,
                                'absolute_area_tolerance': options.critical_area_tolerance,
                                'status_counts': {status: sum(value.status == status for value in critical)
                                                  for status in sorted({value.status for value in critical})},
                                'reference_result': next((asdict(value) for value in critical
                                                          if value.pressure_pa == problem.reference_pressure_pa), None)},
        'pdf_diagnostics': pdf_diagnostics,
    }
    (output / 'diagnostics.json').write_text(json.dumps(summary, indent=2, allow_nan=False) + '\n', encoding='utf-8')
    if plot:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        fig, axes = plt.subplots(1, 4, figsize=(16, 4), layout='constrained')
        if problem.description:
            fig.suptitle(textwrap.fill(problem.description, width=130), fontsize=10)
        axes[0].plot(problem.psd.log10_q, problem.psd.log10_c, '.-')
        axes[0].set(xlabel='log10(q / m$^{-1}$)', ylabel='log10(C / m$^4$)', title='Input PSD')
        axes[1].semilogx(table.pressures_pa / 1e6, table.contact_fractions, 'o-')
        axes[1].set(xlabel='Nominal dry pressure (MPa)', ylabel='Relative contact area A/A0', title='Relative contact area')
        axes[2].semilogx(table.pressures_pa / 1e6, table.mean_gaps_m * 1e6, 'o-')
        axes[2].set(xlabel='Nominal dry pressure (MPa)', ylabel=r'Mean gap ($\mu$m)', title='Mean separation')
        axes[3].semilogx(table.pressures_pa/1e6,
                        [value.critical_separation_m*1e6 if value.critical_separation_m is not None else np.nan
                         for value in critical], 'o-')
        axes[3].set(xlabel='Nominal dry pressure (MPa)', ylabel=r'$u_c=u_1$ ($\mu$m)',
                    title=f'Critical separation (A/A0 ~ {options.critical_relative_contact_area:g})')
        if problem.reference_pressure_pa is not None:
            for ax in axes[1:]:
                ax.axvline(problem.reference_pressure_pa/1e6, color='grey', linestyle=':', label='Reference pressure')
                ax.legend()
        for ax in axes:
            ax.grid(alpha=0.25)
        fig.savefig(output / 'dry_contact.png', dpi=180)
        plt.close(fig)
    print(f'接触表与诊断已保存到 {output.resolve()}')
    print(f'压力范围 {table.pressures_pa[0]:g}–{table.pressures_pa[-1]:g} Pa，'
          f'共 {len(table.pressures_pa)} 点；A/A0 为 relative contact area。')
    if compute_gap_pdf:
        print('间隙分布仍是实验版，请查看概率与第一矩诊断。')
    return summary


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('config', type=Path, help='JSON case file; SI units')
    parser.add_argument('--output', type=Path, default=Path('results/dry_contact'))
    parser.add_argument('--plot', action='store_true', help='Save a PNG; requires matplotlib')
    parser.add_argument('--no-pdf', action='store_true', help='Skip experimental gap distributions')
    args = parser.parse_args(argv)
    problem, options = load_problem(args.config)
    config_bytes = args.config.read_bytes()
    config_raw = json.loads(config_bytes.decode('utf-8-sig'))
    psd_path = Path(config_raw['psd']['path'])
    if not psd_path.is_absolute():
        psd_path = args.config.parent / psd_path
    run_contact(problem, options, args.output, plot=args.plot, compute_gap_pdf=not args.no_pdf,
                provenance={'config': str(args.config.resolve()),
                            'config_sha256': hashlib.sha256(config_bytes).hexdigest(),
                            'psd_sha256': hashlib.sha256(psd_path.read_bytes()).hexdigest(),
                            'material_input': config_raw['material']})
    return 0
