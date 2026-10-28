"""Run a dry-contact case from an explicit JSON configuration."""
import argparse
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


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('config', type=Path, help='JSON case file; SI units')
    parser.add_argument('--output', type=Path, default=Path('results/dry_contact'))
    parser.add_argument('--plot', action='store_true', help='Save a PNG; requires matplotlib')
    args = parser.parse_args(argv)
    problem, options = load_problem(args.config)
    config_bytes = args.config.read_bytes()
    config_raw = json.loads(config_bytes.decode('utf-8-sig'))
    psd_path = Path(config_raw['psd']['path'])
    if not psd_path.is_absolute():
        psd_path = args.config.parent / psd_path
    model = PerssonModel(problem.psd, problem.effective_modulus_pa,
                         gamma=options.gamma, n_q=options.n_q,
                         pdf_afact=options.pdf_afact)
    table = build_contact_table(model, problem.pressures_pa)
    args.output.mkdir(parents=True, exist_ok=True)
    np.savetxt(args.output / 'contact_table.csv',
               np.column_stack((table.pressures_pa, table.mean_gaps_m, table.contact_fractions)),
               delimiter=',', header='nominal_contact_pressure_Pa,mean_gap_m,contact_area_fraction', comments='')
    hrms = float(problem.psd.hrms)
    smax = max(20 * hrms, 10 * float(np.max(table.mean_gaps_m)))
    smin = max(hrms * 1e-7, np.finfo(float).tiny)
    gaps = np.r_[0.0, np.geomspace(smin, smax, options.n_gap - 1)]
    pdf_diagnostics = []
    for index, pressure in enumerate(table.pressures_pa):
        if table.contact_fractions[index] == 1.0:
            pdf_diagnostics.append({'pressure_pa': float(pressure), 'status': 'full_contact_to_floating_point_precision; conditional positive-gap density undefined'})
            continue
        distribution = model.gap_distribution(float(pressure), gaps)
        filename = f'gap_pdf_{index:03d}.csv'
        # The Dirac atom is separate, and arbitrarily narrow positive-gap
        # components make density(0) unsuitable as a tabulated point sample.
        positive = gaps > 0
        np.savetxt(args.output / filename,
                   np.column_stack((gaps[positive], distribution.continuous_density_per_m[positive],
                                    distribution.positive_density_per_m[positive],
                                    distribution.continuous_cdf[positive], distribution.positive_cdf[positive])),
                   delimiter=',', header='microgap_m,nominal_area_density_per_m,conditional_positive_density_per_m,resolved_nominal_area_cdf,resolved_conditional_positive_cdf', comments='')
        bin_filename = f'gap_bins_{index:03d}.csv'
        np.savetxt(args.output / bin_filename,
                   np.column_stack((gaps[:-1], gaps[1:], np.diff(distribution.continuous_cdf),
                                    np.diff(distribution.positive_cdf),
                                    np.diff(distribution.continuous_cdf) / np.diff(gaps),
                                    np.diff(distribution.positive_cdf) / np.diff(gaps))),
                   delimiter=',', header='left_microgap_m,right_microgap_m,nominal_area_probability,conditional_positive_probability,bin_average_nominal_density_per_m,bin_average_conditional_density_per_m', comments='')
        pdf_diagnostics.append({'pressure_pa': float(pressure), 'contact_fraction': float(distribution.contact_fraction),
                                'file': filename, 'bin_file': bin_filename,
                                'diagnostics': _jsonable(distribution.diagnostics)})
    summary = {
        'config': str(args.config.resolve()),
        'description': problem.description,
        'config_sha256': hashlib.sha256(config_bytes).hexdigest(),
        'psd_sha256': hashlib.sha256(psd_path.read_bytes()).hexdigest(),
        'package_version': __version__,
        'model': 'standard_Persson_area_gamma_energy_reflected_Gaussian_FV',
        'scope': 'nonadhesive, linearly elastic, isotropic 2D PSD; local half-space response',
        'validation_limit': 'Original Fortran has not been run; gap PDF is an experimental source-inspired approximation, with explicit mass and first-moment diagnostics.',
        'probability_representation': 'Use CDF-difference gap_bins probabilities for integration. Pointwise density CSV excludes s=0 and can under-resolve narrow peaks; sampled/trapezoid diagnostics are not probability-normalization certificates.',
        'effective_modulus_pa': float(problem.effective_modulus_pa),
        'psd_hrms_m': hrms, 'psd_rms_slope': float(problem.psd.rms_slope),
        'psd_cutoffs_m_inverse': {'q0': problem.psd.q0, 'qr': problem.psd.qr, 'q1': problem.psd.q1},
        'tail_is_extrapolated': problem.psd_tail_is_extrapolated,
        'options': _jsonable(options),
        'valid_mean_gap_range_m': table.valid_gap_range_m,
        'pdf_diagnostics': pdf_diagnostics,
    }
    (args.output / 'diagnostics.json').write_text(json.dumps(summary, indent=2, allow_nan=False) + '\n', encoding='utf-8')
    if args.plot:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        fig, axes = plt.subplots(1, 3, figsize=(13, 4), layout='constrained')
        if problem.description:
            fig.suptitle(textwrap.fill(problem.description, width=130), fontsize=10)
        axes[0].plot(problem.psd.log10_q, problem.psd.log10_c, '.-')
        axes[0].set(xlabel='log10(q / m$^{-1}$)', ylabel='log10(C / m$^4$)', title='Input PSD')
        axes[1].semilogx(table.pressures_pa / 1e6, table.contact_fractions, 'o-')
        axes[1].set(xlabel='Nominal dry pressure (MPa)', ylabel='A/A0', title='Standard Persson area')
        axes[2].semilogx(table.pressures_pa / 1e6, table.mean_gaps_m * 1e6, 'o-')
        axes[2].set(xlabel='Nominal dry pressure (MPa)', ylabel=r'Mean gap ($\mu$m)', title='Mean separation')
        for ax in axes:
            ax.grid(alpha=0.25)
        fig.savefig(args.output / 'dry_contact.png', dpi=180)
        plt.close(fig)
    print(f'Contact table and diagnostics saved to {args.output.resolve()}')
    print('Material and numerical values in the bundled case are illustrative, not calibrated research parameters.')
    print('Inspect gap-PDF moment diagnostics; Fortran comparison is pending.')
    return 0
