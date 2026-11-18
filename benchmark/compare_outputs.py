"""Compare lower-plane primary outputs; keep native output defects visible."""
import argparse
import csv
import json
import math
import re
from pathlib import Path
import numpy as np
from openpyxl import Workbook

HERE = Path(__file__).resolve().parent

def native_number(value):
    if isinstance(value, bytes):
        value = value.decode('ascii')
    value = value.replace('D', 'E')
    if set(value) == {'*'}:
        return float('nan')
    # E11.4 drops the E for a three-digit exponent, e.g. -0.4346-114.
    if 'E' not in value.upper():
        value = re.sub(r'(?<=\d)([+-]\d{3})$', r'E\1', value)
    return float(value)

def one_file(directory, pattern):
    matches = list(directory.glob(pattern))
    if len(matches) != 1:
        raise ValueError(f'Expected one {pattern}, found {len(matches)}')
    return matches[0]

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--python', type=Path, required=True)
    parser.add_argument('--fortran', type=Path, required=True)
    parser.add_argument('--fortran-kind', choices=['pristine', 'real-erf-diagnostic'], required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    cfg = json.loads((HERE / 'manifest.json').read_text())
    nx = int((args.fortran / 'IN.MIXEDLUBRICATION.mathematical').read_text().splitlines()[1].split()[0])
    tau = float((args.fortran / 'IN.MIXEDLUBRICATION').read_text().splitlines()[4].split()[0])
    if not math.isclose(tau, cfg['physical']['tau_pa']):
        raise ValueError('Fortran tau differs from the benchmark')
    f = np.loadtxt(one_file(args.fortran, 'out.MSLUBE.LOW.WITH.FrictionFlowfactor.*'), ndmin=2, converters=native_number)
    profiles = np.loadtxt(one_file(args.fortran, 'out.MSLUBE.1time.2vel.3x.*'), ndmin=2, converters=native_number)
    gaps = np.loadtxt(one_file(args.fortran, 'out.MSLUBE.1time.2vel.3minimumSeparation.*'), ndmin=2, converters=native_number)
    if f.shape[1] != 16 or not np.isfinite(f[:, :8]).all():
        raise ValueError('Fortran primary columns missing/nonfinite')
    with (args.python / 'mixed_cases.csv').open(newline='') as file:
        p = list(csv.DictReader(file))
    if len(p) != len(f):
        raise ValueError('Speed counts differ; partial outputs are not a complete comparison')
    headers = ['speed_m_s', 'Ws_prime_n_m', 'Wf_prime_n_m', 'Ar_over_B_m',
               'hmin_effective_profile_m', 'pfluid_max_pa', 'mu_solid_resisting',
               'mu_fluid_sliding_resisting', 'mu_fluid_pressure_resisting', 'mu_total_resisting']
    wb = Workbook()
    fp = wb.active
    fp.title = 'Python'
    ff = wb.create_sheet('Fortran_' + ('pristine' if args.fortran_kind == 'pristine' else 'diagnostic'))
    diff = wb.create_sheet('相对差异')
    for ws in (fp, ff):
        ws.append(headers)
    diff.append(['speed_m_s'] + ['Python_over_Fortran_minus_1__' + k for k in headers[1:]])
    issues = wb.create_sheet('原版导出问题')
    issues.append(['speed_m_s', 'native_hmin_file_m', 'recovered_first_Nx_hmin_m', 'nonfinite_auxiliary_columns'])
    records = []
    for py, fr in zip(p, f):
        v, load = float(py['lower_wall_velocity_m_s']), fr[1]
        if v <= 0 or not math.isclose(v, fr[0], rel_tol=1e-3):
            raise ValueError('This converter requires matched positive lower-plane velocities')
        # Native elastic low-speed loop computes Nx samples, while its export writes 2*Nx+1.
        block = profiles[np.isclose(profiles[:, 0], v, rtol=1e-3, atol=0)][:nx]
        if len(block) != nx or not np.isfinite(block).all() or not (block[:, 2] > 0).all():
            raise ValueError('Invalid physical-domain native profile; inspect original output')
        pyload = float(py['contact_line_load_n_m']) + float(py['fluid_line_load_n_m'])
        pr = [v, float(py['contact_line_load_n_m']), float(py['fluid_line_load_n_m']),
              float(py['real_contact_area_per_width_m']), float(py['minimum_mean_gap_m']),
              float(py['maximum_fluid_pressure_pa']),
              *(-float(py[k]) / pyload for k in ('lower_solid_force_n_m',
                 'lower_fluid_sliding_force_n_m', 'lower_fluid_pressure_flow_force_n_m')),
              float(py['lower_resisting_friction_coefficient'])]
        nr = [v, fr[2], fr[3], abs(fr[5]) * load / tau, float(block[:, 2].min()),
              float(block[:, 4].max()), -fr[5], -fr[6], -fr[7], -float(fr[5:8].sum())]
        fp.append(pr)
        ff.append(nr)
        diff.append([v] + [x / y - 1 if y else None for x, y in zip(pr[1:], nr[1:])])
        gap = gaps[np.isclose(gaps[:, 0], v, rtol=1e-3, atol=0)]
        issues.append([v, float(gap[0, 1]) if len(gap) == 1 else None, nr[4],
                       ','.join(str(i+1) for i, x in enumerate(fr) if not np.isfinite(x))])
        records.append({'speed_m_s': v, 'python': dict(zip(headers, pr)),
                        'fortran': dict(zip(headers, nr))})
    notes = wb.create_sheet('口径与限制')
    notes.append(['Fortran source kind', args.fortran_kind])
    notes.append(['Python snapshot', cfg['python_snapshot']])
    notes.append(['Fortran expected pristine source SHA256', cfg['fortran_source_sha256']])
    notes.append(['Area definition', 'Python integral of local A(pc); Fortran native mean-pressure estimate.'])
    notes.append(['hmin', 'Fortran recovered from first Nx physical samples; native hmin retained separately.'])
    notes.append(['Scope', 'Selected moving lower plane only. Auxiliary native NaNs are not valid forces.'])
    notes.append(['Validation', 'Output comparison does not by itself prove mesh/domain/constitutive convergence.'])
    for ws in (fp, ff, diff, issues):
        ws.freeze_panes = 'A2'
        ws.auto_filter.ref = ws.dimensions
    args.output.mkdir(parents=True, exist_ok=False)
    wb.save(args.output / 'comparison.xlsx')
    (args.output / 'comparison.json').write_text(json.dumps({
        'fortran_kind': args.fortran_kind, 'convergence_certified': False,
        'headers': headers, 'records': records}, indent=2) + '\n')
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(2, 2, figsize=(10, 7))
    for ax, key, ylabel in zip(axes.flat, ['mu_total_resisting', 'Ar_over_B_m', 'Wf_prime_n_m',
                                         'hmin_effective_profile_m'],
                               ['Lower-plane total mu', 'Real area / B [m]', 'Fluid load [N/m]', 'Mean gap min [m]']):
        for source, style in [('python', 'o-'), ('fortran', 's--')]:
            ax.semilogx([r['speed_m_s'] for r in records], [r[source][key] for r in records],
                        style, label=source + (' (' + args.fortran_kind + ')' if source == 'fortran' else ''))
        ax.set_xlabel('Speed [m/s]')
        ax.set_ylabel(ylabel)
        ax.grid(alpha=.3)
        ax.legend(fontsize=8)
    fig.suptitle('Common synthetic PSD — preliminary comparison, no fitted tau')
    fig.tight_layout()
    fig.savefig(args.output / 'comparison.png', dpi=160)
    print(args.output / 'comparison.xlsx')

if __name__ == '__main__':
    main()
