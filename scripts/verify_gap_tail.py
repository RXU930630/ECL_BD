#!/usr/bin/env python3
"""Verify repaired PDF tails against completed, unchanged Fortran benchmarks.

Requires raw reference outputs (not shipped). No reference simulation is rerun.
The compact fixture and exported plot samples retain real Fortran values.
"""
import argparse
from dataclasses import asdict
import json
from pathlib import Path
from types import SimpleNamespace

import numpy as np
from scipy.integrate import cumulative_trapezoid, trapezoid

import compare_fortran_reference as reference
from persson_contact.core import PerssonModel
from persson_contact.psd import read_psd


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--matrix', required=True)
    parser.add_argument('--psd', required=True)
    parser.add_argument('--output', required=True)
    parser.add_argument('--fixture', help='Write compact, provenance-stamped regression data')
    args = parser.parse_args()
    output = Path(args.output)
    output.mkdir(parents=True, exist_ok=False)
    psd = read_psd(args.psd)
    records, curves, fixture = [], [], []
    for manifest_path in sorted(Path(args.matrix).glob('*/manifest.json')):
        manifest = json.loads(manifest_path.read_text())
        case = manifest.get('case', {})
        if case.get('mode') != 3 or manifest.get('status') != 'validated':
            continue
        model = PerssonModel(psd, case['effective_modulus_pa'], .4,
                             pdf_log_step=case['delzeta'])
        provenance = reference.completed_reference(manifest_path.parent, model, SimpleNamespace(
            effective_modulus_pa=case['effective_modulus_pa'], gamma=.4, pdf_afact=1., scale_c=1.))
        pressure, data = reference.read_numeric(manifest_path.parent/reference.PDF_NAME,
            2, pressure_header=True, minimum_rows=40000)
        if pressure != case['pressure_pa']:
            raise ValueError('Pressure header differs from manifest')
        gaps, density = data.T
        python = model.gap_distribution(pressure, gaps)
        calculated = python.continuous_density_per_m
        cdf = cumulative_trapezoid(density, gaps, initial=0.)
        pcdf = python.continuous_cdf - python.continuous_cdf[0]
        tail = gaps >= gaps[np.argmax(gaps*density)]
        nonzero = density > 0
        relative = np.abs(calculated[nonzero]/density[nonzero]-1)
        summary = reference.read_numeric(manifest_path.parent/reference.SUMMARY_NAME, 5)[0]
        worst = int(np.flatnonzero(nonzero)[np.argmax(relative)])
        record = {'worst_density_error_gap_m': float(gaps[worst]),
                  'worst_density_error_reference_per_m': float(density[worst]),
                  'tail_definition': 'gaps >= argmax(u * reference_density)',
                  'name': manifest_path.parent.name, **case, 'provenance': provenance,
                  'points': len(gaps), 'max_pointwise_relative_density_error': float(max(relative)),
                  'tail_max_pointwise_relative_density_error': float(np.max(abs(calculated[tail]/density[tail]-1))),
                  'cdf_max_absolute_error': float(max(abs(pcdf-cdf))),
                  'fortran_first_moment_m': float(summary[1]),
                  'python_full_axis_first_moment_m': python.diagnostics.analytic_first_moment_m,
                  'full_first_moment_relative_error': float(python.diagnostics.analytic_first_moment_m/summary[1]-1),
                  'common_first_moment_relative_error': float(trapezoid(gaps*calculated,gaps)/trapezoid(gaps*density,gaps)-1),
                  'min_reference_density_per_m': float(min(density)),
                  'diagnostics': asdict(python.diagnostics)}
        if manifest_path.parent.name == 'E3_p10k_mode3_d01':
            record['diagnostic_source_u1_control'] = source_u1_control(
                manifest_path.parent, model, pressure, gaps, density)
        records.append(record)
        indices = np.unique(np.r_[np.linspace(0,len(gaps)-1,450).astype(int),
                                 np.arange(len(gaps)-800,len(gaps),8)])
        reference.csv(output, record['name']+'.csv',
            ['microgap_m','fortran_density_per_m','python_density_per_m','fortran_common_cdf','python_common_cdf'],
            [gaps[indices],density[indices],calculated[indices],cdf[indices],pcdf[indices]])
        curves.append((case, gaps[indices], density[indices], calculated[indices]))
        take = np.unique(np.r_[np.linspace(0,len(gaps)-1,40).astype(int),
                               np.arange(len(gaps)-1000,len(gaps),50),len(gaps)-1])
        fixture.append({'name': record['name'], 'case': case,
                        'source_sha256': manifest['source']['original_sha256'],
                        'reference_file_sha256': reference.sha256(manifest_path.parent/reference.PDF_NAME),
                        'gaps_m': gaps[take].tolist(), 'density_per_m': density[take].tolist(),
                        'first_moment_m':float(summary[1])})
        print(record['name'], 'max PDF relative error',record['max_pointwise_relative_density_error'],
              'first moment',record['full_first_moment_relative_error'],flush=True)
    if len(records) != 10:
        raise ValueError(f'Expected all 10 completed PDF cases, found {len(records)}')
    (output/'comparison.json').write_text(json.dumps({'case_count': len(records), 'cases': records,
        'psd_sha256':reference.sha256(args.psd),
        'python_core_sha256':reference.sha256(Path(__file__).resolve().parents[1]/'src/persson_contact/core.py'),
        'reference_policy':'Reuse completed Fortran runs; preserve raw output hashes. No fitting or normalization.',
        'pdf_policy':'Method-1 scale quadrature and cutoffs; accurate separation kernel; match delzeta per case.'},
        indent=2,allow_nan=False)+'\n')
    if args.fixture:
        target=Path(args.fixture);target.parent.mkdir(parents=True,exist_ok=True)
        target.write_text(json.dumps({'psd_sha256':reference.sha256(args.psd),
                                     'cases':fixture},indent=2)+'\n')
    plot(output,curves)


def source_u1_control(directory, model, pressure, gaps, reference_density):
    """Diagnostic only: isolate u1 from weights/widths/Gaussian evaluation.

    The production model does not read these reference u1 values. Replacing
    ONLY its centers in this control locates the remaining density difference.
    """
    from scipy.special import erf
    q, c, increments, i3, variance = model._pdf_geometry
    area = erf(pressure/(model.effective_modulus_pa*np.sqrt(np.pi*i3)))
    centers = reference.read_numeric(directory/'out.1u1.2u1min', 2)
    centers = centers[centers[:,0] >= model.pdf_min_u1_m, 0]
    scales = reference.read_numeric(directory/
        'out.1zeta.2du.over.dmu.3Area.4dA.over.dmu.5dlogA.over.dmu.6separation.7separation1',7)
    indices = np.rint(np.log(scales[:,0])/model.pdf_log_step).astype(int)-1
    if len(indices) != len(centers):
        raise ValueError('Diagnostic source u1 and scale output lengths differ')
    weights = (area[indices-1]-area[indices+1])/2
    widths = 1/np.sqrt(model.pdf_afact/centers**2+1/variance[indices])
    density = np.zeros_like(gaps)
    for weight, center, width in zip(weights, centers, widths):
        density += weight * (np.exp(-.5*((gaps-center)/width)**2)
                             + np.exp(-.5*((gaps+center)/width)**2))/(np.sqrt(2*np.pi)*width)
    return {'description':'Diagnostic using full-precision Fortran u1 only; NOT production output.',
            'max_relative_density_error':float(max(abs(density/reference_density-1))),
            'source_u1_file_sha256':reference.sha256(directory/'out.1u1.2u1min')}


def plot(output,curves):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig,axes=plt.subplots(2,3,figsize=(14,8),layout='constrained')
    selected=sorted((item for item in curves if item[0]['young_modulus_pa']==3e6
                     and item[0]['delzeta']==.01),key=lambda item:item[0]['pressure_pa'])
    for column,(case,gaps,reference_density,calculated) in enumerate(selected):
        axes[0,column].loglog(gaps*1e6,reference_density,label='Fortran')
        axes[0,column].loglog(gaps*1e6,calculated,'--',label='Python (corrected)')
        axes[0,column].set(title=f"p = {case['pressure_pa']/1e6:g} MPa",ylabel='Nominal-area density (1/m)')
        axes[1,column].semilogx(gaps*1e6,100*(calculated/reference_density-1))
        axes[1,column].set(ylabel='Density relative difference (%)')
    for ax in axes.flat:ax.set_xlabel('Microgap (um)');ax.grid(alpha=.2)
    for ax in axes[0]:ax.legend()
    fig.savefig(output/'pdf_tail_corrected.png',dpi=180);plt.close(fig)


if __name__ == '__main__':
    main()
