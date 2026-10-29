#!/usr/bin/env python3
"""Compare a completed benchmark matrix without modifying model parameters.

Keep native integrals and common-range integrals separate. Large relative
point-density errors in nearly zero tails are not used as agreement metrics.
"""
import argparse
from dataclasses import asdict
import json
from pathlib import Path

import numpy as np
from scipy.integrate import cumulative_trapezoid, trapezoid

import compare_fortran_reference as compare
from persson_contact.core import PerssonModel
from persson_contact.psd import read_psd

SCALE_FILE = 'out.1pressure.2zeta.3logzeta.4relArea.5separation.6hrms.zeta'
U1_FILE = 'out.1zeta.2du.over.dmu.3Area.4dA.over.dmu.5dlogA.over.dmu.6separation.7separation1'


def scalar_error(reference, python, floor=0.):
    reference, python = np.asarray(reference), np.asarray(python)
    mask = np.abs(reference) >= floor if floor else reference != 0
    return {'max_absolute': float(np.max(abs(reference-python))),
            'max_relative': float(np.max(abs((python[mask]-reference[mask])/reference[mask]))) if np.any(mask) else None,
            'relative_rows': int(np.sum(mask)), 'absolute_floor': floor}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--matrix', required=True)
    parser.add_argument('--psd', required=True)
    parser.add_argument('--output', required=True)
    parser.add_argument('--anchor')
    parser.add_argument('--exclude-case', action='append', default=[],
                        help='Explicitly exclude a named incomplete/cancelled run; record it in the report')
    args = parser.parse_args()
    root, output = Path(args.matrix), Path(args.output)
    all_manifests = sorted(root.glob('*/manifest.json'))
    excluded = [{'name': path.parent.name, 'status':json.loads(path.read_text()).get('status'),
                 'reason':json.loads(path.read_text()).get('cancellation_reason')}
                for path in all_manifests if path.parent.name in args.exclude_case]
    manifests = [path for path in all_manifests if path.parent.name not in args.exclude_case]
    if not manifests:
        parser.error('No case manifests found')
    for path in manifests:
        if json.loads(path.read_text()).get('status') != 'validated':
            parser.error(f'Case not complete: {path}; no comparison output created')
    output.mkdir(parents=True, exist_ok=False)
    psd = read_psd(args.psd)
    models, results = {}, []
    geometry_cache = {}
    curves, distributions, scale_curves = [], [], []
    for manifest_path in manifests:
        manifest = json.loads(manifest_path.read_text())
        if manifest.get('status') != 'validated':
            raise RuntimeError(f'Case not complete: {manifest_path}')
        case = manifest['case']
        e = case['effective_modulus_pa']
        model_key = (e, case['delzeta'])
        if model_key not in models:
            models[model_key] = PerssonModel(psd, e, .4, n_q=256, pdf_afact=1.,
                                               pdf_log_step=case['delzeta'])
            # Geometry-only quadrature depends on the PSD, not on E or p.
            models[model_key]._quadrature_cache = geometry_cache
        model = models[model_key]
        directory = manifest_path.parent
        options = argparse.Namespace(effective_modulus_pa=e, gamma=.4, pdf_afact=1., scale_c=1.)
        provenance = compare.completed_reference(directory, model, options)
        record = {'name': directory.name, **case, 'provenance': provenance}
        area = compare.read_numeric(directory/compare.AREA_NAME, 3, minimum_rows=1000)
        record['area_pressure_curve'] = scalar_error(area[:,1], model.area(area[:,0]))
        if case['mode'] == 2:
            data = compare.read_numeric(directory/compare.MEAN_NAME, 4, minimum_rows=108)
            h, a = model.mean_gap(data[:,0]), model.area(data[:,0])
            record['mean_gap_pressure_curve'] = scalar_error(data[:,1], h, 1e-12)
            record['relative_contact_area_pressure_curve'] = scalar_error(data[:,2], a)
            compare.csv(output, directory.name+'_pressure.csv',
                        ['pressure_pa','fortran_relative_contact_area','python_relative_contact_area',
                         'fortran_mean_gap_m','python_mean_gap_m'],
                        [data[:,0],data[:,2],a,data[:,1],h])
            curves.append((case, data, h, a))
        else:
            pressure = case['pressure_pa']
            scale = compare.read_numeric(directory/SCALE_FILE,6,minimum_rows=100)
            area_exact = model.area(pressure, psd.q0*scale[:,1])
            record['relative_contact_area_scale_curve'] = scalar_error(scale[:,3],area_exact)
            # Sample the full log-scale curve at 180 points and retain the
            # actual reference values/indices. No interpolated Fortran baseline.
            indices = np.unique(np.linspace(0,len(scale)-1,180).astype(int))
            selected = scale[indices]
            q = psd.q0*selected[:,1]
            python_a = model.area(pressure,q)
            python_h = np.array([model._mean_gap_above(pressure/ai,float(qi)) for qi,ai in zip(q,python_a)])
            record['mean_gap_scale_curve'] = scalar_error(selected[:,4],python_h,1e-12)
            compare.csv(output,directory.name+'_scale.csv',
                        ['magnification','fortran_relative_contact_area','python_relative_contact_area',
                         'fortran_ubar_at_scale_m','python_ubar_at_scale_m'],
                        [selected[:,1],selected[:,3],python_a,selected[:,4],python_h])
            scale_curves.append((case,selected,python_a,python_h))
            critical = model.critical_separation(pressure)
            u1_data = compare.read_numeric(directory/U1_FILE,7,minimum_rows=1)
            if critical.status == 'found' and np.min(scale[:,3]) <= .42:
                index = int(np.argmin(abs(u1_data[:,2]-.42)))
                source_u1 = float(u1_data[index,6])
                record['critical'] = {'fortran_nearest_relative_contact_area':float(u1_data[index,2]),
                                      'fortran_nearest_magnification':float(u1_data[index,0]),
                                      'fortran_uc_m':source_u1,'python':asdict(critical),
                                      'native_uc_relative_difference':float((critical.u1_m-source_u1)/source_u1)}
            else:
                record['critical'] = {'python':asdict(critical),'fortran_threshold_reached':bool(np.min(scale[:,3])<=.42)}
            _, pdf = compare.read_numeric(directory/compare.PDF_NAME,2,pressure_header=True,minimum_rows=40000)
            gaps,density = pdf.T
            python = model.gap_distribution(pressure,gaps)
            source_cdf = cumulative_trapezoid(density,gaps,initial=0.)
            python_cdf = python.continuous_cdf-python.continuous_cdf[0]
            # Integration by parts over the exact Gaussian CDF also captures
            # narrow components that point-density quadrature could miss.
            python_common_first = float(gaps[-1]*python_cdf[-1]-trapezoid(python_cdf,gaps))
            summary = compare.read_numeric(directory/compare.SUMMARY_NAME,5)[0]
            record['pdf'] = {
                'rows':len(pdf),'common_min_gap_m':float(gaps[0]),'common_max_gap_m':float(gaps[-1]),
                'common_cdf_max_absolute_difference':float(np.max(abs(source_cdf-python_cdf))),
                'common_density_L1_difference':float(trapezoid(abs(density-python.continuous_density_per_m),gaps)),
                'fortran_common_mass':float(source_cdf[-1]),'python_common_mass':float(python_cdf[-1]),
                'fortran_native_first_moment_m':float(summary[1]),
                'fortran_energy_mean_gap_m':float(summary[2]),
                'fortran_first_moment_relative_error':float(summary[1]/summary[2]-1),
                'python_native_first_moment_m':python.diagnostics.analytic_first_moment_m,
                'python_energy_mean_gap_m':python.mean_gap_m,
                'python_first_moment_relative_error':python.diagnostics.first_moment_relative_error,
                'common_first_moment_relative_difference':float(trapezoid(gaps*python.continuous_density_per_m,gaps)/trapezoid(gaps*density,gaps)-1),
                'native_first_moment_relative_difference':float(python.diagnostics.analytic_first_moment_m/summary[1]-1),
                'python_common_first_moment_from_cdf_m':python_common_first,
                'python_first_moment_outside_common_range_m':float(python.diagnostics.analytic_first_moment_m-python_common_first),
                'common_first_moment_cdf_relative_difference':float(python_common_first/summary[1]-1),
                'fortran_native_continuous_mass':float(1-summary[3]),
                'fortran_standard_relative_contact_area':float(summary[4]),
                'python_resolved_continuous_mass':python.diagnostics.resolved_continuous_mass,
                'python_unresolved_continuous_mass':python.diagnostics.unresolved_continuous_mass,
            }
            # Plot/export a thinned copy only after integrals over ALL 40000
            # points have been computed. Internal mass is never renormalized.
            take=np.unique(np.linspace(0,len(gaps)-1,1000).astype(int))
            compare.csv(output,directory.name+'_pdf.csv',
                        ['microgap_m','fortran_density_per_m','python_density_per_m','fortran_common_cdf','python_common_cdf'],
                        [gaps[take],density[take],python.continuous_density_per_m[take],source_cdf[take],python_cdf[take]])
            distributions.append((case,gaps[take],density[take],python.continuous_density_per_m[take],source_cdf[take],python_cdf[take]))
        results.append(record)
        print('Compared',directory.name,flush=True)
    report = {'case_count':len(results),'cases':results,'excluded_runs':excluded,
              'scope':'Supplied isotropic 2D PSD; nonadhesive elastic cases only; no leakage or mixed-friction validation.',
              'comparison_policy':'No fitted/rescaled curves. Native truncation and common-range metrics reported separately.',
              'scalar_relative_error_floor_m':1e-12,
              'batch_source_variant':'real-erf+dense-pdf-output; original contact equations unchanged'}
    if args.anchor:
        anchor=Path(args.anchor)
        manifest=json.loads((anchor/'manifest.json').read_text())
        report['original_anchor_status']=manifest['status']
        if manifest['status']=='validated':
            batch=root/'E3_p100k_mode3_d01'
            source_summary=compare.read_numeric(anchor/compare.SUMMARY_NAME,5)
            accelerated_summary=compare.read_numeric(batch/compare.SUMMARY_NAME,5)
            report['original_anchor_summary']=source_summary.tolist()
            report['original_anchor_source']=manifest['source']
            report['original_anchor_runtime_seconds']=manifest['run']['elapsed_seconds']
            report['original_anchor_provenance']={key:manifest[key] for key in
                ('compiler','compile_flags','binary','source','run','validation','used_input_hashes')}
            report['anchor_vs_accelerated_summary_max_absolute_difference']=float(np.max(abs(source_summary-accelerated_summary)))
            for name,columns in ((compare.AREA_NAME,3),(SCALE_FILE,6),(U1_FILE,7)):
                pristine=compare.read_numeric(anchor/name,columns)
                fast=compare.read_numeric(batch/name,columns)
                if pristine.shape != fast.shape:raise ValueError('Anchor output shapes differ')
                report.setdefault('anchor_output_differences',{})[name]=float(np.max(abs(pristine-fast)))
            _, sparse=compare.read_numeric(anchor/compare.PDF_NAME,2,pressure_header=True)
            _, dense=compare.read_numeric(batch/compare.PDF_NAME,2,pressure_header=True)
            selected=dense[np.r_[0,np.arange(99,len(dense),100)]]
            report['anchor_pdf_relative_density_difference']=scalar_error(sparse[:,1],selected[:,1])
            report['anchor_pdf_grid_relative_difference']=scalar_error(sparse[:,0],selected[:,0])
    (output/'comparison.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
    plot(output,curves,scale_curves,distributions)
    refinement_plot(output,curves,scale_curves)


def refinement_plot(output, curves, scales):
    import matplotlib.pyplot as plt
    fig,axes=plt.subplots(1,2,figsize=(12,4),layout='constrained')
    for case,data,h,a in curves:
        if case['young_modulus_pa']==3e6:
            axes[0].semilogx(data[:,0]/1e6,100*(h/data[:,1]-1),label=f"Fortran step {case['delzeta']:g}")
    for case,data,a,h in scales:
        if case['young_modulus_pa']==3e6 and case['pressure_pa']==1e5:
            keep=data[:,4]>=1e-12
            axes[1].semilogx(data[keep,1],100*(h[keep]/data[keep,4]-1),label=f"Fortran step {case['delzeta']:g}")
    axes[0].set(xlabel='Nominal pressure (MPa)',ylabel='Python / Fortran - 1 (%)',title='Mean-gap pressure curve')
    axes[1].set(xlabel='Magnification',ylabel='Python / Fortran - 1 (%)',title='Scale mean gap (reference >= 1e-12 m)')
    for ax in axes:ax.grid(alpha=.2);ax.legend()
    fig.savefig(output/'refinement.png',dpi=180);plt.close(fig)


def plot(output,curves,scales,pdfs):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    primary=[item for item in curves if item[0]['delzeta']==.01]
    fig,axes=plt.subplots(2,3,figsize=(14,8),layout='constrained')
    for column,(case,data,h,a) in enumerate(primary):
        axes[0,column].semilogx(data[:,0]/1e6,data[:,2],'-',label='Fortran')
        axes[0,column].semilogx(data[:,0]/1e6,a,'--',label='Python')
        axes[1,column].loglog(data[:,0]/1e6,data[:,1]*1e6,'-',label='Fortran')
        axes[1,column].loglog(data[:,0]/1e6,h*1e6,'--',label='Python')
        axes[0,column].set(title=f'E = {case["young_modulus_pa"]/1e6:g} MPa',ylabel='Relative contact area A/A0')
        axes[1,column].set(ylabel='Mean separation (um)')
    for ax in axes.flat:ax.set_xlabel('Nominal pressure (MPa)');ax.grid(alpha=.2);ax.legend()
    fig.savefig(output/'pressure_curves.png',dpi=180);plt.close(fig)
    primary_scales=[s for s in scales if s[0]['delzeta']==.01 and s[0]['young_modulus_pa']==3e6]
    fig,axes=plt.subplots(2,3,figsize=(14,8),layout='constrained')
    for column,(case,data,a,h) in enumerate(sorted(primary_scales,key=lambda x:x[0]['pressure_pa'])):
        axes[0,column].semilogx(data[:,1],data[:,3],'-',label='Fortran')
        axes[0,column].semilogx(data[:,1],a,'--',label='Python')
        axes[0,column].axhline(.42,color='grey',linestyle=':')
        axes[0,column].set(title=f'p = {case["pressure_pa"]/1e6:g} MPa',ylabel='Relative contact area A/A0')
        axes[1,column].loglog(data[:,1],np.where(data[:,4]>0,data[:,4],np.nan)*1e6,'-',label='Fortran')
        axes[1,column].loglog(data[:,1],np.where(h>0,h,np.nan)*1e6,'--',label='Python')
        axes[1,column].set(ylabel='Apparent-contact mean separation (um)',ylim=(1e-6,10.))
    for ax in axes.flat:ax.set_xlabel('Magnification');ax.grid(alpha=.2);ax.legend()
    fig.savefig(output/'scale_curves.png',dpi=180);plt.close(fig)
    primary_pdf=[p for p in pdfs if p[0]['delzeta']==.01 and p[0]['young_modulus_pa']==3e6]
    fig,axes=plt.subplots(2,3,figsize=(14,8),layout='constrained')
    for column,(case,u,f,p,fc,pc) in enumerate(sorted(primary_pdf,key=lambda x:x[0]['pressure_pa'])):
        axes[0,column].loglog(u*1e6,np.maximum(f,1e-30),'-',label='Fortran')
        axes[0,column].loglog(u*1e6,np.maximum(p,1e-30),'--',label='Python')
        axes[0,column].set(title=f'p = {case["pressure_pa"]/1e6:g} MPa',ylabel='Nominal-area gap density (1/m)')
        axes[1,column].semilogx(u*1e6,fc,'-',label='Fortran')
        axes[1,column].semilogx(u*1e6,pc,'--',label='Python')
        axes[1,column].set(ylabel='Cumulative probability on common range')
    for ax in axes.flat:ax.set_xlabel('Microgap (um)');ax.grid(alpha=.2);ax.legend()
    fig.savefig(output/'pdf_curves.png',dpi=180);plt.close(fig)


if __name__=='__main__':main()
