"""Output shear factors and uniform-state pressure/motion flux decomposition."""
import csv
import json
from pathlib import Path

from .shear import shear_length


def run_shear_flow(law,pressure_rows,output,*,plot=False):
    output=Path(output);rows=[]
    sign=1. if law.options.rough_surface=='upper' else -1.
    for r in pressure_rows:
        h=r['mean_gap_m'];length=float(shear_length(h,law.rms_height_m,law.percolation_gap_m))
        phi=length/law.rms_height_m if law.rms_height_m>0 else 0.
        c0,c1=law.coefficients(h)
        motion=float(law.source(h))
        pressure=-r['hydraulic_conductivity_m3_pa_s']*law.options.pressure_gradient_pa_m
        rows.append({'contact_pressure_pa':r['pressure_pa'],'mean_gap_m':h,
            'pressure_flow_factor':r['pressure_flow_factor'],
            'shear_flow_factor':phi,'signed_shear_flow_factor':sign*phi,
            'shear_length_m':length,'lower_velocity_coefficient_m':float(c0),
            'upper_velocity_coefficient_m':float(c1),
            'pressure_flux_per_width_m2_s':pressure,
            'motion_flux_per_width_m2_s':motion,
            'combined_fixed_state_flux_per_width_m2_s':pressure+motion,
            'pressure_conductivity_status':r['status']})
    path=output/'shear_flow_table.csv'
    with path.open('w',encoding='utf-8',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=list(rows[0]),lineterminator='\n');writer.writeheader();writer.writerows(rows)
    summary={**law.metadata,'table_file':path.name,'results':rows,
        'table_flux_scope':'Each row is a uniform fixed contact state at the given test pressure gradient; not the spatial seal solution.',
        'blocked_state_scope':'A nonzero motion source in a pressure-blocked state is not evidence of through-seal leakage; spatial moving flow requires positive mobility.',
        'uncertainty_scope':'Shear interpolation does not supply PDF-derived confidence intervals. Pressure PDF bounds do not bound the combined moving-wall flux.'}
    (output/'shear_diagnostics.json').write_text(json.dumps(summary,indent=2,allow_nan=False)+'\n')
    if plot:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        p=[r['contact_pressure_pa']/1e6 for r in rows]
        fig,axes=plt.subplots(1,2,figsize=(10,4),layout='constrained')
        axes[0].semilogx(p,[r['shear_flow_factor'] for r in rows])
        axes[0].set(xlabel='Contact pressure (MPa)',ylabel='Shear factor amplitude',title='Single rough surface interpolation')
        for key,label in [('pressure_flux_per_width_m2_s','Pressure'),('motion_flux_per_width_m2_s','Motion'),
                          ('combined_fixed_state_flux_per_width_m2_s','Combined')]:
            axes[1].semilogx(p,[r[key] for r in rows],label=label)
        axes[1].set(xlabel='Contact pressure (MPa)',ylabel='Fixed-state flux per width (m^2/s)')
        axes[1].legend()
        for ax in axes:ax.grid(alpha=.25)
        fig.savefig(output/'shear_flow.png',dpi=180);plt.close(fig)
    return summary
