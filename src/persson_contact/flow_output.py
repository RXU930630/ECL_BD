"""Flow-table outputs for the Python workflow; independent physics in flow.py."""
import csv
from dataclasses import asdict
import json
from pathlib import Path

from .applications import run_leakage
from .flow import pressure_flow


def _csv(path, rows):
    with Path(path).open('w',encoding='utf-8',newline='') as stream:
        writer=csv.DictWriter(stream,fieldnames=list(rows[0]))
        writer.writeheader();writer.writerows(rows)


def run_pressure_flow(model, pressures, interface, options, output, *, percolation,
                      leakage=None, plot=False):
    output=Path(output)
    results=[pressure_flow(model,float(p),interface,options,
                          percolation_relative_contact_area=percolation) for p in pressures]
    rows=[value.as_dict() for value in results]
    _csv(output/'pressure_flow_table.csv',rows)
    leaks=[run_leakage(value,leakage,interface) for value in results] if leakage is not None else []
    if leaks:
        _csv(output/'leakage_table.csv',leaks)
    summary={'scope':'Isotropic pressure flow; incompressible Newtonian no-slip fluid; fixed dry contact.',
        'interface':{'rheology':'newtonian','wall_slip':'no_slip','viscosity_pa_s':interface.viscosity_pa_s},
        'options':asdict(options),'percolation_relative_contact_area':percolation,
        'convention':'K [m^3]; hydraulic sigma=K/(12 eta) [m^3/(Pa s)]; phi_p=K/h^3.',
        'missing_pdf_policy':'Unrepresented PDF mass is zero hydraulic conductance in the default/lower result, '
          'not additional real contact. The upper result assigns it infinite conductance while true contact remains zero.',
        'reference_scope':'Sealing3 isotropic Bruggeman branch and area-dependent dimension; '
          'direct PDF integration only. Large-gap asymptotic switching, continuity correction, '
          'anisotropy, boundary films and shear-flow factors are not implemented.',
        'leakage_parameters':asdict(leakage) if leakage is not None else None,
        'leakage_scope':'Uniform fixed contact; no fluid load feedback, cavitation or cylinder mechanics.' if leaks else None,
        'status_counts':{name:sum(row['status']==name for row in rows) for name in sorted({row['status'] for row in rows})},
        'results':rows,'leakage_results':leaks}
    (output/'flow_diagnostics.json').write_text(json.dumps(summary,indent=2,allow_nan=False)+'\n')
    if plot:
        _plot(output,results,leaks)
    return summary


def _plot(output,results,leaks):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import numpy as np
    p=np.array([r.pressure_pa for r in results])/1e6
    factor=np.array([r.pressure_flow_factor for r in results])
    upper=np.array([r.upper_pressure_flow_factor if r.upper_pressure_flow_factor is not None else np.nan for r in results])
    fig,axes=plt.subplots(1,3,figsize=(15,4),layout='constrained')
    axes[0].semilogx(p,factor,label='Zero-conductance missing phase')
    axes[0].semilogx(p,upper,'--',label='Unknown-PDF upper bound')
    axes[0].set(ylabel='Pressure flow factor',title='Isotropic pressure flow')
    axes[0].legend(fontsize=8)
    conductance=np.array([r.geometric_conductivity_m3 for r in results])
    axes[1].loglog(p,np.where(conductance>0,conductance,np.nan))
    axes[1].set(ylabel='Geometric conductivity (m^3)',title='Positive conductance (blocked = 0)')
    if leaks:
        axes[2].semilogx(p,[r['flow_rate_m3_s'] for r in leaks])
        axes[2].set(ylabel='Signed flow rate (m^3/s)',title='Fixed-contact leakage')
    else:
        axes[2].semilogx(p,[r.unrepresented_noncontact_probability for r in results])
        axes[2].set(ylabel='Nominal-area probability',title='Unrepresented noncontact probability')
    for ax in axes:ax.set_xlabel('Nominal dry contact pressure (MPa)');ax.grid(alpha=.25)
    fig.savefig(output/'pressure_flow.png',dpi=180);plt.close(fig)
