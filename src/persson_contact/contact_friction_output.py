"""Dry contact-friction configuration, velocity sweep and scientific outputs."""
from dataclasses import dataclass
import csv
import json
from pathlib import Path

import numpy as np

from .contact_friction import finite,make_shear_strength,evaluate_contact_friction


@dataclass(frozen=True)
class ContactFrictionPlan:
    law: object
    relative_velocities_m_s: tuple
    temperature_k: float
    transverse_width_m: float | None


def prepare_contact_friction(raw,base_directory):
    allowed={'model','model_parameters','relative_velocities_m_s','temperature_k','transverse_width_m'}
    if not isinstance(raw,dict) or set(raw)-allowed:raise ValueError('Unknown contact_friction parameter')
    velocities=raw.get('relative_velocities_m_s')
    if not isinstance(velocities,(list,tuple)) or not velocities:
        raise ValueError('Provide nonempty relative_velocities_m_s')
    speeds=tuple(finite(v,'relative velocity') for v in velocities)
    temperature=finite(raw.get('temperature_k',298.15),'temperature_k',minimum=0,strict=True)
    width=raw.get('transverse_width_m')
    if width is not None:width=finite(width,'transverse_width_m',minimum=0,strict=True)
    return ContactFrictionPlan(make_shear_strength(raw,base_directory),speeds,temperature,width)


def calculate_friction_sweep(normal,plan):
    return [evaluate_contact_friction(normal.cell_edges_m,normal.relative_contact_area,normal.pressure_pa,
        normal.mean_gap_m,v,plan.law,temperature_k=plan.temperature_k,transverse_width_m=plan.transverse_width_m)
        for v in plan.relative_velocities_m_s]


def _csv(path,rows):
    with path.open('w',newline='',encoding='utf-8') as out:
        writer=csv.DictWriter(out,fieldnames=list(rows[0]),lineterminator='\n');writer.writeheader();writer.writerows(rows)


def write_contact_friction(out,plan,results,*,save_plots=True):
    out=Path(out);cases=[];profiles=[]
    for index,result in enumerate(results):
        name=f'contact_friction_profile_{index:03d}.csv';profiles.append(name)
        _csv(out/name,[{key:values[i] for key,values in result.profile.items()} for i in range(len(result.profile['x_m']))])
        cases.append({'case_index':index,**result.summary,'profile_file':name})
    _csv(out/'contact_friction_cases.csv',cases)
    diagnostics={'scope':'Solid-contact sliding friction postprocessed on the dry normal solution only.',
        'model':plan.law.metadata,'temperature_k':plan.temperature_k,'transverse_width_m':plan.transverse_width_m,
        'stress_convention':'tau_c is strength magnitude on real area; nominal traction magnitude=a*tau_c. p_real=p_c/a is the mean contact-spot pressure.',
        'direction_convention':'V=U_upper-U_lower along +x. Upper traction=-sign(V)*a*tau_c; lower traction is equal and opposite.',
        'zero_velocity_policy':'Applied sliding traction and dissipation are zero. Static friction and its limit are not solved; kinetic coefficient is blank.',
        'noncontact_policy':'No strength law is evaluated where a=0; strength/real-pressure profile placeholders are zero there.',
        'units':'Cell traction Pa; real area per width m; forces per width N/m; optional total force N; dissipation W/m or W.',
        'normal_state_reused':True,'velocity_affects_normal_contact':False,
        'limitations':['No fluid friction or viscoelastic hysteresis contribution.',
            'No tangential compliance, partial stick, rolling/spin, heating or friction/normal feedback.',
            'Constant/log-Gaussian laws have no built-in pressure or temperature dependence; custom callbacks receive both.',
            'Experimental strength values must be supplied; no strength is inferred from elastic modulus or PSD.'],
        'cases':cases,'outputs':['contact_friction_cases.csv','contact_friction_diagnostics.json',*profiles]}
    if save_plots:
        plot_contact_friction(out,results)
        diagnostics['outputs'].append('contact_friction.png')
    (out/'contact_friction_diagnostics.json').write_text(json.dumps(diagnostics,ensure_ascii=False,indent=2,allow_nan=False)+'\n')
    print(f'固体接触摩擦已计算：{len(cases)} 个相对速度工况；剪切模型 {plan.law.metadata["model"]}。')
    return diagnostics


def plot_contact_friction(out,results):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.ticker import ScalarFormatter
    moving=[r for r in results if r.summary['speed_m_s']>0]
    selected=moving[-1] if moving else results[-1]
    fig,axes=plt.subplots(1,3,figsize=(13,4),layout='constrained')
    x=selected.profile['x_m']*1e3
    axes[0].plot(x,selected.profile['nominal_sliding_strength_pa']/1e6,label='a * tau_c')
    axes[0].plot(x,selected.profile['real_contact_shear_strength_pa']/1e6,'--',label='tau_c (real contact)')
    axes[0].set(xlabel='x [mm]',ylabel='Strength magnitude [MPa]',title=f'Real vs nominal area: |V|={selected.summary["speed_m_s"]:g} m/s')
    axes[0].legend(fontsize=8)
    if moving:
        unique={r.summary['speed_m_s']:r for r in moving}
        speed=sorted(unique)
        axes[1].semilogx(speed,[unique[v].summary['kinetic_contact_friction_coefficient'] for v in speed],'o-')
        coefficients=[unique[v].summary['kinetic_contact_friction_coefficient'] for v in speed]
        coefficients=[v for v in coefficients if v is not None]
        if coefficients:axes[1].set_ylim(0,max(max(coefficients)*1.2,.01))
    axes[1].set(xlabel='Relative speed |V| [m/s]',ylabel='Contact kinetic coefficient',title='Dry solid-contact contribution')
    signed=sorted(results,key=lambda r:r.summary['relative_velocity_m_s'])
    axes[2].plot([r.summary['relative_velocity_m_s'] for r in signed],
        [r.summary['upper_friction_force_per_width_n_m'] for r in signed],'o')
    positive=[r.summary['speed_m_s'] for r in moving]
    if positive:
        axes[2].set_xscale('symlog',linthresh=min(positive)*.1)
        values=[r.summary['relative_velocity_m_s'] for r in results]
        ticks=sorted(set([min(values),0.,min(positive),max(values)]))
        axes[2].set_xticks(ticks,labels=[f'{v:g}' for v in ticks])
    axes[2].set(xlabel='Signed relative V [m/s]',ylabel='Upper contact force [N/m]',title='Signed force on upper body\nV=0: no static solve')
    for ax in axes:
        ax.grid(alpha=.25);ax.yaxis.set_major_formatter(ScalarFormatter(useOffset=False))
    fig.savefig(Path(out)/'contact_friction.png',dpi=170);plt.close(fig)
