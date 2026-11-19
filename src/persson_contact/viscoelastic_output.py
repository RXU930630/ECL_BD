"""Independent dry Persson viscoelastic friction. Never silently modifies mixed lubrication."""
import hashlib,json
from pathlib import Path
import numpy as np
from .material import parse_material
from .psd import read_psd
from .viscoelastic_material import MasterCurve
from .viscoelastic_friction import PerssonViscoelasticFriction
from .macro_output import _csv
from .experiment_excel import _workbook,_sheet


def run_viscoelastic_project(parameters,*,base_directory,save_plots=True):
    base=Path(base_directory).resolve();raw=parameters.get('viscoelastic',{})
    allowed={'reference_temperature_c','temperature_c','shift_factors_file','shift_factors','modulus_type',
        'frequency_policy','pressures_pa','speeds_m_s','q_intervals','angular_points','kappa','relative_tolerance'}
    if not isinstance(raw,dict) or set(raw)-allowed:raise ValueError('Unknown viscoelastic settings')
    required={'reference_temperature_c','temperature_c','pressures_pa','speeds_m_s'}
    if not required<=set(raw):raise ValueError('Missing viscoelastic settings: '+', '.join(sorted(required-set(raw))))
    material=parse_material(parameters['material'],base)
    if material.model!='viscoelastic':raise ValueError('Select material.model=viscoelastic for this task')
    body=parameters.get('macro_contact',{}).get('deformable_body')
    if body not in ('upper','lower'):raise ValueError('Explicit deformable_body upper/lower required for viscoelastic task')
    rough=parameters.get('mixed_friction',{}).get('rough_surface','upper')
    if rough==body:
        raise NotImplementedError('This kernel requires roughness on the rigid body. Soft-body roughness does not receive the same swept-asperity loss; macro viscoelastic loss is not implemented here')
    def path(value):
        if value is None:return None
        p=Path(value);return p if p.is_absolute() else base/p
    curve=MasterCurve(material.master_curve_file,material.poisson_ratio,
        reference_temperature_c=raw['reference_temperature_c'],temperature_c=raw['temperature_c'],
        shift_factors_file=path(raw.get('shift_factors_file')),shift_factors=raw.get('shift_factors'),
        modulus_type=raw.get('modulus_type','young'),frequency_policy=raw.get('frequency_policy','error'))
    psd_raw=parameters['psd'];psd_path=path(psd_raw['path']);psd=read_psd(psd_path,psd_raw.get('scale_c',1.))
    pressures=raw.get('pressures_pa');speeds=raw.get('speeds_m_s')
    for values in (pressures,speeds):
        if (not isinstance(values,(list,tuple)) or not values or any(isinstance(v,bool) or not np.isfinite(v) or v<=0 for v in values)):
            raise ValueError('Provide nonempty positive pressures_pa and speeds_m_s; static friction is excluded')
    strength=parameters.get('contact_friction',{})
    if strength.get('model','constant_shear_strength')!='constant_shear_strength':raise NotImplementedError('First dry viscoelastic task supports constant interface shear strength only')
    tau=strength.get('model_parameters',{}).get('constant_shear_strength',{}).get('shear_strength_pa')
    if isinstance(tau,bool) or tau is None or not np.isfinite(tau) or tau<0:raise ValueError('Provide nonnegative interface shear_strength_pa; use explicit zero to isolate hysteresis')
    qn=raw.get('q_intervals',512);an=raw.get('angular_points',64);kappa=raw.get('kappa',.5)
    tolerance=raw.get('relative_tolerance',.002)
    if isinstance(tolerance,bool) or not np.isfinite(tolerance) or not 0<tolerance<.1:raise ValueError('Invalid relative_tolerance')
    if qn>10000 or an>512:raise ValueError('Leave room for independent doubled q/angular convergence checks')
    model=PerssonViscoelasticFriction(psd,curve,q_intervals=qn,angular_points=an,kappa=kappa)
    name=parameters.get('output_directory')
    if not isinstance(name,str) or not name.strip():raise ValueError('Provide output_directory')
    out=path(name)
    if out.exists() and (not out.is_dir() or any(out.iterdir())):raise FileExistsError(f'输出目录非空：{out}')
    cases=[]
    for p in pressures:
        for v in speeds:
            cases.append(model.evaluate(p,v,tau));print(f'黏弹干摩擦：p={p:g} Pa, v={v:g} m/s, mu_visco={cases[-1]["mu_viscoelastic"]:.6g}',flush=True)
    checks=[]
    for label,probe in [('q',PerssonViscoelasticFriction(psd,curve,q_intervals=2*qn,angular_points=an,kappa=kappa)),
                        ('angle',PerssonViscoelasticFriction(psd,curve,q_intervals=qn,angular_points=2*an,kappa=kappa))]:
        for c in cases:
            fine=probe.evaluate(c['pressure_pa'],c['speed_m_s'],tau)
            for field in ('mu_viscoelastic','relative_contact_area'):
                difference=abs(fine[field]-c[field]);limit=1e-10+tolerance*max(abs(fine[field]),abs(c[field]))
                checks.append({'check':label,'pressure_pa':c['pressure_pa'],'speed_m_s':c['speed_m_s'],'quantity':field,
                    'baseline':c[field],'refined':fine[field],'absolute_difference':difference,'allowed_difference':limit,'passed':bool(difference<=limit)})
    summary={'workflow_status':'complete','scope':'Dry isothermal uniform-pressure rigid-rough/smooth-viscoelastic half-space; no macro geometry or liquid simulation',
        'numerical_status':'passed' if all(c['passed'] for c in checks) else 'needs_attention',
        'contact_configuration':{'rough_surface':rough,'deformable_body':body,'rigid_body':'lower' if body=='upper' else 'upper'},
        'material':curve.diagnostics,'psd_sha256':hashlib.sha256(psd_path.read_bytes()).hexdigest(),
        'q0_m_inv':psd.q0,'q1_m_inv':psd.q1,'q_intervals':qn,'angular_points':an,'kappa':kappa,
        'shear_strength_pa':float(tau),'parameters_snapshot':parameters,'cases':cases,'checks':checks,
        'limitations':['No macro viscoelastic deformation/rolling loss, fluid coupling, adhesion force, flash heating or nonlinear softening.',
            'mu_adhesive means interface shear tau*A/F, not a normal adhesive-force model.',
            'mu_fluid is zero because this task is dry. Mixed-lubrication and fitting tasks still reject viscoelastic materials.',
            'Finite-band interpolated DMA curves are not automatically causal; convergence is not experimental validation.'],
        'theory_reference':'https://doi.org/10.1063/5.0293616, equations 17-20'}
    out.mkdir(parents=True,exist_ok=True);_csv(out/'viscoelastic_cases.csv',cases);_csv(out/'viscoelastic_checks.csv',checks)
    (out/'viscoelastic_diagnostics.json').write_text(json.dumps(summary,ensure_ascii=False,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    book=_workbook();_sheet(book,'摩擦贡献',list(cases[0]),[[r[k] for k in cases[0]] for r in cases]);_sheet(book,'数值核查',list(checks[0]),[[r[k] for k in checks[0]] for r in checks]);_sheet(book,'计算说明',('项目','内容'),[(k,json.dumps(v,ensure_ascii=False) if isinstance(v,(dict,list)) else v) for k,v in summary.items() if k not in ('cases','checks')]);book.save(out/'viscoelastic_results.xlsx');book.close()
    if save_plots:
        import matplotlib.pyplot as plt
        fig,ax=plt.subplots()
        for p in pressures:
            rows=[r for r in cases if r['pressure_pa']==p];v=[r['speed_m_s'] for r in rows]
            for field,style in [('mu_total','-'),('mu_viscoelastic','--'),('mu_adhesive',':')]:ax.semilogx(v,[r[field] for r in rows],style,label=f'{p:g} Pa {field}')
        ax.set(xlabel='Sliding speed [m/s]',ylabel='Friction coefficient',title='Dry Persson viscoelastic friction');ax.legend();ax.grid(alpha=.2);fig.tight_layout();fig.savefig(out/'viscoelastic_friction.png',dpi=160);plt.close(fig)
    return summary
