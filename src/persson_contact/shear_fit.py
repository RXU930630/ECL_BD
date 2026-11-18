"""Constant interface shear-strength inverse problem, independent of the solver.

mu = mu_fluid + mu_background + tau * (real area / TOTAL normal load).
No intercept, pressure adjustment, contact-area rescaling or hidden hysteresis.
"""
from dataclasses import asdict,dataclass
import csv
import math
from pathlib import Path

import numpy as np

from .contact_friction import finite


@dataclass(frozen=True)
class Experiment:
    row_number: int
    speed_m_s: float
    line_load_n_m: float
    friction_coefficient: float
    background_friction_coefficient: float = 0.
    standard_uncertainty_mu: float | None = None
    use_for_fit: bool = True

    def __post_init__(self):
        if type(self.row_number) is not int or self.row_number<1 or type(self.use_for_fit) is not bool:
            raise ValueError('Invalid experiment row number or use_for_fit flag')
        for key in ('speed_m_s','line_load_n_m'):
            object.__setattr__(self,key,finite(getattr(self,key),key,minimum=0,strict=True))
        for key in ('friction_coefficient','background_friction_coefficient'):
            object.__setattr__(self,key,finite(getattr(self,key),key))
        if self.standard_uncertainty_mu is not None:
            object.__setattr__(self,'standard_uncertainty_mu',finite(self.standard_uncertainty_mu,'standard_uncertainty_mu',minimum=0,strict=True))


def read_experiments(path,*,default_line_load_n_m,transverse_width_m=None):
    """Excel .xlsx or UTF-8/BOM CSV; same units and validation for both.

    Optional load columns: line_load_n_m OR normal_force_n (requires width).
    If both are given they must agree. Blank load uses the common macro load.
    Blank use_for_fit means True; false reserves a point for validation.
    """
    default=finite(default_line_load_n_m,'default line load',minimum=0,strict=True)
    width=None if transverse_width_m is None else finite(transverse_width_m,'width',minimum=0,strict=True)
    required={'speed_m_s','friction_coefficient'}
    allowed=required|{'line_load_n_m','normal_force_n','background_friction_coefficient','standard_uncertainty_mu','use_for_fit'}
    rows=[];path=Path(path)
    if path.suffix.lower()=='.xlsx':
        from .experiment_excel import read_excel_table
        header,table=read_excel_table(path)
    elif path.suffix.lower()=='.csv':
        with path.open(encoding='utf-8-sig',newline='') as stream:
            reader=csv.reader(stream);header=next(reader,[])
            table=[(number,row) for number,row in enumerate(reader,2) if row]
    else:raise ValueError('实验文件请选择.xlsx（Excel工作簿）或.csv；旧.xls文件请在Excel中另存为.xlsx。')
    header=['' if v is None else str(v).strip() for v in header]
    if (len(set(header))!=len(header) or not required<=set(header) or set(header)-allowed):
        raise ValueError('实验表必须包含speed_m_s,friction_coefficient；请按模板使用列名，不允许重复列名。')
    for number,cells in table:
        if len(cells)!=len(header):raise ValueError(f'实验表第{number}行列数不匹配')
        row={key:'' if value is None else str(value).strip() for key,value in zip(header,cells)}
        if not any(v.strip() for v in row.values()):continue
        def optional(key,*,positive=False):
            value=row.get(key,'').strip()
            return None if not value else finite(value,f'实验表第{number}行: {key}',minimum=0 if positive else None,strict=positive)
        try:
            speed=finite(row['speed_m_s'],f'实验表第{number}行: speed_m_s',minimum=0,strict=True)
            mu=finite(row['friction_coefficient'],f'实验表第{number}行: friction_coefficient')
            load=optional('line_load_n_m',positive=True);force=optional('normal_force_n',positive=True)
            if force is not None:
                if width is None:raise ValueError('normal_force_n需要几何轴向长度或contact_friction.transverse_width_m')
                converted=force/width
                if load is not None and not math.isclose(load,converted,rel_tol=1e-8):
                    raise ValueError('normal_force_n与line_load_n_m、横向宽度不一致')
                load=converted
            flag=row.get('use_for_fit','').strip().lower()
            if flag not in ('','true','false','1','0'):raise ValueError('use_for_fit只能为true/false或1/0')
            background=optional('background_friction_coefficient')
            uncertainty=optional('standard_uncertainty_mu',positive=True)
            rows.append(Experiment(number,speed,default if load is None else load,mu,
                0. if background is None else background,uncertainty,flag not in ('false','0')))
        except (ValueError,TypeError) as exc:raise ValueError(f'实验表第{number}行：{exc}') from exc
    if not rows:raise ValueError('实验表没有数据；请在模板中填写恒速段数据。')
    return rows


def fit_constant_strength(experiments,area_over_load_pa_inv,fluid_coefficients,*,
                          reference_strength_pa=1e6,weighting='equal',minimum_solid_coefficient_at_1mpa=1e-5,
                          speed_range_m_s=(None,None),diagnostics=None):
    """Nonnegative weighted least squares, one slope and no fitted intercept.

    Reference strength only scales the design variable; physical tau is
    invariant. Conditional statistical uncertainty excludes all model errors.
    At the nonnegative boundary, no symmetric standard error is reported.
    """
    reference=finite(reference_strength_pa,'reference_strength_pa',minimum=0,strict=True)
    threshold=finite(minimum_solid_coefficient_at_1mpa,'minimum_solid_coefficient_at_1mpa',minimum=0)
    if weighting not in ('equal','standard_uncertainty_mu'):raise ValueError('weighting must be equal or standard_uncertainty_mu')
    if not isinstance(speed_range_m_s,(list,tuple)) or len(speed_range_m_s)!=2:raise ValueError('speed_range_m_s must contain [minimum,maximum], either may be None')
    limits={'maximum_relative_rmse':.1,'maximum_tau_ratio_across_loads':2.,'maximum_reduced_chi_square':4.}
    diagnostics={} if diagnostics is None else diagnostics
    if not isinstance(diagnostics,dict) or set(diagnostics)-set(limits):raise ValueError('Unknown fit diagnostics settings')
    limits.update(diagnostics)
    for key,value in limits.items():
        if (isinstance(value,bool) or not isinstance(value,(int,float,np.integer,np.floating))
                or not np.isfinite(value) or value<=0):raise ValueError('Invalid fit diagnostic limit: '+key)
    lo,hi=(None if v is None else finite(v,'fit speed bound',minimum=0,strict=True) for v in speed_range_m_s)
    if lo is not None and hi is not None and lo>hi:raise ValueError('Invalid fit speed range')
    n=len(experiments);c=np.asarray(area_over_load_pa_inv,float);fluid=np.asarray(fluid_coefficients,float)
    if (n==0 or c.shape!=(n,) or fluid.shape!=(n,) or np.any(~np.isfinite(c)) or np.any(c<0)
            or np.any(~np.isfinite(fluid))):raise ValueError('Invalid real-area/load sensitivities or fluid coefficients')
    x=c*reference
    if np.any(~np.isfinite(x)):raise ValueError('Reference shear coefficient overflow')
    measured=np.array([e.friction_coefficient for e in experiments])
    background=np.array([e.background_friction_coefficient for e in experiments])
    y=measured-fluid-background
    reasons=[]
    for i,(e,xi) in enumerate(zip(experiments,x)):
        reason=''
        if not e.use_for_fit:reason='reserved_for_validation'
        elif (lo is not None and e.speed_m_s<lo) or (hi is not None and e.speed_m_s>hi):reason='outside_fit_speed_range'
        elif xi<=0 or c[i]*1e6<threshold:reason='insufficient_solid_contact_sensitivity'
        reasons.append(reason)
    selected=np.array([not v for v in reasons])
    if not np.any(selected):raise ValueError('没有可拟合的数据点：检查use_for_fit、速度范围及固体接触贡献。')
    w=np.ones(n)
    if weighting=='standard_uncertainty_mu':
        for i,e in enumerate(experiments):
            if e.standard_uncertainty_mu is None:
                if selected[i]:raise ValueError('加权拟合的每个参与点都需要standard_uncertainty_mu')
            else:w[i]=finite(e.standard_uncertainty_mu,'standard_uncertainty_mu',minimum=0,strict=True)
        # Normalize uncertainties before squaring to avoid overflow; absolute
        # sigma is retained separately for the conditional uncertainty.
        minimum=min(e.standard_uncertainty_mu for i,e in enumerate(experiments) if selected[i])
        w=(minimum/w)**2
    # Normalize the one-column design before dot products. This protects
    # valid tiny real-area sensitivities without introducing a contact floor.
    scale=float(np.max(x[selected]));z=x/scale
    denominator=float(np.sum(w[selected]*z[selected]**2))
    unconstrained=float(np.sum(w[selected]*z[selected]*y[selected])/denominator/scale)
    alpha=max(0.,unconstrained);tau=alpha*reference
    predicted=fluid+background+alpha*x;residual=measured-predicted
    if not all(np.all(np.isfinite(v)) for v in (tau,predicted,residual)):
        raise ValueError('Fit overflow; check experimental units and identifiability')
    count=int(np.sum(selected));dof=count-1;stderr=None
    if alpha>0:
        if weighting=='standard_uncertainty_mu':stderr=reference*minimum/(scale*math.sqrt(denominator))
        elif dof>0:stderr=reference*math.sqrt(float(np.sum(residual[selected]**2))/dof)/(scale*math.sqrt(denominator))
    warnings=[]
    if alpha==0:warnings.append('Nonnegative optimum is zero: supplied fluid/background may exceed experiment. This is not evidence of a precisely known zero interface strength.')
    if count<2:warnings.append('Only one training point; no residual-based estimate of repeatability or fit consistency.')
    if not np.any(~selected):warnings.append('No held-out points; transferability has not been checked.')
    warnings.append('Statistical standard error is conditional on exact area and fluid/background predictions; it excludes shear-cutoff, elasticity, PSD, cavitation and hysteresis model errors.')
    rmse=float(np.sqrt(np.mean(residual[selected]**2)))
    measured_rms=float(np.sqrt(np.mean(measured[selected]**2)))
    relative_rmse=None if measured_rms==0 else rmse/measured_rms
    issues=[];groups=[];seen=[]
    for load in sorted(set(e.line_load_n_m for e in experiments)):
        if any(np.isclose(load,v,rtol=1e-9,atol=0.) for v in seen):continue
        seen.append(load)
        mask=selected & np.array([np.isclose(e.line_load_n_m,load,rtol=1e-9,atol=0.) for e in experiments])
        count_group=int(np.sum(mask))
        if count_group:
            norm=float(max(x[mask]));design=x[mask]/norm
            with np.errstate(over='ignore',invalid='ignore'):
                group_tau=max(0.,float(np.sum(w[mask]*design*y[mask])/np.sum(w[mask]*design**2))/norm)*reference
            group_tau=group_tau if np.isfinite(group_tau) else None
            local_residual=None if group_tau is None else y[mask]-group_tau*c[mask]
            groups.append({'line_load_n_m':float(load),'training_points':count_group,
                'separately_fitted_strength_pa':group_tau,
                'global_fit_rmse_mu':float(np.sqrt(np.mean(residual[mask]**2))),
                'separate_fit_rmse_mu':None if local_residual is None else float(np.sqrt(np.mean(local_residual**2))),
                'mean_global_residual_mu':float(np.mean(residual[mask]))})
    eligible=[g['separately_fitted_strength_pa'] for g in groups if g['training_points']>=2 and g['separately_fitted_strength_pa'] is not None]
    tau_ratio=None if len(eligible)<2 or min(eligible)==0 else max(eligible)/min(eligible)
    if relative_rmse is not None and relative_rmse>limits['maximum_relative_rmse']:
        issues.append('整体残差偏大：当前常数tau与面积/流体预测不能充分解释实验；请勿直接将拟合值当作已验证界面属性。')
    if len(eligible)>=2 and ((min(eligible)==0<max(eligible)) or (tau_ratio is not None and tau_ratio>limits['maximum_tau_ratio_across_loads'])):
        issues.append('不同载荷分别拟合的tau差异偏大：需要检查工况、接触面积和常数剪切强度假设。')
    reduced_chi=None
    if weighting=='standard_uncertainty_mu' and dof>0:
        with np.errstate(over='ignore'):
            sigma=np.array([e.standard_uncertainty_mu for i,e in enumerate(experiments) if selected[i]])
            chi=float(np.sum((residual[selected]/sigma)**2)/dof)
        if np.isfinite(chi):reduced_chi=chi
        if not np.isfinite(chi) or chi>limits['maximum_reduced_chi_square']:
            issues.append('残差超过给定测量不确定度所能解释的水平；还需考虑模型或实验系统误差。')
    leverage=np.zeros(n);leverage[selected]=w[selected]*z[selected]**2/denominator
    if count>=3 and max(leverage)>.8:
        warnings.append('拟合主要由单个实验点决定（权重/接触灵敏度）；请检查tau的代表性。')
    if any(v=='insufficient_solid_contact_sensitivity' for v in reasons):
        warnings.append('部分点对固体剪切强度不敏感，未参与拟合；这些点的流体摩擦失配不能靠调整tau修正。')
    warnings.extend(issues)
    records=[]
    for i,e in enumerate(experiments):
        with np.errstate(over='ignore',invalid='ignore',divide='ignore'):
            implied=float(y[i]/c[i]) if selected[i] else None
        if implied is not None and not np.isfinite(implied):implied=None
        records.append({**asdict(e),'used_for_fit':bool(selected[i]),'exclusion_reason':reasons[i],
            'area_over_total_load_pa_inv':float(c[i]),'reference_solid_coefficient':float(x[i]),
            'tau_sensitivity_mu_per_mpa':float(c[i]*1e6),'fit_leverage':float(leverage[i]),
            'pointwise_implied_strength_pa':implied,
            'fluid_coefficient':float(fluid[i]),'fitted_solid_coefficient':float(alpha*x[i]),
            'predicted_total_coefficient':float(predicted[i]),'residual_experiment_minus_model':float(residual[i])})
    summary={'model':'constant_shear_strength','fitted_shear_strength_pa':tau,
        'reference_strength_pa':reference,'reference_multiplier':alpha,
        'unconstrained_shear_strength_pa':unconstrained*reference,'nonnegative_boundary_active':alpha==0,
        'conditional_standard_error_pa':stderr,'weighting':weighting,
        'training_points':count,'held_out_or_excluded_points':n-count,'residual_degrees_of_freedom':dof,
        'training_rmse_mu':rmse,'relative_training_rmse':relative_rmse,
        'fit_consistency_status':'needs_attention' if issues else 'no_diagnostic_threshold_exceeded',
        'diagnostic_limits':limits,'reduced_chi_square':reduced_chi,
        'per_load_diagnostics':groups,'separate_tau_max_min_ratio':tau_ratio,
        'interpretation':'Candidate conditional estimate, not a physically validated material property. Diagnostic thresholds are screening criteria, not proof of validity.',
        'held_out_or_excluded_rmse_mu':None if np.all(selected) else float(np.sqrt(np.mean(residual[~selected]**2))),
        'minimum_solid_coefficient_at_1mpa':threshold,'speed_range_m_s':[lo,hi],
        'force_convention':'All coefficients are signed resisting coefficients on the moving lower plane, divided by TOTAL imposed normal load.',
        'background_policy':'Only explicit CSV background_friction_coefficient is subtracted; missing values mean the user assumes zero. No viscoelastic or capillary contribution is calculated.',
        'warnings':warnings}
    return summary,records
