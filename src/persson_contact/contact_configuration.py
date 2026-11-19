"""Single-side roughness selection with explicit, inactive future interfaces."""
from copy import deepcopy
from pathlib import Path


def resolve_configuration(parameters):
    p=deepcopy(parameters)
    rough=p.get('mixed_friction',{}).get('rough_surface','upper')
    if rough=='both':
        raise NotImplementedError('rough_surface=both 是双侧粗糙预留接口，当前尚未实现；不会自动合并两份PSD')
    if rough not in ('upper','lower'):raise ValueError('rough_surface must be upper/lower (both reserved)')
    raw=p.get('roughness_parameters')
    if raw is not None:
        if not isinstance(raw,dict) or set(raw)-{'upper','lower'}:
            raise ValueError('roughness_parameters requires upper/lower blocks')
        active=raw.get(rough,{})
        if not isinstance(active,dict) or set(active)-{'path','scale_c','tail_is_extrapolated'}:
            raise ValueError('Invalid active roughness_parameters block')
        if active.get('path') is not None:
            legacy=p.get('psd',{})
            if legacy.get('path') is not None:
                old={'scale_c':1.,'tail_is_extrapolated':False,**legacy}
                new={'scale_c':1.,'tail_is_extrapolated':False,**active}
                if old!=new:raise ValueError('psd 与选中的 roughness_parameters 冲突；使用新接口时将 psd.path 设为 None')
            p['psd']={'scale_c':1.,'tail_is_extrapolated':False,**active}
    macro=p.get('macro_contact',{})
    body=macro.get('deformable_body')
    if body not in (None,'upper','lower'):raise ValueError('deformable_body must be upper/lower')
    if body=='lower':
        raise NotImplementedError('deformable_body=lower 已预留；真实下壁形变、受力投影及黏弹运动学将在下一阶段实现。当前可使用既有等效弹性表示')
    if body=='upper' and macro.get('body_representation','upper_elastic_lower_rigid')!='upper_elastic_lower_rigid':
        raise ValueError('deformable_body 与旧 body_representation 冲突；只保留一个明确的构型')
    if p.get('material',{}).get('model')=='viscoelastic':
        raise NotImplementedError('viscoelastic 为预留材料接口，当前稳定版尚未实现主曲线读取及黏弹性计算，不会回退为弹性')
    return p
