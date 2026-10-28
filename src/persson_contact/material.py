"""Mutually exclusive elastic constants and future viscoelastic master curves."""
from dataclasses import dataclass
import math
from pathlib import Path

from .psd import _positive_finite


@dataclass(frozen=True)
class Material:
    model: str
    effective_modulus_pa: float | None = None
    poisson_ratio: float | None = None
    master_curve_file: Path | None = None

    def elastic_modulus(self):
        if self.model != 'elastic':
            raise NotImplementedError(
                '黏弹性接口已预留，但主曲线读取及黏弹性接触计算尚未实现。'
                '当前计算请选择 material.model="elastic"；不会自动改用固定杨氏模量。')
        return self.effective_modulus_pa


def parse_material(raw, base_directory=Path('.')):
    """Validate a material specification without loading a future curve format.

    Elastic input accepts E* OR E/nu. Viscoelastic input requires a curve
    filename and nu, and forbids either fixed modulus. No curve is evaluated
    and no fixed-modulus fallback is provided in the current release.
    """
    allowed = {'model', 'effective_modulus_pa', 'young_modulus_pa',
               'poisson_ratio', 'master_curve_file'}
    if not isinstance(raw, dict) or set(raw) - allowed:
        raise ValueError('material 必须为参数字典，且只能包含已定义的材料字段。')
    model = raw.get('model', 'elastic')
    if model not in ('elastic', 'viscoelastic'):
        raise ValueError('material.model 必须为 elastic 或 viscoelastic。')
    values = {key: value for key, value in raw.items() if value is not None and key != 'model'}
    nu = values.get('poisson_ratio')
    if nu is not None:
        if isinstance(nu, bool):
            raise ValueError('poisson_ratio 必须满足 -1 < nu <= 0.5。')
        try:
            nu = float(nu)
        except (ValueError, TypeError) as exc:
            raise ValueError('poisson_ratio 必须为有限数值。') from exc
        if not math.isfinite(nu) or not -1 < nu <= .5:
            raise ValueError('poisson_ratio 必须满足 -1 < nu <= 0.5。')
    if model == 'viscoelastic':
        if 'young_modulus_pa' in values or 'effective_modulus_pa' in values:
            raise ValueError('使用 master curve 时不应同时给定固定杨氏模量或有效模量。')
        curve = values.get('master_curve_file')
        if not isinstance(curve, (str, Path)) or not str(curve).strip():
            raise ValueError('黏弹性材料必须给出 master_curve_file。')
        if nu is None:
            raise ValueError('主曲线不代替泊松比：请给出 poisson_ratio。')
        path = Path(curve)
        if not path.is_absolute():
            path = Path(base_directory) / path
        return Material(model, poisson_ratio=nu, master_curve_file=path.resolve())
    if 'master_curve_file' in values:
        raise ValueError('弹性模式请使用固定模量；主曲线请选择 viscoelastic 模式。')
    if 'effective_modulus_pa' in values:
        if set(values) != {'effective_modulus_pa'}:
            raise ValueError('Provide effective_modulus_pa OR young_modulus_pa and poisson_ratio, not both')
        return Material(model, _positive_finite(values['effective_modulus_pa'], 'effective_modulus_pa'))
    if set(values) != {'young_modulus_pa', 'poisson_ratio'}:
        raise ValueError('弹性材料需要 effective_modulus_pa，或 young_modulus_pa 与 poisson_ratio。')
    young = _positive_finite(values['young_modulus_pa'], 'young_modulus_pa')
    return Material(model, young / (1 - nu**2), nu)
