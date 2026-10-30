# 500气体估算与5000独立ballistic计算

在 `main.py` 中开启 `RUN_LEAKAGE=True`，其余选择在所选参数文件的 `CASE['leakage']` 分组中完成（配置布局见 [参数说明](parameters.md)）：

```python
'fluid': 'liquid',                 # liquid对应500；gas对应5000的ballistic模型
'compute_equivalent_gas': True,    # 仅500生效；False时只输出液体结果
'gas': {
    'temperature_c': 25.0,
    'mean_molecular_speed_m_s': 471.0,
    'viscosity_pa_s': 1.82e-5,
    'ambient_pressure_pa': 1e5,
    'reference_pressure_pa': 1e5,
},
```

更换流体或工况时同时更换 `OUTPUT_DIRECTORY`，例如 `results/liquid` 与 `results/ballistic_gas`；已有非空目录会被拒绝覆盖。旧的 `fortran_mode=500/5000` 及 `transport='liquid'` 配置继续兼容，但新配置只用 `fluid` 选择流体，避免两个开关冲突。输出同时记录fluid与原Fortran编号。

所有入口、出口压力继续使用Pa表压。`ambient_pressure_pa` 给出表压零点对应的绝对压力；`reference_pressure_pa` 是报告气体体积时使用的参考绝对压力。两者分别设置，通常均为100000 Pa。体积参考温度为 `temperature_c`，不是额外指定的标准状态温度。

## 模式500

非零液体Navier滑移时需关闭compute_equivalent_gas；这里的辅助换算采用无滑移液体关系。独立气体模式不使用液体滑移参数。

求解现有不可压缩液体问题，使用局部导流律 `u^3/(12 eta_liquid)` 与原有局部卸载。`compute_equivalent_gas=True` 时，按原Fortran 500的辅助公式增加气体估算。

以总流量Q、密封流向长度L、横向宽度W和压差Δp定义：

```
K3_seal = 12 eta_liquid L Q / (W Δp)
u_eff = K3_seal^(1/3)
p_mean_abs = ambient_pressure + (p_in_gauge + p_out_gauge)/2
Q_diffusive_ref = Q_liquid * (eta_liquid/eta_gas) * p_mean_abs/p_ref
Q_ballistic_est_ref = Q_liquid * 12 eta_liquid vbar / (2 p_ref u_eff)
```

这两列使用液体解的接触状态与等效间隙，没有独立求解可压缩气体或ballistic压力分布。不能把估算列当作5000输出，也不默认将两列相加作为经过验证的过渡区模型。

零压差时两种估算均为零，等效间隙留空；阻断且压差非零时估算为零、等效间隙为零。压差反向时保留有符号流量。

## 模式5000

从同一接触模型与间隙PDF独立求解二次有效介质导流系数K2，不使用K3的幂换算：

```
local molecular mobility = vbar u^2 / (2 kB T)
K2 [m^2], phi_2 = K2/h^2
j_N = -vbar K2/(2 kB T) * dp/dx
d(j_N)/dx = 0
Q_ref = Ndot kB T / p_ref
```

空间压力与局部卸载独立重算；复用现有守恒Newton求解器，其未知通量在5000中是分子数通量。5000当前支持 `effective_medium`，沿用既有几何、`fixed_contact/local_unloading`、`spc`平台或`physical_zero`规则。`compute_equivalent_gas` 在5000中不生效。

该模式对应自由分子、壁面碰撞主导的ballistic极限，未加入黏性/Knudsen过渡插值或非局部FEM柔度。气体黏度用于500估算和平均自由程诊断，不参与5000的分子通量系数。平均热运动速度由用户填写；程序不会根据气体种类或温度自动替换它。

`physical_zero` 可能存在阻断区或需要可行初值的导通分支。确实阻断时输出零流量和未知区间压力；没有求得导通分支时明确失败，不以正导流下限伪造结果。底层 `solve_case(..., initial_pressure_profile=(x,p))` 可接收已收敛的相邻工况作为数值初值，最终仍按原始零导流规则求解。

## 输出字段

两种模式继续使用 `leakage_cases.csv`、`leakage_profile_*.csv` 和 `leakage_diagnostics.json`，由 `fortran_mode` 与 `flow_quantity` 标记物理量。

| 文件/字段 | 含义 |
| --- | --- |
| `pressure_flow_table.csv` | 500的K3表，单位m³ |
| `ballistic_flow_table.csv` | 5000的K2表，单位m²；另有分子迁移率 |
| `flow_rate_m3_s`、`flow_rate_mm3_h` | 总流量；500为液体体积，5000为参考压力下气体体积 |
| `square_unit_flow_m3_s` | 总体积流量乘L/W，可与原Fortran square-unit列对照；参考状态与主模式相同 |
| `estimated_diffusive_gas_reference_flow_*` | 500的扩散气体估算，含m³/s与mm³/h |
| `estimated_ballistic_gas_reference_flow_*` | 500的等效ballistic估算，含m³/s与mm³/h |
| `molecular_flow_rate_per_s` | 5000独立求得的总分子流率 |
| `gas_reference_volume_flow_*` | 5000参考压力下的总体积流率 |
| `gas_reference_pressure_pa`、`gas_temperature_k` | 明确标注体积参考状态 |
| `effective_ballistic_gap_m` | 5000空间等效K2的平方根，用于原程序等效间隙对照 |
| `knudsen_at_reference_pressure` | 参考压力下平均自由程/等效间隙；不是全空间Kn分布 |

环形密封应设置 `seal_width_m=2*pi*R`。例如R=9.525 mm时，宽度约0.05984734 m。总m³/s转总mm³/h只需乘 `3.6e12`；已有总体积流量不再除以L或乘周长。

PDF未解析概率的零/无限导流界限仍分别求解；它们不是统计置信区间，也不包含气体模型或结构近似误差。

## 常数与验证范围

Python使用精确SI玻尔兹曼常数 `1.380649e-23 J/K` 与 `T=temperature_c+273.15`。原程序部分位置使用 `1.380e-23` 与 `+273.0`。对比分子数时应区分这一常数差异；同一计算内换算参考体积时温度因子会抵消，平均热运动速度仍按给定值使用。

单元测试覆盖均匀狭缝、接触/非接触两相解析解、二次与三次均质化不可互换、反向/零压差、阻断、参考压力/温度/速度缩放、500估算对历史输出，以及两个主模式的完整工作流。原Fortran参考驱动保留原始`Sealing3gas`迭代和5000空间循环的字节；对照范围与误差记录在 `docs/gas_validation/report.md`。独立导流检查与共享导流表的空间检查分开记录，不宣称整包所有气体分支已验证。
