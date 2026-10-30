# 参数文件

日常使用只需编辑 `main.py` 和 `parameters.py`。main选择计算任务；parameters选择工况、模型与数值精度。开发复算输入单独保存在tests/data中，不需要在日常入口中切换。

| 位置 | 内容 |
| --- | --- |
| `CASE` | 材料、PSD、流体、几何、压力、接触卸载与渗流规则 |
| `NUMERICS` | 接触数据库节点、模型修正参数、积分精度、插值与空间求解器 |
| `OUTPUT_DIRECTORY` | 输出目录；新工况使用新目录 |
| `PARAMETERS` | 自动合并后的程序输入，通常不用编辑 |

## 液体与气体

在CASE的leakage分组中设置：

```python
'fluid': 'liquid',             # 或'gas'
'compute_equivalent_gas': True,
```

liquid使用三次导流律，compute_equivalent_gas控制附加气体估算。gas独立求解二次导流和ballistic气体压力，估算开关不生效。gas目前仅表示自由分子ballistic模型，尚无可压缩黏性气体或过渡流求解器。

液体黏度在 `CASE['interface']['rheology_parameters']['newtonian']`，气体性质在 `CASE['leakage']['gas']`。气体参考状态、单位与模型范围见 [气体说明](gas_leakage.md)。空间计算需要在main.py开启RUN_LEAKAGE。

## 几何选择

```python
'geometry': 'oring_hertz',
```

只修改这一行选择几何。geometry_parameters中完整列出所有备选参数，程序只读取选中的一块。

| 名称 | 模型 | 所需参数 |
| --- | --- | --- |
| oring_hertz | O-ring的Hertz接触压力 | 峰值接触压力，以及截面直径或接触长度（二选一） |
| square_ring | 均匀接触压力 | 流向接触长度、参考接触压力 |
| prescribed_profile | 导入接触压力曲线 | 压力文件、文件格式、插值方式 |
| gaussian | 高斯接触压力 | 流向计算域长度、峰值压力、高斯分布标准差 |

Hertz截面直径不是密封环整体直径。Square ring当前不从截面尺寸计算变形。压力文件使用CSV表头x_m,contact_pressure_pa，单位为m、Pa；保留原坐标，不自动补零。入口位于小坐标端。

所有几何共用seal_width_m、入口和出口压力、接触卸载规则及渗流规则。环形密封的seal_width_m填2πR。spc平台和physical_zero代表不同导流处理，仍保留为正常计算选项。

## 牛顿流体与壁面滑移

interface中用rheology选择流变，用wall_slip选择壁面边界；所有已支持预设均列出，仅读取所选块。目前流变为newtonian恒黏度；壁面为no_slip或navier，后者分别给lower_slip_length_m、upper_slip_length_m（单位m，允许0）。两侧均为0恢复无滑移。非零滑移时需将leakage.compute_equivalent_gas设为False，避免套用无滑移气体估算。

旧的interface.viscosity_pa_s仍可读取，但不能与新的rheology_parameters同时填写。压温相关黏度和非牛顿流体尚未启用，只有底层接口预留。详见[模块3说明](module3_interface.md)。

## 剪切流动与粗糙侧

main.py中的COMPUTE_SHEAR_FLOW默认关闭。打开后，在CASE的shear_flow块选择rough_surface=lower/upper、两侧速度和测试压力梯度，所有字段均有注释。可以关闭RUN_LEAKAGE只生成因子表；非零速度与空间输运一起计算时，当前需fixed_contact、无滑移、有效介质，并关闭附加气体估算。详见[剪切流动模块](shear_flow.md)。

## PDF积分与精度

NUMERICS.options中的pdf_scale_method选择尺度积分：conservative使用准确谱积分和面积差权重，当前日常参数推荐这一项；fortran保留历史右矩形积分供复算。旧JSON或直接API省略该字段时仍使用fortran。

pdf_log_step控制PDF尺度步长，n_q控制能量积分，pressure_flow.gap_points控制导流分箱，三者需分别检查收敛。导出的n_gap只控制显示网格。当前导流quadrature_relative_difference仅代表分箱收敛，不能证明PDF尺度或模型误差已消失。新工况可用步长0.01、0.005比较；近渗流还需检查遗漏概率及上下界。详细方法和数据见[独立数值审计](numerical_convergence/README.md)。
