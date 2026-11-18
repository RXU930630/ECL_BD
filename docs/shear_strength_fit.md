# 独立任务：由实验摩擦反演常数界面剪切强度

当前任务适用于恒速稳定段、持续供油的一维等温弹性圆柱—平面混合润滑模型。硅油使用恒黏度牛顿模型、无滑移、守恒空化，要求连通导流分支。实验启动／停止、有限供油、热效应、黏弹滞后和毛细作用尚未计算。

## 方程与参考倍率

```text
A_real / F_N = (integral a(x) dx) / Wprime     [Pa^-1]
mu_model = mu_fluid + mu_background + tau * (A_real / F_N)
mu_solid_ref = tau_ref * (A_real / F_N)
mu_model = mu_fluid + mu_background + alpha * mu_solid_ref
tau_fit = alpha * tau_ref, alpha >= 0
```

F_N和Wprime均为总外加法向载荷，绝不只用固体接触承担的载荷。当前法向解不依赖tau；求出面积和流体力后，用一列设计变量作非负最小二乘。多点联合拟合一个常数，无额外自由截距。改变参考值不改变物理拟合结果；其默认1MPa只是计算尺度。

流体摩擦保留其方向，背景项按用户给定值相加，不随alpha缩放。viscoelastic friction是橡胶滞后，fluid friction是硅油黏性摩擦，两者不同。未提供背景项时表示用户暂时假定它为0，程序不会生成或推断橡胶滞后、毛细力或装置阻力。

主系数针对运动下平面的阻力。若实验测的是圆柱受力，需要在稳态、惯性可忽略且端部／外部供油压力作用已一致处理的条件下对应其反作用力；起停惯性或未扣除的端部压力力不能直接拟合成界面剪切强度。

## 实验Excel（推荐）或CSV

推荐用Excel打开`inputs/experiment_template.xlsx`，填写后另存为`inputs/experiment.xlsx`；参数`experiment_file`默认指向此文件。读取名为`experiment`的工作表；没有这个名称时读取第一个工作表。所有载荷放在同一个数据工作表，每个速度/载荷组合一行。首行保留英文表头，中文解释在批注和`填写说明`页。

Excel可以直接填写数字和TRUE/FALSE。公式使用Excel保存的计算值，程序不重算公式；没有保存的计算值或出现Excel错误值时，会提示单元格位置。旧`.xls`需在Excel中另存为`.xlsx`，不能只改扩展名。

也兼容原逗号分隔UTF-8 CSV，模板在`inputs/experiment_template.csv`。两种格式最少需要下面两列：

```csv
speed_m_s,friction_coefficient
```

填写数据后再运行。字段如下：

| 字段 | 必填 | 单位／定义 |
| --- | --- | --- |
| speed_m_s | 是 | 正的速度大小m/s；不接受0，当前不拟合静摩擦 |
| friction_coefficient | 是 | 以总法向载荷归一化的阻力系数；按阻碍运动为正的统一约定 |
| line_load_n_m | 否 | 总线载荷N/m；可逐行变化 |
| normal_force_n | 否 | 总法向力N；有限曲面模式自动使用轴向长度B，旧圆柱模式填写contact_friction.transverse_width_m；Wprime=F_N/B |
| background_friction_coefficient | 否 | 独立已知背景项；空白取0，正负按同一阻力约定 |
| standard_uncertainty_mu | 否 | 摩擦系数的正的标准不确定度；加权拟合时参与点必须填写 |
| use_for_fit | 否 | true/false或1/0；空白为true，false作留出验证 |

两个载荷列均空时使用macro_contact的默认载荷（总力或线载荷）。两个载荷列同时给出时必须与宽度换算一致。有限曲面的输入和边界解释见[几何说明](finite_geometry.md)。摩擦系数不是摩擦力；不要把N、mm/s、MPa数值直接填入其他单位的列。

两种格式均支持重复实验点和空行；拒绝NaN、Infinity、未知／重复列名和列数不匹配。CSV支持UTF-8 BOM。重复实验使用同一个正向状态，但保留为多个测量点。滑动方向在shear_strength_fit.sliding_direction中选择+1或-1，表示圆柱随动坐标系下平面的运动方向，实验表始终填写正速度大小。

## 参数与来源

开启RUN_SHEAR_STRENGTH_FIT时关闭独立泄漏、干接触、干摩擦和混合润滑任务；内部自动依赖所需模块。共用物理参数仍在原分组，不再复制一套“拟合材料”。

| 参数 | 定义 |
| --- | --- |
| model | 当前只拟合constant_shear_strength |
| experiment_file | 相对于项目根目录的.xlsx或.csv实验文件路径 |
| source | simulate或saved_forward |
| source_parameters.simulate | 当前无额外参数，使用共用物理输入 |
| source_parameters.saved_forward.directories | 已完成混合润滑结果目录列表，每个目录含mixed_cases.csv与mixed_diagnostics.json |
| reference_strength_pa | 参考倍率的强度尺度，默认1e6Pa |
| sliding_direction | 下平面方向+1或-1 |
| weighting | equal等权；standard_uncertainty_mu使用测量标准不确定度的倒方差 |
| speed_range_m_s | [最小,最大]；None表示不设该边界 |
| minimum_solid_coefficient_at_1mpa | 接触敏感性筛选；固定以1MPa计算固体系数，与选用的tau_ref无关 |

simulate按实验中所有速度／载荷计算，包括留出点。每个不同载荷生成一次正向速度扫描，参考强度临时替换为常数tau_ref，不修改原参数文件。各载荷首次创建输运数据库可能需要数分钟，具体取决于机器和积分精度。

saved_forward不读取当前PSD或重新计算材料，保存的物理工况是唯一依据。多个目录必须有相同材料、PSD哈希、流体、几何和闭合设置，允许不同线载荷和速度。按速度、方向和线载荷精确匹配，不对不同速度之间的面积／摩擦自动插值。缺少所需状态或出现多个候选时明确报错。

输入缓存还会检查CSV与诊断是否一致、总载荷与质量守恒残差。旧结果仅保存了面积而没有流体力、或尚未完成，均不应拿来拟合。正向计算失败时记录forward_failed；已完成的load目录保留，未生成有效材料参数。

## 输出与解释

| 文件 | 内容 |
| --- | --- |
| fit_results.xlsx | Excel结果：拟合概况（含Pa/MPa）、各实验点、正向响应及所用实验 |
| fit_summary.json | tau（Pa）、参考倍率、参与／排除点数、RMSE、条件标准误差及诊断 |
| fit_cases.csv | 每个实验点的面积、流体项、拟合固体项、显式背景项、合计、残差和是否参与拟合 |
| forward_response.csv | 可复核的面积／总载荷灵敏度和流体系数，不是新的实验数据 |
| experiment_used.csv | 单位统一后的实验输入，可再次作为实验CSV读取 |
| fitted_contact_friction.json | 可供复制的常数剪切模型参数；不自动修改parameters.py |
| fit_diagnostics.json | 参数快照、实验哈希、缓存来源及其物理工况、限制和演示标记 |
| shear_strength_fit.png | 实验与总预测、固体／背景分项及残差图 |
| forward/load_*/ | simulate模式下完整的耦合计算结果 |

留下的点、超出速度范围的点及固体贡献过小的点均不参与求tau，但仍输出预测和排除原因。高流体承载区不适合独立确定固体剪切强度。如果所有参与点都缺乏接触敏感性，程序拒绝给出tau。

非负最优解落在0时会明确标记；这也可能表示流体／背景项已高于实验，并不证明真实剪切强度精确等于0。单点也能反演，但不能通过残差估计一致性。没有留出点时，不宣称参数已经具有可迁移性。

条件标准误差仅是给定面积、流体及背景模型后的统计量：等权时从参与点残差估计；已知测量不确定度加权时按其给定值传播。边界解不输出对称标准误差。它不包含PSD、模量、空化、边界膜、剪切下限等模型误差，不应解释成完整的材料参数置信范围。

fluid_friction.minimum_shear_gap_m默认1nm仅是演示取值。0.5/1/2nm等敏感性需要分别重算流体摩擦并拟合，不能只改变当前参数后复用旧结果。拟合量应先理解为该橡胶—玻璃—硅油体系、指定条件和当前模型下的有效界面剪切强度，再以未参与拟合的数据检验。
