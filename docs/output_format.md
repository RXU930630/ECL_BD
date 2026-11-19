# 固定发布入口与输出格式

从本次编码修复起，GitHub 的 `stable` 分支作为固定的纯净下载入口，下载和解压目录名不再按当次任务更换。版本记录在 `PACKAGE_MANIFEST.json` 与程序内部；“stable”表示固定的发布入口与文件接口，不代表物理模型已完成全部实验验证。

程序仍由 main.py 选择任务，parameters.py 指定工况及 OUTPUT_DIRECTORY。修复不会自动切换你的任务、几何、PSD、载荷、速度或输出目录。正向与拟合是不同任务，各自的文件名称和字段固定；今后的兼容修复不重命名已有文件、列或工作表。格式若确实必须不兼容变更，需要先说明并提供迁移方法。

## 固定结果

正向：`mixed_results.xlsx`、`mixed_cases.csv`、`mixed_boundaries.csv`、`mixed_contact_table.csv`、`mixed_transport_table.csv`、`mixed_diagnostics.json`、`mixed_profile_000.csv`等按速度编号的剖面。绘图开启时输出 `mixed_friction.png`、`mixed_profiles.png`、`mixed_geometry.png`。

数值核查：`mixed_numerical_checks.json`，实际执行比较时另有 `mixed_numerical_checks.csv`；辅助解固定在 `numerical_checks/mesh`、`tables_and_spectral_quadrature`、`domain`。原有 Excel 工作表保持“计算概况”“正向摩擦与面积”“数值核查差异”“核查执行情况”，其中差异页仅在实际比较时存在。

拟合：`fit_results.xlsx`、`fit_cases.csv`、`fit_summary.json`、`fit_diagnostics.json`、`fitted_contact_friction.json`、`forward_response.csv`、`experiment_used.csv`；绘图开启时输出 `shear_strength_fit.png`。正向参考解固定在 `forward/load_000`等按载荷编号的子目录。

拟合 Excel 工作表保持“拟合概况”“各载荷独立诊断”“各实验点”“正向响应”“所用实验”。字段、单位、正负号及计算公式与rc5相同。所有程序文本输出明确采用UTF-8；CSV继续保持原有UTF-8，Excel保持.xlsx，不通过删除中文或修改数据来规避Windows编码问题。

## 更新现有程序，保留自己的设置

先解压新的完整包。在**新包目录**打开终端，执行：

```text
python scripts/update_runtime.py "C:\你的旧项目目录"
```

目标路径必须是包含 main.py、parameters.py、src 的旧项目根目录。脚本先备份旧模块到旧项目的 `runtime_backup`，再仅更新运行模块。不会覆盖你的 main.py、parameters.py、输入文件、实验表和已有结果。备份目录已存在时会停止，避免覆盖备份。以后仍在原项目目录运行 `python main.py`。

## rc5编码错误后的已完成计算

截图中的 `needs_attention` 表示精度核查需要注意，不是拟合失败。最后的 `'charmap' codec can't encode...` 是写中文JSON时使用Windows默认编码导致的程序错误，现已改为明确UTF-8。

更新旧项目的运行模块后，可在新包目录运行：

```text
python scripts/recover_fit.py "C:\你的旧项目目录"
```

脚本读取旧项目现有参数，检查其记录确为失败状态，再复用完整正向解，在临时目录完成拟合。只有全部所需速度和载荷都能匹配且输出成功后，才将固定名称的拟合文件保存回**原输出目录**；旧失败文件备份在 `runtime_backup/fit_before_recovery`。原 main.py、parameters.py、正向解及输入文件不改写。已完成的拟合不会覆盖；恢复失败也不替换原结果。原精度核查失败或缺失时仍明确提示，不把恢复拟合当作通过精度核查。

如果希望手动复用：原输出目录的 `forward/load_000` 等目录中 `mixed_cases.csv` 完整且 `mixed_diagnostics.json` 的 `workflow_status` 为 `complete` 时，将 `shear_strength_fit.source` 设为 `saved_forward`，在对应 `directories` 中填写原正向目录；全部所需载荷、速度都必须有成功解。手动调用普通任务时仍需选择空输出目录，输出文件名和格式保持相同。

编码失败可能留下空的 `mixed_numerical_checks.json`，不能把它视为通过或已完成核查。复用旧解时没有完整数值核查记录，会明确标记未核查。若主 `mixed_diagnostics.json` 也不完整，则不能直接复用，需要更新后重算相应正向工况。

## 自动查表缓存（rc7）

混合润滑正向及 simulate 拟合自动使用项目根目录 `.cache/tables/`，无需开关。
`contact/` 保存接触表；`gap_scales/` 保存指定压力下构建间隙分布所需的
u1、概率权重和剩余高度方差；`transport/` 保存几何导流率和流动/摩擦因子。
完整 PDF 可由尺度数据重建，并非保存所有可能的微观间隙采样网格。

每个目录提供 `CONDITIONS.txt` 供直接查看条件摘要；`table.json` 包含完整条件和数据：PSD 文件名及内容 SHA256、
波数范围、高度 RMS、材料/有效模量、积分精度、压力或间隙网格、适用模型参数、
程序源码指纹、NumPy/SciPy 版本、生成时间和数据校验值。文件名仅作说明，不影响匹配。
黏度在当前牛顿无滑移模型中于求解时施加，缓存的是几何导流率；改变黏度可以复用它。
不同速度、载荷和宏观网格复用微观表；压力范围不足仍要求扩展 contact_table。
接触压力网格的分段点固定为最大表压力的 0.001 倍，不再随 Hertz 载荷移动。
这会使首次计算的离散点与旧版略有不同，应以数值核查评估离散误差。

只使用完整且校验通过的匹配数据，损坏后自动重算；写入采用临时文件原子替换。
计算失败前已完整生成的有效表可以继续复用。尺度严重未解析的数据不保存；
混合润滑会拒绝上界无界的异常输运点，提示减小 pdf_log_step，不能将其当作真实封堵。
程序或计算环境变化保守地使旧缓存失效；核查中的加密参数分别匹配缓存。
可删除整个 `.cache/` 强制重算，既有 results 输出不受影响。
本功能目前接入混合润滑/拟合任务，其他独立任务仍按原流程运行。

## cylinder_plane 自动计算域（rc8）

`macro_contact.domain_half_width_m=None` 或省略该字段时，按实际线载荷和有效模量计算
`a=sqrt(4 Wprime R/(pi E*))`，初始域半宽自动取 `L=6a`。拟合每个载荷分别计算；
总载荷仍需通过实际轴向长度换算成线载荷。显式正数仍优先，不会被程序覆盖。
终端报告 a、L；诊断 geometry 中报告 domain_selection、hertz_half_width_m 和实际区间。
自动域是初始化估计，不是远场收敛证明；混合润滑的扩大域核查使用实际 L 的两倍，
而不是重新使用 None 或自动通过。高速液膜、粗糙接触可能需要更大区间。
finite_cylindrical_pad 仍严格使用实际尺寸，不按 Hertz 估计扩展。

## Excel 简洁表格（rc9）

Excel 输入模板和结果工作簿统一使用普通表格：第一行列名，后续每行一条记录。
不再设置彩色表头、粗体、自动换行、特殊行高、冻结窗格、自动筛选或强制小数位数；
仅按表头长度设置列宽。数字保持数值及 General 格式，不预填格式化空行。
文件名、工作表名、列名、列顺序及输入校验保持原有规则，表头批注和填写说明仍保留。

mixed_cases.csv / mixed_results.xlsx「正向摩擦与面积」：
mu solid = lower_resisting_solid_coefficient；mu fluid = lower_resisting_fluid_coefficient；
mu total = lower_resisting_friction_coefficient。均以总法向载荷归一化，fluid 包含滑动和压力流两项。
带 force_n_m 的列是单位轴向长度的有符号切向力，不是摩擦系数。

fit_cases.csv / fit_results.xlsx「各实验点」：
mu solid = fitted_solid_coefficient；mu fluid = fluid_coefficient；
mu total = predicted_total_coefficient（包含表格明确提供的 background_friction_coefficient）。
若从拟合任务的 forward 目录读取 mixed_cases.csv，则其中固体摩擦按参考 tau 计算，
而不是最后拟合的 tau；应使用 fit_cases.csv 中的上述列取得最终拟合结果。

## 多载荷正向计算及构型预留（rc10）

RUN_MIXED_FRICTION 正向任务中，macro_contact.line_load_n_m 既接受单个N/m数值，
也接受非空、正值、无重复的列表，例如 [310.,910.,1510.,2110.]。或者用
normal_force_n=[31.,91.,151.,211.] 输入总力N，再按轴向长度逐项换算一次。
两种载荷字段不能同时给出。输入顺序即运行顺序；先校验所有工况再开始。
单值输入保持原目录结构；列表输入生成固定的load_000、load_001…子目录，
其中的文件名、列名与单载荷输出一致；每个载荷分别执行速度扫描和数值核查。
根目录 mixed_cases.csv / mixed_results.xlsx 汇总所有工况，增加load_index、
line_load_n_m、normal_force_n；原case_index为各载荷内编号。
根目录 mixed_loads.csv 和Excel「载荷列表」给出载荷、对应子目录和完成状态。
mixed_diagnostics.json记录运行状态与各载荷状态；中途失败不会报告整批complete，
已完成的子目录保留。缓存继续自动复用，自动Hertz初始域按各载荷分别计算。
汇总mixed_friction.png叠加各载荷曲线；空间剖面与几何图保留在各子目录。
后续saved_forward拟合请明确选择这些单载荷子目录，不能将批量汇总目录作为单一载荷。
fit仍按实验表组织载荷，不使用正向载荷列表作为实验来源。

旧psd配置继续有效。新增可选roughness_parameters.upper/lower块，每块包含
path、scale_c、tail_is_extrapolated。仅按mixed_friction.rough_surface读取选中侧，
另一侧不打开文件。新接口选中path不为空时，旧psd.path应为None；同时给出不一致
来源会报错。当前只实现upper/lower，both是双侧粗糙预留选择，会明确报未实现，
不会合并两谱。下表面PSD路径的预留不意味着双侧粗糙模型已经实现。

macro_contact.deformable_body支持upper；lower作为下一阶段真实下壁变形与黏弹
运动学接口预留，选择会报未实现。旧body_representation继续兼容；新模板不要求
重复填写它。material.model=viscoelastic及master_curve_file仍为预留，不会用固定
弹性模量代算。当前实际计算范围保持单侧粗糙、上壁弹性等效表示、牛顿无滑移。
