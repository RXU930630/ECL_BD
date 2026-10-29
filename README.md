# Persson 接触与液体泄漏计算

项目按十模块方案组织，目前包括模块1输入配置、模块2 Persson干接触、模块3恒黏度无滑移界面规律，以及模块4的各向同性压力流动部分。已增加模块6参考压力几何、模块7压力边界和模块8—9稳态守恒求解，可计算空间分布及局部卸载后的液体泄漏。使用 Python；不需要 MATLAB 许可证。

第一版范围：线弹性、非黏附、各向同性二维 PSD、局部弹性半空间响应。当前液体泄漏考虑局部流体压力卸载；非局部结构柔度、剪切流动因子、混合摩擦、有限供油、空化及ballistic气体尚未实现。

## 新增：三种几何的空间液体泄漏

[使用说明、算例及 Fortran 对照](docs/liquid_leakage/report.md)。默认开启泄漏，运行 `python main.py`，结果在 `results/leakage`。只需在 `parameters.py` 中修改 `leakage.geometry` 字符串，参数在 `geometry_parameters` 同名分组：

```python
'geometry': 'oring_hertz',       # O-ring：Hertz分布
# 'geometry': 'square_ring',     # Square ring：恒压
# 'geometry': 'prescribed_profile', # 读取单独压力文件
```

四个几何预设均已放在 parameters.py 中；只读取所选类型。也继续支持直接将 geometry 写成带 type 的字典，用于脚本调用。

CSV表头为 `x_m,contact_pressure_pa`，还可增加 `interface_reference_pressure_pa` 列。原Fortran文件选择 `file_format='fortran'`：名称、参考外部流体压力MPa、随后两列x[m]与接触压力MPa，自动转换单位。默认保留真实坐标，不重新等距采样或补零。

入口压力扫描、周向宽度、有效介质/临界通道、是否卸载以及spc处理都在同一 `leakage` 分组。当前 `legacy_spc_plateau` 明确保留原程序的有限导流平台；选择 `physical_zero` 则保留真正阻断。两者可能给出不同泄漏结论。FEM外侧载荷保持输入参考工况，程序只更新粗糙界面内压力及局部接触，未计算非局部FEM柔度。

64项测试、默认111个数据库压力节点和4个入口工况完成；另完成6组原Fortran-500空间卸载核心、三级加密共18次实际参考运行。最细参考网格下泄漏率最大差约0.164%，原包FEM曲线两工况约0.003%和0.005%。这是共享参考压力与导流数据库的空间求解对照；整个原程序的所有物理分支没有因此完成验证。压力数据库加密另使示例泄漏率变化最高约0.286%。

此前接触/流动核心的对照保存在 [压力流动报告](docs/pressure_flow/report.md)，与新泄漏求解共享同一基础。

## 与 Fortran 的干接触对照结果

[完整报告、曲线及CSV数据](docs/fortran_comparison/report.md)。已完成三个E、多个压力的14次批量参考计算，以及1次完全未改源码的独立基准。面积—压力曲线最大相对差约0.065%，平均分离量—压力曲线约0.603%，有效临界间隙查询约1.724%。

**PDF 尾部已修复，最新结果见 [尾部修复报告](docs/gap_tail_correction/report.md)。** 统一 Fortran method 1 的尺度截断、导数和积分网格后，10组参考分布的第一矩最大差约0.0024%，共同区间CDF最大绝对差约9.70e-6。完整40000点密度的最大相对差约4.36%，出现在约1e-12 m的小间隙端；大间隙尾部最大相对差约1.96%；图中异常长尾和小峰已消失。前一份报告保留修复前9.521%的差异，供追溯。模型本身的高斯展宽及尺度截断仍会使 PDF 第一矩不同于能量平均间隙，不能据此宣称完整润滑模型已经验证。

## 日常使用：只改两个文件

安装下文的依赖后即可运行，无需先安装项目本身。请保留完整目录，尤其是 `src/persson_contact`，不要只复制 `main.py`：

1. 在 **`main.py`** 中选择功能开关。
2. 在 **`parameters.py`** 中修改所有工况、材料、PSD 路径与修正参数。
3. 运行 `python main.py`（Windows 也可使用 `.venv\Scripts\python.exe main.py`，macOS/Linux 使用 `.venv/bin/python main.py`）。结果默认在 `results/leakage`。

| main.py 开关 | 当前作用 |
| --- | --- |
| `RUN_CONTACT = True` | 计算共用 Persson 弹性接触数据库 |
| `RUN_LEAKAGE = True` | 计算所选几何的空间泄漏与局部卸载，自动调用接触和流动基础 |
| `RUN_MIXED_FRICTION = False` | 预留混合摩擦入口；求解器尚未实现 |
| `COMPUTE_PRESSURE_FLOW = True` | 接触计算后计算压力流动系数；泄漏开启时自动计算 |
| `COMPUTE_GAP_PDF = True` | 是否导出 Fortran method 1 间隙分布，不关闭流动所需的内部分布计算 |
| `SAVE_PLOTS = True` | 是否生成 PNG 图 |

混合摩擦、黏弹性、复杂流变和壁面滑移尚未实现；选择后会明确报错。开启泄漏会自动计算所需接触和压力流动，即使 RUN_CONTACT=False 或 COMPUTE_PRESSURE_FLOW=False。三个 RUN 开关全关时不执行计算。要只运行原干接触部分，保持 RUN_CONTACT=True，并关闭 COMPUTE_PRESSURE_FLOW 和两个应用开关。

`parameters.py` 中提供粗糙度、材料、名义干接触压力序列、积分及修正参数、输出目录。`interface` 的黏度取原Fortran输入0.1 Pa·s；`leakage` 提供示例尺寸和压差，均不代表实验标定；`mixed_friction` 仍预留。路径始终相对于项目中的 `parameters.py`，不随启动终端的目录改变。为避免把新工况和旧结果混在一起，非空输出目录会被拒绝，请修改 `output_directory` 后重新运行。

压力扫描现在只需设置参考压力：

```python
'contact_pressure': {
    'reference_pa': 1e5,
    'upper_factor': 2.0,
    'number_of_points': 110,
}
```

本例基础配置生成 10 Pa 到 200000 Pa 的 110 个正压力点。泄漏开启时会插入必要的spc节点，并在参考压力曲线超出范围时自动扩充上限；实际范围与新增点数写入diagnostics.json。默认Hertz算例因此为111点。下限固定为参考压力的 1e-4 倍；上限为参考压力乘以上限倍率。点近似对数分布，替换一个邻近节点以准确保留参考压力和总点数。`upper_factor` 至少为 1。仍可使用旧的 `pressures_pa` 列表，但不能与 `contact_pressure` 同时提供。110 点是默认分辨率，不代表对所有工况都已保证插值精度。

术语统一：**A 为 real contact area；A/A0 为 relative contact area**。没有输入名义面积 A0 时，不输出带面积单位的 A。`contact_table.csv` 第三列现命名为 `relative_contact_area`。

临界分离量定义为 `uc=u1`（修正系数 alpha=1）。`options.critical_relative_contact_area=0.42`，`critical_area_tolerance=0.005` 表示 A/A0 的绝对容差，接受接近 0.42 而非严格相等的尺度。先检查 PSD 的倍率范围确实跨过阈值，再定位近似临界尺度；没有跨过时不以最接近的端点代替交点。

`critical_separation.csv` 对每个压力输出临界倍率、实际 relative contact area、偏差、该尺度下的平均间隙、u1 和 uc，长度单位均为 m。u1 使用 A*ubar 对 A 的局部有限体积导数，自动缩小步长进行比较，并输出 `u1_refinement_relative_difference`。只有 `status=found` 时才提供 uc；未跨过阈值、面积分辨率不足、u1 未解析或导数未收敛时保留状态与空值，不填造数值。此表是实验性接触模型输出，尚非泄漏率结果，已做上述有限工况下的 Fortran 对照。

`dry_contact.png` 现有四幅图：PSD、relative contact area—压力、平均间隙—压力、uc—压力。后面三幅图包含整组压力点，并以竖线标记参考压力；不是单一压力下的分布图。诊断 JSON 同时保存完整压力节点和参考压力的临界结果。

材料输入互斥：

```python
# 当前：弹性常数，自动计算 E*=E/(1-nu**2)
'material': {'model': 'elastic', 'young_modulus_pa': 3e6, 'poisson_ratio': 0.5}

# 或直接指定 E*，不再修正一次泊松比
'material': {'model': 'elastic', 'effective_modulus_pa': 4e6}

# 未来：主曲线提供模量，不要求也不允许同时给固定 E 或 E*
'material': {'model': 'viscoelastic', 'poisson_ratio': 0.5,
             'master_curve_file': 'data/master_curve.txt'}
```

以上三种形式只选一种。主曲线文件接口目前只保留路径与材料模式；文件格式、读取及黏弹性接触计算尚未实现，也不会退回弹性计算。原有 JSON 配置及命令行入口继续支持，可加 `--no-pdf` 跳过间隙分布。

## 下载到自己的电脑

[下载最新版 ZIP](https://github.com/RXU930630/ECL_BD/archive/refs/heads/main.zip)，或者在仓库页面选择 **Code → Download ZIP**。私有仓库需要先登录有访问权限的 GitHub 账户。

这是 Python 源码版，需要安装 Python 3.11 或更新版本；尚未制作 Windows/macOS 独立安装包。代码下载后可在自己的电脑上修改、运行和测试，不依赖本云环境。

解压后进入包含此 README 的目录。Windows PowerShell：

```powershell
py -m venv .venv
.\.venv\Scripts\python.exe -m pip install -r requirements.txt
.\.venv\Scripts\python.exe main.py
```

macOS/Linux 终端：

```bash
python3 -m venv .venv
.venv/bin/python -m pip install -r requirements.txt
.venv/bin/python main.py
```

计算后查看 `results/leakage` 中的表格、诊断 JSON 和 PNG 图。示例参数及概率分布的限制见下文。

如果直接使用终端中的 `python`，应使用同一个解释器安装依赖：

```bash
python -m pip install -r requirements.txt
python main.py
```

只有使用旧的 `python -m persson_contact ...` 入口或在其他项目中导入本库时，才需要先执行 `python -m pip install -e .`。若旧版出现 `No module named 'persson_contact'`，也可用此安装命令修复；最新版 `main.py` 已能自行定位本地源码目录。

## 当前云环境中运行

```bash
cd /workspace/ECL_BD
.venv/bin/python -m unittest discover -s tests -v
MPLCONFIGDIR=/workspace/ECL_BD/.cache/matplotlib XDG_CACHE_HOME=/workspace/ECL_BD/.cache .venv/bin/python -m persson_contact examples/dry_contact.json --output results/example --plot
```

在新机器安装（Python 3.11 或更新版本）：

```bash
python -m venv .venv
.venv/bin/python -m pip install -r requirements.txt
.venv/bin/python -m pip install -e .
```

Windows 下将 `.venv/bin/python` 换为 `.venv\Scripts\python.exe`。使用现有 checkout 即可，每个云任务已经隔离，无需另建 Git worktree。

## 示例与研究参数

`examples/data/roughness_psd.dat` 是用户提供的真实 PSD，字节内容未更改。文件最后一个点是用户确认的高波数外推点，与前一点之间做 log–log 线性插值。

**`examples/dry_contact.json` 中的有效模量、压力序列、γ 和分布参数是演示值，不是已确认的研究参数。** 将来应根据原程序配套控制文件设置；不要把演示曲线当作实验结果。

程序输出：

- `contact_table.csv`：按名义面积平均的干接触压力 Pa、平均间隙 m、relative contact area。
- `critical_separation.csv`：各压力下的近似临界尺度、u1、uc 和查找状态。
- `gap_pdf_*.csv`：微观间隙 m、按名义面积加权的连续密度 1/m、已知间隙大于零的条件密度 1/m。
- `gap_bins_*.csv`：由高斯 CDF 差直接得到各间隙区间的概率，以及概率除以区间宽度得到的平均密度；优先用这些区间概率做积分。
- `diagnostics.json`：材料及数值选项、有效查询范围、分布概率质量和第一矩检查。
- `dry_contact.png`：PSD、relative contact area、平均间隙和 uc 四幅图。
- `pressure_flow_table.csv`：几何导流系数K、液压导流系数、压力流动因子、PDF遗漏概率、上界和计算状态。
- `flow_diagnostics.json`：界面模型、有效介质选项、积分收敛和泄漏参数记录。
- `pressure_flow.png`：压力流动系数、正导流系数和泄漏（未开启泄漏时显示遗漏概率）。
- `leakage_cases.csv`：每个入口压力工况的泄漏率、上下界、状态及收敛检查。
- `leakage_profile_*.csv`：对应入口工况的坐标、流体压力、接触压力、平均间隙和relative contact area。
- `leakage_transport_table.csv`：原始/实际使用的导流系数，显式显示spc平台处理。
- `leakage_diagnostics.json`、`leakage.png`：泄漏模型/输入记录和图像；图中压力与间隙显示列表最后一个入口工况，并标注其入口压力。
- `leakage_table.csv`：仅使用旧四参数固定接触接口时生成，继续兼容。

结果文件夹可重复生成，已从 Git 跟踪中排除。

## 模块接口

```python
from persson_contact.config import load_problem
from persson_contact.core import PerssonModel
from persson_contact.table import build_contact_table, evaluate_contact

problem, options = load_problem("examples/dry_contact.json")
model = PerssonModel(problem.psd, problem.effective_modulus_pa,
                     gamma=options.gamma, n_q=options.n_q,
                     pdf_afact=options.pdf_afact,
                     pdf_log_step=options.pdf_log_step,
                     pdf_max_relative_contact_area=options.pdf_max_relative_contact_area,
                     pdf_min_relative_contact_area=options.pdf_min_relative_contact_area,
                     pdf_min_u1_m=options.pdf_min_u1_m)
table = build_contact_table(model, problem.pressures_pa)
lower, upper = table.valid_gap_range_m
response = evaluate_contact(table, (lower + upper) / 2)
print(response.pressure_pa, response.relative_contact_area)
```

PSD 两列固定为 `log10(q)`、`log10(C)`，分别采用 m⁻¹ 和 m⁴ 的二维各向同性谱约定。前六行为名称、q0、qr、q1、h0、H；RMS 从谱积分计算，不直接采用头部 h0。超出 q0/q1 的查询会报错，不静默外推或修改 PSD。

配置中 `material` 可以给出 `effective_modulus_pa`，或给出 `young_modulus_pa` 与 `poisson_ratio`。后者使用 E*=E/(1−ν²)，假定另一侧刚性；已给有效模量时不再修正。`gamma`、`pdf_afact`、`n_q` 必须显式指定。

接触表至少需要三个递增的正压力点。查询在平均间隙有效区间内对 log(p) 做保形单调插值，面积用同一压力回代物理模型，避免面积与分布采用不同压力。查询分布时还返回 `interpolation_mean_gap_error_m`，表示反查压力回代模型后的平均间隙与查询值之差。超出表范围会报错，需显式扩充数据库。

## 流动模块的范围与缺口处理

界面局部规律为 `sigma_local=s^3/(12 eta)`，以及局部Poiseuille+Couette通量。均质化只实现压力驱动部分；局部Couette规律不等于已实现剪切流动因子。

压力流动沿用原Sealing3的各向同性Bruggeman闭合，渗流相对面积共用 `options.critical_relative_contact_area`。`pressure_flow.gap_points` 为内部分布积分网格，独立于 `options.n_gap` 的导出网格；采用CDF分箱和加密比较，未收敛时明确报错。

默认结果将遗漏PDF概率作为零导流相，并输出该概率无限导流时的上界；两者均保留原接触面积、不归一化PDF。上界不是统计置信区间。状态 `conducting` 为默认闭合可导通，`blocked` 为上下界均阻断，`unresolved_connectivity` 为默认值零但上界非零。极端情况下上界无穷，CSV上界留空并置 `upper_bound_unbounded=True`。

新空间泄漏联立恒定通量与局部接触卸载；旧固定状态接口仍使用 `Q=(W/L)*sigma_eff*DeltaP`。压差反向则流动反向，零压差给Q=0。导流表仍是已验证的各向同性直接PDF分支；大间隙渐近切换、各向异性、边界膜、剪切流动因子和完整混合摩擦尚未实现。

## 概率与干接触验证限制

完整间隙分布包含接触概率质量 a 位于零间隙处，以及正间隙的连续部分。此处 P₊ 是已知 s>0 的条件密度，积分为 1；按总名义面积加权的连续密度积分为 1−a。两者不是同一条曲线。

现在采用原 Fortran method 1 的尺度网格和截断：`pdf_log_step=0.01` 是 ln(zeta) 步长，`pdf_max_relative_contact_area=0.999`、`pdf_min_relative_contact_area=0.0001`、`pdf_min_u1_m=1e-12` 是原程序的积分范围。全部集中在 `parameters.py` 的 `options` 中，旧 JSON 省略时采用这些默认值。`n_q` 继续控制能量平均间隙和临界间隙计算，PDF 的尺度分辨率改由 `pdf_log_step` 单独控制。与 Fortran 对照时必须使用相同的尺度步长和截断参数。

PDF 使用解析 A' 和中心差分 ubar' 计算 u1，按照中心差分面积变化加权，避免旧有限体积公式在 A≈1 时产生巨大 u1。能量核函数仍使用独立准确积分，没有为贴合原程序的近似查表而降低精度。原程序的尺度截断作用于混合分布的组成；程序不会直接剪掉输出间隙轴上的曲线。

`unresolved_continuous_mass` 记录积分范围之外的非接触概率，另分解为 `excluded_high_area_mass`、`excluded_low_area_mass`、`excluded_small_u1_mass` 和 `boundary_unresolved_mass`。不把缺口当成接触，也不强行归一化。PDF 采用右矩形尺度积分，接触面积仍采用独立准确的 PSD 积分，因此需另计有符号的 `quadrature_mass_error`：

```
已解析概率 + 未解析概率 = 1 - relative_contact_area + quadrature_mass_error
```

返回的条件密度仍除以完整非接触概率 1-a，其保留部分的积分不必等于1。最小间隙附近可出现极窄峰，PDF CSV 不输出 s=0 的密度值；区间概率优先使用解析 CDF 差得到的 `gap_bins`，并检查概率缺口、尺度积分误差和第一矩诊断。若诊断提示初始尺度面积小于0.99，需要减小 `pdf_log_step` 或补充低波数 roll-off；不能将欠分辨率的分布作为已验证结果使用。

10组实际 Fortran 输出覆盖 E=1.5/3/9 MPa、p=0.01/0.05/0.1 MPa，以及额外的0.005尺度步长。最新 [曲线和逐工况数据](docs/gap_tail_correction/report.md) 包含极低密度尾部及第一矩比较，回归测试保留原始参考采样和来源哈希。用于流动或摩擦积分时，仍须明确积分域和截断约定，并检查相应高阶矩；与 Fortran 的复现一致不代表近似分布严格满足能量平均间隙等式。

已用GNU Fortran完成原参数gamma=0.4、E=3 MPa、p=0.1 MPa的完全未改源码运行。批量副本仅加速实数误差函数并增加PDF输出密度，源代码差异和与原版锚点的核对保存在上述报告中。详细公式、例程位置和最初输入分析见 [源码分析与模块计划](docs/fortran_review_and_module_plan.md)。超时、失败及主动停止的试算不计为成功。

参考脚本需要另行提供原始Fortran ZIP和编译器；GitHub中不包含完整Fortran包和编译器。`--real-erf`为可选构建副本优化。`scripts/benchmark_fortran_matrix.py`、`scripts/analyze_fortran_matrix.py`用于多E/压力比较；`scripts/verify_gap_tail.py`复核修复后的PDF尾部，命令见相应报告。
