# Persson 干接触计算

项目先实现十模块方案中的模块 1（输入与配置）和模块 2（Persson 干接触），为后续圆柱—平面、一维稳态混合润滑计算提供接触数据库。使用 Python；不需要 MATLAB 许可证。

第一版范围：线弹性、非黏附、各向同性二维 PSD、局部弹性半空间响应。宏观圆柱形变、流体、摩擦、供油及空化尚未接入。

## 日常使用：只改两个文件

安装下文的依赖后即可运行，无需先安装项目本身。请保留完整目录，尤其是 `src/persson_contact`，不要只复制 `main.py`：

1. 在 **`main.py`** 中选择功能开关。
2. 在 **`parameters.py`** 中修改所有工况、材料、PSD 路径与修正参数。
3. 运行 `python main.py`（Windows 也可使用 `.venv\Scripts\python.exe main.py`，macOS/Linux 使用 `.venv/bin/python main.py`）。结果默认在 `results/contact`。

| main.py 开关 | 当前作用 |
| --- | --- |
| `RUN_CONTACT = True` | 计算共用 Persson 弹性接触数据库 |
| `RUN_LEAKAGE = False` | 预留泄漏入口；求解器尚未实现 |
| `RUN_MIXED_FRICTION = False` | 预留混合摩擦入口；求解器尚未实现 |
| `COMPUTE_GAP_PDF = True` | 是否额外计算实验性间隙分布 |
| `SAVE_PLOTS = True` | 是否生成 PNG 图 |

若启用未实现的应用或黏弹性模型，会在计算、创建结果之前明确报错，不会输出伪造的泄漏率或摩擦系数。两个应用开关同时开启也不代表已经实现耦合求解。三个计算开关全部关闭时不执行任何计算。

`parameters.py` 中提供粗糙度、材料、名义干接触压力序列、积分及修正参数、输出目录。`leakage` 与 `mixed_friction` 分组暂仅预留，`None` 表示尚未给定。路径始终相对于项目中的 `parameters.py`，不随启动终端的目录改变。为避免把新工况和旧结果混在一起，非空输出目录会被拒绝，请修改 `output_directory` 后重新运行。

压力扫描现在只需设置参考压力：

```python
'contact_pressure': {
    'reference_pa': 1e5,
    'upper_factor': 2.0,
    'number_of_points': 110,
}
```

本例自动生成 10 Pa 到 200000 Pa 的 110 个正压力点。下限固定为参考压力的 1e-4 倍；上限为参考压力乘以上限倍率。点近似对数分布，替换一个邻近节点以准确保留参考压力和总点数。`upper_factor` 至少为 1。仍可使用旧的 `pressures_pa` 列表，但不能与 `contact_pressure` 同时提供。110 点是默认分辨率，不代表对所有工况都已保证插值精度。

术语统一：**A 为 real contact area；A/A0 为 relative contact area**。没有输入名义面积 A0 时，不输出带面积单位的 A。`contact_table.csv` 第三列现命名为 `relative_contact_area`。

临界分离量定义为 `uc=u1`（修正系数 alpha=1）。`options.critical_relative_contact_area=0.42`，`critical_area_tolerance=0.005` 表示 A/A0 的绝对容差，接受接近 0.42 而非严格相等的尺度。先检查 PSD 的倍率范围确实跨过阈值，再定位近似临界尺度；没有跨过时不以最接近的端点代替交点。

`critical_separation.csv` 对每个压力输出临界倍率、实际 relative contact area、偏差、该尺度下的平均间隙、u1 和 uc，长度单位均为 m。u1 使用 A*ubar 对 A 的局部有限体积导数，自动缩小步长进行比较，并输出 `u1_refinement_relative_difference`。只有 `status=found` 时才提供 uc；未跨过阈值、面积分辨率不足、u1 未解析或导数未收敛时保留状态与空值，不填造数值。此表是实验性接触模型输出，尚非泄漏率结果，也未完成 Fortran 对照。

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

计算后查看 `results/contact` 中的表格、诊断 JSON 和 PNG 图。示例参数及概率分布的限制见下文。

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

结果文件夹可重复生成，已从 Git 跟踪中排除。

## 模块接口

```python
from persson_contact.config import load_problem
from persson_contact.core import PerssonModel
from persson_contact.table import build_contact_table, evaluate_contact

problem, options = load_problem("examples/dry_contact.json")
model = PerssonModel(problem.psd, problem.effective_modulus_pa,
                     gamma=options.gamma, n_q=options.n_q,
                     pdf_afact=options.pdf_afact)
table = build_contact_table(model, problem.pressures_pa)
lower, upper = table.valid_gap_range_m
response = evaluate_contact(table, (lower + upper) / 2)
print(response.pressure_pa, response.relative_contact_area)
```

PSD 两列固定为 `log10(q)`、`log10(C)`，分别采用 m⁻¹ 和 m⁴ 的二维各向同性谱约定。前六行为名称、q0、qr、q1、h0、H；RMS 从谱积分计算，不直接采用头部 h0。超出 q0/q1 的查询会报错，不静默外推或修改 PSD。

配置中 `material` 可以给出 `effective_modulus_pa`，或给出 `young_modulus_pa` 与 `poisson_ratio`。后者使用 E*=E/(1−ν²)，假定另一侧刚性；已给有效模量时不再修正。`gamma`、`pdf_afact`、`n_q` 必须显式指定。

接触表至少需要三个递增的正压力点。查询在平均间隙有效区间内对 log(p) 做保形单调插值，面积用同一压力回代物理模型，避免面积与分布采用不同压力。查询分布时还返回 `interpolation_mean_gap_error_m`，表示反查压力回代模型后的平均间隙与查询值之差。超出表范围会报错，需显式扩充数据库。

## 概率与验证限制

完整间隙分布包含接触概率质量 a 位于零间隙处，以及正间隙的连续部分。此处 P₊ 是已知 s>0 的条件密度，积分为 1；按总名义面积加权的连续密度积分为 1−a。两者不是同一条曲线。

实际离散计算另报告 `unresolved_continuous_mass`（未解析的非接触概率），不把它当成接触，也不把它强行分给已有曲线。返回的密度/CDF 仅包含已解析部分：接触概率 + 已解析概率 + 未解析概率 = 1；存在未解析量时，返回的条件密度积分小于 1。最小间隙附近可出现极窄峰，普通梯形积分点值密度可能严重失真；PDF CSV 不输出 s=0 的密度值，区间概率请使用 `gap_bins`，并同时检查未解析量。

间隙分布采用原 Fortran 启发的倍率混合、反射高斯近似，并使用有限体积方式离散面积变化。**分布目前属于实验性实现，不是逐行移植或已与原程序核对的版本。** 高斯展宽、数值积分和有限取样范围可能导致第一矩与独立平均间隙计算不同；程序保留并报告偏差，不强制重标定分布来消除它。用于后续流动或摩擦积分前，必须检查这些诊断并完成网格/截断及 Fortran 对照验证。

在随附演示参数下，分布解析第一矩偏差最大约 6.35%，未解析非接触概率最大约为名义面积的 2.94×10⁻⁴；这不是研究工况的已验证精度。面积与平均间隙的独立数值检查不代表分布已经通过同等验证。

后续收到完整 Fortran 包后，已用 GNU Fortran 编译，并完成原程序的面积扫描及单独 γ=1 的平均间隙测试；原参数 γ=0.4 的分离运行未完成，已按用户要求停止。不将超时或停止的算例计入成功结果，也尚未验证完整间隙分布与原程序一致。详细公式、例程位置和最初输入分析见 [源码分析与模块计划](docs/fortran_review_and_module_plan.md)。

`scripts/run_fortran_reference.py` 与 `scripts/compare_fortran_reference.py` 用于后续自行运行参考算例，需要另行提供原始 Fortran ZIP 和编译器；GitHub 中不包含用户上传的完整 Fortran 包、编译器或生成结果。`--real-erf` 是可选的构建副本优化，会记录源码差异，尚未完成该优化副本的端到端验证。
