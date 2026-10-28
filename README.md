# Persson 干接触计算

项目先实现十模块方案中的模块 1（输入与配置）和模块 2（Persson 干接触），为后续圆柱—平面、一维稳态混合润滑计算提供接触数据库。使用 Python；不需要 MATLAB 许可证。

第一版范围：线弹性、非黏附、各向同性二维 PSD、局部弹性半空间响应。宏观圆柱形变、流体、摩擦、供油及空化尚未接入。

## 下载到自己的电脑

[下载最新版 ZIP](https://github.com/RXU930630/ECL_BD/archive/refs/heads/main.zip)，或者在仓库页面选择 **Code → Download ZIP**。私有仓库需要先登录有访问权限的 GitHub 账户。

这是 Python 源码版，需要安装 Python 3.11 或更新版本；尚未制作 Windows/macOS 独立安装包。代码下载后可在自己的电脑上修改、运行和测试，不依赖本云环境。

解压后进入包含此 README 的目录。Windows PowerShell：

```powershell
py -m venv .venv
.\.venv\Scripts\python.exe -m pip install -r requirements.txt
.\.venv\Scripts\python.exe -m pip install -e .
.\.venv\Scripts\python.exe -m persson_contact examples/dry_contact.json --output results/example --plot
```

macOS/Linux 终端：

```bash
python3 -m venv .venv
.venv/bin/python -m pip install -r requirements.txt
.venv/bin/python -m pip install -e .
.venv/bin/python -m persson_contact examples/dry_contact.json --output results/example --plot
```

计算后查看 `results/example` 中的表格、诊断 JSON 和 PNG 图。示例参数及概率分布的限制见下文。

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

- `contact_table.csv`：按名义面积平均的干接触压力 Pa、平均间隙 m、接触面积比。
- `gap_pdf_*.csv`：微观间隙 m、按名义面积加权的连续密度 1/m、已知间隙大于零的条件密度 1/m。
- `gap_bins_*.csv`：由高斯 CDF 差直接得到各间隙区间的概率，以及概率除以区间宽度得到的平均密度；优先用这些区间概率做积分。
- `diagnostics.json`：材料及数值选项、有效查询范围、分布概率质量和第一矩检查。
- `dry_contact.png`：使用 `--plot` 时输出 PSD、接触面积和平均间隙图。

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
print(response.pressure_pa, response.contact_fraction)
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
