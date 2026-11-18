# 原始 Fortran / Python 混合润滑对照输入

这是独立的数值核查算例。Python 求解代码冻结在 rc4（`026c8951add03e4440a55234849f236a1bf303a9`），不改写现有模块，也不覆盖日常 `main.py` / `parameters.py`。当前为**收敛预检输入**，不能称为已验证的完整理论复现或实验稳定版。纯净程序和实验输入保留在仓库根目录；本目录的合成谱仅用于验证。

## 共同物理输入

| 量 | 数值 |
|---|---:|
| 线弹性 E / ν / E* | 3 MPa / 0.5 / 4 MPa |
| 圆柱半径 R | 0.005 m |
| 线载荷 W′ | 2 N/m |
| 牛顿黏度 η | 0.1 Pa·s |
| 已知界面剪切强度 τ | 1 MPa；本算例不做拟合 |
| 下平面速度 | 1e-8、1e-7、1e-6、1e-5、1e-4、1e-3 m/s |
| 上圆柱速度 | 0；在圆柱随动坐标中计算 |
| 粗糙度 | 上圆柱粗糙、下平面光滑；各向同性 |
| 两端表压 / 空化表压 | 0 / 0 Pa；充分供油 |
| 滑移、体黏弹滞后、毛细力、非牛顿 | 均关闭 |
| 参考轴向宽度 B | 0.01 m；仅用于将每单位宽度量换算成总量 |
| Hertz 参考半宽 aH | 5.6418958355e-5 m |
| 计算区间 | ±12 aH = ±0.0006770275002573 m |

上表是合成基准，**不是**用户的 17 cm 圆柱、1 Pa·s 硅油实验。选较低载荷是为了先检查连通流动支路。τ=1 MPa 是公开指定的测试常数，得到较大 μ 不代表拟合值。

`benchmark_psd.dat` 与 Fortran 输入中的 `IN_ROADTYPE.q0.qr.q1.h0.H.1logWavevector.2logPowerSpectrum` 字节完全相同：q0=1e3、qr=3e4、q1=3e6 m⁻¹，H=0.8，平台后 C∝q⁻³·⁶，hrms=1 μm。PSD 数据列为 log10(q)、log10(C)，C 的 SI 单位为 m⁴；hrms按谱积分归一化。仅给有限波数范围，不在外侧额外延伸谱。SHA256、源文件校验值、各输入文件校验值见 `manifest.json`。

## 运行 Python

在解压后的**项目根目录**打开终端：

```text
python -m pip install -r requirements.txt
python benchmark/run_python.py
python benchmark/export_python_excel.py benchmark_runs/python
```

直接读取 `benchmark/python_parameters.py`，不用修改日常开关。输出到 `benchmark_runs/python`；重复运行前为旧结果改名。`forward_components.xlsx` 将固体、流体滑动剪切、流体压力驱动剪切和总 μ 分开，并保留原始数据表及 Excel 图。

两版运行完成后，可在同一台电脑汇总原生输出：

```text
python benchmark/compare_outputs.py --python benchmark_runs/python --fortran benchmark_runs/fortran --fortran-kind pristine --output benchmark_runs/comparison
```

生成 `comparison.xlsx`、`comparison.json` 和 `comparison.png`，保留原版 hmin 导出问题以及辅助列 NaN 的记录。如果来源是另行标注的 ERF 诊断加速副本，必须使用 `--fortran-kind real-erf-diagnostic`，不能标作 pristine。

## 运行原版 Fortran

使用本对话**最开始的 ZIP** 中的 `AI.ContactMechanics.14.09.2026.f90`，不是后来泄漏计算的 `GOOD.GNU...f90`。此包提供全部配套 IN 文件，不另外分发原始 Fortran 源码。GNU Fortran 建议 14.x：

```text
python benchmark/run_fortran.py --source /你的路径/AI.ContactMechanics.14.09.2026.f90
```

源文件 SHA256 不匹配会停止。脚本在独立的 `benchmark_runs/fortran` 中编译并运行；不改源码。源代码中的 BoF 查表预计算可持续数十分钟。若服务器上使用已有可执行文件，也可将 `fortran_inputs` 全部文件复制到一个新工作目录，确认源码版本、编译选项一致后运行；请保留原生输出和日志。

关键 Fortran 选择：`iflag011=31`（有效介质压力流动，并由 A′ 链计算逆间隙），`iflag000=0`，`iflag002=1`，其余任务关闭；`irough=0`、`imoving=0`，`loaddivision=0.4`，`bthin=1e30`，`Cuniversal=0`，`uboundary=1e-9`，滑移长度 0，能量修正 γ=0.4。这些数值在两版输入中已填好。

**Fortran 返回码 0 不足以证明成功。** 原版的 STOP 也可能返回 0；需检查有内容的 ERROR 文件、每个速度是否都有有限输出、原生迭代残差以及网格细化。原版会无条件创建几个空 ERROR 文件。附带脚本检查执行及主比较列的完整性，单独报告辅助列的 NaN，不自动声称收敛。

## 先检查哪些输出

Python `mixed_cases.csv` 每行一个速度。固体/流体承载、真实面积 Ar/B、最小平均间隙、最大流体压力、下平面固体力、流体滑动/压力剪切力、总 μ、质量/载荷残差均可分别绘图。`mixed_profile_*.csv` 给出空间分布 h、pc、pf、A/A0、饱和度 θ 和壁面牵引等；`mixed_diagnostics.json` 记录模型与求解诊断。

拟合任务中，`fit_results.xlsx` 才使用拟合后的 τ。其 `forward/load_*` 原始表通常使用 `tau_ref=1 MPa`，不要把该表的原始总 μ 当成拟合后的总 μ。本基准直接指定相同 τ，只运行正向解。

### 原版列序必须按 WRITE 核对

此处针对源 SHA256 指定的**纯弹性** `MixedLubricationCylinder`，不要根据冗长文件名推断列序。

| 原版输出 | 实际列序 |
|---|---|
| `out.MSLUBE.LOW.WITH.FrictionFlowfactor...` | **16 列**：v、W′、Ws′、Wf′、原版相对接触面积、带符号 μs、流体滑动剪切 μ、流体压力剪切 μ、8 个其他物体/斜面受力项。开头**没有 time**。 |
| `out.MSLUBE.1time.2vel.3Fload...` | **7 列**：time、v、W′、Ws′、Wf′、原版相对接触面积、abs(μs)。 |
| `out.MSLUBE.1time.2vel.3minimumSeparation...` | **5 列**：v、hmin、xhmin、baruc、uboundary。开头没有 time。 |
| `out.MSLUBE.1time.2vel.3x.4baru...` | **7 列**：v、x、h、pc、pf−pa、滑动剪切应力、负的压力剪切应力。没有 time。 |

本工况报告**运动下平面**上的正向阻力：Fortran LOW/WITH 文件第 6–8 列先保留符号，再分别取负转换成阻力贡献；μtot=−(c6+c7+c8)。压力剪切可能助推，不能逐项取绝对值。上圆柱的斜面压力力另列，不能随意加到下平面上。

原版用于固体摩擦的面积是 `Ar/B = 2*aH*A(Ws′/(2*aH))`；Python 是 `∫A(pc(x))dx`。原版第 5 列的“相对面积”与全计算域平均面积不是同一定义。可先从 `abs(μs)*W′/τ` 恢复原版原生等效 Ar/B，再单列 Python 积分面积；两者不能直接以同一个 A/A0 标题比较。原版空间导出循环还可能写出计算域外额外行，解析时限定物理 x 区间并检查有效值。

源文件的空间导出和最小间隙搜索循环均写到 `2*Nx+1`，而低速算法本体只计算前 `Nx` 个值。本预检中后续行可能为零，导致原版原生 hmin 文件报告 0；请按每个速度的**前 Nx 行**重算有效 hmin，保留原生 0 作为导出问题记录。此外，LOW/WITH 第 9–12 列的上圆柱斜面受力可能为 NaN；不能将这些值冒充有效的整体受力结果。下平面本次主比较只使用第 1–8 列。应保留全部原生输出，分别检查这些问题。

## 比较范围和预检状态

两版共同物理量相同，但算法并非逐行翻译：Python 使用守恒 PDF 积分、受限插值、强耦合求解、质量守恒空化；原版使用自己的谱/PDF求积、样条与外推、分步松弛和低/高速算法切换。逆间隙截断、φf 下限、远间隙处理也存在差别。

因此比较顺序是：谱矩 → 干接触 p(h)、A(p) → 输运/摩擦因子 → 同速度的 h、pc、pf、承载、面积 → 摩擦分量 → 总 μ。先定位差异，不能只看总曲线接近就宣布通过。

当前采用粗预检网格：Python 129 单元/41 个输运基点；Fortran Nx=128/12 个压力基点；谱对数步长 0.1。原版 Nx 取偶数以包含 x=0，避免极小初始 Hertz 区完全落在节点之间。Python 边缘压力**诊断门槛**为 0.005，较默认宽松，以先比较同一有限区间；未改变控制方程或求解残差容限。后续必须细化网格、表、谱步长并扩大区间，尤其不能把这个门槛当作远场收敛证明。

先前 Nx=129、从 1e-5 m/s 开始的原版预检在第一点报 `ERROR.renormalization.factor.negative`，未产生有效混合润滑输出。此候选输入改为 Nx=128，并从 1e-8 开始；最终原版运行结果和收敛结论需以独立报告为准。
