# 正向混合摩擦：统一谱与原版输入

下载的 mixed-forward rc5 包默认开启 `RUN_MIXED_FRICTION=True`，关闭泄漏、独立干接触和拟合。直接在项目根目录运行 `python main.py`。独立 `COMPUTE_SHEAR_FLOW=False` 正确，混合润滑内部一直使用剪切流动。

## 统一基准

E=3 MPa，ν=0.5，η=0.1 Pa·s，R=0.005 m，W′=2 N/m；B=0.01 m只用于总量换算；tau=1 MPa为公开假设。上圆柱粗糙且固定，下平面沿+x运动；无滑移、充分供油、两端零表压。速度为1e-8、1e-7、1e-6、1e-5、1e-4、1e-3 m/s。

`inputs/mixed_benchmark_psd.dat` 是有限波数合成谱，不是用户实验PSD。其字节与本目录 `fortran_inputs/IN_ROADTYPE.q0.qr.q1.h0.H.1logWavevector.2logPowerSpectrum` 完全相同，SHA256：`7a442e27a66ad810f7cdf934cc219e18180fc5df9cd2ca7ef36f2fe92ec041e1`。

Python 主计算129单元，接触表601点、输运表41基点、n_q=128、pdf_log_step=0.1；Fortran Nx=128、流动表12压力点、谱步长0.1。两版相同的区间为±0.0006770275002573076 m。节点数不要求相同，但需分别核查精度。

该有限区间是上一轮原版对照的固定基准，已知最高速度对扩大区间敏感。rc5会保存相同主结果，并另外执行网格、表/积分、区间检查；默认容限下出现 `needs_attention` 是预期的数值诊断。请读取Excel中的超限项，不能把主计算完成等同于数值或物理验证通过。

## 原版服务器运行

使用最开始ZIP中的 `AI.ContactMechanics.14.09.2026.f90`，SHA256应为 `99d3b2e2c47a663ffe8b65bbdb7415db5eeebb434debd0ef404006600e861854`。本包不分发原始Fortran源码。将本目录 `fortran_inputs` 的全部文件复制到一个新工作目录，在其中编译和运行原版。例如GNU Fortran 14.x：

```bash
gfortran -O3 -mcmodel=large -ffree-line-length-none -std=legacy -fallow-argument-mismatch /你的路径/AI.ContactMechanics.14.09.2026.f90 -o original
./original > run.log 2>&1
```

保留全部原生输出。BoF预计算可能持续数十分钟。主要设置为 `iflag011=31`、`irough=0`、`imoving=0`，载荷切换阈值0.4，恒黏度、无滑移。原版 STOP 也可能返回0；需检查有内容的ERROR文件、完整速度行及残差。

原版 LOW/WITH/FrictionFlowfactor 文件的实际前8列为：v、总线载荷、固体线载荷、流体线载荷、原版相对面积、带符号mu_s、滑动流体mu、压力驱动流体mu；开头没有time。对运动下平面的阻力，前三个摩擦分项分别取负，再相加；不要逐项取绝对值。

原版固体摩擦所用面积是平均接触压力估计，Python积分局部接触面积，两者定义有差别。原版hmin搜索可能误读前Nx有效点以外的数据而报0，且一些上部物体的辅助受力列可能为NaN；不要将这些列用作物理对照目标。六点基准未覆盖完整高速分支。

## 正向输出与实验工况

`results/mixed_benchmark_rc5/mixed_results.xlsx` 含全部速度的固体、流体滑动剪切、流体压力剪切、总mu、承载、面积与膜厚，以及数值核查。CSV、空间剖面和PNG也保留。`numerical_checks` 中是辅助解，主结果不自动替换。

要使用用户的17 cm曲率半径、10 cm轴向、3 cm滑动方向、2 cm厚的胶块与1 Pa·s硅油，可将 `examples/silicone_forward_parameters.py` 的内容复制到根目录 `parameters.py`。该示例采用40N（400N/m）和仓库原有PSD，tau仍是1MPa的明确假设；请按本次与原版对照的实际载荷和剪切强度修改。有限物理尺寸不会被程序扩大。

数值设置、容限、状态和拟合诊断见 `docs/mixed_numerical_checks.md`。空化的物理闭合保持rc4，尚未增加Noël混合物模型、气泡动力学或毛细力。
