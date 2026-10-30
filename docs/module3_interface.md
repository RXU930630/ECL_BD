# 模块3：恒黏度牛顿流体与Navier滑移

本模块提供局部平行间隙内的压力流、运动带液和两侧壁面剪切应力，并将滑移压力导流接入现有液体泄漏工作流。粗糙表面的有效剪切流动因子、混合摩擦、压温相关黏度、非牛顿流体及温度场求解尚未实现。

## 参数入口

main.py继续选择任务。parameters.py的CASE中，流体、流变模型、壁面规律、几何分别使用独立选择器；每个选择器的备选参数放在对应预设块中：

```python
'interface': {
    'rheology': 'newtonian',
    'rheology_parameters': {
        'newtonian': {'viscosity_pa_s': 0.1},
    },
    'wall_slip': 'no_slip',       # 或navier
    'slip_parameters': {
        'no_slip': {},
        'navier': {
            'lower_slip_length_m': 1e-8,
            'upper_slip_length_m': 1e-8,
        },
    },
},
```

只读取所选预设。旧的interface.viscosity_pa_s仍可读取，但不能与rheology_parameters同时给出。未实现的流变或滑移名称会报错，不会自动回退为牛顿无滑移。

上下壁面指膜厚方向z=0与z=h；流向坐标x的小坐标端仍为入口。滑移长度单位m，允许0；两侧均为0时严格恢复无滑移。

切换为非零Navier滑移时，将CASE的leakage.compute_equivalent_gas设为False。原辅助气体估算采用无滑移液体等效间隙关系，不能套用滑移液体结果；程序会在创建输出前检查这一组合。独立ballistic气体路径使用自己的局部规律，液体Navier参数不参与气体求解。

## 局部规律与量纲

设动力黏度eta，间隙h，两侧速度U0、U1，压力梯度G=dp/dx，两侧滑移长度b0、b1。边界条件：

```
u(0)-U0 = b0 du/dz(0)
u(h)-U1 = -b1 du/dz(h)
eta d²u/dz² = G
```

每单位横向宽度的体积通量q，单位m²/s：

```
D = h+b0+b1
K_local = h³[1+3(b0+b1)/D] + 12h² b0 b1/D   [m³]
q = -K_local G/(12 eta)
    + h[(h/2+b1)U0+(h/2+b0)U1]/D
```

两侧等滑移b时，压力导流相对无滑移增强1+6b/h；只有一侧滑移且长度趋于无穷时，增强趋于4。两侧滑移不等时，运动带液项不能直接写成h(U0+U1)/2。

定义A=[eta(U1-U0)-G h(h/2+b1)]/(eta D)，则u(z)=U0+b0 A+A z+Gz²/(2 eta)。输出剪切应力采用“流体作用在壁面上的沿x方向牵引”：

```
tau_lower = eta A
tau_upper = -eta A-Gh
```

因此tau_lower+tau_upper=-Gh。零间隙时q=0，流体速度和剪切应力为NaN；没有流体层，不能由此构造固体接触摩擦。

## 局部API与后续流变扩展

```python
from persson_contact.interface import NewtonianNavierSlip

law = NewtonianNavierSlip(
    viscosity_pa_s=0.1,
    lower_slip_length_m=1e-8,
    upper_slip_length_m=2e-8,
)
response = law.local_response(
    gap_m=1e-6,
    pressure_gradient_pa_m=-1e8,
    lower_wall_velocity_m_s=0.0,
    upper_wall_velocity_m_s=0.1,
    fluid_pressure_pa=1e5,
    temperature_k=298.15,
)
# response.flux_per_width_m2_s
# response.lower_wall_shear_stress_pa
# response.upper_wall_shear_stress_pa
```

viscosity(shear_rate_s_inv, fluid_pressure_pa, temperature_k)已预留剪切速率大小、流体绝对压力与温度参数；本版本始终返回配置的恒定黏度。压力单位Pa、温度K、剪切速率1/s。流体绝对压力与接触数据库中的名义接触压力是不同物理量。

预留参数不代表已经支持压力/温度相关黏度。今后接入Eyring、自定义eta(gamma_dot,p,T)等模型时，还需相应扩展局部流动求解、均质化及宏观方程，不能只替换常数返回值。

可直接运行 `python examples/local_interface.py` 查看局部流量与剪切应力。局部速度输入尚未接入粗糙接触的宏观剪切输运；main中的泄漏目前仍为静止壁面的压力驱动问题。

## 导流与验证

有效介质在各PDF间隙点使用完整K_local求根，保持K的m³单位及phi_p=K/h³定义。低压延拓按完整局部导流比值缩放，临界通道也使用所选滑移规律。空间阻力仍为12eta/K。

测试包括：独立求解平行板边界值问题并积分速度分布；双壁等滑移及单壁极限；不对称滑移；壁面互换、整体速度平移；流体与壁面力平衡；黏性及滑移耗散；零间隙；均匀间隙/接触相混合；低压延拓；临界通道；完整泄漏流程和下载源码后的参数切换。既有液体及ballistic参考回归继续保留。

本次发布前运行完整测试集：91项通过。局部示例已实际运行；默认无滑移和Navier滑移的源码入口、液体与ballistic分支均有工作流覆盖。
