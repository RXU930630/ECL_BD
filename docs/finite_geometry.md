# 直接输入有限圆柱曲面的尺寸

`macro_contact.geometry='finite_cylindrical_pad'`用于已知曲率的有限胶块。泄漏任务仍读取自己的`leakage.geometry`，不受这里影响。宏观干接触、固体摩擦、混合润滑和剪切强度拟合共用此几何。

```python
'macro_contact': {
    'geometry': 'finite_cylindrical_pad',
    'geometry_parameters': {
        'finite_cylindrical_pad': {
            'outer_radius_m': 0.17,
            'axial_length_m': 0.10,
            'sliding_length_m': 0.03,
            'length_definition': 'projected',
            'rubber_thickness_m': 0.02,
        },
    },
    'normal_force_n': 100.,
    # 保留其余elastic_model、solver、contact_table等设置
},
```

曲率半径是粘贴成形后外表面的半径。轴向长度是垂直于滑动方向的长度。`projected`表示流向长度在平面上的投影，本例自动生成`[-0.015,0.015] m`。如果3cm是沿圆弧测量，选择`arc_length`，边界自动改为`±R*sin(length/(2R))`。表面采用真实圆弧矢高`g(x)=x²/(R+sqrt(R²-x²))`。

`rubber_thickness_m`仅记录实验尺寸。弹性仍采用既有平面应变半空间模型，不引入有限厚度、刚性背衬粘接约束、轴向端部效应或贴合预应力。有限表面范围和有限尺寸弹性是不同的近似层次。

不需要填写或估计`domain_half_width_m`；有限曲面模式忽略它。旧`cylinder_plane`模式继续使用局部抛物线和人工指定的远场截断范围。

## 载荷与拟合

`normal_force_n`填总力N，程序用轴向长度B自动换算`W'=F/B`。也可改用`line_load_n_m`直接填N/m；两项不能同时给非空值。`contact_friction.transverse_width_m`在有限模式留`None`即可，若手填则必须等于轴向长度，以免重复宽度产生冲突。

本例40、100、160N分别对应400、1000、1600N/m。一次拟合可以包含全部载荷：每行填写速度m/s、摩擦系数和总力N，例如：

```csv
speed_m_s,friction_coefficient,normal_force_n
0.000001,填写实验mu,40
0.000001,填写实验mu,100
0.000001,填写实验mu,160
```

上面是列格式说明，含占位文字，不能直接计算。正式文件用Excel打开`inputs/experiment_template.xlsx`后填写并另存；`line_load_n_m`可留空。CSV逐行的载荷覆盖默认100N。拟合以CSV速度为准，无需同步修改`mixed_friction.lower_wall_velocities_m_s`。

## 物理边缘与输出

计算范围止于实际曲面边缘，不自动延长胶块。混合润滑在两端**单元面**施加储液压力，流入端为满液；相邻单元中心的压力无需等于边界压力。有限模式仍检查载荷、几何、质量守恒和空化互补条件，但不使用旧模式的远场边缘压力衰减阈值。

`mixed_boundaries.csv`记录每个速度、两侧物理坐标、储液压力与边界流量。`mixed_diagnostics.json`/`macro_diagnostics.json`记录几何解释与弹性近似。变形后的平均轮廓继续输出到几何CSV和图片。

257单元是初始数值设置；试算后可用513单元检查接触面积和拟合剪切强度的变化。有限模式加密网格时保持实际长度不变。改变几何或流体后需重算正向解；不同几何的旧结果不能混入同一次拟合。

## 本次试运行包

带`silicone_cylinder`预设的0.2.0rc3包已填写上述几何、E=1.8MPa、ν=0.5、η=1Pa·s，启用独立拟合，接触剪切强度保持未知。1000cSt是运动黏度；1Pa·s是本次采用的动力黏度近似，之后可据油品密度修正。

包内不含实测摩擦数据或已计算结果。填写`inputs/experiment.xlsx`后运行`python main.py`。无需先单独运行混合润滑；拟合任务会按CSV中的各载荷和速度计算。
