"""集中输入：全部使用 SI 单位；路径相对于本文件所在目录。

材料/修正参数取自所提供 Fortran 包，配合用户的 block PSD。
这是研究模型的输入起点，不代表已经校准的实验工况。
"""

PARAMETERS = {
    # 图注用英文，避免未安装中文绘图字体的电脑出现缺字。
    'description': 'Elastic contact: archive material and supplied PSD; not an experimental fit.',
    'psd': {
        'path': 'examples/data/roughness_psd.dat',
        'scale_c': 1.0,                 # PSD 幅值乘数；不修改原文件
        'tail_is_extrapolated': True,    # 记录已给定的外推尾部，不自动外推
    },
    'material': {
        'model': 'elastic',
        'young_modulus_pa': 3.0e6,       # 杨氏模量 Pa
        'poisson_ratio': 0.5,
    },
    # 如已知 E*，可将整个 material 改为：
    # {'model': 'elastic', 'effective_modulus_pa': 4.0e6}
    # 未来黏弹性输入采用以下互斥形式，不再给固定杨氏/有效模量：
    # {'model': 'viscoelastic', 'poisson_ratio': 0.5,
    #  'master_curve_file': 'data/master_curve.txt'}
    # 目前只校验该接口，尚不读取主曲线，也不进行黏弹性计算。
    'contact_pressure': {
        'reference_pa': 1.0e5,       # 关注的名义接触压力 Pa，保证包含该点
        'upper_factor': 2.0,         # 扫描上限 = 参考压力 × 倍率
        'number_of_points': 110,     # 近似对数分布，自动正下限 = 参考压力 × 1e-4
    },
    'options': {
        'gamma': 0.4,       # 弹性能量修正，不是摩擦系数
        'n_q': 256,         # 平均间隙/临界间隙积分分辨率，整数 >= 32
        'pdf_afact': 1.0,   # 分布展宽修正
        'pdf_log_step': 0.01,  # PDF 的 ln(zeta) 步长；与 Fortran 对照时须相同
        'pdf_max_relative_contact_area': 0.999,  # 原程序上截断，排除 A≈1 的不适定 u1
        'pdf_min_relative_contact_area': 0.0001, # 原程序下截断
        'pdf_min_u1_m': 1e-12,  # 原程序 u1min；不等于强行裁去输出曲线的尾部
        'n_gap': 600,       # 微观间隙输出网格，整数 >= 32
        'critical_relative_contact_area': 0.42,
        'critical_area_tolerance': 0.005,  # A/A0 的绝对容差，非相对百分数
    },
    'output_directory': 'results/contact',  # 非空时拒绝覆盖；新工况请换目录
    # 模块3：目前仅支持恒黏度、无滑移。黏度取原 Fortran 输入，不是实验标定。
    'interface': {
        'rheology': 'newtonian',
        'wall_slip': 'no_slip',
        'viscosity_pa_s': 0.1,
    },
    # 模块4：第一阶段为各向同性压力流动；剪切流动因子稍后实现。
    # 渗流面积沿用 options.critical_relative_contact_area。
    'pressure_flow': {
        'closure': 'fortran_isotropic_bruggeman',
        'gap_points': 2049, # 内部 CDF 积分网格；与 n_gap 导出网格相互独立
        'quadrature_relative_tolerance': 0.002,
        'max_gap_refinements': 2,
        # 缺失PDF概率按零导流处理，并给出该概率无限导流时的上界；不更改真实接触面积。
        'unresolved_flow_policy': 'zero_conductance_with_upper_bound',
    },
    # 泄漏示例尺寸/压差仅为演示。RUN_LEAKAGE=True 后计算固定接触状态下的流量。
    'leakage': {
        'inlet_pressure_pa': 1e4,
        'outlet_pressure_pa': 0.0,
        'seal_length_m': 1e-3,
        'seal_width_m': 1e-2,
    },
    'mixed_friction': {
        'cylinder_radius_m': None,
        'normal_load_per_width_n_m': None,
        'sliding_velocities_m_s': None,
        'fluid_viscosity_pa_s': None,
    },
}
