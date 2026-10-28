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
    'pressures_pa': [1.0e4, 2.0e4, 5.0e4, 1.0e5, 2.0e5],  # 名义干接触压力 Pa
    'options': {
        'gamma': 0.4,       # 弹性能量修正，不是摩擦系数
        'n_q': 256,         # 波数/尺度分辨率，整数 >= 32
        'pdf_afact': 1.0,   # 分布展宽修正
        'n_gap': 600,       # 微观间隙输出网格，整数 >= 32
    },
    'output_directory': 'results/contact',  # 非空时拒绝覆盖；新工况请换目录
    # 以下仅为未来接口，None 表示尚未给定；当前不参与接触计算。
    'leakage': {
        'fluid_viscosity_pa_s': None,
        'inlet_pressure_pa': None,
        'outlet_pressure_pa': None,
        'seal_length_m': None,
        'seal_width_m': None,
    },
    'mixed_friction': {
        'cylinder_radius_m': None,
        'normal_load_per_width_n_m': None,
        'sliding_velocities_m_s': None,
        'fluid_viscosity_pa_s': None,
    },
}
