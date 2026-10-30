"""日常配置：先改工况，必要时再改数值设置。全部使用SI单位。

当前为Hertz示例；功能开关在 main.py。
路径相对于项目目录；末尾自动组装程序读取的 PARAMETERS，无需手改。
"""

# 一、工况：材料、粗糙度、流体、几何与压力
CASE = {
    'description': 'Hertz example.',
    'material': {
        'model': 'elastic',
        'young_modulus_pa': 3.0e6,
        'poisson_ratio': 0.5,
    },
    'psd': {
        'path': 'examples/data/roughness_psd.dat',
        'scale_c': 1.0,
        'tail_is_extrapolated': True,
    },
    # 宏观干接触：main.py开启RUN_MACRO_CONTACT；与当前泄漏任务分开运行
    'macro_contact': {
        'geometry': 'cylinder_plane',
        'geometry_parameters': {
            'cylinder_plane': {
                'radius_m': 5e-3,     # 接触截面曲率半径；不是塞子的整体周向半径
            },
        },
        'elastic_model': 'half_space_plane_strain', # 一维线接触、平面应变、线弹性半空间
        # upper_elastic_lower_rigid：上部弹性圆柱/下部刚性平面，material为上部材料
        # relative_only：material给合并E*，只确定相对变形；显示的上轮廓为等效轮廓
        'body_representation': 'upper_elastic_lower_rigid',
        'contact_model': 'persson',    # persson：模块2+6粗糙干接触；smooth：光滑Hertz检验
        'line_load_n_m': 100.,         # N/m线载荷；总力需除以横向宽度，不能直接填N
        'domain_half_width_m': 1.2e-3, # 计算域[-L,L]；检查边界压力和域加宽收敛
        'reference_length_m': None,   # 2D位移参考长度；None取全域长度，只影响绝对位移与刚体偏移
    },
    'leakage': {
        'mode': 'spatial',
        # 流体：liquid为液体；gas为独立ballistic气体
        'fluid': 'liquid',
        'compute_equivalent_gas': True, # 仅无滑移liquid：附加气体估算；非零滑移时设False

        # 几何选择器：只选一个名称；下面保留全部备选参数，只读取选中的一块
        # oring_hertz         O-ring的Hertz接触压力分布
        # square_ring         给定接触长度、均匀接触压力的方环
        # prescribed_profile  导入FEM或其他接触压力曲线
        # gaussian            给定高斯接触压力分布
        'geometry': 'oring_hertz',
        'geometry_parameters': {

            # O-ring / Hertz：最大压力p_max；由截面直径与E*推算接触长度
            'oring_hertz': {
                'cross_section_diameter_m': 3e-3, # 截面直径，m；不是密封环整体直径
                'peak_pressure_pa': 1e5,          # p_max：峰值接触压力，Pa
                # 若直接给定接触长度，用下面一项替换cross_section_diameter_m
                # 'length_m': 1e-3,               # 流向接触长度，m；两者不能同时给
            },

            # Square ring：当前为均匀接触压力模型，不从截面尺寸计算变形
            'square_ring': {
                'length_m': 1e-3,                 # 沿泄漏方向的接触长度，m
                'reference_contact_pressure_pa': 1e5, # p_0：参考接触压力，Pa
            },

            # FEM / 给定曲线：文件路径相对于项目目录
            'prescribed_profile': {
                'pressure_file': 'examples/data/seal_pressure_profile.csv',
                'file_format': 'csv',             # CSV压力曲线
                # csv：表头x_m,contact_pressure_pa；单位m、Pa
                # 文件坐标按原值读取，不重新排列或自动补零
                'interpolation': 'pchip',         # pchip保形插值 / linear线性插值
            },

            # Gaussian：独立给定计算域和压力分布宽度
            'gaussian': {
                'length_m': 1e-3,                 # 流向计算域总长度，m
                'peak_pressure_pa': 1e5,          # 中心峰值接触压力，Pa
                'sigma_m': 2e-4,                  # 高斯分布标准差，m
            },
        },

        # 所有几何共用：横向宽度、流体压力与接触处理
        'seal_width_m': 1e-2,          # 横向总宽度，m；环形密封填2*pi*R
        'inlet_pressures_pa': [0.0, 1e4, 2e4, 5e4], # Pa表压；小坐标端入口
        'outlet_pressure_pa': 0.0,
        'coupling': 'local_unloading', # 局部压力卸载；fixed_contact为固定接触
        'conductance_model': 'effective_medium', # 有效介质；critical_junction为临界通道
        # gas目前仅支持effective_medium；临界通道可另给critical_alpha，默认1.0
        'percolation_policy': 'legacy_spc_plateau', # 有限导流平台；physical_zero保留阻断
        'spc': 0.8,                    # 仅平台规则使用，会影响物理结果
        'gas': {                      # gas求解或liquid附加估算时使用
            'temperature_c': 25.0,
            'mean_molecular_speed_m_s': 471.0, # 平均分子速率，m/s；不随温度自动更新
            'viscosity_pa_s': 1.82e-5, # 估算及自由程诊断，不控制ballistic通量
            'ambient_pressure_pa': 1e5, # 表压零点对应的绝对压力
            'reference_pressure_pa': 1e5, # 报告气体体积时的参考绝对压力
        },
    },
    # 剪切流动：main.py开启COMPUTE_SHEAR_FLOW后读取；首版限牛顿、无滑移、单侧粗糙
    'shear_flow': {
        'model': 'persson_interpolation', # 当前唯一实现的各向同性剪切因子模型
        'rough_surface': 'lower',         # lower下壁粗糙 / upper上壁粗糙；另一侧光滑
        'lower_wall_velocity_m_s': 0.0,   # 下壁沿+x方向速度，m/s，允许负值
        'upper_wall_velocity_m_s': 0.01,  # 上壁沿+x方向速度，m/s，允许负值
        'pressure_gradient_pa_m': 0.0,    # 仅用于均匀接触状态表的测试梯度，Pa/m
        'pressure_floor_pa': 0.0,         # 运动空间输运的压力检查下限，Pa表压；低于则报错，不裁剪
        # 同时计算运动壁面空间流动时：coupling='fixed_contact'，附加气体估算设False
        # 这一步不求空化、摩擦、载荷分担或膜厚；现有压力卸载仍可独立使用
    },
    'interface': {                     # 仅liquid使用；两个选择器互相独立
        # 流变选择器：目前实现newtonian恒黏度；压温依赖及剪切稀化留待扩展
        'rheology': 'newtonian',
        'rheology_parameters': {
            'newtonian': {
                'viscosity_pa_s': 0.1, # 动力黏度，Pa·s
            },
        },

        # 壁面选择器：no_slip无滑移 / navier两侧独立Navier滑移
        'wall_slip': 'no_slip',
        'slip_parameters': {
            'no_slip': {},            # 无滑移，无额外参数
            'navier': {
                'lower_slip_length_m': 1e-8, # 下壁z=0的滑移长度，m，允许0
                'upper_slip_length_m': 1e-8, # 上壁z=h的滑移长度，m，允许0
                # z是膜厚方向；与流向的小坐标端入口不是同一个方向
            },
        },
    },
}

# 二、模型与数值设置：控制接触模型、积分精度与空间求解
NUMERICS = {
    'macro_contact': {
        'solver': {
            'cells': 257,             # 均匀压力单元数；通常129/257/513检查收敛
            'relative_tolerance': 1e-9,
            'max_iterations': 80,
            'edge_pressure_ratio_tolerance': 1e-4, # 边界压力/峰值过大时报错，需增大计算域
        },
        'contact_table': {            # 仅persson读取；远场极小压力保留正值，不加压力地板
            'minimum_pressure_pa': 1e-100,
            'maximum_pressure_pa': 2e6, # 必须覆盖局部峰值；范围不足明确报错，不外推
            'points': 1201,           # 独立于泄漏数据库；远场尾部与工作区分配不同密度
        },
    },
    'contact_pressure': {             # 接触数据库，不是流体入口压力
        'reference_pa': 1.0e5,         # 参考接触压力，Pa；扫描包含该点
        'upper_factor': 2.0,           # 扫描上限/参考压力；必要时为泄漏自动扩展
        'number_of_points': 110,       # 基础压力节点数，近似对数分布
    },
    'options': {
        'gamma': 0.4,                 # 弹性能量修正
        'pdf_afact': 1.0,             # 分布展宽修正
        'critical_relative_contact_area': 0.42,
        'critical_area_tolerance': 0.005, # 面积分数的绝对容差
        'n_q': 256,                   # 平均/临界间隙积分分辨率
        # conservative：准确谱积分、面积差守恒的尺度积分（推荐）
        # fortran：保留原右矩形和中心差分规则，供历史复算
        'pdf_scale_method': 'conservative',
        'pdf_log_step': 0.01,         # PDF尺度步长；需独立减半检查，导流分箱收敛不能代替它
        'pdf_max_relative_contact_area': 0.999,
        'pdf_min_relative_contact_area': 0.0001,
        'pdf_min_u1_m': 1e-12,
        'n_gap': 600,                 # 仅控制导出的间隙网格
    },
    'pressure_flow': {
        'closure': 'fortran_isotropic_bruggeman',
        'gap_points': 2049,           # 内部导流积分网格，独立于n_gap
        'quadrature_relative_tolerance': 0.002,
        'max_gap_refinements': 2,
        'unresolved_flow_policy': 'zero_conductance_with_upper_bound',
    },
    'leakage': {
        'database_interpolation': 'log_pchip', # 对数保形插值；linear为线性插值
        'compute_upper_bound': True,   # 遗漏PDF概率的导流界限
        'solver': {
            'nodes': 201,             # 初始空间网格点数
            'output_points': 1001,    # 导出的压力曲线点数
            'tolerance': 1e-5,        # 空间求解收敛容差
            'max_nodes': 20001,       # 自适应空间网格点数上限
        },
    },
}

# 三、输出：切换工况或流体后请换目录，已有非空目录不会覆盖
OUTPUT_DIRECTORY = 'results/leakage'

# 自动组装，兼容原有程序入口；上面两块的leakage字段合并
PARAMETERS = {
    **CASE,
    **NUMERICS,
    'leakage': {**CASE['leakage'], **NUMERICS['leakage']},
    'macro_contact': {**CASE['macro_contact'], **NUMERICS['macro_contact']},
    'output_directory': OUTPUT_DIRECTORY,
}
