# 本地运行与候选版本范围

本包版本0.2.0rc3，推荐Python3.12；启动方式与安装命令见根目录README.md。所有计算都在本地完成。src中的源代码、main.py、parameters.py和必要输入必须保留在同一解压目录结构内。有限圆柱曲面的尺寸、自动载荷换算和物理边缘见[几何说明](finite_geometry.md)。

## 任务选择

日常只编辑main.py和parameters.py，不需要编辑src。安装后先运行scripts/check_installation.py，再用examples/fit_demo.py检查绘图和拟合；示例会在运行时生成合成数据，不是用户实验。

| 开关 | 作用 |
| --- | --- |
| RUN_CONTACT | 微观Persson接触与平均间隙基础 |
| RUN_LEAKAGE | 液体或ballistic气体泄漏，使用leakage中的几何与边界 |
| RUN_MACRO_CONTACT | 独立圆柱—平面宏观干接触，smooth或persson |
| RUN_CONTACT_FRICTION | 给定实验tau后的独立干摩擦 |
| RUN_MIXED_FRICTION | 给定tau后的流体—接触—形变耦合与摩擦 |
| RUN_SHEAR_STRENGTH_FIT | 独立实验反演；内部调用混合润滑或复用旧解 |
| COMPUTE_PRESSURE_FLOW、COMPUTE_SHEAR_FLOW | 独立流量因子任务；混合润滑和拟合内部自动计算所需因子 |
| COMPUTE_GAP_PDF | 控制微观间隙分布导出，不关闭求解所需的内部计算 |
| SAVE_PLOTS | 生成PNG图 |

源码保留既有默认泄漏示例，尚未替换为用户的实验配置。使用拟合前请按根目录README切换开关。多个独立应用同时打开会明确报错，避免把不同任务的几何或载荷混用。

## 物理输入

材料、PSD和恒黏度流体在共用分组。参数全部为SI单位：m、Pa、Pa·s、m/s、N或N/m。宏观圆柱曲率半径不是环形密封用于2πR宽度换算的整体半径；总力F_N与线载荷Wprime通过横向宽度B联系。

macro_contact用于圆柱接触、混合润滑及拟合；leakage.geometry用于已有的Hertz参考压力、方环、FEM给定曲线或高斯压力泄漏。两者是不同工作流。PSD文件继续使用仓库中的examples/data/roughness_psd.dat；用户有新实测PSD时通过parameters.py指定路径。

混合润滑当前为一维稳态、线弹性、各向同性单侧粗糙、恒黏度牛顿液体、无滑移、充分供油入口、压力—充满比例守恒空化。不能求完全阻断的封闭液袋、有限供油、瞬态起停、橡胶黏弹滞后、毛细力或热反馈。流体对两壁的摩擦力分开输出，主摩擦系数针对运动的下平面。

界面剪切强度恒定时，固体力由真实接触面积乘tau得到；log10速度高斯和自定义函数可用于普通正向摩擦任务，当前自动反演仅支持常数。非牛顿流变和黏度压温依赖只有接口预留，未实现计算。

## 常见问题

- 缺numpy/scipy/matplotlib/openpyxl：用运行程序的同一个Python执行`-m pip install -r requirements.txt`。
- 输出目录非空：换一个新的OUTPUT_DIRECTORY，避免覆盖前次结果。
- 实验表没有数据：请用Excel打开inputs/experiment_template.xlsx填写后另存为experiment.xlsx，空白模板没有测量值，必须填入恒速段测量值。
- reference tau与实验tau：1MPa只是拟合参考，不是最后结果；普通正向计算则需要给定真实或明确假设的tau。
- normal_force_n报缺宽度：有限曲面模式填写geometry_parameters.finite_cylindrical_pad.axial_length_m；旧圆柱模式在contact_friction.transverse_width_m填写横向宽度，单位m。
- 缓存缺少对应速度／载荷：选择simulate按实验点计算；不自动插值。
- 模型表域或计算域不足：按报错检查压力表域及数值分辨率。有限曲面保持真实尺寸，不为改善收敛而扩大胶块；旧圆柱远场截断才检查域宽收敛。
- 零导流／封闭液袋：这是当前模型适用范围限制，不通过增加spc地板来强行求解。
- 矩阵病态或不收敛：保留报错和参数；检查工况、查表范围及连通分支，不把失败结果用于反演。

## 清洁包与可复现

不包含历史结果或测试夹具。PACKAGE_MANIFEST.json记录每个运行文件的SHA256；纯净ZIP旁的.sha256用于校验传输。脚本build_local_package.py使用明确的文件清单，即使项目中已产生results、用户实验或缓存，再次打包也不会把它们带入。

数值检查与历史基准保留在开发仓库中。本地试运行时请先核对模型假设，再保存自己的输入、正向结果、拟合诊断和程序版本；候选版本名称本身不构成实验精度保证。
