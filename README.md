本包已预填本次硅油实验几何、材料与黏度，并开启独立拟合。请用Excel填写 inputs/experiment.xlsx 后运行；该文件不含你的实测摩擦数据。默认100N仅用于CSV未给载荷的行。详见 docs/finite_geometry.md。

# ECL_BD 0.2.0rc3 — 本地试运行版（候选稳定版）

这是完整的Python运行源码包，供首次本地计算和实验拟合使用。包含现有接触、泄漏、宏观形变、混合润滑和界面剪切强度拟合功能；不包含历史计算结果、原Fortran包、开发测试数据、缓存或实验数据。安装依赖需要联网；计算与拟合在本地运行，不需要云端服务。

`rc3`新增直接读取Excel实验工作簿与Excel拟合结果，继续支持有限圆柱曲面尺寸。已有数值检查，尚待用户真实实验的首次本地试运行，仍为候选版。几何说明见[几何与边界](docs/finite_geometry.md)，使用范围见[本地运行说明](docs/local_run.md)。

## 第一次安装

推荐Python 3.12，支持Python 3.11及以上。先解压完整目录，不要只复制main.py。打开终端并进入该目录。

Windows：

```powershell
py -3.12 -m venv .venv
.venv\Scripts\python.exe -m pip install -r requirements.txt
.venv\Scripts\python.exe scripts\check_installation.py
.venv\Scripts\python.exe examples\fit_demo.py
```

macOS/Linux：

```bash
python3 -m venv .venv
.venv/bin/python -m pip install -r requirements.txt
.venv/bin/python scripts/check_installation.py
.venv/bin/python examples/fit_demo.py
```

若只有其他受支持版本的Python，使用对应解释器创建虚拟环境即可。fit_demo在运行时生成合成数据，已知剪切强度0.3MPa；它用于检查拟合流程和出图，不计算真实混合润滑，也不是历史测试结果。输出目录results/fit_demo已存在且非空时，请用`--output results/fit_demo_2`指定新目录。

## 使用自己的实验

1. 用Excel打开`inputs/experiment_template.xlsx`，在`experiment`工作表填写数据，另存为`inputs/experiment.xlsx`。第1行英文表头保持不变；中文说明在表头批注与`填写说明`工作表中。每个速度/载荷组合填一行，速度单位m/s、总力单位N。无需另存为CSV。
2. 在parameters.py填写实际材料、PSD、几何、总法向载荷或线载荷、硅油黏度与供油边界。默认物理参数只是示例，不是用户实验参数。
3. 在main.py设置：

```python
RUN_LEAKAGE = False
RUN_MACRO_CONTACT = False
RUN_CONTACT_FRICTION = False
RUN_MIXED_FRICTION = False
RUN_SHEAR_STRENGTH_FIT = True
```

4. parameters.py中选择`shear_strength_fit.source='simulate'`，使用新的`OUTPUT_DIRECTORY='results/my_fit'`，运行`python main.py`（使用虚拟环境中的Python）。独立拟合自动计算正向解，RUN_CONTACT和COMPUTE开关无需额外打开。
5. 用Excel打开输出目录中的`fit_results.xlsx`，在`拟合概况`查看剪切强度Pa/MPa，在`各实验点`查看面积、固体/流体分项及残差；图片为`shear_strength_fit.png`。原CSV/JSON输出继续保留。

无需在contact_friction中预先填写未知剪切强度：拟合任务会临时使用参考值1MPa，并从真实面积和流体摩擦求出结果。材料和流体参数不会自动反演。

拟合后，若要做普通混合润滑预测，将得到的剪切强度填回`contact_friction.model_parameters.constant_shear_strength.shear_strength_pa`，关闭拟合开关，开启RUN_MIXED_FRICTION。程序不会自动覆盖你的参数文件。

## 反复拟合不必重算

第一次计算会按载荷在输出目录中保存`forward/load_000`、`forward/load_001`等完整混合润滑结果。调整拟合范围、权重、实验值或留出点时，改为：

```python
'source': 'saved_forward',
'source_parameters': {
    'simulate': {},
    'saved_forward': {'directories': ['results/my_fit/forward/load_000']},
},
```

多载荷时把所需load目录全部列出。输出目录换成新的空目录。复用时不重跑流固耦合，且以保存的物理工况为准；改动当前材料、PSD、硅油或剪切下限不会更新旧结果，须重新选择simulate计算。

Excel和CSV使用相同的列名与SI单位。旧CSV仍可通过`experiment_file`选择；不用改文件扩展名来转换格式。全部列、拟合约定、物理限制和常见错误见[剪切强度拟合说明](docs/shear_strength_fit.md)。

## 包内容与版本

- main.py：选择任务。
- parameters.py：工况与数值设置，注释中列出当前可选模型。
- src/persson_contact：全部运行模块。
- inputs：仅空白实验模板，不含用户实验。
- examples：必要输入和运行时生成数据的快速演示。
- docs：本地使用与拟合说明，不包含历史验证结果。
- scripts：安装自检与纯净包构建脚本。
- PACKAGE_MANIFEST.json：版本、文件SHA256和来源基准提交；不包含它自身的哈希。

非空计算目录不会覆盖。包内没有已计算的results目录；所有结果由本地运行生成。原开发仓库仍保留历史验证材料，与此运行包分开。
