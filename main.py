"""日常入口：在这里开启或关闭功能，然后运行 python main.py。"""
from pathlib import Path
import sys

PROJECT_DIRECTORY = Path(__file__).resolve().parent
# ZIP 源码版可直接运行，不要求先执行 pip install -e .。
sys.path.insert(0, str(PROJECT_DIRECTORY / 'src'))

from parameters import PARAMETERS

# 当前已实现的共用 Persson 弹性接触基础。
RUN_CONTACT = True
# 模块6及2+6：独立宏观干接触；开启时关闭泄漏、压力/剪切流动和混合摩擦开关
RUN_MACRO_CONTACT = False
# 固体接触滑动摩擦：自动调用宏观干接触；在parameters中填写实验剪切强度
# 开启时关闭RUN_LEAKAGE、RUN_MIXED_FRICTION、COMPUTE_PRESSURE_FLOW、COMPUTE_SHEAR_FLOW
RUN_CONTACT_FRICTION = False
# 在参数文件中用fluid='liquid'/'gas'选择液体或ballistic气体。
RUN_LEAKAGE = True
RUN_MIXED_FRICTION = False

COMPUTE_PRESSURE_FLOW = True # 压力流动系数；泄漏开启时会自动计算
COMPUTE_SHEAR_FLOW = False   # 剪切流动因子；速度与粗糙侧在parameters.py的shear_flow中
COMPUTE_GAP_PDF = True      # 是否导出微观间隙分布，不控制流动所需的内部计算
SAVE_PLOTS = True           # 输出 PNG，需 matplotlib


def main():
    if not (PROJECT_DIRECTORY / 'src' / 'persson_contact' / '__init__.py').is_file():
        print('缺少 src/persson_contact。请下载并解压完整项目，不要只复制 main.py。')
        return 2
    try:
        from persson_contact.workflow import run_project
        run_project(PARAMETERS, base_directory=Path(__file__).resolve().parent,
                    run_contact=RUN_CONTACT, run_leakage=RUN_LEAKAGE,
                    run_mixed_friction=RUN_MIXED_FRICTION,
                    compute_gap_pdf=COMPUTE_GAP_PDF, save_plots=SAVE_PLOTS,
                    compute_pressure_flow=COMPUTE_PRESSURE_FLOW,
                    compute_shear_flow=COMPUTE_SHEAR_FLOW,
                    run_macro_contact=RUN_MACRO_CONTACT,
                    run_contact_friction=RUN_CONTACT_FRICTION)
    except ModuleNotFoundError as exc:
        if exc.name not in {'numpy', 'scipy', 'matplotlib'}:
            raise
        print(f'当前 Python 缺少依赖 {exc.name}。请在项目目录执行：')
        print(f'"{sys.executable}" -m pip install -r requirements.txt')
        return 2
    except (ValueError, NotImplementedError, FileNotFoundError, FileExistsError) as exc:
        print(f'无法执行：{exc}')
        return 2
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
