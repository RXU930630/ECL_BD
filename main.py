"""日常入口：在这里开启或关闭功能，然后运行 python main.py。"""
from pathlib import Path
import sys

PROJECT_DIRECTORY = Path(__file__).resolve().parent
# ZIP 源码版可直接运行，不要求先执行 pip install -e .。
sys.path.insert(0, str(PROJECT_DIRECTORY / 'src'))

from parameters import PARAMETERS

# 当前已实现的共用 Persson 弹性接触基础。
RUN_CONTACT = True
# 两个应用接口已预留，求解器尚未实现；开启后会明确提示。
RUN_LEAKAGE = False
RUN_MIXED_FRICTION = False

COMPUTE_GAP_PDF = True      # 实验性微观间隙分布
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
                    compute_gap_pdf=COMPUTE_GAP_PDF, save_plots=SAVE_PLOTS)
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
