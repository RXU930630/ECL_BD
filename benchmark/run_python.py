"""Run the paired forward benchmark without changing main.py or its parameters."""
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'src'))
from python_parameters import PARAMETERS
from persson_contact.workflow import run_project

if __name__ == '__main__':
    run_project(PARAMETERS, base_directory=ROOT, run_contact=True,
                run_mixed_friction=True, run_leakage=False, save_plots=True)
