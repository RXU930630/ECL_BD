"""Exercise the downloaded ZIP layout, not just the editable installation."""
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import unittest


class SourceEntryTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.project = Path(self.temp.name) / 'downloaded project'
        self.project.mkdir()
        source = Path(__file__).resolve().parents[1]
        for name in ('main.py', 'parameters.py'):
            shutil.copyfile(source / name, self.project / name)
        shutil.copytree(source / 'src', self.project / 'src', ignore=shutil.ignore_patterns('__pycache__'))
        shutil.copytree(source / 'examples', self.project / 'examples')
        path = self.project / 'main.py'
        path.write_text(path.read_text(encoding='utf-8').replace('COMPUTE_GAP_PDF = True', 'COMPUTE_GAP_PDF = False')
                        .replace('SAVE_PLOTS = True', 'SAVE_PLOTS = False')
                        .replace('COMPUTE_PRESSURE_FLOW = True', 'COMPUTE_PRESSURE_FLOW = False')
                        .replace('RUN_LEAKAGE = True', 'RUN_LEAKAGE = False'), encoding='utf-8')

    def test_downloaded_source_runs_from_another_working_directory(self):
        driver = '''
import pathlib, runpy, sys
entry = pathlib.Path(sys.argv[1])
sys.path.insert(0, str(entry.parent))
try:
    runpy.run_path(str(entry), run_name='__main__')
except SystemExit as exc:
    import persson_contact
    assert pathlib.Path(persson_contact.__file__).resolve().is_relative_to(entry.parent / 'src')
    raise
'''
        result = subprocess.run([sys.executable, '-I', '-c', driver, str(self.project / 'main.py')],
                                cwd=self.temp.name, capture_output=True, text=True, timeout=60)
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        self.assertTrue((self.project / 'results' / 'leakage' / 'contact_table.csv').is_file())

    def test_downloaded_spatial_leakage_entry_and_automatic_dependencies(self):
        import json
        entry=self.project/'main.py'
        entry.write_text(entry.read_text(encoding='utf-8').replace('RUN_LEAKAGE = False','RUN_LEAKAGE = True'),encoding='utf-8')
        config=self.project/'parameters.py'
        config.write_text(config.read_text(encoding='utf-8').replace("'number_of_points': 110", "'number_of_points': 3")
                          .replace("[0.0, 1e4, 2e4, 5e4]", "[0.0, 1e4]"),encoding='utf-8')
        driver="import runpy,sys,pathlib; p=pathlib.Path(sys.argv[1]); sys.path.insert(0,str(p.parent)); runpy.run_path(str(p),run_name='__main__')"
        result=subprocess.run([sys.executable,'-I','-c',driver,str(entry)],cwd=self.temp.name,
                              capture_output=True,text=True,timeout=60)
        self.assertEqual(result.returncode,0,result.stdout+result.stderr)
        output=self.project/'results/leakage'
        self.assertTrue((output/'leakage_cases.csv').is_file())
        saved=json.loads((output/'diagnostics.json').read_text())
        self.assertEqual(saved['workflow_status'],'complete')
        self.assertEqual(saved['leakage']['case_count'],2)

    def test_missing_dependency_reports_install_command_without_traceback(self):
        driver = '''
import importlib.abc, pathlib, runpy, sys
class MissingNumpy(importlib.abc.MetaPathFinder):
    def find_spec(self, fullname, path=None, target=None):
        if fullname == 'numpy':
            raise ModuleNotFoundError("No module named 'numpy'", name='numpy')
sys.meta_path.insert(0, MissingNumpy())
entry = pathlib.Path(sys.argv[1])
sys.path.insert(0, str(entry.parent))
runpy.run_path(str(entry), run_name='__main__')
'''
        result = subprocess.run([sys.executable, '-I', '-c', driver, str(self.project / 'main.py')],
                                cwd=self.temp.name, capture_output=True, text=True, timeout=30)
        self.assertEqual(result.returncode, 2, result.stdout + result.stderr)
        self.assertIn('pip install -r requirements.txt', result.stdout)
        self.assertNotIn('Traceback', result.stderr)

    def test_downloaded_shear_switch_and_parameter_presets(self):
        import json
        entry=self.project/'main.py'
        entry.write_text(entry.read_text().replace('RUN_LEAKAGE = False','RUN_LEAKAGE = True')
            .replace('COMPUTE_SHEAR_FLOW = False','COMPUTE_SHEAR_FLOW = True'))
        config=self.project/'parameters.py'
        config.write_text(config.read_text().replace("'geometry': 'oring_hertz',","'geometry': 'square_ring',")
            .replace("'coupling': 'local_unloading',","'coupling': 'fixed_contact',")
            .replace("'compute_equivalent_gas': True,","'compute_equivalent_gas': False,")
            .replace("'number_of_points': 110","'number_of_points': 3")
            .replace('[0.0, 1e4, 2e4, 5e4]','[0.0, 1e4]'))
        driver="import runpy,sys,pathlib; p=pathlib.Path(sys.argv[1]); sys.path.insert(0,str(p.parent)); runpy.run_path(str(p),run_name='__main__')"
        result=subprocess.run([sys.executable,'-I','-c',driver,str(entry)],cwd=self.temp.name,
                              capture_output=True,text=True,timeout=60)
        self.assertEqual(result.returncode,0,result.stdout+result.stderr)
        output=self.project/'results/leakage'
        self.assertTrue((output/'shear_flow_table.csv').exists())
        saved=json.loads((output/'leakage_diagnostics.json').read_text())
        self.assertTrue(saved['motion_transport_enabled'])
        self.assertGreater(saved['cases'][0]['flow_rate_m3_s'],0.)

    def test_downloaded_navier_switch_in_parameters(self):
        import json
        entry=self.project/'main.py'
        entry.write_text(entry.read_text().replace('RUN_LEAKAGE = False','RUN_LEAKAGE = True'))
        config=self.project/'parameters.py'
        config.write_text(config.read_text().replace("'wall_slip': 'no_slip',","'wall_slip': 'navier',")
            .replace("'compute_equivalent_gas': True,","'compute_equivalent_gas': False,")
            .replace("'number_of_points': 110","'number_of_points': 3")
            .replace('[0.0, 1e4, 2e4, 5e4]','[0.0, 1e4]'))
        driver="import runpy,sys,pathlib; p=pathlib.Path(sys.argv[1]); sys.path.insert(0,str(p.parent)); runpy.run_path(str(p),run_name='__main__')"
        result=subprocess.run([sys.executable,'-I','-c',driver,str(entry)],cwd=self.temp.name,
                              capture_output=True,text=True,timeout=60)
        self.assertEqual(result.returncode,0,result.stdout+result.stderr)
        saved=json.loads((self.project/'results/leakage/leakage_diagnostics.json').read_text())
        self.assertEqual(saved['liquid_interface']['wall_slip'],'navier')
        self.assertEqual(saved['liquid_interface']['lower_slip_length_m'],1e-8)
        self.assertGreater(saved['cases'][1]['flow_rate_m3_s'],0)

    def test_downloaded_gas_fluid_switch_in_parameters(self):
        import json
        entry=self.project/'main.py'
        entry.write_text(entry.read_text().replace('RUN_LEAKAGE = False','RUN_LEAKAGE = True'))
        config=self.project/'parameters.py'
        config.write_text(config.read_text().replace("'fluid': 'liquid',","'fluid': 'gas',")
            .replace("'number_of_points': 110","'number_of_points': 3")
            .replace('[0.0, 1e4, 2e4, 5e4]','[0.0, 1e4]'))
        driver="import runpy,sys,pathlib; p=pathlib.Path(sys.argv[1]); sys.path.insert(0,str(p.parent)); runpy.run_path(str(p),run_name='__main__')"
        result=subprocess.run([sys.executable,'-I','-c',driver,str(entry)],cwd=self.temp.name,
                              capture_output=True,text=True,timeout=60)
        self.assertEqual(result.returncode,0,result.stdout+result.stderr)
        output=self.project/'results/leakage'
        saved=json.loads((output/'leakage_diagnostics.json').read_text())
        self.assertEqual(saved['fortran_mode'],5000)
        self.assertEqual(saved['fluid'],'gas')
        self.assertTrue((output/'ballistic_flow_table.csv').exists())
        self.assertGreater(saved['cases'][1]['molecular_flow_rate_per_s'],0.)


if __name__ == '__main__':
    unittest.main()
