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
                        .replace('SAVE_PLOTS = True', 'SAVE_PLOTS = False'), encoding='utf-8')

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
        self.assertTrue((self.project / 'results' / 'contact' / 'contact_table.csv').is_file())

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


if __name__ == '__main__':
    unittest.main()
