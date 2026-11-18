"""Build/extract the real clean ZIP and execute its entry points from elsewhere."""
import hashlib
import importlib.util
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
import zipfile


ROOT=Path(__file__).resolve().parents[1]
spec=importlib.util.spec_from_file_location('local_package_builder',ROOT/'scripts/build_local_package.py')
builder=importlib.util.module_from_spec(spec);spec.loader.exec_module(builder)


class LocalPackageTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.temp=tempfile.TemporaryDirectory();cls.addClassCleanup(cls.temp.cleanup)
        cls.workspace=Path(cls.temp.name);cls.zip=cls.workspace/'candidate.zip'
        builder.build_package(ROOT,cls.zip)
        with zipfile.ZipFile(cls.zip) as z:z.extractall(cls.workspace/'extracted')
        cls.project=next((cls.workspace/'extracted').iterdir())

    def test_manifest_complete_runtime_clean_paths_and_deterministic_build(self):
        manifest=json.loads((self.project/'PACKAGE_MANIFEST.json').read_text())
        self.assertEqual(manifest['version'],'0.2.0rc3')
        for name,data in manifest['files'].items():
            self.assertEqual(hashlib.sha256((self.project/name).read_bytes()).hexdigest(),data['sha256'])
            self.assertFalse(any(part in name.split('/') for part in ('tests','results','.git','__pycache__','.venv')))
            self.assertFalse(name.endswith(('.png','.zip','.pyc')))
        for source in (ROOT/'src/persson_contact').glob('*.py'):
            self.assertEqual((self.project/'src/persson_contact'/source.name).read_bytes(),source.read_bytes())
        self.assertFalse((self.project/'inputs/experiment.csv').exists())
        self.assertEqual(len((self.project/'inputs/experiment_template.csv').read_text().splitlines()),1)
        self.assertTrue((self.project/'inputs/experiment_template.xlsx').is_file())
        self.assertFalse((self.project/'inputs/experiment.xlsx').exists())
        second=self.workspace/'same_again.zip';builder.build_package(ROOT,second)
        self.assertEqual(second.read_bytes(),self.zip.read_bytes())
        with self.assertRaises(FileExistsError):builder.build_package(ROOT,self.zip)

    def test_extracted_installation_and_plotted_fit_demo_from_other_cwd(self):
        for script,extra in [('scripts/check_installation.py',[]),('examples/fit_demo.py',['--output',str(self.workspace/'demo_output')])]:
            result=subprocess.run([sys.executable,str(self.project/script),*extra],cwd=self.workspace,
                capture_output=True,text=True,timeout=60)
            self.assertEqual(result.returncode,0,result.stdout+result.stderr)
        d=json.loads((self.workspace/'demo_output/fit_diagnostics.json').read_text())
        self.assertTrue(d['synthetic_demo'])
        self.assertAlmostEqual(d['fit']['fitted_shear_strength_pa'],3e5)
        self.assertTrue((self.workspace/'demo_output/shear_strength_fit.png').is_file())

    def test_extracted_main_fit_switch(self):
        # Use a second extraction, preserving the immutable candidate tested above.
        with zipfile.ZipFile(self.zip) as z:z.extractall(self.workspace/'main_check')
        project=next((self.workspace/'main_check').iterdir())
        sys.path.insert(0,str(ROOT/'examples'))
        try:
            from fit_demo import write_demo_inputs
            inputs=self.workspace/'main_inputs';inputs.mkdir();experiment,forward=write_demo_inputs(inputs)
        finally:sys.path.pop(0)
        p={'macro_contact':{'line_load_n_m':100.},'contact_friction':{'transverse_width_m':.01},
            'shear_strength_fit':{'experiment_file':str(experiment),'source':'saved_forward',
                'source_parameters':{'saved_forward':{'directories':[str(forward)]}}},
            'output_directory':str(self.workspace/'main_result')}
        (project/'parameters.py').write_text('PARAMETERS = '+repr(p)+'\n')
        entry=project/'main.py';entry.write_text(entry.read_text().replace('RUN_LEAKAGE = True','RUN_LEAKAGE = False')
            .replace('RUN_SHEAR_STRENGTH_FIT = False','RUN_SHEAR_STRENGTH_FIT = True').replace('SAVE_PLOTS = True','SAVE_PLOTS = False'))
        result=subprocess.run([sys.executable,str(entry)],cwd=self.workspace,capture_output=True,text=True,timeout=60)
        self.assertEqual(result.returncode,0,result.stdout+result.stderr)
        d=json.loads((self.workspace/'main_result/fit_summary.json').read_text())
        self.assertAlmostEqual(d['fitted_shear_strength_pa'],3e5)

    def test_silicone_trial_preset_has_real_dimensions_and_fit_ready_entry(self):
        import runpy
        target=self.workspace/'silicone.zip'
        builder.build_package(ROOT,target,preset='silicone_cylinder')
        with zipfile.ZipFile(target) as z:z.extractall(self.workspace/'silicone')
        project=next((self.workspace/'silicone').iterdir())
        p=runpy.run_path(str(project/'parameters.py'))['PARAMETERS']
        from persson_contact.macro_output import prepare_macro
        args=prepare_macro(p,project)
        self.assertEqual(args[2].line_load_n_m,1000.)
        self.assertEqual(args[3].edges_m[-1],.015)
        self.assertEqual(p['material']['young_modulus_pa'],1.8e6)
        self.assertEqual(p['interface']['rheology_parameters']['newtonian']['viscosity_pa_s'],1.)
        self.assertIn('RUN_SHEAR_STRENGTH_FIT = True',(project/'main.py').read_text())
        self.assertIn('RUN_LEAKAGE = False',(project/'main.py').read_text())
        self.assertFalse((project/'inputs/experiment.csv').exists())


if __name__=='__main__':unittest.main()
