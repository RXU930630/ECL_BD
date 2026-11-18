"""Actual diagnostics and fitting under a simulated Windows ANSI text default."""
from contextlib import contextmanager
from copy import deepcopy
import importlib.util
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

from openpyxl import load_workbook
import test_mixed_checks as mixed_fixture
import test_shear_fit as fit_fixture
from persson_contact.mixed_checks import run_checks,check_settings,write_mixed_workbook
from persson_contact.workflow import run_project


ROOT=Path(__file__).resolve().parents[1]


@contextmanager
def windows_ansi_defaults():
    read,write=Path.read_text,Path.write_text
    def ansi_write(path,data,encoding=None,errors=None,newline=None):
        return write(path,data,encoding=encoding or 'cp1252',errors=errors,newline=newline)
    def ansi_read(path,encoding=None,errors=None):
        return read(path,encoding=encoding or 'cp1252',errors=errors)
    with patch.object(Path,'write_text',ansi_write),patch.object(Path,'read_text',ansi_read):yield


def load_script(name):
    spec=importlib.util.spec_from_file_location(name,ROOT/'scripts'/f'{name}.py')
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    return module


class WindowsOutputTests(unittest.TestCase):
    def test_numerical_check_chinese_warning_roundtrips_and_excel_names_stay_fixed(self):
        with tempfile.TemporaryDirectory() as temporary:
            root=Path(temporary)/'中文工况';root.mkdir()
            p,s=mixed_fixture.NumericalAuditTests().setup_audit(root)
            with windows_ansi_defaults():
                result=run_checks(p,root,s,{},check_settings({'enabled':True}),
                                 lambda *a,**k:{'cases':[mixed_fixture.case()]},root)
                s.update(numerical_validation=result,solid_friction={})
                write_mixed_workbook(root/'mixed_results.xlsx',s)
            saved=json.loads((root/'mixed_numerical_checks.json').read_text(encoding='utf-8'))
            self.assertEqual(saved['status'],'passed_at_selected_cases')
            self.assertIn('有限胶块',saved['warnings'][0])
            book=load_workbook(root/'mixed_results.xlsx',read_only=True)
            self.assertEqual(book.sheetnames,['计算概况','正向摩擦与面积','数值核查差异','核查执行情况'])
            book.close()

    def test_saved_fit_chinese_json_keeps_rc5_filenames_and_sheets(self):
        with tempfile.TemporaryDirectory() as temporary:
            p=fit_fixture.FitWorkflowTests().configure(temporary)
            p['description']='硅油实验：原参数保留'
            with windows_ansi_defaults():
                result=run_project(p,base_directory=temporary,run_contact=False,
                                   run_shear_strength_fit=True,save_plots=False)
            out=Path(temporary)/'fit'
            self.assertEqual({p.name for p in out.iterdir()},
                {'fit_results.xlsx','fit_cases.csv','fit_summary.json','fit_diagnostics.json',
                 'fitted_contact_friction.json','forward_response.csv','experiment_used.csv'})
            saved=json.loads((out/'fit_diagnostics.json').read_text(encoding='utf-8'))
            self.assertEqual(saved['parameters_snapshot']['description'],p['description'])
            self.assertAlmostEqual(result['fit']['fitted_shear_strength_pa'],4e5)
            book=load_workbook(out/'fit_results.xlsx',read_only=True)
            self.assertEqual(book.sheetnames,['拟合概况','各载荷独立诊断','各实验点','正向响应','所用实验'])
            book.close()

    def test_failure_report_preserves_original_solver_error_instead_of_masking_with_encoding(self):
        with tempfile.TemporaryDirectory() as temporary:
            p=fit_fixture.FitWorkflowTests().configure(temporary)
            p['psd']['path']=str(ROOT/'examples/data/roughness_psd.dat')
            p['shear_strength_fit']['source']='simulate'
            with windows_ansi_defaults(),patch('persson_contact.shear_fit_output.run_mixed_project',
                    side_effect=ValueError('原始求解错误')):
                with self.assertRaisesRegex(ValueError,'原始求解错误'):
                    run_project(p,base_directory=temporary,run_contact=False,run_shear_strength_fit=True,save_plots=False)
            saved=json.loads((Path(temporary)/'fit/fit_diagnostics.json').read_text(encoding='utf-8'))
            self.assertEqual(saved['error'],'原始求解错误')


class SafeUpdateRecoveryTests(unittest.TestCase):
    def test_update_preserves_user_settings_inputs_and_results_with_backup(self):
        module=load_script('update_runtime')
        with tempfile.TemporaryDirectory() as temporary:
            root=Path(temporary);old=root/'old';source=root/'new'
            for p in (old/'src/persson_contact',source/'src/persson_contact',old/'inputs',old/'results'):p.mkdir(parents=True)
            for name in ('main.py','parameters.py','inputs/实验.xlsx','results/原结果.txt'):
                (old/name).write_text('原有内容',encoding='utf-8')
            (old/'src/persson_contact/__init__.py').write_text('old',encoding='utf-8')
            (source/'src/persson_contact/__init__.py').write_text('new',encoding='utf-8')
            with windows_ansi_defaults():module.update_runtime(source,old)
            self.assertEqual((old/'src/persson_contact/__init__.py').read_text(),'new')
            self.assertEqual((old/'runtime_backup/__init__.py').read_text(),'old')
            for name in ('main.py','parameters.py','inputs/实验.xlsx','results/原结果.txt'):
                self.assertEqual((old/name).read_text(encoding='utf-8'),'原有内容')
            with self.assertRaises(FileExistsError):module.update_runtime(source,old)

    def test_recover_failed_fit_in_place_without_forward_recalculation_or_setting_changes(self):
        module=load_script('recover_fit')
        with tempfile.TemporaryDirectory() as temporary:
            project=Path(temporary);p=fit_fixture.FitWorkflowTests().configure(temporary)
            (project/'main.py').write_text('SAVE_PLOTS = False\n',encoding='utf-8')
            p['shear_strength_fit']['source']='simulate';p['description']='编码失败后的拟合'
            (project/'parameters.py').write_text('PARAMETERS = '+repr(p),encoding='utf-8')
            before=(project/'parameters.py').read_bytes()
            out=project/'fit';out.mkdir();(out/'forward').mkdir()
            (project/'forward_0').rename(out/'forward/load_000')
            failure={'workflow_status':'forward_failed','error':'charmap codec'}
            (out/'fit_diagnostics.json').write_text(json.dumps(failure),encoding='utf-8')
            with windows_ansi_defaults(),patch('persson_contact.shear_fit_output.run_mixed_project',
                    side_effect=AssertionError('Saved states must not be recalculated')):
                result=module.recover_fit(project)
            self.assertAlmostEqual(result['fit']['fitted_shear_strength_pa'],4e5)
            self.assertEqual((project/'parameters.py').read_bytes(),before)
            self.assertEqual(json.loads((out/'fit_diagnostics.json').read_text(encoding='utf-8'))['workflow_status'],'complete')
            backup=project/'runtime_backup/fit_before_recovery/fit_diagnostics.json'
            self.assertEqual(json.loads(backup.read_text(encoding='utf-8')),failure)
            with self.assertRaises(ValueError):module.recover_fit(project)


if __name__=='__main__':unittest.main()
