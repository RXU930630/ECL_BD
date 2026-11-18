"""Real workbook input, formula values, validation and fitted Excel output."""
from copy import deepcopy
import csv
from pathlib import Path
import tempfile
import unittest
from zipfile import ZipFile,ZIP_DEFLATED
from xml.etree import ElementTree as ET

from openpyxl import Workbook,load_workbook
from openpyxl.styles import Font

from persson_contact.experiment_excel import HEADERS,write_experiment_workbook
from persson_contact.shear_fit import read_experiments
from persson_contact.workflow import run_project
from test_shear_fit import make_forward


class ExcelExperimentTests(unittest.TestCase):
    def setUp(self):
        self.temp=tempfile.TemporaryDirectory();self.addCleanup(self.temp.cleanup)
        self.root=Path(self.temp.name);self.path=self.root/'experiment.xlsx'

    def read(self):return read_experiments(self.path,default_line_load_n_m=1000.,transverse_width_m=.1)

    def raw(self,header,rows):
        book=Workbook();sheet=book.active;sheet.append(header)
        for row in rows:sheet.append(row)
        book.save(self.path);book.close()

    def test_excel_and_csv_identical_values_flags_and_force_conversion(self):
        data=[dict(zip(HEADERS,[1e-6,.338072,None,40,None,None,True])),
              dict(zip(HEADERS,[.01,.029435,None,100,.01,.002,False]))]
        write_experiment_workbook(self.path,data)
        path=self.root/'same.csv'
        with path.open('w',newline='') as stream:
            writer=csv.DictWriter(stream,fieldnames=HEADERS);writer.writeheader();writer.writerows(data)
        excel=self.read();legacy=read_experiments(path,default_line_load_n_m=1000.,transverse_width_m=.1)
        self.assertEqual(excel,legacy);self.assertEqual(excel[0].line_load_n_m,400.)
        self.assertFalse(excel[1].use_for_fit)

    def test_named_sheet_selected_even_if_help_sheet_first_and_active(self):
        write_experiment_workbook(self.path,[{'speed_m_s':1e-6,'friction_coefficient':.2}])
        book=load_workbook(self.path);book.move_sheet('填写说明',offset=-1);book.active=0;book.save(self.path);book.close()
        self.assertEqual(self.read()[0].speed_m_s,1e-6)

    def test_first_sheet_fallback_blanks_and_formatted_extra_columns(self):
        self.raw(['speed_m_s','friction_coefficient'],[[.001,.2],[],[.01,.1]])
        book=load_workbook(self.path);book.active['Z100'].font=Font(bold=True);book.save(self.path);book.close()
        rows=self.read();self.assertEqual([r.row_number for r in rows],[2,4])

    def test_template_empty_and_existing_destination_rejected(self):
        write_experiment_workbook(self.path)
        with self.assertRaisesRegex(ValueError,'没有数据'):self.read()
        with self.assertRaises(FileExistsError):write_experiment_workbook(self.path)

    def test_invalid_headers_and_extra_nonempty_columns(self):
        for header,row in [(['speed_m_s','speed_m_s','friction_coefficient'],[1,1,.1]),
            (['speed_m_s','friction_coefficient','notes'],[1,.1,'text']),
            (['speed_m_s'],[1]),(['speed_m_s','friction_coefficient'],[1,.1,'extra'])]:
            self.raw(header,[row])
            with self.subTest(header=header,row=row),self.assertRaises(ValueError):self.read()

    def test_formula_without_saved_value_reports_cell(self):
        self.raw(['speed_m_s','friction_coefficient'],[['=1/1000',.2]])
        with self.assertRaisesRegex(ValueError,'A2.*公式'):self.read()

    def test_saved_formula_value_is_read_without_evaluation(self):
        self.raw(['speed_m_s','friction_coefficient'],[['=1/1000',.2]])
        with ZipFile(self.path) as archive:parts={name:archive.read(name) for name in archive.namelist()}
        ns={'s':'http://schemas.openxmlformats.org/spreadsheetml/2006/main'}
        xml=ET.fromstring(parts['xl/worksheets/sheet1.xml'])
        cell=xml.find(".//s:c[@r='A2']/s:v",ns);self.assertIsNotNone(cell);cell.text='0.001'
        parts['xl/worksheets/sheet1.xml']=ET.tostring(xml,encoding='utf-8')
        with ZipFile(self.path,'w',compression=ZIP_DEFLATED) as archive:
            for name,data in parts.items():archive.writestr(name,data)
        self.assertEqual(self.read()[0].speed_m_s,.001)

    def test_excel_error_cells_and_legacy_xls_have_clear_errors(self):
        self.raw(['speed_m_s','friction_coefficient'],[[.001,'#DIV/0!']])
        with self.assertRaisesRegex(ValueError,'B2.*错误'):self.read()
        with self.assertRaisesRegex(ValueError,'另存为.xlsx'):
            read_experiments(self.root/'old.xls',default_line_load_n_m=100.)

    def test_renamed_csv_is_not_treated_as_a_workbook(self):
        self.path.write_text('speed_m_s,friction_coefficient\n.001,.2\n')
        with self.assertRaisesRegex(ValueError,'不能只改扩展名'):self.read()

    def test_invalid_numeric_and_boolean_fields_preserve_excel_row(self):
        for row in [[0,.1,40,True],[.001,.1,-40,True],[.001,.1,40,'maybe']]:
            self.raw(['speed_m_s','friction_coefficient','normal_force_n','use_for_fit'],[row])
            with self.subTest(row=row),self.assertRaisesRegex(ValueError,'第2行'):self.read()

    def test_native_workbook_through_fitting_task_and_excel_results(self):
        from parameters import PARAMETERS
        p=deepcopy(PARAMETERS);p['contact_friction']['transverse_width_m']=.1
        cases=make_forward(self.root/'forward',load=1000.)
        known=4e5;rows=[]
        for case in cases:
            rows.append({'speed_m_s':case['lower_wall_velocity_m_s'],'normal_force_n':100.,
                'friction_coefficient':-case['lower_fluid_force_n_m']/1000+known*case['real_contact_area_per_width_m']/1000,
                'use_for_fit':True})
        write_experiment_workbook(self.path,rows)
        p['shear_strength_fit'].update(experiment_file=str(self.path),source='saved_forward')
        p['shear_strength_fit']['source_parameters']['saved_forward']['directories']=[str(self.root/'forward')]
        p['output_directory']=str(self.root/'fit')
        result=run_project(p,base_directory=self.root,run_shear_strength_fit=True,save_plots=False)
        self.assertAlmostEqual(result['fit']['fitted_shear_strength_pa'],known)
        book=load_workbook(self.root/'fit/fit_results.xlsx',data_only=True)
        self.assertEqual(book.sheetnames,['拟合概况','各载荷独立诊断','各实验点','正向响应','所用实验'])
        self.assertAlmostEqual(book['拟合概况']['B2'].value,known)
        self.assertAlmostEqual(book['拟合概况']['B3'].value,known/1e6)
        self.assertEqual(book['各实验点'].max_row,4);book.close()


if __name__=='__main__':unittest.main()
