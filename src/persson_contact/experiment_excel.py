"""Excel experiment input and human-readable fit output; no equation changes."""
from contextlib import closing
from datetime import datetime
import json
from pathlib import Path
from zipfile import BadZipFile


HEADERS = ('speed_m_s','friction_coefficient','line_load_n_m','normal_force_n',
    'background_friction_coefficient','standard_uncertainty_mu','use_for_fit')


def read_excel_table(path):
    """Read experiment sheet (otherwise first sheet), using saved formula values."""
    from openpyxl import load_workbook
    from openpyxl.utils.exceptions import InvalidFileException
    try:
        with closing(load_workbook(path,read_only=True,data_only=False,keep_links=False)) as formulas, \
             closing(load_workbook(path,read_only=True,data_only=True,keep_links=False)) as values:
            name='experiment' if 'experiment' in formulas.sheetnames else formulas.sheetnames[0]
            source=formulas[name];saved=values[name]
            header=None;rows=[]
            for number,(cells,cached) in enumerate(zip(source.iter_rows(),saved.iter_rows()),1):
                row=[]
                for cell,result in zip(cells,cached):
                    if cell.data_type=='f' and result.value is None:
                        raise ValueError(f'实验Excel {name}!{cell.coordinate}的公式没有保存的计算值；请用Excel重新计算并保存，或复制后粘贴为数值。')
                    if result.data_type=='e':
                        raise ValueError(f'实验Excel {name}!{cell.coordinate}包含错误值{result.value}；请修正该单元格。')
                    row.append(result.value)
                if number==1:
                    while row and row[-1] is None:row.pop()
                    header=row
                else:
                    if not any(value is not None and str(value).strip() for value in row):continue
                    # Excel formatting often creates empty trailing columns.
                    while len(row)>len(header) and row[-1] is None:row.pop()
                    row.extend([None]*max(0,len(header)-len(row)))
                    rows.append((number,row))
            return header or [],rows
    except (BadZipFile,InvalidFileException,KeyError) as exc:
        raise ValueError('无法读取Excel工作簿；请使用Excel另存为真正的.xlsx文件，不能只改扩展名。') from exc


def _workbook():
    from openpyxl import Workbook
    book=Workbook();book.remove(book.active)
    book.properties.creator='ECL_BD'
    book.properties.created=book.properties.modified=datetime(2026,1,1)
    return book


def _sheet(book,name,headers,rows):
    from openpyxl.styles import Alignment,Font,PatternFill
    sheet=book.create_sheet(name);sheet.append(list(headers))
    for row in rows:sheet.append(list(row))
    sheet.freeze_panes='A2';sheet.auto_filter.ref=sheet.dimensions
    for cell in sheet[1]:
        cell.font=Font(bold=True,color='FFFFFF');cell.fill=PatternFill('solid',fgColor='24476B')
        cell.alignment=Alignment(wrap_text=True,vertical='center')
        sheet.column_dimensions[cell.column_letter].width=min(46,max(19,len(str(cell.value))+2))
    sheet.row_dimensions[1].height=34
    return sheet


def write_experiment_workbook(path,rows=()):
    """Blank template or user-supplied values; input columns always use SI units."""
    from openpyxl.comments import Comment
    from openpyxl.styles import Alignment
    from openpyxl.worksheet.datavalidation import DataValidation
    path=Path(path)
    if path.exists():raise FileExistsError(f'不覆盖已有Excel文件：{path}')
    book=_workbook();rows=list(rows)
    sheet=_sheet(book,'experiment',HEADERS,[[row.get(key) for key in HEADERS] for row in rows])
    notes=[
        ('speed_m_s','必填。速度大小，m/s；1微米/秒填0.000001。不要填log速度。'),
        ('friction_coefficient','必填。实验摩擦系数mu，无单位。'),
        ('line_load_n_m','可空。线载荷N/m；填写总力时建议留空此列。'),
        ('normal_force_n','总法向力N；逐行填写40、100、160等。轴向长度由parameters读取。'),
        ('background_friction_coefficient','可空。已知背景摩擦系数；空白表示0。'),
        ('standard_uncertainty_mu','等权拟合可空；按测量不确定度加权时填写正值。'),
        ('use_for_fit','TRUE参与拟合；FALSE留作验证；空白默认TRUE。')]
    for cell,(_,note) in zip(sheet[1],notes):cell.comment=Comment(note,'ECL_BD')
    for row in sheet.iter_rows(min_row=2,max_row=max(28,len(rows)+1),max_col=7):
        row[0].number_format='0.00E+00';row[1].number_format='0.000000'
        for cell in row[2:6]:cell.number_format='0.############'
    choice=DataValidation(type='list',formula1='"TRUE,FALSE"',allow_blank=True)
    sheet.add_data_validation(choice);choice.add('G2:G10001')
    instructions=[('使用方法','在experiment工作表中从第2行填写。保存为.xlsx，放到inputs/experiment.xlsx。'),
        ('多载荷','一个工作表可包含多档载荷；每个速度/载荷组合一行，不分成多个载荷工作表。'),
        ('表头','第1行英文表头保持不变；中文解释见本页或表头批注。'),
        ('公式','读取Excel保存的计算值，不重算公式；公式报错或没有保存计算值时会提示单元格。'),
        ('拟合结果','计算结束打开输出目录中的fit_results.xlsx。剪切强度同时显示Pa和MPa。'),*notes]
    help_sheet=_sheet(book,'填写说明',('项目','说明'),instructions)
    help_sheet.column_dimensions['A'].width=36;help_sheet.column_dimensions['B'].width=100
    for row in help_sheet.iter_rows(min_row=2):row[1].alignment=Alignment(wrap_text=True)
    path.parent.mkdir(parents=True,exist_ok=True);book.save(path);book.close()


def write_fit_workbook(path,summary,records,responses,experiments):
    book=_workbook()
    def value(v):
        if isinstance(v,(dict,list,tuple)):return json.dumps(v,ensure_ascii=False)
        # Preserve external text as literal cells, never as Excel formulas.
        return v
    top=[('拟合剪切强度 tau [Pa]',summary['fitted_shear_strength_pa']),
        ('拟合剪切强度 tau [MPa]',summary['fitted_shear_strength_pa']/1e6)]
    top.extend((key,value(v)) for key,v in summary.items())
    _sheet(book,'拟合概况',('项目','数值'),top)
    groups=summary.get('per_load_diagnostics',[])
    if groups:_sheet(book,'各载荷独立诊断',list(groups[0]),[[row[key] for key in groups[0]] for row in groups])
    for name,rows in [('各实验点',records),('正向响应',responses),('所用实验',experiments)]:
        if rows:_sheet(book,name,tuple(rows[0]),[[value(row[key]) for key in rows[0]] for row in rows])
    for sheet in book:
        for row in sheet:
            for cell in row:
                if isinstance(cell.value,str):cell.data_type='s'
    book.save(path);book.close()
