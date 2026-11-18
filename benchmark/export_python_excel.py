"""Export forward components to Excel; use the frozen model's existing CSV outputs."""
import argparse
import csv
from pathlib import Path
from openpyxl import Workbook
from openpyxl.chart import ScatterChart, Series, Reference

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('directory', type=Path)
    args = parser.parse_args()
    with (args.directory / 'mixed_cases.csv').open(newline='', encoding='utf-8-sig') as f:
        raw = list(csv.DictReader(f))
    if not raw:
        parser.error('mixed_cases.csv has no results')
    wb = Workbook()
    ws = wb.active
    ws.title = '正向摩擦分量'
    ws.append(['speed_m_s', 'mu_solid', 'mu_fluid_sliding', 'mu_fluid_pressure',
               'mu_total', 'real_contact_area_per_width_m', 'fluid_load_fraction',
               'minimum_mean_gap_m', 'load_relative_error', 'mass_relative_residual'])
    for row in raw:
        v = float(row['lower_wall_velocity_m_s'])
        load = float(row['contact_line_load_n_m']) + float(row['fluid_line_load_n_m'])
        # Signed lower-body resistance; do not take the absolute value of each fluid term.
        sign = -1.0 if v > 0 else 1.0
        ws.append([abs(v), *(sign * float(row[key]) / load for key in (
            'lower_solid_force_n_m', 'lower_fluid_sliding_force_n_m',
            'lower_fluid_pressure_flow_force_n_m')),
            float(row['lower_resisting_friction_coefficient']),
            *(float(row[key]) for key in ('real_contact_area_per_width_m', 'fluid_load_fraction',
                                         'minimum_mean_gap_m', 'load_relative_error',
                                         'mass_relative_residual'))])
    ws.freeze_panes = 'A2'
    ws.auto_filter.ref = ws.dimensions
    native = wb.create_sheet('原始正向输出')
    native.append(list(raw[0]))
    for row in raw:
        native.append([_number(value) for value in row.values()])
    note = wb.create_sheet('说明')
    note.append(['本表直接读取正向解；tau为输入的已知值。'])
    note.append(['拟合任务的forward目录通常使用tau_ref=1MPa，不是拟合后的tau。'])
    note.append(['面积列是Ar/B，乘轴向宽度B才得到m²。'])
    chart = ScatterChart()
    chart.title = 'Lower-plane friction components'
    chart.x_axis.title = 'Speed [m/s]'
    chart.y_axis.title = 'Resisting friction coefficient'
    chart.x_axis.scaling.logBase = 10
    for column in range(2, 6):
        chart.series.append(Series(Reference(ws, min_col=column, min_row=1, max_row=ws.max_row),
                                   Reference(ws, min_col=1, min_row=2, max_row=ws.max_row),
                                   title_from_data=True))
    ws.add_chart(chart, 'L2')
    target = args.directory / 'forward_components.xlsx'
    wb.save(target)
    print(target)

def _number(value):
    try:
        return float(value)
    except ValueError:
        return value

if __name__ == '__main__':
    main()
