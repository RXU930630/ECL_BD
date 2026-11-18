#!/usr/bin/env python3
"""Create a blank Excel template, or convert an existing comma/semicolon CSV."""
import argparse
import csv
from pathlib import Path
import sys

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'src'))


def main():
    from persson_contact.experiment_excel import HEADERS,write_experiment_workbook
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--csv',type=Path,default=None)
    args=parser.parse_args();rows=[]
    if args.csv:
        with args.csv.open(encoding='utf-8-sig',newline='') as stream:
            first=stream.readline();stream.seek(0)
            reader=csv.DictReader(stream,delimiter=';' if ';' in first and ',' not in first else ',')
            header=reader.fieldnames or []
            if (len(set(header))!=len(header) or not {'speed_m_s','friction_coefficient'}<=set(header)
                or set(header)-set(HEADERS)):raise ValueError('请使用实验模板的列名')
            for row in reader:
                if None in row or any(v is None for v in row.values()):raise ValueError('CSV列数不匹配')
                converted={}
                for key,value in row.items():
                    value=value.strip()
                    converted[key]=(None if not value else value.lower() in ('true','1')) if key=='use_for_fit' else (None if not value else float(value))
                    if key=='use_for_fit' and value.lower() not in ('','true','false','1','0'):raise ValueError('Invalid use_for_fit')
                rows.append(converted)
    write_experiment_workbook(args.output,rows)
    print(f'Excel文件已保存：{args.output.resolve()}')


if __name__=='__main__':main()
