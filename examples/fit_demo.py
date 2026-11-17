#!/usr/bin/env python3
"""Fast local installation/fit demo; creates SYNTHETIC data at run time.

No mixed physics is simulated here and no experiment is supplied. A prescribed
area and fluid coefficient provide an exact 0.3 MPa inverse-problem benchmark.
"""
import argparse
import csv
import json
from pathlib import Path
import sys
import tempfile

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'src'))

from persson_contact.workflow import run_project


def write_demo_inputs(directory):
    directory=Path(directory);forward=directory/'synthetic_forward';forward.mkdir()
    cases=[];experiment=[]
    for i,(speed,area,fluid,fluid_load) in enumerate(((1e-4,2e-4,.01,10.),(1e-3,1.2e-4,.04,50.),(1e-2,5e-5,.12,90.))):
        cases.append({'lower_wall_velocity_m_s':speed,'real_contact_area_per_width_m':area,
            'lower_fluid_force_n_m':-fluid*100.,'contact_line_load_n_m':100.-fluid_load,
            'fluid_line_load_n_m':fluid_load,'load_relative_error':0.,'geometry_relative_residual':0.,'mass_relative_residual':0.})
        experiment.append({'speed_m_s':speed,'friction_coefficient':fluid+3e5*area/100.,'use_for_fit':i<2})
    for path,rows in ((forward/'mixed_cases.csv',cases),(directory/'experiment.csv',experiment)):
        with path.open('w',encoding='utf-8',newline='') as stream:
            writer=csv.DictWriter(stream,fieldnames=list(rows[0]),lineterminator='\n');writer.writeheader();writer.writerows(rows)
    diag={'workflow_status':'complete','synthetic_demo':True,'target_line_load_n_m':100.,'cases':cases,
        'parameters_snapshot':{'description':'Synthetic mathematical fit benchmark; NOT a mixed-lubrication prediction.'},
        'limitations':['Synthetic area/fluid data. No material, PSD, cavitation or fluid-solid calculation was performed.']}
    (forward/'mixed_diagnostics.json').write_text(json.dumps(diag,indent=2)+'\n')
    return directory/'experiment.csv',forward


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,default=ROOT/'results/fit_demo')
    parser.add_argument('--no-plots',action='store_true')
    args=parser.parse_args()
    with tempfile.TemporaryDirectory(prefix='ecl_fit_demo_') as directory:
        experiment,forward=write_demo_inputs(directory)
        p={'macro_contact':{'line_load_n_m':100.},'contact_friction':{'transverse_width_m':.01},
            'shear_strength_fit':{'experiment_file':str(experiment),'source':'saved_forward',
                'source_parameters':{'saved_forward':{'directories':[str(forward)]}},
                'reference_strength_pa':1e6},'output_directory':str(args.output.resolve())}
        result=run_project(p,base_directory=ROOT,run_contact=False,run_shear_strength_fit=True,save_plots=not args.no_plots)
    actual=result['fit']['fitted_shear_strength_pa']
    if abs(actual/3e5-1)>1e-12:raise AssertionError('Synthetic known-strength recovery failed')
    print('合成数据演示完成：已知0.3MPa，拟合0.3MPa。此演示不代表真实实验或混合润滑物理验证。')


if __name__=='__main__':main()
