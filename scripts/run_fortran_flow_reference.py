#!/usr/bin/env python3
"""Run full flow-factor calculation from the first uploaded Fortran ZIP.

Unmodified source is the default. Optional ERF acceleration and diagnostic-only
print correction are explicit, with exact patches/hashes saved in manifest.json.
"""
import argparse
from datetime import datetime,timezone
import json
from pathlib import Path
import subprocess

import run_fortran_reference as ref

SOURCE_SHA='99d3b2e2c47a663ffe8b65bbdb7415db5eeebb434debd0ef404006600e861854'


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--archive',required=True);p.add_argument('--psd',required=True)
    p.add_argument('--output',required=True);p.add_argument('--compiler',default='gfortran')
    p.add_argument('--step',type=float,default=.01);p.add_argument('--nodes',type=int,default=12)
    p.add_argument('--maximum-pressure-pa',type=float,default=2e5)
    p.add_argument('--pdf-points',type=int,default=512,help='Original nurange input; internal PDF has four times this count')
    p.add_argument('--iterations',type=int,default=6000)
    p.add_argument('--real-erf',action='store_true');p.add_argument('--fix-diagnostic-print',action='store_true')
    p.add_argument('--anchor-pressure-pa',type=float,help='Single full-program iflag008=2 anchor, no flow-table extrapolation')
    p.add_argument('--timeout',type=float,default=1800)
    p.add_argument('--prepare-only',action='store_true',help='Extract and record inputs without compiling or running')
    args=p.parse_args()
    if not (4<=args.nodes<=300 and 32<=args.pdf_points<=10000 and .002<=args.step<=.1
            and args.iterations>0 and args.maximum_pressure_pa>0 and args.timeout>0):
        p.error('Invalid reference resolution or pressure')
    args.mode=3;args.delzeta=args.step
    manifest={'schema_version':1,'created_at':datetime.now(timezone.utc).isoformat()}
    out,build,source=ref.prepare(args,manifest)
    if manifest['source']['original_sha256']!=SOURCE_SHA:raise ValueError('Unexpected original source version')
    def change(name,updates):
        path=out/name;before=path.read_text();rows=before.splitlines()
        for index,value in updates.items():rows[index]=ref.replace_first_token(rows[index],value)
        path.write_text('\n'.join(rows)+'\n')
        manifest.setdefault('flow_control_changes',{})[name]={'before':before,'after':path.read_text()}
    anchor=args.anchor_pressure_pa is not None
    change('IN.flagparameters',{4:0,9:2 if anchor else 0,17:0 if anchor else 1,18:0})
    change('IN.FLOWFACTORS',{0:2,1:.01,2:.001,3:args.nodes,4:args.maximum_pressure_pa,5:0 if anchor else 6,6:1.})
    controls={13:args.pdf_points,18:.9,19:args.iterations}
    if anchor:
        if args.anchor_pressure_pa<=0:p.error('Anchor pressure must be positive')
        controls[2]=1
        change('IN.SEALING',{1:1})
        change('IN.sigma.gammabar.xnu.ElasticModul.yieldstress',{0:args.anchor_pressure_pa})
    change('IN.mathematical',controls)
    if args.fix_diagnostic_print:
        raw=source.read_bytes();old=b'      WRITE(555,8276) LowPressure';new=b'      WRITE(555,8276) DBLE(LowPressure)'
        if raw.count(old)!=1:raise ValueError('Unrecognized diagnostic WRITE')
        source.write_bytes(raw.replace(old,new))
        manifest['source'].update(modified=True,variant=('real-erf+' if args.real_erf else '')+'diagnostic-real-cast',
            printing_fix={'before':old.decode(),'after':new.decode(),'scope':'Only final diagnostic formatting; no physical calculation changed'})
    manifest['source']['compiled_sha256']=ref.sha256(source)
    manifest['purpose']='Complete original-program PSD/contact/PDF/flow-factor chain; not an extracted kernel'
    manifest['used_input_hashes']={f.name:ref.sha256(f) for f in out.glob('IN*') if f.is_file()}
    manifest['status']='prepared';ref.save_manifest(out,manifest)
    if args.prepare_only:
        print(json.dumps({'status':'prepared','output':str(out)}))
        return 0
    manifest['compiler']=subprocess.check_output([args.compiler,'--version'],text=True).splitlines()[0]
    manifest['compile']=ref.run_command([args.compiler,*ref.COMPILE_FLAGS,str(source),'-o',str(build/'original')],
        build,out/'compile.log',180)
    ref.save_manifest(out,manifest)
    if manifest['compile']['returncode']!=0:raise RuntimeError('Fortran compilation failed; see compile.log')
    manifest['binary_sha256']=ref.sha256(build/'original');manifest['status']='running';ref.save_manifest(out,manifest)
    manifest['run']=ref.run_command([str(build/'original')],out,out/'run.log',args.timeout)
    flow=list(out.glob('out.FLUID.FLOWFACTORS.*'))
    expected=1 if anchor else args.nodes+6
    rows=ref.numeric_rows(flow[0],8) if len(flow)==1 else []
    complete=len(rows)==expected
    if complete and not anchor:
        tests=list(out.glob('TESTING.1u.2p.*'))
        complete=len(tests)==1 and tests[0].stat().st_size>0
    manifest['validation']={'table_complete':complete,'expected_rows':expected,'observed_rows':len(rows),
        'passed':complete and manifest['run']['returncode']==0 and not manifest['run']['timed_out']}
    manifest['status']='validated' if manifest['validation']['passed'] else 'failed'
    manifest['output_sha256']={f.name:ref.sha256(f) for f in flow}
    ref.save_manifest(out,manifest)
    print(json.dumps({'status':manifest['status'],'run':manifest['run'],'validation':manifest['validation']}))
    return 0 if manifest['validation']['passed'] else 2


if __name__=='__main__':raise SystemExit(main())
