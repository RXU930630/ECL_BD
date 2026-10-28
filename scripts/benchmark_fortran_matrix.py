#!/usr/bin/env python3
"""Run controlled elastic Fortran cases; keep source variants/provenance explicit.

The batch variant changes only real DERF1 to intrinsic ERF and writes every PDF
point instead of every 100th point. No contact, separation or PDF formula changes.
Use run_fortran_reference.py without --real-erf for the independent pristine anchor.
"""
import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
import json
from pathlib import Path
import shutil
from types import SimpleNamespace

import run_fortran_reference as reference


def prepare_case(args, name, young, pressure, mode, step):
    manifest = {'schema_version': 1, 'status': 'preparing'}
    out, build, source = reference.prepare(SimpleNamespace(
        archive=args.archive, psd=args.psd, output=str(Path(args.output)/name),
        mode=mode, delzeta=step, real_erf=True), manifest)
    material = out/reference.MATERIAL_NAME
    lines = material.read_text().splitlines()
    lines[0], lines[3] = f'{pressure:.17g}', f'{young:.17g}'
    material.write_text('\n'.join(lines)+'\n')
    manifest['controls'][reference.MATERIAL_NAME]['used'] = material.read_text()
    data = source.read_bytes()
    before = data
    first = b'      if (j.eq.1) write(7037,*) uu, Prob'
    thinning = b'      if (j.eq.100*(j/100)) write(7037,*) uu, Prob'
    # This source version has the same output sites in the two separation
    # methods; instrument both, without changing either method's arithmetic.
    if data.count(first) != 2 or data.count(thinning) != 2:
        raise ValueError('Unrecognized PDF write sites; refusing output instrumentation')
    data = data.replace(first, b'!     First point is included in the dense diagnostic output below.')
    data = data.replace(thinning, b'      write(7037,*) uu, Prob')
    source.write_bytes(data)
    manifest['source']['variant'] = 'real-erf+dense-pdf-output'
    manifest['source']['output_instrumentation'] = 'Write every internal PDF point; remove duplicate first-point write; no arithmetic changes.'
    manifest['source']['output_patch_counts'] = {'first_point_site': before.count(first), 'thinning_site': before.count(thinning)}
    manifest['source']['compiled_sha256'] = reference.sha256(source)
    manifest['unchanged_model_parameters'] = 'Poisson ratio and all mathematical controls except requested delzeta and disabled polynomial smoothing; E and nominal pressure explicitly overridden.'
    manifest['case'] = {'young_modulus_pa': young, 'poisson_ratio': .5,
                        'effective_modulus_pa': young/.75, 'pressure_pa': pressure,
                        'mode': mode, 'delzeta': step}
    manifest['used_input_hashes'] = {p.name: reference.sha256(p) for p in out.glob('IN*') if p.is_file()}
    reference.save_manifest(out, manifest)
    return out, build, source, manifest


def run_case(case, binary, source_hash, args):
    out, build, source, manifest = case
    if reference.sha256(source) != source_hash:
        raise ValueError('Reused binary does not match prepared source hash')
    target = build/'contact_mechanics'
    shutil.copy2(binary, target)
    manifest['compiler'] = {'compiled_here': False, 'binary_origin': str(binary),
                            'corresponding_compiled_source_sha256': source_hash}
    manifest['binary'] = {'path': str(target.relative_to(out)), 'sha256': reference.sha256(target)}
    manifest['status'] = 'running'
    reference.save_manifest(out, manifest)
    manifest['run'] = reference.run_command([str(target)], out, out/'run.log', args.timeout)
    manifest['validation'] = reference.validate_outputs(out, manifest['mode'])
    manifest['source']['sha256_after_run'] = reference.sha256(source)
    good = (manifest['run']['returncode'] == 0 and not manifest['run']['timed_out']
            and manifest['validation']['passed'] and reference.sha256(source) == source_hash)
    manifest['status'] = 'validated' if good else 'failed'
    reference.save_manifest(out, manifest)
    return {'case': out.name, 'status': manifest['status'], 'seconds': manifest['run']['elapsed_seconds']}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--archive', required=True)
    parser.add_argument('--psd', required=True)
    parser.add_argument('--output', required=True)
    parser.add_argument('--compiler', default='gfortran')
    parser.add_argument('--jobs', type=int, default=2)
    parser.add_argument('--timeout', type=float, default=1200)
    args = parser.parse_args()
    root = Path(args.output).resolve()
    if root.exists():
        parser.error('Output directory must not already exist')
    root.mkdir(parents=True)
    args.output = str(root)
    jobs = [('E3_p100k_mode2_d01', 3e6, 1e5, 2, .01),
            ('E3_p100k_mode3_d01', 3e6, 1e5, 3, .01)]
    for e, label in ((1.5e6, 'E1p5'), (3e6, 'E3'), (9e6, 'E9')):
        if e != 3e6:
            jobs.append((f'{label}_p100k_mode2_d01', e, 1e5, 2, .01))
        for p, plabel in ((1e4,'p10k'), (5e4,'p50k'), (1e5,'p100k')):
            if e == 3e6 and p == 1e5:
                continue
            jobs.append((f'{label}_{plabel}_mode3_d01', e, p, 3, .01))
    jobs.extend([('E3_p100k_mode2_d005', 3e6, 1e5, 2, .005),
                 ('E3_p100k_mode3_d005', 3e6, 1e5, 3, .005)])
    cases = [prepare_case(args, *job) for job in jobs]
    first = cases[0]
    binary = root/'contact_mechanics_batch'
    compiled = reference.run_command([args.compiler, *reference.COMPILE_FLAGS,
                                       first[2].name, '-o', str(binary)],
                                      first[1], root/'compile.log', 300)
    (root/'compile_manifest.json').write_text(json.dumps({
        'compile': compiled, 'source': first[3]['source'],
        'compiler_flags': reference.COMPILE_FLAGS}, indent=2)+'\n')
    if compiled['returncode'] != 0 or compiled['timed_out']:
        raise RuntimeError('Batch source compilation failed')
    source_hash = first[3]['source']['compiled_sha256']
    results = []
    with ThreadPoolExecutor(max_workers=args.jobs) as pool:
        futures = [pool.submit(run_case, case, binary, source_hash, args) for case in cases]
        for future in as_completed(futures):
            result = future.result()
            results.append(result)
            print(json.dumps(result), flush=True)
            (root/'status.json').write_text(json.dumps(results, indent=2)+'\n')
    return 0 if all(result['status'] == 'validated' for result in results) else 1


if __name__ == '__main__':
    raise SystemExit(main())
