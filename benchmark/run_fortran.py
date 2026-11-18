"""Compile the user's pristine original source and retain every native output."""
import argparse
import hashlib
import json
import math
from pathlib import Path
import shutil
import subprocess
import sys

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source', type=Path, required=True,
                        help='Path to AI.ContactMechanics.14.09.2026.f90 from the first ZIP')
    parser.add_argument('--compiler', default='gfortran')
    parser.add_argument('--output', type=Path, default=ROOT / 'benchmark_runs' / 'fortran')
    args = parser.parse_args()
    source = args.source.resolve()
    manifest = json.loads((HERE / 'manifest.json').read_text(encoding='utf-8'))
    digest = hashlib.sha256(source.read_bytes()).hexdigest()
    if digest != manifest['fortran_source_sha256']:
        parser.error('Source SHA256 differs from the original first uploaded archive.')
    target = args.output.resolve()
    target.mkdir(parents=True, exist_ok=False)
    for item in (HERE / 'fortran_inputs').iterdir():
        shutil.copyfile(item, target / item.name)
    # Keep the include next to the copied, byte-identical source.
    copied_source = target / source.name
    shutil.copyfile(source, copied_source)
    executable = target / ('original.exe' if sys.platform == 'win32' else 'original')
    command = [args.compiler, '-O3', '-mcmodel=large', '-ffree-line-length-none',
               '-std=legacy', '-fallow-argument-mismatch', str(copied_source), '-o', str(executable)]
    with (target / 'compile.log').open('w', encoding='utf-8') as log:
        subprocess.run(command, cwd=target, stdout=log, stderr=subprocess.STDOUT, check=True)
    print('Running pristine Fortran. BoF precomputation can take tens of minutes.', flush=True)
    with (target / 'run.log').open('w', encoding='utf-8') as log:
        run = subprocess.run([str(executable)], cwd=target, stdout=log, stderr=subprocess.STDOUT)
    files = list(target.glob('out.MSLUBE.LOW.WITH.FrictionFlowfactor.*'))
    rows = []
    if len(files) == 1:
        for line in files[0].read_text().splitlines():
            values = [float(x.replace('D', 'E')) for x in line.split()]
            if values:
                rows.append(values)
    # The original opens several ERROR files unconditionally; empty files are normal.
    errors = [p.name for p in target.glob('ERROR*') if p.is_file() and p.stat().st_size]
    expected = manifest['physical']['speeds_m_s']
    rows_ok = len(rows) == len(expected) and all(
        len(row) == 16 and all(math.isfinite(x) for x in row[:8])
        and math.isclose(row[0], speed, rel_tol=1e-3)
        for row, speed in zip(rows, expected))
    auxiliary_nonfinite = sorted({i + 1 for row in rows for i, x in enumerate(row)
                                  if i >= 8 and not math.isfinite(x)})
    record = {'source_sha256': digest, 'command': command,
              'returncode': run.returncode, 'native_rows': len(rows),
              'error_files': errors, 'expected_primary_rows_present_and_finite': rows_ok,
              'nonfinite_auxiliary_columns': auxiliary_nonfinite,
              'convergence_certified': False,
              'note': 'Completion and exit status alone do not establish solver or mesh convergence.'}
    (target / 'execution_record.json').write_text(json.dumps(record, indent=2) + '\n')
    if run.returncode or errors or not rows_ok:
        print('Incomplete or invalid native outputs; inspect execution_record.json and native logs.')
        return 2
    if auxiliary_nonfinite:
        print('WARNING: nonfinite auxiliary body/slope-force columns:', auxiliary_nonfinite)
        print('Only columns 1-8 are usable for the selected lower-plane comparison.')
    print('Native rows present. Check residuals and refinement before treating this as validated.')
    return 0

if __name__ == '__main__':
    raise SystemExit(main())
