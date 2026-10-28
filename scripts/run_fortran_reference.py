#!/usr/bin/env python3
"""Reproduce a controlled Fortran dry-contact reference run.

Only an isolated copy of the IN* controls is changed.  The original archive,
PSD and archived Fortran source remain untouched.  By default the compiled source
is byte-identical.  --real-erf explicitly changes only DERF1 in the build copy and
records the exact patch and both hashes.  This script requires Python's standard
library and a working GNU Fortran compiler; it does not install a compiler.
"""

import argparse
from datetime import datetime, timezone
import difflib
import hashlib
import json
import math
from pathlib import Path, PurePosixPath
import shutil
import stat
import subprocess
import sys
import time
import zipfile


PSD_NAME = "IN_ROADTYPE.q0.qr.q1.h0.H.1logWavevector.2logPowerSpectrum"
FLAG_NAME = "IN.flagparameters"
MATH_NAME = "IN.mathematical"
POLY_NAME = "IN_npolyCq.npolyab.npolyE"
MATERIAL_NAME = "IN.sigma.gammabar.xnu.ElasticModul.yieldstress"
SOURCE_PREFIX = "AI.ContactMechanics."
COMPILE_FLAGS = ["-mcmodel=large", "-ffree-line-length-none", "-std=legacy",
                 "-fallow-argument-mismatch", "-O3"]
MAX_SELECTED_BYTES = 64 * 1024 * 1024

DERF1_ORIGINAL = (b"     real*8 function DERF1(x)\n"
                  b"     real*8 x\n"
                  b"     complex*16 zg,zgg \n"
                  b"!      y=DERF(x)\n"
                  b"      zg=x\n"
                  b"      call cerror(zg,zgg)\n"
                  b"      DERF1=zgg\n"
                  b"     return\n"
                  b"     end\n")
DERF1_REAL_ERF = (b"     real*8 function DERF1(x)\n"
                  b"     real*8 x\n"
                  b"      DERF1=ERF(x)\n"
                  b"     return\n"
                  b"     end\n")


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _zip_path(info):
    path = PurePosixPath(info.filename)
    if (path.is_absolute() or ".." in path.parts or "\\" in info.filename
            or not path.name or ":" in path.parts[0]):
        raise ValueError(f"Unsafe archive path: {info.filename!r}")
    mode = info.external_attr >> 16
    if stat.S_ISLNK(mode):
        raise ValueError(f"Symlink entries are not accepted: {info.filename!r}")
    return path


def archive_members(archive):
    """Select one source and its sibling IN* files, rejecting ambiguous inputs."""
    files = [info for info in archive.infolist() if not info.is_dir()]
    sources = [info for info in files
               if PurePosixPath(info.filename).name.startswith(SOURCE_PREFIX)
               and PurePosixPath(info.filename).suffix.lower() == ".f90"]
    if len(sources) != 1:
        raise ValueError(f"Expected exactly one {SOURCE_PREFIX}*.f90 source; found {len(sources)}")
    source_path = _zip_path(sources[0])
    selected = [sources[0]]
    selected.extend(info for info in files
                    if PurePosixPath(info.filename).parent == source_path.parent
                    and PurePosixPath(info.filename).name.startswith("IN"))
    names = []
    total = 0
    for info in selected:
        path = _zip_path(info)
        names.append(path.name)
        total += info.file_size
    if len(names) != len(set(names)):
        raise ValueError("Archive contains duplicate selected filenames")
    required = {"IN.PARAMETER", FLAG_NAME, MATH_NAME, POLY_NAME, MATERIAL_NAME, PSD_NAME,
                "IN.HERTZ", "IN.MIXEDLUBRICATION", "IN.MIXEDLUBRICATION.mathematical",
                "IN.FLUID.viscosity.and.TrippNumber", "IN.FLOWFACTORS", "IN.SEALING",
                "IN.ContactRadiusAndWear", "IN.ADHESION.mathematical", "IN.CAPILLARY",
                "IN.addCq", "IN.BREAKLOOSE"}
    missing = required - set(names)
    if missing:
        raise ValueError(f"Archive lacks required original inputs: {', '.join(sorted(missing))}")
    if total > MAX_SELECTED_BYTES:
        raise ValueError("Selected source and IN* files exceed the 64 MiB safety limit")
    return sources[0], selected


def replace_first_token(line, value):
    parts = line.split(maxsplit=1)
    if not parts:
        raise ValueError("Expected a numeric control value, found an empty row")
    return str(value) + ("  " + parts[1] if len(parts) == 2 else "")


def _control_text(path):
    return path.read_text(encoding="utf-8-sig")


def modify_controls(output, mode, delzeta):
    names = (FLAG_NAME, MATH_NAME, POLY_NAME, MATERIAL_NAME)
    before = {name: _control_text(output / name) for name in names}
    rows = before[FLAG_NAME].splitlines()
    if len(rows) != 20:
        raise ValueError("IN.flagparameters must have the 20 rows of this reference version")
    # Row 1 is temperature.  All 19 physical flags are zero except PSD input
    # and the requested dry separation mode.  iNonlinear is included here.
    for index in range(1, 20):
        rows[index] = replace_first_token(rows[index], 1 if index == 3 else mode if index == 4 else 0)
    (output / FLAG_NAME).write_text("\n".join(rows) + "\n", encoding="utf-8")
    rows = before[POLY_NAME].splitlines()
    if len(rows) < 3:
        raise ValueError("Polynomial input must contain three control rows")
    rows[0] = replace_first_token(rows[0], 0)
    (output / POLY_NAME).write_text("\n".join(rows) + "\n", encoding="utf-8")
    rows = before[MATH_NAME].splitlines()
    if len(rows) < 22:
        raise ValueError("Mathematical input must contain the complete original control set")
    rows[0] = replace_first_token(rows[0], format(delzeta, ".17g"))
    (output / MATH_NAME).write_text("\n".join(rows) + "\n", encoding="utf-8")
    after = {name: _control_text(output / name) for name in names}
    return {name: {"original": before[name], "used": after[name]} for name in names}


def save_manifest(output, manifest):
    (output / "manifest.json").write_text(
        json.dumps(manifest, ensure_ascii=False, indent=2, allow_nan=False) + "\n", encoding="utf-8")


def apply_real_erf(source):
    """Apply one narrowly recognized real-only DERF1 replacement, or fail closed."""
    data = source.read_bytes()
    # The full exact function must be unique.  Unknown versions are refused;
    # there is no regex/global substitution that could touch CERROR or callers.
    newline = b"\r\n" if b"\r\n" in data else b"\n"
    old = DERF1_ORIGINAL.replace(b"\n", newline)
    new = DERF1_REAL_ERF.replace(b"\n", newline)
    if data.count(old) != 1:
        raise ValueError("--real-erf requires exactly one recognized original DERF1 function; source version differs")
    replacement = data.replace(old, new, 1)
    source.write_bytes(replacement)
    before, after = old.decode("ascii"), new.decode("ascii")
    return {"name": "real-erf", "description": "Replace the real-valued DERF1 wrapper around complex CERROR with standard intrinsic ERF(x); CERROR and all other source bytes remain unchanged",
            "before": before, "after": after,
            "diff": "".join(difflib.unified_diff(before.splitlines(keepends=True),
                                                  after.splitlines(keepends=True),
                                                  fromfile="DERF1.original", tofile="DERF1.real-erf"))}


def prepare(args, manifest):
    archive_path, psd_path = Path(args.archive).resolve(), Path(args.psd).resolve()
    if not archive_path.is_file() or not psd_path.is_file():
        raise ValueError("--archive and --psd must be existing files")
    output = Path(args.output).absolute()
    # Even an empty pre-existing directory is rejected, so outputs cannot be
    # confused with an old run and no user files can be overwritten.
    if output.exists() or output.is_symlink():
        raise FileExistsError(f"Refusing to overwrite existing output: {output}")
    with zipfile.ZipFile(archive_path) as archive:
        source_info, members = archive_members(archive)
        output.mkdir(parents=True, exist_ok=False)
        manifest["working_directory"] = str(output)
        build = output / "build"
        build.mkdir()
        selected = {}
        for info in members:
            name = PurePosixPath(info.filename).name
            destination = build / name if info == source_info else output / name
            destination.write_bytes(archive.read(info))
            selected[name] = {"archive_member": info.filename, "sha256_original": sha256(destination)}
        source = build / PurePosixPath(source_info.filename).name
    manifest.update({
        "archive": {"path": str(archive_path), "sha256": sha256(archive_path)},
        "psd": {"path": str(psd_path), "sha256": sha256(psd_path), "input_name": PSD_NAME},
        "source": {"file": str(source.relative_to(output)), "sha256": sha256(source),
                   "original_sha256": sha256(source), "modified": False, "variant": "original"},
        "selected_original_files": selected,
        "controls": modify_controls(output, args.mode, args.delzeta),
        "mode": args.mode, "delzeta": args.delzeta,
        "compile_flags": COMPILE_FLAGS,
        "working_directory": str(output),
        "changes": ["Use supplied PSD in the fixed input filename", "Set numeric PSD input and dry separation mode only",
                    "Disable polynomial PSD smoothing", "Set requested logarithmic magnification step"],
        "unchanged_model_parameters": "All remaining IN.mathematical values, including gamma, afact and u1min; original material inputs",
    })
    shutil.copyfile(psd_path, output / PSD_NAME)
    shutil.copyfile(output / "IN.PARAMETER", build / "IN.PARAMETER")
    if args.real_erf:
        manifest["source"]["optimization"] = apply_real_erf(source)
        manifest["source"]["modified"] = True
        manifest["source"]["variant"] = "real-erf"
        manifest["changes"].append("Explicit real-erf optimization applied only to DERF1 in build source")
    manifest["source"]["compiled_sha256"] = sha256(source)
    manifest["used_input_hashes"] = {
        path.name: sha256(path) for path in sorted(output.glob("IN*")) if path.is_file()}
    manifest["status"] = "prepared"
    save_manifest(output, manifest)
    return output, build, source


def run_command(command, cwd, logfile, timeout):
    started = time.monotonic()
    timed_out = False
    with logfile.open("wb") as log:
        try:
            completed = subprocess.run(command, cwd=cwd, stdout=log, stderr=subprocess.STDOUT,
                                       timeout=timeout, check=False)
            returncode = completed.returncode
        except subprocess.TimeoutExpired:
            # subprocess.run kills and waits for its child before raising.
            timed_out = True
            returncode = None
    return {"command": [str(part) for part in command], "cwd": str(cwd),
            "log": str(logfile), "returncode": returncode, "timed_out": timed_out,
            "elapsed_seconds": time.monotonic() - started}


def numeric_rows(path, columns):
    rows = []
    if not path.is_file():
        return rows
    for line in path.read_text(errors="replace").splitlines():
        if line.lstrip().startswith(("#", "!")):
            continue
        tokens = line.split()
        if len(tokens) != columns:
            continue
        try:
            values = [float(token.replace("D", "E").replace("d", "e")) for token in tokens]
        except ValueError:
            continue
        if not all(math.isfinite(value) for value in values):
            raise ValueError(f"Non-finite values in {path.name}")
        rows.append(values)
    return rows


def validate_outputs(output, mode):
    specifications = [("out_No.Adhesion.1pressure.2A.3Aexact", 3, 1000)]
    if mode == 2:
        specifications.append(("RealUnits.1pressure.Pa.2average.interfacial.separation.meter.3A.4Aexact", 4, 10))
    else:
        specifications.extend([
            ("out.1u.2ProbabilityInterfacialSeparation", 2, 10),
            ("out.1pressure.2baruFromP.3baruFromEnergy.4AreaFromP.5AreaFromEnergy", 5, 1),
            ("out.1pressure.2zeta.3logzeta.4relArea.5separation.6hrms.zeta", 6, 10)])
    checks, failures = {}, []
    for name, columns, minimum in specifications:
        try:
            rows = numeric_rows(output / name, columns)
            failures_for_file = []
            if len(rows) < minimum:
                failures_for_file.append(f"requires >= {minimum} finite {columns}-column rows; got {len(rows)}")
            if rows:
                if name.startswith("out_No.Adhesion"):
                    if any(row[0] <= 0 or not 0 <= row[1] <= 1 for row in rows):
                        failures_for_file.append("pressure/standard contact area outside physical range")
                elif name.startswith("RealUnits"):
                    if any(row[0] <= 0 or row[1] < 0 or not 0 <= row[2] <= 1 for row in rows):
                        failures_for_file.append("pressure/gap/standard contact area outside physical range")
                elif name == "out.1u.2ProbabilityInterfacialSeparation":
                    if any(row[0] <= 0 or row[1] < 0 for row in rows):
                        failures_for_file.append("gap or density outside physical range")
                elif name.startswith("out.1pressure.2zeta"):
                    if any(row[0] <= 0 or row[1] < 1 or not 0 <= row[3] <= 1
                           or row[4] < 0 or row[5] < 0 for row in rows):
                        failures_for_file.append("scale-resolved pressure/area/gaps outside physical range")
                elif name.startswith("out.1pressure.2baruFromP"):
                    if any(row[0] <= 0 or row[1] < 0 or row[2] < 0 for row in rows):
                        failures_for_file.append("pressure or mean gaps outside physical range")
            checks[name] = {"numeric_rows": len(rows), "required_minimum": minimum,
                            "passed": not failures_for_file, "failures": failures_for_file,
                            "sha256": sha256(output / name) if (output / name).is_file() else None}
            failures.extend(f"{name}: {failure}" for failure in failures_for_file)
        except ValueError as exc:
            checks[name] = {"passed": False, "failures": [str(exc)]}
            failures.append(str(exc))
    error_files = {path.name: {"bytes": path.stat().st_size, "sha256": sha256(path)}
                   for path in sorted(output.iterdir())
                   if path.is_file() and path.name.upper().startswith("ERROR") and path.stat().st_size}
    warning_files = {path.name: {"bytes": path.stat().st_size, "sha256": sha256(path)}
                     for path in sorted(output.iterdir())
                     if path.is_file() and path.name.upper().startswith("WARNING") and path.stat().st_size}
    failures.extend(f"Nonempty error file: {name}" for name in error_files)
    return {"passed": not failures, "outputs": checks, "error_files": error_files,
            "warning_files": warning_files, "failures": failures}


def parser():
    result = argparse.ArgumentParser(description=__doc__)
    result.add_argument("--archive", required=True, help="Original complete ContactMechanics ZIP")
    result.add_argument("--psd", required=True, help="Research PSD to use, without modification")
    result.add_argument("--output", required=True, help="New flat run directory; must not exist")
    result.add_argument("--compiler", default="gfortran", help="GNU Fortran executable or wrapper")
    result.add_argument("--mode", type=int, choices=(2, 3), default=2)
    result.add_argument("--delzeta", type=float, default=0.01)
    result.add_argument("--binary", help="Reuse a caller-supplied compiled executable; records its hash")
    result.add_argument("--real-erf", action="store_true",
                        help="Opt in to replacing real DERF1 with intrinsic ERF in the build copy; records patch and changes source hash")
    result.add_argument("--timeout", type=float, default=600.0, help="Seconds per compile/run command")
    result.add_argument("--prepare-only", action="store_true", help="Extract and record inputs; do not compile or run")
    return result


def main(argv=None):
    args = parser().parse_args(argv)
    if not math.isfinite(args.delzeta) or args.delzeta <= 0:
        parser().error("--delzeta must be finite and positive")
    if not math.isfinite(args.timeout) or args.timeout <= 0:
        parser().error("--timeout must be finite and positive")
    if args.real_erf and args.binary:
        parser().error("--real-erf cannot be combined with --binary: prebuilt patch provenance cannot be established")
    manifest = {"schema_version": 1, "created_at": datetime.now(timezone.utc).isoformat(),
                "status": "preparing", "verification": "not_run"}
    output = None
    try:
        output, build, source = prepare(args, manifest)
        if args.prepare_only:
            print(f"Prepared inputs only [source variant: {manifest['source']['variant']}]; no calculation executed: {output / 'manifest.json'}")
            return 0
        binary = build / "contact_mechanics"
        if args.binary:
            supplied = Path(args.binary).resolve()
            if not supplied.is_file():
                raise ValueError("--binary must be an existing executable")
            shutil.copy2(supplied, binary)
            manifest["compiler"] = {"compiled_here": False, "binary_supplied": str(supplied),
                                    "source_correspondence": "caller-provided; not established by this script"}
        else:
            compiler = shutil.which(args.compiler)
            if compiler is None:
                raise ValueError(f"GNU Fortran compiler not found: {args.compiler}")
            version = subprocess.run([compiler, "--version"], capture_output=True, text=True,
                                     timeout=min(args.timeout, 20), check=False)
            manifest["compiler"] = {"compiled_here": True, "executable": compiler,
                                    "version": version.stdout.splitlines()[0] if version.stdout else "unknown"}
            manifest["compile"] = run_command([compiler, *COMPILE_FLAGS, source.name, "-o", binary.name],
                                              build, output / "compile.log", args.timeout)
            if manifest["compile"]["timed_out"] or manifest["compile"]["returncode"] != 0:
                manifest["status"] = "compile_failed"
                save_manifest(output, manifest)
                print(f"Compilation failed; inspect {output / 'manifest.json'}", file=sys.stderr)
                return 124 if manifest["compile"]["timed_out"] else 1
        manifest["binary"] = {"path": str(binary.relative_to(output)), "sha256": sha256(binary)}
        manifest["status"] = "running"
        save_manifest(output, manifest)
        manifest["run"] = run_command([str(binary)], output, output / "run.log", args.timeout)
        manifest["validation"] = validate_outputs(output, args.mode)
        successful = (manifest["run"]["returncode"] == 0 and not manifest["run"]["timed_out"]
                      and manifest["validation"]["passed"])
        manifest["status"] = "validated" if successful else "failed"
        manifest["verification"] = "meaningful_outputs_checked"
        manifest["source"]["sha256_after_run"] = sha256(source)
        if manifest["source"]["sha256_after_run"] != manifest["source"]["compiled_sha256"]:
            manifest["status"] = "failed"
            manifest["validation"]["failures"].append("Fortran source changed during run")
            manifest["validation"]["passed"] = False
            successful = False
        save_manifest(output, manifest)
        print(f"{manifest['status']} [source variant: {manifest['source']['variant']}]: {output / 'manifest.json'}")
        if manifest["run"]["timed_out"]:
            return 124
        return 0 if successful else 1
    except (OSError, ValueError, zipfile.BadZipFile, subprocess.SubprocessError) as exc:
        manifest["status"] = "failed"
        manifest["failure"] = str(exc)
        if output is None and "working_directory" in manifest:
            # Preparation may fail after this script exclusively created the
            # directory; record that failure without touching pre-existing dirs.
            output = Path(manifest["working_directory"])
        if output is not None:
            save_manifest(output, manifest)
        print(f"Reference run failed: {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
