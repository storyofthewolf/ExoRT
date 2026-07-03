#!/usr/bin/env python3
"""Per-column config acceptance check (Stage E2b: per-column grav/scon).

A batch whose columns carry their own 'grav'/'scon' input variables must
reproduce single-column runs configured with the same values through the
runtime namelist (exo_g / shr_const_scon):

  col A: no per-column vars           == single @ namelist (680, 9.80616)
  col B: grav=3.711                   == single @ namelist (680, 3.711)
  col C: grav=3.711, scon=451.166     == single @ namelist (451.166, 3.711)

Columns A and B must match EXACTLY (max |delta| == 0): the per-column value
follows the identical arithmetic path as the namelist value. Column C is
gated at rel <= 1e-10 per variable: a namelist scon is baked into the
spectrum normalization at init, while a per-column scon rides the msdist
flux factor — same physics, different rounding order.

All three clones share the same P/T/zint grid (cloned from the TS300K
fixture), exactly as the namelist comparison runs do, so the comparison
isolates the config plumbing.

Assumes run/exort.exe is already built.

  USER_FC=gfortran python percol_config_check.py
"""

import os
import shutil
import subprocess
import sys

import numpy as np
import netCDF4

import run_regression as rr

FIXTURE = os.path.join(rr.REPO, "tests", "regression", "fixtures",
                       "RTprofile_in_TS300K.nc")
STACKER = os.path.join(rr.REPO, "tools", "stackColumns.py")
BATCH_INPUT = os.path.join(rr.RUN_DIR, "RTprofile_in.nc")

EARTH_G, MARS_G = 9.80616, 3.711
SUN_SCON, MARS_SCON = 680.0, 451.166
RTOL_SCON = 1e-10

NL = """&exort_config
  solar_file     = 'G2V_SUN_n84.nc',
  shr_const_scon = {scon},
  exo_g          = {g},
/
"""

# (tag, per-column vars for the batch clone, single-run namelist (scon, g))
COLUMNS = [
    ("default",   {},                                       (SUN_SCON, EARTH_G)),
    ("marsG",     {"grav": MARS_G},                         (SUN_SCON, MARS_G)),
    ("marsGscon", {"grav": MARS_G, "scon": MARS_SCON},      (MARS_SCON, MARS_G)),
]


def clone_with_vars(dst, pervars):
    shutil.copyfile(FIXTURE, dst)
    if not pervars:
        return
    with netCDF4.Dataset(dst, "a") as d:
        for name, val in pervars.items():
            v = d.createVariable(name, "f8", ())
            v[...] = val


def run_exe(exe, env, tag):
    if os.path.exists(rr.RUN_OUTPUT):
        os.remove(rr.RUN_OUTPUT)
    proc = subprocess.run(
        [exe], cwd=rr.RUN_DIR, env=env,
        stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True,
    )
    if proc.returncode != 0 or not os.path.isfile(rr.RUN_OUTPUT):
        sys.stderr.write(proc.stdout)
        sys.exit(f"ERROR: run failed ({tag})")
    dst = os.path.join(rr.OUTPUT_DIR, f"RTprofile_out_percol_{tag}.nc")
    shutil.copyfile(rr.RUN_OUTPUT, dst)
    return dst


def main():
    exe = os.path.join(rr.RUN_DIR, "exort.exe")
    if not os.path.isfile(exe):
        sys.exit("ERROR: run/exort.exe not built (cd build && make exort)")
    env = rr.runtime_env()
    os.makedirs(rr.OUTPUT_DIR, exist_ok=True)

    saved_nl = None
    if os.path.exists(rr.NAMELIST):
        with open(rr.NAMELIST) as fh:
            saved_nl = fh.read()

    try:
        # 1. singles: plain fixture, config through the namelist
        single_paths = []
        for tag, _, (scon, g) in COLUMNS:
            print(f"-- single: {tag} (namelist scon={scon}, g={g})")
            shutil.copyfile(FIXTURE, BATCH_INPUT)
            with open(rr.NAMELIST, "w") as fh:
                fh.write(NL.format(scon=scon, g=g))
            single_paths.append(run_exe(exe, env, f"single_{tag}"))

        # 2. batch: per-column vars in the input, default namelist
        clones = []
        for tag, pervars, _ in COLUMNS:
            path = os.path.join(rr.OUTPUT_DIR, f"RTprofile_in_percol_{tag}.nc")
            clone_with_vars(path, pervars)
            clones.append(path)
        subprocess.run([sys.executable, STACKER, *clones, "-o", BATCH_INPUT],
                       check=True)
        print(f"-- batch: ncol={len(COLUMNS)} (per-column grav/scon vars)")
        with open(rr.NAMELIST, "w") as fh:
            fh.write(NL.format(scon=SUN_SCON, g=EARTH_G))
        batch_path = run_exe(exe, env, "batch")
    finally:
        if saved_nl is not None:
            with open(rr.NAMELIST, "w") as fh:
                fh.write(saved_nl)
        elif os.path.exists(rr.NAMELIST):
            os.remove(rr.NAMELIST)

    # 3. compare
    n_fail = 0
    batch = netCDF4.Dataset(batch_path)
    print("\n== batch (per-column config) vs singles (namelist config) ==")
    for ic, ((tag, pervars, _), spath) in enumerate(zip(COLUMNS, single_paths)):
        exact = "scon" not in pervars  # scon path differs in rounding order
        single = netCDF4.Dataset(spath)
        worst, worst_var = 0.0, ""
        for name, bvar in batch.variables.items():
            b = np.asarray(bvar[:], dtype="f8")
            s = np.asarray(single.variables[name][:], dtype="f8").squeeze()
            col = np.asarray(b[ic] if "ncol" in bvar.dimensions
                             else b.squeeze()).squeeze()
            if col.shape != s.shape:
                print(f"  col {ic} {name}: SHAPE MISMATCH {col.shape} vs {s.shape}")
                n_fail += 1
                continue
            d = float(np.max(np.abs(col - s))) if col.size else 0.0
            if not exact:
                d = d / max(float(np.max(np.abs(s))) if s.size else 0.0, 1e-30)
            if d > worst:
                worst, worst_var = d, name
        single.close()
        if exact:
            ok = worst == 0.0
            crit = "max|d|"
        else:
            ok = worst <= RTOL_SCON
            crit = "max rel"
        n_fail += 0 if ok else 1
        print(f"  col {ic} ({tag}): {crit}={worst:.3e}  "
              f"{'PASS' if ok else f'FAIL (worst {worst_var})'}")
    batch.close()

    if n_fail:
        sys.exit(f"FAIL: {n_fail} column(s) mismatched")
    print("PASS")


if __name__ == "__main__":
    main()
