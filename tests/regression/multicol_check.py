#!/usr/bin/env python3
"""
Stage E1 acceptance check: a multi-column batch run must reproduce the
corresponding single-column runs EXACTLY (max |delta| == 0 on every output
variable), since batch and single columns execute the identical per-column
code path (run_one_column) serially.

Flow:
  1. run three single-column cases that share runtime config (star/scon/g):
     TS250K_G2V, TS300K_G2V, TS340K_G2V
  2. stack their fixtures with tools/stackColumns.py -> ncol=3 input
  3. run the batch through run/exort.exe with OMP_NUM_THREADS=1, then again
     with OMP_NUM_THREADS=8 (Stage E2 OpenMP check; default thread stacks —
     the big solver arrays are heap-allocated)
  4. compare every variable of every batch column against the single runs,
     and the threaded batch against the serial batch (both exact)

Reuses the regression harness helpers (env, namelist, run_case). Assumes
run/exort.exe is already built (run run_regression.py first, or make exort).

  USER_FC=gfortran python multicol_check.py
"""

import os
import shutil
import subprocess
import sys

import numpy as np
import netCDF4

import run_regression as rr

CASE_NAMES = ["TS250K_G2V", "TS300K_G2V", "TS340K_G2V"]
STACKER = os.path.join(rr.REPO, "tools", "stackColumns.py")
BATCH_INPUT = os.path.join(rr.RUN_DIR, "RTprofile_in.nc")
BATCH_CAPTURE = os.path.join(rr.OUTPUT_DIR, "RTprofile_out_multicol3.nc")
BATCH_CAPTURE_OMP = os.path.join(rr.OUTPUT_DIR, "RTprofile_out_multicol3_omp8.nc")
OMP_THREADS = "8"


def run_batch(exe, env, capture, tag):
    if os.path.exists(rr.RUN_OUTPUT):
        os.remove(rr.RUN_OUTPUT)
    print(f"-- batch: {tag}")
    proc = subprocess.run(
        [exe], cwd=rr.RUN_DIR, env=env,
        stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True,
    )
    if proc.returncode != 0 or not os.path.isfile(rr.RUN_OUTPUT):
        sys.stderr.write(proc.stdout)
        sys.exit(f"ERROR: batch run ({tag}) failed")
    shutil.copyfile(rr.RUN_OUTPUT, capture)


def main():
    cases = {c["name"]: c for c in rr.build_cases()}
    picked = [cases[n] for n in CASE_NAMES]

    exe = os.path.join(rr.RUN_DIR, "exort.exe")
    if not os.path.isfile(exe):
        sys.exit("ERROR: run/exort.exe not built (cd build && make exort)")
    env = rr.runtime_env()
    os.makedirs(rr.OUTPUT_DIR, exist_ok=True)

    # preserve any pre-existing namelist, as the harness does
    saved_nl = None
    if os.path.exists(rr.NAMELIST):
        with open(rr.NAMELIST) as fh:
            saved_nl = fh.read()

    try:
        # 1. singles
        single_paths = []
        for case in picked:
            print(f"-- single: {case['name']}")
            single_paths.append(rr.run_case(case, exe, env))

        # 2. stacked input (all three share star/scon/g; namelist from case 0)
        fixtures = [rr.fixture_path(c) for c in picked]
        subprocess.run(
            [sys.executable, STACKER, *fixtures, "-o", BATCH_INPUT],
            check=True,
        )

        # 3. batch runs: serial reference, then multi-threaded (Stage E2)
        rr.write_namelist(picked[0])
        run_batch(exe, dict(env, OMP_NUM_THREADS="1"),
                  BATCH_CAPTURE, f"ncol={len(picked)}, OMP_NUM_THREADS=1")
        run_batch(exe, dict(env, OMP_NUM_THREADS=OMP_THREADS),
                  BATCH_CAPTURE_OMP,
                  f"ncol={len(picked)}, OMP_NUM_THREADS={OMP_THREADS}")
    finally:
        if saved_nl is not None:
            with open(rr.NAMELIST, "w") as fh:
                fh.write(saved_nl)
        elif os.path.exists(rr.NAMELIST):
            os.remove(rr.NAMELIST)

    # 4. compare
    batch = netCDF4.Dataset(BATCH_CAPTURE)
    ncol = len(batch.dimensions["ncol"])
    assert ncol == len(picked), f"batch ncol={ncol} != {len(picked)}"

    n_fail = 0
    print(f"\n== batch vs singles (exact, max |delta| must be 0) ==")
    for ic, (case, spath) in enumerate(zip(picked, single_paths)):
        single = netCDF4.Dataset(spath)
        worst = 0.0
        worst_var = ""
        for name, bvar in batch.variables.items():
            svar = single.variables[name]
            b = np.asarray(bvar[:])
            s = np.asarray(svar[:]).squeeze()
            col = b[ic] if "ncol" in bvar.dimensions else b.squeeze()
            col = np.asarray(col).squeeze()
            if col.shape != s.shape:
                print(f"  col {ic} {name}: SHAPE MISMATCH {col.shape} vs {s.shape}")
                n_fail += 1
                continue
            d = float(np.max(np.abs(col - s))) if col.size else 0.0
            if d > worst:
                worst, worst_var = d, name
            if d != 0.0:
                print(f"  col {ic} ({case['name']}) {name}: max|d|={d:.3e}  FAIL")
                n_fail += 1
        status = "PASS" if worst == 0.0 else f"FAIL (worst {worst_var})"
        print(f"  col {ic} ({case['name']}): max|d|={worst:.3e}  {status}")
        single.close()

    # 5. threaded batch vs serial batch (exact)
    print(f"\n== batch OMP_NUM_THREADS={OMP_THREADS} vs 1 (exact) ==")
    omp = netCDF4.Dataset(BATCH_CAPTURE_OMP)
    worst, worst_var = 0.0, ""
    for name, bvar in batch.variables.items():
        d = float(np.max(np.abs(np.asarray(omp.variables[name][:], dtype="f8")
                                - np.asarray(bvar[:], dtype="f8"))))
        if d > worst:
            worst, worst_var = d, name
        if d != 0.0:
            print(f"  {name}: max|d|={d:.3e}  FAIL")
            n_fail += 1
    status = "PASS" if worst == 0.0 else f"FAIL (worst {worst_var})"
    print(f"  threaded vs serial: max|d|={worst:.3e}  {status}")
    omp.close()
    batch.close()

    if n_fail:
        sys.exit(f"FAIL: {n_fail} mismatching variables")
    print("PASS")


if __name__ == "__main__":
    main()
