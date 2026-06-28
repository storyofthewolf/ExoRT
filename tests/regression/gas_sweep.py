#!/usr/bin/env python3
"""
Per-gas k-coefficient sweep: HITRAN-2016 vs HITRAN-2024, grid 68 vs 84.

Purpose (diagnostic, NOT pass/fail vs golden baselines): for each gas whose
HITRAN line list changed between the legacy build and the v2 src.exort build,
run a single-gas fixture (that gas non-zero, N2 background) through three exes
that each encode a known vintage/grid, and tabulate the outgoing longwave
(OLR = LWUP at top-of-model) so two effects separate cleanly:

  n68equiv.exe : HITRAN-2016, 68-band grid   (legacy reference)
  n84equiv.exe : HITRAN-2016, 84-band grid   (legacy reference, wider grid)
  exort.exe    : HITRAN-2024, 84-band grid   (v2 target)

  grid effect      = n84 - n68     (both HITRAN-2016; should be small)
  line-list effect = exort - n84   (HITRAN-2016 -> 2024, same 84 grid)

A large line-list effect flags a suspect HITRAN-2024 k-file (this is exactly how
the 2-bar-CO2 far-wing artifact was caught: n84 h16 -> 92.7, exort h24 -> 47.8).

The three exes carry their vintage via their own (hardcoded) kabs.F90 — no
string swapping. n68equiv/n84equiv kabs are HITRAN-2016; src.exort kabs is the
committed HITRAN-2024. This script BUILDS all three (honoring USER_FC), then runs
each gas fixture. It writes run/user_nl_exort per case and restores it after.

Stars: n68equiv uses the *_n68 stellar file; n84equiv/exort use *_n84. OLR is
LW-only and (in this offline mode) independent of the stellar file, so the star
choice does not bias the OLR comparison — it only needs to load on each grid.

Usage
-----
  python gas_sweep.py                 # build all three, run all gas fixtures
  python gas_sweep.py --gases CH4 NH3 # subset (substring match on fixture tag)
  python gas_sweep.py --no-build      # skip rebuild, use existing run/*.exe
"""

import argparse
import os
import shutil
import subprocess
import sys

import numpy as np

try:
    import netCDF4
except ImportError:
    sys.exit("ERROR: netCDF4 Python package is required (pip install netCDF4).")

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, os.pardir, os.pardir))
RUN_DIR = os.path.join(REPO, "run")
BUILD_DIR = os.path.join(REPO, "build")
FIXTURE_DIR = os.path.join(HERE, "fixtures")
RUN_INPUT = os.path.join(RUN_DIR, "RTprofile_in.nc")
RUN_OUTPUT = os.path.join(RUN_DIR, "RTprofile_out.nc")
NAMELIST = os.path.join(RUN_DIR, "user_nl_exort")

# exe -> (make target, vintage label, stellar file). Order = table column order.
EXES = [
    ("n68equiv.exe", "n68equiv", "HITRAN-2016 / 68-grid", "G2V_SUN_n68.nc"),
    ("n84equiv.exe", "n84equiv", "HITRAN-2016 / 84-grid", "G2V_SUN_n84.nc"),
    ("exort.exe",    "exort",    "HITRAN-2024 / 84-grid", "G2V_SUN_n84.nc"),
]

# Optional 4th column (--with-n84h24): the n84equiv BUNDLE built against the
# HITRAN-2024 k-files (its native h2o/co2/ch4/c2h6 strings swapped h16->h24, same
# as exort reads). Same data (h24) + same 84 grid as exort, but the LEGACY n84
# code path -> n84h24 vs exort isolates pure CODE differences (merge/reader/wiring)
# between the old bundle and the new src.exort refactor, with the line list held
# fixed. Built by string-swapping src.n84equiv/kabs.F90, compiling, copying the
# artifact aside, then restoring kabs and rebuilding plain n84equiv (h16).
N84H24 = ("n84equiv_h24.exe", "n84equiv_h24", "HITRAN-2024 / 84-grid", "G2V_SUN_n84.nc")
N84_KABS = os.path.join(REPO, "source", "src.n84equiv", "kabs.F90")
# The 4 native gases whose n84equiv strings are h16 and must be swapped to h24
# (nh3/co/o2/o3 are already h24/h20 in n84equiv and are left untouched).
N84_SWAP_GASES = ("h2o", "co2", "ch4", "c2h6")

# Single-gas fixtures (built by makeColumn.py; that gas non-zero, N2 background).
# Each entry: tag -> fixture basename. realistic + elevated per gas.
GAS_FIXTURES = [
    ("H2O_realistic",  "RTprofile_in_H2O_realistic.nc"),
    ("H2O_elevated",   "RTprofile_in_H2O_elevated.nc"),
    ("CO2_realistic",  "RTprofile_in_CO2_realistic.nc"),
    ("CO2_elevated",   "RTprofile_in_CO2_elevated.nc"),
    ("CH4_realistic",  "RTprofile_in_CH4_realistic.nc"),
    ("CH4_elevated",   "RTprofile_in_CH4_elevated.nc"),
    ("C2H6_realistic", "RTprofile_in_C2H6_realistic.nc"),
    ("C2H6_elevated",  "RTprofile_in_C2H6_elevated.nc"),
    ("NH3_realistic",  "RTprofile_in_NH3_realistic.nc"),
    ("NH3_elevated",   "RTprofile_in_NH3_elevated.nc"),
    ("CO_realistic",   "RTprofile_in_CO_realistic.nc"),
    ("CO_elevated",    "RTprofile_in_CO_elevated.nc"),
]

# Insolation/gravity: Earth-like, fixed across the sweep (LW is insensitive).
SCON = 680.0
GRAV = 9.80616


def runtime_env():
    env = os.environ.copy()
    lib_dir = env.get("EXORT_LIB_DIR")
    if not lib_dir:
        try:
            prefix = subprocess.run(
                ["nf-config", "--prefix"],
                stdout=subprocess.PIPE, stderr=subprocess.DEVNULL, text=True,
            ).stdout.strip()
            if prefix:
                lib_dir = os.path.join(prefix, "lib")
        except FileNotFoundError:
            lib_dir = None
    if lib_dir and os.path.isdir(lib_dir):
        for key in ("DYLD_LIBRARY_PATH", "LD_LIBRARY_PATH"):
            cur = env.get(key, "")
            env[key] = lib_dir + (os.pathsep + cur if cur else "")
    return env


def build(target, env):
    cmd = ["make", target]
    if env.get("USER_FC"):
        cmd.append(f"USER_FC={env['USER_FC']}")
    proc = subprocess.run(cmd, cwd=BUILD_DIR, env=env,
                          stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
    exe = os.path.join(RUN_DIR, target + ".exe")
    if proc.returncode != 0 or not os.path.isfile(exe):
        sys.stderr.write(proc.stdout)
        raise RuntimeError(f"`make {target}` failed (rc={proc.returncode})")


def build_n84h24(env):
    """Build the n84equiv bundle against HITRAN-2024 native-gas k-files.

    String-swaps the 4 native-gas filenames in src.n84equiv/kabs.F90 from
    hitran16 -> hitran24, builds, copies run/n84equiv.exe -> run/n84equiv_h24.exe,
    then ALWAYS restores kabs.F90 and rebuilds plain n84equiv (h16) so the regular
    n84 column is unaffected. Returns the path to the h24 exe.
    """
    with open(N84_KABS) as fh:
        original = fh.read()
    swapped = original
    for gas in N84_SWAP_GASES:
        old = f"n84_8gpt_{gas}_hitran16_"
        new = f"n84_8gpt_{gas}_hitran24_"
        if old not in swapped:
            raise RuntimeError(f"n84h24: expected '{old}' in {N84_KABS}; "
                               "kabs.F90 layout changed, aborting swap")
        swapped = swapped.replace(old, new)
    h24_exe = os.path.join(RUN_DIR, N84H24[0])
    try:
        with open(N84_KABS, "w") as fh:
            fh.write(swapped)
        build("n84equiv", env)                       # produces run/n84equiv.exe (h24)
        shutil.copyfile(os.path.join(RUN_DIR, "n84equiv.exe"), h24_exe)
        os.chmod(h24_exe, 0o755)
    finally:
        with open(N84_KABS, "w") as fh:
            fh.write(original)                        # restore h16 source
        build("n84equiv", env)                        # rebuild plain n84 (h16)
    return h24_exe


def write_namelist(star):
    with open(NAMELIST, "w") as fh:
        fh.write("&exort_config\n"
                 f"  solar_file     = '{star}',\n"
                 f"  shr_const_scon = {SCON},\n"
                 f"  exo_g          = {GRAV}\n"
                 "/\n")


def run_one(exe_path, fixture, star, env, save_as=None):
    """Run one (exe, fixture); return OLR = LWUP[0], or raise.

    If save_as is given, copy the resulting RTprofile_out.nc there (so the full
    spectral output is preserved per (gas, version) for plotting/analysis).
    """
    write_namelist(star)
    shutil.copyfile(fixture, RUN_INPUT)
    if os.path.exists(RUN_OUTPUT):
        os.remove(RUN_OUTPUT)
    proc = subprocess.run([exe_path], cwd=RUN_DIR, env=env,
                          stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
    if proc.returncode != 0 or not os.path.isfile(RUN_OUTPUT):
        sys.stderr.write(proc.stdout[-2000:])
        raise RuntimeError(f"run failed: {os.path.basename(exe_path)} / "
                           f"{os.path.basename(fixture)} (rc={proc.returncode})")
    if save_as is not None:
        shutil.copyfile(RUN_OUTPUT, save_as)
    with netCDF4.Dataset(RUN_OUTPUT) as ds:
        return float(np.asarray(ds.variables["LWUP"][:])[0])


def main():
    ap = argparse.ArgumentParser(description="Per-gas HITRAN-2016 vs 2024 OLR sweep")
    ap.add_argument("--gases", nargs="+", default=None,
                    help="subset of fixture tags (substring match), e.g. CH4 NH3")
    ap.add_argument("--no-build", action="store_true",
                    help="skip rebuild; use existing run/*.exe")
    ap.add_argument("--with-n84h24", action="store_true",
                    help="add a 4th column: the n84equiv BUNDLE built against the "
                         "HITRAN-2024 native-gas k-files (same data+grid as exort, "
                         "legacy code path). n84h24 vs exort isolates code-only "
                         "differences. With --no-build this reuses an existing "
                         "run/n84equiv_h24.exe; otherwise it is built (and plain "
                         "n84equiv is restored to h16 afterward).")
    ap.add_argument("--save-spectra", metavar="DIR", nargs="?",
                    const=os.path.join(HERE, "spectra_out"), default=None,
                    help="save per-(gas,version) RTprofile_out.nc into DIR "
                         "(default: tests/regression/spectra_out) for spectral "
                         "analysis with plot_gas_sweep_spectra.py")
    args = ap.parse_args()

    out_dir = None
    if args.save_spectra is not None:
        out_dir = os.path.abspath(args.save_spectra)
        os.makedirs(out_dir, exist_ok=True)
        print(f"=== saving spectral output to {out_dir} ===")

    fixtures = GAS_FIXTURES
    if args.gases:
        fixtures = [(t, f) for (t, f) in fixtures
                    if any(g.lower() in t.lower() for g in args.gases)]
        if not fixtures:
            sys.exit(f"ERROR: no fixtures match {args.gases}")

    env = runtime_env()

    # Active column set: base three, optionally + n84h24 (inserted before exort so
    # the table reads n68 | n84 | n84h24 | exort = grid | line-list | code effects).
    exes = list(EXES)
    if args.with_n84h24:
        exes.insert(2, N84H24)

    if not args.no_build:
        for exe_name, target, _, _ in EXES:
            print(f"=== building {target} ===")
            build(target, env)
        if args.with_n84h24:
            print("=== building n84equiv_h24 (n84 bundle + HITRAN-2024 k-files) ===")
            build_n84h24(env)
    elif args.with_n84h24 and not os.path.isfile(os.path.join(RUN_DIR, N84H24[0])):
        print("=== building n84equiv_h24 (missing; --no-build can't skip it) ===")
        build_n84h24(env)

    # Preserve the user's namelist.
    saved_nl = None
    if os.path.isfile(NAMELIST):
        with open(NAMELIST) as fh:
            saved_nl = fh.read()

    # results[tag][label] = OLR
    results = {}
    try:
        for tag, fbase in fixtures:
            fixture = os.path.join(FIXTURE_DIR, fbase)
            if not os.path.isfile(fixture):
                print(f"  SKIP {tag}: fixture missing ({fbase})")
                continue
            results[tag] = {}
            for exe_name, label, _, star in exes:
                exe_path = os.path.join(RUN_DIR, exe_name)
                if not os.path.isfile(exe_path):
                    results[tag][label] = float("nan")
                    continue
                save_as = None
                if out_dir is not None:
                    # exe_name like "n68equiv.exe" -> "n68equiv"
                    ver = os.path.splitext(exe_name)[0]
                    save_as = os.path.join(out_dir, f"out_{tag}_{ver}.nc")
                try:
                    results[tag][label] = run_one(exe_path, fixture, star, env,
                                                  save_as=save_as)
                except Exception as exc:  # noqa: BLE001
                    sys.stderr.write(f"  {tag}/{label}: {exc}\n")
                    results[tag][label] = float("nan")
    finally:
        if saved_nl is not None:
            with open(NAMELIST, "w") as fh:
                fh.write(saved_nl)

    # ----- report -----
    n68 = "n68equiv"; n84 = "n84equiv"; n84h = "n84equiv_h24"; exo = "exort"
    print()
    print("=== per-gas OLR sweep (LWUP top-of-model, W/m2) ===")
    if args.with_n84h24:
        print("  grid = n84-n68 (h16, grid only) | line-list = n84h24-n84 (h16->h24, "
              "same code) | code = exort-n84h24 (h24, code only)")
        print()
        hdr = (f"  {'fixture':<16} {n68:>10} {n84:>10} {n84h:>13} {exo:>10}   "
               f"{'grid Δ':>9} {'lineΔ':>9} {'codeΔ':>9}")
        print(hdr)
        print("  " + "-" * (len(hdr) - 2))
        for tag, _ in fixtures:
            if tag not in results:
                continue
            r = results[tag]
            v68 = r.get(n68, float("nan")); v84 = r.get(n84, float("nan"))
            vh = r.get(n84h, float("nan")); vex = r.get(exo, float("nan"))
            grid = v84 - v68
            ll = vh - v84
            code = vex - vh
            flag = ""
            if code == code and abs(code) > 0.01:  # any nonzero code path delta
                flag = "  <== CODE DIFF"
            elif ll == ll and abs(ll) > 1.0:
                flag = "  <== check"
            print(f"  {tag:<16} {v68:>10.3f} {v84:>10.3f} {vh:>13.3f} {vex:>10.3f}   "
                  f"{grid:>+9.3f} {ll:>+9.3f} {code:>+9.3f}{flag}")
    else:
        print("  grid effect = n84-n68 (both HITRAN-2016) | line-list = exort-n84 (h16->h24)")
        print()
        hdr = f"  {'fixture':<16} {n68:>10} {n84:>10} {exo:>10}   {'grid Δ':>9} {'linelist Δ':>11}"
        print(hdr)
        print("  " + "-" * (len(hdr) - 2))
        for tag, _ in fixtures:
            if tag not in results:
                continue
            r = results[tag]
            v68, v84, vex = r.get(n68, float("nan")), r.get(n84, float("nan")), r.get(exo, float("nan"))
            grid = v84 - v68
            ll = vex - v84
            flag = ""
            if ll == ll and abs(ll) > 1.0:  # >1 W/m2 line-list shift worth a look
                flag = "  <== check"
            print(f"  {tag:<16} {v68:>10.3f} {v84:>10.3f} {vex:>10.3f}   "
                  f"{grid:>+9.3f} {ll:>+11.3f}{flag}")
    print()
    print("  vintages: " + " | ".join(f"{lbl}" for _, _, lbl, _ in exes))
    return 0


if __name__ == "__main__":
    sys.exit(main())
