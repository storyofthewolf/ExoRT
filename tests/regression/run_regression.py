#!/usr/bin/env python3
"""
ExoRT 1-D regression harness.

Runs a set of standard atmospheric profiles through the ExoRT executable and
compares the output flux/heating fields against committed golden baselines with
a numerical tolerance.

Each case (see build_cases()) carries its own fixture, stellar spectrum,
insolation (shr_const_scon) and gravity (exo_g), so heterogeneous planets sit
alongside the Earth-like temperature sequence:
  - Earth-like   TS250K..TS360K x {G2V_SUN, blackbody_3400K}   (12 cases)
  - Mars-like    2barCO2_dry_Mars_G2V                          (1 case)

Per case:
  1. Write run/user_nl_exort   (solar_file, shr_const_scon, exo_g for the case)
  2. Copy fixture            -> run/RTprofile_in.nc
  3. Run the executable      -> run/RTprofile_out.nc
  4. Capture output         -> outputs/RTprofile_out_<name>.nc
  5. Compare vs baselines/   (rel/abs tolerance) across four field groups:
       - integrated fluxes  (TOA / surface scalars from the flux arrays + FSDTOA)
       - per-level fluxes    (full LWUP/LWDN/SWUP/SWDN vectors)
       - spectral fluxes     (LW/SW UP/DN _SPECTRAL arrays)
       - heating rates       (LWHR, SWHR)

Usage
-----
  python run_regression.py                  # run all cases, compare to baselines
  python run_regression.py --list           # list case names + their physics
  python run_regression.py --cases TS300K Mars   # subset by substring match
  python run_regression.py --generate-baselines
  python run_regression.py --rtol 1e-3 --atol 1e-3

Exit status is 0 only if every compared case passes (non-generate mode).
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

# ---------------------------------------------------------------------------
# Paths (this script lives in <repo>/tests/regression/)
# ---------------------------------------------------------------------------
HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, os.pardir, os.pardir))
RUN_DIR = os.path.join(REPO, "run")
FIXTURE_DIR = os.path.join(HERE, "fixtures")
BASELINE_DIR = os.path.join(HERE, "baselines")
OUTPUT_DIR = os.path.join(HERE, "outputs")

EXE_CANDIDATES = ["n68equiv.exe"]          # n68equiv is the reference build
RUN_INPUT = os.path.join(RUN_DIR, "RTprofile_in.nc")
RUN_OUTPUT = os.path.join(RUN_DIR, "RTprofile_out.nc")
NAMELIST = os.path.join(RUN_DIR, "user_nl_exort")

# ---------------------------------------------------------------------------
# Test matrix
# ---------------------------------------------------------------------------
# Each case is a dict:
#   name    unique, filename-safe case id (drives output/baseline filenames)
#   fixture input NetCDF basename in fixtures/
#   star    solar_file written into user_nl_exort
#   scon    shr_const_scon  [W m-2]
#   g       exo_g           [m s-2]
# Homogeneous suites (same physics across stars) are convenient to express as a
# product; heterogeneous cases (different gravity/insolation/fixture) are just
# appended. Both end up as plain entries in CASES.

# Earth-like temperature sequence: 6 profiles x 2 stars, fixed g and insolation.
_TS_PROFILES = ["TS250K", "TS273K", "TS300K", "TS320K", "TS340K", "TS360K"]
_TS_STARS = ["G2V_SUN_n68.nc", "blackbody_3400K_n68.nc"]
_TS_SCON = 680.0
_TS_G = 9.80616


def _ts_cases():
    cases = []
    for profile in _TS_PROFILES:
        for star in _TS_STARS:
            cases.append({
                "name": f"{profile}_{star_tag(star)}",
                "fixture": f"RTprofile_in_{profile}.nc",
                "star": star,
                "scon": _TS_SCON,
                "g": _TS_G,
            })
    return cases


def build_cases():
    cases = _ts_cases()
    # Mars-like 2 bar dry CO2 atmosphere: reduced gravity and insolation, G2V star.
    cases.append({
        "name": "2barCO2_dry_Mars_G2V",
        "fixture": "RTprofile_in_2barCO2_dry_pver300.nc",
        "star": "G2V_SUN_n68.nc",
        "scon": 451.166,
        "g": 3.711,
    })
    return cases


DEFAULT_RTOL = 1e-3
DEFAULT_ATOL = 1e-3

# Output variables grouped by category. Per-level and integrated draw from the
# same four flux arrays; we compare the full vectors (per-level) and separately
# pull TOA (index 0) and surface (index -1) scalars (integrated).
FLUX_VARS = ["LWUP", "LWDN", "SWUP", "SWDN"]
SPECTRAL_VARS = ["LWUP_SPECTRAL", "LWDN_SPECTRAL", "SWUP_SPECTRAL", "SWDN_SPECTRAL"]
HEATING_VARS = ["LWHR", "SWHR"]
INTEGRATED_SCALARS = ["FSDTOA"]

ALL_COMPARE_VARS = FLUX_VARS + SPECTRAL_VARS + HEATING_VARS + INTEGRATED_SCALARS

# First-order glance metrics shown inline in the results table.
# Vertical index convention (confirmed via PINT): index 0 = top-of-model,
# index -1 = surface.
#   LWUP[0]   outgoing longwave at top-of-model  [W m-2]
#   SWDN[-1]  downwelling shortwave at surface    [W m-2]
GLANCE_METRICS = [
    ("LWUP_TOM", "LWUP", 0),
    ("SWDN_SFC", "SWDN", -1),
]


def star_tag(star_file):
    """Short, filename-safe tag for a solar file."""
    base = os.path.splitext(star_file)[0]
    if base.startswith("G2V"):
        return "G2V"
    if base.startswith("blackbody_3400K"):
        return "BB3400"
    return base.replace(".", "_")


def output_path(case):
    return os.path.join(OUTPUT_DIR, f"RTprofile_out_{case['name']}.nc")


def baseline_path(case):
    return os.path.join(BASELINE_DIR, f"RTprofile_out_{case['name']}.nc")


def fixture_path(case):
    return os.path.join(FIXTURE_DIR, case["fixture"])


def runtime_env():
    """
    Build the environment for the executable, ensuring the NetCDF Fortran
    shared library is locatable.

    On macOS the exe is linked against @rpath/libnetcdff.*.dylib but the
    Anaconda lib dir is not on the default loader path, so we prepend the
    directory reported by `nf-config --prefix` to (DY)LD_LIBRARY_PATH. An
    explicit EXORT_LIB_DIR env var overrides the autodetection.
    """
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


def find_exe():
    for name in EXE_CANDIDATES:
        path = os.path.join(RUN_DIR, name)
        if os.path.isfile(path) and os.access(path, os.X_OK):
            return path
    return None


def write_namelist(case):
    """Write run/user_nl_exort with the case's stellar spectrum and physics."""
    text = (
        "&exort_config\n"
        f"  solar_file     = '{case['star']}',\n"
        f"  shr_const_scon = {case['scon']},\n"
        f"  exo_g          = {case['g']}\n"
        "/\n"
    )
    with open(NAMELIST, "w") as fh:
        fh.write(text)


def run_case(case, exe, env, verbose=False):
    """Run one case; return the captured output path."""
    fixture = fixture_path(case)
    if not os.path.isfile(fixture):
        raise FileNotFoundError(f"fixture missing: {fixture}")

    write_namelist(case)
    shutil.copyfile(fixture, RUN_INPUT)
    if os.path.exists(RUN_OUTPUT):
        os.remove(RUN_OUTPUT)

    proc = subprocess.run(
        [exe],
        cwd=RUN_DIR,
        env=env,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
    )
    if proc.returncode != 0 or not os.path.isfile(RUN_OUTPUT):
        sys.stderr.write(proc.stdout)
        raise RuntimeError(
            f"executable failed for {case['name']} "
            f"(rc={proc.returncode}, output written={os.path.isfile(RUN_OUTPUT)})"
        )
    if verbose:
        sys.stdout.write(proc.stdout)

    dest = output_path(case)
    shutil.copyfile(RUN_OUTPUT, dest)
    return dest


def _read_vars(path, names):
    out = {}
    with netCDF4.Dataset(path) as ds:
        for n in names:
            if n in ds.variables:
                out[n] = np.asarray(ds.variables[n][:], dtype=np.float64)
    return out


def glance_metrics(path):
    """Return {label: value} for the first-order glance metrics in a file."""
    data = _read_vars(path, [v for _, v, _ in GLANCE_METRICS])
    out = {}
    for label, var, idx in GLANCE_METRICS:
        out[label] = float(data[var][idx]) if var in data else float("nan")
    return out


def compare_case(new_path, base_path, rtol, atol):
    """
    Compare one output file against its baseline.

    Returns (passed, list_of_failure_strings, glance_dict) where glance_dict
    maps each glance-metric label to (new_value, base_value).
    """
    new_glance = glance_metrics(new_path)
    if not os.path.isfile(base_path):
        return (False,
                [f"baseline missing: {os.path.basename(base_path)} "
                 f"(run with --generate-baselines first)"],
                {label: (new_glance[label], float("nan"))
                 for label, _, _ in GLANCE_METRICS})

    base_glance = glance_metrics(base_path)
    glance = {label: (new_glance[label], base_glance[label])
              for label, _, _ in GLANCE_METRICS}

    new = _read_vars(new_path, ALL_COMPARE_VARS)
    base = _read_vars(base_path, ALL_COMPARE_VARS)

    failures = []

    def check(var, a, b):
        if a.shape != b.shape:
            failures.append(f"{var}: shape {a.shape} != baseline {b.shape}")
            return
        close = np.isclose(a, b, rtol=rtol, atol=atol, equal_nan=True)
        if not close.all():
            diff = np.abs(a - b)
            denom = np.maximum(np.abs(b), 1e-30)
            rel = diff / denom
            idx = np.unravel_index(np.argmax(np.where(close, -1.0, rel)), a.shape)
            failures.append(
                f"{var}: {(~close).sum()}/{a.size} pts exceed tol; "
                f"max |abs|={diff.max():.3e} max rel={rel.max():.3e} "
                f"at {idx} (new={a[idx]:.6e} base={b[idx]:.6e})"
            )

    # Group 1+2: flux arrays (full vectors cover per-level; baseline covers integrated
    # TOA/surface endpoints implicitly since they are elements of the same array).
    for var in FLUX_VARS + SPECTRAL_VARS + HEATING_VARS + INTEGRATED_SCALARS:
        if var not in new:
            failures.append(f"{var}: absent from run output")
            continue
        if var not in base:
            failures.append(f"{var}: absent from baseline")
            continue
        check(var, new[var], base[var])

    return (len(failures) == 0), failures, glance


def main():
    ap = argparse.ArgumentParser(description="ExoRT 1-D regression harness")
    ap.add_argument("--generate-baselines", action="store_true",
                    help="run the matrix and store outputs as new baselines")
    ap.add_argument("--cases", nargs="+", default=None,
                    help="subset of case names to run (default: all). "
                         "Substring match, e.g. --cases TS300K Mars")
    ap.add_argument("--list", action="store_true",
                    help="list all case names and exit")
    ap.add_argument("--rtol", type=float, default=DEFAULT_RTOL)
    ap.add_argument("--atol", type=float, default=DEFAULT_ATOL)
    ap.add_argument("--verbose", action="store_true",
                    help="echo executable stdout for each run")
    args = ap.parse_args()

    cases = build_cases()
    if args.list:
        for c in cases:
            print(f"  {c['name']:<22}  star={c['star']:<24} "
                  f"scon={c['scon']:<8} g={c['g']}")
        return 0
    if args.cases:
        cases = [c for c in cases
                 if any(sub in c["name"] for sub in args.cases)]
        if not cases:
            sys.exit(f"ERROR: no cases match {args.cases}")

    os.makedirs(OUTPUT_DIR, exist_ok=True)
    if args.generate_baselines:
        os.makedirs(BASELINE_DIR, exist_ok=True)

    exe = find_exe()
    if exe is None:
        sys.exit(f"ERROR: no executable found in {RUN_DIR} "
                 f"(expected one of {EXE_CANDIDATES}; build with "
                 f"`cd build && make n68equiv`).")

    env = runtime_env()

    # Preserve any existing namelist so we don't clobber the user's setup.
    saved_nl = None
    if os.path.isfile(NAMELIST):
        with open(NAMELIST) as fh:
            saved_nl = fh.read()

    results = []  # (name, status, details, glance)
    try:
        for case in cases:
            name = case["name"]
            try:
                out = run_case(case, exe, env, verbose=args.verbose)
            except Exception as exc:  # noqa: BLE001
                results.append((name, "ERROR", [str(exc)], {}))
                continue

            if args.generate_baselines:
                shutil.copyfile(out, baseline_path(case))
                g = glance_metrics(out)
                results.append((name, "BASELINE", [],
                                {lbl: (g[lbl], g[lbl])
                                 for lbl, _, _ in GLANCE_METRICS}))
            else:
                ok, fails, glance = compare_case(
                    out, baseline_path(case), args.rtol, args.atol
                )
                results.append((name, "PASS" if ok else "FAIL", fails, glance))
    finally:
        # Restore the user's original namelist.
        if saved_nl is not None:
            with open(NAMELIST, "w") as fh:
                fh.write(saved_nl)

    # ----- report -----
    print()
    if args.generate_baselines:
        print(f"=== baseline generation (rtol={args.rtol}, atol={args.atol} unused) ===")
    else:
        print(f"=== regression results (rtol={args.rtol}, atol={args.atol}) ===")

    width = max(len(n) for n, _, _, _ in results)
    labels = [lbl for lbl, _, _ in GLANCE_METRICS]
    # Metric body "{:9.3f} vs {:9.3f} (Δ{:+8.3f})" is exactly COL chars wide;
    # header labels and data cells share this width so columns line up.
    COL = 34

    def fmt_metric(pair):
        new, base = pair
        if base != base:  # NaN baseline
            return f"{new:9.3f} (no base)".ljust(COL)
        return f"{new:9.3f} vs {base:9.3f} (Δ{new - base:+8.3f})".ljust(COL)

    prefix = f"  {{:<{width}}}  {{:<8}}"  # case + status block, shared by all rows
    header = prefix.format("case", "")
    for lbl in labels:
        header += "  " + f"{lbl} [W/m2]".ljust(COL)
    print(header)

    n_fail = 0
    for name, status, details, glance in results:
        row = prefix.format(name, status)
        for lbl in labels:
            row += "  " + (fmt_metric(glance[lbl]) if lbl in glance else " " * COL)
        print(row)
        if status in ("FAIL", "ERROR"):
            n_fail += 1
            for line in details:
                print(f"      - {line}")

    total = len(results)
    print()
    if args.generate_baselines:
        print(f"Generated {total} baselines in {BASELINE_DIR}")
        return 0

    passed = total - n_fail
    print(f"{passed}/{total} cases passed.")
    return 0 if n_fail == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
