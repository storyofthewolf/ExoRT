"""
verify_lib.py -- Stage D acceptance check for libexort.

Runs a regression fixture through the library API and compares every output
against the committed golden baseline for the same case (produced by
exort.exe via tests/regression/run_regression.py). Also checks that
exort_run_columns and repeated exort_run_column calls reproduce the first
result bit-for-bit (clear-sky path is deterministic).

The baseline NetCDF stores float32 (nf_real), so the file comparison
tolerance is float32-limited (~1e-7 relative); the in-memory determinism
checks are exact (max |delta| = 0). Any wiring error in the library path
shows up orders of magnitude above these thresholds.

Usage (defaults reproduce the TS300K_G2V regression case):

    python verify_lib.py
    python verify_lib.py --case 2barCO2_dry_Mars_G2V
    python verify_lib.py --list
"""

import argparse
import os
import sys

import numpy as np
import netCDF4

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from exort_api import ExoRT, load_rtprofile  # noqa: E402

_REPO_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
_REG = os.path.join(_REPO_ROOT, "tests", "regression")

# clear-sky regression cases (name -> fixture, star, scon, g); cloud/haze
# cases are excluded here: MCICA keeps the cloudy path out of scope for the
# bitwise determinism check, and clear cases already exercise every gas path
# stellar files are the _n84 variants: the committed baselines were generated
# by exort.exe on the n84 grid (see _STAR_N84 in run_regression.py)
_TS = {f"TS{t}K_{s}": (f"RTprofile_in_TS{t}K.nc",
                       "G2V_SUN_n84.nc" if s == "G2V" else "blackbody_3400K_n84.nc",
                       680.0, 9.80616)
       for t in (250, 273, 300, 320, 340, 360) for s in ("G2V", "BB3400")}
CASES = dict(_TS)
CASES["2barCO2_dry_Mars_G2V"] = ("RTprofile_in_2barCO2_dry_pver300.nc",
                                 "G2V_SUN_n84.nc", 451.166, 3.711)

# baseline variable -> (result key, scale applied to library output)
CSEC = 86400.0
COMPARE = {
    "LWUP": ("lw_upflux", 1.0), "LWDN": ("lw_dnflux", 1.0),
    "SWUP": ("sw_upflux", 1.0), "SWDN": ("sw_dnflux", 1.0),
    "LWHR": ("lw_dtdt", CSEC), "SWHR": ("sw_dtdt", CSEC),
    "FSDTOA": ("sol_toa", 1.0),
    "LWUP_SPECTRAL": ("lw_upflux_spectral", 1.0),
    "LWDN_SPECTRAL": ("lw_dnflux_spectral", 1.0),
    "SWUP_SPECTRAL": ("sw_upflux_spectral", 1.0),
    "SWDN_SPECTRAL": ("sw_dnflux_spectral", 1.0),
}

# baseline files are float32: tolerance is storage-precision, not physics
RTOL, ATOL = 1e-5, 1e-4


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    ap.add_argument("--case", default="TS300K_G2V", choices=sorted(CASES),
                    metavar="NAME", help="regression case (default TS300K_G2V)")
    ap.add_argument("--lib", default=None, help="path to libexort")
    ap.add_argument("--list", action="store_true", help="list cases and exit")
    args = ap.parse_args()

    if args.list:
        print("\n".join(sorted(CASES)))
        return 0

    fixture, star, scon, g = CASES[args.case]
    fixture = os.path.join(_REG, "fixtures", fixture)
    baseline = os.path.join(_REG, "baselines", f"RTprofile_out_{args.case}.nc")

    rt = ExoRT(args.lib)
    print(f"libexort: {rt.lib_path}")
    print(f"dims: {rt.dims}")
    rt.init(solar_file=star, scon=scon, g=g)

    state = load_rtprofile(fixture, rt.dims)
    result = rt.run_column(state)

    # --- baseline comparison (float32 storage precision) ---
    print(f"\n=== {args.case}: library vs baseline "
          f"(rtol={RTOL}, atol={ATOL}) ===")
    n_fail = 0
    with netCDF4.Dataset(baseline) as ds:
        for var, (key, scale) in COMPARE.items():
            base = np.asarray(ds.variables[var][:], dtype=np.float64)
            new = np.atleast_1d(np.asarray(result[key]) * scale)
            if base.ndim == 2:
                # file stores (nwave, pverp) [netCDF C dim order]; the
                # library returns (pverp, nwave)
                base = base.T
            base = base.reshape(new.shape)
            adiff = np.abs(new - base)
            ok = np.all(adiff <= ATOL + RTOL * np.abs(base))
            n_fail += 0 if ok else 1
            denom = np.maximum(np.abs(base), 1e-30)
            print(f"  {var:15s} {'PASS' if ok else 'FAIL':4s}  "
                  f"max|d|={adiff.max():.3e}  "
                  f"max rel={np.max(adiff / denom):.3e}")

    # --- determinism: repeat call and batch path must be bit-for-bit ---
    print("\n=== determinism (exact, max |delta| must be 0) ===")
    again = rt.run_column(state)
    batch = rt.run_columns([state, state])
    for label, other in [("repeat run_column", again),
                         ("run_columns[0]", batch[0]),
                         ("run_columns[1]", batch[1])]:
        worst = max(np.max(np.abs(np.asarray(other[k]) - np.asarray(result[k])))
                    for k, _ in COMPARE.values())
        ok = worst == 0.0
        n_fail += 0 if ok else 1
        print(f"  {label:20s} {'PASS' if ok else 'FAIL':4s}  max|d|={worst:.3e}")

    olr = result["lw_upflux"][0]
    swdn_sfc = result["sw_dnflux"][-1]
    print(f"\nLWUP_TOM={olr:.3f} W/m2   SWDN_SFC={swdn_sfc:.3f} W/m2")
    print("PASS" if n_fail == 0 else f"FAIL ({n_fail} checks failed)")
    return 0 if n_fail == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
