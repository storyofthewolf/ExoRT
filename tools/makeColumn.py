"""
makeColumn.py  --  build RTprofile_in.nc for ExoRT 1-D standalone runs.

Usage:
  python makeColumn.py              # uses defaults below
  python makeColumn.py --help

Edit the USER SETTINGS block below to set the output filename, profile,
mixing ratios, gravity, and surface albedos.  Gas mixing ratios not set
(or set to 0.0) are simply omitted from the output file; ExoRT will read
them as zero when input.F90 finds them absent.
"""

import argparse
import sys
import numpy as np
import netCDF4 as nc

from profile_data import get_profile, AVAILABLE_PROFILES

# ===========================================================================
# USER SETTINGS -- edit these for each run
# ===========================================================================

OUTPUT_FILE = "RTprofile_in.nc"
PROFILE_TAG = "TS273K"          # one of: smart_2bar_t250, TS273K, US1976

# Dry volume mixing ratios (VMR).  N2 fills the remainder automatically.
# Set a species to 0.0 to exclude it from the output file entirely.
CO2VMR  = 400e-6    # 400 ppm
CH4VMR  = 0.0
C2H6VMR = 0.0
NH3VMR  = 0.0
COVMR   = 0.0
H2VMR   = 0.0
O2VMR   = 0.0
O3VMR   = 0.0
# N2VMR is derived: 1 - sum of above

# Surface albedos (shortwave direct/diffuse, longwave direct/diffuse)
ASDIR = 0.25
ASDIF = 0.25
ALDIR = 0.25
ALDIF = 0.25

# Cosine of solar zenith angle (only relevant for shortwave runs)
COSZRS = 0.5

# Surface gravity (m/s^2) -- used only for hypsometric height calculation
GRAV = 9.80616   # Earth
# GRAV = 3.7     # Mars

# Gas constant for dry air (J/kg/K) -- used for height calculation
R_DRY = 287.04   # Earth N2/O2
# R_DRY = 188.92 # Mars CO2

# ===========================================================================
# Molecular weights (g/mol) and specific heats (J/kg/K)
# ===========================================================================
MW = {
    "n2":   28.014,
    "co2":  44.010,
    "ch4":  16.043,
    "c2h6": 30.069,
    "nh3":  17.031,
    "co":   28.010,
    "h2":    2.016,
    "o2":   32.000,
    "o3":   48.000,
    "h2o":  18.015,
}

CP = {
    "n2":   1039.0,
    "co2":   751.0,
    "ch4":  2226.0,
    "c2h6": 1750.0,
    "nh3":  2175.0,
    "co":   1040.0,
    "h2":  14320.0,
    "o2":    919.0,
    "o3":    817.0,
}


def make_column(output_file=OUTPUT_FILE, profile_tag=PROFILE_TAG,
                co2vmr=CO2VMR, ch4vmr=CH4VMR, c2h6vmr=C2H6VMR,
                nh3vmr=NH3VMR, covmr=COVMR, h2vmr=H2VMR,
                o2vmr=O2VMR, o3vmr=O3VMR,
                asdir=ASDIR, asdif=ASDIF, aldir=ALDIR, aldif=ALDIF,
                coszrs=COSZRS, grav=GRAV, r_dry=R_DRY,
                zero_h2o=False):

    prof = get_profile(profile_tag)
    pint = prof["pint"]
    pmid = prof["pmid"]
    tint = prof["tint"]
    tmid = prof["tmid"]
    q    = np.zeros_like(prof["q"]) if zero_h2o else prof["q"].copy()

    nilev = len(pint)
    nlev  = len(pmid)
    assert nilev == nlev + 1, "Profile arrays have inconsistent lengths."

    ts = tint[-1]
    ps = pint[-1]

    # Derived N2 VMR (remainder after all other dry gases)
    n2vmr = 1.0 - co2vmr - ch4vmr - c2h6vmr - nh3vmr - covmr - h2vmr - o2vmr - o3vmr
    if n2vmr < 0.0:
        sys.exit(f"ERROR: dry VMRs sum to > 1 (n2vmr = {n2vmr:.4f})")

    # Dry molecular weight and specific heat
    mwdry = (n2vmr  * MW["n2"]   + co2vmr  * MW["co2"]  +
             ch4vmr * MW["ch4"]  + c2h6vmr * MW["c2h6"] +
             nh3vmr * MW["nh3"]  + covmr   * MW["co"]   +
             h2vmr  * MW["h2"]   + o2vmr   * MW["o2"]   +
             o3vmr  * MW["o3"])

    # Mass mixing ratios (dry)
    n2mmr   = n2vmr   * MW["n2"]   / mwdry
    co2mmr  = co2vmr  * MW["co2"]  / mwdry
    ch4mmr  = ch4vmr  * MW["ch4"]  / mwdry
    c2h6mmr = c2h6vmr * MW["c2h6"] / mwdry
    nh3mmr  = nh3vmr  * MW["nh3"]  / mwdry
    commr   = covmr   * MW["co"]   / mwdry
    h2mmr   = h2vmr   * MW["h2"]   / mwdry
    o2mmr   = o2vmr   * MW["o2"]   / mwdry
    o3mmr   = o3vmr   * MW["o3"]   / mwdry

    cpdry = (n2mmr  * CP["n2"]   + co2mmr  * CP["co2"]  +
             ch4mmr * CP["ch4"]  + c2h6mmr * CP["c2h6"] +
             nh3mmr * CP["nh3"]  + commr   * CP["co"]   +
             h2mmr  * CP["h2"]   + o2mmr   * CP["o2"]   +
             o3mmr  * CP["o3"])

    # Layer thickness (Pa)
    pdel = pint[1:] - pint[:-1]

    # Interface heights via hypsometric equation (surface = 0)
    zint = np.zeros(nilev)
    for k in range(nilev - 2, -1, -1):
        zint[k] = zint[k + 1] + r_dry * tmid[k] / grav * np.log(pint[k + 1] / pint[k])

    # Print summary
    print("----------------------------------------------")
    print("makeColumn.py")
    print("----------------------------------------------")
    print(f"output file  : {output_file}")
    print(f"profile tag  : {profile_tag}")
    print(f"Ts           : {ts:.2f} K")
    print(f"Ps           : {ps:.2f} Pa")
    print(f"levels       : {nlev}")
    print(f"mwdry        : {mwdry:.4f} g/mol")
    print(f"cpdry        : {cpdry:.2f} J/kg/K")
    print("-- dry VMRs --")
    print(f"  n2vmr  = {n2vmr:.6g}")
    print(f"  co2vmr = {co2vmr:.6g}")
    if ch4vmr  > 0: print(f"  ch4vmr = {ch4vmr:.6g}")
    if c2h6vmr > 0: print(f"  c2h6vmr= {c2h6vmr:.6g}")
    if nh3vmr  > 0: print(f"  nh3vmr = {nh3vmr:.6g}")
    if covmr   > 0: print(f"  covmr  = {covmr:.6g}")
    if h2vmr   > 0: print(f"  h2vmr  = {h2vmr:.6g}")
    if o2vmr   > 0: print(f"  o2vmr  = {o2vmr:.6g}")
    if o3vmr   > 0: print(f"  o3vmr  = {o3vmr:.6g}")
    print(f"-- H2O specific humidity --")
    print(f"  surface  : {q[-1]:.4e} kg/kg")
    print(f"  zero_h2o : {zero_h2o}")

    # Build uniform-profile arrays
    def uniform(val):
        return np.full(nlev, val)

    # Write NetCDF
    with nc.Dataset(output_file, "w") as ds:
        ds.createDimension("pverp", nilev)
        ds.createDimension("pver",  nlev)
        ds.createDimension("one",   1)

        def var1d(name, dim, data, title, units):
            v = ds.createVariable(name, "f8", (dim,))
            v.title = title
            v.units = units
            v[:] = data

        def var0d(name, data, title, units):
            v = ds.createVariable(name, "f8", ("one",))
            v.title = title
            v.units = units
            v[:] = data

        # Required: P, T, Z
        var0d("ts",   ts,   "Surface temperature",    "K")
        var0d("ps",   ps,   "Surface pressure",       "Pa")
        var1d("tmid", "pver",  tmid, "Mid-layer temperatures",  "K")
        var1d("tint", "pverp", tint, "Interface temperatures",  "K")
        var1d("pmid", "pver",  pmid, "Mid-layer pressures",     "Pa")
        var1d("pdel", "pver",  pdel, "Pressure difference",     "Pa")
        var1d("pint", "pverp", pint, "Interface pressures",     "Pa")
        var1d("zint", "pverp", zint, "Interface heights",       "m")

        # Required: albedos, cosz, mw, cp
        var0d("asdir",  asdir,  "Albedo shortwave direct",       "1")
        var0d("asdif",  asdif,  "Albedo shortwave diffuse",      "1")
        var0d("aldir",  aldir,  "Albedo longwave direct",        "1")
        var0d("aldif",  aldif,  "Albedo longwave diffuse",       "1")
        var0d("coszrs", coszrs, "Cosine of solar zenith angle",  "1")
        var0d("mw",     mwdry,  "Molecular weight of dry air",   "g/mol")
        var0d("cp",     cpdry,  "Specific heat of dry air",      "J/kg/K")

        # H2O specific humidity (always written; zero if profile has none)
        var1d("h2ommr", "pver", q, "H2O specific humidity", "kg/kg")

        # Gases: write only if VMR > 0
        if n2vmr > 0:
            var1d("n2mmr",   "pver", uniform(n2mmr),   "N2 mass mixing ratio",   "kg/kg")
        if co2vmr > 0:
            var1d("co2mmr",  "pver", uniform(co2mmr),  "CO2 mass mixing ratio",  "kg/kg")
        if ch4vmr > 0:
            var1d("ch4mmr",  "pver", uniform(ch4mmr),  "CH4 mass mixing ratio",  "kg/kg")
        if c2h6vmr > 0:
            var1d("c2h6mmr", "pver", uniform(c2h6mmr), "C2H6 mass mixing ratio", "kg/kg")
        if nh3vmr > 0:
            var1d("nh3mmr",  "pver", uniform(nh3mmr),  "NH3 mass mixing ratio",  "kg/kg")
        if covmr > 0:
            var1d("commr",   "pver", uniform(commr),   "CO mass mixing ratio",   "kg/kg")
        if h2vmr > 0:
            var1d("h2mmr",   "pver", uniform(h2mmr),   "H2 mass mixing ratio",   "kg/kg")
        if o2vmr > 0:
            var1d("o2mmr",   "pver", uniform(o2mmr),   "O2 mass mixing ratio",   "kg/kg")
        if o3vmr > 0:
            var1d("o3mmr",   "pver", uniform(o3mmr),   "O3 mass mixing ratio",   "kg/kg")

    print(f"\nWrote {output_file}")
    print("----------------------------------------------")


def print_defaults():
    print("-- default settings in makeColumn.py --")
    print(f"  output   : {OUTPUT_FILE}")
    print(f"  profile  : {PROFILE_TAG}  (available: {', '.join(AVAILABLE_PROFILES)})")
    print("  dry VMRs :")
    print(f"    co2vmr  = {CO2VMR:.6g}")
    print(f"    ch4vmr  = {CH4VMR:.6g}")
    print(f"    c2h6vmr = {C2H6VMR:.6g}")
    print(f"    nh3vmr  = {NH3VMR:.6g}")
    print(f"    covmr   = {COVMR:.6g}")
    print(f"    h2vmr   = {H2VMR:.6g}")
    print(f"    o2vmr   = {O2VMR:.6g}")
    print(f"    o3vmr   = {O3VMR:.6g}")
    print(f"    n2vmr   = (remainder)")
    print("  albedos  :")
    print(f"    asdir   = {ASDIR}")
    print(f"    asdif   = {ASDIF}")
    print(f"    aldir   = {ALDIR}")
    print(f"    aldif   = {ALDIF}")
    print(f"    coszrs  = {COSZRS}")
    print(f"  gravity  : {GRAV} m/s^2")
    print(f"  R_dry    : {R_DRY} J/kg/K")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Build RTprofile_in.nc for ExoRT.")
    parser.add_argument("--defaults", action="store_true",
                        help="Print the current default settings and exit")
    parser.add_argument("--output",  default=OUTPUT_FILE,  help="Output filename")
    parser.add_argument("--profile", default=PROFILE_TAG,
                        choices=AVAILABLE_PROFILES, help="P-T profile to use")
    parser.add_argument("--co2vmr",  type=float, default=CO2VMR)
    parser.add_argument("--ch4vmr",  type=float, default=CH4VMR)
    parser.add_argument("--c2h6vmr", type=float, default=C2H6VMR)
    parser.add_argument("--nh3vmr",  type=float, default=NH3VMR)
    parser.add_argument("--covmr",   type=float, default=COVMR)
    parser.add_argument("--h2vmr",   type=float, default=H2VMR)
    parser.add_argument("--o2vmr",   type=float, default=O2VMR)
    parser.add_argument("--o3vmr",   type=float, default=O3VMR)
    parser.add_argument("--asdir",   type=float, default=ASDIR)
    parser.add_argument("--asdif",   type=float, default=ASDIF)
    parser.add_argument("--aldir",   type=float, default=ALDIR)
    parser.add_argument("--aldif",   type=float, default=ALDIF)
    parser.add_argument("--coszrs",  type=float, default=COSZRS)
    parser.add_argument("--grav",     type=float, default=GRAV)
    parser.add_argument("--zero-h2o", action="store_true",
                        help="Force H2O specific humidity to zero at all levels")
    args = parser.parse_args()

    if args.defaults:
        print_defaults()
        sys.exit(0)

    make_column(
        output_file=args.output,
        profile_tag=args.profile,
        co2vmr=args.co2vmr,
        ch4vmr=args.ch4vmr,
        c2h6vmr=args.c2h6vmr,
        nh3vmr=args.nh3vmr,
        covmr=args.covmr,
        h2vmr=args.h2vmr,
        o2vmr=args.o2vmr,
        o3vmr=args.o3vmr,
        asdir=args.asdir,
        asdif=args.asdif,
        aldir=args.aldir,
        zero_h2o=args.zero_h2o,
        aldif=args.aldif,
        coszrs=args.coszrs,
        grav=args.grav,
    )
