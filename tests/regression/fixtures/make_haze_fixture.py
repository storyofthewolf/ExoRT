#!/usr/bin/env python3
"""Generate the CARMA haze regression fixtures for the Stage C3 aerosol path.

Clones the committed Earth-like TS300K column (RTprofile_in_TS300K.nc,
pver=300) and layers a haze aerosol into the upper atmosphere (layers
40..100, ~35-220 Pa) as a binwise mass mixing ratio carmammr. Mass is placed
in CARMA bins 15 (~0.10 um) and 18 (~0.24 um). Exercising this path requires
do_exo_haze=.true. in ExoRT.

Two fixtures are written, covering the optically thin and thick regimes:

  RTprofile_in_TS300K_haze.nc        visible-band tau ~ 0.5   (thin)
  RTprofile_in_TS300K_hazethick.nc   visible-band tau ~ 8.6   (thick)

The mixing ratios are SOLVED for a target optical depth against the optics
table ExoRT actually reads (data/aerosol/haze_n84_b40_mie.nc), rather than
hardcoded, so the fixtures stay meaningful if the optics are regenerated.
Re-run this script after any change to the haze optics, then rebaseline.

Cloning the existing 300-level fixture keeps the vertical grid identical to
the clear TS300K baseline, so the hazy-vs-clear comparison isolates the
aerosol effect.

Run from this directory:
    python make_haze_fixture.py
"""
import os
import shutil
import sys

import numpy as np
import netCDF4 as nc

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, "..", "..", ".."))
sys.path.insert(0, os.path.join(REPO, "tools"))

import spectral_intervals as si  # noqa: E402

SRC = os.path.join(HERE, "RTprofile_in_TS300K.nc")
OPTICS = os.path.join(REPO, "data", "aerosol", "haze_n84_b40_mie.nc")

NELEM, NBIN = 1, 40
KLO, KHI = 40, 100              # haze slab, 0-based mid-layer indices
GRAV = 9.80616                  # [m s-2], matches the TS300K case config

# Relative mass split between the two populated bins; the pair is scaled as a
# whole to hit each fixture's target optical depth.
BIN_WEIGHTS = {15: 2.0,         # ~0.10 um bin
               18: 1.0}         # ~0.24 um bin

# Target visible-band optical depths -> output filename.
TARGETS = {
    0.5: "RTprofile_in_TS300K_haze.nc",
    8.6: "RTprofile_in_TS300K_hazethick.nc",
}

VIS_WAVENUM = 20000.0           # ~0.5 um, the band the targets refer to


def visible_band_index():
    edges = np.asarray(si.wavenum_edge_n84, dtype=float)
    mids = 0.5 * (edges[:-1] + edges[1:])
    return int(np.argmin(np.abs(mids - VIS_WAVENUM))), 1.0e4 / mids[
        int(np.argmin(np.abs(mids - VIS_WAVENUM)))
    ]


def slab_masspath():
    """Mass path [kg m-2] of the haze slab per unit mass mixing ratio."""
    with nc.Dataset(SRC) as d:
        pint = np.asarray(d.variables["pint"][:], dtype=float)
    return (np.diff(pint) / GRAV)[KLO:KHI].sum()


def write_fixture(out_path, mmr, tau, jvis, wl_um):
    shutil.copyfile(SRC, out_path)
    with nc.Dataset(out_path, "a") as d:
        nlev = len(d.dimensions["pver"])
        arr = np.zeros((nlev, NELEM, NBIN))
        for ib, val in mmr.items():
            arr[KLO:KHI, 0, ib] = val
        d.createDimension("nelements", NELEM)
        d.createDimension("nbins", NBIN)
        # Fortran reads carmammr(pver,nelem,nbin): NetCDF dims fastest-last.
        v = d.createVariable("carmammr", "f8", ("nbins", "nelements", "pver"))
        v.title = "CARMA haze binwise mass mixing ratio"
        v.units = "kg/kg"
        v[:] = arr.T
        d.haze_optical_depth = (
            f"visible-band (band {jvis + 1}, {wl_um:.3f} um) tau = {tau:.4f} "
            f"against {os.path.basename(OPTICS)}"
        )
    print(f"{os.path.basename(out_path)}: layers {KLO}..{KHI}, "
          f"visible tau = {tau:.4f}")
    for ib, val in sorted(mmr.items()):
        print(f"    bin {ib}: {val:.4e} kg/kg")


def main():
    jvis, wl_um = visible_band_index()
    with nc.Dataset(OPTICS) as d:
        kext = np.asarray(d.variables["Kext"][:], dtype=float)

    masspath = slab_masspath()
    # tau per unit scale factor applied to BIN_WEIGHTS.
    tau_unit = sum(w * masspath * kext[jvis, ib, 0]
                   for ib, w in BIN_WEIGHTS.items())

    print(f"visible band {jvis + 1} ({wl_um:.3f} um), "
          f"slab mass path {masspath:.3f} kg/m2 per unit mmr")

    for tau_target, fname in sorted(TARGETS.items()):
        scale = tau_target / tau_unit
        mmr = {ib: w * scale for ib, w in BIN_WEIGHTS.items()}
        write_fixture(os.path.join(HERE, fname), mmr, tau_target, jvis, wl_um)


if __name__ == "__main__":
    main()
