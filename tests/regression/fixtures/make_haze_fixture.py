#!/usr/bin/env python3
"""Generate the CARMA haze regression fixture for the Stage C3 aerosol path.

Clones the committed Earth-like TS300K column (RTprofile_in_TS300K.nc,
pver=300) and layers a fractal-haze aerosol into the upper atmosphere
(layers 40..100, ~35-220 Pa) as a binwise mass mixing ratio carmammr.
Mass is placed in CARMA bins 15 (~0.10 um) and 18 (~0.24 um), giving a
visible-band haze optical depth of ~0.5 against the provisional
haze_n84_b40_fractal_interp.nc optics. Exercising this path requires
do_exo_haze=.true. in ExoRT.

Cloning the existing 300-level fixture keeps the vertical grid identical to
the clear TS300K baseline, so the hazy-vs-clear comparison isolates the
aerosol effect.

Run from this directory:
    python make_haze_fixture.py
writes RTprofile_in_TS300K_haze.nc here.
"""
import os
import shutil

import numpy as np
import netCDF4 as nc

HERE = os.path.dirname(os.path.abspath(__file__))
SRC = os.path.join(HERE, "RTprofile_in_TS300K.nc")
OUT = os.path.join(HERE, "RTprofile_in_TS300K_haze.nc")

NELEM, NBIN = 1, 40
KLO, KHI = 40, 100              # haze slab, 0-based mid-layer indices
BIN_MMR = {15: 2.0e-5,          # ~0.10 um bin [kg/kg]
           18: 1.0e-5}          # ~0.24 um bin [kg/kg]


def main():
    shutil.copyfile(SRC, OUT)
    with nc.Dataset(OUT, "a") as d:
        nlev = len(d.dimensions["pver"])
        arr = np.zeros((nlev, NELEM, NBIN))
        for ib, val in BIN_MMR.items():
            arr[KLO:KHI, 0, ib] = val
        d.createDimension("nelements", NELEM)
        d.createDimension("nbins", NBIN)
        # Fortran reads carmammr(pver,nelem,nbin): NetCDF dims fastest-last.
        v = d.createVariable("carmammr", "f8", ("nbins", "nelements", "pver"))
        v.title = "CARMA haze binwise mass mixing ratio"
        v.units = "kg/kg"
        v[:] = arr.T
    print(f"haze fixture: {OUT}")
    print(f"  haze layers {KLO}..{KHI}: bins "
          + ", ".join(f"{ib} @ {val:g} kg/kg" for ib, val in BIN_MMR.items()))


if __name__ == "__main__":
    main()
