#!/usr/bin/env python3
"""Make an n84 stellar spectrum from an n68 one by appending empty UV bins.

The n84 grid is the n68 grid unchanged (bands 1-68, identical edges) plus 16
UV bins above 42087 cm-1 (wavelengths below 0.238 um). An n68 file carries no
flux information there. For a cool star (an M dwarf photosphere) the true flux
below 0.238 um is negligible, so copying bands 1-68 and setting bins 69-84 to
zero reproduces the n68 stellar forcing exactly. That is what this does.

Do NOT use it for stars with real UV flux (G/K/F stars, active M dwarfs where
chromospheric UV matters, white dwarfs): regenerate those from the SED with
tools/makeStellarSpectrum_fromSED.py instead.

usage: stellar_n68_to_n84.py IN_n68.nc OUT_n84.nc [--grid data/stellar/G2V_SUN_n84.nc]
"""
import argparse
import os
import sys

import numpy as np
from netCDF4 import Dataset

HERE = os.path.dirname(os.path.abspath(__file__))
DEFAULT_GRID = os.path.join(HERE, "..", "data", "stellar", "G2V_SUN_n84.nc")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("src", help="input *_n68.nc stellar file")
    ap.add_argument("dst", help="output *_n84.nc stellar file (must not exist)")
    ap.add_argument("--grid", default=DEFAULT_GRID, help="an existing n84 stellar file supplying the band edges")
    a = ap.parse_args()
    if os.path.exists(a.dst):
        sys.exit(f"refusing to overwrite {a.dst}")

    src, grid = Dataset(a.src), Dataset(a.grid)
    n68, n84 = len(src["solarflux"]), len(grid["solarflux"])
    if (n68, n84) != (68, 84):
        sys.exit(f"expected an n68 input and an n84 grid file, got {n68} and {n84} bins")
    for v in ("wav_low", "wav_high"):
        if not np.array_equal(np.asarray(src[v][:]), np.asarray(grid[v][:68])):
            sys.exit(f"{v}: bands 1-68 of the input do not match the n84 grid")

    out = Dataset(a.dst, "w", format=grid.data_model)
    for name, dim in grid.dimensions.items():
        out.createDimension(name, len(dim))
    for name, var in grid.variables.items():
        w = out.createVariable(name, var.dtype, var.dimensions)
        for k in var.ncattrs():
            w.setncattr(k, var.getncattr(k))
    out["wav_low"][:] = grid["wav_low"][:]
    out["wav_high"][:] = grid["wav_high"][:]
    out["S0"][:] = src["S0"][:]
    flux = np.zeros(84, dtype=grid["solarflux"].dtype)
    flux[:68] = src["solarflux"][:]
    out["solarflux"][:] = flux
    out.setncattr("history", f"n84 version of {os.path.basename(a.src)} made by tools/stellar_n68_to_n84.py: "
                             "bands 1-68 copied unchanged, UV bins 69-84 (>42087 cm-1) set to 0 "
                             "(valid only for stars with negligible flux below 0.238 um).")
    out.close()
    print(f"wrote {a.dst}: sum {flux.sum():.6g} W/m2 (input {np.asarray(src['solarflux'][:], float).sum():.6g})")


if __name__ == "__main__":
    main()
