#!/usr/bin/env python3
"""
Stack single-column RTprofile_in.nc files into one multi-column input file
(Stage E1 'ncol' layout).

The output carries an 'ncol' dimension; every variable gains a leading
'ncol' axis in the file's C-order layout, which the Fortran reader sees as
a trailing column dimension:

  scalars   (ts, ps, coszrs, albedos, mw, cp, srf_emiss)  ->  (ncol)
  mid-layer (gases, clouds)                               ->  (ncol, pver)
  interface (tint, pint, zint)                            ->  (ncol, pverp)
  carmammr                                                ->  (ncol, nbins, nelements, pver)

Optional variables (gases, clouds, haze, srf_emiss) are unioned across the
inputs: a variable present in ANY input appears in the output, zero-filled
(1.0 for srf_emiss) for columns whose file lacks it — matching the model's
absent-variable semantics.

All inputs must share pver/pverp. Columns share the process-level runtime
config (solar_file, scon, exo_g, cloud/haze flags in user_nl_exort) — only
per-column state can vary.

Usage:
  python stackColumns.py col1.nc col2.nc [...] -o RTprofile_in.nc
"""

import argparse
import sys

import numpy as np
import netCDF4

SCALARS = ["ts", "ps", "coszrs", "asdir", "asdif", "aldir", "aldif", "mw", "cp"]
OPT_SCALARS = {"srf_emiss": 1.0}
MID = ["tmid", "pmid", "pdel"]
INT = ["tint", "pint", "zint"]
OPT_MID = ["h2ommr", "co2mmr", "ch4mmr", "c2h6mmr", "nh3mmr", "commr",
           "o2mmr", "o3mmr", "n2mmr", "h2mmr",
           "cicewp", "cliqwp", "rei", "rel", "cicewp_co2", "rei_co2"]
CARMA = "carmammr"


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[1])
    ap.add_argument("inputs", nargs="+", help="single-column RTprofile_in files")
    ap.add_argument("-o", "--output", required=True, help="multi-column output file")
    args = ap.parse_args()

    ncol = len(args.inputs)
    if ncol < 2:
        sys.exit("ERROR: need at least 2 input files to stack")

    dss = [netCDF4.Dataset(p) for p in args.inputs]
    pver = len(dss[0].dimensions["pver"])
    pverp = len(dss[0].dimensions["pverp"])
    for p, ds in zip(args.inputs, dss):
        if len(ds.dimensions["pver"]) != pver or len(ds.dimensions["pverp"]) != pverp:
            sys.exit(f"ERROR: {p} has pver/pverp mismatching {args.inputs[0]}")

    def present(name):
        return [name in ds.variables for ds in dss]

    out = netCDF4.Dataset(args.output, "w")
    out.createDimension("ncol", ncol)
    out.createDimension("pver", pver)
    out.createDimension("pverp", pverp)

    for name in SCALARS:
        v = out.createVariable(name, "f8", ("ncol",))
        v[:] = [float(np.asarray(ds.variables[name][:]).ravel()[0]) for ds in dss]

    for name, default in OPT_SCALARS.items():
        if any(present(name)):
            v = out.createVariable(name, "f8", ("ncol",))
            v[:] = [float(np.asarray(ds.variables[name][:]).ravel()[0])
                    if name in ds.variables else default for ds in dss]

    for name, dim, n in [(m, "pver", pver) for m in MID] + \
                        [(i, "pverp", pverp) for i in INT]:
        v = out.createVariable(name, "f8", ("ncol", dim))
        for ic, ds in enumerate(dss):
            v[ic, :] = np.asarray(ds.variables[name][:]).ravel()

    for name in OPT_MID:
        if not any(present(name)):
            continue
        v = out.createVariable(name, "f8", ("ncol", "pver"))
        for ic, ds in enumerate(dss):
            if name in ds.variables:
                v[ic, :] = np.asarray(ds.variables[name][:]).ravel()
            else:
                v[ic, :] = 0.0

    if any(present(CARMA)):
        # single-column layout is ("nbins","nelements","pver") — see makeColumn.py
        ref = next(ds for ds in dss if CARMA in ds.variables)
        nbins, nelem = ref.variables[CARMA].shape[0], ref.variables[CARMA].shape[1]
        out.createDimension("nbins", nbins)
        out.createDimension("nelements", nelem)
        v = out.createVariable(CARMA, "f8", ("ncol", "nbins", "nelements", "pver"))
        for ic, ds in enumerate(dss):
            if CARMA in ds.variables:
                v[ic, :, :, :] = np.asarray(ds.variables[CARMA][:])
            else:
                v[ic, :, :, :] = 0.0

    out.source = "stackColumns.py: " + " ".join(args.inputs)
    out.close()
    for ds in dss:
        ds.close()
    print(f"wrote {args.output}: ncol={ncol}, pver={pver}")


if __name__ == "__main__":
    main()
