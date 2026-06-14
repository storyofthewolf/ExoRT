#!/usr/bin/env python
"""Plot vertical profiles from an ExoRT RTprofile netCDF file.

Python rewrite of tools/plotprofile.pro (Wolf, E.T., 2020).

Two file formats are supported and auto-detected:

  * RTprofile_out.nc (default) -- uppercase variable names; contains P-T,
    P-constituent, P-flux, and P-heating-rate fields.  A 4-panel figure is
    produced: temperature, constituent mixing ratios, up/down LW+SW fluxes,
    and LW+SW heating rates.

  * RTprofile_in.nc -- lowercase variable names; contains only P-T and
    P-constituent fields.  A 2-panel figure is produced: temperature and
    constituent mixing ratios.

Defaults:
  - vertical coordinate is pressure on a log axis (--height switches to km)
  - constituents are plotted as volume mixing ratios (--mmr switches to mass)
  - temperature is always in Kelvin

For now only a single profile is ingested and plotted at a time.
"""

import argparse
import os

import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset

# Gas constituents: legend label -> (out.nc name, in.nc name).
# H2O is specific humidity (per moist air); the rest are per dry air.
# For VMR we apply the standard MMR -> VMR conversion using the dry-air MW.
GASES = [
    ("H2O", "H2OMMR", "h2ommr"),
    ("CO2", "CO2MMR", "co2mmr"),
    ("CH4", "CH4MMR", "ch4mmr"),
    ("C2H6", "C2H6MMR", "c2h6mmr"),
    ("NH3", "NH3MMR", "nh3mmr"),
    ("CO", "COMMR", "commr"),
    ("O2", "O2MMR", "o2mmr"),
    ("O3", "O3MMR", "o3mmr"),
    ("N2", "N2MMR", "n2mmr"),
    ("H2", "H2MMR", "h2mmr"),
]

# Molecular weights (g/mol) used to convert mass -> volume mixing ratio.
MW_GAS = {
    "H2O": 18.015,
    "CO2": 44.010,
    "CH4": 16.043,
    "C2H6": 30.070,
    "NH3": 17.031,
    "CO": 28.010,
    "O2": 31.999,
    "O3": 47.998,
    "N2": 28.014,
    "H2": 2.016,
}


def _scalar(nc, *names):
    """Return the first scalar variable found among names, else None."""
    for name in names:
        if name in nc.variables:
            return float(np.asarray(nc.variables[name][:]).ravel()[0])
    return None


def read_profile(fname):
    """Read profiles from an ExoRT in.nc or out.nc file.

    Detects the format from the variable names present and returns a dict.
    The 'has_fluxes' key indicates the out.nc format (flux + heating rate
    fields available).
    """
    data = {}
    with Dataset(fname, "r") as nc:
        v = nc.variables
        # Vertical coordinates (uppercase out.nc, lowercase in.nc).
        data["PMID"] = (v["PMID"] if "PMID" in v else v["pmid"])[:]   # Pa, midlayer
        data["PINT"] = (v["PINT"] if "PINT" in v else v["pint"])[:]   # Pa, interfaces
        data["TMID"] = (v["TMID"] if "TMID" in v else v["tmid"])[:]   # K, midlayer
        data["ZINT"] = (v["ZINT"] if "ZINT" in v else v["zint"])[:]   # m, interfaces
        data["MWDRY"] = _scalar(nc, "MWDRY", "mw")

        gases = {}
        for label, up, low in GASES:
            if up in v:
                gases[label] = np.asarray(v[up][:])
            elif low in v:
                gases[label] = np.asarray(v[low][:])
        data["gases"] = gases

        # Flux / heating-rate fields exist only in the out.nc format.
        data["has_fluxes"] = "LWUP" in v
        if data["has_fluxes"]:
            for key in ("LWUP", "LWDN", "SWUP", "SWDN", "LWHR", "SWHR"):
                data[key] = np.asarray(v[key][:])
    return data


def midlayer_height(zint):
    """Interface geopotential heights (m) -> midlayer heights (km)."""
    return 0.5 * (zint[:-1] + zint[1:]) / 1000.0


def vertical_coord(data, use_height, n):
    """Return (yvals, ylabel, invert) for an array of length n (pver or pverp).

    n == len(PMID)  -> midlayer field
    n == len(PINT)  -> interface field
    """
    if use_height:
        if n == len(data["PINT"]):
            return data["ZINT"] / 1000.0, "Height (km)", False
        return midlayer_height(data["ZINT"]), "Height (km)", False
    if n == len(data["PINT"]):
        return data["PINT"] / 100.0, "Pressure (hPa)", True
    return data["PMID"] / 100.0, "Pressure (hPa)", True


def apply_vertical_axis(ax, ylabel, invert, use_height):
    ax.set_ylabel(ylabel)
    if not use_height:
        ax.set_yscale("log")
    if invert:
        # Pressure increases downward; height increases upward.
        ax.set_ylim(*sorted(ax.get_ylim(), reverse=True))


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    default_out = os.path.join(os.path.dirname(__file__), "..", "run", "RTprofile_out.nc")
    p.add_argument("infile", nargs="?", default=default_out,
                   help="ExoRT netCDF file in RTprofile_out.nc format "
                        "(default: ../run/RTprofile_out.nc)")
    p.add_argument("--infile", dest="infile_opt", default=None, metavar="FILE",
                   help="ExoRT netCDF file in RTprofile_in.nc format "
                        "(plots only P-T and P-constituents)")
    p.add_argument("--height", action="store_true",
                   help="use height (km) as the vertical coordinate instead of pressure")
    p.add_argument("--mmr", action="store_true",
                   help="plot mass mixing ratios instead of volume mixing ratios")
    p.add_argument("-o", "--output", default=None,
                   help="save the figure to this path (PNG); still shown unless --no-show")
    p.add_argument("--no-show", action="store_true",
                   help="do not open an interactive window")
    args = p.parse_args()

    fname = args.infile_opt if args.infile_opt is not None else args.infile
    data = read_profile(fname)

    # The figure layout follows the file content: 4 panels for the full output
    # file, 2 panels for the input file (no fluxes / heating rates available).
    full = data["has_fluxes"]

    # Constituent mixing ratios.
    if args.mmr:
        mr = data["gases"]
        gas_xlabel = "Mass mixing ratio (kg/kg)"
    else:
        mr = {label: mmr * (data["MWDRY"] / MW_GAS[label])
              for label, mmr in data["gases"].items()}
        gas_xlabel = "Volume mixing ratio (mol/mol)"

    ncols = 4 if full else 2
    fig, axes = plt.subplots(1, ncols, figsize=(5.5 * ncols, 6), sharey=True)
    axes = np.atleast_1d(axes)

    # --- Panel 1: temperature (midlayer) ---
    ax = axes[0]
    y, ylabel, invert = vertical_coord(data, args.height, len(data["TMID"]))
    ax.plot(data["TMID"], y, color="k", lw=2)
    ax.set_xlabel("Temperature (K)")
    ax.set_title("Temperature")
    ax.grid(True, which="both", alpha=0.3)
    apply_vertical_axis(ax, ylabel, invert, args.height)

    # --- Panel 2: constituent mixing ratios (midlayer) ---
    ax = axes[1]
    y, ylabel, invert = vertical_coord(data, args.height, len(data["TMID"]))
    cmap = plt.get_cmap("turbo")
    labels = list(mr.keys())
    for i, label in enumerate(labels):
        color = cmap(i / max(len(labels) - 1, 1))
        ax.plot(np.clip(mr[label], 1e-30, None), y, label=label, color=color, lw=2)
    ax.set_xscale("log")
    ax.set_xlim(1e-12, 2.0)
    ax.set_xlabel(gas_xlabel)
    ax.set_title("Constituent mixing ratios")
    ax.grid(True, which="both", alpha=0.3)
    ax.legend(fontsize="small", ncol=2)

    if full:
        # --- Panel 3: LW/SW up & down fluxes (interfaces) ---
        ax = axes[2]
        y, ylabel, invert = vertical_coord(data, args.height, len(data["LWUP"]))
        ax.plot(data["SWDN"], y, color="C1", ls="--", lw=2, label="SW down")
        ax.plot(data["SWUP"], y, color="C1", ls="-", lw=2, label="SW up")
        ax.plot(data["LWDN"], y, color="C3", ls="--", lw=2, label="LW down")
        ax.plot(data["LWUP"], y, color="C3", ls="-", lw=2, label="LW up")
        ax.set_xlabel(r"Flux (W m$^{-2}$)")
        ax.set_title("Radiative fluxes")
        ax.grid(True, which="both", alpha=0.3)
        ax.legend(fontsize="small")

        # --- Panel 4: LW/SW heating rates (midlayer) ---
        ax = axes[3]
        y, ylabel, invert = vertical_coord(data, args.height, len(data["LWHR"]))
        ax.plot(data["SWHR"], y, color="C1", ls="-", lw=2, label="SW")
        ax.plot(data["LWHR"], y, color="C3", ls="-", lw=2, label="LW")
        ax.axvline(0.0, color="0.5", lw=0.8)
        ax.set_xlabel("Heating rate (K/s)")
        ax.set_title("Heating/cooling rates")
        ax.grid(True, which="both", alpha=0.3)
        ax.legend(fontsize="small")

    fig.suptitle(os.path.basename(fname))
    fig.tight_layout()

    if args.output:
        fig.savefig(args.output, dpi=150, bbox_inches="tight")
        print(f"saved figure to {args.output}")
    if not args.no_show:
        plt.show()


if __name__ == "__main__":
    main()
