#!/usr/bin/env python3
"""
plot_gas_sweep_spectra.py -- overlay n68equiv / n84equiv / exort spectra per gas.

Companion to tests/regression/gas_sweep.py run with `--save-spectra`. That writes
per-(gas, version) output files named:

    out_<TAG>_<VERSION>.nc      e.g. out_CO2_elevated_exort.nc

into tests/regression/spectra_out/. This script reads the three versions for each
gas fixture and overlays their spectral fluxes so the *spectral* origin of any
OLR difference (which gas_sweep.py only reports as a scalar) is visible:

  n68equiv : HITRAN-2016, 68-band grid   (solid)
  n84equiv : HITRAN-2016, 84-band grid   (dashed)
  exort    : HITRAN-2024, 84-band grid   (dotted, the v2 target)

Because all fluxes are plotted as spectral DENSITY (per-band flux / band width),
the 68- and 84-band grids overlay correctly on a common wavenumber/wavelength
axis even though they have different band counts. In the longwave (<~42000 cm-1)
the n68 and n84 edge arrays are identical, so any LW divergence is purely the
line list (h16 vs h24), not the grid. The line-list effect (exort vs n84) is the
diagnostic of interest -- a clean h24 table sits on top of n84.

For each gas a 2-panel figure (LW + SW) is written. A LW-difference summary panel
across all gases is also written.

Usage
-----
  python plot_gas_sweep_spectra.py                 # all gases found in spectra_out/
  python plot_gas_sweep_spectra.py --gases CO2 H2O # subset
  python plot_gas_sweep_spectra.py --indir DIR --outdir DIR
"""

import argparse
import glob
import os
import re
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import netCDF4 as nc

from spectral_intervals import get_spectral_intervals

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, os.pardir))
DEFAULT_INDIR = os.path.join(REPO, "tests", "regression", "spectra_out")
DEFAULT_OUTDIR = os.path.join(REPO, "tests", "regression", "spectra_out", "plots")

NSPEC_TO_TAG = {28: "n28", 42: "n42", 68: "n68", 84: "n84"}

# version -> (nice label, color, linestyle, linewidth)
VERSIONS = {
    "n68equiv": ("n68 (h16/68)", "tab:blue",   "-",  1.6),
    "n84equiv": ("n84 (h16/84)", "tab:green",  "--", 1.6),
    "exort":    ("exort (h24/84)", "tab:red",  ":",  2.0),
}
VERSION_ORDER = ["n68equiv", "n84equiv", "exort"]

LW_XRANGE_WN = [0, 3000]      # cm-1; thermal IR of interest
SW_XRANGE_WL = [0.2, 5.0]     # um


def read_output(fname):
    d = {}
    with nc.Dataset(fname) as f:
        for v in ["LWUP", "LWDN", "SWUP", "SWDN",
                  "LWUP_SPECTRAL", "LWDN_SPECTRAL",
                  "SWUP_SPECTRAL", "SWDN_SPECTRAL"]:
            if v in f.variables:
                d[v] = np.asarray(f.variables[v][:])
    return d


def grid_for(nspec):
    tag = NSPEC_TO_TAG.get(nspec)
    if tag is None:
        sys.exit(f"unsupported nspec={nspec}")
    nband, wn_edge, wn_mid, wl_edge, wl_mid = get_spectral_intervals(tag)
    wn_diff = np.diff(wn_edge)
    wl_diff = np.abs(np.diff(wl_edge))
    return wn_edge, wn_diff, wl_edge, wl_diff


def histxy(edges, values):
    return np.repeat(edges, 2)[1:-1], np.repeat(values, 2)


def discover_tags(indir):
    """Return ordered list of gas tags present (those with >=1 version file)."""
    tags = {}
    pat = re.compile(r"^out_(.+)_(n68equiv|n84equiv|exort)\.nc$")
    for path in sorted(glob.glob(os.path.join(indir, "out_*.nc"))):
        m = pat.match(os.path.basename(path))
        if m:
            tags.setdefault(m.group(1), set()).add(m.group(2))
    return tags


def load_versions(indir, tag):
    """Return {version: data dict} for all versions present for this tag."""
    out = {}
    for ver in VERSION_ORDER:
        path = os.path.join(indir, f"out_{tag}_{ver}.nc")
        if os.path.isfile(path):
            out[ver] = read_output(path)
    return out


def plot_gas(tag, datasets, outdir):
    """LW + SW overlay of all three versions for one gas fixture."""
    fig, (axlw, axsw) = plt.subplots(1, 2, figsize=(15, 5.2))

    # ---- LW panel: OLR (LW up at TOA), wavenumber axis ----
    for ver in VERSION_ORDER:
        if ver not in datasets:
            continue
        d = datasets[ver]
        nspec = d["LWUP_SPECTRAL"].shape[0]
        wn_edge, wn_diff, _, _ = grid_for(nspec)
        olr = d["LWUP_SPECTRAL"][:, 0] / wn_diff
        label, color, ls, lw = VERSIONS[ver]
        x, y = histxy(wn_edge, olr)
        axlw.plot(x, y, color=color, ls=ls, lw=lw,
                  label=f"{label}  (OLR={d['LWUP'][0]:.2f})")
    axlw.set_xlim(LW_XRANGE_WN)
    axlw.set_ylim(bottom=0)
    axlw.set_xlabel("Wavenumber (cm$^{-1}$)")
    axlw.set_ylabel("Spectral OLR  (W m$^{-2}$ / cm$^{-1}$)")
    axlw.set_title(f"{tag}  —  Longwave (LW↑ TOA)")
    axlw.legend(fontsize=9)
    axlw.xaxis.set_minor_locator(ticker.AutoMinorLocator())

    # ---- SW panel: SW down at surface, wavelength axis ----
    for ver in VERSION_ORDER:
        if ver not in datasets:
            continue
        d = datasets[ver]
        if "SWDN_SPECTRAL" not in d:
            continue
        nspec = d["SWDN_SPECTRAL"].shape[0]
        _, _, wl_edge, wl_diff = grid_for(nspec)
        pverp = d["SWDN"].shape[0]
        swsrf = d["SWDN_SPECTRAL"][:, pverp - 1] / wl_diff
        label, color, ls, lw = VERSIONS[ver]
        x, y = histxy(wl_edge, swsrf)
        axsw.plot(x, y, color=color, ls=ls, lw=lw,
                  label=f"{label}  (SW↓sfc={d['SWDN'][pverp-1]:.2f})")
    axsw.set_xlim(SW_XRANGE_WL)
    axsw.set_ylim(bottom=0)
    axsw.set_xlabel("Wavelength (µm)")
    axsw.set_ylabel("Spectral SW↓ surface  (W m$^{-2}$ / µm)")
    axsw.set_title(f"{tag}  —  Shortwave (SW↓ surface)")
    axsw.legend(fontsize=9)
    axsw.xaxis.set_minor_locator(ticker.AutoMinorLocator())

    fig.tight_layout()
    out = os.path.join(outdir, f"spectra_{tag}.png")
    fig.savefig(out, dpi=140)
    plt.close(fig)
    return out


def plot_linelist_diff(tags_data, outdir):
    """Per-gas LW spectral difference exort - n84 (the line-list effect), and
    grid difference n84 - n68, on a shared wavenumber axis."""
    items = [(t, d) for t, d in tags_data.items()
             if "n84equiv" in d and "exort" in d]
    if not items:
        return None
    n = len(items)
    ncol = 2
    nrow = (n + ncol - 1) // ncol
    fig, axes = plt.subplots(nrow, ncol, figsize=(14, 3.6 * nrow),
                             squeeze=False)
    for idx, (tag, d) in enumerate(items):
        ax = axes[idx // ncol][idx % ncol]
        # n84 and exort share the 84-grid -> element-wise difference is valid
        wn_edge, wn_diff, _, _ = grid_for(d["exort"]["LWUP_SPECTRAL"].shape[0])
        n84 = d["n84equiv"]["LWUP_SPECTRAL"][:, 0] / wn_diff
        exo = d["exort"]["LWUP_SPECTRAL"][:, 0] / wn_diff
        x, y = histxy(wn_edge, exo - n84)
        ax.plot(x, y, color="tab:red", lw=1.6, label="exort − n84 (line list h16→h24)")
        ax.axhline(0, color="0.5", lw=0.8)
        ax.set_xlim(LW_XRANGE_WN)
        ax.set_xlabel("Wavenumber (cm$^{-1}$)")
        ax.set_ylabel("Δ spectral OLR")
        dol = d["exort"]["LWUP"][0] - d["n84equiv"]["LWUP"][0]
        ax.set_title(f"{tag}   ΔOLR(exort−n84) = {dol:+.2f} W/m²")
        ax.legend(fontsize=8)
    # hide any unused axes
    for j in range(n, nrow * ncol):
        axes[j // ncol][j % ncol].axis("off")
    fig.suptitle("Line-list effect on spectral OLR (HITRAN-2016 → 2024)",
                 fontsize=13)
    fig.tight_layout(rect=[0, 0, 1, 0.98])
    out = os.path.join(outdir, "linelist_diff_summary.png")
    fig.savefig(out, dpi=140)
    plt.close(fig)
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--indir", default=DEFAULT_INDIR,
                    help="directory of out_<tag>_<version>.nc files")
    ap.add_argument("--outdir", default=DEFAULT_OUTDIR,
                    help="where to write PNGs")
    ap.add_argument("--gases", nargs="+", default=None,
                    help="substring filter on gas tags (e.g. CO2 H2O)")
    args = ap.parse_args()

    if not os.path.isdir(args.indir):
        sys.exit(f"ERROR: input dir not found: {args.indir}\n"
                 f"Run: python tests/regression/gas_sweep.py --save-spectra")

    tags = discover_tags(args.indir)
    if not tags:
        sys.exit(f"ERROR: no out_*.nc files in {args.indir}")

    selected = sorted(tags)
    if args.gases:
        selected = [t for t in selected
                    if any(g.lower() in t.lower() for g in args.gases)]
        if not selected:
            sys.exit(f"no tags match {args.gases}; available: {sorted(tags)}")

    os.makedirs(args.outdir, exist_ok=True)

    tags_data = {}
    print(f"Reading from {args.indir}")
    for tag in selected:
        ds = load_versions(args.indir, tag)
        if not ds:
            continue
        tags_data[tag] = ds
        present = ",".join(v for v in VERSION_ORDER if v in ds)
        out = plot_gas(tag, ds, args.outdir)
        print(f"  {tag:<16} versions=[{present}]  -> {os.path.relpath(out, REPO)}")

    summ = plot_linelist_diff(tags_data, args.outdir)
    if summ:
        print(f"  summary -> {os.path.relpath(summ, REPO)}")
    print(f"\nDone. {len(tags_data)} gas figure(s) in {args.outdir}")


if __name__ == "__main__":
    sys.exit(main())
