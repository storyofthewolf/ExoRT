"""
plotspectra_1D.py  --  spectral flux plots from ExoRT 1-D output.

Reads RTprofile_out.nc and plots shortwave and longwave spectra as
histogram-style band plots.  Supports optional comparison with a second file.

Usage:
  python plotspectra_1D.py                        # uses defaults below
  python plotspectra_1D.py --help
  python plotspectra_1D.py --save                 # save SW/LW to PNG instead of display
  python plotspectra_1D.py --f2 run2/RTprofile_out.nc --label2 "no NH3"
  python plotspectra_1D.py --lw-xaxis wavelength  # LW on wavelength axis (µm)
  python plotspectra_1D.py --sw-xaxis wavenumber  # SW on wavenumber axis (cm-1)
"""

import argparse
import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import netCDF4 as nc

from spectral_intervals import get_spectral_intervals

# ===========================================================================
# USER SETTINGS
# ===========================================================================

FILE1  = "../run/RTprofile_out_0nh3.nc"
LABEL1 = "run 1"

FILE2  = "../run/RTprofile_out.nc"          # set to a path string to enable comparison
LABEL2 = "run2 "

SW_XAXIS   = "wavelength"   # "wavelength" (µm) or "wavenumber" (cm-1)
LW_XAXIS   = "wavenumber"   # "wavenumber" (cm-1) or "wavelength" (µm)

SW_XRANGE_WL = [0.0, 6.0]    # SW wavelength axis range (µm)
SW_XRANGE_WN = [2000, 50000] # SW wavenumber axis range (cm-1)
SW_YRANGE    = None           # None = auto
LW_XRANGE_WN = [0, 4000]     # LW wavenumber axis range (cm-1)
LW_XRANGE_WL = [2.0, 50.0]  # LW wavelength axis range (µm)
LW_YRANGE    = None           # None = auto

SAVE_FIGS  = False          # True = save PNG; False = display interactively
SW_OUTFILE = "spectra_sw.png"
LW_OUTFILE = "spectra_lw.png"

# ===========================================================================

NSPEC_TO_TAG = {28: "n28", 42: "n42", 68: "n68", 84: "n84"}


def read_output(fname):
    """Read ExoRT output NetCDF; return dict of arrays."""
    d = {}
    with nc.Dataset(fname) as f:
        for v in ["LWUP", "LWDN", "SWUP", "SWDN",
                  "LWUP_SPECTRAL", "LWDN_SPECTRAL",
                  "SWUP_SPECTRAL", "SWDN_SPECTRAL",
                  "LWHR", "SWHR"]:
            if v in f.variables:
                d[v] = f.variables[v][:]
    return d


def make_histogram_xy(edges, values):
    """Convert band edges and per-band values to step-plot x/y pairs."""
    x = np.repeat(edges, 2)[1:-1]   # interleave edges
    y = np.repeat(values, 2)
    return x, y


def get_grid(nspec):
    tag = NSPEC_TO_TAG.get(nspec)
    if tag is None:
        sys.exit(f"ERROR: unsupported number of spectral intervals ({nspec}). "
                 f"Supported: {list(NSPEC_TO_TAG.keys())}")
    nband, wn_edge, wn_mid, wl_edge, wl_mid = get_spectral_intervals(tag)
    wn_diff = np.diff(wn_edge)                  # cm-1 per band
    # wavelength edges run high→low for bands ordered low→high in wavenumber
    # compute |Δλ| per band
    wl_diff = np.abs(np.diff(wl_edge))          # µm per band
    return wn_edge, wn_mid, wn_diff, wl_edge, wl_mid, wl_diff


def print_band_table(wn_edge, lwup_toa, lwup_srf, swdn_toa, swdn_srf):
    nband = len(lwup_toa)
    hdr = (f"{'band':>5}  {'wn1':>10}  {'wn2':>10}  "
           f"{'LWUP_TOA':>12}  {'LWUP_SRF':>12}  "
           f"{'SWDN_TOA':>12}  {'SWDN_SRF':>12}")
    print(hdr)
    print("-" * len(hdr))
    for i in range(nband):
        print(f"{i+1:5d}  {wn_edge[i]:10.2f}  {wn_edge[i+1]:10.2f}  "
              f"{lwup_toa[i]:12.4f}  {lwup_srf[i]:12.4f}  "
              f"{swdn_toa[i]:12.4f}  {swdn_srf[i]:12.4f}")


def plot_sw(data1, wn_edge, wn_diff, wl_edge, wl_diff, label1,
            data2=None, label2=None,
            xaxis="wavelength",
            xrange=None, yrange=None, save=False, outfile=SW_OUTFILE):

    pverp = data1["SWDN"].shape[0]

    if xaxis == "wavenumber":
        edges = wn_edge
        diff  = wn_diff
        xlabel = "Wavenumber (cm⁻¹)"
        ylabel = "W m⁻² cm"
    else:
        edges = wl_edge
        diff  = wl_diff
        xlabel = "Wavelength (µm)"
        ylabel = "W m⁻² µm⁻¹"

    def spectral_density(spec, level):
        return spec[:, level] / diff

    fig, ax = plt.subplots(figsize=(9, 5))

    x, y = make_histogram_xy(edges, spectral_density(data1["SWDN_SPECTRAL"], 0))
    ax.plot(x, y, color="steelblue", lw=1.5, ls="--",
            label=f"{label1} – SW↓ TOA")
    x, y = make_histogram_xy(edges, spectral_density(data1["SWDN_SPECTRAL"], pverp - 1))
    ax.plot(x, y, color="steelblue", lw=1.5, ls="-",
            label=f"{label1} – SW↓ surface")

    if data2 is not None:
        pverp2 = data2["SWDN"].shape[0]
        x, y = make_histogram_xy(edges, spectral_density(data2["SWDN_SPECTRAL"], 0))
        ax.plot(x, y, color="tomato", lw=1.5, ls="--",
                label=f"{label2} – SW↓ TOA")
        x, y = make_histogram_xy(edges, spectral_density(data2["SWDN_SPECTRAL"], pverp2 - 1))
        ax.plot(x, y, color="tomato", lw=1.5, ls="-",
                label=f"{label2} – SW↓ surface")

    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title("Shortwave spectrum")
    if xrange:
        ax.set_xlim(xrange)
    if yrange:
        ax.set_ylim(yrange)
    else:
        ax.set_ylim(bottom=0)
    ax.legend(fontsize=9)
    ax.xaxis.set_minor_locator(ticker.AutoMinorLocator())
    ax.yaxis.set_minor_locator(ticker.AutoMinorLocator())
    fig.tight_layout()

    if save:
        fig.savefig(outfile, dpi=150)
        print(f"Saved {outfile}")
    else:
        plt.show()
    plt.close(fig)


def plot_lw(data1, wn_edge, wn_diff, wl_edge, wl_diff, label1,
            data2=None, label2=None,
            xaxis="wavenumber",
            xrange=None, yrange=None, save=False, outfile=LW_OUTFILE):

    pverp = data1["LWUP"].shape[0]

    if xaxis == "wavelength":
        edges = wl_edge
        diff  = wl_diff
        xlabel = "Wavelength (µm)"
        ylabel = "W m⁻² µm⁻¹"
    else:
        edges = wn_edge
        diff  = wn_diff
        xlabel = "Wavenumber (cm⁻¹)"
        ylabel = "W m⁻² cm"

    def spectral_density(spec, level):
        return spec[:, level] / diff

    fig, ax = plt.subplots(figsize=(9, 5))

    x, y = make_histogram_xy(edges, spectral_density(data1["LWUP_SPECTRAL"], 0))
    ax.plot(x, y, color="tomato", lw=2.0, ls="-",
            label=f"{label1} – LW↑ TOA (OLR)")
    x, y = make_histogram_xy(edges, spectral_density(data1["LWUP_SPECTRAL"], pverp - 1))
    ax.plot(x, y, color="steelblue", lw=2.0, ls="-",
            label=f"{label1} – LW↑ surface")

    if data2 is not None:
        pverp2 = data2["LWUP"].shape[0]
        x, y = make_histogram_xy(edges, spectral_density(data2["LWUP_SPECTRAL"], 0))
        ax.plot(x, y, color="tomato", lw=2.0, ls="--",
                label=f"{label2} – LW↑ TOA (OLR)")
        x, y = make_histogram_xy(edges, spectral_density(data2["LWUP_SPECTRAL"], pverp2 - 1))
        ax.plot(x, y, color="steelblue", lw=2.0, ls="--",
                label=f"{label2} – LW↑ surface")

    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title("Longwave spectrum")
    if xrange:
        ax.set_xlim(xrange)
    if yrange:
        ax.set_ylim(yrange)
    else:
        ax.set_ylim(bottom=0)
    ax.legend(fontsize=9)
    ax.xaxis.set_minor_locator(ticker.AutoMinorLocator())
    ax.yaxis.set_minor_locator(ticker.AutoMinorLocator())
    fig.tight_layout()

    if save:
        fig.savefig(outfile, dpi=150)
        print(f"Saved {outfile}")
    else:
        plt.show()
    plt.close(fig)


def main(file1=FILE1, label1=LABEL1, file2=FILE2, label2=LABEL2,
         sw_xaxis=SW_XAXIS, lw_xaxis=LW_XAXIS,
         sw_xrange=None, sw_yrange=SW_YRANGE,
         lw_xrange=None, lw_yrange=LW_YRANGE,
         save=SAVE_FIGS, sw_out=SW_OUTFILE, lw_out=LW_OUTFILE):

    if sw_xrange is None:
        sw_xrange = SW_XRANGE_WN if sw_xaxis == "wavenumber" else SW_XRANGE_WL
    if lw_xrange is None:
        lw_xrange = LW_XRANGE_WL if lw_xaxis == "wavelength" else LW_XRANGE_WN

    print(f"Reading {file1}")
    d1 = read_output(file1)

    d2 = None
    if file2:
        print(f"Reading {file2}")
        d2 = read_output(file2)

    nspec = d1["LWUP_SPECTRAL"].shape[0]   # shape is (nspec, pverp)
    print(f"Spectral intervals: {nspec}")
    wn_edge, wn_mid, wn_diff, wl_edge, wl_mid, wl_diff = get_grid(nspec)

    pverp = d1["LWUP"].shape[0]

    # Band table — spectral arrays are (nspec, pverp)
    lwup_toa = d1["LWUP_SPECTRAL"][:, 0]
    lwup_srf = d1["LWUP_SPECTRAL"][:, pverp - 1]
    swdn_toa = d1["SWDN_SPECTRAL"][:, 0]
    swdn_srf = d1["SWDN_SPECTRAL"][:, pverp - 1]
    print()
    print_band_table(wn_edge, lwup_toa, lwup_srf, swdn_toa, swdn_srf)

    # Integrated totals
    print()
    print(f"{'':30s}  {'LWUP_TOA':>10}  {'LWUP_SRF':>10}  {'SWDN_TOA':>10}  {'SWDN_SRF':>10}")
    print(f"{'broadband (from flux arrays)':30s}  "
          f"{d1['LWUP'][0]:10.4f}  {d1['LWUP'][pverp-1]:10.4f}  "
          f"{d1['SWDN'][0]:10.4f}  {d1['SWDN'][pverp-1]:10.4f}")
    print(f"{'sum over spectral bands':30s}  "
          f"{lwup_toa.sum():10.4f}  {lwup_srf.sum():10.4f}  "
          f"{swdn_toa.sum():10.4f}  {swdn_srf.sum():10.4f}")

    plot_sw(d1, wn_edge, wn_diff, wl_edge, wl_diff, label1,
            data2=d2, label2=label2,
            xaxis=sw_xaxis,
            xrange=sw_xrange, yrange=sw_yrange,
            save=save, outfile=sw_out)

    plot_lw(d1, wn_edge, wn_diff, wl_edge, wl_diff, label1,
            data2=d2, label2=label2,
            xaxis=lw_xaxis,
            xrange=lw_xrange, yrange=lw_yrange,
            save=save, outfile=lw_out)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Plot ExoRT 1-D spectral output.")
    parser.add_argument("--f1",     default=FILE1,      help="Primary output file")
    parser.add_argument("--label1", default=LABEL1,     help="Label for file 1")
    parser.add_argument("--f2",     default=FILE2,      help="Comparison output file (optional)")
    parser.add_argument("--label2", default=LABEL2,     help="Label for file 2")
    parser.add_argument("--save",   action="store_true", help="Save figures to PNG instead of display")
    parser.add_argument("--sw-out", default=SW_OUTFILE, help="SW output filename")
    parser.add_argument("--lw-out", default=LW_OUTFILE, help="LW output filename")
    parser.add_argument("--sw-xaxis", default=SW_XAXIS,
                        choices=["wavelength", "wavenumber"],
                        help="SW x-axis units (default: wavelength)")
    parser.add_argument("--lw-xaxis", default=LW_XAXIS,
                        choices=["wavelength", "wavenumber"],
                        help="LW x-axis units (default: wavenumber)")
    parser.add_argument("--sw-xmax", type=float, default=None,
                        help="SW x-axis max (µm for wavelength, cm-1 for wavenumber)")
    parser.add_argument("--lw-xmax", type=float, default=None,
                        help="LW x-axis max (cm-1 for wavenumber, µm for wavelength)")
    args = parser.parse_args()

    sw_xrange = SW_XRANGE_WN if args.sw_xaxis == "wavenumber" else SW_XRANGE_WL
    if args.sw_xmax is not None:
        sw_xrange = [sw_xrange[0], args.sw_xmax]

    lw_xrange = LW_XRANGE_WL if args.lw_xaxis == "wavelength" else LW_XRANGE_WN
    if args.lw_xmax is not None:
        lw_xrange = [lw_xrange[0], args.lw_xmax]

    main(
        file1=args.f1,
        label1=args.label1,
        file2=args.f2,
        label2=args.label2,
        sw_xaxis=args.sw_xaxis,
        lw_xaxis=args.lw_xaxis,
        sw_xrange=sw_xrange,
        lw_xrange=lw_xrange,
        save=args.save,
        sw_out=args.sw_out,
        lw_out=args.lw_out,
    )
