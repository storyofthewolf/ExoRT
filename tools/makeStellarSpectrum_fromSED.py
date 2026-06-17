#!/usr/bin/env python
"""makeStellarSpectrum_fromSED.py

Author: Wolf, E.T. (original IDL, makeStellarSpectrum_fromSED.pro, Aug 2020)
Converted to Python.

DESCRIPTION:
    Creates solar/stellar spectral files for ExoRT.

    Reads a "raw" high-resolution stellar spectral energy distribution (SED)
    from a text file, bins the radiance onto one of ExoRT's supported spectral
    grids (n28, n42, n68, n84), normalises the integrated flux to a chosen
    solar constant, and writes a netCDF file usable by ExoRT.  The raw and
    binned spectra are also plotted for confirmation.

NOTES ON RAW SED FILES:
    Raw SED files vary from source to source, and this script is written to
    cope with that:

      1. different total number of rows                 -> read all rows
      2. different number of header lines                -> auto-detected
         (override with --nhead)
      3. occasionally a different number of columns      -> first two numeric
         columns (wavelength, flux) are used; --wcol/--fcol override
      4. wavelength usually in angstroms, sometimes in   -> auto-detected from
         microns (needs a scale factor)                     the value range;
                                                             override with --units
      5. ExoRT spectral intervals are no longer hardcoded here; they come from
         tools/spectral_intervals.py via get_spectral_intervals().

    The raw irradiance/flux units are arbitrary: the binned spectrum is
    renormalised to Snorm (default 1360 W/m2), so only the *shape* of the input
    matters.
"""

import argparse
import os

import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc

from spectral_intervals import get_spectral_intervals


def detect_nhead(path, max_scan=500):
    """Count leading non-data lines in a raw SED text file.

    A line is treated as data once it begins (after optional whitespace) with a
    token that parses as a float.  Comment lines (starting with #, ;, or !) and
    descriptive header text are skipped.  Returns the number of header lines.
    """
    nhead = 0
    with open(path, "r", errors="replace") as fh:
        for i, line in enumerate(fh):
            if i >= max_scan:
                break
            s = line.strip()
            if not s or s[0] in "#;!":
                nhead += 1
                continue
            try:
                float(s.split()[0])
            except (ValueError, IndexError):
                nhead += 1
                continue
            break
    return nhead


def load_sed(path, nhead, wcol, fcol):
    """Read (wavelength, flux) columns from a raw SED text file.

    Rows that cannot be parsed (ragged / stray text) are skipped.  Returns two
    1-D float arrays of equal length.  ``nhead`` lines are skipped first; the
    wavelength and flux columns are ``wcol`` and ``fcol`` (0-based).
    """
    need = max(wcol, fcol) + 1
    wav, flux = [], []
    with open(path, "r", errors="replace") as fh:
        for _ in range(nhead):
            fh.readline()
        for line in fh:
            s = line.strip()
            if not s or s[0] in "#;!":
                continue
            tok = s.replace(",", " ").split()
            if len(tok) < need:
                continue
            try:
                w = float(tok[wcol])
                f = float(tok[fcol])
            except ValueError:
                continue
            wav.append(w)
            flux.append(f)
    if not wav:
        raise ValueError(f"no numeric data rows parsed from {path} "
                         f"(check --nhead/--wcol/--fcol)")
    return np.asarray(wav, dtype=float), np.asarray(flux, dtype=float)


def detect_mu_scale(wav, units):
    """Return the factor that converts the wavelength column to microns.

    With units='auto', angstroms vs microns is inferred from the maximum
    wavelength: optical/IR SEDs span up to ~1e5-1e7 A but only ~1e1-1e3 um, so
    a large maximum implies angstroms.
    """
    if units == "micron":
        return 1.0
    if units == "angstrom":
        return 1.0e-4
    # auto: angstrom values run large; micron values stay small.
    if np.nanmax(wav) > 5.0e3:
        print("  wavelength units auto-detected: angstroms (scaling by 1e-4)")
        return 1.0e-4
    print("  wavelength units auto-detected: microns (no scaling)")
    return 1.0


def bin_spectrum(wnm, sunm, rtwavlow, rtwavhi):
    """Bin raw radiance onto the ExoRT grid.

    For each interval [wn_lw, wn_hi] (cm^-1), average the raw radiance of all
    raw points whose wavenumber falls in the interval, then multiply by the
    interval's wavelength width (microns).  Mirrors the IDL loop exactly.
    Returns (solarflux, counts) with one entry per interval.
    """
    nrtwavl = len(rtwavlow)
    solarflux = np.zeros(nrtwavl)
    counts = np.zeros(nrtwavl, dtype=int)
    nxsum = 0
    for iw in range(nrtwavl):
        wn_lw = rtwavlow[iw]
        wn_hi = rtwavhi[iw]
        dwl = 1.0e4 / wn_lw - 1.0e4 / wn_hi
        mask = (wnm >= wn_lw) & (wnm <= wn_hi)
        nx = int(np.count_nonzero(mask))
        counts[iw] = nx
        nxsum += nx
        ftmp = np.mean(sunm[mask]) * dwl if nx > 0 else 0.0
        if ftmp >= 0.0:
            solarflux[iw] = ftmp
        print(f"  {iw + 1:3d}  {wn_lw:12.3f}  {wn_hi:12.3f}  {nx:8d}  {ftmp:.6e}")
    print(f"  nxsum {nxsum}")
    return solarflux, counts


def write_netcdf(outname, rtwavlow, rtwavhi, S0, solarflux):
    """Write the binned, normalised SED to an ExoRT-format netCDF file."""
    with nc.Dataset(outname, "w", format="NETCDF4") as ds:
        ds.createDimension("nw", len(rtwavlow))
        ds.createDimension("one", 1)

        v1 = ds.createVariable("wav_low", "f4", ("nw",))
        v2 = ds.createVariable("wav_high", "f4", ("nw",))
        v3 = ds.createVariable("S0", "f4", ("one",))
        v4 = ds.createVariable("solarflux", "f4", ("nw",))

        v1.title = "wavenumber edge, low"; v1.units = "cm -1"
        v2.title = "wavenumber edge, high"; v2.units = "cm -1"
        v3.title = "solar constant"; v3.units = "W/m2"
        v4.title = "spectrally integrated solar flux"
        v4.units = "W/m2 in each interval"

        v1[:] = rtwavlow
        v2[:] = rtwavhi
        v3[:] = S0
        v4[:] = solarflux


def plot_spectra(wlgth, sunm, scalefac, rtwavlow, rtwavhi, solarflux,
                 rtwvldel, plot_out, xmax=6.0):
    """Plot the scaled raw spectrum against the binned step spectrum."""
    # Step ("bar") representation of the binned SED, in radiance per micron.
    nrtwavl = len(rtwavlow)
    xbar = np.empty(nrtwavl * 2)
    ybar = np.empty(nrtwavl * 2)
    for i in range(nrtwavl):
        xbar[2 * i] = 1.0e4 / rtwavlow[i]
        xbar[2 * i + 1] = 1.0e4 / rtwavhi[i]
        ybar[2 * i] = solarflux[i] / rtwvldel[i]
        ybar[2 * i + 1] = solarflux[i] / rtwvldel[i]

    fig, ax = plt.subplots(figsize=(8, 6))
    ax.plot(wlgth, sunm * scalefac, color="steelblue", lw=1.0,
            label="raw spectrum (scaled)")
    ax.plot(xbar, ybar, color="darkorange", lw=2.5, label="binned SED")
    ax.set_xlim(0.0, xmax)
    ax.set_xlabel("wavelength (microns)")
    ax.set_ylabel(r"Radiance (W m$^{-2}$ micron$^{-1}$)")
    ax.legend()
    fig.tight_layout()
    fig.savefig(plot_out, dpi=150)
    plt.close(fig)
    print(f"plotting, {plot_out}")


def main():
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("infile", help="raw stellar SED text file")
    p.add_argument("-o", "--outname", default=None,
                   help="output netCDF name "
                        "(default: <infile stem>_<tag>.nc)")
    p.add_argument("-r", "--resolution", default="n68",
                   choices=["n28", "n42", "n68", "n84"],
                   help="ExoRT spectral resolution tag (default: n68)")
    p.add_argument("--snorm", type=float, default=1360.0,
                   help="normalise integrated flux to this value [W/m2] "
                        "(default: 1360)")
    p.add_argument("--units", default="auto",
                   choices=["auto", "angstrom", "micron"],
                   help="wavelength units of the raw file (default: auto-detect)")
    p.add_argument("--nhead", type=int, default=None,
                   help="number of header lines to skip "
                        "(default: auto-detect)")
    p.add_argument("--wcol", type=int, default=0,
                   help="0-based wavelength column index (default: 0)")
    p.add_argument("--fcol", type=int, default=1,
                   help="0-based flux column index (default: 1)")
    p.add_argument("--fscale", type=float, default=1.0,
                   help="multiplicative scale applied to the raw flux column "
                        "(default: 1.0)")
    p.add_argument("--no-write", action="store_true",
                   help="do not write the netCDF file")
    p.add_argument("--no-plot", action="store_true",
                   help="do not write the confirmation plot")
    args = p.parse_args()

    tag = args.resolution
    if args.outname is None:
        stem = os.path.splitext(os.path.basename(args.infile))[0]
        outname = f"{stem}_{tag}.nc"
    else:
        outname = args.outname

    print(f"reading {args.infile}")

    # -- ingest the raw SED (robust to header / row / column variation) --
    nhead = args.nhead if args.nhead is not None else detect_nhead(args.infile)
    print(f"  skipping {nhead} header line(s)")
    wav_raw, flux_raw = load_sed(args.infile, nhead, args.wcol, args.fcol)
    print(f"  read {wav_raw.size} data rows")

    mu_scale = detect_mu_scale(wav_raw, args.units)

    # Convert to wavelength (microns), wavenumber (cm^-1), and scaled radiance.
    wlgth = wav_raw * mu_scale
    wnm = 1.0e4 / wlgth
    sunm = flux_raw * args.fscale

    # Sort by wavenumber so binning masks are well behaved regardless of the
    # input ordering (raw files may run blue->red or red->blue).
    order = np.argsort(wnm)
    wnm = wnm[order]
    sunm = sunm[order]
    wlgth = wlgth[order]

    # -- ExoRT spectral grid (edges -> low/high interval bounds) --
    (nrtwavl, wavenum_edge, _wavenum_mid,
     _wavelength_edge, _wavelength_mid) = get_spectral_intervals(tag)
    # spectral_intervals.py already substitutes the longwave 0 cm^-1 edge with
    # 1.0 cm^-1, so the conversion below stays finite.
    rtwavlow = wavenum_edge[:-1]
    rtwavhi = wavenum_edge[1:]
    rtwvldel = 1.0e4 / rtwavlow - 1.0e4 / rtwavhi
    print(f"using {nrtwavl} spectral intervals ({tag})")

    # -- bin and normalise --
    solarflux, _counts = bin_spectrum(wnm, sunm, rtwavlow, rtwavhi)

    total = np.sum(solarflux)
    print(f"TOTAL SOLAR FLUX (pre-norm): {total}")
    scalefac = args.snorm / total
    solarflux_out = solarflux * scalefac
    S0 = np.sum(solarflux_out)
    print(f"TOTAL SOLAR FLUX (post-norm): {S0}")
    print(f"ScaleFac {scalefac}")

    # -- write netCDF --
    if not args.no_write:
        print(f"writing to {outname}")
        write_netcdf(outname, rtwavlow, rtwavhi, S0, solarflux_out)

    # -- confirmation plot --
    if not args.no_plot:
        plot_out = outname + ".plot.png"
        plot_spectra(wlgth, sunm, scalefac, rtwavlow, rtwavhi,
                     solarflux_out, rtwvldel, plot_out)


if __name__ == "__main__":
    main()
