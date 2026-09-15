#!/usr/bin/env python3
"""Generate CARMA haze optical-property tables on an ExoRT spectral grid.

Python port of the Mie branch of `tools/makeCARMAOptics.pro`, used to produce
`data/aerosol/haze_<grid>_b40_mie.nc`. The fractal branch of the IDL code shells
out to an external `fractaloptics.exe` (mean-field fractal aggregate solver) and
is NOT ported here; see `--help` and the module docstring of the fractal tool.

Pipeline (identical to the IDL original):
  1. Read Khare et al. (1984) Titan tholin refractive indices; k = a*10^-b.
  2. Reverse to ascending wavelength and drop points whose size parameter would
     exceed the Mie solver's limit (2*pi*100/lambda > 12000).
  3. Interpolate n and k onto a 1 cm-1 grid spanning the spectral grid, then
     average into each band -> one (n, k) pair per band.
  4. Build the CARMA mass-ratio bin structure (rmin, rmrat, rho).
  5. Mie-solve every (bin, band) pair at the band-midpoint wavelength.
  6. Kext = 3/4 * Qext / (rho * r), converted from cm2/g to m2/kg.

UNITS NOTE: the written `Kext` is in m2 kg-1. The legacy n68 files carry a
`units = "cm2 g-1"` attribute that does NOT match their own values (which are
m2 kg-1); this writer stamps the correct units.

Usage:
    python tools/makeCARMAOptics.py --grid n84 --out data/aerosol/haze_n84_b40_mie.nc
    python tools/makeCARMAOptics.py --grid n68 --out /tmp/haze_n68_check.nc
"""
import argparse
import datetime
import os
import sys

import numpy as np
import netCDF4 as nc

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

from mie import mie_single  # noqa: E402
import spectral_intervals as si  # noqa: E402

REFRACT_FILE = os.path.join(
    HERE, "..", "data", "aerosol", "refractive_indices", "Khare_haze.txt"
)

# --- CARMA element definition: Titan haze (must match the ExoCAM-CARMA setup) --
RHO_PARTICLE = 0.64      # [g cm-3] tholin particle density
NBINS = 40               # number of CARMA bins
RMIN = 1.0e-7            # [cm] radius of the smallest bin
RMRAT = 2.5              # ratio of bin masses

# Mie-solver size-parameter guard from the IDL original: refractive-index
# points with 2*pi*100/lambda >= 12000 are dropped before interpolation.
XPARAM_LIMIT = 12000.0


def read_khare(path=REFRACT_FILE):
    """Read Khare tholin indices -> (wavelength_um, n, k), ascending in lambda."""
    data = np.loadtxt(path, skiprows=4)
    wl, a, b, n = data[:, 0], data[:, 1], data[:, 2], data[:, 3]
    k = a * 10.0 ** (-b)

    # File is descending in wavelength; reverse to ascending.
    wl, n, k = wl[::-1], n[::-1], k[::-1]

    keep = (2.0 * np.pi * 100.0 / wl) < XPARAM_LIMIT
    return wl[keep], n[keep], k[keep]


def carma_bins():
    """CARMA mass-ratio bin structure -> bin_radius [cm]."""
    rmassmin = 4.0 / 3.0 * np.pi * RMIN ** 3 * RHO_PARTICLE
    ib = np.arange(NBINS)
    bin_mass = rmassmin * RMRAT ** ib
    bin_radius = (3.0 * bin_mass / (4.0 * np.pi * RHO_PARTICLE)) ** (1.0 / 3.0)
    return bin_radius


def band_average_indices(wavenum_edge, wl_ri, n_ri, k_ri):
    """Average n and k into each spectral band.

    Follows the IDL recipe: interpolate onto a uniform 1 cm-1 grid spanning the
    spectral range, then take the plain mean of the points falling inside each
    band. Beyond the range of the refractive-index data the interpolation holds
    the end values constant (np.interp default), which is what the IDL
    `interpol` + clamp produced at the grid extremes.
    """
    e = np.asarray(wavenum_edge, dtype="f8")
    nband = len(e) - 1

    # Uniform 1 cm-1 grid over the full spectral range (skip the 0 cm-1 edge).
    wn_fine = np.arange(max(e[0], 1.0), e[-1] + 1.0, 1.0)
    wl_fine = 1.0e4 / wn_fine

    # np.interp needs ascending x; refractive-index arrays are ascending in
    # wavelength, so interpolate in wavelength directly.
    n_fine = np.interp(wl_fine, wl_ri, n_ri)
    k_fine = np.interp(wl_fine, wl_ri, k_ri)
    k_fine = np.where(k_fine < 0.0, 1.0e-6, k_fine)

    n_band = np.empty(nband)
    k_band = np.empty(nband)
    for j in range(nband):
        sel = (wn_fine >= e[j]) & (wn_fine <= e[j + 1])
        if not np.any(sel):
            # Band narrower than the 1 cm-1 sampling: fall back to its midpoint.
            wl_mid = 1.0e4 / (0.5 * (e[j] + e[j + 1]))
            n_band[j] = np.interp(wl_mid, wl_ri, n_ri)
            k_band[j] = max(np.interp(wl_mid, wl_ri, k_ri), 1.0e-6)
        else:
            n_band[j] = n_fine[sel].mean()
            k_band[j] = k_fine[sel].mean()
    return n_band, k_band


def compute_optics(wavenum_edge, bin_radius, n_band, k_band):
    """Mie-solve every (bin, band) pair. Returns Qext, Kext, W, G."""
    e = np.asarray(wavenum_edge, dtype="f8")
    nband = len(e) - 1
    nbin = len(bin_radius)

    # Band-midpoint wavelength [um], as in the IDL original.
    wn_mid = 0.5 * (e[:-1] + e[1:])
    wl_mid_um = 1.0e4 / wn_mid

    qext = np.zeros((nband, nbin))
    kext = np.zeros((nband, nbin))
    walb = np.zeros((nband, nbin))
    gasy = np.zeros((nband, nbin))

    for ib in range(nbin):
        r_cm = bin_radius[ib]
        for j in range(nband):
            m = complex(n_band[j], -k_band[j])
            x = 2.0 * np.pi * r_cm / (wl_mid_um[j] / 1.0e4)
            qe, qs, g = mie_single(x, m)
            qext[j, ib] = qe
            # cm2/g -> m2/kg (factor 1/10)
            kext[j, ib] = (0.75 * qe / RHO_PARTICLE / r_cm) / 10.0
            walb[j, ib] = qs / qe if qe > 0.0 else 0.0
            gasy[j, ib] = g

    # Same clamps as the IDL original.
    walb = np.minimum(walb, 0.9999)
    gasy = np.minimum(gasy, 0.9999)
    return qext, kext, walb, gasy


def write_netcdf(path, wavenum_edge, bin_radius, qext, kext, walb, gasy,
                 grid_tag):
    nband = len(wavenum_edge) - 1
    nbin = len(bin_radius)

    with nc.Dataset(path, "w") as out:
        out.createDimension("nelements", 1)
        out.createDimension("nbins", nbin)
        out.createDimension("nwavlrng", nband)
        out.createDimension("nwave_edge", nband + 1)

        v = out.createVariable("rmrat", "f4", ("nelements",))
        v.title = "carma ratio of bin masses"
        v.units = "unitless"
        v[:] = [RMRAT]

        v = out.createVariable("rbins", "f4", ("nbins", "nelements"))
        v.title = "carma bin equivalent sphere radii"
        v.units = "cm"
        v[:] = bin_radius.reshape(nbin, 1)

        v = out.createVariable("wvnrng", "f4", ("nwave_edge",))
        v.title = "wavenumber at edges"
        v.units = "cm-1"
        v[:] = wavenum_edge

        dims = ("nwavlrng", "nbins", "nelements")
        for name, arr, title, units in (
            ("Qext", qext, "extinction efficiency", "unitless"),
            ("Kext", kext, "mass extinction efficiency", "m2 kg-1"),
            ("W", walb, "single scattering albedo", "unitless"),
            ("G", gasy, "asymmetry parameter", "unitless"),
        ):
            v = out.createVariable(name, "f4", dims)
            v.title = title
            v.units = units
            v[:] = arr.reshape(nband, nbin, 1)

        out.aerosol = "Titan tholin haze (Khare et al. 1984 optical constants)"
        out.method = "Mie theory, homogeneous spheres"
        out.spectral_grid = grid_tag
        out.rho_particle = f"{RHO_PARTICLE} g cm-3"
        out.carma_bins = (
            f"rmin={RMIN} cm, rmrat={RMRAT}, nbins={nbin}"
        )
        out.kext_units_note = (
            "Kext is in m2 kg-1. The legacy haze_n68_b40_*.nc files carry a "
            "units attribute of 'cm2 g-1' that does not match their values, "
            "which are also m2 kg-1."
        )
        out.history = (
            f"{datetime.date.today().isoformat()}: generated by "
            "tools/makeCARMAOptics.py (port of makeCARMAOptics.pro)"
        )
    print(f"wrote {path}")


def main():
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("--grid", default="n84", choices=("n28", "n42", "n68", "n84"),
                   help="ExoRT spectral grid tag (default: n84)")
    p.add_argument("--out", default=None,
                   help="output NetCDF path (default: data/aerosol/haze_<grid>_b40_mie.nc)")
    args = p.parse_args()

    # Use the raw edge arrays: get_spectral_intervals() rewrites the 0 cm-1
    # longwave edge to 1.0, but the optics files store the grid verbatim.
    wavenum_edge = np.asarray(
        getattr(si, f"wavenum_edge_{args.grid}"), dtype="f8"
    )
    out_path = args.out or os.path.join(
        HERE, "..", "data", "aerosol", f"haze_{args.grid}_b40_mie.nc"
    )

    wl_ri, n_ri, k_ri = read_khare()
    print(f"refractive indices: {len(wl_ri)} points, "
          f"{wl_ri.min():.4g}-{wl_ri.max():.4g} um")

    bin_radius = carma_bins()
    print(f"CARMA bins: {NBINS}, r = {bin_radius[0]:.3e} to "
          f"{bin_radius[-1]:.3e} cm")

    n_band, k_band = band_average_indices(wavenum_edge, wl_ri, n_ri, k_ri)
    nband = len(wavenum_edge) - 1
    print(f"spectral grid {args.grid}: {nband} bands")

    qext, kext, walb, gasy = compute_optics(
        wavenum_edge, bin_radius, n_band, k_band
    )
    write_netcdf(out_path, wavenum_edge, bin_radius, qext, kext, walb, gasy,
                 args.grid)


if __name__ == "__main__":
    main()
