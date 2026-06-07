#!/usr/bin/env python3
"""
regrid_rtprofile.py
-------------------
Regrid the vertical structure of an ExoRT RTprofile_in.nc file to an
arbitrary number of layers (pver), while preserving temperature, pressure,
and constituent mixing-ratio profiles.

Regridding strategy
-------------------
All interpolation is performed in log(pressure) space, which is the natural
vertical coordinate for atmospheric profiles (profiles are smooth in ln(p)
because the atmosphere is nearly in hydrostatic balance with an approximately
constant scale height).  This is analogous to the way ExoRT itself uses
ln(p) as the vertical coordinate internally.

The interface pressure grid is regenerated to span the same TOA and surface
pressures as the original, with the new pver+1 interfaces evenly spaced in
log(p).  All mid-layer quantities (tmid, pmid, pdel, constituent MMRs) are
then derived from the new interface grid or interpolated onto it.

Invariants enforced (same as ExoRT/CAM convention):
  pdel[k] = pint[k+1] - pint[k]          (layer thickness)
  pmid[k] = 0.5*(pint[k] + pint[k+1])    (mid-layer pressure)
  tmid[k] = 0.5*(tint[k] + tint[k+1])    (mid-layer temperature)
  zint[-1] = 0                            (surface height = 0)
  zint[0]  = zint[1]                      (TOA convention)

Usage
-----
  python regrid_rtprofile.py RTprofile_in_TS320K.nc 60 --desc myfoo
  python regrid_rtprofile.py RTprofile_in_TS320K.nc 120
"""

import argparse
import os
import sys
import warnings
import numpy as np
import scipy.io as sio
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d


# ---------------------------------------------------------------------------
# Physical constants
# ---------------------------------------------------------------------------
GRAV = 9.80665       # m s-2, standard gravity
R_UNIV = 8.31446     # J mol-1 K-1, universal gas constant


# ---------------------------------------------------------------------------
# I/O helpers
# ---------------------------------------------------------------------------

def load_rtprofile(path):
    """
    Read an RTprofile_in.nc file and return its contents as a plain dict.

    Returns
    -------
    dict with keys:
      'pver'  : int, number of mid-layers in input file
      'scalars': dict of scalar variables  {name: (value, attrs)}
      'mid'   : dict of mid-layer arrays   {name: (array, attrs)}
      'int'   : dict of interface arrays   {name: (array, attrs)}
    """
    ds = sio.netcdf_file(path, "r", mmap=False)

    pver = ds.dimensions["pver"]

    # Collect variable attributes as plain Python dicts (strings decoded)
    def get_attrs(var):
        attrs = {}
        for k in var._attributes:
            v = getattr(var, k)
            attrs[k] = v.decode() if isinstance(v, bytes) else v
        return attrs

    # Identify variables by their dimension
    scalars, mid_vars, int_vars = {}, {}, {}

    for vname, var in ds.variables.items():
        data = var.data.copy()
        attrs = get_attrs(var)
        dims = var.dimensions

        if "pver" in dims:
            mid_vars[vname] = (data, attrs)
        elif "pverp" in dims:
            int_vars[vname] = (data, attrs)
        else:
            # 'one' dimension — scalar
            scalars[vname] = (data, attrs)

    ds.close()
    return {"pver": pver, "scalars": scalars, "mid": mid_vars, "int": int_vars}


def write_rtprofile(data, path):
    """
    Write a regridded RTprofile dict to a NetCDF file at *path*.
    Dimension names and variable attributes are preserved exactly.
    """
    pver  = data["pver"]
    pverp = pver + 1

    f = sio.netcdf_file(path, "w")
    f.createDimension("pver",  pver)
    f.createDimension("pverp", pverp)
    f.createDimension("one",   1)

    def _write_var(name, arr, attrs, dims):
        # Choose float32 to match typical ExoRT output
        v = f.createVariable(name, "f", dims)
        for k, val in attrs.items():
            setattr(v, k, val)
        v[:] = arr.astype(np.float32)

    for name, (arr, attrs) in data["scalars"].items():
        _write_var(name, arr, attrs, ("one",))

    for name, (arr, attrs) in data["mid"].items():
        _write_var(name, arr, attrs, ("pver",))

    for name, (arr, attrs) in data["int"].items():
        _write_var(name, arr, attrs, ("pverp",))

    f.close()


# ---------------------------------------------------------------------------
# Regridding core
# ---------------------------------------------------------------------------

def make_new_pint(pint_old, pver_new):
    """
    Generate a new interface pressure grid with pver_new+1 levels, spanning
    the same TOA and surface pressures as the original, evenly spaced in
    log(pressure).

    Parameters
    ----------
    pint_old : 1-D array, shape (pver_old+1,)
        Original interface pressures in Pa, ordered TOA -> surface (ascending).
    pver_new : int
        Target number of mid-layers.

    Returns
    -------
    pint_new : 1-D array, shape (pver_new+1,)
    """
    p_toa = pint_old[0]   # top of atmosphere (smallest pressure)
    p_sfc = pint_old[-1]  # surface (largest pressure)

    # Evenly spaced in log(p) — equivalent to geometric progression in linear p.
    # This mirrors the level distribution strategy common in GCM sigma coordinates.
    pint_new = np.exp(np.linspace(np.log(p_toa), np.log(p_sfc), pver_new + 1))
    return pint_new


def regrid_rtprofile(data_in, pver_new):
    """
    Regrid a loaded RTprofile dict to *pver_new* layers.

    Interface temperatures are interpolated from the old pint grid to the new
    pint grid in log(p) space.  Mid-layer quantities are derived from the new
    interface grid or interpolated from old pmid to new pmid.

    Parameters
    ----------
    data_in : dict, as returned by load_rtprofile()
    pver_new : int

    Returns
    -------
    data_out : dict in the same format, with updated arrays and pver.
    """
    pint_old, _ = data_in["int"]["pint"]
    tint_old, _ = data_in["int"]["tint"]

    pmid_old, _ = data_in["mid"]["pmid"]

    # --- 1. New interface pressure grid (log-spaced) -----------------------
    pint_new = make_new_pint(pint_old, pver_new)

    # --- 2. Interpolate tint onto new pint grid ----------------------------
    # Interpolate in log(p).  The temperature profile is smooth in log(p).
    # Extrapolation is set to the boundary values (should not occur since we
    # span the exact same pressure range).
    tint_interp = interp1d(
        np.log(pint_old), tint_old,
        kind="linear",
        bounds_error=False,
        fill_value=(tint_old[0], tint_old[-1]),
    )
    tint_new = tint_interp(np.log(pint_new))

    # --- 3. Derive mid-layer quantities from interfaces --------------------
    pmid_new = 0.5 * (pint_new[:-1] + pint_new[1:])
    pdel_new = pint_new[1:] - pint_new[:-1]
    tmid_new = 0.5 * (tint_new[:-1] + tint_new[1:])

    # --- 4. Interpolate constituent MMRs -----------------------------------
    # MMRs live on mid-layer centers.  Interpolate from old pmid to new pmid
    # in log(p) space.  Clamp to non-negative to avoid interpolation artifacts.
    mmr_names = [n for n in data_in["mid"] if n not in ("pmid", "pdel", "tmid")]

    mmr_new = {}
    for name in mmr_names:
        arr_old, attrs = data_in["mid"][name]
        mmr_interp = interp1d(
            np.log(pmid_old), arr_old,
            kind="linear",
            bounds_error=False,
            fill_value=(arr_old[0], arr_old[-1]),
        )
        regridded = mmr_interp(np.log(pmid_new))
        # Physical constraint: mixing ratios must be >= 0
        regridded = np.maximum(regridded, 0.0)
        mmr_new[name] = (regridded, attrs)

    # --- 5. Reconstruct zint via hydrostatic equation ----------------------
    # Hydrostatic: dz = (R*T_mid / g) * ln(p_int[k+1] / p_int[k])
    # z decreases from TOA (k=0) toward surface (k=pverp-1), zint[-1] = 0.
    # Integrate upward from the surface.
    mw_amu = data_in["scalars"]["mw"][0][0]   # scalar, AMU
    R_spec = R_UNIV / (mw_amu * 1e-3)         # J kg-1 K-1

    dz = R_spec * tmid_new / GRAV * np.log(pint_new[1:] / pint_new[:-1])

    zint_new = np.zeros(pver_new + 1)
    zint_new[-1] = 0.0
    for k in range(pver_new - 1, -1, -1):
        zint_new[k] = zint_new[k + 1] + dz[k]
    # ExoRT convention: zint[0] == zint[1] (top boundary duplicated)
    zint_new[0] = zint_new[1]

    # --- 6. Assemble output dict -------------------------------------------
    data_out = {
        "pver": pver_new,
        "scalars": {k: v for k, v in data_in["scalars"].items()},  # unchanged
        "mid": {},
        "int": {},
    }

    # Interface variables
    _, pint_attrs = data_in["int"]["pint"]
    _, tint_attrs = data_in["int"]["tint"]
    _, zint_attrs = data_in["int"]["zint"]
    data_out["int"]["pint"] = (pint_new, pint_attrs)
    data_out["int"]["tint"] = (tint_new, tint_attrs)
    data_out["int"]["zint"] = (zint_new, zint_attrs)

    # Mid-layer pressure/thickness/temperature
    _, pmid_attrs = data_in["mid"]["pmid"]
    _, pdel_attrs = data_in["mid"]["pdel"]
    _, tmid_attrs = data_in["mid"]["tmid"]
    data_out["mid"]["tmid"] = (tmid_new, tmid_attrs)
    data_out["mid"]["pmid"] = (pmid_new, pmid_attrs)
    data_out["mid"]["pdel"] = (pdel_new, pdel_attrs)

    # Constituent MMRs
    data_out["mid"].update(mmr_new)

    return data_out


# ---------------------------------------------------------------------------
# Plotting
# ---------------------------------------------------------------------------

def plot_comparison(data_old, data_new, out_path):
    """
    Two-panel comparison plot of T(p) and H2O MMR(p) profiles
    for the original and regridded files.

    Pressure axis is logarithmic and inverted (surface at bottom).
    """
    pmid_old, _ = data_old["mid"]["pmid"]
    tmid_old, _ = data_old["mid"]["tmid"]
    h2o_old,  _ = data_old["mid"]["h2ommr"]

    pmid_new, _ = data_new["mid"]["pmid"]
    tmid_new, _ = data_new["mid"]["tmid"]
    h2o_new,  _ = data_new["mid"]["h2ommr"]

    pver_old = data_old["pver"]
    pver_new = data_new["pver"]

    fig, axes = plt.subplots(1, 2, figsize=(10, 6), sharey=True)

    # Pressure axis setup (shared y)
    p_min = min(pmid_old[0], pmid_new[0]) * 0.9
    p_max = max(pmid_old[-1], pmid_new[-1]) * 1.05

    for ax in axes:
        ax.set_yscale("log")
        ax.invert_yaxis()
        ax.set_ylim(p_max, p_min)
        ax.set_ylabel("Pressure (Pa)")
        ax.grid(True, which="both", alpha=0.3, linestyle="--")

    # Panel 1: Temperature
    ax = axes[0]
    ax.plot(tmid_old, pmid_old, color="steelblue",  lw=1.5, label=f"Original (pver={pver_old})")
    ax.plot(tmid_new, pmid_new, color="tomato", lw=1.0, ls="--", label=f"Regridded (pver={pver_new})")
    ax.set_xlabel("Temperature (K)")
    ax.set_title("Temperature Profile")
    ax.legend(fontsize=9)

    # Panel 2: H2O MMR
    ax = axes[1]
    ax.plot(h2o_old, pmid_old, color="steelblue",  lw=1.5, label=f"Original (pver={pver_old})")
    ax.plot(h2o_new, pmid_new, color="tomato", lw=1.0, ls="--", label=f"Regridded (pver={pver_new})")
    ax.set_xscale("log")
    ax.set_xlabel("H$_2$O MMR (kg/kg)")
    ax.set_title("Water Vapor Profile")
    ax.legend(fontsize=9)

    fig.suptitle("RTprofile Vertical Regrid Comparison", fontsize=12, fontweight="bold")
    plt.tight_layout()
    plt.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"  Plot saved: {out_path}")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def parse_args():
    p = argparse.ArgumentParser(
        description="Regrid an RTprofile_in.nc file to a new vertical resolution."
    )
    p.add_argument("input",  help="Path to input RTprofile_in.nc file")
    p.add_argument("pver_new", type=int, help="Target number of mid-layers (pver)")
    p.add_argument(
        "--desc", default=None,
        help="Descriptor string for output filename (default: derived from input name)"
    )
    p.add_argument(
        "--outdir", default=".",
        help="Output directory (default: current directory)"
    )
    return p.parse_args()


def main():
    args = parse_args()

    if not os.path.isfile(args.input):
        sys.exit(f"ERROR: Input file not found: {args.input}")

    if args.pver_new < 2:
        sys.exit("ERROR: pver_new must be >= 2")

    # Derive descriptor from input filename if not provided
    # e.g. 'RTprofile_in_TS320K.nc' -> desc='TS320K'
    if args.desc is None:
        base = os.path.splitext(os.path.basename(args.input))[0]
        # Strip leading 'RTprofile_in_' if present
        prefix = "RTprofile_in_"
        if base.startswith(prefix):
            args.desc = base[len(prefix):]
        else:
            args.desc = base

    pver_new = args.pver_new
    desc     = args.desc

    print(f"Input file : {args.input}")

    # Load
    print("Loading input profile...")
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        data_old = load_rtprofile(args.input)

    pver_old = data_old["pver"]
    pint_old = data_old["int"]["pint"][0]
    print(f"  pver = {pver_old}  (pverp = {pver_old+1})")
    print(f"  pint range: {pint_old[0]:.4f} – {pint_old[-1]:.2f} Pa")

    if pver_new == pver_old:
        print(f"WARNING: pver_new ({pver_new}) == pver_old ({pver_old}). "
              "Output will be a resampled copy on the same log-p grid.")

    # Regrid
    print(f"\nRegridding to pver = {pver_new} ...")
    data_new = regrid_rtprofile(data_old, pver_new)

    pint_new = data_new["int"]["pint"][0]
    print(f"  pint range: {pint_new[0]:.4f} – {pint_new[-1]:.2f} Pa  (check)")
    print(f"  pdel sum = {data_new['mid']['pdel'][0].sum():.2f} Pa  "
          f"(ps = {data_old['scalars']['ps'][0][0]:.2f} Pa)")

    # Write output NetCDF
    os.makedirs(args.outdir, exist_ok=True)
    nc_name = f"RTprofile_in_{desc}_pver{pver_new}.nc"
    nc_path = os.path.join(args.outdir, nc_name)
    print(f"\nWriting: {nc_path}")
    write_rtprofile(data_new, nc_path)

    # Diagnostic: verify round-trip invariants
    _verify(data_new)

    # Plot
    plot_name = f"RTprofile_regrid_{desc}_pver{pver_new}.png"
    plot_path = os.path.join(args.outdir, plot_name)
    print(f"Plotting comparison...")
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        plot_comparison(data_old, data_new, plot_path)

    print("\nDone.")


def _verify(data):
    """Spot-check that internal invariants hold on the regridded data."""
    pint  = data["int"]["pint"][0]
    tint  = data["int"]["tint"][0]
    pmid  = data["mid"]["pmid"][0]
    pdel  = data["mid"]["pdel"][0]
    tmid  = data["mid"]["tmid"][0]

    ok_pdel = np.allclose(pdel, pint[1:] - pint[:-1])
    ok_pmid = np.allclose(pmid, 0.5*(pint[:-1] + pint[1:]))
    ok_tmid = np.allclose(tmid, 0.5*(tint[:-1] + tint[1:]))

    print(f"  Invariant checks: pdel={ok_pdel}, pmid={ok_pmid}, tmid={ok_tmid}")
    if not (ok_pdel and ok_pmid and ok_tmid):
        print("  WARNING: one or more invariants violated — check regrid logic.")


if __name__ == "__main__":
    main()
