#!/usr/bin/env python3
"""
plot_minor_gas_cases.py -- figures for the CH4 / C2H6 / NH3 / CO regression cases.

The regression suite's minor-gas cases (tests/regression, `_gas_cases()`) put
one gas at a time in an N2 background on the TS273K profile, at a realistic and
an elevated abundance, under a G2V and a 3400 K blackbody star. This script
plots each gas against a pure-N2 reference column (the same deck with the gas
zeroed), so the bands each gas acts in are visible:

  * OLR spectral density vs wavenumber (star-independent, so G2V only)
  * surface SW-down spectral density vs wavelength, one panel per star
  * LW heating rate, and SW heating rate per star, vs pressure

Inputs: the committed baselines tests/regression/baselines/RTprofile_out_<GAS>_
<level>_<star>.nc. The N2 reference is not a regression case; it is computed by
running run/exort.exe (build it first: cd build && make exort) and cached in
tests/regression/outputs/ (gitignored).

Outputs (default figures_minor_gases/):
  minor_gas_<GAS>.png      one 2x3 figure per gas
  minor_gas_summary.png    OLR and surface-SW change vs N2, all gases
  minor_gas_data.json      the plotted series
  minor_gas_page.html      self-contained interactive page (open in a browser);
                           tools/build_figure_pages.py builds it and its siblings

Heating rates in RTprofile_out.nc are K/day (io_1D converts from K/s). Baselines
written before 2026-09-23 carry a stale `units = "K/s"` attribute; the values
were always K/day.

Usage
-----
  python tools/plot_minor_gas_cases.py
  python tools/plot_minor_gas_cases.py --gases CH4 NH3 --outdir /tmp/figs
  python tools/plot_minor_gas_cases.py --rerun-reference   # recompute the N2 column
"""

import argparse
import json
import os
import shutil
import subprocess
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import netCDF4 as nc

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, os.pardir))
sys.path.insert(0, HERE)
from spectral_intervals import get_spectral_intervals  # noqa: E402

REG = os.path.join(REPO, "tests", "regression")
BASELINES = os.path.join(REG, "baselines")
FIXTURES = os.path.join(REG, "fixtures")
CACHE = os.path.join(REG, "outputs")
RUN_DIR = os.path.join(REPO, "run")

GASES = ["CH4", "C2H6", "NH3", "CO"]
GAS_LABEL = {"CH4": "CH$_4$", "C2H6": "C$_2$H$_6$", "NH3": "NH$_3$", "CO": "CO"}
GAS_VAR = {"CH4": "ch4mmr", "C2H6": "c2h6mmr", "NH3": "nh3mmr", "CO": "commr"}
# Volume mixing ratios of the committed fixtures (see run_regression.py).
ABUND = {"CH4": ("1.8 ppmv", "1%"), "C2H6": ("1 ppbv", "100 ppmv"),
         "NH3": ("1 ppbv", "100 ppmv"), "CO": ("100 ppbv", "1%")}
STARS = {"G2V": ("G2V_SUN_n84.nc", "G2V (Sun)"),
         "BB3400": ("blackbody_3400K_n84.nc", "3400 K blackbody")}

# Reference palette (dataviz skill), light-mode steps: slot 1 blue, slot 2
# orange; the N2 reference is neutral secondary ink, dashed.
C_REF, C_REAL, C_ELEV = "#52514e", "#2a78d6", "#eb6834"
INK, MUTED, GRID = "#0b0b0b", "#898781", "#e1e0d9"


def nf_env():
    env = os.environ.copy()
    try:
        prefix = subprocess.run(["nf-config", "--prefix"], capture_output=True,
                                text=True).stdout.strip()
    except FileNotFoundError:
        prefix = ""
    if prefix:
        for key in ("DYLD_LIBRARY_PATH", "LD_LIBRARY_PATH"):
            env[key] = os.path.join(prefix, "lib") + os.pathsep + env.get(key, "")
    return env


def reference_path(star):
    return os.path.join(CACHE, f"RTprofile_out_N2only_{star}.nc")


def make_reference(star, rerun=False):
    """Run exort.exe on the TS273K N2-only column (a gas fixture with the gas zeroed)."""
    out = reference_path(star)
    if os.path.exists(out) and not rerun:
        return out
    exe = os.path.join(RUN_DIR, "exort.exe")
    if not os.path.exists(exe):
        sys.exit("run/exort.exe not found: build it first (cd build && make exort)")
    os.makedirs(CACHE, exist_ok=True)
    deck = os.path.join(RUN_DIR, "RTprofile_in.nc")
    namelist = os.path.join(RUN_DIR, "user_nl_exort")
    saved = open(namelist).read() if os.path.exists(namelist) else None
    try:
        shutil.copy(os.path.join(FIXTURES, "RTprofile_in_CH4_realistic.nc"), deck)
        with nc.Dataset(deck, "a") as d:
            d["ch4mmr"][:] = 0.0
            d["n2mmr"][:] = 1.0
        with open(namelist, "w") as f:
            f.write(f"&exort_config\n  solar_file = '{STARS[star][0]}',\n"
                    "  shr_const_scon = 680.0,\n  exo_g = 9.80616,\n/\n")
        p = subprocess.run(["./exort.exe"], cwd=RUN_DIR, env=nf_env(),
                           capture_output=True, text=True)
        if p.returncode != 0:
            sys.exit(f"exort.exe failed for the N2 reference ({star}):\n{p.stdout[-2000:]}")
        shutil.copy(os.path.join(RUN_DIR, "RTprofile_out.nc"), out)
    finally:
        if saved is None:
            if os.path.exists(namelist):
                os.remove(namelist)
        else:
            with open(namelist, "w") as f:
                f.write(saved)
    return out


def read(path):
    with nc.Dataset(path) as d:
        return {v: np.array(d[v][:], dtype=float) for v in
                ("LWUP", "SWDN", "LWUP_SPECTRAL", "SWDN_SPECTRAL", "LWHR", "SWHR", "PMID")}


def load(gases, rerun):
    ref = {s: read(make_reference(s, rerun)) for s in STARS}
    runs = {}
    for g in gases:
        for lvl in ("realistic", "elevated"):
            for s in STARS:
                runs[(g, lvl, s)] = read(os.path.join(
                    BASELINES, f"RTprofile_out_{g}_{lvl}_{s}.nc"))
    return ref, runs


def grid():
    """Band edges and widths. (get_spectral_intervals returns edges and band
    MIDPOINTS; widths come from the edges. Wavelength edges are descending.)"""
    _, wn_edge, _, wl_edge, _ = get_spectral_intervals("n84")
    wn_edge, wl_edge = np.asarray(wn_edge), np.asarray(wl_edge)
    return wn_edge, np.diff(wn_edge), wl_edge, -np.diff(wl_edge)


def check_density(r, g_grid):
    """The spectral densities must integrate back to the broadband fluxes."""
    wn_edge, wn_diff, wl_edge, wl_diff = g_grid
    for dens, width, total, name in (
            (olr_density(r, wn_diff), wn_diff, r["LWUP"][0], "OLR"),
            (sws_density(r, wl_diff), wl_diff, r["SWDN"][-1], "surface SW down")):
        integ = float(np.sum(dens * width))
        if not np.isclose(integ, total, rtol=1e-4, atol=1e-3):
            sys.exit(f"{name}: spectral density integrates to {integ:.4f}, "
                     f"broadband is {total:.4f}")


def olr_density(r, wn_diff):
    return r["LWUP_SPECTRAL"][:, 0] / wn_diff            # W m-2 (cm-1)-1


def sws_density(r, wl_diff):
    return r["SWDN_SPECTRAL"][:, -1] / wl_diff           # W m-2 um-1


def steps(edges, values):
    """Step-plot coordinates for per-band values on band edges."""
    x = np.repeat(edges, 2)[1:-1]
    y = np.repeat(values, 2)
    return x, y


def style(ax, xlabel, ylabel):
    ax.grid(True, color=GRID, lw=0.6)
    ax.set_axisbelow(True)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_color("#c3c2b7")
    ax.tick_params(colors=MUTED, labelcolor="#52514e", labelsize=8)
    ax.set_xlabel(xlabel, color="#52514e", fontsize=9)
    ax.set_ylabel(ylabel, color="#52514e", fontsize=9)


def plot_gas(g, ref, runs, g_grid, outdir):
    wn_edge, wn_diff, wl_edge, wl_diff = g_grid
    real_lab, elev_lab = ABUND[g]
    series = [("N$_2$ only", ref, C_REF, "--"),
              (f"{GAS_LABEL[g]} {real_lab}", None, C_REAL, "-"),
              (f"{GAS_LABEL[g]} {elev_lab}", None, C_ELEV, "-")]

    def pick(i, star):
        if i == 0:
            return ref[star]
        return runs[(g, "realistic" if i == 1 else "elevated", star)]

    fig, axs = plt.subplots(2, 3, figsize=(15, 8.4))
    fig.patch.set_facecolor("#fcfcfb")

    # OLR spectrum (LW is star-independent here, so G2V only)
    ax = axs[0, 0]
    lw = wn_edge <= 3000
    nb = lw.sum() - 1
    for i, (lab, _, c, ls) in enumerate(series):
        x, y = steps(wn_edge[:nb + 1], olr_density(pick(i, "G2V"), wn_diff)[:nb])
        ax.plot(x, y, color=c, ls=ls, lw=1.6, label=lab,
                    zorder=3 if i == 0 else 2)
    style(ax, "wavenumber [cm$^{-1}$]", "OLR [W m$^{-2}$ (cm$^{-1}$)$^{-1}$]")
    ax.set_xlim(0, 3000)
    ax.set_title("Outgoing longwave (top of model)", fontsize=10, color=INK, loc="left")

    # Surface SW-down spectra, one panel per star
    sw = (wl_edge >= 0.2) & (wl_edge <= 5.0)
    for j, star in enumerate(STARS):
        ax = axs[0, 1 + j]
        idx = np.where(sw)[0]
        lo, hi = idx.min(), idx.max()
        for i, (lab, _, c, ls) in enumerate(series):
            vals = sws_density(pick(i, star), wl_diff)
            x, y = steps(wl_edge[lo:hi + 1], vals[lo:hi])
            ax.plot(x, y, color=c, ls=ls, lw=1.6, label=lab,
                    zorder=3 if i == 0 else 2)
        ax.set_xscale("log")
        ax.set_xlim(0.2, 5.0)
        ax.set_xticks([0.2, 0.3, 0.5, 1, 2, 3, 5])
        ax.set_xticklabels(["0.2", "0.3", "0.5", "1", "2", "3", "5"])
        style(ax, "wavelength [µm]", "surface SW↓ [W m$^{-2}$ µm$^{-1}$]")
        ax.set_title(f"Surface shortwave down — {STARS[star][1]}", fontsize=10,
                     color=INK, loc="left")

    # Heating rates vs pressure
    panels = [("LWHR", "G2V", "LW heating rate"),
              ("SWHR", "G2V", f"SW heating rate — {STARS['G2V'][1]}"),
              ("SWHR", "BB3400", f"SW heating rate — {STARS['BB3400'][1]}")]
    for j, (var, star, title) in enumerate(panels):
        ax = axs[1, j]
        for i, (lab, _, c, ls) in enumerate(series):
            r = pick(i, star)
            ax.plot(r[var], r["PMID"] / 100.0, color=c, ls=ls, lw=1.6, label=lab,
                    zorder=3 if i == 0 else 2)
        ax.set_yscale("log")
        ax.invert_yaxis()
        ax.axvline(0, color="#c3c2b7", lw=0.8)
        style(ax, "heating rate [K day$^{-1}$]", "pressure [hPa]")
        ax.set_title(title, fontsize=10, color=INK, loc="left")

    handles, labels = axs[0, 0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="upper right", ncol=3, frameon=False, fontsize=10,
               bbox_to_anchor=(0.985, 0.995))
    fig.suptitle(f"{GAS_LABEL[g]} in N$_2$ — TS273K column, ExoRT exort (HITRAN-2024)",
                 x=0.012, ha="left", y=0.985, fontsize=13, color=INK)
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    path = os.path.join(outdir, f"minor_gas_{g}.png")
    fig.savefig(path, dpi=150, facecolor=fig.get_facecolor())
    plt.close(fig)
    return path


def summary_rows(gases, ref, runs):
    rows = []
    for g in gases:
        for i, lvl in enumerate(("realistic", "elevated")):
            row = {"gas": g, "level": lvl, "vmr": ABUND[g][i]}
            r = runs[(g, lvl, "G2V")]
            row["olr"] = float(r["LWUP"][0])
            row["d_olr"] = float(r["LWUP"][0] - ref["G2V"]["LWUP"][0])
            for s in STARS:
                r = runs[(g, lvl, s)]
                row[f"sws_{s}"] = float(r["SWDN"][-1])
                row[f"d_sws_{s}"] = float(r["SWDN"][-1] - ref[s]["SWDN"][-1])
            rows.append(row)
    return rows


def plot_summary(rows, outdir):
    gases = list(dict.fromkeys(r["gas"] for r in rows))
    fig, axs = plt.subplots(1, 3, figsize=(15, 4.2), sharey=True)
    fig.patch.set_facecolor("#fcfcfb")
    y = np.arange(len(gases))
    h = 0.36
    panels = [("d_olr", "ΔOLR vs N$_2$ [W m$^{-2}$]"),
              ("d_sws_G2V", "Δ surface SW↓ vs N$_2$, G2V [W m$^{-2}$]"),
              ("d_sws_BB3400", "Δ surface SW↓ vs N$_2$, 3400 K [W m$^{-2}$]")]
    for ax, (key, title) in zip(axs, panels):
        for k, (lvl, c) in enumerate((("realistic", C_REAL), ("elevated", C_ELEV))):
            vals = [next(r[key] for r in rows if r["gas"] == g and r["level"] == lvl)
                    for g in gases]
            ax.barh(y + (k - 0.5) * h, vals, height=h * 0.9, color=c,
                    label=lvl, edgecolor="#fcfcfb", linewidth=1)
            for yy, v in zip(y + (k - 0.5) * h, vals):
                ax.text(v, yy, f" {v:+.2f} ", va="center", fontsize=7.5, color="#52514e",
                        ha="left" if v >= 0 else "right")
        ax.axvline(0, color="#c3c2b7", lw=0.8)
        style(ax, "W m$^{-2}$", "")
        ax.set_title(title, fontsize=10, color=INK, loc="left")
        lo, hi = ax.get_xlim()
        ax.set_xlim(lo * 1.25, max(hi, 0) + 0.05 * (hi - lo))
    axs[0].set_yticks(y)
    axs[0].set_yticklabels([GAS_LABEL[g] for g in gases], fontsize=10)
    axs[0].invert_yaxis()
    axs[-1].legend(frameon=False, fontsize=9, loc="lower left")
    fig.suptitle("Minor-gas regression cases: change relative to a pure-N$_2$ column (TS273K)",
                 x=0.012, ha="left", fontsize=12, color=INK)
    fig.tight_layout(rect=(0, 0, 1, 0.93))
    path = os.path.join(outdir, "minor_gas_summary.png")
    fig.savefig(path, dpi=150, facecolor=fig.get_facecolor())
    plt.close(fig)
    return path


def write_json(gases, ref, runs, rows, g_grid, outdir):
    wn_edge, wn_diff, wl_edge, wl_diff = g_grid
    rnd = lambda a, n=5: [float(f"{v:.{n}g}") for v in a]  # noqa: E731
    data = {"wn_edge": rnd(wn_edge), "wl_edge": rnd(wl_edge), "summary": rows,
            "abund": ABUND, "stars": {k: v[1] for k, v in STARS.items()}, "gases": {}}
    data["p_hpa"] = rnd(ref["G2V"]["PMID"] / 100.0)

    def pack(r):
        return {"olr": rnd(olr_density(r, wn_diff)),
                "sws": rnd(sws_density(r, wl_diff)),
                "lwhr": rnd(r["LWHR"], 4), "swhr": rnd(r["SWHR"], 4)}

    data["reference"] = {s: pack(ref[s]) for s in STARS}
    for g in gases:
        data["gases"][g] = {lvl: {s: pack(runs[(g, lvl, s)]) for s in STARS}
                            for lvl in ("realistic", "elevated")}
    path = os.path.join(outdir, "minor_gas_data.json")
    with open(path, "w") as f:
        json.dump(data, f, separators=(",", ":"))
    return path


def write_page(outdir):
    """Rebuild the interactive page (shared builder: tools/build_figure_pages.py)."""
    if os.path.abspath(outdir) != os.path.join(REPO, "figures_minor_gases"):
        return None   # the page is built from the repo copy of the data only
    import build_figure_pages as bfp
    links = {k: f"../{d}/{stem}.html" for k, (d, stem, _) in bfp.PAGES.items()}
    return bfp.build("minor", links, None)


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[1])
    ap.add_argument("--gases", nargs="+", default=GASES, choices=GASES)
    ap.add_argument("--outdir", default=os.path.join(REPO, "figures_minor_gases"))
    ap.add_argument("--rerun-reference", action="store_true",
                    help="recompute the cached N2-only reference outputs")
    args = ap.parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    ref, runs = load(args.gases, args.rerun_reference)
    g_grid = grid()
    for r in list(ref.values()) + list(runs.values()):
        check_density(r, g_grid)
    for g in args.gases:
        print("wrote", plot_gas(g, ref, runs, g_grid, args.outdir))
    rows = summary_rows(args.gases, ref, runs)
    print("wrote", plot_summary(rows, args.outdir))
    print("wrote", write_json(args.gases, ref, runs, rows, g_grid, args.outdir))
    page = write_page(args.outdir)
    if page:
        print("wrote", page)
    print(f"\n{'gas':5s} {'level':9s} {'VMR':>9s} {'ΔOLR':>8s} {'ΔSW↓ G2V':>9s} {'ΔSW↓ 3400K':>11s}")
    for r in rows:
        print(f"{r['gas']:5s} {r['level']:9s} {r['vmr']:>9s} {r['d_olr']:+8.2f} "
              f"{r['d_sws_G2V']:+9.2f} {r['d_sws_BB3400']:+11.2f}")


if __name__ == "__main__":
    main()
