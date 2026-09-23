#!/usr/bin/env python3
"""
build_figure_pages.py -- build the self-contained interactive figure pages.

Each page is a template (<dir>/<name>.tmpl.html) plus a data file
(<dir>/<name>_data.json). This script inlines the data, the shared styles and
chart code (figures_common/page.css, charts.js) and the cross-page nav, and
writes <dir>/<name>.html, which opens in any browser with no server.

  figures_regression/suite_page       v2 test-suite demonstration (data built here
                                      from tests/regression/baselines)
  figures_h2o_h16fix/hitran_page      H2O line-list progression and the HITRAN-2016
                                      temperature-index fix (data committed; its
                                      source runs predate this script)
  figures_minor_gases/minor_gas_page  CH4/C2H6/NH3/CO cases (data from
                                      tools/plot_minor_gas_cases.py)

Pages link to each other. By default the links are relative file paths; pass
--links links.json ({"suite": url, "hitran": url, "minor": url}) to point them at
published copies instead, and --outdir to write those copies elsewhere.

Usage
-----
  python tools/build_figure_pages.py
  python tools/build_figure_pages.py --links links.json --outdir /tmp/pages
"""

import argparse
import json
import os
import sys

import numpy as np
import netCDF4 as nc

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, os.pardir))
sys.path.insert(0, HERE)
from spectral_intervals import get_spectral_intervals  # noqa: E402

BASE = os.path.join(REPO, "tests", "regression", "baselines")

PAGES = {  # key: (directory, file stem, page title used in the nav)
    "suite": ("figures_regression", "suite_page", "v2 test suite"),
    "hitran": ("figures_h2o_h16fix", "hitran_page", "H₂O line lists"),
    "minor": ("figures_minor_gases", "minor_gas_page", "Minor gases"),
}


def rnd(a, n=5):
    return [float(f"{v:.{n}g}") for v in np.asarray(a, dtype=float).ravel()]


def read(case):
    with nc.Dataset(os.path.join(BASE, f"RTprofile_out_{case}.nc")) as d:
        return {v: np.array(d[v][:], dtype=float) for v in
                ("LWUP", "SWDN", "SWUP", "LWUP_SPECTRAL", "SWDN_SPECTRAL",
                 "LWHR", "SWHR", "PMID")}


def suite_data():
    """Data for the v2 test-suite page, from the committed baselines."""
    _, wn, _, wl, _ = get_spectral_intervals("n84")
    wn, wl = np.asarray(wn), np.asarray(wl)
    dwn, dwl = np.diff(wn), -np.diff(wl)

    def pack(r):
        olr = r["LWUP_SPECTRAL"][:, 0] / dwn
        sws = r["SWDN_SPECTRAL"][:, -1] / dwl
        for dens, width, total in ((olr, dwn, r["LWUP"][0]), (sws, dwl, r["SWDN"][-1])):
            if not np.isclose(np.sum(dens * width), total, rtol=1e-4, atol=1e-3):
                sys.exit("spectral density does not integrate to the broadband flux")
        return {"olr": rnd(olr), "sws": rnd(sws), "lwhr": rnd(r["LWHR"], 4),
                "swhr": rnd(r["SWHR"], 4), "p": rnd(r["PMID"] / 100.0),
                "OLR": float(r["LWUP"][0]), "SWS": float(r["SWDN"][-1]),
                "ALB": float(r["SWUP"][0] / r["SWDN"][0])}

    ts = [250, 273, 300, 320, 340, 360]
    seq = {star: {"OLR": [], "SWS": [], "ALB": []} for star in ("G2V", "BB3400")}
    for t in ts:
        for star in seq:
            r = read(f"TS{t}K_{star}")
            seq[star]["OLR"].append(float(r["LWUP"][0]))
            seq[star]["SWS"].append(float(r["SWDN"][-1]))
            seq[star]["ALB"].append(float(r["SWUP"][0] / r["SWDN"][0]))

    cases = sorted(f[len("RTprofile_out_"):-3] for f in os.listdir(BASE) if f.endswith(".nc"))
    table = []
    for c in cases:
        r = read(c)
        table.append({"case": c, "OLR": float(r["LWUP"][0]), "SWS": float(r["SWDN"][-1]),
                      "ALB": float(r["SWUP"][0] / r["SWDN"][0])})

    return {
        "wn_edge": rnd(wn), "wl_edge": rnd(wl), "ts": ts, "seq": seq, "table": table,
        "mars": {"dry": pack(read("2barCO2_dry_Mars_G2V")),
                 "cloud": pack(read("2barCO2_co2cloud_Mars_G2V"))},
        "haze": {"clear": pack(read("TS300K_G2V")),
                 "thin": pack(read("TS300K_haze_G2V")),
                 "thick": pack(read("TS300K_hazethick_G2V"))},
    }


def nav_html(current, links):
    items = []
    for key, (_, _, title) in PAGES.items():
        if key == current:
            items.append(f'<span aria-current="page">{title}</span>')
        else:
            items.append(f'<a href="{links[key]}">{title}</a>')
    return '<nav class="pages" aria-label="ExoRT figure pages">' + "".join(items) + "</nav>"


def build(key, links, outdir):
    d, stem, _ = PAGES[key]
    src = os.path.join(REPO, d)
    tmpl = os.path.join(src, f"{stem}.tmpl.html")
    data = os.path.join(src, f"{stem}_data.json") if key != "minor" else \
        os.path.join(src, "minor_gas_data.json")
    if not (os.path.exists(tmpl) and os.path.exists(data)):
        print(f"skip {key}: missing {tmpl if not os.path.exists(tmpl) else data}")
        return None
    with open(tmpl) as f:
        page = f.read()
    with open(data) as f:
        page = page.replace("/*DATA*/null", f.read())
    page = page.replace("<!--NAV-->", nav_html(key, links))
    common = os.path.join(REPO, "figures_common")
    with open(os.path.join(common, "page.css")) as f:
        page = page.replace("<!--COMMON_CSS-->", "<style>\n" + f.read() + "</style>")
    with open(os.path.join(common, "charts.js")) as f:
        page = page.replace("/*COMMON_JS*/", f.read())
    dest = os.path.join(outdir or src, f"{stem}.html")
    os.makedirs(os.path.dirname(dest), exist_ok=True)
    with open(dest, "w") as f:
        if outdir:   # publishing copy: the host supplies the document skeleton
            f.write(page)
        else:
            f.write('<!doctype html>\n<html lang="en"><head><meta charset="utf-8">'
                    '<meta name="viewport" content="width=device-width,initial-scale=1">'
                    "</head><body>\n" + page + "\n</body></html>\n")
    return dest


def main():
    ap = argparse.ArgumentParser(description="Build the ExoRT interactive figure pages.")
    ap.add_argument("--links", help="JSON mapping page key -> URL for the cross-links")
    ap.add_argument("--outdir", help="write the pages here instead of next to their data")
    args = ap.parse_args()

    os.makedirs(os.path.join(REPO, "figures_regression"), exist_ok=True)
    with open(os.path.join(REPO, "figures_regression", "suite_page_data.json"), "w") as f:
        json.dump(suite_data(), f, separators=(",", ":"))

    if args.links:
        with open(args.links) as f:
            links = json.load(f)
    else:  # relative links between the repo's figure directories
        links = {k: f"../{d}/{stem}.html" for k, (d, stem, _) in PAGES.items()}
    for key in PAGES:
        path = build(key, links, args.outdir)
        if path:
            print("wrote", path)


if __name__ == "__main__":
    main()
