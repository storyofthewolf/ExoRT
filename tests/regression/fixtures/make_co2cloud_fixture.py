#!/usr/bin/env python3
"""Generate the CO2-ice cloud regression fixture for the Stage C cloud path.

Clones the committed 2-bar pure-CO2 Mars column
(RTprofile_in_2barCO2_dry_pver300.nc, pver=300) and layers a CO2 ice cloud into
the mid-troposphere, plus a sub-unity surface emissivity. This exercises the
C1b CO2-ice cloud path and the C1-emiss srf_emiss input (both require
do_exo_clouds=.true. in ExoRT for the cloud part).

Cloning the existing 300-level fixture (rather than rebuilding from the native
69-level smart_2bar_t250 profile) keeps the vertical grid identical to the
clear Mars baseline, so the cloudy-vs-clear comparison isolates the cloud
effect rather than a regrid difference.

Run from this directory:
    python make_co2cloud_fixture.py
writes RTprofile_in_2barCO2_co2cloud_pver300.nc here.
"""
import os
import shutil

import numpy as np
import netCDF4 as nc

HERE = os.path.dirname(os.path.abspath(__file__))
SRC = os.path.join(HERE, "RTprofile_in_2barCO2_dry_pver300.nc")
OUT = os.path.join(HERE, "RTprofile_in_2barCO2_co2cloud_pver300.nc")

# CO2 ice cloud: a ~20-layer slab in the mid column.
KLO, KHI = 140, 160
CICEWP_CO2 = 5.0    # g/m2 in-cloud ice path per layer
REI_CO2 = 10.0      # micron effective radius
SRF_EMISS = 0.95    # slightly sub-blackbody regolith


def main():
    shutil.copyfile(SRC, OUT)
    with nc.Dataset(OUT, "a") as d:
        nlev = len(d.dimensions["pver"])
        cicewp = np.zeros(nlev, "f4"); cicewp[KLO:KHI] = CICEWP_CO2
        rei = np.zeros(nlev, "f4");    rei[KLO:KHI] = REI_CO2
        for name, arr, title, units in (
            ("cicewp_co2", cicewp, "in-cloud CO2 ice water path", "g/m2"),
            ("rei_co2", rei, "CO2 ice cloud effective radius", "microns"),
        ):
            v = d.variables[name] if name in d.variables \
                else d.createVariable(name, "f4", ("pver",))
            v.title = title; v.units = units
            v[:] = arr
        v = d.variables["srf_emiss"] if "srf_emiss" in d.variables \
            else d.createVariable("srf_emiss", "f4", ("one",) if "one" in d.dimensions else ())
        v.title = "Surface thermal emissivity"; v.units = "1"
        v[...] = SRF_EMISS
    print(f"CO2-cloud fixture: {OUT}")
    print(f"  cloud layers {KLO}..{KHI}: cicewp_co2={CICEWP_CO2} g/m2, "
          f"rei_co2={REI_CO2} um; srf_emiss={SRF_EMISS}")


if __name__ == "__main__":
    main()
