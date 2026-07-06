#!/bin/bash
# All supported CPP combinations of the src.cam.exort bundle.
set -e
HERE=$(cd "$(dirname "$0")" && pwd)
"$HERE/build.sh"
"$HERE/build.sh" -DEXORT_CO2CLD
"$HERE/build.sh" -DEXORT_CARMA
"$HERE/build.sh" -DEXORT_SRF_EMISS
"$HERE/build.sh" -DEXORT_CO2CLD -DEXORT_CARMA -DEXORT_SRF_EMISS
echo "ALL OK"
