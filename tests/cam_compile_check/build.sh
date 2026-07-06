#!/bin/bash
# Compile-check 3dmodels/src.cam.exort against stubbed CESM modules.
#
# This is the local (no-HPC) compile gate for the CAM bundle: it catches
# use-graph breaks, signature mismatches, and CPP-guard errors without a
# CESM checkout. It does NOT link or run anything — real validation is an
# ExoCAM build on the target machine.
#
# Usage:
#   ./build.sh                       # mainline (no feature macros)
#   ./build.sh -DEXORT_CO2CLD -DEXORT_CARMA -DEXORT_SRF_EMISS
#   ./run_all.sh                     # all four supported combinations
#
# exoplanet_mod_stub.F90 mirrors ExoCAM's src.share/exoplanet_mod.F90 and
# doubles as the template for the ExoCAM-side additions the bundle needs
# (do_exo_clouds / do_exo_haze; see the bundle README).
set -e
HERE=$(cd "$(dirname "$0")" && pwd)
RT=$(cd "$HERE/../.." && pwd)
MISC=$RT/source/src.misc
BUN=$RT/3dmodels/src.cam.exort
WORK=$HERE/work$(echo "$@" | tr -cd 'A-Z_')
rm -rf "$WORK"; mkdir -p "$WORK"; cd "$WORK"

FC="gfortran -c -cpp -I$MISC -I$HERE \
    -fallow-invalid-boz -fno-range-check -fallow-argument-mismatch \
    -ffree-line-length-none -w $*"

# infrastructure: real src.misc shims where possible, stubs otherwise
$FC $MISC/shr_kind_mod.F90
$FC $HERE/exoplanet_mod_stub.F90
$FC $MISC/shr_const_mod.F90
$FC $MISC/ppgrid.F90
$FC $MISC/pmgrid.F90
$FC $MISC/spmd_utils.F90
$FC $MISC/infnan.F90
$FC $MISC/physconst.F90
$FC $MISC/ioFileMod.F90
$FC $HERE/cesm_stubs.F90

# the bundle, in dependency order
$FC $BUN/sys_rootdir.F90
$FC $BUN/rayleigh_data.F90
$FC $BUN/radgrid.F90
$FC $BUN/cloud.F90
$FC $BUN/kabs.F90
$FC $BUN/rad_interp_mod.F90
$FC $BUN/mcica_random_numbers.F90
$FC $BUN/mcica.F90
$FC $BUN/planck_mod.F90
$FC $BUN/exo_init_ref.F90
$FC $BUN/model_specific.F90
$FC $BUN/calc_opd_mod.F90
$FC $BUN/exo_radiation_mod.F90
$FC $BUN/spectral_output_cam.F90
$FC $BUN/initialize_rad_mod_cam.F90
$FC $BUN/exo_radiation_cam_intr.F90
echo "OK: src.cam.exort compiled with flags: $*"
