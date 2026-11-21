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
#
# shr_const_mod and physconst are what ExoCAM provides, not ExoRT's 1-D
# shims: shr_const_mod_exocam.F90 is a verbatim copy of mainline ExoCAM's
# (cam_aqua_fv) and physconst_stub.F90 exports only names ExoCAM's
# physconst has. A bundle import that ExoCAM lacks therefore fails here.
# (The src.misc shr_const_mod carries ExoRT-only constants; compiling
# against it hid SHR_CONST_DENSITYCO2ICE until the first real ExoCAM build.)
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
# ExoCAM's real shr_const_mod (vendored), NOT source/src.misc's: the bundle must
# compile against the constants mainline ExoCAM actually provides.
$FC $HERE/shr_const_mod_exocam.F90
$FC $MISC/ppgrid.F90
$FC $MISC/pmgrid.F90
$FC $MISC/spmd_utils.F90
$FC $MISC/infnan.F90
$FC $HERE/physconst_stub.F90      # CAM's physconst names, not ExoRT's 1-D shim
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
