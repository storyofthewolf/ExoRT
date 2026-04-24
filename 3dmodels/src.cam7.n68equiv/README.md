# src.cam7.n68equiv

## Overview
This directory contains the 12 files needed to link the ExoRT radiative transfer package with CAM7. It is a modernized update of the previous `src.cam.n68equiv` package (which worked with CESM1.2.1) and is designed to run with the development alpha version of **CESM3**.

Copies of the files are placed here instead of symbolic links so that one-off changes can be incorporated into 1D or GCM distributions without corrupting the other. These should be mostly identical to those found in `/ExoRT/source`, with a few specific exceptions required for 3D computational efficiency and CESM3 compatibility.

## NCAR CESM3-Planets Integration
ExoRT n68equiv is now natively integrated into the **CESM3-planets framework** as an external component. 

* **GitHub Repository**: [NCAR/CESM3-planets](https://github.com/NCAR/CESM3-planets)
* **Project Wiki**: [CESM3-planets Wiki](https://github.com/NCAR/CESM3-planets/wiki)

The current `CESM3-planets` port is in its early stages but successfully features a 2-degree Mars atmosphere and a CTSM land model running a very basic configuration. **We hope to entrain a larger community to improve and maintain this project!**

## Key Updates for CAM7 / CESM3-Planets
The following core modifications were made to the original ExoRT files to successfully transition from CESM1 to CESM3-planets:

* **`calc_opd_mod.F90`**: 
  * **Vacuum Safety Cutoff**: The CO2 k-distribution table can fail at extreme low pressures (returning unphysical opacities that crash the model). Radiation calculations are now safely skipped for pressures `< 0.05 mb` (5 Pa) to stabilize the model top.
  * Updated physical constants definitions to properly leverage CAM's `physconst` and `shr_const_mod`.
* **`exo_radiation_mod.F90`**: 
  * **Planetary Gravity Overrides**: Refactored the gravity constant imports. To correctly handle variable planetary gravity instead of Earth's default, `SHR_CONST_G` is now explicitly aliased to `gravit` from the `physconst` module, mapping the correct planetary acceleration to the radiation solver.
* **`exo_init_ref.F90`**: 
  * Cleaned up MPI integer/real imports from `spmd_utils` to align with the modernized CAM7 infrastructure.
* **`spectral_output_cam.F90`**: 
  * Updated the `cam_history` module imports, transitioning from the deprecated `phys_decomp` to `horiz_only` to properly handle historical field additions.
* **`mcica.F90`, `mcica_random_numbers.F90` & `model_specific.F90`**:
  * Cleaned up obsolete `kabs` imports and updated Fortran `iso_fortran_env` and integer types for the newer CESM3 build environments.

*(Note: As in previous versions, `calc_opd_mod.F90` and `exo_init_ref.F90` omit certain spectral loops and restrict integration ranges compared to the 1D model to maximize 3D efficiency. Consider these tunable parameters in your build).*

## Integration Instructions

Because ExoRT is now natively tracked by the `CESM3-planets` framework, **manual file copying to `SourceMods` or using the `-usr_src` flag is no longer required.** When you clone the `CESM3-planets` repository, the framework automatically pulls this code in as an external component. 

For complete instructions on checking out the model, synchronizing externals, and creating a Mars run, please refer to the **[CESM3-planets Quick Start Guide](https://github.com/NCAR/CESM3-planets/wiki/Quick-Start)**.
