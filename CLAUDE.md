# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What This Is

ExoRT is a two-stream radiative transfer code for exoplanet climate modeling. It supports both 1-D offline calculations and 3-D integration with CESM/CAM4 (via ExoCAM). The codebase is primarily Fortran 90, with Python and IDL tools for pre/post-processing.

## Build Commands

All builds run from `build/`:

```bash
cd build

make n68equiv        # recommended version for terrestrial planets
make n84equiv        # for F-stars (>6500 K)
make n28archean      # Archean climate simulations
make n42h2o          # H2O-only, 42 bins
make n68h2o          # H2O-only, 68 bins
make n68equiv_exp    # experimental build

make clean           # remove all build artifacts and run/*.exe
```

The compiler defaults to `ifort`; on Apple Silicon Macs use `USER_FC=gfortran make n68equiv` (ifort has no arm64 port). Requires NetCDF4 Fortran library (`nf-config` must be on PATH). Executables are placed in `run/`.

## Running the 1-D Model

```bash
# 1. Set exo_pver in source/exoplanet_mod.F90 to match your input file (compile-time only)
# 2. Build: cd build && USER_FC=gfortran make n68equiv
# 3. Place input file in run/ as RTprofile_in.nc
# 4. Optionally create run/user_nl_exort to override solar_file, shr_const_scon, exo_g
cd run
./n68equiv.exe
# Output written to run/RTprofile_out.nc
```

Template input files are in `iofiles/input_files/`. The input NetCDF must match the vertical level count set in `exoplanet_mod.F90` (`exo_pver`).

### Runtime namelist: `user_nl_exort`

Three parameters that were formerly compile-time constants are now overridable at runtime without a rebuild:

| Variable | Meaning | Default |
|----------|---------|---------|
| `solar_file` | Stellar spectrum filename in `data/solar/` | `'BT_SETTL_interp_t2550_logg5.0_m0.0_a0.0.nc'` |
| `shr_const_scon` | Stellar constant ÷ 2 [W m⁻²] | `680.0` (≈ present Earth) |
| `exo_g` | Surface gravity [m s⁻²] | `0.93 × 9.80616` (≈ 9.12) |

**To use:** copy the template, edit it, and place it in the run directory before invoking the executable.

```bash
cp iofiles/user_nl_exort.template run/user_nl_exort
# edit run/user_nl_exort, then:
cd run && ./n68equiv.exe
```

If `user_nl_exort` is absent the model runs with the defaults above — bit-for-bit identical to the old hardcoded behaviour. The active values are always printed at startup under the `=== exort_config ===` header.

**Namelist syntax** (standard Fortran):
```fortran
&exort_config
  solar_file     = 'G2V_SUN_n68.nc',
  shr_const_scon = 680.0,
  exo_g          = 9.80616,
/
```

Only variables you want to change need to appear. `exo_pver` is **not** in the namelist — it is still compile-time only and requires a rebuild.

## Architecture

### Source Layout

```
source/
  exoplanet_mod.F90          # Top-level config: solar file, gravity, pver, scon
  src.main/                  # Shared drivers (used by all RT versions)
  src.misc/                  # CESM stubs required for standalone operation
  src.n68equiv/              # Version-specific RT code (recommended)
  src.n84equiv/
  src.n28archean/
  src.n42h2o/
  src.n68h2o/
  experimental/
3dmodels/src.cam.*/          # CAM hook-up copies (duplicates of source/, not symlinks)
```

### Key Module Responsibilities

Each RT version (`src.n*`) provides these version-specific files:
- `radgrid.F90` — spectral interval and Gauss point definitions
- `kabs.F90` — k-distribution absorption lookup
- `calc_opd_mod.F90` — optical depth calculation
- `rad_interp_mod.F90` — interpolation in P-T space
- `cloud.F90` — cloud optical properties
- `model_specific.F90` — gas species configuration
- `initialize_rad_mod_1D.F90` — 1-D initialization

Shared across all versions (`src.main/`):
- `exo_radiation_mod.F90` — top-level RT driver (two-stream solver)
- `exo_init_ref.F90` — reference atmosphere initialization
- `main.F90` — 1-D standalone entry point
- `input.F90` / `output.F90` — NetCDF I/O
- `mcica.F90` — Monte Carlo Independent Column Approximation for clouds
- `planck_mod.F90` — Planck function

### Configuration in `exoplanet_mod.F90`

This is the primary file to edit for each run. Key parameters:
- `solar_file` — stellar spectrum NetCDF from `data/solar/`
- `exo_pver` — number of vertical levels (must match input file)
- `exo_g` — surface gravity (m/s²)
- `shr_const_scon` — stellar constant (W/m², divided by 2 for 1-D)
- `do_exo_rt_optimize_bands` — enable/disable band optimization

### 3-D Model Integration (ExoCAM)

`3dmodels/` contains file copies (not symlinks) for CAM integration. To use with ExoCAM:

```bash
xmlchange CAM_CONFIG_OPTS="-usr_src /$MYDIR/ExoRT/3dmodels/src.cam.n68equiv"
```

Files in `3dmodels/` may diverge slightly from `source/` to accommodate 1-D vs. 3-D differences.

## Data Files

| Directory | Contents |
|-----------|----------|
| `data/solar/` | Stellar spectra NetCDF files (named `*_n68.nc`, `*_n42.nc`, etc., matching RT version) |
| `data/kdist/` | Correlated k-distribution tables |
| `data/cia/` | Collision-induced absorption (N₂–N₂, N₂–H₂, H₂–H₂, CO₂–CO₂, etc.) |
| `data/continuum/` | H₂O and CO₂ continuum from MT_CKD/LBLRTM |
| `data/cloud/` | Mie scattering cloud optical properties |
| `data/aerosol/` | CARMA aerosol optical constants |

Solar spectrum filenames encode the RT version (e.g., `G2V_SUN_n68.nc` for n68equiv; `LHS1140_spectra_n42.nc` for n42h2o).

## Tools

`tools/` contains pre/post-processing scripts. Most are IDL (`.pro`); Python equivalents exist for the primary workflow:

### Python tools (preferred for new work)

- **`profile_data.py`** — P-T profile data library. Contains `smart_2bar_t250`, `TS273K`, and `US1976` profiles as numpy arrays. Use `get_profile(tag)` to retrieve a profile dict with keys `pint`, `pmid`, `tint`, `tmid`, `q`. Equivalent to `profile_data.pro`.
- **`makeColumn.py`** — Build `RTprofile_in.nc` from scratch. Edit the `USER SETTINGS` block at the top, or use CLI flags. All gas species are optional; absent species are simply omitted from the output file and ExoRT reads them as zero. Equivalent to `makeColumn.pro`.
  ```bash
  python makeColumn.py                        # use USER SETTINGS block
  python makeColumn.py --defaults             # print current defaults
  python makeColumn.py --profile TS273K --co2vmr 400e-6 --nh3vmr 1e-6 --output RTprofile_in.nc
  python makeColumn.py --zero-h2o            # force H2O to zero
  ```
- **`plotspectra_1D.py`** — Plot shortwave and longwave spectral fluxes from `RTprofile_out.nc`. SW panel uses wavelength axis (µm, W m⁻² µm⁻¹); LW panel uses wavenumber axis (cm⁻¹, W m⁻² cm). Supports two-file comparison. Equivalent to `plotspectra_1D.pro`.
  ```bash
  python plotspectra_1D.py                    # use USER SETTINGS block
  python plotspectra_1D.py --save             # save PNG instead of display
  python plotspectra_1D.py --f2 RTprofile_out_ref.nc --label2 "no NH3"
  ```
- **`spectral_intervals.py`** — Spectral grid edge arrays and Gauss weights for all RT versions. Used by `plotspectra_1D.py`; call `get_spectral_intervals(tag)` with tag `n28`, `n42`, `n68`, or `n84`.
- **`makeStellarSpectrum_blackbody.py`** — Generate stellar spectrum NetCDF from blackbody temperature.
- **`check_kcoeff.py`** — Inspect k-coefficient files.
- **`populate3Dmodels.py`** — Sync files between `source/` and `3dmodels/`.

### IDL tools (legacy)

- `makeColumn.pro` / `profile_data.pro` — original column builder and profile data
- `plotspectra_1D.pro` — original spectral plotting
- `makeStellarSpectrum_blackbody.pro`, `plotprofile.pro`, `plotstellar.pro`, etc.

## Input File Behavior (`input.F90`)

All gas species variables in `RTprofile_in.nc` are **optional**. `input.F90` uses `nf_inq_varid` directly (bypassing `wrap_inq_varid`) for every gas and cloud variable so that missing variables produce a clean diagnostic line rather than a NetCDF error dump. The output at startup looks like:

```
--- gas species in input file ---
  h2ommr:  found
  co2mmr:  found
  ch4mmr:  not found, set to zero
  nh3mmr:  found
  ...
--- cloud variables in input file ---
  cicewp:  not found, set to zero
  ...
```

Only the P/T/Z arrays, albedos, `coszrs`, `mw`, and `cp` are required. When adding a new gas to `input.F90`, always use the optional pattern — never `wrap_inq_varid` for a gas species.

## .gitignore Notes

Patterns in `.gitignore` are relative to the repo root. The prefix `/ExoRT/` does **not** work — it would only match if the repo were checked out inside a directory literally named `ExoRT`. Use bare relative paths:

```
tools/*.nc       # correct
/ExoRT/tools/*.nc  # wrong — silently ignored by git
```

## RT Version Selection Guide

- **Terrestrial planets, most cases:** `n68equiv` (68 bins, H₂O/CO₂/CH₄/C₂H₆/NH₃/CO, HITRAN 2016/2020)
- **F-stars (6500–10000 K):** `n84equiv` (adds UV bins below 0.24 µm)
- **Archean Earth:** `n28archean` (28 bins, includes CH₄ near-IR caution)
- **H₂O/N₂/H₂ atmospheres only:** `n42h2o` or `n68h2o`

## Adding a New Gas Absorber (n68equiv)

NH₃ and CO were added to `n68equiv` on 2026-04-27 as an example. To add another gas, touch these 9 files in this order:

1. **`source/src.misc/physconst.F90`** — add `mwXXX` molecular weight parameter
2. **`source/src.n68equiv/kabs.F90`** — add `dirk_xxx` (data directory) and `k_xxx_file` (filename) parameters
3. **`source/src.n68equiv/radgrid.F90`** — increment `nspecies`, add `iXXX` index parameter, extend `gas_name`, add `k_xxx` array to the individual k-array declaration
4. **`source/src.n68equiv/initialize_rad_mod_1D.F90`** — add a `getfil` + `wrap_open` + `wrap_get_var_realx` block to read the k-coefficient file
5. **`source/src.n68equiv/model_specific.F90`** — add one line each to `k_major_data` and `k_grey_data` assignments
6. **`source/src.main/input.F90`** — add `XXXMMR_in` declaration and zero-init; use `nf_inq_varid` directly (not the wrapper) so the variable is optional: print a clean found/not-found message and call `wrap_get_var_realx` only if the variable is present
7. **`source/src.main/exo_radiation_mod.F90`** — add `ext_XXX` to `aerad_driver` signature and intent declaration; add `qXXX` local array; assign pseudo-layer and mid-layer values; pass to `calc_gasopd`
8. **`source/src.main/main.F90`** — add `XXXMMR_in` to the `aerad_driver` call
9. **`source/src.n68equiv/calc_opd_mod.F90`** — add `mwxxx` to physconst import; add `qxxx` argument and intent declaration; add `xxxvmr` and `u_xxx` local variables; compute VMR and column density in the level loop; extend `ugas` array; add `bilinear_interpK_grey` call and `tau_grey` entry

**Data files:** k-coefficient NetCDF files for new gases go in `data/kdist/n68xxx/` and must follow the same P-T grid as existing species (`kc_npress=61`, `kc_ntemp=17`, `ngauss_8gpt=8`). The variable name inside the NetCDF must be `data`. NH₃ files go in `data/kdist/n68nh3/`; CO files go in `data/kdist/n68co/`. Input profile NetCDFs must include the corresponding `xxxmmr` variable (e.g. `nh3mmr`, `commr`).

## Connecting a New Gas to the 3-D Interface (n68equiv + ExoCAM)

After completing the 1-D steps above, sync the 3-D files and extend the CAM interface. NH₃ and CO were connected on 2026-04-29 as the reference example.

**Files to copy from source to `3dmodels/src.cam.n68equiv/`** (should be identical):

- `source/src.n68equiv/radgrid.F90`
- `source/src.n68equiv/kabs.F90`
- `source/src.n68equiv/model_specific.F90`
- `source/src.n68equiv/calc_opd_mod.F90`
- `source/src.main/rayleigh_data.F90`
- `source/src.main/exo_radiation_mod.F90`

**Files requiring new edits in both `source/src.n68equiv/` and `3dmodels/src.cam.n68equiv/`** (keep in sync):

10. **`source/src.n68equiv/initialize_rad_mod_cam.F90`** — add a PIO file-read block for the new k-file (same pattern as the existing species: `getfil` + `cam_pio_openfile` + `pio_inq_varid` + `pio_get_var` + `pio_closefile`); add a matching `mpibcast` call in the `#if defined SPMD` block.

11. **`source/src.main/exo_radiation_cam_intr.F90`** — add `mwXXX` to the `physconst` import; add `exo_xxxmmr` to the `exoplanet_mod` import; declare a `xxxmmr(pcols,pver)` local array; assign `xxxmmr(:,:) = exo_xxxmmr` in the gas-retrieval block; add `xxxmmr(i,:)` to both `aerad_driver` call sites (clearsky and full-sky).

**Files to update in ExoCAM** (in all active configs under `cesm1.2.1/configs/*/SourceMods/`):

12. **`src.cam/physconst.F90`** — add `mwXXX` molecular weight parameter; add `write(iulog,*)` diagnostic line for VMR/MMR after the last existing gas line.

13. **`src.share/exoplanet_mod.F90`** — add `exo_xxxbar = 0.0_r8` inventory parameter; add `cpXXX` specific heat parameter; derive `exo_xxxvmr` and `exo_xxxmmr` from them; add `exo_xxxmmr*cpXXX` to the `exo_cpdair` formula.

**Sync rule:** After all edits, `diff source/src.n68equiv/<file> 3dmodels/src.cam.n68equiv/<file>` must produce no output for every shared file except `sys_rootdir.F90` (which intentionally differs by machine path).
