# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What This Is

ExoRT is a two-stream radiative transfer code for exoplanet climate modeling. It supports both 1-D offline calculations and 3-D integration with CESM/CAM4 (via ExoCAM). The codebase is primarily Fortran 90, with Python and IDL tools for pre/post-processing.

## Build Commands

All builds run from `build/`:

```bash
cd build

make exort           # v2 single bundle (84-band, HITRAN-2024, NH3/CO) — primary target
make libexort        # Stage D shared library (run/libexort.dylib|.so) with the C API
make n68equiv        # legacy HITRAN-2016 reference (slated for retirement)
make n84equiv        # legacy HITRAN-2016 reference, +UV bins (slated for retirement)

make clean           # remove all build artifacts, run/*.exe, run/libexort.*
```

**Library API (Stage D):** `libexort` exposes `exort_init` / `exort_run_column`
/ `exort_run_columns` / `exort_finalize` / `exort_get_dims`
(`source/src.main/exort_lib_mod.F90`; structs in `exort_column_mod.F90`).
Python binding + usage in `tools/exort_pytools/README.md`; acceptance checks:
`python tools/exort_pytools/verify_lib.py` and `cd tests/lib && make run`.
To reproduce the committed baselines through the library, pass the `_n84`
stellar files (`G2V_SUN_n84.nc`, `blackbody_3400K_n84.nc`) — the regression
harness maps `_n68 → _n84` for exort runs.

**v2 build status (refactor branch):** the legacy `n28archean` / `n42h2o` /
`n68h2o` targets were **removed in v2** (they live in the `v1.0.0` tag). `exort`
is the v2 path forward; `n68equiv`/`n84equiv` remain only as HITRAN-2016
comparison references. Note `n68equiv`/`n84equiv` now have their `kabs.F90` data
paths hand-edited to the flat `data/kdist/<gas>/` layout (test scaffolding) and
will not run against the old `data/kdist/n68*` tree.

✅ **`src.exort` runs on the validated HITRAN-2016 native-gas k-files by default**
(as of 2026-06-28). The structural refactor is decoupled from the HITRAN-2024
line-list upgrade: `kabs.F90` pins H₂O/CO₂/CH₄/C₂H₆ to `hitran16` (NH₃/CO are
`hitran24`-only but proven clean; O₂/O₃ are `hitran20`). The regression suite
baselines against this build. The HITRAN-2024 tables remain **unvalidated**
(CO₂ far-IR χ-factor pipeline bug; H₂O possibly partly real) and are reachable
only via the `run_regression.py --exort h24` side-path. See
`tests/regression/EXORT_H16_EQUIVALENCE.md`, `gas_sweep.py`, and `REFACTOR_LOG.md`.

The compiler defaults to `ifort`; on Apple Silicon Macs use `USER_FC=gfortran make exort` (ifort has no arm64 port). Requires NetCDF4 Fortran library (`nf-config` must be on PATH). Executables are placed in `run/`.

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

Template input files are in `iofiles/input_files/`. The input NetCDF must match the vertical level count set in `exoplanet_mod.F90` (`exo_pver`). The input may carry any number of columns via an optional `ncol` dimension (see "Input File Behavior" below); one invocation then solves all columns.

### Runtime namelist: `user_nl_exort`

Three parameters that were formerly compile-time constants are now overridable at runtime without a rebuild:

| Variable | Meaning | Default |
|----------|---------|---------|
| `solar_file` | Stellar spectrum filename in `data/stellar/` | `'G2V_SUN_n84.nc'` (v2; was a BT_SETTL file in v1) |
| `shr_const_scon` | Stellar constant ÷ 2 [W m⁻²] | `680.0` (≈ present Earth) |
| `exo_g` | Surface gravity [m s⁻²] | `9.80616` (Earth) |
| `do_exo_clouds` | Enable the cloud RT path (H₂O + CO₂ ice; reads `cicewp*`/`rei*` from the input file) | `.false.` |
| `do_exo_haze` | Enable the CARMA haze aerosol RT path (reads `carmammr(pver,nelem,nbin)` from the input file; optics from `data/aerosol/haze_n84_b40_*.nc`) | `.false.` |
| `mcica_percol_seed` | Opt-in per-column MCICA seed: each batch column offsets the stochastic-cloud seed by its column index (column 1 stays bit-identical to legacy). Off = constant seed 9404 for every column. Only affects cloudy H₂O runs; enabling is a rebaseline decision | `.false.` |

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
  do_exo_clouds  = .false.,
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
  src.exort/                 # v2 single RT bundle (primary)
  src.n68equiv/              # legacy HITRAN-2016 reference (slated for retirement)
  src.n84equiv/              # legacy HITRAN-2016 reference, +UV bins (slated for retirement)
3dmodels/src.cam.exort/      # v2 CAM hook-up (byte-synced with source/; CO2 clouds/haze/emiss via -cppdefs)
3dmodels/src.cam.*/          # pre-v2 CAM hook-ups (frozen legacy copies)
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
- `main.F90` — 1-D standalone entry point (serial loop over columns)
- `io_1D.F90` — NetCDF I/O (module `io`; holds no column state since E1)
- `exort_column_mod.F90` — bind(c) `column_state_t`/`column_result_t` structs
- `exort_column_run.F90` — `run_one_column`, the single shared per-column
  solve used by both `main.F90` and the library (`exort_lib_mod.F90`)
- `mcica.F90` — Monte Carlo Independent Column Approximation for clouds
- `planck_mod.F90` — Planck function

### Configuration in `exoplanet_mod.F90`

This is the primary file to edit for each run. Key parameters:
- `solar_file` — stellar spectrum NetCDF from `data/stellar/`
- `exo_pver` — number of vertical levels (must match input file)
- `exo_g` — surface gravity (m/s²)
- `shr_const_scon` — stellar constant (W/m², divided by 2 for 1-D)
- `do_exo_rt_optimize_bands` — enable/disable band optimization

### 3-D Model Integration (ExoCAM)

`3dmodels/` contains file copies (not symlinks) for CAM integration. To use with ExoCAM:

```bash
# v2 (preferred; see 3dmodels/src.cam.exort/README for ExoCAM-side requirements)
xmlchange CAM_CONFIG_OPTS="-usr_src /$MYDIR/ExoRT/3dmodels/src.cam.exort"
# optional physics: -cppdefs '-DEXORT_CO2CLD -DEXORT_CARMA -DEXORT_SRF_EMISS'

# pre-v2 legacy bundles (frozen)
xmlchange CAM_CONFIG_OPTS="-usr_src /$MYDIR/ExoRT/3dmodels/src.cam.n68equiv"
```

`src.cam.exort` is kept byte-identical to `source/` (gate:
`python tools/populate3Dmodels.py check --exort`; local compile gate:
`tests/cam_compile_check/run_all.sh`). The pre-v2 bundles intentionally
diverge from `source/` and must not be regenerated.

## Data Files

| Directory | Contents |
|-----------|----------|
| `data/stellar/` | Stellar spectra NetCDF files (renamed from `data/solar/` in v2; named `*_n84.nc`, `*_n68.nc`, etc., matching the RT grid) |
| `data/kdist/` | Correlated k-distribution tables. **v2 layout: flat per-gas dirs `data/kdist/<gas>/`** (legacy `n68<gas>/` tree removed). Each gas dir holds both `n84_8gpt_<gas>_hitran16…` and `…hitran24…` files keyed by filename prefix. |
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

## Input File Behavior (`io_1D.F90`, module `io`)

All gas species variables in `RTprofile_in.nc` are **optional**. `input_profile` uses `nf_inq_varid` directly (bypassing `wrap_inq_varid`) for every gas and cloud variable so that missing variables produce a clean diagnostic line rather than a NetCDF error dump. The output at startup looks like:

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

Only the P/T/Z arrays, albedos, `coszrs`, `mw`, and `cp` are required. When adding a new gas, always use the optional pattern (the `opt_mid` helper inside `input_profile`) — never `wrap_inq_varid` for a gas species.

**Multi-column input (Stage E1):** `RTprofile_in.nc` may carry an optional `ncol` dimension. Absent = classic single-column file (bit-for-bit the legacy path). Present = every variable carries a trailing column dimension in Fortran order (`tmid(pver,ncol)`, `ts(ncol)`, `carmammr(pver,nelem,nbin,ncol)`), columns are solved in an **OpenMP-parallel loop** (Stage E2; threads spawn only when `ncol > 1`, count from `OMP_NUM_THREADS`, results bitwise independent of thread count), and `RTprofile_out.nc` mirrors the `ncol` dimension. Build multi-column inputs with `python tools/stackColumns.py col1.nc col2.nc -o RTprofile_in.nc`. Columns share the process-level `solar_file` and cloud/haze flags, but **gravity and insolation may vary per column** (Stage E2b): optional scalar input variables `grav`/`scon` override the namelist `exo_g`/`shr_const_scon` for that column (absent or ≤0 = namelist value; `makeColumn.py --write-grav` / `--scon` writes them). H₂O cloud fraction is drivable from the deck via the optional `cfrc` variable (E2; absent = zero). Acceptance checks: `tests/regression/multicol_check.py` (batch == singles exactly; threaded == serial exactly) and `tests/regression/percol_config_check.py` (per-column grav exact; per-column scon ≤1e-10 rel vs namelist-configured singles). User-facing recipe: `MULTICOLUMN_BATCH_HOWTO.md`.

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
6. **`source/src.main/exort_column_mod.F90` + `io_1D.F90`** — (post-E1 form of the old `input.F90` step) add a `xxxmmr(pver)` field to `column_state_t` (and zero it in `zero_column_state`); in `input_profile` read it via the optional `opt_mid` pattern; add the output variable in `output_data`
7. **`source/src.main/exo_radiation_mod.F90`** — add `ext_XXX` to `aerad_driver` signature and intent declaration; add `qXXX` local array; assign pseudo-layer and mid-layer values; pass to `calc_opd_gas`
8. **`source/src.main/exort_column_run.F90`** — (post-E1 form of the old `main.F90` step) add `state%xxxmmr` to the `aerad_driver` call in `run_one_column`
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

## Regression Test Suite

A Python regression harness lives in `tests/regression/` (see its README). It runs
a set of standard profiles through `run/n68equiv.exe` and compares flux/heating
outputs against committed golden baselines with rtol/atol = 1e-3.

```bash
cd tests/regression
python run_regression.py                 # run all cases, compare to baselines
python run_regression.py --list          # list case names + their physics
python run_regression.py --cases TS300K Mars   # subset by substring
python run_regression.py --generate-baselines  # (re)create golden baselines
```

Cases are defined in `build_cases()`; each carries its own fixture, stellar
spectrum, insolation (`shr_const_scon`), and gravity (`exo_g`), so heterogeneous
planets (e.g. the Mars-like `2barCO2_dry_Mars_G2V`, g=3.711) coexist with the
Earth-like TS250K–TS360K × {G2V, blackbody_3400K} sequence. 15 cases total
(12 clear TS + clear Mars + gated CO₂-cloud Mars + gated hazy TS300K), all
`pver=300`. The harness auto-sets the NetCDF lib path for the macOS loader and
preserves/restores any existing `run/user_nl_exort`. Use this to verify any code
change is bit-for-bit (Δ=0) or to gate intended physics changes — see
`REFACTOR_PLAN.md` for the rebaseline workflow.

Two inline glance columns (`LWUP_TOM`, `SWDN_SFC`) print `new vs base (Δ)` for
every case. **`exo_pver=300` is now the standard level count** in
`source/exoplanet_mod.F90`.

### `--exort {h16,h24}` mode

`run_regression.py --exort h16|h24` builds and runs `exort.exe` (swapping the
native-gas HITRAN-16/24 filename strings in `src.exort/kabs.F90`, rebuilding,
then restoring the file) against the **existing n68equiv baselines** — the
"n84 supersedes n68" equivalence check. `_SPECTRAL` arrays are skipped (84-band
vs 68-band rows can't align element-wise).

### `gas_sweep.py` — per-gas HITRAN-2016 vs 2024 sweep

`tests/regression/gas_sweep.py` runs single-gas fixtures (one gas non-zero, N₂
background; realistic + elevated per gas, 12 fixtures) through **n68equiv
(h16/68-grid)**, **n84equiv (h16/84-grid)**, and **exort (h24/84-grid)**, and
tabulates OLR so two effects separate: **grid effect** (n84−n68) and
**line-list effect** (exort−n84). This is the tool that caught the unvalidated
HITRAN-2024 tables (see the data warning in Build Commands). Rerun after a
k-coefficient re-fit:

```bash
python tests/regression/gas_sweep.py                 # all gases
python tests/regression/gas_sweep.py --gases CO2 C2H6
```

It requires `n68equiv`/`n84equiv` `kabs.F90` on the flat `data/kdist/<gas>/`
layout (HITRAN-2016) and `exort` at HITRAN-2024 — the current working state.
`gas_sweep.py` is force-tracked past the `tests/regression/*` gitignore rule.

## Session Handoff (2026-07-06 — 3-D port: src.cam.exort)

**Branch:** `refactor`. Created **`3dmodels/src.cam.exort`**, the v2 CAM
bundle superseding `src.cam.n68equiv` + `n84equiv` + `n68equiv.haze` (all
preserved; `src.cam7.n68equiv` never touched). Pure addition + docs;
regression 15/15 Δ=0.

- **Contract:** `exo_radiation_tend`'s argument list is FROZEN (the
  ExoRT↔ExoCAM linkage) — verified identical to both old bundles.
  Internals use the v2 keyword-tail `aerad_driver` calls; clear-sky calls
  omit condensed-phase args (no `*_zero` fabrication — argument-handling
  Increment 2 is done for this bundle).
- **Optional physics via `-cppdefs`** (one bundle, all ExoCAM configs):
  `-DEXORT_CO2CLD` (CO2 ice clouds; needs the co2condense config's pbuf
  fields + `do_exo_condense_co2`), `-DEXORT_CARMA` (haze; needs
  `carma_exort_mod`, only `carma_exort_get_mmr` is used — ExoRT loads its
  own haze optics via new PIO `initialize_hazeopts`), `-DEXORT_SRF_EMISS`
  (`cam_in%srf_emiss`). No macros = mainline, H2O clouds only.
- **ExoCAM side must add** `do_exo_clouds`/`do_exo_haze` to
  `src.share/exoplanet_mod.F90` (shared driver gates kernels on them) and
  point `exo_solar_file` at `data/stellar/*_n84.nc`. Template:
  `tests/cam_compile_check/exoplanet_mod_stub.F90`; full requirements in
  the bundle README.
- **Sync model:** 13 of 16 files are byte-copies of `source/`;
  `populate3Dmodels.py check --exort` PASSES and is the gate.
  `initialize_rad_mod_cam.F90` (+CO2 cloud optics, +rel/rei/rei_co2 grid
  reads, +hazeopts) and `exo_radiation_cam_intr.F90` (merge) were edited
  in `source/` first, then copied. Old bundles = expected drift, never
  regenerate (documented in the tool).
- **Compile gate without HPC:** `tests/cam_compile_check/run_all.sh` —
  gfortran-compiles the whole bundle against stubbed CESM modules in all
  four macro combinations. All OK at commit.
- **NOT done: real ExoCAM build/run** (no HPC from this machine). First
  validation: clear-sky aquaplanet vs identical `src.cam.n68equiv` case
  (expect LW bit-close, SW <0.4% = n68→n84 stellar UV regrid), then the
  gated physics paths. Known data caveat: `haze_n84_b40_*` optics remain
  provisional above band 68.

## Session Handoff (2026-07-03 — Stage E2 complete)

**Branch:** `refactor`, commits `c176bde` (per-column MCICA seed) + `b47c07c`
(OpenMP + heap arrays) + docs. Every gate green at both commits: regression
15/15 Δ=0, `multicol_check.py` batch==singles exact AND 8-thread==serial
exact, `verify_lib.py` + `tests/lib` PASS. **Stage E is complete (E1+E2);
next is either the HITRAN-2024 k-table re-fit follow-up, the 3-D port
session, or distant E3 (GPU) — see "What is NOT done" in REFACTOR_LOG.md.**

- **OpenMP** (`b47c07c`): `!$omp parallel do schedule(dynamic) if(n>1)` over
  the column loops in `main.F90` and `exort_run_columns`; Makefile appends
  `-fopenmp`/`-qopenmp` to every target (opt out: `make OMPFLAGS=`).
  8-column scaling smoke: 4.4× on 8 threads, bitwise identical to serial.
- **Heap work arrays — the load-bearing discovery:** `-fopenmp` implies
  `-frecursive`, which moved the solver's formerly-STATIC large locals onto
  the stack (~120 MB/solve: aerad_driver 71 MB, refine_lwflux 34 MB, …) →
  instant SIGSEGV/SIGILL on macOS even single-threaded (also explains why
  the code ever ran on an 8 MB main stack: those "automatics" were static).
  All 143 solve-path local arrays ≥200 KB (`exo_radiation_mod`, `mcica`,
  `calc_opd_mod` ×3 bundles) are now `allocatable` + `allocate` at the
  `! Start Code` marker. Default thread stacks now suffice everywhere —
  workers, Python/C library caller threads — **no OMP_STACKSIZE needed**
  (the audit's stacksize plan became unnecessary).
- **Per-column MCICA seed** (`c176bde`, opt-in): namelist flag
  `mcica_percol_seed` (default `.false.` = constant 9404, bit-for-bit; the
  regression gate stays valid). Seed = `get_nstep() + icol - 1` computed in
  `aerad_driver` from optional keyword-tail `ext_icol`; `calc_opd_cld_h2o`
  takes `permuteseed` as an argument in all 3 bundles; flag gating lives in
  `exort_column_run` (1-D-only) so shared driver code never references the
  flag — CAM path untouched. Library: `exort_set_percol_seed(enable)`.
  Functional proof: 3 identical cloudy columns — flag on → col 1 == legacy
  exactly, cols 2/3 decorrelated (max|Δ| ~1e2 W/m² level shifts in slab).
- **Gap found+fixed en route:** the 1-D deck could never drive H₂O cloud
  fraction — `io_1D` didn't read `cfrc` (struct field silently 0 → MCICA
  produced zero cloud from nonzero condensate; `makeColumn.py` writes cfrc
  and documents it). `io_1D` now reads optional `cfrc`; `stackColumns.py`
  stacks it. Also fixed: `calc_opd_cld_h2o` call site aliased `cICE_mcica`
  into both mcica condensate outputs and never passed `cLIQ_mcica`
  (write-only diagnostics; Δ=0).
- **Scope limit CLOSED by E2b** (`934252c`, same day): per-column gravity
  and stellar constant via trailing `grav`/`scon` fields on
  `column_state_t` (append-only ABI; ≤0 = process value) and optional
  `grav`/`scon` input variables. Gravity flows as `ext_grav` on the
  driver keyword tail (the coldens/dzc/masspath block was the only
  solve-path G read — no kernel changes); insolation rides the existing
  `ext_msdist` flux factor (`msdist = scon_process/scon_col`), so a
  per-column scon matches a namelist-configured run to ~1e-15 rel
  (rounding order; gated 1e-10 by the new committed
  `percol_config_check.py` — per-column gravity is gated EXACT, which
  doubles as the detector for any future solve-path `SHR_CONST_G` read).
  Per-column stellar *spectrum* deliberately out of scope
  (maintainer-agreed): mixed-star sets run as an outer loop of
  processes.

## Session Handoff (2026-07-02, session 2 — Stage E through E1)

**Branch:** `refactor`, commits `d94f81d` → `675803a` (5 commits, all
regression-gated 15/15 Δ=0). Stage E is done through E1; **E2 (OpenMP) is the
next task** and everything is at a clean commit boundary.

- **Audit checkpoint** (`d94f81d`): `STAGE_E_AUDIT.md` — module-scope mutable
  inventory of the column solve. Kernel read-only after init except two
  writes, both since fixed. MCICA RNG is call-local (thread-safe) but every
  column gets the constant seed 9404 (`time_manager` stub) — the E2
  per-column seed is a decorrelation feature, not a race fix.
- **Fix 2** (`f5b48ea`): `do camtop=` → `do k=` copy-paste bug in the
  dead-in-1-D `part_in_tshadow` branch (`exo_radiation_mod.F90`); the 3-D
  copies still carry it and get the same line at the 3-D port.
- **Fix 1** (`fb17e2e`): per-column `mwdry`/`cpdry` now flow through
  `aerad_driver`'s optional keyword tail (`ext_mwdry`/`ext_cpdry`; absent →
  physconst module values = the CAM path). `physconst_setgas` **deleted**;
  `calc_opd_gas` takes `mwdry` as an argument in **all three bundles**
  (src.exort + legacy n68/n84equiv, which link the shared driver — 3dmodels
  copies deliberately NOT touched, drift is tracked). After this, a column
  solve writes no module-scope state.
- **E1** (`675803a`): `io_1D.F90` holds no column state — `input_profile`
  fills a `column_state_t` array, `output_data` writes from state/result
  arrays (K/s→K/day moved inside). `run_one_column` lives in ONE shared
  place, new `source/src.main/exort_column_run.F90`, used by both `main.F90`
  (serial column loop) and `exort_lib_mod`. `RTprofile_in.nc` optionally
  carries an `ncol` dimension (absent = bit-for-bit legacy single-column
  layout; present = trailing column dim on every variable, mirrored in the
  output). `tools/stackColumns.py` stacks single-column files;
  `tests/regression/multicol_check.py` proves 3-col batch == 3 singles
  exactly (max|Δ|=0 every variable). Both lib harnesses PASS unchanged.
- **Parallelization decision** (`0b42ea3`, in REFACTOR_PLAN.md Stage E):
  CPU/OpenMP now; GPU deferred as a possible distant Stage E3, gated on the
  (recorded) end goal of `exort_run_columns` with thousands of columns for
  parameter sweeps / emulator training-data generation on NVIDIA hardware.

**Next: E2** — (1) `!$omp parallel do` over the `main.F90` column loop
(and/or `exort_run_columns`); Makefile needs an OpenMP flag path
(`-fopenmp` for gfortran). (2) Opt-in per-column MCICA seed — default MUST
stay constant-9404 (today's behavior; regression Δ=0 gate), enabling it
only changes cloudy/hazy cases and is a rebaseline decision. (3) Document
and set `OMP_STACKSIZE` — the solver carries tens-of-MB automatic arrays
per thread (audit note); worker threads don't inherit the shell ulimit.
Verify: `OMP_NUM_THREADS=1` vs `8` within 1e-10 (use `multicol_check.py`
machinery); regression 15/15 Δ=0 with seed off. **Known scope limit:**
batch columns share the process-level namelist config (star/scon/`exo_g`) —
per-column gravity/insolation is an open API question for the emulator
end goal.

## Session Handoff (2026-07-02)

**Branch:** `refactor`. **Completed: Stage D — `iso_c_binding` library +
Python binding** (pure addition; regression 15/15 Δ=0 before and after).

- New: `source/src.main/exort_column_mod.F90` (bind(c) `column_state_t` /
  `column_result_t`, all-double, field order = ABI contract),
  `source/src.main/exort_lib_mod.F90` (C API: `exort_get_dims`, `exort_init`
  with data_root/solar_file/scon/g/cloud/haze args, `exort_run_column`,
  `exort_run_columns` serial batch, `exort_finalize`), `make libexort`
  target → `run/libexort.dylib|.so`, `tools/exort_pytools/` (cffi binding
  `exort_api.py`, harness `verify_lib.py`), `tests/lib/` (dimension-agnostic
  C test via dlopen + netCDF-C).
- `exort_rootdir` demoted parameter→variable (`src.main/sys_rootdir.F90`
  only; 3-D copies untouched) so `exort_init(data_root)` works at runtime.
- Verified vs committed goldens (TS300K_G2V, TS250K_BB3400, Mars): max rel
  ~6e-8 = float32 storage precision of the baselines (they're nf_real);
  repeat/batch determinism exactly 0. One init per process; `run_column`
  not thread-safe yet (Stage E adds per-column MCICA seed + OpenMP).
- Gotcha: baselines were generated with `_n84` stellar files — library
  callers must pass `G2V_SUN_n84.nc` etc., not the `_n68` names the
  harness case table lists (it maps them).

**Also completed 2026-07-02 — argument-handling Increment 1:** `aerad_driver`
now has a mandatory positional core (thermo/geometry/albedos/all 10 gas MMRs;
gas off = zero MMR) plus an **optional keyword tail after `sol_toa`** (H2O
clouds, CO2 clouds, `carmammr`, `srf_emiss`) — keyword-passed, append-only, so
positional call sites never break when physics args are added. Kernels keep
all-mandatory args; gating stays at the driver call sites (flag on + arg
absent → kernel skipped, zero contribution). `main.F90`/`exort_lib_mod` pass
the tail keyworded unconditionally; the stale source-side
`exo_radiation_cam_intr.F90` was repaired to the new form (clear-sky call
omits cloud args — `*_zero` fabrication deleted). Contract documented in the
`aerad_driver` header and `ARGUMENT_HANDLING.md` (implementation note added).
All 3 targets build; regression 15/15 Δ=0; lib harnesses PASS. The 3dmodels
copies adopt this signature during the 3-D port (they were already stale).

**Next stage: Stage E** (multi-column batch + OpenMP), starting with a
**thread-safety audit as its own committed checkpoint**: enumerate every
module-scope mutable written during a column solve (`exo_radiation_mod`
internals, `calc_opd_mod`, `rad_interp_mod`, MCICA RNG, `physconst_setgas`)
before touching any code; record it in a `STAGE_E_AUDIT.md` or a
`REFACTOR_LOG.md` entry. Then E1 (serial multi-column I/O, optional `ncol`
dim in `io_1D`, single-column path bit-for-bit) and E2 (OpenMP + opt-in
per-column MCICA seed) as separate commits.

**3-D constraints (user decisions, 2026-07-02):** preserve the existing
`3dmodels/` bundles as-is — they may persist as legacy connections to
existing ExoCAM setups; do NOT regenerate during 1-D work (drift is tracked,
currently 11 files, expected). `src.cam7.n68equiv` disposition is an open
question for the maintainer. Argument-handling Increment 2 (keyword-tail
call sites in cam bundles, haze `carmammr_zero` removal) waits for a
dedicated 3-D port session, ideally gated by a real ExoCAM compile. Any
Stage E edit to shared `src.main` files (esp. `mcica.F90`) must keep
CAM-path behavior unchanged by default.

## Session Handoff (2026-07-01)

**Branch:** `refactor` (local commits `cc892f1` → …, on top of `510500d`).

**Completed this session — Stage C3, CARMA haze into `src.exort` (1-D):**

- **C3 code** (`cc892f1`): filled the `calc_aeropd()` stub; new runtime
  namelist flag `do_exo_haze` (default `.false.`); optics loader
  `initialize_hazeopts` (Kext cm²/g → m²/kg, dims checked against
  `nelem_carma=1`/`nbin_carma=40`/`ntot_wavlnrng`); band-indexed kernel
  presums τ, w·τ, g·w·τ over CARMA elements/bins so the `rad_precalc` merge is
  a plain accumulation (same algebra as the published haze bundle); optional
  `carmammr(pver,nelem,nbin)` input via the `nf_inq_varid` pattern. The
  CARMA-module coupling is severed on the 1-D side — `carmammr` comes from the
  deck. **Deliberate deviation:** the published 3-D kernel pairs
  `qcarmammr(ik)` with `pdel(ik)` and never filled the bottom rad level
  (off-by-one); the 1-D kernel instead matches the driver's `coldens`
  convention (rad layer k ← atmos layer k−1). The 3-D kernel in
  `3dmodels/src.cam.n68equiv.haze` was fixed to the same convention on
  2026-07-01 (maintainer-confirmed).
- **C3 data** (`df2b7bb`): `tools/regrid_haze_optics.py` +
  **provisional** `data/aerosol/haze_n84_b40_{mie,fractal_interp}.nc` — bands
  1–68 verbatim from the validated n68 tables (n68 edges are a strict prefix
  of n84, verified), UV bands 69–84 nearest-band-extended (stamped
  `provisional`). Maintainer regenerates properly from
  `data/aerosol/refractive_indices/`, then rebaseline `TS300K_haze_G2V`.
- **C3 test** (`5489019`): `makeColumn.py` haze writer (`add_haze_layer`);
  hazy TS300K fixture (fractal haze slab, layers 40–100, bins 15/18, visible
  τ≈0.5) + generator; `TS300K_haze_G2V` regression case gated via the
  per-case `"haze": True` toggle. **Suite is now 15/15** (14 haze-off Δ=0 +
  1 hazy gated, OLR 266.550). Hazy physics verified sensible & monotonic:
  SWDN_SFC 233.9→194.3 (→131.7 at 5×), OLR 268.7→266.6 (→258.3), SWUP_TOM
  80.1→77.1 (anti-greenhouse + tholin darkening).

- **Legacy build fix** (`02ba2ed`): n68equiv/n84equiv had failed to link since
  C1b (shared `exo_radiation_mod` call sites without legacy symbols). Restored
  with zero-opacity `calc_cldopd_co2`/`calc_aeropd` stubs + loud-stop
  `initialize_{cld,haze}opts` stubs; all legacy targets build, gas_sweep CO₂
  slice reproduces. **Also observed:** `gas_sweep.py`'s exort column no longer
  reflects h24 (the `d77edb9` kabs h16-pinning makes its swap a no-op; linelist
  Δ ≡ 0). Revisit with the h24 re-fit.

**Note:** the `exo_g` compiled default in `source/exoplanet_mod.F90` was
normalized back to `9.80616` (Earth) by the maintainer on 2026-07-01 (it had
been left at Mars `3.711` from cloud testing). `shr_const_mod`'s own default
is still `0.93×9.80616`; both are runtime-overridden by the namelist, which
the harness always writes.

**Also this session:** kernel renames to `calc_opd_gas`/`calc_opd_cld_h2o`/
`calc_opd_cld_co2`/`calc_opd_aero` + uniform `ik=1,pverp` loops over
driver-premapped rad-grid arrays (`masspath` alongside `coldens`); 3-D haze
kernel off-by-one fixed; **`source/experimental/` deleted** (Stage C consumed
it; `n68equiv_exp` target removed). The `_CLD` double-run reference,
experimental `getColumn.pro`, and `_CLD`-aware `plotspectra_1D.pro` are
recoverable at `0e409c3` — concepts recorded in `REFACTOR_LOG.md`'s deletion
entry.

**Still open (decision-gated):** (1) HITRAN-2024 k-table re-fit (CO₂/H₂O);
(2) proper 84-band haze optics regen (provisional tables in place; UV
underestimated); (3) haze/cloud 3-D port + `src.cam.n68equiv.haze`
reconciliation; (4) clear-sky/cloud-forcing `_CLD` double-run (deferred;
reference at `0e409c3`); (5) merge to `main` (maintainer not ready). See
`REFACTOR_LOG.md`.

## Session Handoff (2026-06-28)

**Branch:** `refactor` (pushed to origin, at `76386e5`). `main` and `v1.0.0`
untouched. See `REFACTOR_LOG.md` for the human-readable logbook with per-commit
undo instructions.

**Completed this session (commits `d77edb9` → `76386e5`):**

- **Regression rebaseline to src.exort/h16** (`d77edb9`): the 13-case suite now
  builds+runs `exort.exe` on **HITRAN-2016** native gases as the primary target
  (was `n68equiv.exe`); regenerated golden baselines; `--exort` inverted so `h24`
  is the swap-based side-path. `EXORT_H16_EQUIVALENCE.md` records the one-time
  proof that exort-on-h16 reproduces the retired n68 physics. This decoupled the
  refactor from the unvalidated HITRAN-2024 tables.
- **gas_sweep tooling** (`cad1643`): `--save-spectra`, `--with-n84h24`
  code-isolation column, `tools/plot_gas_sweep_spectra.py`.
- **Stage C — CO₂ ice clouds + surface emissivity** (`03f3007` → `63e27fd`),
  all additive and gated OFF by default:
  - **C1a** restored the 1-D cloud-optics loader (was commented out → 1-D clouds
    never worked) behind a new `do_exo_clouds` flag.
  - **C1b** added the CO₂-ice cloud kernel `calc_cldopd_co2` + solver wiring
    (band-indexed, cloud fraction = 1, no MCICA; condensate is a grid-box mean
    supplied from outside ExoRT). The CO₂-cloud→solver merge was never written in
    the experimental thread — this authored it new.
  - **C1-emiss** added optional `srf_emiss` input (was hardcoded ε=1.0).
  - **C1c** ported the cloud writers from experimental `makeColumn.pro` into
    `tools/makeColumn.py`; committed a CO₂-cloud test fixture + generator.
  - **C1d** made `do_exo_clouds` a **runtime** namelist flag (was compile-time)
    so one build iterates clear/cloudy decks; added the `2barCO2_co2cloud_Mars_G2V`
    regression case. **Suite is now 14/14** (13 clear Δ=0 + 1 cloudy gated).
- **REFACTOR_LOG.md** (`76386e5`): repo-root logbook of the whole v2 refactor.

**Cloud design invariants (Stage C):** clouds do nothing unless BOTH
`do_exo_clouds=.true.` (namelist) AND the input file carries cloud condensate
(`cicewp`/`cliqwp`/`rei`/`rel` for H₂O; `cicewp_co2`/`rei_co2` for CO₂ ice).
Existing decks have neither, so clear-sky physics is bit-for-bit unchanged.
`smart_2bar_t250` is a 69-level profile; the committed Mars fixtures are 300-level
(regridded), so the cloudy fixture is built by *cloning* the 300-level file.

**Still open (decision-gated):** (1) HITRAN-2024 k-table re-fit (CO₂/H₂O); (2)
**Stage C3 — CARMA haze**: `calc_aeropd()` is still an empty stub, blocked on
regenerating 84-band haze optics (`data/aerosol/haze_n84_*.nc`; only
`haze_n68_b40_*` exist); haze enters via a `carmammr` input array, no CARMA-module
coupling on the 1-D side; (3) clear-sky / cloud-forcing `_CLD` double-run
(deferred); (4) merge to `main` (maintainer not ready). See `REFACTOR_LOG.md`.

## Session Handoff (2026-06-17)

**Branch:** `refactor`. v2 single-bundle work landed (commits `77fd3ee` →
`45c1897`, local — not pushed). Builds green: `exort`, `n68equiv`, `n84equiv`,
`n68equiv_exp`.

**Completed (v2 Stages A–B + HITRAN-2024 validation):**
- **Stage A — prune legacy bundles** (`77fd3ee` ≈ `fa013cf`): deleted
  `n28archean`/`n42h2o`/`n68h2o` from `source/`, `3dmodels/`, `build/Makefile`,
  and `populate3Dmodels.py`. Preserved in the `v1.0.0` tag.
- **Stage B — merge to `src.exort`** (`ff7a84e`): n84 grid + NH₃/CO (nspecies 8);
  new `make exort`. n84 supersedes n68 (runtime band optimizer). Verified the
  merge is physics-neutral — LW bit-equivalent to n68 baselines; SW <0.4% is the
  n68→n84 stellar UV regrid, not the k-merge.
- **Data reorg** (`93d9f00`, `51f4b9c`): `data/solar/`→`data/stellar/`;
  `data/kdist/` flattened to per-gas dirs; HITRAN-2024 84-band k-files added;
  `hitran2024`→`hitran24` filename normalization; `blackbody_3400K_n84` added.
- **Wiring** (`a1d7114`): `src.exort/kabs.F90` + `exoplanet_mod.F90` point at the
  new layout; default `solar_file = G2V_SUN_n84.nc`.
- **Test tooling**: `run_regression.py --exort {h16,h24}` (`ed6e80c`);
  `gas_sweep.py` + 12 single-gas fixtures (`a69524e`); NH₃/CO added to n84equiv
  so it's a valid sweep reference (`d04bf7f`).

**⚠️ OPEN ISSUE — HITRAN-2024 k-coefficients (HELIOS-K pipeline):**
The gas sweep found the 84-band HITRAN-2024 tables give non-physical LW for
**H₂O, CO₂, and C₂H₆**. CH₄/NH₃/CO are clean. **The ExoRT code is proven
correct** — a control run (exort pointed at the n84 *HITRAN-2016* files)
reproduces n68/n84 bit-for-bit on all 12 sweep cases, so the defect is entirely
in the k-files. After any re-fit: drop new files into `data/kdist/<gas>/`, rerun
`USER_FC=gfortran python tests/regression/gas_sweep.py` (default `ifort` won't
link on arm64), expect line-list Δ → ~0. Full diagnosis in memory:
`h24-gas-sweep-findings.md`, `h24-co2-farwing-artifact.md`.

Per-gas status (updated 2026-06-17):
- **C₂H₆ — FIXED & committed (`1a536b7`).** Root cause: the original h24 C₂H₆
  table had been generated with the **CH₄ line list** (byte-identical to the CH₄
  h24 file). Re-fit with the correct C₂H₆ list: k-table median h24/h16 ratio
  0.241 → 1.000; gas sweep `C2H6_elevated` line-list Δ +27.0 → −0.039. Done.
  (Will be regenerated again in the full-set rerun, but currently valid.)
- **CO₂ — still broken, fully fingerprinted.** Far-IR rotation-band wing
  (bands 3–8, 24–81 µm) has k ∝ P¹ in h24 vs P⁰ in h16; band-5 k24/k16 grows
  6× → 76,000× from 0.01 mb → 10 bar. = sub-Lorentzian χ-factor not applied to
  the wing *below* band center (< ~440 cm⁻¹). 2-bar CO₂ loses ~48% OLR. A
  partition-file-cleanup rerun came back **byte-identical** → it's a recipe bug,
  invariant to line-data. Fix is the χ-factor application in HELIOS-K.
- **H₂O — possibly partly real; needs a literature check.** Not the "uniform
  +12%" first thought: per-band ratio is spectrally structured (8–25 µm window
  +1.3–1.75×, 2.7/4.3 µm vib +1.9×, a few bands <1) and **pressure-flat** (slopes
  match h16) → a per-band intensity/line-shape offset, *not* continuum or plinth.
  HITRAN 2016→2024 did revise H₂O intensities, so this may be adoptable — verify
  the per-band shape against published 2024 water revisions. Band 63 (0.37 µm,
  112× too strong, sign-flipped P-slope) is a discrete bad UV band = bug
  regardless. Partition rerun also byte-identical.

**Working-tree note:** `n68equiv`/`n84equiv` `kabs.F90` are hand-edited to the
flat `data/kdist/<gas>/` layout (test scaffolding so the legacy codes run against
v2 data) — committed as such. The regression harness `run_regression.py` (13
n68 cases) is currently NOT rebaselined to `exort` (blocked on valid h24 files).

**Next:** (1) maintainer re-fits HELIOS-K k-files; (2) rebaseline `run_regression`
to `src.exort` once h24 is valid; (3) Stage C — fold experimental haze + CO₂
clouds into `src.exort`; (4) retire n68equiv/n84equiv. See `REFACTOR_PLAN.md`.

---

## Session Handoff (2026-06-10)

**Branch:** `refactor` (pushed to origin). Builds green for all 6 targets;
n68equiv regression 13/13 bit-for-bit.

**Completed this session (commits edf9287 → 79c3278):**
- **Regression harness** (`edf9287`, `089dfad`): `tests/regression/run_regression.py`
  + 13 committed fixtures/baselines (6 TS profiles × 2 stars + Mars 2barCO2).
- **regrid gravity bug** (`a97db99`): `tools/regrid_rtprofile.py` hardcoded Earth
  gravity when rebuilding `zint`; `zint` sets the CIA path length in 1-D RT (via
  aerad_driver → zlayer → calc_gasopd `pathlength`), so a wrong g biased OLR.
  `--gravity` is now required. Set `exo_pver=300` as the standard.
- **main.F90 hygiene**: per-phase CPU timing (`8500b4d`); inlined the `h2oint`
  scratch array (`0fd9716`).
- **physconst import fix** (`2107866`): it referenced `shr_const_cpc2h6` but
  borrowed it via a vestigial `use input` (whole-module re-export). Fixed the
  `only:` list, dropped `use input` — severed a phantom physconst→input edge.
- **io_1D consolidation** (`b553be3`): merged `input.F90` + `output.F90` →
  `source/src.main/io_1D.F90` (module `io`); extracted `read_namelist` and
  `print_diagnostics` out of `main.F90`. Unblocked by the physconst fix.
- **calc_gasopd signature fix** (`e0fd01e`): n84/n42/n68h2o/n28archean failed to
  build because the shared `exo_radiation_mod.F90` passes `qnh3`/`qco` to
  `calc_gasopd`, but only n68equiv's `calc_opd_mod.F90` had been updated. Added
  unused dummy args to the other four — all 5 production builds now compile.
- **src.misc cleanup** (`4c12119`, `fed637b`, `6e5093a`, `2729a5a`, `79c3278`):
  removed ~3,775 lines of dead CESM shim. Trimmed `wrap_nf.F90` to the 5 used
  wrappers; deleted 4 unreferenced headers; trimmed `ioFileMod` to `getfil` and
  deleted now-orphan `shr_sys_mod`; deleted the stale orphan `src.misc/mcica.F90`
  (the live one is `src.main/mcica.F90`); **resolved the NetCDF version skew** —
  the build now uses the library's own `netcdf.inc` via
  `-I$(nf-config --includedir)` instead of a stale NetCDF-3 local copy. This
  fixed the v4-upgrade compatibility class of issues *and* repaired the
  previously-broken `n68equiv_exp` build.

**src.misc rule (important):** `source/src.misc/` is a local ExoRT-only shim
("not used in 3dmodels RT builds, as it is already present in CESM" — its README).
It is safe to modify/condense/delete. The files duplicated in `/3dmodels/`
(byte-identical to `source/src.main` and `source/src.n68equiv`) are
**dual-contract with CESM** — their `use`/`include` statements must NOT change.
1-D-local files (e.g. experimental `input.F90`) may be changed freely.

**3dmodels/ status:** Untouched. `src.cam7.n68equiv` and `src.cam.n68equiv.haze`
are special cases whose files do NOT map cleanly to `source/` and must be
preserved.

**Known broken (pre-existing, not from this session):** legacy 3-D-only paths and
the n42h2o/n68h2o/n28archean *physics* (they build now but only n84equiv is
physics-ready to match the regression suite — it needs an `_n84.nc` solar file +
baselines). Disposition of the legacy versions is undecided.

**Next:** see `REFACTOR_PLAN.md`. Remaining `src.misc` files are all live or
vendor headers — the high-confidence dead weight is gone. Larger open work:
source↔3dmodels dedup (Stages 1–2, gated on CESM-ExoCAM workflow decisions) and
unifying experimental haze / CO₂-ice optics.
