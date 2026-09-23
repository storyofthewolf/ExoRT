# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What This Is

ExoRT is a two-stream radiative transfer code for exoplanet climate modeling. It supports both 1-D offline calculations and 3-D integration with CESM/CAM4 (via ExoCAM). The codebase is primarily Fortran 90, with Python and IDL tools for pre/post-processing.

**Branches.** `main` = v1 maintenance (`v1.0.0` published checkpoint, `v1.1.0` = H₂O
T-index fix). `refactor` = v2 development — a single branch, no sub-branches. What v2
changed and why: `CHANGES.md`. What is left: `REFACTOR_PLAN.md`.

## Where to look

| Need | File |
|---|---|
| What v2 changed vs v1, with the measured physics deltas | `CHANGES.md` |
| Open work, priorities, ground rules for changes | `REFACTOR_PLAN.md` |
| Implementation reference (library API, driver contract, k-file format, known limitations) | `DEVELOPER_NOTES.md` |
| Current test record | `tests/regression/REGRESSION_STATUS.md` |
| How-tos | `HAZE_1D_HOWTO.md`, `MULTICOLUMN_BATCH_HOWTO.md`, `tools/exort_pytools/README.md`, `3dmodels/src.cam.exort/README` |
| Performance analysis (not yet applied) | `OPTIMIZATION.md` |
| CESM3 / `src.cam7.n68equiv` assessment | `CESM.CAM7.md` |

## Build Commands

All builds run from `build/`:

```bash
cd build

make exort           # v2 single bundle (84-band, HITRAN-2024, NH3/CO) — primary target
make libexort        # shared library (run/libexort.dylib|.so) with the C API
make n68equiv        # legacy HITRAN-2016 reference (slated for retirement)
make n84equiv        # legacy HITRAN-2016 reference, +UV bins (slated for retirement)

make clean           # remove all build artifacts, run/*.exe, run/libexort.*
```

- **Compiler.** Defaults to `ifort`. On Apple Silicon use
  `USER_FC=gfortran make exort`; ifort has no arm64 port. Requires NetCDF4
  Fortran (`nf-config` on PATH). Executables go to `run/`. All targets build
  with OpenMP; opt out with `make OMPFLAGS= <target>`.
- **Link gotcha (this Mac).** `make exort` fails with
  `library not found for -lnetcdf` unless `LIBRARY_PATH=/opt/homebrew/lib` is
  set, because `nf-config --flibs` omits the netCDF-C directory. The Makefile
  is deliberately unchanged.
- **Library API.** `libexort` exposes `exort_init` / `exort_run_column` /
  `exort_run_columns` / `exort_finalize` / `exort_get_dims` /
  `exort_set_percol_seed` (`source/src.main/exort_lib_mod.F90`; structs in
  `exort_column_mod.F90`). The Python binding is in `tools/exort_pytools/`.
  Acceptance checks: `python tools/exort_pytools/verify_lib.py` and
  `cd tests/lib && make run`. To reproduce the committed baselines through the
  library, pass the `_n84` stellar files; the regression harness maps
  `_n68 → _n84` for exort runs.
- **Legacy bundles.** `n28archean` / `n42h2o` / `n68h2o` were removed in v2
  (they live in `v1.0.0`). `n68equiv`/`n84equiv` remain only as HITRAN-2016
  references for `gas_sweep.py`. Their `kabs.F90` paths are hand-edited to the
  flat `data/kdist/<gas>/` layout.

**Line list.** `src.exort` reads HITRAN-2024 for H₂O/CO₂/CH₄/C₂H₆/NH₃/CO and
HITRAN-2020 for O₂/O₃. The HITRAN-2016 native-gas tables stay in
`data/kdist/<gas>/`, reachable via `run_regression.py --exort h16`.

**HITRAN-2016 H₂O — the 25 K temperature-index error.**
- **What.** The table n68equiv used from 2020-11-05 (`d3d340e`) to 2026-09-23
  held k(T−25 K) in the slot labelled T, about 20 % under-absorbing. Wolf et al.
  2022 Fig. 1 used it.
- **Where it lives now.** It is kept as
  `data/kdist/h2o/n{68,84}_…_grrtm_Tindex-error.nc`; the corrected table has
  the plain `…_grrtm.nc` name. The k-file grid check rejects the broken file.
- **Proof.** The 2020 per-bin files' own Temperature coordinate reads
  100, 100, 125 … 475; sibling HELIOS-K families match the fix.
- **Other branches.** `main` is fixed as `v1.1.0` (`3dd35ba`). The frozen
  `3dmodels/src.cam.n68equiv*`/`.n84equiv` bundles point at the pre-v2
  `data/kdist/n68h2o/hitran2016/`, which does not exist on `refactor`.
- **Deltas:** `CHANGES.md`.

## Running the 1-D Model

```bash
# 1. exo_pver in source/exoplanet_mod.F90 must match the input's level count (compile-time; 300)
# 2. Build: cd build && USER_FC=gfortran make exort
# 3. Place input file in run/ as RTprofile_in.nc
# 4. Optionally create run/user_nl_exort (runtime config, below)
cd run
./exort.exe
# Output written to run/RTprofile_out.nc
```

Template inputs are in `iofiles/input_files/`; the regression fixtures in `tests/regression/fixtures/` are 300-level decks. The input may carry any number of columns via an optional `ncol` dimension (see "Input File Behavior" below); one invocation then solves all columns.

### Runtime namelist: `user_nl_exort`

These formerly compile-time settings are now read at runtime, with no rebuild:

| Variable | Meaning | Default |
|----------|---------|---------|
| `solar_file` | Stellar spectrum filename in `data/stellar/` (must be `*_n84.nc` for exort) | `'G2V_SUN_n84.nc'` |
| `shr_const_scon` | Stellar constant ÷ 2 [W m⁻²] | `680.0` (≈ present Earth) |
| `exo_g` | Surface gravity [m s⁻²] | `9.80616` (Earth) |
| `do_exo_clouds` | Enable the cloud RT path (H₂O + CO₂ ice; reads `cicewp*`/`rei*`/`cfrc` from the input file) | `.false.` |
| `do_exo_haze` | Enable the CARMA haze aerosol RT path (reads `carmammr(pver,nelem,nbin)`; optics from `data/aerosol/haze_n84_b40_mie.nc` — see "Haze optics" below) | `.false.` |
| `mcica_percol_seed` | Opt-in per-column MCICA seed. Each batch column offsets the stochastic-cloud seed by its column index; column 1 stays bit-identical to legacy. Off = constant seed 9404 for every column. Affects cloudy H₂O runs only; enabling it is a rebaseline decision | `.false.` |

```bash
cp iofiles/user_nl_exort.template run/user_nl_exort   # edit, then run
```

If `user_nl_exort` is absent the model runs with the defaults above. The active values are always printed at startup under the `=== exort_config ===` header. Standard Fortran namelist syntax (`&exort_config … /`); only variables you want to change need to appear. `exo_pver` is **not** in the namelist — it is compile-time only.

## Architecture

### Source Layout

```
source/
  exoplanet_mod.F90          # Compiled defaults: solar file, gravity, pver, scon
  src.main/                  # Shared drivers and solver
  src.misc/                  # CESM stubs required for standalone operation
  src.exort/                 # v2 single RT bundle (primary)
  src.n68equiv/              # legacy HITRAN-2016 reference (slated for retirement)
  src.n84equiv/              # legacy HITRAN-2016 reference, +UV bins (slated for retirement)
3dmodels/src.cam.exort/      # v2 CAM hook-up (byte-synced with source/; CO2 clouds/haze/emiss via -cppdefs)
3dmodels/src.cam.*/          # pre-v2 CAM hook-ups (frozen legacy copies)
```

### Key Module Responsibilities

Bundle-specific files (`src.exort/`; the legacy references carry their own copies):
- `radgrid.F90` — spectral intervals, Gauss points, P-T grid (`tgrid`, `pgrid`)
- `kabs.F90` — k-file directories and filenames
- `calc_opd_mod.F90` — optical depth kernels (`calc_opd_gas`, `_cld_h2o`, `_cld_co2`, `_aero`)
- `rad_interp_mod.F90` — interpolation in P-T space
- `cloud.F90` — cloud and haze optics filenames/arrays
- `model_specific.F90` — gas species configuration
- `initialize_rad_mod_1D.F90` / `initialize_rad_mod_cam.F90` — table loading (1-D / CAM), including `check_kfile_grid`

Shared (`src.main/`):
- `exo_radiation_mod.F90` — `aerad_driver` and the two-stream solver (argument contract in its header)
- `exo_init_ref.F90` — reference atmosphere, band optimization
- `main.F90` — 1-D entry point (OpenMP loop over columns)
- `io_1D.F90` — NetCDF I/O (module `io`; holds no column state)
- `exort_column_mod.F90` — bind(c) `column_state_t`/`column_result_t` structs
- `exort_column_run.F90` — `run_one_column`, the single per-column solve used by both `main.F90` and the library
- `exort_lib_mod.F90` — C API
- `exo_radiation_cam_intr.F90` — CAM interface (`exo_radiation_tend`; its argument list is the frozen ExoRT↔ExoCAM contract)
- `mcica.F90` — Monte Carlo Independent Column Approximation for clouds
- `planck_mod.F90` — Planck function

### 3-D Model Integration (ExoCAM)

`3dmodels/` contains file copies (not symlinks) for CAM integration:

```bash
# v2 (see 3dmodels/src.cam.exort/README for ExoCAM-side requirements)
xmlchange CAM_CONFIG_OPTS="-usr_src /$MYDIR/ExoRT/3dmodels/src.cam.exort"
# optional physics: -cppdefs '-DEXORT_CO2CLD -DEXORT_CARMA -DEXORT_SRF_EMISS'

# pre-v2 legacy bundles (frozen; use them from main, not refactor)
xmlchange CAM_CONFIG_OPTS="-usr_src /$MYDIR/ExoRT/3dmodels/src.cam.n68equiv"
```

`src.cam.exort` is kept byte-identical to `source/`. Edit `source/`, then
`python tools/populate3Dmodels.py regenerate --exort`, then `check --exort`,
then the local compile gate `tests/cam_compile_check/run_all.sh`.
`src.cam.exort` has not yet run in a real ExoCAM case. The pre-v2 bundles
intentionally diverge from `source/` and must not be regenerated.

## Data Files

| Directory | Contents |
|-----------|----------|
| `data/stellar/` | Stellar spectra, one per grid (`*_n84.nc` for exort, `*_n68.nc` for n68equiv) |
| `data/kdist/` | Correlated k-distribution tables, **flat per-gas dirs `data/kdist/<gas>/`**. Each gas dir holds its `hitran16…`/`hitran24…` files keyed by filename. |
| `data/cia/` | Collision-induced absorption (N₂–N₂, N₂–H₂, H₂–H₂, CO₂–CO₂, O₂ pairs, …) |
| `data/continuum/` | H₂O and CO₂ continuum from MT_CKD/LBLRTM |
| `data/cloud/` | Mie scattering cloud optical properties (H₂O liquid/ice, CO₂ ice) |
| `data/aerosol/` | CARMA haze optics and refractive indices |

### k-file grid check (`check_kfile_grid`)

Both `src.exort` loaders check every gas k-file at load:
`initialize_rad_mod_1D.F90` (1-D exe + `libexort`, `stop 1`) and
`initialize_rad_mod_cam.F90` (CAM, `endrun`). A file is rejected unless the
`data` dims equal the compiled (bands, gauss, press, temp) sizes and the
coordinates are present, nonzero, and equal to the compiled grid:
- `Temperature` = `tgrid`
- `Pressure` = `pgrid` in **mb** (rel 1e-4)
- `GaussWeights` = g-interval **midpoints**,
  `g_xpos_edge_8gpt + g_weight_8gpt/2`. The variable is misnamed; it never
  held weights.

`SpectralBands` holds only indices 1…N, so bands are checked by count.
**New k-files must carry correct coordinates or they will not load.**

Coordinate repairs made on 2026-09-23 (`data` bit-identical in every file):
- **O₂/O₃/`null`:** Pressure was stored in bar under an `mb` label; rescaled
  ×1000.
- **HITRAN-2016 CO₂/CH₄/C₂H₆:** coordinates had been zeroed by the 2023 merge;
  restored from their per-bin source files.
- **`_Tindex-error` H₂O:** now carry their true labels, 100, 100 … 475.

### Haze optics (`data/aerosol/haze_*`)

`src.exort` reads **`haze_n84_b40_mie.nc`** — Mie spheres, regenerated over
the full 84-band grid from the Khare et al. (1984) tholin indices by
`tools/makeCARMAOptics.py` (2026-09-15). Bands 69–77 rest on real Khare
measurements; 78–84 are extrapolated (Khare ends at 0.0588 µm) but sit in the
EUV where stellar flux is negligible.

`haze_n84_b40_fractal_interp.nc` (fractal aggregates) is **still provisional**
above band 68 — its UV bands are a nearest-band copy of band 68. Regenerating
it requires the external mean-field fractal solver (`fractaloptics.exe` on
Discover), which is why the Mie table is the default. Switching back is a
one-line change in `source/src.exort/cloud.F90` plus a rebaseline.

**Units — important.** `Kext` in all four haze files is **m² kg⁻¹**, and the
loaders apply **no** conversion. `rbins` is in **cm**. (Older attributes
claimed cm² g⁻¹ and microns; they were corrected in place on 2026-09-15, and
a leftover `kcarma*0.1` that made haze τ 10× too small was deleted.)

## Tools

`tools/` contains pre/post-processing scripts. Most are IDL (`.pro`); Python equivalents exist for the primary workflow:

### Python tools (preferred for new work)

- **`profile_data.py`** — P-T profile data library. Contains `smart_2bar_t250`, `TS273K`, and `US1976` profiles as numpy arrays. Use `get_profile(tag)` to retrieve a profile dict with keys `pint`, `pmid`, `tint`, `tmid`, `q`. Equivalent to `profile_data.pro`.
- **`makeColumn.py`** — Build `RTprofile_in.nc` from scratch. Edit the `USER SETTINGS` block at the top, or use CLI flags. All gas species are optional; absent species are simply omitted from the output file and ExoRT reads them as zero. Also writes clouds, haze, `srf_emiss`, and per-column `grav`/`scon`. Equivalent to `makeColumn.pro`.
  ```bash
  python makeColumn.py                        # use USER SETTINGS block
  python makeColumn.py --defaults             # print current defaults
  python makeColumn.py --profile TS273K --co2vmr 400e-6 --nh3vmr 1e-6 --output RTprofile_in.nc
  python makeColumn.py --zero-h2o            # force H2O to zero
  ```
- **`stackColumns.py`** — stack single-column files into one multi-column (`ncol`) input.
- **`plotspectra_1D.py`** — Plot shortwave and longwave spectral fluxes from `RTprofile_out.nc`. SW panel uses wavelength axis (µm, W m⁻² µm⁻¹); LW panel uses wavenumber axis (cm⁻¹, W m⁻² cm). Supports two-file comparison. Equivalent to `plotspectra_1D.pro`.
  ```bash
  python plotspectra_1D.py                    # use USER SETTINGS block
  python plotspectra_1D.py --save             # save PNG instead of display
  python plotspectra_1D.py --f2 RTprofile_out_ref.nc --label2 "no NH3"
  ```
- **`spectral_intervals.py`** — Spectral grid edge arrays and Gauss weights for all RT versions. Call `get_spectral_intervals(tag)` with tag `n28`, `n42`, `n68`, or `n84`; it returns `(nbands, wavenum_edges, …)`.
- **`makeStellarSpectrum_blackbody.py`** — Generate stellar spectrum NetCDF from blackbody temperature.
- **`check_kcoeff.py`** — Inspect k-coefficient files.
- **`mie.py`** — Mie scattering for homogeneous spheres (Bohren & Huffman). Port of `mie_single.pro`; agrees with it to ~4e-6 worst case. `python tools/mie.py` runs a benchmark selftest against published values.
- **`makeCARMAOptics.py`** — Build CARMA aerosol optics tables on any ExoRT grid (`--grid n28|n42|n68|n84`). Port of the **Mie branch** of `makeCARMAOptics.pro`; regenerating the n68 grid reproduces the committed IDL-made table to 1.7e-4 (median 4e-7). The fractal branch is not ported — it shells out to an external solver.
- **`populate3Dmodels.py`** — Sync `source/` → `3dmodels/src.cam.exort` (`check` / `diff` / `regenerate` / `list`).

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

**Multi-column input.** `RTprofile_in.nc` may carry an optional `ncol` dimension.
- **Absent `ncol`:** a classic single-column file, bit-for-bit the legacy path.
- **Layout:** every variable carries a trailing column dimension in Fortran
  order (`tmid(pver,ncol)`, `ts(ncol)`, `carmammr(pver,nelem,nbin,ncol)`), and
  `RTprofile_out.nc` mirrors the `ncol` dimension.
- **Threading:** columns are solved in an OpenMP-parallel loop. Threads spawn
  only when `ncol > 1`; the count comes from `OMP_NUM_THREADS`; results are
  bitwise independent of thread count.
- **Building inputs:** `python tools/stackColumns.py col1.nc col2.nc -o RTprofile_in.nc`.
- **Per-process vs per-column config:** columns share the process-level
  `solar_file` and cloud/haze flags. Optional scalar `grav`/`scon` input
  variables override `exo_g`/`shr_const_scon` for that column (absent or
  ≤0 = namelist value).
- **Cloud fraction:** optional `cfrc` input drives H₂O cloud fraction
  (absent = zero).
- **Checks:** `tests/regression/multicol_check.py` and
  `percol_config_check.py`.
- **Recipe:** `MULTICOLUMN_BATCH_HOWTO.md`.

## .gitignore Notes

Patterns in `.gitignore` are relative to the repo root. The prefix `/ExoRT/` does **not** work — it would only match if the repo were checked out inside a directory literally named `ExoRT`. Use bare relative paths:

```
tools/*.nc       # correct
/ExoRT/tools/*.nc  # wrong — silently ignored by git
```

`tests/regression/*` is ignored except for whitelisted files (`!` lines); a new
committed file there needs a whitelist entry or `git add -f`.

## Adding a New Gas Absorber (`src.exort`)

NH₃ and CO (2026-04-27, then carried into `src.exort`) are the reference example — grep for `nh3` to see every touch point. Touch these files in this order:

1. **`source/src.misc/physconst.F90`** — add `mwXXX` molecular weight parameter
2. **`source/src.exort/kabs.F90`** — add `dirk_xxx` (`'data/kdist/xxx/'`) and `k_xxx_file` (filename) parameters
3. **`source/src.exort/radgrid.F90`** — increment `nspecies`, add `iXXX` index parameter, extend `gas_name`, add `k_xxx` to the individual k-array declaration
4. **`source/src.exort/initialize_rad_mod_1D.F90`** — add a `getfil` + `wrap_open` + `check_kfile_grid` + `wrap_get_var_realx` block (copy the NH₃ block)
5. **`source/src.exort/model_specific.F90`** — add one line each to `k_major_data` and `k_grey_data` assignments
6. **`source/src.main/exort_column_mod.F90` + `io_1D.F90`** — add a `xxxmmr(pver)` field to `column_state_t` and zero it in `zero_column_state`. In `input_profile`, read it via the optional `opt_mid` pattern; add the output variable in `output_data`.
7. **`source/src.main/exo_radiation_mod.F90`** — add the MMR argument to `aerad_driver` (positional core, with the other gases); add a `qXXX` local array; assign pseudo-layer and mid-layer values; pass to `calc_opd_gas`
8. **`source/src.main/exort_column_run.F90`** — add `state%xxxmmr` to the `aerad_driver` call in `run_one_column` (also update the `column_state_t` mirrors in `tools/exort_pytools/exort_api.py` and `tests/lib/`)
9. **`source/src.exort/calc_opd_mod.F90`** — add `mwxxx` to the physconst import; add a `qxxx` argument; add `xxxvmr`/`u_xxx` locals; compute VMR and column density in the level loop; extend `ugas`; add the `bilinear_interpK_grey` call and `tau_grey` entry. Add Rayleigh terms in `src.main/rayleigh_data.F90` if needed.

**Data files:** the k-table goes in `data/kdist/xxx/` on the same grid as the other species (`NTemp=17`, `NPress=61`, `NGauss=8`, `NBins=84`), with variable `data` **and correct `Temperature`/`Pressure`[mb]/`GaussWeights` coordinates** — the grid check rejects it otherwise. Input decks carry the matching `xxxmmr`. Add a single-gas regression case (see `_gas_cases()` in `run_regression.py`).

### Connecting it to the 3-D interface (`src.cam.exort` + ExoCAM)

10. **`source/src.exort/initialize_rad_mod_cam.F90`** — add a PIO read block (`getfil` + `cam_pio_openfile` + `check_kfile_grid` + `pio_inq_varid` + `pio_get_var` + `pio_closefile`) and a matching `mpibcast` in the `#if defined SPMD` block.
11. **`source/src.main/exo_radiation_cam_intr.F90`** — import `mwXXX` from `physconst` and `exo_xxxmmr` from `exoplanet_mod`. Declare a `xxxmmr(pcols,pver)` local, set `xxxmmr(:,:) = exo_xxxmmr`, and pass `xxxmmr(i,:)` at both `aerad_driver` call sites (clear-sky and full-sky).
12. **Sync:** `populate3Dmodels.py regenerate --exort` → `check --exort` → `tests/cam_compile_check/run_all.sh` (add `exo_xxxmmr` to `exoplanet_mod_stub.F90` there).
13. **ExoCAM** (all active configs under `cesm1.2.1/configs/*/SourceMods/`):
    - `src.cam/physconst.F90`: `mwXXX` plus a VMR/MMR `write(iulog,*)` line.
    - `src.share/exoplanet_mod.F90`: `exo_xxxbar = 0.0_r8`, `cpXXX`, the derived `exo_xxxvmr`/`exo_xxxmmr`, and `exo_xxxmmr*cpXXX` in `exo_cpdair`.

## Regression Test Suite

A Python regression harness lives in `tests/regression/` (see its README). It builds and runs `run/exort.exe` over 32 standard cases and compares flux/heating/spectral outputs against committed golden baselines with rtol/atol = 1e-3.

```bash
cd tests/regression
python run_regression.py                 # run all cases, compare to baselines
python run_regression.py --list          # list case names + their physics
python run_regression.py --cases TS300K Mars   # subset by substring
python run_regression.py --generate-baselines  # (re)create golden baselines
python run_regression.py --exort h16     # HITRAN-2016 side-path (comparison, not a gate)
```

- **Cases** are defined in `build_cases()`; each carries its own fixture,
  stellar spectrum, insolation and gravity. The 32:
  - 12 Earth-like TS250K–TS360K × {G2V, blackbody_3400K};
  - clear and CO₂-cloudy 2-bar Mars (g = 3.711);
  - thin and thick haze on TS300K;
  - 16 minor-gas cases, CH₄/CO/NH₃/C₂H₆ × {realistic, elevated} × 2 stars,
    single gas in N₂ on TS273K.

  All are `pver=300`.
- **Environment.** The harness sets the NetCDF lib path for the macOS loader
  and preserves/restores any existing `run/user_nl_exort`.
- **Output.** Two inline glance columns (`LWUP_TOM`, `SWDN_SFC`) print
  `new vs base (Δ)` for every case.
- **When to use it.** A change that should not move physics must show Δ=0.
  An intended change regenerates baselines in the same commit and refreshes
  `REGRESSION_STATUS.md` (rules in `REFACTOR_PLAN.md`).
- **`--exort h16`** temporarily swaps the four native-gas filenames in
  `src.exort/kabs.F90` to HITRAN-2016, rebuilds, runs against the h24
  baselines, then restores the file. It cannot be combined with
  `--generate-baselines`.

**Other gates** (all must pass before committing shared-code changes):
`multicol_check.py`, `percol_config_check.py`, `verify_lib.py`, `tests/lib`
(`make run`), `populate3Dmodels.py check --exort`, and
`tests/cam_compile_check/run_all.sh`.

### `gas_sweep.py` — per-gas grid vs line-list sweep

`tests/regression/gas_sweep.py` runs the single-gas fixtures through three
builds and tabulates OLR:
- **n68equiv** (HITRAN-2016, 68-band grid);
- **n84equiv** (HITRAN-2016, 84-band grid);
- **exort** (HITRAN-2024, 84-band grid).

That separates the **grid effect** (n84−n68) from the **line-list effect**
(exort−n84). Rerun it after any k-coefficient re-fit:

```bash
python tests/regression/gas_sweep.py                 # all gases
python tests/regression/gas_sweep.py --gases CO2 C2H6
```

## Current state (2026-09-23)

`refactor`, all gates green (`REGRESSION_STATUS.md`).
- **Done:** stages A–E, the Stage C physics, the `src.cam.exort` code, the
  HITRAN-2024 default, the H₂O T-index fix, minor-gas regression cases, and
  the k-file grid check in both loaders.
- **Next:** in priority order in `REFACTOR_PLAN.md`. The top item is the first
  real ExoCAM run of `src.cam.exort` (HPC).
- **Working constraints:**
  - Never regenerate the pre-v2 `3dmodels/` bundles.
  - Shared `src.main` edits must keep CAM-path behaviour unchanged by default.
  - Keep steps small and mechanically verifiable, so a later session on a
    smaller model can pick them up.
