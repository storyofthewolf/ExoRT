# DEVELOPER_NOTES.md

Implementation reference for ExoRT (v2 / `refactor` branch). Granular detail for
code review — complements `README.md` (user-facing) and `CLAUDE.md` (architecture
+ session handoff). See `CLAUDE.md` for the "add a gas" 9-file pattern, the
source↔3dmodels sync rule, and the current HITRAN-2024 validation status.

---

## Build targets (`build/Makefile`)

| Target | Produces | Notes |
|--------|----------|-------|
| `make exort` | `run/exort.exe` | v2 single bundle (84-band, HITRAN-2024, `nspecies=8`) — primary |
| `make libexort` | `run/libexort.dylib` (`.so` on Linux) | Stage D shared library: same physics objects as `exort.exe` with `main.F90` replaced by the C-API modules |
| `make n68equiv` | `run/n68equiv.exe` | legacy HITRAN-2016 reference, 68-band |
| `make n84equiv` | `run/n84equiv.exe` | legacy HITRAN-2016 reference, 84-band + UV |
| `make clean` | — | removes build artifacts, `run/*.exe`, `run/libexort.*` |

- Default compiler `ifort`; on Apple Silicon add `USER_FC=gfortran` to every
  `make`. Requires `nf-config` on PATH (NetCDF4 Fortran).
- Legacy `n28archean` / `n42h2o` / `n68h2o` targets were removed in v2 (present
  in the `v1.0.0` tag).

---

## Library API (`libexort`, Stage D)

C entry points in `source/src.main/exort_lib_mod.F90` (bind(c) structs in
`exort_column_mod.F90`; field order is the ABI contract):

| Function | Purpose |
|----------|---------|
| `exort_get_dims(pver, pverp, nwave, nelem, nbin)` | query compiled dims (exo_pver etc.) |
| `exort_init(data_root, solar_file, scon, g, do_clouds, do_haze)` | load all tables; once per process; `''`/non-positive → compiled default |
| `exort_run_column(state, result)` | one column |
| `exort_run_columns(n, states, results)` | serial batch |
| `exort_finalize()` | refuse further calls (no re-init in Stage D) |

Return codes: 0 OK, 1 already-init, 2 not-init, 3 bad args. Heating rates
come back in K s⁻¹ (×86400 for the file convention). Consumers:
`tools/exort_pytools/` (cffi binding + `verify_lib.py`; see its README) and
`tests/lib/` (`make run`; dlopen + netCDF-C, dimension-agnostic).

- **To reproduce the committed regression baselines through the library, pass
  the `_n84` stellar files** (`G2V_SUN_n84.nc`, `blackbody_3400K_n84.nc`) —
  the harness maps `_n68 → _n84` for exort runs.
- Tables are read-only after init and a column solve writes no module-scope
  state, so concurrent column solves are safe. `exort_run_columns` is itself
  OpenMP-parallel (results bitwise independent of thread count).
- `exort_set_percol_seed(enable)` switches on the per-column MCICA seed (see
  Multi-column I/O below).
- The per-column solve (`run_one_column`) lives in
  `source/src.main/exort_column_run.F90`, shared with the executable's column
  loop in `main.F90` — the two entry points cannot drift (Stage E1).

## `aerad_driver` argument contract (since Increment 1, 2026-07-02)

Mandatory positional core: thermo state, geometry, albedos, all 10 gas MMRs
(gas off = zero MMR; `calc_opd_gas` short-circuits zero-abundance gases).
Optional keyword tail after `sol_toa`: H₂O clouds (`ext_cicewp/cliqwp/cfrc/
rei/rel`, all five or none), CO₂ ice clouds (`ext_cicewp_co2/rei_co2`, both
or neither), CARMA haze (`ext_carmammr`), `ext_srf_emiss`, and per-column
dry-air properties `ext_mwdry`/`ext_cpdry` (added Stage E; absent → the
physconst module values, i.e. the CAM path — `physconst_setgas` is deleted
and 1-D callers always pass them). **Pass by keyword only; append new
optionals at the end** so positional call sites never break.
Flag on + arg absent → kernel not called (zero contribution). Kernels
(`calc_opd_*`) keep all-mandatory args. Full contract in the driver header
(`source/src.main/exo_radiation_mod.F90`).

---

## Runtime namelist `user_nl_exort` (`&exort_config`)

Overrides three formerly-compile-time constants without a rebuild. Place in the
run directory before invoking the executable; absent ⇒ compiled-in defaults.

| Variable | Meaning | v2 default (`source/exoplanet_mod.F90`) |
|----------|---------|------------------------------------------|
| `solar_file` | stellar spectrum NetCDF in `data/stellar/` | `'G2V_SUN_n84.nc'` |
| `shr_const_scon` | stellar constant ÷ 2 [W m⁻²] | `680.0` (≈ present Earth) |
| `exo_g` | surface gravity [m s⁻²] | `9.80616` (Earth) |
| `do_exo_clouds` | enable the cloud RT path (H₂O + CO₂ ice) | `.false.` |
| `do_exo_haze` | enable the CARMA haze aerosol RT path | `.false.` |
| `mcica_percol_seed` | offset the MCICA seed by column index (cloudy H₂O batches only) | `.false.` |

`do_exo_clouds` (added Stage C, 2026-06-28) became a runtime flag — when `.true.`
ExoRT loads the cloud Mie optics at init and reads condensate fields from the
input file; when `.false.` (default) the cloud path is skipped (bit-for-bit
cloud-free). It is read before the `initialize_cldopts` gate in `main.F90`, so the
namelist value applies. `do_exo_haze` (added Stage C3, 2026-07-01) works the same
way: when `.true.` ExoRT loads `data/aerosol/haze_n84_b40_mie.nc` at
init and reads `carmammr(pver,nelem,nbin)` from the input file; when `.false.`
(default) the aerosol path is skipped (bit-for-bit haze-free). `exo_pver` is **compile-time only** (not in the namelist);
active value `300`
(`integer, parameter :: exo_pver = 300`). The active config prints at startup
under `=== exort_config ===`.

---

## `tools/makeColumn.py` — build `RTprofile_in.nc`

All gas VMR flags default to the module-level `USER SETTINGS` block; absent/zero
species are omitted from the output (ExoRT reads them as zero). `n2vmr` is the
remainder (`1 - sum of others`).

```bash
python makeColumn.py --defaults                       # print current defaults
python makeColumn.py --profile TS273K --co2vmr 400e-6 --output RTprofile_in.nc
python makeColumn.py --profile TS273K --zero-h2o --co2vmr 0 --ch4vmr 1e-2 \
       --output fixture.nc                             # single-gas (CH4) fixture
```

Flags (all present in argparse): `--defaults`, `--output`, `--profile`
{`smart_2bar_t250`,`TS273K`,`US1976`}, `--co2vmr`, `--ch4vmr`, `--c2h6vmr`,
`--nh3vmr`, `--covmr`, `--h2vmr`, `--o2vmr`, `--o3vmr`, `--asdir`, `--asdif`,
`--aldir`, `--aldif`, `--coszrs`, `--grav`, `--srf-emiss`, `--zero-h2o`.

- H₂O is NOT a VMR flag — it comes from the profile's specific-humidity `q`
  array. `--zero-h2o` zeros it. `TS273K` is the 300-level (pverp=301) Earth
  moist-adiabat base.
- `--co2vmr 0` omits `co2mmr` entirely (needed to isolate non-CO₂ gases, since
  the default block carries 400 ppm CO₂).

**Clouds + surface emissivity (Stage C, optional, written only if set):**
- `--srf-emiss <float>` — broadband surface thermal emissivity scalar; written
  only if `!= 1.0` (ExoRT defaults to 1.0 when absent).
- Per-level cloud fields are API/USER-block driven (not CLI flags, since they are
  arrays): pass a `clouds` dict to `make_column()` or set the `CLOUDS` global.
  Supported fields (`CLOUD_FIELDS`): `cicewp`, `cliqwp`, `rei`, `rel`, `cfrc`
  (H₂O); `cicewp_co2`, `rei_co2` (CO₂ ice). The `add_cloud_layer(clouds, nlev,
  klo, khi, {field: value})` helper sets a field over a layer range. All-zero
  fields are skipped on write. Requires `do_exo_clouds=.true.` in ExoRT to take
  effect. See `tests/regression/fixtures/make_co2cloud_fixture.py` for a worked
  example.

---

## Multi-column I/O (`ncol`, Stage E1)

`RTprofile_in.nc` may carry an optional `ncol` dimension; `main.F90` then
solves every column in one invocation. **The column loop is OpenMP-parallel
(Stage E2)**: threads spawn only when `ncol > 1` (`if(ncol>1)` clause), the
thread count comes from `OMP_NUM_THREADS`, and threaded results are bitwise
identical to serial (gated by `multicol_check.py`). The large solver work
arrays are heap-allocated (E2 converted them from automatics when `-fopenmp`'s
implied `-frecursive` would have put ~120 MB per solve on the stack), so
default thread stacks suffice — no `OMP_STACKSIZE` needed. Build without
OpenMP via `make OMPFLAGS= <target>`.

- **Layout (Fortran order = trailing column dim):** scalars `ts/ps/coszrs/
  asdir/asdif/aldir/aldif/mw/cp[/srf_emiss]` → `(ncol)`; mid-layer and
  interface profiles → `(pver,ncol)` / `(pverp,ncol)`; `carmammr` →
  `(pver,nelem,nbin,ncol)`. In Python/C order that is `('ncol','pver')` etc.
- **Absent `ncol` = classic single-column file**, bit-for-bit the pre-E1 path;
  the output keeps the legacy layout (no `ncol`, scalars on `ONE`).
- `RTprofile_out.nc` mirrors the input: with `ncol`, every output variable
  gains the trailing column dimension (`LWUP(pverp,ncol)`, `FSDTOA(ncol)`, …).
- **`tools/stackColumns.py col1.nc col2.nc [...] -o out.nc`** stacks
  single-column files: optional variables are unioned across inputs and
  zero-filled (1.0 for `srf_emiss`) where a file lacks them; all inputs must
  share `pver`/`pverp`.
- **Columns share the process-level runtime config** (`solar_file`,
  cloud/haze flags) — but gravity and insolation may vary per column
  (Stage E2b): optional scalar input variables `grav` [m s⁻²] and `scon`
  [W m⁻², true value ÷ 2] override the namelist `exo_g`/`shr_const_scon`
  for that column (absent or ≤0 = namelist value; write them with
  `makeColumn.py --write-grav` / `--scon`). Per-column gravity is exact
  (same arithmetic as the namelist path); per-column `scon` rides the
  driver's `msdist` flux factor and agrees with a namelist-configured run
  to ~1e-15 relative (rounding order only) — gated at 1e-10 by
  `percol_config_check.py`. The stellar *spectrum shape* is per-process by
  design; run mixed-star sets as an outer loop of processes.
- **Per-column MCICA seed (Stage E2, opt-in):** by default every column's
  stochastic cloud generator uses the same constant seed (9404), bit-for-bit
  the legacy behavior. Setting `mcica_percol_seed = .true.` in `user_nl_exort`
  offsets the seed by the column index (column 1 stays bit-identical to
  legacy), decorrelating the MCICA subcolumns across batch columns. Only
  affects cloudy H₂O runs; enabling it is a rebaseline decision. Library
  equivalent: `exort_set_percol_seed(1)`.
- **Cloud fraction input:** `cfrc(pver[,ncol])` is now an optional input
  variable (E2; previously the 1-D deck could not drive H₂O cloud fraction —
  it was silently zero, so H₂O condensate produced no MCICA cloud). Absent →
  zero, as before.
- Acceptance check: `tests/regression/multicol_check.py` runs three
  same-config cases (TS250K/TS300K/TS340K × G2V) singly and as an `ncol=3`
  batch, requires exact equality (max |Δ| = 0 on every output variable), and
  re-runs the batch at `OMP_NUM_THREADS=8` requiring exact equality with the
  serial batch. Requires `run/exort.exe` already built.

---

## Regression suite (`tests/regression/run_regression.py`)

Builds and runs **`run/exort.exe`** (84-band grid, HITRAN-2024 native gases)
and compares flux/heating/spectral outputs to committed golden baselines at
`rtol = atol = 1e-3` (`DEFAULT_RTOL`, `DEFAULT_ATOL`). **32 cases**, all
`pver=300`, defined in `build_cases()`; current record in
`tests/regression/REGRESSION_STATUS.md`. Use `USER_FC=gfortran` on Apple
Silicon.

```bash
USER_FC=gfortran python run_regression.py            # all cases vs baselines (32/32)
python run_regression.py --list                      # list case names + physics
python run_regression.py --cases TS300K Mars
USER_FC=gfortran python run_regression.py --generate-baselines
USER_FC=gfortran python run_regression.py --exort h16  # HITRAN-2016 side-path
```

Flags: `--generate-baselines`, `--cases`, `--list`, `--rtol`, `--atol`,
`--verbose`, `--exort {h16,h24}`. Default (no `--exort`) is h24, the baselined
line list. `--exort h16` cannot be combined with `--generate-baselines`.

Cases:
- 12 Earth-like: TS250K–TS360K × {G2V, blackbody_3400K}.
- Mars: `2barCO2_dry_Mars_G2V` and `2barCO2_co2cloud_Mars_G2V` (cloudy).
- Haze: `TS300K_haze_G2V` and `TS300K_hazethick_G2V`.
- 16 minor-gas: CH₄/CO/NH₃/C₂H₆ × {realistic, elevated} × 2 stars.

A case dict with `"clouds": True` / `"haze": True` makes `write_namelist()`
emit `do_exo_clouds = .true.` / `do_exo_haze = .true.`. Both are toggled per
case at runtime, with no rebuild.

Compared variables:
- `FLUX_VARS = [LWUP, LWDN, SWUP, SWDN]`
- `SPECTRAL_VARS = [LWUP_SPECTRAL, LWDN_SPECTRAL, SWUP_SPECTRAL, SWDN_SPECTRAL]`
- `HEATING_VARS = [LWHR, SWHR]`
- `INTEGRATED_SCALARS = [FSDTOA]`

### `--exort h16` side-path mechanism

The default build needs no source edit (`kabs.F90` commits the h24 native
gases). `--exort h16` makes `build_exort()` do the following:
- string-swap the four native-gas filenames h24→h16 via `H16_TO_H24`;
- run `make exort` (honours `USER_FC`) with `_n84` stars (`N68_TO_N84_STAR`);
- compare against the h24 baselines;
- **restore** `kabs.F90` and `user_nl_exort` in a `finally` block.

Baselines are never written in this mode.

> **Maintenance:** the `H16_TO_H24` values are exact h24 filenames. If a
> k-file is renamed, update both `kabs.F90` and this map, or the swap
> silently no-ops for that gas.

---

## Per-gas sweep (`tests/regression/gas_sweep.py`)

Diagnostic (not pass/fail). Runs single-gas fixtures through three exes that each
encode a vintage/grid and tabulates OLR (`LWUP[0]`):

| exe | vintage / grid | star |
|-----|----------------|------|
| `n68equiv.exe` | HITRAN-2016 / 68-band | `G2V_SUN_n68.nc` |
| `n84equiv.exe` | HITRAN-2016 / 84-band | `G2V_SUN_n84.nc` |
| `exort.exe` | HITRAN-2024 / 84-band | `G2V_SUN_n84.nc` |

Reports **grid effect** `= n84 − n68` (both h16, expect ~0) and **line-list
effect** `= exort − n84` (h16→h24); flags `|line-list| > 1 W/m²`.

```bash
python gas_sweep.py                  # build all three, run all 12 fixtures
python gas_sweep.py --gases CO2 C2H6 # subset (substring on fixture tag)
python gas_sweep.py --no-build       # reuse existing run/*.exe
```

Flags: `--gases`, `--no-build`. Constants: `SCON = 680.0`, `GRAV = 9.80616`
(LW is insolation-independent in offline mode, so the star/scon don't bias OLR).
`GAS_FIXTURES` lists the 12 fixtures (H2O/CO2/CH4/C2H6/NH3/CO × realistic+elevated).

- Requires `n68equiv`/`n84equiv` `kabs.F90` on the flat `data/kdist/<gas>/`
  layout (h16) and `exort` at committed h24.
- `gas_sweep.py` and the `RTprofile_in_*` fixtures are force-tracked past the
  `tests/regression/*` gitignore rule (`git add -f`).
- NH₃/CO have only h24 fits (no h16). n84equiv was given NH₃/CO (`nspecies` 6→8)
  so all three codes are comparable; for NH₃/CO all three agree exactly.

---

## k-coefficient file format (`data/kdist/<gas>/*.nc`)

Flat per-gas layout (v2). Each file holds correlated-k tables on a fixed grid:

| dimension | size | source constant (`radgrid.F90`) |
|-----------|------|----------------------------------|
| NTemp | 17 | `kc_ntemp` (100–500 K) |
| NPress | 61 | `kc_npress` (0.01–10⁴ mb) |
| NGauss | 8 | `ngauss_8gpt` |
| NBins | 84 | `ntot_wavlnrng` |

- Variable read by ExoRT: **`data`** with shape `(NTemp, NPress, NGauss, NBins)`.
- **Coordinates are checked at load** (`check_kfile_grid`, 1-D and CAM). The
  run stops (`stop 1` / `endrun`) unless the `data` dims match and:
  - `Temperature` = `tgrid`;
  - `Pressure` = `pgrid` in **mb**;
  - `GaussWeights` = the g-interval midpoints.

  A missing or all-zero coordinate is also fatal, so new files must carry
  correct coordinates. `SpectralBands` holds only indices 1…N and is checked
  by count.
- ⚠️ The file variable named `GaussWeights` actually holds g-interval
  **midpoints** (`g_xpos_edge_8gpt + g_weight_8gpt/2`), not the integration
  weights. The integration weights live in
  `radgrid.F90` (`g_weight_8gpt = [0.30192, 0.27379, 0.22012, 0.14595, 0.04712,
  0.00686, 0.00363, 0.00061]`) and are correct.
- Filenames encode vintage/lineshape, e.g.
  `n84_8gpt_<gas>_hitran24_Nnu1e4_c<cutoff>_<lineshape>_q<n>_grrtm.nc`. Each gas
  dir holds both `hitran16` and `hitran24` files; `kabs.F90` selects by name.

---

## Gas species in `src.exort` (`radgrid.F90`)

`nspecies = 8`: indices `iH2O=1, iCO2=2, iCH4=3, iC2H6=4, iO3=5, iO2=6, iNH3=7,
iCO=8`. `kabs.F90` points H₂O/CO₂/CH₄/C₂H₆/NH₃/CO at **HITRAN-2024** (the
default since 2026-09-23; NH₃/CO never had an h16 table); O₂/O₃ are
HITRAN-2020. Same indexing in `n84equiv` (NH₃/CO added 2026-06-17 as a sweep
reference).

---

## Known limitations

Open work is tracked in `REFACTOR_PLAN.md`; the ones a user or reviewer is
most likely to hit are listed here.

- **`src.cam.exort` has not run in a real ExoCAM case.** It compiles against
  stubbed CESM in every CPP combination; the first HPC run is pending.
- **CO k-table edge bands 29/31** disagreed with an ad-hoc line-by-line
  check. The result is unconfirmed and the impact is small: 1 % CO moves OLR
  by −1.1 W/m² in total.
- **The fractal haze table** (`haze_n84_b40_fractal_interp.nc`) is
  provisional above band 68. The default Mie table covers all 84 bands.
- **No clear-sky / cloud-forcing `_CLD` double run in 1-D.** A cloudy run
  gives full-sky fluxes only.
- **Batch columns share one stellar spectrum shape.** Gravity and insolation
  may vary per column, but `solar_file` is per process by design. Run
  mixed-star sets as separate processes.
- **`libexort` is single-init.** One `exort_init` per process, with no
  re-init after `exort_finalize`.
- **Fatal data errors stop the process.** In the library that includes the
  caller. Only the k-file grid check exits nonzero (`stop 1`); the older
  loaders' bare `stop` exits 0.
- **Legacy references.** `n68equiv`/`n84equiv` `kabs.F90` are hand-edited to
  the flat `data/kdist/<gas>/` layout, so they won't run against the old
  `data/kdist/n68*` tree. They exist only as HITRAN-2016 references for
  `gas_sweep.py`.
- **Pre-v2 `3dmodels/` bundles are frozen** (`src.cam.n68equiv`, `.haze`,
  `.n84equiv`, `src.cam7.n68equiv`). They drift from `source/` by design and
  must not be regenerated.
