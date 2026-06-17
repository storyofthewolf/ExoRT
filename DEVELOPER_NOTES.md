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
| `make n68equiv` | `run/n68equiv.exe` | legacy HITRAN-2016 reference, 68-band |
| `make n84equiv` | `run/n84equiv.exe` | legacy HITRAN-2016 reference, 84-band + UV |
| `make n68equiv_exp` | `run/n68equiv_exp.exe` | experimental |
| `make clean` | — | removes build artifacts and `run/*.exe` |

- Default compiler `ifort`; on Apple Silicon add `USER_FC=gfortran` to every
  `make`. Requires `nf-config` on PATH (NetCDF4 Fortran).
- Legacy `n28archean` / `n42h2o` / `n68h2o` targets were removed in v2 (present
  in the `v1.0.0` tag).

---

## Runtime namelist `user_nl_exort` (`&exort_config`)

Overrides three formerly-compile-time constants without a rebuild. Place in the
run directory before invoking the executable; absent ⇒ compiled-in defaults.

| Variable | Meaning | v2 default (`source/exoplanet_mod.F90`) |
|----------|---------|------------------------------------------|
| `solar_file` | stellar spectrum NetCDF in `data/stellar/` | `'G2V_SUN_n84.nc'` |
| `shr_const_scon` | stellar constant ÷ 2 [W m⁻²] | `680.0` (≈ present Earth) |
| `exo_g` | surface gravity [m s⁻²] | currently `3.711` (Mars; set during testing) |

`exo_pver` is **compile-time only** (not in the namelist); active value `300`
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
`--aldir`, `--aldif`, `--coszrs`, `--grav`, `--zero-h2o`.

- H₂O is NOT a VMR flag — it comes from the profile's specific-humidity `q`
  array. `--zero-h2o` zeros it. `TS273K` is the 300-level (pverp=301) Earth
  moist-adiabat base.
- `--co2vmr 0` omits `co2mmr` entirely (needed to isolate non-CO₂ gases, since
  the default block carries 400 ppm CO₂).

---

## Regression suite (`tests/regression/run_regression.py`)

Runs standard profiles through `run/n68equiv.exe`, compares flux/heating/spectral
outputs to committed golden baselines at `rtol = atol = 1e-3` (`DEFAULT_RTOL`,
`DEFAULT_ATOL`). 13 cases, all `pver=300`; defined in `build_cases()`.

```bash
python run_regression.py                  # all cases vs baselines
python run_regression.py --list           # list case names + physics
python run_regression.py --cases TS300K Mars
python run_regression.py --generate-baselines
python run_regression.py --exort h16      # build/run exort vs n68 baselines (equivalence)
python run_regression.py --exort h24      # build/run exort vs n68 baselines (line-list delta)
```

Flags: `--generate-baselines`, `--cases`, `--list`, `--rtol`, `--atol`,
`--verbose`, `--exort {h16,h24}` (mutually exclusive with `--generate-baselines`).

Compared variables:
- `FLUX_VARS = [LWUP, LWDN, SWUP, SWDN]`
- `SPECTRAL_VARS = [LWUP_SPECTRAL, LWDN_SPECTRAL, SWUP_SPECTRAL, SWDN_SPECTRAL]`
- `HEATING_VARS = [LWHR, SWHR]`
- `INTEGRATED_SCALARS = [FSDTOA]`

`compare_case(..., skip_spectral=)` drops `SPECTRAL_VARS` — set automatically in
`--exort` mode because 84-band rows can't align element-wise with 68-band
baselines.

### `--exort {h16,h24}` mechanism

`build_exort()` edits `source/src.exort/kabs.F90`: for `h16` it string-swaps the
four native-gas filenames h24→h16 via `H24_TO_H16`; for `h24` it leaves the
committed file. Then `make exort` (honors `USER_FC`), runs cases with `_n84`
stars (`N68_TO_N84_STAR`), and **restores** `kabs.F90` + `user_nl_exort` in a
`finally` block. Baselines are never written in this mode.

> **Maintenance:** `H24_TO_H16` keys are exact h24 filenames. If a k-file is
> renamed (e.g. the 2026-06-17 `hitran2024`→`hitran24` rename), update both
> `kabs.F90` and this map or the h16 swap silently no-ops for that gas.

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
  The reader uses only `data`; the embedded `Temperature`/`Pressure`/
  `SpectralBands`/`GaussWeights` arrays are NOT consumed.
- ⚠️ The file variable named `GaussWeights` actually holds g-interval
  **midpoints**, not the integration weights. The integration weights live in
  `radgrid.F90` (`g_weight_8gpt = [0.30192, 0.27379, 0.22012, 0.14595, 0.04712,
  0.00686, 0.00363, 0.00061]`) and are correct.
- Filenames encode vintage/lineshape, e.g.
  `n84_8gpt_<gas>_hitran24_Nnu1e4_c<cutoff>_<lineshape>_q<n>_grrtm.nc`. Each gas
  dir holds both `hitran16` and `hitran24` files; `kabs.F90` selects by name.

---

## Gas species in `src.exort` (`radgrid.F90`)

`nspecies = 8`: indices `iH2O=1, iCO2=2, iCH4=3, iC2H6=4, iO3=5, iO2=6, iNH3=7,
iCO=8`. Native gases (H₂O/CO₂/CH₄/C₂H₆/NH₃/CO) on HITRAN-2024; O₂/O₃ on
HITRAN-2020. Same indexing now in `n84equiv` (NH₃/CO added 2026-06-17 as a sweep
reference; its native gases stay HITRAN-2016).

---

## Known limitations

- **HITRAN-2024 k-coefficients are not validated** (2026-06-17). The h24 tables
  give non-physical LW for H₂O (~+12%), CO₂ (far-IR wing; 2-bar CO₂ loses ~48%
  OLR), and C₂H₆ (~4× too weak). CH₄/NH₃/CO are clean. The ExoRT code is proven
  correct (exort with h16 files reproduces n68/n84 bit-for-bit). Defect is in the
  offline HELIOS-K generation; under investigation. Treat `src.exort` LW results
  involving H₂O/CO₂/C₂H₆ as provisional until re-fit.
- `run_regression.py` is **not yet rebaselined** to `src.exort` (blocked on valid
  h24 files); it still runs `n68equiv.exe` for the 13-case golden comparison.
- `n68equiv`/`n84equiv` `kabs.F90` are hand-edited to the flat `data/kdist/<gas>/`
  layout — test scaffolding so the legacy codes run against v2 data; they will
  not run against the old `data/kdist/n68*` tree.
- Experimental haze / CO₂-cloud optics are not yet folded into `src.exort`
  (Stage C, planned).
