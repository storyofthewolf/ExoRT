# ExoRT 1-D regression harness

Runs a matrix of standard atmospheric profiles through the `n68equiv`
executable under multiple stellar spectra and compares the resulting
flux / heating-rate fields against committed golden baselines.

## Cases

Each case carries its own fixture, stellar spectrum, insolation, and gravity, so
heterogeneous planets sit alongside the Earth-like suite. Cases are defined in
`build_cases()`; run `python run_regression.py --list` to see them.

**Earth-like temperature sequence** (12 cases) — `shr_const_scon = 680.0`,
`exo_g = 9.80616`, each profile run under both stars:

| Profiles | Stars |
|---|---|
| `TS250K`, `TS273K`, `TS300K`, `TS320K`, `TS340K`, `TS360K` | `G2V_SUN_n68.nc`, `blackbody_3400K_n68.nc` |

**Mars-like** (1 case) — `2barCO2_dry_Mars_G2V`: 2 bar dry CO₂,
`shr_const_scon = 451.166`, `exo_g = 3.711`, `G2V_SUN_n68.nc`.

**13 cases total.** All fixtures are `pver = 300` and require a build with
`exo_pver = 300` in `source/exoplanet_mod.F90`.

### Adding a case

Append a dict to `build_cases()` with keys `name`, `fixture`, `star`, `scon`,
`g`; drop the input NetCDF in `fixtures/`; then
`python run_regression.py --cases <name> --generate-baselines`.

## Layout

```
tests/regression/
  run_regression.py   driver + comparator
  fixtures/           input profiles   (committed)
  baselines/          golden outputs   (committed)
  outputs/            fresh run outputs (gitignored)
```

## Usage

Build first (`cd build && USER_FC=gfortran make n68equiv` on Apple Silicon),
so `run/n68equiv.exe` exists. Then:

```bash
cd tests/regression
python run_regression.py                  # run all cases, compare to baselines
python run_regression.py --list           # list case names + their physics
python run_regression.py --cases TS300K Mars   # subset by substring match
python run_regression.py --generate-baselines  # (re)create golden baselines
```

Exit status is non-zero if any compared case fails.

## Comparison

Tolerance is `rtol = 1e-3`, `atol = 1e-3` (override with `--rtol/--atol`),
applied with `numpy.isclose` over four field groups:

- **Integrated fluxes** — `FSDTOA`, plus TOA/surface endpoints of the flux arrays
- **Per-level fluxes** — full `LWUP`, `LWDN`, `SWUP`, `SWDN` vectors
- **Spectral fluxes** — `LWUP_SPECTRAL`, `LWDN_SPECTRAL`, `SWUP_SPECTRAL`, `SWDN_SPECTRAL`
- **Heating rates** — `LWHR`, `SWHR`

The harness preserves and restores any existing `run/user_nl_exort`.

### Inline glance metrics

Every row also prints two first-order diagnostics as `new vs base (Δ)`:

- **`LWUP_TOM`** — outgoing longwave at top-of-model (`LWUP[0]`)
- **`SWDN_SFC`** — downwelling shortwave at the surface (`SWDN[-1]`)

(Vertical index convention: 0 = top-of-model, −1 = surface.) These show at a
glance how far a run has moved from its baseline even when the case passes.

## Regenerating baselines

Only regenerate when an *intended* physics/code change alters the outputs.
Review the diff carefully — `--generate-baselines` overwrites the golden files
and a careless run will mask a real regression.
