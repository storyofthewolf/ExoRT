# ExoRT 1-D regression harness

Runs a matrix of standard atmospheric profiles through the `n68equiv`
executable under multiple stellar spectra and compares the resulting
flux / heating-rate fields against committed golden baselines.

## Matrix

| | |
|---|---|
| Profiles | `TS250K`, `TS273K`, `TS300K`, `TS320K`, `TS340K`, `TS360K` |
| Stars | `G2V_SUN_n68.nc`, `blackbody_3400K_n68.nc` |
| Cases | 6 × 2 = **12** |

All cases hold `shr_const_scon = 680.0` and `exo_g = 9.80616` fixed; only the
stellar spectrum varies. The fixtures are `pver = 300` and require a build with
`exo_pver = 300` in `source/exoplanet_mod.F90`.

## Layout

```
tests/regression/
  run_regression.py   driver + comparator
  fixtures/           6 input profiles  (committed)
  baselines/          12 golden outputs (committed)
  outputs/            fresh run outputs (gitignored)
```

## Usage

Build first (`cd build && USER_FC=gfortran make n68equiv` on Apple Silicon),
so `run/n68equiv.exe` exists. Then:

```bash
cd tests/regression
python run_regression.py                  # run all 12, compare to baselines
python run_regression.py --profiles TS300K --stars G2V_SUN_n68.nc
python run_regression.py --generate-baselines   # (re)create golden baselines
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
