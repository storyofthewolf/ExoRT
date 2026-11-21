# ExoRT v2 — changes since v1.0.0

What changed between `v1.0.0` (the published Wolf et al. 2022 checkpoint, on
`main`) and v2 (the `refactor` branch), and the evidence for each
physics-affecting change. v1 is untouched: the `v1.0.0` tag still builds every
legacy bundle. Commit SHAs are on `refactor`; the full step-by-step narrative
lives in the commit messages.

Physics numbers below are 1-D, `pver=300`, from the regression suite
(`tests/regression/`). OLR = `LWUP` at the top of the model; SFC SW↓ =
`SWDN` at the surface.

---

## 1. Structure

- **One RT bundle, `src.exort`** (`make exort`). It replaces the five v1
  bundles. The legacy `n28archean`, `n42h2o` and `n68h2o` are gone from v2
  (`77fd3ee`); they remain in `v1.0.0`. `n68equiv` and `n84equiv` are merged
  into `src.exort` (`ff7a84e`): the 84-band grid (n68 + UV bins below
  0.24 µm) plus the NH₃ and CO absorbers, `nspecies = 8`. `n68equiv` and
  `n84equiv` survive only as HITRAN-2016 comparison references.
- **The merge is physics-neutral.** On HITRAN-2016 data, `exort` reproduced
  the `n68equiv` baselines on all 13 cases of the time. Longwave was
  bit-for-bit identical everywhere. Shortwave differed only on the G2V cases,
  by up to 0.37 % (SFC SW↓ −0.08 … −0.21 W/m²). That is the n68→n84 regrid of
  the structured solar UV spectrum; blackbody cases were identical
  (`d77edb9`, 2026-06-28).
- **Data layout.** `data/solar/` → `data/stellar/`; k-tables are flat per gas,
  `data/kdist/<gas>/` (`93d9f00`). Stellar files are per grid, `*_n84.nc`.
- **Runtime namelist `user_nl_exort`** (`&exort_config`): `solar_file`,
  `shr_const_scon`, `exo_g`, `do_exo_clouds`, `do_exo_haze` and
  `mcica_percol_seed`, formerly compile-time, are now set without a rebuild.
  `exo_pver` stays compile-time, and the standard is now 300.
- **Input files.** Every gas and cloud variable in `RTprofile_in.nc` is
  optional; an absent one is read as zero.

## 2. New capabilities

- **CO₂ ice clouds** (`calc_opd_cld_co2`) and **surface emissivity**
  (`srf_emiss`, was hard-wired to 1). The 1-D H₂O cloud optics loader, which
  had been commented out, was restored. Stage C, 2026-06-28.
- **CARMA haze** (`calc_opd_aero`, `carmammr` input). Optics are Mie spheres
  from the Khare et al. (1984) tholin indices over all 84 bands (`6c17346`,
  `d4ede97`). See `HAZE_1D_HOWTO.md`.
- **Library API** `libexort` with C and Python bindings: `exort_init` /
  `exort_run_column(s)` / `exort_finalize` (`c00ee1a`). See
  `tools/exort_pytools/README.md`.
- **Multi-column batches.** Optional `ncol` input dimension (`675803a`);
  OpenMP column loop whose results are bitwise independent of thread count
  (`b47c07c`); opt-in per-column MCICA seed (`c176bde`); per-column gravity
  and insolation (`934252c`). See `MULTICOLUMN_BATCH_HOWTO.md`.
- **The 1-D deck can drive H₂O cloud fraction** via an optional `cfrc` input
  (`c176bde`). Previously cloud fraction was silently zero in 1-D, so H₂O
  condensate produced no cloud.
- **v2 CAM bundle `3dmodels/src.cam.exort`** (`cbc15d6`). CO₂ clouds, haze
  and surface emissivity are enabled by `-cppdefs`. It is compile-checked
  only and has not yet run in a real ExoCAM case.
- **k-file grid check** (`f09b153`, `f82532c`). The 1-D and CAM loaders
  reject a k-table whose dimensions or Temperature/Pressure/g coordinates do
  not match the compiled grid.

## 3. Physics and data corrections

Each of these moved the regression numbers on purpose; the baselines were
regenerated in the same commit.

| Change | Effect | Commit |
|---|---|---|
| **HITRAN-2016 H₂O temperature-index error fixed.** The table used since 2020-11-05 (`d3d340e`) held k(T−25 K) in the slot labelled T, ~20 % under-absorbing. Its own 2020 Temperature coordinate reads 100, 100, 125 … 475; clean sibling HELIOS-K runs agree with the corrected table at the same T. Also fixed on `main` as `v1.1.0` (`3dd35ba`). Wolf et al. 2022 Fig. 1 used the shifted table. | Earth-like cases: OLR −1.5 … −4.5 W/m², SFC SW↓ −0.3 … −2.3. Mars unchanged. vs Yang et al. 2016 LBLRTM, OLR error goes from −0.8 … +2.8 to −0.6 … −3.9 W/m² | `5eb0202`, `184fef7` |
| **HITRAN-2024 becomes the default line list** for H₂O/CO₂/CH₄/C₂H₆ (NH₃/CO were only ever fitted to HITRAN-2024; O₂/O₃ stay HITRAN-2020). HITRAN-2016 remains reachable via `run_regression.py --exort h16`. | vs corrected HITRAN-2016: OLR +0.04 … +0.19 W/m², SFC SW↓ −0.02 … −1.12. The largest shift, TS360K G2V, is partly real near-UV H₂O lines: the HITRAN-2016 list stops at 25,711 cm⁻¹; the 2024 list reaches 42,000 cm⁻¹ | `2082527` |
| **HITRAN-2024 CO₂ re-fit.** The first fit omitted the sub-Lorentzian χ-factor below the band centre. | 2-bar CO₂ Mars OLR 92.94 (the first fit gave 47.8) vs 92.74 on HITRAN-2016 | `a85b473` |
| **HITRAN-2024 C₂H₆ re-fit.** The first table had been generated from the CH₄ line list. | Median k ratio vs HITRAN-2016: 0.24 → 1.00 | `1a536b7` |
| **Haze extinction units.** A leftover `×0.1` made haze optical depth 10× too small. | Haze cases only | `d4ede97` |
| **3-D haze kernel off-by-one.** The published haze bundle paired layer k−1 mass with `pdel(k)` and never filled the bottom level. The 1-D kernel and `src.cam.exort` use the correct mapping. | Haze, 3-D | 2026-07-01 |
| **O₂/O₃ k-table Pressure coordinate** was stored in bar under an `mb` label. Metadata only; the data are unchanged. | None | `f09b153` |

## 4. Tests

- **Regression suite** `tests/regression/run_regression.py`: 32 cases,
  rtol = atol = 1e-3. They cover the Earth-like TS sequence under two stars,
  Mars dry and CO₂-cloudy, thin and thick haze, and the CH₄/CO/NH₃/C₂H₆
  single-gas cases. The current record is in
  `tests/regression/REGRESSION_STATUS.md`.
- **Other gates:**
  - `multicol_check.py`: batch ≡ singles, threaded ≡ serial, exactly.
  - `percol_config_check.py`: per-column gravity exact, scon ≤ 1e-10.
  - `verify_lib.py` and `tests/lib`: the library reproduces the baselines.
  - `populate3Dmodels.py check --exort`: the CAM bundle matches `source/`.
  - `tests/cam_compile_check/run_all.sh`: the CAM bundle compiles against stubbed CESM.
- **`gas_sweep.py`** separates the grid effect (n84 − n68) from the
  line-list effect (HITRAN-2024 − HITRAN-2016) per gas.

## 5. Performance

Δ = 0 throughout: the regression suite is bit-identical to its baselines
after each change. Measurements come from `exort-bench`: 3312 real columns
from the ExoCAM case `exocam_ML_grp4_pt04` (`pver = 51`), with the production
Intel flags on AMD Milan.

- **Fused spectral-block solver** (`669af08`).
  - **Problem.** On a loaded node ExoRT was memory-bound, not compute-bound.
    The per-column solver working set (~64 arrays of `ntot_gpt × pverp`,
    15–18 MB) overflows a core's share of L3, so per-core cost rose 4–7×
    from one idle core to a full 126-core Milan node.
  - **Change.** The two-stream chain now runs one spectral block at a time on
    cache-resident arrays (`DEVELOPER_NOTES.md`). It also skips the night-side
    SW solve, whose result was discarded.
  - **Result.** Full Milan node: n68equiv clear sky 91.0 → 16.3 s (5.6×),
    exort cloudy 121.1 → 32.0 s (3.8×).
  - **In ExoCAM** (4×5 aquaplanet, 126 ranks), with the same change ported to
    `src.cam.n68equiv`: 22.27 → 42.52 simulated years/day. Day-30 restarts
    were bit-identical.
- **Cached cloud radius-grid bounds** (`26274c0`).
  - **Problem.** `interpolate_cld` recomputed `minval`/`maxval` of the
    300-point ice grid on every call.
  - **Result.** Cloudy bench 7.96 → 4.38 s on one core.
  - **In ExoCAM:** cloud optics 11.1 → 7.8 s per 30 days, and 43.38 simulated
    years/day overall (1.95× the original). Real CAM clouds are sparser than
    the bench's monthly means.
- **Optics only for the bands the solver uses.**
  - **Problem.** Band optimization limits the solver to the bands that carry
    flux, but the gas, cloud and haze optics still computed every band. For
    an M dwarf the solver uses bands 1–61; exort's other 23 bands (27 %) had
    their optics computed and never read.
  - **Change.** `calc_opd_gas`, `calc_opd_cld_h2o`, `calc_opd_cld_co2` and
    `calc_opd_aero` take an optional band range (default: all bands), and
    the driver passes the solver's range. MCICA still draws every random
    number in order, so the cloud subcolumns are unchanged; only the overlap
    and property assignment are limited to the range.
  - **Result.** Δ = 0 in the regression suite (`-O0`). With
    `-ffp-contract=off` at `-O2` the output is bit-identical. With FMA
    contraction on, 65 of ~14 M spectral values differ by one float32 ulp,
    from different vectorization (the same class as `REFACTOR_PLAN.md` §3
    10.5). Cloudy bench: 4.35 → 4.10 s.
  - **n68 vs exort.** The optimized n68equiv and exort runs now do the same
    work (488 solver g-points; optics over bands 1–61). The legacy bundles
    accept the range and still compute every band.
- **Two leftovers found by the first ExoCAM run of `src.cam.exort`.** Found
  with GPTL timers against the fused n68 port.
  - **`rad_precalc` added the CO₂-ice and haze optics to every g-point and
    level even when those optics are off.** Their arrays are then all zero.
    Now guarded by per-column flags: adding exact zeros changes nothing, so
    the result is identical. Solver 27.3 → 25.8 s per 30 days.
  - **MCICA drew random numbers for every g-point, including those past the
    solver's range.** Draws are subcolumn-major, so stopping at the last
    used g-point leaves every earlier draw unchanged. MCICA 9.0 → 7.6 s.
    Outputs outside the range are no longer zero-filled; nothing reads them.
  - Bit-identical: regression 32/32; exort-bench cloudy run identical to the
    previous commit, serial and 6 threads.
  - **Continuum and CIA interpolators** (`interpH2Omtckd_ng` and the ten
    CIA pairs) computed every band internally. They now take the band
    range. The redundant full-array zero-fills of `tau_gas` (driver and
    `calc_opd_gas`) and of the H2O cloud optics (driver and
    `calc_opd_cld_h2o`) are dropped: in-range values are always fully
    assigned, and the driver zeroes the cloud arrays only for a column
    without the H2O cloud path. Shrinking those fills to the range was
    tried first; it was slower (strided section fills instead of one
    `memset`).
  - **Still open:** exort v2 radiation is ~9 s per 30 days (~20 %) slower
    than the fused n68 port. Cloud optics outside MCICA account for ~3–5 s
    and gas optics for ~1.4 s. Unattributed.
- **Where ExoCAM's time goes now** (atmosphere, 30 days):
  - dynamics 45 %
  - radiation 29 % (solver 24 s, gas optics 12 s, cloud optics 8 s)
  - stratiform 17 %
