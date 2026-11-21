# ExoRT performance work, November 2026

What was done to make ExoRT fast inside ExoCAM, why it worked, how it was
measured, and what is left. Commit-level detail is in `CHANGES.md` §5 and the
commit messages. Every change here is Δ = 0 (bit-identical) against the
regression baselines.

## Result

ExoCAM `exocam_ML_grp4_pt04`: 4×5 aquaplanet, M-dwarf star (Teff 2850 K),
N₂ + CO₂, 51 levels. 30-day branch run from year 191 on one 126-core AMD Milan
node (Discover), ifort `-O2 -fp-model source`.

| Build | Throughput (sim yr/day) | Radiation per 30 days |
|---|---|---|
| Original (`main` `src.cam.n68equiv`) | 22.27 | 198 s |
| Same physics, fused solver + cloud-bounds cache (hand port) | **43.4** | 45 s |
| `refactor` `src.cam.exort` at `bc17591` | **41.7** | 46 s |

The model runs about 1.9× faster. Radiation fell from 64 % to 28 % of the
atmosphere's time. The FV dynamical core (`stepon_run1`, ~72 s, 44 %) is now
the largest cost.

## Diagnosis: ExoRT was memory-bound on a loaded node

The assumption in `OPTIMIZATION.md` was that the cost lay in the k-table
interpolation. Measurement said otherwise.

- **Offline bench.** 3312 real ExoCAM columns, run through the 1-D code. It
  reproduces production's per-column cost (~30 ms per column on a full node).
- **Load slowdown.** On a full Milan node, per-core cost was 4–7× higher than
  on one idle core. Sharing one table copy (OpenMP) vs private copies (MPI)
  made no difference.
- **Cause.** A cache-placement test isolated it. 8 processes packed on one
  8-core CCX (32 MB L3) ran 2.8× slower; spread one per CCX, 1.0×.
  - The two-stream solver's per-column working set (~64 arrays of
    `ntot_gpt × pverp` doubles, 15–18 MB) does not fit a core's share of L3.
  - The k-tables were ~4 % of single-core time.

## Changes (all on `refactor`)

| Commit | Change | Effect |
|---|---|---|
| `669af08` | **Fused spectral-block solver.** The whole two-stream chain runs one block of 16 g-points at a time, on block-sized arrays declared with global bounds `(gb0:gb1, …)` so loop bodies are unchanged. Night-side SW solve skipped. | Full Milan node: n68 5.6×, exort 3.8× faster |
| `26274c0` | **Cached cloud radius-grid bounds.** `interpolate_cld` recomputed `minval`/`maxval` of the 300-point ice grid on every call. | 46 % of the cloudy bench profile, gone |
| `ea2d7e5` | **Optics only over the solver's optimized band range** (gas, cloud, haze, MCICA assignment). | exort does the same work as n68 |
| `ce60fa3` | Skip the all-zero CO₂-ice/haze adds in `rad_precalc`. MCICA stops drawing after the last used g-point (draws are subcolumn-major, so earlier numbers are unchanged). | ExoCAM: solver −1.5 s, MCICA −1.4 s |
| `bc17591` | Band range passed into the continuum/CIA interpolators. Redundant full-array zero-fills of `tau_gas` and the cloud optics dropped. | **ExoCAM radiation 54.8 → 46.3 s** |
| `f06d2a7` | Fix: CO₂-ice density owned by ExoRT. The CAM compile gate now builds against ExoCAM's real `shr_const_mod`. | `src.cam.exort` builds in mainline ExoCAM |

## Lessons

- **Measure on a loaded node.**
  - An idle laptop core showed ~1 % for `bc17591`; ExoCAM on 126 ranks
    showed 8.5 s (18 % of radiation).
  - Cache-footprint effects (zero-fills, larger arrays) only bite when every
    core competes for L3 and DRAM bandwidth.
  - The fully loaded Milan bench (126 processes) is the right proxy.
- **Remove redundant memory traffic rather than shrinking it.** Resizing a
  whole-array zero-fill to a strided section was ~11 % *slower* (many small
  fills instead of one `memset`).
- **"Same code, more time" means footprint.** MCICA's code did not change in
  `bc17591`, yet it got 1.8 s faster.
- **Noise.** Two identical 30-day runs differ by up to ~2 s in total
  radiation and ~1 s per sub-timer. Repeat anything smaller before trusting it.
- **Ruled out for the v2-vs-n68 gap:** per-call allocations (allocate-once
  gave nothing on Linux), the Mersenne Twister, band ranges, cloud inputs and
  cloud amount, and compile flags.
- **Bit-identity under FMA.**
  - Moving code can change how the compiler contracts multiply-adds.
  - gfortran: compare with `-ffp-contract=off`.
  - Intel `-fp-model source` stayed bit-identical throughout.

## How it was measured

- **Offline bench:** `~/Research/models/exort-bench` (not in this repo).
  - Its `README.md` has every table.
  - `make_columns.py` builds the 3312-column input from a pt04 `cam.h0`.
  - `patch_fused.py`, `patch_cldbounds.py` are the mechanical patch
    generators. `patch_fused.py` also produced the `src.cam.n68equiv` hand
    port.
  - `discover/*.sh` are the Milan job scripts.
- **ExoCAM sandbox cases on Discover** (exocam-casemgr matrices
  `sbx_exort_timers.yaml`, `sbx_exort_v2.yaml`):
  - `sbx_exort_timers_n1`: the n68equiv hand port with GPTL timers.
  - `sbx_exort_v2_n1`: `src.cam.exort` with GPTL timers (`exo_gasopd`,
    `exo_cldopd`, `exo_mcica`, `exo_solve`, `exo_rt_fullsky`) in the case's
    SourceMods copy, and `profile_detail_limit = 1` in `user_nl_cpl`.
- **Star.** pt04's star exists only on the n68 grid.
  `tools/stellar_n68_to_n84.py` made the n84 copy: bands 1–68 identical, UV
  bins zero. This is valid for cool stars only.

## Left open

- **The remaining v2-vs-n68 gap,** ~1.6 s radiation per 30 days, spread
  across solver/MCICA/clouds and near noise. Plausibly the 84-vs-68-band
  footprint.
- **Band-major k-table:** deferred.
  - The 8 g-points of a band are 5.4 KB apart in `k_major_data`, so 32 cache
    lines per lookup.
  - It only touches gas optics (~11 s per 30 days, ~7 % of the atmosphere),
    so the ceiling is small.
  - If done, reorder in `model_specific.F90`'s existing startup copy, not in
    the data files.
- **`OPTIMIZATION.md` interpolation fixes:** not applied. Gas optics is now at
  parity with n68.
- **Other REFACTOR_PLAN §1.1 checks:**
  - the OpenMP NTHRDS 1-vs-2 identity check in CAM;
  - the clear-sky n68 vs exort physics comparison;
  - the gated physics (`-DEXORT_CO2CLD`, `-DEXORT_CARMA`, `-DEXORT_SRF_EMISS`)
    in a real case.
- **Production ExoCAM** still builds `main`'s `src.cam.n68equiv`, which has
  none of this. Moving a case to `refactor`'s `src.cam.exort` picks it all up.
  The n68 hand port lives only in exort-bench and the sandbox case.
