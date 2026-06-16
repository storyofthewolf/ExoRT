# ExoRT Refactor Plan — v2.0.0

> **2026-06-16 rewrite.** This plan was substantially re-scoped after `v1.0.0`
> was frozen on `main`. The original multi-version refactor plan (preserved in
> git history at commit `3745405` and earlier) assumed all five RT versions
> (`n28archean`, `n42h2o`, `n68h2o`, `n68equiv`, `n84equiv`) would be carried
> forward and deduplicated. **That assumption is dropped.** v2 collapses to a
> single RT bundle. See "v2 strategy" below for the new direction and "What
> changed from the original plan" for the diff against prior intent.

## Versioning model

- **`main` = `v1.0.0`** — frozen at the published Wolf et al. (2022) checkpoint.
  Tagged `v1.0.0`. Bugfix-only maintenance. **All five legacy RT bundles live
  here.** Anyone who needs `n28archean` / `n42h2o` / `n68h2o` uses v1.
- **`refactor` = v2.0.0 development** — single-bundle, modernized. Breaking
  changes are allowed (new namelist, parameter→runtime demotions, I/O renames).

### v1 legacy regression checkpoint — MAINTAINER-OWNED PREREQUISITE

The regression harness (`tests/regression/`) exists only on `refactor`, never on
`v1.0.0`. Before the legacy bundles are deleted from `refactor`, the maintainer
will **manually** run a last formal checkpoint of the legacy bundles
(`n28archean`, `n42h2o`, `n68h2o`) against their published outputs and store the
results **in v1 (`main`)**, where those bundles belong — *not* back-port the
harness into v1.

- **Why manual:** the v2 harness depends on the new namelist and on
  parameter→float/int demotions that don't exist in v1; back-porting it is more
  work than a ~30-minute manual run-and-save, the same way the original v1
  outputs were captured.
- **Ordering (relaxed 2026-06-16):** the prune happens only on `refactor`; it
  does not touch `main`, and the legacy bundles remain fully recoverable from the
  `v1.0.0` tag and git history. So the v1 checkpoint is *not* a hard predecessor
  of the prune — the maintainer can capture it in `main` at any time. (Stage A
  was executed before the checkpoint on this basis.)
- This is the only piece of the prune that is a maintainer action rather than a
  code task in this repo.

---

## v2 strategy — one bundle to rule them all

v2's core thesis: **`n84equiv` cleanly supersedes `n68equiv`.** The only physics
distinction is the *input* band count (68 vs 84 spectral intervals — n84 adds UV
bins below 0.24 µm), and that distinction never reaches the working RT loop
because `optimize_band_sw` / `optimize_band_lw` in
`source/src.main/exo_init_ref.F90` already condense the full band set at runtime.
68 and 84 bands are never used in full; both are optimized down. So a single
wider-grid bundle covers every case the two equiv versions covered.

**Target: a single RT bundle `src.exort`, built by a single target `make exort`.**
It is:

```
n84equiv band grid (the wider superset)
  + NH3 / CO absorbers ported forward from n68equiv  (nspecies 6 → 8)
  + haze (CARMA aerosol) folded in from experimental
  + CO2 clouds folded in from experimental
  = src.exort   (one bundle to maintain)
```

`make n68equiv`, `make n84equiv`, `make n28archean`, `make n42h2o`,
`make n68h2o`, `make n68equiv_exp` are all **removed** from `build/Makefile`.
The single new target is `make exort`. No compat alias is kept — v2 is a clean
break.

### Grounded merge surface (diff `src.n68equiv` vs `src.n84equiv`, 2026-06-16)

Same 10 files in each version. Per-file divergence and what it actually is:

| File | Δ lines | Nature | v2 action |
|---|---:|---|---|
| `radgrid.F90` | 80 | `nspecies` 8 (n68: +NH3,CO) vs 6 (n84); band counts | **Real port.** Base on n84 grid, add NH3/CO index+name+arrays |
| `kabs.F90` | 43 | NH3/CO `dirk_*`/`k_*_file` params; n84 filenames | **Real port.** Add NH3/CO file params onto n84 base |
| `calc_opd_mod.F90` | 56 | NH3/CO VMR + column density + tau_grey; band dims | **Real port.** Add NH3/CO opacity path onto n84 base |
| `model_specific.F90` | 10 | NH3/CO entries in `k_major_data`/`k_grey_data` | **Real port.** Add 2 lines |
| `initialize_rad_mod_1D.F90` | 14 | NH3/CO k-file read blocks | **Real port.** Add 2 read blocks |
| `rad_interp_mod.F90` | 340 | **whitespace/indentation only** (verified: `interpO2O2cia` reindent) | Reconcile; expect ~0 logic delta |
| `cloud.F90` | 4 | filename suffix | Trivial |
| `initialize_rad_mod_cam.F90` | 18 | NH3/CO + band dims; **CAM-only** | 3-D path only (not `make exort`) |
| `spectral_output_cam.F90` | 548 | **band-count `addfld` blocks** (intervals 69–84); **CAM-only** | 3-D path only (not `make exort`) |

**Takeaway:** the substantive work is the NH3/CO replay onto the n84 grid —
exactly the 9-file "add a gas" pattern documented in `CLAUDE.md`. The two largest
raw diffs (`rad_interp_mod` 340, `spectral_output_cam` 548) are **whitespace and
CAM-only band extension** — noise, not work. NH3/CO is expected to port clean
because n84 predates that work; the merge replays a known delta onto a wider grid.

---

## What changed from the original plan

| Original plan element | v2 disposition |
|---|---|
| Carry all 5 RT versions; dedup `cloud.F90`/`kabs.F90` (Stage 5) | **Dropped.** Legacy 3 pruned; equiv 2 merged. Nothing left to dedup across versions. |
| `source/shared/` + compat-shim rewrite (Stage 3a/3b) | **Dropped as an architecture.** The conservative `src.misc` cleanup already removed the dead CESM glue, fixed the netCDF-v4 skew, and left `src.misc` at 718 lines of *live* shims. The full shared/+shim rewrite buys mostly cosmetics + double churn (shims added now, removed in old Stage 8) for foundation-wide risk. The two genuinely useful pieces (CPP-overridable `pcols`; per-column MCICA seed) are extracted surgically when the multi-column stage needs them — not via a global shim layer. |
| Multi-version regression matrix | **Narrowed** to the single `src.exort` bundle. |
| `ARGUMENT_HANDLING.md` cross-version signature guarding | **Superseded.** With N=1 there are no cross-version signatures to keep consistent. The one surviving idea — `optional` args for the cloud/aerosol merge so `main.F90` need not fabricate zero arrays — is folded into Stage C below. The doc carries a superseded banner. |
| Stage 4 library (iso_c_binding + Python) | **Kept**, rescoped to `src.exort`. |
| Stage 6 multi-column + OpenMP | **Kept**, rescoped to `src.exort`. |
| Stage 7 IDL → Python | **Kept** as an independent parallel track. |
| 3dmodels symlink dedup (Stages 1–2) + `populate3Dmodels` | **Kept** but rescoped: only the surviving cam bundles matter. `populate3Dmodels.py` was already rewritten (`146de33`, `48f3baa`) with `check`/`diff`/`regenerate`/`list` subcommands. |

---

## Done so far on `refactor` (state at 2026-06-16)

- **Regression harness** (`tests/regression/run_regression.py`): Python/numpy,
  rtol/atol = 1e-3, **n68equiv only**, 13 cases. Green. In v2 this re-targets
  `src.exort`.
- **`src.misc` conservative cleanup**: ~3,775 lines of dead CESM glue removed;
  `wrap_nf` trimmed to 5 used wrappers; build switched to the netCDF library's
  own `netcdf.inc` (fixed v4 skew, repaired the experimental build). `src.misc`
  now = 718 lines of live shims (14 files).
- **`io_1D.F90` consolidation**: merged `input.F90` + `output.F90`; extracted
  `read_namelist` / `print_diagnostics` from `main.F90`.
- **`populate3Dmodels.py` rewrite** with subcommands + `diff`; **regenerated all
  3dmodels bundles** from source (fixed a real stale-n84 bug:
  `k_major_data(iC2H6)=k_ch4`, restored `nspecies=6`, removed a non-UTF-8 byte).
- **physconst import fix**, **calc_gasopd signature fix** (unused `qnh3`/`qco`
  dummies in the 4 non-n68 versions — these dummies **disappear** when those
  versions are pruned/merged), **`exo_pver=300`** standard.
- **`make exort` does not exist yet.** Current targets are still the 6 legacy
  ones (`build/Makefile:159-164`).

---

## v2 stages (priority order, as confirmed)

Each stage is independently shippable and gated by the regression suite
(`src.exort` must reproduce reference within tolerance).

### Stage A — Prune legacy bundles  *(DONE 2026-06-16)*

**v1 checkpoint:** not a predecessor (see relaxed ordering above) — the prune is
`refactor`-only and the bundles survive in `v1.0.0`. Maintainer captures the v1
checkpoint in `main` on their own schedule.

**Delete from `refactor`:**
- `source/src.n28archean/`, `source/src.n42h2o/`, `source/src.n68h2o/`
- `3dmodels/src.cam.n28archean/`, `3dmodels/src.cam.n42h2o/`,
  `3dmodels/src.cam.n68h2o/`
- `build/Makefile`: `OBJS_N28ARCHEAN`, `OBJS_N42H2O`, `OBJS_N68H2O` and their
  `.exe` targets (lines 102–164 region)
- `tests/regression/`: any legacy-version cases (currently none — harness is
  n68equiv-only, so this is a no-op for the harness)
- Remove the now-dead unused `qnh3`/`qco` dummy args from any file that only had
  them to satisfy the pruned versions (the `e0fd01e` fix becomes unnecessary once
  only the equiv lineage remains).

**Verify:** remaining builds (`n68equiv`, `n84equiv`, experimental) still green;
n68equiv regression 13/13. `populate3Dmodels.py check` green on the reduced tree.

**Risk:** Low — deletion of preserved-elsewhere code.

### Stage B — Merge `n84equiv` + `n68equiv` → `src.exort`  *(the core v2 task)*

**Goal:** one bundle, one `make exort` target, n84 band grid + NH3/CO, runtime
optimizer condenses bands as today.

**Steps:**
1. Create `source/src.exort/` from `source/src.n84equiv/` (the wider grid is the
   base).
2. Replay the NH3/CO additions onto it using the 9-file "add a gas" pattern in
   `CLAUDE.md` (`radgrid`, `kabs`, `calc_opd_mod`, `model_specific`,
   `initialize_rad_mod_1D` — plus the CAM-side `initialize_rad_mod_cam`,
   `spectral_output_cam` for the 3-D path). nspecies 6 → 8.
3. Reconcile `rad_interp_mod.F90` (expected whitespace-only; confirm no logic
   delta before adopting either copy).
4. Add `OBJS_EXORT` + `exort : exort.exe` to `build/Makefile`; remove
   `n68equiv` / `n84equiv` targets.
5. Point the regression harness at `src.exort` (the n84-grid output for the
   existing 13 cases becomes the new baseline — **rebaseline expected**, since
   the band grid changes vs the n68 reference; this is an intended physics change,
   not a regression. Capture new goldens and note the rebaseline in the commit).
6. Update `3dmodels/`: collapse `src.cam.n68equiv` / `src.cam.n84equiv` into the
   `src.exort`-derived cam bundle; update `populate3Dmodels` MANIFEST/config.
   `src.cam7.n68equiv` (CESM3-planets external) and `src.cam.n68equiv.haze` are
   special — handle per their existing fork status, don't blindly collapse.

**Verify:** `make exort` green; harness green against rebaselined goldens;
`populate3Dmodels check` green.

**Risk:** Medium — the NH3/CO port and the band rebaseline. Bounded by the
harness and by the fact that the delta is a known, documented gas-add pattern.

### Stage C — Fold in haze + CO2 clouds  *(additive)*

**Goal:** absorb the two experimental physics threads
(`source/experimental/src.n68equiv_exp/{shr,shr2}`) into `src.exort`:
- **CO2 clouds** (`shr/`): the extended `aerad_driver` signature
  (`CICEWP_CO2_in`, `REI_CO2_in`, `SRF_EMISS_in`), `initialize_cldopts()`,
  `do_exo_clouds` flag, clear-sky + cloudy-sky double run; extended `output_data`.
- **Haze** (`shr2/`): CARMA aerosol/haze hookup (merges with
  `src.cam.n68equiv.haze`).

**Surviving idea from `ARGUMENT_HANDLING.md`:** make the new cloud/aerosol args
**`optional`** (legal today — `calc_*opd` are module procedures with explicit
interfaces). This lets `main.F90` and any non-cloud caller simply not pass them,
rather than fabricating zero arrays. This is now a clean-code choice *within* the
single bundle, not a cross-version drift guard (that problem is gone with N=1).

**Verify:** `make exort` with clouds off reproduces Stage B baseline (within
MCICA stochastic envelope where applicable); clouds-on path produces physically
sensible output; experimental `.F90~` backups resolved/committed before merge.

**Risk:** Medium — genuine new physics integration. The experimental threads are
in-flight; coordinate with their state before merging.

### Stage D — `iso_c_binding` library + Python binding  *(kept, rescoped)*

As the original Stage 4, but only for `src.exort`:
`source/src.main/exort_lib_mod.F90` (`exort_init` / `exort_run_column` /
`exort_run_columns` / `exort_finalize`, `data_root` string API),
`exort_column_mod.F90` (`column_state_t` / `column_result_t`),
`libexort.{so,dylib}`, `tools/exort_pytools/` cffi binding. Thread-safety
contract: tables read-only after init.

**Verify:** standalone C test + Python harness reproduce single-column reference
within 1e-10.

**Risk:** Low — pure addition; `aerad_driver` unchanged.

### Stage E — Multi-column batch + OpenMP  *(kept, rescoped)*

As the original Stage 6, for `src.exort`. Decouple `io_1D.F90` module-scope
state into `column_state_t`; multi-column NetCDF read/write with optional `ncol`
dim; OpenMP loop in `main.F90`; per-column MCICA seed (`set_nstep(9404+col)`).
**This is where the surgical `exort_grid.F90` (`pcols` CPP-overridable) and the
per-column seed — the only load-bearing pieces of the dropped `shared/` stage —
are introduced, scoped to exactly what this stage needs.**

**Verify:** 1-column matches reference; 64-column batch == 64 single runs within
1e-10; `OMP_NUM_THREADS=1` vs `8` within 1e-10.

**Risk:** Medium — touches `io_1D.F90` / `main.F90`; bounded by the legacy
single-column path + harness.

### Stage F — IDL → Python migration  *(independent parallel track)*

Unchanged from the original Stage 7 (plotting/utils → Mie/optics → CIA/MTCKD),
except validation targets only `src.exort` data. Move `.pro` files to
`tools/legacy_idl/`; delete after data regeneration is validated.

### Stage G — Cleanup  *(last)*

- Drop any remaining transitional scaffolding.
- I/O file naming normalization (`RTprofile_in.nc`/`RTprofile_out.nc` → new
  symmetric, consistently-cased names) — **v2-only breaking change**, one
  coherent commit with the harness proving equivalence before/after. Touches
  `io_1D.F90`, `main.F90`, README/CLAUDE, `tools/`, regression fixtures.
- Update `README.md` and `CLAUDE.md` to the single-bundle layout.

---

## Preserved unchanged through v2 (in-flight / external contracts)

- `source/experimental/src.n68equiv_exp/` — until Stage C consumes it. Do not
  rename/restructure mid-flight; resolve `.F90~` backups with the maintainer
  before merging.
- `3dmodels/src.cam7.n68equiv/` — CESM3-planets external; forked files stay real
  files. Stage B updates *source*, then regenerates/reconciles cam bundles via
  `populate3Dmodels`; cam7 forks are declared, not blindly collapsed.
- `3dmodels/src.cam.n68equiv.haze/` — merges in Stage C, not before.

## Reproducibility

Numerical equivalence within 1e-10 for refactors that should not change physics
(Stages D, E, F, G). **Intended physics changes that rebaseline the harness**
(Stage B band-grid change; Stage C clouds-on) are gated by capturing new goldens
in the same commit, with the rebaseline called out explicitly. No byte-identical
NetCDF requirement.

## Session handoff protocol

Each session ends by updating the "Done so far" section above (or `CLAUDE.md`'s
Session Handoff) with: what completed (file paths), what tests pass, branch/commit
state, surprises, and where the next session starts. A new session reads this
plan, then `git log` the last few commits, then runs the regression suite to
confirm a known-green start before stacking new work.
