# CESM.CAM7.md — `src.cam7.n68equiv` assessment

Findings from analyzing `3dmodels/src.cam7.n68equiv` (the CESM3 bundle) against
`3dmodels/src.cam.n68equiv` (the CESM1.2.1 bundle) and the canonical `source/`
tree. Two separate concerns: (A) how to bring cam7 onto a maintainable footing
inside ExoRT, and (B) what cam7 improvements are worth back-porting to
`cam.n68equiv` / `source/`.

> Date of assessment: 2026-06-13. Branch: `refactor`.
>
> **Update 2026-06-13:** the `pdel` item in Part B has since been investigated
> and resolved — `pdel` was a vestigial (unused) argument and has been removed
> from `calc_gasopd` across all five 1-D RT versions and re-synced into the five
> non-haze `3dmodels/src.cam.*` CAM1 bundles. See "Resolved: the `pdel`
> argument" below. `src.cam7.n68equiv`, `src.cam.n68equiv.haze`, and
> `source/experimental/src.n68equiv_exp` were deliberately left untouched.

---

## Background: the two bundles

Both `3dmodels/src.cam*.n68equiv` dirs are **file copies** (not symlinks) of the
ExoRT RT code, staged for a GCM build.

- **`src.cam.n68equiv`** — links to **CESM1.2.1**. Its files are intended to be
  byte-identical to `source/src.main` + `source/src.n68equiv` (plus
  `sys_rootdir.F90`, which differs by machine path). Today the shared physics
  files (`radgrid.F90`, `exo_radiation_mod.F90`, `calc_opd_mod.F90`,
  `model_specific.F90`, …) are **0-diff against `source/`**. This is the
  single-source-of-truth invariant that makes it easy to maintain (and is
  enforced by the `populate3Dmodels.py` sync rule in CLAUDE.md).

- **`src.cam7.n68equiv`** — links to **CESM3 (development alpha)**. Authored by
  an NCAR CESM software engineer whose goal was importing ExoRT *into CESM3*,
  not interfacing cleanly with the ExoRT `source/` conventions. As a result it
  has diverged on three independent axes at once, which is why it is hard to
  keep in sync with the rest of the code.

### File-inventory differences

| | cam.n68equiv (CESM1) | cam7.n68equiv (CESM3) |
|---|---|---|
| `kabs.F90` | present | **dropped** (`use kabs` commented out everywhere) |
| `initialize_rad_mod_cam.F90` (k-file reader) | present | **dropped** (glue moved to CESM3-planets external) |
| `exo_radiation_cam_intr.F90` (CAM↔ExoRT glue) | present | **dropped** (glue moved to CESM3-planets external) |
| `cloud.F90` | present | dropped from bundle |
| `sys_rootdir.F90` | present | dropped |
| `exoplanet_mod.F90` | inherited from `source/` | **added** as a bundle-local pared-down copy |
| `README` | plain text | `README.md` (documents CESM3-planets integration) |

The k-coefficient arrays (`k_h2o`, `k_co2`, …) are still **declared** in cam7's
`radgrid.F90`, but the module that **reads** them
(`initialize_rad_mod_cam.F90`) was removed from the bundle — per the cam7
README, that responsibility now lives in the **NCAR/CESM3-planets** external
repository, which pulls ExoRT in as a native external component (no `-usr_src`
flag or `SourceMods` copy required).

### Content divergence of shared files (cam1 → cam7)

`rad_interp_mod.F90`, `planck_mod.F90` are identical. The rest differ:

| file | diff size | nature of divergence |
|---|---|---|
| `calc_opd_mod.F90` | ~497 lines | robustness guards (portable) + missing NH₃/CO (regression) + CESM3 imports |
| `spectral_output_cam.F90` | ~2406 lines | `cam_history` API: `horiz_only` replaces deprecated `phys_decomp` |
| `mcica_random_numbers.F90` | ~240 lines | `iso_fortran_env` integer types for CESM3 build env |
| `mcica.F90` | ~141 lines | obsolete `kabs` import cleanup, integer types |
| `exo_radiation_mod.F90` | ~33 lines | `SHR_CONST_G => gravit` gravity aliasing + missing NH₃/CO + `pdel` note |
| `radgrid.F90` | ~12 lines | `nspecies=6` vs `8` (missing NH₃/CO) |
| `rayleigh_data.F90` | ~30 lines | missing CO/NH₃ depolarization + A/B constants (formatting) |
| `model_specific.F90` | ~9 lines | missing NH₃/CO k-array assignments; `use kabs` commented |

The divergences sort into **three categories**:

1. **Legitimate CESM3 infrastructure** — must be preserved, must NOT go to cam1.
2. **Numerical-robustness improvements** — physics-version-agnostic, worth
   back-porting (see Part B).
3. **Stale physics** — cam7 predates the 2026-04-27 NH₃/CO gas additions; this
   is cam7 *behind* `source/`, and flows `source/` → cam7, not back.

---

## Part A — Incorporating cam7 into a maintainable pathway

**Root cause of the maintenance pain:** cam7 broke the single-source-of-truth
invariant that cam.n68equiv obeys, on all three axes — drifted content, missing
physics, and a *restructured file/ownership boundary* (glue moved into the
CESM3-planets external). The normal "copy from `source/`, diff must be empty"
sync rule cannot be applied to it as-is.

**Immediate constraint (current priority):** keep `src.cam7.n68equiv`
**buildable for CESM3 exactly as it is**. Do not edit it yet. A full disposition
is **gated on cloning `NCAR/CESM3-planets`** and inspecting what the external
actually supplies vs. expects ExoRT to provide.

**Target end-state:** make cam7 obey a sync discipline like cam.n68equiv's, where
the *only* tracked divergence is the deliberate CESM3 infrastructure delta. Get
there by separating the three categories so only category 1 remains:

1. **Re-baseline cam7 physics onto current `source/`** (eliminates category 3).
   Port the NH₃/CO additions into cam7 so `nspecies=8` and the
   `k_nh3`/`k_co`, Rayleigh, and `model_specific` entries match `source/`. This
   removes the largest part of the diff and is required anyway for NH₃/CO
   science in CESM3.

2. **Isolate the CESM3 infrastructure behind CPP guards** (makes category 1
   mergeable). Wrap the gravity aliasing (`SHR_CONST_G => gravit`), import
   lists, vacuum cutoff, and `horiz_only`/`phys_decomp` swap in
   `#if defined CESM3 … #else … #endif` inside the **single `source/` file**, so
   both bundles regenerate from one source of truth with different `-D` flags.

3. **Decide the ownership boundary for the dropped glue files** — the one real
   open decision, **deferred** pending the CESM3-planets clone:
   - **(A) CESM3-planets external owns the glue** — formalize the README's
     claim; ExoRT ships only version-specific physics; the sync gate checks only
     those files. Less for ExoRT to carry, but the CAM interface is no longer
     testable from this repo.
   - **(B) Bundle owns the glue (symmetric with cam1)** — re-add CESM3 versions
     of the glue into both cam7 and `source/`, CPP-guarded. ExoRT stays
     self-contained; cost is maintaining CESM3 glue not compilable standalone.

   Note: cam7's added `exoplanet_mod.F90` is the standalone-radiation subset of
   ExoCAM's module — it should resolve to `source/exoplanet_mod.F90` (or a
   documented CESM3 variant), not a third independent copy.

4. **Extend the sync tooling.** Teach `tools/populate3Dmodels.py` a `cam7`
   target with a per-file rule table: "diff-clean against source" for physics
   files, "CPP-variant of source" for the infra files, "ignored /
   externally-owned" for whatever step 3 excludes. Wire it into the deferred
   drift gate (see the `3dmodels-sync-gate` memory) so cam7 cannot silently rot
   again.

**Open / gating items:** clone & inspect `NCAR/CESM3-planets` to (a) confirm
exactly which files the external supplies vs. expects from ExoRT, which
empirically resolves the step-3 ownership question; (b) decide the
CPP-vs-separate-copy strategy for the infra deltas.

---

## Part B — Back-porting cam7 features to `cam.n68equiv` / `source/`

cam7 contains three genuinely valuable **numerical-robustness** changes that are
physics-version-agnostic and would harden the 1-D and CESM1 builds. They all
live in `calc_opd_mod.F90` (which is shared: `source/src.n68equiv` ≡ cam1), so
they belong in `source/` and propagate to cam1 via the sync tool.

### Tier 1 — Back-port these (real robustness fixes, portable as-is)

1. **Vacuum safety cutoff (model-top stability).** Two related guards:
   - Skip the entire layer when `pmid(ik) < 0.05 mb` (5 Pa): zero `tau_gas`,
     `tau_ray`, and **all ten CIA terms**, then `cycle`.
   - Skip when `pressure < log10pgrid(1)` (below the k-table pressure floor)
     instead of clamping vacuum onto "thick atmosphere" physics.
   - *Why it matters here:* the CO₂ k-table returns unphysical opacities
     (τ > 20,000) at very low pressure. The `exo_pver=300` standard profile
     reaches a fine model top that exercises exactly this regime. **Highest-value
     item to lift.**

2. **Out-of-range temperature guards on the continuum + every CIA pair.** cam1
   currently clamps T to the table edge and computes anyway. cam7 wraps the
   MT_CKD water continuum and all five CIA blocks (N2-N2, H2-H2, N2-H2,
   CO2-CO2 LW+SW, CO2-CH4) in
   `if (T >= tgrid(1) .and. T <= tgrid(ntemp)) … else ans_cia = 0`, preventing
   cold-top extrapolation blow-ups.

3. **MT_CKD dry/cold skip with `itc` counter preserved.** cam7 disables the
   water-continuum loop when `u_h2o < 1e16 molec/cm²` or T is out of table
   range, but still advances `itc` so downstream band indexing stays aligned —
   a careful, correct optimization + safety change.

### Tier 2 — Do NOT back-port (CESM3-specific; would break CESM1)

`SHR_CONST_G => gravit` aliasing, `physconst` import restructuring, `horiz_only`
vs `phys_decomp`, `iso_fortran_env` integer types. Keep isolated to cam7.

### Tier 3 — Not improvements (cam7 is *behind*; these flow source → cam7)

Missing NH₃ + CO throughout (`nspecies=6`, no `qnh3/qco`, no `u_nh3/u_co`,
Rayleigh missing CO/NH₃ terms). This is cam7 lacking ExoRT's 2026-04-27 work,
plus whitespace/formatting noise and a trivial `SHR_CONST_PI,SHR_CONST_PI`
duplicate typo in cam1 worth fixing in passing.

### Resolved: the `pdel` argument (was "one item to verify") — 2026-06-13

cam7's `exo_radiation_mod.F90` carries a `pdel` note —
`!jt Indexing error, Need to pass pseudo layer version of pdel` — and passes a
pseudo-layer `pdel/100.0` to `calc_gasopd` where cam1 (≡ `source/`) passed
`ext_pdel/100.0`. The concern was that `source/` had a latent indexing bug.

**Verdict: the indexing is genuinely inconsistent but completely inert, because
`pdel` was never read inside `calc_gasopd`.** Investigation trace:

- `pdel` originates from `state%pdel` (3-D) or the input file (1-D), and reaches
  `aerad_driver` (`exo_radiation_mod.F90`). There it is used **only** to build
  `coldens`/`coldens_dry` and `dzc` (mass column density). That work is fully
  discharged before the `calc_gasopd` call.
- Inside `calc_gasopd`, `pdel` was a dummy argument with **zero reads in the
  body** in n68equiv / n84equiv / n42h2o / n68h2o. In n42h2o/n68h2o the only
  other reference was a *commented-out* `pdelik = pdel(ik-1)` block. In
  n28archean `pdel` fed a local `pdelik` that was itself **assigned but never
  read**.
- The continuum/CIA path-length term that historically may have used a
  `pdel`-derived thickness is now driven entirely by `pathlength` (the
  `zint`/`zlayer` geometric thickness), not `pdel`. `pdelik` is the orphaned
  fossil of that older scheme.

**Action taken (bit-for-bit):** removed `pdel` from the `calc_gasopd` signature
and declaration in all five `source/src.n*/calc_opd_mod.F90`; removed the
orphaned `pdelik` declaration/assignments (n28archean) and the dead commented
`pdelik` blocks (n42h2o, n68h2o); dropped the `ext_pdel/100.0` actual argument
from the single shared call site in `source/src.main/exo_radiation_mod.F90:644`.
`ext_pdel` itself remains in `aerad_driver` (still needed for `coldens`/`dzc`).
Re-synced the edited `calc_opd_mod.F90` + `exo_radiation_mod.F90` into the five
non-haze `3dmodels/src.cam.*` CAM1 bundles (all now diff-clean against
`source/`).

**Verification:** all 5 production 1-D builds compile green
(`USER_FC=gfortran`); the n68equiv regression suite passes **13/13 bit-for-bit
(Δ=0.000)** on every case, confirming `pdel` was inert.

**Deliberately NOT touched:** `src.cam7.n68equiv` (keep buildable for CESM3 as
is; its `!jt` note is now moot once it is eventually re-baselined),
`src.cam.n68equiv.haze`, and `source/experimental/src.n68equiv_exp` (its
`calc_gasopd` is called from within its own `calc_opd_mod.F90`, decoupled from
the shared `src.main` caller). These lingering bundles still carry their own
matched `pdel` call+signature pairs and remain internally consistent; fold the
`pdel` removal into them when they are merged.

### Recommended back-port sequencing

1. Port Tier-1 guards into **`source/src.n68equiv/calc_opd_mod.F90`**,
   **NH₃/CO-aware** — merge the guards *around* the existing 8-species code, not
   cam7's 6-species version.
2. Run `tests/regression/run_regression.py`; confirm Δ=0 on the 13 cases
   (guards are no-ops within table range). The Mars `2barCO2_dry` case at fine
   model top is the one that may legitimately change — that is the guard working.
3. Investigate the `pdel` note separately; rebaseline if it is a genuine fix.
4. `populate3Dmodels.py` re-propagates the hardened `calc_opd_mod.F90` to cam1
   (diff stays clean).

This gives `cam.n68equiv` and the 1-D builds all of cam7's robustness wins
without importing any CESM3 coupling.
