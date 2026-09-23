# CESM.CAM7.md — `src.cam7.n68equiv` assessment

Findings from analyzing `3dmodels/src.cam7.n68equiv` (the CESM3 bundle) against
`3dmodels/src.cam.n68equiv` (the CESM1.2.1 bundle) and the canonical `source/`
tree. Two separate concerns: (A) how to bring cam7 onto a maintainable footing
inside ExoRT, and (B) what cam7 improvements are worth back-porting to
`cam.n68equiv` / `source/`.

> Date of assessment: 2026-06-13. Branch: `refactor`.
>
> **v2 context (refreshed 2026-09-23).** `src.n68equiv` and `src.cam.n68equiv`
> are now legacy. The canonical source is `source/src.exort` + `source/src.main`,
> and the CESM1 bundle that obeys the sync discipline is `3dmodels/src.cam.exort`
> (`populate3Dmodels.py check --exort`). Read "cam.n68equiv / cam1" below as the
> *pattern* cam7 should follow; the target for back-ports and for re-baselining
> cam7 is `src.exort`. None of Part A or Part B has been done yet; both are
> tracked in `REFACTOR_PLAN.md` (§3.2, §4.4).
>
> The `pdel` note in cam7's `exo_radiation_mod.F90` was resolved 2026-06-13:
> `pdel` was a vestigial argument and was removed from ExoRT's `calc_gasopd`.
> cam7 itself was left untouched.

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

### Recommended back-port sequencing

1. Port the Tier-1 guards into **`source/src.exort/calc_opd_mod.F90`**. Merge
   them *around* the existing 8-species code, not cam7's 6-species version.
2. Run `tests/regression/run_regression.py` and confirm Δ = 0 on all cases;
   the guards are no-ops within table range. The Mars cases at the fine model
   top are the ones that may legitimately change — that is the guard working,
   and it rebaselines explicitly.
3. `populate3Dmodels.py regenerate --exort` carries the hardened file into
   `src.cam.exort`; `check --exort` and the CAM compile gate must pass.

This gives `src.cam.exort` and the 1-D builds all of cam7's robustness wins
without importing any CESM3 coupling.
