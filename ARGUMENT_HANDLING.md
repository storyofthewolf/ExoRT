# ARGUMENT_HANDLING.md — aerad_driver / calc_gasopd signature-guard design

> **SUPERSEDED (2026-06-16).** This document addresses how to keep the
> `aerad_driver` / `calc_gasopd` signatures consistent **across the five RT
> versions** (unused-dummy-args vs `optional`-args, the stale-n84 bug, the
> cross-version interface module). v2.0.0 collapses ExoRT to a **single** RT
> bundle (`src.exort`) — see `REFACTOR_PLAN.md`. With N=1 there are no
> cross-version signatures to keep consistent, so the core problem here is
> **dissolved**: D1, D4, the rollout fan-out (steps 2–4), and the entire "Open
> questions" section no longer apply, and the `qnh3`/`qco` unused dummies get
> deleted with the pruned versions.
>
> **What survives:** (1) the `optional`-args idea (D2) keeps independent value
> for folding the haze / CO2-cloud args into the single bundle so `main.F90`
> need not fabricate zero arrays — folded into Stage C of `REFACTOR_PLAN.md`;
> (2) the `populate3Dmodels.py check` sync-gate (step 5) still guards
> source↔3dmodels parity, orthogonal to version count. The rest is retained
> below for historical context only.
>
> **IMPLEMENTED (2026-07-02).** D2 landed in revised form on `refactor`:
> `aerad_driver` now has a **mandatory positional core** (thermo/geometry/
> surface/all 10 gases — a gas is switched off by zero MMR, not by an absent
> argument) plus an **optional keyword tail** placed after the outputs (H2O
> clouds, CO2 clouds, CARMA haze, srf_emiss; keyword-passed, append-only, so
> positional call sites never break when physics is added). Kernels
> (`calc_opd_*`) keep all-mandatory arguments and are gated at the call site.
> This revises D2's gas split (qH2O/qCO2/qN2-mandatory was a five-version-era
> choice). Full contract in the `aerad_driver` header
> (`source/src.main/exo_radiation_mod.F90`) and the Increment-1 entry in
> `REFACTOR_LOG.md`; 3dmodels call sites adopt it during the 3-D port.

---


Session handoff covering: (1) `populate3Dmodels.py` `diff` subcommand, (2) a
full source↔3dmodels drift diagnosis + regeneration, and (3) the design
discussion for guarding the `aerad_driver` / `calc_opd` signatures against
future drift. Items (1) and (2) are **done and committed**; item (3) is
**design-only, nothing written yet.**

Branch: `refactor`. Last commits this session:
- `48f3baa` tools(populate3Dmodels): add diff subcommand; switch verbs to subcommands
- `7601537` 3dmodels: regenerate n28archean/n42h2o/n68h2o/n84equiv bundles from source

---

## Goal

Make the radiative-transfer call chain robust against signature drift as new
`aerad_driver` arguments are added over time — concretely, the *upcoming*
refactor adds **CO₂-cloud** and **aerosol** arguments that go **only into
n68equiv** (consumed by `calc_cldopd` / `calc_aeropd`), while the legacy
versions (n28archean, n42h2o, n68h2o) must remain runnable without growing
matching plumbing. The user wants something **more rigorous** than the current
"mandatory unused dummy argument" convention, and wants future signature
changes to be **easier / harder to get wrong**.

---

## The call chain (established facts)

```
main.F90  ────────────┐
                       ├─► aerad_driver (SHARED, src.main/exo_radiation_mod.F90)
exo_radiation_cam_intr ┘        │  builds pverp-sized q* arrays from ext_* mid-layer inputs
   (SHARED, src.main)           ▼
                          calc_gasopd   (VERSION-SPECIFIC, src.<ver>/calc_opd_mod.F90)
                          calc_cldopd   (VERSION-SPECIFIC)
                          calc_aeropd   (VERSION-SPECIFIC, currently empty stub everywhere)
```

Two distinct signature boundaries:
- **External:** `aerad_driver`'s `ext_*` args — the cross-component contract
  called by `main.F90` and `exo_radiation_cam_intr.F90`. `aerad_driver` lives
  in **shared** `src.main`, so its signature is identical in every bundle.
- **Internal:** `calc_gasopd`'s `q*` args — version-specific. `aerad_driver`
  always builds the full `q*` arrays internally regardless of version.

Key code facts verified this session:
- `exo_radiation_mod` reaches `calc_gasopd` via `use calc_opd_mod` — i.e. as a
  **module procedure**, which already carries a compiler-generated **explicit
  interface**. ⇒ `optional` arguments are *already legal here today* with no new
  interface file.
- Source-side signatures are **already unified**: every production version's
  `calc_gasopd` takes the full arg list `(... qh2o, qco2, qch4, qc2h6, qnh3,
  qco, qo2, qo3, qh2, qn2, ...)`; legacy versions accept the extra gases as
  unused dummies (the `e0fd01e` fix). `calc_cldopd` also already has a uniform
  signature across all 5 versions; `calc_aeropd( )` is a uniform empty stub.
- In n68equiv's `calc_gasopd`, each situational gas is dereferenced in exactly
  **one** place — a single VMR computation per level (lines ~191–196). So a
  `present()` guard is a trivial, localized change:
  `xvmr = 0; if (present(qx)) xvmr = qx(ik)*mwdry/mwx`.

---

## Key decisions made (and why)

### D1. Keep one universal `aerad_driver` signature; do NOT split it per version
`aerad_driver` lives in shared `src.main`. Splitting it per-version would
un-share the very file whose sharing was just restored by the regeneration,
multiplying the drift surface. The universal-signature design is sound; the
only failure mode is the *physical N× duplication* of the shared file going
stale.

### D2. `optional` arguments are the right primitive (user's insight)
Convert *situational* gas args to `optional` (keep `qN2, qH2O, qCO2` mandatory —
present in virtually every sim). Why this is better than mandatory dummies:
- Encodes physical truth in the type system (always-present vs situational).
- Removes the fake zero-array burden in `main.F90` / legacy callers — an absent
  gas is simply not passed, not passed as a zero array.
- New absorbers/clouds/aerosols become **purely additive** optional args:
  existing call sites that don't pass them keep compiling unchanged. This is the
  biggest "make future drift easier" win.
- Folds the future CO₂-cloud/aerosol (n68equiv-only) case into the *same*
  mechanism — they become optional args too; legacy `calc_cldopd`/`calc_aeropd`
  never reference them. No special-casing.

### D3. A standalone abstract-interface module is NOT required for `optional`
**Mid-discussion correction (important):** I initially claimed an
`abstract interface` module (`rt_opd_interface.F90`) was needed and would have
caught the stale-n84 bug. Both are wrong:
- `optional` is already legal because `calc_gasopd` is a module procedure (D-fact
  above) — no new interface file needed to use it.
- The n84 bug got through because the **whole module pair** (`exo_radiation_mod`
  + `calc_opd_mod`) was copied stale *together*, so the call and the body
  matched each other at the old vintage. An interface wouldn't have caught a
  self-consistent stale pair. **The sync-gate would have.**

An abstract-interface module still has *some* value — it enforces that the five
versions' `calc_gasopd` conform to one contract **across independent builds**
(belt-and-suspenders). But it is not a prerequisite for `optional` and is the
larger/riskier change.

### D4. Rollout scope: n68equiv + src.main first, prove bit-for-bit, then fan out
User chose this over an all-five-at-once pass. Lower blast radius; prove the
mechanism end-to-end before touching legacy versions.

---

## Approaches considered and rejected

| Approach | Verdict | Why |
|---|---|---|
| Mandatory unused dummy args (status quo, NH₃/CO pattern) | superseded by D2 | Forces 1-D/legacy callers to fabricate & pass placeholder zero arrays; doesn't document non-universality. |
| Per-version `aerad_driver` signatures | rejected (D1) | Un-shares src.main; multiplies drift surface. |
| `rt_opd_interface.F90` as the *primary* guard / prerequisite | rejected as framed (D3) | Not needed for `optional`; would not have caught the stale-module-pair bug. Demoted to optional belt-and-suspenders. |
| All-five-versions in one pass | rejected (D4) | Larger blast radius before first green build/regression. |
| Contract-comment-only ("required dummy, do not remove") | rejected | Convention, not enforcement. |

---

## Open questions

1. **First increment, precise scope** (the unanswered AskUserQuestion the user
   paused to clarify). Candidates:
   - (a) Make only n68equiv's `calc_gasopd` situational args `optional` (with
     `present()` guards); leave `aerad_driver` + the call site passing all of
     them. Zero blast radius to the other 4; bit-for-bit trivially.
   - (b) Also add `rt_opd_interface.F90` now (forces all 5 to compile against it
     immediately — larger first step).
   - (c) Wire the sync-gate first (it's what would actually have caught n84).
2. **Which boundary to make optional first** — internal (`calc_gasopd` `q*`) vs
   external (`aerad_driver` `ext_*`). Internal is lower-risk; external is what
   lets `main.F90` drop fake zero arrays.
3. **What "rigorous" means to the user** — compile-time enforcement *within* a
   build, cross-version uniformity *across* independent builds, or guaranteed
   source↔3dmodels parity. These need different mechanisms.
4. **Decided-later (user):** whether future CO₂-cloud/aerosol args are `optional`
   vs mandatory dummies — defer until that work starts. (D2 strongly implies
   `optional`.)

---

## Concrete next steps (proposed increments, in order)

0. **Resolve open question 1** (pick the first-increment scope).
1. **n68equiv `calc_gasopd` → optional** for the 6 situational gases
   (`qch4, qc2h6, qnh3, qco, qo2, qo3`); keep `qh2o, qco2, qn2` mandatory.
   Guard the 6 VMR computations with `present()`. Leave `aerad_driver` + the
   call site passing all args (optional-ness invisible to the other 4 versions).
   **Gate: build n68equiv + regression 13/13 bit-for-bit (Δ=0).**
2. **(separate)** `aerad_driver` `ext_*` → optional; drop fake zero-array
   plumbing in `main.F90` / legacy paths.
3. **(separate)** Fan the `optional` pattern to the other four `calc_gasopd`.
4. **(separate, optional)** Add a shared `abstract interface` module for
   cross-version `calc_gasopd`/`calc_cldopd`/`calc_aeropd` conformance.
5. **(separate)** Wire `populate3Dmodels.py check` into pre-commit/CI — the
   real guard against the stale-module-pair drift class. See
   memory `3dmodels-sync-gate` (designed, deferred, enable scoped to n68equiv).

**Invariant for every increment:** physics unchanged ⇒ n68equiv regression
must stay Δ=0. Any shared-file edit must keep `populate3Dmodels.py check`
green (regenerate bundles after).

---

## Side note: tooling done this session (context, not pending)
- `populate3Dmodels.py` now uses **subcommands** (`check`/`diff`/`regenerate`/
  `list`) instead of `--`flags; added `diff` (unified diff of drifted files,
  exit 1 on drift) with lenient UTF-8 decode.
- Regenerated the four drifted bundles. This **fixed a real bug** baked into the
  stale n84equiv bundle: `k_major_data(iC2H6,...) = k_ch4` (copy-paste; source
  had since fixed to `k_c2h6`), plus restored `nspecies=6`/O₂-O₃ and removed a
  `0xc3` non-UTF-8 byte present in all four stale `exo_radiation_mod.F90` copies.
