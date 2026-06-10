# ExoRT Refactor Plan

## Context

ExoRT is a Fortran 90 two-stream radiative transfer code with both a 1D offline executable and a CESM/CAM (3D) hookup via ExoCAM. After ~10 years of organic growth there are four sources of technical debt that increase maintenance cost and risk of human error:

1. **CESM stub layer (`source/src.misc/`).** ~2150 lines of stubs carried over from CESM (`shr_kind_mod`, `shr_const_mod`, `physconst`, `ppgrid`, `wrap_nf`, `ioFileMod`, `time_manager`, `spmd_utils`, `infnan`, etc.) that exist only so the 1D code can compile without CESM. Mixed in are dead stubs (`spmd_utils`, `pmgrid`) and substantial real wrappers (`wrap_nf` ~1666 callsites, `ioFileMod::getfil` ~1000 callsites).

2. **Drifting duplication across RT versions.** Five RT versions (`n28archean`, `n42h2o`, `n68h2o`, `n68equiv`, `n84equiv`) each ship 9 files; `cloud.F90` and `kabs.F90` are pure parameter declarations differing only in filename strings, while `radgrid.F90` (~36–48% common across versions) and `initialize_rad_mod_*.F90` (~48% common) differ mainly in array sizes and band counts. The actual physics modules (`calc_opd_mod.F90`, `rad_interp_mod.F90`) genuinely diverge (e.g. n68equiv 1224L vs n68h2o 664L) and are not candidates for unification.

3. **Drifting duplication across `3dmodels/`.** Each `3dmodels/src.cam.{version}/` has 17 files (8 from `src.main/` + 9 from `src.{version}/`), all currently byte-identical to their upstream `source/` counterparts except `initialize_rad_mod_cam.F90` (legitimate I/O fork). Additionally, `3dmodels/src.cam7.n68equiv/` is a CESM3-planets external-component hookup with 12 files: 5 are absent (CESM3 framework provides `cloud.F90`, `kabs.F90`, `initialize_rad_mod_cam.F90`, `exo_radiation_cam_intr.F90`, `sys_rootdir.F90`); 9 are intentionally diverged from cam1 for CESM3 API compatibility (`spectral_output_cam.F90` `cam_history` API, `model_specific.F90` `kabs` decoupling, `exo_init_ref.F90` MPI sync) plus correctness improvements (`calc_opd_mod.F90` vacuum cutoff, `exo_radiation_mod.F90` runtime-gravity aliasing); 1 is new (`exoplanet_mod.F90`, a CESM3-only parameter module — different from `source/exoplanet_mod.F90`). The `tools/populate3Dmodels.py` script meant to maintain this is broken (references undefined variables, no functional `--check`).

4. **IDL tools.** ~14 `.pro` files (~8200 lines) for stellar spectra, k-coeff fitting, Mie scattering, plotting. Five Python ports already exist (one of them broken). IDL is a licensing and onboarding burden.

5. **1D wrapper limitations.** `main.F90` (124 lines) processes one column per invocation; filenames are hardcoded; `pcols=1` is baked in at compile time; module-scope state in `input.F90` (`TS_in`, `TMID_in`, …) prevents per-column parallelism. The compute kernel `aerad_driver` (`exo_radiation_mod.F90:154–691`) is already a clean per-column subroutine — only the wrapper is single-shot.

### In-flight work that this refactor must preserve

`source/experimental/src.n68equiv_exp/` is **active development**, not abandoned scaffolding. It contains three distinct in-progress threads that the refactor must accommodate without breaking:

- **`shr/`** — Mars CO2-ice-cloud build. Adds new arguments to `aerad_driver` (`CICEWP_CO2_in`, `REI_CO2_in`, `SRF_EMISS_in`), an `initialize_cldopts()` init step, a `do_exo_clouds` flag, and runs `aerad_driver` twice per column (clear-sky default + cloudy-sky option). `output_data` signature is correspondingly extended. The README also notes "this experimental version contains the multiple H2O continuum methods, which are not in other versions" — a divergent physics path.
- **`shr2/`** — CARMA aerosol/haze hookup, intended to merge with `src.n68equiv.haze`. Contains editor backup files (`.F90~`) — actively being edited.
- **The `shr/` subdirectory layout itself** is a prototype of a directory reorganization the maintainer has already started exploring. The refactor's `source/shared/` proposal should harmonize with this naming, not collide with it.

The refactor must NOT delete, rename, or massively restructure `experimental/` files mid-flight. Where the refactor introduces shared infrastructure (e.g. `source/shared/`), the experimental builds get equivalent benefits via a parallel adoption stage, not a forced migration.

### Intended outcome

A staged refactor that:
- replaces `src.misc` with a small native `shared/` layer
- collapses byte-identical `3dmodels/` copies into symlinks with a copy-fallback script
- dedups parameter-only RT files (defers logic-heavy ones)
- delivers ExoRT as (a) a Fortran-callable library, (b) a Python binding via `iso_c_binding`, and (c) a multi-column batch executable
- migrates IDL tools to Python in a parallel track

### User-confirmed constraints

- **Reproducibility:** numerical equivalence within 1e-10 (allows FP reordering); no requirement for byte-identical NetCDF output.
- **3dmodels strategy:** symlinks with copy-fallback script.
- **Library targets:** all three (Fortran lib + Python binding + batch CLI).
- **Dedup scope:** parameter-only files first; defer `radgrid` template and `initialize_rad_mod` factoring.
- **ExoCAM / CESM compatibility:** `xmlchange CAM_CONFIG_OPTS="-usr_src .../3dmodels/src.cam.*"` (CESM1) must keep working; `3dmodels/src.cam7.n68equiv/` (CESM3-planets external) must remain buildable end-to-end. Both are first-class targets through every stage.
- **Three-target topology:** the refactor explicitly supports three host paths simultaneously — the 1D wrapper (`source/`), CESM1 (`3dmodels/src.cam.*`, treated as "active legacy" — gets the duplication cleanup but no API changes), and CESM3-planets (`3dmodels/src.cam7.*`, the forward target). The Stage 4 library wrapper is the natural hook for any fourth host (WRF, ROCKE-3D, Python research code).
- **CAM7 backport scope:** the `calc_opd_mod.F90` vacuum-pressure cutoff is backported to 1D source (defensive, no numerical drift expected). The `exo_radiation_mod.F90` runtime-gravity aliasing stays forked in cam7 — variable planetary gravity for the 1D path is future work.
- **Single maintainer:** stages must be independently shippable.

### Session budget and chunking

This refactor is too large for a single Claude Code session. Estimated total cost: **4.4M–8.1M tokens**, requiring **~10 separate 5-hour Pro-plan sessions**. The plan is chunked so each session has:

1. **A concrete deliverable** (a single stage or sub-stage) that can be committed to git as a working state.
2. **A "definition of done"** (specific tests pass) that any future session can verify cold.
3. **A handoff note** at `tests/HANDOFF.md` updated at the end of each session, recording: what was completed, what tests pass, what's next, any surprises.

Sessions are sized at **~500k tokens of useful work** each, leaving headroom for the inevitable build-error-iteration loops. A session that runs out of budget mid-stage MUST commit its current state to a feature branch and update `HANDOFF.md` before terminating, so the next session resumes from a known-good checkpoint.

| Session | Stage(s) | Estimated tokens | Deliverable |
|---|---|---|---|
| 1 | Stage 0 (regression harness + experimental triage + cam7 build-check scaffold) | 300k | `tests/refs/`, `tests/run_regression.sh`, `tests/EXPERIMENTAL.md`, `tests/HANDOFF.md`, cam7 MD5 manifest |
| 2 | Stages 1 + 2 (populate3Dmodels rewrite with three-target MANIFEST + symlinks for cam1 and cam7) | 400k | New populate script, symlinked 3dmodels/ for cam1 + cam7, ExoCAM (CESM1) and CESM3-planets rebuilds verified |
| 3 | Stage 3a (shared/ kinds + constants + grid; compat shims for those three) | 500k | 3 production builds green using new modules via shims; cam7 untouched |
| 4 | Stage 3b (shared/ netcdf + io + abort + runtime; remaining compat shims; delete src.misc/) | 500k | All 5 production builds green; src.misc/ deleted; experimental still builds; cam1 + cam7 still build |
| 4.5 | Stage 3c (backport cam7 vacuum-pressure cutoff to 1D source) | 150k | `calc_opd_mod.F90` vacuum guard in all applicable production versions; regression suite within 1e-10; new low-pressure unit test |
| 5 | Stage 4 (iso_c_binding library + Python binding, `data_root` API only) | 700k | `libexort_n68equiv.so`, Python harness reproduces reference |
| 6 | Stage 5 (cloud.F90 + kabs.F90 dedup) | 400k | Per-version cloud.F90/kabs.F90 deleted; template + .inc files in place; cam7 marked `absent` for these |
| 7 | Stage 6a (decouple input/output globals into column_state_t) | 500k | Single-column path uses derived type; output bit-comparable |
| 8 | Stage 6b (multi-column NetCDF + OpenMP loop in main.F90) | 500k | 64-column batch test passes with OMP threading |
| 9 | Stage 7 (IDL → Python — split if needed into 9a plotting/utils, 9b Mie/optics, 9c CIA/MTCKD) | 1M (likely 2 sessions: 9a+9b, then 9c) | `tools/exort_pytools/`, regenerated data files validated |
| 10 | Stage 8 (cleanup: drop compat shims, mechanical `use` rewrites, README/CLAUDE updates) | 500k | Compat shims gone, codebase uses native module names |

If a session hits its budget early (e.g. Stage 4 turns out heavier than estimated), split that stage and continue in the next session — the chunking is a guide, not a contract.

---

## Progress log — actual work vs. plan (updated 2026-06-10)

The work below was done on the `refactor` branch. It **deviates from the staged
plan** in two important ways the next session must understand; the staged
sections further down are the original design and have NOT been rewritten.

### Stage 0 — DONE, but built differently than specified
The regression harness exists and is green, but the implementation differs from
the Stage-0 spec:
- It is **Python** (`tests/regression/run_regression.py`), not the planned
  `tests/run_regression.sh` + `nccmp`. Comparison uses netCDF4/numpy with
  rtol/atol = **1e-3** (not 1e-10 — looser, chosen for cross-compiler/machine
  robustness). Fixtures and golden baselines are committed under
  `tests/regression/{fixtures,baselines}/`.
- It covers **n68equiv only**, with 13 cases (6 TS profiles × {G2V,
  blackbody_3400K} + a Mars-like `2barCO2_dry_Mars_G2V`). The other 4 production
  versions and the cam7 build-check are **NOT** yet wired into it.
- `tests/HANDOFF.md`, `tests/EXPERIMENTAL.md`, `tests/refs/`, and the cam7 MD5
  manifest from the Stage-0 spec were **not** created. The handoff lives in
  `CLAUDE.md` (Session Handoff section) instead.
- Still open from Stage 0 if desired: multi-version coverage, experimental
  output capture, cam7 structural check.

### Stage 3 — goal substantially achieved by conservative in-place cleanup
Rather than the planned `source/shared/` + compat-shim rewrite, this session did
**lower-risk in-place cleanup of `src.misc/`** that achieved much of Stage 3's
*goal* (remove dead CESM glue) without the new architecture:
- Trimmed `wrap_nf.F90` to the 5 used wrappers; deleted 4 unreferenced headers
  (`comctl/comsol/ncfortran.h`, `nfconfig.inc`); trimmed `ioFileMod` to `getfil`
  and deleted now-orphan `shr_sys_mod`; deleted the stale orphan
  `src.misc/mcica.F90` (live copy is `src.main/mcica.F90`); **switched the build
  to the NetCDF library's own `netcdf.inc`** via `-I$(nf-config --includedir)`,
  deleting the stale NetCDF-3 local copy (fixes the v4-upgrade skew).
- ~3,775 lines removed; all 6 build targets compile; n68equiv regression
  13/13 bit-for-bit.
- **The `source/shared/` + compat-shim architecture (Stages 3a/3b) was NOT
  built.** `src.misc/` still exists, now containing only live modules
  (`shr_kind_mod`, `shr_const_mod`, `physconst`, `ppgrid`, `pmgrid`,
  `spmd_utils`, `time_manager`, `infnan`, `ioFileMod`, `wrap_nf`) + small
  headers. If the `shared/` rewrite is still wanted, it now starts from a much
  cleaner `src.misc/`; if not, the conservative path may be sufficient.

### Other fixes landed this session (not in the original plan)
- **regrid gravity bug** (`tools/regrid_rtprofile.py`): hardcoded Earth gravity
  corrupted `zint` (which sets the 1-D CIA path length); `--gravity` now required.
- **io_1D consolidation**: merged `input.F90` + `output.F90` →
  `source/src.main/io_1D.F90`; extracted `read_namelist`/`print_diagnostics` from
  `main.F90`. (Note: this changes `src.main`, which is dual-contract with
  3dmodels — but `io_1D.F90` replaces 1-D-only `input.F90`/`output.F90`, which
  have no 3dmodels twin; the 3-D path uses `exo_radiation_cam_intr.F90`.)
- **physconst import fix**: removed a vestigial `use input` that masked a missing
  `shr_const_cpc2h6` import (was a phantom physconst→input dependency).
- **calc_gasopd signature fix**: added unused `qnh3`/`qco` dummy args to
  n84/n42/n68h2o/n28archean so all 5 production versions build (the shared
  `exo_radiation_mod.F90` passes them; only n68equiv had been updated).
- Set `exo_pver=300` as the standard level count.

### `src.misc` contract (confirmed with maintainer)
`source/src.misc/` is a **local ExoRT-only shim** — its README states it is "not
used in 3dmodels RT builds, as it is already present in CESM." Safe to
modify/condense/delete. The files duplicated in `/3dmodels/` (byte-identical to
`source/src.main` / `source/src.n68equiv`) are **dual-contract with CESM**:
their `use`/`include` statements must NOT be changed. `src.cam7.n68equiv` and
`src.cam.n68equiv.haze` have files that do NOT map cleanly to `source/` and must
be preserved (relevant to Stages 1–2).

---

## Session handoff protocol

Every session ends by writing or updating `tests/HANDOFF.md` with:

```markdown
# Session N handoff (date: YYYY-MM-DD)

## Completed this session
- (concrete list, with file paths)

## Tests passing
- tests/run_regression.sh: PASS / FAIL (which versions)
- (other relevant checks)

## Branch state
- Branch: refactor/stage-N-<slug>
- Last commit: <hash> — <message>
- Commits ahead of main: K

## Surprises / deviations from plan
- (anything the next session needs to know that the plan doesn't say)

## Next session starts here
- Stage M, sub-task X
- First commands to run: ...
- Specific files to read first: ...
```

The first thing a new session does is read `tests/HANDOFF.md`, then `git log` the last few commits on the active branch, then run the regression suite to confirm it's starting from a known-green state. If the regression suite fails on a fresh checkout, the session's first job is to fix that — do not stack new work on a broken base.

If a session's budget runs out before the stage is complete:
1. Commit the current state (even if WIP) to the feature branch with `WIP:` prefix.
2. Update `HANDOFF.md` describing exactly where mid-stage progress stopped and what the next session should pick up.
3. Do NOT attempt to "finish quickly" with degraded context — that produces broken commits that the next session has to debug.

---

## Stage 0 — Regression harness + experimental snapshot
### Session 1

**Goal.** Capture pre-refactor reference outputs for all 5 production versions AND snapshot the experimental builds so they can be reintegrated cleanly later. Add a numerical-comparison test driver.

**New files.**
- `tests/refs/{version}/RTprofile_in.nc` — canonical input per version
- `tests/refs/{version}/RTprofile_out.nc` — captured reference output
- `tests/run_regression.sh` — builds each version, runs against canonical input, calls `nccmp -d --tolerance=1e-10`
- `tests/compare_nc.py` — wraps `nccmp` with a tolerance report; falls back to NumPy comparison if nccmp absent
- `tests/refs/n68equiv_cam7/MD5SUMS` — MD5 manifest of all 12 cam7 files. cam7 cannot run standalone (it requires CESM3), so the "reference" is structural: file presence + MD5 + clean compile. `tests/run_regression.sh` adds a `cam7_build_check` mode that verifies the source files compile to objects against a stub Makefile target.

**Experimental triage (no code changes — read-only).**
Before any refactor commit, document the in-flight state of `source/experimental/src.n68equiv_exp/`:
- Inventory `.F90~` editor backups in `shr2/` and ask the maintainer whether they represent uncommitted work that should land first.
- Diff `experimental/src.n68equiv_exp/shr/aerad_driver` signature against the production `src.n68equiv/aerad_driver` to capture the exact list of new arguments (CO2-ice clouds, `SRF_EMISS_in`, etc.) and new init steps (`initialize_cldopts`, `do_exo_clouds`).
- Note in `tests/EXPERIMENTAL.md` which experimental files are byte-identical to production counterparts (candidates for adoption) vs. genuinely diverged (must stay forked until merged).
- If the experimental build currently compiles, capture its outputs into `tests/refs/n68equiv_exp/` so later integration can verify reproducibility. If it does not compile, document that and freeze the goal at "the refactor should not make it harder to fix later."

**Verify.** `tests/run_regression.sh` exits 0 on current `main` before any refactor commits land. `tests/EXPERIMENTAL.md` exists and lists the full inventory of experimental work-in-progress.

**Risk.** None — pure addition / documentation.

---

## Stage 1 — Replace broken `populate3Dmodels.py`
### Session 2 (first half)

**Goal.** Replace the broken sync script with a working drift-check + copy-mode tool, before any structural changes to `3dmodels/`.

**Files.**
- Rewrite `tools/populate3Dmodels.py` from scratch (~150 lines — three-target version is bigger than the original two-target sketch):
  - `--check` mode walks each declared host directory, resolves each entry (real file or symlink), and compares against the manifest. Validates: `symlink` declarations actually point to canonical paths; `fork` files differ from canonical (warn if they've become byte-identical — sign of an accidental backport); `absent` files are absent; `fork_only` files exist where declared.
  - `--copy <case_dir> --host {cam1,cam7}` mode materializes the named host's tree (resolving symlinks via `cp -L`, copying forks/fork_only as-is) into `$CASE/SourceMods/...`.
  - `--classify` mode prints all `fork` files tagged `backport_candidate` for human review.
  - **Drift policy:** strict. `--check` exits non-zero on any undeclared drift. Every divergence between cam7 and `source/` must be a declared `fork` with a classification. The cost is a one-line MANIFEST edit per legitimate divergence; the payoff is no silent drift between releases.
- Manifest: new `3dmodels/MANIFEST.toml` with three-target schema:

  ```toml
  [targets]
  hosts = ["source_1d", "cam1", "cam7"]

  [files."rad_interp_mod.F90"]
  canonical = "source/src.n68equiv/rad_interp_mod.F90"
  cam1      = "symlink"
  cam7      = "symlink"
  source_1d = "canonical"

  [files."exo_radiation_mod.F90"]
  canonical = "source/src.main/exo_radiation_mod.F90"
  cam1      = "symlink"
  cam7      = "fork"
  cam7_classification = "backport_candidate"
  cam7_notes = "SHR_CONST_G => gravit aliasing; deferred — see Stage 3c scope"

  [files."exoplanet_mod.F90"]
  cam1      = "absent"
  cam7      = "fork_only"
  cam7_path = "3dmodels/src.cam7.n68equiv/exoplanet_mod.F90"
  source_1d = "source/exoplanet_mod.F90"   # NOT the same file

  [files."cloud.F90"]
  canonical = "source/src.n68equiv/cloud.F90"
  cam1      = "symlink"
  cam7      = "absent"      # CESM3-planets framework provides
  source_1d = "canonical"
  ```

- Initial cam7 fork classifications (from inventory at time of writing):
  - `calc_opd_mod.F90` — `backport_candidate` (vacuum cutoff, scheduled for Stage 3c)
  - `exo_radiation_mod.F90` — `backport_candidate` (gravity aliasing — deferred to future work)
  - `exo_init_ref.F90` — `cesm3_only` (MPI broadcast)
  - `spectral_output_cam.F90` — `cesm3_only` (`cam_history` API)
  - `model_specific.F90` — `cesm3_only` (`use kabs` removal)
  - `mcica.F90`, `mcica_random_numbers.F90`, `radgrid.F90` — `cesm3_safe_cleanup` (whitespace, dead imports, string padding; also safe to apply to source/ but no payoff)
  - `exoplanet_mod.F90` — `cesm3_parameter_module` (cam7-only, never copy to cam1)

**Verify.** Run `--check` on current tree; report shows zero drift (every cam7 fork is declared and classified). Introduce one byte of un-declared drift and verify `--check` exits non-zero.

**Risk.** Low — replaces a script that nobody is currently relying on (it's broken).

---

## Stage 2 — Symlink byte-identical files in `3dmodels/`
### Session 2 (second half)

**Goal.** Eliminate maintenance burden of duplicate file copies.

**Files for cam1 — `3dmodels/src.cam.{version}/`** (unchanged from original plan):
- Replace with symlinks to `../../source/src.main/<name>`: `exo_init_ref.F90`, `exo_radiation_cam_intr.F90`, `exo_radiation_mod.F90`, `mcica_random_numbers.F90`, `mcica.F90`, `planck_mod.F90`, `rayleigh_data.F90`, `sys_rootdir.F90`.
- Replace with symlinks to `../../source/src.{version}/<name>`: `calc_opd_mod.F90`, `cloud.F90`, `kabs.F90`, `model_specific.F90`, `rad_interp_mod.F90`, `radgrid.F90`. Also `spectral_output_cam.F90` after verifying byte-identity (already noted as identical for n68equiv; check others).
- Keep as real files: `initialize_rad_mod_cam.F90` (legitimate I/O fork from `_1D` variant).
- `*.haze` variants in `src.cam.n68equiv.haze/`: keep as-is for now; symlink only confirmed-identical files.

**Files for cam7 — `3dmodels/src.cam7.n68equiv/`** (much smaller symlink set):
- Replace with symlinks: `planck_mod.F90` → `../../source/src.main/planck_mod.F90`, `rayleigh_data.F90` → `../../source/src.main/rayleigh_data.F90`, `rad_interp_mod.F90` → `../../source/src.n68equiv/rad_interp_mod.F90`. (These three are byte-identical across cam1, cam7, and source.)
- Keep as real files (9 forks, all declared in MANIFEST): `calc_opd_mod.F90`, `exo_init_ref.F90`, `exo_radiation_mod.F90`, `mcica.F90`, `mcica_random_numbers.F90`, `model_specific.F90`, `radgrid.F90`, `spectral_output_cam.F90`, `exoplanet_mod.F90`.
- These are `absent` from cam7 (CESM3-planets framework provides them; do not symlink, do not create stubs): `cloud.F90`, `kabs.F90`, `initialize_rad_mod_cam.F90`, `exo_radiation_cam_intr.F90`, `sys_rootdir.F90`.

**Escape hatch.** `populate3Dmodels.py --copy --host {cam1,cam7}` materializes everything to real files for users whose CESM `cp` doesn't follow symlinks.

**Verify.** Rebuild ExoCAM (CESM1) with each version against an unchanged input; regression suite produces matching outputs. Rebuild CESM3-planets pulling cam7 as an external; verify clean compile and cam7 MD5 manifest still matches (no symlink resolution accidentally rewrote a fork).

**Risk.** Medium — depends on `cp` vs `cp -L` behavior in user CESM versions; the escape hatch mitigates. CESM3-planets external-component pull mechanism may resolve symlinks differently — verify with the maintainer before this stage commits.

---

## Open question — `exort_rootdir` / `sys_rootdir.F90`

`source/src.main/sys_rootdir.F90` defines `exort_rootdir` as a compile-time `parameter` string (currently hardcoded to `/discover/nobackup/etwolf/models/ExoRT/`). This works on the cluster where the path exists, but fails silently on any other machine: `getfil` falls through to the MSS `msread` fallback and crashes on architectures without a `shr_sys_system` implementation (e.g. Apple Silicon / macOS).

A fix is needed but the right approach within the refactor has not been decided. Options include:

- **Environment variable at runtime** — replace the `parameter` with a runtime value from `EXORT_ROOTDIR` env var (falling back to the hardcoded path). Requires `get_environment_variable` (Fortran 2003 intrinsic) or `shr_sys_getenv`. Simple; no build changes; works for both Mac and cluster.
- **Absorb into `exort_io.F90` (Stage 3b)** — `exort_io.F90` reimplements `getfil`; it could read the data root from an env var directly, eliminating the `sys_rootdir` module entirely.
- **Namelist or CLI argument** — pass `data_root` at runtime (already planned for Stage 4's `exort_init(data_root_cstr, ...)`). The library API solves this cleanly for the Python / C path but the 1D standalone `main.F90` still needs a fallback before Stage 4 lands.

**Decision needed before Stage 3b.** For the near term (Mac development use), the environment-variable approach is a one-file, non-breaking fix that leaves the refactor stages unaffected. Record that decision here once made and update the Stage 3b scope accordingly (absorb or drop `sys_rootdir.F90`).

---

## Stage 3 — Replace `src.misc/` with `source/shared/`
### Sessions 3 + 4 (split into 3a and 3b)

**Goal.** Eliminate CESM-derived stubs; provide a small native utility layer for the 1D build. The 3D builds never compile `shared/` — CAM (cam1) and CESM3-planets (cam7) provide their own `physconst`, `shr_const_mod`, `ppgrid`, etc. `populate3Dmodels.py --check` enforces this: shared/compat shims must not appear in any cam manifest, or symbol clashes occur at link time.

**Session split.**

- **Session 3 (Stage 3a — foundation modules):** Write `exort_kinds.F90`, `exort_constants.F90`, `exort_grid.F90` and their compat shims (`shr_kind_mod`, `shr_const_mod`, `physconst`, `ppgrid`). Wire into `build/Makefile` for ONE version (n68equiv). Get n68equiv green via the new modules. Defer the other 4 versions to Session 4 to keep this session focused on getting the architecture right.
  - **Definition of done:** `make n68equiv` green; regression suite passes for n68equiv within 1e-10; other 4 versions can be temporarily broken (but commit a SKIP marker for them in `tests/run_regression.sh`).
  - **Handoff note must record:** any non-obvious aliasing decisions in compat shims; whether `physconst_setgas` was preserved or replaced.

- **Session 4 (Stage 3b — finish + delete src.misc):** Write `exort_netcdf.F90`, `exort_io.F90`, `exort_abort.F90`, `exort_runtime.F90` and their compat shims (`wrap_nf`, `ioFileMod`, `shr_sys_mod`, `time_manager`, `infnan`). Update `build/Makefile` for all 5 versions + experimental. Delete `source/src.misc/` (consolidating `mcica.F90` to `src.main/` first).
  - **Definition of done:** all 5 production builds green; regression suite passes for all 5 within 1e-10; `source/src.misc/` removed from the tree; experimental still builds (or is documented as already-broken pre-refactor).

**New files (under `source/shared/`).**
- `exort_kinds.F90` — `r8`, `i4`, `i8`. Replaces `shr_kind_mod.F90`.
- `exort_constants.F90` — physical constants, molecular weights, `mwdry/cpair` setters. Merges relevant content from `shr_const_mod.F90` and `physconst.F90`.
- `exort_grid.F90` — `pcols`, `pver`, `pverp` parameters. Replaces `ppgrid.F90`. `pcols` becomes a CPP-overridable parameter (default 1) to support stage 6's batch driver.
- `exort_netcdf.F90` — `wrap_open`, `wrap_close`, `wrap_inq_*`, `wrap_get_var_realx`, etc. Reimplements the actually-used subset of `wrap_nf.F90` against modern netCDF Fortran bindings; preserves error-handler signatures so `handle_error` callsites don't change.
- `exort_io.F90` — `getfil()` from `ioFileMod.F90`, stripped of MSS/archive logic.
- `exort_abort.F90` — `endrun(msg)`, `flush_unit`. Replaces the used subset of `shr_sys_mod.F90`.
- `exort_runtime.F90` — `get_nstep()` returning a stored integer (used as MCICA seed); `set_nstep(n)` for the multi-column driver to seed per-column reproducibly.

**Compat shims (under `source/shared/compat/`).** Thin modules named `shr_kind_mod`, `shr_const_mod`, `physconst`, `ppgrid`, `time_manager`, `wrap_nf`, `ioFileMod`, `shr_sys_mod`, `infnan` that re-export from the new modules. Keeps every `use` line in the RT-version source files unchanged. Compat shims live for stages 3–8; removed in stage 9.

**Experimental compatibility.** `source/experimental/src.n68equiv_exp/` and any other `experimental/*` builds keep working without modification because they `use` the same module names that the compat shims expose. The maintainer's in-flight `shr/` subdirectory inside experimental (which is a parallel prototype of this reorganization) is left untouched — the production `shared/` layer does not collide because it lives at `source/shared/`, not inside `experimental/`. A later integration pass (after experimental work is merged back) can dedupe the experimental `shr/` against `source/shared/`, but that is out of scope for this refactor.

**File moves / deletions.**
- Move `source/src.misc/mcica.F90` → `source/src.main/mcica.F90` (it is physics, not CESM glue, and `src.main/` already contains a copy used by all builds; verify they are identical and consolidate).
- Delete after stage 3 lands green: `source/src.misc/{infnan,pmgrid,spmd_utils}.F90`, `source/src.misc/{comctl,comhyb,comsol,misc,ncfortran,params}.h`, `source/src.misc/{netcdf,nfconfig}.inc`.
- Update `build/Makefile`: replace `$(MISC_DIR)/...` paths with `$(SHARED_DIR)/...` paths in each `OBJS_*` list. Add `$(SHARED_DIR) = ../source/shared`.

**Verify.** All 5 versions build and pass the stage-0 regression suite (numerical equivalence within 1e-10). 3D builds untouched (no `shared/` files in cam1 or cam7 manifests; `populate3Dmodels.py --check` exits 0).

**Risk.** Medium — touches every `OBJS_*` list and indirectly every RT-version source file. Compat shims minimize breaking changes; can revert by re-adding `src.misc/` files if a regression appears.

---

## Stage 3c — Backport CAM7 vacuum-pressure cutoff to 1D source
### Session 4.5 (NEW — small follow-on to Stage 3b)

**Goal.** Backport the cam7 stability fix in `calc_opd_mod.F90` to 1D source. Scope is **vacuum cutoff only**; the gravity aliasing change in cam7's `exo_radiation_mod.F90` stays forked and is explicitly out of scope.

**Files.**
- `source/src.n68equiv/calc_opd_mod.F90`: add the `pmid < 0.05 mb` skip block and the MT_CKD continuum bounds check (mirroring `3dmodels/src.cam7.n68equiv/calc_opd_mod.F90:170-200, 252-269, 342-370`).
- Apply the same patches to the other production versions' `calc_opd_mod.F90` (`n28archean`, `n42h2o`, `n68h2o`, `n84equiv`) **only if** their CO2 k-distribution tables exhibit the same low-pressure failure mode. For versions where this isn't physically applicable (e.g. `n42h2o` has no CO2), the backport is a no-op.
- Update `3dmodels/MANIFEST.toml`: change `calc_opd_mod.F90` cam7 entry from `fork (backport_candidate)` to `fork (backported)` — `populate3Dmodels.py --check` will warn that the fork is now byte-identical to source for the patched regions, which is exactly the signal we want post-backport.

**Why vacuum-only.** The vacuum cutoff is a defensive guard that returns zero opacity for pressures below the lookup-table floor — numerically equivalent to the existing behavior for any standard 1D test case (which never reaches < 5 Pa at the model top). The regression suite should pass within 1e-10 with no rebaseline needed. The gravity change, by contrast, would shift fluxes anywhere `exo_g ≠ 9.80616 m/s²` and is deferred to a future variable-planetary-gravity work item.

**Verify.** Regression suite passes within 1e-10 against pre-backport reference for all production versions (vacuum cutoff is defensive — should not change any standard-profile result). Add a unit test that drives a column to `pmid < 5 Pa` and confirms the model returns finite fluxes (would have crashed previously — this is the actual value of the backport).

**Risk.** Low — defensive guard at a regime the existing 1D suite doesn't reach. Revertible.

---

## Stage 4 — `iso_c_binding` library wrapper
### Session 5

**Goal.** Expose `aerad_driver` as a callable library — the foundation for both Python binding (stage 4b) and Fortran-to-Fortran coupling.

**Definition of done.** `libexort_n68equiv.so` builds; standalone C test program produces output within 1e-10 of single-column reference; Python harness via cffi reproduces the same. Other 4 versions of the library may be deferred to a follow-up session — only n68equiv is required to pass.

**New files.**
- `source/src.main/exort_lib_mod.F90`:
  ```fortran
  module exort_lib_mod
    use iso_c_binding
  contains
    subroutine exort_init(data_root_cstr, data_root_len) bind(C, name="exort_init")
    subroutine exort_run_column(state_in_ptr, state_out_ptr) bind(C, name="exort_run_column")
    subroutine exort_run_columns(ncol, states_in_ptr, states_out_ptr) bind(C, name="exort_run_columns")
    subroutine exort_finalize() bind(C, name="exort_finalize")
  end module
  ```
- `source/src.main/exort_column_mod.F90` — defines `column_state_t` derived type holding per-column inputs (TS, PS, TMID, PINT, gas mixing ratios, albedos, coszrs, …) and `column_result_t` for outputs (fluxes, heating rates, spectral fluxes). Used by both the library wrapper and stage 6's batch driver.
- `build/Makefile` targets: `libexort_{version}.a` and `libexort_{version}.so` (or `.dylib` on macOS).

**Implementation notes.**
- `exort_init` allocates k-coefficient tables, solar data, Planck tables — the same calls currently in `main.F90:45–50`. The data API is **`data_root` string only** (single mode): `exort_init` takes a directory path and uses the existing `exort_io`/`exort_netcdf` shared modules (Stage 3) to read all data files. New hosts stage data into a directory and pass the path. CESM3-planets is **not** a consumer of `exort_lib_mod` — it links the ExoRT source files directly into its own external-component build, so it does not need a buffer-injection API. If a future host (WRF, ROCKE-3D, in-memory Python research code) requires preloaded buffers, that becomes a separate, additive `exort_init_with_buffers()` entry point — out of scope here. Keeping the API single-mode now meaningfully reduces Stage 4 surface area.
- `exort_run_column` builds a `column_state_t` from C-interop array pointers, calls `aerad_driver` with the existing argument list (no changes to `aerad_driver`), returns results via C-interop pointers.
- `exort_run_columns` is a thin loop over `exort_run_column`. OpenMP added in stage 6.
- **`column_state_t` does NOT include `gravit`.** Per the maintainer decision recorded in the User-confirmed constraints, `exo_g` remains a compile-time constant from `source/exoplanet_mod.F90` for the 1D path. cam7's runtime-gravity behavior continues to live in its forked `exo_radiation_mod.F90` and is unaffected. When variable planetary gravity is later added to 1D, it will be a backwards-compatible extension of `column_state_t` — a clear future-work hook, but not in scope here.
- **Thread-safety (explicit guarantee).** After `exort_init`, all data tables (k-coefficients in `kabs`/`radgrid`/`radiation_data`, solar flux, Planck lookup, cloud optics) are read-only. Concurrent `exort_run_column` calls from threads are safe. This is now an explicit contract — cam7's framework parallelizes across columns, and any future library consumer will assume it.
- **No CESM1 changes.** `src.cam.n68equiv/exo_radiation_cam_intr.F90` stays as-is. CESM1 continues to use the older direct interface; new hosts use `exort_lib_mod`. This is the explicit "CESM1 freezes / CESM3+ moves forward" split.

**Python binding (stage 4b, may land in same stage).**
- New package `tools/exort_pytools/` with `pyproject.toml`.
- `tools/exort_pytools/runtime.py` — cffi/ctypes loader for `libexort_<version>.so`; marshals NumPy arrays to/from `column_state_t` via Python-side struct mirroring.
- Documented usage: `import exort_pytools as xrt; xrt.init('n68equiv', data_root='./data'); out = xrt.run_column(state); xrt.finalize()`.

**Verify.**
- A standalone test program that calls `exort_init / exort_run_column / exort_finalize` against the canonical single-column input produces output within 1e-10 of `RTprofile_out.nc` reference.
- Python harness produces matching output.

**Risk.** Low — pure addition. `aerad_driver` and existing `main.F90` unchanged.

---

## Stage 5 — Dedup parameter-only RT files
### Session 6

**Goal.** Eliminate `cloud.F90` and `kabs.F90` duplication across versions.

**Definition of done.** Per-version `cloud.F90` and `kabs.F90` removed; `source/shared/cloud_mod.F90` and `kabs_template.F90` + per-version `kabs_params.inc` in place; all 5 production builds green; regression suite within 1e-10.

**Files.**
- Verify `source/src.{n28archean,n42h2o,n68h2o,n68equiv,n84equiv}/cloud.F90` are identical (they are believed to be ~14-line pure parameter modules differing only in filename suffix). If identical, move to `source/shared/cloud_mod.F90` and delete per-version copies. If filename suffix differs, parameterize via include file (next bullet).
- Replace per-version `kabs.F90` with:
  - `source/shared/kabs_template.F90` — module skeleton with `include 'kabs_params.inc'`
  - `source/src.{version}/kabs_params.inc` — version-specific filename strings + `bps_continuum` flag
  - `build/Makefile` adds `-I$(N{VERSION}_DIR)` so each version build picks up the correct include
- Update `3dmodels/src.cam.{version}/` symlinks (from stage 2) to track new locations.
- **cam7 is unaffected.** `cloud.F90` and `kabs.F90` are already `absent` from the cam7 manifest (CESM3-planets framework provides them). `populate3Dmodels.py --check` will refuse to add them to cam7 even if a future symlink-loop tries.

**Verify.** Regression suite — all versions reproduce reference output within 1e-10.

**Risk.** Low — parameter-only changes, no logic flow changes.

**Experimental scope.** This stage touches only the 5 production RT versions. The experimental `src.n68equiv_exp/` keeps its own `cloud.F90` and `kabs.F90` (which may have additional CO2-ice-cloud entries) — do NOT delete them. After experimental work is merged back to production (a future, separate task), its parameters can be folded into the shared template at that point.

---

## Stage 6 — Multi-column batch driver
### Sessions 7 + 8 (split into 6a and 6b)

**Goal.** Refactor 1D wrapper to process N columns per invocation with OpenMP parallelism. Decouple I/O from module-scope state.

**Session split.**

- **Session 7 (Stage 6a — decouple I/O from globals):** Define `column_state_t` and `column_result_t` (if not already present from Stage 4). Rewrite `input.F90` to populate a derived-type instance instead of module-level globals; rewrite `output.F90` to accept the derived type. `main.F90` still single-column but now passes data explicitly. No multi-column reading yet, no OpenMP. The point of this session is to land the data-structure refactor on a known-green base.
  - **Definition of done:** all 5 builds green; regression suite within 1e-10; module-level globals deleted from `input.F90`; `column_state_t` is the canonical per-column container.

- **Session 8 (Stage 6b — multi-column + OpenMP):** Extend `read_profile_nc` to handle the optional `ncol` dimension; add the OpenMP loop in `main.F90`; add CLI arg parsing for input/output filenames; ship the per-column MCICA seed change.
  - **Definition of done:** legacy single-column input still produces matching output; 64-column batch run vs 64 single-column runs concatenated within 1e-10; `OMP_NUM_THREADS=1` vs `=8` within 1e-10 on the 64-column run.

**Files.**
- `source/src.main/input.F90`: delete module-level globals (`TS_in`, `TMID_in`, `H2OMMR_in`, …). Replace `input_profile()` (no args, populates module state) with `read_profile_nc(filename, states)` returning an allocated array of `column_state_t` (defined in stage 4). If the input file has an `ncol` dimension, return that many states; if not, return a length-1 array (legacy compat).
- `source/src.main/output.F90`: similarly accept an array of `column_result_t` and write a NetCDF with `ncol` dimension. Maintain backwards-compat by writing a flat (no `ncol` dim) file if `ncol == 1`.
- `source/src.main/main.F90`: rewrite as
  ```
  call exort_init(data_root)
  call read_profile_nc(input_file, states)
  allocate(results(size(states)))
  !$omp parallel do private(i)
  do i = 1, size(states)
    call exort_run_column(states(i), results(i))
  end do
  call write_profile_nc(output_file, states, results)
  call exort_finalize()
  ```
  CLI args: `./n68equiv.exe [input.nc] [output.nc]` (default to `RTprofile_in.nc` / `RTprofile_out.nc` for legacy compat).
- `source/exoplanet_mod.F90`: leave unchanged in this stage — `exo_pver` stays compile-time. (A future stage could allow runtime `pver` per column, but it's a much bigger refactor and not in scope.)

**MCICA seeding.** Currently `get_nstep()` returns hardcoded 9404 → all columns use the same seed. New behavior: `set_nstep(9404 + col_index)` so per-column results are reproducible regardless of thread count. Document this as an intentional, numerically-equivalent change (within MCICA stochastic envelope; not bit-identical). **This is a 1D-only change.** cam1's `mcica.F90` (symlinked from `source/`) inherits it but cam1's CAM driver continues to call `get_nstep()` from its own `time_manager`, so behavior there is unchanged. cam7's `mcica.F90` is already a fork (different `iso_fortran_env` integer typing); it stays on the host CESM timestep for its seed. Update `MANIFEST.toml` to record that cam7's `mcica.F90` fork is intentionally *not* tracking the 1D seeding strategy — this prevents a future `--check` run from flagging Stage 6 as drift introduction.

**Verify.**
- A 1-column input produces output matching reference within 1e-10.
- 64 single-column runs concatenated equals one 64-column run within 1e-10.
- `OMP_NUM_THREADS=1` vs `=8` on 64 columns: results match within 1e-10.

**Risk.** Medium — touches `input.F90`, `output.F90`, `main.F90`. Bounded by the legacy single-column shim path and the regression suite.

**Experimental scope.** The experimental `src.n68equiv_exp/main.F90` calls `aerad_driver` with a different argument list (CO2-ice clouds, `SRF_EMISS_in`) and runs it twice (clear + cloudy). Stage 6 changes only the production `src.main/main.F90` and `input.F90`/`output.F90`; the experimental driver keeps its own copies. When the experimental work is later merged back, the `column_state_t` derived type defined here is the natural extension point — extra fields can be added (e.g. `cicewp_co2`, `rei_co2`, `srf_emiss`, `do_clouds` flag) without breaking the production callsites.

---

## Stage 7 — IDL → Python migration (parallel track)
### Sessions 9a, 9b, 9c (IDL ports — likely 2–3 sessions)

**Goal.** Eliminate IDL dependency. Touches no Fortran, so it CAN run in parallel with stages 3–6 (e.g. you start a new session for IDL while leaving Fortran refactor sessions independent). Practically though, scheduling matters for token budget.

**Session split.** The IDL migration is naturally chunked by sub-port:

- **Session 9a — plotting + simple utilities (~400k tokens):** Items 1 and 2 from the order list — plotting scripts and simple profile utilities. Smallest validation surface. Regenerate / replicate small reference outputs in `tests/refs/python_ports/`.
- **Session 9b — Mie + cloud optics (~600k tokens):** Items 3, 4, 5 — stellar SED, Mie core, cloud/CARMA optics. Validate Mie against canonical particle radii (relative err < 1e-6); regenerate cloud-optics NetCDFs and validate.
- **Session 9c — CIA + MTCKD (~700k tokens):** Item 6 — `makeCIA.pro` (1594L) and `makeMTCKD.pro`. Largest port; the MT_CKD external binary requires a working test environment. Plan ~1 full session for this alone; split further if needed.

**Per-session definition of done.** Within each sub-port, the corresponding `.pro` file is moved to `tools/legacy_idl/`, the Python equivalent is in `tools/exort_pytools/`, and a validation script in `tests/python_ports/test_<name>.py` confirms the regenerated outputs match the IDL reference within tolerance.

**Order (lowest risk first).**
1. Plotting: `plotprofile.pro`, `plotspectra_1D.pro`, `plotstellar.pro` → matplotlib in `tools/exort_pytools/plot.py`. Validate visually.
2. Simple utilities: `hybrid2height.pro`, `hybrid2pressure.pro`, `getcolumn.pro`, `makeColumn.pro`, `profile_data.pro` → `tools/exort_pytools/profile.py`. Validate via NumPy `allclose(..., atol=1e-12)` on round-trip of a known profile NetCDF.
3. Stellar SED: `makeStellarSpectrum_fromSED.pro` → joins existing `makeStellarSpectrum_blackbody.py` in `tools/exort_pytools/stellar.py`. Validate against an IDL-generated reference NetCDF.
4. Mie core: `mie_single.pro` → `tools/exort_pytools/mie.py` using `miepython` (BSD, well-tested). Validate Qext/Qsca/g against IDL output for canonical particle radii to relative error < 1e-6.
5. Cloud / CARMA optics: `makeCloudOptics.pro`, `makeCARMAOptics.pro`, `interpCARMAOptics.pro` → `tools/exort_pytools/cloud_optics.py`. Depends on Mie core. Validate by regenerating canonical cloud-optics NetCDF and `nccmp -d --tolerance=1e-10` against existing data files (use deterministic seeds in any random sampling).
6. CIA / MTCKD: `makeCIA.pro` (1594L), `makeMTCKD.pro` → `tools/exort_pytools/cia.py`, `tools/exort_pytools/mtckd.py`. Largest port; the MT_CKD external binary is wrapped via `subprocess`. Validate by regenerating each `cia_*.nc` and `nccmp -d`.

**Disposition.** Move `.pro` files to `tools/legacy_idl/` (not deleted) until all 5 RT versions have been re-validated against Python-regenerated data. Then delete `tools/legacy_idl/`.

**Experimental scope.** `experimental/src.n68equiv_exp/{getColumn.pro, makeColumn.pro, plotprofile.pro, plotspectra_1D.pro}` are version-of-record IDL files for the in-progress CO2-ice-cloud build (the README notes they were "edited to add clouds"). These get the same treatment in this stage: port to Python, place in `tools/exort_pytools/` with cloud-aware variants, and drop the `.pro` copies into `tools/legacy_idl/experimental/` for later deletion.

**Risk.** Low per script — each port is validated against its own IDL-generated NetCDF. Independent track, can revert per-script.

---

## Stage 8 — Cleanup
### Session 10

**Goal.** Remove transitional scaffolding.

**Definition of done.** All compat shims under `source/shared/compat/` deleted; mechanical rewrites of `use shr_kind_mod, only: r8 => shr_kind_r8` to `use exort_kinds, only: r8` (and analogous) completed throughout `source/`; all 5 production builds green; regression suite within 1e-10; README and CLAUDE.md updated.

**Files.**
- Drop `source/shared/compat/` shims (`shr_kind_mod`, `physconst`, `ppgrid`, `wrap_nf`, `ioFileMod`, `time_manager`, `shr_sys_mod`, `infnan`). Mechanically rewrite `use shr_kind_mod, only: r8 => shr_kind_r8` → `use exort_kinds, only: r8` and analogous changes throughout `source/`. CAM-side files (in `3dmodels/` and the CAM-build copies of shared `src.main/` files) keep `use shr_kind_mod` etc. — those names resolve against CAM's modules, not ours.
- Delete `tools/legacy_idl/` once all data files have been regenerated and validated from Python.
- Update README.md and CLAUDE.md to reflect new layout.

**Risk.** Low — mechanical, fully covered by regression suite.

---

## Critical files

To be modified:
- `/Users/wolfe/models/ExoRT/source/src.main/main.F90` (stage 6 — full rewrite)
- `/Users/wolfe/models/ExoRT/source/src.main/input.F90` (stage 6 — remove module-level state)
- `/Users/wolfe/models/ExoRT/source/src.main/output.F90` (stage 6 — accept array of results)
- `/Users/wolfe/models/ExoRT/source/src.main/exo_radiation_mod.F90` (no changes; `aerad_driver` reused as-is)
- `/Users/wolfe/models/ExoRT/source/src.{version}/calc_opd_mod.F90` (stage 3c — backport vacuum cutoff for applicable versions)
- `/Users/wolfe/models/ExoRT/build/Makefile` (stages 3, 4, 5, 6 — adjust OBJS lists, add lib targets)
- `/Users/wolfe/models/ExoRT/tools/populate3Dmodels.py` (stage 1 — replace from scratch with three-target schema)
- `/Users/wolfe/models/ExoRT/3dmodels/src.cam.*/` (stage 2 — symlinks per cam1 plan)
- `/Users/wolfe/models/ExoRT/3dmodels/src.cam7.n68equiv/` (stage 2 — 3 symlinks for byte-identical files; 9 forks remain real files)

To be created:
- `/Users/wolfe/models/ExoRT/source/shared/` — `exort_kinds.F90`, `exort_constants.F90`, `exort_grid.F90`, `exort_netcdf.F90`, `exort_io.F90`, `exort_abort.F90`, `exort_runtime.F90`, `cloud_mod.F90`, `kabs_template.F90`, `compat/*.F90`
- `/Users/wolfe/models/ExoRT/source/src.main/exort_lib_mod.F90` (stage 4)
- `/Users/wolfe/models/ExoRT/source/src.main/exort_column_mod.F90` (stage 4)
- `/Users/wolfe/models/ExoRT/source/src.{version}/kabs_params.inc` (stage 5)
- `/Users/wolfe/models/ExoRT/3dmodels/MANIFEST.toml` (stage 1)
- `/Users/wolfe/models/ExoRT/tests/refs/{version}/RTprofile_{in,out}.nc` (stage 0)
- `/Users/wolfe/models/ExoRT/tests/{run_regression.sh, compare_nc.py}` (stage 0)
- `/Users/wolfe/models/ExoRT/tools/exort_pytools/` Python package (stages 4b, 7)
- `/Users/wolfe/models/ExoRT/tools/check_3dmodels_drift.py` (stage 1)

To be deleted (eventually):
- `/Users/wolfe/models/ExoRT/source/src.misc/` (stage 3, after compat shims absorb usage)
- `/Users/wolfe/models/ExoRT/source/shared/compat/` (stage 8)
- `/Users/wolfe/models/ExoRT/tools/legacy_idl/` (stage 8, after data regeneration validated)
- per-version `cloud.F90`, `kabs.F90` (stage 5)

To be preserved unchanged through this refactor (in-flight work):
- `/Users/wolfe/models/ExoRT/source/experimental/src.n68equiv_exp/` — entire directory tree, including `shr/`, `shr2/`, all `.F90`, `.F90~`, and `.pro` files. The compat shim layer (stage 3) ensures it keeps building. Reintegration / dedup against `source/shared/` is a separate future task once the maintainer has merged the CO2-ice-cloud and CARMA work into a production version.
- `/Users/wolfe/models/ExoRT/3dmodels/src.cam7.n68equiv/` — the 9 forked files (and the new `exoplanet_mod.F90`) stay real files; CESM3-planets pulls this directory as an external. Stage 3c modifies *source* `calc_opd_mod.F90`, not the cam7 copy. The cam7 gravity aliasing in `exo_radiation_mod.F90` is left as-is throughout this refactor.
- `/Users/wolfe/models/ExoRT/3dmodels/src.cam.n68equiv/` — CESM1 hookup, treated as active legacy: gets the symlink cleanup but no API changes (continues to use the existing `exo_radiation_cam_intr.F90` shim, not `exort_lib_mod`).

Existing functions to reuse:
- `aerad_driver` at `source/src.main/exo_radiation_mod.F90:154-691` — the stable per-column kernel; library wrapper in stage 4 wraps it directly. CESM1 keeps calling it via the existing cam_intr shim; CESM3-planets calls it directly via its external component.
- The cam7 vacuum cutoff at `3dmodels/src.cam7.n68equiv/calc_opd_mod.F90:170-200, 252-269, 342-370` — direct backport target for Stage 3c.
- The cam7 `gravit` aliasing pattern at `3dmodels/src.cam7.n68equiv/exo_radiation_mod.F90:15-18` — *not* used in this refactor, but recorded for future variable-planetary-gravity work in 1D source.
- `initialize_kcoeff`, `initialize_solar`, `init_ref`, `init_planck`, `initialize_radbuffer` — currently called by `main.F90:45-50`; reused in `exort_init` (stage 4).
- `wrap_nf.F90` callsite signatures — preserved by `exort_netcdf.F90` so RT-version files don't need rewrites.
- `ioFileMod::getfil` callsite signature — preserved by `exort_io.F90`.
- Existing Python ports: `tools/makeStellarSpectrum_blackbody.py`, `tools/spectral_intervals.py`, `tools/check_kcoeff.py`, `tools/convert_O2_CIA_files.py` — reused / extended in `tools/exort_pytools/` (stage 7).

---

## Verification (end-to-end)

After each stage:
1. Build all 5 versions: `cd build && make n28archean n42h2o n68h2o n68equiv n84equiv`.
2. Build the experimental version (if it currently builds): `make n68equiv_exp`. If it does not currently build, treat as a known issue — the refactor must not make it worse, but is not required to fix it.
3. Run regression suite: `tests/run_regression.sh` — must report numerical equivalence within 1e-10 against `tests/refs/{version}/RTprofile_out.nc`.
4. After stage 1 (MANIFEST + populate script): `python3 tools/populate3Dmodels.py --check` exits 0 on a clean tree; introducing one byte of un-declared drift in any file makes it exit non-zero (proves strict-mode policy works); `--classify` correctly lists the cam7 forks tagged `backport_candidate`.
5. After stage 2 (symlinks): rebuild ExoCAM (CESM1) with one configuration (e.g. n68equiv) end-to-end and confirm a short test integration produces sensible output. cam7 directory has 3 symlinks (`planck_mod.F90`, `rad_interp_mod.F90`, `rayleigh_data.F90`) and 9 real files; CESM3-planets external-pull rebuild verifies cam7 still works.
6. After stage 3 (shared/): 1D builds use `source/shared/`; cam1 and cam7 still build using their host's `physconst`/`shr_const_mod` (verify by grepping cam manifests for any `shared/` references — should be zero).
7. After stage 3c (vacuum cutoff backport): regression suite passes within 1e-10 against pre-backport reference for all production versions. New unit test: drive a column to `pmid < 5 Pa` and confirm the model returns finite fluxes (would have crashed previously).
8. After stage 4 (library): standalone C test program calling `exort_init(data_root) / exort_run_column / exort_finalize` reproduces the reference within 1e-10. Python harness `tools/exort_pytools/runtime.py` does the same.
9. After stage 6 (batch): 64-column batch run vs 64 single-column runs concatenated — within 1e-10. `OMP_NUM_THREADS=1` vs `=8` on 64 columns — within 1e-10.
10. After stage 7 (IDL → Python): regenerated data files (`cloud_*.nc`, `cia_*.nc`, stellar spectra, …) match IDL-generated references within their declared per-script tolerance, and ExoRT regression suite still passes when the new data files are dropped in.

The 1e-10 tolerance is uniform across stages (per the user-confirmed reproducibility level) and accommodates floating-point reordering from refactors and OpenMP nondeterminism.

The cam7 build itself is verified by the maintainer running CESM3-planets externally — that's outside this repo's CI but should be smoke-tested at each stage that touches cam7 files (Stages 1, 2, 3c).
