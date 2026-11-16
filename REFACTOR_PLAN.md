# ExoRT v2 — plan: what is left

Forward-looking only. What has already landed is summarised in `CHANGES.md`;
the detail is in `git log`. Rewritten 2026-09-23; stages A–E, the Stage C
physics, the 3-D port code and the HITRAN-2024 switch are all done.

## Ground rules (apply to every item)

- **Branches.** `main` is the v1 maintenance line (`v1.0.0`, `v1.1.0` = H₂O
  T-index fix). All v2 work stays on `refactor` until the maintainer decides to
  merge.
- **Gate every change on the regression suite.** A change that should not move
  physics must give Δ = 0 on all cases. An intended physics change regenerates
  the baselines in the same commit and says so in the message, then refreshes
  `tests/regression/REGRESSION_STATUS.md`.
- **Keep the CAM bundle in sync.** `3dmodels/src.cam.exort` stays
  byte-identical to `source/`: edit `source/`, then
  `populate3Dmodels.py regenerate --exort`, then `check`, then
  `tests/cam_compile_check/run_all.sh`. The pre-v2 bundles
  (`src.cam.n68equiv`, `.n68equiv.haze`, `.n84equiv`, `src.cam7.n68equiv`) are
  frozen; never regenerate them.
- **Keep changes small.** One purpose per commit, so `git revert` is always a
  clean undo.

---

## 1. Validation (blocks a v2 release)

1. **First real ExoCAM run of `src.cam.exort`** (HPC). Code-complete and
   compile-checked only.
   - **Setup:** add `do_exo_clouds`/`do_exo_haze` to ExoCAM's
     `src.share/exoplanet_mod.F90` (template:
     `tests/cam_compile_check/exoplanet_mod_stub.F90`); point `exo_solar_file`
     at `data/stellar/*_n84.nc`; use an ExoRT checkout whose `data/kdist` is
     current (the k-file grid check will `endrun` on old O₂/O₃ files).
   - **Order:** a clear-sky aquaplanet against an identical `src.cam.n68equiv`
     case, expecting LW close and SW within ~0.4 % (UV regrid) plus the H₂O
     T-fix and HITRAN-2024 deltas in `CHANGES.md`. Then the gated physics per
     config: `-DEXORT_CO2CLD`, `-DEXORT_CARMA`, `-DEXORT_SRF_EMISS`.
   - This first run is also the first real-PIO test of the k-file grid check.
   - **PR #7 (per-thread solver work arrays, merged 2026-11-16)** has not
     yet run in CAM; the same run also checks it:
     - Build with Intel and OpenMP (`NTHRDS_ATM > 1`). The compile alone
       confirms that the CESM compiler accepts `!$omp threadprivate` on
       SAVEd allocatable locals (OpenMP 3.0). Only ifx and gfortran have
       been tried, through the stub harness.
     - Run NTHRDS 1 vs 2 with the same executable: `cam.h0` must be
       identical. A PR build vs a parent build may differ at roundoff
       (about 1e-9 relative at `-O2`; see §3 item 10), which is not a
       regression.
     - Compare the radiation timer in the CESM timing file against the
       parent build. The 2.4× gain was measured for ExoColumn on
       glibc; CAM's heap may never hit the trim pattern, so the gain there
       may be close to zero.
     - Memory: each thread now keeps about 20–30 MB at CAM level counts for
       the whole run (e.g. 28 tasks per node ≈ 0.6–0.8 GB per node, held,
       not transient).
2. **Drop the three `EXORT_*` CPP macros** once item 1 passes. The deferred
   plan: `do_exo_condense_co2` everywhere, `srf_emiss` standard in
   `camsrfexch`, and CARMA gated by its own build directive.
3. **ExoCAM impact of the H₂O T-index fix.** The shifted table was in every
   n68equiv CAM run since 2020-11-09 (`da04522`). Which published runs need a
   note or a rerun is the maintainer's call.
4. **v1 legacy checkpoint** (maintainer action, on `main`): a final
   run-and-save of `n28archean`/`n42h2o`/`n68h2o` against their published
   outputs. Not a blocker for anything on `refactor`.

## 2. Physics and data

1. **CO edge bands 29 (1950–2050 cm⁻¹) and 31 (2200–2397 cm⁻¹).** A
   hand-rolled line-by-line check disagreed at every temperature offset.
   Suspected: line-wing or cutoff handling. The n68 CO table is `subL`, the
   n84 one `voigt`. **Unconfirmed:** re-check with a trusted LBL code before
   re-fitting. Low impact: 1 % CO moves OLR by only −1.1 W/m² in total. A
   re-fit will intentionally fail the four `CO_*` regression cases; rebaseline
   them on purpose.
2. **Fractal haze optics above band 68.** `haze_n84_b40_fractal_interp.nc`
   UV bands are a nearest-band copy of band 68. Regenerating it needs
   `fractaloptics.exe` on Discover. Mie stays the default until then.
3. **k-table grid expansion** (science-driven): P to 100 bar (71 levels), T
   to 3000 K+ (68 levels, or fewer if non-uniform). This means new HELIOS-K
   tables and a rebaseline. Land `OPTIMIZATION.md` Fix 1 (binary search) and
   Fix 2 (free duplicate k arrays) first, so the search change isn't
   entangled with the rebaseline.
4. **Band edges in k-files.** `SpectralBands` holds only 1…N. Have the
   `heliosk2netcdf` converter write the wavenumber edges (e.g. `BandEdges`,
   N+1 values); then `check_kfile_grid` can compare them against
   `wavenum_edge`.

## 3. Code

1. **Clear-sky / cloud-forcing `_CLD` double run** (1-D). The reference
   implementation is at `0e409c3` under `source/experimental/src.n68equiv_exp/`
   (`main.F90` calls `aerad_driver` twice; `output.F90` writes `*_CLD`).
   Convention: in 1-D, clear sky is the default and cloudy is the extra — the
   reverse of the 3-D model.
2. **Robustness back-ports from `src.cam7.n68equiv`** — see `CESM.CAM7.md`
   Tier 1: the model-top vacuum cutoff (`pmid < 0.05 mb`), out-of-range-T
   guards on MT_CKD and every CIA pair, and the MT_CKD dry/cold skip.
   Expected Δ = 0 except at the very model top.
3. **Performance** — `OPTIMIZATION.md`. Profile first; then the Δ = 0 wins
   (skip zero gases, hoist interpolation factors, binary-search CIA T);
   then evaluate the k-array layout change.
4. **Input checks for the other tables.** The MT_CKD continuum and CIA
   files load with no grid or dimension checks. Extend the
   `check_kfile_grid` idea to them.
5. **Fatal-error exit codes.** Only the k-file checker uses `stop 1`. The
   other loaders use bare `stop` (exit 0), so scripts see success. Convert
   them all.
6. **Small dead-code cleanups:** `exo_radiation_mod` module-scope
   `openstatus` (unused; the local of the same name is the live one);
   `radgrid` `ntopcld`/`nlevsRT` (never read or written); `physconst`
   `gravit`/`rga`/`cpvir`/`zvir` (unused on the 1-D path).
7. **Retire `n68equiv`/`n84equiv`** — `source/`, the Makefile and
   `populate3Dmodels` mappings. `gas_sweep.py` uses them as HITRAN-2016
   references, so decide first whether that diagnostic is still needed or
   can compare `exort --exort h16` against `exort` instead.
8. **I/O naming normalisation** (v2-only breaking change):
   `RTprofile_in.nc`/`RTprofile_out.nc` → consistently cased, symmetric
   names. Do it in one commit, with the regression suite proving equivalence.
9. **GPU (distant, "E3").** Only if a real many-column workload lands on
   NVIDIA hardware: a directive-based (OpenACC or OpenMP `target`) port of
   the batched kernel. Nothing in the current design blocks it.
10. **Solver work-array follow-ups to PR #7** (merged 2026-11-16). #7 made
    the 81 solver work arrays in `exo_radiation_mod`
    `allocatable, save` + `!$omp threadprivate`, allocated once per thread:
    1. **macOS slowdown: pass the scratch arrays as plain arguments.**
       - **Problem.** Homebrew gfortran on macOS has no native thread-local
         storage and emulates it. Every reference to a threadprivate array
         calls `___emutls_get_address` to find that thread's copy. Measured
         on this Mac (88-column batch): serial 7.0 → 8.2 s at `-O0`
         (+17 %), 2.33 → 2.52 s at `-O2` (+8 %); 4 threads 2.29 → 2.66 s.
         The Mac never had the glibc page-fault problem #7 fixes, so there it
         is cost only. Linux uses native TLS and is not affected.
       - **Where the cost is.** Arrays reached through a dummy argument
         carry no TLS lookup; the callee just sees an address. The
         lookups happen where a threadprivate array is referenced by name.
         `aerad_driver`'s own 60 arrays are named only in the allocate
         guard, the whole-array zero-fill and the call sites, so they cost
         about one lookup per array per call. The hot cost is the 21
         threadprivate locals of the kernels, used element by element in
         the inner loops:
         - `add_txrad`: `CPP CM EE3 EL3 AS DF DS XK`
         - `refine_lwflux`: `uintent dintent Y1 Y2 Y4 Y8 Y5 A1 A2 A3 A4 A7`
         - `dplanck`: `ptemp_midlayer`
       - **Change.** Hoist those 21 arrays into `aerad_driver`'s threadprivate
         set: same shapes, same allocate-once guard. Remove `save` and
         `threadprivate` from the kernels, and give each kernel the arrays as
         extra explicit-shape dummy arguments,
         e.g. `real(r8), intent(inout) :: CPP(ntot_gpt,pverp)`. Every kernel
         is `private` and called only from `aerad_driver` (`dplanck` via
         `rad_precalc`), so no public interface changes; `rad_precalc`
         gains the `ptemp_midlayer` pass-through. `add_txrad` is called
         twice per column (solar, then IR), one after the other, so both
         calls can share one set of scratch arrays.
       - **Equivalent alternative.** Make `aerad_driver` a thin wrapper that
         allocates and then calls an inner `aerad_solve` with every array as
         an argument. It is a bigger diff for no extra gain, since the
         driver's own arrays are already lookup-light.
       - **Verify.**
         - Regression Δ = 0, `multicol_check`, `percol_config_check`,
           `verify_lib`, sync + CAM compile gate.
         - `nm`/`otool` on the macOS object: no `___emutls_get_address`
           left in the kernels.
         - Mac timing back to the pre-#7 numbers; Linux timing unchanged.
    2. **Release routine for long-lived processes.** Each thread holds its
       arrays until it exits, about 165 MB per thread at pver = 300 (16
       threads: peak 2.5 → 3.1 GB). Add `exort_release_workspace` (a
       parallel region that deallocates each thread's copies) to the
       library API, for Python/`libexort` callers that want the memory back.
    3. **Drop redundant zero-fills** in `aerad_driver`. After #7 the
       per-call zero-fill is about 16 % of ExoColumn's time.
       - Redundant: `TAUL`/`OPD` (re-zeroed in `rad_precalc`), `tau_gas`
         (in `calc_opd_gas`), `Y3` (fully overwritten); the `*ir` arrays
         are read only over the LW range.
       - Must stay: the `*sol` arrays (night columns; `rad_postcalc` reads
         them for `vis/nir_dif`) and the cloud/haze optics (clear-sky
         companion call).
       - Separate commit, Δ = 0.
    4. **Remaining per-call allocations**, same treatment:
       - `calc_opd_gas` (11 CIA arrays)
       - `calc_opd_cld_h2o` (9 arrays, ~2 MB at pver = 70)
       - `mcica_subcol` (2)

       These live in three bundle copies of `calc_opd_mod.F90` plus
       `mcica.F90`.
    5. **Bit-identity holds only at `-O0`.** #7 is byte-identical to its
       parent at the Makefile's `-O0`. At `-O2` (arm64) outputs differ at
       roundoff, heating rate up to 1e-9 relative, less than the
       `-O0 → -O2` difference itself. The likely cause is heap alignment of
       reused vs fresh arrays changing how loops are vectorized and summed.
       Expect the same in optimized ExoCAM builds.
    6. **Non-OpenMP builds are not thread-safe.** With
       `make OMPFLAGS=` the directive is a comment and the arrays are
       shared `SAVE` storage, as in v1. Documented in the `aerad_driver`
       header, `exort_lib_mod.F90` and `DEVELOPER_NOTES.md`; no default
       target is affected. CAM without OpenMP runs physics serially per
       task, so it is safe.
11. **Data root without editing a tracked file.** `exort_rootdir` is a
    hard-coded absolute path in `source/src.main/sys_rootdir.F90`, so every
    machine or checkout must edit a tracked file before building. The edit
    then shows as a change waiting to be committed, and a commit from one
    machine breaks the next (2026-11-16: the committed
    `/Users/wolfe/models/ExoRT-refactor/` made every case fail on a Mac whose
    home is `/Users/wolf`).
    - **Change.** Resolve the root at runtime and stop editing the file:
      1. the environment variable `EXORT_ROOT`, if set;
      2. else `exort_rootdir` in the `&exort_config` namelist of
         `run/user_nl_exort`, which is already gitignored and already read at
         startup;
      3. else a compiled relative default, `'../'`. `exort.exe` always runs
         from `run/`, so this works in any checkout with no configuration.

      Append a missing trailing slash, as `exort_init` already does, and print
      the resolved root in the `=== exort_config ===` header.
    - **Library.** `exort_init(data_root)` already overrides the root.
      Keep it; an empty `data_root` should fall back to the same order.
    - **CAM.** `3dmodels/src.cam.exort/sys_rootdir.F90` is a separate,
      unsynced file (a `parameter` set to the Discover path). Take the root
      from an ExoCAM namelist variable in `exoplanet_mod`, alongside
      `exo_solar_file`, so a CAM case needs no source edit either.
    - **Clean up.** Delete the per-machine comment block in
      `sys_rootdir.F90`. Update the `GETFIL` hint in `ioFileMod.F90`, which
      says "Check exort_rootdir in source/src.main/sys_rootdir.F90",
      `iofiles/user_nl_exort.template`, the CLAUDE.md run instructions, and
      the `src.cam.exort` README step 1.
    - **Verify.** Regression Δ = 0. Then run from a second checkout with
      `EXORT_ROOT` unset and no namelist entry: it must find its data through
      `'../'`.

## 4. Tooling and tests

1. **Parameter sweeps and training data.** The regression suite now covers
   each minor gas at two abundances. Next steps:
   - Commit a generator for the single-gas fixtures (built from
     `makeColumn.py`; today only the `.nc` files are committed).
   - Add a temperature axis.
   - Build a sweep driver on `exort_run_columns` / multi-column input for
     broad training sweeps. Keep it separate from the Δ = 0 regression suite.
     Promote individual sweep points into the suite only once they are
     checked against LBL.
2. **Pre-commit / CI sync gate.** Run `populate3Dmodels.py check --exort`
   automatically, so source ↔ CAM drift can't be committed.
3. **IDL → Python** (independent track). Remaining: `getColumn.py`,
   porting from the experimental `getColumn.pro` at `0e409c3`, which reads
   the CO₂-cloud, `SRF_EMISS` and CARMA fields. Also plotting and utilities,
   then CIA/MT_CKD generation. Move `.pro` files to `tools/legacy_idl/` and
   delete each once its replacement is validated.
4. **`src.cam7.n68equiv` disposition** (CESM3). Gated on cloning
   `NCAR/CESM3-planets` — see `CESM.CAM7.md` Part A.

## 5. Release

1. Items in §1 pass.
2. Refresh `README.md`: once `src.cam.exort` is validated it becomes the
   recommended CAM bundle.
3. Merge `refactor` → `main` and tag `v2.0.0` (maintainer's decision;
   `CHANGES.md` is the release-notes draft).
