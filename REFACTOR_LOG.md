# ExoRT v2 Refactor — Logbook

A human-readable record of what the v2 refactor changed, **why**, and **how to
undo any of it**. This complements `REFACTOR_PLAN.md` (the forward plan) and the
per-stage docs (`STAGE_C_INTEGRATION_MAP.md`, `tests/regression/EXORT_H16_EQUIVALENCE.md`).
If you read one file to re-ground yourself in the refactor, read this one.

## Safety invariants (true at every step)

1. **Your published code is untouched.** The citable, career-foundation ExoRT is
   the **`v1.0.0` tag** and the **`main` branch**. All refactor work lives on the
   **`refactor` branch** and has NOT been merged to `main`. Nothing here changes
   `main` until you explicitly decide to merge.
2. **Clear-sky physics is provably unchanged.** A 13-case regression suite
   (`tests/regression/run_regression.py`) asserts the standard profiles are
   **bit-for-bit identical (Δ=0)** to committed golden baselines. Every
   "should not change physics" commit is verified Δ=0 before it lands.
3. **Every commit is small, single-purpose, and revertible.** `git revert <sha>`
   undoes exactly one step with no entanglement. Nothing was rewritten in place
   that can't be walked back.
4. **Deleted code is preserved in history.** Retired RT bundles live in the
   `v1.0.0` tag; deleted files remain in git history.

To re-verify invariant #2 yourself at any time:
```bash
cd build && USER_FC=gfortran make exort        # gfortran on Apple Silicon
cd ../tests/regression && USER_FC=gfortran python run_regression.py
# expect: "14/14 cases passed", the 13 non-cloud cases all (Δ +0.000)
```

---

## What v2 is

v2 collapses the five historical RT versions into ONE bundle, `src.exort`
(84-band grid, the n84 grid superseding n68; NH3/CO included). The legacy
`n28archean` / `n42h2o` / `n68h2o` bundles were removed from `source/` (preserved
in `v1.0.0`). `n68equiv` / `n84equiv` remain only as HITRAN-2016 comparison
references.

Two SEPARATE efforts were deliberately decoupled so they don't confound each other:
- **The structural refactor** (5 bundles → 1) — done, proven physics-neutral.
- **The HITRAN-2016 → 2024 line-list upgrade** — gated on k-table fixes, kept as
  a separate later effort. `src.exort` currently runs on the *validated*
  HITRAN-2016 native-gas k-files; HITRAN-2024 is a `--exort h24` side-path.

---

## Stage log (most recent first)

Each entry: what changed, why, the commit(s), and how to undo.

### 3-D port — `3dmodels/src.cam.exort` CAM bundle (2026-07-06)

The v2 bundle for CESM1.2.1/ExoCAM, superseding `src.cam.n68equiv`,
`src.cam.n84equiv`, AND `src.cam.n68equiv.haze` (all preserved untouched
as legacy connections; `src.cam7.n68equiv` never touched). 16 files: 13
byte-copies of `source/` (guarded by `populate3Dmodels.py check --exort`,
which PASSES; the pre-v2 bundles are now documented as expected-drift and
must NOT be regenerated), plus bundle-local `sys_rootdir.F90`/`README`,
plus the two CAM-only files extended in `source/` first and copied in:

- **`exo_radiation_cam_intr.F90`** (src.main): the ExoRT↔CESM linkage.
  `exo_radiation_tend`'s argument list is the frozen contract — verified
  token-identical to both old bundles. Internals merged three optional
  paths behind CPP macros (opt-in via `-cppdefs` in CAM_CONFIG_OPTS):
  `-DEXORT_CO2CLD` (CO2 ice clouds: `do_exo_condense_co2` flag +
  CICEWP_CO2/REI_CO2 pbuf fields from the ExoCAM co2condense config;
  wiring recovered from the purged `source/experimental` intr at
  `938a51a^`), `-DEXORT_CARMA` (haze: `carma_exort_get_mmr` +
  NELEM/NBIN==1/40 endrun check; the old haze bundle's
  `carma_exort_optics_init` is NOT used — ExoRT loads its own n84 optics),
  `-DEXORT_SRF_EMISS` (`cam_in%srf_emiss` → `ext_srf_emiss`, clear-sky
  call included). No macros = mainline configs, H2O clouds only.
  Clear-sky companion calls omit condensed-phase args (the `*_zero`
  fabrication of the old bundles is gone = argument-handling Increment 2).
- **`initialize_rad_mod_cam.F90`** (src.exort): PIO loads added for CO2
  ice cloud optics (`radii`→`rei_co2_grid` has no DATA fallback, so the
  read is load-bearing), `rel_grid`/`rei_grid`, and a new PIO
  `initialize_hazeopts` (dims checked, Kext cm²/g→m²/kg) + mpibcasts.

**ExoCAM-side requirement (documented in the bundle README with a
template):** `src.share/exoplanet_mod.F90` must add `do_exo_clouds` and
`do_exo_haze` (the shared driver gates its kernels on them), and
`exo_solar_file` must move to `data/stellar/*_n84.nc`. Inherited fixes
the old bundles still carry: the `part_in_tshadow` `do camtop=` bug and
the haze-kernel bottom-layer off-by-one are both fixed here.

**Verified (no HPC available):** new committed harness
`tests/cam_compile_check/` compiles the full bundle with gfortran against
stubbed CESM modules — mainline + each macro + all three together, all
OK (`run_all.sh`); the `exoplanet_mod` stub doubles as the ExoCAM update
template. Regression 15/15 Δ=0 (the two edited `source/` files are
CAM-only, not in any local target). **NOT yet compiled or run in a real
ExoCAM case** — first HPC step: clear-sky aquaplanet vs an identical
`src.cam.n68equiv` case (expect LW bit-close, SW <0.4% n68→n84 UV
regrid). **Undo:** revert the commit (pure addition + doc edits).

### Stage E2b — per-column gravity and stellar constant (2026-07-03)

Closes the Stage E scope limit: batch columns may now carry their own
gravity and insolation (`934252c`). `column_state_t` gains trailing
`grav`/`scon` fields (append-only ABI; ≤0 or the zero default = "use the
process-level namelist value" — zero-filled states from existing callers
are bit-for-bit unchanged). Gravity: optional `ext_grav` on the
`aerad_driver` keyword tail; the coldens/dzc/masspath block was the ONLY
solve-path gravity read (mcica/calc_opd import `SHR_CONST_G` but never use
it) so no kernel signatures changed. Insolation: no new RT plumbing — a
per-column `scon` rides the existing `ext_msdist` TOA-flux factor
(`msdist = scon_process/scon_col`); uniform spectral scaling only, the
spectrum *shape* stays per-process by design (mixed-star sets = outer loop
of processes, maintainer-agreed). I/O: optional `grav`/`scon` scalar input
variables (`io_1D`, `stackColumns.py`, `makeColumn.py --write-grav`/
`--scon`); cffi + C-test structs extended in lockstep.

**Verified:** regression 15/15 Δ=0; multicol + threaded stages exact; both
library harnesses PASS; new committed gate `percol_config_check.py` — a
3-column batch with per-column (grav, scon) reproduces namelist-configured
singles: gravity EXACT (max|Δ|=0, doubles as a detector for any future
solve-path `SHR_CONST_G` read), scon at 5.5e-16 rel (gate 1e-10; the
namelist path bakes scon into the init normalization, the per-column path
multiplies the flux — same physics, different rounding order). **Undo:**
revert the commit.

### Stage E2 — OpenMP column loop + opt-in per-column MCICA seed (2026-07-03)

Two commits, each gated regression 15/15 Δ=0 + multicol exact + both library
harnesses PASS.

- **Per-column MCICA seed** (`c176bde`): the permute seed moved out of the
  kernels — `calc_opd_cld_h2o` takes `permuteseed` as an argument (all three
  bundles), computed in `aerad_driver` as `get_nstep() + (ext_icol-1)` from a
  new optional keyword-tail arg `ext_icol`; `run_one_column` forwards the
  column index only when the new namelist flag `mcica_percol_seed` is on.
  **Default `.false.` = constant seed 9404 for every column, bit-for-bit**
  (column 1 is bit-identical even with the flag on). Library:
  `exort_set_percol_seed(enable)` (additive C API + cffi method). Also fixed
  in passing: (a) the `calc_opd_cld_h2o` call site passed `cICE_mcica` twice
  and `cLIQ_mcica` never (aliased intent(out) actuals; write-only diagnostics,
  Δ=0); (b) `io_1D` now reads an optional `cfrc` input variable and
  `stackColumns.py` stacks it — previously the 1-D deck could NOT drive H₂O
  cloud fraction (silently zero → H₂O condensate produced no MCICA cloud;
  only library callers could set it). Functional proof: 3 identical
  H₂O-cloud columns → flag off: identical outputs; flag on: col 1 == legacy
  exactly, cols 2/3 decorrelated. **Undo:** revert the commit.
- **OpenMP + heap work arrays** (`b47c07c`): `!$omp parallel do
  schedule(dynamic) if(n>1)` over the column loop in `main.F90` and
  `exort_run_columns`; Makefile appends `-fopenmp`/`-qopenmp` (opt out:
  `make OMPFLAGS=`); `main.F90` timing switched `cpu_time` → `system_clock`
  (wall). The enabler: `-fopenmp` implies `-frecursive`, which moved the
  solver's formerly-STATIC large locals to the stack (~120 MB/solve —
  instant SIGSEGV on macOS; this is also why the code historically ran fine
  on an 8 MB main stack). All 143 solve-path local arrays ≥200 KB
  (`exo_radiation_mod`, `mcica`, `calc_opd_mod` ×3) became `allocatable` —
  so default thread stacks suffice everywhere, including library caller
  threads; no `OMP_STACKSIZE` needed. `multicol_check.py` gained the
  threads-vs-serial stage. Verified: 8-thread batch == serial batch
  max|Δ|=0; 8-column scaling 4.4× on 8 threads, bitwise identical.
  **Undo:** revert the commit (restores automatics and the serial loop).

### Stage E1 — serial multi-column I/O on column structs (2026-07-02)

The io module (`io_1D.F90`) no longer holds ANY column state. `input_profile`
allocates and fills an array of `column_state_t`; `output_data` writes from
`column_state_t`/`column_result_t` arrays (K/s → K/day conversion moved
inside); the ~50 module-scope `*_in`/`*_zero` variables and
`initialize_to_zero` are deleted. The per-column solve now lives in ONE
shared place — new `exort_column_run.F90` (`run_one_column`, moved verbatim
from `exort_lib_mod`, which now uses it) — and `main.F90` loops it serially
over columns (E2 parallelizes this loop).

**File format:** `RTprofile_in.nc` may carry an optional `ncol` dimension.
Absent = classic single-column file, bit-for-bit the pre-E1 path (regression
15/15 Δ=0; output layout unchanged, no `ncol` dim, scalars stay on `ONE`).
Present = every variable carries a trailing column dimension (Fortran order:
`tmid(pver,ncol)`, `ts(ncol)`, `carmammr(pver,nelem,nbin,ncol)`) and the
output gains a matching `ncol` on every variable. `tools/stackColumns.py`
builds multi-column inputs from single-column files (optional variables
unioned, zero-filled where absent). Columns share the process-level runtime
config (namelist); only per-column state varies.

**Verified:** all 4 targets build; regression 15/15 Δ=0;
`tests/regression/multicol_check.py` (new, whitelisted in .gitignore) proves
a 3-column batch (TS250K/TS300K/TS340K × G2V) reproduces the three single
runs EXACTLY (max |Δ| = 0 on every output variable); both library harnesses
PASS unchanged. **Undo:** revert the commit (restores module-state io and
the single-column-only driver).

### Stage E pre-fixes — audit findings 1 & 2 resolved (2026-07-02)

Both fixes from `STAGE_E_AUDIT.md`, each verified regression 15/15 Δ=0 and
(for finding 1) library harnesses PASS. After these, a column solve writes
**no module-scope state** — the prerequisite for the E2 OpenMP loop.

- **Finding 2** (`f5b48ea`): `do camtop=` → `do k=` in the dead
  `part_in_tshadow` heating branch of `exo_radiation_mod.F90` (copy-paste
  bug that used the radgrid module variable as a DO index). 3dmodels copies
  get the same line at the 3-D port. **Undo:** revert the commit.
- **Finding 1**: per-column dry-air properties now flow through the driver.
  `aerad_driver` gained `ext_mwdry`/`ext_cpdry` appended to the optional
  keyword tail (absent → physconst module values = the CAM path, unchanged);
  the coldens block uses `mwdry_col`, `rad_postcalc` takes `cpair_col` as an
  argument, and `calc_opd_gas` takes `mwdry` as an argument in all three
  bundles (src.exort + the two legacy references, which link against the
  shared driver). `physconst_setgas` **deleted** (both callers — `main.F90`
  and `exort_lib_mod` — now pass the values by keyword); the shim's
  `mwdry`/`cpair` keep Earth-air defaults and exist only as the CAM
  fallback. **Undo:** revert the commit (restores setgas and the module
  writes).

### Stage E audit — thread-safety inventory of the column solve (2026-07-02)

Docs-only checkpoint (`STAGE_E_AUDIT.md`); no code changed. Enumerated every
module-scope mutable reachable from one column solve (declaration sweep +
assignment grep + save/data/initialized-local/allocatable scans over the
`OBJS_EXORT_LIB` file list). Result: the kernel is read-only after init
except for exactly two writes — (1) `physconst_setgas` rewriting
`mwdry`/`cpair` per column in the library path (the primary OpenMP blocker;
fix by passing them through `aerad_driver` as column state), and (2) a
latent `do camtop=` copy-paste bug at `exo_radiation_mod.F90:1928` inside
the dead-in-1-D `part_in_tshadow` branch (writes the `radgrid` module var;
body indexes with `k` — one-line `do k=` fix). MCICA RNG state is call-local
(thread-safe) but seeded with the constant 9404 for every column — the E2
per-column seed is a decorrelation feature, not a race fix. `io_1D`'s
single-column `*_in` module buffers are confirmed as E1's decoupling target.
No allocatables, no saved locals, no solve-time file I/O anywhere on the
path. **Undo:** delete `STAGE_E_AUDIT.md` and this entry.

### Argument-handling Increment 1 — `aerad_driver` optional keyword tail (2026-07-02)

Implements the surviving idea from `ARGUMENT_HANDLING.md` (D2), revised for
the single-bundle world. `aerad_driver`'s signature is now:

- **Mandatory positional core:** thermodynamic state, geometry, surface
  albedos, and **all 10 gas MMRs**. A gas is turned off by zero MMR
  (`calc_opd_gas` short-circuits zero-abundance gases), never by omitting an
  argument.
- **Optional keyword tail (after `sol_toa`):** H2O clouds
  (`ext_cicewp/cliqwp/cfrc/rei/rel`, all five or none), CO2 ice clouds
  (`ext_cicewp_co2/rei_co2`, both or neither), CARMA haze (`ext_carmammr`),
  surface emissivity (`ext_srf_emiss`). **Pass by keyword only; new optional
  args are appended, never inserted** — so positional call sites (1-D,
  library, 3dmodels, ExoCAM SourceMods) never break when physics is added.
  This placement fixes the Stage C mid-list optionals, which had already
  silently broken the source-side `exo_radiation_cam_intr.F90` call sites
  (positional args landed on the wrong dummies — uncompilable, unnoticed
  because no 1-D target compiles that file).
- **Semantics (documented in the driver header):** flag on + arg present →
  physics active; flag on + arg absent → kernel not called, zero
  contribution (this is how a clear-sky companion call works); flag off →
  arg ignored. Kernels (`calc_opd_*`) keep all-mandatory arguments and no
  `present()` logic — gating happens at the driver's call sites.

Call sites updated: `main.F90` and `exort_lib_mod.F90` pass the tail by
keyword unconditionally (their arrays always exist, zero-filled → behavior
identical); `exo_radiation_cam_intr.F90` (source copy) repaired to the new
form — the clear-sky call now **omits** the cloud args (deleting the
`cicewp_zero/cliqwp_zero/cfrc_zero` fabrication), the full-sky call passes
H2O clouds by keyword. CO2-cloud/haze/emissivity plumbing through the CAM
interface remains for the 3-D port, which also updates the 3dmodels copies.

**Verified:** exort/n68equiv/n84equiv all build; regression 15/15 Δ=0
(including gated CO2-cloud and haze cases); libexort rebuilt, Python + C
harnesses PASS (float32-limited agreement unchanged, determinism exact).

**Undo:** revert the single Increment-1 commit.

### Stage D — `iso_c_binding` library + Python binding (2026-07-02)

`src.exort` is now callable in-process from C and Python. Pure addition: the
executable path is untouched (15/15 regression Δ=0 after the change).

**What was added:**

- `source/src.main/exort_column_mod.F90` — `column_state_t` / `column_result_t`,
  bind(c) all-double structs carrying one column's inputs/outputs. Field order
  is the ABI contract shared with the Python and C consumers.
- `source/src.main/exort_lib_mod.F90` — the C API:
  `exort_get_dims` (query compiled pver/pverp/nwave/CARMA dims),
  `exort_init(data_root, solar_file, scon, g, do_clouds, do_haze)` (one per
  process; empty string / non-positive value keeps the compiled default),
  `exort_run_column`, `exort_run_columns` (serial batch loop),
  `exort_finalize` (bookkeeping; re-init not supported in Stage D).
  `run_one_column` mirrors main.F90's per-column sequence exactly
  (physconst_setgas, dry-pressure derivation, terrain trivia, aerad_driver).
- `build/Makefile` — `make libexort` → `run/libexort.dylib|.so` (same object
  list as exort.exe with main.F90 swapped for the two lib modules).
- `tools/exort_pytools/` — cffi ABI-mode binding (`exort_api.py`, dims queried
  at load so nothing is hardcoded) + acceptance harness (`verify_lib.py`).
- `tests/lib/` — dimension-agnostic standalone C test (`test_exort_c.c`,
  dlopen + netCDF-C) with its own Makefile (`make run`).
- `source/src.main/sys_rootdir.F90` — `exort_rootdir` demoted
  parameter → variable so `exort_init`'s data_root can override it at runtime.
  Compiled default unchanged; the 3-D copies of sys_rootdir are untouched.

**Verification:** Python and C harnesses both reproduce the committed golden
baselines (TS300K_G2V, TS250K_BB3400, 2barCO2_dry_Mars_G2V) with max relative
difference ~6e-8 on every output — the float32 storage precision of the
baseline files (`output_data` writes nf_real), i.e. the library's doubles agree
with the executable's to at least that. In-memory repeat-call and batch-vs-
single determinism are exactly 0. Gotcha found on the way: the baselines were
generated with the `_n84` stellar files (the harness maps `_n68 → _n84` for
exort runs); a library caller must pass `G2V_SUN_n84.nc` etc. to reproduce them.

**Thread-safety contract:** tables are read-only after `exort_init`;
`exort_run_column` itself is NOT yet thread-safe (module scratch state in the
RT kernels). Per-column MCICA seeding and OpenMP arrive in Stage E.

**Undo:** revert the Stage D commit (delete the two lib modules, the Makefile
target, `tools/exort_pytools/`, `tests/lib/`; restore the `exort_rootdir`
parameter). Nothing else depends on it.

### Delete `source/experimental/` (Stage C consumed it)

Stage C absorbed the experimental thread: CO2-cloud kernel (C1b), srf_emiss
(C1-emiss), makeColumn.pro cloud writers → Python (C1c), haze (C3, kernel taken
from the published 3-D bundle instead of the shr2 sketch). The `n68equiv_exp`
build target was removed from `build/Makefile` with it. The canonical
`Khare_haze.txt` lives on at `data/aerosol/refractive_indices/`.

**Recovery: the last commit containing the directory is `0e409c3`**
(`git show 0e409c3:source/experimental/src.n68equiv_exp/main.F90`, or
`git checkout 0e409c3 -- source/experimental` to resurrect it all).

**Concepts preserved only in history — needed for future stages:**

1. **`_CLD` clear-sky / cloud-forcing double-run** (deferred stage, open item
   below). Reference implementation: experimental `main.F90` calls
   `aerad_driver` twice — clear-sky first, then cloudy — and `output.F90`
   writes the cloudy fluxes as parallel `*_CLD` variables alongside the
   clear-sky set, so cloud forcing = difference of the two in one output file.
   **Convention note from its header: the 1-D code treats CLEAR-sky as the
   default run and cloudy as the `_CLD` extra — the opposite of the 3-D model,
   where cloudy is default and clearsky is the diagnostic.**
2. **`getColumn.pro` (experimental variant)** — extracts a 1-D column deck
   from 3-D Mars GCM output including the Stage C fields (`CLDICE_CO2`,
   `REI_CO2`, `SRF_EMISS`, CARMA arrays), which `tools/getcolumn.pro` does not
   read. Port from this when writing the future `getColumn.py`.
3. **`plotspectra_1D.pro` (experimental variant)** — reads the `*_CLD`
   outputs (59 references); it documents the expected output-variable naming
   for the double-run stage and the cloud-forcing plots that consume it.

### Stage C3 — CARMA haze aerosols (the haze work)

Additive, gated OFF by default. Haze does nothing unless BOTH the runtime flag
`do_exo_haze=.true.` AND the input file carries `carmammr(pver,nelem,nbin)`.
Existing decks have neither, so haze-off physics is bit-for-bit unchanged
(verified: 14/14 pre-existing cases Δ=0).

- **C3 code `cc892f1`** — Filled the `calc_aeropd()` stub and wired haze opacity
  end-to-end in 1-D: `do_exo_haze` namelist flag, optics loader
  (`initialize_hazeopts`), band-indexed kernel (τ = q·Δp/g·Kext presummed over
  CARMA elements/bins), solver merge at the aerosol placeholder in
  `rad_precalc`, optional `carmammr` input. The CARMA-module coupling is severed
  on the 1-D side (haze formation lives outside ExoRT; `carmammr` comes from the
  deck), per the maintainer-confirmed Stage C design boundary. The published 3-D
  haze kernel's off-by-one (layer ik−1 mass paired with pdel(ik), bottom rad
  level never filled) was not reproduced — the 1-D kernel matches the driver's
  coldens convention — and was subsequently **fixed in the 3-D bundle too**
  (maintainer-confirmed; see the layer-indexing fix commit below). *Undo:*
  `git revert cc892f1`.
- **C3 data `df2b7bb`** — `tools/regrid_haze_optics.py` +
  `data/aerosol/haze_n84_b40_{mie,fractal_interp}.nc`. Bands 1–68 are verbatim
  copies of the validated n68 tables (identical band edges, verified); the 16 UV
  bands are a nearest-band extension and are **PROVISIONAL** (global attr says
  so). Replace via a proper 84-band recomputation from
  `data/aerosol/refractive_indices/`; UV haze extinction is underestimated until
  then. *Undo:* `git revert df2b7bb`.
- **C3 test `5489019`** — `makeColumn.py` haze writer (`add_haze_layer` +
  `carmammr`), committed hazy TS300K fixture + generator, `TS300K_haze_G2V`
  regression case + baseline. **Suite is now 15/15** (14 haze-off Δ=0 + 1 hazy
  gated). *Undo:* `git revert 5489019`.
- **Kernel renames** — the four optical-depth kernels are now named after the
  module: `calc_opd_gas`, `calc_opd_cld_h2o`, `calc_opd_cld_co2`,
  `calc_opd_aero` (were `calc_gasopd` / `calc_cldopd` / `calc_cldopd_co2` /
  `calc_aeropd`; H₂O clouds now tagged explicitly). Renamed in `src.exort`, the
  shared `src.main/exo_radiation_mod.F90`, and the legacy-bundle stubs; the 3-D
  bundles and experimental copies keep the old names until their respective
  ports. Pure rename, 15/15 Δ=0. Historical log entries below keep the names
  current at the time.
- **Kernel uniformity** — all four `calc_*opd` kernels now loop `ik=1,pverp`
  over driver-premapped radiative-grid arrays with no index arithmetic:
  `aerad_driver` maps the air mass path to the rad grid (`masspath(1)=pint(1)/g`
  = buffer layer to space, `masspath(k)=pdel(k-1)/g`, same mapping as
  `coldens`) and `calc_aeropd` consumes it directly. Bit-for-bit (15/15 Δ=0,
  hazy case included). The 3-D haze bundle kernel keeps its (correct)
  `ik=2,pverp` + `pdel(ik-1)` form until the 3-D port adopts the driver-mapped
  style.
- **Legacy build fix `02ba2ed`** — n68equiv/n84equiv had failed to **link since
  C1b** (verified by building at `510500d`): the shared `exo_radiation_mod`
  gained `calc_cldopd_co2` (C1b) and `calc_aeropd` (C3) call sites that only
  src.exort carried. Restored with minimal zero-opacity stubs + loud-stop
  `initialize_{cld,haze}opts` stubs (Stage C physics is exort-only). All legacy
  targets build; gas_sweep CO2 slice reproduces. *Undo:* `git revert 02ba2ed`
  (re-breaks the legacy builds). **Side observation:** since the `d77edb9`
  rebaseline pinned `src.exort/kabs.F90` to h16, `gas_sweep.py`'s exort column
  no longer reflects HITRAN-2024 (its h24 swap now no-ops; the column duplicates
  n84/h16, linelist Δ ≡ 0). Revisit when the h24 re-fit lands.

**Stage C3 verification:** one `exort.exe`, haze toggled by namelist — clear
TS300K OLR 268.726 / SWDN_SFC 233.9; hazy (visible τ≈0.5) OLR 266.550 /
SWDN_SFC 194.3 / SWUP_TOM 80.1→77.1 (haze dims the surface, dark tholin lowers
planetary albedo, high haze mildly reduces OLR — the expected anti-greenhouse
signature); 5× loading is monotonic (OLR 258.3, SWDN_SFC 131.7). The hazy
baseline depends on the provisional optics; **rebaseline `TS300K_haze_G2V` when
the maintainer regenerates the 84-band haze tables.**

### Regression rebaseline to src.exort (HITRAN-2016)  — `d77edb9`
- **What:** Re-pointed the regression suite to validate `exort.exe` on the
  HITRAN-2016 native-gas k-files instead of the legacy `n68equiv.exe`. Regenerated
  the 13 golden baselines with exort.
- **Why:** Anchor the suite to the actual v2 code path on trusted data, and
  decouple the refactor from the unvalidated HITRAN-2024 tables.
- **Proof:** `tests/regression/EXORT_H16_EQUIVALENCE.md` records that exort-on-h16
  reproduces the retired n68 physics (LW bit-for-bit; SW differs only by the known
  <0.4% n68→n84 stellar-UV regrid, G2V cases only).
- **Undo:** `git revert d77edb9` (restores n68equiv baselining).

### gas_sweep tooling  — `cad1643`
- **What:** Added `--save-spectra` and a CO2/code-isolation `--with-n84h24` column
  to `gas_sweep.py`, plus `tools/plot_gas_sweep_spectra.py`.
- **Why:** Diagnostic tooling for the HITRAN h16-vs-h24 comparison. No RT change.
- **Undo:** `git revert cad1643`.

### Stage C — CO2 ice clouds + surface emissivity (the cloud work)

All additive, all gated OFF by default. Clouds do nothing unless BOTH the runtime
flag `do_exo_clouds=.true.` AND the input file contains cloud condensate fields.
Your existing decks have neither, so they are unaffected.

- **C0 `03f3007`** — Integration map (`STAGE_C_INTEGRATION_MAP.md`); deleted 3
  stale `shr2/*.F90~` editor backups. *Undo:* `git revert 03f3007`.
- **C1a `edbdb6f`** — Restored the 1-D cloud-optics loader (it had been commented
  out, so 1-D clouds never worked) behind a new off-by-default `do_exo_clouds`
  flag. Loads the H2O liquid/ice Mie tables at init. *Undo:* `git revert edbdb6f`.
- **C1b `8dd204e`** — Added the CO2-ice cloud kernel (`calc_cldopd_co2`) and wired
  it into the two-stream solver (band-indexed, cloud fraction = 1, no MCICA — the
  condensate is supplied as a grid-box mean from outside ExoRT). *Undo:*
  `git revert 8dd204e`.
- **C1-emiss `778469e`** — Optional surface thermal emissivity input (`srf_emiss`);
  was hardcoded ε=1.0. Absent → ε=1.0, so Δ=0 preserved. *Undo:* `git revert 778469e`.
- **C1c `cce7247`** — Python `makeColumn.py` cloud + `srf_emiss` writers (ports the
  cloud logic from the experimental `makeColumn.pro`); committed CO2-cloud test
  fixture + its generator. *Undo:* `git revert cce7247`.
- **C1d `63e27fd`** — Made `do_exo_clouds` a RUNTIME namelist flag (was compile-time)
  so one build iterates clear/cloudy decks without recompiling; added the
  `2barCO2_co2cloud_Mars_G2V` regression case (suite now 14/14). *Undo:*
  `git revert 63e27fd`.

**Stage C verification:** one `exort.exe`, both modes via namelist — clear Mars
OLR 92.744, `do_exo_clouds=.true.` on the cloudy fixture OLR 33.919 (CO2 ice cloud
strongly reduces OLR, the expected early-Mars greenhouse behavior). The 13 clear
cases stay Δ=0.

**To undo all of Stage C at once** (keeping the rebaseline): revert the 6 cloud
commits `03f3007 edbdb6f 8dd204e 778469e cce7247 63e27fd` (newest-first), or
branch from `cad1643` (the commit just before Stage C began).

---

## What is NOT done (open, decision-gated)

- **HITRAN-2024 k-tables** — CO2 (far-IR χ-factor pipeline bug) and possibly H2O
  need an offline HELIOS-K re-fit before h24 can be adopted. C2H6 is fixed.
  src.exort stays on h16 until then.
- **84-band haze optics regen** — the committed `haze_n84_b40_*.nc` are
  provisional (UV bands 69–84 are a nearest-band extension of band 68). The
  maintainer regenerates them properly from
  `data/aerosol/refractive_indices/`, then reruns
  `tools/regrid_haze_optics.py`-free validation and rebaselines
  `TS300K_haze_G2V`. The C3 kernel/wiring is done (see Stage C3 above).
- **3-D validation on HPC** — `src.cam.exort` (2026-07-06) is code-complete
  and compile-checked against stubbed CESM, but has never been built inside
  a real ExoCAM case. Needs: the small ExoCAM `exoplanet_mod` additions
  (template in the bundle README / `tests/cam_compile_check/
  exoplanet_mod_stub.F90`), an `_n84` `exo_solar_file`, then a clear-sky
  aquaplanet comparison against `src.cam.n68equiv`, then gated CO2-cloud /
  haze / srf-emiss cases per config.
- **Clear-sky / cloud-forcing double-run** (`_CLD` outputs) — deferred to its own
  stage. The reference implementation (dual `aerad_driver` call, `*_CLD` output
  variables, clear-sky-default convention) lives in git history at `0e409c3`
  under `source/experimental/` — see the deletion entry above for the details.
- **Merge to `main`** — not done; the maintainer is not ready. Refactor stays on
  the single `refactor` branch (no sub-branches, by preference).

---

## How to convince yourself, from scratch

```bash
# 1. Your trusted v1 is right there, unchanged:
git show v1.0.0 --stat | head
git log main -1            # main is your maintenance line

# 2. The v2 bundle builds and reproduces the standard cases:
cd build && USER_FC=gfortran make exort
cd ../tests/regression && USER_FC=gfortran python run_regression.py   # 15/15

# 3. Run a known column yourself and eyeball it:
cd ../../run
cp ../tests/regression/fixtures/RTprofile_in_TS300K.nc RTprofile_in.nc
./exort.exe        # inspect RTprofile_out.nc against your expectations
```
