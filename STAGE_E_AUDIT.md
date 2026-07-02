# Stage E Audit — Thread-Safety of the Column Solve (2026-07-02)

Pre-implementation audit for Stage E (multi-column batch + OpenMP), per
`REFACTOR_PLAN.md`. **No code is changed by this audit.** It enumerates every
module-scope mutable touched during a column solve so E1 (serial multi-column)
and E2 (OpenMP + per-column MCICA seed) can proceed against a known state
inventory instead of discovering races mid-implementation.

## Scope and method

Audited call path — everything reachable from one column solve:

- library: `exort_run_column` → `run_one_column` → `physconst_setgas` +
  `aerad_driver` (`exort_lib_mod.F90`)
- executable: `main.F90` → `read_input` (io module) → `physconst_setgas` →
  `aerad_driver`
- kernels: `exo_radiation_mod` (driver, two-stream, fluxes/heating),
  `calc_opd_mod`, `rad_interp_mod`, `mcica` + `mcica_random_numbers`,
  `planck_mod`, and every module they `use`.

Method: (1) extracted all non-`parameter` module-scope declarations from each
module on the path; (2) grepped every procedure body on the path for
assignments to each of those names (including DO-index writes, which a naive
"`x =`" grep misses); (3) scanned for `save`d locals, initialized locals
(implicit `save`), `data` statements, and `allocatable`/`pointer` module data.
Files audited are the `OBJS_EXORT` / `OBJS_EXORT_LIB` lists in
`build/Makefile` (the CAM-only `spectral_output_cam.F90` and
`initialize_rad_mod_cam.F90` / `exo_radiation_cam_intr.F90` are not in the
1-D builds and are out of scope; see "CAM constraint" below).

## Verdict

The RT kernel is **very nearly read-only after init**. All k-tables, CIA,
continuum, cloud/haze optics, spectral grids, band-integration limits, and the
Planck table are filled during the init chain and only read during a solve.
There is no dynamic allocation anywhere on the path (all module data is
fixed-size static), no `save`d or initialized locals in any solve-path
procedure, and the MCICA RNG state is call-local. Exactly **two** module-scope
writes occur during (or can occur during) a column solve — findings 1 and 2 —
plus the known single-column I/O buffers in the io module (E1's target).

## State inventory

### Read-only after init (safe to share across threads once init returns)

| Module | Mutable module-scope state | Written by (init only) |
|---|---|---|
| `radgrid` | `wavenum_edge/mid`, Gauss `g_*` arrays, `ngauss_pts`, `solarflux`, `S0`, `g_weight`, P/T grids (`pgrid`, `log10pgrid`, `tgrid`, `tgrid_*`), k-tables (`k_h2o`…`k_co`, `k_major_data`, `k_grey_data`), MT_CKD arrays, CIA tables (`kh2h2`…`ko2co2`), cloud optics (`Qcld*`, `Wcld*`, `Gcld*`, `rel/rei/rei_co2_grid`), haze optics (`kcarma`, `wcarma`), `taumax/ssamax_co2cld`, `camtop`, `ngb` (`data` stmt) | `initialize_kcoeff/solar/cldopts/hazeopts`, `init_model_specific`, `initialize_radbuffer` (`camtop=1`) — but see finding 2 for `camtop` |
| `exo_init_ref` | `lw_iwbeg/iwend/ipbeg/ipend`, `sw_iwbeg/iwend/ipbeg/ipend`, `gw_solflux`, `solflux`, `g_ang_weight`, `g_angle` | `init_ref` (incl. `optimize_bands_sw/lw` — init-time only, despite operating per stellar spectrum) |
| `exo_radiation_mod` | `ptemp_itp` (Planck table) | `init_planck` |
| `physconst` | `scon` | init (`read_namelist` / `exort_init`) |
| `shr_const_mod` | `SHR_CONST_G` | init |
| `sys_rootdir` | `exort_rootdir` | init (library only) |
| `exoplanet_mod` | `solar_file`, `shr_const_scon`, `exo_g`, `do_exo_clouds`, `do_exo_haze` | init |
| `exort_lib_mod` | `initialized`, `finalized` | `exort_init` / `exort_finalize` (read-only during runs) |

Stateless on the solve path: `calc_opd_mod`, `rad_interp_mod`, `mcica`,
`planck_mod`, `kabs` (all parameters), `rayleigh_data` (parameters),
`time_manager` (`get_nstep()` returns the constant 9404), `spmd_utils`
(`masterproc` is a parameter), `mcica_random_numbers` (type definition only —
no module-level instance). No file I/O occurs during a solve; all NetCDF reads
happen in the init chain.

### Written during a column solve — the actual findings

**Finding 1 (primary blocker): `physconst%mwdry` / `physconst%cpair`,
written by `physconst_setgas` per column.**
`exort_lib_mod.F90` `run_one_column` calls
`physconst_setgas(state%mwdry, state%cpdry)` for **every column**; the values
are read throughout the solve (`calc_opd_mod` VMR conversions via `mwdry`,
heating rates via `cpair` in `exo_radiation_mod`). Under OpenMP this is a
write-then-read race: columns with different `mwdry`/`cpdry` would silently
use each other's values, and even identical values are a formal data race.
The executable path calls it once before its single solve (`main.F90:61`), so
it is unaffected. **Stage E fix direction:** carry `mwdry`/`cpair` as
per-column arguments through `aerad_driver` into the kernels (they are column
state, not constants — `column_state_t` already carries them), or as a
fallback declare them `!$omp threadprivate`. Argument-passing is preferred:
it matches the Increment-1 philosophy (state flows through the signature) and
keeps the CAM path unchanged (CAM sets them once at startup).

**Finding 2 (latent race + pre-existing bug): `radgrid%camtop` used as a DO
index at `source/src.main/exo_radiation_mod.F90:1928`.**
`do camtop=swcut,pverp` inside the `if(part_in_tshadow)` branch of the
flux/heating routine writes the module-scope `camtop` (no local shadows it).
Two problems: (a) any thread entering the branch would corrupt `camtop` for
all other threads mid-solve; (b) it is plainly a copy-paste bug — the loop
body indexes with `k` (`dzc(k-1)`, `sw_dTdt(k-1)`…), which the loop never
advances, so the branch is wrong even single-threaded. In practice the branch
is dead in 1-D: `part_in_tshadow` requires `sw_cutoff > 0`, and with the
zeroed horizon/obstruction inputs both `main.F90` and the library pass
(`TCz_obstruct=0` → `sw_cutoff=0`), it stays `.false.` even though
`tshadow_tog=1`. **Stage E fix:** change `do camtop=` to `do k=` (restoring
the obvious intent, same as the two sibling loops). This is a behavior-neutral
fix on every reachable path; the 3-D copies have the same line and should get
the same one-line fix at the 3-D port.

### Call-local state worth recording (not races)

**MCICA RNG — already thread-safe, but the seed is degenerate.**
The MersenneTwister sequence object is a *local* of
`generate_stochastic_clouds` (`mcica.F90`), constructed fresh each call from
`changeSeed` = `permuteseed` = `get_nstep()` — the constant **9404** in the
1-D `time_manager` stub (`calc_opd_mod.F90:835`). So there is no shared RNG
state, but every column draws the **identical** subcolumn sequence. That is
exactly what the planned opt-in per-column seed (`9404 + col`) addresses in
E2 — as a *decorrelation feature with a rebaseline decision*, not a race fix.
The alternate kissvec path (`irnd=0`) seeds from bottom-layer `pmid` digits
and is also call-local; it is hard-disabled (`irnd=1`).

**io module (`io_1D.F90`) single-column buffers — E1's known target.**
The entire input column lives as module-scope state (`TS_in`, `PMID_in`, ten
`*MMR_in` arrays, cloud/haze arrays, `*_zero` companions, `MWDRY_in`, …),
written by `read_input`, read by `main.F90`. This is executable-path-only
(the library bypasses io entirely) and single-column by construction. E1
decouples it: read into per-column storage (reusing `column_state_t`), loop
`aerad_driver` over columns, optional `ncol` dim in the NetCDF — with the
existing no-`ncol` file layout remaining bit-for-bit.

### Minor observations (cleanup candidates, no action required for E)

- `exo_radiation_mod.F90:54` — module-scope `integer :: openstatus` is dead
  (a local of the same name at line 1605 is the one used).
- `radgrid` `ntopcld` / `nlevsRT` are declared, never written, never read
  (all writes commented out); only `camtop` is live.
- `physconst` `gravit`/`rga`/`cpvir`/`zvir` are never written on the 1-D path
  (`cpvir`/`zvir` stay `inf`); nothing on the solve path reads them.

## OpenMP implementation notes for E2 (from this audit)

1. Resolve findings 1 and 2 **before** adding any `!$omp` directive; after
   that, `aerad_driver` is a pure function of its arguments + read-only
   tables, and the E2 loop can be a plain `parallel do` over columns.
2. **Stack:** the solver uses large automatic arrays (e.g.
   `ntot_gpt×ngangles×pverp` in the LW refinement — tens of MB per thread at
   pver=300). Worker threads get `OMP_STACKSIZE`, not the shell's `ulimit`;
   document a required `OMP_STACKSIZE` (or switch the biggest arrays to heap)
   as part of E2, and make the E2 harness set it.
3. **CAM constraint (per the 2026-07-02 handoff):** any edit to shared
   `src.main` files must keep CAM-path behavior unchanged by default —
   per-column seeding stays opt-in (default = today's constant-9404 behavior),
   and the finding-1 fix must not alter how CAM initializes `mwdry`/`cpair`.
4. Verification gates (from `REFACTOR_PLAN.md`): 1 column == reference;
   64-column batch == 64 singles within 1e-10; `OMP_NUM_THREADS=1` vs `8`
   within 1e-10; regression suite 15/15 Δ=0 with the per-column seed OFF.
