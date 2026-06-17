# Interpolation Scheme — Cost Analysis & Planned Fixes

**Status:** Analysis complete; two fixes designed but **NOT yet applied**. To be
scheduled against `REFACTOR_PLAN.md` (see "Sequencing" below).

**Date:** 2026-06-17
**Bundle:** `source/src.exort` (n84 grid: `ntot_wavlnrng=84`, `nspecies=8`,
`ngauss_8gpt=8`).

---

## Motivation

Planned expansion of the k-coefficient P–T reference grid:

| Dimension | Current | Proposed | Driver |
|-----------|---------|----------|--------|
| Pressure (`kc_npress`) | 61 | 71 | raise P_max to 100 bar |
| Temperature (`kc_ntemp`) | 17 | 68 | raise T_max from 500 K → 3000 K+ at uniform 25 K steps |

Concern raised: **computational cost** of this expansion.

---

## What the interpolation actually does

Code paths examined:
- `source/src.exort/calc_opd_mod.F90:231–267` — bracketing search for the P and
  T reference indices (once per layer per call).
- `source/src.exort/rad_interp_mod.F90` `bilinear_interpK_grey` (lines 162–273)
  and `bilinear_interpK_8gpt_major_gptvec` (lines 46–158) — the bilinear blend.

**Key finding: bilinear interpolation cost is independent of grid size.** Both
routines receive `p_ref_index` / `t_ref_index` already located, then read a fixed
number of corner values from `Kdata`:
- `bilinear_interpK_grey`: 4 corners, one blend.
- `bilinear_interpK_8gpt_major_gptvec`: 4 corners × 8 g-points = 32 reads.

No arithmetic in these routines scales with `kc_npress` or `kc_ntemp`.

The **inner loop** in `calc_opd_mod.F90` (`do sp = 1, ntot_wavlnrng`, line 277)
calls 8 grey interps + 1 major interp per band — this dominates runtime and is
also grid-size-independent.

The **only** grid-dependent work is the bracketing search (next section).

---

## Cost verdict

| Concern | Verdict |
|---------|---------|
| **Per-call interp CPU** | No meaningful change. Fixed corner count regardless of grid. |
| **Bracketing search** | Scales with grid; T-scan 17→68 ≈ 4× more compares. Once per layer (×300), not per band/g-point. Small absolute cost, but see asymmetry below. |
| **Memory** | **The real cost.** ~4.6× growth (61×17 → 71×68), ~0.95 GB → ~4.4 GB resident. Multiplies per MPI rank in 3-D (ExoCAM). |

### Memory footprint detail

`k_major_data` = `nspecies(8) × ntot_wavlnrng(84) × ngauss_8gpt(8) × kc_npress × kc_ntemp × 8 bytes`

| Grid | k_major_data | k_grey_data | 8× individual `k_h2o…k_co` | **Total resident** |
|------|-------------|-------------|----------------------------|--------------------|
| 61×17 (current) | 446 MB | 56 MB | 446 MB | **~0.95 GB** |
| 71×68 (proposed) | 2.08 GB | 260 MB | 2.08 GB | **~4.4 GB** |

The individual `k_h2o, k_co2, … k_co` arrays (`radgrid.F90:372`) **duplicate** the
data packed into `k_major_data`/`k_grey_data`. Verified they are read **only at
init** in `model_specific.F90` (lines 54–61 fill `k_major_data`; lines 97–104
build `k_grey_data`); never touched in the timestep loop. They are pure init-time
scratch carried for the whole run.

---

## The temperature-search asymmetry (the real reason to fix the search)

The T-bracketing in `calc_opd_mod.F90:251–261` starts at `t_ref_index = kc_ntemp`
(the **maximum**) and walks **down** until `tgrid(t_ref_index) <= t_kgas`.

Most atmospheric layers are cold (≪ 500 K), so today they already sit near the
bottom of a 17-point grid — the scan is short. Push `T_max` to 3000 K with 68
points and a 250 K layer now sits ~60 indices from the start of the scan: **every
layer pays ~60 comparisons it did not before.** This is the scaling cost the
expansion introduces. Binary search bounds it to `log2(68) ≈ 7` regardless of
`T_max`.

(Pressure search, lines 234–239, has the same downward-from-max structure but the
grid grows only 61→71, so the effect there is minor.)

---

## Fix 1 — Replace linear bracketing with binary search

**Where:** `calc_opd_mod.F90:231–267`.

**Contract the search must preserve (verified against the interp routines):**
find the largest index `i` such that `grid(i) <= val`, clamped to `[1, n-1]`;
values below `grid(1)` → `i=1` with `val` snapped to `grid(1)`; values at/above
`grid(n)` → `i=n-1` (interp routines bracket `i, i+1`). For T above max, snap
`t_kgas = tgrid(kc_ntemp)` so the interp `temp` factor = 1.0 (pure clamp, matches
old behavior). For P above max, return `kc_npress-1` so the interp extrapolates
via its `else` branch (matches old reliance on the `p_ref_index==kc_npress`
branch).

**Helper to add** (in the `contains` section of `calc_opd_mod.F90`, or alongside
the interp routines in `rad_interp_mod.F90`):

```fortran
  ! Return largest index i in [1, n-1] such that grid(i) <= val.
  ! grid must be monotonically increasing. Values below grid(1) clamp to 1;
  ! values at/above grid(n) clamp to n-1 (the interp routines bracket i,i+1).
  pure function find_ref_index(grid, n, val) result(idx)
    integer,  intent(in) :: n
    real(r8), intent(in) :: grid(n)
    real(r8), intent(in) :: val
    integer :: idx
    integer :: lo, hi, mid

    if (val <= grid(1)) then
      idx = 1
      return
    else if (val >= grid(n)) then
      idx = n - 1
      return
    endif

    lo = 1
    hi = n
    do while (hi - lo > 1)
      mid = (lo + hi) / 2
      if (grid(mid) <= val) then
        lo = mid
      else
        hi = mid
      endif
    enddo
    idx = lo
  end function find_ref_index
```

**Replacement for lines 231–267:**

```fortran
      ! Find pressure coordinates for k-coefficients
      pressure = log10(pmid(ik))                       ! log pressure
      p_ref_index = find_ref_index(log10pgrid, kc_npress, pressure)
      if (pressure <= log10pgrid(1)) pressure = log10pgrid(1)   ! clamp low-P

      ! Find temperature coordinates for k-coefficients
      temperature = tmid(ik)                           ! actual temperature [K]
      t_kgas = temperature
      t_ref_index = find_ref_index(tgrid, kc_ntemp, t_kgas)
      if (t_kgas < tgrid(1))        t_kgas = tgrid(1)         ! clamp below grid
      if (t_kgas > tgrid(kc_ntemp)) t_kgas = tgrid(kc_ntemp)  ! clamp above grid
```

**Verification:** This changes index-finding only, not arithmetic. A correct swap
is **Δ=0 on all 13 regression cases at the current 61×17 grid** —
`cd tests/regression && python run_regression.py`. Achieve Δ=0 at 17 temps
*before* regenerating k-tables on the bigger grid. Any Δ≠0 means an edge case
(exact grid-point hit, or above-max-T) is not reproducing the old bracketing.

**Honest payoff:** The search runs once per layer (×300), not per band (×84) or
g-point (×8), so the inner interp loop still dominates wall-clock by orders of
magnitude. Do not expect a measurable speedup. The real win: removes the
*scaling-with-T_max* property, so cost is bounded no matter how far the grid is
pushed; and the code is clearer.

---

## Fix 2 — Deallocate the duplicate individual k arrays after init

**Where:** `radgrid.F90:372` (declarations); init in
`initialize_rad_mod_1D.F90` / `initialize_rad_mod_cam.F90`; consumed in
`model_specific.F90:54–61, 97–104`.

**Change:** make `k_h2o, k_co2, k_ch4, k_c2h6, k_o3, k_o2, k_nh3, k_co`
`allocatable`; allocate in the init routine, deallocate immediately after
`model_specific` packs them into `k_major_data` / `k_grey_data`.

**Payoff:** reclaims a full redundant copy — ~446 MB today, **~2.08 GB at 71×68**
— before any RT computation runs. Roughly halves steady-state footprint at the
big grid (~4.4 GB → ~2.3 GB). Verified safe: these arrays are init-only scratch,
never read in the timestep loop. Multiplies in 3-D (per-rank).

**Note for 3-D:** the `_cam` init path and the `mpibcast` block in
`initialize_rad_mod_cam.F90` must be kept consistent if the arrays become
allocatable. `model_specific.F90`, `calc_opd_mod.F90`, and `radgrid.F90` are
shared/synced with `3dmodels/src.cam.*` (see CLAUDE.md sync rule) — any edit must
be mirrored and `diff`-clean.

---

## Additional recommendation (not a code fix — a data/grid choice)

**Do not use a uniform 25 K grid to 3000 K.** That uniform spacing is the source
of ~95% of the memory blowup (T-dim 51→68 of the proposed expansion).
k-coefficients vary smoothly and slowly at high T. A **non-uniform `tgrid`**
(25 K up to ~500 K, then 50–100 K steps above) would roughly halve the T
dimension with negligible accuracy loss. **The interpolation code is already
grid-agnostic** — it reads `tgrid(t_ref_index)` directly (`rad_interp_mod.F90:122,
239`), so a non-uniform grid needs **zero code changes**, only a different `tgrid`
array and matching k-tables. (Fix 1's binary search also requires only
monotonicity, not uniform spacing.)

---

## Sequencing — assess against `REFACTOR_PLAN.md`

Both fixes are deferred pending a decision on *when*. Considerations for that
assessment:

- **Both fixes are scoped entirely to `src.exort`** and its interp/init path —
  the same surface that Stages C–G operate on. They do not conflict with the
  v2 single-bundle direction; they harden the bundle that survives.
- **Fix 2 touches files shared with `3dmodels/src.cam.*`** (`radgrid.F90`,
  `model_specific.F90`, `calc_opd_mod.F90`) and the `_cam` init + `mpibcast`
  path. Per the CLAUDE.md sync rule, mirror edits and keep `diff`-clean. Cleanest
  to land *after* any source↔3dmodels dedup work settles, or coordinated with it.
- **The grid expansion itself is a physics/data change** (new k-tables on a new
  P–T grid) and forces a **regression rebaseline** (golden baselines change).
  That is a deliberate, gated event — natural to pair the grid expansion + both
  fixes + rebaseline into one coherent change rather than spreading them.
- **Recommended order when the time comes:**
  1. Fix 1 (binary search) at the *current* 61×17 grid → prove **Δ=0** on all 13
     cases. Lands as a pure no-op refactor commit.
  2. Fix 2 (deallocate duplicates) → still Δ=0 (memory only, no arithmetic).
  3. Expand the grid + regenerate k-tables (prefer non-uniform `tgrid`) →
     intended physics change → **rebaseline** the regression suite.
- **Candidate hook in the plan:** the search/memory fixes are a good fit
  alongside **Stage E** (multi-column batch + OpenMP), which already decouples
  module-scope state and is the point where per-column / memory concerns get
  attention. The grid expansion is independent of the staged refactor and can be
  timed by science need — but should not precede Fix 1 (so the rebaseline isn't
  conflated with a search-logic change).

---

## Parallelism & memory — what the footprint actually limits

**Core distinction:** the k-tables (`k_major_data`, `k_grey_data`, the duplicate
`k_*`) are **shared, read-only reference data**, populated once at init and never
written in the timestep loop (plan thread-safety contract, `REFACTOR_PLAN.md:260`:
"tables read-only after init"). Per-column working state (the `tau_grey`, `ugas`,
`ans` vectors, the column profile) is tiny by comparison. So the memory question
reduces to: **how many copies of the tables exist?** That is an architecture
choice, not a fixed cost. Rule of thumb: **memory = footprint × number of
processes. Threads are free; processes are not.**

| Scenario | Parallelism model | Table copies | Does Fix 2 → more parallelism / speed? |
|----------|-------------------|--------------|----------------------------------------|
| **1-D, OpenMP** (Stage E plan) | Shared-memory threads | **1 per process** | No — tables already shared across threads; adding threads costs only tiny per-column scratch. Fix 2 lowers absolute footprint only. |
| **1-D, separate processes** | N independent `.exe` / MPI | 1 per process | **Yes, directly.** `instances = RAM / footprint`; Fix 2 roughly doubles how many fit. This is the case the "memory limits parallelism" intuition fits. |
| **3-D, CESM (126 MPI ranks, 3312 cols)** | Distributed-memory MPI | **1 per rank** (all 126) | **No speedup, and not a concern on Milan.** |

### 1-D parallel runs (future Stage E)

The lever is *how* you parallelize. Stage E specifies **OpenMP** ("OpenMP loop in
`main.F90`") — shared-memory threads sharing one copy of the read-only tables.
There, the table memory is paid **once per process**, not per thread, so memory
does **not** gate thread count; Fix 2 only reduces the absolute footprint. If
instead you run **N independent processes** (separate `.exe` or MPI ranks), each
carries its own full copy and `instances = RAM / per-instance-footprint` — here
Fix 2 (~2.3 GB vs ~4.4 GB at 71×68) nearly doubles the packing density. OpenMP is
the memory-efficient path and the one the plan takes.

### 3-D climate model (CESM, 3312 columns / 126 cores)

The hypothesis "less memory → more parallel columns → faster" does **not** hold
here:

- **Parallelism is fixed by the decomposition, not by memory.** CESM distributes
  3312 columns across 126 MPI ranks (~26 cols/rank) at setup — a pure
  distributed-memory model. Each rank is a separate process with its **own full
  copy** of the tables (every individual `k_*` is `mpibcast` to all ranks,
  `initialize_rad_mod_cam.F90:242–249`). Rank count is set by the job submission;
  freeing RAM does not add ranks.
- **Time goes to arithmetic + MPI comms**, not memory. The tables far exceed cache
  (446 MB vs ~MB L3) regardless of Fix 2, so the cache story is unchanged.
- **Memory is not a constraint on the target hardware.** On AMD Milan nodes the
  user has 100s of GB of RAM available, so even the 71×68 grid × 126 ranks fits
  comfortably. Fix 2 is therefore **neither a speed nor a feasibility lever in
  3-D** on this hardware — it remains worthwhile only for footprint hygiene and
  to cut wasted init-time `mpibcast` traffic of the duplicate arrays.

**Takeaway:** Fix 2's parallelism payoff is real **only** for the multi-*process*
1-D case. It is neutral for OpenMP 1-D (tables already shared) and effectively
moot for 3-D on Milan (memory-abundant). Speedups in 3-D must come from the
arithmetic, not the footprint — see "Runtime performance opportunities" below.

---

## Runtime performance opportunities (physics-preserving)

Investigation of the hot path (`calc_opd_mod.F90` layer loop, lines 172–304) for
algorithmic/procedural speedups that **do not change physics** (spectral
intervals, Gauss points, gas set, k-data, two-stream solver all untouched —
bit-for-bit identical results required, verified by `run_regression.py`). Ordered
by expected payoff.

> **Scheduling decision (2026-06-17):** these kernel-internal optimizations
> (Opp. 1–6) are **deferred until after the external driver/structure refactor is
> complete** (the `REFACTOR_PLAN.md` stages — bundle consolidation, `iso_c_binding`
> library + Python binding, multi-column/OpenMP, IDL→Python, cleanup). Rationale:
> the refactor restructures the very surface these touch (the inner loop, column
> state, and `aerad_driver`/`calc_opd_mod` call path); optimizing it first would
> create churn and force re-validation. Revisit Opp. 1–6 as a dedicated
> performance pass once the structure is stable. Fixes 1–2 (binary search, k-array
> deallocation) and the grid expansion are tracked separately above.

### Opp. 1 — Hoist the P/T interpolation factors out of the band loop *(highest value)*

**Observation.** For a given layer `ik`, the reference indices `p_ref_index`,
`t_ref_index` and the interpolation factors `press`, `temp` (and `onemp`,
`onemt`) depend **only on the layer's P and T — not on the band `sp` or the
g-point.** Yet `bilinear_interpK_grey` and `bilinear_interpK_8gpt_major_gptvec`
**recompute `press` and `temp` from scratch on every call**
(`rad_interp_mod.F90:103–125, 220–242`): the edge tests, the two divisions, the
`1.-press` / `1.-temp`. With 8 grey calls + 1 major call per band × 84 bands ×
~300 layers, that is **~756 redundant recomputations of the same two factors per
layer**, each with a divide.

**Fix (procedural, results-identical).** Compute `press`, `temp`, `onemp`,
`onemt` **once per layer** (right after the index search at line 267) and pass
them into the interp routines instead of `pressure`/`temperature`. The interp
routines drop their edge-test + divide blocks and just do the 4-corner blend. The
arithmetic that produces `ans` is unchanged → bit-for-bit identical. Two divides
per layer instead of ~1500.

**Risk:** Low. Pure hoist; touches the two interp signatures + their callers in
`calc_opd_mod.F90`. Mirror to `3dmodels/src.cam.*` (synced files).

### Opp. 2 — Collapse the 8 grey calls into one call returning a vector

**Observation.** Lines 279–286 make 8 separate `bilinear_interpK_grey` calls
(one per gas), each re-deriving the same `press`/`temp`, each fetching 4 corners
from `k_grey_data(igas, sp, …)`. The 8 gases are the **leading dimension** of
`k_grey_data` (`nspecies, ntot_wavlnrng, kc_npress, kc_ntemp`) — so the 8 needed
values at a given `(sp, p, t)` corner are **strided far apart** in memory
(stride = `ntot_wavlnrng*kc_npress*kc_ntemp`), the worst cache layout.

**Fix.** A single routine that loops the 8 gases internally and returns
`ans_grey(nspecies)`, combined with Opp. 1 (factors passed in). Removes 7×
call overhead and 7× factor recompute per band. **Bigger win if paired with a
data-layout change** (see Opp. 4) so the 8 gas values at a corner are contiguous.
Results identical (same blend, same inputs).

**Risk:** Low–medium. Changes a call signature and an inner loop; covered by
regression Δ=0.

### Opp. 3 — Replace the per-CIA-species linear T-search with binary search (and hoist)

**Observation.** Beyond the main P/T search (Fix 1), the layer loop repeats the
**same downward-from-max linear T-search** for *each* continuum/CIA species —
MT_CKD (lines 313–329), N2-N2 (348–364), H2-H2 (375+), and the others (CO2-CO2
sw/lw, CO2-H2, CO2-CH4, O2-O2, O2-N2, O2-CO2). Each is an independent linear scan
over its own small T-grid, once per layer. These grids are small today, but the
pattern is the same scaling liability Fix 1 addresses, repeated ~10×.

**Fix.** Reuse the `find_ref_index` helper from Fix 1 for every CIA/continuum
T-lookup. Since all use the **same layer `temperature`**, the index into each grid
can be found once per layer per grid (not recomputed). Low absolute cost today,
but cheap to fix while touching this code and it removes the high-T scaling
liability everywhere, not just the main grid.

**Risk:** Low. Same contract/verification as Fix 1.

### Opp. 4 — Reorder k-data array dimensions for cache-friendly inner access *(largest potential, most invasive)*

**This is the highest-ceiling optimization and the only one that attacks the
likely *true* bottleneck (memory latency). It also has the lowest confidence and
is the most invasive — it must be measured, not assumed.**

#### Background: Fortran is column-major

The **leftmost** array index varies fastest in memory. For the current
declarations (`radgrid.F90:366, 369`):

```fortran
k_major_data(nspecies, ntot_wavlnrng, ngauss_8gpt, kc_npress, kc_ntemp)  ! (8,84,8,61,17)
k_grey_data (nspecies, ntot_wavlnrng,              kc_npress, kc_ntemp)  ! (8,84,  61,17)
```

For `k_grey_data` the per-index memory stride is:

| Index | Position | Stride (elements) | Stride (bytes) |
|-------|----------|-------------------|----------------|
| `igas` | leftmost (fastest) | 1 | 8 |
| `sp`   |          | 8 | 64 |
| `p`    |          | 672 | 5,376 |
| `t`    | rightmost (slowest) | 40,992 | ~328 KB |

#### The access pattern (the problem)

For a **fixed layer**, `p_ref_index` and `t_ref_index` are **constant across the
entire band loop**. The loop varies `sp` (84) and `igas` (8) — which are the two
*leftmost, smallest-stride* indices, so that part is fine. The problem is the
**4 bilinear corners**, which differ in `p` and `t`
(`rad_interp_mod.F90:138–141, 253–256`):

```fortran
kdata(igas, sp, p_ref_index,   t_ref_index)
kdata(igas, sp, p_ref_index+1, t_ref_index)     ! 672 elem (5,376 B) from corner 1
kdata(igas, sp, p_ref_index,   t_ref_index+1)   ! 40,992 elem (~328 KB) from corner 1
kdata(igas, sp, p_ref_index+1, t_ref_index+1)
```

A cache line is 64 B (8 doubles). So **each of the 4 corners lands on a separate
cache line**, and the `t` vs `t+1` pair are ~328 KB apart — far beyond any
hardware prefetcher's reach. Every bilinear interp touches **4 unpredictable
cache lines.**

#### Why this dominates runtime

The bilinear blend is ~7 flops (4 mults, 3 adds) — a few cycles. Fetching 4
cache lines that miss is ~4 × 80–300 cycles of DRAM latency. **The arithmetic is
trivial next to the memory access** — the textbook signature of a
**memory-latency-bound** kernel: the core stalls waiting on the corner loads
rather than computing. This fetch happens ~400,000×/column. It is also the reason
*removing arithmetic (Opp. 1) helps little if the kernel is latency-bound* — the
CPU was stalled on memory anyway — whereas fixing the layout attacks the actual
stall.

#### The fix: make the 4 corners adjacent

Move `p` and `t` to the **leftmost (fastest)** positions so that, for a fixed
`(igas, sp)`, the 2×2 P-T corner block is contiguous:

```fortran
! current
k_grey_data(nspecies, ntot_wavlnrng, kc_npress, kc_ntemp)   ! p,t SLOWEST
!            igas      sp             p          t

! proposed
k_grey_data(kc_npress, kc_ntemp, nspecies, ntot_wavlnrng)   ! p,t FASTEST
!            p          t          igas      sp
```

With the proposed layout the 4 corners for a fixed `(igas, sp)`:

- `(p, t)`, `(p+1, t)` → **1 element apart** (same cache line)
- `(p, t+1)`, `(p+1, t+1)` → `kc_npress` (61) elements away (~488 B, ~8 lines)

The 4 corners collapse from **4 scattered lines (one ~328 KB away)** to **~2
nearby lines**, and consecutive `(igas, sp)` iterations now stride predictably so
the prefetcher can help.

For `k_major_data`, include `ngauss_8gpt` in the reorder. The 8 g-points for a
fixed band are read together in `bilinear_interpK_8gpt_major_gptvec`, so a layout
like `(kc_npress, kc_ntemp, ngauss_8gpt, nspecies, ntot_wavlnrng)` keeps **both**
the g-point loop and the P-T corners local.

#### Free at runtime, bit-identical

The reorder is done **once at init**, when `model_specific.F90` packs the
individual `k_*` arrays into `k_major_data`/`k_grey_data` (lines 54–61, 97–104) —
the copy becomes a transposing copy. The timestep loop reads the better layout
and pays **nothing**. Results are **bit-for-bit identical** because only *storage
order* changes, not any value or computation → passes regression at Δ=0 while
(potentially) running much faster.

#### Honest expectations & what can deflate the win

Cache-layout fixes on a latency-bound inner kernel can deliver **1.5–3× on the
affected routine**, occasionally more. Two things can mute it:

- **Working set may already fit in L2.** The full arrays (56 MB grey, 446 MB
  major) far exceed cache, but the *actively touched* slice for one layer (one
  `(p,t)` neighborhood across 84 bands × 8 gases) may be only a few hundred KB and
  partly L2-resident already.
- **The core may already hide the latency** via out-of-order execution + prefetch
  across the band loop.

Only profiling tells you. **Measure with cache-miss counters, not just
wall-clock:** `perf stat -e L1-dcache-load-misses,LLC-load-misses` before/after to
confirm the memory behavior actually moved.

#### Files touched (all must agree on the new order)

1. `radgrid.F90:366, 369, 372` — array declarations.
2. `model_specific.F90:54–104` — packing becomes a transposing copy.
3. `rad_interp_mod.F90` — `kdata(...)` index order in both interp routines
   (lines 138–141, 253–256).
4. `initialize_rad_mod_1D.F90` / `_cam.F90` — readers fill `k_h2o` etc. in the
   NetCDF native order; cleanest to keep the **individual** arrays in file-order
   and transpose **only** when packing into `k_major_data`/`k_grey_data` (step 2),
   so the readers and on-disk data are untouched.
5. `mpibcast` counts (`_cam.F90:242+`) — total element count unchanged, so the
   broadcast still works; only the in-memory layout differs.
6. Mirror all of the above to `3dmodels/src.cam.*` per the CLAUDE.md sync rule.

**Risk:** Medium–high (most invasive Opp.). Validate Δ=0 *and* profile
before/after — layout wins are empirical.

**Sequencing note.** Best done as the **capstone**, after Opp. 1 and Opp. 2:
hoisting (Opp. 1) makes the kernel *more* memory-bound, which **raises** this
payoff; vectorizing the grey call (Opp. 2) lets all 8 gases' corners be fetched in
one contiguous sweep if `igas` is laid out adjacently. Measure at each step.

### Opp. 5 — Skip absent gases / zero columns

**Observation.** `tau_grey` is computed for all 8 gases every band even when a
gas's column amount `ugas(igas)` is **zero** (common: many runs have no O3, no
CO, no C2H6). The grey interp still runs, then multiplies by zero (lines
287–294). The CIA blocks **already** guard on `if (u_n2 > 0)` etc. (line 347) —
the major/grey path does not.

**Fix.** Skip the grey interp for any gas with `ugas(igas) == 0` (set
`tau_grey(igas)=0` directly). Since a zero-column gas can never be the `maxloc`
major absorber and contributes zero to the minor sum, results are identical. Cuts
interp calls proportional to the number of absent gases — substantial for typical
sparse gas mixes.

**Risk:** Very low. Pure guard; Δ=0 by construction.

### Opp. 6 — `maxloc` on an 8-element array every band

**Observation.** `imajor = maxloc(tau_grey)` (line 295) runs once per band per
layer over 8 elements. `maxloc` returns a rank-1 array and has overhead beyond a
hand-rolled scan. Minor, but it is in the innermost loop.

**Fix.** Replace with a simple inline max-index scan over the (possibly reduced,
post-Opp.5) gas set. Marginal; do only if profiling flags it.

### What was checked and is **not** worth changing

- **The two-stream solver / Planck / MCICA** — physics-bearing; out of scope by
  the hard constraint.
- **The bilinear blend arithmetic itself** — already minimal (4 mults + 3 adds);
  nothing to gain without changing results.
- **Band loop count, g-point count, gas set** — fixed by the physics constraint.

### Recommended order & method

1. **Profile first** (`gprof`/`perf` on a representative 1-D run, or CAM timers in
   3-D) to confirm the hot path is `calc_opd_mod` interpolation as expected — do
   not optimize on assumption.
2. Land the **low-risk, results-identical** wins first: Opp. 5 (skip zero gases),
   Opp. 1 (hoist factors), Opp. 3 (CIA binary search). Each gated at **Δ=0** on
   the 13 regression cases.
3. Then Opp. 2 (vector grey call), then evaluate **Opp. 4 (layout)** with
   before/after profiling — it is the biggest potential 3-D win but the most
   invasive and must be measured, not assumed.
4. All edits to files shared with `3dmodels/src.cam.*` must be mirrored and
   `diff`-clean per the CLAUDE.md sync rule.

---

## Expected speedup — honest estimates

**Caveat up front: these are reasoned estimates, not measurements.** Reliable
percentages require profiling (see "Recommended order", step 1). The numbers below
are bounded by call-count arithmetic and CPU behavior, and flagged where
data-dependent or bottleneck-dependent.

### Hot-path arithmetic (per column, per radiation step)

- Layers ~300 (`pverp`); bands 84; per band 8 grey interps + 1 major interp (8 g-pts).
- Per layer: 84 × (8 grey + 1 major) ≈ **1,344 interp evaluations**.
- Per column: ~300 × 1,344 ≈ **~400,000 interp evaluations.**

This ~400k is the baseline every Opp. is measured against.

### Per-Opp. expectation

| Opp. | Mechanism | Estimated total-RT speedup | Confidence |
|------|-----------|----------------------------|------------|
| **5 — skip zero gases** | Skips grey interp for absent gases (grey = 8 of 9 calls). **Data-dependent.** | **0% (all 8 gases active) → ~20%+ (sparse mix, e.g. 3–4 active)** | Med — depends on typical gas loadout |
| **1 — hoist P/T factors** | Removes ~750 divides/layer (~225k/column → ~600). Real work removed, but smaller if memory-bound. | **~2% → ~8%** | Med |
| **4 — data layout** | Targets memory bandwidth / cache (strided access defeats prefetch). | **~0% → could be the routine's single biggest factor (1.5–3× *on the routine*)** | Low — must be measured |
| **2 — vectorize grey** | Removes 7× call overhead; compiler may already inline most. Main value = vehicle for Opp. 4. | **low single digits standalone** | Low–Med |
| **3 — CIA binary search** | CIA T-grids are tiny; value is high-T scaling guarantee, not speed. | **<1%** | High (small) |
| **6 — replace `maxloc`** | Innermost-loop micro-opt. | **noise** | High (small) |

### Aggregate (typical sparse-gas run, pre-profiling)

| | Pessimistic | Optimistic |
|---|-------------|------------|
| Opp. 5 (sparse gases) | ~5% | ~20% |
| Opp. 1 (hoist) | ~2% | ~8% |
| Opp. 4 (layout) | ~0% | routine's biggest single factor |
| Opp. 2/3/6 | ~0% | ~3% |

**Plausible combined: ~10–30% off total RT runtime, with Opp. 4 a wildcard that
could push higher — or contribute little.** Whole-model (CAM) speedup is smaller,
since RT is only part of a timestep.

### Two caveats that shape the above

1. **They do not simply add.** Removing the divides (Opp. 1) makes the routine
   *more* memory-bound, which **raises** the relative payoff of Opp. 4 and
   **lowers** the marginal benefit of further compute trimming. They interact.
2. **The premise needs verification.** All of this assumes `calc_opd_mod`
   interpolation is the dominant cost. In a coupled CAM run the two-stream solver,
   MCICA, or I/O may dominate — in which case even clean interp optimizations move
   the total little. **This is why step 1 is "profile first."**

### Cheapest path to real numbers

Build with `-pg`, run one representative 1-D profile, `gprof` it. A single
measurement replaces every estimate above with fact and tells you whether Opp. 4
justifies its invasiveness.

---

## Quick reference — files involved

| File | Role |
|------|------|
| `source/src.exort/calc_opd_mod.F90:231–267` | bracketing search (Fix 1 target) |
| `source/src.exort/rad_interp_mod.F90:46–273` | bilinear interp routines (grid-agnostic; read `tgrid`/`log10pgrid` directly) |
| `source/src.exort/radgrid.F90:216–217` | `kc_npress=61`, `kc_ntemp=17` |
| `source/src.exort/radgrid.F90:366–372` | `k_major_data`, `k_grey_data`, individual `k_*` arrays (Fix 2 target: line 372) |
| `source/src.exort/model_specific.F90:54–61, 97–104` | packs individual k arrays into major/grey (init-only) |
| `source/src.exort/initialize_rad_mod_1D.F90` / `_cam.F90` | k-file readers (Fix 2: allocate/deallocate here) |
| `tests/regression/run_regression.py` | Δ=0 verification / rebaseline |
