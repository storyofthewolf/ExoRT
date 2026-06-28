# Stage C integration map — fold experimental CO₂ clouds + haze into `src.exort`

**Status:** survey / map only. No integration code written yet. Drafted
2026-06-28 on `refactor` at the user's request ("map it first, don't merge yet").

Source threads live in `source/experimental/src.n68equiv_exp/`:
- `shr/`  — **CO₂ ice clouds** (the stated live priority: "troubleshoot CO₂ ice
  cloud radiation to get the Mars model running")
- `shr2/` — **CARMA haze** hookup (built *on top of* the cloud thread)

Both are near-complete forks of the whole RT bundle (~14 files each), built on
the **n68equiv / HITRAN-2016 / 68-band** base — they **predate** the n84/exort
merge entirely. So this is a *rebase onto 84-band + current driver*, not a copy.

---

## 1. The `.F90~` backups — RESOLVED

`shr2/` carried three editor backups: `calc_opd_mod.F90~`,
`exo_radiation_mod.F90~`, `exo_radiation_cam_intr.F90~`.

**Finding: the live `.F90` files are NEWER and strictly supersede the `~`
backups.** Diffing each `~` (pre-edit) against its live `.F90`:

- `exo_radiation_mod.F90`: live **adds** the CARMA haze hookup the backup lacks —
  `use carma_model_mod (NELEM,NBIN)`, `use carma_exort_mod`, a new
  `ext_carmammr(pver,nelem,nbin)` arg, `bext/tau/taueff_aer_carma` arrays, the
  `if(do_carma_rt) call calc_aeropd(...)` hook, and CARMA optics threaded into
  the two-stream solver.
- `exo_radiation_cam_intr.F90`: live adds `do_carma_exort`,
  `carma_exort_optics_init`, `carma_exort_get_mmr`, and `carmammr`/`carmammr_zero`
  arrays threaded into both clearsky and full-sky `aerad_driver` calls.
- `calc_opd_mod.F90`: live has a `kc_ntemp` bugfix + one comment the backup
  lacks. Trivial; live is correct.

**Disposition: discard the three `.F90~` files.** They contain nothing the live
versions lack. This clears the REFACTOR_PLAN gate ("`.F90~` backups
resolved/committed before merge"). *(Recommend deleting them in the first
Stage-C-prep commit so the tree is clean before any merge work.)*

---

## 2. Architecture: RT is self-contained; formation physics lives outside

**Maintainer-confirmed design boundary (2026-06-28):** the cloud/haze *formation*
physics (computing CO₂-ice condensate; CARMA aerosol microphysics) lives
**outside** the RT code. ExoRT's job is purely **radiative** — it takes the
condensate/bin arrays as *given inputs* and computes their optical contribution.
It never computes them.

The `carma_exort_get_mmr` / `pbuf_get_index` calls seen in
`exo_radiation_cam_intr.F90` are the **3-D CAM glue** that *fetches* those arrays
from the CESM state. In the **1-D path they are not needed**: the arrays enter
directly from `RTprofile_in.nc`:

- CO₂-ice cloud: `cicewp_co2`, `rei_co2` (per-level) — optional input vars.
- Haze: `carmammr(pver, nelem, nbin)` — optional input var (binwise aerosol MMR).

**Implication (corrects an earlier draft):** **both** CO₂ clouds *and* haze are
1-D-testable, because both consume input arrays rather than CARMA modules. The
1-D RT path must therefore read `carmammr` from the deck and **sever the
`carma_model_mod`/`carma_exort_mod` coupling** (those stay in the CAM-only
`exo_radiation_cam_intr.F90`). `nelem`/`nbin` become small parameters (or
input-file dims) on the 1-D side, not CARMA-module imports.

Haze is a **mature, published** model (maintainer), so it is a committed Stage C
goal — not speculative. It is sequenced after CO₂ clouds only because its 84-band
optics data still needs regenerating (§4), not because of uncertainty.

## 2b. `src.exort` already has the scaffolding

`src.exort` is not greenfield for this — the stubs exist:
- `calc_opd_mod.F90` already defines `calc_cldopd` (H₂O clouds) **and an empty
  `calc_aeropd( )` stub** — the haze hook is a shell waiting to be filled.
- `cloud.F90` already loads `cloudoptics_h2o_{liquid,ice}_mie_n84.nc` — adding a
  CO₂-ice path follows the identical parameter pattern.
- `aerad_driver` (`src.main/exo_radiation_mod.F90:154`) and its `main.F90:86` call
  site are the known seams; new cloud/haze args go here as **`optional`**.

---

## 3. Where the physics actually lives (the real integration surface)

Whole-file diff counts are inflated by the 68→84 grid + h16/h24 base difference
and are **not** a measure of cloud/haze complexity. By keyword, the actual hooks:

### CO₂ clouds (`shr/`)
| File | What it adds |
|------|--------------|
| `calc_opd_mod.F90` | `subroutine calc_cldopd_co2(ext_cICE_co2, ext_REI_CO2, tau_cld_co2, singscat_cld_co2, asym_cld_co2)` — the CO₂-ice cloud optical-depth kernel. ~12 cloud-hook lines. |
| `cloud.F90` | CO₂-ice cloud optical-property file paths/data (only ~6 lines differ from exort — small). |
| `exo_radiation_mod.F90` | threads cloud optics into the two-stream solver; clear-sky vs cloudy double-run scaffolding. |
| `exo_radiation_cam_intr.F90` | 3-D only: pbuf reads, `do_exo_condense_co2`, clearsky forcing diagnostics (QRLC/FLDSC/etc.). |
| `mcica.F90` | cloud-overlap subcolumn generator (the "isCloudy" logic) — needed for fractional cloud cover. ~2 lines vs current. |

### Haze (`shr2/`, on top of `shr/`)
| File | What it adds |
|------|--------------|
| `calc_opd_mod.F90` | `calc_aeropd(...)` CARMA aerosol optical-depth kernel; ~25 haze-hook lines. |
| `exo_radiation_mod.F90` | `ext_carmammr` arg + CARMA optics arrays + `if(do_carma_rt) call calc_aeropd` (the live-vs-backup delta above). |
| `exo_radiation_cam_intr.F90` | `do_carma_exort`, `carma_exort_*` calls, binwise `carmammr` threading. All 3-D. |

---

## 4. The 68→84 optics-data state — RESOLVED

The cloud/haze kernels read **pre-tabulated, band-specific Mie optics** at init;
the kernel bodies are otherwise grid-agnostic (loop `ntot_wavlnrng`, scale by
particle size/density). So only the *data file* is band-keyed. State on disk:

| Feature | 84-band file | Status |
|---------|--------------|--------|
| **CO₂ ice (C1)** | `data/cloud/cloudoptics_co2_ice_mie_n84.nc` | ✅ **exists** — C1 unblocked |
| **Haze (C3)** | `data/aerosol/haze_n84_*.nc` | ❌ **only `haze_n68_b40_{mie,fractal_interp}.nc`** — needs 84-band regen offline before C3 can validate |

`calc_cldopd_co2` interpolates `Qcldice_co2/Wcldice_co2/Gcldice_co2` (loaded from
the `cloudoptics_co2_ice_mie_n84.nc` file) via `interpolate_cld_co2`, then scales
`tau = 3·Q/(4·ρ_ice·r)·cICE`. Confirmed file-based + 84-band data present.

**Net:** CO₂ clouds (C1) have everything they need. Haze (C3) has the kernel but
needs the 84-band aerosol optics regenerated offline (the one outstanding
artifact — analogous to the HITRAN-2024 k-tables, but smaller and uncontroversial
since the haze model itself is published/mature).

---

## 5. Recommended Stage C plan (revised from REFACTOR_PLAN)

Split the monolithic "Stage C" into reviewable sub-stages:

- **C0 (prep):** delete the three stale `shr2/*.F90~` backups; commit. Clears the
  plan gate.
- **C1 (CO₂ clouds, 1-D):** rebase `calc_cldopd_co2` + cloud.F90 optics onto the
  84-band `src.exort`; make the new cloud args to `calc_*opd`/`aerad_driver`
  **`optional`** (per the surviving ARGUMENT_HANDLING idea) so `main.F90` and the
  regression cases pass nothing and stay bit-identical. Add `cicewp_co2`/`rei_co2`
  as optional input-file vars (the established `nf_inq_varid` pattern). **Verify:
  clouds-OFF reproduces the current 13/13 exort-h16 baselines exactly (Δ=0);**
  clouds-ON on a Mars CO₂-ice fixture gives sensible OLR change.
- **C2 (CO₂ clouds, 3-D):** port the `exo_radiation_cam_intr.F90` pbuf path +
  clearsky forcing diagnostics; reconcile with `3dmodels/`. Verify via
  `populate3Dmodels.py diff`.
- **C3 (haze, 1-D):** fill the `calc_aeropd()` stub with the CARMA aerosol
  kernel; read `carmammr(pver,nelem,nbin)` as an optional input var; thread it as
  an `optional` arg through `aerad_driver`/`main.F90`. **Sever the
  `carma_model_mod`/`carma_exort_mod` coupling on the 1-D side** — `nelem`/`nbin`
  are parameters/input dims, the array comes from the deck. **Blocked only on**
  regenerating 84-band haze optics (`data/aerosol/haze_n84_*.nc`; §4) — kernel and
  wiring can be done first, validation waits on the data.
- **C2 / 3-D (either feature):** port the `exo_radiation_cam_intr.F90` pbuf +
  `carma_exort_*` paths; reconcile with `3dmodels/` (incl. `src.cam.n68equiv.haze`)
  via `populate3Dmodels.py diff`. The CARMA-module imports live **here only**.
  Done after the 1-D features land.

**Gating principle (from REFACTOR_PLAN "Reproducibility"):** clouds-OFF / haze-OFF
must be Δ=0 vs current baselines; clouds-ON / haze-ON is an intended physics
change that rebaselines the harness in the same commit, called out explicitly.

---

## Decisions (resolved with maintainer 2026-06-28)

1. **Scope:** C0 → C1 (CO₂ clouds, 1-D) first. ✓
2. **Haze:** mature & published — a committed Stage C goal (C3), not speculative.
   Sequenced after C1 only due to the 84-band optics regen.
3. **Optics grid:** `calc_cldopd_co2` is file-tabulated; the **84-band CO₂-ice
   file already exists** (`cloudoptics_co2_ice_mie_n84.nc`). Haze 84-band optics
   still need regenerating offline (§4).
4. **Inputs:** cloud/haze fields enter via `RTprofile_in.nc` (optional vars,
   `nf_inq_varid` pattern). RT is self-contained — formation physics is external.

## Remaining setup item

- **No CO₂-ice / haze input decks exist** to validate clouds-ON / haze-ON. A
  cloudy fixture must be built — port the cloud-related edits from the
  experimental `makeColumn.pro`/`getColumn.pro` into the Python `makeColumn.py`
  (add `cicewp_co2`/`rei_co2` and optionally `carmammr` writers). Until then, C1
  can only assert the **clouds-OFF Δ=0** half of the gate; the clouds-ON sanity
  check needs this fixture. (Part of C1's deliverable.)
