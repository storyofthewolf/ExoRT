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
  haze kernel's apparent off-by-one (layer ik−1 mass paired with pdel(ik),
  bottom rad level never filled) was deliberately not reproduced — the 1-D
  kernel matches the driver's coldens convention. *Undo:* `git revert cc892f1`.
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
- **Haze 3-D port** — reconcile the 1-D haze path with
  `3dmodels/src.cam.n68equiv.haze` (the CARMA-module imports live there only);
  includes deciding whether to adopt the 1-D kernel's layer-indexing fix in the
  3-D kernel (see the off-by-one note in Stage C3).
- **Clear-sky / cloud-forcing double-run** (`_CLD` outputs) — deferred to its own
  stage; the experimental `plotspectra_1D.pro` expects it.
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
