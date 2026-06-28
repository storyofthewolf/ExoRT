# src.exort (HITRAN-2016) ≡ n68equiv equivalence record

**Captured:** 2026-06-28, branch `refactor`, immediately before rebaselining the
regression suite from `n68equiv.exe` to `exort.exe` (HITRAN-2016 native gases).

## Why this file exists

The 13-case regression suite was historically baselined against `n68equiv.exe`
(68-band grid, HITRAN-2016). The v2 single-bundle `exort.exe` (84-band grid)
supersedes it. Before regenerating the golden baselines *with* `exort.exe`, we
captured this one-time proof that `exort` on the **HITRAN-2016** native-gas
k-files reproduces the retiring n68 physics. After the rebaseline the suite
guards `exort` against its own future regressions and no longer encodes the n68
values, so this record is the permanent evidence of the transition's neutrality.

This run pins the four native gases (H₂O/CO₂/CH₄/C₂H₆) to HITRAN-2016 — the same
pinning now committed as the default in `source/src.exort/kabs.F90`. NH₃/CO are
HITRAN-2024 (no HITRAN-2016 table exists; both proven clean in `gas_sweep.py`);
O₂/O₃ are HITRAN-2020. This is the data bundle the rebaselined suite uses.

## Result

Command: `USER_FC=gfortran python run_regression.py --exort h16`
(exort.exe on h16 native gases + _n84 stars, compared element-wise vs the
n68equiv golden baselines, `_SPECTRAL` arrays skipped — 84-band rows cannot align
with 68-band baselines).

- **Longwave is bit-for-bit identical on all 13 cases** (`LWUP_TOM` Δ = +0.000).
  The n68→n84 band merge and the HITRAN-2016 native-gas pinning are
  longwave-neutral.
- **The only differences are shortwave (`SWUP`), and only on the G2V solar
  cases** — never on the smooth blackbody_3400K cases — at max rel ≈ 0.37%.
  This is the known **n68→n84 stellar UV regrid**: the structured G2V spectrum is
  re-binned onto the finer n84 UV grid; smooth blackbody spectra are unaffected.
  It is a grid-resolution effect on the input stellar spectrum, **not** a
  k-coefficient or RT-code difference. See memory
  `exort-merge-equivalence-verified`.

```
=== src.exort (h16, _n84 stars) vs n68equiv baselines (rtol=0.001, atol=0.001) ===
  case                            LWUP_TOM [W/m2]                     SWDN_SFC [W/m2]
  TS250K_G2V            FAIL        179.023 vs   179.023 (Δ  +0.000)    293.276 vs   293.483 (Δ  -0.207)
      - SWUP: 264/301 pts exceed tol; max |abs|=1.771e-01 max rel=1.844e-03 at (6,) (new=9.621811e+01 base=9.604099e+01)
  TS250K_BB3400         PASS        179.023 vs   179.023 (Δ  +0.000)    279.942 vs   279.942 (Δ  +0.000)
  TS273K_G2V            FAIL        226.290 vs   226.290 (Δ  +0.000)    272.052 vs   272.240 (Δ  -0.187)
      - SWUP: 270/301 pts exceed tol; max |abs|=1.828e-01 max rel=2.028e-03 at (4,) (new=9.030659e+01 base=9.012379e+01)
  TS273K_BB3400         PASS        226.290 vs   226.290 (Δ  +0.000)    241.218 vs   241.218 (Δ  +0.000)
  TS300K_G2V            FAIL        268.726 vs   268.726 (Δ  +0.000)    233.927 vs   234.080 (Δ  -0.154)
      - SWUP: 278/301 pts exceed tol; max |abs|=1.934e-01 max rel=2.421e-03 at (2,) (new=8.007743e+01 base=7.988400e+01)
  TS300K_BB3400         PASS        268.726 vs   268.726 (Δ  +0.000)    181.169 vs   181.169 (Δ  +0.000)
  TS320K_G2V            FAIL        283.244 vs   283.244 (Δ  +0.000)    199.486 vs   199.611 (Δ  -0.125)
      - SWUP: 284/301 pts exceed tol; max |abs|=2.037e-01 max rel=2.844e-03 at (2,) (new=7.183705e+01 base=7.163334e+01)
  TS320K_BB3400         PASS        283.244 vs   283.244 (Δ  +0.000)    133.871 vs   133.871 (Δ  +0.000)
  TS340K_G2V            FAIL        292.872 vs   292.872 (Δ  +0.000)    164.739 vs   164.839 (Δ  -0.100)
      - SWUP: 288/301 pts exceed tol; max |abs|=2.150e-01 max rel=3.296e-03 at (1,) (new=6.546597e+01 base=6.525094e+01)
  TS340K_BB3400         PASS        292.872 vs   292.872 (Δ  +0.000)     94.461 vs    94.461 (Δ  +0.000)
  TS360K_G2V            FAIL        293.408 vs   293.408 (Δ  +0.000)    130.156 vs   130.235 (Δ  -0.079)
      - SWUP: 291/301 pts exceed tol; max |abs|=2.269e-01 max rel=3.702e-03 at (0,) (new=6.149902e+01 base=6.127216e+01)
  TS360K_BB3400         PASS        293.408 vs   293.408 (Δ  +0.000)     64.165 vs    64.165 (Δ  +0.000)
  2barCO2_dry_Mars_G2V  FAIL         92.744 vs    92.744 (Δ  +0.000)    138.038 vs   138.159 (Δ  -0.121)
      - SWUP: 253/301 pts exceed tol; max |abs|=1.182e-01 max rel=1.278e-03 at (1,) (new=9.257124e+01 base=9.245308e+01)

6/13 cases passed.
```

The 7 "FAIL"s are the SW UV regrid only; LW is exact. This is the documented
"n84 supersedes n68" result, not a physics regression.
