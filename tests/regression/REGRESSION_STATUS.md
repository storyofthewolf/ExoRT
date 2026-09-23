# Regression status — current green record

**All gates pass.** This is the live record of the test suite in its passing
state: what was run, on what, and what came out. Regenerate it whenever the
baselines are regenerated, so it never disagrees with `baselines/`.

| | |
|---|---|
| **Commit** | HITRAN-2024 default rebaseline (see `git log -- tests/regression/baselines`) |
| **Branch** | `refactor` |
| **Captured** | 2026-09-23 |
| **Compiler** | GNU Fortran (Homebrew GCC 16.1.0) 16.1.0 |
| **NetCDF Fortran** | 4.6.3 |
| **Platform** | macOS 15.7.9, arm64 |
| **Line list** | HITRAN-2024 native gases (H₂O/CO₂/CH₄/C₂H₆) and NH₃/CO; O₂/O₃ HITRAN-2020. HITRAN-2016 via `--exort h16` |
| **Haze optics** | `data/aerosol/haze_n84_b40_mie.nc` |

> Looking for `EXORT_H16_N68vN84_GRID.md` and wondering why it shows failures?
> That file is a **historical** record from 2026-06-28 comparing the 84-band
> `exort.exe` against the retired 68-band `n68equiv` baselines. Its "FAIL"s are
> the expected stellar-UV regrid difference, documented in that file. **This**
> file is the current state.

---

## 1. Regression suite — 16/16, Δ = 0

```
cd tests/regression && python run_regression.py
```

```
=== regression results: exort.exe (h24) (rtol=0.001, atol=0.001) ===
  case                                 LWUP_TOM [W/m2]                     SWDN_SFC [W/m2]                   
  TS250K_G2V                 PASS        177.595 vs   177.595 (Δ  +0.000)    292.857 vs   292.857 (Δ  +0.000)
  TS250K_BB3400              PASS        177.595 vs   177.595 (Δ  +0.000)    278.939 vs   278.939 (Δ  +0.000)
  TS273K_G2V                 PASS        223.354 vs   223.354 (Δ  +0.000)    271.339 vs   271.339 (Δ  +0.000)
  TS273K_BB3400              PASS        223.354 vs   223.354 (Δ  +0.000)    239.651 vs   239.651 (Δ  +0.000)
  TS300K_G2V                 PASS        264.242 vs   264.242 (Δ  +0.000)    232.653 vs   232.653 (Δ  +0.000)
  TS300K_BB3400              PASS        264.242 vs   264.242 (Δ  +0.000)    178.851 vs   178.851 (Δ  +0.000)
  TS320K_G2V                 PASS        279.050 vs   279.050 (Δ  +0.000)    197.965 vs   197.965 (Δ  +0.000)
  TS320K_BB3400              PASS        279.050 vs   279.050 (Δ  +0.000)    131.674 vs   131.674 (Δ  +0.000)
  TS340K_G2V                 PASS        289.444 vs   289.444 (Δ  +0.000)    162.742 vs   162.742 (Δ  +0.000)
  TS340K_BB3400              PASS        289.444 vs   289.444 (Δ  +0.000)     92.157 vs    92.157 (Δ  +0.000)
  TS360K_G2V                 PASS        290.370 vs   290.370 (Δ  +0.000)    127.705 vs   127.705 (Δ  +0.000)
  TS360K_BB3400              PASS        290.370 vs   290.370 (Δ  +0.000)     62.115 vs    62.115 (Δ  +0.000)
  2barCO2_dry_Mars_G2V       PASS         92.937 vs    92.937 (Δ  +0.000)    138.011 vs   138.011 (Δ  +0.000)
  2barCO2_co2cloud_Mars_G2V  PASS         33.955 vs    33.955 (Δ  +0.000)     71.891 vs    71.891 (Δ  +0.000)
  TS300K_haze_G2V            PASS        263.000 vs   263.000 (Δ  +0.000)    194.278 vs   194.278 (Δ  +0.000)
  TS300K_hazethick_G2V       PASS        244.676 vs   244.676 (Δ  +0.000)     63.514 vs    63.514 (Δ  +0.000)

16/16 cases passed.
```

12 clear TS × {G2V, blackbody_3400K}, 1 clear Mars, 1 CO₂-ice cloud Mars, and
2 haze loadings (thin τ≈0.5, thick τ≈8.6). All `pver=300`.

## 2. Multi-column — batch ≡ singles, threaded ≡ serial

```
cd tests/regression && python multicol_check.py
```

```
== batch vs singles (exact, max |delta| must be 0) ==
  col 0 (TS250K_G2V): max|d|=0.000e+00  PASS
  col 1 (TS300K_G2V): max|d|=0.000e+00  PASS
  col 2 (TS340K_G2V): max|d|=0.000e+00  PASS

== batch OMP_NUM_THREADS=8 vs 1 (exact) ==
  threaded vs serial: max|d|=0.000e+00  PASS
PASS
```

## 3. Per-column gravity / insolation

```
cd tests/regression && python percol_config_check.py
```

```
== batch (per-column config) vs singles (namelist config) ==
  col 0 (default): max|d|=0.000e+00  PASS
  col 1 (marsG): max|d|=0.000e+00  PASS
  col 2 (marsGscon): max rel=0.000e+00  PASS
PASS
```

Per-column gravity is exact; per-column stellar constant differs only by
floating-point rounding order (gate: 1e-10).

## 4. Library API — Python binding

```
python tools/exort_pytools/verify_lib.py
```

```
=== determinism (exact, max |delta| must be 0) ===
  repeat run_column    PASS  max|d|=0.000e+00
  run_columns[0]       PASS  max|d|=0.000e+00
  run_columns[1]       PASS  max|d|=0.000e+00

LWUP_TOM=264.186 W/m2   SWDN_SFC=232.701 W/m2
PASS
```

Residuals vs the committed baselines are ~6e-8 relative — the float32 storage
precision of the baseline files, not a solver difference.

## 5. Library API — C harness

```
cd tests/lib && make run
```

```
=== batch determinism (exact) ===
  run_columns[0]           PASS  max|d|=0.000e+00
  run_columns[1]           PASS  max|d|=0.000e+00

LWUP_TOM=264.186 W/m2  SWDN_SFC=232.701 W/m2
PASS
```

## 6. 3-D bundle sync

```
python tools/populate3Dmodels.py check --exort
```

```
PASS: all selected bundles match source.
```

`3dmodels/src.cam.exort` is byte-identical to `source/` for all shared files.
The pre-v2 bundles intentionally diverge and are not checked.

## 7. CAM compile gate

```
bash tests/cam_compile_check/run_all.sh
```

```
OK: src.cam.exort compiled with flags:
OK: src.cam.exort compiled with flags: -DEXORT_CO2CLD
OK: src.cam.exort compiled with flags: -DEXORT_CARMA
OK: src.cam.exort compiled with flags: -DEXORT_SRF_EMISS
OK: src.cam.exort compiled with flags: -DEXORT_CO2CLD -DEXORT_CARMA -DEXORT_SRF_EMISS
ALL OK
```

Compiles the bundle against stubbed CESM modules in all macro combinations.
**This is a compile gate only** — `src.cam.exort` has never had a real ExoCAM
build or run.

---

## Reproducing this

```bash
cd build && USER_FC=gfortran make exort     # gfortran needed on arm64
cd ../tests/regression && python run_regression.py
```

If the link fails with `ld: library not found for -lnetcdf`, your
`nf-config --flibs` names the C library without giving its directory:

```bash
env LIBRARY_PATH=$(nc-config --libdir) USER_FC=gfortran make exort
```

## Keeping this file honest

Regenerate it in the **same commit** that regenerates baselines. A physics
change that moves numbers should update the transcripts here alongside
`baselines/*.nc`, so the two can never disagree.

If you change code and the suite still shows 16/16 Δ=0, the change is
bit-for-bit neutral. If it moves numbers, that is either a bug or an intended
physics change — and an intended one gets its baselines regenerated in the same
commit, called out explicitly in the message. See `REFACTOR_PLAN.md`.
