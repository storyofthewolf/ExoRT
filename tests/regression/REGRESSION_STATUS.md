# Regression status — current green record

**All gates pass.** This is the live record of the test suite in its passing
state: what was run, on what, and what came out. Regenerate it whenever the
baselines are regenerated, so it never disagrees with `baselines/`.

| | |
|---|---|
| **Commit** | parent of the commit that last touched this file (`git log -- tests/regression/REGRESSION_STATUS.md`) |
| **Branch** | `refactor` |
| **Captured** | 2026-09-23 |
| **Compiler** | GNU Fortran (Homebrew GCC 16.1.0) 16.1.0 |
| **NetCDF Fortran** | 4.6.3 |
| **Platform** | macOS 15.7.9, arm64 |
| **Line list** | HITRAN-2024 for H₂O/CO₂/CH₄/C₂H₆/NH₃/CO; O₂/O₃ HITRAN-2020. HITRAN-2016 via `--exort h16` |
| **Haze optics** | `data/aerosol/haze_n84_b40_mie.nc` |

How the baselines got here (n68equiv → exort, HITRAN-2016 → 2024, the H₂O
temperature-index fix) is summarised in `CHANGES.md`.

---

## 1. Regression suite — 32/32, Δ = 0

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
  CH4_realistic_G2V          PASS        308.390 vs   308.390 (Δ  +0.000)    313.875 vs   313.875 (Δ  +0.000)
  CH4_realistic_BB3400       PASS        308.390 vs   308.390 (Δ  +0.000)    328.411 vs   328.411 (Δ  +0.000)
  CH4_elevated_G2V           PASS        280.718 vs   280.718 (Δ  +0.000)    244.325 vs   244.325 (Δ  +0.000)
  CH4_elevated_BB3400        PASS        280.718 vs   280.718 (Δ  +0.000)    190.439 vs   190.439 (Δ  +0.000)
  CO_realistic_G2V           PASS        311.769 vs   311.769 (Δ  +0.000)    315.046 vs   315.046 (Δ  +0.000)
  CO_realistic_BB3400        PASS        311.769 vs   311.769 (Δ  +0.000)    331.995 vs   331.995 (Δ  +0.000)
  CO_elevated_G2V            PASS        310.694 vs   310.694 (Δ  +0.000)    312.183 vs   312.183 (Δ  +0.000)
  CO_elevated_BB3400         PASS        310.694 vs   310.694 (Δ  +0.000)    323.603 vs   323.603 (Δ  +0.000)
  NH3_realistic_G2V          PASS        311.707 vs   311.707 (Δ  +0.000)    315.081 vs   315.081 (Δ  +0.000)
  NH3_realistic_BB3400       PASS        311.707 vs   311.707 (Δ  +0.000)    332.149 vs   332.149 (Δ  +0.000)
  NH3_elevated_G2V           PASS        220.547 vs   220.547 (Δ  +0.000)    291.237 vs   291.237 (Δ  +0.000)
  NH3_elevated_BB3400        PASS        220.547 vs   220.547 (Δ  +0.000)    275.524 vs   275.524 (Δ  +0.000)
  C2H6_realistic_G2V         PASS        311.808 vs   311.808 (Δ  +0.000)    315.080 vs   315.080 (Δ  +0.000)
  C2H6_realistic_BB3400      PASS        311.808 vs   311.808 (Δ  +0.000)    332.144 vs   332.144 (Δ  +0.000)
  C2H6_elevated_G2V          PASS        270.616 vs   270.616 (Δ  +0.000)    313.419 vs   313.419 (Δ  +0.000)
  C2H6_elevated_BB3400       PASS        270.616 vs   270.616 (Δ  +0.000)    325.279 vs   325.279 (Δ  +0.000)

32/32 cases passed.
```

- 12 clear TS × {G2V, blackbody_3400K};
- clear Mars and CO₂-ice-cloud Mars;
- two haze loadings (thin τ≈0.5, thick τ≈8.6);
- 16 minor-gas cases: CH₄/CO/NH₃/C₂H₆ × {realistic, elevated} × 2 stars, single gas in N₂.

All `pver=300`.

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
=== TS300K_G2V: library vs baseline (rtol=1e-05, atol=0.0001) ===
  LWUP            PASS  max|d|=1.522e-05  max rel=5.755e-08
  LWDN            PASS  max|d|=1.347e-05  max rel=5.717e-08
  SWUP            PASS  max|d|=3.795e-06  max rel=5.860e-08
  SWDN            PASS  max|d|=1.526e-05  max rel=5.563e-08
  LWHR            PASS  max|d|=1.049e-07  max rel=5.792e-08
  SWHR            PASS  max|d|=5.337e-08  max rel=5.738e-08
  FSDTOA          PASS  max|d|=1.025e-05  max rel=3.015e-08
  LWUP_SPECTRAL   PASS  max|d|=1.802e-06  max rel=5.880e-08
  LWDN_SPECTRAL   PASS  max|d|=1.742e-06  max rel=5.901e-08
  SWUP_SPECTRAL   PASS  max|d|=4.763e-07  max rel=5.935e-08
  SWDN_SPECTRAL   PASS  max|d|=1.906e-06  max rel=5.950e-08

=== determinism (exact, max |delta| must be 0) ===
  repeat run_column    PASS  max|d|=0.000e+00
  run_columns[0]       PASS  max|d|=0.000e+00
  run_columns[1]       PASS  max|d|=0.000e+00

LWUP_TOM=264.242 W/m2   SWDN_SFC=232.653 W/m2
PASS
          46   18.769521675747711        18.788832606351178
          47   17.656206095697179        17.674371597009614
          48   10.835751511854070        10.846899821816184
          49   27.255744879613989        27.283786819467796
          50   29.660444200369945        29.690960203362255
          51   23.176187468943358        23.200032176101985
          52   33.266315393725478        33.300541279664351
          53   18.328235589593447        18.347092505164820
          54   41.043594278329799        41.085821779622485
          55   18.311832391384041        18.330672430605141
          56   33.064845972741374        33.098864577851543
          57   16.964925732088403        16.982380012978719
          58   34.927036215601994        34.962970726038442
          59   15.814528433667148        15.830799133921790
          60   8.9092321217831394        8.9183983416903967
          61   88.320571722913712        88.411439911291879
          62   47.584766292845771        47.633723646982865
          63   25.686269708800243        25.712696901798427
          64   5.7181658630539847        5.7240489700435058
          65   8.9164743245274565        8.9256479955395971
          66   5.6020569689378723        5.6078206178593906
          67   1.1259871698765906        1.1271456362708610
          68   2.2287270978257614        2.2310201128031113
          69  0.62870240085846552       0.62934923824958022
          70  0.13756915899441385        0.0000000000000000
          71   0.0000000000000000        0.0000000000000000
          72   0.0000000000000000        0.0000000000000000
          73   0.0000000000000000        0.0000000000000000
          74   0.0000000000000000        0.0000000000000000
          75   0.0000000000000000        0.0000000000000000
          76   0.0000000000000000        0.0000000000000000
          77   0.0000000000000000        0.0000000000000000
          78   0.0000000000000000        0.0000000000000000
          79   0.0000000000000000        0.0000000000000000
          80   0.0000000000000000        0.0000000000000000
          81   0.0000000000000000        0.0000000000000000
          82   0.0000000000000000        0.0000000000000000
          83   0.0000000000000000        0.0000000000000000
          84   0.0000000000000000        0.0000000000000000
 before/after scaling   680.00000000000011        680.00000000000000
 before/after scaling   679.99999999999932        679.99999999999932
 TOTAL SOLAR FLUX:   680.00000000000000        679.99999999999932
 optimizing longwave radiation bands
 optimized to   0.99900001287460327        400.00000000000000       K
 ---------------------------------------
 SW intervals           20          69
 SW gauss pts          153         552
 LW intervals            1          36
 LW gauss pts            1         288
 ---------------------------------------
 INIT_PLANCK: CREATING PLANCK FUNCTION TABLE
```

Residuals vs the committed baselines are ~6e-8 relative — the float32 storage
precision of the baseline files, not a solver difference.

## 5. Library API — C harness

```
cd tests/lib && make run
```

```
=== library vs baseline RTprofile_out_TS300K_G2V.nc (rtol=1e-05, atol=0.0001) ===
  LWDN                     PASS  max|d|=1.347e-05
  LWUP                     PASS  max|d|=1.522e-05
  SWDN                     PASS  max|d|=1.526e-05
  SWUP                     PASS  max|d|=3.795e-06
  LWHR                     PASS  max|d|=1.049e-07
  SWHR                     PASS  max|d|=5.337e-08
  FSDTOA                   PASS  max|d|=1.025e-05
  LWDN_SPECTRAL            PASS  max|d|=1.742e-06
  LWUP_SPECTRAL            PASS  max|d|=1.802e-06
  SWDN_SPECTRAL            PASS  max|d|=1.906e-06
  SWUP_SPECTRAL            PASS  max|d|=4.763e-07

=== batch determinism (exact) ===
  run_columns[0]           PASS  max|d|=0.000e+00
  run_columns[1]           PASS  max|d|=0.000e+00

LWUP_TOM=264.242 W/m2  SWDN_SFC=232.653 W/m2
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

If you change code and the suite still shows 32/32 Δ=0, the change is
bit-for-bit neutral. If it moves numbers, that is either a bug or an intended
physics change — and an intended one gets its baselines regenerated in the same
commit, called out explicitly in the message. See `REFACTOR_PLAN.md`.
