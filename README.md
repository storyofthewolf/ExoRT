# ExoRT

**A two-stream radiative transfer code for 3D climate models and offline 1-D calculations.**

📦 [GitHub Repository](https://github.com/storyofthewolf/ExoRT)  
✉️ Eric T. Wolf — eric.wolf@colorado.edu  

---

## Overview

ExoRT is a flexible two-stream radiative transfer code designed for use with 3D
climate models. It builds as a 1-D offline executable, as a shared library with
C and Python bindings, and as source bundles for CESM1.2.1-ExoCAM and
CESM3-planets.

> **Branches.** `main` carries v1 (`v1.0.0`, the Wolf et al. 2022 release;
> `v1.1.0` adds the HITRAN-2016 H₂O temperature-index fix). This `refactor`
> branch is **v2**, which collapses the RT versions into one bundle,
> `src.exort`. `CHANGES.md` lists everything v2 changes and the measured
> effect of each physics change; `REFACTOR_PLAN.md` lists what is left before
> a v2 release.

## Integration with 3D models

### NCAR CESM1.2.1 (ExoCAM)
- See the **ExoCAM repository** for instructions: [ExoCAM](https://github.com/storyofthewolf/ExoCAM).
- **Production runs:** use `3dmodels/src.cam.n68equiv` from **`main`**
  (`v1.1.0`). The pre-v2 bundles on this branch are frozen and point at the v1
  data layout, which v2 no longer has.
- **v2:** `3dmodels/src.cam.exort` has CO₂ ice clouds, CARMA haze and surface
  emissivity, enabled by `-cppdefs`. It compiles against stubbed CESM in every
  configuration but **has not yet run in a real ExoCAM case**. Setup and
  requirements are in its `README`.

### NCAR CESM3-Planets
- ExoRT n68equiv is integrated into the **CESM3-planets framework** as an
  external component: `3dmodels/src.cam7.n68equiv` (April 2026).
- [NCAR/CESM3-planets](https://github.com/NCAR/CESM3-planets) ·
  [Wiki](https://github.com/NCAR/CESM3-planets/wiki)
- How to bring it onto the v2 source: `CESM.CAM7.md`.

---

## Directory Structure

```
ExoRT/
├── source/
│   ├── src.main/        # Drivers for offline calculation and shared radiation routines
│   ├── src.misc/        # Miscellaneous files and stubs from CESM origin (needed for offline runs)
│   ├── src.exort/       # v2 single RT bundle (84-band, HITRAN-2024, NH3/CO, clouds, haze)
│   ├── src.n68equiv/    # legacy HITRAN-2016 reference (slated for retirement)
│   └── src.n84equiv/    # legacy HITRAN-2016 reference, +UV bins (slated for retirement)
├── data/
│   ├── aerosol/         # CARMA haze optics
│   ├── cia/             # Collision-induced absorption data
│   ├── cloud/           # Cloud optical properties (Mie)
│   ├── continuum/       # CO2, H2O continuum coefficients (from MTCKD / LBLRTM)
│   ├── kdist/           # Correlated k-distributions (flat per-gas dirs: kdist/<gas>/)
│   └── stellar/         # Stellar spectra (renamed from data/solar/ in v2)
├── iofiles/             # Input/output templates for the 1-D model
├── build/               # Build directory
├── run/                 # Run directory for the 1-D model
├── 3dmodels/            # File sets to be linked with CESM
├── tests/               # Regression suite, library tests, CAM compile check
└── tools/               # Pre/post-processing scripts (Python and IDL)
```

---

## Radiative Transfer Versions

### `src.exort` ⭐ (v2)
- One bundle on the **84-band grid**: the 68-band grid plus UV bins shortward
  of 0.24 µm, so it covers F stars (6500 K < T < ~10,000 K) as well. The
  runtime band optimizer condenses to the working set, so 84 bins cost no more
  than 68.
- **Species:** H₂O, CO₂, CH₄, C₂H₆, O₃, O₂, NH₃, CO (`nspecies = 8`).
- **Correlated-k from HELIOS-K** (Grimm et al. 2015), 8 Gauss points.
  **HITRAN-2024** for H₂O/CO₂/CH₄/C₂H₆/NH₃/CO; HITRAN-2020 for O₂/O₃.
  HITRAN-2016 tables remain for comparison (`run_regression.py --exort h16`).
  - **H₂O:** Voigt, 25 cm⁻¹ cutoff, plinth removed; self/foreign continuum
    from MT_CKD 3.3 fit to the Gauss points.
  - **CO₂:** Perrin & Hartmann (1989) sub-Lorentzian line shape, 500 cm⁻¹
    cutoff, with CO₂–CO₂ CIA.
  - **CH₄, C₂H₆, NH₃, CO:** Voigt, 25 cm⁻¹ cutoff.
- **CIA:**
  - N₂–N₂, N₂–H₂, H₂–H₂ from HITRAN.
  - CO₂–H₂ and CO₂–CH₄ from Turbet et al. (2020).
  - O₂–O₂, O₂–N₂, O₂–CO₂.
- **Gas overlap:** equivalent extinction (Amundsen et al. 2016). The major gas
  gets the full 8-point correlated-k; minor species are added as grey
  absorbers, selected on the fly.
- **Clouds and aerosols** (runtime flags):
  - H₂O clouds (MCICA);
  - CO₂ ice clouds;
  - CARMA haze, with Mie optics from Khare et al. (1984) tholin indices over
    all 84 bands (`HAZE_1D_HOWTO.md`).
- **Surface:** optional thermal emissivity.
- **Table range:** pressure 0.01 mb – 10 bar; temperature 100 K – 500 K.
  k-files are checked against this grid at load.
- Reference for the n68equiv lineage: [Wolf et al., PSJ 3:7 (2022)](https://doi.org/10.3847/PSJ/ac3f3d).
  That paper's validation used the HITRAN-2016 H₂O table with the 25 K
  temperature-index error, since fixed; see `CHANGES.md`.

### Legacy and v1-only versions

| Version | Status | Notes / reference |
|---|---|---|
| `n68equiv` | v2 comparison reference only (HITRAN-2016) | The v1 recommended terrestrial version (Sept 2020 – v1); physics as above on 68 bands. [Wolf et al. 2022](https://doi.org/10.3847/PSJ/ac3f3d) |
| `n84equiv` | v2 comparison reference only (HITRAN-2016) | `n68equiv` + UV bins, for F stars |
| `n28archean` | `v1.0.0` tag only | Archean: H₂O, CO₂, CH₄, N₂, H₂; HITRAN 2004, 28 bins; LBLRTM mixed-gas k; misses CH₄ near-IR and overestimates H₂O near-IR around M dwarfs. [Wolf & Toon 2013](https://doi.org/10.1089/ast.2012.0936) |
| `n42h2o` | `v1.0.0` tag only | H₂O/N₂/H₂, HITRAN 2012, 42 bins, HELIOS-K. [Kopparapu et al. 2017](https://doi.org/10.3847/1538-4357/aa7cf8) |
| `n68h2o` | `v1.0.0` tag only | `n42h2o` extended to 68 bins |

---

## Building the Model

```bash
cd ExoRT/build

make exort       # v2 single bundle            -> run/exort.exe
make libexort    # shared library, C/Python API -> run/libexort.dylib|.so
make n68equiv    # legacy HITRAN-2016 references (comparison only)
make n84equiv
```

`libexort` exposes the same single-column physics as `exort.exe` to C and
Python callers (in-process, no NetCDF profile round-trips) — see
`tools/exort_pytools/README.md`. All targets build with OpenMP; `make
OMPFLAGS= <target>` builds without it.

> **macOS users:** the default compiler is `ifort`, which Intel discontinued and
> never ported to Apple Silicon (arm64). On any modern (M-series) Mac you **must**
> build with gfortran by adding `USER_FC=gfortran` to every `make` command, e.g.
> `make USER_FC=gfortran exort`. See [macOS (Apple Silicon)](#macos-apple-silicon) below.

### Linux / Discover (default)

The default compiler is `ifort`. Ensure `nf-config` is on your PATH
(loaded automatically via the NetCDF module on Discover):

```bash
module load netcdf
make exort
```

### macOS (Apple Silicon)

ifort is not available on Apple Silicon — Intel's compiler was never ported
to arm64 and was discontinued in 2023. Use gfortran instead:

```bash
make USER_FC=gfortran exort
```

If the link fails with `library not found for -lnetcdf`, `nf-config --flibs`
is naming the netCDF-C library without its directory. Add it:

```bash
env LIBRARY_PATH=$(nc-config --libdir) USER_FC=gfortran make exort
```

**Dependencies:** NetCDF-Fortran, e.g. from Homebrew (`brew install
netcdf-fortran`) or conda (`conda install -c conda-forge netcdf-fortran`).
Verify `nf-config` is on your PATH before building (`which nf-config`).

**Runtime environment:** the NetCDF libraries must be findable at runtime —
point `DYLD_LIBRARY_PATH` at `$(nf-config --prefix)/lib`:

csh/tcsh:
```csh
setenv DYLD_LIBRARY_PATH `nf-config --prefix`/lib
```

bash:
```bash
export DYLD_LIBRARY_PATH=$(nf-config --prefix)/lib
```

**Compiler warnings:** gfortran warns on some CESM-heritage infrastructure
files (`infnan.F90`, `wrap_nf.F90`) that ifort accepts silently. They are
suppressed with `-w` in the Makefile and are benign.

---

## Running the Model

1. Place your input file `RTprofile_in.nc` in `ExoRT/run/`. Build one with
   `tools/makeColumn.py`, or start from the fixtures in
   `tests/regression/fixtures/`. Its level count must match `exo_pver` in
   `source/exoplanet_mod.F90` (compile-time; 300 by default).

2. Optionally set the star, insolation, gravity and physics switches at
   runtime, with no rebuild:
   ```bash
   cp iofiles/user_nl_exort.template run/user_nl_exort   # then edit
   ```
   The stellar file must match the grid (`*_n84.nc` for `exort`).

3. Run:
   ```bash
   cd ExoRT/run
   ./exort.exe
   ```

Output is written to `RTprofile_out.nc`. An input may hold many columns (an
`ncol` dimension), which are solved in parallel — see
`MULTICOLUMN_BATCH_HOWTO.md`. For haze, see `HAZE_1D_HOWTO.md`.

---

## Input File: `RTprofile_in.nc`

> Note: `pverp = pver + 1` (interface levels = midpoint levels + 1)

Gas species variables are **optional**: if a variable is absent from the input file, ExoRT sets that species to zero and prints a clean diagnostic rather than an error. Only `TS`, `PS`, the pressure/temperature/height arrays, the albedos, `coszrs`, `mw`, and `cp` are required.

| Variable | Dimension | Required | Description |
|----------|-----------|----------|-------------|
| `ts` | (1) | yes | Surface temperature (K) |
| `ps` | (1) | yes | Surface pressure (Pa) |
| `tmid` | (pver) | yes | Temperature at layer midpoints (K) |
| `tint` | (pverp) | yes | Temperature at layer interfaces (K) |
| `pdel` | (pver) | yes | Pressure thickness of each layer (Pa) |
| `pint` | (pverp) | yes | Pressure at interface levels (Pa) |
| `zint` | (pverp) | yes | Height at interfaces (m) |
| `asdir` | (1) | yes | SW albedo, direct |
| `asdif` | (1) | yes | SW albedo, diffuse |
| `aldir` | (1) | yes | Near-IR albedo, direct |
| `aldif` | (1) | yes | Near-IR albedo, diffuse |
| `coszrs` | (1) | yes | Cosine of solar zenith angle |
| `mw` | (1) | yes | Molecular weight of dry air (g/mol) |
| `cp` | (1) | yes | Specific heat of dry air (J/kg/K) |
| `h2ommr` | (pver) | optional | H₂O specific humidity, kg(wv)/kg(air) |
| `co2mmr` | (pver) | optional | CO₂ mass mixing ratio (dry) |
| `ch4mmr` | (pver) | optional | CH₄ mass mixing ratio (dry) |
| `c2h6mmr` | (pver) | optional | C₂H₆ mass mixing ratio (dry) |
| `nh3mmr` | (pver) | optional | NH₃ mass mixing ratio (dry) |
| `commr` | (pver) | optional | CO mass mixing ratio (dry) |
| `o2mmr` | (pver) | optional | O₂ mass mixing ratio (dry) |
| `o3mmr` | (pver) | optional | O₃ mass mixing ratio (dry) |
| `h2mmr` | (pver) | optional | H₂ mass mixing ratio (dry) |
| `n2mmr` | (pver) | optional | N₂ mass mixing ratio (dry) |
| `grav` | (1) | optional | Per-column surface gravity (m s⁻²); absent or ≤0 → namelist `exo_g` |
| `scon` | (1) | optional | Per-column stellar constant ÷ 2 (W m⁻²); absent or ≤0 → namelist `shr_const_scon` |

Use `tools/makeColumn.py` to generate input files. Species with VMR set to 0.0 are automatically omitted from the output file.

**Multi-column input:** the file may carry an optional `ncol` dimension, in
which case every variable gains a column axis (`ts(ncol)`,
`tmid(ncol, pver)` in C/Python dimension order) and one invocation solves
all columns. Stack single-column files with
`python tools/stackColumns.py col1.nc col2.nc -o RTprofile_in.nc`. All
columns share the runtime namelist config and the stellar *spectrum*, but
each column may carry its own gravity and insolation via the optional
`grav`/`scon` variables (write them with `makeColumn.py --write-grav` /
`--scon`); mixed-star batches are done as an outer loop of runs instead.
Without `ncol` the file behaves exactly as before. Multi-column batches
solve in parallel with OpenMP — set `OMP_NUM_THREADS` to control the
thread count (results are bitwise independent of it). See
`MULTICOLUMN_BATCH_HOWTO.md` for a step-by-step recipe (executable and
Python-library paths).

---

## Output File: `RTprofile_out.nc`

> Contains all input variables plus the following radiative flux and heating rate fields.

| Variable | Dimension | Description |
|----------|-----------|-------------|
| `LWUP` | (pverp) | Longwave upwelling flux (W m⁻²) |
| `LWDN` | (pverp) | Longwave downwelling flux (W m⁻²) |
| `SWUP` | (pverp) | Shortwave upwelling flux (W m⁻²) |
| `SWDN` | (pverp) | Shortwave downwelling flux (W m⁻²) |
| `LWUP_SPECTRAL` | (ntot_wavlnrng, pverp) | Spectral LW upwelling flux (W m⁻² per interval) |
| `LWDN_SPECTRAL` | (ntot_wavlnrng, pverp) | Spectral LW downwelling flux (W m⁻² per interval) |
| `SWUP_SPECTRAL` | (ntot_wavlnrng, pverp) | Spectral SW upwelling flux (W m⁻² per interval) |
| `SWDN_SPECTRAL` | (ntot_wavlnrng, pverp) | Spectral SW downwelling flux (W m⁻² per interval) |
| `LWHR` | (pver) | Longwave heating rate (K/Earth day) |
| `SWHR` | (pver) | Shortwave heating rate (K/Earth day) |

For a multi-column input (`ncol` dimension), every output variable gains the
matching column axis; single-column runs keep exactly the layout above.

---

## Data Directory Contents

| Directory | Contents |
|-----------|----------|
| `data/aerosol/` | Aerosol optical constants and CARMA production rates |
| `data/cia/` | Collision-induced absorption data |
| `data/continuum/` | CO₂ and H₂O continuum coefficients from MTCKD (LBLRTM) |
| `data/kdist/` | Correlated k-distributions, one directory per gas (`kdist/<gas>/`) |
| `data/cloud/` | Cloud optical properties (Mie scattering) |
| `data/stellar/` | Stellar spectra, one file per grid (`*_n84.nc` for `exort`) |

---

## Tools

Python scripts in `tools/` handle the primary 1-D pre/post-processing workflow. Legacy IDL scripts (`.pro`) remain for reference.

### Building input profiles

**`tools/makeColumn.py`** — generates `RTprofile_in.nc` from a P-T profile and user-specified gas mixing ratios. Edit the `USER SETTINGS` block at the top or pass CLI flags:

```bash
cd tools/
python makeColumn.py --defaults                        # show current settings
python makeColumn.py --output RTprofile_in.nc          # use USER SETTINGS
python makeColumn.py --profile TS273K \
    --co2vmr 400e-6 --nh3vmr 1e-6 \
    --output RTprofile_in_nh3.nc
python makeColumn.py --zero-h2o --output RTprofile_in_dry.nc
```

Gas species with VMR = 0.0 are omitted from the file; ExoRT reads absent species as zero. Available P-T profiles (from `profile_data.py`):

| Tag | Description |
|-----|-------------|
| `TS273K` | Moist-adiabat Earth, Ts = 273 K, 300 levels |
| `US1976` | US Standard Atmosphere 1976, 49 levels |
| `smart_2bar_t250` | 2-bar dry CO₂ Mars atmosphere, Ts ≈ 250 K, 69 levels |

### Plotting spectral output

**`tools/plotspectra_1D.py`** — plots SW and LW spectra from `RTprofile_out.nc`. The spectral grid is auto-detected from the file. Supports side-by-side comparison of two runs.

```bash
cd tools/
python plotspectra_1D.py                               # interactive display
python plotspectra_1D.py --save                        # save spectra_sw.png / spectra_lw.png
python plotspectra_1D.py --f2 ../run/RTprofile_out_ref.nc --label2 "no NH3"
```

SW panel: wavelength axis (µm), W m⁻² µm⁻¹. LW panel: wavenumber axis (cm⁻¹), W m⁻² cm. Also prints a band-by-band flux table to stdout.

---

## References

**ExoRT applications**
- Wolf, E.T. & Toon, O.B. (2013). Hospitable Archean Climates Simulated by a General Circulation Model. *Astrobiology*, 13(7), 1–18. https://doi.org/10.1089/ast.2012.0936
- Kopparapu, R.K. et al. (2017). Habitable Moist Atmospheres on Terrestrial Planets near the Inner Edge of the Habitable Zone around M Dwarfs. *ApJ*, 845:5. https://doi.org/10.3847/1538-4357/aa7cf8

**Correlated-k / line-by-line**
- Grimm, S.L. & Heng, K. (2015). HELIOS-K: An Ultrafast, Open-source Opacity Calculator for Radiative Transfer. *ApJ*, 808:182. https://doi.org/10.1088/0004-637X/808/2/182
- Mlawer, E.J. et al. (1997). Radiative transfer for inhomogeneous atmospheres: RRTM, a validated correlated-k model for the longwave. *JGR Atmospheres*, 102(D14), 16663–16682. https://doi.org/10.1029/97JD00237
- Shi, G. et al. (2009). An extension of the correlated-k distribution method to non-homogeneous atmospheres. *JQSRT*, 110(18), 2268–2278. https://doi.org/10.1016/j.jqsrt.2009.05.005

**Gas overlap**
- Amundsen, D.S. et al. (2016). Accurate treatment of transmission, emission, and scattering in parameterised radiation schemes for a general circulation model. *A&A*, 564, A59. https://doi.org/10.1051/0004-6361/201323169

**Line shapes / CO₂ broadening**
- Perrin, M.Y. & Hartmann, J.M. (1989). Temperature-dependent measurements and modeling of absorption by CO₂–N₂ mixtures in the far line-wings of the 4.3 µm CO₂ band. *JQSRT*, 42(4), 311–317. https://doi.org/10.1016/0022-4073(89)90077-0

**Continuum absorption**
- Halevy, I., Pierrehumbert, R.T. & Schrag, D.P. (2009). Radiative implications of CO₂ cloud formation in a Neoproterozoic snowball Earth. *JGR Atmospheres*, 114, D18112. https://doi.org/10.1029/2009JD011915

**CIA data**
- Turbet, M. et al. (2020). CO₂ condensation is a serious limit to the deglaciation of Earth-like planets. *Earth and Planetary Science Letters*, 531, 115759. https://doi.org/10.1016/j.epsl.2019.115959
- Gordon, I.E. et al. (HITRAN collaboration) — CIA data retrieved from the [HITRAN database](https://hitran.org/).
