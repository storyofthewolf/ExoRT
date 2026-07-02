# ExoRT

**A two-stream radiative transfer code for 3D climate models and offline 1-D calculations.**

📦 [GitHub Repository](https://github.com/storyofthewolf/ExoRT)  
✉️ Eric T. Wolf — eric.wolf@colorado.edu  

---

## Overview

ExoRT is a flexible two-stream radiative transfer code designed for use with 3D climate models. It includes builds for a 1-D offline version and for direct interfacing with CESM1.2.1-ExoCAM and now with CESM3-planets.

## Integration with 3D models

### NCAR CESM1.2.1 Integration
- see the **ExoCAM Github Repository** for details and instructions: [ExoCAM](https://github.com/storyofthewolf/ExoCAM) 
- **Recommended version for CESM1.2.1-ExoCAM (as of September 2020):** `3dmodels/src.cam.n68equiv`
- Files in `3dmodels/src.cam.n68equiv/` are kept **identical** to their counterparts in `source/src.n68equiv/` and `source/src.main/`. The only intentional exception is `sys_rootdir.F90`, which encodes a machine-specific data path. When updating gas species in `source/`, always sync the corresponding files to `3dmodels/` in the same commit.

### NCAR CESM3-Planets Integration
- ExoRT n68equiv is now natively integrated into the **CESM3-planets framework** as an external component.
- **NEW for CESM3-planets (as of April 2026):** `3dmodels/src.cam7.n68equiv` 
- **GitHub Repository**: [NCAR/CESM3-planets](https://github.com/NCAR/CESM3-planets)
- **Project Wiki**: [CESM3-planets Wiki](https://github.com/NCAR/CESM3-planets/wiki)

---

## Directory Structure

```
ExoRT/
├── source/
│   ├── src.main/        # Drivers for offline calculation and shared radiation routines
│   ├── src.misc/        # Miscellaneous files and stubs from CESM origin (needed for offline runs)
│   ├── src.exort/       # v2 single RT bundle (84-band, HITRAN-2016 default, NH3/CO, CO2 clouds)
│   ├── src.n68equiv/    # legacy HITRAN-2016 reference (slated for retirement)
│   └── src.n84equiv/    # legacy HITRAN-2016 reference, +UV bins (slated for retirement)
├── data/
│   ├── cia/             # Collision-induced absorption data
│   ├── cloud/           # Cloud optical properties (Mie)
│   ├── continuum/       # CO2, H2O continuum coefficients (from MTCKD / LBLRTM)
│   ├── kdist/           # Correlated k-distributions (flat per-gas dirs: kdist/<gas>/)
│   └── stellar/         # Stellar spectra (renamed from data/solar/ in v2)
├── iofiles/             # Input/output files for the 1-D model
├── build/               # Build directory for the 1-D model
├── run/                 # Run directory for the 1-D model
├── 3dmodels/            # File sets to be linked with CESM
└── tools/               # Pre/post-processing scripts (Python and IDL)
```

---

## Radiative Transfer Versions

> **v2.0.0 (refactor branch) — single-bundle direction.** v2 collapses the RT
> versions into one bundle, **`src.exort`** (`make exort`), built on the 84-band
> grid with NH₃/CO. It runs on the validated **HITRAN-2016** native-gas line list
> by default (the HITRAN-2024 upgrade is a decoupled, still-unvalidated effort).
> CO₂ ice clouds + optional surface emissivity are folded in (Stage C; enable via
> `do_exo_clouds` in `user_nl_exort`); CARMA haze is folded in on the 1-D side
> (Stage C3; enable via `do_exo_haze` + a `carmammr` input array — the 84-band
> haze optics are provisional above 42,087 cm⁻¹ pending offline regen). The legacy
> H₂O-only / Archean bundles (`n28archean`, `n42h2o`, `n68h2o`) were **removed in
> v2** and live on in the `v1.0.0` tag. `n68equiv` and `n84equiv` remain in the
> tree for now as comparison references but are slated for retirement once
> `src.exort` is fully validated. See `REFACTOR_PLAN.md` and `REFACTOR_LOG.md`.
>
> ⚠️ **v2 status (2026-06-17): HITRAN-2024 k-coefficients are NOT yet validated.**
> The 84-band HITRAN-2024 tables in `data/kdist/` produce non-physical results
> for **H₂O, CO₂, and C₂H₆** (e.g. a 2-bar CO₂ column loses ~48% of its OLR;
> C₂H₆ is ~4× too weak; H₂O ~12% too strong). CH₄/NH₃/CO are clean, and the
> ExoRT code itself is verified correct (it reproduces the HITRAN-2016 results
> bit-for-bit when fed the h16 files). The defect is in the **HELIOS-K
> k-coefficient generation** for those gases and is under investigation. Until
> the tables are re-fit, treat `src.exort` LW results involving H₂O/CO₂/C₂H₆ as
> provisional. Use `tests/regression/gas_sweep.py` to re-check after a re-fit.

### `src.exort` ⭐ (v2, in validation)
- Single v2 bundle: 84-band grid (supersedes n68equiv + n84equiv; the runtime
  band optimizer condenses to the working set, so 84 bins cost no more than 68)
- Species: H₂O, CO₂, CH₄, C₂H₆, O₃, O₂, NH₃, CO (`nspecies = 8`)
- Correlated-k from HELIOS-K (Grimm et al. 2015); native gases on HITRAN-2024
  (O₂/O₃ HITRAN-2020), 8 Gauss points
- See the v2 status warning above re: HITRAN-2024 validation.

### `src.n68equiv` (legacy reference, HITRAN-2016)
> Was the recommended terrestrial version September 2020 – v1; kept in v2 only as
> a HITRAN-2016 comparison reference.

- Correlated-k coefficients from HELIOS-K (Grimm et al. 2015)
- **H₂O:** HITRAN 2016, Voigt lineshape, 25 cm⁻¹ cutoff, plinth removed; self/foreign continuum from MT_CKDv3.3 fit to Gauss points
- **CO₂:** HITRAN 2016, Perrin & Hartmann (1989) sub-Lorentzian lineshape, 500 cm⁻¹ cutoff, with CO₂–CO₂ CIA
- **CH₄:** HITRAN 2016, Voigt lineshape, 25 cm⁻¹ cutoff
- **C₂H₆:** HITRAN 2016, Voigt lineshape, 25 cm⁻¹ cutoff
- **NH₃:** HITRAN 2024, Voigt lineshape, 25 cm⁻¹ cutoff (added 2026-04-27; 3-D interface 2026-04-29)
- **CO:** HITRAN 2024, Voigt lineshape, 25 cm⁻¹ cutoff (added 2026-04-27; 3-D interface 2026-04-29)
- **CIA:** N₂–N₂, N₂–H₂, H₂–H₂ from HITRAN; CO₂–H₂ and CO₂–CH₄ from Turbet et al. (2020)
- 68 spectral intervals, 8 Gauss points
- Gas overlap via equivalent extinction absorption method (Amundsen et al. 2016): major gas treated with full 8-point correlated-k; minor species added as grey absorbers, selected on the fly
- Pressure range: 10 bar – 0.01 mb
- Temperature range: 100 K – 500 K
- Reference: [Wolf et al., PSJ 3:7 (2022)](https://doi.org/10.3847/PSJ/ac3f3d)

---

### `src.n84equiv`
- Same as `n68equiv`, with additional bins shortward of 0.24 µm
- Use for F stars (6500 K < T < ~10,000 K)

---

### `src.n28archean` — v1.0.0 only (removed from v2)
> Designed for Archean climate simulations. **Removed in v2; available in the
> `v1.0.0` tag.**

- Species: H₂O, CO₂, CH₄, N₂, H₂; HITRAN 2004, 28 spectral bins
- H₂O and CO₂ continuum from MT_CKD2.5 (see Halevy et al. 2009)
- Mixed gas k-distributions via LBLRTM (Mlawer et al. 1997; Shi et al. 2009)
- CO₂ up to multi-bar; reasonable agreement at 2-bar dry CO₂
- CH₄ up to 0.01 bar ⚠️ *Older line list misses CH₄ near-IR absorption*
- H₂O ⚠️ *Overestimates near-IR absorption around M-dwarfs due to coarse bands*
- Up to 100 bar total pressure; N₂–N₂, N₂–H₂, H₂–H₂ CIA
- Reference: [Wolf & Toon, Astrobiology 13(7), 1–18 (2013)](https://doi.org/10.1089/ast.2012.0936)

---

### `src.n42h2o` — v1.0.0 only (removed from v2)
> **Removed in v2; available in the `v1.0.0` tag.**

- Species: H₂O, N₂, H₂; 42 spectral bins; HITRAN 2012
- H₂O single-gas k-distributions via HELIOS-K (Grimm et al. 2015)
- Up to 10 bar total pressure; N₂–N₂, N₂–H₂, H₂–H₂ CIA
- Reference: [Kopparapu et al., ApJ 845:5 (2017)](https://doi.org/10.3847/1538-4357/aa7cf8)

---

### `src.n68h2o` — v1.0.0 only (removed from v2)
> **Removed in v2; available in the `v1.0.0` tag.**

- Same as `n42h2o`, extended to 68 spectral bins

---

## Building the Model

```bash
cd ../ExoRT/build

# v2 single bundle (recommended)
make exort

# legacy comparison references (HITRAN-2016)
make n68equiv
make n84equiv
```

> **macOS users:** the default compiler is `ifort`, which Intel discontinued and
> never ported to Apple Silicon (arm64). On any modern (M-series) Mac you **must**
> build with gfortran by adding `USER_FC=gfortran` to every `make` command, e.g.
> `make USER_FC=gfortran n68equiv`. See [macOS (Apple Silicon)](#macos-apple-silicon) below.

The executable is copied to `../ExoRT/run/`.

### Linux / Discover (default)

The default compiler is `ifort`. Ensure `nf-config` is on your PATH
(loaded automatically via the NetCDF module on Discover):

```bash
module load netcdf
make n68equiv
```

### macOS (Apple Silicon)

ifort is not available on Apple Silicon — Intel's compiler was never ported
to arm64 and was discontinued in 2023. Use gfortran instead:

```bash
make USER_FC=gfortran n68equiv
```

**Dependencies:** NetCDF-Fortran via Anaconda is recommended:

```bash
conda install -c conda-forge netcdf-fortran
```

Verify `nf-config` is on your PATH before building:

```bash
which nf-config
```

**Runtime environment:** The NetCDF libraries must be findable at runtime.

csh/tcsh:
```csh
setenv DYLD_LIBRARY_PATH /opt/anaconda3/lib
```

bash:
```bash
export DYLD_LIBRARY_PATH=/opt/anaconda3/lib
```

Add the appropriate line to your `~/.cshrc` or `~/.bashrc`.

**Note on compiler warnings:** The gfortran build produces warnings from
CESM-heritage infrastructure files (`infnan.F90`, `wrap_nf.F90`,
`shr_sys_mod.F90`) that ifort accepts silently. These are suppressed with
`-w` in the Makefile and are benign for the current codebase, but are
candidates for cleanup during refactoring of `src.misc`.

---

## Running the Model

1. Place your input file `RTprofile_in.nc` in the `../ExoRT/run/` directory.  
   Template inputs are in `../ExoRT/iofiles/`.

2. Set the solar spectrum and number of vertical levels in:
   ```
   source/src.main/exoplanet_mod.F90
   ```

3. Run the executable:
   ```bash
   cd ../ExoRT/run
   ./n42h2o.exe
   ```

Output is written to `RTprofile_out.nc`.

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

Use `tools/makeColumn.py` to generate input files. Species with VMR set to 0.0 are automatically omitted from the output file.

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

---

## Data Directory Contents

| Directory | Contents |
|-----------|----------|
| `data/aerosol/` | Aerosol optical constants and CARMA production rates |
| `data/cia/` | Collision-induced absorption data |
| `data/continuum/` | CO₂ and H₂O continuum coefficients from MTCKD (LBLRTM) |
| `data/kdist/` | Correlated k-distributions |
| `data/cloud/` | Cloud optical properties (Mie scattering) |
| `data/solar/` | Stellar spectra |

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
