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
│   ├── src.n*/          # Radiative transfer versions (n## = number of spectral intervals + descriptor)
│   └── experimental/    # RT builds under development or for testing
├── data/
│   ├── cia/             # Collision-induced absorption data
│   ├── cloud/           # Cloud optical properties (Mie)
│   ├── continuum/       # CO2, H2O continuum coefficients (from MTCKD / LBLRTM)
│   ├── kdist/           # Correlated k-distributions
│   └── solar/           # Stellar spectra
├── iofiles/             # Input/output files for the 1-D model
├── build/               # Build directory for the 1-D model
├── run/                 # Run directory for the 1-D model
├── 3dmodels/            # File sets to be linked with CESM
└── tools/               # IDL tools for file manipulation
```

---

## Radiative Transfer Versions

### `src.n68equiv` ⭐ Recommended
> Recommended for all terrestrial planet cases as of September 2020.

- Correlated-k coefficients from HELIOS-K (Grimm et al. 2015)
- **H₂O:** HITRAN 2016, Voigt lineshape, 25 cm⁻¹ cutoff, plinth removed; self/foreign continuum from MT_CKDv3.3 fit to Gauss points
- **CO₂:** HITRAN 2016, Perrin & Hartmann (1989) sub-Lorentzian lineshape, 500 cm⁻¹ cutoff, with CO₂–CO₂ CIA
- **CH₄:** HITRAN 2016, Voigt lineshape, 25 cm⁻¹ cutoff
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

### `src.n28archean`
> Designed for Archean climate simulations.

- Species: H₂O, CO₂, CH₄, N₂, H₂; HITRAN 2004, 28 spectral bins
- H₂O and CO₂ continuum from MT_CKD2.5 (see Halevy et al. 2009)
- Mixed gas k-distributions via LBLRTM (Mlawer et al. 1997; Shi et al. 2009)
- CO₂ up to multi-bar; reasonable agreement at 2-bar dry CO₂
- CH₄ up to 0.01 bar ⚠️ *Older line list misses CH₄ near-IR absorption*
- H₂O ⚠️ *Overestimates near-IR absorption around M-dwarfs due to coarse bands*
- Up to 100 bar total pressure; N₂–N₂, N₂–H₂, H₂–H₂ CIA
- Reference: [Wolf & Toon, Astrobiology 13(7), 1–18 (2013)](https://doi.org/10.1089/ast.2012.0936)

---

### `src.n42h2o`
- Species: H₂O, N₂, H₂; 42 spectral bins; HITRAN 2012
- H₂O single-gas k-distributions via HELIOS-K (Grimm et al. 2015)
- Up to 10 bar total pressure; N₂–N₂, N₂–H₂, H₂–H₂ CIA
- Reference: [Kopparapu et al., ApJ 845:5 (2017)](https://doi.org/10.3847/1538-4357/aa7cf8)

---

### `src.n68h2o`
- Same as `n42h2o`, extended to 68 spectral bins

---

## Building the Model

```bash
cd ../ExoRT/build

# Build specific version(s)
make n42h2o
make n68equiv
```

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

| Variable | Dimension | Description |
|----------|-----------|-------------|
| `TS` | (1) | Surface temperature |
| `PS` | (1) | Surface pressure |
| `TMID` | (pver) | Temperature at layer midpoints |
| `TINT` | (pverp) | Temperature at layer interfaces |
| `PDEL` | (pver) | Pressure thickness of each layer |
| `PINT` | (pverp) | Pressure at interface levels |
| `ZINT` | (pverp) | Height at interfaces |
| `H2OMMR` | (pver) | Water vapor mass mixing ratio (specific humidity, wet) |
| `CO2MMR` | (pver) | CO₂ mass mixing ratio (dry) |
| `CH4MMR` | (pver) | CH₄ mass mixing ratio (dry) |
| `O2MMR` | (pver) | O₂ mass mixing ratio (dry) |
| `O3MMR` | (pver) | O₃ mass mixing ratio (dry) |
| `H2MMR` | (pver) | H₂ mass mixing ratio (dry) |
| `N2MMR` | (pver) | N₂ mass mixing ratio (dry) |
| `asdir` | (1) | SW albedo (0.2–0.7 µm), direct |
| `asdif` | (1) | SW albedo (0.2–0.7 µm), diffuse |
| `aldir` | (1) | Near-IR albedo (0.7–4.0 µm), direct |
| `aldif` | (1) | Near-IR albedo (0.7–4.0 µm), diffuse |
| `coszrs` | (1) | Cosine of solar zenith angle |
| `mw` | (1) | Molecular weight of dry air |
| `cp` | (1) | Specific heat of dry air |

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
