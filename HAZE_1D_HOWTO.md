# Running 1-D Haze Calculations in ExoRT v2

A walkthrough for setting up, running, and checking a 1-D CARMA haze
calculation on the `refactor` branch.

This branch is a **work in progress**. Things are labelled where they are
rough, and there is a list of known warts and open contribution ideas at the
end. If something here doesn't match what you see, trust the code and tell us —
the docs drift faster than the source.

---

## 0. What you're working with

ExoRT v2 collapses the old five RT bundles into a single `src.exort` target on
an 84-band spectral grid. Haze enters as a **CARMA binwise aerosol**: you supply
a mass mixing ratio per size bin per level in the input file, and ExoRT looks up
pre-tabulated optical properties on exactly that bin grid.

The 1-D path has **no CARMA coupling** — you are not running a microphysics
model. You hand ExoRT a haze distribution and it does radiative transfer on it.
That's the whole contract.

**Gas absorption is HITRAN-2024** for H₂O/CO₂/CH₄/C₂H₆/NH₃/CO and
HITRAN-2020 for O₂/O₃. This is the default and the regression-verified line
list; no action needed. The HITRAN-2016 tables (H₂O corrected for its old
temperature-index error) are still in `data/kdist/`, reachable via
`run_regression.py --exort h16`. That mode swaps filenames temporarily and
restores them afterwards. `CHANGES.md` lists what moved between them.

---

## 1. One-time setup

### 1.1 Point the code at your checkout

`source/src.main/sys_rootdir.F90` hardcodes an absolute path used to find
`data/`. **You must change this** or every data file read will fail:

```fortran
character(len=256) :: exort_rootdir = '/your/path/to/ExoRT-refactor/'
```

Keep the trailing slash. (Yes, this should be an environment variable. See
"Contribution ideas".)

### 1.2 Build

Requires a NetCDF Fortran library with `nf-config` on your `PATH`.

```bash
cd build
make exort              # the v2 target; produces run/exort.exe
```

On Apple Silicon use `USER_FC=gfortran make exort` — the default is `ifort`,
which has no arm64 port.

If the link fails with `ld: library not found for -lnetcdf`, your
`nf-config --flibs` names the C library but doesn't give its directory
(common when netcdf-c and netcdf-fortran are separate packages). Work around
it without editing the Makefile:

```bash
env LIBRARY_PATH=$(nc-config --libdir) USER_FC=gfortran make exort
```

At runtime you may likewise need `DYLD_LIBRARY_PATH` (macOS) or
`LD_LIBRARY_PATH` (Linux) pointing at both NetCDF lib directories.

### 1.3 Check the vertical level count

`exo_pver` in `source/exoplanet_mod.F90` is **compile-time** and must match
your input file's level count exactly. It currently ships as:

```fortran
integer, parameter :: exo_pver = 300
```

Everything below assumes 300 levels. A different profile length means editing
this and rebuilding.

---

## 2. Build an input column with haze

`RTprofile_in.nc` needs an optional variable `carmammr` with shape
`(nlev, nelem, nbin)` = `(300, 1, 40)` in Fortran order. Two ways to make one.

### Route A — clone an existing fixture (recommended)

Best when you want a specific surface temperature. The committed fixtures cover
TS250K/273K/300K/320K/340K/360K, all at 300 levels.

```python
import shutil
import numpy as np
import netCDF4 as nc

REPO = '/your/path/to/ExoRT-refactor'
SRC  = f'{REPO}/tests/regression/fixtures/RTprofile_in_TS300K.nc'
OUT  = 'RTprofile_in.nc'

NELEM, NBIN = 1, 40

shutil.copyfile(SRC, OUT)
with nc.Dataset(OUT, 'a') as d:
    nlev = len(d.dimensions['pver'])
    arr = np.zeros((nlev, NELEM, NBIN))

    # Haze slab: layers 40..100 (0-based, top = 0), bins 15 and 18
    arr[40:100, 0, 15] = 1.16e-06     # kg/kg
    arr[40:100, 0, 18] = 5.80e-07     # kg/kg

    d.createDimension('nelements', NELEM)
    d.createDimension('nbins', NBIN)
    # Fortran reads carmammr(pver,nelem,nbin) => NetCDF dims fastest-last
    v = d.createVariable('carmammr', 'f8', ('nbins', 'nelements', 'pver'))
    v.title = 'CARMA haze binwise mass mixing ratio'
    v.units = 'kg/kg'
    v[:] = arr.T
```

Note the `.T` and the reversed dimension order — this trips everyone once.

### Route B — `makeColumn.py` from scratch

Full control over gas abundances, albedos, gravity etc.

```python
import sys
sys.path.insert(0, '/your/path/to/ExoRT-refactor/tools')
from makeColumn import make_column, add_haze_layer

NLEV = 300
carmammr = add_haze_layer(None, NLEV, 40, 100,
                          {15: 1.16e-06, 18: 5.80e-07})

make_column(output_file='RTprofile_in.nc',
            profile_tag='TS273K',
            co2vmr=400e-6,
            carmammr=carmammr)
```

`add_haze_layer(arr, nlev, klo, khi, {bin: mmr})` sets a slab over layers
`[klo, khi)` and can be called repeatedly to stack slabs at different
altitudes or sizes. Pass `None` the first time to allocate.

> **Gotcha:** `profile_tag` accepts only `TS273K` (300 levels),
> `smart_2bar_t250` (69), and `US1976` (49). The other TS profiles exist as
> *fixture files* but were never transcribed into `tools/profile_data.py` —
> an incomplete piece of the IDL→Python port. For TS300K etc., use Route A.

### Choosing bins and mixing ratios

40 bins, mass ratio 2.5, starting at rmin = 1 nm. Radii:

| bin | radius | bin | radius | bin | radius |
|----:|-------:|----:|-------:|----:|-------:|
| 10 | 0.021 µm | 16 | 0.133 µm | 22 | 0.828 µm |
| 12 | 0.039 µm | 18 | 0.244 µm | 24 | 1.526 µm |
| 14 | 0.072 µm | 20 | 0.450 µm | 26 | 2.811 µm |
| 15 | 0.098 µm | 21 | 0.610 µm | 28 | 5.177 µm |

Bins are **0-based** in the Python tools (matching the NetCDF array), so bin 15
is the 16th bin. Titan-like tholin haze lives around bins 12–20.

To hit a **target optical depth** rather than guessing mixing ratios, copy
`tests/regression/fixtures/make_haze_fixture.py`. It solves for the mixing
ratio given a target τ against the live optics table:

```
tau = sum over bins of ( mmr * masspath * Kext[band, bin] )
```

where `masspath = Δp/g` [kg m⁻²] summed over the slab. That script is the
cleanest worked example in the repo.

---

## 3. Configure the run

Put `RTprofile_in.nc` in `run/`, then create `run/user_nl_exort`:

```fortran
&exort_config
  solar_file     = 'G2V_SUN_n84.nc',
  shr_const_scon = 680.0,
  exo_g          = 9.80616,
  do_exo_haze    = .true.
/
```

| Variable | Meaning | Note |
|---|---|---|
| `solar_file` | stellar spectrum in `data/stellar/` | **must be an `_n84` file** |
| `shr_const_scon` | stellar constant ÷ 2 [W m⁻²] | 680.0 ≈ present Earth |
| `exo_g` | surface gravity [m s⁻²] | |
| `do_exo_haze` | enable the haze path | **required** |
| `do_exo_clouds` | enable H₂O/CO₂ clouds | independent of haze |

Omitted variables keep their compiled defaults. Active values are echoed at
startup under `=== exort_config ===`.

> **Haze is silently inert unless BOTH** `do_exo_haze = .true.` **and**
> `carmammr` is present in the input file. This is the single most common way
> to get a confusing null result.
>
`iofiles/user_nl_exort.template` is the canonical starting point and documents
every variable with its default:

```bash
cp iofiles/user_nl_exort.template run/user_nl_exort
# then edit run/user_nl_exort — set do_exo_haze = .true.
```

`run/` is a scratch directory: `user_nl_exort`, `*.nc`, `*.exe` and
`libexort.*` there are all gitignored, so your run config never shows up as a
repo modification.

---

## 4. Run

```bash
cd run
./exort.exe
```

Output goes to `run/RTprofile_out.nc`. Confirm the haze was picked up — startup
prints an inventory:

```
--- aerosol variables in input file ---
   carmammr:  found
```

`not found, set to zero` means your variable name, shape, or file is wrong.

Runtime is well under a second for a single 300-level column.

---

## 5. Check the answer

### Sanity: run it twice

The cheapest meaningful test is haze on vs. off with an otherwise identical
deck — just flip `do_exo_haze`:

```
do_exo_haze = .true.    SW DN SURF  194.29    LW UP TOA  263.00
do_exo_haze = .false.   SW DN SURF  232.65    LW UP TOA  264.24
```

Tholin haze should **darken the surface** (SW down falls) and **cool via
anti-greenhouse** (OLR falls slightly). If SW at the surface didn't move, the
haze isn't active.

### Reference points

Committed regression cases, TS300K / G2V / `scon=680` / `g=9.80616`:

| case | visible τ | OLR [W m⁻²] | SWDN_SFC [W m⁻²] |
|---|---:|---:|---:|
| clear (`TS300K_G2V`) | 0 | 264.242 | 232.653 |
| `TS300K_haze_G2V` | 0.5 | 263.000 | 194.278 |
| `TS300K_hazethick_G2V` | 8.6 | 244.676 | 63.514 |

The Route A snippet in §2 reproduces the middle row exactly.

### Plot spectra

```bash
cd tools
python plotspectra_1D.py --f1 ../run/RTprofile_out.nc --save
python plotspectra_1D.py --f1 ../run/RTprofile_out.nc \
                         --f2 ../run/RTprofile_out_clear.nc \
                         --label1 "haze" --label2 "clear" --save
```

### Full regression suite

Before and after any code change:

```bash
cd tests/regression
python run_regression.py            # expect 32/32, Δ = 0.000
```

---

## 6. Haze optics — what's solid and what isn't

ExoRT reads `data/aerosol/haze_n84_b40_mie.nc`, set in
`source/src.exort/cloud.F90`.

**Mie spheres, regenerated 2026-09-15** from the Khare et al. (1984) Titan
tholin optical constants by `tools/makeCARMAOptics.py`, across all 84 bands.
Bands 69–77 rest on real Khare measurements. Bands 78–84 are extrapolated —
Khare stops at 0.0588 µm — but that's EUV, negligible for a G or K star. For a
hot or active star with real EUV output, those bands are **not** constrained by
measurement.

**A fractal-aggregate table also exists** (`haze_n84_b40_fractal_interp.nc`)
and is physically more appropriate for Titan-like aggregate hazes. It is
**still provisional above band 68** — its UV bands are a placeholder copy of
band 68 — because regenerating it needs an external mean-field fractal solver
(`fractaloptics.exe`, on Discover) that isn't in this repo. Switching is a
one-line change in `cloud.F90` plus a rebaseline. Mie is the default because
it's the one that's honestly computed end-to-end.

### Mie vs. fractal: what the comparison actually is

Both tables are indexed by the **same CARMA bins**, and a CARMA bin is a bin of
constant *mass*. `rbins` is the equivalent-sphere radius — the radius a sphere
of that bin's mass would have — and it is byte-identical between the two files.
So comparing the two tables at the same bin index is already the physically
meaningful comparison: **a Mie sphere against an equal-mass fractal
aggregate**, not against a sphere of some other size.

The aggregate is built from 50 nm monomers with `nmon = (r_eq/r_mon)³`, a
size-dependent fractal dimension `df = 2.4 − 0.9·exp(−nmon/500)`, and a
resulting aggregate radius `R_f` much larger than `r_eq` at fixed mass:

| bin | r_eq | nmon | df | R_f | R_f/r_eq |
|----:|-----:|-----:|---:|----:|---------:|
| 12 | 0.039 µm | 1 | 3.00 | 0.039 µm | 1.00 |
| 15 | 0.098 µm | 7 | 1.51 | 0.189 µm | 1.93 |
| 18 | 0.244 µm | 116 | 1.69 | 0.839 µm | 3.44 |
| 21 | 0.610 µm | 1819 | 2.38 | 1.177 µm | 1.93 |
| 24 | 1.526 µm | 28422 | 2.40 | 3.586 µm | 2.35 |

Bins at or below the monomer size (bin ≤ 12) are single monomers, so the
fractal treatment reduces exactly to Mie — the two tables agree there to 0.4%,
which is a good consistency check on both.

**There is no single ratio between the tables.** `Kext_mie / Kext_fractal` at
equal mass depends strongly on both wavelength and bin:

| λ | bin 12 | bin 15 | bin 18 | bin 21 | bin 24 |
|---|-------:|-------:|-------:|-------:|-------:|
| 32.8 µm | 1.00 | 1.00 | 1.00 | 1.01 | 1.08 |
| 8.08 µm | 1.00 | 1.00 | 1.03 | 1.22 | 3.06 |
| 1.77 µm | 1.00 | 1.28 | 5.23 | 3.87 | 0.69 |
| 0.645 µm | 1.00 | 1.80 | 4.05 | 0.34 | 0.16 |
| 0.217 µm † | 1.48 | 0.75 | 0.25 | 0.14 | 0.06 |
| 0.073 µm † | 3.07 | 0.55 | 0.22 | 0.13 | 0.06 |

† These two rows are above band 68, where the *fractal* table is still the
provisional band-68 copy. The ratio there is therefore partly an artifact of
the placeholder, not a clean Mie-vs-aggregate result. Rows at 0.645 µm and
longer are both-real.

The structure is physical. In the long-wavelength/Rayleigh limit (x ≪ 1)
extinction depends on mass, not shape, so the ratio → 1. Where the aggregate's
*open structure* lets its monomers act more independently than the compact
sphere's interior — the sphere being partly self-shadowed — the aggregate wins
and the ratio drops well below 1 (large bins in the UV/visible). In between,
the compact sphere can be the more efficient scatterer at its own size
parameter and the ratio exceeds 1.

**Practical consequence:** haze optical depth is not portable between the two
tables by any scalar factor. If you switch tables, recompute τ from the table
you are actually using — and quote which one with any τ you report.

### Units, because this bit was wrong until recently

`Kext` in **all** haze files is **m² kg⁻¹**, and the loaders apply **no**
conversion. `rbins` is in **cm**.

Older files carried `units = "cm2 g-1"` and `"microns"` attributes that
contradicted their own values — a stale label from an old IDL writer. The
attributes were corrected in place on 2026-09-15 (values verified
bit-identical), and a matching `kcarma*0.1` was removed from
`initialize_hazeopts`.

**If you ran 1-D haze on this branch before 2026-09-15, your haze optical
depth was 10× too small.** Rerun. Clear-sky and cloud results were never
affected. `main` never had the bug.

To check units on any new aerosol file, compare against the geometric-optics
limit `Kext = 3Q/(4ρr)` for the largest bin, where Q → 2. A 10× discrepancy
means someone handed you cgs.

---

## 7. Known warts

Things that will look broken because they partly are.

- **`exort_rootdir` is a hardcoded absolute path** (§1.1). Must be edited per
  machine; not overridable at runtime.
- **`tools/profile_data.py` is a partial port** — only 3 of the TS profiles
  made it from `profile_data.pro`.
- **`tools/plotspectra_1D.py` has no haze-specific diagnostics** — it plots
  fluxes, not aerosol optical depth. There is no built-in way to dump the τ
  your deck actually produced; compute it yourself (see
  `make_haze_fixture.py`).
- **No CLI flag for haze in `makeColumn.py`** — Python API only.
- **`carmammr` shape is not validated on read.** `opt_carma` in `io_1D.F90`
  reads straight into a fixed `(pver,1,40)` buffer. A wrong bin count won't
  produce a clean error. The *optics file* dimensions are checked and hard-stop
  on mismatch; the input file's are not.
- **40 bins / 1 element are compile-time constants** in
  `source/src.exort/radgrid.F90` (`nbin_carma`, `nelem_carma`). They match the
  ExoCAM-CARMA grid, and the optics are tabulated on exactly that grid with no
  radius interpolation — so changing them means editing source, rebuilding,
  regenerating optics *and* every input deck, and breaking 3-D compatibility.
  Not a casual knob.
- **The 3-D bundle (`3dmodels/src.cam.exort`) has never had a real ExoCAM
  build or run.** It compiles against stubs only.
- **IDL tools are being retired** in favour of the Python ones in `tools/`, but
  the `.pro` files remain the reference implementation and are still the only
  route for some things (fractal optics, most data-prep).

---

## 8. Contribution ideas

If you're inclined, these are self-contained and genuinely useful:

1. **Backfill `profile_data.py`** with the missing TS profiles so
   `make_column(profile_tag='TS300K')` works. Mechanical; arrays can be read
   out of the committed fixtures.
2. **Validate `carmammr` dimensions on read** in `io_1D.F90` — mirror the
   dimension checks `initialize_hazeopts` already does, and fail with a clear
   message.
3. **Make `exort_rootdir` runtime-settable** (environment variable or
   namelist), falling back to the compiled default.
4. **Write the haze optical depth to `RTprofile_out.nc`** as a diagnostic.
   Right now nothing in the output tells you what τ your deck produced, so it
   has to be recomputed by hand outside the model.
5. **Add a `--haze` CLI path to `makeColumn.py`** so a haze deck can be built
   without writing Python.
6. **Regenerate the fractal optics** on the 84-band grid — needs
   `fractaloptics.exe` and a machine that has it. This is the biggest
   outstanding item on the haze path.

Please run `python run_regression.py` (expect **32/32, Δ = 0**) before and
after any change. If a change is *intended* to move the physics, regenerate
baselines in the same commit and say so explicitly in the message — that's the
project convention.

---

## Quick reference

```bash
# build (macOS/arm64)
cd build && USER_FC=gfortran make exort

# run
cd run && ./exort.exe                 # reads RTprofile_in.nc + user_nl_exort
                                      # writes RTprofile_out.nc

# check
cd tests/regression && python run_regression.py       # expect 32/32 Δ=0

# plot
cd tools && python plotspectra_1D.py --f1 ../run/RTprofile_out.nc --save
```

| What | Where |
|---|---|
| Haze optics table | `data/aerosol/haze_n84_b40_mie.nc` |
| Which table is used | `source/src.exort/cloud.F90` |
| Bin/element counts | `source/src.exort/radgrid.F90` |
| Gas k-file selection | `source/src.exort/kabs.F90` |
| Level count `exo_pver` | `source/exoplanet_mod.F90` |
| Data root path | `source/src.main/sys_rootdir.F90` |
| Haze input read | `source/src.main/io_1D.F90` (`opt_carma`) |
| Optics loader | `source/src.exort/initialize_rad_mod_1D.F90` |
| Worked τ-targeting example | `tests/regression/fixtures/make_haze_fixture.py` |
| Optics generator | `tools/makeCARMAOptics.py`, `tools/mie.py` |
