# Running a Multithreaded, Multi-Column Batch

Step-by-step recipe for solving many columns in one ExoRT invocation with
OpenMP threading (Stage E1/E2/E2b machinery). Two paths are covered: the
standalone executable (`exort.exe` + NetCDF files) and the library
(`libexort` + Python), which share the identical per-column solve
(`run_one_column`).

Results are **bitwise independent of the thread count** — a batch gives the
same answers at `OMP_NUM_THREADS=1` and `=8`, and the same answers as running
each column as its own single-column job (gated by
`tests/regression/multicol_check.py` and `percol_config_check.py`).

---

## What is shared vs. per-column

One batch = one process-level configuration + N independent column states.

| Scope | Quantity | Set by |
|-------|----------|--------|
| **Process (all columns)** | Stellar spectrum *shape* (`solar_file`) | namelist |
| | Default stellar constant (`shr_const_scon`) | namelist |
| | Default gravity (`exo_g`) | namelist |
| | Cloud / haze RT paths on or off (`do_exo_clouds`, `do_exo_haze`) | namelist |
| | Per-column MCICA seed toggle (`mcica_percol_seed`) | namelist |
| | Vertical level count (`exo_pver`) | compile time |
| **Per column** | P/T/z profile, `ts`, `ps`, all gas MMRs | input file |
| | Albedos, `coszrs`, `mw`, `cp`, `srf_emiss` | input file |
| | Cloud fields (`cicewp`, `cliqwp`, `cfrc`, `rei`, `rel`, `cicewp_co2`, `rei_co2`) | input file |
| | Haze `carmammr(pver,nelem,nbin)` | input file |
| | Gravity override (`grav`; absent or ≤0 → namelist `exo_g`) | input file |
| | Stellar constant override (`scon`; absent or ≤0 → namelist `shr_const_scon`) | input file |

To sweep over *stars* (different spectrum shapes), run one batch per star as
an outer loop — per-column `scon` only rescales the process star's spectrum
uniformly.

---

## Path A — standalone executable

### 1. Build (OpenMP is on by default)

```bash
cd build
USER_FC=gfortran make exort        # gfortran (Apple Silicon); ifort is the Linux default
```

The Makefile appends `-fopenmp` (gfortran) / `-qopenmp` (ifort) automatically.
To build without OpenMP: `make OMPFLAGS= exort`. All columns must have the
level count compiled in — check `exo_pver` in `source/exoplanet_mod.F90` and
rebuild if your profiles differ.

### 2. Make the single-column input files

Each column starts life as a normal single-column `RTprofile_in`-style file
(any existing deck works). To generate them:

```bash
cd tools
python makeColumn.py --profile TS300K --co2vmr 400e-6 --output col1.nc
python makeColumn.py --profile TS250K --co2vmr 0.01   --output col2.nc

# per-column gravity / insolation overrides (Stage E2b):
python makeColumn.py --profile TS300K --grav 3.711 --write-grav \
                     --scon 451.166 --output col3_mars.nc
```

`--write-grav` stores the same gravity used for the hypsometric `zint`, so the
column stays self-consistent. Columns without `grav`/`scon` variables use the
namelist values. Cloud/haze fields are added via the `CLOUDS` dict /
`add_haze_layer` helpers in `makeColumn.py` (see its docstring).

All files must share the same `pver` (= the compiled `exo_pver`).

### 3. Stack them into one multi-column input

```bash
python tools/stackColumns.py col1.nc col2.nc col3_mars.nc -o run/RTprofile_in.nc
```

This creates the `ncol` dimension and unions optional variables across inputs
(a file lacking a gas/cloud field contributes zeros; lacking `grav`/`scon`
contributes "use the namelist value"). Column order in the output file = the
argument order.

### 4. Write the runtime namelist

```bash
cp iofiles/user_nl_exort.template run/user_nl_exort
```

Edit `run/user_nl_exort` — this is the shared, process-level config:

```fortran
&exort_config
  solar_file     = 'G2V_SUN_n84.nc',
  shr_const_scon = 680.0,          ! default for columns without 'scon'
  exo_g          = 9.80616,        ! default for columns without 'grav'
  do_exo_clouds  = .false.,
  do_exo_haze    = .false.,
  mcica_percol_seed = .false.,     ! .true. decorrelates MCICA across cloudy columns
/
```

Note on `mcica_percol_seed`: only relevant when `do_exo_clouds=.true.` and
columns carry H₂O condensate + `cfrc`. Off (default), every column draws the
same stochastic subcolumns; on, each column's seed is offset by its column
index (column 1 is unchanged either way).

### 5. Run with threads

```csh
cd run
setenv OMP_NUM_THREADS 8        # tcsh;  bash/zsh: export OMP_NUM_THREADS=8
./exort.exe
```

- Threads spawn **only when `ncol > 1`** — single-column runs always execute
  serially on the master thread, so a threaded build changes nothing for them.
- Unset, `OMP_NUM_THREADS` defaults to all cores. No `OMP_STACKSIZE` is
  needed (the large solver work arrays are heap-allocated as of Stage E2).
- macOS loader note: if the executable can't find the NetCDF dylibs, point the
  loader at them, e.g. `setenv DYLD_LIBRARY_PATH /opt/anaconda3/lib` (use the
  directory `nf-config --flibs` refers to).

The startup banner echoes the active namelist and per-column variables found;
a multi-column threaded run prints:

```
 OpenMP: solving 3 columns on 8 threads
```

and the closing block reports wall-clock timing (`aerad_driver` is the column
solve; expect it to shrink roughly with the thread count — an 8-column batch
measured 4.4× on 8 threads).

### 6. Read the output

`run/RTprofile_out.nc` mirrors the input's `ncol` dimension: every variable
gains a trailing column axis (`LWUP(pverp,ncol)` in Fortran order =
`(ncol, pverp)` from Python), in the same column order as the stack. Without
`ncol` in the input, the output keeps the classic single-column layout.

---

## Path B — library (Python), for programmatic sweeps

```bash
cd build && USER_FC=gfortran make libexort     # -> run/libexort.dylib|.so
```

```csh
setenv OMP_NUM_THREADS 8     # BEFORE starting Python — the OpenMP runtime reads it once, at load
```

```python
import sys
sys.path.insert(0, "tools/exort_pytools")
from exort_api import ExoRT, load_rtprofile

rt = ExoRT()
rt.init(solar_file="G2V_SUN_n84.nc", scon=680.0, g=9.80616)   # once per process

states = [load_rtprofile(f, rt.dims) for f in ("col1.nc", "col2.nc", "col3_mars.nc")]
states[2]["grav"], states[2]["scon"] = 3.711, 451.166   # per-column overrides (0 = process value)

# rt.set_percol_seed(True)          # optional: decorrelate MCICA across cloudy columns
results = rt.run_columns(states)    # OpenMP-parallel over columns
```

`exort_run_columns` threads internally exactly like the executable
(`if(ncols > 1)`); each entry of `results` is a dict of numpy arrays
(`lw_upflux`, `sw_dtdt`, spectral fluxes, …). One `init()` per process; use
the `_n84` stellar filenames. See `tools/exort_pytools/README.md` for the
full API, and `verify_lib.py` / `tests/lib` for worked examples.

---

## Verifying a setup

```bash
cd tests/regression
USER_FC=gfortran python multicol_check.py        # batch == singles exactly; 8 threads == 1 thread exactly
USER_FC=gfortran python percol_config_check.py   # per-column grav/scon == namelist-configured singles
```

Both assume `run/exort.exe` is already built and finish in under a minute.

## Troubleshooting

- **"levels" mismatch / immediate stop at input:** the input `pver` differs
  from the compiled `exo_pver` (`source/exoplanet_mod.F90`) — rebuild.
- **A gas/cloud field seems ignored:** check the startup banner — every
  optional variable prints `found` or `not found, set to zero`. Common cause:
  the variable was absent from one of the stacked files (stacking zero-fills
  it for that column only).
- **H₂O clouds do nothing:** `do_exo_clouds=.true.` is required AND the
  columns need `cfrc` > 0 alongside the condensate fields — condensate with
  zero cloud fraction produces no cloud.
- **Batch slower than expected:** columns are distributed
  `schedule(dynamic)`, so a few expensive columns (cloudy/hazy) can dominate;
  more columns per batch amortize better than few.
