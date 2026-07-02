# exort_pytools — Python binding for libexort (Stage D)

Run ExoRT columns in-process from Python, without `RTprofile_in.nc` round-trips
through the executable.

## Build the library

```bash
cd build
USER_FC=gfortran make libexort        # -> run/libexort.dylib (macOS) / .so (Linux)
```

`exo_pver` (vertical levels) is still compile-time, set in
`source/exoplanet_mod.F90`; the binding queries it at load time.

## Use

```python
from exort_api import ExoRT, load_rtprofile

rt = ExoRT()                                   # finds run/libexort.*
rt.init(solar_file="G2V_SUN_n84.nc",           # '' keeps compiled default
        scon=680.0, g=9.80616,                 # <=0 keeps compiled default
        do_clouds=False, do_haze=False)        # table loading is gated here

state = load_rtprofile("RTprofile_in_TS300K.nc", rt.dims)   # dict of arrays
result = rt.run_column(state)                  # dict of numpy arrays

result["lw_upflux"][0]        # OLR at top of model [W m-2]
result["sw_dnflux"][-1]       # surface SW down [W m-2]
result["lw_dtdt"] * 86400     # LW heating [K day-1], returned as K s-1
```

`state` can also be built by hand — any dict with the `STATE_FIELDS` keys
(`exort_api.py`); missing gas/cloud/haze fields default to zero, `srf_emiss`
to 1.0. `rt.run_columns([s1, s2, ...])` runs a batch through one call.

Constraints (Stage D):

- **one `init()` per process** (Fortran module tables are global); re-init
  and post-`finalize()` re-init are refused.
- `run_column` is serial — tables are read-only after init, but the RT
  kernels keep module scratch state. OpenMP/multi-column parallelism is
  Stage E.
- fatal data errors (missing k-files, wrong-shape inputs) `stop` inside the
  legacy Fortran readers and take the process down, exactly like the
  executable.

## Verify

```bash
python verify_lib.py                    # TS300K_G2V vs committed baseline
python verify_lib.py --case 2barCO2_dry_Mars_G2V
python verify_lib.py --list
```

Compares every output (fluxes, heating rates, per-band spectra) against the
committed regression baseline for the case, plus exact bitwise determinism
checks of repeat and batch calls. Baselines store float32, so the file
comparison is storage-limited (~1e-7 relative); the underlying doubles come
from the same compiled kernels the executable uses.

The equivalent C-level check is `tests/lib/test_exort_c.c`
(`cd tests/lib && make run`).

Requires `cffi`, `numpy`, `netCDF4`.
