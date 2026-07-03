"""
exort_api.py -- Python binding for libexort (ExoRT Stage D library API).

Wraps the C API exported by source/src.main/exort_lib_mod.F90 via cffi in
ABI mode (dlopen). The struct layout is the field order of column_state_t /
column_result_t in source/src.main/exort_column_mod.F90; the array extents
are queried from the library at load time (exort_get_dims), so one binding
works for any compiled exo_pver / band grid.

Typical use:

    from exort_api import ExoRT, load_rtprofile

    rt = ExoRT()                                  # finds run/libexort.dylib
    rt.init(scon=680.0, g=9.80616, solar_file="G2V_SUN_n68.nc")
    state = load_rtprofile("RTprofile_in.nc", rt.dims)
    result = rt.run_column(state)
    print(result["lw_upflux"][0])                 # OLR at top of model

States are plain dicts of floats / numpy arrays keyed by the Fortran field
names; results are dicts of numpy arrays (spectral arrays have shape
(pverp, nwave)). Heating rates are returned in K s-1 (multiply by 86400
for the K day-1 convention of RTprofile_out.nc).

Requires: cffi, numpy (netCDF4 only for load_rtprofile).
"""

import os
import sys

import numpy as np

try:
    import cffi
except ImportError:
    sys.exit("exort_api requires the 'cffi' package (pip install cffi)")

# repo root = two levels up from tools/exort_pytools/
_REPO_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))

_BOOT_CDEF = """
int exort_get_dims(int *pver, int *pverp, int *nwave, int *nelem, int *nbin);
int exort_init(const char *data_root, const char *solar_file,
               double scon, double g, int do_clouds, int do_haze);
int exort_set_percol_seed(int enable);
int exort_finalize(void);
"""

# error codes from exort_lib_mod.F90
EXORT_OK = 0
EXORT_ALREADY_INIT = 1
EXORT_NOT_INITIALIZED = 2
EXORT_BAD_ARGS = 3

# (name, kind) in exact Fortran declaration order; kind is 'scalar', 'pver',
# 'pverp', 'carma' (pver*nelem*nbin) or 'spectral' (pverp*nwave)
STATE_FIELDS = [
    ("ts", "scalar"), ("ps", "scalar"), ("coszrs", "scalar"),
    ("mwdry", "scalar"), ("cpdry", "scalar"), ("srf_emiss", "scalar"),
    ("asdir", "scalar"), ("asdif", "scalar"),
    ("aldir", "scalar"), ("aldif", "scalar"),
    ("tmid", "pver"), ("tint", "pverp"),
    ("pmid", "pver"), ("pdel", "pver"), ("pint", "pverp"),
    ("zint", "pverp"),
    ("h2ommr", "pver"), ("co2mmr", "pver"), ("ch4mmr", "pver"),
    ("c2h6mmr", "pver"), ("nh3mmr", "pver"), ("commr", "pver"),
    ("o2mmr", "pver"), ("o3mmr", "pver"), ("n2mmr", "pver"),
    ("h2mmr", "pver"),
    ("cicewp", "pver"), ("cliqwp", "pver"), ("cfrc", "pver"),
    ("rei", "pver"), ("rel", "pver"),
    ("cicewp_co2", "pver"), ("rei_co2", "pver"),
    ("carmammr", "carma"),
    # per-column planetary/orbital config; <= 0 (incl. the zero default)
    # means "use the process-level exort_init/namelist value"
    ("grav", "scalar"), ("scon", "scalar"),
]

RESULT_FIELDS = [
    ("sol_toa", "scalar"), ("vis_dir", "scalar"), ("vis_dif", "scalar"),
    ("nir_dir", "scalar"), ("nir_dif", "scalar"),
    ("sw_dtdt", "pver"), ("lw_dtdt", "pver"),
    ("lw_dnflux", "pverp"), ("lw_upflux", "pverp"),
    ("sw_dnflux", "pverp"), ("sw_upflux", "pverp"),
    ("lw_dnflux_spectral", "spectral"), ("lw_upflux_spectral", "spectral"),
    ("sw_dnflux_spectral", "spectral"), ("sw_upflux_spectral", "spectral"),
]


def default_lib_path():
    ext = ".dylib" if sys.platform == "darwin" else ".so"
    return os.path.join(_REPO_ROOT, "run", "libexort" + ext)


class ExoRT:
    """Handle on one loaded libexort. One init per process (Fortran module
    tables are global); create the instance, call init() once, then
    run_column / run_columns as often as needed."""

    def __init__(self, lib_path=None):
        lib_path = lib_path or default_lib_path()
        if not os.path.exists(lib_path):
            raise FileNotFoundError(
                f"{lib_path} not found -- build it with "
                "'cd build && USER_FC=gfortran make libexort'")
        boot = cffi.FFI()
        boot.cdef(_BOOT_CDEF)
        blib = boot.dlopen(lib_path)
        p = [boot.new("int *") for _ in range(5)]
        _check(blib.exort_get_dims(*p), "exort_get_dims")
        self.dims = dict(zip(["pver", "pverp", "nwave", "nelem", "nbin"],
                             (x[0] for x in p)))

        self.ffi = cffi.FFI()
        self.ffi.cdef(_BOOT_CDEF + self._struct_cdef())
        self.lib = self.ffi.dlopen(lib_path)
        self.lib_path = lib_path
        self._initialized = False

    # ---- API -------------------------------------------------------------

    def init(self, data_root=None, solar_file="", scon=-1.0, g=-1.0,
             do_clouds=False, do_haze=False):
        """Load all tables. Empty solar_file / non-positive scon,g keep the
        compiled defaults; data_root defaults to this repo checkout."""
        if data_root is None:
            data_root = _REPO_ROOT
        _check(self.lib.exort_init(data_root.encode(), solar_file.encode(),
                                   float(scon), float(g),
                                   int(bool(do_clouds)), int(bool(do_haze))),
               "exort_init")
        self._initialized = True

    def run_column(self, state):
        """Run one column. `state` is a dict keyed by STATE_FIELDS names;
        missing gas/cloud/haze fields default to zero, srf_emiss to 1.0.
        Returns a dict of numpy arrays keyed by RESULT_FIELDS names."""
        c_state = self._pack_state(state)
        c_result = self.ffi.new("column_result_t *")
        _check(self.lib.exort_run_column(c_state, c_result), "exort_run_column")
        return self._unpack_result(c_result)

    def run_columns(self, states):
        """Run a list of column states through the batch entry point."""
        n = len(states)
        c_states = self.ffi.new(f"column_state_t[{n}]")
        for i, s in enumerate(states):
            self._pack_state(s, into=self.ffi.addressof(c_states, i))
        c_results = self.ffi.new(f"column_result_t[{n}]")
        _check(self.lib.exort_run_columns(n, c_states, c_results),
               "exort_run_columns")
        return [self._unpack_result(self.ffi.addressof(c_results, i))
                for i in range(n)]

    def set_percol_seed(self, enable):
        """Opt-in per-column MCICA seed: when enabled, each column of a
        run_columns batch draws decorrelated stochastic cloud subcolumns
        (seed offset by column index). Disabled (default) = every column
        uses the same constant seed, the legacy behavior. Only affects
        cloudy columns."""
        _check(self.lib.exort_set_percol_seed(int(bool(enable))),
               "exort_set_percol_seed")

    def finalize(self):
        _check(self.lib.exort_finalize(), "exort_finalize")
        self._initialized = False

    # ---- internals ---------------------------------------------------------

    def _extent(self, kind):
        d = self.dims
        return {"scalar": 1, "pver": d["pver"], "pverp": d["pverp"],
                "carma": d["pver"] * d["nelem"] * d["nbin"],
                "spectral": d["pverp"] * d["nwave"]}[kind]

    def _struct_cdef(self):
        def fields(spec):
            out = []
            for name, kind in spec:
                n = self._extent(kind)
                out.append(f"  double {name};" if kind == "scalar"
                           else f"  double {name}[{n}];")
            return "\n".join(out)
        return (f"typedef struct {{\n{fields(STATE_FIELDS)}\n}} column_state_t;\n"
                f"typedef struct {{\n{fields(RESULT_FIELDS)}\n}} column_result_t;\n"
                "int exort_run_column(const column_state_t *state,"
                " column_result_t *result);\n"
                "int exort_run_columns(int n, const column_state_t *states,"
                " column_result_t *results);\n")

    def _pack_state(self, state, into=None):
        defaults = {"srf_emiss": 1.0}
        c_state = into if into is not None else self.ffi.new("column_state_t *")
        unknown = set(state) - {name for name, _ in STATE_FIELDS}
        if unknown:
            raise KeyError(f"unknown state fields: {sorted(unknown)}")
        for name, kind in STATE_FIELDS:
            if kind == "scalar":
                setattr(c_state, name,
                        float(state.get(name, defaults.get(name, 0.0))))
            else:
                n = self._extent(kind)
                val = state.get(name)
                if val is None:
                    arr = np.zeros(n)
                else:
                    # ravel in C order: matches the byte sequence the Fortran
                    # side reads into its column-major array (see netCDF
                    # Fortran/C dim-order reversal)
                    arr = np.ascontiguousarray(val, dtype=np.float64).ravel()
                    if arr.size != n:
                        raise ValueError(
                            f"state field '{name}' has {arr.size} elements, "
                            f"expected {n}")
                self.ffi.memmove(getattr(c_state, name), arr, n * 8)
        return c_state

    def _unpack_result(self, c_result):
        out = {}
        for name, kind in RESULT_FIELDS:
            if kind == "scalar":
                out[name] = float(getattr(c_result, name))
            else:
                n = self._extent(kind)
                arr = np.frombuffer(
                    self.ffi.buffer(getattr(c_result, name), n * 8),
                    dtype=np.float64).copy()
                if kind == "spectral":
                    # Fortran (pverp, nwave) column-major -> numpy order='F'
                    arr = arr.reshape((self.dims["pverp"], self.dims["nwave"]),
                                      order="F")
                out[name] = arr
        return out


def _check(ierr, what):
    if ierr != EXORT_OK:
        msg = {EXORT_ALREADY_INIT: "already initialized",
               EXORT_NOT_INITIALIZED: "not initialized",
               EXORT_BAD_ARGS: "bad arguments"}.get(ierr, f"error code {ierr}")
        raise RuntimeError(f"{what} failed: {msg}")


def load_rtprofile(path, dims):
    """Read an RTprofile_in.nc into a column-state dict, mirroring the
    required/optional variable handling of input_profile in io_1D.F90."""
    import netCDF4

    required = ["ts", "ps", "tmid", "tint", "pmid", "pdel", "pint", "zint",
                "asdir", "asdif", "aldir", "aldif", "coszrs"]
    optional = ["srf_emiss", "h2ommr", "co2mmr", "ch4mmr", "c2h6mmr",
                "nh3mmr", "commr", "o2mmr", "o3mmr", "n2mmr", "h2mmr",
                "cicewp", "cliqwp", "rei", "rel", "cicewp_co2", "rei_co2",
                "carmammr"]
    state = {}
    with netCDF4.Dataset(path) as ds:
        if len(ds.dimensions["pver"]) != dims["pver"]:
            raise ValueError(
                f"{path} has pver={len(ds.dimensions['pver'])}, but the "
                f"library was compiled with exo_pver={dims['pver']}")
        def read(name):
            arr = np.asarray(ds.variables[name][:], dtype=np.float64)
            return float(arr.ravel()[0]) if arr.size == 1 else arr

        for name in required:
            state[name] = read(name)
        # file names 'mw'/'cp' map to state fields mwdry/cpdry
        state["mwdry"] = read("mw")
        state["cpdry"] = read("cp")
        for name in optional:
            if name in ds.variables:
                state[name] = read(name)
    return state
