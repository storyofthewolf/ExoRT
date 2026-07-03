/*
 * test_exort_c.c -- standalone C smoke/regression test for libexort
 * (ExoRT Stage D library API).
 *
 * Loads the shared library with dlopen, reads a regression fixture with
 * the netCDF C API, runs one column, and compares every output against
 * the committed golden baseline for the same case (TS300K_G2V). Also
 * checks that the batch entry point (exort_run_columns) reproduces the
 * single-column result bit-for-bit.
 *
 * The column_state_t / column_result_t structs are all-double, so this
 * test addresses them as flat double arrays with offsets computed from
 * exort_get_dims() -- it needs no recompile when exo_pver or the band
 * grid changes. Field order must match exort_column_mod.F90.
 *
 * The baseline NetCDF stores float32, so the comparison tolerance is
 * storage-limited (rtol 1e-5, atol 1e-4); any wiring error is orders of
 * magnitude larger.
 *
 * Build and run (see Makefile in this directory):
 *     cd tests/lib && make && make run
 */

#include <dlfcn.h>
#include <math.h>
#include <netcdf.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define DEFAULT_LIB      "../../run/libexort" SHLIB_EXT
#define DEFAULT_ROOT     "../../"
#define DEFAULT_FIXTURE  "../regression/fixtures/RTprofile_in_TS300K.nc"
#define DEFAULT_BASELINE "../regression/baselines/RTprofile_out_TS300K_G2V.nc"
#define SOLAR_FILE       "G2V_SUN_n84.nc"  /* baselines are n84-grid runs */
#define SCON             680.0
#define GRAV             9.80616

#define RTOL 1e-5
#define ATOL 1e-4
#define CSEC 86400.0

typedef int (*fn_get_dims)(int *, int *, int *, int *, int *);
typedef int (*fn_init)(const char *, const char *, double, double, int, int);
typedef int (*fn_run_column)(const double *, double *);
typedef int (*fn_run_columns)(int, const double *, double *);

static int n_fail = 0;

static void check(int ok, const char *what, double maxd) {
    printf("  %-24s %s  max|d|=%.3e\n", what, ok ? "PASS" : "FAIL", maxd);
    if (!ok) n_fail++;
}

/* read a netCDF variable into buf; if absent and !required, fill with fill */
static void read_var(int ncid, const char *name, double *buf, size_t n,
                     int required, double fill) {
    int varid;
    if (nc_inq_varid(ncid, name, &varid) != NC_NOERR) {
        if (required) {
            fprintf(stderr, "missing required variable '%s'\n", name);
            exit(2);
        }
        for (size_t i = 0; i < n; i++) buf[i] = fill;
        return;
    }
    if (nc_get_var_double(ncid, varid, buf) != NC_NOERR) {
        fprintf(stderr, "error reading variable '%s'\n", name);
        exit(2);
    }
}

/* |a-b| <= ATOL + RTOL*|b| elementwise; returns max |a-b| */
static double maxdiff(const double *a, const double *b, size_t n, int *ok) {
    double m = 0.0;
    *ok = 1;
    for (size_t i = 0; i < n; i++) {
        double d = fabs(a[i] - b[i]);
        if (d > m) m = d;
        if (d > ATOL + RTOL * fabs(b[i])) *ok = 0;
    }
    return m;
}

int main(int argc, char **argv) {
    const char *libpath  = argc > 1 ? argv[1] : DEFAULT_LIB;
    const char *root     = argc > 2 ? argv[2] : DEFAULT_ROOT;
    const char *fixture  = argc > 3 ? argv[3] : DEFAULT_FIXTURE;
    const char *baseline = argc > 4 ? argv[4] : DEFAULT_BASELINE;

    void *h = dlopen(libpath, RTLD_NOW);
    if (!h) { fprintf(stderr, "dlopen: %s\n", dlerror()); return 2; }
    fn_get_dims    get_dims    = (fn_get_dims)dlsym(h, "exort_get_dims");
    fn_init        init        = (fn_init)dlsym(h, "exort_init");
    fn_run_column  run_column  = (fn_run_column)dlsym(h, "exort_run_column");
    fn_run_columns run_columns = (fn_run_columns)dlsym(h, "exort_run_columns");
    if (!get_dims || !init || !run_column || !run_columns) {
        fprintf(stderr, "dlsym: %s\n", dlerror());
        return 2;
    }

    int pver, pverp, nwave, nelem, nbin;
    if (get_dims(&pver, &pverp, &nwave, &nelem, &nbin) != 0) return 2;
    printf("dims: pver=%d pverp=%d nwave=%d nelem=%d nbin=%d\n",
           pver, pverp, nwave, nelem, nbin);

    /* struct sizes in doubles (field order = exort_column_mod.F90) */
    size_t ncarma = (size_t)pver * nelem * nbin;
    /* 10 leading scalars + trailing grav/scon (appended, Stage E2b) */
    size_t nstate = 12 + 20 * (size_t)pver + 3 * (size_t)pverp + ncarma;
    size_t nresult = 5 + 2 * (size_t)pver + 4 * (size_t)pverp
                     + 4 * (size_t)pverp * nwave;
    double *state = calloc(nstate, sizeof(double));
    double *res   = calloc(nresult, sizeof(double));
    double *res2  = calloc(2 * nresult, sizeof(double));
    double *st2   = calloc(2 * nstate, sizeof(double));
    if (!state || !res || !res2 || !st2) return 2;

    /* ---- fill the state from the fixture ---- */
    int ncid;
    if (nc_open(fixture, NC_NOWRITE, &ncid) != NC_NOERR) {
        fprintf(stderr, "cannot open fixture %s\n", fixture);
        return 2;
    }
    double *p = state;
    read_var(ncid, "ts",        p++, 1, 1, 0);
    read_var(ncid, "ps",        p++, 1, 1, 0);
    read_var(ncid, "coszrs",    p++, 1, 1, 0);
    read_var(ncid, "mw",        p++, 1, 1, 0);
    read_var(ncid, "cp",        p++, 1, 1, 0);
    read_var(ncid, "srf_emiss", p++, 1, 0, 1.0);
    read_var(ncid, "asdir",     p++, 1, 1, 0);
    read_var(ncid, "asdif",     p++, 1, 1, 0);
    read_var(ncid, "aldir",     p++, 1, 1, 0);
    read_var(ncid, "aldif",     p++, 1, 1, 0);
    read_var(ncid, "tmid", p, pver,  1, 0); p += pver;
    read_var(ncid, "tint", p, pverp, 1, 0); p += pverp;
    read_var(ncid, "pmid", p, pver,  1, 0); p += pver;
    read_var(ncid, "pdel", p, pver,  1, 0); p += pver;
    read_var(ncid, "pint", p, pverp, 1, 0); p += pverp;
    read_var(ncid, "zint", p, pverp, 1, 0); p += pverp;
    const char *gases[] = {"h2ommr", "co2mmr", "ch4mmr", "c2h6mmr", "nh3mmr",
                           "commr", "o2mmr", "o3mmr", "n2mmr", "h2mmr"};
    for (int i = 0; i < 10; i++) { read_var(ncid, gases[i], p, pver, 0, 0); p += pver; }
    const char *clouds[] = {"cicewp", "cliqwp", "cfrc", "rei", "rel",
                            "cicewp_co2", "rei_co2"};
    for (int i = 0; i < 7; i++) { read_var(ncid, clouds[i], p, pver, 0, 0); p += pver; }
    read_var(ncid, "carmammr", p, ncarma, 0, 0); p += ncarma;
    /* per-column gravity / stellar constant; 0 = use the exort_init values */
    read_var(ncid, "grav", p++, 1, 0, 0.0);
    read_var(ncid, "scon", p++, 1, 0, 0.0);
    nc_close(ncid);
    if ((size_t)(p - state) != nstate) {
        fprintf(stderr, "state packing error\n");
        return 2;
    }

    /* ---- init + run ---- */
    if (init(root, SOLAR_FILE, SCON, GRAV, 0, 0) != 0) {
        fprintf(stderr, "exort_init failed\n");
        return 2;
    }
    if (run_column(state, res) != 0) {
        fprintf(stderr, "exort_run_column failed\n");
        return 2;
    }

    /* ---- compare against the baseline ---- */
    if (nc_open(baseline, NC_NOWRITE, &ncid) != NC_NOERR) {
        fprintf(stderr, "cannot open baseline %s\n", baseline);
        return 2;
    }
    printf("\n=== library vs baseline %s (rtol=%g, atol=%g) ===\n",
           baseline, RTOL, ATOL);

    /* result offsets (field order = exort_column_mod.F90) */
    double *sw_dtdt = res + 5;
    double *lw_dtdt = sw_dtdt + pver;
    double *lw_dn = lw_dtdt + pver;
    double *lw_up = lw_dn + pverp;
    double *sw_dn = lw_up + pverp;
    double *sw_up = sw_dn + pverp;
    double *spec[4]; /* lw_dn, lw_up, sw_dn, sw_up spectral */
    spec[0] = sw_up + pverp;
    for (int i = 1; i < 4; i++) spec[i] = spec[i - 1] + (size_t)pverp * nwave;

    size_t nbuf = (size_t)pverp * nwave;
    double *buf = malloc(nbuf * sizeof(double));
    double *tmp = malloc(nbuf * sizeof(double));
    int ok;

    struct { const char *var; double *val; size_t n; double scale; } flux[] = {
        {"LWDN", lw_dn, (size_t)pverp, 1.0}, {"LWUP", lw_up, (size_t)pverp, 1.0},
        {"SWDN", sw_dn, (size_t)pverp, 1.0}, {"SWUP", sw_up, (size_t)pverp, 1.0},
        {"LWHR", lw_dtdt, (size_t)pver, CSEC}, {"SWHR", sw_dtdt, (size_t)pver, CSEC},
        {"FSDTOA", res, 1, 1.0},
    };
    for (size_t i = 0; i < sizeof(flux) / sizeof(flux[0]); i++) {
        read_var(ncid, flux[i].var, buf, flux[i].n, 1, 0);
        for (size_t k = 0; k < flux[i].n; k++) tmp[k] = flux[i].val[k] * flux[i].scale;
        double m = maxdiff(tmp, buf, flux[i].n, &ok);
        check(ok, flux[i].var, m);
    }

    /* spectral: baseline is C-order (nwave, pverp); library result is
       Fortran-order (pverp, nwave) -> element (k,w) sits at k + pverp*w
       in the result and at w*pverp + k in the file. Same flat layout. */
    const char *svars[] = {"LWDN_SPECTRAL", "LWUP_SPECTRAL",
                           "SWDN_SPECTRAL", "SWUP_SPECTRAL"};
    for (int i = 0; i < 4; i++) {
        read_var(ncid, svars[i], buf, nbuf, 1, 0);
        double m = maxdiff(spec[i], buf, nbuf, &ok);
        check(ok, svars[i], m);
    }
    nc_close(ncid);

    /* ---- batch path must reproduce the single run bit-for-bit ---- */
    memcpy(st2, state, nstate * sizeof(double));
    memcpy(st2 + nstate, state, nstate * sizeof(double));
    if (run_columns(2, st2, res2) != 0) {
        fprintf(stderr, "exort_run_columns failed\n");
        return 2;
    }
    printf("\n=== batch determinism (exact) ===\n");
    for (int i = 0; i < 2; i++) {
        double m = 0.0;
        for (size_t k = 0; k < nresult; k++) {
            double d = fabs(res2[i * nresult + k] - res[k]);
            if (d > m) m = d;
        }
        char label[32];
        snprintf(label, sizeof(label), "run_columns[%d]", i);
        check(m == 0.0, label, m);
    }

    printf("\nLWUP_TOM=%.3f W/m2  SWDN_SFC=%.3f W/m2\n", lw_up[0], sw_dn[pverp - 1]);
    printf("%s\n", n_fail == 0 ? "PASS" : "FAIL");
    return n_fail == 0 ? 0 : 1;
}
