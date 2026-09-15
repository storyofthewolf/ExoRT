#!/usr/bin/env python3
"""Mie scattering for homogeneous spheres — Python port of tools/mie_single.pro.

Implements the Bohren & Huffman (1983) algorithm with downward recurrence for
the logarithmic derivative, matching the EODG `mie_single` IDL routine that
generated the original CARMA optics tables.

Conventions (identical to mie_single.pro):
  * size parameter  x = 2*pi*r/lambda
  * refractive index m = n - i*k   (NEGATIVE imaginary part for absorption)

Validated against the Wiscombe (1979) MIEV0 benchmark cases; see
`python tools/mie.py --selftest`.
"""
import numpy as np


def mie_single(x, m):
    """Mie efficiencies for one size parameter and one complex refractive index.

    Parameters
    ----------
    x : float
        Size parameter 2*pi*r/lambda.
    m : complex
        Refractive index as n - 1j*k (absorption => negative imaginary part).

    Returns
    -------
    qext, qsca, g : float
        Extinction efficiency, scattering efficiency, asymmetry parameter.
    """
    x = float(x)
    if x <= 0.0:
        return 0.0, 0.0, 0.0

    # Work in the B&H convention m = n + i*k internally.
    mm = complex(m.real, abs(m.imag))

    # Series truncation (Wiscombe 1980), with a few extra terms so that very
    # large size parameters (x > ~1000, reached by the biggest CARMA bins in
    # the far UV) stay converged to round-off.
    nstop = int(round(x + 4.0 * x ** (1.0 / 3.0) + 2.0)) + 10
    nmx = int(round(max(nstop, abs(mm * x)) + 32.0))

    # Downward recurrence for D_n(mx) = psi'_n(mx)/psi_n(mx).
    # d[n] holds D_n, so term n below reads d[n] (not d[n-1]).
    d = np.zeros(nmx + 1, dtype=complex)
    y = mm * x
    for n in range(nmx, 0, -1):
        d[n - 1] = (n / y) - 1.0 / (d[n] + n / y)

    # Riccati-Bessel functions by upward recurrence.
    # psi_prev/chi_prev hold order n-1, psi_prev2/chi_prev2 order n-2.
    psi_prev2 = np.cos(x)    # psi_{-1}
    psi_prev = np.sin(x)     # psi_0
    chi_prev2 = -np.sin(x)   # chi_{-1}
    chi_prev = np.cos(x)     # chi_0

    qsca = 0.0
    qext = 0.0
    gsum = 0.0
    an1 = 0.0 + 0.0j   # a_{n-1}
    bn1 = 0.0 + 0.0j   # b_{n-1}

    for n in range(1, nstop + 1):
        psi = (2.0 * n - 1.0) * psi_prev / x - psi_prev2
        chi = (2.0 * n - 1.0) * chi_prev / x - chi_prev2
        xi = complex(psi, -chi)
        xi_prev = complex(psi_prev, -chi_prev)

        dn = d[n]
        ta = dn / mm + n / x
        tb = dn * mm + n / x
        an = (ta * psi - psi_prev) / (ta * xi - xi_prev)
        bn = (tb * psi - psi_prev) / (tb * xi - xi_prev)

        tn = 2.0 * n + 1.0
        qsca += tn * (abs(an) ** 2 + abs(bn) ** 2)
        qext += tn * (an.real + bn.real)

        # Asymmetry parameter (B&H eq. 4.62): cross term uses a_{n-1}, b_{n-1}.
        if n > 1:
            nm = n - 1.0
            gsum += (nm * (nm + 2.0) / (nm + 1.0)) * (
                an1 * np.conj(an) + bn1 * np.conj(bn)
            ).real
        gsum += (tn / (n * (n + 1.0))) * (an * np.conj(bn)).real

        an1, bn1 = an, bn
        psi_prev2, psi_prev = psi_prev, psi
        chi_prev2, chi_prev = chi_prev, chi

    qsca *= 2.0 / (x * x)
    qext *= 2.0 / (x * x)

    g = 0.0 if qsca <= 0.0 else (4.0 / (x * x)) * gsum / qsca

    # Guard against round-off excursions outside physical bounds.
    if qsca > qext:
        qsca = qext
    g = min(max(g, -1.0), 1.0)

    return qext, qsca, g


def mie_arrays(x_arr, m):
    """Convenience wrapper over `mie_single` for an array of size parameters."""
    x_arr = np.atleast_1d(np.asarray(x_arr, dtype=float))
    qext = np.empty(x_arr.shape)
    qsca = np.empty(x_arr.shape)
    g = np.empty(x_arr.shape)
    for i in range(x_arr.size):
        qext.flat[i], qsca.flat[i], g.flat[i] = mie_single(x_arr.flat[i], m)
    return qext, qsca, g


# Benchmark cases spanning the Rayleigh, resonance and geometric-optics
# regimes. Reference values independently computed with miepython 3.2.0.
# (x, m, qext, qsca, g)
_BENCHMARKS = [
    (0.10, complex(1.5, 0.0), 0.00002313, 0.00002313, 0.00198157),
    (0.50, complex(1.7, -0.3), 0.32861082, 0.03046713, 0.05346096),
    (1.00, complex(1.5, 0.0), 0.21509761, 0.21509761, 0.19894237),
    (5.00, complex(1.33, 0.0), 3.59103285, 3.59103285, 0.84534027),
    (10.00, complex(1.5, 0.0), 2.88199906, 2.88199906, 0.74291310),
    (10.00, complex(1.5, -0.1), 2.45979049, 1.23514423, 0.92234953),
    (50.00, complex(2.1, -0.5), 2.14597246, 1.25394728, 0.86820318),
    (100.00, complex(1.5, 0.0), 2.09438784, 2.09438784, 0.81824583),
    (100.00, complex(1.5, -0.1), 2.08982179, 1.13213400, 0.95039208),
    (1000.00, complex(1.6, -0.01), 2.01984414, 1.11949418, 0.94054224),
]


def _selftest():
    ok = True
    print(f"{'x':>8} {'m':>14} {'Qext':>12} {'ref':>12} {'Qsca':>12} {'ref':>12}"
          f" {'g':>10} {'ref':>10}")
    for x, m, qe, qs, gg in _BENCHMARKS:
        a, b, c = mie_single(x, m)
        print(f"{x:8.2f} {str(m):>14} {a:12.6f} {qe:12.6f} {b:12.6f} {qs:12.6f}"
              f" {c:10.6f} {gg:10.6f}")
        if not (abs(a - qe) < 2e-5 and abs(b - qs) < 2e-5 and abs(c - gg) < 2e-5):
            ok = False
    print("SELFTEST", "PASS" if ok else "FAIL")
    return 0 if ok else 1


if __name__ == "__main__":
    import sys
    sys.exit(_selftest())
