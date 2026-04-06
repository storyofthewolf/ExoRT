"""
makeStellarSpectrum_blackbody.py

Author: Wolf, E.T. (original IDL)
Converted to Python using Claude (Sonnet 4.6)
I/O added by Jacob Haqq-Misra

DESCRIPTION:
    Creates solar spectral files for ExoRT.
    Creates blackbody spectrum using a set temperature.
"""

import numpy as np
from scipy.interpolate import UnivariateSpline
import netCDF4 as nc
import matplotlib.pyplot as plt
import matplotlib.backends.backend_pdf as pdf
import argparse

# input arguments and options
parser = argparse.ArgumentParser(description="makeStellarSpectrum_blackbody.py: Creates solar spectral files for ExoRT. Creates blackbody spectrum using a set temperature.",
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('-t'          , type=int,   default='3000', help='Stellar effective temperature [K]')
parser.add_argument('-r'          , type=int,   default='68', help='Spectral resolution [choices: 28, 42, 68*, 73, 84]')
parser.add_argument('--plot',     action='store_true', default=True, help='plot blackbody curve and binned SED')
parser.add_argument('--netcdf',     action='store_true', default=True, help='write netcdf SED file for ExoRT')
parser.add_argument('--text',     action='store_true', default=True, help='write text data file of blackbody curve')
args = parser.parse_args()

# -------------------------
# -- stellar temperature --
# -------------------------
temperature     = int(args.t)
# -------------------------

print(f"Making blackbody spectrum: {temperature} K")

radTEMP = temperature
outname = f"blackbody_{radTEMP}K_n68.nc"

# -- choose one and only one --
# -- spectral resolution --
resolution = int(args.r)
do_n28 = False
do_n42 = False
do_n68 = False
do_n73 = False
do_n84 = False
if resolution==28:
	do_n28 = True
elif resolution==42:
	do_n42 = True
elif resolution==68:
	do_n68 = True
elif resolution==73:
	do_n73 = True
elif resolution==84:
	do_n84 = True
else:
	print(f"Unknown resolution mode, defaulting to n68")
	resolution=68
	do_n68 = True
print("Using spectral resolution n" + str(resolution))

do_plot          = args.plot   # plot blackbody curve and binned SED
do_write_netcdf  = args.netcdf   # write netcdf SED file for ExoRT
do_write_bb_dat  = args.text  # write text data file of blackbody curve


# -- supported spectral resolutions --

# rt = 28
if do_n28:
    rtwavlow = np.array([
        10.,   350.,  500.,  630.,  700.,   820.,   980.,   1100., 1180., 1390.,
        1480., 1800., 2080., 2200., 2380.,  2600.,  3250.,  4000., 4650.,
        5150., 6150., 7700., 8050., 12850., 16000., 22650., 29000., 38000.
    ])
    rtwavhi = np.array([
        350.,   500.,  630.,  700.,  820.,   980.,   1100.,  1180., 1390.,
        1480.,  1800., 2080., 2200., 2380.,  2600.,  3250.,  4000., 4650.,
        5150.,  6150., 7700., 8050., 12850., 16000., 22650., 29000., 38000., 50000.
    ])

# rt = 42
if do_n42:
    rtwavlow = np.array([
        10.0,    200.0,   350.0,   425.0,
        500.0,   630.0,   700.0,   820.0,
        980.0,   1100.0,  1180.0,  1390.0,
        1480.0,  1800.0,  2080.0,  2200.0,
        2380.0,  2600.0,  3250.0,  4000.0,
        4650.0,  4900.0,  5150.0,  5650.0,
        6150.0,  6650.0,  7000.0,  7700.0,
        8050.0,  9100.0,  10000.0, 11000.0,
        11800.0, 12850.0, 13450.0, 14450.0,
        15150.0, 16000.0, 19300.0, 22650.0,
        29000.0, 38000.0
    ])
    rtwavhi = np.array([
        200.0,   350.0,   425.0,
        500.0,   630.0,   700.0,   820.0,
        980.0,   1100.0,  1180.0,  1390.0,
        1480.0,  1800.0,  2080.0,  2200.0,
        2380.0,  2600.0,  3250.0,  4000.0,
        4650.0,  4900.0,  5150.0,  5650.0,
        6150.0,  6650.0,  7000.0,  7700.0,
        8050.0,  9100.0,  10000.0, 11000.0,
        11800.0, 12850.0, 13450.0, 14450.0,
        15150.0, 16000.0, 19300.0, 22650.0,
        29000.0, 38000.0, 50000.0
    ])

# rt = 68
if do_n68:
    rtwavlow = np.array([
        1.000,   40.00000,  100.0000,  160.0000,
        220.0000,  280.0000,  330.0000,  380.0000,
        440.0000,  495.0000,  545.0000,  617.0000,
        667.0000,  720.0000,  800.0000,  875.0000,
        940.0000,  1000.000,  1065.000,  1108.000,
        1200.000,  1275.000,  1350.000,  1450.000,
        1550.000,  1650.000,  1750.000,  1850.000,
        1950.000,  2050.000,  2200.000,  2397.000,
        2494.000,  2796.000,  3087.000,  3425.000,
        3760.000,  4030.000,  4540.000,  4950.000,
        5370.000,  5925.000,  6390.000,  6990.000,
        7650.000,  8315.000,  8850.000,  9350.000,
        9650.000,  10400.00,  11220.00,  11870.00,
        12790.00,  13300.00,  14470.00,  15000.00,
        16000.00,  16528.00,  17649.00,  18198.00,
        18518.00,  22222.00,  25641.00,  29308.00,
        30376.00,  32562.00,  35087.00,  36363.00
    ])
    rtwavhi = np.array([
        40.00000,  100.0000,  160.0000,
        220.0000,  280.0000,  330.0000,  380.0000,
        440.0000,  495.0000,  545.0000,  617.0000,
        667.0000,  720.0000,  800.0000,  875.0000,
        940.0000,  1000.000,  1065.000,  1108.000,
        1200.000,  1275.000,  1350.000,  1450.000,
        1550.000,  1650.000,  1750.000,  1850.000,
        1950.000,  2050.000,  2200.000,  2397.000,
        2494.000,  2796.000,  3087.000,  3425.000,
        3760.000,  4030.000,  4540.000,  4950.000,
        5370.000,  5925.000,  6390.000,  6990.000,
        7650.000,  8315.000,  8850.000,  9350.000,
        9650.000,  10400.00,  11220.00,  11870.00,
        12790.00,  13300.00,  14470.00,  15000.00,
        16000.00,  16528.00,  17649.00,  18198.00,
        18518.00,  22222.00,  25641.00,  29308.00,
        30376.00,  32562.00,  35087.00,  36363.00,
        42087.00
    ])

# rt = 73
if do_n73:
    rtwavlow = np.array([
        1.00000, 40.0000, 100.000, 160.000, 220.000, 280.000, 330.000,
        380.000, 440.000, 495.000, 545.000, 617.000, 667.000, 720.000,
        800.000, 875.000, 940.000, 1000.00, 1065.00, 1108.00, 1200.00,
        1275.00, 1350.00, 1450.00, 1550.00, 1650.00, 1750.00, 1850.00,
        1950.00, 2050.00, 2200.00, 2439.02, 2564.10, 2777.78, 3174.60,
        3508.77, 3773.59, 4081.63, 4545.46, 4716.98, 5154.64, 5376.34,
        5555.56, 5952.38, 6172.84, 6578.95, 6711.41, 6849.31, 7042.25,
        7462.69, 7692.31, 8064.52, 8333.33, 8620.69, 8928.57, 9090.91,
        9259.26, 9708.74, 10869.6, 11111.1, 11363.6, 11494.3, 12500.0,
        12820.5, 14492.8, 16393.4, 18181.8, 20000.0, 22222.2, 23809.5,
        25974.0, 28985.5, 33333.3
    ])
    rtwavhi = np.array([
        40.0000, 100.000, 160.000, 220.000, 280.000, 330.000,
        380.000, 440.000, 495.000, 545.000, 617.000, 667.000, 720.000,
        800.000, 875.000, 940.000, 1000.00, 1065.00, 1108.00, 1200.00,
        1275.00, 1350.00, 1450.00, 1550.00, 1650.00, 1750.00, 1850.00,
        1950.00, 2050.00, 2200.00, 2439.02, 2564.10, 2777.78, 3174.60,
        3508.77, 3773.59, 4081.63, 4545.46, 4716.98, 5154.64, 5376.34,
        5555.56, 5952.38, 6172.84, 6578.95, 6711.41, 6849.31, 7042.25,
        7462.69, 7692.31, 8064.52, 8333.33, 8620.69, 8928.57, 9090.91,
        9259.26, 9708.74, 10869.6, 11111.1, 11363.6, 11494.3, 12500.0,
        12820.5, 14492.8, 16393.4, 18181.8, 20000.0, 22222.2, 23809.5,
        25974.0, 28985.5, 33333.3, 50000.0
    ])

# rt = 84
if do_n84:
    rtwavlow = np.array([
        1.000,   40.00000,   100.0000,   160.0000,
        220.0000,   280.0000,   330.0000,   380.0000,
        440.0000,   495.0000,   545.0000,   617.0000,
        667.0000,   720.0000,   800.0000,   875.0000,
        940.0000,   1000.000,   1065.000,   1108.000,
        1200.000,   1275.000,   1350.000,   1450.000,
        1550.000,   1650.000,   1750.000,   1850.000,
        1950.000,   2050.000,   2200.000,   2397.000,
        2494.000,   2796.000,   3087.000,   3425.000,
        3760.000,   4030.000,   4540.000,   4950.000,
        5370.000,   5925.000,   6390.000,   6990.000,
        7650.000,   8315.000,   8850.000,   9350.000,
        9650.000,   10400.00,   11220.00,   11870.00,
        12790.00,   13300.00,   14470.00,   15000.00,
        16000.00,   16528.00,   17649.00,   18198.00,
        18518.00,   22222.00,   25641.00,   29308.00,
        30376.00,   32562.00,   35087.00,   36363.00,
        42087.00,   50000.00,   60000.00,   70000.00,
        80000.00,   90000.00,  100000.00,  125000.00,
        150000.00, 175000.00, 200000.00,   300000.00,
        400000.00, 500000.00, 750000.00,  1000000.00
    ])
    rtwavhi = np.array([
        40.00000,   100.0000,   160.0000,
        220.0000,   280.0000,   330.0000,   380.0000,
        440.0000,   495.0000,   545.0000,   617.0000,
        667.0000,   720.0000,   800.0000,   875.0000,
        940.0000,   1000.000,   1065.000,   1108.000,
        1200.000,   1275.000,   1350.000,   1450.000,
        1550.000,   1650.000,   1750.000,   1850.000,
        1950.000,   2050.000,   2200.000,   2397.000,
        2494.000,   2796.000,   3087.000,   3425.000,
        3760.000,   4030.000,   4540.000,   4950.000,
        5370.000,   5925.000,   6390.000,   6990.000,
        7650.000,   8315.000,   8850.000,   9350.000,
        9650.000,   10400.00,   11220.00,   11870.00,
        12790.00,   13300.00,   14470.00,   15000.00,
        16000.00,   16528.00,   17649.00,   18198.00,
        18518.00,   22222.00,   25641.00,   29308.00,
        30376.00,   32562.00,   35087.00,   36363.00,
        42087.00,   50000.00,   60000.00,   70000.00,
        80000.00,   90000.00,  100000.00,  125000.00,
        150000.00, 175000.00,  200000.00,  300000.00,
        400000.00, 500000.00,  750000.00,  1000000.00,
        1250000.00
    ])

nrtwavl = len(rtwavlow)

# Convert wavenumber (cm^-1) to wavelength (microns)
# hi/low are crossed intentionally
rtwvlhi  = 1.0e4 / rtwavlow
rtwvllow = 1.0e4 / rtwavhi
rtwvlmid = (rtwvlhi + rtwvllow) / 2.0
rtwvldel = 1.0e4 / rtwavlow - 1.0e4 / rtwavhi

Snorm = 1360.0
S0    = Snorm

# Build high-resolution wavenumber and wavelength grids
npts = 60000
dwm  = 1.0
wnm   = 1.0 + np.arange(npts) * dwm   # wavenumber [cm^-1]
wlgth = 1.0e4 / wnm                   # wavelength [microns]

# Blackbody computation (Planck function per unit wavelength)
# c1 [W m^-2 sr^-1 um^-4],  c2 [K um]
c1 = 1.1911e8
c2 = 1.4388e4
sunm = c1 / (wlgth**5 * (np.exp(c2 / (wlgth * temperature)) - 1.0))

# Build a spline interpolator once (replaces IDL interpol /spline inside loop)
spline_interp = UnivariateSpline(wnm, sunm, s=0, k=3, ext=3)

# Integrate blackbody over each spectral band
dwn  = 1.0
solarflux = np.zeros(nrtwavl)

for iw in range(nrtwavl):
    wn_lw = rtwavlow[iw]
    wn_hi = rtwavhi[iw]
    # wavenumber sample points inside the band (same stepping as IDL loop)
    w_arr = np.arange(wn_lw, wn_hi, dwn)
    # dlambda for each wavenumber step
    dwl_arr = 1.0e4 / w_arr - 1.0e4 / (w_arr + dwn)
    solarflux[iw] = np.sum(spline_interp(w_arr) * dwl_arr)

# Normalise to solar constant
ScaleFac = Snorm / np.sum(solarflux)
solarflux *= ScaleFac

# -- Plot --
if do_plot:
    print("plotting blackbody_spectrum.png")

    # Build bar arrays (step plot of binned SED)
    xbar1 = np.empty(nrtwavl * 2)
    ybar1 = np.empty(nrtwavl * 2)
    for i in range(nrtwavl):
        xbar1[2*i]   = 1.0e4 / rtwavlow[i]
        xbar1[2*i+1] = 1.0e4 / rtwavhi[i]
        ybar1[2*i]   = solarflux[i] / rtwvldel[i]
        ybar1[2*i+1] = solarflux[i] / rtwvldel[i]

    fig, ax = plt.subplots(figsize=(11, 7))
    ax.plot(wlgth, sunm * ScaleFac, color='steelblue', linewidth=1.5,
            label='Blackbody curve')
    ax.plot(xbar1, ybar1, color='darkorange', linewidth=2.5,
            label='Binned SED')
    ax.set_xlim(0.0, 5.0)
    ax.set_ylim(0, 1200)
    ax.set_xlabel('Wavelength (microns)', fontsize=11)
    ax.set_ylabel('Radiance (W m$^{-2}$ micron$^{-1}$)', fontsize=11)
    ax.set_title(f'Blackbody spectrum — {temperature} K', fontsize=12)
    ax.legend()
    fig.tight_layout()
    fig.savefig('blackbody_spectrum.png', dpi=150)
    fig.savefig('blackbody_spectrum.eps')
    plt.close(fig)

# -- Write text data file --
if do_write_bb_dat:
    outname_dat = f"blackbody_{radTEMP}K.dat"
    print(f"writing blackbody as text file output, {outname_dat}")
    with open(outname_dat, 'w') as fh:
        fh.write("wavelength(um)   radiance(W/m2/um)\n")
        for bb in range(npts):
            fh.write(f"{wlgth[bb]:15.8f}{sunm[bb]*ScaleFac:15.8f}\n")

# -- Write NetCDF --
if do_write_netcdf:
    print(f"writing output to {outname}")
    with nc.Dataset(outname, 'w', format='NETCDF4') as ds:
        ds.createDimension('nw',  nrtwavl)
        ds.createDimension('one', 1)

        v1 = ds.createVariable('wav_low',   'f4', ('nw',))
        v2 = ds.createVariable('wav_high',  'f4', ('nw',))
        v3 = ds.createVariable('S0',        'f4', ('one',))
        v4 = ds.createVariable('solarflux', 'f4', ('nw',))

        v1.title = "wavenumber edge, low";  v1.units = "cm -1"
        v2.title = "wavenumber edge, high"; v2.units = "cm -1"
        v3.title = "solar constant";        v3.units = "W/m2"
        v4.title = "spectrally integrated solar flux"; v4.units = "W/m2"

        v1[:] = rtwavlow
        v2[:] = rtwavhi
        v3[:] = S0
        v4[:] = solarflux
