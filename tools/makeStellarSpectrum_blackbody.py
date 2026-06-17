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
import argparse

from spectral_intervals import get_spectral_intervals

# input arguments and options
parser = argparse.ArgumentParser(description="makeStellarSpectrum_blackbody.py: Creates solar spectral files for ExoRT. Creates blackbody spectrum using a set temperature.",
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('-t'          , type=int,   default='3400', help='Stellar effective temperature [K]')
parser.add_argument('-r'          , type=int,   default='84', help='Spectral resolution [choices: 28, 42, 68, 84]')
parser.add_argument('--plot',     action='store_true', default=True, help='plot blackbody curve and binned SED')
parser.add_argument('--netcdf',     action='store_true', default=True, help='write netcdf SED file for ExoRT')
parser.add_argument('--text',     action='store_true', default=False, help='write text data file of blackbody curve')
args = parser.parse_args()

# -------------------------
# -- stellar temperature --
# -------------------------
temperature     = int(args.t)
# -------------------------

print(f"Making blackbody spectrum: {temperature} K")

radTEMP = temperature
resolution = int(args.r)
outname = f"blackbody_{radTEMP}K_n{resolution}.nc"

# -- spectral resolution --
if resolution not in (28, 42, 68, 84):
    print("Unknown resolution mode, defaulting to n68")
    resolution = 68
    outname = f"blackbody_{radTEMP}K_n{resolution}.nc"
tag = f"n{resolution}"
print("Using spectral resolution " + tag)

do_plot          = args.plot   # plot blackbody curve and binned SED
do_write_netcdf  = args.netcdf   # write netcdf SED file for ExoRT
do_write_bb_dat  = args.text  # write text data file of blackbody curve

# -- ExoRT spectral grid (edges -> low/high interval bounds) --
# Grid definitions are centralized in spectral_intervals.py, which already
# substitutes the longwave 0 cm-1 edge with 1.0 cm-1 (finite wavelength).
(nrtwavl, wavenum_edge, _wavenum_mid,
 _wavelength_edge, _wavelength_mid) = get_spectral_intervals(tag)
rtwavlow = wavenum_edge[:-1]
rtwavhi  = wavenum_edge[1:]

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
