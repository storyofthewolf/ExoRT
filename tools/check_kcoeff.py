###  open and read all k-coeff data for a given gas.
###  check which bins have zero absorption
###  create grey absorption plots

import netCDF4 as nc
import numpy   as np
import os
import argparse
import sys
from spectral_intervals import get_spectral_intervals, get_g_weights
import matplotlib.pyplot as plt


# input arguments and options                                                                                                                          
parser = argparse.ArgumentParser()
parser.add_argument('gasname'     , type=str,   nargs=1, default=' ',  help='gas name')
parser.add_argument('desc'          , type=str, nargs=1, default=' ',  help='descriptors [see kcoefficent file names]')
parser.add_argument('-n'          , type=int,   default='68',          help='number of spectral intervals [28,42,68,84]')
parser.add_argument('-g'          , type=str,   default='rrtm8',       help='gauss point intervals [rrtm8, rrmt16, gqaud8]')
parser.add_argument('-l'          , type=str,   default='hitran2016',  help='line list [hitran2012, hitran2016]')
parser.add_argument('--plots',      action='store_true', help='do not do lineplots at end of sequence')
args = parser.parse_args()

# Users set parameters
gasname  = str(args.gasname[0])
desc     = str(args.desc[0])
nbins    = int(args.n) 
gtag     = str(args.g)
linelist = str(args.l)

tag      = "n" + str(nbins)
fdir     = "/discover/nobackup/etwolf/models/ExoRT/data/kdist/" + tag + gasname + "/" + linelist + "/"


# get spectral intervals
nw, wvn_e, wvn_m, wvl_e, wvl_m = get_spectral_intervals(tag)
# get gauss points
g_weights = get_g_weights(gtag)


# open first file to get dimensions
j=1
filename = "n" + str(nbins) + "_8gpt_" + "bin0" + str(j) + "_" + gasname + "_" + desc + ".nc"
openfile = fdir + filename
ncid = nc.Dataset(openfile,'r')
T     = ncid.variables['Temperature'][:]   ; T    = np.squeeze(T)
P     = ncid.variables['Pressure'][:]      ; P    = np.squeeze(P)
G     = ncid.variables['GaussWeights'][:]  ; G    = np.squeeze(G)
ncid.close()
nT = T.size
nP = P.size
nG = G.size

# matrix of grey absorption data
grey = np.zeros([nT,nP,nbins], dtype=float)

# loop over each bin
for i in range(nbins):
    j = i + 1
    if (j < 10):  
        filename = "n" + str(nbins) + "_8gpt_" + "bin0" + str(j) + "_" + gasname + "_" + desc + ".nc"
    elif (j>=10):
        filename = "n" + str(nbins) + "_8gpt_" + "bin" + str(j) + "_" + gasname + "_" + desc + ".nc"
    openfile = fdir + filename
    ncid = nc.Dataset(openfile,'r')
    Kabs  = ncid.variables['data'][:]          ; Kabs = np.squeeze(Kabs)
    ncid.close()
    x = 0.0
    for tt in range(nT):
      for pp in range(nP):
          xtemp = max(Kabs[tt,pp,:])
          x = max(x, xtemp)
          for gg in range(nG):
              grey[tt,pp,i] = grey[tt,pp,i]+Kabs[tt,pp,gg]*g_weights[gg]

    fs = "{:.3f}"
    fs1 = "{:3}"
    fs2 = "{:10.3f}"
    print(fs1.format(j), fs2.format(wvl_e[i]), fs2.format(wvl_e[j]), gasname, "max_val", x)  

    # compile grey spectral information

if args.plots == True:
    print("plotting grey band values")
    x    = wvn_m
    y1    = grey[10,10,:]
    y2    = grey[1,1,:]
    ybot = max(0.0, np.min(y1), np.min(y2))
    ytop = max(np.max(y1), np.max(y2))
    xbot = 1
    xtop = 50000
    plt.loglog(x, y1, linestyle='-', color='b', label=gasname)
    plt.loglog(x, y2, linestyle='-', color='b', label=gasname)
    plt.xlim([xbot, xtop])
    plt.ylim([ybot, ytop])
    plt.legend()
    plt.show()

    
