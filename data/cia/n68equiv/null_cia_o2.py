import numpy as np
import netCDF4 as nc
import pathlib
import shutil

files = ['O2-O2_cia_n68.nc','O2-N2_cia_n68.nc','O2-CO2_cia_n68.nc']

for ifile in files:
    file_in = pathlib.Path(ifile)
    file_out = pathlib.Path('null_'+ifile)

    shutil.copy(file_in,file_out)

    f = nc.Dataset(str(file_out),'r+')
    f['sigma'][...] = 0.0
    f.close()
