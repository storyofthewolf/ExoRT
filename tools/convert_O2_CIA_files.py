import numpy as np
import matplotlib.pyplot as plt
import scipy.interpolate as sint
from scipy.integrate import simpson
import netCDF4 as nc

#setup wavenumber x temperature grid used in netcdf files
#T_ext = np.array([100, 125, 150, 175, 200, 225, 250, 275, 300, 325, 350, 375, 400, 425, 450, 475, 500])
bin_edges =  np.loadtxt('../HELIOS-K/HELIOS-K_n68equiv/bins68.dat',unpack=True)
nu_per_bin = 10000
nu_ref = np.zeros((len(bin_edges)-1)*nu_per_bin)
for ibin in np.arange(len(bin_edges)-1):
    spacing = (bin_edges[ibin+1] - bin_edges[ibin])/nu_per_bin
    nu_ref[ibin*nu_per_bin:(ibin+1)*nu_per_bin] = np.linspace(bin_edges[ibin],bin_edges[ibin+1]-spacing,nu_per_bin)

loschmidt = 2.6867774e19  #loschmidt number in molecules cm^-3
#hitran cia data should have units of cm^5 molecule^-2

def band_average(bins,nu,T,k):
    k_band = np.zeros([len(bins)-1,len(T)])
    for i in np.arange(len(bins)-1):
        for j in np.arange(len(T)):
            waven_range = np.logical_and(nu>=bins[i],nu<=bins[i+1])
            #k_band[i,j] = np.median(k[waven_range,j])
            #k_band[i,j] = 10**(np.mean(np.log10(k[waven_range,j])))

            z = k[waven_range,j]
            z_avg = simpson(z,x=nu[waven_range])/(bins[i+1]-bins[i])
            k_band[i,j] = z_avg

            # z = np.log10(k[waven_range,j])
            # dx = np.zeros_like(z) + (bins[i+1]-bins[i])/(len(z)-1)
            # #import pdb; pdb.set_trace()
            # z_avg = np.nansum(z*dx)/(bins[i+1]-bins[i])
            # k_band[i,j] = 10**z_avg

            z = k[waven_range,j]
            x = nu[waven_range]
            x_left = np.hstack([x[0]-0.5*(x[1]-x[0]),0.5*(x[1:]+x[:-1])])
            x_right = np.hstack([0.5*(x[1:]+x[:-1]),x[-1]+0.5*(x[-1]-x[-2])])
            dx = x_right - x_left

            if (z==0).all():
                k_band[i,j] = 0.0
            elif (z!=0).all():
                k_band[i,j] = 10**(np.sum(np.log10(z)*dx) / np.sum(dx))
                # k_band[i,j] = 10**np.mean(np.log10(z))
            else:
                z_nz = z[z!=0]  #nonzero part of z
                dx_nz = dx[z!=0]
                z_avg_nz = 10**(np.sum(np.log10(z_nz)*dx_nz) / np.sum(dx_nz))

                # dx_zr = dx[z==0]
                k_band[i,j] = np.sum(dx_nz)/(bins[i+1]-bins[i])*z_avg_nz
                # import pdb; pdb.set_trace()


    return k_band

def read_file(filename, species, nsets, ax):
    #from readme file, starting index of each wavenumber group (except the first)
    count = 0
    wavestart = np.array([])
    waveend = np.array([])
    Temps = np.array([])
    ns = []
    ks = []
    chunk = 0
    iset = 0
    group = []
    for line in open(filename,'r'):
        if species in line.split()[0]:
            nlines = int(line.split()[3])
            wavestart = np.concatenate([wavestart,np.array([float(line.split()[1])])])
            waveend = np.concatenate([waveend,np.array([float(line.split()[2])])])
            Temps = np.concatenate([Temps,np.array([float(line.split()[4])])])

            n = np.zeros(nlines)
            k = np.zeros(nlines)
        else:
            n[count] = float(line.split()[0])
            k[count] = float(line.split()[1])
            count += 1

            if count == nlines:
                ax.semilogy(n,k)
                # ax.plot(n,k)
                ns.append(n)
                ks.append(k)
                count = 0
                #if chunk < nsets[iset]:
                if chunk == nsets[iset]:
                    iset += 1
                group.append(iset)

                #if nsets[iset] >= 24:

                #   axes[0].plot(np.arange(len(n)),n,linestyle='-',marker='.')
                chunk += 1

    return Temps, ns, ks

def create_nc_file(filename,bin_edges,temps,k_band):
    fout = nc.Dataset(filename,'w')

    #set dimensions
    wbins = fout.createDimension("wbins",len(bin_edges)-1)
    tbins = fout.createDimension("tbins",len(temps))

    #setup variables
    nu_low = fout.createVariable("nu_low","f8",("wbins",))
    nu_low.title = "Wavenumber grid low"
    nu_low.units = "cm-1"
    nu_low[:] = bin_edges[:-1]

    nu_high = fout.createVariable("nu_high","f8",("wbins",))
    nu_high.title = "Wavenumber grid high"
    nu_high.units = "cm-1"
    nu_high[:] = bin_edges[1:]

    temp = fout.createVariable("temp","f8",("tbins",))
    temp.title = "temperature grid"
    temp.units = "K"
    temp[:] = temps

    sigma = fout.createVariable("sigma","f8",("tbins","wbins"))
    sigma.title = "band averaged CIA absorption"
    sigma.units = "cm-1 amagat-2"
    sigma[:,:] = k_band.T

    fout.close()

#---------O2-O2-----------------------------------------------------------------
#start with o2-o2 from hitran
file_in = 'O2-O2_2018b.cia'
nsets = [15, 16, 17, 18, 19, 20, 24, 29]  #from HITRAN documentation

fig, axes = plt.subplots(ncols=1,nrows=3,figsize=(7.5,9))

Temps, ns, ks = read_file(file_in, 'O2-O2', nsets, axes[0])

# for ins in np.arange(len(ns)):
#     if (ks[ins]<0).any():
#         print(f'Negative values present in O2-O2 chunk {ins}')

T_ext = Temps[:nsets[0]]
#destination k values, grid over nu_ref, T_ext
k_final = np.zeros([len(nu_ref),len(T_ext)])

axes[0].set(xlabel='Wavenumber (cm$^{-1}$)',ylabel='CIA cross section',xlim=(0,30000),ylim=(1e-50,1e-43))


############################################### Group 1
#destination k values, grid over nu_ref, Temps for this group
#k_tmp = np.zeros([len(nu_ref),nsets[0]])
#id valid wavenumber range (same for all sets in this group)
waven_range = np.logical_and(nu_ref>=ns[0][0],nu_ref<=ns[0][-1])

for i in np.arange(1,nsets[0]):
    #first verify wavenumber grid is the same in each set
    # print(f"Group 1, set {i}, waven grid matches set 0: {(ns[0]==ns[i]).all()}")

    k_final[waven_range,i] = sint.interp1d(ns[i],ks[i])(nu_ref[waven_range])

#now interpolate onto T_ext grid
# for iT in np.arange(len(T_ext)):
#     #extrapolate as constant below and above measured values
#     if T_ext[iT] < np.min(Temps[:nsets[0]]):
#         k_final[waven_range,iT] = k_tmp[waven_range,0]
#
#     elif T_ext[iT] > np.max(Temps[:nsets[0]]):
#         k_final[waven_range,iT] = k_tmp[waven_range,-1]
#
#     #interpolate within the measured range
#     else:
#         k_final[waven_range,iT] = sint.interp1d(Temps[:nsets[0]],k_tmp[waven_range,:],axis=1)(T_ext[iT])


############################################### Group 2 - 6
#the next five blocks have only one temperature
for i in np.arange(nsets[0],nsets[5]):
    waven_range = np.logical_and(nu_ref>=ns[i][0],nu_ref<=ns[i][-1])

    k_tmp = sint.interp1d(ns[i],ks[i])(nu_ref[waven_range])
    k_final[waven_range,:] = k_tmp[:,None]


############################################### Group 7
#destination k values, grid over nu_ref, Temps for this group
k_tmp = np.zeros([len(nu_ref),nsets[6]-nsets[5]])
#id valid wavenumber range (one set has slightly diff range in this group)
waven_range = np.logical_and(nu_ref>=ns[20][0],nu_ref<=ns[20][-1])

for i in np.arange(nsets[5],nsets[6]):
    #since one set has slightly different range, we will fill with zeros
    k_tmp[waven_range,i-nsets[5]] = sint.interp1d(ns[i],ks[i],fill_value=0,bounds_error=False)(nu_ref[waven_range])

#now interpolate onto T_ext grid
for iT in np.arange(len(T_ext)):
    #extrapolate as constant below and above measured values
    if T_ext[iT] < np.min(Temps[nsets[5]:nsets[6]]):
        k_final[waven_range,iT] = k_tmp[waven_range,0]

    elif T_ext[iT] > np.max(Temps[nsets[5]:nsets[6]]):
        k_final[waven_range,iT] = k_tmp[waven_range,-1]

    #interpolate within the measured range
    else:
        k_final[waven_range,iT] = sint.interp1d(Temps[nsets[5]:nsets[6]],k_tmp[waven_range,:],axis=1)(T_ext[iT])


############################################### Group 8
#destination k values, grid over nu_ref, Temps for this group
k_tmp = np.zeros([len(nu_ref),nsets[7]-nsets[6]])
#id valid wavenumber range (PITA cuz these all have slightly different ranges)
min = np.min(ns[nsets[6]])
max = np.max(ns[nsets[6]])
for i in np.arange(nsets[6]+1,nsets[7]):
    if np.min(ns[i]) < min:
        min = np.min(ns[i])
    if np.max(ns[i]) > max:
        max = np.max(ns[i])
waven_range = np.logical_and(nu_ref>=min,nu_ref<=max)

for i in np.arange(nsets[6],nsets[7]):
    #since one set has slightly different range, we will fill with zeros
    k_tmp[waven_range,i-nsets[7]] = sint.interp1d(ns[i],ks[i],fill_value=0,bounds_error=False)(nu_ref[waven_range])

#now interpolate onto T_ext grid
for iT in np.arange(len(T_ext)):
    #extrapolate as constant below and above measured values
    if T_ext[iT] < np.min(Temps[nsets[6]:nsets[7]]):
        k_final[waven_range,iT] = k_tmp[waven_range,0]

    elif T_ext[iT] > np.max(Temps[nsets[6]:nsets[7]]):
        k_final[waven_range,iT] = k_tmp[waven_range,-1]

    #interpolate within the measured range
    else:
        k_final[waven_range,iT] = sint.interp1d(Temps[nsets[6]:nsets[7]],k_tmp[waven_range,:],axis=1)(T_ext[iT])

#Looking at Thalman & Volkamer 2013, there are gaps in the data for group 8
#These have been interpolated over in the HITRAN data but most are small.
#2 gaps are significant: 386.1 - 433.5 nm and 487.6 - 511.5 nm.
#I am nulling these here because there is a lot of spurious absortance in the gaps
#Who knows if it is the right thing to do
gap1_edges = 1.0/(np.array([433.5,386.1])/1e7)
gap1 = np.logical_and(nu_ref>np.min(gap1_edges),nu_ref<np.max(gap1_edges))
k_final[gap1,:] = 0.0

gap2_edges = 1.0/(np.array([511.5,487.6])/1e7)
gap2 = np.logical_and(nu_ref>np.min(gap2_edges),nu_ref<np.max(gap2_edges))
k_final[gap2,:] = 0.0


for iT in np.arange(len(T_ext)):
    axes[1].semilogy(nu_ref,k_final[:,iT])

axes[1].set(xlabel='Wavenumber (cm$^{-1}$)',ylabel='CIA cross section',xlim=(0,30000),ylim=(1e-50,1e-43))


# set weird negative values to 0 -- should check papers to see if this makes sense!!
k_final[k_final<0] = 0.0

k_final *= loschmidt**2  #convert to units that are used in ExoRT
#k_final[k_final == 0] = 1e-20

#do band averaging
k_band = band_average(bin_edges,nu_ref,T_ext,k_final)

for iT in np.arange(len(T_ext)):
    axes[2].semilogy(nu_ref,k_final[:,iT])
    for ibin in np.arange(len(bin_edges)-1):
        axes[2].plot(np.array([bin_edges[ibin],bin_edges[ibin+1]]),np.array([k_band[ibin,iT],k_band[ibin,iT]]),'k-',zorder=100)

axes[2].set(xlabel='Wavenumber (cm$^{-1}$)',ylabel='CIA cross section (scaled)',xlim=(0,30000))
ylims = axes[2].get_ylim()
axes[2].vlines(bin_edges,np.min(ylims),np.max(ylims),colors='0.7',zorder=-100)

plt.tight_layout()
plt.savefig('O2-O2_cia_plots.png',dpi=300)
# plt.show()
plt.close()

#create the netcdf file
create_nc_file('O2-O2_cia_68bin.nc',bin_edges,T_ext,k_band)

#---------O2-O2 end-------------------------------------------------------------



#---------O2-N2-----------------------------------------------------------------
file_in = 'O2-N2_2021.cia'
nsets = [7, 8, 13, 18, 19]  #from HITRAN documentation

fig, axes = plt.subplots(ncols=1,nrows=3,figsize=(7.5,9))

Temps, ns, ks = read_file(file_in, 'O2', nsets, axes[0])

# for ins in np.arange(len(ns)):
#     if (ks[ins]<0).any():
#         print(f'Negative values present in O2-N2 chunk {ins}')

T_ext = Temps[:nsets[0]]  #doesn't cover full temperature range but is close enough

#destination k values, grid over nu_ref, T_ext
k_final = np.zeros([len(nu_ref),len(T_ext)])

axes[0].set(xlabel='Wavenumber (cm$^{-1}$)',ylabel='CIA cross section',xlim=(0,14000),ylim=(1e-50,1e-43))

############################################### Group 1
#destination k values, grid over nu_ref, Temps for this group
#k_tmp = np.zeros([len(nu_ref),nsets[0]])
#id valid wavenumber range (not the same for all sets in this group!)
waven_range = np.logical_and(nu_ref>=ns[0][0],nu_ref<=ns[0][-1])

for i in np.arange(1,nsets[0]):
    k_final[waven_range,i] = sint.interp1d(ns[i],ks[i],fill_value=0,bounds_error=False)(nu_ref[waven_range])


############################################### Group 2
#only one temperature
waven_range = np.logical_and(nu_ref>=ns[nsets[0]][0],nu_ref<=ns[nsets[0]][-1])

k_tmp = sint.interp1d(ns[nsets[0]],ks[nsets[0]])(nu_ref[waven_range])
k_final[waven_range,:] = k_tmp[:,None]

############################################### Group 3
#this one and the next are overlapping
#Group 3, low temp, smaller waven range
#Group 4, higher temp, extends farther in waven in both directions
#we'll split them at T = 280 K  (index 5/6 in T_ext)
#At least the waven ranges are consistent within each group!
#destination k values, grid over nu_ref, Temps for this group
k_tmp = np.zeros([len(nu_ref),nsets[2]-nsets[1]])
#id valid wavenumber range (one set has slightly diff range in this group)
waven_range = np.logical_and(nu_ref>=ns[nsets[1]][0],nu_ref<=ns[nsets[1]][-1])

for i in np.arange(nsets[1],nsets[2]):
    #since one set has slightly different range, we will fill with zeros
    k_tmp[waven_range,i-nsets[1]] = sint.interp1d(ns[i],ks[i])(nu_ref[waven_range])

#now interpolate onto T_ext grid
for iT in np.arange(5):
    #extrapolate as constant below and above measured values
    if T_ext[iT] < np.min(Temps[nsets[1]:nsets[2]]):
        k_final[waven_range,iT] = k_tmp[waven_range,0]

    elif T_ext[iT] > np.max(Temps[nsets[1]:nsets[2]]):
        k_final[waven_range,iT] = k_tmp[waven_range,-1]

    #interpolate within the measured range
    else:
        k_final[waven_range,iT] = sint.interp1d(Temps[nsets[1]:nsets[2]],k_tmp[waven_range,:],axis=1)(T_ext[iT])

############################################### Group 4
#this one and the previous are overlapping
#see complaining above
#destination k values, grid over nu_ref, Temps for this group
k_tmp = np.zeros([len(nu_ref),nsets[3]-nsets[2]])
#id valid wavenumber range (one set has slightly diff range in this group)
waven_range = np.logical_and(nu_ref>=ns[nsets[2]][0],nu_ref<=ns[nsets[2]][-1])

for i in np.arange(nsets[2],nsets[3]):
    #since one set has slightly different range, we will fill with zeros
    k_tmp[waven_range,i-nsets[2]] = sint.interp1d(ns[i],ks[i])(nu_ref[waven_range])

#now interpolate onto T_ext grid
for iT in np.arange(5,len(T_ext)):
    #extrapolate as constant below and above measured values
    if T_ext[iT] < np.min(Temps[nsets[2]:nsets[3]]):
        k_final[waven_range,iT] = k_tmp[waven_range,0]

    elif T_ext[iT] > np.max(Temps[nsets[2]:nsets[3]]):
        k_final[waven_range,iT] = k_tmp[waven_range,-1]

    #interpolate within the measured range
    else:
        k_final[waven_range,iT] = sint.interp1d(Temps[nsets[2]:nsets[3]],k_tmp[waven_range,:],axis=1)(T_ext[iT])

############################################### Group 5
#only one temperature
waven_range = np.logical_and(nu_ref>=ns[nsets[3]][0],nu_ref<=ns[nsets[3]][-1])

k_tmp = sint.interp1d(ns[nsets[3]],ks[nsets[3]])(nu_ref[waven_range])
k_final[waven_range,:] = k_tmp[:,None]


for iT in np.arange(len(T_ext)):
    axes[1].semilogy(nu_ref,k_final[:,iT])

axes[1].set(xlabel='Wavenumber (cm$^{-1}$)',ylabel='CIA cross section',xlim=(0,14000),ylim=(1e-50,1e-43))


# set weird negative values to 0 -- should check papers to see if this makes sense!!
k_final[k_final<0] = 0.0

k_final *= loschmidt**2  #convert to units that are used in ExoRT
#k_final[k_final == 0] = 1e-20

#do band averaging
k_band = band_average(bin_edges,nu_ref,T_ext,k_final)

for iT in np.arange(len(T_ext)):
    axes[2].semilogy(nu_ref,k_final[:,iT])
    for ibin in np.arange(len(bin_edges)-1):
        axes[2].plot(np.array([bin_edges[ibin],bin_edges[ibin+1]]),np.array([k_band[ibin,iT],k_band[ibin,iT]]),'k-',zorder=100)

axes[2].set(xlabel='Wavenumber (cm$^{-1}$)',ylabel='CIA cross section (scaled)',xlim=(0,14000))
ylims = axes[2].get_ylim()
axes[2].vlines(bin_edges,np.min(ylims),np.max(ylims),colors='0.7',zorder=-100)


plt.tight_layout()
plt.savefig('O2-N2_cia_plots.png',dpi=300)
plt.close()

#create the netcdf file
create_nc_file('O2-N2_cia_68bin.nc',bin_edges,T_ext,k_band)

#---------O2-N2 end-------------------------------------------------------------



#---------O2-CO2----------------------------------------------------------------
file_in = 'O2-CO2_2011.cia'
nsets = [1]  #from HITRAN documentation

fig, axes = plt.subplots(ncols=1,nrows=3,figsize=(7.5,9))

Temps, ns, ks = read_file(file_in, 'O2', nsets, axes[0])
axes[0].set(xlabel='Wavenumber (cm$^{-1}$)',ylabel='CIA cross section',xlim=(12000,14000),ylim=(1e-50,1e-43))

T_ext = np.array([Temps[0]])
k_final = np.zeros([len(nu_ref),len(T_ext)])

############################################### Group 2
#only one temperature
waven_range = np.logical_and(nu_ref>=ns[0][0],nu_ref<=ns[0][-1])

k_tmp = sint.interp1d(ns[0],ks[0])(nu_ref[waven_range])
k_final[waven_range,:] = k_tmp[:,None]

for iT in np.arange(len(T_ext)):
    axes[1].semilogy(nu_ref,k_final[:,iT])

axes[1].set(xlabel='Wavenumber (cm$^{-1}$)',ylabel='CIA cross section',xlim=(12000,14000),ylim=(1e-50,1e-43))

# set weird negative values to 0 -- should check papers to see if this makes sense!!
k_final[k_final<0] = 0.0

k_final *= loschmidt**2  #convert to units that are used in ExoRT
#k_final[k_final == 0] = 1e-20

#do band averaging
k_band = band_average(bin_edges,nu_ref,T_ext,k_final)

for iT in np.arange(len(T_ext)):
    axes[2].semilogy(nu_ref,k_final[:,iT])
    for ibin in np.arange(len(bin_edges)-1):
        axes[2].plot(np.array([bin_edges[ibin],bin_edges[ibin+1]]),np.array([k_band[ibin,iT],k_band[ibin,iT]]),'k-',zorder=100)

axes[2].set(xlabel='Wavenumber (cm$^{-1}$)',ylabel='CIA cross section (scaled)',xlim=(0,14000))
ylims = axes[2].get_ylim()
axes[2].vlines(bin_edges,np.min(ylims),np.max(ylims),colors='0.7',zorder=-100)


plt.tight_layout()
plt.savefig('O2-CO2_cia_plots.png',dpi=300)
plt.close()

#create the netcdf file
create_nc_file('O2-CO2_cia_68bin.nc',bin_edges,T_ext,k_band)
