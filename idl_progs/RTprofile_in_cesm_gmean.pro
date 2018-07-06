pro RTprofile_in_cesm_gmean
;----- 2/19/14 -----
;make input profile for RT standalone

GRAV = 9.80616

no_csky=0
;get file
fname_in ='/scratch/summit/wolfet/archive/t3300_s1400_p22.1392_0.25bar/atm/hist/t3300_s1400_p22.1392_0.25bar.cam.h0.0048-12.nc'                       
;fname_in ='/scratch/summit/wolfet/archive/trappist1e_1barN2_0.4barCO2_aqua/atm/hist/trappist1e_1barN2_0.4barCO2_aqua.cam.h0.0061-01-01-01800.nc'
;fname_in ='/scratch/summit/wolfet/archive/trappist1e_1barN2_1barCO2_aqua/atm/hist/trappist1e_1barN2_1barCO2_aqua.cam.h0.0061-01-01-01800.nc'   
;fname_in ='/scratch/summit/wolfet/archive/rtbuild_test/atm/hist/rtbuild_test.cam.h0.0061-01-01-01800.nc'  
;fname_in ='/scratch/summit/wolfet/archive/rtbuild_test/atm/hist/rtbuild_test_ORIG.cam.h0.0061-01-01-01800.nc'                                                          
;fname_in ='/scratch/summit/wolfet/archive/trappist1e_1barN2_0.4barCO2_aqua/atm/hist/trappist1e_1barN2_0.4barCO2_aqua.cam.h0.0061-01-01-01800.nc'                       
;fname_in ='/scratch/summit/wolfet/archive/trappist1e_1barCO2_aqua/atm/hist/trappist1e_1barCO2_aqua.cam.h0.0061-01-01-01800.nc'                                         
;fname_in ='/projects/wolfet/trappist-1/simulations/trappist1e_1barN2_0.01barCO2_aqua.cam.h0.avg.nc'                                                                    
;fname_in ='/projects/wolfet/trappist-1/simulations/trappist1e_1barN2_2barCO2_aqua.cam.h0.avg.nc'                                                                       
;fname_in ='/projects/wolfet/EXO_ANALYSIS/co2_forcing_co2_forcing_3.2bar_s75.cam.h0.0001-01-01-00000.nc'           
;fname_in = '/projects/wolfet/trappist-1/simulations/trappist1e_1barN2_0.01barCO2_aqua.cam.h0.avg.nc'

ncid=ncdf_open(fname_in, /nowrite)
ncdf_varget,ncid,'lon',lon 
ncdf_varget,ncid,'lat',lat 
ncdf_varget,ncid,'lev',lev 
ncdf_varget,ncid,'ilev',ilev
ncdf_varget,ncid,'hyai',hyai 
ncdf_varget,ncid,'hybi',hybi 
ncdf_varget,ncid,'hyam',hyam 
ncdf_varget,ncid,'hybm',hybm 
ncdf_varget,ncid,'ASDIR',ASDIR
ncdf_varget,ncid,'ASDIF',ASDIF
ncdf_varget,ncid,'ALDIR',ALDIR
ncdf_varget,ncid,'ALDIF',ALDIF
ncdf_varget,ncid,'PS',PS
ncdf_varget,ncid,'P0',P0
ncdf_varget,ncid,'T',T       
ncdf_varget,ncid,'TS',TS     
if (no_csky ne 1) then begin
ncdf_varget,ncid,'QRSC',QRSC    ;shortwave hearting rate, K/s
ncdf_varget,ncid,'QRLC',QRLC    ;longwave heating rate, K/s  
ncdf_varget,ncid,'FDLC',FDLC  ;clearsky downwelling longwave   
ncdf_varget,ncid,'FULC',FULC  ;clearsky upwelling longwave
ncdf_varget,ncid,'FDSC',FDSC  ;clearsky downwelling shortwave
ncdf_varget,ncid,'FUSC',FUSC  ;clearsky upwelling shortwave 
endif
ncdf_varget,ncid,'FDL',FDL  ;clearsky downwelling longwave   
ncdf_varget,ncid,'FUL',FUL  ;clearsky upwelling longwave
ncdf_varget,ncid,'FDS',FDS  ;clearsky downwelling shortwave
ncdf_varget,ncid,'FUS',FUS  ;clearsky upwelling shortwave 
ncdf_varget,ncid,'QRS',QRS  ;shortwave heatng rate
ncdf_varget,ncid,'QRL',QRL  ;longwave heatng rate
ncdf_varget,ncid,'Z3',Z3 ;geopotential height, m
ncdf_varget,ncid,'Q',Q   ;kg/kg, specific humidity
ncdf_varget,ncid,'m_Q_c',m_Q_c  ;kg/m2, mass per layer
ncdf_varget,ncid,'m_CH4_c',m_CH4_c  ;kg/m2, mass per layer
ncdf_varget,ncid,'m_CO2_c',m_CO2_c  ;kg/m2, mass per layer
;ncdf_varget,ncid,'rad_CO2',rad_CO2  ;kg/kg
;ncdf_varget,ncid,'rad_CH4',rad_CH4  ;kg/kg
;ncdf_varget,ncid,'rad_Q',rad_Q  ;kg/kg
ncdf_varget,ncid,'co2vmr',co2vmr
ncdf_varget,ncid,'ch4vmr',ch4vmr
;--- need in cloud ice water path
;--- need in cloud liquid water path
;--- need in cloud fraction
;    ---ext_rei, ext_rel --- need to add these to output data sets if
;                            this is what I want to do.

ncdf_close,ncid

nlon=n_elements(lon)
nlat=n_elements(lat)
nlev=n_elements(lev)
nilev=n_elements(ilev)

if (no_csky eq 1) then begin
FULC = fltarr(nlon,nlat,nlev+1) & FULC(*,*,*) = 0.0
FDLC = fltarr(nlon,nlat,nlev+1) & FDLC(*,*,*) = 0.0
FUSC = fltarr(nlon,nlat,nlev+1) & FUSC(*,*,*) = 0.0
FDSC = fltarr(nlon,nlat,nlev+1) & FDSC(*,*,*) = 0.0
QRLC = fltarr(nlon,nlat,nlev) & QRLC(*,*,*) = 0.0
QRSC = fltarr(nlon,nlat,nlev) & QRSC(*,*,*) = 0.0
endif
;kludge
QRLC = fltarr(nlon,nlat,nlev) & QRLC(*,*,*) = 0.0
QRSC = fltarr(nlon,nlat,nlev) & QRSC(*,*,*) = 0.0

levP=fltarr(nlon,nlat,nlev)
ilevP=fltarr(nlon,nlat,nilev)
levZ=fltarr(nlon,nlat,nlev)
ilevZ=fltarr(nlon,nlat,nilev)

hybrid2pressure,nlon,nlat,nlev,PS,P0,hyam,hybm,hyai,hybi,levP,ilevP
hybrid2height,nlon,nlat,nilev,PS,P0,hyai,hybi,hyam,hybm,T,levZ,ilevz

; calculate mixing ratios, etc
; mass of the atmosphere
;are these supposed to be dry mixing ratios???
mwn2 = 28.
mwar = 40.
mwco2 = 44.
mwch4 = 16.
mwh2o = 18.

cpn2 = 1.039e3
cpco2 = 0.846e3
cpch4 = 2.226e3

n2vmr=1.0-co2vmr-ch4vmr

mwdry = n2vmr*mwn2 + co2vmr*mwco2 + ch4vmr*mwch4 ;+ arvmr*mwar 
cpdry = n2vmr*cpn2 + co2vmr*cpco2 + ch4vmr*cpch4 ;+ arvmr*cpar


mass_atm = fltarr(nlev)
co2mmr_temp = fltarr(nlev)
ch4mmr_temp = fltarr(nlev)
n2mmr_temp = fltarr(nlev)
for i=0, nlev-1 do begin
  ;mass_atm(i) = pdel_temp(i)/grav
  ;co2mmr_temp(i) = m_CO2_c(lonVC,latVC,i)/mass_atm(i)
  ;ch4mmr_temp(i) = m_CH4_c(lonVC,latVC,i)/mass_atm(i)
  ;n2mmr_temp(i) = 1.0 - CO2mmr_temp(i)-CH4mmr_temp(i)
  ;co2mmr_temp(i) = rad_CO2(lonVC,latVC,i)  ;dry
  ;ch4mmr_temp(i) = rad_CH4(lonVC,latVC,i)  ;dry
  ;n2mmr_temp(i) = 1.0 - CO2mmr_temp(i)-CH4mmr_temp(i) ;dry
  co2mmr_temp(i) = co2vmr*mwco2/mwdry
  ch4mmr_temp(i) = ch4vmr*mwch4/mwdry
  n2mmr_temp(i) = 1.0 - CO2mmr_temp(i)-CH4mmr_temp(i) ;dry
endfor

;is this correct?
h2ommr_temp = fltarr(nlev) & h2ommr_temp(*) = Q(lonVC,latVC,*)  ;specific humidity
;h2ommr_temp = fltarr(nlev) & h2ommr_temp(*) = rad_Q(lonVC,latVC,*)

;calculate TINT
tint_temp = fltarr(nilev)
tint_temp(nilev-1) = TS(lonVC,latVC)
tint_temp(0) = T(lonVC,latVC,0) 
for k=1,nlev-1 do begin
  dy = ( alog10(ilevP(lonVC,latVC,k)) - alog10(levP(lonVC,latVC,k)) ) / $
       ( alog10(levP(lonVC,latVC,k-1)) - alog10(levP(lonVC,latVC,k))) 
  tint_temp(k) = T(lonVC,latVC,k) - dy*(T(lonVC,latVC,k) - T(lonVC,latVC,k-1))
  t_print = T(lonVC,latVC,k)
  tint_print = tint_temp(k)
endfor
;tint_temp(0) = T(lonVC,latVC,0) + (T(lonVC,latVC,0) - T(lonVC,latVC,1))*0.5

;-- set final output values --
TS_out = TS(lonVC,latVC)
PS_out = PS(lonVC,latVC)
TMID_out = fltarr(nlev)   & TMID_out(*) = T(lonVC,latVC,*)
TINT_out = fltarr(nilev) & TINT_out(*) = tint_temp(*)
PMID_out = fltarr(nlev)   & PMID_out(*) = levP(lonVC, latVC, *)
PDEL_out = fltarr(nlev)   & PDEL_OUT(*) = pdel_temp(*)
PINT_OUT = fltarr(nilev)  &  PINT_out(*) = ilevP(lonVC, latVC, *)
ZINT_OUT = fltarr(nilev)  & ZINT_out(*) = ilevZ(lonVC, latVC, *)
H2OMMR_out = fltarr(nlev) & H2OMMR_out(*) = h2ommr_temp(*)
CO2MMR_out = fltarr(nlev) & CO2MMR_out(*) = co2mmr_temp(*)
CH4MMR_out = fltarr(nlev) & CH4MMR_out(*) = ch4mmr_temp(*)
O2MMR_out = fltarr(nlev)  & O2MMR_out(*) = 0.0
O3MMR_out = fltarr(nlev)  & O3MMR_out(*) = 0.0
H2MMR_out = fltarr(nlev)  & H2MMR_out(*) = 0.0
N2MMR_out = fltarr(nlev)  & N2MMR_out(*) = n2mmr_temp(*)
;ALDIR_out = ALDIR(lonVC, latVC)
;ALDIF_out = ALDIF(lonVC, latVC)
;ASDIR_out = ASDIR(lonVC, latVC)
;ASDIF_out = ASDIF(lonVC, latVC)
ALDIR_out = 0.0
ALDIF_out = 0.0
ASDIR_out = 0.0
ASDIF_out = 0.0
MW_OUT = mwdry
CP_OUT = cpdry
;; diagnostic outputs for comparison later
FULC_out =fltarr(nilev) & FULC_out(*) = FULC(lonVC, latVC, *)
FDLC_out =fltarr(nilev) & FDLC_out(*) = FDLC(lonVC, latVC, *)
FUSC_out =fltarr(nilev) & FUSC_out(*) = FUSC(lonVC, latVC, *)
FDSC_out =fltarr(nilev) & FDSC_out(*) = FDSC(lonVC, latVC, *)
QRSC_out = fltarr(nlev) & QRSC_out(*)  = QRSC(lonVC, latVC, *)
QRLC_out = fltarr(nlev) & QRLC_out(*)  = QRLC(lonVC, latVC, *)
FUL_out =fltarr(nilev) & FUL_out(*) = FUL(lonVC, latVC, *)
FDL_out =fltarr(nilev) & FDL_out(*) = FDL(lonVC, latVC, *)
FUS_out =fltarr(nilev) & FUS_out(*) = FUS(lonVC, latVC, *)
FDS_out =fltarr(nilev) & FDS_out(*) = FDS(lonVC, latVC, *)
QRS_out = fltarr(nlev) & QRS_out(*)  = QRS(lonVC, latVC, *)  ;*60.*60.*24.
QRL_out = fltarr(nlev) & QRL_out(*)  = QRL(lonVC, latVC, *)  ;*60.*60.*24.

COSZRS_out = 1.0 ;; only matters for a solar computation
                     ;; does not matter for longwave computation



print, tint_out

plot, tmid_out,pmid_out/100., /ylog, yrange=[2000.,0.01], psym=4,symsize=1.0, ystyle=1, xrange=[100,400], xstyle=1
oplot, tmid_out,pmid_out/100., linestyle=1
oplot, tint_out,pint_out/100., psym=1, symsize=1.0
oplot, tint_out,pint_out/100., linestyle=0
stop


; ----- create RTprofle_in.nc --------

filename = 'RTprofile_in.nc'
id = NCDF_CREATE(filename,  /CLOBBER)

dim1 = NCDF_DIMDEF(id, 'pverp', nilev)
dim2 = NCDF_DIMDEF(id, 'pver', nlev)
dim3 = NCDF_DIMDEF(id, 'one', 1)

varid1 = NCDF_VARDEF(id,'ts',dim3,/float)
varid2 = NCDF_VARDEF(id,'ps',dim3,/float)
varid3 = NCDF_VARDEF(id,'tmid',dim2,/float)
varid4 = NCDF_VARDEF(id,'tint',dim1,/float)
varid5 = NCDF_VARDEF(id,'pmid',dim2,/float)
varid6 = NCDF_VARDEF(id,'pdel',dim2,/float)
varid7 = NCDF_VARDEF(id,'pint',dim1,/float)
varid8 = NCDF_VARDEF(id,'zint',dim1,/float)
varid9 = NCDF_VARDEF(id,'h2ommr',dim2,/float)
varid10 = NCDF_VARDEF(id,'co2mmr',dim2,/float)
varid11 = NCDF_VARDEF(id,'ch4mmr',dim2,/float)
varid12 = NCDF_VARDEF(id,'o2mmr',dim2,/float)
varid13 = NCDF_VARDEF(id,'o3mmr',dim2,/float)
varid14 = NCDF_VARDEF(id,'n2mmr',dim2,/float)
varid15 = NCDF_VARDEF(id,'h2mmr',dim2,/float)
varid16 = NCDF_VARDEF(id,'asdir',dim3,/float)
varid17 = NCDF_VARDEF(id,'asdif',dim3,/float)
varid18 = NCDF_VARDEF(id,'aldir',dim3,/float)
varid19 = NCDF_VARDEF(id,'aldif',dim3,/float)
varid20 = NCDF_VARDEF(id,'coszrs',dim3,/float)
varid21 = NCDF_VARDEF(id,'mw',dim3,/float)
varid22 = NCDF_VARDEF(id,'cp',dim3,/float)
varid23 = NCDF_VARDEF(id,'FUL',dim1,/float)
varid24 = NCDF_VARDEF(id,'FDL',dim1,/float)
varid25 = NCDF_VARDEF(id,'FUS',dim1,/float)
varid26 = NCDF_VARDEF(id,'FDS',dim1,/float)
varid27 = NCDF_VARDEF(id,'QRS',dim2,/float)
varid28 = NCDF_VARDEF(id,'QRL',dim2,/float)
varid29 = NCDF_VARDEF(id,'FULC',dim1,/float)
varid30 = NCDF_VARDEF(id,'FDLC',dim1,/float)
varid31 = NCDF_VARDEF(id,'FUSC',dim1,/float)
varid32 = NCDF_VARDEF(id,'FDSC',dim1,/float)
varid33 = NCDF_VARDEF(id,'QRSC',dim2,/float)
varid34 = NCDF_VARDEF(id,'QRLC',dim2,/float)


NCDF_ATTPUT, id, varid1, "title", "Surface temperature"           & NCDF_ATTPUT, id, varid1, "units", "K"
NCDF_ATTPUT, id, varid2, "title", "Surface pressure"              & NCDF_ATTPUT, id, varid2, "units", "Pa"
NCDF_ATTPUT, id, varid3, "title", "Mid-layer temperatures"        & NCDF_ATTPUT, id, varid3, "units", "K"
NCDF_ATTPUT, id, varid4, "title", "Interface temperatures"        & NCDF_ATTPUT, id, varid4, "units", "K"
NCDF_ATTPUT, id, varid5, "title", "Mid-layer pressures"           & NCDF_ATTPUT, id, varid5, "units", "Pa"
NCDF_ATTPUT, id, varid6, "title", "Pressures difference"          & NCDF_ATTPUT, id, varid6, "units", "Pa"
NCDF_ATTPUT, id, varid7, "title", "Interface pressures"           & NCDF_ATTPUT, id, varid7, "units", "Pa"
NCDF_ATTPUT, id, varid8, "title", "Interface heights"             & NCDF_ATTPUT, id, varid8, "units", "m"
NCDF_ATTPUT, id, varid9, "title", "H2O MMR"                       & NCDF_ATTPUT, id, varid9, "units", "mmr"
NCDF_ATTPUT, id, varid10, "title", "CO2 MMR"                      & NCDF_ATTPUT, id, varid10, "units", "mmr"
NCDF_ATTPUT, id, varid11, "title", "CH4 MMR"                      & NCDF_ATTPUT, id, varid11, "units", "mmr"
NCDF_ATTPUT, id, varid12, "title", "O2 MMR"                       & NCDF_ATTPUT, id, varid12, "units", "mmr"
NCDF_ATTPUT, id, varid13, "title", "O3 MMR"                       & NCDF_ATTPUT, id, varid13, "units", "mmr"
NCDF_ATTPUT, id, varid14, "title", "N2 MMR"                       & NCDF_ATTPUT, id, varid14, "units", "mmr"
NCDF_ATTPUT, id, varid15, "title", "H2 MMR"                       & NCDF_ATTPUT, id, varid15, "units", "mmr"
NCDF_ATTPUT, id, varid16, "title", "albedo shortwave direct"      & NCDF_ATTPUT, id, varid16, "units", "albedo"
NCDF_ATTPUT, id, varid17, "title", "albedo shortwave diffuse"     & NCDF_ATTPUT, id, varid17, "units", "albedo"
NCDF_ATTPUT, id, varid18, "title", "albedo longwave direct"       & NCDF_ATTPUT, id, varid18, "units", "albedo"
NCDF_ATTPUT, id, varid19, "title", "albedo longwave diffuse"      & NCDF_ATTPUT, id, varid19, "units", "albedo"
NCDF_ATTPUT, id, varid20, "title", "cosine of zenith angle"       & NCDF_ATTPUT, id, varid20, "units", "none"
NCDF_ATTPUT, id, varid21, "title", "molecular weight of dry air"  & NCDF_ATTPUT, id, varid21, "units", "AMU"
NCDF_ATTPUT, id, varid22, "title", "specific heat of dry air"     & NCDF_ATTPUT, id, varid22, "units", "J/kg K"
NCDF_ATTPUT, id, varid23, "title", "diag :: FUL"     & NCDF_ATTPUT, id, varid23, "units", "W m-2"
NCDF_ATTPUT, id, varid24, "title", "diag :: FDL"     & NCDF_ATTPUT, id, varid24, "units", "W m-2"
NCDF_ATTPUT, id, varid25, "title", "diag :: FUS"     & NCDF_ATTPUT, id, varid25, "units", "W m-2"
NCDF_ATTPUT, id, varid26, "title", "diag :: FDS"     & NCDF_ATTPUT, id, varid26, "units", "W m-2"
NCDF_ATTPUT, id, varid27, "title", "diag :: QRS"      & NCDF_ATTPUT, id, varid27, "units", "K/day"
NCDF_ATTPUT, id, varid28, "title", "diag :: QRL"      & NCDF_ATTPUT, id, varid28, "units", "K/day"
NCDF_ATTPUT, id, varid29, "title", "diag :: FULC"     & NCDF_ATTPUT, id, varid29, "units", "W m-2"
NCDF_ATTPUT, id, varid30, "title", "diag :: FDLC"     & NCDF_ATTPUT, id, varid30, "units", "W m-2"
NCDF_ATTPUT, id, varid31, "title", "diag :: FUSC"     & NCDF_ATTPUT, id, varid31, "units", "W m-2"
NCDF_ATTPUT, id, varid32, "title", "diag :: FDSC"     & NCDF_ATTPUT, id, varid32, "units", "W m-2"
NCDF_ATTPUT, id, varid33, "title", "diag :: QRSC"      & NCDF_ATTPUT, id, varid33, "units", "K/day"
NCDF_ATTPUT, id, varid34, "title", "diag :: QRLC"      & NCDF_ATTPUT, id, varid34, "units", "K/day"

NCDF_CONTROL, id, /ENDEF

NCDF_VARPUT, id, varid1,  TS_OUT
NCDF_VARPUT, id, varid2,  PS_OUT
NCDF_VARPUT, id, varid3,  TMID_OUT
NCDF_VARPUT, id, varid4,  TINT_OUT
NCDF_VARPUT, id, varid5,  PMID_OUT
NCDF_VARPUT, id, varid6,  PDEL_OUT
NCDF_VARPUT, id, varid7,  PINT_OUT
NCDF_VARPUT, id, varid8,  ZINT_OUT
NCDF_VARPUT, id, varid9,  H2OMMR_OUT
NCDF_VARPUT, id, varid10, CO2MMR_OUT
NCDF_VARPUT, id, varid11, CH4MMR_OUT
NCDF_VARPUT, id, varid12, O2MMR_OUT
NCDF_VARPUT, id, varid13, O3MMR_OUT
NCDF_VARPUT, id, varid14, N2MMR_OUT
NCDF_VARPUT, id, varid15, H2MMR_OUT
NCDF_VARPUT, id, varid16, ASDIR_OUT
NCDF_VARPUT, id, varid17, ASDIF_OUT
NCDF_VARPUT, id, varid18, ALDIR_OUT
NCDF_VARPUT, id, varid19, ALDIF_OUT
NCDF_VARPUT, id, varid20, COSZRS_OUT
NCDF_VARPUT, id, varid21, MW_OUT
NCDF_VARPUT, id, varid22, CP_OUT
NCDF_VARPUT, id, varid23, FUL_OUT
NCDF_VARPUT, id, varid24, FDL_OUT
NCDF_VARPUT, id, varid25, FUS_OUT
NCDF_VARPUT, id, varid26, FDS_OUT
NCDF_VARPUT, id, varid27, QRS_OUT
NCDF_VARPUT, id, varid28, QRL_OUT
NCDF_VARPUT, id, varid29, FULC_OUT
NCDF_VARPUT, id, varid30, FDLC_OUT
NCDF_VARPUT, id, varid31, FUSC_OUT
NCDF_VARPUT, id, varid32, FDSC_OUT
NCDF_VARPUT, id, varid33, QRSC_OUT
NCDF_VARPUT, id, varid34, QRLC_OUT

NCDF_CLOSE, id

FINISH:

end
