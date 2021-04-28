pro makecolumn
;------------------------------------------------
;Author: Wolf, E.T.
;Revision history
;6/30/2020 make input profile for RT standalone

;
;DESCRIPTION:  
;Manually create a 1D column to feed into exoRT
;RTprofile_in.nc
;
;------------------------------------------------

;set gravity accordingly for your planet!
;Earth
;GRAV = 9.80616
;R = 287.04

; Mars
GRAV = 3.7
R = 188.92

do_plot = 1
do_override = 1


;;================ 2 bar CO2 dry column ========================
pint_smart_2bar_t250      =  [  5.65200e-05,  7.21300e-05,  9.18000e-05,  0.000116500,  0.000147500, $
                                0.000186200,  0.000234400,  0.000294400,  0.000368700,  0.000460500,  0.000573500, $
                                0.000712400,  0.000882600,  0.00109100,   0.00134400,   0.00165100,   0.00202400, $
                                0.00247400,   0.00301600,   0.00366600,   0.00444500,   0.00537500,   0.00648200, $
                                0.00779600,   0.00935000,   0.0111800,    0.0133400,    0.0158800,    0.0188400, $
                                0.0223000,    0.0263100,    0.0309700,    0.0363600,    0.0425600,    0.0497000, $
                                0.0578700,    0.0672000,    0.0778300,    0.0899000,    0.103600,     0.119000, $
                                0.136300,     0.155800,     0.177500,     0.201700,     0.228700,     0.258500, $
                                0.291400,     0.327700,     0.367400,     0.410900,     0.458300,     0.509700, $
                                0.565500,     0.625600,     0.690200,     0.759500,     0.833500,     0.912200, $
                                0.995700,     1.08400,      1.17700,      1.27400,      1.37600,      1.48200, $
                                1.59100,      1.70400,      1.82100,      1.94000,      2.00100 ] *1.0e5 ; Pascals  

tint_smart_2bar_t250 =       [  167.000,      167.000,      167.000,      167.000,      167.000, $
                                167.000,      167.000,      167.000,      167.000,      167.000,      167.000, $
                                167.000,      167.000,      167.000,      167.000,      167.000,      167.000, $
                                167.000,      167.000,      167.000,      167.000,      167.000,      167.000, $
                                167.000,      167.000,      167.000,      167.000,      167.000,      167.000, $
                                167.000,      167.000,      167.000,      167.000,      167.000,      167.000, $
                                167.000,      167.000,      167.100,      168.400,      169.700,      171.000, $
                                172.300,      173.500,      174.800,      176.100,      177.300,      178.600, $
                                179.800,      181.000,      182.200,      183.400,      184.600,      185.800, $
                                187.000,      188.100,      189.200,      190.300,      194.300,      199.400, $
                                204.500,      209.600,      214.600,      219.600,      224.500,      229.300, $
                                234.100,      238.800,      243.300,      247.800,      250.000 ] 

pmid_smart_2bar_t250 =       [  6.43250e-05,  8.19650e-05,  0.000104150,  0.000132000, $
                                0.000166850,  0.000210300,  0.000264400,  0.000331550,  0.000414600,  0.000517000, $
                                0.000642950,  0.000797500,  0.000986800,  0.00121750,   0.00149750,   0.00183750, $
                                0.00224900,   0.00274500,   0.00334100,   0.00405550,   0.00491000,   0.00592850, $
                                0.00713900,   0.00857300,   0.0102650,    0.0122600,    0.0146100,    0.0173600, $
                                0.0205700,    0.0243050,    0.0286400,    0.0336650,    0.0394600,    0.0461300, $
                                0.0537850,    0.0625350,    0.0725150,    0.0838650,    0.0967500,    0.111300, $
                                0.127650,     0.146050,     0.166650,     0.189600,     0.215200,     0.243600, $
                                0.274950,     0.309550,     0.347550,     0.389150,     0.434600,     0.484000, $
                                0.537600,     0.595550,     0.657900,     0.724850,     0.796500,     0.872850, $
                                0.953950,     1.03980,      1.13050,      1.22550,      1.32500,      1.42900, $
                                1.53650,      1.64750,      1.76250,      1.88050,      1.97050 ] *1.0e5 ;Pascal s  

tmid_smart_2bar_t250       = [  167.000,      167.000,      167.000,      167.000, $
                                167.000,      167.000,      167.000,      167.000,      167.000,      167.000, $
                                167.000,      167.000,      167.000,      167.000,      167.000,      167.000, $
                                167.000,      167.000,      167.000,      167.000,      167.000,      167.000, $
                                167.000,      167.000,      167.000,      167.000,      167.000,      167.000, $
                                167.000,      167.000,      167.000,      167.000,      167.000,      167.000, $
                                167.000,      167.000,      167.050,      167.750,      169.050,      170.350, $
                                171.650,      172.900,      174.150,      175.450,      176.700,      177.950, $
                                179.200,      180.400,      181.600,      182.800,      184.000,      185.200, $
                                186.400,      187.550,      188.650,      189.750,      192.300,      196.850, $
                                201.950,      207.050,      212.100,      217.100,      222.050,      226.900, $
                                231.700,      236.450,      241.050,      245.550,      248.900 ] 

nlev = n_elements(tmid_smart_2bar_t250)
nilev = n_elements(tint_smart_2bar_t250)

pint_in = pint_smart_2bar_t250
pmid_in = pmid_smart_2bar_t250
tint_in = tint_smart_2bar_t250
tmid_in = tmid_smart_2bar_t250
ts_in   = tint_smart_2bar_t250(nilev-1)
ps_in   = pint_smart_2bar_t250(nilev-1)
print, "TS: ", ts_in
print, "PS: ", ps_in
print, "nlevs: ", nlev

pdel_in = fltarr(nlev)
pdel_in(*) = 0.0
; calculate pdel
for j=0, nlev-1 do begin
  pdel_in(j) = pint_in(j+1)-pint_in(j)
endfor

zint_in = fltarr(nilev)
zint_in(*) = 0.
; calculate zint
for z=nilev-1, 1,-1 do begin
  delta_Z=R*tmid_in(z-1)/GRAV*alog(pint_in(z)/pint_in(z-1)) 
  zint_in(z-1) = zint_in(z) + delta_Z
endfor

;set via volume mixing ratios
co2vmr_o = 1.0
ch4vmr_o = 0.0
h2vmr_o = 0.0
n2vmr_o = 0.0 
o2vmr_o = 0.0
o3vmr_o = 0.0
h2ovmr_o = 0.0

; calculate mixing ratios, etc
; mass of the atmosphere
;are these supposed to be dry mixing ratios???
mwn2 = 28.
mwar = 40.
mwco2 = 44.01
mwch4 = 16.
mwh2o = 18.
mwh2 = 2.

cpn2 = 1.039e3
;cpco2 = 0.846e3
cpco2 = 0.751e3
cpch4 = 2.226e3
cph2 = 14.32e3


n2vmr_temp  = fltarr(nlev)  & n2vmr_temp(*) = n2vmr_o
co2vmr_temp = fltarr(nlev)  & co2vmr_temp(*) = co2vmr_o
ch4vmr_temp = fltarr(nlev)  & ch4vmr_temp(*) = ch4vmr_o
n2vmr_temp  = fltarr(nlev)   & n2vmr_temp(*) = n2vmr_o
h2vmr_temp  = fltarr(nlev)   & h2vmr_temp(*) = h2vmr_o
o2vmr_temp  = fltarr(nlev)   & o2vmr_temp(*) = o2vmr_o
o3vmr_temp  = fltarr(nlev)   & o3vmr_temp(*) = o3vmr_o
h2ovmr_temp  = fltarr(nlev)  & h2ovmr_temp(*) = h2ovmr_o

mwdry = n2vmr_o*mwn2 + co2vmr_o*mwco2 + ch4vmr_o*mwch4 +  h2vmr_o*mwh2  ;+ o2vmr_o*mwo2 + o3vmr_o*mwo3   

n2mmr_temp  = fltarr(nlev) & n2mmr_temp(*)  = n2vmr_temp(*)*mwn2/mwdry
h2mmr_temp  = fltarr(nlev) & h2mmr_temp(*)  = h2vmr_temp(*)*mwh2/mwdry
co2mmr_temp  = fltarr(nlev) & co2mmr_temp(*)  = co2vmr_temp(*)*mwco2/mwdry
ch4mmr_temp  = fltarr(nlev) & ch4mmr_temp(*)  = ch4vmr_temp(*)*mwch4/mwdry
h2ommr_temp  = fltarr(nlev) & h2ommr_temp(*)  = h2ovmr_temp(*)*mwh2o/mwdry

  
cpdry = n2mmr_temp(0)*cpn2 + co2mmr_temp(0)*cpco2 + ch4mmr_temp(0)*cpch4 +  h2mmr_temp(0)*cph2  ;+ o2vmr_o*cpo2 + o3vmr_o*cpo3   


;-- set final output values --
TS_out = ts_in
PS_out = ps_in
TMID_out = fltarr(nlev)   & TMID_out = tmid_in
TINT_out = fltarr(nilev)  & TINT_out = tint_in
PMID_out = fltarr(nlev)   & PMID_out = pmid_in 
PDEL_out = fltarr(nlev)   & PDEL_OUT = pdel_in
PINT_out = fltarr(nilev)  & PINT_OUT = pint_in
ZINT_OUT = fltarr(nilev)  & ZINT_out = zint_in
H2OMMR_out = fltarr(nlev) & H2OMMR_out(*) = h2ommr_temp(*)
CO2MMR_out = fltarr(nlev) & CO2MMR_out(*) = co2mmr_temp(*)
CH4MMR_out = fltarr(nlev) & CH4MMR_out(*) = ch4mmr_temp(*)
O2MMR_out = fltarr(nlev)  & O2MMR_out(*) = 0.0
O3MMR_out = fltarr(nlev)  & O3MMR_out(*) = 0.0
H2MMR_out = fltarr(nlev)  & H2MMR_out(*) = h2mmr_temp(*)
N2MMR_out = fltarr(nlev)  & N2MMR_out(*) = n2mmr_temp(*)
;ALDIR_out = ALDIR(lonVC, latVC)
;ALDIF_out = ALDIF(lonVC, latVC)
;ASDIR_out = ASDIR(lonVC, latVC)
;ASDIF_out = ASDIF(lonVC, latVC)
ALDIR_out = 0.25
ALDIF_out = 0.25
ASDIR_out = 0.25
ASDIF_out = 0.25
MW_OUT = mwdry
CP_OUT = cpdry

COSZRS_out = 0.5 ;; only matters for a solar computation
                     ;; does not matter for longwave computation


if (do_plot eq 1) then begin
  plot, tmid_out,pmid_out/100., /ylog, yrange=[5000.,0.01], psym=4,symsize=1.0, ystyle=1, xrange=[100,400], xstyle=1
  oplot, tmid_out,pmid_out/100., linestyle=1
  oplot, tint_out,pint_out/100., psym=1, symsize=1.0
  oplot, tint_out,pint_out/100., linestyle=0
endif


for k=0, nilev-1 do print, k+1, pint_in(k), tint_in(k), zint_in(k)
print, "---------------------------------------------------------"
for k=0, nlev-1 do print, k+1, pmid_in(k), tmid_in(k), pdel_in(k)
print, "---------------------------------------------------------"
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
;varid23 = NCDF_VARDEF(id,'FUL',dim1,/float)
;varid24 = NCDF_VARDEF(id,'FDL',dim1,/float)
;varid25 = NCDF_VARDEF(id,'FUS',dim1,/float)
;varid26 = NCDF_VARDEF(id,'FDS',dim1,/float)
;varid27 = NCDF_VARDEF(id,'QRS',dim2,/float)
;varid28 = NCDF_VARDEF(id,'QRL',dim2,/float)
;varid29 = NCDF_VARDEF(id,'FULC',dim1,/float)
;varid30 = NCDF_VARDEF(id,'FDLC',dim1,/float)
;varid31 = NCDF_VARDEF(id,'FUSC',dim1,/float)
;varid32 = NCDF_VARDEF(id,'FDSC',dim1,/float)
;varid33 = NCDF_VARDEF(id,'QRSC',dim2,/float)
;varid34 = NCDF_VARDEF(id,'QRLC',dim2,/float)


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
;NCDF_ATTPUT, id, varid23, "title", "diag :: FUL"     & NCDF_ATTPUT, id, varid23, "units", "W m-2"
;NCDF_ATTPUT, id, varid24, "title", "diag :: FDL"     & NCDF_ATTPUT, id, varid24, "units", "W m-2"
;NCDF_ATTPUT, id, varid25, "title", "diag :: FUS"     & NCDF_ATTPUT, id, varid25, "units", "W m-2"
;NCDF_ATTPUT, id, varid26, "title", "diag :: FDS"     & NCDF_ATTPUT, id, varid26, "units", "W m-2"
;NCDF_ATTPUT, id, varid27, "title", "diag :: QRS"      & NCDF_ATTPUT, id, varid27, "units", "K/day"
;NCDF_ATTPUT, id, varid28, "title", "diag :: QRL"      & NCDF_ATTPUT, id, varid28, "units", "K/day"
;NCDF_ATTPUT, id, varid29, "title", "diag :: FULC"     & NCDF_ATTPUT, id, varid29, "units", "W m-2"
;NCDF_ATTPUT, id, varid30, "title", "diag :: FDLC"     & NCDF_ATTPUT, id, varid30, "units", "W m-2"
;NCDF_ATTPUT, id, varid31, "title", "diag :: FUSC"     & NCDF_ATTPUT, id, varid31, "units", "W m-2"
;NCDF_ATTPUT, id, varid32, "title", "diag :: FDSC"     & NCDF_ATTPUT, id, varid32, "units", "W m-2"
;NCDF_ATTPUT, id, varid33, "title", "diag :: QRSC"      & NCDF_ATTPUT, id, varid33, "units", "K/day"
;NCDF_ATTPUT, id, varid34, "title", "diag :: QRLC"      & NCDF_ATTPUT, id, varid34, "units", "K/day"

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
;NCDF_VARPUT, id, varid23, FUL_OUT
;NCDF_VARPUT, id, varid24, FDL_OUT
;NCDF_VARPUT, id, varid25, FUS_OUT
;NCDF_VARPUT, id, varid26, FDS_OUT
;NCDF_VARPUT, id, varid27, QRS_OUT
;NCDF_VARPUT, id, varid28, QRL_OUT
;NCDF_VARPUT, id, varid29, FULC_OUT
;NCDF_VARPUT, id, varid30, FDLC_OUT
;NCDF_VARPUT, id, varid31, FUSC_OUT
;NCDF_VARPUT, id, varid32, FDSC_OUT
;NCDF_VARPUT, id, varid33, QRSC_OUT
;NCDF_VARPUT, id, varid34, QRLC_OUT

NCDF_CLOSE, id

FINISH:

end
