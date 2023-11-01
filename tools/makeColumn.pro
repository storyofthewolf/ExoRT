pro makeColumn
;------------------------------------------------
; Author: Wolf, E.T.
; Revision history
; 6/30/2020 make input profile for RT standalone
;  1/11/2022 P-T input profiles included in seperate IDL module
; 
; DESCRIPTION:  make RTprofile_in.nc from scratch.
; make sure to set your P-T input profiles from exteneral sources
;
;
;------------------------------------------------


filename = 'RTprofile_in_TS273_10000ppm_c2h6.nc'
do_T_plot = 0
do_verbose = 0

;;=================== Specification of Profile ===========================

;---------------------------------
; USER SET MIXING RATIOS HERE
;---------------------------------
; USE DRY MIXING RATIOS for non-H2O species
co2vmr_o  = 0.01
ch4vmr_o  = 0.001
c2h6vmr_o = 1.0e-2
h2vmr_o   = 0.0
o2vmr_o   = 0.0
o3vmr_o   = 0.0
n2vmr_o   = 1.0 - co2vmr_o - ch4vmr_o - c2h6vmr_o

;---------------------------------------------------------------
; USER SET PROFILES SELECTION IN profile_data.pro
; Profile input to this subroutine here
;---------------------------------------------------------------
profile_data, pint_in, pmid_in, tint_in, tmid_in, q_in, nlev, nilev, tag_in

;check inouts
;print, "pint_in", pint_in
;print, "pmid_in", pmid_in
;print, "tint_in", tint_in
;print, "tmid_in", tmid_in
;print, "nlev_out", nlev_out
;print, "nilev_out", nilev_out

ts_in   = tint_in(nilev-1)
ps_in   = pint_in(nilev-1)

print, "----------------------------------------------"
print, "entering makeColumn.pro"
print, "----------------------------------------------"
print, "filename: ", filename
print, "profile_data tag: ", tag_in
print, "TS: ", ts_in
print, "PS: ", ps_in
print, "nlevs: ", nlev

print, "-- user input dry mixing ratios --"
print, "co2vmr ",co2vmr_o
print, "ch4vmr ",ch4vmr_o
print, "c2h6   ",c2h6vmr_o
print, "n2vmr  ",n2vmr_o

;-----------------------------------------
; USER SET GRAVITY HERE  (it matters!)
;-----------------------------------------
; Earth
GRAV = 9.80616
R = 287.04

; Mars
;GRAV = 3.7
;R = 188.92

;-----------------------------------------
; PDEL and ZINT automatically calculated
;-----------------------------------------
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


; SET SPECIFIC HUMIDITY (Kg water / Kg total air mass)
; otherwise read from templates profile_data.pro
;h2o_spchum_o = 0.000

; calculate mass mixing ratios, etc
; mass of the atmosphere
mwn2   = 28.
mwar   = 40.
mwco2  = 44.01
mwch4  = 16.
mwc2h6 = 30.
mwh2o  = 18.
mwh2   = 2.

cpn2   = 1.039e3
;cpco2 = 0.846e3
cpco2  = 0.751e3
cpch4  = 2.226e3
cpc2h6 = 1.75e3
cph2   = 14.32e3

; set profiles with uniform values
n2vmr_temp   = fltarr(nlev)   & n2vmr_temp(*)   = n2vmr_o
co2vmr_temp  = fltarr(nlev)   & co2vmr_temp(*)  = co2vmr_o
ch4vmr_temp  = fltarr(nlev)   & ch4vmr_temp(*)  = ch4vmr_o
c2h6vmr_temp = fltarr(nlev)   & c2h6vmr_temp(*) = c2h6vmr_o
n2vmr_temp   = fltarr(nlev)   & n2vmr_temp(*)   = n2vmr_o
h2vmr_temp   = fltarr(nlev)   & h2vmr_temp(*)   = h2vmr_o
o2vmr_temp   = fltarr(nlev)   & o2vmr_temp(*)   = o2vmr_o
o3vmr_temp   = fltarr(nlev)   & o3vmr_temp(*)   = o3vmr_o

; calculate dry molecular weight of air
mwdry = n2vmr_o*mwn2 + co2vmr_o*mwco2 + ch4vmr_o*mwch4 +  h2vmr_o*mwh2 + c2h6vmr_o*mwc2h6  ;+ o2vmr_o*mwo2 + o3vmr_o*mwo3   

; set mass mixing ratios of dry components
n2mmr_temp    = fltarr(nlev) & n2mmr_temp(*)   = n2vmr_temp(*)*mwn2/mwdry
h2mmr_temp    = fltarr(nlev) & h2mmr_temp(*)   = h2vmr_temp(*)*mwh2/mwdry
co2mmr_temp   = fltarr(nlev) & co2mmr_temp(*)  = co2vmr_temp(*)*mwco2/mwdry
ch4mmr_temp   = fltarr(nlev) & ch4mmr_temp(*)  = ch4vmr_temp(*)*mwch4/mwdry
c2h6mmr_temp  = fltarr(nlev) & c2h6mmr_temp(*) = c2h6vmr_temp(*)*mwc2h6/mwdry

; calculate dry specific heat of air
cpdry = n2mmr_temp(0)*cpn2 + co2mmr_temp(0)*cpco2 + ch4mmr_temp(0)*cpch4 + h2mmr_temp(0)*cph2 + c2h6mmr_temp(0)*cpc2h6  ;+ o2vmr_o*cpo2 + o3vmr_o*cpo3   

; set specific humidity
; to fixed value
;h2ommr_temp  = fltarr(nlev) & h2ommr_temp(*)  = h2o_spchum_o(*)
; to input profile
h2ommr_temp  = fltarr(nlev) & h2ommr_temp(*)  = q_in(*)


;----------------------------------
;-- USER SET SURFACE ALBEDOS
;----------------------------------
aldir_o = 0.25
aldif_o = 0.25
asdir_o = 0.25
asdif_o = 0.25


;----------------------------------
;-- FINAL OUTPUT VALUES SET HERE 
;----------------------------------
TS_out      = ts_in
PS_out      = ps_in
TMID_out    = fltarr(nlev)   & TMID_out       = tmid_in
TINT_out    = fltarr(nilev)  & TINT_out       = tint_in
PMID_out    = fltarr(nlev)   & PMID_out       = pmid_in 
PDEL_out    = fltarr(nlev)   & PDEL_OUT       = pdel_in
PINT_out    = fltarr(nilev)  & PINT_OUT       = pint_in
ZINT_OUT    = fltarr(nilev)  & ZINT_out       = zint_in
H2OMMR_out  = fltarr(nlev)   & H2OMMR_out(*)  = h2ommr_temp(*)
CO2MMR_out  = fltarr(nlev)   & CO2MMR_out(*)  = co2mmr_temp(*)
CH4MMR_out  = fltarr(nlev)   & CH4MMR_out(*)  = ch4mmr_temp(*)
C2H6MMR_out = fltarr(nlev)   & C2H6MMR_out(*) = c2h6mmr_temp(*)
O2MMR_out   = fltarr(nlev)   & O2MMR_out(*)   = 0.0
O3MMR_out   = fltarr(nlev)   & O3MMR_out(*)   = 0.0
H2MMR_out   = fltarr(nlev)   & H2MMR_out(*)   = h2mmr_temp(*)
N2MMR_out   = fltarr(nlev)   & N2MMR_out(*)   = n2mmr_temp(*)
ALDIR_out   = aldir_o
ALDIF_out   = aldif_o
ASDIR_out   = asdir_o
ASDIF_out   = asdif_o
MW_OUT      = mwdry
CP_OUT      = cpdry
COSZRS_out  = 0.5 ;; only matters for a solar computation
                  ;; does not matter for longwave computation

if (do_T_plot eq 1) then begin
  plot, tmid_out,pmid_out/100., /ylog, yrange=[5000.,0.01], psym=4,symsize=1.0, ystyle=1, xrange=[100,400], xstyle=1
  oplot, tmid_out,pmid_out/100., linestyle=1
  oplot, tint_out,pint_out/100., psym=1, symsize=1.0
  oplot, tint_out,pint_out/100., linestyle=0
endif

if (do_verbose eq 1) then begin
  print, "interface layers"
  print, "---------------------------------------------------------"
  for k=0, nilev-1 do print, k+1, pint_in(k), tint_in(k), zint_in(k)
  print, "mid layers"
  print, "---------------------------------------------------------"
  for k=0, nlev-1 do print, k+1, pmid_in(k), tmid_in(k), pdel_in(k)
endif
; ----- create RTprofle_in.nc --------

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
varid9 = NCDF_VARDEF(id,'asdir',dim3,/float)
varid10 = NCDF_VARDEF(id,'asdif',dim3,/float)
varid11 = NCDF_VARDEF(id,'aldir',dim3,/float)
varid12 = NCDF_VARDEF(id,'aldif',dim3,/float)
varid13 = NCDF_VARDEF(id,'coszrs',dim3,/float)
varid14 = NCDF_VARDEF(id,'mw',dim3,/float)
varid15 = NCDF_VARDEF(id,'cp',dim3,/float)

varid16 = NCDF_VARDEF(id,'n2mmr',dim2,/float)
varid17 = NCDF_VARDEF(id,'h2mmr',dim2,/float)
varid18 = NCDF_VARDEF(id,'o2mmr',dim2,/float)
varid19 = NCDF_VARDEF(id,'o3mmr',dim2,/float)
varid20 = NCDF_VARDEF(id,'h2ommr',dim2,/float)
varid21 = NCDF_VARDEF(id,'co2mmr',dim2,/float)
varid22 = NCDF_VARDEF(id,'ch4mmr',dim2,/float)
varid23 = NCDF_VARDEF(id,'c2h6mmr',dim2,/float)

NCDF_ATTPUT, id, varid1, "title", "Surface temperature"           & NCDF_ATTPUT, id, varid1, "units", "K"
NCDF_ATTPUT, id, varid2, "title", "Surface pressure"              & NCDF_ATTPUT, id, varid2, "units", "Pa"
NCDF_ATTPUT, id, varid3, "title", "Mid-layer temperatures"        & NCDF_ATTPUT, id, varid3, "units", "K"
NCDF_ATTPUT, id, varid4, "title", "Interface temperatures"        & NCDF_ATTPUT, id, varid4, "units", "K"
NCDF_ATTPUT, id, varid5, "title", "Mid-layer pressures"           & NCDF_ATTPUT, id, varid5, "units", "Pa"
NCDF_ATTPUT, id, varid6, "title", "Pressures difference"          & NCDF_ATTPUT, id, varid6, "units", "Pa"
NCDF_ATTPUT, id, varid7, "title", "Interface pressures"           & NCDF_ATTPUT, id, varid7, "units", "Pa"
NCDF_ATTPUT, id, varid8, "title", "Interface heights"             & NCDF_ATTPUT, id, varid8, "units", "m"
NCDF_ATTPUT, id, varid9, "title", "albedo shortwave direct"       & NCDF_ATTPUT, id, varid9, "units", "albedo"
NCDF_ATTPUT, id, varid10, "title", "albedo shortwave diffuse"     & NCDF_ATTPUT, id, varid10, "units", "albedo"
NCDF_ATTPUT, id, varid11, "title", "albedo longwave direct"       & NCDF_ATTPUT, id, varid11, "units", "albedo"
NCDF_ATTPUT, id, varid12, "title", "albedo longwave diffuse"      & NCDF_ATTPUT, id, varid12, "units", "albedo"
NCDF_ATTPUT, id, varid13, "title", "cosine of zenith angle"       & NCDF_ATTPUT, id, varid13, "units", "none"
NCDF_ATTPUT, id, varid14, "title", "molecular weight of dry air"  & NCDF_ATTPUT, id, varid14, "units", "AMU"
NCDF_ATTPUT, id, varid15, "title", "specific heat of dry air"     & NCDF_ATTPUT, id, varid15, "units", "J/kg K"
NCDF_ATTPUT, id, varid16, "title", "N2 MMR"                       & NCDF_ATTPUT, id, varid16, "units", "mmr"
NCDF_ATTPUT, id, varid17, "title", "H2 MMR"                       & NCDF_ATTPUT, id, varid17, "units", "mmr"
NCDF_ATTPUT, id, varid18, "title", "O2 MMR"                       & NCDF_ATTPUT, id, varid18, "units", "mmr"
NCDF_ATTPUT, id, varid19, "title", "O3 MMR"                       & NCDF_ATTPUT, id, varid19, "units", "mmr"
NCDF_ATTPUT, id, varid20, "title", "H2O MMR"                      & NCDF_ATTPUT, id, varid20, "units", "mmr"
NCDF_ATTPUT, id, varid21, "title", "CO2 MMR"                      & NCDF_ATTPUT, id, varid21, "units", "mmr"
NCDF_ATTPUT, id, varid22, "title", "CH4 MMR"                      & NCDF_ATTPUT, id, varid22, "units", "mmr"
NCDF_ATTPUT, id, varid23, "title", "C2H6 MMR"                     & NCDF_ATTPUT, id, varid23, "units", "mmr"

NCDF_CONTROL, id, /ENDEF

NCDF_VARPUT, id, varid1,  TS_OUT
NCDF_VARPUT, id, varid2,  PS_OUT
NCDF_VARPUT, id, varid3,  TMID_OUT
NCDF_VARPUT, id, varid4,  TINT_OUT
NCDF_VARPUT, id, varid5,  PMID_OUT
NCDF_VARPUT, id, varid6,  PDEL_OUT
NCDF_VARPUT, id, varid7,  PINT_OUT
NCDF_VARPUT, id, varid8,  ZINT_OUT
NCDF_VARPUT, id, varid9,  ASDIR_OUT
NCDF_VARPUT, id, varid10, ASDIF_OUT
NCDF_VARPUT, id, varid11, ALDIR_OUT
NCDF_VARPUT, id, varid12, ALDIF_OUT
NCDF_VARPUT, id, varid13, COSZRS_OUT
NCDF_VARPUT, id, varid14, MW_OUT
NCDF_VARPUT, id, varid15, CP_OUT
NCDF_VARPUT, id, varid16, N2MMR_OUT
NCDF_VARPUT, id, varid17, H2MMR_OUT
NCDF_VARPUT, id, varid18, O2MMR_OUT
NCDF_VARPUT, id, varid19, O3MMR_OUT
NCDF_VARPUT, id, varid20, H2OMMR_OUT
NCDF_VARPUT, id, varid21, CO2MMR_OUT
NCDF_VARPUT, id, varid22, CH4MMR_OUT
NCDF_VARPUT, id, varid23, C2H6MMR_OUT

NCDF_CLOSE, id

FINISH:

print, "----------------------------------------------"
print, "exiting makeColumn.pro"
print, "----------------------------------------------"

end
