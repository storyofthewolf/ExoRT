pro makeColumn
;------------------------------------------------
;Author: Wolf, E.T.
;Revision history
;6/30/2020 make input profile for RT standalone
; 2022 --- included clouds in profile
;
;DESCRIPTION:  
;Manually create a 1D column to feed into exoRT
;RTprofile_in.nc
;
;------------------------------------------------

filename = 'RTprofile_in.nc'

do_plot = 0
do_override = 1
do_clouds = 1
do_copy_profile_to_rtrun = 1
do_carma = 0

;;=================== Specification of Profile ===========================

;---------------------------------------------------------------
; USER SET PROFILES SELECTION IN profile_data.pro
; Profile input to this subroutine here
;---------------------------------------------------------------
profile_data, pint_in, pmid_in, tint_in, tmid_in, q_in, nlev, nilev

;check inouts
;print, "pint_in", pint_in
;print, "pmid_in", pmid_in
;print, "tint_in", tint_in
;print, "tmid_in", tmid_in
;print, "nlev_out", nlev_out
;print, "nilev_out", nilev_out

ts_in   = tint_in(nilev-1)
ps_in   = pint_in(nilev-1)

print, "----------------"
print, "TS: ", ts_in
print, "PS: ", ps_in
print, "nlevs: ", nlev


;-----------------------------------------
; USER SET GRAVITY HERE  (it matters!)
;-----------------------------------------
; Earth
;GRAV = 9.80616
;R = 287.04

; Mars
GRAV = 3.7
R = 188.92


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


;---------------------------------
;USER SET MIXING RATIOS HERE
;---------------------------------

; USE DRY MIXING RATIOS for non-H2O species
co2vmr_o = 1.0
ch4vmr_o = 0.0
h2vmr_o = 0.0
o2vmr_o = 0.0
o3vmr_o = 0.0
n2vmr_o = 0.0 

; SET SPECIFIC HUMIDITY (Kg water / Kg total air mass)
h2o_spchum_o = 0.0

; calculate mixing ratios, etc
; mass of the atmosphere
mwn2 = 28.
mwar = 40.
mwco2 = 44.01
mwch4 = 16.
mwh2o = 18.
mwh2 = 2.

cpn2 = 1.039e3
;cpco2 = 0.846e3
cpco2 = 0.751e3  ;0.735e3_r8 
cpch4 = 2.226e3
cph2 = 14.32e3

; set profiles with uniform values
n2vmr_temp  = fltarr(nlev)  & n2vmr_temp(*) = n2vmr_o
co2vmr_temp = fltarr(nlev)  & co2vmr_temp(*) = co2vmr_o
ch4vmr_temp = fltarr(nlev)  & ch4vmr_temp(*) = ch4vmr_o
n2vmr_temp  = fltarr(nlev)   & n2vmr_temp(*) = n2vmr_o
h2vmr_temp  = fltarr(nlev)   & h2vmr_temp(*) = h2vmr_o
o2vmr_temp  = fltarr(nlev)   & o2vmr_temp(*) = o2vmr_o
o3vmr_temp  = fltarr(nlev)   & o3vmr_temp(*) = o3vmr_o


; calculate dry molecular weight of air
mwdry = n2vmr_o*mwn2 + co2vmr_o*mwco2 + ch4vmr_o*mwch4 +  h2vmr_o*mwh2  ;+ o2vmr_o*mwo2 + o3vmr_o*mwo3   

; set mass mixing ratios of dry components
n2mmr_temp  = fltarr(nlev) & n2mmr_temp(*)  = n2vmr_temp(*)*mwn2/mwdry
h2mmr_temp  = fltarr(nlev) & h2mmr_temp(*)  = h2vmr_temp(*)*mwh2/mwdry
co2mmr_temp  = fltarr(nlev) & co2mmr_temp(*)  = co2vmr_temp(*)*mwco2/mwdry
ch4mmr_temp  = fltarr(nlev) & ch4mmr_temp(*)  = ch4vmr_temp(*)*mwch4/mwdry

; calculate dry specific heat of air 
cpdry = n2mmr_temp(0)*cpn2 + co2mmr_temp(0)*cpco2 + ch4mmr_temp(0)*cpch4 +  h2mmr_temp(0)*cph2  ;+ o2vmr_o*cpo2 + o3vmr_o*cpo3   

;------------------------------------
; SET specific humidity                                                                                                 
;------------------------------------
; to fixed value                                                                                                        
;h2ommr_temp  = fltarr(nlev) & h2ommr_temp(*)  = h2o_spchum_o(*)                                                        
; to input profile                                                                                                      
h2ommr_temp  = fltarr(nlev) & h2ommr_temp(*)  = q_in(*)


;------------------------------------
;-- USER SET CLOUDS 
;------------------------------------
cliqwp_temp      = fltarr(nlev)  & cliqwp_temp(*) = 0.0
cicewp_temp      = fltarr(nlev)  & cicewp_temp(*) = 0.0
cicewp_co2_temp  = fltarr(nlev)  & cicewp_co2_temp(*) = 0.0
rel_temp      = fltarr(nlev)  & rel_temp(*) = 0.0
rei_temp      = fltarr(nlev)  & rei_temp(*) = 0.0
rei_co2_temp  = fltarr(nlev)  & rei_co2_temp(*) = 0.0

;52 ~ 500 mb

; centered at 25 Km, ~5 km thick using 2bar CO2 smart profile
cicewp_co2_temp(41) = 0.1  ; gm-2
cicewp_co2_temp(42) = 1.0  ; gm-2
cicewp_co2_temp(43) = 10.0  ; gm-2
cicewp_co2_temp(44) = 1.0  ; gm-2
cicewp_co2_temp(45) = 0.1  ; gm-2

rei_co2_temp(41:48) = 169.   ; microns

cicewp_co2_temp(*) = cicewp_co2_temp(*) * 100.


;cicewp_co2_temp(40:51) = [2.111749e-16, 3.564572e-12, 4.869048e-08, 0.001089489, 32.9725, 215.8229, $
;                          292.2275, 1420.496, 33.85143, 0.001785828, 8.611823e-06, 2.9916e-08]

;rei_co2_temp(40:51) =  [6.213649, 6.246476, 6.002151, 6.01229, 6.109071, $ 
;   5.980229, 6.200058, 6.768147, 9.319643, 14.55946, 23.35444, 30.94548]

;------------------------------------
;-- USER SET CARMA AEROSOLS
;------------------------------------
;NELEM = 1  ; number of elements
;NBIN  = 1   ; number of bins
;carmammr_temp = fltarr(nlev, nlem, nibn)  & carmammr_temp(*,*,*) = 0.0

;---------------------------------------------
;-- USER SET SURFACE ALBEDOS and EMISSIVITY
;---------------------------------------------
; Remember these are broadband integrated quantitites
; applied to each radiation stream seperately.
; Thus a + e need not equal 1.

aldir_o = 0.25      ; visible direct
aldif_o = 0.25     ; visible diffuse
asdir_o = 0.25      ; visible direct
asdif_o = 0.25     ; visible diffuse
srf_emiss_o = 1.0 ; thermal emissivity


;----------------------------------
;-- FINAL OUTPUT VALUES SET HERE
;----------------------------------
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
ALDIR_out = aldir_o
ALDIF_out = aldif_o
ASDIR_out = asdir_o
ASDIF_out = asdif_o
SRF_EMISS_out = srf_emiss_o
MW_OUT = mwdry
CP_OUT = cpdry
if (do_clouds eq 1) then begin
  CLIQWP_OUT = fltarr(nlev)  & CLIQWP_OUT(*) = cliqwp_temp(*)
  CICEWP_OUT = fltarr(nlev)  & CICEWP_OUT(*) = cicewp_temp(*)
  CICEWP_CO2_OUT = fltarr(nlev)  & CICEWP_CO2_OUT(*) = cicewp_co2_temp(*)
  REL_OUT = fltarr(nlev)  & REL_OUT(*) = rel_temp(*)
  REI_OUT = fltarr(nlev)  & REI_OUT(*) = rei_temp(*)
  REI_CO2_OUT = fltarr(nlev)  & REI_CO2_OUT(*) = rei_co2_temp(*)
endif
;if (do_carma eq 1) then begin
;  CARMAMMR_out = carmammr_temp
;endif
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


;-------------------------------------------
;----- CREATE RTprofle_in.nc 
;-------------------------------------------

print, "creating ", filename
id = NCDF_CREATE(filename,  /CLOBBER)

dim1 = NCDF_DIMDEF(id, 'pverp', nilev)
dim2 = NCDF_DIMDEF(id, 'pver', nlev)
dim3 = NCDF_DIMDEF(id, 'one', 1)
if (do_carma eq 1) then begin
  dim4 = NCDF_DIMDEF(id, 'nelem', nelem)
  dim5 = NCDF_DIMDEF(id, 'nbin', nbin)
endif
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
varid20 = NCDF_VARDEF(id,'srf_emiss',dim3,/float)
varid21 = NCDF_VARDEF(id,'coszrs',dim3,/float)
varid22 = NCDF_VARDEF(id,'mw',dim3,/float)
varid23 = NCDF_VARDEF(id,'cp',dim3,/float)
if (do_clouds) then begin
  varid24 = NCDF_VARDEF(id,'cliqwp',dim2,/float)
  varid25 = NCDF_VARDEF(id,'cicewp',dim2,/float)
  varid26 = NCDF_VARDEF(id,'cicewp_co2',dim2,/float)
  varid27 = NCDF_VARDEF(id,'rel',dim2,/float)
  varid28 = NCDF_VARDEF(id,'rei',dim2,/float)
  varid29 = NCDF_VARDEF(id,'rei_co2',dim2,/float)
endif
;if (do_carma) then begin
;  varid30 = NCDF_VARDEF(id,'carmammr',dim2,dim4,dim5,/float)
;endif

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
NCDF_ATTPUT, id, varid20, "title", "surface emissivity"           & NCDF_ATTPUT, id, varid19, "units", "emissivity"
NCDF_ATTPUT, id, varid21, "title", "cosine of zenith angle"       & NCDF_ATTPUT, id, varid20, "units", "none"
NCDF_ATTPUT, id, varid22, "title", "molecular weight of dry air"  & NCDF_ATTPUT, id, varid21, "units", "AMU"
NCDF_ATTPUT, id, varid23, "title", "specific heat of dry air"     & NCDF_ATTPUT, id, varid22, "units", "J/kg K"
if (do_clouds eq 1) then begin 
  NCDF_ATTPUT, id, varid24, "title", "in-cloud liquid water"      & NCDF_ATTPUT, id, varid24, "units", "g/m2"
  NCDF_ATTPUT, id, varid25, "title", "in-cloud ice water"         & NCDF_ATTPUT, id, varid25, "units", "g/m2"
  NCDF_ATTPUT, id, varid26, "title", "in-cloud CO2 ice"           & NCDF_ATTPUT, id, varid26, "units", "g/m2"
  NCDF_ATTPUT, id, varid27, "title", "liquid water cloud Reff"    & NCDF_ATTPUT, id, varid27, "units", "microns"
  NCDF_ATTPUT, id, varid28, "title", "ice water cloud Reff"       & NCDF_ATTPUT, id, varid28, "units", "microns"
  NCDF_ATTPUT, id, varid29, "title", "CO2 ice cloud Reff"         & NCDF_ATTPUT, id, varid29, "units", "microns"
endif
if (do_carma) then begin
  NCDF_ATTPUT, id, varid30, "title", "carma aerosol mixing ratios" & NCDF_ATTPUT, id, varid30, "units", "mmr"
endif

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
NCDF_VARPUT, id, varid20, SRF_EMISS_OUT
NCDF_VARPUT, id, varid21, COSZRS_OUT
NCDF_VARPUT, id, varid22, MW_OUT
NCDF_VARPUT, id, varid23, CP_OUT
if (do_clouds eq 1) then begin
  NCDF_VARPUT, id, varid24, CLIQWP_OUT
  NCDF_VARPUT, id, varid25, CICEWP_OUT
  NCDF_VARPUT, id, varid26, CICEWP_CO2_OUT
  NCDF_VARPUT, id, varid27, REL_OUT
  NCDF_VARPUT, id, varid28, REI_OUT
  NCDF_VARPUT, id, varid29, REI_CO2_OUT
endif
if (do_carma) then begin
  NCDF_VARPUT, id, varid30, CARMAMMR_OUT
endif

NCDF_CLOSE, id

if (do_copy_profile_to_rtrun eq 1) then spawn, "cp RTprofile_in.nc /gpfsm/dnb53/etwolf/models/ExoRT/run"


FINISH:



end
