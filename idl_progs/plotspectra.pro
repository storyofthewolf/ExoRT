pro plotspectra
;------------------------------------------------
;Author: Wolf, E.T.
;Created: August 2020
;Revision history: adpated from older code
;                                                         
;DESCRIPTION:
; Creates spectral histrogram style plots of the 
; shortwave and longwave calculations of ExoRT
;------------------------------------------------





do_28 = 0
do_42 = 0
do_68 = 1

plot_ps = 1 ; if eq 0, then plot to x windows
            ; if eq 1, then plot to postscript
if (plot_ps eq 1) then begin
  print, "plotting to postscript"
endif else begin
  print, "plotting to x-windows"
endelse
wait, 5

loadct, 33
nfiles = 1

fname = strarr(nfiles)
fname_short = strarr(nfiles)

fname_short(0) = " "
fname(0) =  "/gpfsm/dnb53/etwolf/models/ExoRT/run/RTprofile_out.nc"
;fname(0) =  "/gpfsm/dnb53/etwolf/models/ExoRT/run/RTprofile_out_TS340_H2Oonly_mtckd.nc"
;fname(0) =  "/gpfsm/dnb53/etwolf/models/ExoRT/run/RTprofile_out_TS340_H2Oonly_nocont.nc"


line_index = [0,0,0,0,0,0,2,2,2,2,2,2]
ncolors=6
colorscale = floor(255/ncolors)
color_index =(floor(findgen(ncolors)*255/ncolors))



if( do_28 eq 1) then begin
;; for 28 bin 
ntot_wavlnrng = 28
wavenum_edge =  [  $          ; all wavenumber edges [cm-1]
            10., 350., 500., 630., 700., 820., 980., 1100., 1180., 1390., $
            1480., 1800., 2080., 2200., 2380., 2600., 3250., 4000., 4650., $
            5150., 6150., 7700., 8050., 12850., 16000., 22650., 29000., $
            38000., 50000. ]
endif

if( do_42 eq 1) then begin
;; for 42 bin
ntot_wavlnrng=42                                                                                                 
wavenum_edge = [      10.0,       200.0,         350.0,        425.0, $
                500.0,       630.0,         700.0,        820.0, $ 
                980.0,       1100.0,        1180.0,       1390.0, $
               1480.0,       1800.0,        2080.0,       2200.0, $
               2380.0,       2600.0,        3250.0,       4000.0, $
               4650.0,       4900.0,        5150.0,       5650.0, $                                              
               6150.0,       6650.0,        7000.0,       7700.0, $
               8050.0,       9100.0,        10000.0,      11000.0, $
              11800.0,       12850.0,       13450.0,      14450.0, $
              15150.0,       16000.0,       19300.0,      22650.0, $
              29000.0,       38000.0,       50000.0 ]
endif


if( do_68 eq 1) then begin
;; for 68 bin
ntot_wavlnrng=68                                                                                                 
wavenum_edge = [0.00E+00,   40.00000,   100.0000,   160.0000, $
                220.0000,   280.0000,   330.0000,   380.0000, $
                440.0000,   495.0000,   545.0000,   617.0000, $
                667.0000,   720.0000,   800.0000,   875.0000, $
                940.0000,   1000.000,   1065.000,   1108.000, $
                1200.000,   1275.000,   1350.000,   1450.000, $
                1550.000,   1650.000,   1750.000,   1850.000, $
                1950.000,   2050.000,   2200.000,   2397.000, $
                2494.000,   2796.000,   3087.000,   3425.000, $
                3760.000,   4030.000,   4540.000,   4950.000, $
                5370.000,   5925.000,   6390.000,   6990.000, $
                7650.000,   8315.000,   8850.000,   9350.000, $
                9650.000,   10400.00,   11220.00,   11870.00, $
                12790.00,   13300.00,   14470.00,   15000.00, $
                16000.00,   16528.00,   17649.00,   18198.00, $
                18518.00,   22222.00,   25641.00,   29308.00, $
                30376.00,   32562.00,   35087.00,   36363.00,  42087.00 ]

endif



wavln_edge = 1.0e4/wavenum_edge
wavln_diff = wavln_edge(0:ntot_wavlnrng-1)-wavln_edge(1:ntot_wavlnrng)

wavenum_mid = (wavenum_edge(0:ntot_wavlnrng-2)+wavenum_edge(1:ntot_wavlnrng-1))/2.
wavln_mid = 1.0e4/wavenum_mid

wavenum_diff = (wavenum_edge(0:ntot_wavlnrng-1)-wavenum_edge(1:ntot_wavlnrng)) * (-1.)

SWABS_MATRIX = fltarr(nfiles)
ALBEDO_MATRIX = fltarr(nfiles)



ncid=ncdf_open(fname(0), /nowrite)
ncdf_varget,ncid,'TMID',TMID_IN
ncdf_varget,ncid,'TINT',TINT_IN
ncdf_varget,ncid,'PMID',PMID_IN       ; wet pressures, includes water vapor
ncdf_varget,ncid,'PINT',PINT_IN
ncdf_varget,ncid,'ZINT',ZINT_IN
ncdf_varget,ncid,'H2OMMR',H2OMMR_IN   ; specific humidity kg(wv) / (kg(dryair)+kg(wv)) 
ncdf_varget,ncid,'CO2MMR',CO2MMR_IN   ; dry mass mixing ratio kg(co2) / kg(dry air)    
ncdf_varget,ncid,'CH4MMR',CH4MMR_IN   ; dry mass mixing ratio kg(ch4) / kg(dry air)    
ncdf_varget,ncid,'O2MMR',O2MMR_IN     ; dry mass mixing ratio kg(o2) / kg(dry air)     
ncdf_varget,ncid,'O3MMR',O3MMR_IN     ; dry mass mixing ratio kg(o3) / kg(dry air)     
ncdf_varget,ncid,'N2MMR',N2MMR_IN     ; dry mass mixing ratio kg(n2) / kg(dry air)     
ncdf_varget,ncid,'H2MMR',H2MMR_IN     ; dry mass mixing ratio kg(h2) / kg(dry air)     
ncdf_varget,ncid,'LWUP',LWUP_IN
ncdf_varget,ncid,'LWDN',LWDN_IN
ncdf_varget,ncid,'SWUP',SWUP_IN
ncdf_varget,ncid,'SWDN',SWDN_IN
ncdf_varget,ncid,'LWUP_SPECTRAL',LWUP_SPECTRAL_IN
ncdf_varget,ncid,'LWDN_SPECTRAL',LWDN_SPECTRAL_IN
ncdf_varget,ncid,'SWUP_SPECTRAL',SWUP_SPECTRAL_IN
ncdf_varget,ncid,'SWDN_SPECTRAL',SWDN_SPECTRAL_IN
ncdf_varget,ncid,'LWHR',LWHR_IN
ncdf_varget,ncid,'SWHR',SWHR_IN
ncdf_varget,ncid,'MWDRY',MWDRY_IN
ncdf_varget,ncid,'CPDRY',CPDRY_IN
ncdf_close,ncid

pverp = n_elements(LWUP_IN)

for z=0, pverp-1 do begin
  lwup_level_sum = total(LWUP_SPECTRAL_IN(z,*))
  lwdn_level_sum = total(LWDN_SPECTRAL_IN(z,*))
  swup_level_sum = total(SWUP_SPECTRAL_IN(z,*))
  swdn_level_sum = total(SWDN_SPECTRAL_IN(z,*))
endfor
;print, fname(i)
;print, "LWUP TOA: ", total(LWUP_SPECTRAL_IN(0,*)), LWUP_IN(0)
;print, "LWDN TOA: ", total(LWDN_SPECTRAL_IN(0,*)), LWDN_IN(0)
;print, "SWUP TOA: ", total(SWUP_SPECTRAL_IN(0,*)), SWUP_IN(0)
;print, "SWDN TOA: ", total(SWDN_SPECTRAL_IN(0,*)), SWDN_IN(0)

;--- cumulative fluxes --

LWUP_TOA_CUM = fltarr(ntot_wavlnrng) & LWUP_TOA_CUM(ntot_wavlnrng-1) = LWUP_SPECTRAL_IN(0,ntot_wavlnrng-1)
LWDN_TOA_CUM = fltarr(ntot_wavlnrng) & LWDN_TOA_CUM(ntot_wavlnrng-1) = LWDN_SPECTRAL_IN(0,ntot_wavlnrng-1)
LWUP_SURFACE_CUM = fltarr(ntot_wavlnrng) & LWUP_SURFACE_CUM(ntot_wavlnrng-1) = LWUP_SPECTRAL_IN(pverp-1,ntot_wavlnrng-1)
LWDN_SURFACE_CUM = fltarr(ntot_wavlnrng) & LWDN_SURFACE_CUM(ntot_wavlnrng-1) = LWDN_SPECTRAL_IN(pverp-1,ntot_wavlnrng-1)



SWUP_TOA_CUM = fltarr(ntot_wavlnrng) & SWUP_TOA_CUM(0) = SWUP_SPECTRAL_IN(0,0)
SWDN_TOA_CUM = fltarr(ntot_wavlnrng) & SWDN_TOA_CUM(0) = SWDN_SPECTRAL_IN(0,0)
SWUP_SURFACE_CUM = fltarr(ntot_wavlnrng) & SWUP_SURFACE_CUM(0) = SWUP_SPECTRAL_IN(pverp-1,0)
SWDN_SURFACE_CUM = fltarr(ntot_wavlnrng) & SWDN_SURFACE_CUM(0) = SWDN_SPECTRAL_IN(pverp-1,0)



;print, "zcum, LWUP_TOA_CUM(zcum), LWUP_SURFACE_CUM(zcum), SWDN_TOA_CUM(zcum), SWDN_TOA_CUM(zcum), SWUP_SURFACE_CUM(zcum)"
;zcum=0
;print, zcum, LWUP_TOA_CUM(zcum), LWUP_SURFACE_CUM(zcum), SWDN_TOA_CUM(zcum), SWDN_TOA_CUM(zcum)-SWUP_SURFACE_CUM(zcum)

zrev=ntot_wavlnrng-1
for zcum=1, ntot_wavlnrng-1 do begin
  ; longwave
  zrev = zrev-1
  LWUP_TOA_CUM(zrev) = LWUP_TOA_CUM(zrev+1)+LWUP_SPECTRAL_IN(0,zrev)
  LWDN_TOA_CUM(zrev) = LWUP_TOA_CUM(zrev+1)+LWDN_SPECTRAL_IN(0,zrev)
  LWUP_SURFACE_CUM(zrev) = LWUP_SURFACE_CUM(zrev+1)+LWUP_SPECTRAL_IN(pverp-1,zrev)
  LWDN_SURFACE_CUM(zrev) = LWUP_SURFACE_CUM(zrev+1)+LWDN_SPECTRAL_IN(pverp-1,zrev)

  ;shortwave
  SWUP_TOA_CUM(zcum) = SWUP_TOA_CUM(zcum-1)+SWUP_SPECTRAL_IN(0,zcum)
  SWDN_TOA_CUM(zcum) = SWDN_TOA_CUM(zcum-1)+SWDN_SPECTRAL_IN(0,zcum)
  SWUP_SURFACE_CUM(zcum) = SWUP_SURFACE_CUM(zcum-1)+SWUP_SPECTRAL_IN(pverp-1,zcum)
  SWDN_SURFACE_CUM(zcum) = SWDN_SURFACE_CUM(zcum-1)+SWDN_SPECTRAL_IN(pverp-1,zcum)

endfor



print, "     LWUP TOA,   LWUP SURFACE,    SWDN TOA,   SW_ABS "
print, "------------------------------------------------------"
for z=0, ntot_wavlnrng-1 do begin
print, z,wavenum_edge(z), wavenum_edge(z+1),1.0e4/wavenum_edge(z), 1.0e4/wavenum_edge(z+1), LWUP_TOA_CUM(z), LWUP_SURFACE_CUM(z), SWDN_TOA_CUM(z), SWDN_TOA_CUM(z)-SWUP_TOA_CUM(z)
endfor


; cut offs for reduced RT
;lw_cut=17  ;(10 to 4000 cm-1)  (2.5 to 1000 µm)
;sw_cut=6   ;(820 to 50000 cm-1) (12.1915 to 0.2 µm)

lw_cut=ntot_wavlnrng  ;(10 to 4000 cm-1)  (2.5 to 1000 µm)
sw_cut=0   ;(820 to 50000 cm-1) (12.1915 to 0.2 µm)

print, "     LWUP TOA,   LWUP SURFACE,    SWDN TOA,   SW ABS "
print, LWUP_IN(0), LWUP_IN(pverp-1), SWDN_IN(0), SWDN_IN(0)-SWUP_IN(0)
print, total(LWUP_SPECTRAL_IN(0,0:lw_cut-1)), total(LWUP_SPECTRAL_IN(pverp-1,0:lw_cut-1)), total(SWDN_SPECTRAL_IN(0,sw_cut-1:ntot_wavlnrng-1)), total(SWDN_SPECTRAL_IN(0,sw_cut-1:ntot_wavlnrng-1)-SWUP_SPECTRAL_IN(0,sw_cut-1:ntot_wavlnrng-1))
print, total(LWUP_SPECTRAL_IN(0,0:lw_cut-1))/LWUP_IN(0), total(LWUP_SPECTRAL_IN(pverp-1,0:lw_cut-1))/LWUP_IN(pverp-1), total(SWDN_SPECTRAL_IN(0,sw_cut-1:ntot_wavlnrng-1))/SWDN_IN(0), total(SWDN_SPECTRAL_IN(0,sw_cut-1:ntot_wavlnrng-1)-SWUP_SPECTRAL_IN(0,sw_cut-1:ntot_wavlnrng-1))/(SWDN_IN(0)-SWUP_IN(0))


;----------------    plotting   ----------------
; shortwave plot
;-----------------------------------------------
!P.font=0


if (plot_ps eq 1) then begin 
  set_plot,'PS'
  device,file='idl_sw.ps'
  device,/color,BITS=8 ;, /ENCAPSULATED ;, /CMYK       
  device,xsize=18.5,ysize=15,/CM
  device, set_font='Helvetica-Oblique', FONT_INDEX=20
  device, set_font='Helvetica-Bold', FONT_INDEX=19
  device, set_font='helvetica',FONT_INDEX=18
endif else begin
  set_plot,'x'
endelse

xr = [0,5.0]
yr = [0.011,75]
plot, wavln_mid, SWDN_SPECTRAL_IN(pverp-1,*)/wavln_diff,/nodata, $
         xrange=xr, xstyle=1, yrange=yr, ystyle=1, xthick=3.0, ythick=3.0, $
         xtitle="wavelength (microns)", $
         ytitle="W m!U-2!N micron!U-1!N" ;, $
;         ytitle= "Cumumlative"


;create artificial bar
q=0
xbar = fltarr(ntot_wavlnrng*2)
ybar1 = fltarr(ntot_wavlnrng*2)
ybar2 = fltarr(ntot_wavlnrng*2)

for nw=0,ntot_wavlnrng-1 do begin

  ;x coordinateds
  xbar(q) = 1.0e4/wavenum_edge(nw)
  xbar(q+1) = 1.0e4/wavenum_edge(nw+1)
;  xbar(q) = wavenum_edge(nw)
;  xbar(q+1) = wavenum_edge(nw+1)

  ; short wave plots
  ; downwelling TOA
  ybar1(q) = SWDN_SPECTRAL_IN(0,nw)/wavln_diff(nw)
  ybar1(q+1) =  SWDN_SPECTRAL_IN(0,nw)/wavln_diff(nw)

  ; absorded solar
  ;ybar1(q) = SWDN_SPECTRAL_IN(0,nw)/wavln_diff(nw)-SWUP_SPECTRAL_IN(0,nw)/wavln_diff(nw)
  ;ybar1(q+1) =  SWDN_SPECTRAL_IN(0,nw)/wavln_diff(nw)-SWUP_SPECTRAL_IN(0,nw)/wavln_diff(nw)

  ; upwelling TOA
  ybar2(q) = SWUP_SPECTRAL_IN(0,nw)/wavln_diff(nw)
  ybar2(q+1) =  SWUP_SPECTRAL_IN(0,nw)/wavln_diff(nw)

  ; downwelling surface
;  ybar2(q) = SWDN_SPECTRAL_IN(pverp-1,nw)/wavln_diff(nw)
;  ybar2(q+1) =  SWDN_SPECTRAL_IN(pverp-1,nw)/wavln_diff(nw)


  q=q+2

endfor
  ci=254
  ; histogram plots
  oplot, xbar, ybar1, color=ci, thick=3.0, linestyle=2
  oplot, xbar, ybar1, color=ci, psym=4, symsize=0.5, thick=2.0

  oplot, xbar, ybar2, color=ci, thick=3.0, linestyle=0
  oplot, xbar, ybar2, color=ci, psym=4, symsize=0.5, thick=2.0



  ;=== longwave plots ===
  ;oplot, wavln_mid, LWUP_SPECTRAL_IN(0,*)/wavln_diff(*), color=color_index(i)
  ;oplot, wavln_mid, LWUP_SPECTRAL_IN(0,*)/wavln_diff(*), psym=4, color=color_index(i)

  ;oplot, wavln_mid, LWUP_SPECTRAL_IN(pverp-1,*)/wavln_diff(*), color=color_index(i), linestyle=1
  ;oplot, wavln_mid, LWUP_SPECTRAL_IN(pverp-1,*)/wavln_diff(*), psym=4, color=color_index(i)

  ;cumulative plots
;  oplot, wavln_mid, (SWDN_TOA_CUM(*)-SWUP_TOA_CUM(*))/total(SWDN_TOA_CUM(*)-SWUP_TOA_CUM(*)), color=color_index(i), linestyle=line_index(i), thick=3.0
;  oplot, wavln_mid, SWDN_TOA_CUM(*)-SWUP_TOA_CUM(*)/total(SWDN_TOA_CUM(*)-SWUP_TOA_CUM(*)), color=color_index(i), psym=4, symsize=0.5, thick=2.0
;  oplot, wavln_mid, (SWDN_TOA_CUM(*)-SWUP_TOA_CUM(*))/total(SWDN_IN(0,*)-SWUP_IN(0,*)), color=color_index(i), linestyle=line_index(i), thick=3.0
;  oplot, wavln_mid, (SWDN_TOA_CUM(*)-SWUP_TOA_CUM(*))/total(SWDN_IN(0,*)-SWUP_IN(0,*)), color=color_index(i), psym=4, symsize=0.5, thick=2.0


  ;oplot, wavln_mid, LWUP_TOA_CUM(*)/total(LWUP_IN(0,*)), color=color_index(i)
  ;oplot, wavln_mid, LWUP_TOA_CUM(*)/total(LWUP_IN(0,*)), psym=4, color=color_index(i)
  ;oplot, wavln_mid, LWUP_SURFACE_CUM(*)/total(LWUP_IN(pverp-1,*)), color=color_index(i)
  ;oplot, wavln_mid, LWUP_SURFACE_CUM(*)/total(LWUP_IN(pverp-1,*)), psym=4, color=color_index(i)

if  (plot_ps eq 1) then begin 
  device, /close
endif else begin
  stop
endelse
;----------------    plotting   ----------------
; longwave plot
;-----------------------------------------------
  !P.font=0

if (plot_ps eq 1) then begin
  set_plot,'PS'
  device,file='idl_lw.ps'
  device,/color,BITS=8 ;, /ENCAPSULATED ;, /CMYK       
  device,xsize=18.5,ysize=15,/CM
  device, set_font='Helvetica-Oblique', FONT_INDEX=20
  device, set_font='Helvetica-Bold', FONT_INDEX=19
  device, set_font='helvetica',FONT_INDEX=18
endif else begin
  set_plot,'x'
endelse

  xr = [0,2000]
  yr = [0,0.3]
  plot, wavln_mid, LWDN_SPECTRAL_IN(pverp-1,*)/wavln_diff,/nodata, $
         xrange=xr, xstyle=1, yrange=yr, ystyle=1, xthick=3.0, ythick=3.0, $
;         xtitle="wavenumber (cm-1)", $
         xtitle="wavelength (microns)", $
         ytitle="W m!U-2!N micron!U-1!N" ;, $
;         ytitle= "Cumumlative"

;create artificial bar
q=0
xbar = fltarr(ntot_wavlnrng*2)
ybar1 = fltarr(ntot_wavlnrng*2)
ybar2 = fltarr(ntot_wavlnrng*2)

for nw=0,ntot_wavlnrng-1 do begin

  ;x coordinateds
;  xbar(q) = 1.0e4/wavenum_edge(nw)
;  xbar(q+1) = 1.0e4/wavenum_edge(nw+1)
  xbar(q) = wavenum_edge(nw)
  xbar(q+1) = wavenum_edge(nw+1)


  ; long wave plots
  ybar1(q) = LWUP_SPECTRAL_IN(0,nw)/wavenum_diff(nw)
  ybar1(q+1) =  LWUP_SPECTRAL_IN(0,nw)/wavenum_diff(nw) ;wavln_diff(nw)

  ybar2(q) = LWUP_SPECTRAL_IN(pverp-1,nw)/wavenum_diff(nw)
  ybar2(q+1) =  LWUP_SPECTRAL_IN(pverp-1,nw)/wavenum_diff(nw)


  q=q+2

endfor
  ci=254
  ; histogram plots
  oplot, xbar, ybar1, color=ci, thick=3.0, linestyle=2
  oplot, xbar, ybar1, color=ci, psym=4, symsize=0.5, thick=2.0

  oplot, xbar, ybar2, color=ci, thick=3.0, linestyle=0
  oplot, xbar, ybar2, color=ci, psym=4, symsize=0.5, thick=2.0

oplot, wavenum_mid, LWUP_SPECTRAL_IN(pverp-1,*)/wavenum_diff(*)
oplot, wavenum_mid, LWUP_SPECTRAL_IN(0,*)/wavenum_diff(*)

;== solar plots ===
  ;downwelling surface
 ; oplot, wavln_mid, SWDN_SPECTRAL_IN(pverp-1,*)/wavln_diff, color=color_index(i), thick=3.0, linestyle=line_index(i)
 ; oplot, wavln_mid, SWDN_SPECTRAL_IN(pverp-1,*)/wavln_diff, color=color_index(i), psym=4, symsize=0.5, thick=2.0

  ;=== absorbed solar flux ===
  ;oplot, wavln_mid, SWDN_SPECTRAL_IN(0,*)/wavln_diff-SWUP_SPECTRAL_IN(0,*)/wavln_diff, color=color_index(i), linestyle=line_index(i), thick=3.0
  ;oplot, wavln_mid, SWDN_SPECTRAL_IN(0,*)/wavln_diff-SWUP_SPECTRAL_IN(0,*)/wavln_diff, color=color_index(i), psym=4, symsize=0.5, thick=2.0

  ;=== TOA downewlling ===
;  oplot, wavln_mid, SWDN_SPECTRAL_IN(0,*)/wavln_diff, color=color_index(i), linestyle=1, thick=3.0
;  oplot, wavln_mid, SWDN_SPECTRAL_IN(0,*)/wavln_diff, color=color_index(i), psym=4, symsize=0.5, thick=2.0

  ;cumulative plots
;  oplot, wavln_mid, (SWDN_TOA_CUM(*)-SWUP_TOA_CUM(*))/total(SWDN_TOA_CUM(*)-SWUP_TOA_CUM(*)), color=color_index(i), linestyle=line_index(i), thick=3.0
;  oplot, wavln_mid, SWDN_TOA_CUM(*)-SWUP_TOA_CUM(*)/total(SWDN_TOA_CUM(*)-SWUP_TOA_CUM(*)), color=color_index(i), psym=4, symsize=0.5, thick=2.0
;  oplot, wavln_mid, (SWDN_TOA_CUM(*)-SWUP_TOA_CUM(*))/total(SWDN_IN(0,*)-SWUP_IN(0,*)), color=color_index(i), linestyle=line_index(i), thick=3.0
;  oplot, wavln_mid, (SWDN_TOA_CUM(*)-SWUP_TOA_CUM(*))/total(SWDN_IN(0,*)-SWUP_IN(0,*)), color=color_index(i), psym=4, symsize=0.5, thick=2.0

if  (plot_ps eq 1) then device, /close



end

