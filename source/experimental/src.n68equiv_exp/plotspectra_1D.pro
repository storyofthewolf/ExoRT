pro plotspectra_1D
;------------------------------------------------
;Author: Wolf, E.T.
;Created: August 2020
;Revision history: adpated from older code
;                                                         
;DESCRIPTION:
; Creates spectral histrogram style plots of the 
; shortwave and longwave spectra from 
; 1D offline ExoRT calculations
;------------------------------------------------

do_clouds = 1
nfiles = 1
fname = strarr(nfiles)
fname_short = strarr(nfiles)
fname_short(0) = " "

level=45
; filename of radiation output
fname(0) =  "/gpfsm/dnb53/etwolf/models/ExoRT/run/RTprofile_out.nc"
;fname(0) =  "/gpfsm/dnb53/etwolf/models/ExoRT/run/RTprofile_out_TS340_H2Oonly_mtckd.nc"
;fname(0) =  "/gpfsm/dnb53/etwolf/models/ExoRT/run/RTprofile_out_TS340_H2Oonly_nocont.nc"

; number of spectral intervals in radiation output file
; choose one and only one
do_n28 = 0
do_n42 = 0
do_n68 = 1
do_n73 = 0
do_n84 = 0

; select whether to plot in x windows or postscript
plot_ps = 0 ; if eq 0, then plot to x windows
            ; if eq 1, then plot to postscript
if (plot_ps eq 1) then begin
  print, "plotting to postscript"
endif else begin
  print, "plotting to x-windows"
endelse
wait, 5

; plotting options
loadct, 33
line_index = [0,0,0,0,0,0,2,2,2,2,2,2]
ncolors=6
colorscale = floor(255/ncolors)
color_index =(floor(findgen(ncolors)*255/ncolors))


;------- define spectral interval grids -------
if( do_n28 eq 1) then begin
;; for 28 bin 
ntot_wavlnrng = 28
wavenum_edge =  [  $          ; all wavenumber edges [cm-1]
            10., 350., 500., 630., 700., 820., 980., 1100., 1180., 1390., $
            1480., 1800., 2080., 2200., 2380., 2600., 3250., 4000., 4650., $
            5150., 6150., 7700., 8050., 12850., 16000., 22650., 29000., $
            38000., 50000. ]
endif

if( do_n42 eq 1) then begin
;; for 42 bin
ntot_wavlnrng=42                                                                                                 
wavenum_edge = [  $   ; all wavenumber edges [cm-1] 
                 10.0,       200.0,         350.0,        425.0, $
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


if( do_n68 eq 1) then begin
;; for 68 bin
ntot_wavlnrng=68                                                                                                 
wavenum_edge = [  $  ; all wavenumber edges [cm-1]
                0.00E+00,   40.00000,   100.0000,   160.0000, $
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

if( do_n73 eq 1) then begin
;; for 73 bin
ntot_wavlnrng=73
wavenum_edge  = [  $   ; all wavenumber edges [cm-1]
                 0.00000, 40.0000, 100.000, 160.000, 220.000, 280.000, 330.000, $
                 380.000, 440.000, 495.000, 545.000, 617.000, 667.000, 720.000, $
                 800.000, 875.000, 940.000, 1000.00, 1065.00, 1108.00, 1200.00, $
                 1275.00, 1350.00, 1450.00, 1550.00, 1650.00, 1750.00, 1850.00, $
                 1950.00, 2050.00, 2200.00, 2439.02, 2564.10, 2777.78, 3174.60, $
                 3508.77, 3773.59, 4081.63, 4545.46, 4716.98, 5154.64, 5376.34, $
                 5555.56, 5952.38, 6172.84, 6578.95, 6711.41, 6849.31, 7042.25, $
                 7462.69, 7692.31, 8064.52, 8333.33, 8620.69, 8928.57, 9090.91, $
                 9259.26, 9708.74, 10869.6, 11111.1, 11363.6, 11494.3, 12500.0, $
                 12820.5, 14492.8, 16393.4, 18181.8, 20000.0, 22222.2, 23809.5, $
                 25974.0, 28985.5, 33333.3, 50000.0  ]
endif


if (do_n84 eq 1) then begin
;; for 84 bin
ntot_wavlnrng=84
wavenum_edge =  [  $  ; all wavenumber edges [cm-1]
                 0.00E+00,   40.00000,   100.0000,   160.0000, $
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
                 30376.00,   32562.00,   35087.00,   36363.00, $
                 42087.00,   50000.00,   60000.00,   70000.00, $
                 80000.00,   90000.00,  100000.00,  125000.00, $
                 150000.00, 175000.00,  200000.00,  300000.00, $
                 400000.00, 500000.00, 750000.00,  1000000.00, $
                 1250000.00 ]
endif


; cut offs for reduced RT
;lw_cut=17  ;(10 to 4000 cm-1)  (2.5 to 1000 µm)
;sw_cut=6   ;(820 to 50000 cm-1) (12.1915 to 0.2 µm)

lw_cut=ntot_wavlnrng  ;(10 to 4000 cm-1)  (2.5 to 1000 µm)
sw_cut=0   ;(820 to 50000 cm-1) (12.1915 to 0.2 µm)



wavln_edge = 1.0e4/wavenum_edge
wavln_diff = wavln_edge(0:ntot_wavlnrng-1)-wavln_edge(1:ntot_wavlnrng)

wavenum_mid = (wavenum_edge(0:ntot_wavlnrng-2)+wavenum_edge(1:ntot_wavlnrng-1))/2.
wavln_mid = 1.0e4/wavenum_mid

wavenum_diff = (wavenum_edge(0:ntot_wavlnrng-1)-wavenum_edge(1:ntot_wavlnrng)) * (-1.)

SWABS_MATRIX = fltarr(nfiles)
ALBEDO_MATRIX = fltarr(nfiles)


; read in radiation output netcdf file
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
if (do_clouds eq 1) then begin
  ncdf_varget,ncid,'LWUP_CLD',LWUP_CLD_IN
  ncdf_varget,ncid,'LWDN_CLD',LWDN_CLD_IN
  ncdf_varget,ncid,'SWUP_CLD',SWUP_CLD_IN
  ncdf_varget,ncid,'SWDN_CLD',SWDN_CLD_IN
  ncdf_varget,ncid,'LWUP_SPECTRAL_CLD',LWUP_SPECTRAL_CLD_IN
  ncdf_varget,ncid,'LWDN_SPECTRAL_CLD',LWDN_SPECTRAL_CLD_IN
  ncdf_varget,ncid,'SWUP_SPECTRAL_CLD',SWUP_SPECTRAL_CLD_IN
  ncdf_varget,ncid,'SWDN_SPECTRAL_CLD',SWDN_SPECTRAL_CLD_IN
  ncdf_varget,ncid,'LWHR_CLD',LWHR_CLD_IN
  ncdf_varget,ncid,'SWHR_CLD',SWHR_CLD_IN
endif

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


if (do_clouds eq 1) then begin
  LWUP_CLD_TOA_CUM = fltarr(ntot_wavlnrng) & LWUP_CLD_TOA_CUM(ntot_wavlnrng-1) = LWUP_SPECTRAL_CLD_IN(0,ntot_wavlnrng-1)
  LWDN_CLD_TOA_CUM = fltarr(ntot_wavlnrng) & LWDN_CLD_TOA_CUM(ntot_wavlnrng-1) = LWDN_SPECTRAL_CLD_IN(0,ntot_wavlnrng-1)
  LWUP_CLD_SURFACE_CUM = fltarr(ntot_wavlnrng) & LWUP_CLD_SURFACE_CUM(ntot_wavlnrng-1) = LWUP_SPECTRAL_CLD_IN(pverp-1,ntot_wavlnrng-1)
  LWDN_CLD_SURFACE_CUM = fltarr(ntot_wavlnrng) & LWDN_CLD_SURFACE_CUM(ntot_wavlnrng-1) = LWDN_SPECTRAL_CLD_IN(pverp-1,ntot_wavlnrng-1)

  SWUP_CLD_TOA_CUM = fltarr(ntot_wavlnrng) & SWUP_CLD_TOA_CUM(0) = SWUP_SPECTRAL_CLD_IN(0,0)
  SWDN_CLD_TOA_CUM = fltarr(ntot_wavlnrng) & SWDN_CLD_TOA_CUM(0) = SWDN_SPECTRAL_CLD_IN(0,0)
  SWUP_CLD_SURFACE_CUM = fltarr(ntot_wavlnrng) & SWUP_CLD_SURFACE_CUM(0) = SWUP_SPECTRAL_CLD_IN(pverp-1,0)
  SWDN_CLD_SURFACE_CUM = fltarr(ntot_wavlnrng) & SWDN_CLD_SURFACE_CUM(0) = SWDN_SPECTRAL_CLD_IN(pverp-1,0)
endif


;print, "zcum, LWUP_TOA_CUM(zcum), LWUP_SURFACE_CUM(zcum), SWDN_TOA_CUM(zcum), SWDN_TOA_CUM(zcum), SWUP_SURFACE_CUM(zcum)"
;zcum=0
;print, zcum, LWUP_TOA_CUM(zcum), LWUP_SURFACE_CUM(zcum), SWDN_TOA_CUM(zcum), SWDN_TOA_CUM(zcum)-SWUP_SURFACE_CUM(zcum)

zrev=ntot_wavlnrng-1
for zcum=1, ntot_wavlnrng-1 do begin

  zrev = zrev-1

  ; longwave
  LWUP_TOA_CUM(zrev) = LWUP_TOA_CUM(zrev+1)+LWUP_SPECTRAL_IN(0,zrev)
  LWDN_TOA_CUM(zrev) = LWUP_TOA_CUM(zrev+1)+LWDN_SPECTRAL_IN(0,zrev)
  LWUP_SURFACE_CUM(zrev) = LWUP_SURFACE_CUM(zrev+1)+LWUP_SPECTRAL_IN(pverp-1,zrev)
  LWDN_SURFACE_CUM(zrev) = LWUP_SURFACE_CUM(zrev+1)+LWDN_SPECTRAL_IN(pverp-1,zrev)

  ;shortwave
  SWUP_TOA_CUM(zcum) = SWUP_TOA_CUM(zcum-1)+SWUP_SPECTRAL_IN(0,zcum)
  SWDN_TOA_CUM(zcum) = SWDN_TOA_CUM(zcum-1)+SWDN_SPECTRAL_IN(0,zcum)
  SWUP_SURFACE_CUM(zcum) = SWUP_SURFACE_CUM(zcum-1)+SWUP_SPECTRAL_IN(pverp-1,zcum)
  SWDN_SURFACE_CUM(zcum) = SWDN_SURFACE_CUM(zcum-1)+SWDN_SPECTRAL_IN(pverp-1,zcum)

  if (do_clouds eq 1) then begin
    ;longwave  
    LWUP_CLD_TOA_CUM(zrev) = LWUP_CLD_TOA_CUM(zrev+1)+LWUP_SPECTRAL_CLD_IN(0,zrev)
    LWDN_CLD_TOA_CUM(zrev) = LWUP_CLD_TOA_CUM(zrev+1)+LWDN_SPECTRAL_CLD_IN(0,zrev)
    LWUP_CLD_SURFACE_CUM(zrev) = LWUP_CLD_SURFACE_CUM(zrev+1)+LWUP_SPECTRAL_CLD_IN(pverp-1,zrev)
    LWDN_CLD_SURFACE_CUM(zrev) = LWUP_CLD_SURFACE_CUM(zrev+1)+LWDN_SPECTRAL_CLD_IN(pverp-1,zrev)
  
    ;shortwave
    SWUP_CLD_TOA_CUM(zcum) = SWUP_CLD_TOA_CUM(zcum-1)+SWUP_SPECTRAL_CLD_IN(0,zcum)
    SWDN_CLD_TOA_CUM(zcum) = SWDN_CLD_TOA_CUM(zcum-1)+SWDN_SPECTRAL_CLD_IN(0,zcum)
    SWUP_CLD_SURFACE_CUM(zcum) = SWUP_CLD_SURFACE_CUM(zcum-1)+SWUP_SPECTRAL_CLD_IN(pverp-1,zcum)
    SWDN_CLD_SURFACE_CUM(zcum) = SWDN_CLD_SURFACE_CUM(zcum-1)+SWDN_SPECTRAL_CLD_IN(pverp-1,zcum)
  endif


endfor


print, "------------------------------------------------------"
print, "== clearsky fluxes"
print, "     LWUP TOA,   LWUP SURFACE,    SWDN TOA,   SW_ABS "
print, "------------------------------------------------------"
for z=0, ntot_wavlnrng-1 do begin
  print, z,wavenum_edge(z), wavenum_edge(z+1),1.0e4/wavenum_edge(z), 1.0e4/wavenum_edge(z+1), LWUP_TOA_CUM(z), LWUP_SURFACE_CUM(z), SWDN_TOA_CUM(z), SWDN_TOA_CUM(z)-SWUP_TOA_CUM(z)
endfor
print, "     LWUP TOA,   LWUP SURFACE,    SWDN TOA,   SW ABS "
print, LWUP_IN(0), LWUP_IN(pverp-1), SWDN_IN(0), SWDN_IN(0)-SWUP_IN(0)
print, total(LWUP_SPECTRAL_IN(0,0:lw_cut-1)), total(LWUP_SPECTRAL_IN(pverp-1,0:lw_cut-1)), total(SWDN_SPECTRAL_IN(0,sw_cut-1:ntot_wavlnrng-1)), total(SWDN_SPECTRAL_IN(0,sw_cut-1:ntot_wavlnrng-1)-SWUP_SPECTRAL_IN(0,sw_cut-1:ntot_wavlnrng-1))
print, total(LWUP_SPECTRAL_IN(0,0:lw_cut-1))/LWUP_IN(0), total(LWUP_SPECTRAL_IN(pverp-1,0:lw_cut-1))/LWUP_IN(pverp-1), total(SWDN_SPECTRAL_IN(0,sw_cut-1:ntot_wavlnrng-1))/SWDN_IN(0), total(SWDN_SPECTRAL_IN(0,sw_cut-1:ntot_wavlnrng-1)-SWUP_SPECTRAL_IN(0,sw_cut-1:ntot_wavlnrng-1))/(SWDN_IN(0)-SWUP_IN(0))

if (do_clouds eq 1) then begin
  print, "------------------------------------------------------"
  print, "== cloudy sky fluxes"
  print, "     LWUP TOA,   LWUP SURFACE,    SWDN TOA,   SW_ABS "
  print, "------------------------------------------------------"
  for z=0, ntot_wavlnrng-1 do begin
    print, z,wavenum_edge(z), wavenum_edge(z+1),1.0e4/wavenum_edge(z), 1.0e4/wavenum_edge(z+1), LWUP_CLD_TOA_CUM(z), LWUP_CLD_SURFACE_CUM(z), SWDN_CLD_TOA_CUM(z), SWDN_CLD_TOA_CUM(z)-SWUP_CLD_TOA_CUM(z)
  endfor
endif
print, "     LWUP TOA,   LWUP SURFACE,    SWDN TOA,   SW ABS "
print, LWUP_CLD_IN(0), LWUP_CLD_IN(pverp-1), SWDN_CLD_IN(0), SWDN_CLD_IN(0)-SWUP_CLD_IN(0)
print, total(LWUP_SPECTRAL_CLD_IN(0,0:lw_cut-1)), total(LWUP_SPECTRAL_CLD_IN(pverp-1,0:lw_cut-1)), total(SWDN_SPECTRAL_CLD_IN(0,sw_cut-1:ntot_wavlnrng-1)), total(SWDN_SPECTRAL_CLD_IN(0,sw_cut-1:ntot_wavlnrng-1)-SWUP_SPECTRAL_CLD_IN(0,sw_cut-1:ntot_wavlnrng-1))
print, total(LWUP_SPECTRAL_CLD_IN(0,0:lw_cut-1))/LWUP_CLD_IN(0), total(LWUP_SPECTRAL_CLD_IN(pverp-1,0:lw_cut-1))/LWUP_CLD_IN(pverp-1), total(SWDN_SPECTRAL_CLD_IN(0,sw_cut-1:ntot_wavlnrng-1))/SWDN_CLD_IN(0), total(SWDN_SPECTRAL_CLD_IN(0,sw_cut-1:ntot_wavlnrng-1)-SWUP_SPECTRAL_CLD_IN(0,sw_cut-1:ntot_wavlnrng-1))/(SWDN_CLD_IN(0)-SWUP_CLD_IN(0))




;print individual level diagnostics"
print, "------------------------------"
print, "level, ", level
for z=0, ntot_wavlnrng-1 do begin
  print, z, SWDN_SPECTRAL_IN(level,z), SWDN_SPECTRAL_CLD_IN(level,z)
endfor
print,total(SWDN_SPECTRAL_IN(level,*)), total(SWDN_SPECTRAL_CLD_IN(level,*))



;----------------    plotting   ----------------
; shortwave plot
;-----------------------------------------------
!P.font=0


if (plot_ps eq 1) then begin 
  set_plot,'PS'
  device,file='idl_sw.eps'
  device,/color,BITS=8, /ENCAPSULATED ;, /CMYK       
  device,xsize=18.5,ysize=15,/CM
  device, set_font='Helvetica-Oblique', FONT_INDEX=20
  device, set_font='Helvetica-Bold', FONT_INDEX=19
  device, set_font='helvetica',FONT_INDEX=18
endif else begin
  set_plot,'x'
endelse

xr = [0.0,5]
yr = [0.0,300]
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
ybar1_cld = fltarr(ntot_wavlnrng*2)
ybar2_cld = fltarr(ntot_wavlnrng*2)



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
  ;ybar2(q) = SWDN_SPECTRAL_IN(pverp-1,nw)/wavln_diff(nw)
  ;ybar2(q+1) =  SWDN_SPECTRAL_IN(pverp-1,nw)/wavln_diff(nw)

  if (do_clouds eq 1) then begin
    ; downwelling TOA
    ybar1_cld(q) = SWDN_SPECTRAL_CLD_IN(0,nw)/wavln_diff(nw)
    ybar1_cld(q+1) =  SWDN_SPECTRAL_CLD_IN(0,nw)/wavln_diff(nw)

    ; downwelling
    ybar2_cld(q) = SWDN_SPECTRAL_CLD_IN(level,nw)/wavln_diff(nw)
    ybar2_cld(q+1) =  SWDN_SPECTRAL_CLD_IN(level,nw)/wavln_diff(nw)

    ; absorded solar
    ;ybar1_cld(q) = SWDN_SPECTRAL_CLD_IN(0,nw)/wavln_diff(nw)-SWUP_SPECTRAL_CLD_IN(0,nw)/wavln_diff(nw)
    ;ybar1_Cld(q+1) =  SWDN_SPECTRAL_CLD_IN(0,nw)/wavln_diff(nw)-SWUP_SPECTRAL_CLD_IN(0,nw)/wavln_diff(nw)

    ; upwelling TOA
    ;ybar2_cld(q) = SWUP_SPECTRAL_CLD_IN(0,nw)/wavln_diff(nw)
    ;ybar2_cld(q+1) =  SWUP_SPECTRAL_CLD_IN(0,nw)/wavln_diff(nw)

    ; downwelling surface
    ;ybar2_cld(q) = SWDN_SPECTRAL_CLD_IN(pverp-1,nw)/wavln_diff(nw)
    ;ybar2_cld(q+1) =  SWDN_SPECTRAL_CLD_IN(pverp-1,nw)/wavln_diff(nw)
   endif


  q=q+2

endfor
  ci=254
  ; histogram plots
  oplot, xbar, ybar1, color=ci, thick=3.0, linestyle=0
  oplot, xbar, ybar1, color=ci, psym=4, symsize=0.5, thick=2.0

  oplot, xbar, ybar2, color=ci, thick=3.0, linestyle=0
  oplot, xbar, ybar2, color=ci, psym=4, symsize=0.5, thick=2.0
  
  if (do_clouds eq 1) then begin
    oplot, xbar, ybar1_cld, color=ci, thick=3.0, linestyle=2
    oplot, xbar, ybar1_cld, color=ci, psym=4, symsize=0.5, thick=2.0

    oplot, xbar, ybar2_cld, color=ci, thick=3.0, linestyle=2
    oplot, xbar, ybar2_cld, color=ci, psym=4, symsize=0.5, thick=2.0
  endif


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
  device,file='idl_lw.eps'
  device,/color,BITS=8, /ENCAPSULATED ;, /CMYK       
  device,xsize=18.5,ysize=15,/CM
  device, set_font='Helvetica-Oblique', FONT_INDEX=20
  device, set_font='Helvetica-Bold', FONT_INDEX=19
  device, set_font='helvetica',FONT_INDEX=18
endif else begin
  set_plot,'x'
endelse

  xr = [0,2500]
  yr = [0,0.6]
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

  if (do_clouds eq 1) then begin
    ybar1_cld(q) = LWUP_SPECTRAL_CLD_IN(level,nw)/wavenum_diff(nw)
    ybar1_cld(q+1) =  LWUP_SPECTRAL_CLD_IN(level,nw)/wavenum_diff(nw)
    ybar2_cld(q) = LWUP_SPECTRAL_CLD_IN(pverp-1,nw)/wavenum_diff(nw)
    ybar2_cld(q+1) =  LWUP_SPECTRAL_CLD_IN(pverp-1,nw)/wavenum_diff(nw)
  endif

  q=q+2

endfor
print, "LWUP, level", level, TOTAL(LWUP_SPECTRAL_CLD_IN(level,*))

  cib=90
  cir=240

  ; histogram plots
  oplot, xbar, ybar1, color=cir, thick=4.0, linestyle=0
  oplot, xbar, ybar1, color=cir, psym=4, symsize=0.5, thick=4.0

  oplot, xbar, ybar2, color=cib, thick=4.0, linestyle=0
  oplot, xbar, ybar2, color=cib, psym=4, symsize=0.5, thick=4.0

;  oplot, wavenum_mid, LWUP_SPECTRAL_CLD_IN(pverp-1,*)/wavenum_diff(*), color=cib, thick=4
;  oplot, wavenum_mid, LWUP_SPECTRAL_CLD_IN(0,*)/wavenum_diff(*), color=cir, thick=4

  if (do_clouds eq 1) then begin
    oplot, xbar, ybar1_cld, color=cir, thick=4.0, linestyle=2
    oplot, xbar, ybar1_cld, color=cir, psym=4, symsize=0.5, thick=4.0

    oplot, xbar, ybar2_cld, color=cib, thick=4.0, linestyle=2
    oplot, xbar, ybar2_cld, color=cib, psym=4, symsize=0.5, thick=4.0
  endif

;  xyouts, 0.55, 0.9, '2 bar CO!D2!N atmosphere', /normal, charsize=1, color=0
;  xyouts, 0.55, 0.86, '250 K blackbody: 221.5 Wm!U-2!N', /normal, charsize=1, color=cib
;  xyouts, 0.55, 0.82, 'OLR 96.1319 Wm!U-2!N', /norma, charsize=1, color=cir


if  (plot_ps eq 1) then device, /close



end

