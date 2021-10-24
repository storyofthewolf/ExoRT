pro makeStellarSpectrum_blackbody
;------------------------------------------------
;Author: Wolf, E.T.
;Created: October 2021
;Revision history: adapted from old code
;
;DESCRIPTION:
; Creates solar spectral files for ExoRT.
;  Creates blackbody spectrum using a set temperature
;  This is a separate script because I was too lazy to
;  merge this script with makeStellarSpectrum_fromSED.pro
;  to create a single script
;
;
;------------------------------------------------

;-------------------------
;-- stellar temperature --
;-------------------------
temperature = 3000 ;[K]
;-------------------------

print, "Making blackbody spectrum: ", temperature, " K"

radTEMP = temperature
outname = "blackbody_"+string(radTEMP)+"K_n68.nc"
outname=STRJOIN(STRSPLIT(outname, /EXTRACT), '')

;-- choose one and only one --
;-- spectral resolution --
do_n28 = 0
do_n35 = 0
do_n42 = 0
do_n68 = 1
do_n73 = 0
do_n84 = 0

do_plot = 1              ; plot blackbody curve and binned SED
do_write_netcdf = 0    ; write netcdf SED file for ExoRT
do_write_bb_dat = 1   ; write text data file of blackbody curve


;; -- supported spectral resolutions ---
;=  rt = 28
if (do_n28 eq 1) then begin
  nrtwavl = 28
  rtwavlow = dblarr(nrtwavl)
  rtwavhi = dblarr(nrtwavl)
  rtwavlow = [10.,   350.,  500.,  630.,  700.,   820.,   980.,   1100., 1180., 1390., $
              1480., 1800., 2080., 2200., 2380.,  2600.,  3250.,  4000., 4650., $
              5150., 6150., 7700., 8050., 12850., 16000., 22650., 29000., $
              38000. ]
  rtwavhi  = [350.,   500.,  630.,  700.,  820.,   980.,   1100.,  1180., 1390., $
              1480.,  1800., 2080., 2200., 2380.,  2600.,  3250.,  4000., 4650., $
              5150.,  6150., 7700., 8050., 12850., 16000., 22650., 29000., $
              38000., 50000. ]
endif

;=  rt = 42 
if (do_n42 eq 1) then begin
  nrtwavl = 42
  rtwavlow = dblarr(nrtwavl)
  rtwavhi = dblarr(nrtwavl)
  rtwavlow =  [  10.0,       200.0,         350.0,        425.0, $
                 500.0,      630.0,         700.0,        820.0, $
                 980.0,      1100.0,        1180.0,       1390.0, $
                 1480.0,     1800.0,        2080.0,       2200.0, $
                 2380.0,     2600.0,        3250.0,       4000.0, $
                 4650.0,     4900.0,        5150.0,       5650.0, $
                 6150.0,     6650.0,        7000.0,       7700.0, $
                 8050.0,     9100.0,        10000.0,      11000.0, $
                 11800.0,    12850.0,       13450.0,      14450.0, $
                 15150.0,    16000.0,       19300.0,      22650.0, $
                 29000.0,    38000.0 ]
  rtwavhi =   [  200.0,      350.0,        425.0, $
                 500.0,      630.0,        700.0,        820.0, $
                 980.0,      1100.0,       1180.0,       1390.0, $
                 1480.0,     1800.0,       2080.0,       2200.0, $
                 2380.0,     2600.0,       3250.0,       4000.0, $
                 4650.0,     4900.0,       5150.0,       5650.0, $
                 6150.0,     6650.0,       7000.0,       7700.0, $
                 8050.0,     9100.0,       10000.0,      11000.0, $
                 11800.0,    12850.0,      13450.0,      14450.0, $
                 15150.0,    16000.0,      19300.0,      22650.0, $
                 29000.0,    38000.0,      50000.0 ]
endif


;=  rt = 68
if (do_n68 eq 1) then begin 
  nrtwavl = 68
  rtwavlow = dblarr(nrtwavl)
  rtwavhi = dblarr(nrtwavl) 
  rtwavlow =  [1.000,   40.00000,   100.0000,   160.0000, $
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
            30376.00,   32562.00,   35087.00,   36363.00 ]
  rtwavhi =   [ 40.00000,   100.0000,   160.0000, $
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
                42087.00 ]
endif

; = rt =  73 bin
if (do_n73 eq 1) then begin
  nrtwavl = 73
  rtwavlow = dblarr(nrtwavl)
  rtwavhi = dblarr(nrtwavl) 
  rtwavlow  = [ 1.00000, 40.0000, 100.000, 160.000, 220.000, 280.000, 330.000, $
                380.000, 440.000, 495.000, 545.000, 617.000, 667.000, 720.000, $
                800.000, 875.000, 940.000, 1000.00, 1065.00, 1108.00, 1200.00, $
                1275.00, 1350.00, 1450.00, 1550.00, 1650.00, 1750.00, 1850.00, $
                1950.00, 2050.00, 2200.00, 2439.02, 2564.10, 2777.78, 3174.60, $
                3508.77, 3773.59, 4081.63, 4545.46, 4716.98, 5154.64, 5376.34, $
                5555.56, 5952.38, 6172.84, 6578.95, 6711.41, 6849.31, 7042.25, $
                7462.69, 7692.31, 8064.52, 8333.33, 8620.69, 8928.57, 9090.91, $
                9259.26, 9708.74, 10869.6, 11111.1, 11363.6, 11494.3, 12500.0, $
                12820.5, 14492.8, 16393.4, 18181.8, 20000.0, 22222.2, 23809.5, $
                25974.0, 28985.5, 33333.3 ]
  rtwavhi  = [ 40.0000, 100.000, 160.000, 220.000, 280.000, 330.000, $
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

; = rt =  84 bin
if (do_n84 eq 1) then begin
  nrtwavl = 84
  rtwavlow = dblarr(nrtwavl)
  rtwavhi = dblarr(nrtwavl) 
  rtwavlow =  [1.000,   40.00000,   100.0000,   160.0000, $
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
            150000.00, 175000.00, 200000.00,   300000.00, $
            400000.00, 500000.00, 750000.00,  1000000.00 ]
rtwavhi =   [ 40.00000,   100.0000,   160.0000, $
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
              400000.00, 500000.00,  750000.00,  1000000.00, $
             1250000.00 ]
endif

;convert to micronds
;hi/low are crossed intentionally
rtwvlhi = 1.0e4/rtwavlow(*)  
rtwvllow = 1.0e4/rtwavhi(*)  
rtwvlmid = (rtwvlhi+rtwvllow)/2.0
rtwvldel = (1.0e4/rtwavlow(*) - 1.0e4/rtwavhi(*))

SunSolidAngle = 6.87e-5 ; as the Sun appears from Earth today
Snorm = 1360.
S0 = Snorm
Rsun = 6.957e8

npts=60000.
wnm=dblarr(npts)
wlgth=dblarr(npts)
dwm =1.0
sunm=dblarr(npts)

for i=0, npts-1 do begin
  wnm(i) = 1.+i*dwm
  wlgth(i) = 1.0e4/wnm(i)
endfor

wm1=wnm(0)
wm2=wnm(npts-1)

; Blackbody computation
c1 = 1.1911E8 ; [W m^-2 sr-1 um-4]            
c2 = 1.4388E4 ; [K um]
for i=0,npts-1 do begin
  Blambda1 = c1/(wlgth(i)^5*(exp(c2/wlgth(i)/temperature)-1))
  sunm(i)=Blambda1
endfor

solarflux = dblarr(nrtwavl)

dwn = 1.0
dwni = 1
for iw=0, nrtwavl-1 do begin
  ftmp=0.
  wn_lw=rtwavlow(iw)
  wn_hi=rtwavhi(iw)
  for w=wn_lw,wn_hi-dwn,dwni do begin   ; integrate over each  spectral interval in steps of 10 cm-1  
    dwl=1.0e4/w - 1.0e4/(w+dwn)
    ftmp=ftmp+interpol(sunm(*),wnm,w, /spline)*dwl
 endfor
 solarflux(iw)=ftmp
endfor

ScaleFac = Snorm/total(solarflux(*))
solarflux(*) = solarflux(*)*ScaleFac

if (do_plot eq 1) then begin
  print, "plotting blackbody_spectrum.eps"
  loadct,40
  !P.font=0
  set_plot,'PS'
  device,xsize=11.0,ysize=7.0,xoff=1.0,yoff=3.0,/cm
  device,file='blackbody_spectrum.eps'
  device,/color,BITS=8,/ENCAPSULATED, /CMYK   
  device, set_font='Helvetica-Oblique', FONT_INDEX=20
  device, set_font='Helvetica-Bold', FONT_INDEX=19
  device, set_font='helvetica',FONT_INDEX=18

  ;create artifical bar
  q=0
  xbar1 = fltarr(nrtwavl*2)
  ybar1 = fltarr(nrtwavl*2)
  for i=0,nrtwavl-1 do begin
    xbar1(q) = 1.0e4/rtwavlow(i)        
    xbar1(q+1) = 1.0e4/rtwavhi(i)
    ybar1(q) = solarflux(i)/rtwvldel(i)
    ybar1(q+1) = solarflux(i)/rtwvldel(i)
    q=q+2
  endfor

  lthk=3.0
  lthk2=6.0

  loadct,40

  plot, wlgth, sunm(*), xtitle="!18wavenlength (microns)", xrange=[0.0,5.0], xstyle=1, yrange=[0,1200], $
                          xthick=3.0, ythick=3.0, /nodata, ytitle="Radiance (W m!U-2!N micron!U-1!N)", charsize=0.7
  oplot, wlgth, sunm(*)*ScaleFac, linestyle=0, thick=lthk, color=60
  oplot, xbar1, ybar1, linestyle=0, thick=lthk2, color=230

  device, /close
  set_plot,'X'

endif

if (do_write_bb_dat eq 1) then begin

  radTEMP = temperature
  outname = "blackbody_"+string(radTEMP)+"K.dat"
  outname=STRJOIN(STRSPLIT(outname, /EXTRACT), '')
  print, "writing blackbody as text file output, ", outname
  OPENW,lun,outname, /Get_Lun
  printf,lun, format = '("wavelength(um)   radiance(W/m2/um)")'
  for bb=0, npts-1 do begin
    printf,lun, format = '(2F15.8)', wlgth(bb), sunm(bb)*ScaleFac
  endfor
endif


if (do_write_netcdf eq 1) then begin
  print, "writing output to ", outname
  solarflux_out = solarflux
  ;; write output file
  one=1
  id = NCDF_CREATE(outname, /CLOBBER)
  dim1 = NCDF_DIMDEF(id,'nw',nrtwavl)
  dim2 = NCDF_DIMDEF(id,'one',one)

  varid1 = NCDF_VARDEF(id, 'wav_low',[dim1], /float)
  varid2 = NCDF_VARDEF(id, 'wav_high', [dim1], /float)
  varid3 = NCDF_VARDEF(id, 'S0', [dim2], /float)
  varid4 = NCDF_VARDEF(id, 'solarflux', [dim1], /float)
  NCDF_ATTPUT, id, varid1, "title", "wavenumber edge, low"
  NCDF_ATTPUT, id, varid1, "units", "cm -1"
  NCDF_ATTPUT, id, varid2, "title", "wavenumber edge, high"
  NCDF_ATTPUT, id, varid2, "units", "cm -1"
  NCDF_ATTPUT, id, varid3, "title", "solar constant"
  NCDF_ATTPUT, id, varid3, "units", "W/m2"
  NCDF_ATTPUT, id, varid4, "title", "spectrally integrated solar flux"
  NCDF_ATTPUT, id, varid4, "units", "W/m2"
  NCDF_CONTROL, id, /ENDEF

  NCDF_VARPUT, id, varid1, rtwavlow
  NCDF_VARPUT, id, varid2, rtwavhi
  NCDF_VARPUT, id, varid3, S0
  NCDF_VARPUT, id, varid4, solarflux_out

  NCDF_CLOSE, id

endif

end
