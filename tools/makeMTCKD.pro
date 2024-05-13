pro makeMTCKD
;------------------------------------------------
;Author: Wolf, E.T. 
;Created: sometime long ago                                     
;
;DESCRIPTION:
; Program to make input files containing MT_CKD H2O self 
;  continuum absorption coefficients for ExoRT.
;  Requires the users to have the aer mt_ckd continuum
;  model downloaded on their machine, with binary remnamed
;  "mtckdrun".
;   http://rtweb.aer.com/continuum_frame.html
;  on discover, run in /discover/nobackup/etwolf/models/mt_ckd_3.3/cntnm/rundir
;-----------------------------------------------





do_write = 1  ;if eq 1, then write the output file

do_band_mean    = 0
do_band_median  = 0
do_ksort_fit     = 1

outfile = "/discover/nobackup/etwolf/models/ExoRT/data/continuum/KH2O_MTCKD3.3_SELF.FRGN_n68_ngauss.nc"
;outfile = "/discover/nobackup/etwolf/models/ExoRT/data/continuum/KH2O_MTCKD3.3_SELF.FRGN_n68_avg.nc"

;-- choose one and only one --
;-- spectral resolution
do_n28 = 0
do_n35 = 0
do_n42 = 0
do_n68 = 1
do_n73 = 0
do_n84 = 0



;; 28 bin
if (do_n28 eq 1) then begin
  print, "using 28 spectral intervals"
  NU_LOW = [ 10,   350,  500,  630,  700,  820,  980,  1100, $
             1180, 1390, 1480, 1800, 2080, $
             2200, 2380, 2600,  3250,  4000,  4650,  5150, 6150, $
             7700, 8050, 12850, 16000, 22650, 29000, 38000 ]
  NU_HIGH = [ 350,  500,  630,  700,  820,  980,  1100, $
              1180, 1390, 1480, 1800, 2080, 2200, $
              2380, 2600,  3250,  4000,  4650,  5150, 6150, $
              7700, 8050, 12850, 16000, 22650, 29000, 38000, 50000 ]
endif


;;35 bin
if (do_n35 eq 1) then begin
  print, "using 35 spectral intervals"
  NU_LOW = [   10.0,          200.0,         350.0,        425.0,         500.0,   $
               630.0,         700.0,         820.0,        980.0,         1100.0,  $
              1180.0,        1390.0,        1480.0,       1800.0,        2080.0,  $
              2200.0,        2380.0,        2600.0,       3250.0,        4000.0,  $
              4650.0,        5150.0,        6150.0,       6650.0,        7700.0,  $
              8050.0,        9650.0,        11650.0,      12850.0,       14800.0, $
              16000.0,       19300.0,       22650.0,      29000.0,       38000.0 ]
  NU_HIGH = [  200.0,         350.0,        425.0,         500.0,   $
               630.0,         700.0,         820.0,        980.0,         1100.0,  $
              1180.0,        1390.0,        1480.0,       1800.0,        2080.0,  $
              2200.0,        2380.0,        2600.0,       3250.0,        4000.0,  $
              4650.0,        5150.0,        6150.0,       6650.0,        7700.0,  $
              8050.0,        9650.0,        11650.0,      12850.0,       14800.0, $
              16000.0,       19300.0,       22650.0,      29000.0,       38000.0, 50000.0 ]
endif


;42 bin
if (do_n42 eq 1) then begin
  print, "using 42 spectral intervals"
  NU_LOW = [     10.0,       200.0,         350.0,        425.0, $
                500.0,       630.0,         700.0,        820.0, $
                980.0,       1100.0,        1180.0,       1390.0, $
               1480.0,       1800.0,        2080.0,       2200.0, $
               2380.0,       2600.0,        3250.0,       4000.0, $
               4650.0,       4900.0,        5150.0,       5650.0, $
               6150.0,       6650.0,        7000.0,       7700.0, $
               8050.0,       9100.0,        10000.0,      11000.0, $
              11800.0,       12850.0,       13450.0,      14450.0, $
              15150.0,       16000.0,       19300.0,      22650.0, $
              29000.0,       38000.0 ]
  NU_HIGH = [   200.0,         350.0,        425.0, $
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


;; 68 bin
if (do_n68 eq 1) then begin
  print, "using 68 spectral intervals"
  NU_LOW =  [0.00E+00,   40.00000,   100.0000,   160.0000, $
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
  NU_HIGH =[40.00000,   100.0000,   160.0000, $
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


;; 73 bin 
if (do_n73 eq 1) then begin
  print, "using 73 spectral intervals"
  NU_LOW  = [ 0.00000, 40.0000, 100.000, 160.000, 220.000, 280.000, 330.000, $
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
  NU_HIGH= [ 40.0000, 100.000, 160.000, 220.000, 280.000, 330.000, $
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


;;84 bin (shortwave extension of n68)
if (do_n84 eq 1) then begin
  print, "using 84 spectral intervals (shortwave extension of n68)"
  NU_LOW =  [0.00E+00,   40.00000,   100.0000,   160.0000, $
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
             400000.00, 500000.00, 750000.00,  1000000.00 ]
  NU_HIGH = [40.00000,   100.0000,   160.0000, $
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

;; LMD co2h2ovar 63 bin
;fname="/Users/wolfe/Models/RT/RT_offline/absorption_data/co2_h2ovar/CO2_H2Ovar_MERGED.nc"
;ncid=ncdf_open(fname, /nowrite)
;ncdf_varget,ncid,'NU_LOW',NU_LOW
;ncdf_varget,ncid,'NU_HIGH',NU_HIGH
;ncdf_close, ncid


; end spectral interval defintions
;-------------
NU_MID = (NU_HIGH(*)+NU_LOW(*))/2.



PRESSURES = [1013.]
;TEMPERATURES = [296]
TEMPERATURES = [100, 110, 120, 130, 140, 150, 160, 170, 180, 190, 200, $
                     210, 220, 230, 240, 250, 260, 270, 280, 290, 300, $
                     310, 320, 330, 340, 350, 360, 370, 380, 390, 400, $
                     410, 420, 430, 440, 450, 460, 470, 480, 490, 500]

;TEMPERATURES = [ 100,150, 200, 250, 300, 350, 400, 450, 500]

;TEMPERATURES = [ 200, 250, 300, 350, 400]


H2O = [1.0]
CO2 = [0.0]
O2  = [0.0]

V1 = 0.0
V2 = 20000.
DV = 10.


nwvl = n_elements(NU_LOW)
np = n_elements(PRESSURES)
nt = n_elements(TEMPERATURES)
nw = n_elements(WEIGHTS)

pathlength = 1.0

inputfile="INPUT"
readfile="WATER.COEF"
header = strarr(101) 
;NDATA=1001   ;0 to 10,000 cm-1
NDATA=2001  ;0 to 20,000 cm-1
;NDATA=5001  ;0 to 50,000 cm-1
data = fltarr(5,NDATA)
wavenumber_vec = fltarr(NDATA)
kself_raw_vec = fltarr(NDATA)
kfrgn_raw_vec = fltarr(NDATA)
kco2cont_raw_vec = fltarr(NDATA)
kself = fltarr(nwvl)
kfrgn = fltarr(nwvl)
kco2cont = fltarr(nwvl)

kself_out = fltarr(nwvl,nt)
kfrgn_out = fltarr(nwvl,nt)
kco2cont_out = fltarr(nwvl,nt)

if (do_ksort_fit eq 1) then begin
  ng = 8
  npoints=10001
  npts_vec = findgen(npoints)
  interp_wav=fltarr(npoints)

  kself_ng_out = fltarr(ng,nwvl,nt)
  kfrgn_ng_out = fltarr(ng,nwvl,nt)
  
; 8 point RRTMG
  g_xpos_min = fltarr(8)
  g_xpos_min = [0.00000, 0.30192, 0.57571, 0.79583, $
                0.94178, 0.98890, 0.99576, 0.99939 ]

  g_xpos_max = fltarr(8)
  g_xpos_max = [0.30192, 0.57571, 0.79583, 0.94178, $
               0.98890, 0.99576, 0.99939, 1.00000 ]

  g_xpos = fltarr(8)
  g_xpos = (g_xpos_max + g_xpos_min)/2.0

  g_weights = [0.30192, 0.27379, 0.22012, 0.14595, $
               0.04712, 0.00686, 0.00363, 0.00061 ]
endif



for b=0,nt-1 do begin
  print, "TEMPERATURE:", TEMPERATURES(b)
     
  print, "----begin step ----"
  spawn, "date"
  ;--- create input file ---
  OPENW,lun,inputfile, /Get_Lun
  printf,lun, "Pressure (mb), Temperature (K), Path Length (cm),    VMR H2O,     VMRCO2,      VMRO2,  V1, V2, DV"
  printf,lun,format = '(f13.6, f17.4, f18.4, f12.8, f12.8, f12.8, f8.1, f8.1, f8.1)', $
              PRESSURES, TEMPERATURES(b), pathlength, H2O, CO2, O2, V1, V2, DV
  FREE_LUN,lun
  ;--- run ---
  spawn, "mtckdrun"
  ;--- read ---
  OPENR,lun,readfile
  READF,lun, header
  READF,lun, data
  FREE_LUN,lun
  wavenumber_vec(*) = data(0,*)
  kself_raw_vec(*) = data(3,*) 
  kfrgn_raw_vec(*) = data(4,*) 
  
  ; wavelength plotting
 ; if (b eq 0) then plot, 1.0e4/wavenumber_vec,kself_raw_vec, /ylog, xrange=[0,20000]
 ; if (b gt 0) then oplot,1.0e4/wavenumber_vec,kself_raw_vec

  ; wavenumber plotting
 ; if (b eq 0) then plot, wavenumber_vec,kself_raw_vec, /ylog, xrange=[0,20000]
 ; if (b gt 0) then oplot,wavenumber_vec,kself_raw_vec, linestyle=0, color=255
;  if (b gt 0) then oplot,wavenumber_vec,kfrgn_raw_vec, color=200, linestyle=0
   
  ct=0
  for d=0,nwvl-1 do begin
    print, "-------------------------------------------------------------"
    print, "BAND: ", d
    print, "LOW/HI: ",NU_LOW(d), NU_HIGH(d)
    x=where(wavenumber_vec le NU_HIGH(d) and wavenumber_vec ge NU_LOW(d),num)
;       print, wavenumber_vec(ct), wavenumber_vec(ct+num-1)
    if (num gt 0) then begin
 ;     oplot, wavenumber_vec(x), kself_raw_vec(x), psym=4
;      print, "wavenumbers ", wavenumber_vec(x)
;      print, "kself_raw ", kself_raw_vec(x)      
;      print, "min, avg, max ", min(kself_raw_vec(x)), total(kself_raw_vec(x))/num, max(kself_raw_vec(x))

      if (do_ksort_fit) then begin
      ;-------- K-sort ------
        temp_self=kself_raw_vec(x)
        temp_frgn=kfrgn_raw_vec(x)

;        print, "temp(sort(temp))", temp(sort(temp))

;        plot, wavenumber_vec(x), temp, /ylog
;        oplot, wavenumber_vec(x), temp, psym=4
;stop
;        plot, wavenumber_vec(x), temp(sort(temp)), /ylog
;        oplot, wavenumber_vec(x), temp(sort(temp)), psym=4
;stop
;        print, num, wavenumber_vec(ct),wavenumber_vec(ct+num)
;        oplot, [wavenumber_vec(ct),wavenumber_vec(ct+num)], [total(kself_raw_vec(x))/num,total(kself_raw_vec(x))/num]
;        oplot, [wavenumber_vec(ct),wavenumber_vec(ct+num)], [MEDIAN(kself_raw_vec(x)), MEDIAN(kself_raw_vec(x))], linestyle=2

        delta=(wavenumber_vec(ct+num-1)-wavenumber_vec(ct))/npoints
        interp_wav(*) = 0.0
        FOR nkb=0, npoints-1 DO interp_wav(nkb) = min(wavenumber_vec(x)) + delta*nkb

        gfunc_self = interpol(temp_self(sort(temp_self)), wavenumber_vec(x), interp_wav, /lsquadratic)
        gfunc_frgn = interpol(temp_frgn(sort(temp_frgn)), wavenumber_vec(x), interp_wav, /lsquadratic)
;print, "gfunc", min(gfunc), max(gfunc)
;stop
;        oplot, interp_wav, gfunc, linestyle=0, color=250, thick=4
;stop
        ;print, g_xpos*npoints, round(g_xpos*npoints), interp_wav(round(g_xpos*npoints))
print, d,"Gauss points SELF:", gfunc_self(round(g_xpos*npoints))
print, d,"Gauss points FOREIGN:", gfunc_frgn(round(g_xpos*npoints))
;        plot, npts_vec/npoints, gfunc, linestyle=0, thick=2, xrange=[0.0,1], /ylog
;        oplot, npts_vec(round(g_xpos*npoints))/npoints, gfunc(round(g_xpos*npoints)), psym=4, thick=5, color=250
        ct=ct+(num-1)

        kself_ng_out(*,d,b) =  gfunc_self(round(g_xpos*npoints))
        kfrgn_ng_out(*,d,b) =  gfunc_frgn(round(g_xpos*npoints))
 ;       stop
       ;\------- K-sort ------
     endif else begin
       if (do_band_mean) then begin
         ; use average of all values within band
         kself(d) = total(kself_raw_vec(x))/num 
         kfrgn(d) = total(kfrgn_raw_vec(x))/num 
       endif
       if (do_band_median) then begin
         ; use median of all values within band
         kself(d) = MEDIAN(kself_raw_vec(x))
         kfrgn(d) = MEDIAN(kfrgn_raw_vec(x))
       endif
       ; finalize outputs for netcdf
       kself_out(d,b) = kself(d)
       kfrgn_out(d,b) = kfrgn(d)
  ;     print,"self",  d, 1.0e4/NU_MID(d),total(kself_raw_vec(x))/num, MEDIAN(kself_raw_vec(x)), 10.^(mean(temps))
  ;     print,"frgn",  d, 1.0e4/NU_MID(d),total(kfrgn_raw_vec(x))/num, MEDIAN(kfrgn_raw_vec(x)), 10.^(mean(tempf))
      endelse
    endif else begin
; set to zero if no continuum found in band
      kself_out(d,b) = 0.0
      kfrgn_out(d,b) = 0.0
;      kco2cont_out(d,b) = 0.0
    endelse
  endfor
  
  ;wavelength plotting
  ;oplot, 1.0e4/NU_MID,kself_out(*,b), psym=2, color=200
  ;oplot, 1.0e4/NU_MID,kself, psym=2

  ;wavenumber plotting
  oplot, NU_MID,kself_out(*,b), psym=2, color=255
;  oplot, NU_MID,kfrgn_out(*,b), psym=2, color=200
;  oplot, NU_MID,kself, psym=2

;  print, "h2O self: ", kself_out
;  print, "co2comt: ",kco2cont
  print, "---- end step ----"
  spawn, "date"
;stop
endfor

;stop  ; to see plot





;--- write to netcdf ---
if (do_write eq 1) then begin

id = NCDF_CREATE(outfile, /CLOBBER)
dim1 = NCDF_DIMDEF(id, 'nbands',nwvl)
dim3 = NCDF_DIMDEF(id, 'ntemp',nt)
if (do_ksort_fit eq 1) then dim4 = NCDF_DIMDEF(id, 'ngauss',ng)


varid1 = NCDF_VARDEF(id, 'NU_LOW', dim1, /float)
varid2 = NCDF_VARDEF(id, 'NU_HIGH', dim1, /float)
varid4 = NCDF_VARDEF(id, 'TEMPERATURES', dim3, /float)

if (do_ksort_fit eq 1) then begin
  varid6 = NCDF_VARDEF(id, 'KSELF',[dim4, dim1, dim3], /float)
  varid7 = NCDF_VARDEF(id, 'KFRGN',[dim4, dim1, dim3], /float)

endif else begin
  varid6 = NCDF_VARDEF(id, 'KSELF',[dim1,dim3], /float)
  varid7 = NCDF_VARDEF(id, 'KFRGN',[dim1,dim3], /float)
endelse

NCDF_ATTPUT, id, varid1, "title", "Wavenumber grid low"
NCDF_ATTPUT, id, varid1, "units", "cm-1"
NCDF_ATTPUT, id, varid2, "title", "Wavenumber grid high"
NCDF_ATTPUT, id, varid2, "units", "cm-1"
NCDF_ATTPUT, id, varid4, "title", "temperature grid"
NCDF_ATTPUT, id, varid4, "units", "K"
NCDF_ATTPUT, id, varid6, "title", "h2o self continuum coefficients"
NCDF_ATTPUT, id, varid6, "units", "cm2 / molecule"
NCDF_ATTPUT, id, varid7, "title", "h2o foreign continuum coefficients"
NCDF_ATTPUT, id, varid7, "units", "cm2 / molecule"

NCDF_CONTROL,id, /ENDEF

NCDF_VARPUT, id, varid1, NU_LOW
NCDF_VARPUT, id, varid2, NU_HIGH
NCDF_VARPUT, id, varid4, TEMPERATURES
if (do_ksort_fit eq 1) then begin
  NCDF_VARPUT, id, varid6, kself_ng_out
  NCDF_VARPUT, id, varid7, kfrgn_ng_out
endif else begin
  NCDF_VARPUT, id, varid6, kself_out
  NCDF_VARPUT, id, varid7, kfrgn_out
endelse

NCDF_CLOSE, id

endif

end
