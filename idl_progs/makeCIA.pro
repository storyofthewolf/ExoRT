pro makeCIA
;=====================================================
;Purpose :: construct CIA netcdf files for input into 
; ExoRT using *.cia files downloaded from HITRAN,
; Author: Wolf, E.T.

do_write = 0

do_n2n2 = 0
do_h2h2 = 0
do_n2h2 = 0
do_co2co2 = 1
do_co2ch4 = 0
do_co2h2 = 0

do_h2oh2o = 0 ; experimental
do_h2on2 = 0  ; experimental

outfile_n2n2 = "N2-N2_cia_n73.nc"
outfile_h2h2 = "H2-H2_cia_n73.nc"
outfile_n2h2 = "N2-H2_cia_n73.nc"

outfile_co2co2_lw= "CO2-CO2_cia_lw_n73.nc"
outfile_co2co2_sw= "CO2-CO2_cia_sw_n73.nc"
outfile_co2ch4 = "CO2-CH4_cia_n73.nc"
outfile_co2h2 = "CO2-H2_cia_n73.nc"

outfile_h2oh2o = "H2O-H2O_cia_68bin.nc"  ; experimental
outfile_h2on2 = "H2O-N2_cia_68bin.nc"    ; experimental

loschmidt = 2.6867774e19 ;molecule cm-3 at STP

;; 28 bin
;NU_LOW = [ 10,   350,  500,  630,  700,  820,  980,  1100, $
;               1180, 1390, 1480, 1800, 2080, $
;               2200, 2380, 2600,  3250,  4000,  4650,  5150, 6150, $
;               7700, 8050, 12850, 16000, 22650, 29000, 38000 ]

;NU_HIGH = [ 350,  500,  630,  700,  820,  980,  1100, $
;                 1180, 1390, 1480, 1800, 2080, 2200, $
;                 2380, 2600,  3250,  4000,  4650,  5150, 6150, $
;                 7700, 8050, 12850, 16000, 22650, 29000, 38000, 50000 ]


;; 68 bin
;NU_LOW =  [0.00E+00,   40.00000,   100.0000,   160.0000, $
;           220.0000,   280.0000,   330.0000,   380.0000, $
;           440.0000,   495.0000,   545.0000,   617.0000, $
;           667.0000,   720.0000,   800.0000,   875.0000, $
;           940.0000,   1000.000,   1065.000,   1108.000, $
;           1200.000,   1275.000,   1350.000,   1450.000, $
;           1550.000,   1650.000,   1750.000,   1850.000, $
;           1950.000,   2050.000,   2200.000,   2397.000, $
;           2494.000,   2796.000,   3087.000,   3425.000, $
;           3760.000,   4030.000,   4540.000,   4950.000, $
;           5370.000,   5925.000,   6390.000,   6990.000, $
;           7650.000,   8315.000,   8850.000,   9350.000, $
;           9650.000,   10400.00,   11220.00,   11870.00, $
;           12790.00,   13300.00,   14470.00,   15000.00, $
;           16000.00,   16528.00,   17649.00,   18198.00, $
;           18518.00,   22222.00,   25641.00,   29308.00, $
;           30376.00,   32562.00,   35087.00,   36363.00 ]

;NU_HIGH =   [ 40.00000,   100.0000,   160.0000, $
;           220.0000,   280.0000,   330.0000,   380.0000, $
;           440.0000,   495.0000,   545.0000,   617.0000, $
;           667.0000,   720.0000,   800.0000,   875.0000, $
;           940.0000,   1000.000,   1065.000,   1108.000, $
;           1200.000,   1275.000,   1350.000,   1450.000, $
;           1550.000,   1650.000,   1750.000,   1850.000, $
;           1950.000,   2050.000,   2200.000,   2397.000, $
;           2494.000,   2796.000,   3087.000,   3425.000, $
;           3760.000,   4030.000,   4540.000,   4950.000, $
;           5370.000,   5925.000,   6390.000,   6990.000, $
;           7650.000,   8315.000,   8850.000,   9350.000, $
;           9650.000,   10400.00,   11220.00,   11870.00, $
;           12790.00,   13300.00,   14470.00,   15000.00, $
;           16000.00,   16528.00,   17649.00,   18198.00, $
;           18518.00,   22222.00,   25641.00,   29308.00, $
;           30376.00,   32562.00,   35087.00,   36363.00, $
;           42087.00 ]

;42 bin
;NU_LOW = [      10.0,       200.0,         350.0,        425.0, $
;                500.0,       630.0,         700.0,        820.0, $
;                980.0,       1100.0,        1180.0,       1390.0, $
;               1480.0,       1800.0,        2080.0,       2200.0, $
;               2380.0,       2600.0,        3250.0,       4000.0, $
;               4650.0,       4900.0,        5150.0,       5650.0, $
;               6150.0,       6650.0,        7000.0,       7700.0, $
;               8050.0,       9100.0,        10000.0,      11000.0, $
;              11800.0,       12850.0,       13450.0,      14450.0, $
;              15150.0,       16000.0,       19300.0,      22650.0, $
;              29000.0,       38000.0 ]

;NU_HIGH = [      200.0,         350.0,        425.0, $
;                500.0,       630.0,         700.0,        820.0, $
;                980.0,       1100.0,        1180.0,       1390.0, $
;               1480.0,       1800.0,        2080.0,       2200.0, $
;               2380.0,       2600.0,        3250.0,       4000.0, $
;               4650.0,       4900.0,        5150.0,       5650.0, $
;               6150.0,       6650.0,        7000.0,       7700.0, $
;               8050.0,       9100.0,        10000.0,      11000.0, $
;              11800.0,       12850.0,       13450.0,      14450.0, $
;              15150.0,       16000.0,       19300.0,      22650.0, $
;              29000.0,       38000.0,       50000.0 ]

;;35 bin
;NU_LOW = [   10.0,          200.0,         350.0,        425.0,         500.0,   $
;               630.0,         700.0,         820.0,        980.0,         1100.0,  $
;              1180.0,        1390.0,        1480.0,       1800.0,        2080.0,  $
;              2200.0,        2380.0,        2600.0,       3250.0,        4000.0,  $
;              4650.0,        5150.0,        6150.0,       6650.0,        7700.0,  $
;              8050.0,        9650.0,        11650.0,      12850.0,       14800.0, $
;              16000.0,       19300.0,       22650.0,      29000.0,       38000.0 ]

;NU_HIGH = [    200.0,         350.0,        425.0,         500.0,   $
;               630.0,         700.0,         820.0,        980.0,         1100.0,  $
;              1180.0,        1390.0,        1480.0,       1800.0,        2080.0,  $
;              2200.0,        2380.0,        2600.0,       3250.0,        4000.0,  $
;              4650.0,        5150.0,        6150.0,       6650.0,        7700.0,  $
;              8050.0,        9650.0,        11650.0,      12850.0,       14800.0, $
;              16000.0,       19300.0,       22650.0,      29000.0,       38000.0, 50000.0 ]




; 73 bin
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

NU_HIGH  = [ 40.0000, 100.000, 160.000, 220.000, 280.000, 330.000, $
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




;; LMD co2h2ovar 63 bin
;fname="/Users/wolfe/Models/RT/RT_offline/absorption_data/co2_h2ovar/CO2_H2Ovar_MERGED.nc"
;ncid=ncdf_open(fname, /nowrite)
;ncdf_varget,ncid,'NU_LOW',NU_LOW
;ncdf_varget,ncid,'NU_HIGH',NU_HIGH
;ncdf_close, ncid






;-------------------------------------
;----- n2-n2 cia ---------------------
;-------------------------------------
if (do_n2n2 eq 1 ) then begin
print, "====== calculating N2-N2 CIA ======"

;;;  --- read n2-n2 cia cia file  ---- ;;;;
inputfile = "/gpfsm/dnb53/etwolf/models/ExoRT/data/cia/orig/N2-N2_2011.cia"
OPENR,lun,inputfile, /Get_Lun

chemsym = ' '
; read first block
; spectral range    :: 0.020   554.000  cm-1
; temperature range :: 40, 51.7, 66.7, 86.2, 11.3, 143.8, 185.7, 239.8, 309.7, 400.0
data1 = dblarr(10,2,582)  & data1t = dblarr(2,582)
temp1 = dblarr(10)
for t=0, 9 do begin
  readf,lun,format = '(A20, 2F10, I7, F7, E10, F10)', $
                     chemsym, wmax, wmin, npts, temp, maxcia, res
  readf,lun, data1t
  temp1(t) = temp
  data1(t,*,*) = data1t(*,*)
endfor

; read second block
; spectral range    :: 1850.004  3000.094  cm-1
; temperature range :: 300.9, 323.6, 343.5, 355.3, 362.5
data2 = dblarr(5,2,9543)  & data2t = dblarr(2,9543)
temp2 = dblarr(5)
for t=0, 4 do begin
  readf,lun,format = '(A20, 2F10, I7, F7, E10, F10)', $
                     chemsym, wmax, wmin, npts, temp, maxcia, res
  readf,lun, data2t
  temp2(t) = temp
  data2(t,*,*) = data2t(*,*)
endfor

; read third block
; spectral range    :: 1999.900 2697.900 cm-1
; temperature range :: 228.2, 233.7, 243.2, 253.1, 272.1
data3 = dblarr(5,3,2806)  & data3t = dblarr(3,2806)
data3t_2 = dblarr(2,2806)
temp3 = dblarr(5)
for t=0, 1 do begin
  readf,lun,format = '(A20, 2F10, I7, F7, E10, F10)', $
                     chemsym, wmax, wmin, npts, temp, maxcia, res
  readf,lun, data3t
  temp3(t) = temp
  data3(t,*,*) = data3t(*,*)
endfor
for t=2, 2 do begin
  readf,lun,format = '(A20, 2F10, I7, F7, E10, F10)', $
                     chemsym, wmax, wmin, npts, temp, maxcia, res
  readf,lun, data3t_2
  temp3(t) = temp
  data3(t,0:1,*) = data3t_2(*,*)
endfor
for t=3, 4 do begin
  readf,lun,format = '(A20, 2F10, I7, F7, E10, F10)', $
                     chemsym, wmax, wmin, npts, temp, maxcia, res
  readf,lun, data3t
  temp3(t) = temp
  data3(t,*,*) = data3t(*,*)
endfor
FREE_LUN,lun

;set_plot,'PS'
;loadct,40
;device,xsize=7.0,ysize=4,xoff=0.5,yoff=5.5, /inches, /encapsulated
;device,/color,bits=8
;device, file="n2-n2.ps"

;plot block 1
for i=0,9 do begin
  wave = data1(i,0,*)
  abs = data1(i,1,*)
  abs = loschmidt*abs*loschmidt
  if i eq 0 then plot, wave, abs, /ylog, xrange=[0, 3100], /nodata, xstyle=1
  oplot, wave, abs, color=175
endfor

;plot block 2
for i=0,4 do begin
  wave = data2(i,0,*)
  abs = data2(i,1,*)
;  if i eq 0 then plot, wave, abs
  abs = loschmidt*abs*loschmidt
  oplot, wave, abs, color=150
endfor

;plot block 3
for i=0,4 do begin
  wave = data3(i,0,*)
  abs = data3(i,1,*)
;  if i eq 0 then plot, wave, abs
  abs = loschmidt*abs*loschmidt
  oplot, wave, abs, color=255
endfor

;device, /close
;set_plot,'X'

; interp
nwl = n_elements(nu_low)
k1 = fltarr(nwl, 10)
k2 = fltarr(nwl, 5)
k3 = fltarr(nwl, 5)
for d=0, nwl-1 do begin

   for t=0, n_elements(temp1)-1 do begin
     wavenumber_vec = data1(t,0,*)  
     abs_vec = loschmidt*data1(t,1,*)*loschmidt
  ;   abs_vec = data1(t,1,*)
     x1=where(wavenumber_vec le NU_HIGH(d) and wavenumber_vec ge NU_LOW(d),num)
     if (num gt 0) then begin
;       k1(d,t) = total(abs_vec(x1))/num
       k1(d,t) = MEDIAN(abs_vec(x1))
       if (k1(d,t) lt 0.0) then k1(d,t) = 0.0
     endif else begin
       k1(d,t) = 0.0
     endelse
   endfor

   for t=0, n_elements(temp2)-1 do begin
     wavenumber_vec = data2(t,0,*)  
     abs_vec = loschmidt*data2(t,1,*)*loschmidt
 ;    abs_vec = data2(t,1,*)
     x2=where(wavenumber_vec le NU_HIGH(d) and wavenumber_vec ge NU_LOW(d),num)
     if (num gt 0) then begin
;       k2(d,t) = total(abs_vec(x2))/num
       k2(d,t) = MEDIAN(abs_vec(x2))
       if (k2(d,t) lt 0.0) then k2(d,t) = 0.0
     endif else begin
       k2(d,t) = 0.0
     endelse
   endfor


   for t=0, n_elements(temp3)-1 do begin
     wavenumber_vec = data3(t,0,*)  
     abs_vec = loschmidt*data3(t,1,*)*loschmidt
;     abs_vec = data3(t,1,*)
     x3=where(wavenumber_vec le NU_HIGH(d) and wavenumber_vec ge NU_LOW(d),num)
     if (num gt 0) then begin
;       k3(d,t) = total(abs_vec(x3))/num
       k3(d,t) = MEDIAN(abs_vec(x3))
       if (k3(d,t) lt 0.0) then k3(d,t) = 0.0
     endif else begin
       k3(d,t) = 0.0
     endelse
   endfor

endfor

;
for h=0, n_elements(temp1)-1 do begin
   oplot, (NU_LOW(*) + nu_HIGH(*))/2., k1(*,h), color=200, thick=4, psym=4, symsize=0.9
endfor
for h=0, n_elements(temp2)-1 do begin
   oplot, (NU_LOW(*) + nu_HIGH(*))/2., k2(*,h), color=200, thick=4, psym=4, symsize=0.9
endfor
for h=0, n_elements(temp3)-1 do begin
   oplot, (NU_LOW(*) + nu_HIGH(*))/2., k3(*,h), color=200, thick=4, psym=4, symsize=0.9
endfor

;create interp array
karr = fltarr(nwl, 10)
kmaster_temp = fltarr(nwl, 10)
kmaster = fltarr(nwl, 10)

; combine blocks 2 and 3, since they cover the same approximate range
tarr = fltarr(10)
for k=0, 4 do begin
  karr(*,k) = k3(*,k)
  tarr(k) = temp3(k)
endfor
for k=5, 9 do begin
  karr(*,k) = k2(*,k-5)
  tarr(k) = temp2(k-5)
endfor
print, "temp1", temp1 
print, "temp2", temp2 
print, "temp3", temp3 
print, "tarr", tarr

temp_master = temp1
print, "temp_master", temp_master
;stop
for d=0,nwl-1 do begin
    kmaster_temp(d,*) = interpol(karr(d,*), tarr ,temp_master) ; + ktemp1(d,*)
endfor


for t=0, 9 do begin
  ; if temperature grid point is less than block 2 and 3 combined array
  if (temp_master(t) lt tarr(0)) then kmaster(*,t) = k1(*,t) + k3(*,0)  
  if (temp_master(t) gt tarr(9)) then kmaster(*,t) = k1(*,t) + k2(*,4)
  if (temp_master(t) gt tarr(0) and (temp_master(t) lt tarr(9))) then kmaster(*,t) = k1(*,t) + kmaster_temp(*,t)
endfor





;for h=0,4 do begin
;   oplot, (nu_low(*) + nu_high(*))/2., ktemp2(*,h), color=200, thick=4, psym=4
;   oplot, (nu_low(*) + nu_high(*))/2., ktemp3(*,h), color=200, thick=4, psym=4
;endfor
for h=0, 9 do begin
   oplot, (nu_low(*) + nu_high(*))/2., kmaster(*,h), psym=4, thick=4
;   oplot, (nu_low(*) + nu_high(*))/2., ktemp1(*,h), color=200, thick=4, psym=4
  print, temp_master(h), " / ",kmaster(*,h)
endfor



;print, temp_master
;print, kmaster
;print, "--"
;print, temp1
;print, ktemp1
;print, temp2
;print, ktemp2
;print, temp3
;print, ktemp3
;print, tarr
;stop

nwvl = nwl
nt = n_elements(temp_master)
TEMPERATURES = temp_master
k_out = kmaster

;--- write to netcdf ---
if (do_write eq 1) then begin
outfile = outfile_n2n2 
id = NCDF_CREATE(outfile, /CLOBBER)
dim1 = NCDF_DIMDEF(id, 'wbins',nwvl)
dim3 = NCDF_DIMDEF(id, 'tbins',nt)

varid1 = NCDF_VARDEF(id, 'nu_low', dim1, /float)
varid2 = NCDF_VARDEF(id, 'nu_high', dim1, /float)
varid4 = NCDF_VARDEF(id, 'temp', dim3, /float)
varid6 = NCDF_VARDEF(id, 'sigma',[dim1,dim3], /float)


NCDF_ATTPUT, id, varid1, "title", "Wavenumber grid low"
NCDF_ATTPUT, id, varid1, "units", "cm-1"
NCDF_ATTPUT, id, varid2, "title", "Wavenumber grid high"
NCDF_ATTPUT, id, varid2, "units", "cm-1"
NCDF_ATTPUT, id, varid4, "title", "temperature grid"
NCDF_ATTPUT, id, varid4, "units", "K"
NCDF_ATTPUT, id, varid6, "title", "band averaged CIA absorption"
NCDF_ATTPUT, id, varid6, "units", "cm^-1 amagat^-2"

NCDF_CONTROL,id, /ENDEF

NCDF_VARPUT, id, varid1, NU_LOW
NCDF_VARPUT, id, varid2, NU_HIGH
NCDF_VARPUT, id, varid4, TEMPERATURES
NCDF_VARPUT, id, varid6, k_out
NCDF_CLOSE, id

endif

endif




;-------------------------------------
;----- h2-h2 cia ---------------------
;-------------------------------------
if (do_h2h2 eq 1 ) then begin
print, "====== calculating H2-H2 CIA ======"

;;;  --- read h2-h2 cia cia file  ---- ;;;;
inputfile = "/gpfsm/dnb53/etwolf/models/ExoRT/data/cia/orig/H2-H2_2011.cia"
OPENR,lun,inputfile, /Get_Lun

chemsym = ' '
; read block

data1 = dblarr(113,2,9981)  & data1t = dblarr(2,9981)
temp1 = dblarr(113)

for t=0, 112 do begin
  readf,lun,format = '(A20, 2F10, I7, F7, E10, F10)', $
                     chemsym, wmax, wmin, npts, temp, maxcia, res
  readf,lun, data1t
  temp1(t) = temp
  data1(t,*,*) = data1t(*,*)
endfor
FREE_LUN,lun

;print, temp1
nt = n_elements(temp1)
;print, nt
;inputfile = "/Users/wolfe/Models/RT/ExoRT/data/cia/orig/H2-H2_norm_2011.cia"
;OPENR,lun,inputfile, /Get_Lun
;chemsym = ' '
;; read block
;data1n = dblarr(10,2,2428)  & data1tn = dblarr(2,2428)
;temp1n = dblarr(10)
;for t=0, 9 do begin
;  readf,lun,format = '(A20, 2F10, I7, F7, E10, F10)', $
;                     chemsym, wmax, wmin, npts, temp, maxcia, res
;  readf,lun, data1tn
;  temp1n(t) = temp
;  data1n(t,*,*) = data1tn(*,*)
;endfor
;FREE_LUN,lun



;plot block 1
;spectral range    :: 20.20   10000.000  cm-1                                                                              
; temperature range :: 200 - 3000, 25 
for i=0,nt-1 do begin
  wave = data1(i,0,*)
  abs = data1(i,1,*)
  abs = loschmidt*abs*loschmidt
  if i eq 0 then plot, wave, abs, /ylog, xrange=[0, 10000]
  oplot, wave, abs
endfor

;plot block 2
;for i=5,9 do begin
;  wave = data1n(i,0,*)
;  abs = data1n(i,1,*)
;  abs = loschmidt*abs*loschmidt
;;  if i eq 0 then plot, wave, abs, /ylog, xrange=[0, 10000]
;  oplot, wave, abs, color=200
;endfor


;; interp to grid ;; 
nwl = n_elements(NU_LOW)
kabs_out = fltarr(nwl, nt)
for t=0, nt-1 do begin
  wavenumber_vec = data1(t,0,*)
  abs_vec = loschmidt*data1(t,1,*)*loschmidt
 ; abs_vec = data1(t,1,*)
;print, abs_vec
  for d=0,nwl-1 do begin
     x=where(wavenumber_vec le NU_HIGH(d) and wavenumber_vec ge NU_LOW(d),num)
     if (num gt 0) then begin
;       kabs_out(d,t) = double(total(abs_vec(x))/num)
       kabs_out(d,t) = MEDIAN(abs_vec(x))
     endif else begin
       kabs_out(d,t) = 0.0
     endelse
  endfor
endfor
;print, temp1
;print, kabs_out


nwvl = nwl
nt = n_elements(temp1)
TEMPERATURES = temp1
k_out = kabs_out


for h=0, nt-1 do begin
;   oplot, (NU_LOW(*) + NU_HIGH(*))/2., kmaster(*,h), color=200, psym=4, thick=4
   oplot, (NU_LOW(*) + nu_HIGH(*))/2., k_out(*,h), color=200, thick=4, psym=4
endfor

;--- write to netcdf ---
if (do_write eq 1) then begin
outfile = outfile_h2h2
id = NCDF_CREATE(outfile, /CLOBBER)
dim1 = NCDF_DIMDEF(id, 'wbins',nwvl)
dim3 = NCDF_DIMDEF(id, 'tbins',nt)

varid1 = NCDF_VARDEF(id, 'nu_low', dim1, /float)
varid2 = NCDF_VARDEF(id, 'nu_high', dim1, /float)
varid4 = NCDF_VARDEF(id, 'temp', dim3, /float)
varid6 = NCDF_VARDEF(id, 'sigma',[dim1,dim3], /float)


NCDF_ATTPUT, id, varid1, "title", "Wavenumber grid low"
NCDF_ATTPUT, id, varid1, "units", "cm-1"
NCDF_ATTPUT, id, varid2, "title", "Wavenumber grid high"
NCDF_ATTPUT, id, varid2, "units", "cm-1"
NCDF_ATTPUT, id, varid4, "title", "temperature grid"
NCDF_ATTPUT, id, varid4, "units", "K"
NCDF_ATTPUT, id, varid6, "title", "band averaged CIA absorption"
NCDF_ATTPUT, id, varid6, "units", "cm^-1 amagat^-2"

NCDF_CONTROL,id, /ENDEF

NCDF_VARPUT, id, varid1, NU_LOW
NCDF_VARPUT, id, varid2, NU_HIGH
NCDF_VARPUT, id, varid4, TEMPERATURES
NCDF_VARPUT, id, varid6, k_out
NCDF_CLOSE, id

endif

endif





;-------------------------------------
;----- n2-h2 cia ---------------------
;-------------------------------------
if (do_n2h2 eq 1 ) then begin
print, "====== calculating N2-H2 CIA ======"

;;;  --- read n2-h2 cia cia file  ---- ;;;;
inputfile = "/gpfsm/dnb53/etwolf/models/ExoRT/data/cia/orig/N2-H2_2011.cia"
OPENR,lun,inputfile, /Get_Lun

chemsym = ' '
; read block
data1 = dblarr(10,2,1914)  & data1t = dblarr(2,1914)
temp1 = dblarr(10)
for t=0, 9 do begin
  readf,lun,format = '(A20, 2F10, I7, F7, E10, F10)', $
                     chemsym, wmax, wmin, npts, temp, maxcia, res
  readf,lun, data1t
  temp1(t) = temp
  data1(t,*,*) = data1t(*,*)
endfor
FREE_LUN,lun

;plot block 1
for i=0,9 do begin
  wave = data1(i,0,*)
  abs = data1(i,1,*)
  abs = loschmidt*abs*loschmidt
  if i eq 0 then plot, wave, abs, /ylog, xrange=[0, 3200]
  oplot, wave, abs
endfor


;; interp to grid ;; 
nwl = n_elements(nu_low)
kabs_out = fltarr(nwl, 10)
for t=0, 9 do begin
  wavenumber_vec = data1(t,0,*)
  abs_vec = loschmidt*data1(t,1,*)*loschmidt
  for d=0,nwl-1 do begin
     x=where(wavenumber_vec le NU_HIGH(d) and wavenumber_vec ge NU_LOW(d),num)
     ;print, d, nu_high(d), nu_low(d), num
     if (num gt 0) then begin
;       kabs_out(d,t) = total(abs_vec(x))/num
       kabs_out(d,t) = MEDIAN(abs_vec(x))
     endif else begin
       kabs_out(d,t) = 0.0
     endelse
  endfor
endfor
;print, temp1
;print, kabs_out


nwvl = nwl
nt = n_elements(temp1)
TEMPERATURES = temp1
k_out = kabs_out


for h=0, nt-1 do begin
;   oplot, (NU_LOW(*) + NU_HIGH(*))/2., kmaster(*,h), color=200, psym=4,
;   thick=4                                                                     
   oplot, (NU_LOW(*) + nu_HIGH(*))/2., k_out(*,h), color=200, thick=4, psym=4
endfor


;--- write to netcdf ---
if (do_write eq 1) then begin
outfile = outfile_n2h2
id = NCDF_CREATE(outfile, /CLOBBER)
dim1 = NCDF_DIMDEF(id, 'wbins',nwvl)
dim3 = NCDF_DIMDEF(id, 'tbins',nt)

varid1 = NCDF_VARDEF(id, 'nu_low', dim1, /float)
varid2 = NCDF_VARDEF(id, 'nu_high', dim1, /float)
varid4 = NCDF_VARDEF(id, 'temp', dim3, /float)
varid6 = NCDF_VARDEF(id, 'sigma',[dim1,dim3], /float)


NCDF_ATTPUT, id, varid1, "title", "Wavenumber grid low"
NCDF_ATTPUT, id, varid1, "units", "cm-1"
NCDF_ATTPUT, id, varid2, "title", "Wavenumber grid high"
NCDF_ATTPUT, id, varid2, "units", "cm-1"
NCDF_ATTPUT, id, varid4, "title", "temperature grid"
NCDF_ATTPUT, id, varid4, "units", "K"
NCDF_ATTPUT, id, varid6, "title", "band averaged CIA absorption"
NCDF_ATTPUT, id, varid6, "units", "cm-1 amagat-2"

NCDF_CONTROL,id, /ENDEF

NCDF_VARPUT, id, varid1, NU_LOW
NCDF_VARPUT, id, varid2, NU_HIGH
NCDF_VARPUT, id, varid4, TEMPERATURES
NCDF_VARPUT, id, varid6, k_out
NCDF_CLOSE, id

endif

endif

;-------------------------------------
;----- co2-co2 cia ---------------------
;-------------------------------------
if (do_co2co2 eq 1 ) then begin
print, "====== calculating CO2-CO2 CIA ======"

;;;  --- read co2-co2 cia file  ---- ;;;;
inputfile = "/gpfsm/dnb53/etwolf/models/ExoRT/data/cia/orig/CO2-CO2_2018.cia"
OPENR,lun,inputfile, /Get_Lun

chemsym = ' '
; read block 1
data1 = dblarr(10,2,750)  & data1t = dblarr(2,750)
temp1 = dblarr(10)
for t=0, 9 do begin
  readf,lun,format = '(A20, 2F10, I7, F7, E10, F10)', $
                     chemsym, wmax, wmin, npts, temp, maxcia, res
  readf,lun, data1t
  temp1(t) = temp
  data1(t,*,*) = data1t(*,*)
endfor

; read block 2
data2 = dblarr(6,2,801)  & data2t = dblarr(2,801)
temp2 = dblarr(6)
for t=0, 5 do begin
  readf,lun,format = '(A20, 2F10, I7, F7, E10, F10)', $
                     chemsym, wmax, wmin, npts, temp, maxcia, res
  readf,lun, data2t
  temp2(t) = temp
  data2(t,*,*) = data2t(*,*)
endfor

; read block 3
data3 = dblarr(3,2,1361)  & data3t = dblarr(2,1361)
temp3 = dblarr(3)
for t=0, 1 do begin
  readf,lun,format = '(A20, 2F10, I7, F7, E10, F10)', $
                     chemsym, wmax, wmin, npts, temp, maxcia, res
  readf,lun, data3t
  temp3(t) = temp
  data3(t,*,*) = data3t(*,*)
endfor

; read block 4
data4 = dblarr(1,2,1411)  & data4t = dblarr(2,1411)
temp4 = dblarr(1)
for t=0, 0 do begin
  readf,lun,format = '(A20, 2F10, I7, F7, E10, F10)', $
                     chemsym, wmax, wmin, npts, temp;, maxcia, res
  readf,lun, data4t
  temp4(t) = temp
  data4(t,*,*) = data4t(*,*)
endfor

; read block 5
data5 = dblarr(1,2,3530)  & data5t = dblarr(2,3530)
temp5 = dblarr(1)
for t=0, 0 do begin
  readf,lun,format = '(A20, 2F10, I7, F7, E10, F10)', $
                     chemsym, wmax, wmin, npts, temp;, maxcia, res
  readf,lun, data5t
  temp5(t) = temp
  data5(t,*,*) = data5t(*,*)
endfor
FREE_LUN,lun


; interpolat block 3 and 4 together, because they cover the same spectral region
MG = interpol(data4(0,1,*),data4(0,0,*),data3(0,0,*))
temp3(2) = temp4(0) 
print, temp3
data3(2,1,*) = MG(*)
data3(2,0,*) = data3(0,0,*)
;print, data3(2,1,*)
;stop

;interpolate blocks 3 and 5 to the block 1 temperature grid


;plot block 1
for i=0,9 do begin
  wave = data1(i,0,*)
  abs = data1(i,1,*)
  abs = loschmidt*abs*loschmidt
  if i eq 0 then plot, wave, abs, /ylog, xrange=[0, 2000], yrange=[1.e-12, 1.e-3], xstyle=1, ystyle=1
  oplot, wave, abs 
  print, temp1(i)

;help, abs
;  abs2 = data1(i,1,*)/1.385277e-39
;  oplot, wave, abs2
;  for l=0, 100 do begin
;  print, abs(l), abs2(l)
;  endfor

endfor

;plot block 2
for i=0,5 do begin
  wave = data2(i,0,*)
  abs = data2(i,1,*)
  abs = loschmidt*abs*loschmidt
;  if i eq 0 then plot, wave, abs, /ylog, xrange=[0, 3000]
  oplot, wave, abs 
  print, temp2(i)

endfor

;plot block 3
for i=0,2 do begin
  wave = data3(i,0,*)
  abs = data3(i,1,*)
  abs = loschmidt*abs*loschmidt
;  if i eq 0 then plot, wave, abs, /ylog, xrange=[0, 3000]
  oplot, wave, abs 
endfor


;plot block 4
for i=0,0 do begin
  wave = data4(i,0,*)
  abs = data4(i,1,*)
  abs = loschmidt*abs*loschmidt
;  if i eq 0 then p
  oplot, wave, abs ;* 100.
endfor


;plot block 5
for i=0,0 do begin
  wave = data5(i,0,*)
  abs = data5(i,1,*)
  abs = loschmidt*abs*loschmidt
;  if i eq 0 then plot, wave, abs, /ylog, xrange=[0, 3000]
  oplot, wave, abs 
endfor


;; interp to grid ;; 
; interp
nwl = n_elements(nu_low)
k1 = fltarr(nwl, 10)
k2 = fltarr(nwl, 6)
k3 = fltarr(nwl, 3)
k4 = fltarr(nwl, 1)
k5 = fltarr(nwl, 1)

for d=0, nwl-1 do begin

   for t=0, n_elements(temp1)-1 do begin
     wavenumber_vec = data1(t,0,*)  
     abs_vec = loschmidt*data1(t,1,*)*loschmidt
     x1=where(wavenumber_vec le NU_HIGH(d) and wavenumber_vec ge NU_LOW(d),num)
     abs_vec_tot = 0.0
     if (num gt 0) then begin

;; use planck weighting
;      abs_vec_temp = abs_vec(x1)
;      wavenumber_vec_temp = wavenumber_vec(x1)
;      planck1 = fltarr(num)
;      for l=0,num-1 do begin 
;         planck1(l) = planck_nu(wavenumber_vec_temp(l),temp1(t))
;       abs_vec_tot  = abs_vec_tot + abs_vec_temp(l)*planck1(l)
;       endfor
;       k1(d,t) = abs_vec_tot/total(planck1)  ; planck weighting in interval
;\ planck weighting

; use interval average
       k1(d,t) = total(abs_vec(x1))/num 

; use interval median
;       k1(d,t) = MEDIAN(abs_vec(x1))    

; use logarithmic average
;     tempk1=alog10(abs_vec(x1))
;     k1(d,t) =  10.^(mean(tempk1))


;print, d,t, total(abs_vec(x1))/num, MEDIAN(abs_vec(x1)), 10.^(mean(tempk1))


       if (k1(d,t) lt 0.0) then k1(d,t) = 0.0
     endif else begin
       k1(d,t) = 0.0
     endelse
   endfor

   for t=0, n_elements(temp2)-1 do begin
     wavenumber_vec = data2(t,0,*)  
     abs_vec = loschmidt*data2(t,1,*)*loschmidt
     x2=where(wavenumber_vec le NU_HIGH(d) and wavenumber_vec ge NU_LOW(d),num)
     abs_vec_tot = 0.0
     if (num gt 0) then begin

; use planck weighting
;      abs_vec_temp = abs_vec(x2)
;       wavenumber_vec_temp = wavenumber_vec(x2)
;       planck2 = fltarr(num)
;       for l=0,num-1 do begin 
;         planck2(l) = planck_nu(wavenumber_vec_temp(l),temp1(t))
;         abs_vec_tot  = abs_vec_tot + abs_vec_temp(l)*planck2(l)
;       endfor
;       k2(d,t) = abs_vec_tot/total(planck2)
;/ planck weighting

; use interval average
;       k2(d,t) = total(abs_vec(x2))/num       ;interval avergae

; use interval median
;       k2(d,t) = MEDIAN(abs_vec(x2))          ;median value in interval

; use logarithmic average
;     tempk2=alog10(abs_vec(x2))
;     k2(d,t) =  10.^(mean(tempk2))

;print, abs_vec_tot/total(planck2), total(abs_vec(x2))/num,MEDIAN(abs_vec(x2)), 10.^(mean(tempk2))
;print, total(abs_vec(x2))/num,MEDIAN(abs_vec(x2)), 10.^(mean(tempk2))

       if (k2(d,t) lt 0.0) then k2(d,t) = 0.0
     endif else begin
       k2(d,t) = 0.0
     endelse
   endfor

   for t=0, n_elements(temp3)-1 do begin
     wavenumber_vec = data3(t,0,*)  
     abs_vec = loschmidt*data3(t,1,*)*loschmidt
     x3=where(wavenumber_vec le NU_HIGH(d) and wavenumber_vec ge NU_LOW(d),num)
     if (num gt 0) then begin
        print, num
       k3(d,t) = total(abs_vec(x3))/num
;       k3(d,t) = MEDIAN(abs_vec(x3))
;       k3(d,t) = max(abs_vec(x3))
       if (k3(d,t) lt 0.0) then k3(d,t) = 0.0
     endif else begin
       k3(d,t) = 0.0
     endelse
   endfor

;   for t=0, n_elements(temp4)-1 do begin
;     wavenumber_vec = data4(t,0,*)  
;     abs_vec = loschmidt*data4(t,1,*)*loschmidt
;     x4=where(wavenumber_vec le NU_HIGH(d) and wavenumber_vec ge NU_LOW(d),num)
;     if (num gt 0) then begin
;;       k4(d,t) = total(abs_vec(x4))/num
;       k4(d,t) = MEDIAN(abs_vec(x4))
;       if (k4(d,t) lt 0.0) then k4(d,t) = 0.0
;     endif else begin
;       k4(d,t) = 0.0
;     endelse
;   endfor

   for t=0, n_elements(temp5)-1 do begin
     wavenumber_vec = data5(t,0,*)  
     abs_vec = loschmidt*data5(t,1,*)*loschmidt
     x5=where(wavenumber_vec le NU_HIGH(d) and wavenumber_vec ge NU_LOW(d),num)
     if (num gt 0) then begin
        k5(d,t) = total(abs_vec(x5))/num
;       k5(d,t) = MEDIAN(abs_vec(x5))
;       k5(d,t) = max(abs_vec(x5))
       if (k5(d,t) lt 0.0) then k5(d,t) = 0.0
     endif else begin
       k5(d,t) = 0.0
     endelse
   endfor

endfor

;for h=0, n_elements(temp1)-1 do begin
;   oplot, (NU_LOW(*) + nu_HIGH(*))/2., k1(*,h), color=200, thick=4, psym=4, symsize=0.9
;endfor
;for h=0, n_elements(temp2)-1 do begin
;   oplot, (NU_LOW(*) + nu_HIGH(*))/2., k2(*,h), color=200, thick=4, psym=4, symsize=0.9
;endfor
;for h=0, n_elements(temp3+1)-1 do begin
;   oplot, (NU_LOW(*) + nu_HIGH(*))/2., k3(*,h), color=200, thick=4, psym=4, symsize=0.9
;endfor
;;for h=0, n_elements(temp4)-1 do begin
;;   oplot, (NU_LOW(*) + nu_HIGH(*))/2., k4(*,h), color=200, thick=4, psym=4, symsize=0.9
;;endfor
;for h=0, n_elements(temp5)-1 do begin
;   oplot, (NU_LOW(*) + nu_HIGH(*))/2., k5(*,h), color=200, thick=4, psym=4, symsize=0.9
;endfor





;kmaster_lw = fltarr(nwl, 10)
kmaster_lw = fltarr(nwl, 12)
kmaster_sw = fltarr(nwl, 3)

;temp_master_lw = temp1
; add/interpolate to lower temperatures in longwave grid
temp_master_lw=fltarr(12)
temp_master_lw(0) = 100.  & temp_master_lw(1) = 150. 
temp_master_lw(2:11) = temp1(*) 

temp_master_sw = temp3
print, "temp1", temp1
print, "temp2", temp2
print, "temp3", temp3
print, "temp5", temp5
print, "temp_master_lw", temp_master_lw
print, "temp_master_sw", temp_master_sw





;grid up longwave CO2-CO2 CIA
for t=0, 9 do begin
; if temperature grid point is less than block 2 and 3 combined array 
; if (temp_master(t) lt tarr(0)) then kmaster(*,t) = k1(*,t) + k3(*,0)
 kmaster_lw(*,t+2) = k1(*,t)  ;k1 maps to all temperatures
 ;k2 maps to temperatures 1-6.  set the higher temperatures to the max value
 if (temp_master_lw(t+2) le temp_master_lw(5)) then kmaster_lw(*,t+2) = kmaster_lw(*,t+2) + k2(*,t+2)
 if (temp_master_lw(t+2) gt temp_master_lw(5)) then kmaster_lw(*,t+2) = kmaster_lw(*,t+2) + k2(*,5)
endfor


for jj=0, nwl-1 do begin
  kmaster_lw(jj,0) = interpol(kmaster_lw(jj,2:5),temp_master_lw(2:5),100.)
  kmaster_lw(jj,1) = interpol(kmaster_lw(jj,2:5),temp_master_lw(2:5),150.)
  if (kmaster_lw(jj,0) lt 0.0) then kmaster_lw(jj,0) = 0.0
  if (kmaster_lw(jj,1) lt 0.0) then kmaster_lw(jj,1) = 0.0
endfor

;grid up shortwave CO2-CO2 CIA
for t=0, 2 do begin
; if temperature grid point is less than block 2 and 3 combined array 
; if (temp_master(t) lt tarr(0)) then kmaster(*,t) = k1(*,t) + k3(*,0)
 kmaster_sw(*,t) = k3(*,t)+ k5(*,0)  ;k1 maps to all temperatures
 ;k2 maps to temperatures 1-6.  set the higher temperatures to the max value
; if (temp_master_sw(t) le temp_master_sw(5)) then kmaster_sw(*,t) = kmaster_sw(*,t) + k5(*,0)
endfor


; plot

;for h=0,4 do begin
;   oplot, (nu_low(*) + nu_high(*))/2., ktemp2(*,h), color=200, thick=4, psym=4
;   oplot, (nu_low(*) + nu_high(*))/2., ktemp3(*,h), color=200, thick=4, psym=4
;endfor
for h=0, n_elements(temp_master_lw)-1 do begin
   oplot, (nu_low(*) + nu_high(*))/2., kmaster_lw(*,h), psym=4, thick=2

print, temp_master_lw(h)
;stop
endfor

for h=0, 2 do begin
   oplot, (nu_low(*) + nu_high(*))/2., kmaster_sw(*,h), psym=4, thick=4
endfor

;print, temp_master
;print, kmaster
;print, "--"
;print, temp1
;print, ktemp1
;print, temp2
;print, ktemp2
;print, temp3
;print, ktemp3
;print, tarr
;stop


nwvl = nwl
nt_lw = n_elements(temp_master_lw)
TEMPERATURES_lw = temp_master_lw
k_out_lw = kmaster_lw

nt_sw = n_elements(temp_master_sw)
TEMPERATURES_sw = temp_master_sw
k_out_sw = kmaster_sw


;--- write to netcdf ---
if (do_write eq 1) then begin
outfile = outfile_co2co2_lw
id = NCDF_CREATE(outfile, /CLOBBER)
dim1 = NCDF_DIMDEF(id, 'wbins',nwvl)
dim3 = NCDF_DIMDEF(id, 'tbins',nt_lw)

varid1 = NCDF_VARDEF(id, 'nu_low', dim1, /float)
varid2 = NCDF_VARDEF(id, 'nu_high', dim1, /float)
varid4 = NCDF_VARDEF(id, 'temp', dim3, /float)
varid6 = NCDF_VARDEF(id, 'sigma',[dim1,dim3], /float)


NCDF_ATTPUT, id, varid1, "title", "Wavenumber grid low"
NCDF_ATTPUT, id, varid1, "units", "cm-1"
NCDF_ATTPUT, id, varid2, "title", "Wavenumber grid high"
NCDF_ATTPUT, id, varid2, "units", "cm-1"
NCDF_ATTPUT, id, varid4, "title", "temperature grid"
NCDF_ATTPUT, id, varid4, "units", "K"
NCDF_ATTPUT, id, varid6, "title", "band averaged CIA absorption"
NCDF_ATTPUT, id, varid6, "units", "cm-1 amagat-2"

NCDF_CONTROL,id, /ENDEF

NCDF_VARPUT, id, varid1, NU_LOW
NCDF_VARPUT, id, varid2, NU_HIGH
NCDF_VARPUT, id, varid4, TEMPERATURES_lw
NCDF_VARPUT, id, varid6, k_out_lw
NCDF_CLOSE, id


outfile = outfile_co2co2_sw
id = NCDF_CREATE(outfile, /CLOBBER)
dim1 = NCDF_DIMDEF(id, 'wbins',nwvl)
dim3 = NCDF_DIMDEF(id, 'tbins',nt_sw)

varid1 = NCDF_VARDEF(id, 'nu_low', dim1, /float)
varid2 = NCDF_VARDEF(id, 'nu_high', dim1, /float)
varid4 = NCDF_VARDEF(id, 'temp', dim3, /float)
varid6 = NCDF_VARDEF(id, 'sigma',[dim1,dim3], /float)


NCDF_ATTPUT, id, varid1, "title", "Wavenumber grid low"
NCDF_ATTPUT, id, varid1, "units", "cm-1"
NCDF_ATTPUT, id, varid2, "title", "Wavenumber grid high"
NCDF_ATTPUT, id, varid2, "units", "cm-1"
NCDF_ATTPUT, id, varid4, "title", "temperature grid"
NCDF_ATTPUT, id, varid4, "units", "K"
NCDF_ATTPUT, id, varid6, "title", "band averaged CIA absorption"
NCDF_ATTPUT, id, varid6, "units", "cm-1 amagat-2"

NCDF_CONTROL,id, /ENDEF

NCDF_VARPUT, id, varid1, NU_LOW
NCDF_VARPUT, id, varid2, NU_HIGH
NCDF_VARPUT, id, varid4, TEMPERATURES_sw
NCDF_VARPUT, id, varid6, k_out_sw
NCDF_CLOSE, id




endif

endif



;-------------------------------------
;----- co2-h2 cia ---------------------
;-------------------------------------
if (do_co2h2 eq 1 ) then begin
print, "====== calculating CO2-H2 CIA ======"

;;;  --- read co2-h2 cia cia file  ---- ;;;;
inputfile = "/gpfsm/dnb53/etwolf/models/ExoRT/data/cia/orig/CO2-H2_TURBET2020_cm5mol-2.cia"
OPENR,lun,inputfile, /Get_Lun

chemsym = ' '
; read block

data1 = dblarr(6,2,300)  & data1t = dblarr(2,300)
temp1 = dblarr(6)

for t=0, 6-1 do begin
  readf,lun,format = '(A20, 2F10, I7, F7)', $
                     chemsym, wmax, wmin, npts, temp 
  readf,lun, data1t
  temp1(t) = temp
  data1(t,*,*) = data1t(*,*)
endfor
FREE_LUN,lun

;print, temp1
nt = n_elements(temp1)

;plot block 1
;spectral range    :: 20.20   10000.000  cm-1                                                                              
; temperature range :: 200 - 3000, 25 
for i=0,nt-1 do begin
  wave = data1(i,0,*)
  abs = data1(i,1,*)
  abs = loschmidt*abs*loschmidt
  if i eq 0 then plot, wave, abs, /ylog, xrange=[0, 2000]
  oplot, wave, abs
endfor

;; interp to grid ;; 
nwl = n_elements(NU_LOW)
kabs_out = fltarr(nwl, nt)
for t=0, nt-1 do begin
  wavenumber_vec = data1(t,0,*)
  abs_vec = loschmidt*data1(t,1,*)*loschmidt
 ; abs_vec = data1(t,1,*)
  for d=0,nwl-1 do begin
     x=where(wavenumber_vec le NU_HIGH(d) and wavenumber_vec ge NU_LOW(d),num)
     if (num gt 0) then begin
;       kabs_out(d,t) = double(total(abs_vec(x))/num)
       kabs_out(d,t) = MEDIAN(abs_vec(x))
       print, num, nu_low(d), temp1(t), kabs_out(d,t)
     endif else begin
       kabs_out(d,t) = 0.0
     endelse
  endfor
endfor
;print, temp1
;print, kabs_out


nwvl = nwl
nt = n_elements(temp1)
TEMPERATURES = temp1
k_out = kabs_out


for h=0, nt-1 do begin
;   oplot, (NU_LOW(*) + NU_HIGH(*))/2., kmaster(*,h), color=200, psym=4, thick=4
   oplot, (NU_LOW(*) + nu_HIGH(*))/2., k_out(*,h), color=200, thick=4, psym=4
endfor

;--- write to netcdf ---
if (do_write eq 1) then begin
outfile = outfile_co2h2
id = NCDF_CREATE(outfile, /CLOBBER)
dim1 = NCDF_DIMDEF(id, 'wbins',nwvl)
dim3 = NCDF_DIMDEF(id, 'tbins',nt)

varid1 = NCDF_VARDEF(id, 'nu_low', dim1, /float)
varid2 = NCDF_VARDEF(id, 'nu_high', dim1, /float)
varid4 = NCDF_VARDEF(id, 'temp', dim3, /float)
varid6 = NCDF_VARDEF(id, 'sigma',[dim1,dim3], /float)


NCDF_ATTPUT, id, varid1, "title", "Wavenumber grid low"
NCDF_ATTPUT, id, varid1, "units", "cm-1"
NCDF_ATTPUT, id, varid2, "title", "Wavenumber grid high"
NCDF_ATTPUT, id, varid2, "units", "cm-1"
NCDF_ATTPUT, id, varid4, "title", "temperature grid"
NCDF_ATTPUT, id, varid4, "units", "K"
NCDF_ATTPUT, id, varid6, "title", "band averaged CIA absorption"
NCDF_ATTPUT, id, varid6, "units", "cm^-1 amagat^-2"

NCDF_CONTROL,id, /ENDEF

NCDF_VARPUT, id, varid1, NU_LOW
NCDF_VARPUT, id, varid2, NU_HIGH
NCDF_VARPUT, id, varid4, TEMPERATURES
NCDF_VARPUT, id, varid6, k_out
NCDF_CLOSE, id

endif

endif


;-------------------------------------
;----- co2-ch4 cia ---------------------
;-------------------------------------
if (do_co2ch4 eq 1 ) then begin
print, "====== calculating CO2-CH4 CIA ======"

;;;  --- read co2-ch4 cia cia file  ---- ;;;;
inputfile = "/gpfsm/dnb53/etwolf/models/ExoRT/data/cia/orig/CO2-CH4_TURBET2020_cm5mol-2.cia"
OPENR,lun,inputfile, /Get_Lun

chemsym = ' '
; read block

data1 = dblarr(6,2,240)  & data1t = dblarr(2,240)
temp1 = dblarr(6)

for t=0, 6-1 do begin
  readf,lun,format = '(A20, 2F10, I7, F7)', $
                     chemsym, wmax, wmin, npts, temp 
  readf,lun, data1t
  temp1(t) = temp
  data1(t,*,*) = data1t(*,*)
endfor
FREE_LUN,lun

;print, temp1
nt = n_elements(temp1)

;plot block 1
;spectral range    :: 20.20   10000.000  cm-1                                                                              
; temperature range :: 200 - 3000, 25 
for i=0,nt-1 do begin
  wave = data1(i,0,*)
  abs = data1(i,1,*)
  abs = loschmidt*abs*loschmidt
  if i eq 0 then plot, wave, abs, /ylog, xrange=[0, 2000]
  oplot, wave, abs
endfor

;; interp to grid ;; 
nwl = n_elements(NU_LOW)
kabs_out = fltarr(nwl, nt)
for t=0, nt-1 do begin
  wavenumber_vec = data1(t,0,*)
  abs_vec = loschmidt*data1(t,1,*)*loschmidt
 ; abs_vec = data1(t,1,*)
  for d=0,nwl-1 do begin
     x=where(wavenumber_vec le NU_HIGH(d) and wavenumber_vec ge NU_LOW(d),num)
     if (num gt 0) then begin
;       kabs_out(d,t) = double(total(abs_vec(x))/num)
       kabs_out(d,t) = MEDIAN(abs_vec(x))
       if (d eq 20) then  kabs_out(d,t) = 0.0
       print, d, num, nu_low(d), temp1(t), kabs_out(d,t)
     endif else begin
       kabs_out(d,t) = 0.0
     endelse
  endfor
endfor
;print, temp1
;print, kabs_out


nwvl = nwl
nt = n_elements(temp1)
TEMPERATURES = temp1
k_out = kabs_out


for h=0, nt-1 do begin
;   oplot, (NU_LOW(*) + NU_HIGH(*))/2., kmaster(*,h), color=200, psym=4, thick=4
   oplot, (NU_LOW(*) + nu_HIGH(*))/2., k_out(*,h), color=200, thick=4, psym=4
endfor

;--- write to netcdf ---
if (do_write eq 1) then begin
outfile = outfile_co2ch4
id = NCDF_CREATE(outfile, /CLOBBER)
dim1 = NCDF_DIMDEF(id, 'wbins',nwvl)
dim3 = NCDF_DIMDEF(id, 'tbins',nt)

varid1 = NCDF_VARDEF(id, 'nu_low', dim1, /float)
varid2 = NCDF_VARDEF(id, 'nu_high', dim1, /float)
varid4 = NCDF_VARDEF(id, 'temp', dim3, /float)
varid6 = NCDF_VARDEF(id, 'sigma',[dim1,dim3], /float)


NCDF_ATTPUT, id, varid1, "title", "Wavenumber grid low"
NCDF_ATTPUT, id, varid1, "units", "cm-1"
NCDF_ATTPUT, id, varid2, "title", "Wavenumber grid high"
NCDF_ATTPUT, id, varid2, "units", "cm-1"
NCDF_ATTPUT, id, varid4, "title", "temperature grid"
NCDF_ATTPUT, id, varid4, "units", "K"
NCDF_ATTPUT, id, varid6, "title", "band averaged CIA absorption"
NCDF_ATTPUT, id, varid6, "units", "cm^-1 amagat^-2"

NCDF_CONTROL,id, /ENDEF

NCDF_VARPUT, id, varid1, NU_LOW
NCDF_VARPUT, id, varid2, NU_HIGH
NCDF_VARPUT, id, varid4, TEMPERATURES
NCDF_VARPUT, id, varid6, k_out
NCDF_CLOSE, id

endif

endif




;-------------------------------------
;----- h2o-h2o cia ---------------------
;-------------------------------------
if (do_h2oh2o eq 1 ) then begin
print, "====== calculating H2O-H2O CIA ======"

;experimental H2O-H2O CIA derived by Villanueva and Kaufmann to replace MT_CKD model

;;;  --- read h2o-h2o cia cia file  ---- ;;;;
inputfile = "/gpfsm/dnb53/etwolf/models/ExoRT/data/cia/orig/H2O-H2O_2020.cia"
OPENR,lun,inputfile, /Get_Lun

chemsym = ' '
; read block

data1 = dblarr(21,2,2001)  & data1t = dblarr(2,2001)
temp1 = dblarr(21)

for t=0, 21-1 do begin
  readf,lun,format = '(A20, 2F10, I7, F7)', $
                     chemsym, wmax, wmin, npts, temp 
  readf,lun, data1t
  temp1(t) = temp
  data1(t,*,*) = data1t(*,*)
endfor
FREE_LUN,lun

;print, temp1
nt = n_elements(temp1)

;plot block 1
;spectral range    :: 0.0  20000.000  cm-1                                                                              
; temperature range :: 250 - 450, 21 
for i=0,nt-1 do begin
  wave = data1(i,0,*)
  abs = data1(i,1,*)
  if i eq 0 then plot, wave, abs, /ylog, xrange=[0, 20000]
  oplot, wave, abs
  abs = loschmidt*abs*loschmidt
endfor

;; interp to grid ;; 
nwl = n_elements(NU_LOW)
kabs_out = fltarr(nwl, nt)
for t=0, nt-1 do begin
  wavenumber_vec = data1(t,0,*)
  abs_vec = loschmidt*data1(t,1,*)*loschmidt
 ; abs_vec = data1(t,1,*)
  for d=0,nwl-1 do begin
     x=where(wavenumber_vec le NU_HIGH(d) and wavenumber_vec ge NU_LOW(d),num)
     if (num gt 0) then begin
       kabs_out(d,t) = double(total(abs_vec(x))/num)
;       kabs_out(d,t) = MEDIAN(abs_vec(x))
       if (d eq 20) then  kabs_out(d,t) = 0.0
       print, d, num, nu_low(d), temp1(t), kabs_out(d,t)
     endif else begin
       kabs_out(d,t) = 0.0
     endelse
  endfor
endfor
;print, temp1
;print, kabs_out


nwvl = nwl
nt = n_elements(temp1)
TEMPERATURES = temp1
k_out = kabs_out


for h=0, nt-1 do begin
;   oplot, (NU_LOW(*) + NU_HIGH(*))/2., kmaster(*,h), color=200, psym=4, thick=4
   oplot, (NU_LOW(*) + nu_HIGH(*))/2., k_out(*,h), color=200, thick=4, psym=4
endfor

;--- write to netcdf ---
if (do_write eq 1) then begin
outfile = outfile_h2oh2o
id = NCDF_CREATE(outfile, /CLOBBER)
dim1 = NCDF_DIMDEF(id, 'wbins',nwvl)
dim3 = NCDF_DIMDEF(id, 'tbins',nt)

varid1 = NCDF_VARDEF(id, 'nu_low', dim1, /float)
varid2 = NCDF_VARDEF(id, 'nu_high', dim1, /float)
varid4 = NCDF_VARDEF(id, 'temp', dim3, /float)
varid6 = NCDF_VARDEF(id, 'sigma',[dim1,dim3], /float)


NCDF_ATTPUT, id, varid1, "title", "Wavenumber grid low"
NCDF_ATTPUT, id, varid1, "units", "cm-1"
NCDF_ATTPUT, id, varid2, "title", "Wavenumber grid high"
NCDF_ATTPUT, id, varid2, "units", "cm-1"
NCDF_ATTPUT, id, varid4, "title", "temperature grid"
NCDF_ATTPUT, id, varid4, "units", "K"
NCDF_ATTPUT, id, varid6, "title", "band averaged CIA absorption"
NCDF_ATTPUT, id, varid6, "units", "cm^-1 amagat^-2"

NCDF_CONTROL,id, /ENDEF

NCDF_VARPUT, id, varid1, NU_LOW
NCDF_VARPUT, id, varid2, NU_HIGH
NCDF_VARPUT, id, varid4, TEMPERATURES
NCDF_VARPUT, id, varid6, k_out
NCDF_CLOSE, id

endif

endif


;-------------------------------------
;----- h2o-n2 cia ---------------------
;-------------------------------------
if (do_h2on2 eq 1 ) then begin
print, "====== calculating H2O-N2 CIA ======"

;experimental H2O-N2 CIA derived by Villanueva and Kaufmann to replace MT_CKD model

;;;  --- read h2o-n2 cia cia file  ---- ;;;;
inputfile = "/gpfsm/dnb53/etwolf/models/ExoRT/data/cia/orig/H2O-N2_2020.cia"
OPENR,lun,inputfile, /Get_Lun

chemsym = ' '
; read block

data1 = dblarr(21,2,2001)  & data1t = dblarr(2,2001)
temp1 = dblarr(21)

for t=0, 21-1 do begin
  readf,lun,format = '(A20, 2F10, I7, F7)', $
                     chemsym, wmax, wmin, npts, temp 
  readf,lun, data1t
  temp1(t) = temp
  data1(t,*,*) = data1t(*,*)
endfor
FREE_LUN,lun

;print, temp1
nt = n_elements(temp1)

;plot block 1
;spectral range    :: 0.0  20000.000  cm-1                                                                              
; temperature range :: 250 - 450, 21 
for i=0,nt-1 do begin
  wave = data1(i,0,*)
  abs = data1(i,1,*)
  if i eq 0 then plot, wave, abs, /ylog, xrange=[0, 20000]
  oplot, wave, abs
  abs = loschmidt*abs*loschmidt
endfor

;; interp to grid ;; 
nwl = n_elements(NU_LOW)
kabs_out = fltarr(nwl, nt)
for t=0, nt-1 do begin
  wavenumber_vec = data1(t,0,*)
  abs_vec = loschmidt*data1(t,1,*)*loschmidt
 ; abs_vec = data1(t,1,*)
  for d=0,nwl-1 do begin
     x=where(wavenumber_vec le NU_HIGH(d) and wavenumber_vec ge NU_LOW(d),num)
     if (num gt 0) then begin
;       kabs_out(d,t) = double(total(abs_vec(x))/num)
       kabs_out(d,t) = MEDIAN(abs_vec(x))
       if (d eq 20) then  kabs_out(d,t) = 0.0
       print, d, num, nu_low(d), temp1(t), kabs_out(d,t)
     endif else begin
       kabs_out(d,t) = 0.0
     endelse
  endfor
endfor
;print, temp1
;print, kabs_out


nwvl = nwl
nt = n_elements(temp1)
TEMPERATURES = temp1
k_out = kabs_out


for h=0, nt-1 do begin
;   oplot, (NU_LOW(*) + NU_HIGH(*))/2., kmaster(*,h), color=200, psym=4, thick=4
   oplot, (NU_LOW(*) + nu_HIGH(*))/2., k_out(*,h), color=200, thick=4, psym=4
endfor

;--- write to netcdf ---
if (do_write eq 1) then begin
outfile = outfile_h2oh2o
id = NCDF_CREATE(outfile, /CLOBBER)
dim1 = NCDF_DIMDEF(id, 'wbins',nwvl)
dim3 = NCDF_DIMDEF(id, 'tbins',nt)

varid1 = NCDF_VARDEF(id, 'nu_low', dim1, /float)
varid2 = NCDF_VARDEF(id, 'nu_high', dim1, /float)
varid4 = NCDF_VARDEF(id, 'temp', dim3, /float)
varid6 = NCDF_VARDEF(id, 'sigma',[dim1,dim3], /float)


NCDF_ATTPUT, id, varid1, "title", "Wavenumber grid low"
NCDF_ATTPUT, id, varid1, "units", "cm-1"
NCDF_ATTPUT, id, varid2, "title", "Wavenumber grid high"
NCDF_ATTPUT, id, varid2, "units", "cm-1"
NCDF_ATTPUT, id, varid4, "title", "temperature grid"
NCDF_ATTPUT, id, varid4, "units", "K"
NCDF_ATTPUT, id, varid6, "title", "band averaged CIA absorption"
NCDF_ATTPUT, id, varid6, "units", "cm^-1 amagat^-2"

NCDF_CONTROL,id, /ENDEF

NCDF_VARPUT, id, varid1, NU_LOW
NCDF_VARPUT, id, varid2, NU_HIGH
NCDF_VARPUT, id, varid4, TEMPERATURES
NCDF_VARPUT, id, varid6, k_out
NCDF_CLOSE, id

endif

endif


end
