pro makeCIA

do_write = 1
do_n2n2 = 1
do_h2h2 = 1
do_n2h2 = 1


outfile_n2n2 = "N2-N2_cia_68bin.nc"
outfile_h2h2 = "H2-H2_cia_68bin.nc"
outfile_n2h2 = "N2-H2_cia_68bin.nc"

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

NU_HIGH =   [ 40.00000,   100.0000,   160.0000, $
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
inputfile = "/projects/wolfet/models/ExoRT/data/cia/orig/N2-N2_2011.cia"
OPENR,lun,inputfile, /Get_Lun

chemsym = ' '
; read first block
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

;plot block 1
for i=0,9 do begin
  wave = data1(i,0,*)
  abs = data1(i,1,*)
  abs = loschmidt*abs*loschmidt
  if i eq 0 then plot, wave, abs, /ylog, xrange=[0, 3000], /nodata
  oplot, wave, abs
endfor

;plot block 2
for i=0,4 do begin
  wave = data2(i,0,*)
  abs = data2(i,1,*)
;  if i eq 0 then plot, wave, abs
  abs = loschmidt*abs*loschmidt
  oplot, wave, abs
endfor

;plot block 3
for i=0,4 do begin
  wave = data3(i,0,*)
  abs = data3(i,1,*)
;  if i eq 0 then plot, wave, abs
  abs = loschmidt*abs*loschmidt
  oplot, wave, abs
endfor



; interp
nwl = n_elements(nu_low)
ktemp1 = fltarr(nwl, 10)
ktemp2 = fltarr(nwl, 5)
ktemp3 = fltarr(nwl, 5)
for d=0, nwl-1 do begin

   for t=0, n_elements(temp1)-1 do begin
     wavenumber_vec = data1(t,0,*)  
     abs_vec = loschmidt*data1(t,1,*)*loschmidt
     x1=where(wavenumber_vec le NU_HIGH(d) and wavenumber_vec ge NU_LOW(d),num)
     if (num gt 0) then begin
       ktemp1(d,t) = total(abs_vec(x1))/num
       if (ktemp1(d,t) lt 0.0) then ktemp3(d,t) = 0.0
     endif else begin
       ktemp1(d,t) = 0.0
     endelse
   endfor

   for t=0, n_elements(temp2)-1 do begin
     wavenumber_vec = data2(t,0,*)  
     abs_vec = loschmidt*data2(t,1,*)*loschmidt
     x2=where(wavenumber_vec le NU_HIGH(d) and wavenumber_vec ge NU_LOW(d),num)
     if (num gt 0) then begin
       ktemp2(d,t) = total(abs_vec(x2))/num
       if (ktemp2(d,t) lt 0.0) then ktemp2(d,t) = 0.0
     endif else begin
       ktemp2(d,t) = 0.0
     endelse
   endfor


   for t=0, n_elements(temp3)-1 do begin
     wavenumber_vec = data3(t,0,*)  
     abs_vec = loschmidt*data3(t,1,*)*loschmidt
     x3=where(wavenumber_vec le NU_HIGH(d) and wavenumber_vec ge NU_LOW(d),num)
     if (num gt 0) then begin
       ktemp3(d,t) = total(abs_vec(x3))/num
       if (ktemp3(d,t) lt 0.0) then ktemp3(d,t) = 0.0
     endif else begin
       ktemp3(d,t) = 0.0
     endelse
   endfor

endfor

;create interp array
karr = fltarr(nwl, 10)
kmaster_temp = fltarr(nwl, 10)
kmaster = fltarr(nwl, 10)
tarr = fltarr(10)
for k=0, 4 do begin
  karr(*,k) = ktemp3(*,k)
  tarr(k) = temp3(k)
endfor
for k=5, 9 do begin
  karr(*,k) = ktemp2(*,k-5)
  tarr(k) = temp2(k-5)
endfor


;print, tarr, karr

temp_master = temp1
for d=0,nwl-1 do begin
    kmaster_temp(d,*) = interpol(karr(d,*), tarr ,temp_master) ; + ktemp1(d,*)
endfor

for t=0, 9 do begin
  if (temp_master(t) lt tarr(0)) then kmaster(*,t) = ktemp1(*,t) + ktemp3(*,0)
  if (temp_master(t) gt tarr(9)) then kmaster(*,t) = ktemp1(*,t) + ktemp2(*,4)
  if (temp_master(t) gt tarr(0) and (temp_master(t) lt tarr(9))) then kmaster(*,t) = ktemp1(*,t) + kmaster_temp(*,t)
endfor



;for h=0,4 do begin
;   oplot, (nu_low(*) + nu_high(*))/2., ktemp2(*,h), color=200, thick=4, psym=4
;   oplot, (nu_low(*) + nu_high(*))/2., ktemp3(*,h), color=200, thick=4, psym=4
;endfor
for h=0, 9 do begin
   oplot, (nu_low(*) + nu_high(*))/2., kmaster(*,h), color=200, psym=4, thick=4
   oplot, (nu_low(*) + nu_high(*))/2., ktemp1(*,h), color=200, thick=4, psym=4
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
;----- h2-h2 cia ---------------------
;-------------------------------------
if (do_h2h2 eq 1 ) then begin
print, "====== calculating H2-H2 CIA ======"

;;;  --- read h2-h2 cia cia file  ---- ;;;;
inputfile = "/projects/wolfet/models/ExoRT/data/cia/orig/H2-H2_2011.cia"
OPENR,lun,inputfile, /Get_Lun

chemsym = ' '
; read block

data1 = dblarr(11,2,9981)  & data1t = dblarr(2,9981)
temp1 = dblarr(11)

for t=0, 10 do begin
  readf,lun,format = '(A20, 2F10, I7, F7, E10, F10)', $
                     chemsym, wmax, wmin, npts, temp, maxcia, res
  readf,lun, data1t
  temp1(t) = temp
  data1(t,*,*) = data1t(*,*)
endfor
FREE_LUN,lun


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
for i=0,10 do begin
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
nwl = n_elements(nu_low)
kabs_out = fltarr(nwl, 11)
for t=0, 10 do begin
  wavenumber_vec = data1(t,0,*)
  abs_vec = loschmidt*data1(t,1,*)*loschmidt
  for d=0,nwl-1 do begin
     x=where(wavenumber_vec le NU_HIGH(d) and wavenumber_vec ge NU_LOW(d),num)
     ;print, d, nu_high(d), nu_low(d), num
     if (num gt 0) then begin
       kabs_out(d,t) = total(abs_vec(x))/num
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
;----- n2-h2 cia ---------------------
;-------------------------------------
if (do_n2h2 eq 1 ) then begin
print, "====== calculating N2-H2 CIA ======"

;;;  --- read h2-h2 cia cia file  ---- ;;;;
inputfile = "/projects/wolfet/models/ExoRT/data/cia/orig/N2-H2_2011.cia"
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
  if i eq 0 then plot, wave, abs, /ylog, xrange=[0, 3000]
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
       kabs_out(d,t) = total(abs_vec(x))/num
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



end
