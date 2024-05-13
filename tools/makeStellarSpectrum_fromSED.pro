pro makeStellarSpectrum_fromSED
;------------------------------------------------
;Author: Wolf, E.T.
;Created: August 2020
;Revision history: adapted from old code 
;                                                         
;DESCRIPTION:
; Creates solar spectral files for ExoRT.  
;  Reads data read in from SED text files and sorts.
;
;------------------------------------------------
;NOTES: 
;
;-- choose one and only one -- 
;-- spectral resolution --
do_n28 = 0
do_n42 = 0
do_n68 = 1
do_n73 = 0
do_n84 = 0


do_write=1   ; do I want to write the file to netcdf?
Snorm=1360.0 ; normalization of SED in W/m2

;-- paths to raw stellar spectra --
spath1 = '/discover/nobackup/etwolf/models/ExoRT/data/solar/raw'
spath2 = '/discover/nobackup/projects/p54/users/etwolf/raw_stellar_spectra/models_bt-settl_1577677702'
spath3 = '/discover/nobackup/projects/p54/users/etwolf/raw_stellar_spectra'
spath4 = '/discover/nobackup/etwolf/models/ExoRT/data/solar/raw/white_dwarfs'

;--
;-- raw stellar spectral files --
;--
;nhead = number of lines of header information in the raw stellar spectral file
;ncol = number of columns in the raw stellar spectral file
;npts = number of points/rows in the raw stellar spectral files

;file = spath1 + '/g2v_sun.txt' & npts = 100000  & outname = "G2V_SUN_n28_test.nc"  & nhead=1 & ncol = 2
;file = spath1 + '/adleo_dat.txt' & npts = 4838 & header=strarr(175) & outname=  "M4.5_adleo_n68.nc" & nhead = 175 & ncol = 5 ;M4.5 3390K
;file = spath1 + '/trappist1_sed.txt' & npts = 83302-18   & outname = 'trappist1_lincowski2018_n68.nc' & nhead=18 & ncol = 2
;file = spath1 + '/hd128167um.txt' & npts = 1737 & outname = "F2V_hd128167.nc" ;F2V, 6595 K
;file = spath1 + '/hd22049um.txt' & npts = 1268 & outname = "K2V_hd22049.nc" ;K2V 5084 K
;file = spath1 + '/gj644_dat.txt' & npts = 5483 & header=strarr(98) & outname = "M3.5v_gj644.nc"
;file = spath1 + '/hd114710um.txt' & npts = 1268 & outname = "G0V_hd114710.nc" ;G0V low activity 5860 K
;file = spath1 + '/hd206860um.txt' & npts = 1268 & outname = "G0V_hd206860.nc" ;G0V high activity 5957 K
;file = spath1 + '/lhs1140_bt-settl-interp_3216K_logg5_meta-0.24.txt' & npts =389369 & outname = "LHS1140"
;file = spath1 + '/wolf1069.txt' & npts = 392387  & outname = "wolf1069_n68.nc"  & nhead=2 & ncol = 2
;file = spath1 + '/kepler62.txt' & npts = 51588  & outname = "kepler62_n84.nc"  & nhead=2 & ncol = 2
file = spath1 + '/TOI700_SED_HST.txt' & npts = 97656  & outname = "toi700_hst.nc"  & nhead=1 & ncol = 2

; used in Kopparapu et al. 2017
;file = spath1 + '/bt-settl_2600_logg4.5_FeH0.txt'  & npts = 390628  & nhead = 9 & ncol = 2
;file = spath1 + '/bt-settl_3000_logg4.5_FeH0.txt'  & npts = 394036  & nhead = 9 & ncol = 2
;file = spath1 + '/bt-settl_3300_logg4.5_FeH0.txt'  & npts = 396483  & nhead = 9 & ncol = 2
;file = spath1 + '/bt-settl_3700_logg4.5_FeH0.txt'  & npts = 398395  & nhead = 9 & ncol = 2
;file = spath1 + '/bt-settl_4000_logg4.5_FeH0.txt'  & npts = 398760  & nhead = 9 & ncol = 2
;file = spath1 + '/bt-settl_4500_logg4.5_FeH0.txt'  & npts = 398558  & nhead = 9 & ncol = 2
;file = spath1 + '/lte032-5.0-0.0a+0.0.BT-NextGen.7.dat.txt' & npts = 390737

;high temperature BT_Settl
;file = spath3 + '/lte055-4.5-0.0a+0.0.BT-NextGen.7.dat.txt'  & npts = 395327  & nhead = 21  & ncol = 2
;file = spath3 + '/lte060-4.5-0.0a+0.0.BT-NextGen.7.dat.txt'  & npts = 395272  & nhead = 21  & ncol = 2
;file = spath3 + '/lte065-4.5-0.0a+0.0.BT-NextGen.7.dat.txt'  & npts = 395343  & nhead = 21  & ncol = 2
;file = spath3 + '/lte070-4.5-0.0a+0.0.BT-NextGen.7.dat.txt'  & npts = 395279  & nhead = 21  & ncol = 2
;file = spath3 + '/lte080-4.5-0.0a+0.0.BT-NextGen.7.dat.txt'  & npts = 395704  & nhead = 8   & ncol = 2
;file = spath3 + '/lte090-4.5-0.0a+0.0.BT-NextGen.7.dat.txt'  & npts = 395645  & nhead = 8   & ncol = 2
;file = spath3 + '/lte100-4.5-0.0a+0.0.BT-NextGen.7.dat.txt'  & npts = 395572  & nhead = 8   & ncol = 2
;file = spath3 + '/lte200-4.5-0.0a+0.0.BT-NextGen.7.dat.txt'  & npts = 394974  & nhead = 8   & ncol = 2
;file = spath3 + '/lte300-4.5-0.0a+0.0.BT-NextGen.7.dat.txt'  & npts = 394921  & nhead = 8   & ncol = 2
;file = spath3 + '/lte400-4.5-0.0a+0.0.BT-NextGen.7.dat.txt'  & npts = 394759  & nhead = 8   & ncol = 2
;file = spath3 + '/lte500-4.5-0.0a+0.0.BT-NextGen.7.dat.txt'  & npts = 394629  & nhead = 8   & ncol = 2
;file = spath3 + '/lte600-4.5-0.0a+0.0.BT-NextGen.7.dat.txt'  & npts = 394698  & nhead = 8   & ncol = 2
;file = spath3 + '/lte700-4.5-0.0a+0.0.BT-NextGen.7.dat.txt'  & npts = 394565  & nhead = 8   & ncol = 2

; white dwarf spectra (from A.Shields)
;file = spath4 + '/WD_3000_spec_new_wider.txt'  & npts = 100000  & nhead = 1  & ncol = 2 
;file = spath4 + '/WD_5000_spec_new_wider.txt'  & npts = 100000  & nhead = 1  & ncol = 2 
;file = spath4 + '/WD_8000_spec_new_wider.txt'  & npts = 100000  & nhead = 1  & ncol = 2 
;file = spath4 + '/WD_10000_spec_new_wider.txt'  & npts = 100000  & nhead = 1  & ncol = 2 
;file = spath4 + '/WD_20000_spec_new_wider.txt'  & npts = 100000  & nhead = 1  & ncol = 2 
;file = spath4 + '/WD_40000_spec_new_wider.txt'  & npts = 100000  & nhead = 1  & ncol = 2 
;file = spath4 + '/WD_80000_spec_new_wider.txt'  & npts = 100000  & nhead = 1  & ncol = 2 
;file = spath1 + '/M_star_spec_new_wider.txt'   & npts =  100000  & nhead = 1  & ncol = 2
;file = spath1 + '/G_star_spec_new_wider.txt'   & npts =  100000  & nhead = 1  & ncol = 2


; output file name
;outname = "WD_5000K_n68.nc"

; set scaler to convert input wavelengths to microns
;mu_scale = 1.0      ; input already in microns
mu_scale = 1.0e-4   ; input in angstroms, convert to microns

filename=STRJOIN(STRSPLIT(file,/EXTRACT,' '))
print, filename


SunSolidAngle = 6.87e-5 ; as the Sun appears from Earth today


; define arrays
; check that nhead, ncol, and npts are set correctly
header=strarr(nhead)
data=dblarr(ncol,npts)
wnm=dblarr(npts)
wlgth=dblarr(npts)
;dwm = 20.0
sunm=dblarr(npts)
dwm = dblarr(npts)




; read text data
OPENR,lun,filename, /GET_LUN
READF,lun,header
READF,lun, data
FREE_LUN,lun


for i=0, npts-1 do begin

  ;check your stellar flux raw input file for its particular 
  ;input units.  Convert appropriately to wavelength in microns, 
  ;wavenumber in cm-1.  The irradiance units can be arbitrary 
  ;because they are renormalized to 1360 Wm-2 anyways.
  
  ; set wavelength to microns
  wlgth(i) = data(0,i) * mu_scale 

  ; wavenumber
  wnm(i) = 1.0e4/wlgth(i)  ; wavenumber
  sunm(i) = data(1,i)

  ;print, wlgth(i), sunm(i)
endfor

;for i=0, npts-2 do begin
;  dwm(i) = wnm(i+1)-wnm(i)
;  print, dwm(i)
;endfor

wm1=wnm(0)
wm2=wnm(npts-1)

set_plot,'PS'
!P.Multi=[0,1,0]
loadct,40
!P.font=0
device,xsize=7.0,ysize=5.0,xoff=1.0,yoff=3.0,/inches
device,file='idl.ps'
device,/color,BITS=8 ;, /ENCAPSULATED, /CMYK                                                  
;device,xsize=8.4,ysize=9,xoff=1.0,yoff=1.0,/CM                                               

device, set_font='Helvetica-Oblique', FONT_INDEX=20
device, set_font='Helvetica-Bold', FONT_INDEX=19
device, set_font='helvetica',FONT_INDEX=18


;sunm(*) = sunm(*) * SunSolidAngle 
plot, wlgth, sunm, xtitle = "!18wavelength (microns)", xrange=[0,100]

device, /close
set_plot,'X'

;; -- supported spectral resolutions ---
;=  rt = 28
if (do_n28 eq 1) then begin
  print, "using 28 spectral intervals"
  nrtwavl = 28
  rtwavlow = dblarr(nrtwavl)
  rtwavhi = dblarr(nrtwavl)
  rtwavlow = [10., 350., 500., 630., 700., 820., 980., 1100., 1180., 1390., $
              1480., 1800., 2080., 2200., 2380., 2600., 3250., 4000., 4650., $
              5150., 6150., 7700., 8050., 12850., 16000., 22650., 29000., $
              38000. ]
  rtwavhi = [350., 500., 630., 700., 820., 980., 1100., 1180., 1390., $
            1480., 1800., 2080., 2200., 2380., 2600., 3250., 4000., 4650., $
            5150., 6150., 7700., 8050., 12850., 16000., 22650., 29000., $
            38000., 50000. ]
endif

; = rt = 42
if (do_n42 eq 1) then begin
  print, "using 43 spectral intervals"
  nrtwavl = 42                                                             
  rtwavlow = dblarr(nrtwavl)                                               
  rtwavhi = dblarr(nrtwavl)                                                
  rtwavlow =  [  10.0,       200.0,         350.0,        425.0, $         
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
  rtwavhi =   [ 200.0,         350.0,        425.0, $                      
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


; = rt = 68
if (do_n68 eq 1) then begin
  print, "using 68 spectral intervals"
  nrtwavl = 68                                                             
  rtwavlow = dblarr(nrtwavl) 
  rtwavhi = dblarr(nrtwavl)  
  rtwavlow =  [1.0,   40.00000,   100.0000,   160.0000, $
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
  rtwavhi =   [40.00000,   100.0000,   160.0000, $
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


; = rt = 73
if (do_n73 eq 1) then begin
  print, "using 73 spectral intervals"
  nrtwavl = 73
  rtwavlow = dblarr(nrtwavl)
  rtwavhi = dblarr(nrtwavl)   
  rtwavlow  = [ 0.00000, 40.0000, 100.000, 160.000, 220.000, 280.000, 330.000, $
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
  rtwavhi =[ 0.00000, 40.0000, 100.000, 160.000, 220.000, 280.000, 330.000, $
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
  print, "using 84 spectral intervals (shortwave extension of n68)"
  nrtwavl = 84
  rtwavlow = dblarr(nrtwavl)
  rtwavhi = dblarr(nrtwavl)
  rtwavlow =  [1.0E+00,   40.00000,   100.0000,   160.0000, $
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
  rtwavhi = [40.00000,   100.0000,   160.0000, $
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



rtwavmid = (rtwavhi(*) + rtwavlow(*) )/2.0
rtwvlmid = (1.0e4/rtwavlow(*) + 1.0e4/rtwavhi(*))/2.0
rtwvldel = (1.0e4/rtwavlow(*) - 1.0e4/rtwavhi(*))
solarflux = dblarr(nrtwavl)
solarflux(*) = 0.0 

for k=0, nrtwavl-1 do begin
  print, rtwavhi(k), rtwavlow(k), rtwavmid(k)
endfor

;-- bin raw stellar into GCM grids ---
nxsum=0.0
for iw=0, nrtwavl-1 do begin
  ftmp=0.
  wn_lw=rtwavlow(iw)
  wn_hi=rtwavhi(iw)
; bin averages
  dwl = 1.e04/wn_lw - 1.0e4/wn_hi
  x=where(wnm ge wn_lw and wnm le wn_hi, nx)
  nxsum=nxsum+nx
  if (nx gt 0) then begin
    ftmp = mean(sunm(x))*dwl
  endif else begin
    ftmp = 0.0
  endelse
  print, iw+1,  wn_lw, wn_hi, nx, ftmp
  if(ftmp ge 0.0) then solarflux(iw)=ftmp
endfor
print, "nxsum", nxsum

S0 = total(solarflux)

print, "TOTAL SOLAR FLUX: ", total(solarflux) 
ScaleFac = Snorm/total(solarflux(*))
solarflux_out = solarflux(*) * ScaleFac
print, "TOTAL SOLAR FLUX: ", total(solarflux_out) 
S0 = total(solarflux_out)
print, "ScaleFac", ScaleFac

;-- write stellar spectra to netcdf ---
if (do_write eq 1) then begin
print, "writing to ",outname

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
NCDF_ATTPUT, id, varid4, "units", "W/m2 in each interval"

NCDF_CONTROL, id, /ENDEF

NCDF_VARPUT, id, varid1, rtwavlow
NCDF_VARPUT, id, varid2, rtwavhi
NCDF_VARPUT, id, varid3, S0
NCDF_VARPUT, id, varid4, solarflux_out

NCDF_CLOSE, id

endif


;-- plot stellar spectra, raw and binned
set_plot,'PS'
!P.Multi=[0,1,0]
loadct,40
!P.font=0
device,xsize=7.0,ysize=5.0,xoff=1.0,yoff=3.0,/inches
device,file='StellarSpectrum.ps'
device,/color,BITS=8            ;, /ENCAPSULATED, /CMYK 
;device,xsize=8.4,ysize=9,xoff=1.0,yoff=1.0,/CM   
device, set_font='Helvetica-Oblique', FONT_INDEX=20
device, set_font='Helvetica-Bold', FONT_INDEX=19
device, set_font='helvetica',FONT_INDEX=18

;create artifical barq=0
xbar1 = fltarr(nrtwavl*2)
ybar1 = fltarr(nrtwavl*2)
xbar2 = fltarr(nrtwavl*2)
ybar2 = fltarr(nrtwavl*2)
q=0
for i=0,nrtwavl-1 do begin
  xbar1(q) = 1.0e4/rtwavlow(i)
  xbar1(q+1) = 1.0e4/rtwavhi(i)
  ybar1(q) = solarflux_out(i)/rtwvldel(i)
  ybar1(q+1) = solarflux_out(i)/rtwvldel(i)
  q=q+2
endfor

loadct,40
plot, wlgth, sunm*ScaleFac, xtitle="!18wavenlength (microns)", xrange=[0.0,5.0], xstyle=1, $
             ytitle="!18Radiance (W m!U-2!N !M"+string("155B)+"!3m)"

;oplot, 1.0e4/rtwavmid, solarflux/rtwvldel, linestyle=90, psym=4,
;symsize=0.7                                                                      
;oplot, 1.0e4/rtwavmid, solarflux2/rtwvldel, linestyle=254, psym=4,
;symsize=0.7                                                                    
oplot, xbar1,ybar1, linestyle=0, thick=8.0, color=90

print, "plotting, StellarSpectrum.ps"
device, /close
set_plot,'X'


end
