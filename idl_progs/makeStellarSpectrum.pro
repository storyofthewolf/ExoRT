pro makeStellarSpectrum
;------------------------------------------------
;Author: Wolf, E.T.
;Created: August 2020
;Revision history: adapted from old code 
;                                                         
;DESCRIPTION:
; Creates solar spectral files for ExoRT.  
;  Reads data read in from text files and sorts.
;
;------------------------------------------------
;NOTES: 
; need to add a blackbody option
;



do_write=1
Snorm=1360.0

spath1 = '/gpfsm/dnb53/etwolf/models/ExoRT/data/solar/raw'
spath2 = '/gpfsm/dnb02/projects/p54/users/etwolf/raw_stellar_spectra/models_bt-settl_1577677702'

spath = spath2


;nhead = number of lines of header information in the raw stellar spectral file
;ncol = number of columns in the raw stellar spectral file
;npts = number of points/rows in the raw stellar spectral files

;file = spath1 + '/adleo_dat.txt' & npts = 4838 & header=strarr(175) & outname;= "M4.5_adleo_n68.nc" & nhead = 175 & ncol = 5 ;M4.5 3390K
file = spath1 + '/trappist1_sed.txt' & npts = 83302-20   & outname = 'trappist1_lincowski2018_n68.nc' & nhead=20 & ncol = 2
;file = spath1 + '/hd128167um.txt' & npts = 1737 & outname = "F2V_hd128167.nc" ;F2V, 6595 K
;file = spath1 + '/hd22049um.txt' & npts = 1268 & outname = "K2V_hd22049.nc" ;K2V 5084 K
;file = spath1 + '/gj644_dat.txt' & npts = 5483 & header=strarr(98) & outname = "M3.5v_gj644.nc"
;file = spath1 + '/hd114710um.txt' & npts = 1268 & outname = "G0V_hd114710.nc" ;G0V low activity 5860 K
;file = spath1 + '/hd206860um.txt' & npts = 1268 & outname = "G0V_hd206860.nc" ;G0V high activity 5957 K
;file = spath1 + '/lhs1140_bt-settl-interp_3216K_logg5_meta-0.24.txt' & npts = 389369 & outname = "LHS1140"

; used in Kopparapu et al. 2017
;file = spath + '/bt-settl_2600_logg4.5_FeH0.txt'  & npts = 390628  & nhead = 9 & ncol = 2
;file = spath + '/bt-settl_3000_logg4.5_FeH0.txt'  & npts = 394036  & nhead = 9 & ncol = 2
;file = spath + '/bt-settl_3300_logg4.5_FeH0.txt'  & npts = 396483  & nhead = 9 & ncol = 2
;file = spath + '/bt-settl_3700_logg4.5_FeH0.txt'  & npts = 398395  & nhead = 9 & ncol = 2
;file = spath + '/bt-settl_4000_logg4.5_FeH0.txt'  & npts = 398760  & nhead = 9 & ncol = 2
;file = spath + '/bt-settl_4500_logg4.5_FeH0.txt'  & npts = 398558  & nhead = 9 & ncol = 2
;file = spath + '/lte032-5.0-0.0a+0.0.BT-NextGen.7.dat.txt' & npts = 390737
;outname = "M4.5_AdLeo_n68.nc"


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
  
  ;wavelength in microns
  ;wlgth(i) = data(0,i)/1.0e4 ;microns
  wlgth(i) = data(0,i)  ; microns

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


!P.Multi=[0,1,0]
loadct,40
!P.font=0
set_plot,'PS'
device,xsize=7.0,ysize=5.0,xoff=1.0,yoff=3.0,/inches
device,file='idl.ps'
device,/color,BITS=8 ;, /ENCAPSULATED, /CMYK                                                  
;device,xsize=8.4,ysize=9,xoff=1.0,yoff=1.0,/CM                                               

device, set_font='Helvetica-Oblique', FONT_INDEX=20
device, set_font='Helvetica-Bold', FONT_INDEX=19
device, set_font='helvetica',FONT_INDEX=18


;sunm(*) = sunm(*) * SunSolidAngle 
plot, wlgth, sunm, xtitle = "!18wavenumber (cm!U-1!N)", xrange=[0,2.5]

device, /close
set_plot,'X'

;=  rt = 28
;nrtwavl = 28
;rtwavlow = dblarr(nrtwavl)
;rtwavhi = dblarr(nrtwavl)
;rtwavlow = [10., 350., 500., 630., 700., 820., 980., 1100., 1180., 1390., $
;            1480., 1800., 2080., 2200., 2380., 2600., 3250., 4000., 4650., $
;            5150., 6150., 7700., 8050., 12850., 16000., 22650., 29000., $
;            38000. ]

;rtwavhi = [350., 500., 630., 700., 820., 980., 1100., 1180., 1390., $
;            1480., 1800., 2080., 2200., 2380., 2600., 3250., 4000., 4650., $
;            5150., 6150., 7700., 8050., 12850., 16000., 22650., 29000., $
;            38000., 50000. ]


; = rt = 42
;nrtwavl = 42                                                             
;rtwavlow = dblarr(nrtwavl)                                               
;rtwavhi = dblarr(nrtwavl)                                                
;rtwavlow =  [  10.0,       200.0,         350.0,        425.0, $         
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
;rtwavhi =   [ 200.0,         350.0,        425.0, $                      
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



; = rt = 68
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


; = rt = 63
;nrtwavl = 63
;rtwavlow = dblarr(nrtwavl)
;rtwavhi = dblarr(nrtwavl)  
;rtwavlow =  [10.,    100.,   160.,   220.,   $
;             280.,   330.,   380.,   440.,   $
;             495.,   545.,   617.,   667.,   $ 
;             720.,   800.,   875.,   940.,   $ 
;             1000.,  1065.,  1108.,  1200.,  $  
;             1275.,  1350.,  1450.,  1550.,  $ 
;             1650.,  1750.,  1850.,  1950.,  $
;             2050.,  2200.,  2397.,  2494.,  $
;             2796.,  3087.,  3425.,  3760.,  $ 
;             4030.,  4540.,  4950.,  5370.,  $ 
;             5925.,  6390.,  6990.,  7650.,  $
;             8315.,  8850.,  9350.,  9650.,  $
;             10400., 11220., 11870., 12750., $
;             13300., 14470., 15000., 15385., $
;             16667., 18182., 20000., 22222., $
;             25000., 28571., 33333. ]    
;rtwavhi =   [100.,   160.,   220.,   $
;             280.,   330.,   380.,   440.,   $
;             495.,   545.,   617.,   667.,   $ 
;             720.,   800.,   875.,   940.,   $ 
;             1000.,  1065.,  1108.,  1200.,  $  
;             1275.,  1350.,  1450.,  1550.,  $ 
;             1650.,  1750.,  1850.,  1950.,  $
;             2050.,  2200.,  2397.,  2494.,  $
;             2796.,  3087.,  3425.,  3760.,  $ 
;             4030.,  4540.,  4950.,  5370.,  $ 
;             5925.,  6390.,  6990.,  7650.,  $
;             8315.,  8850.,  9350.,  9650.,  $
;             10400., 11220., 11870., 12750., $
;             13300., 14470., 15000., 15385., $
;             16667., 18182., 20000., 22222., $
;             25000., 28571., 33333., 40000.  ]

; = rt = 73
;nrtwavl = 73
;rtwavlow = dblarr(nrtwavl)
;rtwavhi = dblarr(nrtwavl)   
;rtwavlow  = [ 0.00000, 40.0000, 100.000, 160.000, 220.000, 280.000, 330.000, $
;            380.000, 440.000, 495.000, 545.000, 617.000, 667.000, 720.000, $
;            800.000, 875.000, 940.000, 1000.00, 1065.00, 1108.00, 1200.00, $
;            1275.00, 1350.00, 1450.00, 1550.00, 1650.00, 1750.00, 1850.00, $
;            1950.00, 2050.00, 2200.00, 2439.02, 2564.10, 2777.78, 3174.60, $
;            3508.77, 3773.59, 4081.63, 4545.46, 4716.98, 5154.64, 5376.34, $
;            5555.56, 5952.38, 6172.84, 6578.95, 6711.41, 6849.31, 7042.25, $
;            7462.69, 7692.31, 8064.52, 8333.33, 8620.69, 8928.57, 9090.91, $
;            9259.26, 9708.74, 10869.6, 11111.1, 11363.6, 11494.3, 12500.0, $
;            12820.5, 14492.8, 16393.4, 18181.8, 20000.0, 22222.2, 23809.5, $
;            25974.0, 28985.5, 33333.3, 50000.0  ]

;rtwavhi  = [ 0.00000, 40.0000, 100.000, 160.000, 220.000, 280.000, 330.000, $
;             380.000, 440.000, 495.000, 545.000, 617.000, 667.000, 720.000, $
;             800.000, 875.000, 940.000, 1000.00, 1065.00, 1108.00, 1200.00, $
;             1275.00, 1350.00, 1450.00, 1550.00, 1650.00, 1750.00, 1850.00, $
;             1950.00, 2050.00, 2200.00, 2439.02, 2564.10, 2777.78, 3174.60, $
;             3508.77, 3773.59, 4081.63, 4545.46, 4716.98, 5154.64, 5376.34, $
;             5555.56, 5952.38, 6172.84, 6578.95, 6711.41, 6849.31, 7042.25, $
;             7462.69, 7692.31, 8064.52, 8333.33, 8620.69, 8928.57, 9090.91, $
;             9259.26, 9708.74, 10869.6, 11111.1, 11363.6, 11494.3, 12500.0, $
;             12820.5, 14492.8, 16393.4, 18181.8, 20000.0, 22222.2, 23809.5, $
;             25974.0, 28985.5, 33333.3, 50000.0  ]



rtwavmid = (rtwavhi(*) + rtwavlow(*) )/2.0
rtwvlmid = (1.0e4/rtwavlow(*) + 1.0e4/rtwavhi(*))/2.0
rtwvldel = (1.0e4/rtwavlow(*) - 1.0e4/rtwavhi(*))
solarflux = dblarr(nrtwavl)
solarflux(*) = 0.0 

;print, rtwavmid, rtwvldel

;new way                                                                                                                             
dwn = 1.0
dwni = 1

for iw=0, nrtwavl-1 do begin
  ftmp=0.
  wn_lw=rtwavlow(iw)
  wn_hi=rtwavhi(iw)
  print, iw,  wn_lw, wn_hi
  for w=wn_lw,wn_hi-dwn,dwni do begin   ; integrate over each  spectral interval in steps of 10 cm-1   
    dwl=1.0e4/w - 1.0e4/(w+dwn)    

    ftmp=ftmp+interpol(sunm,wnm,w)*dwl
  endfor
  if(ftmp ge 0.0) then solarflux(iw)=ftmp

endfor

S0 = total(solarflux)

print, "TOTAL SOLAR FLUX: ", total(solarflux) 
ScaleFac = Snorm/total(solarflux(*))
solarflux_out = solarflux(*) * ScaleFac
print, "TOTAL SOLAR FLUX: ", total(solarflux_out) 
S0 = total(solarflux_out)

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
NCDF_ATTPUT, id, varid4, "units", "W/m2/"

NCDF_CONTROL, id, /ENDEF

NCDF_VARPUT, id, varid1, rtwavlow
NCDF_VARPUT, id, varid2, rtwavhi
NCDF_VARPUT, id, varid3, S0
NCDF_VARPUT, id, varid4, solarflux_out

NCDF_CLOSE, id

endif

!P.Multi=[0,1,0]
loadct,40
!P.font=0
set_plot,'PS'
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
plot, wlgth, sunm*ScaleFac, xtitle="!18wavenlength (micron)", xrange=[0.0,3.0], xstyle=1

;oplot, 1.0e4/rtwavmid, solarflux/rtwvldel, linestyle=90, psym=4,
;symsize=0.7                                                                      
;oplot, 1.0e4/rtwavmid, solarflux2/rtwvldel, linestyle=254, psym=4,
;symsize=0.7                                                                    
oplot, xbar1,ybar1, linestyle=0, thick=3.0, color=90

device, /close
set_plot,'X'


end
