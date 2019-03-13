pro makeSolarSpectrum
;! Creates solar spectral files for ExoRT
;! interval.  Data is read in from text files
;! 



do_write=1
Snorm=1360.0

spath = '/gpfsm/dnb02/projects/p54/users/etwolf/raw_stellar_spectra'
outname = "test_n63.nc"

;file = spath + '/hd128167um.txt' & npts = 1737 & outname = "F2V_hd128167.nc" ;F2V, 6595 K
;file = spath + '/hd22049um.txt' & npts = 1268 & outname = "K2V_hd22049.nc" ;K2V 5084 K
;file = spath + '/adleo_dat.txt' & npts = 4838 & header=strarr(175) & outname = "M4.5_adleo.nc" ;M4.5 3390K
;file = spath + '/gj644_dat.txt' & npts = 5483 & header=strarr(98) & outname = "M3.5v_gj644.nc"
;file = spath + '/hd114710um.txt' & npts = 1268 & outname = "G0V_hd114710.nc" ;G0V low activity 5860 K
;file = spath + '/hd206860um.txt' & npts = 1268 & outname = "G0V_hd206860.nc" ;G0V high activity 5957 K
;file = spath + '/lhs1140_bt-settl-interp_3216K_logg5_meta-0.24.txt' & npts = 389369 & outname = "LHS1140"

; used in Kopparapu et al. 2017
;file = spath + '/bt-settl_2600_logg4.5_FeH0.txt'  & npts = 390628
;file = spath + '/bt-settl_3000_logg4.5_FeH0.txt'  & npts = 394036
;file = spath + '/bt-settl_3300_logg4.5_FeH0.txt'  & npts = 396483
;file = spath + '/bt-settl_3700_logg4.5_FeH0.txt'  & npts = 398395
;file = spath + '/bt-settl_4000_logg4.5_FeH0.txt'  & npts = 398760
file = spath + '/bt-settl_4500_logg4.5_FeH0.txt'  & npts = 398558


filename=STRJOIN(STRSPLIT(file,/EXTRACT,' '))
print, filename


SunSolidAngle = 6.87e-5 ; as the Sun appears from Earth today


; define arrays
data=dblarr(2,npts)
wnm=dblarr(npts)
wlgth=dblarr(npts)
;dwm = 20.0
sunm=dblarr(npts)
dwm = dblarr(npts)

nhead = 9

header=strarr(nhead)

; read text data
OPENR,lun,filename, /GET_LUN
READF,lun,header
READF,lun, data
FREE_LUN,lun


for i=0, npts-1 do begin
  wlgth(i) = data(0,i)/1.0e4 ;microns
  wnm(i) = 1.0e4/wlgth(i)
  sunm(i) = data(1,i)
endfor

for i=0, npts-2 do begin
  dwm(i) = wnm(i)-wnm(i+1)
;  print, dwm
endfor

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
;nrtwavl = 68                                                             
;rtwavlow = dblarr(nrtwavl) 
;rtwavhi = dblarr(nrtwavl)  
;rtwavlow =  [1.0,   40.00000,   100.0000,   160.0000, $
;             220.0000,   280.0000,   330.0000,   380.0000, $
;             440.0000,   495.0000,   545.0000,   617.0000, $
;             667.0000,   720.0000,   800.0000,   875.0000, $
;             940.0000,   1000.000,   1065.000,   1108.000, $
;             1200.000,   1275.000,   1350.000,   1450.000, $
;             1550.000,   1650.000,   1750.000,   1850.000, $
;             1950.000,   2050.000,   2200.000,   2397.000, $
;             2494.000,   2796.000,   3087.000,   3425.000, $
;             3760.000,   4030.000,   4540.000,   4950.000, $
;             5370.000,   5925.000,   6390.000,   6990.000, $
;             7650.000,   8315.000,   8850.000,   9350.000, $
;             9650.000,   10400.00,   11220.00,   11870.00, $
;             12790.00,   13300.00,   14470.00,   15000.00, $
;             16000.00,   16528.00,   17649.00,   18198.00, $
;             18518.00,   22222.00,   25641.00,   29308.00, $
;             30376.00,   32562.00,   35087.00,   36363.00 ]
;rtwavhi =   [40.00000,   100.0000,   160.0000, $
;             220.0000,   280.0000,   330.0000,   380.0000, $
;             440.0000,   495.0000,   545.0000,   617.0000, $
;             667.0000,   720.0000,   800.0000,   875.0000, $
;             940.0000,   1000.000,   1065.000,   1108.000, $
;             1200.000,   1275.000,   1350.000,   1450.000, $
;             1550.000,   1650.000,   1750.000,   1850.000, $
;             1950.000,   2050.000,   2200.000,   2397.000, $
;             2494.000,   2796.000,   3087.000,   3425.000, $
;             3760.000,   4030.000,   4540.000,   4950.000, $
;             5370.000,   5925.000,   6390.000,   6990.000, $
;             7650.000,   8315.000,   8850.000,   9350.000, $
;             9650.000,   10400.00,   11220.00,   11870.00, $
;             12790.00,   13300.00,   14470.00,   15000.00, $
;             16000.00,   16528.00,   17649.00,   18198.00, $
;             18518.00,   22222.00,   25641.00,   29308.00, $
;             30376.00,   32562.00,   35087.00,   36363.00, $
;             42087.00 ]


; = rt = 63
nrtwavl = 63
rtwavlow = dblarr(nrtwavl)
rtwavhi = dblarr(nrtwavl)  
rtwavlow =  [10.,    100.,   160.,   220.,   $
             280.,   330.,   380.,   440.,   $
             495.,   545.,   617.,   667.,   $ 
             720.,   800.,   875.,   940.,   $ 
             1000.,  1065.,  1108.,  1200.,  $  
             1275.,  1350.,  1450.,  1550.,  $ 
             1650.,  1750.,  1850.,  1950.,  $
             2050.,  2200.,  2397.,  2494.,  $
             2796.,  3087.,  3425.,  3760.,  $ 
             4030.,  4540.,  4950.,  5370.,  $ 
             5925.,  6390.,  6990.,  7650.,  $
             8315.,  8850.,  9350.,  9650.,  $
             10400., 11220., 11870., 12750., $
             13300., 14470., 15000., 15385., $
             16667., 18182., 20000., 22222., $
             25000., 28571., 33333. ]    
rtwavhi =   [100.,   160.,   220.,   $
             280.,   330.,   380.,   440.,   $
             495.,   545.,   617.,   667.,   $ 
             720.,   800.,   875.,   940.,   $ 
             1000.,  1065.,  1108.,  1200.,  $  
             1275.,  1350.,  1450.,  1550.,  $ 
             1650.,  1750.,  1850.,  1950.,  $
             2050.,  2200.,  2397.,  2494.,  $
             2796.,  3087.,  3425.,  3760.,  $ 
             4030.,  4540.,  4950.,  5370.,  $ 
             5925.,  6390.,  6990.,  7650.,  $
             8315.,  8850.,  9350.,  9650.,  $
             10400., 11220., 11870., 12750., $
             13300., 14470., 15000., 15385., $
             16667., 18182., 20000., 22222., $
             25000., 28571., 33333., 40000.  ]




rtwavmid = (rtwavhi(*) + rtwavlow(*) )/2.0
rtwvlmid = (1.0e4/rtwavlow(*) + 1.0e4/rtwavhi(*))/2.0
rtwvldel = (1.0e4/rtwavlow(*) - 1.0e4/rtwavhi(*))
solarflux = dblarr(nrtwavl)

print, rtwavmid, rtwvldel

;new way                                                                                                                             
dwn = 1.0
dwni = 1

for iw=0, nrtwavl-1 do begin
  ftmp=0.
  wn_lw=rtwavlow(iw)
  wn_hi=rtwavhi(iw)
  for w=wn_lw,wn_hi-dwn,dwni do begin   ; integrate over each  spectral interval in steps of 10 cm-1                                             
    dwl=1.0e4/w - 1.0e4/(w+dwn)    
    ftmp=ftmp+interpol(sunm,wnm,w)*dwl
  endfor
  solarflux(iw)=ftmp
;  print, ftmp
endfor

S0 = total(solarflux)

print, "TOTAL SOLAR FLUX: ", total(solarflux) 
ScaleFac = Snorm/total(solarflux(*))
solarflux_out = solarflux(*) * ScaleFac
print, "TOTAL SOLAR FLUX: ", total(solarflux_out) 
S0 = total(solarflux_out)

if (do_write eq 1) then begin
print, "writing to ",outname
z
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
device,file='idl.ps'
device,/color,BITS=8            ;, /ENCAPSULATED, /CMYK                                                                                                       
;device,xsize=8.4,ysize=9,xoff=1.0,yoff=1.0,/CM                                                                                                    
device, set_font='Helvetica-Oblique', FONT_INDEX=20
device, set_font='Helvetica-Bold', FONT_INDEX=19
device, set_font='helvetica',FONT_INDEX=18

;create artifical bar                                                                                                                              
q=0
xbar1 = fltarr(nrtwavl*2)
ybar1 = fltarr(nrtwavl*2)
xbar2 = fltarr(nrtwavl*2)
ybar2 = fltarr(nrtwavl*2)
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
