pro makeCloudOptics
;--------------------------------------------------
;AUTHOR: Wolf E.T.
;
;DESCRIPTION:
; Produce optical property netcdf files for
; liquid cloud drops and ice cloud particles.
; Included are various spectral grids used by ExoCAM
; along with various plotting options.
; Choose wisely :-)
;---------------------------------------------------
;---------------------------------------------------

;;=============  options & files =================
do_h2o_liq = 0 ; flag, compute optical properties of H2O liquid cloud droplets
do_h2o_ice = 0 ; flag, compute optical properties of H2O ice particles
do_co2_ice = 1 ; flag, compute optical properties of CO2 ice particles

plot_ps = 1              ; if eq 0, then plot to x windows 
do_plot_refract = 1      ; plot raw refractive indices

do_mie = 1               ; do mie calculations and bin to spectral grid
do_plot_qwg = 1          ; plot extinction, single scattering, and asymmetry parameter  ;Currently not operable
do_write_netcdf = 1      ; flag to write netcdf outputs
 
h2o_liq_filename = 'cloudoptics_h2o_liquid_mie_n68.nc'
h2o_ice_filename = 'cloudoptics_h2o_ice_mie_n68.nc'
co2_ice_filename = 'cloudoptics_co2_ice_n68_r1000.nc'

if (do_h2o_liq eq 1) then print, "calculating H2O liquid clouds"
if (do_h2o_ice eq 1) then print, "calculating H2O ice clouds"
if (do_co2_ice eq 1) then print, "calculating CO2 ice clouds"

;-- choose one and only one --
;   -- spectral resolution 
do_n28 = 0
do_n35 = 0
do_n42 = 0
do_n68 = 1
do_n73 = 0
do_n84 = 0

colortable = 40 ; ;load standard color table

;-------- input file list ----------
root = '/discover/nobackup/etwolf/models'
H2O_LIQUID_REFRACT_FILE = root + '/ExoRT/data/cloud/refractive_indices/Segelstein_h2o_liquid.txt'
H2O_ICE_REFRACT_FILE    = root + '/ExoRT/data/cloud/refractive_indices/Warren_h2o_ice.txt'
CO2_ICE_REFRACT_FILE    = root + '/ExoRT/data/cloud/refractive_indices/Hansen_co2_ice.txt'



;;=============  read in refractive indice files =================
;-------- read h2o liquid files ----------
if (do_h2o_liq eq 1) then begin

  nw_h2o_liq = 1261
  header_h2o_liq = strarr(4)
  data_h2o_liq = fltarr(3,nw_h2o_liq)

  wavelength_h2o_liq = fltarr(nw_h2o_liq)
  real_h2o_liq = fltarr(nw_h2o_liq)
  imaginary_h2o_liq = fltarr(nw_h2o_liq)

  OPENR,lun,H2O_LIQUID_REFRACT_FILE,/Get_lun
  READF,lun,header_h2o_liq
  READF,lun,data_h2o_liq

  wavelength_h2o_liq(*) = data_h2o_liq(0,*)
  real_h2o_liq(*) = data_h2o_liq(1,*)
  imaginary_h2o_liq(*) = data_h2o_liq(2,*)

  FREE_LUN,lun

  ;sort out where size parameter too big, only relevant for far UV
  ;wl = where(2*!pi*30.0/wavelength_liq lt 12000.0, nwliq)
  wl = where(2*!pi*100.0/wavelength_h2o_liq lt 12000.0, nw_h2o_liq)
  wavelength_h2o_liq = wavelength_h2o_liq(wl)
  real_h2o_liq = real_h2o_liq(wl)
  imaginary_h2o_liq = imaginary_h2o_liq(wl)

endif

;-------- read h2o ice files ----------
if (do_h2o_ice eq 1) then begin
  
  nw_h2o_ice = 486
  data_h2o_ice = fltarr(3,nw_h2o_ice)

  wavelength_h2o_ice = fltarr(nw_h2o_ice)
  real_h2o_ice = fltarr(nw_h2o_ice)
  imaginary_h2o_ice = fltarr(nw_h2o_ice)

  OPENR,lun,H2O_ICE_REFRACT_FILE,/Get_lun
  READF,lun,data_h2o_ice

  wavelength_h2o_ice(*) = data_h2o_ice(0,*)
  real_h2o_ice(*) = data_h2o_ice(1,*)
  imaginary_h2o_ice(*) = data_h2o_ice(2,*)

  FREE_LUN,lun

  ;sort out where size parameter too big, only relevant for far UV
  wi = where(2*!pi*300.0/wavelength_h2o_ice lt 12000.0, nw_h2o_ice)
  wavelength_h2o_ice = wavelength_h2o_ice(wi)
  real_h2o_ice = real_h2o_ice(wi)
  imaginary_h2o_ice = imaginary_h2o_ice(wi)

endif

;-------- read co2 ice files ----------
if (do_co2_ice eq 1) then begin
  
  nw_co2_ice = 15355
  header_co2_ice = strarr(5)
  data_co2_ice = fltarr(3,nw_co2_ice)

  wavelength_co2_ice = fltarr(nw_co2_ice)
  real_co2_ice = fltarr(nw_co2_ice)
  imaginary_co2_ice = fltarr(nw_co2_ice)

  OPENR,lun,CO2_ICE_REFRACT_FILE,/Get_lun
  READF,lun,header_co2_ice
  READF,lun,data_co2_ice

  wavelength_co2_ice(*) = data_co2_ice(0,*)
  real_co2_ice(*) = data_co2_ice(1,*)
  imaginary_co2_ice(*) = data_co2_ice(2,*)

  FREE_LUN,lun

  ;sort out where size parameter too big, only relevant for far UV
  wi = where(2*!pi*300.0/wavelength_co2_ice lt 12000.0, nw_co2_ice)
  wavelength_co2_ice = wavelength_co2_ice(wi)
  real_co2_ice = real_co2_ice(wi)
  imaginary_co2_ice = imaginary_co2_ice(wi)

endif



;;=============  plot refractive indices ================= 
if (do_plot_refract eq 1) then  begin

  ;--------------------------
  ; plot real refactive index
  ;--------------------------
  if (plot_ps eq 1) then begin
    print, "plotting to postscript: idl_cld_refrac_real.eps"
    set_plot,'PS'
    device,file='idl_cld_refrac_real.eps'
    device,/color,BITS=8, /ENCAPSULATED ;, /CMYK
    device,xsize=18.5,ysize=11,/CM
    device, set_font='Helvetica-Oblique', FONT_INDEX=20
    device, set_font='Helvetica-Bold', FONT_INDEX=19
    device, set_font='helvetica',FONT_INDEX=18
    loadct, colortable
  endif else begin
    print, "plotting to x-windows"
    set_plot,'x'
  endelse

  xr=[0.1,100]
  yr=[0,4]

  xplaceholder=findgen(10) & yplaceholder=findgen(10)
  plot, xplaceholder, yplaceholder,/nodata, $
         xrange=xr, xstyle=1, yrange=yr, ystyle=1, xthick=3.0, ythick=3.0, /xlog, $
         xtitle="wavelength (microns)", $
         ytitle="refractive indices (n)" ;, $                    

  if (do_h2o_ice eq 1) then begin
    oplot, wavelength_h2o_ice, real_h2o_ice, color=90, thick=2, linestyle=0
  endif

  if (do_h2o_liq eq 1) then begin
    oplot, wavelength_h2o_liq, real_h2o_liq, color=250, thick=2, linestyle=0
   endif

  if (do_co2_ice eq 1) then begin
    oplot, wavelength_co2_ice, real_co2_ice, color=150, thick=2, linestyle=0
  endif

  if  (plot_ps eq 1) then begin 
    device, /close
  endif else begin
    print, "code stopped for plot viewing, press .c to continue"
    stop
  endelse

  ;--------------------------
  ; plot imaginary refactive index
  ;--------------------------
  if (plot_ps eq 1) then begin
    print, "plotting to postscript: idl_cld_refrac_imaginary.eps"
    set_plot,'PS'
    device,file='idl_cld_refrac_imaginary.eps'
    device,/color,BITS=8, /ENCAPSULATED ;, /CMYK
    device,xsize=18.5,ysize=11,/CM
    device, set_font='Helvetica-Oblique', FONT_INDEX=20
    device, set_font='Helvetica-Bold', FONT_INDEX=19
    device, set_font='helvetica',FONT_INDEX=18
    loadct, colortable
  endif else begin
    print, "plotting to x-windows"
    set_plot,'x'
  endelse

  xr=[0.1,100]
  yr=[1.0e-10,1.0e4]

  xplaceholder=findgen(10) & yplaceholder=findgen(10)
  plot, xplaceholder, yplaceholder,/nodata, $
         xrange=xr, xstyle=1, yrange=yr, ystyle=1, xthick=3.0, ythick=3.0, /xlog, /ylog, $
         xtitle="wavelength (microns)", $
         ytitle="refractive indices (i)" ;, $                    

  if (do_h2o_ice eq 1) then begin
    oplot, wavelength_h2o_ice, imaginary_h2o_ice, color=90, thick=2, linestyle=2
  endif

  if (do_h2o_liq eq 1) then begin
    oplot, wavelength_h2o_liq, imaginary_h2o_liq, color=250, thick=2, linestyle=2
   endif

  if (do_co2_ice eq 1) then begin
    oplot, wavelength_co2_ice, imaginary_co2_ice, color=150, thick=2, linestyle=2
 endif

  if  (plot_ps eq 1) then begin 
    device, /close
  endif else begin
    print, "code stopped for plot viewing, press .c to continue"
    stop
  endelse


endif



;;=============  Definitions of spectral interval grid =================
if (do_mie eq 1) then begin
;--------------------------

;wavenumber interval edges [cm-1]
;28 bin
if (do_n28 eq 1) then begin
  wavenum_edge = [ 10., 350., 500., 630., 700., 820., 980., 1080., 1180., 1390., $
                   1480., 1800., 2080., 2250., 2380., 2600., 3250., 4000., 4650., $
                   5150., 6150., 7700., 8050., 12850., 16000., 22650., 29000., $
                   38000., 50000. ]
endif

;42 bin
if (do_n42 eq 1) then begin
  wavenum_edge  = [      10.0,       200.0,         350.0,        425.0, $
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

;n68
if (do_n68 eq 1) then begin
  wavenum_edge =  [0.00E+00,   40.00000,   100.0000,   160.0000, $
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

;n73 
if (do_n73 eq 1) then begin
  wavenum_edge  = [ 0.00000, 40.0000, 100.000, 160.000, 220.000, 280.000, 330.000,$
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

;n84 (shortwave extension of the n68)
if (do_n84 eq 1) then begin
  wavenum_edge =  [0.00E+00,   40.00000,   100.0000,   160.0000, $
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
                   400000.00, 500000.00, 750000.00,  1000000.00, $
                  1250000.00 ]
endif

nwgrid = n_elements(wavenum_edge)
wavelength_edge = fltarr(nwgrid)
wavelength_mid = fltarr(nwgrid-1)

;convert wavenum_edge to wavelength microns
wavelength_edge(*) = 1.0e4/wavenum_edge(*)

for nn=0, nwgrid-2 do begin
  wavelength_mid(nn) = (wavelength_edge(nn)+wavelength_edge(nn+1))/2.
endfor



;;=============  Define effective radii grids =================

;effective radii for h2o liquid drop cloud particles
rel_h2o_grid = [ 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, $
             11, 12, 13, 14, 15, 16, 17, 18, 19, 20, $
             21, 22, 23, 24, 25, 26, 27, 28, 29, 30 ]
rel_h2o_grid = findgen(30) + 1
rel_h2o_grid = [14]  ; single value array, uncomment for fast plotting only 
nrel_h2o = n_elements(rel_h2o_grid)

;effective radiii for h2o ice particles
rei_h2o_grid = findgen(300) + 1
rei_h2o_grid = [100]  ; single value array, uncomment for fast plotting only 
nrei_h2o = n_elements(rei_h2o_grid)
dge_grid = fltarr(nrei_h2o)   ; generalized effective size ice particles

;effective radiii for co2 ice particles
rei_co2_grid = [1, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 60, 70, 80, 90, 100, 125, 150, 175, 200, 250, 300, 350, 400, 450, 500, 600, 700, 800, 900, 1000]
;rei_co2_grid = findgen(1000) + 1
;rei_co2_grid = [1]  ; single value array, uncomment for fast plotting only 
nrei_co2 = n_elements(rei_co2_grid)

;optical constants averaged over spectral intervals
Qext_h2o_liq_out = fltarr(nrel_h2o, nwgrid-1)
W_h2o_liq_out = fltarr(nrel_h2o, nwgrid-1)
G_h2o_liq_out = fltarr(nrel_h2o, nwgrid-1)

Qext_h2o_ice_out = fltarr(nrei_h2o, nwgrid-1)
W_h2o_ice_out = fltarr(nrei_h2o, nwgrid-1)
G_h2o_ice_out = fltarr(nrei_h2o, nwgrid-1)

Qext_co2_ice_out = fltarr(nrei_co2, nwgrid-1)
W_co2_ice_out = fltarr(nrei_co2, nwgrid-1)
G_co2_ice_out = fltarr(nrei_co2, nwgrid-1)



;;=============  Mie computations start here =================  

if (do_h2o_liq eq 1) then begin   

  ;h2o liquid droplet optical constants on Segelstein grid
  Qext_h2o_liq = fltarr(nrel_h2o, nw_h2o_liq)
  W_h2o_liq = fltarr(nrel_h2o, nw_h2o_liq)
  G_h2o_liq = fltarr(nrel_h2o, nw_h2o_liq)

  ; calculate optical properties of liquid droplets on Segelstein grid
  for i = 0, nrel_h2o-1 do begin
    print, i, " size=", rel_h2o_grid(i)
    for j = 0, nw_h2o_liq-1 do begin
      cm = complex(real_h2o_liq(j), -imaginary_h2o_liq(j))
      alpha = 2.0 * !pi * rel_h2o_grid(i) / wavelength_h2o_liq(j)   
      mie_single, alpha, cm, dqext, dqscat, dqbk, dg
      Qext_h2o_liq(i,j) = dqext
      W_h2o_liq(i,j) = dqscat/dqext
      G_h2o_liq(i,j) = dg 
    endfor
  endfor

  ; interpolate/extrapolate the optical properties onto a finer grid 
  ; based on the selected spectral grid.  For n68 sometimes no
  ; optical properties exist in a band.  n84 extends beyond data
  n_newgrid = (wavenum_edge(n_elements(wavenum_edge)-1)-wavenum_edge(0))+1
  ;print, n_newgrid, wavenum_edge(n_elements(wavenum_edge)-1), wavenum_edge(0)
  newgrid = findgen(n_newgrid+1)+wavenum_edge(0)
  newgrid = 1.0e4/newgrid
  Qext_h2o_liq_intex = fltarr(nrel_h2o, n_newgrid+1)
  W_h2o_liq_intex = fltarr(nrel_h2o, n_newgrid+1)
  G_h2o_liq_intex = fltarr(nrel_h2o, n_newgrid+1)
  for i = 0, nrel_h2o-1 do begin
    Qext_h2o_liq_intex(i,*) = interpol(Qext_h2o_liq(i,*), wavelength_h2o_liq(*),  newgrid)
    W_h2o_liq_intex(i,*) = interpol(W_h2o_liq(i,*), wavelength_h2o_liq(*),  newgrid)
    G_h2o_liq_intex(i,*) = interpol(G_h2o_liq(i,*), wavelength_h2o_liq(*),  newgrid)
  endfor

  ; average optical properties over all Segelstein points in spectral interval
  for i = 0, nrel_h2o-1 do begin
    for j = 0, nwgrid-2 do begin
      ; orig
      qx=where(wavelength_h2o_liq le wavelength_edge(j) and wavelength_h2o_liq ge wavelength_edge(j+1), qcount)
      if (qcount gt 0) then begin ; raw refractive indices exist in the bin, use them
        ;print, j, qcount, MEAN(Qext_h2o_liq(i,qx)), MEAN(W_h2o_liq(i,qx)), MEAN(G_h2o_liq(i,qx))
        Qext_h2o_liq_out(i,j) = MEAN(Qext_h2o_liq(i,qx))
        W_h2o_liq_out(i,j)    = MEAN(W_h2o_liq(i,qx))
        G_h2o_liq_out(i,j)    = MEAN(G_h2o_liq(i,qx))
      endif else begin  ; no raw refractive indices exist in the bin, interpolate/extrapolate
        qx=where(newgrid le wavelength_edge(j) and newgrid ge wavelength_edge(j+1), qcount)
        ;print, j, qcount, MEAN(Qext_h2o_liq_intex(i,qx)), MEAN(W_h2o_liq_intex(i,qx)), MEAN(G_h2o_liq_intex(i,qx))
        Qext_h2o_liq_out(i,j) = MEAN(Qext_h2o_liq_intex(i,qx))
        W_h2o_liq_out(i,j)    = MEAN(W_h2o_liq_intex(i,qx))
        G_h2o_liq_out(i,j)    = MEAN(G_h2o_liq_intex(i,qx))
      endelse
    endfor
  endfor

  ;print, "qext_h2o_liq",Qext_h2o_liq_out(0,0), Qext_h2o_liq_out(1,0), Qext_h2o_liq_out(2,0)
  ;print, "W_h2o_liq",W_h2o_liq_out(0,0), W_h2o_liq_out(1,0), W_h2o_liq_out(2,0)
  ;print, "G_h2o_liq",G_h2o_liq_out(0,0), G_h2o_liq_out(1,0), G_h2o_liq_out(2,0)

  if (do_write_netcdf) then begin
    ;write to netcdf file
    id = NCDF_CREATE(h2o_liq_filename, /CLOBBER)
    dim1 = NCDF_DIMDEF(id, 'rel_bins', nrel_h2o)
    dim2 = NCDF_DIMDEF(id, 'nwavlrng', nwgrid-1)
    dim3 = NCDF_DIMDEF(id, 'nwave_edge',nwgrid) 

    varid1 = NCDF_VARDEF(id, 'rel_grid', [dim1], /float)
    varid2 = NCDF_VARDEF(id, 'wvnrng', [dim3], /float)
    varid3 = NCDF_VARDEF(id, 'Qext_liq',[dim1,dim2], /float)
    varid4 = NCDF_VARDEF(id, 'W_liq', [dim1,dim2], /float)
    varid5 = NCDF_VARDEF(id, 'G_liq', [dim1,dim2], /float)

    NCDF_ATTPUT, id, varid1, "title", "effective radii grid, liquid cloud water drops"
    NCDF_ATTPUT, id, varid1, "units", "microns"

    NCDF_ATTPUT, id, varid2, "title", "wavenumber at edges" 
    NCDF_ATTPUT, id, varid2, "units", "cm-1"

    NCDF_ATTPUT, id, varid3, "title", "extinction efficiency"
    NCDF_ATTPUT, id, varid3, "units", "unitless"

    NCDF_ATTPUT, id, varid4, "title", "single scattering albedo"
    NCDF_ATTPUT, id, varid4, "units", "unitless"

    NCDF_ATTPUT, id, varid5, "title", "asymmetry parameter"
    NCDF_ATTPUT, id, varid5, "units", "unitless"
  
    NCDF_CONTROL, id, /ENDEF

    NCDF_VARPUT, id, varid1, rel_h2o_grid
    NCDF_VARPUT, id, varid2, wavenum_edge
    NCDF_VARPUT, id, varid3, Qext_h2o_liq_out
    NCDF_VARPUT, id, varid4, W_h2o_liq_out
    NCDF_VARPUT, id, varid5, G_h2o_liq_out

    NCDF_CLOSE, id

    print, "wrote liquid cloud optical constants to, ",h2o_liq_filename
 endif ;do_write_netcdf
endif                   ; do_liq

;----------------------------------------------------------------------
if (do_h2o_ice eq 1) then begin   

  ;h2o ice optical constants on Warren grid
  Qext_h2o_ice = fltarr(nrei_h2o, nw_h2o_ice)
  W_h2o_ice = fltarr(nrei_h2o, nw_h2o_ice)
  G_h2o_ice = fltarr(nrei_h2o, nw_h2o_ice)

   ; calculate optical properties for h2o ice particles on Warren grid
   for i = 0, nrei_h2o-1 do begin
    print, i, " size=", rei_h2o_grid(i)
    for j = 0, nw_h2o_ice-1 do begin
      cm = complex(real_h2o_ice(j), -imaginary_h2o_ice(j))
      alpha = 2.0 * !pi * rei_h2o_grid(i) / wavelength_h2o_ice(j)   
      mie_single, alpha, cm, dqext, dqscat, dqbk, dg
      Qext_h2o_ice(i,j) = dqext
      W_h2o_ice(i,j) = dqscat/dqext
      G_h2o_ice(i,j) = dg 
    endfor       
  endfor

  ; interpolate/extrapolate the optical properties onto a finer grid 
  ; based on the selected spectral grid.  For n68 sometimes no
  ; optical properties exist in a band.  n84 extends beyond data
  n_newgrid = (wavenum_edge(n_elements(wavenum_edge)-1)-wavenum_edge(0))+1
  ;print, n_newgrid, wavenum_edge(n_elements(wavenum_edge)-1), wavenum_edge(0)
  newgrid = findgen(n_newgrid+1)+wavenum_edge(0)
  newgrid = 1.0e4/newgrid
  Qext_h2o_ice_intex = fltarr(nrei_h2o, n_newgrid+1)
  W_h2o_ice_intex = fltarr(nrei_h2o, n_newgrid+1)
  G_h2o_ice_intex = fltarr(nrei_h2o, n_newgrid+1)
  for i = 0, nrei_h2o-1 do begin
    Qext_h2o_ice_intex(i,*) = interpol(Qext_h2o_ice(i,*), wavelength_h2o_ice(*),  newgrid)
    W_h2o_ice_intex(i,*) = interpol(W_h2o_ice(i,*), wavelength_h2o_ice(*),  newgrid)
    G_h2o_ice_intex(i,*) = interpol(G_h2o_ice(i,*), wavelength_h2o_ice(*),  newgrid)
  endfor

  ; average optical properties over all Warren points in spectral interval
  for i = 0, nrei_h2o-1 do begin
     for j = 0, nwgrid-2 do begin
      qx=where(wavelength_h2o_ice le wavelength_edge(j) and wavelength_h2o_ice ge wavelength_edge(j+1), qcount)
      if (qcount gt 0) then begin ; raw refractive indices exist in the bin, use them
        ;print, j, qcount, MEAN(Qext_h2o_ice(i,qx)), MEAN(W_h2o_ice(i,qx)), MEAN(G_h2o_ice(i,qx))
        Qext_h2o_ice_out(i,j) = MEAN(Qext_h2o_ice(i,qx))
        W_h2o_ice_out(i,j)    = MEAN(W_h2o_ice(i,qx))
        G_h2o_ice_out(i,j)    = MEAN(G_h2o_ice(i,qx))
      endif else begin  ; no raw refractive indices exist in the bin, interpolate/extrapolate 
        qx=where(newgrid le wavelength_edge(j) and newgrid ge wavelength_edge(j+1), qcount)
        ;print, j, qcount, MEAN(Qext_h2o_ice_intex(i,qx)), MEAN(W_h2o_ice_intex(i,qx)), MEAN(G_h2o_ice_intex(i,qx))
        Qext_h2o_ice_out(i,j) = MEAN(Qext_h2o_ice_intex(i,qx))
        W_h2o_ice_out(i,j)    = MEAN(W_h2o_ice_intex(i,qx))
        G_h2o_ice_out(i,j)    = MEAN(G_h2o_ice_intex(i,qx))
      endelse
    endfor
  endfor

;  print, "qext_h2o_ice",Qext_h2o_ice_out(0,0), Qext_h2o_ice_out(1,0), Qext_h2o_ice_out(2,0)
;  print, "W_h2o_ice",W_h2o_ice_out(0,0), W_h2o_ice_out(1,0), W_h2o_ice_out(2,0)
;  print, "G_h2o_ice",G_h2o_ice_out(0,0), G_h2o_ice_out(1,0), G_h2o_ice_out(2,0)

  if (do_write_netcdf eq 1) then begin
    ;write to netcdf file
    id = NCDF_CREATE(h2o_ice_filename, /CLOBBER)
    dim1 = NCDF_DIMDEF(id, 'rei_bins', nrei_h2o)
    dim2 = NCDF_DIMDEF(id, 'nwavlrng', nwgrid-1)
    dim3 = NCDF_DIMDEF(id, 'nwave_edge',nwgrid) 

    varid1 = NCDF_VARDEF(id, 'rei_grid', [dim1], /float)
    varid2 = NCDF_VARDEF(id, 'wvnrng', [dim3], /float)
    varid3 = NCDF_VARDEF(id, 'Qext_ice',[dim1,dim2], /float)
    varid4 = NCDF_VARDEF(id, 'W_ice', [dim1,dim2], /float)
    varid5 = NCDF_VARDEF(id, 'G_ice', [dim1,dim2], /float)

    NCDF_ATTPUT, id, varid1, "title", "effective radii grid, ice cloud particles"
    NCDF_ATTPUT, id, varid1, "units", "microns"

    NCDF_ATTPUT, id, varid2, "title", "wavenumber at edges"
    NCDF_ATTPUT, id, varid2, "units", "cm-1"

    NCDF_ATTPUT, id, varid3, "title", "extinction efficiency"
    NCDF_ATTPUT, id, varid3, "units", "unitless"

    NCDF_ATTPUT, id, varid4, "title", "single scattering albedo"
    NCDF_ATTPUT, id, varid4, "units", "unitless"

    NCDF_ATTPUT, id, varid5, "title", "asymmetry parameter"
    NCDF_ATTPUT, id, varid5, "units", "unitless"
  
    NCDF_CONTROL, id, /ENDEF

    NCDF_VARPUT, id, varid1, rei_h2o_grid
    NCDF_VARPUT, id, varid2, wavenum_edge
    NCDF_VARPUT, id, varid3, Qext_h2o_ice_out
    NCDF_VARPUT, id, varid4, W_h2o_ice_out
    NCDF_VARPUT, id, varid5, G_h2o_ice_out

    NCDF_CLOSE, id

    print, "wrote ice cloud optical constants to, ",h2o_ice_filename
 endif ; do_write_netcdf

endif  ; do_h2o_ice

;----------------------------------------------------------------------
if (do_co2_ice eq 1) then begin   

  ;CO2 ice optical constants on Hansen grid
  Qext_co2_ice = fltarr(nrei_co2, nw_co2_ice)
  W_co2_ice = fltarr(nrei_co2, nw_co2_ice)
  G_co2_ice = fltarr(nrei_co2, nw_co2_ice)

   ; calculate optical properties for co2 ice particles on Hansen grid
   for i = 0, nrei_co2-1 do begin
    print, i, " size=", rei_co2_grid(i)
    for j = 0, nw_co2_ice-1 do begin
      cm = complex(real_co2_ice(j), -imaginary_co2_ice(j))
      alpha = 2.0 * !pi * rei_co2_grid(i) / wavelength_co2_ice(j)   
      mie_single, alpha, cm, dqext, dqscat, dqbk, dg
      Qext_co2_ice(i,j) = dqext
      W_co2_ice(i,j) = dqscat/dqext
      G_co2_ice(i,j) = dg 
    endfor       
  endfor

  ; interpolate/extrapolate the optical properties onto a finer grid 
  ; based on the selected spectral grid.  For n68 sometimes no
  ; optical properties exist in a band.  n84 extends beyond data
  n_newgrid = (wavenum_edge(n_elements(wavenum_edge)-1)-wavenum_edge(0))+1
  ;print, n_newgrid, wavenum_edge(n_elements(wavenum_edge)-1), wavenum_edge(0)
  newgrid = findgen(n_newgrid+1)+wavenum_edge(0)
  newgrid = 1.0e4/newgrid
  Qext_co2_ice_intex = fltarr(nrei_co2, n_newgrid+1)
  W_co2_ice_intex = fltarr(nrei_co2, n_newgrid+1)
  G_co2_ice_intex = fltarr(nrei_co2, n_newgrid+1)
  for i = 0, nrei_co2-1 do begin
    Qext_co2_ice_intex(i,*) = interpol(Qext_co2_ice(i,*), wavelength_co2_ice(*),  newgrid)
    W_co2_ice_intex(i,*) = interpol(W_co2_ice(i,*), wavelength_co2_ice(*),  newgrid)
    G_co2_ice_intex(i,*) = interpol(G_co2_ice(i,*), wavelength_co2_ice(*),  newgrid)
  endfor

  ; average optical properties over all Hansen points in spectral interval
  for i = 0, nrei_co2-1 do begin
    for j = 0, nwgrid-2 do begin
      qx=where(wavelength_co2_ice le wavelength_edge(j) and wavelength_co2_ice ge wavelength_edge(j+1), qcount)
      if (qcount gt 0) then begin ; raw refractive indices exist in the bin, use them
        ;print, j, qcount, MEAN(Qext_co2_ice(i,qx)), MEAN(W_co2_ice(i,qx)), MEAN(G_co2_ice(i,qx))
        Qext_co2_ice_out(i,j) = MEAN(Qext_co2_ice(i,qx))
        W_co2_ice_out(i,j)    = MEAN(W_co2_ice(i,qx))
        G_co2_ice_out(i,j)    = MEAN(G_co2_ice(i,qx))
      endif else begin  ; no raw refractive indices exist in the bin, interpolate/extrapolate
        qx=where(newgrid le wavelength_edge(j) and newgrid ge wavelength_edge(j+1), qcount)
        ;print, j, qcount, MEAN(Qext_co2_ice_intex(i,qx)), MEAN(W_co2_ice_intex(i,qx)), MEAN(G_co2_ice_intex(i,qx))
        Qext_co2_ice_out(i,j) = MEAN(Qext_co2_ice_intex(i,qx))
        W_co2_ice_out(i,j)    = MEAN(W_co2_ice_intex(i,qx))
        G_co2_ice_out(i,j)    = MEAN(G_co2_ice_intex(i,qx))
      endelse
    endfor
  endfor

;  print, "qext_co2_ice",Qext_co2_ice_out(0,0), Qext_co2_ice_out(1,0), Qext_co2_ice_out(2,0)
;  print, "W_co2_ice",W_co2_ice_out(0,0), W_co2_ice_out(1,0), W_co2_ice_out(2,0)
;  print, "G_co2_ice",G_co2_ice_out(0,0), G_co2_ice_out(1,0), G_co2_ice_out(2,0)

  if (do_write_netcdf eq 1) then begin
    ;write to netcdf file
    id = NCDF_CREATE(co2_ice_filename, /CLOBBER)
    dim1 = NCDF_DIMDEF(id, 'nbins', nrei_co2)
    dim2 = NCDF_DIMDEF(id, 'nwavlrng', nwgrid-1)
    dim3 = NCDF_DIMDEF(id, 'nwave_edge',nwgrid) 

    varid1 = NCDF_VARDEF(id, 'radii', [dim1], /float)
    varid2 = NCDF_VARDEF(id, 'wvnrng', [dim3], /float)
    varid3 = NCDF_VARDEF(id, 'Qext',[dim1,dim2], /float)
    varid4 = NCDF_VARDEF(id, 'W', [dim1,dim2], /float)
    varid5 = NCDF_VARDEF(id, 'G', [dim1,dim2], /float)

    NCDF_ATTPUT, id, varid1, "title", "effective radii grid, CO2 ice cloud particles"
    NCDF_ATTPUT, id, varid1, "units", "microns"

    NCDF_ATTPUT, id, varid2, "title", "wavenumber at edges"
    NCDF_ATTPUT, id, varid2, "units", "cm-1"

    NCDF_ATTPUT, id, varid3, "title", "extinction efficiency"
    NCDF_ATTPUT, id, varid3, "units", "unitless"

    NCDF_ATTPUT, id, varid4, "title", "single scattering albedo"
    NCDF_ATTPUT, id, varid4, "units", "unitless"

    NCDF_ATTPUT, id, varid5, "title", "asymmetry parameter"
    NCDF_ATTPUT, id, varid5, "units", "unitless"
  
    NCDF_CONTROL, id, /ENDEF

    NCDF_VARPUT, id, varid1, rei_co2_grid
    NCDF_VARPUT, id, varid2, wavenum_edge
    NCDF_VARPUT, id, varid3, Qext_co2_ice_out
    NCDF_VARPUT, id, varid4, W_co2_ice_out
    NCDF_VARPUT, id, varid5, G_co2_ice_out

    NCDF_CLOSE, id

    print, "wrote CO2 ice cloud optical constants to, ", co2_ice_filename
 endif ; do_write_netcdf

endif  ; do_co2_ice


;;=============  plotting of Q, W, and G from mie computations  =================  

if (do_plot_qwg eq 1) then  begin

  ; select index of radii to plot data
  rei_h2o_select = 0 ;99
  rel_h2o_select = 0 ;13
  rei_co2_select = 0 ;10
  

  ;-------  plot QEXT ----------
  if (plot_ps eq 1) then begin
    print, "plotting to postscript: idl_cld_q.eps"
    set_plot,'PS'
    device,file='idl_cld_q.eps'
    device,/color,BITS=8, /ENCAPSULATED ;, /CMYK
    device,xsize=18.5,ysize=11,/CM
    device, set_font='Helvetica-Oblique', FONT_INDEX=20
    device, set_font='Helvetica-Bold', FONT_INDEX=19
    device, set_font='helvetica',FONT_INDEX=18
    loadct, colortable
   endif else begin
    print, "plotting to x-windows"
    set_plot,'x'
   endelse

   xr=[0.01,1000]
   yr=[0,4]

   xplaceholder=findgen(10) & yplaceholder=findgen(10)
   plot, xplaceholder, yplaceholder,/nodata, $
         xrange=xr, xstyle=1, yrange=yr, ystyle=1, xthick=3.0, ythick=3.0, /xlog, $
         xtitle="wavelength (microns)", $
         ytitle="Qext" ;, $                    

   if (do_h2o_ice eq 1) then begin
     print, "rei_h2o: ", rei_h2o_grid(rei_h2o_select)
     oplot, wavelength_mid, Qext_h2o_ice_out(rei_h2o_select,*), color=240, psym=1, symsize=1, thick=2
     oplot, wavelength_h2o_ice, Qext_h2o_ice(rei_h2o_select,*), color=90, thick=2, linestyle=0
     xyouts, 0.65, 0.95, "rei = "+string(rei_h2o_grid(rei_h2o_select))+" microns", /normal, size=1
   endif

   if (do_h2o_liq eq 1) then begin
     print, "rel_h2o: ", rel_h2o_grid(rel_h2o_select)
     oplot, wavelength_mid, Qext_h2o_liq_out(rel_h2o_select,*), color=240, psym=1, symsize=1, thick=2
     oplot, wavelength_h2o_liq, Qext_h2o_liq(rel_h2o_select,*), color=90, thick=2, linestyle=0
     xyouts, 0.65, 0.95, "rel = "+string(rel_h2o_grid(rel_h2o_select))+" microns", /normal, size=1
   endif

   if (do_co2_ice eq 1) then begin
     print, "rei_co2: ", rei_co2_grid(rei_co2_select)
     oplot, wavelength_mid, Qext_co2_ice_out(rei_co2_select,*), color=240, psym=1, symsize=1, thick=2
     oplot, wavelength_co2_ice, Qext_co2_ice(rei_co2_select,*), color=90, thick=2, linestyle=0
     xyouts, 0.65, 0.95, "rei = "+string(rei_co2_grid(rei_co2_select))+" microns", /normal, size=1
   endif

   if  (plot_ps eq 1) then begin 
     device, /close
   endif else begin
    print, "code stopped for plot viewing, press .c to continue"
    stop
   endelse

  ;-------  plot W (omega) ----------
  if (plot_ps eq 1) then begin
    print, "plotting to postscript: idl_cld_w.eps"
    set_plot,'PS'
    device,file='idl_cld_w.eps'
    device,/color,BITS=8, /ENCAPSULATED ;, /CMYK
    device,xsize=18.5,ysize=11,/CM
    device, set_font='Helvetica-Oblique', FONT_INDEX=20
    device, set_font='Helvetica-Bold', FONT_INDEX=19
    device, set_font='helvetica',FONT_INDEX=18
    loadct, colortable
  endif else begin
    print, "plotting to x-windows"
    set_plot,'x'
  endelse

   xr=[0.01,1000]
   yr=[-0.1,1.1]

   xplaceholder=findgen(10) & yplaceholder=findgen(10)
   plot, xplaceholder, yplaceholder,/nodata, $
         xrange=xr, xstyle=1, yrange=yr, ystyle=1, xthick=3.0, ythick=3.0, /xlog, $
         xtitle="wavelength(microns)", $
         ytitle="single scattering albedo" ;, $                    

   if (do_h2o_ice eq 1) then begin
     print, "rei_h2o: ", rei_h2o_grid(rei_h2o_select)
     oplot, wavelength_mid, W_h2o_ice_out(rei_h2o_select,*), color=240, psym=1, symsize=1, thick=2
     oplot, wavelength_h2o_ice, W_h2o_ice(rei_h2o_select,*), color=90, thick=2, linestyle=0
     xyouts, 0.65, 0.95, "rei = "+string(rei_h2o_grid(rei_h2o_select))+" microns", /normal, size=1
   endif

   if (do_h2o_liq eq 1) then begin
     print, "rel_h2o: ", rel_h2o_grid(rel_h2o_select)
     oplot, wavelength_mid, W_h2o_liq_out(rel_h2o_select,*), color=240, psym=1, symsize=1, thick=2
     oplot, wavelength_h2o_liq, W_h2o_liq(rel_h2o_select,*), color=90, thick=2, linestyle=0
     xyouts, 0.65, 0.95, "rel = "+string(rel_h2o_grid(rel_h2o_select))+" microns", /normal, size=1
   endif

   if (do_co2_ice eq 1) then begin
     print, "rei_co2: ", rei_co2_grid(rei_co2_select)
     oplot, wavelength_mid, W_co2_ice_out(rei_co2_select,*), color=240, psym=1, symsize=1, thick=2
     oplot, wavelength_co2_ice, W_co2_ice(rei_co2_select,*), color=90, thick=2, linestyle=0
     xyouts, 0.65, 0.95, "rei = "+string(rei_co2_grid(rei_co2_select))+" microns", /normal, size=1
   endif

   if  (plot_ps eq 1) then begin 
     device, /close
   endif else begin
    print, "code stopped for plot viewing, press .c to continue"
    stop
   endelse


  ;-------  plot G (asymmetry) ----------
  if (plot_ps eq 1) then begin
    print, "plotting to postscript: idl_cld_g.eps"
    set_plot,'PS'
    device,file='idl_cld_g.eps'
    device,/color,BITS=8, /ENCAPSULATED ;, /CMYK
    device,xsize=18.5,ysize=11,/CM
    device, set_font='Helvetica-Oblique', FONT_INDEX=20
    device, set_font='Helvetica-Bold', FONT_INDEX=19
    device, set_font='helvetica',FONT_INDEX=18
    loadct, colortable
  endif else begin
    print, "plotting to x-windows"
    set_plot,'x'
  endelse

   xr=[0.01,1000]
   yr=[-1.1,1.1]

   xplaceholder=findgen(10) & yplaceholder=findgen(10)
   plot, xplaceholder, yplaceholder,/nodata, $
         xrange=xr, xstyle=1, yrange=yr, ystyle=1, xthick=3.0, ythick=3.0, /xlog, $
         xtitle="wavelength (microns)", $
         ytitle="asymmetry parameter" ;, $                    

   if (do_h2o_ice eq 1) then begin
     print, "rei_h2o: ", rei_h2o_grid(rei_h2o_select)
     oplot, wavelength_mid, G_h2o_ice_out(rei_h2o_select,*), color=240, psym=1, symsize=1, thick=2
     oplot, wavelength_h2o_ice, G_h2o_ice(rei_h2o_select,*), color=90, thick=2, linestyle=0
     xyouts, 0.65, 0.95, "rei = "+string(rei_h2o_grid(rei_h2o_select))+" microns", /normal, size=1
   endif

   if (do_h2o_liq eq 1) then begin
     print, "rel_h2o: ", rel_h2o_grid(rel_h2o_select)
     oplot, wavelength_mid, G_h2o_liq_out(rel_h2o_select,*), color=240, psym=1, symsize=1, thick=2
     oplot, wavelength_h2o_liq, G_h2o_liq(rel_h2o_select,*), color=90, thick=2, linestyle=0
     xyouts, 0.65, 0.95, "rel = "+string(rel_h2o_grid(rel_h2o_select))+" microns", /normal, size=1
   endif

   if (do_co2_ice eq 1) then begin
     print, "rei_co2: ", rei_co2_grid(rei_co2_select)
     oplot, wavelength_mid, G_co2_ice_out(rei_co2_select,*), color=240, psym=1, symsize=1, thick=2
     oplot, wavelength_co2_ice, G_co2_ice(rei_co2_select,*), color=90, thick=2, linestyle=0
     xyouts, 0.65, 0.95, "rei = "+string(rei_co2_grid(rei_co2_select))+" microns", /normal, size=1
   endif

   if  (plot_ps eq 1) then begin 
     device, /close
   endif else begin
    print, "code stopped for plot viewing, press .c to continue"
    stop
   endelse


endif ;do_plot_qwg

endif ; do_mie

end
