pro cloudoptics_mie
;--------------------------------------------------
;AUTHOR: Wolf E.T.
;
;DESCRIPTION:
; Produce optical property netcdf files for
; liquid cloud drops and ice cloud particles
;---------------------------------------------------
;NOTES: mie scattering for ice and water clouds

;---------------------------------------------------

do_liq = 1 ; flag, compute optical properties of liquid cloud droplets
do_ice = 1 ; flag, compute optical properties of ice particles
 
liq_filename = 'cloudoptics_liquid_mie_n73.nc'
ice_filename = 'cloudoptics_ice_mie_n73.nc'

;-------- input file list ----------
root = '/gpfsm/dnb53/etwolf/models'
LIQUID_MIE_FILE = root + '/ExoRT/data/cloud/refractive_indices/Segelstein_liquid.txt'
ICE_MIE_FILE    = root + '/ExoRT/data/cloud/refractive_indices/Warren_ice.txt'

;-------- read liquid files ----------

nwliq = 1261
header_liq = strarr(4)
data_liq = fltarr(3,nwliq)

wavelength_liq = fltarr(nwliq)
real_liq = fltarr(nwliq)
imaginary_liq = fltarr(nwliq)

;LIQUID_MIE_FILE

OPENR,lun,LIQUID_MIE_FILE,/Get_lun
READF,lun,header_liq
READF,lun,data_liq

wavelength_liq(*) = data_liq(0,*)
real_liq(*) = data_liq(1,*)
imaginary_liq(*) = data_liq(2,*)

FREE_LUN,lun

;sort out where size parameter too big, only relevant for far UV
;wl = where(2*!pi*30.0/wavelength_liq lt 12000.0, nwliq)
wl = where(2*!pi*100.0/wavelength_liq lt 12000.0, nwliq)
wavelength_liq = wavelength_liq(wl)
real_liq = real_liq(wl)
imaginary_liq = imaginary_liq(wl)

;-------- read ice files ----------

nwice = 486
data_ice = fltarr(3,nwice)

wavelength_ice = fltarr(nwice)
real_ice = fltarr(nwice)
imaginary_ice = fltarr(nwice)

;ICE_MIE_FILE

OPENR,lun,ICE_MIE_FILE,/Get_lun
READF,lun,data_ice

wavelength_ice(*) = data_ice(0,*)
real_ice(*) = data_ice(1,*)
imaginary_ice(*) = data_ice(2,*)

FREE_LUN,lun

;sort out where size parameter too big, only relevant for far UV
wi = where(2*!pi*300.0/wavelength_ice lt 12000.0, nwice)
wavelength_ice = wavelength_ice(wi)
real_ice = real_ice(wi)
imaginary_ice = imaginary_ice(wi)

;---------- define spectral intervals from model ----------------

;wavenumber interval edges [cm-1]
;28 bin
;wavenum_edge = [ 10., 350., 500., 630., 700., 820., 980., 1080., 1180., 1390., $
;                 1480., 1800., 2080., 2250., 2380., 2600., 3250., 4000., 4650., $
;                 5150., 6150., 7700., 8050., 12850., 16000., 22650., 29000., $
;                 38000., 50000. ]

;42 bin
;wavenum_edge  = [      10.0,       200.0,         350.0,        425.0, $
;                      500.0,       630.0,         700.0,        820.0, $
;                      980.0,       1100.0,        1180.0,       1390.0, $
;                     1480.0,       1800.0,        2080.0,       2200.0, $
;                     2380.0,       2600.0,        3250.0,       4000.0, $
;                     4650.0,       4900.0,        5150.0,       5650.0, $
;                     6150.0,       6650.0,        7000.0,       7700.0, $
;                     8050.0,       9100.0,        10000.0,      11000.0, $
;                    11800.0,       12850.0,       13450.0,      14450.0, $
;                    15150.0,       16000.0,       19300.0,      22650.0, $
;                    29000.0,       38000.0,       50000.0 ]

;63 bin LMD CO2 grid
;wavenum_edge =   [  10.,    100.,   160.,   220.,   280.,   330.,   380.,   440.,   495.,   545.,  617.,  667., 720., $
;                    800.,   875.,   940.,   1000.,  1065.,  1108.,  1200.,  1275.,  1350.,  1450., 1550., 1650.,  $
;                    1750.,  1850.,  1950.,  2050.,  2200.,  2397.,  2494.,  2796.,  3087.,  3425., 3760., 4030.,  $ 
;                    4540.,  4950.,  5370.,  5925.,  6390.,  6990.,  7650.,  8315.,  8850.,  9350., 9650., 10400., $
;                    11220., 11870., 12750., 13300., 14470., 15000., 15385., 16667., 18182., 20000., $
;                    22222., 25000., 28571., 33333., 40000. ]


;n68
;wavenum_edge =  [0.00E+00,   40.00000,   100.0000,   160.0000, $
;                220.0000,   280.0000,   330.0000,   380.0000, $
;                440.0000,   495.0000,   545.0000,   617.0000, $
;                667.0000,   720.0000,   800.0000,   875.0000, $
;                940.0000,   1000.000,   1065.000,   1108.000, $
;                1200.000,   1275.000,   1350.000,   1450.000, $
;                1550.000,   1650.000,   1750.000,   1850.000, $
;                1950.000,   2050.000,   2200.000,   2397.000, $
;                2494.000,   2796.000,   3087.000,   3425.000, $
;                3760.000,   4030.000,   4540.000,   4950.000, $
;                5370.000,   5925.000,   6390.000,   6990.000, $
;                7650.000,   8315.000,   8850.000,   9350.000, $
;                9650.000,   10400.00,   11220.00,   11870.00, $
;                12790.00,   13300.00,   14470.00,   15000.00, $
;                16000.00,   16528.00,   17649.00,   18198.00, $
;                18518.00,   22222.00,   25641.00,   29308.00, $
;                30376.00,   32562.00,   35087.00,   36363.00,  42087.00 ]

;n73 
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


nwgrid = n_elements(wavenum_edge)
wavelength_edge = fltarr(nwgrid)

;convert wavenum_edge to wavelength microns

wavelength_edge(*) = 1.0e4/wavenum_edge(*)


;effective radii for liquid drop cloud particles

rel_grid = [ 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, $
             11, 12, 13, 14, 15, 16, 17, 18, 19, 20, $
             21, 22, 23, 24, 25, 26, 27, 28, 29, 30 ]
rel_grid = findgen(30) + 1

;rel_grid = fltarr(1)
;rel_grid(0) = 30

nrel = n_elements(rel_grid)

;effective radiii for ice particles
;rei_grid= [ 5, 25, 45, 65., 85, 105, 125, 145, 185, 225, 265 ]

;rei_grid = [ 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, $
;             55, 60, 65, 70, 75, 80, 85, 90, 95, 100, $
;             105, 110, 115, 120, 125, 130, 135, 140, 145, 150, $
;             155, 160, 165, 170, 175, 180, 185, 190, 195, 200, $
;             205, 210, 215, 220, 225, 230, 235, 240, 245, 250, $
;             255, 260, 265, 270, 275, 280, 285, 290, 295, 300 ]
rei_grid = findgen(300) + 1

nrei = n_elements(rei_grid)

dge_grid = fltarr(nrei)   ; generalized effective size ice particles

;water optical constants on Segelstein grid
Qext_liq = fltarr(nrel, nwliq)
W_liq = fltarr(nrel, nwliq)
G_liq = fltarr(nrel, nwliq)

;ice optical constants on Warren grid
Qext_ice = fltarr(nrei, nwice)
W_ice = fltarr(nrei, nwice)
G_ice = fltarr(nrei, nwice)

;optical constants averaged over spectral intervals
Qext_liq_out = fltarr(nrel, nwgrid-1)
W_liq_out = fltarr(nrel, nwgrid-1)
G_liq_out = fltarr(nrel, nwgrid-1)

Qext_ice_out = fltarr(nrei, nwgrid-1)
W_ice_out = fltarr(nrei, nwgrid-1)
G_ice_out = fltarr(nrei, nwgrid-1)


if (do_liq eq 1) then begin   

  ; calculate optical properties of liquid droplets on Segelstein grid
  for i = 0, nrel-1 do begin
    for j = 0, nwliq-1 do begin
      cm = complex(real_liq(j), -imaginary_liq(j))
      alpha = 2.0 * !pi * rel_grid(i) / wavelength_liq(j)   
      mie_single, alpha, cm, dqext, dqscat, dqbk, dg
      Qext_liq(i,j) = dqext
      W_liq(i,j) = dqscat/dqext
      G_liq(i,j) = dg 
    endfor
  endfor

  ; average optical properties over all Segelstein points in spectral interval
  for i = 0, nrel-1 do begin
    for j = 0, nwgrid-2 do begin
      qx=where(wavelength_liq le wavelength_edge(j) and wavelength_liq ge wavelength_edge(j+1), qcount)
print, i,j,qcount
      Qext_liq_out(i,j) = MEAN(Qext_liq(i,qx))
      W_liq_out(i,j) = MEAN(W_liq(i,qx))
      G_liq_out(i,j) = MEAN(G_liq(i,qx))
    endfor
  endfor
stop
  ;print, "qext_liq",Qext_liq_out(0,0), Qext_liq_out(1,0), Qext_liq_out(2,0)
  ;print, "W_liq",W_liq_out(0,0), W_liq_out(1,0), W_liq_out(2,0)
  ;print, "G_liq",G_liq_out(0,0), G_liq_out(1,0), G_liq_out(2,0)

  ;write to netcdf file
  id = NCDF_CREATE(liq_filename, /CLOBBER)
  dim1 = NCDF_DIMDEF(id, 'rel_bins', nrel)
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

  NCDF_VARPUT, id, varid1, rel_grid
  NCDF_VARPUT, id, varid2, wavenum_edge
  NCDF_VARPUT, id, varid3, Qext_liq_out
  NCDF_VARPUT, id, varid4, W_liq_out
  NCDF_VARPUT, id, varid5, G_liq_out

  NCDF_CLOSE, id

  print, "wrote liquid cloud optical constants to, ",liq_filename

endif                   ; do_liq

if (do_ice eq 1) then begin   

   ; calculate optical properties for ice particles on Warren grid
   for i = 0, nrei-1 do begin
    for j = 0, nwice-1 do begin
      cm = complex(real_ice(j), -imaginary_ice(j))
      alpha = 2.0 * !pi * rei_grid(i) / wavelength_ice(j)   
      mie_single, alpha, cm, dqext, dqscat, dqbk, dg
      Qext_ice(i,j) = dqext
      W_ice(i,j) = dqscat/dqext
      G_ice(i,j) = dg 
    endfor       
  endfor

  ; average optical properties over all Warren points in spectral interval
  for i = 0, nrei-1 do begin
    for j = 0, nwgrid-2 do begin
      qx=where(wavelength_ice le wavelength_edge(j) and wavelength_ice ge wavelength_edge(j+1), qcount)
      Qext_ice_out(i,j) = MEAN(Qext_ice(i,qx))
      W_ice_out(i,j) = MEAN(W_ice(i,qx))
      G_ice_out(i,j) = MEAN(G_ice(i,qx))
    endfor
  endfor

  print, "qext_ice",Qext_ice_out(0,0), Qext_ice_out(1,0), Qext_ice_out(2,0)
  print, "W_ice",W_ice_out(0,0), W_ice_out(1,0), W_ice_out(2,0)
  print, "G_ice",G_ice_out(0,0), G_ice_out(1,0), G_ice_out(2,0)

  ;write to netcdf file
  id = NCDF_CREATE(ice_filename, /CLOBBER)
  dim1 = NCDF_DIMDEF(id, 'rei_bins', nrei)
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

  NCDF_VARPUT, id, varid1, rei_grid
  NCDF_VARPUT, id, varid2, wavenum_edge
  NCDF_VARPUT, id, varid3, Qext_ice_out
  NCDF_VARPUT, id, varid4, W_ice_out
  NCDF_VARPUT, id, varid5, G_ice_out

  NCDF_CLOSE, id

  print, "wrote ice cloud optical constants to, ",ice_filename

endif  ; do_ice


end
