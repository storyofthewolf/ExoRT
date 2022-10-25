pro makeCARMAOptics
;--------------------------------------------------
;AUTHOR: Wolf E.T.
;
;DESCRIPTION:
; Produce optical property netcdf files for CARMA cloud and aerosol grids
; Included are various spectral grids used by ExoCAM along with various plotting options.
; 
; The refractive indices are interpolated to the radiation grid
; BEFORE doing MIE or Fractal optics calculations
;---------------------------------------------------
;---------------------------------------------------

;;=============  options & files =================
;choose from pre-constructed CARMA elements
do_haze = 1 ; flag, compute optical properties of titan haze particles

plot_ps = 1              ; if eq 0, then plot to x windows 
do_plot_refract = 1      ; plot raw refractive indices

do_optical_calc = 1
do_plot_qwg = 1          ; plot extinction, single scattering, and asymmetry parameter  ;Currently not operable
do_write_netcdf = 1      ; flag to write netcdf outputs
 
carma_output_filename = 'haze_n68.nc'

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
root = '/gpfsm/dnb53/etwolf/models'
HAZE_REFRACT_FILE = root + '/ExoRT/data/aerosol/refractive_indices/Khare_haze.txt'

if (do_haze eq 1) then print, "calculating titan haze optical constants"

;;=============  read in refractive indice files =================
;-------- read haze files ----------
if (do_haze eq 1) then begin

  nw_haze = 90
  header_haze = strarr(4)
  data_haze = fltarr(4,nw_haze)

  wavelength_haze = fltarr(nw_haze)
  real_haze = fltarr(nw_haze)
  a = fltarr(nw_haze)
  b = fltarr(nw_haze)
  ;k = a*10^-6
  imaginary_haze = fltarr(nw_haze)

  OPENR,lun,HAZE_REFRACT_FILE,/Get_lun
  READF,lun,header_haze
  READF,lun,data_haze
  FREE_LUN,lun

  wavelength_haze(*) = data_haze(0,*)
  a(*) = data_haze(1,*)
  b(*) = data_haze(2,*)
  real_haze(*) = data_haze(3,*)
  imaginary_haze(*) = a(*)*10.0^(-b(*))

  wavelength_haze = reverse(wavelength_haze)
  real_haze = reverse(real_haze)
  imaginary_haze = reverse(imaginary_haze)


  ;sort out where size parameter too big, only relevant for far UV
  wl = where(2*!pi*100.0/wavelength_haze lt 12000.0, nw_haze)
  wavelength_haze = wavelength_haze(wl)
  real_haze = real_haze(wl)
  imaginary_haze = imaginary_haze(wl)

endif


;;----------------------------------------------------------------------
;;=============  Definitions of spectral interval grid =================
;;----------------------------------------------------------------------

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
wavelength_mid  = fltarr(nwgrid-1)
wavenum_mid     = fltarr(nwgrid-1)

;convert wavenum_edge to wavelength microns
wavelength_edge(*) = 1.0e4/wavenum_edge(*)

for nn=0, nwgrid-2 do begin
  wavenum_mid(nn) = (wavenum_edge(nn)+wavenum_edge(nn+1))/2.
endfor
wavelength_mid(*) = 1.0e4/wavenum_mid(*)

;print, "Spectral grid edges"
;for n=0, nwgrid-1 do print, n+1, wavenum_edge(n), wavelength_edge(n)
;print, " "
;print, "Spectral mid points"
;for n=0, nwgrid-2 do print, n+1, wavenum_mid(n), wavelength_mid(n)


;;-------------------------------------------------------------------
;;=============  Define the CARMA elements and bins =================
;;-------------------------------------------------------------------
;;!! WARNING  grid setup must match with ExoCAM-CARMA implementation
;;            inorder to have correct results

;; definte number of elements used
ncarma_elems = 1

;; quantities defined once for each element
is_fractal = intarr(ncarma_elems)
rho_particle = fltarr(ncarma_elems)
ncarma_bins = fltarr(ncarma_elems)
rmin = fltarr(ncarma_elems)
rmrat = fltarr(ncarma_elems)
rmassmin = fltarr(ncarma_elems)
vrfact = fltarr(ncarma_elems)
;fractal stuff
rmon  = fltarr(ncarma_elems)
alpha = fltarr(ncarma_elems)
xv    = fltarr(ncarma_elems)
angle = fltarr(ncarma_elems)

;;---  define elements ---''

;; Titan haze particles ;;
if (do_haze eq 1) then begin
  e1=0
  is_fractal(e1) = 1
  rho_particle(e1) = 0.64     ; [g cm-3] haze particle density
  ncarma_bins(e1) = 40        ; number of carma bins
  rmin(e1) = 1.0e-7           ; minimum particle size [cm] 
  rmrat(e1) = 2.0             ; ratio of bin masses
  rmassmin(e1) = 4.0/3.0*!pi*(rmin(e1)^3.0)*rho_particle(e1)    ;[g] mass of smallest bin
  vrfact(e1) = ((3.0/2.0/!pi/(rmrat(e1)+1.0))^(1.0/3.0))*(rmrat(e1)^(1.0/3.0)-1.0) ;; volume ratio factor

  ; define carma bin level properties
  dr = fltarr(ncarma_elems, ncarma_bins)    ;[um] vector containing the width of each bin in um 
  bin_radius = fltarr(ncarma_elems, ncarma_bins)  ;[cm] equivalent spherical mass
  bin_mass = fltarr(ncarma_elems, ncarma_bins)    ;[g] bin mass

  ;; quantities defined for elems x bins
  for ib=0, ncarma_bins(e1)-1 do begin
    bin_mass(e1,ib)=rmassmin(e1)*rmrat(e1)^(ib)
    bin_radius(e1,ib)=((3.0*bin_mass(ib))/(4.0*!pi*rho_particle(e1)))^(1.0/3)     ;[cm]
    dr(e1,ib)=vrfact(e1)*(bin_mass(e1,ib)/rho_particle(e1))^(1.0/3.0)*(1.0e4) ;[um] bin width, 1.0e4 converts cm to um 
  endfor

  if (is_fractal eq 1) then begin
    nmon  = fltarr(ncarma_elems, ncarma_bins)
    df    = fltarr(ncarma_elems, ncarma_bins)
    rf    = fltarr(ncarma_elems, ncarma_bins)
    alpha = fltarr(ncarma_elems, ncarma_bins)

    ;; fractal assumptions
    rmon(e1) = 50.0e-7 ; monomer size
    alpha(e1) = 1.0    ; fractal packing coefficient
    xv(e1) = 1.0
    angle(e1) = 0.0
    N = 1
    do_miess = "FALSE"  & do_miess_bool = 0 ; do coreshell?

    print, "ib+1, bin_radius(e1,ib), nmon(e1,ib), df(e1,ib), rf(e1,ib), rmon(e1)"
    for ib=0, ncarma_bins(e1)-1 do begin
      if (bin_radius(e1,ib) lt rmon) then begin  ; monomer particles
        nmon(e1,ib) = 1.0
        df(e1,ib) = 3.0
        rf(e1,ib) = bin_radius(e1,ib)
      endif else begin  ; fractal aggregates
        nmon(e1,ib) = (bin_radius(e1,ib)/rmon(e1))^3
        df(e1,ib) = 2.4 - 0.9*exp(-nmon(ib)/500.0)
        rf(e1,ib) = (1.0/alpha(e1))^(1.0/df(e1,ib))*bin_radius(e1,ib)^(3.0/df(e1,ib))*rmon(e1)^(1.0-3.0/df(e1,ib))
      endelse
      print, ib+1, bin_radius(e1,ib), nmon(e1,ib), df(e1,ib), rf(e1,ib), rmon(e1)
    endfor
  endif
endif
;; /Titan haze particles


;;----------------------------------------------------------------------
;;=============  Set Refractive Indices to ExoRT grid  =================
;;---------------------------------------------------------------------


ie=0
; Interpolate optical constants onto finer spectral grid.
; Sometimes no optical constant data exists in spectral bands 
; Then average optical constants into the spectral intervals
n_newgrid = (wavenum_edge(n_elements(wavenum_edge)-1)-wavenum_edge(0))+1
newgrid = findgen(n_newgrid)+wavenum_edge(0)+1
newgrid = 1.0e4/newgrid

;refractive indices on interp grid
rfi_real_intex = fltarr(ncarma_elems, n_elements(newgrid))
rfi_imag_intex = fltarr(ncarma_elems, n_elements(newgrid))

;refractive indices on ExoRT spectral intervals
rfi_real = fltarr(ncarma_elems, n_elements(newgrid))
rfi_imag = fltarr(ncarma_elems, n_elements(newgrid))

rfi_real_intex(ie,*) = interpol(real_haze,      wavelength_haze,  newgrid)
rfi_imag_intex(ie,*) = interpol(imaginary_haze, wavelength_haze,  newgrid)

qsum=0
for j = 0, nwgrid-2 do begin
  qx=where(newgrid le wavelength_edge(j) and newgrid ge wavelength_edge(j+1), qcount)
  qsum=qsum+qcount
  rfi_real(ie,j) = MEAN(rfi_real_intex(ie,qx)) 
  rfi_imag(ie,j) = MEAN(rfi_imag_intex(ie,qx)) 
endfor



;; -------------------------------------------------------------- ;;
;; =============  Optical Property Computations ================= ;;
;; -------------------------------------------------------------- ;;
if (do_optical_calc eq 1) then begin

  ;-------- Define Optical constant arrays ---------------

  ;optical constants averaged over ExoRT spectral intervals
  Kext_out = fltarr(ncarma_elems, ncarma_bins, nwgrid-1)
  Qext_out = fltarr(ncarma_elems, ncarma_bins, nwgrid-1)
  W_out    = fltarr(ncarma_elems, ncarma_bins, nwgrid-1)
  G_out    = fltarr(ncarma_elems, ncarma_bins, nwgrid-1)

  ;optical properties on input optical constant grid
  Kext = fltarr(ncarma_elems, ncarma_bins, nw_haze)
  Qext = fltarr(ncarma_elems, ncarma_bins, nw_haze)
  W    = fltarr(ncarma_elems, ncarma_bins, nw_haze)
  G    = fltarr(ncarma_elems, ncarma_bins, nw_haze)

  ; calculate optical properties
  ; loop over elements
  for ie=0, ncarma_elems-1 do begin

    ; MIE Calculations
    if (is_fractal(ie) eq 0) then begin   
      print, "element ",ie," ... doing mie calculation ..."
      for ib = 0, ncarma_bins(ie)-1 do begin
        print, ie,ib, " size=", bin_radius(ib)
        for j=0, nwgrid-2 do begin 
          cm = complex(rfi_real(ie,j), -rfi_imag(ie,j))
          alpha = 2.0 * !pi * bin_radius(ie,ib) / (wavelength_mid(j)/1.0e4)
          mie_single, alpha, cm, dqext, dqscat, dqbk, dg
          Qext(ie,ib,j) = dqext
          Kext(ie,ib,j) = 3.0/4.0*Qext(ie,ib,j)/rho_particle(ie)/bin_radius(ie,ib)
          W(ie,ib,j) = dqscat/dqext
          G(ie,ib,j) = dg 
          ;print, "q,k,w,g", Qext(i,j), Kext(i,j), W(i,j), G(i,j)
        endfor      
      endfor
    endif 

    ;FRACTAL Calculations
    if (is_fractal eq 1) then begin
      print, "element ",ie," ... doing fractal optics calculation ..."

      ; create_input file
      inputfile="INPUT"
;     fractal_model_path="/gpfsm/dnb53/etwolf/models/fractal_optics_coreshell/"
;     inputfile= fractal_model_path + inputfile
;     inptutfile=STRJOIN(STRSPLIT(inputfile,/EXTRACT,' '))

      for ib = 0, ncarma_bins(ie)-1 do begin
        print, ie,ib, " size=", bin_radius(ib)
        for j=0, nwgrid-2 do begin
          OPENW,lun,inputfile,/Get_LUN
          printf,lun, format = '(I5,A7)', N, do_miess
          if (do_miess_bool eq 1 ) then begin
            ; write ouput file with coreshell data, see file "INPUT_CORESHELL" 
            printf,lun, format = '("WVL(um)   K         N         nmon        a    df   rmon      xv   ang  rcore     KSHELL    NSHELL")'
            for i=0,N-1 do begin
              printf,lun,format = "(f10.4,  f10.5,  f10.5,  f12.1,   f5.2,     f5.2,  f10.5,   f5.1,  f5.1,     f10.5,    f10.5,    f10.5)", $
                                    wavelength_mid(j)*1.0e4, $
                                    rfi_imag(ie,j), $
                                    rfi_real(ie,j), $
                                    nmon(ie,ib), $
                                    alpha(ie), $
                                    df(ie,ib), $
                                    rmon(ie), $ 
                                    xv(ie), $
                                    angle(ie), $ 
                                    rcore(ie), $
                                    rfiSH(ie), $
                                    rfrSH(ie)
            endfor
          endif else begin
            ; write input file for homogeneous monomer data, see file "INPUT_MONOMER"
            printf,lun, format = '("WVL(um)   K         N         nmon        a     df    rmon     xv     ang")'  
            for i=0,N-1 do begin
              printf,lun,format = "(f10.4,    f10.5,    f10.5,    f12.1,      f5.2, f5.2, f10.5,   f5.1,  f5.1    )", $
                                    wavelength_mid(j), $
                                    rfi_imag(ie,j), $
                                    rfi_real(ie,j), $
                                    nmon(ie,ib), $
                                    alpha(ie), $
                                    df(ie,ib), $
                                    rmon(ie)*1.0e4, $ ; converts from cm to microns
                                    xv(ie), $
                                    angle(ie)
            endfor
          endelse
          FREE_LUN,lun
    
          ;---------- run fractal optics code -------------
          fractal_executable = "fractaloptics.exe"
          fractal_model_path="/gpfsm/dnb53/etwolf/models/fractal_optics_coreshell/"
          run_string = fractal_model_path + fractal_executable
          run_string = STRJOIN(STRSPLIT(run_string,/EXTRACT,' '))
          spawn, '/gpfsm/dnb53/etwolf/models/fractal_optics_coreshell/fractaloptics.exe'
          ;-----------------------------------------------
    
          ;---------- read fractal optics output -------------
          fractal_output = 'OUTPUT'
          header=strarr(1)
          OPENR,lun,fractal_output,/GET_LUN
          READF,lun,nrows
          data = fltarr(nrows, 5) 
          READF,lun,header
          READF,lun,data
          FREE_LUN, lun
          ;-----------------------------------------------

          Qext(ie,ib,j) = dqext
          Kext(ie,ib,j) = 3.0/4.0*Qext(ie,ib,j)/rho_particle(ie)/bin_radius(ie,ib)
          W(ie,ib,j) = dqscat/dqext
          G(ie,ib,j) = dg 
        endfor
      endfor
    endif
  endfor

  ; print outputs
  for ie=0,ncarma_elems-1 do begin
    for ib=0, ncarma_bins(ie)-1 do begin
      for j=0,nwgrid-2 do begin
        Qext_out(ie,ib,j) = Qext(ie,ib,j)
        Kext_out(ie,ib,j) = Kext(ie,ib,j)
        W_out(ie,ib,j)    = W(ie,ib,j)
        G_out(ie,ib,j)    = G(ie,ib,j)
      endfor
    endfor
  endfor

; print, "qext_haze",Qext_haze_out(0,0), Qext_haze_out(1,0), Qext_haze_out(2,0)
; print, "W_haze",W_haze_out(0,0), W_haze_out(1,0), W_haze_out(2,0)
; print, "G_haze",G_haze_out(0,0), G_haze_out(1,0), G_haze_out(2,0)

endif



if (do_write_netcdf) then begin
  ;write to netcdf file
  id = NCDF_CREATE(carma_output_filename, /CLOBBER)
  dim1 = NCDF_DIMDEF(id, 'nelements', ncarma_elems)
  dim2 = NCDF_DIMDEF(id, 'nbins', ncarma_bins)
  dim3 = NCDF_DIMDEF(id, 'nwavlrng', nwgrid-1)
  dim4 = NCDF_DIMDEF(id, 'nwave_edge',nwgrid) 

  varid1 = NCDF_VARDEF(id, 'rbins', [dim1,dim2], /float)
  varid2 = NCDF_VARDEF(id, 'wvnrng', [dim3], /float)
  varid3 = NCDF_VARDEF(id, 'Qext',[dim1,dim2,dim3], /float)
  varid4 = NCDF_VARDEF(id, 'Kext',[dim1,dim2,dim3], /float)
  varid5 = NCDF_VARDEF(id, 'W_liq', [dim1,dim2,dim3], /float)
  varid6 = NCDF_VARDEF(id, 'G_liq', [dim1,dim2,dim3], /float)

  NCDF_ATTPUT, id, varid1, "title", "carma bin equivalent sphere radii"
  NCDF_ATTPUT, id, varid1, "units", "microns"

  NCDF_ATTPUT, id, varid2, "title", "wavenumber at edges" 
  NCDF_ATTPUT, id, varid2, "units", "cm-1"

  NCDF_ATTPUT, id, varid3, "title", "extinction efficiency"
  NCDF_ATTPUT, id, varid3, "units", "unitless"

  NCDF_ATTPUT, id, varid4, "title", "mass extinction efficiency"
  NCDF_ATTPUT, id, varid4, "units", "cm2 g-1"

  NCDF_ATTPUT, id, varid5, "title", "single scattering albedo"
  NCDF_ATTPUT, id, varid5, "units", "unitless"

  NCDF_ATTPUT, id, varid6, "title", "asymmetry parameter"
  NCDF_ATTPUT, id, varid6, "units", "unitless"
  
  NCDF_CONTROL, id, /ENDEF

  NCDF_VARPUT, id, varid1, bin_radius
  NCDF_VARPUT, id, varid2, wavenum_edge
  NCDF_VARPUT, id, varid3, Qext_out
  NCDF_VARPUT, id, varid4, Kext_out
  NCDF_VARPUT, id, varid5, W_out
  NCDF_VARPUT, id, varid6, G_out
  NCDF_CLOSE, id

  print, "wrote carma optical properties to, ",carma_output_filename
endif ;do_write_netcdf






;----------------------------------------------------------------------
;============================ PLOTTING ================================
;----------------------------------------------------------------------

;;=============  plot refractive indices ================= 
if (do_plot_refract eq 1) then  begin

  ;--------------------------
  ; plot real refactive index
  ;--------------------------
  if (plot_ps eq 1) then begin
    print, "plotting to postscript: idl_carma_refrac_real.eps"
    set_plot,'PS'
    device,file='idl_carma_refrac_real.eps'
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

  if (do_haze eq 1) then begin
    oplot, wavelength_haze, real_haze, color=90, thick=2, linestyle=0
    oplot, wavelength_haze, real_haze, color=90, psym=1, symsize=1, thick=2
    oplot, wavelength_mid(*), rfi_real(0,*), color=240, psym=4, symsize=2, thick=2
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
    print, "plotting to postscript: idl_carma_refrac_imaginary.eps"
    set_plot,'PS'
    device,file='idl_carma_refrac_imaginary.eps'
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

  if (do_haze eq 1) then begin
    oplot, wavelength_haze, imaginary_haze, color=90, thick=2, linestyle=2
    oplot, wavelength_haze, imaginary_haze, color=90, psym=1, symsize=1, thick=2
    oplot, wavelength_mid(*), rfi_imag(0,*), color=240, psym=4, symsize=2, thick=2
  endif





  if  (plot_ps eq 1) then begin 
    device, /close
  endif else begin
    print, "code stopped for plot viewing, press .c to continue"
    stop
  endelse

endif


;;=============  plotting of Q, W, and G =================  

if (do_plot_qwg eq 1) then  begin

  ; select index of radii to plot data
  rcarma_select = 35 ;99

  ;-------  plot QEXT ----------
  if (plot_ps eq 1) then begin
    print, "plotting to postscript: idl_carma_q.eps"
    set_plot,'PS'
    device,file='idl_carma_q.eps'
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

   if (do_haze eq 1) then begin
     print, "rhaze: ", bin_radius(rcarma_select)
     ie=0
     oplot, wavelength_mid(*),  Qext_out(ie, rcarma_select, *), color=240, psym=1, symsize=1, thick=2
     oplot, wavelength_mid(*),  Qext_out(ie, rcarma_select, *), color=90, thick=2, linestyle=0
     xyouts, 0.65, 0.95, "rcarma = "+string(bin_radius(rcarma_select))+" microns", /normal, size=1
   endif

   if  (plot_ps eq 1) then begin 
     device, /close
   endif else begin
    print, "code stopped for plot viewing, press .c to continue"
    stop
   endelse

  ;-------  plot W (omega) ----------
  if (plot_ps eq 1) then begin
    print, "plotting to postscript: idl_carma_w.eps"
    set_plot,'PS'
    device,file='idl_carma_w.eps'
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

   if (do_haze eq 1) then begin
     print, "rhaze: ", bin_radius(rcarma_select)
     ie=0
     oplot, wavelength_mid(*), W_out(ie, rcarma_select, *), color=240, psym=1, symsize=1, thick=2
     oplot, wavelength_mid(*), W_out(ie, rcarma_select, *), color=90, thick=2, linestyle=0
     xyouts, 0.65, 0.95, "rcarma = "+string(bin_radius(rcarma_select))+" microns", /normal, size=1
   endif

   if  (plot_ps eq 1) then begin 
     device, /close
   endif else begin
    print, "code stopped for plot viewing, press .c to continue"
    stop
   endelse


  ;-------  plot G (asymmetry) ----------
  if (plot_ps eq 1) then begin
    print, "plotting to postscript: idl_carma_g.eps"
    set_plot,'PS'
    device,file='idl_carma_g.eps'
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

   if (do_haze eq 1) then begin
     print, "rhaze: ", bin_radius(rcarma_select)
     oplot, wavelength_mid, G_out(ie, rcarma_select, *), color=240, psym=1, symsize=1, thick=2
     oplot, wavelength_mid, G_out(ie, rcarma_select, *), color=90, thick=2, linestyle=0
     xyouts, 0.65, 0.95, "rhaze = "+string(bin_radius(rcarma_select))+" microns", /normal, size=1
   endif

   if  (plot_ps eq 1) then begin 
     device, /close
   endif else begin
    print, "code stopped for plot viewing, press .c to continue"
    stop
   endelse


endif ;do_plot_qwg



end
