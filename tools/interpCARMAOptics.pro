pro interpCARMAOptics
;-----------------------------
; Author: Wolf, E.T>
; Date: January 12, 2023
; Sometimes the fractal mean field code will crash, for extremely large
; particles and short wavelengths (<0.5 microns).  Generally, these particles
; exist only as a long tail of the particle size distribution.  This routine
; fills in the errant values via interpolation.
;

do_write_netcdf = 0
carma_output_filename = "haze_n68_b40_fractal_interp_fixed.nc"

plot_ps = 0
if (plot_ps eq 1) then begin
  print, "plotting to postscript"
endif else begin
  print, "plotting to x-windows"
endelse


filename = "/gpfsm/dnb53/etwolf/models/ExoRT/idl_progs/haze_n68_b40_fractal.nc"
;filename = "/gpfsm/dnb53/etwolf/models/ExoRT/idl_progs/haze_n68_b40_mie.nc"
ncid=ncdf_open(filename, /nowrite)
ncdf_varget,ncid,'rmrat' , rmrat
ncdf_varget,ncid,'rbins' , rbins
ncdf_varget,ncid,'wvnrng', wvnrng
ncdf_varget,ncid,'Qext', Qext
ncdf_varget,ncid,'Kext', Kext
ncdf_varget,ncid,'W', W
ncdf_varget,ncid,'G', G
ncdf_close,ncid

help, rmrat
help, rbins
help, wvnrng
help, Qext
help, Kext
help, W
help, G



; bad values determined by visual inspection
Qext(0,37, 58) = 0.0


nb = n_elements(Qext(0,*,0))
ns = n_elements(Qext(0,0,*))
; for output netcdf
ncarma_elems = n_elements(Qext(*,0,0))
ncarma_bins  = n_elements(rbins)
nwgrid = n_elements(wvnrng)
bin_radius = rbins
wavenum_edge = wvnrng

;rbins(*) = rbins(*) * 1.0e4

print, wvnrng
wvnmid=fltarr(ns)
wvlmid=fltarr(ns)
for i=0, ns-1 do begin
  wvnmid(i) = (wvnrng(i) + wvnrng(i+1))/2.
  wvlmid(i) = 1.0e4/wvnmid(i)
  print, i+1, wvlmid(i), wvnmid(i)
endfor


print, nb, ns


loadct,40
ytit= "Qext"
;xtit="!18Wavelength (!M"+string("155B)+"!3m)"
xtit="wvl"

xr=[0.2, 100.0] ; wvl

yr=[1.0e-6, 1e1] ; mie
yr=[1.0e-5, 1.0e5]

plot, wvlmid, Qext(0,0,*), xtitle=xtit, $                                                  
      xrange=xr, /xlog, xstyle=1, $
      /ylog, yrange=yr, ystyle=1, $
      ytitle=ytit, $
      charsize=0.7, xthick=3, ythick=3, /nodata


for nbi = 0, nb-1 do begin 
   x=where(Qext(0,nbi,*) eq 0.0, nx)
   xx=where(Qext(0,nbi,*) gt 0.0, nxx)

  print, nbi+1, rbins(nbi), ":", nx
  oplot, wvlmid, Qext(0,nbi,*), psym=4, symsize=2.0, thick=2.5
  oplot, wvlmid, Qext(0,nbi,*), linestyle=0, thick=2.5
  if (nx gt 0) then begin
    Qext_short = Qext(0,nbi,xx)
    Kext_short = Kext(0,nbi,xx)
    G_short    = G(0,nbi,xx)
    W_short    = W(0,nbi,xx)
    wvl_short  = wvlmid(xx)
    Qext_new   = interpol(Qext_short, wvl_short, wvlmid)
    Kext_new   = interpol(Kext_short, wvl_short, wvlmid)
    W_new      = interpol(W_short, wvl_short, wvlmid)
    G_new      = interpol(G_short, wvl_short, wvlmid)
;    oplot, wvlmid, Qext_new, psym=4, symsize=2.0, color=200, thick=3.0
;    oplot, wvlmid, Qext_new, linestyle=0, color=200, thick=3.0
    Qext(0,nbi,*) = Qext_new(*)
    Kext(0,nbi,*) = Kext_new(*)
    G(0,nbi,*) = G_new(*)
    W(0,nbi,*) = W_new(*)
  endif
;  oplot, wvlmid, Qext(0,nbi,*), linestyle=0, color=15

endfor


if (do_write_netcdf) then begin
  
  id = NCDF_CREATE(carma_output_filename, /CLOBBER)
  dim1 = NCDF_DIMDEF(id, 'nelements', ncarma_elems)
  dim2 = NCDF_DIMDEF(id, 'nbins', ncarma_bins)
  dim3 = NCDF_DIMDEF(id, 'nwavlrng', nwgrid-1)
  dim4 = NCDF_DIMDEF(id, 'nwave_edge',nwgrid)

  varid1 = NCDF_VARDEF(id, 'rmrat', [dim1], /float)
  varid2 = NCDF_VARDEF(id, 'rbins', [dim1,dim2], /float)
  varid3 = NCDF_VARDEF(id, 'wvnrng', [dim4], /float)
  varid4 = NCDF_VARDEF(id, 'Qext',[dim1,dim2,dim3], /float)
  varid5 = NCDF_VARDEF(id, 'Kext',[dim1,dim2,dim3], /float)
  varid6 = NCDF_VARDEF(id, 'W', [dim1,dim2,dim3], /float)
  varid7 = NCDF_VARDEF(id, 'G', [dim1,dim2,dim3], /float)

  NCDF_ATTPUT, id, varid1, "title", "carma ratio of bin masses"
  NCDF_ATTPUT, id, varid1, "units", "unitless"

  NCDF_ATTPUT, id, varid2, "title", "carma bin equivalent sphere radii"
  NCDF_ATTPUT, id, varid2, "units", "microns"

  NCDF_ATTPUT, id, varid3, "title", "wavenumber at edges"
  NCDF_ATTPUT, id, varid3, "units", "cm-1"

  NCDF_ATTPUT, id, varid4, "title", "extinction efficiency"
  NCDF_ATTPUT, id, varid4, "units", "unitless"

  NCDF_ATTPUT, id, varid5, "title", "mass extinction efficiency"
  NCDF_ATTPUT, id, varid5, "units", "m2 kg-1"

  NCDF_ATTPUT, id, varid6, "title", "single scattering albedo"
  NCDF_ATTPUT, id, varid6, "units", "unitless"

  NCDF_ATTPUT, id, varid7, "title", "asymmetry parameter"
  NCDF_ATTPUT, id, varid7, "units", "unitless"

  NCDF_CONTROL, id, /ENDEF
  NCDF_VARPUT, id, varid1, rmrat
  NCDF_VARPUT, id, varid2, bin_radius
  NCDF_VARPUT, id, varid3, wavenum_edge
  NCDF_VARPUT, id, varid4, Qext
  NCDF_VARPUT, id, varid5, Kext
  NCDF_VARPUT, id, varid6, W
  NCDF_VARPUT, id, varid7, G
  NCDF_CLOSE, id

  print, "wrote carma optical properties to, ",carma_output_filename

endif


end
