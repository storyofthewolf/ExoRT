pro plotprofile
;------------------------------------------------
;Author: Wolf, E.T. 
;Revision history 
;Created: August 2020
;
;DESCRIPTION: 
;Plot vertical profiles of relevant content contained
; in the ExoRT output file RTprofile_out.nc.
;------------------------------------------------

loadct, 33
nfiles = 1

fname = strarr(nfiles)
fname_short = strarr(nfiles)

fname_short(0) = " "
fname(0) =  "/discover/nobackup/etwolf/models/ExoRT/run/RTprofile_out.nc"
;fname(0) =  "/discover/nobackup/etwolf/models/ExoRT/run/RTprofile_out_TS340_H2Oonly_mtckd.nc"
;fname(0) =  "/discover/nobackup/etwolf/models/ExoRT/run/RTprofile_out_TS340_H2Oonly_nocont.nc"


line_index = [0,0,0,0,0,0,2,2,2,2,2,2]
ncolors=6
colorscale = floor(255/ncolors)
color_index =(floor(findgen(ncolors)*255/ncolors))

ncid=ncdf_open(fname(0), /nowrite)
ncdf_varget,ncid,'TMID',TMID_IN
ncdf_varget,ncid,'TINT',TINT_IN
ncdf_varget,ncid,'PMID',PMID_IN       ; wet pressures, includes water vapor
ncdf_varget,ncid,'PINT',PINT_IN
ncdf_varget,ncid,'ZINT',ZINT_IN
ncdf_varget,ncid,'H2OMMR',H2OMMR_IN   ; specific humidity kg(wv) / (kg(dryair)+kg(wv)) 
ncdf_varget,ncid,'CO2MMR',CO2MMR_IN   ; dry mass mixing ratio kg(co2) / kg(dry air)    
ncdf_varget,ncid,'CH4MMR',CH4MMR_IN   ; dry mass mixing ratio kg(ch4) / kg(dry air)    
ncdf_varget,ncid,'O2MMR',O2MMR_IN     ; dry mass mixing ratio kg(o2) / kg(dry air)     
ncdf_varget,ncid,'O3MMR',O3MMR_IN     ; dry mass mixing ratio kg(o3) / kg(dry air)     
ncdf_varget,ncid,'N2MMR',N2MMR_IN     ; dry mass mixing ratio kg(n2) / kg(dry air)     
ncdf_varget,ncid,'H2MMR',H2MMR_IN     ; dry mass mixing ratio kg(h2) / kg(dry air)     
ncdf_varget,ncid,'LWUP',LWUP_IN
ncdf_varget,ncid,'LWDN',LWDN_IN
ncdf_varget,ncid,'SWUP',SWUP_IN
ncdf_varget,ncid,'SWDN',SWDN_IN
ncdf_varget,ncid,'LWUP_SPECTRAL',LWUP_SPECTRAL_IN
ncdf_varget,ncid,'LWDN_SPECTRAL',LWDN_SPECTRAL_IN
ncdf_varget,ncid,'SWUP_SPECTRAL',SWUP_SPECTRAL_IN
ncdf_varget,ncid,'SWDN_SPECTRAL',SWDN_SPECTRAL_IN
ncdf_varget,ncid,'LWHR',LWHR_IN
ncdf_varget,ncid,'SWHR',SWHR_IN
ncdf_varget,ncid,'MWDRY',MWDRY_IN
ncdf_varget,ncid,'CPDRY',CPDRY_IN
ncdf_close,ncid

pverp = n_elements(LWUP_IN)

for z=0, pverp-1 do begin
  lwup_level_sum = total(LWUP_SPECTRAL_IN(z,*))
  lwdn_level_sum = total(LWDN_SPECTRAL_IN(z,*))
  swup_level_sum = total(SWUP_SPECTRAL_IN(z,*))
  swdn_level_sum = total(SWDN_SPECTRAL_IN(z,*))
endfor
print, "LWUP TOA: ", total(LWUP_SPECTRAL_IN(0,*)), LWUP_IN(0)
print, "LWDN TOA: ", total(LWDN_SPECTRAL_IN(0,*)), LWDN_IN(0)
print, "SWUP TOA: ", total(SWUP_SPECTRAL_IN(0,*)), SWUP_IN(0)
print, "SWDN TOA: ", total(SWDN_SPECTRAL_IN(0,*)), SWDN_IN(0)


;----------------    plotting   ----------------
; shortwave fluxes
;-----------------------------------------------
!P.font=0
set_plot,'x'
;  set_plot,'PS'
;  device,file='idl.ps'
;  device,/color,BITS=8 ;, /ENCAPSULATED ;, /CMYK       
;  device,xsize=18.5,ysize=15,/CM
;  device, set_font='Helvetica-Oblique', FONT_INDEX=20
;  device, set_font='Helvetica-Bold', FONT_INDEX=19
;  device, set_font='helvetica',FONT_INDEX=18

xr = [-5,400]
yr = [3.0e3,0.01]
;xr=[110,150]
;yr = [1000,300]
plot, SWDN_IN, PINT_IN/100.,/nodata, /ylog, $
         xrange=xr, xstyle=1, yrange=yr, ystyle=1, xthick=3.0, ythick=3.0, $
         ytitle="Pressure (Pa)", $
         xtitle="Flux (W m!U-2!N)", title="Shortwave fluxes"

oplot, SWDN_IN, PINT_IN/100., linestyle=2, thick=2
;oplot, SWDN_IN, PINT_IN/100., psym=2, thick=2
oplot, SWUP_IN, PINT_IN/100., linestyle=0, thick=2
;oplot, SWUP_IN, PINT_IN/100., psym=2, thick=2

stop

;----------------    plotting   ----------------
; longwave fluxes
;-----------------------------------------------
!P.font=0
set_plot,'x'
;  set_plot,'PS'
;  device,file='idl.ps'
;  device,/color,BITS=8 ;, /ENCAPSULATED ;, /CMYK       
;  device,xsize=18.5,ysize=15,/CM
;  device, set_font='Helvetica-Oblique', FONT_INDEX=20
;  device, set_font='Helvetica-Bold', FONT_INDEX=19
;  device, set_font='helvetica',FONT_INDEX=18

xr = [-5,250]
yr = [3.0e3,0.01]
;xr = [600,700]
;yr = [1000,300]
plot, LWDN_IN, PINT_IN/100.,/nodata, /ylog, $
         xrange=xr, xstyle=1, yrange=yr, ystyle=1, xthick=3.0, ythick=3.0, $
         ytitle="Pressure (mb)", $
         xtitle="Flux (W m!U-2!N)", title="Longwave fluxes"

oplot, LWDN_IN, PINT_IN/100, linestyle=2, thick=2
;oplot, LWDN_IN, PINT_IN/100, psym=2, thick=2
oplot, LWUP_IN, PINT_IN/100, linestyle=0, thick=2
;oplot, LWUP_IN, PINT_IN/100, psym=0, thick=2
stop

;----------------    plotting   ----------------
; heating/cooling rates
;-----------------------------------------------
!P.font=0
set_plot,'x'
;  set_plot,'PS'
;  device,file='idl.ps'
;  device,/color,BITS=8 ;, /ENCAPSULATED ;, /CMYK       
;  device,xsize=18.5,ysize=15,/CM
;  device, set_font='Helvetica-Oblique', FONT_INDEX=20
;  device, set_font='Helvetica-Bold', FONT_INDEX=19
;  device, set_font='helvetica',FONT_INDEX=18

xr = [-100,100]
yr = [3.0e3,0.01]
;xr = [-5,5]
;yr = [1000,300]
plot, LWDN_IN, PINT_IN/100.,/nodata, /ylog, $
         xrange=xr, xstyle=1, yrange=yr, ystyle=1, xthick=3.0, ythick=3.0, $
         ytitle="Pressure (mb)", $
         xtitle="Heating/Cooling Rate (K/day)", title="Heating/Cooling Rates)"

oplot, LWHR_IN, PMID_IN/100, linestyle=2, thick=2
;oplot, LWHR_IN, PMID_IN/100, psym=2, thick=2
oplot, SWHR_IN, PMID_IN/100, linestyle=0, thick=2
;oplot, LWHR_IN, PMID_IN/100, psym=2, thick=2
stop

;----------------    plotting   ----------------
; temperature
;-----------------------------------------------
!P.font=0
set_plot,'x'
;  set_plot,'PS'
;  device,file='idl.ps'
;  device,/color,BITS=8 ;, /ENCAPSULATED ;, /CMYK       
;  device,xsize=18.5,ysize=15,/CM
;  device, set_font='Helvetica-Oblique', FONT_INDEX=20
;  device, set_font='Helvetica-Bold', FONT_INDEX=19
;  device, set_font='helvetica',FONT_INDEX=18

xr = [150,400]
yr = [3.0e4,0.01]
;yr = [1000,300]


plot, TINT_IN, PINT_IN/100.,/nodata, /ylog, $
         xrange=xr, xstyle=1, yrange=yr, ystyle=1, xthick=3.0, ythick=3.0, $
         ytitle="Pressure (mb)", $
         xtitle="Flux (W m!U-2!N)",title="Temperature"

oplot, TINT_IN, PINT_IN/100, linestyle=0, thick=2
;oplot, TINT_IN, PINT_IN/100, psym=2, thick=2

  


;device,/close

;device,/close




end

