pro plotstellar
;--------------------------------------------------
;Author: Wolf, E.T.
; plot stellar spectrum based on ExoRT ready files
;--------------------------------------------------

nfiles = 2
nrtmax=100
filenames = strarr(nfiles)
;filenames(1) = "../data/solar/WD_5000K_n68.nc"
;filenames(0) = "../data/solar/WD_5000K_n84.nc"
filenames(0) = "../data/solar/toi700_n68.nc"
filenames(1) = "../data/solar/toi700_hst.nc"
;filenames(2) = "../data/solar/WD_3000K_n84.nc"
;filenames(3) = "../data/solar/bt-settl_50000_logg4.5_FeH0_n84.nc"


; select whether to plot in x windows or postscript
plot_ps = 1    ; if eq 0, then plot to x windows
               ; if eq 1, then plot to postscript
if (plot_ps eq 1) then begin
  print, "plotting to postscript"
endif else begin
  print, "plotting to x-windows"
endelse
wait, 5



color_index = [70,250,0,0]
line_index = [0,0,0,1]
lthk = [3,3,1,1]

solarflux_arr = fltarr(nfiles,nrtmax)
wavln_low_arr = fltarr(nfiles,nrtmax)
wavln_high_arr = fltarr(nfiles,nrtmax)
wavln_diff_arr = fltarr(nfiles,nrtmax)
nrt_arr = fltarr(nfiles)
for i=0, nfiles-1 do begin
  ncid=ncdf_open(filenames(i), /nowrite)
  ncdf_varget,ncid,'wav_low',wav_low
  ncdf_varget,ncid,'wav_high',wav_high
  ncdf_varget,ncid,'S0',S0
  ncdf_varget,ncid,'solarflux',solarflux
  nrt_arr(i) = n_elements(wav_low)
  solarflux_arr(i,0:nrt_arr(i)-1) = solarflux(0:nrt_arr(i)-1)
;print, "====="
;print, solarflux_arr(i,0:nrt_arr(i)-1)
;print, total(solarflux_arr(i,0:nrt_arr(i)-1))
  wavln_low_arr(i,0:nrt_arr(i)-1) = 1.0e4/wav_low(0:nrt_arr(i)-1)
  wavln_high_arr(i,0:nrt_arr(i)-1) = 1.0e4/wav_high(0:nrt_arr(i)-1)
  wavln_diff_arr(i,0:nrt_arr(i)-1) = 1.0e4/wav_low(0:nrt_arr(i)-1) - 1.0e4/wav_high(0:nrt_arr(i)-1)
;print, "-===="
;print, wavln_low_arr(i,0:nrt_arr(i)-1) 
;print, wavln_high_arr(i,0:nrt_arr(i)-1) 
;print, wavln_diff_arr(i,0:nrt_arr(i)-1)
  ncdf_close,ncid
endfor

print, nrt_arr

!P.Multi=[0,1,0]
loadct,40
!P.font=0

if (plot_ps eq 1) then begin
  set_plot,'PS'
  device,file='plotstellar.eps'
  device,/color,BITS=8, /ENCAPSULATED ;, /CMYK
  device,xsize=8.7,ysize=6,xoff=1.0,yoff=1.0,/CM
  device, set_font='Helvetica-Oblique', FONT_INDEX=20
  device, set_font='Helvetica-Bold', FONT_INDEX=19
  device, set_font='helvetica',FONT_INDEX=18
endif else begin
  set_plot, 'x'
endelse

;create artifical bar
xbar_arr = fltarr(nfiles,nrtmax*2)
ybar_arr = fltarr(nfiles,nrtmax*2)
nrt_maxplot = fltarr(nfiles)
for i=0, nfiles-1 do begin
  q=0
  for j=0,nrt_arr(i)-1 do begin
    xbar_arr(i,q) = wavln_low_arr(i,j)
    xbar_arr(i,q+1) = wavln_high_arr(i,j)
    ybar_arr(i,q) = solarflux_arr(i,j)/wavln_diff_arr(i,j)
    ybar_arr(i,q+1) = solarflux_arr(i,j)/wavln_diff_arr(i,j)
    q=q+2
  endfor
  nrt_maxplot(i) = q-1
;print, "q=",q
;print, "------"
;print, ybar_arr(i,*)
endfor



loadct,40
plot, xbar_arr(0,*), ybar_arr(0,*), xtitle="!18Wavelength (!M"+string("155B)+"!3m)", $
       xrange=[-0.1,4.0], xstyle=1, yrange=[0.0, 1200], ystyle=1, $
       ytitle="!18Radiance (W m!U-2!N !M"+string("155B)+"!3m)", $
       charsize=0.7, xthick=3, ythick=3, /nodata

for i=0, nfiles-1 do begin
  oplot, xbar_arr(i,0:nrt_maxplot(i)), ybar_arr(i,0:nrt_maxplot(i)), linestyle=line_index(i), thick=lthk(i), color=color_index(i)
endfor

;label outputs if desired
;xyouts, 0.333, 0.77, 'White Dwarf, 3000 K', color=100, charsize=0.8, /normal
;xyouts, 0.446, 0.405, "Ad Leo,3390 K ", color=250, charsize=0.8, /normal

;xyouts, 0.380, 0.45, 'White Dwarf, 5000 K', color=100, charsize=0.8, /normal
;xyouts, 0.196, 0.58, "Sun", color=0, charsize=0.8, /normal

xyouts, 0.61, 0.85, 'interpolated', color=color_index(0), charsize=0.6, /normal
xyouts, 0.61, 0.81, 'HST', color=color_index(1), charsize=0.6, /normal


if  (plot_ps eq 1) then begin
  device, /close
endif else begin
  stop
endelse


end
