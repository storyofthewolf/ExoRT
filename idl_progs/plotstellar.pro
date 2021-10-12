pro plotstellar
;--------------------------------------------------
;Author: Wolf, E.T.
; plot stellar spectrum based on ExoRT ready files
;--------------------------------------------------

nfiles = 4
nrtmax=100
filenames = strarr(nfiles)
filenames(0) = "../data/solar/bt-settl_4500_logg4.5_FeH0_n84.nc"
filenames(1) = "../data/solar/bt-settl_7000_logg4.5_FeH0_n84.nc"
filenames(2) = "../data/solar/bt-settl_10000_logg4.5_FeH0_n84.nc"
filenames(3) = "../data/solar/bt-settl_50000_logg4.5_FeH0_n84.nc"


color_index = [250,0,90,120]
line_index = [0,0,0,0]
lthk = 2

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
set_plot,'PS'

device,file='plotstellar.eps'
device,/color,BITS=8, /ENCAPSULATED ;, /CMYK                                                                                           \
device,xsize=14.4,ysize=6,xoff=1.0,yoff=1.0,/CM

                                                                                                                                                   
device, set_font='Helvetica-Oblique', FONT_INDEX=20
device, set_font='Helvetica-Bold', FONT_INDEX=19
device, set_font='helvetica',FONT_INDEX=18

;create artifical bar
;\
                                                                                                                                                   

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
       xrange=[-0.5,3.5], xstyle=1, yrange=[0.0, 5000], ystyle=1, $
       ytitle="!18Radiance (W m!U-2!N !M"+string("155B)+"!3m)", $
charsize=0.7, xthick=3, ythick=3

for i=0, nfiles-1 do begin
  oplot, xbar_arr(i,0:nrt_maxplot(i)), ybar_arr(i,0:nrt_maxplot(i)), linestyle=line_index(i), thick=lthk, color=color_index(i)
endfor

print, "plotting, plotstellar.eps"
device, /close
set_plot,'X'

end
