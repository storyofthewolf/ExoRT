pro corrk2netcdf
;--------------------------------------------------
;AUTHOR: Wolf E.T.
;
;DESCRIPTION: 
; Convert HELIOS-K formatted correlated-k
; output data into netcdf files for ExoRT
; 
;--------------------------------------------------
;NOTES:
; Not happy with ExoRT's accuracy, range, gases?
; Make your own k-coefficients!
;  https://github.com/exoclime/HELIOS-K
;--------------------------------------------------

do_write_netcdf = 1

;\\ =====  define variable dimensions  ====
; Standard Definition resolutions for correlated-K sets
NTemp = 17
NPress = 61
NGauss = 8
Nbins = 68

; h2o_kv
;NTemp = 10
;NPress = 13
;NGauss = 8
;Nbins = 68
;// =====  define variable dimensions  ====  


atm2mb = 0.000986923267 ; conversion factor atm to mb

;name = 't10000'
;dir = '/gpfsm/dnb53/etwolf/models/HELIOS-K/run/'

name = "h2o"
dir = '/gpfsm/dnb53/etwolf/models/ExoRT/data/kdist/n68h2o/'

; define output variables
T  = dblarr(NTemp)
P  = dblarr(NPress)
;K  = dblarr(NTemp, NPress, NGauss)
K  = dblarr(NGauss, NPress, NTemp)
GW = dblarr(NGauss)

npoints = NTemp*NPress*NGauss
data_in = dblarr(5,npoints)
;row order: gauss weight, K, T, P, gauss#

for b=0, Nbins-1 do begin
  if (b lt 10) then  openfile = dir + "Out_" + name + "_bin000" + string(b) + ".dat"
  if (b ge 10) then  openfile = dir + "Out_" + name + "_bin00" + string(b) + ".dat"
  openfile = STRJOIN(STRSPLIT(openfile, /EXTRACT), '')

  ;open and read text files
  OPENR,lun,openfile, /Get_Lun
  READF, lun,data_in
  FREE_LUN, lun
  data_in(3,*) = data_in(3,*)/atm2mb
  ;sort input 
  for a=0,npoints-1 do begin
    if (a lt NGauss) then  GW(a) = data_in(0,a)
    if (a MOD NGauss eq 0 and a/NGauss lt NPress) then  P((a/NGauss) MOD (NGauss * NPRESS)) = data_in(3,a)
    if (a MOD NGauss*NPress eq 0) then T(a/(NGauss*NPress)) = data_in(2,a)
;    K((a/(NGauss*Npress)),(a/NGauss MOD NPress),(a MOD NGauss)) = data_in(1,a)
    K((a MOD NGauss),(a/NGauss MOD NPress),(a/(NGauss*Npress))) = data_in(1,a)
  endfor

  ;print, "-----"
  ;print, GW
  ;print, "-----"
  ;print, P
  ;print, "-----"
  ;print, T

  ; write netcdf
  name_root = "n68_8gpt_bin"
;  name_append = "_co2_hitran16_1e4Nnu_subL_c1000_gclima_qalpha1.0.nc"
  name_append = "_h2o_hitran16_Nnu1e4_noplinth_c25_gclima.nc"
;  name_append = "_ch4_hitran16_Nnu1e5_gclima.nc"

  if (do_write_netcdf eq 1) then begin
    if (b lt 9) then  writefile = name_root + '0' + string(b+1) + name_append
    if (b ge 9) then  writefile = name_root + string(b+1) + name_append
    writefile = STRJOIN(STRSPLIT(writefile, /EXTRACT), '')

    print, "writing to ",writefile
    id = NCDF_CREATE(writefile, /CLOBBER)
    dim1 = NCDF_DIMDEF(id,'NTemp',NTemp)
    dim2 = NCDF_DIMDEF(id,'NPress',NPress)
    dim3 = NCDF_DIMDEF(id,'NGauss',NGauss)

    varid1 = NCDF_VARDEF(id, 'Temperature',[dim1], /float)
    varid2 = NCDF_VARDEF(id, 'Pressure', [dim2], /float)
    varid3 = NCDF_VARDEF(id, 'GaussWeights', [dim3], /float)
    varid4 = NCDF_VARDEF(id, 'data', [dim3,dim2,dim1], /float)

    NCDF_ATTPUT, id, varid1, "title", "Temperature"
    NCDF_ATTPUT, id, varid1, "units", "K"
    NCDF_ATTPUT, id, varid2, "title", "Pressure"
    NCDF_ATTPUT, id, varid2, "units", "mb"
    NCDF_ATTPUT, id, varid2, "title", "GaussWeights"
    NCDF_ATTPUT, id, varid2, "units", "data"
    NCDF_ATTPUT, id, varid4, "title", "absorption coefficient"
    NCDF_ATTPUT, id, varid4, "units", "cm^2 molecule^-1"

    NCDF_CONTROL, id, /ENDEF

    NCDF_VARPUT, id, varid1, T
    NCDF_VARPUT, id, varid2, P
    NCDF_VARPUT, id, varid3, GW
    NCDF_VARPUT, id, varid4, K

    NCDF_CLOSE, id

  endif


endfor



end
