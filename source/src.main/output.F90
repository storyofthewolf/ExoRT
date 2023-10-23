module output


use shr_kind_mod,       only: r8 => shr_kind_r8
use shr_const_mod
use ppgrid
use ioFileMod
use physconst
use radgrid

implicit none
public

contains

subroutine output_data ( sw_dTdt, lw_dTdt, &
                         lw_dnflux, lw_upflux, &
                         sw_dnflux, sw_upflux, &
                         lw_dnflux_spectral, lw_upflux_spectral, &
                         sw_dnflux_spectral, sw_upflux_spectral, &
                         sol_toa, &
                         pmid, pint, tmid, &
                         tint, zint, &
                         h2ommr, co2mmr, &
                         ch4mmr, c2h6mmr,  &
                         o2mmr,  o3mmr,  &
                         n2mmr, h2mmr, &
                         mw, cp)

implicit none
include 'netcdf.inc'

real(r8), intent(in), dimension(pver)  :: sw_dTdt
real(r8), intent(in), dimension(pver)  :: lw_dTdt
real(r8), intent(in), dimension(pverp) :: lw_dnflux
real(r8), intent(in), dimension(pverp) :: lw_upflux
real(r8), intent(in), dimension(pverp) :: sw_dnflux
real(r8), intent(in), dimension(pverp) :: sw_upflux
real(r8), intent(in), dimension(pverp,ntot_wavlnrng) :: lw_dnflux_spectral
real(r8), intent(in), dimension(pverp,ntot_wavlnrng) :: lw_upflux_spectral
real(r8), intent(in), dimension(pverp,ntot_wavlnrng) :: sw_dnflux_spectral
real(r8), intent(in), dimension(pverp,ntot_wavlnrng) :: sw_upflux_spectral
real(r8), intent(in)                   :: sol_toa
real(r8), intent(in), dimension(pver)  :: pmid
real(r8), intent(in), dimension(pverp) :: pint
real(r8), intent(in), dimension(pver)  :: tmid
real(r8), intent(in), dimension(pverp) :: tint
real(r8), intent(in), dimension(pverp) :: zint
real(r8), intent(in), dimension(pver)  :: h2ommr
real(r8), intent(in), dimension(pver)  :: co2mmr
real(r8), intent(in), dimension(pver)  :: ch4mmr
real(r8), intent(in), dimension(pver)  :: c2h6mmr
real(r8), intent(in), dimension(pver)  :: o2mmr
real(r8), intent(in), dimension(pver)  :: o3mmr
real(r8), intent(in), dimension(pver)  :: n2mmr
real(r8), intent(in), dimension(pver)  :: h2mmr
real(r8), intent(in)                   :: mw
real(r8), intent(in)                   :: cp

integer :: ncid, status, pverp_id, pver_id, one_id, nw_id
integer :: lwup_id, lwdn_id, swup_id, swdn_id, fsd_id
integer :: lwup_spectral_id, lwdn_spectral_id, swup_spectral_id, swdn_spectral_id
integer :: lwhr_id, swhr_id
integer :: ps_id, pmid_id, pint_id, tmid_id
integer :: tint_id, zint_id
integer :: h2o_id, co2_id, ch4_id, c2h6_id, o2_id, o3_id, n2_id, h2_id, mw_id, cp_id

! output concentrations
real(r8), dimension(pver)  :: h2o_out
real(r8), dimension(pver)  :: co2_out
real(r8), dimension(pver)  :: ch4_out
real(r8), dimension(pver)  :: c2h6_out
real(r8), dimension(pver)  :: o2_out
real(r8), dimension(pver)  :: o3_out
real(r8), dimension(pver)  :: h2_out
real(r8), dimension(pver)  :: n2_out


!--- specific humidity and dry mass mixing ratios

h2o_out(:)  = h2ommr(:) !*mwdry/mwh2o
co2_out(:)  = co2mmr(:) !*mwdry/mwco2
ch4_out(:)  = ch4mmr(:) !*mwdry/mwch4
c2h6_out(:) = c2h6mmr(:) !*mwdry/mwch4
o2_out(:)   = o2mmr(:)  !*mwdry/mwo2
o3_out(:)   = o3mmr(:)  !*mwdry/mwo3
h2_out(:)   = h2mmr(:)  !*mwdry/mwh2
n2_out(:)   = n2mmr(:)  !*mwdry/mwn2


write(*,*)  "making netcdf file RTprofile_out.nc"
status = NF_CREATE("RTprofile_out.nc",nf_clobber,ncid)
if (status /= nf_noerr) call handle_err(status)

! define dimension
status = NF_DEF_DIM(ncid,"pverp", pverp, pverp_id)
if (status /= nf_noerr) call handle_err(status)
status = NF_DEF_DIM(ncid,"pver", pver, pver_id)
if (status /= nf_noerr) call handle_err(status)
status = NF_DEF_DIM(ncid,"ONE", 1, one_id)
if (status /= nf_noerr) call handle_err(status)

status = NF_DEF_DIM(ncid,"ntot_wavlnrng", ntot_wavlnrng, nw_id)
if (status /= nf_noerr) call handle_err(status)


! defining variables
status = NF_DEF_VAR(ncid,"LWUP", nf_real, 1, pverp_id, lwup_id)
if (status /= nf_noerr) call handle_err(status)
status = NF_DEF_VAR(ncid,"LWDN", nf_real, 1, pverp_id, lwdn_id)
if (status /= nf_noerr) call handle_err(status)
status = NF_DEF_VAR(ncid,"SWUP", nf_real, 1, pverp_id, swup_id)
if (status /= nf_noerr) call handle_err(status)
status = NF_DEF_VAR(ncid,"SWDN", nf_real, 1, pverp_id, swdn_id)
if (status /= nf_noerr) call handle_err(status)
status = NF_DEF_VAR(ncid,"LWUP_SPECTRAL", nf_real, 2, [pverp_id,nw_id], lwup_spectral_id)
if (status /= nf_noerr) call handle_err(status)
status = NF_DEF_VAR(ncid,"LWDN_SPECTRAL", nf_real, 2, [pverp_id,nw_id], lwdn_spectral_id)
if (status /= nf_noerr) call handle_err(status)
status = NF_DEF_VAR(ncid,"SWUP_SPECTRAL", nf_real, 2, [pverp_id,nw_id], swup_spectral_id)
if (status /= nf_noerr) call handle_err(status)
status = NF_DEF_VAR(ncid,"SWDN_SPECTRAL", nf_real, 2, [pverp_id,nw_id], swdn_spectral_id)
if (status /= nf_noerr) call handle_err(status)
status = NF_DEF_VAR(ncid,"FSDTOA", nf_real, 1, one_id, fsd_id)
if (status /= nf_noerr) call handle_err(status)
status = NF_DEF_VAR(ncid,"LWHR", nf_real, 1, pver_id, lwhr_id)
if (status /= nf_noerr) call handle_err(status)
status = NF_DEF_VAR(ncid,"SWHR", nf_real, 1, pver_id, swhr_id)
if (status /= nf_noerr) call handle_err(status)
status = NF_DEF_VAR(ncid,"PMID", nf_real, 1, pver_id, pmid_id)
if (status /= nf_noerr) call handle_err(status)
status = NF_DEF_VAR(ncid,"PINT", nf_real, 1, pverp_id, pint_id)
if (status /= nf_noerr) call handle_err(status)
status = NF_DEF_VAR(ncid,"TMID", nf_real, 1, pver_id, tmid_id)
if (status /= nf_noerr) call handle_err(status)
status = NF_DEF_VAR(ncid,"TINT", nf_real, 1, pverp_id, tint_id)
if (status /= nf_noerr) call handle_err(status)
status = NF_DEF_VAR(ncid,"ZINT", nf_real, 1, pverp_id, zint_id)
if (status /= nf_noerr) call handle_err(status)
status = NF_DEF_VAR(ncid,"H2OMMR", nf_real, 1, pver_id, h2o_id)
if (status /= nf_noerr) call handle_err(status)
status = NF_DEF_VAR(ncid,"CO2MMR", nf_real, 1, pver_id, co2_id)
if (status /= nf_noerr) call handle_err(status)
status = NF_DEF_VAR(ncid,"CH4MMR", nf_real, 1, pver_id, ch4_id)
if (status /= nf_noerr) call handle_err(status)
status = NF_DEF_VAR(ncid,"C2H6MMR", nf_real, 1, pver_id, c2h6_id)
if (status /= nf_noerr) call handle_err(status)
status = NF_DEF_VAR(ncid,"O2MMR", nf_real, 1, pver_id, o2_id)
if (status /= nf_noerr) call handle_err(status)
status = NF_DEF_VAR(ncid,"O3MMR", nf_real, 1, pver_id, o3_id)
if (status /= nf_noerr) call handle_err(status)
status = NF_DEF_VAR(ncid,"N2MMR", nf_real, 1, pver_id, n2_id)
if (status /= nf_noerr) call handle_err(status)
status = NF_DEF_VAR(ncid,"H2MMR", nf_real, 1, pver_id, h2_id)
if (status /= nf_noerr) call handle_err(status)
status = NF_DEF_VAR(ncid,"MWDRY", nf_real, 1, one_id, mw_id)
if (status /= nf_noerr) call handle_err(status)
status = NF_DEF_VAR(ncid,"CPDRY", nf_real, 1, one_id, cp_id)
if (status /= nf_noerr) call handle_err(status)


! end definitions
status = NF_ENDDEF(ncid)
if (status /= nf_noerr) call handle_err(status)

! put variables
status = NF_PUT_VAR_double(ncid,lwup_id,lw_upflux)
if (status /= nf_noerr) call handle_err(status)
status = NF_PUT_VAR_double(ncid,lwdn_id,lw_dnflux)
if (status /= nf_noerr) call handle_err(status)
status = NF_PUT_VAR_double(ncid,swup_id,sw_upflux)
if (status /= nf_noerr) call handle_err(status)
status = NF_PUT_VAR_double(ncid,swdn_id,sw_dnflux)
if (status /= nf_noerr) call handle_err(status)
status = NF_PUT_VAR_double(ncid,lwup_spectral_id,lw_upflux_spectral)
if (status /= nf_noerr) call handle_err(status)
status = NF_PUT_VAR_double(ncid,lwdn_spectral_id,lw_dnflux_spectral)
if (status /= nf_noerr) call handle_err(status)
status = NF_PUT_VAR_double(ncid,swup_spectral_id,sw_upflux_spectral)
if (status /= nf_noerr) call handle_err(status)
status = NF_PUT_VAR_double(ncid,swdn_spectral_id,sw_dnflux_spectral)
if (status /= nf_noerr) call handle_err(status)
status = NF_PUT_VAR_double(ncid,fsd_id,sol_toa)
if (status /= nf_noerr) call handle_err(status)
status = NF_PUT_VAR_double(ncid,lwhr_id,lw_dTdt)
if (status /= nf_noerr) call handle_err(status)
status = NF_PUT_VAR_double(ncid,swhr_id,sw_dTdt)
if (status /= nf_noerr) call handle_err(status)
status = NF_PUT_VAR_double(ncid,pmid_id,pmid)
if (status /= nf_noerr) call handle_err(status)
status = NF_PUT_VAR_double(ncid,pint_id,pint)
if (status /= nf_noerr) call handle_err(status)
status = NF_PUT_VAR_double(ncid,tmid_id,tmid)
if (status /= nf_noerr) call handle_err(status)
status = NF_PUT_VAR_double(ncid,tint_id,tint)
if (status /= nf_noerr) call handle_err(status)
status = NF_PUT_VAR_double(ncid,zint_id,zint)
if (status /= nf_noerr) call handle_err(status)
status = NF_PUT_VAR_double(ncid,h2o_id,h2o_out)
if (status /= nf_noerr) call handle_err(status)
status = NF_PUT_VAR_double(ncid,co2_id,co2_out)
if (status /= nf_noerr) call handle_err(status)
status = NF_PUT_VAR_double(ncid,ch4_id,ch4_out)
if (status /= nf_noerr) call handle_err(status)
status = NF_PUT_VAR_double(ncid,c2h6_id,c2h6_out)
if (status /= nf_noerr) call handle_err(status)
status = NF_PUT_VAR_double(ncid,o2_id,o2_out)
if (status /= nf_noerr) call handle_err(status)
status = NF_PUT_VAR_double(ncid,o3_id,o3_out)
if (status /= nf_noerr) call handle_err(status)
status = NF_PUT_VAR_double(ncid,n2_id,n2_out)
if (status /= nf_noerr) call handle_err(status)
status = NF_PUT_VAR_double(ncid,h2_id,h2_out)
if (status /= nf_noerr) call handle_err(status)
status = NF_PUT_VAR_double(ncid,mw_id,mw)
if (status /= nf_noerr) call handle_err(status)
status = NF_PUT_VAR_double(ncid,cp_id,cp)
if (status /= nf_noerr) call handle_err(status)

! close netcdf
status = NF_CLOSE(ncid)
if (status /= nf_noerr) call handle_err(status)

end subroutine output_data

subroutine handle_err(cdfout)
  include 'netcdf.inc'
  integer, intent(in) :: cdfout
  if(cdfout /= nf_noerr) then
             print *, nf_strerror(cdfout)
             stop "Stopped"
    end if
end subroutine handle_err


end module
