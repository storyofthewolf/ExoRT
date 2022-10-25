module input

use shr_kind_mod,       only: r8 => shr_kind_r8
use shr_const_mod
use ppgrid
use ioFileMod
use exoplanet_mod,   only:  do_exo_clouds

implicit none
public  ! all public data used throughout code

integer, parameter :: ext_nazm_tshadow = 1
real(r8) :: ext_msdist_in
real(r8) :: ext_rtgt_in
real(r8) :: ext_solar_azm_ang_in
real(r8) :: ext_tazm_ang_in
real(r8) :: ext_tslope_ang_in
integer :: ext_tslas_tog_in
integer :: ext_tshadow_tog_in
real(r8), dimension(ext_nazm_tshadow) :: ext_cosz_horizon_in
real(r8), dimension(ext_nazm_tshadow) :: ext_TCx_obstruct_in
real(r8), dimension(ext_nazm_tshadow) :: ext_TCz_obstruct_in

! P,T,Z
real(r8) :: TS_in
real(r8) :: PS_in
real(r8), dimension(pverp) :: TINT_in
real(r8), dimension(pver) :: TMID_in
real(r8), dimension(pver) :: PMID_in
real(r8), dimension(pver) :: PDEL_in
real(r8), dimension(pver) :: PDELDRY_in
real(r8), dimension(pverp) :: PINT_in
real(r8), dimension(pverp) :: PINTDRY_in
real(r8), dimension(pverp) :: ZINT_in
! Absorbing gases
real(r8), dimension(pver) :: H2OMMR_in   ! specific humidity
real(r8), dimension(pver) :: CO2MMR_in   ! dry mass mixing ratio
real(r8), dimension(pver) :: CH4MMR_in   ! dry mass mixing ratio
real(r8), dimension(pver) :: O2MMR_in    ! dry mass mixing ratio
real(r8), dimension(pver) :: O3MMR_in    ! dry mass mixing ratio
real(r8), dimension(pver) :: N2MMR_in    ! dry mass mixing ratio
real(r8), dimension(pver) :: H2MMR_in    ! dry mass mixing ratio
! Surface Albedo and Emissivity
real(r8) :: ASDIR_in  
real(r8) :: ASDIF_in
real(r8) :: ALDIR_in
real(r8) :: ALDIF_in
real(r8) :: SRF_EMISS_in
! Cloud
real(r8), dimension(pver)  :: CICEWP_in,     CICEWP_zero
real(r8), dimension(pver)  :: CLIQWP_in,     CLIQWP_zero
real(r8), dimension(pver)  :: CFRC_in,       CFRC_zero
real(r8), dimension(pver)  :: CICEWP_CO2_in, CICEWP_CO2_zero
real(r8), dimension(pver)  :: REI_in,        REI_zero
real(r8), dimension(pver)  :: REL_in,        REL_zero
real(r8), dimension(pver)  :: REI_CO2_in,    REI_CO2_zero
! CARMA 
!real(r8), dimension(pver,nelem,nbin) :: CARMAMMR_in, CARMAMMR_zero
! Cosine of Zentih angle
real(r8) :: COSZRS_in
real(r8) :: MWDRY_in
real(r8) :: CPDRY_in


contains

subroutine initialize_to_zero
! P,T,Z
  TS_in = 0.
  PS_in = 0.
  TINT_in(:) = 0.
  TMID_in(:) = 0.
  PMID_in(:) = 0.
  PDEL_in(:) = 0.
  PDELDRY_in(:) = 0.
  PINTDRY_in(:) = 0.
  PINT_in(:) = 0.
  ZINT_in(:) = 0.
  ! Absorbing gases
  H2OMMR_in(:) = 0.   ! specific humidity
  CO2MMR_in(:) = 0.   ! dry mass mixing ratio
  CH4MMR_in(:) = 0.   ! dry mass mixing ratio
  O2MMR_in(:) = 0.    ! dry mass mixing ratio
  O3MMR_in(:) = 0.    ! dry mass mixing ratio
  N2MMR_in(:) = 0.    ! dry mass mixing ratio
  H2MMR_in(:) = 0.    ! dry mass mixing ratio
  ! Albedo
  ASDIR_in = 0.
  ASDIF_in = 0.
  ALDIR_in = 0.
  ALDIF_in = 0.
  SRF_EMISS_in = 1.   ! default to 1.
  ! Cloud
  CICEWP_in(:) = 0.       ;  CICEWP_zero(:) = 0.
  CLIQWP_in(:) = 0.       ;  CLIQWP_zero(:) = 0.
  CFRC_in(:) = 0.         ;  CFRC_zero(:) = 0.
  CICEWP_CO2_in(:) = 0.   ;  CICEWP_CO2_zero(:) = 0.
  REI_in(:) = 0.          ;  REI_zero(:) = 0.
  REL_in(:) = 0.          ;  REL_zero(:) = 0.
  REI_CO2_in(:) = 0.      ;  REI_CO2_zero(:) = 0.
  ! CARMA
 ! CARMAMMR_in(:) = 0.0    ;  CARMAMMR_zero(:) = 0.0
  ! Cosine of Zentih angle, mw, cp
  COSZRS_in = 0.
  MWDRY_in = 0.
  CPDRY_in = 0. 

end subroutine

subroutine input_profile

implicit none
include '../source/src.misc/netcdf.inc'

character(len=256) :: input_file
character(len=256) :: locfn
integer :: ncid
integer :: pverp_id, pver_id, npverp, npver 
integer :: ts_id, ps_id
integer :: tmid_id, tint_id, pmid_id, pdel_id, pint_id, zint_id
integer :: h2ommr_id, co2mmr_id, ch4mmr_id, o2mmr_id, o3mmr_id, h2mmr_id, n2mmr_id
integer :: cicewp_id, cliqwp_id, cicewp_co2_id, carmammr_id
integer :: rei_id, rel_id, rei_co2_id
integer :: asdir_id, asdif_id, aldir_id, aldif_id, srf_emiss_id
integer :: coszrs_id
integer :: mw_id, cp_id

write(6,*) "GET INPUT PROFILE DATA"
input_file = "RTprofile_in.nc"

call getfil(input_file, locfn, 0)
call wrap_open(locfn, 0, ncid)

call wrap_inq_dimid(ncid, 'pverp', pverp_id)
call wrap_inq_dimid(ncid, 'pver', pver_id)

call wrap_inq_dimlen(ncid, pverp_id, npverp)
call wrap_inq_dimlen(ncid, pver_id, npver)

write(*,*) "levels: ",pver
write(*,*) "interfaces: ",pverp

call wrap_inq_varid(ncid, 'ts', ts_id)
call wrap_inq_varid(ncid, 'ps', ps_id)
call wrap_inq_varid(ncid, 'tmid', tmid_id)
call wrap_inq_varid(ncid, 'tint', tint_id)
call wrap_inq_varid(ncid, 'pmid', pmid_id)
call wrap_inq_varid(ncid, 'pdel', pdel_id)
call wrap_inq_varid(ncid, 'pint', pint_id)
call wrap_inq_varid(ncid, 'zint', zint_id)
call wrap_inq_varid(ncid, 'h2ommr', h2ommr_id)
call wrap_inq_varid(ncid, 'co2mmr', co2mmr_id)
call wrap_inq_varid(ncid, 'ch4mmr', ch4mmr_id)
call wrap_inq_varid(ncid, 'o2mmr', o2mmr_id)
call wrap_inq_varid(ncid, 'o3mmr', o3mmr_id)
call wrap_inq_varid(ncid, 'n2mmr', n2mmr_id)
call wrap_inq_varid(ncid, 'h2mmr', h2mmr_id)
if (do_exo_clouds) then 
  call wrap_inq_varid(ncid, 'cicewp', cicewp_id)
  call wrap_inq_varid(ncid, 'cliqwp', cliqwp_id)
  call wrap_inq_varid(ncid, 'cicewp_co2', cicewp_co2_id)
  call wrap_inq_varid(ncid, 'rei', rei_id)
  call wrap_inq_varid(ncid, 'rel', rel_id)
  call wrap_inq_varid(ncid, 'rei_co2', rei_co2_id)
endif
!if (do_carma_exort) then
!  call wrap_inq_varid(ncid, 'carmammr', carmammr_id)
!endif
call wrap_inq_varid(ncid, 'asdir', asdir_id)
call wrap_inq_varid(ncid, 'asdif', asdif_id)
call wrap_inq_varid(ncid, 'aldir', aldir_id)
call wrap_inq_varid(ncid, 'aldif', aldif_id)
call wrap_inq_varid(ncid, 'srf_emiss', srf_emiss_id)
call wrap_inq_varid(ncid, 'coszrs', coszrs_id)
call wrap_inq_varid(ncid, 'mw', mw_id)
call wrap_inq_varid(ncid, 'cp', cp_id)

call wrap_get_var_realx(ncid, ts_id, TS_in)
call wrap_get_var_realx(ncid, ps_id, PS_in)
call wrap_get_var_realx(ncid, tmid_id, TMID_in)
call wrap_get_var_realx(ncid, tint_id, TINT_in)
call wrap_get_var_realx(ncid, pmid_id, PMID_in)
call wrap_get_var_realx(ncid, pdel_id, PDEL_in)
call wrap_get_var_realx(ncid, pint_id, PINT_in)
call wrap_get_var_realx(ncid, zint_id, ZINT_in)
call wrap_get_var_realx(ncid, h2ommr_id, H2OMMR_in)
call wrap_get_var_realx(ncid, co2mmr_id, CO2MMR_in)
call wrap_get_var_realx(ncid, ch4mmr_id, CH4MMR_in)
call wrap_get_var_realx(ncid, o2mmr_id, O2MMR_in)
call wrap_get_var_realx(ncid, o3mmr_id, O3MMR_in)
call wrap_get_var_realx(ncid, n2mmr_id, N2MMR_in)
call wrap_get_var_realx(ncid, h2mmr_id, H2MMR_in)
if (do_exo_clouds) then 
  call wrap_get_var_realx(ncid, cicewp_id, CICEWP_in)
  call wrap_get_var_realx(ncid, cliqwp_id, CLIQWP_in)
  call wrap_get_var_realx(ncid, cicewp_co2_id, CICEWP_CO2_in)
  call wrap_get_var_realx(ncid, rei_id, REI_in)
  call wrap_get_var_realx(ncid, rel_id, REL_in)
  call wrap_get_var_realx(ncid, rei_co2_id, REI_CO2_in)
endif
!if (do_carma_exort) then
!  call wrap_inq_varid(ncid, carmammr_id, CARMAMMR_in)
!endif
call wrap_get_var_realx(ncid, asdir_id, ASDIR_in)
call wrap_get_var_realx(ncid, asdif_id, ASDIF_in)
call wrap_get_var_realx(ncid, aldir_id, ALDIR_in)
call wrap_get_var_realx(ncid, aldif_id, ALDIF_in)
call wrap_get_var_realx(ncid, srf_emiss_id, SRF_EMISS_in)
call wrap_get_var_realx(ncid, coszrs_id, COSZRS_in)
call wrap_get_var_realx(ncid, mw_id, MWDRY_in)
call wrap_get_var_realx(ncid, cp_id, CPDRY_in)

end subroutine input_profile

end module input
