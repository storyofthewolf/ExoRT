module io
!----------------------------------------------------------------------
! Consolidated 1-D input/output module for ExoRT (standalone driver).
! Combines the former input.F90 and output.F90:
!   read_namelist      - read user_nl_exort runtime config (1-D only)
!   initialize_to_zero - zero the input arrays
!   input_profile      - read RTprofile_in.nc
!   print_diagnostics  - echo primary fluxes to stdout (1-D only)
!   output_data        - write RTprofile_out.nc
!   handle_err         - NetCDF error handler
! Not used by the CESM/3-D path (which has its own I/O).
!----------------------------------------------------------------------

use shr_kind_mod,       only: r8 => shr_kind_r8
use shr_const_mod
use ppgrid
use ioFileMod
use physconst
use radgrid
use exoplanet_mod,      only: solar_file, exo_g, shr_const_scon, do_exo_clouds, do_exo_haze

implicit none
public  ! all public data used throughout code

! Namelist for 1-D runtime config; defaults are the module-variable
! initializers in exoplanet_mod. Read by read_namelist() below.
namelist /exort_config/ solar_file, shr_const_scon, exo_g, do_exo_clouds, do_exo_haze


integer, parameter :: ext_nazm_tshadow = 1
real(r8) :: ext_msdist_in
real(r8) :: ext_rtgt_in
real(r8) :: ext_solar_azm_ang_in
real(r8) :: ext_tazm_ang_in
real(r8) :: ext_tslope_ang_in
integer  :: ext_tslas_tog_in
integer  :: ext_tshadow_tog_in
real(r8), dimension(ext_nazm_tshadow) :: ext_cosz_horizon_in
real(r8), dimension(ext_nazm_tshadow) :: ext_TCx_obstruct_in
real(r8), dimension(ext_nazm_tshadow) :: ext_TCz_obstruct_in

! P,T,Z
real(r8) :: TS_in
real(r8) :: PS_in
real(r8), dimension(pverp) :: TINT_in
real(r8), dimension(pver)  :: TMID_in
real(r8), dimension(pver)  :: PMID_in
real(r8), dimension(pver)  :: PDEL_in
real(r8), dimension(pver)  :: PDELDRY_in
real(r8), dimension(pverp) :: PINT_in
real(r8), dimension(pverp) :: PINTDRY_in
real(r8), dimension(pverp) :: ZINT_in
! Absorbing gases
real(r8), dimension(pver) :: H2OMMR_in   ! specific humidity
real(r8), dimension(pver) :: CO2MMR_in   ! dry mass mixing ratio
real(r8), dimension(pver) :: CH4MMR_in   ! dry mass mixing ratio
real(r8), dimension(pver) :: C2H6MMR_in  ! dry mass mixing ratio
real(r8), dimension(pver) :: NH3MMR_in   ! dry mass mixing ratio
real(r8), dimension(pver) :: COMMR_in    ! dry mass mixing ratio
real(r8), dimension(pver) :: O2MMR_in    ! dry mass mixing ratio
real(r8), dimension(pver) :: O3MMR_in    ! dry mass mixing ratio
real(r8), dimension(pver) :: N2MMR_in    ! dry mass mixing ratio
real(r8), dimension(pver) :: H2MMR_in    ! dry mass mixing ratio
! Albedo
real(r8) :: ASDIR_in
real(r8) :: ASDIF_in
real(r8) :: ALDIR_in
real(r8) :: ALDIF_in
! Surface thermal emissivity (broadband scalar; default 1.0)
real(r8) :: SRF_EMISS_in
! Cloud
real(r8), dimension(pver)  :: CICEWP_in, CICEWP_zero
real(r8), dimension(pver)  :: CLIQWP_in, CLIQWP_zero
real(r8), dimension(pver)  :: CFRC_in, CFRC_zero
real(r8), dimension(pver)  :: REI_in, REI_zero
real(r8), dimension(pver)  :: REL_in, REL_zero
! CO2 ice cloud
real(r8), dimension(pver)  :: CICEWP_CO2_in, CICEWP_CO2_zero
real(r8), dimension(pver)  :: REI_CO2_in, REI_CO2_zero
! CARMA haze aerosols (binwise mass mixing ratio; zero if absent from input)
real(r8), dimension(pver,nelem_carma,nbin_carma) :: CARMAMMR_in
! Cosine of Zentih angle
real(r8) :: COSZRS_in
real(r8) :: MWDRY_in
real(r8) :: CPDRY_in



contains

subroutine read_namelist
!----------------------------------------------------------------------
! Read the 1-D runtime config namelist user_nl_exort if present, else
! keep the compile-time defaults from exoplanet_mod. Propagate the
! (possibly overridden) values into physconst (scon) and shr_const_mod
! (SHR_CONST_G), and echo the active config. 1-D only; the CESM path
! sets these via exoplanet_mod directly.
!----------------------------------------------------------------------
  integer :: nl_unit, nl_iostat
  logical :: nl_exists

  inquire(file='user_nl_exort', exist=nl_exists)
  if (nl_exists) then
    nl_unit = 10
    open(unit=nl_unit, file='user_nl_exort', status='old', action='read', iostat=nl_iostat)
    if (nl_iostat == 0) then
      read(nl_unit, nml=exort_config, iostat=nl_iostat)
      close(nl_unit)
      if (nl_iostat /= 0) then
        write(*,*) "WARNING: error reading user_nl_exort namelist (iostat=", nl_iostat, "); using defaults"
      endif
    endif
  endif
  ! Propagate possibly-overridden values to physconst and shr_const_mod
  scon        = shr_const_scon
  SHR_CONST_G = exo_g

  write(*,*) '=== exort_config ==========================='
  if (nl_exists) then
    write(*,*) 'namelist file: user_nl_exort (found)'
  else
    write(*,*) 'namelist file: user_nl_exort (not found, using defaults)'
  endif
  write(*,*) 'solar_file    = ', trim(solar_file)
  write(*,*) 'shr_const_scon= ', shr_const_scon, ' W m-2'
  write(*,*) 'exo_g         = ', exo_g, ' m s-2'
  write(*,*) 'do_exo_clouds = ', do_exo_clouds
  write(*,*) 'do_exo_haze   = ', do_exo_haze
  write(*,*) '============================================'

end subroutine read_namelist


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
  C2H6MMR_in(:) = 0.   ! dry mass mixing ratio
  NH3MMR_in(:) = 0.   ! dry mass mixing ratio
  COMMR_in(:) = 0.    ! dry mass mixing ratio
  O2MMR_in(:) = 0.    ! dry mass mixing ratio
  O3MMR_in(:) = 0.    ! dry mass mixing ratio
  N2MMR_in(:) = 0.    ! dry mass mixing ratio
  H2MMR_in(:) = 0.    ! dry mass mixing ratio
  ! Albedo
  ASDIR_in = 0.
  ASDIF_in = 0.
  ALDIR_in = 0.
  ALDIF_in = 0.
  ! Cloud
  CICEWP_in(:) = 0.
  CLIQWP_in(:) = 0.
  CFRC_in(:) = 0.
  REI_in(:) = 0.
  REL_in(:) = 0.
  CICEWP_CO2_in(:) = 0.
  REI_CO2_in(:) = 0.
  ! CARMA haze
  CARMAMMR_in(:,:,:) = 0.
  ! Cosine of Zentih angle, mw, cp
  COSZRS_in = 0.
  MWDRY_in = 0.
  CPDRY_in = 0. 

end subroutine

subroutine input_profile

  implicit none
    include 'netcdf.inc'

  character(len=256) :: input_file
  character(len=256) :: locfn
  integer :: ncid
  integer :: pverp_id, pver_id, npverp, npver
  integer :: ts_id, ps_id
  integer :: tmid_id, tint_id, pmid_id, pdel_id, pint_id, zint_id
  integer :: h2ommr_id, co2mmr_id, ch4mmr_id, c2h6mmr_id, nh3mmr_id, commr_id
  integer :: o2mmr_id, o3mmr_id, h2mmr_id, n2mmr_id
  integer :: cicewp_id, cliqwp_id
  integer :: rei_id, rel_id
  integer :: cicewp_co2_id, rei_co2_id
  integer :: carmammr_id
  integer :: asdir_id, asdif_id, aldir_id, aldif_id
  integer :: srf_emiss_id
  integer :: coszrs_id
  integer :: mw_id, cp_id
  integer :: ret

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

  ! Required: P, T, Z
  call wrap_inq_varid(ncid, 'ts', ts_id)
  call wrap_inq_varid(ncid, 'ps', ps_id)
  call wrap_inq_varid(ncid, 'tmid', tmid_id)
  call wrap_inq_varid(ncid, 'tint', tint_id)
  call wrap_inq_varid(ncid, 'pmid', pmid_id)
  call wrap_inq_varid(ncid, 'pdel', pdel_id)
  call wrap_inq_varid(ncid, 'pint', pint_id)
  call wrap_inq_varid(ncid, 'zint', zint_id)

  call wrap_get_var_realx(ncid, ts_id, TS_in)
  call wrap_get_var_realx(ncid, ps_id, PS_in)
  call wrap_get_var_realx(ncid, tmid_id, TMID_in)
  call wrap_get_var_realx(ncid, tint_id, TINT_in)
  call wrap_get_var_realx(ncid, pmid_id, PMID_in)
  call wrap_get_var_realx(ncid, pdel_id, PDEL_in)
  call wrap_get_var_realx(ncid, pint_id, PINT_in)
  call wrap_get_var_realx(ncid, zint_id, ZINT_in)

  ! Required: albedos, cosz, mw, cp
  call wrap_inq_varid(ncid, 'asdir', asdir_id)
  call wrap_inq_varid(ncid, 'asdif', asdif_id)
  call wrap_inq_varid(ncid, 'aldir', aldir_id)
  call wrap_inq_varid(ncid, 'aldif', aldif_id)
  call wrap_inq_varid(ncid, 'coszrs', coszrs_id)
  call wrap_inq_varid(ncid, 'mw', mw_id)
  call wrap_inq_varid(ncid, 'cp', cp_id)

  call wrap_get_var_realx(ncid, asdir_id, ASDIR_in)
  call wrap_get_var_realx(ncid, asdif_id, ASDIF_in)
  call wrap_get_var_realx(ncid, aldir_id, ALDIR_in)
  call wrap_get_var_realx(ncid, aldif_id, ALDIF_in)
  call wrap_get_var_realx(ncid, coszrs_id, COSZRS_in)
  call wrap_get_var_realx(ncid, mw_id, MWDRY_in)
  call wrap_get_var_realx(ncid, cp_id, CPDRY_in)

  ! Optional: surface thermal emissivity (default 1.0 if absent)
  SRF_EMISS_in = 1.0
  ret = nf_inq_varid(ncid, 'srf_emiss', srf_emiss_id)
  if (ret == NF_NOERR) then
    write(6,*) "  srf_emiss:  found"
    call wrap_get_var_realx(ncid, srf_emiss_id, SRF_EMISS_in)
  else
    write(6,*) "  srf_emiss:  not found, set to 1.0"
  end if

  ! Optional gases: zero-initialized above; read only if present
  write(6,*) "--- gas species in input file ---"
  ret = nf_inq_varid(ncid, 'h2ommr', h2ommr_id)
  if (ret == NF_NOERR) then
    write(6,*) "  h2ommr:  found"
    call wrap_get_var_realx(ncid, h2ommr_id, H2OMMR_in)
  else
    write(6,*) "  h2ommr:  not found, set to zero"
  end if

  ret = nf_inq_varid(ncid, 'co2mmr', co2mmr_id)
  if (ret == NF_NOERR) then
    write(6,*) "  co2mmr:  found"
    call wrap_get_var_realx(ncid, co2mmr_id, CO2MMR_in)
  else
    write(6,*) "  co2mmr:  not found, set to zero"
  end if

  ret = nf_inq_varid(ncid, 'ch4mmr', ch4mmr_id)
  if (ret == NF_NOERR) then
    write(6,*) "  ch4mmr:  found"
    call wrap_get_var_realx(ncid, ch4mmr_id, CH4MMR_in)
  else
    write(6,*) "  ch4mmr:  not found, set to zero"
  end if

  ret = nf_inq_varid(ncid, 'c2h6mmr', c2h6mmr_id)
  if (ret == NF_NOERR) then
    write(6,*) "  c2h6mmr: found"
    call wrap_get_var_realx(ncid, c2h6mmr_id, C2H6MMR_in)
  else
    write(6,*) "  c2h6mmr: not found, set to zero"
  end if

  ret = nf_inq_varid(ncid, 'nh3mmr', nh3mmr_id)
  if (ret == NF_NOERR) then
    write(6,*) "  nh3mmr:  found"
    call wrap_get_var_realx(ncid, nh3mmr_id, NH3MMR_in)
  else
    write(6,*) "  nh3mmr:  not found, set to zero"
  end if

  ret = nf_inq_varid(ncid, 'commr', commr_id)
  if (ret == NF_NOERR) then
    write(6,*) "  commr:   found"
    call wrap_get_var_realx(ncid, commr_id, COMMR_in)
  else
    write(6,*) "  commr:   not found, set to zero"
  end if

  ret = nf_inq_varid(ncid, 'o2mmr', o2mmr_id)
  if (ret == NF_NOERR) then
    write(6,*) "  o2mmr:   found"
    call wrap_get_var_realx(ncid, o2mmr_id, O2MMR_in)
  else
    write(6,*) "  o2mmr:   not found, set to zero"
  end if

  ret = nf_inq_varid(ncid, 'o3mmr', o3mmr_id)
  if (ret == NF_NOERR) then
    write(6,*) "  o3mmr:   found"
    call wrap_get_var_realx(ncid, o3mmr_id, O3MMR_in)
  else
    write(6,*) "  o3mmr:   not found, set to zero"
  end if

  ret = nf_inq_varid(ncid, 'n2mmr', n2mmr_id)
  if (ret == NF_NOERR) then
    write(6,*) "  n2mmr:   found"
    call wrap_get_var_realx(ncid, n2mmr_id, N2MMR_in)
  else
    write(6,*) "  n2mmr:   not found, set to zero"
  end if

  ret = nf_inq_varid(ncid, 'h2mmr', h2mmr_id)
  if (ret == NF_NOERR) then
    write(6,*) "  h2mmr:   found"
    call wrap_get_var_realx(ncid, h2mmr_id, H2MMR_in)
  else
    write(6,*) "  h2mmr:   not found, set to zero"
  end if

  ! Optional: cloud properties (zero-initialized above; no clouds if absent)
  write(6,*) "--- cloud variables in input file ---"
  ret = nf_inq_varid(ncid, 'cicewp', cicewp_id)
  if (ret == NF_NOERR) then
    write(6,*) "  cicewp:  found"
    call wrap_get_var_realx(ncid, cicewp_id, CICEWP_in)
  else
    write(6,*) "  cicewp:  not found, set to zero"
  end if

  ret = nf_inq_varid(ncid, 'cliqwp', cliqwp_id)
  if (ret == NF_NOERR) then
    write(6,*) "  cliqwp:  found"
    call wrap_get_var_realx(ncid, cliqwp_id, CLIQWP_in)
  else
    write(6,*) "  cliqwp:  not found, set to zero"
  end if

  ret = nf_inq_varid(ncid, 'rei', rei_id)
  if (ret == NF_NOERR) then
    write(6,*) "  rei:     found"
    call wrap_get_var_realx(ncid, rei_id, REI_in)
  else
    write(6,*) "  rei:     not found, set to zero"
  end if

  ret = nf_inq_varid(ncid, 'rel', rel_id)
  if (ret == NF_NOERR) then
    write(6,*) "  rel:     found"
    call wrap_get_var_realx(ncid, rel_id, REL_in)
  else
    write(6,*) "  rel:     not found, set to zero"
  end if

  ret = nf_inq_varid(ncid, 'cicewp_co2', cicewp_co2_id)
  if (ret == NF_NOERR) then
    write(6,*) "  cicewp_co2:  found"
    call wrap_get_var_realx(ncid, cicewp_co2_id, CICEWP_CO2_in)
  else
    write(6,*) "  cicewp_co2:  not found, set to zero"
  end if

  ret = nf_inq_varid(ncid, 'rei_co2', rei_co2_id)
  if (ret == NF_NOERR) then
    write(6,*) "  rei_co2:     found"
    call wrap_get_var_realx(ncid, rei_co2_id, REI_CO2_in)
  else
    write(6,*) "  rei_co2:     not found, set to zero"
  end if

  ! Optional: CARMA haze binwise mass mixing ratio (zero-initialized above)
  write(6,*) "--- aerosol variables in input file ---"
  ret = nf_inq_varid(ncid, 'carmammr', carmammr_id)
  if (ret == NF_NOERR) then
    write(6,*) "  carmammr:  found"
    call wrap_get_var_realx(ncid, carmammr_id, CARMAMMR_in)
  else
    write(6,*) "  carmammr:  not found, set to zero"
  end if

  write(6,*) "---------------------------------"

end subroutine input_profile


subroutine print_diagnostics ( sol_toa, vis_dir, vis_dif, nir_dir, nir_dif, &
                               sw_dnflux, sw_upflux, lw_dnflux, lw_upflux )
!----------------------------------------------------------------------
! Echo the primary 1-D diagnostic fluxes to stdout (TOA/surface, direct/
! diffuse). 1-D only; not part of the CESM output path.
!----------------------------------------------------------------------
  real(r8), intent(in) :: sol_toa
  real(r8), intent(in) :: vis_dir, vis_dif, nir_dir, nir_dif
  real(r8), intent(in), dimension(pverp) :: sw_dnflux, sw_upflux
  real(r8), intent(in), dimension(pverp) :: lw_dnflux, lw_upflux

    write(*,*) "------------------------------------"
    write(*,*) "Top-Model Downwelling Stellar"
    write(*,*) 'sol_toa', sol_toa
    write(*,*) "Surface downwelling fluxes"
    write(*,*) "vis_dir", vis_dir
    write(*,*) "vis_dif", vis_dif
    write(*,*) "nir_dir", nir_dir
    write(*,*) "nir_dif", nir_dif
    write(*,*) "total direct", vis_dir+nir_dir
    write(*,*) "total diffuse", vis_dif+nir_dif
    write(*,*) "TOA and Surface Fluxes"
    write(*,*) "SW DN TOA/SURF", sw_dnflux(1), sw_dnflux(pverp)
    write(*,*) "SW UP TOA/SURF", sw_upflux(1), sw_upflux(pverp)
    write(*,*) "LW DN TOA/SURF", lw_dnflux(1), lw_dnflux(pverp)
    write(*,*) "LW UP TOA/SURF", lw_upflux(1), lw_upflux(pverp)

end subroutine print_diagnostics


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
                         nh3mmr, commr, &
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
real(r8), intent(in), dimension(pver)  :: nh3mmr
real(r8), intent(in), dimension(pver)  :: commr
real(r8), intent(in), dimension(pver)  :: h2mmr
real(r8), intent(in)                   :: mw
real(r8), intent(in)                   :: cp

integer :: ncid, status, pverp_id, pver_id, one_id, nw_id
integer :: lwup_id, lwdn_id, swup_id, swdn_id, fsd_id
integer :: lwup_spectral_id, lwdn_spectral_id, swup_spectral_id, swdn_spectral_id
integer :: lwhr_id, swhr_id
integer :: ps_id, pmid_id, pint_id, tmid_id
integer :: tint_id, zint_id
integer :: h2o_id, co2_id, ch4_id, c2h6_id, nh3_id, co_id, o2_id, o3_id, n2_id, h2_id, mw_id, cp_id

! output concentrations
real(r8), dimension(pver)  :: h2o_out
real(r8), dimension(pver)  :: co2_out
real(r8), dimension(pver)  :: ch4_out
real(r8), dimension(pver)  :: c2h6_out
real(r8), dimension(pver)  :: nh3_out
real(r8), dimension(pver)  :: co_out
real(r8), dimension(pver)  :: o2_out
real(r8), dimension(pver)  :: o3_out
real(r8), dimension(pver)  :: h2_out
real(r8), dimension(pver)  :: n2_out


!--- specific humidity and dry mass mixing ratios

h2o_out(:)  = h2ommr(:) !*mwdry/mwh2o
co2_out(:)  = co2mmr(:) !*mwdry/mwco2
ch4_out(:)  = ch4mmr(:) !*mwdry/mwch4
c2h6_out(:) = c2h6mmr(:) !*mwdry/mwch4
nh3_out(:)  = nh3mmr(:)
co_out(:)   = commr(:)
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
status = NF_PUT_ATT_TEXT (ncid,lwup_id,'long_name', 24, 'longwave upwelling flux')
status = NF_PUT_ATT_TEXT (ncid,lwup_id,'units', 4, 'W/m2')
if (status /= nf_noerr) call handle_err(status)

status = NF_DEF_VAR(ncid,"LWDN", nf_real, 1, pverp_id, lwdn_id)
status = NF_PUT_ATT_TEXT (ncid,lwdn_id,'long_name', 26, 'longwave downwelling flux')
status = NF_PUT_ATT_TEXT (ncid,lwdn_id,'units', 4, 'W/m2')
if (status /= nf_noerr) call handle_err(status)

status = NF_DEF_VAR(ncid,"SWUP", nf_real, 1, pverp_id, swup_id)
status = NF_PUT_ATT_TEXT (ncid,swup_id,'long_name', 25, 'shortwave upwelling flux')
status = NF_PUT_ATT_TEXT (ncid,swup_id,'units', 4, 'W/m2')
if (status /= nf_noerr) call handle_err(status)

status = NF_DEF_VAR(ncid,"SWDN", nf_real, 1, pverp_id, swdn_id)
status = NF_PUT_ATT_TEXT (ncid,swdn_id,'long_name', 27, 'shortwave downwelling flux')
status = NF_PUT_ATT_TEXT (ncid,swdn_id,'units', 4, 'W/m2')
if (status /= nf_noerr) call handle_err(status)

status = NF_DEF_VAR(ncid,"LWUP_SPECTRAL", nf_real, 2, [pverp_id,nw_id], lwup_spectral_id)
status = NF_PUT_ATT_TEXT (ncid,lwup_spectral_id,'long_name', 34, 'longwave upwelling spectral fluxes')
status = NF_PUT_ATT_TEXT (ncid,lwup_spectral_id,'units', 4, 'W/m2')
if (status /= nf_noerr) call handle_err(status)

status = NF_DEF_VAR(ncid,"LWDN_SPECTRAL", nf_real, 2, [pverp_id,nw_id], lwdn_spectral_id)
status = NF_PUT_ATT_TEXT (ncid,lwdn_spectral_id,'long_name', 36, 'longwave downwelling spectral fluxes')
status = NF_PUT_ATT_TEXT (ncid,lwdn_spectral_id,'units', 4, 'W/m2')
if (status /= nf_noerr) call handle_err(status)

status = NF_DEF_VAR(ncid,"SWUP_SPECTRAL", nf_real, 2, [pverp_id,nw_id], swup_spectral_id)
status = NF_PUT_ATT_TEXT (ncid,swup_spectral_id,'long_name', 35, 'shortwave upwelling spectral fluxes')
status = NF_PUT_ATT_TEXT (ncid,swup_spectral_id,'units', 4, 'W/m2')
if (status /= nf_noerr) call handle_err(status)

status = NF_DEF_VAR(ncid,"SWDN_SPECTRAL", nf_real, 2, [pverp_id,nw_id], swdn_spectral_id)
status = NF_PUT_ATT_TEXT (ncid,swdn_spectral_id,'long_name', 37, 'shortwave downwelling spectral fluxes')
status = NF_PUT_ATT_TEXT (ncid,swdn_spectral_id,'units', 4, 'W/m2')
if (status /= nf_noerr) call handle_err(status)

status = NF_DEF_VAR(ncid,"FSDTOA", nf_real, 1, one_id, fsd_id)
status = NF_PUT_ATT_TEXT (ncid,fsd_id,'long_name', 25, 'toa incident stellar flux')
status = NF_PUT_ATT_TEXT (ncid,fsd_id,'units', 4, 'W/m2')
if (status /= nf_noerr) call handle_err(status)

status = NF_DEF_VAR(ncid,"LWHR", nf_real, 1, pver_id, lwhr_id)
status = NF_PUT_ATT_TEXT (ncid,lwhr_id,'long_name', 21, 'longwave heating rate')
status = NF_PUT_ATT_TEXT (ncid,lwhr_id,'units',3, 'K/s')
if (status /= nf_noerr) call handle_err(status)

status = NF_DEF_VAR(ncid,"SWHR", nf_real, 1, pver_id, swhr_id)
status = NF_PUT_ATT_TEXT (ncid,swhr_id,'long_name', 22, 'shortwave heating rate')
status = NF_PUT_ATT_TEXT (ncid,swhr_id,'units',3, 'K/s')
if (status /= nf_noerr) call handle_err(status)

status = NF_DEF_VAR(ncid,"PMID", nf_real, 1, pver_id, pmid_id)
status = NF_PUT_ATT_TEXT (ncid,pmid_id,'long_name', 18, 'midlayer pressures')
status = NF_PUT_ATT_TEXT (ncid,pmid_id,'units',2, 'Pa')
if (status /= nf_noerr) call handle_err(status)

status = NF_DEF_VAR(ncid,"PINT", nf_real, 1, pverp_id, pint_id)
status = NF_PUT_ATT_TEXT (ncid,pint_id,'long_name', 19, 'interface pressures')
status = NF_PUT_ATT_TEXT (ncid,pint_id,'units',2, 'Pa')
if (status /= nf_noerr) call handle_err(status)

status = NF_DEF_VAR(ncid,"TMID", nf_real, 1, pver_id, tmid_id)
status = NF_PUT_ATT_TEXT (ncid,tmid_id,'long_name', 22, 'midlayer temperatures')
status = NF_PUT_ATT_TEXT (ncid,tmid_id,'units',1, 'K')
if (status /= nf_noerr) call handle_err(status)

status = NF_DEF_VAR(ncid,"TINT", nf_real, 1, pverp_id, tint_id)
status = NF_PUT_ATT_TEXT (ncid,tint_id,'long_name', 23, 'interface temperatures')
status = NF_PUT_ATT_TEXT (ncid,tint_id,'units',1, 'K')
if (status /= nf_noerr) call handle_err(status)

status = NF_DEF_VAR(ncid,"ZINT", nf_real, 1, pverp_id, zint_id)
status = NF_PUT_ATT_TEXT (ncid,zint_id,'long_name', 30, 'interface geopotential heights')
status = NF_PUT_ATT_TEXT (ncid,zint_id,'units',1, 'm')
if (status /= nf_noerr) call handle_err(status)

status = NF_DEF_VAR(ncid,"H2OMMR", nf_real, 1, pver_id, h2o_id)
status = NF_PUT_ATT_TEXT (ncid,h2o_id,'long_name', 21, 'H2O specific humidity')
status = NF_PUT_ATT_TEXT (ncid,h2o_id,'units',15, 'kg/kg_moist_air')
if (status /= nf_noerr) call handle_err(status)

status = NF_DEF_VAR(ncid,"CO2MMR", nf_real, 1, pver_id, co2_id)
status = NF_PUT_ATT_TEXT (ncid,co2_id,'long_name', 21, 'CO2 mass mixing ratio')
status = NF_PUT_ATT_TEXT (ncid,co2_id,'units',13, 'kg/kg_dry_air')
if (status /= nf_noerr) call handle_err(status)

status = NF_DEF_VAR(ncid,"CH4MMR", nf_real, 1, pver_id, ch4_id)
status = NF_PUT_ATT_TEXT (ncid,ch4_id,'long_name', 21, 'CH4 mass mixing ratio')
status = NF_PUT_ATT_TEXT (ncid,ch4_id,'units',13, 'kg/kg_dry_air')
if (status /= nf_noerr) call handle_err(status)

status = NF_DEF_VAR(ncid,"C2H6MMR", nf_real, 1, pver_id, c2h6_id)
status = NF_PUT_ATT_TEXT (ncid,c2h6_id,'long_name', 22, 'C2H6 mass mixing ratio')
status = NF_PUT_ATT_TEXT (ncid,c2h6_id,'units',13, 'kg/kg_dry_air')
if (status /= nf_noerr) call handle_err(status)

status = NF_DEF_VAR(ncid,"NH3MMR", nf_real, 1, pver_id, nh3_id)
status = NF_PUT_ATT_TEXT (ncid,nh3_id,'long_name', 21, 'NH3 mass mixing ratio')
status = NF_PUT_ATT_TEXT (ncid,nh3_id,'units',13, 'kg/kg_dry_air')
if (status /= nf_noerr) call handle_err(status)

status = NF_DEF_VAR(ncid,"COMMR", nf_real, 1, pver_id, co_id)
status = NF_PUT_ATT_TEXT (ncid,co_id,'long_name', 20, 'CO mass mixing ratio')
status = NF_PUT_ATT_TEXT (ncid,co_id,'units',13, 'kg/kg_dry_air')
if (status /= nf_noerr) call handle_err(status)

status = NF_DEF_VAR(ncid,"O2MMR", nf_real, 1, pver_id, o2_id)
status = NF_PUT_ATT_TEXT (ncid,o2_id,'long_name', 20, 'O2 mass mixing ratio')
status = NF_PUT_ATT_TEXT (ncid,o2_id,'units',13, 'kg/kg_dry_air')
if (status /= nf_noerr) call handle_err(status)

status = NF_DEF_VAR(ncid,"O3MMR", nf_real, 1, pver_id, o3_id)
status = NF_PUT_ATT_TEXT (ncid,o3_id,'long_name', 20, 'O3 mass mixing ratio')
status = NF_PUT_ATT_TEXT (ncid,o3_id,'units',13, 'kg/kg_dry_air')
if (status /= nf_noerr) call handle_err(status)

status = NF_DEF_VAR(ncid,"N2MMR", nf_real, 1, pver_id, n2_id)
status = NF_PUT_ATT_TEXT (ncid,n2_id,'long_name', 20, 'N2 mass mixing ratio')
status = NF_PUT_ATT_TEXT (ncid,n2_id,'units',13, 'kg/kg_dry_air')
if (status /= nf_noerr) call handle_err(status)

status = NF_DEF_VAR(ncid,"H2MMR", nf_real, 1, pver_id, h2_id)
status = NF_PUT_ATT_TEXT (ncid,h2_id,'long_name', 20, 'H2 mass mixing ratio')
status = NF_PUT_ATT_TEXT (ncid,h2_id,'units',13, 'kg/kg_dry_air')
if (status /= nf_noerr) call handle_err(status)

status = NF_DEF_VAR(ncid,"MWDRY", nf_real, 1, one_id, mw_id)
status = NF_PUT_ATT_TEXT (ncid,mw_id,'long_name', 27, 'molecular weight of dry air')
status = NF_PUT_ATT_TEXT (ncid,mw_id,'units',3, 'AMU')
if (status /= nf_noerr) call handle_err(status)

status = NF_DEF_VAR(ncid,"CPDRY", nf_real, 1, one_id, cp_id)
status = NF_PUT_ATT_TEXT (ncid,cp_id,'long_name', 21, 'specific heat of dry air')
status = NF_PUT_ATT_TEXT (ncid,cp_id,'units',6, 'J/kg K')
if (status /= nf_noerr) call handle_err(status)


! global attributes: provenance stamp for the three runtime-config values
status = NF_PUT_ATT_TEXT  (ncid, NF_GLOBAL, 'solar_file',     len_trim(solar_file), trim(solar_file))
if (status /= nf_noerr) call handle_err(status)
status = NF_PUT_ATT_DOUBLE(ncid, NF_GLOBAL, 'shr_const_scon', nf_double, 1, [scon])
if (status /= nf_noerr) call handle_err(status)
status = NF_PUT_ATT_DOUBLE(ncid, NF_GLOBAL, 'exo_g',          nf_double, 1, [exo_g])
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
status = NF_PUT_VAR_double(ncid,nh3_id,nh3_out)
if (status /= nf_noerr) call handle_err(status)
status = NF_PUT_VAR_double(ncid,co_id,co_out)
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



end module io
