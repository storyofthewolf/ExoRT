module io
!----------------------------------------------------------------------
! Consolidated 1-D input/output module for ExoRT (standalone driver).
!   read_namelist      - read user_nl_exort runtime config (1-D only)
!   input_profile      - read RTprofile_in.nc into column_state_t array
!   print_diagnostics  - echo primary fluxes to stdout (1-D only)
!   output_data        - write RTprofile_out.nc from state/result arrays
!   handle_err         - NetCDF error handler
! Not used by the CESM/3-D path (which has its own I/O).
!
! Stage E1: this module holds NO column state. input_profile allocates
! and fills an array of column_state_t; output_data writes from arrays
! of column_state_t/column_result_t. RTprofile_in.nc may carry an
! optional 'ncol' dimension: absent = classic single-column file
! (bit-for-bit the pre-E1 path); present = every variable carries a
! trailing column dimension (Fortran order — e.g. tmid(pver,ncol),
! ts(ncol), carmammr(pver,nelem,nbin,ncol)) and RTprofile_out.nc gains
! the matching 'ncol' dimension on every variable. tools/stackColumns.py
! builds multi-column inputs from single-column files.
!----------------------------------------------------------------------

use shr_kind_mod,       only: r8 => shr_kind_r8
use shr_const_mod
use ppgrid
use ioFileMod
use physconst
use radgrid
use exoplanet_mod,      only: solar_file, exo_g, shr_const_scon, do_exo_clouds, do_exo_haze, &
                              mcica_percol_seed
use exort_column_mod,   only: column_state_t, column_result_t, zero_column_state

implicit none
public

! Namelist for 1-D runtime config; defaults are the module-variable
! initializers in exoplanet_mod. Read by read_namelist() below.
namelist /exort_config/ solar_file, shr_const_scon, exo_g, do_exo_clouds, do_exo_haze, &
                        mcica_percol_seed

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
  write(*,*) 'mcica_percol_seed = ', mcica_percol_seed
  write(*,*) '============================================'

end subroutine read_namelist


subroutine input_profile(states)
!----------------------------------------------------------------------
! Read RTprofile_in.nc into an array of column_state_t (allocated here
! to the file's column count; 1 when the 'ncol' dimension is absent).
! Gas/cloud/haze variables remain optional: absent variables stay zero
! with a clean diagnostic line, never a NetCDF error.
!----------------------------------------------------------------------
  type(column_state_t), allocatable, intent(out) :: states(:)

  include 'netcdf.inc'

  character(len=256) :: input_file
  character(len=256) :: locfn
  integer :: ncid
  integer :: pverp_id, pver_id, ncol_id, npverp, npver
  integer :: ncol, ic, ret

  ! whole-variable read buffers (trailing column dimension)
  real(r8), allocatable :: b0(:)        ! scalars        (ncol)
  real(r8), allocatable :: bm(:,:)      ! mid-layer      (pver,ncol)
  real(r8), allocatable :: bi(:,:)      ! interface      (pverp,ncol)
  real(r8), allocatable :: bc(:,:,:,:)  ! CARMA          (pver,nelem,nbin,ncol)

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

  ! Optional column dimension (absent = single-column file)
  ret = nf_inq_dimid(ncid, 'ncol', ncol_id)
  if (ret == NF_NOERR) then
    call wrap_inq_dimlen(ncid, ncol_id, ncol)
    write(*,*) "columns: ", ncol
  else
    ncol = 1
  endif

  allocate(states(ncol))
  allocate(b0(ncol))
  allocate(bm(pver,ncol))
  allocate(bi(pverp,ncol))
  allocate(bc(pver,nelem_carma,nbin_carma,ncol))

  call zero_column_state(states)

  ! Required: P, T, Z
  call req_scalar('ts');    states(:)%ts = b0(:)
  call req_scalar('ps');    states(:)%ps = b0(:)
  call req_mid('tmid');     do ic=1,ncol; states(ic)%tmid = bm(:,ic); enddo
  call req_int('tint');     do ic=1,ncol; states(ic)%tint = bi(:,ic); enddo
  call req_mid('pmid');     do ic=1,ncol; states(ic)%pmid = bm(:,ic); enddo
  call req_mid('pdel');     do ic=1,ncol; states(ic)%pdel = bm(:,ic); enddo
  call req_int('pint');     do ic=1,ncol; states(ic)%pint = bi(:,ic); enddo
  call req_int('zint');     do ic=1,ncol; states(ic)%zint = bi(:,ic); enddo

  ! Required: albedos, cosz, mw, cp
  call req_scalar('asdir');  states(:)%asdir  = b0(:)
  call req_scalar('asdif');  states(:)%asdif  = b0(:)
  call req_scalar('aldir');  states(:)%aldir  = b0(:)
  call req_scalar('aldif');  states(:)%aldif  = b0(:)
  call req_scalar('coszrs'); states(:)%coszrs = b0(:)
  call req_scalar('mw');     states(:)%mwdry  = b0(:)
  call req_scalar('cp');     states(:)%cpdry  = b0(:)

  ! Optional: surface thermal emissivity (default 1.0 if absent)
  states(:)%srf_emiss = 1.0
  if (opt_scalar('srf_emiss')) then
    write(6,*) "  srf_emiss:  found"
    states(:)%srf_emiss = b0(:)
  else
    write(6,*) "  srf_emiss:  not found, set to 1.0"
  end if

  ! Optional gases: zero-initialized above; read only if present
  write(6,*) "--- gas species in input file ---"
  if (opt_mid('h2ommr', 'h2ommr:  ')) then
    do ic=1,ncol; states(ic)%h2ommr = bm(:,ic); enddo
  endif
  if (opt_mid('co2mmr', 'co2mmr:  ')) then
    do ic=1,ncol; states(ic)%co2mmr = bm(:,ic); enddo
  endif
  if (opt_mid('ch4mmr', 'ch4mmr:  ')) then
    do ic=1,ncol; states(ic)%ch4mmr = bm(:,ic); enddo
  endif
  if (opt_mid('c2h6mmr', 'c2h6mmr: ')) then
    do ic=1,ncol; states(ic)%c2h6mmr = bm(:,ic); enddo
  endif
  if (opt_mid('nh3mmr', 'nh3mmr:  ')) then
    do ic=1,ncol; states(ic)%nh3mmr = bm(:,ic); enddo
  endif
  if (opt_mid('commr', 'commr:   ')) then
    do ic=1,ncol; states(ic)%commr = bm(:,ic); enddo
  endif
  if (opt_mid('o2mmr', 'o2mmr:   ')) then
    do ic=1,ncol; states(ic)%o2mmr = bm(:,ic); enddo
  endif
  if (opt_mid('o3mmr', 'o3mmr:   ')) then
    do ic=1,ncol; states(ic)%o3mmr = bm(:,ic); enddo
  endif
  if (opt_mid('n2mmr', 'n2mmr:   ')) then
    do ic=1,ncol; states(ic)%n2mmr = bm(:,ic); enddo
  endif
  if (opt_mid('h2mmr', 'h2mmr:   ')) then
    do ic=1,ncol; states(ic)%h2mmr = bm(:,ic); enddo
  endif

  ! Optional: cloud properties (zero-initialized above; no clouds if absent)
  write(6,*) "--- cloud variables in input file ---"
  if (opt_mid('cicewp', 'cicewp:  ')) then
    do ic=1,ncol; states(ic)%cicewp = bm(:,ic); enddo
  endif
  if (opt_mid('cliqwp', 'cliqwp:  ')) then
    do ic=1,ncol; states(ic)%cliqwp = bm(:,ic); enddo
  endif
  if (opt_mid('rei', 'rei:     ')) then
    do ic=1,ncol; states(ic)%rei = bm(:,ic); enddo
  endif
  if (opt_mid('rel', 'rel:     ')) then
    do ic=1,ncol; states(ic)%rel = bm(:,ic); enddo
  endif
  if (opt_mid('cfrc', 'cfrc:    ')) then
    do ic=1,ncol; states(ic)%cfrc = bm(:,ic); enddo
  endif
  if (opt_mid('cicewp_co2', 'cicewp_co2:  ')) then
    do ic=1,ncol; states(ic)%cicewp_co2 = bm(:,ic); enddo
  endif
  if (opt_mid('rei_co2', 'rei_co2:     ')) then
    do ic=1,ncol; states(ic)%rei_co2 = bm(:,ic); enddo
  endif

  ! Optional: CARMA haze binwise mass mixing ratio (zero-initialized above)
  write(6,*) "--- aerosol variables in input file ---"
  if (opt_carma('carmammr', 'carmammr:  ')) then
    do ic=1,ncol; states(ic)%carmammr = bc(:,:,:,ic); enddo
  endif

  write(6,*) "---------------------------------"

contains

  subroutine req_scalar(name)
    character(len=*), intent(in) :: name
    integer :: vid
    call wrap_inq_varid(ncid, name, vid)
    call wrap_get_var_realx(ncid, vid, b0)
  end subroutine req_scalar

  subroutine req_mid(name)
    character(len=*), intent(in) :: name
    integer :: vid
    call wrap_inq_varid(ncid, name, vid)
    call wrap_get_var_realx(ncid, vid, bm)
  end subroutine req_mid

  subroutine req_int(name)
    character(len=*), intent(in) :: name
    integer :: vid
    call wrap_inq_varid(ncid, name, vid)
    call wrap_get_var_realx(ncid, vid, bi)
  end subroutine req_int

  logical function opt_scalar(name)
    character(len=*), intent(in) :: name
    integer :: vid, ret2
    ret2 = nf_inq_varid(ncid, name, vid)
    opt_scalar = (ret2 == NF_NOERR)
    if (opt_scalar) call wrap_get_var_realx(ncid, vid, b0)
  end function opt_scalar

  logical function opt_mid(name, label)
    character(len=*), intent(in) :: name, label
    integer :: vid, ret2
    ret2 = nf_inq_varid(ncid, name, vid)
    opt_mid = (ret2 == NF_NOERR)
    if (opt_mid) then
      write(6,*) "  "//label//"found"
      call wrap_get_var_realx(ncid, vid, bm)
    else
      write(6,*) "  "//label//"not found, set to zero"
    endif
  end function opt_mid

  logical function opt_carma(name, label)
    character(len=*), intent(in) :: name, label
    integer :: vid, ret2
    ret2 = nf_inq_varid(ncid, name, vid)
    opt_carma = (ret2 == NF_NOERR)
    if (opt_carma) then
      write(6,*) "  "//label//"found"
      call wrap_get_var_realx(ncid, vid, bc)
    else
      write(6,*) "  "//label//"not found, set to zero"
    endif
  end function opt_carma

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


subroutine output_data ( states, results )
!----------------------------------------------------------------------
! Write RTprofile_out.nc. Single column (size(states)==1): identical
! layout to the pre-E1 file — no 'ncol' dimension, scalars on the ONE
! dimension. Multi-column: every variable gains a trailing 'ncol'
! dimension and the ONE-dimension scalars become ncol vectors.
! Heating rates are converted K/s -> K/day here (SHR_CONST_CSEC).
!----------------------------------------------------------------------
type(column_state_t),  intent(in) :: states(:)
type(column_result_t), intent(in) :: results(:)

include 'netcdf.inc'

integer :: ncol, ic
integer :: ncid, status, pverp_id, pver_id, one_id, nw_id, ncol_id
integer :: lwup_id, lwdn_id, swup_id, swdn_id, fsd_id
integer :: lwup_spectral_id, lwdn_spectral_id, swup_spectral_id, swdn_spectral_id
integer :: lwhr_id, swhr_id
integer :: pmid_id, pint_id, tmid_id
integer :: tint_id, zint_id
integer :: h2o_id, co2_id, ch4_id, c2h6_id, nh3_id, co_id, o2_id, o3_id, n2_id, h2_id, mw_id, cp_id

! dimension signatures: (count, dimids) per variable class
integer :: nd_scl, nd_mid, nd_int, nd_spec
integer :: d_scl(2), d_mid(2), d_int(2), d_spec(3)

! staging buffers (trailing column dimension)
real(r8), allocatable :: w0(:)        ! scalars   (ncol)
real(r8), allocatable :: wm(:,:)      ! mid-layer (pver,ncol)
real(r8), allocatable :: wi(:,:)      ! interface (pverp,ncol)
real(r8), allocatable :: ws(:,:,:)    ! spectral  (pverp,ntot_wavlnrng,ncol)

ncol = size(states)
allocate(w0(ncol))
allocate(wm(pver,ncol))
allocate(wi(pverp,ncol))
allocate(ws(pverp,ntot_wavlnrng,ncol))

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

if (ncol > 1) then
  status = NF_DEF_DIM(ncid,"ncol", ncol, ncol_id)
  if (status /= nf_noerr) call handle_err(status)
  nd_scl = 1;  d_scl = [ncol_id, 0]
  nd_mid = 2;  d_mid = [pver_id, ncol_id]
  nd_int = 2;  d_int = [pverp_id, ncol_id]
  nd_spec = 3; d_spec = [pverp_id, nw_id, ncol_id]
else
  nd_scl = 1;  d_scl = [one_id, 0]
  nd_mid = 1;  d_mid = [pver_id, 0]
  nd_int = 1;  d_int = [pverp_id, 0]
  nd_spec = 2; d_spec = [pverp_id, nw_id, 0]
endif

! defining variables
status = NF_DEF_VAR(ncid,"LWUP", nf_real, nd_int, d_int, lwup_id)
status = NF_PUT_ATT_TEXT (ncid,lwup_id,'long_name', 24, 'longwave upwelling flux')
status = NF_PUT_ATT_TEXT (ncid,lwup_id,'units', 4, 'W/m2')
if (status /= nf_noerr) call handle_err(status)

status = NF_DEF_VAR(ncid,"LWDN", nf_real, nd_int, d_int, lwdn_id)
status = NF_PUT_ATT_TEXT (ncid,lwdn_id,'long_name', 26, 'longwave downwelling flux')
status = NF_PUT_ATT_TEXT (ncid,lwdn_id,'units', 4, 'W/m2')
if (status /= nf_noerr) call handle_err(status)

status = NF_DEF_VAR(ncid,"SWUP", nf_real, nd_int, d_int, swup_id)
status = NF_PUT_ATT_TEXT (ncid,swup_id,'long_name', 25, 'shortwave upwelling flux')
status = NF_PUT_ATT_TEXT (ncid,swup_id,'units', 4, 'W/m2')
if (status /= nf_noerr) call handle_err(status)

status = NF_DEF_VAR(ncid,"SWDN", nf_real, nd_int, d_int, swdn_id)
status = NF_PUT_ATT_TEXT (ncid,swdn_id,'long_name', 27, 'shortwave downwelling flux')
status = NF_PUT_ATT_TEXT (ncid,swdn_id,'units', 4, 'W/m2')
if (status /= nf_noerr) call handle_err(status)

status = NF_DEF_VAR(ncid,"LWUP_SPECTRAL", nf_real, nd_spec, d_spec, lwup_spectral_id)
status = NF_PUT_ATT_TEXT (ncid,lwup_spectral_id,'long_name', 34, 'longwave upwelling spectral fluxes')
status = NF_PUT_ATT_TEXT (ncid,lwup_spectral_id,'units', 4, 'W/m2')
if (status /= nf_noerr) call handle_err(status)

status = NF_DEF_VAR(ncid,"LWDN_SPECTRAL", nf_real, nd_spec, d_spec, lwdn_spectral_id)
status = NF_PUT_ATT_TEXT (ncid,lwdn_spectral_id,'long_name', 36, 'longwave downwelling spectral fluxes')
status = NF_PUT_ATT_TEXT (ncid,lwdn_spectral_id,'units', 4, 'W/m2')
if (status /= nf_noerr) call handle_err(status)

status = NF_DEF_VAR(ncid,"SWUP_SPECTRAL", nf_real, nd_spec, d_spec, swup_spectral_id)
status = NF_PUT_ATT_TEXT (ncid,swup_spectral_id,'long_name', 35, 'shortwave upwelling spectral fluxes')
status = NF_PUT_ATT_TEXT (ncid,swup_spectral_id,'units', 4, 'W/m2')
if (status /= nf_noerr) call handle_err(status)

status = NF_DEF_VAR(ncid,"SWDN_SPECTRAL", nf_real, nd_spec, d_spec, swdn_spectral_id)
status = NF_PUT_ATT_TEXT (ncid,swdn_spectral_id,'long_name', 37, 'shortwave downwelling spectral fluxes')
status = NF_PUT_ATT_TEXT (ncid,swdn_spectral_id,'units', 4, 'W/m2')
if (status /= nf_noerr) call handle_err(status)

status = NF_DEF_VAR(ncid,"FSDTOA", nf_real, nd_scl, d_scl, fsd_id)
status = NF_PUT_ATT_TEXT (ncid,fsd_id,'long_name', 25, 'toa incident stellar flux')
status = NF_PUT_ATT_TEXT (ncid,fsd_id,'units', 4, 'W/m2')
if (status /= nf_noerr) call handle_err(status)

status = NF_DEF_VAR(ncid,"LWHR", nf_real, nd_mid, d_mid, lwhr_id)
status = NF_PUT_ATT_TEXT (ncid,lwhr_id,'long_name', 21, 'longwave heating rate')
status = NF_PUT_ATT_TEXT (ncid,lwhr_id,'units',3, 'K/s')
if (status /= nf_noerr) call handle_err(status)

status = NF_DEF_VAR(ncid,"SWHR", nf_real, nd_mid, d_mid, swhr_id)
status = NF_PUT_ATT_TEXT (ncid,swhr_id,'long_name', 22, 'shortwave heating rate')
status = NF_PUT_ATT_TEXT (ncid,swhr_id,'units',3, 'K/s')
if (status /= nf_noerr) call handle_err(status)

status = NF_DEF_VAR(ncid,"PMID", nf_real, nd_mid, d_mid, pmid_id)
status = NF_PUT_ATT_TEXT (ncid,pmid_id,'long_name', 18, 'midlayer pressures')
status = NF_PUT_ATT_TEXT (ncid,pmid_id,'units',2, 'Pa')
if (status /= nf_noerr) call handle_err(status)

status = NF_DEF_VAR(ncid,"PINT", nf_real, nd_int, d_int, pint_id)
status = NF_PUT_ATT_TEXT (ncid,pint_id,'long_name', 19, 'interface pressures')
status = NF_PUT_ATT_TEXT (ncid,pint_id,'units',2, 'Pa')
if (status /= nf_noerr) call handle_err(status)

status = NF_DEF_VAR(ncid,"TMID", nf_real, nd_mid, d_mid, tmid_id)
status = NF_PUT_ATT_TEXT (ncid,tmid_id,'long_name', 22, 'midlayer temperatures')
status = NF_PUT_ATT_TEXT (ncid,tmid_id,'units',1, 'K')
if (status /= nf_noerr) call handle_err(status)

status = NF_DEF_VAR(ncid,"TINT", nf_real, nd_int, d_int, tint_id)
status = NF_PUT_ATT_TEXT (ncid,tint_id,'long_name', 23, 'interface temperatures')
status = NF_PUT_ATT_TEXT (ncid,tint_id,'units',1, 'K')
if (status /= nf_noerr) call handle_err(status)

status = NF_DEF_VAR(ncid,"ZINT", nf_real, nd_int, d_int, zint_id)
status = NF_PUT_ATT_TEXT (ncid,zint_id,'long_name', 30, 'interface geopotential heights')
status = NF_PUT_ATT_TEXT (ncid,zint_id,'units',1, 'm')
if (status /= nf_noerr) call handle_err(status)

status = NF_DEF_VAR(ncid,"H2OMMR", nf_real, nd_mid, d_mid, h2o_id)
status = NF_PUT_ATT_TEXT (ncid,h2o_id,'long_name', 21, 'H2O specific humidity')
status = NF_PUT_ATT_TEXT (ncid,h2o_id,'units',15, 'kg/kg_moist_air')
if (status /= nf_noerr) call handle_err(status)

status = NF_DEF_VAR(ncid,"CO2MMR", nf_real, nd_mid, d_mid, co2_id)
status = NF_PUT_ATT_TEXT (ncid,co2_id,'long_name', 21, 'CO2 mass mixing ratio')
status = NF_PUT_ATT_TEXT (ncid,co2_id,'units',13, 'kg/kg_dry_air')
if (status /= nf_noerr) call handle_err(status)

status = NF_DEF_VAR(ncid,"CH4MMR", nf_real, nd_mid, d_mid, ch4_id)
status = NF_PUT_ATT_TEXT (ncid,ch4_id,'long_name', 21, 'CH4 mass mixing ratio')
status = NF_PUT_ATT_TEXT (ncid,ch4_id,'units',13, 'kg/kg_dry_air')
if (status /= nf_noerr) call handle_err(status)

status = NF_DEF_VAR(ncid,"C2H6MMR", nf_real, nd_mid, d_mid, c2h6_id)
status = NF_PUT_ATT_TEXT (ncid,c2h6_id,'long_name', 22, 'C2H6 mass mixing ratio')
status = NF_PUT_ATT_TEXT (ncid,c2h6_id,'units',13, 'kg/kg_dry_air')
if (status /= nf_noerr) call handle_err(status)

status = NF_DEF_VAR(ncid,"NH3MMR", nf_real, nd_mid, d_mid, nh3_id)
status = NF_PUT_ATT_TEXT (ncid,nh3_id,'long_name', 21, 'NH3 mass mixing ratio')
status = NF_PUT_ATT_TEXT (ncid,nh3_id,'units',13, 'kg/kg_dry_air')
if (status /= nf_noerr) call handle_err(status)

status = NF_DEF_VAR(ncid,"COMMR", nf_real, nd_mid, d_mid, co_id)
status = NF_PUT_ATT_TEXT (ncid,co_id,'long_name', 20, 'CO mass mixing ratio')
status = NF_PUT_ATT_TEXT (ncid,co_id,'units',13, 'kg/kg_dry_air')
if (status /= nf_noerr) call handle_err(status)

status = NF_DEF_VAR(ncid,"O2MMR", nf_real, nd_mid, d_mid, o2_id)
status = NF_PUT_ATT_TEXT (ncid,o2_id,'long_name', 20, 'O2 mass mixing ratio')
status = NF_PUT_ATT_TEXT (ncid,o2_id,'units',13, 'kg/kg_dry_air')
if (status /= nf_noerr) call handle_err(status)

status = NF_DEF_VAR(ncid,"O3MMR", nf_real, nd_mid, d_mid, o3_id)
status = NF_PUT_ATT_TEXT (ncid,o3_id,'long_name', 20, 'O3 mass mixing ratio')
status = NF_PUT_ATT_TEXT (ncid,o3_id,'units',13, 'kg/kg_dry_air')
if (status /= nf_noerr) call handle_err(status)

status = NF_DEF_VAR(ncid,"N2MMR", nf_real, nd_mid, d_mid, n2_id)
status = NF_PUT_ATT_TEXT (ncid,n2_id,'long_name', 20, 'N2 mass mixing ratio')
status = NF_PUT_ATT_TEXT (ncid,n2_id,'units',13, 'kg/kg_dry_air')
if (status /= nf_noerr) call handle_err(status)

status = NF_DEF_VAR(ncid,"H2MMR", nf_real, nd_mid, d_mid, h2_id)
status = NF_PUT_ATT_TEXT (ncid,h2_id,'long_name', 20, 'H2 mass mixing ratio')
status = NF_PUT_ATT_TEXT (ncid,h2_id,'units',13, 'kg/kg_dry_air')
if (status /= nf_noerr) call handle_err(status)

status = NF_DEF_VAR(ncid,"MWDRY", nf_real, nd_scl, d_scl, mw_id)
status = NF_PUT_ATT_TEXT (ncid,mw_id,'long_name', 27, 'molecular weight of dry air')
status = NF_PUT_ATT_TEXT (ncid,mw_id,'units',3, 'AMU')
if (status /= nf_noerr) call handle_err(status)

status = NF_DEF_VAR(ncid,"CPDRY", nf_real, nd_scl, d_scl, cp_id)
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

! put variables (stage each field across columns, then one whole-var put)
do ic=1,ncol; wi(:,ic) = results(ic)%lw_upflux; enddo
status = NF_PUT_VAR_double(ncid,lwup_id,wi)
if (status /= nf_noerr) call handle_err(status)
do ic=1,ncol; wi(:,ic) = results(ic)%lw_dnflux; enddo
status = NF_PUT_VAR_double(ncid,lwdn_id,wi)
if (status /= nf_noerr) call handle_err(status)
do ic=1,ncol; wi(:,ic) = results(ic)%sw_upflux; enddo
status = NF_PUT_VAR_double(ncid,swup_id,wi)
if (status /= nf_noerr) call handle_err(status)
do ic=1,ncol; wi(:,ic) = results(ic)%sw_dnflux; enddo
status = NF_PUT_VAR_double(ncid,swdn_id,wi)
if (status /= nf_noerr) call handle_err(status)
do ic=1,ncol; ws(:,:,ic) = results(ic)%lw_upflux_spectral; enddo
status = NF_PUT_VAR_double(ncid,lwup_spectral_id,ws)
if (status /= nf_noerr) call handle_err(status)
do ic=1,ncol; ws(:,:,ic) = results(ic)%lw_dnflux_spectral; enddo
status = NF_PUT_VAR_double(ncid,lwdn_spectral_id,ws)
if (status /= nf_noerr) call handle_err(status)
do ic=1,ncol; ws(:,:,ic) = results(ic)%sw_upflux_spectral; enddo
status = NF_PUT_VAR_double(ncid,swup_spectral_id,ws)
if (status /= nf_noerr) call handle_err(status)
do ic=1,ncol; ws(:,:,ic) = results(ic)%sw_dnflux_spectral; enddo
status = NF_PUT_VAR_double(ncid,swdn_spectral_id,ws)
if (status /= nf_noerr) call handle_err(status)
w0(:) = results(:)%sol_toa
status = NF_PUT_VAR_double(ncid,fsd_id,w0)
if (status /= nf_noerr) call handle_err(status)
do ic=1,ncol; wm(:,ic) = results(ic)%lw_dtdt * SHR_CONST_CSEC; enddo
status = NF_PUT_VAR_double(ncid,lwhr_id,wm)
if (status /= nf_noerr) call handle_err(status)
do ic=1,ncol; wm(:,ic) = results(ic)%sw_dtdt * SHR_CONST_CSEC; enddo
status = NF_PUT_VAR_double(ncid,swhr_id,wm)
if (status /= nf_noerr) call handle_err(status)
do ic=1,ncol; wm(:,ic) = states(ic)%pmid; enddo
status = NF_PUT_VAR_double(ncid,pmid_id,wm)
if (status /= nf_noerr) call handle_err(status)
do ic=1,ncol; wi(:,ic) = states(ic)%pint; enddo
status = NF_PUT_VAR_double(ncid,pint_id,wi)
if (status /= nf_noerr) call handle_err(status)
do ic=1,ncol; wm(:,ic) = states(ic)%tmid; enddo
status = NF_PUT_VAR_double(ncid,tmid_id,wm)
if (status /= nf_noerr) call handle_err(status)
do ic=1,ncol; wi(:,ic) = states(ic)%tint; enddo
status = NF_PUT_VAR_double(ncid,tint_id,wi)
if (status /= nf_noerr) call handle_err(status)
do ic=1,ncol; wi(:,ic) = states(ic)%zint; enddo
status = NF_PUT_VAR_double(ncid,zint_id,wi)
if (status /= nf_noerr) call handle_err(status)
do ic=1,ncol; wm(:,ic) = states(ic)%h2ommr; enddo
status = NF_PUT_VAR_double(ncid,h2o_id,wm)
if (status /= nf_noerr) call handle_err(status)
do ic=1,ncol; wm(:,ic) = states(ic)%co2mmr; enddo
status = NF_PUT_VAR_double(ncid,co2_id,wm)
if (status /= nf_noerr) call handle_err(status)
do ic=1,ncol; wm(:,ic) = states(ic)%ch4mmr; enddo
status = NF_PUT_VAR_double(ncid,ch4_id,wm)
if (status /= nf_noerr) call handle_err(status)
do ic=1,ncol; wm(:,ic) = states(ic)%c2h6mmr; enddo
status = NF_PUT_VAR_double(ncid,c2h6_id,wm)
if (status /= nf_noerr) call handle_err(status)
do ic=1,ncol; wm(:,ic) = states(ic)%nh3mmr; enddo
status = NF_PUT_VAR_double(ncid,nh3_id,wm)
if (status /= nf_noerr) call handle_err(status)
do ic=1,ncol; wm(:,ic) = states(ic)%commr; enddo
status = NF_PUT_VAR_double(ncid,co_id,wm)
if (status /= nf_noerr) call handle_err(status)
do ic=1,ncol; wm(:,ic) = states(ic)%o2mmr; enddo
status = NF_PUT_VAR_double(ncid,o2_id,wm)
if (status /= nf_noerr) call handle_err(status)
do ic=1,ncol; wm(:,ic) = states(ic)%o3mmr; enddo
status = NF_PUT_VAR_double(ncid,o3_id,wm)
if (status /= nf_noerr) call handle_err(status)
do ic=1,ncol; wm(:,ic) = states(ic)%n2mmr; enddo
status = NF_PUT_VAR_double(ncid,n2_id,wm)
if (status /= nf_noerr) call handle_err(status)
do ic=1,ncol; wm(:,ic) = states(ic)%h2mmr; enddo
status = NF_PUT_VAR_double(ncid,h2_id,wm)
if (status /= nf_noerr) call handle_err(status)
w0(:) = states(:)%mwdry
status = NF_PUT_VAR_double(ncid,mw_id,w0)
if (status /= nf_noerr) call handle_err(status)
w0(:) = states(:)%cpdry
status = NF_PUT_VAR_double(ncid,cp_id,w0)
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
