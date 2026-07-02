
module exort_lib_mod
!----------------------------------------------------------------------
! C-callable library API for ExoRT (Stage D). Exposes the same physics
! as the standalone 1-D driver (main.F90) behind four bind(c) entry
! points, so C and Python callers can run columns in-process without
! NetCDF profile files:
!
!   exort_get_dims(pver, pverp, nwave, nelem, nbin)   dims for array sizing
!   exort_init(data_root, solar_file, scon, g,
!              do_clouds, do_haze)                     load all tables (once)
!   exort_run_column(state, result)                    one column
!   exort_run_columns(n, states, results)              n columns, serial loop
!   exort_finalize()                                   refuse further runs
!
! All functions return 0 on success (see error codes below). Fatal data
! errors (missing k-files etc.) still `stop` inside the legacy readers
! and terminate the process — same behavior as the executable.
!
! Configuration mirrors the user_nl_exort namelist: pass an empty string
! for data_root/solar_file or a non-positive scon/g to keep the compiled
! defaults. exo_pver remains compile-time; states carry pver levels.
!
! Thread-safety contract: after exort_init returns, all tables are
! read-only. exort_run_column itself is NOT yet thread-safe (module
! scratch state in the RT kernels; per-column MCICA seeding lands in
! Stage E) — call it from one thread at a time.
!----------------------------------------------------------------------

use iso_c_binding
use shr_kind_mod,          only: r8 => shr_kind_r8
use ppgrid,                only: pver, pverp
use radgrid,               only: ntot_wavlnrng, nelem_carma, nbin_carma
use exoplanet_mod,         only: solar_file, exo_g, shr_const_scon, &
                                 do_exo_clouds, do_exo_haze
use sys_rootdir,           only: exort_rootdir
use physconst,             only: scon, physconst_setgas
use shr_const_mod,         only: SHR_CONST_G
use initialize_rad_mod_1D, only: initialize_kcoeff, initialize_solar, &
                                 initialize_cldopts, initialize_hazeopts, &
                                 initialize_radbuffer
use exo_init_ref,          only: init_ref
use exo_model_specific,    only: init_model_specific
use exo_radiation_mod,     only: init_planck, aerad_driver
use exort_column_mod,      only: column_state_t, column_result_t

implicit none
private

! error codes
integer(c_int), parameter :: EXORT_OK              = 0
integer(c_int), parameter :: EXORT_ALREADY_INIT    = 1
integer(c_int), parameter :: EXORT_NOT_INITIALIZED = 2
integer(c_int), parameter :: EXORT_BAD_ARGS        = 3

logical :: initialized = .false.
logical :: finalized   = .false.

public :: exort_get_dims, exort_init, exort_run_column, &
          exort_run_columns, exort_finalize

contains

!============================================================================

function exort_get_dims(n_pver, n_pverp, n_wavlnrng, n_elem_carma, &
                        n_bin_carma) bind(c, name='exort_get_dims') result(ierr)
  integer(c_int), intent(out) :: n_pver, n_pverp, n_wavlnrng
  integer(c_int), intent(out) :: n_elem_carma, n_bin_carma
  integer(c_int) :: ierr

  n_pver       = pver
  n_pverp      = pverp
  n_wavlnrng   = ntot_wavlnrng
  n_elem_carma = nelem_carma
  n_bin_carma  = nbin_carma
  ierr = EXORT_OK

end function exort_get_dims

!============================================================================

function exort_init(data_root, solar_file_in, scon_in, g_in, &
                    do_clouds, do_haze) bind(c, name='exort_init') result(ierr)
!
! Load every table (k-coefficients, CIA, continuum, stellar spectrum,
! optional cloud/haze optics) and build the reference/Planck/band-
! optimizer state. One call per process.
!
  character(kind=c_char), intent(in) :: data_root(*)      ! null-terminated; '' = compiled default
  character(kind=c_char), intent(in) :: solar_file_in(*)  ! null-terminated; '' = compiled default
  real(c_double),  value, intent(in) :: scon_in           ! stellar constant / 2 [W m-2]; <=0 = default
  real(c_double),  value, intent(in) :: g_in              ! surface gravity [m s-2]; <=0 = default
  integer(c_int),  value, intent(in) :: do_clouds         ! nonzero enables the cloud RT path
  integer(c_int),  value, intent(in) :: do_haze           ! nonzero enables the CARMA haze RT path
  integer(c_int) :: ierr

  character(len=256) :: root_f, solar_f

  if (initialized .or. finalized) then
    ierr = EXORT_ALREADY_INIT
    return
  endif

  call c_to_f_string(data_root, root_f)
  call c_to_f_string(solar_file_in, solar_f)

  if (len_trim(root_f) > 0) then
    ! table readers concatenate exort_rootdir//subdir; require trailing slash
    if (root_f(len_trim(root_f):len_trim(root_f)) /= '/') then
      root_f = trim(root_f)//'/'
    endif
    exort_rootdir = root_f
  endif
  if (len_trim(solar_f) > 0) solar_file = solar_f
  if (scon_in > 0.0_r8) shr_const_scon = scon_in
  if (g_in    > 0.0_r8) exo_g = g_in
  do_exo_clouds = (do_clouds /= 0)
  do_exo_haze   = (do_haze /= 0)

  ! propagate runtime config exactly as read_namelist does for the executable
  scon        = shr_const_scon
  SHR_CONST_G = exo_g

  write(*,*) '=== exort_config (library) ================='
  write(*,*) 'data_root     = ', trim(exort_rootdir)
  write(*,*) 'solar_file    = ', trim(solar_file)
  write(*,*) 'shr_const_scon= ', shr_const_scon, ' W m-2'
  write(*,*) 'exo_g         = ', exo_g, ' m s-2'
  write(*,*) 'do_exo_clouds = ', do_exo_clouds
  write(*,*) 'do_exo_haze   = ', do_exo_haze
  write(*,*) '============================================'

  ! same initialization chain as main.F90
  call initialize_kcoeff
  call initialize_solar
  if (do_exo_clouds) call initialize_cldopts
  if (do_exo_haze) call initialize_hazeopts
  call init_ref
  call init_model_specific
  call init_planck
  call initialize_radbuffer

  initialized = .true.
  ierr = EXORT_OK

end function exort_init

!============================================================================

function exort_run_column(state, result_out) &
         bind(c, name='exort_run_column') result(ierr)
  type(column_state_t),  intent(in)  :: state
  type(column_result_t), intent(out) :: result_out
  integer(c_int) :: ierr

  if (.not. initialized .or. finalized) then
    ierr = EXORT_NOT_INITIALIZED
    return
  endif

  call run_one_column(state, result_out)
  ierr = EXORT_OK

end function exort_run_column

!============================================================================

function exort_run_columns(ncols, states, results) &
         bind(c, name='exort_run_columns') result(ierr)
  integer(c_int), value, intent(in)  :: ncols
  type(column_state_t),  intent(in)  :: states(ncols)
  type(column_result_t), intent(out) :: results(ncols)
  integer(c_int) :: ierr

  integer :: i

  if (.not. initialized .or. finalized) then
    ierr = EXORT_NOT_INITIALIZED
    return
  endif
  if (ncols < 1) then
    ierr = EXORT_BAD_ARGS
    return
  endif

  do i = 1, ncols
    call run_one_column(states(i), results(i))
  end do
  ierr = EXORT_OK

end function exort_run_columns

!============================================================================

function exort_finalize() bind(c, name='exort_finalize') result(ierr)
!
! Bookkeeping close-out: further run/init calls are refused. Tables are
! static module data, so nothing is deallocated; re-initialization in
! the same process is not supported in Stage D.
!
  integer(c_int) :: ierr

  finalized = .true.
  ierr = EXORT_OK

end function exort_finalize

!============================================================================

subroutine run_one_column(state, res)
!
! Mirror of the per-column sequence in main.F90: derive dry pressures,
! set the (unused-in-1D) terrain/orbit inputs, run aerad_driver.
!
  type(column_state_t),  intent(in)  :: state
  type(column_result_t), intent(out) :: res

  integer, parameter :: nazm_tshadow = 1
  real(r8) :: msdist, rtgt, solar_azm_ang, tazm_ang, tslope_ang
  integer  :: tslas_tog, tshadow_tog
  real(r8), dimension(nazm_tshadow) :: cosz_horizon, TCx_obstruct, TCz_obstruct
  real(r8), dimension(pver)  :: pdeldry
  real(r8), dimension(pverp) :: pintdry

  call physconst_setgas(state%mwdry, state%cpdry)

  msdist = 1.0
  rtgt = 1.0
  solar_azm_ang = 0.0
  tazm_ang = 0.0
  tslope_ang = 0.0
  tslas_tog = 0
  tshadow_tog = 1
  cosz_horizon(:) = 0.0
  TCx_obstruct(:) = 0.0
  TCz_obstruct(:) = 0.0

  ! define dry as wet*(1-H2OMMR); interfaces map mid-layer H2OMMR the
  ! CESM way (top interface takes the top layer, each lower interface
  ! takes the mid-layer below it)
  pdeldry(:) = state%pdel(:)*(1.-state%h2ommr(:))
  pintdry(1)       = state%pint(1)*(1.-state%h2ommr(1))
  pintdry(2:pverp) = state%pint(2:pverp)*(1.-state%h2ommr(:))

  call aerad_driver(state%h2ommr, state%co2mmr, &
                    state%ch4mmr, state%c2h6mmr, &
                    state%nh3mmr, state%commr, &
                    state%h2mmr,  state%n2mmr, state%o3mmr, state%o2mmr, &
                    state%cicewp, state%cliqwp, state%cfrc, &
                    state%rei, state%rel, &
                    state%cicewp_co2, state%rei_co2, &
                    state%carmammr, &
                    state%ts, state%ps, state%pmid, &
                    state%pdel, pdeldry, state%tmid, state%pint, pintdry, &
                    state%coszrs, msdist, &
                    state%asdir, state%aldir, &
                    state%asdif, state%aldif, &
                    state%srf_emiss, &
                    rtgt, solar_azm_ang, tazm_ang, tslope_ang, &
                    tslas_tog, tshadow_tog, nazm_tshadow, cosz_horizon, &
                    TCx_obstruct, TCz_obstruct, state%zint, &
                    res%sw_dtdt, res%lw_dtdt, &
                    res%lw_dnflux, res%lw_upflux, res%sw_upflux, res%sw_dnflux, &
                    res%lw_dnflux_spectral, res%lw_upflux_spectral, &
                    res%sw_upflux_spectral, res%sw_dnflux_spectral, &
                    res%vis_dir, res%vis_dif, res%nir_dir, res%nir_dif, &
                    res%sol_toa )

end subroutine run_one_column

!============================================================================

subroutine c_to_f_string(cstr, fstr)
  character(kind=c_char), intent(in)  :: cstr(*)
  character(len=*),       intent(out) :: fstr

  integer :: i

  fstr = ' '
  do i = 1, len(fstr)
    if (cstr(i) == c_null_char) exit
    fstr(i:i) = cstr(i)
  end do

end subroutine c_to_f_string

end module exort_lib_mod
