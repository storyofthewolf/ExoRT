
module exort_column_run
!----------------------------------------------------------------------
! Single-column solve on the column_state_t / column_result_t structs
! (Stage E1). run_one_column is the ONE place the per-column call into
! aerad_driver is assembled; both the standalone executable (main.F90
! column loop) and the C library (exort_lib_mod) use it, so the two
! entry points cannot drift.
!
! Thread-safety: run_one_column writes no module-scope state (see
! STAGE_E_AUDIT.md) — the Stage E2 OpenMP loop parallelizes over calls
! to it.
!----------------------------------------------------------------------

use shr_kind_mod,      only: r8 => shr_kind_r8
use ppgrid,            only: pver, pverp
use exo_radiation_mod, only: aerad_driver
use exort_column_mod,  only: column_state_t, column_result_t

implicit none
private

public :: run_one_column

contains

!============================================================================

subroutine run_one_column(state, res)
!
! Per-column sequence: derive dry pressures, set the (unused-in-1D)
! terrain/orbit inputs, run aerad_driver.
!
  type(column_state_t),  intent(in)  :: state
  type(column_result_t), intent(out) :: res

  integer, parameter :: nazm_tshadow = 1
  real(r8) :: msdist, rtgt, solar_azm_ang, tazm_ang, tslope_ang
  integer  :: tslas_tog, tshadow_tog
  real(r8), dimension(nazm_tshadow) :: cosz_horizon, TCx_obstruct, TCz_obstruct
  real(r8), dimension(pver)  :: pdeldry
  real(r8), dimension(pverp) :: pintdry

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

  ! Optional condensed-phase/surface inputs are passed by keyword (the
  ! state struct always carries the fields, zero-filled by callers that
  ! don't use them, so passing them unconditionally preserves behavior).
  call aerad_driver(state%h2ommr, state%co2mmr, &
                    state%ch4mmr, state%c2h6mmr, &
                    state%nh3mmr, state%commr, &
                    state%h2mmr,  state%n2mmr, state%o3mmr, state%o2mmr, &
                    state%ts, state%ps, state%pmid, &
                    state%pdel, pdeldry, state%tmid, state%pint, pintdry, &
                    state%coszrs, msdist, &
                    state%asdir, state%aldir, &
                    state%asdif, state%aldif, &
                    rtgt, solar_azm_ang, tazm_ang, tslope_ang, &
                    tslas_tog, tshadow_tog, nazm_tshadow, cosz_horizon, &
                    TCx_obstruct, TCz_obstruct, state%zint, &
                    res%sw_dtdt, res%lw_dtdt, &
                    res%lw_dnflux, res%lw_upflux, res%sw_upflux, res%sw_dnflux, &
                    res%lw_dnflux_spectral, res%lw_upflux_spectral, &
                    res%sw_upflux_spectral, res%sw_dnflux_spectral, &
                    res%vis_dir, res%vis_dif, res%nir_dir, res%nir_dif, &
                    res%sol_toa, &
                    ext_cicewp=state%cicewp, ext_cliqwp=state%cliqwp, &
                    ext_cfrc=state%cfrc, ext_rei=state%rei, ext_rel=state%rel, &
                    ext_cicewp_co2=state%cicewp_co2, ext_rei_co2=state%rei_co2, &
                    ext_carmammr=state%carmammr, &
                    ext_srf_emiss=state%srf_emiss, &
                    ext_mwdry=state%mwdry, ext_cpdry=state%cpdry )

end subroutine run_one_column

end module exort_column_run
