module physconst
! Harness stand-in for CAM/ExoCAM's physconst (ExoCAM: configs/*/SourceMods/
! src.cam/physconst.F90). It exports ONLY the names src.cam.exort imports,
! each of which exists in mainline ExoCAM's physconst (cam_aqua_fv,
! cam_aqua_se, cam_land_fv, cam_mixed_fv). ExoRT's own 1-D shim
! (source/src.misc/physconst.F90) is not used here: it reads ExoRT-only
! shr_const_mod constants that ExoCAM's shr_const_mod does not have.
! Note: ExoCAM's circumbinary config lacks mwc2h6, mwco and mwnh3.
  use shr_kind_mod, only: r8 => shr_kind_r8
  implicit none
  private
  save
  real(r8), public, parameter :: scon   = 1.361e3_r8
  real(r8), public, parameter :: cpair  = 1.00464e3_r8
  real(r8), public, parameter :: cappa  = 0.2857_r8
  real(r8), public, parameter :: mwdry  = 28.966_r8
  real(r8), public, parameter :: mwn2   = 28.0_r8
  real(r8), public, parameter :: mwco2  = 44.0_r8
  real(r8), public, parameter :: mwch4  = 16.0_r8
  real(r8), public, parameter :: mwc2h6 = 30.0_r8
  real(r8), public, parameter :: mwnh3  = 17.0_r8
  real(r8), public, parameter :: mwco   = 28.0_r8
  real(r8), public, parameter :: mwh2o  = 18.0_r8
  real(r8), public, parameter :: mwo2   = 32.0_r8
  real(r8), public, parameter :: mwh2   = 2.0_r8
  real(r8), public, parameter :: mwo3   = 48.0_r8
end module physconst
