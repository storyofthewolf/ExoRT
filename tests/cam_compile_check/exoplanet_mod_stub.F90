module exoplanet_mod
! Harness stand-in for ExoCAM's SourceMods/src.share/exoplanet_mod.F90.
! Doubles as the TEMPLATE for the ExoCAM-side update required by
! 3dmodels/src.cam.exort: the three "v2 additions" below must be added to
! any ExoCAM config's exoplanet_mod that builds against src.cam.exort.
  use shr_kind_mod, only: r8 => shr_kind_r8
  implicit none
  public

  ! (harness-only: src.misc/ppgrid derives pver from this; real CESM does not)
  integer, public, parameter :: exo_pver = 40
  real(r8), public, parameter :: exo_pstd = 101325.0_r8

  ! --- existing ExoCAM symbols consumed by the bundle ---
  integer,  public, parameter :: exo_rad_step = 3
  logical,  public, parameter :: do_exo_rt_clearsky = .true.
  logical,  public, parameter :: do_exo_rt_spectral = .true.
  logical,  public, parameter :: do_exo_rt_optimize_bands = .true.
  real(r8), public, parameter :: Tmax = 400.
  real(r8), public, parameter :: swFluxLimit = 0.999
  real(r8), public, parameter :: lwFluxLimit = 0.999
  logical,  public, parameter :: do_carma_exort = .true.
  character(len=256), public, parameter :: exo_solar_file = &
      '/discover/nobackup/etwolf/models/ExoRT/data/stellar/G2V_SUN_n84.nc'
  real(r8), public, parameter :: exo_n2mmr   = 0.7_r8
  real(r8), public, parameter :: exo_h2mmr   = 0.0_r8
  real(r8), public, parameter :: exo_co2mmr  = 6.0e-4_r8
  real(r8), public, parameter :: exo_ch4mmr  = 0.0_r8
  real(r8), public, parameter :: exo_c2h6mmr = 0.0_r8
  real(r8), public, parameter :: exo_nh3mmr  = 0.0_r8
  real(r8), public, parameter :: exo_commr   = 0.0_r8

  ! --- v2 additions REQUIRED by src.cam.exort ---
  ! The shared driver (exo_radiation_mod.F90) gates its cloud and haze
  ! kernels on these; they must exist in every exoplanet_mod that links
  ! the v2 bundle.
  logical,  public, parameter :: do_exo_clouds = .true.            ! 3-D always wants H2O clouds
  logical,  public, parameter :: do_exo_haze   = do_carma_exort    ! haze RT iff CARMA runs

  ! --- co2condense-config-only (required when built with -DEXORT_CO2CLD) ---
  logical,  public, parameter :: do_exo_condense_co2 = .true.

end module exoplanet_mod
