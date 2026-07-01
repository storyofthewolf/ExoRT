
module exoplanet_mod
  ! This is paired down version of exoplanet_mod.F90, the primary version of which is
  ! the core module for changing geophysical and stellar parameters in ExoCAM.
  ! This module has been simplified to remove items not explicitly needed in the
  ! standalone radiation.

  ! set basic atmosphere assumptions here.

  use shr_kind_mod,  only: r8 => shr_kind_r8
  implicit none
  public

  ! stellar spectrum file directory (v2: data/solar/ renamed to data/stellar/)
  character(len=256), parameter :: dirsol = 'data/stellar/'

  ! Cloud radiative transfer (Stage C).
  ! When .false. (default) the cloud optical-depth path is skipped entirely and
  ! the model is bit-for-bit identical to the cloud-free build. When .true. the
  ! cloud optics tables are loaded at init and condensate fields from the input
  ! file (cicewp/cliqwp/rei/rel for H2O; cicewp_co2/rei_co2 for CO2 ice) drive
  ! cloud opacity. Runtime-overridable via the user_nl_exort namelist so one
  ! build can iterate clear/cloudy decks without recompiling.
  logical :: do_exo_clouds = .false.

  ! CARMA haze radiative transfer (Stage C3).
  ! When .false. (default) the aerosol optical-depth path is skipped entirely
  ! and the model is bit-for-bit identical to the haze-free build. When .true.
  ! the haze optics table is loaded at init and the binwise aerosol mass mixing
  ! ratio carmammr(pver,nelem,nbin) from the input file drives haze opacity.
  ! Haze formation physics (CARMA microphysics) lives outside ExoRT; the 1-D
  ! path consumes carmammr as a given input array. Runtime-overridable via the
  ! user_nl_exort namelist.
  logical :: do_exo_haze = .false.

  ! Radiation Spectral Band Optimization
  logical, parameter :: do_exo_rt_optimize_bands = .true.
  real(r8), parameter :: Tmax = 400.          !! Maximum expected temperature for thermal band optimization
  real(r8), parameter :: swFluxLimit = 0.999  !! Fraction of stellar flux captured in bands, rescaled
  real(r8), parameter :: lwFluxLimit = 0.999  !! Fraction of thermal flux captured in bands, not rescaled

  ! solar spectral file names
  ! Default: runtime-overridable via user_nl_exort namelist.
  !character(len=256) :: solar_file = 'WD_5000K_n84.nc'
  !character(len=256) :: solar_file = 'trappist1_lincowski2018_n28.nc'
  character(len=256) :: solar_file = 'G2V_SUN_n84.nc'   ! v2 default (n84 grid); override via user_nl_exort
  !character(len=256) :: solar_file = 'blackbody_3400K_n28.nc'
  !character(len=256) :: solar_file = 'G2V_SUN_n68.nc'
  !character(len=256) :: solar_file = 'LHS1140_spectra_n42.nc'
  !character(len=256) :: solar_file = 'bt-settl_2600_logg4.5_FeH0_n68.nc'

  ! solar constant, 1D kludge
  ! in 1D code, we use the true solar constant divided by 2
  ! a second factor of 2 comes with a 0.5 zenith angle.
  ! Thus 1360 Wm-2 --> 680 Wm-2 --> 340 Wm-2 in the standard SW profile
  ! Default: runtime-overridable via user_nl_exort namelist.
  real(r8) :: shr_const_scon = 680.0_r8  !Wm-2, Current Earth (approximate)
  !real(r8) :: shr_const_scon = 451.166_r8 !Wm-2, Early Mars

  ! number of levels
  !integer, parameter :: exo_pver   = 40  !CESM
  !integer, parameter :: exo_pver   = 45  !CESM
  !integer, parameter :: exo_pver   = 66  !WACCM
  !integer, parameter :: exo_pver   = 49  !US1976
  !integer, parameter :: exo_pver   = 69   !2 bar CO2
  !integer, parameter :: exo_pver   = 51   !
  integer, parameter :: exo_pver   = 300 !ExoMIP

  ! must set gravity!!!!
  ! Default: runtime-overridable via user_nl_exort namelist.
  !real(r8) :: exo_g = 0.93_r8*9.80616_r8 !Earth
  !real(r8) :: exo_g = 7.22925_r8 !Trappist-1e
  real(r8) :: exo_g = 3.711_r8 ! Mars

  ! Reference pressures (should not effect final answer)
  real(r8), parameter :: exo_pstd = 100000.  !Pascals

end module exoplanet_mod
