
module cloud

implicit none
public

  ! directory
  character(len=256), parameter :: dircld = 'data/cloud/'

  ! Cloud mie data
  character(len=256), parameter :: cldoptsL_file = 'cloudoptics_h2o_liquid_mie_n84.nc'
  character(len=256), parameter :: cldoptsI_file = 'cloudoptics_h2o_ice_mie_n84.nc'
  character(len=256), parameter :: cldoptsICO2_file = 'cloudoptics_co2_ice_mie_n84.nc'

  ! CARMA haze aerosol optics (fractal aggregate particles; a Mie-sphere
  ! table haze_n84_b40_mie.nc is also available)
  character(len=256), parameter :: diraer = 'data/aerosol/'
  character(len=256), parameter :: hazeopts_file = 'haze_n84_b40_fractal_interp.nc'

end module cloud
