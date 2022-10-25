
module cloud

implicit none
public

  ! directory
  character(len=256), parameter :: dircld = 'data/cloud/'

  ! Cloud mie data
  character(len=256), parameter :: cldoptsL_file = 'cloudoptics_h2o_liquid_mie_n68.nc'
  character(len=256), parameter :: cldoptsI_file = 'cloudoptics_h2o_ice_mie_n68.nc'
!  character(len=256), parameter :: cldoptsICO2_file = 'cloudoptics_co2_ice_mie_n68.nc'
  character(len=256), parameter :: cldoptsICO2_file = 'cloudoptics_co2_ice_mie_n68_r200.nc'

end module cloud
