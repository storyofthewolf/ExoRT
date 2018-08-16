
module cloud

implicit none
public

  ! Cloud mie data
  character(len=256), parameter :: cldoptsL_file = '/projects/wolfet/models/ExoRT/data/cloud/cloudoptics_liquid_mie_n68.nc'
  character(len=256), parameter :: cldoptsI_file = '/projects/wolfet/models/ExoRT/data/cloud/cloudoptics_ice_mie_n68.nc'

end module cloud
