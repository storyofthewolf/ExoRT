#include <params.h>
#include <misc.h>

module solar

implicit none
public

  ! solar spectral file names
  !character(len=256), parameter :: solar_file = '/Users/wolfe/Models/RT/RT_offline/solar_data/G2V_SUN_n35.nc'
  character(len=256), parameter :: solar_file = '/Users/wolfe/Models/RT/RT_offline/solar_data/blackbody_3400K_n35.nc'

end module solar