#include <params.h>
#include <misc.h>

module kabs
! version.co2_h2ovar

  implicit none
  public

  ! K coefficient file names
  character(len=256), parameter :: k_file = '/Users/wolfe/Models/RT/RT_offline/absorption_data/co2_h2ovar/CO2_H2Ovar_MERGED.nc'

  ! H2O Self and Foreign continuum data
  character(len=256), parameter :: kh2ocont_file = '/Users/wolfe/Models/RT/RT_offline/absorption_data/continuum/h2o_cont_n63.nc'

end module kabs