#include <params.h>
#include <misc.h>

module kabs
! version.co2_h2ovar

  implicit none
  public

  ! K coefficient file names
  character(len=256), parameter :: dirk = '/projects/wolfet/models/ExoRT/data/kdist/n63co2h2ovar/'
  character(len=256), parameter :: k_file = 'CO2_H2Ovar_MERGED.nc'

  ! H2O Self and Foreign continuum data
  character(len=256), parameter :: dirct = '/projects/wolfet/models/ExoRT/data/continuum/'
  character(len=256), parameter :: kh2ocont_file = 'h2o_cont_n63.nc'

end module kabs
