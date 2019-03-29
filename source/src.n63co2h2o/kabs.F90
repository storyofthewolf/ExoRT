#include <params.h>
#include <misc.h>

module kabs
! version.n63co2h2o

  implicit none
  public

  ! K coefficient file names
  character(len=256), parameter :: dirk = '/projects/wolfet/models/ExoRT/data/kdist/n63co2h2o/'
  character(len=256), parameter :: k_file = 'n63_CO2_H2O_MERGED.nc'

  ! H2O Self and Foreign continuum data
  character(len=256), parameter :: dirct = '/projects/wolfet/models/ExoRT/data/continuum/'
  character(len=256), parameter :: kh2ocont_file = 'h2o_cont_n63.nc'

end module kabs
