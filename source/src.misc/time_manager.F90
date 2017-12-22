#include <params.h>
#include <misc.h>

module time_manager

implicit none
public

contains

integer function get_nstep()
  !Kludge, for cloud random number generator hook-up
  get_nstep=9404
end function get_nstep

end module time_manager
