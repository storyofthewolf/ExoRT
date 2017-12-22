#include <misc.h>
#include <params.h>

module pmgrid

!----------------------------------------------------------------------- 
! 
! Purpose:  Stub for 1D RT hook-up
! 
!-----------------------------------------------------------------------

   use ppgrid

   integer, parameter :: plev   = pver       ! number of vertical levels
   integer, parameter :: plevp  = plev+1     ! plev + 1

end module pmgrid
