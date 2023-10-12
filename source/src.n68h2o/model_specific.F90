
module exo_model_specific

!---------------------------------------------------------------------       
! Purpose:                                                                   

  use shr_kind_mod,     only: r8 => shr_kind_r8
  use radgrid
  use spmd_utils,       only: masterproc 

  implicit none
  private
  save
  
  public :: init_model_specific


!============================================================================
contains
!============================================================================


!============================================================================
  subroutine init_model_specific
!------------------------------------------------------------------------

    !! nothing to do

  end subroutine init_model_specific

end module exo_init_ref
