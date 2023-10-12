
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

    call map_co2cont_gpt

 end subroutine init_model_specific


!============================================================================
  subroutine map_co2cont_gpt
!------------------------------------------------------------------------
! PURPOSE:  CO2 Continuum data sets are on 8 point gauss interval bin.  Adjust
!           to match number of gauss intervals used for major absorbing gases 
!           Required to populated the CO2 continuum from file!!
!------------------------------------------------------------------------

    implicit none

!------------------------------------------------------------------------
!                                                                                                               
! Local Variables                                                                                               
!                                                                                                               
  integer :: iw
  integer :: ig
  integer :: itc

!------------------------------------------------------------------------                                       
! Start Code                                                                                                    
!                                                                                                               

    ! Initialize reduced continuum k coefficient arrays
    kco2cont(:) = 0.0
    !
    ! Gauss point adjustment for CO2 continuum
    !
    itc = 0
    do iw=1, ntot_wavlnrng
      if (ngauss_pts(iw) .eq. 8) then ! no adjustment needed
        do ig=1, ngauss_pts(iw)
          itc = itc + 1
          kco2cont(itc) = kco2cont_8gpt(ig,iw)
        enddo
      endif
      if (ngauss_pts(iw) .eq. 16) then
        do ig=1, ngauss_pts(iw)
          itc = itc + 1
          kco2cont(itc) = kco2cont_8gpt(map8to16gpt(ig),iw)
        enddo
      endif
    enddo

  end subroutine map_co2cont_gpt

end module exo_model_specific
