
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

    call setup_major_gas_matrix
    call setup_gray_gas_matrix

 end subroutine init_model_specific


!============================================================================
  subroutine setup_major_gas_matrix
!------------------------------------------------------------------------
! PURPOSE:  Load full k-distribution from files into major gas absorber
!           matrix.
!------------------------------------------------------------------------ 

    use kabs
    implicit none

!------------------------------------------------------------------------ 
!
! Local Variables
!
  integer :: iw
!------------------------------------------------------------------------ 
! Start Code
!

    ! Initialize all bins to zero
    k_major_data(:,:,:,:,:) = 0.0
    k_major_data(iH2O,:,:,:,:) = k_h2o(:,:,:,:) 
    k_major_data(iCO2,:,:,:,:) = k_co2(:,:,:,:)
    k_major_data(iCH4,:,:,:,:) = k_ch4(:,:,:,:)
    k_major_data(iC2H6,:,:,:,:) = k_ch4(:,:,:,:)

  
  end subroutine setup_major_gas_matrix


!============================================================================ 
  subroutine setup_gray_gas_matrix
!------------------------------------------------------------------------
! PURPOSE:  Calculate gray gas coefficients to use for minor species.  
!           Gray gases are a gauss-point weighted averaged in each band.
!------------------------------------------------------------------------ 

    use kabs
    use shr_const_mod, only: SHR_CONST_G, SHR_CONST_PSTD, SHR_CONST_AVOGAD
    use rad_interp_mod, only:  bilinear_interpK_grey

    implicit none

!------------------------------------------------------------------------ 
!
! Local Variables
!
  integer :: iq, ig, iw


!------------------------------------------------------------------------ 
! Start Code
!

   ! Create grey gas absorption matrix
   iq=0
   k_grey_data(:,:,:,:) = 0.0
   do iw=1, ntot_wavlnrng
     do ig=1, ngauss_pts(1)
       iq = iq + 1
       k_grey_data(iH2O,iw,:,:)  = k_grey_data(iH2O,iw,:,:)  + k_h2o(iw,ig,:,:)  * g_weight(iq) 
       k_grey_data(iCO2,iw,:,:)  = k_grey_data(iCO2,iw,:,:)  + k_co2(iw,ig,:,:)  * g_weight(iq)
       k_grey_data(iCH4,iw,:,:)  = k_grey_data(iCH4,iw,:,:)  + k_ch4(iw,ig,:,:)  * g_weight(iq)
       k_grey_data(iC2H6,iw,:,:) = k_grey_data(iC2H6,iw,:,:) + k_c2h6(iw,ig,:,:) * g_weight(iq)
     enddo
   enddo

  end subroutine setup_gray_gas_matrix

end module exo_model_specific
