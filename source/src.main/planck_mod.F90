
module planck_mod

use shr_kind_mod,      only: r8 => shr_kind_r8

implicit none
save

public :: PLANCKf

contains

!============================================================================      

  function PLANCKf(e,t1)

!------------------------------------------------------------------------          
!                                                                                  
! Purpose: Computes integral of the Planck function between zero and a             
!          given wavelength                                                        
!                                                                                  
!------------------------------------------------------------------------          
!     ******************************************************                       
!     *  Purpose             :  Calculate Planck Function  *                       
!     *  Subroutines Called  :  None                       *                       
!     *  Input               :  WAVE, TEMP                 *                       
!     *  Output              :  PLANK                      *                       
!     * ****************************************************           
!                                                                                  
!  THIS SUBROUTINE COMPUTES THE INTEGRAL OF THE PLANCK FUNCTION BETWEEN            
!  ZERO AND THE SPECIFIED VALUE OF LAMBDA.  THUS (USING XL AS LAMBDA)              
!  WE WANT TO INTEGRATE                                                            
!  R = INTEGRAL(XL=0 TO XL=XLSPEC) ( C1*XL**-5* / (EXP(C2/XL*T)-1) )*              
!  SUBSTITUTING U=C2/(XL*T), THE INTEGRAL BECOMES                                  
!  R = A CONSTANT TIMES INTEGRAL (USPEC TO INFINITY) OF                            
!            ( U**3 / (EXP(U) - 1) )*DU                                            
!  THE APPROXIMATIONS SHOWN HERE ARE ON PAGE 998 OF ABRAMOWITZ AND                 
!  UNDER THE HEADING OF DEBYE FUNCTIONS.  C2 IS THE PRODUCT OF PLANCK'S            
!  CONSTANT AND THE SPEED OF LIGHT DIVIDED BY BOLTZMANN'S CONSTANT.                
!  C2 = 14390 WHEN LAMBDA IS IN MICRONS.                                           
!  THE FACTOR 0.15399 IS THE RECIPROCAL OF SIX TIMES                               
!  THE SUM OF (1/N**2) FOR ALL N FROM ONE TO INFINITY.  IT IS CHOSEN TO            
!  NORMALIZE THE INTEGRAL TO A MAXIMUM VALUE OF UNITY.                             
!  RADIATION IN REAL UNITS IS OBTAINED BY MULTIPLYING THE INTEGRAL BY              
!  THE STEFAN-BOLTZMANN CONSTANT TIMES T**4.                                       

    implicit none
!------------------------------------------------------------------------          
!                                                                                  
! Input Arguments                                                                  
!   
                                                                                   
    real(r8), intent(in)  :: e                                                     
    real(r8), intent(in)  :: t1                                                    
                                                                                   
                                                                                   
!------------------------------------------------------------------------          
!                                                                                  
! Local Variables                                                                  
!                                                                                  

    real(r8) :: PLANCKf
    real(r8), dimension(5) :: am
    real(r8) :: d
    real(r8) :: v1
    real(r8) :: a
    integer :: m

!------------------------------------------------------------------------          
!                                                                                  
! Start Code                                                                       
!                                             
                                                                                  
    d = 0.d0                                                                       
    v1 = e/t1                                                                      
                                                                                   
    if(v1 <= 1.d0) then                                                            
      d = 1.0-0.15399*V1**3 *  &                                                   
      (1./3.-v1/8.+v1**2/60.-v1**4/5040.+v1**6/272160.-V1**8/13305600.)            
    endif                                                                          
                                                                                   
    if(v1 > 1.d0 .and. v1 <= 50.d0) then                                           
      do m=1,5                                                                     
        a = dble(m)*v1                                                             
        am(m) = 0.15399*exp(-a)/m**4*(((a+3.)*a+6.)*a+6.)                          
      enddo                                                                        
                                                                                   
      d = am(1)+am(2)+am(3)+am(4)+am(5)                                            
    endif                                                                          
                                                                                   
    PLANCKf = d*t1**4                                                              
                                                                                   
    return                                                                         
                                                                                   
  end function PLANCKf                                                                                                    

end module planck_mod
