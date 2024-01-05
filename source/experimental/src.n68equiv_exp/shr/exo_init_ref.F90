
module exo_init_ref

!---------------------------------------------------------------------       
! Purpose:                                                                   

  use shr_kind_mod,     only: r8 => shr_kind_r8
  use shr_const_mod,    only: SHR_CONST_PI
  use physconst,        only: scon, mwdry
  use radgrid
  use spmd_utils,       only: masterproc 

  implicit none
  public  

  ! Approximate smallest double precision floating point difference 
  !real(r8), parameter :: SMALLd = 1.0d-12
  !real(r8), parameter :: SMALLe = 1.0e-12
  real(r8), parameter :: SMALLd = 1.0d-8
  real(r8), parameter :: SMALLe = 1.0e-8
  !real(r8), parameter :: SMALLd = 1.0d-6
  !real(r8), parameter :: SMALLe = 1.0e-6

  real(r8), parameter :: sqrt3 = 1.732050808d0      ! square root of 3
  real(r8), parameter :: mb_to_atm = 9.869233e-4    ! convert pressure from Pa to atm

  !------------------------------------------------------------------------------  
  ! Radiative transfer model variable/array declarations                           
  !
  ! Assign beginning and end wavelength range and point indices for each         
  !  wavelength group                                               

  ! full spectral range
  !integer, parameter  :: lw_iwbeg = 1     ! thermal band wvl integration limits                  
  !integer, parameter  :: lw_iwend = ntot_wavlnrng
  !integer, parameter  :: sw_iwbeg = 1     ! solar band wvl integration limits                    
  !integer, parameter  :: sw_iwend = ntot_wavlnrng
  !integer, parameter  :: lw_ipbeg = 1     ! thermal band gpt integration limits                  
  !integer, parameter  :: lw_ipend = ntot_gpt
  !integer, parameter  :: sw_ipbeg = 1     ! solar band gpt integration limits                    
  !integer, parameter  :: sw_ipend = ntot_gpt

  ! reduced integration limits for 3dmodel efficiency
  ! modifiable to suit stellar spectra and planet emission temperatures
  integer, parameter  :: lw_iwbeg = 1     ! thermal band wvl integration limits                  
  integer, parameter  :: lw_iwend = 35 
  integer, parameter  :: sw_iwbeg = 16    ! solar band wvl integration limits                    
  integer, parameter  :: sw_iwend = 68
  integer, parameter  :: lw_ipbeg = 1     ! thermal band gpt integration limits                  
  integer, parameter  :: lw_ipend = 280
  integer, parameter  :: sw_ipbeg = 121   ! solar band gpt integration limits                    
  integer, parameter  :: sw_ipend = 544

  !                                                                              
  ! set two-stream model coefficients                                            
  !                                           
  ! for solar stream (quadrature)    
  real(r8), parameter :: U1Isol  = sqrt3                             ! 2*PI / mu1 factors 
  real(r8), parameter :: U1I2sol = 0.5d0*sqrt3                       ! 2*PI / mu1 factors 
  real(r8), parameter :: U1Ssol  = 2.0*SHR_CONST_PI/U1Isol  ! mu1 factors
  ! for thermal stream (hemispsheric mean)      
  real(r8), parameter :: U1Iir   = 2.d0                              ! mu1 factors 
  real(r8), parameter :: U1I2ir  = 1.d0                              ! mu1 factors 
  real(r8), parameter :: U1Sir   = 2.0*SHR_CONST_PI/U1Iir            ! mu1 factors 

  real(r8), dimension(ntot_gpt) :: gw_solflux
  real(r8), dimension(ntot_wavlnrng) :: solflux

  ! Define constants for Gauss quadratrue calculations
  integer, parameter  :: ngangles = 3               ! # of Gauss angles to use

  ! Gauss weight vector and zenith angle
  real(r8), dimension(ntot_gpt) :: g_weight
  real(r8), dimension(ngangles) :: g_ang_weight
  real(r8), dimension(ngangles) :: g_angle

  ! Weights for Gaussian quadrature (for each probability interval) over a
  ! hemisphere subdivided into 'ngangles' equal angles  [none] (there are
  ! 'ngangles' of these):

  ! For ngangles = 3:
  data g_angle / 0.21234054, 0.59053314, 0.91141204 /
  data g_ang_weight / 0.06982698, 0.22924111, 0.20093191 /


!============================================================================
contains
!============================================================================


  subroutine init_ref
!------------------------------------------------------------------------      
! Purpose: Initial reference value arrays 
!------------------------------------------------------------------------ 
  implicit none
!------------------------------------------------------------------------ 
! Local Variables                                                         
!
integer :: iq
integer :: iw
integer :: ig
integer :: ip

!------------------------------------------------------------------------      
!                                                                              
! Start Code                                                                   
!                                                                              
    !write(*,*)   "INIT_REF: INITIALIZE GUASS POINT ARRAYS"
    ! Arrange g_weight(ntot_pt) array 
    iq = 0
    do iw=1,ntot_wavlnrng
      do ig=1, ngauss_pts(iw)
        iq = iq + 1
        if (ngauss_pts(iw) .eq. ngauss_8gpt) g_weight(iq) =  g_weight_8gpt(ig)
        if (ngauss_pts(iw) .eq. ngauss_16gpt) g_weight(iq) =  g_weight_16gpt(ig)
!!        if (iw .eq. 27) g_weight(iq) = solar_gweight_sp27(ig)   !!  decpreciated from old O2, O3 version         
!!        if (iw .eq. 28) g_weight(iq) = solar_gweight_sp28(ig)   !!  decpreciated from old O2, O3 version         
      enddo
    enddo

     ! Scale solar constant to namelist value
    solarflux(:) = solarflux(:)*scon/S0

    !
    ! Calculate the "average" wavenumber (1/wavelength) <wavenum()> for each    
    !  wavelength interval in each wavelength group [1/cm]:                     
    !                                                                           
    iq = 0
    ip = lw_ipbeg-1
    do iw=1,ntot_wavlnrng  ! "avg" wavenumber over each band                    
      do ig=1,ngauss_pts(iw)
        ip = ip+1
        ! Gauss-weighted solar flux in each probability interval:               
        gw_solflux(ip) = solarflux(iw)*g_weight(ip)
      enddo
    enddo

    if (masterproc) then
      write(*,*) "INIT_REF: total solar irradiance scaled to ",scon, "W m-2"
      write(*,*) "solar flux [W m-2] in each spectral interval"
      do iw=lw_iwbeg, sw_iwend
        write(*,*) iw, solarflux(iw)
      enddo
      write(*,*) "TOTAL SOLAR FLUX:", SUM(solarflux)
    endif

    call setup_major_gas_matrix
    call setup_gray_gas_matrix


 end subroutine init_ref


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

!------------------------------------------------------------------------ 
! Start Code
!


!=============================================================


   ! Create major gas absorption matrix
   !=====  interval 1:   =====!
   k01_major_data(iH2O,:,:,:) = k01_h2o(:,:,:) 
   k01_major_data(iCO2,:,:,:) = k01_co2(:,:,:)
   k01_major_data(iCH4,:,:,:) = k01_ch4(:,:,:)
   !=====  interval 2:   =====!
   k02_major_data(iH2O,:,:,:) = k02_h2o(:,:,:) 
   k02_major_data(iCO2,:,:,:) = k02_co2(:,:,:)
   k02_major_data(iCH4,:,:,:) = k02_ch4(:,:,:)
   !=====  interval 3:   =====!
   k03_major_data(iH2O,:,:,:) = k03_h2o(:,:,:) 
   k03_major_data(iCO2,:,:,:) = k03_co2(:,:,:)
   k03_major_data(iCH4,:,:,:) = k03_ch4(:,:,:)
   !=====  interval 4:   =====!
   k04_major_data(iH2O,:,:,:) = k04_h2o(:,:,:) 
   k04_major_data(iCO2,:,:,:) = k04_co2(:,:,:)
   k04_major_data(iCH4,:,:,:) = k04_ch4(:,:,:)
   !=====  interval 5:   =====!
   k05_major_data(iH2O,:,:,:) = k05_h2o(:,:,:) 
   k05_major_data(iCO2,:,:,:) = k05_co2(:,:,:)
   k05_major_data(iCH4,:,:,:) = k05_ch4(:,:,:)
   !=====  interval 6:   =====!
   k06_major_data(iH2O,:,:,:) = k06_h2o(:,:,:) 
   k06_major_data(iCO2,:,:,:) = k06_co2(:,:,:)
   k06_major_data(iCH4,:,:,:) = k06_ch4(:,:,:)
   !=====  interval 7:   =====!
   k07_major_data(iH2O,:,:,:) = k07_h2o(:,:,:) 
   k07_major_data(iCO2,:,:,:) = k07_co2(:,:,:)
   k07_major_data(iCH4,:,:,:) = k07_ch4(:,:,:)
   !=====  interval 8:   =====!
   k08_major_data(iH2O,:,:,:) = k08_h2o(:,:,:) 
   k08_major_data(iCO2,:,:,:) = k08_co2(:,:,:)
   k08_major_data(iCH4,:,:,:) = k08_ch4(:,:,:)
   !=====  interval 9:   =====!
   k09_major_data(iH2O,:,:,:) = k09_h2o(:,:,:) 
   k09_major_data(iCO2,:,:,:) = k09_co2(:,:,:)
   k09_major_data(iCH4,:,:,:) = k09_ch4(:,:,:)
   !=====  interval 10:   =====!
   k10_major_data(iH2O,:,:,:) = k10_h2o(:,:,:) 
   k10_major_data(iCO2,:,:,:) = k10_co2(:,:,:)
   k10_major_data(iCH4,:,:,:) = k10_ch4(:,:,:)
   !=====  interval 11:   =====!
   k11_major_data(iH2O,:,:,:) = k11_h2o(:,:,:) 
   k11_major_data(iCO2,:,:,:) = k11_co2(:,:,:)
   k11_major_data(iCH4,:,:,:) = k11_ch4(:,:,:)
   !=====  interval 12:   =====!
   k12_major_data(iH2O,:,:,:) = k12_h2o(:,:,:) 
   k12_major_data(iCO2,:,:,:) = k12_co2(:,:,:)
   k12_major_data(iCH4,:,:,:) = k12_ch4(:,:,:)
   !=====  interval 13:   =====!
   k13_major_data(iH2O,:,:,:) = k13_h2o(:,:,:) 
   k13_major_data(iCO2,:,:,:) = k13_co2(:,:,:)
   k13_major_data(iCH4,:,:,:) = k13_ch4(:,:,:)
   !=====  interval 14:   =====!
   k14_major_data(iH2O,:,:,:) = k14_h2o(:,:,:) 
   k14_major_data(iCO2,:,:,:) = k14_co2(:,:,:)
   k14_major_data(iCH4,:,:,:) = k14_ch4(:,:,:)
   !=====  interval 15:   =====!
   k15_major_data(iH2O,:,:,:) = k15_h2o(:,:,:) 
   k15_major_data(iCO2,:,:,:) = k15_co2(:,:,:)
   k15_major_data(iCH4,:,:,:) = k15_ch4(:,:,:)
   !=====  interval 16:   =====!
   k16_major_data(iH2O,:,:,:) = k16_h2o(:,:,:) 
   k16_major_data(iCO2,:,:,:) = k16_co2(:,:,:)
   k16_major_data(iCH4,:,:,:) = k16_ch4(:,:,:)
   !=====  interval 17:   =====!
   k17_major_data(iH2O,:,:,:) = k17_h2o(:,:,:) 
   k17_major_data(iCO2,:,:,:) = k17_co2(:,:,:)
   k17_major_data(iCH4,:,:,:) = k17_ch4(:,:,:)
   !=====  interval 18:   =====!
   k18_major_data(iH2O,:,:,:) = k18_h2o(:,:,:) 
   k18_major_data(iCO2,:,:,:) = k18_co2(:,:,:)
   k18_major_data(iCH4,:,:,:) = k18_ch4(:,:,:)
   !=====  interval 19:   =====!
   k19_major_data(iH2O,:,:,:) = k19_h2o(:,:,:) 
   k19_major_data(iCO2,:,:,:) = k19_co2(:,:,:)
   k19_major_data(iCH4,:,:,:) = k19_ch4(:,:,:)
   !=====  interval 20:   =====!
   k20_major_data(iH2O,:,:,:) = k20_h2o(:,:,:) 
   k20_major_data(iCO2,:,:,:) = k20_co2(:,:,:)
   k20_major_data(iCH4,:,:,:) = k20_ch4(:,:,:)
   !=====  interval 21:   =====!
   k21_major_data(iH2O,:,:,:) = k21_h2o(:,:,:) 
   k21_major_data(iCO2,:,:,:) = k21_co2(:,:,:)
   k21_major_data(iCH4,:,:,:) = k21_ch4(:,:,:)
   !=====  interval 22:   =====!
   k22_major_data(iH2O,:,:,:) = k22_h2o(:,:,:) 
   k22_major_data(iCO2,:,:,:) = k22_co2(:,:,:)
   k22_major_data(iCH4,:,:,:) = k22_ch4(:,:,:)
   !=====  interval 23:   =====!
   k23_major_data(iH2O,:,:,:) = k23_h2o(:,:,:) 
   k23_major_data(iCO2,:,:,:) = k23_co2(:,:,:)
   k23_major_data(iCH4,:,:,:) = k23_ch4(:,:,:)
   !=====  interval 24:   =====!
   k24_major_data(iH2O,:,:,:) = k24_h2o(:,:,:) 
   k24_major_data(iCO2,:,:,:) = k24_co2(:,:,:)
   k24_major_data(iCH4,:,:,:) = k24_ch4(:,:,:)
   !=====  interval 25:   =====!
   k25_major_data(iH2O,:,:,:) = k25_h2o(:,:,:) 
   k25_major_data(iCO2,:,:,:) = k25_co2(:,:,:)
   k25_major_data(iCH4,:,:,:) = k25_ch4(:,:,:)
   !=====  interval 26:   =====!
   k26_major_data(iH2O,:,:,:) = k26_h2o(:,:,:) 
   k26_major_data(iCO2,:,:,:) = k26_co2(:,:,:)
   k26_major_data(iCH4,:,:,:) = k26_ch4(:,:,:)
   !=====  interval 27:   =====!
   k27_major_data(iH2O,:,:,:) = k27_h2o(:,:,:) 
   k27_major_data(iCO2,:,:,:) = k27_co2(:,:,:)
   k27_major_data(iCH4,:,:,:) = k27_ch4(:,:,:)
   !=====  interval 28:   =====!
   k28_major_data(iH2O,:,:,:) = k28_h2o(:,:,:) 
   k28_major_data(iCO2,:,:,:) = k28_co2(:,:,:)
   k28_major_data(iCH4,:,:,:) = k28_ch4(:,:,:)
   !=====  interval 29:   =====!
   k29_major_data(iH2O,:,:,:) = k29_h2o(:,:,:) 
   k29_major_data(iCO2,:,:,:) = k29_co2(:,:,:)
   k29_major_data(iCH4,:,:,:) = k29_ch4(:,:,:)
   !=====  interval 30:   =====!
   k30_major_data(iH2O,:,:,:) = k30_h2o(:,:,:) 
   k30_major_data(iCO2,:,:,:) = k30_co2(:,:,:)
   k30_major_data(iCH4,:,:,:) = k30_ch4(:,:,:)
   !=====  interval 31:   =====!
   k31_major_data(iH2O,:,:,:) = k31_h2o(:,:,:) 
   k31_major_data(iCO2,:,:,:) = k31_co2(:,:,:)
   k31_major_data(iCH4,:,:,:) = k31_ch4(:,:,:)
   !=====  interval 32:   =====!
   k32_major_data(iH2O,:,:,:) = k32_h2o(:,:,:) 
   k32_major_data(iCO2,:,:,:) = k32_co2(:,:,:)
   k32_major_data(iCH4,:,:,:) = k32_ch4(:,:,:)
   !=====  interval 33:   =====!
   k33_major_data(iH2O,:,:,:) = k33_h2o(:,:,:) 
   k33_major_data(iCO2,:,:,:) = k33_co2(:,:,:)
   k33_major_data(iCH4,:,:,:) = k33_ch4(:,:,:)
   !=====  interval 34:   =====!
   k34_major_data(iH2O,:,:,:) = k34_h2o(:,:,:) 
   k34_major_data(iCO2,:,:,:) = k34_co2(:,:,:)
   k34_major_data(iCH4,:,:,:) = k34_ch4(:,:,:)
   !=====  interval 35:   =====!
   k35_major_data(iH2O,:,:,:) = k35_h2o(:,:,:) 
   k35_major_data(iCO2,:,:,:) = k35_co2(:,:,:)
   k35_major_data(iCH4,:,:,:) = k35_ch4(:,:,:)
   !=====  interval 36:   =====!
   k36_major_data(iH2O,:,:,:) = k36_h2o(:,:,:) 
   k36_major_data(iCO2,:,:,:) = k36_co2(:,:,:)
   k36_major_data(iCH4,:,:,:) = k36_ch4(:,:,:)
   !=====  interval 37:   =====!
   k37_major_data(iH2O,:,:,:) = k37_h2o(:,:,:) 
   k37_major_data(iCO2,:,:,:) = k37_co2(:,:,:)
   k37_major_data(iCH4,:,:,:) = k37_ch4(:,:,:)
   !=====  interval 38:   =====!
   k38_major_data(iH2O,:,:,:) = k38_h2o(:,:,:) 
   k38_major_data(iCO2,:,:,:) = k38_co2(:,:,:)
   k38_major_data(iCH4,:,:,:) = k38_ch4(:,:,:)
   !=====  interval 39:   =====!
   k39_major_data(iH2O,:,:,:) = k39_h2o(:,:,:) 
   k39_major_data(iCO2,:,:,:) = k39_co2(:,:,:)
   k39_major_data(iCH4,:,:,:) = k39_ch4(:,:,:)
   !=====  interval 40:   =====!
   k40_major_data(iH2O,:,:,:) = k40_h2o(:,:,:) 
   k40_major_data(iCO2,:,:,:) = k40_co2(:,:,:)
   k40_major_data(iCH4,:,:,:) = k40_ch4(:,:,:)
   !=====  interval 41:   =====!
   k41_major_data(iH2O,:,:,:) = k41_h2o(:,:,:) 
   k41_major_data(iCO2,:,:,:) = k41_co2(:,:,:)
   k41_major_data(iCH4,:,:,:) = k41_ch4(:,:,:)
   !=====  interval 42:   =====!
   k42_major_data(iH2O,:,:,:) = k42_h2o(:,:,:) 
   k42_major_data(iCO2,:,:,:) = k42_co2(:,:,:)
   k42_major_data(iCH4,:,:,:) = k42_ch4(:,:,:)
   !=====  interval 43:   =====!
   k43_major_data(iH2O,:,:,:) = k43_h2o(:,:,:) 
   k43_major_data(iCO2,:,:,:) = k43_co2(:,:,:)
   k43_major_data(iCH4,:,:,:) = k43_ch4(:,:,:)
   !=====  interval 44:   =====!
   k44_major_data(iH2O,:,:,:) = k44_h2o(:,:,:) 
   k44_major_data(iCO2,:,:,:) = k44_co2(:,:,:)
   k44_major_data(iCH4,:,:,:) = k44_ch4(:,:,:)
   !=====  interval 45:   =====!
   k45_major_data(iH2O,:,:,:) = k45_h2o(:,:,:) 
   k45_major_data(iCO2,:,:,:) = k45_co2(:,:,:)
   k45_major_data(iCH4,:,:,:) = k45_ch4(:,:,:)
   !=====  interval 46:   =====!
   k46_major_data(iH2O,:,:,:) = k46_h2o(:,:,:) 
   k46_major_data(iCO2,:,:,:) = k46_co2(:,:,:)
   k46_major_data(iCH4,:,:,:) = k46_ch4(:,:,:)
   !=====  interval 47:   =====!
   k47_major_data(iH2O,:,:,:) = k47_h2o(:,:,:) 
   k47_major_data(iCO2,:,:,:) = k47_co2(:,:,:)
   k47_major_data(iCH4,:,:,:) = k47_ch4(:,:,:)
   !=====  interval 48:   =====!
   k48_major_data(iH2O,:,:,:) = k48_h2o(:,:,:) 
   k48_major_data(iCO2,:,:,:) = k48_co2(:,:,:)
   k48_major_data(iCH4,:,:,:) = k48_ch4(:,:,:)
   !=====  interval 49:   =====!
   k49_major_data(iH2O,:,:,:) = k49_h2o(:,:,:) 
   k49_major_data(iCO2,:,:,:) = k49_co2(:,:,:)
   k49_major_data(iCH4,:,:,:) = k49_ch4(:,:,:)
   !=====  interval 50:   =====!
   k50_major_data(iH2O,:,:,:) = k50_h2o(:,:,:) 
   k50_major_data(iCO2,:,:,:) = k50_co2(:,:,:)
   k50_major_data(iCH4,:,:,:) = k50_ch4(:,:,:)
   !=====  interval 51:   =====!
   k51_major_data(iH2O,:,:,:) = k51_h2o(:,:,:) 
   k51_major_data(iCO2,:,:,:) = k51_co2(:,:,:)
   k51_major_data(iCH4,:,:,:) = k51_ch4(:,:,:)
   !=====  interval 52:   =====!
   k52_major_data(iH2O,:,:,:) = k52_h2o(:,:,:) 
   k52_major_data(iCO2,:,:,:) = k52_co2(:,:,:)
   k52_major_data(iCH4,:,:,:) = k52_ch4(:,:,:)
   !=====  interval 53:   =====!
   k53_major_data(iH2O,:,:,:) = k53_h2o(:,:,:) 
   k53_major_data(iCO2,:,:,:) = k53_co2(:,:,:)
   k53_major_data(iCH4,:,:,:) = k53_ch4(:,:,:)
   !=====  interval 54:   =====!
   k54_major_data(iH2O,:,:,:) = k54_h2o(:,:,:) 
   k54_major_data(iCO2,:,:,:) = k54_co2(:,:,:)
   k54_major_data(iCH4,:,:,:) = k54_ch4(:,:,:)
   !=====  interval 55:   =====!
   k55_major_data(iH2O,:,:,:) = k55_h2o(:,:,:) 
   k55_major_data(iCO2,:,:,:) = k55_co2(:,:,:)
   k55_major_data(iCH4,:,:,:) = k55_ch4(:,:,:)
   !=====  interval 56:   =====!
   k56_major_data(iH2O,:,:,:) = k56_h2o(:,:,:) 
   k56_major_data(iCO2,:,:,:) = k56_co2(:,:,:)
   k56_major_data(iCH4,:,:,:) = k56_ch4(:,:,:)
   !=====  interval 57:   =====!
   k57_major_data(iH2O,:,:,:) = k57_h2o(:,:,:) 
   k57_major_data(iCO2,:,:,:) = k57_co2(:,:,:)
   k57_major_data(iCH4,:,:,:) = k57_ch4(:,:,:)
   !=====  interval 58:   =====!
   k58_major_data(iH2O,:,:,:) = k58_h2o(:,:,:) 
   k58_major_data(iCO2,:,:,:) = k58_co2(:,:,:)
   k58_major_data(iCH4,:,:,:) = k58_ch4(:,:,:)
   !=====  interval 59:   =====!
   k59_major_data(iH2O,:,:,:) = k59_h2o(:,:,:) 
   k59_major_data(iCO2,:,:,:) = k59_co2(:,:,:)
   k59_major_data(iCH4,:,:,:) = k59_ch4(:,:,:)
   !=====  interval 60:   =====!
   k60_major_data(iH2O,:,:,:) = k60_h2o(:,:,:) 
   k60_major_data(iCO2,:,:,:) = k60_co2(:,:,:)
   k60_major_data(iCH4,:,:,:) = k60_ch4(:,:,:)
   !=====  interval 61:   =====!
   k61_major_data(iH2O,:,:,:) = k61_h2o(:,:,:) 
   k61_major_data(iCO2,:,:,:) = k61_co2(:,:,:)
   k61_major_data(iCH4,:,:,:) = k61_ch4(:,:,:)
   !=====  interval 62:   =====!
   k62_major_data(iH2O,:,:,:) = k62_h2o(:,:,:) 
   k62_major_data(iCO2,:,:,:) = k62_co2(:,:,:)
   k62_major_data(iCH4,:,:,:) = k62_ch4(:,:,:)
   !=====  interval 63:   =====!
   k63_major_data(iH2O,:,:,:) = k63_h2o(:,:,:) 
   k63_major_data(iCO2,:,:,:) = k63_co2(:,:,:)
   k63_major_data(iCH4,:,:,:) = k63_ch4(:,:,:)
   !=====  interval 64:   =====!
   k64_major_data(iH2O,:,:,:) = k64_h2o(:,:,:) 
   k64_major_data(iCO2,:,:,:) = k64_co2(:,:,:)
   k64_major_data(iCH4,:,:,:) = k64_ch4(:,:,:)
   !=====  interval 65:   =====!
   k65_major_data(iH2O,:,:,:) = k65_h2o(:,:,:) 
   k65_major_data(iCO2,:,:,:) = k65_co2(:,:,:)
   k65_major_data(iCH4,:,:,:) = k65_ch4(:,:,:)
   !=====  interval 66:   =====!
   k66_major_data(iH2O,:,:,:) = k66_h2o(:,:,:) 
   k66_major_data(iCO2,:,:,:) = k66_co2(:,:,:)
   k66_major_data(iCH4,:,:,:) = k66_ch4(:,:,:)
   !=====  interval 67:   =====!
   k67_major_data(iH2O,:,:,:) = k67_h2o(:,:,:) 
   k67_major_data(iCO2,:,:,:) = k67_co2(:,:,:)
   k67_major_data(iCH4,:,:,:) = k67_ch4(:,:,:)
   !=====  interval 68:   =====!
   k68_major_data(iH2O,:,:,:) = k68_h2o(:,:,:) 
   k68_major_data(iCO2,:,:,:) = k68_co2(:,:,:)
   k68_major_data(iCH4,:,:,:) = k68_ch4(:,:,:)

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
!    use exoplanet_mod  ! this marks a divergence between 1d and 3d implementations
    ! i want to read mixing ratios from here. 
    implicit none

!------------------------------------------------------------------------ 
!
! Local Variables
!
  integer :: iq, ig


!------------------------------------------------------------------------ 
! Start Code
!

   ! Create grey gas absorption matrix
   iq=0
   !=====  interval 1:   =====!
   k01_grey_data(:,:,:) = 0.0
   do ig=1, ngauss_pts(1)
     iq = iq + 1
     k01_grey_data(iH2O,:,:) = k01_grey_data(iH2O,:,:) + k01_h2o(ig,:,:) * g_weight(iq) 
     k01_grey_data(iCO2,:,:) = k01_grey_data(iCO2,:,:) + k01_co2(ig,:,:) * g_weight(iq)
     k01_grey_data(iCH4,:,:) = k01_grey_data(iCH4,:,:) + k01_ch4(ig,:,:) * g_weight(iq)
   enddo
   !=====  interval 2:   =====!
   k02_grey_data(:,:,:) = 0.0
   do ig=1, ngauss_pts(2)
     iq = iq + 1
     k02_grey_data(iH2O,:,:) = k02_grey_data(iH2O,:,:) + k02_h2o(ig,:,:) * g_weight(iq)
     k02_grey_data(iCO2,:,:) = k02_grey_data(iCO2,:,:) + k02_co2(ig,:,:) * g_weight(iq)
     k02_grey_data(iCH4,:,:) = k02_grey_data(iCH4,:,:) + k02_ch4(ig,:,:) * g_weight(iq)
   enddo
   !=====  interval 3:   =====!
   k03_grey_data(:,:,:) = 0.0
   do ig=1, ngauss_pts(3)
     iq = iq + 1
     k03_grey_data(iH2O,:,:) = k03_grey_data(iH2O,:,:) + k03_h2o(ig,:,:) * g_weight(iq)
     k03_grey_data(iCO2,:,:) = k03_grey_data(iCO2,:,:) + k03_co2(ig,:,:) * g_weight(iq)
     k03_grey_data(iCH4,:,:) = k03_grey_data(iCH4,:,:) + k03_ch4(ig,:,:) * g_weight(iq)
   enddo
   !=====  interval 4:   =====!
   k04_grey_data(:,:,:) = 0.0
   do ig=1, ngauss_pts(4)
     iq = iq + 1
     k04_grey_data(iH2O,:,:) = k04_grey_data(iH2O,:,:) + k04_h2o(ig,:,:) * g_weight(iq)
     k04_grey_data(iCO2,:,:) = k04_grey_data(iCO2,:,:) + k04_co2(ig,:,:) * g_weight(iq)
     k04_grey_data(iCH4,:,:) = k04_grey_data(iCH4,:,:) + k04_ch4(ig,:,:) * g_weight(iq)
   enddo
   !=====  interval 5:   =====!
   k05_grey_data(:,:,:) = 0.0
   do ig=1, ngauss_pts(5)
     iq = iq + 1
     k05_grey_data(iH2O,:,:) = k05_grey_data(iH2O,:,:) + k05_h2o(ig,:,:) * g_weight(iq)
     k05_grey_data(iCO2,:,:) = k05_grey_data(iCO2,:,:) + k05_co2(ig,:,:) * g_weight(iq)
     k05_grey_data(iCH4,:,:) = k05_grey_data(iCH4,:,:) + k05_ch4(ig,:,:) * g_weight(iq)
   enddo
   !=====  interval 6:   =====!
   k06_grey_data(:,:,:) = 0.0
   do ig=1, ngauss_pts(6)
     iq = iq + 1
     k06_grey_data(iH2O,:,:) = k06_grey_data(iH2O,:,:) + k06_h2o(ig,:,:) * g_weight(iq)
     k06_grey_data(iCO2,:,:) = k06_grey_data(iCO2,:,:) + k06_co2(ig,:,:) * g_weight(iq)
     k06_grey_data(iCH4,:,:) = k06_grey_data(iCH4,:,:) + k06_ch4(ig,:,:) * g_weight(iq)
   enddo
   !=====  interval 7:   =====!
   k07_grey_data(:,:,:) = 0.0
   do ig=1, ngauss_pts(7)
     iq = iq + 1
     k07_grey_data(iH2O,:,:) = k07_grey_data(iH2O,:,:) + k07_h2o(ig,:,:) * g_weight(iq)
     k07_grey_data(iCO2,:,:) = k07_grey_data(iCO2,:,:) + k07_co2(ig,:,:) * g_weight(iq)
     k07_grey_data(iCH4,:,:) = k07_grey_data(iCH4,:,:) + k07_ch4(ig,:,:) * g_weight(iq)
   enddo
   !=====  interval 8:   =====!
   k08_grey_data(:,:,:) = 0.0
   do ig=1, ngauss_pts(8)
     iq = iq + 1
     k08_grey_data(iH2O,:,:) = k08_grey_data(iH2O,:,:) + k08_h2o(ig,:,:) * g_weight(iq)
     k08_grey_data(iCO2,:,:) = k08_grey_data(iCO2,:,:) + k08_co2(ig,:,:) * g_weight(iq)
     k08_grey_data(iCH4,:,:) = k08_grey_data(iCH4,:,:) + k08_ch4(ig,:,:) * g_weight(iq)
   enddo
   !=====  interval 9:   =====!
   k09_grey_data(:,:,:) = 0.0
   do ig=1, ngauss_pts(9)
     iq = iq + 1
     k09_grey_data(iH2O,:,:) = k09_grey_data(iH2O,:,:) + k09_h2o(ig,:,:) * g_weight(iq)
     k09_grey_data(iCO2,:,:) = k09_grey_data(iCO2,:,:) + k09_co2(ig,:,:) * g_weight(iq)
     k09_grey_data(iCH4,:,:) = k09_grey_data(iCH4,:,:) + k09_ch4(ig,:,:) * g_weight(iq)
   enddo
   !=====  interval 10:   =====!
   k10_grey_data(:,:,:) = 0.0
   do ig=1, ngauss_pts(10)
     iq = iq + 1
     k10_grey_data(iH2O,:,:) = k10_grey_data(iH2O,:,:) + k10_h2o(ig,:,:) * g_weight(iq)
     k10_grey_data(iCO2,:,:) = k10_grey_data(iCO2,:,:) + k10_co2(ig,:,:) * g_weight(iq)
     k10_grey_data(iCH4,:,:) = k10_grey_data(iCH4,:,:) + k10_ch4(ig,:,:) * g_weight(iq)
   enddo
   !=====  interval 11:   =====!
   k11_grey_data(:,:,:) = 0.0
   do ig=1, ngauss_pts(11)
     iq = iq + 1
     k11_grey_data(iH2O,:,:) = k11_grey_data(iH2O,:,:) + k11_h2o(ig,:,:) * g_weight(iq)
     k11_grey_data(iCO2,:,:) = k11_grey_data(iCO2,:,:) + k11_co2(ig,:,:) * g_weight(iq)
     k11_grey_data(iCH4,:,:) = k11_grey_data(iCH4,:,:) + k11_ch4(ig,:,:) * g_weight(iq)
   enddo
   !=====  interval 12:   =====!
   k12_grey_data(:,:,:) = 0.0
   do ig=1, ngauss_pts(12)
     iq = iq + 1
     k12_grey_data(iH2O,:,:) = k12_grey_data(iH2O,:,:) + k12_h2o(ig,:,:) * g_weight(iq)
     k12_grey_data(iCO2,:,:) = k12_grey_data(iCO2,:,:) + k12_co2(ig,:,:) * g_weight(iq)
     k12_grey_data(iCH4,:,:) = k12_grey_data(iCH4,:,:) + k12_ch4(ig,:,:) * g_weight(iq)
   enddo
   !=====  interval 13:   =====!
   k13_grey_data(:,:,:) = 0.0
   do ig=1, ngauss_pts(13)
     iq = iq + 1
     k13_grey_data(iH2O,:,:) = k13_grey_data(iH2O,:,:) + k13_h2o(ig,:,:) * g_weight(iq)
     k13_grey_data(iCO2,:,:) = k13_grey_data(iCO2,:,:) + k13_co2(ig,:,:) * g_weight(iq)
     k13_grey_data(iCH4,:,:) = k13_grey_data(iCH4,:,:) + k13_ch4(ig,:,:) * g_weight(iq)
   enddo
   !=====  interval 14:   =====!
   k14_grey_data(:,:,:) = 0.0
   do ig=1, ngauss_pts(14)
     iq = iq + 1
     k14_grey_data(iH2O,:,:) = k14_grey_data(iH2O,:,:) + k14_h2o(ig,:,:) * g_weight(iq)
     k14_grey_data(iCO2,:,:) = k14_grey_data(iCO2,:,:) + k14_co2(ig,:,:) * g_weight(iq)
     k14_grey_data(iCH4,:,:) = k14_grey_data(iCH4,:,:) + k14_ch4(ig,:,:) * g_weight(iq)
   enddo
   !=====  interval 15:   =====!
   k15_grey_data(:,:,:) = 0.0
   do ig=1, ngauss_pts(15)
     iq = iq + 1
     k15_grey_data(iH2O,:,:) = k15_grey_data(iH2O,:,:) + k15_h2o(ig,:,:) * g_weight(iq)
     k15_grey_data(iCO2,:,:) = k15_grey_data(iCO2,:,:) + k15_co2(ig,:,:) * g_weight(iq)
     k15_grey_data(iCH4,:,:) = k15_grey_data(iCH4,:,:) + k15_ch4(ig,:,:) * g_weight(iq)
   enddo
   !=====  interval 16:   =====!
   k16_grey_data(:,:,:) = 0.0
   do ig=1, ngauss_pts(16)
     iq = iq + 1
     k16_grey_data(iH2O,:,:) = k16_grey_data(iH2O,:,:) + k16_h2o(ig,:,:) * g_weight(iq)
     k16_grey_data(iCO2,:,:) = k16_grey_data(iCO2,:,:) + k16_co2(ig,:,:) * g_weight(iq)
     k16_grey_data(iCH4,:,:) = k16_grey_data(iCH4,:,:) + k16_ch4(ig,:,:) * g_weight(iq)
   enddo
   !=====  interval 17:   =====!
   k17_grey_data(:,:,:) = 0.0
   do ig=1, ngauss_pts(17)
     iq = iq + 1
     k17_grey_data(iH2O,:,:) = k17_grey_data(iH2O,:,:) + k17_h2o(ig,:,:) * g_weight(iq)
     k17_grey_data(iCO2,:,:) = k17_grey_data(iCO2,:,:) + k17_co2(ig,:,:) * g_weight(iq)
     k17_grey_data(iCH4,:,:) = k17_grey_data(iCH4,:,:) + k17_ch4(ig,:,:) * g_weight(iq)
   enddo
   !=====  interval 18:   =====!
   k18_grey_data(:,:,:) = 0.0
   do ig=1, ngauss_pts(18)
     iq = iq + 1
     k18_grey_data(iH2O,:,:) = k18_grey_data(iH2O,:,:) + k18_h2o(ig,:,:) * g_weight(iq)
     k18_grey_data(iCO2,:,:) = k18_grey_data(iCO2,:,:) + k18_co2(ig,:,:) * g_weight(iq)
     k18_grey_data(iCH4,:,:) = k18_grey_data(iCH4,:,:) + k18_ch4(ig,:,:) * g_weight(iq)
   enddo
   !=====  interval 19:   =====!
   k19_grey_data(:,:,:) = 0.0
   do ig=1, ngauss_pts(19)
     iq = iq + 1
     k19_grey_data(iH2O,:,:) = k19_grey_data(iH2O,:,:) + k19_h2o(ig,:,:) * g_weight(iq)
     k19_grey_data(iCO2,:,:) = k19_grey_data(iCO2,:,:) + k19_co2(ig,:,:) * g_weight(iq)
     k19_grey_data(iCH4,:,:) = k19_grey_data(iCH4,:,:) + k19_ch4(ig,:,:) * g_weight(iq)
   enddo
   !=====  interval 20:   =====!
   k20_grey_data(:,:,:) = 0.0
   do ig=1, ngauss_pts(20)
     iq = iq + 1
     k20_grey_data(iH2O,:,:) = k20_grey_data(iH2O,:,:) + k20_h2o(ig,:,:) * g_weight(iq)
     k20_grey_data(iCO2,:,:) = k20_grey_data(iCO2,:,:) + k20_co2(ig,:,:) * g_weight(iq)
     k20_grey_data(iCH4,:,:) = k20_grey_data(iCH4,:,:) + k20_ch4(ig,:,:) * g_weight(iq)
   enddo
   !=====  interval 21:   =====!
   k21_grey_data(:,:,:) = 0.0
   do ig=1, ngauss_pts(21)
     iq = iq + 1
     k21_grey_data(iH2O,:,:) = k21_grey_data(iH2O,:,:) + k21_h2o(ig,:,:) * g_weight(iq)
     k21_grey_data(iCO2,:,:) = k21_grey_data(iCO2,:,:) + k21_co2(ig,:,:) * g_weight(iq)
     k21_grey_data(iCH4,:,:) = k21_grey_data(iCH4,:,:) + k21_ch4(ig,:,:) * g_weight(iq)
   enddo
   !=====  interval 22:   =====!
   k22_grey_data(:,:,:) = 0.0
   do ig=1, ngauss_pts(22)
     iq = iq + 1
     k22_grey_data(iH2O,:,:) = k22_grey_data(iH2O,:,:) + k22_h2o(ig,:,:) * g_weight(iq)
     k22_grey_data(iCO2,:,:) = k22_grey_data(iCO2,:,:) + k22_co2(ig,:,:) * g_weight(iq)
     k22_grey_data(iCH4,:,:) = k22_grey_data(iCH4,:,:) + k22_ch4(ig,:,:) * g_weight(iq)
   enddo
   !=====  interval 23:   =====!
   k23_grey_data(:,:,:) = 0.0
   do ig=1, ngauss_pts(23)
     iq = iq + 1
     k23_grey_data(iH2O,:,:) = k23_grey_data(iH2O,:,:) + k23_h2o(ig,:,:) * g_weight(iq)
     k23_grey_data(iCO2,:,:) = k23_grey_data(iCO2,:,:) + k23_co2(ig,:,:) * g_weight(iq)
     k23_grey_data(iCH4,:,:) = k23_grey_data(iCH4,:,:) + k23_ch4(ig,:,:) * g_weight(iq)
   enddo
   !=====  interval 24:   =====!
   k24_grey_data(:,:,:) = 0.0
   do ig=1, ngauss_pts(24)
     iq = iq + 1
     k24_grey_data(iH2O,:,:) = k24_grey_data(iH2O,:,:) + k24_h2o(ig,:,:) * g_weight(iq)
     k24_grey_data(iCO2,:,:) = k24_grey_data(iCO2,:,:) + k24_co2(ig,:,:) * g_weight(iq)
     k24_grey_data(iCH4,:,:) = k24_grey_data(iCH4,:,:) + k24_ch4(ig,:,:) * g_weight(iq)
   enddo
   !=====  interval 25:   =====!
   k25_grey_data(:,:,:) = 0.0
   do ig=1, ngauss_pts(25)
     iq = iq + 1
     k25_grey_data(iH2O,:,:) = k25_grey_data(iH2O,:,:) + k25_h2o(ig,:,:) * g_weight(iq)
     k25_grey_data(iCO2,:,:) = k25_grey_data(iCO2,:,:) + k25_co2(ig,:,:) * g_weight(iq)
     k25_grey_data(iCH4,:,:) = k25_grey_data(iCH4,:,:) + k25_ch4(ig,:,:) * g_weight(iq)
   enddo
   !=====  interval 26:   =====!
   k26_grey_data(:,:,:) = 0.0
   do ig=1, ngauss_pts(26)
     iq = iq + 1
     k26_grey_data(iH2O,:,:) = k26_grey_data(iH2O,:,:) + k26_h2o(ig,:,:) * g_weight(iq)
     k26_grey_data(iCO2,:,:) = k26_grey_data(iCO2,:,:) + k26_co2(ig,:,:) * g_weight(iq)
     k26_grey_data(iCH4,:,:) = k26_grey_data(iCH4,:,:) + k26_ch4(ig,:,:) * g_weight(iq)
   enddo
   !=====  interval 27:   =====!
   k27_grey_data(:,:,:) = 0.0
   do ig=1, ngauss_pts(27)
     iq = iq + 1
     k27_grey_data(iH2O,:,:) = k27_grey_data(iH2O,:,:) + k27_h2o(ig,:,:) * g_weight(iq)
     k27_grey_data(iCO2,:,:) = k27_grey_data(iCO2,:,:) + k27_co2(ig,:,:) * g_weight(iq)
     k27_grey_data(iCH4,:,:) = k27_grey_data(iCH4,:,:) + k27_ch4(ig,:,:) * g_weight(iq)
   enddo
   !=====  interval 28:   =====!
   k28_grey_data(:,:,:) = 0.0
   do ig=1, ngauss_pts(28)
     iq = iq + 1
     k28_grey_data(iH2O,:,:) = k28_grey_data(iH2O,:,:) + k28_h2o(ig,:,:) * g_weight(iq)
     k28_grey_data(iCO2,:,:) = k28_grey_data(iCO2,:,:) + k28_co2(ig,:,:) * g_weight(iq)
     k28_grey_data(iCH4,:,:) = k28_grey_data(iCH4,:,:) + k28_ch4(ig,:,:) * g_weight(iq)
   enddo
   !=====  interval 29:   =====!
   k29_grey_data(:,:,:) = 0.0
   do ig=1, ngauss_pts(29)
     iq = iq + 1
     k29_grey_data(iH2O,:,:) = k29_grey_data(iH2O,:,:) + k29_h2o(ig,:,:) * g_weight(iq)
     k29_grey_data(iCO2,:,:) = k29_grey_data(iCO2,:,:) + k29_co2(ig,:,:) * g_weight(iq)
     k29_grey_data(iCH4,:,:) = k29_grey_data(iCH4,:,:) + k29_ch4(ig,:,:) * g_weight(iq)
   enddo
   !=====  interval 30:   =====!
   k30_grey_data(:,:,:) = 0.0
   do ig=1, ngauss_pts(30)
     iq = iq + 1
     k30_grey_data(iH2O,:,:) = k30_grey_data(iH2O,:,:) + k30_h2o(ig,:,:) * g_weight(iq)
     k30_grey_data(iCO2,:,:) = k30_grey_data(iCO2,:,:) + k30_co2(ig,:,:) * g_weight(iq)
     k30_grey_data(iCH4,:,:) = k30_grey_data(iCH4,:,:) + k30_ch4(ig,:,:) * g_weight(iq)
   enddo
   !=====  interval 31:   =====!
   k31_grey_data(:,:,:) = 0.0
   do ig=1, ngauss_pts(31)
     iq = iq + 1
     k31_grey_data(iH2O,:,:) = k31_grey_data(iH2O,:,:) + k31_h2o(ig,:,:) * g_weight(iq)
     k31_grey_data(iCO2,:,:) = k31_grey_data(iCO2,:,:) + k31_co2(ig,:,:) * g_weight(iq)
     k31_grey_data(iCH4,:,:) = k31_grey_data(iCH4,:,:) + k31_ch4(ig,:,:) * g_weight(iq)
   enddo
   !=====  interval 32:   =====!
   k32_grey_data(:,:,:) = 0.0
   do ig=1, ngauss_pts(32)
     iq = iq + 1
     k32_grey_data(iH2O,:,:) = k32_grey_data(iH2O,:,:) + k32_h2o(ig,:,:) * g_weight(iq)
     k32_grey_data(iCO2,:,:) = k32_grey_data(iCO2,:,:) + k32_co2(ig,:,:) * g_weight(iq)
     k32_grey_data(iCH4,:,:) = k32_grey_data(iCH4,:,:) + k32_ch4(ig,:,:) * g_weight(iq)
   enddo
   !=====  interval 33:   =====!
   k33_grey_data(:,:,:) = 0.0
   do ig=1, ngauss_pts(33)
     iq = iq + 1
     k33_grey_data(iH2O,:,:) = k33_grey_data(iH2O,:,:) + k33_h2o(ig,:,:) * g_weight(iq)
     k33_grey_data(iCO2,:,:) = k33_grey_data(iCO2,:,:) + k33_co2(ig,:,:) * g_weight(iq)
     k33_grey_data(iCH4,:,:) = k33_grey_data(iCH4,:,:) + k33_ch4(ig,:,:) * g_weight(iq)
   enddo
   !=====  interval 34:   =====!
   k34_grey_data(:,:,:) = 0.0
   do ig=1, ngauss_pts(34)
     iq = iq + 1
     k34_grey_data(iH2O,:,:) = k34_grey_data(iH2O,:,:) + k34_h2o(ig,:,:) * g_weight(iq)
     k34_grey_data(iCO2,:,:) = k34_grey_data(iCO2,:,:) + k34_co2(ig,:,:) * g_weight(iq)
     k34_grey_data(iCH4,:,:) = k34_grey_data(iCH4,:,:) + k34_ch4(ig,:,:) * g_weight(iq)
   enddo
   !=====  interval 35:   =====!
   k35_grey_data(:,:,:) = 0.0
   do ig=1, ngauss_pts(35)
     iq = iq + 1
     k35_grey_data(iH2O,:,:) = k35_grey_data(iH2O,:,:) + k35_h2o(ig,:,:) * g_weight(iq)
     k35_grey_data(iCO2,:,:) = k35_grey_data(iCO2,:,:) + k35_co2(ig,:,:) * g_weight(iq)
     k35_grey_data(iCH4,:,:) = k35_grey_data(iCH4,:,:) + k35_ch4(ig,:,:) * g_weight(iq)
   enddo
   !=====  interval 36:   =====!
   k36_grey_data(:,:,:) = 0.0
   do ig=1, ngauss_pts(36)
     iq = iq + 1
     k36_grey_data(iH2O,:,:) = k36_grey_data(iH2O,:,:) + k36_h2o(ig,:,:) * g_weight(iq)
     k36_grey_data(iCO2,:,:) = k36_grey_data(iCO2,:,:) + k36_co2(ig,:,:) * g_weight(iq)
     k36_grey_data(iCH4,:,:) = k36_grey_data(iCH4,:,:) + k36_ch4(ig,:,:) * g_weight(iq)
   enddo
   !=====  interval 37:   =====!
   k37_grey_data(:,:,:) = 0.0
   do ig=1, ngauss_pts(37)
     iq = iq + 1
     k37_grey_data(iH2O,:,:) = k37_grey_data(iH2O,:,:) + k37_h2o(ig,:,:) * g_weight(iq)
     k37_grey_data(iCO2,:,:) = k37_grey_data(iCO2,:,:) + k37_co2(ig,:,:) * g_weight(iq)
     k37_grey_data(iCH4,:,:) = k37_grey_data(iCH4,:,:) + k37_ch4(ig,:,:) * g_weight(iq)
   enddo
   !=====  interval 38:   =====!
   k38_grey_data(:,:,:) = 0.0
   do ig=1, ngauss_pts(38)
     iq = iq + 1
     k38_grey_data(iH2O,:,:) = k38_grey_data(iH2O,:,:) + k38_h2o(ig,:,:) * g_weight(iq)
     k38_grey_data(iCO2,:,:) = k38_grey_data(iCO2,:,:) + k38_co2(ig,:,:) * g_weight(iq)
     k38_grey_data(iCH4,:,:) = k38_grey_data(iCH4,:,:) + k38_ch4(ig,:,:) * g_weight(iq)
   enddo
   !=====  interval 39:   =====!
   k39_grey_data(:,:,:) = 0.0
   do ig=1, ngauss_pts(39)
     iq = iq + 1
     k39_grey_data(iH2O,:,:) = k39_grey_data(iH2O,:,:) + k39_h2o(ig,:,:) * g_weight(iq)
     k39_grey_data(iCO2,:,:) = k39_grey_data(iCO2,:,:) + k39_co2(ig,:,:) * g_weight(iq)
     k39_grey_data(iCH4,:,:) = k39_grey_data(iCH4,:,:) + k39_ch4(ig,:,:) * g_weight(iq)
   enddo
   !=====  interval 40:   =====!
   k40_grey_data(:,:,:) = 0.0
   do ig=1, ngauss_pts(40)
     iq = iq + 1
     k40_grey_data(iH2O,:,:) = k40_grey_data(iH2O,:,:) + k40_h2o(ig,:,:) * g_weight(iq)
     k40_grey_data(iCO2,:,:) = k40_grey_data(iCO2,:,:) + k40_co2(ig,:,:) * g_weight(iq)
     k40_grey_data(iCH4,:,:) = k40_grey_data(iCH4,:,:) + k40_ch4(ig,:,:) * g_weight(iq)
   enddo
   !=====  interval 41:   =====!
   k41_grey_data(:,:,:) = 0.0
   do ig=1, ngauss_pts(41)
     iq = iq + 1
     k41_grey_data(iH2O,:,:) = k41_grey_data(iH2O,:,:) + k41_h2o(ig,:,:) * g_weight(iq)
     k41_grey_data(iCO2,:,:) = k41_grey_data(iCO2,:,:) + k41_co2(ig,:,:) * g_weight(iq)
     k41_grey_data(iCH4,:,:) = k41_grey_data(iCH4,:,:) + k41_ch4(ig,:,:) * g_weight(iq)
   enddo
   !=====  interval 42:   =====!
   k42_grey_data(:,:,:) = 0.0
   do ig=1, ngauss_pts(42)
     iq = iq + 1
     k42_grey_data(iH2O,:,:) = k42_grey_data(iH2O,:,:) + k42_h2o(ig,:,:) * g_weight(iq)
     k42_grey_data(iCO2,:,:) = k42_grey_data(iCO2,:,:) + k42_co2(ig,:,:) * g_weight(iq)
     k42_grey_data(iCH4,:,:) = k42_grey_data(iCH4,:,:) + k42_ch4(ig,:,:) * g_weight(iq)
   enddo
   !=====  interval 43:   =====!
   k43_grey_data(:,:,:) = 0.0
   do ig=1, ngauss_pts(43)
     iq = iq + 1
     k43_grey_data(iH2O,:,:) = k43_grey_data(iH2O,:,:) + k43_h2o(ig,:,:) * g_weight(iq)
     k43_grey_data(iCO2,:,:) = k43_grey_data(iCO2,:,:) + k43_co2(ig,:,:) * g_weight(iq)
     k43_grey_data(iCH4,:,:) = k43_grey_data(iCH4,:,:) + k43_ch4(ig,:,:) * g_weight(iq)
   enddo
   !=====  interval 44:   =====!
   k44_grey_data(:,:,:) = 0.0
   do ig=1, ngauss_pts(44)
     iq = iq + 1
     k44_grey_data(iH2O,:,:) = k44_grey_data(iH2O,:,:) + k44_h2o(ig,:,:) * g_weight(iq)
     k44_grey_data(iCO2,:,:) = k44_grey_data(iCO2,:,:) + k44_co2(ig,:,:) * g_weight(iq)
     k44_grey_data(iCH4,:,:) = k44_grey_data(iCH4,:,:) + k44_ch4(ig,:,:) * g_weight(iq)
   enddo
   !=====  interval 45:  =====!
   k45_grey_data(:,:,:) = 0.0
   do ig=1, ngauss_pts(45)
     iq = iq + 1
     k45_grey_data(iH2O,:,:) = k45_grey_data(iH2O,:,:) + k45_h2o(ig,:,:) * g_weight(iq)
     k45_grey_data(iCO2,:,:) = k45_grey_data(iCO2,:,:) + k45_co2(ig,:,:) * g_weight(iq)
     k45_grey_data(iCH4,:,:) = k45_grey_data(iCH4,:,:) + k45_ch4(ig,:,:) * g_weight(iq)
   enddo
   !=====  interval 46: 6 =====!
   k46_grey_data(:,:,:) = 0.0
   do ig=1, ngauss_pts(40)
     iq = iq + 1
     k46_grey_data(iH2O,:,:) = k46_grey_data(iH2O,:,:) + k46_h2o(ig,:,:) * g_weight(iq)
     k46_grey_data(iCO2,:,:) = k46_grey_data(iCO2,:,:) + k46_co2(ig,:,:) * g_weight(iq)
     k46_grey_data(iCH4,:,:) = k46_grey_data(iCH4,:,:) + k46_ch4(ig,:,:) * g_weight(iq)
   enddo
   !=====  interval 47:   =====!
   k47_grey_data(:,:,:) = 0.0
   do ig=1, ngauss_pts(47)
     iq = iq + 1
     k47_grey_data(iH2O,:,:) = k47_grey_data(iH2O,:,:) + k47_h2o(ig,:,:) * g_weight(iq)
     k47_grey_data(iCO2,:,:) = k47_grey_data(iCO2,:,:) + k47_co2(ig,:,:) * g_weight(iq)
     k47_grey_data(iCH4,:,:) = k47_grey_data(iCH4,:,:) + k47_ch4(ig,:,:) * g_weight(iq)
   enddo
   !=====  interval 48:   =====!
   k48_grey_data(:,:,:) = 0.0
   do ig=1, ngauss_pts(48)
     iq = iq + 1
     k48_grey_data(iH2O,:,:) = k48_grey_data(iH2O,:,:) + k48_h2o(ig,:,:) * g_weight(iq)
     k48_grey_data(iCO2,:,:) = k48_grey_data(iCO2,:,:) + k48_co2(ig,:,:) * g_weight(iq)
     k48_grey_data(iCH4,:,:) = k48_grey_data(iCH4,:,:) + k48_ch4(ig,:,:) * g_weight(iq)
   enddo
   !=====  interval 49:   =====!
   k49_grey_data(:,:,:) = 0.0
   do ig=1, ngauss_pts(49)
     iq = iq + 1
     k49_grey_data(iH2O,:,:) = k49_grey_data(iH2O,:,:) + k49_h2o(ig,:,:) * g_weight(iq)
     k49_grey_data(iCO2,:,:) = k49_grey_data(iCO2,:,:) + k49_co2(ig,:,:) * g_weight(iq)
     k49_grey_data(iCH4,:,:) = k49_grey_data(iCH4,:,:) + k49_ch4(ig,:,:) * g_weight(iq)
   enddo
   !=====  interval 50:   =====!
   k50_grey_data(:,:,:) = 0.0
   do ig=1, ngauss_pts(50)
     iq = iq + 1
     k50_grey_data(iH2O,:,:) = k50_grey_data(iH2O,:,:) + k50_h2o(ig,:,:) * g_weight(iq)
     k50_grey_data(iCO2,:,:) = k50_grey_data(iCO2,:,:) + k50_co2(ig,:,:) * g_weight(iq)
     k50_grey_data(iCH4,:,:) = k50_grey_data(iCH4,:,:) + k50_ch4(ig,:,:) * g_weight(iq)
   enddo
   !=====  interval 51:   =====!
   k51_grey_data(:,:,:) = 0.0
   do ig=1, ngauss_pts(51)
     iq = iq + 1
     k51_grey_data(iH2O,:,:) = k51_grey_data(iH2O,:,:) + k51_h2o(ig,:,:) * g_weight(iq)
     k51_grey_data(iCO2,:,:) = k51_grey_data(iCO2,:,:) + k51_co2(ig,:,:) * g_weight(iq)
     k51_grey_data(iCH4,:,:) = k51_grey_data(iCH4,:,:) + k51_ch4(ig,:,:) * g_weight(iq)
   enddo
   !=====  interval 52:   =====!
   k52_grey_data(:,:,:) = 0.0
   do ig=1, ngauss_pts(52)
     iq = iq + 1
     k52_grey_data(iH2O,:,:) = k52_grey_data(iH2O,:,:) + k52_h2o(ig,:,:) * g_weight(iq)
     k52_grey_data(iCO2,:,:) = k52_grey_data(iCO2,:,:) + k52_co2(ig,:,:) * g_weight(iq)
     k52_grey_data(iCH4,:,:) = k52_grey_data(iCH4,:,:) + k52_ch4(ig,:,:) * g_weight(iq)
   enddo
   !=====  interval 53:   =====!
   k53_grey_data(:,:,:) = 0.0
   do ig=1, ngauss_pts(53)
     iq = iq + 1
     k53_grey_data(iH2O,:,:) = k53_grey_data(iH2O,:,:) + k53_h2o(ig,:,:) * g_weight(iq)
     k53_grey_data(iCO2,:,:) = k53_grey_data(iCO2,:,:) + k53_co2(ig,:,:) * g_weight(iq)
     k53_grey_data(iCH4,:,:) = k53_grey_data(iCH4,:,:) + k53_ch4(ig,:,:) * g_weight(iq)
   enddo
   !=====  interval 54:   =====!
   k54_grey_data(:,:,:) = 0.0
   do ig=1, ngauss_pts(54)
     iq = iq + 1
     k54_grey_data(iH2O,:,:) = k54_grey_data(iH2O,:,:) + k54_h2o(ig,:,:) * g_weight(iq)
     k54_grey_data(iCO2,:,:) = k54_grey_data(iCO2,:,:) + k54_co2(ig,:,:) * g_weight(iq)
     k54_grey_data(iCH4,:,:) = k54_grey_data(iCH4,:,:) + k54_ch4(ig,:,:) * g_weight(iq)
   enddo
   !=====  interval 55:   =====!
   k55_grey_data(:,:,:) = 0.0
   do ig=1, ngauss_pts(55)
     iq = iq + 1
     k55_grey_data(iH2O,:,:) = k55_grey_data(iH2O,:,:) + k55_h2o(ig,:,:) * g_weight(iq)
     k55_grey_data(iCO2,:,:) = k55_grey_data(iCO2,:,:) + k55_co2(ig,:,:) * g_weight(iq)
     k55_grey_data(iCH4,:,:) = k55_grey_data(iCH4,:,:) + k55_ch4(ig,:,:) * g_weight(iq)
   enddo
   !=====  interval 56:   =====!
   k56_grey_data(:,:,:) = 0.0
   do ig=1, ngauss_pts(56)
     iq = iq + 1
     k56_grey_data(iH2O,:,:) = k56_grey_data(iH2O,:,:) + k56_h2o(ig,:,:) * g_weight(iq)
     k56_grey_data(iCO2,:,:) = k56_grey_data(iCO2,:,:) + k56_co2(ig,:,:) * g_weight(iq)
     k56_grey_data(iCH4,:,:) = k56_grey_data(iCH4,:,:) + k56_ch4(ig,:,:) * g_weight(iq)
   enddo
   !=====  interval 57:   =====!
   k57_grey_data(:,:,:) = 0.0
   do ig=1, ngauss_pts(57)
     iq = iq + 1
     k57_grey_data(iH2O,:,:) = k57_grey_data(iH2O,:,:) + k57_h2o(ig,:,:) * g_weight(iq)
     k57_grey_data(iCO2,:,:) = k57_grey_data(iCO2,:,:) + k57_co2(ig,:,:) * g_weight(iq)
     k57_grey_data(iCH4,:,:) = k57_grey_data(iCH4,:,:) + k57_ch4(ig,:,:) * g_weight(iq)
   enddo
   !=====  interval 58:   =====!
   k58_grey_data(:,:,:) = 0.0
   do ig=1, ngauss_pts(58)
     iq = iq + 1
     k58_grey_data(iH2O,:,:) = k58_grey_data(iH2O,:,:) + k58_h2o(ig,:,:) * g_weight(iq)
     k58_grey_data(iCO2,:,:) = k58_grey_data(iCO2,:,:) + k58_co2(ig,:,:) * g_weight(iq)
     k58_grey_data(iCH4,:,:) = k58_grey_data(iCH4,:,:) + k58_ch4(ig,:,:) * g_weight(iq)
   enddo
   !=====  interval 59:   =====!
   k59_grey_data(:,:,:) = 0.0
   do ig=1, ngauss_pts(59)
     iq = iq + 1
     k59_grey_data(iH2O,:,:) = k59_grey_data(iH2O,:,:) + k59_h2o(ig,:,:) * g_weight(iq)
     k59_grey_data(iCO2,:,:) = k59_grey_data(iCO2,:,:) + k59_co2(ig,:,:) * g_weight(iq)
     k59_grey_data(iCH4,:,:) = k59_grey_data(iCH4,:,:) + k59_ch4(ig,:,:) * g_weight(iq)
   enddo
   !=====  interval 60:   =====!
   k60_grey_data(:,:,:) = 0.0
   do ig=1, ngauss_pts(60)
     iq = iq + 1
     k60_grey_data(iH2O,:,:) = k60_grey_data(iH2O,:,:) + k60_h2o(ig,:,:) * g_weight(iq)
     k60_grey_data(iCO2,:,:) = k60_grey_data(iCO2,:,:) + k60_co2(ig,:,:) * g_weight(iq)
     k60_grey_data(iCH4,:,:) = k60_grey_data(iCH4,:,:) + k60_ch4(ig,:,:) * g_weight(iq)
   enddo
   !=====  interval 61:   =====!
   k61_grey_data(:,:,:) = 0.0
   do ig=1, ngauss_pts(61)
     iq = iq + 1
     k61_grey_data(iH2O,:,:) = k61_grey_data(iH2O,:,:) + k61_h2o(ig,:,:) * g_weight(iq)
     k61_grey_data(iCO2,:,:) = k61_grey_data(iCO2,:,:) + k61_co2(ig,:,:) * g_weight(iq)
     k61_grey_data(iCH4,:,:) = k61_grey_data(iCH4,:,:) + k61_ch4(ig,:,:) * g_weight(iq)
   enddo
   !=====  interval 62:   =====!
   k62_grey_data(:,:,:) = 0.0
   do ig=1, ngauss_pts(62)
     iq = iq + 1
     k62_grey_data(iH2O,:,:) = k62_grey_data(iH2O,:,:) + k62_h2o(ig,:,:) * g_weight(iq)
     k62_grey_data(iCO2,:,:) = k62_grey_data(iCO2,:,:) + k62_co2(ig,:,:) * g_weight(iq)
     k62_grey_data(iCH4,:,:) = k62_grey_data(iCH4,:,:) + k62_ch4(ig,:,:) * g_weight(iq)
   enddo
   !=====  interval 63:   =====!
   k63_grey_data(:,:,:) = 0.0
   do ig=1, ngauss_pts(63)
     iq = iq + 1
     k63_grey_data(iH2O,:,:) = k63_grey_data(iH2O,:,:) + k63_h2o(ig,:,:) * g_weight(iq)
     k63_grey_data(iCO2,:,:) = k63_grey_data(iCO2,:,:) + k63_co2(ig,:,:) * g_weight(iq)
     k63_grey_data(iCH4,:,:) = k63_grey_data(iCH4,:,:) + k63_ch4(ig,:,:) * g_weight(iq)
   enddo
   !=====  interval 64:   =====!
   k64_grey_data(:,:,:) = 0.0
   do ig=1, ngauss_pts(64)
     iq = iq + 1
     k64_grey_data(iH2O,:,:) = k64_grey_data(iH2O,:,:) + k64_h2o(ig,:,:) * g_weight(iq)
     k64_grey_data(iCO2,:,:) = k64_grey_data(iCO2,:,:) + k64_co2(ig,:,:) * g_weight(iq)
     k64_grey_data(iCH4,:,:) = k64_grey_data(iCH4,:,:) + k64_ch4(ig,:,:) * g_weight(iq)
   enddo
   !=====  interval 65:   =====!
   k65_grey_data(:,:,:) = 0.0
   do ig=1, ngauss_pts(65)
     iq = iq + 1
     k65_grey_data(iH2O,:,:) = k65_grey_data(iH2O,:,:) + k65_h2o(ig,:,:) * g_weight(iq)
     k65_grey_data(iCO2,:,:) = k65_grey_data(iCO2,:,:) + k65_co2(ig,:,:) * g_weight(iq)
     k65_grey_data(iCH4,:,:) = k65_grey_data(iCH4,:,:) + k65_ch4(ig,:,:) * g_weight(iq)
   enddo
   !=====  interval 66:   =====!
   k66_grey_data(:,:,:) = 0.0
   do ig=1, ngauss_pts(66)
     iq = iq + 1
     k66_grey_data(iH2O,:,:) = k66_grey_data(iH2O,:,:) + k66_h2o(ig,:,:) * g_weight(iq)
     k66_grey_data(iCO2,:,:) = k66_grey_data(iCO2,:,:) + k66_co2(ig,:,:) * g_weight(iq)
     k66_grey_data(iCH4,:,:) = k66_grey_data(iCH4,:,:) + k66_ch4(ig,:,:) * g_weight(iq)
   enddo
   !=====  interval 67:   =====!
   k67_grey_data(:,:,:) = 0.0
   do ig=1, ngauss_pts(67)
     iq = iq + 1
     k67_grey_data(iH2O,:,:) = k67_grey_data(iH2O,:,:) + k67_h2o(ig,:,:) * g_weight(iq)
     k67_grey_data(iCO2,:,:) = k67_grey_data(iCO2,:,:) + k67_co2(ig,:,:) * g_weight(iq)
     k67_grey_data(iCH4,:,:) = k67_grey_data(iCH4,:,:) + k67_ch4(ig,:,:) * g_weight(iq)
   enddo
   !=====  interval 68:   =====!
   k68_grey_data(:,:,:) = 0.0
   do ig=1, ngauss_pts(68)
     iq = iq + 1
     k68_grey_data(iH2O,:,:) = k68_grey_data(iH2O,:,:) + k68_h2o(ig,:,:) * g_weight(iq)
     k68_grey_data(iCO2,:,:) = k68_grey_data(iCO2,:,:) + k68_co2(ig,:,:) * g_weight(iq)
     k68_grey_data(iCH4,:,:) = k68_grey_data(iCH4,:,:) + k68_ch4(ig,:,:) * g_weight(iq)
   enddo

  end subroutine setup_gray_gas_matrix

end module exo_init_ref
