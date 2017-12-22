#include <params.h>
#include <misc.h>

module radiation

!version. 1D standalone RT general
!---------------------------------------------------------------------
! Purpose:
!
! Provides the CAM interface to the radiation code
!
! Adapts correlated K radiative transfer code for use with CAM/WACCM.
! Thie code uses the two stream radiative transfer method described in 
! Toon et al (1989).  Quadrature is used for shortwave, hemispheric mean
! is used for longwave.  Gas phase optical depths are calculate using a 
! correlated K-distribution method (Mlawer, 1997), with overlapping bands 
! treated via an amount weighted scheme (Shi et al, 2009). 
!
! Cloud optics treated using mie scattering for both liquid and ice clouds.
! Cloud overlap is treated using Monte Carlo Independent Column Approximation.
!
! Water vapor and CO2 continuum from MTCKD
!
!
! Revision history
! September 2010, E. T. Wolf, R. Urata
! March     2014, E. T. Wolf --- decoupled solar and IR streams
! December  2015  E. T. Wolf --- descope calc_gasopd for buildtime
!                                inclusion of multiple k-set options
!                            --- incorporated DIML, direct integration
!                                using mid-layer planck functions
!---------------------------------------------------------------------
!  
  use shr_kind_mod,     only: r8 => shr_kind_r8
  use shr_const_mod,    only: SHR_CONST_PI,SHR_CONST_PI, SHR_CONST_G, &
                              SHR_CONST_RGAS, SHR_CONST_AVOGAD, &
                              SHR_CONST_STEBOL, SHR_CONST_CSEC, &
                              SHR_CONST_MSDIST2, SHR_CONST_BOLTZ, &
                              SHR_CONST_RHOFW, SHR_CONST_RHOICE, &
 			      SHR_CONST_LOSCHMIDT
  use physconst,        only: mwn2, mwco2, mwch4, mwh2o, mwo2, mwh2, mwo3, mwdry, cpair
  use ppgrid            ! pver, pverp is here
  !use pmgrid            ! masterproc is here
  !use abortutils,       only: endrun
  !use rad_interp_mod    
  use radgrid
  use kabs
  use solar
  use calc_opd_mod
  !use phys_grid,        only: get_lat_p

  ! This defines the parameters and variables (in multiple named common blocks)
  ! that are provided by the CARMA code. CARMA is F77 code and uses common
  ! blocks instead of modules.
  !include 'aerad.h'

  implicit none
  private
  save

!------------------------------------------------------------------------
!
! Public interfaces
!  
  !public :: radiation_register 
  !public :: radiation_defaultopts
  !public :: radiation_setopts      
  public :: radiation_init        
  !public :: radiation_tend
  !public :: radiation_do
  public :: aerad_driver
  public :: init_ref



!------------------------------------------------------------------------
!
! private data
!
  
  ! Default values for namelist variables

  integer :: iradsw = 1      ! freq. of shortwave radiation calc in time steps (positive)
                             ! or hours (negative).
  integer :: iradlw = 1      ! frequency of longwave rad. calc. in time steps (positive)
                             ! or hours (negative).
  integer :: iradae = -12    ! frequency of absorp/emis calc in time steps (positive)
                             ! or hours (negative).
  integer :: irad_always = 0 ! Specifies length of time in timesteps (positive)
                             ! or hours (negative) SW/LW radiation will be
                             ! run continuously from the start of an
                             ! initial run

  ! At present this is used to scale the total strength of the input solar spectrum
  !! Kludge
  !! THIS OVERRIDES S0 from input file!!! 
  !real(r8) :: scon = 1360.0  ! Default Solar Constant, [w m-2]
  !real(r8) :: scon = 1371.734  ! Default Solar Constant, [w m-2]
  ! this is the value used in the 1d profiles
  real(r8) :: scon = 1360./2.  ! div 2 for day/night split for 1D column
  	      	     	       ! other factor of two from coszrs
  			       
  
  ! K coefficient file names (no more splitting sw/lw)
  !character(len=256) :: k01_file
  !character(len=256) :: k02_file
  !character(len=256) :: k03_file
  !character(len=256) :: k04_file
  !character(len=256) :: k05_file
  !character(len=256) :: k06_file
  !character(len=256) :: k07_file
  !character(len=256) :: k08_file
  !character(len=256) :: k09_file
  !character(len=256) :: k10_file
  !character(len=256) :: k11_file
  !character(len=256) :: k12_file
  !character(len=256) :: k13_file
  !character(len=256) :: k14_file 
  !character(len=256) :: k15_file 
  !character(len=256) :: k16_file
  !character(len=256) :: k17_file 
  !character(len=256) :: k18_file 
  !character(len=256) :: k19_file 
  !character(len=256) :: k20_file 
  !character(len=256) :: k21_file 
  !character(len=256) :: k22_file 
  !character(len=256) :: k23_file 
  !character(len=256) :: k24_file 
  !character(len=256) :: k25_file 
  !character(len=256) :: k26_file     
  !character(len=256) :: k27_file     
  !character(len=256) :: k28_file     

  ! K coefficients for continuum files
  !character(len=256) :: kh2oself_file
  !character(len=256) :: kco2cont_file
  !character(len=256) :: kh2h2cia_file
  !character(len=256) :: kh2n2cia_file
   
  ! Cloud optical constant files
  !character(len=256) :: cldoptsL_file ! name of water cloud optics file
  !character(len=256) :: cldoptsI_file ! name of ice cloud optics file

  integer :: openstatus

  ! Approximate smallest double precision floating point difference
  real(r8), parameter :: SMALLd = 1.0d-12    ! was 1.0d12
  real(r8), parameter :: SMALLe = 1.0e-12  ! was 1.0e12

  real(r8), parameter :: sqrt3 = 1.732050808d0	    ! square root of 3
  real(r8), parameter :: mb_to_atm = 9.869233e-4    ! convert pressure from Pa to atm

  ! Define constants for Gauss quadratrue calculations
  integer, parameter  :: ngangles = 3               ! # of Gauss angles to use


  ! Weights for Gaussian quadrature (for each probability interval) over a
  ! hemisphere subdivided into 'ngangles' equal angles  [none] (there are
  ! 'ngangles' of these):  
  ! Gauss weight vector and zenith angle
  real(r8), dimension(ntot_gpt) :: g_weight
  real(r8), dimension(ngangles) :: g_ang_weight 
  real(r8), dimension(ngangles) :: g_angle 

! For ngangles = 3:
  data g_angle / 0.21234054, 0.59053314, 0.91141204 / 
  data g_ang_weight / 0.06982698, 0.22924111, 0.20093191 /
  
  
  !
  ! CO2 continuum information, self and foreign
  !
  !!integer, parameter   :: ngCO2 = 8                           ! # of gauss points per spectral interval
  !!real(r8), parameter :: pref_co2 = 1013.0
  !!real(r8), parameter :: vref_co2 = 0.1 ! reference CO2 for continuum calculation
  !!real(r8), dimension(ngCO2,ntot_wavlnrng) :: kco2cont_8gpt  ! CO2 continuum data from file
  !!real(r8), dimension(ntot_gpt) :: kco2cont  ! CO2 continuum data, reduced to gauss pt grid
  !!!!real(r8), dimension(ntot_wavlnrng, ks_ntemp) :: kco2cont

  ! 8 to 16 gauss point mapping
  !!integer, dimension(16) :: map8to16gpt
  
  !!data map8to16gpt / 1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6, 7, 7, 8, 8 /

  ! CIA absorption data from HITRAN
  !!real(r8), dimension(ntot_wavlnrng,kh2h2_ntemp) :: kh2h2  ! H2-H2 CIA data [cm-1 amagat-2]
  !!real(r8), dimension(ntot_wavlnrng,kh2n2_ntemp) :: kh2n2  ! H2-N2 CIA data [cm-1 amagat-2]
  !!real(r8), dimension(ntot_wavlnrng,kn2n2_ntemp) :: kn2n2  ! N2-N2 CIA data [cm-1 amagat-2]
  
  !
  ! Planck function interpolation table temperatures, etc. [K]:
  !
  integer, parameter  :: tpft_inc = 10                      ! increment [# pts/K]
  integer, parameter  :: tpft_beg = 1                       ! start temp (1 K)
  integer, parameter  :: tpft_end = 1000*tpft_inc           ! end temp (1000 K)
  integer, parameter  :: tpft_nt = tpft_end                 ! table temp dimension
  real(r8), parameter :: tpft_invfinc = 1.d0/dble(tpft_inc) ! factor dble(1/tpft_inc)
  real(r8), parameter :: tpft_dinc = dble(tpft_inc)         ! factor dble(tpft_inc)  
  real(r8), parameter :: tpft_finc = real(tpft_inc)         ! factor real(tpft_inc)
  real(r8), dimension(tpft_nt,ntot_gpt) :: ptemp_itp        ! table [tpft_nt,ntot_gpt]

  ! data structures containing cloud optics data from file
  !!real(r8) :: Qcldliq(nrel, ntot_wavlnrng)
  !!real(r8) :: Wcldliq(nrel, ntot_wavlnrng)
  !!real(r8) :: Gcldliq(nrel, ntot_wavlnrng)

  !!real(r8) :: Qcldice(nrei, ntot_wavlnrng)
  !!real(r8) :: Wcldice(nrei, ntot_wavlnrng)
  !!real(r8) :: Gcldice(nrei, ntot_wavlnrng)
  
  !! top level indices for longwave and cloud calculations
  !!NOTES: these currently not used
  !!integer :: camtop  ! top level to solve using CAM CKD RT.  
  !!integer :: ntopcld ! top level to solve for cloud overlap  
  !!integer :: nlevsRT ! number of levels to calculate RT

  !------------------------------------------------------------------------------
  ! Radiative transfer model variable/array declarations

  ! defined in INIT_REF
  integer  :: lw_iwbeg
  integer  :: lw_iwend 
  integer  :: sw_iwbeg 
  integer  :: sw_iwend 
  integer  :: lw_ipbeg 
  integer  :: lw_ipend 
  integer  :: sw_ipbeg 
  integer  :: sw_ipend

  ! defined in INIT_REF  
  real(r8) :: U1Ssol      ! 2*PI / mu1 factors 
  real(r8) :: U1Sir         ! 2*PI / mu1 factors 
  real(r8) :: U1Isol      ! mu1 factors
  real(r8) :: U1Iir         ! mu1 factors
  real(r8) :: U1I2sol    ! mu1 factors
  real(r8) :: U1I2ir       ! mu1 factors
  real(r8), dimension(ntot_gpt) :: gw_solflux
  real(r8), dimension(ntot_wavlnrng) :: solflux


  !NOTES: THESE AREN"T HOOKED UP YET
  ! Aerosol optical properties
  !real(r8), dimension(naerspc_rad,ntot_wavlnrng,pverp) ::  singscat_aer 
  !real(r8), dimension(naerspc_rad,ntot_wavlnrng,pverp) ::  asym_aer     
  !real(r8), dimension(naerspc_rad,ntot_wavlnrng,pverp) ::  tau_aer      
 
  ! Solar gauss weight corrections for ozone UV absorption region
!!  real(r8), dimension(8) :: solar_gweight_sp27
!!  real(r8), dimension(8) :: solar_gweight_sp28

!!  data solar_gweight_sp27 / 0.5048931837, 0.3053007424, 0.1143089905, 0.0489797853, &
!!                            0.0219137613, 0.0025616873, 0.0018370311, 0.0002052875 /
!!  data solar_gweight_sp28 / 0.1023708880, 0.2846772969, 0.2338085473, 0.2813689709, &
!!                            0.0869001597, 0.0059007024, 0.0037554756, 0.0012206517 /


!============================================================================
contains
!============================================================================

!============================================================================
!
! Public subroutines 
!
!============================================================================

!  subroutine radiation_register
!-----------------------------------------------------------------------
!
! Register radiation fields in the physics buffer
!
!-----------------------------------------------------------------------
!
!    !use phys_buffer,  only: pbuf_times, pbuf_add
!
!    integer idx
!    !call pbuf_add('QRS' , 'global', 1,pver,1, idx) ! shortwave radiative heating rate
!    !call pbuf_add('QRL' , 'global', 1,pver,1, idx) ! longwave  radiative heating rate
!
!  end subroutine radiation_register


!============================================================================

!  subroutine radiation_defaultopts(iradsw_out, iradlw_out, iradae_out, irad_always_out, &
!                                   k01_file_out, k02_file_out, k03_file_out, k04_file_out, &   
!                                   k05_file_out, k06_file_out, k07_file_out, k08_file_out, &   
!                                   k09_file_out, k10_file_out, k11_file_out, k12_file_out, &   
!                                   k13_file_out, k14_file_out, k15_file_out, k16_file_out, &
!                                   k17_file_out, k18_file_out, k19_file_out, k20_file_out, &
!                                   k21_file_out, k22_file_out, k23_file_out, k24_file_out, &
!                                   k25_file_out, k27_file_out, k28_file_out, &
!                                   kh2oself_file_out, kco2cont_file_out, &
!                                   cldoptsL_file_out, cldoptsI_file_out, scon_out )
!-----------------------------------------------------------------------
!
!  Purpose: Return default runtime options
!  NOTES: iradsw, iradlw are not hooked up.  iradae can be completely eliminated.
!
!------------------------------------------------------------------------
!
! Input Arguments
!
!    integer, intent(out), optional :: iradsw_out
!    integer, intent(out), optional :: iradlw_out
!    integer, intent(out), optional :: iradae_out
!    integer, intent(out), optional :: irad_always_out
!    character(len=*), intent(out), optional :: k01_file_out
!    character(len=*), intent(out), optional :: k02_file_out
!    character(len=*), intent(out), optional :: k03_file_out
!    character(len=*), intent(out), optional :: k04_file_out
!    character(len=*), intent(out), optional :: k05_file_out
!    character(len=*), intent(out), optional :: k06_file_out
!    character(len=*), intent(out), optional :: k07_file_out
!    character(len=*), intent(out), optional :: k08_file_out
!    character(len=*), intent(out), optional :: k09_file_out
!    character(len=*), intent(out), optional :: k10_file_out
!    character(len=*), intent(out), optional :: k11_file_out
!    character(len=*), intent(out), optional :: k12_file_out
!    character(len=*), intent(out), optional :: k13_file_out
!    character(len=*), intent(out), optional :: k14_file_out
!    character(len=*), intent(out), optional :: k15_file_out
!    character(len=*), intent(out), optional :: k16_file_out
!    character(len=*), intent(out), optional :: k17_file_out
!    character(len=*), intent(out), optional :: k18_file_out
!    character(len=*), intent(out), optional :: k19_file_out
!    character(len=*), intent(out), optional :: k20_file_out
!    character(len=*), intent(out), optional :: k21_file_out
!    character(len=*), intent(out), optional :: k22_file_out
!    character(len=*), intent(out), optional :: k23_file_out
!    character(len=*), intent(out), optional :: k24_file_out
!    character(len=*), intent(out), optional :: k25_file_out
!    character(len=*), intent(out), optional :: k27_file_out
!    character(len=*), intent(out), optional :: k28_file_out
!    character(len=*), intent(out), optional :: kh2oself_file_out
!    character(len=*), intent(out), optional :: kco2cont_file_out
!    character(len=*), intent(out), optional :: cldoptsL_file_out
!    character(len=*), intent(out), optional :: cldoptsI_file_out
!    real(r8), intent(out), optional :: scon_out

!------------------------------------------------------------------------
!
! Start Code
!
!   if ( present(iradsw_out) )      iradsw_out = iradsw
!   if ( present(iradlw_out) )      iradlw_out = iradlw
!   if ( present(iradae_out) )      iradae_out = iradae
!   if ( present(irad_always_out) ) irad_always_out = irad_always
!   if ( present(k01_file_out) ) k01_file_out = k01_file
!   if ( present(k02_file_out) ) k02_file_out = k02_file
!   if ( present(k03_file_out) ) k03_file_out = k03_file
!   if ( present(k04_file_out) ) k04_file_out = k04_file
!   if ( present(k05_file_out) ) k05_file_out = k05_file
!   if ( present(k06_file_out) ) k06_file_out = k06_file
!   if ( present(k07_file_out) ) k07_file_out = k07_file
!   if ( present(k08_file_out) ) k08_file_out = k08_file
!   if ( present(k09_file_out) ) k09_file_out = k09_file
!   if ( present(k10_file_out) ) k10_file_out = k10_file
!   if ( present(k11_file_out) ) k11_file_out = k11_file
!   if ( present(k12_file_out) ) k12_file_out = k12_file
!   if ( present(k13_file_out) ) k13_file_out = k13_file
!   if ( present(k14_file_out) ) k14_file_out = k14_file
!   if ( present(k15_file_out) ) k15_file_out = k15_file
!   if ( present(k16_file_out) ) k16_file_out = k16_file
!   if ( present(k17_file_out) ) k17_file_out = k17_file
!   if ( present(k18_file_out) ) k18_file_out = k18_file
!   if ( present(k19_file_out) ) k19_file_out = k19_file
!   if ( present(k20_file_out) ) k20_file_out = k20_file
!   if ( present(k21_file_out) ) k21_file_out = k21_file
!   if ( present(k22_file_out) ) k22_file_out = k22_file
!   if ( present(k23_file_out) ) k23_file_out = k23_file
!   if ( present(k24_file_out) ) k24_file_out = k24_file
!   if ( present(k25_file_out) ) k25_file_out = k25_file
!   if ( present(k27_file_out) ) k27_file_out = k27_file
!   if ( present(k28_file_out) ) k28_file_out = k28_file
!   if ( present(kh2oself_file_out) ) kh2oself_file_out = kh2oself_file
!   if ( present(kco2cont_file_out) ) kco2cont_file_out = kco2cont_file
!   if ( present(cldoptsL_file_out) ) cldoptsL_file_out = cldoptsL_file
!   if ( present(cldoptsI_file_out) ) cldoptsI_file_out = cldoptsI_file
!   if ( present(scon_out) ) scon_out = scon
!
!  end subroutine radiation_defaultopts


!============================================================================

!  subroutine radiation_setopts (dtime, nhtfrq, iradsw_in, iradlw_in, iradae_in, irad_always_in, &
!                                k01_file_in, k02_file_in, k03_file_in, k04_file_in, &
!                                k05_file_in, k06_file_in, k07_file_in, k08_file_in, &
!                                k09_file_in, k10_file_in, k11_file_in, k12_file_in, &
!                                k13_file_in, k14_file_in, k15_file_in, k16_file_in, &
!                                k17_file_in, k18_file_in, k19_file_in, k20_file_in, &
!                                k21_file_in, k22_file_in, k23_file_in, k24_file_in, &
!                                k25_file_in, k27_file_in, k28_file_in, &
!                                kh2oself_file_in, kco2cont_file_in, &
!                                cldoptsL_file_in, cldoptsI_file_in, scon_in )                                
!-----------------------------------------------------------------------
! Purpose: Set runtime options
!*** NOTE *** This routine needs information about dtime (init by dycore)
!              and nhtfrq (init by history) to do its checking.  Being called
!              from runtime_opts provides these values possibly before they
!              have been set in the modules responsible for them.
!-----------------------------------------------------------------------
!
! Input Arguments
! 
!    integer, intent(in)           :: dtime           ! timestep size (s)
!    integer, intent(in)           :: nhtfrq          ! output frequency of primary history file
!    integer, intent(in), optional :: iradsw_in
!    integer, intent(in), optional :: iradlw_in
!    integer, intent(in), optional :: iradae_in
!    integer, intent(in), optional :: irad_always_in
!    character(len=*), intent(in), optional :: k01_file_in
!    character(len=*), intent(in), optional :: k02_file_in
!    character(len=*), intent(in), optional :: k03_file_in
!    character(len=*), intent(in), optional :: k04_file_in
!    character(len=*), intent(in), optional :: k05_file_in
!    character(len=*), intent(in), optional :: k06_file_in
!    character(len=*), intent(in), optional :: k07_file_in
!    character(len=*), intent(in), optional :: k08_file_in
!    character(len=*), intent(in), optional :: k09_file_in
!    character(len=*), intent(in), optional :: k10_file_in
!    character(len=*), intent(in), optional :: k11_file_in
!    character(len=*), intent(in), optional :: k12_file_in
!    character(len=*), intent(in), optional :: k13_file_in
!    character(len=*), intent(in), optional :: k14_file_in
!    character(len=*), intent(in), optional :: k15_file_in
!    character(len=*), intent(in), optional :: k16_file_in
!    character(len=*), intent(in), optional :: k17_file_in
!    character(len=*), intent(in), optional :: k18_file_in
!    character(len=*), intent(in), optional :: k19_file_in
!    character(len=*), intent(in), optional :: k20_file_in
!    character(len=*), intent(in), optional :: k21_file_in
!    character(len=*), intent(in), optional :: k22_file_in
!    character(len=*), intent(in), optional :: k23_file_in
!    character(len=*), intent(in), optional :: k24_file_in
!    character(len=*), intent(in), optional :: k25_file_in
!    character(len=*), intent(in), optional :: k27_file_in
!    character(len=*), intent(in), optional :: k28_file_in
!    character(len=*), intent(in), optional :: kh2oself_file_in
!    character(len=*), intent(in), optional :: kco2cont_file_in
!    character(len=*), intent(in), optional :: cldoptsL_file_in
!    character(len=*), intent(in), optional :: cldoptsI_file_in
!    real(r8), intent(in), optional :: scon_in
!
!-----------------------------------------------------------------------
!
! Local Variables
!
!    integer :: ntspdy    ! number of timesteps per day
!    integer :: nhtfrq1   ! local copt of input arg nhtfrq
!-----------------------------------------------------------------------
!
! Start Code
!
!   if ( present(iradsw_in) )      iradsw = iradsw_in
!   if ( present(iradlw_in) )      iradlw = iradlw_in
!   if ( present(iradae_in) )      iradae = iradae_in
!   if ( present(irad_always_in) ) irad_always = irad_always_in
!   if ( present(k01_file_in) ) k01_file = k01_file_in
!   if ( present(k01_file_in) ) k01_file = k01_file_in
!   if ( present(k02_file_in) ) k02_file = k02_file_in
!   if ( present(k03_file_in) ) k03_file = k03_file_in
!   if ( present(k04_file_in) ) k04_file = k04_file_in
!   if ( present(k05_file_in) ) k05_file = k05_file_in
!   if ( present(k06_file_in) ) k06_file = k06_file_in
!   if ( present(k07_file_in) ) k07_file = k07_file_in
!   if ( present(k08_file_in) ) k08_file = k08_file_in
!   if ( present(k09_file_in) ) k09_file = k09_file_in
!   if ( present(k10_file_in) ) k10_file = k10_file_in
!   if ( present(k11_file_in) ) k11_file = k11_file_in
!   if ( present(k12_file_in) ) k12_file = k12_file_in
!   if ( present(k13_file_in) ) k13_file = k13_file_in
!   if ( present(k14_file_in) ) k14_file = k14_file_in
!   if ( present(k15_file_in) ) k15_file = k15_file_in
!   if ( present(k16_file_in) ) k16_file = k16_file_in
!   if ( present(k17_file_in) ) k17_file = k17_file_in
!   if ( present(k18_file_in) ) k18_file = k18_file_in
!   if ( present(k19_file_in) ) k19_file = k19_file_in
!   if ( present(k20_file_in) ) k20_file = k20_file_in
!   if ( present(k21_file_in) ) k21_file = k21_file_in
!   if ( present(k22_file_in) ) k22_file = k22_file_in
!   if ( present(k23_file_in) ) k23_file = k23_file_in
!   if ( present(k24_file_in) ) k24_file = k24_file_in
!   if ( present(k25_file_in) ) k25_file = k25_file_in
!   if ( present(k27_file_in) ) k27_file = k27_file_in
!   if ( present(k28_file_in) ) k28_file = k28_file_in
!   if ( present(kh2oself_file_in) ) kh2oself_file = kh2oself_file_in
!   if ( present(kco2cont_file_in) ) kco2cont_file = kco2cont_file_in
!   if ( present(cldoptsL_file_in) ) cldoptsL_file = cldoptsL_file_in
!   if ( present(cldoptsI_file_in) ) cldoptsI_file = cldoptsI_file_in
!   if ( present(scon_in) ) scon = scon_in

!   ! Convert iradsw, iradlw and irad_always from hours to timesteps if necessary
!   if (iradsw      < 0) iradsw      = nint((-iradsw     *3600.)/dtime)
!   if (iradlw      < 0) iradlw      = nint((-iradlw     *3600.)/dtime)
!   if (irad_always < 0) irad_always = nint((-irad_always*3600.)/dtime)

   ! Convert iradae from hours to timesteps if necessary and check that
   ! iradae must be an even multiple of iradlw
!   if (iradae < 0) iradae = nint((-iradae*3600.)/dtime)
!   if (mod(iradae,iradlw)/=0) then
!      write(6,*)'radiation_setopts: iradae must be an even multiple of iradlw.'
!      write(6,*)'     iradae = ',iradae,', iradlw = ',iradlw
!      call endrun('radiation_setopts: iradae must be an even multiple of iradlw.')
!   end if

   ! Do absorptivities/emissivities have to go on a restart dataset?
   ! The value of nhtfrq from the namelist may need conversion.
!   nhtfrq1 = nhtfrq
!   if (nhtfrq1 < 0) nhtfrq1 = nint((-nhtfrq1*3600.)/dtime)
!   ntspdy = nint(SHR_CONST_CSEC/dtime)
!   if (nhtfrq1 /= 0) then
!      if (masterproc .and. mod(nhtfrq1,iradae)/=0) then
!         write(6,*)'radiation_setopts: *** NOTE: Extra overhead invoked putting',  &
!            ' a/e numbers on restart dataset. ***   ',         &
!            ' To avoid, make mod(nhtfrq,iradae) = 0'
!      end if
!   else
!      if (masterproc) then
!         if (mod(ntspdy,iradae) /= 0 .or. iradae > ntspdy) then
!            write(6,*)'radiation_setopts: *** NOTE: Extra overhead invoked',  &
!               ' putting a/e numbers on restart dataset. ***'
!            write(6,*)' To avoid, make mod(timesteps per day,iradae)= 0'
!         end if
!      end if
!   end if

!  end subroutine radiation_setopts


!============================================================================

!  function radiation_do(op)
!-----------------------------------------------------------------------
!
! Purpose: Returns true if the specified operation is done this timestep.
!
!-----------------------------------------------------------------------
!    use time_manager,   only: get_nstep
!------------------------------------------------------------------------
!
! Arguments
!
!
!    character(len=*), intent(in) :: op             ! name of operation
!    logical                      :: radiation_do   ! return value

!------------------------------------------------------------------------
!
! Local Variables
!

!    integer :: nstep             ! current timestep number

!------------------------------------------------------------------------
!
! Start Code
!
! NOTES: FIX THIS SECTION
!    
!    nstep = get_nstep()
!
!    select case (op)
!
!    case ('sw') ! do a shortwave heating calc this timestep?
!      radiation_do = nstep == 0  .or.  iradsw == 1                     &
!                    .or. (mod(nstep-1,iradsw) == 0  .and.  nstep /= 1) &
!                    .or. nstep <= irad_always
!
!    case ('lw') ! do a longwave heating calc this timestep?
!      radiation_do = nstep == 0  .or.  iradlw == 1                     &
!                    .or. (mod(nstep-1,iradlw) == 0  .and.  nstep /= 1) &
!                    .or. nstep <= irad_always
!
!    case ('absems') ! do an absorptivity/emissivity calculation this timestep?
!      ! no abs/ems calc in correlated-k
!      radiation_do = .false.
!
!
!    case ('aeres') ! write absorptivity/emissivity to restart file this timestep?
!      ! no abs/ems calc in correlated-k 
!      radiation_do = .false. 
!
!    case default
!     call endrun('radiation_do: unknown operation:'//op)
!
!    end select
!
!  end function radiation_do


!============================================================================

  subroutine radiation_printopts(unit)
!-----------------------------------------------------------------------
!
! Purpose: Print runtime options to log.
!
!-----------------------------------------------------------------------
!NOTES: change this subroutine to print options relevent to my correlated-K
!model

!------------------------------------------------------------------------
!
! Arguments
! 
    integer, intent(in) :: unit  ! Fortran unit for log output

!------------------------------------------------------------------------
!
! Start Code
! 

    if(irad_always /= 0) write(unit,10) irad_always
   write(unit,20) iradsw,iradlw,iradae
10 format(' Execute SW/LW radiation continuously until timestep ',i5)
20 format(' Frequency of Shortwave Radiation calc. (IRADSW)     ',i5/, &
          ' Frequency of Longwave Radiation calc. (IRADLW)      ',i5/,  &
          ' Frequency of Absorptivity/Emissivity calc. (IRADAE) ',i5)

  end subroutine radiation_printopts


!============================================================================

  subroutine radiation_get(iradsw_out, iradlw_out, iradae_out, irad_always_out)
!-----------------------------------------------------------------------
!
! Purpose: Provide access to private module data.  (This should be eliminated.)
!
!-----------------------------------------------------------------------

!------------------------------------------------------------------------
!
! Arguments
!
    integer, intent(out), optional :: iradsw_out
    integer, intent(out), optional :: iradlw_out
    integer, intent(out), optional :: iradae_out
    integer, intent(out), optional :: irad_always_out

!------------------------------------------------------------------------
!
! Start Code
!

    if ( present(iradsw_out) )      iradsw_out = iradsw
    if ( present(iradlw_out) )      iradlw_out = iradlw
    if ( present(iradae_out) )      iradae_out = iradae
    if ( present(irad_always_out) ) irad_always_out = irad_always

  end subroutine radiation_get


!============================================================================

  subroutine radiation_init()
!-----------------------------------------------------------------------
!
! Purpose: Initialize the radiation parameterization, add fields to the 
!          history buffer
!
!-----------------------------------------------------------------------

!    use history,     only: addfld, add_default, phys_decomp

!#if ( defined WACCM_MOZART )
!      use solvar_interface, only: solvar_init
!#endif

!------------------------------------------------------------------------
!
! Local Variables
!
    integer :: iv
    integer :: ig
    
!------------------------------------------------------------------------
!
! Start Code
!

    ! Fill some general reference arrays:
    call init_ref

    ! Add Shortwave radiation fields
    !call addfld ('SOLIN   ','W/m2    ',1,    'A','Solar insolation',phys_decomp, sampling_seq='rad_lwsw')
    !call addfld ('SOLL    ','W/m2    ',1,    'A','Solar downward near infrared direct  to surface',phys_decomp, &
    !                                                                                       sampling_seq='rad_lwsw')
    !call addfld ('SOLS    ','W/m2    ',1,    'A','Solar downward visible direct  to surface',phys_decomp, &
    !                                                                                         sampling_seq='rad_lwsw')
    !call addfld ('SOLLD   ','W/m2    ',1,    'A','Solar downward near infrared diffuse to surface',phys_decomp, &
    !                                                                                        sampling_seq='rad_lwsw')
    !call addfld ('SOLSD   ','W/m2    ',1,    'A','Solar downward visible diffuse to surface',phys_decomp, sampling_seq='rad_lwsw')
    !call addfld ('QRS     ','K/day     ',pver, 'A','Solar heating rate',phys_decomp, sampling_seq='rad_lwsw')
    !call addfld ('FSNS    ','W/m2    ',1,    'A','Net solar flux at surface',phys_decomp, sampling_seq='rad_lwsw')
    !call addfld ('FSNT    ','W/m2    ',1,    'A','Net solar flux at top of model',phys_decomp, sampling_seq='rad_lwsw')
    !call addfld ('FSNTOA  ','W/m2    ',1,    'A','Net solar flux at top of atmosphere',phys_decomp, sampling_seq='rad_lwsw')
    !call addfld ('FSNTOAC ','W/m2    ',1,    'A','Clearsky net solar flux at top of atmosphere',phys_decomp, &
    !                                                                                          sampling_seq='rad_lwsw')
    !call addfld ('FSN200  ','W/m2    ',1,    'A','Net shortwave flux at 200 mb',phys_decomp, sampling_seq='rad_lwsw')
    !call addfld ('FSN200C ','W/m2    ',1,    'A','Clearsky net shortwave flux at 200 mb',phys_decomp, sampling_seq='rad_lwsw')
    !call addfld ('FSNTC   ','W/m2    ',1,    'A','Clearsky net solar flux at top of model',phys_decomp, sampling_seq='rad_lwsw')
    !call addfld ('FSNSC   ','W/m2    ',1,    'A','Clearsky net solar flux at surface',phys_decomp, sampling_seq='rad_lwsw')
    !call addfld ('FSDSC   ','W/m2    ',1,    'A','Clearsky downwelling solar flux at surface',phys_decomp, &
    !                                                                                               sampling_seq='rad_lwsw')
    !call addfld ('FSDS    ','W/m2    ',1,    'A','Downwelling solar flux at surface',phys_decomp, sampling_seq='rad_lwsw')
    !call addfld ('FUS     ','W/m2    ',pverp,'I','Shortwave upward flux',phys_decomp)
    !call addfld ('FDS     ','W/m2    ',pverp,'I','Shortwave downward flux',phys_decomp)
    !call addfld ('FUSC    ','W/m2    ',pverp,'I','Shortwave clear-sky upward flux',phys_decomp)
    !call addfld ('FDSC    ','W/m2    ',pverp,'I','Shortwave clear-sky downward flux',phys_decomp)
    !call addfld ('FSNIRTOA','W/m2    ',1,    'A','Net near-infrared flux (Nimbus-7 WFOV) at top of atmosphere', &
    !                                                                           phys_decomp, sampling_seq='rad_lwsw')
    !call addfld ('FSNRTOAC','W/m2    ',1,    'A','Clearsky net near-infrared flux (Nimbus-7 WFOV) at top of atmosphere', &
    !                                                                             phys_decomp, sampling_seq='rad_lwsw')
    !call addfld ('FSNRTOAS','W/m2    ',1,    'A','Net near-infrared flux (>= 0.7 microns) at top of atmosphere',phys_decomp, &
    !                                                                                          sampling_seq='rad_lwsw')
    !call addfld ('SWCF    ','W/m2    ',1,    'A','Shortwave cloud forcing',phys_decomp, sampling_seq='rad_lwsw')
    !
    !call add_default ('QRS     ', 1, ' ')
    !call add_default ('FSNS    ', 1, ' ')
    !call add_default ('FSNT    ', 1, ' ')
    !call add_default ('FSNTOA  ', 1, ' ')
    !call add_default ('FSNTOAC ', 1, ' ')
    !call add_default ('FSNTC   ', 1, ' ')
    !call add_default ('FSNSC   ', 1, ' ')
    !call add_default ('FSDSC   ', 1, ' ')
    !call add_default ('FSDS    ', 1, ' ')
    !call add_default ('SWCF    ', 1, ' ')

    ! aerosol forcing-only calculations
    !call addfld ('FSNT_RF ','W/m^2   ',1, 'I','Total column absorbed solar flux (radforce)' ,phys_decomp)
    !call addfld ('FSNTC_RF','W/m^2   ',1, 'I','Clear sky total column absorbed solar flux (radforce)' ,phys_decomp)
    !call addfld ('FSNS_RF ','W/m^2   ',1, 'I','Surface absorbed solar flux (radforce)' ,phys_decomp)
    !call addfld ('FSNSC_RF','W/m^2   ',1, 'I','Clear sky surface absorbed solar flux (radforce)' ,phys_decomp)
    !call addfld ('QRS_RF  ','K/s     ',pver, 'I','Solar heating rate (radforce)' ,phys_decomp)

    ! Longwave radiation
    !call addfld ('QRL     ','K/day     ',pver, 'A','Longwave heating rate',phys_decomp, sampling_seq='rad_lwsw')
    !call addfld ('FLNS    ','W/m2    ',1,    'A','Net longwave flux at surface',phys_decomp, sampling_seq='rad_lwsw')
    !call addfld ('FLNT    ','W/m2    ',1,    'A','Net longwave flux at top of model',phys_decomp, sampling_seq='rad_lwsw')
    !call addfld ('FLUT    ','W/m2    ',1,    'A','Upwelling longwave flux at top of model',phys_decomp, sampling_seq='rad_lwsw')
    !call addfld ('FLUTC   ','W/m2    ',1,    'A','Clearsky upwelling longwave flux at top of model',phys_decomp, &
    !                                                                                               sampling_seq='rad_lwsw')
    !call addfld ('FLNTC   ','W/m2    ',1,    'A','Clearsky net longwave flux at top of model',phys_decomp, sampling_seq='rad_lwsw')
    !call addfld ('FLN200  ','W/m2    ',1,    'A','Net longwave flux at 200 mb',phys_decomp, sampling_seq='rad_lwsw')
    !call addfld ('FLN200C ','W/m2    ',1,    'A','Clearsky net longwave flux at 200 mb',phys_decomp, sampling_seq='rad_lwsw')
    !call addfld ('FLNSC   ','W/m2    ',1,    'A','Clearsky net longwave flux at surface',phys_decomp, sampling_seq='rad_lwsw')
    !call addfld ('LWCF    ','W/m2    ',1,    'A','Longwave cloud forcing',phys_decomp, sampling_seq='rad_lwsw')
    !call addfld ('FUL     ','W/m2    ',pverp,'I','Longwave upward flux',phys_decomp)
    !call addfld ('FDL     ','W/m2    ',pverp,'I','Longwave downward flux',phys_decomp)
    !call addfld ('FULC    ','W/m2    ',pverp,'I','Longwave clear-sky upward flux',phys_decomp)
    !call addfld ('FDLC    ','W/m2    ',pverp,'I','Longwave clear-sky downward flux',phys_decomp)
    !call add_default ('QRL     ', 1, ' ')
    !call add_default ('FLNS    ', 1, ' ')
    !call add_default ('FLNT    ', 1, ' ')
    !call add_default ('FLUT    ', 1, ' ')
    !call add_default ('FLUTC   ', 1, ' ')
    !call add_default ('FLNTC   ', 1, ' ')
    !call add_default ('FLNSC   ', 1, ' ')
    !call add_default ('LWCF    ', 1, ' ')

    ! Heating rate needed for d(theta)/dt computation
    !call addfld ('HR      ','K/s     ',pver, 'A','Heating rate needed for d(theta)/dt computation',phys_decomp)     

  end subroutine radiation_init


!============================================================================
!
!  subroutine radiation_tend(state, ptend, pbuf, surface_state2d, &
!                            srfflx_state2d, landfrac, landm, icefrac, snowh,&
!                            fsns, fsnt, flns, flnt, fsds, concld, net_flx)
!-----------------------------------------------------------------------
! !!! THIS SUBROUTINE ISN"T CALLED IN STANDALONE COLUMN VERSION !!!
! Purpose: Driver for correlated K radiation computation.  Uses delta eddington 
!          two stream at short wavelengths and hemispheric mean two sream 
!          at IR wavelengths.
!
! Revision history:
! 2009?     R. Urata - adapted code from MRAMS/CARMA Mars model
! September 2010: E.T.Wolf     
!
!-----------------------------------------------------------------------
! 
!    !use phys_buffer,      only: pbuf_size_max, pbuf_fld, pbuf_old_tim_idx, pbuf_get_fld_idx
!    !use phys_grid,        only: get_rlat_all_p, get_rlon_all_p
!    !use param_cldoptics,  only: param_cldoptics_calc
!    !use physics_types,    only: physics_state, physics_ptend
!    !use time_manager,     only: get_curr_calday,get_nstep, get_curr_calday_rotation
!    !use comsrf,           only: surface_state, srfflx_state
!    !use history,          only: outfld
!    use radheat,          only: radheat_tend
!    use pmgrid,           only: plev, plevp
!    use pspect
!    use rad_constituents, only: rad_constituents_get
!    use carma,		  only: carma_is_active
!    use shr_orb_mod
!
!    implicit none
!
!!#if ( defined SCAM )
!      #include <max.h>
!      use scamMod, only: switch,have_cld,cldobs,have_clwp,clwpobs,have_tg,tground
!#endif
!                  
!#include <comctl.h>
!#include <comsol.h>
!#include <comhyb.h>
!
!
!------------------------------------------------------------------------
!
! Input Arguments
!
!    real(r8), intent(in), dimension(pcols) :: landfrac         ! land fraction
!    real(r8), intent(in), dimension(pcols) :: landm            ! land fraction ramp
!    real(r8), intent(in), dimension(pcols) :: icefrac          ! ice fraction
!    real(r8), intent(in), dimension(pcols) :: snowh            ! Snow depth (liquid water equivalent)
!    real(r8), intent(inout), dimension(pcols) :: fsns          ! Net solar flux at surface 
!    real(r8), intent(inout), dimension(pcols) :: fsnt          ! Net column abs solar flux at model top
!    real(r8), intent(inout), dimension(pcols) :: flns          ! Srf longwave cooling (up-down) flux
!    real(r8), intent(inout), dimension(pcols) :: flnt          ! Net outgoing lw flux at model top
!    real(r8), intent(out), dimension(pcols) :: fsds            ! Surface solar down flux
!    real(r8), intent(in), dimension(pcols,pver) :: concld     ! should be pbuf
!    real(r8), intent(inout), dimension(pcols) :: net_flx
!
!    !type(physics_state), intent(in), target :: state
!    !type(physics_ptend), intent(out)        :: ptend
!    !type(pbuf_fld),      intent(inout), dimension(pbuf_size_max) :: pbuf
!    type(surface_state), intent(inout)      :: surface_state2d
!    type(srfflx_state),  intent(in)         :: srfflx_state2d

!------------------------------------------------------------------------
!
!  Local Variables
!
!
!    integer, dimension(pcols) :: nmxrgn            ! Number of maximally overlapped regions
!    real(r8), dimension(pcols,pverp) :: pmxrgn     ! Maximum values of pressure for each
!                                                        ! maximally overlapped region.
!                                                        ! 0->pmxrgn(i,1) is range of pressure for
!                                                        !    1st region,pmxrgn(i,1)->pmxrgn(i,2) for
!                                                        !    2nd region, etc
!    real(r8), dimension(pcols,pver) :: cldemis     ! Cloud longwave emissivity
!    real(r8), dimension(pcols,pver) :: cicewp      ! in-cloud cloud ice water path (from param_cldoptics_calc)
!    real(r8), dimension(pcols,pver) :: cliqwp      ! in-cloud cloud liquid water path (from param_cldoptics_calc)
!    real(r8), dimension(pcols) ::  cltot           ! Diagnostic total cloud cover
!    real(r8), dimension(pcols) ::  cllow           !       "     low  cloud cover
!    real(r8), dimension(pcols) ::  clmed           !       "     mid  cloud cover
!    real(r8), dimension(pcols) ::  clhgh           !       "     hgh  cloud cover
!    real(r8), dimension(pcols,pver) :: ftem        ! Temporary workspace for outfld variables
!    real(r8), dimension(pcols,pver) :: ftem2	   ! Second temp workspace
!
!    integer itim, ifld
!    real(r8), pointer, dimension(:,:) :: rel     ! liquid effective drop radius (microns)
!    real(r8), pointer, dimension(:,:) :: rei     ! ice effective drop size (microns)
!    real(r8), pointer, dimension(:,:) :: cfrc  ! cloud fraction
!    real(r8), pointer, dimension(:,:) :: qrs     ! shortwave radiative heating rate 
!    real(r8), pointer, dimension(:,:) :: qrl     ! longwave  radiative heating rate 
!
!    integer lchnk, ncol
!    !real(r8) :: calday                        ! current calendar day
!    real(r8), dimension(pcols) :: clat        ! current latitudes(radians)
!    real(r8), dimension(pcols) :: clon        ! current longitudes(radians)
!    real(r8), dimension(pcols) ::  coszrs     ! Cosine solar zenith angle
!    !integer  :: nstep                         ! current timestep number
!    logical  :: conserve_energy = .true.      ! flag to carry (QRS,QRL)*dp across time steps
!
!    !
!    ! Local variables from radctl
!    !

 !   integer naer_groups       ! Num of aerosol groups for optical diagnostics
 !   parameter ( naer_groups = 7 )    ! current groupings are sul, sslt, all carbons, all dust, background, and all aerosols

!    real(r8), dimension(pcols,pverp) :: pmxrgnrf        ! temporary copy of pmxrgn
!    integer, dimension(pcols)  :: nmxrgnrf              ! temporary copy of nmxrgn
!    integer :: i, k ,ik             
!    integer :: istat
!  
!    real(r8) :: vis_dir
!    real(r8) :: vis_dif
!    real(r8) :: nir_dir
!    real(r8) :: nir_dif
!    real(r8), dimension(pver) ::  sw_dTdt     
!    real(r8), dimension(pver) ::  lw_dTdt     
!    real(r8), dimension(pverp) ::  sw_upflux   
!    real(r8), dimension(pverp) ::  sw_dnflux   
!    real(r8), dimension(pverp) ::  lw_upflux   
!    real(r8), dimension(pverp) ::  lw_dnflux   
!    real(r8), dimension(pcols) :: fsntoa        ! Net solar flux at TOA
!    real(r8), dimension(pcols) :: fsntoac       ! Clear sky net solar flux at TOA
!    real(r8), dimension(pcols) :: fsnirt        ! Near-IR flux absorbed at toa
!    real(r8), dimension(pcols) :: fsnrtc        ! Clear sky near-IR flux absorbed at toa
!    real(r8), dimension(pcols) :: fsnirtsq      ! Near-IR flux absorbed at toa >= 0.7 microns
!    real(r8), dimension(pcols) :: fsntc         ! Clear sky total column abs solar flux
!    real(r8), dimension(pcols) :: fsnsc         ! Clear sky surface abs solar flux
!    real(r8), dimension(pcols) :: fsdsc         ! Clear sky surface downwelling solar flux
!    real(r8), dimension(pcols) :: flut          ! Upward flux at top of model
!    real(r8), dimension(pcols) :: lwcf          ! longwave cloud forcing
!    real(r8), dimension(pcols) :: swcf          ! shortwave cloud forcing
!    real(r8), dimension(pcols) :: flutc         ! Upward Clear Sky flux at top of model
!    real(r8), dimension(pcols) :: flntc         ! Clear sky lw flux at model top
!    real(r8), dimension(pcols) :: flnsc         ! Clear sky lw flux at srf (up-down)
!    real(r8), dimension(pcols) :: fln200        ! net longwave flux interpolated to 200 mb
!    real(r8), dimension(pcols) :: fln200c       ! net clearsky longwave flux interpolated to 200 mb
!    real(r8), dimension(pcols,pverp) :: fns     ! net shortwave flux
!    real(r8), dimension(pcols,pverp) :: fcns    ! net clear-sky shortwave flux
!    real(r8), dimension(pcols) :: fsn200        ! fns interpolated to 200 mb
!    real(r8), dimension(pcols) :: fsn200c       ! fcns interpolated to 200 mb
!    real(r8), dimension(pcols,pverp) :: fsn     ! net shortwave flux
!    real(r8), dimension(pcols,pverp) :: fln     ! net longwave flux
!    real(r8), dimension(pcols,pverp) :: fnl     ! net longwave flux
!    real(r8), dimension(pcols,pverp) :: fcnl    ! net clear-sky longwave flux
!
!    real(r8) :: eccf                 ! Earth/sun distance factor
!    real(r8) :: delta        ! Solar declination angle
!
!    real(r8), pointer, dimension(:,:) :: co2mmr   ! co2   mass mixing ratio
!    real(r8), pointer, dimension(:,:) :: ch4mmr   ! ch4   mass mixing ratio
!    real(r8), pointer, dimension(:,:) :: o2mmr    ! o2    mass mixing ratio
!    real(r8), pointer, dimension(:,:) :: o3mmr    ! o3    mass mixing ratio
!
!    ! RT variables ------------------------------------------------
!    !the following would have dimension nazm_tshadow if shadows were taken into effect
!    integer, parameter :: ext_nazm_tshadow = 1                  ! take shadows into effect
!    real(r8), dimension(ext_nazm_tshadow) :: ext_cosz_horizon   ! cos of zenith angle of horizon 
!    real(r8), dimension(ext_nazm_tshadow) :: ext_TCx_obstruct
!    real(r8), dimension(ext_nazm_tshadow) :: ext_TCz_obstruct
!    !------------------------------------
!    integer :: ext_tslas_tog
!    integer :: ext_tshadow_tog
!  
!    real(r8) :: ext_solar_azm_ang
!    real(r8) :: ext_tazm_ang
!    real(r8) :: ext_tslope_ang
!    !real(r8), dimension(pver) :: ext_zc_static
!    !real(r8), dimension(pver) :: scale_height  
!
!    real(r8) :: ext_msdist
!    real(r8) :: ext_rtgt
!    ! ***NOTE: the following have CGS units:
!
!    real(r8), dimension(pcols,pverp) :: lwup_rad 
!    real(r8), dimension(pcols,pverp) :: lwdown_rad 
!    real(r8), dimension(pcols,pverp) :: swup_rad
!    real(r8), dimension(pcols,pverp) :: swdown_rad
!    real(r8) :: frac_day
!    real(r8) :: day_in_year
!
!------------------------------------------------------------------------
!
! Start Code
!
!    call t_startf ('radiation_tend')
!
!    lchnk = state%lchnk
!    ncol = state%ncol
!
!    !calday = get_curr_calday()
!
!    !itim = pbuf_old_tim_idx()
!    !ifld = pbuf_get_fld_idx('CLD')
!    !cfrc => pbuf(ifld)%fld_ptr(1,1:pcols,1:pver,lchnk,itim)
!    !ifld = pbuf_get_fld_idx('QRS')
!    !qrs => pbuf(ifld)%fld_ptr(1,1:pcols,1:pver,lchnk,1)
!    !ifld = pbuf_get_fld_idx('QRL')
!    !qrl => pbuf(ifld)%fld_ptr(1,1:pcols,1:pver,lchnk,1)
!    !ifld = pbuf_get_fld_idx('REL')
!    !rel  => pbuf(ifld)%fld_ptr(1,1:pcols,1:pver,lchnk, 1)
!    !ifld = pbuf_get_fld_idx('REI')
!    !rei  => pbuf(ifld)%fld_ptr(1,1:pcols,1:pver,lchnk, 1)
!    
!#if ( defined SCAM )
!      !  For CRM, make cloud equal to input observations:
!      if(switch(CRM_SW+1).and.have_cld) then
!        do k = 1,pver
!          cld(:ncol,k)= cldobs(k)
!        enddo
!      endif
!#endif
!
!    !
!    ! Cosine solar zenith angle for current time step
!    !
!    !call get_rlat_all_p(lchnk, ncol, clat)
!    !call get_rlon_all_p(lchnk, ncol, clon)
!    
!    !call get_curr_calday_rotation(frac_day, day_in_year)
!    !call zenith_rotation (frac_day, calday, clat, clon, coszrs, ncol)
!    
!    !
!    !  Compute cloud water/ice paths and optical properties for input to radiation
!    !  
!    !call t_startf('cldoptics')
!    !call param_cldoptics_calc(state, cfrc, landfrac, landm, icefrac, cicewp, cliqwp, &
!                              cldemis, rel, rei, pmxrgn, nmxrgn, snowh)   
!    !call t_stopf('cldoptics')
!
!    !call shr_orb_decl(calday, eccen, mvelpp, lambm0, obliqr, delta, eccf)
!    !ext_msdist=SHR_CONST_MSDIST2/eccf
!										
!    !-------for now, set these topography angles to 0 (ie, no topography blocking sun)
!    !-------should be moved to inside do loop if we want to use them
!    !-------azimuth can be calculated from cos(azim) = (cos(hr_ang)*cos(dec)*sin(lat)-sin(dec)*cos(lat))/cos(elev)
!    !-------the azimuth is the angle east of south when the hour angle, h, is negative (morning)
!    !-------and the angle west of south when the hour angle, h, is positive (afternoon).

!    !NOTES:  Why are all these set to 0?  What are all these variables?

!    ext_solar_azm_ang = 0.     !solar azimuthal angle? should not be zero
!    ext_tazm_ang = 0.
!    ext_tslope_ang = 0.
!    ext_tslas_tog = 0          
!    ext_tshadow_tog = 1        ! toggle shadowing 
!    ext_cosz_horizon = 0.
!    ext_TCx_obstruct = 0.
!    ext_TCz_obstruct = 0.
!    ext_rtgt = 1.  

!    sw_dTdt(:) = 0.    ! Initialize heating rate arrays
!    lw_dTdt(:) = 0.    !
 
!    lw_dnflux(:) = 0.   !
!    lw_upflux(:) = 0.   ! Initialize entire arrays for summing below
!    sw_upflux(:) = 0.   !
!    sw_dnflux(:) = 0.   !

!    qrs(:,:) = 0.
!    qrl(:,:) = 0.
!    lwup_rad(:,:) = 0.
!    lwdown_rad(:,:) = 0.
!    swup_rad(:,:) = 0.
!    swdown_rad(:,:) = 0.

!    !nstep = get_nstep()

!    !call rad_constituents_get('CO2', state, co2mmr)
!    !call rad_constituents_get('CH4', state, ch4mmr) 
!    !call rad_constituents_get('O2', state, o2mmr) 
!    !call rad_constituents_get('O3', state, o3mmr) 

!    do i = 1, ncol

!      ! scale height calculated at each layer midpoint in meters
!      ! scale_height(:) = SHR_CONST_BOLTZ*state%t(i,:)/(SHR_CONST_G*mwdry/SHR_CONST_AVOGAD) 

!      !ext_zc_static(:)=scale_height(:)*(-1000.)*log(state%pmid(i,:)/ps0)

!      call t_startf ('aerad_driver') 

!      call aerad_driver(state%q(i,:,1), co2mmr(i,:), ch4mmr(i,:), o2mmr(i,:), o3mmr(i,:) &
!                     ,cicewp(i,:), cliqwp(i,:), cfrc(i,:) &
!                     ,rei(i,:), rel(i,:) &
!                     ,srfflx_state2d%ts(i), state%ps(i), state%pmid(i,:) &
!                     ,state%pdel(i,:), state%t(i,:), state%pint(i,:) &
!                     ,coszrs(i), ext_msdist &
!                     ,srfflx_state2d%asdir(i), srfflx_state2d%aldir(i) &
!                     ,srfflx_state2d%asdif(i), srfflx_state2d%aldif(i) &
!                     ,ext_rtgt, ext_solar_azm_ang, ext_tazm_ang, ext_tslope_ang  &
!                     ,ext_tslas_tog, ext_tshadow_tog, ext_nazm_tshadow, ext_cosz_horizon  &
!                     ,ext_TCx_obstruct, ext_TCz_obstruct, state%zi(i,:), &
!                     ,sw_dTdt, lw_dTdt, lw_dnflux, lw_upflux, sw_upflux &
!                     ,sw_dnflux, vis_dir, vis_dif, nir_dir, nir_dif )       

!       call t_stopf ('aerad_driver') 

!       ftem(i,:) = sw_dTdt(:)       
!       ftem2(i,:) = lw_dTdt(:)      

!       lwup_rad(i,:) = lw_upflux(:)
!       lwdown_rad(i,:) = lw_dnflux(:)
!       swup_rad(i,:) = sw_upflux(:)
!       swdown_rad(i,:) = sw_dnflux(:)

!       ! Fluxes sent to land model
!       surface_state2d%sols(i) = vis_dir
!       surface_state2d%soll(i) = nir_dir
!       surface_state2d%solsd(i) = vis_dif
!       surface_state2d%solld(i) = nir_dif
!       surface_state2d%flwds(i) = lw_dnflux(pverp)
!       fsns(i) = sw_dnflux(pverp)-sw_upflux(pverp)
!       flns(i) = lw_upflux(pverp)-lw_dnflux(pverp)
       
!    enddo   ! ncol loop

!    fsn(:,:) = swdown_rad(:,:) - swup_rad(:,:)
!    fln(:,:) = lwup_rad(:,:) - lwdown_rad(:,:)
!    fsnsc(:) = fsns(:)
!    flnsc(:) = flns(:)
!    fsnt(:) = fsn(:,1)
!    fsntc(:) = fsnt(:)
!    flnt(:) = fln(:,1)
!    flntc(:) = flnt(:)
!    fsds(:) = swdown_rad(:,pverp)
!    qrs(:ncol,:pver) = ftem(:ncol,:pver)
!    qrl(:ncol,:pver) = ftem2(:ncol,:pver)
 
    !call outfld('QRS     ',qrs*SHR_CONST_CSEC  , pcols,lchnk)    ! [K/day]
    !call outfld('FSDS    ',fsds  ,pcols,lchnk)
    !call outfld('FSNT    ',fsnt  ,pcols,lchnk)
    !call outfld('FSNS    ',fsns  ,pcols,lchnk)
    !call outfld('QRL     ',qrl*SHR_CONST_CSEC   ,pcols,lchnk)    ! [K/day]
    !call outfld('FLNT    ',flnt  ,pcols,lchnk)
    !call outfld('FLUT    ',lwup_rad(:,2)  ,pcols,lchnk)
    !call outfld('FLNS    ',flns  ,pcols,lchnk)
    !call outfld('FUL     ',lwup_rad, pcols, lchnk)
    !call outfld('FDL     ',lwdown_rad, pcols, lchnk)
    !call outfld('FUS     ',swup_rad, pcols, lchnk)
    !call outfld('FDS     ',swdown_rad, pcols, lchnk)
    !call outfld('SOLS    ',surface_state2d%sols  ,pcols,lchnk)
    !call outfld('SOLL    ',surface_state2d%soll  ,pcols,lchnk)
    !call outfld('SOLSD   ',surface_state2d%solsd ,pcols,lchnk)
    !call outfld('SOLLD   ',surface_state2d%solld ,pcols,lchnk)


    ! Cloud cover diagnostics
    !call cldsav (lchnk, ncol, cfrc, state%pmid, cltot, cllow, clmed, clhgh, nmxrgn, pmxrgn)
    !
    ! Dump cloud field information to history tape buffer (diagnostics)
    !
    !call outfld('CLDTOT  ',cltot  ,pcols,lchnk)
    !call outfld('CLDLOW  ',cllow  ,pcols,lchnk)
    !call outfld('CLDMED  ',clmed  ,pcols,lchnk)
    !call outfld('CLDHGH  ',clhgh  ,pcols,lchnk)
    !call outfld('CLOUD   ',cfrc    ,pcols,lchnk) 

    ! Compute net radiative heating tendency

    !call radheat_tend(state, ptend, qrl, qrs, fsns, &
    !                  fsnt, flns, flnt, srfflx_state2d%asdir, net_flx)

    ! Compute heating rate for dtheta/dt
    !do k=1,pver
    !   do i=1,ncol
    !      ftem(i,k) = (qrs(i,k) + qrl(i,k))/cpair * &
    !                  (SHR_CONST_PSTD/state%pmid(i,k))**cappa
    !   enddo
    !enddo
    !call outfld('HR      ',ftem    ,pcols   ,lchnk   )

    ! convert radiative heating rates to Q*dp for energy conservation
    !if (conserve_energy) then
!DIR$ CONCURRENT
    !   do k =1 , pver
!DIR$ CONCURRENT
    !      do i = 1, ncol
    !         qrs(i,k) = qrs(i,k)*state%pdel(i,k)
    !         qrl(i,k) = qrl(i,k)*state%pdel(i,k)
    !      enddo
    !   enddo
    !endif
    ! 
    !call t_stopf ('radiation_tend')
!
!  end subroutine radiation_tend


!============================================================================

   subroutine init_ref
!------------------------------------------------------------------------
!
! Purpose: Initial reference value arrays
!
!------------------------------------------------------------------------

!NOTES: solvar_cambins not implemented yet
!#if ( defined WACCM_MOZART )
!      use solvar_interface, only: solvar_cambins
!#endif  

  implicit none

!------------------------------------------------------------------------
!
! Local Variables
!
    integer :: iq
    integer :: iw
    integer :: ip
    integer :: ig
    integer :: iv
    integer :: it

    real(r8) :: wn1 
    real(r8) :: wn2 
    real(r8) :: tempsol 
    real(r8) :: tempsol_iw
    real(r8) :: w1 
    real(r8) :: w2 
    real(r8) :: tempk
    logical :: found
 
!------------------------------------------------------------------------
!
! Start Code
! 
    
    ! ??move this to radgrid??
    ! Assign beginning and end wavelength range and point indices for each
    !  wavelength group (currently 2):
    !lw_iwbeg = 1                ! "longwave" group   
    !lw_iwend = nlw_wavlnrng     !
    !sw_iwbeg = lw_iwend+1       ! "shortwave" group
    !sw_iwend = ntot_wavlnrng    !

    !lw_ipbeg = 1                ! "longwave" group   
    !lw_ipend = nlw_gpt          !
    !sw_ipbeg = lw_ipend+1       ! "shortwave" group
    !sw_ipend = ntot_gpt         !

    ! Thermal and Solar spectral integration may be varied (i.e. reduced)
    ! to maximize accuracy and efficiency. 
    lw_iwbeg=1              ! thermal wavenumber bin begin, end index
    lw_iwend=ntot_wavlnrng    
    sw_iwbeg=1              ! solar wavenumber bin begin, end index
    sw_iwend=ntot_wavlnrng
    lw_ipbeg=1              ! thermal gauss point begin, end index
    lw_ipend=ntot_gpt
    sw_ipbeg=1              ! solar gauss point begin, end index
    sw_ipend=ntot_gpt

    ! Thermal and Solar spectral integration may be varied (i.e. reduced)
    ! to maximize accuracy and efficiency. 
    !lw_iwbeg=1              ! thermal wavenumber bin begin, end index
    !lw_iwend=17
    !sw_iwbeg=6              ! solar wavenumber bin begin, end index
    !sw_iwend=ntot_wavlnrng
    !lw_ipbeg=1              ! thermal gauss point begin, end index
    !lw_ipend=160
    !sw_ipbeg=65             ! solar gauss point begin, end index
    !sw_ipend=ntot_gpt       ! 248


    !write(*,*)   "INIT_REF: INITIALIZE GUASS POINT ARRAYS"

    ! Arrange g_weight(ntot_pt) array
    iq = 0
    do iw=1,ntot_wavlnrng 
      do ig=1, ngauss_pts(iw)
        iq = iq + 1
        if (ngauss_pts(iw) .eq. ngauss_8gpt) g_weight(iq) =  g_weight_8gpt(ig)
        if (ngauss_pts(iw) .eq. ngauss_16gpt) g_weight(iq) =  g_weight_16gpt(ig)
!!        if (iw .eq. 27) g_weight(iq) = solar_gweight_sp27(ig)
!!        if (iw .eq. 28) g_weight(iq) = solar_gweight_sp28(ig)
      enddo
    enddo

    ! Scale solar constant to namelist value
    ! In the future this data will be read from input file
!    sunm(:) = sunm(:)*scon/S0
    solarflux(:) = solarflux(:)*scon/S0
  
    write(*,*) "INIT_REF: total solar irradiance scaled to ",scon, "W m-2"
    write(*,*) "solar flux [W m-2] in each spectral interval"
    ! Calculate the "average" wavenumber (1/wavelength) <wavenum()> for each
    !  wavelength interval in each wavelength group [1/cm]:

    iq = 0
    ip = lw_ipbeg-1
    do iw=1,ntot_wavlnrng  ! "avg" wavenumber over each band
 !     iq = iq+1
 !     wn1 = wavenum_edge(iq)    ! [cm^-1]
 !     wn2 = wavenum_edge(iq+1)  ! 
 !     ! Integrated solar flux for each wavelength range:
 !     solflux(iw) = dble(solarflux(wn1,wn2,iq))   !W/m^2
      do ig=1,ngauss_pts(iw)
        ip = ip+1
        ! Gauss-weighted solar flux in each probability interval:
        !gw_solflux(ip) = solflux(iw)*g_weight(ip)
!write(*,*) ip !, solarflux(iw) !, g_weight(ip), gw_solflux(ip)
        gw_solflux(ip) = solarflux(iw)*g_weight(ip)

      enddo
    enddo

    !if (masterproc) then
      do iw=1, ntot_wavlnrng
        write(*,*) iw, solarflux(iw)
      enddo
      write(*,*) "TOTAL SOLAR FLUX:", SUM(solarflux)      
    !endif


    ! Calculate Planck function quantities as a function of wavelength point and
    !  temperature: 

    write(*,*)   "INIT_REF: CREATING PLANCK FUNCTION TABLE"
    ip = lw_ipbeg-1

    do iw=1,ntot_wavlnrng   ! Loop over all wavelength intervals

      ! Set these so that difference below is
      !   f(longer_wavlen)-f(shorter_wavlen),
      !   w1 = longer_wavlen, w2 = shorter_wavlen:
      !   h*c*nu/kc 
      w1 = dble(1.439*wavenum_edge(iw))      ! 1.439 ~ ((h*c)/k) ,
      w2 = dble(1.439*wavenum_edge(iw+1))    ! convert nu cm-1 to lambda (um)
 
      do ig=1,ngauss_pts(iw)
        ip = ip+1

        do it=tpft_beg,tpft_end
          tempk = tpft_invfinc*dble(it)   ! Temperature [K]

          !ptemp_itp, spectrally integrated radiance (isotropic) in each longwave guass interval
          !value weighted according to gauss weights, (B outside mapping)
          ptemp_itp(it,ip) = (PLANCKf(w1,tempk)-PLANCKf(w2,tempk))* &
                              g_weight(ip)*SHR_CONST_STEBOL/SHR_CONST_PI

        enddo
        ptemp_itp(1:tpft_beg,ip) = ptemp_itp(tpft_beg,ip)
      enddo
    enddo

    ! set two-stream model coefficients
    ! for solar stream (quadrature)
    U1Isol = sqrt3        
    U1I2sol = 0.5d0*sqrt3
    U1Ssol = 2.0*SHR_CONST_PI/U1Isol
    ! for thermal stream (hemispsheric mean)
    U1Iir = 2.d0    
    U1I2ir = 1.d0   ! 0.5d0*2.d0
    U1Sir = 2.0*SHR_CONST_PI/U1Iir

    !write(*,*) "INIT_REF: ADJUST CONTINUUM GAUSS POINTS"
    call adjust_continuum_gpt

    return

  end subroutine init_ref


!============================================================================

  subroutine adjust_continuum_gpt

!------------------------------------------------------------------------
!
! PURPOSE:  CO2 Continuum data sets are on 8 point gauss interval bin.  Adjust
!           to match number of gauss intervals used for major absorbing gases
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
!
! Start Code 
!

    ! Initialize reduced continuum k coefficient arrays    
    !kh2oself(:,:,:,:) = 0.0
    kco2cont(:) = 0.0
    !
    ! Gauss point adjustment for water vapor self continuum
    !
    itc = 0
    do iw=1, ntot_wavlnrng
      if (ngauss_pts(iw) .eq. 8) then ! no adjustment needed      
        do ig=1, ngauss_pts(iw)
          itc = itc + 1
          !kh2oself(itc,:,:,:) = kh2oself_8gpt(ig,:,:,:,iw)
          kco2cont(itc) = kco2cont_8gpt(ig,iw)
        enddo
      endif
      if (ngauss_pts(iw) .eq. 16) then
        do ig=1, ngauss_pts(iw)  
          itc = itc + 1
          !kh2oself(itc,:,:,:) = kh2oself_8gpt(map8to16gpt(ig),:,:,:,iw)
          kco2cont(itc) = kco2cont_8gpt(map8to16gpt(ig),iw)
        enddo
      endif      
    enddo
    
  end subroutine adjust_continuum_gpt

!============================================================================

  subroutine aerad_driver(ext_H2O, ext_CO2, ext_CH4, ext_O2, ext_O3, ext_H2, ext_N2, &
      ext_cicewp, ext_cliqwp, ext_cfrc, ext_rei, ext_rel, &
      ext_sfcT, ext_sfcP, ext_pmid, ext_pdel, ext_tmid, &
      ext_pint, ext_cosZ, ext_msdist, ext_asdir,  & 
      ext_aldir, ext_asdif, ext_aldif,  &
      ext_rtgt, ext_solar_azm_ang, ext_tazm_ang, ext_tslope_ang,  &
      ext_tslas_tog, ext_tshadow_tog, ext_nazm_tshadow, ext_cosz_horizon,  &
      ext_TCx_obstruct, ext_TCz_obstruct, ext_zint, &
      sw_dTdt, lw_dTdt, &
      lw_dnflux, lw_upflux, sw_upflux, sw_dnflux, &
      lw_dnflux_spectral, lw_upflux_spectral, sw_upflux_spectral, sw_dnflux_spectral, &
      vis_dir, vis_dif, nir_dir, nir_dif )

!------------------------------------------------------------------------
!
! Purpose: Driver for correlated K radiative transfer code.
!          Recieves column data from radiation_tend
!------------------------------------------------------------------------

    implicit none

!------------------------------------------------------------------------
!
! Input Arguments
!
    ! {intent IN}:
   
    integer, intent(in) :: ext_tslas_tog
    integer, intent(in) :: ext_tshadow_tog
    integer, intent(in) :: ext_nazm_tshadow

    real(r8), intent(in) :: ext_msdist
    real(r8), intent(in) :: ext_asdir          ! direct albedo (0.2-0.7 um) (from srfflx_stat2d%asdir)
    real(r8), intent(in) :: ext_aldir          ! direct albedo (0.7-4.0 um) (from srfflx_state2d%aldir)
    real(r8), intent(in) :: ext_asdif          ! diffuse albedo (0.2-0.7 um) (from srfflx_stat2d%asdif)
    real(r8), intent(in) :: ext_aldif          ! diffuse albedo (0.7-4.0 um) (from srfflx_state2d%aldif)
    real(r8), intent(in) :: ext_sfcT           ! surface temperature radiative  (from srfflx_state2d%ts)
    real(r8), intent(in) :: ext_sfcP           ! surface pressre (from state%ps) 
    real(r8), intent(in) :: ext_cosZ           ! cosine of the zenith angle 
    real(r8), intent(in) :: ext_rtgt           ! scaling used by models grid coordinate?
    real(r8), intent(in) :: ext_solar_azm_ang  ! solar azimuthal angle [rad]
    real(r8), intent(in) :: ext_tazm_ang       ! topographic slope and aspect angles from our model
    real(r8), intent(in) :: ext_tslope_ang     
    real(r8), intent(in), dimension(ext_nazm_tshadow) :: ext_cosz_horizon
    real(r8), intent(in), dimension(ext_nazm_tshadow) :: ext_TCx_obstruct
    real(r8), intent(in), dimension(ext_nazm_tshadow) :: ext_TCz_obstruct
    real(r8), intent(in), dimension(pverp) :: ext_zint     ! geopotential height at interfaces [m]
    real(r8), intent(in), dimension(pver) :: ext_pmid      ! pressure at midpoints [Pa]
    real(r8), intent(in), dimension(pver) :: ext_pdel      ! layer thickness [Pa]
    real(r8), intent(in), dimension(pver) :: ext_tmid      ! temperature at midpoints [K]
    real(r8), intent(in), dimension(pverp) :: ext_pint     ! pressure at interfaces

    real(r8), intent(in), dimension(pver) :: ext_H2O       ! specific humidy (from state%q < q) at midlayer [kg/kg]
    real(r8), intent(in), dimension(pver) :: ext_CO2       ! CO2 mass mixing ratio from state%q < co2mmr)   [kg/kg]
    real(r8), intent(in), dimension(pver) :: ext_CH4       ! CH4 mass mixing ratio from state%q < ch4mmr)   [kg/kg]
    real(r8), intent(in), dimension(pver) :: ext_O2        ! O2 mass mixing ratio from state%q < o2mmr)   [kg/kg]
    real(r8), intent(in), dimension(pver) :: ext_O3        ! O3 mass mixing ratio from state%q < o3mmr)   [kg/kg]
    real(r8), intent(in), dimension(pver) :: ext_H2        ! H2 mass mixing ratio from state%q < h2mmr)   [kg/kg]
    real(r8), intent(in), dimension(pver) :: ext_N2        ! N2 mass mixing ratio from state%q < h2mmr)   [kg/kg]
    real(r8), intent(in), dimension(pver) :: ext_cicewp    ! in cloud ice water path at layer midpoints [g/m2]
    real(r8), intent(in), dimension(pver) :: ext_cliqwp    ! in cloud liquid water path at layer midpoints [g/m2]
    real(r8), intent(in), dimension(pver) :: ext_cFRC      ! cloud fraction]
    real(r8), intent(in), dimension(pver) :: ext_rei       ! ice cloud particle effective drop size ice [microns]
    real(r8), intent(in), dimension(pver) :: ext_rel       ! liquid cloud drop effective drop size liquid [micron   

    real(r8), intent(out), dimension(pver) ::  sw_dTdt     
    real(r8), intent(out), dimension(pver) ::  lw_dTdt     

    real(r8), intent(out), dimension(pverp) ::  sw_upflux   
    real(r8), intent(out), dimension(pverp) ::  sw_dnflux   
    real(r8), intent(out), dimension(pverp) ::  lw_upflux   
    real(r8), intent(out), dimension(pverp) ::  lw_dnflux

    real(r8), intent(out), dimension(pverp,ntot_wavlnrng) ::  sw_upflux_spectral   
    real(r8), intent(out), dimension(pverp,ntot_wavlnrng) ::  sw_dnflux_spectral   
    real(r8), intent(out), dimension(pverp,ntot_wavlnrng) ::  lw_upflux_spectral   
    real(r8), intent(out), dimension(pverp,ntot_wavlnrng) ::  lw_dnflux_spectral

    real(r8), intent(out) ::  vis_dir
    real(r8), intent(out) ::  vis_dif
    real(r8), intent(out) ::  nir_dir
    real(r8), intent(out) ::  nir_dif


!------------------------------------------------------------------------
!
! Local Variables       
!
     real(r8), dimension(pverp) :: coldens       ! [molec m-2] wet columen amount per layer
     real(r8), dimension(pverp) :: coldens_dry   ! [molec m-2] dry columen amount per layer
     real(r8), dimension(pverp) :: qH2O          ! [kg/kg] H2O  mass mixing ratio mid layers 
     real(r8), dimension(pverp) :: qCO2          ! [kg/kg] CO2 mass mixing ratio at mid layers   
     real(r8), dimension(pverp) :: qCH4          ! [kg/kg] CH4 mass mixing ratio at mid layers   
     real(r8), dimension(pverp) :: qO2           ! [kg/kg] O2 mass mixing ratio at mid layers   
     real(r8), dimension(pverp) :: qO3           ! [kg/kg] O3 mass mixing ratio at mid layers   
     real(r8), dimension(pverp) :: qH2           ! [kg/kg] H2 mass mixing ratio at mid layers   
     real(r8), dimension(pverp) :: qN2           ! [kg/kg] H2 mass mixing ratio at mid layers   
     real(r8), dimension(pverp) :: cICE          ! [g/m2] in cloud ice water path at mid layers
     real(r8), dimension(pverp) :: cLIQ          ! [g/m2] in cloud liquid water path at mid layers
     real(r8), dimension(pverp) :: cfrc          ! cloud fraction at mid layers
     real(r8), dimension(pverp) :: REI           ! [microns] ice cloud particle effective radii at mid layers
     real(r8), dimension(pverp) :: REL           ! [microns] liquid cloud drop effective radii at mid layers
     real(r8), dimension(pverp) :: zlayer        ! [m] thickness of each vertical layer

     integer  :: swcut
     real(r8) :: tmp 
     real(r8) :: sinz
     real(r8) :: cosai 
     real(r8) :: atmp 
     real(r8) :: aint 
     real(r8) :: cosz_h 
     real(r8) :: x_obst
     real(r8) :: z_obst 
     real(r8) :: sw_cutoff
     real(r8) :: t1
     real(r8) :: t2
     real(r8), external :: tag_CPUtime

     integer :: iv
     integer :: ig
     integer :: ik
     integer :: k 
     integer :: i
     integer :: j
     integer :: iw
     integer :: ip
     integer :: ia0
     integer :: ia1
     integer :: ip_ibeg, ip_iend

     logical :: found
     logical :: horizon_extension 

     ! local variables used for computation (see Toon, 1989) annotate?
     real(r8), dimension(ntot_gpt,pverp) :: CK1sol, CK1ir
     real(r8), dimension(ntot_gpt,pverp) :: CK2sol, CK2ir       
     real(r8), dimension(ntot_gpt,pverp) :: CPBsol, CPBir
     real(r8), dimension(ntot_gpt,pverp) :: CMBsol, CMBir       
     real(r8), dimension(ntot_gpt,ngangles,pverp) :: Y3   
     real(r8), dimension(ntot_gpt,pverp) :: EM1sol, EM1ir
     real(r8), dimension(ntot_gpt,pverp) :: EM2sol, EM2ir       
     real(r8), dimension(ntot_gpt,pverp) :: EL1sol, EL1ir
     real(r8), dimension(ntot_gpt,pverp) :: EL2sol, EL2ir
     real(r8), dimension(ntot_gpt,2*pverp) :: AFsol, AFir
     real(r8), dimension(ntot_gpt,2*pverp) :: BFsol, BFir
     real(r8), dimension(ntot_gpt,2*pverp) :: EFsol, EFir
     real(r8), dimension(ntot_gpt,pverp) :: AKsol, AKir
     real(r8), dimension(ntot_gpt,pverp) :: GAMIsol, GAMIir
     real(r8), dimension(ntot_gpt,pverp) :: EE1sol, EE1ir       
     real(r8), dimension(ntot_gpt,pverp) :: B1sol, B1ir  ! gamma 1  \
     real(r8), dimension(ntot_gpt,pverp) :: B2sol, B2ir  ! gamma 2   two stream parameters
     real(r8), dimension(ntot_gpt,pverp) :: B3sol, B3ir  ! gamma 3  /
     real(r8), dimension(ntot_gpt,pverp) :: DIRECTsol, DIRECTir       
     real(r8), dimension(ntot_gpt,pverp) :: DIRECTU   
     real(r8), dimension(ntot_gpt,pverp) :: DIREC     
     real(r8), dimension(ntot_gpt,pverp) :: TAUL       ! optical depth of each layer
     real(r8), dimension(ntot_gpt,pverp) :: OPD        ! cumulative optical depth (top down)
     real(r8), dimension(ntot_gpt,pverp) ::  tau_gas    ! gas optical depth array
     real(r8), dimension(ntot_wavlnrng,pverp) ::  tau_ray    ! rayleigh optical depth
     real(r8), dimension(ntot_wavlnrng,pverp) ::  tau_n2n2cia  ! N2-N2 CIA optical depth
     real(r8), dimension(ntot_wavlnrng,pverp) ::  tau_h2n2cia  ! H2-N2 CIA optical depth
     real(r8), dimension(ntot_wavlnrng,pverp) ::  tau_h2h2cia  ! H2-H2 CIA optical depth
     real(r8), dimension(ntot_gpt,pverp) :: SOL       
     real(r8), dimension(ntot_gpt,pverp) :: W0         ! single scattering albedo        
     real(r8), dimension(ntot_gpt,pverp) :: G0         ! asymmetry parameter        

     ! Cloud optical properties
     real(r8), dimension(ncld_grp,ntot_gpt,pverp) ::  singscat_cld_mcica
     real(r8), dimension(ncld_grp,ntot_gpt,pverp) ::  asym_cld_mcica       
     real(r8), dimension(ncld_grp,ntot_gpt,pverp) ::  tau_cld_mcica
  
     ! stochastic bulk cloud properties (MCICA)
     real(r8), dimension(ntot_gpt,pverp) :: cFRC_mcica         
     real(r8), dimension(ntot_gpt,pverp) :: cICE_mcica
     real(r8), dimension(ntot_gpt,pverp) :: cLIQ_mcica

     real(r8), dimension(ntot_gpt,pverp) :: PTEMP     ! Planck function evaluated at each level  
     real(r8), dimension(ntot_gpt) :: PTEMPG          ! Planck function evaluated at ground
     real(r8), dimension(ntot_gpt,pverp) :: SLOPE     

     logical  :: part_in_tshadow
     real(r8), dimension(ntot_gpt) :: EMIS       ! Surface emissivity
     real(r8), dimension(ntot_gpt) :: RSFXdir    ! Surface reflectivity, direct radiation, gauss point grid
     real(r8), dimension(ntot_gpt) :: RSFXdif    ! Surface reflectivity, diffuse radiation, gauss point grid
     real(r8), dimension(ntot_wavlnrng) :: sfc_albedo_dir   ! Surface albedo, direct radiation, wavenumber grid
     real(r8), dimension(ntot_wavlnrng) :: sfc_albedo_dif   ! Surface albedo, diffuse radiation, wavenumber grid
     real(r8), dimension(ntot_wavlnrng) :: sfc_emiss
     real(r8) :: sflux_frac
     real(r8) :: sfc_tempk  
     real(r8) :: sfc_press
     real(r8) :: cos_mu           
     logical  :: sw_on     ! switch for togopography and day/night
     logical  :: beamSolar ! switch for beam

     real(r8), dimension(pver) ::  dzc          ! [kg m-2], column amount of mass 
     real(r8), dimension(pverp) ::  tint        ! [K] temperatures at level interfaces 
     real(r8), dimension(pverp) ::  tmid        ! [K] temperatures at level at mid layers + top (isothermal) 
     real(r8), dimension(pverp) ::  pmid        ! [Pa] pressure at level at mid layers + top (isothermal) 

     real(r8) :: dy     
!------------------------------------------------------------------------
!
! Start Code
!

    ! initialize internal RT variables
    CK1sol(:,:) = 0.0    ;    CK1ir(:,:) = 0.0
    CK2sol(:,:) = 0.0    ;    CK2ir(:,:) = 0.0
    CPBsol(:,:) = 0.0    ;    CPBir(:,:) = 0.0
    CMBsol(:,:) = 0.0    ;    CMBir(:,:) = 0.0
    Y3(:,:,:) = 0.0
    EM1sol(:,:) = 0.0    ;    EM1ir(:,:) = 0.0
    EM2sol(:,:) = 0.0    ;    EM2ir(:,:) = 0.0
    EL1sol(:,:) = 0.0    ;    EL1ir(:,:) = 0.0
    EL2sol(:,:) = 0.0    ;    EL2ir(:,:) = 0.0
    AFsol(:,:) = 0.0     ;    AFir(:,:) = 0.0
    BFsol(:,:) = 0.0     ;    BFir(:,:) = 0.0
    EFsol(:,:) = 0.0     ;    EFir(:,:) = 0.0
    AKsol(:,:) = 0.0     ;    AKir(:,:) = 0.0
    GAMIsol(:,:) = 0.0   ;    GAMIir(:,:) = 0.0
    EE1sol(:,:) = 0.0    ;    EE1ir(:,:) = 0.0
    B1sol(:,:) = 0.0     ;    B1ir(:,:) = 0.0
    B2sol(:,:) = 0.0     ;    B2ir(:,:) = 0.0
    B3sol(:,:) = 0.0     ;    B3ir(:,:) = 0.0
    DIRECTsol(:,:) = 0.0 ;    DIRECTir(:,:) = 0.0
    DIRECTU(:,:) = 0.0
    DIREC(:,:) = 0.0
    TAUL(:,:) = 0.0
    OPD(:,:) = 0.0
    tau_gas(:,:) = 0.0
    tau_n2n2cia(:,:) = 0.0
    tau_h2n2cia(:,:) = 0.0
    tau_h2h2cia(:,:) = 0.0
    tau_cld_mcica(:,:,:) = 0.0
    singscat_cld_mcica(:,:,:) = 0.0
    asym_cld_mcica(:,:,:) = 0.0


    
    ! Fraction of the interplanetary solar flux at top of atmosphere:
    sflux_frac = dble(1./ext_msdist)    ! [1/AU^2]
   
    ! Set amount in layer above top defined model boundary
    qH2O(1) = ext_H2O(1)       ! H2O vapor mass concentration (specific humdity) [kg/kg]
    qCO2(1) = ext_CO2(1)       ! CO2 mass mixing ratio [kg/kg]
    qCH4(1) = ext_CH4(1)       ! CH4 mass mixing ratio [kg/kg]
    qO2(1) = ext_O2(1)         ! O2 mass mixing ratio [kg/kg]
    qO3(1) = ext_O3(1)         ! O3 mass mixing ratio [kg/kg]
    qH2(1) = ext_H2(1)         ! H2 mass mixing ratio [kg/kg]
    qN2(1) = ext_N2(1)         ! N2 mass mixing ratio [kg/kg]
    cICE(1) = ext_cicewp(1)    ! in cloud ice water path [g/m2]
    cLIQ(1) = ext_cliqwp(1)    ! in cloud liquid water path [g/m2]
    cFRC(1) = ext_cfrc(1)      ! cloud fraction
    REI(1) = ext_rei(1)        ! ice cloud particle effective radii [microns]
    REL(1) = ext_rel(1)        ! liquid cloud dropeffective radii [microns]
    tmid(1) = ext_tmid(1)      ! temperatures [K]
    pmid(1) = ext_pint(1)      ! pressure [Pa]  
   
    ! Set amount in midlayers elsewhere 
    ! Top layer =2, bottom layer =pverp
    do k=2, pverp
      qH2O(k) = ext_H2O(k-1)
      qCO2(k) = ext_CO2(k-1) 
      qCH4(k) = ext_CH4(k-1) 
      qO2(k) = ext_O2(k-1)
      qO3(k) = ext_O3(k-1)
      qH2(k) = ext_H2(k-1)
      qN2(k) = ext_N2(k-1)
      cICE(k) = ext_cicewp(k-1) 
      cLIQ(k) = ext_cliqwp(k-1) 
      cFRC(k) = ext_cfrc(k-1) 
      REI(k) = ext_rei(k-1)
      REL(k) = ext_rel(k-1)
      tmid(k) = ext_tmid(k-1)
      pmid(k) = ext_pmid(k-1)
    enddo    

    ! Set ground (surface) values:
    sfc_tempk = ext_sfcT   ! [K]
    sfc_press = ext_sfcP   ! [pa]

    ! Set interface temperatures
    tint(pverp) = sfc_tempk
    tint(1) = ext_tmid(1)

    do k = 2, pver 
      dy = (log10(ext_pint(k)) - log10(ext_pmid(k))) / (log10(ext_pmid(k-1)) - log10(ext_pmid(k)))
      tint(k) = ext_tmid(k) - dy * (ext_tmid(k) - ext_tmid(k-1))
      !write(*,*) k,dy
    enddo    

!! DIAGNOSTIC OUTPUT
!    write(*,*) "pint, pmid"
!    !do k=1,pverp
!    !  write(*,*) k, ext_pint(k), pmid(k), pmid(k)*(1.0-qh2o(k))
!    !enddo
!    write(*,*) "tint, tmid"
!    do k=1,pverp
!      write(*,*) k, tint(k), tmid(k)
!    enddo
!! DIAGNOSTIC OUTPUT

    ! Define molecular column density at each layer [molec m-2]  
  
    ! Set column density in layer above top boundary
!    write(*,*) "MWDRY:", mwdry
    coldens(1) = (ext_pint(1)*SHR_CONST_AVOGAD)/(mwdry*SHR_CONST_G)
    coldens_dry(1) = (ext_pint(1)*SHR_CONST_AVOGAD)/(mwdry*SHR_CONST_G)*(1.0-qh2o(1))
    ! Set column density for other mid layers
    do k=2, pverp   
      coldens(k) = (ext_pdel(k-1)*SHR_CONST_AVOGAD)/(mwdry*SHR_CONST_G)  ! defined from wet air mass
      coldens_dry(k) = (ext_pdel(k-1)*SHR_CONST_AVOGAD)/(mwdry*SHR_CONST_G)*(1.0-qh2o(k)) ! defined from dry air mass
    enddo    
    !write(*,*) "COLDENS", coldens

    ! Define mass column density in each layer [kg m-2]
    dzc(:) = ext_pdel(:)/SHR_CONST_G  

    ! Define height of each layer [m]
    zlayer(1) = 0.0   !thickness of layer with lower boundary at model top is zero 
    do k=2, pverp
      zlayer(k-1) = (ext_zint(k-1) - ext_zint(k))    
    enddo

    !NOTES: surface albedo set from srfflx_state2d variable
    !Current implementation uses a gray albedo for shortwave
    !albedo and a gray albedo for longwave. The demarcation 
    !between to the two regimes is  ~5 um radiation.   

    ! Set surface direct albedo, diffuse albedo, and emissivity
    do iw=1,ntot_wavlnrng    ! Loop over relevant wavelength intervals 
      if (wavenum_edge(iw) .le. 2000) then 
        sfc_albedo_dir(iw) = 0.0 
        sfc_albedo_dif(iw) = 0.0
      endif
      if (wavenum_edge(iw) .gt. 2000 .and. wavenum_edge(iw) .le. 13000) then   ! "infrared" :: 2080 cm-1 to 12850 cm-1 (4.80769 um to 0.778210 um)
        sfc_albedo_dir(iw) = ext_aldir       
        sfc_albedo_dif(iw) = ext_aldif
      endif
      if (wavenum_edge(iw) .ge. 13000) then     ! "solar" :: 12850 cm-1 to 50000 cm-1 (0.778210 um to 0.2 um)
        sfc_albedo_dir(iw) = ext_asdir       
        sfc_albedo_dif(iw) = ext_asdif       
      endif
      !sfc_emiss(iw) = 1.-sfc_albedo_dir(iw)      
      sfc_emiss(iw) = 1.0 - sfc_albedo_dir(iw) ! Only used for thermal emission, assume 1 everywhere

!    write(*,*) iw, " sfc_albedo dir/dif", sfc_albedo_dir(iw), sfc_albedo_dif(iw), " sfc_emiss ", sfc_emiss(iw)
    enddo
    


    !*** Set cosine of incident solar angle (zenith angle); take into account
    !     (crudely) any horizon extension:
    cos_mu = dble(ext_cosZ)     ! [none]

    if(ext_tshadow_tog == 1) then   ! Need to check for horizon extension
    
      ! Calculate cosz_horizon array location/interpolative indices (using
      !  solar azimuth values):
    
      atmp = ext_solar_azm_ang/(2.0*SHR_CONST_PI/real(ext_nazm_tshadow)) 
      ia0 = int(atmp)
      aint = atmp-real(ia0)   ! Location crucial to avoid ia0 complications
    
      if(ia0 == 0) then
        ia0 = ext_nazm_tshadow
        ia1 = 1
      elseif(ia0 == ext_nazm_tshadow) then
        ia1 = 1
      else
        ia1 = ia0+1
      endif

      ! Calc horizon cosz, etc. for current solar azimuth via interpolation:
      z_obst = (1.-aint)*ext_TCz_obstruct(ia0)+ aint*ext_TCz_obstruct(ia1)      
      if(z_obst < 0.0) then   ! Horizon extension
        horizon_extension = .TRUE.
        cosz_h = (1.-aint)*ext_cosz_horizon(ia0)+aint*ext_cosz_horizon(ia1)
        if(real(cos_mu) > cosz_h .and. cos_mu < 4.d-2) then
          ! Sun is above extended horizon, but still low in sky;
          !  cosZ=0.04 =~ 87.71deg
          cos_mu = min(cos_mu-dble(cosz_h),4.d-2)
        endif
      else
        horizon_extension = .FALSE.
      endif
    endif

    ! Only do shortwave calculation if sun is sufficiently above horizon;
    !  set probability point interval indices for use in radiative transfer
    !  calcs:

    if(cos_mu > 1.d-6) then
      sw_on = .TRUE.       ! Sun sufficiently high, do  shortwave calculation

      ! Topographic slope/aspect correction:

      if(ext_tslas_tog == 1) then
        sinz = sin(acos(cos_mu))
        ! If cosai < 0., point not sunlit
        cosai = (cos_mu*cos(ext_tslope_ang))+  &
                 (sinz*sin(ext_tslope_ang)*  &
                 cos(ext_solar_azm_ang-ext_tazm_ang))
      endif

      ! Topographic shadowing:

      part_in_tshadow  = .FALSE.   ! Initialize   ! False, no part in shadow

      if(ext_tshadow_tog == 1 .and. .NOT.horizon_extension) then           ! possible shadowing
        cosz_h = (1.-aint)*ext_cosz_horizon(ia0)+aint*ext_cosz_horizon(ia1)
        x_obst = (1.-aint)*ext_TCx_obstruct(ia0)+aint*ext_TCx_obstruct(ia1)
        sw_cutoff = z_obst-(x_obst*tan(asin(real(cos_mu))))
      
        if(sw_cutoff > 0. .and. sw_cutoff < ext_zint(pverp)*ext_rtgt) then
          part_in_tshadow = .TRUE.
          swcut = 0
          atmp = 0.
          shadow_top: do ik=1,pver        ! Find top of shadow
            aint = ext_zint(ik)*ext_rtgt
            if(sw_cutoff < aint) then
              if(ik /= 1) then
                if(abs(aint-sw_cutoff) < abs(atmp-sw_cutoff)) then
                  swcut = pver-ik+1
                else
                  swcut = pver-ik+2
                endif
              else
                swcut = pver-ik+1
              endif
                exit shadow_top
            endif
            atmp = aint
          enddo shadow_top   ! Column sunlit above layer 'swcut'
       
          if(swcut == 0) then
            write(*,*) 'ERROR: unable to find shadow top.'
            stop '*** ERROR: aerad_driver01 ***'
          endif
      
        endif
      endif

    else     ! (cos_mu < 1.d-6) 
      sw_on = .FALSE.       ! Sun below horizon, do only longwave
    endif

    !call t_startf ('calc_gasopd') 
    call calc_gasopd(tmid, pmid/100.0, ext_pdel/100.0, coldens, coldens_dry, qH2O, qCO2, qCH4, qO2, qO3, qH2, qN2, &
                     zlayer*100.0, tau_gas, tau_ray, tau_n2n2cia, tau_h2n2cia, tau_h2h2cia)
    !call t_stopf ('calc_gasopd') 

    !call t_startf ('calc_aeropd')
    !call calc_aeropd( )
    !call t_stopf ('calc_aeropd')

    !call t_startf('calc_cldopd')
    !call calc_cldopd(ext_pint, cICE, cLIQ, REI, REL, cFRC, tau_cld_mcica, singscat_cld_mcica, & 
    !                 asym_cld_mcica, cFRC_mcica, cICE_mcica, cICE_mcica ) 
    !call t_stopf('calc_cldopd')  

    !do k=1,pverp
    !  do iw=lw_ipbeg, sw_ipend
    !    write(*,*) "liq", iw,k,tau_cld_mcica(1,iw,k), singscat_cld_mcica(1,iw,k),asym_cld_mcica(1,iw,k)
    !    write(*,*) "ice", iw,k,tau_cld_mcica(2,iw,k), singscat_cld_mcica(2,iw,k),asym_cld_mcica(2,iw,k)
    !  enddo 
    !enddo

    !call t_startf ('rad_precalc') 
    call rad_precalc(pmid/100.0, tmid, tint, swcut, tau_gas, tau_ray, &
                     tau_n2n2cia, tau_h2n2cia, tau_h2h2cia, &
                     tau_cld_mcica, singscat_cld_mcica, asym_cld_mcica, &
                     part_in_tshadow, sfc_albedo_dir, sfc_albedo_dif, sfc_emiss, & 
                     sflux_frac, sfc_tempk, cos_mu, sw_on, &                 
                     Y3, TAUL, OPD, PTEMP, PTEMPG, SLOPE, SOL, W0, G0, EMIS, RSFXdir, RSFXdif)
    !call t_stopf ('rad_precalc') 

    !call t_startf ('two_stream') 
    !call t_stopf ('two_stream') 
    beamSolar = .true.  ! do solar calculation, all wavenlengths, two-stream quadrature
    ip_ibeg = sw_ipbeg
    ip_iend = sw_ipend
    call two_stream(TAUL, W0, G0, RSFXdir, RSFXdif, beamSolar, ip_ibeg, ip_iend, &
                    EM1sol, EM2sol, EL1sol, EL2sol, &
                    AFsol, BFsol, EFsol, AKsol, &
                    GAMIsol, B1sol, B2sol, EE1sol)

    call add_txrad(EM1sol, EM2sol, EL1sol, EL2sol, &
                   AFsol, BFsol, EFsol, B1sol, B2sol, AKsol, &
                   TAUL, OPD, PTEMP, PTEMPG, SLOPE, SOL, W0, G0, EMIS, RSFXdir, RSFXdif, cos_mu, sw_on, &
                   ip_ibeg, ip_iend, beamSolar, CK1sol, CK2sol, CPBsol, CMBsol, B3sol, DIRECTsol)    


    beamSolar = .false.  ! do thermal calculation, all wavelengths, two-stream hemispheric mean
    ip_ibeg = lw_ipbeg
    ip_iend = lw_ipend
    call two_stream(TAUL, W0, G0, RSFXdir,RSFXdif, beamSolar, ip_ibeg, ip_iend, &
                    EM1ir, EM2ir, EL1ir, EL2ir, &
                    AFir, BFir, EFir, AKir, &
                    GAMIir, B1ir, B2ir, EE1ir)

    call add_txrad(EM1ir, EM2ir, EL1ir, EL2ir, &
                   AFir, BFir, EFir, B1ir, B2ir, AKir, &
                   TAUL, OPD, PTEMP, PTEMPG, SLOPE, SOL, W0, G0, EMIS, RSFXdir, RSFXdif, cos_mu, sw_on, &
                   ip_ibeg, ip_iend, beamSolar, CK1ir, CK2ir, CPBir, CMBir, B3ir, DIRECTir)    

    !call t_startf ('add_txrad') 
    !call t_stopf ('add_txrad') 

    !call t_startf ('refine_lwflux') 
    call refine_lwflux(CK1ir, CK2ir, Y3, AKir, GAMIir, B3ir, EE1ir, &
                       TAUL, PTEMP, PTEMPG, SLOPE, EMIS, RSFXdir, RSFXdif, cos_mu, DIRECTU, DIREC)
    !call t_stopf ('refine_lwflux') 

    ! Calculate final fluxes / heating rates
    !call t_startf ('rad_postcalc') 
    call rad_postcalc(CK1sol, CK2sol, CPBsol, CMBsol, &
                      EM1sol, EM2sol, EL1sol, EL2sol, &
                      DIRECTsol, DIRECTU, DIREC, dzc, swcut, part_in_tshadow, sw_on, &
                      sw_dTdt, lw_dTdt, &
                      lw_dnflux, lw_upflux, sw_upflux, sw_dnflux, &
                      lw_dnflux_spectral, lw_upflux_spectral, sw_upflux_spectral, sw_dnflux_spectral, &
                      vis_dir, vis_dif, nir_dir, nir_dif) 
    !call t_stopf ('rad_postcalc') 

    write(*,*) "Surface downwelling fluxes"
    write(*,*) "vis_dir", vis_dir !/sw_dnflux(pverp)
    write(*,*) "vis_dif", vis_dif !/sw_dnflux(pverp)
    write(*,*) "nir_dir", nir_dir !/sw_dnflux(pverp)
    write(*,*) "nir_dif", nir_dif !/sw_dnflux(pverp)
    write(*,*) "total direct", vis_dir+nir_dir
    write(*,*) "total diffuse", vis_dif+nir_dif
    write(*,*) "SW DN TOA/SURF", sw_dnflux(1), sw_dnflux(pverp)
    write(*,*) "SW UP TOA/SURF", sw_upflux(1), sw_upflux(pverp)
    write(*,*) "LW DN TOA/SURF", lw_dnflux(1), lw_dnflux(pverp)
    write(*,*) "LW UP TOA/SURF", lw_upflux(1), lw_upflux(pverp)
    

    return

  end subroutine aerad_driver


!============================================================================

  subroutine rad_precalc(pmid, tmid, tint, swcut, tau_gas, tau_ray, &
                         tau_n2n2cia, tau_h2n2cia, tau_h2h2cia, &
                         tau_cld_mcica, singscat_cld_mcica, asym_cld_mcica, &
                         part_in_tshadow, sfc_albedo_dir, sfc_albedo_dif, & 
                         sfc_emiss, sflux_frac, sfc_tempk, cos_mu, sw_on, &
                         Y3, TAUL, OPD, PTEMP, PTEMPG, SLOPE, SOL, W0, G0, EMIS, RSFXdir, RSFXdif)

!------------------------------------------------------------------------
!
! Purpose: Calculates quantities needed before the main radiative transfer 
!          calculation can be performed
!                                        
!------------------------------------------------------------------------

    implicit none

!------------------------------------------------------------------------
!
! Arguments
! 
  real(r8), intent(in), dimension(pverp) :: pmid       ! [mb] pressures at mid layers)
  real(r8), intent(in), dimension(pverp) :: tmid       ! [K] temperatures at mid layers + top (isothermal)
  real(r8), intent(in), dimension(pverp) :: tint       ! [K] temperatures at level interfaces 
  integer, intent(in) :: swcut
  real(r8), intent(in), dimension(ntot_gpt,pverp) ::  tau_gas 
  real(r8), intent(in), dimension(ntot_wavlnrng,pverp) ::  tau_ray   
  real(r8), intent(in), dimension(ntot_wavlnrng,pverp) ::  tau_n2n2cia   
  real(r8), intent(in), dimension(ntot_wavlnrng,pverp) ::  tau_h2n2cia   
  real(r8), intent(in), dimension(ntot_wavlnrng,pverp) ::  tau_h2h2cia   
  real(r8), intent(in), dimension(ncld_grp,ntot_gpt,pverp) ::  tau_cld_mcica
  real(r8), intent(in), dimension(ncld_grp,ntot_gpt,pverp) ::  singscat_cld_mcica
  real(r8), intent(in), dimension(ncld_grp,ntot_gpt,pverp) ::  asym_cld_mcica
  logical, intent(in) :: part_in_tshadow 
  real(r8), intent(in), dimension(ntot_wavlnrng) :: sfc_albedo_dir
  real(r8), intent(in), dimension(ntot_wavlnrng) :: sfc_albedo_dif
  real(r8), intent(in), dimension(ntot_wavlnrng) :: sfc_emiss
  real(r8), intent(in) :: sflux_frac
  real(r8), intent(in) :: sfc_tempk
  real(r8), intent(in) :: cos_mu           
  logical, intent(in) :: sw_on

  real(r8), intent(out), dimension(ntot_gpt,ngangles,pverp) :: Y3   
  real(r8), intent(out), dimension(ntot_gpt,pverp) :: TAUL      
  real(r8), intent(out), dimension(ntot_gpt,pverp) :: OPD       
  real(r8), intent(out), dimension(ntot_gpt,pverp) :: PTEMP     
  real(r8), intent(out), dimension(ntot_gpt) :: PTEMPG
  real(r8), intent(out), dimension(ntot_gpt,pverp) :: SLOPE      
  real(r8), intent(out), dimension(ntot_gpt,pverp) :: SOL       
  real(r8), intent(out), dimension(ntot_gpt,pverp) :: W0
  real(r8), intent(out), dimension(ntot_gpt,pverp) :: G0        
  real(r8), intent(out), dimension(ntot_gpt) :: EMIS
  real(r8), intent(out), dimension(ntot_gpt) :: RSFXdir
  real(r8), intent(out), dimension(ntot_gpt) :: RSFXdif

!------------------------------------------------------------------------
!
! Local Variables
! 
    real(r8) :: den
    real(r8) :: x
    real(r8) :: f0_ig
    real(r8) :: w0_ig
    real(r8) :: g0_ig
    real(r8) :: taul_ig
    real(r8) :: tmp

    integer :: it
    integer :: ig
    integer :: k
    integer :: ip
    integer :: iw
    integer :: ib
    integer :: k_1
    integer :: ia

!   The vertical structure used by this radiation code is as follows:
!
!   Rad Index = Layer Boundaries
!   Atmos Index = Layer Centers
!
!   Layer                       Rad Index   Atmos Index
!
!    top -----------------------  j = 1
!         - - - - - - - - - - -   j = 2       k = 1
!        -----------------------
!         - - - - - - - - - - -               k = 2
!        -----------------------
!         - - - - - - - - - - -
!        -----------------------
!         - - - - - - - - - - -               k = NZ-1
!        -----------------------
!         - - - - - - - - - - -   j = 2*NZ-1  k = NZ
! bottom -----------------------  j = 2*NZ
!============================================================
 
!------------------------------------------------------------------------
!
! Start Code
!     

    ! Map/average important parameters to gauss points (it):
    it = 0
    do iw=1,ntot_wavlnrng    ! Loop over wavenumber bands
      do ig=1,ngauss_pts(iw)
        it = it+1
        EMIS(it) = sfc_emiss(iw)
        RSFXdir(it) = sfc_albedo_dir(iw)
        RSFXdif(it) = sfc_albedo_dif(iw)        
      enddo
    enddo    

    if(sw_on) then
      if(part_in_tshadow) then   ! Eliminate direct "shortwave" flux within shadow

        it = sw_ipbeg-1
        do iw=sw_iwbeg,sw_iwend  
          do ig=1,ngauss_pts(iw)
            it = it+1
            SOL(it,1:swcut-1) = gw_solflux(it)*sflux_frac
            SOL(it,swcut:pverp) = 0.d0
          enddo
        enddo        

      else    ! No topographic shadowing in current column
       
        it = sw_ipbeg-1
        do iw=sw_iwbeg,sw_iwend
          do ig=1,ngauss_pts(iw)
            it = it+1
            SOL(it,1:pverp) = gw_solflux(it)*sflux_frac
          enddo
        enddo

      endif
    endif
    
    !if (sw_on) write(6,*) 'sw_on is true in precalc, solin1 is',solin1,'solinl1 is',solinl1
    
  
    OPD(1:ntot_gpt,:) = 0.d0     ! Initialize necessary array portions in
    TAUL(1:ntot_gpt,:) = 0.d0    !  anticipation of summing below

    !open (unit = 15, file = "opdepths.txt", status = "replace", iostat = openstatus)
    !write (15,*) "it k opdepth(it,k) TAUL(it,k) singscat asym"

    do k=1,pverp     ! Loop over all layer BOUNDARIES  
   
      k_1 = max(1,k-1)   ! index for level above (except for k=1, of course)
      it = 0

      do iw=1,ntot_wavlnrng      ! Only necessary wavelength bands
        do ig=1,ngauss_pts(iw)
          it = it+1
 
          ! Combine the aerosol and gas optical parameters to get the total
          !  "effective" layer parameters:

          taul_ig = tau_gas(it,k) + tau_ray(iw,k) + tau_n2n2cia(iw,k) + tau_h2n2cia(iw,k) + tau_h2h2cia(iw,k)
          w0_ig = tau_ray(iw,k)
          g0_ig = 0.                             

          ! Add cloud optical depths
          do ip =1,2
            taul_ig = taul_ig + tau_cld_mcica(ip,it,k)      
            w0_ig = w0_ig + singscat_cld_mcica(ip,it,k) * tau_cld_mcica(ip,it,k)
            g0_ig = g0_ig + asym_cld_mcica(ip,it,k) * singscat_cld_mcica(ip,it,k) * tau_cld_mcica(ip,it,k)
          enddo
         
          ! Add aerosol optical depths here


          if(taul_ig < SMALLd) then   ! Clip if optical depth too small
            taul_ig = SMALLd
          endif

          w0_ig = w0_ig/taul_ig        ! "total" weighted average singscat_albd
          w0_ig = min(1.d0-SMALLd,w0_ig)  !

          if(w0_ig > SMALLd) then
            g0_ig = g0_ig/(w0_ig*taul_ig)  ! "total" weighted average asym_fact
          else
            !?any point to this?    wsum = SMALLe   !? w0_ig = epsilon
            g0_ig = 0.d0           ! "total" weighted average asym_fact
          endif

          ! Apply delta-Eddington scaling to truncate the forward scattering
          !  lobe of the scattering phase function (effectively leave the
          !  energy within the truncated portion in the direct beam), as
          !  oringinally described in Joseph et al. (1976):
 
          f0_ig = g0_ig*g0_ig
          den = 1.d0-w0_ig*f0_ig
          TAUL(it,k) = taul_ig*den              ! "delta-scaled" total opt_depth
          W0(it,k) = (1.d0-f0_ig)*w0_ig/den     ! "delta-scaled" total singscat_albd
          G0(it,k) = g0_ig/(1.d0+g0_ig)         ! "delta-scaled" total asym_fact
                                                ! Actually '(g0_ig-f0_ig)/(1.d0-f0_ig)'

          ! Cumulative "Delta-scaled" optical depth, sums through layer k
          OPD(it,k) = OPD(it,k_1)+TAUL(it,k)    
                                                 
          !write (15, fmt = '(1X, I3, I3, ES15.7, ES15.7, ES15.7, ES15.7)') & 
          !  it,k,OPD(it,k),TAUL(it,k),w0(it,k),g0(it,k) 

        enddo
      enddo

      do ia=1,ngangles
        do it=1,ntot_gpt
          X = TAUL(it,k)/g_angle(ia)
          if(X <= 400.d0) then
            Y3(it,ia,k) = min(1.d0,exp(-X))
          else
            Y3(it,ia,k) = 0.d0    ! would be effectively zero anyway
          endif
        enddo 
      enddo
    enddo

    !close (15)

    ! Calculate Planck function for current temperatures profiles
    !write(6,*)'calling dplanck'
    call dplanck(pmid, tmid, tint, TAUL, sfc_tempk, PTEMP, PTEMPG, SLOPE)

    !write(6,*)'ending rad_precalc'
    return

  end subroutine rad_precalc

!============================================================================

  subroutine dplanck(pmid, tmid, tint, TAUL, sfc_tempk, PTEMP, PTEMPG, SLOPE)

!------------------------------------------------------------------------
!
! Purpose: Calculate Planck function and its derivative at the ground and
!          at all altitudes.
!          Planck temperatures are cacluated at the interface temperatures.
!          Planck slopes are calculated between midlayer temperatures.                                        
!------------------------------------------------------------------------

    implicit none

!------------------------------------------------------------------------
!
! Arguments
!          
    real(r8), intent(in), dimension(pverp) :: pmid       ! [mb] pressures at mid layers
    real(r8), intent(in), dimension(pverp) :: tmid       ! [K] temperatures at mid layers (isothermal top)
    real(r8), intent(in), dimension(pverp) :: tint       ! [K] temperatures at level interfaces 
    real(r8), intent(in), dimension(ntot_gpt,pverp) :: TAUL   ! Layer optical depth      
    real(r8), intent(in) :: sfc_tempk  
    real(r8), intent(out), dimension(ntot_gpt,pverp) :: PTEMP   ! Planck function at each layer (isothermal top)
    real(r8), intent(out), dimension(ntot_gpt) :: PTEMPG        ! Planck function evaluated at ground
    real(r8), intent(out), dimension(ntot_gpt,pverp) :: SLOPE     

!------------------------------------------------------------------------
!
! Local Variables
!          
    real(r8), dimension(ntot_gpt,pverp) :: ptemp_midlayer
    real(r8) :: wl1
    real(r8) :: wl2
    real(r8) :: ft
    integer :: ip,iw,ig
    integer :: k
    integer :: k_1
    integer :: it

!------------------------------------------------------------------------
!
! Start Code
!          

    ! Calculate stuff based on the wavelength-dependent Planck function at the
    !  ground surface and for all vertical levels:

    write(*,*) "surface planck function [W m-2] in each spectral interval"

    do ip=lw_ipbeg,lw_ipend   ! Loop over longwave gauss points

      it = int(sfc_tempk*tpft_finc)      
      ft = dble(sfc_tempk*tpft_finc-real(it))
      
      ! Interpolate between table values:
      PTEMPG(ip) = ptemp_itp(it,ip)*(1.d0-ft)+ptemp_itp(it+1,ip)*ft
      !write(*,*) "sfc_tempk", sfc_tempk, "it", it, "ft", ft, PTEMPG(ip)

    enddo

    !---- Diagnostic output for 1d simulations -----
    ip=0
    do iw=1,ntot_wavlnrng
      do ig=1,ngauss_pts(iw)  
        ip=ip+1        
      enddo
      write(*,*) iw,PTEMPG(ip)/g_weight(ip)*SHR_CONST_PI
    enddo
    write(*,*) "TOTAL PLANCK FUNCTION:", SUM(PTEMPG)*SHR_CONST_PI
    write(*,*) "SURFACE PLANCK FUNCTION:",SHR_CONST_STEBOL*sfc_tempk**4
    !----         ----

    ! Set Planck function at interfaces
    do ip=lw_ipbeg,lw_ipend   ! Loop over all longwave gauss points
      do k=1,pverp      ! Loop over all layers

        ! Interface temperatures for planck 
        it = int(tint(k)*tpft_finc)        
        ft = dble(tint(k)*tpft_finc-real(it))

        ! If temperature exceeds table range, for to be maximum    
        if (tint(k)*tpft_finc > tpft_end) then
           it = tpft_end-1
        endif
    
        if (tint(k)*tpft_finc < tpft_beg) then
           write(*,*) "negative temperature ", tmid(k),k, "tmid(k), k"
           it = tpft_beg-1
        endif

        ! Interpolate between table values:        
        PTEMP(ip,k) = ptemp_itp(it,ip)*(1.d0-ft)+ptemp_itp(it+1,ip)*ft   

        ! Midlayer temperatures Dplanck/Dtau SLOPE calculation 
        it = int(tmid(k)*tpft_finc)
        ft = dble(tmid(k)*tpft_finc-real(it))
        ptemp_midlayer(ip,k) = ptemp_itp(it,ip)*(1.d0-ft)+ptemp_itp(it+1,ip)*ft   
 
        k_1 = max(1,k-1)   ! index for level above (except for k=1, of course)
        if(TAUL(ip,k) > SMALLd) then
!          SLOPE(ip,k) = (PTEMP(ip,k)-PTEMP(ip,k_1))/TAUL(ip,k)
           SLOPE(ip,k) = 2.0 * (ptemp_midlayer(ip,k)-ptemp_midlayer(ip,k_1)) / (TAUL(ip,k)+TAUL(ip,k_1)) ! DIML
        else
          SLOPE(ip,k) = 0.d0          
        endif
      enddo  
    
    enddo
!! Diagnostics
!!   do k=1,pverp
!!    write(*,*) "----------------"
!!    write(*,*) k,tmid(k), tint(k), PTEMP(1,k), SLOPE(1,k)
!!   enddo
!! Diagnostics
      
    return

  end subroutine dplanck

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


!============================================================================

  subroutine two_stream(TAUL, W0, G0, RSFXdir, RSFXdif, beamSolar, ip_ibeg, ip_iend, &
                        EM1, EM2, EL1, EL2, AF, BF, EF, AK, GAMI, B1, B2, EE1)

!------------------------------------------------------------------------
!
! Purpose: Defines matrix properties and sets up coefficients that do 
!          not depend on solar zenith angle or atmospheric temperature
!                                                  
!------------------------------------------------------------------------

    implicit none

!------------------------------------------------------------------------
!
! Arguments

     real(r8), intent(in), dimension(ntot_gpt,pverp) :: TAUL      
     real(r8), intent(in), dimension(ntot_gpt,pverp) :: W0         
     real(r8), intent(in), dimension(ntot_gpt,pverp) :: G0         
     real(r8), intent(in), dimension(ntot_gpt) :: RSFXdir
     real(r8), intent(in), dimension(ntot_gpt) :: RSFXdif
     logical, intent(in) :: beamSolar
     integer, intent(in) :: ip_ibeg
     integer, intent(in) :: ip_iend
     real(r8), intent(out), dimension(ntot_gpt,pverp) :: EM1
     real(r8), intent(out), dimension(ntot_gpt,pverp) :: EM2
     real(r8), intent(out), dimension(ntot_gpt,pverp) :: EL1
     real(r8), intent(out), dimension(ntot_gpt,pverp) :: EL2
     real(r8), intent(out), dimension(ntot_gpt,2*pverp) :: AF        
     real(r8), intent(out), dimension(ntot_gpt,2*pverp) :: BF        
     real(r8), intent(out), dimension(ntot_gpt,2*pverp) :: EF        
     real(r8), intent(out), dimension(ntot_gpt,pverp) :: AK        
     real(r8), intent(out), dimension(ntot_gpt,pverp) :: GAMI      
     real(r8), intent(out), dimension(ntot_gpt,pverp) :: B1        
     real(r8), intent(out), dimension(ntot_gpt,pverp) :: B2  
     real(r8), intent(out), dimension(ntot_gpt,pverp) :: EE1       
!------------------------------------------------------------------------
!
! Local Variables
!                                                  
   
    real(r8) :: x1
    integer :: ip
    integer :: k
    integer :: kd
    real(r8) u1i, u1i_2
!------------------------------------------------------------------------
!
! Start Code
!                                                  

    ! HERE WE DEFINE LAYER PROPERTIES FOLLOWING GENERAL SCHEME
    !  OF MEADOR AND WEAVOR. THEN WE SET UP LAYER PROPERTIES
    !  NEEDED FOR MATRIX:

    if (beamSolar) then 
      u1i_2 = U1I2sol
    else
      u1i_2 = U1I2ir
    endif

    do k=1,pverp
      do ip=ip_ibeg,ip_iend

        ! THESE ARE FOR QUADRATURE AND HEMISPHERIC MEAN
        B1(ip,k) = u1i_2*(2.d0-W0(ip,k)*(1.d0+G0(ip,k)))
        B2(ip,k) = u1i_2*W0(ip,k)*(1.d0-G0(ip,k))
        AK(ip,k) = SQRT(ABS(B1(ip,k)**2-B2(ip,k)**2))
        GAMI(ip,k) = B2(ip,k)/(B1(ip,k)+AK(ip,k))
        x1 = AK(ip,k)*TAUL(ip,k)         !
        if(x1 <= 400.d0) then            !
          EE1(ip,k) = exp(-x1)           !        
        else                             ! TIM add mod
          EE1(ip,k) = 0.d0               !
        endif                            !
        EL1(ip,k) = 1.d0+GAMI(ip,k)*EE1(ip,k)
        EM1(ip,k) = 1.d0-GAMI(ip,k)*EE1(ip,k)
        EL2(ip,k) = GAMI(ip,k)+EE1(ip,k)
        EM2(ip,k) = GAMI(ip,k)-EE1(ip,k)

      enddo
    
    enddo

    ! WE SEEK TO SOLVE AX(L-1)+BX(L)+EX(L+1) = D.
    !  L=2N FOR EVEN L, L=N+1 FOR ODD L. THE MEAN INTENSITY (TMI/4PI)
    !  AND THE NET FLUX (FNET) ARE RELATED TO X'S AS NOTED IN ADD.
    !  FIRST WE SET UP THE COEFFICIENTS THAT ARE INDEPENDENT OF SOLAR
    !  ANGLE OR TEMPERATURE: A(I),B(I),E(I). D(I) IS DEFINED IN ADD:
    
    k = 0
    do kd=2,2*pverp-1,2
      k = k+1
      do ip=ip_ibeg,ip_iend
        !           HERE ARE THE EVEN MATRIX ELEMENTS
        AF(ip,kd) = EM1(ip,k+1)*EL1(ip,k)-EM2(ip,k+1)*EL2(ip,k)
        BF(ip,kd) = EM1(ip,k+1)*EM1(ip,k)-EM2(ip,k+1)*EM2(ip,k)
        EF(ip,kd) = EL1(ip,k+1)*EM2(ip,k+1)-EL2(ip,k+1)*EM1(ip,k+1)
        !           HERE ARE THE ODD MATRIX ELEMENTS EXCEPT FOR THE TOP.
        AF(ip,kd+1) = EM1(ip,k)*EL2(ip,k)-EL1(ip,k)*EM2(ip,k)
        BF(ip,kd+1) = EL1(ip,k+1)*EL1(ip,k)-EL2(ip,k+1)*EL2(ip,k)
        EF(ip,kd+1) = EL2(ip,k)*EM2(ip,k+1)-EL1(ip,k)*EM1(ip,k+1)
      enddo
    enddo

    ! HERE ARE THE TOP AND BOTTOM BOUNDARY CONDITIONS AS WELL AS THE
    !  BEGINNING OF THE TRIDIAGONAL SOLUTION DEFINITIONS. I ASSUME
    !  NO DIFFUSE RADIATION IS INCIDENT AT UPPER BOUNDARY:

    do ip=ip_ibeg,ip_iend
      AF(ip,1) = 0.d0
      BF(ip,1) = EL1(ip,1)
      EF(ip,1) = -EM1(ip,1)
      AF(ip,2*pverp) = EL1(ip,pverp)-RSFXdif(ip)*EL2(ip,pverp)
      BF(ip,2*pverp) = EM1(ip,pverp)-RSFXdif(ip)*EM2(ip,pverp)
      EF(ip,2*pverp) = 0.d0
    enddo
  
    return    

  end subroutine two_stream


!============================================================================

  subroutine add_txrad (EM1, EM2, EL1, EL2, AF, BF, EF, B1, B2, AK, TAUL, OPD, & 
                        PTEMP, PTEMPG, SLOPE, SOL, W0, G0, EMIS, RSFXdir, RSFXdif, cos_mu, sw_on, &
                        ip_ibeg, ip_iend, beamSolar, CK1, CK2, CPB, CMB, B3, DIRECT )

!------------------------------------------------------------------------
!
! Purpose: Defines source terms, forms matrix for multiple layers and solves
!          the tri-diagonal equations to obtain mean intensity and net flux 
!                                                            
!------------------------------------------------------------------------

    implicit none

!------------------------------------------------------------------------
!
! Arguments
!          
  real(r8), intent(in), dimension(ntot_gpt,pverp) :: EM1       
  real(r8), intent(in), dimension(ntot_gpt,pverp) :: EM2       
  real(r8), intent(in), dimension(ntot_gpt,pverp) :: EL1       
  real(r8), intent(in), dimension(ntot_gpt,pverp) :: EL2       
  real(r8), intent(in), dimension(ntot_gpt,2*pverp) :: AF        
  real(r8), intent(in), dimension(ntot_gpt,2*pverp) :: BF
  real(r8), intent(in), dimension(ntot_gpt,2*pverp) :: EF
  real(r8), intent(in), dimension(ntot_gpt,pverp) :: B1        
  real(r8), intent(in), dimension(ntot_gpt,pverp) :: B2        
  real(r8), intent(in), dimension(ntot_gpt,pverp) :: AK        
  real(r8), intent(in), dimension(ntot_gpt,pverp) :: TAUL      
  real(r8), intent(in), dimension(ntot_gpt,pverp) :: OPD       
  real(r8), intent(in), dimension(ntot_gpt,pverp) :: PTEMP
  real(r8), intent(in), dimension(ntot_gpt) :: PTEMPG 
  real(r8), intent(in), dimension(ntot_gpt,pverp) :: SLOPE     
  real(r8), intent(in), dimension(ntot_gpt,pverp) :: SOL       
  real(r8), intent(in), dimension(ntot_gpt,pverp) :: W0  
  real(r8), intent(in), dimension(ntot_gpt,pverp) :: G0 
  real(r8), intent(in), dimension(ntot_gpt) :: EMIS
  real(r8), intent(in), dimension(ntot_gpt) :: RSFXdir
  real(r8), intent(in), dimension(ntot_gpt) :: RSFXdif
  real(r8), intent(in) :: cos_mu           
  logical, intent(in) :: sw_on
  logical, intent(in) :: beamSolar
  integer, intent(in) :: ip_ibeg
  integer, intent(in) :: ip_iend
  real(r8), intent(out), dimension(ntot_gpt,pverp) :: CK1       
  real(r8), intent(out), dimension(ntot_gpt,pverp) :: CK2       
  real(r8), intent(out), dimension(ntot_gpt,pverp) :: CPB       
  real(r8), intent(out), dimension(ntot_gpt,pverp) :: CMB
  real(r8), intent(out), dimension(ntot_gpt,pverp) :: B3        
  real(r8), intent(out), dimension(ntot_gpt,pverp) :: DIRECT       

!------------------------------------------------------------------------
!
! Local Variables
!          
    
    real(r8), dimension(ntot_gpt,pverp) :: CPP
    real(r8), dimension(ntot_gpt,pverp) :: CM
    real(r8), dimension(ntot_gpt,pverp) :: EE3
    real(r8), dimension(ntot_gpt,pverp) :: EL3
    real(r8), dimension(ntot_gpt,2*pverp) :: AS
    real(r8), dimension(ntot_gpt,2*pverp) :: DF
    real(r8), dimension(ntot_gpt,2*pverp) :: DS
    real(r8), dimension(ntot_gpt,2*pverp) :: XK

    real(r8), dimension(ntot_gpt) :: sfcs
    real(r8) :: DUo
    real(r8) :: B4
    real(r8) :: X2
    real(r8) :: X3
!    real(r8) :: X4
    real(r8) :: C1
    real(r8) :: C2
    real(r8) :: CP1
    real(r8) :: CM1
    real(r8) :: X
    real(r8) :: x_lim
    integer :: KINDEX
    integer :: k
    integer :: ip
    integer ::k_1
    integer :: kd

!------------------------------------------------------------------------
!
! Start Code
!        
    ! THIS SUBROUTINE FORMS THE MATRIX FOR THE MULTIPLE LAYERS AND
    !  USES A TRIDIAGONAL ROUTINE TO FIND RADIATION IN THE ENTIRE
    !  ATMOSPHERE.
    !

    if (beamSolar) then
    ! ******************************
    ! *   CALCULATIONS FOR SOLAR   *
    ! ******************************
    ! use delta eddington scaling (Joseph et al. 1976)

      if(.not. sw_on) return  ! if sun is below horizon, skip calculation

      DUo = 1.d0/cos_mu

      do k=1,pverp

        do ip=sw_ipbeg,sw_ipend  
    
          B3(ip,k) = 0.5d0*(1.d0-sqrt3*G0(ip,k)*cos_mu)     ! TIM mod
          B4 = 1.d0-B3(ip,k)
          X2 = TAUL(ip,k)*DUo
          if(X2 > 400.d0) then   !
            X2 = 400.d0          ! TIM add
          endif                  ! 
          X3 = OPD(ip,k)*DUo  
          if(X3 > 400.d0) then   !
            X3 = 400.d0          ! TIM add
          endif                  !         
          EL3(ip,k) = exp(-X3)*SOL(ip,k)        ! Beers law, layer bottom, tau=opd
          EE3(ip,k) = exp(-(X3-X2))*SOL(ip,k)   ! Beers law, layer top, tau=opd-taul
          DIRECT(ip,k) = cos_mu*EL3(ip,k)   
          C1 = B1(ip,k)-DUo
          if(abs(C1) < SMALLd) then
           C1 = sign(SMALLd,C1)
          endif
          C2 = AK(ip,k)*AK(ip,k)-DUo*DUo
          if(abs(C2) <= SMALLd) then
            C2 = SMALLd
          endif
          CP1 = W0(ip,k)*(B3(ip,k)*C1+B4*B2(ip,k))/C2   
          CM1 = (CP1*B2(ip,k)+W0(ip,k)*B4)/C1
          CPB(ip,k) = CP1*EL3(ip,k)     ! C+ lower boundary (tau = opd)
          CPP(ip,k) = CP1*EE3(ip,k)     ! C+ upper boundary (tau = opd-taul)   
          CMB(ip,k) = CM1*EL3(ip,k)     ! C- lower boundary (tau = opd)
          CM(ip,k) = CM1*EE3(ip,k)      ! C- upper boundary (tau = opd-taul) 

        enddo
      enddo

      ! CALCULATE SFCS, SHORTWAVE SOURCE AT THE BOTTOM (REFLECTION):
    
      do ip=sw_ipbeg,sw_ipend
        sfcs(ip) = DIRECT(ip,pverp)*RSFXdir(ip)
      enddo

    else  ! beamSolar

    ! ******************************
    ! * CALCULATIONS FOR INFRARED. *
    ! ******************************
      do k=1,pverp
        kindex = max(1,k-1)
        do ip=lw_ipbeg,lw_ipend
          B3(ip,k) = 1.d0/(B1(ip,k)+B2(ip,k))
          CPP(ip,k) = (PTEMP(ip,KINDEX)+SLOPE(ip,k)*B3(ip,k))*U1Sir  ! C+ upper boundary (B0)
          CPB(ip,k) = CPP(ip,k)+SLOPE(ip,k)*TAUL(ip,k)*U1Sir         ! C+ lower boundary (B+db*tau)
          CM(ip,k) = (PTEMP(ip,KINDEX)-SLOPE(ip,k)*B3(ip,k))*U1Sir   ! C- upper boundary (B0)
          CMB(ip,k) = CM(ip,k)+SLOPE(ip,k)*TAUL(ip,k)*U1Sir          ! C- lower boundary (B+db*tau)
          EL3(ip,k)    = 0.d0
          DIRECT(ip,k) = 0.d0
          EE3(ip,k)    = 0.d0
        enddo
      enddo
 
      do ip=lw_ipbeg,lw_ipend
        sfcs(ip) = EMIS(ip)*PTEMPG(ip)*SHR_CONST_PI     
      enddo
    
    endif !beamSolar

    k = 0
    do kd=2,2*pverp-1,2
      k = k+1
      do ip=ip_ibeg,ip_iend
        ! HERE ARE THE ODD MATRIX ELEMENTS EXCEPT FOR THE TOP.
        DF(ip,kd+1) = EL2(ip,k)*(CPP(ip,k+1)-CPB(ip,k))+  &   ! TIM mod
           EL1(ip,k)*(CMB(ip,k)-CM(ip,k+1))

        ! HERE ARE THE EVEN MATRIX ELEMENTS
        DF(ip,kd) = (CPP(ip,k+1)-CPB(ip,k))*EM1(ip,k+1)-  &   ! TIM mod
           (CM(ip,k+1)-CMB(ip,k))*EM2(ip,k+1)

      enddo
    enddo

    ! HERE ARE THE TOP AND BOTTOM BOUNDARY CONDITIONS AS WELL AS THE
    !  BEGINNING OF THE TRIDIAGONAL SOLUTION DEFINITIONS. I ASSUME NO
    !  DIFFUSE RADIATION IS INCIDENT AT THE TOP:
    do ip=ip_ibeg,ip_iend
      DF(ip,1) = -CM(ip,1)
      DF(ip,2*pverp) = SFCS(ip)+RSFXdif(ip)*CMB(ip,pverp)-CPB(ip,pverp)
      DS(ip,2*pverp) = DF(ip,2*pverp)/BF(ip,2*pverp)
      AS(ip,2*pverp) = AF(ip,2*pverp)/BF(ip,2*pverp)
    enddo

    ! ********************************************
    ! *     WE SOLVE THE TRIDIAGONAL EQUATIONS   *
    ! ********************************************
    x_lim = tiny(x)*1.d100   ! Lower limit for 'x' below        ! TIM add
    do kd=2,2*pverp
      k = 2*pverp+1-kd     ! TIM add
      do ip=ip_ibeg,ip_iend
        x = (BF(ip,k)-EF(ip,k)*AS(ip,k+1))   ! TIM mod
        if(abs(x) < x_lim) then     !Ã
          x = sign(1.d0,x)*x_lim    ! TIM add
        endif                       !
        X = 1.d0/x                  !
        AS(ip,k) = AF(ip,k)*X
        DS(ip,k) = (DF(ip,k)-EF(ip,k)*DS(ip,k+1))*X
      enddo
    enddo

    do ip=ip_ibeg,ip_iend
      XK(ip,1) = DS(ip,1)
    enddo
    
    do kd=2,2*pverp
      do ip=ip_ibeg,ip_iend
        XK(ip,kd) = DS(ip,kd)-AS(ip,kd)*XK(ip,kd-1)
      enddo
    enddo

    ! ***************************************************************
    !    CALCULATE LAYER COEFFICIENTS, NET FLUX AND MEAN INTENSITY
    ! ***************************************************************
    
    do k=1,pverp
      do ip=ip_ibeg,ip_iend
        CK1(ip,k) = XK(ip,2*k-1)
        CK2(ip,k) = XK(ip,2*k)
        !### FNET and TMI not used in this application (similar quantities are instead
        !     computed in `postcalc`)
        !FNET(ip,k)  = CK1(ip,k)  *( EL1(ip,k) -EL2(ip,k))   +  &
        !    CK2(ip,k) *( EM1(ip,k)-EM2(ip,k) ) + CPB(ip,k) -  &
        !    CMB(ip,k) - DIRECT(ip,k)
        !TMI(ip,k)   =  EL3(ip,k) + U1I(ip) *( CK1(ip,k)  *  &
        !   ( EL1(ip,k) + EL2(ip,k))   +  &
        !   CK2(ip,k) *( EM1(ip,k)+EM2(ip,k) ) +  &
        !   CPB(ip,k) + CMB(ip,k) )
      enddo
    enddo
    
    return
  end subroutine add_txrad
                  

!============================================================================

  subroutine refine_lwflux (CK1, CK2, Y3, AK, GAMI, B3, EE1, TAUL, PTEMP, PTEMPG, SLOPE, EMIS, RSFXdir,RSFXdif, cos_mu, DIRECTU, DIREC)

!------------------------------------------------------------------------
!
! Purpose: Calculate upward and downward radiances and intensities using 
!          Gauss Quadrature angles and weights          
!          original                                                  
!------------------------------------------------------------------------

     implicit none

!------------------------------------------------------------------------
!
! Local Variables
!
    real(r8), intent(in), dimension(ntot_gpt,pverp) :: CK1       
    real(r8), intent(in), dimension(ntot_gpt,pverp) :: CK2       
    real(r8), intent(in), dimension(ntot_gpt,ngangles,pverp) :: Y3   
    real(r8), intent(in), dimension(ntot_gpt,pverp) :: AK        
    real(r8), intent(in), dimension(ntot_gpt,pverp) :: GAMI      
    real(r8), intent(in), dimension(ntot_gpt,pverp) :: B3        
    real(r8), intent(in), dimension(ntot_gpt,pverp) :: EE1           
    real(r8), intent(in), dimension(ntot_gpt,pverp) :: TAUL      
    real(r8), intent(in), dimension(ntot_gpt,pverp) :: PTEMP     
    real(r8), intent(in), dimension(ntot_gpt) :: PTEMPG     
    real(r8), intent(in), dimension(ntot_gpt,pverp) :: SLOPE     
    real(r8), intent(in), dimension(ntot_gpt) :: EMIS
    real(r8), intent(in), dimension(ntot_gpt) :: RSFXdir
    real(r8), intent(in), dimension(ntot_gpt) :: RSFXdif
    real(r8), intent(in) :: cos_mu

    real(r8), intent(out), dimension(ntot_gpt,pverp) :: DIRECTU   
    real(r8), intent(out), dimension(ntot_gpt,pverp) :: DIREC    

!
!     *  Input               :  PTEMP, SLOPE, Y3, B3, EE1, EE2     *
!     *  Output              :  Irradiances: DIREC, DIRECTU
!
!------------------------------------------------------------------------
!
! Local Variables
!          
    real(r8), dimension(ntot_gpt,ngangles,pverp) :: uintent  ! [ntot_gpt,ngangles,pverp]
    real(r8), dimension(ntot_gpt,ngangles,pverp) :: dintent  ! [ntot_gpt,ngangles,pverp]

    integer :: openstatus
    real(r8), dimension(ntot_gpt,ngangles,pverp) :: Y1   
    real(r8), dimension(ntot_gpt,ngangles,pverp) :: Y2   
    real(r8), dimension(ntot_gpt,ngangles,pverp) :: Y4   
    real(r8), dimension(ntot_gpt,ngangles,pverp) :: Y8   
    real(r8), dimension(ntot_gpt,pverp) :: Y5        
    real(r8), dimension(ntot_gpt,pverp) :: A1        
    real(r8), dimension(ntot_gpt,pverp) :: A2        
    real(r8), dimension(ntot_gpt,pverp) :: A3        
    real(r8), dimension(ntot_gpt,pverp) :: A4        
    real(r8), dimension(ntot_gpt,pverp) :: A7        
    real(r8) :: X4
    real(r8) :: YA
    real(r8) :: YB
    real(r8) :: CKP
    real(r8) :: CKM
    real(r8) :: cos_mu1
    integer :: KINDEX
    integer :: k
    integer :: kd
    integer :: ip
    integer :: ia
!------------------------------------------------------------------------
!
! Start Code

    do k=1,pverp
      kindex = max(1,k-1)
      do ip=lw_ipbeg,lw_ipend
        ! HERE WE DO NO SCATTERING COEFFICIENTS
        A3(ip,k) = PTEMP(ip,KINDEX)*2.0*SHR_CONST_PI
        A4(ip,k) = 2.0*SHR_CONST_PI*SLOPE(ip,k)
        A7(ip,k) = A3(ip,k)
        Y5(ip,k) = A4(ip,k)*TAUL(ip,k)
      enddo
      ! HERE WE DO SCATTERING
      do ip=lw_ipbeg,lw_ipend
        X4 = SLOPE(ip,k)*(2.0*SHR_CONST_PI*B3(ip,k)-U1Sir)
        A1(ip,k) = U1Iir-AK(ip,k)
        A2(ip,k) = GAMI(ip,k)*(AK(ip,k)+U1Iir)
        A3(ip,k) = A3(ip,k)+X4
        A7(ip,k) = A7(ip,k)-X4
      enddo
    enddo

    ! CALCULATIONS FOR ALL GAUSS POINTS:
    do k=1,pverp   
      do ia=1,ngangles
       
        ! HERE WE DO NO SCATTERING COEFFS

        do ip=lw_ipbeg,lw_ipend
          !%waste CPU%         Y1(ip,ia,k)  =  0.d0
          !%  Y2(ip,ia,k)  =  0.d0
          Y4(ip,ia,k) = A7(ip,k)-A4(ip,k)*g_angle(ia)
          Y8(ip,ia,k) = A3(ip,k)+A4(ip,k)*g_angle(ia)
          !%      enddo

          ! HERE WE DO SCATTERING

          !% do ip=lw_ipbeg,lw_ipend
          YA = A1(ip,k)*(Y3(ip,ia,k)-EE1(ip,k))/(AK(ip,k)*g_angle(ia)-1.d0)
          YB = A2(ip,k)*(1.d0-EE1(ip,k)*Y3(ip,ia,k))/  &
               (AK(ip,k)*g_angle(ia)+1.d0)
          CKP = CK1(ip,k)+CK2(ip,k)
          CKM = CK1(ip,k)-CK2(ip,k)
          Y1(ip,ia,k) = CKP*YB+CKM*YA
          Y2(ip,ia,k) = CKP*YA+CKM*YB

        enddo

      enddo
    enddo

    ! INITIALIZE IRRADIANCES TO ZERO
    do  k=1, pverp
      do ip=lw_ipbeg,lw_ipend
        DIREC(ip,k) = 0.d0
        DIRECTU(ip,k) = 0.d0
      enddo   
    enddo

    ! DIREC IS DOWNWARD IRRADIANCE. DIRECTU IS UPWARD IRRADIANCE.
    !  CALCULATE DINTENT THE DOWNWARD RADIANCE AND DIREC THE DOWNWARD IRRADIANCE
   
    ! BOUNDARY CONDITIONS: DOWNWARD IRRADIANCE, RADIANCE AT TOA (k = camtop) was 1
    write(*,*) "CAMTOP", camtop
    cos_mu1 = max(cos_mu,0.0)

    do ia=1,ngangles
      do ip=lw_ipbeg,lw_ipend
        
        DINTENT(ip,ia,camtop) = (1.d0-Y3(ip,ia,camtop))*Y4(ip,ia,camtop)+Y1(ip,ia,camtop)

        DIREC(ip,camtop) = DIREC(ip,camtop)+DINTENT(ip,ia,camtop)*g_ang_weight(ia)
        
      enddo
    enddo


    ! DINTENT IS DOWNWARD RADIANCE * TwoPI. DIREC IS THE DOWNWARD IRRADIANCE.   
    !  CALCULATE FOR REST OF ATMOSPHERE.
    do k=camtop+1, pverp     !was k=2
      do ia=1,ngangles
        do ip=lw_ipbeg,lw_ipend
          DINTENT(ip,ia,k) = DINTENT(ip,ia,k-1)*Y3(ip,ia,k)  & 
               +Y1(ip,ia,k)+Y5(ip,k)+(1.d0-Y3(ip,ia,k))*Y4(ip,ia,k)             

          DIREC(ip,k) = DIREC(ip,k)+DINTENT(ip,ia,k)*g_ang_weight(ia)
        enddo
      enddo             
    enddo

    ! UINTENT IS THE UPWARD RADIANCE * TwoPI. DIRECTU IS THE UPWARD IRRADIANCE.
    !  ASSUME THAT THE REFLECTIVITY IS LAMBERT.
    
    ! BOUNDARY CONDITIONS: UPWARD IRRADIANCE, RADIANCE AT BOTTOM (k = pverp)

    do ia=1,ngangles
      do ip=lw_ipbeg,lw_ipend

        UINTENT(ip,ia,pverp) = EMIS(ip)*PTEMPG(ip)*2.0*SHR_CONST_PI+RSFXdif(ip)*DIREC(ip,pverp)*2.0d0

        DIRECTU(ip,pverp) = DIRECTU(ip,pverp)+UINTENT(ip,ia,pverp)*g_ang_weight(ia)

      enddo
    enddo

    ! CALCULATE FOR THE REST OF THE ATMOSPHERE
    ! ITERATED FROM BOTTOM UP
    kd = pverp
    do k=camtop+1,pverp   ! was k=2
      kd = kd-1
      do ia=1,ngangles
        do ip=lw_ipbeg,lw_ipend

          UINTENT(ip,ia,kd) = (UINTENT(ip,ia,kd+1)-Y5(ip,kd+1))  &
                              *Y3(ip,ia,kd+1)+Y2(ip,ia,kd+1)+  &
                              (1.d0-Y3(ip,ia,kd+1))*Y8(ip,ia,kd+1)          

          DIRECTU(ip,kd) = DIRECTU(ip,kd)+UINTENT(ip,ia,kd)*g_ang_weight(ia)

        enddo
      enddo
    enddo

    return

  end subroutine refine_lwflux

!============================================================================


  subroutine rad_postcalc (CK1, CK2, &
                           CPB, CMB, &
                           EM1, EM2, EL1, EL2, &
                           DIRECT, DIRECTU, DIREC, dzc, swcut, part_in_tshadow, sw_on, &
                           sw_dTdt, lw_dTdt, & 
                           lw_dnflux, lw_upflux, sw_upflux, sw_dnflux, &
                           lw_dnflux_spectral, lw_upflux_spectral, sw_upflux_spectral, sw_dnflux_spectral, &
                           vis_dir, vis_dif, nir_dir, nir_dif )

!------------------------------------------------------------------------
!
! Purpose: Calculate total radiative fluxes; put in a form suitable for output 
!                                                                      
!------------------------------------------------------------------------
 
    implicit none

!------------------------------------------------------------------------
!
! Arguments
!              

    real(r8), intent(in), dimension(ntot_gpt,pverp) :: CK1
    real(r8), intent(in), dimension(ntot_gpt,pverp) :: CK2
    real(r8), intent(in), dimension(ntot_gpt,pverp) :: CPB
    real(r8), intent(in), dimension(ntot_gpt,pverp) :: CMB
    real(r8), intent(in), dimension(ntot_gpt,pverp) :: EM1
    real(r8), intent(in), dimension(ntot_gpt,pverp) :: EM2
    real(r8), intent(in), dimension(ntot_gpt,pverp) :: EL1
    real(r8), intent(in), dimension(ntot_gpt,pverp) :: EL2
    real(r8), intent(in), dimension(ntot_gpt,pverp) :: DIRECT       
    real(r8), intent(in), dimension(ntot_gpt,pverp) :: DIRECTU   
    real(r8), intent(in), dimension(ntot_gpt,pverp) :: DIREC
    real(r8), intent(in), dimension(pver) ::  dzc          ! [kg m-2], column amount of mass 
    integer, intent(in) :: swcut
    logical, intent(in) :: part_in_tshadow
    logical, intent(in) :: sw_on
 
    real(r8), intent(out), dimension(pver) ::  sw_dTdt     
    real(r8), intent(out), dimension(pver) ::  lw_dTdt     
    real(r8), intent(out), dimension(pverp) ::  sw_upflux   
    real(r8), intent(out), dimension(pverp) ::  sw_dnflux   
    real(r8), intent(out), dimension(pverp) ::  lw_upflux   
    real(r8), intent(out), dimension(pverp) ::  lw_dnflux       
    real(r8), intent(out), dimension(pverp,ntot_wavlnrng) ::  sw_upflux_spectral   
    real(r8), intent(out), dimension(pverp,ntot_wavlnrng) ::  sw_dnflux_spectral   
    real(r8), intent(out), dimension(pverp,ntot_wavlnrng) ::  lw_upflux_spectral
    real(r8), intent(out), dimension(pverp,ntot_wavlnrng) ::  lw_dnflux_spectral
    real(r8), intent(out) ::  vis_dir       
    real(r8), intent(out) ::  vis_dif       
    real(r8), intent(out) ::  nir_dir       
    real(r8), intent(out) ::  nir_dif       


!------------------------------------------------------------------------
!
! Local Variables
!          
    
    integer :: k
    integer :: ip,ipc,iw
    integer :: k_1
    integer :: j, it, ig
    real(r8) :: lyr_mass_fact
    real(r8) :: fdiv_sw
    real(r8) :: fdiv_lw


!------------------------------------------------------------------------
!
! Start Code
!          

    sw_dTdt(:) = 0.    ! Initialize heating rate arrays
    lw_dTdt(:) = 0.    !
 
    lw_dnflux(:) = 0.   !
    lw_upflux(:) = 0.   ! Initialize entire arrays for summing below
    sw_upflux(:) = 0.   !
    sw_dnflux(:) = 0.   !

    lw_dnflux_spectral(:,:) = 0.   !
    lw_upflux_spectral(:,:) = 0.   ! Initialize entire arrays for summing below
    sw_upflux_spectral(:,:) = 0.   !
    sw_dnflux_spectral(:,:) = 0.   !

    vis_dir = 0.     !
    vis_dif = 0.     !  Initialize solar fluxes to surface
    nir_dir = 0.     !  Pass to land model in CESM
    nir_dif = 0.     !
    
    ! Finalize fluxes: 
    do k=camtop,pverp    ! Loop over all layer BOUNDARIES, was k=1

      !Loop for broadband longwave fluxes
      do ip=lw_ipbeg,lw_ipend
        lw_dnflux(k) = lw_dnflux(k)+DIREC(ip,k)
        lw_upflux(k) = lw_upflux(k)+DIRECTU(ip,k)
      enddo


      !Loop for spectral longwave fluxes
      ! THIS FEATURE NOT IN CESM VERSION
      it=0
      do iw=1,ntot_wavlnrng
        do ig=1,ngauss_pts(iw)
	  it=it+1
          lw_dnflux_spectral(k,iw) = lw_dnflux_spectral(k,iw)+DIREC(it,k)
          lw_upflux_spectral(k,iw) = lw_upflux_spectral(k,iw)+DIRECTU(it,k)
        enddo
      enddo

    enddo
   
    if (sw_on) then
      do k=camtop,pverp  ! Loop over all layer BOUNDARIES, was k=1

        !Loop for broadband shortwave fluxes
        do ip=sw_ipbeg,sw_ipend  ! Loop over all "shortwave" wavelength intervals                 
          sw_upflux(k) = sw_upflux(k)+CK1(ip,k)*EL1(ip,k)+  &
                         CK2(ip,k)*EM1(ip,k)+CPB(ip,k)
          sw_dnflux(k) = sw_dnflux(k)+CK1(ip,k)*EL2(ip,k)+  &
                         CK2(ip,k)*EM2(ip,k)+CMB(ip,k)+DIRECT(ip,k)
        enddo

        !Loop for spectral shortwave fluxes
        ! THIS FEATURE NOT IN CESM VERSION
        it=0
        do iw=1,ntot_wavlnrng
          do ig=1,ngauss_pts(iw)
	    it=it+1
            sw_upflux_spectral(k,iw) = sw_upflux_spectral(k,iw)+CK1(it,k)*EL1(it,k)+  &
                                       CK2(it,k)*EM1(it,k)+CPB(it,k)
            sw_dnflux_spectral(k,iw) = sw_dnflux_spectral(k,iw)+CK1(it,k)*EL2(it,k)+  &
                                       CK2(it,k)*EM2(it,k)+CMB(it,k)+DIRECT(it,k)
          enddo
        enddo

      enddo  
    endif
   
    if(sw_on) then      !***** BOTH "longwave" and "shortwave" *****

      ! Calculate atmospheric heating rates (dT/dt) [K/s]:
      ! Column sunlit above layer 'swcut'

      if(part_in_tshadow) then  ! NOTES: set to false always
       
        do k=camtop+1,swcut-1  ! Above shadow, was k=2
      
          lyr_mass_fact = dzc(k-1)*cpair
        
          fdiv_sw = (sw_upflux(k)-sw_dnflux(k))-(sw_upflux(k-1)-sw_dnflux(k-1))
          sw_dTdt(k-1) = fdiv_sw/lyr_mass_fact      ! "shortwave" heating rate [K/s]

          fdiv_lw = (lw_upflux(k)-lw_dnflux(k))-(lw_upflux(k-1)-lw_dnflux(k-1))
          lw_dTdt(k-1) = fdiv_lw/lyr_mass_fact       ! "longwave" heating rate [K/s]
          
        enddo

        do k=swcut,pverp  ! Within shadow, no shortwave calculation

          lyr_mass_fact = dzc(k-1)*cpair
      
          sw_dTdt(k-1) = 0.0   ! "shortwave" heating rate [K/s]
          
          fdiv_lw = (lw_upflux(k)-lw_dnflux(k))-(lw_upflux(k-1)-lw_dnflux(k-1))
          lw_dTdt(k-1) = fdiv_lw/lyr_mass_fact    ! "longwave" heating rate [K/s]

        enddo 

      else   ! No shadow in column
       
        do k=camtop+1,pverp  ! was k=2

          lyr_mass_fact = dzc(k-1)*cpair

          fdiv_sw = (sw_upflux(k)-sw_dnflux(k))-(sw_upflux(k-1)-sw_dnflux(k-1))
          sw_dTdt(k-1) = fdiv_sw/lyr_mass_fact      ! "shortwave" heating rate [K/s]

          fdiv_lw = (lw_upflux(k)-lw_dnflux(k))-(lw_upflux(k-1)-lw_dnflux(k-1))
          lw_dTdt(k-1) = fdiv_lw/lyr_mass_fact       ! "longwave" heating rate [K/s]
  
        enddo

      endif
 
    else   !***** ONLY "longwave" ***** 
     
      do k=camtop+1, pverp   ! was k=2

        lyr_mass_fact = dzc(k-1)*cpair
        fdiv_lw = (lw_upflux(k)-lw_dnflux(k))-(lw_upflux(k-1)-lw_dnflux(k-1))
        lw_dTdt(k-1) = fdiv_lw/lyr_mass_fact       ! "longwave" heating rate [K/s]
      
        !write(*,*) "--------------------------------------------------------"
        !write(*,*) k,"lw_up(k)  lw_up(k-1)  lw_dn(k)  lw_dn(k-1)"
        !write(*,*) lw_upflux(k), lw_upflux(k-1), lw_dnflux(k), lw_dnflux(k-1)
        !write(*,*) fdiv_lw,lyr_mass_fact
        !write(*,*) lw_dTdt(k-1)*24.0*60.*60.0

      enddo

    endif

    !if (camtop > 1) then
    !  lw_dTdt(:camtop-1) = 0.0
    !  sw_dTdt(:camtop-1) = 0.0
    !endif

    ! Calculate surface shortwave fluxes to land model
    ipc=sw_ipbeg 
    do iw=sw_iwbeg,sw_iwend    ! Loop over relevant wavelength intervals
      if (wavenum_edge(iw) .gt. 13000) then 
        do ip=1,ngauss_pts(iw)
          vis_dir = vis_dir + DIRECT(ipc,pverp)
          vis_dif = vis_dif + CK1(ipc,pverp)*EL2(ipc,pverp)+CK2(ipc,pverp)*EM2(ipc,pverp)+CMB(ipc,pverp)
          ipc=ipc+1
        enddo
      else
        do ip=1,ngauss_pts(iw)
          nir_dir = nir_dir + DIRECT(ipc,pverp)
          nir_dif = nir_dif + CK1(ipc,pverp)*EL2(ipc,pverp)+CK2(ipc,pverp)*EM2(ipc,pverp)+CMB(ipc,pverp)
          ipc=ipc+1
        enddo
      endif
    enddo

    return

  end subroutine rad_postcalc

!!============================================================================
!
!  function solarflux(wn_lw, wn_up, ib)
!
!!------------------------------------------------------------------------
!!
!! Purpose: Computes integrated solar intensity (flux) over given 
!!          wavenumber range
!!                                                                      
!!------------------------------------------------------------------------
!
!    implicit none
!
!!------------------------------------------------------------------------
!!
!! Input Arguments
!!          
!    real(r8), intent(in) :: wn_lw     ! wavenumber [cm-1]
!    real(r8), intent(in) :: wn_up     ! wavenumber [cm-1]
!    integer, intent(in) :: ib
!
!!------------------------------------------------------------------------
!!
!! Local Variables
!!          
!    real(r8) :: solarflux
!    real(r8) :: ftmp
!    real(r8) :: dwl
!    real(r8) :: wn
!    real(r8) :: solar
!    real(r8) :: J
!    integer :: itop
!    integer :: ibot
!    integer :: I
!    integer :: iband
!
!!------------------------------------------------------------------------
!!
!! Start Code
!!          
!    ftmp = 0.0
!    iband = 0
!
!    iband = iband + ib
!
!    do I = wn_lw, wn_up, 10
!      if(iband.eq.ntot_wavlnrng) iband = iband - 1
!      J=real(I)
!      if(I.eq.wn_up) then
!        dwl=1.e4/(J-10.)-1.e4/J
!      else
!        dwl=1.e4/J-1.e4/(J+10.)
!      endif
!      wn=I
!      ftmp = ftmp+solarrad(wn)*dwl 
!    enddo
!
!    solarflux = ftmp 
!
!    return
!
!  end function solarflux

!!============================================================================
!
!  function solarrad(w)
!
!!------------------------------------------------------------------------
!!
!! Purpose
!!         
!!------------------------------------------------------------------------
!
!    implicit none
!
!!------------------------------------------------------------------------
!!
!! Input Arguments
!!         
!    real(r8), intent(in) :: w   ! wavenumber cm-1
!
!!------------------------------------------------------------------------
!!
!! Local Variables
!!         
!    real(r8) :: wave
!    real(r8) :: solarrad
!    real(r8) :: i
!    real(r8) :: f
!    integer :: j
!
!    wave = w                  !wavenumber [cm^-1]
!    i=int((wave-wm1)/dwm+1.00001)
!    j=int(i)
!    solarrad=0.
!    if(i.lt.1 .or. i.gt.nw-1) return
!    f = (wave-wm1-dwm*(i-1))/dwm
!    solarrad=(1.-f)*sunm(j)+f*sunm(j+1)
!    return
!  end function solarrad



end module radiation
