module physconst
!
! Physical constants (replaces /comcon/)
!
   use shr_kind_mod,      only: r8 => shr_kind_r8
   use infnan,            only: inf
   use input
   use shr_const_mod,     only: shr_const_g,      shr_const_stebol, shr_const_tkfrz,  &
                                shr_const_mwdair, shr_const_rdair,  shr_const_mwwv,   &
                                shr_const_latice, shr_const_latvap, shr_const_cpdair, &
                                shr_const_rhofw,  shr_const_cpwv,   shr_const_rgas,   &
                                shr_const_karman, &  !shr_const_pstd,   shr_const_rhodair,&
                                shr_const_avogad, shr_const_boltz,  shr_const_cpn2,   &
                                shr_const_cpco2,  shr_const_cpch4,  shr_const_cpar,   &
                                shr_const_cph2,   shr_const_tkfrz

   implicit none
!
! Make module data private by default (so that modules used don't become public)
!
   private
   save

   public :: physconst_defaultopts
   public :: physconst_setopts
   public :: physconst_setgas
 
! Universal constants
   real(r8), public, parameter :: r_universal = shr_const_rgas ! Universal gas constant (J/K/kmol)
   real(r8), public, parameter :: stebol = shr_const_stebol    ! Stefan-Boltzmann's constant

! Constants for Earth
   real(r8), public, parameter :: gravit = shr_const_g      ! gravitational acceleration
   real(r8), public, parameter :: rga = 1./gravit           ! reciprocal of gravit

! Specific Heats
   real(r8), public, parameter :: cpn2 = shr_const_cpn2     ! specific heat, dinitrogen gas
   real(r8), public, parameter :: cpco2 = shr_const_cpco2   ! specific heat, carbon dioxide gas
   real(r8), public, parameter :: cpch4 = shr_const_cpch4   ! specific heat, methane gas
   real(r8), public, parameter :: cpar = shr_const_cpar     ! specific heat, argon gas
   real(r8), public, parameter :: cph2 = shr_const_cph2     ! specific heat of H2

! Standard Pressure
   !real(r8), public, parameter :: pstd = shr_const_pstd     ! Standard pressure Pascals

! Molecular weights
   real(r8), public, parameter :: mwn2  =  28.              ! molecular weight of n2
   real(r8), public, parameter :: mwar  =  40.              ! molecular weight of Ar
   real(r8), public, parameter :: mwco2 =  44.              ! molecular weight co2
   real(r8), public, parameter :: mwh2o =  shr_const_mwwv   ! molecular weight h2o
   real(r8), public, parameter :: mwn2o =  44.              ! molecular weight n2o
   real(r8), public, parameter :: mwch4 =  16.              ! molecular weight ch4
   real(r8), public, parameter :: mwf11 = 136.              ! molecular weight cfc11
   real(r8), public, parameter :: mwf12 = 120.              ! molecular weight cfc12
   real(r8), public, parameter :: mwo3  =  48.              ! molecular weight O3
   real(r8), public, parameter :: mwo2  =  32.		    ! molecular weight O2
   real(r8), public, parameter :: mwh2  =  2.		    ! molecular weight H2   
   
   !real(r8), public :: mwdry = inf                         ! molecular weight of dry air

! Constants dependent on gas mixture
   !real(r8), public, parameter :: mwdry =  shr_const_mwdair
   real(r8), public :: mwdry 
   !real(r8), public, parameter :: &   ! molecular weight of archean air
   !     mwdry = aN2vmr*mwn2 + aCO2vmr*mwco2 + aCH4vmr*mwch4  + aARvmr*mwar !+ aH2vmr*mwh2
   real(r8), public :: cpair 
   !real(r8), public, parameter :: cpair = shr_const_cpdair
   !real(r8), public, parameter :: &
   !     cpair = aN2vmr*cpn2 + aCO2vmr*cpco2 + aCH4vmr*cpch4 + aARvmr*cpar !+ aH2vmr*cph2
    

   !real(r8), public, parameter :: rair = r_universal/mwdry
   !real(r8), public, parameter :: cappa = rair/cpair
   !real(r8), public, parameter :: rhodair = pstd/(rair*shr_const_tkfrz)

   !real(r8), public :: rair = inf                                 ! Gas constant for dry air (J/K/kg)
   !real(r8), public :: cpair = inf                                ! specific heat of dry air (J/K/kg)
   !real(r8), public :: cappa = inf                                ! R/Cp
   !real(r8), public :: rhodair = inf                              ! density of dry air at STP (kg/m3)

   ! real(r8), public, parameter :: rair = shr_const_rdair
   ! real(r8), public, parameter :: cpair = shr_const_cpdair 
   ! real(r8), public, parameter :: cappa = shr_const_rdair/shr_const_cpdair       
   ! real(r8), public, parameter :: rhodair = shr_const_rhodair 

! Constants for water
   real(r8), public, parameter :: cph2o = shr_const_cpwv    ! specific heat of water vapor (J/K/kg)
   real(r8), public, parameter :: latvap = shr_const_latvap ! Latent heat of vaporization
   real(r8), public, parameter :: latice = shr_const_latice ! Latent heat of fusion
   real(r8), public, parameter :: rhoh2o = shr_const_rhofw  ! Density of liquid water (STP)
   real(r8), public, parameter :: tmelt = shr_const_tkfrz   ! Freezing point of water
   real(r8), public, parameter :: epsilo = shr_const_mwwv/shr_const_mwdair ! ratio of h2o to dry air molecular weights    

   !real(r8), public, parameter :: epsilo = shr_const_mwwv/mwdry
   !real(r8), public :: epsilo = inf

! Other molecular constants
   real(r8), public, parameter :: avogad = shr_const_avogad ! Avogadro's number
   real(r8), public, parameter :: boltz  = shr_const_boltz  ! Boltzman's constant

! Turbulence
   real(r8), public, parameter :: karman = shr_const_karman ! VonKarman constant

!
! Defer setting the following variables to their actual values because their value depends on
! whether the run is adiabatic or not
!
   real(r8), public :: rh2o = inf          ! gas constant for water vapor
   real(r8), public :: cpvir = inf         ! cpwv/cpair - 1
   real(r8), public :: cpwv = inf          ! Specific heat of water vapor
   real(r8), public :: zvir = inf          ! rh2o/rair - 1

! Private Data
   logical :: use_archean_constants = .false.
   logical :: no_clm = .false.  !default is to use clm model

!======================================================================
   contains
!======================================================================
   
   subroutine physconst_defaultopts(use_archean_constants_out, no_clm_out)

     implicit none
     logical, intent(out), optional :: use_archean_constants_out
     logical, intent(out), optional :: no_clm_out

     if ( present(use_archean_constants_out) ) use_archean_constants_out = use_archean_constants
     if ( present(no_clm_out) ) no_clm_out = no_clm

   end subroutine physconst_defaultopts

!====================================================

   subroutine physconst_setopts(use_archean_constants_in, no_clm_in)

     implicit none
     logical, intent(in), optional :: use_archean_constants_in
     logical, intent(in), optional :: no_clm_in

     if ( present(use_archean_constants_in) ) use_archean_constants = use_archean_constants_in
     if ( present(no_clm_in) ) no_clm = no_clm_in

   end subroutine physconst_setopts
 
!====================================================

   subroutine physconst_setgas(mwin, cpin )
     real(r8), intent(in) :: mwin
     real(r8), intent(in) :: cpin

     mwdry = mwin
     cpair = cpin
   end subroutine physconst_setgas


end module physconst