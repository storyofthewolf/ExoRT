!===============================================================================
! CVS $Id: shr_const_mod.F90,v 1.2.4.1 2004/01/02 18:50:53 mvr Exp $
! CVS $Source: /fs/cgd/csm/models/CVS.REPOS/shared/csm_share/shr/shr_const_mod.F90,v $
! CVS $Name: cam3_1_9_brnch_waccm $
!===============================================================================

MODULE shr_const_mod

   use shr_kind_mod, only: SHR_KIND_R8

   !----------------------------------------------------------------------------
   ! physical constants (all data public)
   ! Values below are for the current Earth.
   ! Procedure in physconst.F90 changes RDAIR, RHODAIR, CPAIR, MWDAIR for variable 
   ! gas mixtures. 
   !----------------------------------------------------------------------------
   public
   real(SHR_KIND_R8),parameter :: SHR_CONST_PI     = 3.14159265358979323846_SHR_KIND_R8  ! pi
   real(SHR_KIND_R8),parameter :: SHR_CONST_TWOPI  = 2.0*SHR_CONST_PI         ! 2 pi
   real(SHR_KIND_R8),parameter :: SHR_CONST_CSEC   = 86400.0_SHR_KIND_R8      ! sec in calendar day ~ sec
   real(SHR_KIND_R8),parameter :: SHR_CONST_CSECR  = 86400.0                  ! sec in calendar day ~ sec
   integer,          parameter :: SHR_CONST_CSECI  = 86400                    ! sec in calendar day, integer type 
   real(SHR_KIND_R8),parameter :: SHR_CONST_SSEC   = 86163.3_SHR_KIND_R8      ! sec in siderial day ~ sec  = CSEC - CSEC/DAYS
   real(SHR_KIND_R8),parameter :: SHR_CONST_DAYS   = 365.0_SHR_KIND_R8        ! number of days in one year
   real(SHR_KIND_R8),parameter :: SHR_CONST_DAYSL  = 365.25_SHR_KIND_R8       ! number of days in one year (leap year calculations) 
   integer,          parameter :: SHR_CONST_DAYSI  = 365                      ! number of days in one year, integer type    
   real(SHR_KIND_R8),parameter :: SHR_CONST_CDVE   = 80.5                     ! calendar day of vernal equinox (March 21)
   ! Time/rotation rate scaling parameters
   real(SHR_KIND_R8),parameter :: DAY_RATIO = 1.0_SHR_KIND_R8
   real(SHR_KIND_R8),parameter :: YEAR_RATIO = 1.0_SHR_KIND_R8
   !
   real(SHR_KIND_R8),parameter :: SHR_CONST_OMEGA  = 2.0_SHR_KIND_R8*SHR_CONST_PI/SHR_CONST_SSEC/DAY_RATIO ! earth rot ~ rad/sec  
   real(SHR_KIND_R8),parameter :: SHR_CONST_REARTH = 6.37122e6_SHR_KIND_R8    ! radius of earth ~ m
   real(SHR_KIND_R8),parameter :: SHR_CONST_G      = 9.80616_SHR_KIND_R8      ! acceleration of gravity ~ m/s^2
     !real(SHR_KIND_R8),parameter :: SHR_CONST_PSTD   = 101325.0_SHR_KIND_R8     ! standard pressure ~ pascals
   ! standard pressure ~ pascals
    !real(SHR_KIND_R8),parameter :: SHR_CONST_PSTD   = (aN2bar + aARbar + aCO2bar + aCH4bar) * 1.0e5
    !real(SHR_KIND_R8),parameter :: SHR_CONST_PREF   = 100000.0_SHR_KIND_R8    ! reference pressure
    !real(SHR_KIND_R8),parameter :: SHR_CONST_PREF   = SHR_CONST_PSTD           ! reference pressure
   real(SHR_KIND_R8),parameter :: SHR_CONST_MSDIST2= 1.0                      ! semimajor axis squared ~ AU^2
   real(SHR_KIND_R8),parameter :: SHR_CONST_STEBOL = 5.67e-8_SHR_KIND_R8      ! Stefan-Boltzmann constant ~ W/m^2/K^4
   real(SHR_KIND_R8),parameter :: SHR_CONST_BOLTZ  = 1.38065e-23_SHR_KIND_R8  ! Boltzmann's constant ~ J/K/molecule
   real(SHR_KIND_R8),parameter :: SHR_CONST_AVOGAD = 6.02214e26_SHR_KIND_R8   ! Avogadro's number ~ molecules/kmole
   real(SHR_KIND_R8),parameter :: SHR_CONST_RGAS   = SHR_CONST_AVOGAD*SHR_CONST_BOLTZ ! Universal gas constant ~ J/K/kmole
   real(SHR_KIND_R8),parameter :: SHR_CONST_MWDAIR = 28.966_SHR_KIND_R8       ! molecular weight dry air ~ kg/kmole
   real(SHR_KIND_R8),parameter :: SHR_CONST_MWWV   = 18.016_SHR_KIND_R8       ! molecular weight water vapor
   real(SHR_KIND_R8),parameter :: SHR_CONST_RDAIR  = SHR_CONST_RGAS/SHR_CONST_MWDAIR  ! Dry air gas constant ~ J/K/kg
   real(SHR_KIND_R8),parameter :: SHR_CONST_RWV    = SHR_CONST_RGAS/SHR_CONST_MWWV    ! Water vapor gas constant ~ J/K/kg
   real(SHR_KIND_R8),parameter :: SHR_CONST_ZVIR   = (SHR_CONST_RWV/SHR_CONST_RDAIR)-1.0_SHR_KIND_R8   ! RWV/RDAIR - 1.0
   real(SHR_KIND_R8),parameter :: SHR_CONST_KARMAN = 0.4_SHR_KIND_R8          ! Von Karman constant
   real(SHR_KIND_R8),parameter :: SHR_CONST_LOSCHMIDT = 2.6867774e25          ! Loschmidt number ~# m^-3
   real(SHR_KIND_R8),parameter :: SHR_CONST_TKFRZ  = 273.16_SHR_KIND_R8       ! freezing T of fresh water ~ K (intentionally made == to TKTRIP)
   real(SHR_KIND_R8),parameter :: SHR_CONST_TKTRIP = 273.16_SHR_KIND_R8       ! triple point of fresh water ~ K
   ! real(SHR_KIND_R8),parameter :: SHR_CONST_RHODAIR=SHR_CONST_PSTD/ &
   !   (SHR_CONST_RDAIR*SHR_CONST_TKFRZ)         ! density of dry air at STP   ~ kg/m^3
   real(SHR_KIND_R8),parameter :: SHR_CONST_RHOFW  = 1.000e3_SHR_KIND_R8      ! density of fresh water ~ kg/m^3
   real(SHR_KIND_R8),parameter :: SHR_CONST_RHOSW  = 1.026e3_SHR_KIND_R8      ! density of sea water ~ kg/m^3
   real(SHR_KIND_R8),parameter :: SHR_CONST_RHOICE = 0.917e3_SHR_KIND_R8      ! density of ice   ~ kg/m^3
   real(SHR_KIND_R8),parameter :: SHR_CONST_CPDAIR = 1.00464e3_SHR_KIND_R8    ! specific heat of dry air ~ J/kg K
   real(SHR_KIND_R8),parameter :: SHR_CONST_CPFW   = 4.188e3_SHR_KIND_R8      ! specific heat of fresh h2o ~ J/kg K
   real(SHR_KIND_R8),parameter :: SHR_CONST_CPSW   = 3.996e3_SHR_KIND_R8      ! specific heat of sea h2o ~ J/kg K
   real(SHR_KIND_R8),parameter :: SHR_CONST_CPWV   = 1.810e3_SHR_KIND_R8      ! specific heat of water vap ~ J/kg K
   real(SHR_KIND_R8),parameter :: SHR_CONST_CPICE  = 2.11727e3_SHR_KIND_R8    ! specific heat of fresh ice ~ J/kg K
   real(SHR_KIND_R8),parameter :: SHR_CONST_CPN2   = 1.039e3_SHR_KIND_R8      ! specific heat of N2 gas  ~ J/kg K
   real(SHR_KIND_R8),parameter :: SHR_CONST_CPCO2  = 0.846e3_SHR_KIND_R8      ! specific heat of CO2 gas ~ J/kg K
   real(SHR_KIND_R8),parameter :: SHR_CONST_CPCH4  = 2.226e3_SHR_KIND_R8      ! specific heat of CH4 gas ~ J/kg K
   real(SHR_KIND_R8),parameter :: SHR_CONST_CPAR   = 0.520e3_SHR_KIND_R8      ! specific heat of Ar gas ~ J/kg K
   real(SHR_KIND_R8),parameter :: SHR_CONST_CPH2   = 14.30e3_SHR_KIND_R8      ! specific heat of H2 gas ~ J/kg K
   

   real(SHR_KIND_R8),parameter :: SHR_CONST_LATICE = 3.337e5_SHR_KIND_R8      ! latent heat of fusion ~ J/kg
   real(SHR_KIND_R8),parameter :: SHR_CONST_LATVAP = 2.501e6_SHR_KIND_R8      ! latent heat of evaporation ~ J/kg
   real(SHR_KIND_R8),parameter :: SHR_CONST_LATSUB = SHR_CONST_LATICE + SHR_CONST_LATVAP ! latent heat of sublimation ~ J/kg
   real(SHR_KIND_R8),parameter :: SHR_CONST_OCN_REF_SAL = 34.7_SHR_KIND_R8    ! ocn ref salinity (psu)
   real(SHR_KIND_R8),parameter :: SHR_CONST_ICE_REF_SAL =  4.0_SHR_KIND_R8    ! ice ref salinity (psu)

   real(SHR_KIND_R8),parameter :: SHR_CONST_SPVAL       = 1.0e30_SHR_KIND_R8  ! special missing value



END MODULE shr_const_mod
