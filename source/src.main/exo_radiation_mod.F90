
module exo_radiation_mod

!---------------------------------------------------------------------
! Purpose:
!
! Contains radiative transfer algorithm for ExoRT
!
! Revision history
! September 2010, E. T. Wolf, R. Urata CAM3
! March     2014, E. T. Wolf --- decoupled solar and IR streams
!                            --- merged with CESM1.2.1
! February  2015, E.T. Wolf  --- new build
!---------------------------------------------------------------------
!
  use shr_kind_mod,     only: r8 => shr_kind_r8
  use shr_const_mod,    only: SHR_CONST_PI,SHR_CONST_G, SHR_CONST_PSTD, &
                              SHR_CONST_RGAS, SHR_CONST_AVOGAD, &
                              SHR_CONST_STEBOL, SHR_CONST_CDAY, &
                              SHR_CONST_BOLTZ, &
                              SHR_CONST_RHOFW, SHR_CONST_RHOICE, &
			      SHR_CONST_LOSCHMIDT
  use physconst,        only: scon,mwn2, mwco2, mwch4, mwc2h6, mwh2o, mwo2, mwh2, mwdry, cpair, cappa
  use ppgrid            ! pver, pverp is here
  use pmgrid            ! ?masterproc is here?
  use spmd_utils,       only: masterproc
  use rad_interp_mod
  use radgrid
  use kabs
  use time_manager,     only: get_nstep
  use calc_opd_mod
  use exo_init_ref
  use planck_mod
  use exoplanet_mod,    only: do_exo_clouds, do_exo_haze

  implicit none
  private
  save

!------------------------------------------------------------------------
!
! Public interfaces
!
  public :: init_planck
  public :: aerad_driver

!------------------------------------------------------------------------
!
! private data
!

  ! Default values for namelist variables

  integer :: openstatus

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


!============================================================================
  ! fused solver: current spectral block (global g-point bounds and
  ! band range) and the first g-point of each band. Per thread (OpenMP).
  integer :: gb0 = 1, gb1 = 1, bw0 = 1, bw1 = 1
  integer, dimension(ntot_wavlnrng+1) :: gstart = 1
  ! Target g-points per spectral block: the solver's block work arrays
  ! (~54 x nblk_gpt x pverp doubles) stay cache-resident. 8-32 measured within
  ! 10% of each other on Milan and Apple M5; 16 = 2 bands of 8 g-points.
  integer, parameter :: nblk_gpt = 16
!$omp threadprivate(gb0, gb1, bw0, bw1, gstart)

contains
!============================================================================

!============================================================================
!
! Public subroutines
!
!============================================================================

   subroutine init_planck
!------------------------------------------------------------------------
!
! Purpose: Initial reference value arrays
!
!------------------------------------------------------------------------

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

!    real(r8) :: wn1
!    real(r8) :: wn2
!    real(r8) :: tempsol
!    real(r8) :: tempsol_iw
    real(r8) :: w1
    real(r8) :: w2
    real(r8) :: tempk
 !   logical :: found

!------------------------------------------------------------------------
!
! Start Code
!

    !
    ! Calculate Planck function quantities as a function of wavelength point and
    !  temperature:
    !

    if (masterproc) then
      write(*,*)   "INIT_PLANCK: CREATING PLANCK FUNCTION TABLE"
    endif
    ip = lw_ipbeg-1

    do iw=1,ntot_wavlnrng   ! Loop over all wavelength intervals

      ! Set these so that difference below is
      !   f(longer_wavlen)-f(shorter_wavlen),
      !   w1 = longer_wavlen, w2 = shorter_wavlen:
      !   h*c*nu/kc
      w1 = dble(1.439*wavenum_edge(iw))      ! 1.439 ~ ((h*c)/k) ,
      w2 = dble(1.439*wavenum_edge(iw+1))    ! convert nu cm-1 to lambda (um

      do ig=1,ngauss_pts(iw)
        ip = ip+1

        do it=tpft_beg,tpft_end
          tempk = tpft_invfinc*dble(it)   ! Temperature [K]

          !ptemp_itp, spectrally integrated radiance (isotropic) in each longwave guass interval
          !value weighted according to gauss weights, (B outside mapping)
          ptemp_itp(it,ip) = (PLANCKf(w1,tempk)-PLANCKf(w2,tempk))* &
                              g_weight(ip)*SHR_CONST_STEBOL/SHR_CONST_PI
          !write(*,*) "PLANCK", (PLANCKf(w1,tempk)-PLANCKf(w2,tempk)), g_weight(ip)
        enddo
        ptemp_itp(1:tpft_beg,ip) = ptemp_itp(tpft_beg,ip)
      enddo
    enddo

    return

  end subroutine init_planck


!============================================================================

  subroutine aerad_driver(ext_H2O, ext_CO2, &
                          ext_CH4, ext_C2H6, &
                          ext_NH3, ext_CO, &
                          ext_H2, ext_N2, ext_O3, ext_O2, &
                          ext_sfcT, ext_sfcP, ext_pmid, &
                          ext_pdel, ext_pdeldry, ext_tmid, ext_pint, ext_pintdry, &
                          ext_cosZ, ext_msdist, &
                          ext_asdir, ext_aldir, &
                          ext_asdif, ext_aldif,  &
                          ext_rtgt, ext_solar_azm_ang, ext_tazm_ang, ext_tslope_ang,  &
                          ext_tslas_tog, ext_tshadow_tog, ext_nazm_tshadow, ext_cosz_horizon,  &
                          ext_TCx_obstruct, ext_TCz_obstruct, ext_zint, &
                          sw_dTdt, lw_dTdt, lw_dnflux, lw_upflux, sw_upflux, sw_dnflux, &
                          lw_dnflux_spec, lw_upflux_spec, sw_upflux_spec, sw_dnflux_spec, &
                          vis_dir, vis_dif, nir_dir, nir_dif, sol_toa, &
                          ext_cicewp, ext_cliqwp, ext_cfrc, &
                          ext_rei, ext_rel, &
                          ext_cicewp_co2, ext_rei_co2, &
                          ext_carmammr, &
                          ext_srf_emiss, &
                          ext_mwdry, ext_cpdry, ext_icol, ext_grav )


!------------------------------------------------------------------------
!
! Purpose: Driver for correlated K radiative transfer code.
!          Recieves column data from radiation_tend
!
! Argument contract:
!   Mandatory (positional): thermodynamic state, geometry, surface
!   albedos, and all 10 gas mass mixing ratios. The off switch for a gas
!   is a zero MMR (calc_opd_gas short-circuits zero-abundance gases),
!   never an absent argument.
!
!   Optional keyword tail (everything after sol_toa): condensed-phase and
!   situational inputs -- H2O clouds (cicewp/cliqwp/cfrc/rei/rel, all five
!   or none), CO2 ice clouds (cicewp_co2/rei_co2, both or neither), CARMA
!   haze (carmammr), surface emissivity (srf_emiss), per-column dry-air
!   properties (mwdry/cpdry; absent -> the physconst module values, which
!   is the CAM path where CAM owns them), a 1-based batch column
!   index (icol; used ONLY to offset the MCICA permute seed so batch
!   columns decorrelate — absent or 1 reproduces the legacy constant
!   seed), and per-column surface gravity (grav; absent -> the
!   SHR_CONST_G module value = the CAM path). Per-column insolation
!   needs no new argument: ext_msdist already scales the TOA flux
!   (sflux_frac = 1/msdist), so callers map a per-column stellar
!   constant onto it. Callers MUST pass these by
!   keyword and new optional arguments MUST be appended at the end, so
!   existing positional call sites (1-D, library, 3dmodels, ExoCAM
!   SourceMods) never break when physics is added.
!
!   Thread-safety: with mwdry/cpdry passed as arguments, a column solve
!   reads NO module-scope mutable state (all tables are read-only after
!   init). Its full-size optics arrays are per-thread (SAVEd + OpenMP
!   threadprivate) and reused across calls; in a build without OpenMP they
!   are shared, so do not call aerad_driver from concurrent threads in such
!   a build. The solver's block work arrays are automatic (stack): a calling
!   thread needs ~3 MB of stack at pver = 300. The main thread and
!   libgomp's default worker stacks suffice (multicol_check, 8 threads,
!   pver = 300); a caller-created thread with a smaller stack needs more.
!
!   Semantics per input: runtime flag on + argument present -> physics
!   active; flag on + argument absent -> zero contribution (the kernel is
!   not called; this is how a clear-sky companion call works); flag off
!   (do_exo_clouds / do_exo_haze) -> argument ignored, kernel never
!   called.
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
    real(r8), intent(in) :: ext_asdir          ! direct albedo (0.2-0.7 um) (from cam_in%asdir)
    real(r8), intent(in) :: ext_aldir          ! direct albedo (0.7-4.0 um) (from cam_in%aldir)
    real(r8), intent(in) :: ext_asdif          ! diffuse albedo (0.2-0.7 um) (from cam_in%asdif)
    real(r8), intent(in) :: ext_aldif          ! diffuse albedo (0.7-4.0 um) (from cam_in%aldif)
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
    real(r8), intent(in), dimension(pver) :: ext_pdeldry   ! layer thickness dry atmosphere [Pa]
    real(r8), intent(in), dimension(pver) :: ext_tmid      ! temperature at midpoints [K]
    real(r8), intent(in), dimension(pverp) :: ext_pint     ! wet pressure at interfaces
    real(r8), intent(in), dimension(pverp) :: ext_pintdry     ! dry pressure at interfaces

    real(r8), intent(in), dimension(pver) :: ext_H2O       ! specific humidy (from state%q < q) at midlayer [kg/kg]
    real(r8), intent(in), dimension(pver) :: ext_CO2       ! CO2 mass mixing ratio from state%q < co2mmr)   [kg/kg]
    real(r8), intent(in), dimension(pver) :: ext_CH4       ! CH4 mass mixing ratio from state%q < ch4mmr)   [kg/kg]
    real(r8), intent(in), dimension(pver) :: ext_C2H6      ! C2H6 mass mixing ratio from state%q < c2h6mmr)   [kg/kg]
    real(r8), intent(in), dimension(pver) :: ext_NH3       ! NH3 mass mixing ratio from state%q < nh3mmr)   [kg/kg]
    real(r8), intent(in), dimension(pver) :: ext_CO        ! CO mass mixing ratio from state%q < commr)   [kg/kg]
    real(r8), intent(in), dimension(pver) :: ext_H2        ! H2 mass mixing ratio from state%q < h2mmr)   [kg/kg]
    real(r8), intent(in), dimension(pver) :: ext_N2        ! N2 mass mixing ratio from state%q < h2mmr)   [kg/kg]
    real(r8), intent(in), dimension(pver) :: ext_O3        ! O3 mass mixing ratio from state%q < h2mmr)   [kg/kg]
    real(r8), intent(in), dimension(pver) :: ext_O2        ! O2 mass mixing ratio from state%q < h2mmr)   [kg/kg]

    ! {intent IN, optional keyword tail} -- see argument contract above
    real(r8), intent(in), dimension(pver), optional :: ext_cicewp    ! in cloud ice water path at layer midpoints [g/m2]
    real(r8), intent(in), dimension(pver), optional :: ext_cliqwp    ! in cloud liquid water path at layer midpoints [g/m2]
    real(r8), intent(in), dimension(pver), optional :: ext_cFRC      ! cloud fraction]
    real(r8), intent(in), dimension(pver), optional :: ext_rei       ! ice cloud particle effective drop size ice [microns]
    real(r8), intent(in), dimension(pver), optional :: ext_rel       ! liquid cloud drop effective drop size liquid [micron
    real(r8), intent(in), dimension(pver), optional :: ext_cicewp_co2  ! CO2 ice cloud water path at layer midpoints [g/m2]
    real(r8), intent(in), dimension(pver), optional :: ext_rei_co2     ! CO2 ice cloud particle effective radius [microns]
    real(r8), intent(in), dimension(pver,nelem_carma,nbin_carma), optional :: ext_carmammr  ! CARMA haze binwise mass mixing ratio [kg/kg]
    real(r8), intent(in), optional :: ext_srf_emiss  ! broadband surface thermal emissivity (default 1.0 if absent)
    real(r8), intent(in), optional :: ext_mwdry      ! dry-air molecular weight [g/mol] (absent -> physconst mwdry)
    real(r8), intent(in), optional :: ext_cpdry      ! dry-air specific heat [J/kg/K]   (absent -> physconst cpair)
    integer,  intent(in), optional :: ext_icol       ! 1-based batch column index; offsets the MCICA seed (absent or 1 -> legacy constant seed)
    real(r8), intent(in), optional :: ext_grav       ! surface gravity [m s-2] (absent -> SHR_CONST_G, the CAM path)

    real(r8), intent(out), dimension(pver) ::  sw_dTdt
    real(r8), intent(out), dimension(pver) ::  lw_dTdt

    real(r8), intent(out), dimension(pverp) ::  sw_upflux
    real(r8), intent(out), dimension(pverp) ::  sw_dnflux
    real(r8), intent(out), dimension(pverp) ::  lw_upflux
    real(r8), intent(out), dimension(pverp) ::  lw_dnflux

    real(r8), intent(out), dimension(pverp,ntot_wavlnrng) ::  sw_upflux_spec
    real(r8), intent(out), dimension(pverp,ntot_wavlnrng) ::  sw_dnflux_spec
    real(r8), intent(out), dimension(pverp,ntot_wavlnrng) ::  lw_upflux_spec
    real(r8), intent(out), dimension(pverp,ntot_wavlnrng) ::  lw_dnflux_spec

    real(r8), intent(out) ::  vis_dir
    real(r8), intent(out) ::  vis_dif
    real(r8), intent(out) ::  nir_dir
    real(r8), intent(out) ::  nir_dif
    real(r8), intent(out) ::  sol_toa


!------------------------------------------------------------------------
!
! Local Variables
!
     real(r8), dimension(pverp) :: coldens      ! [molec m-2] wet columen amount per layer
     real(r8), dimension(pverp) :: coldens_dry  ! [molec m-2] dry columen amount per layer
     real(r8), dimension(pverp) :: qH2O         ! [kg/kg] H2O  mass mixing ratio mid layers
     real(r8), dimension(pverp) :: qCO2         ! [kg/kg] CO2 mass mixing ratio at mid layers
     real(r8), dimension(pverp) :: qCH4         ! [kg/kg] CH4 mass mixing ratio at mid layers
     real(r8), dimension(pverp) :: qC2H6        ! [kg/kg] C2H6 mass mixing ratio at mid layers
     real(r8), dimension(pverp) :: qNH3         ! [kg/kg] NH3 mass mixing ratio at mid layers
     real(r8), dimension(pverp) :: qCO          ! [kg/kg] CO mass mixing ratio at mid layers
     real(r8), dimension(pverp) :: qO2          ! [kg/kg] O2 mass mixing ratio at mid layers
     real(r8), dimension(pverp) :: qO3          ! [kg/kg] O3 mass mixing ratio at mid layers
     real(r8), dimension(pverp) :: qH2          ! [kg/kg] H2 mass mixing ratio at mid layers
     real(r8), dimension(pverp) :: qN2          ! [kg/kg] H2 mass mixing ratio at mid layers
     real(r8), dimension(pverp) :: cICE         ! [g/m2] in cloud ice water path at mid layers
     real(r8), dimension(pverp) :: cLIQ         ! [g/m2] in cloud liquid water path at mid layers
     real(r8), dimension(pverp) :: cfrc         ! cloud fraction at mid layers
     real(r8), dimension(pverp) :: REI          ! [microns] ice cloud particle effective radii at mid layers
     real(r8), dimension(pverp) :: REL          ! [microns] liquid cloud drop effective radii at mid layers
     real(r8), dimension(pverp) :: cICE_co2     ! [g/m2] CO2 ice cloud water path at mid layers
     real(r8), dimension(pverp) :: REI_co2      ! [microns] CO2 ice cloud particle effective radii at mid layers
     real(r8), dimension(pverp,nelem_carma,nbin_carma) :: qcarma  ! [kg/kg] CARMA haze binwise mass mixing ratio at mid layers
     real(r8), dimension(pverp) :: masspath     ! [kg m-2] air mass path of each radiative layer
     real(r8), dimension(pverp) :: zlayer       ! [m] thickness of each vertical layer

     logical  :: have_cld_h2o  ! all five H2O-cloud optional args supplied
     logical  :: have_cld_co2  ! both CO2-cloud optional args supplied

     integer  :: mcica_seed    ! MCICA permute seed (nstep, offset by ext_icol-1 when supplied)

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
     integer :: iwb, iwlo, iwhi, nbw   ! spectral block loop

     logical :: found
     logical :: horizon_extension

     ! local variables used for computation (see Toon, 1989) annotate?
     !
     ! Two-stream solver work arrays live in solve_block, sized to one
     ! spectral block (nblk_gpt g-points) as automatic arrays: ~54 x nblk_gpt
     ! x pverp doubles plus add_txrad/refine_lwflux scratch, about 3 MB of
     ! thread stack at pver = 300 (0.5 MB at pver = 51).
     !
     ! The full-size optics arrays below stay on the heap (allocatable, so OpenMP
     ! workers and library caller threads run on default stacks, as in E2)
     ! but are SAVEd and OpenMP threadprivate: each thread allocates its own
     ! copy on its first call and reuses it afterwards. Allocating and freeing
     ! them on every call (~30 MB per column at pver=70, with the subroutines'
     ! arrays) can cost as much as the radiative transfer itself: when the
     ! freed block is at the top of the heap, glibc returns it to the OS and
     ! the next call page-faults it back in (ExoColumn: ~60% of its run time).
     ! Every call still initializes them as before (zero-filled below, or
     ! fully overwritten before use), so nothing carries over from one column
     ! to the next. Without OpenMP the directive is a comment (see
     ! Thread-safety above).
     real(r8), dimension(:,:), allocatable, save :: tau_gas       ! gas optical depth array
     real(r8), dimension(:,:), allocatable, save :: tau_ray  ! rayleigh optical depth

     ! Cloud optical properties
     real(r8), dimension(:,:,:), allocatable, save :: singscat_cld_mcica
     real(r8), dimension(:,:,:), allocatable, save :: asym_cld_mcica
     real(r8), dimension(:,:,:), allocatable, save :: tau_cld_mcica

     ! CO2 ice cloud optics, band-indexed (cloud fraction = 1, no MCICA)
     real(r8), dimension(:,:), allocatable, save :: tau_cld_co2
     real(r8), dimension(:,:), allocatable, save :: singscat_cld_co2
     real(r8), dimension(:,:), allocatable, save :: asym_cld_co2

     ! CARMA haze aerosol optics, band-indexed, presummed over elements/bins
     real(r8), dimension(:,:), allocatable, save :: tau_aer    ! sum of tau
     real(r8), dimension(:,:), allocatable, save :: wtau_aer   ! sum of w*tau
     real(r8), dimension(:,:), allocatable, save :: gwtau_aer  ! sum of g*w*tau

     ! stochastic bulk cloud properties (MCICA)
     real(r8), dimension(:,:), allocatable, save :: cFRC_mcica
     real(r8), dimension(:,:), allocatable, save :: cICE_mcica
     real(r8), dimension(:,:), allocatable, save :: cLIQ_mcica

     ! Planck function items

!$omp threadprivate(tau_gas, tau_ray, singscat_cld_mcica, asym_cld_mcica, &
!$omp               tau_cld_mcica, tau_cld_co2, singscat_cld_co2, asym_cld_co2, &
!$omp               tau_aer, wtau_aer, gwtau_aer, cFRC_mcica, cICE_mcica, cLIQ_mcica)

     ! albedo and emissivity
     real(r8) :: srf_emiss_val                   ! broadband surface emissivity (1.0 unless supplied)
     real(r8) :: mwdry_col                       ! per-column dry-air molecular weight [g/mol]
     real(r8) :: cpair_col                       ! per-column dry-air specific heat [J/kg/K]
     real(r8) :: grav_col                        ! per-column surface gravity [m s-2]
     real(r8), dimension(ntot_wavlnrng) :: sfc_albedo_dir   ! Surface albedo, direct radiation, spectral interval grid
     real(r8), dimension(ntot_wavlnrng) :: sfc_albedo_dif   ! Surface albedo, diffuse radiation, spectral interval grid
     real(r8), dimension(ntot_wavlnrng) :: sfc_emiss        ! Surface emissivity, spectral interval grid

     logical  :: part_in_tshadow
     real(r8) :: sflux_frac
     real(r8) :: sfc_tempk
     real(r8) :: sfc_press
     real(r8) :: cos_mu
     logical  :: sw_on     ! switch for togopography and day/night
     logical  :: beamSolar ! switch for beam

     real(r8), dimension(pver) :: dzc          ! [kg m-2], column amount of mass
     real(r8), dimension(pverp) :: tint        ! [K] temperatures at level interfaces
     real(r8), dimension(pverp) :: tmid        ! [K] temperatures at level at mid layers + top (isothermal)
     real(r8), dimension(pverp) :: pmid        ! [Pa] pressure at level at mid layers + top (isothermal)

     real(r8) :: dy

!------------------------------------------------------------------------
!
! Start Code
    ! Allocate this thread's work arrays on its first call only (see above)
    if (.not. allocated(tau_gas)) then
      allocate( tau_gas(ntot_gpt,pverp) )
      allocate( tau_ray(ntot_wavlnrng,pverp) )
      allocate( singscat_cld_mcica(ncld_grp,ntot_gpt,pverp) )
      allocate( asym_cld_mcica(ncld_grp,ntot_gpt,pverp) )
      allocate( tau_cld_mcica(ncld_grp,ntot_gpt,pverp) )
      allocate( tau_cld_co2(ntot_wavlnrng,pverp) )
      allocate( singscat_cld_co2(ntot_wavlnrng,pverp) )
      allocate( asym_cld_co2(ntot_wavlnrng,pverp) )
      allocate( tau_aer(ntot_wavlnrng,pverp) )
      allocate( wtau_aer(ntot_wavlnrng,pverp) )
      allocate( gwtau_aer(ntot_wavlnrng,pverp) )
      allocate( cFRC_mcica(ntot_gpt,pverp) )
      allocate( cICE_mcica(ntot_gpt,pverp) )
      allocate( cLIQ_mcica(ntot_gpt,pverp) )
    endif
!

    ! initialize internal RT variables
    tau_gas(:,:) = 0.0
    tau_cld_mcica(:,:,:) = 0.0
    singscat_cld_mcica(:,:,:) = 0.0
    asym_cld_mcica(:,:,:) = 0.0
    tau_cld_co2(:,:) = 0.0
    singscat_cld_co2(:,:) = 0.0
    asym_cld_co2(:,:) = 0.0
    tau_aer(:,:) = 0.0
    wtau_aer(:,:) = 0.0
    gwtau_aer(:,:) = 0.0
    sfc_emiss(:) = 0.0
    sfc_albedo_dir(:) = 0.0
    sfc_albedo_dif(:) = 0.0


    ! Fraction of the interplanetary solar flux at top of atmosphere:
    sflux_frac = dble(1./ext_msdist)    ! [1/AU^2]

    ! Set amount in pseudo layer, extends above model top to infinity
    ! Set P, T, and gases in psuedo layer equal to the top model layer
    tmid(1)  = ext_tmid(1)      ! temperatures [K]
    pmid(1)  = ext_pint(1)      ! pressure [Pa]
    qH2O(1)  = ext_H2O(1)       ! H2O vapor mass concentration (specific humdity) [kg/kg]
    qCO2(1)  = ext_CO2(1)       ! CO2 mass mixing ratio [kg/kg]
    qCH4(1)  = ext_CH4(1)       ! CH4 mass mixing ratio [kg/kg]
    qC2H6(1) = ext_C2H6(1)      ! C2H6 mass mixing ratio [kg/kg]
    qNH3(1)  = ext_NH3(1)       ! NH3 mass mixing ratio [kg/kg]
    qCO(1)   = ext_CO(1)        ! CO mass mixing ratio [kg/kg]
    qH2(1)   = ext_H2(1)        ! H2 mass mixing ratio [kg/kg]
    qN2(1)   = ext_N2(1)        ! N2 mass mixing ratio [kg/kg]
    qO3(1)   = ext_O3(1)        ! O3 mass mixing ratio [kg/kg]
    qO2(1)   = ext_O2(1)        ! O2 mass mixing ratio [kg/kg]


    ! Per-column dry-air properties: from the keyword tail when supplied
    ! (1-D and library callers always pass them), else the physconst module
    ! values (CAM path, where CAM owns mwdry/cpair)
    if (present(ext_mwdry)) then
      mwdry_col = ext_mwdry
    else
      mwdry_col = mwdry
    endif
    if (present(ext_cpdry)) then
      cpair_col = ext_cpdry
    else
      cpair_col = cpair
    endif
    if (present(ext_grav)) then
      grav_col = ext_grav
    else
      grav_col = SHR_CONST_G
    endif

    ! Optional condensed-phase inputs: default to zero everywhere (includes
    ! the pseudo layer, which never carries cloud or aerosol)
    have_cld_h2o = present(ext_cicewp) .and. present(ext_cliqwp) .and. &
                   present(ext_cfrc) .and. present(ext_rei) .and. present(ext_rel)
    have_cld_co2 = present(ext_cicewp_co2) .and. present(ext_rei_co2)
    cICE(:) = 0.0     ! in cloud ice water path [g/m2]
    cLIQ(:) = 0.0     ! in cloud liquid water path [g/m2]
    cFRC(:) = 0.0     ! cloud fraction
    REI(:)  = 0.0     ! ice cloud particle effective radii [microns]
    REL(:)  = 0.0     ! liquid cloud drop effective radii [microns]
    cICE_co2(:) = 0.0
    REI_co2(:)  = 0.0

    ! Set amounts in midlayers elsewhere
    do k=2, pverp
      qH2O(k)  = ext_H2O(k-1)
      qCO2(k)  = ext_CO2(k-1)
      qCH4(k)  = ext_CH4(k-1)
      qC2H6(k) = ext_C2H6(k-1)
      qNH3(k)  = ext_NH3(k-1)
      qCO(k)   = ext_CO(k-1)
      qH2(k)   = ext_H2(k-1)
      qN2(k)   = ext_N2(k-1)
      qO3(k)   = ext_O3(k-1)
      qO2(k)   = ext_O2(k-1)
      tmid(k)  = ext_tmid(k-1)
      pmid(k)  = ext_pmid(k-1)
    enddo

    ! H2O cloud condensate/radii/fraction, mapped to mid layers if supplied
    if (have_cld_h2o) then
      do k=2, pverp
        cICE(k)  = ext_cicewp(k-1)
        cLIQ(k)  = ext_cliqwp(k-1)
        cFRC(k)  = ext_cfrc(k-1)
        REI(k)   = ext_rei(k-1)
        REL(k)   = ext_rel(k-1)
      enddo
    endif

    ! CO2 ice cloud condensate/radii, mapped to mid layers if supplied
    if (have_cld_co2) then
      do k=2, pverp
        cICE_co2(k) = ext_cicewp_co2(k-1)
        REI_co2(k)  = ext_rei_co2(k-1)
      enddo
    endif

    ! CARMA haze mass mixing ratios, mapped to mid layers if supplied
    ! (pseudo layer contains no aerosol)
    qcarma(:,:,:) = 0.0
    if (present(ext_carmammr)) then
      do k=2, pverp
        qcarma(k,:,:) = ext_carmammr(k-1,:,:)
      enddo
    endif

    ! Set ground (surface) values:
    sfc_tempk = ext_sfcT   ! [K]
    sfc_press = ext_sfcP   ! [pa]

    ! Set interface temperatures
    tint(pverp) = sfc_tempk
    tint(1) = ext_tmid(1)

    do k = 2, pver
      dy = (log10(ext_pint(k)) - log10(ext_pmid(k))) / (log10(ext_pmid(k-1)) - log10(ext_pmid(k)))
      tint(k) = ext_tmid(k) - dy * (ext_tmid(k) - ext_tmid(k-1))
    enddo

    ! Define molecular column density at each layer [molec m-2]
    ! We use coldens_dry for computing column densities.

    ! Set column density in layer above top boundary
    coldens_dry(1) = (ext_pint(1)/grav_col)*(1.0-qh2o(1)) * SHR_CONST_AVOGAD/mwdry_col
    coldens(1) = (ext_pint(1)/grav_col) * SHR_CONST_AVOGAD/mwdry_col

    ! Set column density for other mid layers
    do k=2, pverp
      coldens(k) = (ext_pdel(k-1)/grav_col) * SHR_CONST_AVOGAD/mwdry_col
      !coldens_dry(k) = (ext_pdeldry(k-1)/grav_col) * SHR_CONST_AVOGAD/mwdry_col  !gives identical answer as below
      coldens_dry(k) = (ext_pdel(k-1)/grav_col)*(1.0-qh2o(k)) * SHR_CONST_AVOGAD/mwdry_col
    enddo

    ! Define mass column density in each layer [kg m-2]
    dzc(:) = ext_pdel(:)/grav_col

    ! Air mass path per radiative layer [kg m-2], same mapping as coldens:
    ! index 1 is the buffer layer from the model top to space.
    masspath(1) = ext_pint(1)/grav_col
    do k=2, pverp
      masspath(k) = ext_pdel(k-1)/grav_col
    enddo

    ! Define height of each layer [m]
    zlayer(1) = 0.0   !thickness of layer with lower boundary at model top is zero
! this only affects the CIA  calculations
!    zlayer(1) = 1.0   !thickness of layer with lower boundary at model top is zero
    do k=2, pverp
      zlayer(k) = (ext_zint(k-1) - ext_zint(k))
    enddo

    !=======================================================================================
    ! Surface Albedo and Emissivity Treatment
    !
    ! NOTES: The surface albedo is set from the cam_in variable.  The current implementation
    ! uses a three channel gray albedo/emissiviy scheme, divided between "visible", "near-IR",
    ! and "thermal" bands.
    !
    ! For the Earth-Sun combination, the shortwave and longwave streams are generally
    ! non-overlapping and can be divided cleanly at ~5 microns (2000 cm-1).  However, in
    ! the limits of very hot planets, or those around very red stars, shortwave and longwave
    ! radiative streams begin to overlap.  Then, the assumption of a broadband albedo and
    ! emissivity split at 5 microns begins to break down.
    !
    ! Here, the radiative processes in the atmosphere are more appropriately represented by
    ! by assuming an emissivitiy of 1 across the whole spectra for the longwave stream, while
    ! extending the near-IR albedo through the "thermal" wavelengths for the shortwave stream.
    ! While stricitly speaking we have decoupled the emissivity from the albedo, however, note
    ! that Kirchoff's law applys at a specific wavelenth and not neccessarily across broadband
    ! regions as are used here.
    !
    ! Presently the visible and near-IR direct and diffuse albedos are calculated in the surface
    ! models and then passed through the coupler to the atmosphere model and radiation.  The surface
    ! emissivity currently is not coupled between surface and atmosphere models through the coupler,
    ! and instead is set seperately to ~1 in atmosphere, land, ocean, and ice componenets respectively.
    !

    ! Broadband surface thermal emissivity: 1.0 unless an input value is supplied.
    srf_emiss_val = 1.0
    if (present(ext_srf_emiss)) srf_emiss_val = ext_srf_emiss

    ! Set spectral interval albedos and emissivities from broadband quantities
    do iw=1,ntot_wavlnrng    ! Loop over relevant wavelength intervals
      if (wavenum_edge(iw) .le. 2000) then  ! "thermal"
        sfc_albedo_dir(iw) = ext_aldir
        sfc_albedo_dif(iw) = ext_aldif
        sfc_emiss(iw) = srf_emiss_val
      endif
      if (wavenum_edge(iw) .gt. 2000 .and. wavenum_edge(iw) .lt. 13000) then   ! "near-IR"
        sfc_albedo_dir(iw) = ext_aldir
        sfc_albedo_dif(iw) = ext_aldif
        sfc_emiss(iw) = srf_emiss_val
      endif
      if (wavenum_edge(iw) .ge. 13000) then     ! "visible"
        sfc_albedo_dir(iw) = ext_asdir
        sfc_albedo_dif(iw) = ext_asdif
        sfc_emiss(iw) = srf_emiss_val
      endif
    enddo

    !if (masterproc) then
    !    write(*,*) "sfc_albedo", sfc_albedo_dir, sfc_albedo_dif, "sfc_emiss",sfc_emiss
    !endif

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

      ! On exceedingly rare occurences, zenith angles at the far horizon (89.9<z<90 deg)
      ! will be injested from the orbital calculations and can cause numerical faults.
      cos_mu = max(cos_mu, 1.0e-3)        ! Trim far horizon zenith angles to ~89.9 deg

    else     ! (cos_mu < 1.d-6)
      sw_on = .FALSE.       ! Sun below horizon, do only longwave
    endif


    ! Spectral range the solver uses (union of the optimized SW and LW ranges).
    ! The optics below are computed only for these bands; the rest are never read.
    iwlo = min(sw_iwbeg, lw_iwbeg)
    iwhi = max(sw_iwend, lw_iwend)

    call calc_opd_gas(tmid, pmid/100.0, coldens, coldens_dry, &
                     qH2O, qCO2, qCH4, qC2H6, qNH3, qCO, qO2, qO3, qH2, qN2, &
                     zlayer*100.0, mwdry_col, tau_gas, tau_ray, iwlo=iwlo, iwhi=iwhi)

    ! Cloud optics (tau/ssa/asym already zeroed above; when do_exo_clouds is
    ! .false., or the optional condensate arguments are absent, the cloud
    ! path is skipped entirely and contributes no opacity).
    if (do_exo_clouds) then
      if (have_cld_h2o) then
        ! In 1-D get_nstep() is the constant 9404 (time_manager stub); the
        ! optional ext_icol offset decorrelates the MCICA subcolumns across
        ! batch columns (absent or 1 = the legacy constant seed).
        mcica_seed = get_nstep()
        if (present(ext_icol)) mcica_seed = mcica_seed + ext_icol - 1
        call calc_opd_cld_h2o(ext_pint, cICE, cLIQ, REI, REL, cFRC, tau_cld_mcica, singscat_cld_mcica, &
                         asym_cld_mcica, cFRC_mcica, cICE_mcica, cLIQ_mcica, mcica_seed, &
                         iwlo=iwlo, iwhi=iwhi )
      endif
      if (have_cld_co2) then
        call calc_opd_cld_co2(cICE_co2, REI_co2, tau_cld_co2, singscat_cld_co2, asym_cld_co2, &
                              iwlo=iwlo, iwhi=iwhi)
      endif
    endif

    ! CARMA haze aerosol optics (same gating principle as clouds: skipped
    ! entirely, zero opacity, when do_exo_haze is .false. or carmammr is
    ! not supplied).
    if (do_exo_haze .and. present(ext_carmammr)) then
      call calc_opd_aero(qcarma, masspath, tau_aer, wtau_aer, gwtau_aer, iwlo=iwlo, iwhi=iwhi)
    endif

    ! fused solver: whole two-stream chain per spectral block
    gstart(1) = 1
    do iw=1,ntot_wavlnrng
      gstart(iw+1) = gstart(iw) + ngauss_pts(iw)
    enddo
    lw_dnflux(:) = 0. ; lw_upflux(:) = 0. ; sw_upflux(:) = 0. ; sw_dnflux(:) = 0.
    lw_dnflux_spec(:,:) = 0. ; lw_upflux_spec(:,:) = 0.
    sw_upflux_spec(:,:) = 0. ; sw_dnflux_spec(:,:) = 0.
    vis_dir = 0. ; vis_dif = 0. ; nir_dir = 0. ; nir_dif = 0. ; sol_toa = 0.
    nbw = max(1, nblk_gpt/ngauss_8gpt)
    do iwb = iwlo, iwhi, nbw
      bw0 = iwb
      bw1 = min(iwb+nbw-1, iwhi)
      gb0 = gstart(bw0)
      gb1 = gstart(bw1+1)-1
      call solve_block(pmid/100.0, tmid, tint, swcut, tau_gas, tau_ray, &
                       tau_cld_mcica, singscat_cld_mcica, asym_cld_mcica, &
                       tau_cld_co2, singscat_cld_co2, asym_cld_co2, &
                       tau_aer, wtau_aer, gwtau_aer, &
                       part_in_tshadow, sfc_albedo_dir, sfc_albedo_dif, sfc_emiss, &
                       sflux_frac, sfc_tempk, cos_mu, sw_on, &
                       lw_dnflux, lw_upflux, sw_upflux, sw_dnflux, &
                       lw_dnflux_spec, lw_upflux_spec, sw_upflux_spec, sw_dnflux_spec, &
                       vis_dir, vis_dif, nir_dir, nir_dif, sol_toa)
    enddo
    call rad_postcalc_final(dzc, cpair_col, swcut, part_in_tshadow, sw_on, &
                            sw_dTdt, lw_dTdt, lw_dnflux, lw_upflux, sw_upflux, sw_dnflux)

    return

  end subroutine aerad_driver


!============================================================================

  subroutine rad_precalc(pmid, tmid, tint, swcut, tau_gas, tau_ray, &
                         tau_cld_mcica, singscat_cld_mcica, asym_cld_mcica, &
                         tau_cld_co2, singscat_cld_co2, asym_cld_co2, &
                         tau_aer, wtau_aer, gwtau_aer, &
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
    integer,  intent(in) :: swcut
    real(r8), intent(in), dimension(ntot_gpt,pverp) ::  tau_gas
    real(r8), intent(in), dimension(ntot_wavlnrng,pverp) ::  tau_ray
    real(r8), intent(in), dimension(ncld_grp,ntot_gpt,pverp) ::  tau_cld_mcica
    real(r8), intent(in), dimension(ncld_grp,ntot_gpt,pverp) ::  singscat_cld_mcica
    real(r8), intent(in), dimension(ncld_grp,ntot_gpt,pverp) ::  asym_cld_mcica
    ! CO2 ice cloud optics, band-indexed (cloud fraction = 1, no MCICA overlap).
    real(r8), intent(in), dimension(ntot_wavlnrng,pverp) ::  tau_cld_co2
    real(r8), intent(in), dimension(ntot_wavlnrng,pverp) ::  singscat_cld_co2
    real(r8), intent(in), dimension(ntot_wavlnrng,pverp) ::  asym_cld_co2
    ! CARMA haze aerosol optics, band-indexed, presummed over elements/bins
    ! (tau, w*tau, g*w*tau) so they drop straight into the accumulation below.
    real(r8), intent(in), dimension(ntot_wavlnrng,pverp) ::  tau_aer
    real(r8), intent(in), dimension(ntot_wavlnrng,pverp) ::  wtau_aer
    real(r8), intent(in), dimension(ntot_wavlnrng,pverp) ::  gwtau_aer
    logical,  intent(in) :: part_in_tshadow
    real(r8), intent(in), dimension(ntot_wavlnrng) :: sfc_albedo_dir
    real(r8), intent(in), dimension(ntot_wavlnrng) :: sfc_albedo_dif
    real(r8), intent(in), dimension(ntot_wavlnrng) :: sfc_emiss
    real(r8), intent(in) :: sflux_frac
    real(r8), intent(in) :: sfc_tempk
    real(r8), intent(in) :: cos_mu
    logical,  intent(in) :: sw_on

    real(r8), intent(out), dimension(gb0:gb1,ngangles,pverp) :: Y3
    real(r8), intent(out), dimension(gb0:gb1,pverp) :: TAUL
    real(r8), intent(out), dimension(gb0:gb1,pverp) :: OPD
    real(r8), intent(out), dimension(gb0:gb1,pverp) :: PTEMP
    real(r8), intent(out), dimension(gb0:gb1) :: PTEMPG
    real(r8), intent(out), dimension(gb0:gb1,pverp) :: SLOPE
    real(r8), intent(out), dimension(gb0:gb1,pverp) :: SOL
    real(r8), intent(out), dimension(gb0:gb1,pverp) :: W0
    real(r8), intent(out), dimension(gb0:gb1,pverp) :: G0
    real(r8), intent(out), dimension(gb0:gb1) :: EMIS
    real(r8), intent(out), dimension(gb0:gb1) :: RSFXdir
    real(r8), intent(out), dimension(gb0:gb1) :: RSFXdif

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
    !
    ! Map/average important parameters to gauss points (it):
    !

    ! Map emissivity and surface relfecitivities to gauss point grid
    it = gstart(bw0)-1
    do iw=bw0,bw1    ! block bands
      do ig=1,ngauss_pts(iw)
        it = it+1
        EMIS(it) = sfc_emiss(iw)
        RSFXdir(it) = sfc_albedo_dir(iw)
        RSFXdif(it) = sfc_albedo_dif(iw)
      enddo
    enddo

    if(sw_on) then
      if(part_in_tshadow) then   ! Eliminate direct "shortwave" flux within shadow

        it = gstart(max(sw_iwbeg,bw0))-1
        do iw=max(sw_iwbeg,bw0),min(sw_iwend,bw1)
          do ig=1,ngauss_pts(iw)
            it = it+1
            SOL(it,1:swcut-1) = gw_solflux(it)*sflux_frac
            SOL(it,swcut:pverp) = 0.d0
          enddo
        enddo

      else    ! No topographic shadowing in current column

        it = gstart(max(sw_iwbeg,bw0))-1
        do iw=max(sw_iwbeg,bw0),min(sw_iwend,bw1)
          do ig=1,ngauss_pts(iw)
            it = it+1
            SOL(it,1:pverp) = gw_solflux(it)*sflux_frac
          enddo
        enddo

      endif
    endif
    !if (sw_on) write(6,*) 'sw_on is true in precalc, solin1 is',solin1,'solinl1 is',solinl1

    OPD(:,:) = 0.d0     ! Initialize necessary array portions in
    TAUL(:,:) = 0.d0    !  anticipation of summing below

    do k=1,pverp     ! Loop over all layer BOUNDARIES

      k_1 = max(1,k-1)   ! index for level above (except for k=1, of course)
      it = gstart(bw0)-1

      do iw=bw0,bw1      ! block bands
        do ig=1,ngauss_pts(iw)
          it = it+1

          ! Combine the aerosol and gas optical parameters to get the total
          !  "effective" layer parameters:

          taul_ig = tau_gas(it,k) + tau_ray(iw,k)
          w0_ig = tau_ray(iw,k)
          g0_ig = 0.

          ! Add cloud optical depths
          do ip =1,2
            taul_ig = taul_ig + tau_cld_mcica(ip,it,k)
            w0_ig = w0_ig + singscat_cld_mcica(ip,it,k) * tau_cld_mcica(ip,it,k)
            g0_ig = g0_ig + asym_cld_mcica(ip,it,k) * singscat_cld_mcica(ip,it,k) * tau_cld_mcica(ip,it,k)
          enddo

          ! Add CO2 ice cloud optical depths (band-indexed; cloud fraction = 1,
          ! so the same opacity applies to every gauss point in band iw, no MCICA
          ! subcolumn sampling). Same accumulation algebra as the H2O groups.
          taul_ig = taul_ig + tau_cld_co2(iw,k)
          w0_ig = w0_ig + singscat_cld_co2(iw,k) * tau_cld_co2(iw,k)
          g0_ig = g0_ig + asym_cld_co2(iw,k) * singscat_cld_co2(iw,k) * tau_cld_co2(iw,k)

          ! Add CARMA haze aerosol optical depths (band-indexed; the element/bin
          ! sums of tau, w*tau, and g*w*tau were formed in calc_opd_aero, so the
          ! accumulation algebra matches the cloud groups above).
          taul_ig = taul_ig + tau_aer(iw,k)
          w0_ig = w0_ig + wtau_aer(iw,k)
          g0_ig = g0_ig + gwtau_aer(iw,k)

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
        do it=max(lw_ipbeg,gb0),min(lw_ipend,gb1)   ! Y3 feeds LW only
          X = TAUL(it,k)/g_angle(ia)
          if(X <= 400.d0) then
            Y3(it,ia,k) = min(1.d0,exp(-X))
          else
            Y3(it,ia,k) = 0.d0    ! would be effectively zero anyway
          endif
        enddo
      enddo
    enddo

    ! Calculate Planck function for current temperatures profiles
    call dplanck(pmid, tmid, tint, TAUL, sfc_tempk, PTEMP, PTEMPG, SLOPE)

    return

  end subroutine rad_precalc

!============================================================================

  subroutine dplanck(pmid, tmid, tint, TAUL, sfc_tempk, PTEMP, PTEMPG, SLOPE)

!------------------------------------------------------------------------
!
! Purpose: Calculate Planck function and its derivative at the ground and
!          at all altitudes.
!
!------------------------------------------------------------------------

    implicit none

!------------------------------------------------------------------------
!
! Arguments
!
    real(r8), intent(in), dimension(pverp) :: pmid       ! [mb] pressures at mid layers
    real(r8), intent(in), dimension(pverp) :: tmid       ! [K] temperatures at mid layers (isothermal top)
    real(r8), intent(in), dimension(pverp) :: tint       ! [K] temperatures at level interfaces
    real(r8), intent(in), dimension(gb0:gb1,pverp) :: TAUL   ! Layer optical depth
    real(r8), intent(in) :: sfc_tempk
    real(r8), intent(out), dimension(gb0:gb1,pverp) :: PTEMP   ! Planck function at each layer (isothermal top)
    real(r8), intent(out), dimension(gb0:gb1) :: PTEMPG        ! Planck function evaluated at ground
    real(r8), intent(out), dimension(gb0:gb1,pverp) :: SLOPE

!------------------------------------------------------------------------
!
! Local Variables
!
    real(r8), dimension(gb0:gb1,pverp) :: ptemp_midlayer   ! block-sized
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

    do ip=max(lw_ipbeg,gb0),min(lw_ipend,gb1)   ! block LW

      it = int(sfc_tempk*tpft_finc)
      ft = dble(sfc_tempk*tpft_finc-real(it))

      ! If temperature out of planck table range
      if ( (sfc_tempk*tpft_finc > tpft_end) .and. (masterproc) ) then
        write(*,*) "WARNING: surface temperature exceeds planck table maximum", tint(k),k
      endif
      if ( (sfc_tempk*tpft_finc < tpft_beg) .and. (masterproc) ) then
        write(*,*) "WARNING: surface temperature below planck table minimum", tint(k),k
      endif

      ! Interpolate between table values:
      PTEMPG(ip) = ptemp_itp(it,ip)*(1.d0-ft)+ptemp_itp(it+1,ip)*ft

    enddo

    !---- Diagnostic output -----
    !---- comment out for simulations
!    if (masterproc) then
!      write(*,*) "surface planck function [W m-2] in each spectral interval"
!      ip=0
!      do iw=1,ntot_wavlnrng
!        do ig=1,ngauss_pts(iw)
!          ip=ip+1
!          write(*,*) iw,PTEMPG(ip)/g_weight(ip), PTEMPG(ip)
!        enddo
!      enddo
!      write(*,*) "TOTAL PLANCK FUNCTION:", SUM(PTEMPG)*SHR_CONST_PI
!      write(*,*) "SURFACE PLANCK FUNCTION:",SHR_CONST_STEBOL*sfc_tempk**4
!      write(*,*) "SURFACE TEMPERATURE:",sfc_tempk
!    endif
    !-----------

    ! Set Planck function at interfaces
    do k=1,pverp      ! Loop over all layers
      do ip=max(lw_ipbeg,gb0),min(lw_ipend,gb1)   ! block LW

        it = int(tint(k)*tpft_finc)
        ft = dble(tint(k)*tpft_finc-real(it))

        ! If temperature out of planck table range
        if ( (tint(k)*tpft_finc > tpft_end) .and. (masterproc) ) then
           write(*,*) "WARNING: temperature exceeds planck table maximum", tint(k),k
        endif
        if ( (tint(k)*tpft_finc < tpft_beg) .and. (masterproc) )then
           write(*,*) "WARNING: temperature below planck table minimum", tint(k),k
        endif

        ! Interpolate between table values:
        PTEMP(ip,k) = ptemp_itp(it,ip)*(1.d0-ft)+ptemp_itp(it+1,ip)*ft

        ! Determine midlayer planck function
        ! Slope derived from midlayer planck
        it = int(tmid(k)*tpft_finc)
        ft = dble(tmid(k)*tpft_finc-real(it))
        ptemp_midlayer(ip,k) = ptemp_itp(it,ip)*(1.d0-ft)+ptemp_itp(it+1,ip)*ft

        k_1 = max(1,k-1)   ! index for level above (except for k=1, of course)
        if(TAUL(ip,k) > SMALLd) then
          SLOPE(ip,k) = 2.0 * (ptemp_midlayer(ip,k)-ptemp_midlayer(ip,k_1)) / (TAUL(ip,k)+TAUL(ip,k_1))  ! DIML
        else
          SLOPE(ip,k) = 0.d0
        endif

      enddo

    enddo

    return

  end subroutine dplanck


!============================================================================

  subroutine two_stream(TAUL, W0, G0, EMIS, RSFXdir, RSFXdif, beamSolar, ip_ibeg, ip_iend, &
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

     real(r8), intent(in), dimension(gb0:gb1,pverp) :: TAUL
     real(r8), intent(in), dimension(gb0:gb1,pverp) :: W0
     real(r8), intent(in), dimension(gb0:gb1,pverp) :: G0
     real(r8), intent(in), dimension(gb0:gb1) :: EMIS
     real(r8), intent(in), dimension(gb0:gb1) :: RSFXdir
     real(r8), intent(in), dimension(gb0:gb1) :: RSFXdif
     logical, intent(in) :: beamSolar
     integer, intent(in) :: ip_ibeg
     integer, intent(in) :: ip_iend
     real(r8), intent(out), dimension(gb0:gb1,pverp) :: EM1
     real(r8), intent(out), dimension(gb0:gb1,pverp) :: EM2
     real(r8), intent(out), dimension(gb0:gb1,pverp) :: EL1
     real(r8), intent(out), dimension(gb0:gb1,pverp) :: EL2
     real(r8), intent(out), dimension(gb0:gb1,2*pverp) :: AF
     real(r8), intent(out), dimension(gb0:gb1,2*pverp) :: BF
     real(r8), intent(out), dimension(gb0:gb1,2*pverp) :: EF
     real(r8), intent(out), dimension(gb0:gb1,pverp) :: AK
     real(r8), intent(out), dimension(gb0:gb1,pverp) :: GAMI
     real(r8), intent(out), dimension(gb0:gb1,pverp) :: B1
     real(r8), intent(out), dimension(gb0:gb1,pverp) :: B2
     real(r8), intent(out), dimension(gb0:gb1,pverp) :: EE1

!------------------------------------------------------------------------
!
! Local Variables
!

    real(r8) :: x1
    integer :: ip
    integer :: k
    integer :: kd
    real(r8) u1i, u1i_2
    real(r8), dimension(gb0:gb1) :: srf_reflect_dif
    real(r8), dimension(gb0:gb1) :: srf_reflect_dir  ! not used

!------------------------------------------------------------------------
!
! Start Code
!

    ! HERE WE DEFINE LAYER PROPERTIES FOLLOWING GENERAL SCHEME
    !  OF MEADOR AND WEAVOR. THEN WE SET UP LAYER PROPERTIES
    !  NEEDED FOR MATRIX:

    if (beamSolar) then ! shortwave stream
      u1i_2 = U1I2sol
      srf_reflect_dif(:) = RSFXdif(:)     ! set shortwave surface diffuse reflectivity
    else ! thermal stream
      u1i_2 = U1I2ir
      srf_reflect_dif(:) = 1.0d0-EMIS(:)  ! set longwave surface diffuse reflectivity
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
      AF(ip,2*pverp) = EL1(ip,pverp)-srf_reflect_dif(ip)*EL2(ip,pverp)
      BF(ip,2*pverp) = EM1(ip,pverp)-srf_reflect_dif(ip)*EM2(ip,pverp)
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
  real(r8), intent(in), dimension(gb0:gb1,pverp) :: EM1
  real(r8), intent(in), dimension(gb0:gb1,pverp) :: EM2
  real(r8), intent(in), dimension(gb0:gb1,pverp) :: EL1
  real(r8), intent(in), dimension(gb0:gb1,pverp) :: EL2
  real(r8), intent(in), dimension(gb0:gb1,2*pverp) :: AF
  real(r8), intent(in), dimension(gb0:gb1,2*pverp) :: BF
  real(r8), intent(in), dimension(gb0:gb1,2*pverp) :: EF
  real(r8), intent(in), dimension(gb0:gb1,pverp) :: B1
  real(r8), intent(in), dimension(gb0:gb1,pverp) :: B2
  real(r8), intent(in), dimension(gb0:gb1,pverp) :: AK
  real(r8), intent(in), dimension(gb0:gb1,pverp) :: TAUL
  real(r8), intent(in), dimension(gb0:gb1,pverp) :: OPD
  real(r8), intent(in), dimension(gb0:gb1,pverp) :: PTEMP
  real(r8), intent(in), dimension(gb0:gb1) :: PTEMPG
  real(r8), intent(in), dimension(gb0:gb1,pverp) :: SLOPE
  real(r8), intent(in), dimension(gb0:gb1,pverp) :: SOL
  real(r8), intent(in), dimension(gb0:gb1,pverp) :: W0
  real(r8), intent(in), dimension(gb0:gb1,pverp) :: G0
  real(r8), intent(in), dimension(gb0:gb1) :: EMIS
  real(r8), intent(in), dimension(gb0:gb1) :: RSFXdir
  real(r8), intent(in), dimension(gb0:gb1) :: RSFXdif
  real(r8), intent(in) :: cos_mu
  logical, intent(in) :: sw_on
  logical, intent(in) :: beamSolar
  integer, intent(in) :: ip_ibeg
  integer, intent(in) :: ip_iend
  real(r8), intent(out), dimension(gb0:gb1,pverp) :: CK1
  real(r8), intent(out), dimension(gb0:gb1,pverp) :: CK2
  real(r8), intent(out), dimension(gb0:gb1,pverp) :: CPB
  real(r8), intent(out), dimension(gb0:gb1,pverp) :: CMB
  real(r8), intent(out), dimension(gb0:gb1,pverp) :: B3
  real(r8), intent(out), dimension(gb0:gb1,pverp) :: DIRECT

!------------------------------------------------------------------------
!
! Local Variables
!

    real(r8), dimension(gb0:gb1,pverp) :: CPP
    real(r8), dimension(gb0:gb1,pverp) :: CM
    real(r8), dimension(gb0:gb1,pverp) :: EE3
    real(r8), dimension(gb0:gb1,pverp) :: EL3
    real(r8), dimension(gb0:gb1,2*pverp) :: AS
    real(r8), dimension(gb0:gb1,2*pverp) :: DF
    real(r8), dimension(gb0:gb1,2*pverp) :: DS
    real(r8), dimension(gb0:gb1,2*pverp) :: XK

    real(r8), dimension(gb0:gb1) :: sfcs
    real(r8), dimension(gb0:gb1) :: srf_reflect_dif
    real(r8), dimension(gb0:gb1) :: srf_reflect_dir  ! not used
    real(r8) :: DUo
    real(r8) :: B4
    real(r8) :: X2
    real(r8) :: X3
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

        do ip=ip_ibeg,ip_iend

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

      do ip=ip_ibeg,ip_iend
        sfcs(ip) = DIRECT(ip,pverp)*RSFXdir(ip)
        srf_reflect_dif(ip) = RSFXdif(ip)
      enddo

    else  ! beamSolar

    ! ******************************
    ! * CALCULATIONS FOR INFRARED. *
    ! ******************************
      do k=1,pverp
        kindex = max(1,k-1)
        do ip=ip_ibeg,ip_iend
          B3(ip,k) = 1.d0/(B1(ip,k)+B2(ip,k))
          CPP(ip,k) = (PTEMP(ip,KINDEX)+SLOPE(ip,k)*B3(ip,k))*U1Sir  ! TIM mod
          CPB(ip,k) = CPP(ip,k)+SLOPE(ip,k)*TAUL(ip,k)*U1Sir
          CM(ip,k) = (PTEMP(ip,KINDEX)-SLOPE(ip,k)*B3(ip,k))*U1Sir
          CMB(ip,k) = CM(ip,k)+SLOPE(ip,k)*TAUL(ip,k)*U1Sir
          EL3(ip,k)    = 0.d0
          DIRECT(ip,k) = 0.d0
          EE3(ip,k)    = 0.d0
        enddo
      enddo

      do ip=ip_ibeg,ip_iend
        sfcs(ip) = EMIS(ip)*PTEMPG(ip)*SHR_CONST_PI
        srf_reflect_dif(ip) = 1.0d0-EMIS(ip)
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
      DF(ip,2*pverp) = SFCS(ip)+srf_reflect_dif(ip)*CMB(ip,pverp)-CPB(ip,pverp)
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
        if(abs(x) < x_lim) then     !�
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

  subroutine refine_lwflux (CK1, CK2, Y3, AK, GAMI, B3, EE1, TAUL, PTEMP, PTEMPG, SLOPE, EMIS, RSFXdir,RSFXdif, cos_mu, DIRECTU, DIREC, ip_ibeg, ip_iend)

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
    integer, intent(in) :: ip_ibeg, ip_iend   ! LW g-points of this block
    real(r8), intent(in), dimension(gb0:gb1,pverp) :: CK1
    real(r8), intent(in), dimension(gb0:gb1,pverp) :: CK2
    real(r8), intent(in), dimension(gb0:gb1,ngangles,pverp) :: Y3
    real(r8), intent(in), dimension(gb0:gb1,pverp) :: AK
    real(r8), intent(in), dimension(gb0:gb1,pverp) :: GAMI
    real(r8), intent(in), dimension(gb0:gb1,pverp) :: B3
    real(r8), intent(in), dimension(gb0:gb1,pverp) :: EE1
    real(r8), intent(in), dimension(gb0:gb1,pverp) :: TAUL
    real(r8), intent(in), dimension(gb0:gb1,pverp) :: PTEMP
    real(r8), intent(in), dimension(gb0:gb1) :: PTEMPG
    real(r8), intent(in), dimension(gb0:gb1,pverp) :: SLOPE
    real(r8), intent(in), dimension(gb0:gb1) :: EMIS
    real(r8), intent(in), dimension(gb0:gb1) :: RSFXdir
    real(r8), intent(in), dimension(gb0:gb1) :: RSFXdif
    real(r8), intent(in) :: cos_mu

    real(r8), intent(out), dimension(gb0:gb1,pverp) :: DIRECTU
    real(r8), intent(out), dimension(gb0:gb1,pverp) :: DIREC

!
!     *  Input               :  PTEMP, SLOPE, Y3, B3, EE1, EE2     *
!     *  Output              :  Irradiances: DIREC, DIRECTU
!
!------------------------------------------------------------------------
!
! Local Variables
!
    real(r8), dimension(gb0:gb1,ngangles,pverp) :: uintent
    real(r8), dimension(gb0:gb1,ngangles,pverp) :: dintent

    integer :: openstatus
    real(r8), dimension(gb0:gb1,ngangles,pverp) :: Y1
    real(r8), dimension(gb0:gb1,ngangles,pverp) :: Y2
    real(r8), dimension(gb0:gb1,ngangles,pverp) :: Y4
    real(r8), dimension(gb0:gb1,ngangles,pverp) :: Y8
    real(r8), dimension(gb0:gb1,pverp) :: Y5
    real(r8), dimension(gb0:gb1,pverp) :: A1
    real(r8), dimension(gb0:gb1,pverp) :: A2
    real(r8), dimension(gb0:gb1,pverp) :: A3
    real(r8), dimension(gb0:gb1,pverp) :: A4
    real(r8), dimension(gb0:gb1,pverp) :: A7
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
      do ip=ip_ibeg,ip_iend
        ! HERE WE DO NO SCATTERING COEFFICIENTS
        A3(ip,k) = PTEMP(ip,KINDEX)*2.0*SHR_CONST_PI
        A4(ip,k) = 2.0*SHR_CONST_PI*SLOPE(ip,k)
        A7(ip,k) = A3(ip,k)
        Y5(ip,k) = A4(ip,k)*TAUL(ip,k)
      enddo
      ! HERE WE DO SCATTERING
      do ip=ip_ibeg,ip_iend
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

        do ip=ip_ibeg,ip_iend
          !%waste CPU%         Y1(ip,ia,k)  =  0.d0
          !%  Y2(ip,ia,k)  =  0.d0
          Y4(ip,ia,k) = A7(ip,k)-A4(ip,k)*g_angle(ia)
          Y8(ip,ia,k) = A3(ip,k)+A4(ip,k)*g_angle(ia)
          !%      enddo

          ! HERE WE DO SCATTERING

          !% do ip=ip_ibeg,ip_iend
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
      do ip=ip_ibeg,ip_iend
        DIREC(ip,k) = 0.d0
        DIRECTU(ip,k) = 0.d0
      enddo
    enddo

    ! DIREC IS DOWNWARD IRRADIANCE. DIRECTU IS UPWARD IRRADIANCE.
    !  CALCULATE DINTENT THE DOWNWARD RADIANCE AND DIREC THE DOWNWARD IRRADIANCE

    ! BOUNDARY CONDITIONS: DOWNWARD IRRADIANCE, RADIANCE AT TOA (k = camtop) was 1

    cos_mu1 = max(cos_mu,0.0)

    do ia=1,ngangles
      do ip=ip_ibeg,ip_iend

        DINTENT(ip,ia,camtop) = (1.d0-Y3(ip,ia,camtop))*Y4(ip,ia,camtop)+Y1(ip,ia,camtop)

        DIREC(ip,camtop) = DIREC(ip,camtop)+DINTENT(ip,ia,camtop)*g_ang_weight(ia)

      enddo
    enddo


    ! DINTENT IS DOWNWARD RADIANCE * TwoPI. DIREC IS THE DOWNWARD IRRADIANCE.
    !  CALCULATE FOR REST OF ATMOSPHERE.
    do k=camtop+1, pverp     !was k=2
      do ia=1,ngangles
        do ip=ip_ibeg,ip_iend
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
      do ip=ip_ibeg,ip_iend

        ! ETW:  assume no reflection of downwelling thermal radiation, until a more
        ! comprehensive albedo/emissivity treatment is implemented
        UINTENT(ip,ia,pverp) = EMIS(ip)*PTEMPG(ip)*2.0*SHR_CONST_PI+(1.0-EMIS(ip))*DIREC(ip,pverp)*2.0d0

        DIRECTU(ip,pverp) = DIRECTU(ip,pverp)+UINTENT(ip,ia,pverp)*g_ang_weight(ia)

      enddo
    enddo

    ! CALCULATE FOR THE REST OF THE ATMOSPHERE
    ! ITERATED FROM BOTTOM UP
    kd = pverp
    do k=camtop+1,pverp   ! was k=2
      kd = kd-1
      do ia=1,ngangles
        do ip=ip_ibeg,ip_iend

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


  subroutine rad_postcalc_accum (CK1s, CK2s, CPBs, CMBs, EM1s, EM2s, EL1s, EL2s, &
                                 DIRECTs, DIRECTU, DIREC, SOL, cos_mu, sw_on, &
                                 lw_dnflux, lw_upflux, sw_upflux, sw_dnflux, &
                                 lw_dnflux_spec, lw_upflux_spec, sw_upflux_spec, sw_dnflux_spec, &
                                 vis_dir, vis_dif, nir_dir, nir_dif, sol_toa)
!------------------------------------------------------------------------
! the flux sums of the original rad_postcalc, for the bands of
! the current block only. Blocks are visited in ascending band order, so each
! running sum is formed in exactly the original order.
!------------------------------------------------------------------------
    implicit none
    real(r8), intent(in), dimension(gb0:gb1,pverp) :: CK1s, CK2s, CPBs, CMBs
    real(r8), intent(in), dimension(gb0:gb1,pverp) :: EM1s, EM2s, EL1s, EL2s
    real(r8), intent(in), dimension(gb0:gb1,pverp) :: DIRECTs, DIRECTU, DIREC, SOL
    real(r8), intent(in) :: cos_mu
    logical,  intent(in) :: sw_on
    real(r8), intent(inout), dimension(pverp) :: lw_dnflux, lw_upflux, sw_upflux, sw_dnflux
    real(r8), intent(inout), dimension(pverp,ntot_wavlnrng) :: lw_dnflux_spec, lw_upflux_spec
    real(r8), intent(inout), dimension(pverp,ntot_wavlnrng) :: sw_upflux_spec, sw_dnflux_spec
    real(r8), intent(inout) :: vis_dir, vis_dif, nir_dir, nir_dif, sol_toa
    integer :: k, ip, ipc, iw, ig, w0, w1

    w0 = max(lw_iwbeg,bw0) ; w1 = min(lw_iwend,bw1)
    if (w0 <= w1) then
      do k=camtop,pverp
        ip = gstart(w0)
        do iw=w0,w1
          do ig=1,ngauss_pts(iw)
            lw_dnflux(k) = lw_dnflux(k)+DIREC(ip,k)
            lw_upflux(k) = lw_upflux(k)+DIRECTU(ip,k)
            lw_dnflux_spec(k,iw) = lw_dnflux_spec(k,iw)+DIREC(ip,k)
            lw_upflux_spec(k,iw) = lw_upflux_spec(k,iw)+DIRECTU(ip,k)
            ip = ip+1
          enddo
        enddo
      enddo
    endif

    ! Night side: the original summed exact zeros here (add_txrad left the
    ! SW coefficients zero), so skipping leaves every SW output at zero.
    if (.not. sw_on) return
    w0 = max(sw_iwbeg,bw0) ; w1 = min(sw_iwend,bw1)
    if (w0 > w1) return

    do k=camtop,pverp
      ip = gstart(w0)
      do iw=w0,w1
        do ig=1,ngauss_pts(iw)
          sw_upflux(k) = sw_upflux(k)+CK1s(ip,k)*EL1s(ip,k)+  &
                         CK2s(ip,k)*EM1s(ip,k)+CPBs(ip,k)
          sw_dnflux(k) = sw_dnflux(k)+CK1s(ip,k)*EL2s(ip,k)+  &
                         CK2s(ip,k)*EM2s(ip,k)+CMBs(ip,k)+DIRECTs(ip,k)
          sw_upflux_spec(k,iw) = sw_upflux_spec(k,iw)+CK1s(ip,k)*EL1s(ip,k)+  &
                                 CK2s(ip,k)*EM1s(ip,k)+CPBs(ip,k)
          sw_dnflux_spec(k,iw) = sw_dnflux_spec(k,iw)+CK1s(ip,k)*EL2s(ip,k)+  &
                                 CK2s(ip,k)*EM2s(ip,k)+CMBs(ip,k)+DIRECTs(ip,k)
          ip = ip+1
        enddo
      enddo
    enddo

    ip = gstart(w0)
    do iw=w0,w1
      do ig=1,ngauss_pts(iw)
        sol_toa = sol_toa + SOL(ip,1)*cos_mu
        ip = ip+1
      enddo
    enddo

    ipc = gstart(w0)
    do iw=w0,w1
      if (wavenum_edge(iw) .gt. 13000) then
        do ip=1,ngauss_pts(iw)
          vis_dir = vis_dir + DIRECTs(ipc,pverp)
          vis_dif = vis_dif + CK1s(ipc,pverp)*EL2s(ipc,pverp)+CK2s(ipc,pverp)*EM2s(ipc,pverp)+CMBs(ipc,pverp)
          ipc = ipc+1
        enddo
      else
        do ip=1,ngauss_pts(iw)
          nir_dir = nir_dir + DIRECTs(ipc,pverp)
          nir_dif = nir_dif + CK1s(ipc,pverp)*EL2s(ipc,pverp)+CK2s(ipc,pverp)*EM2s(ipc,pverp)+CMBs(ipc,pverp)
          ipc = ipc+1
        enddo
      endif
    enddo
  end subroutine rad_postcalc_accum

!============================================================================

  subroutine rad_postcalc_final (dzc, cpair_col, swcut, part_in_tshadow, sw_on, &
                                 sw_dTdt, lw_dTdt, lw_dnflux, lw_upflux, sw_upflux, sw_dnflux)
!------------------------------------------------------------------------
! heating rates from the accumulated fluxes (unchanged logic of
! the original rad_postcalc).
!------------------------------------------------------------------------
    implicit none
    real(r8), intent(in), dimension(pver) :: dzc
    real(r8), intent(in) :: cpair_col
    integer,  intent(in) :: swcut
    logical,  intent(in) :: part_in_tshadow
    logical,  intent(in) :: sw_on
    real(r8), intent(out), dimension(pver) :: sw_dTdt, lw_dTdt
    real(r8), intent(in), dimension(pverp) :: lw_dnflux, lw_upflux, sw_upflux, sw_dnflux
    integer :: k
    real(r8) :: lyr_mass_fact, fdiv_sw, fdiv_lw

    sw_dTdt(:) = 0.
    lw_dTdt(:) = 0.
    if (sw_on) then
      if (part_in_tshadow) then
        do k=camtop+1,swcut-1
          lyr_mass_fact = dzc(k-1)*cpair_col
          fdiv_sw = (sw_upflux(k)-sw_dnflux(k))-(sw_upflux(k-1)-sw_dnflux(k-1))
          sw_dTdt(k-1) = fdiv_sw/lyr_mass_fact
          fdiv_lw = (lw_upflux(k)-lw_dnflux(k))-(lw_upflux(k-1)-lw_dnflux(k-1))
          lw_dTdt(k-1) = fdiv_lw/lyr_mass_fact
        enddo
        do k=swcut,pverp
          lyr_mass_fact = dzc(k-1)*cpair_col
          sw_dTdt(k-1) = 0.0
          fdiv_lw = (lw_upflux(k)-lw_dnflux(k))-(lw_upflux(k-1)-lw_dnflux(k-1))
          lw_dTdt(k-1) = fdiv_lw/lyr_mass_fact
        enddo
      else
        do k=camtop+1,pverp
          lyr_mass_fact = dzc(k-1)*cpair_col
          fdiv_sw = (sw_upflux(k)-sw_dnflux(k))-(sw_upflux(k-1)-sw_dnflux(k-1))
          sw_dTdt(k-1) = fdiv_sw/lyr_mass_fact
          fdiv_lw = (lw_upflux(k)-lw_dnflux(k))-(lw_upflux(k-1)-lw_dnflux(k-1))
          lw_dTdt(k-1) = fdiv_lw/lyr_mass_fact
        enddo
      endif
    else
      do k=camtop+1, pverp
        lyr_mass_fact = dzc(k-1)*cpair_col
        fdiv_lw = (lw_upflux(k)-lw_dnflux(k))-(lw_upflux(k-1)-lw_dnflux(k-1))
        lw_dTdt(k-1) = fdiv_lw/lyr_mass_fact
      enddo
    endif
  end subroutine rad_postcalc_final

!============================================================================

  subroutine solve_block(pmid, tmid, tint, swcut, tau_gas, tau_ray, &
                         tau_cld_mcica, singscat_cld_mcica, asym_cld_mcica, &
                         tau_cld_co2, singscat_cld_co2, asym_cld_co2, &
                         tau_aer, wtau_aer, gwtau_aer, &
                         part_in_tshadow, sfc_albedo_dir, sfc_albedo_dif, &
                         sfc_emiss, sflux_frac, sfc_tempk, cos_mu, sw_on, &
                         lw_dnflux, lw_upflux, sw_upflux, sw_dnflux, &
                         lw_dnflux_spec, lw_upflux_spec, sw_upflux_spec, sw_dnflux_spec, &
                         vis_dir, vis_dif, nir_dir, nir_dif, sol_toa)
!------------------------------------------------------------------------
! the whole two-stream chain for one spectral block (bands
! bw0:bw1, g-points gb0:gb1) on block-sized work arrays.
!------------------------------------------------------------------------
    implicit none
    real(r8), intent(in), dimension(pverp) :: pmid, tmid, tint
    integer,  intent(in) :: swcut
    real(r8), intent(in), dimension(ntot_gpt,pverp) :: tau_gas
    real(r8), intent(in), dimension(ntot_wavlnrng,pverp) :: tau_ray
    real(r8), intent(in), dimension(ncld_grp,ntot_gpt,pverp) :: tau_cld_mcica, singscat_cld_mcica, asym_cld_mcica
    real(r8), intent(in), dimension(ntot_wavlnrng,pverp) :: tau_cld_co2, singscat_cld_co2, asym_cld_co2
    real(r8), intent(in), dimension(ntot_wavlnrng,pverp) :: tau_aer, wtau_aer, gwtau_aer
    logical,  intent(in) :: part_in_tshadow
    real(r8), intent(in), dimension(ntot_wavlnrng) :: sfc_albedo_dir, sfc_albedo_dif, sfc_emiss
    real(r8), intent(in) :: sflux_frac, sfc_tempk, cos_mu
    logical,  intent(in) :: sw_on
    real(r8), intent(inout), dimension(pverp) :: lw_dnflux, lw_upflux, sw_upflux, sw_dnflux
    real(r8), intent(inout), dimension(pverp,ntot_wavlnrng) :: lw_dnflux_spec, lw_upflux_spec
    real(r8), intent(inout), dimension(pverp,ntot_wavlnrng) :: sw_upflux_spec, sw_dnflux_spec
    real(r8), intent(inout) :: vis_dir, vis_dif, nir_dir, nir_dif, sol_toa

    ! block-sized work arrays (global g-point bounds)
    real(r8), dimension(gb0:gb1,ngangles,pverp) :: Y3
    real(r8), dimension(gb0:gb1,pverp) :: TAUL, OPD, PTEMP, SLOPE, SOL, W0, G0
    real(r8), dimension(gb0:gb1) :: PTEMPG, EMIS, RSFXdir, RSFXdif
    real(r8), dimension(gb0:gb1,pverp) :: EM1sol, EM2sol, EL1sol, EL2sol, AKsol, GAMIsol, &
                                          B1sol, B2sol, B3sol, EE1sol, CK1sol, CK2sol, &
                                          CPBsol, CMBsol, DIRECTsol
    real(r8), dimension(gb0:gb1,2*pverp) :: AFsol, BFsol, EFsol
    real(r8), dimension(gb0:gb1,pverp) :: EM1ir, EM2ir, EL1ir, EL2ir, AKir, GAMIir, &
                                          B1ir, B2ir, B3ir, EE1ir, CK1ir, CK2ir, &
                                          CPBir, CMBir, DIRECTir, DIRECTU, DIREC
    real(r8), dimension(gb0:gb1,2*pverp) :: AFir, BFir, EFir
    integer :: ip_ibeg, ip_iend
    logical :: beamSolar

    ! same zero state the original gave its full-size arrays on every call
    Y3 = 0.0 ; TAUL = 0.0 ; OPD = 0.0 ; PTEMP = 0.0 ; SLOPE = 0.0 ; SOL = 0.0 ; W0 = 0.0 ; G0 = 0.0
    PTEMPG = 0.0 ; EMIS = 0.0 ; RSFXdir = 0.0 ; RSFXdif = 0.0
    EM1sol = 0.0 ; EM2sol = 0.0 ; EL1sol = 0.0 ; EL2sol = 0.0 ; AKsol = 0.0 ; GAMIsol = 0.0
    B1sol = 0.0 ; B2sol = 0.0 ; B3sol = 0.0 ; EE1sol = 0.0 ; CK1sol = 0.0 ; CK2sol = 0.0
    CPBsol = 0.0 ; CMBsol = 0.0 ; DIRECTsol = 0.0 ; AFsol = 0.0 ; BFsol = 0.0 ; EFsol = 0.0
    EM1ir = 0.0 ; EM2ir = 0.0 ; EL1ir = 0.0 ; EL2ir = 0.0 ; AKir = 0.0 ; GAMIir = 0.0
    B1ir = 0.0 ; B2ir = 0.0 ; B3ir = 0.0 ; EE1ir = 0.0 ; CK1ir = 0.0 ; CK2ir = 0.0
    CPBir = 0.0 ; CMBir = 0.0 ; DIRECTir = 0.0 ; DIRECTU = 0.0 ; DIREC = 0.0
    AFir = 0.0 ; BFir = 0.0 ; EFir = 0.0

    call rad_precalc(pmid, tmid, tint, swcut, tau_gas, tau_ray, &
                     tau_cld_mcica, singscat_cld_mcica, asym_cld_mcica, &
                     tau_cld_co2, singscat_cld_co2, asym_cld_co2, &
                     tau_aer, wtau_aer, gwtau_aer, &
                     part_in_tshadow, sfc_albedo_dir, sfc_albedo_dif, sfc_emiss, &
                     sflux_frac, sfc_tempk, cos_mu, sw_on, &
                     Y3, TAUL, OPD, PTEMP, PTEMPG, SLOPE, SOL, W0, G0, EMIS, RSFXdir, RSFXdif)

    ! Solar pass: skipped on the night side, where its results were discarded
    if (sw_on) then
      ip_ibeg = max(sw_ipbeg,gb0) ; ip_iend = min(sw_ipend,gb1)
      if (ip_ibeg <= ip_iend) then
        beamSolar = .true.
        call two_stream(TAUL, W0, G0, EMIS, RSFXdir, RSFXdif, beamSolar, ip_ibeg, ip_iend, &
                        EM1sol, EM2sol, EL1sol, EL2sol, AFsol, BFsol, EFsol, AKsol, &
                        GAMIsol, B1sol, B2sol, EE1sol)
        call add_txrad(EM1sol, EM2sol, EL1sol, EL2sol, AFsol, BFsol, EFsol, B1sol, B2sol, AKsol, &
                       TAUL, OPD, PTEMP, PTEMPG, SLOPE, SOL, W0, G0, EMIS, RSFXdir, RSFXdif, cos_mu, sw_on, &
                       ip_ibeg, ip_iend, beamSolar, CK1sol, CK2sol, CPBsol, CMBsol, B3sol, DIRECTsol)
      endif
    endif

    ! Thermal pass
    ip_ibeg = max(lw_ipbeg,gb0) ; ip_iend = min(lw_ipend,gb1)
    if (ip_ibeg <= ip_iend) then
      beamSolar = .false.
      call two_stream(TAUL, W0, G0, EMIS, RSFXdir, RSFXdif, beamSolar, ip_ibeg, ip_iend, &
                      EM1ir, EM2ir, EL1ir, EL2ir, AFir, BFir, EFir, AKir, &
                      GAMIir, B1ir, B2ir, EE1ir)
      call add_txrad(EM1ir, EM2ir, EL1ir, EL2ir, AFir, BFir, EFir, B1ir, B2ir, AKir, &
                     TAUL, OPD, PTEMP, PTEMPG, SLOPE, SOL, W0, G0, EMIS, RSFXdir, RSFXdif, cos_mu, sw_on, &
                     ip_ibeg, ip_iend, beamSolar, CK1ir, CK2ir, CPBir, CMBir, B3ir, DIRECTir)
      call refine_lwflux(CK1ir, CK2ir, Y3, AKir, GAMIir, B3ir, EE1ir, &
                         TAUL, PTEMP, PTEMPG, SLOPE, EMIS, RSFXdir, RSFXdif, cos_mu, DIRECTU, DIREC, &
                         ip_ibeg, ip_iend)
    endif

    call rad_postcalc_accum(CK1sol, CK2sol, CPBsol, CMBsol, EM1sol, EM2sol, EL1sol, EL2sol, &
                            DIRECTsol, DIRECTU, DIREC, SOL, cos_mu, sw_on, &
                            lw_dnflux, lw_upflux, sw_upflux, sw_dnflux, &
                            lw_dnflux_spec, lw_upflux_spec, sw_upflux_spec, sw_dnflux_spec, &
                            vis_dir, vis_dif, nir_dir, nir_dif, sol_toa)
  end subroutine solve_block

!============================================================================


end module exo_radiation_mod
