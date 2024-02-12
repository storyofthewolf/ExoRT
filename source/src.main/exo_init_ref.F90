
module exo_init_ref

!---------------------------------------------------------------------
! Purpose:

  use shr_kind_mod,     only: r8 => shr_kind_r8
  use shr_const_mod,    only: SHR_CONST_PI, SHR_CONST_STEBOL
  use physconst,        only: scon, mwdry
  use radgrid
  use spmd_utils,       only: masterproc
  use planck_mod
  use exoplanet_mod

  implicit none

  public

  ! Approximate smallest double precision floating point difference
  !real(r8), parameter :: SMALLd = 1.0d-12
  !real(r8), parameter :: SMALLe = 1.0e-12
  real(r8), parameter :: SMALLd = 1.0d-8
  real(r8), parameter :: SMALLe = 1.0e-8

  real(r8), parameter :: sqrt3 = 1.732050808d0      ! square root of 3
  real(r8), parameter :: mb_to_atm = 9.869233e-4    ! convert pressure from Pa to atm

  !------------------------------------------------------------------------------
  ! Radiative transfer model variable/array declarations
  !

  ! Assign beginning and end wavelength range and point indices for each
  !  wavelength group
  integer :: lw_iwbeg = 1     ! thermal band wvl integration limits
  integer :: lw_iwend = ntot_wavlnrng
  integer :: sw_iwbeg = 1     ! solar band wvl integration limits
  integer :: sw_iwend = ntot_wavlnrng
  integer :: lw_ipbeg = 1     ! thermal band gpt integration limits
  integer :: lw_ipend = ntot_gpt
  integer :: sw_ipbeg = 1     ! solar band gpt integration limits
  integer :: sw_ipend = ntot_gpt

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


!============================================================================
  subroutine init_ref
!------------------------------------------------------------------------
! Purpose: Initial reference value arrays
!------------------------------------------------------------------------
  implicit none
!------------------------------------------------------------------------
! Local Variables
!

  integer :: iq, iw, ig, ip

!------------------------------------------------------------------------
!
! Start Code
!
    ! Arrange g_weight(ntot_pt) array
    iq = 0
    do iw=1,ntot_wavlnrng
      do ig=1, ngauss_pts(iw)
        iq = iq + 1
        if (ngauss_pts(iw) .eq. ngauss_8gpt) g_weight(iq) =  g_weight_8gpt(ig)
        if (ngauss_pts(iw) .eq. ngauss_16gpt) g_weight(iq) =  g_weight_16gpt(ig)
      enddo
    enddo

    ! Scale solar constant to namelist value
    solarflux(:) = solarflux(:)*scon/SUM(solarflux(:))

    !
    ! Calculate the "average" wavenumber (1/wavelength) <wavenum()> for each
    !  wavelength interval in each wavelength group [1/cm]:
    !
    ip = 0
    do iw=1,ntot_wavlnrng  ! "avg" wavenumber over each band
      do ig=1,ngauss_pts(iw)
        ip = ip+1
        ! Gauss-weighted solar flux in each probability interval:
        gw_solflux(ip) = solarflux(iw)*g_weight(ip)
      enddo
    enddo
    gw_solflux(:) = gw_solflux(:)*scon/SUM(gw_solflux(:))

    if (masterproc) then
      write (6, '(2x, a)') '_______________________________________________________'
      write (6, '(2x, a)') '_____________ stellar flux information ________________'
      write (6, '(2x, a)') '_______________________________________________________'
      write(*,*) "INIT_REF: total solar irradiance scaled to ",scon, "W m-2"
      write(*,*) "INPUT: solar flux [W m-2] in each spectral interval"
    endif

    if (do_exo_rt_optimize_bands) then
      call optimize_bands_sw
    else
      if (masterproc) then
        do iw=1, ntot_wavlnrng
          write(*,*) iw, solarflux(iw)
        enddo
      endif
    endif
    if (masterproc) then
      write(*,*) "TOTAL SOLAR FLUX:", SUM(solarflux), SUM(gw_solflux)
    endif
    if (do_exo_rt_optimize_bands) then
      call optimize_bands_lw
    endif

    if (masterproc) then
      write(*,*) "---------------------------------------"
      write(*,*) "SW intervals ", sw_iwbeg, sw_iwend
      write(*,*) "SW gauss pts ", sw_ipbeg, sw_ipend
      write(*,*) "LW intervals ", lw_iwbeg, lw_iwend
      write(*,*) "LW gauss pts ", lw_ipbeg, lw_ipend
      write(*,*) "---------------------------------------"
    endif

 end subroutine init_ref


!============================================================================
  subroutine optimize_bands_sw
!------------------------------------------------------------------------
! Purpose: Reduces the spectral integration limits in the shortwave
!          in order to improve model efficiency.
!          Lost flux rescaled into remaining bins.
!------------------------------------------------------------------------
  implicit none
!------------------------------------------------------------------------
! Local Variables
!
    real(r8) :: total_tsi, scale_ip, scale_iw
    real(r8), dimension(ntot_wavlnrng) :: cumulative_tsi
    real(r8), dimension(ntot_wavlnrng) :: solarflux_scale
    real(r8), dimension(ntot_gpt)      :: gw_solflux_scale
    integer :: iw, iwb, iwe, ip, temp, ig

!------------------------------------------------------------------------
!
! Start Code
!

    if (masterproc) write(*,*) "Optimizing shortwave bands..."

    ! set solar limits
    total_tsi = SUM(solarflux)
    do iw=1, ntot_wavlnrng
      cumulative_tsi(iw) = SUM(solarflux(iw:ntot_wavlnrng))/total_tsi
      !!write(*,*) iw, solarflux(iw), cumulative_tsi(iw)
    enddo

    iwe = ntot_wavlnrng
    do iwb=ntot_wavlnrng, 1, -1
      if (cumulative_tsi(iwb) .lt. 1.0-swFluxLimit)  iwe = iwe - 1
      if (cumulative_tsi(iwb) > swFluxLimit) exit
    enddo

    sw_iwbeg = iwb
    temp = 0
    do ig=1, iwb-1
      temp = temp + ngauss_pts(ig)
    enddo
    sw_ipbeg = temp + 1

    sw_iwend = iwe
    temp = sw_ipbeg - 1
    do ig=iwb+1,iwe+1
      temp = temp + ngauss_pts(ig)
    enddo
    sw_ipend = temp

    scale_ip = SUM(gw_solflux(:)) / SUM(gw_solflux(sw_ipbeg:sw_ipend))
    scale_iw = SUM(solarflux(:))  / SUM(solarflux(sw_iwbeg:sw_iwend))

    do iw=1, ntot_wavlnrng
      solarflux_scale(iw) = solarflux(iw)  * scale_iw
      if (iw .lt. sw_iwbeg) solarflux_scale(iw) = 0.0
      if (iw .gt. sw_iwend) solarflux_scale(iw) = 0.0
    enddo

    do ip=1, ntot_gpt
      gw_solflux_scale(ip) = gw_solflux(ip)  * scale_ip
      if (ip .lt. sw_ipbeg) gw_solflux_scale(ip) = 0.0
      if (ip .gt. sw_ipend) gw_solflux_scale(ip) = 0.0
    enddo

    if (masterproc) then
      write(*,*) "optimizing stellar radiation bands"
      write(*,*) "optimized to ", swfluxLimit
      write(*,*) "SW intervals ", sw_iwbeg, sw_iwend
      write(*,*) "SW gauss pts ", sw_ipbeg, sw_ipend
      write(*,*) "scale ",        scale_iw, scale_ip
      write(*,*) "          iw full                    scaled"
      do iw=1, ntot_wavlnrng
        write(*,*) iw, solarflux(iw), solarflux_scale(iw)
      enddo
      write(*,*) "before/after scaling", SUM(solarflux), SUM(solarflux_scale)
      write(*,*) "before/after scaling", SUM(gw_solflux), SUM(gw_solflux_scale)
    endif

    solarflux(:) = solarflux_scale(:)
    gw_solflux(:) = gw_solflux_scale(:)

  end subroutine optimize_bands_sw



!============================================================================
  subroutine optimize_bands_lw
!------------------------------------------------------------------------
! Purpose: Reduces the spectral integration limits in the longwave
!          in order to improve model efficiency.
!          No rescaling of lost flux.
!------------------------------------------------------------------------

  implicit none

!------------------------------------------------------------------------
! Local Variables
!
    real(r8) :: planck, w1, w2
    real(r8), dimension(ntot_wavlnrng) :: flux, cumulative_flux
    integer :: iw, iwr, ig, temp

!------------------------------------------------------------------------
!
! Start Code
!
    do iw=1,ntot_wavlnrng

       w1 = dble(1.439*wavenum_edge(iw))      ! 1.439 ~ ((h*c)/k)
       w2 = dble(1.439*wavenum_edge(iw+1))    ! convert nu cm-1 to lambda (um)
       planck = (PLANCKf(w1,Tmax)-PLANCKf(w2,Tmax))*SHR_CONST_STEBOL
       flux(iw) = planck
    enddo

    do iw=1, ntot_wavlnrng
      cumulative_flux(iw) = SUM(flux(1:iw))/SUM(flux)
    enddo
    do iwr=1, ntot_wavlnrng
      if (cumulative_flux(iwr) > lwFluxLimit) then
        exit
      endif
    enddo

    lw_iwend = iwr
    temp = 0
    do ig=1, iwr
      temp = temp + ngauss_pts(ig)
    enddo
    lw_ipend = temp

    if (masterproc) then
      write(*,*) "optimizing longwave radiation bands"
      write(*,*) "optimized to ", lwfluxLimit, Tmax, " K"
    endif

  end subroutine optimize_bands_lw


end module exo_init_ref
