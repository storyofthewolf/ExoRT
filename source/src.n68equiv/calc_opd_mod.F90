
module calc_opd_mod
! version n68equiv
!
! Equivalent Extinction:  For every grid box, first compare the grey optical depths.
! Which ever has the largest grey optical depth (k*u), use the k-distribution.
! the other gases will be treated using equivalent extinction.
! This method would allow for dynamic selection of major and minor species
! continuously during the simulation.

  use shr_kind_mod,     only: r8 => shr_kind_r8
  use rayleigh_data
  use shr_const_mod,    only: SHR_CONST_PI,SHR_CONST_PI, SHR_CONST_G, &
                              SHR_CONST_RGAS, SHR_CONST_AVOGAD, &
                              SHR_CONST_STEBOL, &
                              SHR_CONST_BOLTZ, &
                              SHR_CONST_RHOFW, SHR_CONST_RHOICE, &
                              SHR_CONST_LOSCHMIDT
  use physconst,        only: mwn2, mwco2, mwch4, mwc2h6, mwh2o, mwo2, mwh2, mwo3, mwdry, cpair, epsilo
  use radgrid
  use rad_interp_mod
  use ppgrid
  use kabs

  implicit none
  private
  save

  public :: calc_gasopd
  public :: calc_cldopd


!============================================================================
contains
!============================================================================

!============================================================================

  subroutine calc_gasopd(tmid, pmid, pdel, coldens, coldens_dry, qh2o, qco2, qch4, qc2h6, qO2, qO3, qH2, qN2, &
                         pathlength, tau_gas, tau_ray)

!------------------------------------------------------------------------
!
! Purpose: Calculate the optical depths of gases
!          Optical depths stored in tau_gas(g_value,wavelength_band,vertical_level)
!
!------------------------------------------------------------------------
!

    !use time_manager,   only: get_nstep

    implicit none

!------------------------------------------------------------------------
!
! Input Arguments
!
    real(r8), intent(in), dimension(pverp) :: tmid         ! temperatures at mid layers [K]
    real(r8), intent(in), dimension(pverp) :: pmid         ! pressure midlayers [mb]
    real(r8), intent(in), dimension(pver)  :: pdel         ! layer thickness [mb]
    real(r8), intent(in), dimension(pverp) :: coldens      ! Wet Column density profile [molec m-2]
    real(r8), intent(in), dimension(pverp) :: coldens_dry  ! Dry Column density profile [molec m-2]
    real(r8), intent(in), dimension(pverp) :: qh2o         ! mass mixing ratio h2o profile [kg/kg] wet
    real(r8), intent(in), dimension(pverp) :: qco2         ! mass mixing ratio co2 profile [kg/kg] dry
    real(r8), intent(in), dimension(pverp) :: qch4         ! mass mixing ratio ch4 profile [kg/kg] dry
    real(r8), intent(in), dimension(pverp) :: qc2h6        ! mass mixing ratio c2h6 profile [kg/kg] dry
    real(r8), intent(in), dimension(pverp) :: qo2          ! mass mixing ratio o2 profile [kg/kg] dry
    real(r8), intent(in), dimension(pverp) :: qo3          ! mass mixing ratio o3 profile [kg/kg] dry
    real(r8), intent(in), dimension(pverp) :: qh2          ! mass mixing ratio h2 profile [kg/kg] dry
    real(r8), intent(in), dimension(pverp) :: qn2          ! mass mixing ratio n2 profile [kg/kg] dry
    real(r8), intent(in), dimension(pverp) :: pathlength   ! thickness of layer [cm]
    real(r8), intent(out), dimension(ntot_gpt,pverp) ::  tau_gas
    real(r8), intent(out), dimension(ntot_wavlnrng,pverp) ::  tau_ray

!------------------------------------------------------------------------
!
! Local Variables
!
    ! layer pressure, temperature
    ! used for interpolation of k-coefficients
    real(r8) :: pressure
    real(r8) :: temperature

    ! gas volume mixing ratios
    real(r8) :: h2ovmr, co2vmr
    real(r8) :: ch4vmr, c2h6vmr
    real(r8) :: h2vmr, n2vmr
    real(r8) :: o2vmr, o3vmr

    ! indices for interpolation
    integer :: p_ref_index, t_ref_index, t_ref_index_mtckd
    integer :: t_ref_index_h2h2,t_ref_index_n2h2, t_ref_index_n2n2
    integer :: t_ref_index_co2co2_sw, t_ref_index_co2co2_lw
    integer :: t_ref_index_co2ch4, t_ref_index_co2h2
    integer :: ik, iw, ig, itu, itl, itc, sp  ! indices
    integer :: iwbeg, iwend  ! first and last band

    ! absorption coefficient "answers" returned from interpolators
    real(r8) :: ans_kmajor, ans_kgrey_h2o, ans_kgrey_co2, ans_kgrey_ch4, ans_kgrey_c2h6, ans_kgrey_o3, ans_kgrey_o2, ans
    real(r8), dimension(ngpt_max) :: ans_kmajor_gptvec
    real(r8), dimension(ntot_wavlnrng) :: ans_cia
    real(r8), dimension(ngauss_8gpt, ntot_wavlnrng) :: ans_h2os, ans_h2of

    ! Gas quantities
    real(r8) :: u_col, u_h2o, u_co2, u_ch4, u_c2h6, u_h2, u_n2, u_o2, u_o3
    real(r8), dimension(nspecies) :: ugas, tau_grey
    integer, dimension(1) :: imajor

    ! place holder temperatures for interpolation
    real(r8) :: t_kgas, t_n2n2, t_n2h2, t_h2h2, t_h2o_mtckd
    real(r8) :: t_co2ch4, t_co2h2, t_co2co2_sw, t_co2co2_lw

    real(r8) :: wm, wl, wla, r, ns, w

    ! for rayleigh scattering calc, depolarization
    real(r8) :: depolN2, depolCO2, depolH2O, depolCH4, depolO2
    ! for rayleigh scattering calc, Allen (1976) coefficients
    real(r8) :: allenN2, allenCO2, allenCH4, allenO2
    ! rayleigh scattering cross sections [cm2 molecule-1]
    real(r8) :: sigmaRayl, sigmaRaylCO2, sigmaRaylN2, sigmaRaylH2, sigmaRaylH2O, sigmaRaylCH4, sigmaRaylO2
    real(r8) :: kg_sw_minval  !! minimum value to check sw_abs error

    ! partial pressures
    real(r8) :: ppN2, ppH2, ppCO2, ppCH4, ppH2O

    ! amagats
    real(r8) :: amaN2, amaH2, amaCO2, amaCH4, amaH2O, amaFRGN

    ! CIA optical depths
    real(r8), dimension(ntot_wavlnrng,pverp) ::  tau_cia  ! total
    ! individual
    real(r8), dimension(ntot_wavlnrng,pverp) ::  tau_n2n2_cia
    real(r8), dimension(ntot_wavlnrng,pverp) ::  tau_n2h2_cia
    real(r8), dimension(ntot_wavlnrng,pverp) ::  tau_h2h2_cia
    real(r8), dimension(ntot_wavlnrng,pverp) ::  tau_co2co2_lw_cia
    real(r8), dimension(ntot_wavlnrng,pverp) ::  tau_co2co2_sw_cia
    real(r8), dimension(ntot_wavlnrng,pverp) ::  tau_co2h2_cia
    real(r8), dimension(ntot_wavlnrng,pverp) ::  tau_co2ch4_cia


!------------------------------------------------------------------------
!
! Start Code
!

    !
    ! calculate optical depth over all spectral intervals
    iwbeg = 1
    iwend = ntot_wavlnrng

    ! initialize some things to zero, just to be safe
    tau_gas(:,:) = 0.0
    tau_ray(:,:) = 0.0
    tau_cia(:,:) = 0.0
    tau_n2n2_cia(:,:) = 0.0
    tau_n2h2_cia(:,:) = 0.0
    tau_h2h2_cia(:,:) = 0.0
    tau_co2co2_lw_cia(:,:) = 0.0
    tau_co2co2_sw_cia(:,:) = 0.0
    tau_co2h2_cia(:,:) = 0.0
    tau_co2ch4_cia(:,:) = 0.0

    do ik = 1,pverp

      ! optical depths are calculated at mid-layers
      ! level pverp is the atmosphere grid box nearest the surface
      ! level 1 is the psuedo layer above the model top to infinity

!write(*,*) "------- level ------ ", ik, " -----------------------"

      !! Define grid box volume mixing ratios and column densities

      ! qco2, qh2o, qch4 are in mass mixing ratio, multiply by mwdry/mw*** to get volume mixing ratio
      ! factor of 10^-4 converts molecules m-2 to molecules cm-2
      ! qh2o is defined relative to wet atm mass
      ! all others are defined relative to dry atm mass
      ! set all vmr quantities relative to dry mass except where noted
      w = qh2o(ik)/(1.0-qh2o(ik))             ! H2O mass mixing ratio relative to dry air
      h2ovmr  = w*mwdry/mwh2o                 ! H2O dry volume mixing ratio
      co2vmr  = qco2(ik)*mwdry/mwco2          ! CO2 volume mixing ratio dry
      ch4vmr  = qch4(ik)*mwdry/mwch4          ! CH4 volume mixing ratio dry
      c2h6vmr = qc2h6(ik)*mwdry/mwc2h6        ! C2H6 volume mixing ratio dry
      o2vmr   = qo2(ik)*mwdry/mwo2            ! O2 volume mixing ratio dry
      o3vmr   = qo3(ik)*mwdry/mwo3            ! O3 volume mixing ratio dry
      h2vmr   = qh2(ik)*mwdry/mwh2            ! H2 volume mixing ratio dry
      n2vmr   = qn2(ik)*mwdry/mwn2            ! N2 volume mixing ratio dry

      u_h2o = h2ovmr*coldens_dry(ik)/10000.       !   water column amount [ molecules cm-2 ]
      u_co2 = co2vmr*coldens_dry(ik)/10000.       !   co2 column amount [ molecules cm-2 ]
      u_ch4 = ch4vmr*coldens_dry(ik)/10000.       !   ch4 column amount [ molecules cm-2 ]
      u_c2h6 = c2h6vmr*coldens_dry(ik)/10000.     !   c2h6 column amount [ molecules cm-2 ]
      u_o2 = o2vmr*coldens_dry(ik)/10000.         !   o2 column amount [ molecules cm-2 ]
      u_o3 = o3vmr*coldens_dry(ik)/10000.         !   o3 column amount [ molecules cm-2 ]
      u_h2 = h2vmr*coldens_dry(ik)/10000.         !   h2 column amount [ molecules cm-2 ]
      u_n2 = n2vmr*coldens_dry(ik)/10000.         !   n2 column amount [ molecules cm-2 ]

      ! determine partial pressures (mb)
      ppH2O = h2ovmr/(1.+h2ovmr)*pmid(ik)
      ppN2  = n2vmr*(pmid(ik))/(1.+h2ovmr)
      ppH2  = h2vmr*(pmid(ik))/(1.+h2ovmr)
      ppCO2 = co2vmr*(pmid(ik))/(1.+h2ovmr)
      ppCH4 = ch4vmr*(pmid(ik))/(1.+h2ovmr)

      ! calculate amagats of various gases
      amaN2   = (273.15/tmid(ik)) * (ppN2/1013.25)
      amaH2   = (273.15/tmid(ik)) * (ppH2/1013.25)
      amaCO2  = (273.15/tmid(ik)) * (ppCO2/1013.25)
      amaCH4  = (273.15/tmid(ik)) * (ppCH4/1013.25)
      amaH2O  = (273.15/tmid(ik)) * (ppH2O/1013.25)
      amaFRGN = (273.15/tmid(ik)) * ((pmid(ik)-ppH2O)/1013.25)

      ! create array of major gases
      ugas = (/u_h2o, u_co2, u_ch4, u_c2h6, u_o3, u_o2/)

      ! Find pressure coordinates for k-coefficients
      pressure = log10(pmid(ik))       ! log pressure
      ! find the reference pressure value, exit if pressure less than minimum of grid
      p_ref_index = kc_npress          ! index of reference pressure, default is max
      do  ! for K coefficient data sets
        if (p_ref_index .le. 1) exit
        if (log10pgrid(p_ref_index) .le. pressure ) exit
        p_ref_index = p_ref_index - 1
      enddo
      ! if pressure less than minimum of grid, force reference to minimum grid value
      ! force reference pressure for interpolation to minimum pressure in pgrid
      if (pressure .le. log10pgrid(1)) then
        p_ref_index = 1
        pressure = log10pgrid(p_ref_index)
      endif


      ! Find temperature coordinates for k-coefficients
      temperature = tmid(ik)              ! actual temperature [K]
      ! For gas absorption k coefficients
      t_ref_index = kc_ntemp  ! index of reference temperature
      t_kgas = temperature
      do
        if (t_ref_index .le. 1) exit                   ! exit if temperature less than minimum grid
        if (t_kgas .gt. tgrid(kc_ntemp)) then  ! temperature greater than grid max
          t_kgas = tgrid(t_ref_index)   ! set t to max grid value
          exit                                              ! exit
        endif
        if ((tgrid(t_ref_index) .le. t_kgas)) exit ! found reference temperature
        t_ref_index = t_ref_index - 1                   ! increment reference temperature
      enddo
      ! if temperature less than minimum of grid, force reference to minimum grid value
      ! force reference temperature for interpolation to minimum temperature in tgrid
      if (t_ref_index .lt. 1) then
        t_ref_index = 1
        t_kgas = tgrid(t_ref_index)
      endif

      !
      !  Calculate gas absorption using equivalent extinction
      !
      !  First gray optical depths are computed for each gas.  The most optically thick
      !  gas is designated as the "major" absorber, and treated with a full k-distribution.
      !  All other species are "minor" and included using their gray absorption.

      itl = 0
      do sp=1, ntot_wavlnrng
        tau_grey(:) = 0.0
        call bilinear_interpK_grey(k_grey_data, iH2O,   sp, pressure, p_ref_index, t_kgas, t_ref_index, ans_kgrey_h2o)
        call bilinear_interpK_grey(k_grey_data, iCO2,   sp, pressure, p_ref_index, t_kgas, t_ref_index, ans_kgrey_co2)
        call bilinear_interpK_grey(k_grey_data, iCH4,   sp, pressure, p_ref_index, t_kgas, t_ref_index, ans_kgrey_ch4)
        call bilinear_interpK_grey(k_grey_data, iC2H6,  sp, pressure, p_ref_index, t_kgas, t_ref_index, ans_kgrey_c2h6)
        call bilinear_interpK_grey(k_grey_data, iO3,    sp, pressure, p_ref_index, t_kgas, t_ref_index, ans_kgrey_o3)
        call bilinear_interpK_grey(k_grey_data, iO2,    sp, pressure, p_ref_index, t_kgas, t_ref_index, ans_kgrey_o2)
        tau_grey(iH2O)  = ans_kgrey_h2o * ugas(iH2O)
        tau_grey(iCO2)  = ans_kgrey_co2 * ugas(iCO2)
        tau_grey(iCH4)  = ans_kgrey_ch4 * ugas(iCH4)
        tau_grey(iC2H6) = ans_kgrey_c2h6 * ugas(iC2H6)
        tau_grey(iO3)   = ans_kgrey_o3 * ugas(iO3)
        tau_grey(iO2)   = ans_kgrey_o2 * ugas(iO2)
        imajor = maxloc(tau_grey)
        !write(*,*) ik, sp, gas_name(imajor), tau_grey, ans_kgrey_h2o, ans_kgrey_co2, ans_kgrey_ch4, ans_kgrey_c2h6
        ! major gas (k-distribution)
        call bilinear_interpK_8gpt_major_gptvec(k_major_data, imajor(1), sp, pressure, p_ref_index, t_kgas, t_ref_index, ans_kmajor_gptvec)
        !!write(*,*) "kmajor", ans_kmajor_gptvec
        do ig = 1,ngauss_pts(sp)
          itl = itl + 1
          tau_gas(itl,ik) = ans_kmajor_gptvec(ig)*ugas(imajor(1)) + (SUM(tau_grey) - tau_grey(imajor(1)))
        enddo
      enddo

      !
      !  Water Vapor Self Continuum
      !
      !===== MT_CKD =====!

      ! For water vapor self continuum, find temperature index
      ! find the reference temperature value, exit if temperature less than minimum of grid
      t_ref_index_mtckd = kmtckd_ntemp  ! index of reference temperature
      t_h2o_mtckd = temperature
      do
        if (t_ref_index_mtckd .le. 1) exit                   ! exit if temperature less than minimum grid
        if (t_h2o_mtckd .gt. tgrid_mtckd(kmtckd_ntemp)) then  ! temperature greater than grid max
          t_h2o_mtckd = tgrid_mtckd(t_ref_index_mtckd)   ! set t to max grid value
          exit                                              ! exit
        endif
        if ((tgrid_mtckd(t_ref_index_mtckd) .le. t_h2o_mtckd)) exit ! found reference temperature
        t_ref_index_mtckd = t_ref_index_mtckd - 1                   ! increment reference temperature
      enddo
      ! if temperature less than minimum of grid, force reference to minimum grid value
      ! force reference temperature for interpolation to minimum temperature in tgrid
      if (t_ref_index_mtckd .lt. 1) then
        t_ref_index_mtckd = 1
        t_h2o_mtckd = tgrid_mtckd(t_ref_index_mtckd)
      endif

      !!  Do MT_CKD continuum
      itc=0
      call interpH2Omtckd_ng(kh2oself_mtckd, t_h2o_mtckd, t_ref_index_mtckd, ans_h2os)
      call interpH2Omtckd_ng(kh2ofrgn_mtckd, t_h2o_mtckd, t_ref_index_mtckd, ans_h2of)
      do iw = iwbeg,iwend     ! loop over bands
        do ig=1, ngauss_pts(iw)
          itc=itc+1
          tau_gas(itc,ik) = tau_gas(itc,ik) + (ans_h2os(ig,iw)*amaH2O + ans_h2of(ig,iw)*amaFRGN) * u_h2o
         enddo
      enddo    ! close band loop


      !
      !  Collision Induced Absorption
      !
      !!====  Calculate N2-N2 CIA  ====!!
      if (u_n2 .gt. 0) then
        t_ref_index_n2n2 = kn2n2_ntemp
        t_n2n2 = temperature
        do
          if (t_ref_index_n2n2 .le. 1) exit                   ! exit if temperature less than minimum grid
          if (t_n2n2 .gt. tgrid_n2n2(kn2n2_ntemp)) then  ! temperature greater than grid max
            t_n2n2 = tgrid_n2n2(t_ref_index_n2n2)   ! set t to max grid value
            exit                                              ! exit
          endif
          if ((tgrid_n2n2(t_ref_index_n2n2) .le. t_n2n2)) exit ! found reference temperature
          t_ref_index_n2n2 = t_ref_index_n2n2 - 1                   ! increment reference temperature
        enddo
        ! if temperature less than minimum of grid, force reference to minimum grid value
        ! force reference temperature for interpolation to minimum temperature in tgrid
        if (t_ref_index_n2n2 .lt. 1) then
          t_ref_index_n2n2 = 1
          t_n2n2 = tgrid_n2n2(t_ref_index_n2n2)
        endif
        call interpN2N2cia(kn2n2, t_n2n2, t_ref_index_n2n2, ans_cia)
        do iw=iwbeg,iwend      ! loop over bands
          tau_n2n2_cia(iw,ik) = ans_cia(iw) * amaN2 * amaN2 * pathlength(ik)
          !!!write(*,*) "N2-N2 CIA",iw, ans, n2vmr, tau_n2n2cia(iw,ik)
        enddo
      endif


      !!====  Calculate H2-H2 CIA  ====!!
      if (u_h2 .gt. 0) then
        t_ref_index_h2h2 = kh2h2_ntemp
        t_h2h2 = temperature
        do
          if (t_ref_index_h2h2 .le. 1) exit                   ! exit if temperature less than minimum grid
          if (t_h2h2 .gt. tgrid_h2h2(kh2h2_ntemp)) then  ! temperature greater than grid max
            t_h2h2 = tgrid_h2h2(t_ref_index_h2h2)   ! set t to max grid value
            exit                                              ! exit
          endif
          if ((tgrid_h2h2(t_ref_index_h2h2) .le. t_h2h2)) exit ! found reference temperature
          t_ref_index_h2h2 = t_ref_index_h2h2 - 1                   ! increment reference temperature
        enddo
        ! if temperature less than minimum of grid, force reference to minimum grid value
        ! force reference temperature for interpolation to minimum temperature in tgrid
        if (t_ref_index_h2h2 .lt. 1) then
          t_ref_index_h2h2 = 1
          t_h2h2 = tgrid_h2h2(t_ref_index_h2h2)
        endif
        call interpH2H2cia(kh2h2, t_h2h2, t_ref_index_h2h2, ans_cia)
        do iw=iwbeg,iwend      ! loop over bands
           tau_h2h2_cia(iw,ik) = ans_cia(iw) * amaH2 * amaH2 * pathlength(ik)
          !!!write(*,*) "H2-H2 CIA",iw, ans, h2vmr, tau_h2h2cia(iw,ik)
        enddo
      endif


      !!====  Calculate N2-H2 CIA  ====!!
      if (u_h2 .gt. 0 .and. u_n2 .gt. 0) then
        t_ref_index_n2h2 = kn2h2_ntemp
        t_n2h2 = temperature
        do
          if (t_ref_index_n2h2 .le. 1) exit                   ! exit if temperature less than minimum grid
          if (t_n2h2 .gt. tgrid_n2h2(kn2h2_ntemp)) then  ! temperature greater than grid max
            t_n2h2 = tgrid_n2h2(t_ref_index_n2h2)   ! set t to max grid value
            exit                                              ! exit
          endif
          if ((tgrid_n2h2(t_ref_index_n2h2) .le. t_n2h2)) exit ! found reference temperature
          t_ref_index_n2h2 = t_ref_index_n2h2 - 1                   ! increment reference temperature
        enddo
        ! if temperature less than minimum of grid, force reference to minimum grid value
        ! force reference temperature for interpolation to minimum temperature in tgrid
        if (t_ref_index_n2h2 .lt. 1) then
          t_ref_index_n2h2 = 1
          t_n2h2 = tgrid_n2h2(t_ref_index_n2h2)
        endif
        call interpN2H2cia(kn2h2, t_n2h2, t_ref_index_n2h2, ans_cia)
        do iw=iwbeg,iwend      ! loop over bands
          tau_n2h2_cia(iw,ik)  = ans_cia(iw) * amaN2 * amaH2 * pathlength(ik)
          !!!write(*,*) "N2-H2 CIA",iw, ans, h2vmr !tau_n2h2cia(iw,ik)
        enddo
      endif

     !!====  Calculate CO2-CO2 CIA  ====!!
     if (u_co2 .gt. 0) then
       t_ref_index_co2co2_lw = kco2co2_lw_ntemp
       t_co2co2_lw = temperature
       do
         if (t_ref_index_co2co2_lw .le. 1) exit                   ! exit if temperature less than minimum grid
         if (t_co2co2_lw .gt. tgrid_co2co2_lw(kco2co2_lw_ntemp)) then  ! temperature greater than grid max
           t_co2co2_lw = tgrid_co2co2_lw(t_ref_index_co2co2_lw)   ! set t to max grid value
           exit                                              ! exit
         endif
         if ((tgrid_co2co2_lw(t_ref_index_co2co2_lw) .le. t_co2co2_lw)) exit ! found reference temperature
         t_ref_index_co2co2_lw = t_ref_index_co2co2_lw - 1                   ! increment reference temperature
       enddo
       ! if temperature less than minimum of grid, force reference to minimum grid value
       ! force reference temperature for interpolation to minimum temperature in tgrid
       if (t_ref_index_co2co2_lw .lt. 1) then
         t_ref_index_co2co2_lw = 1
         t_co2co2_lw = tgrid_co2co2_lw(t_ref_index_co2co2_lw)
       endif
       call interpCO2CO2cia_lw(kco2co2_lw, t_co2co2_lw, t_ref_index_co2co2_lw, ans_cia)
       do iw=iwbeg,iwend      ! loop over bands
         tau_co2co2_lw_cia(iw,ik) = ans_cia(iw) * amaCO2 * amaCO2 * pathlength(ik)
         !!!write(*,*) "CO2-CO2 CIA",iw, ans_cia(iw), amaCO2, pathlength(ik),  tau_co2co2_lw_cia(iw,ik)
       enddo

       t_ref_index_co2co2_sw = kco2co2_sw_ntemp
       t_co2co2_sw = temperature
       do
         if (t_ref_index_co2co2_sw .le. 1) exit                   ! exit if temperature less than minimum grid
         if (t_co2co2_sw .gt. tgrid_co2co2_sw(kco2co2_sw_ntemp)) then  ! temperature greater than grid max
           t_co2co2_sw = tgrid_co2co2_sw(t_ref_index_co2co2_sw)   ! set t to max grid value
           exit                                              ! exit
         endif
         if ((tgrid_co2co2_sw(t_ref_index_co2co2_sw) .le. t_co2co2_sw)) exit ! found reference temperature
         t_ref_index_co2co2_sw = t_ref_index_co2co2_sw - 1                   ! increment reference temperature
       enddo
       ! if temperature less than minimum of grid, force reference to minimum grid value
       ! force reference temperature for interpolation to minimum temperature in tgrid
       if (t_ref_index_co2co2_sw .lt. 1) then
         t_ref_index_co2co2_sw = 1
         t_co2co2_sw = tgrid_co2co2_sw(t_ref_index_co2co2_sw)
       endif
       call interpCO2CO2cia_sw(kco2co2_sw, t_co2co2_sw, t_ref_index_co2co2_sw, ans_cia)
       do iw=iwbeg,iwend      ! loop over bands
         tau_co2co2_sw_cia(iw,ik) = ans_cia(iw) * amaCO2 * amaCO2 * pathlength(ik)
         !!!write(*,*) "CO2-CO2 CIA",iw, ans, co2vmr, tau_co2co2cia(iw,ik)
       enddo
     endif

     !!====  Calculate CO2-CH4 CIA  ====!!
     if (u_co2 .gt. 0 .and. u_ch4 .gt. 0) then
        t_ref_index_co2ch4 = kco2ch4_ntemp
        t_co2ch4 = temperature
        do
          if (t_ref_index_co2ch4 .le. 1) exit                   ! exit if temperature less than minimum grid
          if (t_co2ch4 .gt. tgrid_co2ch4(kco2ch4_ntemp)) then  ! temperature greater than grid max
            t_co2ch4 = tgrid_co2ch4(t_ref_index_co2ch4)   ! set t to max grid value
            exit                                              ! exit
          endif
          if ((tgrid_co2ch4(t_ref_index_co2ch4) .le. t_co2ch4)) exit ! found reference temperature
          t_ref_index_co2ch4 = t_ref_index_co2ch4 - 1                   ! increment reference temperature
        enddo
        ! if temperature less than minimum of grid, force reference to minimum grid value
        ! force reference temperature for interpolation to minimum temperature in tgrid
        if (t_ref_index_co2ch4 .lt. 1) then
          t_ref_index_co2ch4 = 1
          t_co2ch4 = tgrid_co2ch4(t_ref_index_co2ch4)
        endif
        t_ref_index_co2ch4    = kco2ch4_ntemp
        call interpCO2CH4cia(kco2ch4, t_co2ch4, t_ref_index_co2ch4, ans_cia)
        do iw=iwbeg,iwend      ! loop over bands
          tau_co2ch4_cia(iw,ik) = ans_cia(iw) * amaCO2 * amaCH4 * pathlength(ik)
!          !write(*,*) "CO2-CH4 CIA",iw, ans, amaCO2, pathlength(ik),  tau_co2ch4cia(iw,ik)
        enddo
      endif

      !!====  Calculate CO2-H2 CIA  ====!!
      if (u_co2 .gt. 0 .and. u_h2 .gt. 0) then
        t_ref_index_co2h2 = kco2h2_ntemp
        t_co2h2 = temperature
        do
          if (t_ref_index_co2h2 .le. 1) exit                   ! exit if temperature less than minimum grid
          if (t_co2h2 .gt. tgrid_co2h2(kco2h2_ntemp)) then  ! temperature greater than grid max
            t_co2h2 = tgrid_co2h2(t_ref_index_co2h2)   ! set t to max grid value
            exit                                              ! exit
          endif
          if ((tgrid_co2h2(t_ref_index_co2h2) .le. t_co2h2)) exit ! found reference temperature
          t_ref_index_co2h2 = t_ref_index_co2h2 - 1                   ! increment reference temperature
        enddo
        ! if temperature less than minimum of grid, force reference to minimum grid value
        ! force reference temperature for interpolation to minimum temperature in tgrid
        if (t_ref_index_co2h2 .lt. 1) then
          t_ref_index_co2h2 = 1
          t_co2h2 = tgrid_co2h2(t_ref_index_co2h2)
        endif
        call interpCO2H2cia(kco2h2, t_co2h2, t_ref_index_co2h2, ans_cia)
        do iw=iwbeg,iwend      ! loop over bands
          tau_co2h2_cia(iw,ik) = ans_cia(iw) * amaCO2 * amaH2 * pathlength(ik)
          !!!write(*,*) "CO2-H2 CIA",iw, ans, co2vmr, tau_co2ch2cia(iw,ik)
        enddo
      endif

      !
      ! Add CIA optical depths to total optical depth
      !
      itc = 0
      do iw=iwbeg,iwend      ! loop over bands
         do ig=1, ngauss_pts(iw)
           itc = itc + 1
            tau_gas(itc, ik) = tau_gas(itc, ik) + tau_n2n2_cia(iw,ik) + tau_n2h2_cia(iw,ik) + tau_h2h2_cia(iw,ik)  &
                                               +  tau_co2ch4_cia(iw,ik) + tau_co2h2_cia(iw,ik)  &
                                               +  tau_co2co2_sw_cia(iw,ik) + tau_co2co2_lw_cia(iw,ik)
         enddo
       enddo


      !
      !  Calculate Rayleigh scattering optical depth
      !
      do iw = iwbeg,iwend     ! loop over bands

        wm = (wavenum_edge(iw) + wavenum_edge(iw+1))/2.0    ! wavenumber at bin midpoint
        wl = 1.e4/wm  ! wavelength in microns
        wla = wl*1.0e4 ! wavelength in angstroms

        !
        ! Vardavas and Carver (1984), Allen (1976) N2, CO2 Rayleigh scattering
        !
        ! Rayleigh scattering for CO2, N2
        depolCO2 = (6+3*delCO2)/(6-7*delCO2)
        depolN2  = (6+3*delN2)/(6-7*delN2)
        depolCH4 = (6+3*delCH4)/(6-7*delCH4)
        depolO2  = (6+3*delO2)/(6-7*delO2)
        allenCO2 = (1.0E-5*raylA_CO2*(1.0+1.0E-3*raylB_CO2/wl**2))**2
        allenN2 = (1.0E-5*raylA_N2*(1.0+1.0E-3*raylB_N2/wl**2))**2
        allenO2 = (1.0E-5*raylA_O2*(1.0+1.0E-3*raylB_O2/wl**2))**2
        allenCH4 = ((4869.8 + 4.1023e14/(1.133e10 - wm**2))*1.0e-8)**2   ! He et al. 2021, doi.org/10.5194/acp-21-14927-2021
        sigmaRaylCO2 = 4.577E-21/wl**4*depolCO2*allenCO2
        sigmaRaylN2 = 4.577E-21/wl**4*depolN2*allenN2
        sigmaRaylO2 = 4.577E-21/wl**4*depolO2*allenO2
        sigmaRaylCH4 = 4.577E-21/wl**4*depolCH4*allenCH4
        !  Rayleigh scattering from H2O
        ns = (5791817./(238.0185-(1./wl)**2) + 167909./(57.362-(1./wl)**2))/1.0E8  ! Bucholtz (1995)  !new
        r = 0.85*ns
        depolH2O = (6+3*delH2O)/(6-7*delH2O)
        sigmaRaylH2O = 4.577e-21*depolH2O*(r**2)/(wl**4)  !new
	! Rayleigh scattering from H2, Dalgarno & Williams 1962, ApJ, 136, 690D
      	sigmaRaylH2 = 8.14e-13/(wla**4) + 1.28e-6/(wla**6) + 1.61/(wla**8)
        ! Total Rayleigh scattering
        tau_ray(iw,ik) = sigmaRaylCO2*u_co2 + sigmaRaylN2*u_n2 + sigmaRaylH2O*u_h2o + sigmaRaylH2*u_h2 + sigmaRaylCH4*u_ch4 + sigmaRaylO2*u_o2
      enddo  ! close band loop

    enddo  ! close level loop

!    close(10)
    return

  end subroutine calc_gasopd



!============================================================================

  subroutine calc_cldopd(ext_pmid, ext_cICE, ext_cLIQ, ext_REI, ext_REL, ext_cFRC, &
                         tau_cld_mcica, singscat_cld_mcica, asym_cld_mcica, cFRC_mcica, &
                         cICE_mcica, cLIQ_mcica )

!------------------------------------------------------------------------
!
!  Purpose: Calculate the optical depth and related properties of clouds.
!  Uses MIE optical properties for both ice clouds and liquid water clouds.
!  Uses Monte Carlo Independent Column Approximation (MCICA) to treat cloud
!    overlap. Pincus (2003), Raisanen (2004)
!------------------------------------------------------------------------

   use mcica,            only: mcica_subcol
   use time_manager,     only: get_nstep

   implicit none

!------------------------------------------------------------------------
!
! Input Arguments
!
    real(r8), intent(in), dimension(pverp)  ::  ext_pmid      ! [Pa] interface pressure
    real(r8), intent(in), dimension(pverp)  ::  ext_cICE   ! [g/m2]  ice
    real(r8), intent(in), dimension(pverp)  ::  ext_cLIQ   ! [g/m2]  liquid
    real(r8), intent(in), dimension(pverp)  ::  ext_REI    ! [microns]  ice radii
    real(r8), intent(in), dimension(pverp)  ::  ext_REL    ! [microns]  liquid
    real(r8), intent(in), dimension(pverp)  ::  ext_cFRC   ! cloud fraction

    ! output bulk cloud optical properties after MCICA
    real(r8), intent(out), dimension(ncld_grp, ntot_gpt, pverp)  ::  tau_cld_mcica       ! cloud optical depth
    real(r8), intent(out), dimension(ncld_grp, ntot_gpt, pverp)  ::  singscat_cld_mcica  ! cloud single scattering albedo
    real(r8), intent(out), dimension(ncld_grp, ntot_gpt, pverp)  ::  asym_cld_mcica      ! cloud asymmetry parameter

    real(r8), intent(out), dimension(ntot_gpt, pverp) :: cFRC_mcica  ! stochastic cloud fraction
    real(r8), intent(out), dimension(ntot_gpt, pverp) :: cICE_mcica  ! stochastic ice cloud water content
    real(r8), intent(out), dimension(ntot_gpt, pverp) :: cLIQ_mcica  ! stochastic liquid cloud water content

!------------------------------------------------------------------------
!
! Local Variables
!
    real(r8), dimension(ncld_grp, ntot_wavlnrng, pverp)  ::  tau_cld_temp       ! cloud optical depth temporary array
    real(r8), dimension(ncld_grp, ntot_wavlnrng, pverp)  ::  singscat_cld_temp  ! cloud single scattering albedo temporary array
    real(r8), dimension(ncld_grp, ntot_wavlnrng, pverp)  ::  asym_cld_temp      ! cloud asymmetry parameter temporary array
    real(r8), dimension(ntot_gpt, pverp) :: tau_mcica_ice
    real(r8), dimension(ntot_gpt, pverp) :: ssa_mcica_ice
    real(r8), dimension(ntot_gpt, pverp) :: asym_mcica_ice
    real(r8), dimension(ntot_gpt, pverp) :: tau_mcica_liq
    real(r8), dimension(ntot_gpt, pverp) :: ssa_mcica_liq
    real(r8), dimension(ntot_gpt, pverp) :: asym_mcica_liq

    real(r8) :: Qliq
    real(r8) :: Wliq
    real(r8) :: Gliq

    real(r8) :: Qice
    real(r8) :: Wice
    real(r8) :: Gice

    real(r8) :: rho_liq
    real(r8) :: rho_ice
    real(r8) :: r_liq
    real(r8) :: r_ice

    integer :: iwend
    integer :: iwbeg
    integer :: iw
    integer :: ik
    integer :: ig

    integer :: icldovr
    integer :: permuteseed
    integer :: liqcld
    integer :: icecld

    integer :: nstep

!------------------------------------------------------------------------
!
! Start Code
!

    rho_liq = SHR_CONST_RHOFW*1000.0   ! density of fresh water [g m-3]
    rho_ice = SHR_CONST_RHOICE*1000.0  ! density of ice water [g m-3]

    ! Calculate gas optical depths over all spectral intervals
    iwbeg = 1
    iwend = ntot_wavlnrng

    liqcld = 1
    icecld = 2

    !initialize cloud optics array
    tau_cld_mcica(:,:,:) = 0.0
    singscat_cld_mcica(:,:,:) = 0.0
    asym_cld_mcica(:,:,:) = 0.0
    tau_cld_temp(:,:,:) = 0.0
    singscat_cld_temp(:,:,:) = 0.0
    asym_cld_temp(:,:,:) = 0.0

    !
    ! Water clouds
    !
    do ik=1, pverp
      do iw=iwbeg, iwend

      r_liq = ext_REL(ik) * 1.0e-6   ! liquid cloud drop size [m]

        if (ext_cLIQ(ik) .le. cldmin) then
          tau_cld_temp(liqcld,iw,ik) = 0.0
          singscat_cld_temp(liqcld,iw,ik) = 0.0
          asym_cld_temp(liqcld,iw,ik) = 0.0
        else
          call interpolate_cld(liqcld, iw, r_liq*1.0e6, Qliq, Wliq, Gliq, Qcldliq, Qcldice, Wcldliq, Wcldice, Gcldliq, Gcldice)
          !!!write(*,*) "interpolate_cld, liquid: r,q,w,g ",ext_REL(ik), Qliq, Wliq, Gliq, ext_cLIQ(ik)
          tau_cld_temp(liqcld,iw,ik) = 3.*Qliq  / (4.*rho_liq*r_liq) * ext_cLIQ(ik)
          singscat_cld_temp(liqcld,iw,ik) = Wliq
          asym_cld_temp(liqcld,iw,ik) = Gliq
        endif

      enddo
    enddo

    !
    ! Ice clouds
    !
    do ik=1, pverp
      do iw=iwbeg, iwend

        r_ice = ext_REI(ik) * 1.0e-6   ! ice cloud particle size [m]

        if (ext_cICE(ik) .le. cldmin) then
          tau_cld_temp(icecld,iw,ik) = 0.0
          singscat_cld_temp(icecld,iw,ik) = 0.0
          asym_cld_temp(icecld,iw,ik) = 0.0
        else
          call interpolate_cld(icecld, iw, ext_REI(ik), Qice, Wice, Gice, Qcldliq, Qcldice, Wcldliq, Wcldice, Gcldliq, Gcldice)
          !!!write(*,*) "interpolate_cld, ice: r,q,w,g ",ext_REI(ik), Qice, Wice, Gice, ext_cICE(ik)
          tau_cld_temp(icecld,iw,ik) = 3.*Qice  / (4.*rho_ice*r_ice) * ext_cICE(ik)
          singscat_cld_temp(icecld,iw,ik) = Wice
          asym_cld_temp(icecld,iw,ik) = Gice
        endif

      enddo
    enddo

    ! Call MCICA sub-column generator (Monte Carlo Independent Column Approximation for cloud overlap)

    ! Select cloud overlap approach (1=random, 2=maximum-random, 3=maximum)
    icldovr = 2
    nstep = get_nstep()
    permuteseed = (nstep)

    call mcica_subcol(icldovr, permuteseed, ext_pmid, ext_cFRC, ext_cICE, ext_cLIQ, &
                      tau_cld_temp(icecld,:,:), singscat_cld_temp(icecld,:,:), asym_cld_temp(icecld,:,:), &
                      tau_cld_temp(liqcld,:,:), singscat_cld_temp(liqcld,:,:), asym_cld_temp(liqcld,:,:), &
                      cFRC_mcica, cICE_mcica, cLIQ_mcica, &
                      tau_mcica_ice, ssa_mcica_ice, asym_mcica_ice, &
                      tau_mcica_liq, ssa_mcica_liq, asym_mcica_liq )

    tau_cld_mcica(icecld,:,:) = tau_mcica_ice
    singscat_cld_mcica(icecld,:,:) = ssa_mcica_ice
    asym_cld_mcica(icecld,:,:) = asym_mcica_ice
    tau_cld_mcica(liqcld,:,:) = tau_mcica_liq
    singscat_cld_mcica(liqcld,:,:) = ssa_mcica_liq
    asym_cld_mcica(liqcld,:,:) = asym_mcica_liq

    return

  end subroutine calc_cldopd

!============================================================================

  subroutine calc_aeropd( )

!------------------------------------------------------------------------
!
! Purpose: Calculate the optical depths of aerosols
!          Optical depths stored in calculates the current 'tau_aer(1:ip,wavelength_band,vertical_level)',
!          'singscat_albd(1:ip,wavelength_band,vertical_level)',
!          'asym_fact(1:ip,wavelength_band,vertical_level)' [i.e., layer opacity,
!          single-scattering albedo, and asymmetry factor of aerosol specie 'ip'].
!
!------------------------------------------------------------------------
! NOTES: Incomplete, hook up with CARMA haze aerosols

    implicit none

!------------------------------------------------------------------------
!
! Input Arguments
!
!
! Local Variables
!

!------------------------------------------------------------------------
!
! Start Code
!
   !loop over wavelengths
   !loop oer levels find optical coefficients
   !find optical depths over particles sizes
   ! find single scatter albedo

    return
  end subroutine calc_aeropd



end module calc_opd_mod
