
module calc_opd_mod
! version n68equiv
!
! Idea:  For every grid box, first compare the grey optical depths.
! Which ever has the largest grey optical depth (k*u), use the k-distribution.
! the other gases will be treated using equivalent extinction.
! This method would allow for dynamic selection of major and minor species 
! continuously.

  use shr_kind_mod,     only: r8 => shr_kind_r8
  use rayleigh_data
  use shr_const_mod,    only: SHR_CONST_PI,SHR_CONST_PI, SHR_CONST_G, &
                              SHR_CONST_RGAS, SHR_CONST_AVOGAD, &
                              SHR_CONST_STEBOL, &
                              SHR_CONST_BOLTZ, &
                              SHR_CONST_RHOFW, SHR_CONST_RHOICE, &
                              SHR_CONST_LOSCHMIDT
  use physconst,        only: mwn2, mwco2, mwch4, mwh2o, mwo2, mwh2, mwo3, mwdry, cpair, epsilo
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

  subroutine calc_gasopd(tmid, pmid, pdel, coldens, coldens_dry, qh2o, qco2, qch4, qO2, qO3, qH2, qN2, pathlength, &
                         tau_gas, tau_ray)

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
    real(r8), intent(in), dimension(pver) :: pdel          ! layer thickness [mb]  
    real(r8), intent(in), dimension(pverp) :: coldens      ! Wet Column density profile [molec m-2]
    real(r8), intent(in), dimension(pverp) :: coldens_dry  ! Dry Column density profile [molec m-2]
    real(r8), intent(in), dimension(pverp) :: qh2o         ! mass mixing ratio h2o profile [kg/kg] wet
    real(r8), intent(in), dimension(pverp) :: qco2         ! mass mixing ratio co2 profile [kg/kg] dry
    real(r8), intent(in), dimension(pverp) :: qch4         ! mass mixing ratio ch4 profile [kg/kg] dry
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
    ! layer pressure, temperature place hodler
    real(r8) :: pressure
    real(r8) :: temperature

    ! gas volume mixing ratios
    real(r8) :: h2ovmr, co2vmr, ch4vmr
    real(r8) :: h2vmr, n2vmr
    real(r8) :: o2vmr, o3vmr

    ! indices for interpolation
    integer :: p_ref_index, t_ref_index, t_ref_index_s
    integer :: t_ref_index_h2oh2o, t_ref_index_h2on2
    integer :: t_ref_index_h2h2,t_ref_index_n2h2, t_ref_index_n2n2
    integer :: t_ref_index_co2co2_sw, t_ref_index_co2co2_lw
    integer ::  t_ref_index_co2ch4, t_ref_index_co2h2
    integer :: ik, iw, ig, itu, itl, itc  ! indices
    integer  :: iwbeg, iwend  ! first and last band

    ! absorption coefficient "answers" returned from interpolators
    real(r8) :: ans_kmajor, ans_kgrey_h2o, ans_kgrey_co2, ans_kgrey_ch4, ans
    real(r8), dimension(ngpt_max) :: ans_kmajor_gptvec
    real(r8), dimension(ntot_wavlnrng) :: ans_cia, ans_h2os_avg, ans_h2of_avg
    real(r8), dimension(8, ntot_wavlnrng) :: ans_h2os, ans_h2of

    ! Gas quantities
    real(r8) :: u_col, u_h2o, u_co2, u_ch4, u_h2, u_n2, u_o2, u_o3, u_tot, u_tot_dry
    real(r8), dimension(nspecies) :: ugas, tau_grey
    integer, dimension(1) :: imajor

    ! place holder temperatures for interpolation
    real(r8) :: t_kgas, t_n2n2, t_n2h2, t_h2h2, t_h2os, t_h2oh2o, t_h2on2
    real(r8) :: t_co2ch4, t_co2h2, t_co2co2_sw, t_co2co2_lw

    real(r8) :: wm, wl, wla, r, ns, sp, w

    ! Bps continuum variables
    real(r8) :: arg1, arg2, radfield
    real(r8), dimension(ntot_gpt,pverp) ::  tau_h2oself 
    real(r8), dimension(ntot_gpt,pverp) ::  tau_h2oforeign     

    ! for rayleigh scattering calc, depolarization
    real(r8) :: depolN2, depolCO2, depolH2O      
    ! for rayleigh scattering calc, Allen (1976) coefficients
    real(r8) :: allenN2, allenCO2      
    ! rayleigh scattering cross sections [cm2 molecule-1]
    real(r8) :: sigmaRayl, sigmaRaylCO2, sigmaRaylN2, sigmaRaylH2, sigmaRaylH2O
    real(r8) :: kg_sw_minval  !! minimum value to check sw_abs error

    ! partial pressures
    real(r8) :: ppN2, ppH2, ppCO2, ppCH4, ppH2O

    ! amagats
    real(r8) :: amaN2, amaH2, amaCO2, amaCH4, amaH2O, amaFRGN

    ! CIA optical depths
    real(r8), dimension(ntot_wavlnrng,pverp) ::  tau_cia  ! total
    ! individual
    real(r8), dimension(ntot_wavlnrng,pverp) ::  tau_h2oh2o_cia
    real(r8), dimension(ntot_wavlnrng,pverp) ::  tau_h2on2_cia
    real(r8), dimension(ntot_wavlnrng,pverp) ::  tau_n2n2_cia
    real(r8), dimension(ntot_wavlnrng,pverp) ::  tau_n2h2_cia
    real(r8), dimension(ntot_wavlnrng,pverp) ::  tau_h2h2_cia
    real(r8), dimension(ntot_wavlnrng,pverp) ::  tau_co2co2_lw_cia
    real(r8), dimension(ntot_wavlnrng,pverp) ::  tau_co2co2_sw_cia
    real(r8), dimension(ntot_wavlnrng,pverp) ::  tau_co2h2_cia
    real(r8), dimension(ntot_wavlnrng,pverp) ::  tau_co2ch4_cia

    ! MT_CKD optical depths
    real(r8), dimension(ntot_wavlnrng,pverp) ::  mtckd_self_tau
    real(r8), dimension(ntot_wavlnrng,pverp) ::  mtckd_frgn_tau

 
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
    tau_h2oh2o_cia(:,:) = 0.0
    tau_h2on2_cia(:,:) = 0.0
    tau_n2n2_cia(:,:) = 0.0
    tau_n2h2_cia(:,:) = 0.0
    tau_h2h2_cia(:,:) = 0.0
    tau_co2co2_lw_cia(:,:) = 0.0
    tau_co2co2_sw_cia(:,:) = 0.0
    tau_co2h2_cia(:,:) = 0.0
    tau_co2ch4_cia(:,:) = 0.0

    do ik = 1,pverp
    
!write(*,*) "------- level ------ ", ik, " -----------------------"

      !! Define grid box volume mixing ratios and column densities

      ! qco2, qh2o, qch4 are in mass mixing ratio, multiply by mwdry/mw*** to get volume mixing ratio
      ! factor of 10^-4 converts molecules m-2 to molecules cm-2
      ! qh2o is defined relative to wet atm mass
      ! all others are defined relative to dry atm mass
      ! set all vmr quantities relative to dry mass except where noted
      w = qh2o(ik)/(1.0-qh2o(ik))           ! H2O mass mixing ratio relative to dry air
      h2ovmr = w*mwdry/mwh2o                ! H2O dry volume mixing ratio

      co2vmr = qco2(ik)*mwdry/mwco2         ! CO2 volume mixing ratio dry
      ch4vmr = qch4(ik)*mwdry/mwch4         ! CH4 volume mixing ratio dry
      o2vmr = qo2(ik)*mwdry/mwo2            ! O2 volume mixing ratio dry
      o3vmr = qo3(ik)*mwdry/mwo3            ! O3 volume mixing ratio dry
      h2vmr = qh2(ik)*mwdry/mwh2            ! H2 volume mixing ratio dry
      n2vmr = qn2(ik)*mwdry/mwn2            ! N2 volume mixing ratio dry

! kludge value for experiments
!ch4vmr=1.0e-6
!co2vmr=0.0

      u_h2o = h2ovmr*coldens_dry(ik)/10000.     !   water column amount [ molecules cm-2 ]
      u_co2 = co2vmr*coldens_dry(ik)/10000.     !   co2 column amount [ molecules cm-2 ]
      u_ch4 = ch4vmr*coldens_dry(ik)/10000.     !   ch4 column amount [ molecules cm-2 ]
      u_o2 = o2vmr*coldens_dry(ik)/10000.       !   o2 column amount [ molecules cm-2 ]
      u_o3 = o3vmr*coldens_dry(ik)/10000.       !   o3 column amount [ molecules cm-2 ]
      u_h2 = h2vmr*coldens_dry(ik)/10000.       !   h2 column amount [ molecules cm-2 ]
      u_n2 = n2vmr*coldens_dry(ik)/10000.       !   n2 column amount [ molecules cm-2 ]
      u_tot = u_h2o + u_co2 + u_ch4 + u_o2 + u_o3 + u_h2 + u_n2
      u_tot_dry = u_tot - u_h2o

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
      ugas = (/u_h2o, u_co2, u_ch4/)

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
        if (t_kgas .gt. tgrid(kn2h2_ntemp)) then  ! temperature greater than grid max
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
      !=====  interval 1:   =====!
      !=====  0 - 40 cm-1  
      sp=1
      call bilinear_interpK_grey(k01_grey_data, iH2O, pressure, p_ref_index, t_kgas, t_ref_index, ans_kgrey_h2o) 
      call bilinear_interpK_grey(k01_grey_data, iCO2, pressure, p_ref_index, t_kgas, t_ref_index, ans_kgrey_co2) 
      call bilinear_interpK_grey(k01_grey_data, iCH4, pressure, p_ref_index, t_kgas, t_ref_index, ans_kgrey_ch4)
      tau_grey(iH2O) = ans_kgrey_h2o * ugas(iH2O)
      tau_grey(iCO2) = ans_kgrey_co2 * ugas(iCO2)
      tau_grey(iCH4) = ans_kgrey_ch4 * ugas(iCH4)
      imajor = maxloc(tau_grey)
!write(*,*) ik, sp, gas_name(imajor), tau_grey, ans_kgrey_h2o, ans_kgrey_co2, ans_kgrey_ch4
      ! major gas (k-distribution)
      call bilinear_interpK_8gpt_major_gptvec(k01_major_data, imajor(1), pressure, p_ref_index, t_kgas, t_ref_index, ans_kmajor_gptvec) 
!!write(*,*) "kmajor", ans_kmajor_gptvec
      do ig = 1,ngauss_pts(sp)
        itl = itl + 1
        tau_gas(itl,ik) = ans_kmajor_gptvec(ig)*ugas(imajor(1)) + (tau_grey(iH2O) + tau_grey(iCO2) + tau_grey(iCH4) - tau_grey(imajor(1)))
      enddo

      !=====  interval 2:   =====!
      !=====  40 - 100 cm-1  
      sp=2
      call bilinear_interpK_grey(k02_grey_data, iH2O, pressure, p_ref_index, t_kgas, t_ref_index, ans_kgrey_h2o) 
      call bilinear_interpK_grey(k02_grey_data, iCO2, pressure, p_ref_index, t_kgas, t_ref_index, ans_kgrey_co2) 
      call bilinear_interpK_grey(k02_grey_data, iCH4, pressure, p_ref_index, t_kgas, t_ref_index, ans_kgrey_ch4)
      tau_grey(iH2O) = ans_kgrey_h2o * ugas(iH2O)
      tau_grey(iCO2) = ans_kgrey_co2 * ugas(iCO2)
      tau_grey(iCH4) = ans_kgrey_ch4 * ugas(iCH4)
      imajor = maxloc(tau_grey)
!write(*,*) ik, sp, gas_name(imajor), tau_grey, ans_kgrey_h2o, ans_kgrey_co2, ans_kgrey_ch4
      ! major gas (k-distribution)      
      call bilinear_interpK_8gpt_major_gptvec(k02_major_data, imajor(1), pressure, p_ref_index, t_kgas, t_ref_index, ans_kmajor_gptvec) 
!!write(*,*) "kmajor", ans_kmajor_gptvec
      do ig = 1,ngauss_pts(sp)      
        itl = itl + 1         
        tau_gas(itl,ik) = ans_kmajor_gptvec(ig)*ugas(imajor(1)) + (tau_grey(iH2O) + tau_grey(iCO2) + tau_grey(iCH4) - tau_grey(imajor(1)))
      enddo   

      !=====  interval 3:   =====!
      !=====  100 - 160 cm-1  
      sp=3
      call bilinear_interpK_grey(k03_grey_data, iH2O, pressure, p_ref_index, t_kgas, t_ref_index, ans_kgrey_h2o) 
      call bilinear_interpK_grey(k03_grey_data, iCO2, pressure, p_ref_index, t_kgas, t_ref_index, ans_kgrey_co2) 
      call bilinear_interpK_grey(k03_grey_data, iCH4, pressure, p_ref_index, t_kgas, t_ref_index, ans_kgrey_ch4)
      tau_grey(iH2O) = ans_kgrey_h2o * ugas(iH2O)
      tau_grey(iCO2) = ans_kgrey_co2 * ugas(iCO2)
      tau_grey(iCH4) = ans_kgrey_ch4 * ugas(iCH4)
      imajor = maxloc(tau_grey)
!write(*,*) ik, sp, gas_name(imajor), tau_grey, ans_kgrey_h2o, ans_kgrey_co2, ans_kgrey_ch4
      ! major gas (k-distribution)
      call bilinear_interpK_8gpt_major_gptvec(k03_major_data, imajor(1), pressure, p_ref_index, t_kgas, t_ref_index, ans_kmajor_gptvec) 
!!write(*,*) "kmajor", ans_kmajor_gptvec
      do ig = 1,ngauss_pts(sp)      
        itl = itl + 1         
        tau_gas(itl,ik) = ans_kmajor_gptvec(ig)*ugas(imajor(1)) + (tau_grey(iH2O) + tau_grey(iCO2) + tau_grey(iCH4) - tau_grey(imajor(1)))
      enddo   

      !=====  interval 4:   =====!
      !=====  160 - 220 cm-1  
      sp=4
      call bilinear_interpK_grey(k04_grey_data, iH2O, pressure, p_ref_index, t_kgas, t_ref_index, ans_kgrey_h2o) 
      call bilinear_interpK_grey(k04_grey_data, iCO2, pressure, p_ref_index, t_kgas, t_ref_index, ans_kgrey_co2) 
      call bilinear_interpK_grey(k04_grey_data, iCH4, pressure, p_ref_index, t_kgas, t_ref_index, ans_kgrey_ch4)
      tau_grey(iH2O) = ans_kgrey_h2o * ugas(iH2O)
      tau_grey(iCO2) = ans_kgrey_co2 * ugas(iCO2)
      tau_grey(iCH4) = ans_kgrey_ch4 * ugas(iCH4)
      imajor = maxloc(tau_grey)
!write(*,*) ik, sp, gas_name(imajor), tau_grey, ans_kgrey_h2o, ans_kgrey_co2, ans_kgrey_ch4
      ! major gas (k-distribution)
      call bilinear_interpK_8gpt_major_gptvec(k04_major_data, imajor(1), pressure, p_ref_index, t_kgas, t_ref_index, ans_kmajor_gptvec) 
!!write(*,*) "kmajor", ans_kmajor_gptvec
      do ig = 1,ngauss_pts(sp)      
        itl = itl + 1         
        tau_gas(itl,ik) = ans_kmajor_gptvec(ig)*ugas(imajor(1)) + (tau_grey(iH2O) + tau_grey(iCO2) + tau_grey(iCH4) - tau_grey(imajor(1)))
      enddo   

      !=====  interval 5:   =====!
      !=====  220 - 280 cm-1  
      sp=5
      call bilinear_interpK_grey(k05_grey_data, iH2O, pressure, p_ref_index, t_kgas, t_ref_index, ans_kgrey_h2o) 
      call bilinear_interpK_grey(k05_grey_data, iCO2, pressure, p_ref_index, t_kgas, t_ref_index, ans_kgrey_co2) 
      call bilinear_interpK_grey(k05_grey_data, iCH4, pressure, p_ref_index, t_kgas, t_ref_index, ans_kgrey_ch4)
      tau_grey(iH2O) = ans_kgrey_h2o * ugas(iH2O)
      tau_grey(iCO2) = ans_kgrey_co2 * ugas(iCO2)
      tau_grey(iCH4) = ans_kgrey_ch4 * ugas(iCH4)
      imajor = maxloc(tau_grey)
!write(*,*) ik, sp, gas_name(imajor), tau_grey, ans_kgrey_h2o, ans_kgrey_co2, ans_kgrey_ch4
      ! major gas (k-distribution)
      call bilinear_interpK_8gpt_major_gptvec(k05_major_data, imajor(1), pressure, p_ref_index, t_kgas, t_ref_index, ans_kmajor_gptvec) 
!!write(*,*) "kmajor", ans_kmajor_gptvec
      do ig = 1,ngauss_pts(sp)      
        itl = itl + 1          
        tau_gas(itl,ik) = ans_kmajor_gptvec(ig)*ugas(imajor(1)) + (tau_grey(iH2O) + tau_grey(iCO2) + tau_grey(iCH4) - tau_grey(imajor(1)))
      enddo   

      !=====  interval 6:   =====!
      !=====  280 - 330 cm-1  
      sp=6
      call bilinear_interpK_grey(k06_grey_data, iH2O, pressure, p_ref_index, t_kgas, t_ref_index, ans_kgrey_h2o) 
      call bilinear_interpK_grey(k06_grey_data, iCO2, pressure, p_ref_index, t_kgas, t_ref_index, ans_kgrey_co2) 
      call bilinear_interpK_grey(k06_grey_data, iCH4, pressure, p_ref_index, t_kgas, t_ref_index, ans_kgrey_ch4)
      tau_grey(iH2O) = ans_kgrey_h2o * ugas(iH2O)
      tau_grey(iCO2) = ans_kgrey_co2 * ugas(iCO2)
      tau_grey(iCH4) = ans_kgrey_ch4 * ugas(iCH4)
      imajor = maxloc(tau_grey)
!write(*,*) ik, sp, gas_name(imajor), tau_grey, ans_kgrey_h2o, ans_kgrey_co2, ans_kgrey_ch4
      ! major gas (k-distribution)
      call bilinear_interpK_8gpt_major_gptvec(k06_major_data, imajor(1), pressure, p_ref_index, t_kgas, t_ref_index, ans_kmajor_gptvec) 
!!write(*,*) "kmajor", ans_kmajor_gptvec
      do ig = 1,ngauss_pts(sp)      
        itl = itl + 1        
        tau_gas(itl,ik) = ans_kmajor_gptvec(ig)*ugas(imajor(1)) + (tau_grey(iH2O) + tau_grey(iCO2) + tau_grey(iCH4) - tau_grey(imajor(1)))
      enddo   

      !=====  interval 7:   =====!
      !=====  330 - 380 cm-1  
      sp=7
      call bilinear_interpK_grey(k07_grey_data, iH2O, pressure, p_ref_index, t_kgas, t_ref_index, ans_kgrey_h2o) 
      call bilinear_interpK_grey(k07_grey_data, iCO2, pressure, p_ref_index, t_kgas, t_ref_index, ans_kgrey_co2) 
      call bilinear_interpK_grey(k07_grey_data, iCH4, pressure, p_ref_index, t_kgas, t_ref_index, ans_kgrey_ch4)
      tau_grey(iH2O) = ans_kgrey_h2o * ugas(iH2O)
      tau_grey(iCO2) = ans_kgrey_co2 * ugas(iCO2)
      tau_grey(iCH4) = ans_kgrey_ch4 * ugas(iCH4)
      imajor = maxloc(tau_grey)
!write(*,*) ik, sp, gas_name(imajor), tau_grey, ans_kgrey_h2o, ans_kgrey_co2, ans_kgrey_ch4
      ! major gas (k-distribution)
      call bilinear_interpK_8gpt_major_gptvec(k07_major_data, imajor(1), pressure, p_ref_index, t_kgas, t_ref_index, ans_kmajor_gptvec) 
!!write(*,*) "kmajor", ans_kmajor_gptvec
      do ig = 1,ngauss_pts(sp)      
        itl = itl + 1         
        tau_gas(itl,ik) = ans_kmajor_gptvec(ig)*ugas(imajor(1)) + (tau_grey(iH2O) + tau_grey(iCO2) + tau_grey(iCH4) - tau_grey(imajor(1)))
      enddo   

      !=====  interval 8:   =====!
      !=====  380 - 440 cm-1  
      sp=8
      call bilinear_interpK_grey(k08_grey_data, iH2O, pressure, p_ref_index, t_kgas, t_ref_index, ans_kgrey_h2o) 
      call bilinear_interpK_grey(k08_grey_data, iCO2, pressure, p_ref_index, t_kgas, t_ref_index, ans_kgrey_co2) 
      call bilinear_interpK_grey(k08_grey_data, iCH4, pressure, p_ref_index, t_kgas, t_ref_index, ans_kgrey_ch4)
      tau_grey(iH2O) = ans_kgrey_h2o * ugas(iH2O)
      tau_grey(iCO2) = ans_kgrey_co2 * ugas(iCO2)
      tau_grey(iCH4) = ans_kgrey_ch4 * ugas(iCH4)
      imajor = maxloc(tau_grey)
!write(*,*) ik, sp, gas_name(imajor), tau_grey, ans_kgrey_h2o, ans_kgrey_co2, ans_kgrey_ch4
      ! major gas (k-distribution)
      call bilinear_interpK_8gpt_major_gptvec(k08_major_data, imajor(1), pressure, p_ref_index, t_kgas, t_ref_index, ans_kmajor_gptvec) 
!!write(*,*) "kmajor", ans_kmajor_gptvec
      do ig = 1,ngauss_pts(sp)      
        itl = itl + 1         
        tau_gas(itl,ik) = ans_kmajor_gptvec(ig)*ugas(imajor(1)) + (tau_grey(iH2O) + tau_grey(iCO2) + tau_grey(iCH4) - tau_grey(imajor(1)))
      enddo   

      !=====  interval 9:   =====!
      !=====  440 - 495 cm-1  
      sp=9
      call bilinear_interpK_grey(k09_grey_data, iH2O, pressure, p_ref_index, t_kgas, t_ref_index, ans_kgrey_h2o) 
      call bilinear_interpK_grey(k09_grey_data, iCO2, pressure, p_ref_index, t_kgas, t_ref_index, ans_kgrey_co2) 
      call bilinear_interpK_grey(k09_grey_data, iCH4, pressure, p_ref_index, t_kgas, t_ref_index, ans_kgrey_ch4)
      tau_grey(iH2O) = ans_kgrey_h2o * ugas(iH2O)
      tau_grey(iCO2) = ans_kgrey_co2 * ugas(iCO2)
      tau_grey(iCH4) = ans_kgrey_ch4 * ugas(iCH4)
      imajor = maxloc(tau_grey)
!write(*,*) ik, sp, gas_name(imajor), tau_grey, ans_kgrey_h2o, ans_kgrey_co2, ans_kgrey_ch4
      ! major gas (k-distribution)
      call bilinear_interpK_8gpt_major_gptvec(k09_major_data, imajor(1), pressure, p_ref_index, t_kgas, t_ref_index, ans_kmajor_gptvec) 
!!write(*,*) "kmajor", ans_kmajor_gptvec
      do ig = 1,ngauss_pts(sp)      
        itl = itl + 1         
        tau_gas(itl,ik) = ans_kmajor_gptvec(ig)*ugas(imajor(1)) + (tau_grey(iH2O) + tau_grey(iCO2) + tau_grey(iCH4) - tau_grey(imajor(1)))
      enddo   

      !=====  interval 10:  =====!
      !=====  495 - 545 cm-1  
      sp=10
      call bilinear_interpK_grey(k10_grey_data, iH2O, pressure, p_ref_index, t_kgas, t_ref_index, ans_kgrey_h2o) 
      call bilinear_interpK_grey(k10_grey_data, iCO2, pressure, p_ref_index, t_kgas, t_ref_index, ans_kgrey_co2) 
      call bilinear_interpK_grey(k10_grey_data, iCH4, pressure, p_ref_index, t_kgas, t_ref_index, ans_kgrey_ch4)
      tau_grey(iH2O) = ans_kgrey_h2o * ugas(iH2O)
      tau_grey(iCO2) = ans_kgrey_co2 * ugas(iCO2)
      tau_grey(iCH4) = ans_kgrey_ch4 * ugas(iCH4)
      imajor = maxloc(tau_grey)
!write(*,*) ik, sp, gas_name(imajor), tau_grey, ans_kgrey_h2o, ans_kgrey_co2, ans_kgrey_ch4
      ! major gas (k-distribution)
      call bilinear_interpK_8gpt_major_gptvec(k10_major_data, imajor(1), pressure, p_ref_index, t_kgas, t_ref_index, ans_kmajor_gptvec) 
!!write(*,*) "kmajor", ans_kmajor_gptvec
      do ig = 1,ngauss_pts(sp)      
        itl = itl + 1         
        tau_gas(itl,ik) = ans_kmajor_gptvec(ig)*ugas(imajor(1)) + (tau_grey(iH2O) + tau_grey(iCO2) + tau_grey(iCH4) - tau_grey(imajor(1)))
      enddo   

      !=====  interval 11:  =====!
      !=====  545 - 617 cm-1  
      sp=11
      call bilinear_interpK_grey(k11_grey_data, iH2O, pressure, p_ref_index, t_kgas, t_ref_index, ans_kgrey_h2o) 
      call bilinear_interpK_grey(k11_grey_data, iCO2, pressure, p_ref_index, t_kgas, t_ref_index, ans_kgrey_co2) 
      call bilinear_interpK_grey(k11_grey_data, iCH4, pressure, p_ref_index, t_kgas, t_ref_index, ans_kgrey_ch4)
      tau_grey(iH2O) = ans_kgrey_h2o * ugas(iH2O)
      tau_grey(iCO2) = ans_kgrey_co2 * ugas(iCO2)
      tau_grey(iCH4) = ans_kgrey_ch4 * ugas(iCH4)
      imajor = maxloc(tau_grey)
!write(*,*) ik, sp, gas_name(imajor), tau_grey, ans_kgrey_h2o, ans_kgrey_co2, ans_kgrey_ch4
      ! major gas (k-distribution)
      call bilinear_interpK_8gpt_major_gptvec(k11_major_data, imajor(1), pressure, p_ref_index, t_kgas, t_ref_index, ans_kmajor_gptvec) 
!!write(*,*) "kmajor", ans_kmajor_gptvec
      do ig = 1,ngauss_pts(sp)      
        itl = itl + 1         
        tau_gas(itl,ik) = ans_kmajor_gptvec(ig)*ugas(imajor(1)) + (tau_grey(iH2O) + tau_grey(iCO2) + tau_grey(iCH4) - tau_grey(imajor(1)))
      enddo   

      !=====  interval 12:  =====
      !=====  617 - 667 cm-1 
      sp=12
      call bilinear_interpK_grey(k12_grey_data, iH2O, pressure, p_ref_index, t_kgas, t_ref_index, ans_kgrey_h2o) 
      call bilinear_interpK_grey(k12_grey_data, iCO2, pressure, p_ref_index, t_kgas, t_ref_index, ans_kgrey_co2) 
      call bilinear_interpK_grey(k12_grey_data, iCH4, pressure, p_ref_index, t_kgas, t_ref_index, ans_kgrey_ch4)
      tau_grey(iH2O) = ans_kgrey_h2o * ugas(iH2O)
      tau_grey(iCO2) = ans_kgrey_co2 * ugas(iCO2)
      tau_grey(iCH4) = ans_kgrey_ch4 * ugas(iCH4)
      imajor = maxloc(tau_grey)
!write(*,*) ik, sp, gas_name(imajor), tau_grey, ans_kgrey_h2o, ans_kgrey_co2, ans_kgrey_ch4
      ! major gas (k-distribution)
      call bilinear_interpK_8gpt_major_gptvec(k12_major_data, imajor(1), pressure, p_ref_index, t_kgas, t_ref_index, ans_kmajor_gptvec) 
!!write(*,*) "kmajor", ans_kmajor_gptvec
      do ig = 1,ngauss_pts(sp)      
        itl = itl + 1         
        tau_gas(itl,ik) = ans_kmajor_gptvec(ig)*ugas(imajor(1)) + (tau_grey(iH2O) + tau_grey(iCO2) + tau_grey(iCH4) - tau_grey(imajor(1)))
      enddo   

      !=====  interval 13:  =====!
      !=====  667 - 720 cm-1 
      sp=13
      call bilinear_interpK_grey(k13_grey_data, iH2O, pressure, p_ref_index, t_kgas, t_ref_index, ans_kgrey_h2o) 
      call bilinear_interpK_grey(k13_grey_data, iCO2, pressure, p_ref_index, t_kgas, t_ref_index, ans_kgrey_co2) 
      call bilinear_interpK_grey(k13_grey_data, iCH4, pressure, p_ref_index, t_kgas, t_ref_index, ans_kgrey_ch4)
      tau_grey(iH2O) = ans_kgrey_h2o * ugas(iH2O)
      tau_grey(iCO2) = ans_kgrey_co2 * ugas(iCO2)
      tau_grey(iCH4) = ans_kgrey_ch4 * ugas(iCH4)
      imajor = maxloc(tau_grey)
!write(*,*) ik, sp, gas_name(imajor), tau_grey, ans_kgrey_h2o, ans_kgrey_co2, ans_kgrey_ch4
      ! major gas (k-distribution)
      call bilinear_interpK_8gpt_major_gptvec(k13_major_data, imajor(1), pressure, p_ref_index, t_kgas, t_ref_index, ans_kmajor_gptvec) 
!!write(*,*) "kmajor", ans_kmajor_gptvec
      do ig = 1,ngauss_pts(sp)      
        itl = itl + 1         
        tau_gas(itl,ik) = ans_kmajor_gptvec(ig)*ugas(imajor(1)) + (tau_grey(iH2O) + tau_grey(iCO2) + tau_grey(iCH4) - tau_grey(imajor(1)))
      enddo   

      !=====  interval 14:  =====!
      !=====  720 - 800 cm-1 
      sp=14
      call bilinear_interpK_grey(k14_grey_data, iH2O, pressure, p_ref_index, t_kgas, t_ref_index, ans_kgrey_h2o) 
      call bilinear_interpK_grey(k14_grey_data, iCO2, pressure, p_ref_index, t_kgas, t_ref_index, ans_kgrey_co2) 
      call bilinear_interpK_grey(k14_grey_data, iCH4, pressure, p_ref_index, t_kgas, t_ref_index, ans_kgrey_ch4)
      tau_grey(iH2O) = ans_kgrey_h2o * ugas(iH2O)
      tau_grey(iCO2) = ans_kgrey_co2 * ugas(iCO2)
      tau_grey(iCH4) = ans_kgrey_ch4 * ugas(iCH4)
      imajor = maxloc(tau_grey)
!write(*,*) ik, sp, gas_name(imajor), tau_grey, ans_kgrey_h2o, ans_kgrey_co2, ans_kgrey_ch4
      ! major gas (k-distribution)
      call bilinear_interpK_8gpt_major_gptvec(k14_major_data, imajor(1), pressure, p_ref_index, t_kgas, t_ref_index, ans_kmajor_gptvec) 
!!write(*,*) "kmajor", ans_kmajor_gptvec
      do ig = 1,ngauss_pts(sp)      
        itl = itl + 1         
        tau_gas(itl,ik) = ans_kmajor_gptvec(ig)*ugas(imajor(1)) + (tau_grey(iH2O) + tau_grey(iCO2) + tau_grey(iCH4) - tau_grey(imajor(1)))
      enddo   

      !=====  interval 15:  =====!
      !=====  800 - 875 cm-1 
      sp=15
      call bilinear_interpK_grey(k15_grey_data, iH2O, pressure, p_ref_index, t_kgas, t_ref_index, ans_kgrey_h2o) 
      call bilinear_interpK_grey(k15_grey_data, iCO2, pressure, p_ref_index, t_kgas, t_ref_index, ans_kgrey_co2) 
      call bilinear_interpK_grey(k15_grey_data, iCH4, pressure, p_ref_index, t_kgas, t_ref_index, ans_kgrey_ch4)
      tau_grey(iH2O) = ans_kgrey_h2o * ugas(iH2O)
      tau_grey(iCO2) = ans_kgrey_co2 * ugas(iCO2)
      tau_grey(iCH4) = ans_kgrey_ch4 * ugas(iCH4)
      imajor = maxloc(tau_grey)
!write(*,*) ik, sp, gas_name(imajor), tau_grey, ans_kgrey_h2o, ans_kgrey_co2, ans_kgrey_ch4
      ! major gas (k-distribution)
      call bilinear_interpK_8gpt_major_gptvec(k15_major_data, imajor(1), pressure, p_ref_index, t_kgas, t_ref_index, ans_kmajor_gptvec) 
!!write(*,*) "kmajor", ans_kmajor_gptvec
      do ig = 1,ngauss_pts(sp)      
        itl = itl + 1         
        tau_gas(itl,ik) = ans_kmajor_gptvec(ig)*ugas(imajor(1)) + (tau_grey(iH2O) + tau_grey(iCO2) + tau_grey(iCH4) - tau_grey(imajor(1)))
      enddo   

      !=====  interval 16:  =====!
      !=====  875 - 940 cm-1 
      sp=16
      call bilinear_interpK_grey(k16_grey_data, iH2O, pressure, p_ref_index, t_kgas, t_ref_index, ans_kgrey_h2o) 
      call bilinear_interpK_grey(k16_grey_data, iCO2, pressure, p_ref_index, t_kgas, t_ref_index, ans_kgrey_co2) 
      call bilinear_interpK_grey(k16_grey_data, iCH4, pressure, p_ref_index, t_kgas, t_ref_index, ans_kgrey_ch4)
      tau_grey(iH2O) = ans_kgrey_h2o * ugas(iH2O)
      tau_grey(iCO2) = ans_kgrey_co2 * ugas(iCO2)
      tau_grey(iCH4) = ans_kgrey_ch4 * ugas(iCH4)
      imajor = maxloc(tau_grey)
!write(*,*) ik, sp, gas_name(imajor), tau_grey, ans_kgrey_h2o, ans_kgrey_co2, ans_kgrey_ch4
      ! major gas (k-distribution)
      call bilinear_interpK_8gpt_major_gptvec(k16_major_data, imajor(1), pressure, p_ref_index, t_kgas, t_ref_index, ans_kmajor_gptvec) 
!!write(*,*) "kmajor", ans_kmajor_gptvec
      do ig = 1,ngauss_pts(sp)      
        itl = itl + 1         
        tau_gas(itl,ik) = ans_kmajor_gptvec(ig)*ugas(imajor(1)) + (tau_grey(iH2O) + tau_grey(iCO2) + tau_grey(iCH4) - tau_grey(imajor(1)))
      enddo   

      !=====  interval 17:  =====!
      !=====  940 - 1000 cm-1 
      sp=17
      call bilinear_interpK_grey(k17_grey_data, iH2O, pressure, p_ref_index, t_kgas, t_ref_index, ans_kgrey_h2o) 
      call bilinear_interpK_grey(k17_grey_data, iCO2, pressure, p_ref_index, t_kgas, t_ref_index, ans_kgrey_co2) 
      call bilinear_interpK_grey(k17_grey_data, iCH4, pressure, p_ref_index, t_kgas, t_ref_index, ans_kgrey_ch4)
      tau_grey(iH2O) = ans_kgrey_h2o * ugas(iH2O)
      tau_grey(iCO2) = ans_kgrey_co2 * ugas(iCO2)
      tau_grey(iCH4) = ans_kgrey_ch4 * ugas(iCH4)
      imajor = maxloc(tau_grey)
!write(*,*) ik, sp, gas_name(imajor), tau_grey, ans_kgrey_h2o, ans_kgrey_co2, ans_kgrey_ch4
      ! major gas (k-distribution)
      call bilinear_interpK_8gpt_major_gptvec(k17_major_data, imajor(1), pressure, p_ref_index, t_kgas, t_ref_index, ans_kmajor_gptvec) 
!!write(*,*) "kmajor", ans_kmajor_gptvec
      do ig = 1,ngauss_pts(sp)      
        itl = itl + 1         
        tau_gas(itl,ik) = ans_kmajor_gptvec(ig)*ugas(imajor(1)) + (tau_grey(iH2O) + tau_grey(iCO2) + tau_grey(iCH4) - tau_grey(imajor(1)))
      enddo   

      !=====  interval 18:  =====!
      !=====  1000 - 1065 cm-1 
      sp=18
      call bilinear_interpK_grey(k18_grey_data, iH2O, pressure, p_ref_index, t_kgas, t_ref_index, ans_kgrey_h2o) 
      call bilinear_interpK_grey(k18_grey_data, iCO2, pressure, p_ref_index, t_kgas, t_ref_index, ans_kgrey_co2) 
      call bilinear_interpK_grey(k18_grey_data, iCH4, pressure, p_ref_index, t_kgas, t_ref_index, ans_kgrey_ch4)
      tau_grey(iH2O) = ans_kgrey_h2o * ugas(iH2O)
      tau_grey(iCO2) = ans_kgrey_co2 * ugas(iCO2)
      tau_grey(iCH4) = ans_kgrey_ch4 * ugas(iCH4)
      imajor = maxloc(tau_grey)
!write(*,*) ik, sp, gas_name(imajor), tau_grey, ans_kgrey_h2o, ans_kgrey_co2, ans_kgrey_ch4
      ! major gas (k-distribution)
      call bilinear_interpK_8gpt_major_gptvec(k18_major_data, imajor(1), pressure, p_ref_index, t_kgas, t_ref_index, ans_kmajor_gptvec) 
!!write(*,*) "kmajor", ans_kmajor_gptvec
      do ig = 1,ngauss_pts(sp)      
        itl = itl + 1         
        tau_gas(itl,ik) = ans_kmajor_gptvec(ig)*ugas(imajor(1)) + (tau_grey(iH2O) + tau_grey(iCO2) + tau_grey(iCH4) - tau_grey(imajor(1)))
      enddo   

      !=====  interval 19:  =====!
      !=====  1065 - 1108 cm-1 
      sp=19
      call bilinear_interpK_grey(k19_grey_data, iH2O, pressure, p_ref_index, t_kgas, t_ref_index, ans_kgrey_h2o) 
      call bilinear_interpK_grey(k19_grey_data, iCO2, pressure, p_ref_index, t_kgas, t_ref_index, ans_kgrey_co2) 
      call bilinear_interpK_grey(k19_grey_data, iCH4, pressure, p_ref_index, t_kgas, t_ref_index, ans_kgrey_ch4)
      tau_grey(iH2O) = ans_kgrey_h2o * ugas(iH2O)
      tau_grey(iCO2) = ans_kgrey_co2 * ugas(iCO2)
      tau_grey(iCH4) = ans_kgrey_ch4 * ugas(iCH4)
      imajor = maxloc(tau_grey)
!write(*,*) ik, sp, gas_name(imajor), tau_grey, ans_kgrey_h2o, ans_kgrey_co2, ans_kgrey_ch4
      ! major gas (k-distribution)
      call bilinear_interpK_8gpt_major_gptvec(k19_major_data, imajor(1), pressure, p_ref_index, t_kgas, t_ref_index, ans_kmajor_gptvec) 
!!write(*,*) "kmajor", ans_kmajor_gptvec
      do ig = 1,ngauss_pts(sp)      
        itl = itl + 1         
        tau_gas(itl,ik) = ans_kmajor_gptvec(ig)*ugas(imajor(1)) + (tau_grey(iH2O) + tau_grey(iCO2) + tau_grey(iCH4) - tau_grey(imajor(1)))
      enddo   

      !=====  interval 20:  =====!
      !=====  1108 - 1200 cm-1 
      sp=20
      call bilinear_interpK_grey(k20_grey_data, iH2O, pressure, p_ref_index, t_kgas, t_ref_index, ans_kgrey_h2o) 
      call bilinear_interpK_grey(k20_grey_data, iCO2, pressure, p_ref_index, t_kgas, t_ref_index, ans_kgrey_co2) 
      call bilinear_interpK_grey(k20_grey_data, iCH4, pressure, p_ref_index, t_kgas, t_ref_index, ans_kgrey_ch4)
      tau_grey(iH2O) = ans_kgrey_h2o * ugas(iH2O)
      tau_grey(iCO2) = ans_kgrey_co2 * ugas(iCO2)
      tau_grey(iCH4) = ans_kgrey_ch4 * ugas(iCH4)
      imajor = maxloc(tau_grey)
!write(*,*) ik, sp, gas_name(imajor), tau_grey, ans_kgrey_h2o, ans_kgrey_co2, ans_kgrey_ch4
      ! major gas (k-distribution)
      call bilinear_interpK_8gpt_major_gptvec(k20_major_data, imajor(1), pressure, p_ref_index, t_kgas, t_ref_index, ans_kmajor_gptvec) 
!!write(*,*) "kmajor", ans_kmajor_gptvec
      do ig = 1,ngauss_pts(sp)      
        itl = itl + 1         
        tau_gas(itl,ik) = ans_kmajor_gptvec(ig)*ugas(imajor(1)) + (tau_grey(iH2O) + tau_grey(iCO2) + tau_grey(iCH4) - tau_grey(imajor(1)))
      enddo   

      !=====  interval 21:  =====!
      !=====  1200 - 1275 cm-1 
      sp=21
      call bilinear_interpK_grey(k21_grey_data, iH2O, pressure, p_ref_index, t_kgas, t_ref_index, ans_kgrey_h2o) 
      call bilinear_interpK_grey(k21_grey_data, iCO2, pressure, p_ref_index, t_kgas, t_ref_index, ans_kgrey_co2) 
      call bilinear_interpK_grey(k21_grey_data, iCH4, pressure, p_ref_index, t_kgas, t_ref_index, ans_kgrey_ch4)
      tau_grey(iH2O) = ans_kgrey_h2o * ugas(iH2O)
      tau_grey(iCO2) = ans_kgrey_co2 * ugas(iCO2)
      tau_grey(iCH4) = ans_kgrey_ch4 * ugas(iCH4)
      imajor = maxloc(tau_grey)
!write(*,*) ik, sp, gas_name(imajor), tau_grey, ans_kgrey_h2o, ans_kgrey_co2, ans_kgrey_ch4
      ! major gas (k-distribution)
      call bilinear_interpK_8gpt_major_gptvec(k21_major_data, imajor(1), pressure, p_ref_index, t_kgas, t_ref_index, ans_kmajor_gptvec) 
!!write(*,*) "kmajor", ans_kmajor_gptvec
      do ig = 1,ngauss_pts(sp)      
        itl = itl + 1         
        tau_gas(itl,ik) = ans_kmajor_gptvec(ig)*ugas(imajor(1)) + (tau_grey(iH2O) + tau_grey(iCO2) + tau_grey(iCH4) - tau_grey(imajor(1)))
      enddo   

      !=====  interval 22:  =====!
      !=====  1275 - 1350 cm-1 
      sp=22
      call bilinear_interpK_grey(k22_grey_data, iH2O, pressure, p_ref_index, t_kgas, t_ref_index, ans_kgrey_h2o) 
      call bilinear_interpK_grey(k22_grey_data, iCO2, pressure, p_ref_index, t_kgas, t_ref_index, ans_kgrey_co2) 
      call bilinear_interpK_grey(k22_grey_data, iCH4, pressure, p_ref_index, t_kgas, t_ref_index, ans_kgrey_ch4)
      tau_grey(iH2O) = ans_kgrey_h2o * ugas(iH2O)
      tau_grey(iCO2) = ans_kgrey_co2 * ugas(iCO2)
      tau_grey(iCH4) = ans_kgrey_ch4 * ugas(iCH4)
      imajor = maxloc(tau_grey)
!write(*,*) ik, sp, gas_name(imajor), tau_grey, ans_kgrey_h2o, ans_kgrey_co2, ans_kgrey_ch4
      ! major gas (k-distribution)
      call bilinear_interpK_8gpt_major_gptvec(k22_major_data, imajor(1), pressure, p_ref_index, t_kgas, t_ref_index, ans_kmajor_gptvec) 
!!write(*,*) "kmajor", ans_kmajor_gptvec
      do ig = 1,ngauss_pts(sp)      
        itl = itl + 1         
        tau_gas(itl,ik) = ans_kmajor_gptvec(ig)*ugas(imajor(1)) + (tau_grey(iH2O) + tau_grey(iCO2) + tau_grey(iCH4) - tau_grey(imajor(1)))
      enddo   

      !=====  interval 23:  =====!
      !=====  1350 - 1450 cm-1 
      sp=23
      call bilinear_interpK_grey(k23_grey_data, iH2O, pressure, p_ref_index, t_kgas, t_ref_index, ans_kgrey_h2o) 
      call bilinear_interpK_grey(k23_grey_data, iCO2, pressure, p_ref_index, t_kgas, t_ref_index, ans_kgrey_co2) 
      call bilinear_interpK_grey(k23_grey_data, iCH4, pressure, p_ref_index, t_kgas, t_ref_index, ans_kgrey_ch4)
      tau_grey(iH2O) = ans_kgrey_h2o * ugas(iH2O)
      tau_grey(iCO2) = ans_kgrey_co2 * ugas(iCO2)
      tau_grey(iCH4) = ans_kgrey_ch4 * ugas(iCH4)
      imajor = maxloc(tau_grey)
!write(*,*) ik, sp, gas_name(imajor), tau_grey, ans_kgrey_h2o, ans_kgrey_co2, ans_kgrey_ch4
      ! major gas (k-distribution)
      call bilinear_interpK_8gpt_major_gptvec(k23_major_data, imajor(1), pressure, p_ref_index, t_kgas, t_ref_index, ans_kmajor_gptvec) 
!!write(*,*) "kmajor", ans_kmajor_gptvec
      do ig = 1,ngauss_pts(sp)      
        itl = itl + 1         
        tau_gas(itl,ik) = ans_kmajor_gptvec(ig)*ugas(imajor(1)) + (tau_grey(iH2O) + tau_grey(iCO2) + tau_grey(iCH4) - tau_grey(imajor(1)))
      enddo   

      !=====  interval 24:  =====!
      !=====  1450 - 1550 cm-1 
      sp=24
      call bilinear_interpK_grey(k24_grey_data, iH2O, pressure, p_ref_index, t_kgas, t_ref_index, ans_kgrey_h2o) 
      call bilinear_interpK_grey(k24_grey_data, iCO2, pressure, p_ref_index, t_kgas, t_ref_index, ans_kgrey_co2) 
      call bilinear_interpK_grey(k24_grey_data, iCH4, pressure, p_ref_index, t_kgas, t_ref_index, ans_kgrey_ch4)
      tau_grey(iH2O) = ans_kgrey_h2o * ugas(iH2O)
      tau_grey(iCO2) = ans_kgrey_co2 * ugas(iCO2)
      tau_grey(iCH4) = ans_kgrey_ch4 * ugas(iCH4)
       imajor = maxloc(tau_grey)
!write(*,*) ik, sp, gas_name(imajor), tau_grey, ans_kgrey_h2o, ans_kgrey_co2, ans_kgrey_ch4
      ! major gas (k-distribution)
      call bilinear_interpK_8gpt_major_gptvec(k24_major_data, imajor(1), pressure, p_ref_index, t_kgas, t_ref_index, ans_kmajor_gptvec) 
!!write(*,*) "kmajor", ans_kmajor_gptvec
      do ig = 1,ngauss_pts(sp)      
        itl = itl + 1         
        tau_gas(itl,ik) = ans_kmajor_gptvec(ig)*ugas(imajor(1)) + (tau_grey(iH2O) + tau_grey(iCO2) + tau_grey(iCH4) - tau_grey(imajor(1)))
      enddo   

      !=====  interval 25:  =====!
      !=====  1550 - 1650 cm-1 
      sp=25
      call bilinear_interpK_grey(k25_grey_data, iH2O, pressure, p_ref_index, t_kgas, t_ref_index, ans_kgrey_h2o) 
      call bilinear_interpK_grey(k25_grey_data, iCO2, pressure, p_ref_index, t_kgas, t_ref_index, ans_kgrey_co2) 
      call bilinear_interpK_grey(k25_grey_data, iCH4, pressure, p_ref_index, t_kgas, t_ref_index, ans_kgrey_ch4)
      tau_grey(iH2O) = ans_kgrey_h2o * ugas(iH2O)
      tau_grey(iCO2) = ans_kgrey_co2 * ugas(iCO2)
      tau_grey(iCH4) = ans_kgrey_ch4 * ugas(iCH4)
      imajor = maxloc(tau_grey)
!write(*,*) ik, sp, gas_name(imajor), tau_grey, ans_kgrey_h2o, ans_kgrey_co2, ans_kgrey_ch4
      ! major gas (k-distribution)
      ! kmajor = k25_major_data(imajor(1),:,:,:)
      call bilinear_interpK_8gpt_major_gptvec(k25_major_data, imajor(1), pressure, p_ref_index, t_kgas, t_ref_index, ans_kmajor_gptvec) 
!!write(*,*) "kmajor", ans_kmajor_gptvec
      do ig = 1,ngauss_pts(sp)      
        itl = itl + 1         
        tau_gas(itl,ik) = ans_kmajor_gptvec(ig)*ugas(imajor(1)) + (tau_grey(iH2O) + tau_grey(iCO2) + tau_grey(iCH4) - tau_grey(imajor(1)))
      enddo   

      !=====  interval 26:  =====!
      !=====  1650 - 1750 cm-1 
      sp=26
      call bilinear_interpK_grey(k26_grey_data, iH2O, pressure, p_ref_index, t_kgas, t_ref_index, ans_kgrey_h2o) 
      call bilinear_interpK_grey(k26_grey_data, iCO2, pressure, p_ref_index, t_kgas, t_ref_index, ans_kgrey_co2) 
      call bilinear_interpK_grey(k26_grey_data, iCH4, pressure, p_ref_index, t_kgas, t_ref_index, ans_kgrey_ch4)
      tau_grey(iH2O) = ans_kgrey_h2o * ugas(iH2O)
      tau_grey(iCO2) = ans_kgrey_co2 * ugas(iCO2)
      tau_grey(iCH4) = ans_kgrey_ch4 * ugas(iCH4)
      imajor = maxloc(tau_grey)
!write(*,*) ik, sp, gas_name(imajor), tau_grey, ans_kgrey_h2o, ans_kgrey_co2, ans_kgrey_ch4
      ! major gas (k-distribution)
      call bilinear_interpK_8gpt_major_gptvec(k26_major_data, imajor(1), pressure, p_ref_index, t_kgas, t_ref_index, ans_kmajor_gptvec) 
!!write(*,*) "kmajor", ans_kmajor_gptvec
      do ig = 1,ngauss_pts(sp)      
        itl = itl + 1         
        tau_gas(itl,ik) = ans_kmajor_gptvec(ig)*ugas(imajor(1)) + (tau_grey(iH2O) + tau_grey(iCO2) + tau_grey(iCH4) - tau_grey(imajor(1)))
      enddo   

      !=====  interval 27:  =====!
      !=====  1750 - 1850 cm-1 
      sp=27
      call bilinear_interpK_grey(k27_grey_data, iH2O, pressure, p_ref_index, t_kgas, t_ref_index, ans_kgrey_h2o) 
      call bilinear_interpK_grey(k27_grey_data, iCO2, pressure, p_ref_index, t_kgas, t_ref_index, ans_kgrey_co2) 
      call bilinear_interpK_grey(k27_grey_data, iCH4, pressure, p_ref_index, t_kgas, t_ref_index, ans_kgrey_ch4)
      tau_grey(iH2O) = ans_kgrey_h2o * ugas(iH2O)
      tau_grey(iCO2) = ans_kgrey_co2 * ugas(iCO2)
      tau_grey(iCH4) = ans_kgrey_ch4 * ugas(iCH4)
      imajor = maxloc(tau_grey)
!write(*,*) ik, sp, gas_name(imajor), tau_grey, ans_kgrey_h2o, ans_kgrey_co2, ans_kgrey_ch4
      ! major gas (k-distribution)
      call bilinear_interpK_8gpt_major_gptvec(k27_major_data, imajor(1), pressure, p_ref_index, t_kgas, t_ref_index, ans_kmajor_gptvec) 
!!write(*,*) "kmajor", ans_kmajor_gptvec
      do ig = 1,ngauss_pts(sp)      
        itl = itl + 1         
        tau_gas(itl,ik) = ans_kmajor_gptvec(ig)*ugas(imajor(1)) + (tau_grey(iH2O) + tau_grey(iCO2) + tau_grey(iCH4) - tau_grey(imajor(1)))
      enddo   

      !=====  interval 28:  =====!
      !=====  1850 - 1950 cm-1 
      sp=28
      call bilinear_interpK_grey(k28_grey_data, iH2O, pressure, p_ref_index, t_kgas, t_ref_index, ans_kgrey_h2o) 
      call bilinear_interpK_grey(k28_grey_data, iCO2, pressure, p_ref_index, t_kgas, t_ref_index, ans_kgrey_co2) 
      call bilinear_interpK_grey(k28_grey_data, iCH4, pressure, p_ref_index, t_kgas, t_ref_index, ans_kgrey_ch4)
      tau_grey(iH2O) = ans_kgrey_h2o * ugas(iH2O)
      tau_grey(iCO2) = ans_kgrey_co2 * ugas(iCO2)
      tau_grey(iCH4) = ans_kgrey_ch4 * ugas(iCH4)
      imajor = maxloc(tau_grey)
!write(*,*) ik, sp, gas_name(imajor), tau_grey, ans_kgrey_h2o, ans_kgrey_co2, ans_kgrey_ch4
      ! major gas (k-distribution)
      call bilinear_interpK_8gpt_major_gptvec(k28_major_data, imajor(1), pressure, p_ref_index, t_kgas, t_ref_index, ans_kmajor_gptvec) 
!!write(*,*) "kmajor", ans_kmajor_gptvec
      do ig = 1,ngauss_pts(sp)      
        itl = itl + 1         
        tau_gas(itl,ik) = ans_kmajor_gptvec(ig)*ugas(imajor(1)) + (tau_grey(iH2O) + tau_grey(iCO2) + tau_grey(iCH4) - tau_grey(imajor(1)))
      enddo   

      !=====  interval 29:  =====!
      !=====  1950 - 2050 cm-1 
      sp=29
      call bilinear_interpK_grey(k29_grey_data, iH2O, pressure, p_ref_index, t_kgas, t_ref_index, ans_kgrey_h2o) 
      call bilinear_interpK_grey(k29_grey_data, iCO2, pressure, p_ref_index, t_kgas, t_ref_index, ans_kgrey_co2) 
      call bilinear_interpK_grey(k29_grey_data, iCH4, pressure, p_ref_index, t_kgas, t_ref_index, ans_kgrey_ch4)
      tau_grey(iH2O) = ans_kgrey_h2o * ugas(iH2O)
      tau_grey(iCO2) = ans_kgrey_co2 * ugas(iCO2)
      tau_grey(iCH4) = ans_kgrey_ch4 * ugas(iCH4)
      imajor = maxloc(tau_grey)
!write(*,*) ik, sp, gas_name(imajor), tau_grey, ans_kgrey_h2o, ans_kgrey_co2, ans_kgrey_ch4
      ! major gas (k-distribution)
      call bilinear_interpK_8gpt_major_gptvec(k29_major_data, imajor(1), pressure, p_ref_index, t_kgas, t_ref_index, ans_kmajor_gptvec) 
!!write(*,*) "kmajor", ans_kmajor_gptvec
      do ig = 1,ngauss_pts(sp)      
        itl = itl + 1         
        tau_gas(itl,ik) = ans_kmajor_gptvec(ig)*ugas(imajor(1)) + (tau_grey(iH2O) + tau_grey(iCO2) + tau_grey(iCH4) - tau_grey(imajor(1)))
      enddo   

      !=====  interval 30:  =====!
      !=====  2050 - 2200 cm-1 
      sp=30
      call bilinear_interpK_grey(k30_grey_data, iH2O, pressure, p_ref_index, t_kgas, t_ref_index, ans_kgrey_h2o) 
      call bilinear_interpK_grey(k30_grey_data, iCO2, pressure, p_ref_index, t_kgas, t_ref_index, ans_kgrey_co2) 
      call bilinear_interpK_grey(k30_grey_data, iCH4, pressure, p_ref_index, t_kgas, t_ref_index, ans_kgrey_ch4)
      tau_grey(iH2O) = ans_kgrey_h2o * ugas(iH2O)
      tau_grey(iCO2) = ans_kgrey_co2 * ugas(iCO2)
      tau_grey(iCH4) = ans_kgrey_ch4 * ugas(iCH4)
      imajor = maxloc(tau_grey)
!write(*,*) ik, sp, gas_name(imajor), tau_grey, ans_kgrey_h2o, ans_kgrey_co2, ans_kgrey_ch4
      ! major gas (k-distribution)
      call bilinear_interpK_8gpt_major_gptvec(k30_major_data, imajor(1), pressure, p_ref_index, t_kgas, t_ref_index, ans_kmajor_gptvec) 
!!write(*,*) "kmajor", ans_kmajor_gptvec
      do ig = 1,ngauss_pts(sp)      
        itl = itl + 1         
        tau_gas(itl,ik) = ans_kmajor_gptvec(ig)*ugas(imajor(1)) + (tau_grey(iH2O) + tau_grey(iCO2) + tau_grey(iCH4) - tau_grey(imajor(1)))
      enddo   

      !=====  interval 31:  =====!
      !=====  2200 - 2397 cm-1 
      sp=31
      call bilinear_interpK_grey(k31_grey_data, iH2O, pressure, p_ref_index, t_kgas, t_ref_index, ans_kgrey_h2o) 
      call bilinear_interpK_grey(k31_grey_data, iCO2, pressure, p_ref_index, t_kgas, t_ref_index, ans_kgrey_co2) 
      call bilinear_interpK_grey(k31_grey_data, iCH4, pressure, p_ref_index, t_kgas, t_ref_index, ans_kgrey_ch4)
      tau_grey(iH2O) = ans_kgrey_h2o * ugas(iH2O)
      tau_grey(iCO2) = ans_kgrey_co2 * ugas(iCO2)
      tau_grey(iCH4) = ans_kgrey_ch4 * ugas(iCH4)
      imajor = maxloc(tau_grey)
!write(*,*) ik, sp, gas_name(imajor), tau_grey, ans_kgrey_h2o, ans_kgrey_co2, ans_kgrey_ch4
      ! major gas (k-distribution)
      call bilinear_interpK_8gpt_major_gptvec(k31_major_data, imajor(1), pressure, p_ref_index, t_kgas, t_ref_index, ans_kmajor_gptvec) 
!!write(*,*) "kmajor", ans_kmajor_gptvec
      do ig = 1,ngauss_pts(sp)      
        itl = itl + 1         
        tau_gas(itl,ik) = ans_kmajor_gptvec(ig)*ugas(imajor(1)) + (tau_grey(iH2O) + tau_grey(iCO2) + tau_grey(iCH4) - tau_grey(imajor(1)))
      enddo   

      !=====  interval 32:  =====
      !=====  2397 - 2494 cm-1 
      sp=32
      call bilinear_interpK_grey(k32_grey_data, iH2O, pressure, p_ref_index, t_kgas, t_ref_index, ans_kgrey_h2o) 
      call bilinear_interpK_grey(k32_grey_data, iCO2, pressure, p_ref_index, t_kgas, t_ref_index, ans_kgrey_co2) 
      call bilinear_interpK_grey(k32_grey_data, iCH4, pressure, p_ref_index, t_kgas, t_ref_index, ans_kgrey_ch4)
      tau_grey(iH2O) = ans_kgrey_h2o * ugas(iH2O)
      tau_grey(iCO2) = ans_kgrey_co2 * ugas(iCO2)
      tau_grey(iCH4) = ans_kgrey_ch4 * ugas(iCH4)
      imajor = maxloc(tau_grey)
!write(*,*) ik, sp, gas_name(imajor), tau_grey, ans_kgrey_h2o, ans_kgrey_co2, ans_kgrey_ch4
      ! major gas (k-distribution)
      call bilinear_interpK_8gpt_major_gptvec(k32_major_data, imajor(1), pressure, p_ref_index, t_kgas, t_ref_index, ans_kmajor_gptvec) 
!!write(*,*) "kmajor", ans_kmajor_gptvec
      do ig = 1,ngauss_pts(sp)      
        itl = itl + 1         
        tau_gas(itl,ik) = ans_kmajor_gptvec(ig)*ugas(imajor(1)) + (tau_grey(iH2O) + tau_grey(iCO2) + tau_grey(iCH4) - tau_grey(imajor(1)))
      enddo   

      !=====  interval 33:  =====!
      !=====  2494 - 2796 cm-1 
      sp=33
      call bilinear_interpK_grey(k33_grey_data, iH2O, pressure, p_ref_index, t_kgas, t_ref_index, ans_kgrey_h2o) 
      call bilinear_interpK_grey(k33_grey_data, iCO2, pressure, p_ref_index, t_kgas, t_ref_index, ans_kgrey_co2) 
      call bilinear_interpK_grey(k33_grey_data, iCH4, pressure, p_ref_index, t_kgas, t_ref_index, ans_kgrey_ch4)
      tau_grey(iH2O) = ans_kgrey_h2o * ugas(iH2O)
      tau_grey(iCO2) = ans_kgrey_co2 * ugas(iCO2)
      tau_grey(iCH4) = ans_kgrey_ch4 * ugas(iCH4)
      imajor = maxloc(tau_grey)
!write(*,*) ik, sp, gas_name(imajor), tau_grey, ans_kgrey_h2o, ans_kgrey_co2, ans_kgrey_ch4
      ! major gas (k-distribution)
      call bilinear_interpK_8gpt_major_gptvec(k33_major_data, imajor(1), pressure, p_ref_index, t_kgas, t_ref_index, ans_kmajor_gptvec) 
!!write(*,*) "kmajor", ans_kmajor_gptvec
      do ig = 1,ngauss_pts(sp)      
        itl = itl + 1         
        tau_gas(itl,ik) = ans_kmajor_gptvec(ig)*ugas(imajor(1)) + (tau_grey(iH2O) + tau_grey(iCO2) + tau_grey(iCH4) - tau_grey(imajor(1)))
      enddo   

      !=====  interval 34:  =====!
      !=====  2796 - 3087 cm-1 
      sp=34
      call bilinear_interpK_grey(k34_grey_data, iH2O, pressure, p_ref_index, t_kgas, t_ref_index, ans_kgrey_h2o) 
      call bilinear_interpK_grey(k34_grey_data, iCO2, pressure, p_ref_index, t_kgas, t_ref_index, ans_kgrey_co2) 
      call bilinear_interpK_grey(k34_grey_data, iCH4, pressure, p_ref_index, t_kgas, t_ref_index, ans_kgrey_ch4)
      tau_grey(iH2O) = ans_kgrey_h2o * ugas(iH2O)
      tau_grey(iCO2) = ans_kgrey_co2 * ugas(iCO2)
      tau_grey(iCH4) = ans_kgrey_ch4 * ugas(iCH4)
      imajor = maxloc(tau_grey)
!write(*,*) ik, sp, gas_name(imajor), tau_grey, ans_kgrey_h2o, ans_kgrey_co2, ans_kgrey_ch4
      ! major gas (k-distribution)
      call bilinear_interpK_8gpt_major_gptvec(k34_major_data, imajor(1), pressure, p_ref_index, t_kgas, t_ref_index, ans_kmajor_gptvec) 
!!write(*,*) "kmajor", ans_kmajor_gptvec
      do ig = 1,ngauss_pts(sp)      
        itl = itl + 1         
        tau_gas(itl,ik) = ans_kmajor_gptvec(ig)*ugas(imajor(1)) + (tau_grey(iH2O) + tau_grey(iCO2) + tau_grey(iCH4) - tau_grey(imajor(1)))
      enddo   

      !=====  interval 35:  =====!
      !=====  3087 - 3425 cm-1 
      sp=35
      call bilinear_interpK_grey(k35_grey_data, iH2O, pressure, p_ref_index, t_kgas, t_ref_index, ans_kgrey_h2o) 
      call bilinear_interpK_grey(k35_grey_data, iCO2, pressure, p_ref_index, t_kgas, t_ref_index, ans_kgrey_co2) 
      call bilinear_interpK_grey(k35_grey_data, iCH4, pressure, p_ref_index, t_kgas, t_ref_index, ans_kgrey_ch4)
      tau_grey(iH2O) = ans_kgrey_h2o * ugas(iH2O)
      tau_grey(iCO2) = ans_kgrey_co2 * ugas(iCO2)
      tau_grey(iCH4) = ans_kgrey_ch4 * ugas(iCH4)
      imajor = maxloc(tau_grey)
!write(*,*) ik, sp, gas_name(imajor), tau_grey, ans_kgrey_h2o, ans_kgrey_co2, ans_kgrey_ch4
      ! major gas (k-distribution)
      call bilinear_interpK_8gpt_major_gptvec(k35_major_data, imajor(1), pressure, p_ref_index, t_kgas, t_ref_index, ans_kmajor_gptvec) 
!!write(*,*) "kmajor", ans_kmajor_gptvec
      do ig = 1,ngauss_pts(sp)      
        itl = itl + 1         
        tau_gas(itl,ik) = ans_kmajor_gptvec(ig)*ugas(imajor(1)) + (tau_grey(iH2O) + tau_grey(iCO2) + tau_grey(iCH4) - tau_grey(imajor(1)))
      enddo   
      
      !=====  interval 36:  =====!
      !=====  3425 - 3760 cm-1 
      sp=36
      call bilinear_interpK_grey(k36_grey_data, iH2O, pressure, p_ref_index, t_kgas, t_ref_index, ans_kgrey_h2o) 
      call bilinear_interpK_grey(k36_grey_data, iCO2, pressure, p_ref_index, t_kgas, t_ref_index, ans_kgrey_co2) 
      call bilinear_interpK_grey(k36_grey_data, iCH4, pressure, p_ref_index, t_kgas, t_ref_index, ans_kgrey_ch4)
      tau_grey(iH2O) = ans_kgrey_h2o * ugas(iH2O)
      tau_grey(iCO2) = ans_kgrey_co2 * ugas(iCO2)
      tau_grey(iCH4) = ans_kgrey_ch4 * ugas(iCH4)
      imajor = maxloc(tau_grey)
!write(*,*) ik, sp, gas_name(imajor), tau_grey, ans_kgrey_h2o, ans_kgrey_co2, ans_kgrey_ch4
      ! major gas (k-distribution)
      ! kmajor = k36_major_data(imajor(1),:,:,:)
      call bilinear_interpK_8gpt_major_gptvec(k36_major_data, imajor(1), pressure, p_ref_index, t_kgas, t_ref_index, ans_kmajor_gptvec) 
!!write(*,*) "kmajor", ans_kmajor_gptvec
      do ig = 1,ngauss_pts(sp)      
        itl = itl + 1         
        tau_gas(itl,ik) = ans_kmajor_gptvec(ig)*ugas(imajor(1)) + (tau_grey(iH2O) + tau_grey(iCO2) + tau_grey(iCH4) - tau_grey(imajor(1)))
      enddo   

      !=====  interval 37:  =====!
      !=====  3760 - 4030 cm-1 
      sp=37
      call bilinear_interpK_grey(k37_grey_data, iH2O, pressure, p_ref_index, t_kgas, t_ref_index, ans_kgrey_h2o) 
      call bilinear_interpK_grey(k37_grey_data, iCO2, pressure, p_ref_index, t_kgas, t_ref_index, ans_kgrey_co2) 
      call bilinear_interpK_grey(k37_grey_data, iCH4, pressure, p_ref_index, t_kgas, t_ref_index, ans_kgrey_ch4)
      tau_grey(iH2O) = ans_kgrey_h2o * ugas(iH2O)
      tau_grey(iCO2) = ans_kgrey_co2 * ugas(iCO2)
      tau_grey(iCH4) = ans_kgrey_ch4 * ugas(iCH4)
      imajor = maxloc(tau_grey)
!write(*,*) ik, sp, gas_name(imajor), tau_grey, ans_kgrey_h2o, ans_kgrey_co2, ans_kgrey_ch4
      ! major gas (k-distribution)
      call bilinear_interpK_8gpt_major_gptvec(k37_major_data, imajor(1), pressure, p_ref_index, t_kgas, t_ref_index, ans_kmajor_gptvec) 
!!write(*,*) "kmajor", ans_kmajor_gptvec
      do ig = 1,ngauss_pts(sp)      
        itl = itl + 1         
        tau_gas(itl,ik) = ans_kmajor_gptvec(ig)*ugas(imajor(1)) + (tau_grey(iH2O) + tau_grey(iCO2) + tau_grey(iCH4) - tau_grey(imajor(1)))
      enddo   

      !=====  interval 38:  =====!
      !=====  4030 - 4540 cm-1 
      sp=38
      call bilinear_interpK_grey(k38_grey_data, iH2O, pressure, p_ref_index, t_kgas, t_ref_index, ans_kgrey_h2o) 
      call bilinear_interpK_grey(k38_grey_data, iCO2, pressure, p_ref_index, t_kgas, t_ref_index, ans_kgrey_co2) 
      call bilinear_interpK_grey(k38_grey_data, iCH4, pressure, p_ref_index, t_kgas, t_ref_index, ans_kgrey_ch4)
      tau_grey(iH2O) = ans_kgrey_h2o * ugas(iH2O)
      tau_grey(iCO2) = ans_kgrey_co2 * ugas(iCO2)
      tau_grey(iCH4) = ans_kgrey_ch4 * ugas(iCH4)
      imajor = maxloc(tau_grey)
!write(*,*) ik, sp, gas_name(imajor), tau_grey, ans_kgrey_h2o, ans_kgrey_co2, ans_kgrey_ch4
      ! major gas (k-distribution)
      call bilinear_interpK_8gpt_major_gptvec(k38_major_data, imajor(1), pressure, p_ref_index, t_kgas, t_ref_index, ans_kmajor_gptvec) 
!!write(*,*) "kmajor", ans_kmajor_gptvec
      do ig = 1,ngauss_pts(sp)      
        itl = itl + 1         
        tau_gas(itl,ik) = ans_kmajor_gptvec(ig)*ugas(imajor(1)) + (tau_grey(iH2O) + tau_grey(iCO2) + tau_grey(iCH4) - tau_grey(imajor(1)))
      enddo   

      !=====  interval 39:  =====!
      !=====  4540 - 4950 cm-1 
      sp=39
      call bilinear_interpK_grey(k39_grey_data, iH2O, pressure, p_ref_index, t_kgas, t_ref_index, ans_kgrey_h2o) 
      call bilinear_interpK_grey(k39_grey_data, iCO2, pressure, p_ref_index, t_kgas, t_ref_index, ans_kgrey_co2) 
      call bilinear_interpK_grey(k39_grey_data, iCH4, pressure, p_ref_index, t_kgas, t_ref_index, ans_kgrey_ch4)
      tau_grey(iH2O) = ans_kgrey_h2o * ugas(iH2O)
      tau_grey(iCO2) = ans_kgrey_co2 * ugas(iCO2)
      tau_grey(iCH4) = ans_kgrey_ch4 * ugas(iCH4)
      imajor = maxloc(tau_grey)
!write(*,*) ik, sp, gas_name(imajor), tau_grey, ans_kgrey_h2o, ans_kgrey_co2, ans_kgrey_ch4
      ! major gas (k-distribution)
      call bilinear_interpK_8gpt_major_gptvec(k39_major_data, imajor(1), pressure, p_ref_index, t_kgas, t_ref_index, ans_kmajor_gptvec) 
!!write(*,*) "kmajor", ans_kmajor_gptvec
      do ig = 1,ngauss_pts(sp)      
        itl = itl + 1         
        tau_gas(itl,ik) = ans_kmajor_gptvec(ig)*ugas(imajor(1)) + (tau_grey(iH2O) + tau_grey(iCO2) + tau_grey(iCH4) - tau_grey(imajor(1)))
      enddo   

      !=====  interval 40:  =====!
      !=====  4950 - 5370 cm-1 
      sp=40
      call bilinear_interpK_grey(k40_grey_data, iH2O, pressure, p_ref_index, t_kgas, t_ref_index, ans_kgrey_h2o) 
      call bilinear_interpK_grey(k40_grey_data, iCO2, pressure, p_ref_index, t_kgas, t_ref_index, ans_kgrey_co2) 
      call bilinear_interpK_grey(k40_grey_data, iCH4, pressure, p_ref_index, t_kgas, t_ref_index, ans_kgrey_ch4)
      tau_grey(iH2O) = ans_kgrey_h2o * ugas(iH2O)
      tau_grey(iCO2) = ans_kgrey_co2 * ugas(iCO2)
      tau_grey(iCH4) = ans_kgrey_ch4 * ugas(iCH4)
      imajor = maxloc(tau_grey)
!write(*,*) ik, sp, gas_name(imajor), tau_grey, ans_kgrey_h2o, ans_kgrey_co2, ans_kgrey_ch4
      ! major gas (k-distribution)
      call bilinear_interpK_8gpt_major_gptvec(k40_major_data, imajor(1), pressure, p_ref_index, t_kgas, t_ref_index, ans_kmajor_gptvec) 
!!write(*,*) "kmajor", ans_kmajor_gptvec
      do ig = 1,ngauss_pts(sp)      
        itl = itl + 1         
        tau_gas(itl,ik) = ans_kmajor_gptvec(ig)*ugas(imajor(1)) + (tau_grey(iH2O) + tau_grey(iCO2) + tau_grey(iCH4) - tau_grey(imajor(1)))
      enddo   

      !=====  interval 41:  =====!
      !=====  5370 - 5925 cm-1 
      sp=41
      call bilinear_interpK_grey(k41_grey_data, iH2O, pressure, p_ref_index, t_kgas, t_ref_index, ans_kgrey_h2o) 
      call bilinear_interpK_grey(k41_grey_data, iCO2, pressure, p_ref_index, t_kgas, t_ref_index, ans_kgrey_co2) 
      call bilinear_interpK_grey(k41_grey_data, iCH4, pressure, p_ref_index, t_kgas, t_ref_index, ans_kgrey_ch4)
      tau_grey(iH2O) = ans_kgrey_h2o * ugas(iH2O)
      tau_grey(iCO2) = ans_kgrey_co2 * ugas(iCO2)
      tau_grey(iCH4) = ans_kgrey_ch4 * ugas(iCH4)
      imajor = maxloc(tau_grey)
!write(*,*) ik, sp, gas_name(imajor), tau_grey, ans_kgrey_h2o, ans_kgrey_co2, ans_kgrey_ch4
      ! major gas (k-distribution)
      call bilinear_interpK_8gpt_major_gptvec(k41_major_data, imajor(1), pressure, p_ref_index, t_kgas, t_ref_index, ans_kmajor_gptvec) 
!!write(*,*) "kmajor", ans_kmajor_gptvec
      do ig = 1,ngauss_pts(sp)      
        itl = itl + 1         
        tau_gas(itl,ik) = ans_kmajor_gptvec(ig)*ugas(imajor(1)) + (tau_grey(iH2O) + tau_grey(iCO2) + tau_grey(iCH4) - tau_grey(imajor(1)))
      enddo   

      !=====  interval 42:  =====!
      !=====  5925 - 6390 cm-1 
      sp=42
      call bilinear_interpK_grey(k42_grey_data, iH2O, pressure, p_ref_index, t_kgas, t_ref_index, ans_kgrey_h2o) 
      call bilinear_interpK_grey(k42_grey_data, iCO2, pressure, p_ref_index, t_kgas, t_ref_index, ans_kgrey_co2) 
      call bilinear_interpK_grey(k42_grey_data, iCH4, pressure, p_ref_index, t_kgas, t_ref_index, ans_kgrey_ch4)
      tau_grey(iH2O) = ans_kgrey_h2o * ugas(iH2O)
      tau_grey(iCO2) = ans_kgrey_co2 * ugas(iCO2)
      tau_grey(iCH4) = ans_kgrey_ch4 * ugas(iCH4)
      imajor = maxloc(tau_grey)
!write(*,*) ik, sp, gas_name(imajor), tau_grey, ans_kgrey_h2o, ans_kgrey_co2, ans_kgrey_ch4
      ! major gas (k-distribution)
      call bilinear_interpK_8gpt_major_gptvec(k42_major_data, imajor(1), pressure, p_ref_index, t_kgas, t_ref_index, ans_kmajor_gptvec) 
!!write(*,*) "kmajor", ans_kmajor_gptvec
      do ig = 1,ngauss_pts(sp)      
        itl = itl + 1         
        tau_gas(itl,ik) = ans_kmajor_gptvec(ig)*ugas(imajor(1)) + (tau_grey(iH2O) + tau_grey(iCO2) + tau_grey(iCH4) - tau_grey(imajor(1)))
      enddo   

      !=====  interval 43:  =====!
      !=====  6390 - 6990 cm-1 
      sp=43
      call bilinear_interpK_grey(k43_grey_data, iH2O, pressure, p_ref_index, t_kgas, t_ref_index, ans_kgrey_h2o) 
      call bilinear_interpK_grey(k43_grey_data, iCO2, pressure, p_ref_index, t_kgas, t_ref_index, ans_kgrey_co2) 
      call bilinear_interpK_grey(k43_grey_data, iCH4, pressure, p_ref_index, t_kgas, t_ref_index, ans_kgrey_ch4)
      tau_grey(iH2O) = ans_kgrey_h2o * ugas(iH2O)
      tau_grey(iCO2) = ans_kgrey_co2 * ugas(iCO2)
      tau_grey(iCH4) = ans_kgrey_ch4 * ugas(iCH4)
      imajor = maxloc(tau_grey)
!write(*,*) ik, sp, gas_name(imajor), tau_grey, ans_kgrey_h2o, ans_kgrey_co2, ans_kgrey_ch4
      ! major gas (k-distribution)
      call bilinear_interpK_8gpt_major_gptvec(k43_major_data, imajor(1), pressure, p_ref_index, t_kgas, t_ref_index, ans_kmajor_gptvec) 
!!write(*,*) "kmajor", ans_kmajor_gptvec
      do ig = 1,ngauss_pts(sp)      
        itl = itl + 1         
        tau_gas(itl,ik) = ans_kmajor_gptvec(ig)*ugas(imajor(1)) + (tau_grey(iH2O) + tau_grey(iCO2) + tau_grey(iCH4) - tau_grey(imajor(1)))
      enddo   

      !=====  interval 44:  =====!
      !=====  6990 - 7650 cm-1 
      sp=44
      call bilinear_interpK_grey(k44_grey_data, iH2O, pressure, p_ref_index, t_kgas, t_ref_index, ans_kgrey_h2o) 
      call bilinear_interpK_grey(k44_grey_data, iCO2, pressure, p_ref_index, t_kgas, t_ref_index, ans_kgrey_co2) 
      call bilinear_interpK_grey(k44_grey_data, iCH4, pressure, p_ref_index, t_kgas, t_ref_index, ans_kgrey_ch4)
      tau_grey(iH2O) = ans_kgrey_h2o * ugas(iH2O)
      tau_grey(iCO2) = ans_kgrey_co2 * ugas(iCO2)
      tau_grey(iCH4) = ans_kgrey_ch4 * ugas(iCH4)
      imajor = maxloc(tau_grey)
!write(*,*) ik, sp, gas_name(imajor), tau_grey, ans_kgrey_h2o, ans_kgrey_co2, ans_kgrey_ch4
      ! major gas (k-distribution)
      call bilinear_interpK_8gpt_major_gptvec(k44_major_data, imajor(1), pressure, p_ref_index, t_kgas, t_ref_index, ans_kmajor_gptvec) 
!!write(*,*) "kmajor", ans_kmajor_gptvec
      do ig = 1,ngauss_pts(sp)      
        itl = itl + 1         
        tau_gas(itl,ik) = ans_kmajor_gptvec(ig)*ugas(imajor(1)) + (tau_grey(iH2O) + tau_grey(iCO2) + tau_grey(iCH4) - tau_grey(imajor(1)))
      enddo   

      !=====  interval 45:  =====!
      !=====  7650 - 8315 cm-1 
      sp=45
      call bilinear_interpK_grey(k45_grey_data, iH2O, pressure, p_ref_index, t_kgas, t_ref_index, ans_kgrey_h2o) 
      call bilinear_interpK_grey(k45_grey_data, iCO2, pressure, p_ref_index, t_kgas, t_ref_index, ans_kgrey_co2) 
      call bilinear_interpK_grey(k45_grey_data, iCH4, pressure, p_ref_index, t_kgas, t_ref_index, ans_kgrey_ch4)
      tau_grey(iH2O) = ans_kgrey_h2o * ugas(iH2O)
      tau_grey(iCO2) = ans_kgrey_co2 * ugas(iCO2)
      tau_grey(iCH4) = ans_kgrey_ch4 * ugas(iCH4)
      imajor = maxloc(tau_grey)
!write(*,*) ik, sp, gas_name(imajor), tau_grey, ans_kgrey_h2o, ans_kgrey_co2, ans_kgrey_ch4
      ! major gas (k-distribution)
      call bilinear_interpK_8gpt_major_gptvec(k45_major_data, imajor(1), pressure, p_ref_index, t_kgas, t_ref_index, ans_kmajor_gptvec) 
!!write(*,*) "kmajor", ans_kmajor_gptvec
      do ig = 1,ngauss_pts(sp)      
        itl = itl + 1         
        tau_gas(itl,ik) = ans_kmajor_gptvec(ig)*ugas(imajor(1)) + (tau_grey(iH2O) + tau_grey(iCO2) + tau_grey(iCH4) - tau_grey(imajor(1)))
      enddo   

      !=====  interval 46:  =====!
      !=====  8315 - 8850 cm-1 
      sp=46
      call bilinear_interpK_grey(k46_grey_data, iH2O, pressure, p_ref_index, t_kgas, t_ref_index, ans_kgrey_h2o) 
      call bilinear_interpK_grey(k46_grey_data, iCO2, pressure, p_ref_index, t_kgas, t_ref_index, ans_kgrey_co2) 
      call bilinear_interpK_grey(k46_grey_data, iCH4, pressure, p_ref_index, t_kgas, t_ref_index, ans_kgrey_ch4)
      tau_grey(iH2O) = ans_kgrey_h2o * ugas(iH2O)
      tau_grey(iCO2) = ans_kgrey_co2 * ugas(iCO2)
      tau_grey(iCH4) = ans_kgrey_ch4 * ugas(iCH4)
      imajor = maxloc(tau_grey)
!write(*,*) ik, sp, gas_name(imajor), tau_grey, ans_kgrey_h2o, ans_kgrey_co2, ans_kgrey_ch4
      ! major gas (k-distribution)
      call bilinear_interpK_8gpt_major_gptvec(k46_major_data, imajor(1), pressure, p_ref_index, t_kgas, t_ref_index, ans_kmajor_gptvec) 
!!write(*,*) "kmajor", ans_kmajor_gptvec
      do ig = 1,ngauss_pts(sp)      
        itl = itl + 1         
        tau_gas(itl,ik) = ans_kmajor_gptvec(ig)*ugas(imajor(1)) + (tau_grey(iH2O) + tau_grey(iCO2) + tau_grey(iCH4) - tau_grey(imajor(1)))
      enddo   

      !=====  interval 47:  =====!
      !=====  8850 - 9350 cm-1 
      sp=47
      call bilinear_interpK_grey(k47_grey_data, iH2O, pressure, p_ref_index, t_kgas, t_ref_index, ans_kgrey_h2o) 
      call bilinear_interpK_grey(k47_grey_data, iCO2, pressure, p_ref_index, t_kgas, t_ref_index, ans_kgrey_co2) 
      call bilinear_interpK_grey(k47_grey_data, iCH4, pressure, p_ref_index, t_kgas, t_ref_index, ans_kgrey_ch4)
      tau_grey(iH2O) = ans_kgrey_h2o * ugas(iH2O)
      tau_grey(iCO2) = ans_kgrey_co2 * ugas(iCO2)
      tau_grey(iCH4) = ans_kgrey_ch4 * ugas(iCH4)
      imajor = maxloc(tau_grey)
!write(*,*) ik, sp, gas_name(imajor), tau_grey, ans_kgrey_h2o, ans_kgrey_co2, ans_kgrey_ch4
      ! major gas (k-distribution)
      call bilinear_interpK_8gpt_major_gptvec(k47_major_data, imajor(1), pressure, p_ref_index, t_kgas, t_ref_index, ans_kmajor_gptvec) 
!!write(*,*) "kmajor", ans_kmajor_gptvec
      do ig = 1,ngauss_pts(sp)      
        itl = itl + 1         
        tau_gas(itl,ik) = ans_kmajor_gptvec(ig)*ugas(imajor(1)) + (tau_grey(iH2O) + tau_grey(iCO2) + tau_grey(iCH4) - tau_grey(imajor(1)))
      enddo   

      !=====  interval 48:  =====!
      !=====  9350 - 9650 cm-1 
      sp=48
      call bilinear_interpK_grey(k48_grey_data, iH2O, pressure, p_ref_index, t_kgas, t_ref_index, ans_kgrey_h2o) 
      call bilinear_interpK_grey(k48_grey_data, iCO2, pressure, p_ref_index, t_kgas, t_ref_index, ans_kgrey_co2) 
      call bilinear_interpK_grey(k48_grey_data, iCH4, pressure, p_ref_index, t_kgas, t_ref_index, ans_kgrey_ch4)
      tau_grey(iH2O) = ans_kgrey_h2o * ugas(iH2O)
      tau_grey(iCO2) = ans_kgrey_co2 * ugas(iCO2)
      tau_grey(iCH4) = ans_kgrey_ch4 * ugas(iCH4)
      imajor = maxloc(tau_grey)
!write(*,*) ik, sp, gas_name(imajor), tau_grey, ans_kgrey_h2o, ans_kgrey_co2, ans_kgrey_ch4
      ! major gas (k-distribution)
      call bilinear_interpK_8gpt_major_gptvec(k48_major_data, imajor(1), pressure, p_ref_index, t_kgas, t_ref_index, ans_kmajor_gptvec) 
!!write(*,*) "kmajor", ans_kmajor_gptvec
      do ig = 1,ngauss_pts(sp)      
        itl = itl + 1         
        tau_gas(itl,ik) = ans_kmajor_gptvec(ig)*ugas(imajor(1)) + (tau_grey(iH2O) + tau_grey(iCO2) + tau_grey(iCH4) - tau_grey(imajor(1)))
      enddo   

      !=====  interval 49:  =====!
      !=====  9650 - 10400 cm-1 
      sp=49
      call bilinear_interpK_grey(k49_grey_data, iH2O, pressure, p_ref_index, t_kgas, t_ref_index, ans_kgrey_h2o) 
      call bilinear_interpK_grey(k49_grey_data, iCO2, pressure, p_ref_index, t_kgas, t_ref_index, ans_kgrey_co2) 
      call bilinear_interpK_grey(k49_grey_data, iCH4, pressure, p_ref_index, t_kgas, t_ref_index, ans_kgrey_ch4)
      tau_grey(iH2O) = ans_kgrey_h2o * ugas(iH2O)
      tau_grey(iCO2) = ans_kgrey_co2 * ugas(iCO2)
      tau_grey(iCH4) = ans_kgrey_ch4 * ugas(iCH4)
      imajor = maxloc(tau_grey)
!write(*,*) ik, sp, gas_name(imajor), tau_grey, ans_kgrey_h2o, ans_kgrey_co2, ans_kgrey_ch4
      ! major gas (k-distribution)
      call bilinear_interpK_8gpt_major_gptvec(k49_major_data, imajor(1), pressure, p_ref_index, t_kgas, t_ref_index, ans_kmajor_gptvec) 
!!write(*,*) "kmajor", ans_kmajor_gptvec
      do ig = 1,ngauss_pts(sp)      
        itl = itl + 1         
        tau_gas(itl,ik) = ans_kmajor_gptvec(ig)*ugas(imajor(1)) + (tau_grey(iH2O) + tau_grey(iCO2) + tau_grey(iCH4) - tau_grey(imajor(1)))
      enddo   

      !=====  interval 50:  =====!
      !==== 10400 - 11220 cm-1
      sp=50
      call bilinear_interpK_grey(k50_grey_data, iH2O, pressure, p_ref_index, t_kgas, t_ref_index, ans_kgrey_h2o) 
      call bilinear_interpK_grey(k50_grey_data, iCO2, pressure, p_ref_index, t_kgas, t_ref_index, ans_kgrey_co2) 
      call bilinear_interpK_grey(k50_grey_data, iCH4, pressure, p_ref_index, t_kgas, t_ref_index, ans_kgrey_ch4)
      tau_grey(iH2O) = ans_kgrey_h2o * ugas(iH2O)
      tau_grey(iCO2) = ans_kgrey_co2 * ugas(iCO2)
      tau_grey(iCH4) = ans_kgrey_ch4 * ugas(iCH4)
      imajor = maxloc(tau_grey)
!write(*,*) ik, sp, gas_name(imajor), tau_grey, ans_kgrey_h2o, ans_kgrey_co2, ans_kgrey_ch4
      ! major gas (k-distribution)
      call bilinear_interpK_8gpt_major_gptvec(k50_major_data, imajor(1), pressure, p_ref_index, t_kgas, t_ref_index, ans_kmajor_gptvec) 
!!write(*,*) "kmajor", ans_kmajor_gptvec
      do ig = 1,ngauss_pts(sp)      
        itl = itl + 1         
        tau_gas(itl,ik) = ans_kmajor_gptvec(ig)*ugas(imajor(1)) + (tau_grey(iH2O) + tau_grey(iCO2) + tau_grey(iCH4) - tau_grey(imajor(1)))
      enddo   

      !=====  interval 51:  =====!
      !=====  11220 - 11870 cm-1 
      sp=51
      call bilinear_interpK_grey(k51_grey_data, iH2O, pressure, p_ref_index, t_kgas, t_ref_index, ans_kgrey_h2o) 
      call bilinear_interpK_grey(k51_grey_data, iCO2, pressure, p_ref_index, t_kgas, t_ref_index, ans_kgrey_co2) 
      call bilinear_interpK_grey(k51_grey_data, iCH4, pressure, p_ref_index, t_kgas, t_ref_index, ans_kgrey_ch4)
      tau_grey(iH2O) = ans_kgrey_h2o * ugas(iH2O)
      tau_grey(iCO2) = ans_kgrey_co2 * ugas(iCO2)
      tau_grey(iCH4) = ans_kgrey_ch4 * ugas(iCH4)
      imajor = maxloc(tau_grey)
!write(*,*) ik, sp, gas_name(imajor), tau_grey, ans_kgrey_h2o, ans_kgrey_co2, ans_kgrey_ch4
      ! major gas (k-distribution)
      call bilinear_interpK_8gpt_major_gptvec(k51_major_data, imajor(1), pressure, p_ref_index, t_kgas, t_ref_index, ans_kmajor_gptvec) 
!!write(*,*) "kmajor", ans_kmajor_gptvec
      do ig = 1,ngauss_pts(sp)      
        itl = itl + 1         
        tau_gas(itl,ik) = ans_kmajor_gptvec(ig)*ugas(imajor(1)) + (tau_grey(iH2O) + tau_grey(iCO2) + tau_grey(iCH4) - tau_grey(imajor(1)))
      enddo   

      !=====  interval 52:  =====!
      !=====  11870 - 12790 cm-1 
      sp=52
      call bilinear_interpK_grey(k52_grey_data, iH2O, pressure, p_ref_index, t_kgas, t_ref_index, ans_kgrey_h2o) 
      call bilinear_interpK_grey(k52_grey_data, iCO2, pressure, p_ref_index, t_kgas, t_ref_index, ans_kgrey_co2) 
      call bilinear_interpK_grey(k52_grey_data, iCH4, pressure, p_ref_index, t_kgas, t_ref_index, ans_kgrey_ch4)
      tau_grey(iH2O) = ans_kgrey_h2o * ugas(iH2O)
      tau_grey(iCO2) = ans_kgrey_co2 * ugas(iCO2)
      tau_grey(iCH4) = ans_kgrey_ch4 * ugas(iCH4)
      imajor = maxloc(tau_grey)
!write(*,*) ik, sp, gas_name(imajor), tau_grey, ans_kgrey_h2o, ans_kgrey_co2, ans_kgrey_ch4
      ! major gas (k-distribution)
      call bilinear_interpK_8gpt_major_gptvec(k52_major_data, imajor(1), pressure, p_ref_index, t_kgas, t_ref_index, ans_kmajor_gptvec) 
!!write(*,*) "kmajor", ans_kmajor_gptvec
      do ig = 1,ngauss_pts(sp)      
        itl = itl + 1         
        tau_gas(itl,ik) = ans_kmajor_gptvec(ig)*ugas(imajor(1)) + (tau_grey(iH2O) + tau_grey(iCO2) + tau_grey(iCH4) - tau_grey(imajor(1)))
      enddo   

      !=====  interval 53:  =====!
      !=====  12790 - 13300 cm-1 
      sp=53
      call bilinear_interpK_grey(k53_grey_data, iH2O, pressure, p_ref_index, t_kgas, t_ref_index, ans_kgrey_h2o) 
      call bilinear_interpK_grey(k53_grey_data, iCO2, pressure, p_ref_index, t_kgas, t_ref_index, ans_kgrey_co2) 
      call bilinear_interpK_grey(k53_grey_data, iCH4, pressure, p_ref_index, t_kgas, t_ref_index, ans_kgrey_ch4)
      tau_grey(iH2O) = ans_kgrey_h2o * ugas(iH2O)
      tau_grey(iCO2) = ans_kgrey_co2 * ugas(iCO2)
      tau_grey(iCH4) = ans_kgrey_ch4 * ugas(iCH4)
      imajor = maxloc(tau_grey)
!write(*,*) ik, sp, gas_name(imajor), tau_grey, ans_kgrey_h2o, ans_kgrey_co2, ans_kgrey_ch4
      ! major gas (k-distribution)
      call bilinear_interpK_8gpt_major_gptvec(k53_major_data, imajor(1), pressure, p_ref_index, t_kgas, t_ref_index, ans_kmajor_gptvec) 
!!write(*,*) "kmajor", ans_kmajor_gptvec
      do ig = 1,ngauss_pts(sp)      
        itl = itl + 1         
        tau_gas(itl,ik) = ans_kmajor_gptvec(ig)*ugas(imajor(1)) + (tau_grey(iH2O) + tau_grey(iCO2) + tau_grey(iCH4) - tau_grey(imajor(1)))
      enddo   

      !=====  interval 54:  =====!
      !=====  13300 - 14470 cm-1 
      sp=54
      call bilinear_interpK_grey(k54_grey_data, iH2O, pressure, p_ref_index, t_kgas, t_ref_index, ans_kgrey_h2o) 
      call bilinear_interpK_grey(k54_grey_data, iCO2, pressure, p_ref_index, t_kgas, t_ref_index, ans_kgrey_co2) 
      call bilinear_interpK_grey(k54_grey_data, iCH4, pressure, p_ref_index, t_kgas, t_ref_index, ans_kgrey_ch4)
      tau_grey(iH2O) = ans_kgrey_h2o * ugas(iH2O)
      tau_grey(iCO2) = ans_kgrey_co2 * ugas(iCO2)
      tau_grey(iCH4) = ans_kgrey_ch4 * ugas(iCH4)
      imajor = maxloc(tau_grey)
!write(*,*) ik, sp, gas_name(imajor), tau_grey, ans_kgrey_h2o, ans_kgrey_co2, ans_kgrey_ch4
      ! major gas (k-distribution)
      call bilinear_interpK_8gpt_major_gptvec(k54_major_data, imajor(1), pressure, p_ref_index, t_kgas, t_ref_index, ans_kmajor_gptvec) 
!!write(*,*) "kmajor", ans_kmajor_gptvec
      do ig = 1,ngauss_pts(sp)      
        itl = itl + 1         
        tau_gas(itl,ik) = ans_kmajor_gptvec(ig)*ugas(imajor(1)) + (tau_grey(iH2O) + tau_grey(iCO2) + tau_grey(iCH4) - tau_grey(imajor(1)))
      enddo   

      !=====  interval 55:  =====!
      !=====  14470 - 15000 cm-1 
      sp=55
      call bilinear_interpK_grey(k55_grey_data, iH2O, pressure, p_ref_index, t_kgas, t_ref_index, ans_kgrey_h2o) 
      call bilinear_interpK_grey(k55_grey_data, iCO2, pressure, p_ref_index, t_kgas, t_ref_index, ans_kgrey_co2) 
      call bilinear_interpK_grey(k55_grey_data, iCH4, pressure, p_ref_index, t_kgas, t_ref_index, ans_kgrey_ch4)
      tau_grey(iH2O) = ans_kgrey_h2o * ugas(iH2O)
      tau_grey(iCO2) = ans_kgrey_co2 * ugas(iCO2)
      tau_grey(iCH4) = ans_kgrey_ch4 * ugas(iCH4)
      imajor = maxloc(tau_grey)
!write(*,*) ik, sp, gas_name(imajor), tau_grey, ans_kgrey_h2o, ans_kgrey_co2, ans_kgrey_ch4
      ! major gas (k-distribution)
      call bilinear_interpK_8gpt_major_gptvec(k55_major_data, imajor(1), pressure, p_ref_index, t_kgas, t_ref_index, ans_kmajor_gptvec) 
!!write(*,*) "kmajor", ans_kmajor_gptvec
      do ig = 1,ngauss_pts(sp)      
        itl = itl + 1         
        tau_gas(itl,ik) = ans_kmajor_gptvec(ig)*ugas(imajor(1)) + (tau_grey(iH2O) + tau_grey(iCO2) + tau_grey(iCH4) - tau_grey(imajor(1)))
      enddo   

      !=====  interval 56:  =====!
      !=====  15000 - 16000 cm-1 
      sp=56
      call bilinear_interpK_grey(k56_grey_data, iH2O, pressure, p_ref_index, t_kgas, t_ref_index, ans_kgrey_h2o) 
      call bilinear_interpK_grey(k56_grey_data, iCO2, pressure, p_ref_index, t_kgas, t_ref_index, ans_kgrey_co2) 
      call bilinear_interpK_grey(k56_grey_data, iCH4, pressure, p_ref_index, t_kgas, t_ref_index, ans_kgrey_ch4)
      tau_grey(iH2O) = ans_kgrey_h2o * ugas(iH2O)
      tau_grey(iCO2) = ans_kgrey_co2 * ugas(iCO2)
      tau_grey(iCH4) = ans_kgrey_ch4 * ugas(iCH4)
      imajor = maxloc(tau_grey)
!write(*,*) ik, sp, gas_name(imajor), tau_grey, ans_kgrey_h2o, ans_kgrey_co2, ans_kgrey_ch4
      ! major gas (k-distribution)
      call bilinear_interpK_8gpt_major_gptvec(k56_major_data, imajor(1), pressure, p_ref_index, t_kgas, t_ref_index, ans_kmajor_gptvec) 
!!write(*,*) "kmajor", ans_kmajor_gptvec
      do ig = 1,ngauss_pts(sp)      
        itl = itl + 1         
        tau_gas(itl,ik) = ans_kmajor_gptvec(ig)*ugas(imajor(1)) + (tau_grey(iH2O) + tau_grey(iCO2) + tau_grey(iCH4) - tau_grey(imajor(1)))
      enddo   

      !=====  interval 57:  =====!
      !=====  16000 - 16528 cm-1 
      sp=57
      call bilinear_interpK_grey(k57_grey_data, iH2O, pressure, p_ref_index, t_kgas, t_ref_index, ans_kgrey_h2o) 
      call bilinear_interpK_grey(k57_grey_data, iCO2, pressure, p_ref_index, t_kgas, t_ref_index, ans_kgrey_co2) 
      call bilinear_interpK_grey(k57_grey_data, iCH4, pressure, p_ref_index, t_kgas, t_ref_index, ans_kgrey_ch4)
      tau_grey(iH2O) = ans_kgrey_h2o * ugas(iH2O)
      tau_grey(iCO2) = ans_kgrey_co2 * ugas(iCO2)
      tau_grey(iCH4) = ans_kgrey_ch4 * ugas(iCH4)
      imajor = maxloc(tau_grey)
!write(*,*) ik, sp, gas_name(imajor), tau_grey, ans_kgrey_h2o, ans_kgrey_co2, ans_kgrey_ch4
      ! major gas (k-distribution)
      call bilinear_interpK_8gpt_major_gptvec(k57_major_data, imajor(1), pressure, p_ref_index, t_kgas, t_ref_index, ans_kmajor_gptvec) 
!!write(*,*) "kmajor", ans_kmajor_gptvec
      do ig = 1,ngauss_pts(sp)      
        itl = itl + 1         
        tau_gas(itl,ik) = ans_kmajor_gptvec(ig)*ugas(imajor(1)) + (tau_grey(iH2O) + tau_grey(iCO2) + tau_grey(iCH4) - tau_grey(imajor(1)))
      enddo   

      !=====  interval 58:  =====!
      !=====  16528 - 17649 cm-1 
      sp=58
      call bilinear_interpK_grey(k58_grey_data, iH2O, pressure, p_ref_index, t_kgas, t_ref_index, ans_kgrey_h2o) 
      call bilinear_interpK_grey(k58_grey_data, iCO2, pressure, p_ref_index, t_kgas, t_ref_index, ans_kgrey_co2) 
      call bilinear_interpK_grey(k58_grey_data, iCH4, pressure, p_ref_index, t_kgas, t_ref_index, ans_kgrey_ch4)
      tau_grey(iH2O) = ans_kgrey_h2o * ugas(iH2O)
      tau_grey(iCO2) = ans_kgrey_co2 * ugas(iCO2)
      tau_grey(iCH4) = ans_kgrey_ch4 * ugas(iCH4)
      imajor = maxloc(tau_grey)
!write(*,*) ik, sp, gas_name(imajor), tau_grey, ans_kgrey_h2o, ans_kgrey_co2, ans_kgrey_ch4
      ! major gas (k-distribution)
      call bilinear_interpK_8gpt_major_gptvec(k58_major_data, imajor(1), pressure, p_ref_index, t_kgas, t_ref_index, ans_kmajor_gptvec) 
!!write(*,*) "kmajor", ans_kmajor_gptvec
      do ig = 1,ngauss_pts(sp)      
        itl = itl + 1         
        tau_gas(itl,ik) = ans_kmajor_gptvec(ig)*ugas(imajor(1)) + (tau_grey(iH2O) + tau_grey(iCO2) + tau_grey(iCH4) - tau_grey(imajor(1)))
      enddo   

      !=====  interval 59:  =====!
      !=====  17649 - 18198 cm-1 
      sp=59
      call bilinear_interpK_grey(k59_grey_data, iH2O, pressure, p_ref_index, t_kgas, t_ref_index, ans_kgrey_h2o) 
      call bilinear_interpK_grey(k59_grey_data, iCO2, pressure, p_ref_index, t_kgas, t_ref_index, ans_kgrey_co2) 
      call bilinear_interpK_grey(k59_grey_data, iCH4, pressure, p_ref_index, t_kgas, t_ref_index, ans_kgrey_ch4)
      tau_grey(iH2O) = ans_kgrey_h2o * ugas(iH2O)
      tau_grey(iCO2) = ans_kgrey_co2 * ugas(iCO2)
      tau_grey(iCH4) = ans_kgrey_ch4 * ugas(iCH4)
      imajor = maxloc(tau_grey)
!write(*,*) ik, sp, gas_name(imajor), tau_grey, ans_kgrey_h2o, ans_kgrey_co2, ans_kgrey_ch4
      ! major gas (k-distribution)
      call bilinear_interpK_8gpt_major_gptvec(k59_major_data, imajor(1), pressure, p_ref_index, t_kgas, t_ref_index, ans_kmajor_gptvec) 
!!write(*,*) "kmajor", ans_kmajor_gptvec
      do ig = 1,ngauss_pts(sp)      
        itl = itl + 1         
        tau_gas(itl,ik) = ans_kmajor_gptvec(ig)*ugas(imajor(1)) + (tau_grey(iH2O) + tau_grey(iCO2) + tau_grey(iCH4) - tau_grey(imajor(1)))
      enddo   

      !=====  interval 60:  =====!
      !=====  18198 - 18518 cm-1 
      sp=60
      call bilinear_interpK_grey(k60_grey_data, iH2O, pressure, p_ref_index, t_kgas, t_ref_index, ans_kgrey_h2o) 
      call bilinear_interpK_grey(k60_grey_data, iCO2, pressure, p_ref_index, t_kgas, t_ref_index, ans_kgrey_co2) 
      call bilinear_interpK_grey(k60_grey_data, iCH4, pressure, p_ref_index, t_kgas, t_ref_index, ans_kgrey_ch4)
      tau_grey(iH2O) = ans_kgrey_h2o * ugas(iH2O)
      tau_grey(iCO2) = ans_kgrey_co2 * ugas(iCO2)
      tau_grey(iCH4) = ans_kgrey_ch4 * ugas(iCH4)
      imajor = maxloc(tau_grey)
!write(*,*) ik, sp, gas_name(imajor), tau_grey, ans_kgrey_h2o, ans_kgrey_co2, ans_kgrey_ch4
      ! major gas (k-distribution)
      call bilinear_interpK_8gpt_major_gptvec(k60_major_data, imajor(1), pressure, p_ref_index, t_kgas, t_ref_index, ans_kmajor_gptvec) 
!!write(*,*) "kmajor", ans_kmajor_gptvec
      do ig = 1,ngauss_pts(sp)      
        itl = itl + 1         
        tau_gas(itl,ik) = ans_kmajor_gptvec(ig)*ugas(imajor(1)) + (tau_grey(iH2O) + tau_grey(iCO2) + tau_grey(iCH4) - tau_grey(imajor(1)))
      enddo   

      !=====  interval 61:  =====!
      !=====  18518 - 22222 cm-1 
      sp=61
      call bilinear_interpK_grey(k61_grey_data, iH2O, pressure, p_ref_index, t_kgas, t_ref_index, ans_kgrey_h2o) 
      call bilinear_interpK_grey(k61_grey_data, iCO2, pressure, p_ref_index, t_kgas, t_ref_index, ans_kgrey_co2) 
      call bilinear_interpK_grey(k61_grey_data, iCH4, pressure, p_ref_index, t_kgas, t_ref_index, ans_kgrey_ch4)
      tau_grey(iH2O) = ans_kgrey_h2o * ugas(iH2O)
      tau_grey(iCO2) = ans_kgrey_co2 * ugas(iCO2)
      tau_grey(iCH4) = ans_kgrey_ch4 * ugas(iCH4)
      imajor = maxloc(tau_grey)
!write(*,*) ik, sp, gas_name(imajor), tau_grey, ans_kgrey_h2o, ans_kgrey_co2, ans_kgrey_ch4
      ! major gas (k-distribution)
      call bilinear_interpK_8gpt_major_gptvec(k61_major_data, imajor(1), pressure, p_ref_index, t_kgas, t_ref_index, ans_kmajor_gptvec) 
!!write(*,*) "kmajor", ans_kmajor_gptvec
      do ig = 1,ngauss_pts(sp)      
        itl = itl + 1         
        tau_gas(itl,ik) = ans_kmajor_gptvec(ig)*ugas(imajor(1)) + (tau_grey(iH2O) + tau_grey(iCO2) + tau_grey(iCH4) - tau_grey(imajor(1)))
      enddo   

      !=====  interval 62:  =====!
      !=====  22222 - 25641 cm-1 
      sp=62
      call bilinear_interpK_grey(k62_grey_data, iH2O, pressure, p_ref_index, t_kgas, t_ref_index, ans_kgrey_h2o) 
      call bilinear_interpK_grey(k62_grey_data, iCO2, pressure, p_ref_index, t_kgas, t_ref_index, ans_kgrey_co2) 
      call bilinear_interpK_grey(k62_grey_data, iCH4, pressure, p_ref_index, t_kgas, t_ref_index, ans_kgrey_ch4)
      tau_grey(iH2O) = ans_kgrey_h2o * ugas(iH2O)
      tau_grey(iCO2) = ans_kgrey_co2 * ugas(iCO2)
      tau_grey(iCH4) = ans_kgrey_ch4 * ugas(iCH4)
      imajor = maxloc(tau_grey)
!write(*,*) ik, sp, gas_name(imajor), tau_grey, ans_kgrey_h2o, ans_kgrey_co2, ans_kgrey_ch4
      ! major gas (k-distribution)
      call bilinear_interpK_8gpt_major_gptvec(k62_major_data, imajor(1), pressure, p_ref_index, t_kgas, t_ref_index, ans_kmajor_gptvec) 
!!write(*,*) "kmajor", ans_kmajor_gptvec
      do ig = 1,ngauss_pts(sp)      
        itl = itl + 1         
        tau_gas(itl,ik) = ans_kmajor_gptvec(ig)*ugas(imajor(1)) + (tau_grey(iH2O) + tau_grey(iCO2) + tau_grey(iCH4) - tau_grey(imajor(1)))
      enddo   

      !=====  interval 63:  =====!
      !=====  25641 - 29308 cm-1 
      sp=63
      call bilinear_interpK_grey(k63_grey_data, iH2O, pressure, p_ref_index, t_kgas, t_ref_index, ans_kgrey_h2o) 
      call bilinear_interpK_grey(k63_grey_data, iCO2, pressure, p_ref_index, t_kgas, t_ref_index, ans_kgrey_co2) 
      call bilinear_interpK_grey(k63_grey_data, iCH4, pressure, p_ref_index, t_kgas, t_ref_index, ans_kgrey_ch4)
      tau_grey(iH2O) = ans_kgrey_h2o * ugas(iH2O)
      tau_grey(iCO2) = ans_kgrey_co2 * ugas(iCO2)
      tau_grey(iCH4) = ans_kgrey_ch4 * ugas(iCH4)
      imajor = maxloc(tau_grey)
!write(*,*) ik, sp, gas_name(imajor), tau_grey, ans_kgrey_h2o, ans_kgrey_co2, ans_kgrey_ch4
      ! major gas (k-distribution)
      call bilinear_interpK_8gpt_major_gptvec(k63_major_data, imajor(1), pressure, p_ref_index, t_kgas, t_ref_index, ans_kmajor_gptvec) 
!!write(*,*) "kmajor", ans_kmajor_gptvec
      do ig = 1,ngauss_pts(sp)      
        itl = itl + 1         
        tau_gas(itl,ik) = ans_kmajor_gptvec(ig)*ugas(imajor(1)) + (tau_grey(iH2O) + tau_grey(iCO2) + tau_grey(iCH4) - tau_grey(imajor(1)))
      enddo   

      !=====  interval 64:  =====!
      !=====  29308 - 30376 cm-1 
      sp=64
      call bilinear_interpK_grey(k64_grey_data, iH2O, pressure, p_ref_index, t_kgas, t_ref_index, ans_kgrey_h2o) 
      call bilinear_interpK_grey(k64_grey_data, iCO2, pressure, p_ref_index, t_kgas, t_ref_index, ans_kgrey_co2) 
      call bilinear_interpK_grey(k64_grey_data, iCH4, pressure, p_ref_index, t_kgas, t_ref_index, ans_kgrey_ch4)
      tau_grey(iH2O) = ans_kgrey_h2o * ugas(iH2O)
      tau_grey(iCO2) = ans_kgrey_co2 * ugas(iCO2)
      tau_grey(iCH4) = ans_kgrey_ch4 * ugas(iCH4)
      imajor = maxloc(tau_grey)
!write(*,*) ik, sp, gas_name(imajor), tau_grey, ans_kgrey_h2o, ans_kgrey_co2, ans_kgrey_ch4
      ! major gas (k-distribution)
      call bilinear_interpK_8gpt_major_gptvec(k64_major_data, imajor(1), pressure, p_ref_index, t_kgas, t_ref_index, ans_kmajor_gptvec) 
!!write(*,*) "kmajor", ans_kmajor_gptvec
      do ig = 1,ngauss_pts(sp)      
        itl = itl + 1         
        tau_gas(itl,ik) = ans_kmajor_gptvec(ig)*ugas(imajor(1)) + (tau_grey(iH2O) + tau_grey(iCO2) + tau_grey(iCH4) - tau_grey(imajor(1)))
      enddo   

      !=====  interval 65:  =====!
      !=====  30376 - 32562 cm-1 
      sp=65
      call bilinear_interpK_grey(k65_grey_data, iH2O, pressure, p_ref_index, t_kgas, t_ref_index, ans_kgrey_h2o) 
      call bilinear_interpK_grey(k65_grey_data, iCO2, pressure, p_ref_index, t_kgas, t_ref_index, ans_kgrey_co2) 
      call bilinear_interpK_grey(k65_grey_data, iCH4, pressure, p_ref_index, t_kgas, t_ref_index, ans_kgrey_ch4)
      tau_grey(iH2O) = ans_kgrey_h2o * ugas(iH2O)
      tau_grey(iCO2) = ans_kgrey_co2 * ugas(iCO2)
      tau_grey(iCH4) = ans_kgrey_ch4 * ugas(iCH4)
      imajor = maxloc(tau_grey)
!write(*,*) ik, sp, gas_name(imajor), tau_grey, ans_kgrey_h2o, ans_kgrey_co2, ans_kgrey_ch4
      ! major gas (k-distribution)
      call bilinear_interpK_8gpt_major_gptvec(k65_major_data, imajor(1), pressure, p_ref_index, t_kgas, t_ref_index, ans_kmajor_gptvec) 
!!write(*,*) "kmajor", ans_kmajor_gptvec
      do ig = 1,ngauss_pts(sp)      
        itl = itl + 1         
        tau_gas(itl,ik) = ans_kmajor_gptvec(ig)*ugas(imajor(1)) + (tau_grey(iH2O) + tau_grey(iCO2) + tau_grey(iCH4) - tau_grey(imajor(1)))
      enddo   

      !=====  interval 66:  =====!
      !=====  32562 - 35087 cm-1 
      sp=66
      call bilinear_interpK_grey(k66_grey_data, iH2O, pressure, p_ref_index, t_kgas, t_ref_index, ans_kgrey_h2o) 
      call bilinear_interpK_grey(k66_grey_data, iCO2, pressure, p_ref_index, t_kgas, t_ref_index, ans_kgrey_co2) 
      call bilinear_interpK_grey(k66_grey_data, iCH4, pressure, p_ref_index, t_kgas, t_ref_index, ans_kgrey_ch4)
      tau_grey(iH2O) = ans_kgrey_h2o * ugas(iH2O)
      tau_grey(iCO2) = ans_kgrey_co2 * ugas(iCO2)
      tau_grey(iCH4) = ans_kgrey_ch4 * ugas(iCH4)
      imajor = maxloc(tau_grey)
!write(*,*) ik, sp, gas_name(imajor), tau_grey, ans_kgrey_h2o, ans_kgrey_co2, ans_kgrey_ch4
      ! major gas (k-distribution)
      call bilinear_interpK_8gpt_major_gptvec(k66_major_data, imajor(1), pressure, p_ref_index, t_kgas, t_ref_index, ans_kmajor_gptvec) 
!!write(*,*) "kmajor", ans_kmajor_gptvec
      do ig = 1,ngauss_pts(sp)      
        itl = itl + 1         
        tau_gas(itl,ik) = ans_kmajor_gptvec(ig)*ugas(imajor(1)) + (tau_grey(iH2O) + tau_grey(iCO2) + tau_grey(iCH4) - tau_grey(imajor(1)))
      enddo   

      !=====  interval 67:  =====!
      !=====  35087 - 36363 cm-1 
      sp=67
      call bilinear_interpK_grey(k67_grey_data, iH2O, pressure, p_ref_index, t_kgas, t_ref_index, ans_kgrey_h2o) 
      call bilinear_interpK_grey(k67_grey_data, iCO2, pressure, p_ref_index, t_kgas, t_ref_index, ans_kgrey_co2) 
      call bilinear_interpK_grey(k67_grey_data, iCH4, pressure, p_ref_index, t_kgas, t_ref_index, ans_kgrey_ch4)
      tau_grey(iH2O) = ans_kgrey_h2o * ugas(iH2O)
      tau_grey(iCO2) = ans_kgrey_co2 * ugas(iCO2)
      tau_grey(iCH4) = ans_kgrey_ch4 * ugas(iCH4)
!write(*,*) ik, sp, gas_name(imajor), tau_grey, ans_kgrey_h2o, ans_kgrey_co2, ans_kgrey_ch4
      imajor = maxloc(tau_grey)
      ! major gas (k-distribution)
      call bilinear_interpK_8gpt_major_gptvec(k67_major_data, imajor(1), pressure, p_ref_index, t_kgas, t_ref_index, ans_kmajor_gptvec) 
!!write(*,*) "kmajor", ans_kmajor_gptvec
      do ig = 1,ngauss_pts(sp)      
        itl = itl + 1         
        tau_gas(itl,ik) = ans_kmajor_gptvec(ig)*ugas(imajor(1)) + (tau_grey(iH2O) + tau_grey(iCO2) + tau_grey(iCH4) - tau_grey(imajor(1)))
      enddo   

      !=====  interval 68:  =====!
      !=====  36363 - 42087 cm-1 
      sp=68
      call bilinear_interpK_grey(k68_grey_data, iH2O, pressure, p_ref_index, t_kgas, t_ref_index, ans_kgrey_h2o) 
      call bilinear_interpK_grey(k68_grey_data, iCO2, pressure, p_ref_index, t_kgas, t_ref_index, ans_kgrey_co2) 
      call bilinear_interpK_grey(k68_grey_data, iCH4, pressure, p_ref_index, t_kgas, t_ref_index, ans_kgrey_ch4)
      tau_grey(iH2O) = ans_kgrey_h2o * ugas(iH2O)
      tau_grey(iCO2) = ans_kgrey_co2 * ugas(iCO2)
      tau_grey(iCH4) = ans_kgrey_ch4 * ugas(iCH4)
      imajor = maxloc(tau_grey)
!write(*,*) ik, sp, gas_name(imajor), tau_grey, ans_kgrey_h2o, ans_kgrey_co2, ans_kgrey_ch4
      ! major gas (k-distribution)
!!write(*,*) "kmajor", ans_kmajor_gptvec
      call bilinear_interpK_8gpt_major_gptvec(k68_major_data, imajor(1), pressure, p_ref_index, t_kgas, t_ref_index, ans_kmajor_gptvec) 
      do ig = 1,ngauss_pts(sp)      
        itl = itl + 1         
        tau_gas(itl,ik) = ans_kmajor_gptvec(ig)*ugas(imajor(1)) + (tau_grey(iH2O) + tau_grey(iCO2) + tau_grey(iCH4) - tau_grey(imajor(1)))
      enddo   

      !
      !  Water Vapor Self Continuum 
      !
      !===== MT_CKD =====!

      ! For water vapor self continuum, find temperature index
      ! find the reference temperature value, exit if temperature less than minimum of grid
      t_ref_index_s = kmtckd_ntemp  ! index of reference temperature
      t_h2os = temperature
      do  
        if (t_ref_index_s .le. 1) exit                   ! exit if temperature less than minimum grid
        if (t_h2os .gt. tgrid_mtckd(kmtckd_ntemp)) then  ! temperature greater than grid max
          t_h2os = tgrid_mtckd(t_ref_index_s)   ! set t to max grid value
          exit                                              ! exit
        endif
        if ((tgrid_mtckd(t_ref_index_s) .le. t_h2os)) exit ! found reference temperature
        t_ref_index_s = t_ref_index_s - 1                   ! increment reference temperature
      enddo
      ! if temperature less than minimum of grid, force reference to minimum grid value
      ! force reference temperature for interpolation to minimum temperature in tgrid
      if (t_ref_index_s .lt. 1) then
        t_ref_index_s = 1
        t_h2os = tgrid_mtckd(t_ref_index_s)
      endif

!write(*,*) "t_h2os", t_h2os
      !!  Do MT_CKD self_continuum           
      itc=0
      call interpH2Omtckd(kh2oself_avg_mtckd, t_h2os, t_ref_index_s, ans_h2os_avg)
      call interpH2Omtckd(kh2ofrgn_avg_mtckd, t_h2os, t_ref_index_s, ans_h2of_avg)
      call interpH2Omtckd_ng(kh2oself_mtckd, t_h2os, t_ref_index_s, ans_h2os)
      call interpH2Omtckd_ng(kh2ofrgn_mtckd, t_h2os, t_ref_index_s, ans_h2of)

      do iw = iwbeg,iwend     ! loop over bands      
        ! apply H2O self continuum from infinity down to ~2 microns, 5000 cm-1
        ! extension further results in over-estimation of absorption compared to recent literature
!        if (iw.gt.41) ans_h2os(iw) = 0.


        mtckd_self_tau(iw, ik) = ans_h2os_avg(iw)*u_h2o*(273.15/temperature)*(ppH2O/1013.250)
        mtckd_frgn_tau(iw, ik) = ans_h2of_avg(iw)*u_h2o*(273.15/temperature)*((pmid(ik)-ppH2O)/1013.250)

        do ig=1, ngauss_pts(iw)
          itc=itc+1
! ngauss MTCKD 
!write(*,*) "mtckd, ans_h2os", ig, iw, ans_h2os(ig, iw)
!write(*,*) "mtckd, ans_h2of", ig, iw, ans_h2of(ig, iw)
!          tau_gas(itc,ik) = tau_gas(itc,ik) + (ans_h2os(ig,iw)*amaH2O + ans_h2of(ig,iw)*amaFRGN) * u_h2o
 

! band averaged value        
!          tau_gas(itc,ik) = tau_gas(itc,ik) + mtckd_self_tau(iw,ik) + mtckd_frgn_tau(iw,ik)

        enddo

      enddo    ! close band loop

      !  
      !  Collision Induced Absorption
      !

      !!====  Calculate H2O-H2O CIA  ====!! 
      !! Villaneuva and Koffman H2O-H2O CIA Method
      if (u_h2o .gt. 0) then
        t_ref_index_h2oh2o = kh2oh2o_ntemp
        t_h2oh2o = temperature
        do
          if (t_ref_index_h2oh2o .le. 1) exit                  ! exit if temperature less than minimum grid
          if (t_h2oh2o .gt. tgrid_h2oh2o(kh2oh2o_ntemp)) then  ! temperature greater than grid max 
            t_h2oh2o = tgrid_h2oh2o(t_ref_index_h2oh2o)   ! set t to max grid value
            exit                                              ! exit 
          endif
          if ((tgrid_h2oh2o(t_ref_index_h2oh2o) .le. t_h2oh2o)) exit ! found reference temperature 
          t_ref_index_h2oh2o = t_ref_index_h2oh2o - 1                   ! increment reference temperature
        enddo
        ! if temperature less than minimum of grid, force reference to minimum grid value 
        ! force reference temperature for interpolation to minimum temperature in tgrid  
        if (t_ref_index_h2oh2o .lt. 1) then
          t_ref_index_h2oh2o = 1
          t_h2oh2o = tgrid_h2oh2o(t_ref_index_h2oh2o)
        endif
        call interpH2OH2Ocia(kh2oh2o, t_h2oh2o, t_ref_index_h2oh2o, ans_cia)
        do iw=iwbeg,iwend      ! loop over bands         
! !   if (iw.gt.39) ans_cia(iw) = 0.
          tau_h2oh2o_cia(iw,ik) = ans_cia(iw) * amaH2O * amaH2O * pathlength(ik)
          write(*,*) "H2O-H2O CIA",iw, ans, tau_h2oh2o_cia(iw,ik)
        enddo
      endif

!      !!====  Calculate H2O-N2 CIA  ====!! 
!      !! Villaneuva and Koffman H2O-N2 CIA Method
      if (u_h2o .gt. 0 .and. u_n2 .gt. 0) then
        t_ref_index_h2on2 = kh2on2_ntemp
        t_h2on2 = temperature
        do
          if (t_ref_index_h2on2 .le. 1) exit                  ! exit if temperature less than minimum grid
          if (t_h2on2 .gt. tgrid_h2on2(kh2on2_ntemp)) then  ! temperature greater than grid max 
            t_h2oh2o = tgrid_h2on2(t_ref_index_h2on2)   ! set t to max grid value
            exit                                              ! exit 
          endif
          if ((tgrid_h2on2(t_ref_index_h2on2) .le. t_h2on2)) exit ! found reference temperature 
          t_ref_index_h2on2 = t_ref_index_h2on2 - 1                   ! increment reference temperature
        enddo
        ! if temperature less than minimum of grid, force reference to minimum grid value 
        ! force reference temperature for interpolation to minimum temperature in tgrid  
        if (t_ref_index_h2on2 .lt. 1) then
          t_ref_index_h2on2 = 1
          t_h2on2 = tgrid_h2on2(t_ref_index_h2on2)
        endif
        call interpH2ON2cia(kh2on2, t_h2on2, t_ref_index_h2on2, ans_cia)
        do iw=iwbeg,iwend      ! loop over bands         
   !!if (iw.gt.38) ans_cia(iw) = 0.
          tau_h2on2_cia(iw,ik) = ans_cia(iw) * amaH2O * amaFRGN * pathlength(ik)
          write(*,*) "H2O-N2 CIA",iw, ans, tau_h2on2_cia(iw,ik)
        enddo
      endif

!! CIA and MTCKD comparison
  do iw=iwbeg, iwend
    write(*,*) "---------------------------------------------------"
    write(*,*) "H2O-N2,  FRGN",iw, tau_h2on2_cia(iw,ik), mtckd_frgn_tau(iw,ik)
    write(*,*) "H2O-H2O, SELF",iw, tau_h2oh2o_cia(iw,ik), mtckd_self_tau(iw,ik)
  enddo

!!===========================


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

!          if (iw .lt. 14) ans_cia(iw) = ans_cia(iw)*10.
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
                                               +  tau_co2co2_sw_cia(iw,ik) + tau_co2co2_lw_cia(iw,ik) &
                                               +  tau_h2oh2o_cia(iw,ik) +  tau_h2on2_cia(iw,ik)
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
        ! Vardavas and Carter (1984), Allen (1976) N2, CO2 Rayleigh scattering
        !
        ! Rayleigh scattering for CO2, N2
        depolCO2 = (6+3*delCO2)/(6-7*delCO2)
        depolN2  = (6+3*delN2)/(6-7*delN2)
        allenCO2 = (1.0E-5*raylA_CO2*(1.0+1.0E-3*raylB_CO2/wl**2))**2
        allenN2 = (1.0E-5*raylA_N2*(1.0+1.0E-3*raylB_N2/wl**2))**2
        sigmaRaylCO2 = 4.577E-21/wl**4*depolCO2*allenCO2
        sigmaRaylN2 = 4.577E-21/wl**4*depolN2*allenN2         
        !  Rayleigh scattering from H2O 
        ns = (5791817./(238.0185-(1./wl)**2) + 167909./(57.362-(1./wl)**2))/1.0E8  ! Bucholtz (1995)  !new
        r = 0.85*ns 
        depolH2O = (6+3*delH2O)/(6-7*delH2O)
        sigmaRaylH2O = 4.577e-21*depolH2O*(r**2)/(wl**4)  !new
	! Rayleigh scattering from H2, Dalgarno & Williams 1962, ApJ, 136, 690D
      	sigmaRaylH2 = 8.14e-13/(wla**4) + 1.28e-6/(wla**6) + 1.61/(wla**8)
        ! Total Rayleigh scattering
        tau_ray(iw,ik) = sigmaRaylCO2*u_co2 + sigmaRaylN2*u_n2 + sigmaRaylH2O*u_h2o + sigmaRaylH2*u_h2
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
