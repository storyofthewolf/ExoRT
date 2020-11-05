

module rad_interp_mod
! version n68equiv
! Contains subroutines used for interpolating K coefficient and 
! cloud optical data in radiation.F90

  use shr_kind_mod,   only: r8 => shr_kind_r8
  use radgrid

  implicit none
  private
  save

!------------------------------------------------------------------------
!
!  Public Interfaces
!

  public :: bilinear_interpK_8gpt_major_gptvec
  public :: bilinear_interpK_grey
  public :: interpH2Omtckd
  public :: interpolate_cld
  public :: interpH2H2cia
  public :: interpN2H2cia
  public :: interpN2N2cia
  public :: interpCO2CO2cia_lw
  public :: interpCO2CO2cia_sw
  public :: interpCO2H2cia
  public :: interpCO2CH4cia

!========================================================================
contains
!========================================================================


!============================================================================

  subroutine bilinear_interpK_8gpt_major_gptvec(Kdata, imaj, pressure, p_ref_index, temperature, t_ref_index, ans) 

!------------------------------------------------------------------------
!
! Purpose: Linearly interpolate k-coefficient from reference P,T
!                                        
!------------------------------------------------------------------------
   
    implicit none
   
!------------------------------------------------------------------------
!
! Input Arguments
!
    real(r8), intent(in) ::Kdata(nspecies, 8, kc_npress, kc_ntemp)
!    real(r8), intent(in) :: Kdata(*, *, *)
    integer, intent(in) :: imaj
!    integer, intent(in) :: ig
    integer, intent(inout) :: p_ref_index
    integer, intent(inout) :: t_ref_index

    real(r8), intent(in) :: pressure  
    real(r8), intent(in) :: temperature
    real(r8), intent(out) :: ans(8)

!------------------------------------------------------------------------
!
! Local Variables
!
    integer :: p_ref_indexp1
    integer :: t_ref_indexp1
    integer :: ik, igi
   
    real(r8) :: u_col
    real(r8) :: press
    real(r8) :: temp
    real(r8) :: weight
    real(r8) :: onemp
    real(r8) :: onemt
    real(r8) :: onemw
    real(r8), dimension(4) :: vbi
    real(r8) :: ans1   
 
    logical :: interpp
    logical :: interpt

!------------------------------------------------------------------------
!
! Start Code
!
    !write(*,*) "-------- bilinear_interpK_8gpt ----------"
    p_ref_indexp1 = p_ref_index + 1    
    t_ref_indexp1 = t_ref_index + 1

    vbi(:) = 0.
    ans = 0.

    ! find the interpolating factor for the various parameters.  at the tops of grids, extrapolate
    if (p_ref_index .eq. kc_npress) then
      p_ref_index = p_ref_index - 1
      p_ref_indexp1 = p_ref_index + 1
      press = (pressure - log10pgrid(p_ref_index)) / &
              (log10pgrid(p_ref_indexp1) - log10pgrid(p_ref_index))      
    else
      press = (pressure - log10pgrid(p_ref_index)) / &
              (log10pgrid(p_ref_indexp1) - log10pgrid(p_ref_index))
    endif
   
!    write(*,*) "------------------------------------------------------" 
!    write(*,*) "pressure", pressure
!    write(*,*) "p_ref_index", p_ref_index
!    write(*,*) "reference", log10pgrid(p_ref_index),log10pgrid(p_ref_indexp1)
!    write(*,*) "interp_press", press
  
    if (t_ref_index .eq. kc_ntemp) then
      t_ref_index = t_ref_index - 1
      t_ref_indexp1 = t_ref_index + 1
      temp = (temperature - tgrid(t_ref_index))/(tgrid(t_ref_indexp1) - tgrid(t_ref_index))
    else
      temp = (temperature - tgrid(t_ref_index))/(tgrid(t_ref_indexp1) - tgrid(t_ref_index))
    endif

!    write(*,*) "temperature", temperature
!    write(*,*) "tgrid(t_ref_index)", tgrid(t_ref_index)
!    write(*,*) "tgrid(t_ref_indexp1)", tgrid(t_ref_indexp1)
!    write(*,*) "interp_temp", temp


    ! perform bilinear interpolation between P,T

    ans(:) = 0.0 
    do igi = 1, 8 

      vbi(1) = kdata(imaj, igi, p_ref_index,   t_ref_index)
      vbi(2) = kdata(imaj, igi, p_ref_indexp1, t_ref_index)
      vbi(3) = kdata(imaj, igi, p_ref_index,   t_ref_indexp1)
      vbi(4) = kdata(imaj, igi, p_ref_indexp1, t_ref_indexp1)

      !    write(*,*) "V", vbi(:)
 
      onemp = 1. - press
      onemt = 1. - temp
 
      ! bilinear  
      ans(igi) = vbi(1)*onemp*onemt &
               + vbi(2)*press*onemt &
               + vbi(3)*onemp*temp &
               + vbi(4)*press*temp
    enddo
!    write(*,*) "------------------------------------------------------" 
     !write(*,*) "gptvec", ans
    return

  end subroutine bilinear_interpK_8gpt_major_gptvec

!============================================================================

  subroutine bilinear_interpK_grey(Kdata, igas, pressure, p_ref_index, temperature, t_ref_index, ans) 

!------------------------------------------------------------------------
!
! Purpose: Linearly interpolate kgrey absorption coefficient in each band 
!          from reference P,T
!                                        
!------------------------------------------------------------------------
   
    implicit none
   
!------------------------------------------------------------------------
!
! Input Arguments
!
    real(r8), intent(in) :: Kdata(nspecies, kc_npress, kc_ntemp)
    integer, intent(in) :: igas
    integer, intent(inout) :: p_ref_index
    integer, intent(inout) :: t_ref_index

    real(r8), intent(in) :: pressure  
    real(r8), intent(in) :: temperature
    real(r8), intent(out) :: ans

!------------------------------------------------------------------------
!
! Local Variables
!
    integer :: p_ref_indexp1
    integer :: t_ref_indexp1
    integer :: ik
   
    real(r8) :: u_col
    real(r8) :: press
    real(r8) :: temp
    real(r8) :: weight
    real(r8) :: onemp
    real(r8) :: onemt
    real(r8) :: onemw
    real(r8), dimension(4) :: vbi
    real(r8) :: ans1   
 
    logical :: interpp
    logical :: interpt

!------------------------------------------------------------------------
!
! Start Code
!
!    write(*,*) "-------- bilinear_interpK_grey ----------"
    p_ref_indexp1 = p_ref_index + 1    
    t_ref_indexp1 = t_ref_index + 1

    vbi(:) = 0.
    ans = 0.

    ! find the interpolating factor for the various parameters.  at the tops of grids, extrapolate
    if (p_ref_index .eq. kc_npress) then
      p_ref_index = p_ref_index - 1
      p_ref_indexp1 = p_ref_index + 1
      press = (pressure - log10pgrid(p_ref_index)) / &
              (log10pgrid(p_ref_indexp1) - log10pgrid(p_ref_index))      
    else
      press = (pressure - log10pgrid(p_ref_index)) / &
              (log10pgrid(p_ref_indexp1) - log10pgrid(p_ref_index))
    endif
   
!    write(*,*) "------------------------------------------------------" 
!    write(*,*) "pressure", pressure, 10.**(pressure)
!    write(*,*) "p_ref_index", p_ref_index
!    write(*,*) "reference", log10pgrid(p_ref_index),log10pgrid(p_ref_indexp1)
!    write(*,*) "interp_press", press
  
    if (t_ref_index .eq. kc_ntemp) then
      t_ref_index = t_ref_index - 1
      t_ref_indexp1 = t_ref_index + 1
      temp = (temperature - tgrid(t_ref_index))/(tgrid(t_ref_indexp1) - tgrid(t_ref_index))
    else
      temp = (temperature - tgrid(t_ref_index))/(tgrid(t_ref_indexp1) - tgrid(t_ref_index))
    endif


!    write(*,*) "temperature", temperature
!    write(*,*) "t_ref_index", t_ref_index
!    write(*,*) "reference", tgrid(t_ref_index), tgrid(t_ref_indexp1)
!    write(*,*) "interp_temp", temp


    ! perform bilinear interpolation between P,T

    vbi(1) = kdata(igas, p_ref_index,   t_ref_index)
    vbi(2) = kdata(igas, p_ref_indexp1, t_ref_index)
    vbi(3) = kdata(igas, p_ref_index,   t_ref_indexp1)
    vbi(4) = kdata(igas, p_ref_indexp1, t_ref_indexp1)

!    write(*,*) "V", vbi(:)
 
    onemp = 1. - press
    onemt = 1. - temp
 
    ! bilinear  
    ans = vbi(1)*onemp*onemt &
        + vbi(2)*press*onemt &
        + vbi(3)*onemp*temp &
        + vbi(4)*press*temp
 !   write(*,*) "ans", ans
 !   write(*,*) "------------------------------------------------------" 

    return

  end subroutine bilinear_interpK_grey


!========================================================================

  subroutine interpH2Omtckd(Kdata, temperature, t_ref_index, ans) 

!------------------------------------------------------------------------
!
! Purpose: Linearly interpolate k-coefficient from reference P,T,W
!                                        
!------------------------------------------------------------------------
   
    implicit none
   
!------------------------------------------------------------------------
!
! Input Arguments
!
    real(r8), intent(in) :: Kdata(ntot_wavlnrng, kmtckd_ntemp)
    integer, intent(inout) :: t_ref_index
    real(r8), intent(in) :: temperature
    real(r8), intent(out) :: ans(ntot_wavlnrng)

!------------------------------------------------------------------------
!
! Local Variables
!
    integer :: p_ref_indexp1
    integer :: t_ref_indexp1
    integer :: w_ref_indexp1
    integer :: ik, iwi
   
    real(r8) :: u_col
    real(r8) :: press
    real(r8) :: temp
    real(r8) :: weight
    real(r8) :: onemp
    real(r8) :: onemt
    real(r8) :: onemw
!    real(r8), dimension(8) :: vtri
    
    logical :: interpp
    logical :: interpt

!------------------------------------------------------------------------
!
! Start Code
!

    t_ref_indexp1 = t_ref_index + 1

 !   vtri(:) = 0.
    ans(:) = 0.

    if (t_ref_index .eq. kmtckd_ntemp) then
      t_ref_index = t_ref_index - 1
      t_ref_indexp1 = t_ref_index + 1
      temp = (temperature - tgrid_mtckd(t_ref_index))/(tgrid_mtckd(t_ref_indexp1) - tgrid_mtckd(t_ref_index))
    else
      temp = (temperature - tgrid_mtckd(t_ref_index))/(tgrid_mtckd(t_ref_indexp1) - tgrid_mtckd(t_ref_index))
    endif

    do iwi = 1, ntot_wavlnrng
      ans(iwi) = kdata(iwi,t_ref_index)*(1.d0-temp)+kdata(iwi,t_ref_index+1)*temp 
    enddo

    return

  end subroutine interpH2Omtckd


!============================================================================

  subroutine interpolate_cld(cldgrp, iw_indx, Rcld, Qcld, Wcld, Gcld, & 
                             Qcldliq, Qcldice, Wcldliq, Wcldice, Gcldliq, Gcldice)

!------------------------------------------------------------------------
!
! Purpose: Interpolate cloud optical data to mode effective radius for ice
!          ice particles and liquid cloud drops.
!
!------------------------------------------------------------------------
 
    implicit none

!------------------------------------------------------------------------
!
! Arguments
! 
    integer, intent(in) :: cldgrp    ! cloud group index, 1 => liquid, 2 => ice
    integer, intent(in) :: iw_indx  ! spectral interval index
    real(r8), intent(in) :: Rcld      ! model output radius of cloud particle
    real(r8), intent(out) :: Qcld     ! interpolated value of mass extinction coefficient
    real(r8), intent(out) :: Wcld     ! interpolated value of single scattering albedo
    real(r8), intent(out) :: Gcld     ! interpolated value of asymmetry parameter
    real(r8), intent(in), dimension(nrel,ntot_wavlnrng) :: Qcldliq 
    real(r8), intent(in), dimension(nrei,ntot_wavlnrng) :: Qcldice
    real(r8), intent(in), dimension(nrel,ntot_wavlnrng) :: Wcldliq 
    real(r8), intent(in), dimension(nrei,ntot_wavlnrng) :: Wcldice
    real(r8), intent(in), dimension(nrel,ntot_wavlnrng) :: Gcldliq 
    real(r8), intent(in), dimension(nrei,ntot_wavlnrng) :: Gcldice
 
!------------------------------------------------------------------------
!
! Local Variables
! 

    integer :: Rcld_ref_index
    integer :: ir
    real(r8) :: fr

!------------------------------------------------------------------------
!
! Start Code
! 
    ir = int(Rcld)
    fr = dble(Rcld-real(ir))

    if (cldgrp .eq. 1) then   ! interpolate for liquid cloud droplets

      ! if Rcld less than minimum, force to be minimum grid value 
      if (Rcld .le. minval(rel_grid)) then 
        Qcld = Qcldliq(1,iw_indx)
        Wcld = Wcldliq(1,iw_indx)
        Gcld = gcldliq(1,iw_indx)
      endif

      ! if Rcld greater than maximum, force to be maximum grid value 
      if (Rcld .ge. maxval(rel_grid)) then 
        Qcld = Qcldliq(nrel,iw_indx)
        Wcld = Wcldliq(nrel,iw_indx)
        Gcld = gcldliq(nrel,iw_indx)
      endif

      ! interpolate Rcld to grid
      if (Rcld .gt. minval(rel_grid) .and. Rcld .lt. maxval(rel_grid)) then
        Qcld = Qcldliq(ir,iw_indx)*(1.d0-fr)+Qcldliq(ir+1,iw_indx)*fr 
        Wcld = Wcldliq(ir,iw_indx)*(1.d0-fr)+Wcldliq(ir+1,iw_indx)*fr 
        Gcld = Gcldliq(ir,iw_indx)*(1.d0-fr)+Gcldliq(ir+1,iw_indx)*fr 
      endif   
 
    endif

    if (cldgrp .eq. 2) then ! interpolate for ice cloud droplet
 
      ! if Rcld less than minimum, force to be minimum grid value 
      if (Rcld .le. minval(rei_grid)) then 
        Qcld = Qcldice(1,iw_indx)
        Wcld = Wcldice(1,iw_indx)
        Gcld = gcldice(1,iw_indx)
      endif

      ! if Rcld greater than maximum, force to be maximum grid value 
      if (Rcld .ge. maxval(rei_grid)) then 
        Qcld = Qcldice(nrei,iw_indx)
        Wcld = Wcldice(nrei,iw_indx)
        Gcld = gcldice(nrei,iw_indx)
      endif

      ! interpolate Rcld to grid
      if (Rcld .gt. minval(rei_grid) .and. Rcld .lt. maxval(rei_grid)) then
        Qcld = Qcldice(ir,iw_indx)*(1.d0-fr)+Qcldice(ir+1,iw_indx)*fr 
        Wcld = Wcldice(ir,iw_indx)*(1.d0-fr)+Wcldice(ir+1,iw_indx)*fr 
        Gcld = Gcldice(ir,iw_indx)*(1.d0-fr)+Gcldice(ir+1,iw_indx)*fr  
      endif  

    endif
   
  end subroutine interpolate_cld



!============================================================================

  subroutine interpH2H2cia(Kdata, temperature, t_ref_index, ans )

!------------------------------------------------------------------------
!
! Purpose: Interpolate H2-H2 collision induced absorption coefficients
!
!------------------------------------------------------------------------
 
    implicit none

!------------------------------------------------------------------------
!
! Arguments
! 
    real(r8), intent(in) :: Kdata(ntot_wavlnrng,kh2h2_ntemp)
    real(r8), intent(in) :: temperature
    integer, intent(inout) :: t_ref_index
    real(r8), intent(out) :: ans(ntot_wavlnrng)  

!------------------------------------------------------------------------
!
! Local Variables
!

    integer :: t_ref_indexp1
    integer :: iwi
    real(r8) :: temp
    real(r8), dimension(2) :: vli
    real(r8) :: ydiff

!------------------------------------------------------------------------
!
! Start Code
!
    !write(*,*) "-------- interpH2H2cia ----------"
    t_ref_indexp1 = t_ref_index + 1
    ans = 0.

    if (t_ref_index .eq. kh2h2_ntemp) then
      t_ref_index = t_ref_index - 1
      t_ref_indexp1 = t_ref_index + 1
      temp = (temperature - tgrid_h2h2(t_ref_index))/(tgrid_h2h2(t_ref_indexp1) - tgrid_h2h2(t_ref_index))
    else
      temp = (temperature - tgrid_h2h2(t_ref_index))/(tgrid_h2h2(t_ref_indexp1) - tgrid_h2h2(t_ref_index))
    endif

    !write(*,*) "temperature", temperature
    !write(*,*) "t_ref_index", t_ref_index
    !write(*,*) "reference", tgrid(t_ref_index),tgrid(t_ref_indexp1)
    !write(*,*) "interp_temp", temp

    do iwi=1, ntot_wavlnrng
      ! linear interpolation T
      vli(1) = kdata(iwi,   t_ref_index)
      vli(2) = kdata(iwi,   t_ref_indexp1)
 
      ydiff = vli(2) - vli(1)
      ans(iwi) = vli(1) + ydiff*temp 
    enddo
    !write(*,*) "ans", ans
    !write(*,*) "------------------------------------------------------" 

    return

  end subroutine interpH2H2cia



!============================================================================

  subroutine interpN2H2cia(Kdata, temperature, t_ref_index, ans )

!------------------------------------------------------------------------
!
! Purpose: Interpolate N2-H2 collision induced absorption coefficients
!
!------------------------------------------------------------------------
 
    implicit none

!------------------------------------------------------------------------
!
! Arguments
! 
    real(r8), intent(in) :: Kdata(ntot_wavlnrng,kn2h2_ntemp)
    real(r8), intent(in) :: temperature
    integer, intent(inout) :: t_ref_index
    real(r8), intent(out) :: ans(ntot_wavlnrng)  

!------------------------------------------------------------------------
!
! Local Variables
!

    integer :: t_ref_indexp1
    integer :: iwi
    real(r8) :: temp
    real(r8), dimension(2) :: vli
    real(r8) :: ydiff

!------------------------------------------------------------------------
!
! Start Code
!
    !write(*,*) "-------- interpN2H2cia ----------"
    t_ref_indexp1 = t_ref_index + 1
    ans(:) = 0.

    if (t_ref_index .eq. kn2h2_ntemp) then
      t_ref_index = t_ref_index - 1
      t_ref_indexp1 = t_ref_index + 1
      temp = (temperature - tgrid_n2h2(t_ref_index))/(tgrid_n2h2(t_ref_indexp1) - tgrid_n2h2(t_ref_index))
    else
      temp = (temperature - tgrid_n2h2(t_ref_index))/(tgrid_n2h2(t_ref_indexp1) - tgrid_n2h2(t_ref_index))
    endif

    !write(*,*) "temperature", temperature
    !write(*,*) "t_ref_index", t_ref_index
    !write(*,*) "reference", tgrid(t_ref_index),tgrid(t_ref_indexp1)
    !write(*,*) "interp_temp", temp

    do iwi=1, ntot_wavlnrng
      ! linear interpolation T
      vli(1) = kdata(iwi,   t_ref_index)
      vli(2) = kdata(iwi,   t_ref_indexp1)
      ydiff = vli(2) - vli(1)
      ans(iwi) = vli(1) + ydiff*temp 
    enddo
    !write(*,*) "ans", ans
    !write(*,*) "------------------------------------------------------" 

    return

  end subroutine interpN2H2cia


!============================================================================

  subroutine interpN2N2cia(Kdata, temperature, t_ref_index, ans )


!------------------------------------------------------------------------
!
! Purpose: Interpolate N2-N2 collision induced absorption coefficients
!
!------------------------------------------------------------------------
 
    implicit none

!------------------------------------------------------------------------
!
! Arguments
! 
    real(r8), intent(in) :: Kdata(ntot_wavlnrng,kn2n2_ntemp)
    real(r8), intent(in) :: temperature
    integer, intent(inout) :: t_ref_index
    real(r8), intent(out) :: ans(ntot_wavlnrng)  

!------------------------------------------------------------------------
!
! Local Variables
!

    integer :: t_ref_indexp1
    integer :: iwi
    real(r8) :: temp
    real(r8), dimension(2) :: vli
    real(r8) :: ydiff

!------------------------------------------------------------------------
!
! Start Code
!
    !write(*,*) "-------- interpN2N2cia ----------"
    t_ref_indexp1 = t_ref_index + 1
    ans(:) = 0.

    if (t_ref_index .eq. kn2n2_ntemp) then
      t_ref_index = t_ref_index - 1
      t_ref_indexp1 = t_ref_index + 1
      temp = (temperature - tgrid_n2n2(t_ref_index))/(tgrid_n2n2(t_ref_indexp1) - tgrid_n2n2(t_ref_index))
    else
      temp = (temperature - tgrid_n2n2(t_ref_index))/(tgrid_n2n2(t_ref_indexp1) - tgrid_n2n2(t_ref_index))
    endif

    !write(*,*) "temperature", temperature
    !write(*,*) "t_ref_index", t_ref_index
    !write(*,*) "reference", tgrid(t_ref_index),tgrid(t_ref_indexp1)
    !write(*,*) "interp_temp", temp

    do iwi=1, ntot_wavlnrng
      ! linear interpolation T
      vli(1) = kdata(iwi,   t_ref_index)
      vli(2) = kdata(iwi,   t_ref_indexp1)
 
      ydiff = vli(2) - vli(1)
      ans(iwi) = vli(1) + ydiff*temp 
    enddo
    !write(*,*) "ans", ans
    !write(*,*) "------------------------------------------------------" 

    return

  end subroutine interpN2N2cia



!============================================================================

  subroutine interpCO2CO2cia_sw(Kdata, temperature, t_ref_index, ans )

!------------------------------------------------------------------------
!
! Purpose: Interpolate CO2-CO2 SW collision induced absorption coefficients
!
!------------------------------------------------------------------------
 
    implicit none

!------------------------------------------------------------------------
!
! Arguments
! 
    real(r8), intent(in) :: Kdata(ntot_wavlnrng,kco2co2_sw_ntemp)
    real(r8), intent(in) :: temperature
    integer, intent(inout) :: t_ref_index
    real(r8), intent(out) :: ans(ntot_wavlnrng)  

!------------------------------------------------------------------------
!
! Local Variables
!

    integer :: t_ref_indexp1
    integer :: iwi
    real(r8) :: temp
    real(r8), dimension(2) :: vli
    real(r8) :: ydiff

!------------------------------------------------------------------------
!
! Start Code
!
    !write(*,*) "-------- interpCO2CO2cia_sw ----------"
    t_ref_indexp1 = t_ref_index + 1
    ans(:) = 0.

    if (t_ref_index .eq. kco2co2_sw_ntemp) then
      t_ref_index = t_ref_index - 1
      t_ref_indexp1 = t_ref_index + 1
      temp = (temperature - tgrid_co2co2_sw(t_ref_index))/(tgrid_co2co2_sw(t_ref_indexp1) - tgrid_co2co2_sw(t_ref_index))
    else
      temp = (temperature - tgrid_co2co2_sw(t_ref_index))/(tgrid_co2co2_sw(t_ref_indexp1) - tgrid_co2co2_sw(t_ref_index))
    endif

    !write(*,*) "temperature", temperature
    !write(*,*) "t_ref_index", t_ref_index
    !write(*,*) "reference", tgrid(t_ref_index),tgrid(t_ref_indexp1)
    !write(*,*) "interp_temp", temp

    do iwi=1, ntot_wavlnrng
      ! linear interpolation T
      vli(1) = kdata(iwi,   t_ref_index)
      vli(2) = kdata(iwi,   t_ref_indexp1)
 
      ydiff = vli(2) - vli(1)
      ans(iwi) = vli(1) + ydiff*temp 
    enddo
    !write(*,*) "ans", ans
    !write(*,*) "------------------------------------------------------" 

    return

  end subroutine interpCO2CO2cia_sw


!============================================================================

  subroutine interpCO2CO2cia_lw(Kdata, temperature, t_ref_index, ans )

!------------------------------------------------------------------------
!
! Purpose: Interpolate CO2-CO2 LW collision induced absorption coefficients
!
!------------------------------------------------------------------------
 
    implicit none

!------------------------------------------------------------------------
!
! Arguments
! 
    real(r8), intent(in) :: Kdata(ntot_wavlnrng,kco2co2_lw_ntemp)
    real(r8), intent(in) :: temperature
    integer, intent(inout) :: t_ref_index
    real(r8), intent(out) :: ans(ntot_wavlnrng)  

!------------------------------------------------------------------------
!
! Local Variables
!

    integer :: t_ref_indexp1
    integer :: iwi
    real(r8) :: temp
    real(r8), dimension(2) :: vli
    real(r8) :: ydiff

!------------------------------------------------------------------------
!
! Start Code
!
    !write(*,*) "-------- interpCO2CO2cia_lw ----------"
    t_ref_indexp1 = t_ref_index + 1
    ans(:) = 0.

    if (t_ref_index .eq. kco2co2_lw_ntemp) then
      t_ref_index = t_ref_index - 1
      t_ref_indexp1 = t_ref_index + 1
      temp = (temperature - tgrid_co2co2_lw(t_ref_index))/(tgrid_co2co2_lw(t_ref_indexp1) - tgrid_co2co2_lw(t_ref_index))
    else
      temp = (temperature - tgrid_co2co2_lw(t_ref_index))/(tgrid_co2co2_lw(t_ref_indexp1) - tgrid_co2co2_lw(t_ref_index))
    endif

    !write(*,*) "temperature", temperature
    !write(*,*) "t_ref_index", t_ref_index
    !write(*,*) "reference", tgrid(t_ref_index),tgrid(t_ref_indexp1)
    !write(*,*) "interp_temp", temp

    do iwi=1, ntot_wavlnrng
      ! linear interpolation T
      vli(1) = kdata(iwi,   t_ref_index)
      vli(2) = kdata(iwi,   t_ref_indexp1)
 
      ydiff = vli(2) - vli(1)
      ans(iwi) = vli(1) + ydiff*temp 
    enddo
    !write(*,*) "ans", ans
    !write(*,*) "------------------------------------------------------" 

    return

  end subroutine interpCO2CO2cia_lw


!============================================================================

  subroutine interpCO2H2cia(Kdata, temperature, t_ref_index, ans )

!------------------------------------------------------------------------
!
! Purpose: Interpolate CO2-H2 collision induced absorption coefficients
!
!------------------------------------------------------------------------
 
    implicit none

!------------------------------------------------------------------------
!
! Arguments
! 
    real(r8), intent(in) :: Kdata(ntot_wavlnrng,kco2h2_ntemp)
    real(r8), intent(in) :: temperature
    integer, intent(inout) :: t_ref_index
    real(r8), intent(out) :: ans(ntot_wavlnrng)  

!------------------------------------------------------------------------
!
! Local Variables
!

    integer :: t_ref_indexp1
    integer :: iwi
    real(r8) :: temp
    real(r8), dimension(2) :: vli
    real(r8) :: ydiff

!------------------------------------------------------------------------
!
! Start Code
!
    !write(*,*) "-------- interpCO2H2cia ----------"
    t_ref_indexp1 = t_ref_index + 1
    ans(:) = 0.

    if (t_ref_index .eq. kco2h2_ntemp) then
      t_ref_index = t_ref_index - 1
      t_ref_indexp1 = t_ref_index + 1
      temp = (temperature - tgrid_co2h2(t_ref_index))/(tgrid_co2h2(t_ref_indexp1) - tgrid_co2h2(t_ref_index))
    else
      temp = (temperature - tgrid_co2h2(t_ref_index))/(tgrid_co2h2(t_ref_indexp1) - tgrid_co2h2(t_ref_index))
    endif

    !write(*,*) "temperature", temperature
    !write(*,*) "t_ref_index", t_ref_index
    !write(*,*) "reference", tgrid(t_ref_index),tgrid(t_ref_indexp1)
    !write(*,*) "interp_temp", temp

    do iwi=1, ntot_wavlnrng
      ! linear interpolation T
      vli(1) = kdata(iwi,   t_ref_index)
      vli(2) = kdata(iwi,   t_ref_indexp1)
 
      ydiff = vli(2) - vli(1)
      ans(iwi) = vli(1) + ydiff*temp 
    enddo
    !write(*,*) "ans", ans
    !write(*,*) "------------------------------------------------------" 

    return

  end subroutine interpCO2H2cia


!============================================================================

  subroutine interpCO2CH4cia(Kdata, temperature, t_ref_index, ans )

!------------------------------------------------------------------------
!
! Purpose: Interpolate CO2-CH4 collision induced absorption coefficients
!
!------------------------------------------------------------------------
 
    implicit none

!------------------------------------------------------------------------
!
! Arguments
! 
    real(r8), intent(in) :: Kdata(ntot_wavlnrng,kco2ch4_ntemp)
    real(r8), intent(in) :: temperature
    integer, intent(inout) :: t_ref_index
    real(r8), intent(out) :: ans(ntot_wavlnrng)  

!------------------------------------------------------------------------
!
! Local Variables
!

    integer :: t_ref_indexp1
    integer :: iwi
    real(r8) :: temp
    real(r8), dimension(2) :: vli
    real(r8) :: ydiff

!------------------------------------------------------------------------
!
! Start Code
!
    !write(*,*) "-------- interpCO2CH4cia ----------"
    t_ref_indexp1 = t_ref_index + 1
    ans(:) = 0.

    if (t_ref_index .eq. kco2ch4_ntemp) then
      t_ref_index = t_ref_index - 1
      t_ref_indexp1 = t_ref_index + 1
      temp = (temperature - tgrid_co2ch4(t_ref_index))/(tgrid_co2ch4(t_ref_indexp1) - tgrid_co2ch4(t_ref_index))
    else
      temp = (temperature - tgrid_co2ch4(t_ref_index))/(tgrid_co2ch4(t_ref_indexp1) - tgrid_co2ch4(t_ref_index))
    endif

    !write(*,*) "temperature", temperature
    !write(*,*) "t_ref_index", t_ref_index
    !write(*,*) "reference", tgrid(t_ref_index),tgrid(t_ref_indexp1)
    !write(*,*) "interp_temp", temp

    do iwi=1, ntot_wavlnrng
      ! linear interpolation T
      vli(1) = kdata(iwi,   t_ref_index)
      vli(2) = kdata(iwi,   t_ref_indexp1)
 
      ydiff = vli(2) - vli(1)
      ans(iwi) = vli(1) + ydiff*temp 
    enddo
    !write(*,*) "ans", ans
    !write(*,*) "------------------------------------------------------" 

    return

  end subroutine interpCO2CH4cia



end module rad_interp_mod
