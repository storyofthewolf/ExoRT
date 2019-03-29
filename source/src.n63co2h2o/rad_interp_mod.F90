#include <params.h>
#include <misc.h>

module rad_interp_mod

! version highco2
! Contains subroutines used for interpolating K coefficient and 
! cloud optical data in radiation.F90
! this can probably be moved to src.main

  use shr_kind_mod,   only: r8 => shr_kind_r8
  use radgrid

  implicit none
  private
  save

!------------------------------------------------------------------------
!
!  Public Interfaces
!

  public :: trilinear_interpK_16gpt
  public :: interph2o_cont
  public :: interpolate_cld

!========================================================================
contains
!========================================================================

  subroutine interph2o_cont(Kdata, iw, temperature, t_ref_index, ans) 

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
    real(r8), intent(in) :: Kdata(ntot_wavlnrng, ks_ntemp)
    integer, intent(in) :: iw
    integer, intent(inout) :: t_ref_index
    real(r8), intent(in) :: temperature
    real(r8), intent(out) :: ans

!------------------------------------------------------------------------
!
! Local Variables
!
    integer :: p_ref_indexp1
    integer :: t_ref_indexp1
    integer :: w_ref_indexp1
    integer :: ik
   
    real(r8) :: u_col
    real(r8) :: press
    real(r8) :: temp
    real(r8) :: weight
    real(r8) :: onemp
    real(r8) :: onemt
    real(r8) :: onemw
    real(r8), dimension(8) :: vtri
    
    logical :: interpp
    logical :: interpt

!------------------------------------------------------------------------
!
! Start Code
!

    t_ref_indexp1 = t_ref_index + 1

    vtri(:) = 0.
    ans = 0.

    if (t_ref_index .eq. ks_ntemp) then
      t_ref_index = t_ref_index - 1
      t_ref_indexp1 = t_ref_index + 1
      temp = (temperature - tgrid_self(t_ref_index))/(tgrid_self(t_ref_indexp1) - tgrid_self(t_ref_index))
    else
      temp = (temperature - tgrid_self(t_ref_index))/(tgrid_self(t_ref_indexp1) - tgrid_self(t_ref_index))
    endif


    ans = kdata(iw,t_ref_index)*(1.d0-temp)+kdata(iw,t_ref_index+1)*temp 
!! Experimental
!!EW Klugde, return answer at bottom temperaturem
!    ans = kdata(iw,t_ref_index+1)

!    write(*,*) "--------------------------------------------------------------"
!    write(*,*) temperature, tgrid_self(t_ref_index), tgrid_self(t_ref_indexp1)
!    write(*,*) kdata(iw,t_ref_index)*(1.d0-temp)+kdata(iw,t_ref_index+1)*temp 
!    write(*,*) "interpKself2", kdata(iw,t_ref_index), kdata(iw,t_ref_index+1)
!    write(*,*) "interpKself2", tgrid_self(t_ref_index)/temperature, tgrid_self(t_ref_indexp1)/temperature
!    write(*,*) "interpKself2", kdata(iw,t_ref_index)*tgrid_self(t_ref_index)/temperature, kdata(iw,t_ref_index+1)*tgrid_self(t_ref_indexp1)/temperature

    return

  end subroutine interph2o_cont


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

  subroutine trilinear_interpK_16gpt(Kdata, ig, pressure, p_ref_index, temperature, t_ref_index, &
                           species_vmr, w_ref_index, ans) 

!------------------------------------------------------------------------
!
! Purpose: Linearly interpolate k-coefficient from reference P,T,W
! for co2_h2ovar atmosphere pressure grid                                       
!
!------------------------------------------------------------------------
   
    implicit none
   
!------------------------------------------------------------------------
!
! Input Arguments
!
    real(r8), intent(in) :: Kdata(16, kc_nvar, kc_npress, kc_ntemp)
    integer, intent(in) :: ig
    integer, intent(inout) :: p_ref_index
    integer, intent(inout) :: t_ref_index
    integer, intent(inout) :: w_ref_index

    real(r8), intent(in) :: pressure  
    real(r8), intent(in) :: temperature
    real(r8), intent(in) :: species_vmr
    real(r8), intent(out) :: ans

!------------------------------------------------------------------------
!
! Local Variables
!
    integer :: p_ref_indexp1
    integer :: t_ref_indexp1
    integer :: w_ref_indexp1
    integer :: ik
   
    real(r8) :: u_col
    real(r8) :: press
    real(r8) :: temp
    real(r8) :: weight
    real(r8) :: onemp
    real(r8) :: onemt
    real(r8) :: onemw
    real(r8), dimension(8) :: vtri
    
    logical :: interpp
    logical :: interpt

!------------------------------------------------------------------------
!
! Start Code
!
    !write(*,*) "-------- trilinear_interpK ----------"
    p_ref_indexp1 = p_ref_index + 1    
    t_ref_indexp1 = t_ref_index + 1
    w_ref_indexp1 = w_ref_index + 1

    vtri(:) = 0.
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
   
 !   write(*,*) "pressure", pressure
 !   write(*,*) "p_ref_index", p_ref_index
 !   write(*,*) "reference", log10pgrid(p_ref_index),log10pgrid(p_ref_indexp1)
 !   write(*,*) "interp_press", press
  
    if (t_ref_index .eq. kc_ntemp) then
      t_ref_index = t_ref_index - 1
      t_ref_indexp1 = t_ref_index + 1
      temp = (temperature - tgrid(t_ref_index))/(tgrid(t_ref_indexp1) - tgrid(t_ref_index))
    else
      temp = (temperature - tgrid(t_ref_index))/(tgrid(t_ref_indexp1) - tgrid(t_ref_index))
    endif

 !   write(*,*) "temperature", temperature
 !   write(*,*) "t_ref_index", t_ref_index
 !   write(*,*) "reference", tgrid(t_ref_index),tgrid(t_ref_indexp1)
 !   write(*,*) "interp_temp", temp


    ! find the interpolating factor for the various parameters.  at the tops of grids, extrapolate
    if (w_ref_index .eq. kc_nvar) then
      w_ref_index = w_ref_index - 1
      w_ref_indexp1 = w_ref_index + 1
      weight = (log10(species_vmr) - log10wgrid(w_ref_index)) / &
              (log10wgrid(w_ref_indexp1) - log10wgrid(w_ref_index))           
    else
      weight = (log10(species_vmr) - log10wgrid(w_ref_index)) / &
              (log10wgrid(w_ref_indexp1) - log10wgrid(p_ref_index))  
    endif

    if (w_ref_index .eq. kc_nvar) then
      w_ref_index = w_ref_index - 1
      w_ref_indexp1 = w_ref_index + 1
      weight = (species_vmr - wgrid(w_ref_index))/(wgrid(w_ref_indexp1) - wgrid(w_ref_index))
    else
      weight = (species_vmr - wgrid(w_ref_index))/(wgrid(w_ref_indexp1) - wgrid(w_ref_index))
    endif

!    write(*,*) "---------------------------------------------"
!    write(*,*) "species_vmr", species_vmr, log10(species_vmr)
!    write(*,*) "w_ref_index", w_ref_index
!    write(*,*) "reference", wgrid(w_ref_index),wgrid(w_ref_indexp1)
!    write(*,*) "reference", log10wgrid(w_ref_index),log10wgrid(w_ref_indexp1)
!    write(*,*) "interp_weight", weight
! 
    ! perform trilinear interpolation between P,T,W

!    vtri(1) = kdata(ig, p_ref_index,   t_ref_index,   w_ref_index)
!    vtri(2) = kdata(ig, p_ref_indexp1, t_ref_index,   w_ref_index)
!    vtri(3) = kdata(ig, p_ref_index,   t_ref_indexp1, w_ref_index)
!    vtri(4) = kdata(ig, p_ref_index,   t_ref_index,   w_ref_indexp1)
!    vtri(5) = kdata(ig, p_ref_indexp1, t_ref_index,   w_ref_indexp1)
!    vtri(6) = kdata(ig, p_ref_index,   t_ref_indexp1, w_ref_indexp1)
!    vtri(7) = kdata(ig, p_ref_indexp1, t_ref_indexp1, w_ref_index)
!    vtri(8) = kdata(ig, p_ref_indexp1, t_ref_indexp1, w_ref_indexp1)

    vtri(1) = kdata(ig, w_ref_index,   p_ref_index,   t_ref_index)
    vtri(2) = kdata(ig, w_ref_indexp1, p_ref_index,   t_ref_index)
    vtri(3) = kdata(ig, w_ref_index,   p_ref_indexp1, t_ref_index)
    vtri(4) = kdata(ig, w_ref_index,   p_ref_index,   t_ref_indexp1)
    vtri(5) = kdata(ig, w_ref_indexp1, p_ref_index,   t_ref_indexp1)
    vtri(6) = kdata(ig, w_ref_index,   p_ref_indexp1, t_ref_indexp1)
    vtri(7) = kdata(ig, w_ref_indexp1, p_ref_indexp1, t_ref_index)
    vtri(8) = kdata(ig, w_ref_indexp1, p_ref_indexp1, t_ref_indexp1)

  !  write(*,*) "V", vtri(:)
 
    onemp = 1. - press
    onemt = 1. - temp
    onemw = 1. - weight

    ans = vtri(1)*onemp*onemt*onemw &
        + vtri(2)*press*onemt*onemw &
        + vtri(3)*onemp*temp*onemw &
        + vtri(4)*onemp*onemt*weight &
        + vtri(5)*press*onemt*weight &
        + vtri(6)*onemp*temp*weight &
        + vtri(7)*press*temp*onemw &
        + vtri(8)*press*temp*weight

  !  write(*,*) "ans", ans

    return

  end subroutine trilinear_interpK_16gpt



end module rad_interp_mod
