! Minimal CESM/CAM module stubs — compile-check harness for
! 3dmodels/src.cam.exort. Interfaces mirror CESM1.2.1 usage in the bundle;
! bodies are empty. NOT for execution.

module abortutils
  implicit none
  public
contains
  subroutine endrun(msg)
    character(len=*), intent(in), optional :: msg
    stop 1
  end subroutine endrun
end module abortutils

module time_manager
  use shr_kind_mod, only: r8 => shr_kind_r8
  implicit none
  public
contains
  integer function get_nstep()
    get_nstep = 0
  end function get_nstep
  integer function get_step_size()
    get_step_size = 1800
  end function get_step_size
  real(r8) function get_curr_calday(offset)
    integer, intent(in), optional :: offset
    get_curr_calday = 1.0_r8
  end function get_curr_calday
  subroutine get_curr_calday_rotation(frac_day, day_in_year)
    real(r8), intent(out) :: frac_day, day_in_year
    frac_day = 0.0_r8; day_in_year = 1.0_r8
  end subroutine get_curr_calday_rotation
end module time_manager

module pio
  use shr_kind_mod, only: r8 => shr_kind_r8
  implicit none
  public
  integer, parameter :: pio_nowrite = 0
  type file_desc_t
    integer :: fh = -1
  end type file_desc_t
  type var_desc_t
    integer :: varid = -1
  end type var_desc_t
contains
  integer function pio_inq_varid(ncid, name, varid)
    type(file_desc_t), intent(in) :: ncid
    character(len=*), intent(in) :: name
    integer, intent(out) :: varid
    varid = 1; pio_inq_varid = 0
  end function pio_inq_varid
  integer function pio_inq_dimid(ncid, name, dimid)
    type(file_desc_t), intent(in) :: ncid
    character(len=*), intent(in) :: name
    integer, intent(out) :: dimid
    dimid = 1; pio_inq_dimid = 0
  end function pio_inq_dimid
  integer function pio_inquire_dimension(ncid, dimid, len)
    type(file_desc_t), intent(in) :: ncid
    integer, intent(in) :: dimid
    integer, intent(out), optional :: len
    if (present(len)) len = 0
    pio_inquire_dimension = 0
  end function pio_inquire_dimension
  integer function pio_get_var(ncid, varid, values)
    type(file_desc_t), intent(in) :: ncid
    integer, intent(in) :: varid
    real(r8), intent(out), dimension(..) :: values
    pio_get_var = 0
  end function pio_get_var
  subroutine pio_closefile(ncid)
    type(file_desc_t), intent(inout) :: ncid
  end subroutine pio_closefile
end module pio

module cam_pio_utils
  use pio, only: file_desc_t
  implicit none
  public
contains
  subroutine cam_pio_openfile(file, fname, mode)
    type(file_desc_t), intent(inout) :: file
    character(len=*), intent(in) :: fname
    integer, intent(in) :: mode
  end subroutine cam_pio_openfile
end module cam_pio_utils

module cam_history
  use shr_kind_mod, only: r8 => shr_kind_r8
  implicit none
  public
  integer, parameter :: phys_decomp = 1
contains
  subroutine addfld(fname, units, numlev, avgflag, long_name, decomp_type, &
                    flag_xyfill, flag_isccplev, sampling_seq)
    character(len=*), intent(in) :: fname, units, long_name
    integer, intent(in) :: numlev, decomp_type
    character(len=1), intent(in) :: avgflag
    logical, intent(in), optional :: flag_xyfill, flag_isccplev
    character(len=*), intent(in), optional :: sampling_seq
  end subroutine addfld
  subroutine add_default(name, tindex, flag)
    character(len=*), intent(in) :: name
    integer, intent(in) :: tindex
    character(len=1), intent(in) :: flag
  end subroutine add_default
  subroutine outfld(fname, field, idim, c)
    character(len=*), intent(in) :: fname
    real(r8), intent(in), dimension(..) :: field
    integer, intent(in) :: idim, c
  end subroutine outfld
end module cam_history

module physics_buffer
  use shr_kind_mod, only: r8 => shr_kind_r8
  implicit none
  public
  type physics_buffer_desc
    integer :: dummy = 0
  end type physics_buffer_desc
contains
  integer function pbuf_get_index(name)
    character(len=*), intent(in) :: name
    pbuf_get_index = 1
  end function pbuf_get_index
  integer function pbuf_old_tim_idx()
    pbuf_old_tim_idx = 1
  end function pbuf_old_tim_idx
  subroutine pbuf_get_field(pbuf, index, field, start, kount)
    type(physics_buffer_desc), pointer :: pbuf(:)
    integer, intent(in) :: index
    real(r8), pointer :: field(:,:)
    integer, intent(in), optional :: start(:), kount(:)
    nullify(field)
  end subroutine pbuf_get_field
end module physics_buffer

module physics_types
  use shr_kind_mod, only: r8 => shr_kind_r8
  use ppgrid,       only: pcols, pver, pverp
  implicit none
  public
  type physics_state
    integer :: lchnk = 1
    integer :: ncol = pcols
    real(r8), dimension(pcols) :: ps
    real(r8), dimension(pcols,pver) :: pmid, pdel, pdeldry, t
    real(r8), dimension(pcols,pverp) :: pint, pintdry, zi
    real(r8), dimension(pcols,pver,1) :: q
  end type physics_state
  type physics_ptend
    integer :: dummy = 0
  end type physics_ptend
end module physics_types

module camsrfexch
  use shr_kind_mod, only: r8 => shr_kind_r8
  use ppgrid,       only: pcols
  implicit none
  public
  type cam_out_t
    real(r8), dimension(pcols) :: sols, soll, solsd, solld, flwds
  end type cam_out_t
  type cam_in_t
    real(r8), dimension(pcols) :: ts, asdir, aldir, asdif, aldif
    real(r8), dimension(pcols) :: srf_emiss
  end type cam_in_t
end module camsrfexch

module phys_grid
  use shr_kind_mod, only: r8 => shr_kind_r8
  implicit none
  public
contains
  subroutine get_rlat_all_p(lchnk, ncol, clat)
    integer, intent(in) :: lchnk, ncol
    real(r8), intent(out) :: clat(:)
    clat = 0.0_r8
  end subroutine get_rlat_all_p
  subroutine get_rlon_all_p(lchnk, ncol, clon)
    integer, intent(in) :: lchnk, ncol
    real(r8), intent(out) :: clon(:)
    clon = 0.0_r8
  end subroutine get_rlon_all_p
end module phys_grid

module radheat
  use shr_kind_mod,   only: r8 => shr_kind_r8
  use physics_types,  only: physics_state, physics_ptend
  use physics_buffer, only: physics_buffer_desc
  implicit none
  public
contains
  subroutine radheat_tend(state, pbuf, ptend, qrl, qrs, fsns, fsnt, flns, flnt, asdir, net_flx)
    type(physics_state), intent(in) :: state
    type(physics_buffer_desc), pointer :: pbuf(:)
    type(physics_ptend), intent(out) :: ptend
    real(r8), intent(in) :: qrl(:,:), qrs(:,:)
    real(r8), intent(in) :: fsns(:), fsnt(:), flns(:), flnt(:), asdir(:)
    real(r8), intent(out) :: net_flx(:)
    net_flx = 0.0_r8
  end subroutine radheat_tend
end module radheat

module pspect
  implicit none
  public
end module pspect

module rad_constituents
  use shr_kind_mod,   only: r8 => shr_kind_r8
  use physics_types,  only: physics_state
  use physics_buffer, only: physics_buffer_desc
  implicit none
  public
contains
  subroutine rad_cnst_get_gas(list_idx, gasname, state, pbuf, mmr)
    integer, intent(in) :: list_idx
    character(len=*), intent(in) :: gasname
    type(physics_state), intent(in), target :: state
    type(physics_buffer_desc), pointer :: pbuf(:)
    real(r8), pointer :: mmr(:,:)
    nullify(mmr)
  end subroutine rad_cnst_get_gas
  subroutine rad_cnst_out()
  end subroutine rad_cnst_out
end module rad_constituents

module shr_orb_mod
  use shr_kind_mod, only: r8 => shr_kind_r8
  implicit none
  public
contains
  subroutine shr_orb_decl(calday, eccen, mvelpp, lambm0, obliqr, delta, eccf)
    real(r8), intent(in) :: calday, eccen, mvelpp, lambm0, obliqr
    real(r8), intent(out) :: delta, eccf
    delta = 0.0_r8; eccf = 1.0_r8
  end subroutine shr_orb_decl
end module shr_orb_mod

module cloud_cover_diags
  use shr_kind_mod, only: r8 => shr_kind_r8
  implicit none
  public
contains
  subroutine cloud_cover_diags_out(lchnk, ncol, cld, pmid, nmxrgn, pmxrgn)
    integer, intent(in) :: lchnk, ncol
    real(r8), intent(in) :: cld(:,:), pmid(:,:)
    integer, intent(in) :: nmxrgn(:)
    real(r8), intent(in) :: pmxrgn(:,:)
  end subroutine cloud_cover_diags_out
end module cloud_cover_diags

module cam_control_mod
  use shr_kind_mod, only: r8 => shr_kind_r8
  implicit none
  public
  real(r8) :: lambm0 = 0.0_r8
  real(r8) :: obliqr = 0.0_r8
  real(r8) :: eccen  = 0.0_r8
  real(r8) :: mvelpp = 0.0_r8
end module cam_control_mod

module radiation_data
  use shr_kind_mod,   only: r8 => shr_kind_r8
  use physics_types,  only: physics_state
  use physics_buffer, only: physics_buffer_desc
  use camsrfexch,     only: cam_in_t
  implicit none
  public
contains
  subroutine init_rad_data()
  end subroutine init_rad_data
  subroutine output_rad_data(pbuf, state, cam_in, landm, coszrs)
    type(physics_buffer_desc), pointer :: pbuf(:)
    type(physics_state), intent(in) :: state
    type(cam_in_t), intent(in) :: cam_in
    real(r8), intent(in) :: landm(:)
    real(r8), intent(in) :: coszrs
  end subroutine output_rad_data
end module radiation_data

module phys_control
  implicit none
  public
contains
  subroutine phys_getopts(microp_scheme_out)
    character(len=*), intent(out), optional :: microp_scheme_out
    if (present(microp_scheme_out)) microp_scheme_out = 'RK'
  end subroutine phys_getopts
end module phys_control

module carma_model_mod
  implicit none
  public
  integer, parameter :: NELEM = 1
  integer, parameter :: NBIN  = 40
end module carma_model_mod

module carma_exort_mod
  use shr_kind_mod,    only: r8 => shr_kind_r8
  use ppgrid,          only: pcols, pver
  use carma_model_mod, only: NELEM, NBIN
  use physics_types,   only: physics_state
  implicit none
  public
contains
  subroutine carma_exort_get_mmr(state, carmamix)
    type(physics_state), intent(in) :: state
    real(r8), intent(out), dimension(pcols,pver,NELEM,NBIN) :: carmamix
    carmamix = 0.0_r8
  end subroutine carma_exort_get_mmr
end module carma_exort_mod

! bare external referenced by exo_radiation_cam_intr (lives in ExoCAM's
! SourceMods zenith.F90 in a real build)
subroutine zenith_rotation(frac_day, calday, clat, clon, coszrs, ncol)
  use shr_kind_mod, only: r8 => shr_kind_r8
  implicit none
  real(r8), intent(in) :: frac_day, calday
  real(r8), intent(in) :: clat(:), clon(:)
  real(r8), intent(out) :: coszrs(:)
  integer, intent(in) :: ncol
  coszrs = 0.5_r8
end subroutine zenith_rotation
