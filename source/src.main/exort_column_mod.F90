
module exort_column_mod
!----------------------------------------------------------------------
! C-interoperable column state / result types for the ExoRT library
! API (Stage D). One column_state_t carries everything aerad_driver
! needs for a single column; one column_result_t carries everything it
! returns. Both are bind(c) so the identical struct layout is visible
! from C and Python (see tools/exort_pytools/).
!
! Array extents are compile-time (exo_pver in exoplanet_mod.F90, band
! and CARMA dims in radgrid.F90); callers query them at runtime via
! exort_get_dims() rather than hardcoding.
!
! Field order is the ABI contract with tools/exort_pytools/exort_api.py
! and tests/lib/test_exort_c.c — change all three together.
!----------------------------------------------------------------------

use iso_c_binding, only: c_double
use ppgrid,        only: pver, pverp
use radgrid,       only: ntot_wavlnrng, nelem_carma, nbin_carma

implicit none
public

type, bind(c) :: column_state_t
  ! scalars
  real(c_double) :: ts                       ! surface temperature [K]
  real(c_double) :: ps                       ! surface pressure [Pa]
  real(c_double) :: coszrs                   ! cosine of zenith angle
  real(c_double) :: mwdry                    ! molecular weight of dry air [AMU]
  real(c_double) :: cpdry                    ! specific heat of dry air [J kg-1 K-1]
  real(c_double) :: srf_emiss                ! broadband surface thermal emissivity
  real(c_double) :: asdir, asdif             ! shortwave direct/diffuse surface albedo
  real(c_double) :: aldir, aldif             ! longwave(near-IR) direct/diffuse surface albedo
  ! P, T, Z
  real(c_double) :: tmid(pver)               ! mid-layer temperature [K]
  real(c_double) :: tint(pverp)              ! interface temperature [K]
  real(c_double) :: pmid(pver)               ! mid-layer pressure [Pa]
  real(c_double) :: pdel(pver)               ! layer pressure thickness [Pa]
  real(c_double) :: pint(pverp)              ! interface pressure [Pa]
  real(c_double) :: zint(pverp)              ! interface geopotential height [m]
  ! gases (h2ommr is specific humidity; the rest are dry mass mixing ratios)
  real(c_double) :: h2ommr(pver)
  real(c_double) :: co2mmr(pver)
  real(c_double) :: ch4mmr(pver)
  real(c_double) :: c2h6mmr(pver)
  real(c_double) :: nh3mmr(pver)
  real(c_double) :: commr(pver)
  real(c_double) :: o2mmr(pver)
  real(c_double) :: o3mmr(pver)
  real(c_double) :: n2mmr(pver)
  real(c_double) :: h2mmr(pver)
  ! H2O clouds (used only when do_exo_clouds was set at exort_init)
  real(c_double) :: cicewp(pver)             ! in-cloud ice water path [g m-2]
  real(c_double) :: cliqwp(pver)             ! in-cloud liquid water path [g m-2]
  real(c_double) :: cfrc(pver)               ! cloud fraction
  real(c_double) :: rei(pver)                ! ice effective radius [micron]
  real(c_double) :: rel(pver)                ! liquid effective radius [micron]
  ! CO2 ice clouds
  real(c_double) :: cicewp_co2(pver)
  real(c_double) :: rei_co2(pver)
  ! CARMA haze binwise mass mixing ratio (used only when do_exo_haze was set)
  real(c_double) :: carmammr(pver, nelem_carma, nbin_carma)
end type column_state_t

! zero-init helper for Fortran callers (C/Python callers get zeroed
! structs from calloc/ffi.new); bind(c) types cannot carry default
! component initialization, so this is the one canonical zeroing spot
public :: zero_column_state

type, bind(c) :: column_result_t
  ! scalars
  real(c_double) :: sol_toa                  ! TOA incident stellar flux [W m-2]
  real(c_double) :: vis_dir, vis_dif         ! surface downwelling visible direct/diffuse [W m-2]
  real(c_double) :: nir_dir, nir_dif         ! surface downwelling near-IR direct/diffuse [W m-2]
  ! heating rates [K s-1] (multiply by 86400 for the K day-1 written to RTprofile_out.nc)
  real(c_double) :: sw_dtdt(pver)
  real(c_double) :: lw_dtdt(pver)
  ! broadband fluxes [W m-2]
  real(c_double) :: lw_dnflux(pverp)
  real(c_double) :: lw_upflux(pverp)
  real(c_double) :: sw_dnflux(pverp)
  real(c_double) :: sw_upflux(pverp)
  ! per-band spectral fluxes [W m-2]
  real(c_double) :: lw_dnflux_spectral(pverp, ntot_wavlnrng)
  real(c_double) :: lw_upflux_spectral(pverp, ntot_wavlnrng)
  real(c_double) :: sw_dnflux_spectral(pverp, ntot_wavlnrng)
  real(c_double) :: sw_upflux_spectral(pverp, ntot_wavlnrng)
end type column_result_t

contains

elemental subroutine zero_column_state(s)
  type(column_state_t), intent(out) :: s

  s%ts = 0.0_c_double
  s%ps = 0.0_c_double
  s%coszrs = 0.0_c_double
  s%mwdry = 0.0_c_double
  s%cpdry = 0.0_c_double
  s%srf_emiss = 0.0_c_double
  s%asdir = 0.0_c_double
  s%asdif = 0.0_c_double
  s%aldir = 0.0_c_double
  s%aldif = 0.0_c_double
  s%tmid = 0.0_c_double
  s%tint = 0.0_c_double
  s%pmid = 0.0_c_double
  s%pdel = 0.0_c_double
  s%pint = 0.0_c_double
  s%zint = 0.0_c_double
  s%h2ommr = 0.0_c_double
  s%co2mmr = 0.0_c_double
  s%ch4mmr = 0.0_c_double
  s%c2h6mmr = 0.0_c_double
  s%nh3mmr = 0.0_c_double
  s%commr = 0.0_c_double
  s%o2mmr = 0.0_c_double
  s%o3mmr = 0.0_c_double
  s%n2mmr = 0.0_c_double
  s%h2mmr = 0.0_c_double
  s%cicewp = 0.0_c_double
  s%cliqwp = 0.0_c_double
  s%cfrc = 0.0_c_double
  s%rei = 0.0_c_double
  s%rel = 0.0_c_double
  s%cicewp_co2 = 0.0_c_double
  s%rei_co2 = 0.0_c_double
  s%carmammr = 0.0_c_double

end subroutine zero_column_state

end module exort_column_mod
