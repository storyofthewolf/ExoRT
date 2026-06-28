program main
!----------------------------------------------------
! Driver for 1d offline radiative transfer calculations using
! CAM radiative transfer code written by E.Wolf for Archean
! and exoplanetary atmospheres
!----------------------------------------------------

use shr_kind_mod,       only: r8 => shr_kind_r8
use radgrid
use exo_radiation_mod
use shr_const_mod
use io
use kabs
use exoplanet_mod
use ppgrid
use physconst
use initialize_rad_mod_1D
use exo_init_ref
use exo_model_specific

implicit none

integer :: k
real    :: t0, t1, t_input, t_init, t_kernel, t_output

!---- output variables ----
real(r8), dimension(pver) :: sw_dTdt_out
real(r8), dimension(pver) :: lw_dTdt_out
real(r8), dimension(pverp) :: lw_dnflux_out
real(r8), dimension(pverp) :: lw_upflux_out
real(r8), dimension(pverp) :: sw_upflux_out
real(r8), dimension(pverp) :: sw_dnflux_out
real(r8), dimension(pverp,ntot_wavlnrng) :: lw_dnflux_spectral_out
real(r8), dimension(pverp,ntot_wavlnrng) :: lw_upflux_spectral_out
real(r8), dimension(pverp,ntot_wavlnrng) :: sw_upflux_spectral_out
real(r8), dimension(pverp,ntot_wavlnrng) :: sw_dnflux_spectral_out
real(r8) :: vis_dir_out
real(r8) :: vis_dif_out
real(r8) :: nir_dir_out
real(r8) :: nir_dif_out
real(r8) :: sol_toa_out

! --- Runtime namelist: read user_nl_exort if present, else compile-time defaults ---
call read_namelist

call cpu_time(t0)
call initialize_kcoeff
call initialize_solar
if (do_exo_clouds) call initialize_cldopts
call init_ref
call init_model_specific
call init_planck
call initialize_radbuffer
call initialize_to_zero
call cpu_time(t1)
t_init = t1 - t0

call cpu_time(t0)
call input_profile
call physconst_setgas(MWDRY_in, CPDRY_in)
call cpu_time(t1)
t_input = t1 - t0

! --- random inputs ---
ext_msdist_in = 1.0
ext_rtgt_in = 1.0
ext_solar_azm_ang_in = 0.0
ext_tazm_ang_in = 0.0
ext_tslope_ang_in = 0.0
ext_tslas_tog_in = 0
ext_tshadow_tog_in = 1
ext_cosz_horizon_in(:) = 0.0
ext_TCx_obstruct_in(:) = 0.0
ext_TCz_obstruct_in(:) = 0.0


! define dry  as wet*(1-H2OMMr)
PDELDRY_in(:) = PDEL_in(:)*(1-H2OMMR_in(:))

! define dry interface pressure as wet*(1-H2OMMR), mapping mid-layer H2OMMR to
! interfaces the CESM way: top interface takes the top layer, each lower
! interface takes the mid-layer below it.
PINTDRY_in(1)       = PINT_in(1)*(1.-H2OMMR_in(1))
PINTDRY_in(2:pverp) = PINT_in(2:pverp)*(1.-H2OMMR_in(:))

call cpu_time(t0)
call aerad_driver(H2OMMR_in, CO2MMR_in, &
                  CH4MMR_in, C2H6MMR_in, &
                  NH3MMR_in, COMMR_in, &
                  H2MMR_in,  N2MMR_in, O3MMR_in, O2MMR_in, &
                  CICEWP_in, CLIQWP_in, CFRC_in,  &
                  REI_in, REL_in,  &
                  TS_in, PS_in, PMID_in,  &
                  PDEL_in, PDELDRY_in, TMID_in, PINT_in, PINTDRY_in,  &
                  COSZRS_in, ext_msdist_in,  &
                  ASDIR_in, ALDIR_in,  &
                  ASDIF_in, ALDIF_in,  &
                  ext_rtgt_in, ext_solar_azm_ang_in, ext_tazm_ang_in, ext_tslope_ang_in,   &
                  ext_tslas_tog_in, ext_tshadow_tog_in, ext_nazm_tshadow, ext_cosz_horizon_in,  &
                  ext_TCx_obstruct_in, ext_TCz_obstruct_in, ZINT_in,  &
                  sw_dTdt_out, lw_dTdt_out, &
                  lw_dnflux_out, lw_upflux_out, sw_upflux_out, sw_dnflux_out,  &
                  lw_dnflux_spectral_out, lw_upflux_spectral_out, sw_upflux_spectral_out, sw_dnflux_spectral_out,  &
                  vis_dir_out, vis_dif_out, nir_dir_out, nir_dif_out, sol_toa_out  )
call cpu_time(t1)
t_kernel = t1 - t0

! Print Primary Diagnostic outputs
call print_diagnostics( sol_toa_out, vis_dir_out, vis_dif_out, nir_dir_out, nir_dif_out, &
                        sw_dnflux_out, sw_upflux_out, lw_dnflux_out, lw_upflux_out )


call cpu_time(t0)
call output_data( sw_dTdt_out*SHR_CONST_CSEC, lw_dTdt_out*SHR_CONST_CSEC, &
                  lw_dnflux_out, lw_upflux_out, &
                  sw_dnflux_out, sw_upflux_out, &
                  lw_dnflux_spectral_out, lw_upflux_spectral_out, &
                  sw_dnflux_spectral_out, sw_upflux_spectral_out, &
                  sol_toa_out, &
                  PMID_in, PINT_in, TMID_in, &
                  TINT_in, ZINT_in, &
		  H2OMMR_in, CO2MMR_in, &
                  CH4MMR_in, C2H6MMR_in, &
                  NH3MMR_in, COMMR_in, &
                  O2MMR_in,  O3MMR_in,  N2MMR_in, H2MMR_in, &
                  MWDRY_in, CPDRY_in )
call cpu_time(t1)
t_output = t1 - t0

write(*,*) '=== timing (cpu seconds) ==================='
write(*,'(a,f10.4)') '  initialization : ', t_init
write(*,'(a,f10.4)') '  input          : ', t_input
write(*,'(a,f10.4)') '  aerad_driver   : ', t_kernel
write(*,'(a,f10.4)') '  output         : ', t_output
write(*,'(a,f10.4)') '  total          : ', t_init+t_input+t_kernel+t_output
write(*,*) '============================================'

end program main
