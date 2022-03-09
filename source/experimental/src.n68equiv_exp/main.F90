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
use input
use output
use kabs
use exoplanet_mod
use ppgrid
use physconst
use initialize_rad_mod_1D
use exo_init_ref

implicit none

integer :: k
!---- output variables ----
! note the clear vs. cloudy sky convention here is the opposite of the 3D model.
! in the 3D model the default is cloudy-sky, with a clearsky option.
! Here, the clear-sky computation is the default, and computations with clouds are the option.
real(r8), dimension(pver) :: sw_dTdt_out,     sw_dTdt_cld_out
real(r8), dimension(pver) :: lw_dTdt_out,     lw_dTdt_cld_out
real(r8), dimension(pverp) :: lw_dnflux_out,  lw_dnflux_cld_out
real(r8), dimension(pverp) :: lw_upflux_out,  lw_upflux_cld_out
real(r8), dimension(pverp) :: sw_upflux_out,  sw_upflux_cld_out 
real(r8), dimension(pverp) :: sw_dnflux_out,  sw_dnflux_cld_out
real(r8), dimension(pverp,ntot_wavlnrng) :: lw_dnflux_spectral_out, lw_dnflux_spectral_cld_out
real(r8), dimension(pverp,ntot_wavlnrng) :: lw_upflux_spectral_out, lw_upflux_spectral_cld_out
real(r8), dimension(pverp,ntot_wavlnrng) :: sw_upflux_spectral_out, sw_upflux_spectral_cld_out 
real(r8), dimension(pverp,ntot_wavlnrng) :: sw_dnflux_spectral_out, sw_dnflux_spectral_cld_out
real(r8) :: vis_dir_out, vis_dir_cld_out
real(r8) :: vis_dif_out, vis_dif_cld_out
real(r8) :: nir_dir_out, nir_dir_cld_out
real(r8) :: nir_dif_out, nir_dif_cld_out

! junk variables
real(r8), dimension(pverp) :: h2oint 

call initialize_kcoeff
call initialize_solar
if (do_exo_clouds) call initialize_cldopts
call init_ref
call init_planck
call initialize_radbuffer
call initialize_to_zero  ! inouts initialized to zero
call input_profile
call physconst_setgas(MWDRY_in, CPDRY_in)

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

! initialize output varibales to zero
sw_dTdt_out(:)       = 0.0  ;  sw_dTdt_cld_out(:) = 0.0
lw_dTdt_out(:)       = 0.0  ;  lw_dTdt_cld_out(:) = 0.0
lw_dnflux_out(:)     = 0.0  ;  lw_dnflux_cld_out(:)= 0.0
lw_upflux_out(:)     = 0.0  ;  lw_upflux_cld_out(:)= 0.0
sw_upflux_out(:)     = 0.0  ;  sw_upflux_cld_out(:) = 0.0
sw_dnflux_out(:)     = 0.0  ;  sw_dnflux_cld_out(:)= 0.0
lw_dnflux_spectral_out(:,:) = 0.0  ;  lw_dnflux_spectral_cld_out(:,:) = 0.0
lw_upflux_spectral_out(:,:) = 0.0  ;  lw_upflux_spectral_cld_out(:,:) = 0.0
sw_upflux_spectral_out(:,:) = 0.0  ;  sw_upflux_spectral_cld_out(:,:) = 0.0
sw_dnflux_spectral_out(:,:) = 0.0  ;  sw_dnflux_spectral_cld_out(:,:) = 0.0


! define dry  as wet*(!-H2OMMr)
PDELDRY_in(:) = PDEL_in(:)*(1-H2OMMR_in(:))

h2oint(1) = H2OMMR_in(1)
h2oint(2:pverp) = H2OMMR_in(:)
PINTDRY_in(:) = PINT_in(:)*(1.-h2oint(:))


! clear sky calculation
call aerad_driver(H2OMMR_in, CO2MMR_in, CH4MMR_in, &
                  H2MMR_in,  N2MMR_in, &
                  CICEWP_zero, CLIQWP_zero, CFRC_zero,  &
                  CICEWP_CO2_zero, &
                  REI_zero, REL_zero, REI_CO2_zero,  &
                  TS_in, PS_in, PMID_in,  &
                  PDEL_in, PDELDRY_in, TMID_in, PINT_in, PINTDRY_in,  &
                  COSZRS_in, ext_msdist_in,  &
                  ASDIR_in, ALDIR_in,  &
                  ASDIF_in, ALDIF_in,  &
                  SRF_EMISS_in, &
                  ext_rtgt_in, ext_solar_azm_ang_in, ext_tazm_ang_in, ext_tslope_ang_in,   &
                  ext_tslas_tog_in, ext_tshadow_tog_in, ext_nazm_tshadow, ext_cosz_horizon_in,  &
                  ext_TCx_obstruct_in, ext_TCz_obstruct_in, ZINT_in,  &
                  sw_dTdt_out, lw_dTdt_out, &
                  lw_dnflux_out, lw_upflux_out, sw_upflux_out, sw_dnflux_out,  &
                  lw_dnflux_spectral_out, lw_upflux_spectral_out, sw_upflux_spectral_out, sw_dnflux_spectral_out,  &
                  vis_dir_out, vis_dif_out, nir_dir_out, nir_dif_out  )

! Print Primary Diagnostic outputs
write(*,*) "============= CLEAR SKY ==============="
write(*,*) "Surface downwelling fluxes"
write(*,*) "vis_dir", vis_dir_out 
write(*,*) "vis_dif", vis_dif_out 
write(*,*) "nir_dir", nir_dir_out 
write(*,*) "nir_dif", nir_dif_out 
write(*,*) "total direct", vis_dir_out+nir_dir_out
write(*,*) "total diffuse", vis_dif_out+nir_dif_out
write(*,*) "SW DN TOA/SURF", sw_dnflux_out(1), sw_dnflux_out(pverp)
write(*,*) "SW UP TOA/SURF", sw_upflux_out(1), sw_upflux_out(pverp)
write(*,*) "LW DN TOA/SURF", lw_dnflux_out(1), lw_dnflux_out(pverp)
write(*,*) "LW UP TOA/SURF", lw_upflux_out(1), lw_upflux_out(pverp)



if (do_exo_clouds) then

  ! cloudy sky calculation
  call aerad_driver(H2OMMR_in, CO2MMR_in, CH4MMR_in, &
                    H2MMR_in,  N2MMR_in, &
                    CICEWP_in, CLIQWP_in, CFRC_in,  &
                    CICEWP_CO2_in, &
                    REI_in, REL_in, REI_CO2_in,  &                   
                    TS_in, PS_in, PMID_in,  &
                    PDEL_in, PDELDRY_in, TMID_in, PINT_in, PINTDRY_in,  &
                    COSZRS_in, ext_msdist_in,  &
                    ASDIR_in, ALDIR_in,  &
                    ASDIF_in, ALDIF_in,  &
                    SRF_EMISS_in, &
                    ext_rtgt_in, ext_solar_azm_ang_in, ext_tazm_ang_in, ext_tslope_ang_in,   &
                    ext_tslas_tog_in, ext_tshadow_tog_in, ext_nazm_tshadow, ext_cosz_horizon_in,  &
                    ext_TCx_obstruct_in, ext_TCz_obstruct_in, ZINT_in,  &
                    sw_dTdt_cld_out, lw_dTdt_cld_out, &
                    lw_dnflux_cld_out, lw_upflux_cld_out, sw_upflux_cld_out, sw_dnflux_cld_out,  &
                    lw_dnflux_spectral_cld_out, lw_upflux_spectral_cld_out, sw_upflux_spectral_cld_out, sw_dnflux_spectral_cld_out,  &
                    vis_dir_cld_out, vis_dif_cld_out, nir_dir_cld_out, nir_dif_cld_out  )

   ! Print Primary Diagnostic outputs
    write(*,*) "   "
    write(*,*) "============= CLOUDY SKY ==============="
    write(*,*) "Surface downwelling fluxes"
    write(*,*) "vis_dir", vis_dir_cld_out 
    write(*,*) "vis_dif", vis_dif_cld_out 
    write(*,*) "nir_dir", nir_dir_cld_out 
    write(*,*) "nir_dif", nir_dif_cld_out 
    write(*,*) "total direct", vis_dir_cld_out+nir_dir_cld_out
    write(*,*) "total diffuse", vis_dif_cld_out+nir_dif_cld_out
    write(*,*) "SW DN TOA/SURF", sw_dnflux_cld_out(1), sw_dnflux_cld_out(pverp)
    write(*,*) "SW UP TOA/SURF", sw_upflux_cld_out(1), sw_upflux_cld_out(pverp)
    write(*,*) "LW DN TOA/SURF", lw_dnflux_cld_out(1), lw_dnflux_cld_out(pverp)
    write(*,*) "LW UP TOA/SURF", lw_upflux_cld_out(1), lw_upflux_cld_out(pverp)

endif


call output_data( sw_dTdt_out*SHR_CONST_CSEC, lw_dTdt_out*SHR_CONST_CSEC, &
                  lw_dnflux_out, lw_upflux_out, &
                  sw_dnflux_out, sw_upflux_out, &
                  lw_dnflux_spectral_out, lw_upflux_spectral_out, &
                  sw_dnflux_spectral_out, sw_upflux_spectral_out, &
                  sw_dTdt_cld_out*SHR_CONST_CSEC, lw_dTdt_cld_out*SHR_CONST_CSEC, &
                  lw_dnflux_cld_out, lw_upflux_cld_out, &
                  sw_dnflux_cld_out, sw_upflux_cld_out, &
                  lw_dnflux_spectral_cld_out, lw_upflux_spectral_cld_out, &
                  sw_dnflux_spectral_cld_out, sw_upflux_spectral_cld_out, &
                  PMID_in, PINT_in, TMID_in, &
                  TINT_in, ZINT_in, &
		  H2OMMR_in, CO2MMR_in, CH4MMR_in, &
                  O2MMR_in,  O3MMR_in,  N2MMR_in, H2MMR_in, &
                  CICEWP_in, CLIQWP_in, CFRC_in,  &
                  CICEWP_CO2_in, &
                  REI_in, REL_in, REI_CO2_in, &
                  MWDRY_in, CPDRY_in ) 

end program main
