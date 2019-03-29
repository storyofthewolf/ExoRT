
module initialize_rad_mod_cam
! version.highco2
!
! read in and initialize radiaive transfer grids
!

use kabs
use exoplanet_mod, only: solar_file => exo_solar_file
use cloud
use radgrid
use spmd_utils,   only: masterproc

implicit none
private 
save

!
! Pubic Interfaces
!
  public :: initialize_kcoeff
  public :: initialize_solar
  public :: initialize_cldopts
  public :: initialize_radbuffer


!============================================================================
contains
!============================================================================

!============================================================================
!
! Public subroutines
!
!============================================================================

!============================================================================

  subroutine initialize_kcoeff

!------------------------------------------------------------------------
!
! Purpose:  Initialize k coefficient data from input file.
!
!------------------------------------------------------------------------
!
#if ( defined SPMD)
  use mpishorthand
#endif

    use ioFileMod, only: getfil
    use cam_pio_utils, only: cam_pio_openfile
    use pio,  only: pio_inq_varid, pio_get_var, pio_closefile, pio_nowrite,  &
              file_desc_t, var_desc_t


    implicit none
    include 'netcdf.inc'

!------------------------------------------------------------------------
!
! Local Variables
!
    type(file_desc_t) :: ncid
    integer :: gid
    integer :: pid
    integer :: tid
    integer :: wid
    integer :: nid
    integer :: keff_id_lower
    integer :: keff_id_upper
    integer :: keff_id
    integer :: npress
    integer :: ntemp
    integer :: nweights
    integer :: nbands
    integer :: ierr
    character(len=256) :: locfn, filename

!------------------------------------------------------------------------
!
! Start Code
!
    if ( masterproc ) then

      write (6, '(2x, a)') '_______________________________________________________'
      write (6, '(2x, a)') '_________ initializing gas absorption coeffs __________'
      write (6, '(2x, a)') '_______________________________________________________'

    end if     

    ! Load K coefficients, interval 1  
    call getfil(k_file, locfn, 0)
    call cam_pio_openfile(ncid, locfn, PIO_NOWRITE)
    ierr =  pio_inq_varid(ncid, 'KEFF',   keff_id)
    ierr =  pio_get_var(ncid, keff_id, k_data)
    call pio_closefile(ncid)

    k01_data(:,:,:,:) = k_data(:,1,:,:,:)
    k02_data(:,:,:,:) = k_data(:,2,:,:,:)
    k03_data(:,:,:,:) = k_data(:,3,:,:,:)
    k04_data(:,:,:,:) = k_data(:,4,:,:,:)
    k05_data(:,:,:,:) = k_data(:,5,:,:,:)
    k06_data(:,:,:,:) = k_data(:,6,:,:,:)
    k07_data(:,:,:,:) = k_data(:,7,:,:,:)
    k08_data(:,:,:,:) = k_data(:,8,:,:,:)
    k09_data(:,:,:,:) = k_data(:,9,:,:,:)
    k10_data(:,:,:,:) = k_data(:,10,:,:,:)
    k11_data(:,:,:,:) = k_data(:,11,:,:,:)
    k12_data(:,:,:,:) = k_data(:,12,:,:,:)
    k13_data(:,:,:,:) = k_data(:,13,:,:,:)
    k14_data(:,:,:,:) = k_data(:,14,:,:,:)
    k15_data(:,:,:,:) = k_data(:,15,:,:,:)
    k16_data(:,:,:,:) = k_data(:,16,:,:,:)
    k17_data(:,:,:,:) = k_data(:,17,:,:,:)
    k18_data(:,:,:,:) = k_data(:,18,:,:,:)
    k19_data(:,:,:,:) = k_data(:,19,:,:,:)
    k20_data(:,:,:,:) = k_data(:,20,:,:,:)
    k21_data(:,:,:,:) = k_data(:,21,:,:,:)
    k22_data(:,:,:,:) = k_data(:,22,:,:,:)
    k23_data(:,:,:,:) = k_data(:,23,:,:,:)
    k24_data(:,:,:,:) = k_data(:,24,:,:,:)
    k25_data(:,:,:,:) = k_data(:,25,:,:,:)
    k26_data(:,:,:,:) = k_data(:,26,:,:,:)
    k27_data(:,:,:,:) = k_data(:,27,:,:,:)
    k28_data(:,:,:,:) = k_data(:,28,:,:,:)
    k29_data(:,:,:,:) = k_data(:,29,:,:,:)
    k30_data(:,:,:,:) = k_data(:,30,:,:,:)
    k31_data(:,:,:,:) = k_data(:,31,:,:,:)
    k32_data(:,:,:,:) = k_data(:,32,:,:,:)
    k33_data(:,:,:,:) = k_data(:,33,:,:,:)
    k34_data(:,:,:,:) = k_data(:,34,:,:,:)
    k35_data(:,:,:,:) = k_data(:,35,:,:,:)
    k36_data(:,:,:,:) = k_data(:,36,:,:,:)
    k37_data(:,:,:,:) = k_data(:,37,:,:,:)
    k38_data(:,:,:,:) = k_data(:,38,:,:,:)
    k39_data(:,:,:,:) = k_data(:,39,:,:,:)
    k40_data(:,:,:,:) = k_data(:,40,:,:,:)
    k41_data(:,:,:,:) = k_data(:,41,:,:,:)
    k42_data(:,:,:,:) = k_data(:,42,:,:,:)
    k43_data(:,:,:,:) = k_data(:,43,:,:,:)
    k44_data(:,:,:,:) = k_data(:,44,:,:,:)
    k45_data(:,:,:,:) = k_data(:,45,:,:,:)
    k46_data(:,:,:,:) = k_data(:,46,:,:,:)
    k47_data(:,:,:,:) = k_data(:,47,:,:,:)
    k48_data(:,:,:,:) = k_data(:,48,:,:,:)
    k49_data(:,:,:,:) = k_data(:,49,:,:,:)
    k50_data(:,:,:,:) = k_data(:,50,:,:,:)
    k51_data(:,:,:,:) = k_data(:,51,:,:,:)
    k52_data(:,:,:,:) = k_data(:,52,:,:,:)
    k53_data(:,:,:,:) = k_data(:,53,:,:,:)
    k54_data(:,:,:,:) = k_data(:,54,:,:,:)
    k55_data(:,:,:,:) = k_data(:,55,:,:,:)
    k56_data(:,:,:,:) = k_data(:,56,:,:,:)
    k57_data(:,:,:,:) = k_data(:,57,:,:,:)
    k58_data(:,:,:,:) = k_data(:,58,:,:,:)
    k59_data(:,:,:,:) = k_data(:,59,:,:,:)
    k60_data(:,:,:,:) = k_data(:,60,:,:,:)
    k61_data(:,:,:,:) = k_data(:,61,:,:,:)
    k62_data(:,:,:,:) = k_data(:,62,:,:,:)
    k63_data(:,:,:,:) = k_data(:,63,:,:,:)

    ! Load K coefficients, water vapor self continuum
    filename = trim(dirct)//trim(kh2oself_file)
    call getfil(filename, locfn, 0)
    call cam_pio_openfile(ncid, locfn, PIO_NOWRITE)
    ierr =  pio_inq_varid(ncid, 'KSELF',   keff_id)
    ierr =  pio_get_var(ncid, keff_id, Kh2oself)
    call pio_closefile(ncid)

    !! Load K coefficients, for co2 continuum
    filename = trim(dirct)//trim(kco2cont_file)
    call getfil(filename, locfn, 0)
    call cam_pio_openfile(ncid, locfn, PIO_NOWRITE)
    ierr =  pio_inq_varid(ncid, 'KSELF',   keff_id)
    ierr =  pio_get_var(ncid, keff_id, Kco2cont_8gpt)
    call pio_closefile(ncid)

! broadcast optical constants to all nodes
#if ( defined SPMD )
      call mpibcast(k01_lower, ngauss_pts(1)*kc_npress_lower*kc_ntemp_lower, mpir8, 0, mpicom)
      call mpibcast(k01_upper, ngauss_pts(1)*kc_npress_upper*kc_ntemp_upper, mpir8, 0, mpicom)
      call mpibcast(k02_lower, ngauss_pts(2)*kc_npress_lower*kc_ntemp_lower, mpir8, 0, mpicom)
      call mpibcast(k02_upper, ngauss_pts(2)*kc_npress_upper*kc_ntemp_upper, mpir8, 0, mpicom)
      call mpibcast(k03_lower, ngauss_pts(3)*kc_npress_lower*kc_ntemp_lower*kc_nweight, mpir8, 0, mpicom)
      call mpibcast(k03_upper, ngauss_pts(3)*kc_npress_upper*kc_ntemp_upper*kc_nweight, mpir8, 0, mpicom)
      call mpibcast(k04_lower, ngauss_pts(4)*kc_npress_lower*kc_ntemp_lower*kc_nweight, mpir8, 0, mpicom)
      call mpibcast(k04_upper, ngauss_pts(4)*kc_npress_upper*kc_ntemp_upper*kc_nweight, mpir8, 0, mpicom)
      call mpibcast(k05_lower, ngauss_pts(5)*kc_npress_lower*kc_ntemp_lower*kc_nweight, mpir8, 0, mpicom)
      call mpibcast(k05_upper, ngauss_pts(5)*kc_npress_upper*kc_ntemp_upper*kc_nweight, mpir8, 0, mpicom)
      call mpibcast(k06_lower, ngauss_pts(6)*kc_npress_lower*kc_ntemp_lower*kc_nweight, mpir8, 0, mpicom)
      call mpibcast(k06_upper, ngauss_pts(6)*kc_npress_upper*kc_ntemp_upper*kc_nweight, mpir8, 0, mpicom)
      call mpibcast(k07_lower, ngauss_pts(7)*kc_npress_lower*kc_ntemp_lower*kc_nweight, mpir8, 0, mpicom)
      call mpibcast(k07_upper, ngauss_pts(7)*kc_npress_upper*kc_ntemp_upper*kc_nweight, mpir8, 0, mpicom)
      call mpibcast(k08_lower, ngauss_pts(8)*kc_npress_lower*kc_ntemp_lower*kc_nweight, mpir8, 0, mpicom)
      call mpibcast(k08_upper, ngauss_pts(8)*kc_npress_upper*kc_ntemp_upper*kc_nweight, mpir8, 0, mpicom)
      call mpibcast(k09_lower, ngauss_pts(9)*kc_npress_lower*kc_ntemp_lower*kc_nweight, mpir8, 0, mpicom)
      call mpibcast(k09_upper, ngauss_pts(9)*kc_npress_upper*kc_ntemp_upper*kc_nweight, mpir8, 0, mpicom)
      call mpibcast(k10_lower, ngauss_pts(10)*kc_npress_lower*kc_ntemp_lower*kc_nweight, mpir8, 0, mpicom)
      call mpibcast(k10_upper, ngauss_pts(10)*kc_npress_upper*kc_ntemp_upper*kc_nweight, mpir8, 0, mpicom)
      call mpibcast(k11_lower, ngauss_pts(11)*kc_npress_lower*kc_ntemp_lower*kc_nweight, mpir8, 0, mpicom)
      call mpibcast(k11_upper, ngauss_pts(11)*kc_npress_upper*kc_ntemp_upper*kc_nweight, mpir8, 0, mpicom)
      call mpibcast(k12_lower, ngauss_pts(12)*kc_npress_lower*kc_ntemp_lower*kc_nweight, mpir8, 0, mpicom)
      call mpibcast(k12_upper, ngauss_pts(12)*kc_npress_upper*kc_ntemp_upper*kc_nweight, mpir8, 0, mpicom)
      call mpibcast(k13_lower, ngauss_pts(13)*kc_npress_lower*kc_ntemp_lower*kc_nweight, mpir8, 0, mpicom)
      call mpibcast(k13_upper, ngauss_pts(13)*kc_npress_upper*kc_ntemp_upper*kc_nweight, mpir8, 0, mpicom)
      call mpibcast(k14_lower, ngauss_pts(14)*kc_npress_lower*kc_ntemp_lower, mpir8, 0, mpicom)
      call mpibcast(k14_upper, ngauss_pts(14)*kc_npress_upper*kc_ntemp_upper, mpir8, 0, mpicom)
      call mpibcast(k15_lower, ngauss_pts(15)*kc_npress_lower*kc_ntemp_lower, mpir8, 0, mpicom)
      call mpibcast(k15_upper, ngauss_pts(15)*kc_npress_upper*kc_ntemp_upper, mpir8, 0, mpicom)
      call mpibcast(k16_lower, ngauss_pts(16)*kc_npress_lower*kc_ntemp_lower*kc_nweight, mpir8, 0, mpicom)
      call mpibcast(k16_upper, ngauss_pts(16)*kc_npress_upper*kc_ntemp_upper*kc_nweight, mpir8, 0, mpicom)
      call mpibcast(k17_lower, ngauss_pts(17)*kc_npress_lower*kc_ntemp_lower*kc_nweight, mpir8, 0, mpicom)
      call mpibcast(k17_upper, ngauss_pts(17)*kc_npress_upper*kc_ntemp_upper*kc_nweight, mpir8, 0, mpicom)
      call mpibcast(k18_lower, ngauss_pts(18)*kc_npress_lower*kc_ntemp_lower*kc_nweight, mpir8, 0, mpicom)
      call mpibcast(k18_upper, ngauss_pts(18)*kc_npress_upper*kc_ntemp_upper*kc_nweight, mpir8, 0, mpicom)
      call mpibcast(k19_lower, ngauss_pts(19)*kc_npress_lower*kc_ntemp_lower*kc_nweight, mpir8, 0, mpicom)
      call mpibcast(k19_upper, ngauss_pts(19)*kc_npress_upper*kc_ntemp_upper*kc_nweight, mpir8, 0, mpicom)
      call mpibcast(k20_lower, ngauss_pts(20)*kc_npress_lower*kc_ntemp_lower, mpir8, 0, mpicom)
      call mpibcast(k20_upper, ngauss_pts(20)*kc_npress_upper*kc_ntemp_upper, mpir8, 0, mpicom)
      call mpibcast(k21_lower, ngauss_pts(21)*kc_npress_lower*kc_ntemp_lower*kc_nweight, mpir8, 0, mpicom)
      call mpibcast(k21_upper, ngauss_pts(21)*kc_npress_upper*kc_ntemp_upper*kc_nweight, mpir8, 0, mpicom)
      call mpibcast(k22_lower, ngauss_pts(22)*kc_npress_lower*kc_ntemp_lower, mpir8, 0, mpicom)
      call mpibcast(k22_upper, ngauss_pts(22)*kc_npress_upper*kc_ntemp_upper, mpir8, 0, mpicom)
      call mpibcast(k23_lower, ngauss_pts(23)*kc_npress_lower*kc_ntemp_lower, mpir8, 0, mpicom)
      call mpibcast(k23_upper, ngauss_pts(23)*kc_npress_upper*kc_ntemp_upper, mpir8, 0, mpicom)
      call mpibcast(k24_lower, ngauss_pts(24)*kc_npress_lower*kc_ntemp_lower, mpir8, 0, mpicom)
      call mpibcast(k24_upper, ngauss_pts(24)*kc_npress_upper*kc_ntemp_upper, mpir8, 0, mpicom)
      call mpibcast(k25_lower, ngauss_pts(25)*kc_npress_lower*kc_ntemp_lower, mpir8, 0, mpicom)
      call mpibcast(k25_upper, ngauss_pts(25)*kc_npress_upper*kc_ntemp_upper, mpir8, 0, mpicom)
      call mpibcast(kh2oself, ngH2O*ks_npress*ks_ntemp*ks_nweight*ntot_wavlnrng, mpir8, 0, mpicom)
      call mpibcast(kco2cont, ngCO2*ntot_wavlnrng, mpir8, 0, mpicom)
#endif


  end subroutine initialize_kcoeff


!============================================================================

  subroutine initialize_solar

!------------------------------------------------------------------------
!
! Purpose:  Initialize solar data from input file.
!
!------------------------------------------------------------------------

#if ( defined SPMD)
  use mpishorthand
#endif

    use ioFileMod, only: getfil
    use cam_pio_utils, only: cam_pio_openfile
    use pio,  only: pio_inq_varid, pio_get_var, pio_closefile, pio_nowrite,  &
              file_desc_t, var_desc_t

    implicit none
    include 'netcdf.inc'

!------------------------------------------------------------------------
!
! Local Variables
!

    type(file_desc_t) :: ncid
    character(len=256) :: locfn
    integer :: solarflux_id
    integer :: S0_id
    integer :: ierr

    ! Load solar data
    call getfil(solar_file, locfn, 0)
    call cam_pio_openfile(ncid, locfn, PIO_NOWRITE)
    ierr =  pio_inq_varid(ncid, 'S0',   S0_id)
    ierr =  pio_get_var(ncid, S0_id, S0)
    ierr =  pio_inq_varid(ncid, 'solarflux',   solarflux_id)
    ierr =  pio_get_var(ncid, solarflux_id, solarflux)
    call pio_closefile(ncid)


#if ( defined SPMD )
    call mpibcast(S0, 1, mpir8, 0, mpicom)
    call mpibcast(solarflux, ntot_wavlnrng, mpir8, 0, mpicom)
#endif

  end subroutine initialize_solar


!============================================================================

  subroutine initialize_cldopts

!------------------------------------------------------------------------
!
! Purpose:  Initialize the cloud optical constants from input file.
!
!------------------------------------------------------------------------
!
#if ( defined SPMD)
  use mpishorthand
#endif

    use ioFileMod, only: getfil
    use cam_pio_utils, only: cam_pio_openfile
    use pio,  only: pio_inq_varid, pio_get_var, pio_closefile, pio_nowrite,  &
                    file_desc_t, var_desc_t, pio_inq_dimid, pio_inquire_dimension


    implicit none
    include 'netcdf.inc'

!------------------------------------------------------------------------
!
! Local Variables
!    
    type(file_desc_t) :: ncid
    integer :: bin_id
    integer :: wav_id
    integer :: ncldopt_lbins
    integer :: ncldopt_lwavs
    integer :: ncldopt_ibins
    integer :: ncldopt_iwavs
    integer :: q_id
    integer :: w_id
    integer :: g_id
    character(len=256) :: locfn
    integer :: ierr

!------------------------------------------------------------------------
!
! Start Code
!

    !write(6,*) "CLDOPTS: INITIALIZING WATER CLOUD OPTICAL PROPERTIES"

    ! Load K water cloud optics file
    call getfil(cldoptsL_file, locfn, 0)
    call cam_pio_openfile(ncid, locfn, PIO_NOWRITE)
    ierr =  pio_inq_dimid(ncid, 'rel_bins',   bin_id)
    ierr =  pio_inquire_dimension(ncid, bin_id, len=ncldopt_lbins)
    ierr =  pio_inq_dimid(ncid, 'nwavlrng',   wav_id)
    ierr =  pio_inquire_dimension(ncid, wav_id, len=ncldopt_lwavs)

    if ( masterproc ) then      
      write(6,*) "CLDOPTS: nrel = ",ncldopt_lbins
      write(6,*) "CLDOPTS: nwavlrng = ",ncldopt_lwavs
    end if

    !if (ncldopt_lwavs .ne. ntot_wavlnrng .or. ncldopt_lbins .ne. nrel) then
    !  write(6,*) "CLDOPTS: file size mismatch, liquid" 
    !  call endrun
    !end if

    ierr =  pio_inq_varid(ncid, 'Qext_liq',   q_id)
    ierr =  pio_get_var(ncid, q_id, Qcldliq)
    ierr =  pio_inq_varid(ncid, 'W_liq',   w_id)
    ierr =  pio_get_var(ncid, w_id, Wcldliq)
    ierr =  pio_inq_varid(ncid, 'G_liq',   g_id)
    ierr =  pio_get_var(ncid, g_id, Gcldliq)

    call pio_closefile(ncid)

    !write(6,*) "CLDOPTS: INITIALIZING ICE OPTICAL PROPERTIES"

    ! Load ice cloud optics file
    call getfil(cldoptsI_file, locfn, 0)
    call cam_pio_openfile(ncid, locfn, PIO_NOWRITE)
    ierr =  pio_inq_dimid(ncid, 'rei_bins',   bin_id)
    ierr =  pio_inquire_dimension(ncid, bin_id, len=ncldopt_ibins)
    ierr =  pio_inq_dimid(ncid, 'nwavlrng',   wav_id)
    ierr =  pio_inquire_dimension(ncid, wav_id, len=ncldopt_iwavs)

    if ( masterproc ) then      
      write(6,*) "CLDOPTS: nrei = ",ncldopt_ibins
      write(6,*) "CLDOPTS: nwavlrng = ",ncldopt_iwavs
    endif

    !if (ncldopt_iwavs .ne. ntot_wavlnrng .or. ncldopt_ibins .ne. nrei) then
    !  write(6,*) "CLDOPTS: file size mismatch, ice" 
    !  call endrun
    !end if

    ierr =  pio_inq_varid(ncid, 'Qext_ice',   q_id)
    ierr =  pio_get_var(ncid, q_id, Qcldice)
    ierr =  pio_inq_varid(ncid, 'W_ice',   w_id)
    ierr =  pio_get_var(ncid, w_id, Wcldice)
    ierr =  pio_inq_varid(ncid, 'G_ice',   g_id)
    ierr =  pio_get_var(ncid, g_id, Gcldice)

    call pio_closefile(ncid)



! broadcast water cloud optical constants to all nodes
#if ( defined SPMD )
      call mpibcast(Qcldliq, nrel*ntot_wavlnrng, mpir8, 0, mpicom)
      call mpibcast(Wcldliq, nrel*ntot_wavlnrng, mpir8, 0, mpicom)
      call mpibcast(Gcldliq, nrel*ntot_wavlnrng, mpir8, 0, mpicom)
#endif

! broadcast ice cloud optical constants to all nodes
#if ( defined SPMD )
      call mpibcast(Qcldice, nrei*ntot_wavlnrng, mpir8, 0, mpicom)
      call mpibcast(Wcldice, nrei*ntot_wavlnrng, mpir8, 0, mpicom)
      call mpibcast(Gcldice, nrei*ntot_wavlnrng, mpir8, 0, mpicom)
#endif

  end subroutine initialize_cldopts


!============================================================================

subroutine initialize_radbuffer

!
! Initialize radiation buffer data
!

!#include <comhyb.h>

!   integer :: k

 !If the top model level is above ~90 km (0.1 Pa), set the top level to compute
 !longwave cooling to about 80 km (1 Pa)
 !  if (hypm(1) .lt. 0.1) then
 !     do k = 1, pver
 !        if (hypm(k) .lt. 1) camtop = k 
 !        ! set top of cloud layer for cloud overlap assumption (1 hpa)
 !        !if (hypm(k) .lt. 1.e2) ntopcld  = k
 !     end do
 !  else
      camtop  = 1
 !     ntopcld = 2
 !  end if
 !  nlevsRT = pverp-camtop+1   
 !  if (masterproc) then
 !     write (6,*) 'INITIALIZE_RADBUFFER: camtop =',camtop
 !     write (6,*) 'INITIALIZE_RADBUFFER: pressure:',hypm(camtop)
 !     write (6,*) 'INITIALIZE_RADBUFFER: nlevsRT:',nlevsRT
 !  endif
  return
end subroutine initialize_radbuffer

!====================================================================================



end module initialize_rad_mod_cam
