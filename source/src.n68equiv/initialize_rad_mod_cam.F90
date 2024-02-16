
module initialize_rad_mod_cam

! version h2on68
!
! read in and initialize radiaive transfer grids
!

use kabs
use exoplanet_mod, only: solar_file => exo_solar_file
use cloud
use radgrid
use spmd_utils,       only: masterproc
use sys_rootdir

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
    integer :: keff_id
    character(len=256) :: locfn
    character(len=256) :: filename
    integer :: ierr

!------------------------------------------------------------------------
!
! Start Code
!

    if ( masterproc ) then
      write (6, '(2x, a)') '_______________________________________________________'
      write (6, '(2x, a)') '_________ initializing gas absorption coeffs __________'
      write (6, '(2x, a)') '_______________________________________________________'
    endif

    ! Load K coefficients
    filename = trim(exort_rootdir)//trim(dirk_h2o)//trim(k_h2o_file)
    call getfil(filename, locfn, 0)
    call cam_pio_openfile(ncid, locfn, PIO_NOWRITE)
    ierr =  pio_inq_varid(ncid, 'data',   keff_id)
    ierr =  pio_get_var(ncid, keff_id, k_h2o)
    call pio_closefile(ncid)

    filename = trim(exort_rootdir)//trim(dirk_co2)//trim(k_co2_file)
    call getfil(filename, locfn, 0)
    call cam_pio_openfile(ncid, locfn, PIO_NOWRITE)
    ierr =  pio_inq_varid(ncid, 'data',   keff_id)
    ierr =  pio_get_var(ncid, keff_id, k_co2)
    call pio_closefile(ncid)

    filename = trim(exort_rootdir)//trim(dirk_ch4)//trim(k_ch4_file)
    call getfil(filename, locfn, 0)
    call cam_pio_openfile(ncid, locfn, PIO_NOWRITE)
    ierr =  pio_inq_varid(ncid, 'data',   keff_id)
    ierr =  pio_get_var(ncid, keff_id, k_ch4)
    call pio_closefile(ncid)

    filename = trim(exort_rootdir)//trim(dirk_c2h6)//trim(k_c2h6_file)
    call getfil(filename, locfn, 0)
    call cam_pio_openfile(ncid, locfn, PIO_NOWRITE)
    ierr =  pio_inq_varid(ncid, 'data',   keff_id)
    ierr =  pio_get_var(ncid, keff_id, k_c2h6)
    call pio_closefile(ncid)

    filename = trim(exort_rootdir)//trim(dirk_o3)//trim(k_o3_file)
    call getfil(filename, locfn, 0)
    call cam_pio_openfile(ncid, locfn, PIO_NOWRITE)
    ierr =  pio_inq_varid(ncid, 'data',   keff_id)
    ierr =  pio_get_var(ncid, keff_id, k_o3)
    call pio_closefile(ncid)

    filename = trim(exort_rootdir)//trim(dirk_o2)//trim(k_o2_file)
    call getfil(filename, locfn, 0)
    call cam_pio_openfile(ncid, locfn, PIO_NOWRITE)
    ierr =  pio_inq_varid(ncid, 'data',   keff_id)
    ierr =  pio_get_var(ncid, keff_id, k_o2)
    call pio_closefile(ncid)

    !! Load mtckd h2o continuum
    filename = trim(exort_rootdir)//trim(dirct)//trim(kh2o_mtckd_file)
    call getfil(filename, locfn, 0)
    call cam_pio_openfile(ncid, locfn, PIO_NOWRITE)
    ierr =  pio_inq_varid(ncid, 'KSELF',   keff_id)
    ierr =  pio_get_var(ncid, keff_id, kh2oself_mtckd)
    call pio_closefile(ncid)

    filename = trim(exort_rootdir)//trim(dirct)//trim(kh2o_mtckd_file)
    call getfil(filename, locfn, 0)
    call cam_pio_openfile(ncid, locfn, PIO_NOWRITE)
    ierr =  pio_inq_varid(ncid, 'KFRGN',   keff_id)
    ierr =  pio_get_var(ncid, keff_id, kh2ofrgn_mtckd)
    call pio_closefile(ncid)
    !! mtckd

    ! Load K coefficients, for n2n2 continuum
    filename = trim(exort_rootdir)//trim(dirci)//trim(kn2n2cia_file )
    call getfil(filename, locfn, 0)
    call cam_pio_openfile(ncid, locfn, PIO_NOWRITE)
    ierr =  pio_inq_varid(ncid, 'sigma',   keff_id)
    ierr =  pio_get_var(ncid, keff_id, kn2n2)
    call pio_closefile(ncid)

    ! Load K coefficients, for n2h2 continuum
    filename = trim(exort_rootdir)//trim(dirci)//trim(kn2h2cia_file )
    call getfil(filename, locfn, 0)
    call cam_pio_openfile(ncid, locfn, PIO_NOWRITE)
    ierr =  pio_inq_varid(ncid, 'sigma',   keff_id)
    ierr =  pio_get_var(ncid, keff_id, kn2h2)
    call pio_closefile(ncid)

    ! Load K coefficients, for h2h2 continuum
    filename = trim(exort_rootdir)//trim(dirci)//trim(kh2h2cia_file )
    call getfil(filename, locfn, 0)
    call cam_pio_openfile(ncid, locfn, PIO_NOWRITE)
    ierr =  pio_inq_varid(ncid, 'sigma',   keff_id)
    ierr =  pio_get_var(ncid, keff_id, kh2h2)
    call pio_closefile(ncid)

    ! Load K coefficients, for co2co2 lw continuum
    filename = trim(exort_rootdir)//trim(dirci)//trim(kco2co2cia_lw_file )
    call getfil(filename, locfn, 0)
    call cam_pio_openfile(ncid, locfn, PIO_NOWRITE)
    ierr =  pio_inq_varid(ncid, 'sigma',   keff_id)
    ierr =  pio_get_var(ncid, keff_id, kco2co2_lw)
    call pio_closefile(ncid)

    ! Load K coefficients, for co2co2 sw continuum
    filename = trim(exort_rootdir)//trim(dirci)//trim(kco2co2cia_sw_file )
    call getfil(filename, locfn, 0)
    call cam_pio_openfile(ncid, locfn, PIO_NOWRITE)
    ierr =  pio_inq_varid(ncid, 'sigma',   keff_id)
    ierr =  pio_get_var(ncid, keff_id, kco2co2_sw)
    call pio_closefile(ncid)

    ! Load K coefficients, for co2h2 continuum
    filename = trim(exort_rootdir)//trim(dirci)//trim(kco2h2cia_file )
    call getfil(filename, locfn, 0)
    call cam_pio_openfile(ncid, locfn, PIO_NOWRITE)
    ierr =  pio_inq_varid(ncid, 'sigma',   keff_id)
    ierr =  pio_get_var(ncid, keff_id, kco2h2)
    call pio_closefile(ncid)

    ! Load K coefficients, for co2h4 continuum
    filename = trim(exort_rootdir)//trim(dirci)//trim(kco2ch4cia_file )
    call getfil(filename, locfn, 0)
    call cam_pio_openfile(ncid, locfn, PIO_NOWRITE)
    ierr =  pio_inq_varid(ncid, 'sigma',   keff_id)
    ierr =  pio_get_var(ncid, keff_id, kco2ch4)
    call pio_closefile(ncid)


! broadcast optical constants to all nodes
#if ( defined SPMD )
    call mpibcast(k_h2o,  ntot_wavlnrng*ngauss_8gpt*kc_npress*kc_ntemp, mpir8, 0, mpicom)
    call mpibcast(k_co2,  ntot_wavlnrng*ngauss_8gpt*kc_npress*kc_ntemp, mpir8, 0, mpicom)
    call mpibcast(k_ch4,  ntot_wavlnrng*ngauss_8gpt*kc_npress*kc_ntemp, mpir8, 0, mpicom)
    call mpibcast(k_c2h6, ntot_wavlnrng*ngauss_8gpt*kc_npress*kc_ntemp, mpir8, 0, mpicom)

    call mpibcast(kh2oself_mtckd, ngauss_8gpt*ntot_wavlnrng*kmtckd_ntemp, mpir8, 0, mpicom)
    call mpibcast(kh2ofrgn_mtckd, ngauss_8gpt*ntot_wavlnrng*kmtckd_ntemp, mpir8, 0, mpicom)

    call mpibcast(kn2n2, ntot_wavlnrng*kn2n2_ntemp, mpir8, 0, mpicom)
    call mpibcast(kn2h2, ntot_wavlnrng*kn2h2_ntemp, mpir8, 0, mpicom)
    call mpibcast(kh2h2, ntot_wavlnrng*kh2h2_ntemp, mpir8, 0, mpicom)

    call mpibcast(kco2co2_sw, ntot_wavlnrng*kco2co2_sw_ntemp, mpir8, 0, mpicom)
    call mpibcast(kco2co2_lw, ntot_wavlnrng*kco2co2_lw_ntemp, mpir8, 0, mpicom)
    call mpibcast(kco2h2, ntot_wavlnrng*kco2h2_ntemp, mpir8, 0, mpicom)
    call mpibcast(kco2ch4, ntot_wavlnrng*kco2ch4_ntemp, mpir8, 0, mpicom)

#endif


  end subroutine initialize_kcoeff


!============================================================================

  subroutine initialize_solar

!------------------------------------------------------------------------
!
! Purpose:  Initialize solar data from input file.
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
    character(len=256) :: locfn
    integer :: solarflux_id
    integer :: S0_id
    integer :: ierr

    if (masterproc) then
        write(6,*) "INITIALIZING SOLAR SPECTRAL FILE"
    endif

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
    character(len=256) :: filename
    character(len=256) :: locfn
    integer :: ierr

!------------------------------------------------------------------------
!
! Start Code
!
    if (masterproc) then
      write(6,*) "CLDOPTS: INITIALIZING CLOUD OPTICAL PROPERTIES"
    endif

    ! Load K water cloud optics file
    filename = trim(exort_rootdir)//trim(dircld)//trim(cldoptsL_file)
    call getfil(filename, locfn, 0)
    call cam_pio_openfile(ncid, locfn, PIO_NOWRITE)
!!    ierr =  pio_inq_dimid(ncid, 'rel_bins',   bin_id)
!!    ierr =  pio_inquire_dimension(ncid, bin_id, len=ncldopt_lbins)
    ierr =  pio_inq_varid(ncid, 'Qext_liq',   q_id)
    ierr =  pio_get_var(ncid, q_id, Qcldliq)
    ierr =  pio_inq_varid(ncid, 'W_liq',   w_id)
    ierr =  pio_get_var(ncid, w_id, Wcldliq)
    ierr =  pio_inq_varid(ncid, 'G_liq',   g_id)
    ierr =  pio_get_var(ncid, g_id, Gcldliq)
    call pio_closefile(ncid)


    ! Load ice cloud optics file
    filename = trim(exort_rootdir)//trim(dircld)//trim(cldoptsI_file)
    call getfil(filename, locfn, 0)
    call cam_pio_openfile(ncid, locfn, PIO_NOWRITE)
!!    ierr =  pio_inq_dimid(ncid, 'rei_bins',   bin_id)
!!    ierr =  pio_inquire_dimension(ncid, bin_id, len=ncldopt_ibins)
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
