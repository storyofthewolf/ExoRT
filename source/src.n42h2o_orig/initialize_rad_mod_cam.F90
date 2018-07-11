
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
     
    ! Load K coefficients, interval 1  
    filename = trim(dirk)//trim(k01_file)
    call getfil(filename, locfn, 0)
    call cam_pio_openfile(ncid, locfn, PIO_NOWRITE)
    ierr =  pio_inq_varid(ncid, 'data',   keff_id)
    ierr =  pio_get_var(ncid, keff_id, K01_data)
    call pio_closefile(ncid)

    ! Load K coefficients, interval 2  
    filename = trim(dirk)//trim(k02_file)
    call getfil(filename, locfn, 0)
    call cam_pio_openfile(ncid, locfn, PIO_NOWRITE)
    ierr =  pio_inq_varid(ncid, 'data',   keff_id)
    ierr =  pio_get_var(ncid, keff_id, K02_data)
    call pio_closefile(ncid)

    ! Load K coefficients, interval 3  
    filename = trim(dirk)//trim(k03_file)
    call getfil(filename, locfn, 0)
    call cam_pio_openfile(ncid, locfn, PIO_NOWRITE)
    ierr =  pio_inq_varid(ncid, 'data',   keff_id)
    ierr =  pio_get_var(ncid, keff_id, K03_data)
    call pio_closefile(ncid)

    ! Load K coefficients, interval 4  
    filename = trim(dirk)//trim(k04_file)
    call getfil(filename, locfn, 0)
    call cam_pio_openfile(ncid, locfn, PIO_NOWRITE)
    ierr =  pio_inq_varid(ncid, 'data',   keff_id)
    ierr =  pio_get_var(ncid, keff_id, K04_data)
    call pio_closefile(ncid)

    ! Load K coefficients, interval 5  
    filename = trim(dirk)//trim(k05_file)
    call getfil(filename, locfn, 0)
    call cam_pio_openfile(ncid, locfn, PIO_NOWRITE)
    ierr =  pio_inq_varid(ncid, 'data',   keff_id)
    ierr =  pio_get_var(ncid, keff_id, K05_data)
    call pio_closefile(ncid)

    ! Load K coefficients, interval 6  
    filename = trim(dirk)//trim(k06_file)
    call getfil(filename, locfn, 0)
    call cam_pio_openfile(ncid, locfn, PIO_NOWRITE)
    ierr =  pio_inq_varid(ncid, 'data',   keff_id)
    ierr =  pio_get_var(ncid, keff_id, K06_data)
    call pio_closefile(ncid)

    ! Load K coefficients, interval 7  
    filename = trim(dirk)//trim(k07_file)
    call getfil(filename, locfn, 0)
    call cam_pio_openfile(ncid, locfn, PIO_NOWRITE)
    ierr =  pio_inq_varid(ncid, 'data',   keff_id)
    ierr =  pio_get_var(ncid, keff_id, K07_data)
    call pio_closefile(ncid)

    ! Load K coefficients, interval 8  
    filename = trim(dirk)//trim(k08_file)
    call getfil(filename, locfn, 0)
    call cam_pio_openfile(ncid, locfn, PIO_NOWRITE)
    ierr =  pio_inq_varid(ncid, 'data',   keff_id)
    ierr =  pio_get_var(ncid, keff_id, K08_data)
    call pio_closefile(ncid)

    ! Load K coefficients, interval 9  
    filename = trim(dirk)//trim(k09_file)
    call getfil(filename, locfn, 0)
    call cam_pio_openfile(ncid, locfn, PIO_NOWRITE)
    ierr =  pio_inq_varid(ncid, 'data',   keff_id)
    ierr =  pio_get_var(ncid, keff_id, K09_data)
    call pio_closefile(ncid)

    ! Load K coefficients, interval 10  
    filename = trim(dirk)//trim(k10_file)
    call getfil(filename, locfn, 0)
    call cam_pio_openfile(ncid, locfn, PIO_NOWRITE)
    ierr =  pio_inq_varid(ncid, 'data',   keff_id)
    ierr =  pio_get_var(ncid, keff_id, K10_data)
    call pio_closefile(ncid)

    ! Load K coefficients, interval 11  
    filename = trim(dirk)//trim(k11_file)
    call getfil(filename, locfn, 0)
    call cam_pio_openfile(ncid, locfn, PIO_NOWRITE)
    ierr =  pio_inq_varid(ncid, 'data',   keff_id)
    ierr =  pio_get_var(ncid, keff_id, K11_data)
    call pio_closefile(ncid)

    ! Load K coefficients, interval 12  
    filename = trim(dirk)//trim(k12_file)
    call getfil(filename, locfn, 0)
    call cam_pio_openfile(ncid, locfn, PIO_NOWRITE)
    ierr =  pio_inq_varid(ncid, 'data',   keff_id)
    ierr =  pio_get_var(ncid, keff_id, K12_data)
    call pio_closefile(ncid)

    ! Load K coefficients, interval 13  
    filename = trim(dirk)//trim(k13_file)
    call getfil(filename, locfn, 0)
    call cam_pio_openfile(ncid, locfn, PIO_NOWRITE)
    ierr =  pio_inq_varid(ncid, 'data',   keff_id)
    ierr =  pio_get_var(ncid, keff_id, K13_data)
    call pio_closefile(ncid)

    ! Load K coefficients, interval 14  
    filename = trim(dirk)//trim(k14_file)
    call getfil(filename, locfn, 0)
    call cam_pio_openfile(ncid, locfn, PIO_NOWRITE)
    ierr =  pio_inq_varid(ncid, 'data',   keff_id)
    ierr =  pio_get_var(ncid, keff_id, K14_data)
    call pio_closefile(ncid)

    ! Load K coefficients, interval 15  
    filename = trim(dirk)//trim(k15_file)
    call getfil(filename, locfn, 0)
    call cam_pio_openfile(ncid, locfn, PIO_NOWRITE)
    ierr =  pio_inq_varid(ncid, 'data',   keff_id)
    ierr =  pio_get_var(ncid, keff_id, K15_data)
    call pio_closefile(ncid)

    ! Load K coefficients, interval 16  
    filename = trim(dirk)//trim(k16_file)
    call getfil(filename, locfn, 0)
    call cam_pio_openfile(ncid, locfn, PIO_NOWRITE)
    ierr =  pio_inq_varid(ncid, 'data',   keff_id)
    ierr =  pio_get_var(ncid, keff_id, K16_data)
    call pio_closefile(ncid)

    ! Load K coefficients, interval 17  
    filename = trim(dirk)//trim(k17_file)
    call getfil(filename, locfn, 0)
    call cam_pio_openfile(ncid, locfn, PIO_NOWRITE)
    ierr =  pio_inq_varid(ncid, 'data',   keff_id)
    ierr =  pio_get_var(ncid, keff_id, K17_data)
    call pio_closefile(ncid)

    ! Load K coefficients, interval 18  
    filename = trim(dirk)//trim(k18_file)
    call getfil(filename, locfn, 0)
    call cam_pio_openfile(ncid, locfn, PIO_NOWRITE)
    ierr =  pio_inq_varid(ncid, 'data',   keff_id)
    ierr =  pio_get_var(ncid, keff_id, K18_data)
    call pio_closefile(ncid)

    ! Load K coefficients, interval 19  
    filename = trim(dirk)//trim(k19_file)
    call getfil(filename, locfn, 0)
    call cam_pio_openfile(ncid, locfn, PIO_NOWRITE)
    ierr =  pio_inq_varid(ncid, 'data',   keff_id)
    ierr =  pio_get_var(ncid, keff_id, K19_data)
    call pio_closefile(ncid)

    ! Load K coefficients, interval 20  
    filename = trim(dirk)//trim(k20_file)
    call getfil(filename, locfn, 0)
    call cam_pio_openfile(ncid, locfn, PIO_NOWRITE)
    ierr =  pio_inq_varid(ncid, 'data',   keff_id)
    ierr =  pio_get_var(ncid, keff_id, K20_data)
    call pio_closefile(ncid)

    ! Load K coefficients, interval 21  
    filename = trim(dirk)//trim(k21_file)
    call getfil(filename, locfn, 0)
    call cam_pio_openfile(ncid, locfn, PIO_NOWRITE)
    ierr =  pio_inq_varid(ncid, 'data',   keff_id)
    ierr =  pio_get_var(ncid, keff_id, K21_data)
    call pio_closefile(ncid)

    ! Load K coefficients, interval 22  
    filename = trim(dirk)//trim(k22_file)
    call getfil(filename, locfn, 0)
    call cam_pio_openfile(ncid, locfn, PIO_NOWRITE)
    ierr =  pio_inq_varid(ncid, 'data',   keff_id)
    ierr =  pio_get_var(ncid, keff_id, K22_data)
    call pio_closefile(ncid)

    ! Load K coefficients, interval 23  
    filename = trim(dirk)//trim(k23_file)
    call getfil(filename, locfn, 0)
    call cam_pio_openfile(ncid, locfn, PIO_NOWRITE)
    ierr =  pio_inq_varid(ncid, 'data',   keff_id)
    ierr =  pio_get_var(ncid, keff_id, K23_data)
    call pio_closefile(ncid)

    ! Load K coefficients, interval 24  
    filename = trim(dirk)//trim(k24_file)
    call getfil(filename, locfn, 0)
    call cam_pio_openfile(ncid, locfn, PIO_NOWRITE)
    ierr =  pio_inq_varid(ncid, 'data',   keff_id)
    ierr =  pio_get_var(ncid, keff_id, K24_data)
    call pio_closefile(ncid)

    ! Load K coefficients, interval 25  
    filename = trim(dirk)//trim(k25_file)
    call getfil(filename, locfn, 0)
    call cam_pio_openfile(ncid, locfn, PIO_NOWRITE)
    ierr =  pio_inq_varid(ncid, 'data',   keff_id)
    ierr =  pio_get_var(ncid, keff_id, K25_data)
    call pio_closefile(ncid)

    ! Load K coefficients, interval 26  
    filename = trim(dirk)//trim(k26_file)
    call getfil(filename, locfn, 0)
    call cam_pio_openfile(ncid, locfn, PIO_NOWRITE)
    ierr =  pio_inq_varid(ncid, 'data',   keff_id)
    ierr =  pio_get_var(ncid, keff_id, K26_data)
    call pio_closefile(ncid)

    ! Load K coefficients, interval 27  
    filename = trim(dirk)//trim(k27_file)
    call getfil(filename, locfn, 0)
    call cam_pio_openfile(ncid, locfn, PIO_NOWRITE)
    ierr =  pio_inq_varid(ncid, 'data',   keff_id)
    ierr =  pio_get_var(ncid, keff_id, K27_data)
    call pio_closefile(ncid)

    ! Load K coefficients, interval 28  
    filename = trim(dirk)//trim(k28_file)
    call getfil(filename, locfn, 0)
    call cam_pio_openfile(ncid, locfn, PIO_NOWRITE)
    ierr =  pio_inq_varid(ncid, 'data',   keff_id)
    ierr =  pio_get_var(ncid, keff_id, K28_data)
    call pio_closefile(ncid)

    ! Load K coefficients, interval 29  
    filename = trim(dirk)//trim(k29_file)
    call getfil(filename, locfn, 0)
    call cam_pio_openfile(ncid, locfn, PIO_NOWRITE)
    ierr =  pio_inq_varid(ncid, 'data',   keff_id)
    ierr =  pio_get_var(ncid, keff_id, K29_data)
    call pio_closefile(ncid)

    ! Load K coefficients, interval 30  
    filename = trim(dirk)//trim(k30_file)
    call getfil(filename, locfn, 0)
    call cam_pio_openfile(ncid, locfn, PIO_NOWRITE)
    ierr =  pio_inq_varid(ncid, 'data',   keff_id)
    ierr =  pio_get_var(ncid, keff_id, K30_data)
    call pio_closefile(ncid)

    ! Load K coefficients, interval 31  
    filename = trim(dirk)//trim(k31_file)
    call getfil(filename, locfn, 0)
    call cam_pio_openfile(ncid, locfn, PIO_NOWRITE)
    ierr =  pio_inq_varid(ncid, 'data',   keff_id)
    ierr =  pio_get_var(ncid, keff_id, K31_data)
    call pio_closefile(ncid)

    ! Load K coefficients, interval 32  
    filename = trim(dirk)//trim(k32_file)
    call getfil(filename, locfn, 0)
    call cam_pio_openfile(ncid, locfn, PIO_NOWRITE)
    ierr =  pio_inq_varid(ncid, 'data',   keff_id)
    ierr =  pio_get_var(ncid, keff_id, K32_data)
    call pio_closefile(ncid)

    ! Load K coefficients, interval 33  
    filename = trim(dirk)//trim(k33_file)
    call getfil(filename, locfn, 0)
    call cam_pio_openfile(ncid, locfn, PIO_NOWRITE)
    ierr =  pio_inq_varid(ncid, 'data',   keff_id)
    ierr =  pio_get_var(ncid, keff_id, K33_data)
    call pio_closefile(ncid)

    ! Load K coefficients, interval 34  
    filename = trim(dirk)//trim(k34_file)
    call getfil(filename, locfn, 0)
    call cam_pio_openfile(ncid, locfn, PIO_NOWRITE)
    ierr =  pio_inq_varid(ncid, 'data',   keff_id)
    ierr =  pio_get_var(ncid, keff_id, K34_data)
    call pio_closefile(ncid)

    ! Load K coefficients, interval 35  
    filename = trim(dirk)//trim(k35_file)
    call getfil(filename, locfn, 0)
    call cam_pio_openfile(ncid, locfn, PIO_NOWRITE)
    ierr =  pio_inq_varid(ncid, 'data',   keff_id)
    ierr =  pio_get_var(ncid, keff_id, K35_data)
    call pio_closefile(ncid)

    ! Load K coefficients, interval 36  
    filename = trim(dirk)//trim(k36_file)
    call getfil(filename, locfn, 0)
    call cam_pio_openfile(ncid, locfn, PIO_NOWRITE)
    ierr =  pio_inq_varid(ncid, 'data',   keff_id)
    ierr =  pio_get_var(ncid, keff_id, K36_data)
    call pio_closefile(ncid)

    ! Load K coefficients, interval 37  
    filename = trim(dirk)//trim(k37_file)
    call getfil(filename, locfn, 0)
    call cam_pio_openfile(ncid, locfn, PIO_NOWRITE)
    ierr =  pio_inq_varid(ncid, 'data',   keff_id)
    ierr =  pio_get_var(ncid, keff_id, K37_data)
    call pio_closefile(ncid)

    ! Load K coefficients, interval 38  
    filename = trim(dirk)//trim(k38_file)
    call getfil(filename, locfn, 0)
    call cam_pio_openfile(ncid, locfn, PIO_NOWRITE)
    ierr =  pio_inq_varid(ncid, 'data',   keff_id)
    ierr =  pio_get_var(ncid, keff_id, K38_data)
    call pio_closefile(ncid)

    ! Load K coefficients, interval 39  
    filename = trim(dirk)//trim(k39_file)
    call getfil(filename, locfn, 0)
    call cam_pio_openfile(ncid, locfn, PIO_NOWRITE)
    ierr =  pio_inq_varid(ncid, 'data',   keff_id)
    ierr =  pio_get_var(ncid, keff_id, K39_data)
    call pio_closefile(ncid)

    ! Load K coefficients, interval 40  
    filename = trim(dirk)//trim(k40_file)
    call getfil(filename, locfn, 0)
    call cam_pio_openfile(ncid, locfn, PIO_NOWRITE)
    ierr =  pio_inq_varid(ncid, 'data',   keff_id)
    ierr =  pio_get_var(ncid, keff_id, K40_data)
    call pio_closefile(ncid)

    ! Load K coefficients, interval 41  
    filename = trim(dirk)//trim(k41_file)
    call getfil(filename, locfn, 0)
    call cam_pio_openfile(ncid, locfn, PIO_NOWRITE)
    ierr =  pio_inq_varid(ncid, 'data',   keff_id)
    ierr =  pio_get_var(ncid, keff_id, K41_data)
    call pio_closefile(ncid)

    ! Load K coefficients, interval 42  
    filename = trim(dirk)//trim(k42_file)
    call getfil(filename, locfn, 0)
    call cam_pio_openfile(ncid, locfn, PIO_NOWRITE)
    ierr =  pio_inq_varid(ncid, 'data',   keff_id)
    ierr =  pio_get_var(ncid, keff_id, K42_data)
    call pio_closefile(ncid)

    !! Load mtckd h2o self continuum
    filename = trim(dirct)//trim(kh2oself_mtckd_file)
    call getfil(filename, locfn, 0)
    call cam_pio_openfile(ncid, locfn, PIO_NOWRITE)
    ierr =  pio_inq_varid(ncid, 'KSELF',   keff_id)
    ierr =  pio_get_var(ncid, keff_id, kh2oself_mtckd)
    call pio_closefile(ncid)
    !! mtckd                         

    !! Load bps h2o vapor self continuum
    filename = trim(dirct)//trim(kh2oself_bps_file)
    call getfil(filename, locfn, 0)
    call cam_pio_openfile(ncid, locfn, PIO_NOWRITE)
    ierr =  pio_inq_varid(ncid, 'self',   keff_id)
    ierr =  pio_get_var(ncid, keff_id, self)
    ierr =  pio_inq_varid(ncid, 'foreign',   keff_id)
    ierr =  pio_get_var(ncid, keff_id, foreign)
    ierr =  pio_inq_varid(ncid, 'base_self',   keff_id)
    ierr =  pio_get_var(ncid, keff_id, base_self)
    ierr =  pio_inq_varid(ncid, 'base_foreign',   keff_id)
    ierr =  pio_get_var(ncid, keff_id, base_foreign)
    ierr =  pio_inq_varid(ncid, 'TempCoeff',   keff_id)
    ierr =  pio_get_var(ncid, keff_id, TempCoeff)
    call pio_closefile(ncid)

    ! Load K coefficients, for n2n2 continuum                                                                               
    filename = trim(dirci)//trim(kn2n2cia_file )
    call getfil(filename, locfn, 0)
    call cam_pio_openfile(ncid, locfn, PIO_NOWRITE)
    ierr =  pio_inq_varid(ncid, 'sigma',   keff_id)
    ierr =  pio_get_var(ncid, keff_id, Kn2n2)
    call pio_closefile(ncid)

    ! Load K coefficients, for h2n2 continuum                                                                               
    !filename = trim(dirci)//trim(kh2n2cia_file )
    !call getfil(filename, locfn, 0)
    !call cam_pio_openfile(ncid, locfn, PIO_NOWRITE)
    !ierr =  pio_inq_varid(ncid, 'sigma',   keff_id)
    !ierr =  pio_get_var(ncid, keff_id, Kh2n2)
    !call pio_closefile(ncid)

    ! Load K coefficients, for h2h2 continuum                                                                               
    !filename = trim(dirci)//trim(kh2h2cia_file )
    !call getfil(filename, locfn, 0)
    !call cam_pio_openfile(ncid, locfn, PIO_NOWRITE)
    !ierr =  pio_inq_varid(ncid, 'sigma',   keff_id)
    !ierr =  pio_get_var(ncid, keff_id, Kh2h2)
    !call pio_closefile(ncid)

  
! broadcast optical constants to all nodes
#if ( defined SPMD )
    call mpibcast(k01_data, ngauss_pts(1)*kc_npress*kc_ntemp, mpir8, 0, mpicom)
    call mpibcast(k02_data, ngauss_pts(2)*kc_npress*kc_ntemp, mpir8, 0, mpicom)
    call mpibcast(k03_data, ngauss_pts(3)*kc_npress*kc_ntemp, mpir8, 0, mpicom)
    call mpibcast(k04_data, ngauss_pts(4)*kc_npress*kc_ntemp, mpir8, 0, mpicom)
    call mpibcast(k05_data, ngauss_pts(5)*kc_npress*kc_ntemp, mpir8, 0, mpicom)
    call mpibcast(k06_data, ngauss_pts(6)*kc_npress*kc_ntemp, mpir8, 0, mpicom)
    call mpibcast(k07_data, ngauss_pts(7)*kc_npress*kc_ntemp, mpir8, 0, mpicom)
    call mpibcast(k08_data, ngauss_pts(8)*kc_npress*kc_ntemp, mpir8, 0, mpicom)
    call mpibcast(k09_data, ngauss_pts(9)*kc_npress*kc_ntemp, mpir8, 0, mpicom)
    call mpibcast(k10_data, ngauss_pts(10)*kc_npress*kc_ntemp, mpir8, 0, mpicom)
    call mpibcast(k11_data, ngauss_pts(11)*kc_npress*kc_ntemp, mpir8, 0, mpicom)
    call mpibcast(k12_data, ngauss_pts(12)*kc_npress*kc_ntemp, mpir8, 0, mpicom)
    call mpibcast(k13_data, ngauss_pts(13)*kc_npress*kc_ntemp, mpir8, 0, mpicom)
    call mpibcast(k14_data, ngauss_pts(14)*kc_npress*kc_ntemp, mpir8, 0, mpicom)
    call mpibcast(k15_data, ngauss_pts(15)*kc_npress*kc_ntemp, mpir8, 0, mpicom)
    call mpibcast(k16_data, ngauss_pts(16)*kc_npress*kc_ntemp, mpir8, 0, mpicom)
    call mpibcast(k17_data, ngauss_pts(17)*kc_npress*kc_ntemp, mpir8, 0, mpicom)
    call mpibcast(k18_data, ngauss_pts(18)*kc_npress*kc_ntemp, mpir8, 0, mpicom)
    call mpibcast(k19_data, ngauss_pts(19)*kc_npress*kc_ntemp, mpir8, 0, mpicom)
    call mpibcast(k20_data, ngauss_pts(20)*kc_npress*kc_ntemp, mpir8, 0, mpicom)
    call mpibcast(k21_data, ngauss_pts(21)*kc_npress*kc_ntemp, mpir8, 0, mpicom)
    call mpibcast(k22_data, ngauss_pts(22)*kc_npress*kc_ntemp, mpir8, 0, mpicom)
    call mpibcast(k23_data, ngauss_pts(23)*kc_npress*kc_ntemp, mpir8, 0, mpicom)
    call mpibcast(k24_data, ngauss_pts(24)*kc_npress*kc_ntemp, mpir8, 0, mpicom)
    call mpibcast(k25_data, ngauss_pts(25)*kc_npress*kc_ntemp, mpir8, 0, mpicom)
    call mpibcast(k26_data, ngauss_pts(26)*kc_npress*kc_ntemp, mpir8, 0, mpicom)
    call mpibcast(k27_data, ngauss_pts(27)*kc_npress*kc_ntemp, mpir8, 0, mpicom)
    call mpibcast(k28_data, ngauss_pts(28)*kc_npress*kc_ntemp, mpir8, 0, mpicom)
    call mpibcast(k29_data, ngauss_pts(29)*kc_npress*kc_ntemp, mpir8, 0, mpicom)
    call mpibcast(k30_data, ngauss_pts(30)*kc_npress*kc_ntemp, mpir8, 0, mpicom)
    call mpibcast(k31_data, ngauss_pts(31)*kc_npress*kc_ntemp, mpir8, 0, mpicom)
    call mpibcast(k32_data, ngauss_pts(32)*kc_npress*kc_ntemp, mpir8, 0, mpicom)
    call mpibcast(k33_data, ngauss_pts(33)*kc_npress*kc_ntemp, mpir8, 0, mpicom)
    call mpibcast(k34_data, ngauss_pts(34)*kc_npress*kc_ntemp, mpir8, 0, mpicom)
    call mpibcast(k35_data, ngauss_pts(35)*kc_npress*kc_ntemp, mpir8, 0, mpicom)
    call mpibcast(k36_data, ngauss_pts(36)*kc_npress*kc_ntemp, mpir8, 0, mpicom)
    call mpibcast(k37_data, ngauss_pts(37)*kc_npress*kc_ntemp, mpir8, 0, mpicom)
    call mpibcast(k38_data, ngauss_pts(38)*kc_npress*kc_ntemp, mpir8, 0, mpicom)
    call mpibcast(k39_data, ngauss_pts(39)*kc_npress*kc_ntemp, mpir8, 0, mpicom)
    call mpibcast(k40_data, ngauss_pts(40)*kc_npress*kc_ntemp, mpir8, 0, mpicom)
    call mpibcast(k41_data, ngauss_pts(41)*kc_npress*kc_ntemp, mpir8, 0, mpicom)
    call mpibcast(k42_data, ngauss_pts(42)*kc_npress*kc_ntemp, mpir8, 0, mpicom)

    call mpibcast(self, ntot_wavlnrng, mpir8, 0, mpicom)
    call mpibcast(base_self, ntot_wavlnrng, mpir8, 0, mpicom)
    call mpibcast(foreign, ntot_wavlnrng, mpir8, 0, mpicom)
    call mpibcast(base_foreign, ntot_wavlnrng, mpir8, 0, mpicom)
    call mpibcast(TempCoeff, ntot_wavlnrng, mpir8, 0, mpicom)

    call mpibcast(kh2oself_mtckd, ntot_wavlnrng*ks_ntemp, mpir8, 0, mpicom)

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
    call getfil(cldoptsL_file, locfn, 0)
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
    call getfil(cldoptsI_file, locfn, 0)
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
