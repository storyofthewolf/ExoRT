
module initialize_rad_mod_1D

! version n42h2o
!
! read in and initialize radiaive transfer grids
!

use kabs
use exoplanet_mod, only: solar_file, dirsol
use cloud  
use radgrid
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
!
!#if ( defined SPMD)
!  use mpishorthand
!#endif

    use ioFileMod, only: getfil

    implicit none
    include 'netcdf.inc'

!------------------------------------------------------------------------
!
! Local Variables
!
    integer :: ncid
    integer :: gid
    integer :: pid
    integer :: tid
    integer :: wid
    integer :: nid
    integer :: keff_id
    integer :: npress
    integer :: ntemp
    integer :: nweights
    integer :: nbands
    character(len=256) :: locfn, filename

!------------------------------------------------------------------------
!
! Start Code
!
!    if ( masterproc ) then

      write (6, '(2x, a)') '_______________________________________________________'
      write (6, '(2x, a)') '_________ initializing gas absorption coeffs __________'
      write (6, '(2x, a)') '_______________________________________________________'
     
      !----  Load K coefficients, interval 1  ----
      filename = trim(exort_rootdir)//trim(dirk_h2o)//trim(k01_h2o_file)
      call getfil(filename, locfn, 0)
      call wrap_open(locfn, 0, ncid)
      call wrap_inq_varid(ncid, 'data', keff_id)
      call wrap_get_var_realx(ncid, keff_id, k01_h2o)

      filename = trim(exort_rootdir)//trim(dirk_co2)//trim(k01_co2_file)
      call getfil(filename, locfn, 0)
      call wrap_open(locfn, 0, ncid)
      call wrap_inq_varid(ncid, 'data', keff_id)
      call wrap_get_var_realx(ncid, keff_id, k01_co2)

      filename = trim(exort_rootdir)//trim(dirk_ch4)//trim(k01_ch4_file)
      call getfil(filename, locfn, 0)
      call wrap_open(locfn, 0, ncid)
      call wrap_inq_varid(ncid, 'data', keff_id)
      call wrap_get_var_realx(ncid, keff_id, k01_ch4)


      !----  Load K coefficients, interval 2  ----
      filename = trim(exort_rootdir)//trim(dirk_h2o)//trim(k02_h2o_file)
      call getfil(filename, locfn, 0)
      call wrap_open(locfn, 0, ncid)
      call wrap_inq_varid(ncid, 'data', keff_id)
      call wrap_get_var_realx(ncid, keff_id, k02_h2o)

      filename = trim(exort_rootdir)//trim(dirk_co2)//trim(k02_co2_file)
      call getfil(filename, locfn, 0)
      call wrap_open(locfn, 0, ncid)
      call wrap_inq_varid(ncid, 'data', keff_id)
      call wrap_get_var_realx(ncid, keff_id, k02_co2)

      filename = trim(exort_rootdir)//trim(dirk_ch4)//trim(k02_ch4_file)
      call getfil(filename, locfn, 0)
      call wrap_open(locfn, 0, ncid)
      call wrap_inq_varid(ncid, 'data', keff_id)
      call wrap_get_var_realx(ncid, keff_id, k02_ch4)


      !----  Load K coefficients, interval 3  ----
      filename = trim(exort_rootdir)//trim(dirk_h2o)//trim(k03_h2o_file)
      call getfil(filename, locfn, 0)
      call wrap_open(locfn, 0, ncid)
      call wrap_inq_varid(ncid, 'data', keff_id)
      call wrap_get_var_realx(ncid, keff_id, k03_h2o)

      filename = trim(exort_rootdir)//trim(dirk_co2)//trim(k03_co2_file)
      call getfil(filename, locfn, 0)
      call wrap_open(locfn, 0, ncid)
      call wrap_inq_varid(ncid, 'data', keff_id)
      call wrap_get_var_realx(ncid, keff_id, k03_co2)

      filename = trim(exort_rootdir)//trim(dirk_ch4)//trim(k03_ch4_file)
      call getfil(filename, locfn, 0)
      call wrap_open(locfn, 0, ncid)
      call wrap_inq_varid(ncid, 'data', keff_id)
      call wrap_get_var_realx(ncid, keff_id, k03_ch4)


      !----  Load K coefficients, interval 4  ----
      filename = trim(exort_rootdir)//trim(dirk_h2o)//trim(k04_h2o_file)
      call getfil(filename, locfn, 0)
      call wrap_open(locfn, 0, ncid)
      call wrap_inq_varid(ncid, 'data', keff_id)
      call wrap_get_var_realx(ncid, keff_id, k04_h2o)

      filename = trim(exort_rootdir)//trim(dirk_co2)//trim(k04_co2_file)
      call getfil(filename, locfn, 0)
      call wrap_open(locfn, 0, ncid)
      call wrap_inq_varid(ncid, 'data', keff_id)
      call wrap_get_var_realx(ncid, keff_id, k04_co2)

      filename = trim(exort_rootdir)//trim(dirk_ch4)//trim(k04_ch4_file)
      call getfil(filename, locfn, 0)
      call wrap_open(locfn, 0, ncid)
      call wrap_inq_varid(ncid, 'data', keff_id)
      call wrap_get_var_realx(ncid, keff_id, k04_ch4)


      !----  Load K coefficients, interval 5  ----
      filename = trim(exort_rootdir)//trim(dirk_h2o)//trim(k05_h2o_file)
      call getfil(filename, locfn, 0)
      call wrap_open(locfn, 0, ncid)
      call wrap_inq_varid(ncid, 'data', keff_id)
      call wrap_get_var_realx(ncid, keff_id, k05_h2o)

      filename = trim(exort_rootdir)//trim(dirk_co2)//trim(k05_co2_file)
      call getfil(filename, locfn, 0)
      call wrap_open(locfn, 0, ncid)
      call wrap_inq_varid(ncid, 'data', keff_id)
      call wrap_get_var_realx(ncid, keff_id, k05_co2)

      filename = trim(exort_rootdir)//trim(dirk_ch4)//trim(k05_ch4_file)
      call getfil(filename, locfn, 0)
      call wrap_open(locfn, 0, ncid)
      call wrap_inq_varid(ncid, 'data', keff_id)
      call wrap_get_var_realx(ncid, keff_id, k05_ch4)


      !----  Load K coefficients, interval 6  ----
      filename = trim(exort_rootdir)//trim(dirk_h2o)//trim(k06_h2o_file)
      call getfil(filename, locfn, 0)
      call wrap_open(locfn, 0, ncid)
      call wrap_inq_varid(ncid, 'data', keff_id)
      call wrap_get_var_realx(ncid, keff_id, k06_h2o)

      filename = trim(exort_rootdir)//trim(dirk_co2)//trim(k06_co2_file)
      call getfil(filename, locfn, 0)
      call wrap_open(locfn, 0, ncid)
      call wrap_inq_varid(ncid, 'data', keff_id)
      call wrap_get_var_realx(ncid, keff_id, k06_co2)

      filename = trim(exort_rootdir)//trim(dirk_ch4)//trim(k06_ch4_file)
      call getfil(filename, locfn, 0)
      call wrap_open(locfn, 0, ncid)
      call wrap_inq_varid(ncid, 'data', keff_id)
      call wrap_get_var_realx(ncid, keff_id, k06_ch4)


      !----  Load K coefficients, interval 7  ----
      filename = trim(exort_rootdir)//trim(dirk_h2o)//trim(k07_h2o_file)
      call getfil(filename, locfn, 0)
      call wrap_open(locfn, 0, ncid)
      call wrap_inq_varid(ncid, 'data', keff_id)
      call wrap_get_var_realx(ncid, keff_id, k07_h2o)

      filename = trim(exort_rootdir)//trim(dirk_co2)//trim(k07_co2_file)
      call getfil(filename, locfn, 0)
      call wrap_open(locfn, 0, ncid)
      call wrap_inq_varid(ncid, 'data', keff_id)
      call wrap_get_var_realx(ncid, keff_id, k07_co2)

      filename = trim(exort_rootdir)//trim(dirk_ch4)//trim(k07_ch4_file)
      call getfil(filename, locfn, 0)
      call wrap_open(locfn, 0, ncid)
      call wrap_inq_varid(ncid, 'data', keff_id)
      call wrap_get_var_realx(ncid, keff_id, k07_ch4)


      !----  Load K coefficients, interval 8  ----
      filename = trim(exort_rootdir)//trim(dirk_h2o)//trim(k08_h2o_file)
      call getfil(filename, locfn, 0)
      call wrap_open(locfn, 0, ncid)
      call wrap_inq_varid(ncid, 'data', keff_id)
      call wrap_get_var_realx(ncid, keff_id, k08_h2o)

      filename = trim(exort_rootdir)//trim(dirk_co2)//trim(k08_co2_file)
      call getfil(filename, locfn, 0)
      call wrap_open(locfn, 0, ncid)
      call wrap_inq_varid(ncid, 'data', keff_id)
      call wrap_get_var_realx(ncid, keff_id, k08_co2)

      filename = trim(exort_rootdir)//trim(dirk_ch4)//trim(k08_ch4_file)
      call getfil(filename, locfn, 0)
      call wrap_open(locfn, 0, ncid)
      call wrap_inq_varid(ncid, 'data', keff_id)
      call wrap_get_var_realx(ncid, keff_id, k08_ch4)


      !----  Load K coefficients, interval 9  ----
      filename = trim(exort_rootdir)//trim(dirk_h2o)//trim(k09_h2o_file)
      call getfil(filename, locfn, 0)
      call wrap_open(locfn, 0, ncid)
      call wrap_inq_varid(ncid, 'data', keff_id)
      call wrap_get_var_realx(ncid, keff_id, k09_h2o)

      filename = trim(exort_rootdir)//trim(dirk_co2)//trim(k09_co2_file)
      call getfil(filename, locfn, 0)
      call wrap_open(locfn, 0, ncid)
      call wrap_inq_varid(ncid, 'data', keff_id)
      call wrap_get_var_realx(ncid, keff_id, k09_co2)

      filename = trim(exort_rootdir)//trim(dirk_ch4)//trim(k09_ch4_file)
      call getfil(filename, locfn, 0)
      call wrap_open(locfn, 0, ncid)
      call wrap_inq_varid(ncid, 'data', keff_id)
      call wrap_get_var_realx(ncid, keff_id, k09_ch4)


      !----  Load K coefficients, interval 10  ----
      filename = trim(exort_rootdir)//trim(dirk_h2o)//trim(k10_h2o_file)
      call getfil(filename, locfn, 0)
      call wrap_open(locfn, 0, ncid)
      call wrap_inq_varid(ncid, 'data', keff_id)
      call wrap_get_var_realx(ncid, keff_id, k10_h2o)

      filename = trim(exort_rootdir)//trim(dirk_co2)//trim(k10_co2_file)
      call getfil(filename, locfn, 0)
      call wrap_open(locfn, 0, ncid)
      call wrap_inq_varid(ncid, 'data', keff_id)
      call wrap_get_var_realx(ncid, keff_id, k10_co2)

      filename = trim(exort_rootdir)//trim(dirk_ch4)//trim(k10_ch4_file)
      call getfil(filename, locfn, 0)
      call wrap_open(locfn, 0, ncid)
      call wrap_inq_varid(ncid, 'data', keff_id)
      call wrap_get_var_realx(ncid, keff_id, k10_ch4)


      !----  Load K coefficients, interval 11  ----
      filename = trim(exort_rootdir)//trim(dirk_h2o)//trim(k11_h2o_file)
      call getfil(filename, locfn, 0)
      call wrap_open(locfn, 0, ncid)
      call wrap_inq_varid(ncid, 'data', keff_id)
      call wrap_get_var_realx(ncid, keff_id, k11_h2o)

      filename = trim(exort_rootdir)//trim(dirk_co2)//trim(k11_co2_file)
      call getfil(filename, locfn, 0)
      call wrap_open(locfn, 0, ncid)
      call wrap_inq_varid(ncid, 'data', keff_id)
      call wrap_get_var_realx(ncid, keff_id, k11_co2)

      filename = trim(exort_rootdir)//trim(dirk_ch4)//trim(k11_ch4_file)
      call getfil(filename, locfn, 0)
      call wrap_open(locfn, 0, ncid)
      call wrap_inq_varid(ncid, 'data', keff_id)
      call wrap_get_var_realx(ncid, keff_id, k11_ch4)

      !----  Load K coefficients, interval 12  ----
      filename = trim(exort_rootdir)//trim(dirk_h2o)//trim(k12_h2o_file)
      call getfil(filename, locfn, 0)
      call wrap_open(locfn, 0, ncid)
      call wrap_inq_varid(ncid, 'data', keff_id)
      call wrap_get_var_realx(ncid, keff_id, k12_h2o)

      filename = trim(exort_rootdir)//trim(dirk_co2)//trim(k12_co2_file)
      call getfil(filename, locfn, 0)
      call wrap_open(locfn, 0, ncid)
      call wrap_inq_varid(ncid, 'data', keff_id)
      call wrap_get_var_realx(ncid, keff_id, k12_co2)

      filename = trim(exort_rootdir)//trim(dirk_ch4)//trim(k12_ch4_file)
      call getfil(filename, locfn, 0)
      call wrap_open(locfn, 0, ncid)
      call wrap_inq_varid(ncid, 'data', keff_id)
      call wrap_get_var_realx(ncid, keff_id, k12_ch4)

      !----  Load K coefficients, interval 13  ----
      filename = trim(exort_rootdir)//trim(dirk_h2o)//trim(k13_h2o_file)
      call getfil(filename, locfn, 0)
      call wrap_open(locfn, 0, ncid)
      call wrap_inq_varid(ncid, 'data', keff_id)
      call wrap_get_var_realx(ncid, keff_id, k13_h2o)

      filename = trim(exort_rootdir)//trim(dirk_co2)//trim(k13_co2_file)
      call getfil(filename, locfn, 0)
      call wrap_open(locfn, 0, ncid)
      call wrap_inq_varid(ncid, 'data', keff_id)
      call wrap_get_var_realx(ncid, keff_id, k13_co2)

      filename = trim(exort_rootdir)//trim(dirk_ch4)//trim(k13_ch4_file)
      call getfil(filename, locfn, 0)
      call wrap_open(locfn, 0, ncid)
      call wrap_inq_varid(ncid, 'data', keff_id)
      call wrap_get_var_realx(ncid, keff_id, k13_ch4)


      !----  Load K coefficients, interval 14  ----
      filename = trim(exort_rootdir)//trim(dirk_h2o)//trim(k14_h2o_file)
      call getfil(filename, locfn, 0)
      call wrap_open(locfn, 0, ncid)
      call wrap_inq_varid(ncid, 'data', keff_id)
      call wrap_get_var_realx(ncid, keff_id, k14_h2o)

      filename = trim(exort_rootdir)//trim(dirk_co2)//trim(k14_co2_file)
      call getfil(filename, locfn, 0)
      call wrap_open(locfn, 0, ncid)
      call wrap_inq_varid(ncid, 'data', keff_id)
      call wrap_get_var_realx(ncid, keff_id, k14_co2)

      filename = trim(exort_rootdir)//trim(dirk_ch4)//trim(k14_ch4_file)
      call getfil(filename, locfn, 0)
      call wrap_open(locfn, 0, ncid)
      call wrap_inq_varid(ncid, 'data', keff_id)
      call wrap_get_var_realx(ncid, keff_id, k14_ch4)


      !----  Load K coefficients, interval 15  ----
      filename = trim(exort_rootdir)//trim(dirk_h2o)//trim(k15_h2o_file)
      call getfil(filename, locfn, 0)
      call wrap_open(locfn, 0, ncid)
      call wrap_inq_varid(ncid, 'data', keff_id)
      call wrap_get_var_realx(ncid, keff_id, k15_h2o)

      filename = trim(exort_rootdir)//trim(dirk_co2)//trim(k15_co2_file)
      call getfil(filename, locfn, 0)
      call wrap_open(locfn, 0, ncid)
      call wrap_inq_varid(ncid, 'data', keff_id)
      call wrap_get_var_realx(ncid, keff_id, k15_co2)

      filename = trim(exort_rootdir)//trim(dirk_ch4)//trim(k15_ch4_file)
      call getfil(filename, locfn, 0)
      call wrap_open(locfn, 0, ncid)
      call wrap_inq_varid(ncid, 'data', keff_id)
      call wrap_get_var_realx(ncid, keff_id, k15_ch4)


      !----  Load K coefficients, interval 16  ----
      filename = trim(exort_rootdir)//trim(dirk_h2o)//trim(k16_h2o_file)
      call getfil(filename, locfn, 0)
      call wrap_open(locfn, 0, ncid)
      call wrap_inq_varid(ncid, 'data', keff_id)
      call wrap_get_var_realx(ncid, keff_id, k16_h2o)

      filename = trim(exort_rootdir)//trim(dirk_co2)//trim(k16_co2_file)
      call getfil(filename, locfn, 0)
      call wrap_open(locfn, 0, ncid)
      call wrap_inq_varid(ncid, 'data', keff_id)
      call wrap_get_var_realx(ncid, keff_id, k16_co2)

      filename = trim(exort_rootdir)//trim(dirk_ch4)//trim(k16_ch4_file)
      call getfil(filename, locfn, 0)
      call wrap_open(locfn, 0, ncid)
      call wrap_inq_varid(ncid, 'data', keff_id)
      call wrap_get_var_realx(ncid, keff_id, k16_ch4)


      !----  Load K coefficients, interval 17  ----
      filename = trim(exort_rootdir)//trim(dirk_h2o)//trim(k17_h2o_file)
      call getfil(filename, locfn, 0)
      call wrap_open(locfn, 0, ncid)
      call wrap_inq_varid(ncid, 'data', keff_id)
      call wrap_get_var_realx(ncid, keff_id, k17_h2o)

      filename = trim(exort_rootdir)//trim(dirk_co2)//trim(k17_co2_file)
      call getfil(filename, locfn, 0)
      call wrap_open(locfn, 0, ncid)
      call wrap_inq_varid(ncid, 'data', keff_id)
      call wrap_get_var_realx(ncid, keff_id, k17_co2)

      filename = trim(exort_rootdir)//trim(dirk_ch4)//trim(k17_ch4_file)
      call getfil(filename, locfn, 0)
      call wrap_open(locfn, 0, ncid)
      call wrap_inq_varid(ncid, 'data', keff_id)
      call wrap_get_var_realx(ncid, keff_id, k17_ch4)


      !----  Load K coefficients, interval 18  ----
      filename = trim(exort_rootdir)//trim(dirk_h2o)//trim(k18_h2o_file)
      call getfil(filename, locfn, 0)
      call wrap_open(locfn, 0, ncid)
      call wrap_inq_varid(ncid, 'data', keff_id)
      call wrap_get_var_realx(ncid, keff_id, k18_h2o)

      filename = trim(exort_rootdir)//trim(dirk_co2)//trim(k18_co2_file)
      call getfil(filename, locfn, 0)
      call wrap_open(locfn, 0, ncid)
      call wrap_inq_varid(ncid, 'data', keff_id)
      call wrap_get_var_realx(ncid, keff_id, k18_co2)

      filename = trim(exort_rootdir)//trim(dirk_ch4)//trim(k18_ch4_file)
      call getfil(filename, locfn, 0)
      call wrap_open(locfn, 0, ncid)
      call wrap_inq_varid(ncid, 'data', keff_id)
      call wrap_get_var_realx(ncid, keff_id, k18_ch4)


      !----  Load K coefficients, interval 19  ----
      filename = trim(exort_rootdir)//trim(dirk_h2o)//trim(k19_h2o_file)
      call getfil(filename, locfn, 0)
      call wrap_open(locfn, 0, ncid)
      call wrap_inq_varid(ncid, 'data', keff_id)
      call wrap_get_var_realx(ncid, keff_id, k19_h2o)

      filename = trim(exort_rootdir)//trim(dirk_co2)//trim(k19_co2_file)
      call getfil(filename, locfn, 0)
      call wrap_open(locfn, 0, ncid)
      call wrap_inq_varid(ncid, 'data', keff_id)
      call wrap_get_var_realx(ncid, keff_id, k19_co2)

      filename = trim(exort_rootdir)//trim(dirk_ch4)//trim(k19_ch4_file)
      call getfil(filename, locfn, 0)
      call wrap_open(locfn, 0, ncid)
      call wrap_inq_varid(ncid, 'data', keff_id)
      call wrap_get_var_realx(ncid, keff_id, k19_ch4)


      !----  Load K coefficients, interval 20  ----
      filename = trim(exort_rootdir)//trim(dirk_h2o)//trim(k20_h2o_file)
      call getfil(filename, locfn, 0)
      call wrap_open(locfn, 0, ncid)
      call wrap_inq_varid(ncid, 'data', keff_id)
      call wrap_get_var_realx(ncid, keff_id, k20_h2o)

      filename = trim(exort_rootdir)//trim(dirk_co2)//trim(k20_co2_file)
      call getfil(filename, locfn, 0)
      call wrap_open(locfn, 0, ncid)
      call wrap_inq_varid(ncid, 'data', keff_id)
      call wrap_get_var_realx(ncid, keff_id, k20_co2)

      filename = trim(exort_rootdir)//trim(dirk_ch4)//trim(k20_ch4_file)
      call getfil(filename, locfn, 0)
      call wrap_open(locfn, 0, ncid)
      call wrap_inq_varid(ncid, 'data', keff_id)
      call wrap_get_var_realx(ncid, keff_id, k20_ch4)


      !----  Load K coefficients, interval 21  ----
      filename = trim(exort_rootdir)//trim(dirk_h2o)//trim(k21_h2o_file)
      call getfil(filename, locfn, 0)
      call wrap_open(locfn, 0, ncid)
      call wrap_inq_varid(ncid, 'data', keff_id)
      call wrap_get_var_realx(ncid, keff_id, k21_h2o)

      filename = trim(exort_rootdir)//trim(dirk_co2)//trim(k21_co2_file)
      call getfil(filename, locfn, 0)
      call wrap_open(locfn, 0, ncid)
      call wrap_inq_varid(ncid, 'data', keff_id)
      call wrap_get_var_realx(ncid, keff_id, k21_co2)

      filename = trim(exort_rootdir)//trim(dirk_ch4)//trim(k21_ch4_file)
      call getfil(filename, locfn, 0)
      call wrap_open(locfn, 0, ncid)
      call wrap_inq_varid(ncid, 'data', keff_id)
      call wrap_get_var_realx(ncid, keff_id, k21_ch4)


      !----  Load K coefficients, interval 22  ----
      filename = trim(exort_rootdir)//trim(dirk_h2o)//trim(k22_h2o_file)
      call getfil(filename, locfn, 0)
      call wrap_open(locfn, 0, ncid)
      call wrap_inq_varid(ncid, 'data', keff_id)
      call wrap_get_var_realx(ncid, keff_id, k22_h2o)

      filename = trim(exort_rootdir)//trim(dirk_co2)//trim(k22_co2_file)
      call getfil(filename, locfn, 0)
      call wrap_open(locfn, 0, ncid)
      call wrap_inq_varid(ncid, 'data', keff_id)
      call wrap_get_var_realx(ncid, keff_id, k22_co2)

      filename = trim(exort_rootdir)//trim(dirk_ch4)//trim(k22_ch4_file)
      call getfil(filename, locfn, 0)
      call wrap_open(locfn, 0, ncid)
      call wrap_inq_varid(ncid, 'data', keff_id)
      call wrap_get_var_realx(ncid, keff_id, k22_ch4)


      !----  Load K coefficients, interval 23  ----
      filename = trim(exort_rootdir)//trim(dirk_h2o)//trim(k23_h2o_file)
      call getfil(filename, locfn, 0)
      call wrap_open(locfn, 0, ncid)
      call wrap_inq_varid(ncid, 'data', keff_id)
      call wrap_get_var_realx(ncid, keff_id, k23_h2o)

      filename = trim(exort_rootdir)//trim(dirk_co2)//trim(k23_co2_file)
      call getfil(filename, locfn, 0)
      call wrap_open(locfn, 0, ncid)
      call wrap_inq_varid(ncid, 'data', keff_id)
      call wrap_get_var_realx(ncid, keff_id, k23_co2)

      filename = trim(exort_rootdir)//trim(dirk_ch4)//trim(k23_ch4_file)
      call getfil(filename, locfn, 0)
      call wrap_open(locfn, 0, ncid)
      call wrap_inq_varid(ncid, 'data', keff_id)
      call wrap_get_var_realx(ncid, keff_id, k23_ch4)


      !----  Load K coefficients, interval 24  ----
      filename = trim(exort_rootdir)//trim(dirk_h2o)//trim(k24_h2o_file)
      call getfil(filename, locfn, 0)
      call wrap_open(locfn, 0, ncid)
      call wrap_inq_varid(ncid, 'data', keff_id)
      call wrap_get_var_realx(ncid, keff_id, k24_h2o)

      filename = trim(exort_rootdir)//trim(dirk_co2)//trim(k24_co2_file)
      call getfil(filename, locfn, 0)
      call wrap_open(locfn, 0, ncid)
      call wrap_inq_varid(ncid, 'data', keff_id)
      call wrap_get_var_realx(ncid, keff_id, k24_co2)

      filename = trim(exort_rootdir)//trim(dirk_ch4)//trim(k24_ch4_file)
      call getfil(filename, locfn, 0)
      call wrap_open(locfn, 0, ncid)
      call wrap_inq_varid(ncid, 'data', keff_id)
      call wrap_get_var_realx(ncid, keff_id, k24_ch4)


      !----  Load K coefficients, interval 25  ----
      filename = trim(exort_rootdir)//trim(dirk_h2o)//trim(k25_h2o_file)
      call getfil(filename, locfn, 0)
      call wrap_open(locfn, 0, ncid)
      call wrap_inq_varid(ncid, 'data', keff_id)
      call wrap_get_var_realx(ncid, keff_id, k25_h2o)

      filename = trim(exort_rootdir)//trim(dirk_co2)//trim(k25_co2_file)
      call getfil(filename, locfn, 0)
      call wrap_open(locfn, 0, ncid)
      call wrap_inq_varid(ncid, 'data', keff_id)
      call wrap_get_var_realx(ncid, keff_id, k25_co2)

      filename = trim(exort_rootdir)//trim(dirk_ch4)//trim(k25_ch4_file)
      call getfil(filename, locfn, 0)
      call wrap_open(locfn, 0, ncid)
      call wrap_inq_varid(ncid, 'data', keff_id)
      call wrap_get_var_realx(ncid, keff_id, k25_ch4)


      !----  Load K coefficients, interval 26  ----
      filename = trim(exort_rootdir)//trim(dirk_h2o)//trim(k26_h2o_file)
      call getfil(filename, locfn, 0)
      call wrap_open(locfn, 0, ncid)
      call wrap_inq_varid(ncid, 'data', keff_id)
      call wrap_get_var_realx(ncid, keff_id, k26_h2o)

      filename = trim(exort_rootdir)//trim(dirk_co2)//trim(k26_co2_file)
      call getfil(filename, locfn, 0)
      call wrap_open(locfn, 0, ncid)
      call wrap_inq_varid(ncid, 'data', keff_id)
      call wrap_get_var_realx(ncid, keff_id, k26_co2)

      filename = trim(exort_rootdir)//trim(dirk_ch4)//trim(k26_ch4_file)
      call getfil(filename, locfn, 0)
      call wrap_open(locfn, 0, ncid)
      call wrap_inq_varid(ncid, 'data', keff_id)
      call wrap_get_var_realx(ncid, keff_id, k26_ch4)


      !----  Load K coefficients, interval 27  ----
      filename = trim(exort_rootdir)//trim(dirk_h2o)//trim(k27_h2o_file)
      call getfil(filename, locfn, 0)
      call wrap_open(locfn, 0, ncid)
      call wrap_inq_varid(ncid, 'data', keff_id)
      call wrap_get_var_realx(ncid, keff_id, k27_h2o)

      filename = trim(exort_rootdir)//trim(dirk_co2)//trim(k27_co2_file)
      call getfil(filename, locfn, 0)
      call wrap_open(locfn, 0, ncid)
      call wrap_inq_varid(ncid, 'data', keff_id)
      call wrap_get_var_realx(ncid, keff_id, k27_co2)

      filename = trim(exort_rootdir)//trim(dirk_ch4)//trim(k27_ch4_file)
      call getfil(filename, locfn, 0)
      call wrap_open(locfn, 0, ncid)
      call wrap_inq_varid(ncid, 'data', keff_id)
      call wrap_get_var_realx(ncid, keff_id, k27_ch4)


      !----  Load K coefficients, interval 28  ----
      filename = trim(exort_rootdir)//trim(dirk_h2o)//trim(k28_h2o_file)
      call getfil(filename, locfn, 0)
      call wrap_open(locfn, 0, ncid)
      call wrap_inq_varid(ncid, 'data', keff_id)
      call wrap_get_var_realx(ncid, keff_id, k28_h2o)

      filename = trim(exort_rootdir)//trim(dirk_co2)//trim(k28_co2_file)
      call getfil(filename, locfn, 0)
      call wrap_open(locfn, 0, ncid)
      call wrap_inq_varid(ncid, 'data', keff_id)
      call wrap_get_var_realx(ncid, keff_id, k28_co2)

      filename = trim(exort_rootdir)//trim(dirk_ch4)//trim(k28_ch4_file)
      call getfil(filename, locfn, 0)
      call wrap_open(locfn, 0, ncid)
      call wrap_inq_varid(ncid, 'data', keff_id)
      call wrap_get_var_realx(ncid, keff_id, k28_ch4)


      !----  Load K coefficients, interval 29  ----
      filename = trim(exort_rootdir)//trim(dirk_h2o)//trim(k29_h2o_file)
      call getfil(filename, locfn, 0)
      call wrap_open(locfn, 0, ncid)
      call wrap_inq_varid(ncid, 'data', keff_id)
      call wrap_get_var_realx(ncid, keff_id, k29_h2o)

      filename = trim(exort_rootdir)//trim(dirk_co2)//trim(k29_co2_file)
      call getfil(filename, locfn, 0)
      call wrap_open(locfn, 0, ncid)
      call wrap_inq_varid(ncid, 'data', keff_id)
      call wrap_get_var_realx(ncid, keff_id, k29_co2)

      filename = trim(exort_rootdir)//trim(dirk_ch4)//trim(k29_ch4_file)
      call getfil(filename, locfn, 0)
      call wrap_open(locfn, 0, ncid)
      call wrap_inq_varid(ncid, 'data', keff_id)
      call wrap_get_var_realx(ncid, keff_id, k29_ch4)


      !----  Load K coefficients, interval 30  ----
      filename = trim(exort_rootdir)//trim(dirk_h2o)//trim(k30_h2o_file)
      call getfil(filename, locfn, 0)
      call wrap_open(locfn, 0, ncid)
      call wrap_inq_varid(ncid, 'data', keff_id)
      call wrap_get_var_realx(ncid, keff_id, k30_h2o)

      filename = trim(exort_rootdir)//trim(dirk_co2)//trim(k30_co2_file)
      call getfil(filename, locfn, 0)
      call wrap_open(locfn, 0, ncid)
      call wrap_inq_varid(ncid, 'data', keff_id)
      call wrap_get_var_realx(ncid, keff_id, k30_co2)

      filename = trim(exort_rootdir)//trim(dirk_ch4)//trim(k30_ch4_file)
      call getfil(filename, locfn, 0)
      call wrap_open(locfn, 0, ncid)
      call wrap_inq_varid(ncid, 'data', keff_id)
      call wrap_get_var_realx(ncid, keff_id, k30_ch4)


      !----  Load K coefficients, interval 31  ----
      filename = trim(exort_rootdir)//trim(dirk_h2o)//trim(k31_h2o_file)
      call getfil(filename, locfn, 0)
      call wrap_open(locfn, 0, ncid)
      call wrap_inq_varid(ncid, 'data', keff_id)
      call wrap_get_var_realx(ncid, keff_id, k31_h2o)

      filename = trim(exort_rootdir)//trim(dirk_co2)//trim(k31_co2_file)
      call getfil(filename, locfn, 0)
      call wrap_open(locfn, 0, ncid)
      call wrap_inq_varid(ncid, 'data', keff_id)
      call wrap_get_var_realx(ncid, keff_id, k31_co2)

      filename = trim(exort_rootdir)//trim(dirk_ch4)//trim(k31_ch4_file)
      call getfil(filename, locfn, 0)
      call wrap_open(locfn, 0, ncid)
      call wrap_inq_varid(ncid, 'data', keff_id)
      call wrap_get_var_realx(ncid, keff_id, k31_ch4)


      !----  Load K coefficients, interval 32  ----
      filename = trim(exort_rootdir)//trim(dirk_h2o)//trim(k32_h2o_file)
      call getfil(filename, locfn, 0)
      call wrap_open(locfn, 0, ncid)
      call wrap_inq_varid(ncid, 'data', keff_id)
      call wrap_get_var_realx(ncid, keff_id, k32_h2o)

      filename = trim(exort_rootdir)//trim(dirk_co2)//trim(k32_co2_file)
      call getfil(filename, locfn, 0)
      call wrap_open(locfn, 0, ncid)
      call wrap_inq_varid(ncid, 'data', keff_id)
      call wrap_get_var_realx(ncid, keff_id, k32_co2)

      filename = trim(exort_rootdir)//trim(dirk_ch4)//trim(k32_ch4_file)
      call getfil(filename, locfn, 0)
      call wrap_open(locfn, 0, ncid)
      call wrap_inq_varid(ncid, 'data', keff_id)
      call wrap_get_var_realx(ncid, keff_id, k32_ch4)


      !----  Load K coefficients, interval 33  ----
      filename = trim(exort_rootdir)//trim(dirk_h2o)//trim(k33_h2o_file)
      call getfil(filename, locfn, 0)
      call wrap_open(locfn, 0, ncid)
      call wrap_inq_varid(ncid, 'data', keff_id)
      call wrap_get_var_realx(ncid, keff_id, k33_h2o)

      filename = trim(exort_rootdir)//trim(dirk_co2)//trim(k33_co2_file)
      call getfil(filename, locfn, 0)
      call wrap_open(locfn, 0, ncid)
      call wrap_inq_varid(ncid, 'data', keff_id)
      call wrap_get_var_realx(ncid, keff_id, k33_co2)

      filename = trim(exort_rootdir)//trim(dirk_ch4)//trim(k33_ch4_file)
      call getfil(filename, locfn, 0)
      call wrap_open(locfn, 0, ncid)
      call wrap_inq_varid(ncid, 'data', keff_id)
      call wrap_get_var_realx(ncid, keff_id, k33_ch4)


      !----  Load K coefficients, interval 34  ----
      filename = trim(exort_rootdir)//trim(dirk_h2o)//trim(k34_h2o_file)
      call getfil(filename, locfn, 0)
      call wrap_open(locfn, 0, ncid)
      call wrap_inq_varid(ncid, 'data', keff_id)
      call wrap_get_var_realx(ncid, keff_id, k34_h2o)

      filename = trim(exort_rootdir)//trim(dirk_co2)//trim(k34_co2_file)
      call getfil(filename, locfn, 0)
      call wrap_open(locfn, 0, ncid)
      call wrap_inq_varid(ncid, 'data', keff_id)
      call wrap_get_var_realx(ncid, keff_id, k34_co2)

      filename = trim(exort_rootdir)//trim(dirk_ch4)//trim(k34_ch4_file)
      call getfil(filename, locfn, 0)
      call wrap_open(locfn, 0, ncid)
      call wrap_inq_varid(ncid, 'data', keff_id)
      call wrap_get_var_realx(ncid, keff_id, k34_ch4)


      !----  Load K coefficients, interval 35  ----
      filename = trim(exort_rootdir)//trim(dirk_h2o)//trim(k35_h2o_file)
      call getfil(filename, locfn, 0)
      call wrap_open(locfn, 0, ncid)
      call wrap_inq_varid(ncid, 'data', keff_id)
      call wrap_get_var_realx(ncid, keff_id, k35_h2o)

      filename = trim(exort_rootdir)//trim(dirk_co2)//trim(k35_co2_file)
      call getfil(filename, locfn, 0)
      call wrap_open(locfn, 0, ncid)
      call wrap_inq_varid(ncid, 'data', keff_id)
      call wrap_get_var_realx(ncid, keff_id, k35_co2)

      filename = trim(exort_rootdir)//trim(dirk_ch4)//trim(k35_ch4_file)
      call getfil(filename, locfn, 0)
      call wrap_open(locfn, 0, ncid)
      call wrap_inq_varid(ncid, 'data', keff_id)
      call wrap_get_var_realx(ncid, keff_id, k35_ch4)


      !----  Load K coefficients, interval 36  ----
      filename = trim(exort_rootdir)//trim(dirk_h2o)//trim(k36_h2o_file)
      call getfil(filename, locfn, 0)
      call wrap_open(locfn, 0, ncid)
      call wrap_inq_varid(ncid, 'data', keff_id)
      call wrap_get_var_realx(ncid, keff_id, k36_h2o)

      filename = trim(exort_rootdir)//trim(dirk_co2)//trim(k36_co2_file)
      call getfil(filename, locfn, 0)
      call wrap_open(locfn, 0, ncid)
      call wrap_inq_varid(ncid, 'data', keff_id)
      call wrap_get_var_realx(ncid, keff_id, k36_co2)

      filename = trim(exort_rootdir)//trim(dirk_ch4)//trim(k36_ch4_file)
      call getfil(filename, locfn, 0)
      call wrap_open(locfn, 0, ncid)
      call wrap_inq_varid(ncid, 'data', keff_id)
      call wrap_get_var_realx(ncid, keff_id, k36_ch4)


      !----  Load K coefficients, interval 37  ----
      filename = trim(exort_rootdir)//trim(dirk_h2o)//trim(k37_h2o_file)
      call getfil(filename, locfn, 0)
      call wrap_open(locfn, 0, ncid)
      call wrap_inq_varid(ncid, 'data', keff_id)
      call wrap_get_var_realx(ncid, keff_id, k37_h2o)

      filename = trim(exort_rootdir)//trim(dirk_co2)//trim(k37_co2_file)
      call getfil(filename, locfn, 0)
      call wrap_open(locfn, 0, ncid)
      call wrap_inq_varid(ncid, 'data', keff_id)
      call wrap_get_var_realx(ncid, keff_id, k37_co2)

      filename = trim(exort_rootdir)//trim(dirk_ch4)//trim(k37_ch4_file)
      call getfil(filename, locfn, 0)
      call wrap_open(locfn, 0, ncid)
      call wrap_inq_varid(ncid, 'data', keff_id)
      call wrap_get_var_realx(ncid, keff_id, k37_ch4)


      !----  Load K coefficients, interval 38  ----
      filename = trim(exort_rootdir)//trim(dirk_h2o)//trim(k38_h2o_file)
      call getfil(filename, locfn, 0)
      call wrap_open(locfn, 0, ncid)
      call wrap_inq_varid(ncid, 'data', keff_id)
      call wrap_get_var_realx(ncid, keff_id, k38_h2o)

      filename = trim(exort_rootdir)//trim(dirk_co2)//trim(k38_co2_file)
      call getfil(filename, locfn, 0)
      call wrap_open(locfn, 0, ncid)
      call wrap_inq_varid(ncid, 'data', keff_id)
      call wrap_get_var_realx(ncid, keff_id, k38_co2)

      filename = trim(exort_rootdir)//trim(dirk_ch4)//trim(k38_ch4_file)
      call getfil(filename, locfn, 0)
      call wrap_open(locfn, 0, ncid)
      call wrap_inq_varid(ncid, 'data', keff_id)
      call wrap_get_var_realx(ncid, keff_id, k38_ch4)


      !----  Load K coefficients, interval 39  ----
      filename = trim(exort_rootdir)//trim(dirk_h2o)//trim(k39_h2o_file)
      call getfil(filename, locfn, 0)
      call wrap_open(locfn, 0, ncid)
      call wrap_inq_varid(ncid, 'data', keff_id)
      call wrap_get_var_realx(ncid, keff_id, k39_h2o)

      filename = trim(exort_rootdir)//trim(dirk_co2)//trim(k39_co2_file)
      call getfil(filename, locfn, 0)
      call wrap_open(locfn, 0, ncid)
      call wrap_inq_varid(ncid, 'data', keff_id)
      call wrap_get_var_realx(ncid, keff_id, k39_co2)

      filename = trim(exort_rootdir)//trim(dirk_ch4)//trim(k39_ch4_file)
      call getfil(filename, locfn, 0)
      call wrap_open(locfn, 0, ncid)
      call wrap_inq_varid(ncid, 'data', keff_id)
      call wrap_get_var_realx(ncid, keff_id, k39_ch4)


      !----  Load K coefficients, interval 40  ----
      filename = trim(exort_rootdir)//trim(dirk_h2o)//trim(k40_h2o_file)
      call getfil(filename, locfn, 0)
      call wrap_open(locfn, 0, ncid)
      call wrap_inq_varid(ncid, 'data', keff_id)
      call wrap_get_var_realx(ncid, keff_id, k40_h2o)

      filename = trim(exort_rootdir)//trim(dirk_co2)//trim(k40_co2_file)
      call getfil(filename, locfn, 0)
      call wrap_open(locfn, 0, ncid)
      call wrap_inq_varid(ncid, 'data', keff_id)
      call wrap_get_var_realx(ncid, keff_id, k40_co2)

      filename = trim(exort_rootdir)//trim(dirk_ch4)//trim(k40_ch4_file)
      call getfil(filename, locfn, 0)
      call wrap_open(locfn, 0, ncid)
      call wrap_inq_varid(ncid, 'data', keff_id)
      call wrap_get_var_realx(ncid, keff_id, k40_ch4)


      !----  Load K coefficients, interval 41  ----
      filename = trim(exort_rootdir)//trim(dirk_h2o)//trim(k41_h2o_file)
      call getfil(filename, locfn, 0)
      call wrap_open(locfn, 0, ncid)
      call wrap_inq_varid(ncid, 'data', keff_id)
      call wrap_get_var_realx(ncid, keff_id, k41_h2o)

      filename = trim(exort_rootdir)//trim(dirk_co2)//trim(k41_co2_file)
      call getfil(filename, locfn, 0)
      call wrap_open(locfn, 0, ncid)
      call wrap_inq_varid(ncid, 'data', keff_id)
      call wrap_get_var_realx(ncid, keff_id, k41_co2)

      filename = trim(exort_rootdir)//trim(dirk_ch4)//trim(k41_ch4_file)
      call getfil(filename, locfn, 0)
      call wrap_open(locfn, 0, ncid)
      call wrap_inq_varid(ncid, 'data', keff_id)
      call wrap_get_var_realx(ncid, keff_id, k41_ch4)


      !----  Load K coefficients, interval 42  ----  
      filename = trim(exort_rootdir)//trim(dirk_h2o)//trim(k42_h2o_file)
      call getfil(filename, locfn, 0)
      call wrap_open(locfn, 0, ncid)
      call wrap_inq_varid(ncid, 'data', keff_id)
      call wrap_get_var_realx(ncid, keff_id, k42_h2o)

      filename = trim(exort_rootdir)//trim(dirk_co2)//trim(k42_co2_file)
      call getfil(filename, locfn, 0)
      call wrap_open(locfn, 0, ncid)
      call wrap_inq_varid(ncid, 'data', keff_id)
      call wrap_get_var_realx(ncid, keff_id, k42_co2)

      filename = trim(exort_rootdir)//trim(dirk_ch4)//trim(k42_ch4_file)
      call getfil(filename, locfn, 0)
      call wrap_open(locfn, 0, ncid)
      call wrap_inq_varid(ncid, 'data', keff_id)
      call wrap_get_var_realx(ncid, keff_id, k42_ch4)


      !----  Load K coefficients, interval 43  ----
      filename = trim(exort_rootdir)//trim(dirk_h2o)//trim(k43_h2o_file)
      call getfil(filename, locfn, 0)
      call wrap_open(locfn, 0, ncid)
      call wrap_inq_varid(ncid, 'data', keff_id)
      call wrap_get_var_realx(ncid, keff_id, k43_h2o)

      filename = trim(exort_rootdir)//trim(dirk_co2)//trim(k43_co2_file)
      call getfil(filename, locfn, 0)
      call wrap_open(locfn, 0, ncid)
      call wrap_inq_varid(ncid, 'data', keff_id)
      call wrap_get_var_realx(ncid, keff_id, k43_co2)

      filename = trim(exort_rootdir)//trim(dirk_ch4)//trim(k43_ch4_file)
      call getfil(filename, locfn, 0)
      call wrap_open(locfn, 0, ncid)
      call wrap_inq_varid(ncid, 'data', keff_id)
      call wrap_get_var_realx(ncid, keff_id, k43_ch4)


      !----  Load K coefficients, interval 44  ----
      filename = trim(exort_rootdir)//trim(dirk_h2o)//trim(k44_h2o_file)
      call getfil(filename, locfn, 0)
      call wrap_open(locfn, 0, ncid)
      call wrap_inq_varid(ncid, 'data', keff_id)
      call wrap_get_var_realx(ncid, keff_id, k44_h2o)

      filename = trim(exort_rootdir)//trim(dirk_co2)//trim(k44_co2_file)
      call getfil(filename, locfn, 0)
      call wrap_open(locfn, 0, ncid)
      call wrap_inq_varid(ncid, 'data', keff_id)
      call wrap_get_var_realx(ncid, keff_id, k44_co2)

      filename = trim(exort_rootdir)//trim(dirk_ch4)//trim(k44_ch4_file)
      call getfil(filename, locfn, 0)
      call wrap_open(locfn, 0, ncid)
      call wrap_inq_varid(ncid, 'data', keff_id)
      call wrap_get_var_realx(ncid, keff_id, k44_ch4)


      !----  Load K coefficients, interval 45  ----
      filename = trim(exort_rootdir)//trim(dirk_h2o)//trim(k45_h2o_file)
      call getfil(filename, locfn, 0)
      call wrap_open(locfn, 0, ncid)
      call wrap_inq_varid(ncid, 'data', keff_id)
      call wrap_get_var_realx(ncid, keff_id, k45_h2o)

      filename = trim(exort_rootdir)//trim(dirk_co2)//trim(k45_co2_file)
      call getfil(filename, locfn, 0)
      call wrap_open(locfn, 0, ncid)
      call wrap_inq_varid(ncid, 'data', keff_id)
      call wrap_get_var_realx(ncid, keff_id, k45_co2)

      filename = trim(exort_rootdir)//trim(dirk_ch4)//trim(k45_ch4_file)
      call getfil(filename, locfn, 0)
      call wrap_open(locfn, 0, ncid)
      call wrap_inq_varid(ncid, 'data', keff_id)
      call wrap_get_var_realx(ncid, keff_id, k45_ch4)


      !----  Load K coefficients, interval 46  ----
      filename = trim(exort_rootdir)//trim(dirk_h2o)//trim(k46_h2o_file)
      call getfil(filename, locfn, 0)
      call wrap_open(locfn, 0, ncid)
      call wrap_inq_varid(ncid, 'data', keff_id)
      call wrap_get_var_realx(ncid, keff_id, k46_h2o)

      filename = trim(exort_rootdir)//trim(dirk_co2)//trim(k46_co2_file)
      call getfil(filename, locfn, 0)
      call wrap_open(locfn, 0, ncid)
      call wrap_inq_varid(ncid, 'data', keff_id)
      call wrap_get_var_realx(ncid, keff_id, k46_co2)

      filename = trim(exort_rootdir)//trim(dirk_ch4)//trim(k46_ch4_file)
      call getfil(filename, locfn, 0)
      call wrap_open(locfn, 0, ncid)
      call wrap_inq_varid(ncid, 'data', keff_id)
      call wrap_get_var_realx(ncid, keff_id, k46_ch4)


      !----  Load K coefficients, interval 47  ----
      filename = trim(exort_rootdir)//trim(dirk_h2o)//trim(k47_h2o_file)
      call getfil(filename, locfn, 0)
      call wrap_open(locfn, 0, ncid)
      call wrap_inq_varid(ncid, 'data', keff_id)
      call wrap_get_var_realx(ncid, keff_id, k47_h2o)

      filename = trim(exort_rootdir)//trim(dirk_co2)//trim(k47_co2_file)
      call getfil(filename, locfn, 0)
      call wrap_open(locfn, 0, ncid)
      call wrap_inq_varid(ncid, 'data', keff_id)
      call wrap_get_var_realx(ncid, keff_id, k47_co2)

      filename = trim(exort_rootdir)//trim(dirk_ch4)//trim(k47_ch4_file)
      call getfil(filename, locfn, 0)
      call wrap_open(locfn, 0, ncid)
      call wrap_inq_varid(ncid, 'data', keff_id)
      call wrap_get_var_realx(ncid, keff_id, k47_ch4)


      !----  Load K coefficients, interval 48  -----
      filename = trim(exort_rootdir)//trim(dirk_h2o)//trim(k48_h2o_file)
      call getfil(filename, locfn, 0)
      call wrap_open(locfn, 0, ncid)
      call wrap_inq_varid(ncid, 'data', keff_id)
      call wrap_get_var_realx(ncid, keff_id, k48_h2o)

      filename = trim(exort_rootdir)//trim(dirk_co2)//trim(k48_co2_file)
      call getfil(filename, locfn, 0)
      call wrap_open(locfn, 0, ncid)
      call wrap_inq_varid(ncid, 'data', keff_id)
      call wrap_get_var_realx(ncid, keff_id, k48_co2)

      filename = trim(exort_rootdir)//trim(dirk_ch4)//trim(k48_ch4_file)
      call getfil(filename, locfn, 0)
      call wrap_open(locfn, 0, ncid)
      call wrap_inq_varid(ncid, 'data', keff_id)
      call wrap_get_var_realx(ncid, keff_id, k48_ch4)


      !----  Load K coefficients, interval 49  ----
      filename = trim(exort_rootdir)//trim(dirk_h2o)//trim(k49_h2o_file)
      call getfil(filename, locfn, 0)
      call wrap_open(locfn, 0, ncid)
      call wrap_inq_varid(ncid, 'data', keff_id)
      call wrap_get_var_realx(ncid, keff_id, k49_h2o)

      filename = trim(exort_rootdir)//trim(dirk_co2)//trim(k49_co2_file)
      call getfil(filename, locfn, 0)
      call wrap_open(locfn, 0, ncid)
      call wrap_inq_varid(ncid, 'data', keff_id)
      call wrap_get_var_realx(ncid, keff_id, k49_co2)

      filename = trim(exort_rootdir)//trim(dirk_ch4)//trim(k49_ch4_file)
      call getfil(filename, locfn, 0)
      call wrap_open(locfn, 0, ncid)
      call wrap_inq_varid(ncid, 'data', keff_id)
      call wrap_get_var_realx(ncid, keff_id, k49_ch4)


      !----  Load K coefficients, interval 50  ----
      filename = trim(exort_rootdir)//trim(dirk_h2o)//trim(k50_h2o_file)
      call getfil(filename, locfn, 0)
      call wrap_open(locfn, 0, ncid)
      call wrap_inq_varid(ncid, 'data', keff_id)
      call wrap_get_var_realx(ncid, keff_id, k50_h2o)

      filename = trim(exort_rootdir)//trim(dirk_co2)//trim(k50_co2_file)
      call getfil(filename, locfn, 0)
      call wrap_open(locfn, 0, ncid)
      call wrap_inq_varid(ncid, 'data', keff_id)
      call wrap_get_var_realx(ncid, keff_id, k50_co2)

      filename = trim(exort_rootdir)//trim(dirk_ch4)//trim(k50_ch4_file)
      call getfil(filename, locfn, 0)
      call wrap_open(locfn, 0, ncid)
      call wrap_inq_varid(ncid, 'data', keff_id)
      call wrap_get_var_realx(ncid, keff_id, k50_ch4)


      !----  Load K coefficients, interval 51  ----
      filename = trim(exort_rootdir)//trim(dirk_h2o)//trim(k51_h2o_file)
      call getfil(filename, locfn, 0)
      call wrap_open(locfn, 0, ncid)
      call wrap_inq_varid(ncid, 'data', keff_id)
      call wrap_get_var_realx(ncid, keff_id, k51_h2o)

      filename = trim(exort_rootdir)//trim(dirk_co2)//trim(k51_co2_file)
      call getfil(filename, locfn, 0)
      call wrap_open(locfn, 0, ncid)
      call wrap_inq_varid(ncid, 'data', keff_id)
      call wrap_get_var_realx(ncid, keff_id, k51_co2)

      filename = trim(exort_rootdir)//trim(dirk_ch4)//trim(k51_ch4_file)
      call getfil(filename, locfn, 0)
      call wrap_open(locfn, 0, ncid)
      call wrap_inq_varid(ncid, 'data', keff_id)
      call wrap_get_var_realx(ncid, keff_id, k51_ch4)


      !----  Load K coefficients, interval 52  ----
      filename = trim(exort_rootdir)//trim(dirk_h2o)//trim(k52_h2o_file)
      call getfil(filename, locfn, 0)
      call wrap_open(locfn, 0, ncid)
      call wrap_inq_varid(ncid, 'data', keff_id)
      call wrap_get_var_realx(ncid, keff_id, k52_h2o)

      filename = trim(exort_rootdir)//trim(dirk_co2)//trim(k52_co2_file)
      call getfil(filename, locfn, 0)
      call wrap_open(locfn, 0, ncid)
      call wrap_inq_varid(ncid, 'data', keff_id)
      call wrap_get_var_realx(ncid, keff_id, k52_co2)

      filename = trim(exort_rootdir)//trim(dirk_ch4)//trim(k52_ch4_file)
      call getfil(filename, locfn, 0)
      call wrap_open(locfn, 0, ncid)
      call wrap_inq_varid(ncid, 'data', keff_id)
      call wrap_get_var_realx(ncid, keff_id, k52_ch4)


      !----  Load K coefficients, interval 53  ----
      filename = trim(exort_rootdir)//trim(dirk_h2o)//trim(k53_h2o_file)
      call getfil(filename, locfn, 0)
      call wrap_open(locfn, 0, ncid)
      call wrap_inq_varid(ncid, 'data', keff_id)
      call wrap_get_var_realx(ncid, keff_id, k53_h2o)

      filename = trim(exort_rootdir)//trim(dirk_co2)//trim(k53_co2_file)
      call getfil(filename, locfn, 0)
      call wrap_open(locfn, 0, ncid)
      call wrap_inq_varid(ncid, 'data', keff_id)
      call wrap_get_var_realx(ncid, keff_id, k53_co2)

      filename = trim(exort_rootdir)//trim(dirk_ch4)//trim(k54_ch4_file)
      call getfil(filename, locfn, 0)
      call wrap_open(locfn, 0, ncid)
      call wrap_inq_varid(ncid, 'data', keff_id)
      call wrap_get_var_realx(ncid, keff_id, k54_ch4)


      !----  Load K coefficients, interval 54  ----
      filename = trim(exort_rootdir)//trim(dirk_h2o)//trim(k54_h2o_file)
      call getfil(filename, locfn, 0)
      call wrap_open(locfn, 0, ncid)
      call wrap_inq_varid(ncid, 'data', keff_id)
      call wrap_get_var_realx(ncid, keff_id, k54_h2o)

      filename = trim(exort_rootdir)//trim(dirk_co2)//trim(k54_co2_file)
      call getfil(filename, locfn, 0)
      call wrap_open(locfn, 0, ncid)
      call wrap_inq_varid(ncid, 'data', keff_id)
      call wrap_get_var_realx(ncid, keff_id, k54_co2)

      filename = trim(exort_rootdir)//trim(dirk_ch4)//trim(k54_ch4_file)
      call getfil(filename, locfn, 0)
      call wrap_open(locfn, 0, ncid)
      call wrap_inq_varid(ncid, 'data', keff_id)
      call wrap_get_var_realx(ncid, keff_id, k54_ch4)


      !----  Load K coefficients, interval 55  ----
      filename = trim(exort_rootdir)//trim(dirk_h2o)//trim(k55_h2o_file)
      call getfil(filename, locfn, 0)
      call wrap_open(locfn, 0, ncid)
      call wrap_inq_varid(ncid, 'data', keff_id)
      call wrap_get_var_realx(ncid, keff_id, k55_h2o)

      filename = trim(exort_rootdir)//trim(dirk_co2)//trim(k55_co2_file)
      call getfil(filename, locfn, 0)
      call wrap_open(locfn, 0, ncid)
      call wrap_inq_varid(ncid, 'data', keff_id)
      call wrap_get_var_realx(ncid, keff_id, k55_co2)

      filename = trim(exort_rootdir)//trim(dirk_ch4)//trim(k55_ch4_file)
      call getfil(filename, locfn, 0)
      call wrap_open(locfn, 0, ncid)
      call wrap_inq_varid(ncid, 'data', keff_id)
      call wrap_get_var_realx(ncid, keff_id, k55_ch4)


      !----  Load K coefficients, interval 56  ----
      filename = trim(exort_rootdir)//trim(dirk_h2o)//trim(k56_h2o_file)
      call getfil(filename, locfn, 0)
      call wrap_open(locfn, 0, ncid)
      call wrap_inq_varid(ncid, 'data', keff_id)
      call wrap_get_var_realx(ncid, keff_id, k56_h2o)

      filename = trim(exort_rootdir)//trim(dirk_co2)//trim(k56_co2_file)
      call getfil(filename, locfn, 0)
      call wrap_open(locfn, 0, ncid)
      call wrap_inq_varid(ncid, 'data', keff_id)
      call wrap_get_var_realx(ncid, keff_id, k56_co2)

      filename = trim(exort_rootdir)//trim(dirk_ch4)//trim(k56_ch4_file)
      call getfil(filename, locfn, 0)
      call wrap_open(locfn, 0, ncid)
      call wrap_inq_varid(ncid, 'data', keff_id)
      call wrap_get_var_realx(ncid, keff_id, k56_ch4)


      !----  Load K coefficients, interval 57  ----
      filename = trim(exort_rootdir)//trim(dirk_h2o)//trim(k57_h2o_file)
      call getfil(filename, locfn, 0)
      call wrap_open(locfn, 0, ncid)
      call wrap_inq_varid(ncid, 'data', keff_id)
      call wrap_get_var_realx(ncid, keff_id, k57_h2o)

      filename = trim(exort_rootdir)//trim(dirk_co2)//trim(k57_co2_file)
      call getfil(filename, locfn, 0)
      call wrap_open(locfn, 0, ncid)
      call wrap_inq_varid(ncid, 'data', keff_id)
      call wrap_get_var_realx(ncid, keff_id, k57_co2)

      filename = trim(exort_rootdir)//trim(dirk_ch4)//trim(k57_ch4_file)
      call getfil(filename, locfn, 0)
      call wrap_open(locfn, 0, ncid)
      call wrap_inq_varid(ncid, 'data', keff_id)
      call wrap_get_var_realx(ncid, keff_id, k57_ch4)


      ! Load K coefficients, interval 58  
      filename = trim(exort_rootdir)//trim(dirk_h2o)//trim(k58_h2o_file)
      call getfil(filename, locfn, 0)
      call wrap_open(locfn, 0, ncid)
      call wrap_inq_varid(ncid, 'data', keff_id)
      call wrap_get_var_realx(ncid, keff_id, k58_h2o)

      filename = trim(exort_rootdir)//trim(dirk_co2)//trim(k58_co2_file)
      call getfil(filename, locfn, 0)
      call wrap_open(locfn, 0, ncid)
      call wrap_inq_varid(ncid, 'data', keff_id)
      call wrap_get_var_realx(ncid, keff_id, k58_co2)

      filename = trim(exort_rootdir)//trim(dirk_ch4)//trim(k58_ch4_file)
      call getfil(filename, locfn, 0)
      call wrap_open(locfn, 0, ncid)
      call wrap_inq_varid(ncid, 'data', keff_id)
      call wrap_get_var_realx(ncid, keff_id, k58_ch4)


      ! Load K coefficients, interval 59  
      filename = trim(exort_rootdir)//trim(dirk_h2o)//trim(k59_h2o_file)
      call getfil(filename, locfn, 0)
      call wrap_open(locfn, 0, ncid)
      call wrap_inq_varid(ncid, 'data', keff_id)
      call wrap_get_var_realx(ncid, keff_id, k59_h2o)

      filename = trim(exort_rootdir)//trim(dirk_co2)//trim(k59_co2_file)
      call getfil(filename, locfn, 0)
      call wrap_open(locfn, 0, ncid)
      call wrap_inq_varid(ncid, 'data', keff_id)
      call wrap_get_var_realx(ncid, keff_id, k59_co2)

      filename = trim(exort_rootdir)//trim(dirk_ch4)//trim(k59_ch4_file)
      call getfil(filename, locfn, 0)
      call wrap_open(locfn, 0, ncid)
      call wrap_inq_varid(ncid, 'data', keff_id)
      call wrap_get_var_realx(ncid, keff_id, k59_ch4)


      ! Load K coefficients, interval 60  
      filename = trim(exort_rootdir)//trim(dirk_h2o)//trim(k60_h2o_file)
      call getfil(filename, locfn, 0)
      call wrap_open(locfn, 0, ncid)
      call wrap_inq_varid(ncid, 'data', keff_id)
      call wrap_get_var_realx(ncid, keff_id, k60_h2o)

      filename = trim(exort_rootdir)//trim(dirk_co2)//trim(k60_co2_file)
      call getfil(filename, locfn, 0)
      call wrap_open(locfn, 0, ncid)
      call wrap_inq_varid(ncid, 'data', keff_id)
      call wrap_get_var_realx(ncid, keff_id, k60_co2)

      filename = trim(exort_rootdir)//trim(dirk_ch4)//trim(k60_ch4_file)
      call getfil(filename, locfn, 0)
      call wrap_open(locfn, 0, ncid)
      call wrap_inq_varid(ncid, 'data', keff_id)
      call wrap_get_var_realx(ncid, keff_id, k60_ch4)


      !----  Load K coefficients, interval 61  ----
      filename = trim(exort_rootdir)//trim(dirk_h2o)//trim(k61_h2o_file)
      call getfil(filename, locfn, 0)
      call wrap_open(locfn, 0, ncid)
      call wrap_inq_varid(ncid, 'data', keff_id)
      call wrap_get_var_realx(ncid, keff_id, k61_h2o)

      filename = trim(exort_rootdir)//trim(dirk_co2)//trim(k61_co2_file)
      call getfil(filename, locfn, 0)
      call wrap_open(locfn, 0, ncid)
      call wrap_inq_varid(ncid, 'data', keff_id)
      call wrap_get_var_realx(ncid, keff_id, k61_co2)

      filename = trim(exort_rootdir)//trim(dirk_ch4)//trim(k61_ch4_file)
      call getfil(filename, locfn, 0)
      call wrap_open(locfn, 0, ncid)
      call wrap_inq_varid(ncid, 'data', keff_id)
      call wrap_get_var_realx(ncid, keff_id, k61_ch4)


      !----  Load K coefficients, interval 62  -----
      filename = trim(exort_rootdir)//trim(dirk_h2o)//trim(k62_h2o_file)
      call getfil(filename, locfn, 0)
      call wrap_open(locfn, 0, ncid)
      call wrap_inq_varid(ncid, 'data', keff_id)
      call wrap_get_var_realx(ncid, keff_id, k62_h2o)

      filename = trim(exort_rootdir)//trim(dirk_co2)//trim(k62_co2_file)
      call getfil(filename, locfn, 0)
      call wrap_open(locfn, 0, ncid)
      call wrap_inq_varid(ncid, 'data', keff_id)
      call wrap_get_var_realx(ncid, keff_id, k62_co2)

      filename = trim(exort_rootdir)//trim(dirk_ch4)//trim(k62_ch4_file)
      call getfil(filename, locfn, 0)
      call wrap_open(locfn, 0, ncid)
      call wrap_inq_varid(ncid, 'data', keff_id)
      call wrap_get_var_realx(ncid, keff_id, k62_ch4)


      !---- Load K coefficients, interval 63  ----
      filename = trim(exort_rootdir)//trim(dirk_h2o)//trim(k63_h2o_file)
      call getfil(filename, locfn, 0)
      call wrap_open(locfn, 0, ncid)
      call wrap_inq_varid(ncid, 'data', keff_id)
      call wrap_get_var_realx(ncid, keff_id, k63_h2o)

      filename = trim(exort_rootdir)//trim(dirk_co2)//trim(k63_co2_file)
      call getfil(filename, locfn, 0)
      call wrap_open(locfn, 0, ncid)
      call wrap_inq_varid(ncid, 'data', keff_id)
      call wrap_get_var_realx(ncid, keff_id, k63_co2)

      filename = trim(exort_rootdir)//trim(dirk_ch4)//trim(k63_ch4_file)
      call getfil(filename, locfn, 0)
      call wrap_open(locfn, 0, ncid)
      call wrap_inq_varid(ncid, 'data', keff_id)
      call wrap_get_var_realx(ncid, keff_id, k63_ch4)


      !----  Load K coefficients, interval 64  ----
      filename = trim(exort_rootdir)//trim(dirk_h2o)//trim(k64_h2o_file)
      call getfil(filename, locfn, 0)
      call wrap_open(locfn, 0, ncid)
      call wrap_inq_varid(ncid, 'data', keff_id)
      call wrap_get_var_realx(ncid, keff_id, k64_h2o)

      filename = trim(exort_rootdir)//trim(dirk_co2)//trim(k64_co2_file)
      call getfil(filename, locfn, 0)
      call wrap_open(locfn, 0, ncid)
      call wrap_inq_varid(ncid, 'data', keff_id)
      call wrap_get_var_realx(ncid, keff_id, k64_co2)

      filename = trim(exort_rootdir)//trim(dirk_ch4)//trim(k64_ch4_file)
      call getfil(filename, locfn, 0)
      call wrap_open(locfn, 0, ncid)
      call wrap_inq_varid(ncid, 'data', keff_id)
      call wrap_get_var_realx(ncid, keff_id, k64_ch4)


      !----  Load K coefficients, interval 65  ----
      filename = trim(exort_rootdir)//trim(dirk_h2o)//trim(k65_h2o_file)
      call getfil(filename, locfn, 0)
      call wrap_open(locfn, 0, ncid)
      call wrap_inq_varid(ncid, 'data', keff_id)
      call wrap_get_var_realx(ncid, keff_id, k65_h2o)

      filename = trim(exort_rootdir)//trim(dirk_co2)//trim(k65_co2_file)
      call getfil(filename, locfn, 0)
      call wrap_open(locfn, 0, ncid)
      call wrap_inq_varid(ncid, 'data', keff_id)
      call wrap_get_var_realx(ncid, keff_id, k65_co2)

      filename = trim(exort_rootdir)//trim(dirk_ch4)//trim(k65_ch4_file)
      call getfil(filename, locfn, 0)
      call wrap_open(locfn, 0, ncid)
      call wrap_inq_varid(ncid, 'data', keff_id)
      call wrap_get_var_realx(ncid, keff_id, k65_ch4)


      !----  Load K coefficients, interval 66  ----
      filename = trim(exort_rootdir)//trim(dirk_h2o)//trim(k66_h2o_file)
      call getfil(filename, locfn, 0)
      call wrap_open(locfn, 0, ncid)
      call wrap_inq_varid(ncid, 'data', keff_id)
      call wrap_get_var_realx(ncid, keff_id, k66_h2o)

      filename = trim(exort_rootdir)//trim(dirk_co2)//trim(k66_co2_file)
      call getfil(filename, locfn, 0)
      call wrap_open(locfn, 0, ncid)
      call wrap_inq_varid(ncid, 'data', keff_id)
      call wrap_get_var_realx(ncid, keff_id, k66_co2)

      filename = trim(exort_rootdir)//trim(dirk_ch4)//trim(k66_ch4_file)
      call getfil(filename, locfn, 0)
      call wrap_open(locfn, 0, ncid)
      call wrap_inq_varid(ncid, 'data', keff_id)
      call wrap_get_var_realx(ncid, keff_id, k66_ch4)


      !----  Load K coefficients, interval 67  ----
      filename = trim(exort_rootdir)//trim(dirk_h2o)//trim(k67_h2o_file)
      call getfil(filename, locfn, 0)
      call wrap_open(locfn, 0, ncid)
      call wrap_inq_varid(ncid, 'data', keff_id)
      call wrap_get_var_realx(ncid, keff_id, k67_h2o)

      filename = trim(exort_rootdir)//trim(dirk_co2)//trim(k67_co2_file)
      call getfil(filename, locfn, 0)
      call wrap_open(locfn, 0, ncid)
      call wrap_inq_varid(ncid, 'data', keff_id)
      call wrap_get_var_realx(ncid, keff_id, k67_co2)

      filename = trim(exort_rootdir)//trim(dirk_ch4)//trim(k67_ch4_file)
      call getfil(filename, locfn, 0)
      call wrap_open(locfn, 0, ncid)
      call wrap_inq_varid(ncid, 'data', keff_id)
      call wrap_get_var_realx(ncid, keff_id, k67_ch4)


      !----  Load K coefficients, interval 68  ----
      filename = trim(exort_rootdir)//trim(dirk_h2o)//trim(k68_h2o_file)
      call getfil(filename, locfn, 0)
      call wrap_open(locfn, 0, ncid)
      call wrap_inq_varid(ncid, 'data', keff_id)
      call wrap_get_var_realx(ncid, keff_id, k68_h2o)

      filename = trim(exort_rootdir)//trim(dirk_co2)//trim(k68_co2_file)
      call getfil(filename, locfn, 0)
      call wrap_open(locfn, 0, ncid)
      call wrap_inq_varid(ncid, 'data', keff_id)
      call wrap_get_var_realx(ncid, keff_id, k68_co2)

      filename = trim(exort_rootdir)//trim(dirk_ch4)//trim(k68_ch4_file)
      call getfil(filename, locfn, 0)
      call wrap_open(locfn, 0, ncid)
      call wrap_inq_varid(ncid, 'data', keff_id)
      call wrap_get_var_realx(ncid, keff_id, k68_ch4)


      ! Load water vapor continuum
      !! mtckd, ngauss point file
      filename = trim(exort_rootdir)//trim(dirct)//trim(kh2o_mtckd_file)
      call getfil(filename, locfn, 0)
      call wrap_open(locfn, 0, ncid)
      call wrap_inq_varid(ncid, 'KSELF', keff_id)
      call wrap_get_var_realx(ncid, keff_id, kh2oself_mtckd)

      filename = trim(exort_rootdir)//trim(dirct)//trim(kh2o_mtckd_file)
      call getfil(filename, locfn, 0)
      call wrap_open(locfn, 0, ncid)
      call wrap_inq_varid(ncid, 'KFRGN', keff_id)
      call wrap_get_var_realx(ncid, keff_id, kh2ofrgn_mtckd)
      !! mtckd

      !! mtckd, average file
      filename = trim(exort_rootdir)//trim(dirct)//trim(kh2o_mtckd_file_avg)
      call getfil(filename, locfn, 0)
      call wrap_open(locfn, 0, ncid)
      call wrap_inq_varid(ncid, 'KSELF', keff_id)
      call wrap_get_var_realx(ncid, keff_id, kh2oself_avg_mtckd)

      filename = trim(exort_rootdir)//trim(dirct)//trim(kh2o_mtckd_file_avg)
      call getfil(filename, locfn, 0)
      call wrap_open(locfn, 0, ncid)
      call wrap_inq_varid(ncid, 'KFRGN', keff_id)
      call wrap_get_var_realx(ncid, keff_id, kh2ofrgn_avg_mtckd)
      !! mtckd


      ! Load absorption coefficients, for h2oh2o continuum
      filename = trim(exort_rootdir)//trim(dirci)//trim(kh2oh2ocia_file )
      call getfil(filename, locfn, 0)
      call wrap_open(locfn, 0, ncid)
      call wrap_inq_varid(ncid, 'sigma', keff_id)
      call wrap_get_var_realx(ncid, keff_id, kh2oh2o )

      ! Load absorption coefficients, for h2on2 continuum
      filename = trim(exort_rootdir)//trim(dirci)//trim(kh2on2cia_file )
      call getfil(filename, locfn, 0)
      call wrap_open(locfn, 0, ncid)
      call wrap_inq_varid(ncid, 'sigma', keff_id)
      call wrap_get_var_realx(ncid, keff_id, kh2on2 )

      ! Load absorption coefficients, for n2n2 continuum
      filename = trim(exort_rootdir)//trim(dirci)//trim(kn2n2cia_file )
      call getfil(filename, locfn, 0)
      call wrap_open(locfn, 0, ncid)
      call wrap_inq_varid(ncid, 'sigma', keff_id)
      call wrap_get_var_realx(ncid, keff_id, kn2n2 )

      ! Load absorption coefficients, for n2h2 continuum
      filename = trim(exort_rootdir)//trim(dirci)//trim(kn2h2cia_file )
      call getfil(filename, locfn, 0)
      call wrap_open(locfn, 0, ncid)
      call wrap_inq_varid(ncid, 'sigma', keff_id)
      call wrap_get_var_realx(ncid, keff_id, kn2h2 )

      ! Load absorption coefficients, for h2h2 continuum
      filename = trim(exort_rootdir)//trim(dirci)//trim(kh2h2cia_file )
      call getfil(filename, locfn, 0)
      call wrap_open(locfn, 0, ncid)
      call wrap_inq_varid(ncid, 'sigma', keff_id)
      call wrap_get_var_realx(ncid, keff_id, kh2h2 )

      ! Load absorption coefficients, for co2co2 lw continuum 
      filename = trim(exort_rootdir)//trim(dirci)//trim(kco2co2cia_lw_file )
      call getfil(filename, locfn, 0)
      call wrap_open(locfn, 0, ncid)
      call wrap_inq_varid(ncid, 'sigma', keff_id)
      call wrap_get_var_realx(ncid, keff_id, kco2co2_lw )

      ! Load absorption coefficients, for co2co2 sw continuum 
      filename = trim(exort_rootdir)//trim(dirci)//trim(kco2co2cia_sw_file )
      call getfil(filename, locfn, 0)
      call wrap_open(locfn, 0, ncid)
      call wrap_inq_varid(ncid, 'sigma', keff_id)
      call wrap_get_var_realx(ncid, keff_id, kco2co2_sw )

      ! Load absorption coefficients, for co2ch4 continuum    
      filename = trim(exort_rootdir)//trim(dirci)//trim(kco2ch4cia_file )
      call getfil(filename, locfn, 0)
      call wrap_open(locfn, 0, ncid)
      call wrap_inq_varid(ncid, 'sigma', keff_id)
      call wrap_get_var_realx(ncid, keff_id, kco2ch4 )

      ! Load absorption coefficients, for co2h2 continuum   
      filename = trim(exort_rootdir)//trim(dirci)//trim(kco2h2cia_file )
      call getfil(filename, locfn, 0)
      call wrap_open(locfn, 0, ncid)
      call wrap_inq_varid(ncid, 'sigma', keff_id)
      call wrap_get_var_realx(ncid, keff_id, kco2h2 )


  end subroutine initialize_kcoeff


!============================================================================

  subroutine initialize_solar

!------------------------------------------------------------------------
!
! Purpose:  Initialize solar data from input file.
!
!------------------------------------------------------------------------
!
!#if ( defined SPMD)
!  use mpishorthand
!#endif

    use ioFileMod, only: getfil

    implicit none
    include 'netcdf.inc'

!------------------------------------------------------------------------
!
! Local Variables
!
    integer :: ncid
    character(len=256) :: locfn, filename
!    integer :: sunm_id
    integer :: solarflux_id
    integer :: S0_id
!    integer :: wav_low_id
!    integer :: wav_high_id
!    integer :: dwm_id


    !if ( masterproc ) then

      ! Load solar data
      filename = trim(exort_rootdir)//trim(dirsol)//trim(solar_file)
      call getfil(filename, locfn, 0)
      call wrap_open(locfn, 0, ncid)
!      call wrap_inq_varid(ncid, 'wm1', wm1_id)
!      call wrap_get_var_realx(ncid, wm1_id, wm1_in )
!      call wrap_inq_varid(ncid, 'wm2', wm2_id)
!      call wrap_get_var_realx(ncid, wm2_id, wm2 )
!      !call wrap_inq_varid(ncid, 'dwm', dwm_id)
!      !call wrap_get_var_realx(ncid, dwm_id, dwm )
      call wrap_inq_varid(ncid, 'S0', S0_id)
      call wrap_get_var_realx(ncid, S0_id, S0 )      
      call wrap_inq_varid(ncid, 'solarflux', solarflux_id)
      call wrap_get_var_realx(ncid, solarflux_id, solarflux )

    !endif

!#if ( defined SPMD )
!      call mpibcast(wm1, 1, mpir8, 0, mpicom)
!      call mpibcast(wm2, 1, mpir8, 0, mpicom)
!      call mpibcast(dwm, 1, mpir8, 0, mpicom)
!      call mpibcast(S0, 1, mpir8, 0, mpicom)
!      call mpibcast(sunm, nw, mpir8, 0, mpicom)
!#endif

  end subroutine initialize_solar


!============================================================================

  subroutine initialize_cldopts

!------------------------------------------------------------------------
!
! Purpose:  Initialize the cloud optical constants from input file.
!
!------------------------------------------------------------------------

!#if ( defined SPMD)
!  use mpishorthand
!#endif

    use ioFileMod, only: getfil

!    implicit none
!    include 'netcdf.inc'

!------------------------------------------------------------------------
!
! Local Variables
!    
    integer :: ncid
    integer :: bin_id
    integer :: wav_id
    integer :: ncldopt_lbins
    integer :: ncldopt_lwavs
    integer :: ncldopt_ibins
    integer :: ncldopt_iwavs
    integer :: q_id
    integer :: w_id
    integer :: g_id
    character(len=256) :: locfn, filename

!------------------------------------------------------------------------
!
! Start Code
!
!    !if ( masterproc ) then
      
      write(6,*) "CLDOPTS: INITIALIZING WATER CLOUD OPTICAL PROPERTIES"

      filename = trim(exort_rootdir)//trim(dircld)//trim(cldoptsL_file)
      ! Load K water cloud optics file
      call getfil(filename, locfn, 0)
      call wrap_open(locfn, 0, ncid)
      call wrap_inq_dimid(ncid, 'rel_bins', bin_id)
      call wrap_inq_dimid(ncid, 'nwavlrng', wav_id)

      call wrap_inq_dimlen(ncid, bin_id, ncldopt_lbins) 
      call wrap_inq_dimlen(ncid, wav_id, ncldopt_lwavs) 

      write(6,*) "CLDOPTS: nrel = ",ncldopt_lbins
      write(6,*) "CLDOPTS: nwavlrng = ",ncldopt_lwavs

      if (ncldopt_lwavs .ne. ntot_wavlnrng .or. ncldopt_lbins .ne. nrel) then
        write(6,*) "CLDOPTS: file size mismatch, liquid" 
!        call endrun
      end if

      call wrap_inq_varid(ncid, 'Qext_liq', q_id)
      call wrap_inq_varid(ncid, 'W_liq', w_id)
      call wrap_inq_varid(ncid, 'G_liq', g_id)

      call wrap_get_var_realx(ncid, q_id, Qcldliq)
      call wrap_get_var_realx(ncid, w_id, Wcldliq)
      call wrap_get_var_realx(ncid, g_id, Gcldliq)

!      write(*,*) "Qcldliq", Qcldliq(1,1), Qcldliq(2,1), Qcldliq(3,1)
!      write(*,*) "should be, "
!      write(*,*) "Wcldliq", Wcldliq(1,1), Wcldliq(2,1), Wcldliq(3,1)
!      write(*,*) "should be, "
!      write(*,*) "Gcldliq", Gcldliq(1,1), Gcldliq(2,1), Gcldliq(3,1)
!      write(*,*) "should be, "

      write(6,*) "CLDOPTS: INITIALIZING ICE OPTICAL PROPERTIES"

!      ! Load ice cloud optics file
     filename = trim(exort_rootdir)//trim(dircld)//trim(cldoptsI_file) 
      call getfil(filename, locfn, 0)

      call wrap_open(locfn, 0, ncid)
      call wrap_inq_dimid(ncid, 'rei_bins', bin_id)
      call wrap_inq_dimid(ncid, 'nwavlrng', wav_id)

      call wrap_inq_dimlen(ncid, bin_id, ncldopt_ibins) 
      call wrap_inq_dimlen(ncid, wav_id, ncldopt_iwavs) 

      write(6,*) "CLDOPTS: nrei = ",ncldopt_ibins
      write(6,*) "CLDOPTS: nwavlrng = ",ncldopt_iwavs

      if (ncldopt_iwavs .ne. ntot_wavlnrng .or. ncldopt_ibins .ne. nrei) then
        write(6,*) "CLDOPTS: file size mismatch, ice" 
 !       call endrun
      end if

      call wrap_inq_varid(ncid, 'Qext_ice', q_id)
      call wrap_inq_varid(ncid, 'W_ice', w_id)
      call wrap_inq_varid(ncid, 'G_ice', g_id)

      call wrap_get_var_realx(ncid, q_id, Qcldice)
      call wrap_get_var_realx(ncid, w_id, Wcldice)
      call wrap_get_var_realx(ncid, g_id, Gcldice)

      write(6,*) "CLDOPTS: INITIALIZING CO2 ICE OPTICAL PROPERTIES"

!      ! Load ice cloud optics file
     filename = trim(exort_rootdir)//trim(dircld)//trim(cldoptsICO2_file) 
      call getfil(filename, locfn, 0)

      call wrap_open(locfn, 0, ncid)
      call wrap_inq_dimid(ncid, 'rei_bins', bin_id)
      call wrap_inq_dimid(ncid, 'nwavlrng', wav_id)

      call wrap_inq_dimlen(ncid, bin_id, ncldopt_ibins) 
      call wrap_inq_dimlen(ncid, wav_id, ncldopt_iwavs) 

      write(6,*) "CLDOPTS: nrei = ",ncldopt_ibins
      write(6,*) "CLDOPTS: nwavlrng = ",ncldopt_iwavs

      if (ncldopt_iwavs .ne. ntot_wavlnrng .or. ncldopt_ibins .ne. nrei) then
        write(6,*) "CLDOPTS: file size mismatch, ice" 
  !      call endrun
      end if

      call wrap_inq_varid(ncid, 'Qext', q_id)
      call wrap_inq_varid(ncid, 'W', w_id)
      call wrap_inq_varid(ncid, 'G', g_id)

      call wrap_get_var_realx(ncid, q_id, Qcldice)
      call wrap_get_var_realx(ncid, w_id, Wcldice)
      call wrap_get_var_realx(ncid, g_id, Gcldice)


!      write(*,*) "Qcldice", Qcldice(1,1), Qcldice(2,1), Qcldice(3,1)
!      write(*,*) "should be, "
!      write(*,*) "Wcldice", Wcldice(1,1), Wcldice(2,1), Wcldice(3,1)
!      write(*,*) "should be, "
!      write(*,*) "Gcldice", Gcldice(1,1), Gcldice(2,1), Gcldice(3,1)
!      write(*,*) "should be, "

   !end if ! masterproc

! broadcast water cloud optical constants to all nodes
!#if ( defined SPMD )
!      call mpibcast(Qcldliq, nrel*ntot_wavlnrng, mpir8, 0, mpicom)
!      call mpibcast(Wcldliq, nrel*ntot_wavlnrng, mpir8, 0, mpicom)
!      call mpibcast(Gcldliq, nrel*ntot_wavlnrng, mpir8, 0, mpicom)
!#endif

! broadcast ice cloud optical constants to all nodes
!#if ( defined SPMD )
!      call mpibcast(Qcldice, nrei*ntot_wavlnrng, mpir8, 0, mpicom)
!      call mpibcast(Wcldice, nrei*ntot_wavlnrng, mpir8, 0, mpicom)
!      call mpibcast(Gcldice, nrei*ntot_wavlnrng, mpir8, 0, mpicom)
!#endif
!
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



end module initialize_rad_mod_1D
