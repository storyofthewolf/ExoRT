#include <params.h>
#include <misc.h>

module initialize_rad_mod

! version h2on68
!
! read in and initialize radiaive transfer grids
!

use kabs
use solar
use radgrid

implicit none
private 
save

!
! Pubic Interfaces
!
  public :: initialize_kcoeff
  public :: initialize_solar
  !public :: initialize_cldopts
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
    character(len=256) :: locfn

!------------------------------------------------------------------------
!
! Start Code
!
!    if ( masterproc ) then

      write (6, '(2x, a)') '_______________________________________________________'
      write (6, '(2x, a)') '_________ initializing gas absorption coeffs __________'
      write (6, '(2x, a)') '_______________________________________________________'
     
      ! Load K coefficients, interval 1  
      call getfil(k01_file, locfn, 0)
      call wrap_open(locfn, 0, ncid)
      call wrap_inq_varid(ncid, 'data', keff_id)
      call wrap_get_var_realx(ncid, keff_id, k01_data)

      ! Load K coefficients, interval 2  
      call getfil(k02_file, locfn, 0)
      call wrap_open(locfn, 0, ncid)
      call wrap_inq_varid(ncid, 'data', keff_id)
      call wrap_get_var_realx(ncid, keff_id, k02_data)

      ! Load K coefficients, interval 3  
      call getfil(k03_file, locfn, 0)
      call wrap_open(locfn, 0, ncid)
      call wrap_inq_varid(ncid, 'data', keff_id)
      call wrap_get_var_realx(ncid, keff_id, k03_data)

      ! Load K coefficients, interval 4  
      call getfil(k04_file, locfn, 0)
      call wrap_open(locfn, 0, ncid)
      call wrap_inq_varid(ncid, 'data', keff_id)
      call wrap_get_var_realx(ncid, keff_id, k04_data)

      ! Load K coefficients, interval 5  
      call getfil(k05_file, locfn, 0)
      call wrap_open(locfn, 0, ncid)
      call wrap_inq_varid(ncid, 'data', keff_id)
      call wrap_get_var_realx(ncid, keff_id, k05_data)

      ! Load K coefficients, interval 6  
      call getfil(k06_file, locfn, 0)
      call wrap_open(locfn, 0, ncid)
      call wrap_inq_varid(ncid, 'data', keff_id)
      call wrap_get_var_realx(ncid, keff_id, k06_data)

      ! Load K coefficients, interval 7  
      call getfil(k07_file, locfn, 0)
      call wrap_open(locfn, 0, ncid)
      call wrap_inq_varid(ncid, 'data', keff_id)
      call wrap_get_var_realx(ncid, keff_id, k07_data)

      ! Load K coefficients, interval 8  
      call getfil(k08_file, locfn, 0)
      call wrap_open(locfn, 0, ncid)
      call wrap_inq_varid(ncid, 'data', keff_id)
      call wrap_get_var_realx(ncid, keff_id, k08_data)

      ! Load K coefficients, interval 9  
      call getfil(k09_file, locfn, 0)
      call wrap_open(locfn, 0, ncid)
      call wrap_inq_varid(ncid, 'data', keff_id)
      call wrap_get_var_realx(ncid, keff_id, k09_data)

      ! Load K coefficients, interval 10  
      call getfil(k10_file, locfn, 0)
      call wrap_open(locfn, 0, ncid)
      call wrap_inq_varid(ncid, 'data', keff_id)
      call wrap_get_var_realx(ncid, keff_id, k10_data)

      ! Load K coefficients, interval 11  
      call getfil(k11_file, locfn, 0)
      call wrap_open(locfn, 0, ncid)
      call wrap_inq_varid(ncid, 'data', keff_id)
      call wrap_get_var_realx(ncid, keff_id, k11_data)

      ! Load K coefficients, interval 12  
      call getfil(k12_file, locfn, 0)
      call wrap_open(locfn, 0, ncid)
      call wrap_inq_varid(ncid, 'data', keff_id)
      call wrap_get_var_realx(ncid, keff_id, k12_data)

      ! Load K coefficients, interval 13  
      call getfil(k13_file, locfn, 0)
      call wrap_open(locfn, 0, ncid)
      call wrap_inq_varid(ncid, 'data', keff_id)
      call wrap_get_var_realx(ncid, keff_id, k13_data)

      ! Load K coefficients, interval 14  
      call getfil(k14_file, locfn, 0)
      call wrap_open(locfn, 0, ncid)
      call wrap_inq_varid(ncid, 'data', keff_id)
      call wrap_get_var_realx(ncid, keff_id, k14_data)

      ! Load K coefficients, interval 15  
      call getfil(k15_file, locfn, 0)
      call wrap_open(locfn, 0, ncid)
      call wrap_inq_varid(ncid, 'data', keff_id)
      call wrap_get_var_realx(ncid, keff_id, k15_data)

      ! Load K coefficients, interval 16  
      call getfil(k16_file, locfn, 0)
      call wrap_open(locfn, 0, ncid)
      call wrap_inq_varid(ncid, 'data', keff_id)
      call wrap_get_var_realx(ncid, keff_id, k16_data)

      ! Load K coefficients, interval 17  
      call getfil(k17_file, locfn, 0)
      call wrap_open(locfn, 0, ncid)
      call wrap_inq_varid(ncid, 'data', keff_id)
      call wrap_get_var_realx(ncid, keff_id, k17_data)

      ! Load K coefficients, interval 18  
      call getfil(k18_file, locfn, 0)
      call wrap_open(locfn, 0, ncid)
      call wrap_inq_varid(ncid, 'data', keff_id)
      call wrap_get_var_realx(ncid, keff_id, k18_data)

      ! Load K coefficients, interval 19  
      call getfil(k19_file, locfn, 0)
      call wrap_open(locfn, 0, ncid)
      call wrap_inq_varid(ncid, 'data', keff_id)
      call wrap_get_var_realx(ncid, keff_id, k19_data)

      ! Load K coefficients, interval 20  
      call getfil(k20_file, locfn, 0)
      call wrap_open(locfn, 0, ncid)
      call wrap_inq_varid(ncid, 'data', keff_id)
      call wrap_get_var_realx(ncid, keff_id, k20_data)

      ! Load K coefficients, interval 21  
      call getfil(k21_file, locfn, 0)
      call wrap_open(locfn, 0, ncid)
      call wrap_inq_varid(ncid, 'data', keff_id)
      call wrap_get_var_realx(ncid, keff_id, k21_data)

      ! Load K coefficients, interval 22  
      call getfil(k22_file, locfn, 0)
      call wrap_open(locfn, 0, ncid)
      call wrap_inq_varid(ncid, 'data', keff_id)
      call wrap_get_var_realx(ncid, keff_id, k22_data)

      ! Load K coefficients, interval 23  
      call getfil(k23_file, locfn, 0)
      call wrap_open(locfn, 0, ncid)
      call wrap_inq_varid(ncid, 'data', keff_id)
      call wrap_get_var_realx(ncid, keff_id, k23_data)

      ! Load K coefficients, interval 24  
      call getfil(k24_file, locfn, 0)
      call wrap_open(locfn, 0, ncid)
      call wrap_inq_varid(ncid, 'data', keff_id)
      call wrap_get_var_realx(ncid, keff_id, k24_data)

      ! Load K coefficients, interval 25  
      call getfil(k25_file, locfn, 0)
      call wrap_open(locfn, 0, ncid)
      call wrap_inq_varid(ncid, 'data', keff_id)
      call wrap_get_var_realx(ncid, keff_id, k25_data)

      ! Load K coefficients, interval 26  
      call getfil(k26_file, locfn, 0)
      call wrap_open(locfn, 0, ncid)
      call wrap_inq_varid(ncid, 'data', keff_id)
      call wrap_get_var_realx(ncid, keff_id, k26_data)

      ! Load K coefficients, interval 27  
      call getfil(k27_file, locfn, 0)
      call wrap_open(locfn, 0, ncid)
      call wrap_inq_varid(ncid, 'data', keff_id)
      call wrap_get_var_realx(ncid, keff_id, k27_data)

      ! Load K coefficients, interval 28  
      call getfil(k28_file, locfn, 0)
      call wrap_open(locfn, 0, ncid)
      call wrap_inq_varid(ncid, 'data', keff_id)
      call wrap_get_var_realx(ncid, keff_id, k28_data)

      ! Load K coefficients, interval 29  
      call getfil(k29_file, locfn, 0)
      call wrap_open(locfn, 0, ncid)
      call wrap_inq_varid(ncid, 'data', keff_id)
      call wrap_get_var_realx(ncid, keff_id, k29_data)

      ! Load K coefficients, interval 30  
      call getfil(k30_file, locfn, 0)
      call wrap_open(locfn, 0, ncid)
      call wrap_inq_varid(ncid, 'data', keff_id)
      call wrap_get_var_realx(ncid, keff_id, k30_data)


      ! Load K coefficients, interval 31  
      call getfil(k31_file, locfn, 0)
      call wrap_open(locfn, 0, ncid)
      call wrap_inq_varid(ncid, 'data', keff_id)
      call wrap_get_var_realx(ncid, keff_id, k31_data)

      ! Load K coefficients, interval 32  
      call getfil(k32_file, locfn, 0)
      call wrap_open(locfn, 0, ncid)
      call wrap_inq_varid(ncid, 'data', keff_id)
      call wrap_get_var_realx(ncid, keff_id, k32_data)

      ! Load K coefficients, interval 33  
      call getfil(k33_file, locfn, 0)
      call wrap_open(locfn, 0, ncid)
      call wrap_inq_varid(ncid, 'data', keff_id)
      call wrap_get_var_realx(ncid, keff_id, k33_data)

      ! Load K coefficients, interval 34  
      call getfil(k34_file, locfn, 0)
      call wrap_open(locfn, 0, ncid)
      call wrap_inq_varid(ncid, 'data', keff_id)
      call wrap_get_var_realx(ncid, keff_id, k34_data)

      ! Load K coefficients, interval 35  
      call getfil(k35_file, locfn, 0)
      call wrap_open(locfn, 0, ncid)
      call wrap_inq_varid(ncid, 'data', keff_id)
      call wrap_get_var_realx(ncid, keff_id, k35_data)

      ! Load K coefficients, interval 36  
      call getfil(k36_file, locfn, 0)
      call wrap_open(locfn, 0, ncid)
      call wrap_inq_varid(ncid, 'data', keff_id)
      call wrap_get_var_realx(ncid, keff_id, k36_data)

      ! Load K coefficients, interval 37  
      call getfil(k37_file, locfn, 0)
      call wrap_open(locfn, 0, ncid)
      call wrap_inq_varid(ncid, 'data', keff_id)
      call wrap_get_var_realx(ncid, keff_id, k37_data)

      ! Load K coefficients, interval 38  
      call getfil(k38_file, locfn, 0)
      call wrap_open(locfn, 0, ncid)
      call wrap_inq_varid(ncid, 'data', keff_id)
      call wrap_get_var_realx(ncid, keff_id, k38_data)

      ! Load K coefficients, interval 39  
      call getfil(k39_file, locfn, 0)
      call wrap_open(locfn, 0, ncid)
      call wrap_inq_varid(ncid, 'data', keff_id)
      call wrap_get_var_realx(ncid, keff_id, k39_data)

      ! Load K coefficients, interval 40  
      call getfil(k40_file, locfn, 0)
      call wrap_open(locfn, 0, ncid)
      call wrap_inq_varid(ncid, 'data', keff_id)
      call wrap_get_var_realx(ncid, keff_id, k40_data)

      ! Load K coefficients, interval 41  
      call getfil(k41_file, locfn, 0)
      call wrap_open(locfn, 0, ncid)
      call wrap_inq_varid(ncid, 'data', keff_id)
      call wrap_get_var_realx(ncid, keff_id, k41_data)

      ! Load K coefficients, interval 42  
      call getfil(k42_file, locfn, 0)
      call wrap_open(locfn, 0, ncid)
      call wrap_inq_varid(ncid, 'data', keff_id)
      call wrap_get_var_realx(ncid, keff_id, k42_data)

      ! Load K coefficients, interval 43  
      call getfil(k43_file, locfn, 0)
      call wrap_open(locfn, 0, ncid)
      call wrap_inq_varid(ncid, 'data', keff_id)
      call wrap_get_var_realx(ncid, keff_id, k43_data)

      ! Load K coefficients, interval 44  
      call getfil(k44_file, locfn, 0)
      call wrap_open(locfn, 0, ncid)
      call wrap_inq_varid(ncid, 'data', keff_id)
      call wrap_get_var_realx(ncid, keff_id, k44_data)

      ! Load K coefficients, interval 45  
      call getfil(k45_file, locfn, 0)
      call wrap_open(locfn, 0, ncid)
      call wrap_inq_varid(ncid, 'data', keff_id)
      call wrap_get_var_realx(ncid, keff_id, k45_data)

      ! Load K coefficients, interval 46  
      call getfil(k46_file, locfn, 0)
      call wrap_open(locfn, 0, ncid)
      call wrap_inq_varid(ncid, 'data', keff_id)
      call wrap_get_var_realx(ncid, keff_id, k46_data)

      ! Load K coefficients, interval 47  
      call getfil(k47_file, locfn, 0)
      call wrap_open(locfn, 0, ncid)
      call wrap_inq_varid(ncid, 'data', keff_id)
      call wrap_get_var_realx(ncid, keff_id, k47_data)

      ! Load K coefficients, interval 48  
      call getfil(k48_file, locfn, 0)
      call wrap_open(locfn, 0, ncid)
      call wrap_inq_varid(ncid, 'data', keff_id)
      call wrap_get_var_realx(ncid, keff_id, k48_data)

      ! Load K coefficients, interval 49  
      call getfil(k49_file, locfn, 0)
      call wrap_open(locfn, 0, ncid)
      call wrap_inq_varid(ncid, 'data', keff_id)
      call wrap_get_var_realx(ncid, keff_id, k49_data)

      ! Load K coefficients, interval 50  
      call getfil(k50_file, locfn, 0)
      call wrap_open(locfn, 0, ncid)
      call wrap_inq_varid(ncid, 'data', keff_id)
      call wrap_get_var_realx(ncid, keff_id, k50_data)

      ! Load K coefficients, interval 51  
      call getfil(k51_file, locfn, 0)
      call wrap_open(locfn, 0, ncid)
      call wrap_inq_varid(ncid, 'data', keff_id)
      call wrap_get_var_realx(ncid, keff_id, k51_data)

      ! Load K coefficients, interval 52  
      call getfil(k52_file, locfn, 0)
      call wrap_open(locfn, 0, ncid)
      call wrap_inq_varid(ncid, 'data', keff_id)
      call wrap_get_var_realx(ncid, keff_id, k52_data)

      ! Load K coefficients, interval 53  
      call getfil(k53_file, locfn, 0)
      call wrap_open(locfn, 0, ncid)
      call wrap_inq_varid(ncid, 'data', keff_id)
      call wrap_get_var_realx(ncid, keff_id, k53_data)

      ! Load K coefficients, interval 54  
      call getfil(k54_file, locfn, 0)
      call wrap_open(locfn, 0, ncid)
      call wrap_inq_varid(ncid, 'data', keff_id)
      call wrap_get_var_realx(ncid, keff_id, k54_data)

      ! Load K coefficients, interval 55  
      call getfil(k55_file, locfn, 0)
      call wrap_open(locfn, 0, ncid)
      call wrap_inq_varid(ncid, 'data', keff_id)
      call wrap_get_var_realx(ncid, keff_id, k55_data)

      ! Load K coefficients, interval 56  
      call getfil(k56_file, locfn, 0)
      call wrap_open(locfn, 0, ncid)
      call wrap_inq_varid(ncid, 'data', keff_id)
      call wrap_get_var_realx(ncid, keff_id, k56_data)

      ! Load K coefficients, interval 57  
      call getfil(k57_file, locfn, 0)
      call wrap_open(locfn, 0, ncid)
      call wrap_inq_varid(ncid, 'data', keff_id)
      call wrap_get_var_realx(ncid, keff_id, k57_data)

      ! Load K coefficients, interval 58  
      call getfil(k58_file, locfn, 0)
      call wrap_open(locfn, 0, ncid)
      call wrap_inq_varid(ncid, 'data', keff_id)
      call wrap_get_var_realx(ncid, keff_id, k58_data)

      ! Load K coefficients, interval 59  
      call getfil(k59_file, locfn, 0)
      call wrap_open(locfn, 0, ncid)
      call wrap_inq_varid(ncid, 'data', keff_id)
      call wrap_get_var_realx(ncid, keff_id, k59_data)

      ! Load K coefficients, interval 60  
      call getfil(k60_file, locfn, 0)
      call wrap_open(locfn, 0, ncid)
      call wrap_inq_varid(ncid, 'data', keff_id)
      call wrap_get_var_realx(ncid, keff_id, k60_data)

      ! Load K coefficients, interval 61  
      call getfil(k61_file, locfn, 0)
      call wrap_open(locfn, 0, ncid)
      call wrap_inq_varid(ncid, 'data', keff_id)
      call wrap_get_var_realx(ncid, keff_id, k61_data)

      ! Load K coefficients, interval 62  
      call getfil(k62_file, locfn, 0)
      call wrap_open(locfn, 0, ncid)
      call wrap_inq_varid(ncid, 'data', keff_id)
      call wrap_get_var_realx(ncid, keff_id, k62_data)

      ! Load K coefficients, interval 63  
      call getfil(k63_file, locfn, 0)
      call wrap_open(locfn, 0, ncid)
      call wrap_inq_varid(ncid, 'data', keff_id)
      call wrap_get_var_realx(ncid, keff_id, k63_data)

      ! Load K coefficients, interval 64  
      call getfil(k64_file, locfn, 0)
      call wrap_open(locfn, 0, ncid)
      call wrap_inq_varid(ncid, 'data', keff_id)
      call wrap_get_var_realx(ncid, keff_id, k64_data)

      ! Load K coefficients, interval 65  
      call getfil(k65_file, locfn, 0)
      call wrap_open(locfn, 0, ncid)
      call wrap_inq_varid(ncid, 'data', keff_id)
      call wrap_get_var_realx(ncid, keff_id, k65_data)

      ! Load K coefficients, interval 66  
      call getfil(k66_file, locfn, 0)
      call wrap_open(locfn, 0, ncid)
      call wrap_inq_varid(ncid, 'data', keff_id)
      call wrap_get_var_realx(ncid, keff_id, k66_data)

      ! Load K coefficients, interval 67  
      call getfil(k67_file, locfn, 0)
      call wrap_open(locfn, 0, ncid)
      call wrap_inq_varid(ncid, 'data', keff_id)
      call wrap_get_var_realx(ncid, keff_id, k67_data)

      ! Load K coefficients, interval 68  
      call getfil(k68_file, locfn, 0)
      call wrap_open(locfn, 0, ncid)
      call wrap_inq_varid(ncid, 'data', keff_id)
      call wrap_get_var_realx(ncid, keff_id, k68_data)


      ! Load K coefficients, water vapor self continuum

      !! mtckd
      call getfil(kh2oself_mtckd_file, locfn, 0)
      call wrap_open(locfn, 0, ncid)
      call wrap_inq_varid(ncid, 'KSELF', keff_id)
      call wrap_get_var_realx(ncid, keff_id, kh2oself)
      !! mtckd

      !! bps
      call getfil(kh2oself_bps_file, locfn, 0)
      call wrap_open(locfn, 0, ncid)
      call wrap_inq_varid(ncid, 'self', keff_id)
      call wrap_get_var_realx(ncid, keff_id, self)
      call wrap_inq_varid(ncid, 'foreign', keff_id)
      call wrap_get_var_realx(ncid, keff_id, foreign)
      call wrap_inq_varid(ncid, 'base_self', keff_id)
      call wrap_get_var_realx(ncid, keff_id, base_self)
      call wrap_inq_varid(ncid, 'base_foreign', keff_id)
      call wrap_get_var_realx(ncid, keff_id, base_foreign)
      call wrap_inq_varid(ncid, 'TempCoeff', keff_id)
      call wrap_get_var_realx(ncid, keff_id, TempCoeff)
      !!bps


      ! Load K coefficients, for n2n2 continuum
      call getfil(kn2n2cia_file, locfn, 0)
      call wrap_open(locfn, 0, ncid)
      call wrap_inq_varid(ncid, 'sigma', keff_id)
      call wrap_get_var_realx(ncid, keff_id, kn2n2 ) 
 
      ! Load K coefficients, for h2n2 continuum
      call getfil(kh2n2cia_file, locfn, 0)
      call wrap_open(locfn, 0, ncid)
      call wrap_inq_varid(ncid, 'sigma', keff_id)
      call wrap_get_var_realx(ncid, keff_id, kh2n2 )

      ! Load K coefficients, for h2h2 continuum
      call getfil(kh2h2cia_file, locfn, 0)
      call wrap_open(locfn, 0, ncid)
      call wrap_inq_varid(ncid, 'sigma', keff_id)
      call wrap_get_var_realx(ncid, keff_id, kh2h2 ) 

!   end if ! masterproc

! broadcast optical constants to all nodes
!#if ( defined SPMD )
!      call mpibcast(k01_lower, ngauss_pts(1)*kc_npress_lower*kc_ntemp_lower, mpir8, 0, mpicom)
!      call mpibcast(k01_upper, ngauss_pts(1)*kc_npress_upper*kc_ntemp_upper, mpir8, 0, mpicom)
!      call mpibcast(k02_lower, ngauss_pts(2)*kc_npress_lower*kc_ntemp_lower, mpir8, 0, mpicom)
!      call mpibcast(k02_upper, ngauss_pts(2)*kc_npress_upper*kc_ntemp_upper, mpir8, 0, mpicom)
!      call mpibcast(k03_lower, ngauss_pts(3)*kc_npress_lower*kc_ntemp_lower*kc_nweight, mpir8, 0, mpicom)
!      call mpibcast(k03_upper, ngauss_pts(3)*kc_npress_upper*kc_ntemp_upper*kc_nweight, mpir8, 0, mpicom)
!      call mpibcast(k04_lower, ngauss_pts(4)*kc_npress_lower*kc_ntemp_lower*kc_nweight, mpir8, 0, mpicom)
!      call mpibcast(k04_upper, ngauss_pts(4)*kc_npress_upper*kc_ntemp_upper*kc_nweight, mpir8, 0, mpicom)
!      call mpibcast(k05_lower, ngauss_pts(5)*kc_npress_lower*kc_ntemp_lower*kc_nweight, mpir8, 0, mpicom)
!      call mpibcast(k05_upper, ngauss_pts(5)*kc_npress_upper*kc_ntemp_upper*kc_nweight, mpir8, 0, mpicom)
!      call mpibcast(k06_lower, ngauss_pts(6)*kc_npress_lower*kc_ntemp_lower*kc_nweight, mpir8, 0, mpicom)
!      call mpibcast(k06_upper, ngauss_pts(6)*kc_npress_upper*kc_ntemp_upper*kc_nweight, mpir8, 0, mpicom)
!      call mpibcast(k07_lower, ngauss_pts(7)*kc_npress_lower*kc_ntemp_lower*kc_nweight, mpir8, 0, mpicom)
!      call mpibcast(k07_upper, ngauss_pts(7)*kc_npress_upper*kc_ntemp_upper*kc_nweight, mpir8, 0, mpicom)
!      call mpibcast(k08_lower, ngauss_pts(8)*kc_npress_lower*kc_ntemp_lower*kc_nweight, mpir8, 0, mpicom)
!      call mpibcast(k08_upper, ngauss_pts(8)*kc_npress_upper*kc_ntemp_upper*kc_nweight, mpir8, 0, mpicom)
!      call mpibcast(k09_lower, ngauss_pts(9)*kc_npress_lower*kc_ntemp_lower*kc_nweight, mpir8, 0, mpicom)
!      call mpibcast(k09_upper, ngauss_pts(9)*kc_npress_upper*kc_ntemp_upper*kc_nweight, mpir8, 0, mpicom)
!      call mpibcast(k10_lower, ngauss_pts(10)*kc_npress_lower*kc_ntemp_lower*kc_nweight, mpir8, 0, mpicom)
!      call mpibcast(k10_upper, ngauss_pts(10)*kc_npress_upper*kc_ntemp_upper*kc_nweight, mpir8, 0, mpicom)
!      call mpibcast(k11_lower, ngauss_pts(11)*kc_npress_lower*kc_ntemp_lower*kc_nweight, mpir8, 0, mpicom)
!      call mpibcast(k11_upper, ngauss_pts(11)*kc_npress_upper*kc_ntemp_upper*kc_nweight, mpir8, 0, mpicom)
!      call mpibcast(k12_lower, ngauss_pts(12)*kc_npress_lower*kc_ntemp_lower*kc_nweight, mpir8, 0, mpicom)
!      call mpibcast(k12_upper, ngauss_pts(12)*kc_npress_upper*kc_ntemp_upper*kc_nweight, mpir8, 0, mpicom)
!      call mpibcast(k13_lower, ngauss_pts(13)*kc_npress_lower*kc_ntemp_lower*kc_nweight, mpir8, 0, mpicom)
!      call mpibcast(k13_upper, ngauss_pts(13)*kc_npress_upper*kc_ntemp_upper*kc_nweight, mpir8, 0, mpicom)
!      call mpibcast(k14_lower, ngauss_pts(14)*kc_npress_lower*kc_ntemp_lower, mpir8, 0, mpicom)
!      call mpibcast(k14_upper, ngauss_pts(14)*kc_npress_upper*kc_ntemp_upper, mpir8, 0, mpicom)
!      call mpibcast(k15_lower, ngauss_pts(15)*kc_npress_lower*kc_ntemp_lower, mpir8, 0, mpicom)
!      call mpibcast(k15_upper, ngauss_pts(15)*kc_npress_upper*kc_ntemp_upper, mpir8, 0, mpicom)
!      call mpibcast(k16_lower, ngauss_pts(16)*kc_npress_lower*kc_ntemp_lower*kc_nweight, mpir8, 0, mpicom)
!      call mpibcast(k16_upper, ngauss_pts(16)*kc_npress_upper*kc_ntemp_upper*kc_nweight, mpir8, 0, mpicom)
!      call mpibcast(k17_lower, ngauss_pts(17)*kc_npress_lower*kc_ntemp_lower*kc_nweight, mpir8, 0, mpicom)
!      call mpibcast(k17_upper, ngauss_pts(17)*kc_npress_upper*kc_ntemp_upper*kc_nweight, mpir8, 0, mpicom)
!      call mpibcast(k18_lower, ngauss_pts(18)*kc_npress_lower*kc_ntemp_lower*kc_nweight, mpir8, 0, mpicom)
!      call mpibcast(k18_upper, ngauss_pts(18)*kc_npress_upper*kc_ntemp_upper*kc_nweight, mpir8, 0, mpicom)
!      call mpibcast(k19_lower, ngauss_pts(19)*kc_npress_lower*kc_ntemp_lower*kc_nweight, mpir8, 0, mpicom)
!      call mpibcast(k19_upper, ngauss_pts(19)*kc_npress_upper*kc_ntemp_upper*kc_nweight, mpir8, 0, mpicom)
!      call mpibcast(k20_lower, ngauss_pts(20)*kc_npress_lower*kc_ntemp_lower, mpir8, 0, mpicom)
!      call mpibcast(k20_upper, ngauss_pts(20)*kc_npress_upper*kc_ntemp_upper, mpir8, 0, mpicom)
!      call mpibcast(k21_lower, ngauss_pts(21)*kc_npress_lower*kc_ntemp_lower*kc_nweight, mpir8, 0, mpicom)
!      call mpibcast(k21_upper, ngauss_pts(21)*kc_npress_upper*kc_ntemp_upper*kc_nweight, mpir8, 0, mpicom)
!      call mpibcast(k22_lower, ngauss_pts(22)*kc_npress_lower*kc_ntemp_lower, mpir8, 0, mpicom)
!      call mpibcast(k22_upper, ngauss_pts(22)*kc_npress_upper*kc_ntemp_upper, mpir8, 0, mpicom)
!      call mpibcast(k23_lower, ngauss_pts(23)*kc_npress_lower*kc_ntemp_lower, mpir8, 0, mpicom)
!      call mpibcast(k23_upper, ngauss_pts(23)*kc_npress_upper*kc_ntemp_upper, mpir8, 0, mpicom)
!      call mpibcast(k24_lower, ngauss_pts(24)*kc_npress_lower*kc_ntemp_lower, mpir8, 0, mpicom)
!      call mpibcast(k24_upper, ngauss_pts(24)*kc_npress_upper*kc_ntemp_upper, mpir8, 0, mpicom)
!      call mpibcast(k25_lower, ngauss_pts(25)*kc_npress_lower*kc_ntemp_lower, mpir8, 0, mpicom)
!      call mpibcast(k25_upper, ngauss_pts(25)*kc_npress_upper*kc_ntemp_upper, mpir8, 0, mpicom)
!      call mpibcast(kh2oself_8gpt, ngH2O*ks_npress*ks_ntemp*ks_nweight*ntot_wavlnrng, mpir8, 0, mpicom)
!      call mpibcast(kco2cont_8gpt, ngCO2*ntot_wavlnrng, mpir8, 0, mpicom)
!#endif
!

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
    character(len=256) :: locfn
!    integer :: sunm_id
    integer :: solarflux_id
    integer :: S0_id
!    integer :: wav_low_id
!    integer :: wav_high_id
!    integer :: dwm_id


    !if ( masterproc ) then

      ! Load solar data
      call getfil(solar_file, locfn, 0)
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
!
!  subroutine initialize_cldopts
!
!------------------------------------------------------------------------
!
! Purpose:  Initialize the cloud optical constants from input file.
!
!------------------------------------------------------------------------
!
!#if ( defined SPMD)
!  use mpishorthand
!#endif

!    use ioFileMod, only: getfil

!    implicit none
!    include 'netcdf.inc'

!------------------------------------------------------------------------
!
! Local Variables
!    
!    integer :: ncid
!    integer :: bin_id
!    integer :: wav_id
!    integer :: ncldopt_lbins
!    integer :: ncldopt_lwavs
!    integer :: ncldopt_ibins
!    integer :: ncldopt_iwavs
!    integer :: q_id
!    integer :: w_id
!    integer :: g_id
!    character(len=256) :: locfn

!------------------------------------------------------------------------
!
! Start Code
!
!    !if ( masterproc ) then
      
!      write(6,*) "CLDOPTS: INITIALIZING WATER CLOUD OPTICAL PROPERTIES"

!      ! Load K water cloud optics file
!      call getfil(cldoptsL_file, locfn, 0)

!      call wrap_open(locfn, 0, ncid)

!      call wrap_inq_dimid(ncid, 'rel_bins', bin_id)
!      call wrap_inq_dimid(ncid, 'nwavlrng', wav_id)

!      call wrap_inq_dimlen(ncid, bin_id, ncldopt_lbins) 
!      call wrap_inq_dimlen(ncid, wav_id, ncldopt_lwavs) 

!      write(6,*) "CLDOPTS: nrel = ",ncldopt_lbins
!      write(6,*) "CLDOPTS: nwavlrng = ",ncldopt_lwavs

!      if (ncldopt_lwavs .ne. ntot_wavlnrng .or. ncldopt_lbins .ne. nrel) then
!        write(6,*) "CLDOPTS: file size mismatch, liquid" 
!        call endrun
!      end if

!      call wrap_inq_varid(ncid, 'Qext_liq', q_id)
!      call wrap_inq_varid(ncid, 'W_liq', w_id)
!      call wrap_inq_varid(ncid, 'G_liq', g_id)

!      call wrap_get_var_realx(ncid, q_id, Qcldliq)
!      call wrap_get_var_realx(ncid, w_id, Wcldliq)
!      call wrap_get_var_realx(ncid, g_id, Gcldliq)

!      write(*,*) "Qcldliq", Qcldliq(1,1), Qcldliq(2,1), Qcldliq(3,1)
!      write(*,*) "should be, "
!      write(*,*) "Wcldliq", Wcldliq(1,1), Wcldliq(2,1), Wcldliq(3,1)
!      write(*,*) "should be, "
!      write(*,*) "Gcldliq", Gcldliq(1,1), Gcldliq(2,1), Gcldliq(3,1)
!      write(*,*) "should be, "

!      write(6,*) "CLDOPTS: INITIALIZING ICE OPTICAL PROPERTIES"

!      ! Load ice cloud optics file
!      call getfil(cldoptsI_file, locfn, 0)

!      call wrap_open(locfn, 0, ncid)
!      call wrap_inq_dimid(ncid, 'rei_bins', bin_id)
!      call wrap_inq_dimid(ncid, 'nwavlrng', wav_id)

!      call wrap_inq_dimlen(ncid, bin_id, ncldopt_ibins) 
!      call wrap_inq_dimlen(ncid, wav_id, ncldopt_iwavs) 

!      write(6,*) "CLDOPTS: nrei = ",ncldopt_ibins
!      write(6,*) "CLDOPTS: nwavlrng = ",ncldopt_iwavs

!      if (ncldopt_iwavs .ne. ntot_wavlnrng .or. ncldopt_ibins .ne. nrei) then
!        write(6,*) "CLDOPTS: file size mismatch, ice" 
!        call endrun
!      end if

!      call wrap_inq_varid(ncid, 'Qext_ice', q_id)
!      call wrap_inq_varid(ncid, 'W_ice', w_id)
!      call wrap_inq_varid(ncid, 'G_ice', g_id)

!      call wrap_get_var_realx(ncid, q_id, Qcldice)
!      call wrap_get_var_realx(ncid, w_id, Wcldice)
!      call wrap_get_var_realx(ncid, g_id, Gcldice)

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
!  end subroutine initialize_cldopts


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



end module initialize_rad_mod