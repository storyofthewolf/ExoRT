
module initialize_rad_mod_1D
! version.archean 1D version
!
! read in and initialize radiaive transfer grids
!

use kabs
use exoplanet_mod, only: solar_file
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
    integer :: keff_id_lower
    integer :: keff_id_upper
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
     
      ! Load K coefficients, interval 1  
      filename = trim(dirk)//trim(K01_file)
      call getfil(filename, locfn, 0)
      call wrap_open(locfn, 0, ncid)
      call wrap_inq_varid(ncid, 'KEFF_LOWER', keff_id_lower)
      call wrap_get_var_realx(ncid, keff_id_lower, k01_lower)
      call wrap_inq_varid(ncid, 'KEFF_UPPER', keff_id_upper)
      call wrap_get_var_realx(ncid, keff_id_upper, k01_upper)

      ! Load K coefficients, interval 2  
      filename = trim(dirk)//trim(K02_file)
      call getfil(filename, locfn, 0)
      call wrap_open(locfn, 0, ncid)
      call wrap_inq_varid(ncid, 'KEFF_LOWER', keff_id_lower)
      call wrap_get_var_realx(ncid, keff_id_lower, k02_lower)
      call wrap_inq_varid(ncid, 'KEFF_UPPER', keff_id_upper)
      call wrap_get_var_realx(ncid, keff_id_upper, k02_upper)

      ! Load K coefficients, interval 3
      filename = trim(dirk)//trim(K03_file)
      call getfil(filename, locfn, 0)
      call wrap_open(locfn, 0, ncid)
      call wrap_inq_varid(ncid, 'KEFF_LOWER', keff_id_lower)
      call wrap_get_var_realx(ncid, keff_id_lower, k03_lower)
      call wrap_inq_varid(ncid, 'KEFF_UPPER', keff_id_upper)
      call wrap_get_var_realx(ncid, keff_id_upper, k03_upper)

      ! Load K coefficients, interval 4
      filename = trim(dirk)//trim(K04_file)
      call getfil(filename, locfn, 0)
      call wrap_open(locfn, 0, ncid)
      call wrap_inq_varid(ncid, 'KEFF_LOWER', keff_id_lower)
      call wrap_get_var_realx(ncid, keff_id_lower, k04_lower)
      call wrap_inq_varid(ncid, 'KEFF_UPPER', keff_id_upper)
      call wrap_get_var_realx(ncid, keff_id_upper, k04_upper)

      ! Load K coefficients, interval 5
      filename = trim(dirk)//trim(K05_file)
      call getfil(filename, locfn, 0)
      call wrap_open(locfn, 0, ncid)
      call wrap_inq_varid(ncid, 'KEFF_LOWER', keff_id_lower)
      call wrap_get_var_realx(ncid, keff_id_lower, k05_lower)
      call wrap_inq_varid(ncid, 'KEFF_UPPER', keff_id_upper)
      call wrap_get_var_realx(ncid, keff_id_upper, k05_upper)

      ! Load K coefficients, interval 6
      filename = trim(dirk)//trim(K06_file)
      call getfil(filename, locfn, 0)
      call wrap_open(locfn, 0, ncid)
      call wrap_inq_varid(ncid, 'KEFF_LOWER', keff_id_lower)
      call wrap_get_var_realx(ncid, keff_id_lower, k06_lower)
      call wrap_inq_varid(ncid, 'KEFF_UPPER', keff_id_upper)
      call wrap_get_var_realx(ncid, keff_id_upper, k06_upper)

      ! Load K coefficients, interval 7
      filename = trim(dirk)//trim(K07_file)
      call getfil(filename, locfn, 0)
      call wrap_open(locfn, 0, ncid)
      call wrap_inq_varid(ncid, 'KEFF_LOWER', keff_id_lower)
      call wrap_get_var_realx(ncid, keff_id_lower, k07_lower)
      call wrap_inq_varid(ncid, 'KEFF_UPPER', keff_id_upper)
      call wrap_get_var_realx(ncid, keff_id_upper, k07_upper)

      ! Load K coefficients, interval 8
      filename = trim(dirk)//trim(K08_file)
      call getfil(filename, locfn, 0)
      call wrap_open(locfn, 0, ncid)
      call wrap_inq_varid(ncid, 'KEFF_LOWER', keff_id_lower)
      call wrap_get_var_realx(ncid, keff_id_lower, k08_lower)
      call wrap_inq_varid(ncid, 'KEFF_UPPER', keff_id_upper)
      call wrap_get_var_realx(ncid, keff_id_upper, k08_upper)

      ! Load K coefficients, interval 9
      filename = trim(dirk)//trim(K09_file)
      call getfil(filename, locfn, 0)
      call wrap_open(locfn, 0, ncid)
      call wrap_inq_varid(ncid, 'KEFF_LOWER', keff_id_lower)
      call wrap_get_var_realx(ncid, keff_id_lower, k09_lower)
      call wrap_inq_varid(ncid, 'KEFF_UPPER', keff_id_upper)
      call wrap_get_var_realx(ncid, keff_id_upper, k09_upper)

      ! Load K coefficients, interval 10
      filename = trim(dirk)//trim(K10_file)
      call getfil(filename, locfn, 0)
      call wrap_open(locfn, 0, ncid)
      call wrap_inq_varid(ncid, 'KEFF_LOWER', keff_id_lower)
      call wrap_get_var_realx(ncid, keff_id_lower, k10_lower)
      call wrap_inq_varid(ncid, 'KEFF_UPPER', keff_id_upper)
      call wrap_get_var_realx(ncid, keff_id_upper, k10_upper)

      ! Load K coefficients, interval 11
      filename = trim(dirk)//trim(K11_file)
      call getfil(filename, locfn, 0)
      call wrap_open(locfn, 0, ncid)
      call wrap_inq_varid(ncid, 'KEFF_LOWER', keff_id_lower)
      call wrap_get_var_realx(ncid, keff_id_lower, k11_lower)
      call wrap_inq_varid(ncid, 'KEFF_UPPER', keff_id_upper)
      call wrap_get_var_realx(ncid, keff_id_upper, k11_upper)

      ! Load K coefficients, interval 12
      filename = trim(dirk)//trim(K12_file)
      call getfil(filename, locfn, 0)
      call wrap_open(locfn, 0, ncid)
      call wrap_inq_varid(ncid, 'KEFF_LOWER', keff_id_lower)
      call wrap_get_var_realx(ncid, keff_id_lower, k12_lower)
      call wrap_inq_varid(ncid, 'KEFF_UPPER', keff_id_upper)
      call wrap_get_var_realx(ncid, keff_id_upper, k12_upper)

      ! Load K coefficients, interval 13
      filename = trim(dirk)//trim(K13_file)
      call getfil(filename, locfn, 0)
      call wrap_open(locfn, 0, ncid)
      call wrap_inq_varid(ncid, 'KEFF_LOWER', keff_id_lower)
      call wrap_get_var_realx(ncid, keff_id_lower, k13_lower)
      call wrap_inq_varid(ncid, 'KEFF_UPPER', keff_id_upper)
      call wrap_get_var_realx(ncid, keff_id_upper, k13_upper)

      ! Load K coefficients, interval 14
      filename = trim(dirk)//trim(K14_file)
      call getfil(filename, locfn, 0)
      call wrap_open(locfn, 0, ncid)
      call wrap_inq_varid(ncid, 'KEFF_LOWER', keff_id_lower)
      call wrap_get_var_realx(ncid, keff_id_lower, k14_lower)
      call wrap_inq_varid(ncid, 'KEFF_UPPER', keff_id_upper)
      call wrap_get_var_realx(ncid, keff_id_upper, k14_upper)

      ! Load K coefficients, interval 15
      filename = trim(dirk)//trim(K15_file)
      call getfil(filename, locfn, 0)
      call wrap_open(locfn, 0, ncid)
      call wrap_inq_varid(ncid, 'KEFF_LOWER', keff_id_lower)
      call wrap_get_var_realx(ncid, keff_id_lower, k15_lower)
      call wrap_inq_varid(ncid, 'KEFF_UPPER', keff_id_upper)
      call wrap_get_var_realx(ncid, keff_id_upper, k15_upper)

      ! Load K coefficients, interval 16
      filename = trim(dirk)//trim(K16_file)
      call getfil(filename, locfn, 0)
      call wrap_open(locfn, 0, ncid)
      call wrap_inq_varid(ncid, 'KEFF_LOWER', keff_id_lower)
      call wrap_get_var_realx(ncid, keff_id_lower, k16_lower)
      call wrap_inq_varid(ncid, 'KEFF_UPPER', keff_id_upper)
      call wrap_get_var_realx(ncid, keff_id_upper, k16_upper)

      ! Load K coefficients, interval 17
      filename = trim(dirk)//trim(K17_file)
      call getfil(filename, locfn, 0)
      call wrap_open(locfn, 0, ncid)
      call wrap_inq_varid(ncid, 'KEFF_LOWER', keff_id_lower)
      call wrap_get_var_realx(ncid, keff_id_lower, k17_lower)
      call wrap_inq_varid(ncid, 'KEFF_UPPER', keff_id_upper)
      call wrap_get_var_realx(ncid, keff_id_upper, k17_upper)

      ! Load K coefficients, interval 18
      filename = trim(dirk)//trim(K18_file)
      call getfil(filename, locfn, 0)
      call wrap_open(locfn, 0, ncid)
      call wrap_inq_varid(ncid, 'KEFF_LOWER', keff_id_lower)
      call wrap_get_var_realx(ncid, keff_id_lower, k18_lower)
      call wrap_inq_varid(ncid, 'KEFF_UPPER', keff_id_upper)
      call wrap_get_var_realx(ncid, keff_id_upper, k18_upper)

      ! Load K coefficients, interval 19
      filename = trim(dirk)//trim(K19_file)
      call getfil(filename, locfn, 0)
      call wrap_open(locfn, 0, ncid)
      call wrap_inq_varid(ncid, 'KEFF_LOWER', keff_id_lower)
      call wrap_get_var_realx(ncid, keff_id_lower, k19_lower)
      call wrap_inq_varid(ncid, 'KEFF_UPPER', keff_id_upper)
      call wrap_get_var_realx(ncid, keff_id_upper, k19_upper)

      ! Load K coefficients, interval 20
      filename = trim(dirk)//trim(K20_file)
      call getfil(filename, locfn, 0)
      call wrap_open(locfn, 0, ncid)
      call wrap_inq_varid(ncid, 'KEFF_LOWER', keff_id_lower)
      call wrap_get_var_realx(ncid, keff_id_lower, k20_lower)
      call wrap_inq_varid(ncid, 'KEFF_UPPER', keff_id_upper)
      call wrap_get_var_realx(ncid, keff_id_upper, k20_upper)

      ! Load K coefficients, interval 21
      filename = trim(dirk)//trim(K21_file)
      call getfil(filename, locfn, 0)
      call wrap_open(locfn, 0, ncid)
      call wrap_inq_varid(ncid, 'KEFF_LOWER', keff_id_lower)
      call wrap_get_var_realx(ncid, keff_id_lower, k21_lower)
      call wrap_inq_varid(ncid, 'KEFF_UPPER', keff_id_upper)
      call wrap_get_var_realx(ncid, keff_id_upper, k21_upper)

      ! Load K coefficients, interval 22
      filename = trim(dirk)//trim(K22_file)
      call getfil(filename, locfn, 0)
      call wrap_open(locfn, 0, ncid)
      call wrap_inq_varid(ncid, 'KEFF_LOWER', keff_id_lower)
      call wrap_get_var_realx(ncid, keff_id_lower, k22_lower)
      call wrap_inq_varid(ncid, 'KEFF_UPPER', keff_id_upper)
      call wrap_get_var_realx(ncid, keff_id_upper, k22_upper)

      ! Load K coefficients, interval 23
      filename = trim(dirk)//trim(K23_file)
      call getfil(filename, locfn, 0)
      call wrap_open(locfn, 0, ncid)
      call wrap_inq_varid(ncid, 'KEFF_LOWER', keff_id_lower)
      call wrap_get_var_realx(ncid, keff_id_lower, k23_lower)
      call wrap_inq_varid(ncid, 'KEFF_UPPER', keff_id_upper)
      call wrap_get_var_realx(ncid, keff_id_upper, k23_upper)

      ! Load K coefficients, interval 24
      filename = trim(dirk)//trim(K24_file)
      call getfil(filename, locfn, 0)
      call wrap_open(locfn, 0, ncid)
      call wrap_inq_varid(ncid, 'KEFF_LOWER', keff_id_lower)
      call wrap_get_var_realx(ncid, keff_id_lower, k24_lower)
      call wrap_inq_varid(ncid, 'KEFF_UPPER', keff_id_upper)
      call wrap_get_var_realx(ncid, keff_id_upper, k24_upper)

      ! Load K coefficients, interval 25
      filename = trim(dirk)//trim(K25_file)
      call getfil(filename, locfn, 0)
      call wrap_open(locfn, 0, ncid)
      call wrap_inq_varid(ncid, 'KEFF_LOWER', keff_id_lower)
      call wrap_get_var_realx(ncid, keff_id_lower, k25_lower)
      call wrap_inq_varid(ncid, 'KEFF_UPPER', keff_id_upper)
      call wrap_get_var_realx(ncid, keff_id_upper, k25_upper)

      ! Load K coefficients, water vapor self continuum
      filename = trim(dirct)//trim(kh2oself_file)
      call getfil(filename, locfn, 0)
      call wrap_open(locfn, 0, ncid)
      call wrap_inq_varid(ncid, 'KSELF', keff_id)
      call wrap_get_var_realx(ncid, keff_id, kh2oself)

      ! Load K coefficients, for co2 continuum
      filename = trim(dirct)//trim(kco2cont_file)
      call getfil(filename, locfn, 0)
      call wrap_open(locfn, 0, ncid)
      call wrap_inq_varid(ncid, 'KSELF', keff_id)
      call wrap_get_var_realx(ncid, keff_id, kco2cont_8gpt)

      ! Load absorptions coefficients, for n2n2 continuum
      filename = trim(dirci)//trim(kn2n2cia_file )
      call getfil(filename, locfn, 0)
      call wrap_open(locfn, 0, ncid)
      call wrap_inq_varid(ncid, 'sigma', keff_id)
      call wrap_get_var_realx(ncid, keff_id, kn2n2 ) 
 
      ! Load K coefficients, for h2n2 continuum
      filename = trim(dirci)//trim(kh2n2cia_file )
      call getfil(filename, locfn, 0)
      call wrap_open(locfn, 0, ncid)
      call wrap_inq_varid(ncid, 'sigma', keff_id)
      call wrap_get_var_realx(ncid, keff_id, kh2n2 )

      ! Load K coefficients, for h2h2 continuum
      filename = trim(dirci)//trim(kh2h2cia_file )
      call getfil(filename, locfn, 0)
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



end module initialize_rad_mod_1D
