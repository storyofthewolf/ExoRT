
module initialize_rad_mod_1D

! version n42h2o
!
! read in and initialize radiaive transfer grids
!

use kabs
use exoplanet_mod, only: solar_file, dirsol
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
    character(len=256) :: locfn, filename

!------------------------------------------------------------------------
!
! Start Code
!
!    if ( masterproc ) then

      write (6, '(2x, a)') '_______________________________________________________'
      write (6, '(2x, a)') '_______ read in gas correlated-k coefficients _________'
      write (6, '(2x, a)') '_______________________________________________________'

      !----  Load K coefficients ----
      !----  H2O, CO2, CH4, C2H6  ----
      filename = trim(exort_rootdir)//trim(dirk_h2o)//trim(k_h2o_file)
      call getfil(filename, locfn, 0)
      call wrap_open(locfn, 0, ncid)
      call wrap_inq_varid(ncid, 'data', keff_id)
      call wrap_get_var_realx(ncid, keff_id, k_h2o)

      filename = trim(exort_rootdir)//trim(dirk_co2)//trim(k_co2_file)
      call getfil(filename, locfn, 0)
      call wrap_open(locfn, 0, ncid)
      call wrap_inq_varid(ncid, 'data', keff_id)
      call wrap_get_var_realx(ncid, keff_id, k_co2)

      filename = trim(exort_rootdir)//trim(dirk_ch4)//trim(k_ch4_file)
      call getfil(filename, locfn, 0)
      call wrap_open(locfn, 0, ncid)
      call wrap_inq_varid(ncid, 'data', keff_id)
      call wrap_get_var_realx(ncid, keff_id, k_ch4)

      filename = trim(exort_rootdir)//trim(dirk_c2h6)//trim(k_c2h6_file)
      call getfil(filename, locfn, 0)
      call wrap_open(locfn, 0, ncid)
      call wrap_inq_varid(ncid, 'data', keff_id)
      call wrap_get_var_realx(ncid, keff_id, k_c2h6)

      filename = trim(exort_rootdir)//trim(dirk_o3)//trim(k_o3_file)
      call getfil(filename, locfn, 0)
      call wrap_open(locfn, 0, ncid)
      call wrap_inq_varid(ncid, 'data', keff_id)
      call wrap_get_var_realx(ncid, keff_id, k_o3)

      filename = trim(exort_rootdir)//trim(dirk_o2)//trim(k_o2_file)
      call getfil(filename, locfn, 0)
      call wrap_open(locfn, 0, ncid)
      call wrap_inq_varid(ncid, 'data', keff_id)
      call wrap_get_var_realx(ncid, keff_id, k_o2)



      write (6, '(2x, a)') '_______________________________________________________'
      write (6, '(2x, a)') '____________ read in water vapor continuum ____________'
      write (6, '(2x, a)') '_______________________________________________________'
      ! Load water vapor continuum
      !! mtckd
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
      !! /mtckd

      write (6, '(2x, a)') '_______________________________________________________'
      write (6, '(2x, a)') '_____________________ read in CIA  ____________________'
      write (6, '(2x, a)') '_______________________________________________________'
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

      filename = trim(exort_rootdir)//trim(dirci)//trim(ko2o2cia_file )
      call getfil(filename, locfn, 0)
      call wrap_open(locfn, 0, ncid)
      call wrap_inq_varid(ncid, 'sigma', keff_id)
      call wrap_get_var_realx(ncid, keff_id, ko2o2 )

      filename = trim(exort_rootdir)//trim(dirci)//trim(ko2n2cia_file )
      call getfil(filename, locfn, 0)
      call wrap_open(locfn, 0, ncid)
      call wrap_inq_varid(ncid, 'sigma', keff_id)
      call wrap_get_var_realx(ncid, keff_id, ko2n2 )

      filename = trim(exort_rootdir)//trim(dirci)//trim(ko2co2cia_file )
      call getfil(filename, locfn, 0)
      call wrap_open(locfn, 0, ncid)
      call wrap_inq_varid(ncid, 'sigma', keff_id)
      call wrap_get_var_realx(ncid, keff_id, ko2co2 )

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
