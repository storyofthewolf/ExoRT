
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
  public :: initialize_cldopts
  public :: initialize_hazeopts
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
      call check_kfile_grid(ncid, filename)
      call wrap_inq_varid(ncid, 'data', keff_id)
      call wrap_get_var_realx(ncid, keff_id, k_h2o)

      filename = trim(exort_rootdir)//trim(dirk_co2)//trim(k_co2_file)
      call getfil(filename, locfn, 0)
      call wrap_open(locfn, 0, ncid)
      call check_kfile_grid(ncid, filename)
      call wrap_inq_varid(ncid, 'data', keff_id)
      call wrap_get_var_realx(ncid, keff_id, k_co2)

      filename = trim(exort_rootdir)//trim(dirk_ch4)//trim(k_ch4_file)
      call getfil(filename, locfn, 0)
      call wrap_open(locfn, 0, ncid)
      call check_kfile_grid(ncid, filename)
      call wrap_inq_varid(ncid, 'data', keff_id)
      call wrap_get_var_realx(ncid, keff_id, k_ch4)

      filename = trim(exort_rootdir)//trim(dirk_c2h6)//trim(k_c2h6_file)
      call getfil(filename, locfn, 0)
      call wrap_open(locfn, 0, ncid)
      call check_kfile_grid(ncid, filename)
      call wrap_inq_varid(ncid, 'data', keff_id)
      call wrap_get_var_realx(ncid, keff_id, k_c2h6)

      filename = trim(exort_rootdir)//trim(dirk_nh3)//trim(k_nh3_file)
      call getfil(filename, locfn, 0)
      call wrap_open(locfn, 0, ncid)
      call check_kfile_grid(ncid, filename)
      call wrap_inq_varid(ncid, 'data', keff_id)
      call wrap_get_var_realx(ncid, keff_id, k_nh3)

      filename = trim(exort_rootdir)//trim(dirk_co)//trim(k_co_file)
      call getfil(filename, locfn, 0)
      call wrap_open(locfn, 0, ncid)
      call check_kfile_grid(ncid, filename)
      call wrap_inq_varid(ncid, 'data', keff_id)
      call wrap_get_var_realx(ncid, keff_id, k_co)

      filename = trim(exort_rootdir)//trim(dirk_o3)//trim(k_o3_file)
      call getfil(filename, locfn, 0)
      call wrap_open(locfn, 0, ncid)
      call check_kfile_grid(ncid, filename)
      call wrap_inq_varid(ncid, 'data', keff_id)
      call wrap_get_var_realx(ncid, keff_id, k_o3)

      filename = trim(exort_rootdir)//trim(dirk_o2)//trim(k_o2_file)
      call getfil(filename, locfn, 0)
      call wrap_open(locfn, 0, ncid)
      call check_kfile_grid(ncid, filename)
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

      write (6, '(2x, a)') '_______________________________________________________'
      write (6, '(2x, a)') '_________________ initializing solar file _____________'
      write (6, '(2x, a)') '_______________________________________________________'

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
! Purpose:  Initialize the cloud optical constants from input files.
!   Loads the H2O liquid and ice Mie optical-property tables (extinction
!   efficiency Q, single-scattering albedo W, asymmetry G) and their effective-
!   radius grids into the radgrid module arrays. Called only when do_exo_clouds
!   is enabled; otherwise the cloud optics arrays stay zero and the cloud path
!   contributes no opacity.
!
!------------------------------------------------------------------------

    use ioFileMod, only: getfil
    use cloud,     only: dircld, cldoptsL_file, cldoptsI_file, cldoptsICO2_file

    implicit none
    include 'netcdf.inc'

!------------------------------------------------------------------------
!
! Local Variables
!
    integer :: ncid
    integer :: q_id, w_id, g_id, r_id
    character(len=256) :: locfn, filename

!------------------------------------------------------------------------
!
! Start Code
!
      write (6, '(2x, a)') '_______________________________________________________'
      write (6, '(2x, a)') '___________ read in cloud optical properties __________'
      write (6, '(2x, a)') '_______________________________________________________'

      ! ---- H2O liquid cloud optics ----
      filename = trim(exort_rootdir)//trim(dircld)//trim(cldoptsL_file)
      call getfil(filename, locfn, 0)
      call wrap_open(locfn, 0, ncid)
      call wrap_inq_varid(ncid, 'rel_grid', r_id)
      call wrap_get_var_realx(ncid, r_id, rel_grid)
      call wrap_inq_varid(ncid, 'Qext_liq', q_id)
      call wrap_get_var_realx(ncid, q_id, Qcldliq)
      call wrap_inq_varid(ncid, 'W_liq', w_id)
      call wrap_get_var_realx(ncid, w_id, Wcldliq)
      call wrap_inq_varid(ncid, 'G_liq', g_id)
      call wrap_get_var_realx(ncid, g_id, Gcldliq)

      ! ---- H2O ice cloud optics ----
      filename = trim(exort_rootdir)//trim(dircld)//trim(cldoptsI_file)
      call getfil(filename, locfn, 0)
      call wrap_open(locfn, 0, ncid)
      call wrap_inq_varid(ncid, 'rei_grid', r_id)
      call wrap_get_var_realx(ncid, r_id, rei_grid)
      call wrap_inq_varid(ncid, 'Qext_ice', q_id)
      call wrap_get_var_realx(ncid, q_id, Qcldice)
      call wrap_inq_varid(ncid, 'W_ice', w_id)
      call wrap_get_var_realx(ncid, w_id, Wcldice)
      call wrap_inq_varid(ncid, 'G_ice', g_id)
      call wrap_get_var_realx(ncid, g_id, Gcldice)

      ! ---- CO2 ice cloud optics ----
      ! (variable/dim names differ from the H2O files: radii / Qext / W / G)
      filename = trim(exort_rootdir)//trim(dircld)//trim(cldoptsICO2_file)
      call getfil(filename, locfn, 0)
      call wrap_open(locfn, 0, ncid)
      call wrap_inq_varid(ncid, 'radii', r_id)
      call wrap_get_var_realx(ncid, r_id, rei_co2_grid)
      call wrap_inq_varid(ncid, 'Qext', q_id)
      call wrap_get_var_realx(ncid, q_id, Qcldice_co2)
      call wrap_inq_varid(ncid, 'W', w_id)
      call wrap_get_var_realx(ncid, w_id, Wcldice_co2)
      call wrap_inq_varid(ncid, 'G', g_id)
      call wrap_get_var_realx(ncid, g_id, Gcldice_co2)

  end subroutine initialize_cldopts


!============================================================================

  subroutine initialize_hazeopts

!------------------------------------------------------------------------
!
! Purpose:  Initialize the CARMA haze aerosol optical constants from the
!   input file. Loads the pre-tabulated haze optics (mass extinction Kext,
!   single-scattering albedo W, asymmetry G) on the exact CARMA element/bin
!   grid into the radgrid module arrays. Called only when do_exo_haze is
!   enabled; otherwise the haze optics arrays are never read and the aerosol
!   path contributes no opacity.
!
!------------------------------------------------------------------------

    use ioFileMod, only: getfil
    use cloud,     only: diraer, hazeopts_file

    implicit none
    include 'netcdf.inc'

!------------------------------------------------------------------------
!
! Local Variables
!
    integer :: ncid
    integer :: k_id, w_id, g_id
    integer :: dim_id, n_file
    character(len=256) :: locfn, filename

!------------------------------------------------------------------------
!
! Start Code
!
      write (6, '(2x, a)') '_______________________________________________________'
      write (6, '(2x, a)') '____________ read in haze optical properties __________'
      write (6, '(2x, a)') '_______________________________________________________'

      filename = trim(exort_rootdir)//trim(diraer)//trim(hazeopts_file)
      call getfil(filename, locfn, 0)
      call wrap_open(locfn, 0, ncid)

      ! The optics tables are used on the exact CARMA element/bin grid with no
      ! interpolation, so the file dimensions must match the compiled sizes.
      call wrap_inq_dimid(ncid, 'nelements', dim_id)
      call wrap_inq_dimlen(ncid, dim_id, n_file)
      if (n_file /= nelem_carma) then
        write(6,*) 'initialize_hazeopts: nelements mismatch: file', n_file, 'expected', nelem_carma
        stop
      endif
      call wrap_inq_dimid(ncid, 'nbins', dim_id)
      call wrap_inq_dimlen(ncid, dim_id, n_file)
      if (n_file /= nbin_carma) then
        write(6,*) 'initialize_hazeopts: nbins mismatch: file', n_file, 'expected', nbin_carma
        stop
      endif
      call wrap_inq_dimid(ncid, 'nwavlrng', dim_id)
      call wrap_inq_dimlen(ncid, dim_id, n_file)
      if (n_file /= ntot_wavlnrng) then
        write(6,*) 'initialize_hazeopts: nwavlrng mismatch: file', n_file, 'expected', ntot_wavlnrng
        stop
      endif

      call wrap_inq_varid(ncid, 'Kext', k_id)
      call wrap_get_var_realx(ncid, k_id, kcarma)
      call wrap_inq_varid(ncid, 'W', w_id)
      call wrap_get_var_realx(ncid, w_id, wcarma)
      call wrap_inq_varid(ncid, 'G', g_id)
      call wrap_get_var_realx(ncid, g_id, gcarma)

      ! Kext is stored in [m2 kg-1], which is what the tau kernel wants, so no
      ! unit conversion is applied here. (The generator computes 3Q/(4 rho r)
      ! in cm2 g-1 and divides by 10 before writing.) Some older haze files
      ! carry a stale "cm2 g-1" units attribute that does not match their own
      ! values -- trust the values, not that attribute.

  end subroutine initialize_hazeopts


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

  subroutine check_kfile_grid(ncid, fname)

!------------------------------------------------------------------------
!
! Purpose:  Verify that an open k-coefficient file sits on the compiled
!           grid. The 'data' dimensions must match exactly. The Temperature,
!           Pressure [mb] and GaussWeights (g-interval midpoints) coordinates
!           are compared against tgrid, pgrid and the radgrid g-intervals;
!           any mismatch stops the run. The lookup never reads these
!           coordinates, so a table on the wrong grid is otherwise silent --
!           the 2020-2026 HITRAN-2016 H2O table carried Temperature =
!           100,100,125..475 and was read as 100..500.
!           A coordinate that is absent or all zero cannot be verified and
!           also stops the run. SpectralBands holds band indices only, so
!           bands are checked by count (data dimension).
!
!------------------------------------------------------------------------

    implicit none
    include 'netcdf.inc'

    integer, intent(in) :: ncid
    character(len=*), intent(in) :: fname

    integer :: keff_id, ndims, i, dimlen
    integer, dimension(4) :: dimids, expect
    real(r8), dimension(ngauss_8gpt) :: gmid

    ! Fortran order of the file's (NTemp, NPress, NGauss, NBins)
    expect = (/ ntot_wavlnrng, ngauss_8gpt, kc_npress, kc_ntemp /)
    call wrap_inq_varid(ncid, 'data', keff_id)
    if (nf_inq_varndims(ncid, keff_id, ndims) /= NF_NOERR .or. ndims /= 4) then
      write(6,*) 'check_kfile_grid: ', trim(fname), ': data is not 4-D'
      stop 1
    endif
    if (nf_inq_vardimid(ncid, keff_id, dimids) /= NF_NOERR) then
      write(6,*) 'check_kfile_grid: ', trim(fname), ': cannot read data dimensions'
      stop 1
    endif
    do i = 1, 4
      call wrap_inq_dimlen(ncid, dimids(i), dimlen)
      if (dimlen /= expect(i)) then
        write(6,*) 'check_kfile_grid: ', trim(fname), ': data dimension', i, &
                   'is', dimlen, 'expected', expect(i), '(bands, gauss, press, temp)'
        stop 1
      endif
    enddo

    gmid(:) = g_xpos_edge_8gpt(:) + 0.5_r8*g_weight_8gpt(:)
    call check_kfile_coord(ncid, fname, 'Temperature',  tgrid, 1.0e-3_r8, .false.)
    call check_kfile_coord(ncid, fname, 'Pressure',     pgrid, 1.0e-4_r8, .true.)
    call check_kfile_coord(ncid, fname, 'GaussWeights', gmid,  1.0e-5_r8, .false.)

  end subroutine check_kfile_grid

!====================================================================================

  subroutine check_kfile_coord(ncid, fname, vname, expected, tol, relative)

!------------------------------------------------------------------------
!
! Purpose:  Compare one 1-D coordinate variable of a k-coefficient file
!           against the compiled grid (see check_kfile_grid).
!
!------------------------------------------------------------------------

    implicit none
    include 'netcdf.inc'

    integer, intent(in) :: ncid
    character(len=*), intent(in) :: fname, vname
    real(r8), dimension(:), intent(in) :: expected
    real(r8), intent(in) :: tol
    logical, intent(in) :: relative

    integer :: vid, ndims, dimid, dimlen, i
    real(r8), dimension(size(expected)) :: vals, err

    if (nf_inq_varid(ncid, vname, vid) /= NF_NOERR) then
      write(6,*) 'check_kfile_grid: ', trim(fname), ': no ', vname, ' coordinate, grid cannot be verified'
      stop 1
    endif
    if (nf_inq_varndims(ncid, vid, ndims) /= NF_NOERR .or. ndims /= 1) then
      write(6,*) 'check_kfile_grid: ', trim(fname), ': ', vname, ' is not 1-D'
      stop 1
    endif
    if (nf_inq_vardimid(ncid, vid, dimid) /= NF_NOERR) then
      write(6,*) 'check_kfile_grid: ', trim(fname), ': cannot read ', vname, ' dimension'
      stop 1
    endif
    call wrap_inq_dimlen(ncid, dimid, dimlen)
    if (dimlen /= size(expected)) then
      write(6,*) 'check_kfile_grid: ', trim(fname), ': ', vname, ' has', dimlen, &
                 'values, expected', size(expected)
      stop 1
    endif
    call wrap_get_var_realx(ncid, vid, vals)

    if (all(vals == 0.0_r8)) then
      write(6,*) 'check_kfile_grid: ', trim(fname), ': ', vname, ' coordinate is all zero, grid cannot be verified'
      stop 1
    endif

    err(:) = abs(vals(:) - expected(:))
    if (relative) err(:) = err(:) / abs(expected(:))
    if (any(err > tol)) then
      write(6,*) 'check_kfile_grid: ', trim(fname), ': ', vname, ' does not match the compiled grid'
      do i = 1, size(expected)
        if (err(i) > tol) write(6,*) '   index', i, ' file', vals(i), ' expected', expected(i)
      enddo
      stop 1
    endif

  end subroutine check_kfile_coord

!====================================================================================

end module initialize_rad_mod_1D
