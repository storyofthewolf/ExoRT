
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

    filename = trim(exort_rootdir)//trim(dirk_nh3)//trim(k_nh3_file)
    call getfil(filename, locfn, 0)
    call cam_pio_openfile(ncid, locfn, PIO_NOWRITE)
    ierr =  pio_inq_varid(ncid, 'data',   keff_id)
    ierr =  pio_get_var(ncid, keff_id, k_nh3)
    call pio_closefile(ncid)

    filename = trim(exort_rootdir)//trim(dirk_co)//trim(k_co_file)
    call getfil(filename, locfn, 0)
    call cam_pio_openfile(ncid, locfn, PIO_NOWRITE)
    ierr =  pio_inq_varid(ncid, 'data',   keff_id)
    ierr =  pio_get_var(ncid, keff_id, k_co)
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

    filename = trim(exort_rootdir)//trim(dirci)//trim(ko2o2cia_file )
    call getfil(filename, locfn, 0)
    call cam_pio_openfile(ncid, locfn, PIO_NOWRITE)
    ierr =  pio_inq_varid(ncid, 'sigma',   keff_id)
    ierr =  pio_get_var(ncid, keff_id, ko2o2)
    call pio_closefile(ncid)

    filename = trim(exort_rootdir)//trim(dirci)//trim(ko2n2cia_file )
    call getfil(filename, locfn, 0)
    call cam_pio_openfile(ncid, locfn, PIO_NOWRITE)
    ierr =  pio_inq_varid(ncid, 'sigma',   keff_id)
    ierr =  pio_get_var(ncid, keff_id, ko2n2)
    call pio_closefile(ncid)

    filename = trim(exort_rootdir)//trim(dirci)//trim(ko2co2cia_file )
    call getfil(filename, locfn, 0)
    call cam_pio_openfile(ncid, locfn, PIO_NOWRITE)
    ierr =  pio_inq_varid(ncid, 'sigma',   keff_id)
    ierr =  pio_get_var(ncid, keff_id, ko2co2)
    call pio_closefile(ncid)



! broadcast optical constants to all nodes
#if ( defined SPMD )
    call mpibcast(k_h2o,  ntot_wavlnrng*ngauss_8gpt*kc_npress*kc_ntemp, mpir8, 0, mpicom)
    call mpibcast(k_co2,  ntot_wavlnrng*ngauss_8gpt*kc_npress*kc_ntemp, mpir8, 0, mpicom)
    call mpibcast(k_ch4,  ntot_wavlnrng*ngauss_8gpt*kc_npress*kc_ntemp, mpir8, 0, mpicom)
    call mpibcast(k_c2h6, ntot_wavlnrng*ngauss_8gpt*kc_npress*kc_ntemp, mpir8, 0, mpicom)
    call mpibcast(k_o2,   ntot_wavlnrng*ngauss_8gpt*kc_npress*kc_ntemp, mpir8, 0, mpicom)
    call mpibcast(k_o3,   ntot_wavlnrng*ngauss_8gpt*kc_npress*kc_ntemp, mpir8, 0, mpicom)
    call mpibcast(k_nh3,  ntot_wavlnrng*ngauss_8gpt*kc_npress*kc_ntemp, mpir8, 0, mpicom)
    call mpibcast(k_co,   ntot_wavlnrng*ngauss_8gpt*kc_npress*kc_ntemp, mpir8, 0, mpicom)

    call mpibcast(kh2oself_mtckd, ngauss_8gpt*ntot_wavlnrng*kmtckd_ntemp, mpir8, 0, mpicom)
    call mpibcast(kh2ofrgn_mtckd, ngauss_8gpt*ntot_wavlnrng*kmtckd_ntemp, mpir8, 0, mpicom)

    call mpibcast(kn2n2, ntot_wavlnrng*kn2n2_ntemp, mpir8, 0, mpicom)
    call mpibcast(kn2h2, ntot_wavlnrng*kn2h2_ntemp, mpir8, 0, mpicom)
    call mpibcast(kh2h2, ntot_wavlnrng*kh2h2_ntemp, mpir8, 0, mpicom)

    call mpibcast(kco2co2_sw, ntot_wavlnrng*kco2co2_sw_ntemp, mpir8, 0, mpicom)
    call mpibcast(kco2co2_lw, ntot_wavlnrng*kco2co2_lw_ntemp, mpir8, 0, mpicom)
    call mpibcast(kco2h2, ntot_wavlnrng*kco2h2_ntemp, mpir8, 0, mpicom)
    call mpibcast(kco2ch4, ntot_wavlnrng*kco2ch4_ntemp, mpir8, 0, mpicom)

    call mpibcast(ko2o2,  ntot_wavlnrng*ko2o2_ntemp,  mpir8, 0, mpicom)
    call mpibcast(ko2n2,  ntot_wavlnrng*ko2n2_ntemp,  mpir8, 0, mpicom)
    call mpibcast(ko2co2, ntot_wavlnrng*ko2co2_ntemp, mpir8, 0, mpicom)
    

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
    integer :: r_id
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
    ierr =  pio_inq_varid(ncid, 'rel_grid',   r_id)
    ierr =  pio_get_var(ncid, r_id, rel_grid)
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
    ierr =  pio_inq_varid(ncid, 'rei_grid',   r_id)
    ierr =  pio_get_var(ncid, r_id, rei_grid)
    ierr =  pio_inq_varid(ncid, 'Qext_ice',   q_id)
    ierr =  pio_get_var(ncid, q_id, Qcldice)
    ierr =  pio_inq_varid(ncid, 'W_ice',   w_id)
    ierr =  pio_get_var(ncid, w_id, Wcldice)
    ierr =  pio_inq_varid(ncid, 'G_ice',   g_id)
    ierr =  pio_get_var(ncid, g_id, Gcldice)
    call pio_closefile(ncid)

    ! Load CO2 ice cloud optics file
    ! (variable names differ from the H2O files: radii / Qext / W / G;
    !  rei_co2_grid has no compiled-in DATA fallback, so this read is required)
    filename = trim(exort_rootdir)//trim(dircld)//trim(cldoptsICO2_file)
    call getfil(filename, locfn, 0)
    call cam_pio_openfile(ncid, locfn, PIO_NOWRITE)
    ierr =  pio_inq_varid(ncid, 'radii',   r_id)
    ierr =  pio_get_var(ncid, r_id, rei_co2_grid)
    ierr =  pio_inq_varid(ncid, 'Qext',   q_id)
    ierr =  pio_get_var(ncid, q_id, Qcldice_co2)
    ierr =  pio_inq_varid(ncid, 'W',   w_id)
    ierr =  pio_get_var(ncid, w_id, Wcldice_co2)
    ierr =  pio_inq_varid(ncid, 'G',   g_id)
    ierr =  pio_get_var(ncid, g_id, Gcldice_co2)
    call pio_closefile(ncid)

! broadcast water cloud optical constants to all nodes
#if ( defined SPMD )
    call mpibcast(rel_grid, nrel, mpir8, 0, mpicom)
    call mpibcast(Qcldliq, nrel*ntot_wavlnrng, mpir8, 0, mpicom)
    call mpibcast(Wcldliq, nrel*ntot_wavlnrng, mpir8, 0, mpicom)
    call mpibcast(Gcldliq, nrel*ntot_wavlnrng, mpir8, 0, mpicom)
#endif

! broadcast ice cloud optical constants to all nodes
#if ( defined SPMD )
    call mpibcast(rei_grid, nrei, mpir8, 0, mpicom)
    call mpibcast(Qcldice, nrei*ntot_wavlnrng, mpir8, 0, mpicom)
    call mpibcast(Wcldice, nrei*ntot_wavlnrng, mpir8, 0, mpicom)
    call mpibcast(Gcldice, nrei*ntot_wavlnrng, mpir8, 0, mpicom)
#endif

! broadcast CO2 ice cloud optical constants to all nodes
#if ( defined SPMD )
    call mpibcast(rei_co2_grid, nrei_co2, mpir8, 0, mpicom)
    call mpibcast(Qcldice_co2, nrei_co2*ntot_wavlnrng, mpir8, 0, mpicom)
    call mpibcast(Wcldice_co2, nrei_co2*ntot_wavlnrng, mpir8, 0, mpicom)
    call mpibcast(Gcldice_co2, nrei_co2*ntot_wavlnrng, mpir8, 0, mpicom)
#endif

  end subroutine initialize_cldopts


!============================================================================

  subroutine initialize_hazeopts

!------------------------------------------------------------------------
!
! Purpose:  Initialize the CARMA haze aerosol optical constants from the
!   optics file (PIO twin of the 1-D loader in initialize_rad_mod_1D.F90).
!   Loads the pre-tabulated haze optics (mass extinction Kext,
!   single-scattering albedo W, asymmetry G) on the exact CARMA element/bin
!   grid into the radgrid module arrays (kcarma/wcarma/gcarma). Call only
!   when the haze path is enabled (do_exo_haze); otherwise the haze optics
!   arrays are never read and the aerosol path contributes no opacity.
!
!------------------------------------------------------------------------

#if ( defined SPMD)
  use mpishorthand
#endif

    use ioFileMod, only: getfil
    use cam_pio_utils, only: cam_pio_openfile
    use pio,  only: pio_inq_varid, pio_get_var, pio_closefile, pio_nowrite,  &
                    file_desc_t, var_desc_t, pio_inq_dimid, pio_inquire_dimension
    use abortutils, only: endrun

    implicit none

!------------------------------------------------------------------------
!
! Local Variables
!
    type(file_desc_t) :: ncid
    integer :: k_id, w_id, g_id
    integer :: dim_id, n_elem_file, n_bin_file, n_wav_file
    character(len=256) :: locfn
    character(len=256) :: filename
    integer :: ierr

!------------------------------------------------------------------------
!
! Start Code
!
    if ( masterproc ) then
      write (6, '(2x, a)') '_______________________________________________________'
      write (6, '(2x, a)') '____________ read in haze optical properties __________'
      write (6, '(2x, a)') '_______________________________________________________'
    endif

    filename = trim(exort_rootdir)//trim(diraer)//trim(hazeopts_file)
    call getfil(filename, locfn, 0)
    call cam_pio_openfile(ncid, locfn, PIO_NOWRITE)

    ! The optics tables are used on the exact CARMA element/bin grid with no
    ! interpolation, so the file dimensions must match the compiled sizes.
    ierr = pio_inq_dimid(ncid, 'nelements', dim_id)
    ierr = pio_inquire_dimension(ncid, dim_id, len=n_elem_file)
    ierr = pio_inq_dimid(ncid, 'nbins', dim_id)
    ierr = pio_inquire_dimension(ncid, dim_id, len=n_bin_file)
    ierr = pio_inq_dimid(ncid, 'nwavlrng', dim_id)
    ierr = pio_inquire_dimension(ncid, dim_id, len=n_wav_file)
    if (n_elem_file /= nelem_carma .or. n_bin_file /= nbin_carma &
        .or. n_wav_file /= ntot_wavlnrng) then
      if (masterproc) then
        write(6,*) 'initialize_hazeopts: optics file dims (nelements,nbins,nwavlrng) =', &
                   n_elem_file, n_bin_file, n_wav_file, &
                   ' expected', nelem_carma, nbin_carma, ntot_wavlnrng
      endif
      call endrun('initialize_hazeopts: haze optics file dimensions do not match compiled sizes')
    endif

    ierr =  pio_inq_varid(ncid, 'Kext',   k_id)
    ierr =  pio_get_var(ncid, k_id, kcarma)
    ierr =  pio_inq_varid(ncid, 'W',   w_id)
    ierr =  pio_get_var(ncid, w_id, wcarma)
    ierr =  pio_inq_varid(ncid, 'G',   g_id)
    ierr =  pio_get_var(ncid, g_id, gcarma)
    call pio_closefile(ncid)

    ! Kext is stored in [cm2 g-1]; convert to [m2 kg-1] for the tau kernel
    kcarma(:,:,:) = kcarma(:,:,:)*0.1

! broadcast haze optical constants to all nodes
#if ( defined SPMD )
    call mpibcast(kcarma, nelem_carma*nbin_carma*ntot_wavlnrng, mpir8, 0, mpicom)
    call mpibcast(wcarma, nelem_carma*nbin_carma*ntot_wavlnrng, mpir8, 0, mpicom)
    call mpibcast(gcarma, nelem_carma*nbin_carma*ntot_wavlnrng, mpir8, 0, mpicom)
#endif

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



end module initialize_rad_mod_cam
