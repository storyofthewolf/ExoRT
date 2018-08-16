
module radgrid

! version n68h2o
!-----------------------------------------------------------
! Purpose:  Stores basic radiation module data for use
!           by radiation.F90, mcica.F90, rad_interp_mod.F90
!-----------------------------------------------------------
  
  use shr_kind_mod,      only: r8 => shr_kind_r8

  implicit none
  public


  ! Define number of wavelength intervals
  integer, parameter  :: ntot_wavlnrng = 68        ! total # of wv intervals
  integer, parameter :: ngauss_8gpt = 8
  integer, parameter :: ngauss_16gpt = 16
 
  ! <wavenum_edge_>x refers to the wavenumbers (1/wavelength) at the wavelength
  !  interval edges (note that there are 'ntot_wavlnrng+1' edges for each
  !  wavelength group) [cm^-1]:
  real(r8), dimension(ntot_wavlnrng+1) :: wavenum_edge
  data wavenum_edge /  &          ! all wavenumber edges [cm-1]
            0.00E+00,   40.00000,   100.0000,   160.0000, &
            220.0000,   280.0000,   330.0000,   380.0000, &
            440.0000,   495.0000,   545.0000,   617.0000, &
            667.0000,   720.0000,   800.0000,   875.0000, &
            940.0000,   1000.000,   1065.000,   1108.000, &
            1200.000,   1275.000,   1350.000,   1450.000, &
            1550.000,   1650.000,   1750.000,   1850.000, &
            1950.000,   2050.000,   2200.000,   2397.000, &
            2494.000,   2796.000,   3087.000,   3425.000, &
            3760.000,   4030.000,   4540.000,   4950.000, &
            5370.000,   5925.000,   6390.000,   6990.000, &
            7650.000,   8315.000,   8850.000,   9350.000, &
            9650.000,   10400.00,   11220.00,   11870.00, &
            12790.00,   13300.00,   14470.00,   15000.00, &
            16000.00,   16528.00,   17649.00,   18198.00, &
            18518.00,   22222.00,   25641.00,   29308.00, &
            30376.00,   32562.00,   35087.00,   36363.00, &
            42087.00 /
 
  real(r8), dimension(ntot_wavlnrng) :: wavenum_mid
  data wavenum_mid / &    ! all wavenumber midpoints
       20.,     70.,     130.,      190.,     250.,      305., &
       355.,    410.,    467.5,     520.,     581.,      642., &
       693.5,   760.,    837.5,     907.5,    970.,      1032.5, &
       1086.5,  1154.,   1237.5,    1312.5,   1400.,     1500., &
       1600.,   1700.,   1800.,     1900.,    2000.,     2125., &
       2298.5,  2445.5,  2645.,     2941.5,   3256.,     3592.5, &
       3895.,   4285.,   4745.,     5160.,    5647.5,    6157.5, &
       6690.,   7320.,   7982.5,    8582.50,  9100.,     9500., &
       10025.,  10810.,  11545.,    12330.,   13045.,    13885., &
       14735.,  15500.,  16264.,    17088.5,  17923.5,   18358., &
       20370.,  23931.5, 27474.5,   29842.,   31469.,    33824.5, &
       35725.,  39225. /

  !
  ! guass point information
  !
  real(r8), dimension(8) :: g_xpos_edge_8gpt
  real(r8), dimension(8) :: g_weight_8gpt
  real(r8), dimension(16) :: g_xpos_edge_16gpt
  real(r8), dimension(16) :: g_weight_16gpt

  ! For 8 gauss pts
  ! "x"-positions for Gaussian quadrature within each wavelength band [none]
  !  (there are 'ngauss_pts' of these):
  data g_xpos_edge_8gpt / 0.00000, 0.30192, 0.57571, 0.79583, 0.94178, 0.98890, 0.99576, 0.99939 /

  ! Weights for Gaussian quadrature within each wavelength band [none] (there are
  !  'ngauss_pts' of these):
  data g_weight_8gpt / 0.30192, 0.27379, 0.22012, 0.14595, 0.04712, 0.00686, 0.00363, 0.00061 /

  ! For 16 gauss pts
  data g_xpos_edge_16gpt / 0.00000, 0.15275, 0.30192, 0.44402, 0.57571, 0.6939, 0.79583, 0.87911, &
                          0.94178, 0.98427, 0.9889, 0.99273, 0.99576, 0.99798, 0.99939, 0.99993 /

  data g_weight_16gpt / 0.15275, 0.14917, 0.14210, 0.13169, 0.11819, 0.10193, 0.08328, 0.06267, &
                       0.04249, 0.00463, 0.00383, 0.00303, 0.00222, 0.00141, 0.00054, 0.00007 /



  ! Gauss point gridding
  integer, parameter  :: ng1 = 8
  integer, parameter  :: ng2 = 8
  integer, parameter  :: ng3 = 8
  integer, parameter  :: ng4 = 8
  integer, parameter  :: ng5 = 8
  integer, parameter  :: ng6 = 8
  integer, parameter  :: ng7 = 8
  integer, parameter  :: ng8 = 8
  integer, parameter  :: ng9 = 8
  integer, parameter  :: ng10 = 8
  integer, parameter  :: ng11 = 8
  integer, parameter  :: ng12 = 8
  integer, parameter  :: ng13 = 8
  integer, parameter  :: ng14 = 8
  integer, parameter  :: ng15 = 8
  integer, parameter  :: ng16 = 8
  integer, parameter  :: ng17 = 8
  integer, parameter  :: ng18 = 8
  integer, parameter  :: ng19 = 8
  integer, parameter  :: ng20 = 8
  integer, parameter  :: ng21 = 8
  integer, parameter  :: ng22 = 8
  integer, parameter  :: ng23 = 8
  integer, parameter  :: ng24 = 8
  integer, parameter  :: ng25 = 8
  integer, parameter  :: ng26 = 8
  integer, parameter  :: ng27 = 8
  integer, parameter  :: ng28 = 8
  integer, parameter  :: ng29 = 8
  integer, parameter  :: ng30 = 8
  integer, parameter  :: ng31 = 8
  integer, parameter  :: ng32 = 8
  integer, parameter  :: ng33 = 8
  integer, parameter  :: ng34 = 8
  integer, parameter  :: ng35 = 8
  integer, parameter  :: ng36 = 8
  integer, parameter  :: ng37 = 8
  integer, parameter  :: ng38 = 8
  integer, parameter  :: ng39 = 8
  integer, parameter  :: ng40 = 8
  integer, parameter  :: ng41 = 8
  integer, parameter  :: ng42 = 8
  integer, parameter  :: ng43 = 8
  integer, parameter  :: ng44 = 8
  integer, parameter  :: ng45 = 8
  integer, parameter  :: ng46 = 8
  integer, parameter  :: ng47 = 8
  integer, parameter  :: ng48 = 8
  integer, parameter  :: ng49 = 8
  integer, parameter  :: ng50 = 8
  integer, parameter  :: ng51 = 8
  integer, parameter  :: ng52 = 8
  integer, parameter  :: ng53 = 8
  integer, parameter  :: ng54 = 8
  integer, parameter  :: ng55 = 8
  integer, parameter  :: ng56 = 8
  integer, parameter  :: ng57 = 8
  integer, parameter  :: ng58 = 8
  integer, parameter  :: ng59 = 8
  integer, parameter  :: ng60 = 8
  integer, parameter  :: ng61 = 8
  integer, parameter  :: ng62 = 8
  integer, parameter  :: ng63 = 8
  integer, parameter  :: ng64 = 8
  integer, parameter  :: ng65 = 8
  integer, parameter  :: ng66 = 8
  integer, parameter  :: ng67 = 8
  integer, parameter  :: ng68 = 8

  ! Define number of gauss points in each spectral interval
  integer, dimension(ntot_wavlnrng) :: ngauss_pts           ! # of Gauss quad pts per interval
  data ngauss_pts / ng1,  ng2,  ng3,  ng4,  ng5,  ng6,  ng7,  ng8,  ng9,  ng10, &
                    ng11, ng12, ng13, ng14, ng15, ng16, ng17, ng18, ng19, ng20, &
                    ng21, ng22, ng23, ng24, ng25, ng26, ng27, ng28, ng29, ng30, &
                    ng31, ng32, ng33, ng34, ng35, ng36, ng37, ng38, ng39, ng40, &
                    ng41, ng42, ng43, ng44, ng45, ng46, ng47, ng48, ng49, ng50, &
                    ng51, ng52, ng53, ng54, ng55, ng56, ng57, ng58, ng59, ng60, &
                    ng61, ng62, ng63, ng64, ng65, ng66, ng67, ng68  /


  real(r8),dimension(ntot_wavlnrng) :: solarflux
  real(r8) :: S0


  ! Dimensions for spectral interval gauss points
  integer, parameter  :: ntot_gpt &                 ! total # of probability intervals
                         = ng1 + ng2 + ng3 + ng4 + ng5 + ng6 + ng7 + ng8  &
                         + ng9 + ng10 + ng11 + ng12 + ng13 + ng14 + ng15  &
                         + ng16 + ng17 + ng18 + ng19 + ng20 + ng21 + ng22 &
                         + ng23 + ng24 + ng25 + ng26 + ng27 + ng28 + ng29 &
                         + ng30 + ng31 + ng32 + ng33 + ng34 + ng35 + ng36 + ng37 + ng38 + ng39 &
                         + ng40 + ng41 + ng42 + ng43 + ng44 + ng45 + ng46 + ng47 + ng48 + ng49 &
                         + ng50 + ng51 + ng52 + ng53 + ng54 + ng55 + ng56 + ng57 + ng58 + ng59 &
                         + ng60 + ng61 + ng62 + ng63 + ng64 + ng65 + ng66 + ng67 + ng68 
            

  ! Dimensions of current k-coefficient datasets
  integer, parameter  :: kc_npress = 61     ! # of reference pressure
  integer, parameter  :: kc_ntemp = 17      ! # of reference temps

  ! Pressure, temperature, and species weight grids
  real(r8), dimension(kc_npress) :: log10pgrid
  real(r8), dimension(kc_npress) :: pgrid
  real(r8), dimension(kc_ntemp) :: tgrid

  ! log 10 pressure grid whole atmosphere
  data log10pgrid / -2.0, -1.9, -1.8, -1.7, -1.6, -1.5, -1.4, -1.3, -1.2, -1.1, -1.0, -0.9, -0.8, -0.7, &
                    -0.6, -0.5, -0.4, -0.3, -0.2, -0.1, 0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, &
                     1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.1, 2.2, 2.3, 2.4, 2.5, 2.6, &
                     2.7, 2.8, 2.9, 3.0, 3.1, 3.2, 3.3, 3.4, 3.5, 3.6, 3.7, 3.8, 3.9, 4.0 /


  ! pressure grid [mb] whole atmopshere ; for old grid
  data pgrid / 0.0100000, 0.0125893, 0.0158489, 0.0199526, 0.0251189, 0.0316228, 0.0398107, 0.0501187, &
               0.0630957, 0.0794328, 0.100000,  0.125893,  0.158489,  0.199526,  0.251189,  0.316228,  &
               0.398107,  0.501187,  0.630957,  0.794328,  1.00000,   1.25893,   1.58489,   1.99526,   &
               2.51189,   3.16228,   3.98107,   5.01187,   6.30958,   7.94328,   10.0000,   12.5893,   &
               15.8489,   19.9526,   25.1189,   31.6228,   39.8107,   50.1187,   63.0958,   79.4328,   &
               100.000,   125.893,   158.489,   199.526,   251.189,   316.228,   398.107,   501.188,   &
               630.958,   794.328,   1000.00,   1258.93,   1584.89,   1995.26,   2511.89,   3162.28,   &
               3981.07,   5011.87,   6309.57,   7943.28,   10000.0 /


  ! temperature grids [K]
  data tgrid / 100, 125, 150, 175, 200, 225, 250, 275, 300, 325, 350, 375, 400, 425, 450, 475, 500  /

  !!
  !! Water Vapor Self Continuum parameters  
  !!

  !! ==== mtckd definitions ===
  integer, parameter  :: ks_ntemp = 8       ! # of reference temperatures
  integer, parameter  :: ngH2O = 1         ! # of gauss points per spectral interval in
  integer, parameter  :: ks_npress = 1
  integer, parameter  :: ks_nweight = 1

  ! T, P, W grid for H2O continuum
  real(r8), dimension(ks_npress) :: pgrid_self
  real(r8), dimension(ks_npress) :: log10pgrid_self
  real(r8), dimension(ks_ntemp) :: tgrid_self
  real(r8), dimension(ks_nweight) :: wgrid_self

  ! temperature grid [K], water vapor self continuum
  data tgrid_self / 150, 200, 250, 300, 350, 400, 450, 500 /

  ! h2o amount grid [vmr]
  data wgrid_self / 1.0 /

  ! water vapor self continuum data array for mtckd
  real(r8), dimension(ntot_wavlnrng, ks_ntemp) :: kh2oself_mtckd
  !! ==== end mtckd definitions ===

  !!!=== bps continuum definitions ===
  real(r8), dimension(ntot_wavlnrng) :: self
  real(r8), dimension(ntot_wavlnrng) :: foreign
  real(r8), dimension(ntot_wavlnrng) :: base_self
  real(r8), dimension(ntot_wavlnrng) :: base_foreign
  real(r8), dimension(ntot_wavlnrng) :: TempCoeff
  !!!=== end bps definitions ====


  ! N2-N2 CIA temperature
  integer, parameter  :: kn2n2_ntemp = 10       ! # of reference temperatures
  real(r8), dimension(kn2n2_ntemp) :: tgrid_n2n2
  data tgrid_n2n2 / 40.0, 51.7, 66.7, 86.2, 111.3, 143.8, 185.7, 239.8, 309.7, 400.0 /

  ! H2-N2 CIA temperature
  integer, parameter  :: kh2n2_ntemp = 10       ! # of reference temperatures
  real(r8), dimension(kh2n2_ntemp) :: tgrid_h2n2
  data tgrid_h2n2 / 40.0, 51.7, 66.7, 86.2, 111.3, 143.8, 185.7, 239.8, 309.7, 400.0 /

  ! H2-H2 CIA temperature
  integer, parameter  :: kh2h2_ntemp = 11       ! # of reference temperatures
  real(r8), dimension(kh2h2_ntemp) :: tgrid_h2h2
  data tgrid_h2h2 / 200.0, 225.0, 250.0, 275.0, 300.0, 325.0, 350.0, 375.0, 400.0, 425.0, 450.0 /


  !
  ! K coefficient variables
  ! 

  real(r8), dimension(ng1, kc_npress, kc_ntemp) :: k01_data
  real(r8), dimension(ng2, kc_npress, kc_ntemp) :: k02_data
  real(r8), dimension(ng3, kc_npress, kc_ntemp) :: k03_data
  real(r8), dimension(ng4, kc_npress, kc_ntemp) :: k04_data
  real(r8), dimension(ng5, kc_npress, kc_ntemp) :: k05_data
  real(r8), dimension(ng6, kc_npress, kc_ntemp) :: k06_data
  real(r8), dimension(ng7, kc_npress, kc_ntemp) :: k07_data
  real(r8), dimension(ng8, kc_npress, kc_ntemp) :: k08_data
  real(r8), dimension(ng9, kc_npress, kc_ntemp) :: k09_data
  real(r8), dimension(ng10, kc_npress, kc_ntemp) :: k10_data
  real(r8), dimension(ng11, kc_npress, kc_ntemp) :: k11_data
  real(r8), dimension(ng12, kc_npress, kc_ntemp) :: k12_data
  real(r8), dimension(ng13, kc_npress, kc_ntemp) :: k13_data
  real(r8), dimension(ng14, kc_npress, kc_ntemp) :: k14_data
  real(r8), dimension(ng15, kc_npress, kc_ntemp) :: k15_data
  real(r8), dimension(ng16, kc_npress, kc_ntemp) :: k16_data
  real(r8), dimension(ng17, kc_npress, kc_ntemp) :: k17_data
  real(r8), dimension(ng18, kc_npress, kc_ntemp) :: k18_data
  real(r8), dimension(ng19, kc_npress, kc_ntemp) :: k19_data
  real(r8), dimension(ng20, kc_npress, kc_ntemp) :: k20_data
  real(r8), dimension(ng21, kc_npress, kc_ntemp) :: k21_data
  real(r8), dimension(ng22, kc_npress, kc_ntemp) :: k22_data
  real(r8), dimension(ng23, kc_npress, kc_ntemp) :: k23_data
  real(r8), dimension(ng24, kc_npress, kc_ntemp) :: k24_data
  real(r8), dimension(ng25, kc_npress, kc_ntemp) :: k25_data
  real(r8), dimension(ng26, kc_npress, kc_ntemp) :: k26_data
  real(r8), dimension(ng27, kc_npress, kc_ntemp) :: k27_data
  real(r8), dimension(ng28, kc_npress, kc_ntemp) :: k28_data
  real(r8), dimension(ng29, kc_npress, kc_ntemp) :: k29_data
  real(r8), dimension(ng30, kc_npress, kc_ntemp) :: k30_data
  real(r8), dimension(ng31, kc_npress, kc_ntemp) :: k31_data
  real(r8), dimension(ng32, kc_npress, kc_ntemp) :: k32_data
  real(r8), dimension(ng33, kc_npress, kc_ntemp) :: k33_data
  real(r8), dimension(ng34, kc_npress, kc_ntemp) :: k34_data
  real(r8), dimension(ng35, kc_npress, kc_ntemp) :: k35_data
  real(r8), dimension(ng36, kc_npress, kc_ntemp) :: k36_data
  real(r8), dimension(ng37, kc_npress, kc_ntemp) :: k37_data
  real(r8), dimension(ng38, kc_npress, kc_ntemp) :: k38_data
  real(r8), dimension(ng39, kc_npress, kc_ntemp) :: k39_data
  real(r8), dimension(ng40, kc_npress, kc_ntemp) :: k40_data
  real(r8), dimension(ng41, kc_npress, kc_ntemp) :: k41_data
  real(r8), dimension(ng42, kc_npress, kc_ntemp) :: k42_data
  real(r8), dimension(ng43, kc_npress, kc_ntemp) :: k43_data
  real(r8), dimension(ng44, kc_npress, kc_ntemp) :: k44_data
  real(r8), dimension(ng45, kc_npress, kc_ntemp) :: k45_data
  real(r8), dimension(ng46, kc_npress, kc_ntemp) :: k46_data
  real(r8), dimension(ng47, kc_npress, kc_ntemp) :: k47_data
  real(r8), dimension(ng48, kc_npress, kc_ntemp) :: k48_data
  real(r8), dimension(ng49, kc_npress, kc_ntemp) :: k49_data
  real(r8), dimension(ng50, kc_npress, kc_ntemp) :: k50_data
  real(r8), dimension(ng51, kc_npress, kc_ntemp) :: k51_data
  real(r8), dimension(ng52, kc_npress, kc_ntemp) :: k52_data
  real(r8), dimension(ng53, kc_npress, kc_ntemp) :: k53_data
  real(r8), dimension(ng54, kc_npress, kc_ntemp) :: k54_data
  real(r8), dimension(ng55, kc_npress, kc_ntemp) :: k55_data
  real(r8), dimension(ng56, kc_npress, kc_ntemp) :: k56_data
  real(r8), dimension(ng57, kc_npress, kc_ntemp) :: k57_data
  real(r8), dimension(ng58, kc_npress, kc_ntemp) :: k58_data
  real(r8), dimension(ng59, kc_npress, kc_ntemp) :: k59_data
  real(r8), dimension(ng60, kc_npress, kc_ntemp) :: k60_data
  real(r8), dimension(ng61, kc_npress, kc_ntemp) :: k61_data
  real(r8), dimension(ng62, kc_npress, kc_ntemp) :: k62_data
  real(r8), dimension(ng63, kc_npress, kc_ntemp) :: k63_data
  real(r8), dimension(ng64, kc_npress, kc_ntemp) :: k64_data
  real(r8), dimension(ng65, kc_npress, kc_ntemp) :: k65_data
  real(r8), dimension(ng66, kc_npress, kc_ntemp) :: k66_data
  real(r8), dimension(ng67, kc_npress, kc_ntemp) :: k67_data
  real(r8), dimension(ng68, kc_npress, kc_ntemp) :: k68_data

  
  !! Continuum Data !!
  !
  ! CO2 continuum information, self and foreign
  !
  integer, parameter   :: ngCO2 = 8                           ! # of gauss points per spectral interval
  real(r8), parameter :: pref_co2 = 1013.0
  real(r8), parameter :: vref_co2 = 0.1 ! reference CO2 for continuum calculation
  real(r8), dimension(ngCO2,ntot_wavlnrng) :: kco2cont_8gpt  ! CO2 continuum data from file
  real(r8), dimension(ntot_gpt) :: kco2cont  ! CO2 continuum data, reduced to gauss pt grid
  !!!!real(r8), dimension(ntot_wavlnrng, ks_ntemp) :: kco2cont


  ! CIA absorption data from HITRAN
  real(r8), dimension(ntot_wavlnrng,kh2h2_ntemp) :: kh2h2  ! H2-H2 CIA data [cm-1 amagat-2]
  real(r8), dimension(ntot_wavlnrng,kh2n2_ntemp) :: kh2n2  ! H2-N2 CIA data [cm-1 amagat-2]
  real(r8), dimension(ntot_wavlnrng,kn2n2_ntemp) :: kn2n2  ! N2-N2 CIA data [cm-1 amagat-2]

 ! 8 to 16 gauss point mapping
  integer, dimension(16) :: map8to16gpt

  data map8to16gpt / 1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6, 7, 7, 8, 8 /



  !
  ! Cloud optics parameters  
  !
  integer, parameter :: ncld_grp = 2     ! number of cloud groups, 1=water clouds, 2=ice clouds
  integer, parameter :: nrel = 30        ! number of radii grid points, liquid
  integer, parameter :: nrei = 300        ! number of radii grid points, ice
  real(r8), parameter :: cldmin = 1.0d-80
  real(r8), dimension(nrel) :: rel_grid
  real(r8), dimension(nrei) :: rei_grid

  ! data structures containing cloud optics data from file
  real(r8) :: Qcldliq(nrel, ntot_wavlnrng)
  real(r8) :: Wcldliq(nrel, ntot_wavlnrng)
  real(r8) :: Gcldliq(nrel, ntot_wavlnrng)

  real(r8) :: Qcldice(nrei, ntot_wavlnrng)
  real(r8) :: Wcldice(nrei, ntot_wavlnrng)
  real(r8) :: Gcldice(nrei, ntot_wavlnrng)


  data rel_grid / 1., 2., 3., 4., 5., 6., 7., 8., 9., 10., &
                  11., 12., 13., 14., 15., 16., 17., 18., 19., 20., &
                  21., 22., 23., 24., 25., 26., 27., 28., 29., 30.  /

  data rei_grid / 1., 2., 3., 4., 5., 6., 7., 8., 9., 10., &
                  11., 12., 13., 14., 15., 16., 17., 18., 19., 20., &
                  21., 22., 23., 24., 25., 26., 27., 28., 29., 30., &
                  31., 32., 33., 34., 35., 36., 37., 38., 39., 40., &
                  41., 42., 43., 44., 45., 46., 47., 48., 49., 50., &
                  51., 52., 53., 54., 55., 56., 57., 58., 59., 60., &
                  61., 62., 63., 64., 65., 66., 67., 68., 69., 70., &
                  71., 72., 73., 74., 75., 76., 77., 78., 79., 80., &
                  81., 82., 83., 84., 85., 86., 87., 88., 89., 90., &
                  91., 92., 93., 94., 95., 96., 97., 98., 99., 100., &
                  101., 102., 103., 104., 105., 106., 107., 108., 109., 110., &
                  111., 112., 113., 114., 115., 116., 117., 118., 119., 120., &
                  121., 122., 123., 124., 125., 126., 127., 128., 129., 130., &
                  131., 132., 133., 134., 135., 136., 137., 138., 139., 140., &
                  141., 142., 143., 144., 145., 146., 147., 148., 149., 150., &
                  151., 152., 153., 154., 155., 156., 157., 158., 159., 160., &
                  161., 162., 163., 164., 165., 166., 167., 168., 169., 170., &
                  171., 172., 173., 174., 175., 176., 177., 178., 179., 180., &
                  181., 182., 183., 184., 185., 186., 187., 188., 189., 190., &
                  191., 192., 193., 194., 195., 196., 197., 198., 199., 200., &
                  201., 202., 203., 204., 205., 206., 207., 208., 209., 210., &
                  211., 212., 213., 214., 215., 216., 217., 218., 219., 220., &
                  221., 222., 223., 224., 225., 226., 227., 228., 229., 230., &
                  231., 232., 233., 234., 235., 236., 237., 238., 239., 240., &
                  241., 242., 243., 244., 245., 246., 247., 248., 249., 250., &
                  251., 252., 253., 254., 255., 256., 257., 258., 259., 260., &
                  261., 262., 263., 264., 265., 266., 267., 268., 269., 270., &
                  271., 272., 273., 274., 275., 276., 277., 278., 279., 280., &
                  281., 282., 283., 284., 285., 286., 287., 288., 289., 290., &
                  291., 292., 293., 294., 295., 296., 297., 298., 299., 300. /



  !! top level indices for longwave and cloud calculations
  integer, public :: camtop  ! top level to solve using CAM CKD RT.
  integer, public :: ntopcld ! top level to solve for cloud overlap
  integer, public :: nlevsRT ! number of levels to calculate RT

! For MCICA calculation, must match ng* specified above
  integer, dimension(ntot_gpt) :: ngb
  data ngb / 1,1,1,1,1,1,1,1, &
             2,2,2,2,2,2,2,2, &
             3,3,3,3,3,3,3,3, &
             4,4,4,4,4,4,4,4, &
             5,5,5,5,5,5,5,5, &
             6,6,6,6,6,6,6,6, &
             7,7,7,7,7,7,7,7, &
             8,8,8,8,8,8,8,8, &
             9,9,9,9,9,9,9,9, &
             10,10,10,10,10,10,10,10, &
             11,11,11,11,11,11,11,11, &
             12,12,12,12,12,12,12,12, &
             13,13,13,13,13,13,13,13, &
             14,14,14,14,14,14,14,14, &
             15,15,15,15,15,15,15,15, &
             16,16,16,16,16,16,16,16, &
             17,17,17,17,17,17,17,17, &
             18,18,18,18,18,18,18,18, &
             19,19,19,19,19,19,19,19, &
             20,20,20,20,20,20,20,20, &
             21,21,21,21,21,21,21,21, &
             22,22,22,22,22,22,22,22, &
             23,23,23,23,23,23,23,23, &
             24,24,24,24,24,24,24,24, &
             25,25,25,25,25,25,25,25, &
             26,26,26,26,26,26,26,26, &
             27,27,27,27,27,27,27,27, &
             28,28,28,28,28,28,28,28, &
             29,29,29,29,29,29,29,29, &
             30,30,30,30,30,30,30,30, &
             31,31,31,31,31,31,31,31, &
             32,32,32,32,32,32,32,32, &
             33,33,33,33,33,33,33,33, &
             34,34,34,34,34,34,34,34, &
             35,35,35,35,35,35,35,35, &
             36,36,36,36,36,36,36,36, &
             37,37,37,37,37,37,37,37, &
             38,38,38,38,38,38,38,38, &
             39,39,39,39,39,39,39,39, &
             40,40,40,40,40,40,40,40, &
             41,41,41,41,41,41,41,41, &
             42,42,42,42,42,42,42,42, &
             43,43,43,43,43,43,43,43, &
             44,44,44,44,44,44,44,44, &
             45,45,45,45,45,45,45,45, &
             46,46,46,46,46,46,46,46, &
             47,47,47,47,47,47,47,47, &
             48,48,48,48,48,48,48,48, &
             49,49,49,49,49,49,49,49, &
             50,50,50,50,50,50,50,50, &
             51,51,51,51,51,51,51,51, &
             52,52,52,52,52,52,52,52, &
             53,53,53,53,53,53,53,53, &
             54,54,54,54,54,54,54,54, &
             55,55,55,55,55,55,55,55, &
             56,56,56,56,56,56,56,56, &
             57,57,57,57,57,57,57,57, &
             58,58,58,58,58,58,58,58, &
             59,59,59,59,59,59,59,59, &
             60,60,60,60,60,60,60,60, &
             61,61,61,61,61,61,61,61, &
             62,62,62,62,62,62,62,62, &
             63,63,63,63,63,63,63,63, &
             64,64,64,64,64,64,64,64, &
             65,65,65,65,65,65,65,65, &
             66,66,66,66,66,66,66,66, &
             67,67,67,67,67,67,67,67, &
             68,68,68,68,68,68,68,68 /


end module radgrid
