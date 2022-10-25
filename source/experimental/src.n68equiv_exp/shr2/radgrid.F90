
module radgrid

! version n68h2o
!-----------------------------------------------------------
! Purpose:  Stores basic radiation module data for use
!           by radiation.F90, mcica.F90, rad_interp_mod.F90
!-----------------------------------------------------------
  
  use shr_kind_mod,      only: r8 => shr_kind_r8

  implicit none
  public

  ! set maxmimum values to limit numerical issues
  real(r8) :: taumax=400.0d0
  real(r8) :: taumax_co2cld=40.0d0
  real(r8) :: ssamax_co2cld=1.0d0

  ! Define number of wavelength intervals
  integer, parameter  :: ntot_wavlnrng = 68        ! total # of wv intervals
  integer, parameter :: ngauss_8gpt = 8
  integer, parameter :: ngauss_16gpt = 16
  integer, parameter :: ngpt_max = ngauss_8gpt     ! used for array definitions
 
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

  !Gaussian quadrate gauss bins and weights, same as used in Clima
!  data g_xpos_edge_8gpt /  0.0000000, 0.16523105, 0.47499999, 0.78476894,0.94999999,0.95869637, 0.97500002, 0.99130368 / 
!  data g_weight_8gpt   / 0.16523105, 0.30976894, 0.30976894, 0.16523105, 0.0086963773, 0.016303658, 0.016303658, 0.0086963177 / 


  ! For 16 gauss pts
!  data g_xpos_edge_16gpt / 0.00000, 0.15275, 0.30192, 0.44402, 0.57571, 0.6939, 0.79583, 0.87911, &
!                          0.94178, 0.98427, 0.9889, 0.99273, 0.99576, 0.99798, 0.99939, 0.99993 /

!  data g_weight_16gpt / 0.15275, 0.14917, 0.14210, 0.13169, 0.11819, 0.10193, 0.08328, 0.06267, &
!                       0.04249, 0.00463, 0.00383, 0.00303, 0.00222, 0.00141, 0.00054, 0.00007 /


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
 ! integer, parameter  :: kmtckd_ntemp = 9       ! # of reference temperatures
 integer, parameter  :: kmtckd_ntemp = 41       ! # of reference temperatures
!  integer, parameter  :: ks_npress = 1
!  integer, parameter  :: ks_nweight = 1

  ! T, P, W grid for H2O continuum
  real(r8), dimension(kmtckd_ntemp) :: tgrid_mtckd

  ! temperature grid [K], water vapor self continuum
!  data tgrid_mtckd / 100, 150, 200, 250, 300, 350, 400, 450, 500 /
  data tgrid_mtckd / 100, 110, 120, 130, 140, 150, 160, 170, 180, 190, &
                     200, 210, 220, 230, 240, 250, 260, 270, 280, 290, &
                     300, 310, 320, 330, 340, 350, 360, 370, 380, 390, &
                     400, 410, 420, 430, 440, 450, 460, 470, 480, 490, 500 /
                    


  ! water vapor self continuum data array for mtckd
  real(r8), dimension(ntot_wavlnrng, kmtckd_ntemp) :: kh2oself_avg_mtckd
  real(r8), dimension(ngauss_8gpt, ntot_wavlnrng, kmtckd_ntemp) :: kh2oself_mtckd

  ! water vapor frgn continuum data array for mtckd
  real(r8), dimension(ntot_wavlnrng, kmtckd_ntemp) :: kh2ofrgn_avg_mtckd
  real(r8), dimension(ngauss_8gpt, ntot_wavlnrng, kmtckd_ntemp) :: kh2ofrgn_mtckd
  !! ==== end mtckd definitions ===

  !!!=== bps continuum definitions ===
  real(r8), dimension(ntot_wavlnrng) :: self
  real(r8), dimension(ntot_wavlnrng) :: foreign
  real(r8), dimension(ntot_wavlnrng) :: base_self
  real(r8), dimension(ntot_wavlnrng) :: base_foreign
  real(r8), dimension(ntot_wavlnrng) :: TempCoeff
  !!!=== end bps definitions ====


  ! CIA grids

  ! H2O-H2O CIA temperature grid
  integer, parameter  :: kh2oh2o_ntemp = 21       ! # of reference temperatures 
  real(r8), dimension(kh2oh2o_ntemp) :: tgrid_h2oh2o
  data tgrid_h2oh2o / 250, 260, 270, 280, 290, 300, 310, 320, 330, 340, 350, 360, &
                    370, 380, 390, 400, 410, 420, 430, 440, 450 /

  ! N2-H2O CIA temperature grid
  integer, parameter  :: kh2on2_ntemp = 21       ! # of reference temperatures 
  real(r8), dimension(kh2on2_ntemp) :: tgrid_h2on2
  data tgrid_h2on2 / 250, 260, 270, 280, 290, 300, 310, 320, 330, 340, 350, 360, &
                    370, 380, 390, 400, 410, 420, 430, 440, 450 /


                  
  ! N2-N2 CIA temperature grid
  integer, parameter  :: kn2n2_ntemp = 10       ! # of reference temperatures         
  real(r8), dimension(kn2n2_ntemp) :: tgrid_n2n2
  data tgrid_n2n2 / 40.0, 51.7, 66.7, 86.2, 111.3, 143.8, 185.7, 239.8, 309.7, 400.0 /

  ! H2-H2 CIA temperature grid
  integer, parameter  :: kh2h2_ntemp = 113       ! # of reference temperatures 
  real(r8), dimension(kh2h2_ntemp) :: tgrid_h2h2
  data tgrid_h2h2 / 200.0, 225.0, 250.0, 275.0, 300.0, 325.0, 350.0, 375.0, 400.0, 425.0, 450.0, 475.0,  &
                    500.0, 525.0, 550.0, 575.0, 600.0, 625.0, 650.0, 675.0, 700.0, 725.0, 750.0, 775.0,  &
                    800.0, 825.0, 850.0, 875.0, 900.0, 925.0, 950.0, 975.0, 1000.0,1025.0,1050.0,1075.0,  &
                    1100.0,1125.0,1150.0,1175.0,1200.0,1225.0,1250.0,1275.0,1300.0,1325.0,1350.0,1375.0,  &
                    1400.0,1425.0,1450.0,1475.0,1500.0,1525.0,1550.0,1575.0,1600.0,1625.0,1650.0,1675.0,  &
                    1700.0,1725.0,1750.0,1775.0,1800.0,1825.0,1850.0,1875.0,1900.0,1925.0,1950.0,1975.0,  &
                    2000.0,2025.0,2050.0,2075.0,2100.0,2125.0,2150.0,2175.0,2200.0,2225.0,2250.0,2275.0,  &
                    2300.0,2325.0,2350.0,2375.0,2400.0,2425.0,2450.0,2475.0,2500.0,2525.0,2550.0,2575.0,  &
                    2600.0,2625.0,2650.0,2675.0,2700.0,2725.0,2750.0,2775.0,2800.0,2825.0,2850.0,2875.0,  &
                    2900.0,2925.0,2950.0,2975.0,3000 /

  ! N2-H2 CIA temperature grid
  integer, parameter  :: kn2h2_ntemp = 10       ! # of reference temperatures         
  real(r8), dimension(kn2h2_ntemp) :: tgrid_n2h2
  data tgrid_n2h2 / 40.0, 51.7, 66.7, 86.2, 111.3, 143.8, 185.7, 239.8, 309.7, 400.0 /

  ! CO2-CO2 CIA temperature grid
  ! Note separate grids for lw and sw contributions 
  integer, parameter  :: kco2co2_lw_ntemp = 10       ! # of reference temperatures         
  real(r8), dimension(kco2co2_lw_ntemp) :: tgrid_co2co2_lw
  data tgrid_co2co2_lw /200.0, 224.0, 250.0, 280.0, 313.0, 350.0, 430.0, 529.0, 651.0, 800.0 /

  integer, parameter  :: kco2co2_sw_ntemp = 3       ! # of reference temperatures         
  real(r8), dimension(kco2co2_sw_ntemp) :: tgrid_co2co2_sw
  data tgrid_co2co2_sw / 221.0, 235.0, 297.0 /

  ! CO2-CH4 CIA temperature grid
  integer, parameter  :: kco2ch4_ntemp = 6       ! # of reference temperatures
  real(r8), dimension(kco2ch4_ntemp) :: tgrid_co2ch4
  data tgrid_co2ch4 / 100.0, 200.0, 300.0, 400.0, 500.0, 600.0 /

  ! CO2-H2 CIA temperature grid
  integer, parameter  :: kco2h2_ntemp = 6       ! # of reference temperatures
  real(r8), dimension(kco2h2_ntemp) :: tgrid_co2h2
  data tgrid_co2h2 / 100.0, 200.0, 300.0, 400.0, 500.0, 600.0 /

  ! Gas gases for line absoprtion 
  integer, parameter  :: nspecies = 3  !H2O, CO2, CH4
  ! gas list
  integer, parameter :: iH2O = 1
  integer, parameter :: iCO2 = 2
  integer, parameter :: iCH4 = 3
  character(len=32), dimension(nspecies), parameter :: & 
             gas_name = (/'H2O','CO2','CH4'/)
   
  integer, parameter :: nalpha = 2

  ! Correlated-K coefficient Arrays
  !
  ! These are defined for each spectral interval separetaly to allow
  ! instituting a variable number of species per bin, or different 
  ! numbers of gauss points, etc.
  !

  ! 
  ! "major" gas correlated-k arrays.  Updated for every grid box, at every timestep.  
  real(r8), dimension(nspecies, ng1, kc_npress, kc_ntemp) :: k01_major_data
  real(r8), dimension(nspecies, ng2, kc_npress, kc_ntemp) :: k02_major_data
  real(r8), dimension(nspecies, ng3, kc_npress, kc_ntemp) :: k03_major_data
  real(r8), dimension(nspecies, ng4, kc_npress, kc_ntemp) :: k04_major_data
  real(r8), dimension(nspecies, ng5, kc_npress, kc_ntemp) :: k05_major_data
  real(r8), dimension(nspecies, ng6, kc_npress, kc_ntemp) :: k06_major_data
  real(r8), dimension(nspecies, ng7, kc_npress, kc_ntemp) :: k07_major_data
  real(r8), dimension(nspecies, ng8, kc_npress, kc_ntemp) :: k08_major_data
  real(r8), dimension(nspecies, ng9, kc_npress, kc_ntemp) :: k09_major_data
  real(r8), dimension(nspecies, ng10, kc_npress, kc_ntemp) :: k10_major_data
  real(r8), dimension(nspecies, ng11, kc_npress, kc_ntemp) :: k11_major_data
  real(r8), dimension(nspecies, ng12, kc_npress, kc_ntemp) :: k12_major_data
  real(r8), dimension(nspecies, ng13, kc_npress, kc_ntemp) :: k13_major_data
  real(r8), dimension(nspecies, ng14, kc_npress, kc_ntemp) :: k14_major_data
  real(r8), dimension(nspecies, ng15, kc_npress, kc_ntemp) :: k15_major_data
  real(r8), dimension(nspecies, ng16, kc_npress, kc_ntemp) :: k16_major_data
  real(r8), dimension(nspecies, ng17, kc_npress, kc_ntemp) :: k17_major_data
  real(r8), dimension(nspecies, ng18, kc_npress, kc_ntemp) :: k18_major_data
  real(r8), dimension(nspecies, ng19, kc_npress, kc_ntemp) :: k19_major_data
  real(r8), dimension(nspecies, ng20, kc_npress, kc_ntemp) :: k20_major_data
  real(r8), dimension(nspecies, ng21, kc_npress, kc_ntemp) :: k21_major_data
  real(r8), dimension(nspecies, ng22, kc_npress, kc_ntemp) :: k22_major_data
  real(r8), dimension(nspecies, ng23, kc_npress, kc_ntemp) :: k23_major_data
  real(r8), dimension(nspecies, ng24, kc_npress, kc_ntemp) :: k24_major_data
  real(r8), dimension(nspecies, ng25, kc_npress, kc_ntemp) :: k25_major_data
  real(r8), dimension(nspecies, ng26, kc_npress, kc_ntemp) :: k26_major_data
  real(r8), dimension(nspecies, ng27, kc_npress, kc_ntemp) :: k27_major_data
  real(r8), dimension(nspecies, ng28, kc_npress, kc_ntemp) :: k28_major_data
  real(r8), dimension(nspecies, ng29, kc_npress, kc_ntemp) :: k29_major_data
  real(r8), dimension(nspecies, ng30, kc_npress, kc_ntemp) :: k30_major_data
  real(r8), dimension(nspecies, ng31, kc_npress, kc_ntemp) :: k31_major_data
  real(r8), dimension(nspecies, ng32, kc_npress, kc_ntemp) :: k32_major_data
  real(r8), dimension(nspecies, ng33, kc_npress, kc_ntemp) :: k33_major_data
  real(r8), dimension(nspecies, ng34, kc_npress, kc_ntemp) :: k34_major_data
  real(r8), dimension(nspecies, ng35, kc_npress, kc_ntemp) :: k35_major_data
  real(r8), dimension(nspecies, ng36, kc_npress, kc_ntemp) :: k36_major_data
  real(r8), dimension(nspecies, ng37, kc_npress, kc_ntemp) :: k37_major_data
  real(r8), dimension(nspecies, ng38, kc_npress, kc_ntemp) :: k38_major_data
  real(r8), dimension(nspecies, ng39, kc_npress, kc_ntemp) :: k39_major_data
  real(r8), dimension(nspecies, ng40, kc_npress, kc_ntemp) :: k40_major_data
  real(r8), dimension(nspecies, ng41, kc_npress, kc_ntemp) :: k41_major_data
  real(r8), dimension(nspecies, ng42, kc_npress, kc_ntemp) :: k42_major_data
  real(r8), dimension(nspecies, ng43, kc_npress, kc_ntemp) :: k43_major_data
  real(r8), dimension(nspecies, ng44, kc_npress, kc_ntemp) :: k44_major_data
  real(r8), dimension(nspecies, ng45, kc_npress, kc_ntemp) :: k45_major_data
  real(r8), dimension(nspecies, ng46, kc_npress, kc_ntemp) :: k46_major_data
  real(r8), dimension(nspecies, ng47, kc_npress, kc_ntemp) :: k47_major_data
  real(r8), dimension(nspecies, ng48, kc_npress, kc_ntemp) :: k48_major_data
  real(r8), dimension(nspecies, ng49, kc_npress, kc_ntemp) :: k49_major_data
  real(r8), dimension(nspecies, ng50, kc_npress, kc_ntemp) :: k50_major_data
  real(r8), dimension(nspecies, ng51, kc_npress, kc_ntemp) :: k51_major_data
  real(r8), dimension(nspecies, ng52, kc_npress, kc_ntemp) :: k52_major_data
  real(r8), dimension(nspecies, ng53, kc_npress, kc_ntemp) :: k53_major_data
  real(r8), dimension(nspecies, ng54, kc_npress, kc_ntemp) :: k54_major_data
  real(r8), dimension(nspecies, ng55, kc_npress, kc_ntemp) :: k55_major_data
  real(r8), dimension(nspecies, ng56, kc_npress, kc_ntemp) :: k56_major_data
  real(r8), dimension(nspecies, ng57, kc_npress, kc_ntemp) :: k57_major_data
  real(r8), dimension(nspecies, ng58, kc_npress, kc_ntemp) :: k58_major_data
  real(r8), dimension(nspecies, ng59, kc_npress, kc_ntemp) :: k59_major_data
  real(r8), dimension(nspecies, ng60, kc_npress, kc_ntemp) :: k60_major_data
  real(r8), dimension(nspecies, ng61, kc_npress, kc_ntemp) :: k61_major_data
  real(r8), dimension(nspecies, ng62, kc_npress, kc_ntemp) :: k62_major_data
  real(r8), dimension(nspecies, ng63, kc_npress, kc_ntemp) :: k63_major_data
  real(r8), dimension(nspecies, ng64, kc_npress, kc_ntemp) :: k64_major_data
  real(r8), dimension(nspecies, ng65, kc_npress, kc_ntemp) :: k65_major_data
  real(r8), dimension(nspecies, ng66, kc_npress, kc_ntemp) :: k66_major_data
  real(r8), dimension(nspecies, ng67, kc_npress, kc_ntemp) :: k67_major_data
  real(r8), dimension(nspecies, ng68, kc_npress, kc_ntemp) :: k68_major_data

  ! grey gas correlated-k arrays
  real(r8), dimension(nspecies, kc_npress, kc_ntemp) :: k01_grey_data
  real(r8), dimension(nspecies, kc_npress, kc_ntemp) :: k02_grey_data
  real(r8), dimension(nspecies, kc_npress, kc_ntemp) :: k03_grey_data
  real(r8), dimension(nspecies, kc_npress, kc_ntemp) :: k04_grey_data
  real(r8), dimension(nspecies, kc_npress, kc_ntemp) :: k05_grey_data
  real(r8), dimension(nspecies, kc_npress, kc_ntemp) :: k06_grey_data
  real(r8), dimension(nspecies, kc_npress, kc_ntemp) :: k07_grey_data
  real(r8), dimension(nspecies, kc_npress, kc_ntemp) :: k08_grey_data
  real(r8), dimension(nspecies, kc_npress, kc_ntemp) :: k09_grey_data
  real(r8), dimension(nspecies, kc_npress, kc_ntemp) :: k10_grey_data
  real(r8), dimension(nspecies, kc_npress, kc_ntemp) :: k11_grey_data
  real(r8), dimension(nspecies, kc_npress, kc_ntemp) :: k12_grey_data
  real(r8), dimension(nspecies, kc_npress, kc_ntemp) :: k13_grey_data
  real(r8), dimension(nspecies, kc_npress, kc_ntemp) :: k14_grey_data
  real(r8), dimension(nspecies, kc_npress, kc_ntemp) :: k15_grey_data
  real(r8), dimension(nspecies, kc_npress, kc_ntemp) :: k16_grey_data
  real(r8), dimension(nspecies, kc_npress, kc_ntemp) :: k17_grey_data
  real(r8), dimension(nspecies, kc_npress, kc_ntemp) :: k18_grey_data
  real(r8), dimension(nspecies, kc_npress, kc_ntemp) :: k19_grey_data
  real(r8), dimension(nspecies, kc_npress, kc_ntemp) :: k20_grey_data
  real(r8), dimension(nspecies, kc_npress, kc_ntemp) :: k21_grey_data
  real(r8), dimension(nspecies, kc_npress, kc_ntemp) :: k22_grey_data
  real(r8), dimension(nspecies, kc_npress, kc_ntemp) :: k23_grey_data
  real(r8), dimension(nspecies, kc_npress, kc_ntemp) :: k24_grey_data
  real(r8), dimension(nspecies, kc_npress, kc_ntemp) :: k25_grey_data
  real(r8), dimension(nspecies, kc_npress, kc_ntemp) :: k26_grey_data
  real(r8), dimension(nspecies, kc_npress, kc_ntemp) :: k27_grey_data
  real(r8), dimension(nspecies, kc_npress, kc_ntemp) :: k28_grey_data
  real(r8), dimension(nspecies, kc_npress, kc_ntemp) :: k29_grey_data
  real(r8), dimension(nspecies, kc_npress, kc_ntemp) :: k30_grey_data
  real(r8), dimension(nspecies, kc_npress, kc_ntemp) :: k31_grey_data
  real(r8), dimension(nspecies, kc_npress, kc_ntemp) :: k32_grey_data
  real(r8), dimension(nspecies, kc_npress, kc_ntemp) :: k33_grey_data
  real(r8), dimension(nspecies, kc_npress, kc_ntemp) :: k34_grey_data
  real(r8), dimension(nspecies, kc_npress, kc_ntemp) :: k35_grey_data
  real(r8), dimension(nspecies, kc_npress, kc_ntemp) :: k36_grey_data
  real(r8), dimension(nspecies, kc_npress, kc_ntemp) :: k37_grey_data
  real(r8), dimension(nspecies, kc_npress, kc_ntemp) :: k38_grey_data
  real(r8), dimension(nspecies, kc_npress, kc_ntemp) :: k39_grey_data
  real(r8), dimension(nspecies, kc_npress, kc_ntemp) :: k40_grey_data
  real(r8), dimension(nspecies, kc_npress, kc_ntemp) :: k41_grey_data
  real(r8), dimension(nspecies, kc_npress, kc_ntemp) :: k42_grey_data
  real(r8), dimension(nspecies, kc_npress, kc_ntemp) :: k43_grey_data
  real(r8), dimension(nspecies, kc_npress, kc_ntemp) :: k44_grey_data
  real(r8), dimension(nspecies, kc_npress, kc_ntemp) :: k45_grey_data
  real(r8), dimension(nspecies, kc_npress, kc_ntemp) :: k46_grey_data
  real(r8), dimension(nspecies, kc_npress, kc_ntemp) :: k47_grey_data
  real(r8), dimension(nspecies, kc_npress, kc_ntemp) :: k48_grey_data
  real(r8), dimension(nspecies, kc_npress, kc_ntemp) :: k49_grey_data
  real(r8), dimension(nspecies, kc_npress, kc_ntemp) :: k50_grey_data
  real(r8), dimension(nspecies, kc_npress, kc_ntemp) :: k51_grey_data
  real(r8), dimension(nspecies, kc_npress, kc_ntemp) :: k52_grey_data
  real(r8), dimension(nspecies, kc_npress, kc_ntemp) :: k53_grey_data
  real(r8), dimension(nspecies, kc_npress, kc_ntemp) :: k54_grey_data
  real(r8), dimension(nspecies, kc_npress, kc_ntemp) :: k55_grey_data
  real(r8), dimension(nspecies, kc_npress, kc_ntemp) :: k56_grey_data
  real(r8), dimension(nspecies, kc_npress, kc_ntemp) :: k57_grey_data
  real(r8), dimension(nspecies, kc_npress, kc_ntemp) :: k58_grey_data
  real(r8), dimension(nspecies, kc_npress, kc_ntemp) :: k59_grey_data
  real(r8), dimension(nspecies, kc_npress, kc_ntemp) :: k60_grey_data
  real(r8), dimension(nspecies, kc_npress, kc_ntemp) :: k61_grey_data
  real(r8), dimension(nspecies, kc_npress, kc_ntemp) :: k62_grey_data
  real(r8), dimension(nspecies, kc_npress, kc_ntemp) :: k63_grey_data
  real(r8), dimension(nspecies, kc_npress, kc_ntemp) :: k64_grey_data
  real(r8), dimension(nspecies, kc_npress, kc_ntemp) :: k65_grey_data
  real(r8), dimension(nspecies, kc_npress, kc_ntemp) :: k66_grey_data
  real(r8), dimension(nspecies, kc_npress, kc_ntemp) :: k67_grey_data
  real(r8), dimension(nspecies, kc_npress, kc_ntemp) :: k68_grey_data

  ! individual gas correlated-k arrays.
  real(r8), dimension(ng1, kc_npress, kc_ntemp) :: k01_h2o, k01_co2, k01_ch4
  real(r8), dimension(ng2, kc_npress, kc_ntemp) :: k02_h2o, k02_co2, k02_ch4
  real(r8), dimension(ng3, kc_npress, kc_ntemp) :: k03_h2o, k03_co2, k03_ch4
  real(r8), dimension(ng4, kc_npress, kc_ntemp) :: k04_h2o, k04_co2, k04_ch4
  real(r8), dimension(ng5, kc_npress, kc_ntemp) :: k05_h2o, k05_co2, k05_ch4
  real(r8), dimension(ng6, kc_npress, kc_ntemp) :: k06_h2o, k06_co2, k06_ch4
  real(r8), dimension(ng7, kc_npress, kc_ntemp) :: k07_h2o, k07_co2, k07_ch4
  real(r8), dimension(ng8, kc_npress, kc_ntemp) :: k08_h2o, k08_co2, k08_ch4
  real(r8), dimension(ng9, kc_npress, kc_ntemp) :: k09_h2o, k09_co2, k09_ch4
  real(r8), dimension(ng10, kc_npress, kc_ntemp) :: k10_h2o, k10_co2, k10_ch4
  real(r8), dimension(ng11, kc_npress, kc_ntemp) :: k11_h2o, k11_co2, k11_ch4
  real(r8), dimension(ng12, kc_npress, kc_ntemp) :: k12_h2o, k12_co2, k12_ch4
  real(r8), dimension(ng13, kc_npress, kc_ntemp) :: k13_h2o, k13_co2, k13_ch4
  real(r8), dimension(ng14, kc_npress, kc_ntemp) :: k14_h2o, k14_co2, k14_ch4
  real(r8), dimension(ng15, kc_npress, kc_ntemp) :: k15_h2o, k15_co2, k15_ch4
  real(r8), dimension(ng16, kc_npress, kc_ntemp) :: k16_h2o, k16_co2, k16_ch4
  real(r8), dimension(ng17, kc_npress, kc_ntemp) :: k17_h2o, k17_co2, k17_ch4
  real(r8), dimension(ng18, kc_npress, kc_ntemp) :: k18_h2o, k18_co2, k18_ch4
  real(r8), dimension(ng19, kc_npress, kc_ntemp) :: k19_h2o, k19_co2, k19_ch4
  real(r8), dimension(ng20, kc_npress, kc_ntemp) :: k20_h2o, k20_co2, k20_ch4
  real(r8), dimension(ng21, kc_npress, kc_ntemp) :: k21_h2o, k21_co2, k21_ch4
  real(r8), dimension(ng22, kc_npress, kc_ntemp) :: k22_h2o, k22_co2, k22_ch4
  real(r8), dimension(ng23, kc_npress, kc_ntemp) :: k23_h2o, k23_co2, k23_ch4
  real(r8), dimension(ng24, kc_npress, kc_ntemp) :: k24_h2o, k24_co2, k24_ch4
  real(r8), dimension(ng25, kc_npress, kc_ntemp) :: k25_h2o, k25_co2, k25_ch4
  real(r8), dimension(ng26, kc_npress, kc_ntemp) :: k26_h2o, k26_co2, k26_ch4
  real(r8), dimension(ng27, kc_npress, kc_ntemp) :: k27_h2o, k27_co2, k27_ch4
  real(r8), dimension(ng28, kc_npress, kc_ntemp) :: k28_h2o, k28_co2, k28_ch4
  real(r8), dimension(ng29, kc_npress, kc_ntemp) :: k29_h2o, k29_co2, k29_ch4
  real(r8), dimension(ng30, kc_npress, kc_ntemp) :: k30_h2o, k30_co2, k30_ch4
  real(r8), dimension(ng31, kc_npress, kc_ntemp) :: k31_h2o, k31_co2, k31_ch4
  real(r8), dimension(ng32, kc_npress, kc_ntemp) :: k32_h2o, k32_co2, k32_ch4
  real(r8), dimension(ng33, kc_npress, kc_ntemp) :: k33_h2o, k33_co2, k33_ch4
  real(r8), dimension(ng34, kc_npress, kc_ntemp) :: k34_h2o, k34_co2, k34_ch4
  real(r8), dimension(ng35, kc_npress, kc_ntemp) :: k35_h2o, k35_co2, k35_ch4
  real(r8), dimension(ng36, kc_npress, kc_ntemp) :: k36_h2o, k36_co2, k36_ch4
  real(r8), dimension(ng37, kc_npress, kc_ntemp) :: k37_h2o, k37_co2, k37_ch4
  real(r8), dimension(ng38, kc_npress, kc_ntemp) :: k38_h2o, k38_co2, k38_ch4
  real(r8), dimension(ng39, kc_npress, kc_ntemp) :: k39_h2o, k39_co2, k39_ch4
  real(r8), dimension(ng40, kc_npress, kc_ntemp) :: k40_h2o, k40_co2, k40_ch4
  real(r8), dimension(ng41, kc_npress, kc_ntemp) :: k41_h2o, k41_co2, k41_ch4
  real(r8), dimension(ng42, kc_npress, kc_ntemp) :: k42_h2o, k42_co2, k42_ch4
  real(r8), dimension(ng43, kc_npress, kc_ntemp) :: k43_h2o, k43_co2, k43_ch4
  real(r8), dimension(ng44, kc_npress, kc_ntemp) :: k44_h2o, k44_co2, k44_ch4
  real(r8), dimension(ng45, kc_npress, kc_ntemp) :: k45_h2o, k45_co2, k45_ch4
  real(r8), dimension(ng46, kc_npress, kc_ntemp) :: k46_h2o, k46_co2, k46_ch4
  real(r8), dimension(ng47, kc_npress, kc_ntemp) :: k47_h2o, k47_co2, k47_ch4
  real(r8), dimension(ng48, kc_npress, kc_ntemp) :: k48_h2o, k48_co2, k48_ch4
  real(r8), dimension(ng49, kc_npress, kc_ntemp) :: k49_h2o, k49_co2, k49_ch4
  real(r8), dimension(ng50, kc_npress, kc_ntemp) :: k50_h2o, k50_co2, k50_ch4
  real(r8), dimension(ng51, kc_npress, kc_ntemp) :: k51_h2o, k51_co2, k51_ch4
  real(r8), dimension(ng52, kc_npress, kc_ntemp) :: k52_h2o, k52_co2, k52_ch4
  real(r8), dimension(ng53, kc_npress, kc_ntemp) :: k53_h2o, k53_co2, k53_ch4
  real(r8), dimension(ng54, kc_npress, kc_ntemp) :: k54_h2o, k54_co2, k54_ch4
  real(r8), dimension(ng55, kc_npress, kc_ntemp) :: k55_h2o, k55_co2, k55_ch4
  real(r8), dimension(ng56, kc_npress, kc_ntemp) :: k56_h2o, k56_co2, k56_ch4
  real(r8), dimension(ng57, kc_npress, kc_ntemp) :: k57_h2o, k57_co2, k57_ch4
  real(r8), dimension(ng58, kc_npress, kc_ntemp) :: k58_h2o, k58_co2, k58_ch4
  real(r8), dimension(ng59, kc_npress, kc_ntemp) :: k59_h2o, k59_co2, k59_ch4
  real(r8), dimension(ng60, kc_npress, kc_ntemp) :: k60_h2o, k60_co2, k60_ch4
  real(r8), dimension(ng61, kc_npress, kc_ntemp) :: k61_h2o, k61_co2, k61_ch4
  real(r8), dimension(ng62, kc_npress, kc_ntemp) :: k62_h2o, k62_co2, k62_ch4
  real(r8), dimension(ng63, kc_npress, kc_ntemp) :: k63_h2o, k63_co2, k63_ch4
  real(r8), dimension(ng64, kc_npress, kc_ntemp) :: k64_h2o, k64_co2, k64_ch4
  real(r8), dimension(ng65, kc_npress, kc_ntemp) :: k65_h2o, k65_co2, k65_ch4
  real(r8), dimension(ng66, kc_npress, kc_ntemp) :: k66_h2o, k66_co2, k66_ch4
  real(r8), dimension(ng67, kc_npress, kc_ntemp) :: k67_h2o, k67_co2, k67_ch4
  real(r8), dimension(ng68, kc_npress, kc_ntemp) :: k68_h2o, k68_co2, k68_ch4

  ! CIA absorption data from HITRAN
  real(r8), dimension(ntot_wavlnrng,kh2oh2o_ntemp) :: kh2oh2o  ! H2O-H2O CIA data [cm-1 amagat-2]
  real(r8), dimension(ntot_wavlnrng,kh2on2_ntemp) :: kh2on2  ! H2O-N2 CIA data [cm-1 amagat-2]

  real(r8), dimension(ntot_wavlnrng,kh2h2_ntemp) :: kh2h2  ! H2-H2 CIA data [cm-1 amagat-2]
  real(r8), dimension(ntot_wavlnrng,kn2h2_ntemp) :: kn2h2  ! H2-N2 CIA data [cm-1 amagat-2]
  real(r8), dimension(ntot_wavlnrng,kn2n2_ntemp) :: kn2n2  ! N2-N2 CIA data [cm-1 amagat-2]

  real(r8), dimension(ntot_wavlnrng,kco2co2_lw_ntemp) :: kco2co2_lw  ! CO2-CO2 CIA lw data [cm-1 amagat-2]
  real(r8), dimension(ntot_wavlnrng,kco2co2_sw_ntemp) :: kco2co2_sw  ! CO2-CO2 CIA sw data [cm-1 amagat-2]
  real(r8), dimension(ntot_wavlnrng,kco2h2_ntemp)  :: kco2h2   ! CO2-H2 CIA data [cm-1 amagat-2]
  real(r8), dimension(ntot_wavlnrng,kco2ch4_ntemp) :: kco2ch4  ! CO2-H2 CIA data [cm-1 amagat-2]

 ! 8 to 16 gauss point mapping
!  integer, dimension(16) :: map8to16gpt
!
!  data map8to16gpt / 1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6, 7, 7, 8, 8 /



  !
  ! Cloud optics parameters  
  !
  integer, parameter :: ncld_grp = 2     ! number of cloud groups, 1=water clouds, 2=ice clouds
  integer, parameter :: nrel = 30        ! number of radii grid points, liquid
  integer, parameter :: nrei = 300        ! number of radii grid points, ice
  !integer, parameter :: nrei_co2 = 20        ! number of radii grid points, CO2 ice    
  integer, parameter :: nrei_co2 = 200        ! number of radii grid points, CO2 ice    

  real(r8), parameter :: cldmin = 1.0d-80
  real(r8), dimension(nrel) :: rel_grid
  real(r8), dimension(nrei) :: rei_grid
  real(r8), dimension(nrei_co2) :: rei_co2_grid

  ! data structures containing cloud optics data from file
  real(r8) :: Qcldliq(nrel, ntot_wavlnrng)
  real(r8) :: Wcldliq(nrel, ntot_wavlnrng)
  real(r8) :: Gcldliq(nrel, ntot_wavlnrng)

  real(r8) :: Qcldice(nrei, ntot_wavlnrng)
  real(r8) :: Wcldice(nrei, ntot_wavlnrng)
  real(r8) :: Gcldice(nrei, ntot_wavlnrng)

  real(r8) :: Qcldice_co2(nrei_co2, ntot_wavlnrng)
  real(r8) :: Wcldice_co2(nrei_co2, ntot_wavlnrng)
  real(r8) :: Gcldice_co2(nrei_co2, ntot_wavlnrng)

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

  data rei_co2_grid / 1., 2., 3., 4., 5., 6., 7., 8., 9., 10., &
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
                    191., 192., 193., 194., 195., 196., 197., 198., 199., 200. /
        

!  data rei_co2_grid / 1.,   5.,   10.,  15.,  20.,  25.,  30.,  35.,  40.,  45., &
!                      50.,  60.,  70.,  80.,  90.,  100., 125., 150., 175., 200. /



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
