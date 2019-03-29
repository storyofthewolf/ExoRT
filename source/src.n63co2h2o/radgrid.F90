#include <params.h>
#include <misc.h>

module radgrid

! version.co2_h2ovar
!-----------------------------------------------------------
! Purpose:  Stores basic radiation module data for use
!           by radiation.F90, mcica.F90, rad_interp_mod.F90
!-----------------------------------------------------------
  
  use shr_kind_mod,      only: r8 => shr_kind_r8

  implicit none
  public


  ! Define number of wavelength intervals
  integer, parameter  :: ntot_wavlnrng = 63        ! total # of wv intervals
  integer, parameter :: ngauss_8gpt = 8
  integer, parameter :: ngauss_16gpt = 16
 
  ! <wavenum_edge_>x refers to the wavenumbers (1/wavelength) at the wavelength
  !  interval edges (note that there are 'ntot_wavlnrng+1' edges for each
  !  wavelength group) [cm^-1]:
  real(r8), dimension(ntot_wavlnrng+1) :: wavenum_edge
  data wavenum_edge /  &          ! all wavenumber edges [cm-1]
       10., 100., 160., 220., 280., 330., 380., 440., 495., 545., 617., 667., 720., & 
       800., 875., 940., 1000., 1065., 1108., 1200., 1275., 1350., 1450., 1550., 1650., & 
       1750., 1850., 1950., 2050., 2200., 2397., 2494., 2796., 3087., 3425., 3760., 4030., & 
       4540., 4950., 5370., 5925., 6390., 6990., 7650., 8315., 8850., 9350., 9650., 10400., & 
       11220., 11870., 12750., 13300., 14470., 15000., 15385., 16667., 18182., 20000., &
       22222., 25000., 28571., 33333., 40000. /

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
  !data g_xpos_edge_8gpt / 0.00000, 0.30192, 0.57571, 0.79583, 0.94178, 0.98890, 0.99576, 0.99939 /

  ! Weights for Gaussian quadrature within each wavelength band [none] (there are
  !  'ngauss_pts' of these):
  !data g_weight_8gpt / 0.30192, 0.27379, 0.22012, 0.14595, 0.04712, 0.00686, 0.00363, 0.00061 /

  ! For 16 gauss pts
  !data g_xpos_edge_16gpt / 0.00000, 0.15275, 0.30192, 0.44402, 0.57571, 0.6939, 0.79583, 0.87911, &
  !                        0.94178, 0.98427, 0.9889, 0.99273, 0.99576, 0.99798, 0.99939, 0.99993 /

  data g_weight_16gpt /  8.136830578222341E-2,   0.121394600875259,     0.134425089573857,     0.130900248184922, &
                         0.118165196764835,      0.101187971554512,     8.316617288672230E-2,  6.601992106904840E-2, &
                         5.077923837800617E-2,   3.787915496627484E-2,  2.737772006598929E-2,  1.911193816262080E-2, &
                         1.280532784767848E-2,   8.138921577805653E-3,  4.795479025026687E-3,  2.484713285219481E-3 /


  ! Gauss point gridding
  integer, parameter  :: ng1 = 16
  integer, parameter  :: ng2 = 16
  integer, parameter  :: ng3 = 16
  integer, parameter  :: ng4 = 16
  integer, parameter  :: ng5 = 16
  integer, parameter  :: ng6 = 16
  integer, parameter  :: ng7 = 16
  integer, parameter  :: ng8 = 16
  integer, parameter  :: ng9 = 16
  integer, parameter  :: ng10 = 16
  integer, parameter  :: ng11 = 16
  integer, parameter  :: ng12 = 16
  integer, parameter  :: ng13 = 16
  integer, parameter  :: ng14 = 16
  integer, parameter  :: ng15 = 16
  integer, parameter  :: ng16 = 16
  integer, parameter  :: ng17 = 16
  integer, parameter  :: ng18 = 16
  integer, parameter  :: ng19 = 16
  integer, parameter  :: ng20 = 16
  integer, parameter  :: ng21 = 16
  integer, parameter  :: ng22 = 16
  integer, parameter  :: ng23 = 16
  integer, parameter  :: ng24 = 16
  integer, parameter  :: ng25 = 16
  integer, parameter  :: ng26 = 16
  integer, parameter  :: ng27 = 16
  integer, parameter  :: ng28 = 16
  integer, parameter  :: ng29 = 16
  integer, parameter  :: ng30 = 16
  integer, parameter  :: ng31 = 16
  integer, parameter  :: ng32 = 16
  integer, parameter  :: ng33 = 16
  integer, parameter  :: ng34 = 16
  integer, parameter  :: ng35 = 16
  integer, parameter  :: ng36 = 16
  integer, parameter  :: ng37 = 16
  integer, parameter  :: ng38 = 16
  integer, parameter  :: ng39 = 16
  integer, parameter  :: ng40 = 16
  integer, parameter  :: ng41 = 16
  integer, parameter  :: ng42 = 16
  integer, parameter  :: ng43 = 16
  integer, parameter  :: ng44 = 16
  integer, parameter  :: ng45 = 16
  integer, parameter  :: ng46 = 16
  integer, parameter  :: ng47 = 16
  integer, parameter  :: ng48 = 16
  integer, parameter  :: ng49 = 16
  integer, parameter  :: ng50 = 16
  integer, parameter  :: ng51 = 16
  integer, parameter  :: ng52 = 16
  integer, parameter  :: ng53 = 16
  integer, parameter  :: ng54 = 16
  integer, parameter  :: ng55 = 16
  integer, parameter  :: ng56 = 16
  integer, parameter  :: ng57 = 16
  integer, parameter  :: ng58 = 16
  integer, parameter  :: ng59 = 16
  integer, parameter  :: ng60 = 16
  integer, parameter  :: ng61 = 16
  integer, parameter  :: ng62 = 16
  integer, parameter  :: ng63 = 16


  ! Define number of gauss points in each spectral interval
  integer, dimension(ntot_wavlnrng) :: ngauss_pts           ! # of Gauss quad pts per interval
  data ngauss_pts / ng1,  ng2,  ng3,  ng4,  ng5,  ng6,  ng7,  ng8,  ng9,  ng10, &
                    ng11, ng12, ng13, ng14, ng15, ng16, ng17, ng18, ng19, ng20, &
                    ng21, ng22, ng23, ng24, ng25, ng26, ng27, ng28, ng29, ng30, & 
                    ng31, ng32, ng33, ng34, ng35, ng36, ng37, ng38, ng39, ng40, &
                    ng41, ng42, ng43, ng44, ng45, ng46, ng47, ng48, ng49, ng50, &
                    ng51, ng52, ng53, ng54, ng55, ng56, ng57, ng58, ng59, ng60, & 
                    ng61, ng62, ng63 /


  real(r8),dimension(ntot_wavlnrng) :: solarflux
  real(r8) :: S0


  ! Dimensions for spectral interval gauss points
  integer, parameter  :: ntot_gpt &                 ! total # of probability intervals
                         = ng1 + ng2 + ng3 + ng4 + ng5 + ng6 + ng7 + ng8  &
                         + ng9 + ng10 + ng11 + ng12 + ng13 + ng14 + ng15  &
                         + ng16 + ng17 + ng18 + ng19 + ng20 + ng21 + ng22 &
                         + ng23 + ng24 + ng25 + ng26 + ng27 + ng28 + ng29 &
                         + ng30 + ng31 + ng32 + ng33 + ng34 + ng35 + ng36 &
                         + ng37 + ng38 + ng39 + ng40 + ng41 + ng42 + ng43 &
                         + ng44 + ng45 + ng46 + ng47 + ng48 + ng49 + ng50 &
                         + ng51 + ng52 + ng53 + ng54 + ng55 + ng56 + ng57 &
                         + ng58 + ng59 + ng60 + ng61 + ng62 + ng63 



  ! Dimensions of current k-coefficient datasets
  integer, parameter  :: kc_npress = 9          ! # of reference pressure
  integer, parameter  :: kc_ntemp = 7           ! # of reference temps
  integer, parameter  :: kc_nvar = 7            ! # of h2o volume mixing ratios

  ! Pressure, temperature, and species weight grids
  real(r8), dimension(kc_npress) :: log10pgrid
  real(r8), dimension(kc_npress) :: pgrid
  real(r8), dimension(kc_ntemp)  :: tgrid
  real(r8), dimension(kc_nvar)   :: wgrid
  real(r8), dimension(kc_nvar)   :: log10wgrid

  ! log 10 pressure grid whole atmosphere
  data log10pgrid / -3.0, -2.0, -1.0, 0.0, 1.0, 2.0, 3.0, 4.0, 5.0 /

  ! pressure grid [mb] whole atmopshere
  data pgrid / 0.001,  0.01,  0.1,  1.0,  10.0,  100.0,   1000.0,   10000.0,   100000.0 /

  ! temperature grids [K]
  data tgrid / 100., 150., 200., 250., 300., 350., 400. / 

  ! species weighting grid  
  data wgrid / 1.0e-7, 1.0e-6, 1.0e-5, 1.0e-4, 1.0e-3, 1.0e-2, 1.0e-1 /
  data log10wgrid /-7.0, -6.0, -5.0, -4.0, -3.0, -2.0, -1.0 /

  !! THIS NEEDS TO BE UPDATED 
!!  ! Water Vapor Self Continuum parameters  
  integer, parameter  :: ks_ntemp = 8       ! # of reference temperatures
!  integer, parameter  :: ks_ntemp = 36       ! # of reference temperatures
  integer, parameter  :: ngH2O = 1         ! # of gauss points per spectral interval in
  integer, parameter  :: ks_npress = 1
  integer, parameter  :: ks_nweight = 1

!  ! T, P, W grid for H2O continuum
  real(r8), dimension(ks_npress) :: pgrid_self
  real(r8), dimension(ks_npress) :: log10pgrid_self
  real(r8), dimension(ks_ntemp) :: tgrid_self
  real(r8), dimension(ks_nweight) :: wgrid_self

  ! temperature grid [K], water vapor self continuum
  data tgrid_self / 150, 200, 250, 300, 350, 400, 450, 500 /
!!  data tgrid_self / 150, 160, 170, 180, 190, 200, 210, 220, 230, 240, 250, &
!!                    260, 270, 280, 290, 300, 310, 320, 330, 340, 350, 360, &
!!                    370, 380, 390, 400, 410, 420, 430, 440, 450, 460, 470, &
!!                    480, 490, 500  /

  ! h2o amount grid [vmr]
!  data wgrid_self / 1.0 /
!
!  ! N2-N2 CIA temperature
!  integer, parameter  :: kn2n2_ntemp = 10       ! # of reference temperatures
!  real(r8), dimension(kn2n2_ntemp) :: tgrid_n2n2
!  data tgrid_n2n2 / 40.0, 51.7, 66.7, 86.2, 111.3, 143.8, 185.7, 239.8, 309.7, 400.0 /
!
!  ! H2-N2 CIA temperature
!  integer, parameter  :: kh2n2_ntemp = 10       ! # of reference temperatures
!  real(r8), dimension(kh2n2_ntemp) :: tgrid_h2n2
!  data tgrid_h2n2 / 40.0, 51.7, 66.7, 86.2, 111.3, 143.8, 185.7, 239.8, 309.7, 400.0 /
!
!  ! H2-H2 CIA temperature
!  integer, parameter  :: kh2h2_ntemp = 7       ! # of reference temperatures
!  real(r8), dimension(kh2h2_ntemp) :: tgrid_h2h2
!  data tgrid_h2h2 / 60.0, 100.0, 150.0, 200.0, 250.0, 300.0, 350.0 /
!

  !
  ! K coefficient variables
  ! 
  real(r8), dimension(ng1, kc_nvar, kc_npress, kc_ntemp) :: k01_data
  real(r8), dimension(ng2, kc_nvar, kc_npress, kc_ntemp) :: k02_data
  real(r8), dimension(ng3, kc_nvar, kc_npress, kc_ntemp) :: k03_data
  real(r8), dimension(ng4, kc_nvar, kc_npress, kc_ntemp) :: k04_data
  real(r8), dimension(ng5, kc_nvar, kc_npress, kc_ntemp) :: k05_data
  real(r8), dimension(ng6, kc_nvar, kc_npress, kc_ntemp) :: k06_data
  real(r8), dimension(ng7, kc_nvar, kc_npress, kc_ntemp) :: k07_data
  real(r8), dimension(ng8, kc_nvar, kc_npress, kc_ntemp) :: k08_data
  real(r8), dimension(ng9, kc_nvar, kc_npress, kc_ntemp) :: k09_data
  real(r8), dimension(ng10, kc_nvar, kc_npress, kc_ntemp) :: k10_data
  real(r8), dimension(ng11, kc_nvar, kc_npress, kc_ntemp) :: k11_data
  real(r8), dimension(ng12, kc_nvar, kc_npress, kc_ntemp) :: k12_data
  real(r8), dimension(ng13, kc_nvar, kc_npress, kc_ntemp) :: k13_data
  real(r8), dimension(ng14, kc_nvar, kc_npress, kc_ntemp) :: k14_data
  real(r8), dimension(ng15, kc_nvar, kc_npress, kc_ntemp) :: k15_data
  real(r8), dimension(ng16, kc_nvar, kc_npress, kc_ntemp) :: k16_data
  real(r8), dimension(ng17, kc_nvar, kc_npress, kc_ntemp) :: k17_data
  real(r8), dimension(ng18, kc_nvar, kc_npress, kc_ntemp) :: k18_data
  real(r8), dimension(ng19, kc_nvar, kc_npress, kc_ntemp) :: k19_data
  real(r8), dimension(ng20, kc_nvar, kc_npress, kc_ntemp) :: k20_data
  real(r8), dimension(ng21, kc_nvar, kc_npress, kc_ntemp) :: k21_data
  real(r8), dimension(ng22, kc_nvar, kc_npress, kc_ntemp) :: k22_data
  real(r8), dimension(ng23, kc_nvar, kc_npress, kc_ntemp) :: k23_data
  real(r8), dimension(ng24, kc_nvar, kc_npress, kc_ntemp) :: k24_data
  real(r8), dimension(ng25, kc_nvar, kc_npress, kc_ntemp) :: k25_data
  real(r8), dimension(ng26, kc_nvar, kc_npress, kc_ntemp) :: k26_data
  real(r8), dimension(ng27, kc_nvar, kc_npress, kc_ntemp) :: k27_data
  real(r8), dimension(ng28, kc_nvar, kc_npress, kc_ntemp) :: k28_data
  real(r8), dimension(ng29, kc_nvar, kc_npress, kc_ntemp) :: k29_data
  real(r8), dimension(ng30, kc_nvar, kc_npress, kc_ntemp) :: k30_data
  real(r8), dimension(ng31, kc_nvar, kc_npress, kc_ntemp) :: k31_data
  real(r8), dimension(ng33, kc_nvar, kc_npress, kc_ntemp) :: k32_data
  real(r8), dimension(ng33, kc_nvar, kc_npress, kc_ntemp) :: k33_data
  real(r8), dimension(ng34, kc_nvar, kc_npress, kc_ntemp) :: k34_data
  real(r8), dimension(ng35, kc_nvar, kc_npress, kc_ntemp) :: k35_data
  real(r8), dimension(ng36, kc_nvar, kc_npress, kc_ntemp) :: k36_data
  real(r8), dimension(ng37, kc_nvar, kc_npress, kc_ntemp) :: k37_data
  real(r8), dimension(ng38, kc_nvar, kc_npress, kc_ntemp) :: k38_data
  real(r8), dimension(ng39, kc_nvar, kc_npress, kc_ntemp) :: k39_data
  real(r8), dimension(ng40, kc_nvar, kc_npress, kc_ntemp) :: k40_data
  real(r8), dimension(ng41, kc_nvar, kc_npress, kc_ntemp) :: k41_data
  real(r8), dimension(ng42, kc_nvar, kc_npress, kc_ntemp) :: k42_data
  real(r8), dimension(ng43, kc_nvar, kc_npress, kc_ntemp) :: k43_data
  real(r8), dimension(ng44, kc_nvar, kc_npress, kc_ntemp) :: k44_data
  real(r8), dimension(ng45, kc_nvar, kc_npress, kc_ntemp) :: k45_data
  real(r8), dimension(ng46, kc_nvar, kc_npress, kc_ntemp) :: k46_data
  real(r8), dimension(ng47, kc_nvar, kc_npress, kc_ntemp) :: k47_data
  real(r8), dimension(ng48, kc_nvar, kc_npress, kc_ntemp) :: k48_data
  real(r8), dimension(ng49, kc_nvar, kc_npress, kc_ntemp) :: k49_data
  real(r8), dimension(ng50, kc_nvar, kc_npress, kc_ntemp) :: k50_data
  real(r8), dimension(ng51, kc_nvar, kc_npress, kc_ntemp) :: k51_data
  real(r8), dimension(ng52, kc_nvar, kc_npress, kc_ntemp) :: k52_data
  real(r8), dimension(ng53, kc_nvar, kc_npress, kc_ntemp) :: k53_data
  real(r8), dimension(ng54, kc_nvar, kc_npress, kc_ntemp) :: k54_data
  real(r8), dimension(ng55, kc_nvar, kc_npress, kc_ntemp) :: k55_data
  real(r8), dimension(ng56, kc_nvar, kc_npress, kc_ntemp) :: k56_data
  real(r8), dimension(ng57, kc_nvar, kc_npress, kc_ntemp) :: k57_data
  real(r8), dimension(ng58, kc_nvar, kc_npress, kc_ntemp) :: k58_data
  real(r8), dimension(ng59, kc_nvar, kc_npress, kc_ntemp) :: k59_data
  real(r8), dimension(ng60, kc_nvar, kc_npress, kc_ntemp) :: k60_data
  real(r8), dimension(ng61, kc_nvar, kc_npress, kc_ntemp) :: k61_data
  real(r8), dimension(ng62, kc_nvar, kc_npress, kc_ntemp) :: k62_data
  real(r8), dimension(ng63, kc_nvar, kc_npress, kc_ntemp) :: k63_data

 
  
  !! Continuum Data !!
  !
  ! water vapor self continuum data
  !
  real(r8), dimension(ntot_wavlnrng, ks_ntemp) :: kh2oself
  real(r8), dimension(ntot_wavlnrng, ks_ntemp) :: kh2ofrgn

  !
  ! CO2 continuum information, self and foreign
  !
  !integer, parameter   :: ngCO2 = 8                           ! # of gauss points per spectral interval
  !real(r8), parameter :: pref_co2 = 1013.0
  !real(r8), parameter :: vref_co2 = 0.1 ! reference CO2 for continuum calculation
  !real(r8), dimension(ngCO2,ntot_wavlnrng) :: kco2cont_8gpt  ! CO2 continuum data from file
  !real(r8), dimension(ntot_gpt) :: kco2cont  ! CO2 continuum data, reduced to gauss pt grid
  !!!!!real(r8), dimension(ntot_wavlnrng, ks_ntemp) :: kco2cont

  ! 8 to 16 gauss point mapping
  !integer, dimension(16) :: map8to16gpt

  !data map8to16gpt / 1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6, 7, 7, 8, 8 /

 ! CIA absorption data from HITRAN
 ! real(r8), dimension(ntot_wavlnrng,kh2h2_ntemp) :: kh2h2  ! H2-H2 CIA data [cm-1 amagat-2]
 ! real(r8), dimension(ntot_wavlnrng,kh2n2_ntemp) :: kh2n2  ! H2-N2 CIA data [cm-1 amagat-2]
 ! real(r8), dimension(ntot_wavlnrng,kn2n2_ntemp) :: kn2n2  ! N2-N2 CIA data [cm-1 amagat-2]



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



  ! For MCICA calculation, must match ng* specified above
  integer, dimension(ntot_gpt) :: ngb 
  data ngb / 1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1, &
             2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2, &
             3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3, &
             4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4, &
             5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5, &
             6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6, &
             7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7, &
             8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8, &
             9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9, &
             10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10, &
             11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11, &
             12,12,12,12,12,12,12,12,12,12,12,12,12,12,12,12, &
             13,13,13,13,13,13,13,13,13,13,13,13,13,13,13,13, &
             14,14,14,14,14,14,14,14,14,14,14,14,14,14,14,14, &
             15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15, &
             16,16,16,16,16,16,16,16,16,16,16,16,16,16,16,16, &
             17,17,17,17,17,17,17,17,17,17,17,17,17,17,17,17, &
             18,18,18,18,18,18,18,18,18,18,18,18,18,18,18,18, &
             19,19,19,19,19,19,19,19,19,19,19,19,19,19,19,19, &
             20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20, &
             21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21, &
             22,22,22,22,22,22,22,22,22,22,22,22,22,22,22,22, &
             23,23,23,23,23,23,23,23,23,23,23,23,23,23,23,23, &
             24,24,24,24,24,24,24,24,24,24,24,24,24,24,24,24, &
             25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25, &
             26,26,26,26,26,26,26,26,26,26,26,26,26,26,26,26, &
             27,27,27,27,27,27,27,27,27,27,27,27,27,27,27,27, &
             28,28,28,28,28,28,28,28,28,28,28,28,28,28,28,28, &
             29,29,29,29,29,29,29,29,29,29,29,29,29,29,29,29, &
             30,30,30,30,30,30,30,30,30,30,30,30,30,30,30,30, &
             31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31, &
             32,32,32,32,32,32,32,32,32,32,32,32,32,32,32,32, &
             33,33,33,33,33,33,33,33,33,33,33,33,33,33,33,33, &
             34,34,34,34,34,34,34,34,34,34,34,34,34,34,34,34, &
             35,35,35,35,35,35,35,35,35,35,35,35,35,35,35,35, &
             36,36,36,36,36,36,36,36,36,36,36,36,36,36,36,36, &
             37,37,37,37,37,37,37,37,37,37,37,37,37,37,37,37, &
             38,38,38,38,38,38,38,38,38,38,38,38,38,38,38,38, &
             39,39,39,39,39,39,39,39,39,39,39,39,39,39,39,39, &
             40,40,40,40,40,40,40,40,40,40,40,40,40,40,40,40, &
             41,41,41,41,41,41,41,41,41,41,41,41,41,41,41,41, &
             42,42,42,42,42,42,42,42,42,42,42,42,42,42,42,42, &
             43,43,43,43,43,43,43,43,43,43,43,43,43,43,43,43, &
             44,44,44,44,44,44,44,44,44,44,44,44,44,44,44,44, &
             45,45,45,45,45,45,45,45,45,45,45,45,45,45,45,45, &
             46,46,46,46,46,46,46,46,46,46,46,46,46,46,46,46, &
             47,47,47,47,47,47,47,47,47,47,47,47,47,47,47,47, &
             48,48,48,48,48,48,48,48,48,48,48,48,48,48,48,48, &
             49,49,49,49,49,49,49,49,49,49,49,49,49,49,49,49, &
             50,50,50,50,50,50,50,50,50,50,50,50,50,50,50,40, &
             51,51,51,51,51,51,51,51,51,51,51,51,51,51,51,51, &
             52,52,52,52,52,52,52,52,52,52,52,52,52,52,52,52, &
             53,53,53,53,53,53,53,53,53,53,53,53,53,53,53,53, &
             54,54,54,54,54,54,54,54,54,54,54,54,54,54,54,54, &
             55,55,55,55,55,55,55,55,55,55,55,55,55,55,55,55, &
             56,56,56,56,56,56,56,56,56,56,56,56,56,56,56,56, &
             57,57,57,57,57,57,57,57,57,57,57,57,57,57,57,57, &
             58,58,58,58,58,58,58,58,58,58,58,58,58,58,58,58, &
             59,59,59,59,59,59,59,59,59,59,59,59,59,59,59,59, &
             60,60,60,60,60,60,60,60,60,60,60,60,60,60,60,60, &
             61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61, &
             62,62,62,62,62,62,62,62,62,62,62,62,62,62,62,62, &
             63,63,63,63,63,63,63,63,63,63,63,63,63,63,63,63 /


  !! top level indices for longwave and cloud calculations
  integer, public :: camtop  ! top level to solve using CAM CKD RT.
  integer, public :: ntopcld ! top level to solve for cloud overlap
  integer, public :: nlevsRT ! number of levels to calculate RT



end module radgrid