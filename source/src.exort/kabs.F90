
module kabs
!  version eqv, read in single gas k-coefficients from files

  implicit none
  public

  ! correlated-k coefficients for gas absorption
  ! v2 layout: data/kdist/<gas>/ (flat per-gas dirs; legacy n68<gas>/ names gone).
  ! 84-band grid, HITRAN-2024 line list (HITRAN-2020 for O2/O3, current for those).
  character(len=256), parameter :: dirk_h2o  = 'data/kdist/h2o/'
  character(len=256), parameter :: dirk_co2  = 'data/kdist/co2/'
  character(len=256), parameter :: dirk_ch4  = 'data/kdist/ch4/'
  character(len=256), parameter :: dirk_c2h6 = 'data/kdist/c2h6/'
  character(len=256), parameter :: dirk_nh3  = 'data/kdist/nh3/'
  character(len=256), parameter :: dirk_co   = 'data/kdist/co/'
  character(len=256), parameter :: dirk_o3   = 'data/kdist/o3/'
  character(len=256), parameter :: dirk_o2   = 'data/kdist/o2/'

  character(len=256), parameter :: k_h2o_file  = 'n84_8gpt_h2o_hitran24_Nnu1e4_c25_voigt_noplinth_q0_grrtm.nc'
  character(len=256), parameter :: k_co2_file  = 'n84_8gpt_co2_hitran24_Nnu1e4_c500_subL_q1_grrtm.nc'
  character(len=256), parameter :: k_ch4_file  = 'n84_8gpt_ch4_hitran24_Nnu1e4_c25_voigt_q0_grrtm.nc'
  character(len=256), parameter :: k_c2h6_file = 'n84_8gpt_c2h6_hitran24_Nnu1e4_c25_voigt_q0_grrtm.nc'
  character(len=256), parameter :: k_nh3_file  = 'n84_8gpt_nh3_hitran24_Nnu1e4_c25_voigt_q0_grrtm.nc'
  character(len=256), parameter :: k_co_file   = 'n84_8gpt_co_hitran24_Nnu1e4_c25_voigt_q0_grrtm.nc'
  character(len=256), parameter :: k_o2_file   = 'n84_8gpt_o2_hitran20_Nnu1e4_c25_voigt_q0_grrtm.nc'
  character(len=256), parameter :: k_o3_file   = 'n84_8gpt_o3_hitran20_Nnu1e4_c25_voigt_q0_ubremen_SW_xsect_grrtm.nc'


  ! water vapor continuum
  character(len=256), parameter :: dirct = 'data/continuum/'
  !mtckd
  character(len=256), parameter :: kh2o_mtckd_file = 'KH2O_MTCKD3.3_SELF.FRGN_n84_ngauss.nc'

  ! CIAs
  character(len=256), parameter :: dirci = 'data/cia/n84equiv/'
  character(len=256), parameter :: kn2n2cia_file = 'N2-N2_cia_n84.nc'
  character(len=256), parameter :: kn2h2cia_file = 'N2-H2_cia_n84.nc'
  character(len=256), parameter :: kh2h2cia_file = 'H2-H2_cia_n84.nc'
  character(len=256), parameter :: kco2co2cia_sw_file = 'CO2-CO2_cia_sw_n84.nc'
  character(len=256), parameter :: kco2co2cia_lw_file = 'CO2-CO2_cia_lw_n84.nc'
  character(len=256), parameter :: kco2ch4cia_file = 'CO2-CH4_cia_n84.nc'
  character(len=256), parameter :: kco2h2cia_file = 'CO2-H2_cia_n84.nc'
  character(len=256), parameter :: ko2o2cia_file = 'O2-O2_cia_n84.nc'
  character(len=256), parameter :: ko2n2cia_file = 'O2-N2_cia_n84.nc'
  character(len=256), parameter :: ko2co2cia_file = 'O2-CO2_cia_n84.nc'

end module kabs
