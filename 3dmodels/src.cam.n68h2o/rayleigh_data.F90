
module rayleigh_data

  use shr_kind_mod,      only: r8 => shr_kind_r8
  !
  ! Rayleigh scattering parameters (Vardavas and Carter, 1984)
  !

  ! depolarization factors
  real(r8), parameter :: delN2      = 0.0305  ! 
  real(r8), parameter :: delO2      = 0.054   ! 
  real(r8), parameter :: delCO2     = 0.0805  ! 
  real(r8), parameter :: delH2O     = 0.17    ! Marshall & Smith, (1990)
  real(r8), parameter :: delCH4     = 2.0e-5  ! 
  real(r8), parameter :: delCO      = 0.0048  ! Sneep & Ubachs (2005), doi:10.1016/j.jqsrt.2004.07.025
  real(r8), parameter :: delNH3     = 0.0922  ! Alms, Burnham & Flygare (1975), doi:10.1063/1.431821
  
  ! From Allen (1976) Astrophysical Quantities pg. 92
  ! Polarizability
  real(r8), parameter :: raylA_N2   = 29.06
  real(r8), parameter :: raylB_N2   = 7.70
  real(r8), parameter :: raylA_O2   = 26.63
  real(r8), parameter :: raylB_O2   = 5.07
  real(r8), parameter :: raylA_CO2  = 43.90
  real(r8), parameter :: raylB_CO2  = 6.40
  real(r8), parameter :: raylA_CO   = 32.70   ! Cuthbertson & Cuthbertson (1920), doi:10.1098/rspa.1920.0020
  real(r8), parameter :: raylB_CO   = 8.81
  real(r8), parameter :: raylA_NH3  = 37.0    ! Cuthbertson & Cuthbertson (1914), doi:10.1098/rsta.1914.0001
  real(r8), parameter :: raylB_NH3  = 12.0     

end module rayleigh_data

