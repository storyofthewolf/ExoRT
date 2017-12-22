#include <params.h>
#include <misc.h>

module kabs

implicit none
public

  ! K coefficient file names
  character(len=256), parameter :: k01_file = '/Users/wolfe/desktop/H2O_CorrelatedK/Helios_correlated-K_files_Corrected/Hitran2012_b28/n28_bin01_h2o_hitran12.nc'
  character(len=256), parameter :: k02_file = '/Users/wolfe/desktop/H2O_CorrelatedK/Helios_correlated-K_files_Corrected/Hitran2012_b28/n28_bin02_h2o_hitran12.nc'
  character(len=256), parameter :: k03_file = '/Users/wolfe/desktop/H2O_CorrelatedK/Helios_correlated-K_files_Corrected/Hitran2012_b28/n28_bin03_h2o_hitran12.nc'
  character(len=256), parameter :: k04_file = '/Users/wolfe/desktop/H2O_CorrelatedK/Helios_correlated-K_files_Corrected/Hitran2012_b28/n28_bin04_h2o_hitran12.nc'
  character(len=256), parameter :: k05_file = '/Users/wolfe/desktop/H2O_CorrelatedK/Helios_correlated-K_files_Corrected/Hitran2012_b28/n28_bin05_h2o_hitran12.nc'
  character(len=256), parameter :: k06_file = '/Users/wolfe/desktop/H2O_CorrelatedK/Helios_correlated-K_files_Corrected/Hitran2012_b28/n28_bin06_h2o_hitran12.nc'
  character(len=256), parameter :: k07_file = '/Users/wolfe/desktop/H2O_CorrelatedK/Helios_correlated-K_files_Corrected/Hitran2012_b28/n28_bin07_h2o_hitran12.nc'
  character(len=256), parameter :: k08_file = '/Users/wolfe/desktop/H2O_CorrelatedK/Helios_correlated-K_files_Corrected/Hitran2012_b28/n28_bin08_h2o_hitran12.nc'
  character(len=256), parameter :: k09_file = '/Users/wolfe/desktop/H2O_CorrelatedK/Helios_correlated-K_files_Corrected/Hitran2012_b28/n28_bin09_h2o_hitran12.nc'
  character(len=256), parameter :: k10_file = '/Users/wolfe/desktop/H2O_CorrelatedK/Helios_correlated-K_files_Corrected/Hitran2012_b28/n28_bin10_h2o_hitran12.nc'
  character(len=256), parameter :: k11_file = '/Users/wolfe/desktop/H2O_CorrelatedK/Helios_correlated-K_files_Corrected/Hitran2012_b28/n28_bin11_h2o_hitran12.nc'
  character(len=256), parameter :: k12_file = '/Users/wolfe/desktop/H2O_CorrelatedK/Helios_correlated-K_files_Corrected/Hitran2012_b28/n28_bin12_h2o_hitran12.nc'
  character(len=256), parameter :: k13_file = '/Users/wolfe/desktop/H2O_CorrelatedK/Helios_correlated-K_files_Corrected/Hitran2012_b28/n28_bin13_h2o_hitran12.nc'
  character(len=256), parameter :: k14_file = '/Users/wolfe/desktop/H2O_CorrelatedK/Helios_correlated-K_files_Corrected/Hitran2012_b28/n28_bin14_h2o_hitran12.nc'
  character(len=256), parameter :: k15_file = '/Users/wolfe/desktop/H2O_CorrelatedK/Helios_correlated-K_files_Corrected/Hitran2012_b28/n28_bin15_h2o_hitran12.nc'
  character(len=256), parameter :: k16_file = '/Users/wolfe/desktop/H2O_CorrelatedK/Helios_correlated-K_files_Corrected/Hitran2012_b28/n28_bin16_h2o_hitran12.nc'
  character(len=256), parameter :: k17_file = '/Users/wolfe/desktop/H2O_CorrelatedK/Helios_correlated-K_files_Corrected/Hitran2012_b28/n28_bin17_h2o_hitran12.nc'
  character(len=256), parameter :: k18_file = '/Users/wolfe/desktop/H2O_CorrelatedK/Helios_correlated-K_files_Corrected/Hitran2012_b28/n28_bin18_h2o_hitran12.nc'
  character(len=256), parameter :: k19_file = '/Users/wolfe/desktop/H2O_CorrelatedK/Helios_correlated-K_files_Corrected/Hitran2012_b28/n28_bin19_h2o_hitran12.nc'
  character(len=256), parameter :: k20_file = '/Users/wolfe/desktop/H2O_CorrelatedK/Helios_correlated-K_files_Corrected/Hitran2012_b28/n28_bin20_h2o_hitran12.nc'
  character(len=256), parameter :: k21_file = '/Users/wolfe/desktop/H2O_CorrelatedK/Helios_correlated-K_files_Corrected/Hitran2012_b28/n28_bin21_h2o_hitran12.nc'
  character(len=256), parameter :: k22_file = '/Users/wolfe/desktop/H2O_CorrelatedK/Helios_correlated-K_files_Corrected/Hitran2012_b28/n28_bin22_h2o_hitran12.nc'
  character(len=256), parameter :: k23_file = '/Users/wolfe/desktop/H2O_CorrelatedK/Helios_correlated-K_files_Corrected/Hitran2012_b28/n28_bin23_h2o_hitran12.nc'
  character(len=256), parameter :: k24_file = '/Users/wolfe/desktop/H2O_CorrelatedK/Helios_correlated-K_files_Corrected/Hitran2012_b28/n28_bin24_h2o_hitran12.nc'
  character(len=256), parameter :: k25_file = '/Users/wolfe/desktop/H2O_CorrelatedK/Helios_correlated-K_files_Corrected/Hitran2012_b28/n28_bin25_h2o_hitran12.nc'
  character(len=256), parameter :: k26_file = '/Users/wolfe/desktop/H2O_CorrelatedK/Helios_correlated-K_files_Corrected/Hitran2012_b28/n28_bin26_h2o_hitran12.nc'
  character(len=256), parameter :: k27_file = '/Users/wolfe/desktop/H2O_CorrelatedK/Helios_correlated-K_files_Corrected/Hitran2012_b28/n28_bin27_h2o_hitran12.nc'
  character(len=256), parameter :: k28_file = '/Users/wolfe/desktop/H2O_CorrelatedK/Helios_correlated-K_files_Corrected/Hitran2012_b28/n28_bin28_h2o_hitran12.nc'




  ! K coefficients for continuum files
  logical, parameter :: bps_continuum  = .true.
  !mtckd
  character(len=256), parameter :: kh2oself_mtckd_file = '/Users/wolfe/Desktop/H2O_CorrelatedK/H2O_selfcontinuum/mtckd2.5/KH2OSELF_MTCKD2.5_28bin.nc'
  !bps
  character(len=256), parameter :: kh2oself_bps_file = '/Users/wolfe/Desktop/H2O_CorrelatedK/H2O_selfcontinuum/bps/bps_h20_continuum_n28.nc'
  !
  character(len=256), parameter :: kco2cont_file = '/Users/wolfe/Models/RT/RT_offline/absorption_data/continuum/KCO2CONT.nc'
  character(len=256), parameter :: kn2n2cia_file = '/Users/wolfe/Models/RT/RT_offline/absorption_data/cia/N2-N2_CIA.nc'
  character(len=256), parameter :: kh2n2cia_file = '/Users/wolfe/Models/RT/RT_offline/absorption_data/cia/H2-N2_CIA.nc'
  character(len=256), parameter :: kh2h2cia_file = '/Users/wolfe/Models/RT/RT_offline/absorption_data/cia/H2-H2_CIA.nc'

end module kabs