
module kabs
!  version n42_h2o

  implicit none
  public

  character(len=256), parameter :: dirk = '/projects/wolfet/models/ExoRT/data/kdist/n68h2o/hitran2012/'

  ! K coefficient file names
  character(len=256), parameter :: k01_file = 'n68_bin01_h2o_hitran12.nc'
  character(len=256), parameter :: k02_file = 'n68_bin02_h2o_hitran12.nc'
  character(len=256), parameter :: k03_file = 'n68_bin03_h2o_hitran12.nc'
  character(len=256), parameter :: k04_file = 'n68_bin04_h2o_hitran12.nc'
  character(len=256), parameter :: k05_file = 'n68_bin05_h2o_hitran12.nc'
  character(len=256), parameter :: k06_file = 'n68_bin06_h2o_hitran12.nc'
  character(len=256), parameter :: k07_file = 'n68_bin07_h2o_hitran12.nc'
  character(len=256), parameter :: k08_file = 'n68_bin08_h2o_hitran12.nc'
  character(len=256), parameter :: k09_file = 'n68_bin09_h2o_hitran12.nc'
  character(len=256), parameter :: k10_file = 'n68_bin10_h2o_hitran12.nc'
  character(len=256), parameter :: k11_file = 'n68_bin11_h2o_hitran12.nc'
  character(len=256), parameter :: k12_file = 'n68_bin12_h2o_hitran12.nc'
  character(len=256), parameter :: k13_file = 'n68_bin13_h2o_hitran12.nc'
  character(len=256), parameter :: k14_file = 'n68_bin14_h2o_hitran12.nc'
  character(len=256), parameter :: k15_file = 'n68_bin15_h2o_hitran12.nc'
  character(len=256), parameter :: k16_file = 'n68_bin16_h2o_hitran12.nc'
  character(len=256), parameter :: k17_file = 'n68_bin17_h2o_hitran12.nc'
  character(len=256), parameter :: k18_file = 'n68_bin18_h2o_hitran12.nc'
  character(len=256), parameter :: k19_file = 'n68_bin19_h2o_hitran12.nc'
  character(len=256), parameter :: k20_file = 'n68_bin20_h2o_hitran12.nc'
  character(len=256), parameter :: k21_file = 'n68_bin21_h2o_hitran12.nc'
  character(len=256), parameter :: k22_file = 'n68_bin22_h2o_hitran12.nc'
  character(len=256), parameter :: k23_file = 'n68_bin23_h2o_hitran12.nc'
  character(len=256), parameter :: k24_file = 'n68_bin24_h2o_hitran12.nc'
  character(len=256), parameter :: k25_file = 'n68_bin25_h2o_hitran12.nc'
  character(len=256), parameter :: k26_file = 'n68_bin26_h2o_hitran12.nc'
  character(len=256), parameter :: k27_file = 'n68_bin27_h2o_hitran12.nc'
  character(len=256), parameter :: k28_file = 'n68_bin28_h2o_hitran12.nc'
  character(len=256), parameter :: k29_file = 'n68_bin29_h2o_hitran12.nc'
  character(len=256), parameter :: k30_file = 'n68_bin30_h2o_hitran12.nc'
  character(len=256), parameter :: k31_file = 'n68_bin31_h2o_hitran12.nc'
  character(len=256), parameter :: k32_file = 'n68_bin32_h2o_hitran12.nc'
  character(len=256), parameter :: k33_file = 'n68_bin33_h2o_hitran12.nc'
  character(len=256), parameter :: k34_file = 'n68_bin34_h2o_hitran12.nc'
  character(len=256), parameter :: k35_file = 'n68_bin35_h2o_hitran12.nc'
  character(len=256), parameter :: k36_file = 'n68_bin36_h2o_hitran12.nc'
  character(len=256), parameter :: k37_file = 'n68_bin37_h2o_hitran12.nc'
  character(len=256), parameter :: k38_file = 'n68_bin38_h2o_hitran12.nc'
  character(len=256), parameter :: k39_file = 'n68_bin39_h2o_hitran12.nc'
  character(len=256), parameter :: k40_file = 'n68_bin40_h2o_hitran12.nc'
  character(len=256), parameter :: k41_file = 'n68_bin41_h2o_hitran12.nc'
  character(len=256), parameter :: k42_file = 'n68_bin42_h2o_hitran12.nc'
  character(len=256), parameter :: k43_file = 'n68_bin43_h2o_hitran12.nc'
  character(len=256), parameter :: k44_file = 'n68_bin44_h2o_hitran12.nc'
  character(len=256), parameter :: k45_file = 'n68_bin45_h2o_hitran12.nc'
  character(len=256), parameter :: k46_file = 'n68_bin46_h2o_hitran12.nc'
  character(len=256), parameter :: k47_file = 'n68_bin47_h2o_hitran12.nc'
  character(len=256), parameter :: k48_file = 'n68_bin48_h2o_hitran12.nc'
  character(len=256), parameter :: k49_file = 'n68_bin49_h2o_hitran12.nc'
  character(len=256), parameter :: k50_file = 'n68_bin50_h2o_hitran12.nc'
  character(len=256), parameter :: k51_file = 'n68_bin51_h2o_hitran12.nc'
  character(len=256), parameter :: k52_file = 'n68_bin52_h2o_hitran12.nc'
  character(len=256), parameter :: k53_file = 'n68_bin53_h2o_hitran12.nc'
  character(len=256), parameter :: k54_file = 'n68_bin54_h2o_hitran12.nc'
  character(len=256), parameter :: k55_file = 'n68_bin55_h2o_hitran12.nc'
  character(len=256), parameter :: k56_file = 'n68_bin56_h2o_hitran12.nc'
  character(len=256), parameter :: k57_file = 'n68_bin57_h2o_hitran12.nc'
  character(len=256), parameter :: k58_file = 'n68_bin58_h2o_hitran12.nc'
  character(len=256), parameter :: k59_file = 'n68_bin59_h2o_hitran12.nc'
  character(len=256), parameter :: k60_file = 'n68_bin60_h2o_hitran12.nc'
  character(len=256), parameter :: k61_file = 'n68_bin61_h2o_hitran12.nc'
  character(len=256), parameter :: k62_file = 'n68_bin62_h2o_hitran12.nc'
  character(len=256), parameter :: k63_file = 'n68_bin63_h2o_hitran12.nc'
  character(len=256), parameter :: k64_file = 'n68_bin64_h2o_hitran12.nc'
  character(len=256), parameter :: k65_file = 'n68_bin65_h2o_hitran12.nc'
  character(len=256), parameter :: k66_file = 'n68_bin66_h2o_hitran12.nc'
  character(len=256), parameter :: k67_file = 'n68_bin67_h2o_hitran12.nc'
  character(len=256), parameter :: k68_file = 'n68_bin68_h2o_hitran12.nc'


  ! K coefficients for continuum files
  logical, parameter :: bps_continuum  = .true.

  character(len=256), parameter :: dirct = '/projects/wolfet/models/ExoRT/data/continuum/'
  !mtckd
  character(len=256), parameter :: kh2oself_mtckd_file = 'KH2OSELF_MTCKD2.5_68bin.nc'
  !bps
  character(len=256), parameter :: kh2oself_bps_file = 'bps_h20_continuum_n68.nc'

  character(len=256), parameter :: dirci = '/projects/wolfet/models/ExoRT/data/cia/'
  character(len=256), parameter :: kn2n2cia_file = 'N2-N2_cia_68bin.nc'
  character(len=256), parameter :: kh2n2cia_file = 'N2-H2_cia_68bin.nc'
  character(len=256), parameter :: kh2h2cia_file = 'H2-H2_cia_68bin.nc'

end module kabs

