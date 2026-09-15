
module cloud

implicit none
public

  ! directory
  character(len=256), parameter :: dircld = 'data/cloud/'

  ! Cloud mie data
  character(len=256), parameter :: cldoptsL_file = 'cloudoptics_h2o_liquid_mie_n84.nc'
  character(len=256), parameter :: cldoptsI_file = 'cloudoptics_h2o_ice_mie_n84.nc'
  character(len=256), parameter :: cldoptsICO2_file = 'cloudoptics_co2_ice_mie_n84.nc'

  ! CARMA haze aerosol optics. Mie spheres, regenerated on the full 84-band
  ! grid from the Khare et al. (1984) tholin indices by
  ! tools/makeCARMAOptics.py. The fractal-aggregate table
  ! haze_n84_b40_fractal_interp.nc is also available, but its UV bands
  ! (69-84) are still a provisional extension of band 68 pending a rerun of
  ! the external mean-field fractal solver.
  character(len=256), parameter :: diraer = 'data/aerosol/'
  character(len=256), parameter :: hazeopts_file = 'haze_n84_b40_mie.nc'

end module cloud
