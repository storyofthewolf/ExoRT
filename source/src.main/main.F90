program main
!----------------------------------------------------
! Driver for 1d offline radiative transfer calculations using
! CAM radiative transfer code written by E.Wolf for Archean
! and exoplanetary atmospheres
!
! Stage E1: the input file may carry any number of columns (optional
! 'ncol' dimension; absent = 1). Columns are solved in a serial loop
! over run_one_column — the same per-column path the C library uses.
! The Stage E2 OpenMP work parallelizes this loop.
!----------------------------------------------------

use io
use exoplanet_mod,        only: do_exo_clouds, do_exo_haze
use initialize_rad_mod_1D
use exo_init_ref,         only: init_ref
use exo_model_specific,   only: init_model_specific
use exo_radiation_mod,    only: init_planck
use exort_column_mod,     only: column_state_t, column_result_t
use exort_column_run,     only: run_one_column

implicit none

integer :: icol, ncol
real    :: t0, t1, t_input, t_init, t_kernel, t_output

type(column_state_t),  allocatable :: states(:)
type(column_result_t), allocatable :: results(:)

! --- Runtime namelist: read user_nl_exort if present, else compile-time defaults ---
call read_namelist

call cpu_time(t0)
call initialize_kcoeff
call initialize_solar
if (do_exo_clouds) call initialize_cldopts
if (do_exo_haze) call initialize_hazeopts
call init_ref
call init_model_specific
call init_planck
call initialize_radbuffer
call cpu_time(t1)
t_init = t1 - t0

call cpu_time(t0)
call input_profile(states)
ncol = size(states)
allocate(results(ncol))
call cpu_time(t1)
t_input = t1 - t0

call cpu_time(t0)
do icol = 1, ncol
  call run_one_column(states(icol), results(icol))
enddo
call cpu_time(t1)
t_kernel = t1 - t0

! Print Primary Diagnostic outputs
do icol = 1, ncol
  if (ncol > 1) write(*,'(a,i0,a,i0,a)') ' === column ', icol, ' of ', ncol, ' ==='
  call print_diagnostics( results(icol)%sol_toa, &
                          results(icol)%vis_dir, results(icol)%vis_dif, &
                          results(icol)%nir_dir, results(icol)%nir_dif, &
                          results(icol)%sw_dnflux, results(icol)%sw_upflux, &
                          results(icol)%lw_dnflux, results(icol)%lw_upflux )
enddo

call cpu_time(t0)
call output_data( states, results )
call cpu_time(t1)
t_output = t1 - t0

write(*,*) '=== timing (cpu seconds) ==================='
write(*,'(a,f10.4)') '  initialization : ', t_init
write(*,'(a,f10.4)') '  input          : ', t_input
write(*,'(a,f10.4)') '  aerad_driver   : ', t_kernel
write(*,'(a,f10.4)') '  output         : ', t_output
write(*,'(a,f10.4)') '  total          : ', t_init+t_input+t_kernel+t_output
write(*,*) '============================================'

end program main
