program main
!----------------------------------------------------
! Driver for 1d offline radiative transfer calculations using
! CAM radiative transfer code written by E.Wolf for Archean
! and exoplanetary atmospheres
!
! Stage E1: the input file may carry any number of columns (optional
! 'ncol' dimension; absent = 1). Columns are solved in a loop over
! run_one_column — the same per-column path the C library uses.
!
! Stage E2: the column loop is OpenMP-parallel (columns are independent;
! run_one_column writes no module-scope state, see STAGE_E_AUDIT.md).
! Threads only spawn for multi-column inputs (if(ncol>1)). The large
! solver work arrays are heap-allocated (E2 converted them from
! automatics), so default thread stacks suffice.
!----------------------------------------------------

use io
use exoplanet_mod,        only: do_exo_clouds, do_exo_haze
use initialize_rad_mod_1D
use exo_init_ref,         only: init_ref
use exo_model_specific,   only: init_model_specific
use exo_radiation_mod,    only: init_planck
use exort_column_mod,     only: column_state_t, column_result_t
use exort_column_run,     only: run_one_column
!$ use omp_lib,           only: omp_get_max_threads

implicit none

integer :: icol, ncol
integer(kind=8) :: clk0, clk1, clkrate
real    :: t_input, t_init, t_kernel, t_output

type(column_state_t),  allocatable :: states(:)
type(column_result_t), allocatable :: results(:)

! --- Runtime namelist: read user_nl_exort if present, else compile-time defaults ---
call read_namelist

call system_clock(clk0, clkrate)
call initialize_kcoeff
call initialize_solar
if (do_exo_clouds) call initialize_cldopts
if (do_exo_haze) call initialize_hazeopts
call init_ref
call init_model_specific
call init_planck
call initialize_radbuffer
call system_clock(clk1)
t_init = real(clk1-clk0)/real(clkrate)

call system_clock(clk0)
call input_profile(states)
ncol = size(states)
allocate(results(ncol))
call system_clock(clk1)
t_input = real(clk1-clk0)/real(clkrate)

!$ if (ncol > 1 .and. omp_get_max_threads() > 1) then
!$   write(*,'(a,i0,a,i0,a)') ' OpenMP: solving ', ncol, ' columns on ', &
!$                            omp_get_max_threads(), ' threads'
!$ endif

call system_clock(clk0)
!$omp parallel do schedule(dynamic) if(ncol > 1)
do icol = 1, ncol
  call run_one_column(states(icol), results(icol), icol)
enddo
!$omp end parallel do
call system_clock(clk1)
t_kernel = real(clk1-clk0)/real(clkrate)

! Print Primary Diagnostic outputs
do icol = 1, ncol
  if (ncol > 1) write(*,'(a,i0,a,i0,a)') ' === column ', icol, ' of ', ncol, ' ==='
  call print_diagnostics( results(icol)%sol_toa, &
                          results(icol)%vis_dir, results(icol)%vis_dif, &
                          results(icol)%nir_dir, results(icol)%nir_dif, &
                          results(icol)%sw_dnflux, results(icol)%sw_upflux, &
                          results(icol)%lw_dnflux, results(icol)%lw_upflux )
enddo

call system_clock(clk0)
call output_data( states, results )
call system_clock(clk1)
t_output = real(clk1-clk0)/real(clkrate)

write(*,*) '=== timing (wall seconds) =================='
write(*,'(a,f10.4)') '  initialization : ', t_init
write(*,'(a,f10.4)') '  input          : ', t_input
write(*,'(a,f10.4)') '  aerad_driver   : ', t_kernel
write(*,'(a,f10.4)') '  output         : ', t_output
write(*,'(a,f10.4)') '  total          : ', t_init+t_input+t_kernel+t_output
write(*,*) '============================================'

end program main
