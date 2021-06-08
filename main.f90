! Program: planets
! By: Noah Egger
!
! This program solves the coupled ODEs system governing the motion of 
! two planets (m_1 and m_2) about a central fixed primary mass utilizing
! the 4th order Runge-Kutta method. The masses are confined to move in the 
! xy plane. The program reads a namelist file or default values that specify 
! the intial conditions of position, velocity as well as the mass for all 
! three objects (the primary mass does not include position or velocity as 
! it is fixed and considered the origin) and the total time interval. 
! The program also calculates the total energy and angular momentum of the 
! system for each time step that writes to a new file.
!
!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
program planets

use types
use read_write, only : read_input, write_results
use ode_solver, only : solve_runge_kutta_4
use mechanics, only : calculate_energy, planets_ode, calculate_momentum

implicit none

real(dp) :: work_array(1:3), initial_condition(1:8)
real(dp) :: final_time
integer :: n_steps
real(dp), allocatable :: time(:), solution(:,:), energy(:), ang_momentum(:)
character(len=1024) :: output_file

call read_input(work_array, initial_condition, final_time, n_steps, output_file)

call solve_runge_kutta_4(planets_ode, initial_condition, final_time, work_array, time, solution, n_steps)

call calculate_energy(solution, energy, work_array)

call calculate_momentum(solution, ang_momentum, work_array)

call write_results(time, solution, energy, ang_momentum, output_file)


end program planets