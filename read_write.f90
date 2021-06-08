!-----------------------------------------------------------------------
! Module: read_write
!-----------------------------------------------------------------------
!! By: Noah Egger
!!
!! This module reads as input the namelist file responsible for initializing
!! the starting positions and velocities as well as the masses, total time,
!! and the number of time steps. In addition, this module writes the coupled ODE
!! results to a file containing the positions, velocities, energies, angular
!! momenta, and time values. 
!! 
!!----------------------------------------------------------------------
!! Included subroutines:
!!
!! read_input 
!! write_results
!!----------------------------------------------------------------------
!! Included functions:
!!
!!
!-----------------------------------------------------------------------
module read_write
use types
implicit none

private
public :: read_input, write_results

contains

!-----------------------------------------------------------------------
!! Subroutine: read_input
!-----------------------------------------------------------------------
!! By: Noah Egger
!!
!! This subroutine reads the provided namelist files containing the
!! masses for the system, as well as the initial conditions of position,
!! velocity, and the numbers of steps and total time value.
!!----------------------------------------------------------------------
!! Output:
!! 
!! work_array(:)        real        array containing primary mass, mass 1, and mass 2
!! initial_condition    real        array containing initial x_1,y_1,x_2,y_2,vx_1,vy_1,vx_2,vy_2
!! final_time           real        total run time for orbits
!! n_steps              real        number of sampling points
!!  
!!
!-----------------------------------------------------------------------
subroutine read_input(work_array, initial_condition, final_time, n_steps, output_file)
    implicit none
    real(dp), intent(out) :: work_array(1:3)
    real(dp), intent(out) :: initial_condition(1:8)
    real(dp), intent(out) :: final_time
    integer, intent(out) :: n_steps
    character(len=1024), intent(out) :: output_file
    character(len=1024) :: namelist_file
    real(dp) :: primary_mass, planet_mass_1, planet_mass_2
    real(dp) :: initial_pos_1(1:2), initial_pos_2(1:2)
    real(dp) :: initial_vel_1(1:2), initial_vel_2(1:2)
    logical :: file_exists
    integer :: n_arguments, unit, ierror

    ! Section names in namelist files, as well as element names
    namelist /masses/ primary_mass, planet_mass_1, planet_mass_2
    namelist /initial_conditions/ initial_pos_1, initial_pos_2, initial_vel_1, initial_vel_2
    namelist /solution_parameters/ final_time, n_steps
    namelist /output/ output_file


    ! Set default values (sort of random). Gives planet-moon system
    primary_mass  = 100000.0
    planet_mass_1 = 1000.0
    planet_mass_2 = 2.0
    initial_pos_1 = [50.0_dp, 50.0_dp]
    initial_pos_2 = [25.0_dp, 75.0_dp]
    initial_vel_1 = [0.0_dp, 15.0_dp]
    initial_vel_2 = [1.0_dp, 30.0_dp]
    final_time = 200.0
    n_steps = 50000
    output_file = 'basic_results.dat'

    ! get namelist file name from command line
    n_arguments = command_argument_count()

    ! read namelists
    if (n_arguments == 1) then
        call get_command_argument(1, namelist_file)
        inquire(file=trim(namelist_file), exist = file_exists)
        if (file_exists) then
            open(newunit = unit, file = namelist_file )
            read(unit, nml = masses, iostat = ierror)
            if(ierror /= 0) then
                print*, 'Error reading masses namelist'
                stop
            endif
            read(unit, nml = initial_conditions, iostat = ierror)
            if(ierror /= 0) then
                print*, 'Error reading initial_conditions namelist'
                stop
            endif
            read(unit, nml = solution_parameters, iostat = ierror)
            if(ierror /= 0) then
                print*, 'Error reading solution_parameters namelist'
                stop
            endif
            read(unit, nml = output, iostat = ierror)
            if(ierror /= 0) then
                print*, 'Error reading output namelist'
                stop
            endif
            close(unit)
        else
            print*, 'Argument, ', trim(namelist_file)
            print*, 'does not exist. Ending program'
            stop
        endif
    else if(n_arguments /= 0) then
        print*, 'Incorrect number of arguments'
        print*, 'The program takes either 0 or 1 arguments'
        print*, 'See documentation in README.md for details'
        stop
    endif

    ! Fill arrays. Either filled from namelist file or default values set.
    work_array = [primary_mass, planet_mass_1, planet_mass_2]
    initial_condition = [initial_pos_1, initial_pos_2, initial_vel_1, initial_vel_2]

end subroutine read_input

!-----------------------------------------------------------------------
!! Subroutine: write_results
!-----------------------------------------------------------------------
!! By: Noah Egger
!!
!! This subroutine writes the solutions for the coupled system of ODEs.
!! By solutions, we mean the positions, velocities, energies, and 
!! angular momenta of the mass 1 and mass 2 at each time step.
!!----------------------------------------------------------------------
!! Input:
!!
!! time(:)          real        1D array containing time values
!! solutions(:,:)   real        2D array containing values of generalized dependent variables at each time step
!! ang_momentum(:)  real        1D array containing total angular momentum of system at each time step
!! energy(:)        real        1D array containg total energy of system at each time step
!! output_file      character   file results are written to
!-----------------------------------------------------------------------
subroutine write_results(time, solution, energy, ang_momentum, output_file)
    implicit none
    real(dp), intent(in) :: time(:), solution(:,:), energy(:), ang_momentum(:)
    character(len=*), intent(in) :: output_file
    integer, allocatable :: shape_solutions(:)
    integer :: unit, i, n_steps
    allocate(shape_solutions(1:2))
    shape_solutions = shape(solution)
    n_steps = shape_solutions(1)

    ! Create file and write to it
    open(newunit = unit, file = output_file)
    write(unit,'(7a25)') 'time', 'x_1', 'y_1', 'x_2', 'y_2', 'E', 'L'
    ! Write parameters (columns) to file at each time value (rows)
    do i = 1, n_steps
        write(unit,*) time(i), solution(i,1), solution(i,2), solution(i,3), solution(i,4), energy(i), ang_momentum(i)
    enddo
    close(unit)
    print*, 'The positions, energies, angular momenta, and'
    print*, 'time values were written to the file named'
    print*, trim(output_file)
end subroutine write_results
end module read_write
