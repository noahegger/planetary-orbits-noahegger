!-----------------------------------------------------------------------
! Module: ode_solver
!-----------------------------------------------------------------------
!! By: Noah Egger
!!
!! This module houses the 4th order Runge-Kutta subroutine used
!! for solving a system of 2nd order ODEs.
!!
!!----------------------------------------------------------------------
!! Included subroutines:
!!
!! solve_runge_kutta_4
!!
!!----------------------------------------------------------------------
!! Included functions:
!!
!-----------------------------------------------------------------------
module ode_solver
use types
implicit none
private

public :: solve_runge_kutta_4

interface
    function func(r, t, work) result(f)
        use types, only : dp
        implicit none
        real(dp), intent(in) :: work(:), r(:), t
        real(dp), allocatable :: f(:)
    end function func
end interface

contains

!-----------------------------------------------------------------------
!! Subroutine: solve_runge_kutta_4
!-----------------------------------------------------------------------
!! By: Noah Egger
!!
!! The RK-4 subroutine utilizes a weighted 6-point formula to approximate
!! the solution to an ODE for each given time step. The RK-4 method
!! assumes the system we are solving contains initial values for the functions. 
!!
!!----------------------------------------------------------------------
!! Input:
!!
!! r(:)         real        vector array containing generalized dependent variables (positions in this case)
!! work_array   real        array containing primary mass, mass 1, and mass 2 (all times G)
!! f            real        vector function that depends on generalized variables and time
!! n_steps      real        number of sampling n_points
!! 
!!
!!----------------------------------------------------------------------
!! Output:
!!
!! solutions(:,:)   real    2D array containg generalized dependent variables at each time step
!! time_array(:)    real    1D array containing time values
!-----------------------------------------------------------------------
subroutine solve_runge_kutta_4(f, r, t, work, time_count_array, solutions, n_steps)
    implicit none
    procedure(func) :: f
    real(dp), intent(in) :: r(:), t, work(:)
    integer, intent(in) :: n_steps
    real(dp), intent(out), allocatable :: solutions(:,:), time_count_array(:)
    real(dp), allocatable :: k_1(:), k_2(:), k_3(:), k_4(:)
    real(dp), allocatable :: variables(:)
    integer :: n_variables, i
    real(dp) :: time, h

    ! Number of generalized dependent variables
    n_variables = size(r) 

    ! Allocate arrays to hold k values for each given dependent variable
    allocate(k_1(1:n_variables))
    allocate(k_2(1:n_variables))
    allocate(k_3(1:n_variables))
    allocate(k_4(1:n_variables))
    allocate(variables(1:n_variables))
    allocate(solutions(1:n_steps,1:n_variables))
    allocate(time_count_array(1:n_steps))

    ! Fill variables array with the initial conditions for the dependent variables
    variables = r

    ! Set initial time value to 0 and choose step size (note t = total time)
    time = 0
    h = t/n_steps

    do i = 1, n_steps
        k_1 = h*f(variables, time, work)
        k_2 = h*f(variables + 0.5*k_1, time + 0.5*h, work)
        k_3 = h*f(variables + 0.5*k_2, time + 0.5*h, work)
        k_4 = h*f(variables + k_3, time + h, work)

        ! Calculate next iteration for dependent variables
        variables = variables + (1/6._dp)*(k_1 + 2*k_2 + 2*k_3 + k_4)

        ! List dependent variables (rows) at each time step (column)
        solutions(i,:) = variables

        ! Record time value at each step
        time_count_array(i) = time  

        ! Increase time step
        time = time + h
    enddo
end subroutine solve_runge_kutta_4

    
end module ode_solver