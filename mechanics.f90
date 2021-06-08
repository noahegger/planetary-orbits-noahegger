!-----------------------------------------------------------------------
! Module: mechanics
!-----------------------------------------------------------------------
!! By: Noah Egger
!!
!! This module contains subroutines responsible for defining the coupled
!! planetary system ODEs as well as for calculating the total energy
!! and momentum for each position throughout the orbit.
!!----------------------------------------------------------------------
!! Included subroutines:
!!
!! calculate_energy
!! calculate_momentum
!!----------------------------------------------------------------------
!! Included functions:
!!
!! planets_ode
!-----------------------------------------------------------------------
module mechanics
use types
implicit none
private
public :: planets_ode, calculate_energy, calculate_momentum

contains


!-----------------------------------------------------------------------
!! function: planets_ode
!-----------------------------------------------------------------------
!! By: Noah Egger
!!
!! This is function will be sent to solve_runge_kutta_4 as an argument "f".
!! Here, we define our planetary motion ODEs. To understand how f(5)-f(8)
!! come about, refer to documentation on the coupled system.
!!
!!----------------------------------------------------------------------
!! Input:
!!
!! r(:)     real        array containing all positions and velocities
!! work(:)  real        1D array containing m_p, m_1, and m_2
!! t        real        duration of orbits (total time)
!!----------------------------------------------------------------------
!! Output:
!!
!! f(:)     real        array containing velocities and accelerations (defines coupled system)
!!
!-----------------------------------------------------------------------
function planets_ode(r, t, work) result(f)
    implicit none
    real(dp), intent(in) :: r(:), t, work(:)
    real(dp), allocatable :: f(:)
    real(dp) :: r_1, r_2, r_12
    integer :: n_variables
    n_variables = size(r)
    allocate(f(1:n_variables))

    ! Define distances between different masses
    r_1 = sqrt(r(1)**2 + r(2)**2)
    r_2 = sqrt(r(3)**2 + r(4)**2)
    r_12 = sqrt((r(1) - r(3))**2 + (r(2) - r(4))**2)

    ! System of differential equations 
    ! r(1) = x_1            
    ! r(2) = y_1            
    ! r(3) = x_2            
    ! r(4) = y_2            
    ! r(5) = vx_1 = f(1)
    ! r(6) = vy_1 = f(2)
    ! r(7) = vx_2 = f(3)
    ! r(8) = vy_2 = f(4)
    !               f(5) = ax_1 = d/dt v_x1     ! Full expression in terms of other variables is found in documentation
    !               f(6) = ay_1 = d/dt vy_1     ! Full expression in terms of other variables is found in documentation
    !               f(7) = ax_2 = d/dt v_x2     ! Full expression in terms of other variables is found in documentation
    !               f(8) = ay_2 = d/dt vy_2     ! Full expression in terms of other variables is found in documentation

    f(1) = r(5)
    f(2) = r(6)
    f(3) = r(7)
    f(4) = r(8)
    f(5) = -work(1)*r(1)/(r_1**3) - work(3)*(r(1)-r(3))/(r_12**3)
    f(6) = -work(1)*r(2)/(r_1**3) - work(3)*(r(2)-r(4))/(r_12**3)
    f(7) = -work(1)*r(3)/(r_2**3) - work(2)*(r(3)-r(1))/(r_12**3)
    f(8) = -work(1)*r(4)/(r_2**3) - work(2)*(r(4)-r(2))/(r_12**3)

end function planets_ode

!-----------------------------------------------------------------------
!! Subroutine: calculate_energy
!-----------------------------------------------------------------------
!! By: Noah Egger
!!
!! This subroutine calculates the total energy of the 3-body gravitational
!! system. The total energy is the sum of the kinetic energies of m_1 and
!! m_2 and the gravitational potential between the primary mass and m_1, m_2
!! and the gravitational potential between m_1 and m_2.
!!----------------------------------------------------------------------
!! Input:
!!
!! solutions(:,:)       real        2D array containing dependent variables (rows) as a function of time(columns)
!! work(:)              real        1D array containing M, m_1, and m_2 (all times G)
!!
!!----------------------------------------------------------------------
!! Output:
!!
!! energy(:)            real        array containing system energy at each time step
!-----------------------------------------------------------------------
subroutine calculate_energy(solutions, energy, work)
    implicit none
    real(dp), intent(in) :: work(:)
    real(dp), allocatable, intent(in) :: solutions(:,:)
    real(dp), allocatable, intent(out) :: energy(:)
    integer, allocatable :: shape_solutions(:)
    integer :: n_steps, i 
    real(dp) :: kinetic, r_1, r_2, r_12, sum_1, sum_2
    real(dp) :: potential_1, potential_2, potential_12, m_1, m_2, m_p

    ! Implicit function shape gives number of rows, columns
    ! Hold those values in 1D array of two elements, shape_solutions
    ! First element will be the number of rows (time steps)
    allocate(shape_solutions(1:2))
    shape_solutions = shape(solutions)
    n_steps = shape_solutions(1)

    ! Make sure energy array is not already written with data
    if(allocated(energy)) deallocate(energy)
    allocate(energy(1:n_steps))

    ! Obtain masses from work array so equations are easier to understand
    m_p = work(1)
    m_1 = work(2)
    m_2 = work(3)

    ! Values are different for each point in space
    ! Thus, must work with arrays ---> Calculate for each i step in time

    do i = 1, n_steps
        ! Calculate intermediate terms in kinetic energy (ignoring 1/2*m factor)
        sum_1 = solutions(i,5)**2 + solutions(i,6)**2 ! v_x1^2 + v_y1^2
        sum_2 = solutions(i,7)**2 + solutions(i,8)**2 ! v_x2^2 + v_y2^2

        ! Calculate distances from solutions array
        r_1 = sqrt(solutions(i,1)**2 + solutions(i,2)**2) ! (x_1^2 + y_1^2)^1/2
        r_2 = sqrt(solutions(i,3)**2 + solutions(i,4)**2) ! (x_2^2 + y_2^2)^1/2
        r_12 = sqrt((solutions(i,1) - solutions(i,3))**2 + (solutions(i,2) - solutions(i,4))**2) ! ((x_1 - x_2)^2 + (y_1 - y_2)^2)^1/2

        ! Calulate total kinetic energy and individual potential energies
        kinetic = 0.5_dp*(m_1*sum_1 + m_2*sum_2)
        potential_1 = -m_1*m_p/r_1
        potential_2 = -m_2*m_p/r_2
        potential_12 = -m_1*m_2/r_12

        ! Calculate total energy
        energy(i) = kinetic + potential_1 + potential_2 + potential_12
    enddo
end subroutine calculate_energy

!-----------------------------------------------------------------------
!! Subroutine: calculate_momentum
!-----------------------------------------------------------------------
!! By: Noah Egger
!!
!! This subroutine calculates the the total angular momentum for the system.
!! We compute the angular momentum of m1 about the primary and m2 about
!! the primary individually by doing a cross product r x p where r and p
!! are vectors. Because the solutions(:,:) array contains the velocties
!! and positions of the masses at various time steps, the calculation is
!! relatively straightforward. 
!!----------------------------------------------------------------------
!! Input:
!!
!! solutions(:,:)       real        2D array containing dependent variables (rows) as a function of time(columns)
!! work(:)              real        1D array containing M, m_1, and m_2 (all times G)
!!
!!----------------------------------------------------------------------
!! Output:
!!
!! ang_momentum(:)      real        1D containing system momentum at each time step
!-----------------------------------------------------------------------
subroutine calculate_momentum(solutions, ang_momentum, work)
    implicit none
    real(dp), intent(in) :: work(:)
    real(dp), allocatable, intent(in) :: solutions(:,:)
    real(dp), allocatable, intent(out) :: ang_momentum(:)
    integer, allocatable :: shape_solutions(:)
    real(dp) :: ang_momentum_1, ang_momentum_2
    integer :: i, n_steps
    allocate(shape_solutions(1:2))

    ! Obtain shape of 2D solutions array. The implicit function 'shape' will
    ! return a 1D array of 2 elements containg the number of rows (1st
    ! element) and the number of columns (2nd element)
    shape_solutions = shape(solutions)

    ! The number of rows (time steps) will be n_steps
    n_steps = shape_solutions(1)

    ! In case ang_momentum is allocated, deallocate
    if (allocated(ang_momentum)) deallocate(ang_momentum)
    allocate(ang_momentum(1:n_steps))

    ! Calculate angular momentum at each time step
    ! L_1 = M_1*(v_1y*x1 - v_1x*y1) = r_1 x p_1
    ! L_2 = M_2*(v_2y*x2 - v_2x*y2) = r_2 x p_2
    do i = 1, n_steps
        ang_momentum_1 = work(2)*(solutions(i,6)*solutions(i,1) - solutions(i,5)*solutions(i,2))
        ang_momentum_2 = work(3)*(solutions(i,8)*solutions(i,3) - solutions(i,7)*solutions(i,4))
        ang_momentum(i) = ang_momentum_1 + ang_momentum_2
    enddo

end subroutine calculate_momentum

   
end module mechanics