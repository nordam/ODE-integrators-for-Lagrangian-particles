program run

use parameters,             only: SP, DP, WP, timesteps, tolerances, timesteps_ref
use input_module,           only: read_initial_positions
use currentdata_module,     only: get_current
use interpolator_module,    only: interpolator
use experiment_module,      only: run_all_fixed, run_all_variable, run_all_special
use experiment_module,      only: experiment_fixed
use integrator_module,      only: rk4

implicit none
! Arrays for particles
real(wp), dimension(:,:), allocatable :: X0      ! two-component vectors
! Coordinate arrays
real(WP), dimension(:),     allocatable :: xc, yc, tc
! Velocity x and y components
real(WP), dimension(:,:,:), allocatable :: u, v
! Derived type to evaluate interpolated current data
type(interpolator) :: f
! Time
real(wp) :: t0, tmax
! Loop variables
integer :: i, order
! Orders of interpolation
integer, dimension(3), parameter :: orders = (/ 2, 4, 6 /)
! Initial timestep
real(WP), parameter :: h0 = 600.0_WP
! Variables for filenames
character(len=256) :: inputfilename
character(len=256) :: currentfilename


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!! Setup of experiment !!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Current data and initial positions
currentfilename = '../data/Arctic-20km.nc'
inputfilename   = '../data/initial_positions_arctic20.txt'


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!! Prepare current data, interpolator, initial positions !!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Read initial particle positions from file
call read_initial_positions(inputfilename, X0)
! Get data from NetCDF file
call get_current(currentfilename, u, v, xc, yc, tc)
! Duration of integration (taken from time coordinates of data)
t0   = tc(6)
tmax = tc(6 + 72)


!!!!!!!!!!!!!!!!!!!!!!!!!
!!!! Run simulations !!!!
!!!!!!!!!!!!!!!!!!!!!!!!!

do i = 1, size(orders)
    order = orders(i)
    ! Create interpolator from discrete data !!!!
    call f%init(xc, yc, tc, u, v, order)
    ! Run simulations
    call experiment_fixed(X0, t0, tmax, timesteps_ref, f, rk4, 'arctic_ref_' // trim(suffix(order)))
    ! Deallocate interpolator
    call f%destroy()
enddo

deallocate(u)
deallocate(v)

contains
    function suffix(order)
        integer, intent(in) :: order
        character(len=12)   :: suffix
        if (order == 2) then
            suffix = 'linear'
        else if (order == 3) then
            suffix = 'quadractic'
        else if (order == 4) then
            suffix = 'cubic'
        else if (order == 5) then
            suffix = 'quadic?'
        else if (order == 6) then
            suffix = 'quintic'
        else
            print*, 'Unsupported order: ', order
            stop
        endif
    end function

end program
