
program run

use parameters,             only: SP, DP, WP, timesteps_ref, tolerances_ref, input_folder
use input_module,           only: read_initial_positions
use currentdata_module,     only: get_current
use interpolator_module,    only: interpolator
use integrator_module,      only: rk4, dp87
use experiment_module,      only: experiment_fixed, experiment_variable, experiment_special
use experiment_module,      only: get_output_filename

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
! Note that order = degree + 1, so e.g. cubic is order 4
integer, dimension(3), parameter :: orders = (/ 2, 4, 6 /)
! Variables for filenames
character(len=256) :: currentdata_filename
character(len=256) :: initial_position_filename
character(len=256) :: output_filename

! Name identifying the dataset
character(len=16), parameter :: dataset_name = 'nordic4'


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!! Prepare current data, interpolator, initial positions !!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Current data and initial positions
currentdata_filename = trim(input_folder) // trim(dataset_name) // '.nc'
initial_position_filename   = trim(input_folder) // 'initial_positions_' // trim(dataset_name) // '.txt'

! Read initial particle positions from file
call read_initial_positions(initial_position_filename, X0)
! Get data from NetCDF file
call get_current(currentdata_filename, u, v, xc, yc, tc)

! Duration of integration (taken from time coordinates of data)
t0   = tc(6)
tmax = tc(6 + 72)


!!!!!!!!!!!!!!!!!!!!!!!!!
!!!! Run simulations !!!!
!!!!!!!!!!!!!!!!!!!!!!!!!

do i = 1, size(orders)
    order = orders(i)
    ! Create interpolator of desired order from discrete data
    call f%init(xc, yc, tc, u, v, order)

    ! Run simulations with rk4
    output_filename = get_output_filename('reference', 'rk4', dataset_name, order)
    call experiment_fixed(X0, t0, tmax, timesteps_ref, f, rk4, output_filename)
    ! Run simulations with dp87
    output_filename = get_output_filename('reference', 'dp87', dataset_name, order)
    call experiment_variable(X0, t0, tmax, tolerances_ref, f, dp87, output_filename, h0input = 1.0_WP)
    ! Run simulations with dp87 special
    output_filename = get_output_filename('reference', 'dp87s', dataset_name, order)
    call experiment_special(X0, t0, tmax, tc, tolerances_ref, f, dp87, output_filename, h0input = 1.0_WP)
    ! Deallocate interpolator
    call f%destroy()
enddo

deallocate(u)
deallocate(v)

end program
