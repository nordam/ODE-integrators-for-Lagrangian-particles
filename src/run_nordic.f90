program run

use parameters,             only: SP, DP, WP, timesteps, tolerances, input_folder
use input_module,           only: read_initial_positions
use currentdata_module,     only: get_current
use interpolator_module,    only: interpolator
use experiment_module,      only: run_all_fixed, run_all_variable, run_all_special

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
! Initial timestep for variable-step integrators
real(WP), parameter :: h0 = 600.0_WP
! Variables for filenames
character(len=256) :: currentdata_filename
character(len=256) :: initial_position_filename

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
call get_current(currentdat_filename, u, v, xc, yc, tc)

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
    ! Create a string identifying dataset and order, for the output filename
    output_filename_suffix = trim(dataset_name) // '_' // trim(order_suffix(order))

    ! Run simulations with all four fixed-step integrators
    call run_all_fixed(X0, t0, tmax, timesteps, f, output_filename_suffix)
    ! Run simulations with all three variable-step integrators
    call run_all_variable(X0, t0, tmax, tolerances, f, output_filename_suffix)
    ! Run simulations with all three special-purpose integrators
    call run_all_special(X0, t0, tmax, tc, tolerances, f, output_filename_suffix)

    ! Deallocate interpolator
    call f%destroy()
enddo

! Deallocate current data
deallocate(u)
deallocate(v)

contains
    ! Convenience function to generate output filenames
    function order_suffix(order)
        implicit none
        integer,          intent(in) :: order
        character(len=12)            :: suffix
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
