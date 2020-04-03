program run

use parameters,             only: SP, DP, WP, timesteps_ref, tolerances_ref
use input_module,           only: read_initial_positions
use currentdata_module,     only: get_current
use interpolator_module,    only: interpolator
use experiment_module,      only: experiment_fixed, experiment_variable, experiment_special
use integrator_module,      only: rk4, dp87

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
! Variables for filenames
character(len=256) :: inputfilename
character(len=256) :: currentfilename


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!! Setup of experiment !!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Current data and initial positions
currentfilename = '../data/Nordic-4km.nc'
inputfilename   = '../data/initial_positions_nordic.txt'


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
    ! Run simulations with rk4
!    call experiment_fixed(X0, t0, tmax, timesteps_ref, f, &
!            rk4, 'reference_nordic_rk4_' // trim(suffix(order)) // '.hdf5')
!    ! Run simulations with dp87
!    call experiment_variable(X0, t0, tmax, tolerances_ref, f, &
!            dp87, 'reference_nordic_dp87_' // trim(suffix(order)) // '.hdf5', h0input = 1.0_WP)
    ! Run simulations with dp87 special
    call experiment_special(X0, t0, tmax, tc, tolerances_ref, f, &
            dp87, 'reference_nordic_dp87_special_' // trim(suffix(order)) // '.hdf5', h0input = 1.0_WP)
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
