program run

use parameters,             only: SP, DP, WP
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
! Timesteps and tolerances
real(WP), dimension(5) :: timesteps
real(WP), dimension(5) :: tolerances
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
character(len=256) :: outputfilename
character(len=256) :: currentfilename


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!! Setup of experiment !!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Current data and initial positions
currentfilename = '../data/Norkyst-800m.nc'
inputfilename   = '../data/initial_positions_norkyst.txt'

! Timesteps and tolerances
timesteps = (/ 3600, 1800, 1200, 900, 600  /)
tolerances = (/ 1e-4_WP, 1e-5_WP, 1e-6_WP, 1e-7_WP, 1e-8_WP /)


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
    call run_all_fixed(X0, t0, tmax, timesteps, f, trim(suffix(order)))
    call run_all_variable(X0, t0, tmax, tolerances, f, trim(suffix(order)), h0)
    call run_all_special(X0, t0, tmax, tc, tolerances, f, trim(suffix(order)), h0)
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
