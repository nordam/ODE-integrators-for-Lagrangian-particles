program run

use parameters,             only: SP, DP, WP
use input_module,           only: read_initial_positions
use currentdata_module,     only: get_current
use interpolator_module,    only: interpolator
use integrator_module,      only: rk1, rk2, rk3, rk4
use integrator_module,      only: bs32, dp54, dp87
use experiment_module,      only: experiment_fixed, experiment_variable, experiment_special

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
tolerances = (/ 1e-3_WP, 1e-4_WP, 1e-5_WP, 1e-6_WP, 1e-7_WP /)



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!! Prepare current data, interpolator, initial positions !!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Get data from NetCDF file
call get_current(currentfilename, u, v, xc, yc, tc)
! Create interpolator from discrete data !!!!
call f%init(xc, yc, tc, u, v, 6)
! Deallocate temporary storage arrays
! (data is now copied and stored in interpolator)
deallocate(u)
deallocate(v)

! Read initial particle positions from file
call read_initial_positions(inputfilename, X0)

! Duration of integration
t0   = tc(6)
tmax = tc(18)



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!! Run experiments for fixed timesteps !!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!outputfilename = 'experiment_rk1_linear.hdf5'
!call experiment_fixed(X0, t0, tmax, timesteps, f, rk1, outputfilename)
!
!outputfilename = 'experiment_rk2_linear.hdf5'
!call experiment_fixed(X0, t0, tmax, timesteps, f, rk2, outputfilename)
!
!outputfilename = 'experiment_rk3_linear.hdf5'
!call experiment_fixed(X0, t0, tmax, timesteps, f, rk3, outputfilename)
!
!outputfilename = 'experiment_rk4_linear.hdf5'
!call experiment_fixed(X0, t0, tmax, timesteps, f, rk4, outputfilename)



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!! Run experiments for variable timesteps !!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

outputfilename = 'experiment_bs32_linear.hdf5'
call experiment_variable(X0, t0, tmax, tolerances, f, bs32, outputfilename, h0input = h0)

!outputfilename = 'experiment_dp54_linear.hdf5'
!call experiment_variable(X0, t0, tmax, tolerances, f, dp54, outputfilename, h0input = h0)
!
!outputfilename = 'experiment_dp87_linear.hdf5'
!call experiment_variable(X0, t0, tmax, tolerances, f, dp87, outputfilename, h0input = h0)
!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!! Run experiments with special-purpose integrators !!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Passing list of timestamps for the currentdata as stoptimes

!outputfilename = 'experiment_bs32_special_linear.hdf5'
!call experiment_special(X0, t0, tmax, tc, tolerances, f, bs32, outputfilename, h0input = h0)
!
!outputfilename = 'experiment_dp54_special_linear.hdf5'
!call experiment_special(X0, t0, tmax, tc, tolerances, f, dp54, outputfilename, h0input = h0)
!
!outputfilename = 'experiment_dp87_special_linear.hdf5'
!call experiment_special(X0, t0, tmax, tc, tolerances, f, dp87, outputfilename, h0input = h0)


end program
