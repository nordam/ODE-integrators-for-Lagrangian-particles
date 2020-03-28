module parameters
    implicit none
    ! Precision parameters
    integer,  parameter :: SP = kind(1.0E0)
    integer,  parameter :: DP = kind(1.0D0)
    integer,  parameter :: WP = DP
    integer,  parameter :: WI = kind(1)

    ! Mathematical constants
    real(DP), parameter :: PI = 3.1415926535897932384626433

    ! Parameters controlling the stepsize adjustment
    real(WP), parameter    :: fac    = 0.9_wp
    real(WP), parameter    :: maxfac = 3.0_wp

    ! Timesteps and tolerances
    real(WP), dimension(3),  parameter :: timesteps_ref = (/ 30, 10, 5  /)
    real(WP), dimension(10), parameter :: timesteps = (/ &
        3600, 1800, 1200, 900, 600, 450, 300, 180, 120, 60  /)
    real(WP), dimension(10), parameter :: tolerances = (/ &
        1e-4_WP, 1e-5_WP, 1e-6_WP, 1e-7_WP, 1e-8_WP, &
        1e-9_WP, 1e-10_WP, 1e-11_WP, 1e-12_WP, 1e-13_WP /)
end module

! vim: ai ts=4 sts=4 et sw=4 tw=79 fenc=utf-8
