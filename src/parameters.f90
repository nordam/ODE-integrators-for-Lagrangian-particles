module parameters
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
end module

! vim: ai ts=4 sts=4 et sw=4 tw=79 fenc=utf-8
