module interpolator_module
    implicit none
    integer, parameter :: WP = kind(1.0D0)

    interface
        ! This is the general form of the right-hand side of an ODE
        function rhs(X, t) result( val )
            import :: WP
            real(WP), dimension(:),      intent(in) :: x
            real(WP),                    intent(in) :: t
            real(WP), dimension(size(x))            :: val
        end function
    end interface

    type interpolator_type
        ! This type would in practice store arrays,
        ! of discrete data to be interpolated.
        real(WP) :: stored_data
        procedure(rhs), nopass, pointer :: eval
    contains
        procedure :: init
    endtype

    class(interpolator_type), pointer :: interpolator

contains

    subroutine init( this, stored_data )
        implicit none
        class(interpolator_type), target :: this
        real(WP) :: stored_data
        this % stored_data = stored_data
        this % eval => evaluate
        interpolator => this
    end subroutine

    function evaluate(X, t) result( val )
        implicit none
        real(WP), dimension(:),      intent(in) :: x
        real(WP),                    intent(in) :: t
        real(WP), dimension(size(x))            :: val
        ! This is where interpolation would happen
        val = interpolator % stored_data * x
    end function

end module


program main
    use interpolator_module, only : interpolator_type
    implicit none
    integer, parameter :: WP = kind(1.0D0)
    type(interpolator_type) :: interp
    real(WP), dimension(2) :: X
    real(WP) :: t, h

    ! initialise interpolator with some data
    call interp % init(-0.1_WP)

    X = (/ 2.0_WP, 1.0_WP /)
    t = 0.0_WP
    h = 1.0_WP

    ! Example of calling rk1 with the "type-bound procedure"
    ! which evaluates an interpolator
    call rk4(X, t, h, interp % eval )
    print *, X

    ! Example of calling rk1 with analytical function
    call rk4(X, t, h, f )
    print *, X

    contains

    subroutine rk4(X, t, h, f)
        ! Makes one step with 4th-order Runge-Kutta.
        ! Calculates next position using timestep h.
        implicit none
        real(WP), intent(inout), dimension(:)   :: X
        real(WP), intent(inout)                 :: t
        real(WP), intent(in)                    :: h
        interface
            function f(X, t) result(val)
                import WP
                real(WP), dimension(:),      intent(in) :: x
                real(WP),                    intent(in) :: t
                real(WP), dimension(size(x))            :: val
            end function
        end interface
        ! Local variables
        real(WP), dimension(size(X))            :: k1, k2, k3, k4

        ! Evaluations of f(x, t)
        k1 = f(X,          t)
        k2 = f(X + k1*h/2, t + h/2)
        k3 = f(X + k2*h/2, t + h/2)
        k4 = f(X + k3,     t + h)

        ! Next position
        X  = X + h*(k1 + 2*k2 + 2*k3 + k4)/6
        t  = t + h
    end subroutine

    pure function f(x, t) result(val)
        implicit none
        real(WP), dimension(:),      intent(in) :: x
        real(WP),                    intent(in) :: t
        real(WP), dimension(size(x))            :: val
        val = -0.1_WP*x
    end function

end program
