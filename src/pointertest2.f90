module mymod
    implicit none
    integer, parameter :: WP = kind(1.0D0)

    interface
        function integrand_interface(X, t) result( val )
            import :: WP
            real(WP), dimension(:)       :: X
            real(WP), dimension(size(X)) :: val
            real(WP)                     :: t
        end function
    end interface

    type mytype
        real(WP) :: q
        procedure(integrand_interface), nopass, pointer :: integrand
    contains
        procedure :: init
    endtype

    class(mytype), pointer :: this_
contains

subroutine init( this )
    class(mytype), target :: this
    this % integrand => myfunc
    this_ => this
end subroutine

function myfunc(X, t) result( val )
    real(WP), dimension(:)       :: X
    real(WP), dimension(size(X)) :: val
    real(WP)                     :: t
    val = this_ % q * x
end function

end module


program main
    use mymod
    implicit none
    !integer, parameter :: WP = kind(1.0D0)
    type(mytype) :: mt
    real(WP), dimension(2) :: X
    real(WP) :: t, h

    call mt% init
    mt% q = 100.0

    X = (/ 0.0_WP, 1.0_WP /)
    t = 0.0_WP
    h = 1.0_WP

    call rk1(X, t, h, mt%integrand)
    !call integrator(0.0_WP,  mt% integrand, 1.0, 2.0, ans )
    print *, X

    contains
    subroutine rk1(X, t, h, f)
        ! Makes one step with the forward Euler method.
        ! Calculates next position using timestep h.
        ! See, e.g., Hairer, Nørsett and Wanner (1993, p. 51)
        implicit none
    !    integer, parameter :: WP = kind(1.0D0)
        ! IO
        real(WP), intent(inout), dimension(:)   :: X
        real(WP), intent(inout)                 :: t
        real(WP), intent(in)                    :: h
        interface
            function f(X, t) result(val)
                import WP
                real(WP), dimension(:)       :: X
                real(WP), dimension(size(X)) :: val
                real(WP)                     :: t
            end function
        end interface
        ! Local variables
        real(WP), dimension(size(X))            :: k1

        ! Evaluations of f(x, t)
        k1 = f(X, t)

        ! Next position
        X  = X + h*k1
        t  = t + h
    end subroutine


end program
