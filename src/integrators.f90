module integrator_module

    use parameters,          only: WP
    use interpolator_module, only: interpolator

    implicit none
    private
    public :: integrate_fixed
    public :: rk1


    contains



    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!     Subroutines to integrate ODEs     !!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine integrate_fixed(X, t0, tmax, h0, f, method)
        ! Calculates X(tmax) as defined by the ODE x' = f(x, t),
        ! with initial value, X(0), given by X at input.
        ! Solution is found by repeatedly calling a fixed-step
        ! Runge-Kutta method, with timestep h0.
        ! The routine adjusts the last timestep to stop exactly at tmax,
        ! and returns the last position.
        implicit none
        ! IO
        real(WP), intent(inout), dimension(:)   :: X
        real(WP), intent(in)                    :: t0, tmax, h0
        type(interpolator), intent(inout)       :: f
        interface
            subroutine method(X, t, h, f)
                import :: WP, interpolator
                real(WP), intent(inout), dimension(:)   :: X
                real(WP), intent(inout)                 :: t
                real(WP), intent(in)                    :: h
                type(interpolator), intent(inout)       :: f
            end subroutine method
        end interface
        ! local variables
        real(WP) :: t, h

        t = t0
        h = h0
        ! Loop over timesteps until t == tmax.
        do while (t < tmax)
            ! If remaining time until tmax is smaller than timestep,
            ! adjust h to stop exactly at tmax.
            h = min(h, tmax - t)
            call method(X, t, h, f)
        end do
    end subroutine

    subroutine integrate_variable(X, t0, tmax, h0, f, method, atol, rtol, Naccepted, Nrejected)
        ! Calculates X(tmax) as defined by the ODE x' = f(x, t),
        ! with initial value, X(0), given by X at input.
        ! Solution is found by repeatedly calling a variable-step
        ! Runge-Kutta method, using tolerance parameters atol and rtol
        ! to dynamically adjust the timestep, with initial timestep h0.
        ! The routine adjusts the last timestep to stop exactly at tmax,
        ! and returns the last position.
        ! On output, the variables Naccepted and Nrejected contain the number
        ! of accepted and rejected steps respectively.
        implicit none
        ! IO
        real(WP), intent(inout), dimension(:)   :: X
        real(WP), intent(in)                    :: t0, tmax, h0, atol, rtol
        type(interpolator), intent(inout)       :: f
        integer(WP), intent(out)                :: Naccepted, Nrejected
        interface
            subroutine method(X, t, h, f, k1, atol, rtol, accepted)
                import :: WP, interpolator
                real(WP), intent(inout), dimension(:)   :: X, k1
                real(WP), intent(inout)                 :: t, h
                real(WP), intent(in)                    :: atol, rtol
                type(interpolator), intent(inout)       :: f
                logical, intent(out)                    :: accepted
            end subroutine method
        end interface
        ! local variables
        real(WP), dimension(size(X)) :: k1
        real(WP) :: t, h
        logical  :: accepted

        t = t0
        h = h0
        Naccepted = 0
        Nrejected = 0
        ! Loop over timesteps until t == tmax
        do while (t < tmax)
            ! Adjust last timestep to stop exactly at tmax
            h = min(h, tmax - t)
            call method(X, t, h, f, k1, atol, rtol, accepted)
            if (accepted) then
                Naccepted = Naccepted + 1
            else
                Nrejected = Nrejected + 1
            endif
        end do
    end subroutine



    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!! Fixed step Runge-Kutta methods !!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine rk1(X, t, h, f)
        ! Makes one step with the forward Euler method.
        ! Calculates next position using timestep h.
        ! See, e.g., Hairer, Nørsett and Wanner (1993, p. 51)
        implicit none
        ! IO
        real(WP), intent(inout), dimension(:)   :: X
        real(WP), intent(inout)                 :: t
        real(WP), intent(in)                    :: h
        type(interpolator), intent(inout)       :: f
        ! Local variables
        real(WP), dimension(size(X))            :: k1

        ! Evaluations of f(x, t)
        k1 = f % eval(X, t)

        ! Next position
        X  = X + h*k1
        t  = t + h
    end subroutine

end module
