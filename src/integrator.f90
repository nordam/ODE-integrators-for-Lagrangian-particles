module integrator_module

    use parameters,          only: WP, DP
    use interpolator_module, only: interpolator

    implicit none
    private
    public :: integrate_fixed, integrate_variable
    public :: rk1, bs32


    contains

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!     Subroutines to integrate ODEs     !!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine integrate_fixed(X, t0, tmax, h0, f, method, Nsteps)
        ! Calculates X(t=tmax) as defined by the ODE x' = f(x, t),
        ! with initial value, X(t=0), given by X at input.
        ! Solution is found by repeatedly calling a fixed-step
        ! Runge-Kutta method, with timestep h0.
        ! The routine adjusts the last timestep to stop exactly at tmax,
        ! and returns the last position.
        implicit none
        ! IO
        real(WP), intent(inout), dimension(:)   :: X
        real(WP), intent(in)                    :: t0, tmax, h0
        type(interpolator), intent(inout)       :: f
        integer(DP), intent(out)                :: Nsteps
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
        Nsteps = 0
        ! Loop over timesteps until t == tmax.
        do while (t < tmax)
            ! If remaining time until tmax is smaller than timestep,
            ! adjust h to stop exactly at tmax.
            h = min(h, tmax - t)
            call method(X, t, h, f)
            Nsteps = Nsteps + 1
        end do
    end subroutine

    subroutine integrate_variable(X, t0, tmax, h0, f, method, atol, rtol, Naccepted, Nrejected)
        ! Calculates X(t=tmax) as defined by the ODE x' = f(x, t),
        ! with initial value, X(t=0), given by X at input.
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

        ! Initialise
        t = t0
        h = h0
        Naccepted = 0
        Nrejected = 0
        ! Evaluate k1 for first step (FSAL)
        k1 = f % eval(X, t)

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



    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!! Variable step Runge-Kutta methods !!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine bs32(X, t, h, f, k1, atol, rtol, accepted)
        ! Bogacki-Shampine 3(2) method.
        ! parameters limiting change in timestep
        use parameters, only: fac, maxfac
        implicit none
        ! IO
        real(WP), intent(inout), dimension(:)   :: X, k1
        real(WP), intent(inout)                 :: t, h
        real(WP), intent(in)                    :: atol, rtol
        type(interpolator), intent(inout)       :: f
        logical, intent(out)                    :: accepted
        ! Local variables
        ! variables needed for stepsize control
        real(WP), dimension(size(X)) :: X2, X3, sc
        real(WP)               :: q, err, h_opt
        ! Evaluations of f(x, t)
        real(WP), dimension(size(X)) ::  k2, k3, k4
        ! Nodes and weights for method
        real(wp), parameter :: c2  = 1.0_wp / 2.0_wp
        real(wp), parameter :: c3  = 3.0_wp / 4.0_wp
        real(wp), parameter :: c4  = 1.0_wp
        real(wp), parameter :: a21 = 1.0_wp / 2.0_wp
        real(wp), parameter :: a31 = 0.0_wp
        real(wp), parameter :: a32 = 3.0_wp / 4.0_wp
        real(wp), parameter :: a41 = 2.0_wp / 9.0_wp
        real(wp), parameter :: a42 = 1.0_wp / 3.0_wp
        real(wp), parameter :: a43 = 4.0_wp / 9.0_wp
        ! third order coefficients
        real(wp), parameter :: b31 = 2.0_wp / 9.0_wp
        real(wp), parameter :: b32 = 1.0_wp / 3.0_wp
        real(wp), parameter :: b33 = 4.0_wp / 9.0_wp
        ! second order coefficients
        real(wp), parameter :: b21 = 7.0_wp / 24.0_wp
        real(wp), parameter :: b22 = 1.0_wp / 4.0_wp
        real(wp), parameter :: b23 = 1.0_wp / 3.0_wp
        real(wp), parameter :: b24 = 1.0_wp / 8.0_wp

        accepted = .false.
        ! k1 comes from previous k4 (FSAL)
        ! k1 = f % eval(X,                        t)
        k2  = f % eval(X + a21*k1*h,             t + c2*h)
        k3  = f % eval(X + a31*k1*h + a32*k2*h,  t + c3*h)
        ! Third order prediction
        X3  = X + h*(b31*k1 + b32*k2 + b33*k3)
        ! One more slope, at X3 (FSAL)
        k4  = f % eval(X3, t + h)
        ! Second order prediction
        X2  = X + h*(b21*k1 + b22*k2 + b23*k3 + b24*k4)

        ! Implementing error checking and variable stepsize roughly as in
        ! Hairer, Nørsett and Wanner:
        ! Solving ordinary differential equations I -- Nonstiff problems
        ! pages 167 and 168 in the 2008 printing

        ! these are 2nd and 3rd order methods, q = min(2, 3)
        q      = 2.0_wp
        sc(1)  = max(abs(X2(1)), abs(X3(1))) * rtol + atol
        sc(2)  = max(abs(X2(2)), abs(X3(2))) * rtol + atol
        err    = sqrt( sum( ((X3 - X2) / sc)**2 ) )
        ! For certain equations, it is possible to get X2 == X3,
        ! hence set minimum err to avoid division by zero.
        err    = max(err, 1e-16_wp)
        ! Calculating optimal h
        h_opt  = h * (1.0_wp/err) ** (1.0_wp/(q+1.0_wp))
        ! Check if step is accepted
        if (err < 1.0_wp) then
            ! Step is accepted, use 3rd order result
            accepted = .true.
            X  = X3
            t  = t + h
            ! Return k4, to be used as k1 at next step (FSAL)
            k1 = k4
            ! If the step is rejected,
            ! X, t and k1 remain unchanged
        endif
        ! Adjust steplength in both cases
        h = min(maxfac * h, fac * h_opt)
    end subroutine

end module
