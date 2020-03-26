module integrator_module

    use parameters,          only: WP, DP
    use interpolator_module, only: interpolator

    implicit none
    private
    public :: integrate_fixed, integrate_variable, integrate_special
    public :: rk1, rk2, rk3, rk4
    public :: bs32, dp54, dp87


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
        ! Using these for debugging temporarily
        !real(WP) :: h_old
        !real(WP), dimension(size(X)) :: X0, X_old
        !X0 = X

        ! Initialise
        t = t0
        h = h0
        Naccepted = 0
        Nrejected = 0
        ! Evaluate k1 for first step (FSAL)
        ! Not all integrators can make use of this, but to keep the code simple
        ! they all take this as an argument.
        k1 = f % eval(X, t)

        ! Loop over timesteps until t == tmax
        do while (t < tmax)
            ! Adjust last timestep to stop exactly at tmax
            h = min(h, tmax - t)
            ! Keep track of these things for debugging for now
            !X_old = X
            !h_old = h
            ! Make step
            call method(X, t, h, f, k1, atol, rtol, accepted)
            ! Increment counters
            if (accepted) then
                Naccepted = Naccepted + 1
            else
                Nrejected = Nrejected + 1
            endif
            ! Finding out why x is sometimes out of area
            !if (f%debugflag) then
            !    print*, 'Debugging info:'
            !    print*, 'X0 = ', X0
            !    print*, 'X_old = ', X_old
            !    print*, 'h_old = ', h_old
            !endif
        end do
    end subroutine

    subroutine integrate_special(X, t0, tmax, h0, stoptimes, f, method, atol, rtol, Naccepted, Nrejected)
        ! Calculates X(t=tmax) as defined by the ODE x' = f(x, t),
        ! with initial value, X(t=0), given by X at input.
        ! Solution is found by repeatedly calling a variable-step
        ! Runge-Kutta method, using tolerance parameters atol and rtol
        ! to dynamically adjust the timestep, with initial timestep h0.
        ! The routine adjusts the last timestep to stop exactly at tmax,
        ! and returns the last position.
        ! Additionally, the routines takes a list of times, stoptimes, at which
        ! there are discontinuities in the data. Integration is stopped and
        ! restarted at these times.
        ! On output, the variables Naccepted and Nrejected contain the number
        ! of accepted and rejected steps respectively.
        implicit none
        ! IO
        real(WP), intent(inout), dimension(:)   :: X
        real(WP), intent(in),    dimension(:)   :: stoptimes
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
        real(WP) :: t, h, h_old
        integer  :: istop
        logical  :: accepted

        ! Initialise
        t = t0
        h = h0
        Naccepted = 0
        Nrejected = 0
        ! Set istop to the first index in stoptimes larger than t0.
        do istop = 1, size(stoptimes)
            if (stoptimes(istop) > t0) exit
        enddo
        ! Evaluate k1 for first step (FSAL)
        ! Not all integrators can make use of this, but to keep the code simple
        ! they all take this as an argument.
        k1 = f % eval(X, t)

        ! Loop over timesteps until t == tmax
        do while (t < tmax)
            ! Adjust last timestep to stop exactly at tmax
            h = min(h, tmax - t)
            ! Handle discontinuities
            if ( (t < stoptimes(istop)) .and. (stoptimes(istop) < (t + h)) ) then
                ! About to step across a discontinuity.
                ! Keep track of the current timestep.
                h_old = h
                ! Calculate remaining time until discontinuity.
                h = stoptimes(istop) - t
                ! Make step
                call method(X, t, h, f, k1, atol, rtol, accepted)
                ! Update book-keeping if step was accepted
                if (accepted) then
                    ! Set timestep back to old value
                    h = h_old
                    ! Increment istop to index of next discontinuity
                    istop = istop + 1
                endif
            else
                ! Normal step, away from any discontinuities
                ! Make step
                call method(X, t, h, f, k1, atol, rtol, accepted)
            endif
            ! Increment counters
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

    subroutine rk2(X, t, h, f)
        ! Make one step with the explicit trapezoid method.
        ! Calculates next position using timestep h.
        ! See, e.g., Griffiths (2010, pp. 44--45)
        implicit none
        ! IO
        real(WP), intent(inout), dimension(:)   :: X
        real(WP), intent(inout)                 :: t
        real(WP), intent(in)                    :: h
        type(interpolator), intent(inout)       :: f
        ! Local variables
        real(WP), dimension(size(X))            ::  k1, k2

        ! Evaluations of f(x, t)
        k1 = f % eval(X,        t)
        k2 = f % eval(X + k1*h, t+h)

        ! Next position
        X  = X + h*(k1 + k2)/2
        t  = t + h
    end subroutine

    subroutine rk3(X, t, h, f)
        ! Make one step with Kutta's method.
        ! Calculates next position using timestep h.
        ! See, e.g., Griffiths (2010, p. 131)
        implicit none
        ! IO
        real(WP), intent(inout), dimension(:)   :: X
        real(WP), intent(inout)                 :: t
        real(WP), intent(in)                    :: h
        type(interpolator), intent(inout)       :: f
        ! Local variables
        real(WP), dimension(size(X))            ::  k1, k2, k3

        ! Evaluations of f(x, t)
        k1 = f % eval(X,                   t)
        k2 = f % eval(X + k1*h/2,          t + h/2)
        k3 = f % eval(X - k1*h   + 2*k2*h, t + h)

        ! Next position
        X  = X + h*(k1 + 4*k2 + k3)/6
        t  = t + h
    end subroutine

    subroutine rk4(X, t, h, f)
        ! Make one step with the classic 4th-order Runge-Kutta method.
        ! Calculates next position using timestep h.
        ! See, e.g., Griffiths (2010, p. 131)
        implicit none
        ! IO
        real(WP), intent(inout), dimension(:)   :: X
        real(WP), intent(inout)                 :: t
        real(WP), intent(in)                    :: h
        type(interpolator), intent(inout)       :: f
        ! Local variables
        real(WP), dimension(size(X))            ::  k1, k2, k3, k4

        ! Evaluations of f(x, t)
        k1 = f % eval( X,          t       )
        k2 = f % eval( X + k1*h/2, t + h/2 )
        k3 = f % eval( X + k2*h/2, t + h/2 )
        k4 = f % eval( X + k3*h,   t + h   )

        ! Next position
        X  = X + h*(k1 + 2*k2 + 2*k3 + k4)/6
        t  = t + h
    end subroutine



    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!! Variable step Runge-Kutta methods !!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine adjust_timestep(X, X_, q, atol, rtol, h, h_new, accepted)
        ! Implementing error checking and variable stepsize approximately as in
        ! Hairer, Nørsett and Wanner:
        ! Solving ordinary differential equations I -- Nonstiff problems
        ! pages 167 and 168 in the 2008 printing

        ! parameters limiting change in timestep,
        ! defined in parameters.f90
        use parameters, only: fac, maxfac
        implicit none
        ! input and output
        real(WP), intent(in), dimension(:)  :: X, X_
        real(WP), intent(in)                :: q, atol, rtol, h
        real(WP), intent(out)               :: h_new
        logical,  intent(out)               :: accepted
        ! local variables
        real(WP), dimension(size(X))        :: sc
        real(WP)                            :: err, h_opt

        sc    = max(X, X_) * rtol + atol
        err   = sqrt( sum( ((X - X_) / sc)**2 ) )
        ! In some cases,  we can get X == X_ by chance,
        ! set minimum value to avoid division by zero.
        err   = max(err, 1e-12_wp)
        ! Calculating optimal h
        h_opt = h * (1.0_wp/err) ** (1.0_wp/(q+1.0_wp))
        ! Check if step is accepted
        if (err < 1.0_wp) then
            ! Step is accepted
            accepted = .true.
        else
            ! Step is rejected
            accepted = .false.
        endif
        ! Adjust steplength in both cases
        h_new = min(maxfac * h, fac * h_opt)
    end subroutine

    subroutine bs32(X, t, h, f, k1, atol, rtol, accepted)
        ! Bogacki-Shampine 3(2) method.
        ! See Bogacki & Shampine (1989)
        implicit none
        ! IO
        real(WP), intent(inout), dimension(:)   :: X, k1
        real(WP), intent(inout)                 :: t, h
        real(WP), intent(in)                    :: atol, rtol
        logical,  intent(out)                   :: accepted
        type(interpolator), intent(inout)       :: f
        ! Local variables
        ! variables needed for stepsize control
        real(WP), dimension(size(X))    :: X2, X3
        real(WP)                        :: q, h_new
        ! Evaluations of f(x, t)
        real(WP), dimension(size(X))    ::  k2, k3, k4

        ! Nodes and weights for method
        real(wp), parameter :: c2  = 1.0_wp / 2.0_wp
        real(wp), parameter :: c3  = 3.0_wp / 4.0_wp
        real(wp), parameter :: c4  = 1.0_wp
        ! Weights
        real(wp), parameter :: a21 = 1.0_wp / 2.0_wp
        real(wp), parameter :: a31 = 0.0_wp
        real(wp), parameter :: a32 = 3.0_wp / 4.0_wp
        real(wp), parameter :: a41 = 2.0_wp / 9.0_wp
        real(wp), parameter :: a42 = 1.0_wp / 3.0_wp
        real(wp), parameter :: a43 = 4.0_wp / 9.0_wp
        ! second order coefficients
        real(wp), parameter :: b21 = 7.0_wp / 24.0_wp
        real(wp), parameter :: b22 = 1.0_wp / 4.0_wp
        real(wp), parameter :: b23 = 1.0_wp / 3.0_wp
        real(wp), parameter :: b24 = 1.0_wp / 8.0_wp
        ! third order coefficients
        real(wp), parameter :: b31 = 2.0_wp / 9.0_wp
        real(wp), parameter :: b32 = 1.0_wp / 3.0_wp
        real(wp), parameter :: b33 = 4.0_wp / 9.0_wp

        ! k1 comes from previous k4 (FSAL)
        ! k1 = f % eval(X,                        t)
        k2 = f % eval(X + a21*k1*h,             t + c2*h)
        k3 = f % eval(X + a31*k1*h + a32*k2*h,  t + c3*h)

        ! Third order prediction
        X3 = X + h*(b31*k1 + b32*k2 + b33*k3)
        ! One more evalutation of f(x,t), at X3
        ! (this will be used as k1 at next iteration, due to FSAL)
        k4 = f % eval(X3, t + h)
        ! Second order prediction
        X2 = X + h*(b21*k1 + b22*k2 + b23*k3 + b24*k4)

        ! these are 2nd and 3rd order methods, q = min(2, 3)
        q  = 2.0_wp
        call adjust_timestep(X2, X3, q, atol, rtol, h, h_new, accepted)
        if (accepted) then
            ! Step is accepted, update, X, t and h
            ! (using 3rd order result to proceed)
            X = X3
            t = t + h
            h = h_new
            ! Return k4, to be used as k1 at next step (FSAL)
            k1 = k4
        else
            ! If step is not accepted, X and t remain unchanged,
            ! but h is adjusted (to a smaller value)
            h = h_new
        endif
    end subroutine

    subroutine dp54(X, t, h, f, k1, atol, rtol, accepted)
        ! Dormand-Prince 5(4) method.
        ! See Dormand & Prince (1980, 1986)

        ! parameters limiting change in timestep,
        ! defined in parameters.f90
        use parameters, only: fac, maxfac

        implicit none
        ! IO
        real(WP), intent(inout), dimension(:)   :: X, k1
        real(WP), intent(inout)                 :: t, h
        real(WP), intent(in)                    :: atol, rtol
        logical,  intent(out)                   :: accepted
        type(interpolator), intent(inout)       :: f
        ! Local variables
        ! variables needed for stepsize control
        real(WP), dimension(size(X))    :: X4, X5
        real(WP)                        :: q, h_new
        ! Evaluations of f(x, t)
        real(WP), dimension(size(X))    ::  k2, k3, k4, k5, k6, k7

        ! Nodes and weights for method
        real(WP), parameter ::c1  =  0.0_wp
        real(WP), parameter ::c2  =  1.0_wp / 5.0_wp
        real(WP), parameter ::c3  =  3.0_wp / 10.0_wp
        real(WP), parameter ::c4  =  4.0_wp / 5.0_wp
        real(WP), parameter ::c5  =  8.0_wp / 9.0_wp
        real(WP), parameter ::c6  =  1.0_wp
        real(WP), parameter ::c7  =  1.0_wp
        ! Weights
        real(WP), parameter ::a21 =  1.0_wp     / 5.0_wp
        real(WP), parameter ::a31 =  3.0_wp     / 40.0_wp
        real(WP), parameter ::a32 =  9.0_wp     / 40.0_wp
        real(WP), parameter ::a41 =  44.0_wp    / 45.0_wp
        real(WP), parameter ::a42 = -56.0_wp    / 15.0_wp
        real(WP), parameter ::a43 =  32.0_wp    / 9.0_wp
        real(WP), parameter ::a51 =  19372.0_wp / 6561.0_wp
        real(WP), parameter ::a52 = -25360.0_wp / 2187.0_wp
        real(WP), parameter ::a53 =  64448.0_wp / 6561.0_wp
        real(WP), parameter ::a54 = -212.0_wp   / 729.0_wp
        real(WP), parameter ::a61 =  9017.0_wp  / 3168.0_wp
        real(WP), parameter ::a62 = -355.0_wp   / 33.0_wp
        real(WP), parameter ::a63 =  46732.0_wp / 5247.0_wp
        real(WP), parameter ::a64 =  49.0_wp    / 176.0_wp
        real(WP), parameter ::a65 = -5103.0_wp  / 18656.0_wp
        real(WP), parameter ::a71 =  35.0_wp    / 384.0_wp
        real(WP), parameter ::a72 =  0.0_wp
        real(WP), parameter ::a73 =  500.0_wp   / 1113.0_wp
        real(WP), parameter ::a74 =  125.0_wp   / 192.0_wp
        real(WP), parameter ::a75 = -2187.0_wp  / 6784.0_wp
        real(WP), parameter ::a76 =  11.0_wp    / 84.0_wp
        ! Fourth order coefficients
        real(WP), parameter ::b41 =  5179.0_wp  / 57600.0_wp
        real(WP), parameter ::b42 =  0.0_wp
        real(WP), parameter ::b43 =  7571.0_wp  / 16695.0_wp
        real(WP), parameter ::b44 =  393.0_wp   / 640.0_wp
        real(WP), parameter ::b45 = -92097.0_wp / 339200.0_wp
        real(WP), parameter ::b46 =  187.0_wp   / 2100.0_wp
        real(WP), parameter ::b47 =  1.0_wp     / 40.0_wp
        ! Fifth order coefficients
        real(WP), parameter ::b51 =  35.0_wp    / 384.0_wp
        real(WP), parameter ::b52 =  0.0_wp
        real(WP), parameter ::b53 =  500.0_wp   / 1113.0_wp
        real(WP), parameter ::b54 =  125.0_wp   / 192.0_wp
        real(WP), parameter ::b55 = -2187.0_wp  / 6784.0_wp
        real(WP), parameter ::b56 =  11.0_wp    / 84.0_wp
        real(WP), parameter ::b57 =  0.0_wp

        ! Flag to indicate if step was accepted
        accepted = .false.
        ! k1 comes from previous k4 (FSAL)
        ! k1  = f % eval(X,                                                        t + c1*h )
        k2  = f % eval(X + a21*h*k1,                                             t + c2*h )
        k3  = f % eval(X + a31*h*k1 + a32*h*k2,                                  t + c3*h )
        k4  = f % eval(X + a41*h*k1 + a42*h*k2 + a43*h*k3,                       t + c4*h )
        k5  = f % eval(X + a51*h*k1 + a52*h*k2 + a53*h*k3 + a54*h*k4,            t + c5*h )
        k6  = f % eval(X + a61*h*k1 + a62*h*k2 + a63*h*k3 + a64*h*k4 + a65*h*k5, t + c6*h )

        ! Fifth order prediction (note b57 = 0)
        X5  = X + h*(k1*b51 + k2*b52 + k3*b53 + k4*b54 + k5*b55 + k6*b56)
        ! One more evalutation of f(x,t), at X5
        ! (this will be used as k1 at next iteration, due to FSAL)
        k7 = f % eval(X5, t + h)
        ! Fourth order prediction
        X4  = X + h*(k1*b41 + k2*b42 + k3*b43 + k4*b44 + k5*b45 + k6*b46 + k7*b47)

        ! these are 4th and 5th order methods, q = min(4, 5)
        q  = 4.0_wp
        call adjust_timestep(X4, X5, q, atol, rtol, h, h_new, accepted)
        if (accepted) then
            ! Step is accepted, update, X, t and h
            ! (using 5th order result to proceed)
            X = X5
            t = t + h
            h = h_new
            ! Return k7, to be used as k1 at next step (FSAL)
            k1 = k7
        else
            ! If step is not accepted, X and t remain unchanged,
            ! but h is adjusted (to a smaller value).
            h = h_new
        endif
    end subroutine


    subroutine dp87(X, t, h, f, k1_, atol, rtol, accepted)
        ! Dormand-Prince 8(7) method. Implemented after
        ! See Prince & Dormand (1981)
        implicit none
        ! IO
        real(WP), intent(inout), dimension(:)   :: X
        ! This method does not have FSAL,
        ! thus k1_ is just an unused argument to keep function
        ! signature consistent with the other methods.
        !
        real(WP), intent(inout), dimension(:)   :: k1_
        real(WP), intent(inout)                 :: t, h
        real(WP), intent(in)                    :: atol, rtol
        logical,  intent(out)                   :: accepted
        type(interpolator), intent(inout)       :: f
        ! Local variables
        ! variables needed for stepsize control
        real(WP), dimension(size(X)) :: X7, X8
        real(WP)                     :: q, h_new
        ! Evaluations of f(x, t)
        real(WP), dimension(size(X)) :: k1, k2, k3, k4, k5, k6, k7
        real(WP), dimension(size(X)) :: k8, k9, k10, k11, k12, k13

        ! Nodes
        real(WP), parameter ::c1    =             0.0_wp
        real(WP), parameter ::c2    =             1.0_wp /          18.0_wp
        real(WP), parameter ::c3    =             1.0_wp /          12.0_wp
        real(WP), parameter ::c4    =             1.0_wp /           8.0_wp
        real(WP), parameter ::c5    =             5.0_wp /          16.0_wp
        real(WP), parameter ::c6    =             3.0_wp /           8.0_wp
        real(WP), parameter ::c7    =            59.0_wp /         400.0_wp
        real(WP), parameter ::c8    =            93.0_wp /         200.0_wp
        real(WP), parameter ::c9    =    5490023248.0_wp /  9719169821.0_wp
        real(WP), parameter ::c10   =            13.0_wp /          20.0_wp
        real(WP), parameter ::c11   =    1201146811.0_wp /  1299019798.0_wp
        real(WP), parameter ::c12   =             1.0_wp
        real(WP), parameter ::c13   =             1.0_wp
        ! Weights
        real(WP), parameter ::a21   =             1.0_wp /          18.0_wp
        real(WP), parameter ::a31   =             1.0_wp /          48.0_wp
        real(WP), parameter ::a32   =             1.0_wp /          16.0_wp
        real(WP), parameter ::a41   =             1.0_wp /          32.0_wp
        real(WP), parameter ::a42   =             0.0_wp
        real(WP), parameter ::a43   =             3.0_wp /          32.0_wp
        real(WP), parameter ::a51   =             5.0_wp /          16.0_wp
        real(WP), parameter ::a52   =             0.0_wp
        real(WP), parameter ::a53   =           -75.0_wp /          64.0_wp
        real(WP), parameter ::a54   =            75.0_wp /          64.0_wp
        real(WP), parameter ::a61   =             3.0_wp /          80.0_wp
        real(WP), parameter ::a62   =             0.0_wp
        real(WP), parameter ::a63   =             0.0_wp
        real(WP), parameter ::a64   =             3.0_wp /          16.0_wp
        real(WP), parameter ::a65   =             3.0_wp /          20.0_wp
        real(WP), parameter ::a71   =      29443841.0_wp /   614563906.0_wp
        real(WP), parameter ::a72   =             0.0_wp
        real(WP), parameter ::a73   =             0.0_wp
        real(WP), parameter ::a74   =      77736538.0_wp /   692538347.0_wp
        real(WP), parameter ::a75   =     -28693883.0_wp /  1125000000.0_wp
        real(WP), parameter ::a76   =      23124283.0_wp /  1800000000.0_wp
        real(WP), parameter ::a81   =      16016141.0_wp /   946692911.0_wp
        real(WP), parameter ::a82   =             0.0_wp
        real(WP), parameter ::a83   =             0.0_wp
        real(WP), parameter ::a84   =      61564180.0_wp /   158732637.0_wp
        real(WP), parameter ::a85   =      22789713.0_wp /   633445777.0_wp
        real(WP), parameter ::a86   =     545815736.0_wp /  2771057229.0_wp
        real(WP), parameter ::a87   =    -180193667.0_wp /  1043307555.0_wp
        real(WP), parameter ::a91   =      39632708.0_wp /   573591083.0_wp
        real(WP), parameter ::a92   =             0.0_wp
        real(WP), parameter ::a93   =             0.0_wp
        real(WP), parameter ::a94   =    -433636366.0_wp /   683701615.0_wp
        real(WP), parameter ::a95   =    -421739975.0_wp /  2616292301.0_wp
        real(WP), parameter ::a96   =     100302831.0_wp /   723423059.0_wp
        real(WP), parameter ::a97   =     790204164.0_wp /   839813087.0_wp
        real(WP), parameter ::a98   =     800635310.0_wp /  3783071287.0_wp
        real(WP), parameter ::a101  =     246121993.0_wp /  1340847787.0_wp
        real(WP), parameter ::a102  =             0.0_wp
        real(WP), parameter ::a103  =             0.0_wp
        real(WP), parameter ::a104  =  -37695042795.0_wp / 15268766246.0_wp
        real(WP), parameter ::a105  =    -309121744.0_wp /  1061227803.0_wp
        real(WP), parameter ::a106  =     -12992083.0_wp /   490766935.0_wp
        real(WP), parameter ::a107  =    6005943493.0_wp /  2108947869.0_wp
        real(WP), parameter ::a108  =     393006217.0_wp /  1396673457.0_wp
        real(WP), parameter ::a109  =     123872331.0_wp /  1001029789.0_wp
        real(WP), parameter ::a111  =   -1028468189.0_wp /   846180014.0_wp
        real(WP), parameter ::a112  =             0.0_wp
        real(WP), parameter ::a113  =             0.0_wp
        real(WP), parameter ::a114  =    8478235783.0_wp /   508512852.0_wp
        real(WP), parameter ::a115  =    1311729495.0_wp /  1432422823.0_wp
        real(WP), parameter ::a116  =  -10304129995.0_wp /  1701304382.0_wp
        real(WP), parameter ::a117  =  -48777925059.0_wp /  3047939560.0_wp
        real(WP), parameter ::a118  =   15336726248.0_wp /  1032824649.0_wp
        real(WP), parameter ::a119  =  -45442868181.0_wp /  3398467696.0_wp
        real(WP), parameter ::a1110 =    3065993473.0_wp /   597172653.0_wp
        real(WP), parameter ::a121  =     185892177.0_wp /   718116043.0_wp
        real(WP), parameter ::a122  =             0.0_wp
        real(WP), parameter ::a123  =             0.0_wp
        real(WP), parameter ::a124  =   -3185094517.0_wp /   667107341.0_wp
        real(WP), parameter ::a125  =    -477755414.0_wp /  1098053517.0_wp
        real(WP), parameter ::a126  =    -703635378.0_wp /   230739211.0_wp
        real(WP), parameter ::a127  =    5731566787.0_wp /  1027545527.0_wp
        real(WP), parameter ::a128  =    5232866602.0_wp /   850066563.0_wp
        real(WP), parameter ::a129  =   -4093664535.0_wp /   808688257.0_wp
        real(WP), parameter ::a1210 =    3962137247.0_wp /  1805957418.0_wp
        real(WP), parameter ::a1211 =      65686358.0_wp /   487910083.0_wp
        real(WP), parameter ::a131  =     403863854.0_wp /   491063109.0_wp
        real(WP), parameter ::a132  =             0.0_wp
        real(WP), parameter ::a133  =             0.0_wp
        real(WP), parameter ::a134  =   -5068492393.0_wp /   434740067.0_wp
        real(WP), parameter ::a135  =    -411421997.0_wp /   543043805.0_wp
        real(WP), parameter ::a136  =     652783627.0_wp /   914296604.0_wp
        real(WP), parameter ::a137  =   11173962825.0_wp /   925320556.0_wp
        real(WP), parameter ::a138  =  -13158990841.0_wp /  6184727034.0_wp
        real(WP), parameter ::a139  =    3936647629.0_wp /  1978049680.0_wp
        real(WP), parameter ::a1310 =    -160528059.0_wp /   685178525.0_wp
        real(WP), parameter ::a1311 =     248638103.0_wp /  1413531060.0_wp
        real(WP), parameter ::a1312 =             0.0_wp
        ! Seventh order coefficients
        real(WP), parameter ::b71  =       13451932.0_wp /   455176623.0_wp
        real(WP), parameter ::b72  =              0.0_wp
        real(WP), parameter ::b73  =              0.0_wp
        real(WP), parameter ::b74  =              0.0_wp
        real(WP), parameter ::b75  =              0.0_wp
        real(WP), parameter ::b76  =     -808719846.0_wp /   976000145.0_wp
        real(WP), parameter ::b77  =     1757004468.0_wp /  5645159321.0_wp
        real(WP), parameter ::b78  =      656045339.0_wp /   265891186.0_wp 
        real(WP), parameter ::b79  =    -3867574721.0_wp /  1518517206.0_wp
        real(WP), parameter ::b710 =      465885868.0_wp /   322736535.0_wp
        real(WP), parameter ::b711 =       53011238.0_wp /   667516719.0_wp
        real(WP), parameter ::b712 =              2.0_wp /          45.0_wp
        real(WP), parameter ::b713 =              0.0_wp
        ! Eigth order coefficients
        real(WP), parameter ::b81  =       14005451.0_wp /   335480064.0_wp
        real(WP), parameter ::b82  =              0.0_wp
        real(WP), parameter ::b83  =              0.0_wp
        real(WP), parameter ::b84  =              0.0_wp
        real(WP), parameter ::b85  =              0.0_wp
        real(WP), parameter ::b86  =      -59238493.0_wp /  1068277825.0_wp
        real(WP), parameter ::b87  =      181606767.0_wp /   758867731.0_wp
        real(WP), parameter ::b88  =      561292985.0_wp /   797845732.0_wp
        real(WP), parameter ::b89  =    -1041891430.0_wp /  1371343529.0_wp
        real(WP), parameter ::b810 =      760417239.0_wp /  1151165299.0_wp
        real(WP), parameter ::b811 =      118820643.0_wp /   751138087.0_wp
        real(WP), parameter ::b812 =     -528747749.0_wp /  2220607170.0_wp
        real(WP), parameter ::b813 =              1.0_wp /           4.0_wp

        ! Evaluations of f(x,t)
k1  = f%eval(X,                                                                                                          t+ c1*h)
k2  = f%eval(X+h*( a21*k1),                                                                                              t+ c2*h)
k3  = f%eval(X+h*( a31*k1+ a32*k2),                                                                                      t+ c3*h)
k4  = f%eval(X+h*( a41*k1+ a42*k2+ a43*k3),                                                                              t+ c4*h)
k5  = f%eval(X+h*( a51*k1+ a52*k2+ a53*k3+ a54*k4),                                                                      t+ c5*h)
k6  = f%eval(X+h*( a61*k1+ a62*k2+ a63*k3+ a64*k4+ a65*k5),                                                              t+ c6*h)
k7  = f%eval(X+h*( a71*k1+ a72*k2+ a73*k3+ a74*k4+ a75*k5+ a76*k6),                                                      t+ c7*h)
k8  = f%eval(X+h*( a81*k1+ a82*k2+ a83*k3+ a84*k4+ a85*k5+ a86*k6+ a87*k7),                                              t+ c8*h)
k9  = f%eval(X+h*( a91*k1+ a92*k2+ a93*k3+ a94*k4+ a95*k5+ a96*k6+ a97*k7+ a98*k8),                                      t+ c9*h)
k10 = f%eval(X+h*(a101*k1+a102*k2+a103*k3+a104*k4+a105*k5+a106*k6+a107*k7+a108*k8+a109*k9),                              t+c10*h)
k11 = f%eval(X+h*(a111*k1+a112*k2+a113*k3+a114*k4+a115*k5+a116*k6+a117*k7+a118*k8+a119*k9+a1110*k10),                    t+c11*h)
k12 = f%eval(X+h*(a121*k1+a122*k2+a123*k3+a124*k4+a125*k5+a126*k6+a127*k7+a128*k8+a129*k9+a1210*k10+a1211*k11),          t+c12*h)
k13 = f%eval(X+h*(a131*k1+a132*k2+a133*k3+a134*k4+a135*k5+a136*k6+a137*k7+a138*k8+a139*k9+a1310*k10+a1311*k11+a1312*k12),t+c13*h)

        ! Calculate seventh order prediction (note b72, b73, b74, b75 and b713 are zero )
        X7  = X + h*(k1*b71 + k6*b76 + k7*b77 + k8*b78 + k9*b79 + k10*b710 + k11*b711 + k12*b712)
        ! Calculate eighth order prediction (note b82, b83, b84 and b85 are zero )
        X8  = X + h*(k1*b81 + k6*b86 + k7*b87 + k8*b88 + k9*b89 + k10*b810 + k11*b811 + k12*b812 + k13*b813)

        ! these are 7th and 8th order methods, q = min(7, 8)
        q   = 7.0_wp
        call adjust_timestep(X7, X8, q, atol, rtol, h, h_new, accepted)
        if (accepted) then
            ! Step is accepted, update, X, t and h
            ! (using 8th order result to proceed)
            X = X8
            t = t + h
            h = h_new
        else
            ! If step is not accepted, X and t remain unchanged,
            ! but h is adjusted (to a smaller value)
            h = h_new
        endif
    end subroutine

end module
