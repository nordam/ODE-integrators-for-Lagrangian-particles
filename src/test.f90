program test

    use parameters,          only: WP, DP
    use interpolator_module, only: interpolator
    use integrator_module,   only: rk1, bs32
    use integrator_module,   only: integrate_fixed, integrate_variable
    use experiment_module,   only: experiment_fixed

    implicit none

    real(WP), dimension(11)         :: xc, yc, tc
    real(WP), dimension(11, 11, 11) :: gvx, gvy

    integer(WP) :: i, j, na, nr
    integer(DP) :: Nsteps
    real(WP) :: t0, tmax, h, atol, rtol
    type(interpolator) :: f

    real(WP), dimension(2,10) :: X
    real(WP), dimension(7) :: timesteps

    timesteps = (/ 1.0, 0.5, 0.2, 0.1, 0.05, 0.02, 0.01 /)


    xc = (/ -5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5 /)
    yc = xc
    tc = xc

    do i = 1, 11
        do j = 1, 11
            gvx(:,i,j) = -xc
            gvy(i,:,j) = -yc
        end do
    end do

    call f%init(xc, yc, tc, gvx, gvy, 4)

    X(1,:) =  4
    x(2,:) = -3
    h = 0.001
    t0 = -4
    tmax = 4

    print*, 'Fixed step'
    print*, 'Analytical: ', X(:,1) * exp(-(tmax - t0))
    call integrate_fixed(X(:,1), t0, tmax, h, f, rk1, Nsteps)
    print*, 'Numerical:  ', X(:,1)
    print*, ''

    X(1,:) =  4
    x(2,:) = -3
    h = 1.00001
    atol = 1e-10
    rtol = 1e-10
    print*, 'Variable step'
    print*, 'Analytical: ', X(:,1) * exp(-(tmax - t0))
    call integrate_variable(X(:,1), t0, tmax, h, f, bs32, atol, rtol, Na, Nr)
    print*, 'Numerical:  ', X(:,1)
    print*, 'Accepted steps: ', na, 'Rejected steps: ', nr

    call experiment_fixed(X, t0, tmax, timesteps, f, rk1, 'testfile.hdf5')


end program
