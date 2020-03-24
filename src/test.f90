program test

    use parameters,          only: WP, DP
    use interpolator_module, only: interpolator
    use integrator_module,   only: rk1, bs32
    use integrator_module,   only: integrate_fixed, integrate_variable
    use experiment_module,   only: experiment_fixed, experiment_variable

    implicit none

    real(WP), dimension(11)         :: xc, yc, tc
    real(WP), dimension(11, 11, 11) :: gvx, gvy

    integer(WP) :: i, j, na, nr
    integer(DP) :: Nsteps
    real(WP) :: t0, tmax, h, atol, rtol
    type(interpolator) :: f

    real(WP), dimension(2,10) :: X
    real(WP), dimension(7) :: timesteps
    real(WP), dimension(7) :: tolerances

    timesteps = (/ 1.0_WP, 0.5_WP, 0.2_WP, 0.1_WP, 0.05_WP, 0.02_WP, 0.01_WP /)
    tolerances = (/ 1e-4_WP, 5e-5_WP, 2e-5_WP, 1e-5_WP, 5e-6_WP, 2e-6_WP, 1e-6_WP /)


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

    X(1,:) =  4
    x(2,:) = -3
    call experiment_fixed(X, t0, tmax, timesteps, f, rk1, 'testfile_fixed.hdf5')
    call experiment_variable(X, t0, tmax, tolerances, f, bs32, 'testfile_variable.hdf5')


end program
