program test

    use parameters,          only: WP
    use interpolator_module, only: interpolator
    use integrator_module,   only: integrate_fixed, rk1

    implicit none

    real(WP), dimension(11)         :: xc, yc, tc
    real(WP), dimension(11, 11, 11) :: gvx, gvy

    integer(WP) :: i, j
    real(WP) :: t0, tmax, h
    type(interpolator) :: f

    real(WP), dimension(2) :: X

    xc = (/ -5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5 /)
    yc = xc
    tc = xc

    do i = 1, 11
        do j = 1, 11
            gvx(:,i,j) = -xc
            gvy(i,:,j) = -yc
        end do
    end do

    call f%init(xc, yc, tc, gvx, gvy, 2)

    X = (/ 4, -4 /)
    h = 0.0001
    t0 = -4
    tmax = 4

    print*, 'Analytical: ', X * exp(-(tmax - t0))
    call integrate_fixed(X, t0, tmax, h, f, rk1)
    print*, 'Numerical:  ', X



end program
