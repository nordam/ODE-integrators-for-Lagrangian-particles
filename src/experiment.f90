module experiment_module

    use hdf5
    use h5lt
    use parameters,          only: SP, DP, WP
    use interpolator_module, only: interpolator
    use integrator_module,   only: integrate_fixed, integrate_variable, integrate_special
    use output_module,       only: write_to_hdf5, create_hdf5_file, close_hdf5_file

    implicit none
    private
    public :: experiment_fixed, experiment_variable, experiment_special

    contains

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!! Subroutine to run experiment with fixed-step Runge-Kutta method !!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine experiment_fixed(X0, t0, tmax, timesteps, f, method, outputfile)
        ! This subroutine takes the initial positions at time t0, for a number,
        ! Np, of particles, along with a list of timesteps, an interpolator
        ! object to evaluate the velocity field, and a numerical ODE integrator.
        ! It then obtains the position at time tmax for each
        ! particle, writes to file, and repeats for each timestep.
        implicit none
        ! Input variables
        real(WP), dimension(:,:), intent(in)    :: X0
        real(WP), dimension(:),   intent(in)    :: timesteps
        real(WP),                 intent(in)    :: t0, tmax
        character(len=*),         intent(in)    :: outputfile
        type(interpolator),       intent(inout) :: f
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
        real(WP), dimension(:,:), allocatable   :: X
        integer(hid_t)                          :: file_id
        integer                                 :: n, idt, Np
        integer(DP), dimension(size(X0,2))      :: Nsteps
        real(WP)                                :: h, tic, toc

        ! Assuming that X0 has shape (2, Np) or (3, Np)
        ! depending on two or three dimensions
        Np = size(X0, 2)
        allocate(X( size(X0,1), size(X0,2) ))

        ! Create file for output
        call create_hdf5_file(outputfile, file_id)

        ! Scan through timesteps
        do idt = 1, size(timesteps)
            h = timesteps(idt)
            ! Reset initial positions
            X = X0
            ! Measure computational time
            call cpu_time(tic)
            ! Transport each particle from time t0 to tmax
            do n = 1, Np
                call integrate_fixed(X(:,n), t0, tmax, h, f, method, Nsteps(n))
            end do
            ! Measure computational time
            call cpu_time(toc)
            ! Write end positions to hdf5 file
            call write_to_hdf5(file_id, X, timesteps(idt))
            ! Print information on timing and number of steps:
            ! filename, timestep, runtime, accepted steps, rejected steps
            print*, trim(outputfile), timesteps(idt), toc - tic, sum(Nsteps), 0
            flush(6)
        end do
        ! Clean up
        deallocate(X)
        call close_hdf5_file(file_id)
    end subroutine



    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!! Subroutine to run experiment with variable-step Runge-Kutta method !!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine experiment_variable(X0, t0, tmax, tolerances, f, method, outputfile, h0input)
        ! This subroutine takes the initial positions at time t0, for a number,
        ! Np, of particles, along with a list of timesteps, an interpolator
        ! object to evaluate the velocity field, and a numerical ODE integrator.
        ! It then obtains the position at time tmax for each
        ! particle, writes to file, and repeats for each timestep.
        implicit none
        ! Input variables
        real(WP), dimension(:,:), intent(in)    :: X0
        real(WP), dimension(:),   intent(in)    :: tolerances
        real(WP),                 intent(in)    :: t0, tmax
        real(WP), optional,       intent(in)    :: h0input
        character(len=*),         intent(in)    :: outputfile
        type(interpolator),       intent(inout) :: f
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
        real(WP), dimension(:,:), allocatable   :: X
        integer(hid_t)                          :: file_id
        integer                                 :: n, idt, Np
        integer(DP), dimension(size(X0,2))      :: Naccepted, Nrejected
        real(WP)                                :: h0, tol, tic, toc

        ! Assuming that X0 has shape (2, Np) or (3, Np)
        ! depending on two or three dimensions
        Np = size(X0, 2)
        allocate(X( size(X0,1), size(X0,2) ))

        ! Create file for output
        call create_hdf5_file(outputfile, file_id)

        ! Initial timestep, setting default value if not supplied
        if (present(h0input)) then
            h0 = h0input
        else
            h0 = (Tmax - t0)/100
        endif

        ! Scan through tolerances
        do idt = 1, size(tolerances)
            ! Set tolerances (using atol = rtol for now)
            tol = tolerances(idt)
            ! Reset initial positions
            X = X0
            ! Measure computational time
            call cpu_time(tic)
            ! Transport each particle from time t0 to tmax
            do n = 1, Np
                call integrate_variable(X(:,n), t0, tmax, h0, f, method, &
                        tol, tol, Naccepted(n), Nrejected(n))
            end do
            ! Measure computational time
            call cpu_time(toc)
            ! Write end positions to hdf5 file
            call write_to_hdf5(file_id, X, tolerances(idt))
            ! Print information on timing and number of steps:
            ! filename, tolerance, runtime, accepted steps, rejected steps
            print*, trim(outputfile), tol, toc - tic, sum(Naccepted), sum(Nrejected)
            flush(6)
        end do
        ! Clean up
        deallocate(X)
        call close_hdf5_file(file_id)
    end subroutine


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!! Subroutine to run experiment with variable-step Runge-Kutta method !!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine experiment_special(X0, t0, tmax, stoptimes, tolerances, f, method, outputfile, h0input)
        ! This subroutine takes the initial positions at time t0, for a number,
        ! Np, of particles, along with a list of timesteps, an interpolator
        ! object to evaluate the velocity field, and a numerical ODE integrator.
        ! It then obtains the position at time tmax for each
        ! particle, writes to file, and repeats for each timestep.
        implicit none
        ! Input variables
        real(WP), dimension(:,:), intent(in)    :: X0
        real(WP), dimension(:),   intent(in)    :: stoptimes
        real(WP), dimension(:),   intent(in)    :: tolerances
        real(WP),                 intent(in)    :: t0, tmax
        real(WP), optional,       intent(in)    :: h0input
        character(len=*),         intent(in)    :: outputfile
        type(interpolator),       intent(inout) :: f
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
        real(WP), dimension(:,:), allocatable   :: X
        integer(hid_t)                          :: file_id
        integer                                 :: n, idt, Np
        integer(DP), dimension(size(X0,2))      :: Naccepted, Nrejected
        real(WP)                                :: h0, tol, tic, toc

        ! Assuming that X0 has shape (2, Np) or (3, Np)
        ! depending on two or three dimensions
        Np = size(X0, 2)
        allocate(X( size(X0,1), size(X0,2) ))

        ! Create file for output
        call create_hdf5_file(outputfile, file_id)

        ! Initial timestep, setting default value if not supplied
        if (present(h0input)) then
            h0 = h0input
        else
            h0 = (Tmax - t0)/100
        endif

        ! Scan through tolerances
        do idt = 1, size(tolerances)
            ! Set tolerances (using atol = rtol for now)
            tol = tolerances(idt)
            ! Reset initial positions
            X = X0
            ! Measure computational time
            call cpu_time(tic)
            ! Transport each particle from time t0 to tmax
            do n = 1, Np
                call integrate_special(X(:,n), t0, tmax, h0, stoptimes, f, method, &
                        tol, tol, Naccepted(n), Nrejected(n))
            end do
            ! Measure computational time
            call cpu_time(toc)
            ! Write end positions to hdf5 file
            call write_to_hdf5(file_id, X, tolerances(idt))
            ! Print information on timing and number of steps:
            ! filename, tolerance, runtime, accepted steps, rejected steps
            print*, trim(outputfile), tol, toc - tic, sum(Naccepted), sum(Nrejected)
            flush(6)
        end do
        ! Clean up
        deallocate(X)
        call close_hdf5_file(file_id)
    end subroutine

end module
