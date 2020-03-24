module input_module
    use parameters, only: SP, DP, WP
    implicit none

    contains

        subroutine check(iostatus, filename, lineno, message)
            integer,            intent(in)                  :: iostatus
            character(len=*),   intent(in)                  :: filename
            integer,            intent(in),     optional    :: lineno
            character(len=*),   intent(in),     optional    :: message

            if (iostatus < 0) then
                if (present(lineno)) then
                    print*, 'End of file ', filename, ' reached after line ', lineno
                    if (present(message)) then
                        print*, message
                    endif
                else
                    print*, 'End of file ', filename, ' reached sooner than expected'
                    if (present(message)) then
                        print*, message
                    endif
                endif
                stop
            else if (iostatus > 0) then
                if (present(lineno)) then
                    print*, 'Error reading ', filename, ' after line ', lineno
                    print*, 'iostat = ', iostatus
                    if (present(message)) then
                        print*, message
                    endif
                else
                    print*, 'Error reading ', filename
                    print*, 'iostat = ', iostatus
                    if (present(message)) then
                        print*, message
                    endif
                endif
                stop
            else
                ! in this case, iostat = 0, and all is well in the world
                continue
            endif
        end subroutine


        subroutine read_initial_positions(filename, X0)
            implicit none
            ! Input
            character(len=*),                       intent(in)  :: filename
            ! Output
            real(WP), allocatable, dimension(:,:),  intent(out) :: X0
            ! local variables
            integer :: Np, i, iostatus
            ! Number of dimensions, hardcoding for now
            integer, parameter :: Ndims = 2

            ! Open input file
            open(42, file = trim(filename))
            ! Read number of particle positions from first line in file
            read(42, *, iostat = iostatus) Np
            ! Confirm that reading went well
            call check(iostatus, filename, message = 'Reading Number of particles')
            ! Allocate space to read values from file
            allocate(X0(Ndims, Np))

            do i = 1, Np
                read(42, *, iostat = iostatus) X0(:,i)
                call check(iostatus, filename, lineno = i)
            enddo
        end subroutine

end module
