module currentdata_module
    use netcdf
    use parameters, only: WP, DP, SP
    implicit none
contains

    subroutine get_current(filename, u, v, xc, yc, tc)
        ! This subroutine will read the u and v variables
        ! from a netCDF file.
        ! A number of assumptions are made:
        !   The variables are assumed to be called u and v in the file.
        !   The dimensions are assumed to be called X, Y, depth and time.
        !   The depth dimension is assumed to have length 1.
        implicit none
        ! Input
        character(len=*), intent(in) :: filename
        ! Output
        real(WP), dimension(:),     allocatable :: xc, yc, tc
        real(WP), dimension(:,:,:), allocatable :: u, v
        ! Local variables
        real(SP) :: scale_factor, add_offset, fill_value
        real(SP), allocatable, dimension(:,:,:) :: u_, v_
        ! file and variable ids, etc
        integer :: nx, ny, nt
        integer :: ncid
        integer :: u_id, v_id, t_id, x_id, y_id
        integer :: t_dim_id, x_dim_id, y_dim_id
        integer, dimension(4) :: start, count
        ! variable/dimension names
        character(len=*), parameter  :: x_name = 'X'
        character(len=*), parameter  :: y_name = 'Y'
        character(len=*), parameter  :: v_name = 'v'
        character(len=*), parameter  :: u_name = 'u'
        character(len=*), parameter  :: t_name = 'time'

        ! Open file for reading
        call check( nf90_open(filename, nf90_nowrite, ncid) )

        ! Get the ids of the dimensions.
        call check( nf90_inq_dimid(ncid, x_name, x_dim_id) )
        call check( nf90_inq_dimid(ncid, y_name, y_dim_id) )
        call check( nf90_inq_dimid(ncid, t_name, t_dim_id) )
        ! Get the ids of the variables.
        call check( nf90_inq_varid(ncid, x_name, x_id) )
        call check( nf90_inq_varid(ncid, y_name, y_id) )
        call check( nf90_inq_varid(ncid, t_name, t_id) )
        call check( nf90_inq_varid(ncid, u_name, u_id) )
        call check( nf90_inq_varid(ncid, v_name, v_id) )

        ! Get the lengths of the dimensions
        call check( nf90_inquire_dimension(ncid, x_dim_id, len=nx) )
        call check( nf90_inquire_dimension(ncid, y_dim_id, len=ny) )
        call check( nf90_inquire_dimension(ncid, t_dim_id, len=nt) )

        ! Allocate arrays for dimensions
        allocate(xc(nx))
        allocate(yc(ny))
        allocate(tc(nt))
        ! Allocate arrays for return data
        allocate(u(nx,ny,nt))
        allocate(v(nx,ny,nt))
        ! Allocate temporary arrays for reading data
        allocate(u_(nx,ny,nt))
        allocate(v_(nx,ny,nt))

        ! Start and count for reading variables
        start  = (/  1,  1, 1,  1 /) ! First element
        count  = (/ nx, ny, 1, nt /) ! Elements to read

        ! get x and y arrays
        call check( nf90_get_var(ncid, x_id, xc, (/ 1 /), (/ nx /)))
        call check( nf90_get_var(ncid, y_id, yc, (/ 1 /), (/ ny /)))
        call check( nf90_get_var(ncid, t_id, tc, (/ 1 /), (/ nt /)))

        ! getting values and handling FillValue, scale, offset
        ! read u (x-component of current)
        call check( nf90_get_var(ncid, u_id, u_, start, count))
        call check( nf90_get_att(ncid, u_id, 'scale_factor', scale_factor))
        call check( nf90_get_att(ncid, u_id, 'add_offset', add_offset))
        call check( nf90_get_att(ncid, u_id, '_FillValue', fill_value))
        where (u_ == fill_value)
            u = 0.0_WP ! Set land cells to 0 to avoid messing up interpolation
        elsewhere
            u = u_ * scale_factor + add_offset ! first scale, then offset
        end where

        ! read v (y-component of current)
        call check( nf90_get_var(ncid, v_id, v_, start, count))
        call check( nf90_get_att(ncid, v_id, 'scale_factor', scale_factor))
        call check( nf90_get_att(ncid, v_id, 'add_offset', add_offset))
        call check( nf90_get_att(ncid, v_id, '_FillValue', fill_value))
        where (v_ == fill_value)
            v = 0.0_WP ! Set land cells to 0 to avoid messing up interpolation
        elsewhere
            v = v_ * scale_factor + add_offset ! first scale, then offset
        end where

        ! Close file
        call check( nf90_close(ncid) )

        ! Clean up temporary arrays
        deallocate(u_)
        deallocate(v_)
    end subroutine

    subroutine check(status)
        integer, intent ( in) :: status
        if (status /= nf90_noerr) then
            print *, trim(nf90_strerror(status))
            stop "Stopped"
      end if
    end subroutine check
end module
