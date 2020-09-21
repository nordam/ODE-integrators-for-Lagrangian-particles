module output_module
    use hdf5
    use h5lt
    use parameters, only: SP, DP, WP
    implicit none

    contains


    subroutine write_to_hdf5(file_id, X, h_or_tol)
        implicit none
        ! Input variables
        ! ID of the file currently in use
        integer(hid_t),           intent(in)  :: file_id
        ! Particle positions
        real(WP), dimension(:,:), intent(in)  :: X
        ! Timestep or tolerance, used to give each group a name
        real(WP),                 intent(in)  :: h_or_tol

        ! Internal variables
        integer(SP)                           :: hdferr
        integer(hid_t)                        :: group_id
        character(len = 100)                  :: group_name

        ! Turn timestep or tolerance into a string
        write( group_name, * ) h_or_tol
        ! Creating a group for each timestep length or tolerance
        call open_or_create_h5group(file_id, trim(group_name), group_id, hdferr)

        ! Writing to group, using group_id as target instead of file_id
        call write_2d_array_double('X',  X, group_id, hdferr)

        ! closing group
        call h5gclose_f(group_id, hdferr)
    end subroutine write_to_hdf5

    subroutine create_hdf5_file(filename, file_id)
        ! Input variables
        character(len=*), intent(in) :: filename

        ! Output variables
        ! ID of the file currently in use
        integer(hid_t), intent(out)  :: file_id
        ! Error flag for the hdf5 library
        integer(SP) :: hdferr

        ! Initialize the hdf5 library
        call h5open_f(hdferr)
        ! Create output file
        call h5fcreate_f(trim(filename) // '.hdf5', h5f_acc_trunc_f, file_id, hdferr)
    end subroutine create_hdf5_file

    subroutine close_hdf5_file(file_id)
        ! Input variables
        ! ID of the file currently in use
        integer(hid_t), intent(in) :: file_id
        ! Output variables
        ! Error flag for the hdf5 library
        integer(SP) :: hdferr
        ! Close the datafile
        call h5fclose_f(file_id, hdferr)
        ! Close the hdf5 library
        call h5close_f(hdferr)
    end subroutine close_hdf5_file

    subroutine open_or_create_h5group(file_id, group_name, group_id, hdferr)
        ! Check if group with the given name exists. Create it if it doesn't,
        ! open it if it does.
        implicit none
        integer(hid_t), intent(in) :: file_id
        character(len=*), intent(in) :: group_name
        integer(SP), intent(inout) :: hdferr
        integer(hid_t), intent(out) :: group_id
        ! Local variables
        ! Variable for checking if a group exists or not
        logical :: group_exists

        call h5lexists_f(file_id, group_name, group_exists, hdferr)
        if (group_exists) then
            call h5gopen_f(file_id, group_name, group_id, hdferr)
        else
            call h5gcreate_f(file_id, group_name, group_id, hdferr)
        end if
    end subroutine open_or_create_h5group

    subroutine write_2d_array_double(array_name, array, file_id, hdferr)
        ! convenience function to dump arrays to hdf5
        ! TODO use interface and shit to make it handle single precision arrays
        implicit none

        ! input variables
        character(len=*), intent(in) :: array_name
        real(DP), intent(in), dimension(:, :) :: array
        ! ID of the file currently in use
        integer(hid_t), intent(in) :: file_id

        ! Output variables
        ! Error flag for the hdf5 library
        integer(SP), intent(out) :: hdferr

        ! Local variables
        ! ID of the dataset and dataspace in question
        integer(hid_t) :: dset_id, dspace_id, dcpl
        ! Rank of the R arrays
        integer(SP), parameter :: rank = 2
        ! Dimensions of the R arrays
        integer(hsize_t), dimension(rank) :: dims, chunk
        integer :: level

        ! Create a dataspace
        dims = shape(array)
        call h5screate_simple_f(rank, dims, dspace_id, hdferr)

        ! Create property list and set compression
        chunk = shape(array)
        level = 9
        call h5pcreate_f(H5P_DATASET_CREATE_F, dcpl, hdferr)
        call h5pset_deflate_f(dcpl, level, hdferr)
        call h5pset_chunk_f(dcpl, rank, chunk, hdferr)

        ! Create, write, and close datasets
        call h5dcreate_f(file_id, array_name, h5t_native_double, dspace_id, dset_id, hdferr, dcpl)
        call h5dwrite_f(dset_id, h5t_native_double, array, dims, hdferr)
        call h5dclose_f(dset_id, hdferr)
        call h5pclose_f(dcpl, hdferr)

        ! Close dataspace
        call h5sclose_f(dspace_id, hdferr)
    end subroutine write_2d_array_double

end module

! vim: ai ts=4 sts=4 et sw=4 tw=79 fenc=utf-8
