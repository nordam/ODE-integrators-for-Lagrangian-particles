cmake_minimum_required(VERSION 3.5)
project (VariableTimestep)
enable_language (Fortran)


# make sure that the default is a RELEASE
if (NOT CMAKE_BUILD_TYPE)
  set (CMAKE_BUILD_TYPE RELEASE CACHE STRING
      "Choose the type of build, options are: None Debug Release."
      FORCE)
endif (NOT CMAKE_BUILD_TYPE)


set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/")
set (CMAKE_Fortran_FLAGS "-Wall")
set (CMAKE_Fortran_FLAGS_DEBUG "-O0 -g -pg -fimplicit-none -fbacktrace -fcheck=all")
set (CMAKE_Fortran_FLAGS_RELEASE "-O2 -mtune=native")


# Files with modules and stuff
set (fortran_modules
    "src/currentdata.f90"
    "src/experiment.f90"
    "src/input.f90"
    "src/integrator.f90"
    "src/interpolator.f90"
    "src/output.f90"
    "src/parameters.f90"
    )


# Using hdf5 to store output
find_package (HDF5 COMPONENTS HL Fortran_HL Fortran)
include_directories (${HDF5_INCLUDE_DIRS})
message ("++ HDF5_LIBRARIES:       ${HDF5_LIBRARIES}")
message ("++ HDF5_INCLUDE_DIRS:    ${HDF5_INCLUDE_DIRS}")

# Using netCDF to read current data from files
file (GLOB  NETCDF_LIBRARIES  "${NETCDF_LIBRARY_PATH}/libnetcdf*.so")
include_directories (${NETCDF_INCLUDE_DIRS})
message ("++ NETCDF_LIBRARIES:     ${NETCDF_LIBRARIES}")
message ("++ NETCDF_INCLUDE_DIRS:  ${NETCDF_INCLUDE_DIRS}")

# Need to link to the bspline-fortran library
# https://github.com/jacobwilliams/bspline-fortran
file (GLOB BSPLINE_LIBRARIES "${BSPLINE_LIBRARY_PATH}/*.a")
include_directories (${BSPLINE_LIBRARY_PATH})
message ("++ BSPLINE_LIBRARY_PATH: ${BSPLINE_LIBRARY_PATH}")
message ("++ BSPLINE_LIBRARIES:    ${BSPLINE_LIBRARIES}")



#### Simple test program
add_executable (testprogram "src/test.f90" ${fortran_modules})
target_link_libraries (testprogram ${BSPLINE_LIBRARIES})
target_link_libraries (testprogram ${HDF5_LIBRARIES})
target_link_libraries (testprogram ${NETCDF_LIBRARIES})

#### Run simulations with current data
add_executable (run_norkyst "src/run_norkyst.f90" ${fortran_modules})
target_link_libraries (run_norkyst ${BSPLINE_LIBRARIES})
target_link_libraries (run_norkyst ${HDF5_LIBRARIES})
target_link_libraries (run_norkyst ${NETCDF_LIBRARIES})

#### Run simulations with current data
add_executable (run_nordic "src/run_nordic.f90" ${fortran_modules})
target_link_libraries (run_nordic ${BSPLINE_LIBRARIES})
target_link_libraries (run_nordic ${HDF5_LIBRARIES})
target_link_libraries (run_nordic ${NETCDF_LIBRARIES})

#### Run simulations with current data
add_executable (run_arctic "src/run_arctic.f90" ${fortran_modules})
target_link_libraries (run_arctic ${BSPLINE_LIBRARIES})
target_link_libraries (run_arctic ${HDF5_LIBRARIES})
target_link_libraries (run_arctic ${NETCDF_LIBRARIES})




#### Reference solutions ####


#### Run simulations with current data
add_executable (run_norkyst_ref "src/run_norkyst_ref.f90" ${fortran_modules})
target_link_libraries (run_norkyst_ref ${BSPLINE_LIBRARIES})
target_link_libraries (run_norkyst_ref ${HDF5_LIBRARIES})
target_link_libraries (run_norkyst_ref ${NETCDF_LIBRARIES})

#### Run simulations with current data
add_executable (run_nordic_ref "src/run_nordic_ref.f90" ${fortran_modules})
target_link_libraries (run_nordic_ref ${BSPLINE_LIBRARIES})
target_link_libraries (run_nordic_ref ${HDF5_LIBRARIES})
target_link_libraries (run_nordic_ref ${NETCDF_LIBRARIES})

#### Run simulations with current data
add_executable (run_arctic_ref "src/run_arctic_ref.f90" ${fortran_modules})
target_link_libraries (run_arctic_ref ${BSPLINE_LIBRARIES})
target_link_libraries (run_arctic_ref ${HDF5_LIBRARIES})
target_link_libraries (run_arctic_ref ${NETCDF_LIBRARIES})
