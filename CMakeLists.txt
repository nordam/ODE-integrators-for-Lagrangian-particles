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
set (CMAKE_Fortran_FLAGS_RELEASE "-funroll-loops -O3 -mtune=native")


# Files with modules and stuff
set (fortran_modules
    "src/parameters.f90"
    "src/interpolator.f90"
    "src/integrators.f90"
    )


# Need to link to the bspline-fortran library
# https://github.com/jacobwilliams/bspline-fortran
message ("++ Bspline library path: ${BSPLINE_LIBRARY_PATH}")
file (GLOB BSPLINE_LIBRARIES "${BSPLINE_LIBRARY_PATH}/*.so")
include_directories (${BSPLINE_LIBRARY_PATH})


#### Simple test program
add_executable (testprogram "src/test.f90" ${fortran_modules})
target_link_libraries (testprogram ${BSPLINE_LIBRARIES})
