# ODE-integrators-for-Lagrangian-particles

This repository contains the code used to run simulations and analyse results for [Nordam & Duran (2020)](https://gmd.copernicus.org/preprints/gmd-2020-154/). The simulation code is implemented in Fortran, using three external libraries (netCDF, hdf5 and bspline-fortran). The analysis and plotting of results is done in jupyter notebooks.


## Build instructions for fortran code

* Download and build the bspline-fortran library from https://github.com/jacobwilliams/bspline-fortran
* Make sure you have
  * cmake version 3.5 or higher
  * hdf5, with fortran interface
  * netcdf, with fortran interface
  * a recent version of gfortran
* Create a separate build directory somewhere
* We use cmake to configure the build, and we need to pass information about the location of the different libraries. Currently, the following seems to work (but could definitely be improved upon)
* `PATH$PATH:PATH_TO_HDF5_LIBS:PATH_TO_HDF5_INCLUDES cmake PATH_TO_CMakeLists.txt -DBSPLINE_LIBRARY_PATH=PATH_TO_BSPLINE_LIBS -DNETCF_LIBRARY_PATH=PATH_TO_NETCDF_LIBS -DNETCDF_INCLUDE_DIRS=PATH_TO_NETCF_INCLUDES`
  * `PATH_TO_HDF5_LIBS` should be the folder where `libhdf5_fortran.so` etc. are found
  * `PATH_TO_HDF5_INCLUDES` should be the folder where the hdf5 header files are found, on my system the mod-files are in a subfolder of that folder, called `shared`
  * `PATH_TO_BSPLINE_LIBS` should be the folder where `bspline_module.o` and `bspline_module.mod` etc. are found
  * `PATH_TO_NETCDF_LIBS` should be the folder where `libnetcdf.so` and `libnetcdff.so` found
  * `PATH_TO_NETCDF_INCLUDES` should be the folder where `netcdf.h` and `netcdf.mod` etc. are found
* Then running `make` should build the project. If that doesn't work, try running `make VERBOSE=1`, copy the last command, where the error occured, and see if all the paths are correct.
