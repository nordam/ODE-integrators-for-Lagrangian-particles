# ODE-integrators-for-Lagrangian-particles

This repository contains the code used to run simulations and analyse results for [Nordam & Duran (2020)](https://gmd.copernicus.org/preprints/gmd-2020-154/). The simulation code is implemented in Fortran, using three external libraries (netCDF, hdf5 and bspline-fortran). The analysis and plotting of results is done in jupyter notebooks.

## Structure of the repo

### build

This directory is where you build your simulation code. It initially contains only a .gitignore file ignoring everything but itself, ensuring that files in this folder do not clutter the git status.

### data

This directory contain the ocean current datasets used to run the simulations, as well as the initial positions used for the 10000 particles. For reference, it also contains a map showing the extent of the datasets. See also data/README.md

### notebooks

This directory contains jupyter notebooks that are used to analyse and plot the results of running the simulations.

### src

This directory contains the fortran source code to run the simulations.

## Build instructions for fortran code

* Download and build the bspline-fortran library from https://github.com/jacobwilliams/bspline-fortran
* Make sure you have
  * cmake version 3.5 or higher
  * hdf5, with fortran interface
  * netcdf, with fortran interface
  * a recent version of gfortran
* Go into the build directory (or create your own, somewhere)
* We use cmake to configure the build, and we need to pass information about the location of the different libraries. The following works for me on xubuntu 18.04 with netcdf installed from apt-get. Adjust the paths to fit your system.
* `cmake .. -DBSPLINE_LIBRARY_PATH=PATH_TO_BSPLINE_LIBS -DNETCF_LIBRARY_PATH=PATH_TO_NETCDF_LIBS -DNETCDF_INCLUDE_DIRS=PATH_TO_NETCF_INCLUDES -DCMAKE_BUILD_TYPE=RELEASE`
  * `PATH_TO_BSPLINE_LIBS` should be the folder where `bspline_module.o` and `bspline_module.mod` etc. are found
  * `PATH_TO_NETCDF_LIBS` should be the folder where `libnetcdf.so` and `libnetcdff.so` found
  * `PATH_TO_NETCDF_INCLUDES` should be the folder where `netcdf.h` and `netcdf.mod` etc. are found
* Then running `make` should build the project. If that doesn't work, try running `make VERBOSE=1`, copy the last command, where the error occured, and see if all the paths are correct.

## Instructions for using the jupyter notebooks

* First, build the fortran code and run the simulations.
