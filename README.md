# ODE-integrators-for-Lagrangian-particles

This repository contains the code used to run simulations and analyse results for [Nordam & Duran (2020)](https://gmd.copernicus.org/preprints/gmd-2020-154/). The simulation code is implemented in Fortran, using three external libraries (netCDF, hdf5 and bspline-fortran). The analysis and plotting of results is done in jupyter notebooks.

## Structure of the repo

### build

This directory is where you build your simulation code. It initially contains only a .gitignore file ignoring everything but itself, ensuring that files in this folder do not clutter the git status.

### data

This directory contain the ocean current datasets used to run the simulations, as well as the initial positions used for the 10000 particles. For reference, it also contains a map showing the extent of the datasets. See also data/README.md

### notebooks

This directory contains jupyter notebooks that are used to analyse and plot the results of running the simulations.

### results

This directory is where the simulation results will be stared. It initially contains only a .gitignore file ignoring everything but itself, ensuring that files in this folder do not clutter the git status.

### scripts

This directory contains a tiny shell script which will run all six executables in parallel (and stop them all if the script is interrupted with ctrl-c).

### src

This directory contains the fortran source code to run the simulations:
* `currentdata.f90` uses the netCDF library to read current data from the .nc files in the data folder.
* `experiment.f90` contains some convenience functions to run a series of simulations.
* `input.f90` contains functions to read initial particle positions from the .txt files in the data folder.
* `integrator.f90` contains implementations of the integrators described in [Nordam & Duran (2020)](https://gmd.copernicus.org/preprints/gmd-2020-154/).
* `interpolator.f90` defines a derived type used to conveniently hold interpolator objects from the bspline-fortran library.
* `output.f90` contains functions to create and write hdf5 files with the particle positions after transport.
* `parameters.f90` defines some parameters. Edit this file to changes paths, and expand range of timesteps or tolerances.
* `run_*.f90` are three almost identical files, one for each dataset, to run a set of simulations.
* `run_*_ref.f90` are three almost identical files, one for each dataset, to obtain reference solutions.

## Build instructions for fortran code

* Download and build the bspline-fortran library from https://github.com/jacobwilliams/bspline-fortran
* Make sure you have
  * cmake version 3.5 or higher
  * hdf5, with fortran interface
  * netcdf, with fortran interface
  * a recent version of gfortran
* Go into the build directory (or create your own, somewhere)
* We use cmake to configure the build, and we need to pass information about the location of the different libraries. The following works for me on xubuntu 18.04 with netcdf installed from apt-get. Adjust the paths to fit your system.
* `cmake .. -DBSPLINE_LIBRARY_PATH=bspline_path -DNETCF_LIBRARY_PATH=netcdf_libs -DNETCDF_INCLUDE_DIRS=netcdf_include -DCMAKE_BUILD_TYPE=RELEASE`
  * `bspline_path` should be the folder where `bspline_module.o` and `bspline_module.mod` etc. are found. This depends on where you built bspline-fortran.
  * `netcdf_libs` should be the folder where `libnetcdf.so` and `libnetcdff.so` found. On my system, running xubuntu 18.04 with `libnetcdf-dev` and `libnetcdff-dev` installed from apt-get, these files are found in `/usr/lib/x86_64-linux-gnu`.
  * `netcdf_include` should be the folder where `netcdf.h` and `netcdf.mod` etc. are found. On my system, running xubuntu 18.04 with `libnetcdf-dev` and `libnetcdff-dev` installed from apt-get, these files are found in `/usr/include`.
* Then running `make` should build the project. If that doesn't work, try running `make VERBOSE=1`, copy the last command, where the error occured, and see if all the paths are correct.
* Running `make` will build six executables, two for each of the three different resolution datasets. Running, e.g., `run_norkyst` will run simulations with the NorKyst800m dataset as input, scanning through the 8 different timesteps and 11 different tolerances investigated in the paper. Running `run_norkyst_ref` will run simulations with a selection of shorter timesteps and stricter tolerances, in order to obtain reference solutions (see Appendix A in [Nordam & Duran (2020)](https://gmd.copernicus.org/preprints/gmd-2020-154/) for details).
* Note that running the full set of timesteps and tolerances used in [Nordam & Duran (2020)](https://gmd.copernicus.org/preprints/gmd-2020-154/) takes a good while (20-30 hours on a 3.3 GHz Intel Xeon CPU). Therefore, the default here is to run a reduced set of simulations. Edit `src/parameters.f90` to expand the range of timesteps and tolerances. The simulations all run on a single core and don't require very much memory, so some trivial parallelisation can be achieved by running all six at the same time. See the convenince script `scripts/run_all.sh`.
* Note also that the executables expect to find input data in `../data/`, and expect to be able to write their output to `../results/`. If you want a different setup, adjust paths and edit the code in `src/parameters.f90` as required.

## Instructions for using the jupyter notebooks

* First, build the fortran code and run the simulations.
* Then, stepping through the cells in order should reproduce Figures 5, B1 and B2 from [Nordam & Duran (2020)](https://gmd.copernicus.org/preprints/gmd-2020-154/), except with a somewhat smaller number of datapoints as mentioned above.
