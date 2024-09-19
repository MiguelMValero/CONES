# CWIPI #

**CWIPI** (Coupling With Interpolation Parallel Interface) is a parallel coupling library under LGPL, with interfaces in C, Fortran and Python.

## Documentation  ##

User documentation is deployed on the Gitlab pages server: https://numerics.gitlab-pages.onera.net/coupling/cwipi/cwipi-1.1.0/index.html

## Build and install ##

### Dependencies

General dependencies for building **CWIPI** are:
- a C++ compiler
- [CMake](https://cmake.org/) (version 3.16 or higher)
- an MPI distribution

### Basic Installation

Follow these steps to build **CWIPI** from the sources:

1. `git clone git@gitlab.onera.net:numerics/coupling/cwipi.git` (for ONERA users only)
1. `cd cwipi`
1. `git submodule update --init` (needed for dependencies)
1. `mkdir build`
1. `cd build`
1. `cmake ..`
1. `make`
1. `make install`
1. `./cwp_run` (if you want to run the test cases)


### CMake general options
    cmake -D<option1_name>=<option1_value> ... -D<option2_name>=<option2_value>

#### Installation prefix
    CMAKE_INSTALL_PREFIX=<prefix>

#### Enable Fortran interface
    CWP_ENABLE_Fortran=<ON | OFF> (default : OFF)

#### Enable Python interface
    CWP_ENABLE_PYTHON_BINDINGS=<ON | OFF> (default : OFF)

If a simple autodetection fails, you can use these options to find Python :

    PYTHON_LIBRARY=<path>
    PYTHON_INCLUDE_DIR=<path>

Refer to [FindPython](https://cmake.org/cmake/help/latest/module/FindPython.html) in the CMake documentation for more information.

#### Build shared library
    CWP_ENABLE_SHARED=<ON | OFF> (default : ON)

#### Build static library
    CWP_ENABLE_STATIC=<ON | OFF> (default : ON)

#### Enable the use of [BLAS](https://www.netlib.org/blas/) (linear algebra)
    CWP_ENABLE_BLAS=<ON | OFF> (default : ON)

If a simple autodetection fails, you can use these options to find BLAS :

    BLAS_DIR=<path>      Where to find the base directory of BLAS
    BLAS_INCDIR=<path>   Where to find the header files
    BLAS_LIBDIR=<path>   Where to find the library files

To force the use of a list of libraries, use :

    DBLAS_LIBRARIES="<lib_1> ... <lib_n>"


#### Enable client-server mode
    CWP_ENABLE_CLIENT_SERVER=<ON | OFF> (default : OFF)

#### Enable documentation mode
    CWP_ENABLE_DOCUMENTATION=<ON | OFF> (default : OFF)

Once built, the documentation can be found in `build/doc/sphinx/html` and launch `index.html` file

### Compiler choice

    CC=<C compiler> CXX=<CXX compiler> FC=<Fortran compiler> cmake ...

or use the following CMake options

    CMAKE_C_COMPILER=<C compiler>
    CMAKE_CXX_COMPILER=<CXX compiler>
    CMAKE_Fortran_COMPILER=<Fortran compiler>

### CMake MPI options

    MPI_C_COMPILER=<C MPI wrapper>
    MPI_CXX_COMPILER=<CXX MPI wrapper>
    MPI_Fortran_COMPILER=<Fortran MPI wrapper>

If a simple autodetection fails, you can use these options to find MPI :

    MPI_<language>_LIBRARIES
    MPI_<language>_INCLUDE_PATH

Refer to [FindMPI](https://cmake.org/cmake/help/latest/module/FindMPI.html) in the CMake documentation for more information.

## Issues ##

Issues can be reported directly in the [Issues](https://gitlab.onera.net/numerics/coupling/cwipi/-/issues) section.


## License ##

**CWIPI** is available under the LGPL3 license (https://www.gnu.org/licenses/lgpl-3.0.fr.html).


## Copyright ##

Copyright 2023, ONERA The French Aerospace Lab
