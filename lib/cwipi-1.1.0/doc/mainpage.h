/*!
  \mainpage CWIPI documentation

   \image html under-construction.jpg

   \section mainpages_install Installation

   \subsection basic_install Basic&nbsp;Installation

   mkdir build<br>
   cd build<br>
   cmake ..<br>
   make<br>
   make install<br>

   \subsection cmake_general_options CMake&nbsp;general&nbsp;options

cmake -D\<option1_name\>=\<option1_value\> ... -D\<option_name\>=\<optionn_value\>

Prefix :<br>
    CMAKE_INSTALL_PREFIX=\<prefix\>

Enable fortran interface :<br>
    CWP_ENABLE_Fortran=\<ON | OFF\> (default : OFF)

Enable python interface :<br>
    CWP_ENABLE_PYTHON_BINDINGS=<ON | OFF> (default : OFF)<br>
      If a simple autodetection fails, you can use these options to find Python :<br>
        PYTHON_LIBRARY=\<path\><br>
        PYTHON_INCLUDE_DIR=\<path\><br>
      Refere to FindPython in the CMake documentation for more informations.

Enable shared libraries :<br>
    CWP_ENABLE_SHARED=\<ON | OFF\> (default : ON)

Enable static libraries :<br>
    CWP_ENABLE_STATIC=\<ON | OFF\> (default : ON)

Enable the use of the library BLAS :<br>
    CWP_ENABLE_BLAS=\<ON | OFF\> (default : ON)

      If a simple autodetection fails, you can use these options to find BLAS :<br>
        BLAS_DIR=\<path\>      Where to find the base directory of blas<br>
        BLAS_INCDIR=\<path\>   Where to find the header files<br>
        BLAS_LIBDIR=\<path\>   Where to find the library files

   \subsection cmake_compiler_options CMake&nbsp;compiler&nbsp;options

  CC=\<C compiler\> CXX=\<CXX compiler\> FC=\<Fortran compiler\> cmake ...

  or use the following cmake options<br>
    CMAKE_C_COMPILER=\<C compiler\><br>
    CMAKE_CXX_COMPILER=\<CXX compiler\><br>
    CMAKE_Fortran_COMPILER=\<Fortran compiler\>


   \subsection cmake_mpi_options CMake&nbsp;MPI&nbsp;options

    MPI_C_COMPILER=\<C mpi wrapper\>
    MPI_CXX_COMPILER=\<CXX mpi wrapper\>
    MPI_Fortran_COMPILER=\<Fortran mpi wrapper\>

If a simple autodetection fails, you can use these options to find MPI :<br>
    MPI_\<lang\>_LIBRARIES<br>
    MPI_\<lang\>_INCLUDE_PATH

Refere to FindMPI in the CMake documentation for more informations.



   \section mainpage_examples Examples

*/
