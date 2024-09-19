.. _installation:

Installation
############

Dependencies
============

General dependencies for building **ParaDiGM** are:

  * a C++ compiler
  * `CMake <https://cmake.org/>`_ (version 3.16 or higher)
  * an MPI distribution

Basic Installation
==================

Follow these steps to build **ParaDiGM** from the sources:

.. code-block:: sh

  mkdir build
  cd build
  cmake ..
  make
  make install


If installation fails, use the following CMake options.


CMake general options
=====================

.. code-block:: sh

    cmake -D<option1_name>=<option1_value> ... -D<option2_name>=<option2_value>

**Installation prefix**

.. code-block:: sh

  CMAKE_INSTALL_PREFIX=<prefix>

.. _enable_fortran_interface:

**Enable Fortran interface**

.. code-block:: sh

  PDM_ENABLE_Fortran=<ON | OFF> (default : OFF)

.. _enable_python_interface:

**Enable Python interface**

.. code-block:: sh

  PDM_ENABLE_PYTHON_BINDINGS=<ON | OFF> (default : OFF)

If a simple autodetection fails, you can use these options to find Python :

.. code-block:: sh

    PYTHON_LIBRARY=<path>
    PYTHON_INCLUDE_DIR=<path>

Refer to `FindPython <https://cmake.org/cmake/help/latest/module/FindPython.html>`_ in the CMake documentation for more information.

**Build shared library**

.. code-block:: sh

  PDM_ENABLE_SHARED=<ON | OFF> (default : ON)

**Build static library**

.. code-block:: sh

  PDM_ENABLE_STATIC=<ON | OFF> (default : ON)

.. _parmetis: https://github.com/KarypisLab/ParMETIS
.. |parmetis| replace:: **ParMETIS**

**Enable the use of** |parmetis|_ **(parallel graph partitioning)**

.. code-block:: sh

    PDM_ENABLE_PARMETIS=<ON | OFF> (default : ON)

If a simple autodetection fails, you can use these options to find ParMETIS :

.. code-block:: sh

    PARMETIS_DIR=<path>

| To link shared libraries, ParMETIS must be compiled with the ``-fPIC`` flag.
| CMake looks for

  * ``parmetis.h`` and ``metis.h`` includes
  * ``parmetis`` and ``metis`` libraries


.. _ptscotch: https://gitlab.inria.fr/scotch/scotch
.. |ptscotch| replace:: **PT-Scotch**

**Enable the use of** |ptscotch|_ **(parallel graph partitioning)**

.. code-block:: sh

    PDM_ENABLE_PTSCOTCH=<ON | OFF> (default : ON)

If a simple autodetection fails, you can use these options to find PT-Scotch :

.. code-block:: sh

    PTSCOTCH_DIR=<path>

| To link shared libraries, PT-Scotch must be compiled with the ``-fPIC`` flag.
| CMake looks for

  * ``ptscotch.h`` include file
  * ``scotch``, ``scotcherr``, ``ptscotch``, ``ptscotcherr`` libraries


.. _blas: https://www.netlib.org/blas/
.. |blas| replace:: **BLAS**

.. _lapack: https://www.netlib.org/lapack/
.. |lapack| replace:: **LAPACK**

**Enable the use of** |blas|_ **/** |lapack|_ **(linear algebra)**

.. code-block:: sh

    PDM_ENABLE_BLASLAPACK=<ON | OFF> (default : OFF)

**Enable long global IDs**

.. code-block:: sh

    PDM_ENABLE_LONG_G_NUM=<ON | OFF> (default : ON)

* ``ON``  : ``PDM_g_num_t`` type is ``long int``
* ``OFF`` : ``PDM_g_num_t`` type is ``int``

**Enable documentation compilation**

.. code-block:: sh

    PDM_ENABLE_DOC=<ON | OFF> (default : OFF)

Once built, the documentation can be found in ``build/doc/sphinx/html`` and launch ``index.html`` file


Compiler choice
===============

.. code-block:: sh

    CC=<C compiler> CXX=<CXX compiler> FC=<Fortran compiler> cmake ...

or use the following CMake options

.. code-block:: sh

    CMAKE_C_COMPILER=<C compiler>
    CMAKE_CXX_COMPILER=<CXX compiler>
    CMAKE_Fortran_COMPILER=<Fortran compiler>


CMake MPI options
=================

.. code-block:: sh

    MPI_C_COMPILER=<C MPI wrapper>
    MPI_CXX_COMPILER=<CXX MPI wrapper>
    MPI_Fortran_COMPILER=<Fortran MPI wrapper>

If a simple autodetection fails, you can use these options to find MPI :

.. code-block:: sh

    MPI_<language>_LIBRARIES
    MPI_<language>_INCLUDE_PATH

Refer to `FindMPI <https://cmake.org/cmake/help/latest/module/FindMPI.html>`_ in the CMake documentation for more information.
