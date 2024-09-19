#==============================================================================
#
#             Top level CMakeLists.txt file for ParaDiGM
#
#==============================================================================

cmake_minimum_required(VERSION 3.16)
cmake_policy(VERSION 3.16)


#------------------------------------------------------------------------------
# Internal CMake module (shared with ParaDiGMA or CWIPI)
#------------------------------------------------------------------------------

set(PDM_CMAKE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/cmake")
set(CMAKE_MODULE_PATH "${PDM_CMAKE_DIR}/modules")

include(${CMAKE_MODULE_PATH}/version.cmake)
project(PDM VERSION ${PDM_DEF_VERSION_MAJOR} LANGUAGES C CXX)

set(PDM_VERSION "${PDM_VERSION_MAJOR}.${PDM_VERSION_MINOR}.${PDM_VERSION_PATCH}")

set(CTEST_TEST_TIMEOUT 120)

# ----------------------------------------------------------------------
find_program( MEMORYCHECK_COMMAND valgrind )
set( MEMORYCHECK_COMMAND_OPTIONS "--trace-children=yes --leak-check=full" )
# ----------------------------------------------------------------------

include(CTest)
enable_testing ()

# ----------------------------------------------------------------------
# RPATH protection


#SET(CMAKE_SKIP_BUILD_RPATH  TRUE)
#set(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE)
#set(CMAKE_INSTALL_RPATH_USE_LINK_PATH FALSE)


set(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE)
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")

#------------------------------------------------------------------------------
# Fake dlclose path
#------------------------------------------------------------------------------

set(PDM_BINARY_DLCLOSE_DIR ${CMAKE_BINARY_DIR}/script/asan/fake_dlclose)
set(PDM_SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR})


#------------------------------------------------------------------------------
# User option Definition
#------------------------------------------------------------------------------

# Multipart in Part
# -----------------
set(PDM_USE_MULTIPART 0)

# Extensions
# ----------

option(PDM_ENABLE_EXTENSION_PDMA "Enable advanced functionalities extension" OFF)

# OpenMP
# ------

option(PDM_ENABLE_OPENMP "Enable OpenMP" OFF)

# CUDA
# ------

option(PDM_ENABLE_CUDA "Enable CUDA" OFF)
if (PDM_ENABLE_CUDA)
  enable_language(CUDA)
endif ()

# ANISO_AGGLO (only if pdma)
# ------
if (PDM_ENABLE_EXTENSION_PDMA)
  option(PDM_ENABLE_ANISO_AGGLO "Enable Anisotropic agglomeration for coarse mesh function" ON)

  if (PDM_ENABLE_ANISO_AGGLO)
    if (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
      if(CMAKE_CXX_COMPILER_VERSION VERSION_LESS 4.8)
        message(FATAL_ERROR "Insufficient g++ version for Anisotropic agglomeration : set PDM_ENABLE_ANISO_AGGLO to OFF and re-run")
      endif()
    elseif (CMAKE_CXX_COMPILER_ID STREQUAL "Intel")
      if(CMAKE_CXX_COMPILER_VERSION VERSION_LESS 15)
        message(FATAL_ERROR "Insufficient icpc version for Anisotropic agglomeration : set PDM_ENABLE_ANISO_AGGLO to OFF and re-run")
      endif()
    elseif (CMAKE_CXX_COMPILER_ID STREQUAL "AppleClang")
      if(CMAKE_CXX_COMPILER_VERSION VERSION_LESS 3.2)
        message(FATAL_ERROR "Insufficient clang version for Anisotropic agglomeration : set PDM_ENABLE_ANISO_AGGLO to OFF and re-run")
      endif()
    elseif (CMAKE_CXX_COMPILER_ID MATCHES "XL")
      if(CMAKE_CXX_COMPILER_VERSION VERSION_LESS 13.1.3)
        message(FATAL_ERROR "Insufficient xlc++ version for Anisotropic agglomeration : set PDM_ENABLE_ANISO_AGGLO to OFF and re-run")
      endif()
    endif ()
    set (PDM_HAVE_ANISO_AGGLO 1)
  endif ()
endif()


# Enable long for absolute number
# -------------------------------

option(PDM_ENABLE_LONG_G_NUM "Enable long for global number" ON)

# Pt-scotch
# ---------

option(PDM_ENABLE_PTSCOTCH "Compile with support for SCOTCH." ON)

# ParMetis
# --------

option(PDM_ENABLE_PARMETIS "Compile with support for SCOTCH." ON)

# Shared libraries
# ----------------

option(PDM_ENABLE_SHARED "Build Shared Libraries" OFF)

# Static libraries
# ----------------

option(PDM_ENABLE_STATIC "Build Static Libraries" ON)

# Check MPI Wrapper
# -----------------

option(PDM_ENABLE_MPI_CHECK "Check MPI Wrapper" ON)

# Python Bindings
# ---------------

option(PDM_ENABLE_PYTHON_BINDINGS "Build Python module" OFF)

# BLAS/LAPACK
# -----------

option(PDM_ENABLE_BLASLAPACK "Build with BLAS/LAPACK" OFF)

# Fortran interface
# -----------------

option(PDM_ENABLE_Fortran "Build Fortran Interface" ON)

if (PDM_ENABLE_Fortran)
  enable_language(Fortran)
endif()

if (PDM_ENABLE_STATIC AND PDM_ENABLE_Fortran)
  set(PDM_ENABLE_STATIC_Fortran_INTERFACE ON)
else ()
  set(PDM_ENABLE_STATIC_Fortran_INTERFACE OFF)
endif ()

if (PDM_ENABLE_SHARED AND PDM_ENABLE_Fortran)
  set(PDM_ENABLE_SHARED_Fortran_INTERFACE ON)
else ()
  set(PDM_ENABLE_SHARED_Fortran_INTERFACE OFF)
endif ()

option(PDM_ENABLE_Fortran_MPI_MODULE "Use Fortran MPI module" ON)

# Unit test mpi with doctect
# ------

option(PDM_ENABLE_UNIT_TEST "Enable unit test " OFF)


# Documentation
# -------------

option(PDM_ENABLE_DOC "Enable Sphinx documentation" OFF)


# Training Notebooks
# ------------------

option(PDM_ENABLE_TRAINING "Enable training Jupyter Notebooks" OFF)

#------------------------------------------------------------------------------
# Default build : Release
#------------------------------------------------------------------------------

if ((PDM_ENABLE_PYTHON_BINDINGS OR PDM_ENABLE_SHARED_Fortran_INTERFACE) AND (NOT PDM_ENABLE_SHARED))
  set (PDM_ENABLE_SHARED "ON" CACHE STRING "Build Shared Libraries" FORCE )
endif()

if ((NOT PDM_ENABLE_STATIC) AND (NOT PDM_ENABLE_SHARED))
  message (FATAL_ERROR "No enabled library. Please, enable shared or static library")
endif()

if (NOT CMAKE_BUILD_TYPE)
 set(CMAKE_BUILD_TYPE "Release" CACHE STRING  "Choose the type of build, options are: None(CMAKE_CXX_FLAGS or CMAKE_C_FLAGS used) Debug Release RelWithDebInfo MinSizeRel" FORCE)
endif()

# ADD "-DDEBUG_CLASSE" in Debug
#------------------------------

if (CMAKE_BUILD_TYPE STREQUAL "Debug")
  add_definitions (-DDEBUG_CLASSE)
endif()

# Add "-DHAVE_GETRUSAGE" if function getrusage and gettimeofday exist
#--------------------------------------------------------------------

include(CheckFunctionExists)
CHECK_FUNCTION_EXISTS(getrusage PDM_HAVE_GETRUSAGE)
CHECK_FUNCTION_EXISTS(gettimeofday PDM_HAVE_GETTIMEOFDAY)

#------------------------------------------------------------------------------
# Default compiler flags
#------------------------------------------------------------------------------


include(${CMAKE_MODULE_PATH}/default_flags.cmake)

#------------------------------------------------------------------------------
# Libraries to add to link
#------------------------------------------------------------------------------

set(COMMON_LIBRARIES)

#------------------------------------------------------------------------------
# Check dependecies
#------------------------------------------------------------------------------


# Check MPI (TODO : check MPI wrapper function check )
# -------------

if (PDM_ENABLE_Fortran)
  find_package(MPI REQUIRED COMPONENTS C Fortran)
else()
  find_package(MPI REQUIRED COMPONENTS C)
endif()

if (PDM_ENABLE_Fortran_MPI_MODULE)
  set(CMAKE_Fortran_PREPROCESS ON)
  set (PDM_HAVE_FORTRAN_MPI_MODULE 1)
endif()

mark_as_advanced (MPI_EXTRA_LIBRARY MPI_LIBRARY)

# ---------------------
# Check MPI Wrapper
# ---------------------

if (PDM_ENABLE_MPI_CHECK)
   include(${CMAKE_MODULE_PATH}/CheckMPI.cmake)
endif()

# Check SCOTCH an ParMETIS
# -------------------------

if (PDM_ENABLE_PARMETIS)
   find_package(ParMETIS 4.0.3)
   if (PARMETIS_FOUND)
      set(PDM_HAVE_PARMETIS 1)
   endif()
endif()

if (PDM_ENABLE_PTSCOTCH)
   find_package(PTSCOTCH 6.0.0)
   if (PTSCOTCH_FOUND)
      set(PDM_HAVE_PTSCOTCH 1)
   endif()
endif()


if (PDM_ENABLE_BLASLAPACK)

  if (CMAKE_C_COMPILER_ID STREQUAL "Intel")
    set (BLA_VENDOR "Intel10_64lp")
    find_package(BLAS)
    find_package(LAPACK)
    if (BLAS_FOUND AND LAPACK_FOUND)
      set(PDM_HAVE_MKL 1)
    endif()
    unset (BLA_VENDOR)
  endif()

  if (NOT BLAS_FOUND)
    find_package(BLAS)
  endif()

  if (NOT LAPACK_FOUND)
    find_package(LAPACK)
  endif()

  if (BLAS_FOUND)
    set(PDM_HAVE_BLAS 1)
  endif()

  if (LAPACK_FOUND)
    set(PDM_HAVE_LAPACK 1)
  endif()

endif()

# find_package(mmg CONFIG)

# Check OpenMP (TODO : check MPI wrapper function check )
# -------------------------------------------------------

if (PDM_ENABLE_OPENMP)
  find_package(OpenMP)
  if (OPENMP_FOUND)
    set(PDM_HAVE_OPENMP 1)
  endif()
endif()

# Python bindings
# ---------------

if (PDM_ENABLE_PYTHON_BINDINGS)
  # NB: cython module must use mpicxx to avoid "mpi.h" not found.
  find_package(MPI REQUIRED COMPONENTS CXX)

  find_package(Python REQUIRED COMPONENTS Interpreter Development NumPy)
  find_package(Mpi4Py REQUIRED)

  include(UseCython)

endif()

#------------------------------------------------------------------------------
# Store variable in configure file
#------------------------------------------------------------------------------

if (PDM_ENABLE_LONG_G_NUM)
  set(PDM_LONG_G_NUM 1)
endif ()

configure_file(pdm_config.h.in "${CMAKE_CURRENT_BINARY_DIR}/pdm_config.h")

if (PDM_ENABLE_STATIC_Fortran_INTERFACE OR PDM_ENABLE_SHARED_Fortran_INTERFACE)
  if ("${CMAKE_CURRENT_BINARY_DIR}/pdm_config.h" IS_NEWER_THAN  "${CMAKE_CURRENT_BINARY_DIR}/pdm_configf.h")
    include (FortranizeFile)
    fortranize ("${CMAKE_CURRENT_BINARY_DIR}/pdm_config.h")
  endif ()
endif ()

#------------------------------------------------------------------------------
# Doc
#------------------------------------------------------------------------------


find_package(Doxygen)

#------------------------------------------------------------------------------
# Print summary of found and not found optional packages
#------------------------------------------------------------------------------


include(FeatureSummary)
FEATURE_SUMMARY(WHAT ALL)

#------------------------------------------------------------------------------
# Installation de PDM
#------------------------------------------------------------------------------

# Append the library version information to the library target properties

if (PDM_WITH_LIBRARY_VERSION)
  string(REPLACE "+" "" PDM_LIBRARY_VERSION ${PDM_VERSION})
  # This setting of SOVERSION assumes that any API change
  # will increment either the minor or major version number.
  set(PDM_LIBRARY_PROPERTIES ${PDM_LIBRARY_PROPERTIES}
    VERSION ${PDM_LIBRARY_VERSION}
    SOVERSION ${PDM_VERSION_MAJOR}.${PDM_VERSION_MINOR}
     )
endif()

# Set PDM install sub-directories
set(PDM_BIN_DIR "bin")
set(PDM_LIB_DIR "lib")
set(PDM_INCLUDE_DIR "include")
set(PDM_PKGCONFIG_DIR "lib/pkgconfig")
set(PDM_SHARE_DIR "share/cedre")
set(PDM_MAN_DIR "share/man")
set(PDM_DOC_DIR "${PDM_SHARE_DIR}/doc")
set(PDM_ETC_DIR "etc")

#------------------------------------------------------------------------------
# Ajout des sources de tous les elements du projet
# En fonction des top_targets a installer
#------------------------------------------------------------------------------

set(PDM_LIST_TEST "")
set(PDM_LIST_NRANK "")

add_subdirectory(script)

add_subdirectory(src)

option(PDM_ENABLE_TESTS "Enable tests" ON)

add_subdirectory(extensions)

if (PDM_ENABLE_PYTHON_BINDINGS)
  add_subdirectory(Cython)
endif()

if (PDM_ENABLE_TESTS)
  add_subdirectory(test)
endif()

add_subdirectory(doc)

if (PDM_ENABLE_TRAINING)
  add_subdirectory(training)
endif()

#------------------------------------------------------------------------------
# Create PDM_run file
#------------------------------------------------------------------------------

string(REPLACE ";" " " PDM_LIST_TEST "${PDM_LIST_TEST}")

string(REPLACE ";" " " PDM_LIST_NRANK "${PDM_LIST_NRANK}")

set (PDM_RUN_PRELOAD "")
set (PDM_RUN_LSAN "")
set (PDM_RUN_PYTHONPATH "")

if (CMAKE_BUILD_TYPE STREQUAL "Sanitize")
  execute_process(COMMAND gcc -print-file-name=libasan.so OUTPUT_VARIABLE PRELOAD_ASAN OUTPUT_STRIP_TRAILING_WHITESPACE)
  set(PDM_RUN_PRELOAD "LD_PRELOAD ${PRELOAD_ASAN}:${PROJECT_BINARY_DIR}/script/asan/fake_dlclose/libdlclose.so")
  set(PDM_RUN_PRELOAD_SRUN "LD_PRELOAD=${PRELOAD_ASAN}:${PROJECT_BINARY_DIR}/script/asan/fake_dlclose/libdlclose.so")
  set(PDM_RUN_LSAN "LSAN_OPTIONS=suppressions=${PROJECT_SOURCE_DIR}/script/asan/asan.supp")

  if (PDM_ENABLE_DOC)
    message(FATAL_ERROR "Documention can not be built in Sanitize CMAKE_BUILD_TYPE")
  endif ()
endif()

if(DEFINED ENV{PYTHONPATH})
  set (PDM_RUN_PYTHONPATH "PYTHONPATH=${PROJECT_BINARY_DIR}/Cython:$ENV{PYTHONPATH}")
else()
  set (PDM_RUN_PYTHONPATH "PYTHONPATH=${PROJECT_BINARY_DIR}/Cython")
endif()

configure_file(${PROJECT_SOURCE_DIR}/script/pdm_run.in "${PROJECT_BINARY_DIR}/pdm_run"
               @ONLY)

 #
 # Voir comment mieux faire !

configure_file(${PROJECT_SOURCE_DIR}/script/pdm_run.in "${PROJECT_BINARY_DIR}/test/pdm_run"
               @ONLY)

if(PDM_ENABLE_EXTENSION_PDMA)
 configure_file(${PROJECT_SOURCE_DIR}/script/pdm_run.in "${PROJECT_BINARY_DIR}/extensions/paradigma/test/pdm_run"
               @ONLY)
endif()

configure_file(${PROJECT_SOURCE_DIR}/script/pdm_launch.in "${PROJECT_BINARY_DIR}/pdm_launch"
               @ONLY)

configure_file(${PROJECT_SOURCE_DIR}/script/pdm_launch_python.in "${PROJECT_BINARY_DIR}/pdm_launch_python"
               @ONLY)


#------------------------------------------------------------------------------

if(PDM_ENABLE_UNIT_TEST)
  include(FetchContent)
  FetchContent_Declare(
    doctest
    GIT_REPOSITORY https://github.com/doctest/doctest.git
    GIT_TAG        v2.4.11
  )
  # FetchContent_Declare(
  #   doctest
  #    GIT_REPOSITORY https://github.com/BerengerBerthoul/doctest.git
  #    GIT_TAG        f13a00607463ce211004be23b5aa0bdf24d3ab61
  # )
  FetchContent_GetProperties(doctest)
  if(NOT doctest_POPULATED)
    FetchContent_Populate(doctest)
    add_subdirectory(${doctest_SOURCE_DIR} ${doctest_BINARY_DIR})
  endif()
endif()

#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# Création d'un package pdmConfig.cmake pour importer facilement paradigm
# pour les clients utilisant CMake
#------------------------------------------------------------------------------

configure_file(
  "${CMAKE_CURRENT_SOURCE_DIR}/cmake/templates/pdmConfig.cmake.in"
  pdmConfig.cmake
  @ONLY
)
install(FILES "${CMAKE_CURRENT_BINARY_DIR}/pdmConfig.cmake"
  DESTINATION lib/cmake/pdm
)

configure_file(
  "${CMAKE_CURRENT_SOURCE_DIR}/cmake/templates/pdmConfigVersion.cmake.in"
  pdmConfigVersion.cmake
  @ONLY
)
install(FILES "${CMAKE_CURRENT_BINARY_DIR}/pdmConfigVersion.cmake"
  DESTINATION lib/cmake/pdm
)

#------------------------------------------------------------------------------
# Affichage du message de fin d'installation
#------------------------------------------------------------------------------


add_subdirectory(${PDM_CMAKE_DIR}/post-install)
