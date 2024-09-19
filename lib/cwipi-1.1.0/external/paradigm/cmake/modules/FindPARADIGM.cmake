#
# CMake module for ParaDIGM 
#
#  PARADIGM_FOUND             - system has ParaDiGM
#  PARADIGM_INCLUDE_DIRS      - include directories for ParaDiGM
#  PARADIGM_LIBRARIES         - libraries for ParaDiGM
#  PARADIGM_FORTRAN_LIBRARIES - libraries for Fortran API
#  PARADIGM_NO_MPI_LIBRARIES  - libraries for ParaDiGM without MPI
#  PARADIGMA_LIBRARIES        - libraries for ParaDiGM
#  PARADIGMA_FORTRAN_LIBRARIES- libraries for Fortran API
#  PARADIGM_VERSION           - version for ParaDiGM
#  PARADIGM_BINDING_PYTHON_PATH - Path of python binding
#

include(CMakeFindDependencyMacro)

if (NOT PARADIGM_FOUND)
  set(PARADIGM_DIR "" CACHE PATH "Installation directory of PARADIGM library")
  if (NOT PARADIGM_FIND_QUIETLY)
    message(STATUS "A cache variable, namely PARADIGM_DIR, has been set to specify the install directory of PARADIGM")
  endif()
endif()

if (NOT MPI_FOUND)
  if (PARADIGM_FIND_REQUIRED)
    find_dependency(MPI REQUIRED)
  else()
    find_dependency(MPI)
  endif()
endif()

set(ENV_PARADIGM_DIR "")
if (DEFINED ENV{PARADIGM_ROOT})
  set(ENV_PARADIGM_DIR "$ENV{PARADIGM_ROOT}")
else ()
  if (DEFINED ENV{PARADIGM_DIR})
    set(ENV_PARADIGM_DIR "$ENV{PARADIGM_DIR}")
  else ()
    if (DEFINED ENV{PARADIGM_HOME})
      set(ENV_PARADIGM_DIR "$ENV{PARADIGM_HOME}")
    endif ()  
  endif()
endif()


find_path(PARADIGM_INCLUDE_DIR
          NAMES    pdm.h pdm_version.h
          PATH_SUFFIXES include
          HINTS "${ENV_PARADIGM_DIR}"
)

find_path(PARADIGM_FORTRAN_INCLUDE_DIR
          NAMES    pdm.mod
          PATH_SUFFIXES include
          HINTS "${ENV_PARADIGM_DIR}"
)

find_library(PARADIGM_LIBRARY
             NAMES pdm
             HINTS "${ENV_PARADIGM_DIR}/lib"
)

find_library(PARADIGMA_LIBRARY
             NAMES pdma
             HINTS "${ENV_PARADIGM_DIR}/lib"
)

find_library(PARADIGM_FORTRAN_LIBRARY
             NAMES pdmf
             HINTS "${ENV_PARADIGM_DIR}/lib"
)

find_library(PARADIGMA_FORTRAN_LIBRARY
             NAMES pdmaf
             HINTS "${ENV_PARADIGM_DIR}/lib"
)

find_library(PARADIGM_MPI_LIBRARY
             NAMES pdm_mpi
             HINTS "${ENV_PARADIGM_DIR}/lib"
)

find_library(PARADIGM_NO_MPI_LIBRARY
             NAMES pdm_no_mpi
             HINTS "${ENV_PARADIGM_DIR}/lib"
)

find_library(PARADIGM_IO_LIBRARY
             NAMES pdm_io
             HINTS "${ENV_PARADIGM_DIR}/lib"
)
   
#find_file(PARADIGM_PYTHON_BINDING_LIBRARY_PATH
#             NAMES Pypdm.so
#             HINTS "${ENV_PARADIGM_DIR}"
#             PATHS "${ENV_PARADIGM_DIR}"
#)

#message(${PARADIGM_PYTHON_BINDING_LIBRARY_PATH} ${ENV_PARADIGM_DIR})
#message("toto=${PARADIGM_PYTHON_BINDING_LIBRARY_PATH}")                             

#if (PARADIGM_PYTHON_BINDING_LIBRARY_PATH)
#   string(REPLACE "Pypdm/" "" 
#          PARADIGM_BINDING_PYTHON_PATH 
#          ${PARADIGM_PYTHON_BINDING_LIBRARY_PATH})
#   message("${PARADIGM_PYTHON_BINDING_LIBRARY_PATH}")                             
#   message("${PARADIGM_BINDING_PYTHON_PATH}")                             
#endif()

set(PARADIGM_LIBRARIES   "")
set(PARADIGM_NO_MPI_LIBRARIES   "")
set(PARADIGMA_LIBRARIES   "")
set(PARADIGM_FORTRAN_LIBRARIES   "")
set(PARADIGMA_FORTRAN_LIBRARIES   "")
set(PARADIGM_INCLUDE_DIRS   "")

if (PARADIGM_INCLUDE_DIR AND
    PARADIGM_LIBRARY     AND
    PARADIGM_MPI_LIBRARY)
  set(PARADIGM_LIBRARIES ${PARADIGM_LIBRARY}  ${PARADIGM_IO_LIBRARY} ${PARADIGM_MPI_LIBRARY})
  set(PARADIGM_INCLUDE_DIRS ${PARADIGM_INCLUDE_DIR}) 
endif()

if (PARADIGM_INCLUDE_DIR AND
    PARADIGM_LIBRARY     AND
    PARADIGM_NO_MPI_LIBRARY)
  set(PARADIGM_NO_MPI_LIBRARIES ${PARADIGM_LIBRARY}  ${PARADIGM_IO_LIBRARY} ${PARADIGM_NO_MPI_LIBRARY})
endif()

if (PARADIGMA_LIBRARY)
  set(PARADIGMA_LIBRARIES ${PARADIGMA_LIBRARY})
endif()

if (PARADIGM_FORTRAN_INCLUDE_DIR AND 
    PARADIGM_FORTRAN_LIBRARY)
  set(PARADIGM_FORTRAN_LIBRARIES ${PARADIGM_FORTRAN_LIBRARY})
endif()

if (PARADIGMA_FORTRAN_LIBRARY)
  set(PARADIGMA_FORTRAN_LIBRARIES ${PARADIGMA_FORTRAN_LIBRARY})
endif()

# check and validate the findings
if (PARADIGM_INCLUDE_DIR AND PARADIGM_LIBRARIES)
  set(REQUIRED_LDFLAGS)
  set(REQUIRED_INCDIRS)
  set(REQUIRED_LIBS)

  # PARaDiGM
  if (PARADIGM_INCLUDE_DIR)
    set(REQUIRED_INCDIRS  "${PARADIGM_INCLUDE_DIR}")
  endif()

  set(REQUIRED_LIBS "${PARADIGM_LIBRARIES}")

  # MPI
  if (MPI_FOUND)
    if (MPI_C_INCLUDE_PATH)
      list(APPEND REQUIRED_INCDIRS "${MPI_C_INCLUDE_PATH}")
    endif()
    if (MPI_C_LINK_FLAGS)
      if (${MPI_C_LINK_FLAGS} MATCHES "  -")
        string(REGEX REPLACE " -" "-" MPI_C_LINK_FLAGS ${MPI_C_LINK_FLAGS})
      endif()
      list(APPEND REQUIRED_LDFLAGS "${MPI_C_LINK_FLAGS}")
    endif()
    list(APPEND REQUIRED_LIBS "${MPI_C_LIBRARIES}")
  endif()

  # set required libraries for link
  set(CMAKE_REQUIRED_INCLUDES "${REQUIRED_INCDIRS}")
  set(CMAKE_REQUIRED_LIBRARIES)

  #list(APPEND CMAKE_REQUIRED_LIBRARIES "${REQUIRED_LDFLAGS}")

  list(APPEND CMAKE_REQUIRED_LIBRARIES "${REQUIRED_LIBS}")

  list(APPEND CMAKE_REQUIRED_FLAGS "${REQUIRED_FLAGS} ${REQUIRED_LDFLAGS}")

  string(REGEX REPLACE "^ -" "-" CMAKE_REQUIRED_LIBRARIES "${CMAKE_REQUIRED_LIBRARIES}")

  # Check ParaDiGM version
  set(PARADIGM_CONFIG_TEST_VERSION_C
    "${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/paradigm_config_test_version.c")

  file(WRITE ${PARADIGM_CONFIG_TEST_VERSION_C} "
#include <stdio.h>
#include <stdlib.h>
#include \"pdm.h\"
#include \"pdm_version.h\"

int main() {

  char *version = PDM_version_get();
  printf(\"%s\",version);
  fflush(stdout);
  free(version);
  return 0;
}
")

  try_run(
     PARADIGM_CONFIG_TEST_VERSION_EXITCODE
     PARADIGM_CONFIG_TEST_VERSION_COMPILED
     ${CMAKE_CURRENT_BINARY_DIR}
     ${PARADIGM_CONFIG_TEST_VERSION_C}
     CMAKE_FLAGS
        "-DINCLUDE_DIRECTORIES:STRING=${CMAKE_REQUIRED_INCLUDES}"
        "-DLINK_LIBRARIES:STRING=${CMAKE_REQUIRED_LIBRARIES}"
     COMPILE_OUTPUT_VARIABLE PARADIGM_CONFIG_TEST_VERSION_COMPILE_OUTPUT
     RUN_OUTPUT_VARIABLE PARADIGM_CONFIG_TEST_VERSION_OUTPUT
  )

  if (NOT PARADIGM_CONFIG_TEST_VERSION_COMPILED)
    message("WARNING: Unable to determine ParaDiGM version")
  else ()  
    if (NOT (PARADIGM_CONFIG_TEST_VERSION_EXITCODE EQUAL 0))
      message("WARNING: Unable to determine ParaDiGM version")
    else ()
      set(PARADIGM_VERSION ${PARADIGM_CONFIG_TEST_VERSION_OUTPUT})
      mark_as_advanced(PARADIGM_VERSION)
    endif()
  endif()
endif()


if(Python_EXECUTABLE)

  # Retrieve the Py version
  EXECUTE_PROCESS(COMMAND
    ${Python_EXECUTABLE} -c "import Pypdm"
    OUTPUT_VARIABLE Mpi4Py_VERSION
    ERROR_VARIABLE  Mpi4Py_VERSION_ERROR
    OUTPUT_STRIP_TRAILING_WHITESPACE
    )
endif()

include(FindPackageHandleStandardArgs)

find_package_handle_standard_args (PARADIGM
  DEFAULT_MSG
  PARADIGM_LIBRARIES 
  PARADIGM_NO_MPI_LIBRARIES 
  PARADIGM_FORTRAN_LIBRARIES 
  PARADIGM_INCLUDE_DIRS 
  PARADIGMA_LIBRARIES
  PARADIGMA_FORTRAN_LIBRARIES
  PARADIGM_VERSION)
#  PARADIGM_BINDING_PYTHON_PATH)
 
include(FindPackageHandleStandardArgs)
find_package_check_version(${PARADIGM_VERSION} result
  HANDLE_VERSION_RANGE
  RESULT_MESSAGE_VARIABLE reason
  )  

if (result)
  message (STATUS "${reason}")
else()
  message (FATAL_ERROR "${reason}")
endif()

#  REASON_FAILURE_MESSAGE "ParaDiGM could not be found. Be sure to set PARADIGM_ROOT.")

# set(CMAKE_REQUIRED_INCLUDES ${PARADIGM_INCLUDE_DIR})
# set(CMAKE_REQUIRED_LIBRARIES ${PARADIGM_LIBRARIES})
# find_package_handle_standard_args(PARADIGM REQUIRED_VARS PARADIGM_INCLUDE_DIR)

