# - Try to find PTSCOTCH
# Once done this will define
#
#  PTSCOTCH_FOUND        - system has found PTSCOTCH
#  PTSCOTCH_INCLUDE_DIRS - include directories for PTSCOTCH
#  PTSCOTCH_LIBARIES     - libraries for PTSCOTCH
#  PTSCOTCH_VERSION      - version for PTSCOTCH
#
include(CMakeFindDependencyMacro)

if (NOT PTSCOTCH_FOUND)
  set(PTSCOTCH_DIR "" CACHE PATH "Installation directory of PTSCOTCH library")
  if (NOT PTSCOTCH_FIND_QUIETLY)
    message(STATUS "A cache variable, namely PTSCOTCH_DIR, has been set to specify the install directory of PTSCOTCH")
  endif()
endif()

# PTSCOTCH depends on Threads, try to find it
include(CMakeFindDependencyMacro)
if (NOT THREADS_FOUND)
  if (PTSCOTCH_FIND_REQUIRED)
    find_dependency(Threads REQUIRED)
  else()
    find_dependency(Threads)
  endif()
endif()

# PTSCOTCH depends on MPI, try to find it
if (NOT MPI_FOUND)
  if (PTSCOTCH_FIND_REQUIRED)
    find_dependency(MPI REQUIRED)
  else()
    find_dependency(MPI)
  endif()
endif()

# Try to get environment variables
# SCOTCH_ROOT for Topaze modules
set(ENV_PTSCOTCH_DIR "")
if (DEFINED ENV{SCOTCH_ROOT})
  set(ENV_PTSCOTCH_DIR "$ENV{SCOTCH_ROOT}")
endif()
# else PTSCOTCH_DIR
if (DEFINED ENV{PTSCOTCH_DIR})
  set(ENV_PTSCOTCH_DIR "$ENV{PTSCOTCH_DIR}")
endif()

set(ENV_PTSCOTCH_INCDIR "")
if (DEFINED ENV{SCOTCH_INCDIR})
  set (ENV_PTSCOTCH_INCDIR "$ENV{SCOTCH_INCDIR}")
endif()
if (DEFINED ENV{PTSCOTCH_INCDIR})
  set (ENV_PTSCOTCH_INCDIR "$ENV{PTSCOTCH_INCDIR}")
endif()

unset(_inc_env)
if(ENV_PTSCOTCH_INCDIR)
  list(APPEND _inc_env "${ENV_PTSCOTCH_INCDIR}")
elseif(ENV_PTSCOTCH_DIR)
  list(APPEND _inc_env "${ENV_PTSCOTCH_DIR}")
  list(APPEND _inc_env "${ENV_PTSCOTCH_DIR}/include")
  list(APPEND _inc_env "${ENV_PTSCOTCH_DIR}/include/ptscotch")
else()
  if(WIN32)
    string(REPLACE ":" ";" _inc_env "$ENV{INCLUDE}")
  else()
    string(REPLACE ":" ";" _path_env "$ENV{INCLUDE}")
    list(APPEND _inc_env "${_path_env}")
    string(REPLACE ":" ";" _path_env "$ENV{C_INCLUDE_PATH}")
    list(APPEND _inc_env "${_path_env}")
    string(REPLACE ":" ";" _path_env "$ENV{CPATH}")
    list(APPEND _inc_env "${_path_env}")
    string(REPLACE ":" ";" _path_env "$ENV{INCLUDE_PATH}")
    list(APPEND _inc_env "${_path_env}")
  endif()
endif()
list(APPEND _inc_env "${CMAKE_PLATFORM_IMPLICIT_INCLUDE_DIRECTORIES}")
list(APPEND _inc_env "${CMAKE_C_IMPLICIT_INCLUDE_DIRECTORIES}")
list(REMOVE_DUPLICATES _inc_env)


# Try to find the ptscotch header in the given paths
# -------------------------------------------------

set(PTSCOTCH_hdrs_to_find "ptscotch.h;scotch.h")

# call cmake macro to find the header path
if(PTSCOTCH_INCDIR)
  foreach(ptscotch_hdr ${PTSCOTCH_hdrs_to_find})
    set(PTSCOTCH_${ptscotch_hdr}_DIRS "PTSCOTCH_${ptscotch_hdr}_DIRS-NOTFOUND")
    find_path(PTSCOTCH_${ptscotch_hdr}_DIRS
      NAMES ${ptscotch_hdr}
      HINTS ${PTSCOTCH_INCDIR})
    mark_as_advanced(PTSCOTCH_${ptscotch_hdr}_DIRS)
  endforeach()
else()
  if(PTSCOTCH_DIR)
    foreach(ptscotch_hdr ${PTSCOTCH_hdrs_to_find})
      set(PTSCOTCH_${ptscotch_hdr}_DIRS "PTSCOTCH_${ptscotch_hdr}_DIRS-NOTFOUND")
      find_path(PTSCOTCH_${ptscotch_hdr}_DIRS
        NAMES ${ptscotch_hdr}
        HINTS ${PTSCOTCH_DIR}
        PATH_SUFFIXES "include" "include/scotch")
      mark_as_advanced(PTSCOTCH_${ptscotch_hdr}_DIRS)
    endforeach()
  else()
    foreach(ptscotch_hdr ${PTSCOTCH_hdrs_to_find})
      set(PTSCOTCH_${ptscotch_hdr}_DIRS "PTSCOTCH_${ptscotch_hdr}_DIRS-NOTFOUND")
      find_path(PTSCOTCH_${ptscotch_hdr}_DIRS
        NAMES ${ptscotch_hdr}
        HINTS ${_inc_env}
        PATH_SUFFIXES "scotch")
      mark_as_advanced(PTSCOTCH_${ptscotch_hdr}_DIRS)
    endforeach()
  endif()
endif()

# If found, add path to cmake variable
# ------------------------------------
foreach(ptscotch_hdr ${PTSCOTCH_hdrs_to_find})
  if (PTSCOTCH_${ptscotch_hdr}_DIRS)
    list(APPEND PTSCOTCH_INCLUDE_DIRS "${PTSCOTCH_${ptscotch_hdr}_DIRS}")
  else ()
    if (NOT PTSCOTCH_FIND_QUIETLY)
      message(STATUS "Looking for ptscotch -- ${ptscotch_hdr} not found")
    endif()
  endif()
endforeach()
list(REMOVE_DUPLICATES PTSCOTCH_INCLUDE_DIRS)

# Looking for lib
# ---------------

# Add system library paths to search lib
# --------------------------------------
unset(_lib_env)
set(ENV_PTSCOTCH_LIBDIR "")
if (DEFINED ENV{SCOTCH_LIBDIR})
  set(ENV_PTSCOTCH_LIBDIR "$ENV{SCOTCH_LIBDIR}")
endif()
if (DEFINED ENV{PTSCOTCH_LIBDIR})
  set(ENV_PTSCOTCH_LIBDIR "$ENV{PTSCOTCH_LIBDIR}")
endif()
#
if(ENV_PTSCOTCH_LIBDIR)
  list(APPEND _lib_env "${ENV_PTSCOTCH_LIBDIR}")
elseif(ENV_PTSCOTCH_DIR)
  list(APPEND _lib_env "${ENV_PTSCOTCH_DIR}")
  list(APPEND _lib_env "${ENV_PTSCOTCH_DIR}/lib")
else()
  if(WIN32)
    string(REPLACE ":" ";" _lib_env "$ENV{LIB}")
  else()
    if(APPLE)
      string(REPLACE ":" ";" _lib_env "$ENV{DYLD_LIBRARY_PATH}")
    else()
      string(REPLACE ":" ";" _lib_env "$ENV{LD_LIBRARY_PATH}")
    endif()
    list(APPEND _lib_env "${CMAKE_PLATFORM_IMPLICIT_LINK_DIRECTORIES}")
    list(APPEND _lib_env "${CMAKE_C_IMPLICIT_LINK_DIRECTORIES}")
  endif()
endif()
list(REMOVE_DUPLICATES _lib_env)

# Try to find the ptscotch lib in the given paths
# ----------------------------------------------

set(PTSCOTCH_libs_to_find "ptscotch;ptscotcherr")
list(APPEND PTSCOTCH_libs_to_find "scotch;scotcherr")

# call cmake macro to find the lib path
if(PTSCOTCH_LIBDIR)
  foreach(ptscotch_lib ${PTSCOTCH_libs_to_find})
    set(PTSCOTCH_${ptscotch_lib}_LIBRARY "PTSCOTCH_${ptscotch_lib}_LIBRARY-NOTFOUND")
    find_library(PTSCOTCH_${ptscotch_lib}_LIBRARY
      NAMES ${ptscotch_lib}
      HINTS ${PTSCOTCH_LIBDIR})
  endforeach()
else()
  if(PTSCOTCH_DIR)
    foreach(ptscotch_lib ${PTSCOTCH_libs_to_find})
      set(PTSCOTCH_${ptscotch_lib}_LIBRARY "PTSCOTCH_${ptscotch_lib}_LIBRARY-NOTFOUND")
      find_library(PTSCOTCH_${ptscotch_lib}_LIBRARY
        NAMES ${ptscotch_lib}
        HINTS ${PTSCOTCH_DIR}
        PATH_SUFFIXES lib lib32 lib64)
    endforeach()
  else()
    foreach(ptscotch_lib ${PTSCOTCH_libs_to_find})
      set(PTSCOTCH_${ptscotch_lib}_LIBRARY "PTSCOTCH_${ptscotch_lib}_LIBRARY-NOTFOUND")
      find_library(PTSCOTCH_${ptscotch_lib}_LIBRARY
        NAMES ${ptscotch_lib}
        HINTS ${_lib_env})
    endforeach()
  endif()
endif()


set(PTSCOTCH_LIBRARIES "")
set(PTSCOTCH_LIBRARY_DIRS "")
# If found, add path to cmake variable
# ------------------------------------
foreach(ptscotch_lib ${PTSCOTCH_libs_to_find})

  if (PTSCOTCH_${ptscotch_lib}_LIBRARY)
    get_filename_component(${ptscotch_lib}_lib_path "${PTSCOTCH_${ptscotch_lib}_LIBRARY}" PATH)
    # set cmake variables
    list(APPEND PTSCOTCH_LIBRARIES "${PTSCOTCH_${ptscotch_lib}_LIBRARY}")
    list(APPEND PTSCOTCH_LIBRARY_DIRS "${${ptscotch_lib}_lib_path}")
  else ()
    if (NOT PTSCOTCH_FIND_QUIETLY)
      message(STATUS "Looking for ptscotch -- lib ${ptscotch_lib} not found")
    endif()
  endif ()

  mark_as_advanced(PTSCOTCH_${ptscotch_lib}_LIBRARY)

endforeach()
list(REMOVE_DUPLICATES PTSCOTCH_LIBRARY_DIRS)


# check a function to validate the find
if(PTSCOTCH_LIBRARIES)

  set(REQUIRED_LDFLAGS)
  set(REQUIRED_INCDIRS)
  set(REQUIRED_LIBDIRS)
  set(REQUIRED_LIBS)

  # PTSCOTCH
  if (PTSCOTCH_INCLUDE_DIRS)
    set(REQUIRED_INCDIRS  "${PTSCOTCH_INCLUDE_DIRS}")
  endif()
  if (PTSCOTCH_LIBRARY_DIRS)
    set(REQUIRED_LIBDIRS "${PTSCOTCH_LIBRARY_DIRS}")
  endif()
  set(REQUIRED_LIBS "${PTSCOTCH_LIBRARIES}")
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
  # THREADS
  if(CMAKE_THREAD_LIBS_INIT)
    list(APPEND REQUIRED_LIBS "${CMAKE_THREAD_LIBS_INIT}")
  endif()
  set(Z_LIBRARY "Z_LIBRARY-NOTFOUND")
  find_library(Z_LIBRARY NAMES z)
  mark_as_advanced(Z_LIBRARY)
  if(Z_LIBRARY)
    list(APPEND REQUIRED_LIBS "-lz")
  endif()
  set(M_LIBRARY "M_LIBRARY-NOTFOUND")
  find_library(M_LIBRARY NAMES m)
  mark_as_advanced(M_LIBRARY)
  if(M_LIBRARY)
    list(APPEND REQUIRED_LIBS "-lm")
  endif()
  set(RT_LIBRARY "RT_LIBRARY-NOTFOUND")
  find_library(RT_LIBRARY NAMES rt)
  mark_as_advanced(RT_LIBRARY)
  if(RT_LIBRARY)
    list(APPEND REQUIRED_LIBS "-lrt")
  endif()

  # set required libraries for link
  set(CMAKE_REQUIRED_INCLUDES "${REQUIRED_INCDIRS}")
  set(CMAKE_REQUIRED_LIBRARIES)
  foreach(lib_dir ${REQUIRED_LIBDIRS})
    list(APPEND CMAKE_REQUIRED_LIBRARIES "-L${lib_dir}")
  endforeach()
  list(APPEND CMAKE_REQUIRED_LIBRARIES "${REQUIRED_LIBS}")
  list(APPEND CMAKE_REQUIRED_FLAGS "${REQUIRED_FLAGS} ${REQUIRED_LDFLAGS}")
  string(REGEX REPLACE "^ -" "-" CMAKE_REQUIRED_LIBRARIES "${CMAKE_REQUIRED_LIBRARIES}")

  # test link
  unset(PTSCOTCH_WORKS CACHE)
  include(CheckFunctionExists)
  check_function_exists(SCOTCH_dgraphInit PTSCOTCH_WORKS)
  mark_as_advanced(PTSCOTCH_WORKS)

  if(PTSCOTCH_WORKS)
    # save link with dependencies
    set(PTSCOTCH_LIBRARIES_DEP "${REQUIRED_LIBS}")
    set(PTSCOTCH_LIBRARY_DIRS_DEP "${REQUIRED_LIBDIRS}")
    set(PTSCOTCH_INCLUDE_DIRS_DEP "${REQUIRED_INCDIRS}")
    set(PTSCOTCH_LINKER_FLAGS "${REQUIRED_LDFLAGS}")
    list(REMOVE_DUPLICATES PTSCOTCH_LIBRARY_DIRS_DEP)
    list(REMOVE_DUPLICATES PTSCOTCH_INCLUDE_DIRS_DEP)
    list(REMOVE_DUPLICATES PTSCOTCH_LINKER_FLAGS)
  else()
    if(NOT PTSCOTCH_FIND_QUIETLY)
      message(STATUS "Looking for PTSCOTCH : test of SCOTCH_dgraphInit with PTSCOTCH library fails")
      message(STATUS "CMAKE_REQUIRED_LIBRARIES: ${CMAKE_REQUIRED_LIBRARIES}")
      message(STATUS "CMAKE_REQUIRED_INCLUDES: ${CMAKE_REQUIRED_INCLUDES}")
      message(STATUS "Check in CMakeFiles/CMakeError.log to figure out why it fails")
    endif()
  endif()

  # Run a C test to get version
  set(PTSCOTCH_CONFIG_TEST_VERSION_C
    "${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/scotch_config_test_version.c")

  file(WRITE ${PTSCOTCH_CONFIG_TEST_VERSION_C} "
#include <stdint.h>
#include <stdio.h>
#include <mpi.h>
#include <ptscotch.h>

int main() {
  printf(\"%i.%i.%i\", SCOTCH_VERSION,
	               SCOTCH_RELEASE,
	               SCOTCH_PATCHLEVEL);
  return 0;
}
")

  try_run(
    PTSCOTCH_CONFIG_TEST_VERSION_EXITCODE
    PTSCOTCH_CONFIG_TEST_VERSION_COMPILED
    ${CMAKE_CURRENT_BINARY_DIR}
    ${PTSCOTCH_CONFIG_TEST_VERSION_C}
    CMAKE_FLAGS
      "-DINCLUDE_DIRECTORIES:STRING=${CMAKE_REQUIRED_INCLUDES}"
      "-DLINK_LIBRARIES:STRING=${CMAKE_REQUIRED_LIBRARIES}"
    COMPILE_OUTPUT_VARIABLE COMPILE_OUTPUT
    RUN_OUTPUT_VARIABLE OUTPUT
  )


  if ((NOT PTSCOTCH_CONFIG_TEST_VERSION_COMPILED)
    OR (NOT (PTSCOTCH_CONFIG_TEST_VERSION_EXITCODE EQUAL 0)))
    message( "WARNING: Unable to determine PTSCOTCH version")
    set(PTSCOTCH_VERSION_OK TRUE)
    set(PTSCOTCH_VERSION "??.??.??" CACHE STRING "PTSCOTCH version")
    set(PTSCOTCH_TEST_COMPILE TRUE)
  else ()
    set(PTSCOTCH_VERSION "${OUTPUT}" CACHE STRING "PTSCOTCH version")
    mark_as_advanced(PTSCOTCH_VERSION)

   if (PTSCOTCH_FIND_VERSION)
      # Check if version found is >= required version
      if (NOT "${PTSCOTCH_VERSION}" VERSION_LESS "${PTSCOTCH_FIND_VERSION}")
  	set(PTSCOTCH_VERSION_OK TRUE)
      endif()
   else()
   # No specific version requested
       set(PTSCOTCH_VERSION_OK TRUE)
   endif()
   mark_as_advanced(PTSCOTCH_VERSION_OK)
  endif ()

  set(CMAKE_REQUIRED_INCLUDES)
  set(CMAKE_REQUIRED_FLAGS)
  set(CMAKE_REQUIRED_LIBRARIES)
else()
  set(PTSCOTCH_VERSION_OK FALSE)
  set(PTSCOTCH_VERSION "??.??.??" CACHE STRING "PTSCOTCH version")
endif()


if (PTSCOTCH_LIBRARIES)
  list(GET PTSCOTCH_LIBRARIES 0 first_lib)
  get_filename_component(first_lib_path "${first_lib}" PATH)
  if (${first_lib_path} MATCHES "/lib(32|64)?$")
    string(REGEX REPLACE "/lib(32|64)?$" "" not_cached_dir "${first_lib_path}")
    set(PTSCOTCH_DIR_FOUND "${not_cached_dir}" CACHE PATH "Installation directory of PTSCOTCH library" FORCE)
  else()
    set(PTSCOTCH_DIR_FOUND "${first_lib_path}" CACHE PATH "Installation directory of PTSCOTCH library" FORCE)
  endif()
endif()
mark_as_advanced(PTSCOTCH_DIR)
mark_as_advanced(PTSCOTCH_DIR_FOUND)

#
# Standard package handling
#
INCLUDE(FindPackageHandleStandardArgs)
find_package_handle_standard_args(PTSCOTCH
                                  "PTSCOTCH could not be found. Be sure to set PTSCOTCH_DIR."
                                  PTSCOTCH_LIBRARIES
                                  PTSCOTCH_INCLUDE_DIRS
                                  PTSCOTCH_WORKS
                                  PTSCOTCH_VERSION
                                  PTSCOTCH_VERSION_OK)
