# - Try to find ParMETIS
# Once done this will define
#
#  PARMETIS_FOUND        - system has ParMETIS
#  PARMETIS_INCLUDE_DIRS - include directories for ParMETIS
#  PARMETIS_LIBRARIES    - libraries for ParMETIS
#  PARMETIS_VERSION      - version for ParMETIS
#=============================================================================
if (NOT ParMETIS_FOUND)
  set(PARMETIS_DIR "" CACHE PATH "Installation directory of PARMETIS library")
  if (NOT PARMETIS_FIND_QUIETLY)
	  message(STATUS "A cache variable, namely PARMETIS_DIR, has been set to specify the install directory of PARMETIS")
  endif()
endif()

include(CMakeFindDependencyMacro)
# PARMETIS depends on MPI, try to find it
if (NOT MPI_FOUND)
  if (PARMETIS_FIND_REQUIRED)
    find_dependency(MPI REQUIRED)
  else()
    find_dependency(MPI)
  endif()
endif()

# Try to get environment variables
# PARMETIS_ROOT for Topaze modules
set(ENV_PARMETIS_DIR "")
if (DEFINED ENV{PARMETIS_ROOT})
  set(ENV_PARMETIS_DIR "$ENV{PARMETIS_ROOT}")
endif()
# else PARMETIS_DIR
if (DEFINED ENV{PARMETIS_DIR})
  set(ENV_PARMETIS_DIR "$ENV{PARMETIS_DIR}")
endif()

set(ENV_PARMETIS_INCDIR "")
if (DEFINED ENV{PARMETIS_INCDIR})
  set (ENV_PARMETIS_INCDIR "$ENV{PARMETIS_INCDIR}")
elseif(DEFINED ENV{PARMETIS_INCLUDE_DIR})
  set (ENV_PARMETIS_INCDIR "$ENV{PARMETIS_INCLUDE_DIR}")
endif()
# message("ENV_PARMETIS_DIR" ${ENV_PARMETIS_DIR})
# message("ENV_PARMETIS_INCDIR" ${ENV_PARMETIS_INCDIR})

# Looking for include
# -------------------

# Add system include paths to search include
# ------------------------------------------
unset(_inc_env)

if(ENV_PARMETIS_INCDIR)
  list(APPEND _inc_env "${ENV_PARMETIS_INCDIR}")
elseif(ENV_PARMETIS_DIR)
  list(APPEND _inc_env "${ENV_PARMETIS_DIR}")
  list(APPEND _inc_env "${ENV_PARMETIS_DIR}/include")
  list(APPEND _inc_env "${ENV_PARMETIS_DIR}/include/parmetis")
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

# Try to find the parmetis header in the given paths
# -------------------------------------------------

set(PARMETIS_hdrs_to_find "parmetis.h")

# call cmake macro to find the header path
if(PARMETIS_INCDIR)
  foreach(parmetis_hdr ${PARMETIS_hdrs_to_find})
    set(PARMETIS_${parmetis_hdr}_DIRS "PARMETIS_${parmetis_hdr}_DIRS-NOTFOUND")
    find_path(PARMETIS_${parmetis_hdr}_DIRS
      NAMES ${parmetis_hdr}
      HINTS ${PARMETIS_INCDIR})
    mark_as_advanced(PARMETIS_${parmetis_hdr}_DIRS)
  endforeach()
else()
  if(PARMETIS_DIR)
    foreach(parmetis_hdr ${PARMETIS_hdrs_to_find})
      set(PARMETIS_${parmetis_hdr}_DIRS "PARMETIS_${parmetis_hdr}_DIRS-NOTFOUND")
      find_path(PARMETIS_${parmetis_hdr}_DIRS
        NAMES ${parmetis_hdr}
        HINTS ${PARMETIS_DIR}
        PATH_SUFFIXES "include" "include/parmetis")
      mark_as_advanced(PARMETIS_${parmetis_hdr}_DIRS)
    endforeach()
  else()
    foreach(parmetis_hdr ${PARMETIS_hdrs_to_find})
      set(PARMETIS_${parmetis_hdr}_DIRS "PARMETIS_${parmetis_hdr}_DIRS-NOTFOUND")
      find_path(PARMETIS_${parmetis_hdr}_DIRS
        NAMES ${parmetis_hdr}
        HINTS ${_inc_env}
        PATH_SUFFIXES "parmetis")
      mark_as_advanced(PARMETIS_${parmetis_hdr}_DIRS)
    endforeach()
  endif()
endif()

# If found, add path to cmake variable
# ------------------------------------
foreach(parmetis_hdr ${PARMETIS_hdrs_to_find})
  if (PARMETIS_${parmetis_hdr}_DIRS)
    list(APPEND PARMETIS_INCLUDE_DIRS "${PARMETIS_${parmetis_hdr}_DIRS}")
  else ()
    if (NOT PARMETIS_FIND_QUIETLY)
      message(STATUS "Looking for parmetis -- ${parmetis_hdr} not found")
    endif()
  endif()
endforeach()
list(REMOVE_DUPLICATES PARMETIS_INCLUDE_DIRS)

# Looking for lib
# ---------------

# Add system library paths to search lib
# --------------------------------------
unset(_lib_env)
set(ENV_PARMETIS_LIBDIR "")
if (DEFINED ENV{PARMETIS_LIBDIR})
  set(ENV_PARMETIS_LIBDIR "$ENV{PARMETIS_LIBDIR}")
endif()
#
if(ENV_PARMETIS_LIBDIR)
  list(APPEND _lib_env "${ENV_PARMETIS_LIBDIR}")
elseif(ENV_PARMETIS_DIR)
  list(APPEND _lib_env "${ENV_PARMETIS_DIR}")
  list(APPEND _lib_env "${ENV_PARMETIS_DIR}/lib")
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

# Try to find the parmetis lib in the given paths
# ----------------------------------------------

set(PARMETIS_libs_to_find "parmetis")

# call cmake macro to find the lib path
if(PARMETIS_LIBDIR)
  foreach(parmetis_lib ${PARMETIS_libs_to_find})
    set(PARMETIS_${parmetis_lib}_LIBRARY "PARMETIS_${parmetis_lib}_LIBRARY-NOTFOUND")
    find_library(PARMETIS_${parmetis_lib}_LIBRARY
      NAMES ${parmetis_lib}
      HINTS ${PARMETIS_LIBDIR})
  endforeach()
else()
  if(PARMETIS_DIR)
    foreach(parmetis_lib ${PARMETIS_libs_to_find})
      set(PARMETIS_${parmetis_lib}_LIBRARY "PARMETIS_${parmetis_lib}_LIBRARY-NOTFOUND")
      find_library(PARMETIS_${parmetis_lib}_LIBRARY
        NAMES ${parmetis_lib}
        HINTS ${PARMETIS_DIR}
        PATH_SUFFIXES lib lib32 lib64)
    endforeach()
  else()
    foreach(parmetis_lib ${PARMETIS_libs_to_find})
      set(PARMETIS_${parmetis_lib}_LIBRARY "PARMETIS_${parmetis_lib}_LIBRARY-NOTFOUND")
      find_library(PARMETIS_${parmetis_lib}_LIBRARY
        NAMES ${parmetis_lib}
        HINTS ${_lib_env})
    endforeach()
  endif()
endif()


set(PARMETIS_LIBRARIES "")
set(PARMETIS_LIBRARY_DIRS "")
foreach(parmetis_lib ${PARMETIS_libs_to_find})
  if (PARMETIS_${parmetis_lib}_LIBRARY)
    get_filename_component(${parmetis_lib}_lib_path "${PARMETIS_${parmetis_lib}_LIBRARY}" PATH)
    # set cmake variables
    list(APPEND PARMETIS_LIBRARIES "${PARMETIS_${parmetis_lib}_LIBRARY}")
    list(APPEND PARMETIS_LIBRARY_DIRS "${${parmetis_lib}_lib_path}")
  else ()
    if (NOT PARMETIS_FIND_QUIETLY)
      message(STATUS "Looking for parmetis -- lib ${parmetis_lib} not found")
    endif()
  endif ()

  mark_as_advanced(PARMETIS_${parmetis_lib}_LIBRARY)
endforeach()
list(REMOVE_DUPLICATES PARMETIS_LIBRARY_DIRS)



# Now search for METIS
# =============================================================
# Use environment variables
# METIS_ROOT for Topaze modules
set(ENV_METIS_DIR "${PARMETIS_DIR}")
if (DEFINED ENV{METIS_ROOT})
  set(ENV_METIS_DIR "$ENV{METIS_ROOT}")
endif()
if (DEFINED ENV{METIS_DIR})
  set(ENV_METIS_DIR "$ENV{METIS_DIR}")
endif()

set(ENV_METIS_INCDIR "")
if (DEFINED ENV{METIS_INCDIR})
  set (ENV_METIS_INCDIR "$ENV{METIS_INCDIR}")
endif()


# Looking for include
# -------------------

# Add system include paths to search include
# ------------------------------------------
unset(_inc_env)

if(ENV_METIS_INCDIR)
  list(APPEND _inc_env "${ENV_METIS_INCDIR}")
elseif(ENV_METIS_DIR)
  list(APPEND _inc_env "${ENV_METIS_DIR}")
  list(APPEND _inc_env "${ENV_METIS_DIR}/include")
  list(APPEND _inc_env "${ENV_METIS_DIR}/include/metis")
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

# Try to find the metis header in the given paths
# -------------------------------------------------

set(METIS_hdrs_to_find "metis.h")

# call cmake macro to find the header path
if(METIS_INCDIR)
  foreach(metis_hdr ${METIS_hdrs_to_find})
    set(METIS_${metis_hdr}_DIRS "METIS_${metis_hdr}_DIRS-NOTFOUND")
    find_path(METIS_${metis_hdr}_DIRS
      NAMES ${metis_hdr}
      HINTS ${METIS_INCDIR})
    mark_as_advanced(METIS_${metis_hdr}_DIRS)
  endforeach()
else()
  if(METIS_DIR)
    foreach(metis_hdr ${METIS_hdrs_to_find})
      set(METIS_${metis_hdr}_DIRS "METIS_${metis_hdr}_DIRS-NOTFOUND")
      find_path(METIS_${metis_hdr}_DIRS
        NAMES ${metis_hdr}
        HINTS ${METIS_DIR}
        PATH_SUFFIXES "include" "include/metis")
      mark_as_advanced(METIS_${metis_hdr}_DIRS)
    endforeach()
  else()
    foreach(metis_hdr ${METIS_hdrs_to_find})
      set(METIS_${metis_hdr}_DIRS "METIS_${metis_hdr}_DIRS-NOTFOUND")
      find_path(METIS_${metis_hdr}_DIRS
        NAMES ${metis_hdr}
        HINTS ${_inc_env}
        PATH_SUFFIXES "metis")
      mark_as_advanced(METIS_${metis_hdr}_DIRS)
    endforeach()
  endif()
endif()

# If found, add path to cmake variable
# ------------------------------------
foreach(metis_hdr ${METIS_hdrs_to_find})
  if (METIS_${metis_hdr}_DIRS)
    list(APPEND METIS_INCLUDE_DIRS "${METIS_${metis_hdr}_DIRS}")
  else ()
    if (NOT METIS_FIND_QUIETLY)
      message(STATUS "Looking for metis -- ${metis_hdr} not found")
    endif()
  endif()
endforeach()
list(REMOVE_DUPLICATES METIS_INCLUDE_DIRS)

# Looking for metis lib
# -------------------------

# Add system library paths to search lib
# --------------------------------------
unset(_lib_env)
set(ENV_METIS_LIBDIR "")
if (DEFINED ENV{METIS_LIBDIR})
  set(ENV_METIS_LIBDIR "$ENV{METIS_LIBDIR}")
endif()
#
if(ENV_METIS_LIBDIR)
  list(APPEND _lib_env "${ENV_METIS_LIBDIR}")
elseif(ENV_METIS_DIR)
  list(APPEND _lib_env "${ENV_METIS_DIR}")
  list(APPEND _lib_env "${ENV_METIS_DIR}/lib")
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

# Try to find the metis lib in the given paths
# ----------------------------------------------

set(METIS_libs_to_find "metis")

# call cmake macro to find the lib path
if(METIS_LIBDIR)
  foreach(metis_lib ${METIS_libs_to_find})
    set(METIS_${metis_lib}_LIBRARY "METIS_${metis_lib}_LIBRARY-NOTFOUND")
    find_library(METIS_${metis_lib}_LIBRARY
      NAMES ${metis_lib}
      HINTS ${METIS_LIBDIR})
  endforeach()
else()
  if(METIS_DIR)
    foreach(metis_lib ${METIS_libs_to_find})
      set(METIS_${metis_lib}_LIBRARY "METIS_${metis_lib}_LIBRARY-NOTFOUND")
      find_library(METIS_${metis_lib}_LIBRARY
        NAMES ${metis_lib}
        HINTS ${METIS_DIR}
        PATH_SUFFIXES lib lib32 lib64)
    endforeach()
  else()
    foreach(metis_lib ${METIS_libs_to_find})
      set(METIS_${metis_lib}_LIBRARY "METIS_${metis_lib}_LIBRARY-NOTFOUND")
      find_library(METIS_${metis_lib}_LIBRARY
        NAMES ${metis_lib}
        HINTS ${_lib_env})
    endforeach()
  endif()
endif()

set(METIS_LIBRARIES "")
set(METIS_LIBRARY_DIRS "")
# If found, add path to cmake variable
# ------------------------------------
foreach(metis_lib ${METIS_libs_to_find})

  if (METIS_${metis_lib}_LIBRARY)
    get_filename_component(${metis_lib}_lib_path "${METIS_${metis_lib}_LIBRARY}" PATH)
    # set cmake variables
    list(APPEND METIS_LIBRARIES "${METIS_${metis_lib}_LIBRARY}")
    list(APPEND METIS_LIBRARY_DIRS "${${metis_lib}_lib_path}")
  else ()
    if (NOT METIS_FIND_QUIETLY)
      message(STATUS "Looking for metis -- lib ${metis_lib} not found")
    endif()
  endif ()

  mark_as_advanced(METIS_${metis_lib}_LIBRARY)

endforeach()
list(REMOVE_DUPLICATES METIS_LIBRARY_DIRS)


# check and validate the findings
if (PARMETIS_LIBRARIES)
  set(REQUIRED_LDFLAGS)
  set(REQUIRED_INCDIRS)
  set(REQUIRED_LIBDIRS)
  set(REQUIRED_LIBS)

  # PARMETIS
  if (PARMETIS_INCLUDE_DIRS)
    set(REQUIRED_INCDIRS  "${PARMETIS_INCLUDE_DIRS}")
  endif()
  if (METIS_INCLUDE_DIRS)
    list(APPEND REQUIRED_INCDIRS "${METIS_INCLUDE_DIRS}")
  endif()
  if (PARMETIS_LIBRARY_DIRS)
    set(REQUIRED_LIBDIRS "${PARMETIS_LIBRARY_DIRS}")
  endif()
  set(REQUIRED_LIBS "${PARMETIS_LIBRARIES}")
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
  foreach(lib_dir ${REQUIRED_LIBDIRS})
    list(APPEND CMAKE_REQUIRED_LIBRARIES "-L${lib_dir}")
  endforeach()
  list(APPEND CMAKE_REQUIRED_LIBRARIES "${REQUIRED_LIBS}")
  list(APPEND CMAKE_REQUIRED_FLAGS "${REQUIRED_FLAGS} ${REQUIRED_LDFLAGS}")
  string(REGEX REPLACE "^ -" "-" CMAKE_REQUIRED_LIBRARIES "${CMAKE_REQUIRED_LIBRARIES}")
  # Check ParMETIS version
  set(PARMETIS_CONFIG_TEST_VERSION_C
    "${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/parmetis_config_test_version.c")

  file(WRITE ${PARMETIS_CONFIG_TEST_VERSION_C} "
#include <stdio.h>
#include \"parmetis.h\"

int main() {
#ifdef PARMETIS_SUBMINOR_VERSION
  printf(\"%i.%i.%i\",
         PARMETIS_MAJOR_VERSION,
         PARMETIS_MINOR_VERSION,
         PARMETIS_SUBMINOR_VERSION);
#else
  printf(\"%i.%i\\n\",
         PARMETIS_MAJOR_VERSION,
         PARMETIS_MINOR_VERSION);
#endif
  return 0;
}
")

  try_run(
     PARMETIS_CONFIG_TEST_VERSION_EXITCODE
     PARMETIS_CONFIG_TEST_VERSION_COMPILED
     ${CMAKE_CURRENT_BINARY_DIR}
     ${PARMETIS_CONFIG_TEST_VERSION_C}
     CMAKE_FLAGS
        "-DINCLUDE_DIRECTORIES:STRING=${CMAKE_REQUIRED_INCLUDES}"
        "-DLINK_LIBRARIES:STRING=${CMAKE_REQUIRED_LIBRARIES}"
     COMPILE_OUTPUT_VARIABLE PARMETIS_CONFIG_TEST_VERSION_COMPILE_OUTPUT
     RUN_OUTPUT_VARIABLE PARMETIS_CONFIG_TEST_VERSION_OUTPUT
  )

  if ((NOT PARMETIS_CONFIG_TEST_VERSION_COMPILED)
      OR (NOT (PARMETIS_CONFIG_TEST_VERSION_EXITCODE EQUAL 0)))
      message("WARNING: Unable to determine ParMETIS version")
      set(PARMETIS_VERSION_OK FALSE)
      set(PARMETIS_VERSION "??.??.??" CACHE STRING "ParMETIS Version")

    else ()
      set(PARMETIS_VERSION ${PARMETIS_CONFIG_TEST_VERSION_OUTPUT} CACHE STRING "ParMETIS Version")
      mark_as_advanced(PARMETIS_VERSION)
      if (ParMETIS_FIND_VERSION)
        # Check if version found is >= required version
        if (NOT "${PARMETIS_VERSION}" VERSION_LESS "${ParMETIS_FIND_VERSION}")
          set(PARMETIS_VERSION_OK TRUE)
        endif()
      else()
        # No specific version requested
          set(PARMETIS_VERSION_OK TRUE)
      endif()
  endif()
  mark_as_advanced(PARMETIS_VERSION_OK)
else()
  set(PARMETIS_VERSION_OK FALSE)
  set(PARMETIS_VERSION "??.??.??" CACHE STRING "ParMETIS version")
endif()


mark_as_advanced(PARMETIS_DIR)
mark_as_advanced(PARMETIS_DIR_FOUND)
#
# Standard package handling
#
INCLUDE(FindPackageHandleStandardArgs)
find_package_handle_standard_args(ParMETIS
	                                "ParMETIS could not be found. Be sure to set PARMETIS_DIR."
	                                PARMETIS_LIBRARIES
	                                PARMETIS_INCLUDE_DIRS
	                                METIS_LIBRARIES
	                                METIS_INCLUDE_DIRS
	                                PARMETIS_VERSION
                                  PARMETIS_VERSION_OK)

mark_as_advanced(PARMETIS_LIBRARIES
                PARMETIS_INCLUDE_DIRS
	              PARMETIS_VERSION
	              PARMETIS_VERSION_OK)
