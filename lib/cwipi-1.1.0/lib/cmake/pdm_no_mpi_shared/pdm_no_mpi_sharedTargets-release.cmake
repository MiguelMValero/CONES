#----------------------------------------------------------------
# Generated CMake target import file for configuration "Release".
#----------------------------------------------------------------

# Commands may need to know the format version.
set(CMAKE_IMPORT_FILE_VERSION 1)

# Import target "pdm::pdm_no_mpi_shared" for configuration "Release"
set_property(TARGET pdm::pdm_no_mpi_shared APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
set_target_properties(pdm::pdm_no_mpi_shared PROPERTIES
  IMPORTED_LOCATION_RELEASE "${_IMPORT_PREFIX}/lib/libpdm_no_mpi.so.2.4.1"
  IMPORTED_SONAME_RELEASE "libpdm_no_mpi.so.2.4.1"
  )

list(APPEND _IMPORT_CHECK_TARGETS pdm::pdm_no_mpi_shared )
list(APPEND _IMPORT_CHECK_FILES_FOR_pdm::pdm_no_mpi_shared "${_IMPORT_PREFIX}/lib/libpdm_no_mpi.so.2.4.1" )

# Commands beyond this point should not need to know the version.
set(CMAKE_IMPORT_FILE_VERSION)
