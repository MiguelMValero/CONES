#----------------------------------------------------------------
# Generated CMake target import file for configuration "Release".
#----------------------------------------------------------------

# Commands may need to know the format version.
set(CMAKE_IMPORT_FILE_VERSION 1)

# Import target "pdm::pdm_mpi_static" for configuration "Release"
set_property(TARGET pdm::pdm_mpi_static APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
set_target_properties(pdm::pdm_mpi_static PROPERTIES
  IMPORTED_LINK_INTERFACE_LANGUAGES_RELEASE "C"
  IMPORTED_LOCATION_RELEASE "${_IMPORT_PREFIX}/lib/libpdm_mpi.a"
  )

list(APPEND _IMPORT_CHECK_TARGETS pdm::pdm_mpi_static )
list(APPEND _IMPORT_CHECK_FILES_FOR_pdm::pdm_mpi_static "${_IMPORT_PREFIX}/lib/libpdm_mpi.a" )

# Commands beyond this point should not need to know the version.
set(CMAKE_IMPORT_FILE_VERSION)
