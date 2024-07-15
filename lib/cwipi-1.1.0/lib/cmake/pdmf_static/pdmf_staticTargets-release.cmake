#----------------------------------------------------------------
# Generated CMake target import file for configuration "Release".
#----------------------------------------------------------------

# Commands may need to know the format version.
set(CMAKE_IMPORT_FILE_VERSION 1)

# Import target "pdm::pdmf_static" for configuration "Release"
set_property(TARGET pdm::pdmf_static APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
set_target_properties(pdm::pdmf_static PROPERTIES
  IMPORTED_LINK_INTERFACE_LANGUAGES_RELEASE "Fortran"
  IMPORTED_LOCATION_RELEASE "${_IMPORT_PREFIX}/lib/libpdmf.a"
  )

list(APPEND _IMPORT_CHECK_TARGETS pdm::pdmf_static )
list(APPEND _IMPORT_CHECK_FILES_FOR_pdm::pdmf_static "${_IMPORT_PREFIX}/lib/libpdmf.a" )

# Commands beyond this point should not need to know the version.
set(CMAKE_IMPORT_FILE_VERSION)
