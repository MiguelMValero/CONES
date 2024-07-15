#----------------------------------------------------------------
# Generated CMake target import file for configuration "Release".
#----------------------------------------------------------------

# Commands may need to know the format version.
set(CMAKE_IMPORT_FILE_VERSION 1)

# Import target "cwp::cwp_static" for configuration "Release"
set_property(TARGET cwp::cwp_static APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
set_target_properties(cwp::cwp_static PROPERTIES
  IMPORTED_LINK_INTERFACE_LANGUAGES_RELEASE "C;CXX"
  IMPORTED_LOCATION_RELEASE "${_IMPORT_PREFIX}/lib/libcwp.a"
  )

list(APPEND _IMPORT_CHECK_TARGETS cwp::cwp_static )
list(APPEND _IMPORT_CHECK_FILES_FOR_cwp::cwp_static "${_IMPORT_PREFIX}/lib/libcwp.a" )

# Commands beyond this point should not need to know the version.
set(CMAKE_IMPORT_FILE_VERSION)
