#----------------------------------------------------------------
# Generated CMake target import file for configuration "Release".
#----------------------------------------------------------------

# Commands may need to know the format version.
set(CMAKE_IMPORT_FILE_VERSION 1)

# Import target "cwp::cwp_shared" for configuration "Release"
set_property(TARGET cwp::cwp_shared APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
set_target_properties(cwp::cwp_shared PROPERTIES
  IMPORTED_LOCATION_RELEASE "${_IMPORT_PREFIX}/lib/libcwp.so.1.1.0"
  IMPORTED_SONAME_RELEASE "libcwp.so.1.1.0"
  )

list(APPEND _IMPORT_CHECK_TARGETS cwp::cwp_shared )
list(APPEND _IMPORT_CHECK_FILES_FOR_cwp::cwp_shared "${_IMPORT_PREFIX}/lib/libcwp.so.1.1.0" )

# Commands beyond this point should not need to know the version.
set(CMAKE_IMPORT_FILE_VERSION)
