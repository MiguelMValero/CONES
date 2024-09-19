
if (CMAKE_BUILD_TYPE STREQUAL "Sanitize" AND CWP_ENABLE_EXTERNAL_PDM)
  add_subdirectory(asan)	
endif()

