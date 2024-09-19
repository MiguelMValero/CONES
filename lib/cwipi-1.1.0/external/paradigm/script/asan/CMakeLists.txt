if (CMAKE_BUILD_TYPE STREQUAL "Sanitize")
  add_subdirectory(fake_dlclose)	
endif()
