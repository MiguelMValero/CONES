file (COPY ${CMAKE_CURRENT_SOURCE_DIR}/prepare_ci DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/")

if (CMAKE_BUILD_TYPE STREQUAL "Sanitize")
  add_subdirectory(asan)	
endif()

