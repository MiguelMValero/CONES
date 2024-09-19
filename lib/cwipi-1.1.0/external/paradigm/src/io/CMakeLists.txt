
file(GLOB_RECURSE SOURCES
  pdm_printf.c
  pdm_error.c)

add_library(pdm_io OBJECT ${SOURCES})

target_include_directories(pdm_io
  INTERFACE $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}>
  PRIVATE   $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}>
  PRIVATE   $<BUILD_INTERFACE:${PDM_BINARY_DIR}>
  PRIVATE   $<BUILD_INTERFACE:${CMAKE_BINARY_DIR}>
  PUBLIC    $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}>
  PUBLIC    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
)

unset(LINK_LIBRARIES)
list(APPEND LINK_LIBRARIES ${COMMON_LIBRARIES})

if (PDM_ENABLE_SHARED)
  set_property(TARGET ${pdm_io} PROPERTY POSITION_INDEPENDENT_CODE 1)
  add_library(pdm_io_shared SHARED  $<TARGET_OBJECTS:pdm_io>)
  set_target_properties(pdm_io_shared PROPERTIES OUTPUT_NAME "pdm_io")
  set_target_properties(pdm_io_shared PROPERTIES VERSION ${PDM_VERSION})
  target_link_libraries(pdm_io_shared ${LINK_LIBRARIES})
  set_target_properties(pdm_io_shared PROPERTIES
                            PUBLIC_HEADER "${LIST_HEADER}")

  get_target_property(PDM_INC pdm_io INCLUDE_DIRECTORIES)
  target_include_directories(pdm_io_shared PUBLIC ${PDM_INC})
endif()

if (PDM_ENABLE_STATIC)
  set_property(TARGET ${pdm_io} PROPERTY POSITION_INDEPENDENT_CODE 1)
  add_library(pdm_io_static STATIC  $<TARGET_OBJECTS:pdm_io>)
  set_target_properties(pdm_io_static PROPERTIES OUTPUT_NAME "pdm_io")
  set_target_properties(pdm_io_static PROPERTIES VERSION ${PDM_VERSION})
  target_link_libraries(pdm_io_static ${LINK_LIBRARIES})
  set_target_properties(pdm_io_static PROPERTIES
                            PUBLIC_HEADER "${LIST_HEADER}")

  get_target_property(PDM_INC pdm_io INCLUDE_DIRECTORIES)
  target_include_directories(pdm_io_static PUBLIC ${PDM_INC})
endif()

# Spécification de l'installation
# -------------------------------

set (_targets "")

if (PDM_ENABLE_SHARED)
  list (APPEND _targets pdm_io_shared)
endif()

if (PDM_ENABLE_STATIC)
  list (APPEND _targets pdm_io_static)
endif()

foreach(_target IN ITEMS ${_targets})
  install(TARGETS ${_target} EXPORT ${_target}Targets
    RUNTIME DESTINATION bin
    LIBRARY DESTINATION lib
    PUBLIC_HEADER DESTINATION include
    ARCHIVE DESTINATION lib
  )
  install(EXPORT ${_target}Targets
    FILE ${_target}Targets.cmake
    NAMESPACE pdm::
    DESTINATION lib/cmake/${_target}
  )
  add_library(pdm::${_target} ALIAS ${_target})
endforeach()
