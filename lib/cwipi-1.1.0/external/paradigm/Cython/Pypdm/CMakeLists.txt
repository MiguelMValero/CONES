
if (PDM_ENABLE_EXTENSION_PDMA)
  file(GLOB SOURCE_PXI
       CONFIGURE_DEPENDS 
       *.pxi
       ${PDM_SOURCE_DIR}/extensions/paradigma/Cython/*.pxi)

  file(GLOB SOURCE_PY
       CONFIGURE_DEPENDS 
       *.py
       ${PDM_SOURCE_DIR}/extensions/paradigma/Cython/*.py)
else ()
  file(GLOB SOURCE_PXI
       CONFIGURE_DEPENDS 
       *.pxi)

  file(GLOB SOURCE_PY
       CONFIGURE_DEPENDS 
       *.py)

endif()

foreach(F IN LISTS SOURCE_PXI SOURCE_PY)
    get_filename_component(F_BASENAME "${F}" NAME )    
    configure_file(${F} ${CMAKE_CURRENT_BINARY_DIR}/${F_BASENAME} ) 
endforeach()

configure_file(Pypdm.pyx.in Pypdm.pyx)
configure_file(pdm_part_coarse.pxi.in pdm_part_coarse.pxi)

include(UseCython)
set(_pyx_file "${CMAKE_CURRENT_BINARY_DIR}/Pypdm.pyx")
set_source_files_properties(${_pyx_file} PROPERTIES CYTHON_IS_CXX FALSE)
# set_source_files_properties(${_pyx_file} PROPERTIES CYTHON_IS_CXX TRUE)

cython_add_module (Pypdm
                   PYX_SOURCES "${_pyx_file}"
                   INCLUDE_DIRECTORIES ${include_dirs}
                   DEPENDS ${SOURCE_PXI} pdm_shared)

set_target_properties(Pypdm PROPERTIES PREFIX "") # do not prefix by lib

set_target_properties(Pypdm PROPERTIES LIBRARY_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/")

unset(LINK_LIBRARIES)
list(APPEND LINK_LIBRARIES ${COMMON_LIBRARIES})
if (CMAKE_C_COMPILER_ID STREQUAL "GNU")
  list(APPEND LINK_LIBRARIES "m")
endif()
if (${PARMETIS_FOUND})
  list(APPEND LINK_LIBRARIES ${PARMETIS_LIBRARIES})
  list(APPEND LINK_LIBRARIES ${METIS_LIBRARIES})
endif ()
if (${PTSCOTCH_FOUND})
  list(APPEND LINK_LIBRARIES ${PTSCOTCH_LIBRARIES})
  list(APPEND LINK_LIBRARIES ${SCOTCH_LIBRARIES})
endif ()

if (PDM_ENABLE_EXTENSION_PDMA)
  list(APPEND LINK_LIBRARIES pdma_shared)
endif()

#find_package(MPI REQUIRED COMPONENTS CXX)
target_link_libraries(Pypdm ${LINK_LIBRARIES} pdm_shared pdm_mpi_shared Python::NumPy MPI::MPI_CXX)

get_target_property(PDM_INC pdm INCLUDE_DIRECTORIES)
target_include_directories(Pypdm PUBLIC ${PDM_INC})

if ((NOT MPI_C_COMPILER) AND MPI_C_COMPILE_FLAGS)
   set_target_properties(Pypdm PROPERTIES
                              COMPILE_FLAGS ${MPI_C_COMPILE_FLAGS})
endif()
if ((NOT MPI_CXX_COMPILER) AND MPI_CXX_COMPILE_FLAGS)
   set_target_properties(Pypdm PROPERTIES
                              COMPILE_FLAGS ${MPI_CXX_COMPILE_FLAGS})
endif()

#set(SETUP_PY    "${CMAKE_CURRENT_BINARY_DIR}/setup.py")
#set(DEPS        "${CMAKE_CURRENT_SOURCE_DIR}/__init__.py")
#set(OUTPUT      "${CMAKE_CURRENT_BINARY_DIR}/timestamp")

#add_custom_command(OUTPUT ${OUTPUT}
#                       COMMAND ${Python_EXECUTABLE} ${SETUP_PY}  build
#                       COMMAND ${CMAKE_COMMAND} -E touch ${OUTPUT}
#                       DEPENDS ${DEPS})

set (SITE_PACKAGES_PDM "${CMAKE_INSTALL_PREFIX}/lib/python${Python_VERSION_MAJOR}.${Python_VERSION_MINOR}/site-packages/Pypdm")

#add_custom_target(target ALL DEPENDS ${OUTPUT})


if ((CMAKE_C_COMPILER_ID STREQUAL "GNU") OR
    (CMAKE_C_COMPILER_ID STREQUAL "Intel") OR
    (CMAKE_C_COMPILER_ID STREQUAL "AppleClang")  OR
    (CMAKE_C_COMPILER_ID STREQUAL "IntelLLVM")  OR
    (CMAKE_C_COMPILER_ID STREQUAL "PGI"))
  set_target_properties (Pypdm PROPERTIES COMPILE_FLAGS "-w")
endif()

install(DIRECTORY DESTINATION ${SITE_PACKAGES_PDM})
install(FILES  "${CMAKE_CURRENT_SOURCE_DIR}/__init__.py"
        DESTINATION            "${SITE_PACKAGES_PDM}")

install(PROGRAMS  "${CMAKE_CURRENT_BINARY_DIR}/Pypdm.so"
        DESTINATION            "${SITE_PACKAGES_PDM}")
