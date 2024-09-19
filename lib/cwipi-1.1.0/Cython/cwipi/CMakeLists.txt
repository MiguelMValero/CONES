file(GLOB SOURCE_PXI
     CONFIGURE_DEPENDS 
     *.pxi)

file(GLOB SOURCE_PY
     CONFIGURE_DEPENDS 
     *.py)

foreach(F IN LISTS SOURCE_PXI SOURCE_PY)
    get_filename_component(F_BASENAME "${F}" NAME )    
    configure_file(${F_BASENAME} ${CMAKE_CURRENT_BINARY_DIR}/${F_BASENAME} ) 
endforeach()

configure_file(cwipi.pyx.in cwipi.pyx)

include(UseCython)

cython_add_module (cwipi
                   PYX_SOURCES "${CMAKE_CURRENT_BINARY_DIR}/cwipi.pyx"
                   INCLUDE_DIRECTORIES ${include_dirs}
                   DEPENDS ${SOURCE_PXI} cwp_shared)

set_target_properties(cwipi PROPERTIES PREFIX "") # do not prefix by lib

get_target_property(CWP_INC cwipi INCLUDE_DIRECTORIES)
target_include_directories(cwipi PUBLIC ${CWP_INC})

target_link_libraries(cwipi Python::NumPy MPI::MPI_CXX)

target_link_libraries(cwipi ${LINK_LIBRARIES} pdm::pdm_shared pdm::pdm_mpi_shared)

target_link_libraries(cwipi cwp_shared)
if (CWP_ENABLE_Fortran)
  target_link_libraries(cwipi cwpf_shared)
endif()

if (CWP_ENABLE_SPACE_NODAL_BASIS)
  target_link_libraries(cwipi spaceBasis_shared)
endif()

if ((NOT MPI_C_COMPILER) AND MPI_C_COMPILE_FLAGS)
   set_target_properties(cwipi PROPERTIES
                         COMPILE_FLAGS ${MPI_C_COMPILE_FLAGS})
endif()


#set(SETUP_PY_OLD "${CMAKE_CURRENT_BINARY_DIR}/setup.py")
#set(OUTPUT       "${CMAKE_CURRENT_BINARY_DIR}/timestamp")
#set(DEPS_OLD    "${CMAKE_CURRENT_SOURCE_DIR}/__init__.py")
#add_custom_command(OUTPUT ${OUTPUT}
#                   COMMAND ${Python_EXECUTABLE} ${SETUP_PY_OLD} --quiet build
#                   COMMAND ${CMAKE_COMMAND} -E touch ${OUTPUT}
#                   DEPENDS ${DEPS_OLD})



set (SITE_PACKAGES_CWP_OLD "${CMAKE_CWP_INSTALL_PYTHON_DIR}/lib/python${Python_VERSION_MAJOR}.${Python_VERSION_MINOR}/site-packages/cwipi")

if ((CMAKE_C_COMPILER_ID STREQUAL "GNU") OR
    (CMAKE_C_COMPILER_ID STREQUAL "Intel") OR
    (CMAKE_C_COMPILER_ID STREQUAL "AppleClang") OR
    (CMAKE_C_COMPILER_ID STREQUAL "IntelLLVM")  OR
    (CMAKE_C_COMPILER_ID STREQUAL "PGI"))
  set_target_properties (cwipi PROPERTIES COMPILE_FLAGS "-w")
endif()


install(DIRECTORY DESTINATION ${SITE_PACKAGES_CWP_OLD})
install(FILES  "${CMAKE_CURRENT_SOURCE_DIR}/__init__.py"
        DESTINATION            "${SITE_PACKAGES_CWP_OLD}")
install(PROGRAMS  "${CMAKE_CURRENT_BINARY_DIR}/cwipi.so"
        DESTINATION            "${SITE_PACKAGES_CWP_OLD}")
