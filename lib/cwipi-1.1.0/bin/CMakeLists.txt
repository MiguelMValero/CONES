# ------------
# Build server
# ------------

if (MPI_C_COMPILER)
    set(CMAKE_C_COMPILER ${MPI_C_COMPILER})
endif()

if (MPI_CXX_COMPILER)
    set(CMAKE_CXX_COMPILER ${MPI_CXX_COMPILER})
endif()

unset(LINK_LIBRARIES)

if (NOT MPI_C_COMPILER)
  if (MPI_C_LIBRARIES)
    list(APPEND LINK_LIBRARIES ${MPI_C_LIBRARIES})
  endif()
endif ()

if (NOT MPI_CXX_COMPILER)
  if (MPI_CXX_LIBRARIES)
    list(APPEND LINK_LIBRARIES ${MPI_C_LIBRARIES})
  else ()
    if (MPI_C_LIBRARIES)
     list(APPEND LINK_LIBRARIES ${MPI_C_LIBRARIES})
    endif()
  endif()
endif ()

if (NOT CWP_ENABLE_EXTERNAL_PDM)

    if (CWP_ENABLE_STATIC)
       list(APPEND LINK_LIBRARIES cwp_static)
       list(APPEND LINK_LIBRARIES
              pdm_mpi_static
              pdm_static)
       list(APPEND LINK_LIBRARIES ${COMMON_LIBRARIES})
    else ()
       list(APPEND LINK_LIBRARIES cwp_shared)
       list(APPEND LINK_LIBRARIES
              pdm_mpi_shared
              pdm_shared)
       list(APPEND LINK_LIBRARIES ${COMMON_LIBRARIES})
    endif()

else ()
    if (CWP_ENABLE_STATIC)
       list(APPEND LINK_LIBRARIES cwp_static)
       list(APPEND LINK_LIBRARIES
              pdm::pdm_mpi_static
              pdm::pdm_static)
       list(APPEND LINK_LIBRARIES ${COMMON_LIBRARIES})
    else ()
       list(APPEND LINK_LIBRARIES cwp_shared)
       list(APPEND LINK_LIBRARIES
              pdm::pdm_mpi_shared
              pdm::pdm_shared)
       list(APPEND LINK_LIBRARIES ${COMMON_LIBRARIES})
    endif()

endif()

# add executable

add_executable (cwp_server server_main.c)
target_include_directories(cwp_server PRIVATE ${CMAKE_SOURCE_DIR}
                           PRIVATE ${CMAKE_BINARY_DIR}
                           PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})

target_link_libraries(cwp_server ${LINK_LIBRARIES})
install(TARGETS cwp_server DESTINATION bin)
