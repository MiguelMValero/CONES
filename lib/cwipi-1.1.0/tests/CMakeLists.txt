# ------------
# Build tests
# ------------

if (MPI_C_COMPILER)
    set(CMAKE_C_COMPILER ${MPI_C_COMPILER})
endif()

if (MPI_CXX_COMPILER)
    set(CMAKE_CXX_COMPILER ${MPI_CXX_COMPILER})
endif()

if (MPI_Fortran_COMPILER)
    set(CMAKE_Fortran_COMPILER ${MPI_Fortran_COMPILER})
endif()

unset(LINK_LIBRARIES)
if (CMAKE_C_COMPILER_ID STREQUAL "GNU")
  list(APPEND LINK_LIBRARIES "m")
endif()

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

if (CWP_ENABLE_STATIC)

   if (CWP_ENABLE_Fortran)
     if (MPI_CXX_LIBRARIES)
       list(APPEND LINK_LIBRARIES ${MPI_CXX_LIBRARIES})
     endif()
     list(APPEND LINK_LIBRARIES cwpf_static)
     if (CWP_ENABLE_EXTERNAL_PDM)
       list(APPEND LINK_LIBRARIES pdm::pdmf_static)
     else()
       list(APPEND LINK_LIBRARIES pdmf_static) 
     endif()  
   endif()
   list(APPEND LINK_LIBRARIES cwp_static)
   if (CWP_ENABLE_EXTERNAL_PDM)
     list(APPEND LINK_LIBRARIES
                         pdm::pdm_mpi_static
                         pdm::pdm_static)
   else()
     list(APPEND LINK_LIBRARIES
                         pdm_mpi_static
                         pdm_static)
   endif()

   if (CWP_ENABLE_Fortran)
    list(APPEND LINK_LIBRARIES cwpf_static)
   endif()
   list(APPEND LINK_LIBRARIES ${COMMON_LIBRARIES})

##   if (CWP_ENABLE_SPACE_NODAL_BASIS)
##    list(APPEND LINK_LIBRARIES spaceBasis_static)
##   endif()



else ()
   if (CWP_ENABLE_Fortran AND NOT CWP_NEW_TO_OLD_API)
     list(APPEND LINK_LIBRARIES cwpf_shared)
     if (CWP_ENABLE_EXTERNAL_PDM)
       list(APPEND LINK_LIBRARIES pdm::pdmf_shared)
     else()
       list(APPEND LINK_LIBRARIES pdmf_shared) 
     endif()  
   endif()
   list(APPEND LINK_LIBRARIES cwp_shared)
   if (CWP_ENABLE_EXTERNAL_PDM)
     list(APPEND LINK_LIBRARIES
               pdm::pdm_shared     
               pdm::pdm_mpi_shared)
   else()
     list(APPEND LINK_LIBRARIES
               pdm_shared     
               pdm_mpi_shared)
   endif() 

#   endif ()
   if (CWP_ENABLE_Fortran)
    list(APPEND LINK_LIBRARIES cwpf_shared)
   endif()
   list(APPEND LINK_LIBRARIES ${COMMON_LIBRARIES})

##   if (CWP_ENABLE_SPACE_NODAL_BASIS)
##    list(APPEND LINK_LIBRARIES spaceBasis_shared)
##   endif()

endif()



#if (CWP_ENABLE_Fortran AND NOT CWP_NEW_TO_OLD_API)
#  file(GLOB SOURCESF
#       fortran_surf_PiQj_common.F90)
#
#  add_library(cwp_checkf ${SOURCESF})
#   target_include_directories(cwp_checkf PRIVATE ${CMAKE_SOURCE_DIR}
#                                         PRIVATE ${CMAKE_BINARY_DIR}
#                                         PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
#   if (MPI_Fortran_ADDITIONAL_INCLUDE_DIRS)
#      target_include_directories(cwp_checkf PRIVATE ${MPI_Fortran_ADDITIONAL_INCLUDE_DIRS})
#   endif()
#   target_include_directories(cwp_checkf PRIVATE ${TEST_INC})
#   target_link_libraries(cwp_checkf ${LINK_LIBRARIES})
#   set_target_properties(cwp_checkf PROPERTIES LINKER_LANGUAGE "Fortran")
#endif()

get_target_property(TEST_INC cwp INCLUDE_DIRECTORIES)

file (COPY ${CMAKE_CURRENT_SOURCE_DIR}/meshes DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

# Add test
# --------

if (CWP_ENABLE_PYTHON_BINDINGS AND NOT CWP_NEW_TO_OLD_API)
  SET (CTEST_ENVIRONMENT
       "LD_LIBRARY_PATH=${LD_LIBRARY_PATH}"
       "PYTHONPATH=${CWP_BINARY_DIR}/Cython:${PYTHONPATH}"
       "PATH=${CWP_BINARY_DIR}/bin:${PATH}"
       )
else()
  SET (CTEST_ENVIRONMENT
       "LD_LIBRARY_PATH=${LD_LIBRARY_PATH}"
       "PATH=${CWP_BINARY_DIR}/bin:${PATH}")
endif()

include(TestCreate)

set (_MPIEXEC_MAX_NUMPROCS 2)

  test_c_create (c_new_api 10 CWP_LIST_TEST CWP_LIST_NRANK)
  test_c_create (c_new_api_c_poly_block 2 CWP_LIST_TEST CWP_LIST_NRANK)
  test_c_create (c_new_api_cellface_set_block 2 CWP_LIST_TEST CWP_LIST_NRANK)
  test_c_create (c_new_api_dcube 6 CWP_LIST_TEST CWP_LIST_NRANK)
  test_c_create (c_new_api_disjoint_comms_multiblock 2 CWP_LIST_TEST CWP_LIST_NRANK)
  test_c_create (c_new_api_multiblock 4 CWP_LIST_TEST CWP_LIST_NRANK)
  test_c_create (c_new_api_surf_P1P0_P0P1 2 CWP_LIST_TEST CWP_LIST_NRANK)
  test_c_create (c_new_api_surf_P1P0_P0P1_dynamic 6 CWP_LIST_TEST CWP_LIST_NRANK)
  test_c_create (c_new_api_surf_P1P0_P0P1_without_part 2 CWP_LIST_TEST CWP_LIST_NRANK)
  test_c_create (c_new_api_global_data 4 CWP_LIST_TEST CWP_LIST_NRANK)
  test_c_create (c_new_api_part_data 4 CWP_LIST_TEST CWP_LIST_NRANK)
  test_c_create (c_new_api_part_data2 2 CWP_LIST_TEST CWP_LIST_NRANK)
  test_c_create (c_new_api_part_data3 2 CWP_LIST_TEST CWP_LIST_NRANK)
  test_c_create (c_new_api_surf_P0P0 2 CWP_LIST_TEST CWP_LIST_NRANK)
  test_c_create (c_new_api_surf_P0P0_disjoint_comm 2 CWP_LIST_TEST CWP_LIST_NRANK)
  test_c_create (c_new_api_surf_P0P0_without_part 3 CWP_LIST_TEST CWP_LIST_NRANK)
  test_c_create (c_new_api_guinea_pig 2 CWP_LIST_TEST CWP_LIST_NRANK)
  test_c_create (c_new_vs_old_HO 2 CWP_LIST_TEST CWP_LIST_NRANK)
  test_c_create (c_new_api_wind_turbine_blade 2 CWP_LIST_TEST CWP_LIST_NRANK)
  test_c_create (c_new_api_user_interp_fct 2 CWP_LIST_TEST CWP_LIST_NRANK)
  test_c_create (c_new_api_reverse_knn 2 CWP_LIST_TEST CWP_LIST_NRANK)
  test_c_create (c_new_api_linear_P1P0_P0P1 2 CWP_LIST_TEST CWP_LIST_NRANK)
  test_c_create (c_new_api_multiple 2 CWP_LIST_TEST CWP_LIST_NRANK)
  test_c_create (c_new_api_variable_mesh 2 CWP_LIST_TEST CWP_LIST_NRANK)
  test_c_create (c_new_api_param 2 CWP_LIST_TEST CWP_LIST_NRANK)
  if (CWP_ENABLE_CLIENT_SERVER)
#      test_c_create (c_ip_io 1)
    test_c_create (c_client_server_main 2 CWP_LIST_TEST CWP_LIST_NRANK)
    test_c_create (c_client_server_new_api_cellface_set_block 2 CWP_LIST_TEST CWP_LIST_NRANK)
    test_c_create (c_client_server_new_api_c_poly_block 2 CWP_LIST_TEST CWP_LIST_NRANK)
    test_c_create (c_client_server_new_api_surf_P1P0_P0P1 4 CWP_LIST_TEST CWP_LIST_NRANK)
    test_c_create (c_client_server_new_api_wind_turbine_blade 2 CWP_LIST_TEST CWP_LIST_NRANK)
    test_c_create (c_client_server_new_api_multiple 2 CWP_LIST_TEST CWP_LIST_NRANK)
    test_c_create (c_client_server_new_api_deformable 2 CWP_LIST_TEST CWP_LIST_NRANK)
    test_c_create (c_client_server_new_api_variable_mesh 2 CWP_LIST_TEST CWP_LIST_NRANK)
    add_dependencies(c_client_server_main cwp_server)
    add_dependencies(c_client_server_new_api_cellface_set_block cwp_server)
    add_dependencies(c_client_server_new_api_c_poly_block cwp_server)
    add_dependencies(c_client_server_new_api_surf_P1P0_P0P1 cwp_server)
    add_dependencies(c_client_server_new_api_wind_turbine_blade cwp_server)
    add_dependencies(c_client_server_new_api_multiple cwp_server)
    add_dependencies(c_client_server_new_api_deformable cwp_server)
    add_dependencies(c_client_server_new_api_variable_mesh cwp_server)
  endif()

  if (CWP_ENABLE_Fortran)
    test_fortran_create(fortran_new_api 10 CWP_LIST_TEST CWP_LIST_NRANK)
    test_fortran_create(fortran_new_api_surf 2 CWP_LIST_TEST CWP_LIST_NRANK)
    test_fortran_create(fortran_new_api_surf2 2 CWP_LIST_TEST CWP_LIST_NRANK)
  endif()

  test_c_create (c_new_to_old_api 2 CWP_LIST_TEST CWP_LIST_NRANK)
  test_c_create (c_new_vs_old_vol 2 CWP_LIST_TEST CWP_LIST_NRANK)
  test_c_create (c_new_vs_old_vol_sendrecv 2 CWP_LIST_TEST CWP_LIST_NRANK)
  test_c_create (c_new_to_old_api_surf_P1P1 2 CWP_LIST_TEST CWP_LIST_NRANK)
  # test_c_create (c_new_to_old_api_surf_P1P0_P0P1_poly 2 CWP_LIST_TEST CWP_LIST_NRANK)
  test_c_create (c_new_api_vol 2 CWP_LIST_TEST CWP_LIST_NRANK)

  test_c_create (c_linear_coupling 2 CWP_LIST_TEST CWP_LIST_NRANK)
  test_c_create (c_surf_location_triaP2 2 CWP_LIST_TEST CWP_LIST_NRANK)
  test_c_create (c_surf_location_curvetriaP2 2 CWP_LIST_TEST CWP_LIST_NRANK)
#  test_c_create (c_ho_palm 2)
  test_c_create (c_surf_coupling_P1P1_micron_polygon 2 CWP_LIST_TEST CWP_LIST_NRANK)
  test_c_create (c_surf_coupling_P1P1_polygon 2 CWP_LIST_TEST CWP_LIST_NRANK)
  test_c_create (c_surf_coupling_P1P1_polygon2 2 CWP_LIST_TEST CWP_LIST_NRANK)
  test_c_create (c_surf_coupling_P1P1_polygon_with_external_points 2 CWP_LIST_TEST CWP_LIST_NRANK)
  test_c_create (c_surf_coupling_P1P1 2 CWP_LIST_TEST CWP_LIST_NRANK)
  test_c_create (c_surf_coupling_P1P0_P0P1 2 CWP_LIST_TEST CWP_LIST_NRANK)
  #test_c_create (c_surf_coupling_P1P0_P0P1_surf_gen 2)
  #test_c_create (c_surf_coupling_P1P0_P0P1_curv 2)
  test_c_create (c_surf_coupling_P1P0_P0P1_savedisk 2 CWP_LIST_TEST CWP_LIST_NRANK)
  test_c_create (c_surf_coupling_P1P0_P0P1_savememory 2 CWP_LIST_TEST CWP_LIST_NRANK)
  test_c_create (c_surf_coupling_P1P1_asynchronous 2 CWP_LIST_TEST CWP_LIST_NRANK)
  #test_c_create (c_surf_coupling_P1P1_asynchronous_curv 2)
  test_c_create (c_surf_coupling_user_target_points 2 CWP_LIST_TEST CWP_LIST_NRANK)
  test_c_create (c_surf_coupling_user_interpolation 2 CWP_LIST_TEST CWP_LIST_NRANK)
  #c_surf_coupling_user_interpolation2_SOURCES = c_surf_coupling_user_interpolation2.c
  test_c_create (c_surf_coupling_vector_exchange 2 CWP_LIST_TEST CWP_LIST_NRANK)
  test_c_create (c_surf_coupling_output_error 2 CWP_LIST_TEST CWP_LIST_NRANK)
  test_c_create (c_surf_coupling_simple_location 2 CWP_LIST_TEST CWP_LIST_NRANK)
  test_c_create (c_surf_without_partitioning_coupling 2 CWP_LIST_TEST CWP_LIST_NRANK)
  test_c_create (c_surf_without_partitioning_coupling2 2 CWP_LIST_TEST CWP_LIST_NRANK)
  test_c_create (c_surf_sequential_coupling 2 CWP_LIST_TEST CWP_LIST_NRANK)
  test_c_create (c_vol_coupling_P1P1 2 CWP_LIST_TEST CWP_LIST_NRANK)
  test_c_create (c_vol_coupling_P1P1_with_external_points 2 CWP_LIST_TEST CWP_LIST_NRANK)
  test_c_create (c_vol_coupling_micro_polyhedra_P1P1 2 CWP_LIST_TEST CWP_LIST_NRANK)
  test_c_create (c_vol_coupling_polyhedra_P1P1 2 CWP_LIST_TEST CWP_LIST_NRANK)
  test_c_create (c_vol_coupling_polyhedra_P1P1_with_external_points 2 CWP_LIST_TEST CWP_LIST_NRANK)

if (CWP_ENABLE_Fortran )
  test_fortran_create (fortran_surf_P1P1            2 CWP_LIST_TEST CWP_LIST_NRANK)
  test_fortran_create (fortran_surf_P1P1_savememory 2 CWP_LIST_TEST CWP_LIST_NRANK)
#  test_fortran_create (fortran_surf_PiQj 2 CWP_LIST_TEST CWP_LIST_NRANK)
#  target_link_libraries(fortran_surf_PiQj cwp_checkf)
#  add_dependencies(fortran_surf_PiQj  cwp_checkf)
#  test_fortran_create (fortran_plane_PiQj 2 CWP_LIST_TEST CWP_LIST_NRANK)
#  target_link_libraries(fortran_plane_PiQj cwp_checkf)
#  add_dependencies(fortran_plane_PiQj  cwp_checkf)

#    list(APPEND LINK_LIBRARIES spaceBasis)
#    test_fortran_create  (fortran_surf_TriaP1_PiPj   ${_MPIEXEC_MAX_NUMPROCS})
#    target_link_libraries(fortran_surf_TriaP1_PiPj   spaceBasis)
#    add_dependencies(fortran_surf_TriaP1_PiPj   spaceBasis)

#    test_fortran_create  (fortran_surf_TriaP2_PiPj   ${_MPIEXEC_MAX_NUMPROCS})
#    target_link_libraries(fortran_surf_TriaP2_PiPj   spaceBasis)
#    add_dependencies(fortran_surf_TriaP2_PiPj   spaceBasis)

endif()

if (CWP_ENABLE_PYTHON_BINDINGS)
  test_python_create (python_api 2 CWP_LIST_TEST CWP_LIST_NRANK)
  test_python_create (python_new_api 2 CWP_LIST_TEST CWP_LIST_NRANK)
  test_python_create (python_new_api_surf 2 CWP_LIST_TEST CWP_LIST_NRANK)
  test_python_create (python_new_api_multiple_field_with_callback 2 CWP_LIST_TEST CWP_LIST_NRANK)
  test_python_create (python_new_api_multiblock_2d 2 CWP_LIST_TEST CWP_LIST_NRANK)
  test_python_create (python_new_api_multiblock_3d 2 CWP_LIST_TEST CWP_LIST_NRANK)
  test_python_create (python_new_api_multiple_orders 2 CWP_LIST_TEST CWP_LIST_NRANK)
  if (CWP_ENABLE_CLIENT_SERVER)
    test_python_create (python_new_api_client_server 2 CWP_LIST_TEST CWP_LIST_NRANK)
    add_dependencies(python_new_api_client_server cwp_server)
  endif()
endif()

add_subdirectory(tutorial)

set (CWP_LIST_TEST ${CWP_LIST_TEST} PARENT_SCOPE )
set (CWP_LIST_NRANK ${CWP_LIST_NRANK} PARENT_SCOPE )
