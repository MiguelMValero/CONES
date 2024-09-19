
add_subdirectory(io)
add_subdirectory(mpi_wrapper)


# Liste de tous les sources
# -------------------------

file(GLOB_RECURSE SOURCES
     mesh/*.[ch]
     meshgen/*.[ch]
     pario/*.[ch]
     ppart/*.[ch]
     ppart/*.cxx
     struct/*.[ch]
     struct/*.cxx
     util/*.[ch]
     ext_wrapper/*.[ch]
     mpi_wrapper/pdm_mpi_node_first_rank.[hc]
     io/pdm_logging.[hc]
     io/pdm_vtk.[hc]
     pdm.h
     pdm.c
     pdm_priv.h)

file(GLOB_RECURSE SOURCESF
     *.[fF]
     *.[fF]90)

file(GLOB_RECURSE SOURCESGPU
     gpu/*.cu
     gpu/*.cuh)

add_library(pdm OBJECT ${SOURCES})

if (PDM_ENABLE_STATIC_Fortran_INTERFACE OR PDM_ENABLE_SHARED_Fortran_INTERFACE)
  add_library(pdmf OBJECT ${SOURCESF})
endif()

if (PDM_ENABLE_CUDA)
  add_library(pdm_gpu STATIC ${SOURCESGPU})
  set_target_properties(pdm_gpu PROPERTIES
                            CUDA_SEPARABLE_COMPILATION ON)
  target_compile_options(pdm_gpu PRIVATE $<$<COMPILE_LANGUAGE:CUDA>:
  -maxrregcount=32;--fmad=false;#-G;-g;
  # -gencode=arch=compute_35,code=sm_35;
  # -gencode=arch=compute_37,code=sm_37;
  # -gencode=arch=compute_50,code=sm_50;
  # -gencode=arch=compute_52,code=sm_52;
  # -gencode=arch=compute_60,code=sm_60;
  -gencode=arch=compute_61,code=sm_61;
  # -gencode=arch=compute_70,code=sm_70;
  # -gencode=arch=compute_75,code=sm_75;
  # -gencode=arch=compute_75,code=compute_75;
  >)
  #set(CUDA_PROPAGATE_HOST_FLAGS OFF)
  #list(APPEND CUDA_NVCC_FLAGS -gencode arch=compute_61,code=sm_51;)
  #set(CUDA_NVCC_FLAGS "-gencode arch=compute_60,code=sm_60")
  #set(CUDA_LAUNCH_BLOCKING 1)
  #set(CMAKE_C_FLAGS "-std=c17 -pedantic-errors ${CMAKE_C_FLAGS}")
endif()

target_include_directories(pdm
  INTERFACE $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}>
  PRIVATE   $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}>
  PRIVATE   $<BUILD_INTERFACE:${CMAKE_BINARY_DIR}>
  PRIVATE   $<BUILD_INTERFACE:${PDM_BINARY_DIR}>
  PUBLIC    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/..>
  PUBLIC    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
  PUBLIC    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/pario>
  PUBLIC    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/ppart>
  PUBLIC    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/io>
  PUBLIC    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/mpi_wrapper>
  PUBLIC    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/ext_wrapper>
  PUBLIC    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/mesh>
  PUBLIC    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/meshgen>
  PUBLIC    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/struct>
  PUBLIC    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/gpu>
  PUBLIC    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/adapt>
  PUBLIC    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/util>
  PUBLIC    $<INSTALL_INTERFACE:include/>
)

if (PDM_ENABLE_EXTENSION_PDMA)
  target_include_directories(pdm
                             PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/../extensions/paradigma/src/plugins/renum>
                             )
endif()

if (PDM_ENABLE_STATIC_Fortran_INTERFACE OR PDM_ENABLE_SHARED_Fortran_INTERFACE)
  target_include_directories(pdmf
    INTERFACE $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}>
    PRIVATE   $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}>
    PRIVATE   $<BUILD_INTERFACE:${CMAKE_BINARY_DIR}>
    PRIVATE   $<BUILD_INTERFACE:${PDM_BINARY_DIR}>
    PUBLIC    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/..>
    PUBLIC    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
    PUBLIC    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/pario>
    PUBLIC    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/ppart>
    PUBLIC    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/io>
    PUBLIC    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/mpi_wrapper>
    PUBLIC    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/ext_wrapper>
    PUBLIC    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/mesh>
    PUBLIC    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/meshgen>
    PUBLIC    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/struct>
    PUBLIC    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/util>
    PUBLIC    $<INSTALL_INTERFACE:include/>
  )
    if (PDM_ENABLE_EXTENSION_PDMA)
      target_include_directories(pdmf
                                 PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/../extensions/paradigma/src/plugins/renum>
                                 )
    endif()
endif()

if (PTSCOTCH_FOUND)
    target_include_directories(pdm
                              PRIVATE $<BUILD_INTERFACE:${PTSCOTCH_INCLUDE_DIRS}>)
endif (PTSCOTCH_FOUND)

if (PARMETIS_FOUND)
    target_include_directories(pdm
                               PRIVATE $<BUILD_INTERFACE:${PARMETIS_INCLUDE_DIRS}>)
    target_include_directories(pdm
                               PRIVATE $<BUILD_INTERFACE:${METIS_INCLUDE_DIRS}>)
endif (PARMETIS_FOUND)

if (PDM_ENABLE_CUDA)
    target_include_directories(pdm_gpu
                              INTERFACE ${CMAKE_CURRENT_BINARY_DIR}
                              PRIVATE   ${CMAKE_CURRENT_BINARY_DIR}
                              PRIVATE  ${CMAKE_BINARY_DIR}
                              PRIVATE  ${PDM_BINARY_DIR}
                              PUBLIC  ${CMAKE_CURRENT_SOURCE_DIR}/..
                              PUBLIC  ${CMAKE_CURRENT_SOURCE_DIR}
                              #PUBLIC  /opt/tools/openmpi/4.0.1-gnu-8.3/include/
                              PUBLIC  ${MPI_INCLUDE_PATH}
                              PUBLIC  ${CUDA_INCLUDE_DIRS}
                              PUBLIC  ${CMAKE_CURRENT_SOURCE_DIR}/pario
                              PUBLIC  ${CMAKE_CURRENT_SOURCE_DIR}/ppart
                              PUBLIC  ${CMAKE_CURRENT_SOURCE_DIR}/io
                              PUBLIC  ${CMAKE_CURRENT_SOURCE_DIR}/mpi_wrapper
                              PUBLIC  ${CMAKE_CURRENT_SOURCE_DIR}/ext_wrapper
                              PUBLIC  ${CMAKE_CURRENT_SOURCE_DIR}/mesh
                              PUBLIC  ${CMAKE_CURRENT_SOURCE_DIR}/meshgen
                              PUBLIC  ${CMAKE_CURRENT_SOURCE_DIR}/struct
                              PUBLIC  ${CMAKE_CURRENT_SOURCE_DIR}/util)
endif ()

# Liste des header
# ----------------

set (LIST_HEADER "")
list(APPEND LIST_HEADER "${CMAKE_CURRENT_SOURCE_DIR}/pdm.h")
list(APPEND LIST_HEADER "${CMAKE_CURRENT_SOURCE_DIR}/pdm_priv.h")
list(APPEND LIST_HEADER "${CMAKE_CURRENT_SOURCE_DIR}/util/pdm_version.h")
list(APPEND LIST_HEADER "${CMAKE_CURRENT_SOURCE_DIR}/util/pdm_array.h")
list(APPEND LIST_HEADER "${CMAKE_CURRENT_SOURCE_DIR}/pario/pdm_io.h")
list(APPEND LIST_HEADER "${CMAKE_CURRENT_SOURCE_DIR}/io/pdm_error.h")
list(APPEND LIST_HEADER "${CMAKE_CURRENT_SOURCE_DIR}/io/pdm_printf.h")
list(APPEND LIST_HEADER "${CMAKE_CURRENT_SOURCE_DIR}/io/pdm_logging.h")
list(APPEND LIST_HEADER "${CMAKE_CURRENT_SOURCE_DIR}/ppart/pdm_part_coarse_mesh.h")
list(APPEND LIST_HEADER "${CMAKE_CURRENT_SOURCE_DIR}/ppart/pdm_part.h")
list(APPEND LIST_HEADER "${CMAKE_CURRENT_SOURCE_DIR}/ppart/pdm_part_renum.h")
list(APPEND LIST_HEADER "${CMAKE_CURRENT_SOURCE_DIR}/ppart/pdm_partitioning_algorithm.h")
list(APPEND LIST_HEADER "${CMAKE_CURRENT_SOURCE_DIR}/ppart/pdm_multipart.h")
list(APPEND LIST_HEADER "${CMAKE_CURRENT_SOURCE_DIR}/ppart/pdm_part_extension.h")
list(APPEND LIST_HEADER "${CMAKE_CURRENT_SOURCE_DIR}/pario/pdm_io_tab.h")
list(APPEND LIST_HEADER "${CMAKE_CURRENT_SOURCE_DIR}/pario/pdm_writer.h")
list(APPEND LIST_HEADER "${CMAKE_CURRENT_SOURCE_DIR}/mesh/pdm_mesh_nodal.h")
list(APPEND LIST_HEADER "${CMAKE_CURRENT_SOURCE_DIR}/mesh/pdm_mesh_intersection.h")
list(APPEND LIST_HEADER "${CMAKE_CURRENT_SOURCE_DIR}/mesh/pdm_mesh_location.h")
list(APPEND LIST_HEADER "${CMAKE_CURRENT_SOURCE_DIR}/mesh/pdm_part_mesh_nodal.h")
list(APPEND LIST_HEADER "${CMAKE_CURRENT_SOURCE_DIR}/mesh/pdm_ho_ordering.h")
list(APPEND LIST_HEADER "${CMAKE_CURRENT_SOURCE_DIR}/mesh/pdm_closest_points.h")
list(APPEND LIST_HEADER "${CMAKE_CURRENT_SOURCE_DIR}/mesh/pdm_domain_interface.h")
list(APPEND LIST_HEADER "${CMAKE_CURRENT_SOURCE_DIR}/mesh/pdm_part_domain_interface.h")
list(APPEND LIST_HEADER "${CMAKE_CURRENT_SOURCE_DIR}/struct/pdm_dmesh.h")
list(APPEND LIST_HEADER "${CMAKE_CURRENT_SOURCE_DIR}/struct/pdm_part_mesh.h")
list(APPEND LIST_HEADER "${CMAKE_CURRENT_SOURCE_DIR}/pario/pdm_writer_priv.h")
list(APPEND LIST_HEADER "${CMAKE_CURRENT_SOURCE_DIR}/util/pdm_mpi_node_first_rank.h")
list(APPEND LIST_HEADER "${CMAKE_CURRENT_SOURCE_DIR}/util/pdm_timer.h")
list(APPEND LIST_HEADER "${CMAKE_CURRENT_SOURCE_DIR}/util/pdm_linear_algebra.h")
list(APPEND LIST_HEADER "${CMAKE_CURRENT_SOURCE_DIR}/util/pdm_fortran_to_c_string.h")
list(APPEND LIST_HEADER "${CMAKE_CURRENT_SOURCE_DIR}/meshgen/pdm_dcube_gen.h")
list(APPEND LIST_HEADER "${CMAKE_CURRENT_SOURCE_DIR}/meshgen/pdm_dcube_nodal_gen.h")
list(APPEND LIST_HEADER "${CMAKE_CURRENT_SOURCE_DIR}/meshgen/pdm_poly_surf_gen.h")
list(APPEND LIST_HEADER "${CMAKE_CURRENT_SOURCE_DIR}/meshgen/pdm_sphere_surf_gen.h")
list(APPEND LIST_HEADER "${CMAKE_CURRENT_SOURCE_DIR}/meshgen/pdm_generate_mesh.h")
list(APPEND LIST_HEADER "${CMAKE_CURRENT_SOURCE_DIR}/meshgen/pdm_poly_vol_gen.h")
list(APPEND LIST_HEADER "${CMAKE_CURRENT_SOURCE_DIR}/struct/pdm_gnum.h")
list(APPEND LIST_HEADER "${CMAKE_CURRENT_SOURCE_DIR}/struct/pdm_gnum_location.h")
list(APPEND LIST_HEADER "${CMAKE_CURRENT_SOURCE_DIR}/struct/pdm_distrib.h")
list(APPEND LIST_HEADER "${CMAKE_CURRENT_SOURCE_DIR}/struct/pdm_gnum_from_hash_values.h")
list(APPEND LIST_HEADER "${CMAKE_CURRENT_SOURCE_DIR}/struct/pdm_part_to_block.h")
list(APPEND LIST_HEADER "${CMAKE_CURRENT_SOURCE_DIR}/struct/pdm_block_to_part.h")
list(APPEND LIST_HEADER "${CMAKE_CURRENT_SOURCE_DIR}/struct/pdm_multi_block_to_part.h")
list(APPEND LIST_HEADER "${CMAKE_CURRENT_SOURCE_DIR}/struct/pdm_part_to_part.h")
list(APPEND LIST_HEADER "${CMAKE_CURRENT_SOURCE_DIR}/struct/pdm_interpolate_from_mesh_location.h")
list(APPEND LIST_HEADER "${CMAKE_CURRENT_SOURCE_DIR}/struct/pdm_hash_tab.h")
list(APPEND LIST_HEADER "${CMAKE_CURRENT_SOURCE_DIR}/struct/pdm_handles.h")
list(APPEND LIST_HEADER "${CMAKE_CURRENT_SOURCE_DIR}/ppart/pdm_part_geom.h")
list(APPEND LIST_HEADER "${CMAKE_CURRENT_SOURCE_DIR}/struct/pdm_sort.h")
list(APPEND LIST_HEADER "${CMAKE_CURRENT_SOURCE_DIR}/struct/pdm_binary_search.h")
list(APPEND LIST_HEADER "${CMAKE_CURRENT_SOURCE_DIR}/struct/pdm_global_mean.h")
list(APPEND LIST_HEADER "${CMAKE_CURRENT_SOURCE_DIR}/struct/pdm_radix_sort.h")
list(APPEND LIST_HEADER "${CMAKE_CURRENT_SOURCE_DIR}/struct/pdm_extract_part.h")
list(APPEND LIST_HEADER "${CMAKE_CURRENT_SOURCE_DIR}/mesh/pdm_cellface_orient.h")
list(APPEND LIST_HEADER "${CMAKE_CURRENT_SOURCE_DIR}/mesh/pdm_dmesh_nodal.h")
list(APPEND LIST_HEADER "${CMAKE_CURRENT_SOURCE_DIR}/mesh/pdm_dmesh_nodal_elmts.h")
list(APPEND LIST_HEADER "${CMAKE_CURRENT_SOURCE_DIR}/mesh/pdm_part_mesh_nodal_elmts.h")
list(APPEND LIST_HEADER "${CMAKE_CURRENT_SOURCE_DIR}/mesh/pdm_surf_mesh.h")
list(APPEND LIST_HEADER "${CMAKE_CURRENT_SOURCE_DIR}/mesh/pdm_dist_cloud_surf.h")
list(APPEND LIST_HEADER "${CMAKE_CURRENT_SOURCE_DIR}/mesh/pdm_dmesh_nodal_to_dmesh.h")
list(APPEND LIST_HEADER "${CMAKE_CURRENT_SOURCE_DIR}/mesh/pdm_dconnectivity_transform.h")
list(APPEND LIST_HEADER "${CMAKE_CURRENT_SOURCE_DIR}/mesh/pdm_part_connectivity_transform.h")
list(APPEND LIST_HEADER "${CMAKE_CURRENT_SOURCE_DIR}/mesh/pdm_overlay.h")
list(APPEND LIST_HEADER "${CMAKE_CURRENT_SOURCE_DIR}/mesh/pdm_geom_elem.h")
list(APPEND LIST_HEADER "${CMAKE_CURRENT_SOURCE_DIR}/mesh/pdm_edges_intersect.h")
list(APPEND LIST_HEADER "${CMAKE_CURRENT_SOURCE_DIR}/mesh/pdm_overlayf.h")
list(APPEND LIST_HEADER "${CMAKE_CURRENT_SOURCE_DIR}/mesh/pdm_poly_clipp.h")
list(APPEND LIST_HEADER "${CMAKE_CURRENT_SOURCE_DIR}/io/pdm_vtk.h")
list(APPEND LIST_HEADER "${PDM_BINARY_DIR}/pdm_config.h")

if (PDM_ENABLE_STATIC_Fortran_INTERFACE OR PDM_ENABLE_SHARED_Fortran_INTERFACE)

  set (LIST_HEADERF "")
  list(APPEND LIST_HEADERF "${CMAKE_CURRENT_SOURCE_DIR}/pario/pdm_io.F90")
  list(APPEND LIST_HEADERF "${CMAKE_CURRENT_SOURCE_DIR}/struct/pdm_cf_array.F90")
  list(APPEND LIST_HEADERF "${CMAKE_CURRENT_BINARY_DIR}/pdm_cf_array.mod")
  list(APPEND LIST_HEADERF "${CMAKE_CURRENT_SOURCE_DIR}/pario/pdm_writer.F90")
  list(APPEND LIST_HEADERF "${CMAKE_CURRENT_SOURCE_DIR}/pdm.F90")
  list(APPEND LIST_HEADERF "${CMAKE_CURRENT_SOURCE_DIR}/ppart/pdm_multipart.F90")
  list(APPEND LIST_HEADERF "${CMAKE_CURRENT_BINARY_DIR}/pdm_multipart.mod")
  list(APPEND LIST_HEADERF "${CMAKE_CURRENT_SOURCE_DIR}/ppart/pdm_part.F90")
  list(APPEND LIST_HEADERF "${CMAKE_CURRENT_SOURCE_DIR}/ppart/pdm_partf.h")
  list(APPEND LIST_HEADERF "${CMAKE_CURRENT_SOURCE_DIR}/pdmf.h")
  list(APPEND LIST_HEADERF "${PDM_BINARY_DIR}/pdm_configf.h")
  list(APPEND LIST_HEADERF "${CMAKE_CURRENT_BINARY_DIR}/pdm_io.mod")
  list(APPEND LIST_HEADERF "${CMAKE_CURRENT_BINARY_DIR}/pdm_fortran.mod")
  list(APPEND LIST_HEADERF "${CMAKE_CURRENT_BINARY_DIR}/pdm_writer.mod")
  list(APPEND LIST_HEADERF "${CMAKE_CURRENT_BINARY_DIR}/pdm.mod")
  list(APPEND LIST_HEADERF "${CMAKE_CURRENT_BINARY_DIR}/pdm_gnum_location.mod")
  list(APPEND LIST_HEADERF "${CMAKE_CURRENT_BINARY_DIR}/pdm_part.mod")
  list(APPEND LIST_HEADERF "${CMAKE_CURRENT_BINARY_DIR}/pdm_mesh_nodal.mod")
  list(APPEND LIST_HEADERF "${CMAKE_CURRENT_SOURCE_DIR}/mesh/pdm_mesh_nodal.F90")
  list(APPEND LIST_HEADERF "${CMAKE_CURRENT_BINARY_DIR}/pdm_part_mesh_nodal.mod")
  list(APPEND LIST_HEADERF "${CMAKE_CURRENT_SOURCE_DIR}/mesh/pdm_part_mesh_nodal.F90")
  list(APPEND LIST_HEADERF "${CMAKE_CURRENT_BINARY_DIR}/pdm_dist_cloud_surf.mod")
  list(APPEND LIST_HEADERF "${CMAKE_CURRENT_SOURCE_DIR}/mesh/pdm_dist_cloud_surf.F90")
  list(APPEND LIST_HEADERF "${CMAKE_CURRENT_BINARY_DIR}/pdm_part_to_block.mod")
  list(APPEND LIST_HEADERF "${CMAKE_CURRENT_SOURCE_DIR}/struct/pdm_part_to_block.F90")
  list(APPEND LIST_HEADERF "${CMAKE_CURRENT_BINARY_DIR}/pdm_part_to_part.mod")
  list(APPEND LIST_HEADERF "${CMAKE_CURRENT_SOURCE_DIR}/struct/pdm_part_to_part.F90")
  list(APPEND LIST_HEADERF "${CMAKE_CURRENT_BINARY_DIR}/pdm_gnum.mod")
  list(APPEND LIST_HEADERF "${CMAKE_CURRENT_SOURCE_DIR}/struct/pdm_gnum.F90")
  list(APPEND LIST_HEADERF "${CMAKE_CURRENT_BINARY_DIR}/pdm_extract_part.mod")
  list(APPEND LIST_HEADERF "${CMAKE_CURRENT_SOURCE_DIR}/struct/pdm_extract_part.F90")
  list(APPEND LIST_HEADERF "${CMAKE_CURRENT_BINARY_DIR}/pdm_block_to_part.mod")
  list(APPEND LIST_HEADERF "${CMAKE_CURRENT_SOURCE_DIR}/struct/pdm_block_to_part.F90")
  list(APPEND LIST_HEADERF "${CMAKE_CURRENT_BINARY_DIR}/pdm_global_mean.mod")
  list(APPEND LIST_HEADERF "${CMAKE_CURRENT_SOURCE_DIR}/struct/pdm_global_mean.F90")
  list(APPEND LIST_HEADERF "${CMAKE_CURRENT_BINARY_DIR}/pdm_mesh_location.mod")
  list(APPEND LIST_HEADERF "${CMAKE_CURRENT_SOURCE_DIR}/mesh/pdm_mesh_location.F90")
  list(APPEND LIST_HEADERF "${CMAKE_CURRENT_BINARY_DIR}/pdm_mesh_intersection.mod")
  list(APPEND LIST_HEADERF "${CMAKE_CURRENT_SOURCE_DIR}/mesh/pdm_mesh_intersection.F90")
  list(APPEND LIST_HEADERF "${CMAKE_CURRENT_BINARY_DIR}/pdm_dcube_gen.mod")
  list(APPEND LIST_HEADERF "${CMAKE_CURRENT_SOURCE_DIR}/meshgen/pdm_dcube_gen.F90")
  list(APPEND LIST_HEADERF "${CMAKE_CURRENT_BINARY_DIR}/pdm_closest_points.mod")
  list(APPEND LIST_HEADERF "${CMAKE_CURRENT_SOURCE_DIR}/mesh/pdm_closest_points.F90")
  list(APPEND LIST_HEADERF "${CMAKE_CURRENT_BINARY_DIR}/pdm_overlay.mod")
  list(APPEND LIST_HEADERF "${CMAKE_CURRENT_SOURCE_DIR}/mesh/pdm_overlay.F90")
  list(APPEND LIST_HEADERF "${CMAKE_CURRENT_BINARY_DIR}/pdm_mesh_check.mod")
  list(APPEND LIST_HEADERF "${CMAKE_CURRENT_BINARY_DIR}/pdm_pointer_array.mod")
  list(APPEND LIST_HEADERF "${CMAKE_CURRENT_SOURCE_DIR}/mesh/pdm_mesh_check.F90")
  list(APPEND LIST_HEADERF "${CMAKE_CURRENT_BINARY_DIR}/pdm_sphere_surf_gen.mod")
  list(APPEND LIST_HEADERF "${CMAKE_CURRENT_SOURCE_DIR}/meshgen/pdm_sphere_surf_gen.F90")
  list(APPEND LIST_HEADERF "${CMAKE_CURRENT_BINARY_DIR}/pdm_part_extension.mod")
  list(APPEND LIST_HEADERF "${CMAKE_CURRENT_SOURCE_DIR}/ppart/pdm_part_extension.F90")
  list(APPEND LIST_HEADERF "${CMAKE_CURRENT_BINARY_DIR}/pdm_global_reduce.mod")
  list(APPEND LIST_HEADERF "${CMAKE_CURRENT_SOURCE_DIR}/struct/pdm_global_reduce.F90")
  list(APPEND LIST_HEADERF "${CMAKE_CURRENT_SOURCE_DIR}/meshgen/pdm_generate_mesh.F90")
  list(APPEND LIST_HEADERF "${CMAKE_CURRENT_SOURCE_DIR}/meshgen/pdm_dcube_nodal_gen.F90")
  list(APPEND LIST_HEADERF "${CMAKE_CURRENT_SOURCE_DIR}/mesh/pdm_dmesh_nodal.F90")
  list(APPEND LIST_HEADERF "${CMAKE_CURRENT_BINARY_DIR}/pdm_generate_mesh.mod")

  list(APPEND LIST_HEADERF "${CMAKE_CURRENT_SOURCE_DIR}/util/pdm_writer_wrapper.F90")

endif ()

if (OPENMP_FOUND)
   set_target_properties(pdm PROPERTIES
                         COMPILE_FLAGS ${OpenMP_C_FLAGS})
endif()

# Options du link
# ---------------

unset(LINK_LIBRARIES)
list(APPEND LINK_LIBRARIES ${COMMON_LIBRARIES})

if (CMAKE_C_COMPILER_ID STREQUAL "GNU")
  list(APPEND LINK_LIBRARIES "m")
endif()

if (PDM_ENABLE_SHARED)
  set_property(TARGET ${pdm} PROPERTY POSITION_INDEPENDENT_CODE 1)
  add_library(pdm_shared SHARED  $<TARGET_OBJECTS:pdm>)
  set_target_properties(pdm_shared PROPERTIES OUTPUT_NAME "pdm")
  set_target_properties(pdm_shared PROPERTIES VERSION ${PDM_VERSION})
  if (OPENMP_FOUND)
    set_target_properties(pdm_shared PROPERTIES
                          COMPILE_FLAGS ${OpenMP_C_FLAGS}
                          LINK_FLAGS ${OpenMP_C_FLAGS})
  endif()
  target_link_libraries(pdm_shared ${LINK_LIBRARIES} pdm_io_shared pdm_mpi_shared)
  set_target_properties(pdm_shared PROPERTIES
                            PUBLIC_HEADER "${LIST_HEADER}")

  get_target_property(PDM_INC pdm INCLUDE_DIRECTORIES)
  target_include_directories(pdm_shared PUBLIC ${PDM_INC})
  if (LAPACK_FOUND)
    target_link_libraries(pdm_shared LAPACK::LAPACK)
  else ()
    if (BLAS_FOUND)
      target_link_libraries(pdm_shared LAPACK::LAPACK)
    endif()
  endif()

endif()

if (PDM_ENABLE_STATIC)
  set_property(TARGET ${pdm} PROPERTY POSITION_INDEPENDENT_CODE 1)
  add_library(pdm_static STATIC  $<TARGET_OBJECTS:pdm>)
  set_target_properties(pdm_static PROPERTIES OUTPUT_NAME "pdm")
  set_target_properties(pdm_static PROPERTIES VERSION ${PDM_VERSION})

  if (OPENMP_FOUND)
    set_target_properties(pdm_static PROPERTIES
                          COMPILE_FLAGS ${OpenMP_C_FLAGS}
                          LINK_FLAGS ${OpenMP_C_FLAGS})
  endif()
  target_link_libraries(pdm_static ${LINK_LIBRARIES} pdm_io_static pdm_mpi_static)
  set_target_properties(pdm_static PROPERTIES
                             PUBLIC_HEADER "${LIST_HEADER}")
  get_target_property(PDM_INC pdm INCLUDE_DIRECTORIES)
  target_include_directories(pdm_static PUBLIC ${PDM_INC})

endif()

if (PDM_ENABLE_SHARED_Fortran_INTERFACE)
  set_property(TARGET ${pdmf} PROPERTY POSITION_INDEPENDENT_CODE 1)
  add_library(pdmf_shared SHARED  $<TARGET_OBJECTS:pdmf>)
  set_target_properties(pdmf_shared PROPERTIES OUTPUT_NAME "pdmf")
  set_target_properties(pdmf_shared PROPERTIES VERSION ${PDM_VERSION})
  target_link_libraries(pdmf_shared pdm_shared)
  set_target_properties(pdmf_shared PROPERTIES
                            PUBLIC_HEADER "${LIST_HEADERF}")

  get_target_property(PDM_INC pdmf INCLUDE_DIRECTORIES)
  target_include_directories(pdmf_shared PUBLIC ${PDM_INC})

endif()

if (PDM_ENABLE_STATIC_Fortran_INTERFACE)
  set_property(TARGET ${pdmf} PROPERTY POSITION_INDEPENDENT_CODE 1)
  add_library(pdmf_static STATIC  $<TARGET_OBJECTS:pdmf>)
  set_target_properties(pdmf_static PROPERTIES OUTPUT_NAME "pdmf")
  set_target_properties(pdmf_static PROPERTIES VERSION ${PDM_VERSION})

  target_link_libraries(pdmf_static pdm_static)
  set_target_properties(pdmf_static PROPERTIES
                             PUBLIC_HEADER "${LIST_HEADERF}")

  get_target_property(PDM_INC pdmf INCLUDE_DIRECTORIES)
  target_include_directories(pdmf_static PUBLIC ${PDM_INC})

endif()

if (PDM_ENABLE_CUDA)
  set_property(TARGET ${pdm_gpu} PROPERTY POSITION_INDEPENDENT_CODE 1)

  if (PDM_ENABLE_STATIC)
    target_link_libraries(pdm_gpu pdm_static)
  elseif (PDM_ENABLE_SHARED)
    target_link_libraries(pdm_gpu pdm_shared)
  endif()
  set_target_properties(pdm_gpu PROPERTIES
                             PUBLIC_HEADER "${LIST_HEADER}")

  get_target_property(PDM_INC pdm_gpu INCLUDE_DIRECTORIES)
  target_include_directories(pdm_gpu PUBLIC ${PDM_INC})
endif()

# Spécification de l'installation
# -------------------------------

set (_targets "")

list (APPEND _targets pdm)
if (PDM_ENABLE_Fortran)
  list (APPEND _targets pdmf)
endif()
if (PDM_ENABLE_SHARED)
  list (APPEND _targets pdm_shared)
endif()

if (PDM_ENABLE_STATIC)
  list (APPEND _targets pdm_static)
endif()

if (PDM_ENABLE_SHARED_Fortran_INTERFACE)
  list (APPEND _targets pdmf_shared)
endif()

if (PDM_ENABLE_STATIC_Fortran_INTERFACE)
  list (APPEND _targets pdmf_static)
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

set (PDM_LIST_TEST ${PDM_LIST_TEST} PARENT_SCOPE )
set (PDM_LIST_NRANK ${PDM_LIST_NRANK} PARENT_SCOPE )

