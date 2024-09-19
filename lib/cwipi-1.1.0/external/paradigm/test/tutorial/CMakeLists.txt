
add_subdirectory(block_to_part/ex1)
add_subdirectory(part_to_block/ex1)
add_subdirectory(gnum/ex1)
add_subdirectory(part_to_block_to_part)
add_subdirectory(mesh)
add_subdirectory(partitioning)
add_subdirectory(mesh_location)

set (PDM_LIST_TEST ${PDM_LIST_TEST} PARENT_SCOPE )
set (PDM_LIST_NRANK ${PDM_LIST_NRANK} PARENT_SCOPE )
