test_c_create (pdm_t_laplacian_smoothing_seq_ex 1 PDM_LIST_TEST PDM_LIST_NRANK)
test_c_create (pdm_t_laplacian_smoothing_seq_sol 1 PDM_LIST_TEST PDM_LIST_NRANK)

test_c_create (pdm_t_laplacian_smoothing_par_ex 1  PDM_LIST_TEST PDM_LIST_NRANK)
test_c_create (pdm_t_laplacian_smoothing_part_to_part_par_ex 1 PDM_LIST_TEST PDM_LIST_NRANK)

#TODO: Make part_extension for 2D meshes 
#test_c_create (pdm_t_laplacian_smoothing_part_extension_ex 1 PDM_LIST_TEST PDM_LIST_NRANK)
test_c_create (pdm_t_laplacian_smoothing_part_extension_3d_ex 1 PDM_LIST_TEST PDM_LIST_NRANK)

test_c_create (pdm_t_cell_extents_ex 1 PDM_LIST_TEST PDM_LIST_NRANK)
test_c_create (pdm_t_cell_extents_sol 1 PDM_LIST_TEST PDM_LIST_NRANK)

test_c_create (pdm_t_tag_cells_ex 1 PDM_LIST_TEST PDM_LIST_NRANK)
test_c_create (pdm_t_tag_cells_sol 1 PDM_LIST_TEST PDM_LIST_NRANK)

test_c_create (pdm_t_groups_ex 1 PDM_LIST_TEST PDM_LIST_NRANK)
test_c_create (pdm_t_groups_sol 1 PDM_LIST_TEST PDM_LIST_NRANK)

test_c_create (pdm_t_ordered_neighbours_2d_triangles_seq 1 PDM_LIST_TEST PDM_LIST_NRANK)
test_c_create (pdm_t_ordered_neighbours_2d_seq 1 PDM_LIST_TEST PDM_LIST_NRANK)

set (PDM_LIST_TEST ${PDM_LIST_TEST} PARENT_SCOPE )
set (PDM_LIST_NRANK ${PDM_LIST_NRANK} PARENT_SCOPE )
