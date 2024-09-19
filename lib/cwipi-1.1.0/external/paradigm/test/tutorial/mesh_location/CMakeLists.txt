test_c_create (pdm_t_mesh_location_sol 2 PDM_LIST_TEST PDM_LIST_NRANK)

if (PDM_ENABLE_STATIC_Fortran_INTERFACE OR PDM_ENABLE_SHARED_Fortran_INTERFACE)
  test_fortran_create (pdm_t_mesh_location_sol_f 2 PDM_LIST_TEST PDM_LIST_NRANK)
endif()

if (PDM_ENABLE_PYTHON_BINDINGS)
  test_python_create (pdm_t_mesh_location_sol_p 2 PDM_LIST_TEST PDM_LIST_NRANK)
  add_dependencies (pdm_t_mesh_location_sol_p Pypdm)
endif()

set (PDM_LIST_TEST ${PDM_LIST_TEST} PARENT_SCOPE )
set (PDM_LIST_NRANK ${PDM_LIST_NRANK} PARENT_SCOPE )
