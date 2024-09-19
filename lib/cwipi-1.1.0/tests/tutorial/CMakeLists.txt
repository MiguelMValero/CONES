
# C

include(TestCreate)

if (CWP_ENABLE_Fortran)

  list (APPEND exe_tuto_ex c_new_api_deformable_ex.c fortran_new_api_deformable_ex.F90)
  list (APPEND exe_tuto_sol c_new_api_deformable_sol.c fortran_new_api_deformable_sol.F90)
  list (APPEND nproc_tuto 2 2) 

  test_create ("${exe_tuto_ex}" "${nproc_tuto}" CWP_LIST_TEST CWP_LIST_NRANK)
  test_create ("${exe_tuto_sol}" "${nproc_tuto}" CWP_LIST_TEST CWP_LIST_NRANK)

endif()

test_c_create (c_new_api_polygon_sol     2 CWP_LIST_TEST CWP_LIST_NRANK)
test_c_create (c_new_api_polygon_ex      2 CWP_LIST_TEST CWP_LIST_NRANK)

# fortran
if (CWP_ENABLE_Fortran)
  test_fortran_create (fortran_new_api_polygon_sol    2 CWP_LIST_TEST CWP_LIST_NRANK)
  test_fortran_create (fortran_new_api_polygon_ex     2 CWP_LIST_TEST CWP_LIST_NRANK)
endif()

# python
if (CWP_ENABLE_PYTHON_BINDINGS)
  test_python_create (python_new_api_polygon_sol 2 CWP_LIST_TEST CWP_LIST_NRANK)
  test_python_create (python_new_api_polygon_ex  2 CWP_LIST_TEST CWP_LIST_NRANK)

  if (CWP_ENABLE_Fortran)

    list (APPEND call_back_tuto_sol c_new_api_callback_sendrecv_sol.c fortran_new_api_callback_sendrecv_sol.F90 python_new_api_callback_sendrecv_sol.py)
    list (APPEND call_back_nproc_tuto 2 2 2)

    test_create ("${call_back_tuto_sol}" "${call_back_nproc_tuto}" CWP_LIST_TEST CWP_LIST_NRANK)
 
  endif()
endif()

set (CWP_LIST_TEST  ${CWP_LIST_TEST} PARENT_SCOPE )
set (CWP_LIST_NRANK ${CWP_LIST_NRANK} PARENT_SCOPE )
