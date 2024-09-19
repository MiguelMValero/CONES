class_name="pdm_t_gnum_n=30"
test_name="pdm_t_gnum_n=1_n=30"
test_n_proc="1"
execute_test 'mpirun -np 1 ./test/pdm_t_gnum  -n 30'
class_name="pdm_t_gnum_n=60"
test_name="pdm_t_gnum_n=1_n=60"
test_n_proc="1"
execute_test 'mpirun -np 1 ./test/pdm_t_gnum  -n 60'
class_name="pdm_t_gnum_n=30"
test_name="pdm_t_gnum_n=2_n=30"
test_n_proc="2"
execute_test 'mpirun -np 2 ./test/pdm_t_gnum  -n 30'
class_name="pdm_t_gnum_n=60"
test_name="pdm_t_gnum_n=2_n=60"
test_n_proc="2"
execute_test 'mpirun -np 2 ./test/pdm_t_gnum  -n 60'
class_name="pdm_t_gnum_n=30"
test_name="pdm_t_gnum_n=3_n=30"
test_n_proc="3"
execute_test 'mpirun -np 3 ./test/pdm_t_gnum  -n 30'
class_name="pdm_t_gnum_n=60"
test_name="pdm_t_gnum_n=3_n=60"
test_n_proc="3"
execute_test 'mpirun -np 3 ./test/pdm_t_gnum  -n 60'
class_name="pdm_t_gnum_n=30"
test_name="pdm_t_gnum_n=4_n=30"
test_n_proc="4"
execute_test 'mpirun -np 4 ./test/pdm_t_gnum  -n 30'
class_name="pdm_t_gnum_n=60"
test_name="pdm_t_gnum_n=4_n=60"
test_n_proc="4"
execute_test 'mpirun -np 4 ./test/pdm_t_gnum  -n 60'
