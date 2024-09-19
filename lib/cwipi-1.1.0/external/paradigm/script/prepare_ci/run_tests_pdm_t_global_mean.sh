class_name="pdm_t_global_mean_n=30_parmetis"
test_name="pdm_t_global_mean_n=1_n=30_parmetis"
test_n_proc="1"
execute_test 'mpirun -np 1 ./test/pdm_t_global_mean  -n 30 -parmetis'
class_name="pdm_t_global_mean_n=60_parmetis"
test_name="pdm_t_global_mean_n=1_n=60_parmetis"
test_n_proc="1"
execute_test 'mpirun -np 1 ./test/pdm_t_global_mean  -n 60 -parmetis'
class_name="pdm_t_global_mean_n=30_parmetis"
test_name="pdm_t_global_mean_n=2_n=30_parmetis"
test_n_proc="2"
execute_test 'mpirun -np 2 ./test/pdm_t_global_mean  -n 30 -parmetis'
class_name="pdm_t_global_mean_n=60_parmetis"
test_name="pdm_t_global_mean_n=2_n=60_parmetis"
test_n_proc="2"
execute_test 'mpirun -np 2 ./test/pdm_t_global_mean  -n 60 -parmetis'
class_name="pdm_t_global_mean_n=30_parmetis"
test_name="pdm_t_global_mean_n=3_n=30_parmetis"
test_n_proc="3"
execute_test 'mpirun -np 3 ./test/pdm_t_global_mean  -n 30 -parmetis'
class_name="pdm_t_global_mean_n=60_parmetis"
test_name="pdm_t_global_mean_n=3_n=60_parmetis"
test_n_proc="3"
execute_test 'mpirun -np 3 ./test/pdm_t_global_mean  -n 60 -parmetis'
class_name="pdm_t_global_mean_n=30_parmetis"
test_name="pdm_t_global_mean_n=4_n=30_parmetis"
test_n_proc="4"
execute_test 'mpirun -np 4 ./test/pdm_t_global_mean  -n 30 -parmetis'
class_name="pdm_t_global_mean_n=60_parmetis"
test_name="pdm_t_global_mean_n=4_n=60_parmetis"
test_n_proc="4"
execute_test 'mpirun -np 4 ./test/pdm_t_global_mean  -n 60 -parmetis'
