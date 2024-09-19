class_name="pdm_t_part_to_block_geom_s=10000"
test_name="pdm_t_part_to_block_geom_n=1_s=10000"
test_n_proc="1"
execute_test 'mpirun -np 1 ./test/pdm_t_part_to_block_geom  -s 10000'
class_name="pdm_t_part_to_block_geom_s=100000"
test_name="pdm_t_part_to_block_geom_n=1_s=100000"
test_n_proc="1"
execute_test 'mpirun -np 1 ./test/pdm_t_part_to_block_geom  -s 100000'
class_name="pdm_t_part_to_block_geom_s=10000"
test_name="pdm_t_part_to_block_geom_n=2_s=10000"
test_n_proc="2"
execute_test 'mpirun -np 2 ./test/pdm_t_part_to_block_geom  -s 10000'
class_name="pdm_t_part_to_block_geom_s=100000"
test_name="pdm_t_part_to_block_geom_n=2_s=100000"
test_n_proc="2"
execute_test 'mpirun -np 2 ./test/pdm_t_part_to_block_geom  -s 100000'
class_name="pdm_t_part_to_block_geom_s=10000"
test_name="pdm_t_part_to_block_geom_n=3_s=10000"
test_n_proc="3"
execute_test 'mpirun -np 3 ./test/pdm_t_part_to_block_geom  -s 10000'
class_name="pdm_t_part_to_block_geom_s=100000"
test_name="pdm_t_part_to_block_geom_n=3_s=100000"
test_n_proc="3"
execute_test 'mpirun -np 3 ./test/pdm_t_part_to_block_geom  -s 100000'
class_name="pdm_t_part_to_block_geom_s=10000"
test_name="pdm_t_part_to_block_geom_n=4_s=10000"
test_n_proc="4"
execute_test 'mpirun -np 4 ./test/pdm_t_part_to_block_geom  -s 10000'
class_name="pdm_t_part_to_block_geom_s=100000"
test_name="pdm_t_part_to_block_geom_n=4_s=100000"
test_n_proc="4"
execute_test 'mpirun -np 4 ./test/pdm_t_part_to_block_geom  -s 100000'
