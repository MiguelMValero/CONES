#
# Exemple of use :
#  - ./script/prepare_ci/run_top_test_base.sh --log_file_name tot
#  - ./script/prepare_ci/run_top_test_base.sh --log_file_name tot -k pdm_t_part_to_block_geom
#

while [ $# -gt 0 ]; do
   if [[ $1 == *"--"* ]]; then
        param="${1/--/}"
        declare $param="$2"
        echo $1 $2 # // Optional to see the parameter:value result
   fi
   shift
done


execute_test () {
  stdout_log_name=out_${test_name}.log
  stderr_log_name=err_${test_name}.log
  test_cmd=$*
  echo "Run : ", $test_name, "with command ", $*
  start=`date +%s.%N`;eval $* 2> $stderr_log_name 1> $stdout_log_name ; status=$?; end=`date +%s.%N`;
  runtime=$( echo "$end - $start" | bc -l )

  ((n_cases = n_cases + 1))
  ((n_cases_tot = n_cases_tot + 1))
  if [ "${status}" -ne "0" ]; then
    echo "${test_cmd} : Failed" >> ${log_file_name}
    echo "Failed"
    ((n_fail = n_fail + 1))
    ((n_fail_tot = n_fail_tot + 1))

    printf "        <testcase classname=\"%s%s\" name=\"%s\" time=\"%s\" status=\"fail\" > \n" $item $matrix_key $test_name $runtime >> $output_test_suite

    if [ "$(find -name '*core*' | wc -l)" ]; then
      rm *core*
    fi

    echo '<failure message="' >> $output_test_suite
    cat $stderr_log_name >> $output_test_suite
    echo '">' >> $output_test_suite
    echo '</failure>' >> $output_test_suite

  else
    echo "${test_cmd} : OK" >> ${log_file_name}
    printf "        <testcase classname=\"%s%s\" name=\"%s\" time=\"%s\" status=\"run\" > \n" $item $matrix_key $test_name $runtime >> $output_test_suite
    echo "OK"
  fi

  # echo '<system-out>' >> $output_test_suite
  # cat $stdout_log_name >> $output_test_suite
  # echo '</system-out>' >> $output_test_suite
  echo '        </testcase>' >> $output_test_suite
}

# A faire : <system-out>
# A faire : failure message=""


if [ -n "$log_file_name" ]; then
  echo "Log in file " $log_file_name
else
  log_file_name="output_test_log"
fi
rm $log_file_name

if [ -n "$matrix_key" ]; then
  echo "matrix_key " ${matrix_key}
else
  matrix_key="test_pdm"
fi

if [ -n "$output_xml" ]; then
  echo "output_xml " ${output_xml}
else
  output_xml="reports/paradigm_alltest.xml"
fi

# echo $k

n_cases_tot=0
n_fail_tot=0
if [ -n "$k" ]; then
  echo "Execute test : " $k
  output_test_suite=out_test_suite_$k
  rm $output_test_suite
  n_cases=0
  n_fail=0
  source script/prepare_ci/run_tests_$k.sh
  printf "    <testsuite errors=\"%i\" failures=\"%i\" id=\"0\" name=\"%s\" tests=\"%i\">\n" $n_fail $n_fail $k $n_cases >> reports/tmp_paradigm_alltest.xml
  cat $output_test_suite >> reports/tmp_paradigm_alltest.xml
  echo '    </testsuite>' >> reports/tmp_paradigm_alltest.xml
else
  # list=("pdm_t_dcube_nodal_gen" "pdm_t_part_to_block_geom" "pdm_t_closest_points" "pdm_t_dist" "pdm_t_dist_strip" "pdm_t_gnum" "pdm_t_global_mean" "pdm_t_gnum_location" "pdm_t_mesh_location_dcube" "pdm_t_dmesh_nodal_to_dmesh" "pdm_t_knn_cube")
  list=("pdm_t_dcube_nodal_gen" "pdm_t_part_to_block_geom" "pdm_t_closest_points" "pdm_t_dist" "pdm_t_dist_strip" "pdm_t_gnum" "pdm_t_global_mean" "pdm_t_gnum_location" "pdm_t_dmesh_nodal_to_dmesh" "pdm_t_knn_cube" "pdm_t_part_extension" "pdm_t_part_dcube")
  n_suite=0
  for item in ${list[*]}; do
    echo $item;
    output_test_suite=out_test_suite_$item
    rm $output_test_suite
    n_cases=0
    n_fail=0
    source script/prepare_ci/run_tests_$item.sh
    printf "    <testsuite errors=\"%i\" failures=\"%i\" id=\"%i\" name=\"%s\" tests=\"%i\">\n" $n_fail $n_fail $n_suite $item $n_cases >> reports/tmp_paradigm_alltest.xml
    cat $output_test_suite >> reports/tmp_paradigm_alltest.xml
    echo '    </testsuite>' >> reports/tmp_paradigm_alltest.xml
    ((n_suite = n_suite + 1))
  done
fi

rm $output_xml
echo '<?xml version="1.0" encoding="UTF-8"?>' > $output_xml
# echo '<testsuites>' >> reports/tmp_paradigm_alltest.xml
printf "    <testsuites id=\"0\" name=\"%s\"  tests=\"%i\" failures=\"%i\" >\n" $matrix_key $n_cases_tot $n_fail_tot  >> $output_xml
cat reports/tmp_paradigm_alltest.xml >> $output_xml
echo '</testsuites>' >> $output_xml


if [[ $n_fail_tot > 0 ]]; then
  echo "One or multiples test fail \"$n_fail_tot\""
  exit 1
fi

# echo '<?xml version="1.0" encoding="UTF-8"?>' > reports/tmp_paradigm_alltest.xml
# echo '<testsuites>' >> reports/tmp_paradigm_alltest.xml
# echo '    <testsuite errors="0" failures="0" id="0" name="my test suite" tests="1">' >> reports/tmp_paradigm_alltest.xml
# echo '        <testcase classname="some.class.name" name="Test1" time="123.345000"/>' >> reports/tmp_paradigm_alltest.xml
# echo '    </testsuite>' >> reports/tmp_paradigm_alltest.xml
# echo '</testsuites>' >> reports/tmp_paradigm_alltest.xml
