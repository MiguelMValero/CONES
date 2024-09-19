#ifndef __PDM_DOCTEST_H__
#define __PDM_DOCTEST_H__


#define CHECK_EQ_C_ARRAY(array, array_expected, size) \
  for(int i = 0; i < size; ++i){ CHECK(array[i] == array_expected[i]); }

#define MPI_CHECK_EQ_C_ARRAY(rank_to_test, array, array_expected, size)        \
  if(rank_to_test == test_rank){                                           \
    for(int i = 0; i < size; ++i){ CHECK(array[i] == array_expected[i]); } \
  }

#endif  /* __PDM_H__ */
