cdef extern from "pdm_distrib.h":
  void PDM_distrib_weight(int            sampling_factor,
                          int            n_active_ranks,
                          int            n_part,
                          int           *n_elmts,
                          PDM_g_num_t  **ln_to_gn,
                          double       **weight,
                          int            n_iter_max,
                          double         tolerance,
                          PDM_MPI_Comm   comm,
                          PDM_g_num_t  **rank_index)


# ===================================================================================

def compute_weighted_distribution(list ln_to_gn_list, list weights, MPI.Comm comm):

  cdef int n_part = len(ln_to_gn_list)
  assert len(ln_to_gn_list) == len(weights)
  for i in range(n_part):
    assert_single_dim_np(ln_to_gn_list[i], npy_pdm_gnum_dtype)
    assert_single_dim_np(weights[i], NPY.double)
    assert ln_to_gn_list[i].size == weights[i].size

  cdef MPI.MPI_Comm c_comm = comm.ob_mpi
  cdef PDM_MPI_Comm PDMC   = PDM_MPI_mpi_2_pdm_mpi_comm(&c_comm)

  cdef int n_active_ranks = comm.Get_size()
  # Those parameters are hard codded in pdm_part_to_block so re-use it
  cdef int sampling_factor = 2
  cdef int n_iter_max = 5
  cdef double tolerance = 0.5

  cdef int* n_elmts = list_to_int_pointer([ln_to_gn.size for ln_to_gn in ln_to_gn_list])
  cdef PDM_g_num_t** _ln_to_gn = np_list_to_gnum_pointers(ln_to_gn_list)
  cdef double**       _weights = np_list_to_double_pointers(weights)

  cdef PDM_g_num_t* distrib = NULL
  PDM_distrib_weight(sampling_factor,
                     n_active_ranks,
                     n_part,
                     n_elmts,
                     _ln_to_gn,
                     _weights,
                     n_iter_max,
                     tolerance,
                     PDMC,
                    &distrib)
  
  free(n_elmts)
  free(_ln_to_gn)
  free(_weights)

  return create_numpy_g(distrib, n_active_ranks+1)
