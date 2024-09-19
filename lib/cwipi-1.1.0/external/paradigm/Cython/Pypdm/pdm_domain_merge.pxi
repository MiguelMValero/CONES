cdef extern from "pdm_multi_block_merge.h":
  ctypedef struct PDM_multi_block_merge_t:
      pass
  PDM_multi_block_merge_t* PDM_multi_block_merge_create(PDM_g_num_t  **block_distrib_idx,
                                                 const int            n_block,
                                                       int           *n_selected,
                                                       PDM_g_num_t  **selected_g_num,
                                                       int            graph_size,
                                                       int           *dmerge_idx,
                                                       int           *dmerge_block_id,
                                                       PDM_g_num_t   *dmerge_g_num,
                                                       PDM_MPI_Comm   comm)

  int PDM_multi_block_merge_get_n_block(PDM_multi_block_merge_t* mbm)
  PDM_g_num_t* PDM_multi_block_merge_get_distrib(PDM_multi_block_merge_t* mbm)
  PDM_g_num_t* PDM_multi_block_merge_get_multi_distrib(PDM_multi_block_merge_t* mbm)

  void PDM_multi_block_merge_get_old_to_new(PDM_multi_block_merge_t  *mbm,
                                            PDM_g_num_t             **old_distrib,
                                            int                     **dold_to_new_idx,
                                            PDM_g_num_t             **dold_to_new)

  void PDM_multi_block_merge_exch(PDM_multi_block_merge_t     *mbm,
                                  size_t                       s_data,
                                  PDM_stride_t                 t_stride,
                                  int                        **block_stride,
                                  void                       **block_data,
                                  int                        **merge_block_stride,
                                  void                       **merge_block_data)

  void PDM_multi_block_merge_exch_and_update(PDM_multi_block_merge_t  *mbm_cur,
                                             PDM_multi_block_merge_t  *mbm_for_update,
                                             PDM_stride_t              t_stride,
                                             int                     **block_stride,
                                             PDM_g_num_t             **block_data,
                                             int                     **block_domain,
                                             int                     **merge_block_stride,
                                             PDM_g_num_t             **merge_block_data)

  void PDM_multi_block_merge_free(PDM_multi_block_merge_t* mbm)

# ===================================================================================

cdef class MultiBlockMerge:

  # Class attributes
  cdef PDM_multi_block_merge_t* m_block_merge
  cdef int                      n_domain
  cdef int                      i_rank
  cdef int                      n_rank

  def  __cinit__(self,
                 int n_domain,
                 list block_distri_l,
                 list selected_per_block,
                 dict interface_graph,
                 MPI.Comm comm):
    """
    Create a pdm structure multi_block_merge.
    Arguments :
      - n_domain            -> number of domains to merge
      - block_distri_l      -> distribution of each block to merge
      - selected_per_blocks -> list of ids (global to each domain) to include in the merge
      - interface_graph     -> description of domains interfaces under graph format :
                               dict must contain graph_idx, graph_ids and graph_dom keys
      - comm                -> MPI communicator 

    """

    # Generic declarations
    cdef NPY.ndarray[NPY.int32_t,    ndim=1, mode='c'] numpy_int
    cdef NPY.ndarray[npy_pdm_gnum_t, ndim=1, mode='c'] numpy_gnum

    #Some checks
    assert len(block_distri_l) == len(selected_per_block) == n_domain
    for i in range(n_domain):
      assert_single_dim_np(block_distri_l[i],     npy_pdm_gnum_dtype)
      assert_single_dim_np(selected_per_block[i], npy_pdm_gnum_dtype)
    for key in ['graph_idx', 'graph_ids', 'graph_dom']:
      assert key in interface_graph
    assert_single_dim_np(interface_graph['graph_idx'], NPY.int32)
    assert_single_dim_np(interface_graph['graph_ids'], npy_pdm_gnum_dtype)
    assert_single_dim_np(interface_graph['graph_dom'], NPY.int32)

    cdef MPI.MPI_Comm c_comm = comm.ob_mpi
    cdef PDM_MPI_Comm PDMC   = PDM_MPI_mpi_2_pdm_mpi_comm(&c_comm)

    cdef int*          _n_selected         = list_to_int_pointer([array.size for array in selected_per_block])
    cdef PDM_g_num_t** _distri_per_block   = np_list_to_gnum_pointers(block_distri_l)
    cdef PDM_g_num_t** _selected_per_block = np_list_to_gnum_pointers(selected_per_block)

    numpy_int = interface_graph['graph_idx']
    cdef int* _graph_idx = <int *> numpy_int.data
    numpy_int = interface_graph['graph_dom']
    cdef int* _graph_dom = <int *> numpy_int.data
    numpy_gnum = interface_graph['graph_ids']
    cdef PDM_g_num_t* _graph_ids = <PDM_g_num_t *> numpy_gnum.data

    self.m_block_merge = PDM_multi_block_merge_create(_distri_per_block,
                                                      n_domain,
                                                      _n_selected,
                                                      _selected_per_block,
                                                <int> interface_graph['graph_idx'].size - 1,
                                                      _graph_idx,
                                                      _graph_dom,
                                                      _graph_ids,
                                                      PDMC)
    self.n_domain = n_domain
    self.i_rank = comm.Get_rank()
    self.n_rank = comm.Get_size()

    free(_n_selected)
    free(_distri_per_block)
    free(_selected_per_block)

  def get_merged_dn(self):
    """
    Return the (distributed) number of elements in the merged block
    """
    return PDM_multi_block_merge_get_n_block(self.m_block_merge)

  def get_merged_distri(self):
    """
    Return the distribution of the merged block
    """
    cdef PDM_g_num_t* _merged_distri = PDM_multi_block_merge_get_distrib(self.m_block_merge)
    merged_distri = create_numpy_g(_merged_distri, self.n_rank+1, flag_owndata=False)
    return NPY.copy(merged_distri)

  def get_blocks_distri(self):
    """
    Return the (shifted) number of entiy in each original block
    """
    cdef PDM_g_num_t* _blocks_distri = PDM_multi_block_merge_get_multi_distrib(self.m_block_merge)
    blocks_distri = create_numpy_g(_blocks_distri, self.n_domain+1, flag_owndata=False)
    return NPY.copy(blocks_distri)

  def get_old_to_new(self):
    """
    Return old to new ordering for the initial blocks. This ordering is distributed over
    the processes following old_distrib distribution. Entity are implicitly ordered 
    following initial blocks order.
    """
    cdef PDM_g_num_t* _old_distrib;
    cdef int*         _old_to_new_idx
    cdef PDM_g_num_t* _old_to_new
    PDM_multi_block_merge_get_old_to_new(self.m_block_merge,
                                         &_old_distrib,
                                         &_old_to_new_idx,
                                         &_old_to_new)
    old_distrib    = create_numpy_g(_old_distrib, self.n_rank+1, flag_owndata=False)
    cdef int dn_olds = _old_distrib[self.i_rank+1] - _old_distrib[self.i_rank]
    old_to_new_idx = create_numpy_i(_old_to_new_idx, dn_olds+1, flag_owndata=False)
    old_to_new     = create_numpy_g(_old_to_new, _old_to_new_idx[dn_olds], flag_owndata=False)
    return NPY.copy(old_distrib), NPY.copy(old_to_new_idx), NPY.copy(old_to_new)

  def merge_field(self, list block_data, list block_stride=None):
    """
    Transfer a field supported by (all the entities) of the orignal blocks to the merged block
    """
    var_stride = bool(block_stride is not None)

    assert len(block_data) == self.n_domain
    if var_stride:
      assert len(block_stride) == self.n_domain

    cdef size_t s_data = block_data[0].dtype.itemsize

    cdef int  **_blocks_stride = NULL
    if var_stride:
      _blocks_stride = np_list_to_int_pointers(block_stride)
    else:
      _blocks_stride = <int **> malloc(self.n_domain * sizeof(int *))
      for i in range(self.n_domain):
        _blocks_stride[i] = <int *> malloc(1*sizeof(int))
        _blocks_stride[i][0] = 1
    cdef void **_blocks_data = np_list_to_void_pointers(block_data)

    cdef PDM_stride_t pdm_stride_t = PDM_STRIDE_VAR_INTERLACED if var_stride else PDM_STRIDE_CST_INTERLACED

    cdef int*  _merged_block_stride = NULL
    cdef void* _merged_block_data = NULL
    PDM_multi_block_merge_exch(self.m_block_merge,
                               s_data,
                               pdm_stride_t,
                               _blocks_stride,
                               _blocks_data,
                              &_merged_block_stride,
                              &_merged_block_data)
    cdef NPY.npy_intp dim_np
    cdef int dn_merged = PDM_multi_block_merge_get_n_block(self.m_block_merge)
    if var_stride:
      merged_block_stride = create_numpy_i(_merged_block_stride, dn_merged)
      dim_np = <NPY.npy_intp> merged_block_stride.sum()
    else:
      assert _merged_block_stride == NULL
      dim_np = <NPY.npy_intp> dn_merged

    #We use first encoutered array to infer numpy type
    merged_block_data = create_numpy(_merged_block_data, block_data[0].dtype.num, dim_np)

    if not var_stride:
      for i in range(self.n_domain):
        free(_blocks_stride[i])
    free(_blocks_stride)
    free(_blocks_data)
  
    if var_stride:
      return merged_block_stride, merged_block_data
    else:
      return merged_block_data


  def merge_and_update(self, MultiBlockMerge other, list block_data, list block_stride=None, list block_domain=None):
    """
    Transfer a field of gnums supported by (all the entities) of the orignal blocks to the merged block.
    Then update it using other MultiBlockMerge object old to new order.
    If block_domain is None, gnums are supposed to belong to the current domain. Otherwise, block_domain
    must have same layout than block_data and indicate the reference domain of each entity before merge.
    """
    var_stride = bool(block_stride is not None)

    assert len(block_data) == self.n_domain
    if var_stride:
      assert len(block_stride) == self.n_domain
    if block_domain is not None:
      assert len(block_domain) == self.n_domain
      for i,bd in enumerate(block_domain): 
        assert_single_dim_np(bd, NPY.int32, size=block_data[i].size)

    cdef size_t s_data = block_data[0].dtype.itemsize

    #Get data

    cdef int**         _blocks_stride = NULL
    if var_stride:
      _blocks_stride = np_list_to_int_pointers(block_stride)
    else:
      _blocks_stride = <int **> malloc(self.n_domain * sizeof(int *))
      for i in range(self.n_domain):
        _blocks_stride[i] = <int *> malloc(1*sizeof(int))
        _blocks_stride[i][0] = 1
    cdef PDM_g_num_t** _blocks_data = np_list_to_gnum_pointers(block_data)
    cdef int**         _blocks_domain = NULL
    if block_domain is not None:
      _blocks_domain = np_list_to_int_pointers(block_domain)



    cdef PDM_stride_t pdm_stride_t = PDM_STRIDE_VAR_INTERLACED if var_stride else PDM_STRIDE_CST_INTERLACED

    cdef int*  _merged_block_stride = NULL
    cdef PDM_g_num_t* _merged_block_data = NULL

    # Run function
    PDM_multi_block_merge_exch_and_update(self.m_block_merge,
                                          other.m_block_merge,
                                          pdm_stride_t,
                                          _blocks_stride,
                                          _blocks_data,
                                          _blocks_domain,
                                         &_merged_block_stride,
                                         &_merged_block_data)

    cdef NPY.npy_intp dim_np
    cdef int dn_merged = PDM_multi_block_merge_get_n_block(self.m_block_merge)
    if var_stride:
      merged_block_stride = create_numpy_i(_merged_block_stride, dn_merged)
      dim_np = <NPY.npy_intp> merged_block_stride.sum()
    else:
      #assert _merged_block_stride == NULL #TODO use it to compute new size (manage suppression)
      dim_np = <NPY.npy_intp> dn_merged

    #We use first encoutered array to infer numpy type
    merged_block_data = create_numpy(_merged_block_data, block_data[0].dtype.num, dim_np)

    if not var_stride:
      for i in range(self.n_domain):
        free(_blocks_stride[i])
    if block_domain is not None:
      free(_blocks_domain)
    free(_blocks_stride)
    free(_blocks_data)

    if var_stride:
      return merged_block_stride, merged_block_data
    else:
      return merged_block_data

  def __dealloc__(self):
    """
    Free pdm structure multi_block_merge
    """
    PDM_multi_block_merge_free(self.m_block_merge)

