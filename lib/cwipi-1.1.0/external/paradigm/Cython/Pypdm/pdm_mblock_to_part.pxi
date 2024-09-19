cdef extern from "pdm_multi_block_to_part.h":
    # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    # > Wrapping of Ppart Structure
    ctypedef struct PDM_multi_block_to_part_t:
      pass
    # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    # > Wrapping of function

    PDM_multi_block_to_part_t* PDM_multi_block_to_part_create(const PDM_g_num_t   *multi_distrib_idx,
                                                              const int            n_block,
                                                              const PDM_g_num_t  **block_distrib_idx,
                                                              const PDM_g_num_t  **gnum_elt,
                                                              const int           *n_elt,
                                                              const int            n_part,
                                                              const PDM_MPI_Comm   comm)
    
    void PDM_multi_block_to_part_exch2(PDM_multi_block_to_part_t   *mbtp,
                                       size_t                       s_data,
                                       PDM_stride_t                 t_stride,
                                       int                        **block_stride,
                                       void                       **block_data,
                                       int                       ***part_stride,
                                       void                      ***part_data)



    PDM_multi_block_to_part_t* PDM_multi_block_to_part_free(PDM_multi_block_to_part_t *mbtp)


# ------------------------------------------------------------------
cdef class MultiBlockToPart:
    """
    """
    # ************************************************************************
    # > Class attributes
    cdef PDM_multi_block_to_part_t* MBTP
    cdef int n_block
    cdef int n_part
    cdef int*  dn_elt
    cdef int*  pn_elt
    # ************************************************************************
    # ------------------------------------------------------------------------
    def __cinit__(self, list distrib_list, list ln_to_gn_list, MPI.Comm comm):
        
        rank = comm.Get_rank()
        self.n_block = len(distrib_list)
        self.n_part  = len(ln_to_gn_list)
        self.dn_elt = list_to_int_pointer([distrib[rank+1] - distrib[rank] for distrib in distrib_list])
        
        # > Some checks
        for ln_to_gn in ln_to_gn_list:
          assert_single_dim_np(ln_to_gn, npy_pdm_gnum_dtype)
        for distrib in distrib_list:
          assert_single_dim_np(distrib, npy_pdm_gnum_dtype, comm.Get_size()+1)

        # > Convert input data
        cdef MPI.MPI_Comm c_comm = comm.ob_mpi
        cdef PDM_MPI_Comm PDMC   = PDM_MPI_mpi_2_pdm_mpi_comm(&c_comm)
        
        self.pn_elt = list_to_int_pointer([lngn.size for lngn in ln_to_gn_list])

        cdef PDM_g_num_t** _distrib  = np_list_to_gnum_pointers(distrib_list)
        cdef PDM_g_num_t** _ln_to_gn = np_list_to_gnum_pointers(ln_to_gn_list)

        cdef PDM_g_num_t* _multi_distrib_idx = <PDM_g_num_t*> malloc((self.n_block+1) * sizeof(PDM_g_num_t))
        _multi_distrib_idx[0] = 0
        for i in range(self.n_block):
          _multi_distrib_idx[i+1] = _multi_distrib_idx[i] + distrib_list[i][comm.Get_size()]

        # > Create PDM structure
        self.MBTP = PDM_multi_block_to_part_create(_multi_distrib_idx,
                                                   self.n_block,
                                                   _distrib,
                                                   _ln_to_gn,
                                                   self.pn_elt,
                                                   self.n_part,
                                                   PDMC)
        # Free working arrays
        free(_ln_to_gn)
        free(_distrib)
        free(_multi_distrib_idx)
    
    # ------------------------------------------------------------------
    def exchange_field(self, list blocks_data, block_stride=1, bint interlaced_str=True):
      
      cdef NPY.ndarray[NPY.int32_t, ndim=1, mode='c'] numpy_int
      cdef PDM_stride_t _stride_t
      cdef int**         _blocks_stride

      assert len(blocks_data) == self.n_block
      dtype = blocks_data[0].dtype

      if isinstance(block_stride, int):
        _stride_t = PDM_STRIDE_CST_INTERLACED if interlaced_str else PDM_STRIDE_CST_INTERLEAVED
        _blocks_stride = <int **> malloc(self.n_block * sizeof(int *))
        for i in range(self.n_block):
          _blocks_stride[i]    = <int*> malloc(sizeof(int))
          _blocks_stride[i][0] = block_stride
        for i,block_data in enumerate(blocks_data):
          assert_single_dim_np(block_data, dtype, block_stride*self.dn_elt[i])
      else:
        raise ValueError("Invalid stride in MBtP exchange")

      cdef size_t s_data = dtype.itemsize

      cdef int**  _part_stride = NULL
      cdef void** _part_data   = NULL
    
      cdef void** _blocks_data = np_list_to_void_pointers(blocks_data)
    
      PDM_multi_block_to_part_exch2(self.MBTP,
                                    s_data,
                                    _stride_t,
                                    _blocks_stride,
                                    _blocks_data,
                                    &_part_stride,
                                    &_part_data)

                              
      for i in range(self.n_block):
        free(_blocks_stride[i])
      free(_blocks_stride)
      free(_blocks_data)


      part_stride = []
      part_data   = []
      cdef NPY.npy_intp dim_np

      for i_part in range(self.n_part):

        if _stride_t == PDM_STRIDE_VAR_INTERLACED:
          stride = create_numpy_i(_part_stride[i_part], self.pn_elt[i_part])
          part_stride.append(stride)
          dim_np = <NPY.npy_intp> stride.sum()
        else:
          dim_np = <NPY.npy_intp> block_stride * self.pn_elt[i_part]

        part_data.append(create_numpy(_part_data[i_part], dtype.num, dim_np))

      if _stride_t == PDM_STRIDE_VAR_INTERLACED:
        free(_part_stride)
      else:
        part_stride = None
      free(_part_data)

      return part_stride, part_data


    # ------------------------------------------------------------------
    def __dealloc__(self):
      """
      Deallocate all the array
      """
      # > Free PDM Structure
      PDM_multi_block_to_part_free(self.MBTP)
      free(self.dn_elt)
      free(self.pn_elt)

