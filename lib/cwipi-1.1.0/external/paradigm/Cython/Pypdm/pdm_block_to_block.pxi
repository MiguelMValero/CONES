import warnings

cdef extern from "pdm_block_to_block.h":
    # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    # > Wrapping of Ppart Structure
    ctypedef struct PDM_block_to_block_t:
      pass
    # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    # > Wrapping of function
    PDM_block_to_block_t *PDM_block_to_block_create(PDM_g_num_t   *blockDistribIniIdx,
                                                    PDM_g_num_t   *blockDistribEndIdx,
                                                    PDM_MPI_Comm   comm)

    int PDM_block_to_block_exch(PDM_block_to_block_t  *btb,
                                size_t                s_data,
                                PDM_stride_t          t_stride,
                                int                   cst_stride,
                                int                  *block_stride_ini,
                                void                 *block_data_ini,
                                int                  *block_stride_end,
                                void                **block_data_end)

    PDM_block_to_block_t *PDM_block_to_block_free(PDM_block_to_block_t *btb)
    # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::


# ------------------------------------------------------------------
cdef class BlockToBlock:
    """
       BlockToPart: Interface for block_to_block.c
    """
    # ************************************************************************
    # > Class attributes
    cdef PDM_block_to_block_t* BTB
    cdef int                   dn_elt_ini
    cdef MPI.Comm              py_comm
    # ************************************************************************
    # ------------------------------------------------------------------------
    def __cinit__(self, NPY.ndarray[npy_pdm_gnum_t, ndim=1, mode='c'] DistribIni,
                        NPY.ndarray[npy_pdm_gnum_t, ndim=1, mode='c'] DistribEnd,
                        MPI.Comm comm):
        """
        Constructor of BlockToBlock object : Python wrapping of PDM library (E. Quémerais)

            :param comm:        MPI Communicator (Caution MPI Comm is a mpi4py object )
            :param DistribIni:  Distribution of distribute array in initiale frame (Size = nRank+1)
            :param DistribEnd:  Distribution of distribute array in finale frame (Size = nRank+1)

        """
        # > Some checks
        assert_single_dim_np(DistribIni, npy_pdm_gnum_dtype, comm.Get_size()+1)
        assert_single_dim_np(DistribEnd, npy_pdm_gnum_dtype, comm.Get_size()+1)

        # > Store class parameters
        self.dn_elt_ini  = DistribIni[comm.Get_rank() + 1] - DistribIni[comm.Get_rank()]
        self.py_comm = comm

        # > Convert input data
        cdef MPI.MPI_Comm c_comm = comm.ob_mpi
        cdef PDM_MPI_Comm PDMC   = PDM_MPI_mpi_2_pdm_mpi_comm(&c_comm)

        # > Create PDM structure
        self.BTB = PDM_block_to_block_create(<PDM_g_num_t *> DistribIni.data,
                                             <PDM_g_num_t *> DistribEnd.data,
                                             PDMC)

    # ------------------------------------------------------------------
    def exchange_field(self, NPY.ndarray block_data, block_stride=1, bint interlaced_str=True):
      """
      Wrapping for PDM_block_to_block_exch: transfert a distributed data from the initial distribution
      to the destination distribution

      :param self:         BlockToPart object
      :param block_data:   Distributed data array, 1 dimensional and with same datatype for each rank
      :param block_stride: Stride for distributed array. Can be either an array of size dn_elt (variable
                           stride will be used) or an integer (cst stride will be used)
      :param interlaced_str: indicate if data are interlaced (True) or interleaved 
      """

      cdef NPY.ndarray[NPY.int32_t, ndim=1, mode='c'] numpy_int
      cdef int   _block_stride_cst = 0
      cdef int*  _block_stride = NULL
      if isinstance(block_stride, int):
        _stride_t = PDM_STRIDE_CST_INTERLACED if interlaced_str else PDM_STRIDE_CST_INTERLEAVED
        _block_stride_cst = block_stride
        assert_single_dim_np(block_data, block_data.dtype, block_stride*self.dn_elt_ini)
      elif isinstance(block_stride, NPY.ndarray):
        raise NotImplementedError("Variable stride for BTB not implemented")
#        _stride_t = PDM_STRIDE_VAR_INTERLACED
#        assert_single_dim_np(block_stride, NPY.int32, self.dn_elt_ini)
#        assert_single_dim_np(block_data, block_data.dtype, block_stride.sum())
#        numpy_int = block_stride
#        _block_stride = <int *> numpy_int.data
      else:
        raise ValueError("Invalid stride in BtP exchange")

      cdef size_t s_data = block_data.dtype.itemsize #Should be the same for all procs
      assert self.py_comm.allreduce(s_data, op=MPI.MIN) == self.py_comm.allreduce(s_data, op=MPI.MAX) == s_data

      cdef int*  _block_stride_end = NULL
      cdef void* _block_data_end   = NULL
      c_size = PDM_block_to_block_exch(self.BTB,
                                       s_data,
                                       _stride_t,
                                       _block_stride_cst,
                                       _block_stride,
                              <void *> block_data.data,
                                       _block_stride_end,
                                      &_block_data_end)

      return create_numpy(_block_data_end, block_data.dtype.num, c_size)

    # ------------------------------------------------------------------
    def BlockToBlock_Exchange(self, dict dFieldIni, 
                                    dict dFieldEnd,
                                    PDM_stride_t t_stride = <PDM_stride_t> -1,
                                    blkStrid=1,
                                    bint interlaced_str=True):
      """ Shortcut to exchange multiple fieds stored in dict """
      if t_stride != -1:
        warnings.warn("Parameter t_stride is deprecated and will be removed in further release",
          DeprecationWarning, stacklevel=2)
      for field_name, field in dFieldIni.items():
        block_data_end = self.exchange_field(field, blkStrid, interlaced_str)
        dFieldEnd[field_name] = block_data_end

    # ------------------------------------------------------------------
    def __dealloc__(self):
      """ Free structure """
      cdef PDM_block_to_block_t* none = PDM_block_to_block_free(self.BTB)
      assert none == NULL

