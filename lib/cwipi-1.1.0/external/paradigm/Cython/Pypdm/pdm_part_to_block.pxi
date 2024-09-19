cdef extern from "pdm_part_to_block.h":
    # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    # > Wrapping of Ppart Structure
    ctypedef struct PDM_part_to_block_t:
        pass

    ctypedef enum PDM_part_to_block_distrib_t:
        pass

    ctypedef enum PDM_part_to_block_post_t:
        pass

    # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    # > Wrapping of function
    PDM_part_to_block_t *PDM_part_to_block_create(PDM_part_to_block_distrib_t   t_distrib,
                                                  PDM_part_to_block_post_t      t_post,
                                                  double                         partActiveNode,
                                                  PDM_g_num_t                 **gnum_elt,
                                                  double                       **weight,
                                                  int                          *n_elt,
                                                  int                           n_part,
                                                  PDM_MPI_Comm                  comm)

    PDM_part_to_block_t *PDM_part_to_block_create_from_distrib(PDM_part_to_block_distrib_t   t_distrib,
                                                   PDM_part_to_block_post_t      t_post,
                                                   double                        partActiveNode,
                                                   PDM_g_num_t                 **gnum_elt,
                                                   PDM_g_num_t                  *dataDistribIndex,
                                                   int                          *n_elt,
                                                   int                           n_part,
                                                   PDM_MPI_Comm                  comm)

    int PDM_part_to_block_n_active_ranks_get(PDM_part_to_block_t *ptb)

    int *PDM_part_to_block_active_ranks_get(PDM_part_to_block_t *ptb)

    int PDM_part_to_block_is_active_rank(PDM_part_to_block_t *ptb)

    int PDM_part_to_block_n_elt_block_get(PDM_part_to_block_t *ptb)

    PDM_g_num_t *PDM_part_to_block_block_gnum_get(PDM_part_to_block_t *ptb)
    int         *PDM_part_to_block_block_gnum_count_get(PDM_part_to_block_t *ptb)

    int PDM_part_to_block_exch(PDM_part_to_block_t       *ptb,
                               size_t                     s_data,
                               PDM_stride_t               t_stride,
                               int                        cst_stride,
                               int                      **part_stride,
                               void                     **part_data,
                               int                      **block_stride,
                               void                     **block_data)

    PDM_part_to_block_t *PDM_part_to_block_free(PDM_part_to_block_t *ptb)

    PDM_g_num_t *PDM_part_to_block_distrib_index_get(PDM_part_to_block_t *ptb)

    PDM_l_num_t *PDM_part_to_block_destination_get(PDM_part_to_block_t *ptb)
    # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::


# ------------------------------------------------------------------
cdef class PartToBlock:
    """
       PartToBlock: Interface for block_to_part.c
    """
    # ************************************************************************
    # > Class attributes
    cdef PDM_part_to_block_t* PTB
    cdef int                  n_part
    cdef int*                 pn_elt
    cdef MPI.Comm             py_comm
    # ************************************************************************
    # ------------------------------------------------------------------------
    def __cinit__(self, MPI.Comm comm, list pLNToGN, list pWeight, int partN,
                        PDM_part_to_block_distrib_t t_distrib = <PDM_part_to_block_distrib_t> (0),
                        PDM_part_to_block_post_t    t_post    = <PDM_part_to_block_post_t   > (0),
                        double partActiveNode = 1.,
                        NPY.ndarray[npy_pdm_gnum_t, mode='c', ndim=1] userDistribution=None):
        """
        Constructor of PartToBlock object : Python wrapping of PDM library (E. Quémerais)
        """
        # > Some checks
        assert(len(pLNToGN) == partN)
        for i in range(partN):
          assert_single_dim_np(pLNToGN[i], npy_pdm_gnum_dtype)
        if (pWeight is not None):
          assert(len(pWeight) == partN)
          for i in range(partN):
            assert_single_dim_np(pWeight[i], NPY.double, pLNToGN[i].size)
        if (userDistribution is not None):
          assert_single_dim_np(userDistribution, npy_pdm_gnum_dtype, comm.Get_size()+1)

        # > Store class parameters
        self.n_part  = partN
        self.pn_elt  = list_to_int_pointer([array.size for array in pLNToGN])
        self.py_comm = comm

        # > Convert input data
        cdef MPI.MPI_Comm c_comm = comm.ob_mpi
        cdef PDM_MPI_Comm PDMC   = PDM_MPI_mpi_2_pdm_mpi_comm(&c_comm)

        cdef PDM_g_num_t** _ln_to_gn = np_list_to_gnum_pointers(pLNToGN)
        cdef double** _weight = NULL
        if (pWeight is not None):
          _weight = np_list_to_double_pointers(pWeight)

        # > Create PDM structure
        if userDistribution is None:
          self.PTB = PDM_part_to_block_create(t_distrib,
                                              t_post,
                                              partActiveNode,
                                              _ln_to_gn,
                                              _weight,
                                              self.pn_elt,
                                              self.n_part,
                                              PDMC)
        else:
          self.PTB = PDM_part_to_block_create_from_distrib(t_distrib,
                                               t_post,
                                               partActiveNode,
                                               _ln_to_gn,
                               <PDM_g_num_t *> userDistribution.data,
                                               self.pn_elt,
                                               self.n_part,
                                               PDMC)
        # Free working arrays
        free(_ln_to_gn)
        if _weight != NULL:
          free(_weight)

    # ------------------------------------------------------------------------
    def exchange_field(self, list part_data, part_stride=1, bint interlaced_str=True):
      """
      Wrapping for PDM_part_to_block_exch : transfert partioned data fields to the
      distribution, allocate and return the distributed array.

      :param self:        PartToBlock object
      :param part_data:   List of partitioned data arrays, each beeing 1 dimensional and with same datatype
      :param part_stride: Stride for partitioned arrays. Can be either a list of n_part array, each element beeing of size
                          pn_elt[i_part] (variable stride will be used) or an integer (cst stride will be used)
      :param interlaced_str: indicate if data are interlaced (True) or interleaved 
      """
      cdef PDM_stride_t _stride_t
      
      cdef int   _part_stride_cst = 0
      cdef int** _part_stride = NULL
      if isinstance(part_stride, int):
        _stride_t = PDM_STRIDE_CST_INTERLACED if interlaced_str else PDM_STRIDE_CST_INTERLEAVED
        _part_stride_cst = part_stride
      elif isinstance(part_stride, list):
        _stride_t = PDM_STRIDE_VAR_INTERLACED
        assert len(part_stride) == self.n_part
        for i_part in range(self.n_part):
          assert_single_dim_np(part_stride[i_part], NPY.int32, self.pn_elt[i_part])
        _part_stride = np_list_to_int_pointers(part_stride)
      else:
        raise ValueError("Invalid stride in PtB exchange")

      assert len(part_data) == self.n_part
      for i_part, p_data in enumerate(part_data):
        expt_size = part_stride[i_part].sum() if _stride_t == PDM_STRIDE_VAR_INTERLACED else part_stride*self.pn_elt[i_part]
        assert_single_dim_np(p_data, part_data[0].dtype, expt_size) #Dtype must be the same for all parts

      cdef void** _part_data   = np_list_to_void_pointers(part_data)

      # Retrieve size of data if proc holds no partition
      ref_dtype = recover_dtype(part_data, self.py_comm)
      cdef size_t s_data         = ref_dtype.itemsize
      cdef size_t dtype_data_num = ref_dtype.num


      cdef int*  _block_stride  = NULL
      cdef void* _block_data    = NULL
      c_size = PDM_part_to_block_exch(self.PTB,
                                      s_data,
                                      _stride_t,
                                      _part_stride_cst,
                                      _part_stride,
                                      _part_data,
                                     &_block_stride,
                                     &_block_data)

      block_data = create_numpy(_block_data, dtype_data_num, c_size)

      if(_stride_t == PDM_STRIDE_VAR_INTERLACED):
        block_stride = create_numpy_i(_block_stride, PDM_part_to_block_n_elt_block_get(self.PTB))
      else:
        block_stride = None

      if _stride_t == PDM_STRIDE_VAR_INTERLACED:
        free(_part_stride)
      free(_part_data)

      return block_stride, block_data
    # ------------------------------------------------------------------------
    def PartToBlock_Exchange(self, dict dField, dict pField, pStrid=1, bint interlaced_str=True):
      """ Shortcut to exchange multiple fieds stored in dict """
      for field_name, part_data in pField.items():
        block_stride, block_data = self.exchange_field(part_data, pStrid, interlaced_str)
        dField[field_name] = block_data
        if block_stride is not None:
          dField[field_name + "#PDM_Stride"] = block_stride

    # ------------------------------------------------------------------------
    def getBlockGnumCopy(self):
      """ Return a copy of the global numbers """
      BlockGnumNPY = create_numpy_g(PDM_part_to_block_block_gnum_get(self.PTB),
                                    PDM_part_to_block_n_elt_block_get(self.PTB),
                                    flag_owndata=False)
      return NPY.copy(BlockGnumNPY)

    # ------------------------------------------------------------------------
    def getBlockGnumCountCopy(self):
      """ Return a copy of the number of occurence of each element """
      BlockGnumCountNPY = create_numpy_i(PDM_part_to_block_block_gnum_count_get(self.PTB),
                                         PDM_part_to_block_n_elt_block_get(self.PTB),
                                         flag_owndata=False)
      return NPY.copy(BlockGnumCountNPY)
    # ------------------------------------------------------------------------
    def getDistributionCopy(self):
      """ Return a copy of the distribution array """
      DistribNPY = create_numpy_g(PDM_part_to_block_distrib_index_get(self.PTB),
                                  self.py_comm.Get_size()+1,
                                  flag_owndata=False)
      return NPY.copy(DistribNPY)

    # ------------------------------------------------------------------------
    def getBeginNbEntryAndGlob(self):
      """ Short cut to get the bounds and size of distribution """
      cdef PDM_g_num_t* Distrib = PDM_part_to_block_distrib_index_get(self.PTB)

      Beg = Distrib[self.py_comm.Get_rank()]
      NbE = Distrib[self.py_comm.Get_rank()+1]-Distrib[self.py_comm.Get_rank()]
      GlB = Distrib[self.py_comm.Get_size()]

      return (Beg, NbE, GlB)

    # ------------------------------------------------------------------------
    def __dealloc__(self):
      """ Free structure """
      # > Free Ppart Structure
      cdef PDM_part_to_block_t* none = PDM_part_to_block_free(self.PTB)
      assert none == NULL
      # > Free allocated array
      free(self.pn_elt)

