cdef extern from "pdm_block_to_part.h":
    # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    # > Wrapping of Ppart Structure
    ctypedef struct PDM_block_to_part_t:
      pass
    # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    # > Wrapping of function
    PDM_block_to_part_t *PDM_block_to_part_create(PDM_g_num_t   *blockDistribIdx,
                                                  PDM_g_num_t  **gnum_elt,
                                                  int           *n_elt,
                                                  int            n_part,
                                                  PDM_MPI_Comm   comm)

    void PDM_block_to_part_exch_in_place(PDM_block_to_part_t  *btp,
                                size_t                s_data,
                                PDM_stride_t          t_stride,
                                int                  *block_stride,
                                void                 *block_data,
                                int                 **part_stride,
                                void                **part_data)

    void PDM_block_to_part_exch(PDM_block_to_part_t  *btp,
                                size_t                s_data,
                                PDM_stride_t          t_stride,
                                int                  *block_stride,
                                void                 *block_data,
                                int                 ***part_stride,
                                void                ***part_data)

    PDM_block_to_part_t *PDM_block_to_part_free(PDM_block_to_part_t *btp)
    # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::


# ------------------------------------------------------------------
cdef class BlockToPart:
    """
       BlockToPart: Interface for block_to_part.c
    """
    # ************************************************************************
    # > Class attributes
    cdef PDM_block_to_part_t* BTP
    cdef int                  n_part
    cdef int                  dn_elt
    cdef int*                 pn_elt
    cdef MPI.Comm             py_comm
    # ************************************************************************
    # ------------------------------------------------------------------------
    def __cinit__(self, NPY.ndarray[npy_pdm_gnum_t, ndim=1, mode='c'] Distrib,
                        MPI.Comm comm,
                        list     pLNToGN,
                        int      partN):
        """
        Constructor of BlockToPart object : Python wrapping of PDM library (E. Quémerais)

            :param comm:     MPI Communicator (Caution MPI Comm is a mpi4py object )
            :param Distrib:  Distribution of distribute array (Size = nRank+1)
            :param pLNToGN:  Part list containaing numpy on LNToGN for each partition (len = partN)
            :param partN:    Number of partitions

        """
        # > Some checks
        assert(len(pLNToGN) == partN)
        for i in range(partN):
          assert_single_dim_np(pLNToGN[i], npy_pdm_gnum_dtype)
        assert Distrib.size == comm.Get_size() + 1

        # > Store class parameters
        self.n_part  = partN
        self.dn_elt  = Distrib[comm.Get_rank() + 1] - Distrib[comm.Get_rank()]
        self.pn_elt = list_to_int_pointer([array.size for array in pLNToGN])
        self.py_comm = comm

        # > Convert input data
        cdef MPI.MPI_Comm c_comm = comm.ob_mpi
        cdef PDM_MPI_Comm PDMC   = PDM_MPI_mpi_2_pdm_mpi_comm(&c_comm)

        cdef PDM_g_num_t*  _distrib  = <PDM_g_num_t *> Distrib.data
        cdef PDM_g_num_t** _ln_to_gn = np_list_to_gnum_pointers(pLNToGN)

        # > Create PDM structure
        self.BTP = PDM_block_to_part_create(_distrib,
                     <const PDM_g_num_t **> _ln_to_gn,
                                            self.pn_elt,
                                            self.n_part,
                                            PDMC)
        # Free working arrays
        free(_ln_to_gn)

    # ------------------------------------------------------------------
    def exchange_field(self, NPY.ndarray block_data, block_stride=1, bint interlaced_str=True):
      """
      Wrapping for PDM_block_to_part_exch : transfert a distributed data field to the
      partitions, allocate and return the partitionned array list.

      :param self:         BlockToPart object
      :param block_data:   Distributed data array, 1 dimensional and with same datatype for each rank
      :param block_stride: Stride for distributed array. Can be either an array of size dn_elt (variable
                           stride will be used) or an integer (cst stride will be used)
      :param interlaced_str: indicate if data are interlaced (True) or interleaved 
      """

      cdef NPY.ndarray[NPY.int32_t, ndim=1, mode='c'] numpy_int
      cdef PDM_stride_t _stride_t
      cdef int*         _block_stride

      if isinstance(block_stride, int):
        _stride_t = PDM_STRIDE_CST_INTERLACED if interlaced_str else PDM_STRIDE_CST_INTERLEAVED
        _block_stride = <int *> malloc(1 * sizeof(int))
        _block_stride[0] = block_stride
        assert_single_dim_np(block_data, block_data.dtype, block_stride*self.dn_elt)
      elif isinstance(block_stride, NPY.ndarray):
        _stride_t = PDM_STRIDE_VAR_INTERLACED
        assert_single_dim_np(block_stride, NPY.int32, self.dn_elt)
        assert_single_dim_np(block_data, block_data.dtype, block_stride.sum())
        numpy_int = block_stride
        _block_stride = <int *> numpy_int.data
      else:
        raise ValueError("Invalid stride in BtP exchange")

      cdef size_t s_data = block_data.dtype.itemsize #Should be the same for all procs
      assert self.py_comm.allreduce(s_data, op=MPI.MIN) == self.py_comm.allreduce(s_data, op=MPI.MAX) == s_data

      cdef int**  _part_stride = NULL
      cdef void** _part_data   = NULL
      PDM_block_to_part_exch(self.BTP,
                              s_data,
                              _stride_t,
                             _block_stride,
                     <void *> block_data.data,
                              &_part_stride,
                              &_part_data)


      part_stride = []
      part_data   = []
      cdef NPY.npy_intp dim_np

      for i_part in range(self.n_part):

        if _stride_t == PDM_STRIDE_VAR_INTERLACED:
          stride = create_numpy_i(_part_stride[i_part], self.pn_elt[i_part])
          part_stride.append(stride)
          dim_np = <NPY.npy_intp> stride.sum()
        else:
          dim_np = <NPY.npy_intp> _block_stride[0] * self.pn_elt[i_part]

        part_data.append(create_numpy(_part_data[i_part], block_data.dtype.num, dim_np))

      if _stride_t == PDM_STRIDE_VAR_INTERLACED:
        free(_part_stride)
      else:
        free(_block_stride)
        part_stride = None
      free(_part_data)

      return part_stride, part_data

    # ------------------------------------------------------------------
    def exchange_field_inplace(self, NPY.ndarray block_data, list part_data, 
      block_stride=1, list part_stride=None, bint interlaced_str=True):
      """
      Wrapping for PDM_block_to_part_exch_in_place : transfert a distributed data field to the
      partitions. Fill the pre-allocated partitionned arrays

      :param self:         BlockToPart object
      :param block_data:   Distributed data array, 1 dimensional and with same datatype for each rank
      :param part_data:    List of the partN pre allocated partitionned data array, each one beeing 1 dimensional 
                           and with same datatype than block_data
      :param block_stride: Stride for distributed array. Can be either an array of size dn_elt (variable
                           stride will be used) or an integer (cst stride will be used)
      :param part_stride:  List of the partN pre allocated parititioned data strides
      :param interlaced_str: indicate if data are interlaced (True) or interleaved 
      """
      cdef NPY.ndarray[NPY.int32_t, ndim=1, mode='c'] numpy_int
      cdef PDM_stride_t _stride_t
      cdef int*         _block_stride
      cdef int**        _part_stride = NULL

      assert len(part_data) == self.n_part
      if isinstance(block_stride, int):
        _stride_t = PDM_STRIDE_CST_INTERLACED if interlaced_str else PDM_STRIDE_CST_INTERLEAVED
        _block_stride = <int *> malloc(1 * sizeof(int *))
        _block_stride[0] = block_stride
        assert_single_dim_np(block_data, block_data.dtype, block_stride*self.dn_elt)
      elif isinstance(block_stride, NPY.ndarray):
        _stride_t = PDM_STRIDE_VAR_INTERLACED
        assert_single_dim_np(block_stride, NPY.int32, self.dn_elt)
        assert_single_dim_np(block_data, block_data.dtype, block_stride.sum())
        numpy_int = block_stride
        _block_stride = <int *> numpy_int.data
        assert len(part_stride) == self.n_part
        for i_part in range(self.n_part):
          assert_single_dim_np(part_stride[i_part], NPY.int32, self.pn_elt[i_part])
        _part_stride = np_list_to_int_pointers(part_stride)
      else:
        raise ValueError("Invalid stride in BtP exchange")

      cdef size_t s_data = block_data.dtype.itemsize #Should be the same for all procs
      assert self.py_comm.allreduce(s_data, op=MPI.MIN) == self.py_comm.allreduce(s_data, op=MPI.MAX) == s_data

      for i_part in range(self.n_part):
        expt_size = part_stride[i_part].sum() if _stride_t == PDM_STRIDE_VAR_INTERLACED else block_stride*self.pn_elt[i_part]
        assert_single_dim_np(part_data[i_part], block_data.dtype, expt_size)

      cdef void** _part_data   = np_list_to_void_pointers(part_data)
      PDM_block_to_part_exch_in_place(self.BTP,
                             s_data,
                             _stride_t,
                            _block_stride,
                    <void *> block_data.data,
                             _part_stride,
                             _part_data)
      if _stride_t != PDM_STRIDE_VAR_INTERLACED:
        free(_block_stride)
      free(_part_data)

    # ------------------------------------------------------------------
    def BlockToPart_Exchange(self, dict         dField,
                                   dict         pField,
                                   BlkStride = 1,
                                   bint interlaced_str=True):
      """ Shortcut to exchange multiple fieds stored in dict """
      dField = {key:data for key,data in dField.items() if not key.endswith("#PDM_Stride")}
      for field_name in dField:
        PartStride = None
        if field_name + '#PDM_Stride' in pField:
          PartStride = pField[field_name + '#PDM_Stride']
        self.exchange_field_inplace(dField[field_name], pField[field_name], BlkStride, PartStride, interlaced_str)

    # ------------------------------------------------------------------
    def BlockToPart_Exchange2(self, dict         dField,
                                    dict         pField,
                                    BlkStride = 1,
                                    bint interlaced_str=True):

      """ Shortcut to exchange multiple fieds stored in dict """
      dField = {key:data for key,data in dField.items() if not key.endswith("#PDM_Stride")}
      for field_name, field in dField.items():
        part_stride, part_data = self.exchange_field(field, BlkStride, interlaced_str)
        pField[field_name] = part_data
        if part_stride is not None:
          pField[field_name + "#PDM_Stride"] = part_stride


    # ------------------------------------------------------------------
    def __dealloc__(self):
      """
      Deallocate all the array
      """
      # > Free PDM Structure
      cdef PDM_block_to_part_t *none = PDM_block_to_part_free(self.BTP)
      assert none == NULL

      # > Free allocated array
      free(self.pn_elt)

