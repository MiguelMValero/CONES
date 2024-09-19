
cdef extern from "pdm_distant_neighbor.h":

    # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    # > Wrapping of Ppart Structure
    ctypedef struct PDM_distant_neighbor_t:
      pass
    # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    # > Wrapping of function
    PDM_distant_neighbor_t* PDM_distant_neighbor_create(PDM_MPI_Comm   comm,
                                                        int            n_part,
                                                        int           *n_entity,
                                                        int          **neighbor_idx,
                                                        int          **neighbor_desc)

    void PDM_distant_neighbor_free(PDM_distant_neighbor_t *dn)

    void PDM_distant_neighbor_exch(PDM_distant_neighbor_t   *dn,
                                   size_t                    s_data,
                                   PDM_stride_t              t_stride,
                                   int                       cst_stride,
                                   int                     **send_entity_stride,
                                   void                    **send_entity_data,
                                   int                    ***recv_entity_stride,
                                   void                   ***recv_entity_data)
    # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::


# ------------------------------------------------------------------
cdef class DistantNeighbor:
    """
       DistantNeighbor: Interface for pdm_distant_neighbor.c
    """
    # ************************************************************************
    # > Class attributes
    cdef PDM_distant_neighbor_t *_dn
    cdef int   _n_part
    cdef int  *_n_entity
    cdef int **_neighbor_idx
    cdef int **_neighbor_desc
    # ************************************************************************
    # ------------------------------------------------------------------------
    def __cinit__(self, MPI.Comm comm,
                        int      n_part,
                        list     neighbor_idx,
                        list     neighbor_desc):
        """
        TODOUX
        """
        # ************************************************************************
        # > Declaration
        cdef NPY.ndarray[NPY.int32_t, ndim=1, mode='fortran'] pneighbor_idx
        cdef NPY.ndarray[NPY.int32_t, ndim=1, mode='fortran'] pneighbor_desc
        # ************************************************************************

        # ::::::::::::::::::::::::::::::::::::::::::::::::::
        # > Some verification
        assert(len(neighbor_idx) == len(neighbor_desc) == n_part)
        # ::::::::::::::::::::::::::::::::::::::::::::::::::

        # ::::::::::::::::::::::::::::::::::::::::::::::::::
        # > Store partN and parameter
        self._n_part        = n_part
        self._neighbor_idx  = <int **> malloc(self._n_part * sizeof(int **))
        self._neighbor_desc = <int **> malloc(self._n_part * sizeof(int **))
        self._n_entity      = <int *>  malloc(self._n_part * sizeof(int * ))
        # ::::::::::::::::::::::::::::::::::::::::::::::::::

        # ::::::::::::::::::::::::::::::::::::::::::::::::::
        # > Convert mpi4py -> PDM_MPI
        cdef MPI.MPI_Comm c_comm = comm.ob_mpi
        cdef PDM_MPI_Comm PDMC   = PDM_MPI_mpi_2_pdm_mpi_comm(&c_comm)
        # ::::::::::::::::::::::::::::::::::::::::::::::::::

        # ::::::::::::::::::::::::::::::::::::::::::::::::::
        # > Prepare
        for i_part in range(self._n_part):

          # ------------------------------------------------
          # > Get shape of array
          n_entity       = neighbor_idx[i_part].shape[0]-1
          pneighbor_idx  = neighbor_idx[i_part]
          pneighbor_desc = neighbor_desc[i_part]

          self._n_entity[i_part]      = <int  > n_entity
          self._neighbor_idx[i_part]  = <int *> pneighbor_idx.data
          self._neighbor_desc[i_part] = <int *> pneighbor_desc.data
          # ------------------------------------------------

        # ::::::::::::::::::::::::::::::::::::::::::::::::::
        # > Create
        self._dn = PDM_distant_neighbor_create(PDMC,
                                                 n_part,
                                                 self._n_entity,
                                                 self._neighbor_idx,
                                                 self._neighbor_desc)
        # ::::::::::::::::::::::::::::::::::::::::::::::::::

    # ------------------------------------------------------------------------
    def DistantNeighbor_Exchange(self, list         l_send_entity_data,
                                       list         l_recv_entity_data,
                                       int          cst_stride = 1,
                                       list         l_send_entity_stri = None,
                                       list         l_recv_entity_stri = None):
        """
        """
        # ************************************************************************
        # > Declaration
        cdef int           i
        cdef int           idx
        cdef int           i_part

        cdef NPY.ndarray psend_entity_data
        cdef NPY.ndarray[NPY.int32_t, ndim=1, mode='fortran'] psend_entity_stri
        cdef NPY.ndarray precv_entity_data
        cdef NPY.ndarray[NPY.int32_t, ndim=1, mode='fortran'] precv_entity_stri

        # > For PDM
        cdef int         **send_entity_stri
        cdef void        **send_entity_data
        cdef int         **recv_entity_stri
        cdef void        **recv_entity_data
        cdef size_t        s_data
        cdef int           stride_one
        cdef int           size_data
        cdef int           npyflags=-1
        cdef NPY.ndarray   tmp_data
        cdef PDM_stride_t  t_stride
        # ************************************************************************

        # ::::::::::::::::::::::::::::::::::::::::::::::::::
        assert(len(l_send_entity_data) == self._n_part)
        assert(len(l_recv_entity_data) == 0           )
        t_stride = PDM_STRIDE_CST_INTERLACED
        if(l_send_entity_stri is not None):
          t_stride = PDM_STRIDE_VAR_INTERLACED
        # ::::::::::::::::::::::::::::::::::::::::::::::::::

        # ::::::::::::::::::::::::::::::::::::::::::::::::::
        # print(" ----> flag 1")
        send_entity_stri  = <int ** > malloc(self._n_part * sizeof(int  **))
        send_entity_data  = <void **> malloc(self._n_part * sizeof(void **))
        # ::::::::::::::::::::::::::::::::::::::::::::::::::

        # ::::::::::::::::::::::::::::::::::::::::::::::::::
        # print(" ----> flag 2")
        for i_part in range(self._n_part):
          psend_entity_data = l_send_entity_data[i_part]
          send_entity_data[i_part] = <void*> psend_entity_data.data
          if(l_send_entity_stri is not None):
            psend_entity_stri = l_send_entity_stri[i_part]
            send_entity_stri[i_part] = <int*> psend_entity_stri.data
          s_data     = psend_entity_data.dtype.itemsize
          dtype_data = psend_entity_data.dtype.num
          dtypep     = psend_entity_data.dtype
        # ::::::::::::::::::::::::::::::::::::::::::::::::::

        #  = <PDM_stride_t> (0),
        # ::::::::::::::::::::::::::::::::::::::::::::::::::
        # print(" ----> flag 3 :: PDM_distant_neighbor_exch")
        PDM_distant_neighbor_exch(self._dn, s_data, t_stride,
                                  cst_stride,
                                  send_entity_stri,
                                  send_entity_data,
                                  &recv_entity_stri,
                                  &recv_entity_data)
        # print(" ----> flag 3 :: PDM_distant_neighbor_exch end ")
        # ::::::::::::::::::::::::::::::::::::::::::::::::::

        # ::::::::::::::::::::::::::::::::::::::::::::::::::
        # > Take owership with numpy and create
        for i_part in range(self._n_part):
          if(l_send_entity_stri is not None):
            size_data = 0
            for i in xrange(self._n_entity[i_part]):
              size_data += recv_entity_stri[i_part][i]
            tmp_data = create_numpy_i(recv_entity_stri[i_part], self._n_entity[i_part])
            l_recv_entity_stri.append(tmp_data)
          else:
            size_data = cst_stride * self._neighbor_idx[i_part][self._n_entity[i_part]]
          tmp_data = create_numpy(recv_entity_data[i_part], dtype_data, size_data)
          l_recv_entity_data.append(tmp_data)
        # ::::::::::::::::::::::::::::::::::::::::::::::::::

        # ::::::::::::::::::::::::::::::::::::::::::::::::::
        # print(" ----> flag 4 :: free")
        free(send_entity_stri)
        free(send_entity_data)
        if(l_send_entity_stri is not None):
          free(recv_entity_stri)
        free(recv_entity_data)
        # print(" ----> flag 4 :: free end ")
        # ::::::::::::::::::::::::::::::::::::::::::::::::::

    # ------------------------------------------------------------------------
    def __dealloc__(self):
      """
         Use the free method of PDM Lib
      """
      # ************************************************************************
      # > Declaration
      # ************************************************************************

      # > Free Ppart Structure
      PDM_distant_neighbor_free(self._dn)

      free(self._n_entity)
      free(self._neighbor_idx)
      free(self._neighbor_desc)
