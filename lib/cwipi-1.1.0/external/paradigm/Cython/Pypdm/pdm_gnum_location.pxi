# ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# > Wrapping of functions
cdef extern from "pdm_gnum_location.h":

    ctypedef struct PDM_gnum_location_t:
        pass

    PDM_gnum_location_t* PDM_gnum_location_create(const int          n_part_in,
                                                  const int          n_part_out,
                                                  const PDM_MPI_Comm comm,
                                                  PDM_ownership_t owner)

    void           PDM_gnum_location_elements_set(      PDM_gnum_location_t *gnum_loc,
                                                  const int                  i_part_in,
                                                  const int                  n_elts_in,
                                                  const PDM_g_num_t         *gnum_in)

    void PDM_gnum_location_requested_elements_set(      PDM_gnum_location_t *gnum_loc,
                                                  const int                  i_part_out,
                                                  const int                  n_elts_out,
                                                  const PDM_g_num_t         *gnum_out)

    void                PDM_gnum_location_compute(PDM_gnum_location_t *gnum_loc)

    void                    PDM_gnum_location_get(      PDM_gnum_location_t  *gnum_loc,
                                                  const int                   i_part_out,
                                                        int                 **location_idx,
                                                        int                 **location)

    void                   PDM_gnum_location_free(      PDM_gnum_location_t *gnum_loc)

# ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

# ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# > Class definition
cdef class GlobalNumberingLocation:
  """
     GlobalNumberingLocation : Interface for pdm_gnum_location.c
  """
  # --------------------------------------------------------------------------
  # > Class attributesint _id
  cdef PDM_gnum_location_t* _gnum_loc
  cdef NPY.int32_t[:] _n_elmts_in_part
  cdef NPY.int32_t[:] _n_elmts_out_part
  # --------------------------------------------------------------------------

  # --------------------------------------------------------------------------
  def __init__(self, int n_part_in, int n_part_out, MPI.Comm comm):
    """
        Init a gnum location structure
    """
    # ************************************************************************
    # > Declaration
    cdef MPI.MPI_Comm c_comm = comm.ob_mpi
    # ************************************************************************
    # > Init private array storing partition sizes
    self._n_elmts_in_part  = NPY.zeros(n_part_in,  dtype=NPY.int32)
    self._n_elmts_out_part = NPY.zeros(n_part_out, dtype=NPY.int32)
    # ************************************************************************

    # ************************************************************************
    # > PDM call
    self._gnum_loc = PDM_gnum_location_create(n_part_in,
                                              n_part_out,
                                              PDM_MPI_mpi_2_pdm_mpi_comm (<void *> &c_comm),
                                              PDM_OWNERSHIP_UNGET_RESULT_IS_FREE)
    # ************************************************************************

  # --------------------------------------------------------------------------
  def gnum_location_elements_set(self,
                                 int i_part_in,
                                 int n_elmts_in,
                                 NPY.ndarray[npy_pdm_gnum_t, mode='c', ndim=1] gnum_in):
    """
       Calls set method for elements location from PDM_gnum_location
    """
    # ************************************************************************
    # > Declaration
    # ************************************************************************

    # ************************************************************************
    # > PDM call
    PDM_gnum_location_elements_set(self._gnum_loc,
                                   i_part_in,
                                   n_elmts_in,
                                   <PDM_g_num_t*> gnum_in.data)
    # ************************************************************************

    # ************************************************************************
    self._n_elmts_in_part[i_part_in] = n_elmts_in
    # ************************************************************************

  # --------------------------------------------------------------------------
  def gnum_location_requested_elements_set(self,
                                           int i_part_out,
                                           int n_elmts_out,
                                           NPY.ndarray[npy_pdm_gnum_t, mode='c', ndim=1] gnum_out):
    """
       Calls set method for requested elements location from PDM_gnum_location
    """
    # ************************************************************************
    # > Declaration
    # ************************************************************************

    # ************************************************************************
    # > PDM call
    PDM_gnum_location_requested_elements_set(self._gnum_loc,
                                             i_part_out,
                                             n_elmts_out,
                              <PDM_g_num_t*> gnum_out.data)
    # ************************************************************************

    # ************************************************************************
    self._n_elmts_out_part[i_part_out] = n_elmts_out
    # ************************************************************************

  # --------------------------------------------------------------------------
  def gnum_location_compute(self):
    """
       Calls compute method from PDM_gnum_location
    """
    # ************************************************************************
    # > Declaration
    # ************************************************************************

    # ************************************************************************
    # > PDM call
    PDM_gnum_location_compute(self._gnum_loc)
    # ************************************************************************

  # --------------------------------------------------------------------------
  def gnum_location_get(self,
                        int i_part_out):
    """
       Calls get method from PDM_gnum_location
    """
    # ************************************************************************
    # > Declaration
    # ************************************************************************
    cdef int          *location_idx
    cdef int          *location
    cdef NPY.npy_intp  dim
    # ************************************************************************
    # > PDM call
    PDM_gnum_location_get(self._gnum_loc,
                          i_part_out,
                          &location_idx,
                          &location)
    # ************************************************************************

    # ************************************************************************
    if(location_idx == NULL):
      locationIdx = None
    else:
      locationIdx = create_numpy_i(location_idx, self._n_elmts_out_part[i_part_out]+1, flag_owndata=False)
    if(location == NULL):
      locationArr = None
    else:
      dim = <NPY.npy_intp> (location_idx[self._n_elmts_out_part[i_part_out]])
      locationArr = create_numpy_i(location, dim, flag_owndata=False)
    # ************************************************************************

    # ************************************************************************
    return (locationIdx, locationArr)
    # ************************************************************************

  # --------------------------------------------------------------------------
  def __dealloc__(self):
    """
       Calls the free method of PDM_gnum_location
    """
    # ************************************************************************
    # > Declaration
    # ************************************************************************

    # ************************************************************************
    # > PDM call
    # Todo : tenir compte du partial ?
    PDM_gnum_location_free(self._gnum_loc)
    # ************************************************************************

# ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
