import sys

cdef extern from "pdm_gnum.h":
  # ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  # > Wrapping of structure
  ctypedef struct PDM_gen_gnum_t:
    pass
  # ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

  # ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  # > Wrapping of functions
  PDM_gen_gnum_t * PDM_gnum_create(const int             dim,
                                   const int             n_part,
                                   const PDM_bool_t      merge,
                                   const double          tolerance,
                                   const PDM_MPI_Comm    comm,
                                   const PDM_ownership_t owner)

  void PDM_gnum_set_from_coords(PDM_gen_gnum_t *gen_gnum,
                                const int       i_part,
                                const int       n_elts,
                                const double   *coords,
                                const double   *char_length)

  void PDM_gnum_set_from_parents(PDM_gen_gnum_t    *gen_gnum,
                                 const int          i_part,
                                 const int          n_elts,
                                 const PDM_g_num_t *parent_gnum)

  void PDM_gnum_set_parents_nuplet(PDM_gen_gnum_t  *gen_gnum,
                                   const int        nuplet)


  void PDM_gnum_compute(PDM_gen_gnum_t *gen_gnum)

  PDM_g_num_t * PDM_gnum_get(PDM_gen_gnum_t *gen_gnum,
                             const int       i_part)

  void PDM_gnum_free(PDM_gen_gnum_t *gen_gnum)
# ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

# ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# > Class definition
cdef class GlobalNumbering:
  """

  """
  # --------------------------------------------------------------------------
  # > Class attributes
  cdef PDM_gen_gnum_t* _gen_gnum
  cdef NPY.npy_intp[:] _n_elem_per_part
  # --------------------------------------------------------------------------

  # --------------------------------------------------------------------------
  def __init__(self,
               int        dim,
               int        n_part,
               PDM_bool_t merge,
               double     tolerance,
               MPI.Comm   comm):
    """
    __init__(dim, n_part, merge, tolerance, comm)
    Create a structure to build a global numbering

    Parameters:
      dim       (int)        : Spatial dimension
      n_part    (int)        : Number of local partitions
      merge     (PDM_bool_t) : Merge double points or not
      tolerance (double)     : Geometric tolerance (used if ``merge`` is set to ``PDM_TRUE``)
      comm      (MPI.Comm)   : MPI communicator
    """
    # ************************************************************************
    # > Declaration
    cdef MPI.MPI_Comm c_comm = comm.ob_mpi
    # ************************************************************************

    # ************************************************************************
    # > Init private array storing partition sizes
    self._n_elem_per_part = NPY.zeros(n_part, dtype=NPY.intp)
    # ************************************************************************

    # ************************************************************************
    # > PDM call
    self._gen_gnum = PDM_gnum_create(dim,
                                     n_part,
                                     merge,
                                     tolerance,
                                     PDM_MPI_mpi_2_pdm_mpi_comm (<void *> &c_comm),
                                     PDM_OWNERSHIP_USER) # Python takes ownership);
    # ************************************************************************

  # --------------------------------------------------------------------------
  def set_from_coords(self,
                      int i_part,
                      NPY.ndarray[NPY.double_t  , mode='c', ndim=1] coords not None,
                      NPY.ndarray[NPY.double_t  , mode='c', ndim=1] char_length):
    """
    set_from_coords(i_part, coords, char_length)
    Set from coordinates

    Note:
      The ordering is based on the `Morton space-filling curve <https://en.wikipedia.org/wiki/Z-order_curve>`_.
      Elements are expected to be unique (i.e. not duplicated on 2 or more ranks).
      If two elements share the same Morton code, their global
      id will be determined by lexicographical ordering of coordinates.

    Warning:
      Coordinates are assumed to be 3-dimensional, even if a lower dimension was specified when creating the ``GlobalNumbering`` instance.

    Parameters:
      i_part      (int)                        : Current partition
      coords      (np.ndarray[np.double_t]   ) : Coordinates (size = 3 * *n_elts*)
      char_length (np.ndarray[npy_pdm_gnum_t]) : Characteristic length (or *None*, used only if ``merge`` is set to ``PDM_TRUE``)
    """
    # ************************************************************************
    # > Declaration
    cdef double *caracteristic_length_data
    cdef double *coords_data
    # ************************************************************************

    # ************************************************************************
    coords_data = <double *> coords.data
    if (char_length is None):
      caracteristic_length_data = NULL
    else:
      caracteristic_length_data = <double *> char_length.data
    # ************************************************************************

    # ************************************************************************
    # > Store size to use it in the get
    if len(coords)%3 != 0:
      print(f"gnum_set_from_coords: Error: coordinates must be provided in 3 dimensions",
            flush=True)
      sys.exit(1)
    cdef int n_elts = len(coords)//3
    self._n_elem_per_part[i_part] = n_elts
    # ************************************************************************

    # ************************************************************************
    # > PDM call
    PDM_gnum_set_from_coords(self._gen_gnum,
                             i_part,
                             n_elts,
                             coords_data,
                             caracteristic_length_data)
    # ************************************************************************

  # --------------------------------------------------------------------------
  def set_from_parent(self,
                      int i_part,
                      NPY.ndarray[npy_pdm_gnum_t  , mode='c', ndim=1] parent_gnum not None):
    """
    set_from_parent(i_part, n_elts, parent_gnum)

    Set parent global numbering

    Parameters:
      i_part      (int)                    : Current partition
      parent_gnum (np.ndarray[np.double_t) : Parent global ids
    """
    # ************************************************************************
    # > Declaration
    # ************************************************************************

    # ************************************************************************
    # > Store size to use it in the get
    cdef int n_elts = len(parent_gnum)
    self._n_elem_per_part[i_part] = n_elts
    # ************************************************************************

    # ************************************************************************
    # > PDM call
    PDM_gnum_set_from_parents(self._gen_gnum,
                              i_part,
                              n_elts,
               <PDM_g_num_t*> parent_gnum.data)
    # ************************************************************************


  # --------------------------------------------------------------------------
  def set_parents_nuplet(self,
                         int nuplet):
    """
    set_parents_nuplet(nuplet)

    Set size of tuple for nuplet

    Parameters:
      nuplet (int) : Size of tuple
    """
    PDM_gnum_set_parents_nuplet(self._gen_gnum,
                                nuplet)

  # --------------------------------------------------------------------------
  def compute(self):
    """
    Build global numbering
    """
    PDM_gnum_compute(self._gen_gnum)

  # --------------------------------------------------------------------------
  def get(self, int i_part):
    """
    get(i_part)

    Get global ids for a given partition

    Parameters:
      i_part (int) : Current partition

    Returns:
      Array of global ids (`np.ndarray[npy_pdm_gnum_t]`)
    """
    # ************************************************************************
    # > Declaration
    cdef PDM_g_num_t *gnum_array
    # ************************************************************************

    # ************************************************************************
    # > PDM call
    gnum_array = PDM_gnum_get(self._gen_gnum, i_part)
    # ************************************************************************

    # ************************************************************************
    if (gnum_array == NULL):
      return None
    else:
      np_gnum_array = create_numpy_g(gnum_array, self._n_elem_per_part[i_part])
    return np_gnum_array
    # ************************************************************************

  # --------------------------------------------------------------------------
  def __dealloc__(self):
    """
    Free a ``GlobalNumbering`` instance
    """
    PDM_gnum_free(self._gen_gnum);

# ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
