
cdef extern from "pdm_closest_points.h":
  # > Wrapping of Ppart Structures
  ctypedef struct PDM_closest_point_t:
    pass
  # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

  PDM_closest_point_t* PDM_closest_points_create(PDM_MPI_Comm    comm,
                                int             n_closest,
                                PDM_ownership_t owner)

  void PDM_closest_points_n_part_cloud_set(PDM_closest_point_t* cls,
                                           int                  n_part_cloud_src,
                                           int                  n_part_cloud_tgt)

  void PDM_closest_points_tgt_cloud_set(PDM_closest_point_t *cls,
                                        int                  i_part,
                                        int                  n_points,
                                        double              *coords,
                                        PDM_g_num_t         *gnum)

  void PDM_closest_points_src_cloud_set(PDM_closest_point_t *cls,
                                        int                  i_part,
                                        int                  n_points,
                                        double              *coords,
                                        PDM_g_num_t         *gnum)

  void PDM_closest_points_compute(PDM_closest_point_t  *cls)

  void PDM_closest_points_get(PDM_closest_point_t  *cls,
                              int                   i_part_tgt,
                              PDM_g_num_t         **closest_src_gnum,
                              double              **closest_src_distance)

  void PDM_closest_points_tgt_in_src_get(PDM_closest_point_t  *cls,
                                         int                   i_part_src,
                                         int                 **tgt_in_src_idx,
                                         PDM_g_num_t         **tgt_in_src)

  void PDM_closest_points_tgt_in_src_dist_get(PDM_closest_point_t  *cls,
                                              int                   i_part_src,
                                              int                 **tgt_in_src_idx,
                                              double              **tgt_in_src_dist)

  void PDM_closest_points_free(PDM_closest_point_t *cls)

  void PDM_closest_points_dump_times(PDM_closest_point_t *cls)

  void PDM_closest_points_part_to_part_get(PDM_closest_point_t  *cls,
                                           PDM_part_to_part_t  **ptp,
                                           PDM_ownership_t       ownership);


# ------------------------------------------------------------------
cdef class ClosestPoints:
  """
  """
  # ************************************************************************
  # > Class attributes
  cdef PDM_closest_point_t *_cls
  cdef int _size
  cdef int _rank
  cdef int* src_n_points
  cdef int* tgt_n_points
  cdef int n_closest
  cdef MPI.Comm py_comm
  # ************************************************************************

  # ------------------------------------------------------------------------
  def __init__(self, MPI.Comm comm,
                     int      n_closest):
    """
    __init__(comm, n_closest)

    Create the structure

    Parameters:
      comm      (MPI.Comm) : MPI communicator
      n_closest (int)      : Number of closest source points to find for each target point
    """
    # ::::::::::::::::::::::::::::::::::::::::::::::::::
    self.n_closest = n_closest
    # ::::::::::::::::::::::::::::::::::::::::::::::::::

    # ::::::::::::::::::::::::::::::::::::::::::::::::::
    self.py_comm = comm
    self._rank = comm.Get_rank()
    self._size = comm.Get_size()
    # ::::::::::::::::::::::::::::::::::::::::::::::::::

    # ::::::::::::::::::::::::::::::::::::::::::::::::::
    # > Convert mpi4py -> PDM_MPI
    cdef MPI.MPI_Comm c_comm = comm.ob_mpi
    cdef PDM_MPI_Comm PDMC   = PDM_MPI_mpi_2_pdm_mpi_comm(<void *> &c_comm)
    # ::::::::::::::::::::::::::::::::::::::::::::::::::

    # ::::::::::::::::::::::::::::::::::::::::::::::::::
    self.src_n_points = NULL
    self.tgt_n_points = NULL
    # ::::::::::::::::::::::::::::::::::::::::::::::::::

    # ::::::::::::::::::::::::::::::::::::::::::::::::::
    self._cls = PDM_closest_points_create(PDMC, n_closest, PDM_OWNERSHIP_USER) # Python take ownership
    # ::::::::::::::::::::::::::::::::::::::::::::::::::

  # ------------------------------------------------------------------------
  def n_part_cloud_set(self, int n_part_cloud_src,
                             int n_part_cloud_tgt):
    """
    n_part_cloud_set(n_part_cloud_src, n_part_cloud_tgt)

    Set the number of partitions of both point clouds

    Parameters:
      n_part_cloud_src (int) : Number of partitions of the source cloud
      n_part_cloud_tgt (int) : Number of partitions of the target cloud
    """
    assert(self.src_n_points  == NULL)
    assert(self.tgt_n_points  == NULL)
    self.src_n_points = <int *> malloc(sizeof(int *) * n_part_cloud_src )
    self.tgt_n_points = <int *> malloc(sizeof(int *) * n_part_cloud_tgt )
    PDM_closest_points_n_part_cloud_set(self._cls,
                                        n_part_cloud_src,
                                        n_part_cloud_tgt)

  # ------------------------------------------------------------------------
  def tgt_cloud_set(self, int                                           i_part,
                          int                                           n_points,
                          NPY.ndarray[NPY.double_t  , mode='c', ndim=1] coords,
                          NPY.ndarray[npy_pdm_gnum_t, mode='c', ndim=1] gnum):
    """
    tgt_cloud_set(i_part, n_points, coords, gnum)

    Set the target point cloud

    Parameters:
      i_part   (int)                        : Partition identifier
      n_points (int)                        : Number of points (to remove)
      coords   (np.ndarray[np.double_t])    : Point coordinates
      gnum     (np.ndarray[npy_pdm_gnum_t]) : Point global ids
    """
    self.tgt_n_points[i_part] = n_points
    PDM_closest_points_tgt_cloud_set(self._cls,
                                     i_part,
                                     n_points,
                           <double*> coords.data,
                      <PDM_g_num_t*> gnum.data)


  # ------------------------------------------------------------------------
  def src_cloud_set(self, int                                           i_part,
                          int                                           n_points,
                          NPY.ndarray[NPY.double_t  , mode='c', ndim=1] coords,
                          NPY.ndarray[npy_pdm_gnum_t, mode='c', ndim=1] gnum):
    """
    src_cloud_set(i_part, n_points, coords, gnum)

    Set the source point cloud

    Parameters:
      i_part   (int)                        : Partition identifier
      n_points (int)                        : Number of points (to remove)
      coords   (np.ndarray[np.double_t])    : Point coordinates
      gnum     (np.ndarray[npy_pdm_gnum_t]) : Point global ids
    """
    self.src_n_points[i_part] = n_points
    PDM_closest_points_src_cloud_set(self._cls,
                                     i_part,
                                     n_points,
                           <double*> coords.data,
                      <PDM_g_num_t*> gnum.data)

  # ------------------------------------------------------------------------
  def compute(self):
    """
    Look for closest points
    """
    PDM_closest_points_compute(self._cls)

  # ------------------------------------------------------------------------
  def points_get(self, int i_part_tgt):
    """
    points_get(i_part_tgt)

    Get closest source points global ids and (squared) distance

    Parameters:
      i_part_tgt (int) : Partition identifier

    Returns:
      Dictionary
        - ``"closest_src_gnum"``     (`np.ndarray[npy_pdm_gnum_t]`) : Global ids of the closest source points
        - ``"closest_src_distance"`` (`np.ndarray[np.double_t]`)    : Squared distance from closest source points
    """
    # ************************************************************************
    # > Declaration
    cdef PDM_g_num_t *closest_src_gnum
    cdef double      *closest_src_distance
    # ************************************************************************

    # > Get
    PDM_closest_points_get(self._cls,
                           i_part_tgt,
                           &closest_src_gnum,
                           &closest_src_distance)

    dim = self.tgt_n_points[i_part_tgt] * self.n_closest
    
    return {'closest_src_gnum'     : create_numpy_g(closest_src_gnum, dim),
            'closest_src_distance' : create_numpy_d(closest_src_distance, dim)}

  # ------------------------------------------------------------------------
  def tgt_in_src_get(self, int i_part_src):
    """
    tgt_in_src_get(i_part_src)

    Get source-to-target mapping and distances

    Parameters:
      i_part_src (int) : Partition identifier

    Returns:
      Dictionary
        - ``"tgt_in_src_idx"``   (`np.ndarray[np.int32_t]`)     : Index for source-to-target mapping
        - ``"tgt_in_src"``       (`np.ndarray[npy_pdm_gnum_t]`) : Source-to-target mapping (global ids)
        - ``"tgt_in_src_dist2"`` (`np.ndarray[np.double_t]`)    : Squared source-to-target distance
    """
    # ************************************************************************
    # > Declaration
    cdef int *tgt_in_src_idx
    cdef PDM_g_num_t *tgt_in_src
    cdef double *tgt_in_src_dist
    # ************************************************************************

    # > Get
    PDM_closest_points_tgt_in_src_get(self._cls,
                                      i_part_src,
                                      &tgt_in_src_idx,
                                      &tgt_in_src)

    PDM_closest_points_tgt_in_src_dist_get(self._cls,
                                           i_part_src,
                                           &tgt_in_src_idx,
                                           &tgt_in_src_dist)

    src_n_pts = self.src_n_points[i_part_src]
    return {
            'tgt_in_src_idx'  : create_numpy_i(tgt_in_src_idx, src_n_pts+1),
            'tgt_in_src'      : create_numpy_g(tgt_in_src, tgt_in_src_idx[src_n_pts]),
            'tgt_in_src_dist2': create_numpy_d(tgt_in_src_dist, tgt_in_src_idx[src_n_pts])
            }

  # ------------------------------------------------------------------------
  def part_to_part_get(self):
    """
    Get the PartToPart object to exchange data between
    the source and target point clouds

    Returns:
      PartToPart object (:py:class:`PartToPart`)
    """
    cdef PDM_part_to_part_t *ptpc
    PDM_closest_points_part_to_part_get(self._cls,
                                        &ptpc,
                                        PDM_OWNERSHIP_USER)

    py_caps = PyCapsule_New(ptpc, NULL, NULL)
    return PartToPartCapsule(py_caps, self.py_comm) # The free is inside the class

  # ------------------------------------------------------------------------
  def dump_times(self):
    """
    Dump elapsed and CPU times
    """
    PDM_closest_points_dump_times(self._cls)

  # ------------------------------------------------------------------------
  def __dealloc__(self):
    """
    """
    free(self.tgt_n_points)
    free(self.src_n_points)
    PDM_closest_points_free(self._cls)
