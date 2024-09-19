
cdef extern from "pdm_mesh_location.h":
  # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  # > Wrapping of Ppart Structure
  ctypedef struct PDM_mesh_location_t:
    pass
  # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

  # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  # > Wrapping for C
  ctypedef enum PDM_interpolate_kind_t:
    PDM_INTERPOLATE_KIND_FROM_CENTER  = 0
    PDM_INTERPOLATE_KIND_FROM_NODE = 1
  # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

  # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  PDM_mesh_location_t* PDM_mesh_location_create(int               _n_point_cloud,
                                                PDM_MPI_Comm      comm,
                                                PDM_ownership_t owner)
  # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

  # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  void PDM_mesh_location_n_part_cloud_set(PDM_mesh_location_t* ml,
                                          int                  i_point_cloud,
                                          int                  n_part)
  # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

  # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  void PDM_mesh_location_cloud_set(PDM_mesh_location_t *ml,
                                   int                  i_point_cloud,
                                   int                  i_part,
                                   int                  n_points,
                                   double              *coords,
                                   PDM_g_num_t         *gnum)

  void PDM_mesh_location_cloud_get (PDM_mesh_location_t  *ml,
                                    int                   i_point_cloud,
                                    int                   i_part,
                                    int                  *n_points,
                                    double              **coords,
                                    PDM_g_num_t         **gnum)

  # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

  # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  # void PDM_mesh_location_shared_nodal_mesh_set(int  id, PDM_Mesh_nodal_t *mesh_nodal)
  void PDM_mesh_location_mesh_n_part_set (PDM_mesh_location_t  *ml, int  n_part)
  # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

  # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  void PDM_mesh_location_part_set(PDM_mesh_location_t  *ml,
                                  int                  i_part,
                                  int                  n_cell,
                                  int                 *cell_face_idx,
                                  int                 *cell_face,
                                  PDM_g_num_t         *cell_ln_to_gn,
                                  int                  n_face,
                                  int                 *face_vtx_idx,
                                  int                 *face_vtx,
                                  PDM_g_num_t         *face_ln_to_gn,
                                  int                  n_vtx,
                                  double              *coords,
                                  PDM_g_num_t         *vtx_ln_to_gn)
  # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

  # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  void PDM_mesh_location_nodal_part_set(      PDM_mesh_location_t *ml,
                                        const int                  i_part,
                                        const int                  n_cell,
                                        const int                 *cell_vtx_idx,
                                        const int                 *cell_vtx,
                                        const PDM_g_num_t         *cell_ln_to_gn,
                                        const int                  n_vtx,
                                        const double              *coords,
                                        const PDM_g_num_t         *vtx_ln_to_gn)
  # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

  # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  void PDM_mesh_location_part_set_2d(PDM_mesh_location_t *ml,
                                     int                  i_part,
                                     int                  n_face,
                                     int                 *face_edge_idx,
                                     int                 *face_edge,
                                     PDM_g_num_t         *face_ln_to_gn,
                                     int                  n_edge,
                                     int                 *edge_vtx,
                                     int                  n_vtx,
                                     double              *coords,
                                     PDM_g_num_t         *vtx_ln_to_gn)
  # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

  # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  void PDM_mesh_location_nodal_part_set_2d(      PDM_mesh_location_t *ml,
                                           const int                  i_part,
                                           const int                  n_face,
                                           const int                 *face_vtx_idx,
                                           const int                 *face_vtx,
                                           const PDM_g_num_t         *face_ln_to_gn,
                                           const int                  n_vtx,
                                           const double              *coords,
                                           const PDM_g_num_t         *vtx_ln_to_gn)
  # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

  # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  void PDM_mesh_location_tolerance_set(PDM_mesh_location_t *ml, double tol)
  # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

  # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  void PDM_mesh_location_method_set(PDM_mesh_location_t *ml, PDM_mesh_location_method_t method)
  # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

  # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  void PDM_mesh_location_compute(PDM_mesh_location_t *ml)
  # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  int PDM_mesh_location_n_located_get (PDM_mesh_location_t *ml,
                                       int                  i_point_cloud,
                                       int                  i_part)

  int PDM_mesh_location_n_unlocated_get (PDM_mesh_location_t *ml,
                                         int                  i_point_cloud,
                                         int                  i_part)

  int *PDM_mesh_location_unlocated_get (PDM_mesh_location_t *ml,
                                        int                  i_point_cloud,
                                        int                  i_part)


  int *PDM_mesh_location_located_get (PDM_mesh_location_t *ml,
                                      int                  i_point_cloud,
                                      int                  i_part)

  # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  void PDM_mesh_location_point_location_get(PDM_mesh_location_t  *ml,
                                            int                   i_point_cloud,
                                            int                   i_part,
                                            PDM_g_num_t         **location,
                                            double              **dist2,
                                            double              **projected_coord)
  # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  void PDM_mesh_location_points_in_elt_get (PDM_mesh_location_t  *ml,
                                            int                   i_part,
                                            int                   i_point_cloud,
                                            int                 **elt_pts_inside_idx,
                                            PDM_g_num_t         **points_gnum,
                                            double              **points_coords,
                                            double              **points_uvw,
                                            int                 **points_weights_idx,
                                            double              **points_weights,
                                            double              **points_dist2,
                                            double              **points_projected_coords)

  # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  void PDM_mesh_location_free(PDM_mesh_location_t  *ml)
  # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

  # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  void PDM_mesh_location_dump_times(PDM_mesh_location_t  *ml)
  # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

  # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  void PDM_mesh_location_mesh_nodal_id_get(PDM_mesh_location_t  *ml)
  # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

  # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  int PDM_mesh_location_n_cell_get (PDM_mesh_location_t  *ml, int i_part)
  # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

  # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  void PDM_mesh_location_cell_vertex_get (PDM_mesh_location_t  *ml, int i_part, int **cell_vtx_idx, int **cell_vtx)
  # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

  # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  void PDM_mesh_location_part_to_part_get(PDM_mesh_location_t  *ml,
                                          const int             icloud,
                                          PDM_part_to_part_t  **ptp,
                                          PDM_ownership_t       ownership)
  # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

# ------------------------------------------------------------------
cdef class MeshLocation:
  # ************************************************************************
  # > Class attributes
  cdef PDM_mesh_location_t* _ml
  cdef MPI.Comm py_comm
  cdef int _n_point_cloud
  cdef int _n_src_part
  cdef list _n_tgt_part_per_cloud
  
  cdef double _tolerance
  cdef PDM_mesh_location_method_t _method

  cdef list _np_located, _np_unlocated, _dic_location, _dic_points_in_elt, _dic_cell_vertex

  #PDM_mesh_location_method_t
  OCTREE = _PDM_MESH_LOCATION_OCTREE
  DBBTREE = _PDM_MESH_LOCATION_DBBTREE
  LOCATE_ALL_TGT = _PDM_MESH_LOCATION_LOCATE_ALL_TGT

  # ************************************************************************

  # ------------------------------------------------------------------------
  def __init__(self, int               n_point_cloud,
                     MPI.Comm          comm):
    """
    __init__(n_point_cloud, comm)

    Create the structure

    Parameters:
      n_point_cloud (int)               : Number of point clouds
      comm          (MPI.Comm)          : MPI communicator
    """

    # ::::::::::::::::::::::::::::::::::::::::::::::::::
    self.py_comm = comm
    self._n_point_cloud = n_point_cloud
    self._n_tgt_part_per_cloud = [0 for i in range(n_point_cloud)]
    # ::::::::::::::::::::::::::::::::::::::::::::::::::

    # ::::::::::::::::::::::::::::::::::::::::::::::::::
    # > Convert mpi4py -> PDM_MPI
    cdef MPI.MPI_Comm c_comm = comm.ob_mpi
    cdef PDM_MPI_Comm PDMC   = PDM_MPI_mpi_2_pdm_mpi_comm(<void *> &c_comm)
    # ::::::::::::::::::::::::::::::::::::::::::::::::::

    # ::::::::::::::::::::::::::::::::::::::::::::::::::
    self._ml = PDM_mesh_location_create(n_point_cloud, PDMC, PDM_OWNERSHIP_UNGET_RESULT_IS_FREE)

    # ::::::::::::::::::::::::::::::::::::::::::::::::::
    # Set defaults to ensure consistency with doc
    self.tolerance = 0
    self.method = MeshLocation.OCTREE

  # ------------------------------------------------------------------------
  def n_part_cloud_set(self, int i_point_cloud,
                             int n_part):
    """
    n_part_cloud_set(i_point_cloud, n_part)

    Set the number of partitions of the specified point cloud

    Parameters:
      i_point_cloud (int) : Point cloud identifier
      n_part        (int) : Number of partitions for this point cloud
    """
    self._n_tgt_part_per_cloud[i_point_cloud] = n_part
    PDM_mesh_location_n_part_cloud_set(self._ml,
                                       i_point_cloud,
                                       n_part)

  # ------------------------------------------------------------------------
  def cloud_set(self, int i_point_cloud,
                      int i_part,
                      NPY.ndarray[NPY.double_t  , mode='c', ndim=1] coords,
                      NPY.ndarray[npy_pdm_gnum_t, mode='c', ndim=1] gnum,
                      ):
    """
    cloud_set(i_point_cloud, i_part, coords, gnum)

    Set a partition of the specified point cloud

    Parameters:
      i_point_cloud (int)                        : Point cloud identifier
      i_part        (int)                        : Partition identifier
      coords        (np.ndarray[np.double_t])    : Point coordinates
      gnum          (np.ndarray[npy_pdm_gnum_t]) : Point global ids
    """
    cdef int n_points = len(gnum)
    assert coords.size == 3*n_points
    PDM_mesh_location_cloud_set(self._ml,
                                i_point_cloud,
                                i_part,
                                n_points,
                 <double*>      coords.data,
                 <PDM_g_num_t*> gnum.data)


  # ------------------------------------------------------------------------
  def mesh_n_part_set(self, int n_part):
    """
    mesh_n_part_set(n_part)

    Set the number of partitions of the mesh

    Parameters:
      n_part (int) : Number of partitions
    """
    self._n_src_part = n_part
    PDM_mesh_location_mesh_n_part_set(self._ml, n_part)

  # ------------------------------------------------------------------------
  def part_set(self, int i_part,
                     NPY.ndarray[NPY.int32_t   , mode='c', ndim=1] cell_face_idx,
                     NPY.ndarray[NPY.int32_t   , mode='c', ndim=1] cell_face,
                     NPY.ndarray[npy_pdm_gnum_t, mode='c', ndim=1] cell_ln_to_gn,
                     NPY.ndarray[NPY.int32_t   , mode='c', ndim=1] face_vtx_idx,
                     NPY.ndarray[NPY.int32_t   , mode='c', ndim=1] face_vtx,
                     NPY.ndarray[npy_pdm_gnum_t, mode='c', ndim=1] face_ln_to_gn,
                     NPY.ndarray[NPY.double_t  , mode='c', ndim=1] coords,
                     NPY.ndarray[npy_pdm_gnum_t, mode='c', ndim=1] vtx_ln_to_gn):
    """
    part_set(i_part, cell_face_idx, cell_face, cell_ln_to_gn, face_vtx_idx, face_vtx, face_ln_to_gn, coords, vtx_ln_to_gn)

    Set a *volume* mesh partition

    Parameters:
      i_part        (int)                        : Partition identifier
      cell_face_idx (np.ndarray[np.int32_t])     : Index for cell -> face connectivity
      cell_face     (np.ndarray[np.int32_t])     : Cell -> face connectivity
      cell_ln_to_gn (np.ndarray[npy_pdm_gnum_t]) : Cell global ids
      face_vtx_idx  (np.ndarray[np.int32_t])     : Index for face -> vertex connectivity
      face_vtx      (np.ndarray[np.int32_t])     : Face -> vertex connectivity
      face_ln_to_gn (np.ndarray[npy_pdm_gnum_t]) : Face global ids
      coords        (np.ndarray[np.double_t])    : Vertex coordinates
      vtx_ln_to_gn  (np.ndarray[npy_pdm_gnum_t]) : Vertex global ids
    """
    cdef int n_cell = len(cell_face_idx) - 1
    cdef int n_face = len(face_vtx_idx)  - 1
    cdef int n_vtx  = len(vtx_ln_to_gn)
    
    assert cell_ln_to_gn.size == n_cell
    assert face_ln_to_gn.size == n_face
    assert coords.size == 3*n_vtx
    
    PDM_mesh_location_part_set(self._ml,
                               i_part,
                               n_cell,
                <int*>         cell_face_idx.data,
                <int*>         cell_face.data,
                <PDM_g_num_t*> cell_ln_to_gn.data,
                               n_face,
                <int*>         face_vtx_idx.data,
                <int*>         face_vtx.data,
                <PDM_g_num_t*> face_ln_to_gn.data,
                               n_vtx,
                <double*>      coords.data,
                <PDM_g_num_t*> vtx_ln_to_gn.data)

  # ------------------------------------------------------------------------
  def nodal_part_set(self, int i_part,
                     NPY.ndarray[NPY.int32_t   , mode='c', ndim=1] cell_vtx_idx,
                     NPY.ndarray[NPY.int32_t   , mode='c', ndim=1] cell_vtx,
                     NPY.ndarray[npy_pdm_gnum_t, mode='c', ndim=1] cell_ln_to_gn,
                     NPY.ndarray[NPY.double_t  , mode='c', ndim=1] coords,
                     NPY.ndarray[npy_pdm_gnum_t, mode='c', ndim=1] vtx_ln_to_gn):
    """
    nodal_part_set(i_part, cell_vtx_idx, cell_vtx, cell_ln_to_gn, coords, vtx_ln_to_gn)

    Set a *volume* mesh partition defined by nodal connectivity

    The mesh is assumed to contain only standard elements
    (tetrahedra, pyramids, prisms, hexahedra).

    Parameters:
      i_part        (int)                        : Partition identifier
      cell_vtx_idx  (np.ndarray[np.int32_t])     : Index for cell -> vertex connectivity
      cell_vtx      (np.ndarray[np.int32_t])     : Cell -> vertex connectivity
      cell_ln_to_gn (np.ndarray[npy_pdm_gnum_t]) : Cell global ids
      coords        (np.ndarray[np.double_t])    : Vertex coordinates
      vtx_ln_to_gn  (np.ndarray[npy_pdm_gnum_t]) : Vertex global ids
    """
    cdef int n_cell = len(cell_vtx_idx) - 1
    cdef int n_vtx  = len(vtx_ln_to_gn)

    assert cell_ln_to_gn.size == n_cell 
    assert coords.size == 3*n_vtx

    PDM_mesh_location_nodal_part_set(self._ml,
                                     i_part,
                                     n_cell,
                      <int*>         cell_vtx_idx.data,
                      <int*>         cell_vtx.data,
                      <PDM_g_num_t*> cell_ln_to_gn.data,
                                     n_vtx,
                      <double*>      coords.data,
                      <PDM_g_num_t*> vtx_ln_to_gn.data)

  # ------------------------------------------------------------------------
  def part_set_2d(self, int i_part,
                  NPY.ndarray[NPY.int32_t   , mode='c', ndim=1] face_edge_idx,
                  NPY.ndarray[NPY.int32_t   , mode='c', ndim=1] face_edge,
                  NPY.ndarray[npy_pdm_gnum_t, mode='c', ndim=1] face_ln_to_gn,
                  NPY.ndarray[NPY.int32_t   , mode='c', ndim=1] edge_vtx,
                  NPY.ndarray[NPY.double_t  , mode='c', ndim=1] coords,
                  NPY.ndarray[npy_pdm_gnum_t, mode='c', ndim=1] vtx_ln_to_gn):
    """
    part_set_2d(i_part, face_edge_idx, face_edge, face_ln_to_gn, edge_vtx, coords, vtx_ln_to_gn)

    Set a *surface* mesh partition

    Parameters:
      i_part        (int)                        : Partition identifier
      face_edge_idx (np.ndarray[np.int32_t])     : Index for face -> edge connectivity
      face_edge     (np.ndarray[np.int32_t])     : Face -> edge connectivity
      face_ln_to_gn (np.ndarray[npy_pdm_gnum_t]) : Face global ids
      edge_vtx      (np.ndarray[np.int32_t])     : Edge -> vertex connectivity
      coords        (np.ndarray[np.double_t])    : Vertex coordinates
      vtx_ln_to_gn  (np.ndarray[npy_pdm_gnum_t]) : Vertex global ids
    """
    cdef int n_vtx  = len(vtx_ln_to_gn)
    cdef int n_edge = len(edge_vtx)//2
    cdef int n_face = len(face_edge_idx) - 1

    assert face_ln_to_gn.size == n_face
    assert coords.size == 3*n_vtx

    PDM_mesh_location_part_set_2d(self._ml,
                                  i_part,
                                  n_face,
                   <int*>         face_edge_idx.data,
                   <int*>         face_edge.data,
                   <PDM_g_num_t*> face_ln_to_gn.data,
                                  n_edge,
                   <int*>         edge_vtx.data,
                                  n_vtx,
                   <double*>      coords.data,
                   <PDM_g_num_t*> vtx_ln_to_gn.data)

  # ------------------------------------------------------------------------
  def nodal_part_set_2d(self, int i_part,
                        NPY.ndarray[NPY.int32_t   , mode='c', ndim=1] face_vtx_idx,
                        NPY.ndarray[NPY.int32_t   , mode='c', ndim=1] face_vtx,
                        NPY.ndarray[npy_pdm_gnum_t, mode='c', ndim=1] face_ln_to_gn,
                        NPY.ndarray[NPY.double_t  , mode='c', ndim=1] coords,
                        NPY.ndarray[npy_pdm_gnum_t, mode='c', ndim=1] vtx_ln_to_gn):
    """
    nodal_part_set_2d(i_part, face_vtx_idx, face_vtx, face_ln_to_gn, coords, vtx_ln_to_gn)

    Set a *surface* mesh partition with nodal connectivity

    Parameters:
      i_part        (int)                        : Partition identifier
      face_vtx_idx  (np.ndarray[np.int32_t])     : Index for face -> vertex connectivity
      face_vtx      (np.ndarray[np.int32_t])     : Face -> vertex connectivity
      face_ln_to_gn (np.ndarray[npy_pdm_gnum_t]) : Face global ids
      coords        (np.ndarray[np.double_t])    : Vertex coordinates
      vtx_ln_to_gn  (np.ndarray[npy_pdm_gnum_t]) : Vertex global ids
    """
    cdef int n_vtx  = len(vtx_ln_to_gn)
    cdef int n_face = len(face_vtx_idx) - 1

    assert face_ln_to_gn.size == n_face
    assert coords.size == 3*n_vtx

    PDM_mesh_location_nodal_part_set_2d(self._ml,
                                        i_part,
                                        n_face,
                         <int*>         face_vtx_idx.data,
                         <int*>         face_vtx.data,
                         <PDM_g_num_t*> face_ln_to_gn.data,
                                        n_vtx,
                         <double*>      coords.data,
                         <PDM_g_num_t*> vtx_ln_to_gn.data)

  # ------------------------------------------------------------------------
  @property
  def tolerance(self):
    """
    Relative tolerance for bounding boxes. Default value is 0.
    """
    return self._tolerance

  @tolerance.setter
  def tolerance(self, double tolerance):
    """ Tolerance setter """
    PDM_mesh_location_tolerance_set(self._ml, tolerance)
    self._tolerance = tolerance

  # ------------------------------------------------------------------------
  @property
  def method(self):
    """
    Method used for the preconditioning stage of computing location

    Admissible values are : 
      - :py:attr:`MeshLocation.OCTREE` : Use point octree (default method)
      - :py:attr:`MeshLocation.DBBTREE` : Use bounding-box tree
      - :py:attr:`MeshLocation.LOCATE_ALL_TGT` : All target points are guaranteed to be located
    """
    return self._method

  @method.setter
  def method(self, PDM_mesh_location_method_t method):
    """ Method setter """
    PDM_mesh_location_method_set(self._ml, method)
    self._method = method

  # ------------------------------------------------------------------------
  def __location_get(self, int i_point_cloud, int i_part):
    """
    """
    # ************************************************************************
    # > Declaration
    cdef int           n_points
    cdef int           n_located
    cdef double       *coords
    cdef PDM_g_num_t  *gnum
    cdef PDM_g_num_t  *location
    cdef double       *dist2
    cdef double       *p_proj_coord
    # ************************************************************************

    # Attention : Mettre les fonction n_located_ et located !!!!
    PDM_mesh_location_cloud_get (self._ml,
                                 i_point_cloud,
                                 i_part,
                                 &n_points,
                                 &coords,
                                 &gnum)

    n_located = PDM_mesh_location_n_located_get(self._ml, i_point_cloud, i_part)

    PDM_mesh_location_point_location_get(self._ml,
                                         i_point_cloud,
                                         i_part,
                                         &location,
                                         &dist2,
                                         &p_proj_coord)

    return {
            'g_num'        : create_numpy_g (gnum,         n_points, flag_owndata=False),
            'location'     : create_numpy_g (location,     n_located),
            'dist2'        : create_numpy_d (dist2,        n_located),
            'p_proj_coord' : create_numpy_d (p_proj_coord, 3*n_located)
           }

  def location_get(self, int i_point_cloud, int i_part):
    """
    location_get(i_point_cloud, i_part)

    Get location data on the target side for the specified point cloud and partition.

    .. note::
      The results are related to located points only, whose ids can be accessed with function :py:func:`located_get`.
    

    Parameters:
      i_point_cloud (int) : Point cloud identifier
      i_part        (int) : Partition identifier

    Returns:
      Dictionary
        - ``"g_num"``        (`np.ndarray[npy_pdm_gnum_t]`) : Point global ids
        - ``"location"``     (`np.ndarray[npy_pdm_gnum_t]`) : Global id of nearest mesh element
        - ``"dist2"``        (`np.ndarray[np.double_t]`)    : Signed squared distance from nearest element (negative if the point is located inside that element)
        - ``"p_proj_coord"`` (`np.ndarray[np.double_t]`)    : Cartesian coordinates of projection onto the nearest element (identity if the point is located inside that element)
    """
    return self._dic_location[i_point_cloud][i_part]


  def __cell_vertex_get (self, int i_part):

    cdef int *cell_vtx_idx
    cdef int *cell_vtx

    cdef int n_elts = PDM_mesh_location_n_cell_get(self._ml, i_part)
    PDM_mesh_location_cell_vertex_get(self._ml, i_part, &cell_vtx_idx, &cell_vtx)

    return {'cell_vtx_idx'      : create_numpy_i(cell_vtx_idx, n_elts+1),
            'cell_vtx'          : create_numpy_i(cell_vtx, cell_vtx_idx[n_elts])
           }
  def cell_vertex_get (self, int i_part):
    """
    cell_vertex_get(i_part)

    Get the cell->vertex connectivity used for internal computations

    .. note::
      For non-standard elements, this connectivity is built by ParaDiGM and is necessary to associate the ``points_weights`` array (returned by :py:meth:`points_in_elt_get`) to the appropriate mesh vertices.

    Parameters:
      i_part (int) : Partition identifier

    Returns:
      Dictionary
        - ``"cell_vtx_idx"`` (`np.ndarray[np.int32_t]`) : Index for cell -> vertex connectivity
        - ``"cell_vtx"``     (`np.ndarray[np.int32_t]`) : Cell -> vertex connectivity
    """
    return self._dic_cell_vertex[i_part]



  def __points_in_elt_get(self, int i_point_cloud, int i_part):
    """
    """
    # ************************************************************************
    # > Declaration
    cdef int          *elt_pts_inside_idx
    cdef PDM_g_num_t  *points_gnum
    cdef double       *points_coords
    cdef double       *points_uvw
    cdef int          *points_weights_idx
    cdef double       *points_weights
    cdef double       *points_dist2
    cdef double       *points_projected_coords
    # ************************************************************************

    PDM_mesh_location_points_in_elt_get(self._ml,
                                        i_point_cloud,
                                        i_part,
                                        &elt_pts_inside_idx,
                                        &points_gnum,
                                        &points_coords,
                                        &points_uvw,
                                        &points_weights_idx,
                                        &points_weights,
                                        &points_dist2,
                                        &points_projected_coords)

    cdef int n_elts =  PDM_mesh_location_n_cell_get(self._ml, i_part)
    cdef int s_loc  = elt_pts_inside_idx[n_elts]
    cdef int s_wei  = points_weights_idx[s_loc]

    return {'elt_pts_inside_idx'      : create_numpy_i (elt_pts_inside_idx,      n_elts+1),
            'points_gnum'             : create_numpy_g (points_gnum,             s_loc   ),
            'points_coords'           : create_numpy_d (points_coords,           3*s_loc ),
            'points_uvw'              : create_numpy_d (points_uvw,              3*s_loc ),
            'points_weights_idx'      : create_numpy_i (points_weights_idx,      s_loc+1 ),
            'points_weights'          : create_numpy_d (points_weights,          s_wei   ),
            'points_dist2'            : create_numpy_d (points_dist2,            s_loc   ),
            'points_projected_coords' : create_numpy_d (points_projected_coords, 3*s_loc )
            }

  def points_in_elt_get(self, int i_point_cloud, int i_part):
    """
    points_in_elt_get(i_point_cloud, i_part)

    Get location data for the mesh cells of the specified partition identifier.

    Parameters:
      i_point_cloud (int) : Point cloud identifier
      i_part        (int) : Partition identifier

    Returns:
      Dictionary
        - ``"elt_pts_inside_idx"``      (`np.ndarray[np.int32_t]`)     : Index for element -> points mapping
        - ``"points_gnum"``             (`np.ndarray[npy_pdm_gnum_t]`) : Located points global ids
        - ``"points_coords"``           (`np.ndarray[np.double_t]`)    : Located points cartesian coordinates
        - ``"points_uvw"``              (`np.ndarray[np.double_t]`)    : Located points parametric coordinates
        - ``"points_weights_idx"``      (`np.ndarray[np.int32_t]`)     : Index for interpolation weights
        - ``"points_weights"``          (`np.ndarray[np.double_t]`)    : Interpolation weights
        - ``"points_dist2"``            (`np.ndarray[np.double_t]`)    : Signed squared distance element-points
        - ``"points_projected_coords"`` (`np.ndarray[np.double_t]`)    : Cartesian coordinates of projection on element
    """
    return self._dic_points_in_elt[i_point_cloud][i_part]

  def __located_get(self, int i_point_cloud, int i_part):
    """
    """
    cdef int n_located = PDM_mesh_location_n_located_get(self._ml, i_point_cloud, i_part)
    cdef int* located = PDM_mesh_location_located_get(self._ml, i_point_cloud, i_part)

    return create_numpy_i(located, n_located)

  def located_get(self, int i_point_cloud, int i_part):
    """
    located_get(i_point_cloud, i_part)

    Get the list of located points

    Parameters:
      i_point_cloud (int) : Point cloud identifier
      i_part        (int) : Partition identifier

    Returns:
      List of located target points (1-based ids) (`np.ndarray[np.int32_t]`)
    """
    return self._np_located[i_point_cloud][i_part]

  def __unlocated_get(self, int i_point_cloud, int i_part):
    """
    """
    cdef int n_unlocated = PDM_mesh_location_n_unlocated_get(self._ml, i_point_cloud, i_part)
    cdef int* unlocated = PDM_mesh_location_unlocated_get(self._ml, i_point_cloud, i_part)

    return create_numpy_i(unlocated, n_unlocated)

  def unlocated_get(self, int i_point_cloud, int i_part):
    """
    unlocated_get(i_point_cloud, i_part)

    Get the list of unlocated points

    Parameters:
      i_point_cloud (int) : Point cloud identifier
      i_part        (int) : Partition identifier

    Returns:
      List of unlocated target points (1-based ids) (`np.ndarray[np.int32_t]`)
    """
    return self._np_unlocated[i_point_cloud][i_part]

  # ------------------------------------------------------------------------
  def compute(self):
    """
    Compute point location
    """
    # ************************************************************************
    # > Declaration
    # ************************************************************************
    PDM_mesh_location_compute(self._ml)

    #Take ownership of all computed data to avoid memory leaks
    self._np_unlocated      = []
    self._np_located        = []
    self._dic_location      = []
    self._dic_points_in_elt = []
    for i_pt_cloud in range(self._n_point_cloud):
      #Target related data
      n_tgt_part = self._n_tgt_part_per_cloud[i_pt_cloud]
      self._np_unlocated.append([self.__unlocated_get(i_pt_cloud, i_part) for i_part in range(n_tgt_part)])
      self._np_located  .append([self.__located_get  (i_pt_cloud, i_part) for i_part in range(n_tgt_part)])
      self._dic_location.append([self.__location_get (i_pt_cloud, i_part) for i_part in range(n_tgt_part)])
      #Source related data
      self._dic_points_in_elt.append([self.__points_in_elt_get(i_pt_cloud, i_part) for i_part in range(self._n_src_part)])
    #Source related data
    self._dic_cell_vertex = [self.__cell_vertex_get(i_part) for i_part in range(self._n_src_part)]

  # ------------------------------------------------------------------------
  def part_to_part_get(self, int i_point_cloud):
    """
    part_to_part_get(i_point_cloud)

    Get the PartToPart object to exchange data between
    the source mesh and a target point cloud

    Parameters:
      i_point_cloud (int) : Point cloud identifier

    Returns:
      PartToPart object (:py:class:`PartToPart`)
    """
    cdef PDM_part_to_part_t *ptpc
    PDM_mesh_location_part_to_part_get(self._ml,
                                       i_point_cloud,
                                       &ptpc,
                                       PDM_OWNERSHIP_USER)

    py_caps = PyCapsule_New(ptpc, NULL, NULL)
    return PartToPartCapsule(py_caps, self.py_comm) # The free is inside the class

  # ------------------------------------------------------------------------
  def dump_times(self):
    """
    Dump elapsed and CPU times
    """
    # ************************************************************************
    # > Declaration
    # ************************************************************************
    PDM_mesh_location_dump_times(self._ml)

  # ------------------------------------------------------------------------
  def __dealloc__(self):
    """
    Free a mesh location structure
    """
    # ************************************************************************
    # > Declaration
    # ************************************************************************
    PDM_mesh_location_free(self._ml)
