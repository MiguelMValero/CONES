cdef extern from "pdm_generate_mesh.h":
  void PDM_generate_mesh_rectangle_simplified(PDM_MPI_Comm   comm,
                                              PDM_g_num_t    n_vtx_seg,
                                              int           *n_vtx,
                                              int           *n_elt,
                                              double       **coords,
                                              int          **elt_vtx_idx,
                                              int          **elt_vtx);

  void PDM_generate_mesh_rectangle_ngon(const PDM_MPI_Comm            comm,
                                        const PDM_Mesh_nodal_elt_t    elt_type,
                                        const double                  xmin,
                                        const double                  ymin,
                                        const double                  zmin,
                                        const double                  lengthx,
                                        const double                  lengthy,
                                        const PDM_g_num_t             n_x,
                                        const PDM_g_num_t             n_y,
                                        const int                     n_part,
                                        const PDM_split_dual_t        part_method,
                                        const double                  random_factor,
                                        int                         **pn_vtx,
                                        int                         **pn_edge,
                                        int                         **pn_face,
                                        double                     ***pvtx_coord,
                                        int                        ***pedge_vtx,
                                        int                        ***pface_edge_idx,
                                        int                        ***pface_edge,
                                        int                        ***pface_vtx,
                                        PDM_g_num_t                ***pvtx_ln_to_gn,
                                        PDM_g_num_t                ***pedge_ln_to_gn,
                                        PDM_g_num_t                ***pface_ln_to_gn);

# ------------------------------------------------------------------

cdef extern from "pdm_mpi.h":
    PDM_MPI_Comm PDM_MPI_mpi_2_pdm_mpi_comm (void *mpi_comm)

# ------------------------------------------------------------------

def generate_mesh_rectangle_simplified(MPI.Comm      comm,
                                       PDM_g_num_t   n_vtx_seg):

  """
  generate_mesh_rectangle_simplified(comm, n_vtx_seg)

  Create a simple partitioned rectangle mesh (2D).

  Parameters:
    comm      (MPI.Comm) : MPI communicator
    n_vtx_seg (int)      : Number of vertices along each side of the rectangle

  Returns:
    Dictionary
      - ``"n_vtx"``       (`int`)                     : Number of vertices
      - ``"n_elt"``       (`int`)                     : Number of elements
      - ``"coords"``      (`np.ndarray[np.double_t]`) : Coordinates
      - ``"elt_vtx_idx"`` (`np.ndarray[np.int32_t]`)  : Index of element-vertex connectivity
      - ``"elt_vtx"``     (`np.ndarray[np.int32_t]`)  : Element-vertex connectivity
  """

  # MPI communicator
  cdef MPI.MPI_Comm c_comm = comm.ob_mpi

  # Get
  cdef int n_vtx = -1
  cdef int n_elt = -1
  cdef double *coords
  cdef int    *elt_vtx_idx
  cdef int    *elt_vtx

  PDM_generate_mesh_rectangle_simplified(PDM_MPI_mpi_2_pdm_mpi_comm (<void *> &c_comm),
                                         n_vtx_seg,
                                         &n_vtx,
                                         &n_elt,
                                         &coords,
                                         &elt_vtx_idx,
                                         &elt_vtx)

  return {
          'n_vtx'       : n_vtx,
          'n_elt'       : n_elt,
          'coords'      : create_numpy_d(coords, 3*n_vtx),
          'elt_vtx_idx' : create_numpy_i(elt_vtx_idx, n_elt+1),
          'elt_vtx'     : create_numpy_i(elt_vtx, elt_vtx_idx[n_elt])
         }

# ------------------------------------------------------------------------
def generate_mesh_rectangle_ngon(MPI.Comm             comm,
                                 PDM_Mesh_nodal_elt_t elt_type,
                                 double               xmin,
                                 double               ymin,
                                 double               zmin,
                                 double               lengthx,
                                 double               lengthy,
                                 PDM_g_num_t          n_x,
                                 PDM_g_num_t          n_y,
                                 int                  n_part,
                                 PDM_split_dual_t     part_method,
                                 double               random_factor=0):
  """
  generate_mesh_rectangle_ngon(comm, elt_type, xmin, ymin, zmin, lengthx, lengthy, n_x, n_y, n_part, part_method, random_factor=0)

  Create a partitioned rectangular mesh (2D) with descending connectivities

  Parameters:
    comm          (MPI.Comm)         : MPI communicator
    elt_type      (int)              : Element type
    xmin          (double)           : Minimal x-coordinate
    ymin          (double)           : Minimal y-coordinate
    zmin          (double)           : Minimal z-coordinate
    lengthx       (double)           : Length of the rectangle in the x-direction
    lengthy       (double)           : Length of the rectangle in the y-direction
    n_x           (int)              : Number of points in the x-direction
    n_y           (int)              : Number of points in the y-direction
    n_part        (int)              : Number of partitions
    part_method   (int)              : Partitioning method
    random_factor (double, optional) : Randomization factor (between 0 and 1, default = 0)

  Returns:
    Dictionary
      - ``"pn_vtx"``         (`list` of `int`)                        : Number of vertices
      - ``"pn_edge"``        (`list` of `int`)                        : Number of edges
      - ``"pn_face"``        (`list` of `int`)                        : Number of faces
      - ``"pvtx_coord"``     (`list` of `np.ndarray[np.double_t]`)    : Vertex coordinates
      - ``"pedge_vtx"``      (`list` of `np.ndarray[np.int32_t]`)     : Edge->vertex connectivity
      - ``"pface_edge_idx"`` (`list` of `np.ndarray[np.int32_t]`)     : Index of face->edge connectivity
      - ``"pface_edge"``     (`list` of `np.ndarray[np.int32_t]`)     : Face->edge connectivity
      - ``"pface_vtx"``      (`list` of `np.ndarray[np.int32_t]`)     : Face->vertex connectivity
      - ``"pvtx_ln_to_gn"``  (`list` of `np.ndarray[npy_pdm_gnum_t]`) : Vertex global ids
      - ``"pedge_ln_to_gn"`` (`list` of `np.ndarray[npy_pdm_gnum_t]`) : Edge global ids
      - ``"pface_ln_to_gn"`` (`list` of `np.ndarray[npy_pdm_gnum_t]`) : Face global ids

  Admissible values for ``elt_type``:
    - 2 : triangles
    - 3 : quadrangles
    - 4 : mixed polygons (triangles, quadrangles and octagons)
  """

  # MPI communicator
  cdef MPI.MPI_Comm c_comm = comm.ob_mpi

  # Get
  cdef int          *_pn_vtx         = NULL
  cdef int          *_pn_edge        = NULL
  cdef int          *_pn_face        = NULL
  cdef double      **_pvtx_coord     = NULL
  cdef int         **_pedge_vtx      = NULL
  cdef int         **_pface_edge_idx = NULL
  cdef int         **_pface_edge     = NULL
  cdef int         **_pface_vtx      = NULL
  cdef PDM_g_num_t **_pvtx_ln_to_gn  = NULL
  cdef PDM_g_num_t **_pedge_ln_to_gn = NULL
  cdef PDM_g_num_t **_pface_ln_to_gn = NULL

  PDM_generate_mesh_rectangle_ngon(PDM_MPI_mpi_2_pdm_mpi_comm (<void *> &c_comm),
                                   elt_type,
                                   xmin,
                                   ymin,
                                   zmin,
                                   lengthx,
                                   lengthy,
                                   n_x,
                                   n_y,
                                   n_part,
                                   part_method,
                                   random_factor,
                                   &_pn_vtx,
                                   &_pn_edge,
                                   &_pn_face,
                                   &_pvtx_coord,
                                   &_pedge_vtx,
                                   &_pface_edge_idx,
                                   &_pface_edge,
                                   &_pface_vtx,
                                   &_pvtx_ln_to_gn,
                                   &_pedge_ln_to_gn,
                                   &_pface_ln_to_gn)

  pn_vtx         = []
  pn_edge        = []
  pn_face        = []
  pvtx_coord     = []
  pedge_vtx      = []
  pface_edge_idx = []
  pface_edge     = []
  pface_vtx      = []
  pvtx_ln_to_gn  = []
  pedge_ln_to_gn = []
  pface_ln_to_gn = []
  for i_part in range(n_part):
    pn_vtx .append(_pn_vtx [i_part])
    pn_edge.append(_pn_edge[i_part])
    pn_face.append(_pn_face[i_part])

    pvtx_coord.append(create_numpy_d(_pvtx_coord[i_part], 3*_pn_vtx[i_part]))

    pedge_vtx.append(create_numpy_i(_pedge_vtx[i_part], 2*_pn_edge[i_part]))

    pface_edge_idx.append(create_numpy_i(_pface_edge_idx[i_part], _pn_face[i_part]+1))
    pface_edge.append(create_numpy_i(_pface_edge[i_part], _pface_edge_idx[i_part][_pn_face[i_part]]))
    pface_vtx.append (create_numpy_i(_pface_vtx[i_part],  _pface_edge_idx[i_part][_pn_face[i_part]]))

    pvtx_ln_to_gn .append(create_numpy_g(_pvtx_ln_to_gn [i_part], _pn_vtx [i_part]))
    pedge_ln_to_gn.append(create_numpy_g(_pedge_ln_to_gn[i_part], _pn_edge[i_part]))
    pface_ln_to_gn.append(create_numpy_g(_pface_ln_to_gn[i_part], _pn_face[i_part]))
  free(_pn_vtx        )
  free(_pn_edge       )
  free(_pn_face       )
  free(_pvtx_coord    )
  free(_pedge_vtx     )
  free(_pface_edge_idx)
  free(_pface_edge    )
  free(_pface_vtx     )
  free(_pvtx_ln_to_gn )
  free(_pedge_ln_to_gn)
  free(_pface_ln_to_gn)

  return {
    "pn_vtx"         : pn_vtx,
    "pn_edge"        : pn_edge,
    "pn_face"        : pn_face,
    "pvtx_coord"     : pvtx_coord,
    "pedge_vtx"      : pedge_vtx,
    "pface_edge_idx" : pface_edge_idx,
    "pface_edge"     : pface_edge,
    "pface_vtx"      : pface_vtx,
    "pvtx_ln_to_gn"  : pvtx_ln_to_gn,
    "pedge_ln_to_gn" : pedge_ln_to_gn,
    "pface_ln_to_gn" : pface_ln_to_gn
  }

