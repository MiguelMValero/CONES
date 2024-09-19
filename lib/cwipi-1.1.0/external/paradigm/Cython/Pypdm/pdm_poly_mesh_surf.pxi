
cdef extern from "pdm_poly_surf_gen.h":
    # ------------------------------------------------------------------
    void PDM_poly_surf_gen(PDM_MPI_Comm  localComm,
                           double        xmin,
                           double        xmax,
                           double        ymin,
                           double        ymax,
                           int           have_random,
                           int           init_random,
                           PDM_g_num_t   nx,
                           PDM_g_num_t   ny,
                           PDM_g_num_t  *n_face,
                           PDM_g_num_t  *n_vtx,
                           PDM_g_num_t  *n_edge,
                           int          *dn_vtx,
                           double      **dvtx_coord,
                           int          *dn_face,
                           int         **dface_vtx_idx,
                           PDM_g_num_t **dface_vtx,
                           PDM_g_num_t **dface_edge,
                           int          *dn_edge,
                           PDM_g_num_t **dedge_vtx,
                           PDM_g_num_t **dedge_face,
                           int          *n_edge_group,
                           int         **dedge_group_idx,
                           PDM_g_num_t **dedge_group)
    # ------------------------------------------------------------------

# ------------------------------------------------------------------------
def PolyMeshSurf(double      xmin,
                 double      xmax,
                 double      ymin,
                 double      ymax,
                 int         have_random,
                 int         init_random,
                 PDM_g_num_t nx,
                 PDM_g_num_t ny,
                 MPI.Comm    comm):
  """
  """
  # ************************************************************************
  # > Declaration
  cdef PDM_g_num_t  n_face
  cdef PDM_g_num_t  n_vtx
  cdef PDM_g_num_t  n_edge
  cdef int          dn_vtx
  cdef double      *dvtx_coord
  cdef int          dn_face
  cdef int         *dface_vtx_idx
  cdef PDM_g_num_t *dface_vtx
  cdef PDM_g_num_t *dface_edge
  cdef int          dn_edge
  cdef PDM_g_num_t *dedge_vtx
  cdef PDM_g_num_t *dedge_face
  cdef int          n_edge_group
  cdef int         *edge_group
  cdef int         *dedge_group_idx
  cdef PDM_g_num_t *dedge_group
  # ************************************************************************

  cdef MPI.MPI_Comm c_comm = comm.ob_mpi
  cdef PDM_MPI_Comm PDMC   = PDM_MPI_mpi_2_pdm_mpi_comm(<void *> &c_comm)

  PDM_poly_surf_gen(PDMC,
                    xmin,
                    xmax,
                    ymin,
                    ymax,
                    have_random,
                    init_random,
                    nx,
                    ny,
                    &n_face,
                    &n_vtx,
                    &n_edge,
                    &dn_vtx,
                    &dvtx_coord,
                    &dn_face,
                    &dface_vtx_idx,
                    &dface_vtx,
                    &dface_edge,
                    &dn_edge,
                    &dedge_vtx,
                    &dedge_face,
                    &n_edge_group,
                    &dedge_group_idx,
                    &dedge_group)

  np_dvtx_coord      = create_numpy_d(dvtx_coord,      3*dn_vtx)
  np_dface_vtx_idx   = create_numpy_i(dface_vtx_idx,   dn_face+1)
  np_dface_vtx       = create_numpy_g(dface_vtx,       np_dface_vtx_idx[np_dface_vtx_idx.shape[0]-1])
  # > In 2D cas number of vertex is the same as number of edge
  np_dface_edge      = create_numpy_g(dface_edge,      np_dface_vtx_idx[np_dface_vtx_idx.shape[0]-1])
  np_dedge_vtx       = create_numpy_g(dedge_vtx,       2*dn_edge)
  np_dedge_face      = create_numpy_g(dedge_face,      2*dn_edge)
  np_dedge_group_idx = create_numpy_i(dedge_group_idx, n_edge_group+1)
  np_dedge_group     = create_numpy_g(dedge_group,     np_dedge_group_idx[np_dedge_group_idx.shape[0]-1])

  return {'n_face'          : n_face,
          'n_vtx'           : n_vtx,
          'n_edge'          : n_edge,
          'dn_face'         : dn_face,
          'dn_vtx'          : dn_vtx,
          'dn_edge'         : dn_edge,
          'n_edge_group'    : n_edge_group,
          'dvtx_coord'      : np_dvtx_coord,
          'dface_vtx_idx'   : np_dface_vtx_idx,
          'dface_vtx'       : np_dface_vtx,
          'dface_edge'      : np_dface_edge,
          'dedge_vtx'       : np_dedge_vtx,
          'dedge_face'      : np_dedge_face,
          'dedge_group_idx' : np_dedge_group_idx,
          'dedge_group'     : np_dedge_group,
          }
