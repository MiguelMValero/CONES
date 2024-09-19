
cdef extern from "pdm_dmesh.h":
    # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    # > Wrapping of structure
    ctypedef struct PDM_dmesh_t:
      pass
    # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    # > Wrapping of function
    PDM_dmesh_t* PDM_dmesh_create(PDM_ownership_t owner,
                                  int             dn_cell,
                                  int             dn_face,
                                  int             dn_edge,
                                  int             dn_vtx,
                                  PDM_MPI_Comm    comm)

    # void PDM_dmesh_set(PDM_dmesh_t  *dm,
    #                    double       *dvtx_coord,
    #                    int          *dface_vtx_idx,
    #                    PDM_g_num_t  *dface_vtx,
    #                    PDM_g_num_t  *dface_cell,
    #                    int          *dface_bound_idx,
    #                    PDM_g_num_t  *dface_bound,
    #                    int          *join_g_dms,
    #                    int          *dface_join_idx,
    #                    PDM_g_num_t  *dface_join)

    # void PDM_dmesh_dims_get(PDM_dmesh_t *dm,
    #                         int         *dn_cell,
    #                         int         *dn_face,
    #                         int         *dn_edge,
    #                         int         *dn_vtx,
    #                         int         *n_bnd,
    #                         int         *n_joins)

    void PDM_dmesh_data_get(PDM_dmesh_t   *dm,
                            double       **dvtx_coord,
                            int          **dface_vtx_idx,
                            PDM_g_num_t  **dface_vtx,
                            PDM_g_num_t  **dface_cell,
                            int          **dface_bound_idx,
                            PDM_g_num_t  **dface_bound,
                            int          **join_g_dms,
                            int          **dface_join_idx,
                            PDM_g_num_t  **dface_join)

    int PDM_dmesh_dn_entity_get(PDM_dmesh_t         *dmesh,
                                PDM_mesh_entities_t  entity_type)
    int PDM_dmesh_connectivity_get(PDM_dmesh_t              *dmesh,
                                   PDM_connectivity_type_t   connectivity_type,
                                   PDM_g_num_t             **connect,
                                   int                     **connect_idx,
                                   PDM_ownership_t           ownership)
    int PDM_dmesh_distrib_get(PDM_dmesh_t              *dmesh,
                               PDM_mesh_entities_t       entity,
                               PDM_g_num_t             **distrib)
    int  PDM_dmesh_bound_get(PDM_dmesh_t       *dmesh,
                             PDM_bound_type_t   bound_type,
                             PDM_g_num_t      **connect,
                             int              **connect_idx,
                             PDM_ownership_t    ownership)
    void PDM_dmesh_vtx_coord_set(PDM_dmesh_t      *dmesh,
                                 double           *dvtx_coord,
                                 PDM_ownership_t   ownership)
    void PDM_dmesh_vtx_coord_get(PDM_dmesh_t      *dmesh,
                                 double          **dvtx_coord,
                                 PDM_ownership_t   ownership)
    void  PDM_dmesh_bound_set(PDM_dmesh_t       *dmesh,
                              PDM_bound_type_t   bound_type,
                              int               n_bound,
                              PDM_g_num_t       *connect,
                              int               *connect_idx,
                              PDM_ownership_t    ownership)
    void PDM_dmesh_connectivity_set(PDM_dmesh_t              *dmesh,
                                    PDM_connectivity_type_t   connectivity_type,
                                    PDM_g_num_t              *connect,
                                    int                      *connect_idx,
                                    PDM_ownership_t           ownership);
    void PDM_dmesh_free(PDM_dmesh_t   *dm)
    # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

# ------------------------------------------------------------------
# cdef void PDM_pydmesh_free(object caps):
#   # print("PDM_pydmesh_free", PyCapsule_GetName(caps))
#   cdef PDM_dmesh_t* dm = <PDM_dmesh_t *> PyCapsule_GetPointer(caps, <const char*> PyCapsule_GetName(caps))
#   PDM_dmesh_free(dm);

# ------------------------------------------------------------------
cdef class DistributedMeshCapsule:
  """
  """
  # ************************************************************************
  # > Class attributes
  cdef PDM_dmesh_t* _dm
  # ************************************************************************
  # ------------------------------------------------------------------------
  def __cinit__(self, object caps):
    """
    """
    # print("DistributedMeshCapsule", PyCapsule_GetName(caps))
    cdef PDM_dmesh_t* dm = <PDM_dmesh_t *> PyCapsule_GetPointer(caps, NULL)
    self._dm = dm;

  # ------------------------------------------------------------------------
  def dmesh_vtx_coord_set(self,
                          NPY.ndarray[NPY.double_t, mode='c', ndim=1] dvtx_coord not None):
      dmesh_vtx_coord_set(self, dvtx_coord)
  def dmesh_vtx_coord_get(self):
      return dmesh_vtx_coord_get(self)
  # ------------------------------------------------------------------------
  def dmesh_connectivity_set(self, PDM_connectivity_type_t connectivity_type,
                             NPY.ndarray[NPY.int32_t   , mode='c', ndim=1] connect_idx,
                             NPY.ndarray[npy_pdm_gnum_t, mode='c', ndim=1] connect    not None):
    """
    """
    dmesh_connectivity_set(self, connectivity_type, connect_idx, connect)

  # ------------------------------------------------------------------------
  def dmesh_connectivity_get(self, PDM_connectivity_type_t connectivity_type):
    """
    """
    return dmesh_connectivity_get(self, connectivity_type)

  # ------------------------------------------------------------------------
  def dmesh_distrib_get(self, PDM_mesh_entities_t entity_type):
    """
    """
    return dmesh_distrib_get(self, entity_type)

  # ------------------------------------------------------------------------
  def dmesh_bound_get(self, PDM_bound_type_t bound_type):
    """
    """
    return dmesh_bound_get(self, bound_type)

  # ------------------------------------------------------------------------
  def dmesh_bound_set(self, PDM_bound_type_t bound_type,
                            NPY.ndarray[NPY.int32_t   , mode='c', ndim=1] connect_idx,
                            NPY.ndarray[npy_pdm_gnum_t, mode='c', ndim=1] connect    not None):
    """
    """
    dmesh_bound_set(self, bound_type, connect_idx, connect)

  # ------------------------------------------------------------------------
  def __dealloc__(self):
    """
       Use the free method of PDM Lib
    """
    PDM_dmesh_free(self._dm)

# ------------------------------------------------------------------
cdef class DistributedMesh:
  """
     DistributedMesh: Distributed mesh structure
  """
  # ************************************************************************
  # > Class attributes
  cdef PDM_dmesh_t* _dm
  cdef int          n_rank
  # ************************************************************************
  # ------------------------------------------------------------------------
  def __cinit__(self, MPI.Comm comm,
                      dn_cell,
                      dn_face,
                      dn_edge,
                      dn_vtx):
    """
    TODOUX
    """
    # ************************************************************************
    # > Declaration
    # ************************************************************************
    cdef MPI.MPI_Comm c_comm = comm.ob_mpi
    self.n_rank = comm.Get_size()
    # ::::::::::::::::::::::::::::::::::::::::::::::::::
    self._dm = PDM_dmesh_create(PDM_OWNERSHIP_UNGET_RESULT_IS_FREE,
                                dn_cell,
                                dn_face,
                                dn_edge,
                                dn_vtx,
                                PDM_MPI_mpi_2_pdm_mpi_comm (<void *> &c_comm))
    # ::::::::::::::::::::::::::::::::::::::::::::::::::

  # ------------------------------------------------------------------------
  # def dmesh_set(self, NPY.ndarray[NPY.double_t  , mode='c', ndim=1] dvtx_coord   not None,
  #                     NPY.ndarray[NPY.int32_t   , mode='c', ndim=1] dface_vtx_idx,
  #                     NPY.ndarray[npy_pdm_gnum_t, mode='c', ndim=1] dface_vtx    not None,
  #                     NPY.ndarray[npy_pdm_gnum_t, mode='c', ndim=1] dface_cell,
  #                     NPY.ndarray[NPY.int32_t   , mode='c', ndim=1] dface_bound_idx,
  #                     NPY.ndarray[npy_pdm_gnum_t, mode='c', ndim=1] dface_bound,
  #                     NPY.ndarray[NPY.int32_t   , mode='c', ndim=1] join_g_dms,
  #                     NPY.ndarray[NPY.int32_t   , mode='c', ndim=1] dface_join_idx,
  #                     NPY.ndarray[npy_pdm_gnum_t, mode='c', ndim=1] dface_join):
  #   """
  #   """
  #   # print("dvtx_coord", dvtx_coord)
  #   # print("dface_vtx_idx", dface_vtx_idx)
  #   # print("dface_vtx", dface_vtx)
  #   # print("dface_cell", dface_cell)
  #   # print("dface_bound_idx", dface_bound_idx)
  #   # print("dface_bound", dface_bound)

  #   PDM_dmesh_set(self._dm,
  #                 <double*>      dvtx_coord.data,
  #                 <int*>         dface_vtx_idx.data,
  #                 <PDM_g_num_t*> dface_vtx.data,
  #                 <PDM_g_num_t*> dface_cell.data,
  #                 <int*>         dface_bound_idx.data,
  #                 <PDM_g_num_t*> dface_bound.data,
  #                 <int*>         join_g_dms.data,
  #                 <int*>         dface_join_idx.data,
  #                 <PDM_g_num_t*> dface_join.data)
  # ------------------------------------------------------------------------
  def dmesh_vtx_coord_set(self,
                          NPY.ndarray[NPY.double_t, mode='c', ndim=1] dvtx_coord not None):
      dmesh_vtx_coord_set(self, dvtx_coord)
  def dmesh_vtx_coord_get(self):
      return dmesh_vtx_coord_get(self)

  # ------------------------------------------------------------------------
  def dmesh_connectivity_get(self, PDM_connectivity_type_t connectivity_type):
    """
    """
    return dmesh_connectivity_get(self, connectivity_type)

  # ------------------------------------------------------------------------
  def dmesh_connectivity_set(self, PDM_connectivity_type_t connectivity_type,
                             NPY.ndarray[NPY.int32_t   , mode='c', ndim=1] connect_idx,
                             NPY.ndarray[npy_pdm_gnum_t, mode='c', ndim=1] connect    not None):
    """
    """
    dmesh_connectivity_set(self, connectivity_type, connect_idx, connect)

  # ------------------------------------------------------------------------
  def dmesh_distrib_get(self, PDM_mesh_entities_t entity_type):
    """
    """
    return dmesh_connectivity_get(self, entity_type)
  # ------------------------------------------------------------------------
  def dmesh_bound_get(self, PDM_bound_type_t bound_type):
    """
    """
    return dmesh_bound_get(self, bound_type)

  # ------------------------------------------------------------------------
  def dmesh_bound_set(self, PDM_bound_type_t bound_type,
                            NPY.ndarray[NPY.int32_t   , mode='c', ndim=1] connect_idx,
                            NPY.ndarray[npy_pdm_gnum_t, mode='c', ndim=1] connect    not None):
    """
    """
    dmesh_bound_set(self, bound_type, connect_idx, connect)

  # ------------------------------------------------------------------------
  def __dealloc__(self):
    """
       Use the free method of PDM Lib
    """
    # ************************************************************************
    # > Declaration
    # ************************************************************************
    # print("DistributedMesh::__dealloc__")
    PDM_dmesh_free(self._dm)
    # print("DistributedMesh::__dealloc__")

ctypedef fused DMesh:
  DistributedMesh
  DistributedMeshCapsule

# ------------------------------------------------------------------------
def dmesh_connectivity_get(DMesh pydm, PDM_connectivity_type_t connectivity_type):
  """
  """
  # ************************************************************************
  # > Declaration
  cdef int          *connect_idx
  cdef PDM_g_num_t  *connect
  cdef int           dn_entity
  cdef NPY.npy_intp  dim
  # ************************************************************************
  dn_entity = PDM_dmesh_connectivity_get(pydm._dm,
                                         connectivity_type,
                                         &connect,
                                         &connect_idx,
                                         PDM_OWNERSHIP_USER)

  np_connect_idx = create_numpy_or_none_i(connect_idx, dn_entity+1)

  if (connect == NULL) :
      np_connect = None
  else :
      if(np_connect_idx is not None):
        dim = <NPY.npy_intp> connect_idx[dn_entity]
      else:
        dim = <NPY.npy_intp> 2 * dn_entity # Face cell
      np_connect = create_numpy_g(connect, dim)
  return (np_connect_idx, np_connect)

# ------------------------------------------------------------------------
def dmesh_distrib_get(DMesh pydm, PDM_mesh_entities_t entity_type):
  """
  """
  # ************************************************************************
  # > Declaration
  cdef PDM_g_num_t  *distrib
  cdef NPY.npy_intp  dim
  cdef int size
  # ************************************************************************
  size = PDM_dmesh_distrib_get(pydm._dm, entity_type, &distrib)
  return create_numpy_or_none_g(distrib, size+1, flag_owndata=False)

# ------------------------------------------------------------------------
def dmesh_bound_get(DMesh pydm, PDM_bound_type_t bound_type):
  """
  """
  # ************************************************************************
  # > Declaration
  cdef int          *connect_idx
  cdef PDM_g_num_t  *connect
  cdef int           dn_entity
  cdef NPY.npy_intp  dim
  cdef int           n_bnd
  # ************************************************************************
  n_bnd = PDM_dmesh_bound_get(pydm._dm,
                              bound_type,
                              &connect,
                              &connect_idx,
                              PDM_OWNERSHIP_USER)

  np_connect_idx = create_numpy_or_none_i(connect_idx, n_bnd+1)

  if (connect == NULL) :
      np_connect = None
  else :
      np_connect = create_numpy_g(connect, connect_idx[n_bnd])
  return (np_connect_idx, np_connect)

# ------------------------------------------------------------------------
def dmesh_bound_set(DMesh pydm,
                    PDM_bound_type_t bound_type,
                    NPY.ndarray[NPY.int32_t   , mode='c', ndim=1] connect_idx,
                    NPY.ndarray[npy_pdm_gnum_t, mode='c', ndim=1] connect    not None
                    ):
  """
  """
  # ************************************************************************
  # > Declaration
  cdef int           n_bnd
  # ************************************************************************

  n_bnd = connect_idx.shape[0]-1

  PDM_dmesh_bound_set(pydm._dm,
                      bound_type,
                      n_bnd,
      <PDM_g_num_t *> connect.data,
              <int *> connect_idx.data,
                      PDM_OWNERSHIP_USER)

# ------------------------------------------------------------------------
def dmesh_connectivity_set(DMesh pydm,
                           PDM_connectivity_type_t entity_type,
                           NPY.ndarray[NPY.int32_t   , mode='c', ndim=1] connect_idx,
                           NPY.ndarray[npy_pdm_gnum_t, mode='c', ndim=1] connect    not None
                           ):
  """
  """
  # ************************************************************************
  # > Declaration
  # ************************************************************************

  cdef int* _connect_idx = NULL
  if connect_idx is not None:
    _connect_idx = <int *> connect_idx.data
  PDM_dmesh_connectivity_set(pydm._dm,
                             entity_type,
             <PDM_g_num_t *> connect.data,
                             _connect_idx,
                             PDM_OWNERSHIP_USER)

# ------------------------------------------------------------------------
def dmesh_vtx_coord_set(DMesh pydm,
                        NPY.ndarray[NPY.double_t, mode='c', ndim=1] dvtx_coord not None):
    PDM_dmesh_vtx_coord_set(pydm._dm,
                 <double *> dvtx_coord.data,
                            PDM_OWNERSHIP_USER)

def dmesh_vtx_coord_get(DMesh pydm):
    dn_vtx = PDM_dmesh_dn_entity_get(pydm._dm, PDM_MESH_ENTITY_VTX)
    cdef double* dvtx_coord = NULL
    PDM_dmesh_vtx_coord_get(pydm._dm, 
                           &dvtx_coord,
                            PDM_OWNERSHIP_USER)
    return create_numpy_d(dvtx_coord, 3*dn_vtx)
