# ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# > Wrapping of functions
cdef extern from "pdm_dmesh_extract.h":
  ctypedef struct PDM_dmesh_extract_t:
    pass

  PDM_dmesh_extract_t* PDM_dmesh_extract_create(int                     dim,
                                                PDM_MPI_Comm            comm);

  void PDM_dmesh_extract_compute(PDM_dmesh_extract_t *dme);

  void PDM_dmesh_extract_selected_gnum_set(PDM_dmesh_extract_t *dme,
                                           PDM_mesh_entities_t  entity_type,
                                           int                  n_selected,
                                           PDM_g_num_t         *selected_gnum);

  void PDM_dmesh_extract_dn_entity_set(PDM_dmesh_extract_t *dme,
                                       PDM_mesh_entities_t  entity_type,
                                       int                  dn_entity);

  void PDM_dmesh_extract_vtx_coord_set(PDM_dmesh_extract_t *dme,
                                       double              *dvtx_coord);

  void PDM_dmesh_extract_dmesh_bound_set(PDM_dmesh_extract_t *dme,
                                         PDM_bound_type_t     bound_type,
                                         int                  n_bound,
                                         PDM_g_num_t         *connect,
                                         int                 *connect_idx);

  void PDM_dmesh_extract_dconnectivity_set(PDM_dmesh_extract_t     *dme,
                                           PDM_connectivity_type_t  connectivity_type,
                                           PDM_g_num_t             *dconnect,
                                           int                     *dconnect_idx);

  void PDM_dmesh_extract_dmesh_get(PDM_dmesh_extract_t     *dme,
                                   PDM_dmesh_t            **dmesh_extract,
                                   PDM_ownership_t          ownership);

  void PDM_dmesh_extract_dmesh_nodal_get(PDM_dmesh_extract_t     *dme,
                                         PDM_dmesh_nodal_t      **dmesh_nodal_extract,
                                         PDM_ownership_t          ownership);

  void PDM_dmesh_extract_parent_gnum_get(PDM_dmesh_extract_t     *dme,
                                         PDM_mesh_entities_t      entity_type,
                                         int                     *dn_entity,
                                         PDM_g_num_t            **parent_gnum,
                                         PDM_ownership_t          ownership);

  void PDM_dmesh_extract_btp_get(PDM_dmesh_extract_t     *dme,
                                 PDM_mesh_entities_t      entity_type,
                                 PDM_block_to_part_t    **btp,
                                 PDM_ownership_t          ownership);

  void PDM_dmesh_extract_btp_group_get(PDM_dmesh_extract_t     *dme,
                                       int                      i_group,
                                       PDM_bound_type_t         bound_type,
                                       PDM_block_to_part_t    **btp,
                                       PDM_ownership_t          ownership);
  void PDM_dmesh_extract_dmesh_set(PDM_dmesh_extract_t     *dme,
                                   PDM_dmesh_t             *dmesh);
  void PDM_dmesh_extract_dmesh_nodal_set(PDM_dmesh_extract_t     *dme,
                                         PDM_dmesh_nodal_t       *dmesh_nodal);

  void PDM_dmesh_extract_free(PDM_dmesh_extract_t  *dme);

# ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# > Class definition
cdef class DMeshExtract:
  """

  """
  # --------------------------------------------------------------------------
  # > Class attributes
  cdef PDM_dmesh_extract_t* _dme
  cdef MPI.Comm py_comm
  cdef dict ptp_objects
  cdef list keep_alive
  # --------------------------------------------------------------------------

  # ------------------------------------------------------------------
  def __cinit__(self,
                int                     dim,
                MPI.Comm                comm):
    """
    Compute the distance from point clouds to a surface
    """
    self.ptp_objects = dict()
    self.keep_alive  = list()

    self.py_comm = comm
    cdef MPI.MPI_Comm c_comm = comm.ob_mpi
    self._dme =  PDM_dmesh_extract_create(dim,
                                          PDM_MPI_mpi_2_pdm_mpi_comm (<void *> &c_comm));

  # ------------------------------------------------------------------
  def compute(self):
    """
    """
    PDM_dmesh_extract_compute(self._dme)

  # ------------------------------------------------------------------
  def register_dmesh(self, DMesh dm):
    """
    """
    PDM_dmesh_extract_dmesh_set(self._dme, dm._dm)

  # ------------------------------------------------------------------
  def register_dmesh_nodal(self, DMeshNodal dmn):
    PDM_dmesh_extract_dmesh_nodal_set(self._dme, dmn.dmn)

  # ------------------------------------------------------------------
  def set_gnum_to_extract(self,
                          entity_type, 
                          NPY.ndarray[npy_pdm_gnum_t, mode='c', ndim=1] selected_gnum):

    PDM_dmesh_extract_selected_gnum_set(self._dme,
                 <PDM_mesh_entities_t>  entity_type,
                                        selected_gnum.size,
                        <PDM_g_num_t*>  selected_gnum.data)

  def get_extract_parent_gnum(self, entity_type):

    cdef int dn_entity = -1
    cdef PDM_g_num_t* parent_gnum = NULL

    PDM_dmesh_extract_parent_gnum_get(self._dme,
                <PDM_mesh_entities_t> entity_type,
                                      &dn_entity,
                                      &parent_gnum,
                                      PDM_OWNERSHIP_USER)

    return create_numpy_g(parent_gnum, dn_entity)

  # ------------------------------------------------------------------
  def get_dmesh(self):
    """
    """
    # ************************************************************************
    # > Declaration
    cdef PDM_dmesh_t* dm
    # ************************************************************************
    PDM_dmesh_extract_dmesh_get(self._dme, &dm, PDM_OWNERSHIP_USER)

    py_caps = PyCapsule_New(dm, NULL, NULL);

    return DistributedMeshCapsule(py_caps) # The free is inside the class

  def get_dmesh_nodal(self):
    cdef PDM_dmesh_nodal_t* dmn
    PDM_dmesh_extract_dmesh_nodal_get(self._dme, &dmn, PDM_OWNERSHIP_USER)
    py_caps = PyCapsule_New(dmn, NULL, NULL);

    return DistributedMeshNodalCapsule(py_caps) # The free is inside the class

  # ------------------------------------------------------------------
  def __dealloc__(self):
    PDM_dmesh_extract_free(self._dme)
