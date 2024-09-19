# ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# > Wrapping of functions
cdef extern from "pdm_extract_part.h":
  ctypedef struct PDM_extract_part_t:
    pass
  PDM_extract_part_t* PDM_extract_part_create(int                     dim,
                                              int                     n_part_in,
                                              int                     n_part_out,
                                              PDM_extract_part_kind_t extract_kind,
                                              PDM_split_dual_t        split_dual_method,
                                              PDM_bool_t              compute_child_gnum,
                                              PDM_ownership_t         ownership,
                                              PDM_MPI_Comm            comm);
  void PDM_extract_part_compute(PDM_extract_part_t        *extrp);
  void PDM_extract_part_selected_lnum_set(PDM_extract_part_t       *extrp,
                                          int                       i_part,
                                          int                       n_extract,
                                          int                      *extract_lnum);

  void PDM_extract_part_target_set(PDM_extract_part_t       *extrp,
                                   int                       i_part,
                                   int                       n_target,
                                   PDM_g_num_t              *target_gnum,
                                   int                      *target_location);

  void PDM_extract_part_part_set(PDM_extract_part_t        *extrp,
                                 int                       i_part,
                                 int                       n_cell,
                                 int                       n_face,
                                 int                       n_edge,
                                 int                       n_vtx,
                                 int                      *cell_face_idx,
                                 int                      *cell_face,
                                 int                      *face_edge_idx,
                                 int                      *face_edge,
                                 int                      *edge_vtx,
                                 int                      *face_vtx_idx,
                                 int                      *face_vtx,
                                 PDM_g_num_t              *cell_ln_to_gn,
                                 PDM_g_num_t              *face_ln_to_gn,
                                 PDM_g_num_t              *edge_ln_to_gn,
                                 PDM_g_num_t              *vtx_ln_to_gn,
                                 double                   *vtx_coord);

  int PDM_extract_part_n_entity_get(PDM_extract_part_t       *extrp,
                                    int                       i_part_out,
                                    PDM_mesh_entities_t       entity_type);

  int PDM_extract_part_connectivity_get(PDM_extract_part_t        *extrp,
                                        int                        i_part_out,
                                        PDM_connectivity_type_t    connectivity_type,
                                        int                      **connect,
                                        int                      **connect_idx,
                                        PDM_ownership_t           ownership);
  int PDM_extract_part_ln_to_gn_get(PDM_extract_part_t        *extrp,
                                    int                        i_part_out,
                                    PDM_mesh_entities_t        entity_type,
                                    PDM_g_num_t              **pentity_ln_to_gn,
                                    PDM_ownership_t            ownership);

  int PDM_extract_part_parent_ln_to_gn_get(PDM_extract_part_t        *extrp,
                                           int                        i_part_out,
                                           PDM_mesh_entities_t        entity_type,
                                           PDM_g_num_t              **parent_entity_ln_to_gn,
                                           PDM_ownership_t            ownership);

  int PDM_extract_part_vtx_coord_get(PDM_extract_part_t         *extrp,
                                     int                        i_part_out,
                                     double                   **pvtx_coord,
                                     PDM_ownership_t            ownership);

  void PDM_extract_part_part_to_part_get(       PDM_extract_part_t   *extrp,
                                         const  PDM_mesh_entities_t   entity_type,
                                                PDM_part_to_part_t  **ptp,
                                                PDM_ownership_t       ownership);

  void PDM_extract_part_free(PDM_extract_part_t  *extrp);

  void PDM_extract_part_n_group_set(PDM_extract_part_t        *extrp,
                                    PDM_bound_type_t           bound_type,
                                    int                        n_group);

  void PDM_extract_part_part_group_set(PDM_extract_part_t        *extrp,
                                       int                       i_part,
                                       int                       i_group,
                                       PDM_bound_type_t          bound_type,
                                       int                       n_group_entity,
                                       int                      *group_entity,
                                       PDM_g_num_t              *group_entity_ln_to_gn);


  void PDM_extract_part_part_to_part_group_get(PDM_extract_part_t   *extrp,
                                               PDM_bound_type_t      bound_type,
                                               int                   i_group,
                                               PDM_part_to_part_t  **ptp,
                                               PDM_ownership_t       ownership);

  void PDM_extract_part_group_get(PDM_extract_part_t   *extrp,
                                  PDM_bound_type_t      bound_type,
                                  int                   i_part,
                                  int                   i_group,
                                  int                  *pn_extract_group_entity,
                                  int                 **pextract_group_entity,
                                  PDM_g_num_t         **pextract_group_entity_ln_to_gn,
                                  PDM_g_num_t         **pextract_group_entity_parent_ln_to_gn,
                                  PDM_ownership_t       ownership);
# ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

# ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# > Class definition
cdef class ExtractPart:
  """

  """
  # --------------------------------------------------------------------------
  # > Class attributes
  cdef PDM_extract_part_t* _extrp
  cdef MPI.Comm py_comm   
  cdef dict ptp_objects
  cdef dict ptp_group_objects
  cdef list keep_alive
  # --------------------------------------------------------------------------

  # ------------------------------------------------------------------
  def __cinit__(self,
                int                     dim,
                int                     n_part_in,
                int                     n_part_out,
                PDM_extract_part_kind_t extract_kind,
                PDM_split_dual_t        split_dual_method,
                PDM_bool_t              compute_child_gnum,
                MPI.Comm                comm):
    """
    Compute the distance from point clouds to a surface
    """
    self.ptp_objects = dict()
    self.ptp_group_objects = dict()
    self.keep_alive  = list()

    self.py_comm = comm
    cdef MPI.MPI_Comm c_comm = comm.ob_mpi
    self._extrp =  PDM_extract_part_create(dim,
                                           n_part_in,
                                           n_part_out,
                                           extract_kind,
                                           split_dual_method,
                                           compute_child_gnum,
                                           PDM_OWNERSHIP_USER,
                                           PDM_MPI_mpi_2_pdm_mpi_comm (<void *> &c_comm));

  # ------------------------------------------------------------------
  def selected_lnum_set(self,
                        int i_part,
                        NPY.ndarray[NPY.int32_t, mode='c', ndim=1] extract_lnum):
    """
    """
    self.keep_alive.append(extract_lnum)
    cdef int n_extract = extract_lnum.shape[0]
    PDM_extract_part_selected_lnum_set(self._extrp,
                                       i_part,
                                       n_extract,
                              <int* >  extract_lnum.data);


  # ------------------------------------------------------------------
  def target_set(self,
                 int i_part,
                 NPY.ndarray[npy_pdm_gnum_t, mode='c', ndim=1] target_gnum,
                 NPY.ndarray[NPY.int32_t   , mode='c', ndim=1] target_location):
    """
    """
    self.keep_alive.append(target_gnum)
    self.keep_alive.append(target_location)

    cdef int n_target = target_gnum.shape[0]
    PDM_extract_part_target_set(self._extrp,
                                i_part,
                                n_target,
                <PDM_g_num_t *> target_gnum.data,
                <int         *> target_location.data);

  # ------------------------------------------------------------------
  def part_set(self,
               int                                           i_part,
               int                                           n_cell,
               int                                           n_face,
               int                                           n_edge,
               int                                           n_vtx,
               NPY.ndarray[NPY.int32_t   , mode='c', ndim=1] cell_face_idx,
               NPY.ndarray[NPY.int32_t   , mode='c', ndim=1] cell_face    ,
               NPY.ndarray[NPY.int32_t   , mode='c', ndim=1] face_edge_idx,
               NPY.ndarray[NPY.int32_t   , mode='c', ndim=1] face_edge    ,
               NPY.ndarray[NPY.int32_t   , mode='c', ndim=1] edge_vtx     ,
               NPY.ndarray[NPY.int32_t   , mode='c', ndim=1] face_vtx_idx ,
               NPY.ndarray[NPY.int32_t   , mode='c', ndim=1] face_vtx     ,
               NPY.ndarray[npy_pdm_gnum_t, mode='c', ndim=1] cell_ln_to_gn,
               NPY.ndarray[npy_pdm_gnum_t, mode='c', ndim=1] face_ln_to_gn,
               NPY.ndarray[npy_pdm_gnum_t, mode='c', ndim=1] edge_ln_to_gn,
               NPY.ndarray[npy_pdm_gnum_t, mode='c', ndim=1] vtx_ln_to_gn ,
               NPY.ndarray[NPY.double_t  , mode='c', ndim=1] coords):
    """
    """
    self.keep_alive.append(cell_face_idx)
    self.keep_alive.append(cell_face)
    self.keep_alive.append(face_edge_idx)
    self.keep_alive.append(face_edge)
    self.keep_alive.append(edge_vtx)
    self.keep_alive.append(face_vtx_idx)
    self.keep_alive.append(face_vtx)
    self.keep_alive.append(cell_ln_to_gn)
    self.keep_alive.append(face_ln_to_gn)
    self.keep_alive.append(edge_ln_to_gn)
    self.keep_alive.append(vtx_ln_to_gn)
    self.keep_alive.append(coords)

    cdef int * _face_edge_idx;
    if (face_edge_idx is None):
        _face_edge_idx = NULL
    else:
        _face_edge_idx = <int *> face_edge_idx.data

    cdef int * _face_edge;
    if (face_edge is None):
        _face_edge = NULL
    else:
        _face_edge = <int *> face_edge.data

    cdef int * _edge_vtx;
    if (edge_vtx is None):
        _edge_vtx = NULL
    else:
        _edge_vtx = <int *> edge_vtx.data

    cdef PDM_g_num_t* _edge_ln_to_gn;
    if (edge_ln_to_gn is None):
        _edge_ln_to_gn = NULL
    else:
        _edge_ln_to_gn = <PDM_g_num_t *> edge_ln_to_gn.data

    PDM_extract_part_part_set(self._extrp,
                              i_part,
                              n_cell,
                              n_face,
                              n_edge,
                              n_vtx,
             <int         *>  cell_face_idx.data,
             <int         *>  cell_face    .data,
                             _face_edge_idx,
                             _face_edge,
                             _edge_vtx,
             <int         *>  face_vtx_idx .data,
             <int         *>  face_vtx     .data,
             <PDM_g_num_t *> cell_ln_to_gn.data,
             <PDM_g_num_t *> face_ln_to_gn.data,
                            _edge_ln_to_gn,
             <PDM_g_num_t *> vtx_ln_to_gn .data,
             <double      *> coords       .data)

  # ------------------------------------------------------------------
  def part_n_group_set(self,
                       PDM_bound_type_t bound_type,
                       int              n_group):
    """
    """
    PDM_extract_part_n_group_set(self._extrp,
                                 bound_type,
                                 n_group)


  # ------------------------------------------------------------------
  def part_group_set(self,
                     int                                           i_part,
                     int                                           i_group,
                     PDM_bound_type_t                              bound_type,
                     NPY.ndarray[NPY.int32_t   , mode='c', ndim=1] np_group_entity,
                     NPY.ndarray[npy_pdm_gnum_t, mode='c', ndim=1] np_group_entity_ln_to_gn):
    """
    """
    self.keep_alive.append(np_group_entity)
    self.keep_alive.append(np_group_entity_ln_to_gn)
    cdef int n_group_entity = np_group_entity.shape[0]
    PDM_extract_part_part_group_set(self._extrp,
                                    i_part,
                                    i_group,
                                    bound_type,
                                    n_group_entity,
                   <int        * >  np_group_entity.data,
                   <PDM_g_num_t* >  np_group_entity_ln_to_gn.data);

  # ------------------------------------------------------------------
  def compute(self):
    """
    """
    PDM_extract_part_compute(self._extrp)

  # ------------------------------------------------------------------
  def n_entity_get(self,
                   int ipart,
                   PDM_mesh_entities_t entity_type):
    """
    """
    return PDM_extract_part_n_entity_get(self._extrp, ipart, entity_type)

  # ------------------------------------------------------------------
  def connectivity_get(self,
                       int ipart,
                       PDM_connectivity_type_t connectivity_type):
    """
    """
    cdef int  n_entity
    cdef int *connect
    cdef int *connect_idx

    n_entity = PDM_extract_part_connectivity_get(self._extrp,
                                                 ipart,
                                                 connectivity_type,
                                                 &connect,
                                                 &connect_idx,
                                                 PDM_OWNERSHIP_USER)

    np_connect_idx = None
    np_connect     = None
    if(connect_idx != NULL):
      np_connect_idx = create_numpy_i(connect_idx, n_entity+1           )
      np_connect     = create_numpy_i(connect    , connect_idx[n_entity])
    else:
      np_connect  = create_numpy_i(connect    , 2 * n_entity)

    # return PDM_extract_part_n_entity_get(self._extrp, ipart, entity_type)
    return np_connect_idx, np_connect

  # ------------------------------------------------------------------
  def ln_to_gn_get(self,
                   int ipart,
                   PDM_mesh_entities_t entity_type):
    """
    """
    cdef int  n_entity
    cdef PDM_g_num_t *entity_ln_to_gn

    n_entity = PDM_extract_part_ln_to_gn_get(self._extrp,
                                             ipart,
                                             entity_type,
                                             &entity_ln_to_gn,
                                             PDM_OWNERSHIP_USER)
    return create_numpy_g(entity_ln_to_gn, n_entity)

  # ------------------------------------------------------------------
  def parent_ln_to_gn_get(self,
                          int ipart,
                          PDM_mesh_entities_t entity_type):
    """
    """
    cdef int  n_entity
    cdef PDM_g_num_t *parent_ln_to_gn

    n_entity = PDM_extract_part_parent_ln_to_gn_get(self._extrp,
                                                    ipart,
                                                    entity_type,
                                                    &parent_ln_to_gn,
                                                    PDM_OWNERSHIP_USER)
    return create_numpy_g(parent_ln_to_gn, n_entity)

  # ------------------------------------------------------------------
  def vtx_coord_get(self,
                    int ipart):
    """
    """
    cdef int     n_vtx
    cdef double *pvtx_coord

    n_vtx = PDM_extract_part_vtx_coord_get(self._extrp, ipart, &pvtx_coord, PDM_OWNERSHIP_USER)

    return create_numpy_d(pvtx_coord, 3 * n_vtx)

  # ------------------------------------------------------------------
  def part_to_part_get(                    self,
                       PDM_mesh_entities_t entity_type):
    """
    """
    cdef PDM_part_to_part_t  *ptpc
    try:
      return self.ptp_objects[entity_type]
    except KeyError:
      PDM_extract_part_part_to_part_get( self._extrp,
                                         entity_type,
                                        &ptpc,
                                         PDM_OWNERSHIP_USER)
      py_caps = PyCapsule_New(ptpc, NULL, NULL);
      self.ptp_objects[entity_type] = PartToPartCapsule(py_caps, self.py_comm) # The free is inside the class
      return self.ptp_objects[entity_type]

  # ------------------------------------------------------------------
  def part_to_part_group_get(                 self,
                             PDM_bound_type_t bound_type,
                             int              i_group):
    """
    """
    cdef PDM_part_to_part_t  *ptpc
    try:
      return self.ptp_group_objects[(i_group, bound_type)]
    except KeyError:
      PDM_extract_part_part_to_part_group_get( self._extrp,
                                               bound_type,
                                               i_group,
                                              &ptpc,
                                               PDM_OWNERSHIP_USER)
      py_caps = PyCapsule_New(ptpc, NULL, NULL);
      self.ptp_group_objects[(i_group, bound_type)] = PartToPartCapsule(py_caps, self.py_comm) # The free is inside the class
      return self.ptp_group_objects[(i_group, bound_type)]

  # ------------------------------------------------------------------
  def extract_part_group_get(self,
                             int              ipart,
                             int              i_group,
                             PDM_bound_type_t bound_type):
    """
    """
    cdef int          pn_extract_group_entity
    cdef int         *pextract_group_entity
    cdef PDM_g_num_t *pextract_group_entity_ln_to_gn
    cdef PDM_g_num_t *pextract_group_entity_parent_ln_to_gn

    PDM_extract_part_group_get(self._extrp,
                               bound_type,
                               ipart,
                               i_group,
                               &pn_extract_group_entity,
                               &pextract_group_entity,
                               &pextract_group_entity_ln_to_gn,
                               &pextract_group_entity_parent_ln_to_gn,
                               PDM_OWNERSHIP_USER)
    return {'group_entity'                 : create_numpy_i(pextract_group_entity                , pn_extract_group_entity),
            'group_entity_ln_to_gn'        : create_numpy_g(pextract_group_entity_ln_to_gn       , pn_extract_group_entity),
            'group_entity_parent_ln_to_gn' : create_numpy_g(pextract_group_entity_parent_ln_to_gn, pn_extract_group_entity)}

  # ------------------------------------------------------------------
  def __dealloc__(self):
      PDM_extract_part_free(self._extrp)
