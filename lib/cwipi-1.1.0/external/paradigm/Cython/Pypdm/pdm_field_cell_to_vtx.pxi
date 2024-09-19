
cdef extern from "pdm_field_cell_to_vtx.h":
    # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    # > Wrapping of Ppart Structure
    ctypedef struct PDM_field_cell_to_vtx_t:
        pass
    # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    # > Wrapping of function
    PDM_field_cell_to_vtx_t* PDM_field_cell_to_vtx_create(int                            n_domain,
                                                          int                           *n_part,
                                                          int                           *n_group,
                                                          PDM_cell_to_vtx_interp_kind_t  interp_kind,
                                                          int                            n_depth,
                                                          PDM_MPI_Comm                   comm);

    void PDM_field_cell_to_vtx_compute(PDM_field_cell_to_vtx_t *part_ext);

    void PDM_field_cell_to_vtx_part_set(PDM_field_cell_to_vtx_t   *mi,
                                       int                       i_domain,
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

    void PDM_field_cell_to_vtx_graph_comm_set(PDM_field_cell_to_vtx_t   *mi,
                                             int                       i_domain,
                                             int                       i_part,
                                             PDM_mesh_entities_t       mesh_entity,
                                             int                      *entity_part_bound_proc_idx,
                                             int                      *entity_part_bound_part_idx,
                                             int                      *entity_part_bound);

    void PDM_field_cell_to_vtx_part_group_set(PDM_field_cell_to_vtx_t   *mi,
                                             int                       i_domain,
                                             int                       i_part,
                                             int                       i_group,
                                             PDM_bound_type_t          bound_type,
                                             int                       n_group_entity,
                                             int                      *group_entity);

    void PDM_field_cell_to_vtx_part_domain_interface_shared_set(PDM_field_cell_to_vtx_t      *mi,
                                                               PDM_part_domain_interface_t *pdi);

    # void PDM_field_cell_to_vtx_part_mesh_nodal_set(PDM_field_cell_to_vtx_t     *mi,
    #                                               PDM_part_mesh_nodal_t      *pmn);

    void PDM_field_cell_to_vtx_exch(PDM_field_cell_to_vtx_t      *mi,
                                   PDM_field_kind_t            field_kind,
                                   double                   ***local_field,
                                   double                  ****bound_field,
                                   double                  ****result_field);

    void PDM_field_cell_to_vtx_free(PDM_field_cell_to_vtx_t *mi);


# ------------------------------------------------------------------
cdef class MeshCellToNode:
    """
       PartToBlock: Interface for block_to_part.c
    """
    # ************************************************************************
    # > Class attributes
    cdef PDM_field_cell_to_vtx_t* mi
    cdef int                     n_part
    cdef int*                    pn_elt
    cdef MPI.Comm                py_comm
    cdef object                  pdi
    # ************************************************************************

    # ------------------------------------------------------------------------
    def __cinit__(self, MPI.Comm comm,
                        int n_domain,
                        NPY.ndarray[NPY.int32_t   , mode='c', ndim=1] n_part,
                        NPY.ndarray[NPY.int32_t   , mode='c', ndim=1] n_group,
                        PDM_cell_to_vtx_interp_kind_t interp_kind,
                        int n_depth = 1):
        """
        Constructor of PartToBlock object
        """
        cdef MPI.MPI_Comm c_comm = comm.ob_mpi
        self.mi = PDM_field_cell_to_vtx_create(n_domain,
                                        <int*> n_part.data,
                                        <int*> n_group.data,
                                               interp_kind,
                                               n_depth,
                                               PDM_MPI_mpi_2_pdm_mpi_comm (<void *> &c_comm))

    # ------------------------------------------------------------------
    def compute(self):
        """
        """
        PDM_field_cell_to_vtx_compute(self.mi)

    # ------------------------------------------------------------------
    def part_set(self,
                 int                                           i_domain,
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

      PDM_field_cell_to_vtx_part_set(self.mi,
                                    i_domain,
                                    i_part,
                                    n_cell,
                                    n_face,
                                    n_edge,
                                    n_vtx,
                   <int         *>  cell_face_idx.data,
                   <int         *>  cell_face    .data,
                   <int         *>  face_edge_idx.data,
                   <int         *>  face_edge    .data,
                   <int         *>  edge_vtx     .data,
                   <int         *>  face_vtx_idx .data,
                   <int         *>  face_vtx     .data,
                   <PDM_g_num_t *> cell_ln_to_gn.data,
                   <PDM_g_num_t *> face_ln_to_gn.data,
                   <PDM_g_num_t *> edge_ln_to_gn.data,
                   <PDM_g_num_t *> vtx_ln_to_gn .data,
                   <double      *> coords       .data)

    # ------------------------------------------------------------------
    def graph_comm(self,
                   int                                        i_domain,
                   int                                        i_part,
                   PDM_mesh_entities_t                        mesh_entity,
                   NPY.ndarray[NPY.int32_t, mode='c', ndim=1] entity_part_bound_proc_idx,
                   NPY.ndarray[NPY.int32_t, mode='c', ndim=1] entity_part_bound_part_idx,
                   NPY.ndarray[NPY.int32_t, mode='c', ndim=1] entity_part_bound):
      """
      """
      PDM_field_cell_to_vtx_graph_comm_set(self.mi,
                                          i_domain,
                                          i_part,
                                          mesh_entity,
                                 <int*>   entity_part_bound_proc_idx.data,
                                 <int*>   entity_part_bound_part_idx.data,
                                 <int*>   entity_part_bound.data)

    # ------------------------------------------------------------------
    def part_group_set(self,
                       int                                        i_domain,
                       int                                        i_part,
                       int                                        i_group,
                       PDM_bound_type_t                           bound_type,
                       int                                        n_group_entity,
                       NPY.ndarray[NPY.int32_t, mode='c', ndim=1] group_entity):
      """
      """
      PDM_field_cell_to_vtx_part_group_set(self.mi,
                                          i_domain,
                                          i_part,
                                          i_group,
                                          bound_type,
                                          n_group_entity,
                                <int*>    group_entity.data)

    # ------------------------------------------------------------------
    def part_domain_interface_shared_set(self, PartDomainInterface pdi):
      """
      """
      self.pdi = pdi # Keep alive
      PDM_field_cell_to_vtx_part_domain_interface_shared_set(self.mi, pdi.pdi)


    # ------------------------------------------------------------------
    def exch(self):
        """
        """
        # PDM_field_cell_to_vtx_exch(PDM_field_cell_to_vtx_t      *mi,
        #                           PDM_field_kind_t            field_kind,
        #                           double                   ***local_field,
        #                           double                  ****bound_field,
        #                           double                  ****result_field);


    # ------------------------------------------------------------------
    def __dealloc__(self):
        PDM_field_cell_to_vtx_free(self.mi)
