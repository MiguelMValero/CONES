
cdef extern from "pdm_mesh_intersection.h":
    # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    # > Wrapping of Ppart Structure
    ctypedef struct PDM_mesh_intersection_t:
        pass
    # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    # > Wrapping of function
    PDM_mesh_intersection_t* PDM_mesh_intersection_create(PDM_mesh_intersection_kind_t intersection_kind,
                                                          int                          dim_mesh_a,
                                                          int                          dim_mesh_b,
                                                          double                       project_coeff,
                                                          PDM_MPI_Comm                 comm,
                                                          PDM_ownership_t              owner);

    void PDM_mesh_intersection_n_part_set(PDM_mesh_intersection_t *mi,
                                          const int                i_mesh,
                                          const int                n_part);

    void PDM_mesh_intersection_part_set(PDM_mesh_intersection_t  *mi,
                                        PDM_ol_mesh_t             i_mesh,
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

    void PDM_mesh_intersection_compute(PDM_mesh_intersection_t  *mi);

    void PDM_mesh_intersection_part_to_part_get(PDM_mesh_intersection_t  *mi,
                                                PDM_part_to_part_t      **ptp,
                                                PDM_ownership_t           ownership);

    void PDM_mesh_intersection_result_from_a_get(PDM_mesh_intersection_t  *mi,
                                                 const int                 ipart,
                                                 int                     **elt_a_elt_b_idx,
                                                 PDM_g_num_t             **elt_a_elt_b,
                                                 double                  **elt_a_elt_b_weight);


    void PDM_mesh_intersection_result_from_b_get(PDM_mesh_intersection_t  *mi,
                                                 const int                 ipart,
                                                 double                  **elt_b_elt_a_weight);
    void PDM_mesh_intersection_free(PDM_mesh_intersection_t* mi);

# ------------------------------------------------------------------
cdef class MeshIntersection:
    """
    Define a method to compute mesh intersections
    """
    # ************************************************************************
    # > Class attributes
    cdef PDM_mesh_intersection_t* _mi
    cdef MPI.Comm py_comm
    cdef dict     ptp_objects
    cdef list     keep_alive
    cdef int      _dim_mesh_a
    cdef int      _dim_mesh_b
    cdef int      _n_part_mesh_a
    cdef int      _n_part_mesh_b
    cdef object   _n_entity_a
    cdef object   _n_entity_b
    # ************************************************************************
    # ------------------------------------------------------------------
    def __cinit__(self,
                  MPI.Comm                     comm,
                  PDM_mesh_intersection_kind_t intersection_kind,
                  int                          dim_mesh_a,
                  int                          dim_mesh_b,
                  int                          n_part_mesh_a,
                  int                          n_part_mesh_b,
                  double                       project_coeff=1.e-6):
        """
        Compute mesh intersections
        """
        self.keep_alive  = list()

        # Convert mpi4py -> PDM_MPI
        self.py_comm = comm
        cdef MPI.MPI_Comm c_comm   = self.py_comm.ob_mpi
        cdef PDM_MPI_Comm pdm_comm = PDM_MPI_mpi_2_pdm_mpi_comm(&c_comm)

        self._mi = PDM_mesh_intersection_create(intersection_kind,
                                                dim_mesh_a,
                                                dim_mesh_b,
                                                project_coeff,
                                                pdm_comm,
                                                PDM_OWNERSHIP_UNGET_RESULT_IS_FREE);

        self._dim_mesh_a    = dim_mesh_a
        self._dim_mesh_b    = dim_mesh_b
        self._n_part_mesh_a = n_part_mesh_a
        self._n_part_mesh_b = n_part_mesh_b
        self._n_entity_a    = NPY.zeros(n_part_mesh_a, dtype='int32', order='C')
        self._n_entity_b    = NPY.zeros(n_part_mesh_b, dtype='int32', order='C')

        PDM_mesh_intersection_n_part_set(self._mi,
                                         0,
                                         n_part_mesh_a)

        PDM_mesh_intersection_n_part_set(self._mi,
                                         1,
                                         n_part_mesh_b)
    # ------------------------------------------------------------------
    def part_set(self,
                 int                                           i_mesh,
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

      if(i_mesh ==  0):
        if(self._dim_mesh_a == 3):
          self._n_entity_a[i_part] = n_cell
        elif(self._dim_mesh_a == 2):
          self._n_entity_a[i_part] = n_face
        elif(self._dim_mesh_a == 1):
          self._n_entity_a[i_part] = n_edge
      else:
        if(self._dim_mesh_b == 3):
          self._n_entity_b[i_part] = n_cell
        elif(self._dim_mesh_b == 2):
          self._n_entity_b[i_part] = n_face
        elif(self._dim_mesh_b == 1):
          self._n_entity_b[i_part] = n_edge

      PDM_mesh_intersection_part_set(self._mi,
                    <PDM_ol_mesh_t>  i_mesh,
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
                    <PDM_g_num_t *>  cell_ln_to_gn.data,
                    <PDM_g_num_t *>  face_ln_to_gn.data,
                    <PDM_g_num_t *>  edge_ln_to_gn.data,
                    <PDM_g_num_t *>  vtx_ln_to_gn .data,
                    <double      *>  coords       .data)


    # ------------------------------------------------------------------
    def compute(self):
        """
        """
        PDM_mesh_intersection_compute(self._mi)

    # ------------------------------------------------------------------------
    def part_to_part_get(self):
      """
      """
      cdef PDM_part_to_part_t *ptpc
      PDM_mesh_intersection_part_to_part_get(self._mi,
                                             &ptpc,
                                             PDM_OWNERSHIP_USER)

      py_caps = PyCapsule_New(ptpc, NULL, NULL)
      return PartToPartCapsule(py_caps, self.py_comm) # The free is inside the class

    # ------------------------------------------------------------------
    def a_to_b_get(self, int i_part):
      """
      """
      # ************************************************************************
      # > Declaration
      cdef int         *a_to_b_idx
      cdef PDM_g_num_t *a_to_b
      cdef double      *a_to_b_weight
      # ************************************************************************
      # > Get
      PDM_mesh_intersection_result_from_a_get(self._mi,
                                              i_part,
                                              &a_to_b_idx,
                                              &a_to_b,
                                              &a_to_b_weight)

      return {
        "a_to_b_idx"    : create_numpy_i(a_to_b_idx,    self._n_entity_a[i_part] +1),
        "a_to_b"        : create_numpy_g(a_to_b,        a_to_b_idx[self._n_entity_a[i_part] ]),
        "a_to_b_weight" : create_numpy_d(a_to_b_weight, a_to_b_idx[self._n_entity_a[i_part] ])
      }

    # ------------------------------------------------------------------
    def b_to_a_get(self, int i_part):
      """
      """
      # ************************************************************************
      # > Declaration
      cdef double *b_to_a_weight
      # ************************************************************************
      # > Get
      PDM_mesh_intersection_result_from_b_get(self._mi,
                                              i_part,
                                              &b_to_a_weight)

      # get nb of elt_b in part #i_part
      n_b = 0#!!!
      # get size of b_to_a_weight (from ptp?)
      s_b_to_a = 0
      return {
        "b_to_a_weight" : create_numpy_d(b_to_a_weight, s_b_to_a)
      }

    # ------------------------------------------------------------------
    def __dealloc__(self):
        """
        Free a mesh_intersection structure
        """
        PDM_mesh_intersection_free(self._mi)
