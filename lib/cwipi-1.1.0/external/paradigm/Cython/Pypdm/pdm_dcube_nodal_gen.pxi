cdef extern from "pdm_dcube_nodal_gen.h":
    # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    # > Wrapping of Dcube Nodal Gen Structure
    ctypedef struct PDM_dcube_nodal_t:
      pass
    # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    # ------------------------------------------------------------------
    PDM_dcube_nodal_t* PDM_dcube_nodal_gen_create(      PDM_MPI_Comm   comm,
                                                  const PDM_g_num_t    nx,
                                                  const PDM_g_num_t    ny,
                                                  const PDM_g_num_t    nz,
                                                  const double         length,
                                                  const double         zero_x,
                                                  const double         zero_y,
                                                  const double         zero_z,
                                                  PDM_Mesh_nodal_elt_t t_elt,
                                                  int                  order,
                                                  PDM_ownership_t      owner)

    # ------------------------------------------------------------------
    void PDM_dcube_nodal_gen_ordering_set(PDM_dcube_nodal_t *dcube,
                                          char              *ordering)

    PDM_dmesh_nodal_t *PDM_dcube_nodal_gen_build(PDM_dcube_nodal_t *dcube)
    void PDM_dcube_nodal_gen_free(PDM_dcube_nodal_t        *pdm_dcube)
    PDM_dmesh_nodal_t* PDM_dcube_nodal_gen_dmesh_nodal_get(PDM_dcube_nodal_t  *dcube)
    # ------------------------------------------------------------------

    # ------------------------------------------------------------------
    void PDM_dcube_nodal_cart_topo(PDM_MPI_Comm              comm,
                                   int                       n_dom_i,
                                   int                       n_dom_j,
                                   int                       n_dom_k,
                                   int                       periodic_i,
                                   int                       periodic_j,
                                   int                       periodic_k,
                                   const PDM_g_num_t         n_vtx_x_in,
                                   const PDM_g_num_t         n_vtx_y_in,
                                   const PDM_g_num_t         n_vtx_z_in,
                                   const double              length,
                                   const double              zero_x,
                                   const double              zero_y,
                                   const double              zero_z,
                                   PDM_Mesh_nodal_elt_t      t_elt,
                                   const int                 order,
                                   PDM_dcube_nodal_t      ***dcube,
                                   PDM_domain_interface_t  **dom_intrf,
                                   PDM_ownership_t           owner);
    # ------------------------------------------------------------------

    # ------------------------------------------------------------------
    void PDM_dcube_nodal_gen_random_factor_set(PDM_dcube_nodal_t *dcube,
                                               double             random_factor)
    # ------------------------------------------------------------------

# ------------------------------------------------------------------
cdef class DCubeNodalGenerator:
    # > For Dcube Nodal Gen
    cdef PDM_dcube_nodal_t* _dcube
    # ------------------------------------------------------------------
    # Fake init (Use only for docstring)
    def __init__(self,
                 npy_pdm_gnum_t                         nx,
                 npy_pdm_gnum_t                         ny,
                 npy_pdm_gnum_t                         nz,
                 NPY.double_t                           length,
                 NPY.double_t                           zero_x,
                 NPY.double_t                           zero_y,
                 NPY.double_t                           zero_z,
                 PDM_Mesh_nodal_elt_t                   t_elmt,
                 int                                    order,
                 MPI.Comm                               comm):
      """
      __init__(nx, ny, nz, length, zero_x, zero_y, zero_z, t_elmt, order, comm)

      Initialize a distributed mesh generator structure

      Parameters:
        nx     (npy_pdm_gnum_t)       : Number of vertices on segments in x-direction
        ny     (npy_pdm_gnum_t)       : Number of vertices on segments in y-direction
        nz     (npy_pdm_gnum_t)       : Number of vertices on segments in z-direction
        length (double)               : Segment length
        zero_x (double)               : x-coordinate of the origin
        zero_y (double)               : y-coordinate of the origin
        zero_z (double)               : z-coordinate of the origin
        t_elmt (PDM_Mesh_nodal_elt_t) : Element type
        order  (int)                  : Element order
        comm   (MPI.Comm)             : MPI communicator
      """
    # ------------------------------------------------------------------
    def __cinit__(self,
                  npy_pdm_gnum_t                         nx,
                  npy_pdm_gnum_t                         ny,
                  npy_pdm_gnum_t                         nz,
                  NPY.double_t                           length,
                  NPY.double_t                           zero_x,
                  NPY.double_t                           zero_y,
                  NPY.double_t                           zero_z,
                  PDM_Mesh_nodal_elt_t                   t_elmt,
                  int                                    order,
                  MPI.Comm                               comm):

        """
        """
        # ~> MPI communicator
        cdef MPI.MPI_Comm c_comm = comm.ob_mpi

        # -> Create dcube
        self._dcube = PDM_dcube_nodal_gen_create(PDM_MPI_mpi_2_pdm_mpi_comm (<void *> &c_comm),
                                                 nx,
                                                 ny,
                                                 nz,
                                                 length,
                                                 zero_x,
                                                 zero_y,
                                                 zero_z,
                                                 t_elmt,
                                                 order,
                                                 PDM_OWNERSHIP_USER) # Python takes ownership



    # ------------------------------------------------------------------
    def __dealloc__(self):
        PDM_dcube_nodal_gen_free(self._dcube)
    # ------------------------------------------------------------------
    def set_ordering(self, char *pdm_ho_ordering):
      """
      set_ordering(pdm_ho_ordering)

      Set the HO-ordering

      Parameters:
        pdm_ho_ordering (str) : Name of the HO-ordering
      """
      PDM_dcube_nodal_gen_ordering_set(self._dcube, pdm_ho_ordering)
    # ------------------------------------------------------------------
    def compute(self):
      """
      Generate a distributed nodal mesh
      """
      PDM_dcube_nodal_gen_build(self._dcube)


    # ------------------------------------------------------------------
    # def dcube_dim_get(self):
    #     """
    #        Get dcube dimensions
    #     """
    #     # ************************************************************************
    #     # > Declaration
    #     cdef int n_face_group
    #     cdef int dn_cell
    #     cdef int dn_face
    #     cdef int dn_vtx
    #     cdef int sface_vtx
    #     cdef int sface_group
    #     # ************************************************************************

    #     PDM_dcube_nodal_gen_dim_get(self._dcube,
    #                                &n_face_group,
    #                                &dn_cell,
    #                                &dn_face,
    #                                &dn_vtx,
    #                                &sface_vtx,
    #                                &sface_group)

    #     return {'n_face_group' : n_face_group,
    #             'dn_cell'      : dn_cell,
    #             'dn_face'      : dn_face,
    #             'dn_vtx'       : dn_vtx,
    #             'sface_vtx'    : sface_vtx,
    #             'sface_group'  : sface_group}

    # ------------------------------------------------------------------
    # def dcube_val_get(self):
    #     """
    #        Get dcube data
    #     """
    #     # ************************************************************************
    #     # > Declaration
    #     cdef NPY.npy_intp dim
    #     cdef PDM_g_num_t  *delmt_vtx
    #     cdef double       *dvtx_coord
    #     cdef int          *dface_group_idx
    #     cdef PDM_g_num_t  *dface_group
    #     # ************************************************************************

    #     dims = self.dcube_dim_get()

    #     PDM_dcube_nodal_gen_data_get(self._dcube,
    #                                  &delmt_vtx,
    #                                  &dvtx_coord,
    #                                  &dface_group_idx,
    #                                  &dface_group)

    #     dim = <NPY.npy_intp> 8*dims['dn_cell'] # Car hexa ...
    #     dim = <NPY.npy_intp> 6*dims['dn_cell'] # Car prismes ...
    #     dim = <NPY.npy_intp> 4*dims['dn_cell'] # Car tetra  ...
    #     np_delmt_vtx = create_numpy_g(delmt_vtx, dim)

    #     np_vtx_coord = create_numpy_d(dvtx_coord, 3*dims['dn_vtx'])
    #     # np_dface_group_idx = create_numpy_i(dface_group_idx, dims['n_face_group'] + 1)
    #     # np_dface_group = create_numpy_g(dface_group, dims['sface_group'])

    #     return {'delmt_vtx'       : np_delmt_vtx,
    #             'dvtx_coord'      : np_vtx_coord}
    #     #         'dface_group_idx' : np_dface_group_idx,
    #     #         'dface_group'     : np_dface_group}
    #     # return {'delmt_vtx'       : np_delmt_vtx,
    #     #         'dvtx_coord'      : np_vtx_coord,
    #     #         'dface_group_idx' : np_dface_group_idx,
    #     #         'dface_group'     : np_dface_group}

    # ------------------------------------------------------------------------
    def get_dmesh_nodal(self):
      """
      Get the associated distributed nodal mesh object

      Returns
        Distributed nodal mesh object (:py:class:`DistributedMeshNodalCaspule`)
      """
      # ************************************************************************
      # > Declaration
      cdef PDM_dmesh_nodal_t* dmn
      # ************************************************************************
      dmn = PDM_dcube_nodal_gen_dmesh_nodal_get(self._dcube)

      py_caps = PyCapsule_New(dmn, NULL, NULL);

      return DistributedMeshNodalCapsule(py_caps) # The free is inside the class

    # ------------------------------------------------------------------
    def set_random_factor(self, double random_factor):
      """
      set_random_factor(random_factor)

      Set randomization factor

      Parameters:
        random_factor (double) : Randomization factor (between 0 and 1)
      """
      PDM_dcube_nodal_gen_random_factor_set(self._dcube, random_factor)

# ------------------------------------------------------------------
cdef class DCubeNodalGeneratorCartTopo:
    """
       DCubeNodalGenerator
    """
    # > For Ppart
    cdef PDM_dcube_nodal_t      **_dcube
    cdef PDM_domain_interface_t  *dom_intrf
    cdef int n_dom
    # ------------------------------------------------------------------
    def __cinit__(self,
                  int                                    n_dom_i,
                  int                                    n_dom_j,
                  int                                    n_dom_k,
                  int                                    periodic_i,
                  int                                    periodic_j,
                  int                                    periodic_k,
                  npy_pdm_gnum_t                         nx,
                  npy_pdm_gnum_t                         ny,
                  npy_pdm_gnum_t                         nz,
                  NPY.double_t                           length,
                  NPY.double_t                           zero_x,
                  NPY.double_t                           zero_y,
                  NPY.double_t                           zero_z,
                  PDM_Mesh_nodal_elt_t                   t_elmt,
                  int                                    order,
                  MPI.Comm                               comm):

        """
        """
        # ~> Communicator Mpi
        cdef MPI.MPI_Comm c_comm = comm.ob_mpi

        self.n_dom = n_dom_i * n_dom_j * n_dom_k;

        # -> Create dcube
        PDM_dcube_nodal_cart_topo(PDM_MPI_mpi_2_pdm_mpi_comm (<void *> &c_comm),
                                  n_dom_i,
                                  n_dom_j,
                                  n_dom_k,
                                  periodic_i,
                                  periodic_j,
                                  periodic_k,
                                  nx,
                                  ny,
                                  nz,
                                  length,
                                  0.,
                                  0.,
                                  0.,
                                  t_elmt,
                                  order,
                                  &self._dcube,
                                  &self.dom_intrf,
                                  PDM_OWNERSHIP_USER); # Python takes ownership

    # ------------------------------------------------------------------------
    def get_n_domain(self):
      """
      """
      return self.n_dom

    # ------------------------------------------------------------------------
    def get_dmesh_nodal(self, int i_domain):
      """
      """
      # ************************************************************************
      # > Declaration
      cdef PDM_dmesh_nodal_t* dmn
      # ************************************************************************
      dmn = PDM_dcube_nodal_gen_dmesh_nodal_get(self._dcube[i_domain])

      py_caps = PyCapsule_New(dmn, NULL, NULL);

      return DistributedMeshNodalCapsule(py_caps) # The free is inside the class

    # ------------------------------------------------------------------------
    def domain_interface_get(self):
      """
      """
      # ************************************************************************
      # > Declaration
      cdef PDM_domain_interface_t* dom_intrf
      # ************************************************************************
      py_caps = PyCapsule_New(self.dom_intrf, NULL, NULL);

      return DomInterfaceCapsule(py_caps) # The free is inside the class

    # ------------------------------------------------------------------
    def __dealloc__(self):
      for i in range(self.n_dom):
        PDM_dcube_nodal_gen_free(self._dcube[i])

      free(self._dcube)

