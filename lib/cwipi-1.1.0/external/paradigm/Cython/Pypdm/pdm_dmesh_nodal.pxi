# cython: c_string_type=str, c_string_encoding=ascii
cdef extern from "pdm_mesh_nodal.h":
  int PDM_Mesh_nodal_n_vertices_element(PDM_Mesh_nodal_elt_t type,
                                        int            order)

cdef extern from "pdm_dmesh_nodal.h":
    # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    # > Wrapping of Ppart Structure
    ctypedef struct PDM_dmesh_nodal_t:
      pass
    # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    # > Wrapping of function
    PDM_dmesh_nodal_t* PDM_DMesh_nodal_create(PDM_MPI_Comm comm,
                                              int          mesh_dimension,
                                              PDM_g_num_t  n_vtx,
                                              PDM_g_num_t  n_cell,
                                              PDM_g_num_t  n_face,
                                              PDM_g_num_t  n_edge)
    void PDM_DMesh_nodal_free(PDM_dmesh_nodal_t* dmn)

    void PDM_DMesh_nodal_coord_set(PDM_dmesh_nodal_t* dmn, int n_vtx, double* coords, PDM_ownership_t    owner)
    void PDM_DMesh_nodal_section_g_dims_get(PDM_dmesh_nodal_t* dmn,
                                            PDM_g_num_t        *n_cell_abs,
                                            PDM_g_num_t        *n_face_abs,
                                            PDM_g_num_t        *n_edge_abs,
                                            PDM_g_num_t        *n_vtx_abs)

    # PDM_g_num_t *PDM_DMesh_nodal_distrib_vtx_get(PDM_dmesh_nodal_t* dmn)
    # PDM_g_num_t *PDM_DMesh_nodal_distrib_section_get(PDM_dmesh_nodal_t* dmn, int id_section)

    int                  PDM_DMesh_nodal_n_vtx_get(PDM_dmesh_nodal_t* dmn)
    int*                 PDM_DMesh_nodal_vtx_tag_get(PDM_dmesh_nodal_t* dmn)

    int*                 PDM_DMesh_nodal_sections_id_get(PDM_dmesh_nodal_t* dmn, PDM_geometry_kind_t geom_kind)

    int                  PDM_DMesh_nodal_n_section_get(PDM_dmesh_nodal_t* dmn, PDM_geometry_kind_t geom_kind)

    PDM_Mesh_nodal_elt_t PDM_DMesh_nodal_section_type_get(PDM_dmesh_nodal_t* dmn, PDM_geometry_kind_t geom_kind, int id_section)

    double*              PDM_DMesh_nodal_vtx_get(PDM_dmesh_nodal_t* dmn)
    PDM_g_num_t*         PDM_DMesh_nodal_section_distri_std_get(PDM_dmesh_nodal_t *dmesh_nodal,
                                                                                   PDM_geometry_kind_t  geom_kind,
                                                                                   int                id_section)
    PDM_g_num_t*         PDM_DMesh_nodal_section_distri_std_copy_get(PDM_dmesh_nodal_t   *dmesh_nodal,
                                                                     int                  id_section)

    PDM_g_num_t*         PDM_dmesh_nodal_vtx_distrib_copy_get(PDM_dmesh_nodal_t *dmesh_nodal)
    int                  PDM_DMesh_nodal_section_add(PDM_dmesh_nodal_t* dmn, PDM_geometry_kind_t geom_kind, PDM_Mesh_nodal_elt_t t_elt)

    void                 PDM_DMesh_nodal_update_ownership(PDM_dmesh_nodal_t* dmn, PDM_ownership_t owner)
    void                 PDM_DMesh_nodal_section_std_set(PDM_dmesh_nodal_t* dmn,
                                                                            PDM_geometry_kind_t geom_kind,
                                                                            int                 id_section,
                                                                            int                 n_elmts,
                                                                            PDM_g_num_t*        connec,
                                                                            PDM_ownership_t     owner)

    PDM_g_num_t* PDM_DMesh_nodal_section_std_get(PDM_dmesh_nodal_t* dmn, PDM_geometry_kind_t geom_kind, int id_section)
    int PDM_DMesh_nodal_section_n_elt_get(PDM_dmesh_nodal_t* dmn, int id_section)

    void PDM_DMesh_nodal_section_poly2d_set(PDM_dmesh_nodal_t   *dmesh_nodal,
                                            PDM_geometry_kind_t  geom_kind,
                                            const int            id_section,
                                            const PDM_l_num_t    n_elt,
                                            PDM_l_num_t         *connec_idx,
                                            PDM_g_num_t         *connec,
                                            PDM_ownership_t      owner)
    void PDM_DMesh_nodal_section_group_elmt_set(PDM_dmesh_nodal_t  *dmesh_nodal,
                                                PDM_geometry_kind_t  geom_kind,
                                                int                 n_group_elmt,
                                                int                *dgroup_elmt_idx,
                                                PDM_g_num_t        *dgroup_elmt,
                                                PDM_ownership_t     owner)
    void PDM_DMesh_nodal_section_group_elmt_get(PDM_dmesh_nodal_t    *dmesh_nodal,
                                                                   PDM_geometry_kind_t   geom_kind,
                                                                   int                  *n_group_elmt,
                                                                   int                 **dgroup_elmt_idx,
                                                                   PDM_g_num_t         **dgroup_elmt)

    void PDM_dmesh_nodal_generate_distribution(PDM_dmesh_nodal_t* dmn)

    void PDM_dmesh_nodal_dump_vtk(PDM_dmesh_nodal_t   *dmn,
                                  PDM_geometry_kind_t  geom_kind,
                                  const char          *filename_patter)
    # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

cdef extern from "pdm_elt_parent_find.h":
    # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    # > Wrapping of Ppart Structure
    void PDM_elt_parent_find_from_distrib(PDM_g_num_t  *elt_distrib,
                                          int          *elt_def_idx,
                                          PDM_g_num_t  *elt_def,
                                          PDM_g_num_t  *elt_to_find_distrib,
                                          int          *elt_to_find_def_idx,
                                          PDM_g_num_t  *elt_to_find_def,
                                          PDM_MPI_Comm  comm,
                                          PDM_g_num_t  *parent)

    void PDM_elt_parent_find(int           dnelt,
                             int          *elt_def_idx,
                             PDM_g_num_t  *elt_def,
                             int           dnelt_to_find,
                             int          *elt_to_find_def_idx,
                             PDM_g_num_t  *elt_to_find_def,
                             PDM_MPI_Comm  comm,
                             PDM_g_num_t  *parent)

cdef extern from "pdm_distrib.h":

    void PDM_distrib_compute(int           dnelt,
                             PDM_g_num_t  *elt_distrib,
                             int           offset,
                             PDM_MPI_Comm  comm)

# :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

# ------------------------------------------------------------------
cdef class DistributedMeshNodal:
    """
       DistributedMeshNodal: Interface to build face from Element->Vtx connectivity
    """
    # ************************************************************************
    # > Class attributes
    cdef PDM_dmesh_nodal_t *dmn
    keep_alive = list()
    # cdef int idmesh
    cdef int n_rank
    # ************************************************************************
    # ------------------------------------------------------------------------
    def __cinit__(self, MPI.Comm    comm,
                        PDM_g_num_t n_vtx,
                        PDM_g_num_t n_cell,
                        PDM_g_num_t n_face = -1,
                        PDM_g_num_t n_edge = -1,
                        int         mesh_dimension = 3):
        """
        TODOUX
        """
        # ************************************************************************
        # > Declaration
        # cdef int      nElts
        # cdef int      idx
        # # > Numpy array
        # cdef NPY.ndarray[npy_pdm_gnum_t, ndim=1, mode='fortran'] partLNToGN
        # ************************************************************************

        # ::::::::::::::::::::::::::::::::::::::::::::::::::
        self.n_rank = comm.Get_size()
        # ::::::::::::::::::::::::::::::::::::::::::::::::::

        # ::::::::::::::::::::::::::::::::::::::::::::::::::
        # > Convert mpi4py -> PDM_MPI
        cdef MPI.MPI_Comm c_comm = comm.ob_mpi
        cdef PDM_MPI_Comm PDMC   = PDM_MPI_mpi_2_pdm_mpi_comm(<void *> &c_comm)
        # ::::::::::::::::::::::::::::::::::::::::::::::::::

        # ::::::::::::::::::::::::::::::::::::::::::::::::::
        self.dmn = PDM_DMesh_nodal_create(PDMC, mesh_dimension, n_vtx, n_cell, n_face, n_edge)
        # ::::::::::::::::::::::::::::::::::::::::::::::::::

    # ------------------------------------------------------------------------
    def set_coordinates(self, NPY.ndarray[NPY.double_t  , mode='c', ndim=1] dvtx_coord):
        """
        """
        # ************************************************************************
        # > Declaration
        cdef int n_vtx
        # ************************************************************************

        # ::::::::::::::::::::::::::::::::::::::::::::::::::
        self.keep_alive.append(dvtx_coord)
        # ::::::::::::::::::::::::::::::::::::::::::::::::::

        # ::::::::::::::::::::::::::::::::::::::::::::::::::
        n_vtx = dvtx_coord.shape[0]//3
        PDM_DMesh_nodal_coord_set(self.dmn, n_vtx, <double *> dvtx_coord.data, PDM_OWNERSHIP_USER)
        # ::::::::::::::::::::::::::::::::::::::::::::::::::

    # ------------------------------------------------------------------------
    def set_sections(self,
                     PDM_geometry_kind_t geom_kind,
                     list elmt_list,
                     NPY.ndarray[NPY.int32_t, mode='c', ndim=1] elmts_type,
                     NPY.ndarray[NPY.int32_t, mode='c', ndim=1] n_elemts):
        """
           TODO : Split function as PDM
        """
        # ************************************************************************
        # > Declaration
        cdef int n_vtx
        cdef NPY.ndarray[npy_pdm_gnum_t, ndim=1, mode='fortran'] connect
        # ************************************************************************

        # ::::::::::::::::::::::::::::::::::::::::::::::::::
        # > Panic assert
        assert(len(elmt_list) == elmts_type.shape[0])
        # ::::::::::::::::::::::::::::::::::::::::::::::::::

        # ::::::::::::::::::::::::::::::::::::::::::::::::::
        self.keep_alive.append(elmt_list)
        # ::::::::::::::::::::::::::::::::::::::::::::::::::

        # ::::::::::::::::::::::::::::::::::::::::::::::::::
        for i_elmt, connect in enumerate(elmt_list):
          id_section = PDM_DMesh_nodal_section_add(self.dmn,
                                                   geom_kind,
                            <PDM_Mesh_nodal_elt_t> elmts_type[i_elmt])
          PDM_DMesh_nodal_section_std_set(self.dmn,
                                          geom_kind,
                                          id_section,
                                          n_elemts[i_elmt],
                                          <PDM_g_num_t *> connect.data,
                                          PDM_OWNERSHIP_USER)
        # ::::::::::::::::::::::::::::::::::::::::::::::::::

    def set_poly2d_section(self, 
                           NPY.ndarray[NPY.int32_t, mode='c', ndim=1] poly_connectivity_idx, 
                           NPY.ndarray[npy_pdm_gnum_t, mode='c', ndim=1] poly_connectivity):
        id_section = PDM_DMesh_nodal_section_add(self.dmn, _PDM_GEOMETRY_KIND_SURFACIC, _PDM_MESH_NODAL_POLY_2D)

        # ::::::::::::::::::::::::::::::::::::::::::::::::::
        self.keep_alive.append(poly_connectivity_idx)
        self.keep_alive.append(poly_connectivity)
        # ::::::::::::::::::::::::::::::::::::::::::::::::::

        PDM_DMesh_nodal_section_poly2d_set(self.dmn,
                                           _PDM_GEOMETRY_KIND_SURFACIC,
                                           id_section,
                                           poly_connectivity_idx.size-1, # n_elts local
                          <PDM_l_num_t*>   poly_connectivity_idx.data,
                          <PDM_g_num_t*>   poly_connectivity.data,
                                           PDM_OWNERSHIP_USER)

    # ------------------------------------------------------------------------
    def set_group_elmt(self,
                       PDM_geometry_kind_t geom_kind,
                       n_group_elmt,
                       NPY.ndarray[NPY.int32_t   , mode='c', ndim=1] dgroup_elmt_idx,
                       NPY.ndarray[npy_pdm_gnum_t, mode='c', ndim=1] dgroup_elmt):
        """
           TODO : Split function as PDM
        """
        if(dgroup_elmt_idx is None):
          PDM_DMesh_nodal_section_group_elmt_set(self.dmn,
                                                 geom_kind,
                                                 n_group_elmt,
                                                 NULL, NULL, PDM_OWNERSHIP_USER)
        else:
          # ::::::::::::::::::::::::::::::::::::::::::::::::::
          self.keep_alive.append(dgroup_elmt_idx)
          self.keep_alive.append(dgroup_elmt)
          # ::::::::::::::::::::::::::::::::::::::::::::::::::
          PDM_DMesh_nodal_section_group_elmt_set(self.dmn,
                                                 geom_kind,
                                                 n_group_elmt,
                                          <int*> dgroup_elmt_idx.data,
                                  <PDM_g_num_t*> dgroup_elmt.data,
                                                 PDM_OWNERSHIP_USER)

    # ------------------------------------------------------------------------
    def generate_distribution(self):
        """
        """
        # ************************************************************************
        # > Declaration
        # ************************************************************************
        PDM_dmesh_nodal_generate_distribution(self.dmn)

    # ------------------------------------------------------------------------
    def __dealloc__(self):
      """
         Use the free method of PDM Lib
      """
      # ************************************************************************
      # > Declaration
      # ************************************************************************
      # > Free Ppart Structure
      PDM_DMesh_nodal_free(self.dmn)

# ------------------------------------------------------------------
cdef class DistributedMeshNodalCapsule:
  """
  """
  # ************************************************************************
  # > Class attributes
  cdef PDM_dmesh_nodal_t* dmn
  # ************************************************************************
  # ------------------------------------------------------------------------
  def __cinit__(self, object caps):
    """
    """
    # print("DistributedMeshNodalCapsule", PyCapsule_GetName(caps))
    cdef PDM_dmesh_nodal_t* caps_dmn = <PDM_dmesh_nodal_t *> PyCapsule_GetPointer(caps, NULL)
    self.dmn = caps_dmn;

  # ------------------------------------------------------------------------
  def dmesh_nodal_get_sections(self, PDM_geometry_kind_t geom_kind, MPI.Comm comm):
    """
    """
    return dmesh_nodal_get_sections(self, geom_kind, comm)

  # ------------------------------------------------------------------------
  def dmesh_nodal_get_vtx(self, MPI.Comm comm):
    """
    """
    return dmesh_nodal_get_vtx(self, comm)

  # ------------------------------------------------------------------------
  def dmesh_nodal_get_group(self, PDM_geometry_kind_t geom_kind):
    """
    """
    return dmesh_nodal_get_group(self, geom_kind)

  # ------------------------------------------------------------------------
  def dmesh_nodal_get_g_dims(self):
    """
    """
    # print("Wrap dmesh_nodal_get_g_dims")
    return dmesh_nodal_get_g_dims(self)

  # ------------------------------------------------------------------------
  def dump_vtk(self,
               PDM_geometry_kind_t  geom_kind,
               char                *filename_pattern):
    """
    dump_vtk(geom_kind, filename_pattern)

    Export in VTK format (ASCII)

    .. note::
      Each rank dumps a file for each nodal section

    Parameters:
      geom_kind        (PDM_geometry_kind_t) : Geometry kind to export
      filename_pattern (str)                 : File name pattern
    """
    PDM_dmesh_nodal_dump_vtk(self.dmn,
                             geom_kind,
                             filename_pattern)

  # ------------------------------------------------------------------------
  def __dealloc__(self):
    """
       Use the free method of PDM Lib
    """
    # print("DistributedMeshNodalCapsule::__dealloc__")
    PDM_DMesh_nodal_free(self.dmn)
    # print("DistributedMeshNodalCapsule::__dealloc__ end z")

ctypedef fused DMeshNodal:
  DistributedMeshNodal
  DistributedMeshNodalCapsule

def generate_distribution(DMeshNodal pydmn):
  """
  """
  PDM_dmesh_nodal_generate_distribution(pydmn.dmn)

def dmesh_nodal_get_g_dims(DMeshNodal pydmn):
  """
  """
  # ************************************************************************
  # > Declaration
  cdef PDM_g_num_t n_cell_abs, n_face_abs, n_edge_abs, n_vtx_abs
  # ************************************************************************

  PDM_DMesh_nodal_section_g_dims_get(pydmn.dmn, &n_cell_abs, &n_face_abs, &n_edge_abs, &n_vtx_abs)
  return {"n_cell_abs" : n_cell_abs,
          "n_face_abs" : n_face_abs,
          "n_edge_abs" : n_edge_abs,
          "n_vtx_abs"  : n_vtx_abs}

def dmesh_nodal_get_vtx(DMeshNodal pydmn, MPI.Comm    comm):
  """
  """
  # ************************************************************************
  # > Declaration
  cdef double               *vtx_coord
  cdef int                  *vtx_tag
  cdef NPY.npy_intp          dim
  # ************************************************************************

  PDM_DMesh_nodal_update_ownership(pydmn.dmn, PDM_OWNERSHIP_USER)
  vtx_distrib = PDM_dmesh_nodal_vtx_distrib_copy_get(pydmn.dmn)
  n_vtx = PDM_DMesh_nodal_n_vtx_get(pydmn.dmn);
  vtx_coord = PDM_DMesh_nodal_vtx_get(pydmn.dmn)
  vtx_tag = PDM_DMesh_nodal_vtx_tag_get(pydmn.dmn)

  return {"np_vtx"         : create_numpy_d(vtx_coord,   3*n_vtx,           False),
          "np_vtx_distrib" : create_numpy_g(vtx_distrib, comm.Get_size()+1, False),
          "np_vtx_tag"     : create_numpy_i(vtx_tag,     n_vtx,             False)}

def dmesh_nodal_get_sections(DMeshNodal          pydmn,
                             PDM_geometry_kind_t geom_kind,
                             MPI.Comm            comm):
  """
  """
  # ************************************************************************
  # > Declaration
  cdef int                   n_section
  cdef int                   n_vtx_per_elmt
  cdef PDM_g_num_t           dn_elmt
  cdef int                  *section_id
  cdef int                  *connect_idx
  cdef PDM_g_num_t          *connect
  cdef PDM_g_num_t          *section_distrib
  cdef PDM_Mesh_nodal_elt_t  t_elmt
  cdef NPY.npy_intp          dim
  # ************************************************************************

  PDM_DMesh_nodal_update_ownership(pydmn.dmn, PDM_OWNERSHIP_USER)
  n_section  = PDM_DMesh_nodal_n_section_get(pydmn.dmn, geom_kind)
  section_id = PDM_DMesh_nodal_sections_id_get(pydmn.dmn, geom_kind)

  # print("n_section : ", n_section)

  sections = []
  for i_section in range(n_section):
    id_section = section_id[i_section]
    t_elmt = PDM_DMesh_nodal_section_type_get(pydmn.dmn, geom_kind, id_section)
    # > For now only use this interfaces for standard element ... (to be refactor with HO and polyhedra elements)
    assert(t_elmt != PDM_MESH_NODAL_POLY_2D)
    assert(t_elmt != PDM_MESH_NODAL_POLY_3D)

    section_distrib = PDM_DMesh_nodal_section_distri_std_get(pydmn.dmn, geom_kind, id_section)
    connect         = PDM_DMesh_nodal_section_std_get(pydmn.dmn, geom_kind, id_section)

    # > Build numpy capsule
    np_distrib_tmp = create_numpy_g(section_distrib, comm.Get_size()+1, flag_owndata=False)
    np_distrib = NPY.copy(np_distrib_tmp)

    # > Build numpy capsule
    dn_elmt = np_distrib[comm.Get_rank()+1] - np_distrib[comm.Get_rank()]
    n_vtx_per_elmt = PDM_Mesh_nodal_n_vertices_element(t_elmt, 1)
    np_connec = create_numpy_g(connect, n_vtx_per_elmt*dn_elmt)

    sections.append({"pdm_type"   : t_elmt,
                     "np_distrib" : np_distrib,
                     "np_connec"  : np_connec})
    # print("t_elmt : ", t_elmt )
    # print("np_distrib : ", np_distrib )
    # print("np_connec : ", np_connec )

  # print("Return dmesh_nodal_get_sections")
  return {"sections" : sections}

def dmesh_nodal_get_group(DMeshNodal pydmn, PDM_geometry_kind_t geom_kind):
  """
  """
  # ************************************************************************
  # > Declaration
  cdef int                   n_group
  cdef int                  *dgroup_elmt_idx
  cdef PDM_g_num_t          *dgroup_elmt
  cdef NPY.npy_intp          dim
  # ************************************************************************

  PDM_DMesh_nodal_update_ownership(pydmn.dmn, PDM_OWNERSHIP_USER)
  PDM_DMesh_nodal_section_group_elmt_get(pydmn.dmn, geom_kind, &n_group, &dgroup_elmt_idx, &dgroup_elmt);

  if n_group == 0:
    return None

  np_dgroup_elmt_idx = create_numpy_i(dgroup_elmt_idx, n_group+1)
  np_dgroup_elmt = create_numpy_g(dgroup_elmt, np_dgroup_elmt_idx[n_group])

  return {"dgroup_elmt_idx" : np_dgroup_elmt_idx,
          "dgroup_elmt"     : np_dgroup_elmt}


# ------------------------------------------------------------------------
def ElementParentFind(int                                           dnelt,
                      NPY.ndarray[NPY.int32_t   , mode='c', ndim=1] elt_def_idx,
                      NPY.ndarray[npy_pdm_gnum_t, mode='c', ndim=1] elt_def,
                      int                                           dnelt_to_find,
                      NPY.ndarray[NPY.int32_t   , mode='c', ndim=1] elt_to_find_def_idx,
                      NPY.ndarray[npy_pdm_gnum_t, mode='c', ndim=1] elt_to_find_def,
                      MPI.Comm    comm,
                      NPY.ndarray[npy_pdm_gnum_t, mode='c', ndim=1] parent):
    """
    """
    # ************************************************************************
    # > Declaration
    # ************************************************************************

    # ::::::::::::::::::::::::::::::::::::::::::::::::::
    # > Convert mpi4py -> PDM_MPI
    cdef MPI.MPI_Comm c_comm = comm.ob_mpi
    cdef PDM_MPI_Comm PDMC   = PDM_MPI_mpi_2_pdm_mpi_comm(<void *> &c_comm)
    # ::::::::::::::::::::::::::::::::::::::::::::::::::

    # ::::::::::::::::::::::::::::::::::::::::::::::::::
    PDM_elt_parent_find(dnelt,
                        <int *>         elt_def_idx.data,
                        <PDM_g_num_t *> elt_def.data,
                        dnelt_to_find,
                        <int *>         elt_to_find_def_idx.data,
                        <PDM_g_num_t *> elt_to_find_def.data,
                        PDMC,
                        <PDM_g_num_t *> parent.data)
    # ::::::::::::::::::::::::::::::::::::::::::::::::::


# ------------------------------------------------------------------------
def ComputeDistributionFromDelmt(int         dnelt,
                                 MPI.Comm    comm,
                                 int         offset=0):
    """
    """
    # ************************************************************************
    # > Declaration
    cdef NPY.ndarray[npy_pdm_gnum_t, ndim=1, mode='c'] elt_distrib = NPY.empty( comm.Get_size() + 1, dtype=npy_pdm_gnum_dtype, order='C')
    # ************************************************************************

    # ::::::::::::::::::::::::::::::::::::::::::::::::::
    # > Convert mpi4py -> PDM_MPI
    cdef MPI.MPI_Comm c_comm = comm.ob_mpi
    cdef PDM_MPI_Comm PDMC   = PDM_MPI_mpi_2_pdm_mpi_comm(<void *> &c_comm)
    # ::::::::::::::::::::::::::::::::::::::::::::::::::

    # ::::::::::::::::::::::::::::::::::::::::::::::::::
    PDM_distrib_compute(dnelt,
                        <PDM_g_num_t *> elt_distrib.data,
                        offset,
                        PDMC)
    # ::::::::::::::::::::::::::::::::::::::::::::::::::

    return elt_distrib

# ------------------------------------------------------------------------


