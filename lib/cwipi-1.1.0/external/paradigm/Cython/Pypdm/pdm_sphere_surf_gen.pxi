cdef extern from "pdm_sphere_surf_gen.h":
    void PDM_sphere_surf_icosphere_gen(const PDM_MPI_Comm   comm,
                                       const PDM_g_num_t    n,
                                       const double         x_center,
                                       const double         y_center,
                                       const double         z_center,
                                       const double         radius,
                                             double       **dvtx_coord,
                                             int          **dface_vtx_idx,
                                             PDM_g_num_t  **dface_vtx,
                                             PDM_g_num_t  **distrib_vtx,
                                             PDM_g_num_t  **distrib_face)

    void PDM_sphere_surf_icosphere_gen_nodal(const PDM_MPI_Comm        comm,
                                             const PDM_g_num_t         n,
                                             const double              x_center,
                                             const double              y_center,
                                             const double              z_center,
                                             const double              radius,
                                             PDM_dmesh_nodal_t **_dmn)

    void PDM_sphere_surf_icosphere_gen_part(const PDM_MPI_Comm        comm,
                                            const PDM_g_num_t         n,
                                            const double              x_center,
                                            const double              y_center,
                                            const double              z_center,
                                            const double              radius,
                                            const int                 n_part,
                                            const PDM_split_dual_t    part_method,
                                                  int               **pn_vtx,
                                                  double           ***pvtx_coord,
                                                  PDM_g_num_t      ***pvtx_ln_to_gn,
                                                  int               **pn_face,
                                                  int              ***pface_vtx_idx,
                                                  int              ***pface_vtx,
                                                  PDM_g_num_t      ***pface_ln_to_gn)

    void PDM_sphere_surf_gen(const PDM_MPI_Comm        comm,
                             const PDM_g_num_t         nu,
                             const PDM_g_num_t         nv,
                             const double              x_center,
                             const double              y_center,
                             const double              z_center,
                             const double              radius,
                                   double            **dvtx_coord,
                                   int               **dface_vtx_idx,
                                   PDM_g_num_t       **dface_vtx,
                                   PDM_g_num_t       **distrib_vtx,
                                   PDM_g_num_t       **distrib_face)


    void PDM_sphere_surf_gen_nodal(const PDM_MPI_Comm        comm,
                                   const PDM_g_num_t         nu,
                                   const PDM_g_num_t         nv,
                                   const double              x_center,
                                   const double              y_center,
                                   const double              z_center,
                                   const double              radius,
                                         PDM_dmesh_nodal_t **dmn)
# ------------------------------------------------------------------------

def sphere_surf_icosphere_gen(MPI.Comm       comm,
                              npy_pdm_gnum_t n,
                              NPY.double_t   x_center,
                              NPY.double_t   y_center,
                              NPY.double_t   z_center,
                              NPY.double_t   radius):
  """
  sphere_surf_icosphere_gen(comm, n, x_center, y_center, z_center, radius)

  Create a sphere surface mesh (icosphere)

  :note: The output mesh is *block-distributed*

  Parameters:
    comm     (MPI.Comm)       : MPI communicator
    n        (npy_pdm_gnum_t) : Icosphere subdivision level
    x_center (double)         : x-coordinate of the sphere center
    y_center (double)         : y-coordinate of the sphere center
    z_center (double)         : z-coordinate of the sphere center
    radius   (double)         : Sphere radius

  Returns:
    Dictionary
      - ``"dvtx_coord"``    : Vertex coordinates (np.ndarray[np.double_t])
      - ``"dface_vtx_idx"`` : Index for face -> vertex connectivity (np.ndarray[np.int32_t])
      - ``"dface_vtx"``     : Face -> vertex connectivity (global ids) (np.ndarray[npy_pdm_gnum_t])
      - ``"distrib_vtx"``   : Vertex distribution (size : *n_rank* + 1) (np.ndarray[npy_pdm_gnum_t])
      - ``"distrib_face"``  : Face distribution (size : *n_rank* + 1) (np.ndarray[npy_pdm_gnum_t])
  """
  i_rank = comm.rank
  n_rank = comm.size

  cdef double      *dvtx_coord    = NULL
  cdef int         *dface_vtx_idx = NULL
  cdef PDM_g_num_t *dface_vtx     = NULL
  cdef PDM_g_num_t *distrib_vtx   = NULL
  cdef PDM_g_num_t *distrib_face  = NULL

  cdef MPI.MPI_Comm c_comm = comm.ob_mpi
  PDM_sphere_surf_icosphere_gen(PDM_MPI_mpi_2_pdm_mpi_comm (<void *> &c_comm),
                                n,
                                x_center,
                                y_center,
                                z_center,
                                radius,
                                &dvtx_coord,
                                &dface_vtx_idx,
                                &dface_vtx,
                                &distrib_vtx,
                                &distrib_face)

  dn_vtx  = distrib_vtx [i_rank+1] - distrib_vtx [i_rank]
  dn_face = distrib_face[i_rank+1] - distrib_face[i_rank]

  return {
  "dvtx_coord"    : create_numpy_d(dvtx_coord,    3*dn_vtx),
  "dface_vtx_idx" : create_numpy_i(dface_vtx_idx, dn_face+1),
  "dface_vtx"     : create_numpy_g(dface_vtx,     dface_vtx_idx[dn_face]),
  "distrib_vtx"   : create_numpy_g(distrib_vtx,   n_rank+1),
  "distrib_face"  : create_numpy_g(distrib_face,  n_rank+1),
  }

# ------------------------------------------------------------------------

def sphere_surf_icosphere_gen_nodal(MPI.Comm       comm,
                                    npy_pdm_gnum_t n,
                                    NPY.double_t   x_center,
                                    NPY.double_t   y_center,
                                    NPY.double_t   z_center,
                                    NPY.double_t   radius):
  """
  sphere_surf_icosphere_gen_nodal(comm, n, x_center, y_center, z_center, radius)

  Create a sphere surface mesh (icosphere)

  :note: The output mesh is *block-distributed* in the form of a :py:class:`DistributedMeshNodalCapsule` object

  Parameters:
    comm     (MPI.Comm)       : MPI communicator
    n        (npy_pdm_gnum_t) : Icosphere subdivision level
    x_center (double)         : x-coordinate of the sphere center
    y_center (double)         : y-coordinate of the sphere center
    z_center (double)         : z-coordinate of the sphere center
    radius   (double)         : Sphere radius

  Returns:
    Distributed nodal icosphere mesh (:py:class:`DistributedMeshNodalCapsule`)
  """
  cdef PDM_dmesh_nodal_t *dmn = NULL
  cdef MPI.MPI_Comm c_comm = comm.ob_mpi

  PDM_sphere_surf_icosphere_gen_nodal(PDM_MPI_mpi_2_pdm_mpi_comm (<void *> &c_comm),
                                      n,
                                      x_center,
                                      y_center,
                                      z_center,
                                      radius,
                                      &dmn)

  py_caps = PyCapsule_New(dmn, NULL, NULL)

  return DistributedMeshNodalCapsule(py_caps)
# ------------------------------------------------------------------------

def sphere_surf_icosphere_gen_part(MPI.Comm         comm,
                                   npy_pdm_gnum_t   n,
                                   NPY.double_t     x_center,
                                   NPY.double_t     y_center,
                                   NPY.double_t     z_center,
                                   NPY.double_t     radius,
                                   int              n_part,
                                   PDM_split_dual_t part_method):
  """
  sphere_surf_icosphere_gen_part(comm, n, x_center, y_center, z_center, radius, n_part, part_method)

  Create a sphere surface mesh (icosphere)

  :note: The output mesh is *partitioned*

  Parameters:
    comm        (MPI.Comm)         : MPI communicator
    n           (npy_pdm_gnum_t)   : Icosphere subdivision level
    x_center    (double)           : x-coordinate of the sphere center
    y_center    (double)           : y-coordinate of the sphere center
    z_center    (double)           : z-coordinate of the sphere center
    radius      (double)           : Sphere radius
    n_part      (int)              : Number of partitions
    part_method (PDM_split_dual_t) : Partitioning method

  Returns:
    Dictionary
      - ``"pn_vtx"``         : Number of vertices (np.ndarray[np.int32_t])
      - ``"pvtx_coord"``     : Vertex coordinates (np.ndarray[np.double_t])
      - ``"pvtx_ln_to_gn"``  : Vertex global ids  (np.ndarray[npy_pdm_gnum_t])
      - ``"pn_face"``        : Number of faces    (np.ndarray[np.int32_t])
      - ``"pface_vtx_idx"``  : Index for face -> vertex connectivity (np.ndarray[np.int32_t])
      - ``"pface_vtx"``      : Face -> vertex connectivity  (np.ndarray[np.int32_t])
      - ``"pface_ln_to_gn"`` : Face global ids  (np.ndarray[npy_pdm_gnum_t])
  """
  cdef int          *_pn_vtx         = NULL
  cdef double      **_pvtx_coord     = NULL
  cdef PDM_g_num_t **_pvtx_ln_to_gn  = NULL
  cdef int          *_pn_face        = NULL
  cdef int         **_pface_vtx_idx  = NULL
  cdef int         **_pface_vtx      = NULL
  cdef PDM_g_num_t **_pface_ln_to_gn = NULL

  cdef MPI.MPI_Comm c_comm = comm.ob_mpi
  PDM_sphere_surf_icosphere_gen_part(PDM_MPI_mpi_2_pdm_mpi_comm (<void *> &c_comm),
                                     n,
                                     x_center,
                                     y_center,
                                     z_center,
                                     radius,
                                     n_part,
                                     part_method,
                                     &_pn_vtx,
                                     &_pvtx_coord,
                                     &_pvtx_ln_to_gn,
                                     &_pn_face,
                                     &_pface_vtx_idx,
                                     &_pface_vtx,
                                     &_pface_ln_to_gn)

  pn_vtx  = create_numpy_i(_pn_vtx,  n_part)
  pn_face = create_numpy_i(_pn_face, n_part)

  pvtx_coord     = []
  pvtx_ln_to_gn  = []
  pface_vtx_idx  = []
  pface_vtx      = []
  pface_ln_to_gn = []
  for i_part in range(n_part):
    _n_vtx  = pn_vtx[i_part]
    _n_face = pn_face[i_part]
    pvtx_coord    .append(create_numpy_d (_pvtx_coord[i_part],     3*_n_vtx))
    pvtx_ln_to_gn .append(create_numpy_g (_pvtx_ln_to_gn[i_part],  _n_vtx))
    pface_vtx_idx .append(create_numpy_i (_pface_vtx_idx[i_part],  _n_face+1))
    pface_vtx     .append(create_numpy_i (_pface_vtx[i_part],      pface_vtx_idx[i_part][_n_face]))
    pface_ln_to_gn.append(create_numpy_g (_pface_ln_to_gn[i_part], _n_face))

  return {
  "pn_vtx"         : pn_vtx,
  "pvtx_coord"     : pvtx_coord,
  "pvtx_ln_to_gn"  : pvtx_ln_to_gn,
  "pn_face"        : pn_face,
  "pface_vtx_idx"  : pface_vtx_idx,
  "pface_vtx"      : pface_vtx,
  "pface_ln_to_gn" : pface_ln_to_gn
  }


# ------------------------------------------------------------------------

def sphere_surf_gen(MPI.Comm       comm,
                    npy_pdm_gnum_t nu,
                    npy_pdm_gnum_t nv,
                    NPY.double_t   x_center,
                    NPY.double_t   y_center,
                    NPY.double_t   z_center,
                    NPY.double_t   radius):

  i_rank = comm.rank
  n_rank = comm.size

  cdef double      *dvtx_coord    = NULL
  cdef int         *dface_vtx_idx = NULL
  cdef PDM_g_num_t *dface_vtx     = NULL
  cdef PDM_g_num_t *distrib_vtx   = NULL
  cdef PDM_g_num_t *distrib_face  = NULL

  cdef MPI.MPI_Comm c_comm = comm.ob_mpi
  PDM_sphere_surf_gen(PDM_MPI_mpi_2_pdm_mpi_comm (<void *> &c_comm),
                      nu,
                      nv,
                      x_center,
                      y_center,
                      z_center,
                      radius,
                      &dvtx_coord,
                      &dface_vtx_idx,
                      &dface_vtx,
                      &distrib_vtx,
                      &distrib_face)

  dn_vtx  = distrib_vtx [i_rank+1] - distrib_vtx [i_rank]
  dn_face = distrib_face[i_rank+1] - distrib_face[i_rank]

  return {
  "dvtx_coord"    : create_numpy_d(dvtx_coord,    3*dn_vtx),
  "dface_vtx_idx" : create_numpy_i(dface_vtx_idx, dn_face+1),
  "dface_vtx"     : create_numpy_g(dface_vtx,     dface_vtx_idx[dn_face]),
  "distrib_vtx"   : create_numpy_g(distrib_vtx,   n_rank+1),
  "distrib_face"  : create_numpy_g(distrib_face,  n_rank+1),
  }

# ------------------------------------------------------------------------

def sphere_surf_gen_nodal(MPI.Comm       comm,
                          npy_pdm_gnum_t nu,
                          npy_pdm_gnum_t nv,
                          NPY.double_t   x_center,
                          NPY.double_t   y_center,
                          NPY.double_t   z_center,
                          NPY.double_t   radius):

  cdef PDM_dmesh_nodal_t *dmn = NULL
  cdef MPI.MPI_Comm c_comm = comm.ob_mpi

  PDM_sphere_surf_gen_nodal(PDM_MPI_mpi_2_pdm_mpi_comm (<void *> &c_comm),
                            nu,
                            nv,
                            x_center,
                            y_center,
                            z_center,
                            radius,
                            &dmn)

  py_caps = PyCapsule_New(dmn, NULL, NULL)

  return DistributedMeshNodalCapsule(py_caps)
