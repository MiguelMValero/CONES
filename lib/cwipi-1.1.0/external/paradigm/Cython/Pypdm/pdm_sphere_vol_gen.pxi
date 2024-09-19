cdef extern from "pdm_sphere_vol_gen.h":
    void PDM_sphere_vol_icosphere_gen_nodal(const PDM_MPI_Comm        comm,
                                             const PDM_g_num_t         n,
                                             const double              x_center,
                                             const double              y_center,
                                             const double              z_center,
                                             const double              radius,
                                             PDM_dmesh_nodal_t **_dmn)

    void PDM_sphere_vol_hollow_gen_nodal(const PDM_MPI_Comm        comm,
                                         const PDM_g_num_t         n,
                                         const PDM_g_num_t         n_layer,
                                         const double              x_center,
                                         const double              y_center,
                                         const double              z_center,
                                         const double              radius_interior,
                                         const double              radius_exterior,
                                         const double              geometric_ratio,
                                         PDM_dmesh_nodal_t       **dmn);

# ------------------------------------------------------------------------
def sphere_vol_icosphere_gen_nodal(MPI.Comm        comm,
                                    npy_pdm_gnum_t n,
                                    NPY.double_t   x_center,
                                    NPY.double_t   y_center,
                                    NPY.double_t   z_center,
                                    NPY.double_t   radius):

  cdef PDM_dmesh_nodal_t *dmn = NULL
  cdef MPI.MPI_Comm c_comm = comm.ob_mpi

  PDM_sphere_vol_icosphere_gen_nodal(PDM_MPI_mpi_2_pdm_mpi_comm (<void *> &c_comm),
                                      n,
                                      x_center,
                                      y_center,
                                      z_center,
                                      radius,
                                      &dmn)

  py_caps = PyCapsule_New(dmn, NULL, NULL)

  return DistributedMeshNodalCapsule(py_caps)


# ------------------------------------------------------------------------
def sphere_vol_hollow_gen_nodal(MPI.Comm        comm,
                                npy_pdm_gnum_t  n,
                                npy_pdm_gnum_t  n_layer,
                                NPY.double_t    x_center,
                                NPY.double_t    y_center,
                                NPY.double_t    z_center,
                                NPY.double_t    radius_int,
                                NPY.double_t    radius_ext,
                                NPY.double_t    geometric_ratio):

  cdef PDM_dmesh_nodal_t *dmn = NULL
  cdef MPI.MPI_Comm c_comm = comm.ob_mpi

  PDM_sphere_vol_hollow_gen_nodal(PDM_MPI_mpi_2_pdm_mpi_comm (<void *> &c_comm),
                                      n,
                                      n_layer,
                                      x_center,
                                      y_center,
                                      z_center,
                                      radius_int,
                                      radius_ext,
                                      geometric_ratio,
                                      &dmn)

  py_caps = PyCapsule_New(dmn, NULL, NULL)

  return DistributedMeshNodalCapsule(py_caps)
