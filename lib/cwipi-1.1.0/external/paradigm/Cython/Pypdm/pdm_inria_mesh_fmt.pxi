
cdef extern from "pdm_reader_gamma.h":
    void PDM_write_meshb(char   *filename,
                         int     n_vtx,
                         int     n_tetra,
                         int     n_tri,
                         int     n_edge,
                         double *vtx_coords,
                         int    *vtx_tags,
                         int    *tetra_vtx,
                         int    *tetra_tag,
                         int    *tria_vtx,
                         int    *tria_tag,
                         int    *edge_vtx,
                         int    *edge_tag)

    PDM_dmesh_nodal_t *PDM_reader_gamma_dmesh_nodal(PDM_MPI_Comm   comm,
                                                    char    *filename,
                                                    int            fix_orientation_2d,
                                                    int            fix_orientation_3d)

    void PDM_write_gamma_sol(char   *filename,
                             int     n_vtx,
                             int     n_field,
                             double    *fields)
    void PDM_read_gamma_sol(char   *filename,
                             int     n_vtx,
                             int     n_field,
                             double    *fields)

    void PDM_write_gamma_matsym(char   *filename,
                                int     n_vtx,
                                double *fields)

def write_meshb(char *filename,
                int n_vtx,
                int n_tetra,
                int n_tri,
                int n_edge,
                NPY.ndarray[NPY.double_t, mode='c', ndim=1] vtx_coords,
                NPY.ndarray[NPY.int32_t , mode='c', ndim=1] vtx_tags,
                NPY.ndarray[NPY.int32_t , mode='c', ndim=1] tetra_vtx,
                NPY.ndarray[NPY.int32_t , mode='c', ndim=1] tetra_tag,
                NPY.ndarray[NPY.int32_t , mode='c', ndim=1] tria_vtx,
                NPY.ndarray[NPY.int32_t , mode='c', ndim=1] tria_tag,
                NPY.ndarray[NPY.int32_t , mode='c', ndim=1] edge_vtx,
                NPY.ndarray[NPY.int32_t , mode='c', ndim=1] edge_tag):
  """
  """
  PDM_write_meshb(filename,
                  n_vtx,
                  n_tetra,
                  n_tri,
                  n_edge,
                  <double *> vtx_coords.data,
                  <int    *> vtx_tags.data,
                  <int    *> tetra_vtx.data,
                  <int    *> tetra_tag.data,
                  <int    *> tria_vtx.data,
                  <int    *> tria_tag.data,
                  <int    *> edge_vtx.data,
                  <int    *> edge_tag.data)

def write_solb(char *filename,
               int n_vtx,
               int n_field,
               NPY.ndarray[NPY.double_t, mode='c', ndim=1] field):
  """
  """
  PDM_write_gamma_sol(filename,
                      n_vtx,
                      n_field,
           <double *> field.data)


def write_matsym_solb(char *filename,
                      int n_vtx,
                      NPY.ndarray[NPY.double_t, mode='c', ndim=1] field):
  """
  """
  PDM_write_gamma_matsym(filename,
                         n_vtx,
              <double *> field.data)


def read_solb(char *filename,
              int n_vtx,
              int n_field,
              NPY.ndarray[NPY.double_t, mode='c', ndim=1] field):
  """
  """
  PDM_read_gamma_sol(filename,
                 n_vtx,
                 n_field,
      <double *> field.data)


def meshb_to_dmesh_nodal(char *filename, MPI.Comm comm, int fix_orientation_2d, int fix_orientation_3d):
  """
  """
  cdef PDM_dmesh_nodal_t* dmn
  cdef MPI.MPI_Comm c_comm   = comm.ob_mpi
  cdef PDM_MPI_Comm pdm_comm = PDM_MPI_mpi_2_pdm_mpi_comm(&c_comm)

  dmn = PDM_reader_gamma_dmesh_nodal(pdm_comm, filename, fix_orientation_2d, fix_orientation_3d)

  py_caps = PyCapsule_New(dmn, NULL, NULL);

  return DistributedMeshNodalCapsule(py_caps) # The free is inside the class

