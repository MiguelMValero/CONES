
cdef extern from "pdm_dmesh_nodal_to_dmesh.h":
    # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    # > Wrapping of structure
    ctypedef struct PDM_dmesh_nodal_to_dmesh_t:
      pass
    # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    ctypedef enum PDM_dmesh_nodal_to_dmesh_transform_t:
      PDM_DMESH_NODAL_TO_DMESH_TRANSFORM_TO_FACE = 0
      PDM_DMESH_NODAL_TO_DMESH_TRANSFORM_TO_EDGE = 1
    # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    ctypedef enum PDM_dmesh_nodal_to_dmesh_translate_group_t:
      PDM_DMESH_NODAL_TO_DMESH_TRANSLATE_GROUP_NONE    = 0
      PDM_DMESH_NODAL_TO_DMESH_TRANSLATE_GROUP_TO_FACE = 1
      PDM_DMESH_NODAL_TO_DMESH_TRANSLATE_GROUP_TO_EDGE = 2
      PDM_DMESH_NODAL_TO_DMESH_TRANSLATE_GROUP_TO_VTX  = 3
    # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    # > Wrapping of function
    PDM_dmesh_nodal_to_dmesh_t* PDM_dmesh_nodal_to_dmesh_create(int             n_mesh,
                                                                PDM_MPI_Comm    comm,
                                                                PDM_ownership_t owner)

    void PDM_dmesh_nodal_to_dmesh_add_dmesh_nodal(PDM_dmesh_nodal_to_dmesh_t *dmn_to_dm,
                                                  int                         i_mesh,
                                                  PDM_dmesh_nodal_t          *dmn)

    void PDM_dmesh_nodal_to_dmesh_compute(PDM_dmesh_nodal_to_dmesh_t                 *dmn_to_dm,
                                          PDM_dmesh_nodal_to_dmesh_transform_t        transform_kind,
                                          PDM_dmesh_nodal_to_dmesh_translate_group_t  transform_group_kind)

    void PDM_dmesh_nodal_to_dmesh_get_dmesh(PDM_dmesh_nodal_to_dmesh_t  *dmn_to_dm,
                                            int                          i_mesh,
                                            PDM_dmesh_t                **dm)

    void PDM_dmesh_nodal_to_dmesh_free(PDM_dmesh_nodal_to_dmesh_t* dmn_to_dm)
    # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::


_PDM_DMESH_NODAL_TO_DMESH_TRANSFORM_TO_FACE = PDM_DMESH_NODAL_TO_DMESH_TRANSFORM_TO_FACE
_PDM_DMESH_NODAL_TO_DMESH_TRANSFORM_TO_EDGE = PDM_DMESH_NODAL_TO_DMESH_TRANSFORM_TO_EDGE

_PDM_DMESH_NODAL_TO_DMESH_TRANSLATE_GROUP_NONE    = PDM_DMESH_NODAL_TO_DMESH_TRANSLATE_GROUP_NONE
_PDM_DMESH_NODAL_TO_DMESH_TRANSLATE_GROUP_TO_FACE = PDM_DMESH_NODAL_TO_DMESH_TRANSLATE_GROUP_TO_FACE
_PDM_DMESH_NODAL_TO_DMESH_TRANSLATE_GROUP_TO_EDGE = PDM_DMESH_NODAL_TO_DMESH_TRANSLATE_GROUP_TO_EDGE
_PDM_DMESH_NODAL_TO_DMESH_TRANSLATE_GROUP_TO_VTX  = PDM_DMESH_NODAL_TO_DMESH_TRANSLATE_GROUP_TO_VTX


# ------------------------------------------------------------------
cdef class DMeshNodalToDMesh:
  """
     DistributedMesh: Distributed mesh structure
  """
  # ************************************************************************
  # > Class attributes
  cdef PDM_dmesh_nodal_to_dmesh_t* dmn_to_dm
  keep_alive = list()
  # ************************************************************************
  # ------------------------------------------------------------------------
  def __cinit__(self, n_mesh,
                      MPI.Comm comm):
    """
    TODOUX
    """
    # ************************************************************************
    # > Declaration
    cdef MPI.MPI_Comm c_comm = comm.ob_mpi
    # ************************************************************************

    self.dmn_to_dm = PDM_dmesh_nodal_to_dmesh_create(n_mesh,
                                                     PDM_MPI_mpi_2_pdm_mpi_comm (<void *> &c_comm),
                                                     PDM_OWNERSHIP_USER) # Python take ownership);

  # ------------------------------------------------------------------------
  def add_dmesh_nodal(self, int i_mesh, DistributedMeshNodal dmn):
    """
    """
    self.keep_alive.append(dmn)
    PDM_dmesh_nodal_to_dmesh_add_dmesh_nodal(self.dmn_to_dm,
                                             i_mesh,
                                             dmn.dmn)

  # ------------------------------------------------------------------------
  def compute(self, PDM_dmesh_nodal_to_dmesh_transform_t       transform_kind,
                    PDM_dmesh_nodal_to_dmesh_translate_group_t transform_group_kind):
    """
    """
    PDM_dmesh_nodal_to_dmesh_compute(self.dmn_to_dm, transform_kind, transform_group_kind)

  # ------------------------------------------------------------------------
  def get_dmesh(self, int i_mesh):
    """
    """
    # ************************************************************************
    # > Declaration
    cdef PDM_dmesh_t* dm
    # ************************************************************************
    PDM_dmesh_nodal_to_dmesh_get_dmesh(self.dmn_to_dm, i_mesh, &dm)

    py_caps = PyCapsule_New(dm, NULL, NULL);

    return DistributedMeshCapsule(py_caps) # The free is inside the class

  # ------------------------------------------------------------------------
  def __dealloc__(self):
    """
       Use the free method of PDM Lib
    """
    # ************************************************************************
    # > Declaration
    # ************************************************************************
    # print("DMeshNodalToDMesh::__dealloc__")
    PDM_dmesh_nodal_to_dmesh_free(self.dmn_to_dm)
    # print("DMeshNodalToDMesh::__dealloc__ end ")
