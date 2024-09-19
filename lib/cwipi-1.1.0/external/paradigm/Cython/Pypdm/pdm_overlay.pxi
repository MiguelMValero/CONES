
cdef extern from "pdm_overlay.h":
  ctypedef struct PDM_ol_t:
      pass

  # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  # > Wrapping of enum
  ctypedef enum PDM_ol_mesh_t:
    PDM_OL_MESH_A = 0
    PDM_OL_MESH_B = 1

  ctypedef enum PDM_ol_parameter_t:
    PDM_OL_CAR_LENGTH_TOL = 0
    PDM_OL_EXTENTS_TOL    = 1
    PDM_OL_SAME_PLANE_TOL = 2

  ctypedef enum PDM_ol_mv_t:
    PDM_OL_MV_TRANSFORMATION = 0
    PDM_OL_MV_UNKNOWN        = 1
  # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

  # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  # > Wrapping of function
  PDM_ol_t* PDM_ol_create(int          n_partMeshA,
                    int          n_partMeshB,
                    double       projectCoeff,
                    PDM_MPI_Comm comm)

  void PDM_ol_parameter_set(PDM_ol_t*          ol,
                            PDM_ol_parameter_t parameter,
                            double             value)

  void PDM_ol_input_mesh_set(PDM_ol_t*          ol,
                             PDM_ol_mesh_t      mesh,
                             int                i_part,
                             int                n_face,
                             int               *face_vtx_idx,
                             int               *face_vtx,
                             PDM_g_num_t       *face_ln_to_gn,
                             int                n_vtx,
                             double            *coords,
                             PDM_g_num_t       *vtx_ln_to_gn)

  void PDM_ol_moving_type_set(PDM_ol_t*     ol,
                              PDM_ol_mesh_t mesh,
                              PDM_ol_mv_t   mv)

  void PDM_ol_translation_set(PDM_ol_t*     ol,
                              double       *vect,
                              double       *center)

  void PDM_ol_rotation_set(PDM_ol_t*     ol,
                           double  *direction,
                           double  *center,
                           double   angle)

  void PDM_ol_compute(PDM_ol_t*     ol)

  void PDM_ol_mesh_dim_get(PDM_ol_t*     ol,
                           PDM_ol_mesh_t   mesh,
                           PDM_g_num_t    *nGOlFace,
                           PDM_g_num_t    *nGOlVtx)

  void PDM_ol_part_mesh_dim_get(PDM_ol_t*     ol,
                                PDM_ol_mesh_t  mesh,
                                int            i_part,
                                int           *nOlFace,
                                int           *nOlLinkedFace,
                                int           *nOlVtx,
                                int           *sOlFaceIniVtx,
                                int           *sOlface_vtx,
                                int           *sInit_to_ol_face)

  void PDM_ol_mesh_entities_get(PDM_ol_t*        ol,
                                PDM_ol_mesh_t    mesh,
                                int              i_part,
                                int            **olFaceIniVtxIdx,
                                int            **olFaceIniVtx,
                                int            **olface_vtx_idx,
                                int            **olface_vtx,
                                int            **olLinkedface_procIdx,
                                int            **olLinkedFace,
                                PDM_g_num_t    **olface_ln_to_gn,
                                double         **olCoords,
                                PDM_g_num_t    **olvtx_ln_to_gn,
                                int            **init_to_ol_face_idx,
                                int            **init_to_ol_face)

  void PDM_ol_del(PDM_ol_t*     ol)

  void PDM_ol_dump_times(PDM_ol_t*     ol)
  # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

# ------------------------------------------------------------------
cdef class Overlay:
  """
  """
  # ************************************************************************
  # > Class attributes
  cdef PDM_ol_t *_ol
  cdef int       _size
  cdef int       _rank
  cdef int       _init_face[2] # Store the init_face for mesh_a / mesh_b
  # ************************************************************************

  # ------------------------------------------------------------------------
  def __init__(self, int          n_partMeshA,
                     int          n_partMeshB,
                     double       projectCoeff,
                     MPI.Comm     comm):
    """
    """

    # ::::::::::::::::::::::::::::::::::::::::::::::::::
    self._rank = comm.Get_rank()
    self._size = comm.Get_size()
    # ::::::::::::::::::::::::::::::::::::::::::::::::::

    # ::::::::::::::::::::::::::::::::::::::::::::::::::
    # > Convert mpi4py -> PDM_MPI
    cdef MPI.MPI_Comm c_comm = comm.ob_mpi
    cdef PDM_MPI_Comm PDMC   = PDM_MPI_mpi_2_pdm_mpi_comm(<void *> &c_comm)
    # ::::::::::::::::::::::::::::::::::::::::::::::::::

    self._ol = PDM_ol_create(n_partMeshA,
                             n_partMeshB,
                             projectCoeff,
                             PDMC)

  # ------------------------------------------------------------------------
  def parameter_set(self, PDM_ol_parameter_t parameter,
                          double             value):
    """
    """
    PDM_ol_parameter_set(self._ol, parameter, value)

  # ------------------------------------------------------------------------
  def input_mesh_set(self, PDM_ol_mesh_t                                 mesh,
                           int                                           i_part,
                           int                                           n_face,
                           NPY.ndarray[NPY.int32_t, mode='c', ndim=1]    face_vtx_idx,
                           NPY.ndarray[NPY.int32_t, mode='c', ndim=1]    face_vtx,
                           NPY.ndarray[npy_pdm_gnum_t, mode='c', ndim=1] face_ln_to_gn,
                           int                                           n_vtx,
                           NPY.ndarray[NPY.double_t  , mode='c', ndim=1] coords,
                           NPY.ndarray[npy_pdm_gnum_t, mode='c', ndim=1] vtx_ln_to_gn):
    """
    """
    self._init_face[<int>mesh] = n_face
    PDM_ol_input_mesh_set(self._ol, mesh,
                          i_part,
                          n_face,
           <int*>         face_vtx_idx.data,
           <int*>         face_vtx.data,
           <PDM_g_num_t*> face_ln_to_gn.data,
                          n_vtx,
           <double*>      coords.data,
           <PDM_g_num_t*> vtx_ln_to_gn.data)

  # ------------------------------------------------------------------------
  def moving_type_set(self, PDM_ol_mesh_t mesh,
                            PDM_ol_mv_t   mv):
    """
    """
    PDM_ol_moving_type_set(self._ol, mesh, mv)

  # ------------------------------------------------------------------------
  def translation_set(self, NPY.ndarray[NPY.double_t  , mode='c', ndim=1] vect,
                            NPY.ndarray[NPY.double_t  , mode='c', ndim=1] center):
    """
    """
    PDM_ol_translation_set(self._ol,
                 <double*> vect.data,
                 <double*> center.data)

  # ------------------------------------------------------------------------
  def rotation_set(self, NPY.ndarray[NPY.double_t  , mode='c', ndim=1] direction,
                         NPY.ndarray[NPY.double_t  , mode='c', ndim=1] center,
                         double angle):
    """
    """
    PDM_ol_rotation_set(self._ol,
              <double*> direction.data,
              <double*> center.data,
                        angle)

  # ------------------------------------------------------------------------
  def compute(self):
    """
    """
    PDM_ol_compute(self._ol)

  # ------------------------------------------------------------------------
  def mesh_dim_get(self, PDM_ol_mesh_t mesh):
    """
    """
    # ************************************************************************
    # > Declaration
    cdef PDM_g_num_t nGOlFace
    cdef PDM_g_num_t nGOlVtx
    # ************************************************************************

    PDM_ol_mesh_dim_get(self._ol, mesh, &nGOlFace, &nGOlVtx)

    return {'nGOlFace' : nGOlFace,
            'nGOlVtx'  : nGOlVtx
            }

  # ------------------------------------------------------------------------
  def part_mesh_dim_get(self, PDM_ol_mesh_t mesh,
                              int           i_part):
    """
    """
    # ************************************************************************
    # > Declaration
    cdef int nOlFace
    cdef int nOlLinkedFace
    cdef int nOlVtx
    cdef int sOlFaceIniVtx
    cdef int sOlface_vtx
    cdef int sInit_to_ol_face
    # ************************************************************************

    PDM_ol_part_mesh_dim_get(self._ol, mesh, i_part,
                             &nOlFace,
                             &nOlLinkedFace,
                             &nOlVtx,
                             &sOlFaceIniVtx,
                             &sOlface_vtx,
                             &sInit_to_ol_face)

    return {'nOlFace'          : nOlFace,
            'nOlLinkedFace'    : nOlLinkedFace,
            'nOlVtx'           : nOlVtx,
            'sOlFaceIniVtx'    : sOlFaceIniVtx,
            'sOlface_vtx'      : sOlface_vtx,
            'sInit_to_ol_face' : sInit_to_ol_face
            }

  # ------------------------------------------------------------------------
  def part_mesh_entities_get(self, PDM_ol_mesh_t mesh,
                                   int           i_part):
    """
    """
    # ************************************************************************
    # > Declaration
    cdef int          nOlFace
    cdef int          nOlLinkedFace
    cdef int          nOlVtx
    cdef int          sOlFaceIniVtx
    cdef int          sOlface_vtx
    cdef int          sInit_to_ol_face
    cdef int         *olFaceIniVtxIdx
    cdef int         *olFaceIniVtx
    cdef int         *olface_vtx_idx
    cdef int         *olface_vtx
    cdef int         *olLinkedface_procIdx
    cdef int         *olLinkedFace
    cdef PDM_g_num_t *olface_ln_to_gn,
    cdef double      *olCoords
    cdef PDM_g_num_t *olvtx_ln_to_gn,
    cdef int         *init_to_ol_face_idx
    cdef int         *init_to_ol_face
    # ************************************************************************

    # > Get dims first
    PDM_ol_part_mesh_dim_get(self._ol, mesh, i_part,
                             &nOlFace,
                             &nOlLinkedFace,
                             &nOlVtx,
                             &sOlFaceIniVtx,
                             &sOlface_vtx,
                             &sInit_to_ol_face)

    PDM_ol_mesh_entities_get(self._ol, mesh, i_part,
                             &olFaceIniVtxIdx,
                             &olFaceIniVtx,
                             &olface_vtx_idx,
                             &olface_vtx,
                             &olLinkedface_procIdx,
                             &olLinkedFace,
                             &olface_ln_to_gn,
                             &olCoords,
                             &olvtx_ln_to_gn,
                             &init_to_ol_face_idx,
                             &init_to_ol_face)

    return {
        "olFaceIniVtxIdx"      : create_numpy_i (olFaceIniVtxIdx,      sInit_to_ol_face+1,           flag_owndata=False),
        "olFaceIniVtx"         : create_numpy_i (olFaceIniVtx,         sInit_to_ol_face,             flag_owndata=False),
        "olface_vtx_idx"       : create_numpy_i (olface_vtx_idx,       nOlFace+1,                    flag_owndata=False),
        "olface_vtx"           : create_numpy_i (olface_vtx,           olface_vtx_idx[nOlFace],      flag_owndata=False),
        "olLinkedface_procIdx" : create_numpy_i (olLinkedface_procIdx, self._size+1,                 flag_owndata=False),
        "olLinkedFace"         : create_numpy_i (olLinkedFace,         4*nOlLinkedFace,              flag_owndata=False),
        "olface_ln_to_gn,"     : create_numpy_g (olface_ln_to_gn,      nOlFace,                      flag_owndata=False),
        "olCoords"             : create_numpy_d (olCoords,             3*nOlVtx,                     flag_owndata=False),
        "olvtx_ln_to_gn,"      : create_numpy_g (olvtx_ln_to_gn,       nOlVtx,                       flag_owndata=False),
        "init_to_ol_face_idx"  : create_numpy_i (init_to_ol_face_idx,  self._init_face[<int>mesh]+1, flag_owndata=False),
        "init_to_ol_face"      : create_numpy_i (init_to_ol_face,      sInit_to_ol_face,             flag_owndata=False),
    }

  # ------------------------------------------------------------------------
  def __dealloc__(self):
    """
       Use the free method of PDM Lib
    """
    # ************************************************************************
    # > Declaration
    # ************************************************************************
    PDM_ol_del(self._ol)
