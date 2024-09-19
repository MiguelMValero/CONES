cdef extern from "pdm_dcube_gen.h":
    # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    # > Wrapping of Ppart Structure
    ctypedef struct PDM_dcube_t:
      pass
    # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    # ------------------------------------------------------------------
    PDM_dcube_t* PDM_dcube_gen_init(PDM_MPI_Comm        comm,
                                    const PDM_g_num_t   n_vtx_seg,
                                    const double        length,
                                    const double        zero_x,
                                    const double        zero_y,
                                    const double        zero_z,
                                    PDM_ownership_t     owner)

    # ------------------------------------------------------------------
    void PDM_dcube_gen_dim_get(PDM_dcube_t        *pdm_dcube,
                               int                *n_face_group,
                               int                *dn_cell,
                               int                *dn_face,
                               int                *dn_vtx,
                               int                *sface_vtx,
                               int                *sface_group)

    # ------------------------------------------------------------------
    void PDM_dcube_gen_data_get(PDM_dcube_t        *pdm_dcube,
                                PDM_g_num_t       **dface_cell,
                                int               **dface_vtx_idx,
                                PDM_g_num_t       **dface_vtx,
                                double            **dvtx_coord,
                                int               **dface_group_idx,
                                PDM_g_num_t       **dface_group)

    # ------------------------------------------------------------------
    void PDM_dcube_gen_free(PDM_dcube_t        *pdm_dcube)
    # ------------------------------------------------------------------

# ------------------------------------------------------------------
cdef class DCubeGenerator:
    """
       DCubeGenerator
    """
    # > For Ppart
    cdef PDM_dcube_t* _dcube
    # ------------------------------------------------------------------
    def __cinit__(self,
                  npy_pdm_gnum_t                         n_vtx_seg,
                  NPY.double_t                           length,
                  NPY.double_t                           zero_x,
                  NPY.double_t                           zero_y,
                  NPY.double_t                           zero_z,
                  MPI.Comm                               comm):

        """
        """
        # ~> Communicator Mpi
        cdef MPI.MPI_Comm c_comm = comm.ob_mpi

        # -> Create dcube
        self._dcube = PDM_dcube_gen_init(PDM_MPI_mpi_2_pdm_mpi_comm (<void *> &c_comm),
                                         n_vtx_seg,
                                         length,
                                         zero_x,
                                         zero_y,
                                         zero_z,
                                         PDM_OWNERSHIP_USER) # Python take owership

    # ------------------------------------------------------------------
    def __dealloc__(self):
        PDM_dcube_gen_free(self._dcube)

    # ------------------------------------------------------------------
    def dcube_dim_get(self):
        """
           Get dcube dimensions
        """
        # ************************************************************************
        # > Declaration
        cdef int n_face_group
        cdef int dn_cell
        cdef int dn_face
        cdef int dn_vtx
        cdef int sface_vtx
        cdef int sface_group
        # ************************************************************************

        PDM_dcube_gen_dim_get(self._dcube,
                              &n_face_group,
                              &dn_cell,
                              &dn_face,
                              &dn_vtx,
                              &sface_vtx,
                              &sface_group)

        return {'n_face_group' : n_face_group,
                'dn_cell'      : dn_cell,
                'dn_face'      : dn_face,
                'dn_vtx'       : dn_vtx,
                'sface_vtx'    : sface_vtx,
                'sface_group'  : sface_group}

    # ------------------------------------------------------------------
    def dcube_val_get(self):
        """
           Get dcube data
        """
        # ************************************************************************
        # > Declaration
        cdef NPY.npy_intp dim
        cdef PDM_g_num_t  *dface_cell
        cdef int          *dface_vtx_idx
        cdef PDM_g_num_t  *dface_vtx
        cdef double       *dvtx_coord
        cdef int          *dface_group_idx
        cdef PDM_g_num_t  *dface_group
        # ************************************************************************

        dims = self.dcube_dim_get()

        PDM_dcube_gen_data_get(self._dcube,
                               &dface_cell,
                               &dface_vtx_idx,
                               &dface_vtx,
                               &dvtx_coord,
                               &dface_group_idx,
                               &dface_group)

        return {'dface_cell'      : create_numpy_g(dface_cell,      2*dims['dn_face']),
                'dface_vtx_idx'   : create_numpy_i(dface_vtx_idx,   dims['dn_face']+1),
                'dface_vtx'       : create_numpy_g(dface_vtx,       dims['sface_vtx']),
                'dvtx_coord'      : create_numpy_d(dvtx_coord,      3*dims['dn_vtx']),
                'dface_group_idx' : create_numpy_i(dface_group_idx, dims['n_face_group']+1),
                'dface_group'     : create_numpy_g(dface_group,     dims['sface_group'])}
