
cdef extern from "pdm_points_merge.h":
    # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    # > Wrapping of structure
    ctypedef struct PDM_points_merge_t:
        pass
    # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    # > Wrapping of function
    PDM_points_merge_t * PDM_points_merge_create(int                   _n_point_cloud,
                                                double                tolerance,
                                                PDM_MPI_Comm          comm,
                                                const PDM_ownership_t owner)
    void PDM_points_merge_free(PDM_points_merge_t *pm)
    void PDM_points_merge_cloud_set(PDM_points_merge_t *pm,
                                    int                i_point_cloud,
                                    int                n_points,
                                    double            *coords,
                                    double            *char_length)
    void PDM_points_merge_process(PDM_points_merge_t *pm)
    void PDM_points_merge_candidates_get(PDM_points_merge_t  *pm,
                                         int                 i_point_cloud,
                                         int               **candidates_idx,
                                         int               **candidates_desc)
    void PDM_points_merge_candidates_size_get(PDM_points_merge_t *pm,
                                              int                i_point_cloud,
                                              int               *_n_point_cloud,
                                              int               *n_candidates_desc)

    void PDM_points_merge_make_interface(PDM_points_merge_t  *pm,
                                         int                 *out_n_g_interface,
                                         int                **out_interface_cloud_pair,
                                         int                **out_dn_vtx_itrf,
                                         PDM_g_num_t       ***out_itrf_gnum_cur,
                                         PDM_g_num_t       ***out_itrf_gnum_opp);
    # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

# ------------------------------------------------------------------
cdef class PointsMerge:
    """
       PointsMerge: Interface to build connection between multiple cloud in parallel
       Useful for get connection between partiton from faces coordinates
    """
    # ************************************************************************
    # > Class attributes
    cdef PDM_points_merge_t* _pm
    cdef int _size
    cdef int _rank
    cdef int n_point_cloud
    # ************************************************************************
    # ------------------------------------------------------------------------
    def __init__(self, MPI.Comm    comm,
                       int         n_point_cloud,
                       double      relative_tolerance):
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
        self._rank = comm.Get_rank()
        self._size = comm.Get_size()
        self.n_point_cloud = n_point_cloud
        # ::::::::::::::::::::::::::::::::::::::::::::::::::

        # ::::::::::::::::::::::::::::::::::::::::::::::::::
        # > Convert mpi4py -> PDM_MPI
        cdef MPI.MPI_Comm c_comm   = comm.ob_mpi
        cdef PDM_MPI_Comm pdm_comm = PDM_MPI_mpi_2_pdm_mpi_comm(<void *> &c_comm)
        # ::::::::::::::::::::::::::::::::::::::::::::::::::

        # ::::::::::::::::::::::::::::::::::::::::::::::::::
        self._pm = PDM_points_merge_create(n_point_cloud,
                                           relative_tolerance,
                                           pdm_comm,
                                           PDM_OWNERSHIP_USER) # Python take ownership
        # ::::::::::::::::::::::::::::::::::::::::::::::::::

    # ------------------------------------------------------------------------
    def cloud_set(self, int i_point_cloud,
                        int n_points,
                        NPY.ndarray[NPY.double_t, mode='c', ndim=1] coords,
                        NPY.ndarray[NPY.double_t, mode='c', ndim=1] char_length):
        """
        """
        # ************************************************************************
        # > Declaration
        # ************************************************************************

        # ::::::::::::::::::::::::::::::::::::::::::::::::::
        PDM_points_merge_cloud_set(self._pm,
                                   i_point_cloud,
                                   n_points,
                                   <double *> coords.data,
                                   <double *> char_length.data)
        # ::::::::::::::::::::::::::::::::::::::::::::::::::

    # ------------------------------------------------------------------------
    def compute(self):
        """
        """
        # ************************************************************************
        # > Declaration
        # ************************************************************************
        PDM_points_merge_process(self._pm)

    # ------------------------------------------------------------------------
    def get_merge_candidates(self, int i_point_cloud):
        """
        """
        # ************************************************************************
        # > Declaration
        cdef int          *candidates_idx
        cdef int          *candidates_desc
        cdef int           n_point_cloud
        cdef int           n_candidates_desc
        cdef NPY.npy_intp  dim
        # ************************************************************************

        PDM_points_merge_candidates_size_get(self._pm,
                                             i_point_cloud,
                                             &n_point_cloud,
                                             &n_candidates_desc)


        # > Get Size
        PDM_points_merge_candidates_get(self._pm,
                                        i_point_cloud,
                                        &candidates_idx,
                                        &candidates_desc)

        return {'candidates_idx'  : create_numpy_i(candidates_idx, n_point_cloud+1),
                'candidates_desc' : create_numpy_i(candidates_desc, 3*candidates_idx[n_point_cloud])
                }

    # ------------------------------------------------------------------------
    def make_interface(self):
        """
        """
        # ************************************************************************
        # > Declaration
        cdef int             n_g_interface
        cdef int            *interface_cloud_pair
        cdef int            *dn_vtx_itrf
        cdef PDM_g_num_t   **itrf_gnum_cur
        cdef PDM_g_num_t   **itrf_gnum_opp
        cdef NPY.npy_intp  dim
        # ************************************************************************

        PDM_points_merge_make_interface(self._pm,
                                        &n_g_interface,
                                        &interface_cloud_pair,
                                        &dn_vtx_itrf,
                                        &itrf_gnum_cur,
                                        &itrf_gnum_opp);

        np_cloud_pair  = create_numpy_i(interface_cloud_pair, 2 * n_g_interface)
        np_dn_vtx_itrf = create_numpy_i(dn_vtx_itrf         ,     n_g_interface)

        litrf_gnum_cur = list()
        litrf_gnum_opp = list()
        for i_cloud in range(n_g_interface):
          litrf_gnum_cur.append(create_numpy_g(itrf_gnum_cur[i_cloud], dn_vtx_itrf[i_cloud]))
          litrf_gnum_opp.append(create_numpy_g(itrf_gnum_opp[i_cloud], dn_vtx_itrf[i_cloud]))

        return {'np_cloud_pair'  : np_cloud_pair,
                'lgnum_cur'      : litrf_gnum_cur,
                'lgnum_opp'      : litrf_gnum_opp}


    # ------------------------------------------------------------------------
    def __dealloc__(self):
      """
         Use the free method of PDM Lib
      """
      PDM_points_merge_free(self._pm)

