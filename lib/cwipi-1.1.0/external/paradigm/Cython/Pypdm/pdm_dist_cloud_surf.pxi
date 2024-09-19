
cdef extern from "pdm_dist_cloud_surf.h":
    # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    # > Wrapping of Ppart Structure
    ctypedef struct PDM_dist_cloud_surf_t:
        pass
    # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    # > Wrapping of function
    PDM_dist_cloud_surf_t* PDM_dist_cloud_surf_create(PDM_mesh_nature_t mesh_nature,
                                                      int               n_point_cloud,
                                                      PDM_MPI_Comm      comm,
                                                      PDM_ownership_t   owner)

    void PDM_dist_cloud_surf_n_part_cloud_set(PDM_dist_cloud_surf_t *dist,
                                              int                    i_point_cloud,
                                              int                    n_part)

    void PDM_dist_cloud_surf_cloud_set(PDM_dist_cloud_surf_t *dist,
                                       int                    i_point_cloud,
                                       int                    i_part,
                                       int                    n_points,
                                       double                *coords,
                                       PDM_g_num_t           *gnum)

#    void PDM_dist_cloud_surf_nodal_mesh_set(PDM_dist_cloud_surf_t *dist,
#                                            int                    mesh_nodal_id)

#    void PDM_dist_cloud_surf_surf_mesh_map(PDM_dist_cloud_surf_t *dist,
#                                           PDM_surf_mesh_t       *surf_mesh)

    void PDM_dist_cloud_surf_surf_mesh_global_data_set(PDM_dist_cloud_surf_t *dist,
                                                       int                    n_part)

    void PDM_dist_cloud_surf_surf_mesh_part_set(PDM_dist_cloud_surf_t *dist,
                                                int                    i_part,
                                                int                     n_face,
                                                int                   *face_vtx_idx,
                                                int                   *face_vtx,
                                                PDM_g_num_t           *face_ln_to_gn,
                                                int                    n_vtx,
                                                double                *coords,
                                                PDM_g_num_t           *vtx_ln_to_gn)

    void PDM_dist_cloud_surf_compute(PDM_dist_cloud_surf_t *dist)

    void PDM_dist_cloud_surf_get(PDM_dist_cloud_surf_t   *dist,
                                 int                      i_point_cloud,
                                 int                      i_part,
                                 double                **closest_elt_distance,
                                 double                **closest_elt_projected,
                                 PDM_g_num_t           **closest_elt_gnum)

    void PDM_dist_cloud_surf_free(PDM_dist_cloud_surf_t   *dist)

    void PDM_dist_cloud_surf_dump_times(PDM_dist_cloud_surf_t   *dist)


    # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

# :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

# ------------------------------------------------------------------
cdef class DistCloudSurf:
    """
    Define a method to compute the distance from a cloud to a surface
    """
    # ************************************************************************
    # > Class attributes
    cdef PDM_dist_cloud_surf_t* _dist
    cdef MPI.Comm               _comm,
    cdef PDM_mesh_nature_t      _mesh_nature
    cdef int                    _n_part_surf
    cdef list                   _nb_pts
    cdef dict                   _results
    # ************************************************************************
    # ------------------------------------------------------------------
    def __cinit__(self,
                  MPI.Comm          comm,
                  PDM_mesh_nature_t mesh_nature,
                  int               n_part_surf=1,
                  list              point_clouds=[]):
        """
        Compute the distance from point clouds to a surface
        """
        #print(f"[{comm.rank}] cinit DistCloudSurf object, n_part_surf = {n_part_surf}, point_clouds = {point_clouds} (len = {len(point_clouds)})")
        self._nb_pts      = []
        self._comm        = comm
        self._mesh_nature = mesh_nature
        self.n_part_surf  = n_part_surf
        #print(f"[{comm.rank}] OK")
        if len(point_clouds) > 0:
            # NB: Create paradigm structure
            self.n_point_cloud = len(point_clouds)
            # Set number of partition per point cloud
            for i_point_cloud, n_part_cloud in enumerate(point_clouds):
                self.set_n_part_cloud(i_point_cloud, n_part_cloud)
            # print("self.n_part_surf = ",self.n_part_surf,", self.n_point_cloud = ",self.n_point_cloud)
            # for i_point_cloud in range(self.n_point_cloud):
                # print("self.get_n_part_cloud(",i_point_cloud,") = ",self.get_n_part_cloud(i_point_cloud))

    def _create_pdm_structure(self):
        # Free memory if it has been previously allocated
        if self._dist:
            PDM_dist_cloud_surf_free(self._dist)

        # Convert mpi4py -> PDM_MPI
        cdef MPI.MPI_Comm c_comm   = self._comm.ob_mpi
        cdef PDM_MPI_Comm pdm_comm = PDM_MPI_mpi_2_pdm_mpi_comm(&c_comm)

        # NB: Python take ownership
        self._dist = PDM_dist_cloud_surf_create(self._mesh_nature,
                                                self.n_point_cloud,
                                                pdm_comm,
                                                PDM_OWNERSHIP_USER)

    # ------------------------------------------------------------------
    def get_n_part_surf(self):
        return self._n_part_surf

    def set_n_part_surf(self, value):
        if value >= 0:
            self._n_part_surf = <int> value
        else:
            raise ValueError("n_part_surf must be >= 1, but {} given here.".format(value))
    n_part_surf = property(get_n_part_surf, set_n_part_surf)

    def get_n_point_cloud(self):
        return len(self._nb_pts)

    def set_n_point_cloud(self, value):
        if value != self.get_n_point_cloud():
            if value >= 1:
                self._nb_pts = [None] * value
                self._create_pdm_structure()
            else:
                raise ValueError("n_point_cloud must be >= 1, but {} given here.".format(value))
    n_point_cloud = property(get_n_point_cloud, set_n_point_cloud)

    def get_n_part_cloud(self,
                         int i_point_cloud):
        if self._nb_pts[i_point_cloud] is not None:
            return self._nb_pts[i_point_cloud].size
        else:
            return None

    def set_n_part_cloud(self,
                         int i_point_cloud,
                         int n_part_cloud):
        """
        Give the number of partitions of a point cloud
        """
        if n_part_cloud >= 0:
            if i_point_cloud >= 0 and i_point_cloud < self.n_point_cloud:
                PDM_dist_cloud_surf_n_part_cloud_set(self._dist, i_point_cloud, n_part_cloud)
                self._nb_pts[i_point_cloud] = NPY.zeros(n_part_cloud, dtype='int32', order='C')
            else:
                raise ValueError("i_point_cloud must be >= 0 and < {}, but {} given here.".format(self.n_point_cloud, i_point_cloud))
        else:
            raise ValueError("n_part_cloud must be >= 1, but {} given here.".format(n_part_cloud))

    # ------------------------------------------------------------------
    def cloud_set(self,
                  int i_point_cloud,
                  int i_part,
                  int n_points,
                  NPY.ndarray[NPY.double_t  , mode='c', ndim=1] coords not None,
                  NPY.ndarray[npy_pdm_gnum_t, mode='c', ndim=1] gnum not None):
        """
        Give the properties of a partition of a point cloud
        """
        PDM_dist_cloud_surf_cloud_set(self._dist,
                                      i_point_cloud,
                                      i_part,
                                      n_points,
                                      <double *> coords.data,
                                      <PDM_g_num_t *> gnum.data)
        self._nb_pts[i_point_cloud][i_part] = n_points

    # ------------------------------------------------------------------
    def surf_mesh_global_data_set(self):
        """
        Give the global properties of the surface mesh
        """
        PDM_dist_cloud_surf_surf_mesh_global_data_set(self._dist,
                                                      self.n_part_surf)


    def surf_mesh_part_set(self,
                           i_part,
                           n_face,
                           NPY.ndarray[NPY.int32_t   , mode='c', ndim=1] face_vtx_idx not None,
                           NPY.ndarray[NPY.int32_t   , mode='c', ndim=1] face_vtx not None,
                           NPY.ndarray[npy_pdm_gnum_t, mode='c', ndim=1] face_ln_to_gn not None,
                           n_vtx,
                           NPY.ndarray[NPY.double_t  , mode='c', ndim=1] coords not None,
                           NPY.ndarray[npy_pdm_gnum_t, mode='c', ndim=1] vtx_ln_to_gn not None):
        """
        Give the properties of a partition of the surface mesh
        """
        PDM_dist_cloud_surf_surf_mesh_part_set (self._dist,
                                                i_part,
                                                n_face,
                                                <int *> face_vtx_idx.data,
                                                <int *> face_vtx.data,
                                                <PDM_g_num_t *> face_ln_to_gn.data,
                                                n_vtx,
                                                <double *> coords.data,
                                                <PDM_g_num_t *> vtx_ln_to_gn.data)

    # ------------------------------------------------------------------
    def compute(self):
        """
        Compute distance
        """
        cdef double *closest_elt_distance
        cdef double *closest_elt_projected
        cdef PDM_g_num_t *closest_elt_gnum

        PDM_dist_cloud_surf_compute(self._dist)

        # Store results in numpy directory
        self._results = {} # [None]*self.n_point_cloud
        for i_point_cloud in range(self.n_point_cloud):
            n_part_cloud = self.get_n_part_cloud(i_point_cloud)
            self._results[i_point_cloud] = {} # [None]*n_part_cloud
            for i_part_cloud in range(n_part_cloud):
                PDM_dist_cloud_surf_get(self._dist,
                                        i_point_cloud,
                                        i_part_cloud,
                                        &closest_elt_distance,
                                        &closest_elt_projected,
                                        &closest_elt_gnum)

                # Encapsulate C array into a numpy array
                n_pts = self._nb_pts[i_point_cloud][i_part_cloud]
                np_closest_elt_distance  = create_numpy_d(closest_elt_distance, n_pts)
                np_closest_elt_projected = create_numpy_d(closest_elt_projected, 3 * n_pts)
                np_closest_elt_gnum      = create_numpy_g(closest_elt_gnum, n_pts)

                dresults = {'ClosestEltDistance'    : np_closest_elt_distance,
                            'ClosestEltProjected'   : np_closest_elt_projected,
                            'ClosestEltGnum'        : np_closest_elt_gnum}
                self._results[i_point_cloud][i_part_cloud] = dresults

    # ------------------------------------------------------------------
    def get(self, int i_point_cloud, int i_part_cloud):
        """
        Return the properties of the closest surface element
        (distance, projected point coordinates and global numbering)
        """
        if self._results is None:
            self.compute()
        if i_point_cloud >= self.n_point_cloud:
            raise KeyError (f"Only {self.n_point_cloud} results are available, query item {i_point_cloud} here.")
        if i_part_cloud >= self.get_n_part_cloud(i_point_cloud):
            raise KeyError (f"Only {self.get_n_part_cloud(i_point_cloud)} results are available, query item {i_part_cloud} here.")
        return self._results[i_point_cloud][i_part_cloud]

    # ------------------------------------------------------------------
    def dump_times(self):
        """
        Print elapsed time
        """
        PDM_dist_cloud_surf_dump_times(self._dist)

    # ------------------------------------------------------------------
    def __dealloc__(self):
        """
        Free a distance mesh structure
        """
        PDM_dist_cloud_surf_free(self._dist)
