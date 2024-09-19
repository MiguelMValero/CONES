
cdef extern from "pdm_part.h":

    ctypedef enum PDM_part_split_t:
        PDM_PART_SPLIT_PARMETIS = 1
        PDM_PART_SPLIT_PTSCOTCH = 2
        PDM_PART_SPLIT_HILBERT  = 3

    ctypedef struct PDM_part_t:
        pass

    # -> PPART bases functions
    # ------------------------------------------------------------------
    # MPI_Comm      comm,
    PDM_part_t * PDM_part_create(PDM_MPI_Comm          comm,
                                 PDM_part_split_t      split_method,
                                 char                  *renum_cell_method,
                                 char                  *renum_face_method,
                                 int                   nPropertyCell,
                                 int                   *renum_properties_cell,
                                 int                   nPropertyFace,
                                 int                   *renum_properties_face,
                                 int                   nPart,
                                 int                   dNCell,
                                 int                   dNFace,
                                 int                   dNVtx,
                                 int                   nFaceGroup,
                                 int                  *dCellFaceIdx,
                                 PDM_g_num_t          *dCellFace,
                                 int                  *dCellTag,
                                 int                  *dCellWeight,
                                 int                   have_dCellPart,
                                 int                  *dCellPart,
                                 PDM_g_num_t          *dFaceCell,
                                 int                  *dFaceVtxIdx,
                                 PDM_g_num_t          *dFaceVtx,
                                 int                  *dFaceTag,
                                 double               *dVtxCoord,
                                 int                  *dVtxTag,
                                 int                  *dFaceGroupIdx,
                                 PDM_g_num_t          *dFaceGroup)

    # ------------------------------------------------------------------
    void PDM_part_part_dim_get(PDM_part_t *ppart,
                               int         ipart,
                               int        *nCell,
                               int        *nFace,
                               int        *nFacePartBound,
                               int        *nVtx,
                               int        *nProc,
                               int        *nTPart,
                               int        *sCellFace,
                               int        *sFaceVtx,
                               int        *sFaceGroup,
                               int        *nFaceGroup)

    # ------------------------------------------------------------------
    void PDM_part_part_val_get(PDM_part_t    *ppart,
                               int            ipart,
                               int          **cellTag,
                               int          **cellFaceIdx,
                               int          **cellFace,
                               PDM_g_num_t **cellLNToGN,
                               int          **faceTag,
                               int          **faceCell,
                               int          **faceVtxIdx,
                               int          **faceVtx,
                               PDM_g_num_t **faceLNToGN,
                               int          **facePartBoundProcIdx,
                               int          **facePartBoundPartIdx,
                               int          **facePartBound,
                               int          **vtxTag,
                               double       **vtx,
                               PDM_g_num_t **vtxLNToGN,
                               int          **faceGroupIdx,
                               int          **faceGroup,
                               PDM_g_num_t **faceGroupLNToGN)

    # ------------------------------------------------------------------
    void PDM_part_part_color_get(PDM_part_t    *ppart,
                                 int            ipart,
                                 int          **cellColor,
                                 int          **faceColor,
                                 int          **threadColor,
                                 int          **hyperPlaneColor)

    # ------------------------------------------------------------------
    void PDM_part_free(PDM_part_t *ppart)
    void PDM_part_partial_free(PDM_part_t *ppart)

    # ------------------------------------------------------------------
    void PDM_part_time_get(PDM_part_t  *ppart,
                           double     **elapsed,
                           double     **cpu,
                           double     **cpu_user,
                           double     **cpu_sys)

    # ------------------------------------------------------------------
    void PDM_part_stat_get(PDM_part_t *ppart,
                           int        *cells_average,
                           int        *cells_median,
                           double     *cells_std_deviation,
                           int        *cells_min,
                           int        *cells_max,
                           int        *bound_part_faces_average,
                           int        *bound_part_faces_median,
                           double     *bound_part_faces_std_deviation,
                           int        *bound_part_faces_min,
                           int        *bound_part_faces_max,
                           int        *bound_part_faces_sum)

# ------------------------------------------------------------------

cdef extern from "pdm_mpi.h":
    PDM_MPI_Comm PDM_MPI_mpi_2_pdm_mpi_comm (void *mpi_comm)

# ------------------------------------------------------------------


cdef class Part:
    """
       Ppart
    """
    # > For Ppart
    cdef PDM_part_t *_ppart
    cdef int _nFaceGroup
    # ------------------------------------------------------------------
    def __cinit__(self,
                 MPI.Comm comm,
                 PDM_part_split_t split_method,
                 char                 *renum_cell_method,
                 char                 *renum_face_method,
                 int                   nPropertyCell,
                 NPY.ndarray[NPY.int32_t   , mode='c', ndim=1] renum_properties_cell,
                 int                   nPropertyFace,
                 NPY.ndarray[NPY.int32_t   , mode='c', ndim=1] renum_properties_face,
                 int nPart,
                 int dNCell,
                 int dNFace,
                 int dNVtx,
                 int have_dCellPart,
                 int nFaceGroup,
                 NPY.ndarray[NPY.int32_t   , mode='c', ndim=1] dCellFaceIdx,
                 NPY.ndarray[npy_pdm_gnum_t, mode='c', ndim=1] dCellFace,
                 NPY.ndarray[NPY.int32_t   , mode='c', ndim=1] dCellTag,
                 NPY.ndarray[NPY.int32_t   , mode='c', ndim=1] dCellWeight,
                 NPY.ndarray[NPY.int32_t   , mode='c', ndim=1] dCellPart,
                 NPY.ndarray[npy_pdm_gnum_t, mode='c', ndim=1] dFaceCell,
                 NPY.ndarray[NPY.int32_t   , mode='c', ndim=1] dFaceVtxIdx not None,
                 NPY.ndarray[npy_pdm_gnum_t, mode='c', ndim=1] dFaceVtx    not None,
                 NPY.ndarray[NPY.int32_t   , mode='c', ndim=1] dFaceTag,
                 NPY.ndarray[NPY.double_t  , mode='c', ndim=1] dVtxCoord   not None,
                 NPY.ndarray[NPY.int32_t   , mode='c', ndim=1] dVtxTag,
                 NPY.ndarray[NPY.int32_t   , mode='c', ndim=1] dFaceGroupIdx,
                 NPY.ndarray[npy_pdm_gnum_t, mode='c', ndim=1] dFaceGroupFace
                 ):

        """
        Create Ppart partition from Ppart Library ( Developed at ONERA by Eric Quemerais )

        :param comm:            MPI Communicator (Caution MPI Comm is a mpi4py object )
        :param split_method:    Partitioning method ( 1: ParMetis / 2: PT-Scotch)
        :param nPart :          Number of partition for the process
        :param dNCell:          Distribute number of Cell
        :param dNFace:          Distribute number of Face
        :param have_dCellPart:  Have a initial partitioning (0:False / 1:True)
        :param nFaceGroup:      Number of Group Face which defined boundary
        :param dCellFaceIdx:    Cell face connectivity ( can be None)


        """

        # ~> Communicator Mpi
        cdef MPI.MPI_Comm c_comm = comm.ob_mpi

        # ~> Set _nFaceGroup
        self._nFaceGroup =  nFaceGroup

        # ~> \param [out]  ppartId        ppart identifier
        # cdef int    _id

        # ~> \param [in]   dCellFaceIdx   Distributed cell face connectivity index or NULL
        cdef int * dCellFaceIdx_data
        if (dCellFaceIdx is None):
            dCellFaceIdx_data = NULL
        else:
            dCellFaceIdx_data = <int *> dCellFaceIdx.data

        # ~> \param [in]   renum_properties_cell
        cdef int * renum_properties_cell_data
        if (renum_properties_cell is None):
            renum_properties_cell_data = NULL
        else:
            renum_properties_cell_data = <int *> renum_properties_cell.data

        # ~> \param [in] renum_properties_face
        cdef int * renum_properties_face_data
        if (renum_properties_face is None):
            renum_properties_face_data = NULL
        else:
            renum_properties_face_data = <int *> renum_properties_face.data

        # ~> \param [in]   dCellFace      Distributed cell face connectivity or NULL
        cdef PDM_g_num_t * dCellFace_data = NULL
        if (dCellFace is None):
           dCellFace_data = NULL
        else:
            dCellFace_data = <PDM_g_num_t *> dCellFace.data

        # \param [in]   dCellTag       Cell tag (size : nCell) or NULL
        cdef int *dCellTag_data
        if (dCellTag is None):
            dCellTag_data = NULL
        else:
            dCellTag_data = <int *> dCellTag.data

        # \param [in]   dCellWeight    Cell weight (size : nCell) or NULL
        cdef int *dCellWeight_data
        if (dCellWeight is None):
            dCellWeight_data = NULL
        else:
            dCellWeight_data = <int *> dCellWeight.data

        # \param [in]   dCellPart      Distributed cell partitioning
        cdef int *dCellPart_data
        if (dCellPart is None):
            dCellPart_data = NULL
        else:
            dCellPart_data = <int *> dCellPart.data

        # \param [in]   dFaceCell      Distributed face cell connectivity or NULL
        cdef PDM_g_num_t * dFaceCell_data
        if (dFaceCell is None):
            dFaceCell_data = NULL
        else:
            dFaceCell_data = <PDM_g_num_t *> dFaceCell.data

        # \param [in]   dFaceTag       Distributed face tag
        cdef int *dFaceTag_data
        if (dFaceTag is None):
            dFaceTag_data = NULL
        else:
            dFaceTag_data = <int *> dFaceTag.data

        # \param [in]   dVtxCoord      Distributed vertex coordinates
        cdef int *dVtxTag_data
        if (dVtxTag is None):
            dVtxTag_data = NULL
        else:
            dVtxTag_data = <int *> dVtxTag.data

        # \param [in]   dFaceGroupIdx  Index of distributed faces list of each group
        cdef int *dFaceGroupIdx_data
        if (dFaceGroupIdx is None):
            dFaceGroupIdx_data = NULL
        else:
            dFaceGroupIdx_data = <int *> dFaceGroupIdx.data

        # \param [in]   dFaceGroup     distributed faces list of each group
        cdef PDM_g_num_t * dFaceGroupFace_data
        if (dFaceGroupFace is None):
            dFaceGroupFace_data = NULL
        else:
            dFaceGroupFace_data = <PDM_g_num_t *> dFaceGroupFace.data

        # print dFaceGroupFace.__array_interface__['data'][0]
        # LOG.info(' '*4 + " --->  nFaceGroup   : {0} ".format(nFaceGroup))
        # LOG.info(' '*4 + " --->  dVtxCoord.data      : {0} ".format(dVtxCoord.__array_interface__['data'][0]) )
        # LOG.info(' '*4 + " --->  dFaceVtxIdx.data    : {0} ".format(dFaceVtxIdx.__array_interface__['data'][0]) )
        # LOG.info(' '*4 + " --->  dFaceVtx.data       : {0} ".format(dFaceVtx.__array_interface__['data'][0]) )
        # LOG.info(' '*4 + " --->  dFaceCell.data      : {0} ".format(dFaceCell.__array_interface__['data'][0]) )
        # LOG.info(' '*4 + " --->  dFaceGroupFace.data : {0} ".format(dFaceGroupFace.__array_interface__['data'][0]) )
        # LOG.info(' '*4 + " --->  dFaceGroupIdx.data  : {0} ".format(dFaceGroupIdx.__array_interface__['data'][0]) )
        # LOG.info(' '*4 + " ---> LibPart.PDM_part_create " )
        # print 'renum_cell_method : ', renum_cell_method
        # print 'renum_face_method : ', renum_face_method
        # -> Create PPART
        self._ppart = PDM_part_create(PDM_MPI_mpi_2_pdm_mpi_comm (<void *> &c_comm),
                                      split_method,
                                      renum_cell_method,
                                      renum_face_method,
                                      nPropertyCell,
                                      renum_properties_cell_data,
                                      nPropertyFace,
                                      renum_properties_face_data,
                                      nPart,
                                      dNCell,
                                      dNFace,
                                      dNVtx,
                                      nFaceGroup,
                                      dCellFaceIdx_data,
                                      dCellFace_data,
                                      dCellTag_data,
                                      dCellWeight_data,
                                      have_dCellPart,    # -> Add this Attention !
                                      dCellPart_data,
                                      dFaceCell_data,
                                      <int *>          dFaceVtxIdx.data,
                                      <PDM_g_num_t *>  dFaceVtx.data,
                                      dFaceTag_data,
                                      <double *>       dVtxCoord.data,
                                      dVtxTag_data,
                                      dFaceGroupIdx_data,
                                      dFaceGroupFace_data)

    # ------------------------------------------------------------------
    def part_partial_free(self):
        PDM_part_partial_free(self._ppart)

    # ------------------------------------------------------------------
    def __dealloc__(self):
        PDM_part_free(self._ppart)

    # ------------------------------------------------------------------
    def part_dim_get(self, int ipart):
        """
           Get partition dimensions
        """
        # ************************************************************************
        # > Declaration
        cdef int nCell
        cdef int nProc
        cdef int nTPart
        cdef int nFace
        cdef int nFacePartBound
        cdef int nVertex
        cdef int sCellFace
        cdef int sFaceVertex
        cdef int sFaceGroup
        cdef int nFaceGroup
        # ************************************************************************

        PDM_part_part_dim_get(self._ppart,
                              ipart,
                              &nCell,
                              &nFace,
                              &nFacePartBound,
                              &nVertex,
                              &nProc,
                              &nTPart,
                              &sCellFace,
                              &sFaceVertex,
                              &sFaceGroup,
                              &nFaceGroup)

        return {'nCell'          :nCell,
                'ipart'          :ipart,
                'nFace'          :nFace,
                'nTPart'         :nTPart,
                'nProc'          :nProc,
                'nFacePartBound' :nFacePartBound,
                'nVertex'        :nVertex,
                'sCellFace'      :sCellFace,
                'sFaceVertex'    :sFaceVertex,
                'sFaceGroup'     :sFaceGroup,
                'nFaceGroup'     :nFaceGroup}

    # ------------------------------------------------------------------
    def part_val_get(self, int ipart):
        """
           Get partition dimensions
        """
        # ************************************************************************
        # > Declaration
        cdef int          *cellTag,
        cdef int          *cellFaceIdx,
        cdef int          *cellFace,
        cdef PDM_g_num_t  *cellLNToGN,
        cdef int          *faceTag,
        cdef int          *faceCell,
        cdef int          *faceVertexIdx,
        cdef int          *faceVertex,
        cdef PDM_g_num_t  *faceLNToGN,
        cdef int          *facePartBound,
        cdef int          *facePartBoundProcIdx,
        cdef int          *facePartBoundPartIdx,
        cdef int          *vertexTag,
        cdef double       *vertex,
        cdef PDM_g_num_t  *vertexLNToGN,
        cdef int          *faceGroupIdx,
        cdef int          *faceGroup,
        cdef PDM_g_num_t  *faceGroupLNToGN
        # ************************************************************************

        # dims = self.part_dim_get(self.id, ipart)
        dims = self.part_dim_get(ipart)

        # -> Call PPART to get info
        PDM_part_part_val_get(self._ppart,
                              ipart,
                              &cellTag,
                              &cellFaceIdx,
                              &cellFace,
                              &cellLNToGN,
                              &faceTag,
                              &faceCell,
                              &faceVertexIdx,
                              &faceVertex,
                              &faceLNToGN,
                              &facePartBoundProcIdx,
                              &facePartBoundPartIdx,
                              &facePartBound,
                              &vertexTag,
                              &vertex,
                              &vertexLNToGN,
                              &faceGroupIdx,
                              &faceGroup,
                              &faceGroupLNToGN)

        return {
            'npCellTag'                  : create_numpy_or_none_i(cellTag,              dims['nCell'],            flag_owndata=False),
            'npCellFaceIdx'              : create_numpy_or_none_i(cellFaceIdx,          dims['nCell']+1,          flag_owndata=False),
            'npCellFace'                 : create_numpy_or_none_i(cellFace,             dims['sCellFace'],        flag_owndata=False),
            'npCellLNToGN'               : create_numpy_or_none_g(cellLNToGN,           dims['nCell'],            flag_owndata=False),
            'npFaceTag'                  : create_numpy_or_none_i(faceTag,              dims['nFace'],            flag_owndata=False),
            'npFaceCell'                 : create_numpy_or_none_i(faceCell,             2*dims['nFace'],          flag_owndata=False),
            'npFaceVertexIdx'            : create_numpy_or_none_i(faceVertexIdx,        dims['nFace']+1,          flag_owndata=False),
            'npFaceVertex'               : create_numpy_or_none_i(faceVertex,           dims['sFaceVertex'],      flag_owndata=False),
            'npFaceLNToGN'               : create_numpy_or_none_g(faceLNToGN,           dims['nFace'],            flag_owndata=False),
            'npfacePartBoundProcIdx'     : create_numpy_or_none_i(facePartBoundProcIdx, dims['nProc']+1,          flag_owndata=False),
            'npfacePartBoundPartIdx'     : create_numpy_or_none_i(facePartBoundPartIdx, dims['nTPart']+1,         flag_owndata=False),
            'npFacePartBound'            : create_numpy_or_none_i(facePartBound,        4*dims['nFacePartBound'], flag_owndata=False),
            'npVertexTag'                : create_numpy_or_none_i(vertexTag,            dims['nVertex'],          flag_owndata=False),
            'npVertex'                   : create_numpy_or_none_d(vertex,               3*dims['nVertex'],        flag_owndata=False),
            'npVertexLNToGN'             : create_numpy_or_none_g(vertexLNToGN,         dims['nVertex'],          flag_owndata=False),
            'npFaceGroupIdx'             : create_numpy_or_none_i(faceGroupIdx,         self._nFaceGroup+1,       flag_owndata=False),
            'npFaceGroup'                : create_numpy_or_none_i(faceGroup,            dims['sFaceGroup'],       flag_owndata=False),
            'npFaceGroupLNToGN'          : create_numpy_or_none_g(faceGroupLNToGN,      dims['sFaceGroup'],       flag_owndata=False),
        }

    # ------------------------------------------------------------------
    def part_color_get(self, int ipart):
        """
           Get partition dimensions
        """
        # ************************************************************************
        # > Declaration
        cdef int          *cellColor,
        cdef int          *faceColor
        cdef int          *threadColor
        cdef int          *hyperPlaneColor
        # ************************************************************************

        # dims = self.part_dim_get(self.id, ipart)
        dims = self.part_dim_get(ipart)

        # -> Call PPART to get info
        PDM_part_part_color_get(self._ppart,
                                ipart,
                                &cellColor,
                                &faceColor,
                                &threadColor,
                                &hyperPlaneColor)
        return {
            'npCellColor'       : create_numpy_or_none_i(cellColor,       dims['nCell'], flag_owndata=False),
            'npFaceColor'       : create_numpy_or_none_i(faceColor,       dims['nFace'], flag_owndata=False),
            'npThreadColor'     : create_numpy_or_none_i(threadColor,     dims['nCell'], flag_owndata=False),
            'npHyperPlaneColor' : create_numpy_or_none_i(hyperPlaneColor, dims['nCell'], flag_owndata=False),
        }

    # ------------------------------------------------------------------
    def part_time_get(self):
        """
        Get times
        """
        # ************************************************************************
        # > Declaration
        cdef double *elapsed
        cdef double *cpu
        cdef double *cpu_user
        cdef double *cpu_sys
        # ************************************************************************

        PDM_part_time_get(self._ppart, &elapsed, &cpu, &cpu_user, &cpu_sys)

        d_elapsed = {'total'              : elapsed[0],
                     'building graph'     : elapsed[1],
                     'splitting graph'    : elapsed[2],
                     'building partitions': elapsed[3]}

        d_cpu     = {'total'              : cpu[0],
                     'building graph'     : cpu[1],
                     'splitting graph'    : cpu[2],
                     'building partitions': cpu[3]}

        d_cpu_user = {'total'              : cpu_user[0],
                      'building graph'     : cpu_user[1],
                      'splitting graph'    : cpu_user[2],
                      'building partitions': cpu_user[3]}

        d_cpu_sys = {'total'              : cpu_sys[0],
                     'building graph'     : cpu_sys[1],
                     'splitting graph'    : cpu_sys[2],
                     'building partitions': cpu_sys[3]}

        return {'elapsed' : d_elapsed, 'cpu' : d_cpu, 'cpu_user' : d_cpu_user,  'cpu_sys' : d_cpu_sys}


    # ------------------------------------------------------------------
    def part_stat_get(self):
        """
        Get statistics
        """
        # ************************************************************************
        # > Declaration
        cdef int      cells_average,
        cdef int      cells_median,
        cdef double   cells_std_deviation,
        cdef int      cells_min,
        cdef int      cells_max,
        cdef int      bound_part_faces_average,
        cdef int      bound_part_faces_median,
        cdef double   bound_part_faces_std_deviation,
        cdef int      bound_part_faces_min,
        cdef int      bound_part_faces_max,
        cdef int      bound_part_faces_sum
        # ************************************************************************

        PDM_part_stat_get(self._ppart,
                          &cells_average,
                          &cells_median,
                          &cells_std_deviation,
                          &cells_min,
                          &cells_max,
                          &bound_part_faces_average,
                          &bound_part_faces_median,
                          &bound_part_faces_std_deviation,
                          &bound_part_faces_min,
                          &bound_part_faces_max,
                          &bound_part_faces_sum)

        return {'cells_average'                  : cells_average,
                'cells_median'                   : cells_median,
                'cells_std_deviation'            : cells_std_deviation,
                'cells_min'                      : cells_min,
                'cells_max'                      : cells_max,
                'bound_part_faces_average'       : bound_part_faces_average,
                'bound_part_faces_median'        : bound_part_faces_median,
                'bound_part_faces_std_deviation' : bound_part_faces_std_deviation,
                'bound_part_faces_min'           : bound_part_faces_min,
                'bound_part_faces_max'           : bound_part_faces_max,
                'bound_part_faces_sum'           : bound_part_faces_sum}
