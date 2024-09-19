
cdef extern from "pdm_dconnectivity_transform.h":
    # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    void PDM_dconnectivity_to_extract_dconnectivity_bis(PDM_MPI_Comm    comm,
                                                    int             n_selected_entity1,
                                                    PDM_g_num_t    *select_entity1,
                                                    PDM_g_num_t    *entity1_distribution,
                                                    int            *dentity1_entity2_idx,
                                                    PDM_g_num_t    *dentity1_entity2,
                                                    PDM_g_num_t   **extract_entity1_distribution,
                                                    PDM_g_num_t   **extract_entity2_distribution,
                                                    int           **dextract_entity1_entity2_idx,
                                                    PDM_g_num_t   **dextract_entity1_entity2,
                                                    PDM_g_num_t   **dparent_entity1_g_num,
                                                    PDM_g_num_t   **dparent_entity2_g_num,
                                                    PDM_g_num_t   **entity1_old_to_new)

    void PDM_dconnectivity_dface_vtx_from_face_and_edge(const PDM_MPI_Comm    comm,
                                                        PDM_g_num_t    *distrib_face,
                                                        PDM_g_num_t    *distrib_edge,
                                                        int            *dface_edge_idx,
                                                        PDM_g_num_t    *dface_edge,
                                                        PDM_g_num_t    *dedge_vtx,
                                                        PDM_g_num_t   **dface_vtx)
    # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

cdef extern from "pdm_distrib.h":
    # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    PDM_g_num_t* PDM_compute_entity_distribution(PDM_MPI_Comm    comm, int dn_entity)

cdef extern from "pdm_dgeom_elem.h":
    # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    void PDM_compute_dface_normal(int          *dface_vtx_idx,
                                  PDM_g_num_t  *dface_vtx,
                                  int           dn_face,
                                  PDM_g_num_t  *dvtx_distrib,
                                  double       *dvtx_coord,
                                  double       *dface_normal,
                                  PDM_MPI_Comm  comm)
    void PDM_compute_vtx_characteristic_length(PDM_MPI_Comm    comm,
                                               int             dn_face,
                                               int             dn_edge,
                                               int             dn_vtx,
                                               int            *dface_vtx_idx,
                                               PDM_g_num_t    *dface_vtx,
                                               PDM_g_num_t    *dedge_vtx,
                                               double         *dvtx_coord,
                                               double        **dchar_length_out);
    # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::


# ===================================================================================
def compute_dface_normal(MPI.Comm                                      comm,
                         NPY.ndarray[NPY.int32_t   , mode='c', ndim=1] dface_vtx_idx,
                         NPY.ndarray[npy_pdm_gnum_t, mode='c', ndim=1] dface_vtx,
                         int                                           dn_face,
                         NPY.ndarray[npy_pdm_gnum_t, mode='c', ndim=1] dvtx_distrib,
                         NPY.ndarray[NPY.double_t  , mode='c', ndim=1] dvtx_coord   not None,
                         NPY.ndarray[NPY.double_t  , mode='c', ndim=1] dface_normal not None):
    # ::::::::::::::::::::::::::::::::::::::::::::::::::
    # > Convert mpi4py -> PDM_MPI
    cdef MPI.MPI_Comm c_comm = comm.ob_mpi
    cdef PDM_MPI_Comm PDMC   = PDM_MPI_mpi_2_pdm_mpi_comm(&c_comm)
    # ::::::::::::::::::::::::::::::::::::::::::::::::::

    PDM_compute_dface_normal(<int         *> dface_vtx_idx.data,
                             <PDM_g_num_t *> dface_vtx.data,
                             dn_face,
                             <PDM_g_num_t *> dvtx_distrib.data,
                             <double      *> dvtx_coord.data,
                             <double      *> dface_normal.data,
                             PDMC)

# ===================================================================================
def compute_vtx_characteristic_length(MPI.Comm                                      comm,
                                      int                                           dn_face,
                                      int                                           dn_edge,
                                      int                                           dn_vtx,
                                      NPY.ndarray[NPY.int32_t   , mode='c', ndim=1] dface_vtx_idx,
                                      NPY.ndarray[npy_pdm_gnum_t, mode='c', ndim=1] dface_vtx,
                                      NPY.ndarray[npy_pdm_gnum_t, mode='c', ndim=1] dedge_vtx,
                                      NPY.ndarray[NPY.double_t  , mode='c', ndim=1] dvtx_coord):
    # ::::::::::::::::::::::::::::::::::::::::::::::::::
    # > Convert mpi4py -> PDM_MPI
    cdef MPI.MPI_Comm c_comm = comm.ob_mpi
    cdef PDM_MPI_Comm PDMC   = PDM_MPI_mpi_2_pdm_mpi_comm(&c_comm)
    cdef double *dchar_lenght
    # ::::::::::::::::::::::::::::::::::::::::::::::::::

    cdef PDM_g_num_t* _dedge_vtx
    if dedge_vtx is None:
        _dedge_vtx = NULL
    else:
        _dedge_vtx = <PDM_g_num_t *> dedge_vtx.data

    PDM_compute_vtx_characteristic_length(PDMC,
                                          dn_face,
                                          dn_edge,
                                          dn_vtx,
                          <int         *> dface_vtx_idx.data,
                          <PDM_g_num_t *> dface_vtx.data,
                                          _dedge_vtx,
                          <double      *> dvtx_coord.data,
                                          &dchar_lenght)

    return create_numpy_d(dchar_lenght, dn_vtx)


# ===================================================================================
def compute_entity_distribution(MPI.Comm comm,
                                int      dn_entity):
    # ::::::::::::::::::::::::::::::::::::::::::::::::::
    # > Convert mpi4py -> PDM_MPI
    cdef MPI.MPI_Comm c_comm = comm.ob_mpi
    cdef PDM_MPI_Comm PDMC   = PDM_MPI_mpi_2_pdm_mpi_comm(&c_comm)
    # ::::::::::::::::::::::::::::::::::::::::::::::::::

    n_rank = comm.Get_size()

    cdef PDM_g_num_t* distrib = PDM_compute_entity_distribution(PDMC, dn_entity)
    return create_numpy_g(distrib, n_rank+1)


# ===================================================================================
def dconnectivity_to_extract_dconnectivity(MPI.Comm                                      comm,
                                           NPY.ndarray[npy_pdm_gnum_t, mode='c', ndim=1] select_entity1,
                                           NPY.ndarray[npy_pdm_gnum_t, mode='c', ndim=1] entity1_distribution,
                                           NPY.ndarray[NPY.int32_t   , mode='c', ndim=1] dentity1_entity2_idx,
                                           NPY.ndarray[npy_pdm_gnum_t, mode='c', ndim=1] dentity1_entity2):

    # ::::::::::::::::::::::::::::::::::::::::::::::::::
    # > Convert mpi4py -> PDM_MPI
    cdef MPI.MPI_Comm c_comm = comm.ob_mpi
    cdef PDM_MPI_Comm PDMC   = PDM_MPI_mpi_2_pdm_mpi_comm(&c_comm)
    # ::::::::::::::::::::::::::::::::::::::::::::::::::

    i_rank = comm.Get_rank()
    n_rank = comm.Get_size()
    assert(entity1_distribution.shape[0] == n_rank+1)

    cdef int n_selected_entity1 = select_entity1.shape[0];

    cdef PDM_g_num_t *extract_entity1_distribution
    cdef PDM_g_num_t *extract_entity2_distribution
    cdef int         *dextract_entity1_entity2_idx
    cdef PDM_g_num_t *dextract_entity1_entity2
    cdef PDM_g_num_t *dparent_entity1_g_num
    cdef PDM_g_num_t *dparent_entity2_g_num
    cdef PDM_g_num_t *entity1_old_to_new

    PDM_dconnectivity_to_extract_dconnectivity_bis(PDMC,
                                                   n_selected_entity1,
                          <PDM_g_num_t *>          select_entity1.data,
                          <PDM_g_num_t *>          entity1_distribution.data,
                          <int *>                  dentity1_entity2_idx.data,
                          <PDM_g_num_t *>          dentity1_entity2.data,
                                                   &extract_entity1_distribution,
                                                   &extract_entity2_distribution,
                                                   &dextract_entity1_entity2_idx,
                                                   &dextract_entity1_entity2,
                                                   &dparent_entity1_g_num,
                                                   &dparent_entity2_g_num,
                                                   &entity1_old_to_new);

    dn_extract_entity1 = extract_entity1_distribution[i_rank+1] - extract_entity1_distribution[i_rank]
    dn_extract_entity2 = extract_entity2_distribution[i_rank+1] - extract_entity2_distribution[i_rank]

    np_extract_entity1_distribution = create_numpy_g(extract_entity1_distribution, n_rank+1)
    np_extract_entity2_distribution = create_numpy_g(extract_entity2_distribution, n_rank+1)
    np_dextract_entity1_entity2_idx = create_numpy_i(dextract_entity1_entity2_idx, dn_extract_entity1+1)
    np_dextract_entity1_entity2 = create_numpy_g(dextract_entity1_entity2, dextract_entity1_entity2_idx[dn_extract_entity1])
    np_dparent_entity1_g_num = create_numpy_g(dparent_entity1_g_num, dn_extract_entity1)
    np_dparent_entity2_g_num = create_numpy_g(dparent_entity2_g_num, dn_extract_entity2)
    np_entity1_old_to_new = create_numpy_g(entity1_old_to_new, n_selected_entity1)

    return (np_extract_entity1_distribution, np_extract_entity2_distribution,
            np_dextract_entity1_entity2_idx, np_dextract_entity1_entity2,
            np_dparent_entity1_g_num, np_dparent_entity2_g_num, np_entity1_old_to_new)

def compute_dfacevtx_from_face_and_edge(MPI.Comm                                      comm,
                                        NPY.ndarray[npy_pdm_gnum_t, mode='c', ndim=1] distrib_face,
                                        NPY.ndarray[npy_pdm_gnum_t, mode='c', ndim=1] distrib_edge,
                                        NPY.ndarray[NPY.int32_t   , mode='c', ndim=1] dface_edge_idx,
                                        NPY.ndarray[npy_pdm_gnum_t, mode='c', ndim=1] dface_edge,
                                        NPY.ndarray[npy_pdm_gnum_t, mode='c', ndim=1] dedge_vtx):

    cdef MPI.MPI_Comm c_comm = comm.ob_mpi
    cdef PDM_MPI_Comm PDMC   = PDM_MPI_mpi_2_pdm_mpi_comm(&c_comm)

    cdef PDM_g_num_t* dface_vtx = NULL
    PDM_dconnectivity_dface_vtx_from_face_and_edge(PDMC,
                                 <PDM_g_num_t*>    distrib_face.data,
                                 <PDM_g_num_t*>    distrib_edge.data,
                                 <int        *>    dface_edge_idx.data,
                                 <PDM_g_num_t*>    dface_edge.data,
                                 <PDM_g_num_t*>    dedge_vtx.data,
                                                  &dface_vtx)
    return create_numpy_g(dface_vtx, dface_edge_idx[-1])
