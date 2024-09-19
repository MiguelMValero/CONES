cdef extern from "pdm_io.h":

  ctypedef enum PDM_io_kind_t:
  
    PDM_IO_KIND_MPIIO_EO   = 0
    PDM_IO_KIND_MPIIO_IP   = 1
    PDM_IO_KIND_MPI_SIMPLE = 2
    PDM_IO_KIND_SEQ        = 3

cdef extern from "pdm_mesh_nodal.h":

  ctypedef struct PDM_Mesh_nodal_t:
    pass

cdef extern from "pdm_writer.h":
  # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  # > Wrapping of Ppart Structure
  ctypedef struct PDM_writer_t:
    pass
  # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

  # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  # > Wrapping for C


  ctypedef enum PDM_writer_status_t:
    PDM_WRITER_OFF = 0
    PDM_WRITER_ON  = 1


  ctypedef enum PDM_writer_elt_geom_t:
    PDM_WRITER_POINT    = PDM_MESH_NODAL_POINT
    PDM_WRITER_BAR2     = PDM_MESH_NODAL_BAR2
    PDM_WRITER_TRIA3    = PDM_MESH_NODAL_TRIA3
    PDM_WRITER_QUAD4    = PDM_MESH_NODAL_QUAD4
    PDM_WRITER_POLY_2D  = PDM_MESH_NODAL_POLY_2D
    PDM_WRITER_TETRA4   = PDM_MESH_NODAL_TETRA4
    PDM_WRITER_PYRAMID5 = PDM_MESH_NODAL_PYRAMID5
    PDM_WRITER_PRISM6   = PDM_MESH_NODAL_PRISM6
    PDM_WRITER_HEXA8    = PDM_MESH_NODAL_HEXA8
    PDM_WRITER_POLY_3D  = PDM_MESH_NODAL_POLY_3D


  ctypedef enum PDM_writer_fmt_fic_t: 
    PDM_WRITER_FMT_BIN   = 0
    PDM_WRITER_FMT_ASCII = 1



  ctypedef enum PDM_writer_var_dim_t:
    PDM_WRITER_VAR_CST           = 0
    PDM_WRITER_VAR_SCALAR        = 1
    PDM_WRITER_VAR_VECTOR        = 3
    PDM_WRITER_VAR_TENSOR_SYM    = 6
    PDM_WRITER_VAR_TENSOR_ASYM   = 9

  ctypedef enum PDM_writer_var_loc_t:
    PDM_WRITER_VAR_VERTICES     = 0
    PDM_WRITER_VAR_ELEMENTS     = 1
    PDM_WRITER_VAR_PARTICLES    = 2

  ctypedef enum PDM_writer_topology_t:
    PDM_WRITER_TOPO_CST         = 0
    PDM_WRITER_TOPO_DEFORMABLE  = 1
    PDM_WRITER_TOPO_VARIABLE    = 2


  # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

  # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  
  PDM_writer_t * PDM_writer_create(char                   *fmt,
                                   PDM_writer_fmt_fic_t    fmt_fic,
                                   PDM_writer_topology_t  topologie,
                                   PDM_writer_status_t     st_reprise,
                                   char                   *rep_sortie,
                                   char                   *nom_sortie,
                                   PDM_MPI_Comm            pdm_mpi_comm,
                                   PDM_io_kind_t          acces,
                                   double                  prop_noeuds_actifs,
                                   char                   *options)
  # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

  void PDM_writer_free(PDM_writer_t *cs)

  void PDM_writer_step_beg(PDM_writer_t  *cs, double   physical_time)

  void PDM_writer_step_end(PDM_writer_t  *cs)

  int PDM_writer_geom_create(PDM_writer_t               *cs,
                             char                 *nom_geom,
                             int                   n_part)

  int PDM_writer_geom_create_from_mesh_nodal(PDM_writer_t        *cs,
                                            char                *nom_geom,
                                            PDM_Mesh_nodal_t    *mesh)


  void PDM_writer_geom_coord_set(PDM_writer_t   *cs,
                                 int             id_geom,
                                 int             id_part,
                                 int             n_som,
                                 PDM_real_t      *coords,
                                 PDM_g_num_t     *numabs,
                                 PDM_ownership_t owner)


  void PDM_writer_geom_coord_from_parent_set(PDM_writer_t *cs,
                                             int           id_geom,
                                             int           id_part,
                                             int           n_som,
                                             int           n_som_parent,
                                             PDM_g_num_t  *numabs,
                                             int          *num_parent,
                                             PDM_real_t   *coords_parent,
                                             PDM_g_num_t  *numabs_parent)


  int PDM_writer_geom_bloc_add(PDM_writer_t                *cs,
                               int                    id_geom,
                               PDM_writer_elt_geom_t  t_elt,
                               PDM_ownership_t        owner)

  void PDM_writer_geom_bloc_std_set(PDM_writer_t  *cs,
                                    int      id_geom,
                                    int      id_bloc,
                                    int      id_part,
                                    int      n_elt,
                                    PDM_l_num_t   *connec,
                                    PDM_g_num_t   *numabs)
  

  void PDM_writer_geom_bloc_poly2d_set(PDM_writer_t  *cs,
                                       int            id_geom,
                                       int            id_bloc,
                                       int            id_part,
                                       PDM_l_num_t    n_elt,
                                       PDM_l_num_t   *connec_idx,
                                       PDM_l_num_t   *connec,
                                       PDM_g_num_t   *numabs)
  
  void PDM_writer_geom_bloc_poly3d_set(PDM_writer_t  *cs,
                                       int            id_geom,
                                       int            id_bloc,
                                       int            id_part,
                                       PDM_l_num_t    n_elt,
                                       PDM_l_num_t    n_face,
                                       PDM_l_num_t   *facsom_idx,
                                       PDM_l_num_t   *facsom,
                                       PDM_l_num_t   *cellfac_idx,
                                       PDM_l_num_t   *cellfac,
                                       PDM_g_num_t   *numabs)
  
  void PDM_writer_geom_cell3d_cellface_add(PDM_writer_t *cs,
                                           int           id_geom,
                                           int           id_part,
                                           int           n_cell,
                                           int           n_face,
                                           PDM_l_num_t  *face_som_idx,
                                           PDM_l_num_t  *face_som_nb,
                                           PDM_l_num_t  *face_som,
                                           PDM_l_num_t  *cell_face_idx,
                                           PDM_l_num_t  *cell_face_nb,
                                           PDM_l_num_t  *cell_face,
                                           PDM_g_num_t  *numabs)
  

  void PDM_writer_geom_cell2d_cellface_add(PDM_writer_t *cs,
                                           int           id_geom,
                                           int           id_part,
                                           int           n_cell,
                                           int           n_face,
                                           PDM_l_num_t  *face_som_idx,
                                           PDM_l_num_t  *face_som_nb,
                                           PDM_l_num_t  *face_som,
                                           PDM_l_num_t  *cell_face_idx,
                                           PDM_l_num_t  *cell_face_nb,
                                           PDM_l_num_t  *cell_face,
                                           PDM_g_num_t  *numabs)

  
  void PDM_writer_geom_faces_facesom_add(PDM_writer_t *cs,
                                         int     id_geom,
                                         int     id_part,
                                         int     n_face,
                                         PDM_l_num_t  *face_som_idx,
                                         PDM_l_num_t  *face_som_nb,
                                         PDM_l_num_t  *face_som,
                                         PDM_g_num_t  *numabs)

  
  void PDM_writer_geom_write(PDM_writer_t *cs, int id_geom)

  
  void PDM_writer_geom_data_free(PDM_writer_t *cs, int id_geom)
  

  void PDM_writer_geom_free(PDM_writer_t *cs, int id_geom)

  
  int PDM_writer_var_create(PDM_writer_t         *cs,
                            PDM_writer_status_t   st_dep_tps,
                            PDM_writer_var_dim_t  dim,
                            PDM_writer_var_loc_t  loc,
                            char                 *nom_var)
  
  void PDM_writer_name_map_add(PDM_writer_t *cs,
                               char   *public_name,
                               char   *private_name)
  
  void PDM_writer_var_write(PDM_writer_t *cs, int id_var)
  
  void PDM_writer_var_set(PDM_writer_t  *cs,
                          int           id_var,
                          int           id_geom,
                          int           id_part,
                          PDM_real_t   *val)
  
  void PDM_writer_var_data_free(PDM_writer_t *cs,
                                int           id_var)
  
  void PDM_writer_var_free(PDM_writer_t *cs, int id_var)
  
  # void PDM_writer_fmt_add(char                  *name,            
  #                         PDM_writer_fct_t       create_fct,      
  #                         PDM_writer_fct_t       free_fct,        
  #                         PDM_writer_fct_t       beg_step_fct,    
  #                         PDM_writer_fct_t       end_step_fct,    
  #                         PDM_writer_geom_fct_t  geom_create_fct, 
  #                         PDM_writer_geom_fct_t  geom_write_fct,  
  #                         PDM_writer_geom_fct_t  geom_free_fct,   
  #                         PDM_writer_var_fct_t   var_create_fct,  
  #                         PDM_writer_var_fct_t   var_write_fct,   
  #                         PDM_writer_var_fct_t   var_free_fct)
  
  # void PDM_writer_fmt_free ()
  
  void PDM_writer_geom_data_reset(PDM_writer_t *cs, int     id_geom)
 
# ------------------------------------------------------------------

cdef class Writer:
  """
  """

  cdef PDM_writer_t* _wt

  def __cinit__(self,                  char  *fmt,
                PDM_writer_fmt_fic_t   fmt_fic,
                PDM_writer_topology_t  topologie,
                PDM_writer_status_t    st_reprise,
                char                  *rep_sortie,
                char                  *nom_sortie,
                MPI.Comm               comm,
                PDM_io_kind_t         acces,
                NPY.double_t           prop_noeuds_actifs,
                char                  *options):
      """
      """

      cdef MPI.MPI_Comm c_comm = comm.ob_mpi
      cdef PDM_MPI_Comm PDMC   = PDM_MPI_mpi_2_pdm_mpi_comm(<void *> &c_comm)

      self._wt = PDM_writer_create(fmt,
                                   fmt_fic,
                                   topologie,
                                   st_reprise,
                                   rep_sortie,
                                   nom_sortie,
                                   PDMC,
                                   acces,
                                   prop_noeuds_actifs,
                                   options)

  def geom_create(self,
                  char  *nom_geom,
                  int    n_part):
      """
      """

      return PDM_writer_geom_create(self._wt,
                                    nom_geom,
                                    n_part)


  # def geom_create_from_mesh_nodal(self,   
  #                                 char                *nom_geom,
  #                                 PDM_Mesh_nodal_t    *mesh):
  #     """
  #     """

  #     return PDM_writer_geom_create_from_mesh_nodal(self._wt,
  #                                                   nom_geom,
  #                                                   mesh)

  def geom_cell2d_cellface_add(self,
                               int id_geom,
                               int id_part,
                               NPY.ndarray[NPY.int32_t  , mode='c', ndim=1] edge_vtx,
                               NPY.ndarray[NPY.int32_t  , mode='c', ndim=1] face_edge_idx,
                               NPY.ndarray[NPY.int32_t  , mode='c', ndim=1] face_edge,
                               NPY.ndarray[npy_pdm_gnum_t, mode='c', ndim=1] numabs):
      """
      """
      cdef int n_edge = len(edge_vtx)//2
      cdef int n_face = len(face_edge_idx) - 1

      PDM_writer_geom_cell2d_cellface_add(self._wt,
                                          id_geom,
                                          id_part,
                                          n_face,
                                          n_edge,
                                          NULL,#<int *> edge_vtx_idx.data,
                                          NULL,#<int *> edge_vtx_nb.data,
                                          <int *> edge_vtx.data,
                                          <int *> face_edge_idx.data,
                                          NULL,#<int *> face_edge_nb.data,
                                          <int *> face_edge.data,
                                          <PDM_g_num_t*> numabs.data)

  def geom_cell3d_cellface_add(self,
                               int id_geom,
                               int id_part,
                               NPY.ndarray[NPY.int32_t  , mode='c', ndim=1] face_vtx_idx,
                               NPY.ndarray[NPY.int32_t  , mode='c', ndim=1] face_vtx,
                               NPY.ndarray[NPY.int32_t  , mode='c', ndim=1] cell_face_idx,
                               NPY.ndarray[NPY.int32_t  , mode='c', ndim=1] cell_face,
                               NPY.ndarray[npy_pdm_gnum_t, mode='c', ndim=1] numabs):
      """
      """
      cdef int n_face = len(face_vtx_idx)  - 1
      cdef int n_cell = len(cell_face_idx) - 1

      PDM_writer_geom_cell3d_cellface_add(self._wt,
                                          id_geom,
                                          id_part,
                                          n_cell,
                                          n_face,
                                          <int *> face_vtx_idx.data,
                                          NULL,#<int *> face_vtx_nb.data,
                                          <int *> face_vtx.data,
                                          <int *> cell_face_idx.data,
                                          NULL,#<int *> cell_face_nb.data,
                                          <int *> cell_face.data,
                                          <PDM_g_num_t*> numabs.data)

  def geom_coord_set(self,
                     int            id_geom,
                     int            id_part,
                     # int            n_vtx,
                     NPY.ndarray[NPY.double_t, mode='c', ndim=1] coords,
                     NPY.ndarray[npy_pdm_gnum_t, mode='c', ndim=1] numabs):
      """
      """
      cdef int n_vtx = len(numabs)

      PDM_writer_geom_coord_set(self._wt,
                               id_geom,
                               id_part,
                               n_vtx,
                               <PDM_real_t *> coords.data,
                               <PDM_g_num_t *> numabs.data,
                               PDM_OWNERSHIP_USER)
  

  def geom_faces_facevtx_add(self,
                             int     id_geom,
                             int     id_part,
                             NPY.ndarray[NPY.int32_t  , mode='c', ndim=1] face_vtx_idx,
                             NPY.ndarray[NPY.int32_t  , mode='c', ndim=1] face_vtx,
                             NPY.ndarray[npy_pdm_gnum_t, mode='c', ndim=1] numabs):
      """
      """
      cdef int n_face = len(face_vtx_idx) - 1
      PDM_writer_geom_faces_facesom_add(self._wt,
                                        id_geom,
                                        id_part,
                                        n_face,
                                        <int *> face_vtx_idx.data,
                                        NULL,#<int *> face_vtx_nb.data,
                                        <int *> face_vtx.data,
                                        <PDM_g_num_t*> numabs.data)

  def geom_block_add(self,
                     int                    id_geom,
                     PDM_writer_elt_geom_t  t_elt,
                     PDM_ownership_t        owner):
    """
    """
    return PDM_writer_geom_bloc_add(self._wt, id_geom, t_elt, owner)


  def geom_block_std_set(self,
                         int id_geom,
                         int id_block,
                         int i_part,
                         NPY.ndarray[NPY.int32_t,    mode='c', ndim=1] elt_vtx,
                         NPY.ndarray[npy_pdm_gnum_t, mode='c', ndim=1] elt_ln_to_gn):
    """
    """
    cdef int n_elt = len(elt_ln_to_gn)

    PDM_writer_geom_bloc_std_set(self._wt,
                                 id_geom,
                                 id_block,
                                 i_part,
                                 n_elt,
                 <int         *> elt_vtx.data,
                 <PDM_g_num_t *> elt_ln_to_gn.data)

  def geom_write(self, int id_geom):
      """
      """

      PDM_writer_geom_write(self._wt, id_geom)

  
  def geom_data_free(self, int id_geom):
      """
      """
  
      PDM_writer_geom_data_free(self._wt, id_geom)


  def geom_free(self, int id_geom):
      """
      """
      PDM_writer_geom_free(self._wt, id_geom)


  def var_create(self,
                 PDM_writer_status_t   st_dep_tps,
                 PDM_writer_var_dim_t  dim,
                 PDM_writer_var_loc_t  loc,
                 char                 *nom_var):
      """
      """
  
      return PDM_writer_var_create(self._wt,
                                   st_dep_tps,
                                   dim,
                                   loc,
                                   nom_var)

  def name_map_add(self,
                   char   *public_name,
                   char   *private_name):
      """
      """
  
      PDM_writer_name_map_add(self._wt,
                              public_name,
                              private_name)

  def var_write(self, int id_var):
      """
      """
  
      PDM_writer_var_write(self._wt, id_var)


  def var_set(self,
              int           id_var,
              int           id_geom,
              int           id_part,
              NPY.ndarray[NPY.double_t  , mode='c', ndim=1] val):
      """
      """
  
      PDM_writer_var_set(self._wt,
                         id_var,
                         id_geom,
                         id_part,
                         <double *> val.data)

  def var_data_free(self,
                    int           id_var):
  
      """
      """
      PDM_writer_var_data_free(self._wt,
                              id_var)


  def var_free(self, int id_var):
      """
      """

      PDM_writer_var_free(self._wt, id_var)

  def step_beg(self, double t):
    PDM_writer_step_beg(self._wt, t)

  def step_end(self):
    PDM_writer_step_end(self._wt);

  # ------------------------------------------------------------------------

  def __dealloc__(self):
    """
       Use the free method of PDM Lib
    """
    # ************************************************************************
    # > Declaration
    # ************************************************************************

    PDM_writer_free(self._wt)


def writer_wrapper(comm,
                   directory,
                   name,
                   vtx_coord,
                   vtx_ln_to_gn,
                   elt_vtx_idx,
                   elt_vtx,
                   elt_ln_to_gn,
                   cell_t=-1,
                   cell_face_idx=None,
                   cell_face=None,
                   format="Ensight",
                   elt_fields=None,
                   vtx_fields=None):
  """
  writer_ez(comm, directory, name, vtx_coord, vtx_ln_to_gn, face_vtx_idx, face_vtx, elt_ln_to_gn, cell_face_idx=None, cell_face=None, format="Ensight", elt_fields=None, vtx_fields=None)

  High-level wrapping for Writer

  Parameters:
    comm          (MPI.Comm)                        : MPI Communicator
    directory     (str)                             : Output files directory
    name          (str)                             : Output files prefix
    vtx_coord     (np.ndarray[np.double_t] list)    : Vertex coordinates
    vtx_ln_to_gn  (np.ndarray[npy_pdm_gnum_t] list) : Vertex global ids
    elt_vtx_idx   (np.ndarray[np.int32_t] list)     : Index for elt->vertex connectivity
    elt_vtx       (np.ndarray[np.int32_t] list)     : Elt->vertex connectivity
    elt_ln_to_gn  (np.ndarray[npy_pdm_gnum_t] list) : Element global ids
    cell_t        (int)                             : Type of cell
    cell_face_idx (np.ndarray[np.int32_t] list)     : Index for cell->face connectivity (only for volume meshes)
    cell_face     (np.ndarray[np.int32_t] list)     : Cell->face connectivity (only for volume meshes)
    format        (str)                             : Output file format (default : "Ensight")
    elt_fields    (dict)                            : Element-based fields (optional)
    vtx_fields    (dict)                            : Vertex-based fields (optional)
  """

  is_3d_nodal = (cell_t != -1)

  is_2d = ((cell_face_idx is None) or (cell_face is None)) and (not is_3d_nodal)

  wrt = Writer(format,
               _PDM_WRITER_FMT_ASCII,
               _PDM_WRITER_TOPO_VARIABLE,
               _PDM_WRITER_OFF,
               directory,
               name,
               comm,
               _PDM_IO_KIND_MPI_SIMPLE,
               1.,
               "")

  n_part = len(elt_ln_to_gn)
  n_elt  = [len(g) for g in elt_ln_to_gn]

  id_geom = wrt.geom_create(name,
                            n_part)

  id_var_part = wrt.var_create(_PDM_WRITER_OFF,
                               _PDM_WRITER_VAR_SCALAR,
                               _PDM_WRITER_VAR_ELEMENTS,
                               "i_part")

  id_var_elt_gnum = wrt.var_create(_PDM_WRITER_OFF,
                                   _PDM_WRITER_VAR_SCALAR,
                                   _PDM_WRITER_VAR_ELEMENTS,
                                   "elt_gnum")

  id_var_elt_fields = dict()
  if elt_fields is not None:
    for name in elt_fields:
      id_var_elt_fields[name] = wrt.var_create(_PDM_WRITER_OFF,
                                               _PDM_WRITER_VAR_SCALAR,
                                               _PDM_WRITER_VAR_ELEMENTS,
                                               name)

  id_var_vtx_fields = dict()
  if vtx_fields is not None:
    for name in vtx_fields:
      id_var_vtx_fields[name] = wrt.var_create(_PDM_WRITER_OFF,
                                               _PDM_WRITER_VAR_SCALAR,
                                               _PDM_WRITER_VAR_VERTICES,
                                               name)

  wrt.step_beg(0.)

  for i_part in range(n_part):
    wrt.geom_coord_set(id_geom,
                       i_part,
                       vtx_coord   [i_part],
                       vtx_ln_to_gn[i_part])

    if is_2d:
      wrt.geom_faces_facevtx_add(id_geom,
                                 i_part,
                                 elt_vtx_idx[i_part],
                                 elt_vtx[i_part],
                                 elt_ln_to_gn[i_part])
    else:
      if is_3d_nodal:
        # mono-section
        id_bloc = wrt.geom_block_add(id_geom,
                                     cell_t,
                                     PDM_OWNERSHIP_USER)
        wrt.geom_block_std_set(id_geom,
                               id_bloc,
                               i_part,
                               elt_vtx[i_part],
                               elt_ln_to_gn[i_part])
      else:
        wrt.geom_cell3d_cellface_add(id_geom,
                                     i_part,
                                     elt_vtx_idx[i_part],
                                     elt_vtx[i_part],
                                     cell_face_idx[i_part],
                                     cell_face[i_part],
                                     elt_ln_to_gn[i_part])


  wrt.geom_write(id_geom)

  for i_part in range(n_part):
    wrt.var_set(id_var_part,
                id_geom,
                i_part,
                (comm.rank*n_part + i_part)*NPY.ones(n_elt[i_part], dtype=NPY.double))

    wrt.var_set(id_var_elt_gnum,
                id_geom,
                i_part,
                elt_ln_to_gn[i_part].astype(NPY.double))

    if elt_fields is not None:
      for name in elt_fields:
        wrt.var_set(id_var_elt_fields[name],
                    id_geom,
                    i_part,
                    elt_fields[name][i_part].astype(NPY.double))

    if vtx_fields is not None:
      for name in vtx_fields:
        wrt.var_set(id_var_vtx_fields[name],
                    id_geom,
                    i_part,
                    vtx_fields[name][i_part].astype(NPY.double))


  wrt.var_write(id_var_part)
  wrt.var_write(id_var_elt_gnum)

  if elt_fields is not None:
    for name in elt_fields:
      wrt.var_write(id_var_elt_fields[name])

  if vtx_fields is not None:
    for name in vtx_fields:
      wrt.var_write(id_var_vtx_fields[name])

  wrt.step_end()
