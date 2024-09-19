cdef extern from "pdm_vtk.h":

  void PDM_vtk_write_polydata(char*        filename,
                              int          n_vtx,
                              double*      vtx_coord,
                              PDM_g_num_t* vtx_g_num,
                              int          n_face,
                              int*         face_vtx_idx,
                              int*         face_vtx,
                              PDM_g_num_t* face_g_num,
                              int*         face_color)

  void PDM_vtk_write_point_cloud(char*        filename,
                                 int          n_vtx,
                                 double*      vtx_coord,
                                 PDM_g_num_t* vtx_g_num,
                                 int*         color)

# ------------------------------------------------------------------

def vtk_write_polydata(char *                                           filename,
                       NPY.ndarray[NPY.double_t, mode='c', ndim=1]      vtx_coord,
                       NPY.ndarray[npy_pdm_gnum_t  , mode='c', ndim=1]  vtx_g_num,
                       NPY.ndarray[NPY.int32_t  , mode='c', ndim=1]     face_vtx_idx,
                       NPY.ndarray[NPY.int32_t  , mode='c', ndim=1]     face_vtx,
                       NPY.ndarray[npy_pdm_gnum_t, mode='c', ndim=1]    face_g_num,
                       NPY.ndarray[NPY.int32_t  , mode='c', ndim=1]     face_color):
    """
    vtk_write_polydata(filename, n_vtx, vtx_coord, vtx_g_num, n_face, face_vtx_idx, face_vtx, face_g_num, face_color)
    Write polygons in a file of vtk format.

    Parameters:
      filename        (str)                           :
      n_vtx           (int)                           :
      vtx_coord       (Numpy array of NPY.double_t)   :
      vtx_g_num       (Numpy array of npy_pdm_gnum_t) :
      n_face          (int)                           :
      face_vtx_idx    (Numpy array of NPY.int32_t)    :
      face_vtx        (Numpy array of  NPY.int32_t)   :
      face_g_num      (Numpy array of npy_pdm_gnum_t) :
      face_color      (Numpy array of  NPY.int32_t)   :
    """

    cdef int n_vtx  = len(vtx_coord) // 3
    cdef int n_face = len(face_vtx_idx) - 1

    cdef PDM_g_num_t *c_vtx_g_num  = NULL
    cdef PDM_g_num_t *c_face_g_num = NULL
    cdef int         *c_face_color = NULL

    if vtx_g_num is not None:
      c_vtx_g_num = <PDM_g_num_t *> vtx_g_num.data

    if face_g_num is not None:
      c_face_g_num = <PDM_g_num_t *> face_g_num.data

    if face_color is not None:
      c_face_color = <int *> face_color.data

    PDM_vtk_write_polydata(filename,
                           n_vtx,
                           <PDM_real_t *> vtx_coord.data,
                           c_vtx_g_num,
                           n_face,
                           <int *> face_vtx_idx.data,
                           <int *> face_vtx.data,
                           c_face_g_num,
                           c_face_color)

def vtk_write_point_cloud(char *                                           filename,
                          NPY.ndarray[NPY.double_t, mode='c', ndim=1]      vtx_coord,
                          NPY.ndarray[npy_pdm_gnum_t  , mode='c', ndim=1]  vtx_g_num,
                          NPY.ndarray[NPY.int32_t  , mode='c', ndim=1]     color):

  """
  vtk_write_point_cloud(filename, n_vtx, vtx_coord, vtx_g_num, color)
  Write point cloud in a file of vtk format.

  Parameters:
    filename        (str)                           :
    n_vtx           (int)                           :
    vtx_coord       (Numpy array of NPY.double_t)   :
    vtx_g_num       (Numpy array of npy_pdm_gnum_t) :
    color           (Numpy array of  NPY.int32_t)   :
  """

  cdef int n_vtx = len(vtx_coord) // 3

  cdef PDM_g_num_t *c_vtx_g_num = NULL
  cdef int         *c_color     = NULL

  if vtx_g_num is not None:
    c_vtx_g_num = <PDM_g_num_t *> vtx_g_num.data

  if color is not None:
    c_color = <int *> color.data

  PDM_vtk_write_point_cloud(filename,
                            n_vtx,
                            <PDM_real_t *> vtx_coord.data,
                            c_vtx_g_num,
                            c_color)
