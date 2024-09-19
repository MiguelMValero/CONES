# =======================================================================================
# ---------------------------------------------------------------------------------------
cdef extern from "pdm_geom_elem.h":
    ctypedef struct PDM_geom_t:
      pass

    void PDM_geom_elem_polygon_properties(int      n_polygon,
                                          int     *connect_idx,
                                          int     *connect,
                                          double  *coordinates,
                                          double  *surface_vector,
                                          double  *center,
                                          double  *char_length,
                                          int     *is_degenerated);

# ---------------------------------------------------------------------------------------
# =======================================================================================



# =======================================================================================
# ---------------------------------------------------------------------------------------
def part_geom_elem_polygon_prop(              int  n_polygon, 
    NPY.ndarray[NPY.int32_t   , mode='c', ndim=1]  connect_idx,
    NPY.ndarray[NPY.int32_t   , mode='c', ndim=1]  connect,
    NPY.ndarray[NPY.double_t  , mode='c', ndim=1]  coordinates,
    NPY.ndarray[NPY.double_t  , mode='c', ndim=1]  surface_vector,
    NPY.ndarray[NPY.double_t  , mode='c', ndim=1]  center,
    NPY.ndarray[NPY.double_t  , mode='c', ndim=1]  char_length,
    NPY.ndarray[NPY.int32_t   , mode='c', ndim=1]  is_degenerated,):
    '''
    '''

    # print("[CYTHON] pdm_geom.pxi::call to PDM_geom_elem_polygon_properties")
    PDM_geom_elem_polygon_properties( n_polygon,
                      <int         *> connect_idx   .data,
                      <int         *> connect       .data,
                      <double      *> coordinates   .data,
                      <double      *> surface_vector.data,
                      <double      *> center        .data,
                      <double      *> char_length   .data,
                      <int         *> is_degenerated.data)
    
    # print("[CYTHON] pdm_geom.pxi::convert PDM_geom_elem_polygon_properties result to NUMPY")
    # np_surface_vector = create_numpy_d(surface_vector, n_polygon)
    # np_center         = create_numpy_d(center        , n_polygon)

    # return np_surface_vector, np_center
    # print("[CYTHON] pdm_geom.pxi::return")
    return surface_vector, center

 

# ---------------------------------------------------------------------------------------
# =======================================================================================