cdef extern from "pdm_global_mean.h":
  ctypedef struct PDM_global_point_mean_t:
    pass

  PDM_global_point_mean_t* PDM_global_mean_create(const int n_part, PDM_MPI_Comm comm)

  void PDM_global_mean_set(PDM_global_point_mean_t *gmean,
                           const int                i_part,
                           const int                n_point,
                           const PDM_g_num_t       *numabs)

  void PDM_global_mean_field_set(PDM_global_point_mean_t *gmean,
                                 const int                i_part,
                                 const int                stride,
                                 const double            *local_field,
                                 const double            *local_weight,
                                 double                  *global_mean_field_ptr)

  void PDM_global_mean_field_compute(PDM_global_point_mean_t *gmean)

  void PDM_global_mean_free(PDM_global_point_mean_t *gmean)




# ------------------------------------------------------------------
cdef class GlobalMean:
  """ Interface for pdm_global_mean.c """
  # ************************************************************************
  cdef PDM_global_point_mean_t* global_mean
  cdef list                     gnum_list # KEEP ALIVE THESE ARRAYS
  cdef list                     n_elem
  cdef int                      n_part
  # ************************************************************************

  # ------------------------------------------------------------------
  def __init__(self, list gnum_list, MPI.Comm comm):

    # Checks
    self.n_part = len(gnum_list)
    for i in range(self.n_part):
      assert_single_dim_np(gnum_list[i], npy_pdm_gnum_dtype)

    self.n_elem = [gnum.size for gnum in gnum_list]
    self.gnum_list = gnum_list # Keep alive until gmean is freed

    # Convert input data
    cdef MPI.MPI_Comm c_comm = comm.ob_mpi
    cdef PDM_MPI_Comm PDMC   = PDM_MPI_mpi_2_pdm_mpi_comm(&c_comm)

    self.global_mean = PDM_global_mean_create(self.n_part, PDMC)

    cdef NPY.ndarray[npy_pdm_gnum_t, ndim=1, mode='c'] np_gnum
    for i_part in range(self.n_part):
      np_gnum = gnum_list[i_part]
      PDM_global_mean_set(self.global_mean,
                          i_part,
                    <int> np_gnum.size, 
           <PDM_g_num_t*> np_gnum.data)


  # ------------------------------------------------------------------
  def compute_field(self, list field_list, list weight_list, int stride=1):
    
    # Declarations
    cdef NPY.ndarray[double, ndim=1, mode='c'] np_field
    cdef NPY.ndarray[double, ndim=1, mode='c'] np_weight
    cdef NPY.ndarray[double, ndim=1, mode='c'] np_global_field

    # Checks
    assert len(field_list) == self.n_part
    for i in range(self.n_part):
      assert_single_dim_np(field_list[i],  NPY.float64, stride*self.n_elem[i])
      assert_single_dim_np(weight_list[i], NPY.float64, self.n_elem[i])

    # Allocate output
    output = [NPY.empty(field.size, NPY.float64) for field in field_list]

    for i_part in range(self.n_part):
      np_field  = field_list[i_part]    #In
      np_weight = weight_list[i_part]   #In
      np_global_field = output[i_part]  #Out

      PDM_global_mean_field_set(self.global_mean,
                                i_part,
                                stride,
                     <double*>  np_field.data,
                     <double*>  np_weight.data,
                     <double*>  np_global_field.data)


    PDM_global_mean_field_compute(self.global_mean)

    return output

  # ------------------------------------------------------------------
  def __dealloc__(self):
    PDM_global_mean_free(self.global_mean)
