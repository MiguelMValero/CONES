
cdef extern from "pdm_part_domain_interface.h":
    # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    # > Wrapping of Ppart Structure
    ctypedef struct PDM_part_domain_interface_t:
      pass
    # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    # > Wrapping of function
    PDM_part_domain_interface_t* PDM_part_domain_interface_create(int                          n_interface,
                                                                  int                          n_domain,
                                                                  int                         *n_part,
                                                                  PDM_domain_interface_mult_t  multidomain_interface,
                                                                  PDM_ownership_t              ownership,
                                                                  PDM_MPI_Comm                 comm);

    void PDM_part_domain_interface_set(PDM_part_domain_interface_t  *dom_intrf,
                                       PDM_bound_type_t              interface_kind,
                                       int                           i_domain,
                                       int                           i_part,
                                       int                           i_interface,
                                       int                           interface_pn,
                                       PDM_g_num_t                  *interface_ln_to_gn,
                                       int                          *interface_sgn,
                                       int                          *interface_sens,
                                       int                          *interface_ids,
                                       int                          *interface_ids_idx,
                                       int                          *interface_dom);

    void PDM_part_domain_interface_get(PDM_part_domain_interface_t   *dom_intrf,
                                       PDM_bound_type_t               interface_kind,
                                       int                            i_domain,
                                       int                            i_part,
                                       int                            i_interface,
                                       int                           *interface_pn,
                                       PDM_g_num_t                  **interface_ln_to_gn,
                                       int                          **interface_sgn,
                                       int                          **interface_sens,
                                       int                          **interface_ids,
                                       int                          **interface_ids_idx,
                                       int                          **interface_dom);

    int PDM_part_domain_interface_n_interface_get(PDM_part_domain_interface_t   *dom_intrf);

    int PDM_part_domain_interface_exist_get(PDM_part_domain_interface_t  *dom_intrf,
                                            PDM_bound_type_t              interface_kind);

    void PDM_part_domain_interface_free(PDM_part_domain_interface_t  *dom_intrf);
    void PDM_part_domain_interface_translation_set(PDM_part_domain_interface_t  *dom_intrf,
                                                   int                           i_interface,
                                                   double                       *vect);

    void PDM_part_domain_interface_rotation_set(PDM_part_domain_interface_t  *dom_intrf,
                                                int                           i_interface,
                                                double                       *direction,
                                                double                       *center,
                                                double                        angle);
    void PDM_part_domain_interface_translation_get(PDM_part_domain_interface_t  *dom_intrf,
                                                   int                           i_interface,
                                                   double                      **vect);

    void PDM_part_domain_interface_rotation_get(PDM_part_domain_interface_t  *dom_intrf,
                                                int                           i_interface,
                                                double                      **direction,
                                                double                      **center,
                                                double                       *angle);


cdef class PartDomainInterface:
  """
  """
  # ************************************************************************
  # > Class attributes
  cdef PDM_part_domain_interface_t *pdi
  cdef object keep_alive
  # ************************************************************************
  def __cinit__(self,
                int                                           n_interface,
                int                                           n_domain,
                NPY.ndarray[NPY.int32_t   , mode='c', ndim=1] n_part,
                bint                                          multidomain_interface,
                MPI.Comm                                      comm):
    """
    """
    cdef MPI.MPI_Comm c_comm = comm.ob_mpi
    cdef PDM_MPI_Comm PDMC   = PDM_MPI_mpi_2_pdm_mpi_comm(<void *> &c_comm)
    self.keep_alive = list()

    cdef _multidomain_interface = PDM_DOMAIN_INTERFACE_MULT_YES if multidomain_interface else PDM_DOMAIN_INTERFACE_MULT_NO
    self.pdi = PDM_part_domain_interface_create(n_interface,
                                                n_domain,
                                     <int*>     n_part.data,
                                                _multidomain_interface,
                                                PDM_OWNERSHIP_USER,
                                                PDMC)

  # ------------------------------------------------------------------------
  def interface_set(self,
                    PDM_bound_type_t                              interface_kind,
                    int                                           i_domain,
                    int                                           i_part,
                    int                                           i_interface,
                    int                                           interface_pn,
                    NPY.ndarray[npy_pdm_gnum_t, mode='c', ndim=1] interface_ln_to_gn,
                    NPY.ndarray[NPY.int32_t   , mode='c', ndim=1] interface_sgn,
                    NPY.ndarray[NPY.int32_t   , mode='c', ndim=1] interface_sens,
                    NPY.ndarray[NPY.int32_t   , mode='c', ndim=1] interface_ids,
                    NPY.ndarray[NPY.int32_t   , mode='c', ndim=1] interface_ids_idx,
                    NPY.ndarray[NPY.int32_t   , mode='c', ndim=1] interface_dom):
    """
    """
    self.keep_alive.append(interface_ln_to_gn)
    self.keep_alive.append(interface_sgn)
    self.keep_alive.append(interface_sens)
    self.keep_alive.append(interface_ids)
    self.keep_alive.append(interface_ids_idx)
    self.keep_alive.append(interface_dom)
    PDM_part_domain_interface_set(self.pdi,
                                  interface_kind,
                                  i_domain,
                                  i_part,
                                  i_interface,
                                  interface_pn,
                  <PDM_g_num_t *> interface_ln_to_gn.data,
                  <int         *> interface_sgn.data,
                  <int         *> interface_sens.data,
                  <int         *> interface_ids.data,
                  <int         *> interface_ids_idx.data,
                  <int         *> interface_dom.data);

  # ------------------------------------------------------------------------
  def translation_set(self,
                      i_interface,
                      NPY.ndarray[NPY.double_t  , mode='c', ndim=1] vect):
    """
    """
    PDM_part_domain_interface_translation_set(self.pdi,
                                              i_interface,
                                   <double *> vect.data)

  # ------------------------------------------------------------------------
  def rotation_set(self,
                   i_interface,
                   NPY.ndarray[NPY.double_t  , mode='c', ndim=1] direction,
                   NPY.ndarray[NPY.double_t  , mode='c', ndim=1] center,
                   double                                        angle):
    """
    """
    PDM_part_domain_interface_rotation_set(self.pdi,
                                           i_interface,
                                <double *> direction.data,
                                <double *> center.data,
                                           angle)

  # ------------------------------------------------------------------------
  def __dealloc__(self):
    """
    """
    PDM_part_domain_interface_free(self.pdi)
