program fortran_new_api_deformable_sol

    use cwp
    use pdm_generate_mesh

    implicit none

    include "mpif.h"

  !--------------------------------------------------------------------
  integer                                     :: ierr
  integer                                     :: i_rank, n_rank
  integer(kind=8), parameter                  :: n_vtx_seg = 10

  integer                                     :: n_code
  character(len = 99),                pointer :: code_names(:)         => null()
  integer                                     :: is_active_rank = CWP_STATUS_ON
  integer,                            pointer :: intra_comms(:)        => null()

  integer                                     :: n_part
  character(len = 99),                pointer :: coupled_code_names(:) => null()
  character(len = 99)                         :: coupling_name

  double precision,   dimension(:,:), pointer :: coords    => null()
  integer(c_long),    dimension(:),   pointer :: vtx_g_num => null()

  integer,            dimension(:),   pointer :: elt_vtx_idx => null()
  integer,            dimension(:),   pointer :: elt_vtx     => null()
  integer(c_long),    dimension(:),   pointer :: elt_g_num   => null()
  integer(c_int)                              :: id_block

  character(len = 99)                         :: field_name
  integer(c_int)                              :: n_components

  double precision                            :: dt
  double precision                            :: degrad
  double precision                            :: x
  double precision                            :: y
  double precision                            :: alpha
  double precision                            :: sina
  double precision                            :: cosa
  double precision                            :: time
  integer                                     :: i, it, itdeb, itend
  integer                                     :: n_vtx, n_elt

  double precision,                   pointer :: field_data(:) => null()
  !--------------------------------------------------------------------

  ! MPI Initialization :
  call MPI_Init(ierr)
  call MPI_Comm_rank(mpi_comm_world, i_rank, ierr)
  call MPI_Comm_size(mpi_comm_world, n_rank, ierr)

  ! Initialize CWIPI :
  n_code = 1

  allocate(code_names(n_code), &
           intra_comms(n_code))

  code_names(1) = "code2"

  call CWP_Init(mpi_comm_world, &
                n_code,         &
                code_names,     &
                is_active_rank, &
                intra_comms)

  ! Create the coupling :
  ! CWP_DYNAMIC_MESH_DEFORMABLE allows us to take into account the modifications
  ! to the mesh over the coupling steps.
  coupling_name     = "coupling"
  allocate(coupled_code_names(n_code))
  coupled_code_names(1) = "code1"
  n_part = 1;
  call CWP_Cpl_create(code_names(1),                                         &
                      coupling_name,                                         &
                      coupled_code_names(1),                                 &
                      CWP_INTERFACE_SURFACE,                                 &
                      CWP_COMM_PAR_WITH_PART,                                &
                      CWP_SPATIAL_INTERP_FROM_LOCATION_MESH_LOCATION_OCTREE, &
                      n_part,                                                &
                      CWP_DYNAMIC_MESH_DEFORMABLE,                           &
                      CWP_TIME_EXCH_USER_CONTROLLED)

  ! Set coupling visualisation:
  call CWP_Visu_set(code_names(1),           &
                    coupling_name,           &
                    1,                       &
                    CWP_VISU_FORMAT_ENSIGHT, &
                    "text")

  ! Create mesh :
  call PDM_generate_mesh_rectangle_simplified(intra_comms(1), &
                                              n_vtx_seg,      &
                                              n_vtx,          &
                                              n_elt,          &
                                              coords,         &
                                              elt_vtx_idx,    &
                                              elt_vtx)

    call CWP_Mesh_interf_vtx_set(code_names(1), &
                               coupling_name, &
                               0,             &
                               n_vtx,         &
                               coords,        &
                               vtx_g_num)

  id_block = CWP_Mesh_interf_block_add(code_names(1),       &
                                       coupling_name,       &
                                       CWP_BLOCK_FACE_POLY)

  call CWP_Mesh_interf_f_poly_block_set(code_names(1), &
                                        coupling_name, &
                                        0,             &
                                        id_block,      &
                                        n_elt,         &
                                        elt_vtx_idx,   &
                                        elt_vtx,       &
                                        elt_g_num)

  call CWP_Mesh_interf_finalize(code_names(1), &
                                coupling_name)

  field_name   = "a super fancy field"
  n_components = 1

  call CWP_Field_create(code_names(1),                &
                        coupling_name,                &
                        field_name,                   &
                        CWP_DOUBLE,                   &
                        CWP_FIELD_STORAGE_INTERLACED, &
                        n_components,                 &
                        CWP_DOF_LOCATION_NODE,        &
                        CWP_FIELD_EXCH_SEND,          &
                        CWP_STATUS_ON)

  allocate(field_data(n_vtx*n_components))

  do i=1, n_vtx
    field_data(i) = coords(1, i)
  enddo

  call CWP_Field_data_set(code_names(1),        &
                          coupling_name,        &
                          field_name,           &
                          0,                    &
                          CWP_FIELD_MAP_SOURCE, &
                          field_data)

  call CWP_Spatial_interp_property_set(code_names(1), &
                                       coupling_name, &
                                       "tolerance",   &
                                       CWP_DOUBLE,    &
                                       "0.001")

  ! Interations :
  ! At each iteration the mesh coordinates and the exchanged fields are modified.
  itdeb = 1
  itend = 10
  time  = 0.0d0
  dt    = 0.1d0

  degrad = acos(-1.0)/180.
  x = 0.0
  y = 0.0
  alpha = 2
  alpha = alpha * degrad
  sina = sin(alpha)
  cosa = cos(alpha)

  do it = itdeb, itend

    time = (it-itdeb)*dt

    ! Begin time step :
    call CWP_Time_step_beg(code_names(1), &
                           time)

    call CWP_Spatial_interp_weights_compute(code_names(1), &
                                            coupling_name)

    call CWP_Field_issend(code_names(1), &
                          coupling_name, &
                          field_name)

    call CWP_Field_wait_issend(code_names(1), &
                               coupling_name, &
                               field_name)

    call CWP_Time_step_end(code_names(1))

  enddo

  ! Delete field :
  call CWP_Field_Del(code_names(1),   &
                     coupling_name,   &
                     field_name)

  ! Delete Mesh :
  call CWP_Mesh_interf_del(code_names(1), &
                           coupling_name)

  ! Delete the coupling :
  call CWP_Cpl_Del(code_names(1), &
                   coupling_name)

  ! free
  deallocate(code_names)
  deallocate(intra_comms)
  deallocate(coupled_code_names)
  deallocate(field_data)

  call pdm_fortran_free_c(c_loc(coords))
  call pdm_fortran_free_c(c_loc(elt_vtx_idx))
  call pdm_fortran_free_c(c_loc(elt_vtx))

  ! Finalize CWIPI :
  call CWP_Finalize()

  ! Finalize MPI :
  call MPI_Finalize(ierr)

end program fortran_new_api_deformable_sol
