#include "cwipi_configf.h"

program fortran_new_api_polygon_sol

    use cwp

    implicit none

    include "mpif.h"

  !--------------------------------------------------------------------
  logical, parameter                      :: verbose = .false.
  logical                                 :: bonus = .false.
  integer                                 :: i
  character(999)                          :: arg

  integer(c_int)                          :: spatial_interp_algorithm
  integer(c_int)                          :: location
  integer                                 :: array_size

  integer, parameter                      :: n_vtx = 11, n_elts = 5

  integer                                 :: ierr
  integer                                 :: i_rank, n_rank

  integer                                 :: n_code
  character(len = 5),            pointer  :: code_names(:)         => null()
  integer                                 :: is_active_rank = CWP_STATUS_ON
  integer,                       pointer  :: intra_comms(:)        => null()

  integer                                 :: n_part
  character(len = 5),            pointer  :: coupled_code_names(:) => null()
  character(len = 99)                     :: coupling_name

  double precision, pointer, dimension(:,:) :: coords => null()
  integer(c_long), pointer, dimension(:)  :: vtx_g_num => null()

  integer, pointer, dimension(:)          :: connec_idx => null()
  integer, pointer, dimension(:)          :: connec => null()
  integer(c_long), pointer, dimension(:)  :: elt_g_num  => null()
  integer(c_int)                          :: id_block

  character(len = 99)                     :: field_name
  integer(c_int)                          :: n_components = 1

  double precision,              pointer  :: field_data(:) => null()
  !--------------------------------------------------------------------

  ! Read command line arguments
  i = 1
  do while (i <= command_argument_count())
    call get_command_argument(i, arg)
    select case(arg)
      case ("-b")
        bonus = .true.
      case default
        print *, "Invalid command argument ", arg
        stop
    end select
    i = i + 1
  enddo

  ! MPI Initialization
  call MPI_Init(ierr)
  call MPI_Comm_rank(mpi_comm_world, i_rank, ierr)
  call MPI_Comm_size(mpi_comm_world, n_rank, ierr)

  n_code = 1

  allocate(code_names(n_code), &
           intra_comms(n_code))

  code_names(1) = "code2"

  call CWP_Init(mpi_comm_world, &
                n_code,         &
                code_names,     &
                is_active_rank, &
                intra_comms)

  coupling_name = "code1_code2";

  allocate(coupled_code_names(n_code))

  coupled_code_names(1) = "code1"

  n_part = 1

  spatial_interp_algorithm = CWP_SPATIAL_INTERP_FROM_LOCATION_MESH_LOCATION_OCTREE
  if (verbose) then
    spatial_interp_algorithm = CWP_SPATIAL_INTERP_FROM_INTERSECTION
  endif

  call CWP_Cpl_create(code_names(1),                                         &
                      coupling_name,                                         &
                      coupled_code_names(1),                                 &
                      CWP_INTERFACE_SURFACE,                                 &
                      CWP_COMM_PAR_WITH_PART,                                &
                      spatial_interp_algorithm, &
                      n_part,                                                &
                      CWP_DYNAMIC_MESH_STATIC,                               &
                      CWP_TIME_EXCH_USER_CONTROLLED)

  call CWP_Visu_set(code_names(1),           &
                    coupling_name,           &
                    1,                       &
                    CWP_VISU_FORMAT_ENSIGHT, &
                    "text")

  allocate(coords(3,n_vtx))
  coords(:, 1) = [0,0,0]
  coords(:, 2) = [1,0,0]
  coords(:, 3) = [2,0,0]
  coords(:, 4) = [3,0,0]
  coords(:, 5) = [0,1,0]
  coords(:, 6) = [2,1,0]
  coords(:, 7) = [3,1,0]
  coords(:, 8) = [1,2,0]
  coords(:, 9) = [0,3,0]
  coords(:,10) = [2,3,0]
  coords(:,11) = [3,3,0]
  call CWP_Mesh_interf_vtx_set(code_names(1), &
                               coupling_name, &
                               0,             &
                               n_vtx,         &
                               coords,        &
                               vtx_g_num)

  id_block = CWP_Mesh_interf_block_add(code_names(1),       &
                                       coupling_name,       &
                                       CWP_BLOCK_FACE_POLY)

  allocate(connec_idx(n_elts+1))
  connec_idx = [0,3,7,11,16,21]
  allocate(connec(21))
  connec = [1,2,5,   3,4,7,6,   5,8,10,9   ,5,2,3,6,8,   6,7,11,10,8]
  call CWP_Mesh_interf_f_poly_block_set(code_names(1), &
                                        coupling_name, &
                                        0,             &
                                        id_block,      &
                                        n_elts,        &
                                        connec_idx,    &
                                        connec,        &
                                        elt_g_num)

  call CWP_Mesh_interf_finalize(code_names(1), &
                                coupling_name)

  field_name   = "a super fancy field"
  n_components = 1

  location = CWP_DOF_LOCATION_NODE
  if (verbose) then
    location = CWP_DOF_LOCATION_CELL_CENTER
  endif

  call CWP_Field_create(code_names(1),                &
                        coupling_name,                &
                        field_name,                   &
                        CWP_DOUBLE,                   &
                        CWP_FIELD_STORAGE_INTERLACED, &
                        n_components,                 &
                        location,        &
                        CWP_FIELD_EXCH_RECV,          &
                        CWP_STATUS_ON)

  array_size = n_vtx * n_components
  if (verbose) then
    array_size = n_elts * n_components
  endif

  allocate(field_data(array_size))

  call CWP_Field_data_set(code_names(1),        &
                          coupling_name,        &
                          field_name,           &
                          0,                    &
                          CWP_FIELD_MAP_TARGET, &
                          field_data)

  call CWP_Time_step_beg(code_names(1), &
                         0.d0)

  call CWP_Spatial_interp_property_set(code_names(1), &
                                       coupling_name, &
                                       "tolerance",   &
                                       CWP_DOUBLE,    &
                                       "0.001")

  call CWP_Spatial_interp_weights_compute(code_names(1), &
                                          coupling_name)

  call CWP_Field_irecv(code_names(1), &
                       coupling_name, &
                       field_name)

  call CWP_Field_wait_irecv(code_names(1), &
                            coupling_name, &
                            field_name)

 call CWP_Time_step_end(code_names(1))

 call CWP_Field_Del(code_names(1),   &
                    coupling_name,   &
                    field_name)

 call CWP_Mesh_interf_del(code_names(1), &
                          coupling_name)

 call CWP_Cpl_Del(code_names(1), &
                  coupling_name)

 deallocate(coords)
 deallocate(connec)
 deallocate(connec_idx)
 deallocate(field_data)

  ! Finalize CWIPI :
  call CWP_Finalize()

  call MPI_Finalize(ierr)


contains



end program fortran_new_api_polygon_sol
