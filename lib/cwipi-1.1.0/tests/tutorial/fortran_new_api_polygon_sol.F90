#include "cwipi_configf.h"

program fortran_new_api_polygon_sol

#ifdef CWP_HAVE_FORTRAN_MPI_MODULE
    use mpi
#endif
    use cwp

    implicit none

#ifndef CWP_HAVE_FORTRAN_MPI_MODULE
    include "mpif.h"
#endif

  !--------------------------------------------------------------------
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
  integer(c_int)                          :: n_components

  integer                                 :: i

  double precision,              pointer  :: send_field_data(:) => null()
  double precision,              pointer  :: recv_field_data(:) => null()

  integer(c_int)                          :: n_uncomputed_tgts
  integer(c_int),                pointer  :: uncomputed_tgts(:) => null()

  logical                                 :: I_am_code1
  !--------------------------------------------------------------------

  ! MPI Initialization
  call MPI_Init(ierr)
  call MPI_Comm_rank(mpi_comm_world, i_rank, ierr)
  call MPI_Comm_size(mpi_comm_world, n_rank, ierr)

  ! This test mimics the coupling between 2 codes running each
  ! on one processor.
  if (n_rank /= 2) then
    print *, " This executable must be executed on 2 processes"
    stop
  endif

  ! Initialize CWIPI :
  ! Here 2 codes are coupled. code1 runs on the processor of
  ! MPI rank 0 and code2 runs on the processor if MPI rank 1.
  ! In this version of CWIPI several codes can execute on the
  ! same MPI rank (here only one code per processor, so n_code = 1).
  ! Therefore, an array of code names is given at initialization.
  ! is_active_rank tells if current ranks will be used
  ! in the CWIPI coupling computations.
  ! intra_comm is an array of MPI communicators
  ! giving the for each code on the processors the communicator
  ! to communicate through the ranks of that code.
  n_code = 1

  allocate(code_names(n_code), &
           intra_comms(n_code))

  ! for code1
  if (i_rank == 0) then
    code_names(1) = "code1"
    I_am_code1    = .true.
  ! for code2
  else
    code_names(1) = "code2"
    I_am_code1    = .false.
  endif

  call CWP_Init(mpi_comm_world, &
                n_code,         &
                code_names,     &
                is_active_rank, &
                intra_comms)

  ! Create the coupling :
  ! One CWIPI context can hold several couplings. Let us set up the
  ! coupling between code1 and code2. CWP_INTERFACE_SURFACE informs
  ! that the geometrical interface of the meshes of the coupled
  ! codes is a surface, still for CWIPI the coordinate system is 3D.
  ! CWP_COMM_PAR_WITH_PART means that each mesh is partitionned
  ! over the processors of its code. Here the mesh does not change
  ! over the coupling, so CWP_DYNAMIC_MESH_STATIC is set.
  ! CWP_TIME_EXCH_USER_CONTROLLED is not used yet.
  coupling_name = "code1_code2";

  allocate(coupled_code_names(n_code))

  ! for code1
  if (I_am_code1) then
    coupled_code_names(1) = "code2"
  ! for code2
  else
    coupled_code_names(1) = "code1"
  endif

  n_part = 1
  call CWP_Cpl_create(code_names(1),                                         &
                      coupling_name,                                         &
                      coupled_code_names(1),                                 &
                      CWP_INTERFACE_SURFACE,                                 &
                      CWP_COMM_PAR_WITH_PART,                                &
                      CWP_SPATIAL_INTERP_FROM_LOCATION_MESH_LOCATION_OCTREE, &
                      n_part,                                                &
                      CWP_DYNAMIC_MESH_STATIC,                               &
                      CWP_TIME_EXCH_USER_CONTROLLED)

  ! Set coupling visualisation:
  ! Output files of the code fields will be written in Ensight ASCII
  ! format (easily readable by paraview) at each iteration.
  call CWP_Visu_set(code_names(1),           &
                    coupling_name,           &
                    1,                       &
                    CWP_VISU_FORMAT_ENSIGHT, &
                    "text")

  ! Create the field
  ! It is possible to operate a bidirectional exchange (see c_new_vs_old_sendrecv).
  ! For sake of simplicity, this example will only send the field
  ! of code1 (CWP_FIELD_EXCH_SEND) to code2 (CWP_FIELD_EXCH_RECV).
  ! On code1 there is a field (CWP_FIELD_MAP_SOURCE) located at
  ! the vertices (CWP_DOF_LOCATION_NODE) with one component (n_components)
  ! which is the x coordinate of the mesh in this test.
  field_name   = "a super fancy field"
  n_components = 1

  ! for code1
  if (I_am_code1) then

    call CWP_Field_create(code_names(1),                &
                          coupling_name,                &
                          field_name,                   &
                          CWP_DOUBLE,                   &
                          CWP_FIELD_STORAGE_INTERLACED, &
                          n_components,                 &
                          CWP_DOF_LOCATION_NODE,        &
                          CWP_FIELD_EXCH_SEND,          &
                          CWP_STATUS_ON)

  ! for code2
  else

    call CWP_Field_create(code_names(1),                &
                          coupling_name,                &
                          field_name,                   &
                          CWP_DOUBLE,                   &
                          CWP_FIELD_STORAGE_INTERLACED, &
                          n_components,                 &
                          CWP_DOF_LOCATION_NODE,        &
                          CWP_FIELD_EXCH_RECV,          &
                          CWP_STATUS_ON)

  endif

  ! Begin time step :
  ! In this example there is only one time step. It is mandatory to create the
  ! coupling and the associated fields before starting the first time step.
  call CWP_Time_step_beg(code_names(1), &
                         0.d0)

  ! Set the mesh vertices coordinates :
  ! The coordinate system in CWIPI is always 3D, so
  ! we allocate an array of the time the number of vertices
  ! (11 here) to set the coordinates in. The coordinates are
  ! interlaced (x0, y0, z0, x1, y1, z1, ..., xn, yn, zn).
  ! vtx_g_num is a null() argument that will be explained later.
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

  ! Set the mesh polygons connectivity :
  ! Let us set a mesh of 5 polygons (CWP_BLOCK_FACE_POLY).
  ! An index array (connec_idx) of size n_elts+1 contains the
  ! information of the number of vertices per polygon. The first
  ! index is always 0, from there we add up the number of vertices
  ! per element. Here one triangle, 2 quadrangles and 2 pentagons.
  ! The connectivity between elements and vertices is an array of
  ! size connec_idx(n_elts+1) (here 21).
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

  ! Finalize mesh :
  ! CWIPI hides the parallelism for users, so it is not
  ! mandatory to give a global numbering for mesh data (the
  ! NULL arguments earlier). If not given this numbering is
  ! generated by CWIPI by the following function,
  ! as well as the underlying mesh data structure
  call CWP_Mesh_interf_finalize(code_names(1), &
                                coupling_name)

  ! Set the field values :
  ! Note that the user has to allocate the array for the
  ! field that will be received by code2 (CWP_FIELD_MAP_TARGET).

  ! for code1
  if (I_am_code1) then

    allocate(send_field_data(n_vtx * n_components))
    do i=1,n_vtx
      send_field_data(i) = coords(1,i)
    end do

    call CWP_Field_data_set(code_names(1),        &
                            coupling_name,        &
                            field_name,           &
                            0,                    &
                            CWP_FIELD_MAP_SOURCE, &
                            send_field_data)
  ! for code2
  else

    allocate(recv_field_data(n_vtx * n_components))

    call CWP_Field_data_set(code_names(1),        &
                            coupling_name,        &
                            field_name,           &
                            0,                    &
                            CWP_FIELD_MAP_TARGET, &
                            recv_field_data)
  endif

  ! Compute interpolation weights :
  ! Set a geometric tolerance of 10% of an element size for
  ! point localisation.
  call CWP_Spatial_interp_property_set(code_names(1), &
                                       coupling_name, &
                                       "tolerance",   &
                                       CWP_DOUBLE,    &
                                       "0.1")

  call CWP_Spatial_interp_weights_compute(code_names(1), &
                                          coupling_name)

  ! Exchange field values between codes :
  ! The field exchange functions mimic the way the associated
  ! MPI functions work, see MPI documentation for more information.

  ! for code1
  if (I_am_code1) then
    call CWP_Field_issend(code_names(1), &
                          coupling_name, &
                          field_name)
  ! for code2
  else
    call CWP_Field_irecv(code_names(1), &
                         coupling_name, &
                         field_name)
  endif

  ! for code1
  if (I_am_code1) then
    call CWP_Field_wait_issend(code_names(1), &
                               coupling_name, &
                               field_name)
  ! for code2
  else
    call CWP_Field_wait_irecv(code_names(1), &
                              coupling_name, &
                              field_name)
  endif

  ! Check interpolation :
  ! These functions allow to know how many and for which target
  ! vertices the interpolation operation has been unsuccessful.
  if (i_rank == 1) then
    n_uncomputed_tgts = CWP_N_uncomputed_tgts_get(code_names(1), &
                                                  coupling_name, &
                                                  field_name,    &
                                                  0)

    allocate(send_field_data(n_uncomputed_tgts))
    uncomputed_tgts => CWP_Uncomputed_tgts_get(code_names(1), &
                                              coupling_name, &
                                              field_name,    &
                                              0)
  endif

  ! End time step :
  call CWP_Time_step_end(code_names(1))

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
  deallocate(coords);
  deallocate(connec);
  deallocate(connec_idx);
  if (I_am_code1) then
    deallocate(send_field_data);
  else
    deallocate(recv_field_data);
  endif

  ! Finalize CWIPI :
  call CWP_Finalize()

  ! Finalize MPI :
  call MPI_Finalize(ierr)

end program fortran_new_api_polygon_sol
