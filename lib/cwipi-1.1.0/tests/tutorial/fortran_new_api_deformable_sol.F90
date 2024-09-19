#include "cwipi_configf.h"

program fortran_new_api_deformable_sol

#ifdef CWP_HAVE_FORTRAN_MPI_MODULE
    use mpi
#endif
    use cwp
    use pdm_generate_mesh

    implicit none

#ifndef CWP_HAVE_FORTRAN_MPI_MODULE
    include "mpif.h"
#endif

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

  double precision,   dimension(:,:), pointer :: xyz_dest  => null()
  integer(c_long),    dimension(:),   pointer :: pts_g_num => null()

  integer,            dimension(:),   pointer :: elt_vtx_idx => null()
  integer,            dimension(:),   pointer :: elt_vtx     => null()
  integer(c_long),    dimension(:),   pointer :: elt_g_num   => null()
  integer(c_int)                              :: id_block

  character(len = 99)                         :: send_field_name
  character(len = 99)                         :: recv_field_name
  integer(c_int)                              :: n_components

  double precision                            :: ampl
  double precision                            :: dt
  double precision                            :: freq
  double precision                            :: omega
  double precision                            :: phi
  double precision                            :: time
  integer                                     :: i, j, it, itdeb, itend
!  integer                                     :: n2, n_partition
  integer                                     :: n_vtx, n_elt

  double precision,                   pointer :: send_field_data(:) => null()
  double precision,                   pointer :: recv_field_data(:) => null()

  integer(c_int)                              :: n_uncomputed_tgts
  integer(c_int),                     pointer :: uncomputed_tgts(:) => null()
  !--------------------------------------------------------------------

  ! MPI Initialization :
  call MPI_Init(ierr)
  call MPI_Comm_rank(mpi_comm_world, i_rank, ierr)
  call MPI_Comm_size(mpi_comm_world, n_rank, ierr)

  ! Initialize CWIPI :
  n_code = 1

  allocate(code_names(n_code), &
           intra_comms(n_code))

  code_names(1) = "code1"

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
  coupled_code_names(1) = "code2"
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

  ! Interations :
  ! At each iteration the mesh coordinates and the exchanged fields are modified.
  itdeb = 1
  itend = 10
  freq  = 0.20d0
  ampl  = 0.012d0
  phi   = 0.1d0
  time  = 0.0d0
  dt    = 0.1d0

  omega = 2.0d0*acos(-1.0d0)*freq

  allocate(send_field_data(n_elt))
  allocate(recv_field_data(n_elt))
  allocate(xyz_dest(3, n_elt))
  allocate(uncomputed_tgts(n_elt))

  send_field_name = "girafe"
  recv_field_name = "chinchilla"

  do it = itdeb, itend

    time = (it-itdeb)*dt

    ! Begin time step :
    call CWP_Time_step_beg(code_names(1), &
                           time)

    ! Deform mesh :

    do i = 1, n_vtx
      coords(3,i) = (coords(1,i)*coords(1,i) + coords(2,i)*coords(2,i)) * ampl*cos(omega*time+phi)
    enddo

    ! Update sent field :

    do  i = 1, n_elt
      send_field_data(i) = 0.
      do j = elt_vtx_idx(i)+1, elt_vtx_idx(i+1)
        send_field_data(i) = send_field_data(i) + coords(3,elt_vtx(j))
      enddo
      send_field_data(i) = send_field_data(i)/ (elt_vtx_idx(i+1) - elt_vtx_idx(i))
    enddo

    ! Update user defined degrees of freedom :

    do i = 1, n_elt
       xyz_dest(:,i) = 0.
       do j = elt_vtx_idx(i)+1, elt_vtx_idx(i+1)
           xyz_dest(:,i) = xyz_dest(:,i) + coords(:,elt_vtx(j))
       enddo
       xyz_dest(:,i) = xyz_dest(:,i) / (elt_vtx_idx(i+1) - elt_vtx_idx(i))
    enddo

    if (it == itdeb) then

      ! Set the mesh vertices coordinates :
      ! If the global numbering array is available, it can be given instead
      ! of the last vtx_g_num => null() argument. If not given, CWIPI will
      ! compute it for you in CWP_Mesh_interf_finalize.
      call CWP_Mesh_interf_vtx_set(code_names(1), &
                                   coupling_name, &
                                   0,             &
                                   n_vtx,         &
                                   coords,        &
                                   vtx_g_num)

      ! Set the mesh polygons connectivity :
      ! Since the mesh elements are triangles, CWP_BLOCK_FACE_TRIA3 could
      ! be used instead of CWP_BLOCK_FACE_POLY.
      id_block = CWP_Mesh_interf_block_add(code_names(1),       &
                                           coupling_name,       &
                                           CWP_BLOCK_FACE_POLY)

      ! CWP_Mesh_interf_from_faceedge_set is not used here since the
      ! mesh is in the form of an face->vertex connectivity. If the mesh
      ! is only available with face->edge and edge->vertex connectivity,
      ! the faceedge function should be used. The input of one does
      ! not distinguish per element type. CWIPI filters that later.
      call CWP_Mesh_interf_f_poly_block_set(code_names(1), &
                                            coupling_name, &
                                            0,             &
                                            id_block,      &
                                            n_elt,         &
                                            elt_vtx_idx,    &
                                            elt_vtx,        &
                                            elt_g_num)

      ! When CWP_DOF_LOCATION_USER, calling CWP_User_tgt_pts_set is mandatory :
      ! This is a point cloud different from the nodes of the mesh that were
      ! previously set.
      call CWP_User_tgt_pts_set(code_names(1), &
                                coupling_name, &
                                0,             &
                                n_elt,         &
                                xyz_dest,      &
                                pts_g_num)

      ! Finalize mesh :
      call CWP_Mesh_interf_finalize(code_names(1), &
                                    coupling_name)

      n_components = 1
      ! Create field :
      call CWP_Field_create(code_names(1),                &
                            coupling_name,                &
                            send_field_name,              &
                            CWP_DOUBLE,                   &
                            CWP_FIELD_STORAGE_INTERLACED, &
                            n_components,                 &
                            CWP_DOF_LOCATION_CELL_CENTER, &
                            CWP_FIELD_EXCH_SEND,          &
                            CWP_STATUS_ON)

      ! Set the field values :
      call CWP_Field_data_set(code_names(1),        &
                              coupling_name,        &
                              send_field_name,      &
                              0,                    &
                              CWP_FIELD_MAP_SOURCE, &
                              send_field_data)

      ! Create field :
      call CWP_Field_create(code_names(1),                &
                            coupling_name,                &
                            recv_field_name,              &
                            CWP_DOUBLE,                   &
                            CWP_FIELD_STORAGE_INTERLACED, &
                            n_components,                 &
                            CWP_DOF_LOCATION_USER,        &
                            CWP_FIELD_EXCH_RECV,          &
                            CWP_STATUS_ON)

      ! Set the field values :
      call CWP_Field_data_set(code_names(1),        &
                              coupling_name,        &
                              recv_field_name,      &
                              0,                    &
                              CWP_FIELD_MAP_TARGET, &
                              recv_field_data)

      ! Set interpolation property :
      call CWP_Spatial_interp_property_set(code_names(1), &
                                           coupling_name, &
                                           "tolerance",   &
                                           CWP_DOUBLE,    &
                                           "0.1")

    endif

    ! Compute interpolation weights :
    call CWP_Spatial_interp_weights_compute(code_names(1), &
                                            coupling_name)

    ! Exchange field values :
    ! If the codes operate a cross exchange, a deadlock could occur if
    ! CWP_Field_wait_issend/CWP_Field_wait_irecv calls are mixed with
    ! CWP_Field_issend/CWP_Field_irecv calls.
    call CWP_Field_irecv(code_names(1), &
                         coupling_name, &
                         recv_field_name)

    call CWP_Field_issend(code_names(1), &
                          coupling_name, &
                          send_field_name)

    call CWP_Field_wait_irecv(code_names(1), &
                              coupling_name, &
                              recv_field_name)

    call CWP_Field_wait_issend(code_names(1), &
                               coupling_name, &
                               send_field_name)

    ! Check interpolation :
    n_uncomputed_tgts = CWP_N_uncomputed_tgts_get(code_names(1),   &
                                                  coupling_name,   &
                                                  recv_field_name, &
                                                  0)

    if (n_uncomputed_tgts /= 0) then
      uncomputed_tgts => CWP_Uncomputed_tgts_get(code_names(1),   &
                                                 coupling_name,   &
                                                 recv_field_name, &
                                                 0)
    endif

    ! End time step :
    call CWP_Time_step_end(code_names(1))

  enddo

  ! Delete field :
  call CWP_Field_Del(code_names(1),   &
                     coupling_name,   &
                     send_field_name)

  call CWP_Field_Del(code_names(1),   &
                     coupling_name,   &
                     recv_field_name)

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
  deallocate(send_field_data)
  deallocate(recv_field_data)
  deallocate(xyz_dest)
  deallocate(uncomputed_tgts)

  call pdm_fortran_free_c(c_loc(coords))
  call pdm_fortran_free_c(c_loc(elt_vtx_idx))
  call pdm_fortran_free_c(c_loc(elt_vtx))

  ! Finalize CWIPI :
  call CWP_Finalize()

  ! Finalize MPI :
  call MPI_Finalize(ierr)

end program fortran_new_api_deformable_sol
