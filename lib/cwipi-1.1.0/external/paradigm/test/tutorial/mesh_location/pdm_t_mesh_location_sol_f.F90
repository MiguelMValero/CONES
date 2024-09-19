#include "pdm_configf.h"

program tp_mesh_location

  use pdm
#ifdef PDM_HAVE_FORTRAN_MPI_MODULE
  use mpi
#endif
  use pdm_pointer_array
  use pdm_generate_mesh
  use pdm_mesh_location
  use pdm_part_to_part
  use pdm_writer_wrapper
  use iso_c_binding

  implicit none

#ifndef PDM_HAVE_FORTRAN_MPI_MODULE
  include "mpif.h"
#endif


  !---------------------------------------------------------------
  integer,              parameter    :: comm = MPI_COMM_WORLD

  ! src mesh generation
  integer(pdm_g_num_s), parameter    :: src_n_vtx_seg = 10
  integer(c_int),       parameter    :: src_n_part = 1

  integer(pdm_l_num_s),      pointer :: src_n_vtx(:)      => null()
  integer(pdm_l_num_s),      pointer :: src_n_edge(:)     => null()
  integer(pdm_l_num_s),      pointer :: src_n_face(:)     => null()
  type(PDM_pointer_array_t), pointer :: src_vtx_coord     => null()
  type(PDM_pointer_array_t), pointer :: src_edge_vtx      => null()
  type(PDM_pointer_array_t), pointer :: src_face_edge_idx => null()
  type(PDM_pointer_array_t), pointer :: src_face_edge     => null()
  type(PDM_pointer_array_t), pointer :: src_face_vtx      => null()
  type(PDM_pointer_array_t), pointer :: src_vtx_ln_to_gn  => null()
  type(PDM_pointer_array_t), pointer :: src_edge_ln_to_gn => null()
  type(PDM_pointer_array_t), pointer :: src_face_ln_to_gn => null()

  ! tgt mesh generation
  integer(pdm_g_num_s), parameter    :: tgt_n_vtx_seg = 8
  integer(c_int),       parameter    :: tgt_n_part = 1
  double precision,     parameter    :: tgt_xmin = 0.3d0
  double precision,     parameter    :: tgt_ymin = 0.3d0

  integer(pdm_l_num_s),      pointer :: tgt_n_vtx(:)      => null()
  integer(pdm_l_num_s),      pointer :: tgt_n_edge(:)     => null()
  integer(pdm_l_num_s),      pointer :: tgt_n_face(:)     => null()
  type(PDM_pointer_array_t), pointer :: tgt_vtx_coord     => null()
  type(PDM_pointer_array_t), pointer :: tgt_edge_vtx      => null()
  type(PDM_pointer_array_t), pointer :: tgt_face_edge_idx => null()
  type(PDM_pointer_array_t), pointer :: tgt_face_edge     => null()
  type(PDM_pointer_array_t), pointer :: tgt_face_vtx      => null()
  type(PDM_pointer_array_t), pointer :: tgt_vtx_ln_to_gn  => null()
  type(PDM_pointer_array_t), pointer :: tgt_edge_ln_to_gn => null()
  type(PDM_pointer_array_t), pointer :: tgt_face_ln_to_gn => null()

  ! mesh location structure
  type(c_ptr)                        :: mesh_loc = C_NULL_PTR

  ! tgt cloud definition
  integer(pdm_g_num_s),      pointer :: vtx_ln_to_gn(:)  => null()
  double precision,          pointer :: vtx_coord(:,:)   => null()

  ! src mesh definition
  logical,              parameter    :: nodal = .false.
  integer(pdm_g_num_s),      pointer :: face_ln_to_gn(:) => null()
  integer(pdm_l_num_s),      pointer :: face_edge_idx(:) => null()
  integer(pdm_l_num_s),      pointer :: face_edge(:)     => null()
  integer(pdm_l_num_s),      pointer :: face_vtx(:)      => null()
  integer(pdm_l_num_s),      pointer :: edge_vtx(:)      => null()

  type(PDM_pointer_array_t), pointer :: tgt_recv_field1  => null()
  type(PDM_pointer_array_t), pointer :: src_send_field2  => null()
  type(PDM_pointer_array_t), pointer :: tgt_recv_field2  => null()
  integer(pdm_g_num_s),      pointer :: field1(:)        => null()
  double precision,          pointer :: field2(:)        => null()

  integer                            :: n_pts, elt_n_vtx
  integer(pdm_l_num_s),      pointer :: src_to_tgt_idx(:)            => null()
  integer(pdm_g_num_s),      pointer :: points_gnum(:)               => null()
  double precision,          pointer :: points_coords(:,:)           => null()
  double precision,          pointer :: points_uvw(:,:)              => null()
  integer(pdm_l_num_s),      pointer :: points_weights_idx(:)        => null()
  double precision,          pointer :: points_weights(:)            => null()
  double precision,          pointer :: points_dist2(:)              => null()
  double precision,          pointer :: points_projected_coords(:,:) => null()
  integer(pdm_l_num_s),      pointer :: cell_vtx_idx(:)              => null()
  integer(pdm_l_num_s),      pointer :: cell_vtx(:)                  => null()

  type(c_ptr)                        :: ptp = C_NULL_PTR
  integer(c_int)                     :: request1, request2

  integer                            :: n_located
  integer(pdm_l_num_s),      pointer :: located(:)   => null()
  integer                            :: n_unlocated
  integer(pdm_l_num_s),      pointer :: unlocated(:) => null()

  double precision                   :: error
  type(my_field_t)                   :: src_visu_field1(1)
  type(my_field_t)                   :: src_visu_field2(1)
  type(my_field_t)                   :: tgt_visu_fields(3)
  double precision,          pointer :: visu_field1(:)        => null()
  double precision,          pointer :: visu_field2(:)        => null()
  double precision,          pointer :: visu_field3(:)        => null()


  integer :: i_rank, ierr
  integer :: i_part, i_elt, i_pt, i_vtx, vtx_id, i
  !---------------------------------------------------------------

  ! Initialize MPI
  call mpi_init(ierr)
  call mpi_comm_rank(comm, i_rank, ierr)

  !
  ! 1) Localization
  !

  ! Generate partitioned source mesh
  call pdm_generate_mesh_rectangle_ngon(comm,                         &
                                        PDM_MESH_NODAL_POLY_2D,       &
                                        0.d0,                         &
                                        0.d0,                         &
                                        0.d0,                         &
                                        1.d0,                         &
                                        1.d0,                         &
                                        src_n_vtx_seg,                &
                                        src_n_vtx_seg,                &
                                        src_n_part,                   &
                                        PDM_SPLIT_DUAL_WITH_PARMETIS, &
                                        src_n_vtx,                    &
                                        src_n_edge,                   &
                                        src_n_face,                   &
                                        src_vtx_coord,                &
                                        src_edge_vtx,                 &
                                        src_face_edge_idx,            &
                                        src_face_edge,                &
                                        src_face_vtx,                 &
                                        src_vtx_ln_to_gn,             &
                                        src_edge_ln_to_gn,            &
                                        src_face_ln_to_gn,            &
                                        random_factor_opt=0.8d0)

  ! Generate partitioned target mesh
  call pdm_generate_mesh_rectangle_ngon(comm,                        &
                                        PDM_MESH_NODAL_QUAD4,        &
                                        tgt_xmin,                    &
                                        tgt_ymin,                    &
                                        0.d0,                        &
                                        1.d0,                        &
                                        1.d0,                        &
                                        tgt_n_vtx_seg,               &
                                        tgt_n_vtx_seg,               &
                                        tgt_n_part,                  &
                                        PDM_SPLIT_DUAL_WITH_HILBERT, &
                                        tgt_n_vtx,                   &
                                        tgt_n_edge,                  &
                                        tgt_n_face,                  &
                                        tgt_vtx_coord,               &
                                        tgt_edge_vtx,                &
                                        tgt_face_edge_idx,           &
                                        tgt_face_edge,               &
                                        tgt_face_vtx,                &
                                        tgt_vtx_ln_to_gn,            &
                                        tgt_edge_ln_to_gn,           &
                                        tgt_face_ln_to_gn)

  ! Create the Mesh Location structure
  call pdm_mesh_location_create(mesh_loc, &
                                1,        &
                                comm,     &
                                PDM_OWNERSHIP_KEEP)


  ! Set target point cloud
  call pdm_mesh_location_n_part_cloud_set(mesh_loc, &
                                          0,        &
                                          tgt_n_part)

  do i_part = 1, tgt_n_part
    call pdm_pointer_array_part_get(tgt_vtx_ln_to_gn, &
                                    i_part-1,         &
                                    vtx_ln_to_gn)

    call pdm_pointer_array_part_get(tgt_vtx_coord,             &
                                    i_part-1,                  &
                                    PDM_STRIDE_CST_INTERLACED, &
                                    3,                         &
                                    vtx_coord)

    call pdm_mesh_location_cloud_set(mesh_loc,          &
                                     0,                 &
                                     i_part-1,          &
                                     tgt_n_vtx(i_part), &
                                     vtx_coord,         &
                                     vtx_ln_to_gn)
  enddo


  ! Set source mesh
  call pdm_mesh_location_mesh_n_part_set(mesh_loc, src_n_part)

  do i_part = 1, src_n_part

    call pdm_pointer_array_part_get(src_vtx_ln_to_gn, &
                                    i_part-1,         &
                                    vtx_ln_to_gn)

    call pdm_pointer_array_part_get(src_vtx_coord,             &
                                    i_part-1,                  &
                                    PDM_STRIDE_CST_INTERLACED, &
                                    3,                         &
                                    vtx_coord)


    call pdm_pointer_array_part_get(src_face_ln_to_gn, &
                                    i_part-1,          &
                                    face_ln_to_gn)

    call pdm_pointer_array_part_get(src_face_edge_idx, &
                                    i_part-1,          &
                                    face_edge_idx)


    if (nodal) then
      call pdm_pointer_array_part_get(src_face_vtx, &
                                      i_part-1,     &
                                      face_vtx)

      call pdm_mesh_location_nodal_part_set_2d(mesh_loc,           &
                                               i_part-1,           &
                                               src_n_face(i_part), &
                                               face_edge_idx,      &
                                               face_vtx,           &
                                               face_ln_to_gn,      &
                                               src_n_vtx(i_part),  &
                                               vtx_coord,          &
                                               vtx_ln_to_gn)
    else
      call pdm_pointer_array_part_get(src_face_edge, &
                                      i_part-1,      &
                                      face_edge)

      call pdm_pointer_array_part_get(src_edge_vtx, &
                                      i_part-1,     &
                                      edge_vtx)

      call pdm_mesh_location_part_set_2d(mesh_loc,           &
                                         i_part-1,           &
                                         src_n_face(i_part), &
                                         face_edge_idx,      &
                                         face_edge,           &
                                         face_ln_to_gn,      &
                                         src_n_edge(i_part), &
                                         edge_vtx,           &
                                         src_n_vtx(i_part),  &
                                         vtx_coord,          &
                                         vtx_ln_to_gn)

    endif
  enddo


  ! Set the location preconditioning method (optional)
  call pdm_mesh_location_method_set(mesh_loc, PDM_MESH_LOCATION_OCTREE)

  ! Set the geometric tolerance (optional)
  call pdm_mesh_location_tolerance_set(mesh_loc, 1.d-6)

  ! Compute location
  call pdm_mesh_location_compute(mesh_loc)

  ! Dump elapsed and CPU times
  call pdm_mesh_location_dump_times(mesh_loc)


  !
  ! 2) Interpolation
  !

  ! Get Part-to-part structure
  call pdm_mesh_location_part_to_part_get(mesh_loc, &
                                          0,        &
                                          ptp,      &
                                          PDM_OWNERSHIP_USER)

  ! Initiate exchange of first field
  request1 = -1
  call pdm_part_to_part_iexch(ptp,                                   &
                              PDM_MPI_COMM_KIND_P2P,                 &
                              PDM_STRIDE_CST_INTERLACED,             &
                              PDM_PART_TO_PART_DATA_DEF_ORDER_PART1, &
                              1,                                     &
                              null(),                                &
                              src_face_ln_to_gn,                     &
                              null(),                                &
                              tgt_recv_field1,                       &
                              request1)

  ! Interpolate second field (node-based)
  call pdm_pointer_array_create(src_send_field2, &
                                src_n_part,      &
                                PDM_TYPE_DOUBLE)
  do i_part = 1, src_n_part

    call pdm_mesh_location_points_in_elt_get(mesh_loc,           &
                                             0,                  &
                                             i_part-1,           &
                                             src_to_tgt_idx,     &
                                             points_gnum,        &
                                             points_coords,      &
                                             points_uvw,         &
                                             points_weights_idx, &
                                             points_weights,     &
                                             points_dist2,       &
                                             points_projected_coords)

    call pdm_pointer_array_part_get(src_vtx_coord,             &
                                    i_part-1,                  &
                                    PDM_STRIDE_CST_INTERLACED, &
                                    3,                         &
                                    vtx_coord)

    call pdm_mesh_location_cell_vertex_get(mesh_loc,     &
                                           i_part-1,     &
                                           cell_vtx_idx, &
                                           cell_vtx)

    n_pts = src_to_tgt_idx(src_n_face(i_part)+1)

    allocate(field2(n_pts))
    do i_elt = 1, src_n_face(i_part)
      do i_pt = src_to_tgt_idx(i_elt)+1, src_to_tgt_idx(i_elt+1)
        field2(i_pt) = 0.d0

        elt_n_vtx = cell_vtx_idx(i_elt+1) - cell_vtx_idx(i_elt)

        if (points_weights_idx(i_pt+1) - points_weights_idx(i_pt) /= elt_n_vtx) then
          print *, "Error elt_n_vtx"
          stop
        endif

        do i_vtx = 1, elt_n_vtx
          vtx_id = cell_vtx(cell_vtx_idx(i_elt) + i_vtx)
          field2(i_pt) = field2(i_pt) + vtx_coord(1,vtx_id) * points_weights(points_weights_idx(i_pt) + i_vtx)
        enddo

      enddo
    enddo

    call pdm_pointer_array_part_set(src_send_field2, &
                                    i_part-1,        &
                                    field2)

  enddo


  ! Initiate exchange of second field
  request2 = -1
  call pdm_part_to_part_iexch(ptp,                                            &
                              PDM_MPI_COMM_KIND_P2P,                          &
                              PDM_STRIDE_CST_INTERLACED,                      &
                              PDM_PART_TO_PART_DATA_DEF_ORDER_PART1_TO_PART2, &
                              1,                                              &
                              null(),                                         &
                              src_send_field2,                                &
                              null(),                                         &
                              tgt_recv_field2,                                &
                              request2)

  ! Finalize both exchanges
  call pdm_part_to_part_iexch_wait(ptp, request1)
  call pdm_part_to_part_iexch_wait(ptp, request2)


  ! Finally, visualize the interpolated target fields.
  ! (Beware of unlocated points!)
  do i = 1, 3
    call pdm_pointer_array_create(tgt_visu_fields(i)%pa, &
                                  tgt_n_part,            &
                                  PDM_TYPE_DOUBLE)
  enddo
  tgt_visu_fields(1)%name = "field1"
  tgt_visu_fields(2)%name = "field2"
  tgt_visu_fields(3)%name = "is_located"


  do i_part = 1, tgt_n_part

    allocate(visu_field1(tgt_n_vtx(i_part)), &
             visu_field2(tgt_n_vtx(i_part)), &
             visu_field3(tgt_n_vtx(i_part)))


    n_located = pdm_mesh_location_n_located_get(mesh_loc, &
                                                0,        &
                                                i_part-1)

    call pdm_mesh_location_located_get(mesh_loc, &
                                       0,        &
                                       i_part-1, &
                                       located)

    n_unlocated = pdm_mesh_location_n_unlocated_get(mesh_loc, &
                                                    0,        &
                                                    i_part-1)

    call pdm_mesh_location_unlocated_get(mesh_loc, &
                                         0,        &
                                         i_part-1, &
                                         unlocated)

    call pdm_pointer_array_part_get(tgt_vtx_ln_to_gn, &
                                    i_part-1,         &
                                    vtx_ln_to_gn)

    call pdm_pointer_array_part_get(tgt_vtx_coord,             &
                                    i_part-1,                  &
                                    PDM_STRIDE_CST_INTERLACED, &
                                    3,                         &
                                    vtx_coord)

    call pdm_pointer_array_part_get(tgt_recv_field1, &
                                    i_part-1,        &
                                    field1)

    call pdm_pointer_array_part_get(tgt_recv_field2, &
                                    i_part-1,        &
                                    field2)

    do i = 1, n_unlocated
      vtx_id = unlocated(i)
      visu_field1(vtx_id) = -1.d0
      visu_field2(vtx_id) = -1.d0
      visu_field3(vtx_id) =  0.d0
    enddo

    do i = 1, n_located
      vtx_id = located(i)
      error  = abs(field2(i) - vtx_coord(1,vtx_id))
      if (error > 1.e-9) then
        print *, "!! error vtx", vtx_ln_to_gn(vtx_id), " :", error
      endif

      visu_field1(vtx_id) = field1(i)
      visu_field2(vtx_id) = field2(i)
      visu_field3(vtx_id) = 1.d0

      call pdm_pointer_array_part_set(tgt_visu_fields(1)%pa, &
                                      i_part-1,              &
                                      visu_field1)

      call pdm_pointer_array_part_set(tgt_visu_fields(2)%pa, &
                                      i_part-1,              &
                                      visu_field2)

      call pdm_pointer_array_part_set(tgt_visu_fields(3)%pa, &
                                      i_part-1,              &
                                      visu_field3)

    enddo

  enddo



  call pdm_pointer_array_create(src_visu_field1(1)%pa, &
                                src_n_part,            &
                                PDM_TYPE_DOUBLE)
  call pdm_pointer_array_create(src_visu_field2(1)%pa, &
                                src_n_part,            &
                                PDM_TYPE_DOUBLE)

  src_visu_field1(1)%name = "field1"
  src_visu_field2(1)%name = "field2"

  do i_part = 1, src_n_part
  call pdm_pointer_array_part_get(src_face_ln_to_gn, &
                                  i_part-1,          &
                                  face_ln_to_gn)
    allocate(field2(src_n_face(i_part)))
    field2(1:src_n_face(i_part)) = face_ln_to_gn(1:src_n_face(i_part))
    call pdm_pointer_array_part_set(src_visu_field1(1)%pa, &
                                    i_part-1,              &
                                    field2)

    call pdm_pointer_array_part_get(src_vtx_coord,             &
                                    i_part-1,                  &
                                    PDM_STRIDE_CST_INTERLACED, &
                                    3,                         &
                                    vtx_coord)
    allocate(field2(src_n_vtx(i_part)))
    field2(1:src_n_vtx(i_part)) = vtx_coord(1,1:src_n_vtx(i_part))
    call pdm_pointer_array_part_set(src_visu_field2(1)%pa, &
                                    i_part-1,              &
                                    field2)
  enddo



  call writer_wrapper(comm,                      &
                      "mesh_location_sol_f",     &
                      "src_mesh",                &
                      src_n_part,                &
                      src_n_vtx,                 &
                      src_vtx_coord,             &
                      src_vtx_ln_to_gn,          &
                      src_n_face,                &
                      src_face_edge_idx,         &
                      src_face_vtx,              &
                      src_face_ln_to_gn,         &
                      elt_field=src_visu_field1, &
                      vtx_field=src_visu_field2)


  call writer_wrapper(comm,                  &
                      "mesh_location_sol_f", &
                      "tgt_mesh",            &
                      tgt_n_part,            &
                      tgt_n_vtx,             &
                      tgt_vtx_coord,         &
                      tgt_vtx_ln_to_gn,      &
                      tgt_n_face,            &
                      tgt_face_edge_idx,     &
                      tgt_face_vtx,          &
                      tgt_face_ln_to_gn,     &
                      vtx_field=tgt_visu_fields)

  ! Free memory
  call pdm_mesh_location_free(mesh_loc)

  call pdm_part_to_part_free(ptp)

  ! TODO: free everything...

  if (i_rank == 0) then
    print *, "The End :D"
  endif

  call mpi_finalize(ierr)


end program tp_mesh_location
