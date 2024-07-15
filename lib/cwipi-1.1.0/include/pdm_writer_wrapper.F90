module pdm_writer_wrapper

  use pdm
  use pdm_writer
  use pdm_io
  use pdm_pointer_array

  implicit none

  type my_field_t
    character(99)                      :: name
    type(PDM_pointer_array_t), pointer :: pa => null()
  end type my_field_t

  contains

  !>
  !! \brief Writter wrapper for training
  !!

  subroutine writer_wrapper (comm,           &
                             folder,         &
                             file,           &
                             n_part,         &
                             n_vtx,          &
                             pcoords,        &
                             pvtx_ln_to_gn,  &
                             n_elt,          &
                             pelt_vtx_idx,   &
                             pelt_vtx,       &
                             pelt_ln_to_gn,  &
                             cell_t,         &
                             n_face,         &
                             pcell_face_idx, &
                             pcell_face,     &
                             format_opt,     &
                             elt_field,      &
                             vtx_field)

    use iso_c_binding
    implicit none

    ! input
    integer,          intent(in)                 :: comm
    character(len = *)                           :: folder
    character(len = *)                           :: file
    integer,          intent(in)                 :: n_part
    integer(pdm_l_num_s),      pointer           :: n_vtx(:)
    type(PDM_pointer_array_t), pointer           :: pcoords
    type(PDM_pointer_array_t), pointer           :: pvtx_ln_to_gn
    integer(pdm_l_num_s),      pointer           :: n_elt(:)
    type(PDM_pointer_array_t), pointer           :: pelt_vtx_idx
    type(PDM_pointer_array_t), pointer           :: pelt_vtx
    type(PDM_pointer_array_t), pointer           :: pelt_ln_to_gn
    integer,          intent(in),       optional :: cell_t
    integer(pdm_l_num_s),      pointer, optional :: n_face(:)
    type(PDM_pointer_array_t), pointer, optional :: pcell_face_idx
    type(PDM_pointer_array_t), pointer, optional :: pcell_face
    character(len = *),                 optional :: format_opt
    type(my_field_t),                   optional :: elt_field(:)
    type(my_field_t),                   optional :: vtx_field(:)

    double precision,          pointer           :: coords(:,:)      => null()
    integer(pdm_g_num_s),      pointer           :: vtx_ln_to_gn(:)  => null()
    integer(pdm_l_num_s),      pointer           :: elt_vtx_idx(:)   => null()
    integer(pdm_l_num_s),      pointer           :: elt_vtx(:)       => null()
    integer(pdm_g_num_s),      pointer           :: elt_ln_to_gn(:)  => null()
    integer(pdm_l_num_s),      pointer           :: cell_face_idx(:) => null()
    integer(pdm_l_num_s),      pointer           :: cell_face(:)     => null()

    ! internal
    type(c_ptr)                                  :: wrt = C_NULL_PTR
    integer, allocatable                         :: id_var_elt_field(:)
    integer, allocatable                         :: id_var_vtx_field(:)
    integer                                      :: i_rank, ierr
    integer                                      :: i_part, i, id_block, id_geom, id_var_part, id_var_elt_gnum
    logical                                      :: is_3d_nodal, is_2d
    integer                                      :: n_vtx_field, n_elt_field
    double precision,          pointer           :: val_part(:)      => null()
    double precision,          pointer           :: val_gnum(:)      => null()
    double precision,          pointer           :: elt_val(:)       => null()
    double precision,          pointer           :: vtx_val(:)       => null()
    character(99)                                :: format

    if (present(cell_t)) then
      is_3d_nodal = cell_t .ne. -1
    else
      is_3d_nodal = .false.
    endif
    is_2d = ((.not. present(pcell_face_idx)) .or. (.not. present(pcell_face))) .and. (.not. is_3d_nodal)

    ! MPI
    call mpi_comm_rank(comm, i_rank, ierr)

    if (present(format_opt)) then
      format(:) = format_opt(:)
    else
      format = "Ensight"
    endif

    call pdm_writer_create(wrt,                      &
                           format,                   &
                           PDM_WRITER_FMT_ASCII,     &
                           PDM_WRITER_TOPO_VARIABLE, &
                           PDM_WRITER_OFF,           &
                           folder,                   &
                           file,                     &
                           comm,                     &
                           PDM_IO_KIND_MPI_SIMPLE,   &
                           1.d0,                     &
                           "")

    call pdm_writer_geom_create(wrt,     &
                                id_geom, &
                                file,    &
                                n_part)

    call pdm_writer_var_create(wrt,                     &
                               id_var_part,             &
                               PDM_WRITER_OFF,          &
                               PDM_WRITER_VAR_SCALAIRE, &
                               PDM_WRITER_VAR_ELEMENTS, &
                               "i_part")

    call pdm_writer_var_create(wrt,                     &
                               id_var_elt_gnum,         &
                               PDM_WRITER_OFF,          &
                               PDM_WRITER_VAR_SCALAIRE, &
                               PDM_WRITER_VAR_ELEMENTS, &
                               "elt_gnum")

    n_elt_field = 0
    if (present(elt_field)) then
      n_elt_field = size(elt_field)

      allocate(id_var_elt_field(n_elt_field))
      do i = 1, n_elt_field
        call pdm_writer_var_create(wrt,                       &
                                   id_var_elt_field(i),       &
                                   PDM_WRITER_OFF,            &
                                   PDM_WRITER_VAR_SCALAIRE,   &
                                   PDM_WRITER_VAR_ELEMENTS,   &
                                   elt_field(i)%name)
      enddo
    endif

    n_vtx_field = 0
    if (present(vtx_field)) then
      n_vtx_field = size(vtx_field)

      allocate(id_var_vtx_field(n_vtx_field))
      do i = 1, n_vtx_field
        call pdm_writer_var_create(wrt,                       &
                                   id_var_vtx_field(i),       &
                                   PDM_WRITER_OFF,            &
                                   PDM_WRITER_VAR_SCALAIRE,   &
                                   PDM_WRITER_VAR_VERTICES,   &
                                   vtx_field(i)%name)
      enddo
    endif

    call pdm_writer_step_beg(wrt, 0.d0)


    do i_part = 1, n_part

      call pdm_pointer_array_part_get(pvtx_ln_to_gn, &
                                      i_part-1,      &
                                      vtx_ln_to_gn)

      call pdm_pointer_array_part_get(pcoords,                   &
                                      i_part-1,                  &
                                      PDM_STRIDE_CST_INTERLACED, &
                                      3,                         &
                                      coords)

      call pdm_writer_geom_coord_set(wrt,            &
                                     id_geom,        &
                                     i_part-1,       &
                                     n_vtx(i_part),  &
                                     coords,         &
                                     vtx_ln_to_gn,   &
                                     PDM_OWNERSHIP_USER)

      if (is_2d) then
         call pdm_pointer_array_part_get(pelt_ln_to_gn, &
                                         i_part-1,       &
                                         elt_ln_to_gn)

         call pdm_pointer_array_part_get(pelt_vtx_idx, &
                                         i_part-1,      &
                                         elt_vtx_idx)

         call pdm_pointer_array_part_get(pelt_vtx, &
                                         i_part-1,  &
                                         elt_vtx)

         call pdm_writer_geom_faces_facesom_add(wrt,           &
                                                id_geom,       &
                                                i_part-1,      &
                                                n_elt(i_part), &
                                                elt_vtx_idx,   &
                                                null(),        &
                                                elt_vtx,       &
                                                elt_ln_to_gn)
      else
        if (is_3d_nodal) then
          call pdm_pointer_array_part_get(pelt_vtx, &
                                          i_part-1, &
                                          elt_vtx)

          call pdm_pointer_array_part_get(pelt_ln_to_gn, &
                                          i_part-1,      &
                                          elt_ln_to_gn)

          call PDM_writer_geom_bloc_add(wrt,                &
                                        id_geom,            &
                                        cell_t,             &
                                        PDM_OWNERSHIP_USER, &
                                        id_block)

          call PDM_writer_geom_bloc_std_set(wrt,           &
                                            id_geom,       &
                                            id_block,      &
                                            i_part-1,      &
                                            n_elt(i_part), &
                                            elt_vtx,       &
                                            elt_ln_to_gn)
        else
          call pdm_pointer_array_part_get(pelt_vtx_idx, &
                                          i_part-1,     &
                                          elt_vtx_idx)

          call pdm_pointer_array_part_get(pelt_vtx, &
                                          i_part-1, &
                                          elt_vtx)

          call pdm_pointer_array_part_get(pcell_face_idx, &
                                          i_part-1,     &
                                          cell_face_idx)

          call pdm_pointer_array_part_get(pcell_face, &
                                          i_part-1, &
                                          cell_face)

          call pdm_pointer_array_part_get(pelt_ln_to_gn, &
                                          i_part-1,      &
                                          elt_ln_to_gn)

          call PDM_writer_geom_cell3d_cellface_add(wrt,            &
                                                   id_geom,        &
                                                   i_part-1,       &
                                                   n_elt(i_part),  &
                                                   n_face(i_part), &
                                                   elt_vtx_idx,    &
                                                   null(),         &
                                                   elt_vtx,        &
                                                   cell_face_idx,  &
                                                   null(),         &
                                                   cell_face,      &
                                                   elt_ln_to_gn)
        endif
      endif

    enddo

    call pdm_writer_geom_write(wrt, id_geom)

    do i_part = 1, n_part

      allocate(val_part(n_elt(i_part)), &
               val_gnum(n_elt(i_part)))

      call pdm_pointer_array_part_get(pelt_ln_to_gn, &
                                      i_part-1,       &
                                      elt_ln_to_gn)

      val_part(1:n_elt(i_part)) = i_rank*n_part + i_part
      val_gnum(1:n_elt(i_part)) = elt_ln_to_gn(1:n_elt(i_part))

      call pdm_writer_var_set(wrt,         &
                              id_var_part, &
                              id_geom,     &
                              i_part-1,    &
                              val_part)

      call pdm_writer_var_set(wrt,             &
                              id_var_elt_gnum, &
                              id_geom,         &
                              i_part-1,        &
                              val_gnum)

      deallocate(val_part, val_gnum)
    enddo

    call pdm_writer_var_write(wrt, id_var_part)
    call pdm_writer_var_write(wrt, id_var_elt_gnum)

    do i = 1, n_elt_field
      do i_part = 1, n_part
        call pdm_pointer_array_part_get(elt_field(i)%pa, &
                                        i_part-1,               &
                                        elt_val)
        call pdm_writer_var_set(wrt,                       &
                                id_var_elt_field(i), &
                                id_geom,                   &
                                i_part-1,                  &
                                elt_val)
      enddo
      call pdm_writer_var_write(wrt, id_var_elt_field(i))
    enddo

    do i = 1, n_vtx_field
      do i_part = 1, n_part
        call pdm_pointer_array_part_get(vtx_field(i)%pa, &
                                        i_part-1,               &
                                        vtx_val)
        call pdm_writer_var_set(wrt,                       &
                                id_var_vtx_field(i), &
                                id_geom,                   &
                                i_part-1,                  &
                                vtx_val)
      enddo
      call pdm_writer_var_write(wrt, id_var_vtx_field(i))
    enddo

    call pdm_writer_step_end(wrt)

    call pdm_writer_free(wrt)

  end subroutine writer_wrapper

end module pdm_writer_wrapper
