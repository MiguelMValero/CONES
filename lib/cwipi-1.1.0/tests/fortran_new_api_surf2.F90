!-----------------------------------------------------------------------------
! This file is part of the CWIPI library.
!
! Copyright (C) 2021-2023  ONERA
!
! This library is free software; you can redistribute it and/or
! modify it under the terms of the GNU Lesser General Public
! License as published by the Free Software Foundation; either
! version 3 of the License, or (at your option) any later version.
!
! This library is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
! Lesser General Public License for more details.
!
! You should have received a copy of the GNU Lesser General Public
! License along with this library. If not, see <http://www.gnu.org/licenses/>.
!-----------------------------------------------------------------------------
#include "cwipi_configf.h"

program testf
#ifdef CWP_HAVE_FORTRAN_MPI_MODULE
  use mpi
#endif
  use cwp
  use pdm_sphere_surf_gen
  use, intrinsic :: iso_c_binding

  implicit none

#ifndef CWP_HAVE_FORTRAN_MPI_MODULE
  include "mpif.h"
#endif

  type my_mesh
    integer(c_int)                     :: n_part
    integer(c_int),            pointer :: pn_vtx(:)      => null()
    type(PDM_pointer_array_t), pointer :: pvtx_coord     => null()
    type(PDM_pointer_array_t), pointer :: pvtx_ln_to_gn  => null()
    integer(c_int),            pointer :: pn_face(:)     => null()
    type(PDM_pointer_array_t), pointer :: pface_vtx_idx  => null()
    type(PDM_pointer_array_t), pointer :: pface_vtx      => null()
    type(PDM_pointer_array_t), pointer :: pface_ln_to_gn => null()
  end type my_mesh

  ! --------------------------------------------------------------------
  integer, parameter                 :: comm = MPI_COMM_WORLD
  integer                            :: ierr
  integer                            :: i_rank
  integer                            :: n_rank

  ! CMD args
  integer(c_int)                     :: all_n_part(2) = [1, 1]
  integer(c_int)                     :: part_method = 3
  logical                            :: disjoint_comms = .false.
  logical                            :: swap_codes     = .false.
  logical                            :: verbose        = .false.
  integer(c_int)                     :: spatial_interp_algo = CWP_SPATIAL_INTERP_FROM_LOCATION_MESH_LOCATION_OCTREE

  integer(c_long)                    :: n_subdiv    = 1
  double precision                   :: x_center    = 0.d0
  double precision                   :: y_center    = 0.d0
  double precision                   :: z_center    = 0.d0
  double precision                   :: radius      = 1.d0

  type(my_mesh), allocatable         :: mesh(:)
  integer(c_int)                     :: id_block

  type(c_ptr)                        :: c_vtx_coord = C_NULL_PTR
  double precision,          pointer :: vtx_coord1(:)    => null()
  double precision,          pointer :: vtx_coord(:,:)   => null()
  integer(c_long),           pointer :: vtx_ln_to_gn(:)  => null()
  integer(c_int),            pointer :: face_vtx_idx(:)  => null()
  integer(c_int),            pointer :: face_vtx(:)      => null()
  integer(c_long),           pointer :: face_ln_to_gn(:) => null()


  logical                            :: has_code(2)
  integer(c_int)                     :: n_code
  integer,                   pointer :: code_id(:)           => null()
  character(len=5),          pointer :: code_name(:)         => null()
  character(len=5),          pointer :: coupled_code_name(:) => null()
  integer(c_int)                     :: is_active_rank = CWP_STATUS_ON
  integer(c_int),            pointer :: intra_comms(:)       => null()
  character(len=99)                  :: coupling_name


  integer(c_int)                     :: visu_status, map_type, exch_type
  double precision,          pointer :: field_data(:) => null(), ptr(:) => null()
  type(PDM_pointer_array_t), pointer :: send_data => null(), recv_data => null()
  type(PDM_pointer_array_t), pointer :: data => null()
  integer(c_int)                     :: stride
  character(len=99)                  :: field_name

  character(len=99)                  :: part_data_name

  ! Global data
  character(len=99)                  :: global_data_name
  integer(c_int)                     :: global_stride, global_n_entity
  integer(c_int),            pointer :: global_data(:,:) => null()


  character                          :: strnum
  integer                            :: i, j, k, l
  integer                            :: iiunit = 13
  character(len=99)                  :: arg
  double precision                   :: expected
  logical                            :: error
  ! --------------------------------------------------------------------

  ! Read command line arguments
  i = 1
  do while (i <= command_argument_count())
    call get_command_argument(i, arg)
    select case(arg)
      case ("-n")
        i = i+1
        call get_command_argument(i, arg)
        read(arg, *) n_subdiv
      case ("-n_part1")
        i = i+1
        call get_command_argument(i, arg)
        read(arg, *) all_n_part(1)
      case ("-n_part2")
        i = i+1
        call get_command_argument(i, arg)
        read(arg, *) all_n_part(2)
      case ("-disjoint")
        disjoint_comms = .true.
      case ("-v")
        verbose = .true.
      case ("-swap_codes")
        swap_codes = .true.
      case ("-location")
        spatial_interp_algo = CWP_SPATIAL_INTERP_FROM_LOCATION_MESH_LOCATION_OCTREE
      case ("-knn")
        spatial_interp_algo = CWP_SPATIAL_INTERP_FROM_NEAREST_SOURCES_LEAST_SQUARES
      case ("-identity")
        spatial_interp_algo = CWP_SPATIAL_INTERP_FROM_IDENTITY
      case ("-algo")
        i = i+1
        call get_command_argument(i, arg)
        read(arg, *) spatial_interp_algo
      case default
        print *, "Invalid command argument ", arg
        stop
    end select

    i = i + 1
  enddo


  ! Initialize MPI
  call MPI_Init(ierr)
  call MPI_Comm_rank(comm, i_rank, ierr)
  call MPI_Comm_size(comm, n_rank, ierr)


  if (verbose) then
    write (strnum, '(i1)') i_rank
    open(unit=iiunit, file='fortran_new_api_surf2_'//strnum//'.log')
  endif


  ! Assign ranks to codes
  has_code(:) = .false.
  if (disjoint_comms) then
    has_code(1) = i_rank < n_rank/2;
    has_code(2) = .not.has_code(1);
  else
    has_code(1) = (i_rank < (2*n_rank) / 3);
    has_code(2) = (i_rank >= n_rank / 3);
  endif

  if (swap_codes) then
    has_code(1:2) = has_code(2:1:-1)
  endif


  n_code = count(has_code)
  allocate(code_id(n_code),           &
           code_name(n_code),         &
           coupled_code_name(n_code), &
           intra_comms(n_code),       &
           mesh(n_code))

  n_code = 0
  do i = 1,2
    if (has_code(i)) then
      n_code = n_code + 1
      code_id(n_code) = i
      write (strnum, '(i1)') i
      code_name        (n_code) = "code" // strnum
      write (strnum, '(i1)') mod(i,2)+1
      coupled_code_name(n_code) = "code" // strnum

      mesh(n_code)%n_part = all_n_part(i)

      if (verbose) then
        write(iiunit, *) "running ", code_name(n_code), ", coupled with ", coupled_code_name(n_code)
      endif
    endif
  enddo


  call CWP_output_fortran_unit_set(iiunit)

  call CWP_Init(comm,            &
                n_code,          &
                code_name,       &
                is_active_rank, &
                intra_comms)

  call MPI_Barrier(MPI_comm_world, ierr)
  if (i_rank == 0) then
    print *, "CWIPI Init OK"
  endif


  coupling_name = "fortran_new_api_surf2"
  do i = 1, n_code
    call CWP_Cpl_create(code_name(i),                  &
                        coupling_name,                 &
                        coupled_code_name(i),          &
                        CWP_INTERFACE_SURFACE,         &
                        CWP_COMM_PAR_WITH_PART,        &
                        spatial_interp_algo,           &
                        mesh(i)%n_part,                &
                        CWP_DYNAMIC_MESH_STATIC,       &
                        CWP_TIME_EXCH_USER_CONTROLLED)
  enddo

  ! this must be performed in 2 separate loops if the intra comms do overlap
  do i = 1, n_code
    call CWP_Visu_set(code_name(i),            &
                      coupling_name,           &
                      1,                       &
                      CWP_VISU_FORMAT_ENSIGHT, &
                      "text")
  enddo


  do i = 1, n_code
    call CWP_Cpl_barrier(code_name(i), &
                         coupling_name)
  enddo


  ! call MPI_Barrier(MPI_comm_world, ierr)
  if (i_rank == 0) then
    print *, "Create coupling OK"
  endif



  ! Define interface mesh
  do i = 1, n_code
    call PDM_sphere_surf_icosphere_gen_part(intra_comms(i),         &
                                            n_subdiv,               &
                                            x_center,               &
                                            y_center,               &
                                            z_center,               &
                                            radius,                 &
                                            mesh(i)%n_part,         &
                                            part_method,            &
                                            mesh(i)%pn_vtx,         &
                                            mesh(i)%pvtx_coord,     &
                                            mesh(i)%pvtx_ln_to_gn,  &
                                            mesh(i)%pn_face,        &
                                            mesh(i)%pface_vtx_idx,  &
                                            mesh(i)%pface_vtx,      &
                                            mesh(i)%pface_ln_to_gn)

    ! id_block = CWP_Mesh_interf_block_add(code_name(i),        &
    !                                      coupling_name,       &
    !                                      CWP_BLOCK_FACE_POLY)

    do j = 1, mesh(i)%n_part
      call PDM_pointer_array_part_get(mesh(i)%pvtx_ln_to_gn, &
                                      j-1,                   &
                                      vtx_ln_to_gn)

      call PDM_pointer_array_part_get(mesh(i)%pvtx_coord, &
                                      j-1,                &
                                      vtx_coord1)

      ! reshape without copy
      c_vtx_coord = c_loc(vtx_coord1)
      call c_f_pointer(c_vtx_coord, vtx_coord, [3, mesh(i)%pn_vtx(j)])

      call CWP_Mesh_interf_vtx_set(code_name(i),      &
                                   coupling_name,     &
                                   j-1,               &
                                   mesh(i)%pn_vtx(j), &
                                   vtx_coord,         &
                                   vtx_ln_to_gn)


      call PDM_pointer_array_part_get(mesh(i)%pface_ln_to_gn, &
                                      j-1,                    &
                                      face_ln_to_gn)

      call PDM_pointer_array_part_get(mesh(i)%pface_vtx_idx, &
                                      j-1,                   &
                                      face_vtx_idx)

      call PDM_pointer_array_part_get(mesh(i)%pface_vtx, &
                                      j-1,               &
                                      face_vtx)

      ! call CWP_Mesh_interf_f_poly_block_set(code_name(i),       &
      !                                       coupling_name,      &
      !                                       j-1,                &
      !                                       id_block,           &
      !                                       mesh(i)%pn_face(j), &
      !                                       face_vtx_idx,       &
      !                                       face_vtx,           &
      !                                       face_ln_to_gn)
      call CWP_Mesh_interf_from_facevtx_set(code_name(i),       &
                                            coupling_name,      &
                                            j-1,                &
                                            mesh(i)%pn_face(j), &
                                            face_vtx_idx,       &
                                            face_vtx,           &
                                            face_ln_to_gn)

    enddo

    call CWP_Mesh_interf_finalize(code_name(i),  &
                                  coupling_name)
  enddo


  call MPI_Barrier(MPI_comm_world, ierr)
  if (i_rank == 0) then
    print *, "Set mesh OK"
  endif



  ! Create fields
  visu_status = CWP_STATUS_ON
  stride = 2

  do i = 1, n_code

    if (code_id(i) == 1) then
      exch_type = CWP_FIELD_EXCH_SEND
      map_type  = CWP_FIELD_MAP_SOURCE

      call PDM_pointer_array_create(send_data,       &
                                    mesh(i)%n_part,  &
                                    PDM_TYPE_DOUBLE)

      do j = 1, mesh(i)%n_part
        allocate(ptr(mesh(i)%pn_face(j) * stride))

        call PDM_pointer_array_part_get(mesh(i)%pface_ln_to_gn, &
                                        j-1,                    &
                                        face_ln_to_gn)

        do k = 1, stride
          ptr(k::stride) = k*face_ln_to_gn(:mesh(i)%pn_face(j))
        enddo

        call PDM_pointer_array_part_set(send_data, &
                                        j-1,       &
                                        ptr)
      enddo

      data => send_data
    else
      exch_type = CWP_FIELD_EXCH_RECV
      map_type  = CWP_FIELD_MAP_TARGET

      call PDM_pointer_array_create(recv_data,       &
                                    mesh(i)%n_part,  &
                                    PDM_TYPE_DOUBLE)

      do j = 1, mesh(i)%n_part
        allocate(ptr(mesh(i)%pn_face(j) * stride))

        call PDM_pointer_array_part_set(recv_data, &
                                        j-1,       &
                                        ptr)
      enddo

      data => recv_data
    endif

    field_name="exchanged_field"
    call CWP_Field_create(code_name(i),                 &
                          coupling_name,                &
                          field_name,                   &
                          CWP_DOUBLE,                   &
                          CWP_FIELD_STORAGE_INTERLACED, &
                          stride,                       &
                          CWP_DOF_LOCATION_CELL_CENTER, &
                          exch_type,                    &
                          visu_status)

    do j = 1, mesh(i)%n_part
      call PDM_pointer_array_part_get(data,       &
                                      j-1,        &
                                      field_data)

      call CWP_Field_data_set(code_name(i),  &
                              coupling_name, &
                              field_name,    &
                              j-1,           &
                              map_type,      &
                              field_data)
    enddo

  enddo

  do i = 1, n_code
    call CWP_Time_step_beg(code_name(i),  &
                           0.d0)
  enddo

  call MPI_Barrier(MPI_comm_world, ierr)
  if (i_rank == 0) then
    print *, "Create fields OK"
  endif



  do i = 1, n_code
    ! Set tolerance for spatial interpolation
    call CWP_Spatial_interp_property_set(code_name(i),  &
                                         coupling_name, &
                                         "tolerance",   &
                                         CWP_DOUBLE,    &
                                         "1e-2")

    ! Compute spatial interpolation weights
    call CWP_Spatial_interp_weights_compute(code_name(i),  &
                                            coupling_name)
  enddo


  call MPI_Barrier(MPI_comm_world, ierr)
  if (i_rank == 0) then
    print *, "Interpolation weights computation OK"
  endif


  ! Exchange interpolated field
  do i = 1, n_code
    if (code_id(i) == 1) then
      call CWP_Field_issend(code_name(i),  &
                            coupling_name, &
                            field_name)
    else
      call CWP_Field_irecv(code_name(i),  &
                           coupling_name, &
                           field_name)
    endif
  enddo


  do i = 1, n_code
    if (code_id(i) == 1) then
      call CWP_Field_wait_issend(code_name(i),  &
                                 coupling_name, &
                                 field_name)
    else
      call CWP_Field_wait_irecv(code_name(i),  &
                                coupling_name, &
                                field_name)
    endif
  enddo

  call MPI_Barrier(MPI_comm_world, ierr)
  if (i_rank == 0) then
    print *, "Exchange interpolated field OK"
  endif


  ! Check recv field
  error = .false.
  do i = 1, n_code
    if (code_id(i) == 2) then
      do j = 1, mesh(i)%n_part
        call PDM_pointer_array_part_get(mesh(i)%pface_ln_to_gn, &
                                        j-1,                    &
                                        face_ln_to_gn)

        call PDM_pointer_array_part_get(recv_data,  &
                                        j-1,        &
                                        field_data)

        do k = 1, mesh(i)%pn_face(j)
          expected = face_ln_to_gn(k)
          if (verbose) then
            write(iiunit, *) face_ln_to_gn(k), " received ", field_data(stride*(k-1)+1:stride*k)
          endif
          do l = 1, stride
            expected = l*face_ln_to_gn(k)
            if (abs(field_data(stride*(k-1)+l) - expected) > 0.d0) then
              error = .true.
              print ("(a1,i0,a2,i0,a10,f10.1,a11,f10.1,a1)"), &
              "[", i_rank, "] ", face_ln_to_gn(k), " received ", &
              field_data(stride*(k-1)+l), " (expected ", expected, ")"
            endif
          enddo
        enddo
      enddo
    endif
  enddo

  if (error) then
    stop
  endif

  call MPI_Barrier(MPI_comm_world, ierr)
  if (i_rank == 0) then
    print *, "Check recv field OK"
  endif


  ! Part data
  if (verbose) then
    write(iiunit,*) "-- Part data --"
  endif
  part_data_name = "part_data"

  do i = 1, n_code
    if (code_id(i) == 1) then
      exch_type = CWP_PARTDATA_SEND
    else
      exch_type = CWP_PARTDATA_RECV

      ! reset recv data
      do j = 1, mesh(i)%n_part
        call PDM_pointer_array_part_get(recv_data,  &
                                        j-1,        &
                                        field_data)

        field_data(:) = -1234.d0
      enddo
    endif

    call CWP_Part_data_create(code_name(i),           &
                              coupling_name,          &
                              part_data_name,         &
                              exch_type,              &
                              mesh(i)%pface_ln_to_gn, &
                              mesh(i)%pn_face,        &
                              mesh(i)%n_part)
  enddo

  call MPI_Barrier(MPI_comm_world, ierr)
  if (i_rank == 0) then
    print *, "Create part data OK"
  endif

  do i = 1, n_code
    if (code_id(i) == 1) then
      call CWP_Part_data_issend(code_name(i),   &
                                coupling_name,  &
                                part_data_name, &
                                0,              &
                                stride,         &
                                send_data)
    else
      call CWP_Part_data_irecv(code_name(i),   &
                               coupling_name,  &
                               part_data_name, &
                               0,              &
                               stride,         &
                               recv_data)
    endif
  enddo


  do i = 1, n_code
    if (code_id(i) == 1) then
      call CWP_Part_data_wait_issend(code_name(i),   &
                                     coupling_name,  &
                                     part_data_name, &
                                     0)
    else
      call CWP_Part_data_wait_irecv(code_name(i),   &
                                    coupling_name,  &
                                    part_data_name, &
                                    0)
    endif
  enddo



  call MPI_Barrier(MPI_comm_world, ierr)
  if (i_rank == 0) then
    print *, "Exchange part data OK"
  endif

  ! Check recv part data
  error = .false.
  do i = 1, n_code
    if (code_id(i) == 2) then
      do j = 1, mesh(i)%n_part
        call PDM_pointer_array_part_get(mesh(i)%pface_ln_to_gn, &
                                        j-1,                    &
                                        face_ln_to_gn)

        call PDM_pointer_array_part_get(recv_data,  &
                                        j-1,        &
                                        field_data)

        do k = 1, mesh(i)%pn_face(j)
          expected = face_ln_to_gn(k)
          if (verbose) then
            write(iiunit, *) face_ln_to_gn(k), " received ", field_data(stride*(k-1)+1:stride*k)
          endif
          do l = 1, stride
            expected = l*face_ln_to_gn(k)
            if (abs(field_data(stride*(k-1)+l) - expected) > 1.d-9) then
              error = .true.
              print ("(a1,i0,a2,i0,a10,f10.1,a11,f10.1,a1)"), &
              "[", i_rank, "] ", face_ln_to_gn(k), " received ", &
              field_data(stride*(k-1)+l), " (expected ", expected, ")"
            endif
          enddo
        enddo
      enddo
    endif
  enddo

  if (error) then
    stop
  endif

  call MPI_Barrier(MPI_comm_world, ierr)
  if (i_rank == 0) then
    print *, "Check recv part data OK"
  endif


  do i = 1, n_code
    if (code_id(i) == 1) then
      exch_type = CWP_PARTDATA_SEND
    else
      exch_type = CWP_PARTDATA_RECV
    endif

    call CWP_Part_data_del(code_name(i),   &
                           coupling_name,  &
                           part_data_name)
  enddo



  ! Global data
  global_data_name = "global_data"
  global_stride   = 3
  global_n_entity = 4
  allocate(global_data(global_stride, global_n_entity))

  do i = 1, n_code
    if (code_id(i) == 1) then
      call CWP_Global_data_irecv(code_name(i),     &
                                 coupling_name,    &
                                 global_data_name, &
                                 global_data)
    else
      do j = 1, global_n_entity
        global_data(:,j) = [(k*j, k=1,global_stride)]
      enddo

      call CWP_Global_data_issend(code_name(i),     &
                                  coupling_name,    &
                                  global_data_name, &
                                  global_data)
    endif
  enddo


  if (verbose) then
    write(iiunit, *) "-- Global Data --"
  endif
  error = .false.
  do i = 1, n_code
    if (code_id(i) == 1) then
      call CWP_Global_data_wait_irecv(code_name(i),     &
                                      coupling_name,    &
                                      global_data_name)

      ! check
      do j = 1, global_n_entity
        if (verbose) then
          write(iiunit, *) j, "received ", global_data(:,j)
        endif
        do k = 1, global_stride
          if (global_data(k,j) /= k*j) then
            error = .true.
            print ("(a1,i0,a2,i0,a10,i0,a11,i0,a1)"), &
              "[", i_rank, "] ", j, " received ", &
              global_data(k,j), " (expected ", k*j, ")"
          endif
        enddo
      enddo
    else
      call CWP_Global_data_wait_issend(code_name(i),     &
                                       coupling_name,    &
                                       global_data_name)
        if (verbose) then
          do j = 1, global_n_entity
            write(iiunit, *) j, "sends ", global_data(:,j)
          enddo
        endif
    endif

  enddo

  if (error) then
    stop
  endif

  call MPI_Barrier(MPI_comm_world, ierr)
  if (i_rank == 0) then
    print *, "Check recv global data OK"
  endif

  deallocate(global_data)


  call MPI_Barrier(MPI_comm_world, ierr)
  if (i_rank == 0) then
    print *, "End"
  endif

  if (verbose) then
    close(iiunit)
  endif

  do i = 1, n_code
    call CWP_Time_step_end(code_name(i))
  enddo

  ! Delete interface mesh
  do i = 1, n_code
    call CWP_Mesh_interf_del(code_name(i),  &
                             coupling_name)

    call free_mesh(mesh(i))
  enddo

  ! Delete coupling
  do i = 1, n_code
    call CWP_Cpl_Del(code_name(i),  &
                     coupling_name)
  enddo


  ! Free memory
  deallocate(code_id,           &
             code_name,         &
             coupled_code_name, &
             intra_comms,       &
             mesh)

  ! Finalize
  call CWP_Finalize()
  call MPI_Finalize(ierr)

contains

  subroutine my_free(pa)
    implicit none
    type(PDM_pointer_array_t) :: pa
    integer                   :: i

    if (associated(pa%cptr)) then
      do i = 1,size(pa%cptr)
        call pdm_fortran_free_c(pa%cptr(i))
      enddo
    endif

  end subroutine my_free

  subroutine free_mesh(mesh)
    implicit none
    type(my_mesh) :: mesh

    call pdm_fortran_free_c(c_loc(mesh%pn_vtx))
    call pdm_fortran_free_c(c_loc(mesh%pn_face))
    ! call my_free(mesh%pvtx_coord)
    ! call my_free(mesh%pvtx_ln_to_gn)
    ! call my_free(mesh%pface_vtx_idx)
    ! call my_free(mesh%pface_vtx)
    ! call my_free(mesh%pface_ln_to_gn)
    call PDM_pointer_array_free(mesh%pvtx_coord)
    call PDM_pointer_array_free(mesh%pvtx_ln_to_gn)
    call PDM_pointer_array_free(mesh%pface_vtx_idx)
    call PDM_pointer_array_free(mesh%pface_vtx)
    call PDM_pointer_array_free(mesh%pface_ln_to_gn)

  end subroutine free_mesh


end program testf



