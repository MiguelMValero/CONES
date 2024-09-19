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

program new_api
#ifdef CWP_HAVE_FORTRAN_MPI_MODULE  
    use mpi
#endif
    use cwp

    implicit none

#ifndef CWP_HAVE_FORTRAN_MPI_MODULE
    include "mpif.h"
#endif  

    integer :: n_code
    integer :: ierr
    integer :: rank, comm_world_size, local_comm_rank, local_comm_size
    integer :: is_active_rank = CWP_STATUS_ON
    character(len = 5), dimension(:), pointer :: code_names
    character(len = 16) :: cpl_id1, cpl_id2, cpl_id3, cpl_id4, cpl_id5, cpl_id6
    real(8), pointer :: coord(:,:)
    integer, pointer :: intra_comms(:)
    integer :: interp_method, block_id, i

    INTEGER(8), POINTER, DIMENSION(:) :: global_num => NULL()

    call MPI_Init(ierr)
    call MPI_Comm_rank(MPI_comm_world, rank, ierr)
    call MPI_Comm_size(MPI_comm_world, comm_world_size, ierr)

    if (comm_world_size < 3) then
      print *, "n_rank must be >= 3"
      stop
    endif

    if (rank == 0 .OR. rank == 3 .OR. rank == 6 .OR. rank == 8) then
        n_code = 1
    else if (rank == 1 .OR. rank == 4 .OR. rank == 5.OR. rank == 9) then
        n_code = 2
    else if (rank == 7) then
        n_code = 3
    else if (rank == 2) then
        n_code = 4
    end if
    allocate(code_names(n_code), intra_comms(n_code))

    if (rank == 0) then
        code_names(1) = "code1";
    else if (rank == 1) then
        code_names(1) = "code1";
        code_names(2) = "code2";
    else if (rank == 2) then
        code_names(1) = "code1";
        code_names(2) = "code2";
        code_names(3) = "code3";
        code_names(4) = "code4";
    else if (rank == 3) then
        code_names(1) = "code3";
    else if (rank == 4) then
        code_names(1) = "code3";
        code_names(2) = "code4";
    else if (rank == 5) then
        code_names(1) = "code1";
        code_names(2) = "code3";
    else if (rank == 6) then
        code_names(1) = "code2";
    else if (rank == 7) then
        code_names(1) = "code1";
        code_names(2) = "code2";
        code_names(3) = "code3";
    else if (rank == 8) then
        code_names(1) = "code4";
    else if (rank == 9) then
        code_names(1) = "code2";
        code_names(2) = "code3";
    end if

    call CWP_Init(MPI_comm_world, n_code, code_names, is_active_rank, intra_comms)

    print *, rank, " CWP_Init OK"

    do i = 1, n_code
        call MPI_Comm_rank(intra_comms(i), local_comm_rank, ierr);
        call MPI_Comm_size(intra_comms(i), local_comm_size, ierr);
    end do

    cpl_id1 = "cpl1_code1_code2";
    cpl_id2 = "cpl2_code1_code3";
    cpl_id3 = "cpl3_code2_code3";
    cpl_id4 = "cpl4_code4_code3";
    cpl_id5 = "cpl5_code1_code4";
    cpl_id6 = "cpl6_code2_code4";

    interp_method = CWP_SPATIAL_INTERP_FROM_NEAREST_SOURCES_LEAST_SQUARES

    ! cpl1
    if (rank == 0 .OR. rank == 1 .OR. rank == 2 .OR. rank == 5 .OR. rank == 7) then
        call CWP_Cpl_create("code1", cpl_id1, "code2", CWP_INTERFACE_SURFACE, &
                CWP_COMM_PAR_WITH_PART, interp_method, 1, CWP_DYNAMIC_MESH_STATIC, CWP_TIME_EXCH_USER_CONTROLLED)
    end if
    if (rank == 1 .OR. rank == 2 .OR. rank == 6 .OR. rank == 7 .OR. rank == 9) then
        call CWP_Cpl_create("code2", cpl_id1, "code1", CWP_INTERFACE_SURFACE, &
                CWP_COMM_PAR_WITH_PART, interp_method, 1, CWP_DYNAMIC_MESH_STATIC, CWP_TIME_EXCH_USER_CONTROLLED)
    end if

    ! cpl2
    if (rank == 0 .OR. rank == 1 .OR. rank == 2 .OR. rank == 5 .OR. rank == 7) then
        call CWP_Cpl_create("code1", cpl_id2, "code3", CWP_INTERFACE_SURFACE, &
        CWP_COMM_PAR_WITH_PART, interp_method, 1, CWP_DYNAMIC_MESH_STATIC, CWP_TIME_EXCH_USER_CONTROLLED)
    end if
    if (rank == 2 .OR. rank == 3 .OR. rank == 4 .OR. rank == 5 .OR. rank == 7 .OR. rank == 9) then
        call CWP_Cpl_create("code3", cpl_id2, "code1", CWP_INTERFACE_SURFACE, &
                CWP_COMM_PAR_WITH_PART, interp_method, 1, CWP_DYNAMIC_MESH_STATIC, CWP_TIME_EXCH_USER_CONTROLLED);
    end if

    ! cpl3
    if (rank == 1 .OR. rank == 2 .OR. rank == 6 .OR. rank == 7 .OR. rank == 9) then
        call CWP_Cpl_create("code2", cpl_id3, "code3", CWP_INTERFACE_SURFACE, &
                CWP_COMM_PAR_WITH_PART, interp_method, 1, CWP_DYNAMIC_MESH_STATIC, CWP_TIME_EXCH_USER_CONTROLLED);
    end if
    if (rank == 2 .OR. rank == 3 .OR. rank == 4 .OR. rank == 5 .OR. rank == 7 .OR. rank == 9) then
        call CWP_Cpl_create("code3", cpl_id3, "code2", CWP_INTERFACE_SURFACE, &
                CWP_COMM_PAR_WITH_PART, interp_method, 1, CWP_DYNAMIC_MESH_STATIC, CWP_TIME_EXCH_USER_CONTROLLED);
    end if

    ! cpl4
    if (rank == 2 .OR. rank == 4 .OR. rank == 8) then
        call CWP_Cpl_create("code4", cpl_id4, "code3", CWP_INTERFACE_SURFACE, &
                CWP_COMM_PAR_WITH_PART, interp_method, 1, CWP_DYNAMIC_MESH_STATIC, CWP_TIME_EXCH_USER_CONTROLLED);
    end if
    if (rank == 2 .OR. rank == 3 .OR. rank == 4 .OR. rank == 5 .OR. rank == 7 .OR. rank == 9) then
        call CWP_Cpl_create("code3", cpl_id4, "code4", CWP_INTERFACE_SURFACE, &
                CWP_COMM_PAR_WITH_PART, interp_method, 1, CWP_DYNAMIC_MESH_STATIC, CWP_TIME_EXCH_USER_CONTROLLED);
    end if

    ! cpl5
    if (rank == 0 .OR. rank == 1 .OR. rank == 2 .OR. rank == 5 .OR. rank == 7) then
        call CWP_Cpl_create("code1", cpl_id5, "code4", CWP_INTERFACE_SURFACE, &
                CWP_COMM_PAR_WITH_PART, interp_method, 1, CWP_DYNAMIC_MESH_STATIC, CWP_TIME_EXCH_USER_CONTROLLED);
    end if
    if (rank == 2 .OR. rank == 4 .OR. rank == 8) then
        call CWP_Cpl_create("code4", cpl_id5, "code1", CWP_INTERFACE_SURFACE, &
                CWP_COMM_PAR_WITH_PART, interp_method, 1, CWP_DYNAMIC_MESH_STATIC, CWP_TIME_EXCH_USER_CONTROLLED);
    end if

    ! cpl6
    if (rank == 1 .OR. rank == 2 .OR. rank == 6 .OR. rank == 7 .OR. rank == 9) then
        call CWP_Cpl_create("code2", cpl_id6, "code4", CWP_INTERFACE_SURFACE, &
                CWP_COMM_PAR_WITH_PART, interp_method, 1, CWP_DYNAMIC_MESH_STATIC, CWP_TIME_EXCH_USER_CONTROLLED);
    end if
    if (rank == 2 .OR. rank == 4 .OR. rank == 8) then
        call CWP_Cpl_create("code4", cpl_id6, "code2", CWP_INTERFACE_SURFACE, &
                CWP_COMM_PAR_WITH_PART, interp_method, 1, CWP_DYNAMIC_MESH_STATIC, CWP_TIME_EXCH_USER_CONTROLLED);
    end if

    if (rank == 0 .OR. rank == 1 .OR. rank == 2 .OR. rank == 5 .OR. rank == 7) then
        call CWP_Visu_set("code1", cpl_id1, 1, CWP_VISU_FORMAT_ENSIGHT, "text")
        call CWP_Time_step_beg("code1", &
                               0.d0)
    end if

    allocate(coord(3,1))
    coord(:,1) = (/9., 4., 2./)
    call CWP_Mesh_interf_vtx_set("code1", cpl_id1, 0, 1, coord, global_num)
    block_id = CWP_Mesh_interf_block_add("code1", cpl_id1, CWP_BLOCK_FACE_QUAD4)

    if (rank == 0 .OR. rank == 1 .OR. rank == 2 .OR. rank == 5 .OR. rank == 7) then
        call CWP_Time_step_end("code1")
    end if

    print *, "All done for rank", rank

    deallocate(code_names, intra_comms, coord)

    call CWP_Finalize()
    call MPI_Finalize(ierr)
end program new_api
