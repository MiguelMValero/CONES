!-----------------------------------------------------------------------------
! This file is part of the ParaDiGM library.
!
! Copyright (C) 2022  ONERA
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

#include "pdm_configf.h"

program testf

  use pdm
#ifdef PDM_HAVE_FORTRAN_MPI_MODULE
  use mpi
#endif

  use pdm_pointer_array
  use pdm_part_to_part
  use iso_c_binding
  use pdm_fortran

  implicit none

#ifndef PDM_HAVE_FORTRAN_MPI_MODULE
  include "mpif.h"
#endif

  type my_type

    integer(pdm_g_num_s), pointer :: g_nums(:)             => null()
    integer(pdm_l_num_s), pointer :: part1_to_part2_idx(:) => null()
    integer(pdm_g_num_s), pointer :: part1_to_part2(:)     => null()

    integer(pdm_l_num_s), pointer :: stride(:)             => null()
    integer(pdm_g_num_s), pointer :: data(:)               => null()

  end type my_type


  !-----------------------------------------------------------
  integer, parameter                    :: comm = MPI_COMM_WORLD

  type(c_ptr)                           :: ptp = C_NULL_PTR

  integer                               :: n_part1
  integer                               :: n_part2

  type(my_type), allocatable            :: my_part1(:)
  type(my_type), allocatable            :: my_part2(:)

  integer(pdm_l_num_s), pointer         :: n_elt1(:) => null()
  type(PDM_pointer_array_t), pointer    :: gnum_elt1  => null()
  type(PDM_pointer_array_t), pointer    :: part1_to_part2_idx  => null()
  type(PDM_pointer_array_t), pointer    :: part1_to_part2  => null()
  integer(pdm_l_num_s), pointer         :: n_elt2(:) => null()
  type(PDM_pointer_array_t), pointer    :: gnum_elt2 => null()

  type(PDM_pointer_array_t), pointer    :: part1_stride  => null()
  type(PDM_pointer_array_t), pointer    :: part1_data  => null()
  type(PDM_pointer_array_t), pointer    :: part2_stride  => null()
  type(PDM_pointer_array_t), pointer    :: part2_data  => null()

  type(PDM_pointer_array_t), pointer    :: part1_stride_r  => null()
  type(PDM_pointer_array_t), pointer    :: part1_data_r  => null()
  type(PDM_pointer_array_t), pointer    :: part2_stride_r  => null()
  type(PDM_pointer_array_t), pointer    :: part2_data_r  => null()

  integer                               :: n_ref_num2
  integer                               :: n_unref_num2
  integer(pdm_l_num_s), pointer         :: ref_num2(:)   => null()
  integer(pdm_l_num_s), pointer         :: unref_num2(:) => null()
  integer(pdm_l_num_s), pointer         :: gnum1_come_from_idx(:) => null()
  integer(pdm_g_num_s), pointer         :: gnum1_come_from(:)     => null()

  integer(pdm_l_num_s), pointer         :: stride(:) => null()
  integer(pdm_g_num_s), pointer         :: data(:)   => null()

  integer                               :: code
  integer                               :: i_rank
  integer                               :: n_rank
  integer                               :: request

  integer                               :: s_part_data

  integer                               :: i, j, k, idx
  character                             :: strnum
  integer                               :: fid = 13
  integer                               :: debug = 0

  !-----------------------------------------------------------

  call mpi_init(code)
  call mpi_comm_rank(comm, i_rank, code)
  call mpi_comm_size(comm, n_rank, code)


  if (n_rank .ne. 2) then
    write(*,*) 'Error : 2 MPI processes are mandatory'
    call mpi_finalize(code)
    stop
  end if

  if (debug .eq. 1) then
    write (strnum, '(i1)') i_rank
    open(unit=fid, file="part_to_part_"//strnum//".log", action='write')
  endif

  !  Define partitions
  n_part1 = 3
  n_part2 = 3

  allocate(n_elt1(n_part1))
  call PDM_pointer_array_create(gnum_elt1,          &
                                n_part1,            &
                                PDM_TYPE_G_NUM)
  call PDM_pointer_array_create(part1_to_part2_idx, &
                                n_part1,            &
                                PDM_TYPE_INT)
  call PDM_pointer_array_create(part1_to_part2,     &
                                n_part1,            &
                                PDM_TYPE_G_NUM)

  allocate(n_elt2(n_part2))
  call PDM_pointer_array_create(gnum_elt2,          &
                                n_part2,            &
                                PDM_TYPE_G_NUM)


  allocate(my_part1(n_part1), my_part2(n_part2))
  if (i_rank .eq. 0) then

    !! Part 1
    n_elt1 = [1, 4, 2]

    allocate(my_part1(1)%g_nums(n_elt1(1)))
    allocate(my_part1(1)%part1_to_part2_idx(n_elt1(1)+1))
    my_part1(1)%g_nums = [9]
    my_part1(1)%part1_to_part2_idx = [0, 0]
    allocate(my_part1(1)%part1_to_part2(my_part1(1)%part1_to_part2_idx(n_elt1(1)+1)))
    ! my_part1(1)%part1_to_part2 = []

    allocate(my_part1(2)%g_nums(n_elt1(2)))
    allocate(my_part1(2)%part1_to_part2_idx(n_elt1(2)+1))
    my_part1(2)%g_nums = [15, 13, 5, 2]
    my_part1(2)%part1_to_part2_idx = [0, 2, 4, 5, 6]
    allocate(my_part1(2)%part1_to_part2(my_part1(2)%part1_to_part2_idx(n_elt1(2)+1)))
    my_part1(2)%part1_to_part2 = [9, 11, 9, 10, 11, 10]

    allocate(my_part1(3)%g_nums(n_elt1(3)))
    allocate(my_part1(3)%part1_to_part2_idx(n_elt1(3)+1))
    my_part1(3)%g_nums = [3, 8]
    my_part1(3)%part1_to_part2_idx = [0, 1, 1]
    allocate(my_part1(3)%part1_to_part2(my_part1(3)%part1_to_part2_idx(n_elt1(3)+1)))
    my_part1(3)%part1_to_part2 = [9]


    !! Part 2
    n_elt2 = [1, 4, 2]

    allocate(my_part2(1)%g_nums(n_elt2(1)))
    my_part2(1)%g_nums = [9]

    allocate(my_part2(2)%g_nums(n_elt2(2)))
    my_part2(2)%g_nums = [15, 13, 5, 2]

    allocate(my_part2(3)%g_nums(n_elt2(3)))
    my_part2(3)%g_nums = [3, 8]


  else

    !! Part 1
    n_elt1 = [2, 6, 2]

    allocate(my_part1(1)%g_nums(n_elt1(1)))
    allocate(my_part1(1)%part1_to_part2_idx(n_elt1(1)+1))
    my_part1(1)%g_nums = [12, 16]
    my_part1(1)%part1_to_part2_idx = [0, 0, 0]
    allocate(my_part1(1)%part1_to_part2(my_part1(1)%part1_to_part2_idx(n_elt1(1)+1)))
    ! my_part1(1)%part1_to_part2 = []

    allocate(my_part1(2)%g_nums(n_elt1(2)))
    allocate(my_part1(2)%part1_to_part2_idx(n_elt1(2)+1))
    my_part1(2)%g_nums = [11, 10, 17, 14, 1, 4]
    my_part1(2)%part1_to_part2_idx = [0, 0, 0, 2, 4, 6, 8]
    allocate(my_part1(2)%part1_to_part2(my_part1(2)%part1_to_part2_idx(n_elt1(2)+1)))
    my_part1(2)%part1_to_part2 = [12, 11, 12, 10, 16, 11, 16, 10]

    allocate(my_part1(3)%g_nums(n_elt1(3)))
    allocate(my_part1(3)%part1_to_part2_idx(n_elt1(3)+1))
    my_part1(3)%g_nums = [7, 6]
    my_part1(3)%part1_to_part2_idx = [0, 1, 2]
    allocate(my_part1(3)%part1_to_part2(my_part1(3)%part1_to_part2_idx(n_elt1(3)+1)))
    my_part1(3)%part1_to_part2 = [12, 16]


    !! Part 2
    n_elt2 = [2, 6, 2]

    allocate(my_part2(1)%g_nums(n_elt2(1)))
    my_part2(1)%g_nums = [12, 16]

    allocate(my_part2(2)%g_nums(n_elt2(2)))
    my_part2(2)%g_nums = [11, 10, 17, 14, 1, 4]

    allocate(my_part2(3)%g_nums(n_elt2(3)))
    my_part2(3)%g_nums = [7, 6]


  end if




  do i = 1, n_part1
    call PDM_pointer_array_part_set(gnum_elt1,                      &
                                    i-1,                            &
                                    my_part1(i)%g_nums)

    call PDM_pointer_array_part_set(part1_to_part2_idx,             &
                                    i-1,                            &
                                    my_part1(i)%part1_to_part2_idx)

    call PDM_pointer_array_part_set(part1_to_part2,                 &
                                    i-1,                            &
                                    my_part1(i)%part1_to_part2)
  end do


  do i = 1, n_part2
    call PDM_pointer_array_part_set(gnum_elt2,                      &
                                    i-1,                            &
                                    my_part2(i)%g_nums)
  end do



  !  Create Part-to-part object
  call PDM_part_to_part_create (ptp,                &
                                gnum_elt1,          &
                                n_elt1,             &
                                n_part1,            &
                                gnum_elt2,          &
                                n_elt2,             &
                                n_part2,            &
                                part1_to_part2_idx, &
                                part1_to_part2,     &
                                comm)

  !  Check Part2 ref/unref/gnum1_come_from
  if (debug .eq. 1) then
    write (fid, *) ""
    write (fid, *) ""
    write (fid, *) "---- Check ref/unref/gnum1_come_from ----"
  endif
  do i = 1, n_part2
    call PDM_part_to_part_ref_lnum2_get (ptp,        &
                                         i-1,        &
                                         n_ref_num2, &
                                         ref_num2)

    call PDM_part_to_part_unref_lnum2_get (ptp,          &
                                           i-1,          &
                                           n_unref_num2, &
                                           unref_num2)

    call PDM_part_to_part_gnum1_come_from_get (ptp,                 &
                                               i-1,                 &
                                               gnum1_come_from_idx, &
                                               gnum1_come_from)

    if (debug .eq. 1) then
      write (fid, *) "part2", i
      write (fid, *) "referenced (l_num) :", ref_num2
      write (fid, *) "unreferenced (g_num) :", my_part2(i)%g_nums(unref_num2)
      write (fid, *) "gnum1_come_from_idx :", gnum1_come_from_idx(1:n_ref_num2+1)
      write (fid, *) "referenced (g_num) -> gnum1_come_from :"
      do j = 1, n_ref_num2
        write (fid, *) my_part2(i)%g_nums(ref_num2(j)), "->", gnum1_come_from(gnum1_come_from_idx(j)+1:gnum1_come_from_idx(j+1))
      end do
      write (fid, *) ""
    endif  
  end do

  !  Prepare exchange from 1 to 2
  call PDM_pointer_array_create (part1_stride,   &
                                 n_part1,        &
                                 PDM_TYPE_INT)

  call PDM_pointer_array_create (part1_data,     &
                                 n_part1,        &
                                 PDM_TYPE_G_NUM)
  do i = 1, n_part1
    allocate(my_part1(i)%stride(n_elt1(i)))
    ! my_part1(i)%stride(:) = 1
    s_part_data = 0
    do j = 1, n_elt1(i)
      my_part1(i)%stride(j) = int(my_part1(i)%g_nums(j))
      s_part_data = s_part_data + my_part1(i)%stride(j)
    end do

    call PDM_pointer_array_part_set (part1_stride,       &
                                     i-1,                &
                                     my_part1(i)%stride)

    ! allocate(my_part1(i)%data(n_elt1(i)))
    ! my_part1(i)%data(:) = my_part1(i)%g_nums(:)
    allocate(my_part1(i)%data(s_part_data))
    idx = 1
    do j = 1, n_elt1(i)
      do k = 1, int(my_part1(i)%g_nums(j))
        my_part1(i)%data(idx) = k
        idx = idx + 1
      end do
    end do

    call PDM_pointer_array_part_set (part1_data,         &
                                     i-1,                &
                                     my_part1(i)%data)
  end do


  print *, "Step 2"


  !  Start non-blocking exchange
  call PDM_part_to_part_iexch (ptp,                                   &
                               PDM_MPI_COMM_KIND_P2P,                 & ! k_comm
                               PDM_STRIDE_VAR_INTERLACED,             & ! t_stride
                               PDM_PART_TO_PART_DATA_DEF_ORDER_PART1, & ! t_part1_data_def
                               1,                                     & ! cst_stride
                               part1_stride,                          &
                               part1_data,                            &
                               part2_stride,                          &
                               part2_data,                            &
                               request)

  print *, "Step 3"


  !  Do stuff to cover MPI communcations
  ! ...

  !  Wait for the exchange to finish
  call PDM_part_to_part_iexch_wait (ptp, &
                                    request)

  print *, "Step 4"

  !  Check data received on Part2
  if (debug .eq. 1) then
    write (fid, *) ""
    write (fid, *) ""
    write (fid, *) "---- Check iexch ----"
  endif  

  do i = 1, n_part2

    call PDM_pointer_array_part_get (part2_stride, &
                                     i-1,          &
                                     stride)

    call PDM_part_to_part_ref_lnum2_get (ptp,        &
                                         i-1,        &
                                         n_ref_num2, &
                                         ref_num2)

    call PDM_part_to_part_gnum1_come_from_get (ptp,                 &
                                               i-1,                 &
                                               gnum1_come_from_idx, &
                                               gnum1_come_from)

    call PDM_pointer_array_part_get (part2_data,   &
                                     i-1,          &
                                     data)


    if (debug .eq. 1) then
      write (fid, *) "part2", i
      write (fid, *) "stride :", stride
      write (fid, *) "referenced (g_num) -> data :"
      idx = 1
      do j = 1, n_ref_num2
        write (fid, *) my_part2(i)%g_nums(ref_num2(j)), " :"
        do k = gnum1_come_from_idx(j)+1, gnum1_come_from_idx(j+1)
          write (fid, *) "    ", gnum1_come_from(k), " -> ", data(idx:idx+stride(k)-1)
          idx = idx + stride(k)
          ! write (fid, *) "    ", gnum1_come_from(k), " -> ", data(k)
        end do
      end do
      write (fid, *) ""
    endif  
  end do

  print *, "Step 5"


  !  Prepare exchange from 2 to 1
  call PDM_pointer_array_create (part2_stride_r, &
                                 n_part2,        &
                                 PDM_TYPE_INT)

  print *, "Step 6"

  call PDM_pointer_array_create (part2_data_r,   &
                                 n_part2,        &
                                 PDM_TYPE_G_NUM)

  print *, "Step 7"

  do i = 1, n_part2
    allocate(my_part2(i)%stride(n_elt2(i)))
    ! my_part2(i)%stride(:) = 1
    s_part_data = 0
    do j = 1, n_elt2(i)
      my_part2(i)%stride(j) = int(my_part2(i)%g_nums(j))
      s_part_data = s_part_data + my_part2(i)%stride(j)
    end do

    call PDM_pointer_array_part_set (part2_stride_r,     &
                                     i-1,                &
                                     my_part2(i)%stride)

    ! allocate(my_part2(i)%data(n_elt2(i)))
    ! my_part2(i)%data(:) = my_part2(i)%g_nums(:)
    allocate(my_part2(i)%data(s_part_data))
    idx = 1
    do j = 1, n_elt2(i)
      do k = 1, int(my_part2(i)%g_nums(j))
        my_part2(i)%data(idx) = k
        idx = idx + 1
      end do
    end do

    call PDM_pointer_array_part_set (part2_data_r,       &
                                     i-1,                &
                                     my_part2(i)%data)
  end do


  print *, "Step 8"


  !  Start reverse non-blocking exchange
  call PDM_part_to_part_reverse_iexch (ptp,                                   &
                                       PDM_MPI_COMM_KIND_P2P,                 & ! k_comm
                                       PDM_STRIDE_VAR_INTERLACED,             & ! t_stride
                                       PDM_PART_TO_PART_DATA_DEF_ORDER_PART2, & ! t_part2_data_def
                                       1,                                     & ! cst_stride
                                       part2_stride_r,                        &
                                       part2_data_r,                          &
                                       part1_stride_r,                        &
                                       part1_data_r,                          &
                                       request)

  !  Do stuff to cover MPI communcations
  ! ...

  print *, "Step 9"

  !  Wait for the exchange to finish
  call PDM_part_to_part_reverse_iexch_wait (ptp, &
                                            request)

  print *, "Step 10"

  !  Check data received on Part1
  if (debug .eq. 1) then
    write (fid, *) ""
    write (fid, *) ""
    write (fid, *) "---- Check reverse iexch ----"
  endif

  do i = 1, n_part1

    call PDM_pointer_array_part_get (part1_stride_r, &
                                     i-1,            &
                                     stride)

    call PDM_pointer_array_part_get (part1_data_r,   &
                                     i-1,            &
                                     data)

    if (debug .eq. 1) then
      write (fid, *) "part1", i
      write (fid, *) "stride :", stride
      write (fid, *) "g_num -> data :"
      idx = 1
      do j = 1, n_elt1(i)
        write (fid, *) my_part1(i)%g_nums(j), " :"
        do k = my_part1(i)%part1_to_part2_idx(j)+1, my_part1(i)%part1_to_part2_idx(j+1)
          write (fid, *) "    ", my_part1(i)%part1_to_part2(k), " -> ", data(idx:idx+stride(k)-1)
          idx = idx + stride(k)
          ! write (fid, *) "    ", my_part1(i)%part1_to_part2(k), " -> ", data(k)
        end do
      end do
      write (fid, *) ""
    endif
  end do







  !  Free memory
  call PDM_part_to_part_free (ptp)

  deallocate(n_elt1, n_elt2)
  call PDM_pointer_array_free (gnum_elt1)
  call PDM_pointer_array_free (part1_to_part2_idx)
  call PDM_pointer_array_free (part1_to_part2)
  call PDM_pointer_array_free (gnum_elt2)
  call PDM_pointer_array_free (part1_stride)
  call PDM_pointer_array_free (part1_data)
  call PDM_pointer_array_free (part2_stride)
  call PDM_pointer_array_free (part2_data)

  call PDM_pointer_array_free (part1_stride_r)
  call PDM_pointer_array_free (part1_data_r)
  call PDM_pointer_array_free (part2_stride_r)
  call PDM_pointer_array_free (part2_data_r)

  do i = 1, n_part1
    call free_my_type(my_part1(i))
  end do
  deallocate(my_part1)

  do i = 1, n_part2
    call free_my_type(my_part2(i))
  end do
  deallocate(my_part2)

  if (debug .eq. 1) then
    close(fid)
  endif  

  if (i_rank .eq. 0) then
    write(*, *) "-- End"
  end if

  call mpi_finalize(code)


contains


  subroutine free_my_type(p)

    implicit none

    type(my_type) :: p

    if (associated(p%g_nums)) then
      deallocate(p%g_nums)
    end if

    if (associated(p%part1_to_part2_idx)) then
      deallocate(p%part1_to_part2_idx)
    end if

    if (associated(p%part1_to_part2)) then
      deallocate(p%part1_to_part2)
    end if

    if (associated(p%stride)) then
      deallocate(p%stride)
    end if

    if (associated(p%data)) then
      deallocate(p%data)
    end if

  end subroutine free_my_type


end program testf
