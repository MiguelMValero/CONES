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

module pdm_part_to_part

  use pdm
  use pdm_pointer_array

  implicit none

  integer, parameter :: PDM_PART_TO_PART_DATA_DEF_ORDER_PART1           = 0
  integer, parameter :: PDM_PART_TO_PART_DATA_DEF_ORDER_PART1_TO_PART2  = 1
  integer, parameter :: PDM_PART_TO_PART_DATA_DEF_ORDER_PART2           = 2
  integer, parameter :: PDM_PART_TO_PART_DATA_DEF_ORDER_GNUM1_COME_FROM = 3


interface

!>
!! \brief Return number of MPI ranks
!!
!! \param [in]   ptp          Part to part structure
!!
!! \return Number of MPI ranks
!!

function PDM_part_to_part_n_ranks_get (ptp) &
result(n_ranks)                              &
bind (c, name = 'PDM_part_to_part_n_ranks_get')
  use iso_c_binding
  implicit none

  type(c_ptr), value :: ptp
  integer(c_int)     :: n_ranks

end function PDM_part_to_part_n_ranks_get


!>
!!
!! \brief Free a part to part structure
!!
!! \param [inout] ptp  Part to part structure
!!
!! \return       NULL
!!

subroutine PDM_part_to_part_free (ptp) &
bind (c, name='PDM_part_to_part_free')
  use iso_c_binding
  implicit none

  type(c_ptr), value :: ptp

end subroutine PDM_part_to_part_free


!>
!!
!! \brief Get number of partitions
!!
!! \param [in]  ptp       Pointer to \ref PDM_part_to_part_t object
!! \param [out] n_part1   Number of partitions on side 1
!! \param [out] n_part2   Number of partitions on side 2
!!
!!

subroutine PDM_part_to_part_n_part_get (ptp,     &
                                        n_part1, &
                                        n_part2) &
bind (c, name='PDM_part_to_part_n_part_get')
  use iso_c_binding
  implicit none

  type(c_ptr), value :: ptp
  integer(c_int)     :: n_part1
  integer(c_int)     :: n_part2

end subroutine PDM_part_to_part_n_part_get


subroutine PDM_part_to_part_iexch_wait (ptp,     &
                                        request) &
bind (c, name='PDM_part_to_part_iexch_wait')
  ! Finalize a non-blocking exchange (Part1→Part2)
  use iso_c_binding
  implicit none

  type(c_ptr),    value :: ptp     ! Part-to-Part instance
  integer(c_int), value :: request ! Request

end subroutine PDM_part_to_part_iexch_wait


subroutine PDM_part_to_part_reverse_iexch_wait (ptp,     &
                                                request) &
bind (c, name='PDM_part_to_part_reverse_iexch_wait')
  ! Finalize a non-blocking exchange (Part2→Part1)
  use iso_c_binding
  implicit none

  type(c_ptr),    value :: ptp     ! Part-to-Part instance
  integer(c_int), value :: request ! Request

end subroutine PDM_part_to_part_reverse_iexch_wait


!>
!!
!! \brief Wait a partial asynchronus send
!!
!! \param [in]  ptp      Part to part structure
!! \param [in]  request  Request
!!
!!

subroutine PDM_part_to_part_issend_wait (ptp,     &
                                         request) &
bind (c, name='PDM_part_to_part_issend_wait')
  use iso_c_binding
  implicit none

  type(c_ptr),    value :: ptp
  integer(c_int), value :: request

end subroutine PDM_part_to_part_issend_wait


!>
!!
!! \brief Wait a partial asynchronus recv
!!
!! \param [in]  ptp      Part to part structure
!! \param [in]  request  Request
!!
!!

subroutine PDM_part_to_part_irecv_wait_raw (ptp,     &
                                            request) &
bind (c, name='PDM_part_to_part_irecv_wait_raw')
  use iso_c_binding
  implicit none

  type(c_ptr),    value :: ptp
  integer(c_int), value :: request

end subroutine PDM_part_to_part_irecv_wait_raw


end interface



contains


subroutine PDM_part_to_part_create (ptp,                &
                                    gnum_elt1,          &
                                    n_elt1,             &
                                    n_part1,            &
                                    gnum_elt2,          &
                                    n_elt2,             &
                                    n_part2,            &
                                    part1_to_part2_idx, &
                                    part1_to_part2,     &
                                    comm)
  ! Create a Partition-to-Partition redistribution from global ids
  use iso_c_binding
  implicit none

  type(c_ptr)                       :: ptp                ! Part-to-part instance
  type(PDM_pointer_array_t), pointer:: gnum_elt1          ! Element global ids in Part1 (size : ``n_part1``)
  integer(pdm_l_num_s), pointer     :: n_elt1(:)          ! Local number of elements in Part1 (size : ``n_part1``)
  integer, intent(in)               :: n_part1            ! Number of partitions in Part1
  type(PDM_pointer_array_t), pointer:: gnum_elt2          ! Element global ids in Part2 (size : ``n_part2``)
  integer(pdm_l_num_s), pointer     :: n_elt2(:)          ! Local number of elements in Part2 (size : ``n_part2``)
  integer, intent(in)               :: n_part2            ! Number of partitions in Part2
  type(PDM_pointer_array_t), pointer:: part1_to_part2_idx ! Index for Part1→Part2 mapping (for each part, size : ``n_elt1`` + 1)
  type(PDM_pointer_array_t), pointer:: part1_to_part2     ! Part1→Part2 mapping (global ids) (for each part, size : ``part1_to_part2_idx(n_elt1 + 1)``)
  integer, intent(in)               :: comm               ! MPI communicator

  integer(c_int)                    :: c_comm

  interface
    function pdm_part_to_part_create_c (gnum_elt1,          &
                                        n_elt1,             &
                                        n_part1,            &
                                        gnum_elt2,          &
                                        n_elt2,             &
                                        n_part2,            &
                                        part1_to_part2_idx, &
                                        part1_to_part2,     &
                                        comm)               &
    result (ptp)                                            &
    bind(c, name='PDM_part_to_part_create')
      use iso_c_binding
      implicit none

      type(c_ptr)           :: ptp
      type(c_ptr),    value :: gnum_elt1
      type(c_ptr),    value :: n_elt1
      integer(c_int), value :: n_part1
      type(c_ptr),    value :: gnum_elt2
      type(c_ptr),    value :: n_elt2
      integer(c_int), value :: n_part2
      type(c_ptr),    value :: part1_to_part2_idx
      type(c_ptr),    value :: part1_to_part2
      integer(c_int), value :: comm

    end function pdm_part_to_part_create_c
  end interface

  c_comm = PDM_MPI_Comm_f2c(comm)

  ptp = pdm_part_to_part_create_c (c_loc(gnum_elt1%cptr),          &
                                   c_loc(n_elt1),                  &
                                   n_part1,                        &
                                   c_loc(gnum_elt2%cptr),          &
                                   c_loc(n_elt2),                  &
                                   n_part2,                        &
                                   c_loc(part1_to_part2_idx%cptr), &
                                   c_loc(part1_to_part2%cptr),     &
                                   c_comm)

end subroutine PDM_part_to_part_create


!<
!!
!! \brief Initialize a asynchronus issend
!!
!! \param [in]   ptp                 Part to part structure
!! \param [in]   s_data              Data size
!! \param [in]   cst_stride          Constant stride
!! \param [in]   part1_to_part2_data Data (order given by part1_to_part2 array)
!! \param [in]   tag                 Tag of the exchange
!! \param [out]  request             Request
!!
!!

subroutine PDM_part_to_part_issend_raw (ptp,        &
                                        s_data,     &
                                        cst_stride, &
                                        raw_buffer, &
                                        tag,        &
                                        request)
  use iso_c_binding
  implicit none

  type(c_ptr), value                :: ptp
  integer, intent(in)               :: cst_stride
  integer, intent(in)               :: s_data
  type(c_ptr), value                :: raw_buffer
  integer, intent(in)               :: tag
  integer, intent(out)              :: request

  interface
    subroutine PDM_part_to_part_issend_c (ptp,        &
                                          s_data,     &
                                          cst_stride, &
                                          raw_buffer, &
                                          tag,        &
                                          request)    &
    bind(c, name='PDM_part_to_part_issend_raw')
      use iso_c_binding
      implicit none

      type(c_ptr),    value :: ptp
      integer(c_int), value :: s_data
      integer(c_int), value :: cst_stride
      type(c_ptr),    value :: raw_buffer
      integer(c_int), value :: tag
      integer(c_int)        :: request

    end subroutine PDM_part_to_part_issend_c
  end interface

  call PDM_part_to_part_issend_c (ptp,        &
                                  s_data,     &
                                  cst_stride, &
                                  raw_buffer, &
                                  tag,        &
                                  request)

end subroutine PDM_part_to_part_issend_raw


!<
!!
!! \brief Initialize a asynchronus irecv
!!
!! \param [in]  ptp           Part to part structure
!! \param [in]  s_data        Data size
!! \param [in]  cst_stride    Constant stride
!! \param [in]  part2_data    Partition 2 data (order given by gnum1_come_from and ref_lnum2 arrays)
!! \param [in]  tag           Tag of the exchange
!! \param [out] request       Request
!!
!!

subroutine PDM_part_to_part_irecv_raw (ptp,        &
                                       cst_stride, &
                                       s_data,     &
                                       raw_buffer, &
                                       tag,        &
                                       request)
  use iso_c_binding
  implicit none

  type(c_ptr), value                :: ptp
  integer, intent(in)               :: cst_stride
  integer, intent(in)               :: s_data
  type(c_ptr), value                :: raw_buffer
  integer, intent(in)               :: tag
  integer, intent(out)              :: request

  interface
    subroutine PDM_part_to_part_irecv_c (ptp,        &
                                         s_data,     &
                                         cst_stride, &
                                         raw_buffer, &
                                         tag,        &
                                         request)    &
    bind(c, name='PDM_part_to_part_irecv_raw')
      use iso_c_binding
      implicit none

      type(c_ptr),    value  :: ptp
      integer(c_int), value  :: s_data
      integer(c_int), value  :: cst_stride
      type(c_ptr),    value  :: raw_buffer
      integer(c_int), value  :: tag
      integer(c_int)         :: request

    end subroutine PDM_part_to_part_irecv_c
  end interface

  call PDM_part_to_part_irecv_c (ptp,        &
                                 s_data,     &
                                 cst_stride, &
                                 raw_buffer, &
                                 tag,        &
                                 request)

end subroutine PDM_part_to_part_irecv_raw


subroutine PDM_part_to_part_iexch (ptp,              &
                                   k_comm,           &
                                   t_stride,         &
                                   t_part1_data_def, &
                                   cst_stride,       &
                                   part1_stride,     &
                                   part1_data,       &
                                   part2_stride,     &
                                   part2_data,       &
                                   request)
  ! Initiate a non-blocking exchange (Part1→Part2)
  use iso_c_binding
  implicit none

  type(c_ptr), value                 :: ptp              ! Part-to-part instance
  integer, intent(in)                :: k_comm           ! Kind of MPI communication
  integer, intent(in)                :: t_stride         ! Kind of stride
  integer, intent(in)                :: t_part1_data_def ! Kind of Part1 data definition
  integer, intent(in)                :: cst_stride       ! Constant stride
  type(PDM_pointer_array_t), pointer :: part1_stride     ! Stride of Part1 data
  type(PDM_pointer_array_t), pointer :: part1_data       ! Part1 data
  type(PDM_pointer_array_t), pointer :: part2_stride     ! Stride of Part2 data
  type(PDM_pointer_array_t), pointer :: part2_data       ! Part2 data
  integer, intent(out)               :: request          ! Request

  integer                           :: n_part1
  integer                           :: n_part2
  type(c_ptr)                       :: c_part1_stride 
  type(c_ptr)                       :: c_part2_stride 
  type(c_ptr)                       :: c_part2_data   
  integer                           :: n_ref
  integer(pdm_l_num_s), pointer     :: ref(:)                 
  integer(pdm_l_num_s), pointer     :: gnum1_come_from_idx(:) 
  integer(pdm_g_num_s), pointer     :: gnum1_come_from(:)     
  integer(pdm_l_num_s), pointer     :: stride(:)              
  integer                           :: s_part2_data
  integer                           :: i, j

  integer, allocatable              :: length_data(:)
  integer, allocatable              :: length_stride(:)
  integer(c_int)                    :: c_s_data

  interface
    subroutine PDM_part_to_part_iexch_c (ptp,              &
                                         k_comm,           &
                                         t_stride,         &
                                         t_part1_data_def, &
                                         cst_stride,       &
                                         s_data,           &
                                         part1_stride,     &
                                         part1_data,       &
                                         part2_stride,     &
                                         part2_data,       &
                                         request)          &
    bind(c, name='PDM_part_to_part_iexch')
      use iso_c_binding
      implicit none

      type(c_ptr),    value :: ptp
      integer(c_int), value :: k_comm
      integer(c_int), value :: t_stride
      integer(c_int), value :: t_part1_data_def
      integer(c_int), value :: cst_stride
      integer(c_int), value :: s_data
      type(c_ptr),    value :: part1_stride
      type(c_ptr),    value :: part1_data
      type(c_ptr)           :: part2_stride
      type(c_ptr)           :: part2_data
      integer(c_int)        :: request

    end subroutine PDM_part_to_part_iexch_c
  end interface

  c_part1_stride = C_NULL_PTR
  c_part2_stride = C_NULL_PTR
  c_part2_data   = C_NULL_PTR
  ref                 => null()
  gnum1_come_from_idx => null()
  gnum1_come_from     => null()
  stride              => null()

  c_s_data = part1_data%s_data

  if (associated (part1_stride)) then
    c_part1_stride = c_loc(part1_stride%cptr)
  endif

  call PDM_part_to_part_iexch_c (ptp,                      &
                                 k_comm,                   &
                                 t_stride,                 &
                                 t_part1_data_def,         &
                                 cst_stride,               &
                                 c_s_data,                 &
                                 c_part1_stride,           &
                                 c_loc(part1_data%cptr),   &
                                 c_part2_stride,           &
                                 c_part2_data,             &
                                 request)

  call PDM_part_to_part_n_part_get (ptp,     &
                                    n_part1, &
                                    n_part2)

  ! Compute lengths

  allocate(length_stride(n_part2))
  allocate(length_data(n_part2))

  if (t_stride .eq. PDM_STRIDE_VAR_INTERLACED) then

    do i = 1, n_part2
      call PDM_part_to_part_ref_lnum2_get (ptp,   &
                                           i-1,   &
                                           n_ref, &
                                           ref)

      call PDM_part_to_part_gnum1_come_from_get (ptp,                 &
                                                 i-1,                 &
                                                 gnum1_come_from_idx, &
                                                 gnum1_come_from)

      length_stride(i) = gnum1_come_from_idx(n_ref+1)

    end do

    call PDM_pointer_array_create (part2_stride, &
                                   n_part2,    &
                                   PDM_TYPE_INT, &
                                   c_part2_stride, &
                                   length_stride, &
                                   PDM_OWNERSHIP_KEEP)


    do i = 1, n_part2
      call PDM_part_to_part_ref_lnum2_get (ptp,   &
                                           i-1,   &
                                           n_ref, &
                                           ref)

      call PDM_part_to_part_gnum1_come_from_get (ptp,                 &
                                                 i-1,                 &
                                                 gnum1_come_from_idx, &
                                                 gnum1_come_from)

      call PDM_pointer_array_part_get (part2_stride, &
                                       i-1,          &
                                       stride)

      s_part2_data = 0
      do j = 1,gnum1_come_from_idx(n_ref+1)
        s_part2_data = s_part2_data + stride(j)
      end do

      length_data(i) = s_part2_data

    end do

  else

    do i = 1, n_part2
      call PDM_part_to_part_ref_lnum2_get (ptp,   &
                                           i-1,   &
                                           n_ref, &
                                           ref)

      call PDM_part_to_part_gnum1_come_from_get (ptp,                 &
                                                 i-1,                 &
                                                 gnum1_come_from_idx, &
                                                 gnum1_come_from)

      s_part2_data = cst_stride * gnum1_come_from_idx(n_ref+1)

      length_data(i) = s_part2_data

    enddo

  endif

  if (part1_data%type .eq. PDM_TYPE_CPTR) then
    call PDM_pointer_array_create (part2_data,      &
                                   n_part2,         &
                                   part1_data%type, &
                                   c_part2_data,    &
                                   length_data,     &
                                   PDM_OWNERSHIP_KEEP)
  else
    call PDM_pointer_array_create (part2_data,         &
                                   n_part2,            &
                                   part1_data%type,    &
                                   c_part2_data,       &
                                   length_data,        &
                                   PDM_OWNERSHIP_KEEP, &
                                   part1_data%s_data)
  endif

  deallocate(length_stride)
  deallocate(length_data)

end subroutine PDM_part_to_part_iexch


subroutine PDM_part_to_part_reverse_iexch (ptp,              &
                                           k_comm,           &
                                           t_stride,         &
                                           t_part2_data_def, &
                                           cst_stride,       &
                                           part2_stride,     &
                                           part2_data,       &
                                           part1_stride,     &
                                           part1_data,       &
                                           request)
  ! Initiate a non-blocking exchange (Part2→Part1)
  use iso_c_binding
  implicit none

  type(c_ptr), value                :: ptp              ! Part-to-part instance
  integer, intent(in)               :: k_comm           ! Kind of MPI communication
  integer, intent(in)               :: t_stride         ! Kind of stride
  integer, intent(in)               :: t_part2_data_def ! Kind of Part2 data definition
  integer, intent(in)               :: cst_stride       ! Constant stride
  type(PDM_pointer_array_t), pointer:: part2_stride     ! Stride of Part1 data
  type(PDM_pointer_array_t), pointer:: part2_data       ! Part1 data
  type(PDM_pointer_array_t), pointer:: part1_stride     ! Stride of Part2 data
  type(PDM_pointer_array_t), pointer:: part1_data       ! Part2 data
  integer, intent(out)              :: request          ! Request

  integer                           :: n_part1
  integer                           :: n_part2
  type(c_ptr)                       :: c_part1_stride 
  type(c_ptr)                       :: c_part1_data   
  integer                           :: n_elt1
  integer(pdm_l_num_s), pointer     :: part1_to_part2_idx(:)
  integer(pdm_l_num_s), pointer     :: stride(:)            
  integer                           :: s_part1_data
  integer                           :: i, j

  integer, allocatable              :: length_data(:)
  integer, allocatable              :: length_stride(:)
  integer(c_int)                    :: c_s_data

  type(c_ptr)                       :: c_part2_stride

  interface
    subroutine PDM_part_to_part_reverse_iexch_c (ptp,              &
                                                 k_comm,           &
                                                 t_stride,         &
                                                 t_part2_data_def, &
                                                 cst_stride,       &
                                                 s_data,           &
                                                 part2_stride,     &
                                                 part2_data,       &
                                                 part1_stride,     &
                                                 part1_data,       &
                                                 request)          &
    bind (c, name='PDM_part_to_part_reverse_iexch')
      use iso_c_binding
      implicit none

      type(c_ptr),    value :: ptp
      integer(c_int), value :: k_comm
      integer(c_int), value :: t_stride
      integer(c_int), value :: t_part2_data_def
      integer(c_int), value :: cst_stride
      integer(c_int), value :: s_data
      type(c_ptr),    value :: part2_stride
      type(c_ptr),    value :: part2_data
      type(c_ptr)           :: part1_stride
      type(c_ptr)           :: part1_data
      integer(c_int)        :: request

    end subroutine PDM_part_to_part_reverse_iexch_c
  end interface

  c_part1_stride = C_NULL_PTR
  c_part1_data   = C_NULL_PTR

  part1_to_part2_idx => null()
  stride             => null()
  c_part2_stride   = C_NULL_PTR

  c_s_data = part2_data%s_data

  if (associated(part2_stride)) then
    c_part2_stride = c_loc(part2_stride%cptr)
  endif

  call PDM_part_to_part_reverse_iexch_c (ptp,                      &
                                         k_comm,                   &
                                         t_stride,                 &
                                         t_part2_data_def,         &
                                         cst_stride,               &
                                         c_s_data,                 &
                                         c_part2_stride,           &
                                         c_loc(part2_data%cptr),   &
                                         c_part1_stride,           &
                                         c_part1_data,             &
                                         request)

  call PDM_part_to_part_n_part_get (ptp,     &
                                    n_part1, &
                                    n_part2)

  allocate(length_stride(n_part2))
  allocate(length_data(n_part2))

  if (t_stride .eq. PDM_STRIDE_VAR_INTERLACED) then
    do i = 1, n_part1
      call PDM_part_to_part_part1_to_part2_idx_get (ptp,                &
                                                i-1,                &
                                                n_elt1,             &
                                                part1_to_part2_idx)


      length_stride(i) = part1_to_part2_idx(n_elt1+1)

    end do

    call PDM_pointer_array_create (part1_stride, &
                                   n_part1,    &
                                   PDM_TYPE_INT, &
                                   c_part1_stride, &
                                   length_stride, &
                                   PDM_OWNERSHIP_KEEP)
  endif

  do i = 1, n_part1
    call PDM_part_to_part_part1_to_part2_idx_get (ptp,                &
                                              i-1,                &
                                              n_elt1,             &
                                              part1_to_part2_idx)

    if (t_stride .eq. PDM_STRIDE_VAR_INTERLACED) then

      call PDM_pointer_array_part_get (part1_stride, &
                                       i-1,          &
                                       stride)

      s_part1_data = 0
      do j = 1,part1_to_part2_idx(n_elt1+1)
        s_part1_data = s_part1_data + stride(j)
      end do

    else

      s_part1_data = cst_stride * part1_to_part2_idx(n_elt1+1)

    end if

    length_data(i) = s_part1_data

  end do

  if (part2_data%type .eq. PDM_TYPE_CPTR) then
    call PDM_pointer_array_create (part1_data,      &
                                   n_part1,         &
                                   part2_data%type, &
                                   c_part1_data,    &
                                   length_data,     &
                                   PDM_OWNERSHIP_KEEP)
  else
    call PDM_pointer_array_create (part1_data,         &
                                   n_part1,            &
                                   part2_data%type,    &
                                   c_part1_data,       &
                                   length_data,        &
                                   PDM_OWNERSHIP_KEEP, &
                                   part2_data%s_data)
  endif

  deallocate(length_stride)
  deallocate(length_data)

end subroutine PDM_part_to_part_reverse_iexch


subroutine PDM_part_to_part_ref_lnum2_get (ptp,         &
                                           i_part,      &
                                           n_ref_lnum2, &
                                           ref_lnum2)
  ! Get referenced Part2 elements
  use iso_c_binding
  implicit none

  type(c_ptr), value            :: ptp          ! Part-to-Part instance
  integer, intent(in)           :: i_part       ! Partition identifier
  integer, intent(out)          :: n_ref_lnum2  ! Number of referenced Part2 elements
  integer(pdm_l_num_s), pointer :: ref_lnum2(:) ! Referenced Part2 elements (1-based local ids)

  type(c_ptr)                   :: c_ref_lnum2

  interface

    subroutine PDM_part_to_part_ref_lnum2_get_c (ptp,         &
                                                 i_part,      &
                                                 n_ref_lnum2, &
                                                 ref_lnum2)   &
    bind (c, name='PDM_part_to_part_ref_lnum2_single_part_get')
      use iso_c_binding
      implicit none

      type(c_ptr),    value :: ptp
      integer(c_int), value :: i_part
      integer(c_int)        :: n_ref_lnum2
      type(c_ptr)           :: ref_lnum2

    end subroutine PDM_part_to_part_ref_lnum2_get_c

  end interface

  c_ref_lnum2 = C_NULL_PTR

  call PDM_part_to_part_ref_lnum2_get_c (ptp,         &
                                         i_part,      &
                                         n_ref_lnum2, &
                                         c_ref_lnum2)

  call c_f_pointer(c_ref_lnum2,   &
                   ref_lnum2,     &
                   [n_ref_lnum2])

end subroutine PDM_part_to_part_ref_lnum2_get


subroutine PDM_part_to_part_unref_lnum2_get (ptp,           &
                                             i_part,        &
                                             n_unref_lnum2, &
                                             unref_lnum2)
  ! Get unreferenced Part2 elements
  use iso_c_binding
  implicit none

  type(c_ptr), value            :: ptp            ! Part-to-Part instance
  integer, intent(in)           :: i_part         ! Partition identifier
  integer, intent(out)          :: n_unref_lnum2  ! Number of unreferenced Part2 elements
  integer(pdm_l_num_s), pointer :: unref_lnum2(:) ! Uneferenced Part2 elements (1-based local ids)

  type(c_ptr)                   :: c_unref_lnum2

  interface

    subroutine PDM_part_to_part_unref_lnum2_get_c (ptp,           &
                                                   i_part,        &
                                                   n_unref_lnum2, &
                                                   unref_lnum2)   &
    bind (c, name='PDM_part_to_part_unref_lnum2_single_part_get')
      use iso_c_binding
      implicit none

      type(c_ptr),    value :: ptp
      integer(c_int), value :: i_part
      integer(c_int)        :: n_unref_lnum2
      type(c_ptr)           :: unref_lnum2

    end subroutine PDM_part_to_part_unref_lnum2_get_c

  end interface

  c_unref_lnum2 = C_NULL_PTR

  call PDM_part_to_part_unref_lnum2_get_c (ptp,           &
                                           i_part,        &
                                           n_unref_lnum2, &
                                           c_unref_lnum2)

  call c_f_pointer(c_unref_lnum2,   &
                   unref_lnum2,     &
                   [n_unref_lnum2])

end subroutine PDM_part_to_part_unref_lnum2_get


subroutine PDM_part_to_part_gnum1_come_from_get (ptp,                 &
                                                 i_part,              &
                                                 gnum1_come_from_idx, &
                                                 gnum1_come_from)
  ! Get Part2→Part1 mapping for referenced Part2 elements
  use iso_c_binding
  implicit none

  type(c_ptr), value            :: ptp                    ! Part-to-Part instance
  integer, intent(in)           :: i_part                 ! Partition identifier
  integer(pdm_l_num_s), pointer :: gnum1_come_from_idx(:) ! Index for Part2→Part1 mapping (size = *n_elt2* + 1)
  integer(pdm_g_num_s), pointer :: gnum1_come_from(:)     ! Part2→Part1 mapping (global ids) (size = ``gnum1_come_from_idx(n_ref_lnum2 + 1)``)

  type(c_ptr)                   :: c_gnum1_come_from_idx
  type(c_ptr)                   :: c_gnum1_come_from    
  integer                       :: n_ref
  integer(pdm_l_num_s), pointer :: ref(:)

  interface
    subroutine PDM_part_to_part_gnum1_come_from_get_c (ptp,                 &
                                                       i_part,              &
                                                       gnum1_come_from_idx, &
                                                       gnum1_come_from)     &
    bind (c, name='PDM_part_to_part_gnum1_come_from_single_part_get')
      use iso_c_binding
      implicit none

      type(c_ptr),    value :: ptp
      integer(c_int), value :: i_part
      type(c_ptr)           :: gnum1_come_from_idx
      type(c_ptr)           :: gnum1_come_from

    end subroutine PDM_part_to_part_gnum1_come_from_get_c
  end interface

  c_gnum1_come_from_idx = C_NULL_PTR
  c_gnum1_come_from     = C_NULL_PTR
  ref => null()

  call PDM_part_to_part_ref_lnum2_get (ptp,    &
                                       i_part, &
                                       n_ref,  &
                                       ref)

  call PDM_part_to_part_gnum1_come_from_get_c (ptp,                   &
                                               i_part,                &
                                               c_gnum1_come_from_idx, &
                                               c_gnum1_come_from)

  call c_f_pointer(c_gnum1_come_from_idx, &
                   gnum1_come_from_idx,   &
                   [n_ref + 1])

  call c_f_pointer(c_gnum1_come_from, &
                   gnum1_come_from,   &
                   [gnum1_come_from_idx(n_ref+1)])

end subroutine PDM_part_to_part_gnum1_come_from_get


!>
!!
!! \brief Get indirection from part1_to_part2 to buffer send (usefull to setup buffer outside ptp )
!!
!! \param [in]   ptp                       Block to part structure
!! \param [out]  gnum1_to_send_buffer_idx  Index of data to send to gnum2 from gnum1
!!                                           (for each part size : \ref n_elt1+1)
!! \param [out]  gnum1_to_send_buffer      For each gnum1 the position in send buffer
!!
!!

subroutine PDM_part_to_part_gnum1_to_send_buffer_get (ptp,                      &
                                                      gnum1_to_send_buffer_idx, &
                                                      gnum1_to_send_buffer)
  use iso_c_binding
  implicit none

  type(c_ptr), value                :: ptp
  type(PDM_pointer_array_t), pointer :: gnum1_to_send_buffer_idx
  type(PDM_pointer_array_t), pointer :: gnum1_to_send_buffer

  integer                       :: i
  integer                       :: n_part1
  integer                       :: n_part2
  integer                       :: n_elt1
  integer(pdm_l_num_s), pointer :: part1_to_part2_idx(:)             
  integer(pdm_l_num_s), pointer :: gnum1_to_send_buffer_idx_ipart(:) 

  type(c_ptr)                   :: c_gnum1_to_send_buffer_idx 
  type(c_ptr)                   :: c_gnum1_to_send_buffer     

  integer, allocatable          :: length_gnum1_to_send_buffer(:)
  integer, allocatable          :: length_gnum1_to_send_buffer_idx(:)

  interface
    subroutine PDM_part_to_part_gnum1_to_send_buffer_get_c (ptp,                      &
                                                            gnum1_to_send_buffer_idx, &
                                                            gnum1_to_send_buffer)     &
    bind (c, name='PDM_part_to_part_gnum1_to_send_buffer_get')
      use iso_c_binding
      implicit none

      type(c_ptr),    value :: ptp
      type(c_ptr)           :: gnum1_to_send_buffer_idx
      type(c_ptr)           :: gnum1_to_send_buffer

    end subroutine PDM_part_to_part_gnum1_to_send_buffer_get_c
  end interface

  part1_to_part2_idx             => null()
  gnum1_to_send_buffer_idx_ipart => null()
  c_gnum1_to_send_buffer_idx = C_NULL_PTR
  c_gnum1_to_send_buffer     = C_NULL_PTR

  call PDM_part_to_part_gnum1_to_send_buffer_get_c (ptp,                        &
                                                    c_gnum1_to_send_buffer_idx, &
                                                    c_gnum1_to_send_buffer)

  call PDM_part_to_part_n_part_get (ptp,     &
                                    n_part1, &
                                    n_part2)

  allocate (length_gnum1_to_send_buffer_idx(n_part1))
  allocate (length_gnum1_to_send_buffer(n_part1) )

  do i = 1, n_part1

    call PDM_part_to_part_part1_to_part2_idx_get (ptp,                &
                                              i-1,                &
                                              n_elt1,             &
                                              part1_to_part2_idx)

    length_gnum1_to_send_buffer_idx(i) = n_elt1 + 1
    length_gnum1_to_send_buffer(i) = part1_to_part2_idx(n_elt1+1)

  end do

  call PDM_pointer_array_create (gnum1_to_send_buffer_idx, &
                                 n_part1,    &
                                 PDM_TYPE_INT, &
                                 c_gnum1_to_send_buffer_idx, &
                                 length_gnum1_to_send_buffer_idx, &
                                 PDM_OWNERSHIP_USER)

  do i = 1, n_part1

    call PDM_pointer_array_part_get(gnum1_to_send_buffer_idx, &
                                    i-1,  &
                                    gnum1_to_send_buffer_idx_ipart)

    length_gnum1_to_send_buffer(i) = gnum1_to_send_buffer_idx_ipart(length_gnum1_to_send_buffer(i)+1)

  end do

  call PDM_pointer_array_create (gnum1_to_send_buffer, &
                                 n_part1,    &
                                 PDM_TYPE_INT, &
                                 c_gnum1_to_send_buffer, &
                                 length_gnum1_to_send_buffer, &
                                 PDM_OWNERSHIP_USER)

  deallocate (length_gnum1_to_send_buffer_idx)
  deallocate (length_gnum1_to_send_buffer)

end subroutine PDM_part_to_part_gnum1_to_send_buffer_get


!>
!!
!! \brief Get indirection from ref_lnum2 to buffer recv (usefull to setup buffer outside ptp )
!!
!! \param [in]   ptp                       Block to part structure
!! \param [out]  recv_buffer_to_ref_lnum2  For each gnum2 the position in recv buffer ( size = gnum1_come_from_idx[n_ref_lnum2])
!!
!!


subroutine PDM_part_to_part_recv_buffer_to_ref_lnum2_get (ptp,                      &
                                                          recv_buffer_to_ref_lnum2)
  use iso_c_binding
  implicit none

  type(c_ptr), value                 :: ptp
  type(PDM_pointer_array_t), pointer :: recv_buffer_to_ref_lnum2

  integer                       :: i
  integer                       :: n_part1
  integer                       :: n_part2
  integer                       :: n_ref
  integer(pdm_l_num_s), pointer :: ref(:)                 
  integer(pdm_l_num_s), pointer :: gnum1_come_from_idx(:) 
  integer(pdm_g_num_s), pointer :: gnum1_come_from(:)     

  type(c_ptr)                   :: c_recv_buffer_to_ref_lnum2
  integer, allocatable          :: length_recv_buffer_to_ref_lnum2(:)

  interface
    subroutine PDM_part_to_part_recv_buffer_to_ref_lnum2_get_c (ptp,                      &
                                                                recv_buffer_to_ref_lnum2) &
    bind (c, name='PDM_part_to_part_recv_buffer_to_ref_lnum2_get')
      use iso_c_binding
      implicit none

      type(c_ptr),    value :: ptp
      type(c_ptr)           :: recv_buffer_to_ref_lnum2

    end subroutine PDM_part_to_part_recv_buffer_to_ref_lnum2_get_c
  end interface

  ref                 => null()
  gnum1_come_from_idx => null()
  gnum1_come_from     => null()

  c_recv_buffer_to_ref_lnum2 = C_NULL_PTR

  call PDM_part_to_part_recv_buffer_to_ref_lnum2_get_c (ptp,                        &
                                                        c_recv_buffer_to_ref_lnum2)

  call PDM_part_to_part_n_part_get (ptp,     &
                                    n_part1, &
                                    n_part2)

  allocate(length_recv_buffer_to_ref_lnum2(n_part2))

  do i = 1, n_part2

    call PDM_part_to_part_ref_lnum2_get (ptp,   &
                                         i-1,   &
                                         n_ref, &
                                         ref)

    call PDM_part_to_part_gnum1_come_from_get (ptp,                 &
                                               i-1,                 &
                                               gnum1_come_from_idx, &
                                               gnum1_come_from)

    length_recv_buffer_to_ref_lnum2(i) = gnum1_come_from_idx(n_ref+1)

  end do

  call PDM_pointer_array_create (recv_buffer_to_ref_lnum2,         &
                                 n_part2,                          &
                                 PDM_TYPE_INT,                     &
                                 c_recv_buffer_to_ref_lnum2,       &
                                 length_recv_buffer_to_ref_lnum2,  &
                                 PDM_OWNERSHIP_USER)

  deallocate( length_recv_buffer_to_ref_lnum2 )

end subroutine PDM_part_to_part_recv_buffer_to_ref_lnum2_get


!>
!!
!!  \brief Get buffer size and stride for send
!!
!!  \param [in]   ptp                       Block to part structure
!!  \param [out]  default_n_send_buffer     Number of entities to send (size = n_rank)
!!  \param [out]  default_i_send_buffer     Index (size = n_rank + 1)
!!
!!

subroutine PDM_part_to_part_default_send_buffer_get (ptp,                   &
                                                     default_n_send_buffer, &
                                                     default_i_send_buffer)
  use iso_c_binding
  implicit none

  type(c_ptr), value            :: ptp
  integer(pdm_l_num_s), pointer :: default_n_send_buffer(:)
  integer(pdm_l_num_s), pointer :: default_i_send_buffer(:)

  integer                       :: n_rank
  type(c_ptr)                   :: c_default_n_send_buffer = C_NULL_PTR
  type(c_ptr)                   :: c_default_i_send_buffer = C_NULL_PTR

interface
  subroutine PDM_part_to_part_default_send_buffer_get_c (ptp,                   &
                                                         default_n_send_buffer, &
                                                         default_i_send_buffer) &
  bind (c, name='PDM_part_to_part_default_send_buffer_get')
    use iso_c_binding
    implicit none

    type(c_ptr),    value :: ptp
    type(c_ptr)           :: default_n_send_buffer
    type(c_ptr)           :: default_i_send_buffer

  end subroutine PDM_part_to_part_default_send_buffer_get_c
end interface

  call PDM_part_to_part_default_send_buffer_get_c (ptp,                     &
                                                   c_default_n_send_buffer, &
                                                   c_default_i_send_buffer)

  n_rank = PDM_part_to_part_n_ranks_get(ptp)

  call c_f_pointer(c_default_n_send_buffer, &
                   default_n_send_buffer,   &
                   [n_rank])

  call c_f_pointer(c_default_i_send_buffer, &
                   default_i_send_buffer,   &
                   [n_rank+1])

end subroutine PDM_part_to_part_default_send_buffer_get


!>
!!
!!  \brief Get buffer size and stride for recv
!!
!!  \param [in]   ptp                       Block to part structure
!!  \param [out]  default_n_recv_buffer     Number of entities to recv (size = n_rank)
!!  \param [out]  default_i_recv_buffer     Index (size = n_rank + 1)
!!
!!

subroutine PDM_part_to_part_default_recv_buffer_get (ptp,                   &
                                                     default_n_recv_buffer, &
                                                     default_i_recv_buffer)
  use iso_c_binding
  implicit none

  type(c_ptr), value            :: ptp
  integer(pdm_l_num_s), pointer :: default_n_recv_buffer(:)
  integer(pdm_l_num_s), pointer :: default_i_recv_buffer(:)

  integer                       :: n_rank
  type(c_ptr)                   :: c_default_n_recv_buffer
  type(c_ptr)                   :: c_default_i_recv_buffer

interface
  subroutine PDM_part_to_part_default_recv_buffer_get_c (ptp,                   &
                                                         default_n_recv_buffer, &
                                                         default_i_recv_buffer) &
  bind (c, name='PDM_part_to_part_default_recv_buffer_get')
    use iso_c_binding
    implicit none

    type(c_ptr),    value :: ptp
    type(c_ptr)           :: default_n_recv_buffer
    type(c_ptr)           :: default_i_recv_buffer

  end subroutine PDM_part_to_part_default_recv_buffer_get_c
end interface

  c_default_n_recv_buffer = C_NULL_PTR  
  c_default_i_recv_buffer = C_NULL_PTR

  call PDM_part_to_part_default_recv_buffer_get_c (ptp,                     &
                                                   c_default_n_recv_buffer, &
                                                   c_default_i_recv_buffer)

  n_rank = PDM_part_to_part_n_ranks_get(ptp)

  call c_f_pointer(c_default_n_recv_buffer, &
                   default_n_recv_buffer,   &
                   [n_rank])

  call c_f_pointer(c_default_i_recv_buffer, &
                   default_i_recv_buffer,   &
                   [n_rank+1])

end subroutine PDM_part_to_part_default_recv_buffer_get


!>
!!
!! \brief
!!
!! \param [in]   ptp                 Pointer to \ref PDM_part_to_part_t object
!! \param [in]   i_part              Id of partition
!! \param [out]  part1_to_part2_idx  Index for part1_to_part2 array
!! \param [out]  part1_to_part2      Gnum2 for each element of part1
!!
!!

subroutine PDM_part_to_part_part1_to_part2_idx_get (ptp,                &
                                                    i_part,             &
                                                    n_elt1,             &
                                                    part1_to_part2_idx)
  use iso_c_binding
  implicit none

  type(c_ptr), value            :: ptp
  integer, intent(in)           :: i_part
  integer, intent(out)          :: n_elt1
  integer(pdm_l_num_s), pointer :: part1_to_part2_idx(:)

  type(c_ptr)                   :: c_part1_to_part2_idx

interface
  subroutine PDM_part_to_part_part1_to_part2_idx_get_c (ptp,                &
                                                        i_part,             &
                                                        n_elt1,             &
                                                        part1_to_part2_idx) &
  bind (c, name='PDM_part_to_part_part1_to_part2_idx_single_part_get')
    use iso_c_binding
    implicit none

    type(c_ptr),    value :: ptp
    integer(c_int), value :: i_part
    integer(c_int)        :: n_elt1
    type(c_ptr)           :: part1_to_part2_idx

  end subroutine PDM_part_to_part_part1_to_part2_idx_get_c
end interface

  c_part1_to_part2_idx = C_NULL_PTR

  call PDM_part_to_part_part1_to_part2_idx_get_c (ptp,                  &
                                                  i_part,               &
                                                  n_elt1,               &
                                                  c_part1_to_part2_idx)

  call c_f_pointer(c_part1_to_part2_idx, &
                   part1_to_part2_idx,   &
                   [n_elt1 + 1])


end subroutine PDM_part_to_part_part1_to_part2_idx_get


end module pdm_part_to_part
