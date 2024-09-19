!-----------------------------------------------------------------------------
! This file is part of the ParaDiGM library.
!
! Copyright (C) 2019  ONERA
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

module pdm_part_to_block

  use pdm
  use pdm_pointer_array

  implicit none

  integer, parameter :: PDM_PART_TO_BLOCK_DISTRIB_ALL_PROC          = 0  ! Distribute block on all processes
  integer, parameter :: PDM_PART_TO_BLOCK_DISTRIB_ONE_PROC_PER_NODE = 1  ! Distribute block on one processe pere node
  integer, parameter :: PDM_PART_TO_BLOCK_DISTRIB_PART_OF_NODE      = 2  ! Distribute block on part of nodes

  integer, parameter :: PDM_PART_TO_BLOCK_POST_NOTHING              = 0  ! No post processing
  integer, parameter :: PDM_PART_TO_BLOCK_POST_CLEANUP              = 1  ! Cleanup multi-elements
  integer, parameter :: PDM_PART_TO_BLOCK_POST_MERGE                = 2  ! Merge multi-elements

interface

!>
!!
!! \brief Free a part to block structure
!!
!! \param [inout] ptb  Part to block structure
!!
!! \return       NULL
!!

subroutine PDM_part_to_block_free (ptb) &
bind (c, name='PDM_part_to_block_free')
  use iso_c_binding
  implicit none

  type(c_ptr), value :: ptb

end subroutine PDM_part_to_block_free


!>
!! \brief Return number of active ranks
!!
!! \param [in]   ptb          Part to block structure
!!
!! \return Number of active ranks
!!

function PDM_part_to_block_n_active_ranks_get (ptb) &
result(n_active_ranks)                              &
bind (c, name = 'PDM_part_to_block_n_active_ranks_get')
  use iso_c_binding
  implicit none

  type(c_ptr), value :: ptb
  integer(c_int)     :: n_active_ranks

end function PDM_part_to_block_n_active_ranks_get


!>
!! \brief Return if current rank is active
!!
!! \param [in]   ptb          Part to block structure
!!
!! \return 1 if current rank is active, 0 else
!!

function PDM_part_to_block_is_active_rank (ptb) &
result(is_active_rank)                          &
bind (c, name = 'PDM_part_to_block_is_active_rank')
  use iso_c_binding
  implicit none

  type(c_ptr), value :: ptb
  integer(c_int)     :: is_active_rank

end function PDM_part_to_block_is_active_rank


!>
!! \brief Return number of element in the current process
!!
!! \param [in]   ptb          Part to block structure
!!
!! \return Number of element in the current process
!!

function PDM_part_to_block_n_elt_block_get (ptb) &
result(n_elt_block)                              &
bind (c, name = 'PDM_part_to_block_n_elt_block_get')
  use iso_c_binding
  implicit none

  type(c_ptr), value :: ptb
  integer(c_int)     :: n_elt_block

end function PDM_part_to_block_n_elt_block_get


!>
!! \brief Return number of MPI ranks
!!
!! \param [in]   ptb          Part to block structure
!!
!! \return Number of MPI ranks
!!

function PDM_part_to_block_n_ranks_get (ptb) &
result(n_ranks)                              &
bind (c, name = 'PDM_part_to_block_n_ranks_get')
  use iso_c_binding
  implicit none

  type(c_ptr), value :: ptb
  integer(c_int)     :: n_ranks

end function PDM_part_to_block_n_ranks_get


!>
!!
!! \brief Return total number of element in the current process (summed over all partitions)
!!
!! \param [in]   ptb          Part to block structure
!!
!! \return Total number of element in the current process
!!

function PDM_part_to_block_n_elt_proc_get (ptb) &
result (n_elt_proc)                             &
bind (c, name='PDM_part_to_block_n_elt_proc_get')
  use iso_c_binding
  implicit none

  type(c_ptr), value :: ptb
  integer(c_int)     :: n_elt_proc

end function PDM_part_to_block_n_elt_proc_get


function PDM_part_to_block_exch_c (ptb,          &
                                   s_data,       &
                                   t_stride,     &
                                   cst_stride,   &
                                   part_stride,  &
                                   part_data,    &
                                   block_stride, &
                                   block_data)   &
result (s_block_data)                            &
bind (c, name='PDM_part_to_block_exch')
  use iso_c_binding
  implicit none

  type(c_ptr),    value :: ptb
  integer(c_int), value :: s_data
  integer(c_int), value :: t_stride
  integer(c_int), value :: cst_stride
  type(c_ptr),    value :: part_stride
  type(c_ptr),    value :: part_data
  type(c_ptr)           :: block_stride
  type(c_ptr)           :: block_data
  integer(c_int)        :: s_block_data

end function PDM_part_to_block_exch_c

end interface


interface PDM_part_to_block_exch
  module procedure PDM_part_to_block_exch_int
#ifdef PDM_LONG_G_NUM
  module procedure PDM_part_to_block_exch_g_num
#endif
  module procedure PDM_part_to_block_exch_double
  module procedure PDM_part_to_block_exch_complex4
  module procedure PDM_part_to_block_exch_complex8
  module procedure PDM_part_to_block_exch_real4
  module procedure PDM_part_to_block_exch_cptr
end interface


contains


!>
!!
!! \brief Create a partitioning to block redistribution
!!
!! \param [out]  ptb             Initialized PDM_part_to_block_t
!! \param [in]   t_distrib       Distribution type
!! \param [in]   t_post          Post processing type
!! \param [in]   partActiveNode  Part of active nodes (\ref PDM_writer_BLOCK_DISTRIB_PART_OF_NODE mode)
!! \param [in]   gnum_elt        Element global number
!! \param [in]   weight          Weight of elements (or NULL)
!! \param [in]   n_elt           Local number of elements
!! \param [in]   n_part          Number of partition
!! \param [in]   comm            MPI communicator
!!

subroutine PDM_part_to_block_create (ptb,            &
                                     t_distrib,      &
                                     t_post,         &
                                     partActiveNode, &
                                     gnum_elt,       &
                                     weight,         &
                                     n_elt,          &
                                     n_part,         &
                                     comm)
  use iso_c_binding
  implicit none

  type(c_ptr)                       :: ptb
  integer,          intent(in)      :: t_distrib
  integer,          intent(in)      :: t_post
  double precision, intent(in)      :: partActiveNode
  type(PDM_pointer_array_t), pointer :: gnum_elt
  type(PDM_pointer_array_t), pointer :: weight
  integer(pdm_l_num_s), pointer     :: n_elt(:)
  integer,          intent(in)      :: n_part
  integer,          intent(in)      :: comm

  integer(c_int)                    :: c_comm
  type(c_ptr)                       :: c_weight

  interface
    function PDM_part_to_block_create_c (t_distrib,      &
                                         t_post,         &
                                         partActiveNode, &
                                         gnum_elt,       &
                                         weight,         &
                                         n_elt,          &
                                         n_part,         &
                                         comm)           &
    result (ptb)                                         &
    bind (c, name='PDM_part_to_block_create')
      use iso_c_binding
      implicit none

      type(c_ptr)           :: ptb
      integer(c_int), value :: t_distrib
      integer(c_int), value :: t_post
      real(c_double), value :: partActiveNode
      type(c_ptr),    value :: gnum_elt
      type(c_ptr),    value :: weight
      type(c_ptr),    value :: n_elt
      integer(c_int), value :: n_part
      integer(c_int), value :: comm

    end function PDM_part_to_block_create_c
  end interface

  c_comm = PDM_MPI_Comm_f2c(comm)
  c_weight = C_NULL_PTR

  if (associated(weight)) then
   c_weight = c_loc(weight%cptr)
  endif

  ptb = PDM_part_to_block_create_c (t_distrib,            &
                                    t_post,               &
                                    partActiveNode,       &
                                    c_loc(gnum_elt%cptr), &
                                    c_weight,   &
                                    c_loc(n_elt),         &
                                    n_part,               &
                                    c_comm)

end subroutine PDM_part_to_block_create


!>
!!
!! \brief Create a partitioning to block redistribution
!!
!! \param [in]   t_distrib       Distribution type
!! \param [in]   t_post          Post processing type
!! \param [in]   partActiveNode  Part of active nodes (\ref PDM_writer_BLOCK_DISTRIB_PART_OF_NODE mode)
!! \param [in]   gnum_elt        Element global number
!! \param [in]   weight          Weight of elements (or NULL)
!! \param [in]   n_elt           Local number of elements
!! \param [in]   n_part          Number of partition
!! \param [in]   comm            MPI communicator
!!
!! \return   Initialized PDM_part_to_block_t
!!
!!

subroutine PDM_part_to_block_create_from_distrib (ptb,              &
                                                  t_distrib,        &
                                                  t_post,           &
                                                  partActiveNode,   &
                                                  gnum_elt,         &
                                                  dataDistribIndex, &
                                                  n_elt,            &
                                                  n_part,           &
                                                  comm)
  use iso_c_binding
  implicit none

  type(c_ptr)                       :: ptb
  integer,          intent(in)      :: t_distrib
  integer,          intent(in)      :: t_post
  double precision, intent(in)      :: partActiveNode
  type(PDM_pointer_array_t), pointer :: gnum_elt
  type(PDM_pointer_array_t), pointer :: dataDistribIndex
  integer(pdm_l_num_s), pointer     :: n_elt(:)
  integer,          intent(in)      :: n_part
  integer,          intent(in)      :: comm

  integer(c_int)                    :: c_comm

  interface
    function PDM_part_to_block_create_from_distrib_c (t_distrib,        &
                                                      t_post,           &
                                                      partActiveNode,   &
                                                      gnum_elt,         &
                                                      dataDistribIndex, &
                                                      n_elt,            &
                                                      n_part,           &
                                                      comm)             &
    result (ptb)                                                        &
    bind (c, name='PDM_part_to_block_create_from_distrib')
      use iso_c_binding
      implicit none

      type(c_ptr)           :: ptb
      integer(c_int), value :: t_distrib
      integer(c_int), value :: t_post
      real(c_double), value :: partActiveNode
      type(c_ptr),    value :: gnum_elt
      type(c_ptr),    value :: dataDistribIndex
      type(c_ptr),    value :: n_elt
      integer(c_int), value :: n_part
      integer(c_int), value :: comm

    end function PDM_part_to_block_create_from_distrib_c
  end interface

  c_comm = PDM_MPI_Comm_f2c(comm)

  ptb = PDM_part_to_block_create_from_distrib_c (t_distrib,                    &
                                                 t_post,                       &
                                                 partActiveNode,               &
                                                 c_loc(gnum_elt%cptr),         &
                                                 c_loc(dataDistribIndex%cptr), &
                                                 c_loc(n_elt),                 &
                                                 n_part,                       &
                                                 c_comm)

end subroutine PDM_part_to_block_create_from_distrib


!>
!! \brief Initialize a data exchange
!!
!! \param [in]   ptb          Part to block structure
!! \param [in]   s_data       Data size
!! \param [in]   t_stride     Stride type
!! \param [in]   var_stride   Variable stride (size = n_part) only for \ref PDM_writer_STRIDE_VAR
!! \param [in]   cst_stride   Stride only for \ref PDM_writer_STRIDE_CST
!!

subroutine PDM_part_to_block_exch_int (ptb,          &
                                       t_stride,     &
                                       cst_stride,   &
                                       part_stride,  &
                                       part_data,    &
                                       block_stride, &
                                       block_data)
  use iso_c_binding
  implicit none

  type(c_ptr), value                :: ptb
  integer, intent(in)               :: t_stride
  integer, intent(in)               :: cst_stride
  type(PDM_pointer_array_t), pointer :: part_stride
  type(PDM_pointer_array_t), pointer :: part_data
  integer(pdm_l_num_s), pointer     :: block_stride(:)
  integer(pdm_l_num_s), pointer     :: block_data(:)

  integer                           :: s_block_data
  type(c_ptr)                       :: c_block_stride 
  type(c_ptr)                       :: c_block_data   
  integer                           :: n_elt_block
  integer, parameter                :: s_data = 4

  type(c_ptr)                       :: c_part_stride 

  c_block_stride = C_NULL_PTR
  c_block_data   = C_NULL_PTR
  c_part_stride  = C_NULL_PTR

  if (part_data%type .ne. PDM_TYPE_int) then
    print *, "PDM_part_to_block_exch_int : wrong type"
    stop
  end if

  if (associated(part_stride)) then 
    c_part_stride = c_loc(part_stride%cptr)    
  endif

  s_block_data = PDM_part_to_block_exch_c (ptb,                     &
                                           s_data,                  &
                                           t_stride,                &
                                           cst_stride,              &
                                           c_part_stride,           &
                                           c_loc(part_data%cptr),   &
                                           c_block_stride,          &
                                           c_block_data)

  n_elt_block = PDM_part_to_block_n_elt_block_get (ptb)

  call c_f_pointer(c_block_stride, &
                   block_stride,   &
                   [n_elt_block])

  call c_f_pointer(c_block_data,   &
                   block_data,     &
                   [s_block_data])

end subroutine PDM_part_to_block_exch_int


#ifdef PDM_LONG_G_NUM
subroutine PDM_part_to_block_exch_g_num (ptb,          &
                                         t_stride,     &
                                         cst_stride,   &
                                         part_stride,  &
                                         part_data,    &
                                         block_stride, &
                                         block_data)
  use iso_c_binding
  implicit none

  type(c_ptr), value                :: ptb
  integer, intent(in)               :: t_stride
  integer, intent(in)               :: cst_stride
  type(PDM_pointer_array_t), pointer :: part_stride
  type(PDM_pointer_array_t), pointer :: part_data
  integer(pdm_l_num_s), pointer     :: block_stride(:)
  integer(pdm_g_num_s), pointer     :: block_data(:)

  integer                           :: s_block_data
  type(c_ptr)                       :: c_block_stride 
  type(c_ptr)                       :: c_block_data   
  integer                           :: n_elt_block
  integer, parameter                :: s_data = 8

  type(c_ptr)                       :: c_part_stride 

  c_block_stride = C_NULL_PTR
  c_block_data   = C_NULL_PTR
  c_part_stride = C_NULL_PTR

  if (part_data%type .ne. PDM_TYPE_G_NUM) then
    print *, "PDM_part_to_block_exch_g_num : wrong type"
    stop
  end if

  if (associated(part_stride)) then 
    c_part_stride = c_loc(part_stride%cptr)    
  endif

  s_block_data = PDM_part_to_block_exch_c (ptb,                     &
                                           s_data,                  &
                                           t_stride,                &
                                           cst_stride,              &
                                           c_part_stride,           &
                                           c_loc(part_data%cptr),   &
                                           c_block_stride,          &
                                           c_block_data)

  n_elt_block = PDM_part_to_block_n_elt_block_get (ptb)

  call c_f_pointer(c_block_stride, &
                   block_stride,   &
                   [n_elt_block])

  call c_f_pointer(c_block_data,   &
                   block_data,     &
                   [s_block_data])

end subroutine PDM_part_to_block_exch_g_num
#endif


subroutine PDM_part_to_block_exch_double (ptb,          &
                                          t_stride,     &
                                          cst_stride,   &
                                          part_stride,  &
                                          part_data,    &
                                          block_stride, &
                                          block_data)
  use iso_c_binding
  implicit none

  type(c_ptr), value                :: ptb
  integer, intent(in)               :: t_stride
  integer, intent(in)               :: cst_stride
  type(PDM_pointer_array_t), pointer :: part_stride
  type(PDM_pointer_array_t), pointer :: part_data
  integer(pdm_l_num_s), pointer     :: block_stride(:)
  double precision,     pointer     :: block_data(:)

  integer                           :: s_block_data
  type(c_ptr)                       :: c_block_stride
  type(c_ptr)                       :: c_block_data  
  integer                           :: n_elt_block
  integer, parameter                :: s_data = 8

  type(c_ptr)                       :: c_part_stride 

  c_block_stride = C_NULL_PTR
  c_block_data   = C_NULL_PTR
  c_part_stride = C_NULL_PTR

  if (part_data%type .ne. PDM_TYPE_DOUBLE) then
    print *, "PDM_part_to_block_exch_double : wrong type"
    stop
  end if

  if (associated(part_stride)) then 
    c_part_stride = c_loc(part_stride%cptr)    
  endif

  s_block_data = PDM_part_to_block_exch_c (ptb,                     &
                                           s_data,                  &
                                           t_stride,                &
                                           cst_stride,              &
                                           c_part_stride,           &
                                           c_loc(part_data%cptr),   &
                                           c_block_stride,          &
                                           c_block_data)

  n_elt_block = PDM_part_to_block_n_elt_block_get (ptb)

  call c_f_pointer(c_block_stride, &
                   block_stride,   &
                   [n_elt_block])

  call c_f_pointer(c_block_data,   &
                   block_data,     &
                   [s_block_data])

end subroutine PDM_part_to_block_exch_double



subroutine PDM_part_to_block_exch_complex4 (ptb,          &
                                          t_stride,     &
                                          cst_stride,   &
                                          part_stride,  &
                                          part_data,    &
                                          block_stride, &
                                          block_data)
  use iso_c_binding
  implicit none

  type(c_ptr), value                :: ptb
  integer, intent(in)               :: t_stride
  integer, intent(in)               :: cst_stride
  type(PDM_pointer_array_t), pointer :: part_stride
  type(PDM_pointer_array_t), pointer :: part_data
  integer(pdm_l_num_s), pointer     :: block_stride(:)
  complex(kind = 4),    pointer     :: block_data(:)

  integer                           :: s_block_data
  type(c_ptr)                       :: c_block_stride 
  type(c_ptr)                       :: c_block_data   
  integer                           :: n_elt_block
  integer, parameter                :: s_data = 8

  type(c_ptr)                       :: c_part_stride

  c_block_stride = C_NULL_PTR
  c_block_data   = C_NULL_PTR
  c_part_stride = C_NULL_PTR

  if (part_data%type .ne. PDM_TYPE_COMPLEX4) then
    print *, "PDM_part_to_block_exch_complex4 : wrong type"
    stop
  end if

  if (associated(part_stride)) then 
    c_part_stride = c_loc(part_stride%cptr)    
  endif

  s_block_data = PDM_part_to_block_exch_c (ptb,                     &
                                           s_data,                  &
                                           t_stride,                &
                                           cst_stride,              &
                                           c_part_stride,           &
                                           c_loc(part_data%cptr),   &
                                           c_block_stride,          &
                                           c_block_data)

  n_elt_block = PDM_part_to_block_n_elt_block_get (ptb)

  call c_f_pointer(c_block_stride, &
                   block_stride,   &
                   [n_elt_block])

  call c_f_pointer(c_block_data,   &
                   block_data,     &
                   [s_block_data])

end subroutine PDM_part_to_block_exch_complex4



subroutine PDM_part_to_block_exch_real4 (ptb,          &
                                          t_stride,     &
                                          cst_stride,   &
                                          part_stride,  &
                                          part_data,    &
                                          block_stride, &
                                          block_data)
  use iso_c_binding
  implicit none

  type(c_ptr), value                :: ptb
  integer, intent(in)               :: t_stride
  integer, intent(in)               :: cst_stride
  type(PDM_pointer_array_t), pointer :: part_stride
  type(PDM_pointer_array_t), pointer :: part_data
  integer(pdm_l_num_s), pointer     :: block_stride(:)
  real(kind = 4),       pointer     :: block_data(:)

  integer                           :: s_block_data
  type(c_ptr)                       :: c_block_stride = C_NULL_PTR
  type(c_ptr)                       :: c_block_data   = C_NULL_PTR
  integer                           :: n_elt_block
  integer, parameter                :: s_data = 4

  type(c_ptr)                       :: c_part_stride = C_NULL_PTR

  c_part_stride  = C_NULL_PTR
  c_block_stride = C_NULL_PTR
  c_block_data   = C_NULL_PTR

  if (part_data%type .ne. PDM_TYPE_REAL4) then
    print *, "PDM_part_to_block_exch_real4 : wrong type"
    stop
  end if

  if (associated(part_stride)) then 
    c_part_stride = c_loc(part_stride%cptr)    
  endif

  s_block_data = PDM_part_to_block_exch_c (ptb,                     &
                                           s_data,                  &
                                           t_stride,                &
                                           cst_stride,              &
                                           c_part_stride,           &
                                           c_loc(part_data%cptr),   &
                                           c_block_stride,          &
                                           c_block_data)

  n_elt_block = PDM_part_to_block_n_elt_block_get (ptb)

  call c_f_pointer(c_block_stride, &
                   block_stride,   &
                   [n_elt_block])

  call c_f_pointer(c_block_data,   &
                   block_data,     &
                   [s_block_data])

end subroutine PDM_part_to_block_exch_real4


subroutine PDM_part_to_block_exch_complex8 (ptb,          &
                                          t_stride,     &
                                          cst_stride,   &
                                          part_stride,  &
                                          part_data,    &
                                          block_stride, &
                                          block_data)
  use iso_c_binding
  implicit none

  type(c_ptr), value                :: ptb
  integer, intent(in)               :: t_stride
  integer, intent(in)               :: cst_stride
  type(PDM_pointer_array_t), pointer :: part_stride
  type(PDM_pointer_array_t), pointer :: part_data
  integer(pdm_l_num_s), pointer     :: block_stride(:)
  complex(kind = 8),    pointer     :: block_data(:)

  integer                           :: s_block_data
  type(c_ptr)                       :: c_block_stride = C_NULL_PTR
  type(c_ptr)                       :: c_block_data   = C_NULL_PTR
  integer                           :: n_elt_block
  integer, parameter                :: s_data = 16

  type(c_ptr)                       :: c_part_stride = C_NULL_PTR

  c_block_stride = C_NULL_PTR
  c_block_data   = C_NULL_PTR
  c_part_stride = C_NULL_PTR

  if (part_data%type .ne. PDM_TYPE_COMPLEX8) then
    print *, "PDM_part_to_block_exch_complex8 : wrong type"
    stop
  end if

  if (associated(part_stride)) then 
    c_part_stride = c_loc(part_stride%cptr)    
  endif

  s_block_data = PDM_part_to_block_exch_c (ptb,                     &
                                           s_data,                  &
                                           t_stride,                &
                                           cst_stride,              &
                                           c_part_stride,           &
                                           c_loc(part_data%cptr),   &
                                           c_block_stride,          &
                                           c_block_data)

  n_elt_block = PDM_part_to_block_n_elt_block_get (ptb)

  call c_f_pointer(c_block_stride, &
                   block_stride,   &
                   [n_elt_block])

  call c_f_pointer(c_block_data,   &
                   block_data,     &
                   [s_block_data])

end subroutine PDM_part_to_block_exch_complex8


subroutine PDM_part_to_block_exch_cptr (ptb,          &
                                        t_stride,     &
                                        cst_stride,   &
                                        part_stride,  &
                                        part_data,    &
                                        block_stride, &
                                        block_data)
  use iso_c_binding
  implicit none

  type(c_ptr), value                :: ptb
  integer, intent(in)               :: t_stride
  integer, intent(in)               :: cst_stride
  type(PDM_pointer_array_t), pointer :: part_stride
  type(PDM_pointer_array_t), pointer :: part_data
  integer(pdm_l_num_s), pointer     :: block_stride(:)
  type(c_ptr)                       :: block_data

  integer                           :: s_block_data
  type(c_ptr)                       :: c_block_stride 
  type(c_ptr)                       :: c_block_data   
  integer                           :: n_elt_block
  integer                           :: s_data

  type(c_ptr)                       :: c_part_stride

  s_data = part_data%s_data

  c_block_stride = C_NULL_PTR
  c_block_data   = C_NULL_PTR
  c_part_stride  = C_NULL_PTR

  if (part_data%type .ne. PDM_TYPE_CPTR) then
    print *, "PDM_part_to_block_exch_double : wrong type"
    stop
  end if

  if (associated(part_stride)) then 
    c_part_stride = c_loc(part_stride%cptr)    
  endif

  s_block_data = PDM_part_to_block_exch_c (ptb,                     &
                                           s_data,                  &
                                           t_stride,                &
                                           cst_stride,              &
                                           c_part_stride,           &
                                           c_loc(part_data%cptr),   &
                                           c_block_stride,          &
                                           c_block_data)

  n_elt_block = PDM_part_to_block_n_elt_block_get (ptb)

  call c_f_pointer(c_block_stride, &
                   block_stride,   &
                   [n_elt_block])

  block_data = c_block_data

end subroutine PDM_part_to_block_exch_cptr

!>
!! \brief Return active ranks
!!
!! \param [in]   ptb           Part to block structure
!! \param [out]  active_ranks  Active ranks
!!

subroutine PDM_part_to_block_active_ranks_get (ptb,          &
                                               active_ranks)
  use iso_c_binding
  implicit none

  type(c_ptr), value            :: ptb
  integer(pdm_l_num_s), pointer :: active_ranks(:)

  type(c_ptr)                   :: c_active_ranks
  integer                       :: n_active_ranks

  interface
    function PDM_part_to_block_active_ranks_get_c (ptb) &
    result(ptr_active_ranks)                            &
    bind (c, name = 'PDM_part_to_block_active_ranks_get')
      use iso_c_binding
      implicit none

      type(c_ptr), value :: ptb
      type(c_ptr)        :: ptr_active_ranks

    end function PDM_part_to_block_active_ranks_get_c
  end interface

  c_active_ranks = PDM_part_to_block_active_ranks_get_c (ptb)

  n_active_ranks = PDM_part_to_block_n_active_ranks_get (ptb)

  call c_f_pointer(c_active_ranks,   &
                   active_ranks,     &
                   [n_active_ranks])

end subroutine PDM_part_to_block_active_ranks_get


!>
!! \brief Return global numbers of element in the current process
!!
!! \param [in]   ptb          Part to block structure
!! \param [out]  g_nums       Global numbers
!!

subroutine PDM_part_to_block_block_gnum_get (ptb,    &
                                               g_nums)
  use iso_c_binding
  implicit none

  type(c_ptr), value            :: ptb
  integer(pdm_g_num_s), pointer :: g_nums(:)

  type(c_ptr)                   :: c_g_nums
  integer                       :: n_elt

  interface
    function PDM_part_to_block_block_gnum_get_c (ptb) &
    result(ptr_block_gnum)                            &
    bind (c, name = 'PDM_part_to_block_block_gnum_get')
      use iso_c_binding
      implicit none

      type(c_ptr), value :: ptb
      type(c_ptr)        :: ptr_block_gnum

    end function PDM_part_to_block_block_gnum_get_c
  end interface

  c_g_nums = PDM_part_to_block_block_gnum_get_c (ptb)

  n_elt = PDM_part_to_block_n_elt_block_get (ptb)

  call c_f_pointer(c_g_nums, &
                   g_nums,   &
                   [n_elt])

end subroutine PDM_part_to_block_block_gnum_get


!>
!! \brief Return block distribution index
!!
!! \param [in] ptb         Part to block structure
!!
!! \return  Distribution (size = communicator size + 1)

subroutine PDM_part_to_block_distrib_index_get (ptb,         &
                                                distrib_idx)
  use iso_c_binding
  implicit none

  type(c_ptr), value            :: ptb
  integer(pdm_g_num_s), pointer :: distrib_idx(:)

  type(c_ptr)                   :: c_distrib_idx
  integer                       :: n_rank

  interface
    function PDM_part_to_block_distrib_index_get_c (ptb) &
    result(ptr_distrib_index)                            &
    bind (c, name = 'PDM_part_to_block_distrib_index_get')
      use iso_c_binding
      implicit none

      type(c_ptr), value :: ptb
      type(c_ptr)        :: ptr_distrib_index

    end function PDM_part_to_block_distrib_index_get_c
  end interface

  c_distrib_idx = PDM_part_to_block_distrib_index_get_c (ptb)

  n_rank = PDM_part_to_block_n_ranks_get (ptb)

  call c_f_pointer(c_distrib_idx, &
                   distrib_idx,   &
                   [n_rank + 1])

end subroutine PDM_part_to_block_distrib_index_get


!>
!! \brief Return processus destination
!!
!! \param [in] ptb         Part to block structure
!!
!! \return  Destination (size = sum of partition elements)
!!

subroutine PDM_part_to_block_destination_get (ptb,         &
                                              destination)
  use iso_c_binding
  implicit none

  type(c_ptr), value            :: ptb
  integer(pdm_l_num_s), pointer :: destination(:)

  type(c_ptr)                   :: c_destination
  integer                       :: n_elt_proc

  interface
    function PDM_part_to_block_destination_get_c (ptb) &
    result(ptr_block_dest)                             &
    bind (c, name = 'PDM_part_to_block_destination_get')
      use iso_c_binding
      implicit none

      type(c_ptr), value :: ptb
      type(c_ptr)        :: ptr_block_dest

    end function PDM_part_to_block_destination_get_c
  end interface

  c_destination = PDM_part_to_block_destination_get_c (ptb)

  n_elt_proc = PDM_part_to_block_n_elt_proc_get (ptb)

  call c_f_pointer(c_destination, &
                   destination,   &
                   [n_elt_proc])

end subroutine PDM_part_to_block_destination_get

end module pdm_part_to_block
