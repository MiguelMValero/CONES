!-----------------------------------------------------------------------------
! This file is part of the ParaDiGM library.
!
! Copyright (C) 2020  ONERA
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

module pdm_block_to_part

  use pdm
  use pdm_pointer_array

  implicit none

interface

!>
!!
!! \brief Free a block to part structure
!!
!! \param [inout] btp  Block to part structure
!!
!! \return       NULL
!!

subroutine PDM_block_to_part_free (btp) &
bind (c, name='PDM_block_to_part_free')
  use iso_c_binding
  implicit none

  type(c_ptr), value :: btp

end subroutine PDM_block_to_part_free


!>
!!
!! \brief Get the number of partitions
!!
!! \param [in] btp         Block to part structure
!!
!! \return  Number of partitions
!!

function PDM_block_to_part_n_part_get (btp) &
result (n_part)                             &
bind (c, name='PDM_block_to_part_n_part_get')
  use iso_c_binding
  implicit none

  type(c_ptr), value :: btp
  integer(c_int)     :: n_part

end function PDM_block_to_part_n_part_get


!>
!!
!! \brief Get the number of elements in a given partition
!!
!! \param [in] btp         Block to part structure
!! \param [in] i_part      Id of current partition
!!
!! \return  Number of element in the current partition
!!

function PDM_block_to_part_n_elt_get (btp,    &
                                      i_part) &
result (n_elt)                                &
bind (c, name='PDM_block_to_part_n_elt_get')
  use iso_c_binding
  implicit none

  type(c_ptr),    value :: btp
  integer(c_int), value :: i_part
  integer(c_int)         :: n_elt

end function PDM_block_to_part_n_elt_get


subroutine PDM_block_to_part_exch_c (btp,          &
                                     s_data,       &
                                     t_stride,     &
                                     block_stride, &
                                     block_data,   &
                                     part_stride,  &
                                     part_data)    &
bind (c, name='PDM_block_to_part_exch')
  use iso_c_binding
  implicit none

  type(c_ptr),    value :: btp
  integer(c_int), value :: s_data
  integer(c_int), value :: t_stride
  type(c_ptr),    value :: block_stride
  type(c_ptr),    value :: block_data
  type(c_ptr)           :: part_stride
  type(c_ptr)           :: part_data

end subroutine PDM_block_to_part_exch_c


subroutine PDM_block_to_part_exch_in_place_c (btp,          &
                                              s_data,       &
                                              t_stride,     &
                                              block_stride, &
                                              block_data,   &
                                              part_stride,  &
                                              part_data)    &
bind (c, name='PDM_block_to_part_exch_in_place')
  use iso_c_binding
  implicit none

  type(c_ptr),    value :: btp
  integer(c_int), value :: s_data
  integer(c_int), value :: t_stride
  type(c_ptr),    value :: block_stride
  type(c_ptr),    value :: block_data
  type(c_ptr),    value :: part_stride
  type(c_ptr),    value :: part_data

end subroutine PDM_block_to_part_exch_in_place_c


end interface


interface PDM_block_to_part_exch
  module procedure PDM_block_to_part_exch_int
#ifdef PDM_LONG_G_NUM
  module procedure PDM_block_to_part_exch_g_num
#endif
  module procedure PDM_block_to_part_exch_double
  module procedure PDM_block_to_part_exch_real4
  module procedure PDM_block_to_part_exch_complex4
  module procedure PDM_block_to_part_exch_complex8
  module procedure PDM_block_to_part_exch_cptr
end interface

interface PDM_block_to_part_exch_in_place
  module procedure PDM_block_to_part_exch_in_place_int
#ifdef PDM_LONG_G_NUM
  module procedure PDM_block_to_part_exch_in_place_g_num
#endif
  module procedure PDM_block_to_part_exch_in_place_double
  module procedure PDM_block_to_part_exch_in_place_real4
  module procedure PDM_block_to_part_exch_in_place_complex4
  module procedure PDM_block_to_part_exch_in_place_complex8
  module procedure PDM_block_to_part_exch_in_place_cptr
end interface

private :: PDM_block_to_part_exch_finalize
private :: PDM_block_to_part_exch_in_place_finalize

contains

!>
!!
!! \brief Create a block to partitions redistribution
!!
!! \param [out]  btp               Initialized \ref PDM_block_to_part instance
!! \param [in]   block_distrib_idx Block distribution (size : \ref size of \ref comm + 1)
!! \param [in]   gnum_elt          Element global number (size : \ref n_part)
!! \param [in]   n_elt             Local number of elements (size : \ref n_part)
!! \param [in]   n_part            Number of partition
!! \param [in]   comm              MPI communicator
!!

subroutine PDM_block_to_part_create (btp,               &
                                     block_distrib_idx, &
                                     gnum_elt,          &
                                     n_elt,             &
                                     n_part,            &
                                     comm)
  use iso_c_binding
  implicit none

  type(c_ptr)                       :: btp
  integer(pdm_g_num_s), pointer     :: block_distrib_idx(:)
  type(PDM_pointer_array_t), pointer :: gnum_elt
  integer(pdm_l_num_s), pointer     :: n_elt(:)
  integer, intent(in)               :: n_part
  integer, intent(in)               :: comm

  integer(c_int)                    :: c_comm

  interface
    function PDM_block_to_part_create_c (block_distrib_idx, &
                                         gnum_elt,          &
                                         n_elt,             &
                                         n_part,            &
                                         comm)              &
    result (btp)                                            &
    bind (c, name='PDM_block_to_part_create')
      use iso_c_binding
      implicit none

      type(c_ptr)           :: btp
      type(c_ptr),    value :: block_distrib_idx
      type(c_ptr),    value :: gnum_elt
      type(c_ptr),    value :: n_elt
      integer(c_int), value :: n_part
      integer(c_int), value :: comm

    end function PDM_block_to_part_create_c
  end interface

  c_comm = PDM_MPI_Comm_f2c(comm)

  btp = PDM_block_to_part_create_c (c_loc(block_distrib_idx), &
                                    c_loc(gnum_elt%cptr),     &
                                    c_loc(n_elt),             &
                                    n_part,                   &
                                    c_comm)

end subroutine PDM_block_to_part_create


!>
!!
!! \brief Initialize an exchange
!!
!! \param [in]   btp          Block to part structure
!! \param [in]   t_stride     Stride type
!! \param [in]   block_stride Stride for each block element for \ref PDM_STRIDE_VAR
!!                            Constant stride for \ref PDM_STRIDE_VAR
!! \param [in]   block_data   Block data
!! \param [out]  part_stride  Partition stride
!! \param [out]  part_data    Partition data
!!

subroutine PDM_block_to_part_exch_in_place_int (btp,          &
                                                t_stride,     &
                                                block_stride, &
                                                block_data,   &
                                                part_stride,  &
                                                part_data)
  use iso_c_binding
  implicit none

  type(c_ptr), value                :: btp
  integer, intent(in)               :: t_stride
  integer(pdm_l_num_s), pointer     :: block_stride(:)
  integer(pdm_l_num_s), pointer     :: block_data(:)
  type(PDM_pointer_array_t), pointer :: part_stride
  type(PDM_pointer_array_t), pointer :: part_data

  integer, parameter                 :: s_data = 4

  type(c_ptr) :: c_block_stride     
  type(c_ptr) :: c_block_data       
  type(c_ptr) :: c_part_stride_cptr 
  type(c_ptr) :: c_part_data_cptr

  if (part_data%type .ne. PDM_TYPE_INT) then
    print *, "PDM_block_to_part_exch_in_place_int : wrong type"
    stop
  end if
  
  c_block_stride = C_NULL_PTR
  if (associated(block_stride)) then
    c_block_stride = c_loc(block_stride)
  endif 
    
  c_block_data = C_NULL_PTR
  if (associated(block_data)) then
    c_block_data = c_loc(block_data)
  endif 
    
  c_part_stride_cptr = C_NULL_PTR
  if (associated(part_stride%cptr)) then
    c_part_stride_cptr = c_loc(part_stride%cptr)
  endif 
    
  c_part_data_cptr = C_NULL_PTR
  if (associated(part_data%cptr)) then
    c_part_data_cptr = c_loc(part_data%cptr)
  endif 

  call PDM_block_to_part_exch_in_place_c (btp,                     &
                                          s_data,                  &
                                          t_stride,                &
                                          c_block_stride,          &
                                          c_block_data,            &
                                          c_part_stride_cptr,      &
                                          c_part_data_cptr)

  call PDM_block_to_part_exch_in_place_finalize (btp,           &
                                                 t_stride,      &
                                                 block_stride,  &
                                                 part_stride,   &
                                                 part_data)

end subroutine PDM_block_to_part_exch_in_place_int


#ifdef PDM_LONG_G_NUM
subroutine PDM_block_to_part_exch_in_place_g_num (btp,          &
                                                  t_stride,     &
                                                  block_stride, &
                                                  block_data,   &
                                                  part_stride,  &
                                                  part_data)
  use iso_c_binding
  implicit none

  type(c_ptr), value                :: btp
  integer, intent(in)               :: t_stride
  integer(pdm_l_num_s), pointer     :: block_stride(:)
  integer(pdm_g_num_s), pointer     :: block_data(:)
  type(PDM_pointer_array_t), pointer :: part_stride
  type(PDM_pointer_array_t), pointer :: part_data

  integer, parameter                 :: s_data = 8

  type(c_ptr) :: c_block_stride     
  type(c_ptr) :: c_block_data       
  type(c_ptr) :: c_part_stride_cptr 
  type(c_ptr) :: c_part_data_cptr

  if (part_data%type .ne. PDM_TYPE_G_NUM) then
    print *, "PDM_block_to_part_exch_in_place_g_num : wrong type"
    stop
  end if

  c_block_stride = C_NULL_PTR
  if (associated(block_stride)) then
    c_block_stride = c_loc(block_stride)
  endif 
    
  c_block_data = C_NULL_PTR
  if (associated(block_data)) then
    c_block_data = c_loc(block_data)
  endif 
    
  c_part_stride_cptr = C_NULL_PTR
  if (associated(part_stride%cptr)) then
    c_part_stride_cptr = c_loc(part_stride%cptr)
  endif 
    
  c_part_data_cptr = C_NULL_PTR
  if (associated(part_data%cptr)) then
    c_part_data_cptr = c_loc(part_data%cptr)
  endif 

  call PDM_block_to_part_exch_in_place_c (btp,                     &
                                          s_data,                  &
                                          t_stride,                &
                                          c_block_stride,     &
                                          c_block_data,       &
                                          c_part_stride_cptr, &
                                          c_part_data_cptr)

  call PDM_block_to_part_exch_in_place_finalize (btp,           &
                                                 t_stride,      &
                                                 block_stride,  &
                                                 part_stride,   &
                                                 part_data)

end subroutine PDM_block_to_part_exch_in_place_g_num
#endif


subroutine PDM_block_to_part_exch_in_place_double (btp,          &
                                                   t_stride,     &
                                                   block_stride, &
                                                   block_data,   &
                                                   part_stride,  &
                                                   part_data)
  use iso_c_binding
  implicit none

  type(c_ptr), value                :: btp
  integer, intent(in)               :: t_stride
  integer(pdm_l_num_s), pointer     :: block_stride(:)
  double precision,     pointer     :: block_data(:)
  type(PDM_pointer_array_t), pointer :: part_stride
  type(PDM_pointer_array_t), pointer :: part_data

  integer, parameter                 :: s_data = 8

  type(c_ptr) :: c_block_stride     
  type(c_ptr) :: c_block_data       
  type(c_ptr) :: c_part_stride_cptr 
  type(c_ptr) :: c_part_data_cptr

  if (part_data%type .ne. PDM_TYPE_DOUBLE) then
    print *, "PDM_block_to_part_exch_in_place_double : wrong type"
    stop
  end if

  c_block_stride = C_NULL_PTR
  if (associated(block_stride)) then
    c_block_stride = c_loc(block_stride)
  endif 
    
  c_block_data = C_NULL_PTR
  if (associated(block_data)) then
    c_block_data = c_loc(block_data)
  endif 
    
  c_part_stride_cptr = C_NULL_PTR
  if (associated(part_stride%cptr)) then
    c_part_stride_cptr = c_loc(part_stride%cptr)
  endif 
    
  c_part_data_cptr = C_NULL_PTR
  if (associated(part_data%cptr)) then
    c_part_data_cptr = c_loc(part_data%cptr)
  endif 

  call PDM_block_to_part_exch_in_place_c (btp,                     &
                                          s_data,                  &
                                          t_stride,                &
                                          c_block_stride,     &
                                          c_block_data,       &
                                          c_part_stride_cptr, &
                                          c_part_data_cptr)

  call PDM_block_to_part_exch_in_place_finalize (btp,           &
                                                 t_stride,      &
                                                 block_stride,  &
                                                 part_stride,   &
                                                 part_data)

end subroutine PDM_block_to_part_exch_in_place_double


subroutine PDM_block_to_part_exch_in_place_real4 (btp,          &
                                                   t_stride,     &
                                                   block_stride, &
                                                   block_data,   &
                                                   part_stride,  &
                                                   part_data)
  use iso_c_binding
  implicit none

  type(c_ptr), value                :: btp
  integer, intent(in)               :: t_stride
  integer(pdm_l_num_s), pointer     :: block_stride(:)
  real (kind = 4),      pointer     :: block_data(:)
  type(PDM_pointer_array_t), pointer :: part_stride
  type(PDM_pointer_array_t), pointer :: part_data

  integer, parameter                 :: s_data = 4

  type(c_ptr) :: c_block_stride     
  type(c_ptr) :: c_block_data       
  type(c_ptr) :: c_part_stride_cptr 
  type(c_ptr) :: c_part_data_cptr

  if (part_data%type .ne. PDM_TYPE_REAL4) then
    print *, "PDM_block_to_part_exch_in_place_real4 : wrong type"
    stop
  end if

  c_block_stride = C_NULL_PTR
  if (associated(block_stride)) then
    c_block_stride = c_loc(block_stride)
  endif 
    
  c_block_data = C_NULL_PTR
  if (associated(block_data)) then
    c_block_data = c_loc(block_data)
  endif 
    
  c_part_stride_cptr = C_NULL_PTR
  if (associated(part_stride%cptr)) then
    c_part_stride_cptr = c_loc(part_stride%cptr)
  endif 
    
  c_part_data_cptr = C_NULL_PTR
  if (associated(part_data%cptr)) then
    c_part_data_cptr = c_loc(part_data%cptr)
  endif 


  call PDM_block_to_part_exch_in_place_c (btp,                     &
                                          s_data,                  &
                                          t_stride,                &
                                          c_block_stride,     &
                                          c_block_data,       &
                                          c_part_stride_cptr, &
                                          c_part_data_cptr)

  call PDM_block_to_part_exch_in_place_finalize (btp,           &
                                                 t_stride,      &
                                                 block_stride,  &
                                                 part_stride,   &
                                                 part_data)

end subroutine PDM_block_to_part_exch_in_place_real4


subroutine PDM_block_to_part_exch_in_place_complex4 (btp,          &
                                                   t_stride,     &
                                                   block_stride, &
                                                   block_data,   &
                                                   part_stride,  &
                                                   part_data)
  use iso_c_binding
  implicit none

  type(c_ptr), value                :: btp
  integer, intent(in)               :: t_stride
  integer(pdm_l_num_s), pointer     :: block_stride(:)
  complex (kind = 4),      pointer :: block_data(:)
  type(PDM_pointer_array_t), pointer :: part_stride
  type(PDM_pointer_array_t), pointer :: part_data

  integer, parameter                 :: s_data = 8

  type(c_ptr) :: c_block_stride     
  type(c_ptr) :: c_block_data       
  type(c_ptr) :: c_part_stride_cptr 
  type(c_ptr) :: c_part_data_cptr

  if (part_data%type .ne. PDM_TYPE_REAL4) then
    print *, "PDM_block_to_part_exch_in_place_complex4 : wrong type"
    stop
  end if

  c_block_stride = C_NULL_PTR
  if (associated(block_stride)) then
    c_block_stride = c_loc(block_stride)
  endif 
    
  c_block_data = C_NULL_PTR
  if (associated(block_data)) then
    c_block_data = c_loc(block_data)
  endif 
    
  c_part_stride_cptr = C_NULL_PTR
  if (associated(part_stride%cptr)) then
    c_part_stride_cptr = c_loc(part_stride%cptr)
  endif 
    
  c_part_data_cptr = C_NULL_PTR
  if (associated(part_data%cptr)) then
    c_part_data_cptr = c_loc(part_data%cptr)
  endif 

  call PDM_block_to_part_exch_in_place_c (btp,                     &
                                          s_data,                  &
                                          t_stride,                &
                                          c_block_stride,     &
                                          c_block_data,       &
                                          c_part_stride_cptr, &
                                          c_part_data_cptr)

  call PDM_block_to_part_exch_in_place_finalize (btp,           &
                                                 t_stride,      &
                                                 block_stride,  &
                                                 part_stride,   &
                                                 part_data)

end subroutine PDM_block_to_part_exch_in_place_complex4


subroutine PDM_block_to_part_exch_in_place_complex8 (btp,          &
                                                   t_stride,     &
                                                   block_stride, &
                                                   block_data,   &
                                                   part_stride,  &
                                                   part_data)
  use iso_c_binding
  implicit none

  type(c_ptr), value                :: btp
  integer, intent(in)               :: t_stride
  integer(pdm_l_num_s), pointer     :: block_stride(:)
  complex (kind = 8),      pointer :: block_data(:)
  type(PDM_pointer_array_t), pointer :: part_stride
  type(PDM_pointer_array_t), pointer :: part_data

  integer, parameter                 :: s_data = 16

  type(c_ptr) :: c_block_stride     
  type(c_ptr) :: c_block_data       
  type(c_ptr) :: c_part_stride_cptr 
  type(c_ptr) :: c_part_data_cptr

  if (part_data%type .ne. PDM_TYPE_REAL4) then
    print *, "PDM_block_to_part_exch_in_place_complex4 : wrong type"
    stop
  end if

  c_block_stride = C_NULL_PTR
  if (associated(block_stride)) then
    c_block_stride = c_loc(block_stride)
  endif 
    
  c_block_data = C_NULL_PTR
  if (associated(block_data)) then
    c_block_data = c_loc(block_data)
  endif 
    
  c_part_stride_cptr = C_NULL_PTR
  if (associated(part_stride%cptr)) then
    c_part_stride_cptr = c_loc(part_stride%cptr)
  endif 
    
  c_part_data_cptr = C_NULL_PTR
  if (associated(part_data%cptr)) then
    c_part_data_cptr = c_loc(part_data%cptr)
  endif 

  call PDM_block_to_part_exch_in_place_c (btp,                     &
                                          s_data,                  &
                                          t_stride,                &
                                          c_block_stride,     &
                                          c_block_data,       &
                                          c_part_stride_cptr, &
                                          c_part_data_cptr)

  call PDM_block_to_part_exch_in_place_finalize (btp,           &
                                                 t_stride,      &
                                                 block_stride,  &
                                                 part_stride,   &
                                                 part_data)

end subroutine PDM_block_to_part_exch_in_place_complex8


subroutine PDM_block_to_part_exch_in_place_cptr (btp,          &
                                                   t_stride,     &
                                                   block_stride, &
                                                   block_data,   &
                                                   part_stride,  &
                                                   part_data)
  use iso_c_binding
  implicit none

  type(c_ptr), value                :: btp
  integer, intent(in)               :: t_stride
  integer(pdm_l_num_s), pointer     :: block_stride(:)
  type (c_ptr),      pointer        :: block_data(:)
  type(PDM_pointer_array_t), pointer :: part_stride
  type(PDM_pointer_array_t), pointer :: part_data

  integer                           :: s_data

  type(c_ptr) :: c_block_stride     
  type(c_ptr) :: c_block_data       
  type(c_ptr) :: c_part_stride_cptr 
  type(c_ptr) :: c_part_data_cptr

  s_data = part_data%s_data

  if (part_data%type .ne. PDM_TYPE_REAL4) then
    print *, "PDM_block_to_part_exch_in_place_complex4 : wrong type"
    stop
  end if

  c_block_stride = C_NULL_PTR
  if (associated(block_stride)) then
    c_block_stride = c_loc(block_stride)
  endif 
    
  c_block_data = C_NULL_PTR
  if (associated(block_data)) then
    c_block_data = c_loc(block_data)
  endif 
    
  c_part_stride_cptr = C_NULL_PTR
  if (associated(part_stride%cptr)) then
    c_part_stride_cptr = c_loc(part_stride%cptr)
  endif 
    
  c_part_data_cptr = C_NULL_PTR
  if (associated(part_data%cptr)) then
    c_part_data_cptr = c_loc(part_data%cptr)
  endif 

  call PDM_block_to_part_exch_in_place_c (btp,                     &
                                          s_data,                  &
                                          t_stride,                &
                                          c_block_stride,     &
                                          c_block_data,       &
                                          c_part_stride_cptr, &
                                          c_part_data_cptr)

  call PDM_block_to_part_exch_in_place_finalize (btp,           &
                                                 t_stride,      &
                                                 block_stride,  &
                                                 part_stride,   &
                                                 part_data)

end subroutine PDM_block_to_part_exch_in_place_cptr


!>
!!
!! \brief Initialize an exchange
!! (part_stride and part_data are allocated in function)
!!
!! \param [in]   btp          Block to part structure
!! \param [in]   t_stride     Stride type
!! \param [in]   block_stride Stride for each block element for \ref PDM_STRIDE_VAR
!!                            Constant stride for \ref PDM_STRIDE_VAR
!! \param [in]   block_data   Block data
!! \param [out]  part_stride  Partition stride
!! \param [out]  part_data    Partition data
!!
!!

subroutine PDM_block_to_part_exch_int (btp,          &
                                       t_stride,     &
                                       block_stride, &
                                       block_data,   &
                                       part_stride,  &
                                       part_data)
  use iso_c_binding
  implicit none

  type(c_ptr), value                :: btp
  integer, intent(in)               :: t_stride
  integer(pdm_l_num_s), pointer     :: block_stride(:)
  integer(pdm_l_num_s), pointer     :: block_data(:)
  type(PDM_pointer_array_t), pointer :: part_stride
  type(PDM_pointer_array_t), pointer :: part_data

  type(c_ptr)                       :: c_block_stride     
  type(c_ptr)                       :: c_block_data       
  type(c_ptr)                       :: c_part_stride
  type(c_ptr)                       :: c_part_data  

  integer, parameter                :: s_data = 4

  c_part_stride = C_NULL_PTR
  c_part_data   = C_NULL_PTR


  c_block_stride = C_NULL_PTR
  if (associated(block_stride)) then
    c_block_stride = c_loc(block_stride)
  endif 
    
  c_block_data = C_NULL_PTR
  if (associated(block_data)) then
    c_block_data = c_loc(block_data)
  endif 

  call PDM_block_to_part_exch_c (btp,                 &
                                 s_data,              &
                                 t_stride,            &
                                 c_block_stride, &
                                 c_block_data,   &
                                 c_part_stride,       &
                                 c_part_data)

  call PDM_block_to_part_exch_finalize (btp,           &
                                        t_stride,      &
                                        block_stride,  &
                                        PDM_TYPE_INT,  &
                                        c_part_stride, &
                                        c_part_data,   &
                                        part_stride,   &
                                        part_data)

end subroutine PDM_block_to_part_exch_int


#ifdef PDM_LONG_G_NUM
subroutine PDM_block_to_part_exch_g_num (btp,          &
                                         t_stride,     &
                                         block_stride, &
                                         block_data,   &
                                         part_stride,  &
                                         part_data)
  use iso_c_binding
  implicit none

  type(c_ptr), value                :: btp
  integer, intent(in)               :: t_stride
  integer(pdm_l_num_s), pointer     :: block_stride(:)
  integer(pdm_g_num_s), pointer     :: block_data(:)
  type(PDM_pointer_array_t), pointer :: part_stride
  type(PDM_pointer_array_t), pointer :: part_data

  type(c_ptr)                       :: c_block_stride     
  type(c_ptr)                       :: c_block_data       
  type(c_ptr)                       :: c_part_stride
  type(c_ptr)                       :: c_part_data  

  integer, parameter                :: s_data = 8

  c_part_stride = C_NULL_PTR
  c_part_data   = C_NULL_PTR

  c_block_stride = C_NULL_PTR
  if (associated(block_stride)) then
    c_block_stride = c_loc(block_stride)
  endif 
    
  c_block_data = C_NULL_PTR
  if (associated(block_data)) then
    c_block_data = c_loc(block_data)
  endif 

  call PDM_block_to_part_exch_c (btp,                 &
                                 s_data,              &
                                 t_stride,            &
                                 c_block_stride, &
                                 c_block_data,   &
                                 c_part_stride,       &
                                 c_part_data)

  call PDM_block_to_part_exch_finalize (btp,            &
                                        t_stride,       &
                                        block_stride,   &
                                        PDM_TYPE_G_NUM, &
                                        c_part_stride,  &
                                        c_part_data,    &
                                        part_stride,    &
                                        part_data)

end subroutine PDM_block_to_part_exch_g_num
#endif

subroutine PDM_block_to_part_exch_double (btp,          &
                                          t_stride,     &
                                          block_stride, &
                                          block_data,   &
                                          part_stride,  &
                                          part_data)
  use iso_c_binding
  implicit none

  type(c_ptr), value                :: btp
  integer, intent(in)               :: t_stride
  integer(pdm_l_num_s), pointer     :: block_stride(:)
  double precision,     pointer     :: block_data(:)
  type(PDM_pointer_array_t), pointer :: part_stride
  type(PDM_pointer_array_t), pointer :: part_data

  type(c_ptr)                       :: c_block_stride     
  type(c_ptr)                       :: c_block_data       
  type(c_ptr)                       :: c_part_stride 
  type(c_ptr)                       :: c_part_data   

  integer, parameter                :: s_data = 8

  c_part_stride = C_NULL_PTR
  c_part_data   = C_NULL_PTR

  c_block_stride = C_NULL_PTR
  if (associated(block_stride)) then
    c_block_stride = c_loc(block_stride)
  endif 
    
  c_block_data = C_NULL_PTR
  if (associated(block_data)) then
    c_block_data = c_loc(block_data)
  endif 

  call PDM_block_to_part_exch_c (btp,                 &
                                 s_data,              &
                                 t_stride,            &
                                 c_block_stride,      &
                                 c_block_data,        &
                                 c_part_stride,       &
                                 c_part_data)

  call PDM_block_to_part_exch_finalize (btp,             &
                                        t_stride,        &
                                        block_stride,    &
                                        PDM_TYPE_DOUBLE, &
                                        c_part_stride,   &
                                        c_part_data,     &
                                        part_stride,     &
                                        part_data)

end subroutine PDM_block_to_part_exch_double

subroutine PDM_block_to_part_exch_real4 (btp,          &
                                          t_stride,     &
                                          block_stride, &
                                          block_data,   &
                                          part_stride,  &
                                          part_data)
  use iso_c_binding
  implicit none

  type(c_ptr), value                :: btp
  integer, intent(in)               :: t_stride
  integer(pdm_l_num_s), pointer     :: block_stride(:)
  real (kind = 4),      pointer     :: block_data(:)
  type(PDM_pointer_array_t), pointer :: part_stride
  type(PDM_pointer_array_t), pointer :: part_data

  type(c_ptr)                       :: c_block_stride     
  type(c_ptr)                       :: c_block_data       
  type(c_ptr)                       :: c_part_stride 
  type(c_ptr)                       :: c_part_data   

  integer, parameter                :: s_data = 4

  c_part_stride = C_NULL_PTR
  c_part_data   = C_NULL_PTR

  c_block_stride = C_NULL_PTR
  if (associated(block_stride)) then
    c_block_stride = c_loc(block_stride)
  endif 
    
  c_block_data = C_NULL_PTR
  if (associated(block_data)) then
    c_block_data = c_loc(block_data)
  endif 
  call PDM_block_to_part_exch_c (btp,                 &
                                 s_data,              &
                                 t_stride,            &
                                 c_block_stride,      &
                                 c_block_data,        &
                                 c_part_stride,       &
                                 c_part_data)

  call PDM_block_to_part_exch_finalize (btp,             &
                                        t_stride,        &
                                        block_stride,    &
                                        PDM_TYPE_REAL4, &
                                        c_part_stride,   &
                                        c_part_data,     &
                                        part_stride,     &
                                        part_data)

end subroutine PDM_block_to_part_exch_real4

subroutine PDM_block_to_part_exch_complex4 (btp,          &
                                          t_stride,     &
                                          block_stride, &
                                          block_data,   &
                                          part_stride,  &
                                          part_data)
  use iso_c_binding
  implicit none

  type(c_ptr), value                :: btp
  integer, intent(in)               :: t_stride
  integer(pdm_l_num_s), pointer     :: block_stride(:)
  complex (kind = 4),     pointer     :: block_data(:)
  type(PDM_pointer_array_t), pointer :: part_stride
  type(PDM_pointer_array_t), pointer :: part_data

  type(c_ptr)                       :: c_block_stride     
  type(c_ptr)                       :: c_block_data       
  type(c_ptr)                       :: c_part_stride 
  type(c_ptr)                       :: c_part_data   

  integer, parameter                :: s_data = 8

  c_part_stride = C_NULL_PTR
  c_part_data   = C_NULL_PTR

  c_block_stride = C_NULL_PTR
  if (associated(block_stride)) then
    c_block_stride = c_loc(block_stride)
  endif 
    
  c_block_data = C_NULL_PTR
  if (associated(block_data)) then
    c_block_data = c_loc(block_data)
  endif 

  call PDM_block_to_part_exch_c (btp,                 &
                                 s_data,              &
                                 t_stride,            &
                                 c_block_stride,      &
                                 c_block_data,        &
                                 c_part_stride,       &
                                 c_part_data)

  call PDM_block_to_part_exch_finalize (btp,             &
                                        t_stride,        &
                                        block_stride,    &
                                        PDM_TYPE_COMPLEX4, &
                                        c_part_stride,   &
                                        c_part_data,     &
                                        part_stride,     &
                                        part_data)

end subroutine PDM_block_to_part_exch_complex4

subroutine PDM_block_to_part_exch_complex8 (btp,          &
                                          t_stride,     &
                                          block_stride, &
                                          block_data,   &
                                          part_stride,  &
                                          part_data)
  use iso_c_binding
  implicit none

  type(c_ptr), value                :: btp
  integer, intent(in)               :: t_stride
  integer(pdm_l_num_s), pointer     :: block_stride(:)
  complex (kind = 8),     pointer     :: block_data(:)
  type(PDM_pointer_array_t), pointer :: part_stride
  type(PDM_pointer_array_t), pointer :: part_data

  type(c_ptr)                       :: c_block_stride     
  type(c_ptr)                       :: c_block_data       
  type(c_ptr)                       :: c_part_stride 
  type(c_ptr)                       :: c_part_data   

  integer, parameter                :: s_data = 16

  c_part_stride = C_NULL_PTR
  c_part_data   = C_NULL_PTR

  c_block_stride = C_NULL_PTR
  if (associated(block_stride)) then
    c_block_stride = c_loc(block_stride)
  endif 
    
  c_block_data = C_NULL_PTR
  if (associated(block_data)) then
    c_block_data = c_loc(block_data)
  endif 

  call PDM_block_to_part_exch_c (btp,                 &
                                 s_data,              &
                                 t_stride,            &
                                 c_block_stride,      &
                                 c_block_data,        &
                                 c_part_stride,       &
                                 c_part_data)

  call PDM_block_to_part_exch_finalize (btp,             &
                                        t_stride,        &
                                        block_stride,    &
                                        PDM_TYPE_COMPLEX8, &
                                        c_part_stride,   &
                                        c_part_data,     &
                                        part_stride,     &
                                        part_data)

end subroutine PDM_block_to_part_exch_complex8

subroutine PDM_block_to_part_exch_cptr (btp,          &
                                        s_data,       &
                                        t_stride,     &
                                        block_stride, &
                                        block_data,   &
                                        part_stride,  &
                                        part_data)
  use iso_c_binding
  implicit none

  type(c_ptr), value                :: btp
  integer, intent(in)               :: s_data
  integer, intent(in)               :: t_stride
  integer(pdm_l_num_s), pointer     :: block_stride(:)
  type (c_ptr),         pointer     :: block_data(:)
  type(PDM_pointer_array_t), pointer :: part_stride
  type(PDM_pointer_array_t), pointer :: part_data

  type(c_ptr)                       :: c_block_stride     
  type(c_ptr)                       :: c_block_data       
  type(c_ptr)                       :: c_part_stride 
  type(c_ptr)                       :: c_part_data   

  c_part_stride = C_NULL_PTR
  c_part_data   = C_NULL_PTR

  c_block_stride = C_NULL_PTR
  if (associated(block_stride)) then
    c_block_stride = c_loc(block_stride)
  endif 
    
  c_block_data = C_NULL_PTR
  if (associated(block_data)) then
    c_block_data = c_loc(block_data)
  endif 

  call PDM_block_to_part_exch_c (btp,                 &
                                 s_data,              &
                                 t_stride,            &
                                 c_block_stride,      &
                                 c_block_data,        &
                                 c_part_stride,       &
                                 c_part_data)

  call PDM_block_to_part_exch_finalize (btp,             &
                                        t_stride,        &
                                        block_stride,    &
                                        PDM_TYPE_CPTR, &
                                        c_part_stride,   &
                                        c_part_data,     &
                                        part_stride,     &
                                        part_data,       &
                                        s_data)

end subroutine PDM_block_to_part_exch_cptr


!>
!!
!! \brief Return index in the block for a gnum
!!
!! \param [in]  ptb         Part to block structure
!! \param [in]  gNum        Global number
!! \param [out] idx         Index
!!

subroutine PDM_block_to_part_gnum_idx_get (btp,  &
                                           gNum, &
                                           idx)
  use iso_c_binding
  implicit none

  type(c_ptr), value                :: btp
  integer(pdm_g_num_s), intent(in)  :: gNum
  integer(pdm_l_num_s), intent(out) :: idx

  interface
    function PDM_block_to_part_gnum_idx_get_c (btp,  &
                                               gNum) &
    result (idx)                                     &
    bind (c, name='PDM_block_to_part_gnum_idx_get')
      use iso_c_binding
      implicit none

      type(c_ptr),     value :: btp
#ifdef PDM_LONG_G_NUM
      integer(c_long), value :: gNum
#else
      integer(c_int),  value :: gNum
#endif
      integer(c_int)         :: idx
    end function PDM_block_to_part_gnum_idx_get_c
  end interface

  idx = PDM_block_to_part_gnum_idx_get_c (btp,  &
                                          gNum)

end subroutine PDM_block_to_part_gnum_idx_get






subroutine PDM_block_to_part_exch_finalize (btp,           &
                                            t_stride,      &
                                            block_stride,  &
                                            data_type,     &
                                            c_part_stride, &
                                            c_part_data,   &
                                            part_stride,   &
                                            part_data,     &
                                            s_data)
  use iso_c_binding
  implicit none

  type(c_ptr), value                :: btp
  integer, intent(in)               :: t_stride
  integer(pdm_l_num_s), pointer     :: block_stride(:)
  integer, intent(in)               :: data_type
  type(c_ptr), value                :: c_part_stride
  type(c_ptr), value                :: c_part_data
  type(PDM_pointer_array_t), pointer:: part_stride
  type(PDM_pointer_array_t), pointer:: part_data
  integer, intent(in), optional     :: s_data

  integer                           :: n_part, n_elt
  integer(pdm_l_num_s), pointer     :: stride(:)
  integer                           :: s_part_data
  integer                           :: i, j

  integer, allocatable              :: length_stride(:)
  integer, allocatable              :: length_data(:)

  n_part = PDM_block_to_part_n_part_get (btp)
  stride => null()

  allocate(length_data(n_part))

  if (t_stride .eq. PDM_STRIDE_VAR_INTERLACED) then

    allocate(length_stride(n_part))

    do i = 1, n_part

      length_stride(i) = PDM_block_to_part_n_elt_get(btp, &
                                              i-1)
    end do

    call  PDM_pointer_array_create (part_stride,     &
                                    n_part, &
                                    PDM_TYPE_INT, &
                                    c_part_stride, &
                                    length_stride, &
                                    PDM_OWNERSHIP_KEEP)

    do i = 1, n_part

      call PDM_pointer_array_part_get (part_stride, &
                                       i-1,         &
                                       stride)

      s_part_data = 0
      do j = 1, length_stride(i)
        s_part_data = s_part_data + stride(j)
      end do

      length_data(i) = s_part_data
    
    enddo  

    if (present(s_data)) then

      call  PDM_pointer_array_create (part_data,     &
                                      n_part, &
                                      data_type, &
                                      c_part_data, &
                                      length_data,      &
                                      PDM_OWNERSHIP_KEEP, &
                                      s_data)

    else 

      if (data_type .eq. PDM_TYPE_CPTR) then
        print *, "Error PDM_block_to_part_exch : s-data is mandatory for PDM_TYPE_CPTR type"
        call exit
      endif

      call  PDM_pointer_array_create (part_data,     &
                                      n_part, &
                                      data_type, &
                                      c_part_data, &
                                      length_data, &
                                      PDM_OWNERSHIP_KEEP)

    endif  

    deallocate(length_stride)

  else

    part_stride => null() 

    do i = 1, n_part

      n_elt = PDM_block_to_part_n_elt_get(btp, &
                                          i-1)

      if (size(block_stride) .gt. 0) then
        s_part_data = block_stride(1) * n_elt
      end if

      length_data(i) = s_part_data

    end do

    if (present(s_data)) then

      call  PDM_pointer_array_create (part_data,     &
                                      n_part, &
                                      data_type, &
                                      c_part_data, &
                                      length_data,      &
                                      PDM_OWNERSHIP_KEEP, &
                                      s_data)

    else 

      if (data_type .eq. PDM_TYPE_CPTR) then
        print *, "Error PDM_block_to_part_exch : s-data is mandatory for PDM_TYPE_CPTR type"
        call exit
      endif

      call  PDM_pointer_array_create (part_data,     &
                                      n_part, &
                                      data_type, &
                                      c_part_data, &
                                      length_data, &
                                      PDM_OWNERSHIP_KEEP);

    endif  

  endif 

  deallocate(length_data)

end subroutine PDM_block_to_part_exch_finalize


subroutine PDM_block_to_part_exch_in_place_finalize (btp,           &
                                                     t_stride,      &
                                                     block_stride,  &
                                                     part_stride,   &
                                                     part_data)
  use iso_c_binding
  implicit none

  type(c_ptr), value                :: btp
  integer, intent(in)               :: t_stride
  integer(pdm_l_num_s), pointer     :: block_stride(:)
  type(PDM_pointer_array_t), pointer :: part_stride
  type(PDM_pointer_array_t), pointer :: part_data

  integer                           :: n_part, n_elt
  integer(pdm_l_num_s), pointer     :: stride(:)
  integer                           :: s_part_data
  integer                           :: i, j

  n_part = PDM_block_to_part_n_part_get (btp)
  stride => null()

  do i = 1, n_part

    n_elt = PDM_block_to_part_n_elt_get(btp, &
                                        i-1)

    if (t_stride .eq. PDM_STRIDE_VAR_INTERLACED) then

      call PDM_pointer_array_part_get (part_stride, &
                                       i-1,         &
                                       stride)

      s_part_data = 0
      do j = 1, n_elt
        s_part_data = s_part_data + stride(j)
      end do

    else

      if (size(block_stride) .gt. 0) then
        s_part_data = block_stride(1) * n_elt
      end if

    end if

    call PDM_pointer_array_part_length_update (part_data, &
                                               i-1, &
                                               s_part_data)

  end do


end subroutine PDM_block_to_part_exch_in_place_finalize

end module pdm_block_to_part
