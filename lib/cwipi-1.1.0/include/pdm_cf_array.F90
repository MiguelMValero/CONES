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

module pdm_cf_array

  use pdm
  use pdm_fortran
  use iso_c_binding

  implicit none

  type PDM_cf_array_t

    integer               :: type      = -1
    integer               :: length    = -1
    type(c_ptr)           :: cptr   
    integer               :: ownership = PDM_OWNERSHIP_KEEP

  end type PDM_cf_array_t

   interface PDM_cf_array_part_get
     module procedure PDM_cf_array_part_get_int_
     module procedure PDM_cf_array_part_get_int_2_
     module procedure PDM_cf_array_part_get_int_3_
#ifdef PDM_LONG_G_NUM
     module procedure PDM_cf_array_part_get_g_num_
     module procedure PDM_cf_array_part_get_g_num_2_
     module procedure PDM_cf_array_part_get_g_num_3_
#endif
     module procedure PDM_cf_array_part_get_double_
     module procedure PDM_cf_array_part_get_double_2_
     module procedure PDM_cf_array_part_get_double_3_
     module procedure PDM_cf_array_part_get_complex8_
     module procedure PDM_cf_array_part_get_complex4_
     module procedure PDM_cf_array_part_get_real4_
   end interface

  interface PDM_cf_array_create
    module procedure PDM_cf_array_create_
  end interface


  private :: &
              PDM_cf_array_part_get_int_, &
              PDM_cf_array_part_get_int_2_, &
              PDM_cf_array_part_get_int_3_, &
#ifdef PDM_LONG_G_NUM
              PDM_cf_array_part_get_g_num_, &
              PDM_cf_array_part_get_g_num_2_, &
              PDM_cf_array_part_get_g_num_3_, &
#endif
              PDM_cf_array_part_get_double_, &
              PDM_cf_array_part_get_double_2_, &
              PDM_cf_array_part_get_double_3_, &
              PDM_cf_array_part_get_complex8_, &
              PDM_cf_array_part_get_complex4_, &
              PDM_cf_array_part_get_real4_, &
              PDM_cf_array_create_ 
  contains


  !>
  !! \brief Initialize a \ref PDM_cf_array_t object
  !!
  !! \param [out]  pa      \ref PDM_cf_array_t object
  !! \param [in]   type    Data type of pointers
  !! \param [in]   length  Length of array
  !!

  subroutine PDM_cf_array_create_ (pa,        &
                               type,      &
                               length,    &
                               cptr,      &
                               ownership)
    use iso_c_binding
    implicit none

    type(PDM_cf_array_t), pointer :: pa
    integer, intent(in)        :: type
    integer, intent(in)        :: length
    type(c_ptr), intent(in)    :: cptr
    integer, intent(in)        :: ownership

    if (associated(pa)) then
      print*, "Error PDM_cf_array_create : pa is already associated ! "
      call exit
    endif

    allocate(pa)

    pa%type      = type
    pa%length    = length
    pa%cptr      = cptr
    pa%ownership = ownership

  end subroutine PDM_cf_array_create_


  !>
  !! \brief Free a \ref PDM_cf_array_t object
  !!
  !! \param [in, out]  pa      \ref PDM_cf_array_t object
  !!

  subroutine PDM_cf_array_free_ (pa)
    use iso_c_binding
    implicit none

    type(PDM_cf_array_t), pointer  :: pa

    if (.not. associated(pa)) then
      print*, "Error PDM_cf_array_free : 'pa' pointer is not associated "
      call exit
    endif

    if (pa%ownership .eq. PDM_OWNERSHIP_KEEP) then
      call pdm_fortran_free_c(pa%cptr)
    endif

    deallocate(pa)
    pa => null()

  end subroutine PDM_cf_array_free_


  !>
  !! \brief Get an array
  !!
  !! Maps a Fortran pointer onto a C pointer
  !!
  !! \param [in]       pa         Array of \ref PDM_cf_array_t
  !! \param [in]       stride1    Dimension 1 of multi-dimension array: Optional
  !! \param [in]       stride2    Dimension 2 of multi-dimension array: Optional
  !! \param [in]       stride3    Dimension 3 of multi-dimension array: Optional
  !! \param [in, out]  pointer_f  Pointer to an integer array
  !!

  subroutine PDM_cf_array_part_get_int_ (pa,        &
                                      pointer_f)
    use iso_c_binding
    implicit none

    type(PDM_cf_array_t),     pointer :: pa
    integer(pdm_l_num_s),  pointer :: pointer_f(:)

    if (.not. associated(pa)) then
      print*, "Error PDM_cf_array_part_get : 'pa' pointer is not associated "
      call exit
    endif

    if (pa%type .ne. PDM_TYPE_INT) then
      print *, "PDM_cf_array_part_get : wrong type for pointer_f"
      stop
    end if

    call c_f_pointer(pa%cptr,      &
                     pointer_f,    &
                     [pa%length])

  end subroutine PDM_cf_array_part_get_int_


  subroutine PDM_cf_array_part_get_int_2_ (pa,        &
                                        stride1,   & 
                                        pointer_f)
    use iso_c_binding
    implicit none

    type(PDM_cf_array_t),     pointer :: pa
    integer, intent(in)            :: stride1
    integer(pdm_l_num_s),  pointer :: pointer_f(:,:)

    if (.not. associated(pa)) then
      print*, "Error PDM_cf_array_part_get : 'pa' pointer is not associated "
      call exit
    endif

    if (pa%type .ne. PDM_TYPE_INT) then
      print *, "PDM_cf_array_part_get : wrong type for pointer_f"
      stop
    end if

    call c_f_pointer(pa%cptr,      &
                     pointer_f,    &
                     [stride1, pa%length/stride1])

  end subroutine PDM_cf_array_part_get_int_2_


  subroutine PDM_cf_array_part_get_int_3_ (pa,        &
                                        stride1,   & 
                                        stride2,   & 
                                        pointer_f)
    use iso_c_binding
    implicit none

    type(PDM_cf_array_t),     pointer :: pa
    integer, intent(in)            :: stride1
    integer, intent(in)            :: stride2
    integer(pdm_l_num_s),  pointer :: pointer_f(:,:,:)

    if (.not. associated(pa)) then
      print*, "Error PDM_cf_array_part_get : 'pa' pointer is not associated "
      call exit
    endif

    if (pa%type .ne. PDM_TYPE_INT) then
      print *, "PDM_cf_array_part_get : wrong type for pointer_f"
      stop
    end if

    call c_f_pointer(pa%cptr,      &
                     pointer_f,    &
                     [stride1, stride2, pa%length/(stride1+stride2)])

  end subroutine PDM_cf_array_part_get_int_3_

#ifdef PDM_LONG_G_NUM

  subroutine PDM_cf_array_part_get_g_num_ (pa,        &
                                      pointer_f)
    use iso_c_binding
    implicit none

    type(PDM_cf_array_t),     pointer :: pa
    integer(pdm_g_num_s),  pointer :: pointer_f(:)

    if (.not. associated(pa)) then
      print*, "Error PDM_cf_array_part_get : 'pa' pointer is not associated "
      call exit
    endif

    if (pa%type .ne. PDM_TYPE_G_NUM) then
      print *, "PDM_cf_array_part_get : wrong type for pointer_f"
      stop
    end if

    call c_f_pointer(pa%cptr,      &
                     pointer_f,    &
                     [pa%length])

  end subroutine PDM_cf_array_part_get_g_num_

  subroutine PDM_cf_array_part_get_g_num_2_ (pa,        &
                                          stride1,   & 
                                          pointer_f)
    use iso_c_binding
    implicit none

    type(PDM_cf_array_t),     pointer :: pa
    integer, intent(in)            :: stride1
    integer(pdm_g_num_s),  pointer :: pointer_f(:,:)

    if (.not. associated(pa)) then
      print*, "Error PDM_cf_array_part_get : 'pa' pointer is not associated "
      call exit
    endif

    if (pa%type .ne. PDM_TYPE_G_NUM) then
      print *, "PDM_cf_array_part_get : wrong type for pointer_f"
      stop
    end if

    call c_f_pointer(pa%cptr,      &
                     pointer_f,    &
                     [stride1, pa%length/stride1])

  end subroutine PDM_cf_array_part_get_g_num_2_

  subroutine PDM_cf_array_part_get_g_num_3_ (pa,        &
                                          stride1,   & 
                                          stride2,   & 
                                          pointer_f)
    use iso_c_binding
    implicit none

    type(PDM_cf_array_t),     pointer :: pa
    integer, intent(in)            :: stride1
    integer, intent(in)            :: stride2
    integer(pdm_g_num_s),  pointer :: pointer_f(:,:,:)

    if (.not. associated(pa)) then
      print*, "Error PDM_cf_array_part_get : 'pa' pointer is not associated "
      call exit
    endif

    if (pa%type .ne. PDM_TYPE_G_NUM) then
      print *, "PDM_cf_array_part_get : wrong type for pointer_f"
      stop
    end if

    call c_f_pointer(pa%cptr,      &
                     pointer_f,    &
                     [stride1, stride2, pa%length/(stride1+stride2)])

  end subroutine PDM_cf_array_part_get_g_num_3_

#endif

  subroutine PDM_cf_array_part_get_double_ (pa,        &
                                         pointer_f)
    use iso_c_binding
    implicit none

    type(PDM_cf_array_t),     pointer :: pa
    double precision,  pointer :: pointer_f(:)

    if (.not. associated(pa)) then
      print*, "Error PDM_cf_array_part_get : 'pa' pointer is not associated "
      call exit
    endif

    if (pa%type .ne. PDM_TYPE_DOUBLE) then
      print *, "PDM_cf_array_part_get : wrong type for pointer_f"
      stop
    end if

    call c_f_pointer(pa%cptr,      &
                     pointer_f,    &
                     [pa%length])

  end subroutine PDM_cf_array_part_get_double_

  subroutine PDM_cf_array_part_get_double_2_ (pa,        &
                                           stride1,   & 
                                           pointer_f)
    use iso_c_binding
    implicit none

    type(PDM_cf_array_t),     pointer :: pa
    integer, intent(in)            :: stride1
    double precision,  pointer :: pointer_f(:,:)

    if (.not. associated(pa)) then
      print*, "Error PDM_cf_array_part_get : 'pa' pointer is not associated "
      call exit
    endif

    if (pa%type .ne. PDM_TYPE_DOUBLE) then
      print *, "PDM_cf_array_part_get : wrong type for pointer_f"
      stop
    end if

    call c_f_pointer(pa%cptr,      &
                     pointer_f,    &
                     [stride1, pa%length/stride1])

  end subroutine PDM_cf_array_part_get_double_2_

  subroutine PDM_cf_array_part_get_double_3_ (pa,        &
                                           stride1,   & 
                                           stride2,   & 
                                           pointer_f)
    use iso_c_binding
    implicit none

    type(PDM_cf_array_t),     pointer :: pa
    integer, intent(in)            :: stride1
    integer, intent(in)            :: stride2
    double precision, pointer :: pointer_f(:,:,:)

    if (.not. associated(pa)) then
      print*, "Error PDM_cf_array_part_get : 'pa' pointer is not associated "
      call exit
    endif

    if (pa%type .ne. PDM_TYPE_DOUBLE) then
      print *, "PDM_cf_array_part_get : wrong type for pointer_f"
      stop
    end if

    call c_f_pointer(pa%cptr,      &
                     pointer_f,    &
                     [stride1, stride2, pa%length/(stride1+stride2)])

  end subroutine PDM_cf_array_part_get_double_3_


  subroutine PDM_cf_array_part_get_complex8_ (pa,        &
                                           pointer_f)
    use iso_c_binding
    implicit none

    type(PDM_cf_array_t),     pointer :: pa
    complex (kind = 8),    pointer :: pointer_f(:)

    if (.not. associated(pa)) then
      print*, "Error PDM_cf_array_part_get : 'pa' pointer is not associated "
      call exit
    endif

    if (pa%type .ne. PDM_TYPE_COMPLEX8) then
      print *, "PDM_cf_array_part_get : wrong type for pointer_f"
      stop
    end if

    call c_f_pointer(pa%cptr,      &
                     pointer_f,    &
                     [pa%length])

  end subroutine PDM_cf_array_part_get_complex8_


  subroutine PDM_cf_array_part_get_complex8_2_ (pa,        &
                                             stride1,   & 
                                             pointer_f)
    use iso_c_binding
    implicit none

    type(PDM_cf_array_t),     pointer :: pa
    integer, intent(in)            :: stride1
    complex (kind = 8),    pointer :: pointer_f(:,:)

    if (.not. associated(pa)) then
      print *, "PDM_cf_array_part_get : wrong type for pointer_f"
      call exit
    endif

    if (pa%type .ne. PDM_TYPE_COMPLEX8) then
      print *, "PDM_cf_array_part_get : wrong type"
      stop
    end if

    call c_f_pointer(pa%cptr,      &
                     pointer_f,    &
                     [stride1, pa%length/stride1])

  end subroutine PDM_cf_array_part_get_complex8_2_

  subroutine PDM_cf_array_part_get_complex8_3_ (pa,        &
                                             stride1,   & 
                                             stride2,   & 
                                             pointer_f)
    use iso_c_binding
    implicit none

    type(PDM_cf_array_t),     pointer :: pa
    integer, intent(in)            :: stride1
    integer, intent(in)            :: stride2
    complex (kind = 8),    pointer :: pointer_f(:,:,:)

    if (.not. associated(pa)) then
      print*, "Error PDM_cf_array_part_get : 'pa' pointer is not associated "
      call exit
    endif

    if (pa%type .ne. PDM_TYPE_COMPLEX8) then
      print *, "PDM_cf_array_part_get : wrong type for pointer_f"
      stop
    end if

    call c_f_pointer(pa%cptr,      &
                     pointer_f,    &
                     [stride1, stride2, pa%length/(stride1+stride2)])

  end subroutine PDM_cf_array_part_get_complex8_3_


  subroutine PDM_cf_array_part_get_complex4_ (pa,        &
                                           pointer_f)
    use iso_c_binding
    implicit none

    type(PDM_cf_array_t),     pointer :: pa
    complex (kind = 4),    pointer :: pointer_f(:)

    if (.not. associated(pa)) then
      print*, "Error PDM_cf_array_part_get : 'pa' pointer is not associated "
      call exit
    endif

    if (pa%type .ne. PDM_TYPE_COMPLEX4) then
      print *, "PDM_cf_array_part_get : wrong type for pointer_f"
      stop
    end if

    call c_f_pointer(pa%cptr,      &
                     pointer_f,    &
                     [pa%length])

  end subroutine PDM_cf_array_part_get_complex4_


  subroutine PDM_cf_array_part_get_complex4_2_ (pa,        &
                                             stride1,   & 
                                             pointer_f)
    use iso_c_binding
    implicit none

    type(PDM_cf_array_t),     pointer :: pa
    integer, intent(in)            :: stride1
    complex (kind = 4),    pointer :: pointer_f(:,:)

    if (.not. associated(pa)) then
      print *, "PDM_cf_array_part_get : wrong type for pointer_f"
      call exit
    endif

    if (pa%type .ne. PDM_TYPE_COMPLEX4) then
      print *, "PDM_cf_array_part_get : wrong type"
      stop
    end if

    call c_f_pointer(pa%cptr,      &
                     pointer_f,    &
                     [stride1, pa%length/stride1])

  end subroutine PDM_cf_array_part_get_complex4_2_

  subroutine PDM_cf_array_part_get_complex4_3_ (pa,        &
                                             stride1,   & 
                                             stride2,   & 
                                             pointer_f)
    use iso_c_binding
    implicit none

    type(PDM_cf_array_t),     pointer :: pa
    integer, intent(in)            :: stride1
    integer, intent(in)            :: stride2
    complex (kind = 4),    pointer :: pointer_f(:,:,:)

    if (.not. associated(pa)) then
      print*, "Error PDM_cf_array_part_get : 'pa' pointer is not associated "
      call exit
    endif

    if (pa%type .ne. PDM_TYPE_COMPLEX4) then
      print *, "PDM_cf_array_part_get : wrong type for pointer_f"
      stop
    end if

    call c_f_pointer(pa%cptr,      &
                     pointer_f,    &
                     [stride1, stride2, pa%length/(stride1+stride2)])

  end subroutine PDM_cf_array_part_get_complex4_3_


  subroutine PDM_cf_array_part_get_real4_ (pa,        &
                                        pointer_f)
    use iso_c_binding
    implicit none

    type(PDM_cf_array_t),     pointer :: pa
    real (kind = 4),       pointer :: pointer_f(:)

    if (.not. associated(pa)) then
      print*, "Error PDM_cf_array_part_get : 'pa' pointer is not associated "
      call exit
    endif

    if (pa%type .ne. PDM_TYPE_REAL4) then
      print *, "PDM_cf_array_part_get : wrong type for pointer_f"
      stop
    end if

    call c_f_pointer(pa%cptr,      &
                     pointer_f,    &
                     [pa%length])

  end subroutine PDM_cf_array_part_get_real4_


  subroutine PDM_cf_array_part_get_real4_2_ (pa,        &
                                          stride1,   & 
                                          pointer_f)
    use iso_c_binding
    implicit none

    type(PDM_cf_array_t),     pointer :: pa
    integer, intent(in)            :: stride1
    real (kind = 4),       pointer :: pointer_f(:,:)

    if (.not. associated(pa)) then
      print *, "PDM_cf_array_part_get : wrong type for pointer_f"
      call exit
    endif

    if (pa%type .ne. PDM_TYPE_REAL4) then
      print *, "PDM_cf_array_part_get : wrong type"
      stop
    end if

    call c_f_pointer(pa%cptr,      &
                     pointer_f,    &
                     [stride1, pa%length/stride1])

  end subroutine PDM_cf_array_part_get_real4_2_

  subroutine PDM_cf_array_part_get_real4_3_ (pa,        &
                                          stride1,   & 
                                          stride2,   & 
                                          pointer_f)
    use iso_c_binding
    implicit none

    type(PDM_cf_array_t),     pointer :: pa
    integer, intent(in)            :: stride1
    integer, intent(in)            :: stride2
    real (kind = 4),       pointer :: pointer_f(:,:,:)

    if (.not. associated(pa)) then
      print*, "Error PDM_cf_array_part_get : 'pa' pointer is not associated "
      call exit
    endif

    if (pa%type .ne. PDM_TYPE_COMPLEX4) then
      print *, "PDM_cf_array_part_get : wrong type for pointer_f"
      stop
    end if

    call c_f_pointer(pa%cptr,      &
                     pointer_f,    &
                     [stride1, stride2, pa%length/(stride1+stride2)])

  end subroutine PDM_cf_array_part_get_real4_3_

end module pdm_cf_array
