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

module pdm_global_mean

  use pdm

  implicit none

  interface

    !>
    !!
    !! \brief Compute the global average field
    !!
    !! \param [in]   gmean         Pointer to \ref PDM_global_mean object
    !!
    !!

    subroutine pdm_global_mean_field_compute (gmean) &
         bind (c, name = 'PDM_global_mean_field_compute')

      use iso_c_binding

      implicit none

      type (c_ptr), value :: gmean

    end subroutine pdm_global_mean_field_compute


    !>
    !!
    !! \brief Free a global point mean structure
    !!
    !! \param [in]   gmean           Pointer to \ref PDM_global_mean object
    !!
    !!

    subroutine pdm_global_mean_free (gmean) &
         bind (c, name = 'PDM_global_mean_free')

      use iso_c_binding

      implicit none

      type (c_ptr), value :: gmean

    end subroutine PDM_global_mean_free

  end interface



  contains

  !>
  !!
  !! \brief Create a structure that compute a global mean
  !!
  !! \param [out]  gmean     Pointer to \ref PDM_global_mean object
  !! \param [in]   n_part    Number of local partitions
  !! \param [in]   comm      PDM_MPI communicator
  !!

  subroutine PDM_global_mean_create (gmean,  &
                                     n_part, &
                                     f_comm)

  use iso_c_binding
  implicit none

  integer, intent(in) :: n_part
  integer, intent(in) :: f_comm
  type (c_ptr)        :: gmean

  integer(c_int)      :: c_n_part
  integer(c_int)      :: c_comm

  interface
    function pdm_global_mean_create_cf (n_part, fComm) &
    result (gmean)                                     &
    bind (c, name = 'PDM_global_mean_create')
      use iso_c_binding
      implicit none

      integer(c_int), value :: n_part
      integer(c_int), value :: fComm
      type (c_ptr)          :: gmean

    end function pdm_global_mean_create_cf
  end interface

  c_comm   = PDM_MPI_Comm_f2c(f_comm)
  c_n_part = n_part

  gmean = pdm_global_mean_create_cf (c_n_part, &
                                     c_comm)

  end subroutine PDM_global_mean_create


  !>
  !!
  !! \brief Set absolute number
  !!
  !! \param [in]   gmean        Pointer to \ref PDM_global_mean object
  !! \param [in]   i_part       Current partition
  !! \param [in]   n_point      Number of points in the partition
  !! \param [in]   numabs       Absolute number of points
  !!
  !!
  subroutine PDM_global_mean_set (gmean,   &
                                  i_part,  &
                                  n_point, &
                                  numabs)
    use iso_c_binding
    implicit none

    type(c_ptr), value            :: gmean
    integer, intent(in)           :: i_part
    integer, intent(in)           :: n_point
    integer(pdm_g_num_s), pointer :: numabs(:)

    interface
      subroutine pdm_global_mean_set_c (gmean,   &
                                        i_part,  &
                                        n_point, &
                                        numabs)  &
      bind (c, name = 'PDM_global_mean_set')
        use iso_c_binding
        implicit none

        type (c_ptr),   value :: gmean
        integer(c_int), value :: i_part
        integer(c_int), value :: n_point
        type (c_ptr),   value :: numabs

      end subroutine pdm_global_mean_set_c
    end interface

    call pdm_global_mean_set_c (gmean,         &
                                i_part,        &
                                n_point,       &
                                c_loc(numabs))

  end subroutine PDM_global_mean_set


  !>
  !!
  !! \brief Set local field and it associated weight
  !!
  !! \param [in]   gmean                 Pointer to \ref PDM_global_mean object
  !! \param [in]   i_part                Current partition
  !! \param [in]   stride                Stride of the field
  !! \param [in]   local_field           Local value of field
  !! \param [in]   local_weight          Local weight used to compute the mean value
  !! \param [in]   global_mean_field_ptr Pointer where global mean field
  !!                                     will be stored after computing
  !!

  subroutine pdm_global_mean_field_set (gmean,                 &
                                        i_part,                &
                                        stride,                &
                                        local_field,           &
                                        local_weight,          &
                                        global_mean_field_ptr)
    use iso_c_binding
    implicit none

    type(c_ptr), value        :: gmean
    integer, intent(in)       :: i_part
    integer, intent(in)       :: stride
    double precision, pointer :: local_field(:)
    double precision, pointer :: local_weight(:)
    double precision, pointer :: global_mean_field_ptr(:)

    interface
      subroutine pdm_global_mean_field_set_c (gmean,                 &
                                              i_part,                &
                                              stride,                &
                                              local_field,           &
                                              local_weight,          &
                                              global_mean_field_ptr) &
      bind (c, name = 'PDM_global_mean_field_set')
        use iso_c_binding
        implicit none

        type (c_ptr),   value :: gmean
        integer(c_int), value :: i_part
        integer(c_int), value :: stride
        type (c_ptr),   value :: local_field
        type (c_ptr),   value :: local_weight
        type (c_ptr),   value :: global_mean_field_ptr

      end subroutine pdm_global_mean_field_set_c
    end interface

    call pdm_global_mean_field_set_c (gmean,                        &
                                      i_part,                       &
                                      stride,                       &
                                      c_loc(local_field),           &
                                      c_loc(local_weight),          &
                                      c_loc(global_mean_field_ptr))

  end subroutine pdm_global_mean_field_set

end module pdm_global_mean
