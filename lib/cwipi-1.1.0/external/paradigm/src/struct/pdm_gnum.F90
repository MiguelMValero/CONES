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

module pdm_gnum

  use pdm

  implicit none

  interface PDM_gnum_create ; module procedure &
  pdm_gnum_create_
  end interface

  interface PDM_gnum_set_from_coords ; module procedure &
  pdm_gnum_set_from_coords_
  end interface

  interface PDM_gnum_set_from_parents ; module procedure &
  pdm_gnum_set_from_parents_
  end interface

  interface PDM_gnum_get ; module procedure &
  pdm_gnum_get_
  end interface

  private :: pdm_gnum_create_
  ! private :: pdm_gnum_set_from_coords_
  private :: pdm_gnum_set_from_parents_
  private :: pdm_gnum_get_

  interface

  !>
  !!
  !! \brief Build a global numbering structure
  !!
  !! \param [in]   dim          Spatial dimension
  !! \param [in]   n_part       Number of local partitions
  !! \param [in]   merge        Merge double points or not
  !! \param [in]   tolerance    Geometric tolerance (if merge double points is activated)
  !! \param [in]   comm         PDM_MPI communicator
  !!
  !! \return     Pointer to \ref PDM_gen_gnum object
  !!

  function PDM_gnum_create_cf (dim,        &
                               n_part,     &
                               merge,      &
                               tolerance,  &
                               comm,       &
                               owner)      &
    result (gen_gnum) &

    bind (c, name = 'PDM_gnum_create')

    use iso_c_binding
    implicit none

    integer(c_int), value :: dim
    integer(c_int), value :: n_part
    integer(c_int), value :: merge
    real(c_double), value :: tolerance
    integer(c_int), value :: comm
    integer(c_int), value :: owner

    type (c_ptr) :: gen_gnum

  end function PDM_gnum_create_cf


  !>
  !!
  !! \brief Set from coordinates
  !!
  !! \param [in]   gen_gnum     Pointer to \ref PDM_gen_gnum object
  !! \param [in]   i_part       Current partition
  !! \param [in]   n_elts       Number of elements
  !! \param [in]   coords       Coordinates (size = 3 * \ref n_elts)
  !! \param [in]   char_length  Characteristic length (or NULL)
  !!                            (used if merge double points is activated)
  !!
  !!

  subroutine PDM_gnum_set_from_coords_cf (gen_gnum,    &
                                          i_part,      &
                                          n_elts,      &
                                          coords,      &
                                          char_length) &
    bind (c, name = 'PDM_gnum_set_from_coords')

    use iso_c_binding
    implicit none

    type(c_ptr),    value :: gen_gnum

    integer(c_int), value :: i_part
    integer(c_int), value :: n_elts
    type(c_ptr),    value :: coords
    type(c_ptr),    value :: char_length

  end subroutine PDM_gnum_set_from_coords_cf


  !>
  !!
  !! \brief Set Parent global numbering
  !!
  !! \param [in]   gen_gnum     Pointer to \ref PDM_gen_gnum object
  !! \param [in]   i_part       Current partition
  !! \param [in]   n_elts       Number of elements
  !! \param [in]   parent_gnum  Parent global numbering (size = \ref n_elts)
  !!
  !!

  subroutine PDM_gnum_set_from_parents_cf (gen_gnum,    &
                                           i_part,      &
                                           n_elts,      &
                                           parent_gnum) &
    bind (c, name = 'PDM_gnum_set_from_parents')

    use iso_c_binding
    implicit none

    type(c_ptr),    value :: gen_gnum

    integer(c_int), value :: i_part
    integer(c_int), value :: n_elts
    type(c_ptr),    value :: parent_gnum

  end subroutine PDM_gnum_set_from_parents_cf


  !>
  !!
  !! \brief Compute
  !!
  !! \param [in]   gen_gnum         Pointer to \ref PDM_gen_gnum object
  !!
  !!

  subroutine PDM_gnum_compute (gen_gnum) &
    bind (c, name = 'PDM_gnum_compute')

    use iso_c_binding
    implicit none

    type(c_ptr), value :: gen_gnum

  end subroutine PDM_gnum_compute


  !>
  !!
  !! \brief Get global ids for a given partition
  !!
  !! \param [in]   gen_gnum     Pointer to \ref PDM_gen_gnum object
  !! \param [in]   i_part       Current partition
  !!
  !! \return     Array of global ids
  !!
  !!

  function PDM_gnum_get_cf (gen_gnum, &
                            i_part)   &
    result (g_nums) &
    bind (c, name = 'PDM_gnum_get')

    use iso_c_binding
    implicit none

    type(c_ptr),    value :: gen_gnum
    integer(c_int), value :: i_part
    type(c_ptr)           :: g_nums

  end function PDM_gnum_get_cf


  !>
  !!
  !! \brief Free
  !!
  !! \param [in]   gen_gnum         Pointer to \ref PDM_gen_gnum object
  !!
  !!

  subroutine PDM_gnum_free (gen_gnum)     &
    bind (c, name = 'PDM_gnum_free')
    use iso_c_binding
    implicit none

    type (c_ptr), value :: gen_gnum

    end subroutine PDM_gnum_free


  !>
  !!
  !! \brief Get number of elements in a partition
  !!
  !! \param [in]   gen_gnum     Pointer to \ref PDM_gen_gnum object
  !! \param [in]   i_part       Current partition
  !! \return     Number of elements
  !!
  !!

  function PDM_gnum_n_elt_get (gen_gnum, &
                               i_part)   &
    result (n_elts)                      &
    bind (c, name='PDM_gnum_n_elt_get')

    use iso_c_binding
    implicit none

    type (c_ptr),   value :: gen_gnum
    integer(c_int), value :: i_part
    integer(c_int)        :: n_elts

  end function PDM_gnum_n_elt_get


  end interface


contains


  subroutine PDM_gnum_create_ (gen_gnum,  &
                               dim,       &
                               n_part,    &
                               merge,     &
                               tolerance, &
                               f_comm,    &
                               owner)
  ! Build a global numbering structure
  use iso_c_binding
  implicit none

  type (c_ptr)                 :: gen_gnum  ! C pointer to PDM_gen_gnum_t object
  integer,          intent(in) :: dim       ! Spatial dimension
  integer,          intent(in) :: n_part    ! Number of local partitions
  integer,          intent(in) :: merge     ! Merge coincident points or not
  double precision, intent(in) :: tolerance ! Geometric tolerance (used only if ``merge`` is enabled)
  integer,          intent(in) :: f_comm    ! Fortran MPI communicator
  integer,          intent(in) :: owner     ! Ownership

  integer(c_int)               :: c_dim
  integer(c_int)               :: c_n_part
  integer(c_int)               :: c_merge
  real(c_double)               :: c_tolerance
  integer(c_int)               :: c_comm
  integer(c_int)               :: c_owner


  c_comm = PDM_MPI_Comm_f2c(f_comm)

  c_dim       = dim
  c_n_part    = n_part
  c_merge     = merge
  c_tolerance = tolerance
  c_owner     = owner

  gen_gnum = PDM_gnum_create_cf (c_dim,       &
                                 c_n_part,    &
                                 c_merge,     &
                                 c_tolerance, &
                                 c_comm,      &
                                 c_owner)

  end subroutine PDM_gnum_create_



  subroutine PDM_gnum_set_from_coords_(gen_gnum,    &
                                       i_part,      &
                                       n_elts,      &
                                       coords,      &
                                       char_length)
    ! Set from coordinates
    use iso_c_binding
    implicit none

    type(c_ptr), value        :: gen_gnum       ! C pointer to PDM_gen_gnum_t object
    integer, intent(in)       :: i_part         ! Current partition
    integer, intent(in)       :: n_elts         ! Number of elements
    double precision, pointer :: coords(:,:)    ! Coordinates (size = 3 * ``n_elts``)
    double precision, pointer :: char_length(:) ! Characteristic length (or *null()*) (used only if ``merge`` was enabled in PDM_gnum_create)

    integer(c_int)            :: c_i_part
    integer(c_int)            :: c_n_elts
    type(c_ptr)               :: c_coords
    type(c_ptr)               :: c_char_length

    c_i_part = i_part
    c_n_elts = n_elts

    c_coords      = c_loc(coords)
    c_char_length = c_loc(char_length)


    call PDM_gnum_set_from_coords_cf (gen_gnum,      &
                                      c_i_part,      &
                                      c_n_elts,      &
                                      c_coords,      &
                                      c_char_length)

  end subroutine PDM_gnum_set_from_coords_



  subroutine PDM_gnum_set_from_parents_ (gen_gnum,    &
                                         i_part,      &
                                         n_elts,      &
                                         parent_gnum)
    ! Set parent global numbering
    use iso_c_binding
    implicit none

    type(c_ptr), value                 :: gen_gnum       ! C pointer to PDM_gen_gnum_t object
    integer, intent(in)                :: i_part         ! Current partition
    integer, intent(in)                :: n_elts         ! Number of elements
    integer(kind=pdm_g_num_s), pointer :: parent_gnum(:) ! Parent global numbering (size = ``n_elts``)

    integer(c_int)                     :: c_i_part
    integer(c_int)                     :: c_n_elts
    type(c_ptr)                        :: c_parent_gnum

    c_i_part = i_part
    c_n_elts = n_elts

    c_parent_gnum = c_loc(parent_gnum)


    call PDM_gnum_set_from_parents_cf (gen_gnum,      &
                                       c_i_part,      &
                                       c_n_elts,      &
                                       c_parent_gnum)

  end subroutine PDM_gnum_set_from_parents_



  subroutine PDM_gnum_get_ (gen_gnum, &
                            i_part,   &
                            g_nums)
    ! Get global ids for a given partition
    use iso_c_binding
    implicit none

    type(c_ptr), value                 :: gen_gnum  ! C pointer to PDM_gen_gnum_t object
    integer, intent(in)                :: i_part    ! Current partition
    integer(kind=pdm_g_num_s), pointer :: g_nums(:) ! Array of global ids

    integer(c_int)                     :: c_i_part
    type(c_ptr)                        :: c_g_nums = C_NULL_PTR
    integer                            :: n_elts

    c_i_part = i_part

    c_g_nums = PDM_gnum_get_cf(gen_gnum, &
                               c_i_part)

    n_elts = PDM_gnum_n_elt_get(gen_gnum, &
                                c_i_part)


    call c_f_pointer(c_g_nums, &
                     g_nums,   &
                     [n_elts])

  end subroutine PDM_gnum_get_

end module pdm_gnum
