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

module pdm_gnum_location

  use pdm

  implicit none


  interface

    !>
    !!
    !! \brief Compute the location (processus, partittion, local number in the partition)
    !!
    !! \param [in]   gloc         Pointer to \ref PDM_gnum_locaion object
    !!
    !!

    subroutine pdm_gnum_location_compute (gloc) &
    bind (c, name = 'PDM_gnum_location_compute')
      use iso_c_binding
      implicit none

      type (c_ptr), value :: gloc

    end subroutine pdm_gnum_location_compute


    !>
    !!
    !! \brief Free
    !!
    !! \param [in]   gloc         Pointer to \ref PDM_gnum_locaion object
    !! \param [in]   partial      If partial = 1, location arrays are not freed
    !!
    !!/

    subroutine pdm_gnum_location_free (gloc) &
    bind (c, name = 'PDM_gnum_location_free')
      use iso_c_binding
      implicit none

      type (c_ptr), value :: gloc

    end subroutine pdm_gnum_location_free


    !>
    !!
    !! \brief Get the number of requested elements in a given partition
    !!
    !! \param [in]  gnum_loc      Pointer to \ref PDM_gnum_locaion object
    !! \param [in]  i_part_out    Id of current partition
    !!
    !! \return  Number of requested elements in the current partition
    !!

    function PDM_gnum_location_n_requested_elt_get (gnum_loc,   &
                                                    i_part_out) &
    result (n_elts_out)                                         &
    bind (c, name='PDM_gnum_location_n_requested_elt_get')
      use iso_c_binding
      implicit none

      type(c_ptr),    value :: gnum_loc
      integer(c_int), value :: i_part_out
      integer(c_int)        :: n_elts_out

    end function PDM_gnum_location_n_requested_elt_get

  end interface


  contains


  !>
  !!
  !! \brief Build a global numbering location structure
  !!
  !! \param [out]  gloc           Pointer to \ref PDM_gnum_locaion object
  !! \param [in]   n_part_in      Number of local partitions for elements
  !! \param [in]   n_part_out     Number of local partitions for requested locations
  !! \param [in]   f_comm         PDM_MPI communicator
  !! \param [in]   owner          Ownership
  !!

  subroutine pdm_gnum_location_create (gloc,       &
                                       n_part_in,  &
                                       n_part_out, &
                                       f_comm, &
                                       owner)
    use iso_c_binding
    implicit none

    integer, intent(in) :: n_part_in
    integer, intent(in) :: n_part_out
    integer, intent(in) :: f_comm
    integer, intent(in) :: owner
    type (c_ptr)        :: gloc

    integer(c_int)      :: c_n_part_in
    integer(c_int)      :: c_n_part_out
    integer(c_int)      :: c_comm
    integer(c_int)      :: c_owner

    interface
      function pdm_gnum_location_create_c (n_part_in,  &
                                           n_part_out, &
                                           comm,       &
                                           owner)      &
      result(gloc)                                     &
      bind (c, name = 'PDM_gnum_location_create')
        use iso_c_binding
        implicit none

        integer(c_int), value :: n_part_in
        integer(c_int), value :: n_part_out
        integer(c_int), value :: comm
        integer(c_int), value :: owner
        type (c_ptr)          :: gloc

      end function pdm_gnum_location_create_c
    end interface

    c_comm       = PDM_MPI_Comm_f2c(f_comm)
    c_n_part_in  = n_part_in
    c_n_part_out = n_part_out
    c_owner      = owner

    gloc = pdm_gnum_location_create_c (c_n_part_in,  &
                                       c_n_part_out, &
                                       c_comm,       &
                                       c_owner)

  end subroutine pdm_gnum_location_create


  !>
  !!
  !! \brief Set global numbering
  !!
  !! \param [in]   gloc        Pointer to \ref PDM_gnum_locaion object
  !! \param [in]   i_part_in   Current partition
  !! \param [in]   n_elts_in   Number of elements
  !! \param [in]   gnum_in     Global numbering
  !!

  subroutine pdm_gnum_location_elements_set (gloc,      &
                                             i_part_in, &
                                             n_elts_in, &
                                             gnum_in)
    use iso_c_binding
    implicit none

    type(c_ptr), value            :: gloc
    integer, intent(in)           :: i_part_in
    integer, intent(in)           :: n_elts_in
    integer(pdm_g_num_s), pointer :: gnum_in(:)

    interface
      subroutine pdm_gnum_location_elements_set_c (gloc,      &
                                                   i_part_in, &
                                                   n_elts_in, &
                                                   gnum_in)   &
      bind (c, name = 'PDM_gnum_location_elements_set')
        use iso_c_binding
        implicit none

        type (c_ptr),   value :: gloc
        integer(c_int), value :: i_part_in
        integer(c_int), value :: n_elts_in
        type(c_ptr),    value :: gnum_in

      end subroutine pdm_gnum_location_elements_set_c
    end interface

    call pdm_gnum_location_elements_set_c (gloc,           &
                                           i_part_in,      &
                                           n_elts_in,      &
                                           c_loc(gnum_in))

  end subroutine pdm_gnum_location_elements_set


  !>
  !!
  !! \brief Set requested elements
  !!
  !! \param [in]   gloc         Pointer to \ref PDM_gnum_locaion object
  !! \param [in]   i_part_out   Current partition
  !! \param [in]   n_elts_out   Number of elements
  !! \param [in]   gnum_out     Global numbering
  !!
  !!

  subroutine pdm_gnum_location_requested_elements_set (gloc,       &
                                                       i_part_out, &
                                                       n_elts_out, &
                                                       gnum_out)
    use iso_c_binding
    implicit none

    type(c_ptr), value            :: gloc
    integer, intent(in)           :: i_part_out
    integer, intent(in)           :: n_elts_out
    integer(pdm_g_num_s), pointer :: gnum_out(:)

    interface
      subroutine pdm_gnum_location_requested_elements_set_c (gloc,       &
                                                             i_part_out, &
                                                             n_elts_out, &
                                                             gnum_out)   &
      bind (c, name = 'PDM_gnum_location_requested_elements_set')
        use iso_c_binding
        implicit none

        type (c_ptr),   value :: gloc
        integer(c_int), value :: i_part_out
        integer(c_int), value :: n_elts_out
        type(c_ptr),    value :: gnum_out

      end subroutine pdm_gnum_location_requested_elements_set_c
    end interface

    call pdm_gnum_location_requested_elements_set_c (gloc,            &
                                                     i_part_out,      &
                                                     n_elts_out,      &
                                                     c_loc(gnum_out))

  end subroutine pdm_gnum_location_requested_elements_set


  !>
  !!
  !! \brief Get localtion
  !!
  !! \param [in]    gloc           Pointer to \ref PDM_gnum_locaion object
  !! \param [in]    i_part_out     Current partition
  !! \param [out]   location_idx   Index in the location arrays (size = \ref n_elts + 1)
  !! \param [out]   location       Locations of each element
  !!                                (Three informations : process, partition, element)
  !!
  subroutine pdm_gnum_location_get (gloc,         &
                                    i_part_out,   &
                                    location_idx, &
                                    location)
    use iso_c_binding
    implicit none

    type(c_ptr), value            :: gloc
    integer, intent(in)           :: i_part_out
    integer(pdm_l_num_s), pointer :: location_idx(:)
    integer(pdm_l_num_s), pointer :: location(:)

    type(c_ptr)                   :: c_location_idx
    type(c_ptr)                   :: c_location    
    integer                       :: n_elt

    interface
      subroutine pdm_gnum_location_get_c (gloc,         &
                                          i_part_out,   &
                                          location_idx, &
                                          location)     &
      bind (c, name = 'PDM_gnum_location_get')
        use iso_c_binding
        implicit none

        type (c_ptr),   value :: gloc
        integer(c_int), value :: i_part_out
        type(c_ptr)           :: location_idx
        type(c_ptr)           :: location

      end subroutine pdm_gnum_location_get_c
    end interface

    c_location_idx = C_NULL_PTR
    c_location     = C_NULL_PTR

    call pdm_gnum_location_get_c (gloc,           &
                                  i_part_out,     &
                                  c_location_idx, &
                                  c_location)

    n_elt = pdm_gnum_location_n_requested_elt_get (gloc,       &
                                                   i_part_out)

    call c_f_pointer(c_location_idx, &
                     location_idx,   &
                     [n_elt+1])

    call c_f_pointer(c_location,               &
                     location,                 &
                     [location_idx(n_elt+1)])

  end subroutine pdm_gnum_location_get

end module pdm_gnum_location
