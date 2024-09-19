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

module pdm_closest_points

  use pdm
  implicit none

  interface PDM_closest_points_create ; module procedure &
  pdm_closest_points_create_
  end interface

  interface PDM_closest_points_tgt_cloud_set ; module procedure &
  pdm_closest_points_tgt_cloud_set_
  end interface

  interface PDM_closest_points_src_cloud_set ; module procedure &
  pdm_closest_points_src_cloud_set_
  end interface

  interface PDM_closest_points_get ; module procedure &
  pdm_closest_points_get_
  end interface

  interface PDM_closest_points_tgt_in_src_get ; module procedure &
  pdm_closest_points_tgt_in_src_get_
  end interface

  private :: pdm_closest_points_create_
  private :: pdm_closest_points_tgt_cloud_set_
  private :: pdm_closest_points_src_cloud_set_
  private :: pdm_closest_points_get_
  private :: pdm_closest_points_tgt_in_src_get_

  interface

    !>
    !!
    !! \brief Create a structure to look for the closest points of a point cloud
    !! (target cloud) in an other point cloud (source cloud)
    !!
    !! \param [in]   comm           MPI communicator
    !! \param [in]   n_closest      Number of closest source points to find for each
    !!                              target point
    !!
    !! \return     Pointer to \ref PDM_closest_points object
    !!
    !!

    function PDM_closest_points_create_cf (comm,      &
                                           n_closest, &
                                           owner)     &
    result(cls)                                       &
    bind (c, name = 'PDM_closest_points_create')

      use iso_c_binding

      implicit none

      integer(c_int), value :: comm
      integer(c_int), value :: n_closest
      integer(c_int), value :: owner

      type(c_ptr)           :: cls

    end function PDM_closest_points_create_cf

    !>
    !!
    !! \brief Set the number of partitions of a point cloud
    !!
    !! \param [in]   cls               Pointer to \ref PDM_closest_points object
    !! \param [in]   n_part_cloud_src  Number of partitions of the source cloud
    !! \param [in]   n_part_cloud_tgt  Number of partitions od the target cloud
    !!
    !!

    subroutine PDM_closest_points_n_part_cloud_set (cls, &
                                                    n_part_cloud_src, &
                                                    n_part_cloud_tgt) &
      bind (c, name = 'PDM_closest_points_n_part_cloud_set')

      use iso_c_binding

      implicit none

      type(c_ptr), value :: cls

      integer(c_int), value :: n_part_cloud_src
      integer(c_int), value :: n_part_cloud_tgt


    end subroutine PDM_closest_points_n_part_cloud_set

    !>
    !!
    !! \brief Set the target point cloud
    !!
    !! \param [in]   cls             Pointer to \ref PDM_closest_points object
    !! \param [in]   i_part          Index of partition
    !! \param [in]   n_points        Number of points
    !! \param [in]   coords          Point coordinates
    !! \param [in]   gnum            Point global number
    !!
    !!

    subroutine PDM_closest_points_tgt_cloud_set_cf (cls, &
                                                    i_part, &
                                                    n_points, &
                                                    coords, &
                                                    gnum) &
      bind (c, name = 'PDM_closest_points_tgt_cloud_set')

      use iso_c_binding

      implicit none

      type(c_ptr), value :: cls

      integer(c_int), value :: i_part
      integer(c_int), value :: n_points

      type(c_ptr), value    :: coords
      type(c_ptr), value    :: gnum


    end subroutine PDM_closest_points_tgt_cloud_set_cf

    !>
    !!
    !! \brief Set the source point cloud
    !!
    !! \param [in]   cls             Pointer to \ref PDM_closest_points object
    !! \param [in]   i_part          Index of partition
    !! \param [in]   n_points        Number of points
    !! \param [in]   coords          Point coordinates
    !! \param [in]   gnum            Point global number
    !!
    !!

    subroutine PDM_closest_points_src_cloud_set_cf (cls, &
                                                    i_part, &
                                                    n_points, &
                                                    coords, &
                                                    gnum) &
      bind (c, name = 'PDM_closest_points_src_cloud_set')

      use iso_c_binding

      implicit none

      type(c_ptr), value :: cls

      integer(c_int), value :: i_part
      integer(c_int), value :: n_points

      type(c_ptr), value    :: coords
      type(c_ptr), value    :: gnum


    end subroutine PDM_closest_points_src_cloud_set_cf

    !>
    !!
    !! \brief Look for closest points
    !!
    !! \param [in]   cls Pointer to \ref PDM_closest_points object
    !!
    !!

    subroutine PDM_closest_points_compute (cls) &
      bind (c, name = 'PDM_closest_points_compute')

      use iso_c_binding

      implicit none

      type(c_ptr), value :: cls

    end subroutine PDM_closest_points_compute

    !>
    !!
    !! \brief Get mesh distance
    !!
    !! \param [in]   cls                   Pointer to \ref PDM_closest_points object
    !! \param [in]   i_part_tgt            Index of partition of the cloud
    !! \param [out]  closest_src_g_num     Global number of the closest element (size = n_closest * n_tgt_points)
    !! \param [out]  closest_src_distance  Distance (size = n_closest * n_tgt_points)
    !!
    !!

    subroutine PDM_closest_points_get_cf (cls, &
                                          i_part_tgt, &
                                          closest_src_gnum, &
                                          closest_src_distance) &
      bind (c, name = 'PDM_closest_points_get')

      use iso_c_binding

      implicit none

      type(c_ptr), value :: cls

      integer(c_int), value :: i_part_tgt

      type(c_ptr)      :: closest_src_gnum
      type(c_ptr)      :: closest_src_distance


    end subroutine PDM_closest_points_get_cf

    !>
    !!
    !! \brief Get Get closest source points global ids and (squared) distance
    !!
    !! \param [in]   cls                   Pointer to \ref PDM_closest_points object
    !! \param [in]   i_part_src            Index of partition of the cloud
    !! \param [out]  tgt_in_src_idx        For each src point the number of target localised  (size = n_src_points )
    !! \param [out]  tgt_in_src            For each src point the globla number of target point located (size = tgt_in_src_idx[n_src_points] )
    !!
    !!

    subroutine PDM_closest_points_tgt_in_src_get_cf (cls, &
                                                     i_part_src, &
                                                     tgt_in_src_idx, &
                                                     tgt_in_src) &
      bind (c, name = 'PDM_closest_points_tgt_in_src_get')

      use iso_c_binding

      implicit none

      type(c_ptr), value :: cls

      integer(c_int), value :: i_part_src

      type(c_ptr)      :: tgt_in_src_idx
      type(c_ptr)      :: tgt_in_src


    end subroutine PDM_closest_points_tgt_in_src_get_cf

    !>
    !!
    !! \brief Free a distance mesh structure
    !!
    !! \param [in]  cls      Pointer to \ref PDM_closest_points object
    !! \param [in]  partial  if partial is equal to 0, all data are removed.
    !!                       Otherwise, results are kept.
    !!
    !!

    subroutine PDM_closest_points_free (cls) &
     bind (c, name = 'PDM_closest_points_free')

      use iso_c_binding

      implicit none

      type(c_ptr), value :: cls

    end subroutine PDM_closest_points_free

    !>
    !!
    !! \brief  Dump elapsed an CPU time
    !!
    !! \param [in]  cls      Pointer to \ref PDM_closest_points object
    !!
    !!

    subroutine PDM_closest_points_dump_times (cls) &
      bind (c, name = 'PDM_closest_points_dump_times')

      use iso_c_binding

      implicit none

      type(c_ptr), value :: cls

    end subroutine PDM_closest_points_dump_times


!>
!!
!! \brief  Get the number of target points in a partition
!!
!! \param [in]  cls     Pointer to \ref PDM_closest_points_t object
!! \param [in]  i_part  Index of partition of the target cloud
!!
!! \return   Number of target point in the partition \ref i_part
!!
!!

function PDM_closest_points_n_tgt_get (cls,     &
                                       i_part)  &
  result (n_tgt)                                &

  bind (c, name='PDM_closest_points_n_tgt_get')

  use iso_c_binding
  implicit none

  type(c_ptr),    value :: cls
  integer(c_int), value :: i_part
  integer(c_int)        :: n_tgt

end function PDM_closest_points_n_tgt_get


!>
!!
!! \brief  Get the number of source points in a partition
!!
!! \param [in]  cls     Pointer to \ref PDM_closest_points_t object
!! \param [in]  i_part  Index of partition of the target cloud
!!
!! \return   Number of source point in the partition \ref i_part
!!
!!

function PDM_closest_points_n_src_get (cls,     &
                                       i_part)  &
  result (n_src)                                &

  bind (c, name='PDM_closest_points_n_src_get')

  use iso_c_binding
  implicit none

  type(c_ptr),    value :: cls
  integer(c_int), value :: i_part
  integer(c_int)        :: n_src

end function PDM_closest_points_n_src_get


!>
!!
!! \brief  Get the number of closest points
!!
!! \param [in]  cls     Pointer to \ref PDM_closest_points_t object
!!
!! \return   Number of closest points
!!
!!

function PDM_closest_points_n_closest_get (cls)  &
  result (n_closest)                             &

  bind (c, name='PDM_closest_points_n_closest_get')

  use iso_c_binding
  implicit none

  type(c_ptr),    value :: cls
  integer(c_int)        :: n_closest

end function PDM_closest_points_n_closest_get




  end interface


  contains


    !>
    !!
    !! \brief Create a structure to look for the closest points of a point cloud
    !! (target cloud) in an other point cloud (source cloud)
    !!
    !! \param [out]  cls         Pointer to \ref PDM_closest_points object
    !! \param [in]   f_comm      MPI communicator
    !! \param [in]   n_closest   Number of closest source points to find for each
    !!                           target point
    !!

    subroutine PDM_closest_points_create_ (cls,       &
                                           f_comm,    &
                                           n_closest, &
                                           owner)

      use iso_c_binding

      implicit none

      integer, intent(in) :: f_comm
      integer, intent(in) :: n_closest
      integer, intent(in) :: owner

      type(c_ptr)    :: cls

      integer(c_int) :: c_comm
      integer(c_int) :: c_n_closest
      integer(c_int) :: c_owner

      c_comm      = PDM_MPI_Comm_f2c(f_comm)
      c_n_closest = n_closest
      c_owner     = owner

      cls = PDM_closest_points_create_cf(c_comm,      &
                                         c_n_closest, &
                                         c_owner)

    end subroutine PDM_closest_points_create_





  !>
  !!
  !! \brief Set the target point cloud
  !!
  !! \param [in]   cls             Pointer to \ref PDM_closest_points object
  !! \param [in]   i_part          Index of partition
  !! \param [in]   n_points        Number of points
  !! \param [in]   coords          Point coordinates
  !! \param [in]   gnum            Point global number
  !!
  !!

  subroutine PDM_closest_points_tgt_cloud_set_ (cls,      &
                                                i_part,   &
                                                n_points, &
                                                coords,   &
                                                gnum)

    use iso_c_binding

    implicit none

    type(c_ptr), value                 :: cls
    integer, intent(in)                :: i_part
    integer, intent(in)                :: n_points
    double precision,          pointer :: coords(:,:)
    integer(kind=pdm_g_num_s), pointer :: gnum(:)

    integer(c_int)                     :: c_i_part
    integer(c_int)                     :: c_n_points
    type(c_ptr)                        :: c_coords
    type(c_ptr)                        :: c_gnum

    c_i_part   = i_part
    c_n_points = n_points

    c_coords = C_NULL_PTR
    if (associated(coords)) then
      c_coords = c_loc(coords)
    end if

    c_gnum = C_NULL_PTR
    if (associated(coords)) then
      c_gnum = c_loc(gnum)
    end if

    call PDM_closest_points_tgt_cloud_set_cf(cls,        &
                                             c_i_part,   &
                                             c_n_points, &
                                             c_coords,   &
                                             c_gnum)

  end subroutine PDM_closest_points_tgt_cloud_set_



  !>
  !!
  !! \brief Set the source point cloud
  !!
  !! \param [in]   cls             Pointer to \ref PDM_closest_points object
  !! \param [in]   i_part          Index of partition
  !! \param [in]   n_points        Number of points
  !! \param [in]   coords          Point coordinates
  !! \param [in]   gnum            Point global number
  !!
  !!

  subroutine PDM_closest_points_src_cloud_set_ (cls, &
                                                i_part, &
                                                n_points, &
                                                coords, &
                                                gnum)

    use iso_c_binding

    type(c_ptr), value                 :: cls
    integer, intent(in)                :: i_part
    integer, intent(in)                :: n_points
    double precision,          pointer :: coords(:,:)
    integer(kind=pdm_g_num_s), pointer :: gnum(:)

    integer(c_int)                     :: c_i_part
    integer(c_int)                     :: c_n_points
    type(c_ptr)                        :: c_coords
    type(c_ptr)                        :: c_gnum

    c_i_part   = i_part
    c_n_points = n_points

    c_coords = C_NULL_PTR
    if (associated(coords)) then
      c_coords = c_loc(coords)
    end if

    c_gnum = C_NULL_PTR
    if (associated(coords)) then
      c_gnum = c_loc(gnum)
    end if

    call PDM_closest_points_src_cloud_set_cf(cls,        &
                                             c_i_part,   &
                                             c_n_points, &
                                             c_coords,   &
                                             c_gnum)


  end subroutine PDM_closest_points_src_cloud_set_



  !>
  !!
  !! \brief Get Get closest source points global ids and (squared) distance
  !!
  !! \param [in]   cls                   Pointer to \ref PDM_closest_points object
  !! \param [in]   i_part_tgt            Index of partition of the cloud
  !! \param [out]  closest_src_g_num     Global number of the closest element (size = n_closest * n_tgt_points)
  !! \param [out]  closest_src_distance  Distance (size = n_closest * n_tgt_points)
  !!
  !!

  subroutine PDM_closest_points_get_ (cls, &
                                      i_part_tgt, &
                                      closest_src_gnum, &
                                      closest_src_distance)

    use iso_c_binding

    implicit none

    type(c_ptr),                 value :: cls
    integer, intent(in)                :: i_part_tgt
    integer(kind=pdm_g_num_s), pointer :: closest_src_gnum(:)
    double precision,          pointer :: closest_src_distance(:)

    integer(c_int)                     :: c_i_part_tgt
    type(c_ptr)                        :: c_closest_src_gnum
    type(c_ptr)                        :: c_closest_src_distance
    integer(c_int)                     :: n_tgt, n_closest

    c_i_part_tgt = i_part_tgt

    c_closest_src_gnum     = C_NULL_PTR
    c_closest_src_distance = C_NULL_PTR
    call PDM_closest_points_get_cf (cls,                    &
                                    c_i_part_tgt,           &
                                    c_closest_src_gnum,     &
                                    c_closest_src_distance)

    n_tgt = pdm_closest_points_n_tgt_get(cls,        &
                                         i_part_tgt)

    n_closest = pdm_closest_points_n_closest_get(cls)

    call c_f_pointer(c_closest_src_gnum, &
                     closest_src_gnum,   &
                     [n_tgt*n_closest])

    call c_f_pointer(c_closest_src_distance, &
                     closest_src_distance,   &
                     [n_tgt*n_closest])

  end subroutine PDM_closest_points_get_

  !>
  !!
  !! \brief Get Get closest source points global ids and (squared) distance
  !!
  !! \param [in]   cls                   Pointer to \ref PDM_closest_points object
  !! \param [in]   i_part_src            Index of partition of the cloud
  !! \param [out]  tgt_in_src_idx        For each src point the number of target localised  (size = n_src_points )
  !! \param [out]  tgt_in_src            For each src point the globla number of target point located (size = tgt_in_src_idx[n_src_points] )
  !!
  !!

  subroutine PDM_closest_points_tgt_in_src_get_ (cls, &
                                                 i_part_src, &
                                                 tgt_in_src_idx, &
                                                 tgt_in_src)

    use iso_c_binding

    implicit none

    type(c_ptr),                 value :: cls
    integer, intent(in)                :: i_part_src
    integer(kind=pdm_l_num_s), pointer :: tgt_in_src_idx(:)
    integer(kind=pdm_g_num_s), pointer :: tgt_in_src(:)

    integer(c_int)                     :: c_i_part_src
    type(c_ptr)                        :: c_tgt_in_src_idx
    type(c_ptr)                        :: c_tgt_in_src
    integer(c_int)                     :: n_src

    c_i_part_src = i_part_src

    c_tgt_in_src_idx = C_NULL_PTR
    c_tgt_in_src     = C_NULL_PTR
    call PDM_closest_points_tgt_in_src_get_cf (cls,                  &
                                               i_part_src,           &
                                               c_tgt_in_src_idx,     &
                                               c_tgt_in_src)

    n_src = pdm_closest_points_n_src_get(cls,        &
                                         i_part_src)

    call c_f_pointer(c_tgt_in_src_idx, &
                     tgt_in_src_idx,   &
                     [n_src+1])

    call c_f_pointer(c_tgt_in_src, &
                     tgt_in_src,   &
                     [tgt_in_src_idx(n_src+1)])

  end subroutine PDM_closest_points_tgt_in_src_get_

end module pdm_closest_points
