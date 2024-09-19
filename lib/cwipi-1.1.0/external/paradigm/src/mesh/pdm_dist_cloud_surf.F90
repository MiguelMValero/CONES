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

module pdm_dist_cloud_surf

  use pdm
  use pdm_pointer_array

  implicit none

  interface PDM_dist_cloud_surf_create ; module procedure &
  pdm_dist_cloud_surf_create_
  end interface

  interface PDM_dist_cloud_surf_get ; module procedure &
  pdm_dist_cloud_surf_get_
  end interface

  interface PDM_dist_cloud_surf_distri_data ; module procedure &
  pdm_dist_cloud_surf_distri_data_
  end interface

  interface PDM_dist_cloud_surf_cloud_set ; module procedure &
  pdm_dist_cloud_surf_cloud_set_
  end interface

  interface PDM_dist_cloud_surf_surf_mesh_part_set ; module procedure &
  pdm_dist_cloud_surf_surf_mesh_part_set_
  end interface

  private :: pdm_dist_cloud_surf_create_
  private :: pdm_dist_cloud_surf_get_
  private :: pdm_dist_cloud_surf_cloud_set_
  private :: pdm_dist_cloud_surf_surf_mesh_part_set_

  interface

    !> \brief Create a structure to compute distance to a mesh nodal
    !!
    !! \param [in]    mesh_nature    Nature of the mesh
    !! \param [in]    n_point_cloud  Number of point cloud
    !! \param [in]    comm           MPI communicator
    !!
    !! \return    Pointer to \ref PDM_dist_cloud_surf object
    !!

    function pdm_dist_cloud_surf_create_cf (mesh_nature,   &
                                            n_point_cloud, &
                                            comm,          &
                                            owner)         &
                                       result(dcs)         &
         bind (c, name = 'PDM_dist_cloud_surf_create')

      use iso_c_binding

      implicit none

      integer(c_int), value :: mesh_nature
      integer(c_int), value :: n_point_cloud
      integer(c_int), value :: comm
      integer(c_int), value :: owner

      type (c_ptr) :: dcs

    end function pdm_dist_cloud_surf_create_cf

    !> \brief Set the number of partitions of a point cloud
    !!
    !! \param [in]   dcs             Pointer to \ref PDM_dist_cloud_surf object
    !! \param [in]   i_point_cloud   Index of point cloud
    !! \param [in]   n_part          Number of partitions
    !!

    subroutine pdm_dist_cloud_surf_n_part_cloud_set (dcs, i_point_cloud, n_part) &
         bind (c, name = 'PDM_dist_cloud_surf_n_part_cloud_set')

      use iso_c_binding

      implicit none

      type (c_ptr), value          :: dcs
      integer(c_int), value        :: i_point_cloud
      integer(c_int), value        :: n_part

    end subroutine pdm_dist_cloud_surf_n_part_cloud_set

    !> \brief Set a point cloud
    !!
    !! \param [in]   dcs             Pointer to \ref PDM_dist_cloud_surf object
    !! \param [in]   i_point_cloud   Index of point cloud
    !! \param [in]   i_part          Index of partition
    !! \param [in]   n_points        Number of points
    !! \param [in]   coords          Point coordinates
    !! \param [in]   gnum            Point global number
    !!

    subroutine pdm_dist_cloud_surf_cloud_set_cf (dcs,           &
                                                 i_point_cloud, &
                                                 i_part,        &
                                                 n_points,      &
                                                 coords,        &
                                                 gnum)          &
         bind (c, name = 'PDM_dist_cloud_surf_cloud_set')

      use iso_c_binding

      implicit none

      type (c_ptr), value          :: dcs
      integer(c_int), value        :: i_point_cloud
      integer(c_int), value        :: i_part
      integer(c_int), value        :: n_points

      type(c_ptr), value        :: coords
      type(c_ptr), value        :: gnum

    end subroutine pdm_dist_cloud_surf_cloud_set_cf

    !> \brief Set the mesh nodal
    !!
    !! \param [in]   dcs            Pointer to \ref PDM_dist_cloud_surf object
    !! \param [in]   mesh_nodal_id  Mesh nodal Pointer to \ref PDM_dist_cloud_surf object
    !!
    !!

    subroutine pdm_dist_cloud_surf_nodal_mesh_set (dcs, mesh_nodal_id) &
         bind (c, name = 'pdm_dist_cloud_surf_nodal_mesh_set')

      use iso_c_binding

      implicit none

      type (c_ptr), value          :: dcs
      integer(c_int), value        :: mesh_nodal_id

    end subroutine pdm_dist_cloud_surf_nodal_mesh_set

    !> \brief Set global data of a surface mesh
    !!
    !! \param [in]   dcs            Pointer to \ref PDM_dist_cloud_surf object
    !! \param [in]   n_part         Number of partition
    !!

    subroutine pdm_dist_cloud_surf_surf_mesh_global_data_set (dcs, n_part) &
         bind (c, name = 'PDM_dist_cloud_surf_surf_mesh_global_data_set')

      use iso_c_binding

      implicit none

      type (c_ptr), value          :: dcs
      integer(c_int), value        :: n_part
    end subroutine pdm_dist_cloud_surf_surf_mesh_global_data_set

    !> \brief Set a part of a surface mesh
    !!
    !! \param [in]   dcs           Pointer to \ref PDM_dist_cloud_surf object
    !! \param [in]   i_part        Partition to define
    !! \param [in]   n_face        Number of faces
    !! \param [in]   face_vtx_idx  Index in the face -> vertex connectivity
    !! \param [in]   face_vtx      face -> vertex connectivity
    !! \param [in]   face_ln_to_gn Local face numbering to global face numbering
    !! \param [in]   n_vtx         Number of vertices
    !! \param [in]   coords        Coordinates
    !! \param [in]   vtx_ln_to_gn  Local vertex numbering to global vertex numbering
    !!

    subroutine pdm_dist_cloud_surf_surf_mesh_part_set_cf (dcs,           &
                                                          i_part,        &
                                                          n_face,        &
                                                          face_vtx_idx,  &
                                                          face_vtx,      &
                                                          face_ln_to_gn, &
                                                          n_vtx,         &
                                                          coords,        &
                                                          vtx_ln_to_gn)  &
      bind (c, name = 'PDM_dist_cloud_surf_surf_mesh_part_set')
      use iso_c_binding

      implicit none

      type (c_ptr),   value  :: dcs
      integer(c_int), value  :: i_part
      integer(c_int), value  :: n_face
      type(c_ptr),    value  :: face_vtx_idx
      type(c_ptr),    value  :: face_vtx
      type(c_ptr),    value  :: face_ln_to_gn
      integer(c_int), value  :: n_vtx
      type(c_ptr),    value  :: coords
      type(c_ptr),    value  :: vtx_ln_to_gn

    end subroutine pdm_dist_cloud_surf_surf_mesh_part_set_cf

    !> \brief Compute distance
    !!
    !! \param [in]   dcs Pointer to \ref PDM_dist_cloud_surf object
    !!

    subroutine pdm_dist_cloud_surf_compute (dcs) &
         bind (c, name = 'PDM_dist_cloud_surf_compute')
      use iso_c_binding

      implicit none

      type (c_ptr), value       :: dcs
    end subroutine pdm_dist_cloud_surf_compute

    !> \brief Get mesh distance
    !!
    !! \param [in]   dcs                   Pointer to \ref PDM_dist_cloud_surf object
    !! \param [in]   i_point_cloud         Current cloud
    !! \param [in]   i_part                Index of partition of the cloud
    !! \param [out]  closest_elt_distance  Distance
    !! \param [out]  closest_elt_projected Projected point coordinates
    !! \param [out]  closest_elt_g_num     Global number of the closest element
    !!

    subroutine pdm_dist_cloud_surf_get_cf (dcs,                   &
                                           i_point_cloud,         &
                                           i_part,                &
                                           closest_elt_distance,  &
                                           closest_elt_projected, &
                                           closest_elt_gnum)      &
      bind (c, name = 'PDM_dist_cloud_surf_get')

      use iso_c_binding

      implicit none

      type (c_ptr), value       :: dcs
      integer(c_int), value     :: i_point_cloud
      integer(c_int), value     :: i_part
      type(c_ptr)               :: closest_elt_distance
      type(c_ptr)               :: closest_elt_projected
      type(c_ptr)               :: closest_elt_gnum

    end subroutine pdm_dist_cloud_surf_get_cf

    !> \brief Get mesh distance
    !!
    !! \param [in]   dcs                   Pointer to \ref PDM_dist_cloud_surf object
    !! \param [in]   i_point_cloud         Current cloud
    !! \param [in]   stride                Stride
    !! \param [in]   surf_data             Data over the surface
    !! \param [out]  cloud_data            Data over the cloud
    !!

    subroutine pdm_dist_cloud_surf_distri_data_cf (dcs,                   &
                                                   i_point_cloud,         &
                                                   stride,                &
                                                   surf_data,             &
                                                   cloud_data)            &
      bind (c, name = 'PDM_dist_cloud_surf_distri_data')

      use iso_c_binding

      implicit none

      type (c_ptr), value       :: dcs
      integer(c_int), value     :: i_point_cloud
      integer(c_int), value     :: stride
      type(c_ptr), value        :: surf_data
      type(c_ptr)               :: cloud_data

    end subroutine pdm_dist_cloud_surf_distri_data_cf

    !> \brief Free a distance mesh structure
    !!
    !! \param [in]  dcs      Pointer to \ref PDM_dist_cloud_surf object
    !! \param [in]  partial  if partial is equal to 0, all data are removed.
    !!                       Otherwise, results are kept.
    !!

    subroutine pdm_dist_cloud_surf_free (dcs, partial) &
       bind (c, name = 'PDM_dist_cloud_surf_free')

      use iso_c_binding

      implicit none

      type (c_ptr), value       :: dcs
      integer(c_int), value     :: partial

    end subroutine pdm_dist_cloud_surf_free

    !> \brief  Dump elapsed an CPU time
    !!
    !! \param [in]  dcs      Pointer to \ref PDM_dist_cloud_surf object
    !!

    subroutine pdm_dist_cloud_surf_dump_times (dcs) &
         bind (c, name = 'PDM_dist_cloud_surf_dump_times')

      use iso_c_binding

      implicit none

      type (c_ptr), value       :: dcs

    end subroutine pdm_dist_cloud_surf_dump_times


  !>
  !!
  !! \brief Get the dimension of a point cloud
  !!
  !! \param [in]   dcs             Pointer to \ref PDM_dist_cloud_surf object
  !! \param [in]   i_point_cloud   Index of point cloud
  !! \param [in]   i_part          Index of partition
  !! \param [out]  n_points        Number of points
  !!
  !!

  subroutine PDM_dist_cloud_surf_cloud_dim_get (dcs,           &
                                                i_point_cloud, &
                                                i_part,        &
                                                n_points)      &

    bind (c, name = 'PDM_dist_cloud_surf_cloud_dim_get')

    use iso_c_binding
    implicit none

    type(c_ptr), value    :: dcs
    integer(c_int), value :: i_point_cloud
    integer(c_int), value :: i_part
    integer(c_int)        :: n_points

  end subroutine PDM_dist_cloud_surf_cloud_dim_get


  end interface


  contains


  !> \brief Create a structure to compute distance to a mesh nodal
  !!
  !! \param [out] dcs            Pointer to \ref PDM_dist_cloud_surf object
  !! \param [in]  mesh_nature    Nature of the mesh
  !! \param [in]  n_point_cloud  Number of point cloud
  !! \param [in]  comm           MPI communicator
  !!

  subroutine pdm_dist_cloud_surf_create_ (dcs,           &
                                          mesh_nature,   &
                                          n_point_cloud, &
                                          f_comm,        &
                                          owner)

    use iso_c_binding

    implicit none

    integer, intent(in) :: mesh_nature
    integer, intent(in) :: n_point_cloud
    integer, intent(in) :: f_comm
    integer, intent(in) :: owner

    type (c_ptr) :: dcs

    integer(c_int) :: c_mesh_nature
    integer(c_int) :: c_n_point_cloud
    integer(c_int) :: c_comm
    integer(c_int) :: c_owner

    c_comm = PDM_MPI_Comm_f2c(f_comm)

    c_mesh_nature   = mesh_nature
    c_n_point_cloud = n_point_cloud
    c_owner         = owner

    dcs = pdm_dist_cloud_surf_create_cf(c_mesh_nature,   &
                                        c_n_point_cloud, &
                                        c_comm,          &
                                        c_owner)

  end subroutine pdm_dist_cloud_surf_create_



  !> \brief Get mesh distance
  !!
  !! \param [in]   dcs                   Pointer to \ref PDM_dist_cloud_surf object
  !! \param [in]   i_point_cloud         Current cloud
  !! \param [in]   i_part                Index of partition of the cloud
  !! \param [out]  closest_elt_distance  Distance
  !! \param [out]  closest_elt_projected Projected point coordinates
  !! \param [out]  closest_elt_g_num     Global number of the closest element
  !!

  subroutine pdm_dist_cloud_surf_get_ (dcs,                   &
                                       i_point_cloud,         &
                                       i_part,                &
                                       closest_elt_distance,  &
                                       closest_elt_projected, &
                                       closest_elt_gnum)

    use iso_c_binding

    implicit none

    type (c_ptr), value                  :: dcs
    integer, intent(in)                  :: i_point_cloud
    integer, intent(in)                  :: i_part
    double precision,            pointer :: closest_elt_distance(:)
    double precision,            pointer :: closest_elt_projected(:,:)
    integer(kind = pdm_g_num_s), pointer :: closest_elt_gnum(:)

    integer(c_int)                       :: c_i_point_cloud
    integer(c_int)                       :: c_i_part
    type(c_ptr)                          :: c_closest_elt_distance
    type(c_ptr)                          :: c_closest_elt_projected
    type(c_ptr)                          :: c_closest_elt_gnum
    integer                              :: n_points

    c_i_point_cloud = i_point_cloud
    c_i_part        = i_part

    call PDM_dist_cloud_surf_cloud_dim_get (dcs,             &
                                            c_i_point_cloud, &
                                            c_i_part,        &
                                            n_points)

    call pdm_dist_cloud_surf_get_cf (dcs,                     &
                                     c_i_point_cloud,         &
                                     c_i_part,                &
                                     c_closest_elt_distance,  &
                                     c_closest_elt_projected, &
                                     c_closest_elt_gnum)

    call c_f_pointer(c_closest_elt_distance, &
                     closest_elt_distance,   &
                     [n_points])

    call c_f_pointer(c_closest_elt_projected, &
                     closest_elt_projected,   &
                     [3,n_points])

    call c_f_pointer(c_closest_elt_gnum, &
                     closest_elt_gnum,   &
                     [n_points])

  end subroutine pdm_dist_cloud_surf_get_



  !> \brief Get mesh distance
  !!
  !! \param [in]   dcs                   Pointer to \ref PDM_dist_cloud_surf object
  !! \param [in]   i_point_cloud         Current cloud
  !! \param [in]   stride                Stride
  !! \param [in]   surf_data             Data over the surface
  !! \param [out]  cloud_data            Data over the cloud
  !!

  subroutine pdm_dist_cloud_surf_distri_data_ (dcs,                   &
                                               i_point_cloud,         &
                                               stride,                &
                                               surf_data,             &
                                               cloud_data)

    use iso_c_binding

    implicit none

    type (c_ptr), value                  :: dcs
    integer, intent(in)                  :: i_point_cloud
    integer, intent(in)                  :: stride
    type(PDM_pointer_array_t), pointer   :: surf_data
    type(PDM_pointer_array_t), pointer   :: cloud_data

    integer(c_int)                       :: c_i_point_cloud
    integer(c_int)                       :: c_i_part
    integer(c_int)                       :: c_stride
    type(c_ptr)                          :: c_cloud_data
    integer                              :: n_points
    integer                              :: n_part
    integer                              :: i_part
    integer, allocatable                 :: length_data(:)

    interface

      function PDM_dist_cloud_surf_cloud_n_part_get_c (dcs, i_point_cloud) result (n_part) &
        bind(c, name='PDM_dist_cloud_surf_cloud_n_part_get')

        use iso_c_binding
        use pdm

        implicit none
        type(c_ptr),                value :: dcs
        integer(c_int), intent(in), value :: i_point_cloud
        integer(c_int)                    :: n_part
      end function

    end interface

    c_i_point_cloud = i_point_cloud
    c_stride        = stride

    n_part = PDM_dist_cloud_surf_cloud_n_part_get_c (dcs, i_point_cloud)

    call pdm_dist_cloud_surf_distri_data_cf (dcs,                    &
                                             c_i_point_cloud,        &
                                             c_stride,               &
                                             c_loc(surf_data%cptr),  &
                                             c_cloud_data)

    allocate( length_data(n_part) )

    do i_part = 1, n_part

      c_i_part = i_part-1

      call PDM_dist_cloud_surf_cloud_dim_get (dcs,             &
                                              c_i_point_cloud, &
                                              c_i_part,        &
                                              n_points)
      length_data(i_part) = n_points*stride

    end do

    call PDM_pointer_array_create (cloud_data,         &
                                   n_part,             &
                                   surf_data%type,     &
                                   c_cloud_data,       &
                                   length_data,        &
                                   PDM_OWNERSHIP_KEEP, &
                                   surf_data%s_data)

    deallocate( length_data )

  end subroutine pdm_dist_cloud_surf_distri_data_



  !> \brief Set a point cloud
  !!
  !! \param [in]   dcs             Pointer to \ref PDM_dist_cloud_surf object
  !! \param [in]   i_point_cloud   Index of point cloud
  !! \param [in]   i_part          Index of partition
  !! \param [in]   n_points        Number of points
  !! \param [in]   coords          Point coordinates
  !! \param [in]   gnum            Point global number
  !!

  subroutine pdm_dist_cloud_surf_cloud_set_ (dcs,           &
                                             i_point_cloud, &
                                             i_part,        &
                                             n_points,      &
                                             coords,        &
                                             gnum)

    use iso_c_binding

    implicit none

    type (c_ptr), value                :: dcs
    integer, intent(in)                :: i_point_cloud
    integer, intent(in)                :: i_part
    integer, intent(in)                :: n_points
    double precision,          pointer :: coords(:,:)
    integer(kind=pdm_g_num_s), pointer :: gnum(:)

    integer(c_int) :: c_i_point_cloud
    integer(c_int) :: c_i_part
    integer(c_int) :: c_n_points

    type(c_ptr)    :: c_coords
    type(c_ptr)    :: c_gnum


    c_i_point_cloud = i_point_cloud
    c_i_part        = i_part
    c_n_points      = n_points

    c_coords = c_loc(coords)
    c_gnum   = c_loc(gnum)

    call pdm_dist_cloud_surf_cloud_set_cf(dcs,             &
                                          c_i_point_cloud, &
                                          c_i_part,        &
                                          c_n_points,      &
                                          c_coords,        &
                                          c_gnum)

  end subroutine pdm_dist_cloud_surf_cloud_set_



  !> \brief Set a part of a surface mesh
  !!
  !! \param [in]   dcs           Pointer to \ref PDM_dist_cloud_surf object
  !! \param [in]   i_part        Partition to define
  !! \param [in]   n_face        Number of faces
  !! \param [in]   face_vtx_idx  Index in the face -> vertex connectivity
  !! \param [in]   face_vtx      face -> vertex connectivity
  !! \param [in]   face_ln_to_gn Local face numbering to global face numbering
  !! \param [in]   n_vtx         Number of vertices
  !! \param [in]   coords        Coordinates
  !! \param [in]   vtx_ln_to_gn  Local vertex numbering to global vertex numbering
  !!

  subroutine pdm_dist_cloud_surf_surf_mesh_part_set_ (dcs,           &
                                                      i_part,        &
                                                      n_face,        &
                                                      face_vtx_idx,  &
                                                      face_vtx,      &
                                                      face_ln_to_gn, &
                                                      n_vtx,         &
                                                      coords,        &
                                                      vtx_ln_to_gn)
    use iso_c_binding

    implicit none

    type (c_ptr)                       :: dcs
    integer, intent(in)                :: i_part
    integer, intent(in)                :: n_face
    integer,                   pointer :: face_vtx_idx(:)
    integer,                   pointer :: face_vtx(:)
    integer(kind=pdm_g_num_s), pointer :: face_ln_to_gn(:)
    integer, intent(in)                :: n_vtx
    double precision,          pointer :: coords(:,:)
    integer(kind=pdm_g_num_s), pointer :: vtx_ln_to_gn(:)

    integer(c_int)     :: c_i_part
    integer(c_int)     :: c_n_face
    type(c_ptr)        :: c_face_vtx_idx
    type(c_ptr)        :: c_face_vtx
    type(c_ptr)        :: c_face_ln_to_gn
    integer(c_int)     :: c_n_vtx
    type(c_ptr)        :: c_coords
    type(c_ptr)        :: c_vtx_ln_to_gn


    c_i_part = i_part
    c_n_face = n_face
    c_n_vtx  = n_vtx

    c_face_vtx_idx  = c_loc(face_vtx_idx)
    c_face_vtx      = c_loc(face_vtx)
    c_face_ln_to_gn = c_loc(face_ln_to_gn)
    c_coords        = c_loc(coords)
    c_vtx_ln_to_gn  = c_loc(vtx_ln_to_gn)

    call pdm_dist_cloud_surf_surf_mesh_part_set_cf(dcs,             &
                                                   c_i_part,        &
                                                   c_n_face,        &
                                                   c_face_vtx_idx,  &
                                                   c_face_vtx,      &
                                                   c_face_ln_to_gn, &
                                                   c_n_vtx,         &
                                                   c_coords,        &
                                                   c_vtx_ln_to_gn)

  end subroutine pdm_dist_cloud_surf_surf_mesh_part_set_

end module pdm_dist_cloud_surf
