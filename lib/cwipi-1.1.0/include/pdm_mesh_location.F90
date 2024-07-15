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

module pdm_mesh_location

  use pdm
  use iso_c_binding

  implicit none

  !!
  !! Enum type PDM_mesh_location_method_t
  !!

  integer(c_int), parameter :: PDM_MESH_LOCATION_OCTREE         = 0 ! Use point octree
  integer(c_int), parameter :: PDM_MESH_LOCATION_DBBTREE        = 1 ! Use bounding-box tree
  integer(c_int), parameter :: PDM_MESH_LOCATION_LOCATE_ALL_TGT = 2 ! Locate all target points
  ! integer(c_int), parameter :: PDM_MESH_LOCATION_DOCTREE        = 3 !


  interface PDM_mesh_location_create ; module procedure &
  pdm_mesh_location_create_
  end interface

  interface PDM_mesh_location_cloud_set ; module procedure &
  pdm_mesh_location_cloud_set_
  end interface

  interface PDM_mesh_location_part_set ; module procedure &
  pdm_mesh_location_part_set_
  end interface

  interface PDM_mesh_location_part_set_2d ; module procedure &
  pdm_mesh_location_part_set_2d_
  end interface

  interface PDM_mesh_location_unlocated_get ; module procedure &
  pdm_mesh_location_unlocated_get_
  end interface

  interface PDM_mesh_location_located_get ; module procedure &
  pdm_mesh_location_located_get_
  end interface

  interface PDM_mesh_location_point_location_get ; module procedure &
  pdm_mesh_location_point_location_get_
  end interface

  interface PDM_mesh_location_points_in_elt_get ; module procedure &
  pdm_mesh_location_points_in_elt_get_
  end interface

  interface PDM_mesh_location_cell_vertex_get;
    module procedure PDM_mesh_location_cell_vertex_get_cptr
    module procedure PDM_mesh_location_cell_vertex_get_f
  end interface


  private :: pdm_mesh_location_create_
  private :: pdm_mesh_location_cloud_set_
  private :: pdm_mesh_location_part_set_
  private :: pdm_mesh_location_part_set_2d_
  private :: pdm_mesh_location_located_get_
  private :: pdm_mesh_location_unlocated_get_
  private :: pdm_mesh_location_point_location_get_
  private :: pdm_mesh_location_points_in_elt_get_
  private :: pdm_mesh_location_cell_vertex_get_cptr
  private :: pdm_mesh_location_cell_vertex_get_f

  interface

    !>
    !!
    !! \brief Create a structure to compute the location of point clouds inta a mesh
    !!
    !! \param [in]   mesh_nature    Nature of the mesh
    !! \param [in]   n_point_cloud  Number of point cloud
    !! \param [in]   comm           MPI communicator
    !!
    !! \return     Pointer to \ref PDM_mesh_location object
    !!
    !!

    function PDM_mesh_location_create_cf (n_point_cloud, &
                                          comm,          &
                                          owner ) &
                                          result(mloc) &
      bind (c, name = 'PDM_mesh_location_create')

      use iso_c_binding

      implicit none

      integer(c_int), value :: n_point_cloud
      integer(c_int), value :: comm
      integer(c_int), value :: owner

      type(c_ptr)           :: mloc

    end function PDM_mesh_location_create_cf


    subroutine PDM_mesh_location_n_part_cloud_set (mloc, &
                                                   i_point_cloud, &
                                                   n_part) &
    bind (c, name = 'PDM_mesh_location_n_part_cloud_set')
      ! Set the number of partitions of a point cloud
      use iso_c_binding

      implicit none

      type (c_ptr),   value :: mloc          ! C pointer to PDM_mesh_location_t object
      integer(c_int), value :: i_point_cloud ! Point cloud identifier
      integer(c_int), value :: n_part        ! Number of partitions

    end subroutine PDM_mesh_location_n_part_cloud_set

    !>
    !!
    !! \brief Set a point cloud
    !!
    !! \param [in]   mloc            Pointer to \ref PDM_mesh_location object
    !! \param [in]   i_point_cloud   Index of point cloud
    !! \param [in]   i_part          Index of partition
    !! \param [in]   n_points        Number of points
    !! \param [in]   coords          Point coordinates
    !! \param [in]   gnum            Point global number
    !!
    !!

    subroutine PDM_mesh_location_cloud_set_cf (mloc, &
                                               i_point_cloud, &
                                               i_part, &
                                               n_points, &
                                               coords, &
                                               gnum) &
     bind (c, name = 'PDM_mesh_location_cloud_set')

      use iso_c_binding

      implicit none


      type (c_ptr), value :: mloc
      integer(c_int), value :: i_point_cloud
      integer(c_int), value :: i_part
      integer(c_int), value :: n_points
      type(c_ptr), value    :: coords
      type(c_ptr), value    :: gnum

    end subroutine PDM_mesh_location_cloud_set_cf

    !>
    !!
    !! \brief Get a point cloud
    !!
    !! \param [in]   mloc            Pointer to \ref PDM_mesh_location object
    !! \param [in]   i_point_cloud   Index of point cloud
    !! \param [in]   i_part          Index of partition
    !! \param [out]  n_points        Number of points
    !! \param [out]  coords          Point coordinates
    !! \param [out]  gnum            Point global number
    !!
    !!

    subroutine PDM_mesh_location_cloud_get (mloc, &
                                            i_point_cloud, &
                                            i_part, &
                                            n_points, &
                                            coords, &
                                            gnum) &
     bind (c, name = 'PDM_mesh_location_cloud_get')

      use iso_c_binding

      implicit none


      type (c_ptr), value :: mloc
      integer(c_int), value :: i_point_cloud
      integer(c_int), value :: i_part
      integer(c_int), value :: n_points
      type(c_ptr)           :: coords
      type(c_ptr)           :: gnum

    end subroutine PDM_mesh_location_cloud_get

    !>
    !!
    !! \brief Set the mesh nodal
    !!
    !! \param [in]   mloc           Pointer to \ref PDM_mesh_location object
    !! \param [in]   mesh_nodal_id  Mesh nodal Pointer to \ref PDM_mesh_location object
    !!
    !!

    subroutine PDM_mesh_location_shared_nodal_mesh_set (mloc, &
                                                        mesh_nodal_id) &
     bind (c, name = 'PDM_mesh_location_shared_nodal_mesh_set')

      use iso_c_binding

      implicit none


      type (c_ptr), value :: mloc
      type (c_ptr), value :: mesh_nodal_id

    end subroutine PDM_mesh_location_shared_nodal_mesh_set

    !>
    !!
    !! \brief Set global data of a mesh
    !!
    !! \param [in]   mloc           Pointer to \ref PDM_mesh_location object
    !! \param [in]   n_part         Number of partition
    !!
    !!

    subroutine PDM_mesh_location_mesh_n_part_set (mloc, &
                                                       n_part) &
     bind (c, name = 'PDM_mesh_location_mesh_n_part_set')

      use iso_c_binding

      implicit none


      type (c_ptr), value :: mloc
      integer(c_int), value :: n_part

    end subroutine PDM_mesh_location_mesh_n_part_set

    !>
    !!
    !! \brief get cell vertex connectivity
    !!
    !! \param [in]   mloc                  Pointer to \ref PDM_mesh_location object
    !! \param [in]   i_part                Index of partition of the cloud
    !! \param [out]  cell_vtx_idx          Index in (size = n_elt + 1)
    !! \param [out]  cell_vtx              Cell vertex connectivity
    !!
    !!

    subroutine PDM_mesh_location_cell_vertex_get_cf(mloc, &
                                                    i_part, &
                                                    cell_vtx_idx, &
                                                    cell_vtx) &
     bind (c, name = 'PDM_mesh_location_cell_vertex_get')

      use iso_c_binding

      implicit none


      type (c_ptr), value :: mloc
      integer(c_int), value :: i_part
      type(c_ptr)           :: cell_vtx_idx
      type(c_ptr)           :: cell_vtx

    end subroutine PDM_mesh_location_cell_vertex_get_cf

    !>
    !!
    !! \brief Set a part of a mesh
    !!
    !! \param [in]   mloc          Pointer to \ref PDM_mesh_location object
    !! \param [in]   i_part        Partition to define
    !! \param [in]   n_cell        Number of cells
    !! \param [in]   cell_face_idx Index in the cell -> face connectivity
    !! \param [in]   cell_face     cell -> face connectivity
    !! \param [in]   cell_ln_to_gn Local cell numbering to global cel numbering
    !! \param [in]   n_face        Number of faces
    !! \param [in]   face_vtx_idx  Index in the face -> vertex connectivity
    !! \param [in]   face_vtx      face -> vertex connectivity
    !! \param [in]   face_ln_to_gn Local face numbering to global face numbering
    !! \param [in]   n_vtx         Number of vertices
    !! \param [in]   coords        Coordinates
    !! \param [in]   vtx_ln_to_gn  Local vertex numbering to global vertex numbering
    !!
    !!

    subroutine PDM_mesh_location_part_set_cf (mloc, &
                                              i_part, &
                                              n_cell, &
                                              cell_face_idx, &
                                              cell_face, &
                                              cell_ln_to_gn, &
                                              n_face, &
                                              face_vtx_idx, &
                                              face_vtx, &
                                              face_ln_to_gn, &
                                              n_vtx, &
                                              coords, &
                                              vtx_ln_to_gn) &
     bind (c, name = 'PDM_mesh_location_part_set')

      use iso_c_binding

      implicit none


      type (c_ptr), value :: mloc
      integer(c_int), value :: i_part
      integer(c_int), value :: n_cell
      type(c_ptr), value    :: cell_face_idx
      type(c_ptr), value    :: cell_face
      type(c_ptr), value    :: cell_ln_to_gn
      integer(c_int), value :: n_face
      type(c_ptr), value    :: face_vtx_idx
      type(c_ptr), value    :: face_vtx
      type(c_ptr), value    :: face_ln_to_gn
      integer(c_int), value :: n_vtx
      type(c_ptr), value    :: coords
      type(c_ptr), value    :: vtx_ln_to_gn

    end subroutine PDM_mesh_location_part_set_cf


    subroutine PDM_mesh_location_nodal_part_set_cf(mloc, &
                                                   i_part, &
                                                   n_cell, &
                                                   cell_vtx_idx, &
                                                   cell_vtx, &
                                                   cell_ln_to_gn, &
                                                   n_vtx, &
                                                   coords, &
                                                   vtx_ln_to_gn) &
     bind (c, name = 'PDM_mesh_location_nodal_part_set')

      use iso_c_binding

      implicit none

      type (c_ptr), value   :: mloc
      integer(c_int), value :: i_part
      integer(c_int), value :: n_cell
      type(c_ptr), value    :: cell_vtx_idx
      type(c_ptr), value    :: cell_vtx
      type(c_ptr), value    :: cell_ln_to_gn
      integer(c_int), value :: n_vtx
      type(c_ptr), value    :: coords
      type(c_ptr), value    :: vtx_ln_to_gn

    end subroutine PDM_mesh_location_nodal_part_set_cf

    !>
    !!
    !! \brief Set a part of a mesh (2d version)
    !!
    !! \param [in]   mloc          Pointer to \ref PDM_mesh_location object
    !! \param [in]   i_part        Partition to define
    !! \param [in]   n_face        Number of faces
    !! \param [in]   face_edge_idx Index in the face -> edge connectivity
    !! \param [in]   face_edge     face -> edge connectivity
    !! \param [in]   face_ln_to_gn Local face numbering to global cel numbering
    !! \param [in]   n_edge        Number of edges
    !! \param [in]   edge_vtx      edge -> vertex connectivity
    !! \param [in]   n_vtx         Number of vertices
    !! \param [in]   coords        Coordinates
    !! \param [in]   vtx_ln_to_gn  Local vertex numbering to global vertex numbering
    !!
    !!

    subroutine PDM_mesh_location_part_set_2d_cf (mloc, &
                                                 i_part, &
                                                 n_face, &
                                                 face_edge_idx, &
                                                 face_edge, &
                                                 face_ln_to_gn, &
                                                 n_edge, &
                                                 edge_vtx, &
                                                 n_vtx, &
                                                 coords, &
                                                 vtx_ln_to_gn) &
     bind (c, name = 'PDM_mesh_location_part_set_2d')

      use iso_c_binding

      implicit none


      type (c_ptr), value  :: mloc
      integer(c_int), value :: i_part
      integer(c_int), value :: n_face
      type(c_ptr), value    :: face_edge_idx
      type(c_ptr), value    :: face_edge
      type(c_ptr), value    :: face_ln_to_gn
      integer(c_int), value :: n_edge
      type(c_ptr), value    :: edge_vtx
      integer(c_int), value :: n_vtx
      type(c_ptr), value    :: coords
      type(c_ptr), value    :: vtx_ln_to_gn

    end subroutine PDM_mesh_location_part_set_2d_cf

    subroutine PDM_mesh_location_nodal_part_set_2d_cf(mloc, &
                                                      i_part, &
                                                      n_face, &
                                                      face_vtx_idx, &
                                                      face_vtx, &
                                                      face_ln_to_gn, &
                                                      n_vtx, &
                                                      coords, &
                                                      vtx_ln_to_gn) &
     bind (c, name = 'PDM_mesh_location_nodal_part_set_2d')

      use iso_c_binding

      implicit none

      type (c_ptr), value   :: mloc
      integer(c_int), value :: i_part
      integer(c_int), value :: n_face
      type(c_ptr), value    :: face_vtx_idx
      type(c_ptr), value    :: face_vtx
      type(c_ptr), value    :: face_ln_to_gn
      integer(c_int), value :: n_vtx
      type(c_ptr), value    :: coords
      type(c_ptr), value    :: vtx_ln_to_gn

    end subroutine PDM_mesh_location_nodal_part_set_2d_cf

    !>
    !!
    !! \brief Set the tolerance for bounding boxes
    !!
    !! \param [in]   mloc            Pointer to \ref PDM_mesh_location object
    !! \param [in]   tol             Tolerance
    !!
    !!

    subroutine PDM_mesh_location_tolerance_set (mloc, &
                                                tol) &
     bind (c, name = 'PDM_mesh_location_tolerance_set')

      use iso_c_binding

      implicit none


      type (c_ptr), value :: mloc
      real(c_double), value :: tol

    end subroutine PDM_mesh_location_tolerance_set

    !>
    !!
    !! \brief Set the method for computing location
    !!
    !! \param [in]   mloc            Pointer to \ref PDM_mesh_location object
    !! \param [in]   method          Method
    !!
    !!

    subroutine PDM_mesh_location_method_set (mloc, &
                                             method) &
    bind (c, name = 'PDM_mesh_location_method_set')
      ! Set the method for computing location (preconditioning stage)
      use iso_c_binding

      implicit none


      type (c_ptr),   value :: mloc   ! C pointer to PDM_mesh_location_t object
      integer(c_int), value :: method ! Preconditioning method

    end subroutine PDM_mesh_location_method_set

    !>
    !!
    !! \brief Compute point location
    !!
    !! \param [in]   mlocPointer to \ref PDM_mesh_location object
    !!
    !!

    subroutine PDM_mesh_location_compute (mloc) &
      bind (c, name = 'PDM_mesh_location_compute')

      use iso_c_binding

      implicit none


      type (c_ptr), value :: mloc

    end subroutine PDM_mesh_location_compute


    !>
    !!
    !! \brief Get the number of located points
    !!
    !! \param [in]   mloc            Pointer to \ref PDM_mesh_location object
    !! \param [in]   i_point_cloud   Index of point cloud
    !! \param [in]   i_part          Index of partition
    !!
    !! \return     The number of located points
    !!

    function PDM_mesh_location_n_located_get (mloc, &
                                              i_point_cloud, &
                                              i_part) &
                                              result(n_located) &
      bind (c, name = 'PDM_mesh_location_n_located_get')

      use iso_c_binding

      implicit none


      type (c_ptr), value :: mloc
      integer(c_int), value :: i_point_cloud
      integer(c_int), value :: i_part
      integer(c_int)        :: n_located

    end function PDM_mesh_location_n_located_get


    !>
    !!
    !! \brief Get the number of cells
    !!
    !! \param [in]  mloc     Pointer to \ref PDM_mesh_location object
    !! \param [in]  i_part   Index of partition of the mesh
    !!
    !! \return Number of cells
    !!

    function PDM_mesh_location_n_cell_get (mloc, &
                                           i_part) &
                                           result(n_cell) &
      bind (c, name = 'PDM_mesh_location_n_cell_get')

      use iso_c_binding

      implicit none


      type (c_ptr), value :: mloc
      integer(c_int), value :: i_part
      integer(c_int)        :: n_cell

    end function PDM_mesh_location_n_cell_get

    !>
    !!
    !! \brief Get the number of unlocated points
    !!
    !! \param [in]   mloc            Pointer to \ref PDM_mesh_location object
    !! \param [in]   i_point_cloud   Index of point cloud
    !! \param [in]   i_part          Index of partition
    !!
    !! \return     The number of unlocated points
    !!

    function PDM_mesh_location_n_unlocated_get (mloc, &
                                                i_point_cloud, &
                                                i_part) &
                                                result(n_unlocated) &
      bind (c, name = 'PDM_mesh_location_n_unlocated_get')

      use iso_c_binding

      implicit none


      type (c_ptr), value :: mloc
      integer(c_int), value :: i_point_cloud
      integer(c_int), value :: i_part
      integer(c_int)        :: n_unlocated

    end function PDM_mesh_location_n_unlocated_get


    !>
    !!
    !! \brief Get the list of unlocated points
    !!
    !! \param [in]   mloc            Pointer to \ref PDM_mesh_location object
    !! \param [in]   i_point_cloud   Index of point cloud
    !! \param [in]   i_part          Index of partition
    !!
    !! \return     The list of unlocated points
    !!
    !!

    function PDM_mesh_location_unlocated_get_cf (mloc, &
                                                 i_point_cloud, &
                                                 i_part) &
    result(unlocated) &
      bind (c, name = 'PDM_mesh_location_unlocated_get')

      use iso_c_binding

      implicit none


      type (c_ptr), value :: mloc
      integer(c_int), value :: i_point_cloud
      integer(c_int), value :: i_part
      type(c_ptr)           :: unlocated

    end function PDM_mesh_location_unlocated_get_cf


    !>
    !!
    !! \brief Get the list of located points
    !!
    !! \param [in]   mloc            Pointer to \ref PDM_mesh_location object
    !! \param [in]   i_point_cloud   Index of point cloud
    !! \param [in]   i_part          Index of partition
    !!
    !! \return     The list of located points
    !!
    !!

    function PDM_mesh_location_located_get_cf (mloc, &
                                               i_point_cloud, &
                                               i_part) &
    result(located) &
      bind (c, name = 'PDM_mesh_location_located_get')

      use iso_c_binding

      implicit none


      type (c_ptr), value :: mloc
      integer(c_int), value :: i_point_cloud
      integer(c_int), value :: i_part
      type(c_ptr)           :: located

    end function PDM_mesh_location_located_get_cf

    !>
    !!
    !! \brief Get point location for located points
    !!
    !! \param [in]   mloc                  Pointer to \ref PDM_mesh_location object
    !! \param [in]   i_point_cloud         Current cloud
    !! \param [in]   i_part                Index of partition of the cloud
    !! \param [out]  n_points              Number of points in point cloud
    !! \param [out]  location              The global number of the closest element for located points
    !! \param [out]  dist2                 Distance to the located element
    !! \param [out]  projected_coord       Projection on the located element
    !!
    !!
    !!

    subroutine PDM_mesh_location_point_location_get_cf (mloc, &
                                                        i_point_cloud, &
                                                        i_part, &
                                                        location, &
                                                        dist2, &
                                                        projected_coords) &
     bind (c, name = 'PDM_mesh_location_point_location_get')

      use iso_c_binding

      implicit none


      type (c_ptr), value :: mloc
      integer(c_int), value :: i_point_cloud
      integer(c_int), value :: i_part
      type(c_ptr)           :: location
      type(c_ptr)           :: dist2
      type(c_ptr)           :: projected_coords


    end subroutine PDM_mesh_location_point_location_get_cf

    !>
    !!
    !! \brief Get point list located in elements
    !!
    !! \param [in]   mloc                    Pointer to \ref PDM_mesh_location object
    !! \param [in]   i_point_cloud           Index of cloud
    !! \param [in]   i_part                  Index of partition of the mesh
    !! \param [out]  elt_pts_inside_idx      Points index (size = n_elt + 1)
    !! \param [out]  points_gnum             Points global number
    !! \param [out]  points_coords           Points coordinates
    !! \param [out]  points_uvw              Points parametric coordinates in elements
    !! \param [out]  points_weights_idx      Interpolation weights index (size = elt_pts_inside_idx[n_elt] + 1)
    !! \param [out]  points_weights          Interpolation weights
    !! \param [out]  points_dist2            Distance element-points (dist < 0 if the point is inside)
    !! \param [out]  points_projected_coords Point projection on element if the point is outside
    !!

    subroutine PDM_mesh_location_points_in_elt_get_cf (mloc, &
                                                       i_point_cloud, &
                                                       i_part, &
                                                       elt_pts_inside_idx, &
                                                       points_gnum, &
                                                       points_coords, &
                                                       points_uvw, &
                                                       points_weights_idx, &
                                                       points_weights, &
                                                       points_dist2, &
                                                       points_projected_coords) &

      bind (c, name = 'PDM_mesh_location_points_in_elt_get')

      use iso_c_binding

      implicit none


      type (c_ptr), value   :: mloc
      integer(c_int), value :: i_point_cloud
      integer(c_int), value :: i_part
      type(c_ptr)           :: elt_pts_inside_idx
      type(c_ptr)           :: points_gnum
      type(c_ptr)           :: points_coords
      type(c_ptr)           :: points_uvw
      type(c_ptr)           :: points_weights_idx
      type(c_ptr)           :: points_weights
      type(c_ptr)           :: points_dist2
      type(c_ptr)           :: points_projected_coords

    end subroutine PDM_mesh_location_points_in_elt_get_cf


    !>
    !!
    !! \brief Free a mesh location structure
    !!
    !! \param [in]  mloc     Pointer to \ref PDM_mesh_location object
    !! \param [in]  partial  if partial is equal to 0, all data are removed.
    !!                       Otherwise, results are kept.
    !!
    !!

    subroutine PDM_mesh_location_free (mloc)&
     bind (c, name = 'PDM_mesh_location_free')

      use iso_c_binding

      implicit none


      type (c_ptr), value :: mloc

    end subroutine PDM_mesh_location_free


    !>
    !!
    !! \brief  Dump elapsed an CPU time
    !!
    !! \param [in]  mloc     Pointer to \ref PDM_mesh_location object
    !!
    !!

    subroutine PDM_mesh_location_dump_times (mloc) &
      bind (c, name = 'PDM_mesh_location_dump_times')

      use iso_c_binding

      implicit none


      type (c_ptr), value :: mloc

    end subroutine PDM_mesh_location_dump_times


    function PDM_mesh_location_mesh_nodal_get (mloc) &
                                                 result(mesh_nodal) &
      bind (c, name = 'PDM_mesh_location_mesh_nodal_get')

      use iso_c_binding

      implicit none


      type (c_ptr), value :: mloc
      type (c_ptr) :: mesh_nodal

    end function PDM_mesh_location_mesh_nodal_get


    subroutine PDM_mesh_location_part_to_part_get(mloc,   &
                                                  icloud, &
                                                  ptp,    &
                                                  owner)  &
    bind (c, name = 'PDM_mesh_location_part_to_part_get')
      ! Get part_to_part object to exchange data between the source mesh and a target point cloud
      use iso_c_binding

      implicit none


      type (c_ptr),   value :: mloc   ! C pointer to PDM_mesh_location_t object
      integer(c_int), value :: icloud ! Point cloud identifier
      type (c_ptr)          :: ptp    ! Pointer to PDM_part_to_part object
      integer(c_int), value :: owner  ! Ownership for ``ptp``

    end subroutine PDM_mesh_location_part_to_part_get

  end interface


  contains


  subroutine PDM_mesh_location_create_ (mloc,          &
                                        n_point_cloud, &
                                        f_comm,        &
                                        owner)
  ! Create a structure to compute the location of point clouds inside a mesh
  use iso_c_binding

  implicit none

  type(c_ptr), intent(out) :: mloc          ! C pointer to PDM_mesh_location_t object
  integer,     intent(in)  :: n_point_cloud ! Number of point clouds
  integer,     intent(in)  :: f_comm        ! Fortran MPI communicator
  integer,     intent(in)  :: owner         ! Ownership


  integer(c_int) :: c_n_point_cloud
  integer(c_int) :: c_comm
  integer(c_int) :: c_owner

  c_comm = PDM_MPI_Comm_f2c(f_comm)

  c_n_point_cloud = n_point_cloud
  c_owner         = owner

  mloc = PDM_mesh_location_create_cf(c_n_point_cloud, &
                                     c_comm,          &
                                     c_owner)

  end subroutine PDM_mesh_location_create_



  subroutine PDM_mesh_location_cloud_set_ (mloc, &
                                           i_point_cloud, &
                                           i_part, &
                                           n_points, &
                                           coords, &
                                           gnum)
    ! Set a point cloud
    use iso_c_binding

    implicit none

    type (c_ptr), intent(in)           :: mloc          ! C pointer to PDM_mesh_location_t object
    integer, intent(in)                :: i_point_cloud ! Point cloud identifier
    integer, intent(in)                :: i_part        ! Partition identifier
    integer, intent(in)                :: n_points      ! Number of points
    real(8),                   pointer :: coords(:,:)   ! Point coordinates (shape = [3, ``n_points``])
    integer(kind=pdm_g_num_s), pointer :: gnum(:)       ! Point global ids (size = ``n_points``)

    integer(c_int)                     :: c_i_point_cloud
    integer(c_int)                     :: c_i_part
    integer(c_int)                     :: c_n_points
    type(c_ptr)                        :: c_coords
    type(c_ptr)                        :: c_gnum

    c_i_point_cloud = i_point_cloud
    c_i_part        = i_part
    c_n_points      = n_points

    c_coords = C_NULL_PTR
    if (associated(coords)) then
      c_coords = c_loc(coords)
    endif
      
    c_gnum = C_NULL_PTR  
    if (associated(gnum)) then
      c_gnum = c_loc(gnum)
    endif  

    call PDM_mesh_location_cloud_set_cf(mloc,            &
                                        c_i_point_cloud, &
                                        c_i_part,        &
                                        c_n_points,      &
                                        c_coords,        &
                                        c_gnum)

  end subroutine PDM_mesh_location_cloud_set_



  subroutine PDM_mesh_location_part_set_(mloc,          &
                                         i_part,        &
                                         n_cell,        &
                                         cell_face_idx, &
                                         cell_face,     &
                                         cell_ln_to_gn, &
                                         n_face,        &
                                         face_vtx_idx,  &
                                         face_vtx,      &
                                         face_ln_to_gn, &
                                         n_vtx,         &
                                         coords,        &
                                         vtx_ln_to_gn)
    ! Set a *volume* mesh partition
    use iso_c_binding

    implicit none

    type (c_ptr), value                :: mloc             ! C pointer to PDM_mesh_location_t object
    integer, intent(in)                :: i_part           ! Partition identifier
    integer, intent(in)                :: n_cell           ! Number of cells
    integer(kind=pdm_l_num_s), pointer :: cell_face_idx(:) ! Index for cell -> face connectivity (size = ``n_cell`` + 1)
    integer(kind=pdm_l_num_s), pointer :: cell_face(:)     ! Cell -> face connectivity (size = ``cell_face_idx(n_cell+1)``)
    integer(kind=pdm_g_num_s), pointer :: cell_ln_to_gn(:) ! Cell global ids (size = ``n_cell``)
    integer, intent(in)                :: n_face           ! Number of faces
    integer(kind=pdm_l_num_s), pointer :: face_vtx_idx(:)  ! Index for face -> vertex connectivity (size = ``n_face`` + 1)
    integer(kind=pdm_l_num_s), pointer :: face_vtx(:)      ! Face -> vertex connectivity (size = ``face_vtx_idx(n_face+1)``)
    integer(kind=pdm_g_num_s), pointer :: face_ln_to_gn(:) ! Face global ids (size = ``n_face``)
    integer, intent(in)                :: n_vtx            ! Number of vertices
    real(8),                   pointer :: coords(:,:)      ! Vertex coordinates (shape = [3, ``n_vtx``])
    integer(kind=pdm_g_num_s), pointer :: vtx_ln_to_gn(:)  ! Vertex global ids (size = ``n_vtx``)

    integer(c_int)                     :: c_i_part
    integer(c_int)                     :: c_n_cell
    type(c_ptr)                        :: c_cell_face_idx
    type(c_ptr)                        :: c_cell_face
    type(c_ptr)                        :: c_cell_ln_to_gn
    integer(c_int)                     :: c_n_face
    type(c_ptr)                        :: c_face_vtx_idx
    type(c_ptr)                        :: c_face_vtx
    type(c_ptr)                        :: c_face_ln_to_gn
    integer(c_int)                     :: c_n_vtx
    type(c_ptr)                        :: c_coords
    type(c_ptr)                        :: c_vtx_ln_to_gn

    c_i_part = i_part
    c_n_cell = n_cell
    c_n_face = n_face
    c_n_vtx  = n_vtx

    c_cell_face_idx = C_NULL_PTR
    if (associated(cell_face_idx)) then
      c_cell_face_idx = c_loc(cell_face_idx)
    endif
      
    c_cell_face = C_NULL_PTR
    if (associated(cell_face)) then
      c_cell_face     = c_loc(cell_face    )
    endif
      
    c_cell_ln_to_gn = C_NULL_PTR
    if (associated(cell_ln_to_gn)) then
      c_cell_ln_to_gn = c_loc(cell_ln_to_gn)
    endif
      
    c_face_vtx_idx = C_NULL_PTR
    if (associated(face_vtx_idx)) then
      c_face_vtx_idx  = c_loc(face_vtx_idx )
    endif
      
    c_face_vtx = C_NULL_PTR
    if (associated(face_vtx)) then
      c_face_vtx      = c_loc(face_vtx     )
    endif
      
    c_face_ln_to_gn = C_NULL_PTR
    if (associated(face_ln_to_gn)) then
      c_face_ln_to_gn = c_loc(face_ln_to_gn)
    endif
      
    c_coords = C_NULL_PTR
    if (associated(coords)) then
      c_coords        = c_loc(coords       )
    endif
      
    c_vtx_ln_to_gn = C_NULL_PTR
    if (associated(vtx_ln_to_gn)) then
      c_vtx_ln_to_gn  = c_loc(vtx_ln_to_gn )
    endif
      

    call PDM_mesh_location_part_set_cf(mloc,            &
                                       c_i_part,        &
                                       c_n_cell,        &
                                       c_cell_face_idx, &
                                       c_cell_face,     &
                                       c_cell_ln_to_gn, &
                                       c_n_face,        &
                                       c_face_vtx_idx,  &
                                       c_face_vtx,      &
                                       c_face_ln_to_gn, &
                                       c_n_vtx,         &
                                       c_coords,        &
                                       c_vtx_ln_to_gn)

  end subroutine PDM_mesh_location_part_set_


  subroutine PDM_mesh_location_nodal_part_set(mloc,          &
                                              i_part,        &
                                              n_cell,        &
                                              cell_vtx_idx,  &
                                              cell_vtx,      &
                                              cell_ln_to_gn, &
                                              n_vtx,         &
                                              coords,        &
                                              vtx_ln_to_gn)
    ! Set a *volume* mesh partition defined by nodal connectivity
    !
    ! The mesh is assumed to contain only standard elements
    ! (tetrahedra, pyramids, prisms, hexahedra).
    use iso_c_binding

    implicit none

    type (c_ptr), value                :: mloc             ! C pointer to PDM_mesh_location_t object
    integer, intent(in)                :: i_part           ! Partition identifier
    integer, intent(in)                :: n_cell           ! Number of cells
    integer(kind=pdm_l_num_s), pointer :: cell_vtx_idx(:)  ! Index for cell -> face connectivity (size = ``n_cell`` + 1)
    integer(kind=pdm_l_num_s), pointer :: cell_vtx(:)      ! Cell -> face connectivity (size = ``cell_face_idx(n_cell+1)``)
    integer(kind=pdm_g_num_s), pointer :: cell_ln_to_gn(:) ! Cell global ids (size = ``n_cell``)
    integer, intent(in)                :: n_vtx            ! Number of vertices
    real(8),                   pointer :: coords(:,:)      ! Vertex coordinates (shape = [3, ``n_vtx``])
    integer(kind=pdm_g_num_s), pointer :: vtx_ln_to_gn(:)  ! Vertex global ids (size = ``n_vtx``)

    integer(c_int)                     :: c_i_part
    integer(c_int)                     :: c_n_cell
    type(c_ptr)                        :: c_cell_vtx_idx
    type(c_ptr)                        :: c_cell_vtx
    type(c_ptr)                        :: c_cell_ln_to_gn
    integer(c_int)                     :: c_n_vtx
    type(c_ptr)                        :: c_coords
    type(c_ptr)                        :: c_vtx_ln_to_gn

    c_i_part = i_part
    c_n_cell = n_cell
    c_n_vtx  = n_vtx

    c_cell_vtx_idx = C_NULL_PTR
    if (associated(cell_vtx_idx)) then
      c_cell_vtx_idx  = c_loc(cell_vtx_idx )
    endif
      
    c_cell_vtx = C_NULL_PTR
    if (associated(cell_vtx)) then
      c_cell_vtx      = c_loc(cell_vtx     )
    endif
      
    c_cell_ln_to_gn = C_NULL_PTR
    if (associated(cell_ln_to_gn)) then
      c_cell_ln_to_gn = c_loc(cell_ln_to_gn)
    endif
      
    c_coords = C_NULL_PTR
    if (associated(coords)) then
      c_coords        = c_loc(coords       )
    endif
      
    c_vtx_ln_to_gn = C_NULL_PTR
    if (associated(vtx_ln_to_gn)) then
      c_vtx_ln_to_gn  = c_loc(vtx_ln_to_gn )
    endif    

    call PDM_mesh_location_nodal_part_set_cf(mloc,            &
                                             c_i_part,        &
                                             c_n_cell,        &
                                             c_cell_vtx_idx,  &
                                             c_cell_vtx,      &
                                             c_cell_ln_to_gn, &
                                             c_n_vtx,         &
                                             c_coords,        &
                                             c_vtx_ln_to_gn)

  end subroutine PDM_mesh_location_nodal_part_set


  subroutine PDM_mesh_location_part_set_2d_(mloc,          &
                                            i_part,        &
                                            n_face,        &
                                            face_edge_idx, &
                                            face_edge,     &
                                            face_ln_to_gn, &
                                            n_edge,        &
                                            edge_vtx,      &
                                            n_vtx,         &
                                            coords,        &
                                            vtx_ln_to_gn)
    ! Set a *surface* mesh partition
    use iso_c_binding

    implicit none

    type (c_ptr), value                :: mloc             ! C pointer to PDM_mesh_location_t object
    integer, intent(in)                :: i_part           ! Partition identifier
    integer, intent(in)                :: n_face           ! Number of faces
    integer(kind=pdm_l_num_s), pointer :: face_edge_idx(:) ! Index for face -> edge connectivity (size = ``n_face`` + 1)
    integer(kind=pdm_l_num_s), pointer :: face_edge(:)     ! Face -> edge connectivity (size = ``face_edge_idx(n_face+1)``)
    integer(kind=pdm_g_num_s), pointer :: face_ln_to_gn(:) ! Face global ids (size = ``n_face``)
    integer, intent(in)                :: n_edge           ! Number of edges
    integer(kind=pdm_l_num_s), pointer :: edge_vtx(:)      ! Edge -> vertex connectivity (size = 2 * ``n_edge``)
    integer, intent(in)                :: n_vtx            ! Number of vertices
    real(8),                   pointer :: coords(:,:)      ! Vertex coordinates (shape = [3, ``n_vtx``])
    integer(kind=pdm_g_num_s), pointer :: vtx_ln_to_gn(:)  ! Vertex global ids (size = ``n_vtx``)


    integer(c_int)                     :: c_i_part
    integer(c_int)                     :: c_n_face
    type(c_ptr)                        :: c_face_edge_idx
    type(c_ptr)                        :: c_face_edge
    type(c_ptr)                        :: c_face_ln_to_gn
    integer(c_int)                     :: c_n_edge
    type(c_ptr)                        :: c_edge_vtx
    integer(c_int)                     :: c_n_vtx
    type(c_ptr)                        :: c_coords
    type(c_ptr)                        :: c_vtx_ln_to_gn

    c_i_part = i_part
    c_n_face = n_face
    c_n_edge = n_edge
    c_n_vtx  = n_vtx

    c_face_edge_idx = C_NULL_PTR
    if (associated (face_edge_idx)) then
      c_face_edge_idx = c_loc(face_edge_idx)
    endif
      
    c_face_edge = C_NULL_PTR
    if (associated (face_edge)) then
      c_face_edge     = c_loc(face_edge)
    endif
      
    c_face_ln_to_gn = C_NULL_PTR
    if (associated (face_ln_to_gn)) then
      c_face_ln_to_gn = c_loc(face_ln_to_gn)
    endif
      
    c_edge_vtx = C_NULL_PTR
    if (associated (edge_vtx)) then
      c_edge_vtx      = c_loc(edge_vtx     )
    endif
      
    c_coords = C_NULL_PTR
    if (associated (coords)) then
      c_coords        = c_loc(coords       )
    endif
      
    c_vtx_ln_to_gn = C_NULL_PTR
    if (associated (vtx_ln_to_gn)) then
      c_vtx_ln_to_gn  = c_loc(vtx_ln_to_gn )
    endif   

    call PDM_mesh_location_part_set_2d_cf(mloc, &
                                          c_i_part, &
                                          c_n_face, &
                                          c_face_edge_idx, &
                                          c_face_edge, &
                                          c_face_ln_to_gn, &
                                          c_n_edge, &
                                          c_edge_vtx, &
                                          c_n_vtx, &
                                          c_coords, &
                                          c_vtx_ln_to_gn)

  end subroutine PDM_mesh_location_part_set_2d_


  subroutine PDM_mesh_location_nodal_part_set_2d(mloc,          &
                                                 i_part,        &
                                                 n_face,        &
                                                 face_vtx_idx,  &
                                                 face_vtx,      &
                                                 face_ln_to_gn, &
                                                 n_vtx,         &
                                                 coords,        &
                                                 vtx_ln_to_gn)
    ! Set a *surface* mesh partition defined by nodal connectivity
    use iso_c_binding

    implicit none

    type (c_ptr), value                :: mloc             ! C pointer to PDM_mesh_location_t object
    integer, intent(in)                :: i_part           ! Partition identifier
    integer, intent(in)                :: n_face           ! Number of faces
    integer(kind=pdm_l_num_s), pointer :: face_vtx_idx(:)  ! Index for face -> vertex connectivity (size = ``n_face`` + 1)
    integer(kind=pdm_l_num_s), pointer :: face_vtx(:)      ! Face -> vertex connectivity (size = ``face_vtx_idx(n_face+1)``)
    integer(kind=pdm_g_num_s), pointer :: face_ln_to_gn(:) ! Face global ids (size = ``n_face``)
    integer, intent(in)                :: n_vtx            ! Number of vertices
    double precision,          pointer :: coords(:,:)      ! Vertex coordinates (shape = [3, ``n_vtx``])
    integer(kind=pdm_g_num_s), pointer :: vtx_ln_to_gn(:)  ! Vertex global ids (size = ``n_vtx``)

    integer(c_int)                     :: c_i_part
    integer(c_int)                     :: c_n_face
    type(c_ptr)                        :: c_face_vtx_idx
    type(c_ptr)                        :: c_face_vtx
    type(c_ptr)                        :: c_face_ln_to_gn
    integer(c_int)                     :: c_n_vtx
    type(c_ptr)                        :: c_coords
    type(c_ptr)                        :: c_vtx_ln_to_gn

    c_i_part = i_part
    c_n_face = n_face
    c_n_vtx  = n_vtx

    c_face_vtx_idx = C_NULL_PTR
    if (associated(face_vtx_idx)) then
      c_face_vtx_idx  = c_loc(face_vtx_idx)
    endif
      
    c_face_vtx = C_NULL_PTR
    if (associated(face_vtx)) then
      c_face_vtx      = c_loc(face_vtx)
    endif
      
    c_face_ln_to_gn = C_NULL_PTR
    if (associated(face_ln_to_gn)) then
      c_face_ln_to_gn = c_loc(face_ln_to_gn)
    endif
      
    c_coords = C_NULL_PTR
    if (associated(coords)) then
      c_coords        = c_loc(coords       )
    endif
      
    c_vtx_ln_to_gn = C_NULL_PTR
    if (associated(vtx_ln_to_gn)) then
      c_vtx_ln_to_gn  = c_loc(vtx_ln_to_gn )
    endif
      

    call PDM_mesh_location_nodal_part_set_2d_cf(mloc, &
                                                c_i_part, &
                                                c_n_face, &
                                                c_face_vtx_idx, &
                                                c_face_vtx, &
                                                c_face_ln_to_gn, &
                                                c_n_vtx, &
                                                c_coords, &
                                                c_vtx_ln_to_gn)

  end subroutine PDM_mesh_location_nodal_part_set_2d



  subroutine PDM_mesh_location_located_get_ (mloc,          &
                                             i_point_cloud, &
                                             i_part,        &
                                             located)
    ! Get the list of located points
    use iso_c_binding

    implicit none

    type (c_ptr), value :: mloc          ! C pointer to PDM_mesh_location_t object
    integer, intent(in) :: i_point_cloud ! Point cloud identifier
    integer, intent(in) :: i_part        ! Partition identifier
    integer, pointer    :: located(:)    ! List of located points

    integer(c_int)      :: c_i_point_cloud
    integer(c_int)      :: c_i_part
    type(c_ptr)         :: c_located
    integer(c_int)      :: n_located

    c_i_point_cloud = i_point_cloud
    c_i_part        = i_part

    n_located = PDM_mesh_location_n_located_get(mloc,            &
                                                c_i_point_cloud, &
                                                c_i_part)

    c_located = PDM_mesh_location_located_get_cf(mloc,            &
                                                 c_i_point_cloud, &
                                                 c_i_part)

    call c_f_pointer(c_located,   &
                     located,     &
                     [n_located])

  end subroutine PDM_mesh_location_located_get_



  subroutine PDM_mesh_location_unlocated_get_ (mloc,          &
                                               i_point_cloud, &
                                               i_part,        &
                                               unlocated)
    ! Get the list of unlocated points
    use iso_c_binding

    implicit none

    type (c_ptr), value :: mloc          ! C pointer to PDM_mesh_location_t object
    integer, intent(in) :: i_point_cloud ! Point cloud identifier
    integer, intent(in) :: i_part        ! Partition identifier
    integer, pointer    :: unlocated(:)  ! List of unlocated points

    integer(c_int)      :: c_i_point_cloud
    integer(c_int)      :: c_i_part
    type(c_ptr)         :: c_unlocated
    integer(c_int)      :: n_unlocated

    c_i_point_cloud = i_point_cloud
    c_i_part        = i_part

    n_unlocated = PDM_mesh_location_n_unlocated_get(mloc,            &
                                                    c_i_point_cloud, &
                                                    c_i_part)

    c_unlocated = PDM_mesh_location_unlocated_get_cf(mloc,            &
                                                     c_i_point_cloud, &
                                                     c_i_part)

    call c_f_pointer(c_unlocated,   &
                     unlocated,     &
                     [n_unlocated])

  end subroutine PDM_mesh_location_unlocated_get_



  subroutine PDM_mesh_location_point_location_get_ (mloc, &
                                                    i_point_cloud, &
                                                    i_part, &
                                                    location, &
                                                    dist2, &
                                                    projected_coords)
    ! Get point location
    !
    ! .. note::
    !   The results are related to located points only
    use iso_c_binding

    implicit none


    type (c_ptr), value                :: mloc                  ! C pointer to PDM_mesh_location_t object
    integer, intent(in)                :: i_point_cloud         ! Point cloud identifier
    integer, intent(in)                :: i_part                ! Partition identifier
    integer(kind=pdm_g_num_s), pointer :: location(:)           ! Global id of nearest mesh element for located points (size = *n_located*)
    real(8),                   pointer :: dist2(:)              ! Signed squared distance from nearest element (negative if the point is located inside that element) (size = *n_located*)
    real(8),                   pointer :: projected_coords(:,:) ! Cartesian coordinates of projection onto the nearest element (identity if the point is located inside that element)  (shape = [3, *n_located*])

    integer(c_int)                     :: c_i_point_cloud
    integer(c_int)                     :: c_i_part
    type(c_ptr)                        :: c_location
    type(c_ptr)                        :: c_dist2
    type(c_ptr)                        :: c_projected_coords
    integer(c_int)                     :: n_located

    c_i_point_cloud = i_point_cloud
    c_i_part        = i_part

    n_located = PDM_mesh_location_n_located_get(mloc,            &
                                                c_i_point_cloud, &
                                                c_i_part)

    c_location         = C_NULL_PTR
    c_dist2            = C_NULL_PTR
    c_projected_coords = C_NULL_PTR

    call PDM_mesh_location_point_location_get_cf(mloc, &
                                                 c_i_point_cloud, &
                                                 c_i_part, &
                                                 c_location, &
                                                 c_dist2, &
                                                 c_projected_coords)

    call c_f_pointer(c_location, &
                     location,   &
                     [n_located])

    call c_f_pointer(c_dist2, &
                     dist2,   &
                     [n_located])

    call c_f_pointer(c_projected_coords, &
                     projected_coords,   &
                     [3,n_located])

  end subroutine PDM_mesh_location_point_location_get_



  subroutine PDM_mesh_location_points_in_elt_get_ (mloc, &
                                                   i_point_cloud, &
                                                   i_part, &
                                                   elt_pts_inside_idx, &
                                                   points_gnum, &
                                                   points_coords, &
                                                   points_uvw, &
                                                   points_weights_idx, &
                                                   points_weights, &
                                                   points_dist2, &
                                                   points_projected_coords)
    ! Get location data for points located in elements
    use iso_c_binding

    implicit none

    type (c_ptr), value                :: mloc                         ! C pointer to PDM_mesh_location_t object
    integer, intent(in)                :: i_point_cloud                ! Point cloud identifier
    integer, intent(in)                :: i_part                       ! Partition identifier
    integer(kind=pdm_l_num_s), pointer :: elt_pts_inside_idx(:)        ! Index for element -> points mapping (size = *n_elt* + 1)
    integer(kind=pdm_g_num_s), pointer :: points_gnum(:)               ! Located points global ids (size = ``elt_pts_inside_idx(n_elt+1)``)
    real(8),                   pointer :: points_coords(:,:)           ! Located points cartesian coordinates (shape = [3, ``elt_pts_inside_idx(n_elt+1)``])
    real(8),                   pointer :: points_uvw(:,:)              ! Located points parametric coordinates (shape = [3, ``elt_pts_inside_idx(n_elt+1)``])
    integer(kind=pdm_l_num_s), pointer :: points_weights_idx(:)        ! Index for interpolation weights (size = ``elt_pts_inside_idx(n_elt+1)`` + 1)
    real(8),                   pointer :: points_weights(:)            ! Interpolation weights (size = ``points_weights_idx(elt_pts_inside_idx(n_elt+1)+1)``)
    real(8),                   pointer :: points_dist2(:)              ! Signed squared distance element-points (< 0 if the point is inside) (size = ``elt_pts_inside_idx(n_elt+1)``)
    real(8),                   pointer :: points_projected_coords(:,:) ! Cartesian coordinates of projection on element (identity if the point is inside) (shape = [3, ``elt_pts_inside_idx(n_elt+1)``])

    integer(c_int)                     :: c_i_part
    integer(c_int)                     :: c_i_point_cloud
    type(c_ptr)                        :: c_elt_pts_inside_idx
    type(c_ptr)                        :: c_points_gnum
    type(c_ptr)                        :: c_points_coords
    type(c_ptr)                        :: c_points_uvw
    type(c_ptr)                        :: c_points_weights_idx
    type(c_ptr)                        :: c_points_weights
    type(c_ptr)                        :: c_points_dist2
    type(c_ptr)                        :: c_points_projected_coords
    integer(c_int)                     :: n_elt
    integer                            :: n_pts_t


    c_i_point_cloud = i_point_cloud
    c_i_part        = i_part

    n_elt = pdm_mesh_location_n_cell_get(mloc,     &
                                         c_i_part)

    c_elt_pts_inside_idx      = C_NULL_PTR
    c_points_gnum             = C_NULL_PTR
    c_points_coords           = C_NULL_PTR
    c_points_uvw              = C_NULL_PTR
    c_points_weights_idx      = C_NULL_PTR
    c_points_weights          = C_NULL_PTR
    c_points_dist2            = C_NULL_PTR
    c_points_projected_coords = C_NULL_PTR
    
    call PDM_mesh_location_points_in_elt_get_cf(mloc, &
                                                c_i_point_cloud, &
                                                c_i_part, &
                                                c_elt_pts_inside_idx, &
                                                c_points_gnum, &
                                                c_points_coords, &
                                                c_points_uvw, &
                                                c_points_weights_idx, &
                                                c_points_weights, &
                                                c_points_dist2, &
                                                c_points_projected_coords)

    call c_f_pointer(c_elt_pts_inside_idx, &
                     elt_pts_inside_idx,   &
                     [n_elt + 1])

    n_pts_t = elt_pts_inside_idx(n_elt+1)

    call c_f_pointer(c_points_gnum, &
                     points_gnum,   &
                     [n_pts_t])

    call c_f_pointer(c_points_coords, &
                     points_coords,   &
                     [3,n_pts_t])

    call c_f_pointer(c_points_uvw, &
                     points_uvw,   &
                     [3,n_pts_t])

    call c_f_pointer(c_points_weights_idx, &
                     points_weights_idx,   &
                     [n_pts_t+1])

    call c_f_pointer(c_points_weights, &
                     points_weights,   &
                     [points_weights_idx(n_pts_t+1)])

    call c_f_pointer(c_points_gnum, &
                     points_gnum,   &
                     [elt_pts_inside_idx(n_elt+1)])

    call c_f_pointer(c_points_dist2, &
                     points_dist2,   &
                     [n_pts_t])

    call c_f_pointer(c_points_projected_coords, &
                     points_projected_coords,   &
                     [3,n_pts_t])

  end subroutine PDM_mesh_location_points_in_elt_get_


  subroutine PDM_mesh_location_cell_vertex_get_cptr(mloc,         &
                                                    i_part,       &
                                                    cell_vtx_idx, &
                                                    cell_vtx)

    use iso_c_binding

    implicit none


    type (c_ptr),        value :: mloc
    integer(c_int), intent(in) :: i_part
    type(c_ptr)                :: cell_vtx_idx
    type(c_ptr)                :: cell_vtx

    call PDM_mesh_location_cell_vertex_get_cf(mloc,         &
                                              i_part,       &
                                              cell_vtx_idx, &
                                              cell_vtx)

  end subroutine PDM_mesh_location_cell_vertex_get_cptr


  subroutine PDM_mesh_location_cell_vertex_get_f(mloc,         &
                                                 i_part,       &
                                                 cell_vtx_idx, &
                                                 cell_vtx)
    ! Get the cell→vertex connectivity used for internal computations
    !
    ! .. note::
    !   For non-standard elements, this connectivity is built by ParaDiGM and is necessary to associate
    !   the `points_weights` array (returned by \ref PDM_mesh_location_points_in_elt_get)
    !   to the appropriate mesh vertices.
    use iso_c_binding

    implicit none


    type (c_ptr),           value :: mloc            ! C pointer to PDM_mesh_location_t object
    integer(c_int),    intent(in) :: i_part          ! Partition identifier
    integer(pdm_l_num_s), pointer :: cell_vtx_idx(:) ! Index for cell -> vertex connectivity
    integer(pdm_l_num_s), pointer :: cell_vtx(:)     ! Cell -> vertex connectivity
    type(c_ptr)                   :: c_cell_vtx_idx
    type(c_ptr)                   :: c_cell_vtx
    integer                       :: n_cell

    c_cell_vtx_idx = C_NULL_PTR
    c_cell_vtx     = C_NULL_PTR

    call PDM_mesh_location_cell_vertex_get_cf(mloc,           &
                                              i_part,         &
                                              c_cell_vtx_idx, &
                                              c_cell_vtx)

    n_cell = pdm_mesh_location_n_cell_get(mloc, i_part)

    call c_f_pointer(c_cell_vtx_idx, cell_vtx_idx, [n_cell+1])
    call c_f_pointer(c_cell_vtx,     cell_vtx,     [cell_vtx_idx(n_cell+1)])

  end subroutine PDM_mesh_location_cell_vertex_get_f

end module pdm_mesh_location
