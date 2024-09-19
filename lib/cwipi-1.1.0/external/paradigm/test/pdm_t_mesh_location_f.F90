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

program testf

  use pdm
#ifdef PDM_HAVE_FORTRAN_MPI_MODULE  
  use mpi
#endif  
  use pdm_mesh_location
  use iso_c_binding

  implicit none

#ifndef PDM_HAVE_FORTRAN_MPI_MODULE  
  include "mpif.h"
#endif  

  !
  ! About point cloud
  !

  integer, parameter :: n_point_cloud = 1
  integer, parameter :: n_part_cloud = 1
  integer, parameter :: i_point_cloud = 0
  integer, parameter :: n_points_into_cloud = 2

  double precision,            pointer :: coords_cloud(:,:) => null()
  integer(kind = pdm_g_num_s), pointer :: gnum_cloud(:)     => null()

  !
  ! About mesh
  !

  integer, parameter :: n_part_mesh = 1

  integer, parameter :: n_cell = 2
  integer(kind = pdm_l_num_s), pointer :: cell_face_idx(:) => null()
  integer(kind = pdm_l_num_s), pointer :: cell_face(:)     => null()
  integer(kind = pdm_g_num_s), pointer :: gnum_cell(:)     => null()


  integer, parameter :: n_face = 11
  integer(kind = pdm_l_num_s), pointer :: face_vtx_idx(:) => null()
  integer(kind = pdm_l_num_s), pointer :: face_vtx(:)     => null()
  integer(kind = pdm_g_num_s), pointer :: gnum_face(:)    => null()


  integer, parameter :: n_vtx = 12
  double precision,            pointer :: coords_vtx(:,:) => null()
  integer(kind = pdm_g_num_s), pointer :: gnum_vtx(:)     => null()

  integer :: i
  integer, parameter :: i_part_cloud = 0
  integer, parameter :: i_part_mesh = 0

  type(c_ptr) :: ml

  ! integer, parameter :: partial = 0 ! Put 1 to keep results when the subroutine closest_points_free is

  integer :: code
  integer :: i_rank
  integer :: n_rank

  !
  ! results
  !

  integer(kind = pdm_g_num_s), pointer :: location(:)           => null()
  double precision,            pointer :: dist2(:)              => null()
  double precision,            pointer :: projected_coords(:,:) => null()

  integer :: n_located
  integer :: n_unlocated
  integer, pointer :: located(:)   => null()
  integer, pointer :: unlocated(:) => null()

  integer(kind = pdm_l_num_s), pointer :: elt_pts_inside_idx(:)        => null()
  integer(kind = pdm_g_num_s), pointer :: points_gnum(:)               => null()
  double precision,            pointer :: points_coords(:,:)           => null()
  double precision,            pointer :: points_uvw(:,:)              => null()
  integer(kind = pdm_l_num_s), pointer :: points_weights_idx(:)        => null()
  double precision,            pointer :: points_weights(:)            => null()
  double precision,            pointer :: points_dist2(:)              => null()
  double precision,            pointer :: points_projected_coords(:,:) => null()

  !
  ! Init
  !

  call mpi_init(code)
  call mpi_comm_rank(mpi_comm_world, i_rank, code)
  call mpi_comm_size(mpi_comm_world, n_rank, code)

  if (n_rank .ne. 1) then
    print *,'Error : 1 MPI process is mandatory'
    call mpi_finalize(code)
    stop
  end if

  !
  ! Set point cloud : 2 points
  !

  allocate(coords_cloud(3,n_points_into_cloud))
  allocate(gnum_cloud(n_points_into_cloud))

  coords_cloud(1,1) = 0.5
  coords_cloud(2,1) = 0.5
  coords_cloud(3,1) = 0.5

  coords_cloud(1,2) = 1.5
  coords_cloud(2,2) = 0.5
  coords_cloud(3,2) = 0.5

  do i = 1, n_points_into_cloud
    gnum_cloud(i) = i
  end do



  !
  ! Set mesh : 2 hexa
  !

  allocate(coords_vtx(3,n_vtx))

  coords_vtx(1,1) = 0.d0
  coords_vtx(2,1) = 0.d0
  coords_vtx(3,1) = 0.d0

  coords_vtx(1,2) = 1.d0
  coords_vtx(2,2) = 0.d0
  coords_vtx(3,2) = 0.d0

  coords_vtx(1,3) = 1.d0
  coords_vtx(2,3) = 1.d0
  coords_vtx(3,3) = 0.d0

  coords_vtx(1,4) = 0.d0
  coords_vtx(2,4) = 1.d0
  coords_vtx(3,4) = 0.d0

  coords_vtx(1,5) = 0.d0
  coords_vtx(2,5) = 0.d0
  coords_vtx(3,5) = 1.d0

  coords_vtx(1,6) = 1.d0
  coords_vtx(2,6) = 0.d0
  coords_vtx(3,6) = 1.d0

  coords_vtx(1,7) = 1.d0
  coords_vtx(2,7) = 1.d0
  coords_vtx(3,7) = 1.d0

  coords_vtx(1,8) = 0.d0
  coords_vtx(2,8) = 1.d0
  coords_vtx(3,8) = 1.d0

  coords_vtx(1,9) = 2.d0
  coords_vtx(2,9) = 0.d0
  coords_vtx(3,9) = 0.d0

  coords_vtx(1,10) = 2.d0
  coords_vtx(2,10) = 1.d0
  coords_vtx(3,10) = 0.d0

  coords_vtx(1,11) = 2.d0
  coords_vtx(2,11) = 0.d0
  coords_vtx(3,11) = 1.d0

  coords_vtx(1,12) = 2.d0
  coords_vtx(2,12) = 1.d0
  coords_vtx(3,12) = 1.d0

  allocate(gnum_vtx(n_vtx))

  do i = 1, n_vtx
    gnum_vtx(i) = i
  end do




  allocate(face_vtx_idx(n_face+1))
  allocate(face_vtx(4*n_face))

  face_vtx_idx(1) = 0
  do i = 1, n_face
    face_vtx_idx(i+1) = face_vtx_idx(i) + 4
  end do

  face_vtx(1) = 1
  face_vtx(2) = 4
  face_vtx(3) = 3
  face_vtx(4) = 2

  face_vtx(5) = 5
  face_vtx(6) = 6
  face_vtx(7) = 7
  face_vtx(8) = 8

  face_vtx(9) = 1
  face_vtx(10) = 5
  face_vtx(11) = 8
  face_vtx(12) = 4

  face_vtx(13) = 6
  face_vtx(14) = 2
  face_vtx(15) = 3
  face_vtx(16) = 7

  face_vtx(17) = 3
  face_vtx(18) = 4
  face_vtx(19) = 8
  face_vtx(20) = 7

  face_vtx(21) = 1
  face_vtx(22) = 2
  face_vtx(23) = 6
  face_vtx(24) = 5

  face_vtx(25) = 11
  face_vtx(26) = 12
  face_vtx(27) = 7
  face_vtx(28) = 6

  face_vtx(29) = 9
  face_vtx(30) = 2
  face_vtx(31) = 3
  face_vtx(32) = 10

  face_vtx(33) = 12
  face_vtx(34) = 10
  face_vtx(35) = 3
  face_vtx(36) = 7

  face_vtx(37) = 11
  face_vtx(38) = 6
  face_vtx(39) = 2
  face_vtx(40) = 9

  face_vtx(41) = 9
  face_vtx(42) = 10
  face_vtx(43) = 12
  face_vtx(44) = 11

  allocate(gnum_face(n_face))

  do i = 1, n_face
    gnum_face(i) = i
  end do




  allocate(cell_face_idx(n_cell+1))
  allocate(cell_face(6*n_cell))

  cell_face_idx(1) = 0
  do i = 1, n_cell
    cell_face_idx(i+1) = cell_face_idx(i) + 6
  end do

  cell_face(1) = 1
  cell_face(2) = 2
  cell_face(3) = 3
  cell_face(4) = 4
  cell_face(5) = 5
  cell_face(6) = 6

  cell_face(7) = 7
  cell_face(8) = 8
  cell_face(9) = 9
  cell_face(10) = 10
  cell_face(11) = 11
  cell_face(12) = -4

  allocate(gnum_cell(n_cell))

  do i = 1, n_cell
    gnum_cell(i) = i
  end do



  !
  ! Create a new PDM_mesh_location structure
  !   The MPI communicator and the number of point cloud are setted
  !

  call PDM_mesh_location_create (ml, &
                                 n_point_cloud, &
                                 MPI_COMM_WORLD, &
                                 PDM_OWNERSHIP_KEEP)

  !
  ! Set the local number partition for any point cloud
  !

  call PDM_mesh_location_n_part_cloud_set (ml, &
                                           i_point_cloud, &
                                           n_part_cloud)

  !
  ! Set point properties for any partition of point cloud
  !

  call PDM_mesh_location_cloud_set (ml, &
                                    i_point_cloud, &
                                    i_part_cloud, &
                                    n_points_into_cloud, &
                                    coords_cloud, &
                                    gnum_cloud)

  !
  ! Set mesh
  !

  call PDM_mesh_location_mesh_n_part_set (ml, &
                                               n_part_mesh)

  call PDM_mesh_location_part_set (ml, &
                                   i_part_mesh, &
                                   n_cell, &
                                   cell_face_idx, &
                                   cell_face, &
                                   gnum_cell, &
                                   n_face, &
                                   face_vtx_idx, &
                                   face_vtx, &
                                   gnum_face, &
                                   n_vtx, &
                                   coords_vtx, &
                                   gnum_vtx)


  !
  ! Set options
  !

  call PDM_mesh_location_tolerance_set (ml, 1.d-4)

  call PDM_mesh_location_method_set (ml, PDM_MESH_LOCATION_OCTREE) ! or PDM_MESH_LOCATION_DBBTREE

  !
  ! Compute
  !

  call PDM_mesh_location_compute (ml)

  !
  ! Get results
  !

  n_unlocated = PDM_mesh_location_n_unlocated_get (ml, &
                                                   i_point_cloud, &
                                                   i_part_cloud)

  n_located = PDM_mesh_location_n_located_get (ml, &
                                               i_point_cloud, &
                                               i_part_cloud)

  call PDM_mesh_location_unlocated_get (ml, &
                                        i_point_cloud, &
                                        i_part_cloud, &
                                        unlocated)

  call PDM_mesh_location_located_get (ml, &
                                      i_point_cloud, &
                                      i_part_cloud, &
                                      located)

  call PDM_mesh_location_point_location_get (ml, &
                                             i_point_cloud, &
                                             i_part_cloud, &
                                             location, &
                                             dist2, &
                                             projected_coords)

  call PDM_mesh_location_points_in_elt_get (ml, &
                                            i_part_mesh, &
                                            i_point_cloud, &
                                            elt_pts_inside_idx, &
                                            points_gnum, &
                                            points_coords, &
                                            points_uvw, &
                                            points_weights_idx, &
                                            points_weights, &
                                            points_dist2, &
                                            points_projected_coords)

  call PDM_mesh_location_dump_times (ml)

  call PDM_mesh_location_free (ml)

  deallocate(coords_cloud)
  deallocate(gnum_cloud)

  deallocate(cell_face_idx)
  deallocate(cell_face)
  deallocate(gnum_cell)

  deallocate(face_vtx_idx)
  deallocate(face_vtx)
  deallocate(gnum_face)

  deallocate(gnum_vtx)
  deallocate(coords_vtx)

  call mpi_finalize(code)


end program testf
