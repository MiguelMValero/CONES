!-----------------------------------------------------------------------------
! This file is part of the CWIPI library. 
!
! Copyright (C) 2011-2018  ONERA
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

subroutine callforthointerpfct(entities_dim, &
                               order, &
                               n_local_vertex, &
                               n_local_element, &
                               n_local_polhyedra, &
                               n_distant_point, &
                               local_coordinates, &
                               local_connectivity_index, &
                               local_connectivity, &
                               local_poly_face_index, &
                               local_poly_cell2face_connec, &
                               local_poly_face_connec_idx, &
                               local_poly_face_connec, &
                               dist_pts_coord, &
                               dist_pts_location, &
                               dist_pts_distance, &
                               dist_pts_barycentric_coord_idx, &
                               dist_pts_barycentric_coord, &
                               dist_uvw, &
                               stride, &
                               solver_type, &
                               local_field, &
                               distant_field, &
                               ptInterpolationFct )

  implicit none

  interface
    subroutine  ptInterpolationFct(entities_dim, &
                                   order, &
                                    n_local_vertex, &
                                    n_local_element, &
                                    n_local_polhyedra, &
                                    n_distant_point, &
                                    local_coordinates, &
                                    local_connectivity_index, &
                                    local_connectivity, &
                                    local_poly_face_index, &
                                    local_poly_cell2face_connec, &
                                    local_poly_face_connec_idx, &
                                    local_poly_face_connec, &
                                    dist_pts_coord, &
                                    dist_pts_location, &
                                    dist_pts_distance, &
                                    dist_pts_barycentric_coord_idx, &
                                    dist_pts_barycentric_coord, &
                                    dist_uvw, &
                                    stride, &
                                    solver_type, &
                                    local_field, &
                                    distant_field)
       integer :: entities_dim
       integer :: order
       integer :: n_local_vertex
       integer :: n_local_element
       integer :: n_local_polhyedra
       integer :: n_distant_point
       double precision, dimension(*) :: local_coordinates
       integer, dimension(*) :: local_connectivity_index
       integer, dimension(*) :: local_connectivity
       integer, dimension(*) :: local_poly_face_index
       integer, dimension(*) :: local_poly_cell2face_connec
       integer, dimension(*) :: local_poly_face_connec_idx
       integer, dimension(*) :: local_poly_face_connec
       double precision, dimension(*) :: dist_pts_coord
       integer, dimension(*) :: dist_pts_location
       real(kind=4), dimension(*)  :: dist_pts_distance
       integer, dimension(*) :: dist_pts_barycentric_coord_idx
       double precision, dimension(*) :: dist_pts_barycentric_coord
       double precision, dimension(*) :: dist_uvw
       integer :: stride
       integer :: solver_type
       double precision, dimension(*) :: local_field
       double precision, dimension(*) :: distant_field
     end subroutine ptInterpolationFct
  end interface

  integer :: entities_dim
  integer :: order
  integer :: n_local_vertex
  integer :: n_local_element
  integer :: n_local_polhyedra
  integer :: n_distant_point
  double precision, dimension(*) :: local_coordinates
  integer, dimension(*) :: local_connectivity_index
  integer, dimension(*) :: local_connectivity
  integer, dimension(*) :: local_poly_face_index
  integer, dimension(*) :: local_poly_cell2face_connec
  integer, dimension(*) :: local_poly_face_connec_idx
  integer, dimension(*) :: local_poly_face_connec
  double precision, dimension(*) :: dist_pts_coord
  integer, dimension(*) :: dist_pts_location
  real(kind=4), dimension(*) :: dist_pts_distance
  integer, dimension(*) :: dist_pts_barycentric_coord_idx
  double precision, dimension(*) :: dist_pts_barycentric_coord
  double precision, dimension(*) :: dist_uvw
  integer :: stride
  integer :: solver_type
  double precision, dimension(*) :: local_field
  double precision, dimension(*) :: distant_field


  call ptInterpolationFct (entities_dim, &
                           order, &
                           n_local_vertex, &
                           n_local_element, &
                           n_local_polhyedra, &
                           n_distant_point, &
                           local_coordinates, &
                           local_connectivity_index, &
                           local_connectivity, &
                           local_poly_face_index, &
                           local_poly_cell2face_connec, &
                           local_poly_face_connec_idx, &
                           local_poly_face_connec, &
                           dist_pts_coord, &
                           dist_pts_location, &
                           dist_pts_distance, &
                           dist_pts_barycentric_coord_idx, &
                           dist_pts_barycentric_coord, &
                           dist_uvw, &
                           stride, &
                           solver_type, &
                           local_field, &
                           distant_field)

end subroutine callforthointerpfct

