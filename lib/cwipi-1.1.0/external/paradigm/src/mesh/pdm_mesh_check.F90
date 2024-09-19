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

module pdm_mesh_check

  use pdm

  implicit none

  contains

  !!
  !! \brief Remove unconnected vertices in a mesh connectivity
  !!
  !! If a vertex is not cited in the face->vtx connectivity, the function
  !! removes it from the mesh to ensure contiguity
  !!
  !! /TODO : PDM_mesh_check_unconnected_vertex complexity is n^2. This function must be optimized
  !!
  !! \param [in, out] n_vtx       Number of vertices
  !! \param [in, out] l_face_vtx  Size of face->vtx connectivity
  !! \param [in, out] face_vtx    Face->vtx connectivity
  !! \param [in, out] coords      Vertices coordinates
  !! \param [in, out] n_holes     Number of holes
  !!
  !!

  subroutine PDM_mesh_check_unconnected_vertex (n_vtx,      &
                                                l_face_vtx, &
                                                face_vtx,   &
                                                coords,     &
                                                n_holes)
    use iso_c_binding
    implicit none

    integer(pdm_g_num_s), intent(inout)          :: n_vtx
    integer(pdm_g_num_s), intent(inout)          :: l_face_vtx
    integer(pdm_g_num_s), intent(inout), pointer :: face_vtx(:)
    double precision,     intent(inout), pointer :: coords(:)
    integer(pdm_l_num_s), intent(inout)          :: n_holes

    type(c_ptr)                                  :: c_face_vtx
    type(c_ptr)                                  :: c_coords

    interface
      subroutine PDM_mesh_check_unconnected_vertex_c (n_vtx,      &
                                                      l_face_vtx, &
                                                      face_vtx,   &
                                                      coords,     &
                                                      n_holes)    &
      bind (c, name='PDM_mesh_check_unconnected_vertex')
        use iso_c_binding
        implicit none

#ifdef PDM_LONG_G_NUM
        integer(c_long)    :: n_vtx
        integer(c_long)    :: l_face_vtx
#else
        integer(c_int)     :: n_vtx
        integer(c_int)     :: l_face_vtx
#endif
        type(c_ptr), value :: face_vtx
        type(c_ptr), value :: coords
        integer(c_int)     :: n_holes

      end subroutine PDM_mesh_check_unconnected_vertex_c
    end interface

    c_face_vtx = c_loc(face_vtx)
    c_coords   = c_loc(coords)

    call PDM_mesh_check_unconnected_vertex_c (n_vtx,      &
                                              l_face_vtx, &
                                              c_face_vtx, &
                                              c_coords,   &
                                              n_holes)

  end subroutine PDM_mesh_check_unconnected_vertex

end module pdm_mesh_check
