!-----------------------------------------------------------------------------
! This file is part of the ParaDiGM library.
!
! Copyright (C) 2023  ONERA
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

module pdm_vtk

  use pdm
  implicit none

  interface PDM_vtk_write_polydata ; module procedure &
  PDM_vtk_write_polydata_
  end interface

  interface PDM_vtk_write_point_cloud ; module procedure &
  PDM_vtk_write_point_cloud_
  end interface

  private :: PDM_vtk_write_polydata_
  private :: PDM_vtk_write_point_cloud_

  interface

    !!>
    !!
    !! \brief PDM_vtk_write_polydata

    subroutine PDM_vtk_write_polydata_cf (filename,     &
                                          l_filename,   &
                                          n_vtx,        &
                                          vtx_coord,    &
                                          vtx_g_num,    &
                                          n_face,       &
                                          face_vtx_idx, &
                                          face_vtx,     &
                                          face_g_num,   &
                                          face_color)   &

    bind (c, name = 'PDM_vtk_write_polydata_cf')

      use iso_c_binding

      implicit none

      character(kind = c_char, len = 1) :: filename

      integer(c_int), value :: l_filename
      integer(c_int), value :: n_vtx
      integer(c_int), value :: n_face

      type(c_ptr), value    :: vtx_coord
      type(c_ptr), value    :: vtx_g_num
      type(c_ptr), value    :: face_vtx_idx
      type(c_ptr), value    :: face_vtx
      type(c_ptr), value    :: face_g_num
      type(c_ptr), value    :: face_color

    end subroutine PDM_vtk_write_polydata_cf

    !>
    !!
    !! \brief PDM_vtk_write_point_cloud

    subroutine PDM_vtk_write_point_cloud_cf (filename,   &
                                             l_filename, &
                                             n_vtx,      &
                                             vtx_coord,  &
                                             vtx_g_num,  &
                                             color)      &

    bind (c, name = 'PDM_vtk_write_point_cloud_cf')

      use iso_c_binding

      implicit none

      character(kind = c_char, len = 1) :: filename

      integer(c_int), value :: l_filename
      integer(c_int), value :: n_vtx

      type(c_ptr), value    :: vtx_coord
      type(c_ptr), value    :: vtx_g_num
      type(c_ptr), value    :: color

    end subroutine PDM_vtk_write_point_cloud_cf

  end interface

  contains

    !!>
    !!
    !! \brief Write polygons to a vtk format file
    !!
    !! \param [in]   filename
    !! \param [in]   n_vtx
    !! \param [in]   vtx_coord
    !! \param [in]   vtx_g_num
    !! \param [in]   n_face
    !! \param [in]   face_vtx_idx
    !! \param [in]   face_vtx
    !! \param [in]   face_g_num
    !! \param [in]   face_color
    !!

    subroutine PDM_vtk_write_polydata_ (filename,     &
                                        n_vtx,        &
                                        vtx_coord,    &
                                        vtx_g_num,    &
                                        n_face,       &
                                        face_vtx_idx, &
                                        face_vtx,     &
                                        face_g_num,   &
                                        face_color)

      use iso_c_binding

      implicit none

      character(kind = c_char, len = *) :: filename

      integer(kind = c_int) :: l_filename

      integer, intent(in) :: n_vtx
      integer, intent(in) :: n_face

      double precision, pointer          :: vtx_coord(:,:)
      integer(kind=pdm_g_num_s), pointer :: vtx_g_num(:)
      integer, pointer                   :: face_vtx_idx(:)
      integer, pointer                   :: face_vtx(:)
      integer(kind=pdm_g_num_s), pointer :: face_g_num(:)
      integer, pointer                   :: face_color(:)

      type(c_ptr)      :: c_vtx_coord
      type(c_ptr)      :: c_vtx_g_num
      type(c_ptr)      :: c_face_vtx_idx
      type(c_ptr)      :: c_face_vtx
      type(c_ptr)      :: c_face_g_num
      type(c_ptr)      :: c_face_color

      integer(c_int) :: c_n_vtx
      integer(c_int) :: c_n_face

      l_filename = len(filename)

      c_n_vtx  = n_vtx
      c_n_face = n_face

      c_vtx_coord    = c_loc(vtx_coord)
      c_vtx_g_num    = c_loc(vtx_g_num)
      c_face_vtx_idx = c_loc(face_vtx_idx)
      c_face_vtx     = c_loc(face_vtx)
      c_face_g_num   = c_loc(face_g_num)
      c_face_color   = c_loc(face_color)

      call PDM_vtk_write_polydata_cf(filename,       &
                                     l_filename,     &
                                     c_n_vtx,        &
                                     c_vtx_coord,    &
                                     c_vtx_g_num,    &
                                     c_n_face,       &
                                     c_face_vtx_idx, &
                                     c_face_vtx,     &
                                     c_face_g_num,   &
                                     c_face_color)

    end subroutine PDM_vtk_write_polydata_

    !!>
    !!
    !! \brief Write point cloud to a vtk format file
    !!
    !! \param [in]   filename
    !! \param [in]   n_vtx
    !! \param [in]   vtx_coord
    !! \param [in]   color
    !!

    subroutine PDM_vtk_write_point_cloud_ (filename,   &
                                           n_vtx,      &
                                           vtx_coord,  &
                                           vtx_g_num,  &
                                           color)

      use iso_c_binding

      implicit none

      character(kind = c_char, len = *) :: filename

      integer(kind = c_int) :: l_filename

      integer, intent(in) :: n_vtx

      double precision, pointer          :: vtx_coord(:,:)
      integer(kind=pdm_g_num_s), pointer :: vtx_g_num(:)
      integer, pointer                   :: color(:)

      type(c_ptr)      :: c_vtx_coord
      type(c_ptr)      :: c_vtx_g_num
      type(c_ptr)      :: c_color

      integer(c_int) :: c_n_vtx

      l_filename = len(filename)

      c_n_vtx  = n_vtx

      c_vtx_coord = c_loc(vtx_coord)
      c_vtx_g_num = c_loc(vtx_g_num)
      c_color     = c_loc(color)

      call PDM_vtk_write_point_cloud_cf(filename,       &
                                        l_filename,     &
                                        c_n_vtx,        &
                                        c_vtx_coord,    &
                                        c_vtx_g_num,    &
                                        c_color)

    end subroutine PDM_vtk_write_point_cloud_

end module pdm_vtk
