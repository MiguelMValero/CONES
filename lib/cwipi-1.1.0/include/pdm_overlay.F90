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

module pdm_overlay

  use pdm
  use iso_c_binding

  implicit none

  !!
  !! \enum PDM_ol_mesh_t
  !! \brief 3 Meshes
  !!
  !!

  integer(c_int), parameter :: PDM_OL_MESH_A =0  !< First mesh to overlay
  integer(c_int), parameter :: PDM_OL_MESH_B =1  !< Second mesh to overlay

  !!
  !! \enum PDM_ol_parameter_t
  !! \brief Parameters for ovelay meshes building
  !!
  !!

  integer(c_int), parameter :: PDM_OL_CAR_LENGTH_TOL = 0   !< Absolute tolerance for caracteristic length
  integer(c_int), parameter :: PDM_OL_EXTENTS_TOL    = 1   !< Absolute tolerance for extents
  integer(c_int), parameter :: PDM_OL_SAME_PLANE_TOL = 2   !< Absolute tolerance for check if 2 surfaces are the same plane surface

  !!
  !! \enum PDM_ol_mv_t
  !! \brief Type of moving mesh
  !!

  integer(c_int), parameter :: PDM_OL_MV_TRANSFORMATION  = 0 !< Moving with combination of geometric transformations
  integer(c_int), parameter :: PDM_OL_MV_UNKNOWN         = 1 !< Unknown moving type


  interface PDM_ol_create ; module procedure &
  pdm_ol_create_
  end interface

  private :: pdm_ol_create_


  interface

  !>
  !! \brief Build and initialize an overlaying object
  !!
  !! This function builds an initializes an overlaying surface meshes object
  !!
  !! \param [in]  n_partMeshA   Number of local partitions of the meshA input
  !! \param [in]  n_partMeshB   Number of local partitions of the meshB input
  !! \param [in]  projectCoeff  Projection coefficient to define the overlay surface
  !!                            projection
  !!                            If value == 0, the surface projection is MeshA
  !!                            If value == 1, the surface projection is MeshB
  !!                            If 0 < value < 1 , the projection surface is an
  !!                            intermediate surface
  !! \param [in]  fComm         MPI communicator.
  !!
  !! return                     Pointer to overlay object.
  !!

  function PDM_ol_create_cf ( &
    n_partMeshA,           &
    n_partMeshB,           &
    projectCoeff,          &
    comm)                  &
  result (ol)              &
  bind (c, name = 'PDM_ol_create')

  use iso_c_binding
  implicit none

  integer(c_int), value :: n_partMeshA
  integer(c_int), value :: n_partMeshB
  real(c_double), value :: projectCoeff
  integer(c_int), value :: comm

  type(c_ptr)           :: ol

end function PDM_ol_create_cf


!>
!! \brief Set an overlay parameter
!!
!! This function sets an overlay parameter
!!
!! \param [in]  ol       Pointer to Overlay object
!! \param [in]  param    Parameter to define :
!!                            - PDM_OL_CAR_LENGTH_TOL
!!                            - PDM_OL_EXTENTS_TOL
!!                            - PDM_OL_SAME_PLANE_TOL
!! \param [in]  val      Parameter value
!!
!!

subroutine PDM_ol_parameter_set ( &
  ol,                             &
  param,                          &
  val)                            &
bind (c, name = 'PDM_ol_parameter_set')

use iso_c_binding
implicit none

type(c_ptr),    value :: ol
integer(c_int), value :: param
real(c_double), value :: val

end subroutine PDM_ol_parameter_set


!>
!! \brief Define input meshes properties
!!
!! This function defines the input meshes properties
!!
!! \param [in]  ol             Pointer to Overlay object
!! \param [in]  mesh           Input mesh to define (\ref PDM_OL_MESH_A or (\ref PDM_OL_MESH_B)
!! \param [in]  i_part         Partition to define
!! \param [in]  n_face         Number of faces
!! \param [in]  face_vtx_idx   Index in the face -> vertex connectivity (Size : n_face + 1, Numbering : 0 to n-1)
!! \param [in]  face_vtx       face -> vertex connectivity (Size : face_vtx_idx(n_face + 1), Numbering : 1 to n)
!! \param [in]  face_ln_to_gn  Local face numbering to global face numbering
!! \param [in]  n_vtx          Number of vertices
!! \param [in]  coords         Coordinates
!! \param [in]  vtx_ln_to_gn   Local vertex numbering to global vertex numbering
!!
!!

subroutine PDM_ol_input_mesh_set ( &
  ol,                              &
  mesh,                            &
  i_part,                          &
  n_face,                          &
  face_vtx_idx,                    &
  face_vtx,                        &
  face_ln_to_gn,                   &
  n_vtx,                           &
  coords,                          &
  vtx_ln_to_gn)                    &
bind (c, name = 'PDM_ol_input_mesh_set')

use iso_c_binding
implicit none

type   (c_ptr), value :: ol
integer(c_int), value :: mesh
integer(c_int), value :: i_part
integer(c_int), value :: n_face
type   (c_ptr), value :: face_vtx_idx
type   (c_ptr), value :: face_vtx
type   (c_ptr), value :: face_ln_to_gn
integer(c_int), value :: n_vtx
type   (c_ptr), value :: coords
type   (c_ptr), value :: vtx_ln_to_gn

end subroutine PDM_ol_input_mesh_set


!>
!! \brief Define the type of a mesh moving
!!
!! This function defines the type of a mesh moving.
!! Only a mesh can move
!!
!! \param [in]  ol       Pointer to Overlay object
!! \param [in]  mesh     Mesh to move (\ref PDM_OL_MESH_A or (\ref PDM_OL_MESH_B)
!! \param [in]  mv       Type of moving :
!!                         - PDM_OL_MV_TRANSFORMATION
!!                         - PDM_OL_MV_UNKNOWN
!!
!!

subroutine PDM_ol_moving_type_set ( &
  ol,                               &
  mesh,                             &
  mv)                               &
bind (c, name = 'PDM_ol_moving_type_set')

use iso_c_binding
implicit none

type   (c_ptr), value :: ol
integer(c_int), value :: mesh
integer(c_int), value :: mv

end subroutine PDM_ol_moving_type_set


!>
!! \brief Define a translation
!!
!! This function defines a translation for the moving mesh
!!
!! \param [in]  ol       Pointer to Overlay object
!! \param [in]  vect     Translation vector
!! \param [in]  center   Translation center
!!
!!

subroutine PDM_ol_translation_set ( &
  ol,                               &
  vect,                             &
  center)                           &
bind (c, name = 'PDM_ol_translation_set')

use iso_c_binding
implicit none

type(c_ptr), value :: ol
type(c_ptr), value :: vect
type(c_ptr), value :: center

end subroutine PDM_ol_translation_set



!>
!! \brief Define a rotation
!!
!! This function defines a rotation for the moving mesh
!!
!! \param [in]  ol        Pointer to Overlay object
!! \param [in]  direction Rotation direction
!! \param [in]  center    Rotation center
!! \param [in]  angle     Rotation center (degrees)
!!
!!

subroutine PDM_ol_rotation_set ( &
  ol,                            &
  direction,                     &
  center,                        &
  angle)                         &
bind (c, name = 'PDM_ol_rotation_set')

use iso_c_binding
implicit none

type(c_ptr),    value :: ol
type(c_ptr),    value :: direction
type(c_ptr),    value :: center
real(c_double), value :: angle

end subroutine PDM_ol_rotation_set


!>
!! \brief Overlaying the input surface meshes
!!
!! This function overlays the input surface meshes
!!
!! \param [in]  ol        Pointer to Overlay object
!!
!!

subroutine PDM_ol_compute (ol) &
  bind (c, name = 'PDM_ol_compute')

  use iso_c_binding
  implicit none

  type(c_ptr), value :: ol

end subroutine PDM_ol_compute



!>
!! \brief Return the entitie sizes of the overlay mesh
!!
!! This function returns the entities sizes of the overlay mesh
!! for each partition of input meshA or input meshB
!!
!! \param [in]  ol        Pointer to Overlay object
!! \param [in]  mesh      Input mesh
!!                        (\ref PDM_OL_MESH_A or (\ref PDM_OL_MESH_B)
!! \param [out] nGOlFace  Global number of faces of the overlay mesh
!! \param [out] nGOlVtx   Global number of vertices of the overlay mesh
!!
!!

subroutine PDM_ol_mesh_dim_get ( &
  ol,                            &
  mesh,                          &
  nGOlFace,                      &
  nGOlVtx)                       &
bind (c, name = 'PDM_ol_mesh_dim_get')

use iso_c_binding
implicit none

type(c_ptr),    value :: ol
integer(c_int), value :: mesh
#ifdef PDM_LONG_G_NUM
integer(c_long)       :: nGOlFace
integer(c_long)       :: nGOlVtx
#else
integer(c_int)        :: nGOlFace
integer(c_int)        :: nGOlVtx
#endif

end subroutine PDM_ol_mesh_dim_get



!>
!! \brief Return the entitie sizes of the overlay mesh
!!
!! This function returns the entities sizes of the overlay mesh
!! for each partition of input meshA or input meshB
!!
!! \param [in]  ol            Pointer to Overlay object
!! \param [in]  mesh          Input mesh
!!                            (\ref PDM_OL_MESH_A or (\ref PDM_OL_MESH_B)
!! \param [in]  i_part         Partition to define
!! \param [out] nGOlFace      Global number of faces of the overlay mesh
!! \param [out] nGOlVtx       Global number of vertices of the overlay mesh
!! \param [out] nOlFace       Number of faces
!!                            (size = partition number of the \ref mesh)
!! \param [out] nOlVtx        Number of faces
!!                            (size = partition number of the \ref mesh)
!! \param [out] sOlface_vtx    Size of olface_vtx
!! \param [out] sInitToOlFace Size of initToOlFace
!!
!!

subroutine PDM_ol_part_mesh_dim_get ( &
  ol,                                 &
  mesh,                               &
  i_part,                             &
  nOlFace,                            &
  nOlLinkedFace,                      &
  nOlVtx,                             &
  sOlFaceIniVtx,                      &
  sOlface_vtx,                        &
  sInitToOlFace)                      &
bind (c, name = 'PDM_ol_part_mesh_dim_get')

use iso_c_binding
implicit none

type(c_ptr),    value :: ol
integer(c_int), value :: mesh
integer(c_int), value :: i_part
integer(c_int)        :: nOlFace
integer(c_int)        :: nOlLinkedFace
integer(c_int)        :: nOlVtx
integer(c_int)        :: sOlFaceIniVtx
integer(c_int)        :: sOlface_vtx
integer(c_int)        :: sInitToOlFace

end subroutine PDM_ol_part_mesh_dim_get


!>
!! \brief Return the entitie of the overlay mesh
!!
!! This function returns the entities of the overlay mesh
!! for each partition of input meshA or input meshB
!!
!! \param [in]  ol              Pointer to Overlay object
!! \param [in]  mesh            Input mesh
!!                              (\ref PDM_OL_MESH_A or (\ref PDM_OL_MESH_B)
!! \param [in]  i_part           Mesh partition identifier
!! \param [out] olface_vtx_idx    Array adress of \ref olface_vtx index
!!                              (size : \ref nOlFace + 1, numbering : 0 to n-1)
!! \param [out] olface_vtx       Array adress of face vertex connectivity
!!                              (size : \ref sOlface_vtx, numbering : 1 to n)
!! \param [out] olFaceLinkFacedd Array adress of linked face in other mesh
!!                              For each face, 3 link properties :
!!                                    - processus number,
!!                                    - part number,
!!                                    - local face number
!!                              (size : \ref 3 * nOlFace)
!! \param [out] olface_ln_to_gn    Array adress of local to global face numbering
!!                              (size : \ref nOlFace)
!! \param [out] olCoords        Array adress of vertex coodinates
!!                              (size : 3 * \ref nOlVtx)
!! \param [out] olvtx_ln_to_gn     Array adress of local to global vertex
!!                              numbering array (size : \ref nOlVtx)
!! \param [out] initToOlFaceIdx Array adress of \ref initToOlFace index
!!                              (size : \ref nOlVtx + 1)
!! \param [out] initToOlFace    Array adress of initial to overlay faces
!! \param [out] initToOlVtx     Array adress of initial to overlay vertices
!!
!!

subroutine PDM_ol_mesh_entities_get ( &
  ol,                                 &
  mesh,                               &
  i_part,                             &
  olFaceIniVtxIdx,                    &
  olFaceIniVtx,                       &
  olface_vtx_idx,                     &
  olface_vtx,                         &
  olLinkedface_procIdx,               &
  olLinkedFace,                       &
  olface_ln_to_gn,                    &
  olCoords,                           &
  olvtx_ln_to_gn,                     &
  initToOlFaceIdx,                    &
  initToOlFace)                       &
bind (c, name = 'PDM_ol_mesh_entities_get')

use iso_c_binding
implicit none

type(c_ptr),    value :: ol
integer(c_int), value :: mesh
integer(c_int), value :: i_part
type(c_ptr)           :: olFaceIniVtxIdx
type(c_ptr)           :: olFaceIniVtx
type(c_ptr)           :: olface_vtx_idx
type(c_ptr)           :: olface_vtx
type(c_ptr)           :: olLinkedface_procIdx
type(c_ptr)           :: olLinkedFace
type(c_ptr)           :: olface_ln_to_gn
type(c_ptr)           :: olCoords
type(c_ptr)           :: olvtx_ln_to_gn
type(c_ptr)           :: initToOlFaceIdx
type(c_ptr)           :: initToOlFace

end subroutine PDM_ol_mesh_entities_get


!>
!! \brief Delete an overlay object
!!
!! This function deletes an overlay object
!!
!! \param [in]  ol              Pointer to Overlay object
!!
!!

subroutine PDM_ol_del (ol) &
  bind (c, name = 'PDM_ol_del')

  use iso_c_binding
  implicit none

  type(c_ptr), value :: ol

end subroutine PDM_ol_del



!>
!! \brief Dump elapsed an CPU time
!!
!! This function dumps the overlay elapsed time
!!
!! \param [in]  ol              Pointer to Overlay object
!!
!!

subroutine PDM_ol_dump_times (ol) &
  bind (c, name = 'PDM_ol_dump_times')

  use iso_c_binding
  implicit none

  type(c_ptr), value :: ol

end subroutine PDM_ol_dump_times


end interface



contains


!>
  !! \brief Build and initialize an overlaying object
  !!
  !! This function builds an initializes an overlaying surface meshes object
  !!
  !! \param [in]  n_partMeshA   Number of local partitions of the meshA input
  !! \param [in]  nGFaceA       Number of global faces of the meshA input
  !! \param [in]  nGVtxA        Number of global vertices of the meshA input
  !! \param [in]  n_partMeshB   Number of local partitions of the meshB input
  !! \param [in]  nGFaceB       Number of global faces of the meshB input
  !! \param [in]  nGVtxB        Number of global vertices of the meshB input
  !! \param [in]  projectCoeff  Projection coefficient to define the overlay surface
  !!                            projection
  !!                            If value == 0, the surface projection is MeshA
  !!                            If value == 1, the surface projection is MeshB
  !!                            If 0 < value < 1 , the projection surface is an
  !!                            intermediate surface
  !! \param [in]  fComm         MPI communicator.
  !!
  !! return                     Pointer to overlay object.
  !!

  subroutine PDM_ol_create_ ( &
    ol,                       &
    n_partMeshA,              &
    n_partMeshB,              &
    projectCoeff,             &
    f_comm)

  use iso_c_binding
  implicit none

  integer,          intent(in) :: n_partMeshA
  integer,          intent(in) :: n_partMeshB
  double precision, intent(in) :: projectCoeff
  integer,          intent(in) :: f_comm

  type(c_ptr)           :: ol


  integer(c_int) :: c_n_partMeshA
  integer(c_int) :: c_n_partMeshB
  real(c_double) :: c_projectCoeff
  integer(c_int) :: c_comm

  c_comm = PDM_MPI_Comm_f2c(f_comm)

  c_n_partMeshA  = n_partMeshA
  c_n_partMeshB  = n_partMeshB
  c_projectCoeff = projectCoeff

  ol = PDM_ol_create_cf(c_n_partMeshA,  &
                        c_n_partMeshB,  &
                        c_projectCoeff, &
                        c_comm)

  end subroutine PDM_ol_create_

end module pdm_overlay
