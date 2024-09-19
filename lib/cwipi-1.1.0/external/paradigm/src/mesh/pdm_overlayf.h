
!===============================================================================
!  Parameters
!===============================================================================


!*******************************************************************************
! Mesh type
!*******************************************************************************

integer, parameter :: PDM_OL_MESH_A       = 0  ! First mesh to overlay
integer, parameter :: PDM_OL_MESH_B       = 1  ! Second mesh to overlay

!*******************************************************************************
! Parameters for ovelay meshes building
!*******************************************************************************

integer, parameter :: PDM_OL_CAR_LENGTH_TOL = 0  ! Absolute tolerance
                                                 ! for caracteristic length
integer, parameter :: PDM_OL_EXTENTS_TOL    = 1  ! Absolute tolerance
                                                 ! for extents
integer, parameter :: PDM_OL_SAME_PLANE_TOL = 2  ! Absolute tolerance
                                                 ! for check if 2 surfaces are
                                                 ! the same plane surface

!*******************************************************************************
! Type of moving mesh
!*******************************************************************************

integer, parameter :: PDM_OL_MV_TRANSFORMATION = 0  ! Moving with combination of
                                                    ! geometric transformations
integer, parameter :: PDM_OL_MV_UNKNOWN        = 1  ! Unknown moving type


!===============================================================================
!  Interface
!===============================================================================

interface


!*******************************************************************************
! \brief Build and initialize an overlaying object
!
! This function builds an initializes an overlaying surface meshes object
!
! \param [in]  n_partMeshA   Number of local partitions of the meshA input
! \param [in]  nGFaceA      Number of global faces of the meshA input
! \param [in]  nGVtxA       Number of global vertices of the meshA input
! \param [in]  n_partMeshB   Number of local partitions of the meshB input
! \param [in]  nGFaceB      Number of global faces of the meshB input
! \param [in]  nGVtxB       Number of global vertices of the meshB input
! \param [in]  projectCoeff Projection coefficient to define the overlay surface
!                           projection
!                           If value == 0, the surface projection is MeshA
!                           If value == 1, the surface projection is MeshB
!                           If 0 < value < 1 , the projection surface is an
!                           intermediate surface
! \param [in]  comm         MPI communicator.
!
! \param [out] id           Overlay object identifier.
!*******************************************************************************

subroutine PDM_ol_create (n_partMeshA, &
                          nGFaceMeshA,&
                          nGVtxMeshA,&
                          n_partMeshB,&
                          nGFaceMeshB,&
                          nGVtxMeshB,&
                          projectCoeff,&
                          comm,&
                          id)
  implicit none

#include "pdmf.h"

  integer                     :: n_partMeshA
  integer (kind = pdm_long_s) :: nGFaceMeshA
  integer (kind = pdm_long_s) :: nGVtxMeshA
  integer                     :: n_partMeshB
  integer (kind = pdm_long_s) :: nGFaceMeshB
  integer (kind = pdm_long_s) :: nGVtxMeshB
  double precision            :: projectCoeff
  integer                     :: comm
  integer                     :: id

end subroutine PDM_ol_create

!*******************************************************************************
! \brief Set an overlay parameter
!
! This function sets en overlay parameter
!
! \param [in]  id       PDM_ol identifier
! \param [in]  param    Parameter to define :
!                            - PDM_OL_CAR_LENGTH_TOL
!                            - PDM_OL_EXTENTS_TOL
!                            - PDM_OL_SAME_PLANE_TOL
! \param [in]  val      Parameter value
!
!*******************************************************************************

subroutine PDM_ol_parameter_set (id, &
                                 param, &
                                 val)

  implicit none

#include "pdmf.h"

  integer          :: id
  integer          :: param
  double precision :: val

end subroutine PDM_ol_parameter_set

!*******************************************************************************
! \brief Define input meshes properties
!
! This function defines the input meshes properties
!
! \param [in]  id          PDM_ol identifier
! \param [in]  mesh        Input mesh to define
!                          (\ref PDM_OL_MESH_A or (\ref PDM_OL_MESH_B)
! \param [in]  i_part       Partition to define
! \param [in]  n_face       Number of faces
! \param [in]  face_vtx_idx  Index in the face -> vertex connectivity
!                          (Size : n_face + 1, Numbering : 0 to n-1)
! \param [in]  face_vtx     face -> vertex connectivity
!                          (Size : face_vtx_idx(n_face + 1), Numbering : 1 to n)
! \param [in]  face_ln_to_gn  Local face numbering to global face numbering
! \param [in]  n_vtx        Number of vertices
! \param [in]  coords      Coordinates
! \param [in]  vtx_ln_to_gn   Local vertex numbering to global vertex numbering
!
!*******************************************************************************

 subroutine PDM_ol_input_mesh_set (id, &
                                   mesh, &
                                   i_part, &
                                   n_face, &
                                   face_vtx_idx, &
                                   face_vtx, &
                                   face_ln_to_gn, &
                                   n_vtx, &
                                   coords, &
                                   vtx_ln_to_gn)
   implicit none

#include "pdmf.h"

   integer                     :: id
   integer                     :: mesh
   integer                     :: i_part
   integer                     :: n_face
   integer                     :: face_vtx_idx
   integer                     :: face_vtx
   integer (kind = pdm_long_s) :: face_ln_to_gn
   integer                     :: n_vtx
   double precision            :: coords
   integer (kind = pdm_long_s) :: vtx_ln_to_gn

 end subroutine PDM_ol_input_mesh_set

!*******************************************************************************
! \brief Define the type of a mesh moving
!
! This function defines the type of a mesh moving.
! Only a mesh can move
!
! \param [in]  id       PDM_ol identifier
! \param [in]  mesh     Mesh to move (\ref PDM_OL_MESH_A or (\ref PDM_OL_MESH_B)
! \param [in]  mv       Type of moving :
!                         - PDM_OL_MV_TRANSFORMATION
!                         - PDM_OL_MV_UNKNOWN
!
!*******************************************************************************

subroutine PDM_ol_moving_type_set (id, &
                                   mesh, &
                                   mv)

  implicit none

#include "pdmf.h"

  integer :: id
  integer :: mesh
  integer :: mv

end subroutine PDM_ol_moving_type_set

!*******************************************************************************
! \brief Define a translation
!
! This function defines a translation for the moving mesh
!
! \param [in]  id       PDM_overlay identifier
! \param [in]  vect     Translation vector
! \param [in]  center   Translation center
!
!*******************************************************************************

subroutine PDM_ol_translation_set (id, &
                                   vect, &
                                   center)
  implicit none

#include "pdmf.h"

  integer          :: id
  double precision :: vect
  double precision :: center

end subroutine PDM_ol_translation_set

!*******************************************************************************
! \brief Define a rotation
!
! This function defines a rotation for the moving mesh
!
! \param [in]  id        PDM_ol identifier
! \param [in]  direction Rotation direction
! \param [in]  center    Rotation center
! \param [in]  angle     Rotation center (degrees)
!
!*******************************************************************************

subroutine PDM_ol_rotation_set (id, &
                                direction, &
                                center, &
                                angle)

  implicit none

#include "pdmf.h"

  integer          :: id
  double precision :: direction
  double precision :: center
  double precision :: angle

end subroutine PDM_ol_rotation_set

!*******************************************************************************
! \brief Overlaying the input surface meshes
!
! This function overlays the input surface meshes
!
! \param [in]  id       PDM_ol identifier
!
!*******************************************************************************

subroutine PDM_ol_compute (id)

  implicit none

#include "pdmf.h"

  integer :: id

end subroutine PDM_ol_compute

!*******************************************************************************
! \brief Return the entitie sizes of the overlay mesh
!
! This function returns the entities sizes of the overlay mesh
! for each partition of input meshA or input meshB
!
! \param [in]  id        PDM_ol identifier
! \param [in]  mesh      Input mesh
!                        (\ref PDM_OL_MESH_A or (\ref PDM_OL_MESH_B)
! \param [out] nGOlFace  Global number of faces of the overlay mesh
! \param [out] nGOlVtx   Global number of vertices of the overlay mesh
!
!*******************************************************************************

subroutine PDM_ol_mesh_dim_get (id, &
                                mesh, &
                                nGOlFace,&
                                nGOlVtx)

  implicit none

#include "pdmf.h"

  integer                     :: id
  integer                     :: mesh
  integer (kind = pdm_long_s) :: nGOlFace
  integer (kind = pdm_long_s) :: nGOlVtx

end subroutine PDM_ol_mesh_dim_get


!*******************************************************************************
! \brief Return the entitie sizes of the overlay mesh
!
! This function returns the entities sizes of the overlay mesh
! for each partition of input meshA or input meshB
!
! \param [in]  id            PDM_ol identifier
! \param [in]  mesh          Input mesh
!                            (\ref PDM_OL_MESH_A or (\ref PDM_OL_MESH_B)
! \param [in]  i_part         Partition to define
! \param [out] nGOlFace      Global number of faces of the overlay mesh
! \param [out] nGOlVtx       Global number of vertices of the overlay mesh
! \param [out] nOlFace       Number of faces
!                            (size = partition number of the \ref mesh)
! \param [out] nOlVtx        Number of faces
!                            (size = partition number of the \ref mesh)
! \param [out] sOlface_vtx    Size of olface_vtx
! \param [out] sInitToOlFace Size of initToOlFace
!
!*******************************************************************************

subroutine PDM_ol_part_mesh_dim_get (id, &
                                     mesh, &
                                     i_part, &
                                     nOlFace, &
                                     nOlVtx, &
                                     sOlface_vtx,&
                                     sInitToOlFace)

  implicit none

#include "pdmf.h"

  integer :: id
  integer :: mesh
  integer :: i_part
  integer :: nOlFace
  integer :: nOlVtx
  integer :: sOlface_vtx
  integer :: sInitToOlFace

end subroutine PDM_ol_part_mesh_dim_get

!*******************************************************************************
! \brief Return the entitie of the overlay mesh
!
! This function returns the entities of the overlay mesh
! for each partition of input meshA or input meshB
!
! \param [in]  id              PDM_ol identifier
! \param [in]  mesh            Input mesh
!                              (\ref PDM_OL_MESH_A or (\ref PDM_OL_MESH_B)
! \param [in]  i_part           Mesh partition identifier
! \param [out] olface_vtx_idx    Array adress of \ref olface_vtx index
!                              (size : \ref nOlFace + 1, numbering : 0 to n-1)
! \param [out] olface_vtx       Array adress of face vertex connectivity
!                              (size : \ref sOlface_vtx, numbering : 1 to n)
! \param [out] olFaceLinkFacedd Array adress of linked face in other mesh
!                              For each face, 3 link properties :
!                                    - processus number,
!                                    - part number,
!                                    - local face number
!                              (size : \ref 3 * nOlFace)
! \param [out] olface_ln_to_gn    Array adress of local to global face numbering
!                              (size : \ref nOlFace)
! \param [out] olCoords        Array adress of vertex coodinates
!                              (size : 3 * \ref nOlVtx)
! \param [out] olvtx_ln_to_gn     Array adress of local to global vertex
!                              numbering array (size : \ref nOlVtx)
! \param [out] initToOlFaceIdx Array adress of \ref initToOlFace index
!                              (size : \ref nOlVtx + 1)
! \param [out] initToOlFace    Array adress of initial to overlay faces
! \param [out] initToOlVtx     Array adress of initial to overlay vertices
!
!*******************************************************************************

subroutine PDM_ol_mesh_entities_get (id, &
                                     mesh, &
                                     i_part, &
                                     olface_vtx_idx, &
                                     olface_vtx, &
                                     olFaceLinkFace, &
                                     olface_ln_to_gn, &
                                     olCoords, &
                                     olvtx_ln_to_gn, &
                                     initToOlFaceIdx, &
                                     initToOlFace, &
                                     initToOlVtx)

  implicit none

#include "pdmf.h"

  integer                     :: id
  integer                     :: mesh
  integer                     :: i_part
  integer                     :: olface_vtx_idx(*)
  integer                     :: olface_vtx(*)
  integer                     :: olFaceLinkFace(*)
  integer (kind = pdm_long_s) :: olface_ln_to_gn(*)
  double precision            :: olCoords(*)
  integer (kind = pdm_long_s) :: olvtx_ln_to_gn(*)
  integer                     :: initToOlFaceIdx(*)
  integer                     :: initToOlFace(*)
  integer                     :: initToOlVtx(*)

end subroutine PDM_ol_mesh_entities_get

!*******************************************************************************
! \brief Delete an overlay object
!
! This function deletes an overlay object
!
! \param [in]  id                PDM_ol identifier.
!
!*******************************************************************************

subroutine PDM_ol_del (id)

  implicit none

#include "pdmf.h"

  integer ::     id

end subroutine PDM_ol_del

end interface
