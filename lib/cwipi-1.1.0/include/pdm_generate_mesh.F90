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


module pdm_generate_mesh

  use pdm
  use pdm_pointer_array

  implicit none

  interface PDM_generate_mesh_rectangle_simplified ; module procedure &
  PDM_generate_mesh_rectangle_simplified_
  end interface

  ! interface PDM_generate_mesh_sphere_simplified ; module procedure &
  ! PDM_generate_mesh_sphere_simplified_
  ! end interface

  ! interface PDM_generate_mesh_ball_simplified ; module procedure &
  ! PDM_generate_mesh_ball_simplified_
  ! end interface

  ! interface PDM_generate_mesh_parallelepiped_simplified ; module procedure &
  ! PDM_generate_mesh_parallelepiped_simplified_
  ! end interface

  ! interface PDM_generate_mesh_sphere ; module procedure &
  ! PDM_generate_mesh_sphere_
  ! end interface

  ! interface PDM_generate_mesh_rectangle ; module procedure &
  ! PDM_generate_mesh_rectangle_
  ! end interface

  ! interface PDM_generate_mesh_ball ; module procedure &
  ! PDM_generate_mesh_ball_
  ! end interface

  ! interface PDM_generate_mesh_parallelepiped ; module procedure &
  ! PDM_generate_mesh_parallelepiped_
  ! end interface

  ! interface PDM_generate_mesh_rectangle_ngon ; module procedure &
  ! PDM_generate_mesh_rectangle_ngon_
  ! end interface

  ! interface PDM_generate_mesh_sphere_ngon ; module procedure &
  ! PDM_generate_mesh_sphere_ngon_
  ! end interface

  ! interface PDM_generate_mesh_ball_ngon ; module procedure &
  ! PDM_generate_mesh_ball_ngon_
  ! end interface

  ! interface PDM_generate_mesh_parallelepiped_ngon ; module procedure &
  ! PDM_generate_mesh_parallelepiped_ngon_
  ! end interface


  private :: PDM_generate_mesh_rectangle_simplified_
  ! private :: PDM_generate_mesh_sphere_simplified
  ! private :: PDM_generate_mesh_ball_simplified
  ! private :: PDM_generate_mesh_parallelepiped_simplified
  ! private :: PDM_generate_mesh_sphere
  ! private :: PDM_generate_mesh_rectangle
  ! private :: PDM_generate_mesh_ball
  ! private :: PDM_generate_mesh_parallelepiped
  ! private :: PDM_generate_mesh_rectangle_ngon
  ! private :: PDM_generate_mesh_sphere_ngon
  ! private :: PDM_generate_mesh_ball_ngon
  ! private :: PDM_generate_mesh_parallelepiped_ngon

  interface

  subroutine PDM_generate_mesh_rectangle_simplified_cf(comm,        &
                                                       n_vtx_seg,   &
                                                       n_vtx,       &
                                                       n_elt,       &
                                                       coords,      &
                                                       elt_vtx_idx, &
                                                       elt_vtx)     &

      bind (c, name = 'PDM_generate_mesh_rectangle_simplified')

      use iso_c_binding
      implicit none

      integer(c_int),  value :: comm
#ifdef PDM_LONG_G_NUM
      integer(c_long), value :: n_vtx_seg
#else
      integer(c_int),  value :: n_vtx_seg
#endif
      integer(c_int)         :: n_vtx
      integer(c_int)         :: n_elt
      type(c_ptr)            :: coords
      type(c_ptr)            :: elt_vtx_idx
      type(c_ptr)            :: elt_vtx

  end subroutine PDM_generate_mesh_rectangle_simplified_cf

  subroutine PDM_generate_mesh_sphere_simplified_cf(comm,        &
                                                    n_vtx,       &
                                                    n_elt,       &
                                                    coords,      &
                                                    elt_vtx_idx, &
                                                    elt_vtx)     &

    bind (c, name = 'PDM_generate_mesh_sphere_simplified')

    use iso_c_binding
    implicit none

    integer(c_int), value  :: comm 
    integer(c_int)         :: n_vtx
    integer(c_int)         :: n_elt
    type(c_ptr)            :: coords
    type(c_ptr)            :: elt_vtx_idx
    type(c_ptr)            :: elt_vtx
  end subroutine PDM_generate_mesh_sphere_simplified_cf

  subroutine PDM_generate_mesh_ball_simplified_cf(comm,        &
                                                  n_vtx,       &
                                                  n_elt,       &
                                                  coords,      &
                                                  elt_vtx_idx, &
                                                  elt_vtx)     &

    bind (c, name = 'PDM_generate_mesh_ball_simplified')

    use iso_c_binding
    implicit none

    integer(c_int), value :: comm 
    integer(c_int)        :: n_vtx
    integer(c_int)        :: n_elt
    type(c_ptr)           :: coords
    type(c_ptr)           :: elt_vtx_idx
    type(c_ptr)           :: elt_vtx
  end subroutine PDM_generate_mesh_ball_simplified_cf

  subroutine PDM_generate_mesh_parallelepiped_simplified_cf(comm,        &
                                                            n_vtx_seg,   &
                                                            n_vtx,       &
                                                            n_elt,       &
                                                            coords,      &
                                                            elt_vtx_idx, &
                                                            elt_vtx)     &

    bind (c, name = 'PDM_generate_mesh_parallelepiped_simplified')

    use iso_c_binding
    implicit none

    integer(c_int), value :: comm 
#ifdef PDM_LONG_G_NUM
    integer(c_long), value :: n_vtx_seg
#else
    integer(c_int),  value :: n_vtx_seg
#endif
    integer(c_int)         :: n_vtx
    integer(c_int)         :: n_elt
    type(c_ptr)            :: coords
    type(c_ptr)            :: elt_vtx_idx
    type(c_ptr)            :: elt_vtx
  end subroutine PDM_generate_mesh_parallelepiped_simplified_cf

  function PDM_generate_mesh_sphere_cf(comm,        &
                                       elt_type,    &
                                       order,       &
                                       ho_ordering, &
                                       radius,      &
                                       center_x,    &
                                       center_y,    &
                                       center_z,    &
                                       n_u,         &
                                       n_v,         &
                                       n_part,      &
                                       part_method) &
    result(mesh_nodal) &
    bind (c, name = 'PDM_generate_mesh_sphere')

    use iso_c_binding
    implicit none

    integer(c_int), value :: comm 
    integer(c_int), value :: elt_type
    integer(c_int), value :: order
    type(c_ptr), value    :: ho_ordering
    real(c_double), value :: radius
    real(c_double), value :: center_x
    real(c_double), value :: center_y
    real(c_double), value :: center_z
#ifdef PDM_LONG_G_NUM
    integer(c_long), value :: n_u
    integer(c_long), value :: n_v
#else
    integer(c_int), value :: n_u
    integer(c_int), value :: n_v
#endif
    integer(c_int), value :: n_part
    integer(c_int), value :: part_method

    type(c_ptr)            :: mesh_nodal
  end function PDM_generate_mesh_sphere_cf

  function PDM_generate_mesh_rectangle_cf (comm,        &
                                           elt_type,    &
                                           order,       &
                                           ho_ordering, &
                                           xmin,        &
                                           ymin,        &
                                           zmin,        &
                                           lengthx,     &
                                           lengthy,     &
                                           n_x,         &
                                           n_y,         &
                                           n_part,      &
                                           part_method) &

    result(mesh_nodal) &
    bind (c, name = 'PDM_generate_mesh_rectangle')

    use iso_c_binding
    implicit none

    integer(c_int), value :: comm
    integer(c_int), value :: elt_type
    integer(c_int), value :: order
    type(c_ptr), value    :: ho_ordering
    real(c_double), value :: xmin
    real(c_double), value :: ymin
    real(c_double), value :: zmin
    real(c_double), value :: lengthx
    real(c_double), value :: lengthy
#ifdef PDM_LONG_G_NUM
    integer(c_long), value :: n_x
    integer(c_long), value :: n_y
#else
    integer(c_int), value :: n_x
    integer(c_int), value :: n_y
#endif
    integer(c_int), value :: n_part
    integer(c_int), value :: part_method

    type(c_ptr)            :: mesh_nodal

  end function PDM_generate_mesh_rectangle_cf

  function PDM_generate_mesh_ball_cf (comm,            &
                                        elt_type,        &
                                        order,           &
                                        ho_ordering,     &
                                        radius,          &
                                        hole_radius,     &
                                        center_x,        &
                                        center_y,        &
                                        center_z,        &
                                        n_x,             &
                                        n_y,             &
                                        n_z,             &
                                        n_layer,         &
                                        geometric_ratio, &
                                        n_part,          &
                                        part_method)     &

    result(mesh_nodal) &
    bind (c, name = 'PDM_generate_mesh_ball')

    use iso_c_binding
    implicit none

    integer(c_int), value :: comm
    integer(c_int), value :: elt_type
    integer(c_int), value :: order
    type(c_ptr), value    :: ho_ordering
    real(c_double), value :: radius
    real(c_double), value :: hole_radius
    real(c_double), value :: center_x
    real(c_double), value :: center_y
    real(c_double), value :: center_z
#ifdef PDM_LONG_G_NUM
    integer(c_long), value :: n_x
    integer(c_long), value :: n_y
    integer(c_long), value :: n_z
    integer(c_long), value :: n_layer
#else
    integer(c_int), value :: n_x
    integer(c_int), value :: n_y
    integer(c_int), value :: n_z
    integer(c_int), value :: n_layer
#endif
    real(c_double), value :: geometric_ratio
    integer(c_int), value :: n_part
    integer(c_int), value :: part_method

    type(c_ptr)            :: mesh_nodal

  end function PDM_generate_mesh_ball_cf

  ! subroutine PDM_generate_mesh_parallelepiped_cf
  !   bind (c, name = 'PDM_generate_mesh_parallelepiped')

  !   use iso_c_binding
  !   implicit none

  ! end subroutine PDM_generate_mesh_parallelepiped_cf

  subroutine PDM_generate_mesh_rectangle_ngon_cf(comm,           &
                                                 elt_type,       &
                                                 xmin,           &
                                                 ymin,           &
                                                 zmin,           &
                                                 lengthx,        &
                                                 lengthy,        &
                                                 n_x,            &
                                                 n_y,            &
                                                 n_part,         &
                                                 part_method,    &
                                                 random_factor,  &
                                                 pn_vtx,         &
                                                 pn_edge,        &
                                                 pn_face,        &
                                                 pvtx_coord,     &
                                                 pedge_vtx,      &
                                                 pface_edge_idx, &
                                                 pface_edge,     &
                                                 pface_vtx,      &
                                                 pvtx_ln_to_gn,  &
                                                 pedge_ln_to_gn, &
                                                 pface_ln_to_gn) &
    bind (c, name = 'PDM_generate_mesh_rectangle_ngon')

    use iso_c_binding
    implicit none

    integer(c_int), value :: comm
    integer(c_int), value :: elt_type
    real(c_double), value :: xmin, ymin, zmin, lengthx, lengthy
#ifdef PDM_LONG_G_NUM
    integer(c_long), value :: n_x, n_y
#else
    integer(c_int),  value :: n_x, n_y
#endif
    integer(c_int),  value :: n_part, part_method
    real(c_double),  value :: random_factor
    type(c_ptr)            :: pn_vtx
    type(c_ptr)            :: pn_edge
    type(c_ptr)            :: pn_face
    type(c_ptr)            :: pvtx_coord
    type(c_ptr)            :: pedge_vtx
    type(c_ptr)            :: pface_edge_idx
    type(c_ptr)            :: pface_edge
    type(c_ptr)            :: pface_vtx
    type(c_ptr)            :: pvtx_ln_to_gn
    type(c_ptr)            :: pedge_ln_to_gn
    type(c_ptr)            :: pface_ln_to_gn
  end subroutine PDM_generate_mesh_rectangle_ngon_cf

  ! subroutine PDM_generate_mesh_sphere_ngon_cf
  !   bind (c, name = 'PDM_generate_mesh_sphere_ngon')

  !   use iso_c_binding
  !   implicit none

  !   integer(c_int), value :: comm 
  ! end subroutine PDM_generate_mesh_sphere_ngon_cf

  ! subroutine PDM_generate_mesh_ball_ngon_cf
  !   bind (c, name = 'PDM_generate_mesh_ball_ngon')

  !   use iso_c_binding
  !   implicit none

  !   integer(c_int), value :: comm 
  !   integer(c_int)         :: n_vtx
  !   integer(c_int)         :: n_elt
  !   type(c_ptr)            :: coords
  !   type(c_ptr)            :: elt_vtx_idx
  !   type(c_ptr)            :: elt_vtx
  ! end subroutine PDM_generate_mesh_ball_ngon_cf

  ! subroutine PDM_generate_mesh_parallelepiped_ngon_cf
  !   bind (c, name = 'PDM_generate_mesh_parallelepiped_ngon')

  !   use iso_c_binding
  !   implicit none

  !   integer(c_int), value :: comm 
  !   integer(c_int)         :: n_vtx
  !   integer(c_int)         :: n_elt
  !   type(c_ptr)            :: coords
  !   type(c_ptr)            :: elt_vtx_idx
  !   type(c_ptr)            :: elt_vtx
  ! end subroutine PDM_generate_mesh_parallelepiped_ngon_cf


  end interface

  contains

  !>
  !!
  !! \brief Create a simple partitionned rectangle mesh (2D).
  !!
  !! \param [in]   comm        MPI communicator
  !! \param [in]   n_vtx_seg   Number of vertices along each side of the rectangle
  !! \param [out]  n_vtx       Number of vertices
  !! \param [out]  n_elt       Number of elements
  !! \param [out]  coords      Array of vertex coordinates
  !! \param [out]  elt_vtx_idx Index array of the element vertex connectivity
  !! \param [out]  elt_vtx     Array of the element vertex connectivity
  !!
  !!

  subroutine PDM_generate_mesh_rectangle_simplified_(comm,        &
                                                     n_vtx_seg,   &
                                                     n_vtx,       &
                                                     n_elt,       &
                                                     coords,      &
                                                     elt_vtx_idx, &
                                                     elt_vtx)

      use iso_c_binding
      implicit none

      integer,                     intent(in) :: comm
      integer(kind=pdm_g_num_s),   intent(in) :: n_vtx_seg
      integer,                    intent(out) :: n_vtx
      integer,                    intent(out) :: n_elt
      double precision,               pointer :: coords(:,:)
      integer(kind=pdm_l_num_s),      pointer :: elt_vtx_idx(:)
      integer(kind=pdm_l_num_s),      pointer :: elt_vtx(:)

      integer(c_int)                          :: c_comm
      integer(c_int)                          :: c_n_vtx
      integer(c_int)                          :: c_n_elt
      type(c_ptr)                             :: c_coords
      type(c_ptr)                             :: c_elt_vtx_idx
      type(c_ptr)                             :: c_elt_vtx

      c_comm = PDM_MPI_Comm_f2c(comm)

      call PDM_generate_mesh_rectangle_simplified_cf(c_comm,        &
                                                     n_vtx_seg,     &
                                                     c_n_vtx,       &
                                                     c_n_elt,       &
                                                     c_coords,      &
                                                     c_elt_vtx_idx, &
                                                     c_elt_vtx)

      n_vtx = c_n_vtx
      n_elt = c_n_elt

      call c_f_pointer(c_coords, &
                       coords,   &
                       [3, n_vtx])

      call c_f_pointer(c_elt_vtx_idx, &
                       elt_vtx_idx,   &
                       [n_elt + 1])

      call c_f_pointer(c_elt_vtx, &
                       elt_vtx,   &
                       [elt_vtx_idx(n_elt + 1)])

  end subroutine PDM_generate_mesh_rectangle_simplified_



  subroutine PDM_generate_mesh_rectangle_ngon(comm,           &
                                              elt_type,       &
                                              xmin,           &
                                              ymin,           &
                                              zmin,           &
                                              lengthx,        &
                                              lengthy,        &
                                              n_x,            &
                                              n_y,            &
                                              n_part,         &
                                              part_method,    &
                                              pn_vtx,         &
                                              pn_edge,        &
                                              pn_face,        &
                                              pvtx_coord,     &
                                              pedge_vtx,      &
                                              pface_edge_idx, &
                                              pface_edge,     &
                                              pface_vtx,      &
                                              pvtx_ln_to_gn,  &
                                              pedge_ln_to_gn, &
                                              pface_ln_to_gn, &
                                              random_factor_opt)
    ! Create a partitioned rectangular mesh (2D) with descending connectivities
    !
    ! Admissible values for ``elt_type``:
    !   - ``PDM_MESH_NODAL_TRIA3``   : triangles
    !   - ``PDM_MESH_NODAL_QUAD4``   : quadrangles
    !   - ``PDM_MESH_NODAL_POLY_2D`` : mixed polygons (triangles, quadrangles and octagons)
    use iso_c_binding
    implicit none

    integer(c_int),       intent(in)           :: comm              ! MPI communicator
    integer(c_int),       intent(in)           :: elt_type          ! Element type
    real(c_double),       intent(in)           :: xmin              ! Minimal x-coordinate
    real(c_double),       intent(in)           :: ymin              ! Minimal y-coordinate
    real(c_double),       intent(in)           :: zmin              ! Minimal z-coordinate
    real(c_double),       intent(in)           :: lengthx           ! Length of the rectangle in the x-direction
    real(c_double),       intent(in)           :: lengthy           ! Length of the rectangle in the y-direction
    integer(pdm_g_num_s), intent(in)           :: n_x               ! Number of points in the x-direction
    integer(pdm_g_num_s), intent(in)           :: n_y               ! Number of points in the y-direction
    integer(c_int),       intent(in)           :: n_part            ! Number of partitions
    integer(c_int),       intent(in)           :: part_method       ! Partitioning method
    integer(pdm_l_num_s),      pointer         :: pn_vtx(:)         ! Number of vertices
    integer(pdm_l_num_s),      pointer         :: pn_edge(:)        ! Number of edges
    integer(pdm_l_num_s),      pointer         :: pn_face(:)        ! Number of faces
    type(PDM_pointer_array_t), pointer         :: pvtx_coord        ! Vertex coordinates
    type(PDM_pointer_array_t), pointer         :: pedge_vtx         ! Edge->vertex connectivity
    type(PDM_pointer_array_t), pointer         :: pface_edge_idx    ! Index of face->edge connectivity
    type(PDM_pointer_array_t), pointer         :: pface_edge        ! Face->edge connectivity
    type(PDM_pointer_array_t), pointer         :: pface_vtx         ! Face->vertex connectivity
    type(PDM_pointer_array_t), pointer         :: pvtx_ln_to_gn     ! Vertex global ids
    type(PDM_pointer_array_t), pointer         :: pedge_ln_to_gn    ! Edge global ids
    type(PDM_pointer_array_t), pointer         :: pface_ln_to_gn    ! Face global ids
    real(c_double),       intent(in), optional :: random_factor_opt ! Randomization factor (between 0 and 1)

    integer(c_int)          :: c_comm
    type(c_ptr)             :: c_pn_vtx
    type(c_ptr)             :: c_pn_edge
    type(c_ptr)             :: c_pn_face
    type(c_ptr)             :: cc_pvtx_coord
    type(c_ptr)             :: cc_pedge_vtx
    type(c_ptr)             :: cc_pface_edge_idx
    type(c_ptr)             :: cc_pface_edge
    type(c_ptr)             :: cc_pface_vtx
    type(c_ptr)             :: cc_pvtx_ln_to_gn
    type(c_ptr)             :: cc_pedge_ln_to_gn
    type(c_ptr)             :: cc_pface_ln_to_gn

    integer                 :: length(n_part)
    type(c_ptr),    pointer :: fptr(:)
    integer(c_int), pointer :: face_edge_idx(:)
    integer                 :: i_part
    real(c_double)          :: random_factor

    c_comm = PDM_MPI_Comm_f2c(comm)

    c_pn_vtx          = C_NULL_PTR
    c_pn_edge         = C_NULL_PTR
    c_pn_face         = C_NULL_PTR
    cc_pvtx_coord     = C_NULL_PTR
    cc_pedge_vtx      = C_NULL_PTR
    cc_pface_edge_idx = C_NULL_PTR
    cc_pface_edge     = C_NULL_PTR
    cc_pface_vtx      = C_NULL_PTR
    cc_pvtx_ln_to_gn  = C_NULL_PTR
    cc_pedge_ln_to_gn = C_NULL_PTR
    cc_pface_ln_to_gn = C_NULL_PTR
    
    fptr          => null() 
    face_edge_idx => null() 
    
    if (present(random_factor_opt)) then
      random_factor = random_factor_opt
    else
      random_factor = 0.d0
    endif

    call PDM_generate_mesh_rectangle_ngon_cf(c_comm,            &
                                             elt_type,          &
                                             xmin,              &
                                             ymin,              &
                                             zmin,              &
                                             lengthx,           &
                                             lengthy,           &
                                             n_x,               &
                                             n_y,               &
                                             n_part,            &
                                             part_method,       &
                                             random_factor,     &
                                             c_pn_vtx,          &
                                             c_pn_edge,         &
                                             c_pn_face,         &
                                             cc_pvtx_coord,     &
                                             cc_pedge_vtx,      &
                                             cc_pface_edge_idx, &
                                             cc_pface_edge,     &
                                             cc_pface_vtx,      &
                                             cc_pvtx_ln_to_gn,  &
                                             cc_pedge_ln_to_gn, &
                                             cc_pface_ln_to_gn)

    call c_f_pointer(c_pn_vtx, &
                     pn_vtx,   &
                     [n_part])

    call c_f_pointer(c_pn_edge, &
                     pn_edge,   &
                     [n_part])

    call c_f_pointer(c_pn_face, &
                     pn_face,   &
                     [n_part])

    ! Vertices
    length(1:n_part) = pn_vtx(1:n_part)
    call PDM_pointer_array_create(pvtx_ln_to_gn,    &
                                  n_part,           &
                                  PDM_TYPE_G_NUM,   &
                                  cc_pvtx_ln_to_gn, &
                                  length,           &
                                  PDM_OWNERSHIP_KEEP)

    length(1:n_part) = pn_vtx(1:n_part) * 3
    call PDM_pointer_array_create(pvtx_coord,      &
                                  n_part,          &
                                  PDM_TYPE_DOUBLE, &
                                  cc_pvtx_coord,   &
                                  length,          &
                                  PDM_OWNERSHIP_KEEP)

    ! Edges
    length(1:n_part) = pn_edge(1:n_part)
    call PDM_pointer_array_create(pedge_ln_to_gn,    &
                                  n_part,            &
                                  PDM_TYPE_G_NUM,    &
                                  cc_pedge_ln_to_gn, &
                                  length,            &
                                  PDM_OWNERSHIP_KEEP)

    length(1:n_part) = pn_edge(1:n_part) * 2
    call PDM_pointer_array_create(pedge_vtx,    &
                                  n_part,       &
                                  PDM_TYPE_INT, &
                                  cc_pedge_vtx, &
                                  length,       &
                                  PDM_OWNERSHIP_KEEP)

    ! Faces
    length(1:n_part) = pn_face(1:n_part)
    call PDM_pointer_array_create(pface_ln_to_gn,    &
                                  n_part,            &
                                  PDM_TYPE_G_NUM,    &
                                  cc_pface_ln_to_gn, &
                                  length,            &
                                  PDM_OWNERSHIP_KEEP)

    length(1:n_part) = pn_face(1:n_part) + 1
    call PDM_pointer_array_create(pface_edge_idx,    &
                                  n_part,            &
                                  PDM_TYPE_INT,      &
                                  cc_pface_edge_idx, &
                                  length,            &
                                  PDM_OWNERSHIP_KEEP)

    do i_part = 1, n_part
      call c_f_pointer(cc_pface_edge_idx, fptr,         [n_part])
      call c_f_pointer(fptr(i_part),      face_edge_idx, [pn_face(i_part) + 1])
      length(i_part) = face_edge_idx(pn_face(i_part))
    enddo
    call PDM_pointer_array_create(pface_vtx,    &
                                  n_part,       &
                                  PDM_TYPE_INT, &
                                  cc_pface_vtx, &
                                  length,       &
                                  PDM_OWNERSHIP_KEEP)

    call PDM_pointer_array_create(pface_edge,    &
                                  n_part,        &
                                  PDM_TYPE_INT,  &
                                  cc_pface_edge, &
                                  length,        &
                                  PDM_OWNERSHIP_KEEP)

  end subroutine PDM_generate_mesh_rectangle_ngon

end module pdm_generate_mesh
