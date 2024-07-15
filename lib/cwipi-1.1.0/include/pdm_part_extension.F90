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

module pdm_part_extension

  use pdm
  use pdm_pointer_array

  implicit none

  integer, parameter :: PDM_EXTEND_FROM_FACE = 0
  integer, parameter :: PDM_EXTEND_FROM_EDGE = 1
  integer, parameter :: PDM_EXTEND_FROM_VTX  = 2

interface

! Compute a part extension structure

subroutine PDM_part_extension_compute (part_ext) &
bind (c, name='PDM_part_extension_compute')
  use iso_c_binding
  implicit none

  type(c_ptr), value :: part_ext ! Part Extension instance

end subroutine PDM_part_extension_compute

! Free a part extension structure

subroutine PDM_part_extension_free (part_ext) &
bind (c, name='PDM_part_extension_free')
  use iso_c_binding
  implicit none

  type(c_ptr), value :: part_ext ! Part Extension instance

end subroutine PDM_part_extension_free



subroutine PDM_part_extension_connectivity_set_c(part_ext, &
                                                 i_domain, &
                                                 i_part, &
                                                 connectivity_type, &
                                                 connect_idx, &
                                                 connect) &
bind(c, name="PDM_part_extension_connectivity_set")
  use iso_c_binding
  implicit none
  type(c_ptr),    value :: part_ext
  integer(c_int), value :: i_domain
  integer(c_int), value :: i_part
  integer(c_int), value :: connectivity_type
  type(c_ptr),    value :: connect_idx
  type(c_ptr),    value :: connect
end subroutine PDM_part_extension_connectivity_set_c

subroutine PDM_part_extension_ln_to_gn_set_c(part_ext, &
                                             i_domain, &
                                             i_part, &
                                             mesh_entity, &
                                             n_entity, &
                                             ln_to_gn) &
bind(c, name="PDM_part_extension_ln_to_gn_set")
  use iso_c_binding
  implicit none
  type(c_ptr),    value :: part_ext
  integer(c_int), value :: i_domain
  integer(c_int), value :: i_part
  integer(c_int), value :: mesh_entity
  integer(c_int), value :: n_entity
  type(c_ptr),    value :: ln_to_gn
end subroutine PDM_part_extension_ln_to_gn_set_c

subroutine PDM_part_extension_vtx_coord_set_c(part_ext, &
                                              i_domain, &
                                              i_part, &
                                              vtx_coord) &
  bind(c, name="PDM_part_extension_vtx_coord_set")
  use iso_c_binding
  implicit none
  type(c_ptr),    value :: part_ext
  integer(c_int), value :: i_domain
  integer(c_int), value :: i_part
  type(c_ptr),    value :: vtx_coord

end subroutine PDM_part_extension_vtx_coord_set_c

subroutine PDM_part_extension_part_bound_graph_set_c(part_ext, &
                                                     i_domain, &
                                                     i_part, &
                                                     entity_type, &
                                                     part_bound_proc_idx, &
                                                     part_bound_part_idx, &
                                                     part_bound) &
bind (c, name="PDM_part_extension_part_bound_graph_set")
  use iso_c_binding
  implicit none
  type(c_ptr),    value :: part_ext
  integer(c_int), value :: i_domain
  integer(c_int), value :: i_part
  integer(c_int), value :: entity_type
  type(c_ptr),    value :: part_bound_proc_idx
  type(c_ptr),    value :: part_bound_part_idx
  type(c_ptr),    value :: part_bound
end subroutine PDM_part_extension_part_bound_graph_set_c

subroutine PDM_part_extension_group_set_c(part_ext, &
                                          i_domain, &
                                          i_part, &
                                          entity_type, &
                                          n_group, &
                                          group_entity_idx, &
                                          group_entity, &
                                          group_entity_ln_to_gn) &
bind (c, name="PDM_part_extension_group_set")
  use iso_c_binding
  implicit none
  type(c_ptr),    value :: part_ext
  integer(c_int), value :: i_domain
  integer(c_int), value :: i_part
  integer(c_int), value :: entity_type
  integer(c_int), value :: n_group
  type(c_ptr),    value :: group_entity_idx
  type(c_ptr),    value :: group_entity
  type(c_ptr),    value :: group_entity_ln_to_gn
end subroutine PDM_part_extension_group_set_c

end interface


contains


subroutine PDM_part_extension_create (part_ext,    &
                                      n_domain,    &
                                      n_part,      &
                                      extend_type, &
                                      depth,       &
                                      comm,        &
                                      owner)
  ! Initialize a part extension structure
  !
  ! Admissible values for ``extend_type``:
  !   - ``PDM_EXTEND_FROM_FACE``
  !   - ``PDM_EXTEND_FROM_EDGE``
  !   - ``PDM_EXTEND_FROM_VTX``
  use iso_c_binding
  implicit none

  type(c_ptr)                   :: part_ext    ! Part Extension instance
  integer, intent(in)           :: n_domain    ! Number of domains
  integer(pdm_l_num_s), pointer :: n_part(:)   ! Number of partitions per domain (size = ``n_domain``)
  integer, intent(in)           :: extend_type ! Type of extension
  integer, intent(in)           :: depth       ! Extension depth
  integer, intent(in)           :: comm        ! MPI communicator
  integer, intent(in)           :: owner       ! Data ownership

  integer(c_int)                :: c_comm

  interface
    function PDM_part_extension_create_c (n_domain,    &
                                          n_part,      &
                                          extend_type, &
                                          depth,       &
                                          comm,        &
                                          owner)       &
    result (part_ext)                                  &
    bind (c, name='PDM_part_extension_create')
      use iso_c_binding
      implicit none

      integer(c_int), value :: n_domain
      type(c_ptr),    value :: n_part
      integer(c_int), value :: extend_type
      integer(c_int), value :: depth
      integer(c_int), value :: comm
      integer(c_int), value :: owner
      type(c_ptr)           :: part_ext

    end function PDM_part_extension_create_c
  end interface

  c_comm = PDM_MPI_Comm_f2c(comm)

  part_ext = PDM_part_extension_create_c (n_domain,      &
                                          c_loc(n_part), &
                                          extend_type,   &
                                          depth,         &
                                          c_comm,        &
                                          owner)

end subroutine PDM_part_extension_create

! Set data to perform the partitioned mesh extension

subroutine PDM_part_extension_set_part (part_ext,                 &
                                        i_domain,                 &
                                        i_part,                   &
                                        n_cell,                   &
                                        n_face,                   &
                                        n_face_part_bound,        &
                                        n_face_group,             &
                                        n_edge,                   &
                                        n_vtx,                    &
                                        cell_face_idx,            &
                                        cell_face,                &
                                        face_cell,                &
                                        face_edge_idx,            &
                                        face_edge,                &
                                        face_vtx_idx,             &
                                        face_vtx,                 &
                                        edge_vtx,                 &
                                        face_bound_idx,           &
                                        face_bound,               &
                                        face_join_idx,            &
                                        face_join,                &
                                        face_part_bound_proc_idx, &
                                        face_part_bound_part_idx, &
                                        face_part_bound,          &
                                        vtx_part_bound_proc_idx,  &
                                        vtx_part_bound_part_idx,  &
                                        vtx_part_bound,           &
                                        cell_ln_to_gn,            &
                                        face_ln_to_gn,            &
                                        edge_ln_to_gn,            &
                                        vtx_ln_to_gn,             &
                                        face_group_ln_to_gn,      &
                                        vtx_coord)
  !
  ! .. warning::
  !   Deprecated: use the individual setters instead
  use iso_c_binding
  implicit none

  type(c_ptr), value            :: part_ext                    ! Part Extension instance
  integer, intent(in)           :: i_domain                    ! Id of the domain
  integer, intent(in)           :: i_part                      ! Id of the partition
  integer, intent(in)           :: n_cell                      ! Number of cells
  integer, intent(in)           :: n_face                      ! Number of faces
  integer, intent(in)           :: n_face_part_bound           ! Number of partition boundary faces
  integer, intent(in)           :: n_face_group                ! Number of face groups
  integer, intent(in)           :: n_edge                      ! Number of edges
  integer, intent(in)           :: n_vtx                       ! Number of vertices
  integer(pdm_l_num_s), pointer :: cell_face_idx(:)            ! Cell-face connectivity index (size = \ref n_cell + 1)
  integer(pdm_l_num_s), pointer :: cell_face(:)                ! Cell-face connectivity (size = \ref cell_face_idx(\ref n_cell + 1))
  integer(pdm_l_num_s), pointer :: face_cell(:)                ! Face-cell connectivity (size = 2 * \ref n_face)
  integer(pdm_l_num_s), pointer :: face_edge_idx(:)            ! Face-edge connectivity index (size = \ref n_face + 1)
  integer(pdm_l_num_s), pointer :: face_edge(:)                ! Face-edge connectivity (size = \ref face_edge_idx(\ref n_face + 1))
  integer(pdm_l_num_s), pointer :: face_vtx_idx(:)             ! Face-vertex connectivity index (size = \ref n_face + 1)
  integer(pdm_l_num_s), pointer :: face_vtx(:)                 ! Face-vertex connectivity (size = \ref face_vtx_idx(\ref n_face + 1))
  integer(pdm_l_num_s), pointer :: edge_vtx(:)                 ! Edge-vertex connectivity (size = 2 * \ref n_edge)
  integer(pdm_l_num_s), pointer :: face_bound_idx(:)           ! Face->group connectivity index (size = \ref n_face_group + 1)
  integer(pdm_l_num_s), pointer :: face_bound(:)               ! Face->group connectivity (size = \ref face_edge_idx(\ref n_face_group + 1))
  integer(pdm_l_num_s), pointer :: face_join_idx(:)            ! Faces connecting domains connectivity index
  integer(pdm_l_num_s), pointer :: face_join(:)                ! Faces connecting domains connectivity
  integer(pdm_l_num_s), pointer :: face_part_bound_proc_idx(:) ! Partitioning boundary faces index from process (size = n_proc + 1)
  integer(pdm_l_num_s), pointer :: face_part_bound_part_idx(:) ! Partitioning boundary faces index from partition (size = n_total_part + 1)
  integer(pdm_l_num_s), pointer :: face_part_bound(:)          ! Partitioning boundary faces (size = 4 * n_face_part_bound)
  integer(pdm_l_num_s), pointer :: vtx_part_bound_proc_idx(:)  ! Partitioning boundary vertices index from process (size = n_proc + 1)
  integer(pdm_l_num_s), pointer :: vtx_part_bound_part_idx(:)  ! Partitioning boundary vertices index from partition (size = n_total_part + 1)
  integer(pdm_l_num_s), pointer :: vtx_part_bound(:)           ! Partitioning boundary vertices (size = 4 * n_vertex_part_bound)
  integer(pdm_g_num_s), pointer :: cell_ln_to_gn(:)            ! Cell global ids (size = \ref n_cell)
  integer(pdm_g_num_s), pointer :: face_ln_to_gn(:)            ! Face global ids (size = \ref n_face)
  integer(pdm_g_num_s), pointer :: edge_ln_to_gn(:)            ! Edge global ids (size = \ref n_edge)
  integer(pdm_g_num_s), pointer :: vtx_ln_to_gn(:)             ! Vertex global ids (size = \ref n_vtx)
  integer(pdm_g_num_s), pointer :: face_group_ln_to_gn(:)      ! Global ids of faces with groups (size = \ref n_face_group)
  double precision,     pointer :: vtx_coord(:,:)              ! Vertex coordinates (shape = [3, n_vtx])

  type(c_ptr) :: c_cell_face_idx
  type(c_ptr) :: c_cell_face
  type(c_ptr) :: c_face_cell
  type(c_ptr) :: c_face_edge_idx
  type(c_ptr) :: c_face_edge
  type(c_ptr) :: c_face_vtx_idx
  type(c_ptr) :: c_face_vtx
  type(c_ptr) :: c_edge_vtx
  type(c_ptr) :: c_face_bound_idx
  type(c_ptr) :: c_face_bound
  type(c_ptr) :: c_face_join_idx
  type(c_ptr) :: c_face_join
  type(c_ptr) :: c_face_part_bound_proc_idx
  type(c_ptr) :: c_face_part_bound_part_idx
  type(c_ptr) :: c_face_part_bound
  type(c_ptr) :: c_vtx_part_bound_proc_idx
  type(c_ptr) :: c_vtx_part_bound_part_idx
  type(c_ptr) :: c_vtx_part_bound
  type(c_ptr) :: c_cell_ln_to_gn
  type(c_ptr) :: c_face_ln_to_gn
  type(c_ptr) :: c_edge_ln_to_gn
  type(c_ptr) :: c_vtx_ln_to_gn
  type(c_ptr) :: c_face_group_ln_to_gn
  type(c_ptr) :: c_vtx_coord

  interface
    subroutine PDM_part_extension_set_part_c (part_ext,                 &
                                              i_domain,                 &
                                              i_part,                   &
                                              n_cell,                   &
                                              n_face,                   &
                                              n_face_part_bound,        &
                                              n_face_group,             &
                                              n_edge,                   &
                                              n_vtx,                    &
                                              cell_face_idx,            &
                                              cell_face,                &
                                              face_cell,                &
                                              face_edge_idx,            &
                                              face_edge,                &
                                              face_vtx_idx,             &
                                              face_vtx,                 &
                                              edge_vtx,                 &
                                              face_bound_idx,           &
                                              face_bound,               &
                                              face_join_idx,            &
                                              face_join,                &
                                              face_part_bound_proc_idx, &
                                              face_part_bound_part_idx, &
                                              face_part_bound,          &
                                              vtx_part_bound_proc_idx,  &
                                              vtx_part_bound_part_idx,  &
                                              vtx_part_bound,           &
                                              cell_ln_to_gn,            &
                                              face_ln_to_gn,            &
                                              edge_ln_to_gn,            &
                                              vtx_ln_to_gn,             &
                                              face_group_ln_to_gn,      &
                                              vtx_coord)                &
    bind (c, name='PDM_part_extension_set_part')
      use iso_c_binding
      implicit none

      type(c_ptr),    value :: part_ext
      integer(c_int), value :: i_domain
      integer(c_int), value :: i_part
      integer(c_int), value :: n_cell
      integer(c_int), value :: n_face
      integer(c_int), value :: n_face_part_bound
      integer(c_int), value :: n_face_group
      integer(c_int), value :: n_edge
      integer(c_int), value :: n_vtx
      type(c_ptr),    value :: cell_face_idx
      type(c_ptr),    value :: cell_face
      type(c_ptr),    value :: face_cell
      type(c_ptr),    value :: face_edge_idx
      type(c_ptr),    value :: face_edge
      type(c_ptr),    value :: face_vtx_idx
      type(c_ptr),    value :: face_vtx
      type(c_ptr),    value :: edge_vtx
      type(c_ptr),    value :: face_bound_idx
      type(c_ptr),    value :: face_bound
      type(c_ptr),    value :: face_join_idx
      type(c_ptr),    value :: face_join
      type(c_ptr),    value :: face_part_bound_proc_idx
      type(c_ptr),    value :: face_part_bound_part_idx
      type(c_ptr),    value :: face_part_bound
      type(c_ptr),    value :: vtx_part_bound_proc_idx
      type(c_ptr),    value :: vtx_part_bound_part_idx
      type(c_ptr),    value :: vtx_part_bound
      type(c_ptr),    value :: cell_ln_to_gn
      type(c_ptr),    value :: face_ln_to_gn
      type(c_ptr),    value :: edge_ln_to_gn
      type(c_ptr),    value :: vtx_ln_to_gn
      type(c_ptr),    value :: face_group_ln_to_gn
      type(c_ptr),    value :: vtx_coord

    end subroutine PDM_part_extension_set_part_c
  end interface

  c_cell_face_idx = C_NULL_PTR
  if (associated(cell_face_idx)) then
    c_cell_face_idx = c_loc(cell_face_idx)
  endif 
    
  c_cell_face = C_NULL_PTR
  if (associated(cell_face)) then
    c_cell_face = c_loc(cell_face)
  endif 
    
  c_face_cell = C_NULL_PTR
  if (associated(face_cell)) then
    c_face_cell = c_loc(face_cell)
  endif 
    
  c_face_edge_idx = C_NULL_PTR
  if (associated(face_edge_idx)) then
    c_face_edge_idx = c_loc(face_edge_idx)
  endif 
    
  c_face_edge = C_NULL_PTR
  if (associated(face_edge)) then
    c_face_edge = c_loc(face_edge)
  endif 
    
  c_face_vtx_idx = C_NULL_PTR
  if (associated(face_vtx_idx)) then
    c_face_vtx_idx = c_loc(face_vtx_idx)
  endif 
    
  c_face_vtx = C_NULL_PTR
  if (associated(face_vtx)) then
    c_face_vtx = c_loc(face_vtx)
  endif 
    
  c_edge_vtx = C_NULL_PTR
  if (associated(edge_vtx)) then
    c_edge_vtx = c_loc(edge_vtx)
  endif 
    
  c_face_bound_idx = C_NULL_PTR
  if (associated(face_bound_idx)) then
    c_face_bound_idx = c_loc(face_bound_idx)
  endif 
    
  c_face_bound = C_NULL_PTR
  if (associated(face_bound)) then
    c_face_bound = c_loc(face_bound)
  endif 
    
  c_face_join_idx = C_NULL_PTR
  if (associated(face_join_idx)) then
    c_face_join_idx = c_loc(face_join_idx)
  endif 
    
  c_face_join = C_NULL_PTR
  if (associated(face_join)) then
    c_face_join = c_loc(face_join)
  endif 
    
  c_face_part_bound_proc_idx = C_NULL_PTR
  if (associated(face_part_bound_proc_idx)) then
    c_face_part_bound_proc_idx = c_loc(face_part_bound_proc_idx)
  endif 
    
  c_face_part_bound_part_idx = C_NULL_PTR
  if (associated(face_part_bound_part_idx)) then
    c_face_part_bound_part_idx = c_loc(face_part_bound_part_idx)
  endif 
    
  c_face_part_bound = C_NULL_PTR
  if (associated(face_part_bound)) then
    c_face_part_bound = c_loc(face_part_bound)
  endif 
    
  c_vtx_part_bound_proc_idx = C_NULL_PTR
  if (associated(vtx_part_bound_proc_idx)) then
    c_vtx_part_bound_proc_idx = c_loc(vtx_part_bound_proc_idx)
  endif 
    
  c_vtx_part_bound_part_idx = C_NULL_PTR
  if (associated(vtx_part_bound_part_idx)) then
    c_vtx_part_bound_part_idx = c_loc(vtx_part_bound_part_idx)
  endif 
    
  c_vtx_part_bound = C_NULL_PTR
  if (associated(vtx_part_bound)) then
    c_vtx_part_bound = c_loc(vtx_part_bound)
  endif 
    
  c_cell_ln_to_gn = C_NULL_PTR
  if (associated(cell_ln_to_gn)) then
    c_cell_ln_to_gn = c_loc(cell_ln_to_gn)
  endif 
    
  c_face_ln_to_gn = C_NULL_PTR
  if (associated(face_ln_to_gn)) then
    c_face_ln_to_gn = c_loc(face_ln_to_gn)
  endif 
    
  c_edge_ln_to_gn = C_NULL_PTR
  if (associated(edge_ln_to_gn)) then
    c_edge_ln_to_gn = c_loc(edge_ln_to_gn)
  endif 
    
  c_vtx_ln_to_gn = C_NULL_PTR
  if (associated(vtx_ln_to_gn)) then
    c_vtx_ln_to_gn = c_loc(vtx_ln_to_gn)
  endif 
    
  c_face_group_ln_to_gn = C_NULL_PTR
  if (associated(face_group_ln_to_gn)) then
    c_face_group_ln_to_gn = c_loc(face_group_ln_to_gn)
  endif 
    
  c_vtx_coord = C_NULL_PTR
  if (associated(vtx_coord)) then
    c_vtx_coord = c_loc(vtx_coord)
  endif 
    

  call PDM_part_extension_set_part_c (part_ext,                        &
                                      i_domain,                        &
                                      i_part,                          &
                                      n_cell,                          &
                                      n_face,                          &
                                      n_face_part_bound,               &
                                      n_face_group,                    &
                                      n_edge,                          &
                                      n_vtx,                           &
                                      c_cell_face_idx,                 &
                                      c_cell_face,                     &
                                      c_face_cell,                     &
                                      c_face_edge_idx,                 &
                                      c_face_edge,                     &
                                      c_face_vtx_idx,                  &
                                      c_face_vtx,                      &
                                      c_edge_vtx,                      &
                                      c_face_bound_idx,                &
                                      c_face_bound,                    &
                                      c_face_join_idx,                 &
                                      c_face_join,                     &
                                      c_face_part_bound_proc_idx,      &
                                      c_face_part_bound_part_idx,      &
                                      c_face_part_bound,               &
                                      c_vtx_part_bound_proc_idx,       &
                                      c_vtx_part_bound_part_idx,       &
                                      c_vtx_part_bound,                &
                                      c_cell_ln_to_gn,                 &
                                      c_face_ln_to_gn,                 &
                                      c_edge_ln_to_gn,                 &
                                      c_vtx_ln_to_gn,                  &
                                      c_face_group_ln_to_gn,           &
                                      c_vtx_coord)

end subroutine PDM_part_extension_set_part


subroutine PDM_part_extension_connectivity_get (part_ext,          &
                                                i_domain,          &
                                                i_part,            &
                                                connectivity_type, &
                                                n_elt,             &
                                                connect_idx,       &
                                                connect)
  ! Get extended connectivity
  use iso_c_binding
  implicit none

  type(c_ptr), value            :: part_ext          ! Part Extension instance
  integer, intent(in)           :: i_domain          ! Id of current domain
  integer, intent(in)           :: i_part            ! Id of current partition
  integer, intent(in)           :: connectivity_type ! Type of mesh entity
  integer, intent(out)          :: n_elt             ! Number of elements
  integer(pdm_l_num_s), pointer :: connect_idx(:)    ! Index for entity->group graph (size = ``n_elt`` + 1)
  integer(pdm_l_num_s), pointer :: connect(:)        ! Entity->group graph (size = ``connect_idx(n_elt)``)

  type(c_ptr)                   :: c_connect
  type(c_ptr)                   :: c_connect_idx

  interface
    function PDM_part_extension_connectivity_get_c (part_ext,          &
                                                    i_domain,          &
                                                    i_part,            &
                                                    connectivity_type, &
                                                    connect_idx,       &
                                                    connect)           &
    result (n_elt)                                                     &
    bind (c, name='PDM_part_extension_connectivity_get')
      use iso_c_binding
      implicit none

      type(c_ptr),    value :: part_ext
      integer(c_int), value :: i_domain
      integer(c_int), value :: i_part
      integer(c_int), value :: connectivity_type
      type(c_ptr)           :: connect_idx
      type(c_ptr)           :: connect
      integer(c_int)        :: n_elt

    end function PDM_part_extension_connectivity_get_c
  end interface

  c_connect     = C_NULL_PTR
  c_connect_idx = C_NULL_PTR
  n_elt =  PDM_part_extension_connectivity_get_c (part_ext,          &
                                                  i_domain,          &
                                                  i_part,            &
                                                  connectivity_type, &
                                                  c_connect_idx,     &
                                                  c_connect)

  call c_f_pointer(c_connect_idx, &
                   connect_idx,   &
                   [n_elt+1])

  call c_f_pointer(c_connect, &
                   connect,   &
                   [connect_idx(n_elt + 1)])

end subroutine PDM_part_extension_connectivity_get


subroutine PDM_part_extension_ln_to_gn_get (part_ext,    &
                                            i_domain,    &
                                            i_part,      &
                                            mesh_entity, &
                                            n_elt,       &
                                            ln_to_gn)
  ! Get global ids of extended entities
  use iso_c_binding
  implicit none

  type(c_ptr), value            :: part_ext    ! Part Extension instance
  integer, intent(in)           :: i_domain    ! Id of current domain
  integer, intent(in)           :: i_part      ! Id of current partition
  integer, intent(in)           :: mesh_entity ! Type of mesh entity
  integer, intent(out)          :: n_elt       ! Number of elements
  integer(pdm_g_num_s), pointer :: ln_to_gn(:) ! Global ids (size = ``n_elt``)

  type(c_ptr)                   :: c_ln_to_gn

  interface
    function PDM_part_extension_ln_to_gn_get_c (part_ext,    &
                                                i_domain,    &
                                                i_part,      &
                                                mesh_entity, &
                                                ln_to_gn)    &
    result (n_elt)                                           &
    bind (c, name='PDM_part_extension_ln_to_gn_get')
      use iso_c_binding
      implicit none

      type(c_ptr),    value :: part_ext
      integer(c_int), value :: i_domain
      integer(c_int), value :: i_part
      integer(c_int), value :: mesh_entity
      type(c_ptr)           :: ln_to_gn
      integer(c_int)        :: n_elt

    end function PDM_part_extension_ln_to_gn_get_c
  end interface

  c_ln_to_gn = C_NULL_PTR
  n_elt =  PDM_part_extension_ln_to_gn_get_c (part_ext,    &
                                              i_domain,    &
                                              i_part,      &
                                              mesh_entity, &
                                              c_ln_to_gn)

  call c_f_pointer(c_ln_to_gn, &
                   ln_to_gn,   &
                   [n_elt])

end subroutine PDM_part_extension_ln_to_gn_get


subroutine PDM_part_extension_group_get (part_ext,         &
                                         i_domain,         &
                                         i_part,           &
                                         mesh_entity,      &
                                         n_group,          &
                                         group_entity_idx, &
                                         group_entity,     &
                                         group_entity_ln_to_gn)
  ! Get groups for extended entities with given type
  use iso_c_binding
  implicit none

  type(c_ptr), value            :: part_ext                 ! Part Extension instance
  integer, intent(in)           :: i_domain                 ! Id of current domain
  integer, intent(in)           :: i_part                   ! Id of current partition
  integer, intent(in)           :: mesh_entity              ! Type of mesh entity
  integer, intent(out)          :: n_group                  ! Number of groups
  integer(pdm_l_num_s), pointer :: group_entity_idx(:)      ! Index for group->entity connectivity (size = ``n_group``)
  integer(pdm_l_num_s), pointer :: group_entity(:)          ! Group->entity connectivity (1-based local ids, size = ``group_entity_idx(n_group)``)
  integer(pdm_g_num_s), pointer :: group_entity_ln_to_gn(:) ! Group->entity connectivity (group-specific global ids, size = ``group_entity_idx(n_group)``)

  type(c_ptr)                   :: c_group_entity          
  type(c_ptr)                   :: c_group_entity_idx      
  type(c_ptr)                   :: c_group_entity_ln_to_gn 

  interface
    function PDM_part_extension_group_get_c (part_ext,              &
                                             i_domain,              &
                                             i_part,                &
                                             entity_type,           &
                                             group_entity_idx,      &
                                             group_entity,          &
                                             group_entity_ln_to_gn) &
    result (n_group)                                                &
    bind (c, name='PDM_part_extension_group_get')
      use iso_c_binding
      implicit none

      type(c_ptr),    value :: part_ext
      integer(c_int), value :: i_domain
      integer(c_int), value :: i_part
      integer(c_int), value :: entity_type
      type(c_ptr)           :: group_entity_idx
      type(c_ptr)           :: group_entity
      type(c_ptr)           :: group_entity_ln_to_gn
      integer(c_int)        :: n_group

    end function PDM_part_extension_group_get_c
  end interface

  c_group_entity          = C_NULL_PTR
  c_group_entity_idx      = C_NULL_PTR
  c_group_entity_ln_to_gn = C_NULL_PTR
  n_group =  PDM_part_extension_group_get_c (part_ext,           &
                                             i_domain,           &
                                             i_part,             &
                                             mesh_entity,        &
                                             c_group_entity_idx, &
                                             c_group_entity,     &
                                             c_group_entity_ln_to_gn)

  call c_f_pointer(c_group_entity_idx, &
                   group_entity_idx,   &
                   [n_group+1])

  call c_f_pointer(c_group_entity, &
                   group_entity,   &
                   [group_entity_idx(n_group + 1)])

  call c_f_pointer(c_group_entity_ln_to_gn, &
                   group_entity_ln_to_gn,   &
                   [group_entity_idx(n_group + 1)])

end subroutine PDM_part_extension_group_get


subroutine PDM_part_extension_vtx_coord_get (part_ext,    &
                                             i_domain,    &
                                             i_part,      &
                                             n_vtx,       &
                                             vtx_coord)
  ! Get coordinates of extended vertices
  use iso_c_binding
  implicit none

  type(c_ptr), value        :: part_ext        ! Part Extension instance
  integer, intent(in)       :: i_domain        ! Id of current domain
  integer, intent(in)       :: i_part          ! Id of current partition
  integer, intent(out)      :: n_vtx           ! Number of vertices
  double precision, pointer :: vtx_coord(:,:)  ! Vertex coordinates (shape = [3, ``n_vtx``])

  type(c_ptr)               :: c_vtx_coord

  interface
    function PDM_part_extension_vtx_coord_get_c (part_ext,    &
                                                 i_domain,    &
                                                 i_part,      &
                                                 vtx_coord)   &
    result (n_vtx)                                            &
    bind (c, name='PDM_part_extension_vtx_coord_get')
      use iso_c_binding
      implicit none

      type(c_ptr),    value :: part_ext
      integer(c_int), value :: i_domain
      integer(c_int), value :: i_part
      type(c_ptr)           :: vtx_coord
      integer(c_int)        :: n_vtx

    end function PDM_part_extension_vtx_coord_get_c
  end interface

  c_vtx_coord = C_NULL_PTR
  n_vtx =  PDM_part_extension_vtx_coord_get_c (part_ext,    &
                                               i_domain,    &
                                               i_part,      &
                                               c_vtx_coord)

  call c_f_pointer(c_vtx_coord, &
                   vtx_coord,   &
                   [3,n_vtx])

end subroutine PDM_part_extension_vtx_coord_get


!>
!! \brief Set a domain interface
!!
!! \param [in]  part_ext        Pointer to \ref PDM_part_extension_t object
!! \param [in]  i_interface     Id of the interface
!! \param [in]  interface_kind  Kind of the interface
!! \param [in]  n_pair          Number of pairs
!! \param [in]  interface_ids   Global ids (size = 2 * \ref n_pair)
!! \param [in]  interface_dom   Ids of the domains (size = 2 * \ref n_pair or 2)
!!

! subroutine PDM_part_extension_domain_interface_set (part_ext,       &
!                                                     i_interface,    &
!                                                     interface_kind, &
!                                                     n_pair,         &
!                                                     interface_ids,  &
!                                                     interface_dom)
!   use iso_c_binding
!   implicit none

!   type(c_ptr), value            :: part_ext
!   integer                       :: i_interface
!   integer                       :: interface_kind
!   integer                       :: n_pair
!   integer(pdm_g_num_s), pointer :: interface_ids(:)
!   integer(pdm_l_num_s), pointer :: interface_dom(:)

!   type(c_ptr)                   :: c_interface_ids = C_NULL_PTR
!   type(c_ptr)                   :: c_interface_dom = C_NULL_PTR

!   interface
!     subroutine PDM_part_extension_domain_interface_set_c (part_ext,       &
!                                                           i_interface,    &
!                                                           interface_kind, &
!                                                           n_pair,         &
!                                                           interface_ids,  &
!                                                           interface_dom)  &
!     bind (c, name='PDM_part_extension_domain_interface_set')
!       use iso_c_binding
!       implicit none

!       type(c_ptr),    value :: part_ext
!       integer(c_int), value :: i_interface
!       integer(c_int), value :: interface_kind
!       integer(c_int), value :: n_pair
!       type(c_ptr),    value :: interface_ids
!       type(c_ptr),    value :: interface_dom

!     end subroutine PDM_part_extension_domain_interface_set_c
!   end interface

!   c_interface_ids = c_loc(interface_ids)
!   c_interface_dom = c_loc(interface_dom)

!   call PDM_part_extension_domain_interface_set_c (part_ext,         &
!                                                   i_interface,      &
!                                                   interface_kind,   &
!                                                   n_pair,           &
!                                                   c_interface_ids,  &
!                                                   c_interface_dom)

! end subroutine PDM_part_extension_domain_interface_set

!>
!! \brief Set the translation vector for a periodic interface
!!
!! \param [in]  part_ext     Pointer to \ref PDM_part_extension_t object
!! \param [in]  i_interface  Id of the interface
!! \param [in]  vect         Translation vector (size = 3)
!!

! subroutine PDM_part_extension_domain_interface_translation_set (part_ext,    &
!                                                                 i_interface, &
!                                                                 vect)
!   use iso_c_binding
!   implicit none

!   type(c_ptr), value        :: part_ext
!   integer                   :: i_interface
!   double precision, pointer :: vect(:)

!   type(c_ptr)               :: c_vect = C_NULL_PTR

!   interface
!     subroutine PDM_part_extension_domain_interface_translation_set_c (part_ext,    &
!                                                                       i_interface, &
!                                                                       vect)        &
!     bind (c, name='PDM_part_extension_domain_interface_translation_set')
!       use iso_c_binding
!       implicit none

!       type(c_ptr),    value :: part_ext
!       integer(c_int), value :: i_interface
!       type(c_ptr),    value :: vect

!     subroutine PDM_part_extension_domain_interface_translation_set_c
!   end interface

!   c_vect = c_loc(vect)

!   call PDM_part_extension_domain_interface_translation_set_c (part_ext,    &
!                                                               i_interface, &
!                                                               c_vect)

! end subroutine PDM_part_extension_domain_interface_translation_set


!>
!! \brief Set the rotation for a periodic interface
!!
!! \param [in]  part_ext     Pointer to \ref PDM_part_extension_t object
!! \param [in]  i_interface  Id of the interface
!! \param [in]  direction    Rotation axis (size = 3)
!! \param [in]  center       Rotation center (size = 3)
!! \param [in]  angle        Rotation angle
!!

! subroutine PDM_part_extension_domain_interface_rotation_set (part_ext,    &
!                                                              i_interface, &
!                                                              direction,   &
!                                                              center,      &
!                                                              angle)
!   use iso_c_binding
!   implicit none

!   type(c_ptr), value        :: part_ext
!   integer                   :: i_interface
!   double precision, pointer :: direction(:)
!   double precision, pointer :: center(:)
!   double precision          :: angle

!   type(c_ptr)               :: c_direction = C_NULL_PTR
!   type(c_ptr)               :: c_center    = C_NULL_PTR

!   interface
!     subroutine PDM_part_extension_domain_interface_rotation_set_c (part_ext,    &
!                                                                    i_interface, &
!                                                                    direction,   &
!                                                                    center,      &
!                                                                    angle)       &
!     bind (c, name='PDM_part_extension_domain_interface_rotation_set')
!       use iso_c_binding
!       implicit none

!       type(c_ptr),    value :: part_ext
!       integer(c_int), value :: i_interface
!       type(c_ptr),    value :: direction
!       type(c_ptr),    value :: center
!       real(c_double), value :: vect

!     subroutine PDM_part_extension_domain_interface_rotation_set_c
!   end interface

!   c_direction = c_loc(direction)
!   c_center    = c_loc(center)

!   call PDM_part_extension_domain_interface_rotation_set_c (part_ext,    &
!                                                            i_interface, &
!                                                            c_direction, &
!                                                            c_center,    &
!                                                            angle)

! end subroutine PDM_part_extension_domain_interface_rotation_set


!>
!!
!! \brief Create part_to_part from interior and extended elements
!!
!! \param [in]   ptp                             Part to part structure
!! \param [in]   n_part                          Number of partitions
!! \param [in]   n_int_cell                      Number of interior elements
!! \param [in]   int_cell_ln_to_gn               gnum of interior elements
!! \param [in]   n_ext_cell                      Number of extended elements
!! \param [in]   ext_cell_ln_to_gn               gnum of extended elements
!! \param [out]  n_selected_cell_to_send         Number of elements selected for send
!! \param [out]  n_selected_cell_to_send_idx     Index of elements selected for send
!! \param [out]  selected_cell_to_send           Local numbering of elements selected for send
!! \param [out]  selected_cell_to_send_ln_to_gn  gnum of elements selected for send
!!

subroutine PDM_part_to_part_create_from_extension (ptp,                            &
                                                   n_part,                         &
                                                   n_int_cell,                     &
                                                   int_cell_ln_to_gn,              &
                                                   n_ext_cell,                     &
                                                   ext_cell_ln_to_gn,              &
                                                   n_selected_cell_to_send,        &
                                                   selected_cell_to_send,          &
                                                   comm)
  use iso_c_binding
  implicit none

  type(c_ptr)                        :: ptp
  integer, intent(in)                :: n_part
  type(PDM_pointer_array_t), pointer :: int_cell_ln_to_gn
  integer(pdm_l_num_s), pointer      :: n_int_cell(:)
  type(PDM_pointer_array_t), pointer :: ext_cell_ln_to_gn
  integer(pdm_l_num_s), pointer      :: n_ext_cell(:)
  type(PDM_pointer_array_t), pointer :: selected_cell_to_send
  integer(pdm_l_num_s), pointer      :: n_selected_cell_to_send(:)
  integer, intent(in)                :: comm

  integer                           :: i
  integer(c_int)                    :: c_comm
  type(c_ptr)                       :: c_selected_cell_to_send_idx
  type(c_ptr)                       :: c_selected_cell_to_send
  type(c_ptr)                       :: c_selected_cell_to_send_ln_to_gn
  type(c_ptr)                       :: c_n_selected_cell_to_send
  integer, allocatable              :: length_selected_cell_to_send(:)

  interface
    subroutine PDM_part_to_part_create_from_extension_c (ptp,                            &
                                                         n_part,                         &
                                                         n_int_cell,                     &
                                                         int_cell_ln_to_gn,              &
                                                         n_ext_cell,                     &
                                                         ext_cell_ln_to_gn,              &
                                                         n_selected_cell_to_send,        &
                                                         selected_cell_to_send,          &
                                                         comm)                           &
    bind(c, name='PDM_part_to_part_create_from_extension')
      use iso_c_binding
      implicit none

      type(c_ptr)           :: ptp
      integer(c_int), value :: n_part
      type(c_ptr),    value :: int_cell_ln_to_gn
      type(c_ptr),    value :: n_int_cell
      type(c_ptr),    value :: ext_cell_ln_to_gn
      type(c_ptr),    value :: n_ext_cell
      type(c_ptr)           :: selected_cell_to_send
      type(c_ptr)           :: n_selected_cell_to_send
      integer(c_int), value :: comm

    end subroutine PDM_part_to_part_create_from_extension_c
  end interface

  c_comm = PDM_MPI_Comm_f2c(comm)

  c_selected_cell_to_send_idx      = C_NULL_PTR
  c_selected_cell_to_send          = C_NULL_PTR
  c_selected_cell_to_send_ln_to_gn = C_NULL_PTR
  c_n_selected_cell_to_send        = C_NULL_PTR
  call PDM_part_to_part_create_from_extension_c (ptp,                              &
                                                 n_part,                           &
                                                 c_loc(n_int_cell),                &
                                                 c_loc(int_cell_ln_to_gn%cptr),    &
                                                 c_loc(n_ext_cell),                &
                                                 c_loc(ext_cell_ln_to_gn%cptr),    &
                                                 c_n_selected_cell_to_send,        &
                                                 c_selected_cell_to_send,          &
                                                 c_comm)

  call c_f_pointer(c_n_selected_cell_to_send,  &
                   n_selected_cell_to_send,    &
                   [n_part])

  allocate( length_selected_cell_to_send(n_part) )

  ! Compute lengths
  do i = 1, n_part
    length_selected_cell_to_send = n_selected_cell_to_send(i)
  end do

  call PDM_pointer_array_create (selected_cell_to_send,             &
                                 n_part,                            &
                                 PDM_TYPE_INT,                      &
                                 c_selected_cell_to_send,           &
                                 length_selected_cell_to_send,      &
                                 PDM_OWNERSHIP_KEEP)

  deallocate( length_selected_cell_to_send )

end subroutine PDM_part_to_part_create_from_extension



subroutine PDM_part_extension_connectivity_set(part_ext, &
                                               i_domain, &
                                               i_part, &
                                               connectivity_type, &
                                               connect_idx, &
                                               connect)
  ! Set connectivity
  use iso_c_binding
  implicit none
  type(c_ptr), intent(in)       :: part_ext            ! Part Extension instance
  integer,     intent(in)       :: i_domain            ! Domain identifier
  integer,     intent(in)       :: i_part              ! Partition identifier
  integer,     intent(in)       :: connectivity_type   ! Type of connectivity
  integer(pdm_l_num_s), pointer :: connect_idx(:)      ! Index for connectivity (can be \p NULL for
  integer(pdm_l_num_s), pointer :: connect(:)          ! Connectivity

  call PDM_part_extension_connectivity_set_c(part_ext, &
                                             i_domain, &
                                             i_part, &
                                             connectivity_type, &
                                             c_loc(connect_idx), &
                                             c_loc(connect))
end subroutine PDM_part_extension_connectivity_set


subroutine PDM_part_extension_ln_to_gn_set(part_ext, &
                                           i_domain, &
                                           i_part, &
                                           mesh_entity, &
                                           n_entity, &
                                           ln_to_gn)
  ! Set global ids
  use iso_c_binding
  implicit none
  type(c_ptr), intent(in)       :: part_ext    ! Part Extension instance
  integer,     intent(in)       :: i_domain    ! Domain identifier
  integer,     intent(in)       :: i_part      ! Partition identifier
  integer,     intent(in)       :: mesh_entity ! Type of mesh entity
  integer,     intent(in)       :: n_entity    ! Local number of entities
  integer(pdm_g_num_s), pointer :: ln_to_gn(:) ! Global ids

  call PDM_part_extension_ln_to_gn_set_c(part_ext, &
                                         i_domain, &
                                         i_part, &
                                         mesh_entity, &
                                         n_entity, &
                                         c_loc(ln_to_gn))
end subroutine PDM_part_extension_ln_to_gn_set


subroutine PDM_part_extension_vtx_coord_set(part_ext, &
                                            i_domain, &
                                            i_part, &
                                            vtx_coord)
  ! Set vertex coordinates
  use iso_c_binding
  implicit none
  type(c_ptr), intent(in)       :: part_ext       ! Part Extension instance
  integer,     intent(in)       :: i_domain       ! Domain identifier
  integer,     intent(in)       :: i_part         ! Partition identifier
  real(8),     pointer          :: vtx_coord(:,:) ! Vertex coordinates (shape = [3, *n_vtx*])

  call PDM_part_extension_vtx_coord_set_c(part_ext, &
                                          i_domain, &
                                          i_part, &
                                          c_loc(vtx_coord))
end subroutine PDM_part_extension_vtx_coord_set


subroutine PDM_part_extension_part_bound_graph_set(part_ext, &
                                                   i_domain, &
                                                   i_part, &
                                                   entity_type, &
                                                   part_bound_proc_idx, &
                                                   part_bound_part_idx, &
                                                   part_bound)
  ! Set the connection graph between partitions for the requested entity type
  use iso_c_binding
  implicit none
  type(c_ptr), intent(in)       :: part_ext               ! Part Extension instance
  integer,     intent(in)       :: i_domain               ! Domain identifier
  integer,     intent(in)       :: i_part                 ! Partition identifier
  integer,     intent(in)       :: entity_type            ! Type of mesh entity
  integer(pdm_l_num_s), pointer :: part_bound_proc_idx(:) ! Partitioning boundary entities index from process (size = *n_rank* + 1)
  integer(pdm_l_num_s), pointer :: part_bound_part_idx(:) ! Partitioning boundary entities index from partition (size = *n_total_part* + 1)
  integer(pdm_l_num_s), pointer :: part_bound(:)          ! Partitioning boundary entities (size = 4 * ``part_bound_proc_idx(n_rank)``)

  call PDM_part_extension_part_bound_graph_set_c(part_ext, &
                                                 i_domain, &
                                                 i_part, &
                                                 entity_type, &
                                                 c_loc(part_bound_proc_idx), &
                                                 c_loc(part_bound_part_idx), &
                                                 c_loc(part_bound))
end subroutine PDM_part_extension_part_bound_graph_set


subroutine PDM_part_extension_group_set(part_ext, &
                                        i_domain, &
                                        i_part, &
                                        entity_type, &
                                        n_group, &
                                        group_entity_idx, &
                                        group_entity, &
                                        group_entity_ln_to_gn)
  ! Set group description
  use iso_c_binding
  implicit none
  type(c_ptr), intent(in)       :: part_ext                 ! Part Extension instance
  integer,     intent(in)       :: i_domain                 ! Domain identifier
  integer,     intent(in)       :: i_part                   ! Partition identifier
  integer,     intent(in)       :: entity_type              ! Type of mesh entity
  integer,     intent(in)       :: n_group                  ! Number of groups
  integer(pdm_l_num_s), pointer :: group_entity_idx(:)      ! Index for group->entity connectivity (size = ``n_group``)
  integer(pdm_l_num_s), pointer :: group_entity(:)          ! Group->entity connectivity (1-based local ids, size = ``group_entity_idx(n_group)``)
  integer(pdm_g_num_s), pointer :: group_entity_ln_to_gn(:) ! Group->entity connectivity (group-specific global ids, size = ``group_entity_idx(``n_group)``)

  call PDM_part_extension_group_set_c(part_ext, &
                                      i_domain, &
                                      i_part, &
                                      entity_type, &
                                      n_group, &
                                      c_loc(group_entity_idx), &
                                      c_loc(group_entity), &
                                      c_loc(group_entity_ln_to_gn))
end subroutine PDM_part_extension_group_set

end module pdm_part_extension
