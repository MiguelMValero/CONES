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

module pdm_multipart

  use pdm

  implicit none

  integer, parameter :: PDM_PART_SIZE_HOMOGENEOUS   = 1
  integer, parameter :: PDM_PART_SIZE_HETEROGENEOUS = 2

  interface PDM_multipart_create ; module procedure  &
    PDM_multipart_create_
  end interface

  interface PDM_multipart_set_reordering_options ; module procedure  &
    PDM_multipart_set_reordering_options_
  end interface

  interface PDM_multipart_set_reordering_options_vtx ; module procedure  &
    PDM_multipart_set_reordering_options_vtx_
  end interface

  interface PDM_multipart_get_part_mesh_nodal ; module procedure  &
    PDM_multipart_get_part_mesh_nodal_
  end interface

  interface PDM_multipart_block_set ; module procedure  &
    PDM_multipart_block_set_
  end interface

  interface PDM_multipart_part_dim_get ; module procedure  &
    PDM_multipart_part_dim_get_
  end interface

  interface PDM_multipart_part_val_get ; module procedure  &
    PDM_multipart_part_val_get_
  end interface

  interface PDM_multipart_part_connectivity_get ; module procedure  &
    PDM_multipart_part_connectivity_get_
  end interface

  interface PDM_multipart_part_ln_to_gn_get ; module procedure  &
    PDM_multipart_part_ln_to_gn_get_
  end interface

  interface PDM_multipart_partition_color_get ; module procedure  &
    PDM_multipart_partition_color_get_
  end interface

  interface PDM_multipart_part_ghost_infomation_get ; module procedure  &
    PDM_multipart_part_ghost_infomation_get_
  end interface

  interface PDM_multipart_part_vtx_coord_get ; module procedure  &
    PDM_multipart_part_vtx_coord_get_
  end interface

  interface PDM_multipart_group_get ; module procedure  &
    PDM_multipart_group_get_
  end interface

  interface PDM_multipart_part_graph_comm_get ; module procedure  &
    PDM_multipart_part_graph_comm_get_
  end interface

interface

  function PDM_multipart_create_c (n_domain, &
                                   n_part, &
                                   merge_domains, &
                                   split_method, &
                                   part_size_method, &
                                   part_fraction, &
                                   comm, &
                                   owner) &

  result(multipart) &
  bind (c, name='PDM_multipart_create')

    use iso_c_binding
    implicit none

    type(c_ptr)            :: multipart
    integer(c_int), value  :: n_domain
    type(c_ptr),    value  :: n_part
    integer(c_int), value  :: merge_domains
    integer(c_int), value  :: split_method
    integer(c_int), value  :: part_size_method
    type(c_ptr),    value  :: part_fraction
    integer(c_int), value  :: comm
    integer(c_int), value  :: owner

  end function PDM_multipart_create_c

  subroutine PDM_multipart_set_reordering_options_c (multipart, &
                                                     i_domain, &
                                                     renum_cell_method, &
                                                     renum_cell_properties, &
                                                     renum_face_method) &
  bind (c, name='PDM_multipart_set_reordering_options')

    use iso_c_binding
    implicit none

    type(c_ptr),    value  :: multipart
    integer(c_int), value  :: i_domain
    character(c_char)      :: renum_cell_method(*)
    type(c_ptr),    value  :: renum_cell_properties
    character(c_char)      :: renum_face_method(*)

  end subroutine PDM_multipart_set_reordering_options_c


  subroutine PDM_multipart_set_reordering_options_vtx_c (multipart, &
                                                         i_domain, &
                                                         renum_vtx_method) &
  bind (c, name='PDM_multipart_set_reordering_options_vtx')

    use iso_c_binding
    implicit none

    type(c_ptr),    value  :: multipart
    integer(c_int), value  :: i_domain
    character(c_char)      :: renum_vtx_method(*)

  end subroutine PDM_multipart_set_reordering_options_vtx_c


  subroutine PDM_multipart_get_part_mesh_nodal_c (multipart, &
                                                  i_domain, &
                                                  pmesh_nodal, &
                                                  ownership) &
  bind (c, name='PDM_multipart_get_part_mesh_nodal')

    use iso_c_binding
    implicit none

    type(c_ptr),    value  :: multipart
    integer(c_int), value  :: i_domain
    type(c_ptr)            :: pmesh_nodal
    integer(c_int), value  :: ownership

  end subroutine PDM_multipart_get_part_mesh_nodal_c


  subroutine PDM_multipart_block_set_c (multipart, &
                                        i_domain, &
                                        dn_cell, &
                                        dn_face, &
                                        dn_vtx, &
                                        n_face_group, &
                                        dcell_face_idx, &
                                        dcell_face, &
                                        dface_cell, &
                                        dface_vtx_idx, &
                                        dface_vtx, &
                                        dvtx_coord, &
                                        dface_group_idx, &
                                        dface_group) &
  bind (c, name='PDM_multipart_block_set')

    use iso_c_binding
    implicit none

    type(c_ptr),    value  :: multipart
    integer(c_int), value  :: i_domain
    integer(c_int), value  :: dn_cell
    integer(c_int), value  :: dn_face
    integer(c_int), value  :: dn_vtx
    integer(c_int), value  :: n_face_group
    type(c_ptr),    value  :: dcell_face_idx
    type(c_ptr),    value  :: dcell_face
    type(c_ptr),    value  :: dface_cell
    type(c_ptr),    value  :: dface_vtx_idx
    type(c_ptr),    value  :: dface_vtx
    type(c_ptr),    value  :: dvtx_coord
    type(c_ptr),    value  :: dface_group_idx
    type(c_ptr),    value  :: dface_group

  end subroutine PDM_multipart_block_set_c


  subroutine PDM_multipart_part_dim_get_c (multipart, &
                                           i_domain, &
                                           i_part, &
                                           n_cell, &
                                           n_face, &
                                           n_face_part_bound, &
                                           n_vtx, &
                                           n_proc, &
                                           n_total_part, &
                                           s_cell_face, &
                                           s_face_vtx, &
                                           s_face_bound, &
                                           n_bound_groups) &
  bind (c, name='PDM_multipart_part_dim_get')

    use iso_c_binding
    implicit none

    type(c_ptr),    value  :: multipart
    integer(c_int), value  :: i_domain
    integer(c_int), value  :: i_part
    integer(c_int)         :: n_cell
    integer(c_int)         :: n_face
    integer(c_int)         :: n_face_part_bound
    integer(c_int)         :: n_vtx
    integer(c_int)         :: n_proc
    integer(c_int)         :: n_total_part
    integer(c_int)         :: s_cell_face
    integer(c_int)         :: s_face_vtx
    integer(c_int)         :: s_face_bound
    integer(c_int)         :: n_bound_groups

  end subroutine PDM_multipart_part_dim_get_c


  subroutine PDM_multipart_part_val_get_c (multipart, &
                                           i_domain, &
                                           i_part, &
                                           cell_face_idx, &
                                           cell_face, &
                                           cell_ln_to_gn, &
                                           face_cell, &
                                           face_vtx_idx, &
                                           face_vtx, &
                                           face_ln_to_gn, &
                                           face_part_bound_proc_idx, &
                                           face_part_bound_part_idx, &
                                           face_part_bound, &
                                           vtx, &
                                           vtx_ln_to_gn, &
                                           face_bound_idx, &
                                           face_bound, &
                                           face_bound_ln_to_gn) &
  bind (c, name='PDM_multipart_part_val_get')

    use iso_c_binding
    implicit none

    type(c_ptr),    value  :: multipart
    integer(c_int), value  :: i_domain
    integer(c_int), value  :: i_part
    type(c_ptr)            :: cell_face_idx
    type(c_ptr)            :: cell_face
    type(c_ptr)            :: cell_ln_to_gn
    type(c_ptr)            :: face_cell
    type(c_ptr)            :: face_vtx_idx
    type(c_ptr)            :: face_vtx
    type(c_ptr)            :: face_ln_to_gn
    type(c_ptr)            :: face_part_bound_proc_idx
    type(c_ptr)            :: face_part_bound_part_idx
    type(c_ptr)            :: face_part_bound
    type(c_ptr)            :: vtx
    type(c_ptr)            :: vtx_ln_to_gn
    type(c_ptr)            :: face_bound_idx
    type(c_ptr)            :: face_bound
    type(c_ptr)            :: face_bound_ln_to_gn

  end subroutine PDM_multipart_part_val_get_c


  function PDM_multipart_part_connectivity_get_c (multipart, &
                                                  i_domain, &
                                                  i_part, &
                                                  connectivity_type, &
                                                  connect_idx, &
                                                  connect, &
                                                  ownership) &
  result (pn_entity) &
  bind (c, name='PDM_multipart_part_connectivity_get')

    use iso_c_binding
    implicit none

    type(c_ptr),    value :: multipart
    integer(c_int), value :: i_domain
    integer(c_int), value :: i_part
    integer(c_int), value :: connectivity_type
    type(c_ptr)           :: connect
    type(c_ptr)           :: connect_idx
    integer(c_int), value :: ownership

    integer(c_int)        :: pn_entity

  end function PDM_multipart_part_connectivity_get_c


  function PDM_multipart_part_ln_to_gn_get_c (multipart, &
                                              i_domain, &
                                              i_part, &
                                              entity_type, &
                                              entity_ln_to_gn, &
                                              ownership) &
  result (pn_entity) &
  bind (c, name='PDM_multipart_part_ln_to_gn_get')

    use iso_c_binding
    implicit none

    type(c_ptr),    value  :: multipart
    integer(c_int), value  :: i_domain
    integer(c_int), value  :: i_part
    integer(c_int), value  :: entity_type
    type(c_ptr)            :: entity_ln_to_gn
    integer(c_int), value  :: ownership

    integer(c_int)         :: pn_entity

  end function PDM_multipart_part_ln_to_gn_get_c


  function PDM_multipart_partition_color_get_c (multipart, &
                                                i_domain, &
                                                i_part, &
                                                entity_type, &
                                                entity_color, &
                                                ownership) &
  result (pn_entity) &
  bind (c, name='PDM_multipart_partition_color_get')

    use iso_c_binding
    implicit none

    type(c_ptr),    value  :: multipart
    integer(c_int), value  :: i_domain
    integer(c_int), value  :: i_part
    integer(c_int), value  :: entity_type
    type(c_ptr)            :: entity_color
    integer(c_int), value  :: ownership

    integer(c_int)         :: pn_entity

  end function PDM_multipart_partition_color_get_c


  subroutine PDM_multipart_part_ghost_infomation_get_c (multipart, &
                                                        i_domain, &
                                                        i_part, &
                                                        vtx_ghost_information) &
  bind (c, name='PDM_multipart_part_ghost_infomation_get')

    use iso_c_binding
    implicit none

    type(c_ptr),    value  :: multipart
    integer(c_int), value  :: i_domain
    integer(c_int), value  :: i_part
    type(c_ptr)            :: vtx_ghost_information

  end subroutine PDM_multipart_part_ghost_infomation_get_c


  function PDM_multipart_part_vtx_coord_get_c (multipart, &
                                               i_domain, &
                                               i_part, &
                                               vtx_coord, &
                                               ownership) &
  result (n_vtx) &
  bind (c, name='PDM_multipart_part_vtx_coord_get')

    use iso_c_binding
    implicit none

    type(c_ptr),    value  :: multipart
    integer(c_int), value  :: i_domain
    integer(c_int), value  :: i_part
    type(c_ptr)            :: vtx_coord
    integer(c_int)         :: ownership
    integer(c_int)         :: n_vtx

  end function PDM_multipart_part_vtx_coord_get_c


  subroutine PDM_multipart_group_get_c (multipart, &
                                        i_domain, &
                                        i_part, &
                                        entity_type, &
                                        n_group, &
                                        group_entity_idx, &
                                        group_entity, &
                                        group_entity_ln_to_gn) &
  bind (c, name='PDM_multipart_group_get')

    use iso_c_binding
    implicit none

    type(c_ptr),    value  :: multipart
    integer(c_int), value  :: i_domain
    integer(c_int), value  :: i_part
    integer(c_int), value  :: entity_type
    integer(c_int)         :: n_group
    type(c_ptr)            :: group_entity_idx
    type(c_ptr)            :: group_entity
    type(c_ptr)            :: group_entity_ln_to_gn

  end subroutine PDM_multipart_group_get_c


  subroutine PDM_multipart_part_graph_comm_get_c (multipart,            &
                                                  i_domain,               &
                                                  i_part,               &
                                                  entity_type,          &
                                                  ppart_bound_proc_idx, &
                                                  ppart_bound_part_idx, &
                                                  ppart_bound,          &
                                                  ownership)            &
  bind (c, name='PDM_multipart_part_graph_comm_get')

    use iso_c_binding
    implicit none

    type(c_ptr),    value  :: multipart
    integer(c_int), value  :: i_domain
    integer(c_int), value  :: i_part
    integer(c_int), value  :: entity_type
    type(c_ptr)            :: ppart_bound_proc_idx
    type(c_ptr)            :: ppart_bound_part_idx
    type(c_ptr)            :: ppart_bound
    integer(c_int), value  :: ownership

  end subroutine PDM_multipart_part_graph_comm_get_c

  !>
  !!
  !! \brief Set distributed mesh data for the input domain
  !!
  !! \param [in]   multipart      Pointer to \p PDM_multipart_t object
  !! \param [in]   domain_id      Global domain id
  !! \param [in]   dmesh          Distributed mesh structure
  !!

  subroutine PDM_multipart_dmesh_set (multipart, &
                                      domain_id,   &
                                      dmesh)     &
  bind (c, name='PDM_multipart_dmesh_set')

    use iso_c_binding
    implicit none

    type(c_ptr),    value  :: multipart
    integer(c_int), value  :: domain_id
    type(c_ptr),    value  :: dmesh

  end subroutine PDM_multipart_dmesh_set

  !>
  !!
  !! \brief Set distributed mesh data for the input domain
  !!
  !! \param [in]   multipart      Pointer to \p PDM_multipart_t object
  !! \param [in]   domain_id      Global domain id
  !! \param [in]   dmesh_nodal    Distributed nodal mesh structure
  !!

  subroutine PDM_multipart_dmesh_nodal_set (multipart,   &
                                            domain_id,     &
                                            dmesh_nodal) &
  bind (c, name='PDM_multipart_dmesh_nodal_set')

    use iso_c_binding
    implicit none

    type(c_ptr),    value  :: multipart
    integer(c_int), value  :: domain_id
    type(c_ptr),    value  :: dmesh_nodal

  end subroutine PDM_multipart_dmesh_nodal_set

  !>
  !!
  !! \brief Construct the partitioned meshes on every domains
  !!
  !! \param [in]   multipart   Pointer to \p PDM_multipart_t object
  !!

  subroutine PDM_multipart_compute (multipart) &
  bind (c, name='PDM_multipart_compute')

    use iso_c_binding
    implicit none

    type(c_ptr), value  :: multipart

  end subroutine PDM_multipart_compute

  !>
  !!
  !! \brief Set the domain interface
  !!
  !! \param [in]   multipart             Pointer to \p PDM_multipart_t object
  !! \param [in]   ditrf                 Domain interface
  !!

  subroutine PDM_multipart_domain_interface_shared_set (multipart, &
                                                        ditrf) &
  bind (c, name='PDM_multipart_domain_interface_shared_set')

    use iso_c_binding
    implicit none

    type(c_ptr), value :: multipart
    type(c_ptr), value :: ditrf

  end subroutine PDM_multipart_domain_interface_shared_set


  subroutine PDM_multipart_free (multipart) &
  bind (c, name='PDM_multipart_free')
    ! Free the structure
    use iso_c_binding
    implicit none

    type(c_ptr), value     :: multipart ! Multipart instance

  end subroutine PDM_multipart_free

end interface

private :: PDM_multipart_create_,&
           PDM_multipart_set_reordering_options_,&
           PDM_multipart_set_reordering_options_vtx_,&
           PDM_multipart_get_part_mesh_nodal_,&
           PDM_multipart_block_set_,&
           PDM_multipart_part_dim_get_,&
           PDM_multipart_part_val_get_,&
           PDM_multipart_part_connectivity_get_,&
           PDM_multipart_part_ln_to_gn_get_,&
           PDM_multipart_partition_color_get_,&
           PDM_multipart_part_ghost_infomation_get_,&
           PDM_multipart_part_vtx_coord_get_,&
           PDM_multipart_group_get_,&
           PDM_multipart_part_graph_comm_get_

contains


  subroutine PDM_multipart_create_ (multipart, &
                                    n_domain, &
                                    n_part, &
                                    merge_domains, &
                                    split_method, &
                                    part_size_method, &
                                    part_fraction, &
                                    comm, &
                                    owner)
    ! Build a Multipart structure instance
    !
    ! Admissible values for ``split_method`` are:
    !   - ``PDM_SPLIT_DUAL_WITH_HILBERT``
    !   - ``PDM_SPLIT_DUAL_WITH_PARMETIS``
    !   - ``PDM_SPLIT_DUAL_WITH_PTSCOTCH``
    !
    ! Admissible values for ``part_size_method`` are:
    !   - ``PDM_PART_SIZE_HOMOGENEOUS``: All requested partition have the same size
    !   - ``PDM_PART_SIZE_HETEROGENEOUS``: Each requested partition can have a portion (within 0. and 1.) of the mesh
    use pdm
    use iso_c_binding
    implicit none

    type(c_ptr)                        :: multipart                     ! Pointer to a new PDM_multipart_t object
    integer(c_int),            value   :: n_domain                      ! Number of domains in the original mesh
    integer(kind=PDM_l_num_s), pointer :: n_part(:)                     ! Number of partition per proc in each domain (size = ``n_domain``)
    integer(c_int),            value   :: merge_domains                 ! Merge or not the domains before splitting
    integer(c_int),            value   :: split_method                  ! Choice of method used to split the mesh
    integer(c_int),            value   :: part_size_method              ! Choice of homogeneous or heterogeneous partitions
    real(8),                   pointer :: part_fraction(:)              ! Weight (in %) of each partition in heterogeneous case
    integer(c_int),            value   :: comm                          ! MPI communicator
    integer(c_int),            value   :: owner                         ! Data ownership
    type(c_ptr)                        :: c_n_part        
    type(c_ptr)                        :: c_part_fraction 
    integer(c_int)                     :: c_comm

    c_n_part = C_NULL_PTR
    if (associated(n_part)) then
      c_n_part        = c_loc(n_part)
    endif
      
    c_part_fraction = C_NULL_PTR
    if (associated(part_fraction)) then
      c_part_fraction = c_loc(part_fraction)
    endif
      
    c_comm = PDM_MPI_Comm_f2c(comm)

    multipart = PDM_multipart_create_c(n_domain, &
                                       c_n_part, &
                                       merge_domains, &
                                       split_method, &
                                       part_size_method, &
                                       c_part_fraction, &
                                       c_comm, &
                                       owner)

  end subroutine PDM_multipart_create_

  ! Set connecting data between all the zones

  ! subroutine PDM_multipart_register_joins_ (multipart, &
  !                                           n_total_joins, &
  !                                           join_to_opposite)

  !   use pdm
  !   use iso_c_binding
  !   implicit none

  !   type(c_ptr),               value   :: multipart                       ! Pointer to \ref PDM_multipart_t object
  !   integer(c_int),            value   :: n_total_joins                   ! Total number of interfaces
  !   integer(kind=PDM_l_num_s), pointer :: join_to_opposite(:)             ! For each global join id, give the global id of the opposite join (size = n_total_joins)
  !   type(c_ptr)                        :: c_join_to_opposite

  !   c_join_to_opposite = C_NULL_PTR
  !   if (associated(join_to_opposite)) then
  !     c_join_to_opposite = c_loc(join_to_opposite)
  !   endif  

  !   call PDM_multipart_register_joins_c(multipart, &
  !                                       n_total_joins, &
  !                                       c_join_to_opposite)

  ! end subroutine PDM_multipart_register_joins_

  ! Set the reordering methods to be used after partitioning

  subroutine PDM_multipart_set_reordering_options_ (multipart, &
                                                    i_domain, &
                                                    renum_cell_method, &
                                                    renum_cell_properties, &
                                                    renum_face_method)
    ! Set the reordering methods to be used after partitioning
    use pdm
    use iso_c_binding
    implicit none

    type(c_ptr),               value   :: multipart                            ! Multipart instance
    integer(c_int),            value   :: i_domain                             ! Id of domain which parameters apply (or -1 for all domains)
    character (len=*)                  :: renum_cell_method                    ! Choice of renumbering method for cells
    integer(kind=PDM_l_num_s), pointer :: renum_cell_properties(:)             ! Parameters used by cache-blocking method : [*n_cell_per_cache_wanted*, *is_asynchronous*, *is_vectorisation*, *n_vect_face*, *split_method*]
    character (len=*)                  :: renum_face_method                    ! Choice of renumbering method for faces
    type(c_ptr)                        :: c_renum_cell_properties

    c_renum_cell_properties = C_NULL_PTR
    if (associated(renum_cell_properties)) then
      c_renum_cell_properties = c_loc(renum_cell_properties)
    endif  
      
    call PDM_multipart_set_reordering_options_c(multipart, &
                                                i_domain, &
                                                trim(renum_cell_method)//C_NULL_CHAR, &
                                                c_renum_cell_properties, &
                                                trim(renum_face_method)//C_NULL_CHAR)

  end subroutine PDM_multipart_set_reordering_options_


  subroutine PDM_multipart_set_reordering_options_vtx_ (multipart, &
                                                        i_domain, &
                                                        renum_vtx_method)
    ! Set the reordering methods to be used after partitioning
    use iso_c_binding
    implicit none

    type(c_ptr),      value  :: multipart        ! Multipart instance
    integer(c_int),   value  :: i_domain         ! Id of domain which parameters apply (or -1 for all domains)
    character (len=*)        :: renum_vtx_method ! Choice of renumbering method for vertices

    call PDM_multipart_set_reordering_options_vtx_c(multipart, &
                                                    i_domain, &
                                                    trim(renum_vtx_method)//C_NULL_CHAR)

  end subroutine PDM_multipart_set_reordering_options_vtx_


  subroutine PDM_multipart_get_part_mesh_nodal_ (multipart, &
                                                 i_domain, &
                                                 pmesh_nodal, &
                                                 ownership)
    ! Retrieve the partitioned nodal mesh
    use iso_c_binding
    implicit none

    type(c_ptr),    value  :: multipart   ! Multipart instance
    integer(c_int), value  :: i_domain    ! Id of domain which parameters apply (or -1 for all domains)
    type(c_ptr)            :: pmesh_nodal ! Partitioned nodal mesh
    integer(c_int), value  :: ownership   ! Data ownership

    call PDM_multipart_get_part_mesh_nodal_c(multipart, &
                                             i_domain, &
                                             pmesh_nodal, &
                                             ownership)

  end subroutine PDM_multipart_get_part_mesh_nodal_


  subroutine PDM_multipart_block_set_ (multipart, &
                                       i_domain, &
                                       dn_cell, &
                                       dn_face, &
                                       dn_vtx, &
                                       n_face_group, &
                                       dcell_face_idx, &
                                       dcell_face, &
                                       dface_cell, &
                                       dface_vtx_idx, &
                                       dface_vtx, &
                                       dvtx_coord, &
                                       dface_group_idx, &
                                       dface_group)
    ! Set block data
    use iso_c_binding
    implicit none

    type(c_ptr)                        :: multipart          ! Multipart instance
    integer, intent(in)                :: i_domain           ! Id of domain
    integer, intent(in)                :: dn_cell            ! Number of distributed cells
    integer, intent(in)                :: dn_face            ! Number of distributed faces
    integer, intent(in)                :: dn_vtx             ! Number of distributed vertices
    integer, intent(in)                :: n_face_group       ! Number of face groups
    integer(kind=PDM_l_num_s), pointer :: dcell_face_idx(:)  ! Distributed cell face connectivity index or *null()* (size : ``dn_cell`` + 1, numbering : 0 to n-1)
    integer(kind=PDM_g_num_s), pointer :: dcell_face(:)      ! Distributed cell face connectivity or *null()* (size : ``dface_vtx_idx(dn_cell+1)``, numbering : 1 to n)
    integer(kind=PDM_g_num_s), pointer :: dface_cell(:)      ! Distributed face cell connectivity
    integer(kind=PDM_l_num_s), pointer :: dface_vtx_idx(:)   ! Distributed face to vertex connectivity index (size : ``dn_face`` + 1, numbering : 0 to n-1)
    integer(kind=PDM_g_num_s), pointer :: dface_vtx(:)       ! Distributed face to vertex connectivity (size : ``dface_vtx_idx(dn_face+1)``, numbering : 1 to n)
    real(8),                   pointer :: dvtx_coord(:,:)    ! Distributed vertex coordinates (shape = [3, ``dn_vtx``])
    integer(kind=PDM_l_num_s), pointer :: dface_group_idx(:) ! Index of distributed faces list of each group (size = ``n_face_group`` + 1) or *null()*
    integer(kind=PDM_g_num_s), pointer :: dface_group(:)     ! Distributed faces list of each group or *null()* (size = ``dface_group_idx(n_face_group+1)``, numbering : 1 to n)

    integer(kind=c_int)                :: c_i_domain
    integer(kind=c_int)                :: c_dn_cell
    integer(kind=c_int)                :: c_dn_face
    integer(kind=c_int)                :: c_dn_vtx
    integer(kind=c_int)                :: c_n_face_group
    type(c_ptr)                        :: c_dcell_face_idx
    type(c_ptr)                        :: c_dcell_face
    type(c_ptr)                        :: c_dface_cell
    type(c_ptr)                        :: c_dface_vtx_idx
    type(c_ptr)                        :: c_dface_vtx
    type(c_ptr)                        :: c_dvtx_coord
    type(c_ptr)                        :: c_dface_group_idx
    type(c_ptr)                        :: c_dface_group

    c_i_domain        = i_domain
    c_dn_cell         = dn_cell
    c_dn_face         = dn_face
    c_dn_vtx          = dn_vtx
    c_n_face_group    = n_face_group

    c_dface_vtx_idx = C_NULL_PTR
    if (associated(dface_vtx_idx)) then
      c_dface_vtx_idx   = c_loc(dface_vtx_idx  )
    endif
      
    c_dface_vtx = C_NULL_PTR
    if (associated(dface_vtx)) then
      c_dface_vtx       = c_loc(dface_vtx      )
    endif
      
    c_dvtx_coord = C_NULL_PTR
    if (associated(dvtx_coord)) then
      c_dvtx_coord      = c_loc(dvtx_coord     )
    endif
      
    c_dface_group_idx = C_NULL_PTR
    if (associated(dface_group_idx)) then
      c_dface_group_idx = c_loc(dface_group_idx)
    endif
      
    c_dface_group = C_NULL_PTR
    if (associated(dface_group)) then
      c_dface_group     = c_loc(dface_group    )
    endif  

    c_dcell_face_idx  = C_NULL_PTR
    c_dcell_face      = C_NULL_PTR
    c_dface_cell      = C_NULL_PTR

    if (associated(dcell_face_idx)) then
      c_dcell_face_idx = c_loc(dcell_face_idx)
    endif
    if (associated(dcell_face)) then
      c_dcell_face = c_loc(dcell_face)
    endif
    if (associated(dface_cell)) then
      c_dface_cell = c_loc(dface_cell)
    endif

    call PDM_multipart_block_set_c (multipart, &
                                    c_i_domain, &
                                    c_dn_cell, &
                                    c_dn_face, &
                                    c_dn_vtx, &
                                    c_n_face_group, &
                                    c_dcell_face_idx, &
                                    c_dcell_face, &
                                    c_dface_cell, &
                                    c_dface_vtx_idx, &
                                    c_dface_vtx, &
                                    c_dvtx_coord, &
                                    c_dface_group_idx, &
                                    c_dface_group)

  end subroutine PDM_multipart_block_set_


  subroutine PDM_multipart_part_dim_get_ (multipart, &
                                          i_domain, &
                                          i_part, &
                                          n_cell, &
                                          n_face, &
                                          n_face_part_bound, &
                                          n_vtx, &
                                          n_proc, &
                                          n_total_part, &
                                          s_cell_face, &
                                          s_face_vtx, &
                                          s_face_bound, &
                                          n_bound_groups)
    ! Returns the dimensions of a given partition
    use pdm
    use iso_c_binding
    implicit none

    type(c_ptr), value       :: multipart         ! Multipart instance
    integer,     intent(in)  :: i_domain          ! Domain identifier
    integer,     intent(in)  :: i_part            ! Partition identifier
    integer,     intent(out) :: n_cell            ! Number of cells
    integer,     intent(out) :: n_face            ! Number of faces
    integer,     intent(out) :: n_face_part_bound ! Number of boundary faces in partition
    integer,     intent(out) :: n_vtx             ! Number of vertices
    integer,     intent(out) :: n_proc            ! Number of processes
    integer,     intent(out) :: n_total_part      ! Total number of partitions
    integer,     intent(out) :: s_cell_face       ! Size of cell->face connectivity
    integer,     intent(out) :: s_face_vtx        ! Size of face->vtx connectivity
    integer,     intent(out) :: s_face_bound      ! Size of boundary faces
    integer,     intent(out) :: n_bound_groups    ! Number of boundary groups

    integer(c_int) :: c_n_cell
    integer(c_int) :: c_n_face
    integer(c_int) :: c_n_face_part_bound
    integer(c_int) :: c_n_vtx
    integer(c_int) :: c_n_proc
    integer(c_int) :: c_n_total_part
    integer(c_int) :: c_s_cell_face
    integer(c_int) :: c_s_face_vtx
    integer(c_int) :: c_s_face_bound
    integer(c_int) :: c_n_bound_groups

    call PDM_multipart_part_dim_get_c(multipart, &
                                      i_domain, &
                                      i_part, &
                                      c_n_cell, &
                                      c_n_face, &
                                      c_n_face_part_bound, &
                                      c_n_vtx, &
                                      c_n_proc, &
                                      c_n_total_part, &
                                      c_s_cell_face, &
                                      c_s_face_vtx, &
                                      c_s_face_bound, &
                                      c_n_bound_groups)

    n_cell            = c_n_cell
    n_face            = c_n_face
    n_face_part_bound = c_n_face_part_bound
    n_vtx             = c_n_vtx
    n_proc            = c_n_proc
    n_total_part      = c_n_total_part
    s_cell_face       = c_s_cell_face
    s_face_vtx        = c_s_face_vtx
    s_face_bound      = c_s_face_bound
    n_bound_groups    = c_n_bound_groups

  end subroutine PDM_multipart_part_dim_get_


  subroutine PDM_multipart_part_val_get_(multipart, &
                                         i_domain, &
                                         i_part, &
                                         cell_face_idx, &
                                         cell_face, &
                                         cell_ln_to_gn, &
                                         face_cell, &
                                         face_vtx_idx, &
                                         face_vtx, &
                                         face_ln_to_gn, &
                                         face_part_bound_proc_idx, &
                                         face_part_bound_part_idx, &
                                         face_part_bound, &
                                         vtx, &
                                         vtx_ln_to_gn, &
                                         face_bound_idx, &
                                         face_bound, &
                                         face_bound_ln_to_gn)
    ! Returns the data arrays of a given partition (Deprecated)
    use pdm
    use iso_c_binding
    use pdm_pointer_array
    implicit none

    type(c_ptr),    value  :: multipart
    integer(c_int), value  :: i_domain
    integer(c_int), value  :: i_part

    integer (kind = PDM_l_num_s), pointer :: cell_face_idx(:)
    integer (kind = PDM_l_num_s), pointer :: cell_face(:)
    integer (kind = PDM_g_num_s), pointer :: cell_ln_to_gn(:)
    integer (kind = PDM_l_num_s), pointer :: face_cell(:)
    integer (kind = PDM_l_num_s), pointer :: face_vtx_idx(:)
    integer (kind = PDM_l_num_s), pointer :: face_vtx(:)
    integer (kind = PDM_g_num_s), pointer :: face_ln_to_gn(:)
    integer (kind = PDM_l_num_s), pointer :: face_part_bound_proc_idx(:)
    integer (kind = PDM_l_num_s), pointer :: face_part_bound_part_idx(:)
    integer (kind = PDM_l_num_s), pointer :: face_part_bound(:)
    real(8),                      pointer :: vtx(:,:)
    integer (kind = PDM_g_num_s), pointer :: vtx_ln_to_gn(:)
    integer (kind = PDM_l_num_s), pointer :: face_bound_idx(:)
    integer (kind = PDM_l_num_s), pointer :: face_bound(:)
    integer (kind = PDM_g_num_s), pointer :: face_bound_ln_to_gn(:)

    type(c_ptr)  :: c_cell_face_idx
    type(c_ptr)  :: c_cell_face
    type(c_ptr)  :: c_cell_ln_to_gn
    type(c_ptr)  :: c_face_cell
    type(c_ptr)  :: c_face_vtx_idx
    type(c_ptr)  :: c_face_vtx
    type(c_ptr)  :: c_face_ln_to_gn
    type(c_ptr)  :: c_face_part_bound_proc_idx
    type(c_ptr)  :: c_face_part_bound_part_idx
    type(c_ptr)  :: c_face_part_bound
    type(c_ptr)  :: c_vtx
    type(c_ptr)  :: c_vtx_ln_to_gn
    type(c_ptr)  :: c_face_bound_idx
    type(c_ptr)  :: c_face_bound
    type(c_ptr)  :: c_face_bound_ln_to_gn

    integer(c_int) :: c_n_cell
    integer(c_int) :: c_n_face
    integer(c_int) :: c_n_face_part_bound
    integer(c_int) :: c_n_vtx
    integer(c_int) :: c_n_proc
    integer(c_int) :: c_n_total_part
    integer(c_int) :: c_s_cell_face
    integer(c_int) :: c_s_face_vtx
    integer(c_int) :: c_s_face_bound
    integer(c_int) :: c_n_bound_groups

    integer :: n_cell
    integer :: n_face
    integer :: n_face_part_bound
    integer :: n_vtx
    integer :: n_bound_groups

    call PDM_multipart_part_dim_get_c(multipart, &
                                      i_domain, &
                                      i_part, &
                                      c_n_cell, &
                                      c_n_face, &
                                      c_n_face_part_bound, &
                                      c_n_vtx, &
                                      c_n_proc, &
                                      c_n_total_part, &
                                      c_s_cell_face, &
                                      c_s_face_vtx, &
                                      c_s_face_bound, &
                                      c_n_bound_groups)

    n_cell             = c_n_cell
    n_face             = c_n_face
    n_face_part_bound  = c_n_face_part_bound
    n_vtx              = c_n_vtx
    n_bound_groups     = c_n_bound_groups
    c_cell_face_idx            = C_NULL_PTR
    c_cell_face                = C_NULL_PTR
    c_cell_ln_to_gn            = C_NULL_PTR
    c_face_cell                = C_NULL_PTR
    c_face_vtx_idx             = C_NULL_PTR
    c_face_vtx                 = C_NULL_PTR
    c_face_ln_to_gn            = C_NULL_PTR
    c_face_part_bound_proc_idx = C_NULL_PTR
    c_face_part_bound_part_idx = C_NULL_PTR
    c_face_part_bound          = C_NULL_PTR
    c_vtx                      = C_NULL_PTR
    c_vtx_ln_to_gn             = C_NULL_PTR
    c_face_bound_idx           = C_NULL_PTR
    c_face_bound               = C_NULL_PTR
    c_face_bound_ln_to_gn      = C_NULL_PTR

    call PDM_multipart_part_val_get_c(multipart, &
                                      i_domain, &
                                      i_part, &
                                      c_cell_face_idx, &
                                      c_cell_face, &
                                      c_cell_ln_to_gn, &
                                      c_face_cell, &
                                      c_face_vtx_idx, &
                                      c_face_vtx, &
                                      c_face_ln_to_gn, &
                                      c_face_part_bound_proc_idx, &
                                      c_face_part_bound_part_idx, &
                                      c_face_part_bound, &
                                      c_vtx, &
                                      c_vtx_ln_to_gn, &
                                      c_face_bound_idx, &
                                      c_face_bound, &
                                      c_face_bound_ln_to_gn)

    call c_f_pointer(c_cell_face_idx, &
                     cell_face_idx,   &
                     [n_cell + 1])

    call c_f_pointer(c_cell_face, &
                     cell_face,   &
                     [cell_face_idx(n_cell+1)])

    call c_f_pointer(c_cell_ln_to_gn, &
                     cell_ln_to_gn,   &
                     [n_cell])

    call c_f_pointer(c_face_cell, &
                     face_cell,   &
                     [2 * n_face])

    call c_f_pointer(c_face_vtx_idx, &
                     face_vtx_idx,   &
                     [n_face + 1])

    call c_f_pointer(c_face_vtx, &
                     face_vtx,   &
                     [face_vtx_idx(n_face+1)])

    call c_f_pointer(c_face_ln_to_gn, &
                     face_ln_to_gn,   &
                     [n_face])

    call c_f_pointer(c_face_part_bound_proc_idx, &
                     face_part_bound_proc_idx,   &
                     [n_face_part_bound + 1])

    call c_f_pointer(c_face_part_bound_part_idx, &
                     face_part_bound_part_idx,   &
                     [n_face_part_bound + 1])

    call c_f_pointer(c_face_part_bound, &
                     face_part_bound,   &
                     [face_part_bound_part_idx(n_face_part_bound+1)])

    call c_f_pointer(c_vtx, &
                     vtx,   &
                     [3, n_vtx])

    call c_f_pointer(c_vtx_ln_to_gn, &
                     vtx_ln_to_gn,   &
                     [n_vtx])

    call c_f_pointer(c_face_bound_idx, &
                     face_bound_idx,   &
                     [n_bound_groups + 1])

    call c_f_pointer(c_face_bound, &
                     face_bound,   &
                     [face_bound_idx(n_bound_groups+1)])

    call c_f_pointer(c_face_bound_ln_to_gn, &
                     face_bound_ln_to_gn,   &
                     [face_bound_idx(n_bound_groups+1)])

  end subroutine PDM_multipart_part_val_get_


  subroutine PDM_multipart_part_connectivity_get_ (multipart, &
                                                   i_domain, &
                                                   i_part, &
                                                   connectivity_type, &
                                                   connect_idx, &
                                                   connect, &
                                                   ownership, &
                                                   n_entity)
    ! Get a partitioned connectivity
    use pdm
    use iso_c_binding
    implicit none

    type(c_ptr),               value   :: multipart                  ! Multipart instance
    integer(c_int),            value   :: i_domain                   ! Domain identifier
    integer(c_int),            value   :: i_part                     ! Partition identifier
    integer(c_int),            value   :: connectivity_type          ! Connectivity type
    integer(kind=PDM_l_num_s), pointer :: connect_idx(:)             ! Connectivity index (size = ``n_entity`` + 1)
    integer(kind=PDM_l_num_s), pointer :: connect(:)                 ! Connectivity (size = ``connect_idx(p_entity+1)``)
    integer(c_int),            value   :: ownership                  ! Data ownership
    integer(c_int)                     :: n_entity                   ! Number of leading entities
    type(c_ptr)                        :: c_connect_idx 
    type(c_ptr)                        :: c_connect     
    integer(c_int)                     :: connec_size

    c_connect = C_NULL_PTR
    c_connect_idx = C_NULL_PTR
    n_entity = PDM_multipart_part_connectivity_get_c(multipart,         &
                                                     i_domain,          &
                                                     i_part,            &
                                                     connectivity_type, &
                                                     c_connect_idx,     &
                                                     c_connect,         &
                                                     ownership)

    connect_idx => null()
    if ( c_associated(c_connect_idx) ) then
      call c_f_pointer(c_connect_idx, &
                       connect_idx,   &
                       [n_entity+1])
    end if

    ! TO DO if c_connect_idx a C_NULL_PTR in other cases than edge_vtx ?
    connec_size = 2 * n_entity
    if ( c_associated(c_connect_idx) ) then
      connec_size = connect_idx(n_entity+1)
    end if

    call c_f_pointer(c_connect, &
                     connect,   &
                     [connec_size])

  end subroutine PDM_multipart_part_connectivity_get_


  subroutine PDM_multipart_part_ln_to_gn_get_ (multipart, &
                                               i_domain, &
                                               i_part, &
                                               entity_type, &
                                               entity_ln_to_gn, &
                                               ownership, &
                                               n_entity)
    ! Get the global ids of entities with given type
    use pdm
    use iso_c_binding
    implicit none

    type(c_ptr),               value   :: multipart                      ! Multipart instance
    integer(c_int),            value   :: i_domain                       ! Domain identifier
    integer(c_int),            value   :: i_part                         ! Partition identifier
    integer(c_int),            value   :: entity_type                    ! Entity type
    integer(kind=PDM_g_num_s), pointer :: entity_ln_to_gn(:)             ! Global ids (size = ``n_entity``)
    integer(c_int),            value   :: ownership                      ! Data ownership
    integer(c_int)                     :: n_entity                       ! Number of entities
    type(c_ptr)                        :: c_entity_ln_to_gn 

    c_entity_ln_to_gn = C_NULL_PTR
    n_entity = PDM_multipart_part_ln_to_gn_get_c(multipart, &
                                                 i_domain, &
                                                 i_part, &
                                                 entity_type, &
                                                 c_entity_ln_to_gn, &
                                                 ownership)

    call c_f_pointer(c_entity_ln_to_gn, &
                     entity_ln_to_gn,   &
                     [n_entity])

  end subroutine PDM_multipart_part_ln_to_gn_get_


  subroutine PDM_multipart_partition_color_get_(multipart, &
                                                i_domain, &
                                                i_part, &
                                                entity_type, &
                                                entity_color, &
                                                ownership, &
                                                n_entity)
    ! Get the color of entities with given type
    use pdm
    use iso_c_binding
    implicit none

    type(c_ptr),               value   :: multipart                   ! Multipart instance
    integer(c_int),            value   :: i_domain                    ! Domain identifier
    integer(c_int),            value   :: i_part                      ! Partition identifier
    integer(c_int),            value   :: entity_type                 ! Entity type
    integer(kind=PDM_l_num_s), pointer :: entity_color(:)             ! Entity color (size = ``n_entity``)
    integer(c_int),            value   :: ownership                   ! Data ownership
    integer(c_int)                     :: n_entity                   ! Number of entities
    type(c_ptr)                        :: c_entity_color

    c_entity_color = C_NULL_PTR
    n_entity = PDM_multipart_partition_color_get_c(multipart, &
                                                   i_domain, &
                                                   i_part, &
                                                   entity_type, &
                                                   c_entity_color, &
                                                   ownership)

    call c_f_pointer(c_entity_color, &
                     entity_color,   &
                     [n_entity])

  end subroutine PDM_multipart_partition_color_get_


  subroutine PDM_multipart_part_ghost_infomation_get_(multipart, &
                                                      i_domain, &
                                                      i_part, &
                                                      vtx_ghost_information)
    ! Get ghost vertex information
    use pdm
    use iso_c_binding
    implicit none

    type(c_ptr),    value              :: multipart                            ! Multipart instance
    integer(c_int), value              :: i_domain                             ! Domain identifier
    integer(c_int), value              :: i_part                               ! Partition identifier
    integer(kind=PDM_l_num_s), pointer :: vtx_ghost_information(:)             ! Integer that gives the current priority of vertices on current partitions
    type(c_ptr)                        :: c_vtx_ghost_information

    integer(c_int) :: c_n_cell
    integer(c_int) :: c_n_face
    integer(c_int) :: c_n_face_part_bound
    integer(c_int) :: c_n_vtx
    integer(c_int) :: c_n_proc
    integer(c_int) :: c_n_total_part
    integer(c_int) :: c_s_cell_face
    integer(c_int) :: c_s_face_vtx
    integer(c_int) :: c_s_face_bound
    integer(c_int) :: c_n_bound_groups


    c_vtx_ghost_information = C_NULL_PTR
    call PDM_multipart_part_dim_get_c(multipart, &
                                      i_domain, &
                                      i_part, &
                                      c_n_cell, &
                                      c_n_face, &
                                      c_n_face_part_bound, &
                                      c_n_vtx, &
                                      c_n_proc, &
                                      c_n_total_part, &
                                      c_s_cell_face, &
                                      c_s_face_vtx, &
                                      c_s_face_bound, &
                                      c_n_bound_groups)

    c_vtx_ghost_information = C_NULL_PTR
    call PDM_multipart_part_ghost_infomation_get_c(multipart, &
                                                   i_domain, &
                                                   i_part, &
                                                   c_vtx_ghost_information)

    call c_f_pointer(c_vtx_ghost_information, &
                     vtx_ghost_information,   &
                     [c_n_vtx])

  end subroutine PDM_multipart_part_ghost_infomation_get_


  subroutine PDM_multipart_part_vtx_coord_get_(multipart, &
                                               i_domain, &
                                               i_part, &
                                               vtx_coord, &
                                               ownership, &
                                               n_vtx)
    ! Get vertex coordinates
    use iso_c_binding
    implicit none

    type(c_ptr),      value     :: multipart                ! Multipart instance
    integer(c_int),   value     :: i_domain                 ! Domain identifier
    integer(c_int),   value     :: i_part                   ! Partition identifier
    real(8),          pointer   :: vtx_coord(:,:)           ! Vertex coordinates (shape = [3, ``n_vtx``])
    integer(c_int),   value     :: ownership                ! Data ownership
    integer(c_int)              :: n_vtx                    ! Number of vertices
    type(c_ptr)                 :: c_vtx_coord

    c_vtx_coord = C_NULL_PTR
    n_vtx = PDM_multipart_part_vtx_coord_get_c(multipart, &
                                               i_domain, &
                                               i_part, &
                                               c_vtx_coord, &
                                               ownership)

    call c_f_pointer(c_vtx_coord, &
                     vtx_coord,   &
                     [3, n_vtx])

  end subroutine PDM_multipart_part_vtx_coord_get_


  subroutine PDM_multipart_group_get_(multipart, &
                                      i_domain, &
                                      i_part, &
                                      entity_type, &
                                      n_group, &
                                      group_entity_idx, &
                                      group_entity, &
                                      group_entity_ln_to_gn)
    ! Get the group description for a given entity
    use pdm
    use iso_c_binding
    implicit none

    type(c_ptr),               value   :: multipart                     ! Multipart instance
    integer(c_int),            value   :: i_domain                      ! Domain identifier
    integer(c_int),            value   :: i_part                        ! Partition identifier
    integer(c_int)                     :: entity_type                   ! Type of mesh entity
    integer(c_int)                     :: n_group                       ! Number of groups
    integer(kind=PDM_l_num_s), pointer :: group_entity_idx(:)           ! Index for group->entity connectivity (size = ``n_group``)
    integer(kind=PDM_l_num_s), pointer :: group_entity(:)               ! Group->entity connectivity (1-based local ids, size = ``group_entity_idx(n_group+1)``)
    integer(kind=PDM_g_num_s), pointer :: group_entity_ln_to_gn(:)      ! Group->entity connectivity (group-specific global ids, size = ``group_entity_idx(n_group+1)``)
    type(c_ptr)                        :: c_group_entity
    type(c_ptr)                        :: c_group_entity_idx
    type(c_ptr)                        :: c_group_entity_ln_to_gn

    c_group_entity          = C_NULL_PTR
    c_group_entity_idx      = C_NULL_PTR
    c_group_entity_ln_to_gn = C_NULL_PTR

    call PDM_multipart_group_get_c(multipart, &
                                   i_domain, &
                                   i_part, &
                                   entity_type, &
                                   n_group, &
                                   c_group_entity_idx, &
                                   c_group_entity, &
                                   c_group_entity_ln_to_gn)

    call c_f_pointer(c_group_entity_idx, &
                     group_entity_idx,   &
                     [n_group + 1])

    call c_f_pointer(c_group_entity, &
                     group_entity,   &
                     [group_entity_idx(n_group+1)])

    call c_f_pointer(c_group_entity_ln_to_gn, &
                     group_entity_ln_to_gn,   &
                     [group_entity_idx(n_group+1)])

  end subroutine PDM_multipart_group_get_


  subroutine PDM_multipart_part_graph_comm_get_(multipart,            &
                                                i_domain,               &
                                                i_part,               &
                                                entity_type,          &
                                                ppart_bound_proc_idx, &
                                                ppart_bound_part_idx, &
                                                ppart_bound,          &
                                                ownership)
    ! Get the connection graph between partition for the requested entity type
    use pdm
    use iso_c_binding
    implicit none

    type(c_ptr),               value   :: multipart                     ! Multipart instance
    integer(c_int),            value   :: i_domain                      ! Domain identifier
    integer(c_int),            value   :: i_part                        ! Partition identifier
    integer(c_int)                     :: entity_type                   ! Type of mesh entity
    integer(kind=PDM_l_num_s), pointer :: ppart_bound_proc_idx(:)       ! Partitioning boundary entities index from process (size = n_proc + 1)
    integer(kind=PDM_l_num_s), pointer :: ppart_bound_part_idx(:)       ! Partitioning boundary entities index from partition (size = n_total_part + 1)
    integer(kind=PDM_l_num_s), pointer :: ppart_bound(:)                ! Partitioning boundary entities (size = 4 * n_entity_part_bound)
    integer(c_int)                     :: ownership                     ! Data ownership

    type(c_ptr)                        :: c_ppart_bound_proc_idx = C_NULL_PTR
    type(c_ptr)                        :: c_ppart_bound_part_idx = C_NULL_PTR
    type(c_ptr)                        :: c_ppart_bound          = C_NULL_PTR

    integer(c_int)                     :: n_cell
    integer(c_int)                     :: n_face
    integer(c_int)                     :: n_face_part_bound
    integer(c_int)                     :: n_vtx
    integer(c_int)                     :: n_proc
    integer(c_int)                     :: n_total_part
    integer(c_int)                     :: s_cell_face
    integer(c_int)                     :: s_face_vtx
    integer(c_int)                     :: s_face_bound
    integer(c_int)                     :: n_bound_groups

    call PDM_multipart_part_graph_comm_get_c(multipart,              &
                                             i_domain,               &
                                             i_part,                 &
                                             entity_type,            &
                                             c_ppart_bound_proc_idx, &
                                             c_ppart_bound_part_idx, &
                                             c_ppart_bound,          &
                                             ownership)

    call PDM_multipart_part_dim_get_ (multipart, &
                                      i_domain, &
                                      i_part, &
                                      n_cell, &
                                      n_face, &
                                      n_face_part_bound, &
                                      n_vtx, &
                                      n_proc, &
                                      n_total_part, &
                                      s_cell_face, &
                                      s_face_vtx, &
                                      s_face_bound, &
                                      n_bound_groups)

    call c_f_pointer(c_ppart_bound_proc_idx, &
                     ppart_bound_proc_idx,   &
                     [n_proc+1])

    call c_f_pointer(c_ppart_bound_part_idx  , &
                     ppart_bound_part_idx  ,   &
                     [n_total_part+1])

    call c_f_pointer(c_ppart_bound, &
                     ppart_bound,   &
                     [4 * ppart_bound_part_idx(n_total_part+1)])

  end subroutine PDM_multipart_part_graph_comm_get_

end module pdm_multipart
