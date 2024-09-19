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

module pdm_extract_part

  use pdm
  use iso_c_binding

  implicit none

  interface

  !>
  !!
  !! \brief Compute extraction
  !!
  !! \param [inout] extrp  PDM_extract_part_t instance
  !!
  !!

  subroutine PDM_extract_part_compute (extrp) &
  bind (c, name='PDM_extract_part_compute')
    use iso_c_binding
    implicit none

    type(c_ptr), value :: extrp

  end subroutine PDM_extract_part_compute

  !>
  !!
  !! \brief Free PDM_extract_part_t instance
  !!
  !! \param [inout] extrp  PDM_extract_part_t instance
  !!
  !! \return       NULL
  !!

  subroutine PDM_extract_part_free (extrp) &
  bind (c, name='PDM_extract_part_free')
    use iso_c_binding
    implicit none

    type(c_ptr), value :: extrp

  end subroutine PDM_extract_part_free

  !>
  !!
  !! \brief Free partially PDM_extract_part_t instance
  !!
  !! \param [inout] extrp  PDM_extract_part_t instance
  !!
  !! \return       NULL
  !!

  subroutine PDM_extract_part_partial_free (extrp) &
  bind (c, name='PDM_extract_part_partial_free')
    use iso_c_binding
    implicit none

    type(c_ptr), value :: extrp

  end subroutine PDM_extract_part_partial_free

  !>
  !!
  !! \brief Set PDM_part_mesh_nodal_elmts_t
  !!
  !! \param [inout] extrp  PDM_extract_part_t instance
  !! \param [in]    pmne   PDM_part_mesh_nodal_elmts_t instance
  !!
  !!

  subroutine PDM_extract_part_part_nodal_set (extrp, pmne) &
  bind (c, name='PDM_extract_part_part_nodal_set')
    use iso_c_binding
    implicit none

    type(c_ptr), value :: extrp
    type(c_ptr), value :: pmne

  end subroutine PDM_extract_part_part_nodal_set

  end interface

  contains

  !>
  !!
  !! \brief Build an extract_part struct
  !!
  !! \param [in]   dim                 Extraction dimension
  !! \param [in]   n_part_in           Number of initial partition
  !! \param [in]   n_part_out          Number of final partition
  !! \param [in]   extract_kind        Extraction kind : (local/requilibrate/from target)
  !! \param [in]   split_dual_method   Split method if requilibrate extract_kind
  !! \param [in]   compute_child_gnum  Yes/No computation of a newest global numbering
  !! \param [in]   ownership           Tell if you want ownership of resulting
  !! \param [in]   comm                MPI communicator
  !!
  !! \return   Initialized \ref PDM_extract_part_t instance
  !!

  subroutine PDM_extract_part_create (extrp,              &
                                      dim,                &
                                      n_part_in,          &
                                      n_part_out,         &
                                      extract_kind,       &
                                      split_dual_method,  &
                                      compute_child_gnum, &
                                      ownership,          &
                                      comm)
    use iso_c_binding
    implicit none

    type(c_ptr)                       :: extrp
    integer, intent(in)               :: dim
    integer, intent(in)               :: n_part_in
    integer, intent(in)               :: n_part_out
    integer, intent(in)               :: extract_kind
    integer, intent(in)               :: split_dual_method
    logical, intent(in)               :: compute_child_gnum
    integer, intent(in)               :: ownership
    integer, intent(in)               :: comm

    integer(c_int)                    :: c_compute_child_gnum
    integer(c_int)                    :: c_comm

    interface
      function pdm_extract_part_create_c (dim,                &
                                          n_part_in,          &
                                          n_part_out,         &
                                          extract_kind,       &
                                          split_dual_method,  &
                                          compute_child_gnum, &
                                          ownership,          &
                                          comm)               &
      result (extrp)                                          &
      bind(c, name='PDM_extract_part_create')
        use iso_c_binding
        implicit none

        type(c_ptr)           :: extrp
        integer(c_int), value :: dim
        integer(c_int), value :: n_part_in
        integer(c_int), value :: n_part_out
        integer(c_int), value :: extract_kind
        integer(c_int), value :: split_dual_method
        integer(c_int), value :: compute_child_gnum
        integer(c_int), value :: ownership
        integer(c_int), value :: comm

      end function pdm_extract_part_create_c
    end interface

    if (compute_child_gnum) then
      c_compute_child_gnum = 1
    else
      c_compute_child_gnum = 0
    endif

    c_comm = PDM_MPI_Comm_f2c(comm)

    extrp = pdm_extract_part_create_c (dim,                  &
                                       n_part_in,            &
                                       n_part_out,           &
                                       extract_kind,         &
                                       split_dual_method,    &
                                       c_compute_child_gnum, &
                                       ownership,            &
                                       c_comm)

  end subroutine PDM_extract_part_create

  !>
  !!
  !! \brief Set partition
  !!
  !! \param [in]   extrp             PDM_extract_part_t
  !! \param [in]   i_part            part identifier
  !!
  !!

  subroutine PDM_extract_part_part_set (extrp,                 &
                                        i_part,                &
                                        n_cell,                &
                                        n_face,                &
                                        n_edge,                &
                                        n_vtx,                 &
                                        cell_face_idx,         &
                                        cell_face,             &
                                        face_edge_idx,         &
                                        face_edge,             &
                                        edge_vtx,              &
                                        face_vtx_idx,          &
                                        face_vtx,              &
                                        cell_ln_to_gn,         &
                                        face_ln_to_gn,         &
                                        edge_ln_to_gn,         &
                                        vtx_ln_to_gn,          &
                                        vtx_coord)

    use iso_c_binding
    implicit none

    type(c_ptr), value            :: extrp
    integer, intent(in)           :: i_part
    integer, intent(in)           :: n_cell
    integer, intent(in)           :: n_face
    integer, intent(in)           :: n_edge
    integer, intent(in)           :: n_vtx
    integer(pdm_l_num_s), pointer :: cell_face_idx(:)
    integer(pdm_l_num_s), pointer :: cell_face(:)
    integer(pdm_l_num_s), pointer :: face_edge_idx(:)
    integer(pdm_l_num_s), pointer :: face_edge(:)
    integer(pdm_l_num_s), pointer :: edge_vtx(:)
    integer(pdm_l_num_s), pointer :: face_vtx_idx(:)
    integer(pdm_l_num_s), pointer :: face_vtx(:)
    integer(pdm_g_num_s), pointer :: cell_ln_to_gn(:)
    integer(pdm_g_num_s), pointer :: face_ln_to_gn(:)
    integer(pdm_g_num_s), pointer :: edge_ln_to_gn(:)
    integer(pdm_g_num_s), pointer :: vtx_ln_to_gn(:)
    double precision,     pointer :: vtx_coord(:,:)

    interface
      subroutine pdm_extract_part_part_set_c (extrp,                    &
                                              i_part,                   &
                                              n_cell,                   &
                                              n_face,                   &
                                              n_edge,                   &
                                              n_vtx,                    &
                                              cell_face_idx,            &
                                              cell_face,                &
                                              face_edge_idx,            &
                                              face_edge,                &
                                              edge_vtx,                 &
                                              face_vtx_idx,             &
                                              face_vtx,                 &
                                              cell_ln_to_gn,            &
                                              face_ln_to_gn,            &
                                              edge_ln_to_gn,            &
                                              vtx_ln_to_gn,             &
                                              vtx_coord)                &
      bind (c, name='PDM_extract_part_part_set')
        use iso_c_binding
        implicit none

        type(c_ptr),    value :: extrp
        integer(c_int), value :: i_part
        integer(c_int), value :: n_cell
        integer(c_int), value :: n_face
        integer(c_int), value :: n_edge
        integer(c_int), value :: n_vtx
        type(c_ptr),    value :: cell_face_idx
        type(c_ptr),    value :: cell_face
        type(c_ptr),    value :: face_edge_idx
        type(c_ptr),    value :: face_edge
        type(c_ptr),    value :: edge_vtx
        type(c_ptr),    value :: face_vtx_idx
        type(c_ptr),    value :: face_vtx
        type(c_ptr),    value :: cell_ln_to_gn
        type(c_ptr),    value :: face_ln_to_gn
        type(c_ptr),    value :: edge_ln_to_gn
        type(c_ptr),    value :: vtx_ln_to_gn
        type(c_ptr),    value :: vtx_coord

      end subroutine pdm_extract_part_part_set_c
    end interface

    call pdm_extract_part_part_set_c (extrp,                           &
                                      i_part,                          &
                                      n_cell,                          &
                                      n_face,                          &
                                      n_edge,                          &
                                      n_vtx,                           &
                                      c_loc(cell_face_idx),            &
                                      c_loc(cell_face),                &
                                      c_loc(face_edge_idx),            &
                                      c_loc(face_edge),                &
                                      c_loc(edge_vtx),                 &
                                      c_loc(face_vtx_idx),             &
                                      c_loc(face_vtx),                 &
                                      c_loc(cell_ln_to_gn),            &
                                      c_loc(face_ln_to_gn),            &
                                      c_loc(edge_ln_to_gn),            &
                                      c_loc(vtx_ln_to_gn),             &
                                      c_loc(vtx_coord))

  end subroutine PDM_extract_part_part_set

  !>
  !!
  !! \brief Set partition group (optional)
  !!
  !! \param [in]   extrp                      PDM_extract_part_t
  !! \param [in]   bound_type                 Kind of group
  !! \param [in]   n_group                    Number of groups
  !!
  !!

  subroutine PDM_extract_part_n_group_set (extrp,      &
                                           bound_type, &
                                           n_group)

    use iso_c_binding
    implicit none

    type(c_ptr), value            :: extrp
    integer, intent(in)           :: bound_type
    integer, intent(in)           :: n_group

    interface
      subroutine pdm_extract_part_n_group_set_c (extrp,      &
                                                 bound_type, &
                                                 n_group)    &
      bind (c, name='PDM_extract_part_n_group_set')
        use iso_c_binding
        implicit none

        type(c_ptr),    value :: extrp
        integer(c_int), value :: bound_type
        integer(c_int), value :: n_group

      end subroutine pdm_extract_part_n_group_set_c
    end interface

    call pdm_extract_part_n_group_set_c (extrp,      &
                                         bound_type, &
                                         n_group)

  end subroutine PDM_extract_part_n_group_set

  !>
  !!
  !! \brief Set partition group (optional)
  !!
  !! \param [in]   extrp                      PDM_extract_part_t
  !! \param [in]   i_part                     part identifier
  !! \param [in]   i_group                    group identifier
  !! \param [in]   bound_type                 Kind of group
  !! \param [in]   n_group_entity             Number of entity in current group
  !! \param [in]   group_entity               List of entity in group (size = n_group_entity)
  !! \param [in]   group_entity_ln_to_gn      Global numbering of entity in group (size = n_group_entity)
  !!
  !!

  subroutine PDM_extract_part_part_group_set (extrp,                 &
                                              i_part,                &
                                              i_group,               &
                                              bound_type,            &
                                              n_group_entity,        &
                                              group_entity,          &
                                              group_entity_ln_to_gn)

    use iso_c_binding
    implicit none

    type(c_ptr), value            :: extrp
    integer, intent(in)           :: i_part
    integer, intent(in)           :: i_group
    integer, intent(in)           :: bound_type
    integer, intent(in)           :: n_group_entity
    integer(pdm_l_num_s), pointer :: group_entity(:)
    integer(pdm_g_num_s), pointer :: group_entity_ln_to_gn(:)

    interface
      subroutine pdm_extract_part_part_group_set_c (extrp,                 &
                                                    i_part,                &
                                                    i_group,               &
                                                    bound_type,            &
                                                    n_group_entity,        &
                                                    group_entity,          &
                                                    group_entity_ln_to_gn) &
      bind (c, name='PDM_extract_part_part_group_set')
        use iso_c_binding
        implicit none

        type(c_ptr),    value :: extrp
        integer(c_int), value :: i_part
        integer(c_int), value :: i_group
        integer(c_int), value :: bound_type
        integer(c_int), value :: n_group_entity
        type(c_ptr),    value :: group_entity
        type(c_ptr),    value :: group_entity_ln_to_gn

      end subroutine pdm_extract_part_part_group_set_c
    end interface

    call pdm_extract_part_part_group_set_c (extrp,                        &
                                            i_part,                       &
                                            i_group,                      &
                                            bound_type,                   &
                                            n_group_entity,               &
                                            c_loc(group_entity),          &
                                            c_loc(group_entity_ln_to_gn))

  end subroutine PDM_extract_part_part_group_set

  !>
  !!
  !! \brief Get partition group (optional)
  !!
  !! \param [in]   extrp                                     PDM_extract_part_t
  !! \param [in]   bound_type                                Kind of group
  !! \param [in]   i_part                                    part identifier
  !! \param [in]   i_group                                   group identifier
  !! \param [out]  n_extract_group_entity                    Number of entity in current group
  !! \param [out]  extract_group_entity                      List of entity in group (size = n_extract_group_entity)
  !! \param [out]  extract_group_entity_ln_to_gn             Global numbering of entity in group (size = n_extract_group_entity)
  !! \param [out]  extract_group_entity_parent_ln_to_gn      Global numbering of entity in group (size = n_extract_group_entity)
  !! \param [in]   ownership                                 Ownership
  !!
  !!

  subroutine PDM_extract_part_group_get (extrp,                                &
                                         bound_type,                           &
                                         i_part,                               &
                                         i_group,                              &
                                         n_extract_group_entity,               &
                                         extract_group_entity,                 &
                                         extract_group_entity_ln_to_gn,        &
                                         extract_group_entity_parent_ln_to_gn, &
                                         ownership)

    use iso_c_binding
    implicit none

    type(c_ptr), value            :: extrp
    integer, intent(in)           :: bound_type
    integer, intent(in)           :: i_part
    integer, intent(in)           :: i_group
    integer                       :: n_extract_group_entity
    integer(pdm_l_num_s), pointer :: extract_group_entity(:)
    integer(pdm_g_num_s), pointer :: extract_group_entity_ln_to_gn(:)
    integer(pdm_g_num_s), pointer :: extract_group_entity_parent_ln_to_gn(:)
    integer, intent(in)           :: ownership

    type(c_ptr)                   :: c_extract_group_entity                 
    type(c_ptr)                   :: c_extract_group_entity_ln_to_gn        
    type(c_ptr)                   :: c_extract_group_entity_parent_ln_to_gn 

    interface
      subroutine pdm_extract_part_group_get_c (extrp,                                &
                                               bound_type,                           &
                                               i_part,                               &
                                               i_group,                              &
                                               n_extract_group_entity,               &
                                               extract_group_entity,                 &
                                               extract_group_entity_ln_to_gn,        &
                                               extract_group_entity_parent_ln_to_gn, &
                                               ownership)                            &
      bind (c, name='PDM_extract_part_group_get')
        use iso_c_binding
        implicit none

        type(c_ptr),    value :: extrp
        integer(c_int), value :: bound_type
        integer(c_int), value :: i_part
        integer(c_int), value :: i_group
        integer(c_int)        :: n_extract_group_entity
        type(c_ptr)           :: extract_group_entity
        type(c_ptr)           :: extract_group_entity_ln_to_gn
        type(c_ptr)           :: extract_group_entity_parent_ln_to_gn
        integer(c_int), value :: ownership

      end subroutine pdm_extract_part_group_get_c
    end interface

    c_extract_group_entity                 = C_NULL_PTR
    c_extract_group_entity_ln_to_gn        = C_NULL_PTR
    c_extract_group_entity_parent_ln_to_gn = C_NULL_PTR

    call pdm_extract_part_group_get_c (extrp,                                  &
                                       bound_type,                             &
                                       i_part,                                 &
                                       i_group,                                &
                                       n_extract_group_entity,                 &
                                       c_extract_group_entity,                 &
                                       c_extract_group_entity_ln_to_gn,        &
                                       c_extract_group_entity_parent_ln_to_gn, &
                                       ownership)

    call c_f_pointer(c_extract_group_entity, &
                     extract_group_entity,   &
                     [n_extract_group_entity])

    call c_f_pointer(c_extract_group_entity_ln_to_gn, &
                     extract_group_entity_ln_to_gn,   &
                     [n_extract_group_entity])

    call c_f_pointer(c_extract_group_entity_parent_ln_to_gn, &
                     extract_group_entity_parent_ln_to_gn,   &
                     [n_extract_group_entity])

  end subroutine PDM_extract_part_group_get

  !>
  !!
  !! \brief Set the extract number
  !!
  !! \param [in]   extrp         PDM_extract_part_t
  !! \param [in]   i_part        part identifier
  !! \param [in]   n_extract     Number of entity to select
  !! \param [in]   extract_lnum  List of id to extract (starting at 0)
  !!
  !!

  subroutine PDM_extract_part_selected_lnum_set (extrp,     &
                                                 i_part_in, &
                                                 n_entity,  &
                                                 extract_lnum)

    use iso_c_binding
    implicit none

    type(c_ptr), value                   :: extrp
    integer, intent(in)                  :: i_part_in
    integer, intent(in)                  :: n_entity
    integer(kind = PDM_l_num_s), pointer :: extract_lnum(:)

    interface
      subroutine pdm_extract_part_selected_lnum_set_c (extrp,        &
                                                       i_part_in,    &
                                                       n_entity,     &
                                                       extract_lnum) &
      bind (c, name='PDM_extract_part_selected_lnum_set')
        use iso_c_binding
        implicit none

        type(c_ptr),    value :: extrp
        integer(c_int), value :: i_part_in
        integer(c_int), value :: n_entity
        type(c_ptr),    value :: extract_lnum

      end subroutine pdm_extract_part_selected_lnum_set_c
    end interface

    call pdm_extract_part_selected_lnum_set_c (extrp,             &
                                               i_part_in,         &
                                               n_entity,          &
                                               c_loc(extract_lnum))

  end subroutine PDM_extract_part_selected_lnum_set

  !>
  !!
  !! \brief Get number of entity in extraction
  !!
  !! \param [in]  extrp         PDM_extract_part_t instance
  !! \param [in]  i_part_out    Number of final partition
  !! \param [in]  entity_type   Type of entity
  !! \param [out] n_entity      Number of entity
  !!
  !!

  subroutine PDM_extract_part_n_entity_get (extrp,       &
                                            i_part_out,  &
                                            entity_type, &
                                            n_entity)

    use iso_c_binding
    implicit none

    type(c_ptr), value                :: extrp
    integer, intent(in)               :: i_part_out
    integer, intent(in)               :: entity_type
    integer                           :: n_entity

    interface
      function pdm_extract_part_n_entity_get_c (extrp,       &
                                                i_part_out,  &
                                                entity_type) &
      result (n_entity)                                      &
      bind (c, name='PDM_extract_part_n_entity_get')
        use iso_c_binding
        implicit none

        type(c_ptr),    value :: extrp
        integer(c_int), value :: i_part_out
        integer(c_int), value :: entity_type
        integer(c_int)        :: n_entity

      end function pdm_extract_part_n_entity_get_c
    end interface

    n_entity = pdm_extract_part_n_entity_get_c (extrp,       &
                                                i_part_out,  &
                                                entity_type)

  end subroutine PDM_extract_part_n_entity_get

  !>
  !!
  !! \brief Get connectivity of entity in extraction
  !!
  !! \param [in]  extrp               PDM_extract_part_t instance
  !! \param [in]  i_part_out          Number of final partition
  !! \param [in]  connectivity_type   Type of connectivity
  !! \param [out] n_entity            Number of entity in connectivity
  !! \param [out] connect             Entity connectivity
  !! \param [out] connect_idx         Entity connectivity index
  !! \param [in]  ownership           Tell if you want ownership of resulting
  !!
  !!

  subroutine PDM_extract_part_connectivity_get (extrp,             &
                                                i_part_out,        &
                                                connectivity_type, &
                                                n_entity,          &
                                                connect,           &
                                                connect_idx,       &
                                                ownership)

    use iso_c_binding
    implicit none

    type(c_ptr), value                   :: extrp
    integer, intent(in)                  :: i_part_out
    integer, intent(in)                  :: connectivity_type
    integer, intent(in)                  :: ownership
    integer                              :: n_entity
    integer(kind = PDM_l_num_s), pointer :: connect(:)
    integer(kind = PDM_l_num_s), pointer :: connect_idx(:)

    type(c_ptr)                          :: c_connect    
    type(c_ptr)                          :: c_connect_idx

    interface
      function pdm_extract_part_connectivity_get_c (extrp,             &
                                                    i_part_out,        &
                                                    connectivity_type, &
                                                    connect,           &
                                                    connect_idx,       &
                                                    ownership)         &
      result (n_entity)                                                &
      bind (c, name='PDM_extract_part_connectivity_get')
        use iso_c_binding
        implicit none

        type(c_ptr),    value :: extrp
        integer(c_int), value :: i_part_out
        integer(c_int), value :: connectivity_type
        integer(c_int), value :: ownership
        integer(c_int)        :: n_entity
        type(c_ptr)           :: connect
        type(c_ptr)           :: connect_idx

      end function pdm_extract_part_connectivity_get_c
    end interface

    c_connect     = C_NULL_PTR
    c_connect_idx = C_NULL_PTR

    n_entity = pdm_extract_part_connectivity_get_c (extrp,             &
                                                    i_part_out,        &
                                                    connectivity_type, &
                                                    c_connect,         &
                                                    c_connect_idx,     &
                                                    ownership)

    call c_f_pointer(c_connect_idx, &
                     connect_idx,   &
                     [n_entity+1])

    call c_f_pointer(c_connect, &
                     connect,   &
                     [connect_idx(n_entity+1)])

  end subroutine PDM_extract_part_connectivity_get

  !>
  !!
  !! \brief Get global numbering of entity in extraction
  !!
  !! \param [in]  extrp               PDM_extract_part_t instance
  !! \param [in]  i_part_out          Number of final partition
  !! \param [in]  entity_type         Type of entity
  !! \param [out] n_entity            Number of entity
  !! \param [out] pentity_ln_to_gn    Entity global numbering
  !! \param [in]  ownership           Tell if you want ownership of resulting
  !!
  !!

  subroutine PDM_extract_part_ln_to_gn_get (extrp,            &
                                            i_part_out,       &
                                            entity_type,      &
                                            n_entity,         &
                                            pentity_ln_to_gn, &
                                            ownership)

    use iso_c_binding
    implicit none

    type(c_ptr), value                   :: extrp
    integer, intent(in)                  :: i_part_out
    integer, intent(in)                  :: entity_type
    integer, intent(in)                  :: ownership
    integer                              :: n_entity
    integer(kind = PDM_g_num_s), pointer :: pentity_ln_to_gn(:)

    type(c_ptr)                          :: c_pentity_ln_to_gn

    interface
      function pdm_extract_part_ln_to_gn_get_c (extrp,            &
                                                i_part_out,       &
                                                entity_type,      &
                                                pentity_ln_to_gn, &
                                                ownership)        &
      result (n_entity)                                           &
      bind (c, name='PDM_extract_part_ln_to_gn_get')
        use iso_c_binding
        implicit none

        type(c_ptr),    value :: extrp
        integer(c_int), value :: i_part_out
        integer(c_int), value :: entity_type
        integer(c_int), value :: ownership
        integer(c_int)        :: n_entity
        type(c_ptr)           :: pentity_ln_to_gn

      end function pdm_extract_part_ln_to_gn_get_c
    end interface

    c_pentity_ln_to_gn = C_NULL_PTR

    n_entity = pdm_extract_part_ln_to_gn_get_c (extrp,              &
                                                i_part_out,         &
                                                entity_type,        &
                                                c_pentity_ln_to_gn, &
                                                ownership)

    call c_f_pointer(c_pentity_ln_to_gn, &
                     pentity_ln_to_gn,   &
                     [n_entity])

  end subroutine PDM_extract_part_ln_to_gn_get


  !>
  !!
  !! \brief Get Parent global numbering of entity 
  !!
  !! \param [in]  extrp               PDM_extract_part_t instance
  !! \param [in]  i_part_out          Number of final partition
  !! \param [in]  entity_type         Type of entity
  !! \param [out] n_entity            Number of entity
  !! \param [out] parent_ln_to_gn     Entity global numbering
  !! \param [in]  ownership           Tell if you want ownership of resulting
  !!
  !!

  subroutine PDM_extract_part_parent_ln_to_gn_get (extrp,            &
                                                   i_part_out,       &
                                                   entity_type,      &
                                                   n_entity,         &
                                                   parent_ln_to_gn, &
                                                   ownership)

    use iso_c_binding
    implicit none

    type(c_ptr), value                   :: extrp
    integer, intent(in)                  :: i_part_out
    integer, intent(in)                  :: entity_type
    integer, intent(in)                  :: ownership
    integer                              :: n_entity
    integer(kind = PDM_g_num_s), pointer :: parent_ln_to_gn(:)

    type(c_ptr)                          :: c_parent_ln_to_gn = C_NULL_PTR

    interface
      function pdm_extract_part_parent_ln_to_gn_get_c (extrp,            &
                                                       i_part_out,       &
                                                       entity_type,      &
                                                       parent_ln_to_gn, &
                                                       ownership)        &
      result (n_entity)                                                  &
      bind (c, name='PDM_extract_part_parent_ln_to_gn_get')
        use iso_c_binding
        implicit none

        type(c_ptr),    value :: extrp
        integer(c_int), value :: i_part_out
        integer(c_int), value :: entity_type
        integer(c_int), value :: ownership
        integer(c_int)        :: n_entity
        type(c_ptr)           :: parent_ln_to_gn

      end function pdm_extract_part_parent_ln_to_gn_get_c
    end interface

    n_entity = pdm_extract_part_parent_ln_to_gn_get_c (extrp,              &
                                                       i_part_out,         &
                                                       entity_type,        &
                                                       c_parent_ln_to_gn,  &
                                                       ownership)

    call c_f_pointer(c_parent_ln_to_gn, &
                     parent_ln_to_gn,   &
                     [n_entity])

  end subroutine PDM_extract_part_parent_ln_to_gn_get

  !>
  !!
  !! \brief Get local numbering of parent entity
  !!
  !! \param [in]  extrp               PDM_extract_part_t instance
  !! \param [in]  i_part_out          Number of final partition
  !! \param [in]  entity_type         Type of entity
  !! \param [out] n_entity            Number of entity
  !! \param [out] parent_entity_lnum  Local numbering of the parent entity
  !! \param [in]  ownership           Tell if you want ownership of resulting
  !!
  !!

  subroutine PDM_extract_part_parent_lnum_get (extrp,              &
                                               i_part_out,         &
                                               entity_type,        &
                                               n_entity,           &
                                               parent_entity_lnum, &
                                               ownership)

    use iso_c_binding
    implicit none

    type(c_ptr), value                   :: extrp
    integer, intent(in)                  :: i_part_out
    integer, intent(in)                  :: entity_type
    integer, intent(in)                  :: ownership
    integer                              :: n_entity
    integer(kind = PDM_l_num_s), pointer :: parent_entity_lnum(:)

    type(c_ptr)                          :: c_parent_entity_lnum

    interface
      function pdm_extract_part_parent_lnum_get_c (extrp,              &
                                                   i_part_out,         &
                                                   entity_type,        &
                                                   parent_entity_lnum, &
                                                   ownership)          &
      result (n_entity)                                                &
      bind (c, name='PDM_extract_part_parent_lnum_get')
        use iso_c_binding
        implicit none

        type(c_ptr),    value :: extrp
        integer(c_int), value :: i_part_out
        integer(c_int), value :: entity_type
        integer(c_int), value :: ownership
        integer(c_int)        :: n_entity
        type(c_ptr)           :: parent_entity_lnum

      end function pdm_extract_part_parent_lnum_get_c
    end interface

    c_parent_entity_lnum = C_NULL_PTR

    n_entity = pdm_extract_part_parent_lnum_get_c (extrp,                &
                                                   i_part_out,           &
                                                   entity_type,          &
                                                   c_parent_entity_lnum, &
                                                   ownership)

    call c_f_pointer(c_parent_entity_lnum, &
                     parent_entity_lnum,   &
                     [n_entity])

  end subroutine PDM_extract_part_parent_lnum_get

  subroutine PDM_extract_part_part_to_part_get (extrp,       &
                                                entity_type, &
                                                ptp,         &
                                                ownership)
    use iso_c_binding
    implicit none

    type(c_ptr)                          :: extrp
    integer, intent(in)                  :: entity_type
    type(c_ptr)                          :: ptp
    integer, intent(in)                  :: ownership


    interface

      subroutine PDM_extract_part_part_to_part_get_c (extrp,       &
                                                       entity_type, &
                                                       ptp,         &
                                                       ownership)   &
        bind (c, name='PDM_extract_part_part_to_part_get')
        use iso_c_binding
        implicit none
    
        type(c_ptr), value           :: extrp
        integer(c_int), value        :: entity_type
        type(c_ptr)                  :: ptp
        integer(c_int), value        :: ownership
      end subroutine PDM_extract_part_part_to_part_get_c

    end interface

    integer(c_int) :: c_entity_type
    integer(c_int) :: c_ownership

    c_entity_type = entity_type 
    c_ownership   = ownership

    call PDM_extract_part_part_to_part_get_c (extrp, c_entity_type, ptp, c_ownership) 

   end subroutine PDM_extract_part_part_to_part_get

  subroutine PDM_extract_part_vtx_coord_get (extrp,            &
                                             i_part_out,       &
                                             n_entity,         &
                                             pvtx_coord,       &
                                             ownership)

    use iso_c_binding
    implicit none

    type(c_ptr), value                   :: extrp
    integer, intent(in)                  :: i_part_out
    integer, intent(in)                  :: ownership
    integer                              :: n_entity
    double precision, pointer            :: pvtx_coord(:,:)

    type(c_ptr)                          :: c_pvtx_coord = C_NULL_PTR

    interface
      function pdm_extract_part_vtx_coord_get_c (extrp,            &
                                                 i_part_out,       &
                                                 pvtx_coord,       &
                                                 ownership)        &
      result (n_entity)                                            &
      bind (c, name='PDM_extract_part_vtx_coord_get')
        use iso_c_binding
        implicit none

        type(c_ptr),    value :: extrp
        integer(c_int), value :: i_part_out
        integer(c_int), value :: ownership
        integer(c_int)        :: n_entity
        type(c_ptr)           :: pvtx_coord

      end function pdm_extract_part_vtx_coord_get_c
    end interface

    n_entity = pdm_extract_part_vtx_coord_get_c (extrp,              &
                                                 i_part_out,         &
                                                 c_pvtx_coord,       &
                                                 ownership)

    call c_f_pointer(c_pvtx_coord, &
                     pvtx_coord,   &
                     [3,n_entity])

  end subroutine PDM_extract_part_vtx_coord_get

end module pdm_extract_part
