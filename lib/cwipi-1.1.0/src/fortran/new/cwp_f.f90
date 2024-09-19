!>
!! \file
!!


!-----------------------------------------------------------------------------
! This file is part of the CWIPI library.
!
! Copyright (C) 2021-2023  ONERA
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

module cwp
    use iso_c_binding
    use pdm_fortran

    ! CWP_Type_t
    enum, bind(c)
        enumerator :: &
                CWP_DOUBLE, &
                CWP_INT, &
                CWP_CHAR
    end enum

    ! CWP_Visu_format_t
    enum, bind(c)
        enumerator :: &
                CWP_VISU_FORMAT_ENSIGHT
    end enum

    ! CWP_Comm_t
    enum, bind(c)
        enumerator :: &
                CWP_COMM_PAR_WITH_PART, &
                CWP_COMM_PAR_WITHOUT_PART
                ! CWP_COMM_INTERNAL
    end enum

    ! CWP_Time_exch_t
    enum, bind(c)
        enumerator :: &
            CWP_TIME_EXCH_USER_CONTROLLED
            ! CWP_TIME_EXCH_EACH_TIME_STEP, &
            ! CWP_TIME_EXCH_N_TIME_STEP, &
            ! CWP_TIME_EXCH_CPL_TIME_STEP, &
            ! CWP_TIME_EXCH_ASYNCHRONOUS, &
            ! CWP_TIME_EXCH_SLAVE, &
            ! CWP_TIME_EXCH_MASTER
    end enum

    ! CWP_Dof_location_t
    enum, bind(c)
        enumerator :: &
                CWP_DOF_LOCATION_UNDEF, &
                CWP_DOF_LOCATION_CELL_CENTER, &
                CWP_DOF_LOCATION_NODE, &
                CWP_DOF_LOCATION_USER
    end enum

    ! CWP_Field_exch_t
    enum, bind(c)
        enumerator :: &
                CWP_FIELD_EXCH_SEND, &
                CWP_FIELD_EXCH_RECV, &
                CWP_FIELD_EXCH_SENDRECV
    end enum

    ! CWP_Field_exch_t
    enum, bind(c)
      enumerator :: &
                CWP_FIELD_MAP_SOURCE, &
                CWP_FIELD_MAP_TARGET
    end enum

    ! CWP_Field_storage_t
    enum, bind(c)
        enumerator :: &
                CWP_FIELD_STORAGE_INTERLACED, &
                CWP_FIELD_STORAGE_INTERLEAVED
    end enum

    ! CWP_Block_t
    enum, bind(c)
        enumerator :: &
                CWP_BLOCK_NODE, &
                CWP_BLOCK_EDGE2, &
                CWP_BLOCK_FACE_TRIA3, &
                CWP_BLOCK_FACE_QUAD4, &
                CWP_BLOCK_FACE_POLY, &
                CWP_BLOCK_CELL_TETRA4, &
                CWP_BLOCK_CELL_HEXA8, &
                CWP_BLOCK_CELL_PRISM6, &
                CWP_BLOCK_CELL_PYRAM5, &
                CWP_BLOCK_CELL_POLY
    end enum

    ! CWP_Spatial_interp_t
    enum, bind(c)
        enumerator :: &
                CWP_SPATIAL_INTERP_FROM_NEAREST_SOURCES_LEAST_SQUARES, &
                CWP_SPATIAL_INTERP_FROM_NEAREST_TARGETS_LEAST_SQUARES, &
                CWP_SPATIAL_INTERP_FROM_INTERSECTION, &
                CWP_SPATIAL_INTERP_FROM_LOCATION_MESH_LOCATION_LOCATE_ALL_TGT, &
                CWP_SPATIAL_INTERP_FROM_LOCATION_MESH_LOCATION_OCTREE, &
                CWP_SPATIAL_INTERP_FROM_LOCATION_MESH_LOCATION_BOXTREE, &
                CWP_SPATIAL_INTERP_FROM_IDENTITY
    end enum

    ! CWP_Interface_t
    enum, bind(c)
        enumerator :: &
                CWP_INTERFACE_POINT, &
                CWP_INTERFACE_LINEAR, &
                CWP_INTERFACE_SURFACE, &
                CWP_INTERFACE_VOLUME
    end enum

    ! CWP_Dynamic_mesh_t
    enum, bind(c)
        enumerator :: &
                CWP_DYNAMIC_MESH_STATIC, &
                CWP_DYNAMIC_MESH_DEFORMABLE, &
                CWP_DYNAMIC_MESH_VARIABLE
    end enum

    ! CWP_Status_t
    enum, bind(c)
        enumerator :: &
                CWP_STATUS_OFF, &
                CWP_STATUS_ON
    end enum

    ! CWP_Err_t
    enum, bind(c)
      enumerator :: &
        CWP_ERR_NO_ERROR, &
        CWP_ERR_DEFAULT
    end enum

    ! CWP_State_t
    enum, bind(c)
      enumerator :: &
        CWP_STATE_IN_PROGRESS, &
        CWP_STATE_END,         &
        CWP_STATE_OUTPUT_ERROR
    end enum

    ! CWP_Op_t
    enum, bind(c)
      enumerator :: &
        CWP_OP_MIN, &
        CWP_OP_MAX, &
        CWP_OP_SUM
    end enum

    ! CWP_PartData_exch_t
    enum, bind(c)
      enumerator :: &
        CWP_PARTDATA_SEND, &
        CWP_PARTDATA_RECV
    end enum

    interface CWP_Param_set; module procedure &
      CWP_Param_set_int_, &
      CWP_Param_set_double_, &
      CWP_Param_set_char_
    end interface CWP_Param_set

    interface CWP_Param_add; module procedure &
      CWP_Param_add_int_, &
      CWP_Param_add_double_, &
      CWP_Param_add_char_
    end interface CWP_Param_add

    interface CWP_Init ; module procedure &
        CWP_Init_
    end interface CWP_Init

    interface CWP_C_to_f_string ; module procedure &
        CWP_C_to_f_string_
    end interface CWP_C_to_f_string

    interface CWP_State_update ; module procedure &
        CWP_State_update_
    end interface CWP_State_update

    interface CWP_Time_step_beg ; module procedure &
        CWP_Time_step_beg_
    end interface CWP_Time_step_beg

    interface CWP_Time_step_end ; module procedure &
        CWP_Time_step_end_
    end interface CWP_Time_step_end

    interface CWP_User_structure_set ; module procedure &
        CWP_User_structure_set_
    end interface CWP_User_structure_set

    interface CWP_User_structure_get ; module procedure &
        CWP_User_structure_get_
    end interface CWP_User_structure_get

    interface CWP_Output_file_set ; module procedure &
        CWP_Output_file_set_
    end interface CWP_Output_file_set

    interface CWP_State_get ; module procedure &
        CWP_State_get_
    end interface CWP_State_get

    interface CWP_Cpl_create ; module procedure &
        CWP_Cpl_create_
    end interface CWP_Cpl_create

    interface CWP_Cpl_barrier ; module procedure &
        CWP_Cpl_barrier_
    end interface CWP_Cpl_barrier

    interface CWP_Cpl_Del ; module procedure &
        CWP_Cpl_Del_
    end interface CWP_Cpl_Del

    interface CWP_Computed_tgts_bcast_enable ; module procedure &
        CWP_Computed_tgts_bcast_enable_
    end interface CWP_Computed_tgts_bcast_enable

    interface CWP_Involved_srcs_bcast_enable ; module procedure &
        CWP_Involved_srcs_bcast_enable_
    end interface CWP_Involved_srcs_bcast_enable

    interface CWP_N_uncomputed_tgts_get ; module procedure &
        CWP_N_uncomputed_tgts_get_
    end interface CWP_N_uncomputed_tgts_get

    interface CWP_Uncomputed_tgts_get ; module procedure &
        CWP_Uncomputed_tgts_get_
    end interface CWP_Uncomputed_tgts_get

    interface CWP_N_computed_tgts_get ; module procedure &
        CWP_N_computed_tgts_get_
    end interface CWP_N_computed_tgts_get

    interface CWP_Computed_tgts_get ; module procedure &
        CWP_Computed_tgts_get_
    end interface CWP_Computed_tgts_get

    interface CWP_N_involved_srcs_get ; module procedure &
        CWP_N_involved_srcs_get_
    end interface CWP_N_involved_srcs_get

    interface CWP_Involved_srcs_get ; module procedure &
        CWP_Involved_srcs_get_
    end interface CWP_Involved_srcs_get

    interface CWP_Spatial_interp_weights_compute ; module procedure &
        CWP_Spatial_interp_weights_compute_
    end interface CWP_Spatial_interp_weights_compute

    interface CWP_Spatial_interp_property_set ; module procedure &
        CWP_Spatial_interp_property_set_
    end interface CWP_Spatial_interp_property_set

    interface CWP_Visu_set ; module procedure &
        CWP_Visu_set_
    end interface CWP_Visu_set

    interface CWP_User_tgt_pts_set ; module procedure &
        CWP_User_tgt_pts_set_
    end interface CWP_User_tgt_pts_set

    interface CWP_Mesh_interf_finalize ; module procedure &
        CWP_Mesh_interf_finalize_
    end interface CWP_Mesh_interf_finalize

    interface CWP_Mesh_interf_vtx_set ; module procedure &
        CWP_Mesh_interf_vtx_set_
    end interface CWP_Mesh_interf_vtx_set

    interface CWP_Mesh_interf_block_add ; module procedure &
        CWP_Mesh_interf_block_add_
    end interface CWP_Mesh_interf_block_add

    interface CWP_Mesh_interf_block_std_set ; module procedure &
        CWP_Mesh_interf_block_std_set_
    end interface CWP_Mesh_interf_block_std_set

    interface CWP_Mesh_interf_block_std_get ; module procedure &
        CWP_Mesh_interf_block_std_get_
    end interface CWP_Mesh_interf_block_std_get

    interface CWP_Mesh_interf_f_poly_block_set ; module procedure &
        CWP_Mesh_interf_f_poly_block_set_
    end interface CWP_Mesh_interf_f_poly_block_set

    interface CWP_Mesh_interf_f_poly_block_get ; module procedure &
        CWP_Mesh_interf_f_poly_block_get_
    end interface CWP_Mesh_interf_f_poly_block_get

    interface CWP_Mesh_interf_c_poly_block_set ; module procedure &
        CWP_Mesh_interf_c_poly_block_set_
    end interface CWP_Mesh_interf_c_poly_block_set

    interface CWP_Mesh_interf_c_poly_block_get ; module procedure &
        CWP_Mesh_interf_c_poly_block_get_
    end interface CWP_Mesh_interf_c_poly_block_get

    interface CWP_Mesh_interf_del ; module procedure &
        CWP_Mesh_interf_del_
    end interface CWP_Mesh_interf_del

    interface CWP_Mesh_interf_from_cellface_set ; module procedure &
        CWP_Mesh_interf_from_cellface_set_
    end interface CWP_Mesh_interf_from_cellface_set

    interface CWP_Mesh_interf_from_faceedge_set ; module procedure &
        CWP_Mesh_interf_from_faceedge_set_
    end interface CWP_Mesh_interf_from_faceedge_set

    interface CWP_Mesh_interf_from_facevtx_set ; module procedure &
        CWP_Mesh_interf_from_facevtx_set_
    end interface CWP_Mesh_interf_from_facevtx_set

    interface CWP_Field_create ; module procedure &
        CWP_Field_create_
    end interface CWP_Field_create

    interface CWP_Field_data_set ; module procedure &
        CWP_Field_data_set_
    end interface CWP_Field_data_set

    interface CWP_Field_dof_location_get ; module procedure &
        CWP_Field_dof_location_get_
    end interface CWP_Field_dof_location_get

    interface CWP_Field_storage_get ; module procedure &
        CWP_Field_storage_get_
    end interface CWP_Field_storage_get

    interface CWP_Field_n_dof_get ; module procedure &
        CWP_Field_n_dof_get_
    end interface CWP_Field_n_dof_get

    interface CWP_Field_del ; module procedure &
        CWP_Field_del_
    end interface CWP_Field_del

    interface CWP_Field_issend ; module procedure &
        CWP_Field_issend_
    end interface CWP_Field_issend

    interface CWP_Field_irecv ; module procedure &
        CWP_Field_irecv_
    end interface CWP_Field_irecv

    interface CWP_Field_wait_issend ; module procedure &
        CWP_Field_wait_issend_
    end interface CWP_Field_wait_issend

    interface CWP_Field_wait_irecv ; module procedure &
        CWP_Field_wait_irecv_
    end interface CWP_Field_wait_irecv

    interface CWP_Field_interp_function_unset ; module procedure &
        CWP_Field_interp_function_unset_
    end interface CWP_Field_interp_function_unset

    interface CWP_Field_interp_function_set ; module procedure &
        CWP_Field_interp_function_set_
    end interface CWP_Field_interp_function_set

    interface CWP_Field_n_components_get ; module procedure &
        CWP_Field_n_components_get_
    end interface CWP_Field_n_components_get

    interface CWP_Field_src_data_properties_get ; module procedure &
        CWP_Field_src_data_properties_get_
    end interface CWP_Field_src_data_properties_get

    interface CWP_Field_tgt_data_properties_get ; module procedure &
        CWP_Field_tgt_data_properties_get_
    end interface CWP_Field_tgt_data_properties_get

    interface CWP_Field_location_weights_get ; module procedure &
        CWP_Field_location_weights_get_
    end interface CWP_Field_location_weights_get

    interface CWP_Field_location_point_data_get ; module procedure &
        CWP_Field_location_point_data_get_
    end interface CWP_Field_location_point_data_get

    interface CWP_Field_location_internal_cell_vtx_get ; module procedure &
        CWP_Field_location_internal_cell_vtx_get_
    end interface CWP_Field_location_internal_cell_vtx_get

    interface CWP_Field_intersection_volumes_get ; module procedure &
        CWP_Field_intersection_volumes_get_
    end interface CWP_Field_intersection_volumes_get

    interface CWP_Field_intersection_tgt_elt_volumes_get ; module procedure &
        CWP_Field_intersection_tgt_elt_volumes_get_
    end interface CWP_Field_intersection_tgt_elt_volumes_get

    interface CWP_Field_nearest_neighbors_distances_get ; module procedure &
        CWP_Field_nearest_neighbors_distances_get_
    end interface CWP_Field_nearest_neighbors_distances_get

    interface CWP_Field_nearest_neighbors_coord_get ; module procedure &
        CWP_Field_nearest_neighbors_coord_get_
    end interface CWP_Field_nearest_neighbors_coord_get

    interface CWP_Param_del ; module procedure &
        CWP_Param_del_
    end interface CWP_Param_del

    interface CWP_Param_n_get ; module procedure &
        CWP_Param_n_get_
    end interface CWP_Param_n_get

    interface CWP_Param_is ; module procedure &
        CWP_Param_is_
    end interface CWP_Param_is

    interface CWP_Param_get
      module procedure CWP_Param_get_int
      module procedure CWP_Param_get_double
      module procedure CWP_Param_get_char
    end interface CWP_Param_get

    interface CWP_Param_reduce
      module procedure CWP_Param_reduce_int
      module procedure CWP_Param_reduce_double
      module procedure CWP_Param_reduce_char
    end interface CWP_Param_reduce

    interface CWP_Param_lock ; module procedure &
        CWP_Param_lock_
    end interface CWP_Param_lock

    interface CWP_Param_unlock ; module procedure &
        CWP_Param_unlock_
    end interface CWP_Param_unlock

    interface CWP_Codes_list_get ; module procedure &
        CWP_Codes_list_get_
    end interface CWP_Codes_list_get

    interface CWP_Loc_codes_list_get ; module procedure &
        CWP_Loc_codes_list_get_
    end interface CWP_Loc_codes_list_get

    interface CWP_Param_list_get ; module procedure &
        CWP_Param_list_get_
    end interface CWP_Param_list_get

    interface CWP_Global_data_issend
      module procedure CWP_Global_data_issend_int
      module procedure CWP_Global_data_issend_long
      module procedure CWP_Global_data_issend_double
      module procedure CWP_Global_data_issend_complex4
      module procedure CWP_Global_data_issend_complex8
      module procedure CWP_Global_data_issend_real4
    end interface CWP_Global_data_issend

    interface CWP_Global_data_irecv
      module procedure CWP_Global_data_irecv_int
      module procedure CWP_Global_data_irecv_long
      module procedure CWP_Global_data_irecv_double
      module procedure CWP_Global_data_irecv_complex4
      module procedure CWP_Global_data_irecv_complex8
      module procedure CWP_Global_data_irecv_real4
    end interface CWP_Global_data_irecv

    interface CWP_Global_data_wait_issend
      module procedure CWP_Global_data_wait_issend_
    end interface CWP_Global_data_wait_issend

    interface CWP_Global_data_wait_irecv
      module procedure CWP_Global_data_wait_irecv_
    end interface CWP_Global_data_wait_irecv

    interface CWP_Part_data_create
      module procedure CWP_Part_data_create_
    end interface CWP_Part_data_create

    interface CWP_Part_data_del
      module procedure CWP_Part_data_del_
    end interface CWP_Part_data_del

    interface CWP_Part_data_issend
      module procedure CWP_Part_data_issend_
    end interface CWP_Part_data_issend

    interface CWP_Part_data_irecv
      module procedure CWP_Part_data_irecv_
    end interface CWP_Part_data_irecv

    interface CWP_Part_data_wait_issend
      module procedure CWP_Part_data_wait_issend_
    end interface CWP_Part_data_wait_issend

    interface CWP_Part_data_wait_irecv
      module procedure CWP_Part_data_wait_irecv_
    end interface CWP_Part_data_wait_irecv

    interface CWP_Cpl_spatial_interp_algo_get
      module procedure CWP_Cpl_spatial_interp_algo_get_
    end interface CWP_Cpl_spatial_interp_algo_get

  !
  ! Private

  private :: c_f_char_array,&
             CWP_Init_ ,&
             CWP_C_to_f_string_ ,&
             CWP_State_update_ ,&
             CWP_Time_step_beg_ ,&
             CWP_Time_step_end_ ,&
             CWP_User_structure_set_ ,&
             CWP_User_structure_get_ ,&
             CWP_Output_file_set_,&
             CWP_State_get_ ,&
             CWP_Cpl_create_ ,&
             CWP_Cpl_barrier_ ,&
             CWP_Cpl_Del_ ,&
             CWP_Computed_tgts_bcast_enable_,&
             CWP_Involved_srcs_bcast_enable_,&
             CWP_N_uncomputed_tgts_get_ ,&
             CWP_Uncomputed_tgts_get_ ,&
             CWP_N_computed_tgts_get_ ,&
             CWP_Computed_tgts_get_ ,&
             CWP_N_involved_srcs_get_ ,&
             CWP_Involved_srcs_get_ ,&
             CWP_Spatial_interp_weights_compute_ ,&
             CWP_Spatial_interp_property_set_ ,&
             CWP_Visu_set_ ,&
             CWP_User_tgt_pts_set_ ,&
             CWP_Mesh_interf_finalize_ ,&
             CWP_Mesh_interf_vtx_set_ ,&
             CWP_Mesh_interf_block_add_ ,&
             CWP_Mesh_interf_block_std_set_ ,&
             CWP_Mesh_interf_block_std_get_ ,&
             CWP_Mesh_interf_f_poly_block_set_ ,&
             CWP_Mesh_interf_f_poly_block_get_ ,&
             CWP_Mesh_interf_c_poly_block_set_ ,&
             CWP_Mesh_interf_c_poly_block_get_ ,&
             CWP_Mesh_interf_del_ ,&
             CWP_Mesh_interf_from_cellface_set_ ,&
             CWP_Mesh_interf_from_faceedge_set_ ,&
             CWP_Mesh_interf_from_facevtx_set_ ,&
             CWP_Field_create_ ,&
             CWP_Field_data_set_ ,&
             CWP_Field_dof_location_get_ ,&
             CWP_Field_storage_get_ ,&
             CWP_Field_n_dof_get_, &
             CWP_Field_del_ ,&
             CWP_Field_issend_ ,&
             CWP_Field_irecv_ ,&
             CWP_Field_wait_issend_ ,&
             CWP_Field_wait_irecv_ ,&
             CWP_Field_interp_function_unset_ ,&
             CWP_Field_interp_function_set_ ,&
             CWP_Field_n_components_get_, &
             CWP_Field_src_data_properties_get_, &
             CWP_Field_tgt_data_properties_get_, &
             CWP_Field_location_weights_get_, &
             CWP_Field_location_point_data_get_, &
             CWP_Field_location_internal_cell_vtx_get_, &
             CWP_Field_intersection_volumes_get_, &
             CWP_Field_intersection_tgt_elt_volumes_get_, &
             CWP_Field_nearest_neighbors_distances_get_, &
             CWP_Field_nearest_neighbors_coord_get_, &
             CWP_Param_add_int_ ,&
             CWP_Param_add_double_ ,&
             CWP_Param_add_char_ ,&
             CWP_Param_set_int_ ,&
             CWP_Param_set_double_ ,&
             CWP_Param_set_char_ ,&
             CWP_Param_del_ ,&
             CWP_Param_n_get_ ,&
             CWP_Param_is_ ,&
             CWP_Param_get_int ,&
             CWP_Param_get_double ,&
             CWP_Param_get_char, &
             CWP_Param_reduce_int,&
             CWP_Param_reduce_double,&
             CWP_Param_reduce_char, &
             CWP_Param_lock_ ,&
             CWP_Param_unlock_,&
             CWP_Codes_list_get_,&
             CWP_Loc_codes_list_get_,&
             CWP_Param_list_get_, &
             CWP_Global_data_issend_int, &
             CWP_Global_data_issend_long, &
             CWP_Global_data_issend_double, &
             CWP_Global_data_issend_complex4, &
             CWP_Global_data_issend_complex8, &
             CWP_Global_data_issend_real4, &
             CWP_Global_data_irecv_int, &
             CWP_Global_data_irecv_long, &
             CWP_Global_data_irecv_double, &
             CWP_Global_data_irecv_complex4, &
             CWP_Global_data_irecv_complex8, &
             CWP_Global_data_irecv_real4, &
             CWP_Global_data_wait_issend_, &
             CWP_Global_data_wait_irecv_, &
             CWP_Part_data_create_, &
             CWP_Part_data_del_, &
             CWP_Part_data_issend_, &
             CWP_Part_data_irecv_, &
             CWP_Part_data_wait_issend_, &
             CWP_Part_data_wait_irecv_, &
             CWP_Cpl_spatial_interp_algo_get_

    interface

    !> \cond DOXYGEN_SHOULD_SKIP_THIS
      subroutine CWP_Init_cf(fcomm, n_code, code_names, l_code_names, is_active_rank, intra_comms) &
              bind(c, name = 'CWP_Init_cf')
        use, intrinsic :: iso_c_binding
        implicit none
        integer(c_int), value :: fcomm
        integer(c_int), value :: n_code
        type(c_ptr),    value :: code_names
        type(c_ptr),    value :: l_code_names
        integer(c_int), value :: is_active_rank
        type(c_ptr),    value :: intra_comms
      end subroutine CWP_Init_cf

      subroutine CWP_State_update_cf(local_code_name, l_local_code_name, state) &
        bind(c, name='CWP_State_update_cf')
        use, intrinsic :: iso_c_binding
        implicit none
        character(kind = c_char, len = 1) :: local_code_name
        integer(c_int), value             :: l_local_code_name, state
      end subroutine CWP_State_update_cf

      subroutine CWP_Time_step_beg_cf(local_code_name, l_local_code_name, current_time) &
        bind(c, name='CWP_Time_step_beg_cf')
        use, intrinsic :: iso_c_binding
        implicit none
        character(kind = c_char, len = 1) :: local_code_name
        integer(c_int), value             :: l_local_code_name
        real(c_double), value             :: current_time
      end subroutine CWP_Time_step_beg_cf

      subroutine CWP_Time_step_end_cf(local_code_name, l_local_code_name) &
        bind(c, name='CWP_Time_step_end_cf')
        use, intrinsic :: iso_c_binding
        implicit none
        character(kind = c_char, len = 1) :: local_code_name
        integer(c_int), value             :: l_local_code_name
      end subroutine CWP_Time_step_end_cf

      subroutine CWP_User_structure_set_cf(local_code_name, l_local_code_name, user_structure) &
        bind (c, name="CWP_User_structure_set_cf")
        use, intrinsic :: iso_c_binding
        implicit none
        character(kind = c_char, len = 1) :: local_code_name
        integer(c_int), value             :: l_local_code_name
        type(c_ptr),    value             :: user_structure
      end subroutine CWP_User_structure_set_cf

      function CWP_User_structure_get_cf(local_code_name, l_local_code_name) &
        result (user_structure)                                              &
        bind (c, name="CWP_User_structure_get_cf")
        use, intrinsic :: iso_c_binding
        implicit none
        character(kind = c_char, len = 1) :: local_code_name
        integer(c_int), value             :: l_local_code_name
        type(c_ptr)                       :: user_structure
      end function CWP_User_structure_get_cf

      subroutine CWP_Output_file_set_cf(f_output_file_name, l_output_file_name) &
        bind (c, name="CWP_Output_file_set_cf")
        use, intrinsic :: iso_c_binding
        implicit none
        character(kind = c_char, len = 1) :: f_output_file_name
        integer(c_int), value             :: l_output_file_name
      end subroutine CWP_Output_file_set_cf

      function CWP_State_get_cf(local_code_name, l_local_code_name) &
        result (state)                                              &
        bind(c, name='CWP_State_get_cf')
        use, intrinsic :: iso_c_binding
        implicit none
        character(kind = c_char, len = 1) :: local_code_name
        integer(c_int), value             :: l_local_code_name
        integer(c_int)                    :: state
      end function CWP_State_get_cf

      subroutine CWP_Cpl_create_cf(local_code_name, l_local_code_name, cpl_id, l_cpl_id, coupled_code_name, l_coupled_code_name, &
        entities_dim, comm_type, spatial_interp, n_part, displacement, freq) &
        bind(c, name = 'CWP_Cpl_create_cf')
        use, intrinsic :: iso_c_binding
        implicit none
        character(kind = c_char, len = 1) :: local_code_name, cpl_id, coupled_code_name
        integer(c_int), value :: l_local_code_name, l_cpl_id, l_coupled_code_name
        integer(c_int), value :: entities_dim
        integer(c_int), value :: comm_type
        integer(c_int), value :: spatial_interp
        integer(c_int), value :: n_part
        integer(c_int), value :: displacement
        integer(c_int), value :: freq
      end subroutine CWP_Cpl_Create_cf

      subroutine CWP_Cpl_barrier_cf(local_code_name, l_local_code_name, cpl_id, l_cpl_id) &
        bind(c, name = 'CWP_Cpl_barrier_cf')
        use, intrinsic :: iso_c_binding
        implicit none
        character(kind = c_char, len = 1) :: local_code_name, cpl_id
        integer(c_int), value             :: l_local_code_name, l_cpl_id
      end subroutine CWP_Cpl_barrier_cf

      subroutine CWP_Cpl_del_cf(local_code_name, l_local_code_name, cpl_id, l_cpl_id) &
        bind(c, name = 'CWP_Cpl_del_cf')
        use, intrinsic :: iso_c_binding
        implicit none
        character(kind = c_char, len = 1) :: local_code_name, cpl_id
        integer(kind = c_int), value :: l_local_code_name, l_cpl_id
      end subroutine CWP_Cpl_del_cf

      subroutine CWP_Computed_tgts_bcast_enable_cf(local_code_name, l_local_code_name, cpl_id, l_cpl_id, field_id, l_field_id) &
        bind(c, name = 'CWP_Computed_tgts_bcast_enable_cf')
        use, intrinsic :: iso_c_binding
        implicit none
        character(kind = c_char, len = 1) :: local_code_name, cpl_id, field_id
        integer(kind = c_int), value :: l_local_code_name, l_cpl_id, l_field_id
      end subroutine CWP_Computed_tgts_bcast_enable_cf

      subroutine CWP_Involved_srcs_bcast_enable_cf(local_code_name, l_local_code_name, cpl_id, l_cpl_id, field_id, l_field_id) &
        bind(c, name = 'CWP_Involved_srcs_bcast_enable_cf')
        use, intrinsic :: iso_c_binding
        implicit none
        character(kind = c_char, len = 1) :: local_code_name, cpl_id, field_id
        integer(kind = c_int), value :: l_local_code_name, l_cpl_id, l_field_id
      end subroutine CWP_Involved_srcs_bcast_enable_cf

      function CWP_N_uncomputed_tgts_get_cf (local_code_name,   &
                                             l_local_code_name, &
                                             cpl_id,            &
                                             l_cpl_id,          &
                                             field_id,          &
                                             l_field_id,        &
                                             i_part)            &
                                             result (n_uncomputed_tgts) &
        bind(c, name = 'CWP_N_uncomputed_tgts_get_cf')
        use, intrinsic :: iso_c_binding
        implicit none

        character(kind = c_char, len = 1) :: local_code_name, cpl_id
        integer(kind = c_int), value :: l_local_code_name, l_cpl_id
        character(kind = c_char, len = 1) :: field_id
        integer(kind = c_int), value :: l_field_id
        integer(c_int), value :: i_part
        integer(kind = c_int) :: n_uncomputed_tgts
      end function CWP_N_uncomputed_tgts_get_cf

      function CWP_N_computed_tgts_get_cf(local_code_name,   &
                                          l_local_code_name, &
                                          cpl_id,            &
                                          l_cpl_id,          &
                                          field_id,          &
                                          l_field_id,        &
                                          i_part)            &
                                          result (n_computed_tgts) &
            bind(c, name = 'CWP_N_computed_tgts_get_cf')
        use, intrinsic :: iso_c_binding
        implicit none
        character(kind = c_char, len = 1) :: local_code_name, cpl_id
        integer(kind = c_int), value :: l_local_code_name, l_cpl_id
        character(kind = c_char, len = 1) :: field_id
        integer(kind = c_int), value :: l_field_id
        integer(c_int), value :: i_part
        integer(kind = c_int) :: n_computed_tgts
      end function CWP_N_computed_tgts_get_cf

      function CWP_N_involved_srcs_get_cf(local_code_name,   &
              l_local_code_name, &
              cpl_id,            &
              l_cpl_id,          &
              field_id,          &
              l_field_id,        &
              i_part)            &
              result (n_involved_srcs) &
                      bind(c, name = 'CWP_N_involved_srcs_get_cf')
        use, intrinsic :: iso_c_binding
        implicit none
        character(kind = c_char, len = 1) :: local_code_name, cpl_id
        integer(kind = c_int), value :: l_local_code_name, l_cpl_id
        character(kind = c_char, len = 1) :: field_id
        integer(kind = c_int), value :: l_field_id
        integer(c_int), value :: i_part
        integer(kind = c_int) :: n_involved_srcs
      end function CWP_N_involved_srcs_get_cf

      function CWP_Computed_tgts_get_cf(local_code_name,   &
                                        l_local_code_name, &
                                        cpl_id,            &
                                        l_cpl_id,          &
                                        field_id,          &
                                        l_field_id,        &
                                        i_part)            &
                                        result (computed_tgts) &
            bind(c, name = 'CWP_Computed_tgts_get_cf')
        use, intrinsic :: iso_c_binding
        implicit none
        character(kind = c_char, len = 1) :: local_code_name, cpl_id
        integer(kind = c_int), value :: l_local_code_name, l_cpl_id
        character(kind = c_char, len = 1) :: field_id
        integer(kind = c_int), value :: l_field_id
        integer(c_int), value :: i_part
        type(c_ptr) :: computed_tgts
      end function CWP_Computed_tgts_get_cf


      function CWP_Uncomputed_tgts_get_cf (local_code_name,   &
                                          l_local_code_name, &
                                          cpl_id,            &
                                          l_cpl_id,          &
                                          field_id,          &
                                          l_field_id,        &
                                          i_part)            &
                                          result (uncomputed_tgts) &
        bind(c, name = 'CWP_Uncomputed_tgts_get_cf')
        use, intrinsic :: iso_c_binding
        implicit none
        character(kind = c_char, len = 1) :: local_code_name, cpl_id
        integer(kind = c_int), value :: l_local_code_name, l_cpl_id
        character(kind = c_char, len = 1) :: field_id
        integer(kind = c_int), value :: l_field_id
        integer(c_int), value :: i_part
        type(c_ptr) :: uncomputed_tgts
      end function CWP_Uncomputed_tgts_get_cf

      function CWP_Involved_srcs_get_cf(local_code_name,   &
              l_local_code_name, &
              cpl_id,            &
              l_cpl_id,          &
              field_id,          &
              l_field_id,        &
              i_part)            &
              result (involved_srcs) &
                      bind(c, name = 'CWP_Involved_srcs_get_cf')
        use, intrinsic :: iso_c_binding
        implicit none
        character(kind = c_char, len = 1) :: local_code_name, cpl_id
        integer(kind = c_int), value :: l_local_code_name, l_cpl_id
        character(kind = c_char, len = 1) :: field_id
        integer(kind = c_int), value :: l_field_id
        integer(c_int), value :: i_part
        type(c_ptr) :: involved_srcs
      end function CWP_Involved_srcs_get_cf

      function CWP_Computed_tgts_dist_to_spatial_interp_get_cf(local_code_name, l_local_code_name, cpl_id, l_cpl_id) &
            result (dists) &
        bind(c, name = 'CWP_Computed_tgts_dist_to_spatial_interp_get_cf')
        use, intrinsic :: iso_c_binding
        implicit none
        character(kind = c_char, len = 1) :: local_code_name, cpl_id
        integer(kind = c_int), value :: l_local_code_name, l_cpl_id
        type(c_ptr) :: dists
      end function CWP_Computed_tgts_dist_to_spatial_interp_get_cf

      subroutine CWP_Spatial_interp_weights_compute_cf(local_code_name, l_local_code_name, cpl_id, l_cpl_id) &
            bind(c, name = 'CWP_Spatial_interp_weights_compute_cf')
        use, intrinsic :: iso_c_binding
        implicit none
        character(kind = c_char, len = 1) :: local_code_name, cpl_id
        integer(kind = c_int), value :: l_local_code_name, l_cpl_id
      end subroutine CWP_Spatial_interp_weights_compute_cf

      subroutine CWP_Visu_set_cf(local_code_name, l_local_code_name, cpl_id, l_cpl_id, freq, &
            format, format_option, l_format_option) &
            bind(c, name = 'CWP_Visu_set_cf')
        use, intrinsic :: iso_c_binding
        implicit none
        character(kind = c_char, len = 1) :: local_code_name, cpl_id, format_option
        integer(c_int), value :: freq
        integer(c_int), value :: format
        integer(kind = c_int), value :: l_local_code_name, l_cpl_id, l_format_option
      end subroutine CWP_Visu_set_cf

      subroutine CWP_User_tgt_pts_set_cf(local_code_name, l_local_code_name, cpl_id, l_cpl_id, i_part, n_pts, coord, global_num) &
              bind(c, name = 'CWP_User_tgt_pts_set_cf')
        use, intrinsic :: iso_c_binding
        implicit none

        character(kind = c_char, len = 1) :: local_code_name, cpl_id
        integer(kind = c_int), value :: l_local_code_name, l_cpl_id
        integer(c_int), value :: i_part, n_pts
        type(c_ptr), value :: coord
        type(c_ptr), value :: global_num
      end subroutine CWP_User_tgt_pts_set_cf

      subroutine CWP_Mesh_interf_finalize_cf(local_code_name, l_local_code_name, cpl_id, l_cpl_id) &
            bind(c, name = 'CWP_Mesh_interf_finalize_cf')
        use, intrinsic :: iso_c_binding
        implicit none
        character(kind = c_char, len = 1) :: local_code_name, cpl_id
        integer(kind = c_int), value :: l_local_code_name, l_cpl_id
      end subroutine CWP_Mesh_interf_finalize_cf

      subroutine CWP_Mesh_interf_vtx_set_cf(local_code_name, l_local_code_name, cpl_id, l_cpl_id, &
            i_part, n_pts, coord, global_num) &
            bind(c, name = 'CWP_Mesh_interf_vtx_set_cf')
        use, intrinsic :: iso_c_binding
        implicit none
        character(kind = c_char, len = 1) :: local_code_name, cpl_id
        integer(c_int), value :: i_part, n_pts
        type(c_ptr), value :: coord
        type(c_ptr), value :: global_num
        integer(kind = c_int), value :: l_local_code_name, l_cpl_id
      end subroutine CWP_Mesh_interf_vtx_set_cf

      function CWP_Mesh_interf_block_add_cf(local_code_name, l_local_code_name, cpl_id, l_cpl_id, block_type) result(block_id) &
        bind(c, name = 'CWP_Mesh_interf_block_add_cf')
        use, intrinsic :: iso_c_binding
        implicit none
        character(kind = c_char, len = 1) :: local_code_name, cpl_id
        integer(c_int), value :: block_type
        integer(c_int) :: block_id
        integer(c_int), value :: l_local_code_name, l_cpl_id
      end function CWP_Mesh_interf_block_add_cf

      subroutine CWP_Mesh_interf_block_std_set_cf(local_code_name, l_local_code_name, cpl_id, l_cpl_id, &
            i_part, block_id, n_elts, connec, global_num) &
        bind(c, name = 'CWP_Mesh_interf_block_std_set_cf')
        use, intrinsic :: iso_c_binding
        implicit none
        character(kind = c_char, len = 1) :: local_code_name, cpl_id
        integer(c_int), value :: i_part, block_id, n_elts
        type(c_ptr), value :: connec, global_num
        integer(kind = c_int), value :: l_local_code_name, l_cpl_id
      end subroutine CWP_Mesh_interf_block_std_set_cf

      subroutine CWP_Mesh_interf_f_poly_block_set_cf(local_code_name, l_local_code_name, cpl_id, l_cpl_id, &
            i_part, block_id, n_elts, connec_idx, connec, global_num) &
        bind(c, name = 'CWP_Mesh_interf_f_poly_block_set_cf')
        use, intrinsic :: iso_c_binding
        implicit none
        character(kind = c_char, len = 1) :: local_code_name, cpl_id
        integer(c_int), value :: i_part, block_id, n_elts
        type(c_ptr), value :: connec_idx, connec, global_num
        integer(kind = c_int), value :: l_local_code_name, l_cpl_id
      end subroutine CWP_Mesh_interf_f_poly_block_set_cf

      subroutine CWP_Mesh_interf_c_poly_block_set_cf(local_code_name, l_local_code_name, cpl_id, l_cpl_id, &
            i_part, block_id, n_elts, n_faces, connec_faces_idx, connec_faces, connec_cells_idx, connec_cells, global_num) &
        bind(c, name = 'CWP_Mesh_interf_c_poly_block_set_cf')
        use, intrinsic :: iso_c_binding
        implicit none
        character(kind = c_char, len = 1) :: local_code_name, cpl_id
        integer(c_int), value :: i_part, block_id, n_elts, n_faces
        type(c_ptr), value :: connec_faces_idx, connec_faces, connec_cells_idx, connec_cells, global_num
        integer(kind = c_int), value :: l_local_code_name, l_cpl_id
      end subroutine CWP_Mesh_interf_c_poly_block_set_cf

      subroutine CWP_Mesh_interf_del_cf(local_code_name, l_local_code_name, cpl_id, l_cpl_id) &
              bind(c, name = 'CWP_Mesh_interf_del_cf')
        use, intrinsic :: iso_c_binding
        implicit none
        character(kind = c_char, len = 1) :: local_code_name, cpl_id
        integer(kind = c_int), value :: l_local_code_name, l_cpl_id
      end subroutine CWP_Mesh_interf_del_cf

      subroutine CWP_Mesh_interf_from_cellface_set_cf(local_code_name, l_local_code_name, cpl_id, l_cpl_id, i_part, n_cells, &
            cell_face_idx, cell_face, n_faces, face_vtx_idx, face_vtx, parent_num) &
            bind(c, name = 'CWP_Mesh_interf_from_cellface_set_cf')
        use, intrinsic :: iso_c_binding
        implicit none
        character(kind = c_char, len = 1) :: local_code_name, cpl_id
        integer(c_int), value :: i_part, n_cells, n_faces
        type(c_ptr), value :: cell_face_idx, cell_face, face_vtx_idx, face_vtx, parent_num
        integer(kind = c_int), value :: l_local_code_name, l_cpl_id
      end subroutine CWP_Mesh_interf_from_cellface_set_cf

      subroutine CWP_Mesh_interf_from_faceedge_set_cf(local_code_name, l_local_code_name, cpl_id, l_cpl_id, i_part, n_faces, &
            face_edge_idx, face_edge, n_edges, edge_vtx, parent_num) &
            bind(c, name = 'CWP_Mesh_interf_from_faceedge_set_cf')
        use, intrinsic :: iso_c_binding
        implicit none
        character(kind = c_char, len = 1) :: local_code_name, cpl_id
        integer(c_int), value :: i_part, n_faces, n_edges
        type(c_ptr), value :: face_edge_idx, face_edge, edge_vtx, parent_num
        integer(kind = c_int), value :: l_local_code_name, l_cpl_id
      end subroutine CWP_Mesh_interf_from_faceedge_set_cf

      subroutine CWP_Mesh_interf_from_facevtx_set_cf(local_code_name, l_local_code_name, cpl_id, l_cpl_id, i_part, n_faces, &
            face_vtx_idx, face_vtx, global_num) &
            bind(c, name = 'CWP_Mesh_interf_from_facevtx_set_cf')
        use, intrinsic :: iso_c_binding
        implicit none
        character(kind = c_char, len = 1) :: local_code_name, cpl_id
        integer(c_int), value :: i_part, n_faces
        type(c_ptr), value :: face_vtx_idx, face_vtx, global_num
        integer(kind = c_int), value :: l_local_code_name, l_cpl_id
      end subroutine CWP_Mesh_interf_from_facevtx_set_cf

      subroutine CWP_Field_create_cf(local_code_name, l_local_code_name, cpl_id, l_cpl_id, field_id, l_field_id, &
          data_type, storage, n_component, target_location, exch_type, visu_status) &
          bind(c, name = 'CWP_Field_create_cf')
        use, intrinsic :: iso_c_binding
        implicit none
        character(kind = c_char, len = 1) :: local_code_name, cpl_id, field_id
        integer(c_int), value :: data_type, storage, n_component, target_location, exch_type, visu_status
        integer(kind = c_int), value :: l_local_code_name, l_cpl_id, l_field_id
      end subroutine

      subroutine CWP_Field_data_set_cf(local_code_name, l_local_code_name, cpl_id, l_cpl_id, field_id, l_field_id, &
            i_part, map_type, data) &
            bind(c, name = 'CWP_Field_data_set_cf')
        use, intrinsic :: iso_c_binding
        implicit none
        character(kind = c_char, len = 1) :: local_code_name, cpl_id, field_id
        integer(c_int), value :: i_part
        integer(c_int), value :: map_type
        type(c_ptr), value :: data
        integer(kind = c_int), value :: l_local_code_name, l_cpl_id, l_field_id
      end subroutine CWP_Field_data_set_cf

      subroutine CWP_Field_issend_cf(local_code_name, l_local_code_name, cpl_id, l_cpl_id, field_id, l_field_id) &
            bind(c, name = 'CWP_Field_issend_cf')
        use, intrinsic :: iso_c_binding
        implicit none
        character(kind = c_char, len = 1) :: local_code_name, cpl_id, field_id
        integer(kind = c_int), value :: l_local_code_name, l_cpl_id, l_field_id
      end subroutine CWP_Field_issend_cf

      subroutine CWP_Field_irecv_cf(local_code_name, l_local_code_name, cpl_id, l_cpl_id, tgt_field_id, l_tgt_field_id) &
            bind(c, name = 'CWP_Field_irecv_cf')
        use, intrinsic :: iso_c_binding
        implicit none
        character(kind = c_char, len = 1) :: local_code_name, cpl_id, tgt_field_id
        integer(kind = c_int), value :: l_local_code_name, l_cpl_id, l_tgt_field_id
      end subroutine CWP_Field_irecv_cf

      subroutine CWP_Field_wait_issend_cf(local_code_name, l_local_code_name, cpl_id, l_cpl_id, field_id, l_field_id) &
            bind(c, name = 'CWP_Field_wait_issend_cf')
        use, intrinsic :: iso_c_binding
        implicit none
        character(kind = c_char, len = 1) :: local_code_name, cpl_id, field_id
        integer(kind = c_int), value :: l_local_code_name, l_cpl_id, l_field_id
      end subroutine CWP_Field_wait_issend_cf

      subroutine CWP_Field_wait_irecv_cf(local_code_name, l_local_code_name, cpl_id, l_cpl_id, tgt_field_id, l_tgt_field_id) &
            bind(c, name = 'CWP_Field_wait_irecv_cf')
        use, intrinsic :: iso_c_binding
        implicit none
        character(kind = c_char, len = 1) :: local_code_name, cpl_id, tgt_field_id
        integer(kind = c_int), value :: l_local_code_name, l_cpl_id, l_tgt_field_id
      end subroutine CWP_Field_wait_irecv_cf

      subroutine CWP_Field_interp_function_set_cf(local_code_name, l_local_code_name, cpl_id, l_cpl_id, &
              field_id, l_field_id, user_interpolation_fct) &
              bind(c, name = 'CWP_Field_interp_function_set_cf')
        use, intrinsic :: iso_c_binding
        implicit none

        character(kind = c_char, len = 1)       :: local_code_name, cpl_id, field_id
        integer(kind = c_int),            value :: l_local_code_name, l_cpl_id, l_field_id
        type(c_funptr),                   value :: user_interpolation_fct
      end subroutine CWP_Field_interp_function_set_cf

      subroutine CWP_Field_interp_function_unset_cf(local_code_name, l_local_code_name, cpl_id, l_cpl_id, &
              field_id, l_field_id) &
              bind(c, name = 'CWP_Field_interp_function_unset_cf')
        use, intrinsic :: iso_c_binding
        implicit none

        character(kind = c_char, len = 1) :: local_code_name, cpl_id, field_id
        integer(kind = c_int), value :: l_local_code_name, l_cpl_id, l_field_id
      end subroutine CWP_Field_interp_function_unset_cf

      function CWP_Field_n_components_get_cf(local_code_name, l_local_code_name, cpl_id, l_cpl_id, &
                                             field_id, l_field_id) &
        result (n_components)                                                          &
        bind(c, name='CWP_Field_n_components_get_cf')
        use, intrinsic :: iso_c_binding
        implicit none
        character(kind = c_char, len = 1) :: local_code_name, cpl_id, field_id
        integer(kind = c_int), value      :: l_local_code_name, l_cpl_id, l_field_id
        integer(c_int)                    :: n_components
      end function CWP_Field_n_components_get_cf

      subroutine CWP_Field_src_data_properties_get_cf(local_code_name,   &
                                                      l_local_code_name, &
                                                      cpl_id,            &
                                                      l_cpl_id,          &
                                                      field_id,          &
                                                      l_field_id,        &
                                                      i_part,            &
                                                      n_elt_src,         &
                                                      c_src_to_tgt_idx)  &
          bind(c, name = 'CWP_Field_src_data_properties_get_cf')
          use, intrinsic :: iso_c_binding
          implicit none
          character(kind = c_char, len = 1) :: local_code_name, cpl_id, field_id
          integer(kind = c_int), value      :: l_local_code_name, l_cpl_id, l_field_id
          integer(c_int), value             :: i_part
          integer(c_int)                    :: n_elt_src
          type(c_ptr)                       :: c_src_to_tgt_idx
      end subroutine CWP_Field_src_data_properties_get_cf

      subroutine CWP_Field_tgt_data_properties_get_cf(local_code_name,         &
                                                      l_local_code_name,       &
                                                      cpl_id,                  &
                                                      l_cpl_id,                &
                                                      field_id,                &
                                                      l_field_id,              &
                                                      i_part,                  &
                                                      n_elt_tgt,               &
                                                      n_referenced_tgt,        &
                                                      c_referenced_tgt,        &
                                                      c_tgt_come_from_src_idx) &
          bind(c, name = 'CWP_Field_tgt_data_properties_get_cf')
          use, intrinsic :: iso_c_binding
          implicit none
          character(kind = c_char, len = 1) :: local_code_name, cpl_id, field_id
          integer(kind = c_int), value      :: l_local_code_name, l_cpl_id, l_field_id
          integer(c_int), value             :: i_part
          integer(c_int)                    :: n_elt_tgt, n_referenced_tgt
          type(c_ptr)                       :: c_referenced_tgt, c_tgt_come_from_src_idx
      end subroutine CWP_Field_tgt_data_properties_get_cf

      subroutine CWP_Field_location_weights_get_cf(local_code_name,          &
                                                   l_local_code_name,        &
                                                   cpl_id,                   &
                                                   l_cpl_id,                 &
                                                   field_id,                 &
                                                   l_field_id,               &
                                                   i_part,                   &
                                                   c_weights,                &
                                                   s_weights)                & !! TO DO: get with coupling
          bind(c, name = 'CWP_Field_location_weights_get_cf')
          use, intrinsic :: iso_c_binding
          implicit none
          character(kind = c_char, len = 1) :: local_code_name, cpl_id, field_id
          integer(kind = c_int), value      :: l_local_code_name, l_cpl_id, l_field_id
          integer(c_int), value             :: i_part
          type(c_ptr)                       :: c_weights
          integer(kind = c_int)             :: s_weights
      end subroutine CWP_Field_location_weights_get_cf

      subroutine CWP_Field_location_point_data_get_cf(local_code_name,           &
                                                       l_local_code_name,        &
                                                       cpl_id,                   &
                                                       l_cpl_id,                 &
                                                       field_id,                 &
                                                       l_field_id,               &
                                                       i_part,                   &
                                                       c_points_coords,          &
                                                       c_points_uvw,             &
                                                       c_points_dist2,           &
                                                       c_points_projected_coords,&
                                                       s_size)&

          bind(c, name = 'CWP_Field_location_point_data_get_cf')
          use, intrinsic :: iso_c_binding
          implicit none
          character(kind = c_char, len = 1) :: local_code_name, cpl_id, field_id
          integer(kind = c_int), value      :: l_local_code_name, l_cpl_id, l_field_id
          integer(c_int), value             :: i_part
          type(c_ptr)                       :: c_points_coords, c_points_uvw, c_points_dist2, c_points_projected_coords
          integer(kind = c_int) :: s_size
      end subroutine CWP_Field_location_point_data_get_cf

      subroutine CWP_Field_intersection_volumes_get_cf(local_code_name,      &
                                                        l_local_code_name,    &
                                                        cpl_id,               &
                                                        l_cpl_id,             &
                                                        field_id,         &
                                                        l_field_id,       &
                                                        i_part,               &
                                                        c_volumes,            &
                                                        s_volumes)            &
          bind(c, name = 'CWP_Field_intersection_volumes_get_cf')
          use, intrinsic :: iso_c_binding
          implicit none
          character(kind = c_char, len = 1) :: local_code_name, cpl_id, field_id
          integer(kind = c_int), value      :: l_local_code_name, l_cpl_id, l_field_id
          integer(c_int), value             :: i_part
          type(c_ptr)                       :: c_volumes
          integer(kind = c_int)             :: s_volumes
      end subroutine CWP_Field_intersection_volumes_get_cf

      subroutine CWP_Field_intersection_tgt_elt_volumes_get_cf(local_code_name,      &
                                                                l_local_code_name,    &
                                                                cpl_id,               &
                                                                l_cpl_id,             &
                                                                field_id,         &
                                                                l_field_id,       &
                                                                i_part,               &
                                                                c_tgt_elt_volumes,    &
                                                                n_elt)                &
          bind(c, name = 'CWP_Field_intersection_tgt_elt_volumes_get_cf')
          use, intrinsic :: iso_c_binding
          implicit none
          character(kind=c_char, len = 1) :: local_code_name, cpl_id, field_id
          integer(c_int), value           :: l_local_code_name, l_cpl_id, l_field_id
          integer(c_int), value           :: i_part
          type(c_ptr)                     :: c_tgt_elt_volumes
          integer(c_int)                  :: n_elt
      end subroutine CWP_Field_intersection_tgt_elt_volumes_get_cf

      subroutine CWP_Field_nearest_neighbors_distances_get_cf(local_code_name,   &
                                                               l_local_code_name, &
                                                               cpl_id,            &
                                                               l_cpl_id,          &
                                                               field_id,          &
                                                               l_field_id,        &
                                                               i_part,            &
                                                               c_distances2,      &
                                                               s_distances2)      &
          bind(c, name = 'CWP_Field_nearest_neighbors_distances_get_cf')
          use, intrinsic :: iso_c_binding
          implicit none
          character(kind = c_char, len = 1) :: local_code_name, cpl_id, field_id
          integer(kind = c_int), value      :: l_local_code_name, l_cpl_id, l_field_id
          integer(c_int), value             :: i_part
          type(c_ptr)                       :: c_distances2
          integer(kind = c_int)             :: s_distances2
      end subroutine CWP_Field_nearest_neighbors_distances_get_cf

      subroutine CWP_Field_nearest_neighbors_coord_get_cf(f_local_code_name,   &
                                                           l_local_code_name,   &
                                                           f_cpl_id,            &
                                                           l_cpl_id,            &
                                                           f_field_id,          &
                                                           l_field_id,          &
                                                           i_part,              &
                                                           c_nearest_src_coord, &
                                                           n_nearest_src_pts)   &
      bind (c, name='CWP_Field_nearest_neighbors_coord_get_cf')
        use, intrinsic :: iso_c_binding
        implicit none
        character(kind=c_char, len = 1) :: f_local_code_name, f_cpl_id, f_field_id
        integer(c_int), value           :: l_local_code_name, l_cpl_id, l_field_id
        integer(c_int), value           :: i_part
        type(c_ptr)                     :: c_nearest_src_coord
        integer(c_int)                  :: n_nearest_src_pts
      end subroutine CWP_Field_nearest_neighbors_coord_get_cf

      subroutine CWP_Field_location_internal_cell_vtx_get_cf(local_code_name,   &
                                                              l_local_code_name, &
                                                              cpl_id,            &
                                                              l_cpl_id,          &
                                                              field_id,          &
                                                              l_field_id,        &
                                                              i_part,            &
                                                              cell_vtx_idx,      &
                                                              n_cell,            &
                                                              cell_vtx)          &

          bind(c, name = 'CWP_Field_location_internal_cell_vtx_get_cf')
          use, intrinsic :: iso_c_binding
          implicit none
          character(kind = c_char, len = 1) :: local_code_name, cpl_id, field_id
          integer(kind = c_int), value      :: l_local_code_name, l_cpl_id, l_field_id
          integer(c_int), value             :: i_part
          type(c_ptr)                       :: cell_vtx_idx, cell_vtx
          integer(kind = c_int) :: n_cell
      end subroutine CWP_Field_location_internal_cell_vtx_get_cf

      subroutine CWP_Spatial_interp_property_set_cf(local_code_name,   &
                                                    l_local_code_name, &
                                                    cpl_id,            &
                                                    l_cpl_id,          &
                                                    property_name,     &
                                                    l_property_name,   &
                                                    property_type,     &
                                                    property_value,    &
                                                    l_property_value)  &
      bind(c, name = 'CWP_Spatial_interp_property_set_cf')
        use, intrinsic :: iso_c_binding
        character(kind = c_char, len = 1) :: local_code_name, cpl_id, property_name, property_value
        integer(kind = c_int), value      :: l_local_code_name, l_cpl_id, l_property_name, property_type, l_property_value
      end subroutine CWP_Spatial_interp_property_set_cf
    !> (DOXYGEN_SHOULD_SKIP_THIS) \endcond

    subroutine CWP_Finalize() &
          bind(c, name = 'CWP_Finalize')
          ! Finalize CWIPI.
    end subroutine CWP_Finalize


    function CWP_Codes_nb_get() &
      result (n_codes)          &
      bind(c, name='CWP_Codes_nb_get')
      ! Return the number of codes known by CWIPI.
      use, intrinsic :: iso_c_binding
      implicit none
      integer(c_int) :: n_codes ! Number of codes
    end function CWP_Codes_nb_get


    function CWP_Loc_codes_nb_get() &
      result (n_local_codes)        &
      bind(c, name='CWP_Loc_codes_nb_get')
      ! Return the number of local codes known by CWIPI.
      use, intrinsic :: iso_c_binding
      implicit none
      integer(c_int) :: n_local_codes ! Number of local codes
    end function CWP_Loc_codes_nb_get


    subroutine CWP_Properties_dump() &
      bind(c, name='CWP_Properties_dump')
      use, intrinsic :: iso_c_binding
      ! Dump code properties.
      implicit none
    end subroutine CWP_Properties_dump

    !> \cond DOXYGEN_SHOULD_SKIP_THIS
    subroutine CWP_Param_add_int_cf(local_code_name,   &
                                    l_local_code_name, &
                                    param_name,        &
                                    l_param_name,      &
                                    initial_value)     &
      bind(c, name='CWP_Param_add_int_cf')
      use, intrinsic :: iso_c_binding
      implicit none
      character(kind = c_char, len = 1) :: local_code_name, param_name
      integer(c_int), value :: l_local_code_name, l_param_name
      integer(c_int), value :: initial_value
    end subroutine CWP_Param_add_int_cf


    subroutine CWP_Param_add_double_cf(local_code_name,   &
                                       l_local_code_name, &
                                       param_name,        &
                                       l_param_name,      &
                                       initial_value)     &
      bind(c, name='CWP_Param_add_double_cf')
      use, intrinsic :: iso_c_binding
      implicit none
      character(kind = c_char, len = 1) :: local_code_name, param_name
      integer(c_int), value :: l_local_code_name, l_param_name
      real(c_double), value :: initial_value
    end subroutine CWP_Param_add_double_cf


    subroutine CWP_Param_add_char_cf(local_code_name,   &
                                     l_local_code_name, &
                                     param_name,        &
                                     l_param_name,      &
                                     initial_value,     &
                                     l_initial_value)   &
      bind(c, name='CWP_Param_add_char_cf')
      use, intrinsic :: iso_c_binding
      implicit none
      character(kind = c_char, len = 1) :: local_code_name, param_name
      integer(c_int), value :: l_local_code_name, l_param_name, l_initial_value
      character(kind = c_char, len = 1) :: initial_value
    end subroutine CWP_Param_add_char_cf

    subroutine CWP_Param_set_int_cf(local_code_name,   &
                                    l_local_code_name, &
                                    param_name,        &
                                    l_param_name,      &
                                    initial_value)     &
      bind(c, name='CWP_Param_set_int_cf')
      use, intrinsic :: iso_c_binding
      implicit none
      character(kind = c_char, len = 1) :: local_code_name, param_name
      integer(c_int), value :: l_local_code_name, l_param_name
      integer(c_int), value :: initial_value
    end subroutine CWP_Param_set_int_cf


    subroutine CWP_Param_set_double_cf(local_code_name,   &
                                       l_local_code_name, &
                                       param_name,        &
                                       l_param_name,      &
                                       initial_value)     &
      bind(c, name='CWP_Param_set_double_cf')
      use, intrinsic :: iso_c_binding
      implicit none
      character(kind = c_char, len = 1) :: local_code_name, param_name
      integer(c_int), value :: l_local_code_name, l_param_name
      real(c_double), value :: initial_value
    end subroutine CWP_Param_set_double_cf


    subroutine CWP_Param_set_char_cf(local_code_name,   &
                                     l_local_code_name, &
                                     param_name,        &
                                     l_param_name,      &
                                     initial_value,     &
                                     l_initial_value)   &
      bind(c, name='CWP_Param_set_char_cf')
      use, intrinsic :: iso_c_binding
      implicit none
      character(kind = c_char, len = 1) :: local_code_name, param_name
      integer(c_int), value :: l_local_code_name, l_param_name, l_initial_value
      character(kind = c_char, len = 1) :: initial_value
    end subroutine CWP_Param_set_char_cf


    subroutine CWP_Param_del_cf(local_code_name,   &
                                l_local_code_name, &
                                param_name,        &
                                l_param_name,      &
                                data_type)         &
      bind(c, name='CWP_Param_del_cf')
      use, intrinsic :: iso_c_binding
      implicit none
      character(kind = c_char, len = 1) :: local_code_name, param_name
      integer(c_int), value :: l_local_code_name, l_param_name
      integer(c_int), value :: data_type
    end subroutine CWP_Param_del_cf


    function CWP_Param_n_get_cf(local_code_name,   &
                                l_local_code_name, &
                                data_type)         &
      result (n_param)                             &
      bind(c, name='CWP_Param_n_get_cf')
      use, intrinsic :: iso_c_binding
      implicit none
      character(kind = c_char, len = 1) :: local_code_name
      integer(c_int), value :: l_local_code_name
      integer(c_int), value :: data_type
      integer(c_int)        :: n_param
    end function CWP_Param_n_get_cf


    function CWP_Param_is_cf(local_code_name,   &
                             l_local_code_name, &
                             param_name,        &
                             l_param_name,      &
                             data_type)         &
      result (is_param)                         &
      bind(c, name='CWP_Param_is_cf')
      use, intrinsic :: iso_c_binding
      implicit none
      character(kind = c_char, len = 1) :: local_code_name, param_name
      integer(c_int), value :: l_local_code_name, l_param_name
      integer(c_int), value :: data_type
      integer(c_int)        :: is_param
    end function CWP_Param_is_cf


    subroutine CWP_Param_get_int_cf(local_code_name,   &
                                    l_local_code_name, &
                                    param_name,        &
                                    l_param_name,      &
                                    value)             &
      bind(c, name='CWP_Param_get_int_cf')
      use, intrinsic :: iso_c_binding
      implicit none
      character(kind = c_char, len = 1) :: local_code_name, param_name
      integer(c_int), value :: l_local_code_name, l_param_name
      integer(c_int)        :: value
    end subroutine CWP_Param_get_int_cf

    subroutine CWP_Param_get_double_cf(local_code_name,   &
                                       l_local_code_name, &
                                       param_name,        &
                                       l_param_name,      &
                                       value)             &
      bind(c, name='CWP_Param_get_double_cf')
      use, intrinsic :: iso_c_binding
      implicit none
      character(kind = c_char, len = 1) :: local_code_name, param_name
      integer(c_int), value :: l_local_code_name, l_param_name
      real(c_double)        :: value
    end subroutine CWP_Param_get_double_cf

    subroutine CWP_Param_get_char_cf(local_code_name,   &
                                     l_local_code_name, &
                                     param_name,        &
                                     l_param_name,      &
                                     val,               &
                                     l_value)           &
      bind(c, name='CWP_Param_get_char_cf')
      use, intrinsic :: iso_c_binding
      implicit none
      character(kind = c_char, len = 1) :: local_code_name, param_name
      integer(c_int), value             :: l_local_code_name, l_param_name
      type(c_ptr)                       :: val
      integer(c_int)                    :: l_value
    end subroutine CWP_Param_get_char_cf


    subroutine CWP_Param_reduce_int_cf(op,           &
                                        param_name,   &
                                        l_param_name, &
                                        res,          &
                                        n_codes,      &
                                        code_names,   &
                                        l_code_names) &

      bind(c, name='CWP_Param_reduce_int_cf')
      use, intrinsic :: iso_c_binding
      implicit none
      integer(c_int), value             :: op, n_codes
      integer(c_int)                    :: res
      character(kind = c_char, len = 1) :: param_name
      integer(c_int), value             :: l_param_name
      type(c_ptr),    value             :: code_names
      type(c_ptr),    value             :: l_code_names
    end subroutine CWP_Param_reduce_int_cf


    subroutine CWP_Param_reduce_double_cf(op,           &
                                          param_name,   &
                                          l_param_name, &
                                          res,          &
                                          n_codes,      &
                                          code_names,   &
                                          l_code_names) &

      bind(c, name='CWP_Param_reduce_double_cf')
      use, intrinsic :: iso_c_binding
      implicit none
      integer(c_int), value             :: op, n_codes
      real(c_double)                    :: res
      character(kind = c_char, len = 1) :: param_name
      integer(c_int), value             :: l_param_name
      type(c_ptr),    value             :: code_names
      type(c_ptr),    value             :: l_code_names
    end subroutine CWP_Param_reduce_double_cf

    subroutine CWP_Param_reduce_char_cf(op,           &
                                        param_name,   &
                                        l_param_name, &
                                        res,          &
                                        l_res,        &
                                        n_codes,      &
                                        code_names,   &
                                        l_code_names) &

      bind(c, name='CWP_Param_reduce_char_cf')
      use, intrinsic :: iso_c_binding
      implicit none
      integer(c_int), value             :: op, n_codes
      type(c_ptr)                       :: res
      integer(c_int)                    :: l_res
      character(kind = c_char, len = 1) :: param_name
      integer(c_int), value             :: l_param_name
      type(c_ptr),    value             :: code_names
      type(c_ptr),    value             :: l_code_names
    end subroutine CWP_Param_reduce_char_cf


    subroutine CWP_Param_lock_cf(local_code_name,   &
                                 l_local_code_name) &
      bind(c, name='CWP_Param_lock_cf')
      use, intrinsic :: iso_c_binding
      implicit none
      character(kind = c_char, len = 1) :: local_code_name
      integer(c_int), value :: l_local_code_name
    end subroutine CWP_Param_lock_cf


    subroutine CWP_Param_unlock_cf(local_code_name,   &
                                   l_local_code_name) &
      bind(c, name='CWP_Param_unlock_cf')
      use, intrinsic :: iso_c_binding
      implicit none
      character(kind = c_char, len = 1) :: local_code_name
      integer(c_int), value :: l_local_code_name
    end subroutine CWP_Param_unlock_cf

    subroutine CWP_Mesh_interf_block_std_get_cf(local_code_name,   &
                                                l_local_code_name, &
                                                cpl_id,            &
                                                l_cpl_id,          &
                                                i_part,            &
                                                block_id,          &
                                                n_elts,            &
                                                c_connec,          &
                                                c_global_num,      &
                                                s_connec)          &
        bind(c, name = 'CWP_Mesh_interf_block_std_get_cf')
        use, intrinsic :: iso_c_binding
        implicit none
        character(kind = c_char, len = 1) :: local_code_name, cpl_id
        integer(c_int), value :: i_part, block_id, n_elts
        type(c_ptr)      :: c_connec, c_global_num
        integer(kind = c_int), value :: l_local_code_name, l_cpl_id
        integer(kind = c_int) :: s_connec
    end subroutine CWP_Mesh_interf_block_std_get_cf


    subroutine CWP_Mesh_interf_f_poly_block_get_cf(local_code_name,   &
                                                   l_local_code_name, &
                                                   cpl_id,            &
                                                   l_cpl_id,          &
                                                   i_part,            &
                                                   block_id,          &
                                                   n_elts,            &
                                                   connec_idx,        &
                                                   connec,            &
                                                   global_num)        &
      bind(c, name='CWP_Mesh_interf_f_poly_block_get_cf')
      use, intrinsic :: iso_c_binding
      implicit none
      character(kind = c_char, len = 1) :: local_code_name, cpl_id
      integer(c_int), value :: l_local_code_name, l_cpl_id
      integer(c_int), value :: i_part, block_id
      integer(c_int)        :: n_elts
      type(c_ptr)           :: connec_idx, connec, global_num
    end subroutine CWP_Mesh_interf_f_poly_block_get_cf


    subroutine CWP_Mesh_interf_c_poly_block_get_cf(local_code_name,   &
                                                   l_local_code_name, &
                                                   cpl_id,            &
                                                   l_cpl_id,          &
                                                   i_part,            &
                                                   block_id,          &
                                                   n_elts,            &
                                                   n_faces,           &
                                                   connec_faces_idx,  &
                                                   connec_faces,      &
                                                   connec_cells_idx,  &
                                                   connec_cells,      &
                                                   global_num)        &
      bind(c, name='CWP_Mesh_interf_c_poly_block_get_cf')
      use, intrinsic :: iso_c_binding
      implicit none
      character(kind = c_char, len = 1) :: local_code_name, cpl_id
      integer(c_int), value :: l_local_code_name, l_cpl_id
      integer(c_int), value :: i_part, block_id
      integer(c_int)        :: n_elts, n_faces
      type(c_ptr)           :: connec_faces_idx, connec_faces, connec_cells_idx, connec_cells, global_num
    end subroutine CWP_Mesh_interf_c_poly_block_get_cf


    function CWP_Field_dof_location_get_cf(local_code_name,   &
                                           l_local_code_name, &
                                           cpl_id,            &
                                           l_cpl_id,          &
                                           field_id,          &
                                           l_field_id)        &
      result (dof_location)                                   &
      bind(c, name='CWP_Field_dof_location_get_cf')
      use, intrinsic :: iso_c_binding
      implicit none
      character(kind = c_char, len = 1) :: local_code_name, cpl_id, field_id
      integer(c_int), value :: l_local_code_name, l_cpl_id, l_field_id
      integer(c_int)        :: dof_location
    end function CWP_Field_dof_location_get_cf


    function CWP_Field_storage_get_cf(local_code_name,   &
                                      l_local_code_name, &
                                      cpl_id,            &
                                      l_cpl_id,          &
                                      field_id,          &
                                      l_field_id)        &
      result (storage)                                   &
      bind(c, name='CWP_Field_storage_get_cf')
      use, intrinsic :: iso_c_binding
      implicit none
      character(kind = c_char, len = 1) :: local_code_name, cpl_id, field_id
      integer(c_int), value :: l_local_code_name, l_cpl_id, l_field_id
      integer(c_int)        :: storage
    end function CWP_Field_storage_get_cf

    function CWP_Field_n_dof_get_cf(local_code_name,   &
                                    l_local_code_name, &
                                    cpl_id,            &
                                    l_cpl_id,          &
                                    field_id,          &
                                    l_field_id,        &
                                    i_part)            &
      result (n_dof)                                     &
      bind(c, name='CWP_Field_n_dof_get_cf')
      use, intrinsic :: iso_c_binding
      implicit none
      character(kind = c_char, len = 1) :: local_code_name, cpl_id, field_id
      integer(c_int), value :: l_local_code_name, l_cpl_id, l_field_id
      integer(c_int), value :: i_part
      integer(c_int)        :: n_dof
    end function CWP_Field_n_dof_get_cf

    subroutine CWP_Field_del_cf(local_code_name,   &
                                l_local_code_name, &
                                cpl_id,            &
                                l_cpl_id,          &
                                field_id,          &
                                l_field_id)        &
      bind(c, name='CWP_Field_del_cf')
      use, intrinsic :: iso_c_binding
      implicit none
      character(kind = c_char, len = 1) :: local_code_name, cpl_id, field_id
      integer(c_int), value :: l_local_code_name, l_cpl_id, l_field_id
    end subroutine CWP_Field_del_cf

    subroutine CWP_Codes_list_get_cf(code_list, code_list_s, n_codes) &
      bind (c, name="CWP_Codes_list_get_cf")
      use, intrinsic :: iso_c_binding
      implicit none
      type(c_ptr)            :: code_list, code_list_s
      integer(c_int)         :: n_codes
    end subroutine CWP_Codes_list_get_cf

    subroutine CWP_Loc_codes_list_get_cf(loc_code_list, loc_code_list_s, n_loc_codes) &
      bind (c, name="CWP_Loc_codes_list_get_cf")
      use, intrinsic :: iso_c_binding
      implicit none
      type(c_ptr)            :: loc_code_list, loc_code_list_s
      integer(c_int)         :: n_loc_codes
    end subroutine CWP_Loc_codes_list_get_cf

    subroutine CWP_Param_list_get_cf(code_name, l_code_name, data_type, n_param, c_param_names, c_param_sizes) &
      bind (c, name="CWP_Param_list_get_cf")
      use, intrinsic :: iso_c_binding
      implicit none
      type(c_ptr)                       :: c_param_names, c_param_sizes
      integer(c_int)                    :: n_param
      character(kind = c_char, len = 1) :: code_name
      integer(c_int), value             :: l_code_name
      integer(c_int), value             :: data_type
    end subroutine CWP_Param_list_get_cf

    subroutine CWP_Global_data_issend_cf(f_local_code_name, &
                                         l_local_code_name, &
                                         f_cpl_id,          &
                                         l_cpl_id,          &
                                         f_global_data_id,  &
                                         l_global_data_id,  &
                                         s_send_entity,     &
                                         send_stride,       &
                                         n_send_entity,     &
                                         send_data)         &
    bind (c, name='CWP_Global_data_issend_cf')
      use, intrinsic :: iso_c_binding
      implicit none
      character(kind=c_char, len=1) :: f_local_code_name, f_cpl_id, f_global_data_id
      integer(c_int),  value        :: l_local_code_name, l_cpl_id, l_global_data_id
      integer(c_long), value        :: s_send_entity
      integer(c_int),  value        :: send_stride
      integer(c_int),  value        :: n_send_entity
      type(c_ptr),     value        :: send_data
    end subroutine CWP_Global_data_issend_cf

    subroutine CWP_Global_data_irecv_cf(f_local_code_name,      &
                                        l_local_code_name,      &
                                        f_cpl_id,               &
                                        l_cpl_id,               &
                                        f_global_data_id,       &
                                        l_global_data_id,       &
                                        s_recv_entity,          &
                                        recv_stride,            &
                                        n_recv_entity,          &
                                        recv_data)              &
    bind (c, name='CWP_Global_data_irecv_cf')
      use, intrinsic :: iso_c_binding
      implicit none
      character(kind=c_char, len=1) :: f_local_code_name, f_cpl_id, f_global_data_id
      integer(c_int),  value        :: l_local_code_name, l_cpl_id, l_global_data_id
      integer(c_long), value        :: s_recv_entity
      integer(c_int),  value        :: recv_stride
      integer(c_int),  value        :: n_recv_entity
      type(c_ptr),     value        :: recv_data
    end subroutine CWP_Global_data_irecv_cf

    subroutine CWP_Global_data_wait_issend_cf(f_local_code_name, &
                                              l_local_code_name, &
                                              f_cpl_id,          &
                                              l_cpl_id,          &
                                              f_global_data_id,  &
                                              l_global_data_id)  &
    bind (c, name='CWP_Global_data_wait_issend_cf')
      use, intrinsic :: iso_c_binding
      implicit none
      character(kind=c_char, len=1) :: f_local_code_name, f_cpl_id, f_global_data_id
      integer(c_int), value         :: l_local_code_name, l_cpl_id, l_global_data_id
    end subroutine CWP_Global_data_wait_issend_cf

    subroutine CWP_Global_data_wait_irecv_cf(f_local_code_name, &
                                             l_local_code_name, &
                                             f_cpl_id,          &
                                             l_cpl_id,          &
                                             f_global_data_id,  &
                                             l_global_data_id)  &
    bind (c, name='CWP_Global_data_wait_irecv_cf')
      use, intrinsic :: iso_c_binding
      implicit none
      character(kind=c_char, len=1) :: f_local_code_name, f_cpl_id, f_global_data_id
      integer(c_int), value         :: l_local_code_name, l_cpl_id, l_global_data_id
    end subroutine CWP_Global_data_wait_irecv_cf

    subroutine CWP_Part_data_create_cf(f_local_code_name, &
                                       l_local_code_name, &
                                       f_cpl_id,          &
                                       l_cpl_id,          &
                                       f_part_data_id,    &
                                       l_part_data_id,    &
                                       exch_type,         &
                                       gnum_elt,          &
                                       n_elt,             &
                                       n_part)            &
    bind (c, name='CWP_Part_data_create_cf')
      use, intrinsic :: iso_c_binding
      implicit none
      character(kind=c_char, len=1) :: f_local_code_name, f_cpl_id, f_part_data_id
      integer(c_int), value         :: l_local_code_name, l_cpl_id, l_part_data_id
      integer(c_int), value         :: exch_type
      type(c_ptr),    value         :: gnum_elt
      type(c_ptr),    value         :: n_elt
      integer(c_int), value         :: n_part
    end subroutine CWP_Part_data_create_cf


    subroutine CWP_Part_data_del_cf(f_local_code_name, &
                                    l_local_code_name, &
                                    f_cpl_id,          &
                                    l_cpl_id,          &
                                    f_part_data_id,    &
                                    l_part_data_id)    &
    bind (c, name='CWP_Part_data_del_cf')
      use, intrinsic :: iso_c_binding
      implicit none
      character(kind=c_char, len=1) :: f_local_code_name, f_cpl_id, f_part_data_id
      integer(c_int), value         :: l_local_code_name, l_cpl_id, l_part_data_id
    end subroutine CWP_Part_data_del_cf

    subroutine CWP_Part_data_issend_cf(f_local_code_name,   &
                                       l_local_code_name,   &
                                       f_cpl_id,            &
                                       l_cpl_id,            &
                                       f_part_data_id,      &
                                       l_part_data_id,      &
                                       exch_id,             &
                                       s_data,              &
                                       n_components,        &
                                       send_data)           &
    bind (c, name='CWP_Part_data_issend_cf')
      use, intrinsic :: iso_c_binding
      implicit none
      character(kind=c_char, len=1) :: f_local_code_name, f_cpl_id, f_part_data_id
      integer(c_int),  value        :: l_local_code_name, l_cpl_id, l_part_data_id
      integer(c_int),  value        :: exch_id
      integer(c_long), value        :: s_data
      integer(c_int),  value        :: n_components
      type(c_ptr),     value        :: send_data
    end subroutine CWP_Part_data_issend_cf

    subroutine CWP_Part_data_irecv_cf(f_local_code_name,   &
                                      l_local_code_name,   &
                                      f_cpl_id,            &
                                      l_cpl_id,            &
                                      f_part_data_id,      &
                                      l_part_data_id,      &
                                      exch_id,             &
                                      s_data,              &
                                      n_components,        &
                                      recv_data)           &
    bind (c, name='CWP_Part_data_irecv_cf')
      use, intrinsic :: iso_c_binding
      implicit none
      character(kind=c_char, len=1) :: f_local_code_name, f_cpl_id, f_part_data_id
      integer(c_int),  value        :: l_local_code_name, l_cpl_id, l_part_data_id
      integer(c_int),  value        :: exch_id
      integer(c_long), value        :: s_data
      integer(c_int),  value        :: n_components
      type(c_ptr),     value        :: recv_data
    end subroutine CWP_Part_data_irecv_cf

    subroutine CWP_Part_data_wait_issend_cf(f_local_code_name,   &
                                            l_local_code_name,   &
                                            f_cpl_id,            &
                                            l_cpl_id,            &
                                            f_part_data_id,      &
                                            l_part_data_id,      &
                                            exch_id)             &
    bind (c, name='CWP_Part_data_wait_issend_cf')
      use, intrinsic :: iso_c_binding
      implicit none
      character(kind=c_char, len=1) :: f_local_code_name, f_cpl_id, f_part_data_id
      integer(c_int), value         :: l_local_code_name, l_cpl_id, l_part_data_id
      integer(c_int),  value        :: exch_id
    end subroutine CWP_Part_data_wait_issend_cf

    subroutine CWP_Part_data_wait_irecv_cf(f_local_code_name,   &
                                           l_local_code_name,   &
                                           f_cpl_id,            &
                                           l_cpl_id,            &
                                           f_part_data_id,      &
                                           l_part_data_id,      &
                                           exch_id)             &
    bind (c, name='CWP_Part_data_wait_irecv_cf')
      use, intrinsic :: iso_c_binding
      implicit none
      character(kind=c_char, len=1) :: f_local_code_name, f_cpl_id, f_part_data_id
      integer(c_int), value         :: l_local_code_name, l_cpl_id, l_part_data_id
      integer(c_int),  value        :: exch_id
    end subroutine CWP_Part_data_wait_irecv_cf

    subroutine CWP_Part_data_n_part_get_cf(f_local_code_name, &
                                           l_local_code_name, &
                                           f_cpl_id,          &
                                           l_cpl_id,          &
                                           f_part_data_id,    &
                                           l_part_data_id,    &
                                           n_part)            &
    bind (c, name='CWP_Part_data_n_part_get_cf')
      use, intrinsic :: iso_c_binding
      implicit none
      character(kind=c_char, len=1) :: f_local_code_name, f_cpl_id, f_part_data_id
      integer(c_int), value         :: l_local_code_name, l_cpl_id, l_part_data_id
      integer(c_int)                :: n_part
    end subroutine CWP_Part_data_n_part_get_cf

    subroutine CWP_Part_data_n_ref_get_cf(f_local_code_name, &
                                          l_local_code_name, &
                                          f_cpl_id,          &
                                          l_cpl_id,          &
                                          f_part_data_id,    &
                                          l_part_data_id,    &
                                          i_part,            &
                                          n_ref)             &
    bind (c, name='CWP_Part_data_n_ref_get_cf')
      use, intrinsic :: iso_c_binding
      implicit none
      character(kind=c_char, len=1) :: f_local_code_name, f_cpl_id, f_part_data_id
      integer(c_int), value         :: l_local_code_name, l_cpl_id, l_part_data_id
      integer(c_int), value         :: i_part
      integer(c_int)                :: n_ref
    end subroutine CWP_Part_data_n_ref_get_cf


    function CWP_Cpl_spatial_interp_algo_get_cf(f_local_code_name,      &
                                                l_local_code_name,      &
                                                f_cpl_id,               &
                                                l_cpl_id) result (algo) &
    bind (c, name='CWP_Cpl_spatial_interp_algo_get_cf')
      use, intrinsic :: iso_c_binding
      implicit none
      character(kind=c_char, len=1) :: f_local_code_name, f_cpl_id
      integer(c_int), value         :: l_local_code_name, l_cpl_id
      integer(c_int)                :: algo
    end function CWP_Cpl_spatial_interp_algo_get_cf
    !> (DOXYGEN_SHOULD_SKIP_THIS) \endcond

  end interface

contains

  !> convert an array of char * in c to an array
  !  of characters of maximum size 256 in FORTRAN

  subroutine c_f_char_array(c_char_array, &
                            c_size_array, &
                            n_chars, &
                            f_char_array, &
                            free_all)

    use, intrinsic :: iso_c_binding
    implicit none

    type(c_ptr) :: c_char_array, c_size_array
    type(c_ptr), pointer :: fptr2(:) => null()
    character(c_char), pointer :: fptr(:) => null()
    integer(c_int), pointer     :: f_size_array(:)
    character(256), allocatable :: f_char_array(:)
    logical, optional           :: free_all
    integer(c_int) :: i, n_chars, strlen

    call c_f_pointer(c_char_array, fptr2, [n_chars])
    call c_f_pointer(c_size_array, f_size_array, [n_chars])

    allocate(f_char_array(n_chars))
    do i = 1, n_chars
      strlen      = f_size_array(i)
      call c_f_pointer(fptr2(i), fptr, [strlen])
      f_char_array(i) = transfer(fptr(1:strlen), f_char_array(i))
      f_char_array(i) = f_char_array(i)(1:strlen)
      if (present(free_all)) then
        if (free_all) then
          call pdm_fortran_free_c(fptr2(i))
        endif
      endif

    end do

    call pdm_fortran_free_c(c_char_array)
    call pdm_fortran_free_c(c_size_array)

  end subroutine c_f_char_array


  subroutine CWP_Init_(f_comm,         &
                       n_code,         &
                       code_names,     &
                       is_active_rank, &
                       intra_comms)
    ! Initialize CWIPI.
    !
    ! This function creates the MPI intra communicators of the codes from
    ! the \p global_comm MPI communicator that contains all code ranks. This
    ! function has to be called from all ranks contained in the ``global_comm``.

    use, intrinsic :: iso_c_binding
    implicit none

    integer(c_int)                                               :: f_comm         ! MPI global communicator
    integer(c_int), intent(in)                                   :: n_code         ! Number of codes on the current rank
    character(kind = c_char, len = *), dimension(n_code), target :: code_names     ! Names of codes on the current rank (size = ``n_code``)
    integer(c_int)                                               :: is_active_rank ! Current rank is available for CWIPI
    integer(c_int), dimension(:), pointer                        :: intra_comms    ! MPI intra communicators of each code (size = ``n_code``)

    integer, dimension(n_code), target :: l_code_names
    integer :: i

    do i = 1, n_code
      l_code_names(i) = len(code_names(i))
    end do

    call CWP_Init_cf(f_comm,              &
                     n_code,              &
                     c_loc(code_names),   &
                     c_loc(l_code_names), &
                     is_active_rank,      &
                     c_loc(intra_comms))

  end subroutine CWP_Init_


  function CWP_C_to_f_string_(c_str)  result(f_str)
    ! Create a Fortran string from a C string
    !
    ! This function creates a Fortran string from a *copy* of the C string.
    use iso_c_binding

    implicit none

    character(kind=c_char,len=1), intent(in) :: c_str(*) ! C string
    character(len=:), pointer                :: f_str    ! Fortran string
    integer i, nchars


    i = 1
    do
       if (c_str(i) == c_null_char) exit
       i = i + 1
    end do
    nchars = i - 1  ! Exclude null character from Fortran string
    allocate(character(len=nchars) :: f_str)

    f_str = transfer(c_str(1:nchars), f_str)
  end function


  subroutine CWP_State_update_(local_code_name, &
                               state)

    ! Update code state.
    use, intrinsic :: iso_c_binding
    implicit none

    character(kind = c_char, len = *) :: local_code_name ! Local code name
    integer(kind = c_int), intent(in) :: state           ! State
    integer(c_int)                    :: l_local_code_name

    l_local_code_name = len(local_code_name)

    call CWP_State_update_cf(local_code_name,   &
                             l_local_code_name, &
                             state)

  end subroutine CWP_State_update_


  subroutine CWP_Time_step_beg_(local_code_name, &
                             current_time)
    ! Begin code time step.
    use, intrinsic :: iso_c_binding
    implicit none

    character(kind = c_char, len = *) :: local_code_name ! Local code name
    real(8), intent(in)               :: current_time    ! Current time
    integer(c_int)                    :: l_local_code_name

    l_local_code_name = len(local_code_name)

    call CWP_Time_step_beg_cf(local_code_name,   &
                            l_local_code_name, &
                            current_time)

  end subroutine CWP_Time_step_beg_


  subroutine CWP_Time_step_end_(local_code_name)
    ! End code time step.
    use, intrinsic :: iso_c_binding
    implicit none

    character(kind = c_char, len = *) :: local_code_name ! Local code name
    integer(c_int)                    :: l_local_code_name

    l_local_code_name = len(local_code_name)

    call CWP_Time_step_end_cf(local_code_name,   &
                              l_local_code_name)

  end subroutine CWP_Time_step_end_

  subroutine CWP_User_structure_set_(local_code_name, &
                                    user_structure)
    ! Define a user structure associated to a code.
    !
    ! This structure can be accessed into a callback.
    use, intrinsic :: iso_c_binding
    implicit none

    character(kind = c_char, len = *) :: local_code_name ! Local code name
    type(c_ptr), value                :: user_structure  ! User structure
    integer(c_int)                    :: l_local_code_name

    l_local_code_name = len(local_code_name)

    call CWP_User_structure_set_cf(local_code_name,   &
                                 l_local_code_name, &
                                 user_structure)

  end subroutine CWP_User_structure_set_


  function CWP_User_structure_get_(local_code_name) &
    result (user_structure)
    ! Return the user structure associated to a code.
    !
    ! This structure can be accessed into a callback.
    use, intrinsic :: iso_c_binding
    implicit none

    character(kind = c_char, len = *) :: local_code_name ! Local code name
    type(c_ptr)                       :: user_structure  ! User structure
    integer(c_int)                    :: l_local_code_name

    l_local_code_name = len(local_code_name)

    user_structure = CWP_User_structure_get_cf(local_code_name,   &
                                               l_local_code_name)

  end function CWP_User_structure_get_


  subroutine CWP_output_fortran_unit_set (outputUnit)
    ! Writing output to Fortran file (shared by fortran and C code).
    !
    ! This function set the file Fortran logical unit for writing output.
    use, intrinsic :: iso_c_binding
    use cwp_printfort

    implicit none

    integer :: outputUnit ! File Fortan logical unit

    ifile = outputUnit

    call CWP_set_output_listing_f(outputUnit)

  end subroutine CWP_output_fortran_unit_set


  subroutine CWP_Output_file_set_ (f_output_file_name)
    ! Define output file (in which only C code writes).
    use, intrinsic :: iso_c_binding
    implicit none

    character(kind = c_char, len = *) :: f_output_file_name ! Output file name
    integer(c_int)                    :: l_output_file_name

    l_output_file_name = len(f_output_file_name)

    call CWP_Output_file_set_cf (f_output_file_name, &
                                 l_output_file_name)
  end subroutine CWP_Output_file_set_


  function CWP_State_get_(local_code_name) &
    result (state)
    ! Return code state.
    use, intrinsic :: iso_c_binding
    implicit none

    character(kind = c_char, len = *) :: local_code_name ! Code name
    integer(c_int)                    :: state           ! Code state
    integer(c_int)                    :: l_local_code_name

    l_local_code_name = len(local_code_name)

    state = CWP_State_get_cf(local_code_name,   &
                             l_local_code_name)

  end function CWP_State_get_


  function CWP_Codes_list_get_() &
    result (fstrings)
    ! Return list of codes known by CWIPI.
    use, intrinsic :: iso_c_binding
    implicit none

    type(c_ptr) :: code_list, code_list_s
    integer(c_int) :: n_codes
    character(256), allocatable :: fstrings(:) ! List of code names

    call CWP_Codes_list_get_cf(code_list, code_list_s, n_codes)

    call c_f_char_array(code_list, code_list_s, n_codes, fstrings)

  end function CWP_Codes_list_get_


  function CWP_Loc_codes_list_get_() &
    result (fstrings)
    ! Return list of local codes known by CWIPI.
    use, intrinsic :: iso_c_binding
    implicit none

    type(c_ptr) :: loc_code_list, loc_code_list_s
    integer(c_int) :: n_loc_codes
    character(256), allocatable :: fstrings(:) ! List of local code names

    call CWP_Loc_codes_list_get_cf(loc_code_list, loc_code_list_s, n_loc_codes)

    call c_f_char_array(loc_code_list, loc_code_list_s, n_loc_codes, fstrings)

  end function CWP_Loc_codes_list_get_


  subroutine CWP_Cpl_create_(local_code_name,   &
                             cpl_id,            &
                             coupled_code_name, &
                             entities_dim,      &
                             comm_type,         &
                             spatial_interp,    &
                             n_part,            &
                             displacement,      &
                             freq)
    ! Create a coupling object and define its properties.
    use, intrinsic :: iso_c_binding
    implicit none

    character(kind = c_char, len = *) :: local_code_name   ! Local code name
    character(kind = c_char, len = *) :: cpl_id            ! Coupling identifier
    character(kind = c_char, len = *) :: coupled_code_name ! Distant or local coupled code name
    integer(kind = c_int), intent(in) :: entities_dim      ! Coupling interface type
    integer(kind = c_int), intent(in) :: comm_type         ! Communication type
    integer(kind = c_int), intent(in) :: spatial_interp    ! :ref:`Spatial interpolation method <spatial_interp>`
    integer(kind = c_int), intent(in) :: n_part            ! Number of interface partition
    integer(kind = c_int), intent(in) :: displacement      ! Mesh moving status
    integer(kind = c_int), intent(in) :: freq              ! Type of receiving frequency
    integer(c_int) :: l_local_code_name, l_cpl_id, l_coupled_code_name

    l_local_code_name   = len(local_code_name)
    l_cpl_id            = len(cpl_id)
    l_coupled_code_name = len(coupled_code_name)

    call CWP_Cpl_create_cf (local_code_name,    &
                            l_local_code_name,  &
                            cpl_id,             &
                            l_cpl_id,           &
                            coupled_code_name,  &
                            l_coupled_code_name,&
                            entities_dim,       &
                            comm_type,          &
                            spatial_interp,     &
                            n_part,             &
                            displacement,       &
                            freq)
  end subroutine CWP_Cpl_Create_


  subroutine CWP_Cpl_barrier_(local_code_name, &
                              cpl_id)
    ! MPI Barrier on the coupling communicator.
    use, intrinsic :: iso_c_binding
    implicit none

    character(kind = c_char, len = *) :: local_code_name ! Local code name
    character(kind = c_char, len = *) :: cpl_id          ! Coupling identifier
    integer(c_int)                    :: l_local_code_name, l_cpl_id

    l_local_code_name   = len(local_code_name)
    l_cpl_id            = len(cpl_id)

    call CWP_Cpl_barrier_cf(local_code_name,   &
                            l_local_code_name, &
                            cpl_id,            &
                            l_cpl_id)
  end subroutine CWP_Cpl_barrier_


  subroutine CWP_Cpl_Del_ (local_code_name, &
                          cpl_id)
    ! Delete a coupling object.
    use, intrinsic :: iso_c_binding
    implicit none

    character(kind = c_char, len = *) :: local_code_name ! Local code name
    character(kind = c_char, len = *) :: cpl_id          ! Coupling identifier
    integer(c_int) :: l_local_code_name, l_cpl_id

    l_local_code_name = len(local_code_name)
    l_cpl_id          = len(cpl_id)

    call CWP_Cpl_del_cf (local_code_name, l_local_code_name, cpl_id, l_cpl_id)
  end subroutine CWP_Cpl_Del_


  subroutine CWP_Computed_tgts_bcast_enable_ (local_code_name, &
                                              cpl_id,          &
                                              field_id)
    ! Enable broadcast of the computed targets ids (in \ref CWP_COMM_PAR_WITHOUT_PART mode).
    !
    ! This function must be called in order for the computed targets to be accessible
    ! on non-root ranks

    use, intrinsic :: iso_c_binding
    implicit none

    character(kind = c_char, len = *) :: local_code_name ! Local code name
    character(kind = c_char, len = *) :: cpl_id          ! Coupling identifier
    character(kind = c_char, len = *) :: field_id        ! Current partition
    integer(c_int) :: l_local_code_name, l_cpl_id, l_field_id

    l_local_code_name = len(local_code_name)
    l_cpl_id          = len(cpl_id)
    l_field_id        = len(field_id)

    call CWP_Computed_tgts_bcast_enable_cf (local_code_name, l_local_code_name, cpl_id, l_cpl_id, field_id, l_field_id)
  end subroutine CWP_Computed_tgts_bcast_enable_


  subroutine CWP_Involved_srcs_bcast_enable_ (local_code_name, &
                                              cpl_id,          &
                                              field_id)
    ! Enable broadcast of the involved sources ids (in \ref CWP_COMM_PAR_WITHOUT_PART mode).
    !
    ! This function must be called in order for the involved sources to be accessible
    ! on non-root ranks
    use, intrinsic :: iso_c_binding
    implicit none

    character(kind = c_char, len = *) :: local_code_name ! Local code name
    character(kind = c_char, len = *) :: cpl_id          ! Coupling identifier
    character(kind = c_char, len = *) :: field_id        ! Current partition
    integer(c_int) :: l_local_code_name, l_cpl_id, l_field_id

    l_local_code_name = len(local_code_name)
    l_cpl_id          = len(cpl_id)
    l_field_id        = len(field_id)

    call CWP_Involved_srcs_bcast_enable_cf (local_code_name, l_local_code_name, cpl_id, l_cpl_id, field_id, l_field_id)
  end subroutine CWP_Involved_srcs_bcast_enable_


  function CWP_N_uncomputed_tgts_get_(local_code_name, &
                                      cpl_id,          &
                                      field_id,        &
                                      i_part)          &
  result (n_uncomputed_tgts)
    ! Return the number of uncomputed targets.
    use, intrinsic :: iso_c_binding
    implicit none

    character(kind = c_char, len = *) :: local_code_name   ! Local code name
    character(kind = c_char, len = *) :: cpl_id            ! Coupling identifier
    character(kind = c_char, len = *) :: field_id          ! Field identifier
    integer(c_int)                    :: i_part            ! Current partition

    integer(c_int)                    :: n_uncomputed_tgts ! Number of uncomputed targets

    integer(c_int) :: l_local_code_name, l_cpl_id, l_field_id

    l_local_code_name = len(local_code_name)
    l_cpl_id = len(cpl_id)
    l_field_id = len(field_id)

    n_uncomputed_tgts = CWP_N_uncomputed_tgts_get_cf (local_code_name,   &
                                                      l_local_code_name, &
                                                      cpl_id,            &
                                                      l_cpl_id,          &
                                                      field_id,          &
                                                      l_field_id,        &
                                                      i_part)
  end function CWP_N_uncomputed_tgts_get_


  function CWP_Uncomputed_tgts_get_(local_code_name, &
                                    cpl_id,          &
                                    field_id,        &
                                    i_part)          &
  result (uncomputed_tgts)
    ! Return uncomputed targets.
    use, intrinsic :: iso_c_binding
    implicit none

    character(kind = c_char, len = *)     :: local_code_name ! Local code name
    character(kind = c_char, len = *)     :: cpl_id          ! Coupling identifier
    character(kind = c_char, len = *)     :: field_id        ! Field identifier
    integer(c_int)                        :: i_part          ! Current partition

    integer(c_int), dimension(:), pointer :: uncomputed_tgts ! Uncomputed targets

    type(c_ptr) :: cptr_uncomputed_tgts
    integer(c_int) :: n_uncomputed_tgts
    integer(c_int) :: l_local_code_name, l_cpl_id, l_field_id

    l_local_code_name = len(local_code_name)
    l_cpl_id = len(cpl_id)
    l_field_id = len(field_id)

    n_uncomputed_tgts = CWP_N_uncomputed_tgts_get_cf (local_code_name,   &
                                                      l_local_code_name, &
                                                      cpl_id,            &
                                                      l_cpl_id,          &
                                                      field_id,          &
                                                      l_field_id,        &
                                                      i_part)

    cptr_uncomputed_tgts = CWP_Uncomputed_tgts_get_cf(local_code_name,   &
                                                      l_local_code_name, &
                                                      cpl_id,            &
                                                      l_cpl_id,          &
                                                      field_id,          &
                                                      l_field_id,        &
                                                      i_part)


    call c_f_pointer (cptr = cptr_uncomputed_tgts, &
                      fptr = uncomputed_tgts  ,    &
                      shape= [ n_uncomputed_tgts])
  end function CWP_Uncomputed_tgts_get_


  function CWP_N_computed_tgts_get_(local_code_name, &
                                    cpl_id,          &
                                    field_id,        &
                                    i_part)          &
                                    result (n_computed_tgts)
    ! Return the number of computed targets.
    use, intrinsic :: iso_c_binding
    implicit none

    character(kind = c_char, len = *) :: local_code_name ! Local code name
    character(kind = c_char, len = *) :: cpl_id          ! Coupling identifier
    character(kind = c_char, len = *) :: field_id        ! Field identifier
    integer(c_int)                    :: i_part          ! Current partition

    integer(c_int)                    :: n_computed_tgts ! Number of computed targets

    integer(c_int) :: l_local_code_name, l_cpl_id, l_field_id

    l_local_code_name = len(local_code_name)
    l_cpl_id = len(cpl_id)
    l_field_id = len(field_id)

    n_computed_tgts = CWP_N_computed_tgts_get_cf (local_code_name,   &
                                                  l_local_code_name, &
                                                  cpl_id,            &
                                                  l_cpl_id,          &
                                                  field_id,          &
                                                  l_field_id,        &
                                                  i_part)
  end function CWP_N_computed_tgts_get_


  function CWP_Computed_tgts_get_(local_code_name, &
                                  cpl_id,          &
                                  field_id,        &
                                  i_part)          &
                                  result (computed_tgts)
    ! Return computed targets.
    use, intrinsic :: iso_c_binding
    implicit none


    character(kind = c_char, len = *)     :: local_code_name ! Local code name
    character(kind = c_char, len = *)     :: cpl_id          ! Coupling identifier
    character(kind = c_char, len = *)     :: field_id        ! Field identifier
    integer(c_int)                        :: i_part          ! Current partition

    integer(c_int), dimension(:), pointer :: computed_tgts   ! Computed targets

    type(c_ptr) :: cptr_computed_tgts
    integer(c_int) :: n_computed_tgts
    integer(c_int) :: l_local_code_name, l_cpl_id, l_field_id

    l_local_code_name = len(local_code_name)
    l_cpl_id = len(cpl_id)
    l_field_id = len(field_id)


    n_computed_tgts = CWP_N_computed_tgts_get_cf (local_code_name,   &
                                                  l_local_code_name, &
                                                  cpl_id,            &
                                                  l_cpl_id,          &
                                                  field_id,          &
                                                  l_field_id,        &
                                                  i_part)

    cptr_computed_tgts = CWP_Computed_tgts_get_cf(local_code_name,   &
                                                  l_local_code_name, &
                                                  cpl_id,            &
                                                  l_cpl_id,          &
                                                  field_id,          &
                                                  l_field_id,        &
                                                  i_part)


    call c_f_pointer (cptr = cptr_computed_tgts, &
                      fptr = computed_tgts  ,    &
                      shape= [n_computed_tgts])
  end function CWP_Computed_tgts_get_


  function CWP_N_involved_srcs_get_(local_code_name, &
                                    cpl_id,          &
                                    field_id,        &
                                    i_part)          &
  result (n_involved_srcs)
    ! Return the number of involved sources.
    use, intrinsic :: iso_c_binding
    implicit none

    character(kind = c_char, len = *) :: local_code_name ! Local code name
    character(kind = c_char, len = *) :: cpl_id          ! Coupling identifier
    character(kind = c_char, len = *) :: field_id        ! Field identifier
    integer(c_int)                    :: i_part          ! Current partition

    integer(c_int)                    :: n_involved_srcs ! Number of involved sources

    integer(c_int) :: l_local_code_name, l_cpl_id, l_field_id

    l_local_code_name = len(local_code_name)
    l_cpl_id = len(cpl_id)
    l_field_id = len(field_id)

    n_involved_srcs = CWP_N_involved_srcs_get_cf (local_code_name,   &
                                                  l_local_code_name, &
                                                  cpl_id,            &
                                                  l_cpl_id,          &
                                                  field_id,          &
                                                  l_field_id,        &
                                                  i_part)
  end function CWP_N_involved_srcs_get_


  function CWP_Involved_srcs_get_(local_code_name, &
                                  cpl_id,          &
                                  field_id,        &
                                  i_part)          &
  result (involved_srcs)
    ! Return involved sources.
    use, intrinsic :: iso_c_binding
    implicit none


    character(kind = c_char, len = *)     :: local_code_name ! Local code name
    character(kind = c_char, len = *)     :: cpl_id          ! Coupling identifier
    character(kind = c_char, len = *)     :: field_id        ! Field identifier
    integer(c_int)                        :: i_part          ! Current partition

    integer(c_int), dimension(:), pointer :: involved_srcs   ! Involved sources

    type(c_ptr) :: cptr_involved_srcs
    integer(c_int) :: n_involved_srcs
    integer(c_int) :: l_local_code_name, l_cpl_id, l_field_id

    l_local_code_name = len(local_code_name)
    l_cpl_id = len(cpl_id)
    l_field_id = len(field_id)


    n_involved_srcs = CWP_N_involved_srcs_get_cf (local_code_name,   &
            l_local_code_name, &
            cpl_id,            &
            l_cpl_id,          &
            field_id,          &
            l_field_id,        &
            i_part)

    cptr_involved_srcs = CWP_Involved_srcs_get_cf(local_code_name,   &
            l_local_code_name, &
            cpl_id,            &
            l_cpl_id,          &
            field_id,          &
            l_field_id,        &
            i_part)


    call c_f_pointer (cptr = cptr_involved_srcs, &
            fptr = involved_srcs  ,    &
            shape= [n_involved_srcs])
  end function CWP_Involved_srcs_get_


  subroutine CWP_Spatial_interp_weights_compute_(local_code_name, &
                                                 cpl_id)
    ! Compute spatial interpolation weights.
    use, intrinsic :: iso_c_binding
    implicit none

    character(kind = c_char, len = *) :: local_code_name ! Local code name
    character(kind = c_char, len = *) :: cpl_id          ! Coupling identifier
    integer(kind = c_int) :: l_local_code_name, l_cpl_id

    l_local_code_name = len(local_code_name)
    l_cpl_id = len(cpl_id)

    call CWP_Spatial_interp_weights_compute_cf(local_code_name, l_local_code_name, cpl_id, l_cpl_id)
  end subroutine CWP_Spatial_interp_weights_compute_


  subroutine CWP_Spatial_interp_property_set_(local_code_name, &
                                              cpl_id,          &
                                              property_name,   &
                                              property_type,   &
                                              property_value)
    ! Set a property of the :ref:`spatial interpolation algorithm <spatial_interp>`.
    use, intrinsic :: iso_c_binding
    implicit none

    character(kind = c_char, len = *) :: local_code_name ! Local code name
    character(kind = c_char, len = *) :: cpl_id          ! Coupling identifier
    character(kind = c_char, len = *) :: property_name   ! Name of the property
    integer  (kind = c_int)           :: property_type   ! Type of the property
    character(kind = c_char, len = *) :: property_value  ! Value of the property

    call CWP_Spatial_interp_property_set_cf(local_code_name,      &
                                            len(local_code_name), &
                                            cpl_id,               &
                                            len(cpl_id),          &
                                            property_name,        &
                                            len(property_name),   &
                                            property_type,        &
                                            property_value,       &
                                            len(property_value))

  end subroutine CWP_Spatial_interp_property_set_


  subroutine CWP_Visu_set_(local_code_name, &
                           cpl_id,          &
                           freq,            &
                           format,          &
                           format_option)
    ! Enable visualization output.
    use, intrinsic :: iso_c_binding
    implicit none

    character(kind = c_char, len = *) :: local_code_name ! Local code name
    character(kind = c_char, len = *) :: cpl_id          ! Coupling identifier
    integer(c_int)                    :: freq            ! Output frequency
    integer(c_int)                    :: format          ! Output format to visualize exchanged fields
    character(kind = c_char, len = *) :: format_option   ! Output options "opt1, opt2, ..."
    integer(kind = c_int) :: l_local_code_name, l_cpl_id, l_format_option

    l_local_code_name = len(local_code_name)
    l_cpl_id = len(cpl_id)
    l_format_option = len(format_option)

    call CWP_Visu_set_cf(local_code_name, l_local_code_name, cpl_id, l_cpl_id, freq, format, format_option, l_format_option)
  end subroutine


  !>
  !!
  !!
  !! \param [in]  local_code_name  Local code name
  !! \param [in]  cpl_id           Coupling identifier
  !! \param [in]  i_part           Current partition
  !! \param [in]  n_pts            Number of points
  !! \param [in]  coord            Coordinates (size = 3 * n_pts)
  !! \param [in]  g_num            global number or NUL (size = n_pts)
  !!

  subroutine CWP_User_tgt_pts_set_(local_code_name, &
                                   cpl_id, &
                                   i_part, &
                                   n_pts, &
                                   coord, &
                                   global_num)
    ! Set a partition of the user target point cloud.
    !
    ! This function must be called if the degrees of freedom locations are
    ! \ref CWP_DOF_LOCATION_USER
    use, intrinsic :: iso_c_binding
    implicit none

    character(kind = c_char, len = *)             :: local_code_name ! Local code name
    character(kind = c_char, len = *)             :: cpl_id          ! Coupling identifier
    integer(kind = c_int)                         :: i_part          ! Current partition
    integer(kind = c_int)                         :: n_pts           ! Number of points
    real(8), dimension(:,:), pointer              :: coord           ! Coordinates (size = 3 * ``n_pts``)
    integer(kind = c_long), dimension(:), pointer :: global_num      ! Global ids (size = ``n_pts`` or  ``null()``)
    integer(kind = c_int) :: l_local_code_name, l_cpl_id
    type(c_ptr) :: c_global_num

    if (associated(global_num)) then
      c_global_num = c_loc(global_num)
    else
      c_global_num = c_null_ptr
    endif

    l_local_code_name = len(local_code_name)
    l_cpl_id = len(cpl_id)

    call CWP_User_tgt_pts_set_cf(local_code_name, l_local_code_name, cpl_id, l_cpl_id, i_part, n_pts, &
            c_loc(coord), c_global_num)
  end subroutine CWP_User_tgt_pts_set_


  subroutine CWP_Mesh_interf_finalize_(local_code_name, &
                                       cpl_id)
    ! Finalize the interface mesh.
    !
    ! This function computes the global ids of mesh entities if they are
    ! not provided.
    use, intrinsic :: iso_c_binding
    implicit none

    character(kind = c_char, len = *) :: local_code_name ! Local code name
    character(kind = c_char, len = *) :: cpl_id          ! Coupling identifier
    integer(kind = c_int) :: l_local_code_name, l_cpl_id

    l_local_code_name = len(local_code_name)
    l_cpl_id = len(cpl_id)

    call CWP_Mesh_interf_finalize_cf(local_code_name, l_local_code_name, cpl_id, l_cpl_id)
  end subroutine CWP_Mesh_interf_finalize_


  subroutine CWP_Mesh_interf_vtx_set_(local_code_name, &
                                      cpl_id,          &
                                      i_part,          &
                                      n_vtx,           &
                                      coord,           &
                                      global_num)
    ! Set the interface mesh vertices.
    use :: iso_c_binding
    implicit none

    character(kind = c_char, len = *)      :: local_code_name ! Local code name
    character(kind = c_char, len = *)      :: cpl_id          ! Coupling identifier
    integer(c_int), intent(in)             :: i_part          ! Current partition
    integer(c_int), intent(in)             :: n_vtx           ! Number of vertices
    real(8), dimension(:,:), pointer       :: coord           ! Coordinates (size = 3 * ``n_vtx``)
    integer(c_long), dimension(:), pointer :: global_num      ! Global vertex ids (size = ``n_vtx`` or ``null()``)
    integer(c_int)                    :: array_size
    integer(kind = c_int) :: l_local_code_name, l_cpl_id

    type(c_ptr) :: c_global_num

    if (associated(global_num)) then
      c_global_num = c_loc(global_num)
    else
      c_global_num = c_null_ptr
    endif

    l_local_code_name = len(local_code_name)
    l_cpl_id = len(cpl_id)

    array_size = size(coord)
    if (modulo(array_size, 3) /= 0) then
        print *, "Error : Length of connectivity array is not a multiple of the 3"
        stop 'error'
    endif


    call CWP_Mesh_interf_vtx_set_cf(local_code_name,   &
                                    l_local_code_name, &
                                    cpl_id,            &
                                    l_cpl_id,          &
                                    i_part,            &
                                    n_vtx,             &
                                    c_loc(coord),      &
                                    c_global_num)

  end subroutine CWP_Mesh_interf_vtx_set_


  function CWP_Mesh_interf_block_add_(local_code_name, &
                                      cpl_id,          &
                                      block_type)      &
  result(block_id)
    ! Add a connectivity block to the interface mesh.
    use, intrinsic :: iso_c_binding
    implicit none

    character(kind = c_char, len = *) :: local_code_name ! Local code name
    character(kind = c_char, len = *) :: cpl_id          ! Coupling identifier
    integer(c_int), intent(in)        :: block_type      ! Block type
    integer(c_int)                    :: block_id        ! Block identifier
    integer(kind = c_int) :: l_local_code_name, l_cpl_id

    l_local_code_name = len(local_code_name)
    l_cpl_id          = len(cpl_id)

    block_id = CWP_Mesh_interf_block_add_cf(local_code_name, l_local_code_name, cpl_id, l_cpl_id, block_type)
  end function CWP_Mesh_interf_block_add_


  subroutine CWP_Mesh_interf_block_std_set_(local_code_name, &
                                            cpl_id,          &
                                            i_part,          &
                                            block_id,        &
                                            n_elts,          &
                                            connec,          &
                                            global_num)
    ! Set a standard block to the interface mesh.
    !
    ! This function adds a connectivity block to the interface mesh.
    ! Refer to CWIPI's :ref:`numbering convention for standard elements <Convention for standard elements>` section.
    !

    use, intrinsic :: iso_c_binding
    implicit none

    character(kind = c_char, len = *)      :: local_code_name ! Local code name
    character(kind = c_char, len = *)      :: cpl_id          ! Coupling identifier
    integer(c_int), intent(in)             :: i_part          ! Partition identifier
    integer(c_int), intent(in)             :: block_id        ! Block identifier
    integer(c_int), intent(in)             :: n_elts          ! Number of elements
    integer(c_int), dimension(:), pointer  :: connec          ! Connectivity (size = *n_vertex_per_elt* * ``n_elts``)
    integer(c_long), dimension(:), pointer :: global_num      ! Global element ids (size = ``n_elts`` or ``null()``)
    ! integer(c_int) :: array_size
    integer(kind = c_int) :: l_local_code_name, l_cpl_id

    type(c_ptr) :: c_global_num

    if (associated(global_num)) then
      c_global_num = c_loc(global_num)
    else
      c_global_num = c_null_ptr
    endif

    l_local_code_name = len(local_code_name)
    l_cpl_id = len(cpl_id)

    ! TODO get stride
    ! array_size = size(connec)
    ! if (modulo(array_size, 3) /= 0) then
    !     print *, "Error : Length of connectivity array is not a multiple of the element size"
    !     stop 'error'
    ! endif

    call CWP_Mesh_interf_block_std_set_cf (local_code_name,   &
                                           l_local_code_name, &
                                           cpl_id,            &
                                           l_cpl_id,          &
                                           i_part,            &
                                           block_id,          &
                                           n_elts,            &
                                           c_loc(connec),     &
                                           c_global_num)
  end subroutine CWP_Mesh_interf_block_std_set_


  subroutine CWP_Mesh_interf_block_std_get_(local_code_name, &
                                            cpl_id,          &
                                            i_part,          &
                                            block_id,        &
                                            n_elts,          &
                                            connec,          &
                                            global_num)
    ! Get the properties of a standard block of the interface mesh.
    use, intrinsic :: iso_c_binding
    implicit none

    character(kind = c_char, len = *)      :: local_code_name ! Local code name
    character(kind = c_char, len = *)      :: cpl_id          ! Coupling identifier
    integer(c_int), intent(in)             :: i_part          ! Partition identifier
    integer(c_int), intent(in)             :: block_id        ! Block identifier
    integer(c_int), intent(out)            :: n_elts          ! Number of elements
    integer(c_int),  dimension(:), pointer :: connec          ! Connectivity (size = `n_vertex_per_elt` * ``n_elts``)
    integer(c_long), dimension(:), pointer :: global_num      ! Global element ids (size = ``n_elts``)
    integer(kind = c_int) :: l_local_code_name, l_cpl_id, s_connec
    type(c_ptr)      :: c_connec, c_global_num

    l_local_code_name = len(local_code_name)
    l_cpl_id          = len(cpl_id)

    call CWP_Mesh_interf_block_std_get_cf(local_code_name,   &
                                          l_local_code_name, &
                                          cpl_id,            &
                                          l_cpl_id,          &
                                          i_part,            &
                                          block_id,          &
                                          n_elts,            &
                                          c_connec,          &
                                          c_global_num,      &
                                          s_connec)

    call c_f_pointer(c_global_num, global_num, [n_elts])
    call c_f_pointer(c_connec,     connec,     [s_connec])

  end subroutine CWP_Mesh_interf_block_std_get_


  subroutine CWP_Mesh_interf_f_poly_block_set_(local_code_name, &
                                               cpl_id,          &
                                               i_part,          &
                                               block_id,        &
                                               n_elts,          &
                                               connec_idx,      &
                                               connec,          &
                                               global_num)
    ! Set the connectivity of a polygon block in an interface mesh partition.
    use, intrinsic :: iso_c_binding
    implicit none

    character(kind = c_char, len = *)      :: local_code_name ! Local code name
    character(kind = c_char, len = *)      :: cpl_id          ! Coupling identifier
    integer(c_int), intent(in)             :: i_part          ! Current partition
    integer(c_int), intent(in)             :: block_id        ! Block identifier
    integer(c_int), intent(in)             :: n_elts          ! Number of elements
    integer(c_int),  dimension(:), pointer :: connec_idx      ! Connectivity index (``connec_idx(0)`` = 0 and size = ``n_elts`` + 1)
    integer(c_int),  dimension(:), pointer :: connec          ! Connectivity (size = ``connec_idx(n_elts+1)``)
    integer(c_long), dimension(:), pointer :: global_num      ! Global element ids (size = ``n_elts`` or ``null()``)
    integer(kind = c_int) :: l_local_code_name, l_cpl_id


    type(c_ptr) :: c_global_num

    if (associated(global_num)) then
      c_global_num = c_loc(global_num)
    else
      c_global_num = c_null_ptr
    endif

    l_local_code_name = len(local_code_name)
    l_cpl_id = len(cpl_id)

    call CWP_Mesh_interf_f_poly_block_set_cf (local_code_name,   &
                                              l_local_code_name, &
                                              cpl_id,            &
                                              l_cpl_id,          &
                                              i_part,            &
                                              block_id,          &
                                              n_elts,            &
                                              c_loc(connec_idx), &
                                              c_loc(connec),     &
                                              c_global_num)
  end subroutine CWP_Mesh_interf_f_poly_block_set_


  subroutine CWP_Mesh_interf_f_poly_block_get_(local_code_name, &
                                               cpl_id,          &
                                               i_part,          &
                                               block_id,        &
                                               n_elts,          &
                                               connec_idx,      &
                                               connec,          &
                                               global_num)
    ! Get the properties of a polygon block of the interface mesh partition.
    use, intrinsic :: iso_c_binding
    implicit none

    character(kind = c_char, len = *)      :: local_code_name ! Local code name
    character(kind = c_char, len = *)      :: cpl_id          ! Coupling identifier
    integer(c_int), intent(in)             :: i_part          ! Partition identifier
    integer(c_int), intent(in)             :: block_id        ! Block identifier
    integer(c_int), intent(out)            :: n_elts          ! Number of elements
    integer(c_int),  dimension(:), pointer :: connec_idx      ! Connectivity index (``connec_idx(0)`` = 0 and size = ``n_elts`` + 1)
    integer(c_int),  dimension(:), pointer :: connec          ! Connectivity (size = ``connec_idx(n_elts+1)``)
    integer(c_long), dimension(:), pointer :: global_num      ! Global element ids (size = ``n_elts``)

    integer(kind = c_int) :: l_local_code_name, l_cpl_id
    type(c_ptr) :: c_connec_idx, c_connec, c_global_num

    l_local_code_name = len(local_code_name)
    l_cpl_id = len(cpl_id)

    call CWP_Mesh_interf_f_poly_block_get_cf(local_code_name,   &
                                             l_local_code_name, &
                                             cpl_id,            &
                                             l_cpl_id,          &
                                             i_part,            &
                                             block_id,          &
                                             n_elts,            &
                                             c_connec_idx,      &
                                             c_connec,          &
                                             c_global_num)

    call c_f_pointer(c_connec_idx, connec_idx, [n_elts+1])
    call c_f_pointer(c_connec,     connec,     [connec_idx(n_elts+1)])
    call c_f_pointer(c_global_num, global_num, [n_elts])

  end subroutine CWP_Mesh_interf_f_poly_block_get_


  subroutine CWP_Mesh_interf_c_poly_block_set_(local_code_name,  &
                                               cpl_id,           &
                                               i_part,           &
                                               block_id,         &
                                               n_cell,           &
                                               n_face,           &
                                               connec_faces_idx, &
                                               connec_faces,     &
                                               connec_cells_idx, &
                                               connec_cells,     &
                                               global_num)
    ! Set the properties of a polyhedron block of the interface mesh partition.
    use, intrinsic :: iso_c_binding
    implicit none

    character(kind = c_char, len = *)      :: local_code_name  ! Local code name
    character(kind = c_char, len = *)      :: cpl_id           ! Coupling identifier
    integer(c_int), intent(in)             :: i_part           ! Partition identifier
    integer(c_int), intent(in)             :: block_id         ! Block identifier
    integer(c_int), intent(in)             :: n_cell           ! Number of polyhedra
    integer(c_int), intent(in)             :: n_face           ! Number of faces
    integer(c_int),  dimension(:), pointer :: connec_faces_idx ! Index for face to vertex connectivity (``connec_faces_idx(0)`` = 0 and size = ``n_face`` + 1)
    integer(c_int),  dimension(:), pointer :: connec_faces     ! Face to vertex connectivity (size = ``connec_faces_idx(n_face+1)``)
    integer(c_int),  dimension(:), pointer :: connec_cells_idx ! Index for polyhedron to face connectivity (``connec_cells_idx(0)`` = 0 and size = ``n_cell`` + 1)
    integer(c_int),  dimension(:), pointer :: connec_cells     ! Polyhedron to face connectivity (size = ``connec_cells_idx(n_cell+1)``)
    integer(c_long), dimension(:), pointer :: global_num       ! Global cell ids (size = ``n_cell`` or ``null()``)
    integer(kind = c_int) :: l_local_code_name, l_cpl_id

    type(c_ptr) :: c_global_num

    if (associated(global_num)) then
      c_global_num = c_loc(global_num)
    else
      c_global_num = c_null_ptr
    endif

    l_local_code_name = len(local_code_name)
    l_cpl_id = len(cpl_id)

    call CWP_Mesh_interf_c_poly_block_set_cf (local_code_name,         &
                                              l_local_code_name,       &
                                              cpl_id,                  &
                                              l_cpl_id,                &
                                              i_part,                  &
                                              block_id,                &
                                              n_cell,                  &
                                              n_face,                  &
                                              c_loc(connec_faces_idx), &
                                              c_loc(connec_faces),     &
                                              c_loc(connec_cells_idx), &
                                              c_loc(connec_cells),     &
                                              c_global_num)
  end subroutine CWP_Mesh_interf_c_poly_block_set_


  subroutine CWP_Mesh_interf_c_poly_block_get_(local_code_name, &
                                              cpl_id,           &
                                              i_part,           &
                                              block_id,         &
                                              n_cell,           &
                                              n_face,           &
                                              connec_faces_idx, &
                                              connec_faces,     &
                                              connec_cells_idx, &
                                              connec_cells,     &
                                              global_num)
    ! Get the properties of a polyhedron block of the interface mesh partition.
    use, intrinsic :: iso_c_binding
    implicit none

    character(kind = c_char, len = *)      :: local_code_name  ! Local code name
    character(kind = c_char, len = *)      :: cpl_id           ! Coupling identifier
    integer(c_int), intent(in)             :: i_part           ! Partition identifier
    integer(c_int), intent(in)             :: block_id         ! Block identifier
    integer(c_int), intent(out)            :: n_cell           ! Number of polyhedra
    integer(c_int), intent(out)            :: n_face           ! Number of faces
    integer(c_int),  dimension(:), pointer :: connec_faces_idx ! Index for face to vertex connectivity (``connec_faces_idx(0)`` = 0 and size = ``n_face`` + 1)
    integer(c_int),  dimension(:), pointer :: connec_faces     ! Face to vertex connectivity (size = ``connec_faces_idx(n_face+1)``)
    integer(c_int),  dimension(:), pointer :: connec_cells_idx ! Index for polyhedron to face connectivity (``connec_cells_idx(0)`` = 0 and size = ``n_cell`` + 1)
    integer(c_int),  dimension(:), pointer :: connec_cells     ! Polyhedron to face connectivity (size = ``connec_cells_idx(n_cell+1)``)
    integer(c_long), dimension(:), pointer :: global_num       ! Global cell ids (size = ``n_cell`` or ``null()``)
    integer(kind = c_int) :: l_local_code_name, l_cpl_id
    type(c_ptr) :: c_connec_faces_idx, c_connec_faces, c_connec_cells_idx, c_connec_cells, c_global_num

    l_local_code_name = len(local_code_name)
    l_cpl_id = len(cpl_id)

    call CWP_Mesh_interf_c_poly_block_get_cf(local_code_name,    &
                                             l_local_code_name,  &
                                             cpl_id,             &
                                             l_cpl_id,           &
                                             i_part,             &
                                             block_id,           &
                                             n_cell,             &
                                             n_face,             &
                                             c_connec_faces_idx, &
                                             c_connec_faces,     &
                                             c_connec_cells_idx, &
                                             c_connec_cells,     &
                                             c_global_num)

    call c_f_pointer(c_connec_faces_idx, connec_faces_idx, [n_face+1])
    call c_f_pointer(c_connec_faces,     connec_faces,     [connec_faces_idx(n_face+1)])
    call c_f_pointer(c_connec_cells_idx, connec_cells_idx, [n_cell+1])
    call c_f_pointer(c_connec_cells,     connec_cells,     [connec_cells_idx(n_cell+1)])
    call c_f_pointer(c_global_num, global_num, [n_cell])

  end subroutine CWP_Mesh_interf_c_poly_block_get_


  subroutine CWP_Mesh_interf_del_ (local_code_name, &
                                  cpl_id)
    ! Delete the interface mesh.
    use, intrinsic :: iso_c_binding
    implicit none

    character(kind = c_char, len = *) :: local_code_name ! Local code name
    character(kind = c_char, len = *) :: cpl_id          ! Coupling identifier
    integer(c_int) :: l_local_code_name, l_cpl_id

    l_local_code_name = len(local_code_name)
    l_cpl_id          = len(cpl_id)

    call CWP_Mesh_interf_del_cf (local_code_name, l_local_code_name, cpl_id, l_cpl_id)
  end subroutine CWP_Mesh_interf_del_


  subroutine CWP_Mesh_interf_from_cellface_set_(local_code_name, &
                                                cpl_id,          &
                                                i_part,          &
                                                n_cells,         &
                                                cell_face_idx,   &
                                                cell_face,       &
                                                n_faces,         &
                                                face_vtx_idx,    &
                                                face_vtx,        &
                                                global_num)
    ! Define the interface mesh from a cell-to-face connectivity.
    use, intrinsic :: iso_c_binding
    implicit none

    character(kind = c_char, len = *)      :: local_code_name ! Local code name
    character(kind = c_char, len = *)      :: cpl_id          ! Coupling identifier
    integer(c_int)                         :: i_part          ! Partition identifier
    integer(c_int)                         :: n_cells         ! Number of cells
    integer(c_int),  dimension(:), pointer :: cell_face_idx   ! Index for cell to face connectivity (``cell_face_idx(0)`` = 0 and size = ``n_cells`` + 1
    integer(c_int),  dimension(:), pointer :: cell_face       ! Cell to face connectivity (size = ``cell_face_idx(n_cells+1)``)
    integer(c_int)                         :: n_faces         ! Number of faces
    integer(c_int),  dimension(:), pointer :: face_vtx_idx    ! Index for face to vertex connectivity (``face_vtx_idx(0)`` = 0 and size = ``n_faces`` + 1)
    integer(c_int),  dimension(:), pointer :: face_vtx        ! Face to vertex connectivity (size = ``face_vtx_idx(n_faces+1)``)
    integer(c_long), dimension(:), pointer :: global_num      ! Global cell ids (size = ``n_cells`` or ``null()``)

    integer(kind = c_int) :: l_local_code_name, l_cpl_id
    type(c_ptr) :: c_global_num

    if (associated(global_num)) then
      c_global_num = c_loc(global_num)
    else
      c_global_num = c_null_ptr
    endif

    l_local_code_name = len(local_code_name)
    l_cpl_id = len(cpl_id)


    call CWP_Mesh_interf_from_cellface_set_cf (local_code_name,      &
                                               l_local_code_name,    &
                                               cpl_id,               &
                                               l_cpl_id,             &
                                               i_part,               &
                                               n_cells,              &
                                               c_loc(cell_face_idx), &
                                               c_loc(cell_face),     &
                                               n_faces,              &
                                               c_loc(face_vtx_idx),  &
                                               c_loc(face_vtx),      &
                                               c_global_num)
  end subroutine CWP_Mesh_interf_from_cellface_set_


  subroutine CWP_Mesh_interf_from_faceedge_set_(local_code_name, &
                                                cpl_id,          &
                                                i_part,          &
                                                n_faces,         &
                                                face_edge_idx,   &
                                                face_edge,       &
                                                n_edges,         &
                                                edge_vtx,        &
                                                global_num)
    ! Define the surface interface mesh from a face-to-edge connectivity.
    use, intrinsic :: iso_c_binding
    implicit none

    character(kind = c_char, len = *)      :: local_code_name ! Local code name
    character(kind = c_char, len = *)      :: cpl_id          ! Coupling identifier
    integer(c_int)                         :: i_part          ! Partition identifier
    integer(c_int)                         :: n_faces         ! Number of faces
    integer(c_int),  dimension(:), pointer :: face_edge_idx   ! Index for face to edge connectivity (``face_edge_idx(0)`` = 0 and size = ``n_faces`` + 1
    integer(c_int),  dimension(:), pointer :: face_edge       ! Face to edge connectivity (size = ``face_edge_idx(n_faces+1)``)
    integer(c_int)                         :: n_edges         ! Number of edges
    integer(c_int),  dimension(:), pointer :: edge_vtx        ! Edge to vertex connectivity (size = 2*``n_edges``)
    integer(c_long), dimension(:), pointer :: global_num      ! Global face ids (size = ``n_faces`` or ``null()``)

    integer(kind = c_int) :: l_local_code_name, l_cpl_id
    type(c_ptr) :: c_global_num

    if (associated(global_num)) then
      c_global_num = c_loc(global_num)
    else
      c_global_num = c_null_ptr
    endif

    l_local_code_name = len(local_code_name)
    l_cpl_id = len(cpl_id)

    call CWP_Mesh_interf_from_faceedge_set_cf (local_code_name,     &
                                               l_local_code_name,   &
                                               cpl_id,              &
                                               l_cpl_id,            &
                                               i_part,              &
                                               n_faces,             &
                                               c_loc(face_edge_idx),&
                                               c_loc(face_edge),    &
                                               n_edges,             &
                                               c_loc(edge_vtx),     &
                                               c_global_num)

  end subroutine CWP_Mesh_interf_from_faceedge_set_


  subroutine CWP_Mesh_interf_from_facevtx_set_(local_code_name, &
                                               cpl_id,          &
                                               i_part,          &
                                               n_faces,         &
                                               face_vtx_idx,    &
                                               face_vtx,        &
                                               global_num)
    ! Define the surface interface mesh from a face-to-vertex connectivity.
    use, intrinsic :: iso_c_binding
    implicit none

    character(kind = c_char, len = *)      :: local_code_name ! Local code name
    character(kind = c_char, len = *)      :: cpl_id          ! Coupling identifier
    integer(c_int)                         :: i_part          ! Partition identifier
    integer(c_int)                         :: n_faces         ! Number of faces
    integer(c_int),  dimension(:), pointer :: face_vtx_idx    ! Index for face to vertex connectivity (``face_vtx_idx(0)`` = 0 and size = ``n_faces`` + 1
    integer(c_int),  dimension(:), pointer :: face_vtx        ! Face to vertex connectivity (size = ``face_vtx_idx(n_faces+1)``)
    integer(c_long), dimension(:), pointer :: global_num      ! Global face ids (size = ``n_faces`` or ``null()``)

    integer(kind = c_int) :: l_local_code_name, l_cpl_id
    type(c_ptr) :: c_global_num

    if (associated(global_num)) then
      c_global_num = c_loc(global_num)
    else
      c_global_num = c_null_ptr
    endif

    l_local_code_name = len(local_code_name)
    l_cpl_id = len(cpl_id)

    call CWP_Mesh_interf_from_facevtx_set_cf(local_code_name,     &
                                             l_local_code_name,   &
                                             cpl_id,              &
                                             l_cpl_id,            &
                                             i_part,              &
                                             n_faces,             &
                                             c_loc(face_vtx_idx), &
                                             c_loc(face_vtx),     &
                                             c_global_num)

  end subroutine CWP_Mesh_interf_from_facevtx_set_


  subroutine CWP_Field_create_(local_code_name,      &
                               cpl_id,               &
                               field_id,             &
                               data_type,            &
                               storage,              &
                               n_component,          &
                               target_location,      &
                               exch_type,            &
                               visu_status)
    ! Create a new field.
    use, intrinsic :: iso_c_binding
    implicit none

    character(kind = c_char, len = *) :: local_code_name ! Local code name
    character(kind = c_char, len = *) :: cpl_id          ! Coupling identifier
    character(kind = c_char, len = *) :: field_id        ! Field id
    integer(c_int)                    :: data_type       ! Data type
    integer(c_int)                    :: storage         ! Storage type
    integer(c_int)                    :: n_component     ! Number of component
    integer(c_int)                    :: target_location ! Target location
    integer(c_int)                    :: exch_type       ! Exchange type
    integer(c_int)                    :: visu_status     ! Visualization status
    integer(kind = c_int) :: l_local_code_name, l_cpl_id, l_field_id

    l_local_code_name = len(local_code_name)
    l_cpl_id = len(cpl_id)
    l_field_id = len(field_id)

    call CWP_Field_create_cf (local_code_name,    &
                              l_local_code_name,  &
                              cpl_id,             &
                              l_cpl_id,           &
                              field_id,           &
                              l_field_id,         &
                              data_type,          &
                              storage,            &
                              n_component,        &
                              target_location,    &
                              exch_type,          &
                              visu_status)
  end subroutine CWP_Field_create_


  subroutine CWP_Field_data_set_(local_code_name, &
                                 cpl_id,          &
                                 field_id,        &
                                 i_part,          &
                                 map_type,        &
                                 data)
    ! Set field data.
    use, intrinsic :: iso_c_binding
    implicit none

    character(kind = c_char, len = *) :: local_code_name ! Local code name
    character(kind = c_char, len = *) :: cpl_id          ! Coupling identifier
    character(kind = c_char, len = *) :: field_id        ! Field id
    integer(c_int)                    :: i_part          ! Current partition
    integer(c_int)                    :: map_type        ! Choice if data is set for the source or the target
    real(8), dimension(:), pointer    :: data            ! Storage array (Mapping)
    integer(c_int) :: n_dof, n_components, array_size, spatial_interp_algorithm
    integer(kind = c_int) :: l_local_code_name, l_cpl_id, l_field_id

    l_local_code_name = len(local_code_name)
    l_cpl_id = len(cpl_id)
    l_field_id = len(field_id)

    ! check data array size
    array_size = size(data)
    n_dof = CWP_Field_n_dof_get(local_code_name,   &
                                cpl_id,            &
                                field_id,          &
                                i_part)
    n_components = CWP_Field_n_components_get(local_code_name,   &
                                              cpl_id,            &
                                              field_id)
    spatial_interp_algorithm = CWP_Cpl_spatial_interp_algo_get(local_code_name,   &
                                                               cpl_id)

    if (spatial_interp_algorithm .eq. CWP_SPATIAL_INTERP_FROM_NEAREST_SOURCES_LEAST_SQUARES .or. &
        spatial_interp_algorithm .eq. CWP_SPATIAL_INTERP_FROM_NEAREST_TARGETS_LEAST_SQUARES .or. &
        spatial_interp_algorithm .eq. CWP_SPATIAL_INTERP_FROM_INTERSECTION .or.                  &
        spatial_interp_algorithm .eq. CWP_SPATIAL_INTERP_FROM_IDENTITY) then

        if (map_type .eq. CWP_FIELD_MAP_TARGET) then

            if (array_size < n_dof * n_components) then
                print *, "Error : Data array size is inconsistent (too small) with the number of components &
                          & and degree of freedom location"
                stop
            endif

        endif

    else if (spatial_interp_algorithm .eq. CWP_SPATIAL_INTERP_FROM_LOCATION_MESH_LOCATION_LOCATE_ALL_TGT .or. &
             spatial_interp_algorithm .eq. CWP_SPATIAL_INTERP_FROM_LOCATION_MESH_LOCATION_OCTREE .or.         &
             spatial_interp_algorithm .eq. CWP_SPATIAL_INTERP_FROM_LOCATION_MESH_LOCATION_BOXTREE) then

        if (map_type .eq. CWP_FIELD_MAP_SOURCE) then

            if (array_size < n_dof * n_components) then
                print *, "Error : Data array size is inconsistent (too small) with the number of components &
                          & and degree of freedom location"
                stop
            endif

        endif

    endif

    call CWP_Field_data_set_cf (local_code_name,   &
                                l_local_code_name, &
                                cpl_id,            &
                                l_cpl_id,          &
                                field_id,          &
                                l_field_id,        &
                                i_part,            &
                                map_type,          &
                                c_loc(data))
  end subroutine CWP_Field_data_set_


  function CWP_Field_dof_location_get_(local_code_name, &
                                       cpl_id,          &
                                       field_id)        &
    result (dof_location)
    ! Get target degrees of freedom location.
    use, intrinsic :: iso_c_binding
    implicit none

    character(kind = c_char, len = *) :: local_code_name ! Local code name
    character(kind = c_char, len = *) :: cpl_id          ! Coupling identifier
    character(kind = c_char, len = *) :: field_id        ! Field id
    integer(c_int)                    :: dof_location    ! Location of degrees of freedom
    integer(kind = c_int) :: l_local_code_name, l_cpl_id, l_field_id

    l_local_code_name = len(local_code_name)
    l_cpl_id          = len(cpl_id)
    l_field_id        = len(field_id)

    dof_location = CWP_Field_dof_location_get_cf(local_code_name,   &
                                                 l_local_code_name, &
                                                 cpl_id,            &
                                                 l_cpl_id,          &
                                                 field_id,          &
                                                 l_field_id)

  end function CWP_Field_dof_location_get_


  function CWP_Field_storage_get_(local_code_name, &
                                  cpl_id,          &
                                  field_id)        &
    result (storage)
    ! Get field storage type.
    use, intrinsic :: iso_c_binding
    implicit none

    character(kind = c_char, len = *) :: local_code_name ! Local code name
    character(kind = c_char, len = *) :: cpl_id          ! Coupling identifier
    character(kind = c_char, len = *) :: field_id        ! Field id
    integer(c_int)                    :: storage         ! Field storage type
    integer(kind = c_int) :: l_local_code_name, l_cpl_id, l_field_id

    l_local_code_name = len(local_code_name)
    l_cpl_id          = len(cpl_id)
    l_field_id        = len(field_id)

    storage = CWP_Field_storage_get_cf(local_code_name,   &
                                       l_local_code_name, &
                                       cpl_id,            &
                                       l_cpl_id,          &
                                       field_id,          &
                                       l_field_id)

  end function CWP_Field_storage_get_


  function CWP_Field_n_dof_get_(local_code_name, &
                                cpl_id,          &
                                field_id,        &
                                i_part)          &
    result (n_dof)
    ! Get field number of degrees of freedom.
    use, intrinsic :: iso_c_binding
    implicit none

    character(kind = c_char, len = *) :: local_code_name ! Local code name
    character(kind = c_char, len = *) :: cpl_id          ! Coupling identifier
    character(kind = c_char, len = *) :: field_id        ! Field id
    integer(c_int)                    :: i_part          ! Current partition
    integer(c_int)                    :: n_dof           ! Field number of degrees of freedom
    integer(kind = c_int) :: l_local_code_name, l_cpl_id, l_field_id

    l_local_code_name = len(local_code_name)
    l_cpl_id          = len(cpl_id)
    l_field_id        = len(field_id)

    n_dof = CWP_Field_n_dof_get_cf(local_code_name,   &
                                   l_local_code_name, &
                                   cpl_id,            &
                                   l_cpl_id,          &
                                   field_id,          &
                                   l_field_id,        &
                                   i_part)

  end function CWP_Field_n_dof_get_


  subroutine CWP_Field_del_(local_code_name, &
                            cpl_id,          &
                            field_id)
    ! Delete a field.
    use, intrinsic :: iso_c_binding
    implicit none

    character(kind = c_char, len = *) :: local_code_name ! Local code name
    character(kind = c_char, len = *) :: cpl_id          ! Coupling identifier
    character(kind = c_char, len = *) :: field_id        ! Field id
    integer(kind = c_int) :: l_local_code_name, l_cpl_id, l_field_id

    l_local_code_name = len(local_code_name)
    l_cpl_id          = len(cpl_id)
    l_field_id        = len(field_id)

    call CWP_Field_del_cf(local_code_name,   &
                          l_local_code_name, &
                          cpl_id,            &
                          l_cpl_id,          &
                          field_id,          &
                          l_field_id)

  end subroutine CWP_Field_del_


  subroutine CWP_Field_issend_(local_code_name, &
                               cpl_id,          &
                               field_id)
    ! Send a spatially interpolated field to the coupled code with
    ! non-blocking communications.
    !
    ! This function is independent of \ref CWP_Time_exch_t mode. The user has to
    ! manually check the consistency of the exchanges.

    use, intrinsic :: iso_c_binding
    implicit none

    character(kind = c_char, len = *) :: local_code_name ! Local code name
    character(kind = c_char, len = *) :: cpl_id          ! Coupling identifier
    character(kind = c_char, len = *) :: field_id        ! Field id
    integer(kind = c_int) :: l_local_code_name, l_cpl_id, l_field_id

    l_local_code_name = len(local_code_name)
    l_cpl_id = len(cpl_id)
    l_field_id = len(field_id)

    call CWP_Field_issend_cf (local_code_name,    &
                              l_local_code_name,  &
                              cpl_id,             &
                              l_cpl_id,           &
                              field_id,       &
                              l_field_id)
  end subroutine CWP_Field_issend_


  subroutine CWP_Field_irecv_(local_code_name, &
                              cpl_id,          &
                              field_id)
    ! Receive a spatially interpolated field from the coupled code
    ! with non-blocking communications.
    !
    ! This function is independent of \ref CWP_Time_exch_t mode. The user has to
    ! manually check the consistency of the exchanges.

    use, intrinsic :: iso_c_binding
    implicit none

    character(kind = c_char, len = *) :: local_code_name ! Local code name
    character(kind = c_char, len = *) :: cpl_id          ! Coupling identifier
    character(kind = c_char, len = *) :: field_id        ! Field id
    integer(kind = c_int) :: l_local_code_name, l_cpl_id, l_field_id

    l_local_code_name = len(local_code_name)
    l_cpl_id          = len(cpl_id)
    l_field_id        = len(field_id)

    call CWP_Field_irecv_cf(local_code_name,   &
                            l_local_code_name, &
                            cpl_id,            &
                            l_cpl_id,          &
                            field_id,          &
                            l_field_id)
  end subroutine CWP_Field_irecv_


  subroutine CWP_Field_wait_issend_(local_code_name, &
                                    cpl_id,          &
                                    field_id)
    ! Wait the end of an exchange related to request from \ref CWP_Field_issend.
    use, intrinsic :: iso_c_binding
    implicit none

    character(kind = c_char, len = *) :: local_code_name ! Local code name
    character(kind = c_char, len = *) :: cpl_id          ! Coupling identifier
    character(kind = c_char, len = *) :: field_id        ! Field id
    integer(kind = c_int) :: l_local_code_name, l_cpl_id, l_field_id

    l_local_code_name = len(local_code_name)
    l_cpl_id          = len(cpl_id)
    l_field_id        = len(field_id)

    call CWP_Field_wait_issend_cf (local_code_name,    &
                                   l_local_code_name,  &
                                   cpl_id,             &
                                   l_cpl_id,           &
                                   field_id,       &
                                   l_field_id)
  end subroutine CWP_Field_wait_issend_


  subroutine CWP_Field_wait_irecv_(local_code_name, &
                                   cpl_id,          &
                                   field_id)
    ! Wait the end of an exchange related to request from \ref CWP_Field_irecv.

    use, intrinsic :: iso_c_binding
    implicit none

    character(kind = c_char, len = *) :: local_code_name ! Local code name
    character(kind = c_char, len = *) :: cpl_id          ! Coupling identifier
    character(kind = c_char, len = *) :: field_id        ! Field id
    integer(kind = c_int) :: l_local_code_name, l_cpl_id, l_field_id

    l_local_code_name = len(local_code_name)
    l_cpl_id          = len(cpl_id)
    l_field_id        = len(field_id)

    call CWP_Field_wait_irecv_cf(local_code_name,     &
                                 l_local_code_name,   &
                                 cpl_id,              &
                                 l_cpl_id,            &
                                 field_id,            &
                                 l_field_id)
  end subroutine CWP_Field_wait_irecv_


  subroutine CWP_Field_interp_function_unset_(local_code_name, &
                                              cpl_id,          &
                                              field_id)
    ! Unset a user-defined spatial interpolation function.
    use, intrinsic :: iso_c_binding
    implicit none

    character(kind = c_char, len = *) :: local_code_name ! Local code name
    character(kind = c_char, len = *) :: cpl_id          ! Coupling identifier
    character(kind = c_char, len = *) :: field_id        ! Field id
    integer(kind = c_int) :: l_local_code_name, l_cpl_id, l_field_id
    l_local_code_name = len(local_code_name)
    l_cpl_id = len(cpl_id)
    l_field_id = len(field_id)

    call CWP_Field_interp_function_unset_cf(local_code_name, l_local_code_name, cpl_id, l_cpl_id, &
            & field_id, l_field_id)

  end subroutine CWP_Field_interp_function_unset_


  subroutine CWP_Field_interp_function_set_(local_code_name,        &
                                            cpl_id,                 &
                                            field_id,               &
                                            user_interpolation_fct)
    ! Set a user-defined spatial interpolation function.
    use, intrinsic :: iso_c_binding
    implicit none

    interface
        subroutine user_interpolation_fct(c_local_code_name, &
                                          c_cpl_id,          &
                                          c_field_id,        &
                                          i_part,            &
                                          c_buffer_in,       &
                                          c_buffer_out)      &
          bind(c)
          ! User-defined spatial interpolation function
          use, intrinsic :: iso_c_binding
          implicit none

          character(kind=c_char,len=1) :: c_local_code_name(*)
          character(kind=c_char,len=1) :: c_cpl_id(*)
          character(kind=c_char,len=1) :: c_field_id(*)
          integer(kind=c_int), value   :: i_part
          type(c_ptr), value           :: c_buffer_in
          type(c_ptr), value           :: c_buffer_out
        end subroutine user_interpolation_fct
    end interface

    character(kind = c_char, len = *) :: local_code_name ! Local code name
    character(kind = c_char, len = *) :: cpl_id          ! Coupling identifier
    character(kind = c_char, len = *) :: field_id        ! Field id
    integer(kind = c_int)             :: l_local_code_name, l_cpl_id, l_field_id
    l_local_code_name = len(local_code_name)
    l_cpl_id          = len(cpl_id)
    l_field_id        = len(field_id)

    call CWP_Field_interp_function_set_cf(local_code_name,   &
                                    l_local_code_name, &
                                    cpl_id,            &
                                    l_cpl_id,          &
                                    field_id,      &
                                    l_field_id,    &
                                    c_funloc(user_interpolation_fct))
  end subroutine CWP_Field_interp_function_set_


  function CWP_Field_n_components_get_(local_code_name, &
                                       cpl_id,          &
                                       field_id)        &
    ! Get number of field components.
    result (n_components)
    use, intrinsic :: iso_c_binding
    implicit none

    character(kind = c_char, len = *) :: local_code_name ! Local code name
    character(kind = c_char, len = *) :: cpl_id          ! Coupling identifier
    character(kind = c_char, len = *) :: field_id        ! Field id
    integer                           :: n_components    ! Number of field components
    integer(kind = c_int) :: l_local_code_name, l_cpl_id, l_field_id
    l_local_code_name = len(local_code_name)
    l_cpl_id          = len(cpl_id)
    l_field_id        = len(field_id)

    n_components = CWP_Field_n_components_get_cf(local_code_name, l_local_code_name, cpl_id, l_cpl_id, &
                                                 & field_id, l_field_id)

  end function CWP_Field_n_components_get_


  subroutine CWP_Field_src_data_properties_get_(local_code_name, &
                                                cpl_id,          &
                                                field_id,        &
                                                i_part,          &
                                                n_src,           &
                                                src_to_tgt_idx)
    ! Get spatial interpolation source data.
    use, intrinsic :: iso_c_binding
    implicit none

    character(kind = c_char, len = *)      :: local_code_name ! Local code name
    character(kind = c_char, len = *)      :: cpl_id          ! Coupling identifier
    character(kind = c_char, len = *)      :: field_id        ! Field id
    integer(c_int), intent(in)             :: i_part          ! Partition identifier
    integer(c_int), intent(out)            :: n_src           ! Number of source dofs
    integer(c_int),  dimension(:), pointer :: src_to_tgt_idx  ! Index for source->target mapping (size = ``n_src`` + 1)
    integer(kind = c_int) :: l_local_code_name, l_cpl_id, l_field_id
    type(c_ptr)      :: c_src_to_tgt_idx

    l_local_code_name = len(local_code_name)
    l_cpl_id = len(cpl_id)
    l_field_id = len(field_id)

    call CWP_Field_src_data_properties_get_cf(local_code_name, &
                                              l_local_code_name, &
                                              cpl_id, &
                                              l_cpl_id, &
                                              field_id, &
                                              l_field_id, &
                                              i_part, &
                                              n_src, &
                                              c_src_to_tgt_idx)

    call c_f_pointer(c_src_to_tgt_idx, src_to_tgt_idx, [n_src+1])

  end subroutine CWP_Field_src_data_properties_get_


  subroutine CWP_Field_tgt_data_properties_get_(local_code_name, &
                                                cpl_id, &
                                                field_id, &
                                                i_part, &
                                                n_tgt, &
                                                n_computed_tgt, &
                                                computed_tgt, &
                                                tgt_to_src_idx)

    use, intrinsic :: iso_c_binding
    implicit none

    character(kind = c_char, len = *)      :: local_code_name ! Local code name
    character(kind = c_char, len = *)      :: cpl_id          ! Coupling identifier
    character(kind = c_char, len = *)      :: field_id        ! Field id
    integer(c_int), intent(in)             :: i_part          ! Partition identifier
    integer(c_int), intent(out)            :: n_tgt           ! Number of target dofs
    integer(c_int), intent(out)            :: n_computed_tgt  ! Number of computed target dofs
    integer(c_int),  dimension(:), pointer :: computed_tgt    ! Computed target dofs (size = ``n_computed_tgt``)
    integer(c_int),  dimension(:), pointer :: tgt_to_src_idx  ! Index for target->source mapping (size = ``n_computed_tgt`` + 1)
    integer(kind = c_int) :: l_local_code_name, l_cpl_id, l_field_id
    type(c_ptr)      :: c_computed_tgt, c_tgt_to_src_idx

    l_local_code_name = len(local_code_name)
    l_cpl_id = len(cpl_id)
    l_field_id = len(field_id)

    call CWP_Field_tgt_data_properties_get_cf(local_code_name, &
                                              l_local_code_name, &
                                              cpl_id, &
                                              l_cpl_id, &
                                              field_id, &
                                              l_field_id, &
                                              i_part, &
                                              n_tgt, &
                                              n_computed_tgt, &
                                              c_computed_tgt, &
                                              c_tgt_to_src_idx)

    call c_f_pointer(c_computed_tgt, computed_tgt, [n_computed_tgt])
    call c_f_pointer(c_tgt_to_src_idx, tgt_to_src_idx, [n_computed_tgt+1])

  end subroutine CWP_Field_tgt_data_properties_get_


  subroutine CWP_Field_location_weights_get_(local_code_name, &
                                             cpl_id,          &
                                             field_id,        &
                                             i_part,          &
                                             weights)
    ! Get spatial interpolation weights (location algorithm).
    use, intrinsic :: iso_c_binding
    implicit none

    character(kind = c_char, len = *) :: local_code_name ! Local code name
    character(kind = c_char, len = *) :: cpl_id          ! Coupling identifier
    character(kind = c_char, len = *) :: field_id        ! Field id
    integer(c_int), intent(in)        :: i_part          ! Partition identifier
    real(8), dimension(:), pointer    :: weights         ! Spatial interpolation weights (barycentric coordinates)
    integer(kind = c_int) :: l_local_code_name, l_cpl_id, l_field_id
    integer(kind = c_int) :: s_weights
    type(c_ptr)           :: c_weights

    l_local_code_name = len(local_code_name)
    l_cpl_id = len(cpl_id)
    l_field_id = len(field_id)

    call CWP_Field_location_weights_get_cf(local_code_name, &
                                           l_local_code_name, &
                                           cpl_id, &
                                           l_cpl_id, &
                                           field_id, &
                                           l_field_id, &
                                           i_part, &
                                           c_weights, &
                                           s_weights)

    call c_f_pointer(c_weights, weights, [s_weights])

  end subroutine CWP_Field_location_weights_get_


  subroutine CWP_Field_location_point_data_get_(local_code_name, &
                                                cpl_id, &
                                                field_id, &
                                                i_part, &
                                                points_coords, &
                                                points_uvw, &
                                                points_dist2, &
                                                points_projected_coords)
    ! Get spatial interpolation point data (location algorithm).
    use, intrinsic :: iso_c_binding
    implicit none

    character(kind = c_char, len = *) :: local_code_name              ! Local code name
    character(kind = c_char, len = *) :: cpl_id                       ! Coupling identifier
    character(kind = c_char, len = *) :: field_id                     ! Field id
    integer(c_int), intent(in)        :: i_part                       ! Partition identifier
    real(8), pointer                  :: points_coords(:,:)           ! Cartesian coordinates of points inside local elements
    real(8), pointer                  :: points_uvw(:,:)              ! Parametric coordinates of points inside local elements
    real(8), pointer                  :: points_dist2(:)              ! Squared distance from points to elements
    real(8), pointer                  :: points_projected_coords(:,:) ! Cartesian coordinates of projection on points on local elements
    integer(kind = c_int) :: l_local_code_name, l_cpl_id, l_field_id
    type(c_ptr)      :: c_points_coords, c_points_uvw, c_points_dist2, c_points_projected_coords
    integer(kind = c_int) :: s_size

    l_local_code_name = len(local_code_name)
    l_cpl_id = len(cpl_id)
    l_field_id = len(field_id)

    call CWP_Field_location_point_data_get_cf(local_code_name, &
                                              l_local_code_name, &
                                              cpl_id, &
                                              l_cpl_id, &
                                              field_id, &
                                              l_field_id, &
                                              i_part, &
                                              c_points_coords, &
                                              c_points_uvw, &
                                              c_points_dist2, &
                                              c_points_projected_coords, &
                                              s_size)

    call c_f_pointer(c_points_coords, points_coords, [3, s_size])
    call c_f_pointer(c_points_uvw, points_uvw, [3, s_size])
    call c_f_pointer(c_points_dist2, points_dist2, [s_size])
    call c_f_pointer(c_points_projected_coords, points_projected_coords, [3, s_size])

  end subroutine CWP_Field_location_point_data_get_


  subroutine CWP_Field_location_internal_cell_vtx_get_(local_code_name, &
                                                       cpl_id, &
                                                       field_id, &
                                                       i_part, &
                                                       cell_vtx_idx, &
                                                       cell_vtx)
    ! Get spatial interpolation internal cell->vertex connectivity (location algorithm).
    use, intrinsic :: iso_c_binding
    implicit none

    character(kind = c_char, len = *)      :: local_code_name ! Local code name
    character(kind = c_char, len = *)      :: cpl_id          ! Coupling identifier
    character(kind = c_char, len = *)      :: field_id        ! Field id
    integer(c_int), intent(in)             :: i_part          ! Partition identifier
    integer(c_int),  dimension(:), pointer :: cell_vtx_idx    ! Index for local cell->vertex connectivity
    integer(c_int),  dimension(:), pointer :: cell_vtx        ! Local cell->vertex connectivity
    integer(kind = c_int) :: l_local_code_name, l_cpl_id, l_field_id
    type(c_ptr)      :: c_cell_vtx_idx, c_cell_vtx
    integer(kind = c_int) :: n_cell

    l_local_code_name = len(local_code_name)
    l_cpl_id = len(cpl_id)
    l_field_id = len(field_id)

    call CWP_Field_location_internal_cell_vtx_get_cf(local_code_name, &
                                                     l_local_code_name, &
                                                     cpl_id, &
                                                     l_cpl_id, &
                                                     field_id, &
                                                     l_field_id, &
                                                     i_part, &
                                                     c_cell_vtx_idx, &
                                                     n_cell, &
                                                     c_cell_vtx)

    call c_f_pointer(c_cell_vtx_idx, cell_vtx_idx, [n_cell+1])
    call c_f_pointer(c_cell_vtx, cell_vtx, [cell_vtx_idx(n_cell+1)])

  end subroutine CWP_Field_location_internal_cell_vtx_get_


  subroutine CWP_Field_intersection_volumes_get_(local_code_name, &
                                                 cpl_id, &
                                                 field_id, &
                                                 i_part, &
                                                 volumes)
    ! Get spatial interpolation volumes (intersection algorithm).
    use, intrinsic :: iso_c_binding
    implicit none

    character(kind = c_char, len = *) :: local_code_name ! Local code name
    character(kind = c_char, len = *) :: cpl_id          ! Coupling identifier
    character(kind = c_char, len = *) :: field_id        ! Field id
    integer(c_int), intent(in)        :: i_part          ! Partition identifier
    real(8),  dimension(:), pointer   :: volumes         ! Volumes of intersection polyhedra
    integer(kind = c_int) :: l_local_code_name, l_cpl_id, l_field_id
    integer(kind = c_int) :: s_volumes
    type(c_ptr)           :: c_volumes

    l_local_code_name = len(local_code_name)
    l_cpl_id = len(cpl_id)
    l_field_id = len(field_id)

    call CWP_Field_intersection_volumes_get_cf(local_code_name, &
                                               l_local_code_name, &
                                               cpl_id, &
                                               l_cpl_id, &
                                               field_id, &
                                               l_field_id, &
                                               i_part, &
                                               c_volumes, &
                                               s_volumes)

    call c_f_pointer(c_volumes, volumes, [s_volumes])

  end subroutine CWP_Field_intersection_volumes_get_


  subroutine CWP_Field_intersection_tgt_elt_volumes_get_(local_code_name, &
                                                         cpl_id,          &
                                                         field_id,        &
                                                         i_part,          &
                                                         tgt_elt_volumes)
    ! Get spatial local target elements volumes (intersection algorithm).
    use, intrinsic :: iso_c_binding
    implicit none

    character(kind = c_char, len = *) :: local_code_name ! Local code name
    character(kind = c_char, len = *) :: cpl_id          ! Coupling identifier
    character(kind = c_char, len = *) :: field_id        ! Field id
    integer(c_int), intent(in)        :: i_part          ! Partition identifier
    real(8),  dimension(:), pointer   :: tgt_elt_volumes ! Volumes of local target elements
    integer(kind = c_int) :: l_local_code_name, l_cpl_id, l_field_id
    type(c_ptr)           :: c_tgt_elt_volumes = C_NULL_PTR
    integer(c_int)        :: n_elt

    l_local_code_name = len(local_code_name)
    l_cpl_id          = len(cpl_id)
    l_field_id        = len(field_id)

    call CWP_Field_intersection_tgt_elt_volumes_get_cf(local_code_name,   &
                                                       l_local_code_name, &
                                                       cpl_id,            &
                                                       l_cpl_id,          &
                                                       field_id,          &
                                                       l_field_id,        &
                                                       i_part,            &
                                                       c_tgt_elt_volumes, &
                                                       n_elt)

    call c_f_pointer(c_tgt_elt_volumes, tgt_elt_volumes, [n_elt])

  end subroutine CWP_Field_intersection_tgt_elt_volumes_get_


  subroutine CWP_Field_nearest_neighbors_distances_get_(local_code_name, &
                                                        cpl_id,          &
                                                        field_id,        &
                                                        i_part,          &
                                                        distances2)
    ! Get spatial interpolation distances (nearest neighbors algorithm).
    use, intrinsic :: iso_c_binding
    implicit none

    character(kind = c_char, len = *) :: local_code_name ! Local code name
    character(kind = c_char, len = *) :: cpl_id          ! Coupling identifier
    character(kind = c_char, len = *) :: field_id        ! Field id
    integer(c_int), intent(in)        :: i_part          ! Partition identifier
    real(8),  dimension(:), pointer   :: distances2      ! Squared distances from nearest source points
    integer(kind = c_int) :: l_local_code_name, l_cpl_id, l_field_id
    integer(kind = c_int) :: s_distances2
    type(c_ptr)           :: c_distances2

    l_local_code_name = len(local_code_name)
    l_cpl_id = len(cpl_id)
    l_field_id = len(field_id)

    call CWP_Field_nearest_neighbors_distances_get_cf(local_code_name, &
                                                      l_local_code_name, &
                                                      cpl_id, &
                                                      l_cpl_id, &
                                                      field_id, &
                                                      l_field_id, &
                                                      i_part, &
                                                      c_distances2, &
                                                      s_distances2)

    call c_f_pointer(c_distances2, distances2, [s_distances2])

  end subroutine CWP_Field_nearest_neighbors_distances_get_


  subroutine CWP_Field_nearest_neighbors_coord_get_(local_code_name,   &
                                                    cpl_id,            &
                                                    field_id,          &
                                                    i_part,            &
                                                    nearest_src_coord)
    ! Get coordinates of nearest source points (nearest neighbors algorithm).
    use, intrinsic :: iso_c_binding
    implicit none

    character(kind = c_char, len = *) :: local_code_name        ! Local code name
    character(kind = c_char, len = *) :: cpl_id                 ! Coupling identifier
    character(kind = c_char, len = *) :: field_id               ! Field id
    integer(c_int), intent(in)        :: i_part                 ! Partition identifier
    real(8), pointer                  :: nearest_src_coord(:,:) ! Coordinates of nearest source points
    integer(c_int) :: l_local_code_name, l_cpl_id, l_field_id

    type(c_ptr)    :: c_nearest_src_coord = C_NULL_PTR
    integer(c_int) :: n_nearest_src_pts

    l_local_code_name = len(local_code_name)
    l_cpl_id          = len(cpl_id)
    l_field_id        = len(field_id)

    call CWP_Field_nearest_neighbors_coord_get_cf(local_code_name,     &
                                                  l_local_code_name,   &
                                                  cpl_id,              &
                                                  l_cpl_id,            &
                                                  field_id,            &
                                                  l_field_id,          &
                                                  i_part,              &
                                                  c_nearest_src_coord, &
                                                  n_nearest_src_pts)

    call c_f_pointer(c_nearest_src_coord, nearest_src_coord, [3, n_nearest_src_pts])

  end subroutine CWP_Field_nearest_neighbors_coord_get_

! /*----------------------------------------------------------------------------*
!  * Functions about code parameters                                            *
!  *----------------------------------------------------------------------------*/


  !>
  !!
  !! \brief Add a new parameter and initialize it.
  !!
  !! \param [in] local_code_name  Local code name
  !! \param [in] param_name       Parameter name
  !! \param [in] data_type        Parameter type
  !! \param [in] initial_value    Initial value
  !!
  !!

  subroutine CWP_Param_add_int_(local_code_name, &
                                param_name,      &
                                initial_value)

    use, intrinsic :: iso_c_binding
    implicit none

    character(kind = c_char, len = *) :: local_code_name
    character(kind = c_char, len = *) :: param_name
    integer(kind = c_int), intent(in) :: initial_value
    integer(kind = c_int)             :: l_local_code_name, l_param_name

    l_local_code_name = len(local_code_name)
    l_param_name      = len(param_name)

    call CWP_Param_add_int_cf(local_code_name,   &
                              l_local_code_name, &
                              param_name,        &
                              l_param_name,      &
                              initial_value)

  end subroutine CWP_Param_add_int_


  subroutine CWP_Param_add_double_(local_code_name, &
                                   param_name,      &
                                   initial_value)

    use, intrinsic :: iso_c_binding
    implicit none

    character(kind = c_char, len = *) :: local_code_name
    character(kind = c_char, len = *) :: param_name
    real(kind = c_double), intent(in) :: initial_value
    integer(kind = c_int)             :: l_local_code_name, l_param_name

    l_local_code_name = len(local_code_name)
    l_param_name      = len(param_name)

    call CWP_Param_add_double_cf(local_code_name,   &
                                 l_local_code_name, &
                                 param_name,        &
                                 l_param_name,      &
                                 initial_value)

  end subroutine CWP_Param_add_double_

  subroutine CWP_Param_add_char_(local_code_name, &
                                 param_name,      &
                                 initial_value)

    use, intrinsic :: iso_c_binding
    implicit none

    character(kind = c_char, len = *) :: local_code_name
    character(kind = c_char, len = *) :: param_name
    character(kind = c_char, len = *) :: initial_value
    integer(kind = c_int)             :: l_local_code_name, l_param_name, l_initial_value

    l_local_code_name = len(local_code_name)
    l_param_name      = len(param_name)
    l_initial_value   = len(initial_value)

    call CWP_Param_add_char_cf(local_code_name,   &
                               l_local_code_name, &
                               param_name,        &
                               l_param_name,      &
                               initial_value,     &
                               l_initial_value)

  end subroutine CWP_Param_add_char_

  !>
  !!
  !! \brief Set a parameter.
  !!
  !! \param [in] local_code_name  Local code name
  !! \param [in] param_name       Parameter name
  !! \param [in] data_type        Parameter type
  !! \param [in] value            Value
  !!
  !!

  subroutine CWP_Param_set_int_(local_code_name, &
                                param_name,      &
                                value)

    use, intrinsic :: iso_c_binding
    implicit none

    character(kind = c_char, len = *) :: local_code_name
    character(kind = c_char, len = *) :: param_name
    integer(kind = c_int), intent(in) :: value
    integer(kind = c_int)             :: l_local_code_name, l_param_name

    l_local_code_name = len(local_code_name)
    l_param_name      = len(param_name)

    call CWP_Param_set_int_cf(local_code_name,   &
                              l_local_code_name, &
                              param_name,        &
                              l_param_name,      &
                              value)

  end subroutine CWP_Param_set_int_


  subroutine CWP_Param_set_double_(local_code_name, &
                                   param_name,      &
                                   value)

    use, intrinsic :: iso_c_binding
    implicit none

    character(kind = c_char, len = *) :: local_code_name
    character(kind = c_char, len = *) :: param_name
    real(kind = c_double), intent(in) :: value
    integer(kind = c_int)             :: l_local_code_name, l_param_name

    l_local_code_name = len(local_code_name)
    l_param_name      = len(param_name)

    call CWP_Param_set_double_cf(local_code_name,   &
                                 l_local_code_name, &
                                 param_name,        &
                                 l_param_name,      &
                                 value)

  end subroutine CWP_Param_set_double_

  subroutine CWP_Param_set_char_(local_code_name, &
                                 param_name,      &
                                 value)

    use, intrinsic :: iso_c_binding
    implicit none

    character(kind = c_char, len = *) :: local_code_name
    character(kind = c_char, len = *) :: param_name
    character(kind = c_char, len = *) :: value
    integer(kind = c_int)             :: l_local_code_name, l_param_name, l_value

    l_local_code_name = len(local_code_name)
    l_param_name      = len(param_name)
    l_value           = len(value)

    call CWP_Param_set_char_cf(local_code_name,   &
                               l_local_code_name, &
                               param_name,        &
                               l_param_name,      &
                               value,             &
                               l_value)

  end subroutine CWP_Param_set_char_


  subroutine CWP_Param_del_(local_code_name, &
                            param_name,      &
                            data_type)
    ! Delete a control parameter.
    use, intrinsic :: iso_c_binding
    implicit none

    character(kind = c_char, len = *) :: local_code_name ! Local code name
    character(kind = c_char, len = *) :: param_name      ! Parameter name
    integer, intent(in)               :: data_type       ! Parameter type
    integer(kind = c_int)             :: l_local_code_name, l_param_name

    l_local_code_name = len(local_code_name)
    l_param_name      = len(param_name)

    call CWP_Param_del_cf(local_code_name,   &
                          l_local_code_name, &
                          param_name,        &
                          l_param_name,      &
                          data_type)

  end subroutine CWP_Param_del_


  function CWP_Param_n_get_(code_name, &
                            data_type) &
    result (n_param)
    ! Return the number of control parameters for the code ``code_name``.
    use, intrinsic :: iso_c_binding
    implicit none

    character(kind = c_char, len = *) :: code_name ! Local or distant code name
    integer, intent(in)               :: data_type ! Parameter type
    integer                           :: n_param   ! Number of control parameters
    integer(kind = c_int)             :: l_code_name

    l_code_name = len(code_name)

    n_param = CWP_Param_n_get_cf(code_name,   &
                                 l_code_name, &
                                 data_type)

  end function CWP_Param_n_get_


  subroutine CWP_Param_list_get_(code_name, &
                                 data_type, &
                                 n_param,   &
                                 param_names)
    ! Return the list of control parameters for the code ``code_name``.
    use, intrinsic :: iso_c_binding
    implicit none

    character(kind = c_char, len = *) :: code_name      ! Local or distant code name
    integer                           :: data_type      ! Parameter type
    integer(c_int), intent(out)       :: n_param        ! Number of parameters
    character(256), allocatable       :: param_names(:) ! Parameter names
    integer(kind = c_int)             :: l_code_name
    type(c_ptr)                       :: c_param_names, c_param_sizes

    l_code_name  = len(code_name)

    call CWP_Param_list_get_cf(code_name, l_code_name, data_type, n_param, c_param_names, c_param_sizes)

    call c_f_char_array(c_param_names, c_param_sizes, n_param, param_names, .true.)

  end subroutine CWP_Param_list_get_


  function CWP_Param_is_(code_name,  &
                         param_name, &
                         data_type)  &
    result (is_param)
    ! Is this an existing control parameter for code ``code_name`` ?
    use, intrinsic :: iso_c_binding
    implicit none

    character(kind = c_char, len = *) :: code_name  ! Local or distant code name
    character(kind = c_char, len = *) :: param_name ! Parameter name
    integer, intent(in)               :: data_type  ! Parameter type
    logical                           :: is_param   ! Existing status
    integer(kind = c_int)             :: l_code_name
    integer(kind = c_int)             :: l_param_name
    integer(kind = c_int)             :: c_is_param

    l_code_name  = len(code_name)
    l_param_name = len(param_name)

    c_is_param = CWP_Param_is_cf(code_name,    &
                                 l_code_name,  &
                                 param_name,   &
                                 l_param_name, &
                                 data_type)

    if (c_is_param .eq. 1) then
      is_param = .true.
    else
      is_param = .false.
    endif

  end function CWP_Param_is_


  !>
  !!
  !! \brief Return the parameter value of \p param_name on \p code_name.
  !!
  !! \param [in]  code_name  Local or distant code name
  !! \param [in]  param_name Parameter name
  !! \param [in]  data_type  Parameter type
  !! \param [out] value      Parameter value
  !!
  !!

  subroutine CWP_Param_get_int(code_name,  &
                               param_name, &
                               value)

    use, intrinsic :: iso_c_binding
    implicit none

    character(kind = c_char, len = *) :: code_name
    character(kind = c_char, len = *) :: param_name
    integer(c_int), intent(out)       :: value
    integer(kind = c_int)             :: l_code_name
    integer(kind = c_int)             :: l_param_name

    l_code_name  = len(code_name)
    l_param_name = len(param_name)

    call CWP_Param_get_int_cf(code_name,    &
                              l_code_name,  &
                              param_name,   &
                              l_param_name, &
                              value)

  end subroutine CWP_Param_get_int


  subroutine CWP_Param_get_double(code_name,  &
                                  param_name, &
                                  value)

    use, intrinsic :: iso_c_binding
    implicit none

    character(kind = c_char, len = *) :: code_name
    character(kind = c_char, len = *) :: param_name
    real(c_double), intent(out)       :: value
    integer(kind = c_int)             :: l_code_name
    integer(kind = c_int)             :: l_param_name

    l_code_name  = len(code_name)
    l_param_name = len(param_name)

    call CWP_Param_get_double_cf(code_name,    &
                                 l_code_name,  &
                                 param_name,   &
                                 l_param_name, &
                                 value)

  end subroutine CWP_Param_get_double

  subroutine CWP_Param_get_char(code_name,  &
                                param_name, &
                                val)

    use, intrinsic :: iso_c_binding
    implicit none

    character(kind = c_char, len = *) :: code_name
    character(kind = c_char, len = *) :: param_name
    integer(kind = c_int)             :: l_code_name
    integer(kind = c_int)             :: l_param_name
    character(c_char), pointer        :: val(:)
    type(c_ptr)                       :: cptr = C_NULL_PTR
    integer(c_int)                    :: l_value

    l_code_name  = len(code_name)
    l_param_name = len(param_name)

    call CWP_Param_get_char_cf(code_name,    &
                               l_code_name,  &
                               param_name,   &
                               l_param_name, &
                               cptr,         &
                               l_value)

    call c_f_pointer(cptr, val, [l_value])

  end subroutine CWP_Param_get_char


  !>
  !!
  !! \brief Return the result of a reduce operation about a parameter
  !!
  !! \param [in]  op           Operation
  !! \param [in]  param_name   Parameter name
  !! \param [in]  data_type    Parameter type
  !! \param [out] res          Result
  !! \param [in]  n_codes      Number of codes
  !! \param [in]  code_names   Codes name
  !!

  subroutine CWP_Param_reduce_int(op,         &
                                  param_name, &
                                  res,        &
                                  n_codes,    &
                                  code_names)

    use, intrinsic :: iso_c_binding
    implicit none

    integer, intent(in)                                           :: op
    integer(c_int), intent(out)                                   :: res
    integer(c_int)                                                :: n_codes, i
    character(kind = c_char, len = *)                             :: param_name
    integer(kind = c_int)                                         :: l_param_name
    character(kind = c_char, len = *), dimension(n_codes), target :: code_names
    integer, dimension(n_codes), target                           :: l_code_names

    l_param_name = len(param_name)
    do i=1,n_codes
      l_code_names(i) = len(code_names(i))
    end do

    call CWP_Param_reduce_int_cf(op,                &
                                 param_name,        &
                                 l_param_name,      &
                                 res,               &
                                 n_codes,           &
                                 c_loc(code_names), &
                                 c_loc(l_code_names))

  end subroutine CWP_Param_reduce_int

  subroutine CWP_Param_reduce_double(op,         &
                                     param_name, &
                                     res,        &
                                     n_codes,    &
                                     code_names)

    use, intrinsic :: iso_c_binding
    implicit none

    integer, intent(in)                                           :: op
    real(c_double), intent(out)                                   :: res
    integer(c_int)                                                :: n_codes, i
    character(kind = c_char, len = *)                             :: param_name
    integer(kind = c_int)                                         :: l_param_name
    character(kind = c_char, len = *), dimension(n_codes), target :: code_names
    integer, dimension(n_codes), target                           :: l_code_names

    l_param_name = len(param_name)
    do i=1,n_codes
      l_code_names(i) = len(code_names(i))
    end do

    call CWP_Param_reduce_double_cf(op,                &
                                    param_name,        &
                                    l_param_name,      &
                                    res,               &
                                    n_codes,           &
                                    c_loc(code_names), &
                                    c_loc(l_code_names))

  end subroutine CWP_Param_reduce_double

  subroutine CWP_Param_reduce_char(op,         &
                                   param_name, &
                                   res,        &
                                   n_codes,    &
                                   code_names)

    use, intrinsic :: iso_c_binding
    implicit none

    integer, intent(in)                                           :: op
    character(c_char), pointer                                    :: res(:)
    integer(c_int)                                                :: n_codes, i
    character(kind = c_char, len = *)                             :: param_name
    integer(kind = c_int)                                         :: l_param_name
    character(kind = c_char, len = *), dimension(n_codes), target :: code_names
    integer, dimension(n_codes), target                           :: l_code_names
    type(c_ptr)                                                   :: cptr = C_NULL_PTR
    integer(c_int)                                                :: l_res

    l_param_name = len(param_name)
    do i=1,n_codes
      l_code_names(i) = len(code_names(i))
    end do

    call CWP_Param_reduce_char_cf(op,                &
                                  param_name,        &
                                  l_param_name,      &
                                  cptr,              &
                                  l_res,             &
                                  n_codes,           &
                                  c_loc(code_names), &
                                  c_loc(l_code_names))

    call c_f_pointer(cptr, res, [l_res])

  end subroutine CWP_Param_reduce_char


  subroutine CWP_Param_lock_(code_name)
    ! Lock access to local control parameters from a distant code.
    use, intrinsic :: iso_c_binding
    implicit none

    character(kind = c_char, len = *) :: code_name ! Code to lock
    integer(kind = c_int)             :: l_code_name

    l_code_name  = len(code_name)

    call CWP_Param_lock_cf(code_name,  &
                           l_code_name)

  end subroutine CWP_Param_lock_


  subroutine CWP_Param_unlock_(code_name)
    ! Unlock access to local control parameters from a distant code.
    use, intrinsic :: iso_c_binding
    implicit none

    character(kind = c_char, len = *) :: code_name ! Code to unlock
    integer(kind = c_int)             :: l_code_name

    l_code_name  = len(code_name)

    call CWP_Param_unlock_cf(code_name,   &
                             l_code_name)

  end subroutine CWP_Param_unlock_


  !>
  !! \brief Initiate the sending of a data array.
  !!
  !! \param [in] local_code_name  Local code name
  !! \param [in] cpl_id           Coupling identifier
  !! \param [in] global_data_id   GlobalData identifier
  !! \param [in] send_data        Pointer to data array
  !!
  !!

  subroutine CWP_Global_data_issend_int(local_code_name, &
                                        cpl_id,          &
                                        global_data_id,  &
                                        send_data)
    use, intrinsic :: iso_c_binding
    implicit none
    character(kind=c_char, len=*) :: local_code_name, cpl_id, global_data_id
    integer(c_int), pointer       :: send_data(:,:)

    integer(c_int)                :: l_local_code_name, l_cpl_id, l_global_data_id
    integer(c_long)               :: s_data
    integer(c_int)                :: send_stride
    integer(c_int)                :: n_send_entity

    l_local_code_name = len(local_code_name)
    l_cpl_id          = len(cpl_id)
    l_global_data_id  = len(global_data_id)

    s_data        = 4
    send_stride   = size(send_data, 1)
    n_send_entity = size(send_data, 2)

    call CWP_Global_data_issend_cf(local_code_name,   &
                                   l_local_code_name, &
                                   cpl_id,            &
                                   l_cpl_id,          &
                                   global_data_id,    &
                                   l_global_data_id,  &
                                   s_data,            &
                                   send_stride,       &
                                   n_send_entity,     &
                                   c_loc(send_data))

  end subroutine CWP_Global_data_issend_int

  subroutine CWP_Global_data_issend_long(local_code_name, &
                                         cpl_id,          &
                                         global_data_id,  &
                                         send_data)
    use, intrinsic :: iso_c_binding
    implicit none
    character(kind=c_char, len=*) :: local_code_name, cpl_id, global_data_id
    integer(c_long), pointer      :: send_data(:,:)

    integer(c_int)                :: l_local_code_name, l_cpl_id, l_global_data_id
    integer(c_long)               :: s_data
    integer(c_int)                :: send_stride
    integer(c_int)                :: n_send_entity

    l_local_code_name = len(local_code_name)
    l_cpl_id          = len(cpl_id)
    l_global_data_id  = len(global_data_id)

    s_data        = 8
    send_stride   = size(send_data, 1)
    n_send_entity = size(send_data, 2)

    call CWP_Global_data_issend_cf(local_code_name,   &
                                   l_local_code_name, &
                                   cpl_id,            &
                                   l_cpl_id,          &
                                   global_data_id,    &
                                   l_global_data_id,  &
                                   s_data,            &
                                   send_stride,       &
                                   n_send_entity,     &
                                   c_loc(send_data))

  end subroutine CWP_Global_data_issend_long

  subroutine CWP_Global_data_issend_double(local_code_name, &
                                           cpl_id,          &
                                           global_data_id,  &
                                           send_data)
    use, intrinsic :: iso_c_binding
    implicit none
    character(kind=c_char, len=*) :: local_code_name, cpl_id, global_data_id
    real(8), pointer     :: send_data(:,:)

    integer(c_int)                :: l_local_code_name, l_cpl_id, l_global_data_id
    integer(c_long)               :: s_data
    integer(c_int)                :: send_stride
    integer(c_int)                :: n_send_entity

    l_local_code_name = len(local_code_name)
    l_cpl_id          = len(cpl_id)
    l_global_data_id  = len(global_data_id)

    s_data        = 8
    send_stride   = size(send_data, 1)
    n_send_entity = size(send_data, 2)

    call CWP_Global_data_issend_cf(local_code_name,   &
                                   l_local_code_name, &
                                   cpl_id,            &
                                   l_cpl_id,          &
                                   global_data_id,    &
                                   l_global_data_id,  &
                                   s_data,            &
                                   send_stride,       &
                                   n_send_entity,     &
                                   c_loc(send_data))

  end subroutine CWP_Global_data_issend_double

  subroutine CWP_Global_data_issend_complex4(local_code_name, &
                                             cpl_id,          &
                                             global_data_id,  &
                                             send_data)
    use, intrinsic :: iso_c_binding
    implicit none
    character(kind=c_char, len=*) :: local_code_name, cpl_id, global_data_id
    complex(kind = 4), pointer    :: send_data(:,:)

    integer(c_int)                :: l_local_code_name, l_cpl_id, l_global_data_id
    integer(c_long)               :: s_data
    integer(c_int)                :: send_stride
    integer(c_int)                :: n_send_entity

    l_local_code_name = len(local_code_name)
    l_cpl_id          = len(cpl_id)
    l_global_data_id  = len(global_data_id)

    s_data        = 8
    send_stride   = size(send_data, 1)
    n_send_entity = size(send_data, 2)

    call CWP_Global_data_issend_cf(local_code_name,   &
                                   l_local_code_name, &
                                   cpl_id,            &
                                   l_cpl_id,          &
                                   global_data_id,    &
                                   l_global_data_id,  &
                                   s_data,            &
                                   send_stride,       &
                                   n_send_entity,     &
                                   c_loc(send_data))

  end subroutine CWP_Global_data_issend_complex4

  subroutine CWP_Global_data_issend_complex8(local_code_name, &
                                             cpl_id,          &
                                             global_data_id,  &
                                             send_data)
    use, intrinsic :: iso_c_binding
    implicit none
    character(kind=c_char, len=*) :: local_code_name, cpl_id, global_data_id
    complex(kind = 8), pointer    :: send_data(:,:)

    integer(c_int)                :: l_local_code_name, l_cpl_id, l_global_data_id
    integer(c_long)               :: s_data
    integer(c_int)                :: send_stride
    integer(c_int)                :: n_send_entity

    l_local_code_name = len(local_code_name)
    l_cpl_id          = len(cpl_id)
    l_global_data_id  = len(global_data_id)

    s_data        = 16
    send_stride   = size(send_data, 1)
    n_send_entity = size(send_data, 2)

    call CWP_Global_data_issend_cf(local_code_name,   &
                                   l_local_code_name, &
                                   cpl_id,            &
                                   l_cpl_id,          &
                                   global_data_id,    &
                                   l_global_data_id,  &
                                   s_data,            &
                                   send_stride,       &
                                   n_send_entity,     &
                                   c_loc(send_data))

  end subroutine CWP_Global_data_issend_complex8

  subroutine CWP_Global_data_issend_real4(local_code_name, &
                                          cpl_id,          &
                                          global_data_id,  &
                                          send_data)
    use, intrinsic :: iso_c_binding
    implicit none
    character(kind=c_char, len=*) :: local_code_name, cpl_id, global_data_id
    real(kind = 4), pointer       :: send_data(:,:)

    integer(c_int)                :: l_local_code_name, l_cpl_id, l_global_data_id
    integer(c_long)               :: s_data
    integer(c_int)                :: send_stride
    integer(c_int)                :: n_send_entity

    l_local_code_name = len(local_code_name)
    l_cpl_id          = len(cpl_id)
    l_global_data_id  = len(global_data_id)

    s_data        = 4
    send_stride   = size(send_data, 1)
    n_send_entity = size(send_data, 2)

    call CWP_Global_data_issend_cf(local_code_name,   &
                                   l_local_code_name, &
                                   cpl_id,            &
                                   l_cpl_id,          &
                                   global_data_id,    &
                                   l_global_data_id,  &
                                   s_data,            &
                                   send_stride,       &
                                   n_send_entity,     &
                                   c_loc(send_data))

  end subroutine CWP_Global_data_issend_real4

  !>
  !! \brief Initiate the reception of a data array.
  !!
  !! \param [in] local_code_name  Local code name
  !! \param [in] cpl_id           Coupling identifier
  !! \param [in] global_data_id   GlobalData identifier
  !! \param [in] recv_data        Pointer to data array
  !!
  !!

  subroutine CWP_Global_data_irecv_int(local_code_name, &
                                       cpl_id,          &
                                       global_data_id,  &
                                       recv_data)
    use, intrinsic :: iso_c_binding
    implicit none
    character(kind=c_char, len=*) :: local_code_name, cpl_id, global_data_id
    integer(c_int), pointer       :: recv_data(:,:)

    integer(c_int)                :: l_local_code_name, l_cpl_id, l_global_data_id
    integer(c_long)               :: s_data
    integer(c_int)                :: recv_stride
    integer(c_int)                :: n_recv_entity

    l_local_code_name = len(local_code_name)
    l_cpl_id          = len(cpl_id)
    l_global_data_id  = len(global_data_id)

    s_data = 4
    recv_stride   = size(recv_data, 1)
    n_recv_entity = size(recv_data, 2)

    call CWP_Global_data_irecv_cf(local_code_name,   &
                                  l_local_code_name, &
                                  cpl_id,            &
                                  l_cpl_id,          &
                                  global_data_id,    &
                                  l_global_data_id,  &
                                  s_data,            &
                                  recv_stride,       &
                                  n_recv_entity,     &
                                  c_loc(recv_data))

  end subroutine CWP_Global_data_irecv_int

  subroutine CWP_Global_data_irecv_long(local_code_name, &
                                        cpl_id,          &
                                        global_data_id,  &
                                        recv_data)
    use, intrinsic :: iso_c_binding
    implicit none
    character(kind=c_char, len=*) :: local_code_name, cpl_id, global_data_id
    integer(c_long), pointer      :: recv_data(:,:)

    integer(c_int)                :: l_local_code_name, l_cpl_id, l_global_data_id
    integer(c_long)               :: s_data
    integer(c_int)                :: recv_stride
    integer(c_int)                :: n_recv_entity

    l_local_code_name = len(local_code_name)
    l_cpl_id          = len(cpl_id)
    l_global_data_id  = len(global_data_id)

    s_data = 8
    recv_stride   = size(recv_data, 1)
    n_recv_entity = size(recv_data, 2)

    call CWP_Global_data_irecv_cf(local_code_name,   &
                                  l_local_code_name, &
                                  cpl_id,            &
                                  l_cpl_id,          &
                                  global_data_id,    &
                                  l_global_data_id,  &
                                  s_data,            &
                                  recv_stride,       &
                                  n_recv_entity,     &
                                  c_loc(recv_data))

  end subroutine CWP_Global_data_irecv_long

  subroutine CWP_Global_data_irecv_double(local_code_name, &
                                          cpl_id,          &
                                          global_data_id,  &
                                          recv_data)
    use, intrinsic :: iso_c_binding
    implicit none
    character(kind=c_char, len=*) :: local_code_name, cpl_id, global_data_id
    real(8), pointer     :: recv_data(:,:)

    integer(c_int)                :: l_local_code_name, l_cpl_id, l_global_data_id
    integer(c_long)               :: s_data
    integer(c_int)                :: recv_stride
    integer(c_int)                :: n_recv_entity

    l_local_code_name = len(local_code_name)
    l_cpl_id          = len(cpl_id)
    l_global_data_id  = len(global_data_id)

    s_data = 8
    recv_stride   = size(recv_data, 1)
    n_recv_entity = size(recv_data, 2)

    call CWP_Global_data_irecv_cf(local_code_name,   &
                                  l_local_code_name, &
                                  cpl_id,            &
                                  l_cpl_id,          &
                                  global_data_id,    &
                                  l_global_data_id,  &
                                  s_data,            &
                                  recv_stride,       &
                                  n_recv_entity,     &
                                  c_loc(recv_data))

  end subroutine CWP_Global_data_irecv_double

  subroutine CWP_Global_data_irecv_complex4(local_code_name, &
                                            cpl_id,          &
                                            global_data_id,  &
                                            recv_data)
    use, intrinsic :: iso_c_binding
    implicit none
    character(kind=c_char, len=*) :: local_code_name, cpl_id, global_data_id
    complex(kind=4), pointer      :: recv_data(:,:)

    integer(c_int)                :: l_local_code_name, l_cpl_id, l_global_data_id
    integer(c_long)               :: s_data
    integer(c_int)                :: recv_stride
    integer(c_int)                :: n_recv_entity

    l_local_code_name = len(local_code_name)
    l_cpl_id          = len(cpl_id)
    l_global_data_id  = len(global_data_id)

    s_data = 8
    recv_stride   = size(recv_data, 1)
    n_recv_entity = size(recv_data, 2)

    call CWP_Global_data_irecv_cf(local_code_name, &
                                  l_local_code_name, &
                                  cpl_id,          &
                                  l_cpl_id,          &
                                  global_data_id,  &
                                  l_global_data_id,  &
                                  s_data,            &
                                  recv_stride,       &
                                  n_recv_entity,     &
                                  c_loc(recv_data))

  end subroutine CWP_Global_data_irecv_complex4

  subroutine CWP_Global_data_irecv_complex8(local_code_name, &
                                            cpl_id,          &
                                            global_data_id,  &
                                            recv_data)
    use, intrinsic :: iso_c_binding
    implicit none
    character(kind=c_char, len=*) :: local_code_name, cpl_id, global_data_id
    complex(kind=8), pointer      :: recv_data(:,:)

    integer(c_int)                :: l_local_code_name, l_cpl_id, l_global_data_id
    integer(c_long)               :: s_data
    integer(c_int)                :: recv_stride
    integer(c_int)                :: n_recv_entity

    l_local_code_name = len(local_code_name)
    l_cpl_id          = len(cpl_id)
    l_global_data_id  = len(global_data_id)

    s_data = 16
    recv_stride   = size(recv_data, 1)
    n_recv_entity = size(recv_data, 2)

    call CWP_Global_data_irecv_cf(local_code_name,   &
                                  l_local_code_name, &
                                  cpl_id,            &
                                  l_cpl_id,          &
                                  global_data_id,    &
                                  l_global_data_id,  &
                                  s_data,            &
                                  recv_stride,       &
                                  n_recv_entity,     &
                                  c_loc(recv_data))

  end subroutine CWP_Global_data_irecv_complex8

  subroutine CWP_Global_data_irecv_real4(local_code_name, &
                                         cpl_id,          &
                                         global_data_id,  &
                                         recv_data)
    use, intrinsic :: iso_c_binding
    implicit none
    character(kind=c_char, len=*) :: local_code_name, cpl_id, global_data_id
    real(kind=4), pointer         :: recv_data(:,:)

    integer(c_int)                :: l_local_code_name, l_cpl_id, l_global_data_id
    integer(c_long)               :: s_data
    integer(c_int)                :: recv_stride
    integer(c_int)                :: n_recv_entity

    l_local_code_name = len(local_code_name)
    l_cpl_id          = len(cpl_id)
    l_global_data_id  = len(global_data_id)

    s_data = 8
    recv_stride   = size(recv_data, 1)
    n_recv_entity = size(recv_data, 2)

    call CWP_Global_data_irecv_cf(local_code_name,   &
                                  l_local_code_name, &
                                  cpl_id,            &
                                  l_cpl_id,          &
                                  global_data_id,    &
                                  l_global_data_id,  &
                                  s_data,            &
                                  recv_stride,       &
                                  n_recv_entity,     &
                                  c_loc(recv_data))

  end subroutine CWP_Global_data_irecv_real4


  subroutine CWP_Global_data_wait_issend_(local_code_name, &
                                          cpl_id,          &
                                          global_data_id)
    ! Finalize the sending of a global data array.
    use, intrinsic :: iso_c_binding
    implicit none
    character(kind=c_char, len=*) :: local_code_name ! Local code name
    character(kind=c_char, len=*) :: cpl_id          ! Coupling identifier
    character(kind=c_char, len=*) :: global_data_id  ! GlobalData identifier

    integer(c_int)                :: l_local_code_name, l_cpl_id, l_global_data_id

    l_local_code_name = len(local_code_name)
    l_cpl_id          = len(cpl_id)
    l_global_data_id  = len(global_data_id)

    call CWP_Global_data_wait_issend_cf(local_code_name,   &
                                        l_local_code_name, &
                                        cpl_id,            &
                                        l_cpl_id,          &
                                        global_data_id,    &
                                        l_global_data_id)

  end subroutine CWP_Global_data_wait_issend_


  subroutine CWP_Global_data_wait_irecv_(local_code_name, &
                                         cpl_id,          &
                                         global_data_id)
    ! Finalize the reception of a global data array.
    use, intrinsic :: iso_c_binding
    implicit none
    character(kind=c_char, len=*) :: local_code_name ! Local code name
    character(kind=c_char, len=*) :: cpl_id          ! Coupling identifier
    character(kind=c_char, len=*) :: global_data_id  ! GlobalData identifier

    integer(c_int)                :: l_local_code_name, l_cpl_id, l_global_data_id

    l_local_code_name = len(local_code_name)
    l_cpl_id          = len(cpl_id)
    l_global_data_id  = len(global_data_id)

    call CWP_Global_data_wait_irecv_cf(local_code_name,   &
                                       l_local_code_name, &
                                       cpl_id,            &
                                       l_cpl_id,          &
                                       global_data_id,    &
                                       l_global_data_id)

  end subroutine CWP_Global_data_wait_irecv_


  subroutine CWP_Part_data_create_(local_code_name, &
                                   cpl_id,          &
                                   part_data_id,    &
                                   exch_type,       &
                                   gnum_elt,        &
                                   n_elt,           &
                                   n_part)
    ! Create a partitioned data exchange object.
    use, intrinsic :: iso_c_binding
    use pdm_pointer_array
    implicit none
    character(kind=c_char, len=*)     :: local_code_name ! Local code name
    character(kind=c_char, len=*)     :: cpl_id          ! Coupling identifier
    character(kind=c_char, len=*)     :: part_data_id    ! PartData identifier
    integer(c_int), intent(in)        :: exch_type       ! Exchange type
    type(PDM_pointer_array_t), target :: gnum_elt        ! Global ids
    integer(c_int), pointer           :: n_elt(:)        ! Number of elements in partitions (size = ``n_part``)
    integer,        intent(in)        :: n_part          ! Number of partitions

    integer(c_int)                    :: l_local_code_name, l_cpl_id, l_part_data_id

    l_local_code_name = len(local_code_name)
    l_cpl_id          = len(cpl_id)
    l_part_data_id    = len(part_data_id)

    call CWP_Part_data_create_cf(local_code_name,      &
                                 l_local_code_name,    &
                                 cpl_id,               &
                                 l_cpl_id,             &
                                 part_data_id,         &
                                 l_part_data_id,       &
                                 exch_type,            &
                                 c_loc(gnum_elt%cptr), &
                                 c_loc(n_elt),         &
                                 n_part)

  end subroutine CWP_Part_data_create_


  subroutine CWP_Part_data_del_(local_code_name, &
                                cpl_id,          &
                                part_data_id)
    ! Delete a partitioned data exchange object.
    use, intrinsic :: iso_c_binding
    implicit none
    character(kind=c_char, len=*)     :: local_code_name ! Local code name
    character(kind=c_char, len=*)     :: cpl_id          ! Coupling identifier
    character(kind=c_char, len=*)     :: part_data_id    ! PartData identifier

    integer(c_int)                :: l_local_code_name, l_cpl_id, l_part_data_id

    l_local_code_name = len(local_code_name)
    l_cpl_id          = len(cpl_id)
    l_part_data_id    = len(part_data_id)

    call CWP_Part_data_del_cf(local_code_name,   &
                              l_local_code_name, &
                              cpl_id,            &
                              l_cpl_id,          &
                              part_data_id,      &
                              l_part_data_id)

  end subroutine CWP_Part_data_del_


  subroutine CWP_Part_data_issend_(local_code_name, &
                                   cpl_id,          &
                                   part_data_id,    &
                                   exch_id,         &
                                   n_components,    &
                                   send_data)
    ! Initiate the sending of a partitioned data array.
    use, intrinsic :: iso_c_binding
    use pdm_pointer_array
    implicit none

    character(kind=c_char, len=*)     :: local_code_name ! Local code name
    character(kind=c_char, len=*)     :: cpl_id          ! Coupling identifier
    character(kind=c_char, len=*)     :: part_data_id    ! PartData identifier
    integer(c_int), intent(in)        :: exch_id         ! Exchange identifier
    integer(c_int), intent(in)        :: n_components    ! Number of components
    type(PDM_pointer_array_t), target :: send_data       ! Pointer to data array to send

    integer(c_int)                    :: l_local_code_name, l_cpl_id, l_part_data_id
    integer(c_long)                   :: s_data

    l_local_code_name = len(local_code_name)
    l_cpl_id          = len(cpl_id)
    l_part_data_id    = len(part_data_id)

    s_data = send_data%s_data

    call CWP_Part_data_issend_cf(local_code_name,       &
                                 l_local_code_name,     &
                                 cpl_id,                &
                                 l_cpl_id,              &
                                 part_data_id,          &
                                 l_part_data_id,        &
                                 exch_id,               &
                                 s_data,                &
                                 n_components,          &
                                 c_loc(send_data%cptr))

  end subroutine CWP_Part_data_issend_


  subroutine CWP_Part_data_irecv_(local_code_name, &
                                  cpl_id,          &
                                  part_data_id,    &
                                  exch_id,         &
                                  n_components,    &
                                  recv_data)
    ! Initiate the reception of a partitioned data array.
    use, intrinsic :: iso_c_binding
    use pdm_pointer_array
    implicit none

    character(kind=c_char, len=*)     :: local_code_name ! Local code name
    character(kind=c_char, len=*)     :: cpl_id          ! Coupling identifier
    character(kind=c_char, len=*)     :: part_data_id    ! PartData identifier
    integer(c_int), intent(in)        :: exch_id         ! Exchange identifier
    integer(c_int), intent(in)        :: n_components    ! Number of components
    type(PDM_pointer_array_t), target :: recv_data       ! Pointer to received data

    integer(c_int)                    :: l_local_code_name, l_cpl_id, l_part_data_id
    integer(c_long)                   :: s_data
    integer(c_int)                    :: n_part, n_ref, i

    l_local_code_name = len(local_code_name)
    l_cpl_id          = len(cpl_id)
    l_part_data_id    = len(part_data_id)

    s_data = recv_data%s_data

    call CWP_Part_data_irecv_cf(local_code_name,        &
                                l_local_code_name,      &
                                cpl_id,                 &
                                l_cpl_id,               &
                                part_data_id,           &
                                l_part_data_id,         &
                                exch_id,                &
                                s_data,                 &
                                n_components,           &
                                c_loc(recv_data%cptr))

    call CWP_Part_data_n_part_get_cf(local_code_name,   &
                                     l_local_code_name, &
                                     cpl_id,            &
                                     l_cpl_id,          &
                                     part_data_id,      &
                                     l_part_data_id,    &
                                     n_part)

    ! Compute lengths
    do i = 1, n_part
      call CWP_Part_data_n_ref_get_cf(local_code_name,   &
                                      l_local_code_name, &
                                      cpl_id,            &
                                      l_cpl_id,          &
                                      part_data_id,      &
                                      l_part_data_id,    &
                                      i-1,               &
                                      n_ref)

      recv_data%length(i) = n_components * n_ref
    enddo

  end subroutine CWP_Part_data_irecv_


  subroutine CWP_Part_data_wait_issend_(local_code_name, &
                                        cpl_id,          &
                                        part_data_id,    &
                                        exch_id)
    ! Finalize the sending of a partitioned data array.
    use, intrinsic :: iso_c_binding
    implicit none
    character(kind=c_char, len=*) :: local_code_name ! Local code name
    character(kind=c_char, len=*) :: cpl_id          ! Coupling identifier
    character(kind=c_char, len=*) :: part_data_id    ! PartData identifier
    integer(c_int), intent(in)    :: exch_id         ! Exchange identifier

    integer(c_int)                :: l_local_code_name, l_cpl_id, l_part_data_id

    l_local_code_name = len(local_code_name)
    l_cpl_id          = len(cpl_id)
    l_part_data_id    = len(part_data_id)

    call CWP_Part_data_wait_issend_cf(local_code_name,   &
                                      l_local_code_name, &
                                      cpl_id,            &
                                      l_cpl_id,          &
                                      part_data_id,      &
                                      l_part_data_id,    &
                                      exch_id)

  end subroutine CWP_Part_data_wait_issend_


  subroutine CWP_Part_data_wait_irecv_(local_code_name, &
                                       cpl_id,          &
                                       part_data_id,    &
                                       exch_id)
    ! Finalize the reception of a data array.
    use, intrinsic :: iso_c_binding
    implicit none
    character(kind=c_char, len=*) :: local_code_name ! Local code name
    character(kind=c_char, len=*) :: cpl_id          ! Coupling identifier
    character(kind=c_char, len=*) :: part_data_id    ! PartData identifier
    integer(c_int), intent(in)    :: exch_id         ! Exchange identifier

    integer(c_int)                :: l_local_code_name, l_cpl_id, l_part_data_id

    l_local_code_name = len(local_code_name)
    l_cpl_id          = len(cpl_id)
    l_part_data_id    = len(part_data_id)

    call CWP_Part_data_wait_irecv_cf(local_code_name,   &
                                     l_local_code_name, &
                                     cpl_id,            &
                                     l_cpl_id,          &
                                     part_data_id,      &
                                     l_part_data_id,    &
                                     exch_id)

  end subroutine CWP_Part_data_wait_irecv_


  function CWP_Cpl_spatial_interp_algo_get_(local_code_name, &
                                            cpl_id) result(algo)
    ! Get the coupling spatial interpolation algorithm.
    use, intrinsic :: iso_c_binding
    implicit none
    character(kind=c_char, len=*) :: local_code_name ! Local code name
    character(kind=c_char, len=*) :: cpl_id          ! Coupling identifier
    integer(c_int)                :: algo            ! Spatial interpolation method
    integer(c_int)                :: l_local_code_name, l_cpl_id

    l_local_code_name = len(local_code_name)
    l_cpl_id          = len(cpl_id)

    algo = CWP_Cpl_spatial_interp_algo_get_cf(local_code_name,   &
                                              l_local_code_name, &
                                              cpl_id,            &
                                              l_cpl_id)

  end function CWP_Cpl_spatial_interp_algo_get_


end module cwp
