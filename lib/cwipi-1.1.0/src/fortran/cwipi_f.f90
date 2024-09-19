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

module cwipi
  implicit none

  integer, parameter :: cwipi_int_l = 4

  !
  ! Parameters
  ! ----------
  !
  ! cwipi_nature_t
  ! TODO:  cwipi_nature a supprimer ?
  integer (kind = cwipi_int_l), parameter :: cwipi_nature_element_center = 0
  integer (kind = cwipi_int_l), parameter :: cwipi_nature_node = 1
  !
  ! cwipi_coupling_type_t
  integer (kind = cwipi_int_l), parameter :: cwipi_cpl_parallel_with_part = 0
  integer (kind = cwipi_int_l), parameter :: cwipi_cpl_parallel_without_part = 1
  integer (kind = cwipi_int_l), parameter :: cwipi_cpl_sequential = 2
  !
  ! cwipi_type_t
  integer (kind = cwipi_int_l), parameter :: cwipi_type_float = 0
  integer (kind = cwipi_int_l), parameter :: cwipi_type_double = 1
  !
  ! cwipi_interpolation_t
  integer (kind = cwipi_int_l), parameter :: cwipi_interpolation_standard = 0
  integer (kind = cwipi_int_l), parameter :: cwipi_interpolation_user = 1
  !
  ! mesh type
  integer (kind = cwipi_int_l), parameter :: cwipi_static_mesh = 0
  integer (kind = cwipi_int_l), parameter :: cwipi_mobile_mesh = 1
  integer (kind = cwipi_int_l), parameter :: cwipi_cyclic_mesh = 2
  !
  ! solver type
  integer (kind = cwipi_int_l), parameter :: cwipi_solver_cell_center = 0
  integer (kind = cwipi_int_l), parameter :: cwipi_solver_cell_vertex = 1
  !
  ! exchange status
  integer (kind = cwipi_int_l), parameter :: cwipi_exchange_ok = 0
  integer (kind = cwipi_int_l), parameter :: cwipi_exchange_bad_receiving = 1


  integer (kind = cwipi_int_l), parameter :: cwipi_node = 0
  integer (kind = cwipi_int_l), parameter :: cwipi_edge2 = 1
  integer (kind = cwipi_int_l), parameter :: cwipi_edgeho = 2
  integer (kind = cwipi_int_l), parameter :: cwipi_face_tria3 = 3
  integer (kind = cwipi_int_l), parameter :: cwipi_face_triaho = 4
  integer (kind = cwipi_int_l), parameter :: cwipi_face_quad4 = 5
  integer (kind = cwipi_int_l), parameter :: cwipi_face_quadho = 6
  integer (kind = cwipi_int_l), parameter :: cwipi_face_poly = 7
  integer (kind = cwipi_int_l), parameter :: cwipi_cell_tetra4 = 8
  integer (kind = cwipi_int_l), parameter :: cwipi_cell_tetraho = 9
  integer (kind = cwipi_int_l), parameter :: cwipi_cell_hexa8 = 10
  integer (kind = cwipi_int_l), parameter :: cwipi_cell_hexaho = 11
  integer (kind = cwipi_int_l), parameter :: cwipi_cell_prism6 = 12
  integer (kind = cwipi_int_l), parameter :: cwipi_cell_prismho = 13
  integer (kind = cwipi_int_l), parameter :: cwipi_cell_pyram5 = 14
  integer (kind = cwipi_int_l), parameter :: cwipi_cell_pyramho = 15
  integer (kind = cwipi_int_l), parameter :: cwipi_cell_poly = 16

  !
  ! Public interfaces
  interface cwipi_exchange_f ; module procedure &
    cwipi_exch_without_user_itp_f_, &
    cwipi_exch_with_user_itp_f_, &
    cwipi_exch_with_user_ho_itp_f_
  end interface cwipi_exchange_f

  interface cwipi_send_f     ; module procedure  &
    cwipi_send_without_user_itp_f_, &
    cwipi_send_with_user_itp_f_,&
    cwipi_send_with_user_ho_itp_f_
  end interface

  interface cwipi_issend_f     ; module procedure  &
    cwipi_issend_without_user_itp_f_, &
    cwipi_issend_with_user_itp_f_, &
    cwipi_issend_with_user_ho_itp_f_
  end interface

  interface cwipi_wait_issend_f     ; module procedure  &
    cwipi_wait_issend_f_
  end interface

  interface cwipi_wait_irecv_f     ; module procedure  &
    cwipi_wait_irecv_f_
  end interface

  interface cwipi_init_f                  ; module procedure &
    cwipi_init_f_
  end interface

  interface cwipi_add_loc_int_ctrl_param_f ; module procedure &
    cwipi_add_loc_int_ctrl_param_f_
  end interface

  interface cwipi_add_loc_dbl_ctrl_param_f ; module procedure &
    cwipi_add_loc_dbl_ctrl_param_f_
  end interface

  interface cwipi_add_loc_str_ctrl_param_f ; module procedure &
    cwipi_add_loc_str_ctrl_param_f_
  end interface

  interface cwipi_set_loc_int_ctrl_param_f ; module procedure &
    cwipi_set_loc_int_ctrl_param_f_
  end interface

  interface cwipi_set_loc_dbl_ctrl_param_f ; module procedure &
    cwipi_set_loc_dbl_ctrl_param_f_
  end interface

  interface cwipi_set_loc_str_ctrl_param_f ; module procedure &
    cwipi_set_loc_str_ctrl_param_f_
  end interface

  interface cwipi_get_loc_int_ctrl_param_f ; module procedure &
    cwipi_get_loc_int_ctrl_param_f_
  end interface

  interface cwipi_get_loc_dbl_ctrl_param_f ; module procedure &
    cwipi_get_loc_dbl_ctrl_param_f_
  end interface

  interface cwipi_get_loc_str_ctrl_param_f ; module procedure &
    cwipi_get_loc_str_ctrl_param_f_
  end interface

  interface cwipi_del_loc_int_ctrl_param_f ; module procedure &
    cwipi_del_loc_int_ctrl_param_f_
  end interface

  interface cwipi_del_loc_dbl_ctrl_param_f ; module procedure &
    cwipi_del_loc_dbl_ctrl_param_f_
  end interface

  interface cwipi_del_loc_str_ctrl_param_f ; module procedure &
    cwipi_del_loc_str_ctrl_param_f_
  end interface

  interface cwipi_get_dis_int_ctrl_param_f ; module procedure &
    cwipi_get_dis_int_ctrl_param_f_
  end interface

  interface cwipi_get_dis_dbl_ctrl_param_f ; module procedure &
    cwipi_get_dis_dbl_ctrl_param_f_
  end interface

  interface cwipi_get_dis_str_ctrl_param_f ; module procedure &
    cwipi_get_dis_str_ctrl_param_f_
  end interface

  interface cwipi_get_n_located_dist_pts_f ; module procedure &
    cwipi_get_n_located_dist_pts_f_
  end interface

  interface cwipi_synch_ctrl_param_f       ; module procedure &
    cwipi_synch_ctrl_param_f_
  end interface

  interface cwipi_create_coupling_f        ; module procedure &
   cwipi_create_coupling_f_
  end interface

  interface cwipi_set_points_to_locate_f   ; module procedure &
   cwipi_set_points_to_locate_f_
  end interface

  interface cwipi_define_mesh_f            ; module procedure &
   cwipi_define_mesh_f_
  end interface

  interface cwipi_ho_define_mesh_f ; module procedure &
   cwipi_ho_define_mesh_f_
  end interface

  interface  cwipi_ho_options_set_f ; module procedure &
    cwipi_ho_options_set_f_
  end interface
  
  interface cwipi_ho_ordering_from_IJK_set_f ; module procedure &
   cwipi_ho_ordering_from_IJK_set_f_
  end interface

  interface cwipi_ho_ordering_from_ref_elt_set_f ; module procedure &
   cwipi_ho_ordering_from_ref_elt_set_f_
  end interface

  interface cwipi_add_polyhedra_f          ; module procedure &
    cwipi_add_polyhedra_f_
  end interface

  interface cwipi_locate_f                 ; module procedure &
    cwipi_locate_f_
  end interface

  interface cwipi_update_location_f                 ; module procedure &
    cwipi_update_location_f_
  end interface

  interface cwipi_get_bary_coord_f         ; module procedure &
    cwipi_get_bary_coord_f_
  end interface

  interface cwipi_get_bary_coord_idx_f     ; module procedure &
    cwipi_get_bary_coord_idx_f_
  end interface

  interface cwipi_get_location_f           ; module procedure &
    cwipi_get_location_f_
  end interface

  interface cwipi_get_distance_f           ; module procedure &
    cwipi_get_distance_f_
  end interface

  interface cwipi_get_coord_f           ; module procedure &
    cwipi_get_coord_f_
  end interface

  interface cwipi_receive_f                ; module procedure &
    cwipi_receive_f_
  end interface

  interface cwipi_ireceive_f                ; module procedure &
    cwipi_ireceive_f_
  end interface

  interface cwipi_delete_coupling_f        ; module procedure &
    cwipi_delete_coupling_f_
  end interface

  interface cwipi_get_not_located_pts_f ; module procedure &
    cwipi_get_not_located_pts_f_
  end interface

  interface cwipi_get_n_not_located_pts_f ; module procedure &
    cwipi_get_n_not_located_pts_f_
  end interface

  interface cwipi_get_n_located_pts_f ; module procedure &
    cwipi_get_n_located_pts_f_
  end interface

  interface cwipi_get_located_pts_f ; module procedure &
    cwipi_get_located_pts_f_
  end interface

  interface  cwipi_dump_appli_properties_f ; module procedure &
      cwipi_dump_appli_properties_f_
  end interface

  interface  cwipi_finalize_f; module procedure &
    cwipi_finalize_f_
  end interface

  interface cwipi_dist_located_pts_get_f; module procedure &
    cwipi_dist_located_pts_get_f_  
  end interface cwipi_dist_located_pts_get_f

  interface cwipi_get_n_dis_ranks_f; module procedure &
    cwipi_get_n_dis_ranks_f_
  end interface cwipi_get_n_dis_ranks_f

  interface cwipi_get_dis_distrib_f; module procedure &
    cwipi_get_dis_distrib_f_
  end interface cwipi_get_dis_distrib_f

  interface cwipi_get_loc_pts_distrib_f; module procedure &
    cwipi_get_loc_pts_distrib_f_
  end interface cwipi_get_loc_pts_distrib_f

  interface cwipi_set_location_index_f; module procedure &
    cwipi_set_location_index_f_
  end interface cwipi_set_location_index_f

  interface cwipi_load_location_f; module procedure &
    cwipi_load_location_f_
  end interface cwipi_load_location_f

  interface cwipi_save_location_f; module procedure &
    cwipi_save_location_f_
  end interface cwipi_save_location_f

  interface cwipi_open_location_file_f; module procedure &
    cwipi_open_location_file_f_
  end interface cwipi_open_location_file_f

  interface cwipi_close_location_file_f; module procedure &
    cwipi_close_location_file_f_
  end interface cwipi_close_location_file_f

  interface cwipi_has_int_ctrl_param_f; module procedure &
    cwipi_has_int_ctrl_param_f_
  end interface cwipi_has_int_ctrl_param_f

  interface cwipi_has_dbl_ctrl_param_f; module procedure &
    cwipi_has_dbl_ctrl_param_f_
  end interface cwipi_has_dbl_ctrl_param_f

  interface cwipi_has_str_ctrl_param_f; module procedure &
    cwipi_has_str_ctrl_param_f_
  end interface cwipi_has_str_ctrl_param_f

  interface cwipi_get_n_int_ctrl_param_f; module procedure &
    cwipi_get_n_int_ctrl_param_f_
  end interface cwipi_get_n_int_ctrl_param_f

  interface cwipi_get_n_dbl_ctrl_param_f; module procedure &
    cwipi_get_n_dbl_ctrl_param_f_
  end interface  cwipi_get_n_dbl_ctrl_param_f

  interface cwipi_get_n_str_ctrl_param_f; module procedure &
    cwipi_get_n_str_ctrl_param_f_
  end interface cwipi_get_n_str_ctrl_param_f

  interface cwipi_get_list_int_ctrl_param_f; module procedure &
    cwipi_get_list_int_ctrl_param_f_
  end interface cwipi_get_list_int_ctrl_param_f

  interface cwipi_get_list_dbl_ctrl_param_f; module procedure &
    cwipi_get_list_dbl_ctrl_param_f_
  end interface cwipi_get_list_dbl_ctrl_param_f

  interface cwipi_get_list_str_ctrl_param_f; module procedure &
       cwipi_get_list_str_ctrl_param_f_
  end interface cwipi_get_list_str_ctrl_param_f
  
  abstract interface
    subroutine cwipi_ho_location_fct_t (entities_dim, &
                                        order, &
                                        n_nodes, &
                                        nodes_coords, &
                                        point_coords, &
                                        projected_coords,&
                                        projected_uvw) bind(c)
      use, intrinsic :: ISO_C_BINDING 
      implicit none
      integer (C_INT), value :: entities_dim
      integer (C_INT), value :: order
      integer (C_INT), value :: n_nodes

      type (C_PTR),    value  :: nodes_coords
      type (C_PTR),    value  :: point_coords
      type (C_PTR),    value  :: projected_coords
      type (C_PTR),    value  :: projected_uvw

    end subroutine cwipi_ho_location_fct_t

    subroutine  cwipi_ho_basis_fct_t (entities_dim, &
                                      order, &
                                      n_nodes, &
                                      n_pts, &
                                      uvw, &
                                      weights) bind(c)
      use, intrinsic :: ISO_C_BINDING 
      implicit none
      integer (C_INT), value :: entities_dim
      integer (C_INT), value :: order
      integer (C_INT), value :: n_nodes
      integer (C_INT), value :: n_pts

      type (C_PTR),    value :: uvw
      type (C_PTR),    value :: weights

    end subroutine cwipi_ho_basis_fct_t
  end interface

  interface cwipi_ho_user_elt_set_f; module procedure &
       cwipi_ho_user_elt_set_f_
  end interface cwipi_ho_user_elt_set_f

  !
  ! Private

  private :: cwipi_init_f_,                   &
             cwipi_exch_without_user_itp_f_,  &
             cwipi_exch_with_user_itp_f_,     &
             cwipi_send_without_user_itp_f_,  &
             cwipi_send_with_user_itp_f_,     &
             cwipi_add_loc_int_ctrl_param_f_, &
             cwipi_add_loc_dbl_ctrl_param_f_, &
             cwipi_add_loc_str_ctrl_param_f_, &
             cwipi_set_loc_int_ctrl_param_f_, &
             cwipi_set_loc_dbl_ctrl_param_f_, &
             cwipi_set_loc_str_ctrl_param_f_, &
             cwipi_get_loc_int_ctrl_param_f_, &
             cwipi_get_loc_dbl_ctrl_param_f_, &
             cwipi_get_loc_str_ctrl_param_f_, &
             cwipi_del_loc_int_ctrl_param_f_, &
             cwipi_del_loc_dbl_ctrl_param_f_, &
             cwipi_del_loc_str_ctrl_param_f_, &
             cwipi_get_dis_int_ctrl_param_f_, &
             cwipi_get_dis_dbl_ctrl_param_f_, &
             cwipi_get_dis_str_ctrl_param_f_, &
             cwipi_has_int_ctrl_param_f_, &
             cwipi_has_dbl_ctrl_param_f_, &
             cwipi_has_str_ctrl_param_f_, &
             cwipi_get_n_int_ctrl_param_f_, &
             cwipi_get_n_dbl_ctrl_param_f_, &
             cwipi_get_n_str_ctrl_param_f_, &
             cwipi_get_list_int_ctrl_param_f_, &
             cwipi_get_list_dbl_ctrl_param_f_, &
             cwipi_get_list_str_ctrl_param_f_, &
             cwipi_synch_ctrl_param_f_,       &
             cwipi_create_coupling_f_,        &
             cwipi_set_points_to_locate_f_,   &
             cwipi_ho_define_mesh_f_,          &
             cwipi_ho_options_set_f_,       &
             cwipi_ho_ordering_from_IJK_set_f_,&
             cwipi_ho_ordering_from_ref_elt_set_f_,&
             cwipi_ho_user_elt_set_f_, &
             cwipi_add_polyhedra_f_,          &
             cwipi_locate_f_,                 &
             cwipi_update_location_f_,        &
             cwipi_receive_f_,                &
             cwipi_delete_coupling_f_,        &
             cwipi_dump_appli_properties_f_,  &
             cwipi_finalize_f_,               &
             cwipi_get_location_f_,           &
             cwipi_get_distance_f_,           &
             cwipi_get_bary_coord_f_,         &
             cwipi_get_bary_coord_idx_f_,     &
             cwipi_get_coord_f_,              &
             cwipi_get_not_located_pts_f_,    &
             cwipi_dist_located_pts_get_f_,   &
             cwipi_get_n_dis_ranks_f_,        &
             cwipi_get_dis_distrib_f_,        &
             cwipi_get_loc_pts_distrib_f_,    &
             cwipi_set_location_index_f_,     &
             cwipi_load_location_f_,          &
             cwipi_save_location_f_,          &
             cwipi_open_location_file_f_,     &
             cwipi_close_location_file_f_      

contains

 subroutine cwipi_finalize_f_

    implicit none

    call cwipi_finalize_cf

  end subroutine cwipi_finalize_f_


  subroutine cwipi_dump_appli_properties_f_

    implicit none

    call cwipi_dump_appli_properties_cf

  end subroutine cwipi_dump_appli_properties_f_


!
!*******************************************************************************
!
! cwipi_init_f_
!
!  Initialize the cwipi library.
!  Redirect outputs in a file (Standard output with output_listing = NULL or
!  output_logical_unit = -1)
!  Create the current communicator application from 'common_comm'.
!
!  parameters:
!    globalComm    <-- Common MPI communicator
!    appliName     <-- Current application name
!    appliComm     --> Internal MPI communicator for the current
!                      application
!
!  It is a synchronization point between all applications
!
!*******************************************************************************
!

  subroutine cwipi_init_f_ (globalComm, appliName, appliComm)

    implicit none

    integer (kind = cwipi_int_l) :: globalcomm, applicomm
    character (len = *) :: appliname

    integer (kind = cwipi_int_l) :: l1

    l1 = len(appliname)

    call cwipi_init_cf (globalcomm, appliname, l1, applicomm)

  end subroutine cwipi_init_f_


!
!********************************************************************************
! cwipi_add_loc_int_ctrl_param_f
!
! Add a integer control parameter
!
! parameters
!    name           <-- parameter name
!    initial_value  <-- initial value
!
!********************************************************************************
!

  subroutine cwipi_add_loc_int_ctrl_param_f_ (name, initialvalue)

    implicit none

    character (len = *) :: name
    integer (kind = cwipi_int_l) :: initialvalue

    integer (kind = cwipi_int_l) :: l

    l = len(name)

    call cwipi_add_loc_int_ctrl_param_cf (name, l, initialvalue)

  end subroutine cwipi_add_loc_int_ctrl_param_f_

!
!********************************************************************************
!
! cwipi_add_loc_dbl_ctrl_param_f
!
! Add a double control parameter
!
! parameters
!    name           <-- parameter name
!    initial_value  <-- initial value
!********************************************************************************
!

  subroutine cwipi_add_loc_dbl_ctrl_param_f_ (name, initialvalue)

    implicit none

    character (len = *) ::name
    double precision :: initialvalue

    integer (kind = cwipi_int_l) :: l

    l = len(name)

    call cwipi_add_loc_dbl_ctrl_param_cf (name, l, initialvalue)

  end subroutine cwipi_add_loc_dbl_ctrl_param_f_

!
!********************************************************************************
!
! cwipi_add_loc_str_ctrl_param_f
!
! Add a double control parameter
!
! parameters
!    name           <-- parameter name
!    initial_value  <-- initial value
!********************************************************************************
!

  subroutine cwipi_add_loc_str_ctrl_param_f_(name, initialvalue)

    implicit none

    character (len = *) :: name
    character (len = *) :: initialvalue

    integer (kind = cwipi_int_l) :: l1, l2

    l1 = len(name)
    l2 = len(initialvalue)

    call cwipi_add_loc_str_ctrl_param_cf(name, l1, initialvalue, l2)

  end subroutine cwipi_add_loc_str_ctrl_param_f_

!
!********************************************************************************
!
! cwipi_set_loc_int_ctrl_param_f
!
!  Set a integer control parameter
!
!  parameters
!     name           <-- parameter name
!     value          <-- value
!
!********************************************************************************
!

  subroutine cwipi_set_loc_int_ctrl_param_f_(name, initialvalue)

    implicit none

    character (len = *) ::name
    integer (kind = cwipi_int_l) :: initialvalue

    integer (kind = cwipi_int_l) :: l

    l = len(name)

    call cwipi_set_loc_int_ctrl_param_cf (name, l, initialvalue)

  end subroutine cwipi_set_loc_int_ctrl_param_f_

!
!********************************************************************************
!
! cwipi_set_loc_dbl_ctrl_param_f
!
! Set a double control parameter
!
! parameters
!    name           <-- parameter name
!    value          <-- value
!
!
!********************************************************************************
!

  subroutine cwipi_set_loc_dbl_ctrl_param_f_ (name, initialvalue)

    implicit none

    character (len = *) :: name
    double precision :: initialvalue

    integer (kind = cwipi_int_l) :: l

    l = len(name)

    call cwipi_set_loc_dbl_ctrl_param_cf (name, l, initialvalue)

  end subroutine cwipi_set_loc_dbl_ctrl_param_f_

!
!********************************************************************************
!
! cwipi_set_loc_str_ctrl_param_f
!
! Set a double control parameter
!
! parameters
!    name           <-- parameter name
!    value          <-- value
!
!
!********************************************************************************
!

  subroutine cwipi_set_loc_str_ctrl_param_f_ (name, initialvalue)

    implicit none

    character (len = *) :: name
    character (len = *) :: initialvalue

    integer (kind = cwipi_int_l) :: l1, l2

    l1 = len(name)
    l2 = len(initialvalue)

    call cwipi_set_loc_str_ctrl_param_cf (name, l1, initialvalue, l2)

  end subroutine cwipi_set_loc_str_ctrl_param_f_

!
!********************************************************************************
!
! cwipi_get_loc_int_ctrl_param_f
!
! Get a integer control parameter of the current application
!
! parameters
!    name           <-- parameter name
!
!********************************************************************************
!

  subroutine cwipi_get_loc_int_ctrl_param_f_ (name, value)

    implicit none

    character (len = *) :: name
    integer (kind = cwipi_int_l) ::value

    integer (kind = cwipi_int_l) :: l

    l = len(name)

    call cwipi_get_loc_int_ctrl_param_cf (name, l, value)

  end subroutine cwipi_get_loc_int_ctrl_param_f_

!
!********************************************************************************
!
! cwipi_get_loc_dbl_ctrl_param_f
!
! Get a double control parameter of the current application
!
! parameters
!    name           <-- parameter name
!
!********************************************************************************
!

  subroutine cwipi_get_loc_dbl_ctrl_param_f_ (name, value)

    implicit none

    character (len = *) :: name
    double precision :: value

    integer (kind = cwipi_int_l) :: l

    l = len(name)

    call cwipi_get_loc_dbl_ctrl_param_cf (name, l, value)

  end subroutine cwipi_get_loc_dbl_ctrl_param_f_

!
!********************************************************************************
!
! cwipi_get_loc_str_ctrl_param_f
!
! Get a double control parameter of the current application
!
! parameters
!    name           <-- parameter name
!
!********************************************************************************
!

  subroutine cwipi_get_loc_str_ctrl_param_f_(name, value_str_f)

    implicit none

    character (len = *) :: name
    character (len = *) :: value_str_f

    integer (kind = cwipi_int_l) :: l1, l2

    l1 = len(name)
    l2 = len(value_str_f)

    call cwipi_get_loc_str_ctrl_param_cf(name, l1, value_str_f, l2)

  end subroutine cwipi_get_loc_str_ctrl_param_f_

!
!********************************************************************************
!
! cwipi_del_loc_int_ctrl_param_f
!
! Delete a current application int parameter
!
! parameters
!    name           <-- parameter name
!
!********************************************************************************
!

  subroutine cwipi_del_loc_int_ctrl_param_f_(name)

    implicit none

    character (len = *) :: name
    integer (kind = cwipi_int_l) l

    l = len(name)

    call cwipi_del_loc_int_ctrl_param_cf(name, l)

  end subroutine cwipi_del_loc_int_ctrl_param_f_

!
!********************************************************************************
!
! cwipi_del_loc_dbl_ctrl_param_f
!
! Delete a current application double parameter
!
! parameters
!    name           <-- parameter name
!
!********************************************************************************
!

  subroutine cwipi_del_loc_dbl_ctrl_param_f_ (name)

    implicit none

    character (len = *) :: name
    integer (kind = cwipi_int_l) :: l

    l = len(name)

    call cwipi_del_loc_dbl_ctrl_param_cf (name, l)

  end subroutine cwipi_del_loc_dbl_ctrl_param_f_

!
!********************************************************************************
!
! cwipi_del_loc_str_ctrl_param_f
!
! Delete a current application double parameter
!
! parameters
!    name           <-- parameter name
!
!********************************************************************************
!

  subroutine cwipi_del_loc_str_ctrl_param_f_ (name)

    implicit none

    character (len = *) :: name
    integer (kind = cwipi_int_l) :: l

    l = len(name)

    call cwipi_del_loc_str_ctrl_param_cf (name, l)

  end subroutine cwipi_del_loc_str_ctrl_param_f_

!
!********************************************************************************
!
! cwipi_has_int_ctrl_param_f
!
! Has int control parameter ?
!
! parameters
!    application_name       <-- application name
!    name                   <-- parameter name
!
!********************************************************************************
!

  subroutine cwipi_has_int_ctrl_param_f_(appliName, &
                                         paramName, &
                                         status)

    implicit none

    character (len = *) :: appliName
    character (len = *) :: paramName
    integer (kind = cwipi_int_l) :: status

    integer (kind = cwipi_int_l) :: l1, l2

    l1 = len(appliName)
    l2 = len(paramName)

    call cwipi_has_int_ctrl_param_cf(appliName, &
                                         l1, &
                                         paramName, &
                                         l2, &
                                         status)

  end subroutine cwipi_has_int_ctrl_param_f_

!
!********************************************************************************
!
! cwipi_has_dbl_ctrl_param_f
!
! has double control parameter ?
!
! parameters
!    application_name       <-- application name
!    name                   <-- parameter name
!
!********************************************************************************
!

  subroutine cwipi_has_dbl_ctrl_param_f_(appliName, &
                                             paramName, &
                                             status)

    implicit none

    character (len = *) :: appliName
    character (len = *) :: paramName
    integer (kind = cwipi_int_l) :: status

    integer (kind = cwipi_int_l) :: l1, l2

    l1 = len(appliName)
    l2 = len(paramName)

    call cwipi_has_dbl_ctrl_param_cf(appliName, &
                                     l1, &
                                     paramName, &
                                     l2, &
                                     status)

  end subroutine cwipi_has_dbl_ctrl_param_f_

!
!********************************************************************************
!
! cwipi_has_str_ctrl_param_f
!
! has double control parameter ?
!
! parameters
!    application_name       <-- application name
!    name                   <-- parameter name
!
!********************************************************************************
!

  subroutine cwipi_has_str_ctrl_param_f_(appliName, &
                                         paramName, &
                                         status)

    implicit none

    character (len = *) :: appliName
    character (len = *) :: paramName
    integer (kind = cwipi_int_l) :: status

    integer (kind = cwipi_int_l) :: l1, l2

    l1 = len(appliName)
    l2 = len(paramName)

    call cwipi_has_str_ctrl_param_cf(appliName, &
                                     l1, &
                                     paramName, &
                                     l2, &
                                     status)

  end subroutine cwipi_has_str_ctrl_param_f_

!
!********************************************************************************
!
! cwipi_get_n_int_ctrl_param_f
!
! Get number of int control parameter
!
! parameters
!    application_name       <-- application name
!    name                   <-- parameter name
!
!********************************************************************************
!

  subroutine cwipi_get_n_int_ctrl_param_f_(appliName, &
                                           n_param)

    implicit none

    character (len = *) :: appliName
    integer (kind = cwipi_int_l) :: n_param

    integer (kind = cwipi_int_l) :: l1

    l1 = len(appliName)

    call cwipi_get_n_int_ctrl_param_cf(appliName, &
                                     l1, &
                                     n_param)

  end subroutine cwipi_get_n_int_ctrl_param_f_


!
!********************************************************************************
!
! cwipi_get_n_dbl_ctrl_param_f
!
! Get number of int control parameter
!
! parameters
!    application_name       <-- application name
!    name                   <-- parameter name
!
!********************************************************************************
!

  subroutine cwipi_get_n_dbl_ctrl_param_f_(appliName, &
                                           n_param)

    implicit none

    character (len = *) :: appliName
    integer (kind = cwipi_int_l) :: n_param

    integer (kind = cwipi_int_l) :: l1

    l1 = len(appliName)

    call cwipi_get_n_dbl_ctrl_param_cf(appliName, &
                                       l1, &
                                       n_param)

  end subroutine cwipi_get_n_dbl_ctrl_param_f_

!
!********************************************************************************
!
! cwipi_get_n_str_ctrl_param_f
!
! Get number of int control parameter
!
! parameters
!    application_name       <-- application name
!    name                   <-- parameter name
!
!********************************************************************************
!

  subroutine cwipi_get_n_str_ctrl_param_f_(appliName, &
                                           n_param)

    implicit none

    character (len = *) :: appliName
    integer (kind = cwipi_int_l) :: n_param

    integer (kind = cwipi_int_l) :: l1

    l1 = len(appliName)

    call cwipi_get_n_str_ctrl_param_cf(appliName, &
                                       l1, &
                                       n_param)

  end subroutine cwipi_get_n_str_ctrl_param_f_


!
!********************************************************************************
!
! cwipi_get_list_int_ctrl_param_f
!
! Get int control parameter list
!
! parameters
!    application_name       <-- application name
!    name                   <-- parameter name
!
!********************************************************************************
!

  subroutine cwipi_get_list_int_ctrl_param_f_(appliName, &
                                              params)

    implicit none

    character (len = *) :: appliName
    character (len = *) :: params(*)

    integer (kind = cwipi_int_l) :: l1, l2

    l1 = len(appliName)
    l2 = len(params)

    call cwipi_get_list_int_ctrl_param_cf(appliName, &
                                          l1, &
                                          params, &
                                          l2)

  end subroutine cwipi_get_list_int_ctrl_param_f_


!
!********************************************************************************
!
! cwipi_get_list_dbl_ctrl_param_f
!
! Get double control parameter list
!
! parameters
!    application_name       <-- application name
!    name                   <-- parameter name
!
!********************************************************************************
!

  subroutine cwipi_get_list_dbl_ctrl_param_f_(appliName, &
                                              params)

    implicit none

    character (len = *) :: appliName
    character (len = *) :: params(*)

    integer (kind = cwipi_int_l) :: l1, l2

    l1 = len(appliName)
    l2 = len(params)

    call cwipi_get_list_dbl_ctrl_param_cf(appliName, &
                                          l1, &
                                          params, &
                                          l2)

  end subroutine cwipi_get_list_dbl_ctrl_param_f_



!
!********************************************************************************
!
! cwipi_get_list_dbl_ctrl_param_f
!
! Get double control parameter list
!
! parameters
!    application_name       <-- application name
!    name                   <-- parameter name
!
!********************************************************************************
!

  subroutine cwipi_get_list_str_ctrl_param_f_(appliName, &
                                              params)

    implicit none

    character (len = *) :: appliName
    character (len = *) :: params(*)

    integer (kind = cwipi_int_l) :: l1, l2

    l1 = len(appliName)
    l2 = len(params)

    call cwipi_get_list_str_ctrl_param_cf(appliName, &
                                          l1, &
                                          params, &
                                          l2)

  end subroutine cwipi_get_list_str_ctrl_param_f_


!
!********************************************************************************
!
! cwipi_get_dis_int_ctrl_param_f
!
! Get a integer control parameter of a other application
!
! parameters
!    application_name       <-- application name
!    name                   <-- parameter name
!
!********************************************************************************
!

  subroutine cwipi_get_dis_int_ctrl_param_f_ (appliName, &
                                                             paramName, &
                                                             value)

    implicit none

    character (len = *) :: appliName
    character (len = *) :: paramName
    integer (kind = cwipi_int_l) :: value

    integer (kind = cwipi_int_l) :: l1, l2

    l1 = len(appliName)
    l2 = len(paramName)

    call cwipi_get_dis_int_ctrl_param_cf (appliName, &
                                                         l1, &
                                                         paramName, &
                                                         l2, &
                                                         value)

  end subroutine cwipi_get_dis_int_ctrl_param_f_

!
!********************************************************************************
!
! cwipi_get_dis_dbl_ctrl_param_f
!
! Get a double control parameter of a other application
!
! parameters
!    application_name    <-- application name
!    name                <-- parameter name
!
!********************************************************************************
!

  subroutine cwipi_get_dis_dbl_ctrl_param_f_  (appliName, &
                                                               paramName, &
                                                               value)

    implicit none

    character (len = *) :: appliName
    character (len = *) :: paramName
    double precision :: value

    integer (kind = cwipi_int_l) l1, l2

    l1 = len(appliName)
    l2 = len(paramName)

    call cwipi_get_dis_dbl_ctrl_param_cf (appliName, &
                                                            l1, &
                                                            paramName, &
                                                            l2, &
                                                            value)

  end subroutine cwipi_get_dis_dbl_ctrl_param_f_


!
!********************************************************************************
!
! cwipi_get_dis_dbl_ctrl_param_f
!
! Get a double control parameter of a other application
!
! parameters
!    application_name    <-- application name
!    name                <-- parameter name
!
!********************************************************************************
!

  subroutine cwipi_get_dis_str_ctrl_param_f_(appliName, &
                                             paramName, &
                                             value_str_f)

    implicit none

    character (len = *) :: appliName
    character (len = *) :: paramName
    character (len = *) :: value_str_f

    integer (kind = cwipi_int_l) :: l1, l2, l3

    l1 = len(appliName)
    l2 = len(paramName)
    l3 = len(value_str_f)

    call cwipi_get_dis_str_ctrl_param_cf (appliName, &
                                          l1, &
                                          paramName, &
                                          l2, &
                                          value_str_f, &
                                          l3)

  end subroutine cwipi_get_dis_str_ctrl_param_f_

!
!********************************************************************************
!
! cwipi_synchronize_control_parameter_f
!
! Synchronize local control parameters with an other application.
!  It is a synchronization point with this second application
!
! parameters
!    appliName           <-- application name
!
!********************************************************************************
!

  subroutine cwipi_synch_ctrl_param_f_ (appliName)

    implicit none

    character (len = *) :: appliName

    integer (kind = cwipi_int_l) l

    l = len(appliName)

    call cwipi_synch_ctrl_param_cf (appliName, l)

  end subroutine cwipi_synch_ctrl_param_f_

!
!********************************************************************************
!
! cwipi_dump_application_properties_f (define into cwipi_cf.hxx)
!
! Dump application properties
!
!********************************************************************************
!

!
!********************************************************************************
!
! cwipi_create_coupling_f
!
! Create a coupling object
!
! parameters:
!   couplingName            <-- Coupling identifier
!   couplingType            <-- Coupling type
!   cplAppli                <-- Coupled application name
!   entitiesDim             <-- Mesh entities dimension (1, 2 or 3)
!   tolerance               <-- Geometric tolerance to locate
!   meshT                   <-- CWIPI_STATIC_MESH
!                               CWIPI_MOBILE_MESH (not implemented yet)
!                               CWIPI_CYCLIC_MESH 
!   solverT                 <-- CWIPI_SOLVER_CELL_CENTER
!                               CWIPI_SOLVER_CELL_VERTEX
!   outputFreq              <-- Output frequency
!   outputFmt               <-- Output format to visualize exchanged fields
!                               on the coupled mesh. Choice between :
!                                 - "EnSight Gold"
!                                 - "MED_fichier"
!                                 - "CGNS"
!   outputFmtOpt            <-- Output options
!                             text                output text files
!                             binary              output binary files (default)
!                             big_endian          force binary files
!                                                 to big-endian
!                             discard_polygons    do not output polygons
!                                                 or related values
!                             discard_polyhedra   do not output polyhedra
!                                                 or related values
!                             divide_polygons     tesselate polygons
!                                                 with triangles
!                             divide_polyhedra    tesselate polyhedra
!                                                 with tetrahedra and pyramids
!                                                 (adding a vertex near
!                                                 each polyhedron's center)
!   nbLocations             <-- maximun number of locations (only with 
!                               CWIPI_CYCLIC_MESH) optional default 1
!
!********************************************************************************
!

  subroutine cwipi_create_coupling_f_ (couplingName, &
                                           couplingType, &
                                           cplAppli, &
                                           entitiesDim, &
                                           tolerance, &
                                           meshT, &
                                           solvert, &
                                           outputfreq, &
                                           outputfmt, &
                                           outputfmtopt,&
                                           nbLocations)

    implicit none

    character (len = *) :: couplingName, cplAppli
    integer (kind = cwipi_int_l) :: entitiesDim, meshT, solverT, couplingType
    double precision :: tolerance
    integer (kind = cwipi_int_l) :: outputFreq
    character (len = *) :: outputFmt, outputFmtOpt

    integer (kind = cwipi_int_l) :: lCouplingName, lCplAppli
    integer (kind = cwipi_int_l) :: lOutputFmt, lOutputFmtOpt
    integer (kind = cwipi_int_l),intent(in), optional :: nbLocations
   
    integer (kind = cwipi_int_l) nbLoc

    if(present(nbLocations)) then
       nbLoc = nbLocations
    else
       nbLoc = 1
    endif

    lCouplingName = len(couplingName)
    lCplAppli     = len(cplAppli)
    lOutputFmt    = len(outputFmt)
    lOutputFmtOpt = len(outputFmtOpt)

    call cwipi_create_coupling_cf(couplingName, &
                                  lCouplingName, &
                                  couplingType, &
                                  cplAppli, &
                                  lCplAppli, &
                                  entitiesDim, &
                                  tolerance, &
                                  meshT, &
                                  solverT, &
                                  outputFreq, &
                                  outputFmt, &
                                  lOutputFmt, &
                                  outputFmtOpt, &
                                  lOutputFmtOpt, &
                                  nbLoc)
    
  end subroutine cwipi_create_coupling_f_

!
!********************************************************************************
!
! cwipi_set_points_to_locate_f
!
! Set points to locate. This function must be called if the points to locate
! do not correspond to :
!        - vertices for CELL_VERTEX nature
!        - cell center for CELL_CENTER nature
!
! parameters:
!   couplingName       <-- coupling identifier
!   nPts               <-- number of points to locate
!   coords             <-- coordinates of points to locate (enterlaced)
!
!********************************************************************************
!

  subroutine cwipi_set_points_to_locate_f_ (couplingName, &
                                                nPts, &
                                                coords)

    implicit none

    character (len = *) :: couplingName

    integer (kind = cwipi_int_l) :: nPts
    double precision, dimension(3 * npts) :: coords

    integer (kind = cwipi_int_l) :: lCouplingName

    lCouplingName  = len(couplingName)

    call cwipi_set_points_to_locate_cf(couplingName, &
                                           lCouplingName, &
                                           nPts, &
                                           coords)

  end subroutine cwipi_set_points_to_locate_f_

!
!********************************************************************************
!
! cwipi_define_mesh_f
!
!
! Define the support mesh for a coupling. The connectivity is sorted if
! necessary.
!
! Order definition :
!    1D : edges
!    2D : triangles, quadrangles, polygons
!    3D : tetrahedra, pyramids, prism, hexaedra
!
! Local connectivity for the following element type :
!
!  - edge :
!
!   1 x-------x 2
!
!  - triangle :
!
!   1 x-------x 3
!      \     /
!       \   /
!        \ /
!         x 2
!
!  - quadrangle :
!
!      4 x-------x 3
!       /       /
!      /       /
!   1 x-------x2
!
!   - tetrahedra :
!
!         x 4
!        /|\
!       / | \
!      /  |  \
!   1 x- -|- -x 3
!      \  |  /
!       \ | /
!        \|/
!         x 2
!
!   - pyramid :
!
!          5 x
!           /|\
!          //| \
!         // |  \
!      4 x/--|---x 3
!       //   |  /
!      //    | /
!   1 x-------x 2
!
!  - prism :
!
!   4 x-------x 6
!     |\     /|
!     | \   / |
!   1 x- \-/ -x 3
!      \ 5x  /
!       \ | /
!        \|/
!         x 2
!
!  -  hexahedra :
!
!      8 x-------x 7
!       /|      /|
!      / |     / |
!   5 x-------x6 |
!     | 4x----|--x 3
!     | /     | /
!     |/      |/
!   1 x-------x 2
!
!
!
! parameters:
!   couplingName       <-- coupling name
!   nVertex            <-- number of vertices
!   nElts              <-- number of elements
!   coords             <-- vertex interlaced coordinates
!   connecIndex        <-- element -> vertices index (O to n-1)
!                          size: n_elements + 1
!   connec             <-- element -> vertex connectivity (1 to n)
!                          size: connectivity_index[n_elements]
!
!********************************************************************************
!

  subroutine cwipi_define_mesh_f_ (couplingName, &
                                       nVertex, &
                                       nElts, &
                                       coords, &
                                       connecIndex, &
                                       connec)

    implicit none

    character (len = *) :: couplingName
    integer (kind = cwipi_int_l) :: lCouplingName

    integer (kind = cwipi_int_l) :: nElts, nVertex
    integer (kind = cwipi_int_l), dimension (nelts+1) :: connecIndex(nelts+1)

    integer (kind = cwipi_int_l), dimension (*) ::  connec
    double precision, dimension(3 * nVertex) :: coords

    lCouplingName    = len(couplingName)

    call cwipi_define_mesh_cf(couplingName, &
                                  lCouplingName, &
                                  nVertex, &
                                  nElts, &
                                  coords, &
                                  connecindex, &
                                  connec)

  end subroutine cwipi_define_mesh_f_

  

  subroutine cwipi_ho_define_mesh_f_ (couplingName, &
                                       nVertex, &
                                       nElts, &
                                       order, &
                                       coords, &
                                       connecIndex, &
                                       connec)

    implicit none

    character (len = *) :: couplingName
    integer (kind = cwipi_int_l) :: lCouplingName

    integer (kind = cwipi_int_l) :: nElts, nVertex, order
    integer (kind = cwipi_int_l), dimension (nelts+1) :: connecIndex(nelts+1)

    integer (kind = cwipi_int_l), dimension (*) ::  connec
    double precision, dimension(3 * nVertex) :: coords

    lCouplingName    = len(couplingName)

    call cwipi_ho_define_mesh_cf(couplingName, &
                                  lCouplingName, &
                                  nVertex, &
                                  nElts, &
                                  order, &
                                  coords, &
                                  connecindex, &
                                  connec)

  end subroutine cwipi_ho_define_mesh_f_


 !********************************************************************************
 !
 ! Define specific options for ho elements
 !
 ! parameters:
 !   coupling_id     <-- coupling name
 !   option          <--  option name, Choice between :
 !                          - "opt_bbox_step" 
 !                              * Description : step of discretization used 
 !                                              to compute the optimized element 
 !                                              bounding boxes
 !                                              -1 to deactivate this computation
 !                              * Default     : 10 
 !   opt_value       <-- option value
 !
 !********************************************************************************

  subroutine cwipi_ho_options_set_f_ (couplingName, &
                                      option, &
                                      opt_value)

    implicit none

    character (len = *) :: couplingName
    integer (kind = cwipi_int_l) :: lCouplingName

    character (len = *) :: option
    integer (kind = cwipi_int_l) :: loption

    character (len = *) :: opt_value
    integer (kind = cwipi_int_l) :: lopt_value

    lCouplingName = len(couplingName)
    loption       = len(option)
    lopt_value    = len(opt_value)

    call cwipi_ho_options_set_cf (couplingName, &
                                  lCouplingName, &
                                  option,&
                                  loption,&
                                  opt_value, &
                                  lopt_value)

  end subroutine cwipi_ho_options_set_f_
  

 !********************************************************************************
 !
 ! Define ho element ordering from the location in the (u, v, w) grid
 !
 ! parameters:
 !   coupling_id     <-- coupling name
 !   t_elt           <-- element type
 !   n_nodes         <-- number of nodes
 !   IJK             <-- user ordering to (I, J, K) grid (size = elt_dim * n_nodes)
 !
 !********************************************************************************

  subroutine cwipi_ho_ordering_from_IJK_set_f_ (couplingName, &
                                               tElt,       &
                                               nNodes,     &
                                               IJK)

    implicit none

    character (len = *) :: couplingName
    integer (kind = cwipi_int_l) :: lCouplingName
    integer (kind = cwipi_int_l) :: tElt, nNodes

    integer (kind = cwipi_int_l), dimension (*) :: IJK

    lCouplingName    = len(couplingName)

    call cwipi_ho_ordering_from_IJK_set_cf (couplingName, lCouplingName, tElt, nNodes, IJK)

  end subroutine cwipi_ho_ordering_from_IJK_set_f_
  

 !********************************************************************************
 !
 ! Define a user ho element 
 !
 ! parameters:
 !   t_elt               <-- element type
 !   element_basis       <-- element basis function
 !   location_in_element <-- location in element function
 !
 !********************************************************************************

  subroutine cwipi_ho_user_elt_set_f_ (elt_type, &
                                       element_basis, &
                                       location_in_element)

    use, intrinsic :: ISO_C_BINDING 

    implicit none

    interface
      subroutine cwipi_ho_user_elt_set_c (elt_type, &
                                          element_basis,&
                                          location_in_element) bind (c, name = "cwipi_ho_user_elt_set")
        use, intrinsic :: ISO_C_BINDING 
        implicit none
        
        integer (C_INT), value :: elt_type
        type (c_funptr), value :: element_basis
        type (c_funptr), value :: location_in_element

      end subroutine cwipi_ho_user_elt_set_c
    end interface
    integer (kind = cwipi_int_l) :: elt_type
    procedure (cwipi_ho_basis_fct_t) :: element_basis
    procedure (cwipi_ho_location_fct_t) :: location_in_element

    integer (C_INT) :: elt_type_
    type (c_funptr) :: element_basis_c
    type (c_funptr) :: location_in_element_c

    elt_type_ = elt_type
    element_basis_c = c_funloc(element_basis)
    location_in_element_c = c_funloc(location_in_element)
    
    call cwipi_ho_user_elt_set_c (elt_type_, element_basis_c, location_in_element_c)

  end subroutine cwipi_ho_user_elt_set_f_
  
!********************************************************************************
!
! Define ho element ordering from reference element (definition between 0 - 1)
!
!   couplingId        <-- coupling name
!   tElt              <-- element type
!   nNodes            <-- number of nodes
!   coords            <-- node coordinates of reference element
!                                TODO: decrire ici les elements de reference
!
!********************************************************************************

 subroutine cwipi_ho_ordering_from_ref_elt_set_f_ (couplingName, &
                                                   tElt,       &
                                                   nNodes,     &
                                                   coords)

   implicit none

   character (len = *) :: couplingName
   integer (kind = cwipi_int_l) :: lCouplingName
   integer (kind = cwipi_int_l) :: tElt, nNodes

   double precision, dimension(*) :: coords

   lCouplingName    = len(couplingName)

   call cwipi_ho_ordering_from_ref_elt_set_cf (couplingName, lCouplingName, tElt, nNodes, coords)

 end subroutine cwipi_ho_ordering_from_ref_elt_set_f_

!
!********************************************************************************
!
! cwipi_add_polyhedra_f
!
! parameters:
!   couplingName       <-- Coupling name
!   nElts              <-- Number of elements
!   cellToFaceIdx      <-- Cell -> faces connectivity index (0 to n-1)
!                          size: nElts + 1
!   cellToFace         <-- Cell -> faces connectivity (1 to n)                         
!                          size: cellToFaceIdx(nElts)
!   nFaces             <-- Number of faces
!   faceConnecIdx      <-- Face -> vertex connectivity index (O to n-1)
!                          size: nFaces + 1
!   faceConnec         <-- Face -> vertex connectivity (1 to n)
!                          size: faceConnecIdx[nFaces]
!*******************************************************************************
!

  subroutine cwipi_add_polyhedra_f_ (couplingName, &
                                     nElts, &
                                     cellToFaceIdx, &
                                     cellToFace, &
                                     nFaces, &
                                     faceConnecIdx, &
                                     faceConnec)

    implicit none

    character (len = *) :: couplingName
    integer (kind = cwipi_int_l) :: lCouplingname, nElts, nFaces
    integer (kind = cwipi_int_l), dimension(nelts) :: cellToFaceIdx
    integer (kind = cwipi_int_l), dimension(*) :: cellToFace, faceConnecIdx, faceConnec

    lCouplingName = len(couplingName)

    call cwipi_add_polyhedra_cf (couplingName, &
                                 lCouplingName, &
                                 nElts, &
                                 cellToFaceIdx, &
                                 cellToFace, &
                                 nFaces, &
                                 faceConnecIdx, &
                                 faceConnec)

  end subroutine cwipi_add_polyhedra_f_

!
!********************************************************************************
!
! cwipi_locate_f
!
! Location completion.
! It is a synchronization point with the coupled application
!
! parameters
!   couplingName         <-- Coupling identifier
!
!*******************************************************************************
!

  subroutine cwipi_locate_f_(couplingName)

    implicit none

    character (len = *) :: couplingName
    integer (kind = cwipi_int_l) :: lcouplingname

    lCouplingName = len(couplingName)

    call cwipi_locate_cf(couplingName, lCouplingName)
  end subroutine cwipi_locate_f_

!
!********************************************************************************
!
! cwipi_dist_located_pts_get_f_
!
! Return located points to interface
!
! parameters
!   couplingName         <-- Coupling identifier
!
!*******************************************************************************
!

  subroutine cwipi_dist_located_pts_get_f_(couplingName, distance)
    implicit none

    character (len = *) :: couplingName
    real(kind=4) :: distance(*)
    integer (kind = cwipi_int_l) :: lcouplingname

    lCouplingName = len(couplingName)

    call cwipi_dist_located_pts_get_cf(couplingName, lCouplingName, distance)
  end subroutine cwipi_dist_located_pts_get_f_
!
!********************************************************************************
!
! cwipi_update_location_f
!
! Location completion.
! It is a synchronization point with the coupled application
!
! parameters
!   couplingName         <-- Coupling identifier
!
!*******************************************************************************
!

  subroutine cwipi_update_location_f_(couplingName)

    implicit none

    character (len = *) :: couplingName
    integer (kind = cwipi_int_l) :: lcouplingname

    lCouplingName = len(couplingName)

    call cwipi_update_location_cf(couplingName, lCouplingName)
  end subroutine cwipi_update_location_f_


!
!********************************************************************************
!
! cwipi_set_location_index_f
!
! parameters
!   couplingName          <-- Coupling identifier
!   index                 <-- location index
!
!*******************************************************************************
!

  subroutine cwipi_set_location_index_f_(couplingName, index)

    implicit none

    character (len = *) :: couplingName
    integer (kind = cwipi_int_l) :: lcouplingname
    integer (kind = cwipi_int_l) :: index

    lCouplingName = len(couplingName)

    call cwipi_set_location_index_cf(couplingName, lCouplingName, index)
  end subroutine cwipi_set_location_index_f_


!
!********************************************************************************
!
! cwipi_load_location_f
!
! parameters
!   couplingName          <-- Coupling identifier
!
!*******************************************************************************
!

  subroutine cwipi_load_location_f_(couplingName)

    implicit none

    character (len = *) :: couplingName
    integer (kind = cwipi_int_l) :: lcouplingname

    lCouplingName = len(couplingName)

    call cwipi_load_location_cf(couplingName, lCouplingName)
  end subroutine cwipi_load_location_f_


!
!********************************************************************************
!
! cwipi_save_location_f
!
! parameters
!   couplingName          <-- Coupling identifier
!
!*******************************************************************************
!

  subroutine cwipi_save_location_f_(couplingName)

    implicit none

    character (len = *) :: couplingName
    integer (kind = cwipi_int_l) :: lcouplingname

    lCouplingName = len(couplingName)

    call cwipi_save_location_cf(couplingName, lCouplingName)
  end subroutine cwipi_save_location_f_


!
!********************************************************************************
!
! cwipi_save_location_f
!
! parameters
!   couplingName          <-- Coupling identifier
!
!*******************************************************************************
!

  subroutine cwipi_open_location_file_f_(couplingName, fileName, mode)

    implicit none

    character (len = *) :: couplingName
    character (len = *) :: fileName
    character (len = *) :: mode
    integer (kind = cwipi_int_l) :: lcouplingname
    integer (kind = cwipi_int_l) :: lfileName
    integer (kind = cwipi_int_l) :: lmode

    lCouplingName = len(couplingName)
    lfileName = len(fileName)
    lmode = len(mode)

    call cwipi_open_location_file_cf(couplingName, lCouplingName, &
                                     fileName, lfileName, &
                                     mode, lmode)
 
  end subroutine cwipi_open_location_file_f_



!
!********************************************************************************
!
! cwipi_save_location_f
!
! parameters
!   couplingName          <-- Coupling identifier
!
!*******************************************************************************
!

  subroutine cwipi_close_location_file_f_(couplingName)

    implicit none

    character (len = *) :: couplingName
    integer (kind = cwipi_int_l) :: lcouplingname

    lCouplingName = len(couplingName)

    call cwipi_close_location_file_cf(couplingName, lCouplingName)
 
  end subroutine cwipi_close_location_file_f_

!
!********************************************************************************
!
! cwipi_get_location_f
!
! Get located points location
!
! parameters
!   couplingName         <-- Coupling identifier
!   location             <-- Get located points location
!
!*******************************************************************************
!

  subroutine cwipi_get_location_f_(couplingName, location)

    implicit none

    character (len = *) :: couplingName
    integer (kind = cwipi_int_l) :: lcouplingname
    integer (kind = cwipi_int_l), dimension(*) :: location

    lCouplingName = len(couplingName)

    call cwipi_get_distant_location_cf (couplingName, lCouplingName, location)
  end subroutine cwipi_get_location_f_

!
!********************************************************************************
!
! cwipi_get_location_f
!
! Get located points location
!
! parameters
!   couplingName         <-- Coupling identifier
!   location             <-- Get located points location
!
!*******************************************************************************
!

  subroutine cwipi_get_distance_f_(couplingName, distance)

    implicit none

    character (len = *) :: couplingName
    integer (kind = cwipi_int_l) :: lcouplingname
    real(kind=4), dimension(*) :: distance

    lCouplingName = len(couplingName)

    call cwipi_get_distant_distance_cf (couplingName, lCouplingName, distance)
  end subroutine cwipi_get_distance_f_

!
!********************************************************************************
!
! cwipi_get_coord
!
! Get distant points coordinates
!
! parameters
!   couplingName      <-- Coupling identifier
!   coordinates       <-- Get distant point coordinates
! 
!*******************************************************************************
!

  subroutine cwipi_get_coord_f_(couplingName, &
                                coordinates)

    implicit none

    character (len = *) :: couplingName
    integer (kind = cwipi_int_l) :: lcouplingname
    double precision, dimension(*) :: coordinates

    lCouplingName = len(couplingName)

    call cwipi_get_dis_coord_cf(couplingName, &
                                lCouplingName, &
                                coordinates)
  end subroutine cwipi_get_coord_f_

!
!********************************************************************************
!
! cwipi_get_bary_coord_idx_f
!
! Get located points barycentric coordinates index
!
! parameters
!   couplingName                 <-- Coupling identifier
!   barycentricCoordinatesIndex  <-- Get located points barycentric coordinates
!                                    index
!*******************************************************************************
!

  subroutine cwipi_get_bary_coord_idx_f_(couplingName, &
                                         barycentricCoordinatesIndex)

    implicit none

    character (len = *) :: couplingName
    integer (kind = cwipi_int_l) :: lcouplingname
    integer (kind = cwipi_int_l), dimension(*) :: barycentricCoordinatesIndex

    lCouplingName = len(couplingName)

    call cwipi_get_dis_bary_coord_idx_cf (couplingName, &
                                          lCouplingName, &
                                          barycentricCoordinatesIndex)
  end subroutine cwipi_get_bary_coord_idx_f_

!
!********************************************************************************
!
! cwipi_get_bary_coord_f
!
! Get located points barycentric coordinates
!
! parameters
!   couplingName              <-- Coupling identifier
!   barycentricCoordinates   <-- Get located points barycentric coordinates
!
!*******************************************************************************
!

  subroutine cwipi_get_bary_coord_f_(couplingName, &
                                     barycentricCoordinates)

    implicit none

    character (len = *) :: couplingName
    integer (kind = cwipi_int_l) :: lcouplingname
    double precision, dimension(*) :: barycentricCoordinates

    lCouplingName = len(couplingName)

    call cwipi_get_dis_bary_coord_cf (couplingName, &
                                      lCouplingName, &
                                      barycentricCoordinates)
  end subroutine cwipi_get_bary_coord_f_

!
!********************************************************************************
!
! cwipi_exch_without_user_itp_f_
!
! Exchange data with the coupled application.
! It is a synchronization point with the coupled application
!
! parameters
!   couplingName         <-- Coupling identifier
!   exchangeName         <-- Exchange name
!   stride               <-- Number of interlaced fields
!   timeStep             <-- Time step  (only for visualization)
!   timeValue            <-- Time value (only for visualization)
!   sendingFieldName     <-- Sending field name
!   sendingField         <-- Sending field (NULL -> no sending)
!   receivingFieldName   <-- Receiving field name
!   receivingField       --> Receiving field
!   nNotLocatedPoints    --> Number of not located points
!   status               --> Cwipi exchange status
!
!********************************************************************************
!

  subroutine cwipi_exch_without_user_itp_f_ (couplingName, &
                                             exchangeName, &
                                             stride, &
                                             nStep, &
                                             timeValue, &
                                             sendingFieldName, &
                                             sendingField, &
                                             receivingFieldName, &
                                             receivingField, &
                                             nNotLocatedPoints, &
                                             status)

    implicit none

    character (len = *) :: couplingName, exchangeName, sendingFieldName
    character (len = *) :: receivingFieldName
    integer (kind = cwipi_int_l) :: stride, nStep, status
    integer (kind = cwipi_int_l) :: nNotLocatedPoints
    double precision :: timeValue
    double precision, dimension(*) :: sendingField, receivingField

    integer (kind = cwipi_int_l) :: lCouplingName, lExchangeName, lSendingFieldName
    integer (kind = cwipi_int_l) :: lReceivingFieldName

    lCouplingName       = len(couplingName)
    lExchangeName       = len(exchangeName)
    lSendingFieldName   = len(sendingFieldName)
    lReceivingFieldName = len(receivingFieldName)

    call cwipi_exchange_cf(couplingName, &
                           lCouplingName, &
                           exchangeName, &
                           lExchangeName, &
                           stride, &
                           nStep, &
                           timeValue, &
                           sendingFieldName, &
                           lSendingFieldName, &
                           sendingField, &
                           receivingFieldName, &
                           lReceivingFieldName, &
                           receivingField, &
                           nNotLocatedPoints, &
                           status)

  end subroutine cwipi_exch_without_user_itp_f_

!
!********************************************************************************
!
! cwipi_issend_without_user_itp_f
!
! Send data to the coupled application (only send)
! Non-blocking communication 
!
! parameters
!   couplingName         <-- Coupling identifier
!   exchangeName         <-- Exchange name
!   tag                  <-- Exchange tag
!   stride               <-- Number of interlaced fields
!   timeStep             <-- Time step  (only for visualization)
!   timeValue            <-- Time value (only for visualization)
!   sendingFieldName     <-- Sending field name
!   sendingField         <-- Sending field
!   request              -->Sending request
!
!********************************************************************************
!

  subroutine cwipi_issend_without_user_itp_f_(couplingName, &
                                              exchangeName, &
                                              tag, &
                                              stride, &
                                              nStep, &
                                              timeValue, &
                                              sendingFieldName, &
                                              sendingField, &
                                              request)

    implicit none

    character (len = *) :: couplingName, exchangeName, sendingFieldName
    integer (kind = cwipi_int_l) :: stride, nStep, tag, request
    double precision :: timeValue
    double precision, dimension(*) :: sendingField

    integer (kind = cwipi_int_l) :: lCouplingName, lExchangeName, lSendingFieldName

    lCouplingName       = len(couplingName)
    lExchangeName       = len(exchangeName)
    lSendingFieldName   = len(sendingFieldName)

    call cwipi_issend_cf(couplingName, &
                         lCouplingName, &
                         exchangeName, &
                         lExchangeName, &
                         tag, &
                         stride, &
                         nStep, &
                         timeValue, &
                         sendingFieldName, &
                         lSendingFieldName, &
                         sendingField, &
                         request)

  end subroutine cwipi_issend_without_user_itp_f_

 
!********************************************************************************
!
! cwipi_issend_with_user_itp_f
!
! Exchange data with the coupled application (only send)
! It is a synchronization point with the coupled application
!
! parameters
!   couplingName         <-- Coupling identifier
!   exchangeName         <-- Exchange name
!   stride               <-- Number of interlaced fields
!   timeStep             <-- Time step  (only for visualization)
!   timeValue            <-- Time value (only for visualization)
!   sendingFieldName     <-- Sending field name
!   sendingField         <-- Sending field (NULL -> no sending)
!   ptInterpolationFct   <-- Callback for interpolation
!   request              --> Exchange request
!
!********************************************************************************
!
  subroutine cwipi_issend_with_user_itp_f_ (couplingName, &
                                            exchangeName, &
                                            tag, &
                                            stride, &
                                            nStep, &
                                            timeValue, &
                                            sendingFieldName, &
                                            sendingField, &
                                            ptInterpolationFct, &
                                            request)

    implicit none

    interface
       subroutine  ptInterpolationFct(entitiesDim, &
                                      nLocalVertex, &
                                      nLocalElement, &
                                      nLocalPolyhedra, &
                                      nDistantPoint, &
                                      localCoordinates, &
                                      localConnectivityIndex, &
                                      localConnectivity, &
                                      localPolyFaceIndex, &
                                      localPolyCellToFaceConnec, &
                                      localPolyFaceConnecIdx, &
                                      localPolyFaceConnec, &
                                      disPtsCoordinates, &
                                      disPtsLocation, &
                                      disPtsDistance, &
                                      disPtsBaryCoordIdx, &
                                      disPtsBaryCoord, &
                                      stride, &
                                      solverType, &
                                      localField, &
                                      distantField)
         integer (kind = 4) :: entitiesDim
         integer (kind = 4) :: nLocalVertex
         integer (kind = 4) :: nLocalElement
         integer (kind = 4) :: nLocalPolyhedra
         integer (kind = 4) :: nDistantPoint
         double precision, dimension(*) :: localCoordinates
         integer (kind = 4), dimension(*) :: localConnectivityIndex
         integer (kind = 4), dimension(*) :: localConnectivity
         integer (kind = 4), dimension(*) :: localPolyFaceIndex
         integer (kind = 4), dimension(*) :: localPolyCellToFaceConnec
         integer (kind = 4), dimension(*) :: localPolyFaceConnecIdx
         integer (kind = 4), dimension(*) :: localPolyFaceConnec
         double precision, dimension(*) :: disPtsCoordinates
         integer (kind = 4), dimension(*) :: disPtsLocation
         real(kind=4), dimension(*) :: disPtsDistance
         integer (kind = 4), dimension(*) :: disPtsBaryCoordIdx
         double precision, dimension(*) :: disPtsBaryCoord
         integer (kind = 4) :: stride
         integer (kind = 4) :: solverType
         double precision, dimension(*) :: localField
         double precision, dimension(*) :: distantField
       end subroutine ptInterpolationFct
    end interface
    character (len = *) :: couplingName, exchangeName, sendingFieldName
    integer (kind = cwipi_int_l) :: stride, nStep, tag, request
    double precision :: timeValue
    double precision, dimension(*) :: sendingField

    integer (kind = cwipi_int_l) :: lCouplingName, lExchangeName, lSendingFieldName

    lCouplingName       = len(couplingName)
    lExchangeName       = len(exchangeName)
    lSendingFieldName   = len(sendingFieldName)

    call cwipi_issend_with_user_itp_cf(couplingName, &
                                       lCouplingName, &
                                       exchangeName, &
                                       lExchangeName, &
                                       tag, &
                                       stride, &
                                       nStep, &
                                       timeValue, &
                                       sendingFieldName, &
                                       lSendingFieldName, &
                                       sendingField, &
                                       ptInterpolationFct, &
                                       request)

  end subroutine cwipi_issend_with_user_itp_f_




  

  subroutine cwipi_issend_with_user_ho_itp_f_ (couplingName, &
                                            exchangeName, &
                                            meshOrder,&
                                            tag, &
                                            stride, &
                                            nStep, &
                                            timeValue, &
                                            sendingFieldName, &
                                            sendingField, &
                                            ptHoInterpolationFct, &
                                            request)

    implicit none

    interface
      subroutine  ptHoInterpolationFct(entitiesDim, &
                                      order, &
                                      nLocalVertex, &
                                      nLocalElement, &
                                      nLocalPolyhedra, &
                                      nDistantPoint, &
                                      localCoordinates, &
                                      localConnectivityIndex, &
                                      localConnectivity, &
                                      localPolyFaceIndex, &
                                      localPolyCellToFaceConnec, &
                                      localPolyFaceConnecIdx, &
                                      localPolyFaceConnec, &
                                      disPtsCoordinates, &
                                      disPtsLocation, &
                                      disPtsDistance, &
                                      disPtsBaryCoordIdx, &
                                      disPtsBaryCoord, &
                                      dist_uvw, &
                                      stride, &
                                      solverType, &
                                      localField, &
                                      distantField)
         integer (kind = 4) :: entitiesDim
         integer (kind = 4) :: order
         integer (kind = 4) :: nLocalVertex
         integer (kind = 4) :: nLocalElement
         integer (kind = 4) :: nLocalPolyhedra
         integer (kind = 4) :: nDistantPoint
         double precision, dimension(*) :: localCoordinates
         integer (kind = 4), dimension(*) :: localConnectivityIndex
         integer (kind = 4), dimension(*) :: localConnectivity
         integer (kind = 4), dimension(*) :: localPolyFaceIndex
         integer (kind = 4), dimension(*) :: localPolyCellToFaceConnec
         integer (kind = 4), dimension(*) :: localPolyFaceConnecIdx
         integer (kind = 4), dimension(*) :: localPolyFaceConnec
         double precision, dimension(*) :: disPtsCoordinates
         integer (kind = 4), dimension(*) :: disPtsLocation
         real(kind=4), dimension(*) :: disPtsDistance
         integer (kind = 4), dimension(*) :: disPtsBaryCoordIdx
         double precision, dimension(*) :: disPtsBaryCoord
         double precision, dimension(*) :: dist_uvw
         integer (kind = 4) :: stride
         integer (kind = 4) :: solverType
         double precision, dimension(*) :: localField
         double precision, dimension(*) :: distantField
       end subroutine ptHoInterpolationFct
    end interface
    
    
    integer (kind=cwipi_int_l) :: meshOrder
    character (len = *) :: couplingName, exchangeName, sendingFieldName
    integer (kind = cwipi_int_l) :: stride, nStep, tag, request
    double precision :: timeValue
    double precision, dimension(*) :: sendingField

    integer (kind = cwipi_int_l) :: lCouplingName, lExchangeName, lSendingFieldName

    lCouplingName       = len(couplingName)
    lExchangeName       = len(exchangeName)
    lSendingFieldName   = len(sendingFieldName)

    call cwipi_issend_with_user_itp_cf(couplingName, &
                                       lCouplingName, &
                                       exchangeName, &
                                       lExchangeName, &
                                       tag, &
                                       stride, &
                                       nStep, &
                                       timeValue, &
                                       sendingFieldName, &
                                       lSendingFieldName, &
                                       sendingField, &
                                       ptHoInterpolationFct, &
                                       request)

  end subroutine cwipi_issend_with_user_ho_itp_f_

!
!********************************************************************************
!
! cwipi_send_without_user_itp_f
!
! Exchange data with the coupled application (only send)
! It is a synchronization point with the coupled application
!
! parameters
!   couplingName         <-- Coupling identifier
!   exchangeName         <-- Exchange name
!   stride               <-- Number of interlaced fields
!   timeStep             <-- Time step  (only for visualization)
!   timeValue            <-- Time value (only for visualization)
!   sendingFieldName     <-- Sending field name
!   sendingField         <-- Sending field (NULL -> no sending)
!   status               --> Cwipi exchange status
!
!********************************************************************************
!

  subroutine cwipi_send_without_user_itp_f_ (couplingName, &
                                             exchangeName, &
                                             stride, &
                                             nStep, &
                                             timeValue, &
                                             sendingFieldName, &
                                             sendingField, &
                                             status)

    implicit none

    character (len = *) :: couplingName, exchangeName, sendingFieldName
    integer (kind = cwipi_int_l) :: stride, nStep, status
    double precision :: timeValue
    double precision, dimension(*) :: sendingField

    integer (kind = cwipi_int_l) :: lCouplingName, lExchangeName, lSendingFieldName

    lCouplingName       = len(couplingName)
    lExchangeName       = len(exchangeName)
    lSendingFieldName   = len(sendingFieldName)

    call cwipi_send_cf(couplingName, &
                       lCouplingName, &
                       exchangeName, &
                       lExchangeName, &
                       stride, &
                       nStep, &
                       timeValue, &
                       sendingFieldName, &
                       lSendingFieldName, &
                       sendingField, &
                       status)

  end subroutine cwipi_send_without_user_itp_f_

 
!********************************************************************************
!
! cwipi_send_with_user_itp_f
!
! Exchange data with the coupled application (only send)
! It is a synchronization point with the coupled application
!
! parameters
!   couplingName         <-- Coupling identifier
!   exchangeName         <-- Exchange name
!   stride               <-- Number of interlaced fields
!   timeStep             <-- Time step  (only for visualization)
!   timeValue            <-- Time value (only for visualization)
!   sendingFieldName     <-- Sending field name
!   sendingField         <-- Sending field (NULL -> no sending)
!   ptInterpolationFct   <-- Callback for interpolation
!   status               --> Cwipi exchange status
!
!********************************************************************************
!
  subroutine cwipi_send_with_user_itp_f_ (couplingName, &
                                          exchangeName, &
                                          stride, &
                                          nStep, &
                                          timeValue, &
                                          sendingFieldName, &
                                          sendingField, &
                                          ptInterpolationFct, &
                                          status)

    implicit none

    interface
       subroutine  ptInterpolationFct(entitiesDim, &
                                      nLocalVertex, &
                                      nLocalElement, &
                                      nLocalPolyhedra, &
                                      nDistantPoint, &
                                      localCoordinates, &
                                      localConnectivityIndex, &
                                      localConnectivity, &
                                      localPolyFaceIndex, &
                                      localPolyCellToFaceConnec, &
                                      localPolyFaceConnecIdx, &
                                      localPolyFaceConnec, &
                                      disPtsCoordinates, &
                                      disPtsLocation, &
                                      disPtsDistance, &
                                      disPtsBaryCoordIdx, &
                                      disPtsBaryCoord, &
                                      stride, &
                                      solverType, &
                                      localField, &
                                      distantField)
         integer (kind = 4) :: entitiesDim
         integer (kind = 4) :: nLocalVertex
         integer (kind = 4) :: nLocalElement
         integer (kind = 4) :: nLocalPolyhedra
         integer (kind = 4) :: nDistantPoint
         double precision, dimension(*) :: localCoordinates
         integer (kind = 4), dimension(*) :: localConnectivityIndex
         integer (kind = 4), dimension(*) :: localConnectivity
         integer (kind = 4), dimension(*) :: localPolyFaceIndex
         integer (kind = 4), dimension(*) :: localPolyCellToFaceConnec
         integer (kind = 4), dimension(*) :: localPolyFaceConnecIdx
         integer (kind = 4), dimension(*) :: localPolyFaceConnec
         double precision, dimension(*) :: disPtsCoordinates
         integer (kind = 4), dimension(*) :: disPtsLocation
         real(kind=4), dimension(*) :: disPtsDistance
         integer (kind = 4), dimension(*) :: disPtsBaryCoordIdx
         double precision, dimension(*) :: disPtsBaryCoord
         integer (kind = 4) :: stride
         integer (kind = 4) :: solverType
         double precision, dimension(*) :: localField
         double precision, dimension(*) :: distantField
       end subroutine ptInterpolationFct
    end interface
    character (len = *) :: couplingName, exchangeName, sendingFieldName
    integer (kind = cwipi_int_l) :: stride, nStep, status
    double precision :: timeValue
    double precision, dimension(*) :: sendingField

    integer (kind = cwipi_int_l) :: lCouplingName, lExchangeName, lSendingFieldName

    lCouplingName       = len(couplingName)
    lExchangeName       = len(exchangeName)
    lSendingFieldName   = len(sendingFieldName)

    call cwipi_send_with_user_itp_cf(couplingName, &
                                        lCouplingName, &
                                        exchangeName, &
                                        lExchangeName, &
                                        stride, &
                                        nStep, &
                                        timeValue, &
                                        sendingFieldName, &
                                        lSendingFieldName, &
                                        sendingField, &
                                        ptInterpolationFct, &
                                        status)

  end subroutine cwipi_send_with_user_itp_f_





  subroutine cwipi_send_with_user_ho_itp_f_ (couplingName, &
                                          exchangeName, &
                                          meshOrder,&
                                          stride, &
                                          nStep, &
                                          timeValue, &
                                          sendingFieldName, &
                                          sendingField, &
                                          ptHoInterpolationFct, &
                                          status)

    implicit none

    interface
      subroutine  ptHoInterpolationFct(entitiesDim, &
                                      order, &
                                      nLocalVertex, &
                                      nLocalElement, &
                                      nLocalPolyhedra, &
                                      nDistantPoint, &
                                      localCoordinates, &
                                      localConnectivityIndex, &
                                      localConnectivity, &
                                      localPolyFaceIndex, &
                                      localPolyCellToFaceConnec, &
                                      localPolyFaceConnecIdx, &
                                      localPolyFaceConnec, &
                                      disPtsCoordinates, &
                                      disPtsLocation, &
                                      disPtsDistance, &
                                      disPtsBaryCoordIdx, &
                                      disPtsBaryCoord, &
                                      dist_uvw, &
                                      stride, &
                                      solverType, &
                                      localField, &
                                      distantField)
         integer (kind = 4) :: entitiesDim
         integer (kind = 4) :: order
         integer (kind = 4) :: nLocalVertex
         integer (kind = 4) :: nLocalElement
         integer (kind = 4) :: nLocalPolyhedra
         integer (kind = 4) :: nDistantPoint
         double precision, dimension(*) :: localCoordinates
         integer (kind = 4), dimension(*) :: localConnectivityIndex
         integer (kind = 4), dimension(*) :: localConnectivity
         integer (kind = 4), dimension(*) :: localPolyFaceIndex
         integer (kind = 4), dimension(*) :: localPolyCellToFaceConnec
         integer (kind = 4), dimension(*) :: localPolyFaceConnecIdx
         integer (kind = 4), dimension(*) :: localPolyFaceConnec
         double precision, dimension(*) :: disPtsCoordinates
         integer (kind = 4), dimension(*) :: disPtsLocation
         real(kind=4), dimension(*) :: disPtsDistance
         integer (kind = 4), dimension(*) :: disPtsBaryCoordIdx
         double precision, dimension(*) :: disPtsBaryCoord
         double precision, dimension(*) :: dist_uvw
         integer (kind = 4) :: stride
         integer (kind = 4) :: solverType
         double precision, dimension(*) :: localField
         double precision, dimension(*) :: distantField
       end subroutine ptHoInterpolationFct
    end interface
    
    integer (kind=cwipi_int_l) :: meshOrder    
    character (len = *) :: couplingName, exchangeName, sendingFieldName
    integer (kind = cwipi_int_l) :: stride, nStep, status
    double precision :: timeValue
    double precision, dimension(*) :: sendingField

    integer (kind = cwipi_int_l) :: lCouplingName, lExchangeName, lSendingFieldName

    lCouplingName       = len(couplingName)
    lExchangeName       = len(exchangeName)
    lSendingFieldName   = len(sendingFieldName)

    call cwipi_send_with_user_itp_cf(couplingName, &
                                        lCouplingName, &
                                        exchangeName, &
                                        lExchangeName, &
                                        stride, &
                                        nStep, &
                                        timeValue, &
                                        sendingFieldName, &
                                        lSendingFieldName, &
                                        sendingField, &
                                        ptHoInterpolationFct, &
                                        status)

  end subroutine cwipi_send_with_user_ho_itp_f_


!
!********************************************************************************
!
! cwipi_receive_f
!
! Exchange data with the coupled application. (only receive)
! It is a synchronization point with the coupled application
!
! parameters
!   couplingName         <-- Coupling identifier
!   exchangeName         <-- Exchange name
!   stride               <-- Number of interlaced fields
!   timeStep             <-- Time step  (only for visualization)
!   timeValue            <-- Time value (only for visualization)
!   receivingFieldName   <-- Receiving field name
!   receivingField       --> Receiving field
!   nNotLocatedPoints    --> Number of not located points
!   status               --> Cwipi exchange status
!
!********************************************************************************
!

  subroutine cwipi_receive_f_ (couplingName, &
                                   exchangeName, &
                                   stride, &
                                   nStep, &
                                   timeValue, &
                                   receivingFieldName, &
                                   receivingField, &
                                   nNotLocatedPoints, &
                                   status)

    implicit none

    character (len = *) :: couplingName, exchangeName, receivingFieldName
    integer (kind = cwipi_int_l) :: stride, nStep, status
    integer (kind = cwipi_int_l) :: nNotlocatedPoints
    double precision :: timeValue
    double precision, dimension(*) :: receivingField

    integer (kind = cwipi_int_l) :: lCouplingName, lExchangeName, lReceivingFieldName

    lCouplingName       = len(couplingName)
    lExchangeName       = len(exchangeName)
    lReceivingFieldName = len(receivingFieldName)

    call cwipi_receive_cf(couplingName, &
                              lCouplingName, &
                              exchangeName, &
                              lExchangeName, &
                              stride, &
                              nStep, &
                              timeValue, &
                              receivingFieldName, &
                              lReceivingFieldName, &
                              receivingField, &
                              nNotLocatedPoints, &
                              status)

  end subroutine cwipi_receive_f_


!
!********************************************************************************
!
! cwipi_ireceive_f
!
! Receive data from the coupled application.
! Non-blocking communication. receivingField is fuly updated 
! after cwipi_wait_irecv calling
!
! parameters
!   couplingName         <-- Coupling identifier
!   exchangeName         <-- Exchange name
!   tag                  <-- Exchange tag
!   stride               <-- Number of interlaced fields
!   timeStep             <-- Time step  (only for visualization)
!   timeValue            <-- Time value (only for visualization)
!   receivingFieldName   <-- Receiving field name
!   receivingField       <-- Receiving field
!   request              --> Exchange request
!
!********************************************************************************
!

  subroutine cwipi_ireceive_f_ (couplingName, &
                                exchangeName, &
                                tag, &
                                stride, &
                                nStep, &
                                timeValue, &
                                receivingFieldName, &
                                receivingField, &
                                request)

    implicit none

    character (len = *) :: couplingName, exchangeName, receivingFieldName
    integer (kind = cwipi_int_l) :: stride, nStep, tag, request
    double precision :: timeValue
    double precision, dimension(*) :: receivingField

    integer (kind = cwipi_int_l) :: lCouplingName, lExchangeName, lReceivingFieldName

    lCouplingName       = len(couplingName)
    lExchangeName       = len(exchangeName)
    lReceivingFieldName = len(receivingFieldName)

    call cwipi_ireceive_cf(couplingName, &
                           lCouplingName, &
                           exchangeName, &
                           lExchangeName, &
                           tag, &
                           stride, &
                           nStep, &
                           timeValue, &
                           receivingFieldName, &
                           lReceivingFieldName, &
                           receivingField, &
                           request)

  end subroutine cwipi_ireceive_f_


!
!********************************************************************************
!
! cwipi_wait_issend
!
! Wait for cwipi_issend
!
! parameters
!   couplingName         <-- Coupling identifier
!   request              <-- Exchange request
!
!********************************************************************************
!

  subroutine cwipi_wait_issend_f_(couplingName, &
                                  request)

    implicit none

    character (len = *) :: couplingName
    integer (kind = cwipi_int_l) :: request
    integer (kind = cwipi_int_l) :: lCouplingName

    lCouplingName       = len(couplingName)

    call cwipi_wait_issend_cf(couplingName, lCouplingName, request)

  end subroutine cwipi_wait_issend_f_


!
!********************************************************************************
!
! cwipi_wait_irecv
!
! Wait for cwipi_irecv
!
! parameters
!   couplingName         <-- Coupling identifier
!   request              <-- Exchange request
!
!********************************************************************************
!

  subroutine cwipi_wait_irecv_f_(couplingName, &
                                 request)

    implicit none

    character (len = *) :: couplingName
    integer (kind = cwipi_int_l) :: request
    integer (kind = cwipi_int_l) :: lCouplingName

    lCouplingName       = len(couplingName)

    call cwipi_wait_irecv_cf(couplingName, lCouplingName, request)

  end subroutine cwipi_wait_irecv_f_


!
!********************************************************************************
!
! cwipi_exchange_f_with_user_itp_f
!
! Exchange data with the coupled application.
! It is a synchronization point with the coupled application
!
! parameters
!   couplingName         <-- Coupling identifier
!   exchangeName         <-- Exchange name
!   stride               <-- Number of interlaced fields
!   timeStep             <-- Time step  (only for visualization)
!   timeValue            <-- Time value (only for visualization)
!   sendingFieldName     <-- Sending field name
!   sendingField         <-- Sending field (NULL -> no sending)
!   receivingFieldName   <-- Receiving field name
!   receivingField       --> Receiving field
!   ptInterpolationFct   <-- Callback for interpolation
!   nNotLocatedPoints    --> Number of not located points
!   status               --> Cwipi exchange status
!
!********************************************************************************
!

  subroutine cwipi_exch_with_user_itp_f_ (couplingName, &
                                          exchangeName, &
                                          exchangeDim, &
                                          nStep, &
                                          timeValue, &
                                          sendingFieldName, &
                                          sendingField, &
                                          receivingFieldName, &
                                          receivingField, &
                                          ptInterpolationFct, &
                                          nNotLocatedPoints, &
                                          status)

    implicit none

    interface
       subroutine  ptInterpolationFct(entitiesDim, &
                                      nLocalVertex, &
                                      nLocalElement, &
                                      nLocalPolyhedra, &
                                      nDistantPoint, &
                                      localCoordinates, &
                                      localConnectivityIndex, &
                                      localConnectivity, &
                                      localPolyFaceIndex, &
                                      localPolyCellToFaceConnec, &
                                      localPolyFaceConnecIdx, &
                                      localPolyFaceConnec, &
                                      disPtsCoordinates, &
                                      disPtsLocation, &
                                      disPtsDistance, &
                                      disPtsBaryCoordIdx, &
                                      disPtsBaryCoord, &
                                      stride, &
                                      solverType, &
                                      localField, &
                                      distantField)
         integer (kind = 4) :: entitiesDim
         integer (kind = 4) :: nLocalVertex
         integer (kind = 4) :: nLocalElement
         integer (kind = 4) :: nLocalPolyhedra
         integer (kind = 4) :: nDistantPoint
         double precision, dimension(*) :: localCoordinates
         integer (kind = 4), dimension(*) :: localConnectivityIndex
         integer (kind = 4), dimension(*) :: localConnectivity
         integer (kind = 4), dimension(*) :: localPolyFaceIndex
         integer (kind = 4), dimension(*) :: localPolyCellToFaceConnec
         integer (kind = 4), dimension(*) :: localPolyFaceConnecIdx
         integer (kind = 4), dimension(*) :: localPolyFaceConnec
         double precision, dimension(*) :: disPtsCoordinates
         integer (kind = 4), dimension(*) :: disPtsLocation
         real(kind=4), dimension(*) :: disPtsDistance
         integer (kind = 4), dimension(*) :: disPtsBaryCoordIdx
         double precision, dimension(*) :: disPtsBaryCoord
         integer (kind = 4) :: stride
         integer (kind = 4) :: solverType
         double precision, dimension(*) :: localField
         double precision, dimension(*) :: distantField
       end subroutine ptInterpolationFct
    end interface

    character (len = *) :: couplingName, exchangeName, sendingFieldName
    character (len = *) :: receivingFieldName
    integer (kind = cwipi_int_l) :: exchangeDim, nStep, status
    integer (kind = cwipi_int_l) :: nnotlocatedpoints
    double precision :: timeValue
    double precision, dimension(*) ::  sendingField, receivingField

    integer (kind = cwipi_int_l) :: lCouplingName, lExchangeName, lSendingFieldName
    integer (kind = cwipi_int_l) :: lReceivingFieldName

    lCouplingName       = len(couplingName)
    lExchangeName       = len(exchangeName)
    lSendingFieldName   = len(sendingFieldName)
    lReceivingFieldName = len(receivingFieldName)

    call cwipi_exch_with_user_itp_cf(couplingName, &
                                     lCouplingName, &
                                     exchangeName, &
                                     lExchangeName, &
                                     exchangeDim, &
                                     nStep, &
                                     timeValue, &
                                     sendingFieldName, &
                                     lSendingFieldName, &
                                     sendingField, &
                                     receivingFieldName, &
                                     lReceivingFieldName, &
                                     receivingField, &
                                     ptInterpolationFct, &
                                     nnotlocatedpoints, &
                                     status)

  end subroutine cwipi_exch_with_user_itp_f_



  subroutine cwipi_exch_with_user_ho_itp_f_ (couplingName, &
                                          exchangeName, &
                                          exchangeDim, &
                                          meshOrder,&
                                          nStep, &
                                          timeValue, &
                                          sendingFieldName, &
                                          sendingField, &
                                          receivingFieldName, &
                                          receivingField, &
                                          ptHoInterpolationFct, &
                                          nNotLocatedPoints, &
                                          status)

    implicit none

    interface
      subroutine  ptHoInterpolationFct(entitiesDim, &
                                      order, &
                                      nLocalVertex, &
                                      nLocalElement, &
                                      nLocalPolyhedra, &
                                      nDistantPoint, &
                                      localCoordinates, &
                                      localConnectivityIndex, &
                                      localConnectivity, &
                                      localPolyFaceIndex, &
                                      localPolyCellToFaceConnec, &
                                      localPolyFaceConnecIdx, &
                                      localPolyFaceConnec, &
                                      disPtsCoordinates, &
                                      disPtsLocation, &
                                      disPtsDistance, &
                                      disPtsBaryCoordIdx, &
                                      disPtsBaryCoord, &
                                      dist_uvw, &
                                      stride, &
                                      solverType, &
                                      localField, &
                                      distantField)
         integer (kind = 4) :: entitiesDim
         integer (kind = 4) :: order
         integer (kind = 4) :: nLocalVertex
         integer (kind = 4) :: nLocalElement
         integer (kind = 4) :: nLocalPolyhedra
         integer (kind = 4) :: nDistantPoint
         double precision, dimension(*) :: localCoordinates
         integer (kind = 4), dimension(*) :: localConnectivityIndex
         integer (kind = 4), dimension(*) :: localConnectivity
         integer (kind = 4), dimension(*) :: localPolyFaceIndex
         integer (kind = 4), dimension(*) :: localPolyCellToFaceConnec
         integer (kind = 4), dimension(*) :: localPolyFaceConnecIdx
         integer (kind = 4), dimension(*) :: localPolyFaceConnec
         double precision, dimension(*) :: disPtsCoordinates
         integer (kind = 4), dimension(*) :: disPtsLocation
         real(kind=4), dimension(*) :: disPtsDistance
         integer (kind = 4), dimension(*) :: disPtsBaryCoordIdx
         double precision, dimension(*) :: disPtsBaryCoord
         double precision, dimension(*) :: dist_uvw
         integer (kind = 4) :: stride
         integer (kind = 4) :: solverType
         double precision, dimension(*) :: localField
         double precision, dimension(*) :: distantField
       end subroutine ptHoInterpolationFct
    end interface

    integer (kind=cwipi_int_l) :: meshOrder
    character (len = *) :: couplingName, exchangeName, sendingFieldName
    character (len = *) :: receivingFieldName
    integer (kind = cwipi_int_l) :: exchangeDim, nStep, status
    integer (kind = cwipi_int_l) :: nnotlocatedpoints
    double precision :: timeValue
    double precision, dimension(*) ::  sendingField, receivingField

    integer (kind = cwipi_int_l) :: lCouplingName, lExchangeName, lSendingFieldName
    integer (kind = cwipi_int_l) :: lReceivingFieldName

    lCouplingName       = len(couplingName)
    lExchangeName       = len(exchangeName)
    lSendingFieldName   = len(sendingFieldName)
    lReceivingFieldName = len(receivingFieldName)

    call cwipi_exch_with_user_itp_cf(couplingName, &
                                     lCouplingName, &
                                     exchangeName, &
                                     lExchangeName, &
                                     exchangeDim, &
                                     nStep, &
                                     timeValue, &
                                     sendingFieldName, &
                                     lSendingFieldName, &
                                     sendingField, &
                                     receivingFieldName, &
                                     lReceivingFieldName, &
                                     receivingField, &
                                     ptHoInterpolationFct, &
                                     nnotlocatedpoints, &
                                     status)

  end subroutine cwipi_exch_with_user_ho_itp_f_
  
!
!********************************************************************************
!
! cwipi_delete_coupling_f
!
! Delete a coupling
!
! parameters
!   couplingName         <-- Coupling identifier
!
!
!********************************************************************************
!

  subroutine cwipi_delete_coupling_f_(couplingName)

    implicit none

    character (len = *) :: couplingName
    integer (kind = cwipi_int_l) :: lCouplingName

    lCouplingName       = len(couplingName)

    call cwipi_delete_coupling_cf(couplingName, lCouplingName)

  end subroutine cwipi_delete_coupling_f_

!
!********************************************************************************
!
! cwipi_finalize()  (define into cwipi_cf.hxx)
!
! Finalize cwipi. This is a synchronization point between all applications
!
!********************************************************************************
!

!
!********************************************************************************
!
! cwipi_get_not_located_pts_f
!
! Get located points
!
! parameters
!   couplingName         <-- Coupling identifier
!   notLocatedPoints     --> Not located points (1 to n)
!
!********************************************************************************
!
  subroutine cwipi_get_not_located_pts_f_(couplingName, notLocatedPoints)

    implicit none

    character (len = *) :: couplingName
    integer (kind = cwipi_int_l), dimension(*) :: notLocatedPoints
    integer (kind = cwipi_int_l) :: lCouplingName

    lCouplingName       = len(couplingName)

    call cwipi_get_not_located_pts_cf(couplingName, lCouplingName, notLocatedPoints)
  end subroutine cwipi_get_not_located_pts_f_

!
!********************************************************************************
!
! cwipi_get_n_located_dist_pts_f
!
! Get located points
!
! parameters
!   couplingName             <-- Coupling identifier
!   nLocatedDistantPoints     --> Number oflocated distan points
!
!********************************************************************************
!
  subroutine cwipi_get_n_located_dist_pts_f_(couplingName, nLocatedDistantPoints)

    implicit none

    character (len = *) :: couplingName
    integer (kind = cwipi_int_l) :: nLocatedDistantPoints
    integer (kind = cwipi_int_l) :: lCouplingName

    lCouplingName       = len(couplingName)

    call cwipi_get_n_located_dist_pts_cf(couplingName, &
                                                   lCouplingName, &
                                                   nLocatedDistantPoints)
  end subroutine cwipi_get_n_located_dist_pts_f_



!
!********************************************************************************
!
! Get located points
!
! parameters
!   couplingName         <-- Coupling identifier
!   locatedPoints        --> Located points (1 to n)
!
!********************************************************************************
!

  subroutine cwipi_get_located_pts_f_(couplingName, locatedPoints)

    implicit none

    character (len = *) :: couplingName
    integer (kind = cwipi_int_l), dimension(*) :: locatedPoints
    integer (kind = cwipi_int_l) :: lCouplingName

    lCouplingName       = len(couplingName)

    call cwipi_get_located_pts_cf(couplingName, lCouplingName, locatedPoints)
  end subroutine cwipi_get_located_pts_f_

!
!********************************************************************************
!
! Get number of located points
!
! parameters
!   couplingName         <-- Coupling identifier
!   locatedPoints        --> Located points (1 to n)
!
!********************************************************************************
!

  subroutine cwipi_get_n_located_pts_f_(couplingName, nLocatedPoints)

    implicit none

    character (len = *) :: couplingName
    integer (kind = cwipi_int_l) :: nLocatedPoints
    integer (kind = cwipi_int_l) :: lCouplingName

    lCouplingName       = len(couplingName)

    call cwipi_get_n_located_pts_cf(couplingName, &
                                         lCouplingName, &
                                         nLocatedPoints)

  end subroutine cwipi_get_n_located_pts_f_

!
!********************************************************************************
!
! Get number of not located points
!
! parameters
!   couplingName         <-- Coupling identifier
!   locatedPoints        --> Located points (1 to n)
!
!********************************************************************************
!

  subroutine cwipi_get_n_not_located_pts_f_(couplingName, nNotLocatedPoints)

    implicit none

    character (len = *) :: couplingName
    integer (kind = cwipi_int_l) :: nNotLocatedPoints
    integer (kind = cwipi_int_l) :: lCouplingName

    lCouplingName       = len(couplingName)

    call cwipi_get_n_not_located_pts_cf(couplingName, &
                                        lCouplingName, &
                                        nNotLocatedPoints)

  end subroutine cwipi_get_n_not_located_pts_f_

!
!********************************************************************************
!
! Get number of distant ranks 
!
! parameters
!   coupling_id          <-- Coupling identifier
!
! return
!                        --> Number of distant ranks
!
!********************************************************************************
!

  subroutine cwipi_get_n_dis_ranks_f_(couplingName, n_dis_ranks)
    implicit none

    character (len = *) :: couplingName
    integer (kind = cwipi_int_l) :: lCouplingName
    integer (kind = cwipi_int_l) :: n_dis_ranks

    lCouplingName       = len(couplingName)

    call cwipi_get_n_dis_ranks_cf(couplingName,  lCouplingName , n_dis_ranks)

  end subroutine cwipi_get_n_dis_ranks_f_

!
!********************************************************************************
!
! Get distant point distribution on distant ranks (size = n_distant_rank + 1)
!
! parameters
!   coupling_id          <-- Coupling identifier
!
! return
!                             Distant point distribution on distant ranks
!
!********************************************************************************
!

  subroutine cwipi_get_dis_distrib_f_(couplingName, distrib)
    implicit none

    character (len = *) :: couplingName
    integer (kind = cwipi_int_l) :: lCouplingName
    integer (kind = cwipi_int_l) :: distrib(*)

    lCouplingName       = len(couplingName)

    call cwipi_get_dis_distrib_cf(couplingName,  lCouplingName , distrib)

  end subroutine cwipi_get_dis_distrib_f_

!
!********************************************************************************
!
! Get located points distribution on distant ranks (size = n_distant_rank + 1)
!
! parameters
!   coupling_id          <-- Coupling identifier
!
! return
!                            Located points distribution
!
!********************************************************************************
!

  subroutine cwipi_get_loc_pts_distrib_f_(couplingName, distrib)
    implicit none

    character (len = *) :: couplingName
    integer (kind = cwipi_int_l) :: lCouplingName
    integer (kind = cwipi_int_l) :: distrib(*)

    lCouplingName       = len(couplingName)

    call cwipi_get_loc_pts_distrib_cf(couplingName,  lCouplingName , distrib)

  end subroutine cwipi_get_loc_pts_distrib_f_


end module cwipi
