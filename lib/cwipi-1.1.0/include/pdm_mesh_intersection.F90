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

module pdm_mesh_intersection

  use pdm
  use pdm_fortran
  use iso_c_binding
  use pdm_cf_array
  use pdm_extract_part

  implicit none

  interface PDM_mesh_intersection_create
    module procedure PDM_mesh_intersection_create_
  end interface

  interface PDM_mesh_intersection_n_part_set
    module procedure PDM_mesh_intersection_n_part_set_ 
  end interface 
    
  interface PDM_mesh_intersection_compute
    module procedure PDM_mesh_intersection_compute_
  end interface  

  interface PDM_mesh_intersection_mesh_nodal_set
    module procedure PDM_mesh_intersection_mesh_nodal_set_ 
  end interface 

  interface PDM_mesh_intersection_tolerance_set
    module procedure PDM_mesh_intersection_tolerance_set_
  end interface
    
  interface PDM_mesh_intersection_part_set
    module procedure PDM_mesh_intersection_part_set_
  end interface
    
  interface PDM_mesh_intersection_free
    module procedure PDM_mesh_intersection_free_
  end interface

  interface PDM_mesh_intersection_preprocessing_get
    module procedure PDM_mesh_intersection_preprocessing_get_
  end interface  
    
  ! interface PDM_mesh_intersection_part_to_part_get
  !   module procedure PDM_mesh_intersection_part_to_part_get_
  ! end interface
    
  ! interface PDM_mesh_intersection_result_from_a_get
  !   module procedure PDM_mesh_intersection_result_from_a_get_
  ! end interface
    
  ! interface PDM_mesh_intersection_result_from_b_get
  !   module procedure PDM_mesh_intersection_result_from_b_get_
  ! end interface
    
  ! interface PDM_mesh_intersection_elt_volume_get
  !   module procedure PDM_mesh_intersection_elt_volume_get_
  ! end interface
    
  ! interface PDM_mesh_intersection_extract_part_get
  !   module procedure PDM_mesh_intersection_extract_part_get_
  ! end interface

  integer, parameter :: PDM_MESH_INTERSECTION_KIND_PREPROCESS    = 0
  integer, parameter :: PDM_MESH_INTERSECTION_KIND_WEIGHT        = 1
  integer, parameter :: PDM_MESH_INTERSECTION_KIND_UNMERGED_POLY = 2
  integer, parameter :: PDM_MESH_INTERSECTION_KIND_MESH          = 3

  private ::                                  &
    PDM_mesh_intersection_create_,            &
    PDM_mesh_intersection_n_part_set_,        & 
    PDM_mesh_intersection_compute_,           &
    PDM_mesh_intersection_mesh_nodal_set_,    &
    PDM_mesh_intersection_free_,              &
    PDM_mesh_intersection_tolerance_set_,     &
  !   PDM_mesh_intersection_part_to_part_get_,  &
  !   PDM_mesh_intersection_result_from_a_get_, &
  !   PDM_mesh_intersection_result_from_b_get_, &
  !   PDM_mesh_intersection_elt_volume_get_,    &
  !   PDM_mesh_intersection_extract_part_get_
    PDM_mesh_intersection_part_set_


  interface

    function PDM_mesh_intersection_mesh_dimension_get_cf(mi, i_mesh) &
    result (dim) &
    bind (c, name="PDM_mesh_intersection_mesh_dimension_get")

      use iso_c_binding
      implicit none

      type(c_ptr),    value :: mi
      integer(c_int), value :: i_mesh
      integer(c_int)        :: dim

    end function PDM_mesh_intersection_mesh_dimension_get_cf


  end interface

contains


!>
!!
!! \brief Create a mesh intersection structure
!!
!! \param [in]   dim_mesh_a     Mesh A dimension
!! \param [in]   dim_mesh_b     Mesh B dimension
!! \param [in]   project_coeff  Projection coefficient
!! \param [in]   comm           Associated communicator
!! \param [in]   owner          Results Ownership 
!!
!! \return       mi             Pointer to \ref PDM_mesh_intersection object
!!
!!

subroutine PDM_mesh_intersection_create_ (mi,                &
                                          intersection_kind, &
                                          dim_mesh_a,        &
                                          dim_mesh_b,        &
                                          project_coeff,     &
                                          f_comm,            &
                                          owner)
  
  use iso_c_binding

  implicit none

  integer, intent(in)          :: intersection_kind
  integer, intent(in)          :: dim_mesh_a
  integer, intent(in)          :: dim_mesh_b
  double precision, intent(in) :: project_coeff
  integer, intent(in)          :: f_comm
  integer, intent(in)          :: owner
  
  type(c_ptr)                  :: mi

  interface
    function PDM_mesh_intersection_create_cf (intersection_kind, &
                                              dim_mesh_a,        &
                                              dim_mesh_b,        &
                                              project_coeff,     &
                                              comm,              &
                                              owner )            &
      result(mi) &
      bind (c, name = 'PDM_mesh_intersection_create')

      use iso_c_binding

      implicit none

      integer(c_int), value :: intersection_kind
      integer(c_int), value :: dim_mesh_a
      integer(c_int), value :: dim_mesh_b
      real(c_double), value :: project_coeff
      integer(c_int), value :: comm
      integer(c_int), value :: owner
  
      type(c_ptr)           :: mi

    end function PDM_mesh_intersection_create_cf
  end interface

  integer(c_int) :: c_intersection_kind
  integer(c_int) :: c_dim_mesh_a
  integer(c_int) :: c_dim_mesh_b
  real(c_double) :: c_project_coeff
  integer(c_int) :: c_comm
  integer(c_int) :: c_owner

  c_intersection_kind = intersection_kind
  c_dim_mesh_a        = dim_mesh_a
  c_dim_mesh_b        = dim_mesh_b
  c_project_coeff     = project_coeff
  c_comm              = PDM_MPI_Comm_f2c(f_comm)
  c_owner             = owner

  mi = PDM_mesh_intersection_create_cf (c_intersection_kind, &
                                        c_dim_mesh_a,        &
                                        c_dim_mesh_b,        &
                                        c_project_coeff,     &
                                        c_comm,              &
                                        c_owner )

end subroutine PDM_mesh_intersection_create_


!>
!!
!! \brief Set global data of a mesh
!!
!! \param [in]   mi             Pointer to \ref PDM_mesh_intersection object
!! \param [in]   i_mesh         Mesh identifier
!! \param [in]   n_part         Number of partitions
!!
!!

subroutine PDM_mesh_intersection_n_part_set_ (mi,     &
                                              i_mesh, &
                                              n_part)
 
  use iso_c_binding

  implicit none

  type(c_ptr), intent(in) :: mi
  integer, intent(in)     :: i_mesh
  integer, intent(in)     :: n_part

  interface 
    subroutine PDM_mesh_intersection_n_part_set_cf (mi,     &
                                                    i_mesh, &
                                                    n_part) &
      bind (c, name = "PDM_mesh_intersection_n_part_set")

      use iso_c_binding

      implicit none

        type(c_ptr),     value :: mi
        integer (c_int), value :: i_mesh
        integer (c_int), value :: n_part

    end subroutine PDM_mesh_intersection_n_part_set_cf
  end interface 

  call PDM_mesh_intersection_n_part_set_cf (mi, i_mesh, n_part)

end subroutine PDM_mesh_intersection_n_part_set_   


!>
!!
!! \brief  Compute meshes intersection
!!
!! \param [in]   mi             Pointer to \ref PDM_mesh_intersection object
!!
!!

subroutine PDM_mesh_intersection_compute_ (mi)  
  
  use iso_c_binding

  implicit none

  type(c_ptr), intent(in) :: mi

  interface 
    subroutine PDM_mesh_intersection_compute_cf (mi) &
      bind (c, name = "PDM_mesh_intersection_compute")

      use iso_c_binding

      implicit none

        type(c_ptr),     value :: mi

    end subroutine PDM_mesh_intersection_compute_cf
  end interface 

  call PDM_mesh_intersection_compute_cf (mi)

end subroutine PDM_mesh_intersection_compute_   


!>
!!
!! \brief Set the mesh nodal
!!
!! \param [in]   mi             Pointer to \ref PDM_mesh_intersection object
!! \param [in]   i_mesh         Mesh identifier
!! \param [in]   mesh           Pointer to \ref PDM_part_mesh_nodal object
!!
!!

subroutine PDM_mesh_intersection_mesh_nodal_set_ (mi,     &
                                                  i_mesh, &
                                                  mesh)
  
  
  use iso_c_binding

  implicit none

  type(c_ptr), intent(in) :: mi
  integer, intent(in)     :: i_mesh
  type(c_ptr), intent(in) :: mesh

  interface 
    subroutine PDM_mesh_intersection_mesh_nodal_set_cf (mi,     &
                                                        i_mesh, &
                                                        mesh) &
      bind (c, name = "PDM_mesh_intersection_mesh_nodal_set")

      use iso_c_binding

      implicit none

        type(c_ptr),     value :: mi
        integer (c_int), value :: i_mesh
        type(c_ptr),     value :: mesh

    end subroutine PDM_mesh_intersection_mesh_nodal_set_cf
  end interface 

  call PDM_mesh_intersection_mesh_nodal_set_cf (mi, i_mesh, mesh)

end subroutine PDM_mesh_intersection_mesh_nodal_set_   

!>
!!
!! \brief Set a mesh partition  
!!
!! \param [in]   mi             Pointer to \ref PDM_mesh_intersection object
!! \param [in]   i_mesh         Mesh identifier
!! \param [in]   i_part         Partition identifier
!! \param [in]   n_cell         Number of cells
!! \param [in]   n_face         Number of faces
!! \param [in]   n_edge         Number of edges
!! \param [in]   n_vtx          Number of vertices
!! \param [in]   cell_face_idx  Index in the cell -> face connectivity
!! \param [in]   cell_face      Cell -> face connectivity
!! \param [in]   face_edge_idx  Index in the face -> edge connectivity
!! \param [in]   face_edge      Cell -> face connectivity
!! \param [in]   edge_vtx       Edge -> vertex conectivity 
!! \param [in]   face_vtx_idx   Index in the face -> vertex connectivity
!! \param [in]   face_vtx       Face -> vertex connectivity
!! \param [in]   cell_ln_to_gn  Local cell numbering to global cel numbering
!! \param [in]   face_ln_to_gn  Local face numbering to global face numbering
!! \param [in]   edge_ln_to_gn  Local edge numbering to global edge numbering
!! \param [in]   vtx_ln_to_gn   Local vertex numbering to global vertex numbering
!! \param [in]   vtx_coord      Vertex coordinates
!!

subroutine PDM_mesh_intersection_part_set_ (mi,            &
                                            i_mesh,        &
                                            i_part,        &
                                            n_cell,        &
                                            n_face,        &
                                            n_edge,        &
                                            n_vtx,         &
                                            cell_face_idx, &
                                            cell_face,     &
                                            face_edge_idx, &
                                            face_edge,     &
                                            edge_vtx,      &
                                            face_vtx_idx,  &
                                            face_vtx,      &
                                            cell_ln_to_gn, &
                                            face_ln_to_gn, &
                                            edge_ln_to_gn, &
                                            vtx_ln_to_gn,  &
                                            coords)      
  
  use iso_c_binding

  implicit none

  type(c_ptr), intent(in)            :: mi
  integer, intent(in)                :: i_mesh
  integer, intent(in)                :: i_part
  integer, intent(in)                :: n_cell
  integer, intent(in)                :: n_face
  integer, intent(in)                :: n_edge
  integer, intent(in)                :: n_vtx
  integer(kind=pdm_l_num_s), pointer :: cell_face_idx(:)
  integer(kind=pdm_l_num_s), pointer :: cell_face(:)
  integer(kind=pdm_l_num_s), pointer :: face_edge_idx(:)
  integer(kind=pdm_l_num_s), pointer :: face_edge(:)
  integer(kind=pdm_l_num_s), pointer :: edge_vtx(:)
  integer(kind=pdm_l_num_s), pointer :: face_vtx_idx(:)
  integer(kind=pdm_l_num_s), pointer :: face_vtx(:)
  integer(kind=pdm_g_num_s), pointer :: cell_ln_to_gn(:)
  integer(kind=pdm_g_num_s), pointer :: face_ln_to_gn(:)
  integer(kind=pdm_g_num_s), pointer :: edge_ln_to_gn(:)
  integer(kind=pdm_g_num_s), pointer :: vtx_ln_to_gn(:)
  double precision,          pointer :: coords(:,:)

  type(c_ptr) :: c_cell_face_idx
  type(c_ptr) :: c_cell_face
  type(c_ptr) :: c_face_edge_idx
  type(c_ptr) :: c_face_edge
  type(c_ptr) :: c_edge_vtx
  type(c_ptr) :: c_face_vtx_idx
  type(c_ptr) :: c_face_vtx
  type(c_ptr) :: c_cell_ln_to_gn
  type(c_ptr) :: c_face_ln_to_gn
  type(c_ptr) :: c_edge_ln_to_gn
  type(c_ptr) :: c_vtx_ln_to_gn
  type(c_ptr) :: c_vtx_coord

  interface 
    subroutine PDM_mesh_intersection_part_set_cf (mi,            &
                                                  i_mesh,        &
                                                  i_part,        &
                                                  n_cell,        &
                                                  n_face,        &
                                                  n_edge,        &
                                                  n_vtx,         &
                                                  cell_face_idx, &
                                                  cell_face,     &
                                                  face_edge_idx, &
                                                  face_edge,     &
                                                  edge_vtx,      &
                                                  face_vtx_idx,  &
                                                  face_vtx,      &
                                                  cell_ln_to_gn, &
                                                  face_ln_to_gn, &
                                                  edge_ln_to_gn, &
                                                  vtx_ln_to_gn,  &
                                                  vtx_coord)     &
      bind (c, name = "PDM_mesh_intersection_part_set")

      use iso_c_binding

      implicit none

        type(c_ptr),     value :: mi
        integer (c_int), value :: i_mesh
        integer (c_int), value :: i_part
        integer (c_int), value :: n_cell
        integer (c_int), value :: n_face
        integer (c_int), value :: n_edge
        integer (c_int), value :: n_vtx
        type(c_ptr),     value :: cell_face_idx
        type(c_ptr),     value :: cell_face
        type(c_ptr),     value :: face_edge_idx
        type(c_ptr),     value :: face_edge
        type(c_ptr),     value :: edge_vtx
        type(c_ptr),     value :: face_vtx_idx
        type(c_ptr),     value :: face_vtx
        type(c_ptr),     value :: cell_ln_to_gn
        type(c_ptr),     value :: face_ln_to_gn
        type(c_ptr),     value :: edge_ln_to_gn
        type(c_ptr),     value :: vtx_ln_to_gn
        type(c_ptr),     value :: vtx_coord

    end subroutine PDM_mesh_intersection_part_set_cf
  end interface 

  c_cell_face_idx = C_NULL_PTR
  if (associated(cell_face_idx)) then
    c_cell_face_idx = c_loc(cell_face_idx)
  endif

  c_cell_face = C_NULL_PTR
  if (associated(cell_face)) then
    c_cell_face = c_loc(cell_face)
  endif

  c_face_edge_idx = C_NULL_PTR
  if (associated(face_edge_idx)) then
    c_face_edge_idx = c_loc(face_edge_idx)
  endif

  c_face_edge = C_NULL_PTR
  if (associated(face_edge)) then
    c_face_edge = c_loc(face_edge)
  endif

  c_edge_vtx = C_NULL_PTR
  if (associated(edge_vtx)) then
    c_edge_vtx = c_loc(edge_vtx)
  endif

  c_face_vtx_idx = C_NULL_PTR
  if (associated(face_vtx_idx)) then
    c_face_vtx_idx = c_loc(face_vtx_idx)
  endif

  c_face_vtx = C_NULL_PTR
  if (associated(face_vtx)) then
    c_face_vtx = c_loc(face_vtx)
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

  c_vtx_coord = C_NULL_PTR
  if (associated(coords)) then
    c_vtx_coord = c_loc(coords)
  endif

  call PDM_mesh_intersection_part_set_cf (mi,                   &
                                          i_mesh,               &
                                          i_part,               &
                                          n_cell,               &
                                          n_face,               &
                                          n_edge,               &
                                          n_vtx,                &
                                          c_cell_face_idx, &
                                          c_cell_face,     &
                                          c_face_edge_idx, &
                                          c_face_edge,     &
                                          c_edge_vtx,      &
                                          c_face_vtx_idx,  &
                                          c_face_vtx,      &
                                          c_cell_ln_to_gn, &
                                          c_face_ln_to_gn, &
                                          c_edge_ln_to_gn, &
                                          c_vtx_ln_to_gn,  &
                                          c_vtx_coord)                                          

end subroutine PDM_mesh_intersection_part_set_   

!>
!!
!! \brief A mesh intersection structure  
!!
!! \param [in]   mi             Pointer to \ref PDM_mesh_intersection object
!!

subroutine PDM_mesh_intersection_free_ (mi)  
  
  use iso_c_binding

  implicit none

  type(c_ptr), intent(in) :: mi

  interface 
    subroutine PDM_mesh_intersection_free_cf (mi) &
      bind (c, name = "PDM_mesh_intersection_free")

      use iso_c_binding

      implicit none

        type(c_ptr),     value :: mi

    end subroutine PDM_mesh_intersection_free_cf
  end interface 

  call PDM_mesh_intersection_free_cf (mi)

end subroutine PDM_mesh_intersection_free_   


!>
!! \brief Get part_to_part object to exchange data between the intersected meshes
!!
!! \param [in ] mi         Pointer to \ref PDM_mesh_intersection_t object
!! \param [out] ptp        Pointer to \ref PDM_part_to_part_t object
!! \param [in ] ownership  Ownership for ptp
!!

subroutine PDM_mesh_intersection_part_to_part_get_ (mi,       &
                                                    ptp,      &
                                                    ownership)
  
  use iso_c_binding

  implicit none

  type(c_ptr), intent(in) :: mi
  type(c_ptr)             :: ptp
  integer, intent(in)     :: ownership

  interface 
    subroutine PDM_mesh_intersection_part_to_part_get_cf (mi,        &
                                                          ptp,       &
                                                          ownership) &
      bind (c, name = "PDM_mesh_intersection_part_to_part_get")

      use iso_c_binding

      implicit none

        type(c_ptr),     value :: mi
        type(c_ptr)            :: ptp
        integer (c_int), value :: ownership

    end subroutine PDM_mesh_intersection_part_to_part_get_cf
  end interface 

  call PDM_mesh_intersection_part_to_part_get_cf (mi, ptp, ownership)

end subroutine PDM_mesh_intersection_part_to_part_get_   


!>
!! \brief Get intersection result for the a point of view
!!
!! \param [in ] mi                 Pointer to \ref PDM_mesh_intersection_t object
!! \param [out] ipart              Partition index
!! \param [out] elt_a_elt_b_idx    Index of list of intersected B elements for each A element 
!! \param [out] elt_a_elt_b        List of intersected B elements for each A element 
!! \param [out] elt_a_elt_b_volume Volume of each intersection 
!!


! subroutine PDM_mesh_intersection_result_from_a_get_ (mi,                &
!                                                      elt_a_elt_b_idx,   &
!                                                      elt_a_elt_b,       & 
!                                                      elt_a_elt_b_volume)
  
!   use iso_c_binding

!   implicit none

!   type(c_ptr), intent(in) :: mi
!   type(PDM_cf_array_t), pointer    :: elt_a_elt_b_idx     
!   type(PDM_cf_array_t), pointer    :: elt_a_elt_b         
!   type(PDM_cf_array_t), pointer    :: elt_a_elt_b_volume        

!   type(c_ptr)             :: c_elt_a_elt_b_idx
!   type(c_ptr)             :: c_elt_a_elt_b
!   type(c_ptr)             :: c_elt_a_elt_b_volume

!   interface 
!     subroutine PDM_mesh_intersection_result_from_a_get_cf (mi,                &
!                                                        elt_a_elt_b_idx,   &
!                                                        elt_a_elt_b,       & 
!                                                        elt_a_elt_b_volume)  &                                                 
!       bind (c, name = "PDM_mesh_intersection_result_from_a_get")

!       use iso_c_binding

!       implicit none

!         type(c_ptr),     value :: mi
!         type(c_ptr)            :: elt_a_elt_b_idx
!         type(c_ptr)            :: elt_a_elt_b
!         type(c_ptr)            :: elt_a_elt_b_volume        

!     end subroutine PDM_mesh_intersection_result_from_a_get_cf
!   end interface 

!   elt_a_elt_b_idx => null()
!   elt_a_elt_b => null()
!   elt_a_elt_b_volume => null()

!   call PDM_mesh_intersection_result_from_a_get_cf (mi,                  &
!                                                   c_elt_a_elt_b_idx,   &
!                                                   c_elt_a_elt_b,       &
!                                                   c_elt_a_elt_b_volume)


!   print *," A finir !!!!"
!   !! TODO: creation des pdm_cf_array pour c_elt_a_elt_b_idx c_elt_a_elt_b c_elt_a_elt_b_volume
!   !!

!    ! call PDM_cf_array_create_ (pa,        &
!    !                            ! type,      &
!    !                            ! length,    &
!    !                            ! cptr,      &
!    !                            ! ownership)


! end subroutine PDM_mesh_intersection_result_from_a_get_   

! void
! PDM_mesh_intersection_result_from_a_get
! (
!        PDM_mesh_intersection_t  *mi,
!  const int                       ipart,
!        int                     **elt_a_elt_b_idx,
!        PDM_g_num_t             **elt_a_elt_b,
!        double                  **elt_a_elt_b_volume
! );


!>
!! \brief Get intersection result for the b point of view
!! 
!! TODO: Do as PDM_mesh_intersection_result_from_a_get
!!
!! \param [in ] mi                 Pointer to \ref PDM_mesh_intersection_t object
!! \param [out] ipart              Partition index
!! \param [out] elt_a_elt_b_idx    Index of list of intersected B elements for each A element 
!! \param [out] elt_a_elt_b        List of intersected B elements for each A element 
!! \param [out] elt_a_elt_b_volume Volume of each intersection 
!!

! void
! PDM_mesh_intersection_result_from_b_get
! (
!        PDM_mesh_intersection_t  *mi,
!  const int                       ipart,
!        double                  **elt_b_elt_a_volume
!  );


!>
!! \brief Get intersection result for the b point of view
!! 
!! TODO: Do as PDM_mesh_intersection_result_from_a_get
!!
!! \param [in ] mi                 Pointer to \ref PDM_mesh_intersection_t object
!! \param [out] ipart              Partition index
!! \param [out] elt_a_elt_b_idx    Index of list of intersected B elements for each A element 
!! \param [out] elt_a_elt_b        List of intersected B elements for each A element 
!! \param [out] elt_a_elt_b_volume Volume of each intersection 
!!

! void
! PDM_mesh_intersection_elt_volume_get
! (
!        PDM_mesh_intersection_t  *mi,
!  const int                       imesh,
!  const int                       ipart,
!        double                  **elt_volume
!  );


!>
!!
!! \brief Set the tolerance for bounding boxes
!!
!! \param [in]   mi              Pointer to \ref PDM_mesh_intersection object
!! \param [in]   tol             Tolerance
!!

subroutine PDM_mesh_intersection_tolerance_set_  (mi ,&
                                                 tol)


  use iso_c_binding

  implicit none

  type(c_ptr), intent(in)      :: mi
  double precision, intent(in) :: tol

  interface
    subroutine PDM_mesh_intersection_tolerance_set_cf  (mi ,&
                                                        tol) &
      bind (c, name = "PDM_mesh_intersection_tolerance_set")

      use iso_c_binding

      implicit none

      type(c_ptr), value      :: mi
      real(c_double), value   :: tol

    end subroutine
  end interface

  call PDM_mesh_intersection_tolerance_set_cf (mi, tol)

end subroutine PDM_mesh_intersection_tolerance_set_

!>
!! \brief Get preprocessing results 
!!
!! \param [in ] mi                 Pointer to \ref PDM_mesh_intersection_t object
!! \param [out] elt_a_elt_b_idx    Index of list of intersected B element candidate for each A element
!!                                 in the extr_mesh distribution 
!! \param [out] elt_a_elt_b        List of intersected B element candidate for each A element in the 
!!                                 extr_mesh distribution 
!! \param [out] extr_mesh_a        Redistributed mesh A with only A element candidate  
!! \param [out] extr_mesh_b        Redistributed mesh B with only B element candidate  
!!
!!

subroutine PDM_mesh_intersection_preprocessing_get_(mi,             &
                                                    box_a_box_b_idx, &
                                                    box_a_box_b,     &
                                                    extr_mesh_a,     &
                                                    extr_mesh_b)

  use iso_c_binding

  implicit none

  type(c_ptr), intent(in) :: mi
  integer, pointer        :: box_a_box_b_idx(:)
  integer, pointer        :: box_a_box_b(:)
  type(c_ptr)             :: extr_mesh_a
  type(c_ptr)             :: extr_mesh_b

  integer                 :: n_elt_a, dim_a
  type(c_ptr)             :: c_box_a_box_b_idx
  type(c_ptr)             :: c_box_a_box_b

  integer                 :: i
  interface 
    subroutine PDM_mesh_intersection_preprocessing_get_cf(mi,            &
                                                          box_a_box_b_idx,&
                                                          box_a_box_b,    &
                                                          extr_mesh_a,    &
                                                          extr_mesh_b)    &

      bind (c, name = "PDM_mesh_intersection_preprocessing_get")

      use iso_c_binding

      implicit none
      
      type(c_ptr), value      :: mi
      type(c_ptr)             :: box_a_box_b_idx
      type(c_ptr)             :: box_a_box_b
      type(c_ptr)             :: extr_mesh_a
      type(c_ptr)             :: extr_mesh_b

    end subroutine PDM_mesh_intersection_preprocessing_get_cf
  end interface 


  call PDM_mesh_intersection_preprocessing_get_cf (mi,                &
                                                   c_box_a_box_b_idx, &
                                                   c_box_a_box_b,     &
                                                   extr_mesh_a,       &
                                                   extr_mesh_b)

  ! Get dimension of mesh A
  dim_a = PDM_mesh_intersection_mesh_dimension_get_cf(mi, 0)

  select case (dim_a)
    case (1)
      call PDM_extract_part_n_entity_get(extr_mesh_a,            &
                                         0,                      &  ! ipart (only one)
                                         PDM_MESH_ENTITY_EDGE,   &
                                         n_elt_a)
    case (2)
      call PDM_extract_part_n_entity_get(extr_mesh_a,            &
                                         0,                      &  ! ipart (only one)
                                         PDM_MESH_ENTITY_FACE,   &
                                         n_elt_a)
    case (3)
      call PDM_extract_part_n_entity_get(extr_mesh_a,            &
                                         0,                      &  ! ipart (only one)
                                         PDM_MESH_ENTITY_CELL,   &
                                         n_elt_a)
    case default
      n_elt_a = 0
  end select


  call c_f_pointer (c_box_a_box_b_idx, box_a_box_b_idx, [n_elt_a + 1])
  call c_f_pointer (c_box_a_box_b, box_a_box_b, [box_a_box_b_idx(n_elt_a + 1)])

  do i = 1, box_a_box_b_idx(n_elt_a + 1)
    box_a_box_b(i) = box_a_box_b(i) + 1
  enddo   

  ! 0-based to 1-based !

end subroutine PDM_mesh_intersection_preprocessing_get_

end module pdm_mesh_intersection


