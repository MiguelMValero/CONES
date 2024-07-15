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

module pdm_dcube_gen

  use pdm

  implicit none

  interface PDM_dcube_gen_init ; module procedure &
  pdm_dcube_gen_init_
  end interface

  interface PDM_dcube_gen_dim_get ; module procedure &
  pdm_dcube_gen_dim_get_
  end interface

  interface PDM_dcube_gen_data_get ; module procedure &
  pdm_dcube_gen_data_get_
  end interface

  private :: pdm_dcube_gen_init_
  private :: pdm_dcube_gen_dim_get_
  private :: pdm_dcube_gen_data_get_

  interface


 !>
 !!
 !! \brief Create a distributed cube
 !!
 !! \param [in]   comm           Communicator
 !! \param [in]   n_vtx_seg        Number of vertices in segments
 !! \param [in]   length         Segment length
 !! \param [in]   zero_x         Coordinates of the origin
 !! \param [in]   zero_y         Coordinates of the origin
 !! \param [in]   zero_z         Coordinates of the origin
 !!
 !! \return      Pointer to \ref PDM_dcube_t object
 !!
 !!
    function PDM_dcube_gen_init_cf (comm,          &
                                    n_vtx_seg,     &
                                    length,        &
                                    zero_x,        &
                                    zero_y,        &
                                    zero_z,        &
                                    owner)         &
                                    result (dcube) &
      bind (c, name = 'PDM_dcube_gen_init')

      use iso_c_binding
      implicit none

      integer(c_int), value :: comm
      integer(c_int), value :: owner

#ifdef PDM_LONG_G_NUM
      integer (c_long),  value :: n_vtx_seg
#else
      integer (c_int),  value :: n_vtx_seg
#endif
      real(c_double), value :: length
      real(c_double), value :: zero_x, zero_y, zero_z

      type (c_ptr) :: dcube

    end function PDM_dcube_gen_init_cf


!>
!!
!! \brief Return distributed cube size
!!
!! \param [in]   dcube         Pointer to \ref PDM_dcube_t object
!! \param [out]  n_face_group  Number of faces groups
!! \param [out]  dn_cell       Number of cells stored in this process
!! \param [out]  dn_face       Number of faces stored in this process
!! \param [out]  dn_vtx        Number of vertices stored in this process
!! \param [out]  sface_vtx     Length of dface_vtx array
!! \param [out]  sface_group   Length of dface_group array
!!
!!

subroutine PDM_dcube_gen_dim_get_cf (dcube,        &
                                     n_face_group, &
                                     dn_cell,      &
                                     dn_face,      &
                                     dn_vtx,       &
                                     sface_vtx,    &
                                     sface_group)  &
  bind (c, name = 'PDM_dcube_gen_dim_get')

  use iso_c_binding
  implicit none

  type(c_ptr), value :: dcube

  integer(c_int) :: n_face_group
  integer(c_int) :: dn_cell
  integer(c_int) :: dn_face
  integer(c_int) :: dn_vtx
  integer(c_int) :: sface_vtx
  integer(c_int) :: sface_group

end subroutine PDM_dcube_gen_dim_get_cf


!>
!!
!! \brief Return distributed cube data
!!
!! \param [in]  dcube           Pointer to \ref PDM_dcube_t object
!! \param [out] dface_cell      Faces from cells connectivity (size = 2 * dn_face)
!! \param [out] dface_vtx_idx   Faces from vertices connectivity index (size = dn_face + 1)
!! \param [out] dface_vtx       Faces from vertices connectivity (size = dface_vtxL)
!! \param [out] dvtx_coord      Vertices coordinates (size = 3 * dn_vtx)
!! \param [out] dface_group_idx Faces groups index (size = n_face_group + 1)
!! \param [out] dface_group     Faces groups (size = dFacegroupL)
!!
!!

subroutine PDM_dcube_gen_data_get_cf (dcube,           &
                                      dface_cell,      &
                                      dface_vtx_idx,   &
                                      dface_vtx,       &
                                      dvtx_coord,      &
                                      dface_group_idx, &
                                      dface_group)     &
  bind (c, name = 'PDM_dcube_gen_data_get')

  use iso_c_binding
  implicit none

  type(c_ptr), value :: dcube
  type(c_ptr)        :: dface_cell
  type(c_ptr)        :: dface_vtx_idx
  type(c_ptr)        :: dface_vtx
  type(c_ptr)        :: dvtx_coord
  type(c_ptr)        :: dface_group_idx
  type(c_ptr)        :: dface_group

end subroutine PDM_dcube_gen_data_get_cf


 !>
 !!
 !! \brief Free a distributed cube
 !!
 !! \param [in]  dcube     Pointer to \ref PDM_dcube_t object
 !!
 !!

    subroutine PDM_dcube_gen_free (dcube)     &
      bind (c, name = 'PDM_dcube_gen_free')
      use iso_c_binding
      implicit none

      type (c_ptr)  , value :: dcube

    end subroutine PDM_dcube_gen_free

  end interface


  contains


 !>
 !!
 !! \brief Create a distributed cube
 !!
 !! \param [out]  dcube       Pointer to \ref PDM_dcube_t object
 !! \param [in]   comm        Communicator
 !! \param [in]   n_vtx_seg   Number of vertices in segments
 !! \param [in]   length      Segment length
 !! \param [in]   zero_x      Coordinates of the origin
 !! \param [in]   zero_y      Coordinates of the origin
 !! \param [in]   zero_z      Coordinates of the origin
 !!

  subroutine PDM_dcube_gen_init_ (dcube,         &
                                  f_comm,        &
                                  n_vtx_seg,     &
                                  length,        &
                                  zero_x,        &
                                  zero_y,        &
                                  zero_z,        &
                                  owner)

      use iso_c_binding
      implicit none

      integer,              intent(in) :: f_comm
      integer,              intent(in) :: owner
      integer(pdm_g_num_s), intent(in) :: n_vtx_seg
      double precision,     intent(in) :: length
      double precision,     intent(in) :: zero_x, zero_y, zero_z
      type (c_ptr)                     :: dcube

      integer(c_int)               :: c_comm
      integer(c_int)               :: c_owner
#ifdef PDM_LONG_G_NUM
      integer (c_long)             :: c_n_vtx_seg
#else
      integer (c_int)              :: c_n_vtx_seg
#endif
      real(c_double)               :: c_length
      real(c_double)               :: c_zero_x, c_zero_y, c_zero_z


      c_comm = PDM_MPI_Comm_f2c(f_comm)

      c_owner     = owner
      c_n_vtx_seg = n_vtx_seg
      c_length    = length
      c_zero_x    = zero_x
      c_zero_y    = zero_y
      c_zero_z    = zero_z

      dcube = PDM_dcube_gen_init_cf(c_comm,       &
                                    c_n_vtx_seg,  &
                                    c_length,     &
                                    c_zero_x,     &
                                    c_zero_y,     &
                                    c_zero_z,     &
                                    c_owner)

    end subroutine PDM_dcube_gen_init_




!>
!!
!! \brief Return distributed cube size
!!
!! \param [in]   dcube         Pointer to \ref PDM_dcube_t object
!! \param [out]  n_face_group  Number of faces groups
!! \param [out]  dn_cell       Number of cells stored in this process
!! \param [out]  dn_face       Number of faces stored in this process
!! \param [out]  dn_vtx        Number of vertices stored in this process
!! \param [out]  sface_vtx     Length of dface_vtx array
!! \param [out]  sface_group   Length of dface_group array
!!
!!

subroutine PDM_dcube_gen_dim_get_ (dcube,        &
                                   n_face_group, &
                                   dn_cell,      &
                                   dn_face,      &
                                   dn_vtx,       &
                                   sface_vtx,    &
                                   sface_group)
  use iso_c_binding
  implicit none

  type(c_ptr), value   :: dcube

  integer, intent(out) :: n_face_group
  integer, intent(out) :: dn_cell
  integer, intent(out) :: dn_face
  integer, intent(out) :: dn_vtx
  integer, intent(out) :: sface_vtx
  integer, intent(out) :: sface_group

  integer(c_int)       :: c_n_face_group
  integer(c_int)       :: c_dn_cell
  integer(c_int)       :: c_dn_face
  integer(c_int)       :: c_dn_vtx
  integer(c_int)       :: c_sface_vtx
  integer(c_int)       :: c_sface_group

  call PDM_dcube_gen_dim_get_cf(dcube,          &
                                c_n_face_group, &
                                c_dn_cell,      &
                                c_dn_face,      &
                                c_dn_vtx,       &
                                c_sface_vtx,    &
                                c_sface_group)

  n_face_group = c_n_face_group
  dn_cell      = c_dn_cell
  dn_face      = c_dn_face
  dn_vtx       = c_dn_vtx
  sface_vtx    = c_sface_vtx
  sface_group  = c_sface_group

end subroutine PDM_dcube_gen_dim_get_



  !>
  !!
  !! \brief Return distributed cube data
  !!
  !! \param [in]  dcube           Pointer to \ref PDM_dcube_t object
  !! \param [out] dface_cell      Faces from cells connectivity (size = 2 * dn_face)
  !! \param [out] dface_vtx_idx   Faces from vertices connectivity index (size = dn_face + 1)
  !! \param [out] dface_vtx       Faces from vertices connectivity (size = dface_vtxL)
  !! \param [out] dvtx_coord      Vertices coordinates (size = 3 * dn_vtx)
  !! \param [out] dface_group_idx Faces groups index (size = n_face_group + 1)
  !! \param [out] dface_group     Faces groups (size = dFacegroupL)
  !!
  !!

  subroutine PDM_dcube_gen_data_get_ (dcube,           &
                                      dface_cell,      &
                                      dface_vtx_idx,   &
                                      dface_vtx,       &
                                      dvtx_coord,      &
                                      dface_group_idx, &
                                      dface_group)

    use iso_c_binding
    implicit none

    type(c_ptr), value                 :: dcube
    integer(kind=pdm_g_num_s), pointer :: dface_cell(:)
    integer(kind=pdm_l_num_s), pointer :: dface_vtx_idx(:)
    integer(kind=pdm_g_num_s), pointer :: dface_vtx(:)
    double precision,          pointer :: dvtx_coord(:,:)
    integer(kind=pdm_l_num_s), pointer :: dface_group_idx(:)
    integer(kind=pdm_g_num_s), pointer :: dface_group(:)

    type(c_ptr)                        :: c_dface_cell
    type(c_ptr)                        :: c_dface_vtx_idx
    type(c_ptr)                        :: c_dface_vtx
    type(c_ptr)                        :: c_dvtx_coord
    type(c_ptr)                        :: c_dface_group_idx
    type(c_ptr)                        :: c_dface_group

    integer                            :: n_face_group
    integer                            :: dn_cell
    integer                            :: dn_face
    integer                            :: dn_vtx
    integer                            :: sface_vtx
    integer                            :: sface_group


    call PDM_dcube_gen_dim_get (dcube,           &
                                n_face_group,    &
                                dn_cell,         &
                                dn_face,         &
                                dn_vtx,          &
                                sface_vtx,       &
                                sface_group)

    call PDM_dcube_gen_data_get_cf (dcube,             &
                                    c_dface_cell,      &
                                    c_dface_vtx_idx,   &
                                    c_dface_vtx,       &
                                    c_dvtx_coord,      &
                                    c_dface_group_idx, &
                                    c_dface_group)

    call c_f_pointer(c_dface_cell, &
                     dface_cell,   &
                     [2*dn_face])

    call c_f_pointer(c_dface_vtx_idx, &
                     dface_vtx_idx,   &
                     [dn_face + 1])

    call c_f_pointer(c_dface_vtx, &
                     dface_vtx,   &
                     [sface_vtx])

    call c_f_pointer(c_dvtx_coord, &
                     dvtx_coord,   &
                     [3,dn_vtx])

    call c_f_pointer(c_dface_group_idx, &
                     dface_group_idx,   &
                     [n_face_group+1])

    call c_f_pointer(c_dface_group, &
                     dface_group,   &
                     [sface_group])

  end subroutine PDM_dcube_gen_data_get_

end module pdm_dcube_gen
