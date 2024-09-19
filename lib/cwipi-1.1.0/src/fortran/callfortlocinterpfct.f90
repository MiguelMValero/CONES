!-----------------------------------------------------------------------------
! This file is part of the CWIPI library. 
!
! Copyright (C) 2011-2022  ONERA
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

subroutine callfortlocinterpfct(c_interface_type, &
                                c_code_name, &
                                c_l_code_name, &
                                c_src_n_block, &
                                c_src_blocks_type, &
                                c_src_i_part, &
                                c_src_n_vtx, &
                                c_src_vtx_coords, &
                                c_src_vtx_global_num, &
                                c_src_n_elts, &
                                c_src_id_block, &
                                c_src_elt_in_block, &
                                c_src_elt_vtx_idx, &
                                c_src_elt_vtx, &
                                c_src_elts_global_num, &
                                c_tgt_n_pts, &
                                c_tgt_pts_elt_idx, &
                                c_tgt_pts_coords, &
                                c_tgt_pts_dist, &
                                c_tgt_pts_uvw, &
                                c_tgt_pts_weights_idx, &
                                c_tgt_pts_weights, &
                                c_stride, &
                                c_src_field_dof_location, &
                                src_field, &
                                tgt_field, &
                                ptFortranLocInterpolationFct)

  use iso_c_binding
  implicit none

  interface
     subroutine  ptFortranLocInterpolationFct(interface_type, &
                                              code_name, &
                                              src_n_block, &
                                              src_blocks_type, &
                                              src_i_part, &
                                              src_n_vtx, &
                                              src_vtx_coords, &
                                              src_vtx_global_num, &
                                              src_n_elts, &
                                              src_id_block, &
                                              src_elt_in_block, &
                                              src_elt_vtx_idx, &
                                              src_elt_vtx, &
                                              src_elts_global_num, &
                                              tgt_n_pts, &
                                              tgt_pts_elt_idx, &
                                              tgt_pts_coords, &
                                              tgt_pts_dist, &
                                              tgt_pts_uvw, &
                                              tgt_pts_weights_idx, &
                                              tgt_pts_weights, &
                                              stride, &
                                              src_field_dof_location, &
                                              src_field, &
                                              tgt_field)


      use iso_c_binding
      integer(c_int), intent(in)  :: interface_type
      character (len=*)           :: code_name
      integer(c_int), intent(in)  :: src_n_block
      integer(c_int), pointer     :: src_blocks_type(:)
      integer(c_int), intent(in)  :: src_i_part
      integer(c_int), intent(in)  :: src_n_vtx
      real(kind = 8), pointer     :: src_vtx_coords(:)
      integer(c_long), pointer    :: src_vtx_global_num(:)
      integer(c_int), intent(in)  :: src_n_elts
      integer(c_int), pointer     :: src_id_block(:)
      integer(c_int), pointer     :: src_elt_in_block(:)
      integer(c_int), pointer     :: src_elt_vtx_idx(:)
      integer(c_int), pointer     :: src_elt_vtx(:)
      integer(c_long), pointer    :: src_elts_global_num(:)
      integer(c_int), intent(in)  :: tgt_n_pts
      integer(c_int), pointer     :: tgt_pts_elt_idx(:)
      real(kind = 8), pointer     :: tgt_pts_coords(:)
      real(kind = 8), pointer     :: tgt_pts_dist(:)
      real(kind = 8), pointer     :: tgt_pts_uvw(:)
      integer(c_int), pointer     :: tgt_pts_weights_idx(:)
      real(kind = 8), pointer     :: tgt_pts_weights(:)
      integer(c_int), intent(in)  :: stride
      integer(c_int), intent(in)  :: src_field_dof_location
      type(c_ptr)                 :: src_field
      type(c_ptr)                 :: tgt_field

     end subroutine ptFortranLocInterpolationFct
  end interface


  integer(c_int), value       :: c_interface_type
  type(c_ptr), value          :: c_code_name
  integer(c_int), value       :: c_l_code_name
  integer(c_int), value       :: c_src_n_block
  type(c_ptr), value          :: c_src_blocks_type
  integer(c_int), value       :: c_src_i_part
  integer(c_int), value       :: c_src_n_vtx
  type(c_ptr), value          :: c_src_vtx_coords
  type(c_ptr), value          :: c_src_vtx_global_num
  integer(c_int), value       :: c_src_n_elts
  type(c_ptr), value          :: c_src_id_block
  type(c_ptr), value          :: c_src_elt_in_block
  type(c_ptr), value          :: c_src_elt_vtx_idx
  type(c_ptr), value          :: c_src_elt_vtx
  type(c_ptr), value          :: c_src_elts_global_num
  integer(c_int), value       :: c_tgt_n_pts
  type(c_ptr), value          :: c_tgt_pts_elt_idx
  type(c_ptr), value          :: c_tgt_pts_coords
  type(c_ptr), value          :: c_tgt_pts_dist
  type(c_ptr), value          :: c_tgt_pts_uvw
  type(c_ptr), value          :: c_tgt_pts_weights_idx
  type(c_ptr), value          :: c_tgt_pts_weights
  integer(c_int), value       :: c_stride
  integer(c_int), value       :: c_src_field_dof_location
  type(c_ptr)                 :: src_field
  type(c_ptr)                 :: tgt_field

  character (len=512)         :: code_name
  character, pointer          :: f_code_name(:)
  integer(c_int), pointer     :: src_blocks_type(:)
  real(kind = 8), pointer     :: src_vtx_coords(:)
  integer(c_long), pointer    :: src_vtx_global_num(:)
  integer(c_int), pointer     :: src_id_block(:)
  integer(c_int), pointer     :: src_elt_in_block(:)
  integer(c_int), pointer     :: src_elt_vtx_idx(:)
  integer(c_int), pointer     :: src_elt_vtx(:)
  integer(c_long), pointer    :: src_elts_global_num(:)
  integer(c_int), pointer     :: tgt_pts_elt_idx(:)
  real(kind = 8), pointer     :: tgt_pts_coords(:)
  real(kind = 8), pointer     :: tgt_pts_dist(:)
  real(kind = 8), pointer     :: tgt_pts_uvw(:)
  integer(c_int), pointer     :: tgt_pts_weights_idx(:)
  real(kind = 8), pointer     :: tgt_pts_weights(:)

  integer                     :: i   

  if (c_l_code_name .gt. 512) then
    print *, "callfortlocinterpfct : l_code_name > 512"
    call exit
  endif

  call c_f_pointer(c_code_name, f_code_name, [512])

  do i = 1, c_l_code_name
    code_name(i:i) = f_code_name(i) 
  enddo

  do i = c_l_code_name, 512
    code_name(i:i) = ' ' 
  enddo

  call c_f_pointer(c_src_blocks_type, src_blocks_type, [c_src_n_block])
  call c_f_pointer(c_src_vtx_coords, src_vtx_coords, [3*c_src_n_vtx])
  call c_f_pointer(c_src_vtx_global_num, src_vtx_global_num, [c_src_n_vtx])
  call c_f_pointer(c_src_id_block, src_id_block, [c_src_n_elts])
  call c_f_pointer(c_src_elt_in_block, src_elt_in_block, [c_src_n_elts])
  call c_f_pointer(c_src_elt_vtx_idx, src_elt_vtx_idx, [c_src_n_elts + 1])
  call c_f_pointer(c_src_elt_vtx, src_elt_vtx, [src_elt_vtx_idx(c_src_n_elts + 1)])
  call c_f_pointer(c_src_elts_global_num, src_elts_global_num, [c_src_n_elts])
  call c_f_pointer(c_tgt_pts_elt_idx, tgt_pts_elt_idx, [c_tgt_n_pts])
  call c_f_pointer(c_tgt_pts_coords, tgt_pts_coords, [3 * c_tgt_n_pts])
  call c_f_pointer(c_tgt_pts_dist, tgt_pts_dist, [c_tgt_n_pts])
  call c_f_pointer(c_tgt_pts_uvw, tgt_pts_uvw, [3*c_tgt_n_pts])
  call c_f_pointer(c_tgt_pts_weights_idx, tgt_pts_weights_idx, [c_tgt_n_pts + 1])
  call c_f_pointer(c_tgt_pts_weights, tgt_pts_weights, [tgt_pts_weights_idx(c_tgt_n_pts + 1)])

  call ptFortranLocInterpolationFct (c_interface_type, &
                                     code_name, &
                                     c_src_n_block, &
                                     src_blocks_type, &
                                     c_src_i_part, &
                                     c_src_n_vtx, &
                                     src_vtx_coords, &
                                     src_vtx_global_num, &
                                     c_src_n_elts, &
                                     src_id_block, &
                                     src_elt_in_block, &
                                     src_elt_vtx_idx, &
                                     src_elt_vtx, &
                                     src_elts_global_num, &
                                     c_tgt_n_pts, &
                                     tgt_pts_elt_idx, &
                                     tgt_pts_coords, &
                                     tgt_pts_dist, &
                                     tgt_pts_uvw, &
                                     tgt_pts_weights_idx, &
                                     tgt_pts_weights, &
                                     c_stride, &
                                     c_src_field_dof_location, &
                                     src_field, &
                                     tgt_field)
                                     

end subroutine callfortlocinterpfct

