!-----------------------------------------------------------------------------
! This file is part of the CWIPI library. 
!
! Copyright (C) 2018  ONERA
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

module mod_fvmc_ho_basis
  implicit none
  
  integer, parameter :: fvmc_edge       = 0
  integer, parameter :: fvmc_face_tria  = 1
  integer, parameter :: fvmc_face_quad  = 2
  integer, parameter :: fvmc_face_poly  = 3
  integer, parameter :: fvmc_cell_tetra = 5
  integer, parameter :: fvmc_cell_pyram = 6
  integer, parameter :: fvmc_cell_prism = 7
  integer, parameter :: fvmc_cell_hexa  = 8
  integer, parameter :: fvmc_cell_poly  = 9 

  interface
    subroutine fvmc_ho_basis  (type, order, n_nodes, n_pts, uvw, weights) &
         bind(c, name='FVMC_ho_basis')
      use iso_c_binding
      implicit none
      integer (C_INT), value :: type
      integer (C_INT), value :: order 
      integer (C_INT), value :: n_nodes
      integer (C_INT), value :: n_pts
    
      type (C_PTR),    value  :: uvw
      type (C_PTR),    value  :: weights
    end subroutine fvmc_ho_basis
  
    subroutine fvmc_ho_basis_free () bind(c, name='FVMC_ho_basis_free')
      use iso_c_binding
      implicit none
    end subroutine fvmc_ho_basis_free

  end interface
end module mod_fvmc_ho_basis
  
