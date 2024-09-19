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

program testf

  use pdm
#ifdef PDM_HAVE_FORTRAN_MPI_MODULE
  use mpi
#endif
  use pdm_mesh_check
  use iso_c_binding
  use pdm_fortran

  implicit none

#ifndef PDM_HAVE_FORTRAN_MPI_MODULE
  include "mpif.h"
#endif

  !-----------------------------------------------------------
  integer(pdm_g_num_s)          :: n_vtx
  integer(pdm_g_num_s)          :: l_face_vtx
  integer(pdm_g_num_s), pointer :: face_vtx(:) => null()
  double precision,     pointer :: coords(:) => null()
  integer(pdm_l_num_s)          :: n_holes

  integer                       :: i, j, k
  integer                       :: code
  !-----------------------------------------------------------

  !  Init
  call mpi_init(code)

  n_vtx      = 9
  l_face_vtx = 8


  allocate(face_vtx(l_face_vtx))
  face_vtx(1) = 1
  face_vtx(2) = 2
  face_vtx(3) = 5
  face_vtx(4) = 4
  ! face_vtx(5) = 2
  ! face_vtx(6) = 3
  ! face_vtx(7) = 6
  ! face_vtx(8) = 5
  face_vtx(5) = 4
  face_vtx(6) = 5
  face_vtx(7) = 8
  face_vtx(8) = 7

  allocate(coords(3*n_vtx))
  k = 1
  outer_loop : do j = 1, 3
    do i = 1, 3
      coords(3*(k-1)+1) = i
      coords(3*(k-1)+2) = j
      coords(3*(k-1)+3) = 0.d0
      k = k+1
      if (k >= n_vtx) then
        exit outer_loop
      endif
    enddo
  enddo outer_loop



  call PDM_mesh_check_unconnected_vertex (n_vtx,      &
                                          l_face_vtx, &
                                          face_vtx,   &
                                          coords,     &
                                          n_holes)

  print *, n_vtx, "vertices, ", n_holes, "holes"

  deallocate(face_vtx)
  deallocate(coords)

  call mpi_finalize(code)

end program testf
