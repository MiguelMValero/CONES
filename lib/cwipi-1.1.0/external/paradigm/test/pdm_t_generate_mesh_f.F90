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
  use pdm_generate_mesh
  use iso_c_binding
  use pdm_fortran

  implicit none

#ifndef PDM_HAVE_FORTRAN_MPI_MODULE
  include "mpif.h"
#endif

  !-----------------------------------------------------------
  integer,                      parameter :: comm = MPI_COMM_WORLD
  integer                                 :: code
  integer                                 :: i_rank
  integer                                 :: n_rank
  integer(kind=pdm_g_num_s),    parameter :: n_vtx_seg = 10

  integer                                 :: n_vtx
  integer                                 :: n_elt
  double precision,               pointer :: coords(:,:)    => null()
  integer(kind=pdm_l_num_s),      pointer :: elt_vtx_idx(:) => null()
  integer(kind=pdm_l_num_s),      pointer :: elt_vtx(:)     => null()
  !-----------------------------------------------------------


  call mpi_init(code)
  call mpi_comm_rank(comm, i_rank, code)
  call mpi_comm_size(comm, n_rank, code)

  call PDM_generate_mesh_rectangle_simplified(comm,        &
                                              n_vtx_seg,   &
                                              n_vtx,       &
                                              n_elt,       &
                                              coords,      &
                                              elt_vtx_idx, &
                                              elt_vtx)

  ! write(*,*) "n_vtx       : ", n_vtx
  ! write(*,*) "n_elt       : ", n_elt
  ! write(*,*) "coords      : ", coords
  ! write(*,*) "elt_vtx_idx : ", elt_vtx_idx
  ! write(*,*) "elt_vtx     : ", elt_vtx

  call pdm_fortran_free_c(c_loc(coords))
  call pdm_fortran_free_c(c_loc(elt_vtx_idx))
  call pdm_fortran_free_c(c_loc(elt_vtx))

  call mpi_finalize(code)

end program testf
