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
  use pdm_dcube_gen
  use iso_c_binding
  use pdm_fortran

  implicit none

#ifndef PDM_HAVE_FORTRAN_MPI_MODULE  
  include "mpif.h"
#endif  

  !-----------------------------------------------------------
  integer,              parameter       :: comm = MPI_COMM_WORLD
  integer(pdm_g_num_s), parameter       :: n_vtx_seg = 10
  double precision,     parameter       :: length = 5.
  double precision,     parameter       :: zero_x = 1.
  double precision,     parameter       :: zero_y = 1.
  double precision,     parameter       :: zero_z = 1.

  type(c_ptr)                           :: dcube

  integer                               :: n_face_group
  integer                               :: dn_cell
  integer                               :: dn_face
  integer                               :: dn_vtx
  integer                               :: sface_vtx
  integer                               :: sface_group

  integer (kind = pdm_g_num_s), pointer :: dface_cell(:)      => null()
  integer (kind = pdm_l_num_s), pointer :: dface_vtx_idx(:)   => null()
  integer (kind = pdm_g_num_s), pointer :: dface_vtx(:)       => null()
  double precision,             pointer :: dvtx_coord(:,:)    => null()
  integer (kind = pdm_l_num_s), pointer :: dface_group_idx(:) => null()
  integer (kind = pdm_g_num_s), pointer :: dface_group(:)     => null()

  integer                               :: code
  integer                               :: i_rank
  integer                               :: n_rank
  !-----------------------------------------------------------

  n_face_group = -1
  dn_cell      = -1
  dn_face      = -1
  dn_vtx       = -1
  sface_vtx    = -1
  sface_group  = -1

  call mpi_init(code)
  call mpi_comm_rank(comm, i_rank, code)
  call mpi_comm_size(comm, n_rank, code)

  call pdm_dcube_gen_init(dcube,              &
                          comm,               &
                          n_vtx_seg,          &
                          length,             &
                          zero_x,             &
                          zero_y,             &
                          zero_z,             &
                          PDM_OWNERSHIP_KEEP)

  call pdm_dcube_gen_dim_get (dcube,           &
                              n_face_group,    &
                              dn_cell,         &
                              dn_face,         &
                              dn_vtx,          &
                              sface_vtx,       &
                              sface_group)

  write(*,*) "n_face_group = ", n_face_group
  write(*,*) "dn_cell      = ", dn_cell
  write(*,*) "dn_face      = ", dn_face
  write(*,*) "dn_vtx       = ", dn_vtx
  write(*,*) "sface_vtx    = ", sface_vtx
  write(*,*) "sface_group  = ", sface_group

  call pdm_dcube_gen_data_get (dcube,           &
                               dface_cell,      &
                               dface_vtx_idx,   &
                               dface_vtx,       &
                               dvtx_coord,      &
                               dface_group_idx, &
                               dface_group)


  write(*,*) dface_cell(1), dface_cell(2)

  write(*,*) dface_vtx_idx(1), dface_vtx_idx(2)
  write(*,*) dface_vtx(1), dface_vtx(2)
  write(*,*) dvtx_coord(1,1), dvtx_coord(2,1)

  write(*,*) "dface_group_idx(sface_group+1)  = ", dface_group_idx(n_face_group+1)


  !  Free memory
  call pdm_dcube_gen_free(dcube)


  !  If PDM_OWNERSHIP_USER
  ! call pdm_fortran_free_c(c_loc(dface_cell))
  ! call pdm_fortran_free_c(c_loc(dface_vtx_idx))
  ! call pdm_fortran_free_c(c_loc(dface_vtx))
  ! call pdm_fortran_free_c(c_loc(dvtx_coord))
  ! call pdm_fortran_free_c(c_loc(dface_group_idx))
  ! call pdm_fortran_free_c(c_loc(dface_group))

  call mpi_finalize(code)

end program testf
