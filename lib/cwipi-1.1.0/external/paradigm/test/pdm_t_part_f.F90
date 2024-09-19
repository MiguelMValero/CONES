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
  use pdm_part
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

  type(c_ptr)                           :: dcube = C_NULL_PTR

  integer                               :: n_face_group = -1
  integer                               :: dn_cell      = -1
  integer                               :: dn_face      = -1
  integer                               :: dn_vtx       = -1
  integer                               :: sface_vtx    = -1
  integer                               :: sface_group  = -1

  integer (kind = pdm_g_num_s), pointer :: dface_cell(:)      => null()
  integer (kind = pdm_l_num_s), pointer :: dface_vtx_idx(:)   => null()
  integer (kind = pdm_g_num_s), pointer :: dface_vtx(:)       => null()
  double precision,             pointer :: dvtx_coord(:,:)    => null()
  integer (kind = pdm_l_num_s), pointer :: dface_group_idx(:) => null()
  integer (kind = pdm_g_num_s), pointer :: dface_group(:)     => null()

  type(c_ptr)                           :: ppart = C_NULL_PTR
  integer(c_int)                        :: method
  integer(c_int)                        :: n_part
  integer                               :: n_property_cell
  integer                               :: n_property_face
  integer(kind=PDM_l_num_s), pointer    :: renum_properties_cell(:) => null()
  integer(kind=PDM_l_num_s), pointer    :: renum_properties_face(:) => null()
  integer(kind=PDM_l_num_s), pointer    :: dcell_face_idx(:)        => null()
  integer(kind=PDM_g_num_s), pointer    :: dcell_face(:)            => null()
  integer(kind=PDM_l_num_s), pointer    :: dcell_tag(:)             => null()
  integer(kind=PDM_l_num_s), pointer    :: dcell_weight(:)          => null()
  integer(kind=PDM_l_num_s), pointer    :: dcell_part(:)            => null()
  integer(kind=PDM_l_num_s), pointer    :: dface_tag(:)             => null()
  integer(kind=PDM_l_num_s), pointer    :: dvtx_tag(:)              => null()
  logical                               :: have_dcell_part

  integer                               :: i_part
  integer                               :: n_cell
  integer                               :: n_face
  integer                               :: n_face_part_bound
  integer                               :: n_vtx
  integer                               :: n_proc
  integer                               :: n_total_part
  integer                               :: scell_face
  integer                               :: sface_vtx2
  integer                               :: sface_group2
  integer                               :: n_face_group2

  integer (kind = PDM_l_num_s), pointer :: cell_tag(:)                  => null()
  integer (kind = PDM_l_num_s), pointer :: cell_face_idx(:)             => null()
  integer (kind = PDM_l_num_s), pointer :: cell_face(:)                 => null()
  integer (kind = PDM_g_num_s), pointer :: cell_ln_to_gn(:)             => null()
  integer (kind = PDM_l_num_s), pointer :: face_tag(:)                  => null()
  integer (kind = PDM_l_num_s), pointer :: face_cell(:)                 => null()
  integer (kind = PDM_l_num_s), pointer :: face_vtx_idx(:)              => null()
  integer (kind = PDM_l_num_s), pointer :: face_vtx(:)                  => null()
  integer (kind = PDM_g_num_s), pointer :: face_ln_to_gn(:)             => null()
  integer (kind = PDM_l_num_s), pointer :: face_part_bound_proc_idx(:)  => null()
  integer (kind = PDM_l_num_s), pointer :: face_part_bound_part_idx(:)  => null()
  integer (kind = PDM_l_num_s), pointer :: face_part_bound(:)           => null()
  integer (kind = PDM_l_num_s), pointer :: vtx_tag(:)                   => null()
  double precision,             pointer :: vtx(:,:)                     => null()
  integer (kind = PDM_g_num_s), pointer :: vtx_ln_to_gn(:)              => null()
  integer (kind = PDM_l_num_s), pointer :: face_group_idx(:)            => null()
  integer (kind = PDM_l_num_s), pointer :: face_group(:)                => null()
  integer (kind = PDM_g_num_s), pointer :: face_group_ln_to_gn(:)       => null()

  double precision                      :: elapsed(4)
  double precision                      :: cpu(4)
  double precision                      :: cpu_user(4)
  double precision                      :: cpu_sys(4)

  integer(kind=PDM_l_num_s)             :: cells_average
  integer(kind=PDM_l_num_s)             :: cells_median
  double precision                      :: cells_std_deviation
  integer(kind=PDM_l_num_s)             :: cells_min
  integer(kind=PDM_l_num_s)             :: cells_max
  integer(kind=PDM_l_num_s)             :: bound_part_faces_average
  integer(kind=PDM_l_num_s)             :: bound_part_faces_median
  double precision                      :: bound_part_faces_std_deviation
  integer(kind=PDM_l_num_s)             :: bound_part_faces_min
  integer(kind=PDM_l_num_s)             :: bound_part_faces_max
  integer(kind=PDM_l_num_s)             :: bound_part_faces_sum

  integer                               :: code
  integer                               :: i_rank
  integer                               :: n_rank
  !-----------------------------------------------------------

  call mpi_init(code)
  call mpi_comm_rank(comm, i_rank, code)
  call mpi_comm_size(comm, n_rank, code)


  !  Generate a distributed mesh
  if (i_rank .eq. 0) then
    write(*, *) "-- Generate distributed mesh"
  end if
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

  call pdm_dcube_gen_data_get (dcube,           &
                               dface_cell,      &
                               dface_vtx_idx,   &
                               dface_vtx,       &
                               dvtx_coord,      &
                               dface_group_idx, &
                               dface_group)

  !  Split mesh
  if (i_rank .eq. 0) then
    write(*, *) "-- Split mesh"
  end if

  method = PDM_PART_SPLIT_HILBERT

  renum_properties_cell => null()
  renum_properties_face => null()
  n_property_cell = 0
  n_property_face = 0

  n_part = 1

  dcell_face_idx => null()
  dcell_face     => null()
  dcell_tag      => null()
  dcell_weight   => null()
  ! dcell_part     => null()
  allocate(dcell_part(dn_cell))
  dface_tag      => null()
  dvtx_tag       => null()

  have_dcell_part = .false.

  call pdm_part_create(ppart,                      &
                       comm,                       &
                       method,                     &
                       "PDM_PART_RENUM_CELL_NONE", &
                       "PDM_PART_RENUM_FACE_NONE", &
                       n_property_cell,            &
                       renum_properties_cell,      &
                       n_property_face,            &
                       renum_properties_face,      &
                       n_part,                     &
                       dn_cell,                    &
                       dn_face,                    &
                       dn_vtx,                     &
                       n_face_group,               &
                       dcell_face_idx,             &
                       dcell_face,                 &
                       dcell_tag,                  &
                       dcell_weight,               &
                       have_dcell_part,            &
                       dcell_part,                 &
                       dface_cell,                 &
                       dface_vtx_idx,              &
                       dface_vtx,                  &
                       dface_tag,                  &
                       dvtx_coord,                 &
                       dvtx_tag,                   &
                       dface_group_idx,            &
                       dface_group)

  deallocate(dcell_part)

  call pdm_part_time_get(ppart,    &
                         elapsed,  &
                         cpu,      &
                         cpu_user, &
                         cpu_sys)

  write (*,*) "elapsed  :", elapsed(1:4)
  write (*,*) "cpu      :", cpu(1:4)
  write (*,*) "cpu_user :", cpu_user(1:4)
  write (*,*) "cpu_sys  :", cpu_sys(1:4)

  call pdm_part_stat_get (ppart,                          &
                          cells_average,                  &
                          cells_median,                   &
                          cells_std_deviation,            &
                          cells_min,                      &
                          cells_max,                      &
                          bound_part_faces_average,       &
                          bound_part_faces_median,        &
                          bound_part_faces_std_deviation, &
                          bound_part_faces_min,           &
                          bound_part_faces_max,           &
                          bound_part_faces_sum)

  write (*,*) "cells avg, med, std, min, max :", &
                cells_average,                   &
                cells_median,                    &
                cells_std_deviation,             &
                cells_min,                       &
                cells_max

  !  Get split mesh
  do i_part = 1, n_part

    call pdm_part_part_dim_get(ppart,             &
                               i_part-1,          &
                               n_cell,            &
                               n_face,            &
                               n_face_part_bound, &
                               n_vtx,             &
                               n_proc,            &
                               n_total_part,      &
                               scell_face,        &
                               sface_vtx2,        &
                               sface_group2,      &
                               n_face_group2)

    call pdm_part_part_val_get (ppart,                    &
                                i_part-1,                 &
                                cell_tag,                 &
                                cell_face_idx,            &
                                cell_face,                &
                                cell_ln_to_gn,            &
                                face_tag,                 &
                                face_cell,                &
                                face_vtx_idx,             &
                                face_vtx,                 &
                                face_ln_to_gn,            &
                                face_part_bound_proc_idx, &
                                face_part_bound_part_idx, &
                                face_part_bound,          &
                                vtx_tag,                  &
                                vtx,                      &
                                vtx_ln_to_gn,             &
                                face_group_idx,           &
                                face_group,               &
                                face_group_ln_to_gn)

  enddo

  !  Free memory
  call pdm_part_free(ppart)
  call pdm_dcube_gen_free(dcube)

  !  If PDM_OWNERSHIP_USER
  ! call pdm_fortran_free_c(c_loc(dface_cell))
  ! call pdm_fortran_free_c(c_loc(dface_vtx_idx))
  ! call pdm_fortran_free_c(c_loc(dface_vtx))
  ! call pdm_fortran_free_c(c_loc(dvtx_coord))
  ! call pdm_fortran_free_c(c_loc(dface_group_idx))
  ! call pdm_fortran_free_c(c_loc(dface_group))

  if (i_rank .eq. 0) then
    write(*, *) "-- End"
  end if

  call mpi_finalize(code)

end program testf
