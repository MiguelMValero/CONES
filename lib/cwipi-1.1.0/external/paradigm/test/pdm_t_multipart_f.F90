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
  use pdm_multipart
  use pdm_dcube_gen
  use iso_c_binding
  use pdm_fortran

  implicit none

#ifndef PDM_HAVE_FORTRAN_MPI_MODULE
  include "mpif.h"
#endif

  !-----------------------------------------------------------
  ! MPI
  integer,                  parameter   :: comm = MPI_COMM_WORLD
  integer                               :: code
  integer                               :: i_rank
  integer                               :: n_rank
  ! UTIL
  integer                               :: i
  ! MULTIPART
  type(c_ptr)                           :: multipart = C_NULL_PTR
  integer(c_int)                        :: split_method
  integer(c_int)                        :: n_part = 1
  integer(c_int)                        :: n_domain = 1
  integer(c_int)                        :: i_domain = -1
  integer(kind=PDM_l_num_s), pointer    :: n_part_domains(:)  => null()
  double precision,          pointer    :: part_fraction(:) => null()
  integer(kind=PDM_l_num_s), pointer    :: renum_cell_properties(:) => null()
  ! MESH
  type(c_ptr)                           :: dcube = C_NULL_PTR
  integer(pdm_g_num_s), parameter       :: n_vtx_seg = 10
  double precision,     parameter       :: length = 5.
  double precision,     parameter       :: zero_x = 1.
  double precision,     parameter       :: zero_y = 1.
  double precision,     parameter       :: zero_z = 1.
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
  integer                               :: n_jn = 0
  integer (kind = pdm_l_num_s), pointer :: dface_join_idx(:)  => null()
  ! integer (kind = pdm_l_num_s), pointer :: joins_glob_id(:)   => null()
  ! integer (kind = pdm_g_num_s), pointer :: dface_join(:)      => null()
  ! type(c_ptr)                           :: dm = C_NULL_PTR
  !-----------------------------------------------------------

  call mpi_init(code)
  call mpi_comm_rank(comm, i_rank, code)
  call mpi_comm_size(comm, n_rank, code)

  ! Initialize multipart

  split_method  = PDM_SPLIT_DUAL_WITH_HILBERT;

  allocate(n_part_domains(n_domain), &
           dface_join_idx(n_jn+1))

  do i = 1, n_domain
    n_part_domains(i) = n_part
  end do

  if (i_rank .eq. 0) then
    write(*, *) "PDM_multipart_create"
  end if

  call PDM_multipart_create(multipart, &
                            n_domain, &
                            n_part_domains, &
                            PDM_FALSE, &
                            split_method, &
                            PDM_PART_SIZE_HOMOGENEOUS, &
                            part_fraction, &
                            comm, &
                            PDM_OWNERSHIP_KEEP)

  ! Reordering options
  if (i_rank .eq. 0) then
    write(*, *) "PDM_multipart_set_reordering_options"
  end if

  call PDM_multipart_set_reordering_options(multipart, &
                                            i_domain, &
                                            "PDM_PART_RENUM_CELL_CUTHILL", &
                                            renum_cell_properties, &
                                            "PDM_PART_RENUM_FACE_LEXICOGRAPHIC")

  ! Generate Mesh (case : n_domain = 1)
  ! > dcube
  if (i_rank .eq. 0) then
    write(*, *) "> dcube"
  end if

  call pdm_dcube_gen_init(dcube,              &
                          comm,               &
                          n_vtx_seg,          &
                          length,             &
                          zero_x,             &
                          zero_y,             &
                          zero_z,             &
                          PDM_OWNERSHIP_KEEP)

  call pdm_dcube_gen_dim_get(dcube,           &
                             n_face_group,    &
                             dn_cell,         &
                             dn_face,         &
                             dn_vtx,          &
                             sface_vtx,       &
                             sface_group)

  call pdm_dcube_gen_data_get(dcube,           &
                              dface_cell,      &
                              dface_vtx_idx,   &
                              dface_vtx,       &
                              dvtx_coord,      &
                              dface_group_idx, &
                              dface_group)
  ! > dmesh
  if (i_rank .eq. 0) then
    write(*, *) "> dmesh"
  end if

  ! call PDM_dmesh_create(PDM_OWNERSHIP_KEEP, &
  !                       dn_cell,            &
  !                       dn_face,            &
  !                       -1,                 &
  !                       dn_vtx,             &
  !                       n_face_group,       &
  !                       n_jn,               &
  !                       comm)
  ! dface_join_idx(0) = 0
  ! call PDM_dmesh_set(dm,              &
  !                    dvtx_coord,      &
  !                    dface_vtx_idx,   &
  !                    dface_vtx,       &
  !                    dface_cell,      &
  !                    dface_group_idx, &
  !                    dface_group,     &
  !                    joins_glob_id,   &
  !                    dface_join_idx,  &
  !                    dface_join)

  ! call PDM_multipart_dmesh_set(multipart, i_domain, dm)
  ! call PDM_multipart_joins_set(multipart, n_total_joins, join_to_opposite)

  ! Run
  ! call PDM_multipart_compute(multipart)

  ! Get
  ! call PDM_multipart_part_dim_get(multipart, i_domain, i_part, n_section, n_elt, &
  !                                 n_cell, n_face, n_part_joins, n_vtx, n_proc, tn_part, &
  !                                 scell_face, sface_vtx, sface_bound, n_bounds, sface_join, n_joins)
  ! call PDM_multipart_part_val_get(multipart, i_domain, i_part, elt_vtx_idx, elt_vtx, elt_section_ln_to_gn, &
  !                                 cell_tag, cell_face_idx, cell_face, cell_ln_to_gn, &
  !                                 face_tag, face_cell, face_vtx_idx, face_vtx, face_ln_to_gn, &
  !                                 face_part_bound_proc_idx, face_part_bound_part_idx, face_part_bound, &
  !                                 vtx_tag, vtx, vtx_ln_to_gn, face_bound_idx, face_bound, &
  !                                 face_bound_ln_to_gn, face_join_idx, face_join, face_join_ln_to_gn)

  ! Free
  if (i_rank .eq. 0) then
    write(*, *) "PDM_multipart_free"
  end if

  call PDM_multipart_free(multipart)

  if (i_rank .eq. 0) then
    write(*, *) "-- End"
  end if

  call mpi_finalize(code)

end program testf
