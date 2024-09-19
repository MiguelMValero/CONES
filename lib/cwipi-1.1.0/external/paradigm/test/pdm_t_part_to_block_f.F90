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
  use pdm_pointer_array
  use pdm_part_to_block
  use iso_c_binding
  use pdm_fortran

  implicit none

#ifndef PDM_HAVE_FORTRAN_MPI_MODULE
  include "mpif.h"
#endif

  !-----------------------------------------------------------
  integer, parameter                    :: comm = MPI_COMM_WORLD

  type(c_ptr)                           :: ptb = C_NULL_PTR

  integer                               :: n_part
  type(PDM_pointer_array_t), pointer    :: gnum_elt => null()
  type(PDM_pointer_array_t), pointer    :: weight   => null()
  integer(pdm_l_num_s), pointer         :: n_elt(:)        => null()
  integer(pdm_g_num_s), pointer         :: elt_ln_to_gn(:) => null()

  type(PDM_pointer_array_t), pointer    :: part_stride => null()
  type(PDM_pointer_array_t), pointer    :: part_data   => null()
  integer(pdm_l_num_s), pointer         :: block_stride(:) => null()
  integer(pdm_l_num_s), pointer         :: block_data(:)   => null()

  integer(pdm_l_num_s), pointer      :: stride(:) => null()
  integer(pdm_l_num_s), pointer      :: data(:)   => null()

  integer                               :: code
  integer                               :: i_rank
  integer                               :: n_rank
  !-----------------------------------------------------------

  call mpi_init(code)
  call mpi_comm_rank(comm, i_rank, code)
  call mpi_comm_size(comm, n_rank, code)

  if (n_rank .ne. 1) then
    print *,'Error : 1 MPI processes are mandatory'
    call mpi_finalize(code)
  end if


  !  Define partitions
  n_part = 1

  allocate(n_elt(n_part))
  n_elt(1) = 5
  call PDM_pointer_array_create (gnum_elt,       &
                                 n_part,         &
                                 PDM_TYPE_G_NUM)

  allocate(elt_ln_to_gn(n_elt(1)))
  elt_ln_to_gn = [3, 5, 2, 4, 1]


  call PDM_pointer_array_part_set (gnum_elt,     &
                                   0,            &
                                   elt_ln_to_gn)

  !  Create Part-to-block object
  call PDM_part_to_block_create (ptb,                                &
                                 PDM_PART_TO_BLOCK_DISTRIB_ALL_PROC, & ! t_distrib
                                 PDM_PART_TO_BLOCK_POST_CLEANUP,     & ! t_post
                                 1.d0,                               & ! partActiveNode
                                 gnum_elt,                           &
                                 weight,                             &
                                 n_elt,                              &
                                 n_part,                             &
                                 comm)


  call PDM_pointer_array_create (part_stride,  &
                                 n_part,       &
                                 PDM_TYPE_INT)
  allocate(stride(n_elt(1)))
  stride(:) = 1

  call PDM_pointer_array_part_set (part_stride, &
                                   0,           &
                                   stride)


  call PDM_pointer_array_create (part_data,    &
                                 n_part,       &
                                 PDM_TYPE_INT)
  allocate(data(n_elt(1)))
  data = [13, 15, 12, 14, 11]

  call PDM_pointer_array_part_set (part_data, &
                                   0,         &
                                   data)


  call PDM_part_to_block_exch (ptb,                       &
                               PDM_STRIDE_VAR_INTERLACED, & ! t_stride
                               1,                         & ! cst_stride
                               part_stride,               &
                               part_data,                 &
                               block_stride,              &
                               block_data)

  print *, "block_stride :", block_stride
  print *, "block_data   :", block_data

  !  Free memory
  call PDM_part_to_block_free(ptb)

  deallocate(n_elt,  &
             stride, &
             data)
  call PDM_pointer_array_free(gnum_elt)
  call PDM_pointer_array_free(part_data)
  call PDM_pointer_array_free(part_stride)
  call pdm_fortran_free_c(c_loc(block_stride))
  call pdm_fortran_free_c(c_loc(block_data))

  if (i_rank .eq. 0) then
    write(*, *) "-- End"
  end if

  call mpi_finalize(code)

end program testf
