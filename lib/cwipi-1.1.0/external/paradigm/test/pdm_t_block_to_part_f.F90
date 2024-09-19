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
  use pdm_block_to_part
  use iso_c_binding
  use pdm_fortran

  implicit none

#ifndef PDM_HAVE_FORTRAN_MPI_MODULE  
  include "mpif.h"
#endif  

  !-----------------------------------------------------------
  integer, parameter                    :: comm = MPI_COMM_WORLD

  logical                               :: exch_in_place = .false.

  type(c_ptr)                           :: btp = C_NULL_PTR

  integer(pdm_g_num_s), pointer         :: block_distrib_idx(:) => null()
  integer                               :: n_part
  type(PDM_pointer_array_t), pointer    :: gnum_elt        => null()
  integer(pdm_l_num_s), pointer         :: n_elt(:)        => null()
  integer(pdm_g_num_s), pointer         :: elt_ln_to_gn(:) => null()

  integer(pdm_l_num_s), pointer         :: block_stride(:) => null()
  double precision,     pointer         :: block_data(:)   => null()

  type(PDM_pointer_array_t), pointer    :: part_stride => null()
  type(PDM_pointer_array_t), pointer    :: part_data   => null()

  integer(pdm_l_num_s), pointer         :: stride(:) => null()
  double precision,     pointer         :: data(:)   => null()


  integer                               :: code
  integer                               :: i_rank
  integer                               :: n_rank
  !-----------------------------------------------------------

  call read_args(exch_in_place)

  print *, "exch_in_place :", exch_in_place

  call mpi_init(code)
  call mpi_comm_rank(comm, i_rank, code)
  call mpi_comm_size(comm, n_rank, code)

  if (n_rank .ne. 1) then
    print *,'Error : 1 MPI processes are mandatory'
    call mpi_finalize(code)
    stop
  end if


  !  Define partitions and block distribution
  allocate(block_distrib_idx(2))
  block_distrib_idx(1) = 0
  block_distrib_idx(2) = 5

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

  !  Create Block-to-part object
  call PDM_block_to_part_create (btp,               &
                                 block_distrib_idx, &
                                 gnum_elt,          &
                                 n_elt,             &
                                 n_part,            &
                                 comm)



  allocate(block_stride(n_elt(1)))
  block_stride(:) = 1

  allocate(block_data(n_elt(1)))
  block_data = [10.d0, 20.d0, 30.d0, 40.d0, 50.d0]

  if (exch_in_place) then
    call PDM_pointer_array_create (part_stride,  &
                                   n_part,       &
                                   PDM_TYPE_INT)
    allocate(stride(n_elt(1)))
    call PDM_pointer_array_part_set (part_stride, &
                                     0,           &
                                     stride)

    call PDM_pointer_array_create (part_data,       &
                                   n_part,          &
                                   PDM_TYPE_DOUBLE)
    allocate(data(n_elt(1)))
    call PDM_pointer_array_part_set (part_data, &
                                     0,         &
                                     data)

    call PDM_block_to_part_exch_in_place (btp,                       &
                                          PDM_STRIDE_VAR_INTERLACED, & ! t_stride
                                          block_stride,              &
                                          block_data,                &
                                          part_stride,               &
                                          part_data)
  else
    call PDM_block_to_part_exch (btp,                       &
                                 PDM_STRIDE_VAR_INTERLACED, & ! t_stride
                                 block_stride,              &
                                 block_data,                &
                                 part_stride,               &
                                 part_data)
  end if

  call PDM_block_to_part_free (btp)


  !  Check part data
  if (.not.exch_in_place) then
    call PDM_pointer_array_part_get (part_stride, &
                                     0,           &
                                     stride)

    call PDM_pointer_array_part_get (part_data,   &
                                     0,           &
                                     data)
  end if

  print *, "part_stride :", stride
  print *, "part_data   :", data


  !  Free memory
  deallocate(block_distrib_idx)
  deallocate(n_elt)
  deallocate(block_stride)
  deallocate(block_data)

  call PDM_pointer_array_free(gnum_elt)
  call PDM_pointer_array_free(part_data)
  call PDM_pointer_array_free(part_stride)

  if (i_rank .eq. 0) then
    write(*, *) "-- End"
  end if

  call mpi_finalize(code)


contains

  subroutine usage(code)

    implicit none

    integer, intent(in) :: code

    write(*,*) "Usage :"
    write(*,*) " -inplace     Perform 'in place' exchange."
    write(*,*) " -h           This message."

    stop

  end subroutine usage

  subroutine read_args(exch_in_place)

    implicit none

    logical,              intent(inout) :: exch_in_place
    integer                             :: argc, i, error
    character(999)                      :: arg

    argc = command_argument_count()

    i = 1
    do while (i <= argc)
      call get_command_argument(i, arg, status=error)
      if (error .ne. 0) then
        call usage(error)
      endif
      select case(arg)

        case ('-h')
          call usage(0)

        case ('-inplace')
          exch_in_place = .true.

      end select

      i = i + 1
    end do

  end subroutine read_args


end program testf
