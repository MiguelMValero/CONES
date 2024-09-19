!-----------------------------------------------------------------------------
! This file is part of the ParaDiGM library.
!
! Copyright (C) 2022  ONERA
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

  use pdm_gnum_location
  use iso_c_binding
  use pdm_fortran

  implicit none

#ifndef PDM_HAVE_FORTRAN_MPI_MODULE
  include "mpif.h"
#endif

  type my_type

    integer                       :: n_elt
    integer(pdm_g_num_s), pointer :: elt_gnum(:)       => null()

    integer                       :: n_requested
    integer(pdm_g_num_s), pointer :: requested_gnum(:) => null()

  end type my_type


  !-----------------------------------------------------------
  integer, parameter                    :: comm = MPI_COMM_WORLD

  type(c_ptr)                           :: gloc = C_NULL_PTR

  integer                               :: n_part1
  integer                               :: n_part2

  type(my_type), allocatable            :: my_part1(:)
  type(my_type), allocatable            :: my_part2(:)

  integer(pdm_l_num_s), pointer         :: location_idx(:) => null()
  integer(pdm_l_num_s), pointer         :: location(:)     => null()

  integer                               :: code
  integer                               :: i_rank
  integer                               :: n_rank
  integer                               :: i, j
  character                             :: strnum
  integer                                ::debug = 0
  integer                               :: fid = 13
  !-----------------------------------------------------------

  call mpi_init(code)
  call mpi_comm_rank(comm, i_rank, code)
  call mpi_comm_size(comm, n_rank, code)

  if (n_rank .ne. 2) then
    write(*,*) 'Error : 2 MPI processes are mandatory'
    call mpi_finalize(code)
    stop
  end if

  if (1 .eq. debug) then
    write (strnum, '(i1)') i_rank
    open(unit=fid, file="gnum_location_"//strnum//".log", action='write')
  end if   

  !  Define partitions
  n_part1 = 3
  n_part2 = 3

  allocate(my_part1(n_part1), my_part2(n_part2))
  if (i_rank .eq. debug) then

    !! Part 1
    my_part1(1)%n_requested = 0
    allocate(my_part1(1)%requested_gnum(my_part1(1)%n_requested))
    !my_part1(1)%requested_gnum = []

    my_part1(2)%n_requested = 6
    allocate(my_part1(2)%requested_gnum(my_part1(2)%n_requested))
    my_part1(2)%requested_gnum = [9, 11, 9, 10, 11, 10]

    my_part1(3)%n_requested = 1
    allocate(my_part1(3)%requested_gnum(my_part1(3)%n_requested))
    my_part1(3)%requested_gnum = [9]


    !! Part 2
    my_part2(1)%n_elt = 1
    allocate(my_part2(1)%elt_gnum(my_part2(1)%n_elt))
    my_part2(1)%elt_gnum = [9]

    my_part2(2)%n_elt = 4
    allocate(my_part2(2)%elt_gnum(my_part2(2)%n_elt))
    my_part2(2)%elt_gnum = [15, 13, 5, 2]

    my_part2(3)%n_elt = 2
    allocate(my_part2(3)%elt_gnum(my_part2(3)%n_elt))
    my_part2(3)%elt_gnum = [3, 8]

  else

    !! Part 1
    my_part1(1)%n_requested = 0
    allocate(my_part1(1)%requested_gnum(my_part1(1)%n_requested))
    !my_part1(1)%requested_gnum = []

    my_part1(2)%n_requested = 8
    allocate(my_part1(2)%requested_gnum(my_part1(2)%n_requested))
    my_part1(2)%requested_gnum = [12, 11, 12, 10, 16, 11, 16, 10]

    my_part1(3)%n_requested = 2
    allocate(my_part1(3)%requested_gnum(my_part1(3)%n_requested))
    my_part1(3)%requested_gnum = [12, 16]


    !! Part 2
    my_part2(1)%n_elt = 2
    allocate(my_part2(1)%elt_gnum(my_part2(1)%n_elt))
    my_part2(1)%elt_gnum = [12, 16]

    my_part2(2)%n_elt = 6
    allocate(my_part2(2)%elt_gnum(my_part2(2)%n_elt))
    my_part2(2)%elt_gnum = [11, 10, 17, 14, 1, 4]

    my_part2(3)%n_elt = 2
    allocate(my_part2(3)%elt_gnum(my_part2(3)%n_elt))
    my_part2(3)%elt_gnum = [7, 6]

  end if


  !  Create Gnum location object
  call pdm_gnum_location_create (gloc,               &
                                 n_part2,            &
                                 n_part1,            &
                                 comm,               &
                                 PDM_OWNERSHIP_KEEP)


  !  Set elements
  if (1 .eq. debug) then
    write (fid, *) "---- Elements ----"
  endif
  do i = 1, n_part2
    call pdm_gnum_location_elements_set (gloc,                 &
                                         i-1,                  &
                                         my_part2(i)%n_elt,    &
                                         my_part2(i)%elt_gnum)

    if (1 .eq. debug) then
      write (fid, *) "part", i-1
      write (fid, *) "lnum -> gnum :"
      do j = 1, my_part2(i)%n_elt
        write (fid, *) j, " ->", my_part2(i)%elt_gnum(j)
      end do
    endif
  end do


  !  Set requested elements
  do i = 1, n_part1
    call pdm_gnum_location_requested_elements_set (gloc,                      &
                                                   i-1,                       &
                                                   my_part1(i)%n_requested,   &
                                                   my_part1(i)%requested_gnum)
  end do


  !  Compute location
  call pdm_gnum_location_compute (gloc)


  !  Get location
  if (1 .eq. debug) then
    write (fid, *) ""
    write (fid, *) ""
    write (fid, *) "---- Location ----"
  endif
  do i = 1, n_part1

    call pdm_gnum_location_get (gloc,         &
                                i-1,          &
                                location_idx, &
                                location)

    if (1 .eq. debug) then
      write (fid, *) "part", i-1
      write (fid, *) "location_idx :", location_idx
      write (fid, *) "requested gnum -> location (rank, part, lnum) :"
      do j = 1, my_part1(i)%n_requested
       write (fid, *) my_part1(i)%requested_gnum(j), " ->", location(location_idx(j)+1:location_idx(j+1))
      end do
    endif

  end do


  !  Free memory
  call pdm_gnum_location_free (gloc)


  do i = 1, n_part1
    call free_my_type(my_part1(i))
  end do
  deallocate(my_part1)

  do i = 1, n_part2
    call free_my_type(my_part2(i))
  end do
  deallocate(my_part2)


  if (1 .eq. debug) then
   close(fid)
  endif


  if (i_rank .eq. debug) then
    write(*, *) "-- End"
  end if

  call mpi_finalize(code)


  contains





  subroutine free_my_type(p)

    implicit none

    type(my_type) :: p

    if (associated(p%elt_gnum)) then
      deallocate(p%elt_gnum)
    end if

    if (associated(p%requested_gnum)) then
      deallocate(p%requested_gnum)
    end if

  end subroutine free_my_type


end program testf
