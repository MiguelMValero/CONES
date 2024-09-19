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
  use PDM_global_reduce
  use iso_c_binding
  use pdm_fortran

  implicit none

#ifndef PDM_HAVE_FORTRAN_MPI_MODULE
  include "mpif.h"
#endif

  !-----------------------------------------------------------

  integer :: code
  integer :: i_rank
  integer :: n_rank

  integer(c_int), parameter :: fComm  = MPI_COMM_WORLD
  integer(c_int), parameter :: n_part = 1
  integer(c_int), parameter :: n_pts  = 6

  integer(kind = pdm_g_num_s), pointer :: pts_ln_to_gn(:)      => null()
  double precision,            pointer :: pts_local_field(:)   => null()
  double precision,            pointer :: pts_reduced_field(:) => null()

  type(c_ptr) :: gre

  integer(kind = pdm_g_num_s) :: i, j, k
  integer                     :: n_error

  !-----------------------------------------------------------



  !*
  ! Init
  !*
  call mpi_init(code)
  call mpi_comm_rank(fComm, i_rank, code)
  call mpi_comm_size(fComm, n_rank, code)

  if (n_rank .ne. 2) then
    print *,'Error : 2 MPI processes are mandatory'
    call mpi_finalize(code)
    stop
  end if

  !*
  !  Create 'partitions'
  !*
  allocate(pts_ln_to_gn     (  n_pts))
  allocate(pts_local_field  (2*n_pts))
  allocate(pts_reduced_field(2*n_pts))
  k = 0
  do j = 0,2
    do i = 0,1
      k = k + 1
      pts_ln_to_gn(k) = 3*j + i + 1
      pts_local_field(2*k - 1) = i
      pts_local_field(2*k    ) = j
    end do
  end do

  if (i_rank .eq. 1) then
    pts_ln_to_gn    = pts_ln_to_gn    + 1
    pts_local_field = pts_local_field + 1
  end if

  write(*, *) "rank #", i_rank, "local :", pts_local_field(:)


  !*
  !  Create global reduction object
  !*
  call PDM_global_reduce_create(gre, n_part, fComm)

  call PDM_global_reduce_g_num_set(gre,          & !
                                   0,            & ! i_part
                                   n_pts,        & !
                                   pts_ln_to_gn)

  call PDM_global_reduce_field_set(gre,              & !
                                   0,                & ! i_part
                                   2,                & ! stride
                                   pts_local_field,  & !
                                   pts_reduced_field)

  ! 0: MIN
  ! 1: MAX
  ! 2: SUM
  call PDM_global_reduce_operation_set(gre, & !
                                       0)     ! operation

  !*
  !  Compute global min
  !*
  call PDM_global_reduce_field_compute(gre)

  write(*, *) "rank #", i_rank, "reduced :", pts_reduced_field(:)


  !*
  !  Check result
  !*
  if (i_rank .eq. 1) then
    pts_reduced_field = pts_reduced_field + 1
  end if

  n_error = 0
  do k = 1, n_pts

    j = (pts_ln_to_gn(k) - 1) / 3
    i = pts_ln_to_gn(k) - 1 - 3*j

    if (abs(pts_local_field(2*k - 1) - i) > 1.d-6 .and. &
        abs(pts_local_field(2*k    ) - j) > 1.d-6) then
      write(*,*) "Error for point ", k
    end if
  end do

  write(*, *) "rank #", i_rank, ":", n_error, "error(s) /", n_pts

  !*
  !  Free memory
  !*
  deallocate(pts_ln_to_gn)
  deallocate(pts_local_field)
  deallocate(pts_reduced_field)

  call PDM_global_reduce_free(gre)

  if (i_rank .eq. 0) then
    write(*, *) "-- End"
  end if

  call mpi_finalize(code)

end program testf
