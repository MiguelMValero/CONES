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
  use pdm_closest_points
  use iso_c_binding


  implicit none

#ifndef PDM_HAVE_FORTRAN_MPI_MODULE  
  include "mpif.h"
#endif  

  ! integer (kind = pdm_g_num_s), parameter :: n_g_points_src = 10
  ! integer (kind = pdm_g_num_s), parameter :: n_g_points_tgt = 10

  integer, parameter :: n_part_cloud_src = 1
  integer, parameter :: n_part_cloud_tgt = 1

  integer, parameter :: n_local_points_src = 5
  integer, parameter :: n_local_points_tgt = 5

  integer, parameter :: n_closest = 2

  double precision, pointer :: coords_src(:,:) ! pointer or allocatble, target
  double precision, pointer :: coords_tgt(:,:) ! pointer or allocatble, target

  integer :: code
  integer :: i_rank
  integer :: n_rank

  integer (kind = pdm_g_num_s), pointer :: gnum_src(:) ! pointer or allocatble, target
  integer (kind = pdm_g_num_s), pointer :: gnum_tgt(:) ! pointer or allocatble, target

  integer (kind = pdm_g_num_s), pointer :: closest_src_gnum(:)
  double precision,             pointer :: closest_src_distance(:)

  integer :: i

  type(c_ptr)              :: cls

  call mpi_init(code)
  call mpi_comm_rank(mpi_comm_world, i_rank, code)
  call mpi_comm_size(mpi_comm_world, n_rank, code)

  if (n_rank .ne. 2) then
    print *,'Error : 2 MPI processes are mandatory'
    call mpi_finalize(code)
    stop
  end if

  allocate (coords_src(3,n_local_points_src))
  allocate (coords_tgt(3,n_local_points_tgt))
  allocate (gnum_src(n_local_points_src))
  allocate (gnum_tgt(n_local_points_tgt))

  if (i_rank .eq. 0) then
    do i = 1, n_local_points_src
      coords_src(1,i) = 0. + i - 1
      coords_src(2,i) = 0. + i - 1
      coords_src(3,i) = 0. + i - 1
      gnum_src(i) = 2*(i-1) + 1
    end do
    do i = 1, n_local_points_tgt
      coords_tgt(1,i) = n_local_points_tgt + i - 1
      coords_tgt(2,i) = n_local_points_tgt + i - 1
      coords_tgt(3,i) = n_local_points_tgt + i - 1
      gnum_tgt(i) = 2*(i-1) + 2
    end do
  else
    do i = 1, n_local_points_src
      coords_src(1,i) = n_local_points_tgt + i - 1
      coords_src(2,i) = n_local_points_tgt + i - 1
      coords_src(3,i) = n_local_points_tgt + i - 1
      gnum_src(i) = 2*(i-1) + 2
    end do
    do i = 1, n_local_points_tgt
      coords_tgt(1,i) = 0. + i - 1
      coords_tgt(2,i) = 0. + i - 1
      coords_tgt(3,i) = 0. + i - 1
      gnum_tgt(i) = 2*(i-1) + 1
    end do
  endif


  do i = 1, n_local_points_tgt
    coords_tgt(1,i) = coords_tgt(1,i) / 10.
    coords_tgt(2,i) = coords_tgt(2,i) / 10.
    coords_tgt(3,i) = coords_tgt(3,i) / 10.
  enddo

  do i = 1, n_local_points_src
    coords_src(1,i) = coords_src(1,i) / 10.
    coords_src(2,i) = coords_src(2,i) / 10.
    coords_src(3,i) = coords_src(3,i) / 10.
  enddo

  !
  ! Create a new PDM_closest_points structure
  !   The MPI communicator and the 'n' closest neighbors are setted
  !


  call PDM_closest_points_create (cls,                &
                                  MPI_COMM_WORLD,     &
                                  n_closest,          &
                                  PDM_OWNERSHIP_KEEP)


  !
  ! Set the number of local partition of the source point cloud and of the target point cloud
  !

  call PDM_closest_points_n_part_cloud_set (cls, &
                                            n_part_cloud_src, &
                                            n_part_cloud_tgt)

  !
  ! Set the point coordinates and the global numbering for any partition of the source
  !

  if (n_part_cloud_src .ne. 1) then
    print *, "For this test, n_part_cloud_src must be equal to 1"
    stop
  end if

  do i = 1, n_part_cloud_src
    call PDM_closest_points_src_cloud_set (cls, &
                                           i-1, & !!! ipart : 0 -> n_part-1 !!!
                                           n_local_points_src, &
                                           coords_src, &
                                           gnum_src)
  end do

  !
  ! Set the point coordinates and the global numbering for any partition of the target
  !

  if (n_part_cloud_tgt .ne. 1) then
    print *, "For this test, n_part_cloud_tgt must be equal to 1"
    stop
  end if

  do i = 1, n_part_cloud_tgt
    call PDM_closest_points_tgt_cloud_set (cls, &
                                           i-1, &  !!! ipart : 0 -> n_part-1 !!!
                                           n_local_points_tgt, &
                                           coords_tgt, &
                                           gnum_tgt)
  end do


  !
  ! Compute the 'n' closest neighbors into the source point cloud for any taget point
  !

  call PDM_closest_points_compute (cls)

  !
  ! Dump the time used to compute
  !

  call PDM_closest_points_dump_times (cls)

  !
  ! Get the 'n' closest neighbors into the source point cloud for any taget point
  !


  do i = 1, n_part_cloud_tgt
    call PDM_closest_points_get (cls, &
                                 i-1, & !!! ipart : 0 -> n_part-1 !!!
                                 closest_src_gnum, &
                                 closest_src_distance)
  end do


  do i = 1, n_local_points_tgt
    print *, "Closest points for ", gnum_tgt(i), ":", closest_src_gnum(2*(i-1)+1),&
         "/", closest_src_distance(2*(i-1)+1), " and ", closest_src_gnum(2*(i-1)+2)&
         ,"/", closest_src_distance(2*(i-1)+2)

  end do

  !
  ! Free the current cloest_point structure
  !

  call PDM_closest_points_free (cls)

  deallocate (coords_src)
  deallocate (coords_tgt)
  deallocate (gnum_src)
  deallocate (gnum_tgt)


  call mpi_finalize(code)


end program testf
