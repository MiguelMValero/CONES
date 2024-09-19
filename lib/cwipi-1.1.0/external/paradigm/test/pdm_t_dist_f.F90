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
  use PDM_dist_cloud_surf
  use iso_c_binding
  use pdm_fortran

  implicit none

#ifndef PDM_HAVE_FORTRAN_MPI_MODULE
  include "mpif.h"
#endif

  !-----------------------------------------------------------
  integer,                   parameter :: f_comm  = MPI_COMM_WORLD
  integer,                   parameter :: n_part  = 1
  integer,                   parameter :: n_pts   = 4
  integer,                   parameter :: n_vtx   = 17
  integer,                   parameter :: n_face  = 12

  type(c_ptr)                          :: dist = C_NULL_PTR
  double precision,          pointer   :: pts_coord(:,:)   => null()
  double precision,          pointer   :: vtx_coord(:,:)   => null()
  integer(kind=pdm_g_num_s), pointer   :: pts_ln_to_gn(:)  => null()
  integer(kind=pdm_g_num_s), pointer   :: vtx_ln_to_gn(:)  => null()
  integer(kind=pdm_g_num_s), pointer   :: face_ln_to_gn(:) => null()
  integer,                   pointer   :: face_vtx_idx(:)  => null()
  integer,                   pointer   :: face_vtx(:)      => null()

  double precision,          pointer   :: closest_elt_distance(:)    => null()
  double precision,          pointer   :: closest_elt_projected(:,:) => null()
  integer(kind=pdm_g_num_s), pointer   :: closest_elt_gnum(:)        => null()

  integer                              :: code
  integer                              :: i_rank
  integer                              :: n_rank
  integer                              :: n_error
  integer                              :: i, j, k
  integer(kind=pdm_g_num_s)            :: offset_gnum
  double precision                     :: offset_y
  integer                              :: fid = 13
  character                            :: strnum
  !-----------------------------------------------------------


  !  Init
  call mpi_init(code)
  call mpi_comm_rank(f_comm, i_rank, code)
  call mpi_comm_size(f_comm, n_rank, code)

  if (n_rank .ne. 2) then
  print *,'Error : 2 MPI processes are mandatory'
  call mpi_finalize(code)
  stop
  end if


  !  Create PDM_dist_cloud_surf object
  call pdm_dist_cloud_surf_create(dist,   &
                                  1,      & ! mesh_nature
                                  1,      & ! n_point_cloud
                                  f_comm, &
                                  1)        ! ownership


  !  Set point cloud
  if (i_rank .eq. 0) then
    write(*, *) "-- Set point cloud"
  end if
  allocate(pts_ln_to_gn(n_pts))
  allocate(pts_coord(3,n_pts))
  if (i_rank == 0) then
    offset_gnum = 0
    offset_y    = 0.d0
  else
    offset_gnum = 4
    offset_y    = 1.d0
  endif

  k = 1
  do j = 1,2
    do i = 1,2
      pts_ln_to_gn(k) = offset_gnum + k
      pts_coord(1,k)  = i-1 + 0.5
      pts_coord(2,k)  = offset_y + 0.5
      pts_coord(3,k)  = j-1 + 0.5
      k = k+1
    enddo
  enddo
  ! print *, i_rank, pts_ln_to_gn
  ! print *, i_rank, pts_coord


  call pdm_dist_cloud_surf_n_part_cloud_set(dist,   &
                                            0,      & ! i_point_cloud
                                            n_part)

  call pdm_dist_cloud_surf_cloud_set(dist,         &
                                     0,            & ! i_point_cloud
                                     0,            & ! i_part
                                     n_pts,        &
                                     pts_coord,    &
                                     pts_ln_to_gn)


  !  Set surface mesh
  if (i_rank .eq. 0) then
    write(*, *) "-- Set surface mesh"
  end if
  allocate(face_ln_to_gn(n_face))
  allocate(face_vtx_idx(n_face+1))
  allocate(face_vtx(4*n_face))
  if (i_rank == 0) then
    offset_gnum = 0
  else
    offset_gnum = 12
  endif
  face_vtx_idx(1) = 0
  do i = 1,n_face
    face_ln_to_gn(i) = offset_gnum + i
    face_vtx_idx(i+1) = face_vtx_idx(i) + 4
  enddo

  allocate(vtx_ln_to_gn(n_vtx))
  allocate(vtx_coord(3,n_vtx))
  if (i_rank == 0) then

    do i = 1, n_vtx
      vtx_ln_to_gn(i) = i
    enddo

    vtx_coord(1, 1) = 0.d0; vtx_coord(2, 1) = 0.d0; vtx_coord(3, 1) = 0.d0
    vtx_coord(1, 2) = 1.d0; vtx_coord(2, 2) = 0.d0; vtx_coord(3, 2) = 0.d0
    vtx_coord(1, 3) = 2.d0; vtx_coord(2, 3) = 0.d0; vtx_coord(3, 3) = 0.d0
    vtx_coord(1, 4) = 0.d0; vtx_coord(2, 4) = 1.d0; vtx_coord(3, 4) = 0.d0
    vtx_coord(1, 5) = 1.d0; vtx_coord(2, 5) = 1.d0; vtx_coord(3, 5) = 0.d0
    vtx_coord(1, 6) = 2.d0; vtx_coord(2, 6) = 1.d0; vtx_coord(3, 6) = 0.d0
    vtx_coord(1, 7) = 0.d0; vtx_coord(2, 7) = 2.d0; vtx_coord(3, 7) = 0.d0
    vtx_coord(1, 8) = 1.d0; vtx_coord(2, 8) = 2.d0; vtx_coord(3, 8) = 0.d0
    vtx_coord(1, 9) = 2.d0; vtx_coord(2, 9) = 2.d0; vtx_coord(3, 9) = 0.d0
    vtx_coord(1,10) = 0.d0; vtx_coord(2,10) = 0.d0; vtx_coord(3,10) = 1.d0
    vtx_coord(1,11) = 1.d0; vtx_coord(2,11) = 0.d0; vtx_coord(3,11) = 1.d0
    vtx_coord(1,12) = 2.d0; vtx_coord(2,12) = 0.d0; vtx_coord(3,12) = 1.d0
    vtx_coord(1,13) = 0.d0; vtx_coord(2,13) = 1.d0; vtx_coord(3,13) = 1.d0
    vtx_coord(1,14) = 2.d0; vtx_coord(2,14) = 1.d0; vtx_coord(3,14) = 1.d0
    vtx_coord(1,15) = 0.d0; vtx_coord(2,15) = 2.d0; vtx_coord(3,15) = 1.d0
    vtx_coord(1,16) = 1.d0; vtx_coord(2,16) = 2.d0; vtx_coord(3,16) = 1.d0
    vtx_coord(1,17) = 2.d0; vtx_coord(2,17) = 2.d0; vtx_coord(3,17) = 1.d0

    face_vtx( 1) =  1; face_vtx( 2) =  2; face_vtx( 3) = 11; face_vtx( 4) = 10;
    face_vtx( 5) =  2; face_vtx( 6) =  3; face_vtx( 7) = 12; face_vtx( 8) = 11;
    face_vtx( 9) =  7; face_vtx(10) = 15; face_vtx(11) = 16; face_vtx(12) =  8;
    face_vtx(13) =  8; face_vtx(14) = 16; face_vtx(15) = 17; face_vtx(16) =  9;
    face_vtx(17) =  3; face_vtx(18) =  6; face_vtx(19) = 14; face_vtx(20) = 12;
    face_vtx(21) =  6; face_vtx(22) =  9; face_vtx(23) = 17; face_vtx(24) = 14;
    face_vtx(25) =  1; face_vtx(26) = 10; face_vtx(27) = 13; face_vtx(28) =  4;
    face_vtx(29) =  4; face_vtx(30) = 13; face_vtx(31) = 15; face_vtx(32) =  7;
    face_vtx(33) =  1; face_vtx(34) =  4; face_vtx(35) =  5; face_vtx(36) =  2;
    face_vtx(37) =  2; face_vtx(38) =  5; face_vtx(39) =  6; face_vtx(40) =  3;
    face_vtx(41) =  4; face_vtx(42) =  7; face_vtx(43) =  8; face_vtx(44) =  5;
    face_vtx(45) =  5; face_vtx(46) =  8; face_vtx(47) =  9; face_vtx(48) =  6;


  else

    do i = 1, n_vtx
      vtx_ln_to_gn(i) = i + 9
    enddo

    vtx_coord(1, 1) = 0.d0; vtx_coord(2, 1) = 0.d0; vtx_coord(3, 1) = 1.d0
    vtx_coord(1, 2) = 1.d0; vtx_coord(2, 2) = 0.d0; vtx_coord(3, 2) = 1.d0
    vtx_coord(1, 3) = 2.d0; vtx_coord(2, 3) = 0.d0; vtx_coord(3, 3) = 1.d0
    vtx_coord(1, 4) = 0.d0; vtx_coord(2, 4) = 1.d0; vtx_coord(3, 4) = 1.d0
    vtx_coord(1, 5) = 2.d0; vtx_coord(2, 5) = 1.d0; vtx_coord(3, 5) = 1.d0
    vtx_coord(1, 6) = 0.d0; vtx_coord(2, 6) = 2.d0; vtx_coord(3, 6) = 1.d0
    vtx_coord(1, 7) = 1.d0; vtx_coord(2, 7) = 2.d0; vtx_coord(3, 7) = 1.d0
    vtx_coord(1, 8) = 2.d0; vtx_coord(2, 8) = 2.d0; vtx_coord(3, 8) = 1.d0
    vtx_coord(1, 9) = 0.d0; vtx_coord(2, 9) = 0.d0; vtx_coord(3, 9) = 2.d0
    vtx_coord(1,10) = 1.d0; vtx_coord(2,10) = 0.d0; vtx_coord(3,10) = 2.d0
    vtx_coord(1,11) = 2.d0; vtx_coord(2,11) = 0.d0; vtx_coord(3,11) = 2.d0
    vtx_coord(1,12) = 0.d0; vtx_coord(2,12) = 1.d0; vtx_coord(3,12) = 2.d0
    vtx_coord(1,13) = 1.d0; vtx_coord(2,13) = 1.d0; vtx_coord(3,13) = 2.d0
    vtx_coord(1,14) = 2.d0; vtx_coord(2,14) = 1.d0; vtx_coord(3,14) = 2.d0
    vtx_coord(1,15) = 0.d0; vtx_coord(2,15) = 2.d0; vtx_coord(3,15) = 2.d0
    vtx_coord(1,16) = 1.d0; vtx_coord(2,16) = 2.d0; vtx_coord(3,16) = 2.d0
    vtx_coord(1,17) = 2.d0; vtx_coord(2,17) = 2.d0; vtx_coord(3,17) = 2.d0

    face_vtx( 1) =  1; face_vtx( 2) =  2; face_vtx( 3) = 10; face_vtx( 4) =  9;
    face_vtx( 5) =  2; face_vtx( 6) =  3; face_vtx( 7) = 11; face_vtx( 8) = 10;
    face_vtx( 9) =  6; face_vtx(10) = 15; face_vtx(11) = 16; face_vtx(12) =  7;
    face_vtx(13) =  7; face_vtx(14) = 16; face_vtx(15) = 17; face_vtx(16) =  8;
    face_vtx(17) =  3; face_vtx(18) =  5; face_vtx(19) = 14; face_vtx(20) = 11;
    face_vtx(21) =  5; face_vtx(22) =  8; face_vtx(23) = 17; face_vtx(24) = 14;
    face_vtx(25) =  1; face_vtx(26) =  9; face_vtx(27) = 12; face_vtx(28) =  4;
    face_vtx(29) =  4; face_vtx(30) = 12; face_vtx(31) = 15; face_vtx(32) =  6;
    face_vtx(33) =  9; face_vtx(34) = 10; face_vtx(35) = 13; face_vtx(36) = 12;
    face_vtx(37) = 10; face_vtx(38) = 11; face_vtx(39) = 14; face_vtx(40) = 13;
    face_vtx(41) = 12; face_vtx(42) = 13; face_vtx(43) = 16; face_vtx(44) = 15;
    face_vtx(45) = 13; face_vtx(46) = 14; face_vtx(47) = 17; face_vtx(48) = 16;

  endif

  ! do i = 1, n_vtx
  !   print *, i_rank, i, vtx_ln_to_gn(i), vtx_coord(1,i), vtx_coord(2,i)), vtx_coord(3,i)
  ! enddo


  !   VTK export
  if (.false.) then

    write (strnum, '(i1)') i_rank

    open(unit=fid, file="points_"//strnum//".vtk", action='write')
    write(fid,'(a)') "# vtk DataFile Version 2.0"
    write(fid,'(a)') "point cloud"
    write(fid,'(a)') "ASCII"
    write(fid,'(a)') "DATASET UNSTRUCTURED_GRID"
    write(fid,'(a7,i0,a7)') "POINTS ", n_pts, " double"
    do i = 1, n_pts
      write(fid, *) pts_coord(1,i), pts_coord(2,i), pts_coord(3,i)
    enddo
    write(fid,'(a5,1x,i0,1x,i0)') "CELLS", n_pts, 2*n_pts
    do i = 1, n_pts
      write(fid, *) 1, i-1
    enddo
    write(fid,'(a10,1x,i0)') "CELL_TYPES", n_pts
    do i = 1, n_pts
      write(fid, '(i0)') 1
    enddo
    write(fid,'(a9,1x,i0)') "CELL_DATA", n_pts
    write(fid,'(a)') "SCALARS pts_ln_to_gn long 1"
    write(fid,'(a)') "LOOKUP_TABLE default"
    do i = 1, n_pts
      write(fid, '(i0)') pts_ln_to_gn(i)
    enddo
    close(fid)


    open(unit=fid, file="surfmesh_"//strnum//".vtk", action='write')
    write(fid,'(a)') "# vtk DataFile Version 2.0"
    write(fid,'(a)') "surface mesh"
    write(fid,'(a)') "ASCII"
    write(fid,'(a)') "DATASET POLYDATA"
    write(fid,'(a7,i0,a7)') "POINTS ", n_vtx, " double"
    do i = 1, n_vtx
      write(fid, *) vtx_coord(1,i), vtx_coord(2,i), vtx_coord(3,i)
    enddo
    write(fid,'(a8,1x,i0,1x,i0)') "POLYGONS", n_face, n_face+face_vtx_idx(n_face+1)
    do i = 1, n_face
      write(fid, '(i0,1x)', advance='no') face_vtx_idx(i+1) - face_vtx_idx(i)
      do j = face_vtx_idx(i), face_vtx_idx(i+1)-1
        write(fid, '(i0,1x)', advance='no') face_vtx(j+1)-1
      enddo
      write(fid, *) ""
    enddo
    write(fid,'(a9,1x,i0)') "CELL_DATA", n_face
    write(fid,'(a)') "SCALARS face_ln_to_gn long 1"
    write(fid,'(a)') "LOOKUP_TABLE default"
    do i = 1, n_face
      write(fid, '(i0)') face_ln_to_gn(i)
    enddo
    write(fid,'(a10,1x,i0)') "POINT_DATA", n_vtx
    write(fid,'(a)') "SCALARS vtx_ln_to_gn long 1"
    write(fid,'(a)') "LOOKUP_TABLE default"
    do i = 1, n_vtx
      write(fid, '(i0)') vtx_ln_to_gn(i)
    enddo
    close(fid)

  endif



  call pdm_dist_cloud_surf_surf_mesh_global_data_set(dist,    &
                                                     n_part)

  call pdm_dist_cloud_surf_surf_mesh_part_set(dist,          &
                                              0,             & ! i_part
                                              n_face,        &
                                              face_vtx_idx,  &
                                              face_vtx,      &
                                              face_ln_to_gn, &
                                              n_vtx,         &
                                              vtx_coord,     &
                                              vtx_ln_to_gn)


  !  Compute distance
  if (i_rank .eq. 0) then
    write(*, *) "-- Compute distance"
  end if
  call pdm_dist_cloud_surf_compute(dist)

  call pdm_dist_cloud_surf_dump_times(dist)


  !  Check results
  if (i_rank .eq. 0) then
    write(*, *) "-- Check"
  end if
  n_error = 0
  call pdm_dist_cloud_surf_get(dist,                  &
                               0,                     & ! i_point_cloud
                               0,                     & ! i_part
                               closest_elt_distance,  &
                               closest_elt_projected, &
                               closest_elt_gnum)

  do i = 1, n_pts
    if (abs(sqrt(closest_elt_distance(i)) - 0.5) > 1.d-6) then
      n_error = n_error + 1
      write (*, *) pts_ln_to_gn(i), closest_elt_gnum(i), sqrt(closest_elt_distance(i)), closest_elt_projected(:,i)
    endif
  enddo


  !  Free memory
  deallocate(pts_coord)
  deallocate(pts_ln_to_gn)
  deallocate(vtx_coord)
  deallocate(vtx_ln_to_gn)
  deallocate(face_ln_to_gn)
  deallocate(face_vtx_idx)
  deallocate(face_vtx)

  call PDM_dist_cloud_surf_free(dist, 0)

  write (*, '(a1,i0,a12,i0)') "[", i_rank, "] n_error = ", n_error

  if (i_rank .eq. 0) then
    write(*, *) "-- End"
  end if

  call mpi_finalize(code)

end program testf
