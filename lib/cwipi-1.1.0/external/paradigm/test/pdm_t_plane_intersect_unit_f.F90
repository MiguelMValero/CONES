!-----------------------------------------------------------------------------
! This file is part of the ParaDiGMA library.
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
  use pdm_overlay
  use iso_c_binding
  use pdm_fortran

  implicit none

#ifndef PDM_HAVE_FORTRAN_MPI_MODULE  
  include "mpif.h"
#endif  

  integer :: code
  integer :: i_rank
  integer :: n_rank

  integer(c_int), parameter :: fComm = MPI_COMM_WORLD
  type (c_ptr) :: c_ol

  integer (c_int), parameter ::  n_partMeshA = 1
  integer (c_int), parameter ::  n_partMeshB = 1
  double precision, parameter :: projectCoeff = 0.d0


  integer (c_int), parameter              :: ipartA = 0
  integer (c_int), parameter              :: nFaceA = 1
  integer (c_int), pointer                :: faceVtxIdxA(:)
  integer (c_int), pointer                :: faceVtxA(:)
  integer (kind = pdm_g_num_s), pointer   :: faceLNToGNA(:)
  integer (c_int), parameter              :: nVtxA = 3
  double precision, pointer               :: vtxCoordA(:)
  integer (kind = pdm_g_num_s), pointer   :: vtxLNToGNA(:)

  integer (kind = pdm_g_num_s) :: nGOlFaceA
  integer (kind = pdm_g_num_s) :: nGOlVtxA
#ifdef PDM_LONG_G_NUM
  integer (c_long)             :: c_nGOlFaceA
  integer (c_long)             :: c_nGOlVtxA
#else
  integer (c_int)              :: c_nGOlFaceA
  integer (c_int)              :: c_nGOlVtxA
#endif

  integer (c_int) :: nOlFaceA
  integer (c_int) :: nOlLinkedFaceA
  integer (c_int) :: nOlVtxA
  integer (c_int) :: sOlFaceIniVtxA
  integer (c_int) :: sOlface_vtxA
  integer (c_int) :: sInitToOlFaceA

  integer (c_int), pointer                :: olFaceIniVtxIdxA(:)
  integer (c_int), pointer                :: olFaceIniVtxA(:)
  integer (c_int), pointer                :: olface_vtx_idxA(:)
  integer (c_int), pointer                :: olface_vtxA(:)
  integer (c_int), pointer                :: olLinkedface_procIdxA(:)
  integer (c_int), pointer                :: olLinkedFaceA(:)
  integer (kind = pdm_g_num_s), pointer   :: olface_ln_to_gnA(:)
  double precision, pointer               :: olCoordsA(:)
  integer (kind = pdm_g_num_s), pointer   :: olvtx_ln_to_gnA(:)
  integer (c_int), pointer                :: initToOlFaceIdxA(:)
  integer (c_int), pointer                :: initToOlFaceA(:)


  integer (c_int), parameter              :: ipartB = 0
  integer (c_int), parameter              :: nFaceB = 1
  integer (c_int), pointer                :: faceVtxIdxB(:)
  integer (c_int), pointer                :: faceVtxB(:)
  integer (kind = pdm_g_num_s), pointer   :: faceLNToGNB(:)
  integer (c_int), parameter              :: nVtxB = 3
  double precision, pointer               :: vtxCoordB(:)
  integer (kind = pdm_g_num_s), pointer   :: vtxLNToGNB(:)

  integer (kind = pdm_g_num_s) :: nGOlFaceB
  integer (kind = pdm_g_num_s) :: nGOlVtxB
#ifdef PDM_LONG_G_NUM
  integer (c_long)             :: c_nGOlFaceB
  integer (c_long)             :: c_nGOlVtxB
#else
  integer (c_int)              :: c_nGOlFaceB
  integer (c_int)              :: c_nGOlVtxB
#endif

  integer (c_int) :: nOlFaceB
  integer (c_int) :: nOlLinkedFaceB
  integer (c_int) :: nOlVtxB
  integer (c_int) :: sOlFaceIniVtxB
  integer (c_int) :: sOlface_vtxB
  integer (c_int) :: sInitToOlFaceB

  integer (c_int), pointer                :: olFaceIniVtxIdxB(:)
  integer (c_int), pointer                :: olFaceIniVtxB(:)
  integer (c_int), pointer                :: olface_vtx_idxB(:)
  integer (c_int), pointer                :: olface_vtxB(:)
  integer (c_int), pointer                :: olLinkedface_procIdxB(:)
  integer (c_int), pointer                :: olLinkedFaceB(:)
  integer (kind = pdm_g_num_s), pointer   :: olface_ln_to_gnB(:)
  double precision, pointer               :: olCoordsB(:)
  integer (kind = pdm_g_num_s), pointer   :: olvtx_ln_to_gnB(:)
  integer (c_int), pointer                :: initToOlFaceIdxB(:)
  integer (c_int), pointer                :: initToOlFaceB(:)

  type (c_ptr)                                :: c_faceVtxIdxA
  type (c_ptr)                                :: c_faceVtxA
  type (c_ptr)                                :: c_faceLNToGNA
  type (c_ptr)                                :: c_vtxCoordA
  type (c_ptr)                                :: c_vtxLNToGNA
  type (c_ptr)                                :: c_faceVtxIdxB
  type (c_ptr)                                :: c_faceVtxB
  type (c_ptr)                                :: c_faceLNToGNB
  type (c_ptr)                                :: c_vtxCoordB
  type (c_ptr)                                :: c_vtxLNToGNB
  type (c_ptr)                                :: c_olFaceIniVtxIdxA
  type (c_ptr)                                :: c_olFaceIniVtxA
  type (c_ptr)                                :: c_olface_vtx_idxA
  type (c_ptr)                                :: c_olface_vtxA
  type (c_ptr)                                :: c_olLinkedface_procIdxA
  type (c_ptr)                                :: c_olLinkedFaceA
  type (c_ptr)                                :: c_olface_ln_to_gnA
  type (c_ptr)                                :: c_olCoordsA
  type (c_ptr)                                :: c_olvtx_ln_to_gnA
  type (c_ptr)                                :: c_initToOlFaceIdxA
  type (c_ptr)                                :: c_initToOlFaceA
  type (c_ptr)                                :: c_olFaceIniVtxIdxB
  type (c_ptr)                                :: c_olFaceIniVtxB
  type (c_ptr)                                :: c_olface_vtx_idxB
  type (c_ptr)                                :: c_olface_vtxB
  type (c_ptr)                                :: c_olLinkedface_procIdxB
  type (c_ptr)                                :: c_olLinkedFaceB
  type (c_ptr)                                :: c_olface_ln_to_gnB
  type (c_ptr)                                :: c_olCoordsB
  type (c_ptr)                                :: c_olvtx_ln_to_gnB
  type (c_ptr)                                :: c_initToOlFaceIdxB
  type (c_ptr)                                :: c_initToOlFaceB

  !
  ! Init
  !

  call mpi_init(code)
  call mpi_comm_rank(mpi_comm_world, i_rank, code)
  call mpi_comm_size(mpi_comm_world, n_rank, code)

  if (n_rank .ne. 1) then
    print *,'Error : 1 MPI process is mandatory'
    call mpi_finalize(code)
    stop
  end if

  !
  ! Set MeshA
  !

  allocate (faceVtxIdxA(nFaceA+1))
  allocate (faceVtxA   (nVtxA   ))
  allocate (faceLNToGNA(nFaceA  ))
  allocate (vtxCoordA  (3*nVtxA ))
  allocate (vtxLNToGNA (nVtxA   ))

  faceVtxIdxA(1) = 0
  faceVtxIdxA(2) = 3

  faceVtxA(1) = 1
  faceVtxA(2) = 2
  faceVtxA(3) = 3

  faceLNToGNA(1) = 1

  vtxCoordA(1) = 0.d0
  vtxCoordA(2) = 0.d0
  vtxCoordA(3) = 0.d0

  vtxCoordA(4) = 0.5d0
  vtxCoordA(5) = sqrt(3.d0)/2.d0
  vtxCoordA(6) = 0.d0

  vtxCoordA(7) = 1.d0
  vtxCoordA(8) = 0.d0
  vtxCoordA(9) = 0.d0

  vtxLNToGNA(1) = 1
  vtxLNToGNA(2) = 2
  vtxLNToGNA(3) = 3


  !
  ! Set MeshB
  !

  allocate (faceVtxIdxB(nFaceB+1))
  allocate (faceVtxB   (nVtxB)   )
  allocate (faceLNToGNB(nFaceB)  )
  allocate (vtxCoordB  (3*nVtxB) )
  allocate (vtxLNToGNB (nVtxB)   )

  faceVtxIdxB(1) = 0
  faceVtxIdxB(2) = 3

  faceVtxB(1) = 1
  faceVtxB(2) = 2
  faceVtxB(3) = 3

  faceLNToGNB(1) = 1

  vtxCoordB(1) = 0.5d0
  vtxCoordB(2) = -sqrt(3.d0)/4.d0
  vtxCoordB(3) = 0.d0

  vtxCoordB(4) = 1.d0
  vtxCoordB(5) = sqrt(3.d0)/4.d0
  vtxCoordB(6) = 0.d0

  vtxCoordB(7) = 0.d0
  vtxCoordB(8) = sqrt(3.d0)/4.d0
  vtxCoordB(9) = 0.d0

  vtxLNToGNB(1) = 1
  vtxLNToGNB(2) = 2
  vtxLNToGNB(3) = 3

  !
  ! Create a new PDM_mesh_location structure
  !   The MPI communicator and the number of point cloud are setted
  !

  call PDM_ol_create (c_ol,         &
                      n_partMeshA,  &
                      n_partMeshB,  &
                      projectCoeff, &
                      fComm)

  call PDM_ol_parameter_set (c_ol,                  &
                             PDM_OL_CAR_LENGTH_TOL, &
                             1.d-4)

  call PDM_ol_parameter_set (c_ol,               &
                             PDM_OL_EXTENTS_TOL, &
                             1.d-4)

  c_faceVtxIdxA = c_loc (faceVtxIdxA)
  c_faceVtxA    = c_loc (faceVtxA)
  c_faceLNToGNA = c_loc (faceLNToGNA)
  c_vtxCoordA   = c_loc (vtxCoordA)
  c_vtxLNToGNA  = c_loc (vtxLNToGNA)

  call PDM_ol_input_mesh_set (c_ol,          &
                              PDM_OL_MESH_A, &
                              ipartA,        &
                              nFaceA,        &
                              c_faceVtxIdxA, &
                              c_faceVtxA,    &
                              c_faceLNToGNA, &
                              nVtxA,         &
                              c_vtxCoordA,   &
                              c_vtxLNToGNA)

  c_faceVtxIdxB = c_loc (faceVtxIdxB)
  c_faceVtxB    = c_loc (faceVtxB)
  c_faceLNToGNB = c_loc (faceLNToGNB)
  c_vtxCoordB   = c_loc (vtxCoordB)
  c_vtxLNToGNB  = c_loc (vtxLNToGNB)

  call PDM_ol_input_mesh_set (c_ol,          &
                              PDM_OL_MESH_B, &
                              ipartB,        &
                              nFaceB,        &
                              c_faceVtxIdxB, &
                              c_faceVtxB,    &
                              c_faceLNToGNB, &
                              nVtxB,         &
                              c_vtxCoordB,   &
                              c_vtxLNToGNB)

  call PDM_ol_compute (c_ol);

  call PDM_ol_dump_times (c_ol);

  call PDM_ol_mesh_dim_get (c_ol,          &
                            PDM_OL_MESH_A, &
                            c_nGOlFaceA,   &
                            c_nGOlVtxA)
  nGOlFaceA = c_nGOlFaceA
  nGOlVtxA  = c_nGOlVtxA

  call PDM_ol_part_mesh_dim_get (c_ol,           &
                                 PDM_OL_MESH_A,  &
                                 ipartA,         &
                                 nOlFaceA,       &
                                 nOlLinkedFaceA, &
                                 nOlVtxA,        &
                                 sOlFaceIniVtxA, &
                                 sOlface_vtxA,   &
                                 sInitToOlFaceA);

  call PDM_ol_mesh_entities_get (c_ol,                    &
                                 PDM_OL_MESH_A,           &
                                 ipartA,                  &
                                 c_olFaceIniVtxIdxA,      &
                                 c_olFaceIniVtxA,         &
                                 c_olface_vtx_idxA,       &
                                 c_olface_vtxA,           &
                                 c_olLinkedface_procIdxA, &
                                 c_olLinkedFaceA,         &
                                 c_olface_ln_to_gnA,      &
                                 c_olCoordsA,             &
                                 c_olvtx_ln_to_gnA,       &
                                 c_initToOlFaceIdxA,      &
                                 c_initToOlFaceA)
  call c_f_pointer (c_olFaceIniVtxIdxA,      olFaceIniVtxIdxA,      [nFaceA+1])
  call c_f_pointer (c_olFaceIniVtxA,         olFaceIniVtxA,         [sOlFaceIniVtxA])
  call c_f_pointer (c_olface_vtx_idxA,       olface_vtx_idxA,       [nOlFaceA+1])
  call c_f_pointer (c_olface_vtxA,           olface_vtxA,           [sOlface_vtxA])
  call c_f_pointer (c_olLinkedface_procIdxA, olLinkedface_procIdxA, [n_rank+1])
  call c_f_pointer (c_olLinkedFaceA,         olLinkedFaceA,         [4*nOlLinkedFaceA])
  call c_f_pointer (c_olface_ln_to_gnA,      olface_ln_to_gnA,      [nOlFaceA])
  call c_f_pointer (c_olCoordsA,             olCoordsA,             [3*nOlVtxA])
  call c_f_pointer (c_olvtx_ln_to_gnA,       olvtx_ln_to_gnA,       [nOlVtxA])
  call c_f_pointer (c_initToOlFaceIdxA,      initToOlFaceIdxA,      [nFaceA+1])
  call c_f_pointer (c_initToOlFaceA,         initToOlFaceA,         [sInitToOlFaceA])


  call PDM_ol_mesh_dim_get (c_ol,          &
                            PDM_OL_MESH_B, &
                            c_nGOlFaceB,   &
                            c_nGOlVtxB)
  nGOlFaceB = c_nGOlFaceB
  nGOlVtxB  = c_nGOlVtxB

  call PDM_ol_part_mesh_dim_get (c_ol,           &
                                 PDM_OL_MESH_B,  &
                                 ipartB,         &
                                 nOlFaceB,       &
                                 nOlLinkedFaceB, &
                                 nOlVtxB,        &
                                 sOlFaceIniVtxB, &
                                 sOlface_vtxB,   &
                                 sInitToOlFaceB);

  call PDM_ol_mesh_entities_get (c_ol,                    &
                                 PDM_OL_MESH_B,           &
                                 ipartB,                  &
                                 c_olFaceIniVtxIdxB,      &
                                 c_olFaceIniVtxB,         &
                                 c_olface_vtx_idxB,       &
                                 c_olface_vtxB,           &
                                 c_olLinkedface_procIdxB, &
                                 c_olLinkedFaceB,         &
                                 c_olface_ln_to_gnB,      &
                                 c_olCoordsB,             &
                                 c_olvtx_ln_to_gnB,       &
                                 c_initToOlFaceIdxB,      &
                                 c_initToOlFaceB)
  call c_f_pointer (c_olFaceIniVtxIdxB,      olFaceIniVtxIdxB,      [nFaceB+1])
  call c_f_pointer (c_olFaceIniVtxB,         olFaceIniVtxB,         [sOlFaceIniVtxB])
  call c_f_pointer (c_olface_vtx_idxB,       olface_vtx_idxB,       [nOlFaceB+1])
  call c_f_pointer (c_olface_vtxB,           olface_vtxB,           [sOlface_vtxB])
  call c_f_pointer (c_olLinkedface_procIdxB, olLinkedface_procIdxB, [n_rank+1])
  call c_f_pointer (c_olLinkedFaceB,         olLinkedFaceB,         [4*nOlLinkedFaceB])
  call c_f_pointer (c_olface_ln_to_gnB,      olface_ln_to_gnB,      [nOlFaceB])
  call c_f_pointer (c_olCoordsB,             olCoordsB,             [3*nOlVtxB])
  call c_f_pointer (c_olvtx_ln_to_gnB,       olvtx_ln_to_gnB,       [nOlVtxB])
  call c_f_pointer (c_initToOlFaceIdxB,      initToOlFaceIdxB,      [nFaceB+1])
  call c_f_pointer (c_initToOlFaceB,         initToOlFaceB,         [sInitToOlFaceB])


  call PDM_ol_del (c_ol)

  call pdm_fortran_free_c (c_olFaceIniVtxIdxA)
  call pdm_fortran_free_c (c_olFaceIniVtxA)
  call pdm_fortran_free_c (c_olface_vtx_idxA)
  call pdm_fortran_free_c (c_olface_vtxA)
  call pdm_fortran_free_c (c_olLinkedface_procIdxA)
  call pdm_fortran_free_c (c_olLinkedFaceA)
  call pdm_fortran_free_c (c_olface_ln_to_gnA)
  call pdm_fortran_free_c (c_olCoordsA)
  call pdm_fortran_free_c (c_olvtx_ln_to_gnA)
  call pdm_fortran_free_c (c_initToOlFaceIdxA)
  call pdm_fortran_free_c (c_initToOlFaceA)

  call pdm_fortran_free_c (c_olFaceIniVtxIdxB)
  call pdm_fortran_free_c (c_olFaceIniVtxB)
  call pdm_fortran_free_c (c_olface_vtx_idxB)
  call pdm_fortran_free_c (c_olface_vtxB)
  call pdm_fortran_free_c (c_olLinkedface_procIdxB)
  call pdm_fortran_free_c (c_olLinkedFaceB)
  call pdm_fortran_free_c (c_olface_ln_to_gnB)
  call pdm_fortran_free_c (c_olCoordsB)
  call pdm_fortran_free_c (c_olvtx_ln_to_gnB)
  call pdm_fortran_free_c (c_initToOlFaceIdxB)
  call pdm_fortran_free_c (c_initToOlFaceB)

  deallocate (faceVtxIdxA)
  deallocate (faceVtxA)
  deallocate (faceLNToGNA)
  deallocate (vtxCoordA)
  deallocate (vtxLNToGNA)

  deallocate (faceVtxIdxB)
  deallocate (faceVtxB)
  deallocate (faceLNToGNB)
  deallocate (vtxCoordB)
  deallocate (vtxLNToGNB)

  call mpi_finalize(code)


end program testf
