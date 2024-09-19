!-----------------------------------------------------------------------------
! This file is part of the CWIPI library. 
!
! Copyright (C) 2011  ONERA
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

#include "cwipi_configf.h"

subroutine printStatus(iiunit, status)
  use cwipi
  implicit none
  integer :: status
  integer :: iiunit

  select case (status)
     case(cwipi_exchange_ok)
        write(iiunit,*) "Exchange ok"
     case(cwipi_exchange_bad_receiving)
        write(iiunit,*) "no or bad receiving"
     case default
        write(iiunit,*) "Unknown receiving status"
        stop
  end select

end subroutine printStatus

!
! ------------------------------------------------------------------------------
! Fonction d'interpolation bidon, juste pour voir si c'est bien pris en compte
! ------------------------------------------------------------------------------
!

subroutine  userInterpolation(entitiesDim, &
                              nLocalVertex, &
                              nLocalElement, &
                              nLocalPolyhedra, &
                              nDistantPoint, &
                              localCoordinates, &
                              localConnectivityIndex, &
                              localConnectivity, &
                              localPolyFaceIndex, &
                              localPolyCellToFaceConnec, &
                              localPolyFaceConnecIdx, &
                              localPolyFaceConnec, &
                              disPtsCoordinates, &
                              disPtsLocation, &
                              disPtsDistance, &
                              disPtsBaryCoordIdx, &
                              disPtsBaryCoord, &
                              stride, &
                              solverType, &
                              localField, &
                              distantField)

  use cwipi

  implicit none

  integer :: entitiesDim
  integer :: nLocalVertex
  integer :: nLocalElement
  integer :: nLocalPolyhedra
  integer :: nDistantPoint
  double precision, dimension(*) :: localCoordinates
  integer, dimension(*) :: localConnectivityIndex
  integer, dimension(*) :: localConnectivity
  integer, dimension(*) :: localPolyFaceIndex
  integer, dimension(*) :: localPolyCellToFaceConnec
  integer, dimension(*) :: localPolyFaceConnecIdx
  integer, dimension(*) :: localPolyFaceConnec
  double precision, dimension(*) :: disPtsCoordinates
  integer, dimension(*) :: disPtsLocation
  real (kind=4), dimension(*) :: disPtsDistance
  integer, dimension(*) :: disPtsBaryCoordIdx
  double precision, dimension(*) :: disPtsBaryCoord
  integer :: stride
  integer :: solverType
  double precision, dimension(*) :: localField
  double precision, dimension(*) :: distantField

  integer :: i

  if (solverType .eq. cwipi_solver_cell_center) then
     do i = 1, nDistantPoint
        distantField(i) = localField(disPtsLocation(i))
     enddo
  else
     print*, 'Error in _userInterpolation : bad solver_type'
     stop
  endif

end subroutine userInterpolation

!
! -----------------------------------------------------------------------------
! Programme de tests
! -----------------------------------------------------------------------------
!

program testf

#ifdef CWP_HAVE_FORTRAN_MPI_MODULE  
  use mpi
#endif
  use cwipi

  implicit none

#ifndef CWP_HAVE_FORTRAN_MPI_MODULE  
  include "mpif.h"
#endif  

  interface
       subroutine  userInterpolation(entitiesDim, &
                                      nLocalVertex, &
                                      nLocalElement, &
                                      nLocalPolyhedra, &
                                      nDistantPoint, &
                                      localCoordinates, &
                                      localConnectivityIndex, &
                                      localConnectivity, &
                                      localPolyFaceIndex, &
                                      localPolyCellToFaceConnec, &
                                      localPolyFaceConnecIdx, &
                                      localPolyFaceConnec, &
                                      disPtsCoordinates, &
                                      disPtsLocation, &
                                      disPtsDistance, &
                                      disPtsBaryCoordIdx, &
                                      disPtsBaryCoord, &
                                      stride, &
                                      solverType, &
                                      localField, &
                                      distantField)
         integer :: entitiesDim
         integer :: nLocalVertex
         integer :: nLocalElement
         integer :: nLocalPolyhedra
         integer :: nDistantPoint
         double precision, dimension(*) :: localCoordinates
         integer, dimension(*) :: localConnectivityIndex
         integer, dimension(*) :: localConnectivity
         integer, dimension(*) :: localPolyFaceIndex
         integer, dimension(*) :: localPolyCellToFaceConnec
         integer, dimension(*) :: localPolyFaceConnecIdx
         integer, dimension(*) :: localPolyFaceConnec
         double precision, dimension(*) :: disPtsCoordinates
         integer, dimension(*) :: disPtsLocation
         real (kind=4), dimension(*) :: disPtsDistance
         integer, dimension(*) :: disPtsBaryCoordIdx
         double precision, dimension(*) :: disPtsBaryCoord
         integer :: stride
         integer :: solverType
         double precision, dimension(*) :: localField
         double precision, dimension(*) :: distantField
       end subroutine userInterpolation
    end interface

!  integer, allocatable, dimension(:) :: location
!  integer, allocatable, dimension(:) :: baryCooIdx
!  double precision, allocatable, dimension(:) :: baryCoo
!  double precision, allocatable, dimension(:) :: tmpDbl
  integer :: nNotLocatedPoints

  integer :: localcom
  integer :: irank, currentRank, localcommsize
  character (len = 4) :: proc
  character (len = 5) :: codeName, codeCoupledName
  integer :: code
  integer :: iiunit

  double precision :: xmin = -100.d0
  double precision :: xmax =  100.d0
  double precision :: ymin = -100.d0
  double precision :: ymax =  100.d0

  integer nvertex, nelts

  integer, parameter :: nvertexm = 4000
  integer, parameter :: neltsm = 4000
  integer, parameter :: lconnecindexm = 12000
  integer, parameter :: nptstolocate = 21

  double precision, allocatable, dimension(:) :: coordstolocate

  double precision, allocatable, dimension(:) :: coords
  integer, allocatable, dimension(:) :: connecindex
  integer, allocatable, dimension(:) :: connec

  double precision, allocatable, dimension(:) :: values
  double precision, allocatable, dimension(:) :: localvalues

  real (kind=4), allocatable, dimension(:) :: distLocPts

  integer status

  integer i

  integer :: stride = 1
  integer :: commWorldSize

  integer :: n_partition, n2, codeId

  integer :: nVertexSeg

  integer :: nLocatedPts

  double precision :: randLevel

  call mpi_init(code)
  call mpi_comm_rank(mpi_comm_world, irank, code)
  call mpi_comm_size(mpi_comm_world, commWorldSize, code)

  write(proc,'(i4.4)') irank
  iiunit = 9
  open(unit=iiunit, file='fortran_surf_P1P1_'//proc//'.txt', &
       form='formatted', status='unknown')

  if (irank == 0) then
     print* 
     print*, 'START: fortran_surf_P1P1_savememory'
  endif

  n_partition = 1
  do while ((2 * n_partition**2) < commWorldSize)
     n_partition = n_partition + 1
  enddo

  n2 = 2 * n_partition**2

  if (n2 /= commWorldSize) then
     if (irank == 0) then
        print *, '      Not executed : only available if the number of processus in the form of 2 * n_partition**2'
     endif
     call mpi_finalize(code)
     stop;
  endif

!
! -----------------------------------------
! Initialisation de l'interface de couplage
! -----------------------------------------
!
 
  call cwipi_set_output_listing_f(iiunit)

  if (irank < commWorldSize / 2) then
     codeName = 'code1'
     codeId = 1
     codeCoupledName = 'code2'
  else 
     codeName = 'code2'
     codeId = 2
     codeCoupledName = 'code1'
  endif

  call cwipi_init_f (mpi_comm_world, &
                     codeName, &
                     localcom)

!
! ------------------------------------------------
! Creation du fichier de sortie listing
! (en parallele, un fichier listing par processus)
! ------------------------------------------------
!

  call mpi_comm_rank(localcom, currentRank, code)
  call mpi_comm_size(localcom, localcommsize, code)

  write(iiunit,*)
  write(iiunit,*) "dump apres initialisation"
  write(iiunit,*) "-------------------------"
  write(iiunit,*)

  call cwipi_dump_appli_properties_f

!
! -------------------------------------
! Test de definition des points
! a interpoler
! -------------------------------------
!
  write(iiunit, *)
  write(iiunit, *) "--------------------------------------------------------"
  write(iiunit, *)
  write(iiunit, *) " Test 3"
  write(iiunit, *)

  if (irank == 0) then
     print*, '       Create coupling'  
  endif

  call cwipi_create_coupling_f("test2D_3", &
                               cwipi_cpl_parallel_with_part,&
                               codeCoupledName, &
                               2,     & ! Dimension des entites geometriques
                               0.1d0, & ! Tolerance geometrique
                               cwipi_static_mesh, &
                               cwipi_solver_cell_center, &
                               1, &
                               "Ensight Gold",&
                               "text")

!
! Construction du maillage

  nVertexSeg   = 10
  randLevel    = 0.1d0
  nvertex      = nVertexSeg * nVertexSeg
  nelts        = (nVertexSeg - 1) * (nVertexSeg - 1)

  allocate(coords(3 * nvertex))
  allocate(coordstolocate(3 * nvertex))

  allocate(connecindex(nelts + 1))
  allocate(connec(4 * nelts))

  allocate(values(nelts))
  allocate(localvalues(nvertex))

  if (irank == 0) then
     print*, '       Create mesh'
  endif

  call grid_mesh_f(xmin, &
                   xmax, &
                   ymin, &
                   ymax, &
                   randLevel, &
                   nVertexSeg, &
                   n_partition, & 
                   coords,  &
                   connecindex,&
                   connec,&
                   localcom)

  call cwipi_define_mesh_f("test2D_3", &
                           nvertex, &
                           nelts, &
                           coords, &
                           connecindex, &
                           connec)

!
! Definition des points a localiser

  do i = 1, 3 * nVertex
     coordstolocate(i) = 0.75 * coords(i)
  enddo

  if (irank == 0) then
     print*, '       Set points to locate'
  endif

  call cwipi_set_points_to_locate_f("test2D_3", &
                                     nptstolocate, &
                                     coordstolocate)


!
! Envoi de la coordonnee Y a codeC
! Reception de la coordonnee Y provenant de codec

  do i = 1, nelts
     if (codeId == 1) then
        values(i) = coords(3*(i-1) + 1)
     else
        values(i) = coords(3*(i-1) + 2)
     endif
  enddo

  stride = 1

  if (irank == 0) then
     print*, '       Exchange'
  endif

  call cwipi_exchange_f ("test2D_3", &
                         "echange1", &
                         stride, &
                         1, &
                         0.1d0, &
                         "cooy", &
                         values, &
                         "coox", &
                         localvalues, &
                         userInterpolation, &
                         nNotLocatedPoints, &
                         status)
  
  call cwipi_get_n_located_pts_f("test2D_3", nLocatedPts)

  allocate(distLocPts(nLocatedPts))

  call cwipi_dist_located_pts_get_f("test2D_3", distLocPts)

  call printStatus(iiunit, status)
  write(iiunit,*) "valeurs recues test2D_3"
  write(iiunit,*) (localvalues(i),i=1,nptstolocate)

!
! Suppression de l'objet couplage

  if (irank == 0) then
     print*, '       Delete coupling'
  endif

  call cwipi_delete_coupling_f("test2D_3");
  write(iiunit, *)
  write(iiunit, *) "--------------------------------------------------------"
  write(iiunit, *)

  call cwipi_finalize_f()

  deallocate(distLocPts)
  deallocate(coords)
  deallocate(coordstolocate)
  deallocate(connecindex)
  deallocate(connec)
  deallocate(values)
  deallocate(localvalues)

  call mpi_finalize(code)

end program testf
