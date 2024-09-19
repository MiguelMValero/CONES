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


!  mpirun -n 1 ./fortran_surf_TetraP1_PiPj : -n 1 ./fortran_surf_TetraP1_PiPj
!  mpirun -n 1 tests/fortran_surf_TetraP1_PiPj : -n 1 tests/fortran_surf_TetraP1_PiPj

module variablesCommunes
  logical :: visu=.false.
  integer :: commWorld,rankWorld,sizeWorld
  integer :: commLocal,rankLocal,sizeLocal
  integer :: order
end module variablesCommunes


subroutine  userInterpolation                        ( &
  &           entitiesDim                             ,&
  &           nLocalVertex                            ,&
  &           nLocalElement                           ,&
  &           nLocalPolhyedra                         ,&
  &           nDistantPoint                           ,&
  &                                                    &
  &           localCoordinates                        ,&
  &                                                    &
  &           localConnectivityIndex                  ,&
  &           localConnectivity                       ,&
  &           localPolyFaceIndex                      ,&
  &           localPolyCellToFaceConnec               ,&
  &           localPolyFaceConnecIdx                  ,&
  &           localPolyFaceConnec                     ,&
  &                                                    &
  &           disPtsCoordinates                       ,&
  &           disPtsLocation                          ,&
  &           disPtsDistance                          ,&
  &           disPtsBaryCoordIdx                      ,&
  &           distantPointsBarycentricCoordinates     ,&
  &                                                    &
  &           stride                                  ,&  ! =ker(calc)
  &           solverType                              ,&
  &           localField                              ,&  !   mySolu
  &           distantField                             )  ! linkSolu
  !---
  use cwipi
  use modDeterminant
  use baseSimplex3D
  
  use variablesCommunes
  !---
  implicit none
  !---
  integer :: entitiesDim
  integer :: nLocalVertex
  integer :: nLocalElement
  integer :: nLocalPolhyedra
  integer :: nDistantPoint
  real(8) :: localCoordinates                        (*)
  integer :: localConnectivityIndex                  (*)
  integer :: localConnectivity                       (*)
  integer :: localPolyFaceIndex                      (*)
  integer :: localPolyCellToFaceConnec               (*)
  integer :: localPolyFaceConnecIdx                  (*)
  integer :: localPolyFaceConnec                     (*)
  real(8) :: disPtsCoordinates                       (*)
  integer :: disPtsLocation                          (*)
  real(4) :: disPtsDistance                          (*)
  integer :: disPtsBaryCoordIdx                      (*)
  real(8) :: distantPointsBarycentricCoordinates     (*)
  integer :: stride
  integer :: solverType
  real(8) ::   localField                            (*)
  real(8) :: distantField                            (*)
  !>
  integer          :: i,j,k,iErr
  integer          :: iNod,nNod,iMod,nMod
  real(8), pointer :: uvw (:,:),rst(:,:),a(:),b(:),c(:),vand(:,:)
  integer          :: iDistantPoint
  integer          :: iBary,iVert
  real(8), pointer :: uvwOut(:,:),lagrange(:,:)
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  if( rankWorld==0 )print'(/">>> userInterpolation rankWorld=",i2)',rankWorld
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  if( visu .and. rankWorld==0 )then
    print '()'
    iVert=0
    do i=1,nLocalVertex
      print '("localCoordinates(",i3,")=",3(f12.5,1x))',i,localCoordinates(iVert+1:iVert+3)
      iVert=iVert+3
    enddo
  endif
  
  if( visu .and. rankWorld==0 )then
    nMod=(order+1)*(order+2)*(order+3)/6
    print '()'
    j=0
    do iMod=1,nMod
      print '("iMod=",i3," localField       =",4(f12.5,1x),t100,"@rkw",i3)',iMod,localField(j+1:j+stride),rankWorld
      j=j+stride
    enddo
  endif
  
  if( visu .and. rankWorld==0 )then
    print '()'
    iVert=0
    do iDistantPoint=1,nDistantPoint
      print '("iDis=",i3," disPtsCoordinates=",3(f12.5,1x),t100,"@rkw",i3)',&
      & iDistantPoint,disPtsCoordinates(iVert+1:iVert+3),rankWorld
      iVert=iVert+3
    enddo
  endif
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  !> Calcul des coordonnées barycentriques dans la cellule mitoyenne
  if( visu  .and. rankWorld==0 )then
    print '()'
    block
      real(8) :: det1,mat1(3,3)
      real(8) :: det2,mat2(3,3)
      real(8) :: bi(1:3)
      real(8) :: u,v,w
      
      mat2(:,1)=localCoordinates( 4: 6)-localCoordinates(1:3)
      mat2(:,2)=localCoordinates( 7: 9)-localCoordinates(1:3)
      mat2(:,3)=localCoordinates(10:12)-localCoordinates(1:3)
      det2 = mat2(1,1)*(mat2(2,2)*mat2(3,3)-mat2(2,3)*mat2(3,2)) &
      &     -mat2(1,2)*(mat2(2,1)*mat2(3,3)-mat2(3,1)*mat2(2,3)) &
      &     +mat2(1,3)*(mat2(2,1)*mat2(3,2)-mat2(3,1)*mat2(2,2))
      det2=1d0/det2
      
      iVert=0
      do iDistantPoint=1,nDistantPoint
        
        bi(1:3)=disPtsCoordinates(iVert+1:iVert+3)-localCoordinates(1:3)
        
        mat1=mat2 ; mat1(:,1)=bi(1:3)
        det1 = mat1(1,1)*(mat1(2,2)*mat1(3,3)-mat1(2,3)*mat1(3,2)) &
        &     -mat1(1,2)*(mat1(2,1)*mat1(3,3)-mat1(3,1)*mat1(2,3)) &
        &     +mat1(1,3)*(mat1(2,1)*mat1(3,2)-mat1(3,1)*mat1(2,2))
        u=det1*det2
        
        mat1=mat2 ; mat1(:,2)=bi(1:3)
        det1 = mat1(1,1)*(mat1(2,2)*mat1(3,3)-mat1(2,3)*mat1(3,2)) &
        &     -mat1(1,2)*(mat1(2,1)*mat1(3,3)-mat1(3,1)*mat1(2,3)) &
        &     +mat1(1,3)*(mat1(2,1)*mat1(3,2)-mat1(3,1)*mat1(2,2))
        v=det1*det2
        
        mat1=mat2 ; mat1(:,3)=bi(1:3)
        det1 = mat1(1,1)*(mat1(2,2)*mat1(3,3)-mat1(2,3)*mat1(3,2)) &
        &     -mat1(1,2)*(mat1(2,1)*mat1(3,3)-mat1(3,1)*mat1(2,3)) &
        &     +mat1(1,3)*(mat1(2,1)*mat1(3,2)-mat1(3,1)*mat1(2,2))
        w=det1*det2
        
        print '("iDis=",i3," uvw0             =",3(f12.5,1x),t100,"@rkw",i3)',&
        & iDistantPoint,u,v,w,rankWorld
        
        iVert=iVert+3
      enddo
    end block
    
  endif
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  !> Affectation des coordonées barycentriques uvwOut
  allocate(uvwOut(1:3,1:nDistantPoint))
  iBary=0
  do iDistantPoint=1,nDistantPoint
    uvwOut(1:3,iDistantPoint)=distantPointsBarycentricCoordinates(iBary+2:iBary+4) ! <= Attention 2:4
    iBary=iBary+4
  enddo
  
  if( visu .and. rankWorld==0 )then
    print '()'
    do iDistantPoint=1,nDistantPoint
      print '("iDis=",i3," uvwOut           =",3(f12.5,1x),t100,"@rkw",i3)',&
      & iDistantPoint,uvwOut(1:3,iDistantPoint),rankWorld
    enddo
  endif
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  !> Base fonctionelle d'ordre order
  
  nMod=(order+1)*(order+2)*(order+3)/6
  nNod=size(uvwOut,2)
  
  !> transpose = .true. => lagrange(1:nMod,1:nNod)
  call lagrange3Dv(ord=order,uvwOut=uvwOut,lagrange=lagrange,transpose=.true.)
  
!  !> Points d'interpolation
!  call nodes3D   (ord=order,uvw=uvw,display=.false.)
!  call nodes3Dopt(ord=order,uvw=uvw,display=.false.)
!  !> Calcul de Vand(:,:)
!  call nodes3Duvw2abc(uvw=uvw,a=a,b=b,c=c,display=.false.)
!  call vandermonde3D(ord=order,a=a,b=b,c=c,vand=vand)
!  deallocate(uvw,a,b,c)
!  !> Calcul des polonômes de Lagrange d'ordre order en uvwOut
!  allocate(lagrange(1:nMod,1:nNod))
!  call nodes3Duvw2abc(uvw=uvwOut,a=a,b=b,c=c,display=.false.)
!  call lagrange3Dv(ord=order,vand=vand,a=a,b=b,c=c,lx=lagrange,transpose=.true.)  !> lagrange= Inverse[Transpose[Vand]].Psi[xyzOut] lxOut(nPt,np)
  
  if( visu .and. rankWorld==0 )then
    print '()'
    do iNod=1,nNod
      print '("iNod=",i3," lagrange         =",*(f12.5,1x))',iNod,lagrange(1:nMod,iNod)
    enddo
  endif
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  !> Calcul de distantField
  j=0
  do iDistantPoint=1,nDistantPoint
    distantField(j+1:j+stride)=0d0
    k=0
    do iMod=1,nMod
      distantField(j+1:j+stride)= distantField(j+1:j+stride)                          &
      &                          +lagrange(iMod,iDistantPoint)*localField(k+1:k+stride)
      k=k+stride
    enddo
    j=j+stride
  enddo
  
  !> Visu de distantField
  if( visu .and. rankWorld==0 )then
    print '()'
    j=0
    do iDistantPoint=1,nDistantPoint
      print '("iDis=",i3," distantField     =",4(f12.5,1x),t100,"@rkw",i3)',&
      & iDistantPoint,distantField(j+1:j+stride),rankWorld
      j=j+stride
    enddo
  endif
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  deallocate(lagrange)
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  if( rankWorld==0 )print'("<<< userInterpolation rankWorld=",i2)',rankWorld
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  return
end subroutine userInterpolation


program testf
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  use iso_fortran_env
  
#ifdef PDM_HAVE_FORTRAN_MPI_MODULE  
  use mpi
#endif

  use cwipi
  
  use modDeterminant
  use baseSimplex2D
  use baseSimplex3D
  use table_tet_mesh
  
  use variablesCommunes
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  implicit none

#ifndef PDM_HAVE_FORTRAN_MPI_MODULE  
  include "mpif.h"
#endif  
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  interface
    subroutine  userInterpolation                      ( &
    &           entitiesDim                             ,&
    &           nLocalVertex                            ,&
    &           nLocalElement                           ,&
    &           nLocalPolhyedra                         ,&
    &           nDistantPoint                           ,&
    &                                                    &
    &           localCoordinates                        ,&
    &                                                    &
    &           localConnectivityIndex                  ,&
    &           localConnectivity                       ,&
    &           localPolyFaceIndex                      ,&
    &           localPolyCellToFaceConnec               ,&
    &           localPolyFaceConnecIdx                  ,&
    &           localPolyFaceConnec                     ,&
    &                                                    &
    &           disPtsCoordinates                       ,&
    &           disPtsLocation                          ,&
    &           disPtsDistance                          ,&
    &           disPtsBaryCoordIdx                      ,&
    &           distantPointsBarycentricCoordinates     ,&
    &                                                    &
    &           stride                                  ,&  ! =ker(calc)
    &           solverType                              ,&
    &           localField                              ,&  !   mySolu
    &           distantField                             )  ! linkSolu
    !---
    integer :: entitiesDim
    integer :: nLocalVertex
    integer :: nLocalElement
    integer :: nLocalPolhyedra
    integer :: nDistantPoint
    real(8) :: localCoordinates                        (*)
    integer :: localConnectivityIndex                  (*)
    integer :: localConnectivity                       (*)
    integer :: localPolyFaceIndex                      (*)
    integer :: localPolyCellToFaceConnec               (*)
    integer :: localPolyFaceConnecIdx                  (*)
    integer :: localPolyFaceConnec                     (*)
    real(8) :: disPtsCoordinates                       (*)
    integer :: disPtsLocation                          (*)
    real(4) :: disPtsDistance                          (*)
    integer :: disPtsBaryCoordIdx                      (*)
    real(8) :: distantPointsBarycentricCoordinates     (*)
    integer :: stride
    integer :: solverType
    real(8) ::   localField                            (*)
    real(8) :: distantField                            (*)
    end subroutine  userInterpolation
  end interface
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  character(5)     :: codeName,codeCoupledName
  
  character(64)    :: meshName
  character(80)    :: fileName
  
  integer          :: iVert,nVert
  real(8), pointer :: vertx(:,:)
  integer          :: iTetra,nTetra
  integer, pointer :: tetra(:,:)
  integer          :: iTrian,nTrian
  integer, pointer :: trian(:,:)
  
  integer          :: iNod,j,k
  integer          :: iCell,nCell
  real(8), pointer :: vertices   (:)
  integer, pointer :: connec     (:)
  integer, pointer :: connecIndex(:)
  integer, pointer :: tetraNodes(:,:)
  
  real(8), pointer :: lagrange(:,:)
  real(8)          :: lagrangeTrianP1(1:3)
 !real(8)          :: lagrangeTetraP1(1:4)
  real(8)          :: xyz(1:3)
  
  integer          :: linkVertSize
  real(8), pointer :: linkVert(:)
  integer          :: notLocatedPoints
  
  integer          :: stride
  real(8), pointer ::   myValues(:)
  real(8), pointer :: linkValues(:)
  
  integer          :: nNod
  real(8), pointer :: uvw  (:,:),a(:),b(:),c(:)
  real(8), pointer :: uv   (:,:),rs (:,:)
  real(8), pointer :: vand (:,:)
  
  real(8)          :: node_xy(1:2,1:3) !> Triangle
  
  integer          :: iRank,iErr
  
  real(8)          :: delta,deltaMax

  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  call mpi_init(iErr)
  commWorld=mpi_comm_world
  
  call mpi_comm_rank(commWorld, rankWorld, iErr)
  call mpi_comm_size(commWorld, sizeWorld, iErr)
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  if( rankWorld==0) print '(/"START: fortran_surf_TetraP1_PiPj")'
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  ! Initialisation de l'interface de couplage
  
  select case(rankWorld)
  case(0)
     codeName        = "code1"
     codeCoupledName = "code2"
  case(1)
     codeName        = "code2"
     codeCoupledName = "code1"
  end select
  
  call cwipi_init_f(           &
  &    globalComm=commWorld   ,&
  &    appliName=codeName     ,&
  &    appliComm=commLocal     )
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  call mpi_comm_rank(commLocal,rankLocal,iErr)
  call mpi_comm_size(commLocal,sizeLocal,iErr)
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  select case(rankWorld)
  case(0) ; order=07
  case(1) ; order=10
  end select
  print '("fortran_surf_TetraP1_PiPj : Order=",i2,t100,"@rkw",i3)',order,rankWorld
  call mpi_barrier(commWorld,iErr)
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  ! Create coupling
  
  if( rankWorld==0 )print '(/"Create coupling")'
  
  call cwipi_create_coupling_f(                  &
  &    couplingName="testPiPj"                  ,&
  &    couplingType=cwipi_cpl_parallel_with_part,&
  &    cplAppli=codeCoupledName                 ,&
  &    entitiesDim=3                            ,& !> Nature du couplage
  &    tolerance=1d-1                           ,& !> Tolerance geometrique 1d-1 par defaut
  &    meshT=cwipi_static_mesh                  ,&
  &    solvert=cwipi_solver_cell_vertex         ,&
  &    outputfreq=1                             ,& !> Frequence du post-traitement
  &    outputfmt="Ensight Gold"                 ,&
  &    outputfmtopt="binary"                     )
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  ! Create Geometric Mesh
  
  if( rankWorld==0 )print '(/"Create Geometric Mesh (inria mesh format")'
  
  select case(rankWorld)
  case(0)
    !> Vertices
    nVert=4
    allocate(vertx(1:3,1:nVert))
    vertx(1:3,1)=[0.00d0,0.00d0,0.00d0]
    vertx(1:3,2)=[1.00d0,0.00d0,0.00d0]
    vertx(1:3,3)=[0.00d0,1.00d0,0.00d0]
    vertx(1:3,4)=[0.00d0,0.00d0,1.00d0]
    !> Tetrahedra
    nTetra=1
    allocate(tetra(1:5,1:nTetra)) !> 4 sommets + 1 marquer
    tetra(1:5,1)=[1,2,3,4, 1]
    !>  Triangles
    nTrian=4
    allocate(trian(1:4,1:nTrian)) !> 3 sommets + 1 marquer
    trian(1:4,1)=[2,3,4, 1]
    trian(1:4,2)=[1,4,3, 2]
    trian(1:4,3)=[1,2,4, 3]
    trian(1:4,4)=[1,3,2, 4]
  case(1)
    !> Vertices
    nVert=4
    allocate(vertx(1:3,1:nVert))
    vertx(1:3,1)=[0.73d0,0.73d0,0.73d0]
    vertx(1:3,2)=[1.00d0,0.00d0,0.00d0]
    vertx(1:3,3)=[0.00d0,0.00d0,1.00d0]
    vertx(1:3,4)=[0.00d0,1.00d0,0.00d0]
    !> Tetrahedra
    nTetra=1
    allocate(tetra(1:5,1:nTetra)) !> 4 sommets + 1 marquer
    tetra(1:5,1)=[1,2,3,4, 1]
    !>  Triangles
    nTrian=4
    allocate(trian(1:4,1:nTrian)) !> 3 sommets + 1 marquer
    trian(1:4,1)=[2,3,4, 1]
    trian(1:4,2)=[1,4,3, 2]
    trian(1:4,3)=[1,2,4, 3]
    trian(1:4,4)=[1,3,2, 4]
  end select
  
  
  if( visu )then
    !> Ecriture des maillages au format mesh de l'inria
    if( rankWorld==0)then
      do iRank=0,sizeWorld-1
        print '(/"Writing mesh file: Tetra",i1,".mesh")',iRank
      enddo
    endif
    
    write(meshName,'("Tetra",i1,".mesh")')rankWorld
    open(unit=100,file=trim(meshName),action='write',status='unknown')
    write(100,'("MeshVersionFormatted 1"/)')
    write(100,'("Dimension 3"/)')
    write(100,'("Vertices")')
    write(100,'(i1)')nVert
    do iVert=1,nVert
      write(100,'(3(e22.15,1x),i2)')vertx(1:3,iVert),0
    enddo
    write(100,'(/"Tetrahedra")')
    write(100,'(i1)')nTetra
    do iTetra=1,nTetra
      write(100,'(*(i6,1x))')tetra(1:5,iTetra)
    enddo
    write(100,'(/"Triangles")')
    write(100,'(i1)')nTrian
    do iTrian=1,nTrian
      write(100,'(*(i6,1x))')trian(1:4,iTrian)
    enddo
    write(100,'(/"End")')
    close(100)
  endif
  
  !> Mise au format pour cwipi
  
  nVert=4
  nCell=1
  allocate( vertices   (1:3*nVert)    )  !> sommets
  allocate( connec     (1:4*nCell)    )  !> tetra
  allocate( connecIndex(1:nCell+1)    )  !> tetra
  
  connec     (1:4)=tetra(1:4,1)
  connecIndex(1:2)=[0,4]
  
  j=0
  do iNod=1,4
    vertices(j+1:j+3)=vertx(1:3,tetra(iNod,1))
    j=j+3
  enddo
  
  !> Transmission des maillages à cwipi
  call cwipi_define_mesh_f(     &
  &   couplingName="testPiPj"  ,&
  &   nVertex     =nVert       ,&
  &   nElts       =nCell       ,&
  &   coords      =vertices    ,&
  &   connecIndex =connecIndex ,&
  &   connec      =connec       )
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  !> Points de couplage situés Triangle avec le marqueur de peau 1 => face commune aux deux tetras <=
  
  if( rankWorld==0 ) print'(/"Calcul des coordonnees des points de couplage")'
  
  !> Calcul des coordonnees barycentriques
  call nodes3Dopt_2D(ord=order,uvw=uvw,display=.false.)
  
  !> Visu des coordonnees barycentriques dans le triangle unité
  if( visu )then
    if(  0<=order .and. order< 10 )write(fileName,'("pointInterpolationP0",i1,".eps")')order
    if( 10<=order .and. order<100 )write(fileName,'("pointInterpolationP" ,i2,".eps")')order
   !print '("writing File: ",a)',trim(fileName)
    node_xy(1:2,1)=[0,0]
    node_xy(1:2,2)=[1,0]
    node_xy(1:2,3)=[0,1]
    call trianglePointsPlot(   &
    &    file_name=fileName   ,&
    &    node_xy=node_xy      ,&
    &    node_show=0          ,&
    &    point_num=size(uvw,2),&
    &    point_xy=uvw         ,&
    &    point_show=2          ) !> point_show=2, shows the points and number them
  endif
  
  !> Calculs des coordonnées des points de couplage
  linkVertSize=size(uvw,2)
  allocate(linkVert(1:3*linkVertSize))
  j=0
  do iVert=1,linkVertSize
    !> Fonction
    lagrangeTrianP1(1)=1d0-uvw(1,iVert)-uvw(2,iVert)
    lagrangeTrianP1(2)=    uvw(1,iVert)
    lagrangeTrianP1(3)=                 uvw(2,iVert)
    
    linkVert(j+1:j+3)= lagrangeTrianP1(1)*vertx(1:3,trian(1,1)) &
    &                 +lagrangeTrianP1(2)*vertx(1:3,trian(2,1)) &
    &                 +lagrangeTrianP1(3)*vertx(1:3,trian(3,1))
    j=j+3
  enddo
  deallocate(uvw)
  
  
  !> Visu des coordonnees barycentriques dans le triangle unité
  if( visu )then
    do iRank=0,sizeWorld-1
      if( iRank==rankWorld )then
        print '()'
        j=0
        do iVert=1,linkVertSize
          print '("linkVert(",i2,")=",3(f12.5,1x),t100,"@rkw",i3)',iVert,linkVert(j+1:j+3),rankWorld
          j=j+3
        enddo
      endif
      call mpi_barrier(commWorld,iErr)
    enddo
  endif
  
  
  call cwipi_set_points_to_locate_f( &
  &    couplingName="testPiPj"      ,&
  &    nPts  =linkVertSize          ,&
  &    coords=linkVert               )
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  ! Localisation
  
  if( rankWorld==0 )print '(/"Localisation")'
  
  call cwipi_locate_f(couplingName="testPiPj")
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  !> Initialisation of myValues
  
  if( rankWorld==0 )print '(/"Initialisation of myValues")'
  
  call nodes3D   (ord=order,uvw=uvw,display=.false.)
  if( visu )then
    call driverTetMesh(ord=order,node_xyz=uvw,tetra_node=tetraNodes)
  endif
  
  call nodes3Dopt(ord=order,uvw=uvw,display=.false.)
  if( visu )then
    call saveTetMesh  (ord=order,node_xyz=uvw,tetra_node=tetraNodes)
    deallocate(tetraNodes)
  endif
  
  stride=4 ; nNod=size(uvw,2)
  allocate(myValues(1:stride*nNod))
  
  !> lagrange(1:Mod,1:Nod)
  call lagrange3Dv(ord=1,uvwOut=uvw,lagrange=lagrange,transpose=.true.)
  
  j=0
  do iNod=1,nNod
    xyz(1:3)= lagrange(1,iNod)*vertx(1:3,tetra(1,1)) &
    &        +lagrange(2,iNod)*vertx(1:3,tetra(2,1)) &
    &        +lagrange(3,iNod)*vertx(1:3,tetra(3,1)) &
    &        +lagrange(4,iNod)*vertx(1:3,tetra(4,1))
    
    myValues(j+1:j+stride)=[xyz(1),xyz(2),xyz(3),real(rankWorld,kind=8)]
    j=j+stride
  enddo
  deallocate(uvw,lagrange)
  
  if( visu .and. rankWorld==0 )then
    print '()'
    j=0
    do iNod=1,nNod
      print '("iMod=",i3," myValues         =",4(f12.5,1x),t100,"@rkw",i3)',iNod,myValues(j+1:j+stride),rankWorld
      j=j+stride
    enddo
  endif
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  allocate(linkValues(1:stride*linkVertSize))
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  if( rankWorld==0 )print '(/"cwipi_exchange_f")'
  
  call cwipi_exchange_f(                       &
  &    couplingName=          "testPiPj"      ,&
  &    exchangeName="exch1_"//"testPiPj"      ,&
  &    exchangeDim=stride                     ,&  ! scalar
  &    ptInterpolationFct=userInterpolation   ,&  ! utilisation de la procedure plug
  &                                            &
  &    sendingFieldName="mySolu"              ,&  ! solution calculee localement
  &    sendingField=myValues                  ,&
  &                                            &
  &    receivingFieldName="linkSolu"          ,&
  &    receivingField=linkValues              ,&  ! solution de raccord
  &                                            &
  &    nStep=1                                ,&  ! pas utilisee juste pour visu cwipi
  &    timeValue=0d0                          ,&  ! pas utilisee juste pour visu cwipi
  &    nNotLocatedPoints=notLocatedPoints     ,&
  &    status=iErr                             )
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  if( rankWorld==0 )print '(/"Delete coupling")'
  call cwipi_delete_coupling_f("testPiPj")
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  if( rankWorld==0 )print '(/"Controling Coupling Results")'
  
  deltaMax=-1d0
  j=0
  k=0
  do iVert=1,linkVertSize
    delta=norm2(linkVert(j+1:j+3)-linkValues(k+1:k+3)) !+( real(rankWorld,kind=8)-linkValues(k+4) )**2
    if( deltaMax<delta )deltaMax=delta
   !print'("Delta=",e22.15)',delta
    j=j+3
    k=k+4
  enddo
  delta=deltaMax
  
  call mpi_allreduce(delta,deltaMax,1,mpi_real8,mpi_max,commWorld,iErr)
  
  if( rankWorld==0 )then
    print '(/"deltaMax=",e22.15/)',deltaMax
  endif
  
  if( deltaMax<1d-12 )then
    if( rankWorld==0 )print '(/"SUCCESS: fortran_surf_TetraP1_PiPj"/)'
  else
    if( rankWorld==0 )print '(/"FAILED: fortran_surf_TetraP1_PiPj"/)'
    stop
  endif
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  deallocate(vertx,tetra,trian)
  deallocate(linkVert,linkValues)
  deallocate(vertices,connecIndex,connec,myValues)
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  call cwipi_finalize_f()
  call mpi_finalize(iErr)
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
end program testf
