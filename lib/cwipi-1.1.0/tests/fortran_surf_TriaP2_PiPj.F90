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


!  mpirun -n 1 ./fortran_surf_TetraP2_PiPj : -n 1 ./fortran_surf_TetraP2_PiPj
!  mpirun -n 1 tests/fortran_surf_TetraP2_PiPj : -n 1 tests/fortran_surf_TetraP2_PiPj

module additionnal_Functions

contains
  
  subroutine triangleP2UV_to_TetraP2UVW(iTrian,uv,uvw)
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    !> Cette procedure retourne les coordonnées barycentriques dans le tetraP2 des sommets de la face iTrian du tetraP2
    !> entree uv (1,2,;) triangleP2
    !> sortie uvw(1,3,;) tetraP2
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    use baseSimplex2D, only: setT3MeshBasis_P2
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    integer         , intent(in)  :: iTrian
    real(8), pointer, intent(in)  :: uv (:,:)
    real(8), pointer, intent(out) :: uvw(:,:)
    !>
    integer                       :: iMod,jMod,nMod
    integer                       :: iNod     ,nNod
    integer                       :: nodes(1:6)
    real(8)                       :: TetraP2(1:3,1:10)
    real(8)                       :: TrianP2(1:3,1:06)
    real(8)                       :: u,v
    real(8), pointer              :: ai(:,:)
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    !> uvw du TetraP2
    TetraP2(1:3,01)=[0d+0,0d+0,0d+0]
    TetraP2(1:3,02)=[1d+0,0d+0,0d+0]
    TetraP2(1:3,03)=[0d+0,1d+0,0d+0]
    TetraP2(1:3,04)=[0d+0,0d+0,1d+0]
    !>
    TetraP2(1:3,05)=[5d-1,0d+0,0d+0]
    TetraP2(1:3,06)=[5d-1,5d-1,0d+0]
    TetraP2(1:3,07)=[0d+0,5d-1,0d+0]
    TetraP2(1:3,08)=[0d+0,0d+0,5d-1]
    TetraP2(1:3,09)=[5d-1,0d+0,5d-1]
    TetraP2(1:3,10)=[0d+0,5d-1,5d-1]
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    !> slection des noeuds
    select case(iTrian)
    case(1) ; nodes(1:6)=[02,03,04, 06,10,09]
    case(2) ; nodes(1:6)=[01,03,02, 07,06,05]
    case(3) ; nodes(1:6)=[01,04,03, 08,10,07]
    case(4) ; nodes(1:6)=[01,02,04, 05,09,08]
    case default ; stop "stop @ triangleP2UV_to_TetraP2UVW"
    end select
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    !> uvw du TriangleP2 face iTrian
    do iMod=1,6
      jMod=nodes(iMod)
      TrianP2(1:3,iMod)=TetraP2(1:3,jMod)
    enddo
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    !> Calcul de uvw dans le TetraP2 correspond a uv dans le TriangleP2
    nNod=size(uv,2)
    nMod=6
    
    allocate(ai(1:nMod,1:nNod))
    call setT3MeshBasis_P2(uv=uv,ai=ai)    
    
    allocate(uvw(1:3,1:nNod))
    
    uvw(1:3,1:nNod)=matmul(TrianP2(1:3,1:nMod),ai(1:nMod,1:nNod))
    
    !do iNod=1,size(uv,2)
    !  !u=uv(1,iNod)
    !  !v=uv(2,iNod)
    !  !call setT3MeshBasis_P2(uv=uv(1;2,iVert),ai=ai)
    !  uvw(1:3,iNod)=0d0
    !  do iMod=1,6
    !    uvw(1:3,iNod)=uvw(1:3,iNod)+ai(iMod,iNod)*TrianP2(1:3,iMod)
    !  enddo
    !enddo
    
    
    deallocate(ai)    
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    return
  end subroutine triangleP2UV_to_TetraP2UVW

!   subroutine setT3MeshBasis_P1(u,v,ai)
!     !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!     !> Numerotation des sommets
!     !>   3
!     !>   1 2
!     !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
!     !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!     ! delcaration des variables passees en argument
!     real(8), intent(in)    :: u,v
!     real(8), intent(inout) :: ai(:)
!     !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
!     !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!     ai(1)=1d0-u-v
!     ai(2)=    u
!     ai(3)=      v
!    !write(*,'("u,v=",2(f12.5,1x),"li=",3(f12.5,1x))')u,v,ai(1:3)
!     !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
!     return
!   end subroutine setT3MeshBasis_P1
!   
!   subroutine setT3MeshBasis_P2(u,v,ai)
!     !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!     !> Numerotation des sommets
!     !>   3
!     !>   6 5
!     !>   1 4 2
!     !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
!     !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!     ! delcaration des variables passees en argument
!     real(8), intent(in)    :: u,v
!     real(8), intent(inout) :: ai(:)
!     !>
!     real(8)                :: w
!     real(8)                :: u2,v2,w2
!     !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
!     !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!     w=1d0-u-v
!     u2=2d0*u ; v2=2d0*v ; w2=2d0*w
!     !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
!     !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!     ai(1)=w*(-1d0+w2)    !> (i,j,k)=(0,0,2)
!     ai(2)=u*(-1d0+u2)    !> (i,j,k)=(2,0,0)
!     ai(3)=v*(-1d0+v2)    !> (i,j,k)=(0,2,0)
!     ai(4)=u2*w2          !> (i,j,k)=(1,0,1)
!     ai(5)=u2*v2          !> (i,j,k)=(1,1,0)
!     ai(6)=v2*w2          !> (i,j,k)=(0,1,1)
!    !write(*,'("u,v=",2(f12.5,1x),"li=",6(f12.5,1x))')u,v,ai(1:6)
!     !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
!     return
!   end subroutine setT3MeshBasis_P2

  subroutine setT4MeshBasisP1(u,v,w,ai)
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    !> Numerotation des sommets
    !> 01 03  04
    !> 02
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    ! delcaration des variables passees en argument
    real(8), intent(in)    :: u,v,w
    real(8), intent(inout) :: ai(:)
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    ai(1)=1d0-u-v-w
    ai(2)=    u
    ai(3)=      v
    ai(4)=        w
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    return
  end subroutine setT4MeshBasisP1
  
  subroutine setT4MeshBasisP2(u,v,w,ai)
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    !> Numerotation des sommets
    !> 01 07 03  08 10  04
    !> 05 06     09
    !> 02
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    ! delcaration des variables passees en argument
    real(8), intent(in)    :: u,v,w
    real(8), intent(inout) :: ai(:)
    !>
    real(8)                :: x
    real(8)                :: u2,v2,w2,x2
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    x=1d0-u-v-w
    u2=2d0*u ; v2=2d0*v ; w2=2d0*w ; x2=2d0*x
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>   
    ai(01)=(-1d0+x2)*x  !> (i,j,k)=(000)
    ai(02)=(-1d0+u2)*u  !> (i,j,k)=(200)
    ai(03)=(-1d0+v2)*v  !> (i,j,k)=(020) 
    ai(04)=(-1d0+w2)*w  !> (i,j,k)=(002)
    ai(05)=u2      *x2  !> (i,j,k)=(100)
    ai(06)=u2*v2        !> (i,j,k)=(110)
    ai(07)=   v2   *x2  !> (i,j,k)=(010)
    ai(08)=      w2*x2  !> (i,j,k)=(001)
    ai(09)=u2   *w2     !> (i,j,k)=(101)    
    ai(10)=   v2*w2     !> (i,j,k)=(011) 
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    return
  end subroutine setT4MeshBasisP2  

end module additionnal_Functions


module variablesCommunes
  logical :: visu=.true.
  integer :: commWorld,rankWorld,sizeWorld
  integer :: commLocal,rankLocal,sizeLocal
  integer :: compOrder
  integer :: meshOrder=2
end module variablesCommunes


subroutine  userInterpolation                        ( &
  &           entitiesDim                             ,&
  &           order                                   ,&
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
  &           dist_uvw                                ,& !> new cwipi
  &                                                    &
  &           stride                                  ,&  ! =ker(calc)
  &           solverType                              ,&
  &           localField                              ,&  !   mySolu
  &           distantField                             )  ! linkSolu
  !---
  use iso_c_binding, only: c_loc,c_f_pointer
  use cwipi
  use modDeterminant
  use baseSimplex3D
  
  use additionnal_Functions
  
  use variablesCommunes
  !---
  implicit none
  !---
  integer :: entitiesDim
  integer :: order
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
  real(8) :: dist_uvw                                (*)
  integer :: stride
  integer :: solverType
  real(8) ::   localField                            (*)
  real(8) :: distantField                            (*)
  !>
  integer          :: i,j,k,iRank,iErr
  integer          :: iNod,nNod,iMod,nMod
  real(8), pointer :: uvw (:,:),rst(:,:),a(:),b(:),c(:),vand(:,:)
  integer          :: iDistantPoint
  integer          :: iTrian,iVert,iBary
  real(8), pointer :: uv(:,:),uvwOut(:,:),lagrange(:,:)
  real(8)          :: lagrangeMesh(1:10)
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  call mpi_barrier(commWorld,iErr)
  if( rankWorld==0 )print'(/3x,">>> userInterpolation")'  
  call mpi_barrier(commWorld,iErr)
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  !> Control localCoordinates localConnectivity
  call mpi_barrier(commWorld,iErr)
  if( rankWorld==0 )print'(/3x,"Control: localCoordinates, localConnectivity")'  
  call mpi_barrier(commWorld,iErr)
  if( 0==1 )then
 !if( visu )then
    do iRank=0,sizeWorld-1
      if( iRank==rankWorld )then
        print '(/3x,"meshOrder=",i1," compOrder=",i2,t120,"@rkw",i3)',meshOrder,compOrder,rankWorld
        iVert=0
        do i=1,nLocalVertex
          print '(6x,"localCoordinates (",i2,")=",3(e22.15,1x))',i,localCoordinates(iVert+1:iVert+3)
          iVert=iVert+3
        enddo
        do i=1,nLocalElement
          print '(6x,"localConnectivity(",i2,")=",*(i3,1x))',i,localConnectivity(localConnectivityIndex(i)+1:localConnectivityIndex(i+1))
        enddo
      endif
      call mpi_barrier(commWorld,iErr)
    enddo
  endif
  call mpi_barrier(commWorld,iErr)
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  !> Control localField
  call mpi_barrier(commWorld,iErr)
  if( rankWorld==0 )print'(/3x,"Control: localField")'
  call mpi_barrier(commWorld,iErr)
 !if( visu )then
  if( 0==1 )then
    do iRank=0,sizeWorld-1
      if( iRank==rankWorld )then
        print '(/3x,"meshOrder=",i1," compOrder=",i2,t120,"@rkw",i3)',meshOrder,compOrder,rankWorld
        nMod=(compOrder+1)*(compOrder+2)*(compOrder+3)/6
        j=0
        do iMod=1,nMod
          print '(6x,"localField(",i3,")=",*(e22.15,1x))',iMod,localField(j+1:j+stride)
          j=j+stride
        enddo
      endif
      call mpi_barrier(commWorld,iErr)
    enddo
    call mpi_barrier(commWorld,iErr)
  endif
  call mpi_barrier(commWorld,iErr)
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  !> Control disPtsCoordinates,disPtsLocation  
  call mpi_barrier(commWorld,iErr)
  if( rankWorld==0 )print'(/3x,"Control: disPtsCoordinates,disPtsLocation")'
  call mpi_barrier(commWorld,iErr)
  if( visu )then
    do iRank=0,sizeWorld-1
      if( iRank==rankWorld )then    
        print '(/3x,"meshOrder=",i1," compOrder=",i2,t120,"@rkw",i3)',meshOrder,compOrder,rankWorld
        iVert=0
        do iDistantPoint=1,nDistantPoint
          if(  iRank==1 .and. iDistantPoint==29 )print '(6x,"disPtsCoordinates(",i3,")=",3(e22.15,1x)," inside Cell: ",i3)',&
          & iDistantPoint,disPtsCoordinates(iVert+1:iVert+3),disPtsLocation(iDistantPoint)
          iVert=iVert+3
        enddo
      endif
      call mpi_barrier(commWorld,iErr)
    enddo
    call mpi_barrier(commWorld,iErr)
  endif
  call mpi_barrier(commWorld,iErr)
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  !> dist_uvw(:) -> uv(1:2,:) (sans dupliquer le bloc mémoire)
  !> la commmande c_f_pointer a besoin du use iso_c_binding, only: c_loc,c_f_pointer
  
  call c_f_pointer(cptr=c_loc(dist_uvw), fptr=uv, shape=[2,nDistantPoint])
  
  
  !if( rankWorld==1 )then
  !  print '(6x,"uv (1:2,11)=",2(e22.15,1x), " -> "$)',uv(1:2,11)  
  !  uv(1:2,11)=[0.966503472726928d-01,0.248894314819442d+00]
  !  print '(2(e22.15,1x))',uv(1:2,11)
  !  print '(6x,"uv (1:2,29)=",2(e22.15,1x), " -> "$)',uv(1:2,29)  
  !  uv(1:2,29)=[0.654455337907865d+00,0.248894314819442d+00]
  !  print '(2(e22.15,1x))',uv(1:2,29)
  !endif

  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  !> Control dist_uvw,uv
  call mpi_barrier(commWorld,iErr)
  if( rankWorld==0 )print'(/3x,"Control: dist_uvw,disPtsLocation")'
  call mpi_barrier(commWorld,iErr)
  if( visu )then
    do iRank=0,sizeWorld-1
      if( iRank==rankWorld )then    
        print '(/3x,"meshOrder=",i1," compOrder=",i2," stride_uvw=",i2,t120,"@rkw",i3)',meshOrder,compOrder,entitiesDim,rankWorld
        call mpi_barrier(commWorld,iErr)
        iVert=0
        do iDistantPoint=1,nDistantPoint
         !print '(6x,"dist_uvw(",i3,")=",*(e22.15,1x))',iDistantPoint,dist_uvw(iVert+1:iVert+entitiesDim)
          if( iRank==1 .and. iDistantPoint==29 )print '(6x,"uv (1:2,",i3,")=",2(e22.15,1x))',iDistantPoint,uv(1:2,iDistantPoint)
          iVert=iVert+entitiesDim
        enddo
      endif
      call mpi_barrier(commWorld,iErr)
    enddo
    call mpi_barrier(commWorld,iErr)
  endif
  call mpi_barrier(commWorld,iErr)
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  !> uv TriangleP2 -> uvw tetraP2
  !> On passe de la face TriangleP2 au TetraP2 afin de pouvoir calculer des gradients du champs
  !> cette étape n'est pas nécessaire si aucun gradient est calculé
  
  allocate(uvw(1:3,1:nDistantPoint))
  
  !> rkw0:Triangle3 <=> rkw1:Triangle4
  select case(rankWorld)
  case(0) ; iTrian=3  !> on se couple sur le triangle 3
  case(1) ; iTrian=4  !> on se couple sur le triangle 4
  end select
  
  call triangleP2UV_to_TetraP2UVW(iTrian=iTrian,uv=uv,uvw=uvw)
  
  call mpi_barrier(commWorld,iErr)
  if( rankWorld==0 )print'(/3x,"Calcul: uv TriangleP2 -> uvw tetraP2")'
  call mpi_barrier(commWorld,iErr)
  if( visu )then
    do iRank=0,sizeWorld-1
      if( iRank==rankWorld )then    
        print '(/3x,"meshOrder=",i1," compOrder=",i2,t120,"@rkw",i3)',meshOrder,compOrder,rankWorld
        call mpi_barrier(commWorld,iErr)
        do iDistantPoint=1,nDistantPoint
          if( iRank==1 .and. iDistantPoint==29 )print '(6x,"uvw(",i3,")=",3(e22.15,1x))',iDistantPoint,uvw(1:3,iDistantPoint)
        enddo
      endif
      call mpi_barrier(commWorld,iErr)
    enddo
    call mpi_barrier(commWorld,iErr)
  endif
  call mpi_barrier(commWorld,iErr)        
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  !> Base fonctionelle d'ordre compOrder
  
  nMod=(compOrder+1)*(compOrder+2)*(compOrder+3)/6  !> Tetra
  
  !> transpose = .true. => lagrange(1:nMod,1:nDistantPoint)
  call lagrange3Dv(ord=compOrder,uvwOut=uvw,lagrange=lagrange,transpose=.true.)
  
!  !> Points d'interpolation
!  call nodes3D   (ord=compOrder,uvw=uvw,display=.false.)
!  call nodes3Dopt(ord=compOrder,uvw=uvw,display=.false.)
!  !> Calcul de vand(:,:)
!  call nodes3Duvw2abc(uvw=uvw,a=a,b=b,c=c,display=.false.)
!  call vandermonde3D(ord=compOrder,a=a,b=b,c=c,vand=vand)
!  deallocate(uvw,a,b,c)
!  !> Calcul des polonômes de Lagrange d'ordre compOrder en uvwOut
!  allocate(lagrange(1:nMod,1:nNod))
!  call nodes3Duvw2abc(uvw=uvwOut,a=a,b=b,c=c,display=.false.)
!  call lagrange3Dv(ord=compOrder,vand=vand,a=a,b=b,c=c,lx=lagrange,transpose=.true.)  !> lagrange= Inverse[Transpose[Vand]].Psi[xyzOut] lagrange(nPt,np)  
  
  call mpi_barrier(commWorld,iErr)
  if( rankWorld==0 )print'(/3x,"Calcul: Bases de lagrange")'
  call mpi_barrier(commWorld,iErr)
  if( visu )then
    do iRank=0,sizeWorld-1
      if( iRank==rankWorld )then    
        !print '(/3x,"meshOrder=",i1," compOrder=",i2,t120,"@rkw",i3)',meshOrder,compOrder,rankWorld
        call mpi_barrier(commWorld,iErr)
        do iDistantPoint=1,nDistantPoint
          !print '(6x,"lagrange(",i3,")=",*(e22.15,1x))',iDistantPoint,lagrange(1:nMod,iDistantPoint)
        enddo
      endif
      call mpi_barrier(commWorld,iErr)
    enddo
    call mpi_barrier(commWorld,iErr)
  endif
  call mpi_barrier(commWorld,iErr)  
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
  
  call mpi_barrier(commWorld,iErr)
  if( rankWorld==0 )print'(/3x,"Calcul: distantField")'
  call mpi_barrier(commWorld,iErr)
  if( visu )then
    do iRank=0,sizeWorld-1
      if( iRank==rankWorld )then    
        print '(/3x,"meshOrder=",i1," compOrder=",i2,t120,"@rkw",i3)',meshOrder,compOrder,rankWorld
        call mpi_barrier(commWorld,iErr)
        j=0
        do iDistantPoint=1,nDistantPoint
          if( iRank==1 .and. iDistantPoint==29 )print '(6x,"distantField(",i3,")=",4(e22.15,1x),t120,"@rkw",i3)',&
          & iDistantPoint,distantField(j+1:j+stride),rankWorld
          j=j+stride
        enddo
      endif
      call mpi_barrier(commWorld,iErr)
    enddo
    call mpi_barrier(commWorld,iErr)
  endif
  call mpi_barrier(commWorld,iErr)    
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  deallocate(lagrange)
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  call mpi_barrier(commWorld,iErr)
  call mpi_barrier(commWorld,iErr)
  if( rankWorld==0 )print'("<<< userInterpolation rankWorld=",i2)',rankWorld
  call mpi_barrier(commWorld,iErr)
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  return
end subroutine userInterpolation


program fortran_surf_TriaP2_PiPj
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  use iso_fortran_env
  
#ifdef CWP_HAVE_FORTRAN_MPI_MODULE  
  use mpi
#endif
  use cwipi
  
  use modDeterminant
  use baseSimplex2D
  use baseSimplex3D
  use table_tet_mesh
  
  use additionnal_Functions
  
  use variablesCommunes
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  implicit none

#ifndef CWP_HAVE_FORTRAN_MPI_MODULE  
  include "mpif.h"
#endif  
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  interface
    subroutine  userInterpolation                      ( &
    &           entitiesDim                             ,&
    &           order                                   ,&
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
    &           dist_uvw                                ,&
    &                                                    &
    &           stride                                  ,&  ! =ker(calc)
    &           solverType                              ,&
    &           localField                              ,&  !   mySolu
    &           distantField                             )  ! linkSolu
    !---
    integer :: entitiesDim
    integer :: order
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
    real(8) :: dist_uvw(*)
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
  
  logical          :: droit=.false.
 !logical          :: droit=.true.
  integer          :: iVert,nVert
  real(8), pointer :: vertx(:,:)
  integer          :: iTetra,nTetra
  integer, pointer :: tetra(:,:)
  integer          :: iTrian,nTrian
  integer, pointer :: trian(:,:)
  
  integer          :: j,k
  integer          :: iMod,nMod
  integer          :: iNod,nNod
  integer          :: iCell,nCell
  real(8), pointer :: vertices   (:)
  integer, pointer :: connec     (:)
  integer, pointer :: ijk     (:)
  integer, pointer :: connecIndex(:)
  integer, pointer :: tetraNodes(:,:)
  
  real(8), pointer :: lagrange(:,:)
  real(8)          :: lagrangeMesh(1:10)
  
  real(8)          :: xyz(1:3)
  integer          :: linkVertSize
  real(8), pointer :: linkVert(:)
  integer          :: notLocatedPoints
  
  integer          :: stride
  real(8), pointer ::   myValues(:)
  real(8), pointer :: linkValues(:)
  
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
  if( rankWorld==0) print '(/"START: fortran_surf_TetraP2_PiPj")'
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  !> Initialisation de l'interface de couplage
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
  case(0) ; compOrder=07 !07 !07
  case(1) ; compOrder=10 !07 !10
  end select
  
  call mpi_barrier(commWorld,iErr)
  print '("fortran_surf_PiPj: meshOrder=",i2," compOrder=",i2,t120,"@rkw",i3)',meshOrder,compOrder,rankWorld
  call mpi_barrier(commWorld,iErr)
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  !> Create coupling
  
  call mpi_barrier(commWorld,iErr)
  if( rankWorld==0 )print '(/"Create Coupling")'
  call mpi_barrier(commWorld,iErr)
  
  call cwipi_create_coupling_f(                  &
  &    couplingName="testPiPj"                  ,&
  &    couplingType=cwipi_cpl_parallel_with_part,&
  &    cplAppli=codeCoupledName                 ,&
  &    entitiesDim=2                            ,& !> Nature du couplage surfacique
  &    tolerance=1d-1                           ,& !> Tolerance geometrique 1d-1 par defaut
  &    meshT=cwipi_static_mesh                  ,&
  &    solvert=cwipi_solver_cell_vertex         ,&
  &    outputfreq=-1                            ,& !> Frequence du post-traitement
  &    outputfmt="Ensight Gold"                 ,&
  &    outputfmtopt="binary"                     )
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  !> Create Geometric Mesh
  
  call mpi_barrier(commWorld,iErr)
  if( rankWorld==0 )print '(/"Creation des maillages géométriques (inria mesh format)")'
  call mpi_barrier(commWorld,iErr)
  
  !  04        niveau2
  !
  !  10
  !  08 09     niveau1
  !
  !  03
  !  07 06
  !  01 05 02  niveau0
  
  ! On propose de couple deux tetras en contact suivant le plan y=0
  
  !> rk0-sd3
  !> 03
  !> 07 10 
  !> 01 08 04 
  
  !> rk1-sd4
  !> 04
  !> 08 09
  !> 01 05 03
  
  !> rk0-sd3 <-> rk1-sd4
  !>   01    <->   01 
  !>   04    <->   04
  !>   03    <->   02
  !>   08    <->   08
  !>   10    <->   09
  !>   07    <->   05
  
  
  select case(rankWorld)
  case(0)
    !> Vertices
    nVert=10
    allocate(vertx(1:3,1:nVert))
    !>
    if( droit )then
      vertx(1:3,01)=[0.00d0, 0.00d0, 0.00d0]
      vertx(1:3,02)=[0.00d0,-1.00d0, 0.00d0]
      vertx(1:3,03)=[1.00d0, 0.00d0, 0.00d0]
      vertx(1:3,04)=[0.00d0, 0.00d0, 1.00d0]
      !>
      vertx(1:3,05)=[0.00d0,-0.50d0, 0.00d0]
      vertx(1:3,06)=[0.50d0,-0.50d0, 0.00d0]
      vertx(1:3,07)=[0.50d0, 0.00d0, 0.00d0]
      vertx(1:3,08)=[0.00d0, 0.00d0, 0.50d0]
      vertx(1:3,09)=[0.00d0,-0.50d0, 0.50d0]
      vertx(1:3,10)=[0.50d0, 0.00d0, 0.50d0]
    else
      vertx(1:3,01)=[ 0.00d0, 0.10d0,-0.50d0]  !> 01 <=> 01
      vertx(1:3,02)=[ 0.00d0,-1.00d0, 0.00d0]
      vertx(1:3,03)=[ 1.50d0, 0.00d0, 0.00d0]  !> 03 <=> 02
      vertx(1:3,04)=[-0.50d0,-0.10d0, 1.00d0]  !> 04 <=> 04
      !>
      vertx(1:3,05)=[ 0.00d0,-0.50d0, 0.00d0]
      vertx(1:3,06)=[ 0.50d0,-0.50d0, 0.00d0]
      vertx(1:3,07)=[ 0.50d0, 0.10d0, 0.00d0]  !> 07 <=> 05
     !vertx(1:3,08)=[ 0.00d0, 0.00d0, 0.50d0]  !> 08 <=> 08
      vertx(1:3,08)=[ 0.00d0, 0.10d0, 0.50d0]  !> 08 <=> 08
      vertx(1:3,09)=[ 0.00d0,-0.50d0, 0.50d0]
      vertx(1:3,10)=[ 0.60d0,-0.10d0, 0.60d0]  !> 10 <=> 09
    endif
    !> Tetrahedra
    nTetra=1
    allocate(tetra(1:11,1:nTetra)) !> 10 sommets + 1 marqueur
    tetra(1:11,1)=[01,02,03,04, 05,06,07,08,09,10, 1]
    !>  Triangles
    nTrian=4
    allocate(trian(1:07,1:nTrian)) !> 6 sommets + 1 marqueur
    trian(1:7,1)=[02,03,04, 06,10,09, 1]
    trian(1:7,2)=[01,03,02, 07,06,05, 1]
    trian(1:7,3)=[01,04,03, 08,10,07, 3] !> Couplage
    trian(1:7,4)=[01,02,04, 05,09,08, 1]
  case(1)
    !> Vertices
    nVert=10
    allocate(vertx(1:3,1:nVert))
    if( droit )then
      vertx(1:3,01)=[0.00d0, 0.00d0, 0.00d0]
      vertx(1:3,02)=[1.00d0, 0.00d0, 0.00d0]
      vertx(1:3,03)=[0.00d0, 1.00d0, 0.00d0]
      vertx(1:3,04)=[0.00d0, 0.00d0, 1.00d0]
      !>
      vertx(1:3,05)=[0.50d0, 0.00d0, 0.00d0]
      vertx(1:3,06)=[0.50d0, 0.50d0, 0.00d0]
      vertx(1:3,07)=[0.00d0, 0.50d0, 0.00d0]
      vertx(1:3,08)=[0.00d0, 0.00d0, 0.50d0]
      vertx(1:3,09)=[0.50d0, 0.00d0, 0.50d0]
      vertx(1:3,10)=[0.00d0, 0.50d0, 0.50d0]
    else
      vertx(1:3,01)=[ 0.00d0, 0.10d0,-0.50d0]  !> 01 <=> 01
      vertx(1:3,02)=[ 1.50d0, 0.00d0, 0.00d0]  !> 02 <=> 03
      vertx(1:3,03)=[ 0.00d0, 1.00d0, 0.00d0]
      vertx(1:3,04)=[-0.50d0,-0.10d0, 1.00d0]  !> 04 <=> 04
      !>
      vertx(1:3,05)=[ 0.50d0, 0.10d0, 0.00d0]  !> 05 <=> 07
      vertx(1:3,06)=[ 0.50d0, 0.50d0, 0.00d0]
      vertx(1:3,07)=[ 0.00d0, 0.50d0, 0.00d0]
      vertx(1:3,08)=[ 0.00d0, 0.10d0, 0.50d0]  !> 08 <=> 08
      vertx(1:3,09)=[ 0.60d0,-0.10d0, 0.60d0]  !> 09 <=> 10
      vertx(1:3,10)=[ 0.00d0, 0.50d0, 0.50d0]
    endif
    !> Tetrahedra
    nTetra=1
    allocate(tetra(1:11,1:nTetra)) !> 10 sommets + 1 marquer
    tetra(1:11,1)=[01,02,03,04, 05,06,07,08,09,10, 1]
    !>  Triangles
    nTrian=4
    allocate(trian(1:7,1:nTrian)) !> 6 sommets + 1 marquer
    trian(1:7,1)=[02,03,04, 06,10,09, 1]
    trian(1:7,2)=[01,03,02, 07,06,05, 1]
    trian(1:7,3)=[01,04,03, 08,10,07, 1]
    trian(1:7,4)=[01,02,04, 05,09,08, 3] !> Couplage    
  end select
  
  
  if( visu )then
    !> Ecriture des maillages au format mesh de l'inria
    if( rankWorld==0)then
      print '(/3x,"Ecriture de : Tetra0.mesh")'
      print '( 3x,"Ecriture de : Tetra1.mesh")'
    endif
    !>
    write(meshName,'("Tetra",i1,".mesh")')rankWorld
    open(unit=100,file=trim(meshName),action='write',status='unknown')
    write(100,'("MeshVersionFormatted 2"/)')
    write(100,'("Dimension 3"/)')
    write(100,'("Vertices")')
    write(100,'(i2)')nVert
    do iVert=1,nVert
      write(100,'(3(e22.15,1x),i2)')vertx(1:3,iVert),0
    enddo
    write(100,'(/"TetrahedraP2")')
    write(100,'(i1)')nTetra
    do iTetra=1,nTetra
      write(100,'(*(i6,1x))')tetra(:,iTetra)
    enddo
    write(100,'(/"TrianglesP2")')
    write(100,'(i1)')nTrian
    do iTrian=1,nTrian
      write(100,'(*(i6,1x))')trian(:,iTrian)
    enddo
    write(100,'(/"End")')
    close(100)
  endif
  
  !> On se couple sur le triangle commun aux deux tetraP2 (y=0) qui va servir de maillage pour le couplage
  
  !  03
  !  06 05
  !  01 04 02
  
  !> rkw0:Triangle3 <=> rkw1:Triangle4
  select case(rankWorld)
  case(0) ; iTrian=3  !> on se couple sur le triangle 3
  case(1) ; iTrian=4  !> on se couple sur le triangle 4
  end select
  
  if( 0==1 )then !> Ancien cwipi
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    if( rankWorld==0)print '(/"Using Cwipi"/)'
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    !> On degrade le maillage à l'ordre 1
    nVert=03
    nCell=01
    allocate( vertices   (1:3*nVert)    )  !> sommets
    allocate( connec     (1:3*nCell)    )  !> triangle P1
    allocate( connecIndex(1:nCell+1)    )  !> triangle
    connec(1:3)=[1,2,3]
    connecIndex(1:2)=[0,3]
    
    j=0
    do iNod=1,nVert
      vertices(j+1:j+3)=vertx(1:3,trian(iNod,iTrian))
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
    if( visu )then
      write(meshName,'("Triangle",i1,".mesh")')rankWorld
      open(unit=100,file=trim(meshName),action='write',status='unknown')
      write(100,'("MeshVersionFormatted 2"/)')
      write(100,'("Dimension 3"/)')
      write(100,'("Vertices")')
      write(100,'(i2)')nVert
      j=0
      do iVert=1,nVert
        write(100,'(3(e22.15,1x),i2)')vertices(j+1:j+3),0
        j=j+3
      enddo
      write(100,'(/"Triangles")')
      write(100,'(i1)')nCell
      do iCell=1,nCell
        write(100,'(*(i6,1x))')connec( connecIndex(iCell)+1:connecIndex(iCell+1) ),0
      enddo
      write(100,'(/"End")')
      close(100)
    endif
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
  else !> Nouveau cwipi
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    if( rankWorld==0)print '(/"Using Cwipi High Order"/)'
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    !> On conserve un maillage d'ordre 2
    nVert=06
    nCell=01
    allocate( vertices   (1:3*nVert) )  !> sommets
    allocate( connec     (1:6*nCell) )  !> triangle P2
    allocate( connecIndex(1:nCell+1) )  !> triangle
    connec(1:6)=[1,2,3,4,5,6]
    connecIndex(1:2)=[0,6]
    
    j=0
    do iNod=1,nVert
      vertices(j+1:j+3)=vertx(1:3,trian(iNod,iTrian))
      j=j+3
    enddo
    
    if( visu )then
      do iRank=0,sizeWorld-1
        if( iRank==rankWorld )then
         !print '(3x,"rankWorld=",i1," meshOrder=",i1)',rankWorld,meshOrder    
          print '(3x,"meshOrder=",i1," compOrder=",i2,t120,"@rkw",i3)',meshOrder,compOrder,rankWorld
          do iVert=1,nVert
            print '(3x,"iVert=",i3," xyz=",*(e22.15,1x))',iVert,vertices(3*(iVert-1)+1:3*iVert)
          enddo
          do iCell=1,nCell
            print '(3x,"iCell=",i3," cel=",*(i3,1x))',iCell,connec( connecIndex(iCell)+1:connecIndex(iCell+1) )
          enddo
        endif
        call mpi_barrier(commWorld,iErr)
      enddo
    endif
    
    !> Transmission des maillages à cwipi
    call cwipi_ho_define_mesh_f(  & !> NEW Cwipi
    &   couplingName="testPiPj"  ,&
    &   nVertex     =nVert       ,&
    &   nElts       =nCell       ,&
    &   order       =meshOrder   ,&  
    &   coords      =vertices    ,&
    &   connecIndex =connecIndex ,&
    &   connec      =connec       )
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    !> Definition du Triangle P2
    
    ! 3
    ! 6 5
    ! 1 4 2
    allocate( ijk(1:12) )  !> sommets
    ijk( 1: 2)=[0,0] !> 1
    ijk( 3: 4)=[2,0] !> 2
    ijk( 5: 6)=[0,2] !> 3
    ijk( 7: 8)=[1,0] !> 4
    ijk( 9:10)=[1,1] !> 5
    ijk(11:12)=[0,1] !> 6
    
    call cwipi_ho_ordering_from_IJK_set_f( & !> NEW Cwipi
    &   couplingName ="testPiPj"          ,&
    &   tElt         = CWIPI_FACE_TRIAHO  ,&
    &   nNodes       = 6                  ,&
    &   IJK          = ijk                 )
    
    deallocate (ijk)
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    if( visu )then
      
      if( rankWorld==0)then
        print '(/"Writing coupling mesh file: Triangle0.mesh")'
        print '( "Writing coupling mesh file: Triangle1.mesh")'
      endif
      
      write(meshName,'("Triangle",i1,".mesh")')rankWorld
      open(unit=100,file=trim(meshName),action='write',status='unknown')
      write(100,'("MeshVersionFormatted 1"/)')
      write(100,'("Dimension 3"/)')
      write(100,'("Vertices")')
      write(100,'(i2)')nVert
      j=0
      do iVert=1,nVert
        write(100,'(3(e22.15,1x),i2)')vertices(j+1:j+3),0
        j=j+3
      enddo
      write(100,'(/"TrianglesP2")')
      write(100,'(i1)')nCell
      do iCell=1,nCell
        write(100,'(*(i6,1x))')connec( connecIndex(iCell)+1:connecIndex(iCell+1) ),0
      enddo
      write(100,'(/"End")')
      close(100)
      
    endif
    
  endif !> ordre de maillage 1 ou 2 (OLD/NEW cwipi)
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  !> Points de couplage situés Triangle avec le marqueur de peau 3 => face commune aux deux tetras <=
  
  call mpi_barrier(commWorld,iErr)
  if( rankWorld==0 ) print'(/"Calcul des coordonnees des points de couplage")'
  call mpi_barrier(commWorld,iErr)
  
  !> Calcul des coordonnees barycentriques optimisées sur face trianglulaire compOrder
  call nodes3Dopt_2D(ord=compOrder,uvw=uvw,display=.false.) !> ordre du calcul
  
  !> Visu uvw
  if( visu )then
    do iRank=0,sizeWorld-1
      call mpi_barrier(commWorld,iErr)
      if( iRank==rankWorld )then
        print '(/3x,"meshOrder=",i1," compOrder=",i2,t120,"@rkw",i3)',meshOrder,compOrder,rankWorld
        do iVert=1,size(uvw,2)
          if( iRank==0 .and. iVert==29 )print '(6x,"uvw(",i2,")=",3(e22.15,1x),"1-u-v-w=",e22.15)',iVert,uvw(1:3,iVert),1d0-uvw(1,iVert)-uvw(2,iVert)-uvw(3,iVert)
        enddo
      endif
      call mpi_barrier(commWorld,iErr)
    enddo
    call mpi_barrier(commWorld,iErr)
  endif
  
  
  
  !> Calculs des coordonnées correspondantes aux coordonnees barycentriques
  linkVertSize=size(uvw,2)
  allocate(linkVert(1:3*linkVertSize))
  
  nMod=(meshOrder+1)*(meshOrder+2)/2  !> Triangle meshOrder (geometric)
  nNod=size(uvw,2)                    !> Triangle compOrder
  
  
  
  allocate(lagrange(1:nMod,1:linkVertSize))
  call setT3MeshBasis_P2(uv=uvw,ai=lagrange)    
  
  
  
  j=0
  do iVert=1,linkVertSize
    !> Fonction
    !call setT3MeshBasis_P2(u=uvw(1,iVert),v=uvw(2,iVert),ai=lagrangeMesh)
    linkVert(j+1:j+3)=0d0
    do iMod=1,nMod
     !linkVert(j+1:j+3)=linkVert(j+1:j+3)+lagrangeMesh(iMod)*vertx(1:3,trian(iMod,iTrian))
      linkVert(j+1:j+3)=linkVert(j+1:j+3)+lagrange(iMod,iVert)*vertx(1:3,trian(iMod,iTrian))
    enddo
    j=j+3
  enddo
  
  deallocate(lagrange)
  deallocate(uvw)
  
  
  
  
  
  
  !> Visu des coordonnees de couplage
  if( visu )then
    do iRank=0,sizeWorld-1
      call mpi_barrier(commWorld,iErr)
      if( iRank==rankWorld )then
        print '(/3x,"meshOrder=",i1," compOrder=",i2,t120,"@rkw",i3)',meshOrder,compOrder,rankWorld
        j=0
        do iVert=1,linkVertSize
          if( iVert==31 ) then
            print '(6x,"linkVert(",i2,")=",3(e22.15,1x))',iVert,linkVert(j+1:j+3)
!            linkVert(1) = linkVert(j+1)
!            linkVert(2) = linkVert(j+2)
!            linkVert(3) = linkVert(j+3)
          endif
          j=j+3
        enddo
      endif
      call mpi_barrier(commWorld,iErr)
    enddo
    call mpi_barrier(commWorld,iErr)
!    linkVertSize = 1
  endif
  
  
  !> Transmission à cwipi des coordonnees de couplage
  call cwipi_set_points_to_locate_f( &
  &    couplingName="testPiPj"      ,&
  &    nPts  =linkVertSize          ,&
  &    coords=linkVert               )
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  !> Localisation par cwipi des coordonnees de couplage
  call mpi_barrier(commWorld,iErr)
  if( rankWorld==0 )print '(/"Localisation")'
  call mpi_barrier(commWorld,iErr)
  
  call cwipi_locate_f(couplingName="testPiPj")
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  
  
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  !> Initialisation of myValues
  
  !> Liste coordonées barycentriques uvw(1:3,:) du Tetra d'ordre compOrder
  call nodes3D   (ord=compOrder,uvw=uvw,display=.false.)
  if( visu )then
    !> construit une connectivite entre les sommets uvw(1:3,:)
    call driverTetMesh(ord=compOrder,node_xyz=uvw,tetra_node=tetraNodes) 
  endif
  
  !> Optimisation du placement des sommets uvw(1:3,:) (deplace les coordonnées uvw pour réduire phénomène de Runge
  call nodes3Dopt(ord=compOrder,uvw=uvw,display=.false.)
  if( visu )then
    !> sauvegarde maillage (sommets et connectivité)
    call saveTetMesh  (ord=compOrder,node_xyz=uvw,tetra_node=tetraNodes)
    deallocate(tetraNodes)
  endif
  
  stride=4 ; nNod=size(uvw,2)
  allocate(myValues(1:stride*nNod))
  
  nMod=(meshOrder+1)*(meshOrder+2)*(meshOrder+3)/6 !> Tetra meshOrder
  j=0
  do iNod=1,nNod
    call setT4MeshBasisP2(u=uvw(1,iNod),v=uvw(2,iNod),w=uvw(3,iNod),ai=lagrangeMesh)    
    xyz(1:3)=0d0
    do iMod=1,nMod
      xyz(1:3)=xyz(1:3)+lagrangeMesh(iMod)*vertx(1:3,tetra(iMod,1))
    enddo
    myValues(j+1:j+stride)=[xyz(1),xyz(2),xyz(3),real(rankWorld,kind=8)]
    j=j+stride
  enddo
  deallocate(uvw)
  
if(0==1)then
  !> Visu des valeurs de couplage
  call mpi_barrier(commWorld,iErr)
  if( rankWorld==0 )print '(/"Initialisation of myValues (x,y,z,rankWorld)")'
  call mpi_barrier(commWorld,iErr)
  if( visu )then
    do iRank=0,sizeWorld-1
      if( iRank==rankWorld )then
        print '(/3x,"Computing myValues nMod=",i3,2x,"nNod=",i3,t120,"@rkw",i3)',nMod,nNod,rankWorld
        j=0
        do iNod=1,nNod
          print '(6x,"myValues(",i3,")=",4(e22.15,1x))',iNod,myValues(j+1:j+stride)
          j=j+stride
        enddo        
      endif
      call mpi_barrier(commWorld,iErr)
    enddo
    call mpi_barrier(commWorld,iErr)
  endif
endif
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  allocate(linkValues(1:stride*linkVertSize))
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  call mpi_barrier(commWorld,iErr)
  if( rankWorld==0 )print '(/"cwipi_exchange_f")'
  call mpi_barrier(commWorld,iErr)
  
   call cwipi_exchange_f(                       &
   &    couplingName=          "testPiPj"      ,&
   &    exchangeName="exch1_"//"testPiPj"      ,&
   &    meshOrder=meshOrder                    ,&  !> NEW cwipi
   &    exchangeDim=stride                     ,&  !> scalar
   &    ptHoInterpolationFct=userInterpolation ,&  !> utilisation de la procedure plug
   &                                            &
   &    sendingFieldName="mySolu"              ,&  !> solution calculee localement
   &    sendingField=myValues                  ,&
   &                                            &
   &    receivingFieldName="linkSolu"          ,&
   &    receivingField=linkValues              ,&  !> solution de raccord
   &                                            &
   &    nStep=1                                ,&  !> pas utilisee juste pour visu cwipi
   &    timeValue=0d0                          ,&  !> pas utilisee juste pour visu cwipi
   &    nNotLocatedPoints=notLocatedPoints     ,&
   &    status=iErr                             )
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  call mpi_barrier(commWorld,iErr)
  if( rankWorld==0 )print '(/"Delete coupling")'
  call mpi_barrier(commWorld,iErr)
  call cwipi_delete_coupling_f("testPiPj")
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  call mpi_barrier(commWorld,iErr)
  if( rankWorld==0 )print '(/"Controling Coupling Results")'
  call mpi_barrier(commWorld,iErr)
  
  deltaMax=-1d0
  j=0
  k=0
  do iRank=0,sizeWorld-1
    if( iRank==rankWorld )then
      print '(/3x,"compOrder =",i2," controling linkValues - linkVertSize=",i2,t120,"@rkw",i3)',compOrder,linkVertSize,rankWorld
      do iVert=1,linkVertSize
        delta=norm2(linkVert(j+1:j+3)-linkValues(k+1:k+3)) !+( real(rankWorld,kind=8)-linkValues(k+4) )**2
        if( deltaMax<delta )deltaMax=delta
        if( 1d-08<=delta )then
          print'(6x,"iVert=",i2,1x,"linkVert=",3(e22.15,1x),"linkValues=",3(e22.15,1x),"Delta=",e22.15)',&
          & ivert,linkVert(j+1:j+3),linkValues(k+1:k+3),delta
        endif
        j=j+3
        k=k+4
      enddo
      delta=deltaMax
    endif
    call mpi_barrier(commWorld,iErr)
  enddo
  
  call mpi_allreduce(delta,deltaMax,1,mpi_real8,mpi_max,commWorld,iErr)
  
  if( rankWorld==0 )then
    print '(/"deltaMax=",e22.15/)',deltaMax
  endif
  
  if( deltaMax<1d-08 )then
    if( rankWorld==0 )print '(/"SUCCESS: fortran_surf_PiPj"/)'
  else
    if( rankWorld==0 )print '(/"FAILED: fortran_surf_PiPj"/)'
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
  
end program fortran_surf_TriaP2_PiPj
