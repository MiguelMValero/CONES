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

subroutine  interpolationbidon_f(entitiesDim, &
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
                                 disPtsBaryCoordIdx, &
                                 disPtsBaryCoord, &
                                 stride, &
                                 solverType, &
                                 localField, &
                                 distantField)

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
  integer, dimension(*) :: disPtsBaryCoordIdx
  double precision, dimension(*) :: disPtsBaryCoord
  integer :: stride
  integer :: solverType
  double precision, dimension(*) :: localField
  double precision, dimension(*) :: distantField

  integer :: i

  do i = 1, nDistantPoint
     distantField(i) = i
  enddo

end subroutine interpolationbidon_f

!
! -----------------------------------------------------------------------------
! Programme de tests
! -----------------------------------------------------------------------------
!

program testf

  use mpi
  use cwipi

  implicit none

  interface
       subroutine  interpolationbidon_f(entitiesDim, &
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
         integer, dimension(*) :: disPtsBaryCoordIdx
         double precision, dimension(*) :: disPtsBaryCoord
         integer :: stride
         integer :: solverType
         double precision, dimension(*) :: localField
         double precision, dimension(*) :: distantField
       end subroutine interpolationbidon_f
    end interface

  integer, allocatable, dimension(:) :: location
  integer, allocatable, dimension(:) :: baryCooIdx
  double precision, allocatable, dimension(:) :: baryCoo, tmpDbl
  integer :: nLocatedPoints
  integer :: nNotLocatedPoints
  integer :: nDistantPoints

  integer :: localcom, localGroup, p1Group, p1Comm
  integer :: irank, currentRank, localcommsize
  character (len = 4) :: proc
  integer :: code
  integer :: iiunit
  integer :: ivalue
  double precision :: dvalue

 double precision :: xmin = -100.d0
  double precision :: xmax =  100.d0
  double precision :: ymin = -100.d0
  double precision :: ymax =  100.d0
  integer :: nx   = 24
  integer  :: ny   = 28
  integer  :: initrandom = 2

  integer nvertex, nelts, lconnecindex

  integer, parameter :: nvertexm = 4000
  integer, parameter :: neltsm = 4000
  integer, parameter :: lconnecindexm = 12000
  integer, parameter :: nptstolocate = 21
  double precision, dimension(3*nptstolocate) :: coordstolocate

  double precision, dimension(3*nvertexm) :: coords
  integer, dimension(neltsm+1) :: connecindex
  integer, dimension(lconnecindexm) :: connec

  double precision, dimension(3*nvertexm) :: values, localvalues

  integer status

  integer i, order, k

  integer :: vpar = 10
  character (len = 6) :: cpar = "niterf"
  character (len = 30) :: disstr = ""

  integer :: stride = 1
  integer, dimension(1) :: rl
  integer :: dislocalcommsize

  call mpi_init(code)
  call mpi_comm_rank(mpi_comm_world, irank, code)


!
! -----------------------------------------
! Initialisation de l'interface de couplage
! -----------------------------------------
!

  call cwipi_init_f (mpi_comm_world, &
                     "CodeFortran", &
                     localcom)


!
! ------------------------------------------------
! Creation du fichier de sortie listing
! (en parallele, un fichier listing par processus)
! ------------------------------------------------
!

  call mpi_comm_rank(localcom, currentRank, code)
  call mpi_comm_size(localcom, localcommsize, code)

  write(proc,'(i4.4)') currentRank
  iiunit = 9
  open(unit=iiunit, file='listing_test2D_1_c2_'//proc, &
       form='formatted', status='unknown')

  call cwipi_set_output_listing_f(iiunit)

  write(iiunit,*)
  write(iiunit,*) "dump apres initialisation"
  write(iiunit,*) "-------------------------"
  write(iiunit,*)

  call cwipi_dump_appli_properties_f

!
! -------------------------------
! Test des parametres de controle
! -------------------------------
!

!
! Ajout de parametres de controle

  call cwipi_add_loc_int_ctrl_param_f("niterf", 10)
!  call cwipi_add_loc_int_ctrl_param_f(cpar, vpar)

  call cwipi_add_loc_dbl_ctrl_param_f("physicaltimef", 1.123d0)

  write(iiunit,*)
  write(iiunit,*) "dump apres ajout de parametres"
  write(iiunit,*) "------------------------------"
  write(iiunit,*)
  call cwipi_dump_appli_properties_f
!
! Modification des parametres de controle
  call cwipi_get_loc_int_ctrl_param_f("niterf", ivalue)

  ivalue = ivalue + 1

  call cwipi_set_loc_int_ctrl_param_f("niterf", ivalue)

  call cwipi_get_loc_dbl_ctrl_param_f("physicaltimef", dvalue)

  dvalue = dvalue + 0.1d0

  call cwipi_set_loc_dbl_ctrl_param_f("physicaltimef", dvalue)

  write(iiunit,*)
  write(iiunit,*) "dump apres modification de parametres"
  write(iiunit,*) "-------------------------------------"
  write(iiunit,*)
  call cwipi_dump_appli_properties_f
!
! Echange des parametres de controle

  call cwipi_synch_ctrl_param_f("CodeC")
  disstr = ""
  call cwipi_get_dis_str_ctrl_param_f("CodeC","varname", disstr)

  write(iiunit,*) "disstr =",disstr


  write(iiunit,*)
  write(iiunit,*) "dump apres synchronisation"
  write(iiunit,*) "--------------------------"
  write(iiunit,*)
  call cwipi_dump_appli_properties_f
!
! Suppression des parametres de controle
  call cwipi_del_loc_int_ctrl_param_f("niterf")
  call cwipi_del_loc_dbl_ctrl_param_f("physicaltimef")

  write(iiunit,*)
  write(iiunit,*) "dump apres suppression des parametres"
  write(iiunit,*) "-------------------------------------"
  write(iiunit,*)
  call cwipi_dump_appli_properties_f

  call cwipi_add_loc_int_ctrl_param_f("localcommsize", localcommsize)
  call cwipi_synch_ctrl_param_f("CodeC")
  call cwipi_get_dis_int_ctrl_param_f("CodeC","localcommsize", dislocalcommsize)
  call cwipi_dump_appli_properties_f

  if ((localcommsize .eq. 1) .and. (dislocalcommsize .eq. 1)) then
     write(iiunit, *)
     write(iiunit, *) "--------------------------------------------------------"
     write(iiunit, *)
     write(iiunit, *) " Test 0"
     write(iiunit, *)

     call cwipi_create_coupling_f("test2D_0", &
          cwipi_cpl_parallel_with_part,&
          "CodeC", &
          1,     & ! Dimension des entites geometriques
          100.d0, & ! Tolerance geometrique
          cwipi_static_mesh, &
          cwipi_solver_cell_vertex, &
          1, &
          "Ensight Gold",&
          "text")

     coords(1) = 0.d0
     coords(2) = 0.d0
     coords(3) = 0.d0

     coords(4) = 0.5d0
     coords(5) = 0.d0
     coords(6) = 0.d0

     coords(7) = 2.d0
     coords(8) = 0.d0
     coords(9) = 0.d0

     nvertex = 3
     nelts = 2

     connecindex(1) = 0
     connecindex(2) = 2
     connecindex(3) = 4

     connec(1) = 1
     connec(2) = 2

     connec(3) = 2
     connec(4) = 3

     call cwipi_define_mesh_f("test2D_0", &
          nvertex, &
          nelts, &
          coords, &
          connecindex, &
          connec)


     values(1) = coords(1)
     values(2) = coords(4)
     values(3) = coords(7)

     stride = 1
     call cwipi_exchange_f ("test2D_0", &
                             "echange1", &
                             stride, &
                             1, &
                             0.1d0, &
                             "coox", &
                             values, &
                             "coox", &
                             localvalues, &
                             nNotLocatedPoints, &
                             status)
     call cwipi_delete_coupling_f("test2D_0");

  endif
!
! ------------------------
! Test couplage p1 <-> p1
! ------------------------
!
! Construction de "l'objet" couplage

  write(iiunit, *)
  write(iiunit, *) "--------------------------------------------------------"
  write(iiunit, *)
  write(iiunit, *) " Test 1"
  write(iiunit, *)

  call cwipi_create_coupling_f("test2D_1", &
                                   cwipi_cpl_parallel_with_part,&
                                   "CodeC", &
                                   2,     & ! Dimension des entites geometriques
                                   0.1d0, & ! Tolerance geometrique
                                   cwipi_static_mesh, &
                                   cwipi_solver_cell_vertex, &
                                   1, &
                                   "Ensight Gold",&
                                   "text")

!
! Construction du maillage

  nvertex      = nvertexm
  nelts        = neltsm
  lconnecindex = lconnecindexm
  order = 1

  xmin = -1.d-4
  xmax =  1.d-4
  ymin = -1.d-4
  ymax =  1.d-4
  nx = 24
  ny = 28

  call creemaillagepolygone2d_f (order, &
                                 localcom, &
                                 xmin, &
                                 xmax, &
                                 ymin, &
                                 ymax, &
                                 initrandom, &
                                 nx, &
                                 ny, &
                                 nvertex, &
                                 coords, &
                                 nelts, &
                                 lconnecindex, &
                                 connecindex, &
                                 connec)

  call cwipi_define_mesh_f("test2D_1", &
                               nvertex, &
                               nelts, &
                               coords, &
                               connecindex, &
                               connec)

!
! Envoi de la coory a codec
! Reception de la coorx provenant de codec

  do i = 1, nvertex
     values(i) = coords(3*(i-1) + 2)
  enddo

  stride = 1
  call cwipi_exchange_f ("test2D_1", &
                             "echange1", &
                             stride, &
                             1, &
                             0.1d0, &
                             "cooy", &
                             values, &
                             "coox", &
                             localvalues, &
                             nNotLocatedPoints, &
                             status)
  call printStatus(iiunit, status)

!
! Suppression de l'objet couplage "couplingcellvertex"

  call cwipi_delete_coupling_f("test2D_1")
  write(iiunit, *)
  write(iiunit, *) "--------------------------------------------------------"
  write(iiunit, *)

  write(iiunit, *)
  write(iiunit, *) "--------------------------------------------------------"
  write(iiunit, *)
  write(iiunit, *) " Test 1 bis"
  write(iiunit, *)

  call cwipi_create_coupling_f("test2D_01", &
                                   cwipi_cpl_parallel_with_part,&
                                   "CodeC", &
                                   2,     & ! Dimension des entites geometriques
                                   0.1d0, & ! Tolerance geometrique
                                   cwipi_static_mesh, &
                                   cwipi_solver_cell_vertex, &
                                   1, &
                                   "Ensight Gold",&
                                   "text")

!
! Construction du maillage

  nvertex      = nvertexm
  nelts        = neltsm
  lconnecindex = lconnecindexm
  order = 1

  xmin = -100.d0
  xmax =  100.d0
  ymin = -100.d0
  ymax =  100.d0
  nx = 24
  ny = 28


  call creemaillagepolygone2d_f (order, &
                                 localcom, &
                                 xmin, &
                                 xmax, &
                                 ymin, &
                                 ymax, &
                                 initrandom, &
                                 nx, &
                                 ny, &
                                 nvertex, &
                                 coords, &
                                 nelts, &
                                 lconnecindex, &
                                 connecindex, &
                                 connec)

  call cwipi_define_mesh_f("test2D_01", &
                               nvertex, &
                               nelts, &
                               coords, &
                               connecindex, &
                               connec)

!
! Envoi de la coory a codec
! Reception de la coorx provenant de codec

  do i = 1, nvertex
     values(i) = coords(3*(i-1) + 2)
  enddo

  stride = 1
  call cwipi_exchange_f ("test2D_01", &
                             "echange1", &
                             stride, &
                             1, &
                             0.1d0, &
                             "cooy", &
                             values, &
                             "coox", &
                             localvalues, &
                             nNotLocatedPoints, &
                             status)
  call printStatus(iiunit, status)

!
! Suppression de l'objet couplage "couplingcellvertex"

  call cwipi_delete_coupling_f("test2D_01")
  write(iiunit, *)
  write(iiunit, *) "--------------------------------------------------------"
  write(iiunit, *)
!
! -------------------------------------
! Test couplage P1 -> P0 puis P0 -> P1
! -------------------------------------
!

  write(iiunit, *)
  write(iiunit, *) "--------------------------------------------------------"
  write(iiunit, *)
  write(iiunit, *) " Test 2"
  write(iiunit, *)

  call cwipi_create_coupling_f("test2D_2", &
                                   cwipi_cpl_parallel_with_part,&
                                   "CodeC", &
                                   2,     & ! Dimension des entites geometriques
                                   0.1d0, & ! Tolerance geometrique
                                   cwipi_static_mesh, &
                                   cwipi_solver_cell_vertex, &
                                   1, &
                                   "Ensight Gold",&
                                   "text")

!
! Construction du maillage

  nvertex      = nvertexm
  nelts        = neltsm
  lconnecindex = lconnecindexm
  order = 1

  xmin = -100.d0
  xmax =  100.d0
  ymin = -100.d0
  ymax =  100.d0
  nx = 24
  ny = 28

  call creemaillagepolygone2d_f (order, &
                                 localcom, &
                                 xmin, &
                                 xmax, &
                                 ymin, &
                                 ymax, &
                                 initrandom, &
                                 nx, &
                                 ny, &
                                 nvertex, &
                                 coords, &
                                 nelts, &
                                 lconnecindex, &
                                 connecindex, &
                                 connec)

  call cwipi_define_mesh_f("test2D_2", &
                               nvertex, &
                               nelts, &
                               coords, &
                               connecindex, &
                               connec)

!
! Envoi de la coory a codec
! Reception de la coorx provenant de codec

  do i = 1, nvertex
     values(i) = coords(3*(i-1) + 2)
  enddo

  stride = 1
  call cwipi_send_f ("test2D_2", &
                         "echange1", &
                         stride, &
                         1, &
                         0.1d0, &
                         "cooY", &
                         values, &
                         status)
  call printStatus(iiunit, status)

  call cwipi_receive_f ("test2D_2", &
                            "echange2", &
                            stride, &
                            1, &
                            0.1d0, &
                            "cooYY", &
                            values, &
                            nNotLocatedPoints, &
                            status)
  call printStatus(iiunit, status)


!
! Suppression de l'objet couplage "couplingcellvertex"

  call cwipi_delete_coupling_f("test2D_2");
  write(iiunit, *)
  write(iiunit, *) "--------------------------------------------------------"
  write(iiunit, *)

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

  call cwipi_create_coupling_f("test2D_3", &
                                   cwipi_cpl_parallel_with_part,&
                                   "CodeC", &
                                   2,     & ! Dimension des entites geometriques
                                   0.1d0, & ! Tolerance geometrique
                                   cwipi_static_mesh, &
                                   cwipi_solver_cell_vertex, &
                                   1, &
                                   "Ensight Gold",&
                                   "text")

!
! Construction du maillage

  nvertex      = nvertexm
  nelts        = neltsm
  lconnecindex = lconnecindexm
  order = 1

  call creemaillagepolygone2d_f (order, &
                                 localcom, &
                                 xmin, &
                                 xmax, &
                                 ymin, &
                                 ymax, &
                                 initrandom, &
                                 nx, &
                                 ny, &
                                 nvertex, &
                                 coords, &
                                 nelts, &
                                 lconnecindex, &
                                 connecindex, &
                                 connec)

  call cwipi_define_mesh_f("test2D_3", &
                               nvertex, &
                               nelts, &
                               coords, &
                               connecindex, &
                               connec)

!
! Definition des points a localiser

  coordstolocate(1) = -75.d0
  coordstolocate(2) = -75.d0
  coordstolocate(3) = 0.d0

  coordstolocate(4) = -75.d0
  coordstolocate(5) = -50.d0
  coordstolocate(6) = 0.d0

  coordstolocate(7) = -75.d0
  coordstolocate(8) = -25.d0
  coordstolocate(9) = 0.d0

  coordstolocate(10) = -75.d0
  coordstolocate(11) = 0.d0
  coordstolocate(12) = 0.d0

  coordstolocate(13) = -75.d0
  coordstolocate(14) = 25.d0
  coordstolocate(15) = 0.d0

  coordstolocate(16) = -75.d0
  coordstolocate(17) = 50.d0
  coordstolocate(18) = 0.d0

  coordstolocate(19) = -75.d0
  coordstolocate(20) = 75.d0
  coordstolocate(21) = 0.d0

  coordstolocate(22) = -25.d0
  coordstolocate(23) = -75.d0
  coordstolocate(24) = 0.d0

  coordstolocate(25) = -25.d0
  coordstolocate(26) = -50.d0
  coordstolocate(27) = 0.d0

  coordstolocate(28) = -25.d0
  coordstolocate(29) = -25.d0
  coordstolocate(30) = 0.d0

  coordstolocate(31) = -25.d0
  coordstolocate(32) = 0.d0
  coordstolocate(33) = 0.d0

  coordstolocate(34) = -25.d0
  coordstolocate(35) = 25.d0
  coordstolocate(36) = 0.d0

  coordstolocate(37) = -25.d0
  coordstolocate(38) = 50.d0
  coordstolocate(39) = 0.d0

  coordstolocate(40) = -25.d0
  coordstolocate(41) = 75.d0
  coordstolocate(42) = 0.d0

  coordstolocate(43) = 25.d0
  coordstolocate(44) = -75.d0
  coordstolocate(45) = 0.d0

  coordstolocate(46) = 25.d0
  coordstolocate(47) = -50.d0
  coordstolocate(48) = 0.d0

  coordstolocate(49) = 25.d0
  coordstolocate(50) = -25.d0
  coordstolocate(51) = 0.d0

  coordstolocate(52) = 25.d0
  coordstolocate(53) = 0.d0
  coordstolocate(54) = 0.d0

  coordstolocate(55) = 25.d0
  coordstolocate(56) = 25.d0
  coordstolocate(57) = 0.d0

  coordstolocate(58) = 25.d0
  coordstolocate(59) = 50.d0
  coordstolocate(60) = 0.d0

  coordstolocate(61) = 25.d0
  coordstolocate(62) = 75.d0
  coordstolocate(63) = 0.d0

  call cwipi_set_points_to_locate_f ("test2D_3", &
                                         nptstolocate, &
                                         coordstolocate)

!
! Envoi de la coordonnee Y a codeC
! Reception de la coordonnee Y provenant de codec

  do i = 1, nvertex
     values(i) = coords(3*(i-1) + 2)
  enddo

  stride = 1
  call cwipi_exchange_f ("test2D_3", &
                             "echange1", &
                             stride, &
                             1, &
                             0.1d0, &
                             "cooy", &
                             values, &
                             "coox", &
                             localvalues, &
                             nNotLocatedPoints, &
                             status)
  call printStatus(iiunit, status)
  write(iiunit,*) "valeurs recues test2D_3"
  write(iiunit,*) (localvalues(i),i=1,nptstolocate)

!
! Suppression de l'objet couplage

  call cwipi_delete_coupling_f("test2D_3");
  write(iiunit, *)
  write(iiunit, *) "--------------------------------------------------------"
  write(iiunit, *)

!
! -------------------------------------
! Test de definition d'une fonction
! d'interpolation en fortran (callback)
! -------------------------------------
!
  write(iiunit, *)
  write(iiunit, *) "--------------------------------------------------------"
  write(iiunit, *)
  write(iiunit, *) " Test 4"
  write(iiunit, *)

  call cwipi_create_coupling_f("test2D_4", &
                                   cwipi_cpl_parallel_with_part,&
                                   "CodeC", &
                                   2,     & ! Dimension des entites geometriques
                                   0.1d0, & ! Tolerance geometrique
                                   cwipi_static_mesh, &
                                   cwipi_solver_cell_vertex, &
                                   1, &
                                   "Ensight Gold",&
                                   "text")




!
! Construction du maillage

  nvertex      = nvertexm
  nelts        = neltsm
  lconnecindex = lconnecindexm
  order = 1

  call creemaillagepolygone2d_f (order, &
                                 localcom, &
                                 xmin, &
                                 xmax, &
                                 ymin, &
                                 ymax, &
                                 initrandom, &
                                 nx, &
                                 ny, &
                                 nvertex, &
                                 coords, &
                                 nelts, &
                                 lconnecindex, &
                                 connecindex, &
                                 connec)

  call cwipi_define_mesh_f("test2D_4", &
                               nvertex, &
                               nelts, &
                               coords, &
                               connecindex, &
                               connec)

!
! Envoi de la coory a codec
! Reception de la coorx provenant de codec

  do i = 1, nvertex
     values(i) = coords(3*(i-1) + 2)
  enddo

  stride = 1
  call cwipi_exchange_f ("test2D_4", &
                             "echange1", &
                             stride, &
                             1, &
                             0.1d0, &
                             "cooy", &
                             values, &
                             "coox", &
                             localvalues, &
                             interpolationbidon_f, &
                             nNotLocatedPoints, &
                             status)
  call printStatus(iiunit, status)

!
! Suppression de l'objet couplage "couplingcellvertex"

  call cwipi_delete_coupling_f("test2D_4")
  write(iiunit, *)
  write(iiunit, *) "--------------------------------------------------------"
  write(iiunit, *)

!
! -------------------------------------
! test de la transmission d'un vecteur
! -------------------------------------
!
!
! Construction de "l'objet" couplage
  write(iiunit, *)
  write(iiunit, *) "--------------------------------------------------------"
  write(iiunit, *)
  write(iiunit, *) " Test 5"
  write(iiunit, *)

  call cwipi_create_coupling_f("test2D_5", &
                                   cwipi_cpl_parallel_with_part,&
                                   "CodeC", &
                                   2,     & ! Dimension des entites geometriques
                                   0.1d0, & ! Tolerance geometrique
                                   cwipi_static_mesh, &
                                   cwipi_solver_cell_vertex, &
                                   1, &
                                   "Ensight Gold",&
                                   "text")

!
! Construction du maillage

  nvertex      = nvertexm
  nelts        = neltsm
  lconnecindex = lconnecindexm
  order = 1

  call creemaillagepolygone2d_f (order, &
                                 localcom, &
                                 xmin, &
                                 xmax, &
                                 ymin, &
                                 ymax, &
                                 initrandom, &
                                 nx, &
                                 ny, &
                                 nvertex, &
                                 coords, &
                                 nelts, &
                                 lconnecindex, &
                                 connecindex, &
                                 connec)

  call cwipi_define_mesh_f("test2D_5", &
                               nvertex, &
                               nelts, &
                               coords, &
                               connecindex, &
                               connec)

!
! Envoi de la coory a codec
! Reception de la coorx provenant de codec

  stride = 6

!
! Finir

  do i = 1, nVertex
     do k = 1, stride
        if (k < 3) then
          values(stride*(i-1)+k) = coords(3*(i-1)+k);
        else if (k < 6) then
          values(stride*(i-1)+k) = coords(3*(i-1)+k-3);
       endif
     enddo
  enddo

  call cwipi_exchange_f ("test2D_5", &
                             "echange1", &
                             stride, &
                             1, &
                             0.1d0, &
                             "cooy", &
                             coords, &
                             "coox", &
                             localvalues, &
                             nNotLocatedPoints, &
                             status)
  call printStatus(iiunit, status)

!
! Suppression de l'objet couplage "couplingcellvertex"

  call cwipi_delete_coupling_f("test2D_5")
  write(iiunit, *)
  write(iiunit, *) "--------------------------------------------------------"
  write(iiunit, *)

!
! -------------------------------------
! test des sorties d'erreur
! -------------------------------------
!

!
! Construction de "l'objet" couplage
  write(iiunit, *)
  write(iiunit, *) "--------------------------------------------------------"
  write(iiunit, *)
  write(iiunit, *) " Test 6"
  write(iiunit, *)

  call cwipi_create_coupling_f("test2D_6", &
                                   cwipi_cpl_parallel_with_part,&
                                   "CodeC", &
                                   2,     & ! Dimension des entites geometriques
                                   0.1d0, & ! Tolerance geometrique
                                   cwipi_static_mesh, &
                                   cwipi_solver_cell_vertex, &
                                   1, &
                                   "Ensight Gold",&
                                   "text")

!
! Construction du maillage

  nvertex      = nvertexm
  nelts        = neltsm
  lconnecindex = lconnecindexm
  order = 1

  call creemaillagepolygone2d_f (order, &
                                 localcom, &
                                 xmin, &
                                 xmax, &
                                 ymin, &
                                 ymax, &
                                 initrandom, &
                                 nx, &
                                 ny, &
                                 nvertex, &
                                 coords, &
                                 nelts, &
                                 lconnecindex, &
                                 connecindex, &
                                 connec)

  call cwipi_define_mesh_f("test2D_6", &
                               nvertex, &
                               nelts, &
                               coords, &
                               connecindex, &
                               connec)

!
! Envoi de la coory a codec
! Reception de la coorx provenant de codec

  do i = 1, nvertex
     values(i) = coords(3*(i-1) + 2)
  enddo

  stride = 1
  call cwipi_exchange_f ("test2D_6", &
                             "echange1", &
                             stride, &
                             1, &
                             0.1d0, &
                             "cooy", &
                             values, &
                             "coox", &
                             localvalues, &
                             nNotLocatedPoints, &
                             status)
  call printStatus(iiunit, status)

  !
  ! Reception mais pas d'envoi alors que code_C attend quelque chose
  !

  stride = 1

  call cwipi_receive_f ("test2D_6", &
                        "echange2", &
                         stride, &
                         1, &
                         0.1d0, &
                         "cooYY", &
                         values, &
                         nNotLocatedPoints, &
                         status)
  call printStatus(iiunit, status)

!
! Suppression de l'objet couplage "couplingcellvertex"

  call cwipi_delete_coupling_f("test2D_6")
  write(iiunit, *)
  write(iiunit, *) "--------------------------------------------------------"
  write(iiunit, *)

!
! -------------------------------------
! Test simple localisation
! -------------------------------------
!

  write(iiunit, *)
  write(iiunit, *) "--------------------------------------------------------"
  write(iiunit, *)
  write(iiunit, *) " Test 7"
  write(iiunit, *)

  call cwipi_create_coupling_f("test2D_7", &
                                   cwipi_cpl_parallel_with_part,&
                                   "CodeC", &
                                   2,     & ! Dimension des entites geometriques
                                   0.1d0, & ! Tolerance geometrique
                                   cwipi_static_mesh, &
                                   cwipi_solver_cell_vertex, &
                                   1, &
                                   "Ensight Gold",&
                                   "text")

!
! Construction du maillage

  nvertex      = nvertexm
  nelts        = neltsm
  lconnecindex = lconnecindexm
  order = 1

  call creemaillagepolygone2d_f (order, &
                                 localcom, &
                                 xmin, &
                                 xmax, &
                                 ymin, &
                                 ymax, &
                                 initrandom, &
                                 nx, &
                                 ny, &
                                 nvertex, &
                                 coords, &
                                 nelts, &
                                 lconnecindex, &
                                 connecindex, &
                                 connec)

  call cwipi_define_mesh_f("test2D_7", &
                               nvertex, &
                               nelts, &
                               coords, &
                               connecindex, &
                               connec)

  call cwipi_locate_f("test2D_7")

  call cwipi_get_n_located_dist_pts_f("test2D_7", nDistantPoints)

  allocate(location(nDistantPoints))
  allocate(baryCooIdx(nDistantPoints+1))

  call cwipi_get_location_f("test2D_7", location)

  write(iiunit,*) "location",(location(i),i=1,nDistantPoints)

  call cwipi_get_bary_coord_idx_f("test2D_7", baryCooIdx)

  allocate(baryCoo(baryCooIdx(nDistantPoints+1)))

  call cwipi_get_bary_coord_f("test2D_7", baryCoo)

  allocate(tmpDbl(nDistantPoints))

  call cwipi_send_cellcenfd_eltcont_f("test2D_7" ,tmpDbl ,1)
  !call cwipi_send_cellvtxfd_eltcon_f("test2D_7" ,tmpDbl ,1)

  deallocate(location)
  deallocate(baryCooIdx)
  deallocate(baryCoo)
  deallocate(tmpDbl)


!
! Suppression de l'objet couplage "couplingcellvertex"

  call cwipi_delete_coupling_f("test2D_7")
  write(iiunit, *)
  write(iiunit, *) "--------------------------------------------------------"
  write(iiunit, *)

!
! ------------------------
! Test couplage p1 <-> p1
! ------------------------
!
! Construction de "l'objet" couplage

  write(iiunit, *)
  write(iiunit, *) "--------------------------------------------------------"
  write(iiunit, *)
  write(iiunit, *) " Test 8"
  write(iiunit, *)

  call cwipi_create_coupling_f("test2D_8", &
                                   cwipi_cpl_parallel_with_part,&
                                   "CodeC", &
                                   2,     & ! Dimension des entites geometriques
                                   0.1d0, & ! Tolerance geometrique
                                   cwipi_static_mesh, &
                                   cwipi_solver_cell_vertex, &
                                   1, &
                                   "Ensight Gold",&
                                   "text")


!
! Construction du maillage

  nvertex      = nvertexm
  nelts        = neltsm
  lconnecindex = lconnecindexm
  order = 1

  call creemaillagepolygone2d_f (order, &
                                 localcom, &
                                 xmin, &
                                 xmax, &
                                 ymin, &
                                 ymax, &
                                 initrandom, &
                                 nx, &
                                 ny, &
                                 nvertex, &
                                 coords, &
                                 nelts, &
                                 lconnecindex, &
                                 connecindex, &
                                 connec)

  call cwipi_define_mesh_f("test2D_8", &
                               nvertex, &
                               nelts, &
                               coords, &
                               connecindex, &
                               connec)



!
! Envoi de la coory a codec
! Reception de la coorx provenant de codec

  do i = 1, nvertex
     values(i) = coords(3*(i-1) + 2)
  enddo

  stride = 1
  call cwipi_exchange_f ("test2D_8", &
                             "echange1", &
                             stride, &
                             1, &
                             0.1d0, &
                             "cooy", &
                             values, &
                             "coox", &
                             localvalues, &
                             nNotLocatedPoints, &
                             status)

 
  call printStatus(iiunit, status)

!
! Suppression de l'objet couplage "couplingcellvertex"

  call cwipi_delete_coupling_f("test2D_8")
  write(iiunit, *)
  write(iiunit, *) "--------------------------------------------------------"
  write(iiunit, *)

!
! ------------------------
! Test couplage p1 <-> p1
! ------------------------
!
! Construction de "l'objet" couplage

  write(iiunit, *)
  write(iiunit, *) "--------------------------------------------------------"
  write(iiunit, *)
  write(iiunit, *) " Test 9"
  write(iiunit, *)

  call cwipi_create_coupling_f("test2D_9", &
                                   cwipi_cpl_parallel_without_part,&
                                   "CodeC", &
                                   2,     & ! Dimension des entites geometriques
                                   0.1d0, & ! Tolerance geometrique
                                   cwipi_static_mesh, &
                                   cwipi_solver_cell_vertex, &
                                   1, &
                                   "Ensight Gold",&
                                   "text")

!
! Construction du maillage

  nvertex      = nvertexm
  nelts        = neltsm
  lconnecindex = lconnecindexm
  order = 1

  call mpi_comm_group(localcom, localGroup, code)

  rl(1) = 0
  call mpi_group_incl(localGroup, 1, rl, p1Group, code)
  call mpi_comm_create(localcom, p1Group, p1Comm, code)

  if (currentRank == 0) then

     call creemaillagepolygone2d_f (order, &
          p1Comm, &
          xmin, &
          xmax, &
          ymin, &
          ymax, &
          initrandom, &
          nx, &
          ny, &
          nvertex, &
          coords, &
          nelts, &
          lconnecindex, &
          connecindex, &
          connec)

     call cwipi_define_mesh_f("test2D_9", &
                                  nvertex, &
                                  nelts, &
                                  coords, &
                                  connecindex, &
                                  connec)
  endif

!
! Envoi de la coory a codec
! Reception de la coorx provenant de codec
  stride = 1

  if (currentRank == 0) then
     do i = 1, nvertex
        values(i) = coords(3*(i-1) + 2)
     enddo

     call cwipi_exchange_f ("test2D_9", &
                                "echange1", &
                                stride, &
                                1, &
                                0.1d0, &
                                "cooy", &
                                values, &
                                "coox", &
                                localvalues, &
                                nNotLocatedPoints, &
                                status)
  else

     call cwipi_receive_f ("test2D_9", &
                               "echange1", &
                               stride, &
                               1, &
                               0.1d0, &
                               "coox", &
                               values, &
                               nNotLocatedPoints, &
                               status)

  endif

  call printStatus(iiunit, status)

!
! Suppression de l'objet couplage "couplingcellvertex"

  call cwipi_delete_coupling_f("test2D_9")
  write(iiunit, *)
  write(iiunit, *) "--------------------------------------------------------"
  write(iiunit, *)

!
! ------------------------
! Test couplage p1 <-> p1
! ------------------------
!
! Construction de "l'objet" couplage

  write(iiunit, *)
  write(iiunit, *) "--------------------------------------------------------"
  write(iiunit, *)
  write(iiunit, *) " Test 10"
  write(iiunit, *)

  call cwipi_create_coupling_f("test2D_10", &
                                   cwipi_cpl_parallel_without_part,&
                                   "CodeC", &
                                   2,     & ! Dimension des entites geometriques
                                   0.1d0, & ! Tolerance geometrique
                                   cwipi_static_mesh, &
                                   cwipi_solver_cell_vertex, &
                                   1, &
                                   "Ensight Gold",&
                                   "text")

!
! Construction du maillage

  nvertex      = nvertexm
  nelts        = neltsm
  lconnecindex = lconnecindexm
  order = 1

  call mpi_comm_group(localcom, localGroup, code)

  rl(1) = 0
  call mpi_group_incl(localGroup, 1, rl, p1Group, code)
  call mpi_comm_create(localcom, p1Group, p1Comm, code)

  if (currentRank == 0) then

     call creemaillagepolygone2d_f (order, &
          p1Comm, &
          xmin, &
          xmax, &
          ymin, &
          ymax, &
          initrandom, &
          nx, &
          ny, &
          nvertex, &
          coords, &
          nelts, &
          lconnecindex, &
          connecindex, &
          connec)

     call cwipi_define_mesh_f("test2D_10", &
                                  nvertex, &
                                  nelts, &
                                  coords, &
                                  connecindex, &
                                  connec)
  endif

!
! Envoi de la coory a codec
! Reception de la coorx provenant de codec
  stride = 1

  if (currentRank == 0) then
     do i = 1, nvertex
        values(i) = coords(3*(i-1) + 2)
     enddo

     call cwipi_exchange_f ("test2D_10", &
                                "echange1", &
                                stride, &
                                1, &
                                0.1d0, &
                                "cooy", &
                                values, &
                                "coox", &
                                localvalues, &
                                nNotLocatedPoints, &
                                status)
  else

     call cwipi_receive_f ("test2D_10", &
                               "echange1", &
                               stride, &
                               1, &
                               0.1d0, &
                               "coox", &
                               values, &
                               nNotLocatedPoints, &
                               status)

  endif

  call printStatus(iiunit, status)

!
! Suppression de l'objet couplage "couplingcellvertex"

  call cwipi_delete_coupling_f("test2D_10")
  write(iiunit, *)
  write(iiunit, *) "--------------------------------------------------------"
  write(iiunit, *)

  call cwipi_finalize_f()

  call mpi_finalize(code)

end program testf
