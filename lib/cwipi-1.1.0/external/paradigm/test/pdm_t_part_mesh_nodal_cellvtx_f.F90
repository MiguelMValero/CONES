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
  use pdm_fortran
  use pdm_part_mesh_nodal
  use iso_c_binding

  implicit none

#ifndef PDM_HAVE_FORTRAN_MPI_MODULE
  include "mpif.h"
#endif

  !-----------------------------------------------------------
  integer,                   parameter :: f_comm  = MPI_COMM_WORLD

  type(c_ptr)                          :: mesh = C_NULL_PTR
  double precision,          pointer   :: vtx_coord(:,:)      => null()
  integer(kind=pdm_g_num_s), pointer   :: vtx_ln_to_gn(:)     => null()
  integer(kind=pdm_g_num_s), pointer   :: cell_ln_to_gn(:)    => null()
  integer,                   pointer   :: cell_vtx_idx(:)     => null()
  integer,                   pointer   :: cell_vtx(:)         => null()
  integer,                   pointer   :: cell_face_idx(:)    => null()
  integer,                   pointer   :: cell_face_nb(:)     => null()
  integer,                   pointer   :: cell_face(:)        => null()
  integer(kind=pdm_g_num_s), pointer   :: face_ln_to_gn(:)    => null()
  integer,                   pointer   :: face_vtx_idx(:)     => null()
  integer,                   pointer   :: face_vtx_nb(:)      => null()
  integer,                   pointer   :: face_vtx(:)         => null()
  integer,                   pointer   :: tmp_cell_vtx_idx(:) => null()
  integer,                   pointer   :: tmp_cell_vtx(:)     => null()
  integer                              :: code
  integer                              :: i_rank
  integer                              :: n_rank
  integer                              :: n_vtx
  integer                              :: n_face
  integer                              :: n_cell

  integer                              :: i, ifac, isom, tmp(6)
  integer                              :: fid = 13
  character                            :: strnum
  !-----------------------------------------------------------

  ! Init
  call mpi_init(code)
  call mpi_comm_rank(f_comm, i_rank, code)
  call mpi_comm_size(f_comm, n_rank, code)

  if (n_rank .ne. 2) then
    print *,'Error : 2 MPI processes are mandatory'
    call mpi_finalize(code)
    stop
  end if

  ! Read mesh
  if (i_rank .eq. 0) then
    write(*, *) "-- Read mesh"
  end if

  write (strnum, '(i1)') i_rank+1
  fid = fid + i_rank
  open(unit=fid, &
file=&
PDM_MESH_DIR_F&
//"/test/meshes/mixed_elements_cellvtx."//strnum, &
action='read')

  do i = 1,26
    read(fid,*)
  end do

  read(fid,*) n_vtx
  read(fid,*) n_face
  read(fid,*) n_cell
  read(fid,*)

  read(fid,*)
  allocate(vtx_ln_to_gn(n_vtx))
  allocate(vtx_coord(3,n_vtx))
  do i = 1,n_vtx
    read(fid,*) vtx_ln_to_gn(i), vtx_coord(1:3,i)
  end do

  read(fid,*)
  allocate(face_ln_to_gn(n_face))
  allocate(face_vtx_idx(n_face+1))
  allocate(face_vtx_nb(n_face))
  allocate(face_vtx(4*n_face))
  face_vtx_idx(1) = 1
  do i = 1,n_face
    read(fid,*) face_ln_to_gn(i), face_vtx_nb(i), tmp(1:face_vtx_nb(i))
    face_vtx_idx(i+1) = face_vtx_idx(i) + face_vtx_nb(i)
    face_vtx(face_vtx_idx(i):face_vtx_idx(i+1)-1) = tmp(1:face_vtx_nb(i))
    do isom = 1,face_vtx_nb(i)
      face_vtx(face_vtx_idx(i)+isom-1) = minloc(abs(vtx_ln_to_gn-face_vtx(face_vtx_idx(i)+isom-1)),1)
    end do
  end do

  read(fid,*)
  allocate(cell_face_idx(n_cell+1))
  allocate(cell_face_nb(n_cell))
  allocate(cell_face(6*n_cell))
  allocate(cell_ln_to_gn(n_cell))
  cell_face_idx(1) = 1
  do i = 1,n_cell
    read(fid,*) cell_ln_to_gn(i), cell_face_nb(i), tmp(1:cell_face_nb(i))
    cell_face_idx(i+1) = cell_face_idx(i) + cell_face_nb(i)
    cell_face(cell_face_idx(i):cell_face_idx(i+1)-1) = tmp(1:cell_face_nb(i))
    do ifac = 1,cell_face_nb(i)
      cell_face(cell_face_idx(i)+ifac-1) = minloc(abs(face_ln_to_gn-cell_face(cell_face_idx(i)+ifac-1)),1)
    end do
  end do
  close(fid)

  ! Convert Fortran -> C
  face_vtx_idx = face_vtx_idx - 1
  cell_face_idx = cell_face_idx - 1

  if (i_rank .eq. 0) then
    write(*, *) "-- Create nodal mesh"
  end if

  ! Création de l'objet "MAILLAGE NODAL"
  call pdm_part_mesh_nodal_create               & !
  (mesh,                                        & !- IDENTIFICATEUR OBJET MAILLAGE NODAL
   3,                                           & !- DIMENSION DU MAILLAGE NODAL
   1,                                           & !- NOMBRE DE PARTITIONS DU MAILLAGE NODAL SUR LE PROCESSUS COURANT
   f_comm)                                        !- COMMUNICATEUR MPI

  if (i_rank .eq. 0) then
    write(*, *) "-- Set vertices"
  end if

  ! Définition des sommets du maillage nodal
  call pdm_part_mesh_nodal_coord_set            &
  (mesh,                                        & !- IDENTIFICATEUR OBJET MAILLAGE NODAL
   0,                                           & !- INDICE DE PARTITION DU MAILLAGE NODAL
   n_vtx,                                       & !- NOMBRE DE SOMMETS
   vtx_coord,                                   & !- COORDONNEES DES SOMMETS
   vtx_ln_to_gn,                                & !- NUMEROTATION ABSOLUE DES SOMMETS
   PDM_OWNERSHIP_USER)                            !- OWNERSHIP

  if (i_rank .eq. 0) then
    write(*, *) "-- Set cell-face connectivity"
  end if

  ! Définition de la connectivité du maillage 3D
  call pdm_part_mesh_nodal_cell3d_cellface_add  &
  (mesh,                                        & !- IDENTIFICATEUR OBJET MAILLAGE NODAL
   0,                                           & !- INDICE DE PARTITION DU MAILLAGE NODAL
   n_cell,                                      & !- NOMBRE DE CELLULES
   n_face,                                      & !- NOMBRE DE FACES
   face_vtx_idx,                                & !- ADRESSES DES NUMEROS DE FACES PAR CELLULE
   face_vtx,                                    & !- NUMEROS DE SOMMETS PAR FACE
   face_ln_to_gn,                               & !- NUMEROTATION ABSOLUE DES FACES
   cell_face_idx,                               & !- ADRESSES DES NUMEROS DE FACES PAR CELLULE
   cell_face,                                   & !- NUMEROS DE FACES PAR CELLULE
   cell_ln_to_gn,                               & !- NUMEROTATION ABSOLUE DES CELLULES
   PDM_OWNERSHIP_KEEP)                            !- OWNERSHIP

  if (i_rank .eq. 0) then
    write(*, *) "-- Get cell-vertex connectivity"
  end if

  ! Récupération de la connectivité du maillage
  call pdm_part_mesh_nodal_cell_vtx_connect_get & !-
  (mesh,                                        & !- IDENTIFICATEUR OBJET LOCALISATEUR
   0,                                           & !- INDICE DE PARTITION DU MAILLAGE NODAL
   tmp_cell_vtx_idx,                            & !- ADRESSES DES NUMEROS DE SOMMETS PAR CELLULE
   tmp_cell_vtx)                                  !- NUMEROS DE SOMMETS PAR CELLULE

  ! Conversion de la connectivité
  allocate(cell_vtx_idx(n_cell+1))
  allocate(cell_vtx(tmp_cell_vtx_idx(n_cell+1)))
  cell_vtx_idx = tmp_cell_vtx_idx
  cell_vtx = tmp_cell_vtx

  ! Free memory
  call pdm_part_mesh_nodal_free(mesh)
  call pdm_fortran_free_c(c_loc(tmp_cell_vtx_idx))
  call pdm_fortran_free_c(c_loc(tmp_cell_vtx))
  nullify(tmp_cell_vtx_idx)
  nullify(tmp_cell_vtx)

  if (i_rank .eq. 0) then
    write(*, *) "-- Create nodal mesh"
  end if

  ! Création de l'objet "MAILLAGE NODAL"
  call pdm_part_mesh_nodal_create               & !
  (mesh,                                        & !- IDENTIFICATEUR OBJET MAILLAGE NODAL
   3,                                           & !- DIMENSION DU MAILLAGE NODAL
   1,                                           & !- NOMBRE DE PARTITIONS DU MAILLAGE NODAL SUR LE PROCESSUS COURANT
   f_comm)                                        !- COMMUNICATEUR MPI

  if (i_rank .eq. 0) then
    write(*, *) "-- Set vertices"
  end if

  ! Définition des sommets du maillage nodal
  call pdm_part_mesh_nodal_coord_set            &
  (mesh,                                        & !- IDENTIFICATEUR OBJET MAILLAGE NODAL
   0,                                           & !- INDICE DE PARTITION DU MAILLAGE NODAL
   n_vtx,                                       & !- NOMBRE DE SOMMETS
   vtx_coord,                                   & !- COORDONNEES DES SOMMETS
   vtx_ln_to_gn,                                & !- NUMEROTATION ABSOLUE DES SOMMETS
   PDM_OWNERSHIP_USER)                            !- OWNERSHIP

  if (i_rank .eq. 0) then
    write(*, *) "-- Set cell-vertex connectivity"
  end if

  ! Définition de la connectivité du maillage nodal
  call pdm_part_mesh_nodal_cells_cellvtx_add    &
  (mesh,                                        & !- IDENTIFICATEUR OBJET MAILLAGE NODAL
   0,                                           & !- INDICE DE PARTITION DU MAILLAGE NODAL
   n_cell,                                      & !- NOMBRE DE CELLULES
   cell_vtx_idx,                                & !- ADRESSES DES NUMEROS DE SOMMETS PAR CELLULE
   cell_vtx,                                    & !- NUMEROS DE SOMMETS PAR CELLULE
   cell_ln_to_gn,                               & !- NUMEROTATION ABSOLUE DES CELLULES
   PDM_OWNERSHIP_KEEP)                            !- OWNERSHIP

  ! Free memory
  call pdm_part_mesh_nodal_free(mesh)
  deallocate(vtx_coord)
  deallocate(vtx_ln_to_gn)
  deallocate(face_ln_to_gn)
  deallocate(face_vtx_idx)
  deallocate(face_vtx_nb)
  deallocate(face_vtx)
  deallocate(cell_face_idx)
  deallocate(cell_face_nb)
  deallocate(cell_face)
  deallocate(cell_ln_to_gn)
  deallocate(cell_vtx_idx)
  deallocate(cell_vtx)

  if (i_rank .eq. 0) then
    write(*, *) "-- End"
  end if

  call mpi_finalize(code)

end program testf
