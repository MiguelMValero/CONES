!-----------------------------------------------------------------------------
! This file is part of the ParaDiGM library.
!
! Copyright (C) 2023  ONERA
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

program pdm_t_mesh_partitioning_sol_f

  use pdm
#ifdef PDM_HAVE_FORTRAN_MPI_MODULE
  use mpi
#endif
  use pdm_multipart
  use pdm_vtk
  use pdm_dcube_nodal_gen
  use pdm_dmesh_nodal
  use pdm_part_mesh_nodal
  use pdm_mesh_nodal
  use pdm_part_extension
  use iso_c_binding
  use pdm_fortran
  use pdm_writer_wrapper
  use pdm_part_connectivity_transform

  implicit none

#ifndef PDM_HAVE_FORTRAN_MPI_MODULE
  include "mpif.h"
#endif

  !-----------------------------------------------------------
  integer (c_int)                       :: i
  logical                               :: fe   = .false.
  logical                               :: visu = .false.

  ! MPI
  integer                               :: code
  integer                               :: i_rank
  integer                               :: n_rank
  integer, parameter                    :: comm = MPI_COMM_WORLD

  ! PDM_dcube_nodal_gen_create
  integer                               :: n_x, n_y, n_z
  integer                               :: elt_type, order
  double precision                      :: length
  double precision                      :: xmin, ymin, zmin
  type (c_ptr)                          :: dcube

  ! PDM_dcube_nodal_gen_dmesh_nodal_get
  type (c_ptr)                          :: dmn

  ! PDM_multipart_create
  type (c_ptr)                          :: mpart
  integer (c_int)                       :: n_domain = 1
  integer(kind=PDM_l_num_s), pointer    :: n_part(:) => null()
  integer (c_int)                       :: part_method
  double precision,          pointer    :: part_fraction(:) => null()

  integer (c_int)                       :: i_domain    = 0
  integer (c_int)                       :: i_section = 0
  integer (c_int)                       :: i_part    = 0

  ! PDM_multipart_set_reordering_options
  integer(kind=PDM_l_num_s), pointer    :: renum_cell_properties(:) => null()

  ! PDM_multipart_get_part_mesh_nodal
  type (c_ptr)                          :: pmn

  ! PDM_part_mesh_nodal_section_n_elt_get
  integer (c_int)                       :: n_elt = -1

  ! PDM_part_mesh_nodal_section_std_get
  integer(kind=PDM_l_num_s), pointer  :: elt_vtx(:) => null()
  integer (pdm_g_num_s), pointer      :: elt_ln_to_gn(:) => null()
  integer(kind=PDM_l_num_s), pointer  :: parent_num(:) => null()
  integer (pdm_g_num_s), pointer      :: parent_entity_g_num(:) => null()

  ! PDM_multipart_part_vtx_coord_get
  double precision, pointer           :: coords(:,:) => null()
  integer(c_int)                      :: n_vtx = 0

  ! PDM_part_mesh_nodal_vtx_g_num_get
  integer (pdm_g_num_s), pointer      :: vtx_ln_to_gn(:) => null()

  ! FV
  integer(kind=PDM_g_num_s), pointer :: face_ln_to_gn(:) => null()
  integer(c_int)                     :: n_face = 0

  integer(kind=PDM_g_num_s), pointer :: cell_ln_to_gn(:) => null()
  integer(c_int)                     :: n_cell = 0
  integer(kind=PDM_l_num_s), pointer :: cell_face(:) => null()
  integer(kind=PDM_l_num_s), pointer :: cell_face_idx(:) => null()

  ! BONUS
  type(c_ptr)                     :: part_ext = C_NULL_PTR
  integer                         :: extend_type
  integer                         :: depth
  integer(PDM_l_num_s), pointer   :: face_vtx_idx(:)             => null()
  integer(PDM_l_num_s), pointer   :: face_vtx(:)                 => null()
  integer(PDM_l_num_s), pointer   :: vtx_part_bound_proc_idx(:)  => null()
  integer(PDM_l_num_s), pointer   :: vtx_part_bound_part_idx(:)  => null()
  integer(PDM_l_num_s), pointer   :: vtx_part_bound(:)           => null()

  integer                         :: n_cell_ext
  integer(pdm_l_num_s), pointer   :: cell_face_ext(:)      => null()
  integer(pdm_l_num_s), pointer   :: cell_face_ext_idx(:)  => null()
  integer(PDM_g_num_s), pointer   :: cell_ln_to_gn_ext(:)  => null()

  integer                         :: n_face_ext
  integer(pdm_l_num_s), pointer   :: face_vtx_ext(:)      => null()
  integer(pdm_l_num_s), pointer   :: face_vtx_ext_idx(:)  => null()
  integer(PDM_g_num_s), pointer   :: face_ln_to_gn_ext(:)  => null()

  integer                         :: n_vtx_ext
  double precision,     pointer   :: vtx_coord_ext(:,:)   => null()
  integer(PDM_g_num_s), pointer   :: vtx_ln_to_gn_ext(:)  => null()

  ! Visualization
  integer(pdm_l_num_s), pointer                :: pn_vtx(:)
  integer(pdm_l_num_s), pointer                :: pn_elt(:)
  integer(pdm_l_num_s), pointer                :: pn_face(:)

  type(PDM_pointer_array_t), pointer           :: pcoords => null()
  type(PDM_pointer_array_t), pointer           :: pvtx_ln_to_gn => null()
  type(PDM_pointer_array_t), pointer           :: pelt_vtx_idx => null()
  type(PDM_pointer_array_t), pointer           :: pelt_vtx => null()
  type(PDM_pointer_array_t), pointer           :: pelt_ln_to_gn => null()
  type(PDM_pointer_array_t), pointer           :: pcell_face_idx => null()
  type(PDM_pointer_array_t), pointer           :: pcell_face => null()
  !-----------------------------------------------------------

  call read_args(fe, visu)

  ! Initialize MPI environment
  call mpi_init(code)
  call mpi_comm_rank(comm, i_rank, code)
  call mpi_comm_size(comm, n_rank, code)

  ! Generate block-distributed parallelepided mesh
  n_x      = 10
  n_y      = 10
  n_z      = 10
  length   = 1.
  xmin     = 0.
  ymin     = 0.
  zmin     = 0.
  elt_type = PDM_MESH_NODAL_TETRA4
  order    = 1
  call PDM_dcube_nodal_gen_create(dcube,     &
                                  comm,      &
                                  n_x,       &
                                  n_y,       &
                                  n_z,       &
                                  length,    &
                                  xmin,      &
                                  ymin,      &
                                  zmin,      &
                                  elt_type,  &
                                  order,     &
                                  PDM_OWNERSHIP_USER)

  call PDM_dcube_nodal_gen_build(dcube, dmn)

  call PDM_dcube_nodal_gen_dmesh_nodal_get(dcube, dmn)

  call PDM_dmesh_nodal_generate_distribution(dmn)

  call PDM_dcube_nodal_gen_free(dcube)

  ! Create partitioning object
  allocate(n_part(n_domain))

  do i = 1, n_domain
    n_part(i) = 1
  end do

  allocate(part_fraction(n_part(i_domain+1)))

  do i = 1, n_part(i_domain+1)
    part_fraction(i) = 1
  end do

  part_method = PDM_SPLIT_DUAL_WITH_HILBERT
  call PDM_multipart_create(mpart,                     &
                            n_domain,                  &
                            n_part,                    &
                            PDM_FALSE,                 &
                            part_method,               &
                            PDM_PART_SIZE_HOMOGENEOUS, &
                            part_fraction,             &
                            comm,                      &
                            PDM_OWNERSHIP_KEEP)

  call PDM_multipart_set_reordering_options(mpart,                      &
                                            i_domain,                   &
                                            "PDM_PART_RENUM_CELL_NONE", &
                                            renum_cell_properties,      &
                                            "PDM_PART_RENUM_FACE_NONE")

  call PDM_multipart_dmesh_nodal_set(mpart,  &
                                     i_domain, &
                                     dmn)

  call PDM_multipart_compute(mpart)

  ! Get mesh arrrays in FE structure
  if (fe) then
    call PDM_multipart_get_part_mesh_nodal(mpart,    &
                                           i_domain, &
                                           pmn,      &
                                           PDM_OWNERSHIP_USER)

    call PDM_part_mesh_nodal_section_n_elt_get(pmn,       &
                                               i_section, &
                                               i_part,    &
                                               n_elt)

    call PDM_part_mesh_nodal_section_std_get(pmn,                 &
                                             i_section,           &
                                             i_part,              &
                                             elt_vtx,             &
                                             elt_ln_to_gn,        &
                                             parent_num,          &
                                             parent_entity_g_num, &
                                             PDM_OWNERSHIP_KEEP)

    call PDM_multipart_part_vtx_coord_get(mpart,              &
                                          i_domain,           &
                                          i_part,             &
                                          coords,             &
                                          PDM_OWNERSHIP_USER, &
                                          n_vtx)

    call PDM_part_mesh_nodal_vtx_g_num_get(pmn,    &
                                           i_part, &
                                           vtx_ln_to_gn)

    if (visu) then
      allocate(pn_vtx(1), &
               pn_elt(1))

      pn_vtx(1) = n_vtx
      pn_elt(1) = n_elt

      call PDM_pointer_array_create(pcoords,        1, PDM_TYPE_DOUBLE)
      call PDM_pointer_array_create(pvtx_ln_to_gn,  1, PDM_TYPE_G_NUM)
      call PDM_pointer_array_create(pelt_vtx,       1, PDM_TYPE_INT)
      call PDM_pointer_array_create(pelt_ln_to_gn,  1, PDM_TYPE_G_NUM)

      call PDM_pointer_array_part_set(pcoords,       0, coords)
      call PDM_pointer_array_part_set(pvtx_ln_to_gn, 0, vtx_ln_to_gn)
      call PDM_pointer_array_part_set(pelt_vtx,      0, elt_vtx)
      call PDM_pointer_array_part_set(pelt_ln_to_gn, 0, elt_ln_to_gn)

      call writer_wrapper(comm,          &
                          "visu",        &
                          "pmesh",       &
                          1,             &
                          pn_vtx,        &
                          pcoords,       &
                          pvtx_ln_to_gn, &
                          pn_elt,        &
                          pelt_vtx_idx,  &
                          pelt_vtx,      &
                          pelt_ln_to_gn, &
                          cell_t = elt_type)

      call PDM_pointer_array_free(pcoords)
      call PDM_pointer_array_free(pvtx_ln_to_gn)
      call PDM_pointer_array_free(pelt_vtx)
      call PDM_pointer_array_free(pelt_ln_to_gn)

      deallocate(pn_vtx, &
                 pn_elt)
    endif

    call PDM_part_mesh_nodal_free(pmn)
  end if

  ! Get mesh arrrays in FV structure
  if (.not.fe) then
    call PDM_multipart_part_ln_to_gn_get(mpart,               &
                                         i_domain,            &
                                         i_part,              &
                                         PDM_MESH_ENTITY_VTX, &
                                         vtx_ln_to_gn,        &
                                         PDM_OWNERSHIP_KEEP,  &
                                         n_vtx)

    call PDM_multipart_part_vtx_coord_get(mpart,             &
                                         i_domain,           &
                                         i_part,             &
                                         coords,             &
                                         PDM_OWNERSHIP_KEEP, &
                                         n_vtx)

    call PDM_multipart_part_ln_to_gn_get(mpart,                &
                                         i_domain,             &
                                         i_part,               &
                                         PDM_MESH_ENTITY_FACE, &
                                         face_ln_to_gn,        &
                                         PDM_OWNERSHIP_KEEP,   &
                                         n_face)

    call PDM_multipart_part_connectivity_get(mpart,                          &
                                             i_domain,                       &
                                             i_part,                         &
                                             PDM_CONNECTIVITY_TYPE_FACE_VTX, &
                                             face_vtx_idx,                   &
                                             face_vtx,                       &
                                             PDM_OWNERSHIP_KEEP,             &
                                             n_face)

    call PDM_multipart_part_ln_to_gn_get(mpart,                &
                                         i_domain,             &
                                         i_part,               &
                                         PDM_MESH_ENTITY_CELL, &
                                         cell_ln_to_gn,        &
                                         PDM_OWNERSHIP_KEEP,   &
                                         n_cell)

    call PDM_multipart_part_connectivity_get(mpart,                           &
                                             i_domain,                        &
                                             i_part,                          &
                                             PDM_CONNECTIVITY_TYPE_CELL_FACE, &
                                             cell_face_idx,                   &
                                             cell_face,                       &
                                             PDM_OWNERSHIP_KEEP,              &
                                             n_cell)


    if (visu) then
      allocate(pn_vtx(1), &
               pn_elt(1), &
               pn_face(1))

      pn_vtx(1)  = n_vtx
      pn_elt(1)  = n_cell
      pn_face(1) = n_face

      call PDM_pointer_array_create(pcoords,        1, PDM_TYPE_DOUBLE)
      call PDM_pointer_array_create(pvtx_ln_to_gn,  1, PDM_TYPE_G_NUM)
      call PDM_pointer_array_create(pelt_vtx_idx,   1, PDM_TYPE_INT)
      call PDM_pointer_array_create(pelt_vtx,       1, PDM_TYPE_INT)
      call PDM_pointer_array_create(pelt_ln_to_gn,  1, PDM_TYPE_G_NUM)
      call PDM_pointer_array_create(pcell_face_idx, 1, PDM_TYPE_INT)
      call PDM_pointer_array_create(pcell_face,     1, PDM_TYPE_INT)

      call PDM_pointer_array_part_set(pcoords,        0, coords)
      call PDM_pointer_array_part_set(pvtx_ln_to_gn,  0, vtx_ln_to_gn)
      call PDM_pointer_array_part_set(pelt_vtx_idx,   0, face_vtx_idx)
      call PDM_pointer_array_part_set(pelt_vtx,       0, face_vtx)
      call PDM_pointer_array_part_set(pelt_ln_to_gn,  0, cell_ln_to_gn)
      call PDM_pointer_array_part_set(pcell_face_idx, 0, cell_face_idx)
      call PDM_pointer_array_part_set(pcell_face,     0, cell_face)

      call writer_wrapper(comm,                            &
                          "visu",                          &
                          "pmesh",                         &
                          1,                               &
                          pn_vtx,                          &
                          pcoords,                         &
                          pvtx_ln_to_gn,                   &
                          pn_elt,                          &
                          pelt_vtx_idx,                    &
                          pelt_vtx,                        &
                          pelt_ln_to_gn,                   &
                          n_face         = pn_face,        &
                          pcell_face_idx = pcell_face_idx, &
                          pcell_face     = pcell_face)

      call PDM_pointer_array_free(pcoords)
      call PDM_pointer_array_free(pvtx_ln_to_gn)
      call PDM_pointer_array_free(pelt_vtx_idx)
      call PDM_pointer_array_free(pelt_vtx)
      call PDM_pointer_array_free(pelt_ln_to_gn)
      call PDM_pointer_array_free(pcell_face_idx)
      call PDM_pointer_array_free(pcell_face)

      deallocate(pn_vtx,  &
                 pn_elt,  &
                 pn_face)
    end if

  end if

  ! BONUS

  if (.not.fe) then
    ! step 1 : create
    extend_type = PDM_EXTEND_FROM_VTX
    depth       = 1
    call PDM_part_extension_create(part_ext,           &
                                   n_domain,           &
                                   n_part,             &
                                   extend_type,        & ! Extend from which element
                                   depth,              & ! Depth of the extension
                                   comm,               &
                                   PDM_OWNERSHIP_KEEP)

    ! step 2 : set
    call PDM_multipart_part_graph_comm_get(mpart,                   &
                                           i_domain,                &
                                           i_part,                  &
                                           PDM_MESH_ENTITY_VTX,     &
                                           vtx_part_bound_proc_idx, &
                                           vtx_part_bound_part_idx, &
                                           vtx_part_bound,          &
                                           PDM_OWNERSHIP_KEEP)

    call PDM_part_extension_connectivity_set(part_ext,                        &
                                             i_domain,                        &
                                             i_part,                          &
                                             PDM_CONNECTIVITY_TYPE_CELL_FACE, &
                                             cell_face_idx,                   &
                                             cell_face)

    call PDM_part_extension_connectivity_set(part_ext,                       &
                                             i_domain,                       &
                                             i_part,                         &
                                             PDM_CONNECTIVITY_TYPE_FACE_VTX, &
                                             face_vtx_idx,                   &
                                             face_vtx)

    call PDM_part_extension_vtx_coord_set(part_ext, &
                                          i_domain, &
                                          i_part,   &
                                          coords)

    call PDM_part_extension_ln_to_gn_set(part_ext,             &
                                         i_domain,             &
                                         i_part,               &
                                         PDM_MESH_ENTITY_CELL, &
                                         n_cell,               &
                                         cell_ln_to_gn)

    call PDM_part_extension_ln_to_gn_set(part_ext,             &
                                         i_domain,             &
                                         i_part,               &
                                         PDM_MESH_ENTITY_FACE, &
                                         n_face,               &
                                         face_ln_to_gn)

    call PDM_part_extension_ln_to_gn_set(part_ext,            &
                                         i_domain,            &
                                         i_part,              &
                                         PDM_MESH_ENTITY_VTX, &
                                         n_vtx,               &
                                         vtx_ln_to_gn)

    call PDM_part_extension_part_bound_graph_set(part_ext,                &
                                                 i_domain,                &
                                                 i_part,                  &
                                                 PDM_MESH_ENTITY_VTX,     &
                                                 vtx_part_bound_proc_idx, &
                                                 vtx_part_bound_part_idx, &
                                                 vtx_part_bound)

    ! step 3 : compute
    call PDM_part_extension_compute (part_ext)

    ! step 4 : get

    ! Cell
    call PDM_part_extension_ln_to_gn_get (part_ext,             &
                                          i_domain,             &
                                          i_part,               &
                                          PDM_MESH_ENTITY_CELL, &
                                          n_cell_ext,           &
                                          cell_ln_to_gn_ext)

    call PDM_part_extension_connectivity_get (part_ext,                        &
                                              i_domain,                        &
                                              i_part,                          &
                                              PDM_CONNECTIVITY_TYPE_CELL_FACE, &
                                              n_cell_ext,                      &
                                              cell_face_ext_idx,               &
                                              cell_face_ext)

    ! Face
    call PDM_part_extension_ln_to_gn_get (part_ext,             &
                                          i_domain,             &
                                          i_part,               &
                                          PDM_MESH_ENTITY_FACE, &
                                          n_face_ext,           &
                                          face_ln_to_gn_ext)

    call PDM_part_extension_connectivity_get (part_ext,                        &
                                              i_domain,                        &
                                              i_part,                          &
                                              PDM_CONNECTIVITY_TYPE_FACE_VTX, &
                                              n_face_ext,                      &
                                              face_vtx_ext_idx,               &
                                              face_vtx_ext)
    ! Vertices
    call PDM_part_extension_vtx_coord_get (part_ext,      &
                                           i_domain,      &
                                           i_part,        &
                                           n_vtx_ext,     &
                                           vtx_coord_ext)

    call PDM_part_extension_ln_to_gn_get (part_ext,            &
                                          i_domain,            &
                                          i_part,              &
                                          PDM_MESH_ENTITY_VTX, &
                                          n_vtx_ext,           &
                                          vtx_ln_to_gn_ext)

    ! step 5 : free
    call PDM_part_extension_free (part_ext)
  end if

  ! free
  deallocate(n_part, &
             part_fraction)
  call PDM_DMesh_nodal_free(dmn)
  call PDM_multipart_free(mpart)

  ! Finalize MPI environment
  call mpi_finalize(code)


contains

  subroutine usage(code)

    implicit none

    integer, intent(in) :: code

    write(*,*) "Usage :"
    write(*,*) " -fe     Use Finite-Element mesh (i.e. nodal connectivity)."
    write(*,*) " -visu   Output visualization files."
    write(*,*) " -h      This message."

    stop

  end subroutine usage


  subroutine read_args(fe, visu)

    implicit none

    logical,              intent(inout) :: fe
    logical,              intent(inout) :: visu
    integer                             :: argc, i, error
    character(999)                      :: arg

    argc = command_argument_count()

    i = 1
    do while (i <= argc)
      call get_command_argument(i, arg, status=error)
      if (error .ne. 0) then
        call usage(error)
      endif
      select case(arg)

        case ('-h')
          call usage(0)

        case ('-fe')
          fe = .true.

        case ('-visu')
          visu = .true.

        ! case ('-n_part')
        !   i = i + 1
        !   call get_command_argument(i, arg, status=error)
        !   if (error .ne. 0) then
        !     call usage(error)
        !   endif
        !   read(arg, *) n_part

      end select

      i = i + 1
    end do

  end subroutine read_args


end program pdm_t_mesh_partitioning_sol_f
