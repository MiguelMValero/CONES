#include "pdm_configf.h"

program testf

  use iso_c_binding
  use pdm
#ifdef PDM_HAVE_FORTRAN_MPI_MODULE
  use mpi
#endif
  use pdm_io
  use pdm_writer
  use pdm_sphere_surf_gen
  use pdm_fortran

  implicit none

#ifndef PDM_HAVE_FORTRAN_MPI_MODULE
  include "mpif.h"
#endif

  !-----------------------------------------------------------
  integer, parameter                :: comm        = MPI_COMM_WORLD
  character(len=99)                 :: arg
  integer                           :: i
  integer                           :: n_part      = 1
  integer                           :: part_method = 3

  integer(pdm_g_num_s)              :: n        = 5
  double precision                  :: x_center = 1.d0
  double precision                  :: y_center = 2.d0
  double precision                  :: z_center = 3.d0
  double precision                  :: radius   = 4.d0

  integer(pdm_l_num_s), pointer     :: pn_vtx(:)  => null()
  type(PDM_pointer_array_t), pointer:: pvtx_coord => null()
  type(PDM_pointer_array_t), pointer:: pvtx_ln_to_gn => null()
  integer(pdm_l_num_s), pointer     :: pn_face(:) => null()
  type(PDM_pointer_array_t), pointer:: pface_vtx_idx => null()
  type(PDM_pointer_array_t), pointer:: pface_vtx => null()
  type(PDM_pointer_array_t), pointer:: pface_ln_to_gn => null()

  ! double precision,     pointer     :: vtx_coord(:)     => null()
  double precision,     pointer     :: vtx_coord2(:,:)  => null()
  integer(pdm_g_num_s), pointer     :: vtx_ln_to_gn(:)  => null()
  integer(pdm_l_num_s), pointer     :: face_vtx(:)      => null()
  integer(pdm_g_num_s), pointer     :: face_ln_to_gn(:) => null()

  ! Writer
  type(c_ptr)                       :: wrt
  integer                           :: id_geom
  integer                           :: id_block

  ! MPI
  integer                           :: code
  integer                           :: i_rank
  integer                           :: n_rank

  integer                           :: ipart
  !-----------------------------------------------------------

  ! Read command line arguments
  i = 1
  do while (i <= command_argument_count())
    call get_command_argument(i, arg)
    select case(arg)
      case ("-n")
        i = i+1
        call get_command_argument(i, arg)
        read(arg, *) n
      case ("-n_part")
        i = i+1
        call get_command_argument(i, arg)
        read(arg, *) n_part
      case ("-parmetis")
        part_method = PDM_SPLIT_DUAL_WITH_PARMETIS
      case ("-pt-scotch")
        part_method = PDM_SPLIT_DUAL_WITH_PTSCOTCH
      case ("-hilbert")
        part_method = PDM_SPLIT_DUAL_WITH_HILBERT
      case default
        print *, "Invalid command argument ", arg
        stop
    end select

    i = i + 1
  enddo

  call mpi_init(code)
  call mpi_comm_rank(comm, i_rank, code)
  call mpi_comm_size(comm, n_rank, code)

  if (i_rank .eq. 0) then
    write(*, *) "-- Generate mesh"
  end if
  call PDM_sphere_surf_icosphere_gen_part(comm,           &
                                          n,              &
                                          x_center,       &
                                          y_center,       &
                                          z_center,       &
                                          radius,         &
                                          n_part,         &
                                          part_method,    &
                                          pn_vtx,         &
                                          pvtx_coord,     &
                                          pvtx_ln_to_gn,  &
                                          pn_face,        &
                                          pface_vtx_idx,  &
                                          pface_vtx,      &
                                          pface_ln_to_gn)


  if (i_rank .eq. 0) then
    write(*, *) "-- Write mesh"
  end if
  call PDM_writer_create(wrt,                    &
                         "Ensight",              &
                         PDM_WRITER_FMT_ASCII,   &
                         PDM_WRITER_TOPO_CST,    &
                         PDM_WRITER_OFF,         &
                         "sphere_surf",          &
                         "sphere_surf",          &
                         comm,                   &
                         PDM_IO_KIND_MPI_SIMPLE, &
                         1.d0,                   &
                         "")

  call PDM_writer_geom_create(wrt,     &
                              id_geom, &
                              "mesh",  &
                              n_part)

  call PDM_writer_step_beg(wrt, 0.d0)

  !  Write geometry
  do ipart = 1, n_part
    call PDM_pointer_array_part_get(pvtx_ln_to_gn, &
                                    ipart-1,       &
                                    vtx_ln_to_gn)
    ! call PDM_pointer_array_part_get(pvtx_coord, &
    !                                 ipart-1,    &
    !                                 vtx_coord)

    ! Reshape without copy
    ! call c_f_pointer(c_loc(vtx_coord), vtx_coord2, [3, pn_vtx(ipart)])
    ! call PDM_pointer_array_part_get(pvtx_coord, &
    !                                 ipart-1,    &
    !                                 vtx_coord)
    call PDM_pointer_array_part_get(pvtx_coord,                &
                                    ipart-1,                   &
                                    PDM_STRIDE_CST_INTERLACED, &
                                    3,                         &
                                    vtx_coord2)

    call PDM_writer_geom_coord_set(wrt,                &
                                   id_geom,            &
                                   ipart-1,            &
                                   pn_vtx(ipart),      &
                                   vtx_coord2,         &
                                   vtx_ln_to_gn,       &
                                   PDM_OWNERSHIP_USER)


    call PDM_writer_geom_bloc_add(wrt,                  &
                                  id_geom,              &
                                  2,                    & ! PDM_MESH_NODAL_TRIA3
                                  PDM_OWNERSHIP_USER,   &
                                  id_block)

    call PDM_pointer_array_part_get(pface_ln_to_gn, &
                                    ipart-1,        &
                                    face_ln_to_gn)
    call PDM_pointer_array_part_get(pface_vtx, &
                                    ipart-1,   &
                                    face_vtx)

    call PDM_writer_geom_bloc_std_set(wrt,            &
                                      id_geom,        &
                                      id_block,       &
                                      ipart-1,        &
                                      pn_face(ipart), &
                                      face_vtx,       &
                                      face_ln_to_gn)
  end do

  call PDM_writer_geom_write(wrt,     &
                             id_geom)

  call PDM_writer_step_end(wrt)

  call PDM_writer_free(wrt)

  !  Free memory

  call pdm_fortran_free_c(c_loc(pn_vtx))
  call pdm_fortran_free_c(c_loc(pn_face))
  call PDM_pointer_array_free(pvtx_coord)
  call PDM_pointer_array_free(pvtx_ln_to_gn)
  call PDM_pointer_array_free(pface_vtx_idx)
  call PDM_pointer_array_free(pface_vtx)
  call PDM_pointer_array_free(pface_ln_to_gn)


  if (i_rank .eq. 0) then
    write(*, *) "-- End"
  end if

  call mpi_finalize(code)

end program testf
