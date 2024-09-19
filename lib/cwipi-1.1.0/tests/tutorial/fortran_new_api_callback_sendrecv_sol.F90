#include "cwipi_configf.h"

  subroutine my_interpolation(c_local_code_name,        &
                              c_cpl_id,                 &
                              c_field_id,               &
                              i_part,                   &
                              c_buffer_in,              &
                              c_buffer_out)             &
    bind(c)
    use, intrinsic :: iso_c_binding

    use cwp

    implicit none

    character(kind=c_char,len=1)   :: c_local_code_name(*)
    character(kind=c_char,len=1)   :: c_cpl_id(*)
    character(kind=c_char,len=1)   :: c_field_id(*)
    integer(kind=c_int), value     :: i_part
    type(c_ptr), value             :: c_buffer_in
    type(c_ptr), value             :: c_buffer_out

    integer(kind = c_int)          :: n_components
    integer(kind = c_int)          :: n_elt_src
    integer(kind = c_int), pointer :: src_to_tgt_idx(:) => null()
    real(kind = c_double), pointer :: buffer_in(:)      => null()
    real(kind = c_double), pointer :: buffer_out(:)     => null()
    integer                        :: i, j, k

    character(len=:),      pointer :: local_code_name => null()
    character(len=:),      pointer :: cpl_id          => null()
    character(len=:),      pointer :: field_id        => null()

    print *, ">> my_interpolation"

    local_code_name => CWP_C_to_f_string(c_local_code_name)
    cpl_id          => CWP_C_to_f_string(c_cpl_id)
    field_id        => CWP_C_to_f_string(c_field_id)


    n_components = CWP_Field_n_components_get(local_code_name, &
                                              cpl_id,          &
                                              field_id)

    call CWP_Field_src_data_properties_get(local_code_name, &
                                           cpl_id,          &
                                           field_id,        &
                                           i_part,          &
                                           n_elt_src,       &
                                           src_to_tgt_idx)

    call c_f_pointer(c_buffer_in,  buffer_in,  [n_elt_src])
    call c_f_pointer(c_buffer_out, buffer_out, [src_to_tgt_idx(n_elt_src+1)])

    do i = 1, n_elt_src
      print *, "i = ", i
      do j = src_to_tgt_idx(i)+1, src_to_tgt_idx(i+1)
        print *, "  j = ", j
        do k = 1, n_components
          buffer_out(n_components*(j-1)+k) = buffer_in(n_components*(i-1)+k)
        enddo
      enddo
    enddo

    deallocate(local_code_name)
    deallocate(cpl_id)
    deallocate(field_id)
  end subroutine my_interpolation


program fortran_new_api_callback_sendrecv_sol

#ifdef CWP_HAVE_FORTRAN_MPI_MODULE
  use mpi
#endif

  use cwp
  use pdm_generate_mesh

  implicit none

#ifndef CWP_HAVE_FORTRAN_MPI_MODULE
    include "mpif.h"
#endif


  !--------------------------------------------------------------------
  integer                                     :: ierr
  integer                                     :: i_rank, n_rank
  integer(c_long), parameter                  :: n_vtx_seg = 4
  integer(c_int),  parameter                  :: n_part    = 1

  ! Coupling
  integer                                     :: n_code
  character(len = 99),                pointer :: code_names(:)         => null()
  integer                                     :: is_active_rank = CWP_STATUS_ON
  integer,                            pointer :: intra_comms(:)        => null()
  character(len = 99),                pointer :: coupled_code_names(:) => null()
  character(len = 99)                         :: coupling_name

  ! Mesh
  integer(c_int)                              :: n_vtx = 0
  double precision,                   pointer :: vtx_coord(:,:) => null()
  integer(c_long),                    pointer :: vtx_g_num(:)   => null()
  integer(c_int)                              :: n_elt = 0
  integer(c_int),                     pointer :: elt_vtx_idx(:) => null()
  integer(c_int),                     pointer :: elt_vtx(:)     => null()
  integer(c_long),                    pointer :: elt_g_num(:)   => null()
  integer(c_int)                              :: id_block

  ! Fields
  character(len = 99)                         :: field_name
  double precision,                   pointer :: send_field_data(:) => null()
  double precision,                   pointer :: recv_field_data(:) => null()

  integer                                     :: i, j
  character :: strnum
  !--------------------------------------------------------------------

  interface

    subroutine my_interpolation(c_local_code_name, &
                                c_cpl_id,          &
                                c_field_id,        &
                                i_part,            &
                                c_buffer_in,       &
                                c_buffer_out) bind(c)

      use, intrinsic :: iso_c_binding

      implicit none

      character(kind=c_char,len=1) :: c_local_code_name(*)
      character(kind=c_char,len=1) :: c_cpl_id(*)
      character(kind=c_char,len=1) :: c_field_id(*)
      integer(kind=c_int), value   :: i_part
      type(c_ptr), value           :: c_buffer_in
      type(c_ptr), value           :: c_buffer_out

    end subroutine
  end interface



  ! Initialize MPI
  call MPI_Init(ierr)
  call MPI_Comm_rank(mpi_comm_world, i_rank, ierr)
  call MPI_Comm_size(mpi_comm_world, n_rank, ierr)

  ! Initialize CWIPI
  n_code = 1

  allocate(code_names (n_code), &
           intra_comms(n_code))

  code_names(1) = "codeFortran"

  call CWP_Init(mpi_comm_world, &
                n_code,         &
                code_names,     &
                is_active_rank, &
                intra_comms)

  ! Create the coupling
  coupling_name = "coupling_Python_Fortran"
  allocate(coupled_code_names(n_code))
  coupled_code_names(1) = "codePython"

  call CWP_Cpl_create(code_names(1),                                         &
                      coupling_name,                                         &
                      coupled_code_names(1),                                 &
                      CWP_INTERFACE_SURFACE,                                 &
                      CWP_COMM_PAR_WITH_PART,                                &
                      CWP_SPATIAL_INTERP_FROM_LOCATION_MESH_LOCATION_OCTREE, &
                      n_part,                                                &
                      CWP_DYNAMIC_MESH_STATIC,                               &
                      CWP_TIME_EXCH_USER_CONTROLLED)

  ! Set coupling visualisation
  call CWP_Visu_set(code_names(1),           &
                    coupling_name,           &
                    1,                       &
                    CWP_VISU_FORMAT_ENSIGHT, &
                    "text")

  ! Generate mesh
  call PDM_generate_mesh_rectangle_simplified(intra_comms(1), &
                                              n_vtx_seg,      &
                                              n_vtx,          &
                                              n_elt,          &
                                              vtx_coord,      &
                                              elt_vtx_idx,    &
                                              elt_vtx)

  ! Set mesh vertices
  call CWP_Mesh_interf_vtx_set(code_names(1), &
                               coupling_name, &
                               0,             &
                               n_vtx,         &
                               vtx_coord,     &
                               vtx_g_num)

  ! Set mesh elements
  id_block = CWP_Mesh_interf_block_add(code_names(1),       &
                                       coupling_name,       &
                                       CWP_BLOCK_FACE_POLY)

  call CWP_Mesh_interf_f_poly_block_set(code_names(1), &
                                        coupling_name, &
                                        0,             &
                                        id_block,      &
                                        n_elt,         &
                                        elt_vtx_idx,   &
                                        elt_vtx,       &
                                        elt_g_num)

  ! Finalize mesh
  call CWP_Mesh_interf_finalize(code_names(1), &
                                coupling_name)


  ! Define field
  field_name = "coord_x"

  call CWP_Field_create(code_names(1),                &
                        coupling_name,                &
                        field_name,                   &
                        CWP_DOUBLE,                   &
                        CWP_FIELD_STORAGE_INTERLACED, &
                        1,                            &
                        CWP_DOF_LOCATION_CELL_CENTER, &
                        CWP_FIELD_EXCH_SENDRECV,      &
                        CWP_STATUS_ON)

  ! Begin time step :
  call CWP_Time_step_beg(code_names(1), &
                         0.d0)

  allocate(send_field_data(n_elt), &
           recv_field_data(n_elt))

  do i = 1, n_elt
    send_field_data(i) = 0.d0
    do j = elt_vtx_idx(i)+1, elt_vtx_idx(i+1)
      send_field_data(i) = send_field_data(i) + vtx_coord(1,elt_vtx(j))
    enddo
    send_field_data(i) = send_field_data(i) / (elt_vtx_idx(i+1) - elt_vtx_idx(i))
    ! print *, "elt ", i, " : ", send_field_data(i)
  enddo


  call CWP_Field_data_set(code_names(1),        &
                          coupling_name,        &
                          field_name,           &
                          0,                    &
                          CWP_FIELD_MAP_SOURCE, &
                          send_field_data)

  call CWP_Field_data_set(code_names(1),        &
                          coupling_name,        &
                          field_name,           &
                          0,                    &
                          CWP_FIELD_MAP_TARGET, &
                          recv_field_data)

  ! Set user-defined interpolation function
  call CWP_Field_interp_function_set(code_names(1), &
                                     coupling_name, &
                                     field_name,    &
                                     my_interpolation)


  ! Spatial interpolation
  call CWP_Spatial_interp_property_set(code_names(1), &
                                       coupling_name, &
                                       "tolerance",   &
                                       CWP_DOUBLE,    &
                                       "0.1")


  call CWP_Spatial_interp_weights_compute(code_names(1), &
                                          coupling_name)



  ! Exchange interpolated fields
  call CWP_Field_irecv(code_names(1), &
                       coupling_name, &
                       field_name)
  call CWP_Field_issend(code_names(1), &
                        coupling_name, &
                        field_name)

  call CWP_Field_wait_irecv(code_names(1), &
                            coupling_name, &
                            field_name)
  call CWP_Field_wait_issend(code_names(1), &
                             coupling_name, &
                             field_name)

  write (strnum, '(i1)') i_rank
  if (.false.) then
    call visu("check_Fortran_"//strnum//".vtk", &
              n_vtx,                            &
              vtx_coord,                        &
              n_elt,                            &
              elt_vtx_idx,                      &
              elt_vtx,                          &
              send_field_data,                  &
              recv_field_data)
  endif

  ! End time step :
  call CWP_Time_step_end(code_names(1))

  ! Delete field
  call CWP_Field_Del(code_names(1), &
                     coupling_name, &
                     field_name)

  ! Delete Mesh
  call CWP_Mesh_interf_del(code_names(1), &
                           coupling_name)

  ! Delete the coupling
  call CWP_Cpl_Del(code_names(1), &
                   coupling_name)


  ! Free memory
  deallocate(code_names,         &
             intra_comms,        &
             coupled_code_names, &
             send_field_data,    &
             recv_field_data)

  call pdm_fortran_free_c(c_loc(vtx_coord))
  call pdm_fortran_free_c(c_loc(elt_vtx_idx))
  call pdm_fortran_free_c(c_loc(elt_vtx))

  ! Finalize CWIPI
  call CWP_Finalize()

  print *, "Fortran ", i_rank, " FINISHED :D"

  ! Finalize MPI
  call MPI_Finalize(ierr)


contains

  subroutine visu(filename,    &
                  n_vtx,       &
                  vtx_coord,   &
                  n_elt,       &
                  elt_vtx_idx, &
                  elt_vtx,     &
                  send_field,  &
                  recv_field)
    implicit none
    character(*),  intent(in) :: filename
    integer,       intent(in) :: n_vtx, n_elt
    double precision, pointer :: vtx_coord(:,:), send_field(:), recv_field(:)
    integer(c_int),   pointer :: elt_vtx_idx(:), elt_vtx(:)

    integer                   :: fid = 13, i, j

    open(unit=fid, file=filename, action='write')
    write (fid,'(A26)') '# vtk DataFile Version 2.0'
    write (fid,*) filename
    write (fid,'(A5)') 'ASCII'
    write (fid,'(A16)') 'DATASET POLYDATA'

    write (fid,'(A6,1x,I0,1x,A6)') 'POINTS', n_vtx, 'double'
    do i = 1, n_vtx
      write (fid, *) vtx_coord(:,i)
    enddo

    write (fid,'(A8,1x,I0,1x,I0)') 'POLYGONS', n_elt, n_elt + elt_vtx_idx(n_elt+1)
    do i = 1, n_elt
      write (fid,'(I0,1x)', advance="no") elt_vtx_idx(i+1) - elt_vtx_idx(i)
      do j = elt_vtx_idx(i)+1,elt_vtx_idx(i+1)
        write (fid,'(I0,1x)', advance="no") elt_vtx(j) - 1
      enddo
      write (fid, *) ""
    enddo

    write (fid,'(A9,1x,I0)') 'CELL_DATA', n_elt
    write (fid,'(A17)') 'FIELD elt_field 2'
    write (fid,'(a12,1x,i0,1x,a6)') 'send_field 1', n_elt, 'double'
    do i = 1, n_elt
      write (fid, '(f20.16)') send_field(i)
    enddo
    write (fid,'(a12,1x,i0,1x,a6)') 'recv_field 1', n_elt, 'double'
    do i = 1, n_elt
      write (fid, '(f20.16)') recv_field(i)
    enddo

    close(fid)

  end subroutine visu

end program fortran_new_api_callback_sendrecv_sol
