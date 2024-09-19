#include "cwipi_configf.h"

program fortran_new_api_polygon_ex

#ifdef CWP_HAVE_FORTRAN_MPI_MODULE
    use mpi
#endif
    use cwp

    implicit none

#ifndef CWP_HAVE_FORTRAN_MPI_MODULE
    include "mpif.h"
#endif

  !--------------------------------------------------------------------
  integer, parameter                      :: n_vtx = 11, n_elts = 5

  integer                                 :: ierr
  integer                                 :: i_rank, n_rank

  integer                                 :: n_code
  character(len = 5),            pointer  :: code_names(:)         => null()

  integer                                 :: n_part
  character(len = 5),            pointer  :: coupled_code_names(:) => null()
  character(len = 99)                     :: coupling_name

  double precision, dimension(:,:), pointer :: coords => null()
!  integer(c_long), pointer, dimension(:)  :: vtx_g_num => null()

  integer, pointer, dimension(:)          :: connec_idx => null()
  integer, pointer, dimension(:)          :: connec => null()
!  integer(c_long), pointer, dimension(:)  :: elt_g_num  => null()
!  integer(c_int)                          :: id_block

  character(len = 99)                     :: field_name
  integer(c_int)                          :: n_components

  integer                                 :: i

  double precision,              pointer  :: send_field_data(:) => null()
  double precision,              pointer  :: recv_field_data(:) => null()

!  integer(c_int)                          :: n_uncomputed_tgts
!  integer(c_int),                pointer  :: uncomputed_tgts(:) => null()

  logical                                 :: I_am_code1
  !--------------------------------------------------------------------

  ! MPI Initialization
  call MPI_Init(ierr)
  call MPI_Comm_rank(mpi_comm_world, i_rank, ierr)
  call MPI_Comm_size(mpi_comm_world, n_rank, ierr)

  ! This test mimics the coupling between 2 codes running each
  ! on one processor.
  if (n_rank /= 2) then
    print *, " This executable must be executed on 2 processes"
    stop
  endif

  ! Initialize CWIPI :
  ! Use CWP_Init for code1 running on MPI rank 0 and code2
  ! running on MPU rank 1.
  ! ------------------------------------------------------- To fill in
  n_code = 1

  allocate(code_names(n_code))

  ! for code1
  if (i_rank == 0) then
    code_names(1) = "code1"
    I_am_code1    = .true.
  ! for code2
  else
    code_names(1) = "code2"
    I_am_code1    = .false.
  endif

  ! ---------------------------------------------------- End To fill in

  ! Create the coupling :
  ! Use CWP_Cpl_create to couple code1 with code2 on a surface
  ! interface. Operate the localization with the octree method.
  ! In this tutorial, one coupling iteration is done on a mesh
  ! partitionned over all processors. This induces that the mesh
  ! is static.
  ! ------------------------------------------------------- To fill in
  coupling_name     = "code1_code2";

  allocate(coupled_code_names(n_code))

  ! for code1
  if (I_am_code1) then
    coupled_code_names(1) = "code2"
  ! for code2
  else
    coupled_code_names(1) = "code1"
  endif

  n_part = 1

  ! ---------------------------------------------------- End To fill in

  ! Set coupling visualisation:
  ! Use CWP_Visu_set to output ASCII Ensight format files
  ! at each iteration step.
  ! ------------------------------------------------------- To fill in

  ! ---------------------------------------------------- End To fill in

  ! Set the mesh vertices coordinates :
  ! Use CWP_Mesh_interf_vtx_set to set the mesh vertex coordinates,
  ! no global numbering of the vertices will be given. In this
  ! simple setting, there is only one partition per processor.
  ! ------------------------------------------------------- To fill in
  allocate(coords(3,n_vtx))
  coords(:, 1) = [0,0,0]
  coords(:, 2) = [1,0,0]
  coords(:, 3) = [2,0,0]
  coords(:, 4) = [3,0,0]
  coords(:, 5) = [0,1,0]
  coords(:, 6) = [2,1,0]
  coords(:, 7) = [3,1,0]
  coords(:, 8) = [1,2,0]
  coords(:, 9) = [0,3,0]
  coords(:,10) = [2,3,0]
  coords(:,11) = [3,3,0]

  ! ---------------------------------------------------- End To fill in

  ! Set the mesh polygons connectivity :
  ! Use CWP_Mesh_interf_block_add to create a block of
  ! of polygons. Choose the correct CWIPI function
  ! to set a polygonal mesh, no need to give the elements
  ! global numbering.
  ! ------------------------------------------------------- To fill in

  allocate(connec_idx(n_elts+1))
  connec_idx = [0,3,7,11,16,21]
  allocate(connec(21))
  connec = [1,2,5,   3,4,7,6,   5,8,10,9   ,5,2,3,6,8,   6,7,11,10,8]

  ! ---------------------------------------------------- End To fill in

  ! Finalize mesh :
  ! Use the correct CWIPI function to generate the
  ! mesh global numbering and the underlying mesh data structure.
  ! ------------------------------------------------------- To fill in

  ! ---------------------------------------------------- End To fill in

  ! Create and set the field values :
  ! Use CWP_Field_create and CWP_Field_data_set to create and set
  ! a field onto the mesh. code1 will send its field which code2
  ! will receive. The field is located at the mesh nodes.
  ! There is only one mesh partition in this tutorial. Activate
  ! visualization for this field if you wish it to be in the
  ! Ensight file.
  ! Do not forget to begin the time step AFTER creating the fields,
  ! but BEFORE setting the fields data!
  ! ------------------------------------------------------- To fill in
  field_name   = "a super fancy field"
  n_components = 1

  ! Create the field :
  ! ...

  ! Begin the time step :
  ! ...

  ! Set the field pointers :
  ! for code1
  if (I_am_code1) then

    allocate(send_field_data(n_vtx * n_components))
    do i=1,n_vtx
      send_field_data(i) = coords(1,i)
    end do

    ! ...

  ! for code2
  else

    allocate(recv_field_data(n_vtx * n_components))

    ! ...

  endif
  ! ---------------------------------------------------- End To fill in

  ! Compute interpolation weights :
  ! Choose the two correct CWIPI functions to set the geometric
  ! tolerance to 10% of an element size for point localisation
  ! and to compute the interpolation weights.
  ! ------------------------------------------------------- To fill in

  ! ---------------------------------------------------- End To fill in

  ! Exchange field values between codes :
  ! Use the CWIPI exchange functions similar to the MPI ones
  ! for code1 to send a field and code2 to receive that field.
  ! The reason for the double if statements is that the exchange
  ! operation is non-blocking. That means that work can be done
  ! before calling the wait function where the field need to be
  ! used or changed.
  ! ------------------------------------------------------- To fill in

  ! for code1
  if (I_am_code1) then

  ! for code2
  else

  endif

  ! for code1
  if (I_am_code1) then

  ! for code2
  else

  endif
  ! ---------------------------------------------------- End To fill in

  ! End the time step :
  ! ------------------------------------------------------- To fill in
  ! ---------------------------------------------------- End To fill in

  ! Check interpolation :
  ! The field that has been sent will be interpolated on the vertices
  ! of the mesh of the receiving code. Depending on the geometric
  ! tolerance and the interpolation algorithm used, some points might
  ! not be located. Users often want all points on the receiving mesh
  ! to have a field value. Thus, for the receiving code, check if there
  ! are vertices for which the interpolation has been unsuccessful.
  ! ------------------------------------------------------- To fill in

  ! ---------------------------------------------------- End To fill in

  ! Delete field :
  ! ------------------------------------------------------- To fill in

  ! ---------------------------------------------------- End To fill in

  ! Delete Mesh :
  ! ------------------------------------------------------- To fill in

  ! ---------------------------------------------------- End To fill in

  ! Delete the coupling :
  ! ------------------------------------------------------- To fill in

  ! ---------------------------------------------------- End To fill in

  ! free
  deallocate(coords);
  deallocate(connec);
  deallocate(connec_idx);
  if (I_am_code1) then
    deallocate(send_field_data);
  else
    deallocate(recv_field_data);
  endif

  ! Finalize CWIPI :
  ! ------------------------------------------------------- To fill in

  ! ---------------------------------------------------- End To fill in

  ! Finalize MPI :
  call MPI_Finalize(ierr)

end program fortran_new_api_polygon_ex
