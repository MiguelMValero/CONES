#include "cwipi_configf.h"

program fortran_new_api_deformable_sol

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
  integer(kind=8), parameter                  :: n_vtx_seg = 10

  integer                                     :: n_code
  character(len = 99),                pointer :: code_names(:)         => null()
!  integer                                     :: is_active_rank = CWP_STATUS_ON
  integer,                            pointer :: intra_comms(:)        => null()

  integer                                     :: n_part
  character(len = 99),                pointer :: coupled_code_names(:) => null()
  character(len = 99)                         :: coupling_name

  double precision,   dimension(:,:), pointer :: coords    => null()
!  integer(c_long),    dimension(:),   pointer :: vtx_g_num => null()

  double precision,   dimension(:,:), pointer :: xyz_dest  => null()
!  integer(c_long),    dimension(:),   pointer :: pts_g_num => null()

  integer,            dimension(:),   pointer :: elt_vtx_idx => null()
  integer,            dimension(:),   pointer :: elt_vtx     => null()
!  integer(c_long),    dimension(:),   pointer :: elt_g_num   => null()
!  integer(c_int)                              :: id_block

  character(len = 99)                         :: send_field_name
  character(len = 99)                         :: recv_field_name
  integer(c_int)                              :: n_components

  double precision                            :: ampl
  double precision                            :: dt
  double precision                            :: freq
  double precision                            :: omega
  double precision                            :: phi
  double precision                            :: time
  integer                                     :: i, j, it, itdeb, itend
  integer                                     :: n_vtx, n_elt

  double precision,                   pointer :: send_field_data(:) => null()
  double precision,                   pointer :: recv_field_data(:) => null()

  integer(c_int)                              :: n_uncomputed_tgts
  integer(c_int),                     pointer :: uncomputed_tgts(:) => null()
  !--------------------------------------------------------------------

  ! MPI Initialization :

  call MPI_Init(ierr)
  call MPI_Comm_rank(mpi_comm_world, i_rank, ierr)
  call MPI_Comm_size(mpi_comm_world, n_rank, ierr)

  ! Initialize CWIPI :
  ! Use CWP_Init for code2 written in C.
  ! ------------------------------------------------------- To fill in
  n_code = 1

  allocate(code_names(n_code))

  code_names(1) = "code1"

  ! ---------------------------------------------------- End To fill in

  ! Create the coupling :
  ! In this tutorial the mesh will be deformed over the coupling
  ! itreations. Those iterations mimic solver steps. Use CWP_Cpl_create
  ! to couple code1 with code2. A rectangular mesh is used with triangle
  ! elements. Operate the localization with the octree method.
  ! ------------------------------------------------------- To fill in
  coupling_name = "coupling"
  allocate(coupled_code_names(n_code))
  coupled_code_names(1) = "code2"
  n_part = 1;

  ! ---------------------------------------------------- End To fill in

  ! Set coupling visualisation:
  ! Use CWP_Visu_set to output ASCII Ensight format files
  ! at each iteration step. Open the CHR.case file in
  ! paraview to visualize the deformed mesh and exchanged fields.
  ! ------------------------------------------------------- To fill in

  ! ---------------------------------------------------- End To fill in

  ! Create mesh :
  ! A ParaDiGM library function is used to generate the rectangular mesh.
  ! PDM_MPI_mpi_2_pdm_mpi_comm is used to convert the MPI communicator into
  ! the expected ParaDiGM format. It is the users responsability to
  ! free arrays from _simplified mesh generation functions. In a real life
  ! coupling case here a user generated mesh would be read/loaded/given.
  ! TO UNCOMMENT -->>
  n_vtx = 0
  n_elt = 0
  ! call PDM_generate_mesh_rectangle_simplified(intra_comms(1), &
  !                                             n_vtx_seg,      &
  !                                             n_vtx,          &
  !                                             n_elt,          &
  !                                             coords,         &
  !                                             elt_vtx_idx,    &
  !                                             elt_vtx)
  ! <<--

  ! Interations :
  ! At each iteration the mesh coordinates and the exchanged fields are modified.
  itdeb = 1
  itend = 10
  freq  = 0.20d0
  ampl  = 0.012d0
  phi   = 0.1d0
  time  = 0.0d0
  dt    = 0.1d0

  omega = 2.0d0*acos(-1.0d0)*freq

  allocate(send_field_data(n_elt))
  allocate(recv_field_data(n_elt))
  allocate(xyz_dest(3, n_elt))
  allocate(uncomputed_tgts(n_elt))

  send_field_name = "girafe"
  recv_field_name = "chinchilla"

  do it = itdeb, itend

    time = (it-itdeb)*dt

    ! Begin time step :
    ! ------------------------------------------------------- To fill in

    ! ---------------------------------------------------- End To fill in

    ! Deform mesh :
    ! The z-coordinate of the mesh is modified with some kind of cosine at
    ! each iteration.
    do i = 1, n_vtx
      coords(3,i) = (coords(1,i)*coords(1,i) + coords(2,i)*coords(2,i)) * ampl*cos(omega*time+phi)
    enddo

    ! Update sent field :
    ! The send field is updated at each iteration. It is an average of
    ! the z-coordoninates of the nodes around the cell center.
    do  i = 1, n_elt
      send_field_data(i) = 0.
      do j = elt_vtx_idx(i)+1, elt_vtx_idx(i+1)
        send_field_data(i) = send_field_data(i) + coords(3,elt_vtx(j))
      enddo
      send_field_data(i) = send_field_data(i)/ (elt_vtx_idx(i+1) - elt_vtx_idx(i))
    enddo

    ! User defined targets :
    ! Targets at a different position than the mesh nodes are used for code1.
    ! Here it is chosen for them to be at the cell center but they could be
    ! at random places. Since here they depend on the mesh node position
    ! they need to be updated when the mesh coordinates are updated.
    do i = 1, n_elt
       xyz_dest(:,i) = 0.
       do j = elt_vtx_idx(i)+1, elt_vtx_idx(i+1)
           xyz_dest(:,i) = xyz_dest(:,i) + coords(:,elt_vtx(j))
       enddo
       xyz_dest(:,i) = xyz_dest(:,i) / (elt_vtx_idx(i+1) - elt_vtx_idx(i))
    enddo

    if (it == itdeb) then

      ! Set the mesh :
      ! Here the appropriate CWIPI functions should be called to
      ! set the arrays difining the mesh.
      ! ------------------------------------------------------- To fill in

      ! ---------------------------------------------------- End To fill in

      ! Set user defined targets :
      ! As explained, the xyz_dest array is filled with user targets (a set of
      ! vertices that are not located on the mesh nodes). Use CWP_User_tgt_pts_set
      ! to set those.
      ! ------------------------------------------------------- To fill in

      ! ---------------------------------------------------- End To fill in

      ! Finalize the mesh :
      ! ------------------------------------------------------- To fill in

      ! ---------------------------------------------------- End To fill in

      n_components = 1
      ! Set fields :
      ! code1 sends the "girafe" field and receives "chinchilla". The field
      ! that is sent is not on the user targets but on the mesh. The field
      ! that is received is located on the user targets.
      ! ------------------------------------------------------- To fill in

      ! ---------------------------------------------------- End To fill in


      ! Set tolerance to 10%.
      ! ------------------------------------------------------- To fill in

      ! ---------------------------------------------------- End To fill in

    else

      ! Find the correct CWIPI function to inform CWIPI that the
      ! the time step has been updated.
      ! ------------------------------------------------------- To fill in

      ! ---------------------------------------------------- End To fill in

    endif

    ! Compute interpolation weights :
    ! ------------------------------------------------------- To fill in

    ! ---------------------------------------------------- End To fill in

    ! Exchange field values between codes :
    ! Use the CWIPI exchange functions similar to the MPI ones
    ! for code1 to send "girafe" and receive "chinchilla".
    ! ------------------------------------------------------- To fill in

    ! ---------------------------------------------------- End To fill in

    ! Check interpolation :
    ! For the receiving field "girafe", check if all  points have been located.
    ! ------------------------------------------------------- To fill in
    n_uncomputed_tgts = 0

    if (n_uncomputed_tgts /= 0) then

    endif
    ! ---------------------------------------------------- End To fill in

    ! End time step :
    ! ------------------------------------------------------- To fill in

    ! ---------------------------------------------------- End To fill in


  enddo

  ! Delete fields :
  ! ------------------------------------------------------- To fill in

  ! ---------------------------------------------------- End To fill in

  ! Delete Mesh :
  ! ------------------------------------------------------- To fill in

  ! ---------------------------------------------------- End To fill in

  ! Delete the coupling :
  ! ------------------------------------------------------- To fill in

  ! ---------------------------------------------------- End To fill in

  ! free
  deallocate(code_names)
  if (associated(intra_comms)) then
    deallocate(intra_comms)
  endif
  deallocate(coupled_code_names)
  deallocate(send_field_data)
  deallocate(recv_field_data)
  deallocate(xyz_dest)
  deallocate(uncomputed_tgts)

  call pdm_fortran_free_c(c_loc(coords))
  call pdm_fortran_free_c(c_loc(elt_vtx_idx))
  call pdm_fortran_free_c(c_loc(elt_vtx))

  ! Finalize CWIPI :
  ! ------------------------------------------------------- To fill in

  ! ---------------------------------------------------- End To fill in

  ! Finalize MPI :
  call MPI_Finalize(ierr)

end program fortran_new_api_deformable_sol
