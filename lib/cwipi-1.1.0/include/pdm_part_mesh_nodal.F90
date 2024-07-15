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

module PDM_part_mesh_nodal

  use pdm
  use pdm_mesh_nodal
  use iso_c_binding

  implicit none

  interface PDM_part_mesh_nodal_section_n_elt_get ; module procedure &
  PDM_part_mesh_nodal_section_n_elt_get_
  end interface

  interface PDM_part_mesh_nodal_section_std_get ; module procedure &
  PDM_part_mesh_nodal_section_std_get_
  end interface

  interface PDM_part_mesh_nodal_vtx_g_num_get ; module procedure &
  PDM_part_mesh_nodal_vtx_g_num_get_
  end interface

  interface PDM_part_mesh_nodal_n_vtx_get ; module procedure &
  PDM_part_mesh_nodal_n_vtx_get_
  end interface

  interface PDM_part_mesh_nodal_section_elt_type_get ; module procedure &
  PDM_part_mesh_nodal_section_elt_type_get_
  end interface

  private :: PDM_part_mesh_nodal_section_n_elt_get_
  private :: PDM_part_mesh_nodal_section_std_get_
  private :: PDM_part_mesh_nodal_vtx_g_num_get_
  private :: PDM_part_mesh_nodal_n_vtx_get_
  private :: PDM_part_mesh_nodal_section_elt_type_get_

  interface

    !>
    !!
    !! \brief Get number of section elements
    !!
    !! \param [in]  pmn        Pointer to \ref PDM_part_mesh_nodal_t object
    !! \param [in]  i_section  Section identifier
    !! \param [in]  i_part     Partition identifier
    !!
    !! \return      Number of elements
    !!

    function PDM_part_mesh_nodal_section_n_elt_get_cf(pmn,       &
                                                      i_section, &
                                                      i_part)    &

      result(n_elt) &
      bind (c, name = 'PDM_part_mesh_nodal_section_n_elt_get')

      use iso_c_binding
      implicit none

      type (c_ptr), value   :: pmn
      integer(c_int), value :: i_section, i_part

      integer(c_int)        :: n_elt

    end function PDM_part_mesh_nodal_section_n_elt_get_cf

    !>
    !!
    !! \brief Return standard section description
    !!
    !! \param [in]  pmn                     Pointer to \ref PDM_part_mesh_nodal_t object
    !! \param [in]  i_section               Section identifier
    !! \param [in]  i_part                  Partition identifier
    !! \param [out] connec                  Connectivity
    !! \param [out] numabs                  Global numbering
    !! \param [out] parent_num              Parent numbering or NULL
    !! \param [out] parent_entity_g_num     Parent global numbering or NULL
    !!

    subroutine PDM_part_mesh_nodal_section_std_get_cf(pmn,                 &
                                                      i_section,           &
                                                      i_part,              &
                                                      connec,              &
                                                      numabs,              &
                                                      parent_num,          &
                                                      parent_entity_g_num, &
                                                      ownership)           &

      bind (c, name = 'PDM_part_mesh_nodal_section_std_get')

      use iso_c_binding
      implicit none

      type (c_ptr),   value :: pmn
      integer(c_int), value :: i_section, i_part, ownership

      type (c_ptr)          :: numabs, parent_entity_g_num
      type (c_ptr)          :: connec, parent_num

    end subroutine PDM_part_mesh_nodal_section_std_get_cf

    !>
    !!
    !! \brief  Return global ids of vertices
    !!
    !! \param [in]  pmn       Pointer to \ref PDM_part_mesh_nodal_t object
    !! \param [in]  id_part   Partition identifier
    !!
    !! \return  Golbal ids of vertices
    !!

    function PDM_part_mesh_nodal_vtx_g_num_get_cf(pmn,    &
                                                  i_part) &

      result(vtx_ln_to_gn) &
      bind (c, name = 'PDM_part_mesh_nodal_vtx_g_num_get')

      use iso_c_binding
      implicit none

      type (c_ptr),   value :: pmn
      integer(c_int), value :: i_part

      type (c_ptr)          :: vtx_ln_to_gn

    end function PDM_part_mesh_nodal_vtx_g_num_get_cf

    !>
    !!
    !! \brief  Return number of vertices
    !!
    !! \param [in]  pmn       Pointer to \ref PDM_part_mesh_nodal_t object
    !! \param [in]  id_part   Partition identifier
    !!
    !! \return  Number of vertices
    !!

    function PDM_part_mesh_nodal_n_vtx_get_cf(pmn,    &
                                              i_part) &

      result(n_vtx) &
      bind (c, name = 'PDM_part_mesh_nodal_n_vtx_get')

      use iso_c_binding
      implicit none

      type (c_ptr),   value :: pmn
      integer(c_int), value :: i_part

      integer(c_int)        :: n_vtx

    end function PDM_part_mesh_nodal_n_vtx_get_cf

    !>
    !!
    !! \brief  Return type of section
    !!
    !! \param [in]  pmn          Pointer to \ref PDM_part_mesh_nodal_t object
    !! \param [in]  i_section    Section identifier
    !!
    !! \return  Type of section
    !!

    function PDM_part_mesh_nodal_section_elt_type_get_cf(pmn,      &
                                                         i_section) &

      result(elt_t) &
      bind (c, name = 'PDM_part_mesh_nodal_section_elt_type_get')

      use iso_c_binding
      implicit none

      type (c_ptr),   value :: pmn
      integer(c_int), value :: i_section

      integer(c_int)        :: elt_t

    end function PDM_part_mesh_nodal_section_elt_type_get_cf

  end interface

  contains



    subroutine PDM_part_mesh_nodal_section_n_elt_get_(pmn,       &
                                                      i_section, &
                                                      i_part,    &
                                                      n_elt)
      ! Get number of section elements
      use iso_c_binding
      implicit none

      type(c_ptr), value   :: pmn       ! Pointer to PDM_part_mesh_nodal_t object
      integer, intent(in)  :: i_section ! Section identifier
      integer, intent(in)  :: i_part    ! Partition identifier

      integer, intent(out) :: n_elt     ! Number of elements

      n_elt = PDM_part_mesh_nodal_section_n_elt_get_cf(pmn,       &
                                                       i_section, &
                                                       i_part)

    end subroutine PDM_part_mesh_nodal_section_n_elt_get_



    subroutine PDM_part_mesh_nodal_section_std_get_(pmn,                 &
                                                    i_section,           &
                                                    i_part,              &
                                                    connec,              &
                                                    numabs,              &
                                                    parent_num,          &
                                                    parent_entity_g_num, &
                                                    ownership)
      ! Return standard section description
      use iso_c_binding
      implicit none

      type (c_ptr),   value          :: pmn                    ! Pointer to PDM_part_mesh_nodal_t object
      integer, intent(in)            :: i_section              ! Section identifier
      integer, intent(in)            :: i_part                 ! Partition identifier
      integer, intent(in)            :: ownership              ! Data ownership

      integer(pdm_l_num_s), pointer  :: connec(:)              ! Connectivity
      integer(pdm_g_num_s), pointer  :: numabs(:)              ! Global ids
      integer(pdm_l_num_s), pointer  :: parent_num(:)          ! Parent local ids or *null()*
      integer(pdm_g_num_s), pointer  :: parent_entity_g_num(:) ! Parent global ids or *null()*

      integer                        :: n_elt, elt_t, n_vtx_per_elt
      integer                        :: order = 1

      type (c_ptr)   :: c_connec
      type (c_ptr)   :: c_numabs
      type (c_ptr)   :: c_parent_entity_g_num
      type (c_ptr)   :: c_parent_num

      c_connec = C_NULL_PTR
      c_numabs = C_NULL_PTR
      c_parent_entity_g_num = C_NULL_PTR
      c_parent_num = C_NULL_PTR
      call PDM_part_mesh_nodal_section_std_get_cf(pmn,                   &
                                                  i_section,             &
                                                  i_part,                &
                                                  c_connec,              &
                                                  c_numabs,              &
                                                  c_parent_num,          &
                                                  c_parent_entity_g_num, &
                                                  ownership)

      n_elt = PDM_part_mesh_nodal_section_n_elt_get_cf(pmn,       &
                                                       i_section, &
                                                       i_part)

      call PDM_part_mesh_nodal_section_elt_type_get(pmn,       &
                                                    i_section, &
                                                    elt_t)

      call PDM_Mesh_nodal_n_vtx_elt_get(elt_t, &
                                        order, &
                                        n_vtx_per_elt)

      call c_f_pointer(c_connec, &
                       connec,   &
                       [n_elt * n_vtx_per_elt])

      call c_f_pointer(c_numabs, &
                       numabs,   &
                       [n_elt])

      parent_num => null()
      if ( c_associated(c_parent_num) ) then
        call c_f_pointer(c_parent_num, &
                         parent_num,   &
                         [n_elt])
      end if

      parent_entity_g_num => null()
      if ( c_associated(c_parent_entity_g_num) ) then
        call c_f_pointer(c_parent_entity_g_num, &
                         parent_entity_g_num,   &
                         [n_elt])
      end if

    end subroutine PDM_part_mesh_nodal_section_std_get_



    subroutine PDM_part_mesh_nodal_vtx_g_num_get_(pmn,          &
                                                  i_part,       &
                                                  vtx_ln_to_gn)
      ! Return global ids of vertices
      use iso_c_binding
      implicit none

      type (c_ptr),   value          :: pmn             ! Pointer to PDM_part_mesh_nodal_t object

      integer, intent(in)            :: i_part          ! Partition identifier
      integer (pdm_g_num_s), pointer :: vtx_ln_to_gn(:) ! Golbal ids of vertices (size = ``n_vtx``)

      type (c_ptr)    :: c_vtx_ln_to_gn = C_NULL_PTR
      integer         :: n_vtx


      c_vtx_ln_to_gn = PDM_part_mesh_nodal_vtx_g_num_get_cf(pmn, &
                                                            i_part)

      n_vtx = PDM_part_mesh_nodal_n_vtx_get_cf(pmn, &
                                               i_part)

      call c_f_pointer(c_vtx_ln_to_gn, &
                       vtx_ln_to_gn,   &
                       [n_vtx])

    end subroutine PDM_part_mesh_nodal_vtx_g_num_get_


    subroutine PDM_part_mesh_nodal_n_vtx_get_(pmn,    &
                                              i_part, &
                                              n_vtx)
      ! Return number of vertices
      use iso_c_binding
      implicit none

      type (c_ptr),   value :: pmn    ! Pointer to PDM_part_mesh_nodal_t object
      integer, intent(in)   :: i_part ! Partition identifier
      integer, intent(out)  :: n_vtx  ! Number of vertices

      n_vtx = PDM_part_mesh_nodal_n_vtx_get_cf(pmn, &
                                                 i_part)

    end subroutine PDM_part_mesh_nodal_n_vtx_get_


    subroutine PDM_part_mesh_nodal_section_elt_type_get_(pmn,       &
                                                         i_section, &
                                                         elt_t)
      ! Return type of section
      use iso_c_binding
      implicit none

      type (c_ptr),   value :: pmn       ! Pointer to PDM_part_mesh_nodal_t object
      integer, intent(in)   :: i_section ! Section identifier
      integer, intent(out)  :: elt_t     ! Type of section

      elt_t = PDM_part_mesh_nodal_section_elt_type_get_cf(pmn, &
                                                          i_section)

    end subroutine PDM_part_mesh_nodal_section_elt_type_get_



  subroutine PDM_part_mesh_nodal_create (mesh, mesh_dimension, n_part, f_comm)
    ! Create a PDM_part_mesh_nodal structure
    use iso_c_binding
    implicit none

    type(c_ptr)          :: mesh           ! Pointer to PDM_part_mesh_nodal object
    integer, intent(in)  :: mesh_dimension ! Mesh dimension
    integer, intent(in)  :: n_part         ! Number of partition on the current process
    integer, intent(in)  :: f_comm         ! MPI communicator

    integer(c_int)       :: c_comm

    interface
      function PDM_part_mesh_nodal_create_c (mesh_dimension, n_part, c_comm) result(mesh) bind(c, name='PDM_part_mesh_nodal_create')
        use iso_c_binding
        implicit none
        integer(c_int), value :: mesh_dimension
        integer(c_int), value :: n_part
        integer(c_int), value :: c_comm
        type(c_ptr)           :: mesh
      end function PDM_part_mesh_nodal_create_c
    end interface

    c_comm = PDM_MPI_Comm_f2c(f_comm)

    mesh = PDM_part_mesh_nodal_create_c (mesh_dimension, n_part, c_comm)

  end subroutine PDM_part_mesh_nodal_create


  subroutine PDM_part_mesh_nodal_partial_free (mesh)
    ! Free partially a PDM_part_mesh_nodal structure
    use iso_c_binding

    implicit none

    type(c_ptr), value     :: mesh ! Pointer to PDM_part_mesh_nodal object

    interface
      subroutine PDM_part_mesh_nodal_partial_free_c (mesh) &
        bind(c, name='PDM_part_mesh_nodal_partial_free')

        use iso_c_binding

        implicit none

        type(c_ptr), value     :: mesh
      end subroutine PDM_part_mesh_nodal_partial_free_c
    end interface

    call PDM_part_mesh_nodal_partial_free_c (mesh)

  end subroutine PDM_part_mesh_nodal_partial_free


  subroutine PDM_part_mesh_nodal_free (mesh)
    ! Free a PDM_part_mesh_nodal structure
    use iso_c_binding

    implicit none

    type(c_ptr), value :: mesh ! Pointer to PDM_part_mesh_nodal object

    interface
      subroutine PDM_part_mesh_nodal_free_c (mesh) &
        bind(c, name='PDM_part_mesh_nodal_free')

        use iso_c_binding

        implicit none

        type(c_ptr), value :: mesh
      end subroutine PDM_part_mesh_nodal_free_c
    end interface

    call PDM_part_mesh_nodal_free_c (mesh)

  end subroutine PDM_part_mesh_nodal_free


  subroutine PDM_part_mesh_nodal_coord_set (mesh, id_part, n_vtx, coords, numabs, owner)
    ! Define partition vertices
    use iso_c_binding

    implicit none

    type(c_ptr), value                  :: mesh        ! Pointer to PDM_part_mesh_nodal object
    integer, intent(in)                 :: id_part     ! Partition identifier
    integer, intent(in)                 :: n_vtx       ! Number of vertices
    integer, intent(in)                 :: owner       ! Data ownership
    double precision, pointer           :: coords(:,:) ! Interlaced coordinates (shape = [3, n_vtx])
    integer (pdm_g_num_s), pointer      :: numabs(:)   ! Global numbering
    type(c_ptr) :: c_coords
    type(c_ptr) :: c_numabs

    interface
      subroutine PDM_part_mesh_nodal_coord_set_c (mesh, id_part, n_vtx, coords, numabs, owner) &
        bind(c, name='PDM_part_mesh_nodal_coord_set')

        use iso_c_binding
        use pdm

        implicit none

        type(c_ptr),                 value :: mesh
        integer (c_int), intent(in), value :: id_part
        integer (c_int), intent(in), value :: n_vtx
        integer (c_int), intent(in), value :: owner
        type (c_ptr),                value :: coords
        type (c_ptr),                value :: numabs
      end subroutine PDM_part_mesh_nodal_coord_set_c
    end interface

    c_coords = C_NULL_PTR
    if (associated(coords)) then
      c_coords = c_loc (coords)
    endif
      
    c_numabs = C_NULL_PTR
    if (associated(numabs)) then
      c_numabs = c_loc (numabs)
    endif    

    call  PDM_part_mesh_nodal_coord_set_c (mesh, id_part, n_vtx, c_coords, c_numabs, owner)

  end subroutine PDM_part_mesh_nodal_coord_set


  subroutine PDM_part_mesh_nodal_cells_cellvtx_add (mesh, id_part, n_cell, cell_vtx_idx, cell_vtx, numabs, owner)
    ! Define standard 3D cells by cell-vertex connectivity
    use iso_c_binding

    implicit none

    type(c_ptr), value                  :: mesh            ! Pointer to PDM_part_mesh_nodal object
    integer, intent(in)                 :: id_part         ! Partition identifier
    integer, intent(in)                 :: n_cell          ! Number of cells
    integer, intent(in)                 :: owner           ! Data ownership
    integer (pdm_l_num_s), pointer      :: cell_vtx_idx(:) ! Index of cell vertex connectivity (size = n_cell + 1)
    integer (pdm_l_num_s), pointer      :: cell_vtx(:)     ! Cell vertex connectivity (size = cell_vtx_idx(n_cell+1))
    integer (pdm_g_num_s), pointer      :: numabs(:)       ! Global numbering
    type(c_ptr) :: c_cell_vtx_idx
    type(c_ptr) :: c_cell_vtx
    type(c_ptr) :: c_numabs

    interface
      subroutine PDM_part_mesh_nodal_cells_cellvtx_add_c (mesh, id_part, n_cell, cell_vtx_idx, cell_vtx, &
        numabs, owner) bind(c, name='PDM_part_mesh_nodal_cells_cellvtx_add')

        use iso_c_binding
        use pdm

        implicit none

        type(c_ptr),                 value :: mesh
        integer (c_int), intent(in), value :: id_part
        integer (c_int), intent(in), value :: n_cell
        integer (c_int), intent(in), value :: owner
        type (c_ptr),                value :: cell_vtx_idx
        type (c_ptr),                value :: cell_vtx
        type (c_ptr),                value :: numabs
      end subroutine PDM_part_mesh_nodal_cells_cellvtx_add_c
    end interface

    c_cell_vtx_idx = C_NULL_PTR
    if (associated (cell_vtx_idx)) then
      c_cell_vtx_idx = c_loc (cell_vtx_idx)
    endif

    c_cell_vtx = C_NULL_PTR
    if (associated (cell_vtx)) then
      c_cell_vtx = c_loc (cell_vtx)
    endif

    c_numabs = C_NULL_PTR
    if (associated (numabs)) then
      c_numabs = c_loc (numabs)
    endif


    call  PDM_part_mesh_nodal_cells_cellvtx_add_c (mesh, id_part, n_cell, c_cell_vtx_idx, c_cell_vtx, &
                                                   c_numabs, owner)

  end subroutine PDM_part_mesh_nodal_cells_cellvtx_add


  subroutine PDM_part_mesh_nodal_faces_facevtx_add (mesh, id_part, n_face, face_vtx_idx, face_vtx, numabs, owner)
    ! Define faces by face-vertex connectivity
    use iso_c_binding

    implicit none

    type(c_ptr), value                  :: mesh            ! Pointer to PDM_part_mesh_nodal object
    integer, intent(in)                 :: id_part         ! Partition identifier
    integer, intent(in)                 :: n_face          ! Number of faces
    integer, intent(in)                 :: owner           ! Data ownership
    integer (pdm_l_num_s), pointer      :: face_vtx_idx(:) ! Index of face vertex connectivity
    integer (pdm_l_num_s), pointer      :: face_vtx(:)     ! Face vertex connectivity
    integer (pdm_g_num_s), pointer      :: numabs(:)       ! Global numbering
    type(c_ptr) :: c_face_vtx_idx
    type(c_ptr) :: c_face_vtx
    type(c_ptr) :: c_numabs

    interface
      subroutine PDM_part_mesh_nodal_faces_facevtx_add_c (mesh, id_part, n_face, face_vtx_idx, face_vtx, &
        numabs, owner) bind(c, name='PDM_part_mesh_nodal_faces_facevtx_add')

        use iso_c_binding
        use pdm

        implicit none

        type(c_ptr),                 value :: mesh
        integer (c_int), intent(in), value :: id_part
        integer (c_int), intent(in), value :: n_face
        integer (c_int), intent(in), value :: owner
        type (c_ptr),                value :: face_vtx_idx
        type (c_ptr),                value :: face_vtx
        type (c_ptr),                value :: numabs
      end subroutine PDM_part_mesh_nodal_faces_facevtx_add_c
    end interface

    c_face_vtx_idx = C_NULL_PTR
    if (associated(face_vtx_idx)) then
      c_face_vtx_idx = c_loc (face_vtx_idx)
    endif
      
    c_face_vtx = C_NULL_PTR
    if (associated(face_vtx)) then
      c_face_vtx = c_loc (face_vtx)
    endif
      
    c_numabs = C_NULL_PTR
    if (associated(numabs)) then
      c_numabs = c_loc (numabs)
    endif
      

    call  PDM_part_mesh_nodal_faces_facevtx_add_c (mesh, id_part, n_face, c_face_vtx_idx, c_face_vtx, &
                                                   c_numabs, owner)

  end subroutine PDM_part_mesh_nodal_faces_facevtx_add


  subroutine PDM_part_mesh_nodal_face2d_faceedge_add (mesh, id_part, n_elt, n_edge, edge_vtx, &
                                                      face_edge_idx, face_edge, numabs, owner)
    ! Define 2D faces by face-edge connectivity
    use iso_c_binding

    implicit none

    type(c_ptr), value                  :: mesh             ! Pointer to PDM_part_mesh_nodal object
    integer, intent(in)                 :: id_part          ! Partition identifier
    integer, intent(in)                 :: n_elt            ! Number of polyhedra
    integer, intent(in)                 :: n_edge           ! Number of edges used to describe polyhedra
    integer, intent(in)                 :: owner            ! Data ownership
    integer (pdm_l_num_s), pointer      :: edge_vtx(:)      ! Edge vertex connectivity
    integer (pdm_l_num_s), pointer      :: face_edge_idx(:) ! Index of face edge connectivity
    integer (pdm_l_num_s), pointer      :: face_edge(:)     ! Face edge connectivity
    integer (pdm_g_num_s), pointer      :: numabs(:)        ! Global numbering
    type(c_ptr) :: c_edge_vtx
    type(c_ptr) :: c_face_edge_idx
    type(c_ptr) :: c_face_edge
    type(c_ptr) :: c_numabs

    interface
      subroutine PDM_part_mesh_nodal_face2d_faceedge_add_c (mesh, id_part, n_elt, n_edge, edge_vtx, &
                                                            face_edge_idx, face_edge, numabs, owner) &
        bind(c, name='PDM_part_mesh_nodal_face2d_faceedge_add')

        use iso_c_binding
        use pdm

        implicit none

        type(c_ptr),                 value :: mesh
        integer (c_int), intent(in), value :: id_part
        integer (c_int), intent(in), value :: n_elt
        integer (c_int), intent(in), value :: n_edge
        integer (c_int), intent(in), value :: owner
        type (c_ptr),                value :: edge_vtx
        type (c_ptr),                value :: face_edge_idx
        type (c_ptr),                value :: face_edge
        type (c_ptr),                value :: numabs
      end subroutine PDM_part_mesh_nodal_face2d_faceedge_add_c
    end interface

    c_edge_vtx = C_NULL_PTR
    if (associated(edge_vtx)) then
      c_edge_vtx = c_loc (edge_vtx)
    endif
      
    c_face_edge_idx = C_NULL_PTR
    if (associated(face_edge_idx)) then
      c_face_edge_idx = c_loc (face_edge_idx)
    endif
      
    c_face_edge = C_NULL_PTR
    if (associated(face_edge)) then
      c_face_edge = c_loc (face_edge)
    endif
      
    c_numabs = C_NULL_PTR
    if (associated(numabs)) then
      c_numabs = c_loc (numabs)
    endif
      

    call  PDM_part_mesh_nodal_face2d_faceedge_add_c (mesh, id_part, n_elt, n_edge, c_edge_vtx, &
                                                     c_face_edge_idx, c_face_edge, c_numabs, owner)

  end subroutine PDM_part_mesh_nodal_face2d_faceedge_add


  subroutine PDM_part_mesh_nodal_cell3d_cellface_add (mesh, id_part, n_elt, n_face, face_vtx_idx, face_vtx, &
                                                      face_ln_to_gn, cell_face_idx, cell_face, numabs, owner)
    ! 3D cells by cell-face connectivity
    use iso_c_binding

    implicit none

    type(c_ptr), value                  :: mesh             ! Pointer to PDM_part_mesh_nodal object
    integer, intent(in)                 :: id_part          ! Partition identifier
    integer, intent(in)                 :: n_elt            ! Number of polyhedra
    integer, intent(in)                 :: n_face           ! Number of faces used to describe polyhedra
    integer, intent(in)                 :: owner            ! Data ownership
    integer (pdm_l_num_s), pointer      :: face_vtx_idx(:)  ! Index of face vertex connectivity
    integer (pdm_l_num_s), pointer      :: face_vtx(:)      ! Face vertex connectivity
    integer (pdm_g_num_s), pointer      :: face_ln_to_gn(:) ! Global face numbering
    integer (pdm_l_num_s), pointer      :: cell_face_idx(:) ! Index of cell face connectivity
    integer (pdm_l_num_s), pointer      :: cell_face(:)     ! Cell face connectivity
    integer (pdm_g_num_s), pointer      :: numabs(:)        ! Global numbering
    type(c_ptr) :: c_face_vtx_idx
    type(c_ptr) :: c_face_vtx
    type(c_ptr) :: c_face_ln_to_gn
    type(c_ptr) :: c_cell_face_idx
    type(c_ptr) :: c_cell_face
    type(c_ptr) :: c_numabs

    interface
      subroutine PDM_part_mesh_nodal_cell3d_cellface_add_c (mesh, id_part, n_elt, n_face, face_vtx_idx, face_vtx, &
                                                            face_ln_to_gn, cell_face_idx, cell_face, numabs, owner) &
        bind(c, name='PDM_part_mesh_nodal_cell3d_cellface_add')

        use iso_c_binding
        use pdm

        implicit none

        type(c_ptr),                 value :: mesh
        integer (c_int), intent(in), value :: id_part
        integer (c_int), intent(in), value :: n_elt
        integer (c_int), intent(in), value :: n_face
        integer (c_int), intent(in), value :: owner
        type (c_ptr),                value :: face_vtx_idx
        type (c_ptr),                value :: face_vtx
        type (c_ptr),                value :: face_ln_to_gn
        type (c_ptr),                value :: cell_face_idx
        type (c_ptr),                value :: cell_face
        type (c_ptr),                value :: numabs
      end subroutine PDM_part_mesh_nodal_cell3d_cellface_add_c
    end interface

    c_face_vtx_idx = C_NULL_PTR
    if (associated(face_vtx_idx)) then
      c_face_vtx_idx = c_loc (face_vtx_idx)
    endif 
      
    c_face_vtx = C_NULL_PTR
    if (associated(face_vtx)) then
      c_face_vtx = c_loc (face_vtx)
    endif 
      
    c_face_ln_to_gn = C_NULL_PTR
    if (associated(face_ln_to_gn)) then
      c_face_ln_to_gn = c_loc (face_ln_to_gn)
    endif 
      
    c_cell_face_idx = C_NULL_PTR
    if (associated(cell_face_idx)) then
      c_cell_face_idx = c_loc (cell_face_idx)
    endif 
      
    c_cell_face = C_NULL_PTR
    if (associated(cell_face)) then
      c_cell_face = c_loc (cell_face)
    endif 
      
    c_numabs = C_NULL_PTR
    if (associated(numabs)) then
      c_numabs = c_loc (numabs)
    endif 
      

    call  PDM_part_mesh_nodal_cell3d_cellface_add_c (mesh, id_part, n_elt, n_face, c_face_vtx_idx, c_face_vtx, &
                                                     c_face_ln_to_gn, c_cell_face_idx, c_cell_face, c_numabs, owner)

  end subroutine PDM_part_mesh_nodal_cell3d_cellface_add


  subroutine PDM_part_mesh_nodal_cell_vtx_connect_get (mesh, id_part, cell_vtx_idx, cell_vtx)
    ! Get cell-vertex connectivity
    use iso_c_binding

    implicit none

    type(c_ptr), value             :: mesh            ! Pointer to PDM_part_mesh_nodal object
    integer, intent(in)            :: id_part         ! Partition identifier
    integer (pdm_l_num_s), pointer :: cell_vtx_idx(:) ! Index of cell vertex connectivity
    integer (pdm_l_num_s), pointer :: cell_vtx(:)     ! Cell vertex connectivity
    type(c_ptr)    :: c_cell_vtx_idx
    type(c_ptr)    :: c_cell_vtx
    integer(c_int) :: geom_kind
    integer(c_int) :: n_elt

    interface

      subroutine PDM_part_mesh_nodal_cell_vtx_connect_get_c (mesh, geom_kind, id_part, cell_vtx_idx, cell_vtx) &
        bind(c, name='PDM_part_mesh_nodal_cell_vtx_connect_get')

        use iso_c_binding
        use pdm

        implicit none

        type(c_ptr),                 value :: mesh
        integer (c_int), intent(in), value :: geom_kind
        integer (c_int), intent(in), value :: id_part
        type (c_ptr)                       :: cell_vtx_idx
        type (c_ptr)                       :: cell_vtx
      end subroutine

      function PDM_part_mesh_nodal_principal_geom_kind_get_c (mesh) result (geom_kind) &
        bind(c, name='PDM_part_mesh_nodal_principal_geom_kind_get')

        use iso_c_binding
        use pdm

        implicit none

        type(c_ptr),    value :: mesh
        integer (c_int)       :: geom_kind
      end function

      function PDM_part_mesh_nodal_n_elmts_get_c (mesh, geom_kind, id_part) result (n_elt) &
        bind(c, name='PDM_part_mesh_nodal_n_elmts_get')

        use iso_c_binding
        use pdm

        implicit none

        type(c_ptr),                 value :: mesh
        integer (c_int), intent(in), value :: geom_kind
        integer (c_int), intent(in), value :: id_part
        integer (c_int)                    :: n_elt
      end function

    end interface

    geom_kind = PDM_part_mesh_nodal_principal_geom_kind_get_c (mesh)
    n_elt = PDM_part_mesh_nodal_n_elmts_get_c (mesh, geom_kind, id_part)

    c_cell_vtx_idx = C_NULL_PTR
    c_cell_vtx = C_NULL_PTR
    call  PDM_part_mesh_nodal_cell_vtx_connect_get_c (mesh, geom_kind, id_part, c_cell_vtx_idx, c_cell_vtx)

    call c_f_pointer(c_cell_vtx_idx,   &
                     cell_vtx_idx,     &
                     [n_elt+1])

    call c_f_pointer(c_cell_vtx,   &
                     cell_vtx,     &
                     [cell_vtx_idx(n_elt+1)])

  end subroutine PDM_part_mesh_nodal_cell_vtx_connect_get

end module PDM_part_mesh_nodal
