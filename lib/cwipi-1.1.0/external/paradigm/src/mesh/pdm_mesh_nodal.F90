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

module pdm_mesh_nodal

  use pdm
  use iso_c_binding

  implicit none

  interface PDM_Mesh_nodal_n_vtx_elt_get ; module procedure &
  PDM_Mesh_nodal_n_vtx_elt_get_
  end interface

  private :: PDM_Mesh_nodal_n_vtx_elt_get_

  interface

    !>
    !!
    !! \brief Get the number of vertices of an element type
    !!
    !! \param [in]   type     Element type
    !! \param [in]   comm     Element order
    !!
    !! \return       Number of vertices
    !!

    function PDM_Mesh_nodal_n_vtx_elt_get_cf(elt_t, &
                                             order) &

      result(n_vtx_per_elt) &
      bind (c, name = 'PDM_Mesh_nodal_n_vtx_elt_get')

      use iso_c_binding
      implicit none

      integer(c_int), value :: elt_t, order
      integer(c_int)        :: n_vtx_per_elt

    end function PDM_Mesh_nodal_n_vtx_elt_get_cf

  end interface

  contains

    !>
    !!
    !! \brief Get the number of vertices of an element type
    !!
    !! \param [in]   type           Element type
    !! \param [in]   comm           Element order
    !!
    !! \param [out]  n_vtx_per_elt  Number of vertices per element
    !!

    subroutine PDM_Mesh_nodal_n_vtx_elt_get_(elt_t, &
                                             order, &
                                             n_vtx_per_elt)

      use iso_c_binding
      implicit none

      integer, intent(in)  :: elt_t, order
      integer, intent(out) :: n_vtx_per_elt

      integer :: c_elt_t, c_order, c_n_vtx_per_elt

      c_elt_t = elt_t
      c_order = order

      c_n_vtx_per_elt = PDM_Mesh_nodal_n_vtx_elt_get_cf(c_elt_t, &
                                                        c_order)

      n_vtx_per_elt = c_n_vtx_per_elt

    end subroutine PDM_Mesh_nodal_n_vtx_elt_get_


!> Create a Mesh nodal structure
!!
!! @param[in]   n_part   Number of partition on the current process
!! @param[in]   f_comm   MPI communicator
!! @param[out]  mesh     Pointer to \ref PDM_Mesh_nodal object
!!

  subroutine PDM_mesh_nodal_create (mesh, n_part, f_comm)
    use iso_c_binding

    implicit none

    type(c_ptr)          :: mesh
    integer, intent(in)  :: n_part
    integer, intent(in)  :: f_comm

    integer(c_int)       :: c_comm

    interface
      function PDM_mesh_nodal_create_c (n_part, c_comm) result(mesh) bind(c, name='PDM_Mesh_nodal_create')
        use iso_c_binding
        implicit none
        integer(c_int), value :: n_part
        integer(c_int), value :: c_comm
        type(c_ptr)           :: mesh
      end function PDM_mesh_nodal_create_c
    end interface

    c_comm = PDM_MPI_Comm_f2c(f_comm)

    mesh = PDM_mesh_nodal_create_c (n_part, c_comm)

  end subroutine PDM_mesh_nodal_create


!> \brief Free partially a nodal mesh structure
!!
!! @param[in]  mesh       Pointer to \ref PDM_Mesh_nodal object
!!

  subroutine PDM_mesh_nodal_partial_free (mesh)
    use iso_c_binding

    implicit none

    type(c_ptr), value     :: mesh

    interface
      subroutine PDM_mesh_nodal_partial_free_c (mesh) &
        bind(c, name='PDM_Mesh_nodal_partial_free')

        use iso_c_binding

        implicit none

        type(c_ptr), value     :: mesh
      end subroutine PDM_mesh_nodal_partial_free_c
    end interface

    call PDM_mesh_nodal_partial_free_c (mesh)

  end subroutine PDM_mesh_nodal_partial_free


!> Free a nodal mesh structure
!!
!! @param[in]  mesh       Pointer to \ref PDM_Mesh_nodal object
!!

  subroutine PDM_mesh_nodal_free (mesh)
    use iso_c_binding

    implicit none

    type(c_ptr), value :: mesh

    interface
      subroutine PDM_mesh_nodal_free_c (mesh) &
        bind(c, name='PDM_Mesh_nodal_free')

        use iso_c_binding

        implicit none

        type(c_ptr), value :: mesh
      end subroutine PDM_mesh_nodal_free_c
    end interface

    call PDM_mesh_nodal_free_c (mesh)

  end subroutine PDM_mesh_nodal_free


!> Define partition vertices
!!
!! @param[in]  mesh     Pointer to \ref PDM_Mesh_nodal object
!! @param[in]  id_part  Partition identifier
!! @param[in]  n_vtx    Number of vertices
!! @param[in]  coords   Interlaced coordinates (size = 3 * \ref n_vtx)
!! @param[in]  numabs   Global numbering
!!

  subroutine PDM_mesh_nodal_coord_set (mesh, id_part, n_vtx, coords, numabs, owner)
    use iso_c_binding

    implicit none

    type(c_ptr), value                  :: mesh
    integer, intent(in)                 :: id_part
    integer, intent(in)                 :: n_vtx
    integer, intent(in)                 :: owner
    double precision, pointer           :: coords(:,:)
    integer (pdm_g_num_s), pointer      :: numabs(:)
    type(c_ptr) :: c_coords
    type(c_ptr) :: c_numabs

    interface
      subroutine PDM_mesh_nodal_coord_set_c (mesh, id_part, n_vtx, coords, numabs, owner) &
        bind(c, name='PDM_Mesh_nodal_coord_set')

        use iso_c_binding
        use pdm

        implicit none

        type(c_ptr),                 value :: mesh
        integer (c_int), intent(in), value :: id_part
        integer (c_int), intent(in), value :: n_vtx
        integer (c_int), intent(in), value :: owner
        type (c_ptr),                value :: coords
        type (c_ptr),                value :: numabs
      end subroutine PDM_mesh_nodal_coord_set_c
    end interface

    c_coords = C_NULL_PTR
    if (associated(coords)) then
      c_coords = c_loc (coords)
    endif
      
    c_numabs = C_NULL_PTR
    if (associated(numabs)) then
      c_numabs = c_loc (numabs)
    endif
      

    call  PDM_mesh_nodal_coord_set_c (mesh, id_part, n_vtx, c_coords, c_numabs, owner)

  end subroutine

!> Define standard 3D cells by cell-vertex connectivity
!!
!! @param[in]  mesh          Pointer to \ref PDM_Mesh_nodal object
!! @param[in]  id_part       Partition identifier
!! @param[in]  n_cell        Number of cells
!! @param[in]  cell_vtx_idx  Index of cell vertex connectivity
!! @param[in]  cell_vtx_nb   Number of vertices for each cell
!! @param[in]  cell_vtx      Cell vertex connectivity
!! @param[in]  numabs        Global numbering
!! @param[in]  owner         Ownership

  subroutine PDM_Mesh_nodal_cells_cellvtx_add (mesh, id_part, n_cell, cell_vtx_idx, cell_vtx_nb, cell_vtx, numabs, owner)
    use iso_c_binding

    implicit none

    type(c_ptr), value                  :: mesh
    integer, intent(in)                 :: id_part
    integer, intent(in)                 :: n_cell
    integer, intent(in)                 :: owner
    integer (pdm_l_num_s), pointer      :: cell_vtx_idx(:)
    integer (pdm_l_num_s), pointer      :: cell_vtx_nb(:)
    integer (pdm_l_num_s), pointer      :: cell_vtx(:)
    integer (pdm_g_num_s), pointer      :: numabs(:)
    type(c_ptr) :: c_cell_vtx_idx
    type(c_ptr) :: c_cell_vtx_nb
    type(c_ptr) :: c_cell_vtx
    type(c_ptr) :: c_numabs

    interface
      subroutine PDM_Mesh_nodal_cells_cellvtx_add_c (mesh, id_part, n_cell, cell_vtx_idx, cell_vtx_nb, cell_vtx, numabs, owner) &
        bind(c, name='PDM_Mesh_nodal_cells_cellvtx_add')

        use iso_c_binding
        use pdm

        implicit none

        type(c_ptr),                 value :: mesh
        integer (c_int), intent(in), value :: id_part
        integer (c_int), intent(in), value :: n_cell
        integer (c_int), intent(in), value :: owner
        type (c_ptr),                value :: cell_vtx_idx
        type (c_ptr),                value :: cell_vtx_nb
        type (c_ptr),                value :: cell_vtx
        type (c_ptr),                value :: numabs
      end subroutine PDM_Mesh_nodal_cells_cellvtx_add_c
    end interface

    c_cell_vtx_idx = C_NULL_PTR
    if (associated(cell_vtx_idx)) then
      c_cell_vtx_idx = c_loc (cell_vtx_idx)
    endif
      
    c_cell_vtx_nb = C_NULL_PTR
    if (associated(cell_vtx_nb)) then
      c_cell_vtx_nb = c_loc (cell_vtx_nb)
    endif
      
    c_cell_vtx = C_NULL_PTR
    if (associated(cell_vtx)) then
      c_cell_vtx = c_loc (cell_vtx)
    endif
      
    c_numabs = C_NULL_PTR
    if (associated(numabs)) then
      c_numabs = c_loc (numabs)
    endif  

    call  PDM_Mesh_nodal_cells_cellvtx_add_c (mesh, id_part, n_cell, c_cell_vtx_idx, c_cell_vtx_nb, c_cell_vtx, c_numabs, owner)

  end subroutine

!> Define faces by face-vertex connectivity
!!
!! @param[in]  mesh          Pointer to \ref PDM_Mesh_nodal object
!! @param[in]  id_part       Partition identifier
!! @param[in]  n_face        Number of faces
!! @param[in]  face_vtx_idx  Index of face vertex connectivity
!! @param[in]  face_vtx_nb   Number of vertices for each face
!! @param[in]  face_vtx      Face vertex connectivity
!! @param[in]  numabs        Global numbering
!! @param[in]  owner         Ownership

  subroutine PDM_Mesh_nodal_faces_facevtx_add (mesh, id_part, n_face, face_vtx_idx, face_vtx_nb, face_vtx, numabs, owner)
    use iso_c_binding

    implicit none

    type(c_ptr), value                  :: mesh
    integer, intent(in)                 :: id_part
    integer, intent(in)                 :: n_face
    integer, intent(in)                 :: owner
    integer (pdm_l_num_s), pointer      :: face_vtx_idx(:)
    integer (pdm_l_num_s), pointer      :: face_vtx_nb(:)
    integer (pdm_l_num_s), pointer      :: face_vtx(:)
    integer (pdm_g_num_s), pointer      :: numabs(:)
    type(c_ptr) :: c_face_vtx_idx
    type(c_ptr) :: c_face_vtx_nb
    type(c_ptr) :: c_face_vtx
    type(c_ptr) :: c_numabs

    interface
      subroutine PDM_Mesh_nodal_faces_facevtx_add_c (mesh, id_part, n_face, face_vtx_idx, face_vtx_nb, face_vtx, numabs, owner) &
        bind(c, name='PDM_Mesh_nodal_faces_facevtx_add')

        use iso_c_binding
        use pdm

        implicit none

        type(c_ptr),                 value :: mesh
        integer (c_int), intent(in), value :: id_part
        integer (c_int), intent(in), value :: n_face
        integer (c_int), intent(in), value :: owner
        type (c_ptr),                value :: face_vtx_idx
        type (c_ptr),                value :: face_vtx_nb
        type (c_ptr),                value :: face_vtx
        type (c_ptr),                value :: numabs
      end subroutine PDM_Mesh_nodal_faces_facevtx_add_c
    end interface

    c_face_vtx_idx = C_NULL_PTR
    if (associated(face_vtx_idx)) then
      c_face_vtx_idx = c_loc (face_vtx_idx)
    endif
      
    c_face_vtx_nb = C_NULL_PTR
    if (associated(face_vtx_nb)) then
      c_face_vtx_nb = c_loc (face_vtx_nb)
    endif
      
    c_face_vtx = C_NULL_PTR
    if (associated(face_vtx)) then
      c_face_vtx = c_loc (face_vtx)
    endif
      
    c_numabs = C_NULL_PTR
    if (associated(numabs)) then
      c_numabs = c_loc (numabs)
    endif
      

    call  PDM_Mesh_nodal_faces_facevtx_add_c (mesh, id_part, n_face, c_face_vtx_idx, c_face_vtx_nb, c_face_vtx, c_numabs, owner)

  end subroutine


!> Define 2D cells by cell-edge connectivity
!!
!! @param[in]  mesh          Pointer to \ref PDM_Mesh_nodal object
!! @param[in]  id_part       Partition identifier
!! @param[in]  n_elt         Number of polyhedra
!! @param[in]  n_edge        Number of edges used to describe polyhedra
!! @param[in]  edge_vtx_idx  Index of edge vertex connectivity
!! @param[in]  edge_vtx_nb   Number of vertices for each edge
!! @param[in]  edge_vtx      Edge vertex connectivity
!! @param[in]  cell_edge_idx Index of cell edge connectivity
!! @param[in]  cell_edge_nb  Number of edges for each cell
!! @param[in]  cell_edge     Cell edge connectivity
!! @param[in]  numabs        Global numbering
!! @param[in]  owner         Ownership

subroutine PDM_Mesh_nodal_cell2d_celledge_add (mesh, id_part, n_elt, n_edge, edge_vtx_idx, edge_vtx_nb, edge_vtx, &
                                               cell_edge_idx, cell_edge_nb, cell_edge, numabs, owner)
    use iso_c_binding

    implicit none

    type(c_ptr), value                  :: mesh
    integer, intent(in)                 :: id_part
    integer, intent(in)                 :: n_elt
    integer, intent(in)                 :: n_edge
    integer, intent(in)                 :: owner
    integer (pdm_l_num_s), pointer      :: edge_vtx_idx(:)
    integer (pdm_l_num_s), pointer      :: edge_vtx_nb(:)
    integer (pdm_l_num_s), pointer      :: edge_vtx(:)
    integer (pdm_l_num_s), pointer      :: cell_edge_idx(:)
    integer (pdm_l_num_s), pointer      :: cell_edge_nb(:)
    integer (pdm_l_num_s), pointer      :: cell_edge(:)
    integer (pdm_g_num_s), pointer      :: numabs(:)
    type(c_ptr) :: c_edge_vtx_idx
    type(c_ptr) :: c_edge_vtx_nb
    type(c_ptr) :: c_edge_vtx
    type(c_ptr) :: c_cell_edge_idx
    type(c_ptr) :: c_cell_edge_nb
    type(c_ptr) :: c_cell_edge
    type(c_ptr) :: c_numabs

    interface
      subroutine PDM_Mesh_nodal_cell2d_celledge_add_c (mesh, id_part, n_elt, n_edge, edge_vtx_idx, edge_vtx_nb, edge_vtx, &
                                                       cell_edge_idx, cell_edge_nb, cell_edge, numabs, owner) &
        bind(c, name='PDM_Mesh_nodal_cell2d_celledge_add')

        use iso_c_binding
        use pdm

        implicit none

        type(c_ptr),                 value :: mesh
        integer (c_int), intent(in), value :: id_part
        integer (c_int), intent(in), value :: n_elt
        integer (c_int), intent(in), value :: n_edge
        integer (c_int), intent(in), value :: owner
        type (c_ptr),                value :: edge_vtx_idx
        type (c_ptr),                value :: edge_vtx_nb
        type (c_ptr),                value :: edge_vtx
        type (c_ptr),                value :: cell_edge_idx
        type (c_ptr),                value :: cell_edge_nb
        type (c_ptr),                value :: cell_edge
        type (c_ptr),                value :: numabs
      end subroutine PDM_Mesh_nodal_cell2d_celledge_add_c
    end interface

    c_edge_vtx_idx = C_NULL_PTR
    if (associated(edge_vtx_idx)) then
      c_edge_vtx_idx = c_loc (edge_vtx_idx)
    endif
      
    c_edge_vtx_nb = C_NULL_PTR
    if (associated(edge_vtx_nb)) then
      c_edge_vtx_nb = c_loc (edge_vtx_nb)
    endif
      
    c_edge_vtx = C_NULL_PTR
    if (associated(edge_vtx)) then
      c_edge_vtx = c_loc (edge_vtx)
    endif
      
    c_cell_edge_idx = C_NULL_PTR
    if (associated(cell_edge_idx)) then
      c_cell_edge_idx = c_loc (cell_edge_idx)
    endif
      
    c_cell_edge_nb = C_NULL_PTR
    if (associated(cell_edge_nb)) then
      c_cell_edge_nb = c_loc (cell_edge_nb)
    endif
      
    c_cell_edge = C_NULL_PTR
    if (associated(cell_edge)) then
      c_cell_edge = c_loc (cell_edge)
    endif
      
    c_numabs = C_NULL_PTR
    if (associated(numabs)) then
      c_numabs = c_loc (numabs)
    endif
      

    call  PDM_Mesh_nodal_cell2d_celledge_add_c (mesh, id_part, n_elt, n_edge, c_edge_vtx_idx, c_edge_vtx_nb, c_edge_vtx, &
                                                c_cell_edge_idx, c_cell_edge_nb, c_cell_edge, c_numabs, owner)

  end subroutine

!> Define 3D cells by cell-face connectivity
!!
!! @param[in]  mesh          Pointer to \ref PDM_Mesh_nodal object
!! @param[in]  id_part       Partition identifier
!! @param[in]  n_elt         Number of polyhedra
!! @param[in]  n_face        Number of faces used to describe polyhedra
!! @param[in]  face_vtx_idx  Index of face vertex connectivity
!! @param[in]  face_vtx_nb   Number of vertices for each face
!! @param[in]  face_vtx      Face vertex connectivity
!! @param[in]  face_ln_to_gn Global face numbering
!! @param[in]  cell_face_idx Index of cell face connectivity
!! @param[in]  cell_face_nb  Number of faces for each cell
!! @param[in]  cell_face     Cell face connectivity
!! @param[in]  numabs        Global numbering
!! @param[in]  owner         Ownership

subroutine PDM_Mesh_nodal_cell3d_cellface_add (mesh, id_part, n_elt, n_face, face_vtx_idx, face_vtx_nb, face_vtx, &
                                               face_ln_to_gn, cell_face_idx, cell_face_nb, cell_face, numabs, owner)
    use iso_c_binding

    implicit none

    type(c_ptr), value                  :: mesh
    integer, intent(in)                 :: id_part
    integer, intent(in)                 :: n_elt
    integer, intent(in)                 :: n_face
    integer, intent(in)                 :: owner
    integer (pdm_l_num_s), pointer      :: face_vtx_idx(:)
    integer (pdm_l_num_s), pointer      :: face_vtx_nb(:)
    integer (pdm_l_num_s), pointer      :: face_vtx(:)
    integer (pdm_g_num_s), pointer      :: face_ln_to_gn(:)
    integer (pdm_l_num_s), pointer      :: cell_face_idx(:)
    integer (pdm_l_num_s), pointer      :: cell_face_nb(:)
    integer (pdm_l_num_s), pointer      :: cell_face(:)
    integer (pdm_g_num_s), pointer      :: numabs(:)
    type(c_ptr) :: c_face_vtx_idx
    type(c_ptr) :: c_face_vtx_nb
    type(c_ptr) :: c_face_vtx
    type(c_ptr) :: c_face_ln_to_gn
    type(c_ptr) :: c_cell_face_idx
    type(c_ptr) :: c_cell_face_nb
    type(c_ptr) :: c_cell_face
    type(c_ptr) :: c_numabs

    interface
      subroutine PDM_Mesh_nodal_cell3d_cellface_add_c (mesh, id_part, n_elt, n_face, face_vtx_idx, face_vtx_nb, face_vtx, &
                                                       face_ln_to_gn, cell_face_idx, cell_face_nb, cell_face, numabs, owner) &
        bind(c, name='PDM_Mesh_nodal_cell3d_cellface_add')

        use iso_c_binding
        use pdm

        implicit none

        type(c_ptr),                 value :: mesh
        integer (c_int), intent(in), value :: id_part
        integer (c_int), intent(in), value :: n_elt
        integer (c_int), intent(in), value :: n_face
        integer (c_int), intent(in), value :: owner
        type (c_ptr),                value :: face_vtx_idx
        type (c_ptr),                value :: face_vtx_nb
        type (c_ptr),                value :: face_vtx
        type (c_ptr),                value :: face_ln_to_gn
        type (c_ptr),                value :: cell_face_idx
        type (c_ptr),                value :: cell_face_nb
        type (c_ptr),                value :: cell_face
        type (c_ptr),                value :: numabs
      end subroutine PDM_Mesh_nodal_cell3d_cellface_add_c
    end interface

    c_face_vtx_idx = C_NULL_PTR
    if (associated(face_vtx_idx)) then
      c_face_vtx_idx = c_loc (face_vtx_idx)
    endif 
      
    c_face_vtx_nb = C_NULL_PTR
    if (associated(face_vtx_nb)) then
      c_face_vtx_nb = c_loc (face_vtx_nb)
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
      
    c_cell_face_nb = C_NULL_PTR
    if (associated(cell_face_nb)) then
      c_cell_face_nb = c_loc (cell_face_nb)
    endif 
      
    c_cell_face = C_NULL_PTR
    if (associated(cell_face)) then
      c_cell_face = c_loc (cell_face)
    endif 
      
    c_numabs = C_NULL_PTR
    if (associated(numabs)) then
      c_numabs = c_loc (numabs)
    endif 
      

    call  PDM_Mesh_nodal_cell3d_cellface_add_c (mesh, id_part, n_elt, n_face, c_face_vtx_idx, c_face_vtx_nb, c_face_vtx, &
                                                c_face_ln_to_gn, c_cell_face_idx, c_cell_face_nb, c_cell_face, c_numabs, owner)

  end subroutine

!> Get cell-vertex connectivity
!!
!! @param[in]  mesh          Pointer to \ref PDM_Mesh_nodal object
!! @param[in]  id_part       Partition identifier
!! @param[out] cell_vtx_idx  Index of cell vertex connectivity
!! @param[out] cell_face     Cell vertex connectivity

  subroutine PDM_Mesh_nodal_cell_vtx_connectivity_get (mesh, id_part, cell_vtx_idx, cell_vtx)
    use iso_c_binding

    implicit none

    type(c_ptr), value             :: mesh
    integer, intent(in)            :: id_part
    integer (pdm_l_num_s), pointer :: cell_vtx_idx(:)
    integer (pdm_l_num_s), pointer :: cell_vtx(:)
    type(c_ptr)    :: c_cell_vtx_idx
    type(c_ptr)    :: c_cell_vtx
    integer(c_int) :: n_elt

    interface

      subroutine PDM_Mesh_nodal_cell_vtx_connectivity_get_c (mesh, id_part, cell_vtx_idx, cell_vtx) &
        bind(c, name='PDM_Mesh_nodal_cell_vtx_connectivity_get')

        use iso_c_binding
        use pdm

        implicit none

        type(c_ptr),                 value :: mesh
        integer (c_int), intent(in), value :: id_part
        type (c_ptr)                       :: cell_vtx_idx
        type (c_ptr)                       :: cell_vtx
      end subroutine

      function PDM_Mesh_nodal_n_cell_get_c (mesh, id_part) result (n_elt) &
        bind(c, name='PDM_Mesh_nodal_n_cell_get')

        use iso_c_binding
        use pdm

        implicit none

        type(c_ptr),                 value :: mesh
        integer (c_int), intent(in), value :: id_part
        integer(c_int)                     :: n_elt
      end function

    end interface

    n_elt = PDM_Mesh_nodal_n_cell_get_c(mesh, id_part)

    c_cell_vtx_idx = C_NULL_PTR
    c_cell_vtx = C_NULL_PTR
    call  PDM_Mesh_nodal_cell_vtx_connectivity_get_c (mesh, id_part, c_cell_vtx_idx, c_cell_vtx)

    call c_f_pointer(c_cell_vtx_idx,   &
                     cell_vtx_idx,     &
                     [n_elt+1])

    call c_f_pointer(c_cell_vtx,   &
                     cell_vtx,     &
                     [cell_vtx_idx(n_elt+1)])

  end subroutine

end module PDM_mesh_nodal
