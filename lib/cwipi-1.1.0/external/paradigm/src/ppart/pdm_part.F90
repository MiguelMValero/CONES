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

module pdm_part

  use pdm

  implicit none

  integer, parameter :: PDM_PART_RENUM_FACE_RANDOM        = 1
  integer, parameter :: PDM_PART_RENUM_FACE_NONE          = 2
  integer, parameter :: PDM_PART_RENUM_FACE_LEXICOGRAPHIC = 3

  integer, parameter :: PDM_PART_RENUM_CELL_HILBERT = 1
  integer, parameter :: PDM_PART_RENUM_CELL_RANDOM  = 2
  integer, parameter :: PDM_PART_RENUM_CELL_NONE    = 3
  integer, parameter :: PDM_PART_RENUM_CELL_CUTHILL = 4

  interface PDM_part_create ; module procedure  &
    PDM_part_create_
  end interface

  interface PDM_part_part_dim_get ; module procedure  &
    PDM_part_part_dim_get_
  end interface

  interface PDM_part_part_val_get ; module procedure  &
  PDM_part_part_val_get_
  end interface

  interface PDM_part_time_get ; module procedure  &
  PDM_part_time_get_
  end interface

  interface PDM_part_stat_get ; module procedure  &
  PDM_part_stat_get_
  end interface

  interface PDM_part_part_color_get ; module procedure  &
  PDM_part_part_color_get_
  end interface

  interface PDM_part_coarse_mesh_create ; module procedure  &
    PDM_part_coarse_mesh_create_
  end interface

  interface PDM_part_coarse_mesh_input ; module procedure  &
    PDM_part_coarse_mesh_input_
  end interface

  interface PDM_part_coarse_mesh_part_dim_get ; module procedure  &
    PDM_part_coarse_mesh_part_dim_get_
  end interface

  interface PDM_part_coarse_mesh_part_get ; module procedure  &
    PDM_part_coarse_mesh_part_get_
  end interface

  interface PDM_part_coarse_color_get ; module procedure  &
  PDM_part_coarse_color_get_
  end interface

  interface PDM_part_coarse_mesh_time_get ; module procedure  &
  PDM_part_coarse_mesh_time_get_
  end interface

interface

  !>
  !!
  !! \brief Build a initial partitioning
  !!  Build a initial partitioning from :
  !!      - Cell block distribution with implicit global numbering
  !!         (the first cell is the first cell of the first process and
  !!          the latest cell is the latest cell of the latest process)
  !!      - Face block distribution with implicit global numbering
  !!      - Vertex block distribution with implicit global numbering
  !!  To repart an existing partition use \ref PDM_part_repart function
  !!
  !! \param [in]   comm                   MPI Comminicator
  !! \param [in]   split_method           Split method
  !! \param [in]   renum_cell_method      Cell renumbering method
  !! \param [in]   renum_face_method      Face renumbering method
  !! \param [in]   n_property_cell        Number of cell properties
  !! \param [in]   renum_properties_cell  For cache blocking [ n_cell_per_cache_wanted, isAsynchrone, isVectorisation ] \ref PDM_renum_cacheblocking
  !! \param [in]   n_property_face        Number of face properties
  !! \param [in]   renum_properties_face  NOT USED
  !! \param [in]   n_part                 Number of partition to build on this process
  !! \param [in]   dn_cell                Number of distributed cells
  !! \param [in]   dn_face                Number of distributed faces
  !! \param [in]   dn_vtx                 Number of distributed vertices
  !! \param [in]   n_face_group           Number of face groups
  !! \param [in]   dcell_face_idx         Distributed cell face connectivity index or NULL
  !!                                      (size : dn_cell + 1, numbering : 0 to n-1)
  !! \param [in]   dcell_face             Distributed cell face connectivity or NULL
  !!                                      (size : dface_vtx_idx[dn_cell], numbering : 1 to n)
  !! \param [in]   dcell_tag              Cell tag (size : n_cell) or NULL
  !! \param [in]   dcell_weight           Cell weight (size : n_cell) or NULL
  !! \param [in]   have_dcell_part        Presence of an array of cell part id
  !! \param [in]   dcell_part             Distributed cell partitioning
  !!                                      (size = dn_cell) or NULL (No partitioning if != NULL)
  !! \param [in]   dface_cell             Distributed face cell connectivity or NULL
  !!                                      (size : 2 * dn_face, numbering : 1 to n)
  !! \param [in]   dface_vtx_idx          Distributed face to vertex connectivity index
  !!                                      (size : dn_face + 1, numbering : 0 to n-1)
  !! \param [in]   dface_vtx              Distributed face to vertex connectivity
  !!                                      (size : dface_vtx_idx[dn_face], numbering : 1 to n)
  !! \param [in]   dface_tag              Distributed face tag (size : dn_face)
  !!                                      or NULL
  !! \param [in]   dvtx_coord             Distributed vertex coordinates
  !!                                      (size : 3*dn_vtx)
  !! \param [in]   dvtx_tag               Distributed vertex tag (size : dn_vtx) or NULL
  !! \param [in]   dface_group_idx        Index of distributed faces list of each group
  !!                                      (size = n_face_group + 1) or NULL
  !! \param [in]   dface_group            Distributed faces list of each group
  !!                                      (size = dface_group[dface_group_idx[n_face_group]], numbering : 1 to n)
  !!                                      or NULL
  !! \return    Pointer to \ref PDM_part object
  !!
  !!

  function PDM_part_create_c (comm, &
                              split_method, &
                              renum_cell_method, &
                              renum_face_method, &
                              n_property_cell, &
                              renum_properties_cell, &
                              n_property_face, &
                              renum_properties_face, &
                              n_part, &
                              dn_cell, &
                              dn_face, &
                              dn_vtx, &
                              n_face_group, &
                              dcell_faceIdx, &
                              dcell_face, &
                              dcell_tag, &
                              dcell_weight, &
                              have_dcell_part, &
                              dcell_part, &
                              dface_cell, &
                              dface_vtx_idx, &
                              dface_vtx, &
                              dface_tag, &
                              dvtx_coord, &
                              dvtx_tag, &
                              dface_group_idx, &
                              dface_group) &
  result (ppart) &
  bind (c, name='PDM_part_create')

    use iso_c_binding
    implicit none

    type(c_ptr)            :: ppart
    integer(c_int), value  :: comm
    integer(c_int), value  :: split_method
    character(kind=c_char) :: renum_cell_method(*)
    character(kind=c_char) :: renum_face_method(*)
    integer(c_int), value  :: n_property_cell
    type(c_ptr),    value  :: renum_properties_cell
    integer(c_int), value  :: n_property_face
    type(c_ptr),    value  :: renum_properties_face
    integer(c_int), value  :: n_part
    integer(c_int), value  :: dn_cell
    integer(c_int), value  :: dn_face
    integer(c_int), value  :: dn_vtx
    integer(c_int), value  :: n_face_group
    type(c_ptr),    value  :: dcell_faceIdx
    type(c_ptr),    value  :: dcell_face
    type(c_ptr),    value  :: dcell_tag
    type(c_ptr),    value  :: dcell_weight
    integer(c_int), value  :: have_dcell_part
    type(c_ptr),    value  :: dcell_part
    type(c_ptr),    value  :: dface_cell
    type(c_ptr),    value  :: dface_vtx_idx
    type(c_ptr),    value  :: dface_vtx
    type(c_ptr),    value  :: dface_tag
    type(c_ptr),    value  :: dvtx_coord
    type(c_ptr),    value  :: dvtx_tag
    type(c_ptr),    value  :: dface_group_idx
    type(c_ptr),    value  :: dface_group

  end function PDM_part_create_c


  !>
  !!
  !! \brief Return a mesh partition dimensions
  !!
  !! \param [in]   ppart               Pointer to \ref PDM_part object
  !! \param [in]   i_part              Current partition
  !! \param [out]  n_cell              Number of cells
  !! \param [out]  n_face              Number of faces
  !! \param [out]  n_face_part_bound   Number of partitioning boundary faces
  !! \param [out]  n_vtx               Number of vertices
  !! \param [out]  n_proc              Number of processus
  !! \param [out]  n_total_part        Number of partitions
  !! \param [out]  scell_face          Size of cell-face connectivity
  !! \param [out]  sface_vtx           Size of face-vertex connectivity
  !! \param [out]  sFacePartBound      Size of face_part_bound array
  !! \param [out]  sface_group         Size of face_group array
  !! \param [out]  n_face_group        Number of face groups
  !!

  subroutine PDM_part_part_dim_get_c (ppart, &
                                      i_part, &
                                      n_cell, &
                                      n_face, &
                                      n_face_part_bound, &
                                      n_vtx, &
                                      n_proc, &
                                      n_total_part, &
                                      scell_face, &
                                      sface_vtx, &
                                      sface_group, &
                                      n_face_group) &
  bind (c, name="PDM_part_part_dim_get")
    use pdm
    use iso_c_binding

    implicit none

    type(c_ptr),    value :: ppart
    integer(c_int), value :: i_part
    integer(c_int)        :: n_cell
    integer(c_int)        :: n_face
    integer(c_int)        :: n_face_part_bound
    integer(c_int)        :: n_vtx
    integer(c_int)        :: n_proc
    integer(c_int)        :: n_total_part
    integer(c_int)        :: scell_face
    integer(c_int)        :: sface_vtx
    integer(c_int)        :: sface_group
    integer(c_int)        :: n_face_group

  end subroutine PDM_part_part_dim_get_c


  !>
  !!
  !! \brief Return a mesh partition
  !!
  !! \param [in]   ppart                     Pointer to \ref PDM_part object
  !! \param [in]   i_part                    Current partition
  !! \param [out]  cell_tag                  Cell tag (size = n_cell)
  !! \param [out]  cell_face_idx             Cell to face connectivity index (size = n_cell + 1, numbering : 0 to n-1)
  !! \param [out]  cell_face                 Cell to face connectivity (size = cell_face_idx[n_cell] = lcell_face
  !!                                                                   numbering : 1 to n)
  !! \param [out]  cell_ln_to_gn             Cell local numbering to global numbering (size = n_cell, numbering : 1 to n)
  !! \param [out]  face_tag                  Face tag (size = n_face)
  !! \param [out]  face_cell                 Face to cell connectivity  (size = 2 * n_face, numbering : 1 to n)
  !! \param [out]  face_vtx_idx              Face to Vertex connectivity index (size = n_face + 1, numbering : 0 to n-1)
  !! \param [out]  face_vtx                  Face to Vertex connectivity (size = faceVertexIdx[n_face], numbering : 1 to n)
  !! \param [out]  face_ln_to_gn             Face local numbering to global numbering (size = n_face, numbering : 1 to n)
  !! \param [out]  face_part_bound_proc_idx  Partitioning boundary faces block distribution from processus (size = n_proc + 1)
  !! \param [out]  face_part_bound_part_idx  Partitioning boundary faces block distribution from partition (size = n_total_part + 1)
  !! \param [out]  face_part_bound           Partitioning boundary faces (size = 4 * n_face_part_bound)
  !!                                          sorted by processus, sorted by partition in each processus, and
  !!                                          sorted by absolute face number in each partition
  !!                                         For each face :
  !!                                           - Face local number (numbering : 1 to n)
  !!                                           - Connected process (numbering : 0 to n-1)
  !!                                           - Connected Partition
  !!                                             on the connected process (numbering :1 to n)
  !!                                           - Connected face local number
  !!                                             in the connected partition (numbering :1 to n)
  !! \param [out]  vtx_tag                   Vertex tag (size = nVertex)
  !! \param [out]  vtx                       Vertex coordinates (size = 3 * nVertex)
  !! \param [out]  vtx_ln_to_gn              Vertex local numbering to global numbering (size = n_vtx, numbering : 1 to n)
  !! \param [out]  face_group_idx            Face group index (size = n_face_group + 1, numbering : 1 to n-1)
  !! \param [out]  face_group                Faces for each group (size = face_group_idx[n_face_group] = lFaceGroup, numbering : 1 to n)
  !! \param [out]  face_group_ln_to_gn       Faces global numbering for each group
  !!                                         (size = face_group_idx[n_face_group] = lFaceGroup, numbering : 1 to n)
  !!

  subroutine PDM_part_part_val_get_c (ppart, &
                                      i_part, &
                                      cell_tag, &
                                      cell_face_idx, &
                                      cell_face, &
                                      cell_ln_to_gn, &
                                      face_tag, &
                                      face_cell, &
                                      face_vtx_idx, &
                                      face_vtx, &
                                      face_ln_to_gn, &
                                      face_part_bound_proc_idx, &
                                      face_part_bound_part_idx, &
                                      face_part_bound, &
                                      vtx_tag, &
                                      vtx, &
                                      vtx_ln_to_gn, &
                                      face_group_idx, &
                                      face_group, &
                                      face_group_ln_to_gn) &
  bind (c, name='PDM_part_part_val_get')

    use pdm
    use iso_c_binding

    implicit none

    type(c_ptr),    value :: ppart
    integer(c_int), value :: i_part
    type(c_ptr)           :: cell_tag
    type(c_ptr)           :: cell_face_idx
    type(c_ptr)           :: cell_face
    type(c_ptr)           :: cell_ln_to_gn
    type(c_ptr)           :: face_tag
    type(c_ptr)           :: face_cell
    type(c_ptr)           :: face_vtx_idx
    type(c_ptr)           :: face_vtx
    type(c_ptr)           :: face_ln_to_gn
    type(c_ptr)           :: face_part_bound_proc_idx
    type(c_ptr)           :: face_part_bound_part_idx
    type(c_ptr)           :: face_part_bound
    type(c_ptr)           :: vtx_tag
    type(c_ptr)           :: vtx
    type(c_ptr)           :: vtx_ln_to_gn
    type(c_ptr)           :: face_group_idx
    type(c_ptr)           :: face_group
    type(c_ptr)           :: face_group_ln_to_gn

  end subroutine PDM_part_part_val_get_c


  !>
  !!
  !! \brief Free ppart
  !!
  !! \param [in]   ppart               Pointer to \ref PDM_part object
  !!
  !!

  subroutine PDM_part_free (ppart) &
  bind (c, name='PDM_part_free')

    use iso_c_binding
    implicit none

    type(c_ptr), value :: ppart

  end subroutine PDM_part_free


  !>
  !!
  !! \brief Return times
  !!
  !! \param [in]   ppart       Pointer to \ref PDM_part object
  !! \param [out]  elapsed     elapsed times (size = 4)
  !! \param [out]  cpu         cpu times (size = 4)
  !! \param [out]  cpu_user    user cpu times (size = 4)
  !! \param [out]  cpu_sys     system cpu times (size = 4)
  !!
  !!

  subroutine PDM_part_time_get_c (ppart,    &
                                  elapsed,  &
                                  cpu,      &
                                  cpu_user, &
                                  cpu_sys)  &
  bind (c, name='PDM_part_time_get')

    use iso_c_binding
    implicit none

    type(c_ptr), value :: ppart
    type(c_ptr)        :: elapsed
    type(c_ptr)        :: cpu
    type(c_ptr)        :: cpu_user
    type(c_ptr)        :: cpu_sys

  end subroutine PDM_part_time_get_c


  !>
  !!
  !! \brief Return statistics
  !!
  !! \param [in]   ppart                          Pointer to \ref PDM_part object
  !! \param [out]  cells_average                  average of cells number
  !! \param [out]  cells_median                   median of cells number
  !! \param [out]  cells_std_deviation            standard deviation of cells number
  !! \param [out]  cells_min                      minimum of cells nummber
  !! \param [out]  cells_max                      maximum of cells nummber
  !! \param [out]  bound_part_faces_average       average of partitioning boundary faces
  !! \param [out]  bound_part_faces_median        median of partitioning boundary faces
  !! \param [out]  bound_part_faces_std_deviation standard deviation of partitioning boundary faces
  !! \param [out]  bound_part_faces_min           minimum of partitioning boundary faces
  !! \param [out]  bound_part_faces_max           maximum of partitioning boundary faces
  !!
  !!

  subroutine PDM_part_stat_get_c (ppart,                          &
                                  cells_average,                  &
                                  cells_median,                   &
                                  cells_std_deviation,            &
                                  cells_min,                      &
                                  cells_max,                      &
                                  bound_part_faces_average,       &
                                  bound_part_faces_median,        &
                                  bound_part_faces_std_deviation, &
                                  bound_part_faces_min,           &
                                  bound_part_faces_max,           &
                                  bound_part_faces_sum)           &
  bind (c, name='PDM_part_stat_get')

    use pdm
    use iso_c_binding
    implicit none

    type(c_ptr), value :: ppart
    integer(c_int)     :: cells_average
    integer(c_int)     :: cells_median
    real(c_double)     :: cells_std_deviation
    integer(c_int)     :: cells_min
    integer(c_int)     :: cells_max
    integer(c_int)     :: bound_part_faces_average
    integer(c_int)     :: bound_part_faces_median
    real(c_double)     :: bound_part_faces_std_deviation
    integer(c_int)     :: bound_part_faces_min
    integer(c_int)     :: bound_part_faces_max
    integer(c_int)     :: bound_part_faces_sum

  end subroutine PDM_part_stat_get_c



  !>
  !!
  !! \brief Return the coloring of a mesh partition
  !!
  !! \param [in]   ppart               Pointer to \ref PDM_part object
  !! \param [in]   i_part              Current partition
  !! \param [out]  cell_color          Cell color (size = n_cell)
  !! \param [out]  face_color          Face color (size = n_face)
  !! \param [out]  thread_color        Thread color (size = n_cell)
  !! \param [out]  hyperplane_color    Hyperplane color (size = n_cell)
  !!
  !!

  subroutine PDM_part_part_color_get_c (ppart, &
                                        i_part, &
                                        cell_color, &
                                        face_color, &
                                        thread_color, &
                                        hyperplane_color) &
  bind (c, name='PDM_part_part_color_get')

    use iso_c_binding
    implicit none

    type(c_ptr), value    :: ppart
    integer(c_int),value  :: i_part
    type(c_ptr)           :: cell_color
    type(c_ptr)           :: face_color
    type(c_ptr)           :: thread_color
    type(c_ptr)           :: hyperplane_color

  end subroutine PDM_part_part_color_get_c




  !>
  !!
  !! \brief Return an initialized \ref PDM_coarse_mesh object
  !!
  !! \param [in]   comm                   MPI Communicator
  !! \param [in]   method                 Split method
  !! \param [in]   renum_cell_method      Cell renumbering method
  !! \param [in]   renum_face_method      Face renumbering method
  !! \param [in]   n_property_cell        Number of cell properties
  !! \param [in]   renum_properties_cell  For cache blocking [ n_cell_per_cache_wanted, isAsynchrone, isVectorisation ] \ref PDM_renum_cacheblocking
  !! \param [in]   n_property_face        Number of face properties
  !! \param [in]   renum_properties_face  NOT USED?
  !! \param [in]   n_part                 Number of partitions
  !! \param [in]   n_total_part           Total number of partitions
  !! \param [in]   have_cell_tag          Presence of an array of cell tags
  !! \param [in]   have_face_tag          Presence of an array of face tags
  !! \param [in]   have_vtx_tag           Presence of an array of vertex tags
  !! \param [in]   have_cell_weight       Presence of an array of cell weights
  !! \param [in]   have_face_weight       Presence of an array of face weights
  !! \param [in]   have_face_group        Presence of an array of faces groups
  !!
  !! \return       Pointer to \ref PDM_coarse_mesh
  !!
  !!

  function PDM_part_coarse_mesh_create_c (comm,                  &
                                          method,                &
                                          renum_cell_method,     &
                                          renum_face_method,     &
                                          n_property_cell,       &
                                          renum_properties_cell, &
                                          n_property_face,       &
                                          renum_properties_face, &
                                          n_part,                &
                                          n_total_part,          &
                                          n_face_group,          &
                                          have_cell_tag,         &
                                          have_face_tag,         &
                                          have_vtx_tag,          &
                                          have_cell_weight,      &
                                          have_face_weight,      &
                                          have_face_group)       &
  result (cm) &
  bind (c, name='PDM_part_coarse_mesh_create')

    use iso_c_binding
    implicit none

    type(c_ptr)           :: cm
    integer(c_int), value :: comm
    character(c_char)     :: method(*)
    character(c_char)     :: renum_cell_method(*)
    character(c_char)     :: renum_face_method(*)
    integer(c_int), value :: n_property_cell
    type(c_ptr),    value :: renum_properties_cell
    integer(c_int), value :: n_property_face
    type(c_ptr),    value :: renum_properties_face
    integer(c_int), value :: n_part
    integer(c_int), value :: n_total_part
    integer(c_int), value :: n_face_group
    integer(c_int), value :: have_cell_tag
    integer(c_int), value :: have_face_tag
    integer(c_int), value :: have_vtx_tag
    integer(c_int), value :: have_cell_weight
    integer(c_int), value :: have_face_weight
    integer(c_int), value :: have_face_group

  end function PDM_part_coarse_mesh_create_c


  !>
  !!
  !! \brief Build a coarse mesh
  !!
  !! \param [in]  cm                        Pointer to \ref PDM_coarse_mesh
  !! \param [in]  i_part                    Partition identifier
  !! \param [in]  n_coarse_cell             Number of cells in the coarse grid
  !! \param [in]  n_cell                    Number of cells
  !! \param [in]  n_face                    Number of faces
  !! \param [in]  n_vtx                     Number of vertices
  !! \param [in]  n_face_group              Number of face groups
  !! \param [in]  n_face_part_bound         Number of partitioning boundary faces
  !! \param [in]  cell_face_idx             Cell to face connectivity index (size = n_cell + 1, numbering : 0 to n-1)
  !! \param [in]  cell_face                 Cell to face connectivity (size = cell_face_idx[n_cell] = lcell_face, numbering : 1 to n)
  !! \param [in]  cell_tag                  Cell tag (size = n_cell)
  !! \param [in]  cell_weight               Cell weight (size = n_cell)
  !! \param [in]  face_weight               Face weight (size = n_face)
  !! \param [in]  cell_ln_to_gn             Cell local numbering to global numbering (size = n_cell, numbering : 1 to n)
  !! \param [in]  face_cell                 Face to cell connectivity  (size = 2 * n_face, numbering : 1 to n)
  !! \param [in]  face_vtx_idx              Face to Vertex connectivity index (size = n_face + 1, numbering : 0 to n-1)
  !! \param [in]  face_vtx                  Face to Vertex connectivity (size = faceVertexIdx[n_face], numbering : 1 to n)
  !! \param [in]  face_tag                  Face tag (size = n_face)
  !! \param [in]  face_ln_to_gn             Face local numbering to global numbering (size = n_face, numbering : 1 to n)
  !! \param [in]  vtxCoord                  Vertex coordinates (size = 3 * nVertex)
  !! \param [in]  vtx_tag                   Vertex tag (size = nVertex)
  !! \param [in]  vtx_ln_to_gn              Vertex local numbering to global numbering (size = n_vtx, numbering : 1 to n)
  !! \param [in]  face_group_idx            Face group index (size = n_face_group + 1, numbering : 1 to n-1)
  !! \param [in]  face_group                Faces for each group (size = face_group_idx[n_face_group] = lFaceGroup, numbering : 1 to n)
  !! \param [in]  face_group_ln_to_gn       Faces global numbering for each group (size = face_group_idx[n_face_group] = lFaceGroup, numbering : 1 to n)
  !! \param [in]  face_part_bound_proc_idx  Partitioning boundary faces block distribution from processus (size = n_proc + 1)
  !! \param [in]  face_part_bound_part_idx  Partitioning boundary faces block distribution from partition (size = n_total_part + 1)
  !! \param [in]  face_part_bound           Partitioning boundary faces (size = 4 * n_face_part_bound)
  !!                                        sorted by processus, sorted by partition in each processus, and
  !!                                        sorted by absolute face number in each partition
  !!                                        For each face :
  !!                                          - Face local number (numbering : 1 to n)
  !!                                          - Connected process (numbering : 0 to n-1)
  !!                                          - Connected Partition
  !!                                            on the connected process (numbering :1 to n)
  !!                                          - Connected face local number
  !!                                            in the connected partition (numbering :1 to n)
  !!

  subroutine PDM_part_coarse_mesh_input_c (cm,                       &
                                            i_part,                   &
                                            n_coarse_cell,            &
                                            n_cell,                   &
                                            n_face,                   &
                                            n_vtx,                    &
                                            n_face_group,             &
                                            n_face_part_bound,        &
                                            cell_face_idx,            &
                                            cell_face,                &
                                            cell_tag,                 &
                                            cell_weight,              &
                                            face_weight,              &
                                            cell_ln_to_gn,            &
                                            face_cell,                &
                                            face_vtx_idx,             &
                                            face_vtx,                 &
                                            face_tag,                 &
                                            face_ln_to_gn,            &
                                            vtxCoord,                 &
                                            vtx_tag,                  &
                                            vtx_ln_to_gn,             &
                                            face_group_idx,           &
                                            face_group,               &
                                            face_group_ln_to_gn,      &
                                            face_part_bound_proc_idx, &
                                            face_part_bound_part_idx, &
                                            face_part_bound)          &
  bind (c, name='PDM_part_coarse_mesh_input')

    use iso_c_binding
    implicit none

    type(c_ptr),    value :: cm
    integer(c_int), value :: i_part
    integer(c_int), value :: n_coarse_cell
    integer(c_int), value :: n_cell
    integer(c_int), value :: n_face
    integer(c_int), value :: n_vtx
    integer(c_int), value :: n_face_group
    integer(c_int), value :: n_face_part_bound
    type(c_ptr),    value :: cell_face_idx
    type(c_ptr),    value :: cell_face
    type(c_ptr),    value :: cell_tag
    type(c_ptr),    value :: cell_weight
    type(c_ptr),    value :: face_weight
    type(c_ptr),    value :: cell_ln_to_gn
    type(c_ptr),    value :: face_cell
    type(c_ptr),    value :: face_vtx_idx
    type(c_ptr),    value :: face_vtx
    type(c_ptr),    value :: face_tag
    type(c_ptr),    value :: face_ln_to_gn
    type(c_ptr),    value :: vtxCoord
    type(c_ptr),    value :: vtx_tag
    type(c_ptr),    value :: vtx_ln_to_gn
    type(c_ptr),    value :: face_group_idx
    type(c_ptr),    value :: face_group
    type(c_ptr),    value :: face_group_ln_to_gn
    type(c_ptr),    value :: face_part_bound_proc_idx
    type(c_ptr),    value :: face_part_bound_part_idx
    type(c_ptr),    value :: face_part_bound

  end subroutine PDM_part_coarse_mesh_input_c


  !>
  !!
  !! \brief Updates all the arrays dealing with MPI exchanges
  !!
  !! \param [in] cm            Pointer to \ref PDM_coarse_mesh
  !!
  !!

  subroutine PDM_part_coarse_mesh_compute (cm) &
  bind (c, name='PDM_part_coarse_mesh_compute')

    use iso_c_binding
    implicit none

    type (c_ptr), value :: cm

  end subroutine PDM_part_coarse_mesh_compute


  !>
  !!
  !! \brief Return a coarse mesh partition dimensions
  !!
  !! \param [in]   cm                     Pointer to \ref PDM_coarse_mesh
  !! \param [in]   i_part                 Current partition
  !! \param [out]  n_cell                 Number of cells
  !! \param [out]  n_face                 Number of faces
  !! \param [out]  n_face_part_bound      Number of partitioning boundary faces
  !! \param [out]  n_vtx                  Number of vertices
  !! \param [out]  n_proc                 Number of processus
  !! \param [out]  n_total_part           Number of partitions
  !! \param [out]  n_face_group           Number of face groups
  !! \param [out]  scell_face             Size of cell-face connectivity
  !! \param [out]  sface_vtx              Size of face-vertex connectivity
  !! \param [out]  sface_group            Size of face_group array
  !! \param [out]  sCoarseCellToFineCell  Size of coarseCellToFineCell array
  !!

  subroutine PDM_part_coarse_mesh_part_dim_get_c (cm,                    &
                                                   i_part,                &
                                                   n_cell,                &
                                                   n_face,                &
                                                   n_face_part_bound,     &
                                                   n_vtx,                 &
                                                   n_proc,                &
                                                   n_total_part,          &
                                                   n_face_group,          &
                                                   scell_face,            &
                                                   sface_vtx,             &
                                                   sface_group,           &
                                                   sCoarseCellToFineCell) &
  bind(c, name='PDM_part_coarse_mesh_part_dim_get')

    use iso_c_binding
    implicit none

    type (c_ptr),   value :: cm
    integer(c_int), value :: i_part
    integer(c_int)        :: n_cell
    integer(c_int)        :: n_face
    integer(c_int)        :: n_face_part_bound
    integer(c_int)        :: n_vtx
    integer(c_int)        :: n_proc
    integer(c_int)        :: n_total_part
    integer(c_int)        :: n_face_group
    integer(c_int)        :: scell_face
    integer(c_int)        :: sface_vtx
    integer(c_int)        :: sface_group
    integer(c_int)        :: sCoarseCellToFineCell

  end subroutine PDM_part_coarse_mesh_part_dim_get_c


  !>
  !! \brief Return a mesh partition
  !!
  !! \param [in]   cm                        Pointer to \ref PDM_coarse_mesh
  !! \param [in]   i_part                    Current partition
  !! \param [out]  cell_face_idx             Cell to face connectivity index (size = n_cell + 1, numbering : 0 to n-1)
  !! \param [out]  cell_face                 Cell to face connectivity (size = cell_face_idx[n_cell] = lcell_face, numbering : 1 to n)
  !! \param [out]  cell_tag                  Cell tag (size = n_cell)
  !! \param [out]  cell_ln_to_gn             Cell local numbering to global numbering (size = n_cell, numbering : 1 to n)
  !! \param [out]  cellInitCellIdx           Array of indexes of the connected partitions (size : n_coarse_cell + 1)
  !! \param [out]  cellInitCell              Partitioning array (size : cellInitCellIdx[n_coarse_cell])
  !! \param [out]  face_cell                 Face to cell connectivity  (size = 2 * n_face, numbering : 1 to n)
  !! \param [out]  face_vtx_idx              Face to Vertex connectivity index (size = n_face + 1, numbering : 0 to n-1)
  !! \param [out]  face_vtx                  Face to Vertex connectivity (size = faceVertexIdx[n_face], numbering : 1 to n)
  !! \param [out]  face_tag                  Face tag (size = n_face)
  !! \param [out]  face_ln_to_gn             Face local numbering to global numbering (size = n_face, numbering : 1 to n)
  !! \param [out]  faceGroupInitFaceGroup    Coarse face group - fine face group connectivity (size = face_group_idx[n_face_group])
  !! \param [out]  faceInitFace              Coarse face - fine face connectivity (size = nCoarseFace)
  !! \param [out]  vtxCoord                  Vertex coordinates (size = 3 * n_vtx)
  !! \param [out]  vtx_tag                   Vertex tag (size = n_vtx)
  !! \param [out]  vtx_ln_to_gn              Vertex local numbering to global numbering (size = n_vtx, numbering : 1 to n)
  !! \param [out]  vtxInitVtx                Coarse vertex - fine vertex connectivity (size = nCoarseVtx)
  !! \param [out]  face_group_idx            Face group index (size = n_face_group + 1, numbering : 1 to n-1)
  !! \param [out]  face_group                Faces for each group (size = face_group_idx[n_face_group] = lFaceGroup, numbering : 1 to n)
  !! \param [out]  face_group_ln_to_gn       Faces global numbering for each group (size = face_group_idx[n_face_group] = lFaceGroup, numbering : 1 to n)
  !! \param [out]  face_part_bound_proc_idx  Partitioning boundary faces block distribution from processus (size = n_proc + 1)
  !! \param [out]  face_part_bound_part_idx  Partitioning boundary faces block distribution from partition (size = n_total_part + 1)
  !! \param [out]  face_part_bound           Partitioning boundary faces (size = 4 * n_face_part_bound)
  !!                                         sorted by processus, sorted by partition in each processus, and
  !!                                         sorted by absolute face number in each partition
  !!                                         For each face :
  !!                                           - Face local number (numbering : 1 to n)
  !!                                           - Connected process (numbering : 0 to n-1)
  !!                                           - Connected Partition
  !!                                             on the connected process (numbering :1 to n)
  !!                                           - Connected face local number
  !!                                             in the connected partition (numbering :1 to n)
  !!

  subroutine PDM_part_coarse_mesh_part_get_c (cm,                       &
                                               i_part,                   &
                                               cell_face_idx,            &
                                               cell_face,                &
                                               cell_tag,                 &
                                               cell_ln_to_gn,            &
                                               cellInitCellIdx,          &
                                               cellInitCell,             &
                                               face_cell,                &
                                               face_vtx_idx,             &
                                               face_vtx,                 &
                                               face_tag,                 &
                                               face_ln_to_gn,            &
                                               faceGroupInitFaceGroup,   &
                                               faceInitFace,             &
                                               vtxCoord,                 &
                                               vtx_tag,                  &
                                               vtx_ln_to_gn,             &
                                               vtxInitVtx,               &
                                               face_group_idx,           &
                                               face_group,               &
                                               face_group_ln_to_gn,      &
                                               face_part_bound_proc_idx, &
                                               face_part_bound_part_idx, &
                                               face_part_bound)          &
  bind (c, name='PDM_part_coarse_mesh_part_get')

    use iso_c_binding
    implicit none

    type(c_ptr),    value :: cm
    integer(c_int), value :: i_part
    type(c_ptr)           :: cell_face_idx
    type(c_ptr)           :: cell_face
    type(c_ptr)           :: cell_tag
    type(c_ptr)           :: cell_ln_to_gn
    type(c_ptr)           :: cellInitCellIdx
    type(c_ptr)           :: cellInitCell
    type(c_ptr)           :: face_cell
    type(c_ptr)           :: face_vtx_idx
    type(c_ptr)           :: face_vtx
    type(c_ptr)           :: face_tag
    type(c_ptr)           :: face_ln_to_gn
    type(c_ptr)           :: faceGroupInitFaceGroup
    type(c_ptr)           :: faceInitFace
    type(c_ptr)           :: vtxCoord
    type(c_ptr)           :: vtx_tag
    type(c_ptr)           :: vtx_ln_to_gn
    type(c_ptr)           :: vtxInitVtx
    type(c_ptr)           :: face_group_idx
    type(c_ptr)           :: face_group
    type(c_ptr)           :: face_group_ln_to_gn
    type(c_ptr)           :: face_part_bound_proc_idx
    type(c_ptr)           :: face_part_bound_part_idx
    type(c_ptr)           :: face_part_bound

  end subroutine PDM_part_coarse_mesh_part_get_c


  !>
  !!
  !! \brief Return the coloring of a coarse mesh
  !!
  !! \param [in]   cm                  Pointer to \ref PDM_coarse_mesh
  !! \param [in]   i_part              Current partition
  !! \param [out]  cell_color          Cell color (size = n_cell)
  !! \param [out]  face_color          Face color (size = n_face)
  !! \param [out]  thread_color        Thread color (size = n_cell)
  !! \param [out]  hyperplane_color    Hyperplane color (size = n_cell)
  !!
  !!

  subroutine PDM_part_coarse_color_get_c (cm, &
                                         i_part, &
                                         cell_color, &
                                         face_color, &
                                         thread_color, &
                                         hyperplane_color) &
  bind (c, name='PDM_part_coarse_color_get')

    use iso_c_binding
    implicit none

    type(c_ptr), value :: cm
    integer(c_int)     :: i_part
    type(c_ptr)        :: cell_color
    type(c_ptr)        :: face_color
    type(c_ptr)        :: thread_color
    type(c_ptr)        :: hyperplane_color

  end subroutine PDM_part_coarse_color_get_c



  !>
  !!
  !! \brief Free coarse mesh
  !!
  !! \param [in]   cm                  Pointer to \ref PDM_coarse_mesh
  !!
  !!

  subroutine PDM_part_coarse_mesh_free (cm) &
  bind (c, name='PDM_part_coarse_mesh_free')

    use iso_c_binding
    implicit none

    type(c_ptr), value :: cm

  end subroutine PDM_part_coarse_mesh_free



  !>
  !!
  !! \brief Return times
  !!
  !! \param [in]   cm          Pointer to \ref PDM_coarse_mesh
  !! \param [out]  elapsed     Elapsed times (size = 18)
  !! \param [out]  cpu         Cpu times (size = 18)
  !! \param [out]  cpu_user    User cpu times (size = 18)
  !! \param [out]  cpu_sys     System cpu times (size = 18)
  !!
  !!

  subroutine PDM_part_coarse_mesh_time_get_c (cm,       &
                                               elapsed,  &
                                               cpu,      &
                                               cpu_user, &
                                               cpu_sys)  &
  bind (c, name='PDM_part_coarse_mesh_time_get')

    use iso_c_binding
    implicit none

    type(c_ptr), value :: cm
    type(c_ptr)        :: elapsed
    type(c_ptr)        :: cpu
    type(c_ptr)        :: cpu_user
    type(c_ptr)        :: cpu_sys

  end subroutine PDM_part_coarse_mesh_time_get_c



  !>
  !!
  !! \brief Displays all the arrays of a coarse mesh
  !!
  !! \param [in]   cm         Pointer to \ref PDM_coarse_mesh
  !!
  !!

  subroutine PDM_part_coarse_mesh_display (cm) &
  bind (c, name='PDM_part_coarse_mesh_display')

    use iso_c_binding
    implicit none

    type (c_ptr), value :: cm

  end subroutine PDM_part_coarse_mesh_display


end interface

private :: PDM_part_create_ ,&
           PDM_part_part_dim_get_ ,&
           PDM_part_part_val_get_ ,&
           PDM_part_time_get_ ,&
           PDM_part_stat_get_ ,&
           PDM_part_part_color_get_ ,&
           PDM_part_coarse_mesh_create_ ,&
           PDM_part_coarse_mesh_input_ ,&
           PDM_part_coarse_mesh_part_dim_get_ ,&
           PDM_part_coarse_mesh_part_get_ ,&
           PDM_part_coarse_color_get_ ,&
           PDM_part_coarse_mesh_time_get_

contains

  !>
  !!
  !! \brief Build a initial partitioning
  !!  Build a initial partitioning from :
  !!      - Cell block distribution with implicit global numbering
  !!         (the first cell is the first cell of the first process and
  !!          the latest cell is the latest cell of the latest process)
  !!      - Face block distribution with implicit global numbering
  !!      - Vertex block distribution with implicit global numbering
  !!  To repart an existing partition use \ref PDM_part_repart function
  !!
  !! \param [out]  ppart                  Pointer to \ref PDM_part object
  !! \param [in]   comm                   MPI Comminicator
  !! \param [in]   split_method           Split method
  !! \param [in]   renum_cell_method      Cell renumbering method
  !! \param [in]   renum_face_method      Face renumbering method
  !! \param [in]   n_property_cell        Number of cell properties
  !! \param [in]   renum_properties_cell  For cache blocking [ n_cell_per_cache_wanted, isAsynchrone, isVectorisation ] \ref PDM_renum_cacheblocking
  !! \param [in]   n_property_face        Number of face properties
  !! \param [in]   renum_properties_face  NOT USED
  !! \param [in]   n_part                 Number of partition to build on this process
  !! \param [in]   dn_cell                Number of distributed cells
  !! \param [in]   dn_face                Number of distributed faces
  !! \param [in]   dn_vtx                 Number of distributed vertices
  !! \param [in]   n_face_group           Number of face groups
  !! \param [in]   dcell_face_idx         Distributed cell face connectivity index or NULL
  !!                                      (size : dn_cell + 1, numbering : 0 to n-1)
  !! \param [in]   dcell_face             Distributed cell face connectivity or NULL
  !!                                      (size : dface_vtx_idx[dn_cell], numbering : 1 to n)
  !! \param [in]   dcell_tag              Cell tag (size : n_cell) or NULL
  !! \param [in]   dcell_weight           Cell weight (size : n_cell) or NULL
  !! \param [in]   have_dcell_part        Presence of an array of cell part id
  !! \param [in]   dcell_part             Distributed cell partitioning
  !!                                      (size = dn_cell) or NULL (No partitioning if != NULL)
  !! \param [in]   dface_cell             Distributed face cell connectivity or NULL
  !!                                      (size : 2 * dn_face, numbering : 1 to n)
  !! \param [in]   dface_vtx_idx          Distributed face to vertex connectivity index
  !!                                      (size : dn_face + 1, numbering : 0 to n-1)
  !! \param [in]   dface_vtx              Distributed face to vertex connectivity
  !!                                      (size : dface_vtx_idx[dn_face], numbering : 1 to n)
  !! \param [in]   dface_tag              Distributed face tag (size : dn_face)
  !!                                      or NULL
  !! \param [in]   dvtx_coord             Distributed vertex coordinates
  !!                                      (size : 3*dn_vtx)
  !! \param [in]   dvtx_tag               Distributed vertex tag (size : dn_vtx) or NULL
  !! \param [in]   dface_group_idx        Index of distributed faces list of each group
  !!                                      (size = n_face_group + 1) or NULL
  !! \param [in]   dface_group            Distributed faces list of each group
  !!                                      (size = dface_group[dface_group_idx[n_face_group]], numbering : 1 to n)
  !!                                      or NULL
  !!
  !!

  subroutine PDM_part_create_(ppart, &
                              f_comm, &
                              split_method,  &
                              renum_cell_method, &
                              renum_face_method, &
                              nPropertyCell, &
                              renum_properties_cell, &
                              nPropertyFace, &
                              renum_properties_face, &
                              nPart, &
                              dNCell, &
                              dNFace, &
                              dNVtx,&
                              nFaceGroup, &
                              dCellFaceIdx,&
                              dCellFace,&
                              dCellTag, &
                              dCellWeight,&
                              have_dCellPart,&
                              dCellPart,&
                              dFaceCell,&
                              dFaceVtxIdx,&
                              dFaceVtx,&
                              dFaceTag,&
                              dVtxCoord,&
                              dVtxTag,&
                              dFaceGroupIdx,&
                              dFaceGroup)

    use pdm
    use iso_c_binding

    implicit none

    type(c_ptr)                        :: ppart
    integer, intent(in)                :: f_comm
    integer, intent(in)                :: split_method
    character (len=*)                  :: renum_cell_method
    character (len=*)                  :: renum_face_method
    integer, intent(in)                :: nPropertyCell
    integer(kind=PDM_l_num_s), pointer :: renum_properties_cell(:)
    integer, intent(in)                :: nPropertyFace
    integer(kind=PDM_l_num_s), pointer :: renum_properties_face(:)
    integer, intent(in)                :: nPart
    integer, intent(in)                :: dNCell
    integer, intent(in)                :: dNFace
    integer, intent(in)                :: dNVtx
    integer, intent(in)                :: nFaceGroup
    integer(kind=PDM_l_num_s), pointer :: dCellFaceIdx(:)
    integer(kind=PDM_g_num_s), pointer :: dCellFace(:)
    integer(kind=PDM_l_num_s), pointer :: dCellTag(:)
    integer(kind=PDM_l_num_s), pointer :: dCellWeight(:)
    logical, intent(in)                :: have_dCellPart
    integer(kind=PDM_l_num_s), pointer :: dCellPart(:)
    integer(kind=PDM_g_num_s), pointer :: dFaceCell(:)
    integer(kind=PDM_l_num_s), pointer :: dFaceVtxIdx(:)
    integer(kind=PDM_g_num_s), pointer :: dFaceVtx(:)
    integer(kind=PDM_l_num_s), pointer :: dFaceTag(:)
    double precision,          pointer :: dVtxCoord(:,:)
    integer(kind=PDM_l_num_s), pointer :: dVtxTag(:)
    integer(kind=PDM_l_num_s), pointer :: dFaceGroupIdx(:)
    integer(kind=PDM_g_num_s), pointer :: dFaceGroup(:)

    integer(kind=c_int)                :: c_split_method
    integer(kind=c_int)                :: c_n_property_cell
    type(c_ptr)                        :: c_renum_properties_cell
    integer(kind=c_int)                :: c_n_property_face
    type(c_ptr)                        :: c_renum_properties_face
    integer(kind=c_int)                :: c_n_part
    integer(kind=c_int)                :: c_dn_cell
    integer(kind=c_int)                :: c_dn_face
    integer(kind=c_int)                :: c_dn_vtx
    integer(kind=c_int)                :: c_n_face_group
    type(c_ptr)                        :: c_dcell_faceIdx
    type(c_ptr)                        :: c_dcell_face
    type(c_ptr)                        :: c_dcell_tag
    type(c_ptr)                        :: c_dcell_weight
    integer(c_int)                     :: c_have_dcell_part
    type(c_ptr)                        :: c_dcell_part
    type(c_ptr)                        :: c_dface_cell
    type(c_ptr)                        :: c_dface_vtx_idx
    type(c_ptr)                        :: c_dface_vtx
    type(c_ptr)                        :: c_dface_tag
    type(c_ptr)                        :: c_dvtx_coord
    type(c_ptr)                        :: c_dvtx_tag
    type(c_ptr)                        :: c_dface_group_idx
    type(c_ptr)                        :: c_dface_group
    integer(c_int)                     :: c_comm


    if (have_dCellPart .and. associated(dCellPart)) then
      c_have_dcell_part = 1
    else
      c_have_dcell_part = 0
    endif

    c_dcell_tag       = C_NULL_PTR
    if (associated(dCellTag)) then
      c_dcell_tag = c_loc(dCellTag)
    endif

    c_dcell_part = C_NULL_PTR
    if (associated(dCellPart)) then
      c_dcell_part = c_loc(dCellPart)
    endif

    c_dface_cell = C_NULL_PTR
    if (associated(dFaceCell)) then
      c_dface_cell = c_loc(dFaceCell)
    endif

    c_dface_tag = C_NULL_PTR
    if (associated(dFaceTag)) then
      c_dface_tag = c_loc(dFaceTag)
    endif

    c_dvtx_tag = C_NULL_PTR
    if (associated(dVtxTag)) then
      c_dvtx_tag = c_loc(dVtxTag)
    endif

    c_dcell_weight = C_NULL_PTR
    if (associated(dCellWeight)) then
      c_dcell_weight = c_loc(dCellWeight)
    endif

    c_renum_properties_face = C_NULL_PTR
    if (associated(renum_properties_face)) then
      c_renum_properties_face = c_loc(renum_properties_face)
    endif

    c_dcell_faceIdx = C_NULL_PTR
    if (associated(dCellFaceIdx)) then
      c_dcell_faceIdx = c_loc(dCellFaceIdx)
    endif

    c_dcell_face = C_NULL_PTR
    if (associated(dCellFace)) then
      c_dcell_face = c_loc(dCellFace)
    endif

    c_split_method    = split_method
    c_n_property_cell = nPropertyCell
    c_n_property_face = nPropertyFace
    c_n_part          = nPart
    c_dn_cell         = dnCell
    c_dn_face         = dnFace
    c_dn_vtx          = dnVtx
    c_n_face_group    = nFaceGroup

    c_renum_properties_cell = C_NULL_PTR
    if (associated (renum_properties_cell)) then 
      c_renum_properties_cell = c_loc(renum_properties_cell)
    endif 

    c_dface_vtx_idx         = C_NULL_PTR
    if (associated (dFaceVtxIdx)) then 
      c_dface_vtx_idx         = c_loc(dFaceVtxIdx          )
    endif 

    c_dface_vtx             = C_NULL_PTR
    if (associated (dFaceVtx)) then 
      c_dface_vtx             = c_loc(dFaceVtx             )
    endif 

    c_dvtx_coord            = C_NULL_PTR
    if (associated (dVtxCoord)) then 
      c_dvtx_coord            = c_loc(dVtxCoord            )
    endif 

    c_dface_group_idx       = C_NULL_PTR
    if (associated (dFaceGroupIdx)) then 
      c_dface_group_idx       = c_loc(dFaceGroupIdx        )
    endif 

    c_dface_group           = C_NULL_PTR
    if (associated (dFaceGroup)) then 
      c_dface_group           = c_loc(dFaceGroup           )
    endif 

    c_comm = PDM_MPI_Comm_f2c(f_comm)

    ppart = PDM_part_create_c (c_comm, &
                                c_split_method, &
                                trim(renum_cell_method)//C_NULL_CHAR, &
                                trim(renum_face_method)//C_NULL_CHAR, &
                                c_n_property_cell, &
                                c_renum_properties_cell, &
                                c_n_property_face, &
                                c_renum_properties_face, &
                                c_n_part, &
                                c_dn_cell, &
                                c_dn_face, &
                                c_dn_vtx, &
                                c_n_face_group, &
                                c_dcell_faceIdx, &
                                c_dcell_face, &
                                c_dcell_tag, &
                                c_dcell_weight, &
                                c_have_dcell_part, &
                                c_dcell_part, &
                                c_dface_cell, &
                                c_dface_vtx_idx, &
                                c_dface_vtx, &
                                c_dface_tag, &
                                c_dvtx_coord, &
                                c_dvtx_tag, &
                                c_dface_group_idx, &
                                c_dface_group)

  end subroutine PDM_part_create_



  !>
  !!
  !! \brief Return a mesh partition dimensions
  !!
  !! \param [in]   ppart               Pointer to \ref PDM_part object
  !! \param [in]   i_part              Current partition
  !! \param [out]  n_cell              Number of cells
  !! \param [out]  n_face              Number of faces
  !! \param [out]  n_face_part_bound   Number of partitioning boundary faces
  !! \param [out]  n_vtx               Number of vertices
  !! \param [out]  n_proc              Number of processus
  !! \param [out]  n_total_part        Number of partitions
  !! \param [out]  scell_face          Size of cell-face connectivity
  !! \param [out]  sface_vtx           Size of face-vertex connectivity
  !! \param [out]  sFacePartBound      Size of face_part_bound array
  !! \param [out]  sface_group         Size of face_group array
  !! \param [out]  n_face_group        Number of face groups
  !!

  subroutine PDM_part_part_dim_get_ (ppart, &
                                     i_part, &
                                     n_cell, &
                                     n_face, &
                                     n_face_part_bound, &
                                     n_vtx, &
                                     n_proc, &
                                     n_total_part, &
                                     scell_face, &
                                     sface_vtx, &
                                     sface_group, &
                                     n_face_group)
    use pdm
    use iso_c_binding

    implicit none

    type(c_ptr), value   :: ppart
    integer, intent(in)  :: i_part
    integer, intent(out) :: n_cell
    integer, intent(out) :: n_face
    integer, intent(out) :: n_face_part_bound
    integer, intent(out) :: n_vtx
    integer, intent(out) :: n_proc
    integer, intent(out) :: n_total_part
    integer, intent(out) :: scell_face
    integer, intent(out) :: sface_vtx
    integer, intent(out) :: sface_group
    integer, intent(out) :: n_face_group

    integer(c_int)       :: c_i_part
    integer(c_int)       :: c_n_cell
    integer(c_int)       :: c_n_face
    integer(c_int)       :: c_n_face_part_bound
    integer(c_int)       :: c_n_vtx
    integer(c_int)       :: c_n_proc
    integer(c_int)       :: c_n_total_part
    integer(c_int)       :: c_scell_face
    integer(c_int)       :: c_sface_vtx
    integer(c_int)       :: c_sface_group
    integer(c_int)       :: c_n_face_group

    c_i_part = i_part

    call PDM_part_part_dim_get_c(ppart, &
                                 c_i_part, &
                                 c_n_cell, &
                                 c_n_face, &
                                 c_n_face_part_bound, &
                                 c_n_vtx, &
                                 c_n_proc, &
                                 c_n_total_part, &
                                 c_scell_face, &
                                 c_sface_vtx, &
                                 c_sface_group, &
                                 c_n_face_group)

    n_cell            = c_n_cell
    n_face            = c_n_face
    n_face_part_bound = c_n_face_part_bound
    n_vtx             = c_n_vtx
    n_proc            = c_n_proc
    n_total_part      = c_n_total_part
    scell_face        = c_scell_face
    sface_vtx         = c_sface_vtx
    sface_group       = c_sface_group
    n_face_group      = c_n_face_group

  end subroutine PDM_part_part_dim_get_


  !>
  !!
  !! \brief Return a mesh partition
  !!
  !! \param [in]   ppart                     Pointer to \ref PDM_part object
  !! \param [in]   i_part                    Current partition
  !! \param [out]  cell_tag                  Cell tag (size = n_cell)
  !! \param [out]  cell_face_idx             Cell to face connectivity index (size = n_cell + 1, numbering : 0 to n-1)
  !! \param [out]  cell_face                 Cell to face connectivity (size = cell_face_idx[n_cell] = lcell_face
  !!                                                                   numbering : 1 to n)
  !! \param [out]  cell_ln_to_gn             Cell local numbering to global numbering (size = n_cell, numbering : 1 to n)
  !! \param [out]  face_tag                  Face tag (size = n_face)
  !! \param [out]  face_cell                 Face to cell connectivity  (size = 2 * n_face, numbering : 1 to n)
  !! \param [out]  face_vtx_idx              Face to Vertex connectivity index (size = n_face + 1, numbering : 0 to n-1)
  !! \param [out]  face_vtx                  Face to Vertex connectivity (size = faceVertexIdx[n_face], numbering : 1 to n)
  !! \param [out]  face_ln_to_gn             Face local numbering to global numbering (size = n_face, numbering : 1 to n)
  !! \param [out]  face_part_bound_proc_idx  Partitioning boundary faces block distribution from processus (size = n_proc + 1)
  !! \param [out]  face_part_bound_part_idx  Partitioning boundary faces block distribution from partition (size = n_total_part + 1)
  !! \param [out]  face_part_bound           Partitioning boundary faces (size = 4 * n_face_part_bound)
  !!                                          sorted by processus, sorted by partition in each processus, and
  !!                                          sorted by absolute face number in each partition
  !!                                         For each face :
  !!                                           - Face local number (numbering : 1 to n)
  !!                                           - Connected process (numbering : 0 to n-1)
  !!                                           - Connected Partition
  !!                                             on the connected process (numbering :1 to n)
  !!                                           - Connected face local number
  !!                                             in the connected partition (numbering :1 to n)
  !! \param [out]  vtx_tag                   Vertex tag (size = nVertex)
  !! \param [out]  vtx                       Vertex coordinates (size = 3 * nVertex)
  !! \param [out]  vtx_ln_to_gn              Vertex local numbering to global numbering (size = n_vtx, numbering : 1 to n)
  !! \param [out]  face_group_idx            Face group index (size = n_face_group + 1, numbering : 1 to n-1)
  !! \param [out]  face_group                Faces for each group (size = face_group_idx[n_face_group] = lFaceGroup, numbering : 1 to n)
  !! \param [out]  face_group_ln_to_gn       Faces global numbering for each group
  !!                                         (size = face_group_idx[n_face_group] = lFaceGroup, numbering : 1 to n)
  !!

  subroutine PDM_part_part_val_get_ (ppart,                    &
                                     i_part,                   &
                                     cell_tag,                 &
                                     cell_face_idx,            &
                                     cell_face,                &
                                     cell_ln_to_gn,            &
                                     face_tag,                 &
                                     face_cell,                &
                                     face_vtx_idx,             &
                                     face_vtx,                 &
                                     face_ln_to_gn,            &
                                     face_part_bound_proc_idx, &
                                     face_part_bound_part_idx, &
                                     face_part_bound,          &
                                     vtx_tag,                  &
                                     vtx,                      &
                                     vtx_ln_to_gn,             &
                                     face_group_idx,           &
                                     face_group,               &
                                     face_group_ln_to_gn)

    use pdm
    use iso_c_binding

    implicit none

    type (c_ptr), value                   :: ppart
    integer, intent(in)                   :: i_part
    integer (kind = PDM_l_num_s), pointer :: cell_tag(:)
    integer (kind = PDM_l_num_s), pointer :: cell_face_idx(:)
    integer (kind = PDM_l_num_s), pointer :: cell_face(:)
    integer (kind = PDM_g_num_s), pointer :: cell_ln_to_gn(:)
    integer (kind = PDM_l_num_s), pointer :: face_tag(:)
    integer (kind = PDM_l_num_s), pointer :: face_cell(:)
    integer (kind = PDM_l_num_s), pointer :: face_vtx_idx(:)
    integer (kind = PDM_l_num_s), pointer :: face_vtx(:)
    integer (kind = PDM_g_num_s), pointer :: face_ln_to_gn(:)
    integer (kind = PDM_l_num_s), pointer :: face_part_bound_proc_idx(:)
    integer (kind = PDM_l_num_s), pointer :: face_part_bound_part_idx(:)
    integer (kind = PDM_l_num_s), pointer :: face_part_bound(:)
    integer (kind = PDM_l_num_s), pointer :: vtx_tag(:)
    double precision,             pointer :: vtx(:,:)
    integer (kind = PDM_g_num_s), pointer :: vtx_ln_to_gn(:)
    integer (kind = PDM_l_num_s), pointer :: face_group_idx(:)
    integer (kind = PDM_l_num_s), pointer :: face_group(:)
    integer (kind = PDM_g_num_s), pointer :: face_group_ln_to_gn(:)

    integer(c_int)                        :: c_i_part
    type(c_ptr)                           :: c_cell_tag
    type(c_ptr)                           :: c_cell_face_idx
    type(c_ptr)                           :: c_cell_face
    type(c_ptr)                           :: c_cell_ln_to_gn
    type(c_ptr)                           :: c_face_tag
    type(c_ptr)                           :: c_face_cell
    type(c_ptr)                           :: c_face_vtx_idx
    type(c_ptr)                           :: c_face_vtx
    type(c_ptr)                           :: c_face_ln_to_gn
    type(c_ptr)                           :: c_face_part_bound_proc_idx
    type(c_ptr)                           :: c_face_part_bound_part_idx
    type(c_ptr)                           :: c_face_part_bound
    type(c_ptr)                           :: c_vtx_tag
    type(c_ptr)                           :: c_vtx
    type(c_ptr)                           :: c_vtx_ln_to_gn
    type(c_ptr)                           :: c_face_group_idx
    type(c_ptr)                           :: c_face_group
    type(c_ptr)                           :: c_face_group_ln_to_gn
    integer(c_int)                        :: n_cell
    integer(c_int)                        :: n_face
    integer(c_int)                        :: n_face_part_bound
    integer(c_int)                        :: n_vtx
    integer(c_int)                        :: n_proc
    integer(c_int)                        :: n_total_part
    integer(c_int)                        :: scell_face
    integer(c_int)                        :: sface_vtx
    integer(c_int)                        :: sface_group
    integer(c_int)                        :: n_face_group

    c_i_part = i_part

    call PDM_part_part_dim_get_c(ppart,             &
                                 c_i_part,          &
                                 n_cell,            &
                                 n_face,            &
                                 n_face_part_bound, &
                                 n_vtx,             &
                                 n_proc,            &
                                 n_total_part,      &
                                 scell_face,        &
                                 sface_vtx,         &
                                 sface_group,       &
                                 n_face_group)

    c_cell_tag                 = C_NULL_PTR
    c_cell_face_idx            = C_NULL_PTR
    c_cell_face                = C_NULL_PTR
    c_cell_ln_to_gn            = C_NULL_PTR
    c_face_tag                 = C_NULL_PTR
    c_face_cell                = C_NULL_PTR
    c_face_vtx_idx             = C_NULL_PTR
    c_face_vtx                 = C_NULL_PTR
    c_face_ln_to_gn            = C_NULL_PTR
    c_face_part_bound_proc_idx = C_NULL_PTR
    c_face_part_bound_part_idx = C_NULL_PTR
    c_face_part_bound          = C_NULL_PTR
    c_vtx_tag                  = C_NULL_PTR
    c_vtx                      = C_NULL_PTR
    c_vtx_ln_to_gn             = C_NULL_PTR
    c_face_group_idx           = C_NULL_PTR
    c_face_group               = C_NULL_PTR
    c_face_group_ln_to_gn      = C_NULL_PTR
    call PDM_part_part_val_get_c(ppart,                     &
                                 c_i_part,                   &
                                 c_cell_tag,                 &
                                 c_cell_face_idx,            &
                                 c_cell_face,                &
                                 c_cell_ln_to_gn,            &
                                 c_face_tag,                 &
                                 c_face_cell,                &
                                 c_face_vtx_idx,             &
                                 c_face_vtx,                 &
                                 c_face_ln_to_gn,            &
                                 c_face_part_bound_proc_idx, &
                                 c_face_part_bound_part_idx, &
                                 c_face_part_bound,          &
                                 c_vtx_tag,                  &
                                 c_vtx,                      &
                                 c_vtx_ln_to_gn,             &
                                 c_face_group_idx,           &
                                 c_face_group,               &
                                 c_face_group_ln_to_gn)

    call c_f_pointer(c_cell_tag, &
                     cell_tag,   &
                     [n_cell])

    call c_f_pointer(c_cell_face_idx, &
                     cell_face_idx,   &
                     [n_cell+1])

    call c_f_pointer(c_cell_face, &
                     cell_face,   &
                     [scell_face])

    call c_f_pointer(c_cell_ln_to_gn, &
                     cell_ln_to_gn,   &
                     [n_cell])

    call c_f_pointer(c_face_tag, &
                     face_tag,   &
                     [n_face])

    call c_f_pointer(c_face_cell, &
                     face_cell,   &
                     [2*n_face])

    call c_f_pointer(c_face_vtx_idx, &
                     face_vtx_idx,   &
                     [n_face+1])

    call c_f_pointer(c_face_vtx, &
                     face_vtx,   &
                     [sface_vtx])

    call c_f_pointer(c_face_ln_to_gn, &
                     face_ln_to_gn,   &
                     [n_face])

    call c_f_pointer(c_face_part_bound_proc_idx, &
                     face_part_bound_proc_idx,   &
                     [n_proc+1])

    call c_f_pointer(c_face_part_bound_part_idx, &
                     face_part_bound_part_idx,   &
                     [n_total_part+1])

    call c_f_pointer(c_face_part_bound, &
                     face_part_bound,   &
                     [4*n_face_part_bound])

    call c_f_pointer(c_vtx_tag, &
                     vtx_tag,   &
                     [n_vtx])

    call c_f_pointer(c_vtx, &
                     vtx,   &
                     [3, n_vtx])

    call c_f_pointer(c_vtx_ln_to_gn, &
                     vtx_ln_to_gn,   &
                     [n_vtx])

    call c_f_pointer(c_face_group_idx, &
                     face_group_idx,   &
                     [n_face_group+1])

    call c_f_pointer(c_face_group, &
                     face_group,   &
                     [sface_group])

    call c_f_pointer(c_face_group_ln_to_gn, &
                     face_group_ln_to_gn,   &
                     [sface_group])

  end subroutine PDM_part_part_val_get_



  !>
  !!
  !! \brief Return times
  !!
  !! \param [in]   ppart       Pointer to \ref PDM_part object
  !! \param [out]  elapsed     Elapsed times (size = 4)
  !! \param [out]  cpu         Cpu times (size = 4)
  !! \param [out]  cpu_user    User cpu times (size = 4)
  !! \param [out]  cpu_sys     System cpu times (size = 4)
  !!
  !!

  subroutine PDM_part_time_get_ (ppart,    &
                                 elapsed,  &
                                 cpu,      &
                                 cpu_user, &
                                 cpu_sys)

    use iso_c_binding
    implicit none

    type(c_ptr), value            :: ppart
    double precision, intent(out) :: elapsed(4)
    double precision, intent(out) :: cpu(4)
    double precision, intent(out) :: cpu_user(4)
    double precision, intent(out) :: cpu_sys(4)

    type(c_ptr)                   :: c_elapsed
    type(c_ptr)                   :: c_cpu
    type(c_ptr)                   :: c_cpu_user
    type(c_ptr)                   :: c_cpu_sys
    double precision, pointer     :: ptr(:)


    c_elapsed  = C_NULL_PTR
    c_cpu      = C_NULL_PTR
    c_cpu_user = C_NULL_PTR
    c_cpu_sys  = C_NULL_PTR
    ptr => null()
    
    call PDM_part_time_get_c (ppart,      &
                              c_elapsed,  &
                              c_cpu,      &
                              c_cpu_user, &
                              c_cpu_sys)

    call c_f_pointer(c_elapsed, &
                     ptr,       &
                     [4])
    elapsed(1:4) = ptr(1:4)

    call c_f_pointer(c_cpu, &
                     ptr,   &
                     [4])
    cpu(1:4) = ptr(1:4)

    call c_f_pointer(c_cpu_user, &
                     ptr,        &
                     [4])
    cpu_user(1:4) = ptr(1:4)

    call c_f_pointer(c_cpu_sys, &
                     ptr,       &
                     [4])
    cpu_sys(1:4) = ptr(1:4)

  end subroutine PDM_part_time_get_



  !>
  !!
  !! \brief Return statistics
  !!
  !! \param [in]   ppart                          Pointer to \ref PDM_part object
  !! \param [out]  cells_average                  average of cells number
  !! \param [out]  cells_median                   median of cells number
  !! \param [out]  cells_std_deviation            standard deviation of cells number
  !! \param [out]  cells_min                      minimum of cells nummber
  !! \param [out]  cells_max                      maximum of cells nummber
  !! \param [out]  bound_part_faces_average       average of partitioning boundary faces
  !! \param [out]  bound_part_faces_median        median of partitioning boundary faces
  !! \param [out]  bound_part_faces_std_deviation standard deviation of partitioning boundary faces
  !! \param [out]  bound_part_faces_min           minimum of partitioning boundary faces
  !! \param [out]  bound_part_faces_max           maximum of partitioning boundary faces
  !!
  !!

  subroutine PDM_part_stat_get_ (ppart,                          &
                                 cells_average,                  &
                                 cells_median,                   &
                                 cells_std_deviation,            &
                                 cells_min,                      &
                                 cells_max,                      &
                                 bound_part_faces_average,       &
                                 bound_part_faces_median,        &
                                 bound_part_faces_std_deviation, &
                                 bound_part_faces_min,           &
                                 bound_part_faces_max,           &
                                 bound_part_faces_sum)

    use pdm
    use iso_c_binding
    implicit none

    type(c_ptr), value                :: ppart
    integer(pdm_l_num_s), intent(out) :: cells_average
    integer(pdm_l_num_s), intent(out) :: cells_median
    double precision,     intent(out) :: cells_std_deviation
    integer(pdm_l_num_s), intent(out) :: cells_min
    integer(pdm_l_num_s), intent(out) :: cells_max
    integer(pdm_l_num_s), intent(out) :: bound_part_faces_average
    integer(pdm_l_num_s), intent(out) :: bound_part_faces_median
    double precision,     intent(out) :: bound_part_faces_std_deviation
    integer(pdm_l_num_s), intent(out) :: bound_part_faces_min
    integer(pdm_l_num_s), intent(out) :: bound_part_faces_max
    integer(pdm_l_num_s), intent(out) :: bound_part_faces_sum

    integer(c_int)                    :: c_cells_average
    integer(c_int)                    :: c_cells_median
    real(c_double)                    :: c_cells_std_deviation
    integer(c_int)                    :: c_cells_min
    integer(c_int)                    :: c_cells_max
    integer(c_int)                    :: c_bound_part_faces_average
    integer(c_int)                    :: c_bound_part_faces_median
    real(c_double)                    :: c_bound_part_faces_std_deviation
    integer(c_int)                    :: c_bound_part_faces_min
    integer(c_int)                    :: c_bound_part_faces_max
    integer(c_int)                    :: c_bound_part_faces_sum


    call PDM_part_stat_get_c(ppart,                            &
                             c_cells_average,                  &
                             c_cells_median,                   &
                             c_cells_std_deviation,            &
                             c_cells_min,                      &
                             c_cells_max,                      &
                             c_bound_part_faces_average,       &
                             c_bound_part_faces_median,        &
                             c_bound_part_faces_std_deviation, &
                             c_bound_part_faces_min,           &
                             c_bound_part_faces_max,           &
                             c_bound_part_faces_sum)

    cells_average                  = c_cells_average
    cells_median                   = c_cells_median
    cells_std_deviation            = c_cells_std_deviation
    cells_min                      = c_cells_min
    cells_max                      = c_cells_max
    bound_part_faces_average       = c_bound_part_faces_average
    bound_part_faces_median        = c_bound_part_faces_median
    bound_part_faces_std_deviation = c_bound_part_faces_std_deviation
    bound_part_faces_min           = c_bound_part_faces_min
    bound_part_faces_max           = c_bound_part_faces_max
    bound_part_faces_sum           = c_bound_part_faces_sum

  end subroutine PDM_part_stat_get_



  !>
  !!
  !! \brief Return the coloring of a mesh partition
  !!
  !! \param [in]   ppart               Pointer to \ref PDM_part object
  !! \param [in]   i_part              Current partition
  !! \param [out]  cell_color          Cell color (size = n_cell)
  !! \param [out]  face_color          Face color (size = n_face)
  !! \param [out]  thread_color        Thread color (size = n_cell)
  !! \param [out]  hyperplane_color    Hyperplane color (size = n_cell)
  !!
  !!

  subroutine PDM_part_part_color_get_ (ppart, &
                                       i_part, &
                                       cell_color, &
                                       face_color, &
                                       thread_color, &
                                       hyperplane_color)

    use iso_c_binding
    implicit none

    type(c_ptr), value                 :: ppart
    integer, intent(in)                :: i_part
    integer(kind=PDM_l_num_s), pointer :: cell_color(:)
    integer(kind=PDM_l_num_s), pointer :: face_color(:)
    integer(kind=PDM_l_num_s), pointer :: thread_color(:)
    integer(kind=PDM_l_num_s), pointer :: hyperplane_color(:)

    integer(c_int)                     :: c_i_part
    type(c_ptr)                        :: c_cell_color
    type(c_ptr)                        :: c_face_color
    type(c_ptr)                        :: c_thread_color
    type(c_ptr)                        :: c_hyperplane_color
    integer(c_int)                     :: n_cell
    integer(c_int)                     :: n_face
    integer(c_int)                     :: n_face_part_bound
    integer(c_int)                     :: n_vtx
    integer(c_int)                     :: n_proc
    integer(c_int)                     :: n_total_part
    integer(c_int)                     :: scell_face
    integer(c_int)                     :: sface_vtx
    integer(c_int)                     :: sface_group
    integer(c_int)                     :: n_face_group


    c_i_part = i_part

    call PDM_part_part_dim_get_c(ppart,             &
                                 c_i_part,          &
                                 n_cell,            &
                                 n_face,            &
                                 n_face_part_bound, &
                                 n_vtx,             &
                                 n_proc,            &
                                 n_total_part,      &
                                 scell_face,        &
                                 sface_vtx,         &
                                 sface_group,       &
                                 n_face_group)

    call PDM_part_part_color_get_c(ppart,              &
                                   c_i_part,           &
                                   c_cell_color,       &
                                   c_face_color,       &
                                   c_thread_color,     &
                                   c_hyperplane_color)

    call c_f_pointer(c_cell_color, &
                     cell_color,   &
                     [n_cell])

    call c_f_pointer(c_face_color, &
                     face_color,   &
                     [n_face])

    call c_f_pointer(c_thread_color, &
                     thread_color,   &
                     [n_cell])

    call c_f_pointer(c_hyperplane_color, &
                     hyperplane_color,   &
                     [n_cell])

  end subroutine PDM_part_part_color_get_


  !>
  !!
  !! \brief Return an initialized \ref PDM_coarse_mesh object
  !!
  !! \param [out]  cm                     Pointer to \ref PDM_coarse_mesh
  !! \param [in]   f_comm                 MPI Communicator
  !! \param [in]   method                 Split method
  !! \param [in]   renum_cell_method      Cell renumbering method
  !! \param [in]   renum_face_method      Face renumbering method
  !! \param [in]   n_property_cell        Number of cell properties
  !! \param [in]   renum_properties_cell  For cache blocking [ n_cell_per_cache_wanted, isAsynchrone, isVectorisation ] \ref PDM_renum_cacheblocking
  !! \param [in]   n_property_face        Number of face properties
  !! \param [in]   renum_properties_face  NOT USED?
  !! \param [in]   n_part                 Number of partitions
  !! \param [in]   n_total_part           Total number of partitions
  !! \param [in]   have_cell_tag          Presence of an array of cell tags
  !! \param [in]   have_face_tag          Presence of an array of face tags
  !! \param [in]   have_vtx_tag           Presence of an array of vertex tags
  !! \param [in]   have_cell_weight       Presence of an array of cell weights
  !! \param [in]   have_face_weight       Presence of an array of face weights
  !! \param [in]   have_face_group        Presence of an array of faces groups
  !!
  !!

  subroutine PDM_part_coarse_mesh_create_ (cm,                    &
                                           f_comm,                &
                                           method,                &
                                           renum_cell_method,     &
                                           renum_face_method,     &
                                           n_property_cell,       &
                                           renum_properties_cell, &
                                           n_property_face,       &
                                           renum_properties_face, &
                                           n_part,                &
                                           n_total_part,          &
                                           n_face_group,          &
                                           have_cell_tag,         &
                                           have_face_tag,         &
                                           have_vtx_tag,          &
                                           have_cell_weight,      &
                                           have_face_weight,      &
                                           have_face_group)

    use iso_c_binding
    implicit none

    type(c_ptr)                        :: cm
    integer, intent(in)                :: f_comm
    character(len = *)                 :: method
    character(len = *)                 :: renum_cell_method
    character(len = *)                 :: renum_face_method
    integer, intent(in)                :: n_property_cell
    integer(kind=PDM_l_num_s), pointer :: renum_properties_cell(:)
    integer, intent(in)                :: n_property_face
    integer(kind=PDM_l_num_s), pointer :: renum_properties_face(:)
    integer, intent(in)                :: n_part
    integer, intent(in)                :: n_total_part
    integer, intent(in)                :: n_face_group
    logical, intent(in)                :: have_cell_tag
    logical, intent(in)                :: have_face_tag
    logical, intent(in)                :: have_vtx_tag
    logical, intent(in)                :: have_cell_weight
    logical, intent(in)                :: have_face_weight
    logical, intent(in)                :: have_face_group

    integer(c_int)                     :: c_comm
    integer(c_int)                     :: c_n_property_cell
    type(c_ptr)                        :: c_renum_properties_cell
    integer(c_int)                     :: c_n_property_face
    type(c_ptr)                        :: c_renum_properties_face
    integer(c_int)                     :: c_n_part
    integer(c_int)                     :: c_n_total_part
    integer(c_int)                     :: c_n_face_group
    integer(c_int)                     :: c_have_cell_tag
    integer(c_int)                     :: c_have_face_tag
    integer(c_int)                     :: c_have_vtx_tag
    integer(c_int)                     :: c_have_cell_weight
    integer(c_int)                     :: c_have_face_weight
    integer(c_int)                     :: c_have_face_group


    c_comm = PDM_MPI_Comm_f2c(f_comm)

    c_n_property_cell = n_property_cell
    c_n_property_face = n_property_face
    c_n_part          = n_part
    c_n_total_part    = n_total_part
    c_n_face_group    = n_face_group

    c_renum_properties_cell = C_NULL_PTR
    if (n_property_cell > 0 .and. associated(renum_properties_cell)) then
      c_renum_properties_cell = c_loc(renum_properties_cell)
    endif

    c_renum_properties_face = C_NULL_PTR
    if (n_property_face > 0 .and. associated(renum_properties_face)) then
      c_renum_properties_face = c_loc(renum_properties_face)
    endif
    
    if (have_cell_tag) then
      c_have_cell_tag = 1
    else
      c_have_cell_tag = 0
    endif

    if (have_face_tag) then
      c_have_face_tag = 1
    else
      c_have_face_tag = 0
    endif

    if (have_vtx_tag) then
      c_have_vtx_tag = 1
    else
      c_have_vtx_tag = 0
    endif

    if (have_cell_weight) then
      c_have_cell_weight = 1
    else
      c_have_cell_weight = 0
    endif

    if (have_face_weight) then
      c_have_face_weight = 1
    else
      c_have_face_weight = 0
    endif

    if (have_face_group) then
      c_have_face_group = 1
    else
      c_have_face_group = 0
    endif

    cm = PDM_part_coarse_mesh_create_c (c_comm,                               &
                                        trim(method)//C_NULL_CHAR,             &
                                        trim(renum_cell_method)//C_NULL_CHAR, &
                                        trim(renum_face_method)//C_NULL_CHAR, &
                                        c_n_property_cell,                    &
                                        c_renum_properties_cell,              &
                                        c_n_property_face,                    &
                                        c_renum_properties_face,              &
                                        c_n_part,                             &
                                        c_n_total_part,                       &
                                        c_n_face_group,                       &
                                        c_have_cell_tag,                      &
                                        c_have_face_tag,                      &
                                        c_have_vtx_tag,                       &
                                        c_have_cell_weight,                   &
                                        c_have_face_weight,                   &
                                        c_have_face_group)

  end subroutine PDM_part_coarse_mesh_create_



  !>
  !!
  !! \brief Build a coarse mesh
  !!
  !! \param [in]  cm                        Pointer to \ref PDM_coarse_mesh
  !! \param [in]  i_part                    Partition identifier
  !! \param [in]  n_coarse_cell             Number of cells in the coarse grid
  !! \param [in]  n_cell                    Number of cells
  !! \param [in]  n_face                    Number of faces
  !! \param [in]  n_vtx                     Number of vertices
  !! \param [in]  n_face_group              Number of face groups
  !! \param [in]  n_face_part_bound         Number of partitioning boundary faces
  !! \param [in]  cell_face_idx             Cell to face connectivity index (size = n_cell + 1, numbering : 0 to n-1)
  !! \param [in]  cell_face                 Cell to face connectivity (size = cell_face_idx[n_cell] = lcell_face, numbering : 1 to n)
  !! \param [in]  cell_tag                  Cell tag (size = n_cell)
  !! \param [in]  cell_weight               Cell weight (size = n_cell)
  !! \param [in]  face_weight               Face weight (size = n_face)
  !! \param [in]  cell_ln_to_gn             Cell local numbering to global numbering (size = n_cell, numbering : 1 to n)
  !! \param [in]  face_cell                 Face to cell connectivity  (size = 2 * n_face, numbering : 1 to n)
  !! \param [in]  face_vtx_idx              Face to Vertex connectivity index (size = n_face + 1, numbering : 0 to n-1)
  !! \param [in]  face_vtx                  Face to Vertex connectivity (size = faceVertexIdx[n_face], numbering : 1 to n)
  !! \param [in]  face_tag                  Face tag (size = n_face)
  !! \param [in]  face_ln_to_gn             Face local numbering to global numbering (size = n_face, numbering : 1 to n)
  !! \param [in]  vtxCoord                  Vertex coordinates (size = 3 * nVertex)
  !! \param [in]  vtx_tag                   Vertex tag (size = nVertex)
  !! \param [in]  vtx_ln_to_gn              Vertex local numbering to global numbering (size = n_vtx, numbering : 1 to n)
  !! \param [in]  face_group_idx            Face group index (size = n_face_group + 1, numbering : 1 to n-1)
  !! \param [in]  face_group                Faces for each group (size = face_group_idx[n_face_group] = lFaceGroup, numbering : 1 to n)
  !! \param [in]  face_group_ln_to_gn       Faces global numbering for each group (size = face_group_idx[n_face_group] = lFaceGroup, numbering : 1 to n)
  !! \param [in]  face_part_bound_proc_idx  Partitioning boundary faces block distribution from processus (size = n_proc + 1)
  !! \param [in]  face_part_bound_part_idx  Partitioning boundary faces block distribution from partition (size = n_total_part + 1)
  !! \param [in]  face_part_bound           Partitioning boundary faces (size = 4 * n_face_part_bound)
  !!                                        sorted by processus, sorted by partition in each processus, and
  !!                                        sorted by absolute face number in each partition
  !!                                        For each face :
  !!                                          - Face local number (numbering : 1 to n)
  !!                                          - Connected process (numbering : 0 to n-1)
  !!                                          - Connected Partition
  !!                                            on the connected process (numbering :1 to n)
  !!                                          - Connected face local number
  !!                                            in the connected partition (numbering :1 to n)
  !!

  subroutine PDM_part_coarse_mesh_input_ (cm,                       &
                                          i_part,                   &
                                          n_coarse_cell,            &
                                          n_cell,                   &
                                          n_face,                   &
                                          n_vtx,                    &
                                          n_face_group,             &
                                          n_face_part_bound,        &
                                          cell_face_idx,            &
                                          cell_face,                &
                                          cell_tag,                 &
                                          cell_weight,              &
                                          face_weight,              &
                                          cell_ln_to_gn,            &
                                          face_cell,                &
                                          face_vtx_idx,             &
                                          face_vtx,                 &
                                          face_tag,                 &
                                          face_ln_to_gn,            &
                                          vtxCoord,                 &
                                          vtx_tag,                  &
                                          vtx_ln_to_gn,             &
                                          face_group_idx,           &
                                          face_group,               &
                                          face_group_ln_to_gn,      &
                                          face_part_bound_proc_idx, &
                                          face_part_bound_part_idx, &
                                          face_part_bound)

    use iso_c_binding
    implicit none

    type(c_ptr), value                 :: cm
    integer, intent(in)                :: i_part
    integer, intent(in)                :: n_coarse_cell
    integer, intent(in)                :: n_cell
    integer, intent(in)                :: n_face
    integer, intent(in)                :: n_vtx
    integer, intent(in)                :: n_face_group
    integer, intent(in)                :: n_face_part_bound
    integer(kind=PDM_l_num_s), pointer :: cell_face_idx(:)
    integer(kind=PDM_l_num_s), pointer :: cell_face(:)
    integer(kind=PDM_l_num_s), pointer :: cell_tag(:)
    integer(kind=PDM_l_num_s), pointer :: cell_weight(:)
    integer(kind=PDM_l_num_s), pointer :: face_weight(:)
    integer(kind=PDM_g_num_s), pointer :: cell_ln_to_gn(:)
    integer(kind=PDM_l_num_s), pointer :: face_cell(:)
    integer(kind=PDM_l_num_s), pointer :: face_vtx_idx(:)
    integer(kind=PDM_l_num_s), pointer :: face_vtx(:)
    integer(kind=PDM_l_num_s), pointer :: face_tag(:)
    integer(kind=PDM_g_num_s), pointer :: face_ln_to_gn(:)
    double precision,          pointer :: vtxCoord(:,:)
    integer(kind=PDM_l_num_s), pointer :: vtx_tag(:)
    integer(kind=PDM_g_num_s), pointer :: vtx_ln_to_gn(:)
    integer(kind=PDM_l_num_s), pointer :: face_group_idx(:)
    integer(kind=PDM_l_num_s), pointer :: face_group(:)
    integer(kind=PDM_g_num_s), pointer :: face_group_ln_to_gn(:)
    integer(kind=PDM_l_num_s), pointer :: face_part_bound_proc_idx(:)
    integer(kind=PDM_l_num_s), pointer :: face_part_bound_part_idx(:)
    integer(kind=PDM_l_num_s), pointer :: face_part_bound(:)

    integer(c_int)                     :: c_i_part
    integer(c_int)                     :: c_n_coarse_cell
    integer(c_int)                     :: c_n_cell
    integer(c_int)                     :: c_n_face
    integer(c_int)                     :: c_n_vtx
    integer(c_int)                     :: c_n_face_group
    integer(c_int)                     :: c_n_face_part_bound
    type(c_ptr)                        :: c_cell_face_idx
    type(c_ptr)                        :: c_cell_face
    type(c_ptr)                        :: c_cell_tag
    type(c_ptr)                        :: c_cell_weight
    type(c_ptr)                        :: c_face_weight
    type(c_ptr)                        :: c_cell_ln_to_gn
    type(c_ptr)                        :: c_face_cell
    type(c_ptr)                        :: c_face_vtx_idx
    type(c_ptr)                        :: c_face_vtx
    type(c_ptr)                        :: c_face_tag
    type(c_ptr)                        :: c_face_ln_to_gn
    type(c_ptr)                        :: c_vtxCoord
    type(c_ptr)                        :: c_vtx_tag
    type(c_ptr)                        :: c_vtx_ln_to_gn
    type(c_ptr)                        :: c_face_group_idx
    type(c_ptr)                        :: c_face_group
    type(c_ptr)                        :: c_face_group_ln_to_gn
    type(c_ptr)                        :: c_face_part_bound_proc_idx
    type(c_ptr)                        :: c_face_part_bound_part_idx
    type(c_ptr)                        :: c_face_part_bound


    c_i_part            = i_part
    c_n_coarse_cell     = n_coarse_cell
    c_n_cell            = n_cell
    c_n_face            = n_face
    c_n_vtx             = n_vtx
    c_n_face_group      = n_face_group
    c_n_face_part_bound = n_face_part_bound

    c_cell_face_idx = C_NULL_PTR
    if (associated (cell_face_idx)) then
      c_cell_face_idx            = c_loc(cell_face_idx           )
    endif

    c_cell_face = C_NULL_PTR
    if (associated (cell_face)) then
      c_cell_face                = c_loc(cell_face               )
    endif

    c_cell_ln_to_gn = C_NULL_PTR
    if (associated (cell_ln_to_gn)) then
      c_cell_ln_to_gn            = c_loc(cell_ln_to_gn           )
    endif

    c_face_cell = C_NULL_PTR
    if (associated (face_cell)) then
      c_face_cell                = c_loc(face_cell               )
    endif

    c_face_vtx_idx = C_NULL_PTR
    if (associated (face_vtx_idx)) then
      c_face_vtx_idx             = c_loc(face_vtx_idx            )
    endif

    c_face_vtx = C_NULL_PTR
    if (associated (face_vtx)) then
      c_face_vtx                 = c_loc(face_vtx                )
    endif

    c_face_ln_to_gn = C_NULL_PTR
    if (associated (face_ln_to_gn)) then
      c_face_ln_to_gn            = c_loc(face_ln_to_gn           )
    endif

    c_vtxCoord = C_NULL_PTR
    if (associated (vtxCoord)) then
      c_vtxCoord                 = c_loc(vtxCoord                )
    endif

    c_vtx_ln_to_gn = C_NULL_PTR
    if (associated (vtx_ln_to_gn)) then
      c_vtx_ln_to_gn             = c_loc(vtx_ln_to_gn            )
    endif

    c_face_part_bound_proc_idx = C_NULL_PTR
    if (associated (face_part_bound_proc_idx)) then
      c_face_part_bound_proc_idx = c_loc(face_part_bound_proc_idx)
    endif

    c_face_part_bound_part_idx = C_NULL_PTR
    if (associated (face_part_bound_part_idx)) then
      c_face_part_bound_part_idx = c_loc(face_part_bound_part_idx)
    endif  

    c_vtx_tag = C_NULL_PTR
    if (associated(vtx_tag)) then
      c_vtx_tag  = c_loc(vtx_tag)
    endif

    c_face_tag = C_NULL_PTR
    if (associated(face_tag)) then
      c_face_tag = c_loc(face_tag)
    endif

    c_cell_tag = C_NULL_PTR
    if (associated(cell_tag)) then
      c_cell_tag = c_loc(cell_tag)
    endif

    c_cell_weight = C_NULL_PTR
    if (associated(cell_weight)) then
      c_cell_weight = c_loc(cell_weight)
    endif

    c_face_weight = C_NULL_PTR
    if (associated(face_weight)) then
      c_face_weight = c_loc(face_weight)
    endif

    c_face_group_idx = C_NULL_PTR
    if (associated(face_group_idx)) then
      c_face_group_idx = c_loc(face_group_idx)
    endif

    c_face_group = C_NULL_PTR
    if (associated(face_group)) then
      c_face_group = c_loc(face_group)
    endif
    
    c_face_group_ln_to_gn = C_NULL_PTR
    if (associated(face_group_ln_to_gn)) then
      c_face_group_ln_to_gn = c_loc(face_group_ln_to_gn)
    endif
    
    c_face_part_bound = C_NULL_PTR
    if (associated(face_part_bound)) then
      c_face_part_bound = c_loc(face_part_bound)
    endif

    call PDM_part_coarse_mesh_input_c (cm,                         &
                                       c_i_part,                   &
                                       c_n_coarse_cell,            &
                                       c_n_cell,                   &
                                       c_n_face,                   &
                                       c_n_vtx,                    &
                                       c_n_face_group,             &
                                       c_n_face_part_bound,        &
                                       c_cell_face_idx,            &
                                       c_cell_face,                &
                                       c_cell_tag,                 &
                                       c_cell_weight,              &
                                       c_face_weight,              &
                                       c_cell_ln_to_gn,            &
                                       c_face_cell,                &
                                       c_face_vtx_idx,             &
                                       c_face_vtx,                 &
                                       c_face_tag,                 &
                                       c_face_ln_to_gn,            &
                                       c_vtxCoord,                 &
                                       c_vtx_tag,                  &
                                       c_vtx_ln_to_gn,             &
                                       c_face_group_idx,           &
                                       c_face_group,               &
                                       c_face_group_ln_to_gn,      &
                                       c_face_part_bound_proc_idx, &
                                       c_face_part_bound_part_idx, &
                                       c_face_part_bound)

  end subroutine PDM_part_coarse_mesh_input_



  !>
  !!
  !! \brief Return a coarse mesh partition dimensions
  !!
  !! \param [in]   cm                     Pointer to \ref PDM_coarse_mesh
  !! \param [in]   i_part                 Current partition
  !! \param [out]  n_cell                 Number of cells
  !! \param [out]  n_face                 Number of faces
  !! \param [out]  n_face_part_bound      Number of partitioning boundary faces
  !! \param [out]  n_vtx                  Number of vertices
  !! \param [out]  n_proc                 Number of processus
  !! \param [out]  n_total_part           Number of partitions
  !! \param [out]  n_face_group           Number of face groups
  !! \param [out]  scell_face             Size of cell-face connectivity
  !! \param [out]  sface_vtx              Size of face-vertex connectivity
  !! \param [out]  sface_group            Size of face_group array
  !! \param [out]  sCoarseCellToFineCell  Size of coarseCellToFineCell array
  !!

  subroutine PDM_part_coarse_mesh_part_dim_get_ (cm,                    &
                                                 i_part,                &
                                                 n_cell,                &
                                                 n_face,                &
                                                 n_face_part_bound,     &
                                                 n_vtx,                 &
                                                 n_proc,                &
                                                 n_total_part,          &
                                                 n_face_group,          &
                                                 scell_face,            &
                                                 sface_vtx,             &
                                                 sface_group,           &
                                                 sCoarseCellToFineCell)

    use iso_c_binding
    implicit none

    type (c_ptr), value  :: cm
    integer, intent(in)  :: i_part
    integer, intent(out) :: n_cell
    integer, intent(out) :: n_face
    integer, intent(out) :: n_face_part_bound
    integer, intent(out) :: n_vtx
    integer, intent(out) :: n_proc
    integer, intent(out) :: n_total_part
    integer, intent(out) :: n_face_group
    integer, intent(out) :: scell_face
    integer, intent(out) :: sface_vtx
    integer, intent(out) :: sface_group
    integer, intent(out) :: sCoarseCellToFineCell

    integer(c_int)       :: c_i_part
    integer(c_int)       :: c_n_cell
    integer(c_int)       :: c_n_face
    integer(c_int)       :: c_n_face_part_bound
    integer(c_int)       :: c_n_vtx
    integer(c_int)       :: c_n_proc
    integer(c_int)       :: c_n_total_part
    integer(c_int)       :: c_n_face_group
    integer(c_int)       :: c_scell_face
    integer(c_int)       :: c_sface_vtx
    integer(c_int)       :: c_sface_group
    integer(c_int)       :: c_sCoarseCellToFineCell

    c_i_part = i_part

    call PDM_part_coarse_mesh_part_dim_get_c (cm,                      &
                                              c_i_part,                &
                                              c_n_cell,                &
                                              c_n_face,                &
                                              c_n_face_part_bound,     &
                                              c_n_vtx,                 &
                                              c_n_proc,                &
                                              c_n_total_part,          &
                                              c_n_face_group,          &
                                              c_scell_face,            &
                                              c_sface_vtx,             &
                                              c_sface_group,           &
                                              c_sCoarseCellToFineCell)

    n_cell                = c_n_cell
    n_face                = c_n_face
    n_face_part_bound     = c_n_face_part_bound
    n_vtx                 = c_n_vtx
    n_proc                = c_n_proc
    n_total_part          = c_n_total_part
    n_face_group          = c_n_face_group
    scell_face            = c_scell_face
    sface_vtx             = c_sface_vtx
    sface_group           = c_sface_group
    sCoarseCellToFineCell = c_sCoarseCellToFineCell

  end subroutine PDM_part_coarse_mesh_part_dim_get_



  !>
  !! \brief Return a mesh partition
  !!
  !! \param [in]   cm                        Pointer to \ref PDM_coarse_mesh
  !! \param [in]   i_part                    Current partition
  !! \param [out]  cell_face_idx             Cell to face connectivity index (size = n_cell + 1, numbering : 0 to n-1)
  !! \param [out]  cell_face                 Cell to face connectivity (size = cell_face_idx[n_cell] = lcell_face, numbering : 1 to n)
  !! \param [out]  cell_tag                  Cell tag (size = n_cell)
  !! \param [out]  cell_ln_to_gn             Cell local numbering to global numbering (size = n_cell, numbering : 1 to n)
  !! \param [out]  cellInitCellIdx           Array of indexes of the connected partitions (size : n_coarse_cell + 1)
  !! \param [out]  cellInitCell              Partitioning array (size : cellInitCellIdx[n_coarse_cell])
  !! \param [out]  face_cell                 Face to cell connectivity  (size = 2 * n_face, numbering : 1 to n)
  !! \param [out]  face_vtx_idx              Face to Vertex connectivity index (size = n_face + 1, numbering : 0 to n-1)
  !! \param [out]  face_vtx                  Face to Vertex connectivity (size = faceVertexIdx[n_face], numbering : 1 to n)
  !! \param [out]  face_tag                  Face tag (size = n_face)
  !! \param [out]  face_ln_to_gn             Face local numbering to global numbering (size = n_face, numbering : 1 to n)
  !! \param [out]  faceGroupInitFaceGroup    Coarse face group - fine face group connectivity (size = face_group_idx[n_face_group])
  !! \param [out]  faceInitFace              Coarse face - fine face connectivity (size = nCoarseFace)
  !! \param [out]  vtxCoord                  Vertex coordinates (size = 3 * n_vtx)
  !! \param [out]  vtx_tag                   Vertex tag (size = n_vtx)
  !! \param [out]  vtx_ln_to_gn              Vertex local numbering to global numbering (size = n_vtx, numbering : 1 to n)
  !! \param [out]  vtxInitVtx                Coarse vertex - fine vertex connectivity (size = nCoarseVtx)
  !! \param [out]  face_group_idx            Face group index (size = n_face_group + 1, numbering : 1 to n-1)
  !! \param [out]  face_group                Faces for each group (size = face_group_idx[n_face_group] = lFaceGroup, numbering : 1 to n)
  !! \param [out]  face_group_ln_to_gn       Faces global numbering for each group (size = face_group_idx[n_face_group] = lFaceGroup, numbering : 1 to n)
  !! \param [out]  face_part_bound_proc_idx  Partitioning boundary faces block distribution from processus (size = n_proc + 1)
  !! \param [out]  face_part_bound_part_idx  Partitioning boundary faces block distribution from partition (size = n_total_part + 1)
  !! \param [out]  face_part_bound           Partitioning boundary faces (size = 4 * n_face_part_bound)
  !!                                         sorted by processus, sorted by partition in each processus, and
  !!                                         sorted by absolute face number in each partition
  !!                                         For each face :
  !!                                           - Face local number (numbering : 1 to n)
  !!                                           - Connected process (numbering : 0 to n-1)
  !!                                           - Connected Partition
  !!                                             on the connected process (numbering :1 to n)
  !!                                           - Connected face local number
  !!                                             in the connected partition (numbering :1 to n)
  !!

  subroutine PDM_part_coarse_mesh_part_get_ (cm,                       &
                                             i_part,                   &
                                             cell_face_idx,            &
                                             cell_face,                &
                                             cell_tag,                 &
                                             cell_ln_to_gn,            &
                                             cellInitCellIdx,          &
                                             cellInitCell,             &
                                             face_cell,                &
                                             face_vtx_idx,             &
                                             face_vtx,                 &
                                             face_tag,                 &
                                             face_ln_to_gn,            &
                                             faceGroupInitFaceGroup,   &
                                             faceInitFace,             &
                                             vtxCoord,                 &
                                             vtx_tag,                  &
                                             vtx_ln_to_gn,             &
                                             vtxInitVtx,               &
                                             face_group_idx,           &
                                             face_group,               &
                                             face_group_ln_to_gn,      &
                                             face_part_bound_proc_idx, &
                                             face_part_bound_part_idx, &
                                             face_part_bound)

    use iso_c_binding
    implicit none

    type(c_ptr), value                 :: cm
    integer, intent(in)                :: i_part
    integer(kind=PDM_l_num_s), pointer :: cell_face_idx(:)
    integer(kind=PDM_l_num_s), pointer :: cell_face(:)
    integer(kind=PDM_l_num_s), pointer :: cell_tag(:)
    integer(kind=PDM_g_num_s), pointer :: cell_ln_to_gn(:)
    integer(kind=PDM_l_num_s), pointer :: cellInitCellIdx(:)
    integer(kind=PDM_l_num_s), pointer :: cellInitCell(:)
    integer(kind=PDM_l_num_s), pointer :: face_cell(:)
    integer(kind=PDM_l_num_s), pointer :: face_vtx_idx(:)
    integer(kind=PDM_l_num_s), pointer :: face_vtx(:)
    integer(kind=PDM_l_num_s), pointer :: face_tag(:)
    integer(kind=PDM_g_num_s), pointer :: face_ln_to_gn(:)
    integer(kind=PDM_l_num_s), pointer :: faceGroupInitFaceGroup(:)
    integer(kind=PDM_l_num_s), pointer :: faceInitFace(:)
    double precision,          pointer :: vtxCoord(:,:)
    integer(kind=PDM_l_num_s), pointer :: vtx_tag(:)
    integer(kind=PDM_g_num_s), pointer :: vtx_ln_to_gn(:)
    integer(kind=PDM_l_num_s), pointer :: vtxInitVtx(:)
    integer(kind=PDM_l_num_s), pointer :: face_group_idx(:)
    integer(kind=PDM_l_num_s), pointer :: face_group(:)
    integer(kind=PDM_g_num_s), pointer :: face_group_ln_to_gn(:)
    integer(kind=PDM_l_num_s), pointer :: face_part_bound_proc_idx(:)
    integer(kind=PDM_l_num_s), pointer :: face_part_bound_part_idx(:)
    integer(kind=PDM_l_num_s), pointer :: face_part_bound(:)

    integer(c_int)                     :: c_i_part
    type(c_ptr)                        :: c_cell_face_idx
    type(c_ptr)                        :: c_cell_face
    type(c_ptr)                        :: c_cell_tag
    type(c_ptr)                        :: c_cell_ln_to_gn
    type(c_ptr)                        :: c_cellInitCellIdx
    type(c_ptr)                        :: c_cellInitCell
    type(c_ptr)                        :: c_face_cell
    type(c_ptr)                        :: c_face_vtx_idx
    type(c_ptr)                        :: c_face_vtx
    type(c_ptr)                        :: c_face_tag
    type(c_ptr)                        :: c_face_ln_to_gn
    type(c_ptr)                        :: c_faceGroupInitFaceGroup
    type(c_ptr)                        :: c_faceInitFace
    type(c_ptr)                        :: c_vtxCoord
    type(c_ptr)                        :: c_vtx_tag
    type(c_ptr)                        :: c_vtx_ln_to_gn
    type(c_ptr)                        :: c_vtxInitVtx
    type(c_ptr)                        :: c_face_group_idx
    type(c_ptr)                        :: c_face_group
    type(c_ptr)                        :: c_face_group_ln_to_gn
    type(c_ptr)                        :: c_face_part_bound_proc_idx
    type(c_ptr)                        :: c_face_part_bound_part_idx
    type(c_ptr)                        :: c_face_part_bound
    integer(c_int)                     :: n_cell
    integer(c_int)                     :: n_face
    integer(c_int)                     :: n_face_part_bound
    integer(c_int)                     :: n_vtx
    integer(c_int)                     :: n_proc
    integer(c_int)                     :: n_total_part
    integer(c_int)                     :: n_face_group
    integer(c_int)                     :: scell_face
    integer(c_int)                     :: sface_vtx
    integer(c_int)                     :: sface_group
    integer(c_int)                     :: sCoarseCellToFineCell

    c_i_part = i_part

    call PDM_part_coarse_mesh_part_dim_get_c (cm,                    &
                                              i_part,                &
                                              n_cell,                &
                                              n_face,                &
                                              n_face_part_bound,     &
                                              n_vtx,                 &
                                              n_proc,                &
                                              n_total_part,          &
                                              n_face_group,          &
                                              scell_face,            &
                                              sface_vtx,             &
                                              sface_group,           &
                                              sCoarseCellToFineCell)

    c_cell_face_idx            = C_NULL_PTR
    c_cell_face                = C_NULL_PTR
    c_cell_tag                 = C_NULL_PTR
    c_cell_ln_to_gn            = C_NULL_PTR
    c_cellInitCellIdx          = C_NULL_PTR
    c_cellInitCell             = C_NULL_PTR
    c_face_cell                = C_NULL_PTR
    c_face_vtx_idx             = C_NULL_PTR
    c_face_vtx                 = C_NULL_PTR
    c_face_tag                 = C_NULL_PTR
    c_face_ln_to_gn            = C_NULL_PTR
    c_faceGroupInitFaceGroup   = C_NULL_PTR
    c_faceInitFace             = C_NULL_PTR
    c_vtxCoord                 = C_NULL_PTR
    c_vtx_tag                  = C_NULL_PTR
    c_vtx_ln_to_gn             = C_NULL_PTR
    c_vtxInitVtx               = C_NULL_PTR
    c_face_group_idx           = C_NULL_PTR
    c_face_group               = C_NULL_PTR
    c_face_group_ln_to_gn      = C_NULL_PTR
    c_face_part_bound_proc_idx = C_NULL_PTR
    c_face_part_bound_part_idx = C_NULL_PTR
    c_face_part_bound          = C_NULL_PTR
    call PDM_part_coarse_mesh_part_get_c (cm,                         &
                                          c_i_part,                   &
                                          c_cell_face_idx,            &
                                          c_cell_face,                &
                                          c_cell_tag,                 &
                                          c_cell_ln_to_gn,            &
                                          c_cellInitCellIdx,          &
                                          c_cellInitCell,             &
                                          c_face_cell,                &
                                          c_face_vtx_idx,             &
                                          c_face_vtx,                 &
                                          c_face_tag,                 &
                                          c_face_ln_to_gn,            &
                                          c_faceGroupInitFaceGroup,   &
                                          c_faceInitFace,             &
                                          c_vtxCoord,                 &
                                          c_vtx_tag,                  &
                                          c_vtx_ln_to_gn,             &
                                          c_vtxInitVtx,               &
                                          c_face_group_idx,           &
                                          c_face_group,               &
                                          c_face_group_ln_to_gn,      &
                                          c_face_part_bound_proc_idx, &
                                          c_face_part_bound_part_idx, &
                                          c_face_part_bound)

    call c_f_pointer(c_cell_face_idx, &
                     cell_face_idx,   &
                     [n_cell+1])

    call c_f_pointer(c_cell_face, &
                     cell_face,   &
                     [scell_face])

    call c_f_pointer(c_cell_tag, &
                     cell_tag,   &
                     [n_cell])

    call c_f_pointer(c_cell_ln_to_gn, &
                     cell_ln_to_gn,   &
                     [n_cell])

    call c_f_pointer(c_cellInitCellIdx, &
                     cellInitCellIdx,   &
                     [n_cell+1])

    call c_f_pointer(c_cellInitCell, &
                     cellInitCell,   &
                     [sCoarseCellToFineCell])

    call c_f_pointer(c_face_cell, &
                     face_cell,   &
                     [2*n_face])

    call c_f_pointer(c_face_vtx_idx, &
                     face_vtx_idx,   &
                     [n_face+1])

    call c_f_pointer(c_face_vtx, &
                     face_vtx,   &
                     [sface_vtx])

    call c_f_pointer(c_face_tag, &
                     face_tag,   &
                     [n_face])

    call c_f_pointer(c_face_ln_to_gn, &
                     face_ln_to_gn,   &
                     [n_face])

    call c_f_pointer(c_faceGroupInitFaceGroup, &
                     faceGroupInitFaceGroup,   &
                     [sface_group])

    call c_f_pointer(c_faceInitFace, &
                     faceInitFace,   &
                     [n_face])

    call c_f_pointer(c_vtxCoord, &
                     vtxCoord,   &
                     [3,n_vtx])

    call c_f_pointer(c_vtx_tag, &
                     vtx_tag,   &
                     [n_vtx])

    call c_f_pointer(c_vtx_ln_to_gn, &
                     vtx_ln_to_gn,   &
                     [n_vtx])

    call c_f_pointer(c_vtxInitVtx, &
                     vtxInitVtx,   &
                     [n_vtx])

    call c_f_pointer(c_face_group_idx, &
                     face_group_idx,   &
                     [n_face_group+1])

    call c_f_pointer(c_face_group, &
                     face_group,   &
                     [sface_group])

    call c_f_pointer(c_face_group_ln_to_gn, &
                     face_group_ln_to_gn,   &
                     [sface_group])

    call c_f_pointer(c_face_part_bound_proc_idx, &
                     face_part_bound_proc_idx,   &
                     [n_proc+1])

    call c_f_pointer(c_face_part_bound_part_idx, &
                     face_part_bound_part_idx,   &
                     [n_total_part+1])

    call c_f_pointer(c_face_part_bound, &
                     face_part_bound,   &
                     [n_face_part_bound])

  end subroutine PDM_part_coarse_mesh_part_get_



  !>
  !!
  !! \brief Return the coloring of a coarse mesh
  !!
  !! \param [in]   cm                  Pointer to \ref PDM_coarse_mesh
  !! \param [in]   i_part              Current partition
  !! \param [out]  cell_color          Cell color (size = n_cell)
  !! \param [out]  face_color          Face color (size = n_face)
  !! \param [out]  thread_color        Thread color (size = n_cell)
  !! \param [out]  hyperplane_color    Hyperplane color (size = n_cell)
  !!
  !!

  subroutine PDM_part_coarse_color_get_ (cm, &
                                         i_part, &
                                         cell_color, &
                                         face_color, &
                                         thread_color, &
                                         hyperplane_color)

    use iso_c_binding
    implicit none

    type(c_ptr), value                 :: cm
    integer, intent(in)                :: i_part
    integer(kind=PDM_l_num_s), pointer :: cell_color(:)
    integer(kind=PDM_l_num_s), pointer :: face_color(:)
    integer(kind=PDM_l_num_s), pointer :: thread_color(:)
    integer(kind=PDM_l_num_s), pointer :: hyperplane_color(:)

    integer(c_int)                     :: c_i_part
    type(c_ptr)                        :: c_cell_color
    type(c_ptr)                        :: c_face_color
    type(c_ptr)                        :: c_thread_color
    type(c_ptr)                        :: c_hyperplane_color
    integer(c_int)                     :: n_cell
    integer(c_int)                     :: n_face
    integer(c_int)                     :: n_face_part_bound
    integer(c_int)                     :: n_vtx
    integer(c_int)                     :: n_proc
    integer(c_int)                     :: n_total_part
    integer(c_int)                     :: n_face_group
    integer(c_int)                     :: scell_face
    integer(c_int)                     :: sface_vtx
    integer(c_int)                     :: sface_group
    integer(c_int)                     :: sCoarseCellToFineCell

    c_i_part = i_part

    call PDM_part_coarse_mesh_part_dim_get_c (cm,                    &
                                              i_part,                &
                                              n_cell,                &
                                              n_face,                &
                                              n_face_part_bound,     &
                                              n_vtx,                 &
                                              n_proc,                &
                                              n_total_part,          &
                                              n_face_group,          &
                                              scell_face,            &
                                              sface_vtx,             &
                                              sface_group,           &
                                              sCoarseCellToFineCell)


    call PDM_part_coarse_color_get_c(cm,                 &
                                     c_i_part,           &
                                     c_cell_color,       &
                                     c_face_color,       &
                                     c_thread_color,     &
                                     c_hyperplane_color)

    call c_f_pointer(c_cell_color, &
                     cell_color,   &
                     [n_cell])

    call c_f_pointer(c_face_color, &
                     face_color,   &
                     [n_face])

    call c_f_pointer(c_thread_color, &
                     thread_color,   &
                     [n_cell])

    call c_f_pointer(c_hyperplane_color, &
                     hyperplane_color,   &
                     [n_cell])

  end subroutine PDM_part_coarse_color_get_



  !>
  !!
  !! \brief Return times
  !!
  !! \param [in]   cm          Pointer to \ref PDM_coarse_mesh
  !! \param [out]  elapsed     Elapsed times (size = 18)
  !! \param [out]  cpu         Cpu times (size = 18)
  !! \param [out]  cpu_user    User cpu times (size = 18)
  !! \param [out]  cpu_sys     System cpu times (size = 18)
  !!
  !!

  subroutine PDM_part_coarse_mesh_time_get_ (cm,       &
                                             elapsed,  &
                                             cpu,      &
                                             cpu_user, &
                                             cpu_sys)

    use iso_c_binding
    implicit none

    type(c_ptr), value            :: cm
    double precision, intent(out) :: elapsed(18)
    double precision, intent(out) :: cpu(18)
    double precision, intent(out) :: cpu_user(18)
    double precision, intent(out) :: cpu_sys(18)

    type(c_ptr)                   :: c_elapsed  
    type(c_ptr)                   :: c_cpu      
    type(c_ptr)                   :: c_cpu_user 
    type(c_ptr)                   :: c_cpu_sys  
    double precision, pointer     :: ptr(:)     


    c_elapsed  = C_NULL_PTR
    c_cpu      = C_NULL_PTR
    c_cpu_user = C_NULL_PTR
    c_cpu_sys  = C_NULL_PTR
    ptr     => null()
    call PDM_part_coarse_mesh_time_get_c (cm,         &
                                          c_elapsed,  &
                                          c_cpu,      &
                                          c_cpu_user, &
                                          c_cpu_sys)

    call c_f_pointer(c_elapsed, &
                     ptr,       &
                     [18])
    elapsed(1:18) = ptr(1:18)

    call c_f_pointer(c_cpu, &
                     ptr,   &
                     [18])
    cpu(1:18) = ptr(1:18)

    call c_f_pointer(c_cpu_user, &
                     ptr,        &
                     [18])
    cpu_user(1:18) = ptr(1:18)

    call c_f_pointer(c_cpu_sys, &
                     ptr,       &
                     [18])
    cpu_sys(1:18) = ptr(1:18)

  end subroutine PDM_part_coarse_mesh_time_get_


  !>
  !!
  !! \brief Get index of a coarse mesh method from it's name
  !!
  !! \param [in]   name   Name of the method
  !! \param [out]  idx    Index of method -1 otherwise
  !!
  !!

  subroutine PDM_coarse_mesh_method_idx_get (name, &
                                             idx)
    use iso_c_binding
    implicit none

    character (len=*), intent(in)  :: name
    integer,           intent(out) :: idx

    interface
      subroutine PDM_coarse_mesh_method_idx_get_c (name, &
                                                   idx)  &
      bind (c, name='PDM_coarse_mesh_method_idx_get')
        use iso_c_binding
        implicit none

        character(c_char) :: name(*)
        integer(c_int)    :: idx

      end subroutine PDM_coarse_mesh_method_idx_get_c
    end interface

    call PDM_coarse_mesh_method_idx_get_c (trim(name)//C_NULL_CHAR, &
                                           idx)

  end subroutine PDM_coarse_mesh_method_idx_get


  !>
  !!
  !! \brief Get the number of coarse mesh method
  !!
  !! \param [out]  n_method  Number of methods
  !!
  !!

  subroutine PDM_coarse_mesh_method_n_get (n_method)

    use iso_c_binding
    implicit none

    integer, intent(out) :: n_method

    interface
      function PDM_coarse_mesh_method_n_get_c () &
      result (n_method)                              &
      bind (c, name='PDM_coarse_mesh_method_n_get')
        use iso_c_binding
        implicit none

        integer(c_int) :: n_method

      end function PDM_coarse_mesh_method_n_get_c
    end interface

    n_method = PDM_coarse_mesh_method_n_get_c ()

  end subroutine PDM_coarse_mesh_method_n_get


  !>
  !!
  !! \brief Get name of a coarse mesh method from it's index
  !!
  !! \param [in]  idx     Index of the method
  !! \param [out] name    Name of the method, '' otherwize
  !!
  !!

  subroutine PDM_coarse_mesh_method_name_get (idx, &
                                              name)
    use iso_c_binding
    implicit none

    integer, intent(in) :: idx
    character(len=*)    :: name

    integer(c_int)      :: l_name

    interface
      subroutine PDM_coarse_mesh_method_name_get_c (idx,    &
                                                    name,   &
                                                    l_name) &
      bind (c, name='PDM_coarse_mesh_method_name_get_cf')
        use iso_c_binding
        implicit none

        integer(c_int), value :: idx
        character(c_char)     :: name(*)
        integer(c_int)        :: l_name

      end subroutine PDM_coarse_mesh_method_name_get_c
    end interface

    call PDM_coarse_mesh_method_name_get_c (idx,    &
                                            name,   &
                                            l_name)

  end subroutine PDM_coarse_mesh_method_name_get


  !>
  !!
  !! \brief Get index of a renumbering face method
  !!
  !! \param [in]   name   Name of the method
  !! \param [out]  idx    Index of method -1 otherwise
  !!
  !!

  subroutine PDM_part_renum_method_face_idx_get (name, &
                                                 idx)
    use iso_c_binding
    implicit none

    character (len=*), intent(in)  :: name
    integer,           intent(out) :: idx

    interface
      subroutine PDM_part_renum_method_face_idx_get_c (name, &
                                                       idx)  &
      bind (c, name='PDM_part_renum_method_face_idx_get')
        use iso_c_binding
        implicit none

        character(c_char) :: name(*)
        integer(c_int)    :: idx

      end subroutine PDM_part_renum_method_face_idx_get_c
    end interface

    call PDM_part_renum_method_face_idx_get_c (trim(name)//C_NULL_CHAR, &
                                               idx)

  end subroutine PDM_part_renum_method_face_idx_get


  !>
  !!
  !! \brief Get index of a renumbering cell method
  !!
  !! \param [in]   name   Name of the method
  !! \param [out]  idx    Index of method -1 otherwise
  !!
  !!

  subroutine PDM_part_renum_method_cell_idx_get (name, &
                                                  idx)
     use iso_c_binding
     implicit none

     character (len=*), intent(in)  :: name
     integer,           intent(out) :: idx

     interface
       subroutine PDM_part_renum_method_cell_idx_get_c (name, &
                                                        idx)  &
       bind (c, name='PDM_part_renum_method_cell_idx_get')
         use iso_c_binding
         implicit none

         character(c_char) :: name(*)
         integer(c_int)    :: idx

       end subroutine PDM_part_renum_method_cell_idx_get_c
     end interface

     call PDM_part_renum_method_cell_idx_get_c (trim(name)//C_NULL_CHAR, &
                                                idx)

  end subroutine PDM_part_renum_method_cell_idx_get


  !>
  !!
  !! \brief Get index of a renumbering edge method
  !!
  !! \param [in]   name   Name of the method
  !! \param [out]  idx    Index of method -1 otherwise
  !!
  !!

  subroutine PDM_part_renum_method_edge_idx_get (name, &
                                                  idx)
     use iso_c_binding
     implicit none

     character (len=*), intent(in)  :: name
     integer,           intent(out) :: idx

     interface
       subroutine PDM_part_renum_method_edge_idx_get_c (name, &
                                                        idx)  &
       bind (c, name='PDM_part_renum_method_edge_idx_get')
         use iso_c_binding
         implicit none

         character(c_char) :: name(*)
         integer(c_int)    :: idx

       end subroutine PDM_part_renum_method_edge_idx_get_c
     end interface

     call PDM_part_renum_method_edge_idx_get_c (trim(name)//C_NULL_CHAR, &
                                                idx)

  end subroutine PDM_part_renum_method_edge_idx_get


  !>
  !!
  !! \brief Get index of a renumbering vtx method
  !!
  !! \param [in]   name   Name of the method
  !! \param [out]  idx    Index of method -1 otherwise
  !!
  !!

  subroutine PDM_part_renum_method_vtx_idx_get (name, &
                                                  idx)
     use iso_c_binding
     implicit none

     character (len=*), intent(in)  :: name
     integer,           intent(out) :: idx

     interface
       subroutine PDM_part_renum_method_vtx_idx_get_c (name, &
                                                        idx)  &
       bind (c, name='PDM_part_renum_method_vtx_idx_get')
         use iso_c_binding
         implicit none

         character(c_char) :: name(*)
         integer(c_int)    :: idx

       end subroutine PDM_part_renum_method_vtx_idx_get_c
     end interface

     call PDM_part_renum_method_vtx_idx_get_c (trim(name)//C_NULL_CHAR, &
                                                idx)

  end subroutine PDM_part_renum_method_vtx_idx_get


  !>
  !!
  !! \brief Get the number of renumbering cell methods
  !!
  !! \param [out] n_method   Number of methods
  !!
  !!

  subroutine PDM_part_n_renum_method_cell_get (n_method)

    use iso_c_binding
    implicit none

    integer, intent(out) :: n_method

    interface
      function PDM_part_n_renum_method_cell_get_c () &
      result (n_method)                              &
      bind (c, name='PDM_part_n_renum_method_cell_get')
        use iso_c_binding
        implicit none

        integer(c_int) :: n_method

      end function PDM_part_n_renum_method_cell_get_c
    end interface

    n_method = PDM_part_n_renum_method_cell_get_c ()

  end subroutine PDM_part_n_renum_method_cell_get


  !>
  !!
  !! \brief Get the number of renumbering face methods
  !!
  !! \param [out]    Number of methods
  !!
  !!

  subroutine PDM_part_n_renum_method_face_get (n_method)

    use iso_c_binding
    implicit none

    integer, intent(out) :: n_method

    interface
      function PDM_part_n_renum_method_face_get_c () &
      result (n_method)                              &
      bind (c, name='PDM_part_n_renum_method_face_get')
        use iso_c_binding
        implicit none

        integer(c_int) :: n_method

      end function PDM_part_n_renum_method_face_get_c
    end interface

    n_method = PDM_part_n_renum_method_face_get_c ()

  end subroutine PDM_part_n_renum_method_face_get


  !>
  !!
  !! \brief Get name of the face renumbering method
  !!
  !! \param [in]  idx     Index of the method
  !! \param [out] name    Name of the method, '' otherwize
  !!
  !!

  subroutine PDM_part_renum_method_face_name_get (idx, &
                                                  name)
    use iso_c_binding
    implicit none

    integer, intent(in) :: idx
    character(len=*)    :: name

    integer(c_int)      :: l_name

    interface
      subroutine PDM_part_renum_method_face_name_get_c (idx,    &
                                                        name,   &
                                                        l_name) &
      bind (c, name='PDM_part_renum_method_face_name_get_cf')
        use iso_c_binding
        implicit none

        integer(c_int), value :: idx
        character(c_char)     :: name(*)
        integer(c_int)        :: l_name

      end subroutine PDM_part_renum_method_face_name_get_c
    end interface

    call PDM_part_renum_method_face_name_get_c (idx,    &
                                                name,   &
                                                l_name)

  end subroutine PDM_part_renum_method_face_name_get


  !>
  !!
  !! \brief Get name of the cell renumbering method
  !!
  !! \param [in]  idx     Index of the method
  !! \param [out] name    Name of the method, '' otherwize
  !!
  !!

  subroutine PDM_part_renum_method_cell_name_get (idx, &
                                                  name)
    use iso_c_binding
    implicit none

    integer, intent(in) :: idx
    character(len=*)    :: name
    integer(c_int)      :: l_name

    interface
      subroutine PDM_part_renum_method_cell_name_get_c (idx,    &
                                                        name,   &
                                                        l_name) &
      bind (c, name='PDM_part_renum_method_cell_name_get_cf')
        use iso_c_binding
        implicit none

        integer(c_int), value :: idx
        character(c_char)     :: name(*)
        integer(c_int)        :: l_name

      end subroutine PDM_part_renum_method_cell_name_get_c
    end interface

    call PDM_part_renum_method_cell_name_get_c (idx,    &
                                                name,   &
                                                l_name)

  end subroutine PDM_part_renum_method_cell_name_get


  !>
  !!
  !! \brief Add option for anisotropic mesh agglomeration
  !!
  !! \param [out]  cm                         Pointer to \ref PDM_coarse_mesh_t object
  !! \param [in]   anisotropic_option
  !!

#ifdef PDM_HAVE_ANISO_AGGLO

  subroutine PDM_part_coarse_mesh_add_option_anisotropic (cm,                 &
                                                          anisotropic_option)
    use iso_c_binding
    implicit none

    type(c_ptr), value            :: cm
    integer(pdm_l_num_s), pointer :: anisotropic_option(:)

    interface
      subroutine PDM_part_coarse_mesh_add_option_anisotropic_f (cm,                 &
                                                                anisotropic_option) &
      bind (c, name='PDM_part_coarse_mesh_add_option_anisotropic')
        use iso_c_binding
        implicit none

        type(c_ptr), value  :: cm
        type(c_ptr), value  :: anisotropic_option

      end subroutine PDM_part_coarse_mesh_add_option_anisotropic_f
    end interface

    call PDM_part_coarse_mesh_add_option_anisotropic_f (cm, c_loc(anisotropic_option))

  end subroutine PDM_part_coarse_mesh_add_option_anisotropic

  !>
  !!
  !! \brief Add isotropic array to current coarse mesh
  !!
  !! \param [in]   cm                           Pointer to \ref PDM_coarse_mesh_t object
  !! \param [in]   i_part                       Current partition
  !! \param [out]  agglomeration_lines_init
  !! \param [out]  agglomeration_lines_init_idx
  !! \param [out]  is_on_fine_bnd_init
  !!

  subroutine PDM_part_coarse_mesh_part_set_anisotropic_info (cm,                                &
                                                             i_part,                            &
                                                             agglomeration_lines_init,          &
                                                             agglomeration_lines_init_idx,      &
                                                             agglomeration_lines_init_idx_size, &
                                                             is_on_fine_bnd_init)
    use iso_c_binding
    implicit none

    type(c_ptr), value            :: cm
    integer, intent(in)           :: i_part
    integer(pdm_l_num_s), pointer :: agglomeration_lines_init(:)
    integer(pdm_l_num_s), pointer :: agglomeration_lines_init_idx(:)
    integer, intent(in)           :: agglomeration_lines_init_idx_size
    integer(pdm_l_num_s), pointer :: is_on_fine_bnd_init(:)

    interface
      subroutine PDM_part_coarse_mesh_part_set_anisotropic_info_f (cm,                                &
                                                                   i_part,                            &
                                                                   agglomeration_lines_init,          &
                                                                   agglomeration_lines_init_idx,      &
                                                                   agglomeration_lines_init_idx_size, &
                                                                   is_on_fine_bnd_init)               &
      bind (c, name='PDM_part_coarse_mesh_part_set_anisotropic_info')

        use iso_c_binding
        implicit none

        type(c_ptr),    value  :: cm
        integer(c_int), value  :: i_part
        type(c_ptr),    value  :: agglomeration_lines_init
        type(c_ptr),    value  :: agglomeration_lines_init_idx
        integer(c_int), value  :: agglomeration_lines_init_idx_size
        type(c_ptr),    value  :: is_on_fine_bnd_init

      end subroutine PDM_part_coarse_mesh_part_set_anisotropic_info_f
    end interface

    call PDM_part_coarse_mesh_part_set_anisotropic_info_f (cm,                                  &
                                                           i_part,                              &
                                                           c_loc(agglomeration_lines_init),     &
                                                           c_loc(agglomeration_lines_init_idx), &
                                                           agglomeration_lines_init_idx_size,   &
                                                           c_loc(is_on_fine_bnd_init))

  end subroutine PDM_part_coarse_mesh_part_set_anisotropic_info


  !>
  !!
  !! \brief Return a mesh partition
  !!
  !! \param [in]   cm                         Pointer to \ref PDM_coarse_mesh_t object
  !! \param [in]   i_part                     Current partition
  !! \param [out]  agglomeration_lines
  !! \param [out]  agglomeration_lines_idx
  !! \param [out]  is_on_fine_bnd
  !!

  subroutine PDM_part_coarse_mesh_part_get_anisotropic_info (cm,                           &
                                                             i_part,                       &
                                                             agglomeration_lines,          &
                                                             agglomeration_lines_idx,      &
                                                             agglomeration_lines_idx_size, &
                                                             is_on_fine_bnd)
    use iso_c_binding
    implicit none

    type(c_ptr), value            :: cm
    integer, intent(in)           :: i_part
    integer(pdm_l_num_s), pointer :: agglomeration_lines(:)
    integer(pdm_l_num_s), pointer :: agglomeration_lines_idx(:)
    integer, intent(out)          :: agglomeration_lines_idx_size
    integer(pdm_l_num_s), pointer :: is_on_fine_bnd(:)

    type(c_ptr)                   :: c_agglomeration_lines          
    type(c_ptr)                   :: c_agglomeration_lines_idx      
    integer(c_int)                :: c_agglomeration_lines_idx_size
    type(c_ptr)                   :: c_is_on_fine_bnd               
    integer                       :: n_cell
    integer                       :: n_face
    integer                       :: n_face_part_bound
    integer                       :: n_vtx
    integer                       :: n_proc
    integer                       :: n_total_part
    integer                       :: n_face_group
    integer                       :: scell_face
    integer                       :: sface_vtx
    integer                       :: sface_group
    integer                       :: sCoarseCellToFineCell

    interface
      subroutine PDM_part_coarse_mesh_part_get_anisotropic_info_f (cm,                           &
                                                                   i_part,                       &
                                                                   agglomeration_lines,          &
                                                                   agglomeration_lines_idx,      &
                                                                   agglomeration_lines_idx_size, &
                                                                   is_on_fine_bnd)               &
      bind (c, name='PDM_part_coarse_mesh_part_get_anisotropic_info')

        use iso_c_binding
        implicit none

        type(c_ptr),    value  :: cm
        integer(c_int), value  :: i_part
        type(c_ptr)            :: agglomeration_lines
        type(c_ptr)            :: agglomeration_lines_idx
        integer(c_int)         :: agglomeration_lines_idx_size
        type(c_ptr)            :: is_on_fine_bnd

      end subroutine PDM_part_coarse_mesh_part_get_anisotropic_info_f
    end interface

    c_agglomeration_lines          = C_NULL_PTR
    c_agglomeration_lines_idx      = C_NULL_PTR
    c_is_on_fine_bnd               = C_NULL_PTR
    call PDM_part_coarse_mesh_part_get_anisotropic_info_f (cm,                             &
                                                           i_part,                         &
                                                           c_agglomeration_lines,          &
                                                           c_agglomeration_lines_idx,      &
                                                           c_agglomeration_lines_idx_size, &
                                                           c_is_on_fine_bnd)

    agglomeration_lines_idx_size = c_agglomeration_lines_idx_size

    call c_f_pointer(c_agglomeration_lines_idx,        &
                     agglomeration_lines_idx,          &
                     [c_agglomeration_lines_idx_size])

    call c_f_pointer(c_agglomeration_lines,                                     &
                     agglomeration_lines,                                       &
                     [agglomeration_lines_idx(c_agglomeration_lines_idx_size)])


    call PDM_part_coarse_mesh_part_dim_get_ (cm,                    &
                                             i_part,                &
                                             n_cell,                &
                                             n_face,                &
                                             n_face_part_bound,     &
                                             n_vtx,                 &
                                             n_proc,                &
                                             n_total_part,          &
                                             n_face_group,          &
                                             scell_face,            &
                                             sface_vtx,             &
                                             sface_group,           &
                                             sCoarseCellToFineCell)

    call c_f_pointer(c_is_on_fine_bnd, &
                     is_on_fine_bnd,   &
                     [n_cell])

  end subroutine PDM_part_coarse_mesh_part_get_anisotropic_info

#endif

end module pdm_part
