#include "pdm_configf.h"
#ifdef PDM_LONG_G_NUM
  integer, parameter :: pdm_g_num_s = 8
#else
  integer, parameter :: pdm_g_num_s = 4
#endif
  integer, parameter :: pdm_l_num_s = 4

  integer, parameter :: pdm_max_char_length = 100

  integer, parameter :: PDM_STRIDE_CST_INTERLACED  = 0
  integer, parameter :: PDM_STRIDE_CST_INTERLEAVED = 1
  integer, parameter :: PDM_STRIDE_VAR_INTERLACED  = 2

  integer, parameter :: PDM_OWNERSHIP_KEEP                 = 0
  integer, parameter :: PDM_OWNERSHIP_USER                 = 1
  integer, parameter :: PDM_OWNERSHIP_UNGET_RESULT_IS_FREE = 2

  integer, parameter :: PDM_MESH_ENTITY_CELL = 0  !< Cell entity
  integer, parameter :: PDM_MESH_ENTITY_FACE = 1  !< Face entity
  integer, parameter :: PDM_MESH_ENTITY_EDGE = 2  !< Edge entity
  integer, parameter :: PDM_MESH_ENTITY_VTX  = 3  !< Vertex entity

  integer, parameter :: PDM_MESH_NATURE_NODAL_SHARED   = 0  !< PDM_mesh_nodal
  integer, parameter :: PDM_MESH_NATURE_MESH_SETTED    = 1  !< PDm_surface_mesh

  integer, parameter :: PDM_MPI_COMM_KIND_P2P                                = 0
  integer, parameter :: PDM_MPI_COMM_KIND_COLLECTIVE                         = 1
  integer, parameter :: PDM_MPI_COMM_KIND_NEIGHBOR_COLLECTIVE                = 2
  integer, parameter :: PDM_MPI_COMM_KIND_WIN_SHARED_AND_P2P                 = 3
  integer, parameter :: PDM_MPI_COMM_KIND_WIN_SHARED_AND_COLLECTIVE          = 4
  integer, parameter :: PDM_MPI_COMM_KIND_WIN_SHARED_AND_NEIGHBOR_COLLECTIVE = 5
  integer, parameter :: PDM_MPI_COMM_KIND_WIN_RMA                            = 6

  integer, parameter :: PDM_CONNECTIVITY_TYPE_CELL_ELMT = 0
  integer, parameter :: PDM_CONNECTIVITY_TYPE_CELL_CELL = 1
  integer, parameter :: PDM_CONNECTIVITY_TYPE_CELL_FACE = 2
  integer, parameter :: PDM_CONNECTIVITY_TYPE_CELL_EDGE = 3
  integer, parameter :: PDM_CONNECTIVITY_TYPE_CELL_VTX  = 4
  integer, parameter :: PDM_CONNECTIVITY_TYPE_FACE_ELMT = 5
  integer, parameter :: PDM_CONNECTIVITY_TYPE_FACE_CELL = 6
  integer, parameter :: PDM_CONNECTIVITY_TYPE_FACE_FACE = 7
  integer, parameter :: PDM_CONNECTIVITY_TYPE_FACE_EDGE = 8
  integer, parameter :: PDM_CONNECTIVITY_TYPE_FACE_VTX  = 9
  integer, parameter :: PDM_CONNECTIVITY_TYPE_EDGE_ELMT = 10
  integer, parameter :: PDM_CONNECTIVITY_TYPE_EDGE_CELL = 11
  integer, parameter :: PDM_CONNECTIVITY_TYPE_EDGE_FACE = 12
  integer, parameter :: PDM_CONNECTIVITY_TYPE_EDGE_EDGE = 13
  integer, parameter :: PDM_CONNECTIVITY_TYPE_EDGE_VTX  = 14
  integer, parameter :: PDM_CONNECTIVITY_TYPE_VTX_ELMT  = 15
  integer, parameter :: PDM_CONNECTIVITY_TYPE_VTX_CELL  = 16
  integer, parameter :: PDM_CONNECTIVITY_TYPE_VTX_FACE  = 17
  integer, parameter :: PDM_CONNECTIVITY_TYPE_VTX_EDGE  = 18
  integer, parameter :: PDM_CONNECTIVITY_TYPE_VTX_VTX   = 19
  integer, parameter :: PDM_CONNECTIVITY_TYPE_ELMT_CELL = 20
  integer, parameter :: PDM_CONNECTIVITY_TYPE_ELMT_FACE = 21
  integer, parameter :: PDM_CONNECTIVITY_TYPE_ELMT_EDGE = 22
  integer, parameter :: PDM_CONNECTIVITY_TYPE_ELMT_VTX  = 23

  integer, parameter :: PDM_EXTRACT_PART_KIND_LOCAL         = 0
  integer, parameter :: PDM_EXTRACT_PART_KIND_REEQUILIBRATE = 1
  integer, parameter :: PDM_EXTRACT_PART_KIND_FROM_TARGET   = 2

  integer, parameter :: PDM_BOUND_TYPE_ELMT = 0
  integer, parameter :: PDM_BOUND_TYPE_CELL = 1
  integer, parameter :: PDM_BOUND_TYPE_FACE = 2
  integer, parameter :: PDM_BOUND_TYPE_EDGE = 3
  integer, parameter :: PDM_BOUND_TYPE_VTX  = 4
  integer, parameter :: PDM_BOUND_TYPE_MAX  = 5

  integer, parameter :: PDM_MESH_NODAL_POINT           = 0
  integer, parameter :: PDM_MESH_NODAL_BAR2            = 1
  integer, parameter :: PDM_MESH_NODAL_TRIA3           = 2
  integer, parameter :: PDM_MESH_NODAL_QUAD4           = 3
  integer, parameter :: PDM_MESH_NODAL_POLY_2D         = 4
  integer, parameter :: PDM_MESH_NODAL_TETRA4          = 5
  integer, parameter :: PDM_MESH_NODAL_PYRAMID5        = 6
  integer, parameter :: PDM_MESH_NODAL_PRISM6          = 7
  integer, parameter :: PDM_MESH_NODAL_HEXA8           = 8
  integer, parameter :: PDM_MESH_NODAL_POLY_3D         = 9
  integer, parameter :: PDM_MESH_NODAL_BARHO           = 10
  integer, parameter :: PDM_MESH_NODAL_TRIAHO          = 11
  integer, parameter :: PDM_MESH_NODAL_QUADHO          = 12
  integer, parameter :: PDM_MESH_NODAL_TETRAHO         = 13
  integer, parameter :: PDM_MESH_NODAL_PYRAMIDHO       = 14
  integer, parameter :: PDM_MESH_NODAL_PRISMHO         = 15
  integer, parameter :: PDM_MESH_NODAL_HEXAHO          = 16
  integer, parameter :: PDM_MESH_NODAL_BARHO_BEZIER    = 17
  integer, parameter :: PDM_MESH_NODAL_TRIAHO_BEZIER   = 18
  integer, parameter :: PDM_MESH_NODAL_N_ELEMENT_TYPES = 19

  integer, parameter :: PDM_PART_SPLIT_PARMETIS        = 1
  integer, parameter :: PDM_PART_SPLIT_PTSCOTCH        = 2
  integer, parameter :: PDM_PART_SPLIT_HILBERT         = 3

  integer, parameter :: PDM_SPLIT_DUAL_WITH_PARMETIS   = 1
  integer, parameter :: PDM_SPLIT_DUAL_WITH_PTSCOTCH   = 2
  integer, parameter :: PDM_SPLIT_DUAL_WITH_HILBERT    = 3
  integer, parameter :: PDM_SPLIT_DUAL_WITH_IMPLICIT   = 4
