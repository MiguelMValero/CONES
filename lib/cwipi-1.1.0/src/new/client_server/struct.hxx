#ifndef __STRUCT_H__
#define __STRUCT_H__
/*
  This file is part of the CWIPI library.

  Copyright (C) 2022-2023  ONERA

  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License as published by the Free Software Foundation; either
  version 3 of the License, or (at your option) any later version.

  This library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public
  License along with this library. If not, see <http://www.gnu.org/licenses/>.
*/

/*
  This file is inspired from OpenPALM.
  OpenPALM is a free software under the GNU Lesser General Public License.
  See: https://www.cerfacs.fr/globc/PALM_WEB/
*/

/*----------------------------------------------------------------------------
 * Standard C++ library headers
 *----------------------------------------------------------------------------*/

#include <string>
#include <map>
#include <tuple>

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include "mpi.h"
#include "cwp.h"
#include "pdm_logging.h"

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force back Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */

/*=============================================================================
 * Macro definitions
 *============================================================================*/

#define CWP_FLAG_VERBOSE    1
#define CWP_FLAG_NOVERBOSE  0

/*============================================================================
 * Types definition
 *============================================================================*/

struct t_field
{
  t_field() {
    data       = NULL;
    srcs       = NULL;
    c_tgts     = NULL;
    u_tgts     = NULL;
    n_entities = NULL;
    size       = NULL;
  }

  // Involved_srcs_get
  int *srcs;

  // Computed_tgts_get
  int *c_tgts;

  // Uncomputed_tgts_get
  int *u_tgts;

  // Field_data_set
  double **data;
  CWP_Dof_location_t location;

  // Field_wait_irecv
  CWP_Field_map_t map_type;
  int             n_component;
  int            *n_entities;
  int            *size;
};

struct t_global_data
{
  t_global_data() {
    send_data = NULL;
    recv_data = NULL;
  }

  // Receive data
  size_t s_recv_entity;
  int recv_stride;
  int n_recv_entity;
  void *recv_data;

  // Send data
  void *send_data;
};

struct t_part_data
{
  t_part_data() {
    exch_type = (CWP_PartData_exch_t) 123;
    gnum_elt  = NULL;
    n_elt     = NULL;
    n_part    = 0;
  }

  CWP_PartData_exch_t     exch_type;

  CWP_g_num_t            **gnum_elt;
  int                     *n_elt;
  int                      n_part;
  std::map<int, void **>   data;
  std::map<int, int>       s_unit;
};

struct t_property
{
  t_property() {
    property_name  = NULL;
    property_type  = (CWP_Type_t) -1;
    property_value = NULL;
  }

  char       *property_name;
  CWP_Type_t  property_type;
  char       *property_value;

};

struct t_block
{
  t_block() {
    connec_faces_idx  = NULL;
    connec_faces      = NULL;
    connec_cells_idx  = NULL;
    connec_cells      = NULL;
    cell_global_num   = NULL;
    connec_idx        = NULL;
    connec            = NULL;
    elt_global_num    = NULL;
    std_connec        = NULL;
    std_global_num    = NULL;
    ho_std_connec     = NULL;
    ho_std_global_num = NULL;
  }

  // Mesh_interf_c_poly_block_set
  int         **connec_faces_idx;
  int         **connec_faces;
  int         **connec_cells_idx;
  int         **connec_cells;
  CWP_g_num_t **cell_global_num;

  // Mesh_interf_f_poly_block_set
  int         **connec_idx;
  int         **connec;
  CWP_g_num_t **elt_global_num;

  // Mesh_interf_block_std_set
  int         **std_connec;
  CWP_g_num_t **std_global_num;

  // Mesh_interf_block_ho_set
  int         **ho_std_connec;
  CWP_g_num_t **ho_std_global_num;
};


struct t_coupling
{
  t_coupling() {
    vtx_coord       = NULL;
    vtx_global_num  = NULL;
    usr_coord       = NULL;
    usr_global_num  = NULL;
    face_edge_idx   = NULL;
    face_edge       = NULL;
    edge_vtx        = NULL;
    face_global_num = NULL;
    cell_face_idx   = NULL;
    cell_face       = NULL;
    face_vtx_idx    = NULL;
    face_vtx        = NULL;
    cell_global_num = NULL;
    ijk_grid        = NULL;
    n_vtx           = NULL;
    n_user_vtx      = NULL;
    n_elt           = NULL;
  }

  // mesh
  int *n_elt;

  // Mesh_interf_vtx_set
  int          *n_vtx;
  double      **vtx_coord;
  CWP_g_num_t **vtx_global_num;

  // User_tgt_pts_set
  int          *n_user_vtx;
  double      **usr_coord;
  CWP_g_num_t **usr_global_num;

  // Mesh_interf_from_faceedge_set
  int         **face_edge_idx;
  int         **face_edge;
  int         **edge_vtx;
  CWP_g_num_t **face_global_num;

  // Mesh_interf_from_cellface_set
  int         **cell_face_idx;
  int         **cell_face;
  int         **face_vtx_idx;
  int         **face_vtx;
  CWP_g_num_t **cell_global_num;

  // Mesh_interf_ho_ordering_from_IJK_set
  int *ijk_grid;

  // Cpl_create
  int                n_part;
  CWP_Dynamic_mesh_t mesh_dynamic;

  // Spatial_interp_property_set
  std::map<std::string, t_property> property;

  // field
  std::map<std::string, t_field> field;

  // global_data
  std::map<std::string, t_global_data> global_data;

  // part_data
  std::map<std::string, t_part_data> part_data;

  // block_id
  std::map<int, t_block> block;
};

struct t_cwp
{

  // Codes_list_get
  int          n_code_names;
  char       **code_names;

  // Loc_codes_list_get
  int          n_loc_code_names;
  char       **loc_code_names;

  // Param_list_get, Param_get
  int                             n_param_names;
  char                          **param_names;
  std::map<std::string, char *>   char_param_value;

  // Output_file_set
  FILE *output_file;

  // cpl_id
  std::map<std::string, t_coupling> coupling;
};

struct t_server_mpi
{
  MPI_Comm  global_comm;
  MPI_Comm *intra_comms;
};

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* __STRUCT_H__ */
