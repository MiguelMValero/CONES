#ifndef __CLIENT_H__
#define __CLIENT_H__
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
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "cwp.h"


#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force back Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */

/*============================================================================
 * Types definition
 *============================================================================*/

/**
 * \struct t_client
 * \brief Client storage structure
 */

typedef struct t_client
{
  MPI_Comm comm;             /**< comm */
  int      i_rank;           /**< i_rank */
  char     *code_name;       /**< code_name */
  int      server_port;      /**< server_port */
  int      flags;            /**< flags */
  int      socket;           /**< socket */
  int      max_msg_size;     /**< max_msg_size */
  int      listen_socket;    /**< listen_socket */
  int      connected_socket; /**< connected_socket */
  int      client_endianess; /**< client_endianess */
  int      server_endianess; /**< server_endianess */
  char     server_name[256]; /**< server_name */

}t_client,*p_client;

/*=============================================================================
 * Client CWIPI function interfaces
 *============================================================================*/

/**
 * \brief Initialize CWIPI.
 *
 * \param [in]  comm           MPI intra communicators of each code
 * \param [in]  config         Configuration file name
 * \param [in]  code_name      Name of the code on current rank
 * \param [in]  is_active_rank Does current rank have to be used by CWIPI
 *
 */

void
CWP_client_Init
(
        MPI_Comm           comm,
        char              *config,
  const char              *code_name,
  const CWP_Status_t       is_active_rank
);

/**
 *
 * \brief Finalize CWIPI.
 *
 */

void
CWP_client_Finalize
(
 void
);

/**
 *
 * \brief Param_lock CWIPI.
 *
 * \param [in]  code_name  Code to lock
 *
 */

void
CWP_client_Param_lock
(
const char *code_name
);

/**
 *
 * \brief Param_unlock CWIPI.
 *
 * \param [in]  code_name  Code to unlock
 *
 */


void
CWP_client_Param_unlock
(
const char *code_name
);

/**
 *
 * \brief Param_add CWIPI.
 *
 * \param [in] local_code_name  Local code name
 * \param [in] param_name       Parameter name
 * \param [in] data_type        Parameter type
 * \param [in] initial_value    Initial value
 *
 */

void
CWP_client_Param_add
(
 const char        *local_code_name,
 const char        *param_name,
 const CWP_Type_t  data_type,
 void              *initial_value
);

/**
 *
 * \brief Param_get CWIPI.
 *
 * \param [in]  code_name  Local or distant code name
 * \param [in]  param_name Parameter name
 * \param [in]  data_type  Parameter type
 * \param [out] value      Parameter value
 *
 */

void
CWP_client_Param_get
(
 const char       *code_name,
 const char       *param_name,
 const CWP_Type_t  data_type,
 void             *value
);

/**
 *
 * \brief Param_set CWIPI.
 *
 * \param [in] local_code_name  Local code name
 * \param [in] param_name       Parameter name
 * \param [in] data_type        Parameter type
 * \param [in] value            Value
 *
 */

void
CWP_client_Param_set
(
 const char             *local_code_name,
 const char             *param_name,
 const CWP_Type_t        data_type,
 void                   *value
);

/**
 *
 * \brief Param_del CWIPI.
 *
 * \param [in] local_code_name  Local code name
 * \param [in] param_name       Parameter name
 * \param [in] data_type        Parameter type
 *
 */

void
CWP_client_Param_del
(
 const char       *local_code_name,
 const char       *param_name,
 const CWP_Type_t  data_type
);

/**
 *
 * \brief Param_n_get CWIPI.
 *
 * \param [in] code_name       Local or distant code name
 * \param [in] data_type       Parameter type,
 *
 * return  Number of parameters
 *
 */

int
CWP_client_Param_n_get
(
 const char             *code_name,
 const CWP_Type_t        data_type
);

/**
 *
 * \brief Param_list_get CWIPI.
 *
 * It is up to the user to free the retreived pointer array
 *
 * \param [in]  code_name      Local or distant code name
 * \param [in]  data_type      Parameter type,
 * \param [out] nParam         Number of parameters
 * \param [out] paramNames     Parameter names
 *
 */

void
CWP_client_Param_list_get
(
 const char             *code_name,
 const CWP_Type_t        data_type,
 int                    *nParam,
 char                 ***paramNames
);

/**
 *
 * \brief Param_is CWIPI.
 *
 * \param [in] code_name      Local or distant code name
 * \param [in] param_name     Parameter name
 * \param [in] data_type      Parameter type
 *
 * return  1 : true / 0 : false
 *
 */

int
CWP_client_Param_is
(
 const char             *code_name,
 const char             *param_name,
 const CWP_Type_t        data_type
);

/**
 *
 * \brief Param_reduce CWIPI.
 *
 * The parameter name has to be the same for all codes.
 *
 * \param [in]  op           Operation
 * \param [in]  param_name   Parameter name
 * \param [in]  data_type    Parameter type,
 * \param [out] res          Result
 * \param [in]  nCode        Number of codes
 * \param       code_names   Codes name
 *
 */

void
CWP_client_Param_reduce
(
 const CWP_Op_t    op,
 const char       *param_name,
 const CWP_Type_t  data_type,
 void             *res,
 const int         nCode,
 const char      **code_names
);

/**
 * \brief Cpl_create CWIPI.
 *
 * \param [in]  local_code_name     Local code name
 * \param [in]  cpl_id              Coupling identifier
 * \param [in]  coupled_code_name   Distant or local coupled code name
 * \param [in]  entities_dim        Coupling interface
 * \param [in]  comm_type           Communication type
 * \param [in]  spatial_interp      Spatial interpolation method
 * \param [in]  n_part              Number of interface partition
 * \param [in]  displacement        Mesh moving status
 * \param [in]  recv_freq_type      Type of receiving frequency
 *
 */

void
CWP_client_Cpl_create
(
 const char                *local_code_name,
 const char                *cpl_id,
 const char                *coupled_code_name,
 CWP_Interface_t            entities_dim,
 const CWP_Comm_t           comm_type,
 const CWP_Spatial_interp_t spatial_interp,
 const int                  n_part,
 const CWP_Dynamic_mesh_t   displacement,
 const CWP_Time_exch_t      recv_freq_type
);

/**
 *
 * \brief CWP_Cpl_barrier CWIPI.
 *
 * \param [in] local_code_name  Local code name
 * \param [in] cpl_id           Coupling identifier
 *
 */

void
CWP_client_Cpl_barrier
(
 const char *local_code_name,
 const char *cpl_id
);

/**
 *
 * \brief Cpl_del CWIPI.
 *
 * \param [in] local_code_name  Local code name
 * \param [in] cpl_id           Coupling identifier
 *
 */

void
CWP_client_Cpl_del
(
 const char *local_code_name,
 const char *cpl_id
);

/**
 * \brief State_update CWIPI.
 *
 * \param [in] local_code_name  Local code name
 * \param [in] state            State
 *
 */

void
CWP_client_State_update
(
 const char* local_code_name,
 const CWP_State_t state
);

/**
 * \brief Time_step_beg CWIPI.
 *
 * \param [in] local_code_name  Local code name
 * \param [in] current_time Current time
 *
 */

void
CWP_client_Time_step_beg
(
 const char* local_code_name,
 const double current_time
);

/**
 * \brief Time_step_end CWIPI.
 *
 * \param [in] local_code_name  Local code name
 *
 */

void
CWP_client_Time_step_end
(
 const char* local_code_name
);

/**
 * \brief Output_file_set CWIPI.
 *
 * /!\ one file per code
 *
 * \param [in] output_file    Output file directory or name
 *
 */

void
CWP_client_Output_file_set
(
  FILE *output_file
);


/**
 * \brief User_structure_set CWIPI.
 *
 * This structure can be called into a callback
 *
 * \param [in] local_code_name  Local code name
 * \param [in] user_structure   User structure
 *
 */

void
CWP_client_User_structure_set
(
 const char* local_code_name,
       void* user_structure
);


/**
 * \brief User_structure_get CWIPI.
 *
 * This structure can be called into a callback
 *
 * \param [in] local_code_name  Local code name
 *
 * \return  User structure
 *
 */

void *
CWP_client_User_structure_get
(
 const char* local_code_name
);

/**
 * \brief State_get CWIPI.
 *
 * \param [in]  code_name    Code name
 *
 * \return      Code state
 */

CWP_State_t
CWP_client_State_get
(
 const char    *code_name
);


/**
 * \brief Codes_nb_get CWIPI.
 *
 * \return Number of codes
 *
 */

int
CWP_client_Codes_nb_get
(
 void
);

/**
 * \brief Codes_list_get CWIPI.
 *
 * \return list of codes.
 */

const char **
CWP_client_Codes_list_get
(
void
);


/**
 * \brief Loc_codes_nb_get CWIPI.
 *
 * \return number of local codes.
 */

int
CWP_client_Loc_codes_nb_get
(
 void
);


/**
 * \brief Loc_codes_list_get CWIPI.
 *
 * \return list of local codes.
 */

const char **
CWP_client_Loc_codes_list_get
(
 void
);

/**
 * \brief Properties_dump CWIPI.
 *
 */

void
CWP_client_Properties_dump
(
void
);

/**
 *
 * \brief Computed_tgts_bcast_enable CWIPI.
 *
 * \param [in] local_code_name  Local code name
 * \param [in] cpl_id           Coupling identifier
 * \param [in] field_id         Field identifier
 */

void
CWP_client_Computed_tgts_bcast_enable
(
 const char *local_code_name,
 const char *cpl_id,
 const char *field_id
);

/**
 *
 * \brief Involved_srcs_bcast_enable CWIPI.
 *
 * \param [in] local_code_name  Local code name
 * \param [in] cpl_id           Coupling identifier
 * \param [in] field_id         Field identifier
 */

void
CWP_client_Involved_srcs_bcast_enable
(
 const char *local_code_name,
 const char *cpl_id,
 const char *field_id
);

/**
 *
 * \brief N_uncomputed_tgts_get CWIPI.
 *
 * \param [in] local_code_name  Local code name
 * \param [in] cpl_id           Coupling identifier
 * \param [in] field_id         Field identifier
 * \param [in] i_part           Current partition
 *
 * \return                Number of uncomputed targets
 */

int
CWP_client_N_uncomputed_tgts_get
(
 const char *local_code_name,
 const char *cpl_id,
 const char *field_id,
 const int   i_part
);


/**
 *
 * \brief Uncomputed_tgts_get CWIPI.
 *
 * \param [in] local_code_name  Local code name
 * \param [in] cpl_id           Coupling identifier
 * \param [in] field_id         Field identifier
 * \param [in] i_part           Current partition
 *
 * \return                Uncomputed targets
 */

const int *
CWP_client_Uncomputed_tgts_get
(
 const char *local_code_name,
 const char *cpl_id,
 const char *field_id,
 const int   i_part
);

/**
 *
 * \brief N_computed_tgts_get CWIPI.
 *
 * \param [in] local_code_name  Local code name
 * \param [in] cpl_id           Coupling identifier
 * \param [in] field_id         Field identifier
 * \param [in] i_part           Current partition
 *
 * \return                Number of computed targets
 */

int
CWP_client_N_computed_tgts_get
(
 const char *local_code_name,
 const char *cpl_id,
 const char *field_id,
 const int   i_part
);

/**
 *
 * \brief Computed_tgts_get CWIPI.
 *
 * \param [in] local_code_name  Local code name
 * \param [in] cpl_id           Coupling identifier
 * \param [in] field_id         Field identifier
 * \param [in] i_part           Current partition
 *
 * \return                Computed targets
 */

const int *
CWP_client_Computed_tgts_get
(
 const char *local_code_name,
 const char *cpl_id,
 const char *field_id,
 const int   i_part
);


/**
 *
 * \brief N_involved_srcs_get CWIPI.
 *
 * \param [in] local_code_name  Local code name
 * \param [in] cpl_id           Coupling identifier
 * \param [in] field_id         Field identifier
 * \param [in] i_part           Current partition
 *
 * \return               Number of involved sources
 */

int
CWP_client_N_involved_srcs_get
(
 const char *local_code_name,
 const char *cpl_id,
 const char *field_id,
 const int   i_part
);


/**
 *
 * \brief Involved_srcs_get CWIPI.
 *
 * \param [in] local_code_name  Local code name
 * \param [in] cpl_id           Coupling identifier
 * \param [in] field_id         Field identifier
 * \param [in] i_part           Current partition
 *
 * \return              Involved sources
 */

const int *
CWP_client_Involved_srcs_get
(
 const char *local_code_name,
 const char *cpl_id,
 const char *field_id,
 const int   i_part
);

/**
 * \brief Spatial_interp_weights_compute CWIPI.
 *
 * \param [in]  local_code_name     Local code name
 * \param [in]  cpl_id              Coupling identifier
 *
 */

void
CWP_client_Spatial_interp_weights_compute
(
 const char     *local_code_name,
 const char     *cpl_id
);


/**
 * \brief Set a property of the spatial interpolation algorithm.
 *
 * Use "n_neighbors" and "polyfit_degree" for the nearest neighbors
 * algorithms. Use "tolerance" for the location algorithm.
 *
 * \param [in]  local_code_name  Local code name
 * \param [in]  cpl_id           Coupling identifier
 * \param [in]  property_name    Name of the property
 * \param [in]  property_type    Type of the property
 * \param [in]  property_value   Value of the property
 *
 */

void
CWP_client_Spatial_interp_property_set
(
 const char       *local_code_name,
 const char       *cpl_id,
 const char       *property_name,
 const CWP_Type_t  property_type,
 const char       *property_value
);

/**
 * \brief Visu_set CWIPI.
 *
 * \param [in]  local_code_name  Local code name
 * \param [in]  cpl_id           Coupling identifier
 * \param [in]  freq             Output frequency
 * \param [in]  format           Output format to visualize exchanged fieldsDouble
 *                               on the coupled mesh. Choice between :
 *                               - "EnSight Gold"
 * \param [in]  format_option   Output options "opt1, opt2, ..."
 *                               - text : output text files
 *                               - binary : output binary files (default)
 */

void
CWP_client_Visu_set
(
 const char                 *local_code_name,
 const char                 *cpl_id,
 const int                   freq,
 const CWP_Visu_format_t     format,
 const char                 *format_option
);

/**
 * \brief SUser_tgt_pts_set CWIPI.
 *
 * This function must be called if the degrees of freedom locations are
 * \ref CWP_DOF_LOCATION_USER
 *
 * \param [in]  local_code_name  Local code name
 * \param [in]  cpl_id           Coupling identifier
 * \param [in]  i_part           Current partition
 * \param [in]  n_pts            Number of points
 * \param [in]  coord            Coordinates (size = 3 * n_pts)
 * \param [in]  global_num       Global number or NUL (size = n_pts)
 *
 */

void
CWP_client_User_tgt_pts_set
(
 const char    *local_code_name,
 const char    *cpl_id,
 const int      i_part,
 const int      n_pts,
 double         coord[],
 CWP_g_num_t    global_num[]
);

/**
 * \brief Mesh_interf_finalize CWIPI.
 *
 * This function computes the global numbers of mesh entities if they are
 * not provided.
 *
 * \param [in]  local_code_name  Local code name
 * \param [in]  cpl_id           Coupling identifier
 *
 */

void
CWP_client_Mesh_interf_finalize
(
 const char           *local_code_name,
 const char           *cpl_id
);


/**
 * \brief Mesh_interf_vtx_set CWIPI.
 *
 * \param [in]  local_code_name  Local code name
 * \param [in]  cpl_id           Coupling identifier
 * \param [in]  i_part           Current partition
 * \param [in]  n_pts            Number of points
 * \param [in]  coord            Coordinates (size = 3 * \p n_pts)
 * \param [in]  global_num       Pointer to parent element number (or NULL)
 *
 */

void
CWP_client_Mesh_interf_vtx_set
(
 const char           *local_code_name,
 const char           *cpl_id,
 const int             i_part,
 const int             n_pts,
 double                coord[],
 CWP_g_num_t           global_num[]
);


/**
 * \brief Mesh_interf_block_add CWIPI.
 *
 * \param [in]  local_code_name  Local code name
 * \param [in]  cpl_id           Coupling identifier
 * \param [in]  block_type       Block type
 *
 * \return block identifier
 */

int
CWP_client_Mesh_interf_block_add
(
 const char           *local_code_name,
 const char           *cpl_id,
 const CWP_Block_t     block_type
);


/**
 * \brief Mesh_interf_block_std_set CWIPI.
 *
 * \param [in]  local_code_name  Local code name
 * \param [in]  cpl_id           Coupling identifier
 * \param [in]  i_part           Partition identifier
 * \param [in]  block_id         Block identifier
 * \param [in]  n_elts           Number of elements
 * \param [in]  connec           Connectivity (size = n_vertex_elt * n_elts)
 * \param [in]  global_num       Pointer to global element number (or NULL)
 *
 */

void
CWP_client_Mesh_interf_block_std_set
(
 const char        *local_code_name,
 const char        *cpl_id,
 const int          i_part,
 const int          block_id,
 const int          n_elts,
 int                connec[],
 CWP_g_num_t        global_num[]
);


/**
 * \brief Mesh_interf_block_std_get CWIPI.
 *
 * /!\ Pointers given here are freed by `CWP_client_Cpl_del`
 *
 * \param [in]   local_code_name  Local code name
 * \param [in]   cpl_id           Coupling identifier
 * \param [in]   i_part           Partition identifier
 * \param [in]   block_id         Block identifier
 * \param [out]  n_elts           Number of elements
 * \param [out]  connec           Connectivity (size = n_vertex_elt * n_elts)
 * \param [out]  global_num       Pointer to global element number (or NULL)
 */

void
CWP_client_Mesh_interf_block_std_get
(
 const char        *local_code_name,
 const char        *cpl_id,
 const int          i_part,
 const int          block_id,
 int               *n_elts,
 int              **connec,
 CWP_g_num_t      **global_num
);

/**
  * \brief Get the standard block type
  *
  * \param [in]  local_code_name  Local code name
  * \param [in]  cpl_id           Coupling identifier
  * \param [in]  block_id    Block identifier
  *
  * \return block type
  */

CWP_Block_t
CWP_client_std_block_type_get
(
 const char             *local_code_name,
 const char             *cpl_id,
 const int               block_id
);

/**
 * \brief Mesh_interf_f_poly_block_set CWIPI.
 *
 * \param [in]  local_code_name  Local code name
 * \param [in]  cpl_id           Coupling identifier
 * \param [in]  i_part           Current partition
 * \param [in]  block_id         Block identifier
 * \param [in]  n_elts           Number of elements
 * \param [in]  connec_idx       Connectivity index (\p connec_idx[0] = 0 and
 *                               size = \p n_elts + 1)
 * \param [in]  connec           Connectivity (size = \p connec_idx[\p n_elts])
 * \param [in]  global_num       Pointer to global element number (or NULL)
 *
 */

void
CWP_client_Mesh_interf_f_poly_block_set
(
 const char             *local_code_name,
 const char             *cpl_id,
 const int               i_part,
 const int               block_id,
 const int               n_elts,
 int                     connec_idx[],
 int                     connec[],
 CWP_g_num_t             global_num[]
);


/**
 * \brief Mesh_interf_f_poly_block_get CWIPI.
 *
 * /!\ Pointers given here are freed by `CWP_client_Cpl_del`
 *
 * \param [in]  local_code_name  Local code name
 * \param [in]  cpl_id           Coupling identifier
 * \param [in]  i_part           Current partition
 * \param [in]  block_id         Block identifier
 * \param [out]  n_elts           Number of elements
 * \param [out]  connec_idx       Connectivity index (\p connec_idx[0] = 0 and
 *                               size = \p n_elts + 1)
 * \param [out]  connec           Connectivity (size = \p connec_idx[\p n_elts])
 * \param [out]  global_num       Pointer to global element number (or NULL)
 *
 */

void
CWP_client_Mesh_interf_f_poly_block_get
(
 const char             *local_code_name,
 const char             *cpl_id,
 const int               i_part,
 const int               block_id,
 int                    *n_elts,
 int                   **connec_idx,
 int                   **connec,
 CWP_g_num_t           **global_num
);


/**
 * \brief Mesh_interf_c_poly_block_set CWIPI.
 *
 * \param [in]  local_code_name   Local code name
 * \param [in]  cpl_id            Coupling identifier
 * \param [in]  i_part            Current partition
 * \param [in]  block_id          Block identifier
 * \param [in]  n_elts            Number of elements
 * \param [in]  connec_cells_idx  Polyhedron to face index
 *                                (\p src_poly_cell_face_idx[0] = 0 and
 *                                 size = \p n_elts + 1)
 * \param [in]  connec_cells      Polyhedron to face connectivity
 *                                (size = \p cell_face_idx[\p n_elts])
 * \param [in]  n_faces           Number of faces
 * \param [in]  connec_faces_idx  Polyhedron face to vertex index
 *                                (\p face_vertex_idx[0] = 0 and
 *                                 size = max(\p cell_face_connec) + 1)
 * \param [in]  connec_faces      Polyhedron face to vertex connectivity
 *                                (size = \p face_vertex_idx[\p n_elts])
 * \param [in]  global_num        Pointer to global element number (or NULL)
 *
 */

void
CWP_client_Mesh_interf_c_poly_block_set
(
 const char           *local_code_name,
 const char           *cpl_id,
 const int             i_part,
 const int             block_id,
 const int             n_elts,
 const int             n_faces,
 int                   connec_faces_idx[],
 int                   connec_faces[],
 int                   connec_cells_idx[],
 int                   connec_cells[],
 CWP_g_num_t           global_num[]
);


/**
 * \brief Mesh_interf_c_poly_block_get CWIPI.
 *
 * /!\ Pointers given here are freed by `CWP_client_Cpl_del`
 *
 * \param [in]  local_code_name   Local code name
 * \param [in]  cpl_id            Coupling identifier
 * \param [in]  i_part            Current partition
 * \param [in]  block_id          Block identifier
 * \param [out]  n_elts            Number of elements
 * \param [out]  connec_cells_idx  Polyhedron to face index
 *                                (\p src_poly_cell_face_idx[0] = 0 and
 *                                 size = \p n_elts + 1)
 * \param [out]  connec_cells      Polyhedron to face connectivity
 *                                (size = \p cell_face_idx[\p n_elts])
 * \param [out]  n_faces           Number of faces
 * \param [out]  connec_faces_idx  Polyhedron face to vertex index
 *                                (\p face_vertex_idx[0] = 0 and
 *                                 size = max(\p cell_face_connec) + 1)
 * \param [out]  connec_faces      Polyhedron face to vertex connectivity
 *                                (size = \p face_vertex_idx[\p n_elts])
 * \param [out]  global_num        Pointer to global element number (or NULL)
 *
 */

void
CWP_client_Mesh_interf_c_poly_block_get
(
 const char           *local_code_name,
 const char           *cpl_id,
 const int             i_part,
 const int             block_id,
 int                  *n_elts,
 int                  *n_faces,
 int                 **connec_faces_idx,
 int                 **connec_faces,
 int                 **connec_cells_idx,
 int                 **connec_cells,
 CWP_g_num_t         **global_num
);


/**
 * \brief Mesh_interf_del CWIPI.
 *
 * \param [in] local_code_name  Local code name
 * \param [in] cpl_id           Coupling identifier
 *
 */

void
CWP_client_Mesh_interf_del
(
 const char *local_code_name,
 const char *cpl_id
);


/**
 * \brief Mesh_interf_from_cellface_set CWIPI.
 *
 * \param [in]  local_code_name   Local code name
 * \param [in]  cpl_id            Coupling identifier
 * \param [in]  i_part            Current partition
 * \param [in]  n_cells           Number of cells
 * \param [in]  cell_face_idx     Polyhedron to face index
 *                                (\p src_poly_cell_face_idx[0] = 0 and
 *                                 size = \p n_elts + 1)
 * \param [in]  cell_face         Cell to face connectivity
 *                                (size = \p cell_face_idx[\p n_elts])
 * \param [in]  n_faces           Number of faces
 * \param [in]  face_vtx_idx      Polyhedron face to vertex index
 *                                (\p face_vtx_idx[0] = 0 and
 *                                 size = \p n_faces + 1)
 * \param [in]  face_vtx          Face to vertex connectivity
 *                                (size = \p face_vtx_idx[\p n_elts])
 * \param [in]  global_num        Pointer to parent element number (or NULL)
 *
 */

void
CWP_client_Mesh_interf_from_cellface_set
(
 const char           *local_code_name,
 const char           *cpl_id,
 const int             i_part,
 const int             n_cells,
 int                   cell_face_idx[],
 int                   cell_face[],
 const int             n_faces,
 int                   face_vtx_idx[],
 int                   face_vtx[],
 CWP_g_num_t           global_num[]
);


/**
 * \brief Mesh_interf_from_faceedge_set CWIPI.
 *
 * \param [in]  local_code_name   Local code name
 * \param [in]  cpl_id            Coupling identifier
 * \param [in]  i_part            Current partition
 * \param [in]  n_faces           Number of cells
 * \param [in]  face_edge_idx     Polygon to edge index
 *                                (\p face_edge_idx[0] = 0 and
 *                                 size =  \p n_faces + 1)
 * \param [in]  face_edge         Face to edge connectivity
 *                                (size = \p face_edge_idx[\p n_faces])
 * \param [in]  n_edges           Number of faces
 * \param [in]  edge_vtx          Edge to vertex connectivity
 *                                (size = 2 * \p n_edges)
 * \param [in]  global_num        Pointer to parent element number (or NULL)
 *
 */

void
CWP_client_Mesh_interf_from_faceedge_set
(
 const char           *local_code_name,
 const char           *cpl_id,
 const int             i_part,
 const int             n_faces,
       int             face_edge_idx[],
       int             face_edge[],
 const int             n_edges,
       int             edge_vtx[],
       CWP_g_num_t     global_num[]
);

/**
 *
 * \brief Field_create CWIPI.
 *
 * \param [in]  local_code_name Local code name
 * \param [in]  cpl_id          Coupling identifier
 * \param [in]  field_id        Field id
 * \param [in]  data_type       Data type
 * \param [in]  storage         Storage type
 * \param [in]  n_component     Number of component
 * \param [in]  target_location Target location
 * \param [in]  exch_type       Exchange type
 * \param [in]  visu_status     Visualization status
 *
 */

void
CWP_client_Field_create
(
 const char                  *local_code_name,
 const char                  *cpl_id,
 const char                  *field_id,
 const CWP_Type_t             data_type,
 const CWP_Field_storage_t    storage,
 const int                    n_component,
 const CWP_Dof_location_t     target_location,
 const CWP_Field_exch_t       exch_type,
 const CWP_Status_t           visu_status
);


/**
 *
 * \brief Field_data_set CWIPI.
 *
 * \param [in] local_code_name   Local code name
 * \param [in] cpl_id            Coupling identifier
 * \param [in] field_id          Field identifier
 * \param [in] i_part            Current partition
 * \param [in] map_type          Choice if data is setted for the source or the target
 * \param [in] data              Storage array (Mapping)
 *
 */

void
CWP_client_Field_data_set
(
 const char              *local_code_name,
 const char              *cpl_id,
 const char              *field_id,
 const int                i_part,
 const CWP_Field_map_t    map_type,
 double                  *data
);

/**
 *
 * \brief Field_dof_location_get CWIPI.
 *
 * \param [in] local_code_name  Local code name
 * \param [in] cpl_id           Coupling identifier
 * \param [in] field_id         Field identifier
 *
 * \return                      Location of degrees of freedom
 *
 */

CWP_Dof_location_t
CWP_client_Field_dof_location_get
(
 const char      *local_code_name,
 const char      *cpl_id,
 const char      *field_id
);

/**
 *
 * \brief Field_storage_get CWIPI.
 *
 * \param [in] local_code_name  Local code name
 * \param [in] cpl_id           Coupling identifier
 * \param [in] field_id         Field identifier
 *
 * \return                      Field storage type
 */

CWP_Field_storage_t
CWP_client_Field_storage_get
(
 const char      *local_code_name,
 const char      *cpl_id         ,
 const char      *field_id
);

/**
 *
 * \brief Field_n_components_get CWIPI.
 *
 * \param [in]  local_code_name           Local code name
 * \param [in]  cpl_id                    Coupling identifier
 * \param [in]  field_id                  Field identifier
 *
 */

int
CWP_client_Field_n_components_get
(
 const char             *local_code_name,
 const char             *cpl_id,
 const char             *field_id
);

/**
 *
 * \brief Field_n_dof_get CWIPI.
 *
 * \param [in]  local_code_name           Local code name
 * \param [in]  cpl_id                    Coupling identifier
 * \param [in]  field_id                  Field identifier
 * \param [in]  i_part                    Partition identifier
 *
 */

int
CWP_client_Field_n_dof_get
(
 const char             *local_code_name,
 const char             *cpl_id,
 const char             *field_id,
 int                     i_part
);

/**
 * \brief Field_del CWIPI.
 *
 * \param [in] local_code_name Local code name
 * \param [in]  cpl_id         Coupling identifier
 * \param [in]  field_id       Field identifier
 *
 */

void
CWP_client_Field_del
(
 const char      *local_code_name,
 const char      *cpl_id         ,
 const char      *field_id
);

/**
 * \brief Field_issend CWIPI.
 *
 * This function is independant of \ref CWP_Time_exch_t mode. The user has to
 * manually check the consistency of the exchanges.
 *
 * \param [in]  local_code_name  Local code name
 * \param [in]  cpl_id           Coupling identifier
 * \param [in]  field_id         Field identifier
 *
 *
 */

void
CWP_client_Field_issend
(
 const char     *local_code_name,
 const char     *cpl_id,
 const char     *field_id
);

/**
 *
 * \brief Field_irecv CWIPI.
 *
 * This function is independant of \ref CWP_Time_exch_t mode. The user has to
 * manually check the consistency of the exchanges.
 *
 * \param [in] local_code_name  Local code name
 * \param [in]  cpl_id          Coupling identifier
 * \param [in]  tgt_field_id    Target field id
 *
 *
 */

void
CWP_client_Field_irecv
(
 const char        *local_code_name,
 const char        *cpl_id,
 const char        *tgt_field_id
);

/**
 *
 * \brief Field_wait_issend CWIPI.
 *
 * \param [in] local_code_name  Local code name
 * \param [in] cpl_id           Coupling identifier
 * \param [in] field_id         Field identifier
 *
 */

void
CWP_client_Field_wait_issend
(
 const char  *local_code_name,
 const char  *cpl_id,
 const char  *field_id
);


/**
 *
 * \brief Field_wait_irecv CWIPI.
 *
 * This function waits the end of exchange related to request
 * from \ref CWP_Field_irecv
 *
 * \param [in] local_code_name  Local code name
 * \param [in] cpl_id           Coupling identifier
 * \param [in] tgt_field_id     Target field id
 *
 */

void
CWP_client_Field_wait_irecv
(
 const char              *local_code_name,
 const char              *cpl_id,
 const char              *tgt_field_id
);

/**
 *
 * \brief Interp_function_unset CWIPI.
 *
 * \param [in] local_code_name  Local code name
 * \param [in] cpl_id           Coupling identifier
 * \param [in] field_id         Field identifier
 *
 */

void
CWP_client_Interp_function_unset
(
 const char                 *local_code_name,
 const char                 *cpl_id,
 const char                 *field_id
);


/**
 *
 * \brief Interp_function_set CWIPI.
 *
 * This function takes into account an user interpolation function written with
 * void (*\ref CWP_Interp_function_t) interface.
 *
 * \param [in] local_code_name  Local code name
 * \param [in] cpl_id           Coupling identifier
 * \param [in] field_id         Field identifier
 * \param [in] fct              Function
 *
 */

void
CWP_client_Interp_function_set
(
 const char                 *local_code_name,
 const char                 *cpl_id,
 const char                 *field_id,
 CWP_Interp_function_t       fct
);

// NB: no need to do getters since this mode does not exist in client-server mode

/**
 * \brief Mesh_interf_block_ho_set CWIPI.
 *
 * \param [in]  local_code_name  Local code name
 * \param [in]  cpl_id           Coupling identifier
 * \param [in]  i_part           Partition identifier
 * \param [in]  block_id         Block identifier
 * \param [in]  n_elts           Number of elements
 * \param [in]  order            Element order
 * \param [in]  connec           Connectivity (size = n_vertex_elt * n_elts)
 * \param [in]  global_num       Pointer to global element number (or NULL)
 *
 */

void
CWP_client_Mesh_interf_block_ho_set
(
 const char        *local_code_name,
 const char        *cpl_id,
 const int          i_part,
 const int          block_id,
 const int          n_elts,
 const int          order,
 int                connec[],
 CWP_g_num_t        global_num[]
);


/**
 * \brief Mesh_interf_block_ho_get CWIPI.
 *
 * \param [in]  local_code_name  Local code name
 * \param [in]  cpl_id           Coupling identifier
 * \param [in]  i_part           Partition identifier
 * \param [in]  block_id         Block identifier
 * \param [out]  n_elts           Number of elements
 * \param [out]  order            Element order
 * \param [out]  connec           Connectivity (size = n_vertex_elt * n_elts)
 * \param [out]  global_num       Pointer to global element number (or NULL)
 *
 */

void
CWP_client_Mesh_interf_block_ho_get
(
 const char        *local_code_name,
 const char        *cpl_id,
 const int          i_part,
 const int          block_id,
 int               *n_elts,
 int               *order,
 int              **connec,
 CWP_g_num_t      **global_num
);

/**
 *
 * \brief Mesh_interf_ho_ordering_from_IJK_set CWIPI.
 *
 * \param [in]  local_code_name  Local code name
 * \param [in]  cpl_id           Coupling identifier
 * \param [in]  block_type       Block type
 * \param [in]  order            Element order
 * \param [in]  n_nodes          Number of nodes
 * \param [in]  ijk_grid         User ordering to (u, v, w) grid (size = elt_dim * n_nodes)
 *
 */

void
CWP_client_Mesh_interf_ho_ordering_from_IJK_set
(
 const char        *local_code_name,
 const char        *cpl_id,
 const CWP_Block_t  block_type,
 const int          order,
 const int          n_nodes,
 const int         *ijk_grid
);

/**
 *
 * \brief Cpl_spatial_interp_algo_get CWIPI.
 *
 * \param [in] local_code_name  Local code name
 * \param [in] cpl_id           Coupling identifier
 *
 * \return                      Spatial interpolation algorithm
 */

CWP_Spatial_interp_t
CWP_client_Cpl_spatial_interp_algo_get
(
 const char *local_code_name,
 const char *cpl_id
);

/**
 * \brief Global_data_issend CWIPI.
 *
 * \param [in] local_code_name  Local code name
 * \param [in] cpl_id           Coupling identifier
 * \param [in] global_data_id
 * \param [in] s_send_entity
 * \param [in] send_stride
 * \param [in] n_send_entity
 * \param [in] send_data
 *
 */

void
CWP_client_Global_data_issend
(
 const char     *local_code_name,
 const char     *cpl_id,
 const char     *global_data_id,
 size_t          s_send_entity,
 int             send_stride,
 int             n_send_entity,
 void           *send_data
);

/**
 * \brief Global_data_irecv CWIPI.
 *
 * \param [in]  local_code_name  Local code name
 * \param [in]  cpl_id           Coupling identifier
 * \param [in]  global_data_id   Global data identifier
 * \param [in]  s_recv_entity    Data size
 * \param [in]  recv_stride      Constant stride value
 * \param [in]  n_recv_entity    Number of entities
 * \param [out] recv_data        Pointer to data array
 *
 */

void
CWP_client_Global_data_irecv
(
 const char     *local_code_name,
 const char     *cpl_id,
 const char     *global_data_id,
       size_t    s_recv_entity,
       int       recv_stride,
       int       n_recv_entity,
       void     *recv_data
);

/**
 * \brief Global_data_wait_issend CWIPI.
 *
 * \param [in] local_code_name  Local code name
 * \param [in] cpl_id           Coupling identifier
 * \param [in] global_data_id
 *
 */

void
CWP_client_Global_data_wait_issend
(
 const char     *local_code_name,
 const char     *cpl_id,
 const char     *global_data_id
);

/**
 * \brief Global_data_wait_irecv CWIPI.
 *
 * \param [in]  local_code_name  Local code name
 * \param [in]  cpl_id           Coupling identifier
 * \param [in]  global_data_id   Global data identifier
 *
 */

void
CWP_client_Global_data_wait_irecv
(
 const char     *local_code_name,
 const char     *cpl_id,
 const char     *global_data_id
);

/**
 * \brief Part_data_create CWIPI.
 *
 * \param [in] local_code_name  Local code name
 * \param [in] cpl_id           Coupling identifier
 * \param [in] part_data_id     Part data identifier
 * \param [in] exch_type        Send or receive
 * \param [in] gnum_elt         Element global number array
 * \param [in] n_elt            Number of elements per partition
 * \param [in] n_part           Number of partitions
 *
 */

void
CWP_client_Part_data_create
(
 const char           *local_code_name,
 const char           *cpl_id,
 const char           *part_data_id,
 CWP_PartData_exch_t   exch_type,
 CWP_g_num_t         **gnum_elt,
 int                  *n_elt,
 int                   n_part
 );

/**
 * \brief Part_data_del CWIPI.
 *
 * \param [in] local_code_name  Local code name
 * \param [in] cpl_id           Coupling identifier
 * \param [in] part_data_id     Part data identifier
 *
 */

void
CWP_client_Part_data_del
(
 const char          *local_code_name,
 const char          *cpl_id,
 const char          *part_data_id
);


/**
 * \brief Part_data_issend CWIPI.
 *
 * \param [in] local_code_name     Local code name
 * \param [in] cpl_id              Coupling identifier
 * \param [in] part_data_id        Partitioned data identifier
 * \param [in] exch_id             Exchange identifier
 * \param [in] s_data              Data size
 * \param [in] n_components        Number of components
 * \param [in] send_data           Pointer to send data
 *
 */

void
CWP_client_Part_data_issend
(
 const char    *local_code_name,
 const char    *cpl_id,
 const char    *part_data_id,
 const int      exch_id,
       size_t   s_data,
       int      n_components,
       void   **send_data
);

/**
 * \brief Part_data_irecv CWIPI.
 *
 * \param [in] local_code_name     Local code name
 * \param [in] cpl_id              Coupling identifier
 * \param [in] part_data_id        Partitioned data identifier
 * \param [in] exch_id             Exchange identifier
 * \param [in] s_data              Data size
 * \param [in] n_components        Number of components
 * \param [in] recv_data           Pointer to received data
 *
 */

void
CWP_client_Part_data_irecv
(
 const char    *local_code_name,
 const char    *cpl_id,
 const char    *part_data_id,
 const int      exch_id,
       size_t   s_data,
       int      n_components,
       void   **recv_data
);

/**
 * \brief Part_data_wait_issend CWIPI.
 *
 * \param [in] local_code_name  Local code name
 * \param [in] cpl_id           Coupling identifier
 * \param [in] part_data_id     Partitioned data identifier
 * \param [in] exch_id          Exchange identifier
 *
 */

void
CWP_client_Part_data_wait_issend
(
 const char    *local_code_name,
 const char    *cpl_id,
 const char    *part_data_id,
 const int      exch_id
);

/**
 * \brief Part_data_wait_irecv CWIPI.
 *
 * \param [in] local_code_name  Local code name
 * \param [in] cpl_id           Coupling identifier
 * \param [in] part_data_id     Partitioned data identifier
 * \param [in] exch_id          Exchange identifier
 *
 */

void
CWP_client_Part_data_wait_irecv
(
 const char    *local_code_name,
 const char    *cpl_id,
 const char    *part_data_id,
 const int      exch_id
);

/*=============================================================================
 * Client function interfaces
 *============================================================================*/

/* Connect to a server */

int
CWP_client_connect
(
 MPI_Comm  comm,
 const char* server_name,
 int server_port,
 int flags
);

/* Disconnect */

int
CWP_client_disconnect
(
 void
);

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* __CLIENT_H__ */
