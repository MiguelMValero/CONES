#ifndef __SERVER_H__
#define __SERVER_H__
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
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <mpi.h>
#include "cwp.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "message.h"

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force back Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */

/*=============================================================================
 * Macro definitions
 *============================================================================*/
/* server status */
#define CWP_SVRSTATE_WAITCONN      0
#define CWP_SVRSTATE_LISTENINGMSG  1
#define CWP_SVRSTATE_RECVPPUTDATA  2
#define CWP_SVRSTATE_SENDPGETDATA  3
#define CWP_SVRSTATE_TERMINATING   4

/*============================================================================
 * Types definition
 *============================================================================*/

typedef struct t_server
{
  uint16_t     port;
  int          state;
  int          flags;
  int          max_msg_size;
  int          listen_socket;
  int          connected_socket;
  int          client_endianess;
  int          server_endianess;
  char         host_name[256];
}t_server,*p_server;

/*=============================================================================
 * Server CWIPI function interfaces
 *============================================================================*/

/**
 * \brief Initialize CWIPI.
 *
 * \param [in]  p_server       Pointer on server data structure
 *
 */

void
CWP_server_Init
(
  p_server                 svr
);

/**
 *
 * \brief Finalize CWIPI.
 *
 */

void
CWP_server_Finalize
(
 p_server                 svr
);

/**
 *
 * \brief Param_lock CWIPI.
 *
 * \param [in]  p_server       Pointer on server data structure
 *
 */

void
CWP_server_Param_lock
(
 p_server                 svr
);

/**
 *
 * \brief Param_unlock CWIPI.
 *
 * \param [in]  p_server       Pointer on server data structure
 *
 */

void
CWP_server_Param_unlock
(
 p_server                 svr
);

/**
 *
 * \brief Param_add CWIPI.
 *
 * \param [in]  p_server       Pointer on server data structure
 *
 */

void
CWP_server_Param_add
(
 p_server                 svr
);

/**
 *
 * \brief Param_get CWIPI.
 *
 * \param [in]  p_server       Pointer on server data structure
 *
 */

void
CWP_server_Param_get
(
 p_server                 svr
);

/**
 *
 * \brief Param_set CWIPI.
 *
 * \param [in]  p_server       Pointer on server data structure
 *
 */

void
CWP_server_Param_set
(
 p_server                 svr
);


/**
 *
 * \brief Param_del CWIPI.
 *
 * \param [in]  p_server       Pointer on server data structure
 *
 */

void
CWP_server_Param_del
(
 p_server                 svr
);


/**
 *
 * \brief Param_n_get CWIPI.
 *
 * \param [in]  p_server       Pointer on server data structure
 *
 */

void
CWP_server_Param_n_get
(
 p_server                 svr
);


/**
 *
 * \brief Param_list_get CWIPI.
 *
 * \param [in]  p_server       Pointer on server data structure
 *
 */

void
CWP_server_Param_list_get
(
 p_server                 svr
);


/**
 *
 * \brief Param_is CWIPI.
 *
 * \param [in]  p_server       Pointer on server data structure
 *
 */

void
CWP_server_Param_is
(
 p_server                 svr
);


/**
 *
 * \brief Param_reduce CWIPI.
 *
 * \param [in]  p_server       Pointer on server data structure
 *
 */

void
CWP_server_Param_reduce
(
 p_server                 svr
);


/**
 *
 * \brief Cpl_create CWIPI.
 *
 * \param [in]  p_server       Pointer on server data structure
 *
 */

void
CWP_server_Cpl_create
(
 p_server                 svr
);

/**
 *
 * \brief CWP_Cpl_barrier CWIPI.
 *
 * \param [in]  p_server       Pointer on server data structure
 *
 */

void
CWP_server_Cpl_barrier
(
 p_server                 svr
);

/**
 *
 * \brief Cpl_del CWIPI.
 *
 * \param [in]  p_server       Pointer on server data structure
 *
 */

void
CWP_server_Cpl_del
(
 p_server                 svr
);

/**
 * \brief State_update CWIPI.
 *
 * \param [in]  p_server       Pointer on server data structure
 *
 */

void
CWP_server_State_update
(
  p_server                 svr
);

/**
 * \brief Time_step_beg CWIPI.
 *
 * \param [in]  p_server       Pointer on server data structure
 *
 */

void
CWP_server_Time_step_beg
(
  p_server                 svr
);

/**
 * \brief Time_step_end CWIPI.
 *
 * \param [in]  p_server       Pointer on server data structure
 *
 */

void
CWP_server_Time_step_end
(
  p_server                 svr
);

/**
 * \brief Output_file_set CWIPI.
 *
 * \param [in]  p_server       Pointer on server data structure
 *
 */

void
CWP_server_Output_file_set
(
  p_server                 svr
);


/**
 * \brief User_structure_set CWIPI.
 *
 * \param [in]  p_server       Pointer on server data structure
 *
 */

void
CWP_server_User_structure_set
(
  p_server                 svr
);


/**
 * \brief User_structure_get CWIPI.
 *
 * \param [in]  p_server       Pointer on server data structure
 *
 */

void
CWP_server_User_structure_get
(
  p_server                 svr
);

/**
 * \brief State_get CWIPI.
 *
 * \param [in]  p_server       Pointer on server data structure
 *
 */

void
CWP_server_State_get
(
  p_server                 svr
);


/**
 * \brief Codes_nb_get CWIPI.
 *
 * \param [in]  p_server       Pointer on server data structure
 *
 */

void
CWP_server_Codes_nb_get
(
  p_server                 svr
);

/**
 * \brief Codes_list_get CWIPI.
 *
 * \param [in]  p_server       Pointer on server data structure
 *
 */

void
CWP_server_Codes_list_get
(
  p_server                 svr
);


/**
 * \brief Loc_codes_nb_get CWIPI.
 *
 * \param [in]  p_server       Pointer on server data structure
 *
 */

void
CWP_server_Loc_codes_nb_get
(
  p_server                 svr
);


/**
 * \brief Loc_codes_list_get CWIPI.
 *
 * \param [in]  p_server       Pointer on server data structure
 *
 */

void
CWP_server_Loc_codes_list_get
(
  p_server                 svr
);

/**
 * \brief Properties_dump CWIPI.
 *
 * \param [in]  p_server       Pointer on server data structure
 *
 */

void
CWP_server_Properties_dump
(
  p_server                 svr
);

/**
 *
 * \brief Computed_tgts_bcast_enable CWIPI.
 *
 * \param [in]  p_server       Pointer on server data structure
 */

void
CWP_server_Computed_tgts_bcast_enable
(
  p_server                 svr
);

/**
 *
 * \brief Involved_srcs_bcast_enable CWIPI.
 *
 * \param [in]  p_server       Pointer on server data structure
 */

void
CWP_server_Involved_srcs_bcast_enable
(
  p_server                 svr
);

/**
 *
 * \brief N_uncomputed_tgts_get CWIPI.
 *
 * \param [in]  p_server       Pointer on server data structure
 *
 */

void
CWP_server_N_uncomputed_tgts_get
(
  p_server                 svr
);


/**
 *
 * \brief Uncomputed_tgts_get CWIPI.
 *
 * \param [in]  p_server       Pointer on server data structure
 *
 */

void
CWP_server_Uncomputed_tgts_get
(
  p_server                 svr
);

/**
 *
 * \brief N_computed_tgts_get CWIPI.
 *
 * \param [in]  p_server       Pointer on server data structure
 *
 */

void
CWP_server_N_computed_tgts_get
(
  p_server                 svr
);

/**
 *
 * \brief Computed_tgts_get CWIPI.
 *
 * \param [in]  p_server       Pointer on server data structure
 *
 */

void
CWP_server_Computed_tgts_get
(
  p_server                 svr
);


/**
 *
 * \brief N_involved_srcs_get CWIPI.
 *
 * \param [in]  p_server       Pointer on server data structure
 *
 */

void
CWP_server_N_involved_srcs_get
(
  p_server                 svr
);


/**
 *
 * \brief Involved_srcs_get CWIPI.
 *
 * \param [in]  p_server       Pointer on server data structure
 *
 */

void
CWP_server_Involved_srcs_get
(
  p_server                 svr
);

/**
 * \brief Spatial_interp_weights_compute CWIPI.
 *
 * \param [in]  p_server       Pointer on server data structure
 *
 */

void
CWP_server_Spatial_interp_weights_compute
(
  p_server                 svr
);

/**
 * \brief Spatial_interp_property_set CWIPI.
 *
 * \param [in]  p_server       Pointer on server data structure
 *
 */

void
CWP_server_Spatial_interp_property_set
(
  p_server                 svr
);

/**
 * \brief Visu_set CWIPI.
 *
 * \param [in]  p_server       Pointer on server data structure
 *
 */

void
CWP_server_Visu_set
(
  p_server                 svr
);

/**
 * \brief SUser_tgt_pts_set CWIPI.
 *
 * \param [in]  p_server       Pointer on server data structure
 *
 */

void
CWP_server_User_tgt_pts_set
(
  p_server                 svr
);

/**
 * \brief Mesh_interf_finalize CWIPI.
 *
 * \param [in]  p_server       Pointer on server data structure
 *
 */

void
CWP_server_Mesh_interf_finalize
(
  p_server                 svr
);


/**
 * \brief Mesh_interf_vtx_set CWIPI.
 *
 * \param [in]  p_server       Pointer on server data structure
 *
 */

void
CWP_server_Mesh_interf_vtx_set
(
  p_server                 svr
);


/**
 * \brief Mesh_interf_block_add CWIPI.
 *
 * \param [in]  p_server       Pointer on server data structure
 *
 */

void
CWP_server_Mesh_interf_block_add
(
  p_server                 svr
);


/**
 * \brief Mesh_interf_block_std_set CWIPI.
 *
 * \param [in]  p_server       Pointer on server data structure
 *
 */

void
CWP_server_Mesh_interf_block_std_set
(
  p_server                 svr
);


/**
 * \brief Mesh_interf_block_std_get CWIPI.
 *
 * \param [in]  p_server       Pointer on server data structure
 *
 */

void
CWP_server_Mesh_interf_block_std_get
(
  p_server                 svr
);

/**
 * \brief Get the standard block type.
 *
 * \param [in]  p_server       Pointer on server data structure
 *
 */

void
CWP_server_std_block_type_get
(
 p_server                 svr
);

/**
 * \brief Mesh_interf_f_poly_block_set CWIPI.
 *
 * \param [in]  p_server       Pointer on server data structure
 *
 */

void
CWP_server_Mesh_interf_f_poly_block_set
(
  p_server                 svr
);


/**
 * \brief Mesh_interf_f_poly_block_get CWIPI.
 *
 * \param [in]  p_server       Pointer on server data structure
 *
 */

void
CWP_server_Mesh_interf_f_poly_block_get
(
  p_server                 svr
);


/**
 * \brief Mesh_interf_c_poly_block_set CWIPI.
 *
 * \param [in]  p_server       Pointer on server data structure
 *
 */

void
CWP_server_Mesh_interf_c_poly_block_set
(
  p_server                 svr
);


/**
 * \brief Mesh_interf_c_poly_block_get CWIPI.
 *
 * \param [in]  p_server       Pointer on server data structure
 *
 */

void
CWP_server_Mesh_interf_c_poly_block_get
(
  p_server                 svr
);


/**
 * \brief Mesh_interf_del CWIPI.
 *
 * \param [in]  p_server       Pointer on server data structure
 *
 */

void
CWP_server_Mesh_interf_del
(
  p_server                 svr
);


/**
 * \brief Mesh_interf_from_cellface_set CWIPI.
 *
 * \param [in]  p_server       Pointer on server data structure
 *
 */

void
CWP_server_Mesh_interf_from_cellface_set
(
  p_server                 svr
);


/**
 * \brief Mesh_interf_from_faceedge_set CWIPI.
 *
 * \param [in]  p_server       Pointer on server data structure
 *
 */

void
CWP_server_Mesh_interf_from_faceedge_set
(
  p_server                 svr
);

/**
 *
 * \brief Field_create CWIPI.
 *
 * \param [in]  p_server       Pointer on server data structure
 *
 */

void
CWP_server_Field_create
(
  p_server                 svr
);


/**
 *
 * \brief Field_data_set CWIPI.
 *
 * \param [in]  p_server       Pointer on server data structure
 *
 */

void
CWP_server_Field_data_set
(
  p_server                 svr
);

/**
 *
 * \brief Field_n_component_get CWIPI.
 *
 * \param [in]  p_server       Pointer on server data structure
 *
 */

void
CWP_server_Field_n_component_get
(
  p_server                 svr
);

/**
 *
 * \brief Field_dof_location_get CWIPI.
 *
 * \param [in]  p_server       Pointer on server data structure
 *
 */

void
CWP_server_Field_dof_location_get
(
  p_server                 svr
);

/**
 *
 * \brief Field_storage_get CWIPI.
 *
 * \param [in]  p_server       Pointer on server data structure
 *
 */

void
CWP_server_Field_storage_get
(
  p_server                 svr
);

/**
 * \brief Field_n_components_get CWIPI.
 *
 * \param [in]  p_server       Pointer on server data structure
 *
 */

void
CWP_server_Field_n_components_get
(
  p_server                 svr
);

/**
 *
 * \brief Field_n_dof_get CWIPI.
 *
 * \param [in]  p_server       Pointer on server data structure
 *
 */

void
CWP_server_Field_n_dof_get
(
  p_server                 svr
);

/**
 * \brief Field_del CWIPI.
 *
 * \param [in]  p_server       Pointer on server data structure
 *
 */

void
CWP_server_Field_del
(
  p_server                 svr
);

/**
 * \brief Field_issend CWIPI.
 *
 * \param [in]  p_server       Pointer on server data structure
 *
 */

void
CWP_server_Field_issend
(
  p_server                 svr
);

/**
 *
 * \brief Field_irecv CWIPI.
 *
 * \param [in]  p_server       Pointer on server data structure
 *
 */

void
CWP_server_Field_irecv
(
  p_server                 svr
);

/**
 *
 * \brief Field_wait_issend CWIPI.
 *
 * \param [in]  p_server       Pointer on server data structure
 *
 */

void
CWP_server_Field_wait_issend
(
  p_server                 svr
);


/**
 *
 * \brief Field_wait_irecv CWIPI.
 *
 * \param [in]  p_server       Pointer on server data structure
 *
 */

void
CWP_server_Field_wait_irecv
(
  p_server                 svr
);

/**
 *
 * \brief Interp_function_unset CWIPI.
 *
 * \param [in]  p_server       Pointer on server data structure
 *
 */

void
CWP_server_Interp_function_unset
(
  p_server                 svr
);


/**
 *
 * \brief Interp_function_set CWIPI.
 *
 * \param [in]  p_server       Pointer on server data structure
 *
 */

void
CWP_server_Interp_function_set
(
  p_server                 svr
);

/**
 * \brief Mesh_interf_h_order_block_set CWIPI.
 *
 * \param [in]  p_server       Pointer on server data structure
 *
 */

void
CWP_server_Mesh_interf_block_ho_set
(
  p_server                 svr
);


/**
 * \brief Mesh_interf_block_ho_get CWIPI.
 *
 * \param [in]  p_server       Pointer on server data structure
 *
 */

void
CWP_server_Mesh_interf_block_ho_get
(
  p_server                 svr
);

/**
 *
 * \brief Mesh_interf_ho_ordering_from_IJK_set CWIPI.
 *
 * \param [in]  p_server       Pointer on server data structure
 *
 */

void
CWP_server_Mesh_interf_ho_ordering_from_IJK_set
(
  p_server                 svr
);

/**
 *
 * \brief Cpl_spatial_interp_algo_get CWIPI.
 *
 * \param [in]  p_server       Pointer on server data structure
 *
 */

void
CWP_server_Cpl_spatial_interp_algo_get
(
  p_server                 svr
);

/**
 * \brief Global_data_issend CWIPI.
 *
 * \param [in]  p_server       Pointer on server data structure
 *
 */

void
CWP_server_Global_data_issend
(
  p_server                 svr
);

/**
 * \brief Global_data_irecv CWIPI.
 *
 * \param [in]  p_server       Pointer on server data structure
 *
 */

void
CWP_server_Global_data_irecv
(
  p_server                 svr
);

/**
 * \brief Global_data_wait_issend CWIPI.
 *
 * \param [in]  p_server       Pointer on server data structure
 *
 */

void
CWP_server_Global_data_wait_issend
(
  p_server                 svr
);

/**
 * \brief Global_data_wait_irecv CWIPI.
 *
 * \param [in]  p_server       Pointer on server data structure
 *
 */

void
CWP_server_Global_data_wait_irecv
(
  p_server                 svr
);

/**
 * \brief Part_data_create CWIPI.
 *
 * \param [in]  p_server       Pointer on server data structure
 *
 */

void
CWP_server_Part_data_create
(
  p_server                 svr
);

/**
 * \brief Part_data_del CWIPI.
 *
 * \param [in]  p_server       Pointer on server data structure
 *
 */

void
CWP_server_Part_data_del
(
  p_server                 svr
);


/**
 * \brief Part_data_issend CWIPI.
 *
 * \param [in]  p_server       Pointer on server data structure
 *
 */

void
CWP_server_Part_data_issend
(
  p_server                 svr
);

/**
 * \brief Part_data_irecv CWIPI.
 *
 * \param [in]  p_server       Pointer on server data structure
 *
 */

void
CWP_server_Part_data_irecv
(
  p_server                 svr
);

/**
 * \brief Part_data_wait_issend CWIPI.
 *
 * \param [in]  p_server       Pointer on server data structure
 *
 */

void
CWP_server_Part_data_wait_issend
(
  p_server                 svr
);

/**
 * \brief Part_data_wait_irecv CWIPI.
 *
 * \param [in]  p_server       Pointer on server data structure
 *
 */

void
CWP_server_Part_data_wait_irecv
(
  p_server                 svr
);

/*=============================================================================
 * Server CWIPI function interfaces that are not implemented yet
 *============================================================================*/

/**
 * \brief Mesh_interf_h_order_block_set CWIPI. <b>(Not implemented yet)</b>
 *
 * \param [in]  p_server       Pointer on server data structure
 *
 */

void
CWP_server_Mesh_interf_h_order_block_set
(
  p_server                 svr
);


/**
 * \brief Mesh_interf_h_order_block_get CWIPI. <b>(Not implemented yet)</b>
 *
 * \param [in]  p_server       Pointer on server data structure
 *
 */

void
CWP_server_Mesh_interf_h_order_block_get
(
  p_server                 svr
);


/**
 * \brief Cpl_trans_init CWIPI. <b>(Not implemented yet)</b>
 *
 * \param [in]  p_server       Pointer on server data structure
 *
 */

void
CWP_server_Cpl_trans_init
(
  p_server                 svr
);


/**
 * \brief Cpl_trans_update CWIPI. <b>(Not implemented yet)</b>
  *
 * \param [in]  p_server       Pointer on server data structure
 *
 */

void
CWP_server_Cpl_trans_update
(
  p_server                 svr
);


/**
 * \brief Cpl_rotation_init CWIPI. <b>(Not implemented yet)</b>
 *
 * \param [in]  p_server       Pointer on server data structure
 *
 */

void
CWP_server_Cpl_rotation_init
(
  p_server                 svr
);


/**
 * \brief Cpl_rotation_update CWIPI. <b>(Not implemented yet)</b>
 *
 * \param [in]  p_server       Pointer on server data structure
 *
 */

void
CWP_server_Cpl_rotation_update
(
  p_server                 svr
);


/**
 * \brief Cpl_storage_properties_set CWIPI. <b>(Not implemented yet)</b>
 *
 * \param [in]  p_server       Pointer on server data structure
 *
 */

void
CWP_server_Cpl_storage_properties_set
(
  p_server                 svr
);


/**
 *
 * \brief Interp_from_intersect_set CWIPI. <b>(Not implemented yet)</b>
 *
 * \param [in]  p_server       Pointer on server data structure
 *
 */

void
CWP_server_Interp_from_intersect_set
(
  p_server                 svr
);

/**
 *
 * \brief Interp_from_closest_pts_set CWIPI. <b>(Not implemented yet)</b>
 *
 * \param [in]  p_server       Pointer on server data structure
 *
 */

void
CWP_server_Interp_from_closest_pts_set
(
  p_server                 svr
);



/**
 * \brief Computed_tgts_dist_to_spatial_interp_get CWIPI. <b>(Not implemented yet)</b>
 *
 * \param [in]  p_server       Pointer on server data structure
 *
 */

void
CWP_server_Computed_tgts_dist_to_spatial_interp_get
(
  p_server                 svr
);

/**
 * \brief Recv_freq_set CWIPI. <b>(Not implemented yet)</b>
 *
 * \param [in]  p_server       Pointer on server data structure
 *
 */

void
CWP_server_Recv_freq_set
(
  p_server                 svr
);

/**
 * \brief next_recv_time_set CWIPI. <b>(Not implemented yet)</b>
 *
 * \param [in]  p_server       Pointer on server data structure
 *
 */

void
CWP_server_next_recv_time_set
(
  p_server                 svr
);

/**
 * \brief Cpl_time_step_set CWIPI. <b>(Not implemented yet)</b>
 *
 * \param [in]  p_server       Pointer on server data structure
 *
 */

void
CWP_server_Cpl_time_step_set
(
  p_server                 svr
);


/**
 * \brief Field_exch CWIPI. <b>(Not implemented yet)</b>
 *
 * \param [in]  p_server       Pointer on server data structure
 *
 */

void
CWP_server_Field_exch
(
  p_server                 svr
);


/**
 * \brief Mesh_interf_shared_pdm_nodal CWIPI. <b>(Not implemented yet)</b>
 *
 * \param [in]  p_server       Pointer on server data structure
 *
 */


void
CWP_server_Mesh_interf_shared_pdm_nodal
(
  p_server                 svr
);


/**
 *
 * \brief Field_data_type_get CWIPI.  <b>(Not implemented yet)</b>
 *
 * \param [in]  p_server       Pointer on server data structure
 *
 */

void
CWP_server_Field_data_type_get
(
  p_server                 svr
);


/**
 *
 * \brief Field_gradient_data_set CWIPI. <b>(Not implemented yet)</b>
 *
 * \param [in]  p_server       Pointer on server data structure
 *
 */

void
CWP_server_Field_gradient_data_set
(
  p_server                 svr
);

/*=============================================================================
 * Server function interfaces
 *============================================================================*/

/* Create a server */

int
CWP_server_create
(
 MPI_Comm global_comm,
 uint16_t server_port,
 int flags,
 p_server svr
);

/* Kill a server */

int
CWP_server_kill
(
 p_server svr
);

/* Message handler */

int
CWP_server_msg_handler
(
 p_server svr,
 p_message msg
);

/* Run a server */

int
CWP_server_run
(
 p_server svr
);

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* __SERVER_H__ */
