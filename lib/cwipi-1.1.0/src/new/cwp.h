/*
 * \file
 */

#ifndef __CWP_H__
#define __CWP_H__

/*
  This file is part of the CWIPI library.

  Copyright (C) 2021-2023  ONERA

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


/** \file cwp.h
  * \brief CWIPI new API header file
  *
  */


#include <mpi.h>
#include <stdio.h>


/*=============================================================================
 * Macro definitions
 *============================================================================*/


/**
 * \cond
 */

#if !defined (__hpux) && !defined (_AIX)
#define PROCF(x, y) x##_
#else
#define PROCF(x, y) x
#endif


/**
 * \endcond
 */

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force back Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */

/*============================================================================
 * Type
 *============================================================================*/

/**
 * \typedef CWP_g_num_t
 * \brief Long int in CWIPI
 *
 */

typedef long CWP_g_num_t;

/*============================================================================
 * Enumeration definitions
 *============================================================================*/

/**
 * \enum CWP_Type_t
 * \brief  Type of data.
 */

typedef enum {

  CWP_DOUBLE,      /*!< Double precision type */
  CWP_INT,         /*!< Integer type */
  CWP_CHAR         /*!< String type */

} CWP_Type_t;

/**
 * \enum CWP_Visu_format_t
 * \brief  List of available formats for Visualization.
 *
 */

typedef enum {

  CWP_VISU_FORMAT_ENSIGHT      /*!< Ensight visualization format */

} CWP_Visu_format_t;


/**
 * \enum CWP_Comm_t
 * \brief Communication mode
 *
 * CWP_Comm_t gives the different communication mode
 */

typedef enum {

  CWP_COMM_PAR_WITH_PART,    /*!< Parallel communication
                                  on partitioned source mesh */
  CWP_COMM_PAR_WITHOUT_PART, /*!< Parallel communication
                                  on unpartitioned source mesh defined on
                                  all processes */
//  CWP_COMM_INTERNAL         /*!< Internal communication within a process */

} CWP_Comm_t;

/**
 * \enum CWP_Time_exch_t
 * \brief  Modes of time exchange used to determine when to trigger interpolation/exchange.
 *
 */


typedef enum {

  CWP_TIME_EXCH_USER_CONTROLLED      /*!< Exchanges are used controlled */
  // CWP_TIME_EXCH_EACH_TIME_STEP,      /*!< Exchange at each time step */
  // CWP_TIME_EXCH_N_TIME_STEP,         /*!< Exchange every <EM> n </EM> time steps  */
  // CWP_TIME_EXCH_CPL_TIME_STEP,       /*!< Coupling time step        */
  // CWP_TIME_EXCH_ASYNCHRONOUS,        /*!< Exchanges are asynchronous with temporal interpolation */
  // CWP_TIME_EXCH_SLAVE,               /*!< Give a converged state    */
  // CWP_TIME_EXCH_MASTER               /*!< Request a converged state */
} CWP_Time_exch_t;


/**
 * \enum CWP_Dof_location_t
 * \brief Modes of degrees of freedom location.
 *
 */

typedef enum {
  CWP_DOF_LOCATION_UNDEF,   /*!< Location is undefined  */
  CWP_DOF_LOCATION_CELL_CENTER,   /*!< Field defined in cell point  */
  CWP_DOF_LOCATION_NODE,     /*!< Cell vertex field */
  CWP_DOF_LOCATION_USER           /*!< User defined field */
} CWP_Dof_location_t ;


/**
 * \enum CWP_Field_exch_t
 * \brief Modes of field exchange.
 */

typedef enum {

  CWP_FIELD_EXCH_SEND,        /*!< Send */
  CWP_FIELD_EXCH_RECV,        /*!< Receive */
  CWP_FIELD_EXCH_SENDRECV     /*!< Send and receive */

} CWP_Field_exch_t ;

/**
 * \enum CWP_PartData_exch_t
 * \brief Modes of field exchange.
 */

typedef enum {

  CWP_PARTDATA_SEND,        /*!< Send */
  CWP_PARTDATA_RECV,        /*!< Receive */

} CWP_PartData_exch_t ;


/**
 * \enum CWP_Field_map_t
 * \brief Modes of field exchange.
 */

typedef enum {

  CWP_FIELD_MAP_SOURCE,        /*!< Data array used to send field */
  CWP_FIELD_MAP_TARGET,        /*!< Data array used to receive interpolated field*/

} CWP_Field_map_t ;


/**
 * \enum CWP_Field_storage_t
 * \brief Modes of field storage.
 */

typedef enum {

  CWP_FIELD_STORAGE_INTERLACED,  /*!< Interlaced storage (x1, y1, z1, ... , xn, yn, zn) */
  CWP_FIELD_STORAGE_INTERLEAVED  /*!< Interleaved storage (x1, ... xn, y1, ..., yn, z1, ...zn) */

} CWP_Field_storage_t ;


/**
 * \enum CWP_Status_t
 * \brief on/off status
 *
 */

typedef enum {

  CWP_STATUS_OFF,           /*!< OFF */
  CWP_STATUS_ON             /*!< ON */

} CWP_Status_t;

/**
 * \enum CWP_Err_t
 * \brief Error codes.
 *
 */

typedef enum {

  CWP_ERR_NO_ERROR,       /*!< Output without error */
  CWP_ERR_DEFAULT         /*!< Output with default error */

} CWP_Err_t;

/**
 * \enum CWP_Block_t
 * \brief Mesh elements supported by CWIPI
 *
 * (alias to PDM_Mesh_nodal_t)
 *
 * CWP_Block_t defines all supported mesh elements
 */

typedef enum {

  CWP_BLOCK_NODE,          /*!< Node */
  CWP_BLOCK_EDGE2,         /*!< Edge with two nodes */
  CWP_BLOCK_FACE_TRIA3,    /*!< Triangle with three nodes */
  CWP_BLOCK_FACE_QUAD4,    /*!< Quadrangle with four nodes */
  CWP_BLOCK_FACE_POLY,     /*!< Generic polygon */
  CWP_BLOCK_CELL_TETRA4,   /*!< Tetrahedron with four nodes */
  CWP_BLOCK_CELL_PYRAM5,   /*!< Pyramid with five nodes */
  CWP_BLOCK_CELL_PRISM6,   /*!< Prism with six nodes */
  CWP_BLOCK_CELL_HEXA8,    /*!< Hexahedron with eight nodes */
  CWP_BLOCK_CELL_POLY,     /*!< Generic polyhedron */
  CWP_BLOCK_EDGEHO,        /*!< High-order Edge */
  CWP_BLOCK_FACE_TRIAHO,   /*!< High-order Triangle */
  CWP_BLOCK_FACE_QUADHO,   /*!< High-order Quadrangle */
  CWP_BLOCK_CELL_TETRAHO,  /*!< High-order Tetrahedron */
  CWP_BLOCK_CELL_PYRAMHO,  /*!< High-order Pyramid */
  CWP_BLOCK_CELL_PRISMHO,  /*!< High-order Prism */
  CWP_BLOCK_CELL_HEXAHO    /*!< High-order Hexahedron */

} CWP_Block_t;


/**
 * \enum CWP_Dynamic_mesh_t
 * \brief Modes of dynamic mesh.
 *
 */

typedef enum {

  CWP_DYNAMIC_MESH_STATIC,      /*!< Static mesh*/
  CWP_DYNAMIC_MESH_DEFORMABLE,  /*!< Deformable mesh with constant topology */
  CWP_DYNAMIC_MESH_VARIABLE     /*!< Variable mesh topology */

} CWP_Dynamic_mesh_t;

/**
 * \enum CWP_Spatial_interp_t
 * \brief List of available spatial interpolation methods.
  */

typedef enum {

  CWP_SPATIAL_INTERP_FROM_NEAREST_SOURCES_LEAST_SQUARES,         /*!< Least squares from nearest sources */
  CWP_SPATIAL_INTERP_FROM_NEAREST_TARGETS_LEAST_SQUARES,         /*!< Least squares from nearest targets */
  CWP_SPATIAL_INTERP_FROM_INTERSECTION,                          /*!< Mesh intersection */
  CWP_SPATIAL_INTERP_FROM_LOCATION_MESH_LOCATION_LOCATE_ALL_TGT, /*!< Location into a mesh (all targets are located) */
  CWP_SPATIAL_INTERP_FROM_LOCATION_MESH_LOCATION_OCTREE,         /*!< Location into a mesh with the octree method*/
  CWP_SPATIAL_INTERP_FROM_LOCATION_MESH_LOCATION_BOXTREE,        /*!< Location into a mesh with the bounding box tree method */
  CWP_SPATIAL_INTERP_FROM_IDENTITY                               /*!< Identity mapping */

} CWP_Spatial_interp_t;

/**
 * \enum CWP_Interface_t
 * \brief Coupling interfaces.
 */

typedef enum{

  CWP_INTERFACE_POINT,    /*!< Point interface */
  CWP_INTERFACE_LINEAR,   /*!< Linear interface */
  CWP_INTERFACE_SURFACE,  /*!< Surface interface */
  CWP_INTERFACE_VOLUME    /*!< Volume interface */

} CWP_Interface_t;

/**
 * \enum CWP_State_t
 * \brief Code state.
 *
 */

typedef enum {

  CWP_STATE_IN_PROGRESS,         /*!< Code running */
  CWP_STATE_END,                 /*!< Computation end */
  CWP_STATE_OUTPUT_ERROR         /*!< Output on error */

} CWP_State_t;

/**
 * \enum CWP_Op_t
 * \brief Reduction operations on control parameters.
 *
 * \warning Reduction operations are only available for control parameters of type ``CWP_INT`` or ``CWP_DOUBLE``.
 *
 */

typedef enum {

  CWP_OP_MIN,            /*!< Minimum */
  CWP_OP_MAX,            /*!< Maximum */
  CWP_OP_SUM             /*!< Sum */

} CWP_Op_t;

/*============================================================================
 * User interpolation type
 *============================================================================*/

/**
 * \typedef void (*CWP_Interp_function_t)
 * \brief User interpolation function interface.
 *
 * void (*CWP_Interp_function_t) defines the user interpolation
 * interface to take into account an user interpolation
 * Use \ref CWP_Field_interp_function_set to activate
 * the function.
 *
 * \param [in]  local_code_name             Local code name
 * \param [in]  cpl_id                      Coupling name
 * \param [in]  field_id                    Field name
 * \param [in]  i_part                      Partition identifier
 * \param [in]  buffer_in                   Input field array
 * \param [out] buffer_out                  Output field array
 *
 */

typedef void (*CWP_Interp_function_t)
(
 const char           *local_code_name,
 const char           *cpl_id,
 const char           *field_id,
 int                   i_part,
 double               *buffer_in,
 double               *buffer_out
);


/*=============================================================================
 * Static global variables
 *============================================================================*/

/*============================================================================*
 *                                                                            *
 *                      Public function prototypes                            *
 *                      --------------------------                            *
 *                                                                            *
 *============================================================================*/

/*----------------------------------------------------------------------------*
 * General functions                                                          *
 *----------------------------------------------------------------------------*/


/**
 * \brief Initialize CWIPI.
 *
 * This function creates the MPI intra communicators of the codes from
 * the \p global_comm MPI communicator that contains all code ranks. This
 * function has to be called from all ranks contained in the \p global_comm.
 *
 * \param [in]  global_comm    MPI global communicator
 * \param [in]  n_code         Number of codes on the current rank
 * \param [in]  code_names     Names of codes on the current rank (size = \p n_code)
 * \param [in]  is_active_rank Is current rank have to be used by CWIPI
 * \param [out] intra_comms    MPI intra communicators of each code (size = \p n_code)
 *
 */

void
CWP_Init
(
 const MPI_Comm           global_comm,
 const int                n_code,
 const char             **code_names,
 const CWP_Status_t       is_active_rank,
 MPI_Comm                *intra_comms
);

/**
 *
 * \brief Finalize CWIPI.
 *
 */

void
CWP_Finalize
(
 void
);

/*----------------------------------------------------------------------------*
 * Functions about current code properties                                    *
 *----------------------------------------------------------------------------*/

/**
 * \brief Update code state.
 *
 * \param [in] local_code_name  Local code name
 * \param [in] state            State
 *
 */

void
CWP_State_update
(
 const char*       local_code_name,
 const CWP_State_t state
);

/**
 * \brief Begin code time step.
 *
 * \param [in] local_code_name  Local code name
 * \param [in] current_time     Current time
 *
 */

void
CWP_Time_step_beg
(
 const char*  local_code_name,
 const double current_time
);

/**
 * \brief End code time step.
 *
 * \param [in] local_code_name  Local code name
 *
 */

void
CWP_Time_step_end
(
 const char* local_code_name
);


/**
 * \brief Define output file.
 *
 * \param [in] output_file    Output file
 *
 */

void
CWP_Output_file_set
(
FILE *output_file
);


/**
 * \brief Define a user structure associated to a code
 * 
 * This structure can be accessed into a callback
 *
 * \param [in] local_code_name  Local code name
 * \param [in] user_structure   User structure
 *
 */

void
CWP_User_structure_set
(
 const char* local_code_name,
       void* user_structure
);


/**
 * \brief Return the user structure associated 
 * 
 * This structure can be called into a callback 
 *
 * \param [in] local_code_name  Local code name
 * 
 * \return  User structure
 *
 */

void *
CWP_User_structure_get
(
 const char* local_code_name
);

//**
// * \brief Writing output to Fortran file.
// *
// * This function set the file Fortran logical unit for writing output.
// *
// * \param [in]  iunit        File Fortran logical unit
// *
// */

//void
//PROCF (cwp_output_fortran_unit_set, CWP_OUTPUT_FORTRAN_UNIT_SET)
//(
// int *iunit
//);

/*----------------------------------------------------------------------------*
 * Functions about other code properties                               *
 *----------------------------------------------------------------------------*/

/**
 * \brief Return code state.
 *
 * \param [in]  code_name    Code name
 *
 * \return      Code state
 */

CWP_State_t
CWP_State_get
(
 const char    *code_name
);


/**
 * \brief Return the number of codes known by CWIPI.
 *
 * \return Number of codes
 *
 */

int
CWP_Codes_nb_get
(
 void
);

/**
 * \brief Return the list of code names known by CWIPI.
 *
 * \return List of codes.
 */

const char **
CWP_Codes_list_get
(
void
);


/**
 * \brief Return the number of local codes known by CWIPI.
 *
 * \return Number of local codes.
 */

int
CWP_Loc_codes_nb_get
(
 void
);


/**
 * \brief Return the list of local code names known by CWIPI.
 *
 * \return List of local codes.
 */

const char **
CWP_Loc_codes_list_get
(
 void
);

/*----------------------------------------------------------------------------*
 * Functions about properties                                                 *
 *----------------------------------------------------------------------------*/

/**
 * \brief Dump code properties.
 *
 */

void
CWP_Properties_dump
(
void
);

/**
 * \brief Dump string of code properties.
 *
 */

int
CWP_Properties_str_dump
(
 char **char_out
);

/*----------------------------------------------------------------------------*
 * General functions about coupling                                           *
 *----------------------------------------------------------------------------*/

/**
 * \brief Create a coupling object and define its properties.
 *
 * \param [in]  local_code_name     Local code name
 * \param [in]  cpl_id              Coupling identifier
 * \param [in]  coupled_code_name   Distant or local coupled code name
 * \param [in]  entities_dim        Coupling interface type
 * \param [in]  comm_type           Communication type
 * \param [in]  spatial_interp      Spatial interpolation method
 * \param [in]  n_part              Number of interface partition
 * \param [in]  displacement        Mesh moving status
 * \param [in]  recv_freq_type      Type of receiving frequency
 *
 */

void
CWP_Cpl_create
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
 * \brief MPI Barrier on the coupling communicator.
 *
 * \param [in] local_code_name  Local code name
 * \param [in] cpl_id           Coupling identifier
 *
 */

void
CWP_Cpl_barrier
(
 const char *local_code_name,
 const char *cpl_id
);


/**
 *
 * \brief Delete a coupling object.
 *
 * \param [in] local_code_name  Local code name
 * \param [in] cpl_id           Coupling identifier
 *
 */

void
CWP_Cpl_del
(
 const char *local_code_name,
 const char *cpl_id
);

/**
 *
 * \brief Enable broadcast of the computed targets ids (in \ref CWP_COMM_PAR_WITHOUT_PART mode)
 *
 * This function must be called in order for the computed targets to be accessible
 * on non-root ranks
 *
 * \param [in] local_code_name  Local code name
 * \param [in] cpl_id           Coupling identifier
 * \param [in] field_id         Field identifier
 *
 */

void
CWP_Computed_tgts_bcast_enable
(
  const char *local_code_name,
  const char *cpl_id,
  const char *field_id
);

/**
 *
 * \brief Return the number of uncomputed targets.
 *
 * \param [in] local_code_name  Local code name
 * \param [in] cpl_id           Coupling identifier
 * \param [in] field_id         Field identifier
 * \param [in] i_part           Current partition
 *
 * \return Number of uncomputed targets
 */

int
CWP_N_uncomputed_tgts_get
(
 const char *local_code_name,
 const char *cpl_id,
 const char *field_id,
 const int   i_part
);


/**
 *
 * \brief Return uncomputed targets.
 *
 * \param [in] local_code_name  Local code name
 * \param [in] cpl_id           Coupling identifier
 * \param [in] field_id         Field identifier
 * \param [in] i_part           Current partition
 *
 * \return Uncomputed targets
 */

const int *
CWP_Uncomputed_tgts_get
(
 const char *local_code_name,
 const char *cpl_id,
 const char *field_id,
 const int   i_part
);

/**
 *
 * \brief Return the number of computed targets.
 *
 * \param [in] local_code_name  Local code name
 * \param [in] cpl_id           Coupling identifier
 * \param [in] field_id         Field identifier
 * \param [in] i_part           Current partition
 *
 * \return Number of computed targets
 */

int
CWP_N_computed_tgts_get
(
 const char *local_code_name,
 const char *cpl_id,
 const char *field_id,
 const int   i_part
);

/**
 *
 * \brief Return computed targets.
 *
 * \param [in] local_code_name  Local code name
 * \param [in] cpl_id           Coupling identifier
 * \param [in] field_id         Field identifier
 * \param [in] i_part           Current partition
 *
 * \return Computed targets
 */

const int *
CWP_Computed_tgts_get
(
 const char *local_code_name,
 const char *cpl_id,
 const char *field_id,
 const int   i_part
);


/**
 *
 * \brief Enable broadcast of the involved sources ids (in \ref CWP_COMM_PAR_WITHOUT_PART mode)
 *
 * This function must be called in order for the involved sources to be accessible
 * on non-root ranks
 *
 * \param [in] local_code_name  Local code name
 * \param [in] cpl_id           Coupling identifier
 * \param [in] field_id         Field identifier
 *
 */

void
CWP_Involved_srcs_bcast_enable
(
  const char *local_code_name,
  const char *cpl_id,
  const char *field_id
);


/**
 *
 * \brief Return number of involved sources.
 *
 * \param [in] local_code_name  Local code name
 * \param [in] cpl_id           Coupling identifier
 * \param [in] field_id         Field identifier
 * \param [in] i_part           Current partition
 *
 * \return Number of involved sources
 */

int
CWP_N_involved_srcs_get
(
 const char *local_code_name,
 const char *cpl_id,
 const char *field_id,
 const int   i_part
);


/**
 *
 * \brief Return involved sources
 *
 * \param [in] local_code_name  Local code name
 * \param [in] cpl_id           Coupling identifier
 * \param [in] field_id         Field identifier
 * \param [in] i_part           Current partition
 *
 * \return Involved sources
 */

const int *
CWP_Involved_srcs_get
(
 const char *local_code_name,
 const char *cpl_id,
 const char *field_id,
 const int   i_part
);

/*----------------------------------------------------------------------------*
 * Functions about spatial interpolation                                      *
 *----------------------------------------------------------------------------*/

/**
 * \brief Compute spatial interpolation weights.
 *
 * \param [in]  local_code_name     Local code name
 * \param [in]  cpl_id              Coupling identifier
 *
 */

void
CWP_Spatial_interp_weights_compute
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
CWP_Spatial_interp_property_set
(
 const char       *local_code_name,
 const char       *cpl_id,
 const char       *property_name,
 const CWP_Type_t  property_type,
 const char       *property_value
);

/*----------------------------------------------------------------------------*
 * Functions about visualization                                              *
 *----------------------------------------------------------------------------*/

/**
 * \brief Enable visualization output.
 *
 * \param [in]  local_code_name  Local code name
 * \param [in]  cpl_id           Coupling identifier
 * \param [in]  freq             Output frequency
 * \param [in]  format           Output format to visualize exchanged fields
 * \param [in]  format_option    Output options "opt1, opt2, ..."
 *                                - text : output text files
 *                                - binary : output binary files (default)
 */

void
CWP_Visu_set
(
 const char                 *local_code_name,
 const char                 *cpl_id,
 const int                   freq,
 const CWP_Visu_format_t     format,
 const char                 *format_option
);

/*----------------------------------------------------------------------------*
 * Functions about User target points                                         *
 *----------------------------------------------------------------------------*/

/**
 * \brief Set a partition of the user target point cloud.
 *
 * This function must be called if the degrees of freedom locations are
 * \ref CWP_DOF_LOCATION_USER
 *
 * \param [in]  local_code_name  Local code name
 * \param [in]  cpl_id           Coupling identifier
 * \param [in]  i_part           Current partition
 * \param [in]  n_pts            Number of points
 * \param [in]  coord            Coordinates (size = 3 * \p n_pts)
 * \param [in]  global_num       Global ids (size = \p n_pts or NULL)
 *
 */

void
CWP_User_tgt_pts_set
(
 const char    *local_code_name,
 const char    *cpl_id,
 const int      i_part,
 const int      n_pts,
 double         coord[],
 CWP_g_num_t    global_num[]
);

/*----------------------------------------------------------------------------*
 * Functions about Mesh                                                    *
 *----------------------------------------------------------------------------*/


/**
 * \brief Finalize the interface mesh.
 *
 * This function computes the global ids of mesh entities if they are
 * not provided.
 *
 * \param [in]  local_code_name  Local code name
 * \param [in]  cpl_id           Coupling identifier
 *
 */

void
CWP_Mesh_interf_finalize
(
 const char           *local_code_name,
 const char           *cpl_id
);


/**
 * \brief Set the interface mesh vertices.
 *
 * \param [in]  local_code_name  Local code name
 * \param [in]  cpl_id           Coupling identifier
 * \param [in]  i_part           Current partition
 * \param [in]  n_vtx            Number of vertices
 * \param [in]  coord            Coordinates (size = 3 * \p n_vtx)
 * \param [in]  global_num       Global vertex ids (size = \p n_vtx or NULL)
 *
 */

void
CWP_Mesh_interf_vtx_set
(
 const char           *local_code_name,
 const char           *cpl_id,
 const int             i_part,
 const int             n_vtx,
 double                coord[],
 CWP_g_num_t           global_num[]
);


/**
 * \brief Add a connectivity block to the interface mesh.
 *
 * \param [in]  local_code_name  Local code name
 * \param [in]  cpl_id           Coupling identifier
 * \param [in]  block_type       Block type
 *
 * \return block identifier
 */

int
CWP_Mesh_interf_block_add
(
 const char           *local_code_name,
 const char           *cpl_id,
 const CWP_Block_t     block_type
);


/**
 * \brief Set a standard block to the interface mesh.
 *
 * This function adds a connectivity block to the interface mesh.
 * Definition of element connectivity is :
 *
 *  - edge (\ref CWP_BLOCK_EDGE2) :
 *
 *   \code
 *       1 x-------x 2
 *   \endcode
 *
 *  - triangle (\ref CWP_BLOCK_FACE_TRIA3):
 *
 *   \code
 *       1 x-------x 3
 *          \     /
 *           \   /
 *            \ /
 *             x 2
 *   \endcode
 *
 *  - quadrangle (\ref CWP_BLOCK_FACE_QUAD4) :
 *
 *   \code
 *          4 x-------x 3
 *           /       /
 *          /       /
 *       1 x-------x2
 *   \endcode
 *
 *  - tetrahedron (\ref CWP_BLOCK_CELL_TETRA4) :
 *
 *   \code
 *             x 4
 *            /|\
 *           / | \
 *          /  |  \
 *       1 x- -|- -x 3
 *          \  |  /
 *           \ | /
 *            \|/
 *             x 2
 *   \endcode
 *
 *  - pyramid (\ref CWP_BLOCK_CELL_PYRAM5) :
 *
 *   \code
 *              5 x
 *               /|\
 *              //| \
 *             // |  \
 *          4 x/--|---x 3
 *           //   |  /
 *          //    | /
 *       1 x-------x 2
 *   \endcode
 *
 *  - prism (\ref CWP_BLOCK_CELL_PRISM6) :
 *
 *   \code
 *       4 x-------x 6
 *         |\     /|
 *         | \   / |
 *       1 x- \-/ -x 3
 *          \ 5x  /
 *           \ | /
 *            \|/
 *             x 2
 *   \endcode
 *
 *  - hexahedron (\ref CWP_BLOCK_CELL_HEXA8) :
 *
 *   \code
 *          8 x-------x 7
 *           /|      /|
 *          / |     / |
 *       5 x-------x6 |
 *         | 4x----|--x 3
 *         | /     | /
 *         |/      |/
 *       1 x-------x 2
 *   \endcode
 *
 * \param [in]  local_code_name  Local code name
 * \param [in]  cpl_id           Coupling identifier
 * \param [in]  i_part           Partition identifier
 * \param [in]  block_id         Block identifier
 * \param [in]  n_elts           Number of elements
 * \param [in]  connec           Connectivity (size = *n_vertex_per_elt* * \p n_elts)
 * \param [in]  global_num       Global element ids (size = \p n_elts or NULL)
 *
 */

void
CWP_Mesh_interf_block_std_set
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
 * \brief Get the properties of a standard block of the interface mesh.
 *
 * \param [in]  local_code_name  Local code name
 * \param [in]  cpl_id           Coupling identifier
 * \param [in]  i_part           Partition identifier
 * \param [in]  block_id         Block identifier
 * \param [out] n_elts           Number of elements
 * \param [out] connec           Connectivity (size = *n_vertex_per_elt* * \p n_elts)
 * \param [out] global_num       Global element ids (size = \p n_elts or NULL)
 */

void
CWP_Mesh_interf_block_std_get
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
  * \param [in]  block_id         Block identifier
  *
  * \return block type
  */

CWP_Block_t
CWP_std_block_type_get
(
 const char             *local_code_name,
 const char             *cpl_id,
 const int               block_id
);


/**
 * \brief Set the connectivity of a polygon block in an interface mesh partition.
 *
 * \param [in]  local_code_name  Local code name
 * \param [in]  cpl_id           Coupling identifier
 * \param [in]  i_part           Current partition
 * \param [in]  block_id         Block identifier
 * \param [in]  n_elts           Number of elements
 * \param [in]  connec_idx       Connectivity index (\p connec_idx[0] = 0 and
 *                               size = \p n_elts + 1)
 * \param [in]  connec           Connectivity (size = \p connec_idx[\p n_elts])
 * \param [in]  global_num       Global element ids (size = \p n_elts or NULL)
 *
 */

void
CWP_Mesh_interf_f_poly_block_set
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
 * \brief Get the properties of a polygon block of the interface mesh partition.
 *
 * \param [in]  local_code_name  Local code name
 * \param [in]  cpl_id           Coupling identifier
 * \param [in]  i_part           Current partition
 * \param [in]  block_id         Block identifier
 * \param [out] n_elts           Number of elements
 * \param [out] connec_idx       Connectivity index (\p connec_idx[0] = 0 and
 *                               size = \p n_elts + 1)
 * \param [out] connec           Connectivity (size = \p connec_idx[\p n_elts])
 * \param [out] global_num       Global element ids (size = \p n_elts or NULL)
 *
 */

void
CWP_Mesh_interf_f_poly_block_get
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
 * \brief Set the properties of a polyhedron block of the interface mesh partition.
 *
 * \param [in]  local_code_name   Local code name
 * \param [in]  cpl_id            Coupling identifier
 * \param [in]  i_part            Current partition
 * \param [in]  block_id          Block identifier
 * \param [in]  n_elts            Number of elements
 * \param [in]  n_faces           Number of faces
 * \param [in]  connec_faces_idx  Polyhedron face to vertex index
 *                                (\p connec_faces_idx[0] = 0 and
 *                                 size = max(\p connec_faces) + 1)
 * \param [in]  connec_faces      Polyhedron face to vertex connectivity
 *                                (size = \p connec_faces_idx[\p n_elts])
 * \param [in]  connec_cells_idx  Polyhedron to face index
 *                                (\p connec_cells_idx[0] = 0 and
 *                                 size = \p n_elts + 1)
 * \param [in]  connec_cells      Polyhedron to face connectivity
 *                                (size = \p connec_cells_idx[\p n_elts])
 * \param [in]  global_num        Global element ids (size = \p n_elts or NULL)
 *
 */

void
CWP_Mesh_interf_c_poly_block_set
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
 * \brief Get the properties of a polyhedron block of the interface mesh partition.
 *
 * \param [in]  local_code_name   Local code name
 * \param [in]  cpl_id            Coupling identifier
 * \param [in]  i_part            Current partition
 * \param [in]  block_id          Block identifier
 * \param [out] n_elts            Number of elements
 * \param [out] n_faces           Number of faces
 * \param [out] connec_faces_idx  Polyhedron face to vertex index
 *                                (\p connec_faces_idx[0] = 0 and
 *                                 size = max(\p connec_cells) + 1)
 * \param [out] connec_faces      Polyhedron face to vertex connectivity
 *                                (size = \p connec_faces_idx[\p n_elts])
 * \param [out] connec_cells_idx  Polyhedron to face index
 *                                (\p connec_cells_idx[0] = 0 and
 *                                 size = \p n_elts + 1)
 * \param [out] connec_cells      Polyhedron to face connectivity
 *                                (size = \p connec_cells_idx[\p n_elts])
 * \param [out] global_num        Global element ids (size = \p n_elts or NULL)
 *
 */

void
CWP_Mesh_interf_c_poly_block_get
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
 * \brief Delete the interface mesh.
 *
 * \param [in] local_code_name  Local code name
 * \param [in] cpl_id           Coupling identifier
 *
 */

void
CWP_Mesh_interf_del
(
 const char *local_code_name,
 const char *cpl_id
);


/**
 * \brief Define the interface mesh from a cell-to-face connectivity.
 *
 * \param [in]  local_code_name   Local code name
 * \param [in]  cpl_id            Coupling identifier
 * \param [in]  i_part            Current partition
 * \param [in]  n_cells           Number of cells
 * \param [in]  cell_face_idx     Polyhedron to face index
 *                                (\p cell_face_idx[0] = 0 and
 *                                 size = \p n_elts + 1)
 * \param [in]  cell_face         Polyhedron to face connectivity
 *                                (size = \p cell_face_idx[\p n_elts])
 * \param [in]  n_faces           Number of faces
 * \param [in]  face_vtx_idx      Face to vertex index
 *                                (\p face_vtx_idx[0] = 0 and
 *                                 size = \p n_faces + 1)
 * \param [in]  face_vtx          Face to vertex connectivity
 *                                (size = \p face_vtx_idx[\p n_elts])
 * \param [in]  global_num        Global polyhedron ids (size = \p n_elts or NULL)
 *
 */

void
CWP_Mesh_interf_from_cellface_set
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
 * \brief Define the surface interface mesh from a face-to-edge connectivity.
 *
 * \param [in]  local_code_name   Local code name
 * \param [in]  cpl_id            Coupling identifier
 * \param [in]  i_part            Current partition
 * \param [in]  n_faces           Number of cells
 * \param [in]  face_edge_idx     Polygon to edge index
 *                                (\p face_edge_idx[0] = 0 and
 *                                 size =  \p n_faces + 1)
 * \param [in]  face_edge         Polygon to edge connectivity
 *                                (size = \p face_edge_idx[\p n_faces])
 * \param [in]  n_edges           Number of faces
 * \param [in]  edge_vtx          Edge to vertex connectivity
 *                                (size = 2 * \p n_edges)
 * \param [in]  global_num        Global polygon ids (size = \p n_faces or NULL)
 *
 */

void
CWP_Mesh_interf_from_faceedge_set
(
 const char        *local_code_name,
 const char        *cpl_id,
 const int          i_part,
 const int          n_faces,
       int          face_edge_idx[],
       int          face_edge[],
 const int          n_edges,
       int          edge_vtx[],
       CWP_g_num_t  global_num[]
);


/**
 * \brief Define the surface interface mesh from a face-to-vertex connectivity.
 *
 * \param [in]  local_code_name   Local code name
 * \param [in]  cpl_id            Coupling identifier
 * \param [in]  i_part            Current partition
 * \param [in]  n_faces           Number of cells
 * \param [in]  face_vtx_idx      Polygon to vertex index
 *                                (\p face_vtx_idx[0] = 0 and
 *                                 size =  \p n_faces + 1)
 * \param [in]  face_vtx          Polygon to vertex connectivity
 *                                (size = \p face_vtx_idx[\p n_faces])
 * \param [in]  global_num        Global polygon ids (size = \p n_faces or NULL)
 *
 */

void
CWP_Mesh_interf_from_facevtx_set
(
 const char        *local_code_name,
 const char        *cpl_id,
 const int          i_part,
 const int          n_faces,
       int          face_vtx_idx[],
       int          face_vtx[],
       CWP_g_num_t  global_num[]
);

/*----------------------------------------------------------------------------*
 * Functions about field                                                      *
 *----------------------------------------------------------------------------*/

/**
 *
 * \brief Create a new field.
 *
 * \param [in]  local_code_name Local code name
 * \param [in]  cpl_id          Coupling identifier
 * \param [in]  field_id        Field id
 * \param [in]  data_type       Data type
 * \param [in]  storage         Storage type
 * \param [in]  n_component     Number of components
 * \param [in]  dof_location    Location of the degrees of freedom
 * \param [in]  exch_type       Exchange type
 * \param [in]  visu_status     Visualization status
 *
 */

void
CWP_Field_create
(
 const char                  *local_code_name,
 const char                  *cpl_id,
 const char                  *field_id,
 const CWP_Type_t             data_type,
 const CWP_Field_storage_t    storage,
 const int                    n_component,
 const CWP_Dof_location_t     dof_location,
 const CWP_Field_exch_t       exch_type,
 const CWP_Status_t           visu_status
);


/**
 *
 * \brief Set field data.
 *
 * \param [in] local_code_name   Local code name
 * \param [in] cpl_id            Coupling identifier
 * \param [in] field_id          Field identifier
 * \param [in] i_part            Current partition
 * \param [in] map_type          Choice if data is set for the source or the target
 * \param [in] data              Storage array (Mapping)
 *
 */

void
CWP_Field_data_set
(
 const char              *local_code_name,
 const char              *cpl_id,
 const char              *field_id,
 const int                i_part,
 const CWP_Field_map_t    map_type,
 double                   data[]
);

/**
 *
 * \brief Get field data.
 *
 * \param [in]  local_code_name   Local code name
 * \param [in]  cpl_id            Coupling identifier
 * \param [in]  field_id          Field identifier
 * \param [in]  i_part            Current partition
 * \param [in]  map_type          Choice if data is set for the source or the target
 * \param [out] data              Storage array (Mapping)
 *
 */

void
CWP_Field_data_get
(
 const char              *local_code_name,
 const char              *cpl_id,
 const char              *field_id,
 const int                i_part,
 const CWP_Field_map_t    map_type,
 double                 **data
);

/**
 *
 * \brief Get the degrees of freedom location.
 *
 * \param [in] local_code_name  Local code name
 * \param [in] cpl_id           Coupling identifier
 * \param [in] field_id         Field identifier
 *
 * \return                      Location of degrees of freedom
 *
 */

CWP_Dof_location_t
CWP_Field_dof_location_get
(
 const char      *local_code_name,
 const char      *cpl_id,
 const char      *field_id
);

/**
 *
 * \brief Get field storage type.
 *
 * \param [in] local_code_name  Local code name
 * \param [in] cpl_id           Coupling identifier
 * \param [in] field_id         Field identifier
 *
 * \return                      Field storage type
 */

CWP_Field_storage_t
CWP_Field_storage_get
(
 const char      *local_code_name,
 const char      *cpl_id         ,
 const char      *field_id
);

/**
 *
 * \brief Get number of field components.
 *
 * \param [in]  local_code_name           Local code name
 * \param [in]  cpl_id                    Coupling identifier
 * \param [in]  field_id                  Field identifier
 *
 * \return Number of field components
 *
 */

int
CWP_Field_n_components_get
(
 const char             *local_code_name,
 const char             *cpl_id,
 const char             *field_id
);

/**
 *
 * \brief Get number of field degrees of freedom.
 *
 * \param [in]  local_code_name           Local code name
 * \param [in]  cpl_id                    Coupling identifier
 * \param [in]  field_id                  Field identifier
 * \param [in]  i_part                    Partition identifier
 *
 */

int
CWP_Field_n_dof_get
(
 const char             *local_code_name,
 const char             *cpl_id,
 const char             *field_id,
 int                     i_part
);

/**
 *
 * \brief Get spatial interpolation source data.
 *
 * \param [in]  local_code_name           Local code name
 * \param [in]  cpl_id                    Coupling identifier
 * \param [in]  field_id                  Field identifier
 * \param [in]  i_part                    Partition identifier
 * \param [out] n_src                     Number of source dofs
 * \param [out] src_to_tgt_idx            Index for source->target mapping (size = \p n_src + 1)
 *
 */

void
CWP_Field_src_data_properties_get
(
 const char   *local_code_name,
 const char   *cpl_id,
 const char   *field_id,
 int           i_part,
 int          *n_src,
 int         **src_to_tgt_idx
);

/**
 *
 * \brief Get spatial interpolation target data.
 *
 * \param [in]  local_code_name Local code name
 * \param [in]  cpl_id          Coupling identifier
 * \param [in]  field_id        Field identifier
 * \param [in]  i_part          Partition identifier
 * \param [out] n_tgt           Number of target dofs
 * \param [out] n_computed_tgt  Number of computed target dofs
 * \param [out] computed_tgt    Computed target dofs (size = \p n_computed_tgt)
 * \param [out] tgt_to_src_idx  Index for target->source mapping (size = \p n_computed_tgt + 1)
 *
 */

void
CWP_Field_tgt_data_properties_get
(
 const char   *local_code_name,
 const char   *cpl_id,
 const char   *field_id,
 int           i_part,
 int          *n_tgt,
 int          *n_computed_tgt,
 int         **computed_tgt,
 int         **tgt_to_src_idx
);

/**
 *
 * \brief Get spatial interpolation weights (location algorithm).
 *
 * \param [in]  local_code_name           Local code name
 * \param [in]  cpl_id                    Coupling identifier
 * \param [in]  field_id                  Field identifier
 * \param [in]  i_part                    Partition identifier
 * \param [out] weights                   Spatial interpolation weights (barycentric coordinates)
 *
 */

void
CWP_Field_location_weights_get
(
 const char            *local_code_name,
 const char            *cpl_id,
 const char            *field_id,
 int                    i_part,
 double               **weights
);

/**
 *
 * \brief Get spatial interpolation point data (location algorithm).
 *
 * \param [in]  local_code_name           Local code name
 * \param [in]  cpl_id                    Coupling identifier
 * \param [in]  field_id                  Field identifier
 * \param [in]  i_part                    Partition identifier
 * \param [out] points_coords             Cartesian coordinates of points inside local elements
 * \param [out] points_uvw                Parametric coordinates of points inside local elements
 * \param [out] points_dist2              Squared distance from points to elements
 * \param [out] points_projected_coords   Cartesian coordinates of projection on points on local elements
 *
 */

void
CWP_Field_location_point_data_get
(
 const char            *local_code_name,
 const char            *cpl_id,
 const char            *field_id,
 int                    i_part,
 double               **points_coords,
 double               **points_uvw,
 double               **points_dist2,
 double               **points_projected_coords
);

/**
 *
 * \brief Get spatial interpolation internal cell->vertex connectivity (location algorithm).
 *
 * \param [in]  local_code_name           Local code name
 * \param [in]  cpl_id                    Coupling identifier
 * \param [in]  field_id                  Field identifier
 * \param [in]  i_part                    Partition identifier
 * \param [out] cell_vtx_idx              Index for local cell->vertex connectivity
 * \param [out] cell_vtx                  Local cell->vertex connectivity
 *
 */

void
CWP_Field_location_internal_cell_vtx_get
(
 const char            *local_code_name,
 const char            *cpl_id,
 const char            *field_id,
 int                    i_part,
 int                  **cell_vtx_idx,
 int                  **cell_vtx
);

/**
 *
 * \brief Get spatial interpolation volumes (intersection algorithm).
 *
 * \param [in]  local_code_name           Local code name
 * \param [in]  cpl_id                    Coupling identifier
 * \param [in]  field_id                  Field identifier
 * \param [in]  i_part                    Partition identifier
 * \param [out] volumes                   Volumes of intersection polyhedra
 *
 */

void
CWP_Field_intersection_volumes_get
(
 const char            *local_code_name,
 const char            *cpl_id,
 const char            *field_id,
 int                    i_part,
 double               **volumes
);

/**
 *
 * \brief Get spatial local target elements volumes (intersection algorithm).
 *
 * \param [in]  local_code_name           Local code name
 * \param [in]  cpl_id                    Coupling identifier
 * \param [in]  field_id                  Field identifier
 * \param [in]  i_part                    Partition identifier
 * \param [out] tgt_elt_volumes           Volumes of local target elements
 *
 */

void
CWP_Field_intersection_tgt_elt_volumes_get
(
 const char            *local_code_name,
 const char            *cpl_id,
 const char            *field_id,
 int                    i_part,
 double               **tgt_elt_volumes
);

/**
 *
 * \brief Get spatial interpolation distances (nearest neighbors algorithm).
 *
 * \param [in]  local_code_name           Local code name
 * \param [in]  cpl_id                    Coupling identifier
 * \param [in]  field_id                  Field identifier
 * \param [in]  i_part                    Partition identifier
 * \param [out] distances2                Squared distances from nearest source points
 *
 */

void
CWP_Field_nearest_neighbors_distances_get
(
 const char            *local_code_name,
 const char            *cpl_id,
 const char            *field_id,
 int                    i_part,
 double               **distances2
);

/**
 *
 * \brief Get coordinates of nearest source points (nearest neighbors algorithm).
 *
 * \param [in]  local_code_name           Local code name
 * \param [in]  cpl_id                    Coupling identifier
 * \param [in]  field_id                  Field identifier
 * \param [in]  i_part                    Partition identifier
 * \param [out] nearest_src_coord         Coordinates of nearest source points
 *
 */

void
CWP_Field_nearest_neighbors_coord_get
(
 const char            *local_code_name,
 const char            *cpl_id,
 const char            *field_id,
 int                    i_part,
 double               **nearest_src_coord
);

/**
 * \brief Delete a field.
 *
 * \param [in]  local_code_name Local code name
 * \param [in]  cpl_id          Coupling identifier
 * \param [in]  field_id        Field identifier
 *
 */

void
CWP_Field_del
(
 const char      *local_code_name,
 const char      *cpl_id         ,
 const char      *field_id
);

/*----------------------------------------------------------------------------*
 * Functions about field exchange                                             *
 *----------------------------------------------------------------------------*/


/**
 * \brief Send a spatially interpolated field to the coupled code with
 *        non-blocking communications.
 *
 * This function is independent of \ref CWP_Time_exch_t mode. The user has to
 * manually check the consistency of the exchanges.
 *
 * \param [in]  local_code_name Local code name
 * \param [in]  cpl_id          Coupling identifier
 * \param [in]  field_id        Field identifier
 *
 *
 */

void
CWP_Field_issend
(
 const char     *local_code_name,
 const char     *cpl_id,
 const char     *field_id
);

/**
 *
 * \brief Receive a spatially interpolated field from the coupled code
 *        with non-blocking communications.
 *
 * This function is independent of \ref CWP_Time_exch_t mode. The user has to
 * manually check the consistency of the exchanges.
 *
 * \param [in]  local_code_name Local code name
 * \param [in]  cpl_id          Coupling identifier
 * \param [in]  tgt_field_id    Target field id
 *
 *
 */

void
CWP_Field_irecv
(
 const char        *local_code_name,
 const char        *cpl_id,
 const char        *tgt_field_id
);

/**
 *
 * \brief Wait the end of an exchange initiated by \ref CWP_Field_issend.
 *
 * \param [in] local_code_name  Local code name
 * \param [in] cpl_id           Coupling identifier
 * \param [in] field_id         Field identifier
 *
 */

void
CWP_Field_wait_issend
(
 const char  *local_code_name,
 const char  *cpl_id,
 const char  *field_id
);


/**
 *
 * \brief Wait the end of an exchange initiated by from \ref CWP_Field_irecv.
 *
 * \param [in] local_code_name  Local code name
 * \param [in] cpl_id           Coupling identifier
 * \param [in] tgt_field_id     Target field id
 *
 */

void
CWP_Field_wait_irecv
(
 const char  *local_code_name,
 const char  *cpl_id,
 const char  *tgt_field_id
);


/*----------------------------------------------------------------------------*
 * Functions about user interpolation                                         *
 *----------------------------------------------------------------------------*/

/**
 *
 * \brief Unset a user-defined spatial interpolation function.
 *
 * \param [in] local_code_name  Local code name
 * \param [in] cpl_id           Coupling identifier
 * \param [in] field_id         Field identifier
 *
 */

void
CWP_Field_interp_function_unset
(
 const char                 *local_code_name,
 const char                 *cpl_id,
 const char                 *field_id
);


/**
 *
 * \brief Set a user-defined spatial interpolation function.
 *
 * \param [in] local_code_name  Local code name
 * \param [in] cpl_id           Coupling identifier
 * \param [in] field_id         Field identifier
 * \param [in] fct              Function
 *
 */

void
CWP_Field_interp_function_set
(
 const char                 *local_code_name,
 const char                 *cpl_id,
 const char                 *field_id,
 CWP_Interp_function_t       fct
);

/*----------------------------------------------------------------------------*
 * Functions about control parameters                                         *
 *----------------------------------------------------------------------------*/

/**
 *
 * \brief Add a new control parameter and initialize it.
 *
 * \param [in] local_code_name  Local code name
 * \param [in] param_name       Parameter name
 * \param [in] data_type        Parameter type
 * \param [in] initial_value    Initial value
 *
 */

void
CWP_Param_add
(
 const char        *local_code_name,
 const char        *param_name,
 const CWP_Type_t   data_type,
 void              *initial_value
);


/**
 *
 * \brief Set a control parameter value.
 *
 * \param [in] local_code_name  Local code name
 * \param [in] param_name       Parameter name
 * \param [in] data_type        Parameter type
 * \param [in] value            Value
 *
 */

void
CWP_Param_set
(
 const char             *local_code_name,
 const char             *param_name,
 const CWP_Type_t        data_type,
 void                   *value
);


/**
 *
 * \brief Delete a control parameter.
 *
 * \param [in] local_code_name  Local code name
 * \param [in] param_name       Parameter name
 * \param [in] data_type        Parameter type
 *
 */

void
CWP_Param_del
(
 const char       *local_code_name,
 const char       *param_name,
 const CWP_Type_t  data_type
);


/**
 *
 * \brief Return the number of control parameters for the code \p code_name.
 *
 * /!\ Still unstable
 *
 * \param [in] code_name       Local or distant code name
 * \param [in] data_type       Parameter type
 *
 * return  Number of control parameters
 *
 */

int
CWP_Param_n_get
(
 const char             *code_name,
 const CWP_Type_t        data_type
);

/**
 *
 * \brief Return the list of control parameters for the code \p code_name.
 *
 * /!\ Still unstable
 *
 * \param [in]  code_name      Local or distant code name
 * \param [in]  data_type      Parameter type
 * \param [out] nParam         Number of parameters
 * \param [out] paramNames     Parameter names
 *
 *
 */

void
CWP_Param_list_get
(
 const char             *code_name,
 const CWP_Type_t        data_type,
 int                    *nParam,
 char                 ***paramNames
);

/**
 *
 * \brief Is this an existing control parameter for code \p code_name ?
 *
 * /!\ Still unstable
 *
 * \param [in] code_name      Local or distant code name
 * \param [in] param_name     Parameter name
 * \param [in] data_type      Parameter type
 *
 * return  1 : true / 0 : false
 *
 */

int
CWP_Param_is
(
 const char             *code_name,
 const char             *param_name,
 const CWP_Type_t        data_type
);

/**
 *
 * \brief Return the value of control parameter \p param_name on \p code_name.
 *
 * /!\ Still unstable
 *
 * \param [in]  code_name  Local or distant code name
 * \param [in]  param_name Parameter name
 * \param [in]  data_type  Parameter type
 * \param [out] value      Parameter value
 *
 */

void
CWP_Param_get
(
 const char       *code_name,
 const char       *param_name,
 const CWP_Type_t  data_type,
 void             *value
);

/**
 *
 * \brief Return the result of a reduce operation on a control parameter shared by multiple codes.
 *
 * /!\ Still unstable
 *
 * The parameter name has to be the same for all the codes involved in the reduction operation.
 *
 * \param [in]  op           Operation
 * \param [in]  param_name   Parameter name
 * \param [in]  data_type    Parameter type
 * \param [out] res          Result
 * \param [in]  n_code       Number of codes involved in the reduction operation
 * \param [in]  code_names   Names of codes involved in the reduction operation
 *
 */

void
CWP_Param_reduce
(
 const CWP_Op_t    op,
 const char       *param_name,
 const CWP_Type_t  data_type,
 void             *res,
 const int         n_code,
 const char      **code_names
);

/**
 *
 * \brief Lock access to local control parameters from a distant code.
 *
 * /!\ Still unstable
 *
 * \param [in]  code_name  Code to lock
 *
 */

void
CWP_Param_lock
(
const char *code_name
);

/**
 *
 * \brief Unlock access to local control parameters from a distant code.
 *
 * /!\ Still unstable
 *
 * \param [in]  code_name  Code to unlock
 *
 */

void
CWP_Param_unlock
(
const char *code_name
);

/*----------------------------------------------------------------------------*
 * Functions about data exchange                                              *
 *----------------------------------------------------------------------------*/

/**
 * \brief Initiate the sending of a global data array.
 *
 * \param [in] local_code_name  Local code name
 * \param [in] cpl_id           Coupling identifier
 * \param [in] global_data_id   Global data identifier
 * \param [in] s_send_entity    Data size
 * \param [in] send_stride      Constant stride value
 * \param [in] n_send_entity    Number of entities
 * \param [in] send_data        Pointer to data array
 *
 */

void
CWP_Global_data_issend
(
 const char     *local_code_name,
 const char     *cpl_id,
 const char     *global_data_id,
       size_t    s_send_entity,
       int       send_stride,
       int       n_send_entity,
       void     *send_data
);

/**
 * \brief Initiate the reception of a global data array.
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
CWP_Global_data_irecv
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
 * \brief Finalize the sending of a global data array.
 *
 * \param [in] local_code_name  Local code name
 * \param [in] cpl_id           Coupling identifier
 * \param [in] global_data_id   Global data identifier
 *
 */

void
CWP_Global_data_wait_issend
(
 const char     *local_code_name,
 const char     *cpl_id,
 const char     *global_data_id
);

/**
 * \brief Finalize the reception of a global data array.
 *
 * \param [in]  local_code_name  Local code name
 * \param [in]  cpl_id           Coupling identifier
 * \param [in]  global_data_id   Global data identifier
 *
 */

void
CWP_Global_data_wait_irecv
(
 const char     *local_code_name,
 const char     *cpl_id,
 const char     *global_data_id
);

/**
 * \brief Create a partitioned data exchange object
 *
 * \param [in] local_code_name  Local code name
 * \param [in] cpl_id           Coupling identifier
 * \param [in] part_data_id     Partitioned data identifier
 * \param [in] exch_type        Exchange type
 * \param [in] gnum_elt         Global ids
 * \param [in] n_elt            Number of elements per partition
 * \param [in] n_part           Number of partitions
 *
 */

void
CWP_Part_data_create
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
 * \brief Delete a partitioned data exchange object
 *
 * \param [in] local_code_name  Local code name
 * \param [in] cpl_id           Coupling identifier
 * \param [in] part_data_id     Partitioned data identifier
 *
 */

void
CWP_Part_data_del
(
 const char          *local_code_name,
 const char          *cpl_id,
 const char          *part_data_id
);

/**
 * \brief Initiate the sending of a partitioned data array.
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
CWP_Part_data_issend
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

 * \brief Initiate the reception of a partitioned data array.
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
CWP_Part_data_irecv
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
 * \brief Finalize the sending of a partitioned data array.
 *
 * \param [in] local_code_name  Local code name
 * \param [in] cpl_id           Coupling identifier
 * \param [in] part_data_id     Partitioned data identifier
 * \param [in] exch_id          Exchange identifier
 *
 */

void
CWP_Part_data_wait_issend
(
 const char    *local_code_name,
 const char    *cpl_id,
 const char    *part_data_id,
 const int      exch_id
);

/**
 * \brief Finalize the reception of a data array.
 *
 * \param [in] local_code_name  Local code name
 * \param [in] cpl_id           Coupling identifier
 * \param [in] part_data_id     Partitioned data identifier
 * \param [in] exch_id          Exchange identifier
 *
 */

void
CWP_Part_data_wait_irecv
(
 const char    *local_code_name,
 const char    *cpl_id,
 const char    *part_data_id,
 const int      exch_id
);


/**
 * \brief Set a generic high-order block to the interface mesh.
 *
 * \param [in]  local_code_name  Local code name
 * \param [in]  cpl_id           Coupling identifier
 * \param [in]  i_part           Partition identifier
 * \param [in]  block_id         Block identifier
 * \param [in]  n_elts           Number of elements
 * \param [in]  order            Element order
 * \param [in]  connec           Connectivity (size = *n_nodes_per_elt* * \p n_elts)
 * \param [in]  global_num       Global element ids (size = \p n_elts or NULL)
 *
 */

void
CWP_Mesh_interf_block_ho_set
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
 * \brief Get the properties of a generic high-order block of the interface mesh.
 *
 * \param [in]   local_code_name  Local code name
 * \param [in]   cpl_id           Coupling identifier
 * \param [in]   i_part           Partition identifier
 * \param [in]   block_id         Block identifier
 * \param [out]  n_elts           Number of elements
 * \param [out]  order            Element order
 * \param [out]  connec           Connectivity (size = *n_nodes_per_elt* * \p n_elts)
 * \param [out]  global_num       Global element ids (size = \p n_elts or NULL)
 *
 */

void
CWP_Mesh_interf_block_ho_get
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
 * \brief Define ho element ordering from the location in the (u, v, w) grid
 *
 * \param [in]  local_code_name  Local code name
 * \param [in]  cpl_id           Coupling identifier
 * \param [in]  block_type       Block type
 * \param [in]  order            Element order
 * \param [in]  n_nodes          Number of nodes per element
 * \param [in]  ijk_grid         User ordering to (u, v, w) grid (size = *elt_dim* * \p n_nodes)
 *
 */

void
CWP_Mesh_interf_ho_ordering_from_IJK_set
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
 * \brief Get the coupling spatial interpolation algorithm.
 *
 * \param [in] local_code_name  Local code name
 * \param [in] cpl_id           Coupling identifier
 *
 * \return                      Spatial interpolation method
 */

CWP_Spatial_interp_t
CWP_Cpl_spatial_interp_algo_get
(
 const char *local_code_name,
 const char *cpl_id
 );


/*****************************************************************************************************
 *                                                                                                   *
 *                              Not yet implemented                                                  *
 *                                                                                                   *
 *****************************************************************************************************/

// /*----------------------------------------------------------------------------
//  *
//  * Define specific options for ho elements
//  *
//  * parameters:
//  *   coupling_id     <-- coupling name
//  *   option          <-- option name, Choice between :
//  *                          - "opt_bbox_step"
//  *                              * Description : step of discretization used
//  *                                              to compute the optimized element
//  *                                              bounding boxes
//  *                                              -1 to deactivate this computation
//  *                              * Default     : 10
//  *   value           <-- option value
//  *
//  *----------------------------------------------------------------------------*/

// void cwipi_ho_options_set (const char *coupling_id,
//                            const char *option,
//                            const char *value);

// /*----------------------------------------------------------------------------
//  *
//  * Define ho element ordering from reference element (definition between 0 - 1)
//  *
//  *   coupling_id        <-- coupling name
//  *   t_elt              <-- element type
//  *   n_nodes            <-- number of nodes
//  *   coords             <-- node coordinates of reference element
//  *                                TODO: decrire ici les elements de reference
//  *
//  *----------------------------------------------------------------------------*/

// PAS GÉRÉ ACTUELLEMENT DANS PDM_HO_ORDERING
// void cwipi_ho_ordering_from_ref_elt_set (const char   *coupling_id,
//                                          const cwipi_element_t t_elt,
//                                          const int n_nodes,
//                                          const double *coords);





// /**
//  * \brief Initialize a translation displacement. <b>(Not implemented yet)</b>
//  *
//  * This function defines the translation direction. The movement is updated
//  * before the end of the current time step by \ref CWP_Cpl_trans_update
//  * Two coupled codes have to define the same properties. The distant code is always
//  * considered as a static interface.
//  *
//  * \param [in]  local_code_name  Local code name
//  * \param [in]  cpl_id           Coupling identifier
//  * \param [in]  vect             Direction
//  *
//  */

// void
// CWP_Cpl_trans_init
// (
//  const char      *local_code_name,
//  const char      *cpl_id,
//  const double     vect[3]
// );


// /**
//  * \brief Define the next time step position in a local mesh translation. <b>(Not implemented yet)</b>
//  *
//  * This function computes the next time step position from a relative distance. If
//  * it is a known position, spatial interpolation weights are not recomputed. Otherwise,
//  * spatial interpolation weights are computed from previous results.
//  *
//  *
//  * \param [in]  local_code_name  Local code name
//  * \param [in]  cpl_id           Coupling identifier
//  * \param [in]  dist             Relative distance from previous displacement
//  *
//  */

// void
// CWP_Cpl_trans_update
// (
//  const char      *local_code_name,
//  const char      *cpl_id,
//  const double     dist
// );


// /**
//  * \brief Initialize a rotation displacement. <b>(Not implemented yet)</b>
//  *
//  * This function defines the rotation properties. The movement is updated
//  * before the end of the current time step by \ref CWP_Cpl_rotation_update.
//  * Two coupled codes have to define the same properties. The distant code is always
//  * considered as a static interface.
//  *
//  * \param [in]  local_code_name  Local code name
//  * \param [in]  cpl_id           Coupling identifier
//  * \param [in]  vect             Direction
//  * \param [in]  center           Center
//  *
//  */

// void
// CWP_Cpl_rotation_init
// (
//  const char      *local_code_name,
//  const char      *cpl_id,
//  const double     vect[3],
//  const double     center[3]
// );


// /**
//  * \brief Define the next time step position in a local mesh rotation. <b>(Not implemented yet)</b>
//  *
//  * This function computes the next time step position from a relative angle. If
//  * it is a known position, the spatial interpolation weights are not reprocessed. Otherwise,
//  * the spatial interpolation weights are computed from previous results.
//  *
//  * \param [in]  cpl_id           Coupling identifier
//  * \param [in]  local_code_name  Local code name
//  * \param [in]  angle            Relative angle from previous displacement
//  *
//  */

// void
// CWP_Cpl_rotation_update
// (
//  const char      *local_code_name,
//  const char      *cpl_id,
//  const double     angle
// );


// /**
//  * \brief Set the storage properties of the spatial interpolation weights. <b>(Not implemented yet)</b>
//  *
//  * This functions activates the storage of the spatial interpolation weights in case of rotation
//  * or translation of the coupling interface.
//  *
//  * \param [in] cpl_id              Coupling identifier
//  * \param [in] local_code_name     Local code name
//  * \param [in] buffer_size         Size of buffer (Mo) on each coupling
//  *                                 communicator rank (same value for each)
//  * \param [in] disk_storage_size   Total size of disk storage when the buffer
//  *                                 is full (Mo) (Same value for each rank)
//  *
//  */

// void
// CWP_Cpl_storage_properties_set
// (
//  const char     *local_code_name,
//  const char     *cpl_id,
//  const int       buffer_size,
//  const int       disk_storage_size
// );

// /**
//  * \brief Return distance from each target to the source interface. <b>(Not implemented yet)</b>
//  *
//  * \param [in]  local_code_name  Local code name
//  * \param [in]  cpl_id           Coupling identifier
//  *
//  * \return               Distance
//  *
//  */

// const double *
// CWP_Computed_tgts_dist_to_spatial_interp_get
// (
//  const char *local_code_name,
//  const char *cpl_id
// );

//----------------------------------------------------------------------------
// Functions about exchange frequency
//----------------------------------------------------------------------------

// /**
//  * \brief Set receiving frequency. <b>(Not implemented yet)</b>
//  *
//  * This function set the receiving frequency. It must be used when
//  * the type of receiving frequency is <b> CWP_TIME_EXCH_N_TIME_STEP </b>
//  *
//  * \param [in]  local_code_name  Local code name
//  * \param [in]  cpl_id           Coupling identifier
//  * \param [in]  n_step           Frequency in steps number
//  *
//  */

// void
// CWP_Recv_freq_set
// (
//  const char      *local_code_name,
//  const char      *cpl_id,
//  const int        n_step
// );

// /**
//  * \brief Set the next receiving time. <b>(Not implemented yet)</b>
//  *
//  * It must be used when the type of receiving frequency is
//  * <b> CWP_TIME_EXCH_ASYNCHRONOUS </b>
//  *
//  * \param [in]  local_code_name  Local code name
//  * \param [in]  cpl_id           Coupling identifier
//  * \param [in]  next_time        Next receiving time
//  *
//  */

// void
// CWP_next_recv_time_set
// (
//  const char      *local_code_name,
//  const char      *cpl_id,
//  const double     next_time
// );


// /**
//  * \brief Set the coupling time step. <b>(Not implemented yet)</b>
//  *
//  * This function sets the coupling time step. It must be used when
//  * the type of receiving frequency is <b> CWP_TIME_EXCH_CPL_TIME_STEP </b>
//  *
//  * \param [in]  local_code_name  Local code name
//  * \param [in]  cpl_id           Coupling identifier
//  * \param [in]  next_time_step   Coupling time step
//  *
//  */

// void
// CWP_Cpl_time_step_set
// (
//  const char      *local_code_name,
//  const char      *cpl_id,
//  const int        next_time_step
// );


// /**
//  * \brief Exchange spatially interpolated fields. <b>(Not implemented yet)</b>
//  *
//  * This function exchanges the interpolated fields for each coupling depending
//  * on mode of time exchange \ref CWP_Time_exch_t.
//  *
//  * \param [in] local_code_name      Local code name
//  * \param [in] cpl_id               Coupling identifier
//  *
//  */

// void
// CWP_Field_exch
// (
//  const char *local_code_name,
//  const char *cpl_id
// );


// /**
//  * \brief Map a PDM mesh nodal as interface mesh. <b>(Not implemented yet)</b>
//  *
//  * \param [in] local_code_name   Local code name
//  * \param [in] cpl_id            Coupling identifier
//  * \param [in] i_part            Current partition
//  * \param [in] pdm_nodal         pdm nodal mesh
//  *
//  */


// void
// CWP_Mesh_interf_shared_pdm_nodal
// (
//  const char   *local_code_name,
//  const char   *cpl_id,
//  const int     i_part,
//  void         *pdm_nodal
// );


// /**
//  *
//  * \brief Get field data type.  <b>(Not implemented yet)</b>
//  *
//  * \param [in] local_code_name  Local code name
//  * \param [in] cpl_id           Coupling identifier
//  * \param [in] field_id         Field identifier
//  *
//  * \return                      Field data type
//  *
//  */

// CWP_Type_t
// CWP_Field_data_type_get
// (
//  const char      *local_code_name,
//  const char      *cpl_id         ,
//  const char      *field_id
// );


// /**
//  *
//  * \brief Set field gradient (optional). <b>(Not implemented yet)</b>
//  *
//  * \param [in] local_code_name Local code name
//  * \param [in] cpl_id          Coupling identifier
//  * \param [in] field_id        Field identifier
//  * \param [in] i_part          Current partition
//  * \param [in] order           Order
//  * \param [in] storage_type    Choice if data is set for the source or the target
//  * \param [in] data            Storage array (Mapping)
//  *
//  */

// void
// CWP_Field_gradient_data_set
// (
//  const char                *local_code_name,
//  const char                *cpl_id,
//  const char                *field_id,
//  const int                  i_part,
//  const int                  order,
//  const CWP_Field_storage_t  storage_type,
//  double                     data[]
// );

#include "fortran/new/cwp_cf.h"

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* __CWP_H__ */
