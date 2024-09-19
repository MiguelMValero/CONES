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

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <cstring>
#include <cstdlib>
#include <cstdarg>

/*----------------------------------------------------------------------------
 * BFT library headers
 *----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------
 * FVM library headers
 *----------------------------------------------------------------------------*/

#include "fvmc_parall.h"
#include "fvmc_nodal.h"
#include "cwp_priv.h"

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "cwp.h"
#include "cwipi_config.h"
#include "factory.hpp"
#include "codeProperties.hxx"
#include "codePropertiesDB.hxx"
#include "codePropertiesDB_i.hxx"
#include "couplingDB.hxx"
#include "couplingDB_i.hxx"
#include "coupling.hxx"
#include "coupling_i.hxx"
#include "commWithPart.hxx"
#include "commWithoutPart.hxx"
#include "field.hxx"
#include "pdm.h"
#include "pdm_printf.h"
#include "pdm_error.h"
#include "pdm_logging.h"
#include "pdm_mesh_nodal.h"
#include "spatialInterpClosestPoint.hxx"
#include "spatialInterpIntersection.hxx"
#include "spatialInterpLocationMeshLocation.hxx"
#include "spatialInterpIdentity.hxx"

#include "mesh.hxx"
#include "block.hxx"
#include "blockStd.hxx"
#include "blockFP.hxx"
#include "blockCP.hxx"
#include <algorithm>
#include <vector>

#include "pdm_mesh_nodal.h"
#include "pdm_ho_ordering.h"

#include <sys/stat.h>
#include <unistd.h>
#include <stdlib.h>

/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force back Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */

using namespace std;

/*============================================================================
 * Fortran function header
 *============================================================================*/

/*============================================================================
 * Local macro definitions
 *============================================================================*/

/*============================================================================
 * Type definitions
 *============================================================================*/

/*============================================================================
 * Global variable
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Output listing File (C printing)
 *----------------------------------------------------------------------------*/

static FILE* _cwipi_output_listing;

/*----------------------------------------------------------------------------
 * Output listing File (C printing)
 *----------------------------------------------------------------------------*/

const char *acceptable_properties[3] = {"n_neighbors", "polyfit_degree", "tolerance"};

/*============================================================================
 * Private function definitions
 *============================================================================*/

static bool
_is_active_rank
(
 const char *local_code_name
 )
{
   cwipi::CodePropertiesDB & propertiesDB = cwipi::CodePropertiesDB::getInstance();
   const cwipi::CodeProperties &properties = propertiesDB.codePropertiesGet(local_code_name);
   return properties.isActiveRank();
}


static cwipi::Coupling&
_cpl_get
(
 const char *local_code_name,
 const char *cpl_id
 )
{
  cwipi::CouplingDB & couplingDB =
    cwipi::CouplingDB::getInstance();

  cwipi::CodePropertiesDB & properties =
    cwipi::CodePropertiesDB::getInstance();

   const string &cpl_name_str = cpl_id;
   return couplingDB.couplingGet (properties.codePropertiesGet(string(local_code_name)),
                                 cpl_name_str);
}


void
_check_interface_spatial_interp_compatibility
(
 const CWP_Interface_t       entities_dim,
 const CWP_Spatial_interp_t  spatial_interp
)
{
  switch (entities_dim) {
    case CWP_INTERFACE_POINT:
    {
      if (spatial_interp != CWP_SPATIAL_INTERP_FROM_NEAREST_SOURCES_LEAST_SQUARES &&
          spatial_interp != CWP_SPATIAL_INTERP_FROM_NEAREST_TARGETS_LEAST_SQUARES &&
          spatial_interp != CWP_SPATIAL_INTERP_FROM_IDENTITY) {
        PDM_error(__FILE__, __LINE__, 0, "Only 'Nearest Neighbors' or 'Identity' spatial interpolation methods can be used with CWP_INTERFACE_POINT\n");
      }
      break;
    }
    case CWP_INTERFACE_LINEAR:
    {
      if (spatial_interp == CWP_SPATIAL_INTERP_FROM_INTERSECTION) {
        PDM_error(__FILE__, __LINE__, 0, "'Intersection' spatial interpolation is not yet implemented for CWP_INTERFACE_LINEAR\n");
      }
      break;
    }
    case CWP_INTERFACE_SURFACE:
    {
      break;
    }
    case CWP_INTERFACE_VOLUME:
    {
      break;
    }
    default:
    {
      PDM_error(__FILE__, __LINE__, 0, "Invalid coupling interface dimension (%d)\n", entities_dim);
    }
  }
}

void
_check_dof_location_spatial_interp_compatibility
(
 const CWP_Dof_location_t    dof_location,
 const CWP_Field_exch_t      exch_type,
 const CWP_Spatial_interp_t  spatial_interp
)
{
  switch (spatial_interp) {
    case CWP_SPATIAL_INTERP_FROM_NEAREST_SOURCES_LEAST_SQUARES:
    case CWP_SPATIAL_INTERP_FROM_NEAREST_TARGETS_LEAST_SQUARES:
    {
      // All dof location are accepted for both source and target fields
      break;
    }
    case CWP_SPATIAL_INTERP_FROM_INTERSECTION:
    {
      if (dof_location == CWP_DOF_LOCATION_USER) {
        PDM_error(__FILE__, __LINE__, 0,
                  "Field dofs must be located at mesh nodes or cell centers when using spatial interpolation from intersection\n");
      }

      if (dof_location == CWP_DOF_LOCATION_NODE) {
        PDM_error(__FILE__, __LINE__, 0,
                  "Spatial interpolation from intersection is not yet implemented for dof located at mesh nodes\n");
      }
      break;
    }
    case CWP_SPATIAL_INTERP_FROM_LOCATION_MESH_LOCATION_LOCATE_ALL_TGT:
    case CWP_SPATIAL_INTERP_FROM_LOCATION_MESH_LOCATION_OCTREE:
    case CWP_SPATIAL_INTERP_FROM_LOCATION_MESH_LOCATION_BOXTREE:
    {
      if (exch_type == CWP_FIELD_EXCH_RECV) {
        // All dof location are accepted for target field
      }
      else {
        if (dof_location == CWP_DOF_LOCATION_USER) {
          PDM_error(__FILE__, __LINE__, 0,
                    "Source field dofs must be located at mesh nodes or cell centers when using spatial interpolation from location\n");
        }
      }
      break;
    }
    case CWP_SPATIAL_INTERP_FROM_IDENTITY:
    {
      // All dof location are accepted for both source and target fields
      break;
    }
    default:
    {
      PDM_error(__FILE__, __LINE__, 0, "Invalid spatial interpolation method (%d)\n", spatial_interp);
    }
  }
}



// --> Fonctions privee qui n'appartiennent pas a l'API

 CWP_g_num_t*
 CWP_GlobalNumGet
 (
  const char  *local_code_name,
  const char  *cpl_id,
  const int    id_block,
  const int    i_part
 )
 {
   cwipi::Coupling& cpl = _cpl_get(local_code_name,cpl_id);

   return cpl.globalNumGet(id_block,i_part);

 }



MPI_Comm
CWP_Connectable_comm_get
(
  char* local_code_name
)
{

  if(_is_active_rank(local_code_name)){
    cwipi::CodePropertiesDB & properties = cwipi::CodePropertiesDB::getInstance();
    const cwipi::CodeProperties & localCodeProperties = properties.codePropertiesGet(string(local_code_name));
    MPI_Comm connecComm = localCodeProperties.connectableCommGet();

    return connecComm;
  }

  return MPI_COMM_NULL;
}

// <-- Fonctions privee qui n'appartiennent pas a l'API

/**
 *
 * \brief Intermediate function to write an output in the C file output
 *
 * \param [in] format     Format (as printf)
 * \param [in] arg_ptr    values to write (as printf)
 *
 * \return                writting status
 */

static int
_cwipi_print_with_c
(
 const char     *const format,
 va_list         arg_ptr
)
{
  return vfprintf(_cwipi_output_listing, format, arg_ptr);
}

/**
 *
 * \brief Intermediate function flush C file output
 *
 */

static int _cwipi_flush_output_listing(void)
{
  return fflush(_cwipi_output_listing);
}

/**
 *
 * \brief Return Coupling instance from it identifier
 *
 * \param [in] cpl_id     Coupling identifier
 * \param [in] fct        Function
 *
 * \return                Coupling instance from it identifier
 */

/*=============================================================================
 * Public function prototypes - with xml data
 *============================================================================*/

/*============================================================================*
 *                                                                            *
 *                     Public function prototypes                             *
 *                     --------------------------                             *
 *                          without xml                                       *
 *                          -----------                                       *
 *                                                                            *
 *============================================================================*/

/*----------------------------------------------------------------------------*
 * General functions                                                          *
 *----------------------------------------------------------------------------*/

/*============================================================================*
 *                                                                            *
 *                     Public function prototypes                             *
 *                     --------------------------                             *
 *                        without xml data                                    *
 *                        ----------------                                    *
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
)
{
  const int n_param_max_default = 100;
  const int str_size_max_default = 80;

  int n_param_max = n_param_max_default;
  int str_size_max = str_size_max_default;

  char* pPath;
  pPath = getenv ("CWP_N_PARAM_MAX");
  if (pPath!=NULL) {
    n_param_max = atoi(pPath);
  }

  if (pPath != NULL) {
   printf ("environment variable 'CWP_N_PARAM_MAX'"
           "is define to: %s\n",pPath);
  }
  pPath = getenv ("CWP_STR_SIZE_MAX");
  if (pPath!=NULL) {
    str_size_max = atoi(pPath);
  }
  if (pPath != NULL) {
   printf ("environment variable 'CWP_STR_SIZE_MAX'"
           "is set to: %s\n",pPath);
  }

  /*
   * Get application properties
   */
  cwipi::CodePropertiesDB & properties =
    cwipi::CodePropertiesDB::getInstance();

  int my_rank;  
  MPI_Comm_rank (global_comm, &my_rank);

  if (my_rank == 0) {

    PDM_printf("\ncwipi " CWP_DEF_VERSION " initializing\n");
    PDM_printf("------------------------\n\n");
    PDM_printf_flush();

  }

  /*
   * Builds application communicator
   */
  properties.init (global_comm,
                   n_code,
                   code_names,
                   is_active_rank,
                   n_param_max,
                   str_size_max,
                   intra_comms);

  /*
   * Create default parameters
   */
  MPI_Barrier(global_comm);

  for (int i = 0; i < n_code; i++) {
    if (is_active_rank == CWP_STATUS_ON) {
      const string &codeNameStr = code_names[i];
      properties.ctrlParamAdd <double> (codeNameStr, "time", 0.0); // WARNING : default first step 0.0
      properties.ctrlParamAdd <int> (codeNameStr, "state", CWP_STATE_IN_PROGRESS);
    }
  }

  MPI_Barrier(global_comm);

  /*
   * Create communication abstract factory
   */

  cwipi::Factory<cwipi::Communication, CWP_Comm_t> &factoryComm =
    cwipi::Factory<cwipi::Communication, CWP_Comm_t>::getInstance();

  factoryComm.Register<cwipi::CommWithPart>(CWP_COMM_PAR_WITH_PART);
  factoryComm.Register<cwipi::CommWithoutPart>(CWP_COMM_PAR_WITHOUT_PART);

  /*
   * Create spatial interpolation abstract factory
   */

  cwipi::Factory<cwipi::SpatialInterp, CWP_Spatial_interp_t> &factorySpatialInterp =
    cwipi::Factory<cwipi::SpatialInterp, CWP_Spatial_interp_t>::getInstance();

  factorySpatialInterp.Register<cwipi::SpatialInterpLocationMeshLocationOctree>(CWP_SPATIAL_INTERP_FROM_LOCATION_MESH_LOCATION_OCTREE);
  factorySpatialInterp.Register<cwipi::SpatialInterpLocationMeshLocationDbbtree>(CWP_SPATIAL_INTERP_FROM_LOCATION_MESH_LOCATION_BOXTREE);
  factorySpatialInterp.Register<cwipi::SpatialInterpLocationMeshLocationLocateAllTgt>(CWP_SPATIAL_INTERP_FROM_LOCATION_MESH_LOCATION_LOCATE_ALL_TGT);
  factorySpatialInterp.Register<cwipi::SpatialInterpIntersection>(CWP_SPATIAL_INTERP_FROM_INTERSECTION);
  factorySpatialInterp.Register<cwipi::SpatialInterpClosestSources>(CWP_SPATIAL_INTERP_FROM_NEAREST_SOURCES_LEAST_SQUARES);
  factorySpatialInterp.Register<cwipi::SpatialInterpClosestTargets>(CWP_SPATIAL_INTERP_FROM_NEAREST_TARGETS_LEAST_SQUARES);
  factorySpatialInterp.Register<cwipi::SpatialInterpIdentity>(CWP_SPATIAL_INTERP_FROM_IDENTITY);

  /*
   * Create block abstract factory
   */

  cwipi::Factory<cwipi::Block, CWP_Block_t> &factoryBlock =
    cwipi::Factory<cwipi::Block, CWP_Block_t>::getInstance();

  factoryBlock.Register<cwipi::BlockStd>(CWP_BLOCK_NODE);
  factoryBlock.Register<cwipi::BlockStd>(CWP_BLOCK_EDGE2);
  factoryBlock.Register<cwipi::BlockStd>(CWP_BLOCK_FACE_TRIA3);
  factoryBlock.Register<cwipi::BlockStd>(CWP_BLOCK_FACE_QUAD4);
  factoryBlock.Register<cwipi::BlockFP >(CWP_BLOCK_FACE_POLY);
  factoryBlock.Register<cwipi::BlockStd>(CWP_BLOCK_CELL_TETRA4);
  factoryBlock.Register<cwipi::BlockStd>(CWP_BLOCK_CELL_HEXA8);
  factoryBlock.Register<cwipi::BlockStd>(CWP_BLOCK_CELL_PRISM6);
  factoryBlock.Register<cwipi::BlockStd>(CWP_BLOCK_CELL_PYRAM5);
  factoryBlock.Register<cwipi::BlockCP >(CWP_BLOCK_CELL_POLY);
  factoryBlock.Register<cwipi::BlockHO >(CWP_BLOCK_FACE_TRIAHO);
  factoryBlock.Register<cwipi::BlockHO >(CWP_BLOCK_FACE_QUADHO);
  factoryBlock.Register<cwipi::BlockHO >(CWP_BLOCK_CELL_TETRAHO);
  factoryBlock.Register<cwipi::BlockHO >(CWP_BLOCK_CELL_HEXAHO);
  factoryBlock.Register<cwipi::BlockHO >(CWP_BLOCK_CELL_PRISMHO);
  factoryBlock.Register<cwipi::BlockHO >(CWP_BLOCK_CELL_PYRAMHO);

  MPI_Barrier(global_comm);

}

/*============================================================================*
 *                                                                            *
 *                   Other public function prototypes                         *
 *                   --------------------------------                         *
 *                                                                            *
 *============================================================================*/

/*----------------------------------------------------------------------------*
 * General functions                                                          *
 *----------------------------------------------------------------------------*/


/**
 *
 * \brief Finalize CWIPI.
 *
 */

void
CWP_Finalize
(
 void
)
{
  int flag = 0;

  MPI_Initialized(&flag);

  cwipi::CodePropertiesDB & properties =
    cwipi::CodePropertiesDB::getInstance();

  const MPI_Comm globalComm = properties.globalCommGet();

 // PDM_printf("CWP_Finalize\n");
  fflush(stdout);
  if (flag != 0) {
    PDM_printf_flush();
    MPI_Barrier(globalComm);
//    MPI_Comm oldFVMComm = fvmc_parall_get_mpi_comm();
  }


  properties.kill();

}

/*----------------------------------------------------------------------------*
 * Functions about current application properties                             *
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
 const char* local_code_name,
 const CWP_State_t state
)
{
  cwipi::CodePropertiesDB & properties =
    cwipi::CodePropertiesDB::getInstance();

  properties.ctrlParamSet<int>(string(local_code_name), "state", state);
}


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
 const char* code_name
)
{
  cwipi::CodePropertiesDB & properties =
    cwipi::CodePropertiesDB::getInstance();

  return static_cast<CWP_State_t> (
      properties.ctrlParamGet<int>(string(code_name), "state"));
}

/**
 * \brief Return the number of codes known by CWIPI.
 *
 * \return Number of codes
 *
 */

int
CWP_Codes_nb_get
(void
)
{
  cwipi::CodePropertiesDB & properties =
    cwipi::CodePropertiesDB::getInstance();

  return properties.codesNbGet();
}

/**
 * \brief Return the list of code names known by CWIPI.
 *
 * \return list of codes.
 */

const char **
CWP_Codes_list_get
(void
)
{
  cwipi::CodePropertiesDB & properties =
    cwipi::CodePropertiesDB::getInstance();
  return properties.codesListGet();
}

/**
 * \brief Return the number of local codes known by CWIPI.
 *
 * \return number of local codes.
 */

int
CWP_Loc_codes_nb_get
(void
)
{
  cwipi::CodePropertiesDB & properties =
    cwipi::CodePropertiesDB::getInstance();

  return properties.localCodesNbGet();
}


/**
 * \brief Return the list of local code names known by CWIPI.
 *
 * \return list of local codes.
 */

const char **
CWP_Loc_codes_list_get
(void
)
{
  cwipi::CodePropertiesDB & properties =
    cwipi::CodePropertiesDB::getInstance();
  return properties.localCodesListGet();
}


/**
 * \brief Begin code time step.
 *
 * \param [in] local_code_name  Local code name
 * \param [in]  current_time Current time
 *
 */

void
CWP_Time_step_beg
(
 const char* local_code_name,
 const double current_time
)
{
  cwipi::CodePropertiesDB & properties = cwipi::CodePropertiesDB::getInstance();

  int is_locked = properties.isLocked(local_code_name);

  if (!is_locked) {
    properties.lock(local_code_name);
  }

  properties.ctrlParamSet<double>(string(local_code_name),"time", current_time);

  if (!is_locked) {
    properties.unLock(local_code_name);
  }

  MPI_Comm connectable_comm = properties.connectableCommGet(local_code_name);
  MPI_Barrier(connectable_comm);

  cwipi::CouplingDB & couplingDB = cwipi::CouplingDB::getInstance();

  couplingDB.time_step_beg(properties.codePropertiesGet(local_code_name),
                           current_time);
}

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
)
{
  cwipi::CodePropertiesDB & properties = cwipi::CodePropertiesDB::getInstance();

  MPI_Comm connectable_comm = properties.connectableCommGet(local_code_name);
  MPI_Barrier(connectable_comm);

  cwipi::CouplingDB & couplingDB = cwipi::CouplingDB::getInstance();

  couplingDB.time_step_end(properties.codePropertiesGet(local_code_name));
}

/**
 * \brief Define a user structure associated to a code
 * 
 * This structure can be called into a callback 
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
)
{
  cwipi::CodePropertiesDB & properties =
  cwipi::CodePropertiesDB::getInstance();
  properties.userStructureSet(string(local_code_name), user_structure);  
}

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
)
{
  cwipi::CodePropertiesDB & properties =
  cwipi::CodePropertiesDB::getInstance();
  return properties.userStructureGet(string(local_code_name));  

}


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
)
{
  _cwipi_output_listing = output_file;
  PDM_printf_proxy_set(_cwipi_print_with_c);
  PDM_printf_flush_proxy_set(_cwipi_flush_output_listing);
}

//
//**
// * \brief Writing output to fortran file.
// *
// * This function set the file fortran logical unit for writing output.
// *
// * \param [in]  iunit        File fortan logical unit
// *
// */
//
//void
//PROCF (cwp_output_fortran_unit_set, CWP_OUTPUT_FORTRAN_UNIT_SET)
//(
// int *iunit
//);

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
)
{
  cwipi::CodePropertiesDB & properties =
    cwipi::CodePropertiesDB::getInstance();
  properties.dump();
}

/**
 * \brief Dump string of code properties.
 *
 */

int
CWP_Properties_str_dump
(
 char **char_out
)
{
  cwipi::CodePropertiesDB & properties =
    cwipi::CodePropertiesDB::getInstance();
  string str_out = properties.str_dump();

  // prepare output
  int size = str_out.length() + 1;
  *char_out = (char *) malloc(sizeof(char) * size);
  strcpy(*char_out, str_out.c_str());

  return size;
}

/*----------------------------------------------------------------------------*
 * General functions about coupling                                           *
 *----------------------------------------------------------------------------*/

/**
 * \brief Create a coupling object and define its properties.
 *
 * \param [in]  local_code_name     Local code name
 * \param [in]  cpl_id              Coupling identifier
 * \param [in]  coupled_code_name   Distant or local coupled code name
 * \param [in]  entities_dim        Dimensions of the entities
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
 const char                 *local_code_name,
 const char                 *cpl_id,
 const char                 *coupled_code_name,
 CWP_Interface_t             entities_dim,
 const CWP_Comm_t            comm_type,
 const CWP_Spatial_interp_t  spatial_interp,
 const int                   n_part,
 const CWP_Dynamic_mesh_t    displacement,
 const CWP_Time_exch_t       recv_freq_type
)
{
  cwipi::CouplingDB & couplingDB =
    cwipi::CouplingDB::getInstance();

  cwipi::CodePropertiesDB & properties =
    cwipi::CodePropertiesDB::getInstance();

  const string &coupling_name_str = cpl_id;
  const string &coupled_application_str = coupled_code_name;
  const string &local_application_str = local_code_name;

  _check_interface_spatial_interp_compatibility(entities_dim, spatial_interp);


 // printf("couplingDB.couplingCreate(proper\n");
  couplingDB.couplingCreate(properties.codePropertiesGet(local_application_str),
                            coupling_name_str,
                            properties.codePropertiesGet(coupled_application_str),
                            entities_dim,
                            comm_type,
                            spatial_interp,
                            n_part,
                            displacement,
                            recv_freq_type);
}


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
 )
{
  if (_is_active_rank(local_code_name)) {
    cwipi::Coupling& cpl = _cpl_get(local_code_name,cpl_id);
    cpl.barrier();
  }
}

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
)
{
  cwipi::CouplingDB & couplingDB =
    cwipi::CouplingDB::getInstance();

  cwipi::CodePropertiesDB & properties =
    cwipi::CodePropertiesDB::getInstance();

  // end (free) visualization
  if(_is_active_rank(local_code_name)){
    cwipi::Coupling& cpl = _cpl_get(local_code_name,cpl_id);
    cpl.visuEnd();
  }

  const string &cpl_id_str = cpl_id;

  couplingDB.couplingDel(properties.codePropertiesGet(string(local_code_name)),
                         cpl_id_str);
}

/**
 * \brief Get coupling communicator and coupling ranks.
 *
 * \param [in]  local_code_name      Local code name
 * \param [in]  cpl_id               Coupling identifier
 * \param [out] cpl_comm             Coupling communicator
 * \param [out] cpl_ranks            Coupling ranks
 *
 * \return Size of \ref cpl_ranks vector
 *
 */
// int
// CWP_Cpl_comm_get
// (
// const char *local_code_name,
// const char *cpl_id,
// MPI_Comm   *cpl_comm,
// int       **cpl_ranks
// )
// {
//   cwipi::Coupling& cpl = _cpl_get(local_code_name,cpl_id);
//   int size = cpl.commGet(cpl_comm,
//                          cpl_ranks);

//   return size;
// }

/**
 * \brief Exchange spatially interpolated fields. <b>(Not implemented yet)</b>
 *
 * This function exchanges the interpolated fields for each coupling depending
 * on mode of time exchange \ref CWP_Time_exch_t.
 *
 * \param [in] local_code_name      Local code name
 * \param [in] cpl_id               Coupling identifier
 *
 */

void
CWP_Field_exch
(
 const char *local_code_name,
 const char *cpl_id
)
{
  //TODO: Voir comment enchainer les appels, voir comment prendre en compte l'interpolation temporelle
  cwipi::Coupling& cpl = _cpl_get(local_code_name,cpl_id);
  CWP_UNUSED(cpl);
  printf("CWP_Field_exch : Not implemented yet\n");
  abort();

  // cpl.exchange();
}


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
)
{
  cwipi::Coupling& cpl = _cpl_get(local_code_name,cpl_id);

  const string &field_name_str = field_id;
  return cpl.computedTargetsBcastEnable(field_name_str);
}


/**
 *
 * \brief Return the number of uncomputed targets.
 *
 * \param [in] local_code_name  Local code name
 * \param [in] cpl_id           Coupling identifier
 * \param [in] field_id         Field identifier
 * \param [in] i_part           Current partition
 *
 * \return                Number of uncomputed targets
 */


int
CWP_N_uncomputed_tgts_get
(
 const char *local_code_name,
 const char *cpl_id,
 const char *field_id,
 const int   i_part
)
{
  cwipi::Coupling& cpl = _cpl_get(local_code_name,cpl_id);

  const string &field_name_str = field_id;
  return cpl.nUncomputedTargetsGet(field_name_str, i_part);
}


/**
 *
 * \brief Return uncomputed targets.
 *
 * \param [in] local_code_name  Local code name
 * \param [in] cpl_id           Coupling identifier
 * \param [in] field_id         Field identifier
 * \param [in] i_part           Current partition
 *
 * \return                Uncomputed targets
 */

const int *
CWP_Uncomputed_tgts_get
(
  const char *local_code_name,
  const char *cpl_id,
  const char *field_id,
  const int   i_part
)
{
  cwipi::Coupling& cpl = _cpl_get(local_code_name,cpl_id);

  const string &field_name_str = field_id;
  return cpl.uncomputedTargetsGet(field_name_str, i_part);
}


/**
 *
 * \brief Return the number of computed targets.
 *
 * \param [in] local_code_name  Local code name
 * \param [in] cpl_id           Coupling identifier
 * \param [in] field_id         Field identifier
 * \param [in] i_part           Current partition
 *
 * \return                Number of computed targets
 */

int
CWP_N_computed_tgts_get
(
  const char *local_code_name,
  const char *cpl_id,
  const char *field_id,
  const int   i_part
)
{
  cwipi::Coupling& cpl = _cpl_get(local_code_name,cpl_id);

  const string &field_name_str = field_id;
  return cpl.nComputedTargetsGet(field_name_str, i_part);
}

/**
 *
 * \brief Return computed targets.
 *
 * \param [in] local_code_name  Local code name
 * \param [in] cpl_id           Coupling identifier
 * \param [in] field_id         Field identifier
 * \param [in] i_part           Current partition
 *
 * \return                Computed targets
 */

const int *
CWP_Computed_tgts_get
(
  const char *local_code_name,
  const char *cpl_id,
  const char *field_id,
  const int   i_part
)
{
  cwipi::Coupling& cpl = _cpl_get(local_code_name,cpl_id);

  const string &field_name_str = field_id;
  return cpl.computedTargetsGet(field_name_str, i_part);
}


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
)
{
  cwipi::Coupling& cpl = _cpl_get(local_code_name,cpl_id);

  const string &field_name_str = field_id;
  return cpl.involvedSourcesBcastEnable(field_name_str);
}

/**
 *
 * \brief Return the number of source elements.
 *
 * \param [in] local_code_name  Local code name
 * \param [in] cpl_id           Coupling identifier
 * \param [in] field_id         Field identifier
 * \param [in] i_part           Current partition
 *
 * \return                Number of source elements
 */

int
CWP_N_involved_srcs_get
(
  const char *local_code_name,
  const char *cpl_id,
  const char *field_id,
  const int   i_part
)
{
  cwipi::Coupling& cpl = _cpl_get(local_code_name,cpl_id);

  const string &field_name_str = field_id;
  return cpl.nInvolvedSourcesGet(field_name_str, i_part);
}

/**
 *
 * \brief Return source elements.
 *
 * \param [in] local_code_name  Local code name
 * \param [in] cpl_id           Coupling identifier
 * \param [in] field_id         Field identifier
 * \param [in] i_part           Current partition
 *
 * \return                Source elements
 */

const int *
CWP_Involved_srcs_get
(
  const char *local_code_name,
  const char *cpl_id,
  const char *field_id,
  const int   i_part
)
{
  cwipi::Coupling& cpl = _cpl_get(local_code_name,cpl_id);

  const string &field_name_str = field_id;
  return cpl.involvedSourcesGet(field_name_str, i_part);
}


/**
 * \brief Return distance from each target to the source interface. <b>(Not implemented yet)</b>
 *
 * \param [in]  local_code_name  Local code name
 * \param [in]  cpl_id           Coupling identifier
 *
 * \return               Distance
 * 
 * A supprimer ???
 *
 */


const double *
CWP_Computed_tgts_dist_to_spatial_interp_get
(
  const char *local_code_name,
  const char *cpl_id
)

{
  CWP_UNUSED (local_code_name);
  CWP_UNUSED (cpl_id);
  printf("CWP_Computed_tgts_dist_to_spatial_interp_get : Not implemented yet\n");
  abort();
  return nullptr;
  // TODO
}


/*----------------------------------------------------------------------------*
 * Functions about exchange frequency                                         *
 *----------------------------------------------------------------------------*/



// void
// CWP_recv_freq_set
// (
//  const char                 *cpl_id,
//  const int                   n_step
// )
// {
//   cwipi::Coupling& cpl = _cpl_get(local_code_name,cpl_id);
//   cpl.recvFreqSet(n_step);
// }


/**
 * \brief Set the next receiving time.
 *
 * It must be used when the type of receiving frequency is
 * \ref CWP_TIME_EXCH_ASYNCHRONOUS
 *
 * \param [in]  local_code_name  Local code name
 * \param [in]  cpl_id           Coupling identifier
 * \param [in]  next_time        Next receiving time
 *
 */

 void
 CWP_next_recv_time_set
 (const char     *local_code_name,
  const char     *cpl_id,
  const double    next_time
 )
 {
   cwipi::Coupling& cpl = _cpl_get(local_code_name,cpl_id);
   // cpl.recvNextTimeSet(next_time);//TO REMOVE
   CWP_UNUSED(cpl);
   CWP_UNUSED(next_time);
 }

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
)
{

  cwipi::Coupling& cpl = _cpl_get(local_code_name,cpl_id);

  cpl.spatialInterpWeightsCompute ();

}


/**
 * \brief Set a property of the spatial interpolation algorithm.
 *
 * Use "n_neighbors" and "polyfit_degree" for the nearest neighbors
 * algorithm. Use "tolerance" for the location algorithm.
 *
 * \param [in]  local_code_name  Local code name
 * \param [in]  cpl_id           Coupling identifier
 * \param [in]  property_name    Name of the property
 * \param [in]  property_type    Type of the property
 * \param [in]  property_value   Value of the property
 *
 */

void
CWP_Spatial_interp_property_set (
 const char       *local_code_name,
 const char       *cpl_id,
 const char       *property_name,
 const CWP_Type_t  property_type,
 const char       *property_value
)
{
  int is_property = 0;

  for (int i = 0 ; i < 3; i++) {
    if (strcmp(property_name, acceptable_properties[i]) == 0) {
      is_property = 1;
      break;
    }
  }

  if (is_property) {

    cwipi::Coupling& cpl = _cpl_get(local_code_name,cpl_id);

    switch (property_type) {
      case CWP_DOUBLE: {
        cpl.spatialInterpPropertyDoubleSet(string(property_name),
                                           atof(property_value));
        break;
      }
      case CWP_INT: {
        cpl.spatialInterpPropertyIntSet(string(property_name),
                                        atoi(property_value));
        break;
      }
      default: {
        printf("Warning: invalid property type %d\n", property_type);
      }
    }

  } else {
    printf("Warning: invalid property name '%s'\n", property_name);
  }
}


// void
// CWP_spatial_interp_update
// (
//  const char     *cpl_id,
//  const int       storage_id
// )
// {
// TODO
// }


// void
// CWP_spatial_interp_properties_set
// (
//  const char     *cpl_id,
//  const char     *fmt,
//                   ...
// )
// {
//   cwipi::Coupling& cpl = _cpl_get(local_code_name,cpl_id);

//   va_list pa;
//   va_start(pa, fmt);
//   cpl.spatialInterpPropertiesSet(fmt, &pa);
//   va_end(pa);
// }

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
)
{
  if(_is_active_rank(local_code_name)){
    cwipi::Coupling& cpl = _cpl_get(local_code_name,cpl_id);
    cpl.visuSet(freq, format, format_option);
  }
}

/*----------------------------------------------------------------------------*
 * Functions about User target points                                         *
 *----------------------------------------------------------------------------*/


/**
 * \brief Setting user target points.
 *
 * This function must be called if the degrees of freedom locations are
 * \ref CWP_DOF_LOCATION_USER
 *
 * \param [in]  local_code_name  Local code name
 * \param [in]  cpl_id           Coupling identifier
 * \param [in]  i_part           Current partition
 * \param [in]  n_pts            Number of points
 * \param [in]  coord            Coordinates (size = 3 * n_pts)
 * \param [in]  g_num            global number or NUL (size = n_pts)
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
)
{
  if(_is_active_rank(local_code_name)){
    cwipi::Coupling& cpl = _cpl_get(local_code_name,cpl_id);
    cpl.userTargetSet(i_part, n_pts, coord, global_num);
  }
}

/*----------------------------------------------------------------------------*
 * Functions about Mesh                                                    *
 *----------------------------------------------------------------------------*/

/**
 * \brief Finalize interface mesh.
 *
 * This function computes the global numbers of mesh entities if they are
 * not provided.
 *
 * \param [in]  local_code_name  Local code name
 * \param [in]  cpl_id           Coupling identifier
 *
 */

void
CWP_Mesh_interf_finalize
(
 const char                 *local_code_name,
 const char                 *cpl_id
)
{
  if(_is_active_rank(local_code_name)){
    cwipi::Coupling& cpl = _cpl_get(local_code_name,cpl_id);
    cpl.meshFinalize();
  }
}



/**
 * \brief Set vertices.
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
CWP_Mesh_interf_vtx_set
(
 const char                 *local_code_name,
 const char                 *cpl_id,
 const int                   i_part,
 const int                   n_pts,
 double                      coord[],
 CWP_g_num_t                 global_num[]
)
{
  if(_is_active_rank(local_code_name)){
    cwipi::Coupling& cpl = _cpl_get(local_code_name,cpl_id);
    cpl.meshVtcsSet(i_part,
                    n_pts,
                    coord,
                    global_num);
  }
}

/**
 * \brief Add a connectivity block to the interface mesh.
 *
 * \param [in]  local_code_name  Local code name
 * \param [in]  cpl_id           Coupling identifier
 * \param [in]  block_type       Block type
 *
 * \return block identifier (-1 if not added)
 */

int
CWP_Mesh_interf_block_add
(
 const char           *local_code_name,
 const char           *cpl_id,
 const CWP_Block_t     block_type
)
{
  int block_id = -1;
  if(_is_active_rank(local_code_name)){
    cwipi::Coupling& cpl = _cpl_get(local_code_name,cpl_id);
    block_id = cpl.meshBlockAdd(block_type);
  }
  return block_id;
}


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
 *   - tetrahedron (\ref CWP_BLOCK_CELL_TETRA4) :
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
 *   - pyramid (\ref CWP_BLOCK_CELL_PYRAM5) :
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
 *  -  hexaedron (\ref CWP_BLOCK_CELL_HEXA8) :
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
 * \param [in]  connec           Connectivity (size = n_vertex_elt * n_elts)
 * \param [in]  global_num       Pointer to global element number (or NULL)
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
 CWP_g_num_t       global_num[]
)
{
  if(_is_active_rank(local_code_name)){
    cwipi::Coupling& cpl = _cpl_get(local_code_name,cpl_id);
    cpl.meshStdBlockSet(i_part,
                        block_id,
                        n_elts,
                        connec,
                        global_num);
  }
}

/**
 * \brief Get the properties of a standard block of the interface mesh.
 *
 * \param [in]  local_code_name  Local code name
 * \param [in]  cpl_id           Coupling identifier
 * \param [in]  i_part           Partition identifier
 * \param [in]  block_id         Block identifier
 * \param [out]  n_elts           Number of elements
 * \param [out]  connec           Connectivity (size = n_vertex_elt * n_elts)
 * \param [out]  global_num       Pointer to global element number (or NULL)
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
)
{
  cwipi::Coupling& cpl = _cpl_get(local_code_name,cpl_id);
  cpl.meshStdBlockGet (i_part,
                       block_id,
                       n_elts,
                       connec,
                       global_num);

}


 /**
  * \brief Get the standard block type
  *
  * \param [in]  block_id    Block identifier
  *
  * \return block type
  */

CWP_Block_t
CWP_std_block_type_get
(
 const char             *local_code_name,
 const char             *cpl_id,
 const int               block_id
)
{
  cwipi::Coupling& cpl = _cpl_get(local_code_name,cpl_id);
  CWP_Block_t block_type = cpl.meshStdBlockTypeGet(block_id);

  return block_type;
}


/**
 * \brief Set the connectivity of a polygon block in a interface mesh partition.
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
)
{

  if(_is_active_rank(local_code_name)){
    cwipi::Coupling& cpl = _cpl_get(local_code_name,cpl_id);
    cpl.meshFPolyBlockSet(i_part,
                          block_id,
                          n_elts,
                          connec_idx,
                          connec,
                          global_num);
  }
}


/**
 * \brief Get the properties of a polygon block of the interface mesh partition.
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
)
{
  cwipi::Coupling& cpl = _cpl_get(local_code_name,cpl_id);
  cpl.meshFPolyBlockGet(i_part,
                        block_id,
                        n_elts,
                        connec_idx,
                        connec,
                        global_num);
  
}


/**
 * \brief Get the properties of a polyhedron block of the interface mesh partition..
 *
 * \param [in]  local_code_name   Local code name
 * \param [in]  cpl_id            Coupling identifier
 * \param [in]  i_part            Current partition
 * \param [in]  block_id          Block identifier
 * \param [in]  n_elts            Number of elements
 * \param [in]  n_faces           Number of faces
 * \param [in]  connec_faces_idx  Polyhedron face to vertex index
 *                                (\p connec_faces_idx[0] = 0 and
 *                                 size = max(\p connec_cells) + 1)
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
)
{
  cwipi::Coupling& cpl = _cpl_get(local_code_name,cpl_id);
  cpl.meshCPolyBlockSet(i_part,
                       block_id,
                       n_elts,
                       n_faces,
                       connec_faces_idx,
                       connec_faces    ,
                       connec_cells_idx,
                       connec_cells    ,
                       global_num);
}


/**
 * \brief Get the properties of a polyhedron block of the interface mesh partition.
 *
 * \param [in]  local_code_name   Local code name
 * \param [in]  cpl_id            Coupling identifier
 * \param [in]  i_part            Current partition
 * \param [in]  block_id          Block identifier
 * \param [in]  n_elts            Number of elements
 * \param [in]  n_faces           Number of faces
 * \param [in]  connec_faces_idx  Polyhedron face to vertex index
 *                                (\p connec_faces_idx[0] = 0 and
 *                                 size = max(\p connec_cells) + 1)
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
)
{
  cwipi::Coupling& cpl = _cpl_get(local_code_name,cpl_id);
  cpl.meshCPolyBlockGet(i_part,
                         block_id,
                         n_elts,
                         n_faces,
                         connec_faces_idx,
                         connec_faces    ,
                         connec_cells_idx,
                         connec_cells    ,
                         global_num);

}


/**
 * \brief Define the interface mesh from a cell to face connectivity.
 *
 * \param [in]  local_code_name   Local code name
 * \param [in]  cpl_id            Coupling identifier
 * \param [in]  i_part            Current partition
 * \param [in]  n_cells           Number of cells
 * \param [in]  cell_face_idx     Polyhedron to face index
 *                                (\p cell_face_idx[0] = 0 and
 *                                 size = \p n_elts + 1)
 * \param [in]  cell_face         Cell to face connectivity
 *                                (size = \p cell_face_idx[\p n_elts])
 * \param [in]  n_faces           Number of faces
 * \param [in]  face_vtx_idx      Polyhedron face to vertex index
 *                                (\p face_vtx_idx[0] = 0 and
 *                                 size = \p n_faces + 1)
 * \param [in]  face_vtx          Face to vertex connectivity
 *                                (size = \p face_vtx_idx[\p n_elts])
 * \param [in]  global_num        Global element ids (size = \p n_elts or NULL)
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
)
{
  cwipi::Coupling& cpl = _cpl_get(local_code_name,cpl_id);
  cpl.meshFromCellFaceSet(i_part,
                          n_cells,
                          cell_face_idx,
                          cell_face,
                          n_faces,
                          face_vtx_idx,
                          face_vtx,
                          global_num);
}


/**
 * \brief Define the surface interface mesh from a face to edge connectivity.
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
 * \param [in]  parent_num        Pointer to parent element number (or NULL)
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
       CWP_g_num_t  parent_num[]
)
{
  cwipi::Coupling& cpl = _cpl_get(local_code_name,cpl_id);
  cpl.meshFromFacesEdgeSet(i_part,
                           n_faces,
                           face_edge_idx,
                           face_edge,
                           n_edges,
                           edge_vtx,
                           parent_num);
}


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
)
{
  cwipi::Coupling& cpl = _cpl_get(local_code_name,cpl_id);
  cpl.meshFromFacesVtxSet(i_part,
                          n_faces,
                          face_vtx_idx,
                          face_vtx,
                          global_num);
}


/**
 * \brief Delete interface mesh.
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
)
{
   cwipi::Coupling& cpl = _cpl_get(local_code_name,cpl_id);
   cpl.meshDel();
}

/*
void
CWP_Mesh_interf_shared_pdm_nodal
(
 const char   *local_code_name,
 const char   *cpl_id,
 const int     i_part,
 fvmc_nodal_t *fvmc_nodal
)
{
   cwipi::Coupling& cpl = _cpl_get(local_code_name,cpl_id);
   cpl.fvmcNodalShared(i_part,
                       fvmc_nodal);
}
*/
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
 * \param [in]  n_component     Number of component
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
)
{
  cwipi::Coupling& cpl = _cpl_get(local_code_name,cpl_id);

  CWP_Spatial_interp_t spatial_interp = CWP_Cpl_spatial_interp_algo_get(local_code_name,
                                                                        cpl_id);

  _check_dof_location_spatial_interp_compatibility(dof_location,
                                                   exch_type,
                                                   spatial_interp);

  cpl.fieldCreate(field_id,
                  data_type,
                  storage,
                  n_component,
                  dof_location,
                  exch_type,
                  visu_status);
}

/*----------------------------------------------------------------------------*
 * Functions about user interpolation                                         *
 *----------------------------------------------------------------------------*/

void
CWP_Field_python_object_set
(
 const char *local_code_name,
 const char *cpl_id,
 const char *field_id,
       void *python_object
 )
 {
  cwipi::Coupling& cpl = _cpl_get(local_code_name,cpl_id);
  cpl.fieldPythonObjectSet(field_id,
                           python_object);
 }

/**
 *
 * \brief Unsetting of an user interpolation.
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
)
{
  cwipi::Coupling& cpl = _cpl_get(local_code_name,cpl_id);
  cpl.interpFunctionUnSet(field_id);
}

// Specific function for Fortran interface

void
CWP_Field_interp_function_f_unset
(
 const char                 *local_code_name,
 const char                 *cpl_id,
 const char                 *field_id
)
{
  cwipi::Coupling& cpl = _cpl_get(local_code_name,cpl_id);
  cpl.interpFunctionFUnSet(field_id);
}

// Specific function for Python interface

void
CWP_Field_interp_function_p_unset
(
 const char                 *local_code_name,
 const char                 *cpl_id,
 const char                 *field_id
 )
{
  cwipi::Coupling& cpl = _cpl_get(local_code_name,cpl_id);
  cpl.interpFunctionPUnset(field_id);
}

/**
 *
 * \brief Setting of an user interpolation from location.
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
CWP_Field_interp_function_set
(
 const char                 *local_code_name,
 const char                 *cpl_id,
 const char                 *field_id,
 CWP_Interp_function_t       fct
 )
{
  cwipi::Coupling& cpl = _cpl_get(local_code_name,cpl_id);
  cpl.interpFunctionSet(field_id,fct);
}

// Specific function for Fortran interface

void
CWP_Field_interp_function_f_set
(
 const char                 *local_code_name,
 const char                 *cpl_id,
 const char                 *field_id,
 CWP_Interp_function_t       fct
 )
{
  cwipi::Coupling& cpl = _cpl_get(local_code_name,cpl_id);
  cpl.interpFunctionFSet(field_id,fct);
}

// Specific function for Python interface

void
CWP_Field_interp_function_p_set
(
 const char                 *local_code_name,
 const char                 *cpl_id,
 const char                 *field_id,
 CWP_Interp_function_p_t     fct
 )
{
  cwipi::Coupling& cpl = _cpl_get(local_code_name,cpl_id);
  cpl.interpFunctionPSet(field_id,fct);
}

/**
 *
 * \brief Set field data.
 *
 * \param [in] local_code_name   Local code name
 * \param [in] cpl_id            Coupling identifier
 * \param [in] field_id          Field identifier
 * \param [in] i_part            Current partition
 * \param [in] data_type         Choice if data is set for the source or the target
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
)
{
  if(_is_active_rank(local_code_name)){
    cwipi::Coupling& cpl = _cpl_get(local_code_name, cpl_id);
    cpl.fieldDataSet(field_id, i_part, map_type, data);
  }
}

/**
 *
 * \brief Get field data.
 *
 * \param [in]  local_code_name   Local code name
 * \param [in]  cpl_id            Coupling identifier
 * \param [in]  field_id          Field identifier
 * \param [in]  i_part            Current partition
 * \param [in]  data_type         Choice if data is set for the source or the target
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
)
{
  if(_is_active_rank(local_code_name)){
    cwipi::Coupling& cpl = _cpl_get(local_code_name, cpl_id);
    cpl.fieldDataGet(field_id, i_part, map_type, (void **) data);
  }
}


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
 const char                  *local_code_name,
 const char                  *cpl_id,
 const char                  *field_id
)
{
  cwipi::Coupling& cpl = _cpl_get(local_code_name,cpl_id);
  return cpl.fieldDofLOcationGet(field_id);
}

/*
 CWP_Type_t
 CWP_Field_type_get
 (
  const char                  *local_code_name,
  const char                  *cpl_id,
  const char                  *field_id
 )
 {
   cwipi::Coupling& cpl = _cpl_get(local_code_name,cpl_id);
   return cpl.fieldTypeGet(field_id);
 }*/


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
 const char                  *local_code_name,
 const char                  *cpl_id,
 const char                  *field_id
)
{
  cwipi::Coupling& cpl = _cpl_get(local_code_name,cpl_id);
  return cpl.fieldStorageGet(field_id);
}

/**
 *
 * \brief Get number of field components.
 *
 * \param [in]  local_code_name           Local code name
 * \param [in]  cpl_id                    Coupling identifier
 * \param [in]  field_id                  Field identifier
 *
 */

int
CWP_Field_n_components_get
(
 const char             *local_code_name,
 const char             *cpl_id,
 const char             *field_id
)
{
  cwipi::Coupling& cpl = _cpl_get(local_code_name,cpl_id);
  return cpl.fieldNComponentGet(field_id);
}

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
)
{
  cwipi::Coupling& cpl = _cpl_get(local_code_name,cpl_id);
  return cpl.fieldNDOFGet(field_id, i_part);
}

/**
 *
 * \brief Get spatial interpolation source data.
 *
 * \param [in]  local_code_name           Local code name
 * \param [in]  cpl_id                    Coupling identifier
 * \param [in]  field_id                  Field identifier
 * \param [in]  i_part                    Partition identifier
 * \param [out] n_elt_src                 Number of source elements
 * \param [out] src_to_tgt_idx            Source elements to target elements index
 *
 */

void
CWP_Field_src_data_properties_get
(
 const char   *local_code_name,
 const char   *cpl_id,
 const char   *field_id,
 int           i_part,
 int          *n_elt_src,
 int         **src_to_tgt_idx
)
{
  cwipi::Coupling& cpl = _cpl_get(local_code_name,cpl_id);

  int           n_part_src = 0;
  int          *tmp_n_elt_src       = NULL;
  int         **tmp_src_to_tgt_idx  = NULL;
  CWP_g_num_t **tmp_src_to_tgt_gnum = NULL;  // UNUSED
  cpl.src_data_get(field_id,
                   &n_part_src,
                   &tmp_n_elt_src,
                   &tmp_src_to_tgt_idx,
                   &tmp_src_to_tgt_gnum);

  *n_elt_src       = tmp_n_elt_src[i_part];
  *src_to_tgt_idx  = tmp_src_to_tgt_idx[i_part];
}

/**
 *
 * \brief Get spatial interpolation target data.
 *
 * \param [in]  local_code_name           Local code name
 * \param [in]  cpl_id                    Coupling identifier
 * \param [in]  field_id                  Field identifier
 * \param [in]  i_part                    Partition identifier
 * \param [out] n_elt_tgt                 Number of target elements
 * \param [out] n_referenced_tgt          Number of referenced target elements
 * \param [out] referenced_tgt            Referenced target elements
 * \param [out] tgt_come_from_src_idx     Target to origin source elements index
 *
 */

void
CWP_Field_tgt_data_properties_get
(
 const char   *local_code_name,
 const char   *cpl_id,
 const char   *field_id,
 int           i_part,
 int          *n_elt_tgt,
 int          *n_referenced_tgt,
 int         **referenced_tgt,
 int         **tgt_come_from_src_idx
)
{
  cwipi::Coupling& cpl = _cpl_get(local_code_name,cpl_id);

  int           n_part_tgt = 0;
  int          *tmp_n_elt_tgt             = NULL;
  int          *tmp_n_referenced_tgt      = NULL;
  int         **tmp_referenced_tgt        = NULL;
  int         **tmp_tgt_come_from_src_idx = NULL;
  CWP_g_num_t **tmp_tgt_come_from_src     = NULL; // UNUSED
  cpl.tgt_data_get(field_id,
                   &n_part_tgt,
                   &tmp_n_elt_tgt,
                   &tmp_n_referenced_tgt,
                   &tmp_referenced_tgt,
                   &tmp_tgt_come_from_src_idx,
                   &tmp_tgt_come_from_src);

  *n_elt_tgt             = tmp_n_elt_tgt[i_part];
  *n_referenced_tgt      = tmp_n_referenced_tgt[i_part];
  *referenced_tgt        = tmp_referenced_tgt[i_part];
  *tgt_come_from_src_idx = tmp_tgt_come_from_src_idx[i_part];
}

/**
 *
 * \brief Get spatial interpolation weights (location algorithm).
 *
 * \param [in]  local_code_name           Local code name
 * \param [in]  cpl_id                    Coupling identifier
 * \param [in]  field_id                  Field identifier
 * \param [in]  i_part                    Partition identifier
 * \param [out] weights                   Spatial interpolation weights
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
)
{
  cwipi::Coupling& cpl = _cpl_get(local_code_name,cpl_id);

  CWP_Spatial_interp_t spatial_interp_algorithm = cpl.spatialInterpAlgoGet();

  if ((spatial_interp_algorithm != CWP_SPATIAL_INTERP_FROM_LOCATION_MESH_LOCATION_LOCATE_ALL_TGT) &&
      (spatial_interp_algorithm != CWP_SPATIAL_INTERP_FROM_LOCATION_MESH_LOCATION_OCTREE)         &&
      (spatial_interp_algorithm != CWP_SPATIAL_INTERP_FROM_LOCATION_MESH_LOCATION_BOXTREE)) {
    PDM_error(__FILE__, __LINE__, 0, "Getter unavailable for spatial interpolation algorithm %d\n", spatial_interp_algorithm);
  }

  int    **tmp_weights_idx = NULL; // UNUSED
  double **tmp_weights     = NULL;
  cpl.weight_get(field_id,
                 &tmp_weights_idx,
                 &tmp_weights);

  *weights     = tmp_weights[i_part];
}

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
)
{
  cwipi::Coupling& cpl = _cpl_get(local_code_name,cpl_id);

  CWP_Spatial_interp_t spatial_interp_algorithm = cpl.spatialInterpAlgoGet();

  if ((spatial_interp_algorithm != CWP_SPATIAL_INTERP_FROM_LOCATION_MESH_LOCATION_LOCATE_ALL_TGT) &&
      (spatial_interp_algorithm != CWP_SPATIAL_INTERP_FROM_LOCATION_MESH_LOCATION_OCTREE)         &&
      (spatial_interp_algorithm != CWP_SPATIAL_INTERP_FROM_LOCATION_MESH_LOCATION_BOXTREE)) {
    PDM_error(__FILE__, __LINE__, 0, "Getter unavailable for spatial interpolation algorithm %d\n", spatial_interp_algorithm);
  }

  double **tmp_points_coords           = NULL;
  double **tmp_points_uvw              = NULL;
  double **tmp_points_dist2            = NULL;
  double **tmp_points_projected_coords = NULL;
  cpl.location_point_data_get(field_id,
                              &tmp_points_coords,
                              &tmp_points_uvw,
                              &tmp_points_dist2,
                              &tmp_points_projected_coords);

  *points_coords           = tmp_points_coords[i_part];
  *points_uvw              = tmp_points_uvw[i_part];
  *points_dist2            = tmp_points_dist2[i_part];
  *points_projected_coords = tmp_points_projected_coords[i_part];
}

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
)
{
  cwipi::Coupling& cpl = _cpl_get(local_code_name,cpl_id);

  CWP_Spatial_interp_t spatial_interp_algorithm = cpl.spatialInterpAlgoGet();

  if ((spatial_interp_algorithm != CWP_SPATIAL_INTERP_FROM_LOCATION_MESH_LOCATION_LOCATE_ALL_TGT) &&
      (spatial_interp_algorithm != CWP_SPATIAL_INTERP_FROM_LOCATION_MESH_LOCATION_OCTREE)         &&
      (spatial_interp_algorithm != CWP_SPATIAL_INTERP_FROM_LOCATION_MESH_LOCATION_BOXTREE)) {
    PDM_error(__FILE__, __LINE__, 0, "Getter unavailable for spatial interpolation algorithm %d\n", spatial_interp_algorithm);
  }

  int **tmp_cell_vtx_idx = NULL;
  int **tmp_cell_vtx     = NULL;
  cpl.location_internal_cell_vtx_get(field_id,
                                     &tmp_cell_vtx_idx,
                                     &tmp_cell_vtx);

  *cell_vtx_idx = tmp_cell_vtx_idx[i_part];
  *cell_vtx     = tmp_cell_vtx[i_part];
}

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
)
{
  cwipi::Coupling& cpl = _cpl_get(local_code_name,cpl_id);

  CWP_Spatial_interp_t spatial_interp_algorithm = cpl.spatialInterpAlgoGet();

  if (spatial_interp_algorithm != CWP_SPATIAL_INTERP_FROM_INTERSECTION) {
    PDM_error(__FILE__, __LINE__, 0, "Getter unavailable for spatial interpolation algorithm %d\n", spatial_interp_algorithm);
  }

  int     **tmp_volumes_idx = NULL; // UNUSED
  double  **tmp_volumes     = NULL;
  cpl.weight_get(field_id,
                 &tmp_volumes_idx,
                 &tmp_volumes);

  *volumes     = tmp_volumes[i_part];
}


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
)
{
  cwipi::Coupling& cpl = _cpl_get(local_code_name,cpl_id);

  CWP_Spatial_interp_t spatial_interp_algorithm = cpl.spatialInterpAlgoGet();

  if (spatial_interp_algorithm != CWP_SPATIAL_INTERP_FROM_INTERSECTION) {
    PDM_error(__FILE__, __LINE__, 0, "Getter unavailable for spatial interpolation algorithm %d\n", spatial_interp_algorithm);
  }

  double **_tgt_elt_volumes = NULL;
  cpl.intersection_tgt_elt_volumes_get(field_id,
                                       &_tgt_elt_volumes);

  *tgt_elt_volumes = _tgt_elt_volumes[i_part];
}

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
)
{
  cwipi::Coupling& cpl = _cpl_get(local_code_name,cpl_id);

  CWP_Spatial_interp_t spatial_interp_algorithm = cpl.spatialInterpAlgoGet();

  if (spatial_interp_algorithm != CWP_SPATIAL_INTERP_FROM_NEAREST_SOURCES_LEAST_SQUARES &&
      spatial_interp_algorithm != CWP_SPATIAL_INTERP_FROM_NEAREST_TARGETS_LEAST_SQUARES) {
    PDM_error(__FILE__, __LINE__, 0, "Getter unavailable for spatial interpolation algorithm %d\n", spatial_interp_algorithm);
  }

  int    **tmp_distances2_idx = NULL; // UNUSED
  double **tmp_distances2     = NULL;
  cpl.weight_get(field_id,
                 &tmp_distances2_idx,
                 &tmp_distances2);

  *distances2     = tmp_distances2[i_part];
}


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
)
{
  cwipi::Coupling& cpl = _cpl_get(local_code_name,cpl_id);

  CWP_Spatial_interp_t spatial_interp_algorithm = cpl.spatialInterpAlgoGet();

  if (spatial_interp_algorithm != CWP_SPATIAL_INTERP_FROM_NEAREST_SOURCES_LEAST_SQUARES &&
      spatial_interp_algorithm != CWP_SPATIAL_INTERP_FROM_NEAREST_TARGETS_LEAST_SQUARES) {
    PDM_error(__FILE__, __LINE__, 0, "Getter unavailable for spatial interpolation algorithm %d\n", spatial_interp_algorithm);
  }

  double **_nearest_src_coord = NULL;
  cpl.closest_point_src_coord_get(field_id,
                                  &_nearest_src_coord);

  *nearest_src_coord = nearest_src_coord[i_part];
}


/**
 * \brief Delete a field.
 *
 * \param [in] local_code_name Local code name
 * \param [in]  cpl_id         Coupling identifier
 * \param [in]  field_id       Field identifier
 *
 */

void
CWP_Field_del
(
 const char                  *local_code_name,
 const char                  *cpl_id,
 const char                  *field_id
)
{
  cwipi::Coupling& cpl = _cpl_get(local_code_name,cpl_id);
  cpl.fieldDel(field_id);
}

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
 * \param [in]  local_code_name  Local code name
 * \param [in]  cpl_id           Coupling identifier
 * \param [in]  field_id        Field identifier
 *
 *
 */

void
CWP_Field_issend
(
 const char        *local_code_name,
 const char        *cpl_id,
 const char        *field_id
)
{
  cwipi::Coupling& cpl = _cpl_get(local_code_name,cpl_id);
  std::string referenceFieldID_str = const_cast<char*>(field_id);

  cpl.issend(referenceFieldID_str);
}


/**
 *
 * \brief Receive a spatially interpolated field from the coupled code
 *        with non-blocking communications.
 *
 * This function is independent of \ref CWP_Time_exch_t mode. The user has to
 * manually check the consistency of the exchanges.
 *
 * \param [in] local_code_name  Local code name
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
)
{
  cwipi::Coupling& cpl = _cpl_get(local_code_name,cpl_id);
  std::string targetFieldID_str = const_cast<char*>(tgt_field_id);
  cpl.irecv(targetFieldID_str);
}


/**
 *
 * \brief Wait the end of an exchange related to request from \ref CWP_Field_issend.
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
)
{
  cwipi::Coupling& cpl = _cpl_get(local_code_name,cpl_id);
  std::string srcFieldID_str = const_cast<char*>(field_id);
  cpl.waitIssend(srcFieldID_str);
}


/**
 *
 * \brief Wait the end of an exchange related to request from \ref CWP_Field_irecv.
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
CWP_Field_wait_irecv
(
 const char  *local_code_name,
 const char  *cpl_id,
 const char  *tgt_field_id
)
{
  cwipi::Coupling& cpl = _cpl_get(local_code_name,cpl_id);
  std::string distantFieldID_str = const_cast<char*>(tgt_field_id);
  cpl.waitIrecv(distantFieldID_str);
}

/*----------------------------------------------------------------------------*
 * Functions about current application control parameters                     *
 *----------------------------------------------------------------------------*/

/**
 *
 * \brief Add a new parameter and initialize it.
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
)
{
  cwipi::CodePropertiesDB & properties =
    cwipi::CodePropertiesDB::getInstance();

  const string &codeNameStr = local_code_name;
  const string &nameStr = param_name;

  switch(data_type) {
  case CWP_INT :
    properties.ctrlParamAdd<int>(codeNameStr,nameStr, *(int *)initial_value);
    break;
  case CWP_DOUBLE :
    properties.ctrlParamAdd<double>(codeNameStr,nameStr, *(double *)initial_value);
    break;
  case CWP_CHAR :
    properties.ctrlParamAdd<char *>(codeNameStr,nameStr, *(char **)initial_value);
    break;
  default :
    PDM_error(__FILE__, __LINE__, 0,
               "Not yet implemented data type\n");
  }
}



/**
 *
 * \brief Set a parameter.
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
)
{
  cwipi::CodePropertiesDB & properties =
    cwipi::CodePropertiesDB::getInstance();

  const string &codeNameStr = local_code_name;
  const string &nameStr = param_name;
  switch(data_type) {
  case CWP_INT :
    properties.ctrlParamSet<int>(codeNameStr, nameStr, * (int *) value);
    break;
  case CWP_DOUBLE :
    properties.ctrlParamSet<double>(codeNameStr, nameStr, * (double *) value);
    break;
  case CWP_CHAR :
    properties.ctrlParamSet<char *>(codeNameStr, nameStr, * (char **) value);
    break;
  default :
    PDM_error(__FILE__, __LINE__, 0,
               "Not yet implemented data type\n");
  }
}



/**
 *
 * \brief Delete a parameter.
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
)
{
  cwipi::CodePropertiesDB & properties =
    cwipi::CodePropertiesDB::getInstance();

  const string &codeNameStr = local_code_name;
  const string &nameStr = param_name;
  switch(data_type) {
  case CWP_INT :
    properties.ctrlParamCancel<int>(codeNameStr, nameStr);
    break;
  case CWP_DOUBLE :
    properties.ctrlParamCancel<double>(codeNameStr, nameStr);
    break;
  case CWP_CHAR :
    properties.ctrlParamCancel<string>(codeNameStr, nameStr);
    break;
  default :
    PDM_error(__FILE__, __LINE__, 0,
               "Not yet implemented data type\n");
  }
}

/*----------------------------------------------------------------------------*
 * Functions about other application control parameters                       *
 *----------------------------------------------------------------------------*/


/**
 *
 * \brief Return the number of parameters for the code \p code_name.
 *
 * \param [in] code_name       Local or distant code name
 * \param [in] data_type       Parameter type,
 *
 * return  Number of parameters
 *
 */

int
CWP_Param_n_get
(
 const char             *code_name,
 const CWP_Type_t        data_type
)
{
  cwipi::CodePropertiesDB & properties =
    cwipi::CodePropertiesDB::getInstance();

  const string &codeNameStr = code_name;
  int nParam;
  switch(data_type) {
  case CWP_INT :
    nParam = properties.ctrlParamNGet<int>(codeNameStr);
    break;
  case CWP_DOUBLE :
    nParam = properties.ctrlParamNGet<double>(codeNameStr);
    break;
  case CWP_CHAR :
    nParam = properties.ctrlParamNGet<string>(codeNameStr);
    break;
  default :
    PDM_error(__FILE__, __LINE__, 0,
                "Not yet implemented data type\n");
  }
  return nParam;
}



/**
 *
 * \brief Return the list of parameters for the code \p code_name.
 *
 * \param [in]  code_name      Local or distant code name
 * \param [in]  data_type      Parameter type,
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
)
{
  cwipi::CodePropertiesDB & properties =
    cwipi::CodePropertiesDB::getInstance();

  const string &codeNameStr = code_name;

  switch(data_type) {
  case CWP_INT :
    properties.ctrlParamListGet<int>(codeNameStr, nParam, paramNames);
    break;
  case CWP_DOUBLE :
    properties.ctrlParamListGet<double>(codeNameStr, nParam, paramNames);
    break;
  case CWP_CHAR :
    properties.ctrlParamListGet<string>(codeNameStr, nParam, paramNames);
    break;
  default :
    PDM_error(__FILE__, __LINE__, 0,
               "Not yet implemented data type\n");
  }
}


/**
 *
 * \brief Is this \p code_name a parameter ?
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
)
{
  cwipi::CodePropertiesDB & properties =
    cwipi::CodePropertiesDB::getInstance();

  const string &codeNameStr = code_name;
  const string &nameStr     = param_name;
  int isParam;
  switch(data_type) {
  case CWP_INT :
    isParam = properties.ctrlParamIs<int>(codeNameStr, nameStr);
    break;
  case CWP_DOUBLE :
    isParam = properties.ctrlParamIs<double>(codeNameStr, nameStr);
    break;
  case CWP_CHAR :
    isParam = properties.ctrlParamIs<string>(codeNameStr, nameStr);
    break;
  default :
    PDM_error(__FILE__, __LINE__, 0,
               "Not yet implemented data type\n");
  }
  return isParam;
}


/**
 *
 * \brief Return the parameter value of \p param_name on \p code_name.
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
)
{
  cwipi::CodePropertiesDB & properties =
    cwipi::CodePropertiesDB::getInstance();

  const string &codeNameStr = code_name;
  const string &nameStr     = param_name;
  switch(data_type) {
  case CWP_INT : {
    int *intValue = (int *) value;
    *intValue = properties.ctrlParamGet<int>(codeNameStr, nameStr);
    break;
  }
  case CWP_DOUBLE : {
    double *dblValue = (double *) value;
    *dblValue = properties.ctrlParamGet<double>(codeNameStr, nameStr);
    break;
  }
  case CWP_CHAR : {
    char **charValue = (char **) value;
    *charValue = properties.ctrlParamGet<char *>(codeNameStr, nameStr);
    break;
  }
  default : {
    PDM_error(__FILE__, __LINE__, 0,
              "Not yet implemented data type\n");
  }
  }
}


/**
 *
 * \brief Return the result of a reduce operation about a parameter.
 *
 * The parameter name has to be the same for all codes.
 *
 * \param [in]  op           Operation
 * \param [in]  param_name   Parameter name
 * \param [in]  data_type    Parameter type,
 * \param [out] res          Result
 * \param [in]  nCode        Number of codes
 * \param [in]  code_names   Codes name
 *
 */

void
CWP_Param_reduce
(
 const CWP_Op_t    op,
 const char       *param_name,
 const CWP_Type_t  data_type,
 void             *res,
 const int         nCode,
 const char      **code_names
)
{
  cwipi::CodePropertiesDB & properties =
    cwipi::CodePropertiesDB::getInstance();

  const string &nameStr = param_name;

  switch(data_type) {
  case CWP_INT : {
    int *intRes = (int *) res;
    properties.ctrlParamReduce<int>(op, nameStr, intRes, nCode, code_names);
    break;
  }
  case CWP_DOUBLE : {
    double *doubleRes = (double *) res;
    properties.ctrlParamReduce<double>(op, nameStr, doubleRes, nCode, code_names);
    break;
  }
  default :
    PDM_error(__FILE__, __LINE__, 0,
               "Not yet implemented data type\n");
  }
}



/**
 *
 * \brief Lock access to local parameters from a distant code.
 *
 * \param [in]  code_name  Code to lock
 *
 */

void
CWP_Param_lock
(
const char *code_name
)
{
  const string &nameStr = code_name;

  cwipi::CodePropertiesDB & properties =
    cwipi::CodePropertiesDB::getInstance();

  properties.lock(nameStr);
}



/**
 *
 * \brief Unlock access to local parameters from a distant code.
 *
 * \param [in]  code_name  Code to unlock
 *
 */


void
CWP_Param_unlock
(
const char *code_name
)
{
  const string &nameStr = code_name;

  cwipi::CodePropertiesDB & properties =
    cwipi::CodePropertiesDB::getInstance();

  properties.unLock(nameStr);
}


/*----------------------------------------------------------------------------*
 * Functions about data exchange                                              *
 *----------------------------------------------------------------------------*/

/**
 * \brief Initiate the sending of a data array.
 *
 * \param [in] local_code_name  Local code name
 * \param [in] cpl_id           Coupling identifier
 * \param [in] global_data_id   Global data identifier
 * \param [in] s_send_data      Data size
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
)
{
  cwipi::Coupling& cpl = _cpl_get(local_code_name,cpl_id);
  cpl.globalDataIssend(global_data_id,
                       s_send_entity,
                       send_stride,
                       n_send_entity,
                       send_data);
}


/**
 * \brief Initiate the reception of a data array.
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
)
{
  assert(recv_data != NULL);

  cwipi::Coupling& cpl = _cpl_get(local_code_name,cpl_id);
  cpl.globalDataIrecv(global_data_id,
                      s_recv_entity,
                      recv_stride,
                      n_recv_entity,
                      recv_data);
}


/**
 * \brief Finalize the sending of a data array.
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
)
{
  cwipi::Coupling& cpl = _cpl_get(local_code_name,cpl_id);
  cpl.globalDataWaitIssend(global_data_id);
}


/**
 * \brief Finalize the reception of a data array.
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
)
{
  cwipi::Coupling& cpl = _cpl_get(local_code_name,cpl_id);
  cpl.globalDataWaitIrecv(global_data_id);
}

/**
 * \brief Create partitioned data exchange object
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
CWP_Part_data_create
(
 const char           *local_code_name,
 const char           *cpl_id,
 const char           *part_data_id,
 CWP_PartData_exch_t   exch_type,
 CWP_g_num_t         **gnum_elt,
 int                  *n_elt,
 int                   n_part
)
{
  cwipi::Coupling& cpl = _cpl_get(local_code_name,cpl_id);
  cpl.partDataCreate(part_data_id,
                     exch_type,
                     gnum_elt,
                     n_elt,
                     n_part);
}


/**
 * \brief Delete partitioned data exchange object
 *
 * \param [in] local_code_name  Local code name
 * \param [in] cpl_id           Coupling identifier
 * \param [in] part_data_id     Part data identifier
 *
 */

void
CWP_Part_data_del
(
 const char          *local_code_name,
 const char          *cpl_id,
 const char          *part_data_id
)
{
  cwipi::Coupling& cpl = _cpl_get(local_code_name,cpl_id);
  cpl.partDataDel(part_data_id);
}


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
)
{
  cwipi::Coupling& cpl = _cpl_get(local_code_name,cpl_id);
  cpl.partDataIssend(part_data_id,
                     exch_id,
                     s_data,
                     n_components,
                     send_data);
}


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
)
{
  cwipi::Coupling& cpl = _cpl_get(local_code_name,cpl_id);
  cpl.partDataIrecv(part_data_id,
                    exch_id,
                    s_data,
                    n_components,
                    recv_data);
}


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
)
{
  cwipi::Coupling& cpl = _cpl_get(local_code_name,cpl_id);
  cpl.partDataWaitIssend(part_data_id,
                         exch_id);
}


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
)
{
  cwipi::Coupling& cpl = _cpl_get(local_code_name,cpl_id);
  cpl.partDataWaitIrecv(part_data_id,
                        exch_id);
}


int
CWP_Part_data_n_part_get
(
 const char    *local_code_name,
 const char    *cpl_id,
 const char    *part_data_id
)
{
  cwipi::Coupling& cpl = _cpl_get(local_code_name,cpl_id);
  return cpl.partDataNPartGet(part_data_id);
}


int
CWP_Part_data_n_ref_get
(
 const char    *local_code_name,
 const char    *cpl_id,
 const char    *part_data_id,
 const int      i_part
)
{
  cwipi::Coupling& cpl = _cpl_get(local_code_name,cpl_id);
  return cpl.partDataNRefGet(part_data_id,
                             i_part);
}


/**
 * \brief Set a generic high order block to the interface mesh. <b>(Not implemented yet)</b>
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
)
{
  if (_is_active_rank(local_code_name)) {
    const char *prefix = "CWP_HO_ORDERING";
    /* Generate an ho_ordering name from local_code_name and cpl_id */
    char *ho_ordering_name = (char *) malloc(sizeof(char) * (strlen(prefix)          +
                                                             strlen(local_code_name) +
                                                             strlen(cpl_id)          +
                                                             3));
    sprintf(ho_ordering_name, "%s_%s_%s", prefix, local_code_name, cpl_id);

    cwipi::Coupling& cpl = _cpl_get(local_code_name,cpl_id);
    cpl.meshHOBlockSet(i_part,
                       block_id,
                       n_elts,
                       connec,
                       global_num,
                       order,
                       ho_ordering_name);

    free(ho_ordering_name);
  }
}


/**
 * \brief Get the properties of a generic high order block of the interface mesh. <b>(Not implemented yet)</b>
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
)
{
  cwipi::Coupling& cpl = _cpl_get(local_code_name,cpl_id);
  cpl.meshHOBlockGet(i_part,
                     block_id,
                     n_elts,
                     order,
                     connec,
                     global_num);
}


/**
 *
 * \brief Define ho element ordering from the location in the (u, v, w) grid
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
CWP_Mesh_interf_ho_ordering_from_IJK_set
(
 const char        *local_code_name,
 const char        *cpl_id,
 const CWP_Block_t  block_type,
 const int          order,
 const int          n_nodes,
 const int         *ijk_grid
 )
{
  const char *prefix = "CWP_HO_ORDERING";
  /* Generate an ho_ordering name from local_code_name and cpl_id */
  char *ho_ordering_name = (char *) malloc(sizeof(char) * (strlen(prefix)          +
                                                           strlen(local_code_name) +
                                                           strlen(cpl_id)          +
                                                           3));
  sprintf(ho_ordering_name, "%s_%s_%s", prefix, local_code_name, cpl_id);

  PDM_Mesh_nodal_elt_t t_elt = (PDM_Mesh_nodal_elt_t) block_type;

  PDM_ho_ordering_user_to_ijk_add(ho_ordering_name,
                                  t_elt,
                                  order,
                                  n_nodes,
                                  ijk_grid);
  free(ho_ordering_name);
}


/**
 *
 * \brief Get the coupling spatial interpolation algorithm.
 *
 * \param [in] local_code_name  Local code name
 * \param [in] cpl_id           Coupling identifier
 *
 * \return                      Spatial interpolation algorithm
 */

CWP_Spatial_interp_t
CWP_Cpl_spatial_interp_algo_get
(
 const char *local_code_name,
 const char *cpl_id
 )
{
  cwipi::Coupling& cpl = _cpl_get(local_code_name, cpl_id);

  return cpl.spatialInterpAlgoGet();
};

/*-----------------------------------------------------------------------------*/

#ifdef __cplusplus
}
#endif /* __cplusplus */
