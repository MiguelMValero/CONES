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
 * Standard C/C++ library headers
 *----------------------------------------------------------------------------*/

#include <mpi.h>

#include <cassert>
#include <cstring>

/*----------------------------------------------------------------------------
 * BFT library headers
 *----------------------------------------------------------------------------*/

#include <bftc_mem.h>
#include <bftc_printf.h>

/*----------------------------------------------------------------------------
 * FVM library headers
 *----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "cwipi_cf.h"
#include "applicationPropertiesDataBase.hxx"
#include "couplingDataBase.hxx"
#include "couplingDataBase_i.hxx"
#include "cwipi.h"
#include "cwipi_config.h"
#include "oldCoupling.hxx"
#include "oldCoupling_i.hxx"

/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force back Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */

/*============================================================================
 * Fortran function header
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Fortran printing
 *----------------------------------------------------------------------------*/

#ifndef CWP_HAVE_NOT_FORTRAN_IN_C
  void PROCF (printfort, PRINTFORT) (char *buf_print_f, int *msgsize);
#endif

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

/*============================================================================
 * Private function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 *
 * Convert a fortran string to a C string
 *
 * parameters:
 *   application_name_f    <-- Fortran string
 *   l_application_name_f  <-- Fortran string length
 *
 * return:
 *   C string
 *
 *----------------------------------------------------------------------------*/

static char *_cwipi_fortran_to_c_string(const char *application_name_f,
                                  const int l_application_name_f)
{
  char *application_name_c = NULL;
  int imin = 0;
  int imax = 0;

  while (imin < l_application_name_f &&
         application_name_f[imin] == ' ')
    imin++;

  while (imax < l_application_name_f &&
         application_name_f[l_application_name_f-imax-1] == ' ')
    imax++;

  imax = l_application_name_f-imax-1;

  assert(imax >= imin);

  if ((imax == l_application_name_f) || (imin == l_application_name_f)) {
    application_name_c =  (char *) malloc (sizeof(char) * (1));
    application_name_c[0] = '\0';
  }
  else {
    int size = imax - imin + 2;
    application_name_c =  (char *) malloc (sizeof(char) * (size));
    int index = 0;
    for (int k = imin; k <= imax; k++)
      application_name_c[index++] = application_name_f[k];
    application_name_c[index] = '\0';
  }

  return application_name_c;
}


/*============================================================================
 * Public function definitions
 *============================================================================*/


/*----------------------------------------------------------------------------
 *
 * Initialize the cwipi library.
 * Redirect outputs in a file (Standard output with output_listing = NULL or
 * output_logical_unit = -1)
 * Create the current communicator application from 'common_comm'.
 *
 * parameters:
 *   common_comm         <-- Common MPI communicator
 *   application_name    <-- Current application name
 *   application_comm    --> Internal MPI communicator for the current
 *                           application
 *
 * This is a synchronization between all applications
 *----------------------------------------------------------------------------*/

void PROCF(cwipi_init_cf, CWIPI_INIT_CF)
  (MPI_Fint  *common_fcomm,
   const char *application_name_f,
   const int  *l_application_name,
   MPI_Fint   *application_fcomm
   ARGF_SUPP_CHAINE)
{
  MPI_Comm common_comm = MPI_Comm_f2c(*common_fcomm);

  MPI_Comm application_comm = MPI_COMM_NULL;

  char *application_name_c = _cwipi_fortran_to_c_string(application_name_f,
                                                            *l_application_name);

  bftc_printf("\ncwipi " CWIPI_VERSION " initializing\n");
  bftc_printf("------------------------\n\n");

  cwipi::ApplicationPropertiesDataBase & properties =
    cwipi::ApplicationPropertiesDataBase::getInstance();

  application_comm = properties.init(application_name_c,
                                     common_comm);

  *application_fcomm = MPI_Comm_c2f(application_comm);


  free ( application_name_c);
}

/*----------------------------------------------------------------------------
 *
 * Add a integer control parameter
 *
 * parameters
 *    name           <-- parameter name
 *    initial_value  <-- initial value
 *
 *----------------------------------------------------------------------------*/

void PROCF(cwipi_add_loc_int_ctrl_param_cf,
           CWIPI_ADD_LOC_INT_CTRL_PARAM_CF)
  (const char *name,
   const int  *l_name,
   int *initial_value
   ARGF_SUPP_CHAINE)
{
  char* nameC = _cwipi_fortran_to_c_string(name, *l_name);
  cwipi_add_local_int_control_parameter(nameC, *initial_value);
  free ( nameC);
}

/*----------------------------------------------------------------------------
 *
 * Add a double control parameter
 *
 * parameters
 *    name           <-- parameter name
 *    initial_value  <-- initial value
 *
 *----------------------------------------------------------------------------*/

void PROCF(cwipi_add_loc_dbl_ctrl_param_cf,
           CWIPI_ADD_LOC_DBL_CTRL_PARAM_CF)
  (const char *name,
   const int  *l_name,
   double *initial_value
   ARGF_SUPP_CHAINE)
{
  char* nameC = _cwipi_fortran_to_c_string(name, *l_name);
  cwipi_add_local_double_control_parameter(nameC, *initial_value);
  free ( nameC);
}

/*----------------------------------------------------------------------------
 *
 * Add a string control parameter
 *
 * parameters
 *    name           <-- parameter name
 *    initial_value  <-- initial value
 *
 *----------------------------------------------------------------------------*/

void PROCF(cwipi_add_loc_str_ctrl_param_cf,
           CWIPI_ADD_LOC_STR_CTRL_PARAM_CF)
  (const char *name,
   const int  *l_name,
   char *initial_value,
   int *l_value
   ARGF_SUPP_CHAINE)
{
  char* nameC = _cwipi_fortran_to_c_string(name, *l_name);
  char* valueC = _cwipi_fortran_to_c_string(initial_value, *l_value);
  cwipi_add_local_string_control_parameter(nameC, valueC);
  free ( nameC);
  free ( valueC);
}

/*----------------------------------------------------------------------------
 *
 * Set a integer control parameter
 *
 * parameters
 *    name           <-- parameter name
 *    value          <-- value
 *
 *----------------------------------------------------------------------------*/

void PROCF(cwipi_set_loc_int_ctrl_param_cf,
           CWIPI_SET_LOC_INT_CTRL_PARAM_CF)
  (const char *name,
   const int  *l_name,
   int *value
   ARGF_SUPP_CHAINE)
{
  char* nameC = _cwipi_fortran_to_c_string(name, *l_name);
  cwipi_set_local_int_control_parameter(nameC, *value);
  free ( nameC);
}

/*----------------------------------------------------------------------------
 *
 * Set a double control parameter
 *
 * parameters
 *    name           <-- parameter name
 *    value          <-- value
 *
 *----------------------------------------------------------------------------*/

void PROCF(cwipi_set_loc_dbl_ctrl_param_cf,
           CWIPI_SET_LOC_DBL_CTRL_PARAM_CF)
  (const char *name,
   const int  *l_name,
   double *value
   ARGF_SUPP_CHAINE)
{
  char* nameC = _cwipi_fortran_to_c_string(name, *l_name);
  cwipi_set_local_double_control_parameter(nameC, *value);
  free ( nameC);
}

/*----------------------------------------------------------------------------
 *
 * Set a string control parameter
 *
 * parameters
 *    name           <-- parameter name
 *    value          <-- value
 *
 *----------------------------------------------------------------------------*/

void PROCF(cwipi_set_loc_str_ctrl_param_cf,
           CWIPI_SET_LOC_STR_CTRL_PARAM_CF)
  (const char *name,
   const int  *l_name,
   char *initial_value,
   int *l_value
   ARGF_SUPP_CHAINE)
{
  char* nameC = _cwipi_fortran_to_c_string(name, *l_name);
  char* valueC = _cwipi_fortran_to_c_string(initial_value, *l_value);
  cwipi_set_local_string_control_parameter(nameC, valueC);
  free ( nameC);
  free ( valueC);
}

/*----------------------------------------------------------------------------
 *
 * Get a integer control parameter of the current application
 *
 * parameters
 *    name           <-- parameter name
 *
 *----------------------------------------------------------------------------*/

void PROCF(cwipi_get_loc_int_ctrl_param_cf,
           CWIPI_GET_LOC_INT_CRTL_PARAM_CF)
  (const char *name,
   const int  *l_name,
   int *value
   ARGF_SUPP_CHAINE)
{
  char* nameC = _cwipi_fortran_to_c_string(name, *l_name);
  *value = cwipi_get_local_int_control_parameter(nameC);
  free ( nameC);
}

/*----------------------------------------------------------------------------
 *
 * Get a double control parameter of the current application
 *
 * parameters
 *    name           <-- parameter name
 *
 *----------------------------------------------------------------------------*/

void PROCF(cwipi_get_loc_dbl_ctrl_param_cf,
           CWIPI_GET_LOC_DBL_CTRL_PARAM_CF)
  (const char *name,
   const int  *l_name,
   double *value
   ARGF_SUPP_CHAINE)
{
  char* nameC = _cwipi_fortran_to_c_string(name, *l_name);
  *value = cwipi_get_local_double_control_parameter(nameC);
  free ( nameC);
}

/*----------------------------------------------------------------------------
 *
 * Get a string control parameter of the current application
 *
 * parameters
 *    name           <-- parameter name
 *
 *----------------------------------------------------------------------------*/

void PROCF(cwipi_get_loc_str_ctrl_param_cf,
           CWIPI_GET_LOC_STR_CTRL_PARAM_CF)
  (const char *name,
   const int  *l_name,
   char *value_str_f,
   int *l_value_str_f
   ARGF_SUPP_CHAINE)
{
  char* nameC = _cwipi_fortran_to_c_string(name, *l_name);
  const char* value_str_c = cwipi_get_local_string_control_parameter(nameC);
  const int l_value_str_c = strlen(value_str_c);

  int i = 0;
  while (i < l_value_str_c && i < *l_value_str_f) {
    value_str_f[i] = value_str_c[i];
    i+=1;
  }

  while (i < *l_value_str_f) {
    value_str_f[i] = ' ';
    i+=1;
  }
  free ( nameC);
}

/*----------------------------------------------------------------------------
 *
 * Delete a current application parameter
 *
 * parameters
 *    name           <-- parameter name
 *
 *----------------------------------------------------------------------------*/

void PROCF(cwipi_del_loc_int_ctrl_param_cf,
           CWIPI_DEL_LOC_INT_CTRL_PARAM_CF)
  (const char *name,
   const int  *l_name
   ARGF_SUPP_CHAINE)
{
  char* nameC = _cwipi_fortran_to_c_string(name, *l_name);
  cwipi_delete_local_int_control_parameter(nameC);
  free ( nameC);
}

/*----------------------------------------------------------------------------
 *
 * Delete a current application parameter
 *
 * parameters
 *    name           <-- parameter name
 *
 *----------------------------------------------------------------------------*/

void PROCF(cwipi_del_loc_dbl_ctrl_param_cf,
           CWIPI_DEL_LOC_DBL_CTRL_PARAM_CF)
  (const char *name,
   const int  *l_name
   ARGF_SUPP_CHAINE)
{
  char* nameC = _cwipi_fortran_to_c_string(name, *l_name);
  cwipi_delete_local_double_control_parameter(nameC);
  free ( nameC);
}

/*----------------------------------------------------------------------------
 *
 * Delete a current application parameter
 *
 * parameters
 *    name           <-- parameter name
 *
 *----------------------------------------------------------------------------*/

void PROCF(cwipi_del_loc_str_ctrl_param_cf,
           CWIPI_DEL_LOC_STR_CTRL_PARAM_CF)
  (const char *name,
   const int  *l_name
   ARGF_SUPP_CHAINE)
{
  char* nameC = _cwipi_fortran_to_c_string(name, *l_name);
  cwipi_delete_local_string_control_parameter(nameC);
  free ( nameC);
}

/*----------------------------------------------------------------------------
 *
 * Get a integer control parameter of a other application
 *
 * parameters
 *    application_name       <-- application name
 *    name                   <-- parameter name
 *
 *----------------------------------------------------------------------------*/

void PROCF(cwipi_get_dis_int_ctrl_param_cf,
           CWIPI_GET_DIS_INT_CTRL_PARAM_CF)
  (const char *application_name,
   const int  *l_application_name,
   const char *name,
   const int  *l_name,
   int *value
   ARGF_SUPP_CHAINE)
{
  char *application_nameC =
    _cwipi_fortran_to_c_string(application_name, *l_application_name);
  char *nameC = _cwipi_fortran_to_c_string(name, *l_name);

  *value = cwipi_get_distant_int_control_parameter(application_nameC,
                                                       nameC);

  free ( nameC);
  free ( application_nameC);
}

/*----------------------------------------------------------------------------
 *
 * Get a double control parameter of a other application
 *
 * parameters
 *    application_name    <-- application name
 *    name                <-- parameter name
 *
 *----------------------------------------------------------------------------*/

void PROCF(cwipi_get_dis_dbl_ctrl_param_cf,
           CWIPI_GET_DIS_DBL_CTRL_PARAM_CF)
  (const char *application_name,
   const int  *l_application_name,
   const char *name,
   const int  *l_name,
   double *value
   ARGF_SUPP_CHAINE)
{
  char *application_nameC =
    _cwipi_fortran_to_c_string(application_name, *l_application_name);
  char *nameC = _cwipi_fortran_to_c_string(name, *l_name);

  *value = cwipi_get_distant_double_control_parameter(application_nameC,
                                                          nameC);

  free ( nameC);
  free ( application_nameC);
}
/*
 * ----------------------------------------------------------------------------
 *
 * Get a string control parameter of a other application
 *
 * parameters
 *    application_name    <-- application name
 *    name                <-- parameter name
 *
 *----------------------------------------------------------------------------*/

void PROCF(cwipi_get_dis_str_ctrl_param_cf,
           CWIPI_GET_DIS_STR_CTRL_PARAM_CF)
  (const char *application_name,
   const int  *l_application_name,
   const char *name,
   const int  *l_name,
   char *value_str_f,
   int  *l_value_str_f
   ARGF_SUPP_CHAINE)

{
  char *application_nameC = _cwipi_fortran_to_c_string(application_name, *l_application_name);
  char *nameC = _cwipi_fortran_to_c_string(name, *l_name);

  const char *value_str_c = cwipi_get_distant_string_control_parameter(application_nameC, nameC);
  const int l_value_str_c = strlen(value_str_c);

  int i = 0;
  while (i < l_value_str_c && i < *l_value_str_f) {
    value_str_f[i] = value_str_c[i];
    i+=1;
  }

  while (i < *l_value_str_f) {
    value_str_f[i] = ' ';
    i+=1;
  }

  free ( nameC);
  free ( application_nameC);
}



/*----------------------------------------------------------------------------
 *
 * Has int control parameter ?
 *
 * parameters
 *    application_name    <-- application name
 *    name                <-- parameter name
 *
 *----------------------------------------------------------------------------*/

void PROCF(cwipi_has_int_ctrl_param_cf,
           CWIPI_HAS_INT_CTRL_PARAM_CF)
  (const char *application_name,
   const int  *l_application_name,
   const char *name,
   const int  *l_name,
   int *status
   ARGF_SUPP_CHAINE)
{
  char *application_nameC = _cwipi_fortran_to_c_string(application_name, *l_application_name);
  char *nameC = _cwipi_fortran_to_c_string(name, *l_name);

  *status = cwipi_has_int_parameter(application_nameC, nameC);
}
/*----------------------------------------------------------------------------
 *
 * Has dbl control parameter ?
 *
 * parameters
 *    application_name    <-- application name
 *    name                <-- parameter name
 *
 *----------------------------------------------------------------------------*/

void PROCF(cwipi_has_dbl_ctrl_param_cf,
           CWIPI_HAS_DBL_CTRL_PARAM_CF)
  (const char *application_name,
   const int  *l_application_name,
   const char *name,
   const int  *l_name,
   int *status
   ARGF_SUPP_CHAINE)
{
  char *application_nameC = _cwipi_fortran_to_c_string(application_name, *l_application_name);
  char *nameC = _cwipi_fortran_to_c_string(name, *l_name);

  *status = cwipi_has_double_parameter(application_nameC, nameC);
}


/*----------------------------------------------------------------------------
 *
 * Has str control parameter ?
 *
 * parameters
 *    application_name    <-- application name
 *    name                <-- parameter name
 *
 *----------------------------------------------------------------------------*/

void PROCF(cwipi_has_str_ctrl_param_cf,
           CWIPI_HAS_STR_CTRL_PARAM_CF)
  (const char *application_name,
   const int  *l_application_name,
   const char *name,
   const int  *l_name,
   int *status
   ARGF_SUPP_CHAINE)
{
  char *application_nameC = _cwipi_fortran_to_c_string(application_name, *l_application_name);
  char *nameC = _cwipi_fortran_to_c_string(name, *l_name);

  *status = cwipi_has_string_parameter(application_nameC, nameC);
}


/*----------------------------------------------------------------------------
 *
 * Get number of int parameters
 *
 * parameters
 *    application_name    <-- application name
 *    name                <-- parameter name
 *
 *----------------------------------------------------------------------------*/

void PROCF(cwipi_get_n_int_ctrl_param_cf,
           CWIPI_GET_N_INT_CTRL_PARAM_CF)
  (const char *application_name,
   const int  *l_application_name,
   int *n_param
   ARGF_SUPP_CHAINE)
{
  char *application_nameC = _cwipi_fortran_to_c_string(application_name, *l_application_name);

  *n_param = cwipi_get_n_int_parameters(application_nameC);
}


/*----------------------------------------------------------------------------
 *
 * Get number of dbl parameters
 *
 * parameters
 *    application_name    <-- application name
 *    name                <-- parameter name
 *
 *----------------------------------------------------------------------------*/

void PROCF(cwipi_get_n_dbl_ctrl_param_cf,
           CWIPI_GET_N_DBL_CTRL_PARAM_CF)
  (const char *application_name,
   const int  *l_application_name,
   int *n_param
   ARGF_SUPP_CHAINE)
{
  char *application_nameC = _cwipi_fortran_to_c_string(application_name, *l_application_name);

  *n_param = cwipi_get_n_double_parameters(application_nameC);
}

/*----------------------------------------------------------------------------
 *
 * Get number of str parameters
 *
 * parameters
 *    application_name    <-- application name
 *    name                <-- parameter name
 *
 *----------------------------------------------------------------------------*/

void PROCF(cwipi_get_n_str_ctrl_param_cf,
           CWIPI_GET_N_STR_CTRL_PARAM_CF)
  (const char *application_name,
   const int  *l_application_name,
   int *n_param
   ARGF_SUPP_CHAINE)
{
  char *application_nameC = _cwipi_fortran_to_c_string(application_name, *l_application_name);

  *n_param = cwipi_get_n_string_parameters(application_nameC);
}


/*----------------------------------------------------------------------------
 *
 * Get int parameters list
 *
 * parameters
 *    application_name    <-- application name
 *    name                <-- parameter name
 *
 *----------------------------------------------------------------------------*/

void PROCF(cwipi_get_list_int_ctrl_param_cf,
           CWIPI_GET_LIST_INT_CTRL_PARAM_CF)
  (const char *application_name,
   const int  *l_application_name,
   char *params,
   const int  *l_param  
   ARGF_SUPP_CHAINE)
{
  char *application_nameC = _cwipi_fortran_to_c_string(application_name, *l_application_name);
  
  int n_param = cwipi_get_n_int_parameters(application_nameC);
  char ** c_params =  cwipi_get_list_int_parameters(application_nameC);
  
  for (int i = 0; i < n_param; i++) {
    for (int j = 0; strlen(c_params[i]); j++)
      params[i * (*l_param) + j] =  c_params[i][j];
    for (int j = strlen(c_params[i]); j < *l_param; j++)
      params[i * (*l_param) + j] = ' ';
    free(c_params[i]);
  }
  free(c_params);
}

/*----------------------------------------------------------------------------
 *
 * Get dbl parameters list
 *
 * parameters
 *    application_name    <-- application name
 *    name                <-- parameter name
 *
 *----------------------------------------------------------------------------*/

void PROCF(cwipi_get_list_dbl_ctrl_param_cf,
           CWIPI_GET_LIST_DBL_CTRL_PARAM_CF)
  (const char *application_name,
   const int  *l_application_name,
   char *params,
   const int  *l_param  
   ARGF_SUPP_CHAINE)
{
  char *application_nameC = _cwipi_fortran_to_c_string(application_name, *l_application_name);
  
  int n_param = cwipi_get_n_double_parameters(application_nameC);
  char ** c_params =  cwipi_get_list_double_parameters(application_nameC);
  
  for (int i = 0; i < n_param; i++) {
    for (int j = 0; strlen(c_params[i]); j++)
      params[i * (*l_param) + j] =  c_params[i][j];
    for (int j = strlen(c_params[i]); j < *l_param; j++)
      params[i * (*l_param) + j] = ' ';
    free(c_params[i]);
  }
  free(c_params);
}

/*----------------------------------------------------------------------------
 *
 * Get str parameters list
 *
 * parameters
 *    application_name    <-- application name
 *    name                <-- parameter name
 *
 *----------------------------------------------------------------------------*/

void PROCF(cwipi_get_list_str_ctrl_param_cf,
           CWIPI_GET_LIST_STR_CTRL_PARAM_CF)
  (const char *application_name,
   const int  *l_application_name,
   char *params,
   const int  *l_param  
   ARGF_SUPP_CHAINE)
{
  char *application_nameC = _cwipi_fortran_to_c_string(application_name, *l_application_name);
  
  int n_param = cwipi_get_n_string_parameters(application_nameC);
  char ** c_params =  cwipi_get_list_string_parameters(application_nameC);
  
  for (int i = 0; i < n_param; i++) {
    for (int j = 0; strlen(c_params[i]); j++)
      params[i * (*l_param) + j] =  c_params[i][j];
    for (int j = strlen(c_params[i]); j < *l_param; j++)
      params[i * (*l_param) + j] = ' ';
    free(c_params[i]);
  }
  free(c_params);
}

/*----------------------------------------------------------------------------
 *
 * Synchronize local control parameters with an other application.
 *  This is a synchronization point with this second application
 *
 * parameters
 *    application_name    <-- application name
 *    name                <22-- parameter name
 *
 *----------------------------------------------------------------------------*/

void PROCF(cwipi_synch_ctrl_param_cf,
           CWIPI_SYNCH_CTRL_PARAM_CF)
  (const char *application_name,
   const int  *l_application_name
   ARGF_SUPP_CHAINE)
{
  char *application_nameC =
    _cwipi_fortran_to_c_string(application_name, *l_application_name);

  cwipi_synchronize_control_parameter(application_nameC);

  free ( application_nameC);
}

/*----------------------------------------------------------------------------
 *
 * Create a coupling object
 *
 * parameters:
 *   couplingName            <-- Coupling identifier
 *   couplingType            <-- Coupling type
 *   cplAppli                <-- Coupled application name
 *   entitiesDim             <-- Mesh entities dimension (1, 2 or 3)
 *   tolerance               <-- Geometric tolerance to locate
 *   meshT                   <-- CWIPI_STATIC_MESH
 *                               CWIPI_MOBILE_MESH (not implemented yet)
 *                               CWIPI_CYCLIC_MESH
 *   solverT                 <-- CWIPI_SOLVER_CELL_CENTER
 *                               CWIPI_SOLVER_CELL_VERTEX
 *   outputFreq              <-- Output frequency
 *   outputFmt               <-- Output format to visualize exchanged fields
 *                               on the coupled mesh. Choice between :
 *                                 - "EnSight Gold"
 *                                 - "MED_fichier"
 *                                 - "CGNS"
 *   outputFmtOpt            <-- Output options
 *                             text                output text files
 *                             binary              output binary files (default)
 *                             big_endian          force binary files
 *                                                 to big-endian
 *                             discard_polygons    do not output polygons
 *                                                 or related values
 *                             discard_polyhedra   do not output polyhedra
 *                                                 or related values
 *                             divide_polygons     tesselate polygons
 *                                                 with triangles
 *                             divide_polyhedra    tesselate polyhedra
 *                                                 with tetrahedra and pyramids
 *                                                 (adding a vertex near
 *                                                 each polyhedron's center)
 *
 *----------------------------------------------------------------------------*/

void PROCF(cwipi_create_coupling_cf,
           CWIPI_CREATE_COUPLING_CF)
( const char *coupling_name,
  const int  *l_coupling_name,
  const int  *coupling_type,
  const char *coupled_application,
  const int  *l_coupled_application,
  const int  *entities_dim,
  const double *tolerance,
  const int *mesh_type,
  const int *solver_type,
  const int  * output_frequency,
  const char  *output_format,
  const int  *l_output_format,
  const char  *output_format_option,
  const int  *l_output_format_option,
  const int *nbLocations
  ARGF_SUPP_CHAINE)

{
  char *coupling_nameC =
    _cwipi_fortran_to_c_string(coupling_name, *l_coupling_name);

  char *coupled_applicationC =
    _cwipi_fortran_to_c_string(coupled_application, *l_coupled_application);

  char *output_formatC =
    _cwipi_fortran_to_c_string(output_format, *l_output_format);

  char *output_format_optionC =
    _cwipi_fortran_to_c_string(output_format_option, *l_output_format_option);

  cwipi_create_coupling(coupling_nameC,
                        (cwipi_coupling_type_t) *coupling_type,
                        coupled_applicationC,
                        *entities_dim,
                        *tolerance,
                        (cwipi_mesh_type_t) *mesh_type,
                        (cwipi_solver_type_t) *solver_type,
                        *output_frequency,
                        output_formatC,
                        output_format_optionC,
                        *nbLocations);

  free ( coupling_nameC);
  free ( coupled_applicationC);
  free ( output_formatC);
  free ( output_format_optionC);
}

/*----------------------------------------------------------------------------
 *
 * Set points to locate. This function must be called if the points to locate
 * do not correspond to :
 *        - vertices for NATURE_NODE nature
 *        - cell center for NATURE_CELL_CENTER nature
 *
 * parameters:
 *   coupling_id        <-- coupling identifier
 *   n_points           <-- number of points to locate
 *   coordinates        <-- coordinates of points to locate (enterlaced)
 *
 *----------------------------------------------------------------------------*/

void PROCF(cwipi_set_points_to_locate_cf,
           CWIPI_SET_POINTS_TO_LOCATE_CF)
  (const char   *coupling_name,
   const int  *l_coupling_name,
   const int    *n_points,
   double *coordinate
   ARGF_SUPP_CHAINE)
{
  char *coupling_nameC =
    _cwipi_fortran_to_c_string(coupling_name, *l_coupling_name);

  cwipi_set_points_to_locate(coupling_nameC,
                                 *n_points,
                                 coordinate);
  free ( coupling_nameC);
}


/*----------------------------------------------------------------------------
 *
 * Define the support mesh for a coupling. The connectivity is ordered if
 * necessary. The connectivity order is :
 *        - 1D : edges
 *        - 2D : triangles, quadrangles, polygons
 *        - 3D : tetrahedra, pyramids, prism, hexaedra, polyhedra
 *
 * parameters:
 *   coupling_id        <-- coupling identifier
 *   dim                <-- space dimension (1, 2 or 3)
 *   n_vertex           <-- number of vertex
 *   n_elements         <-- number of elements
 *   coordinates        <-- vertex enterlaced coordinates
 *   parent_vertex_num  <-- pointer to parent vertex numbers (or NULL)
 *   connectivity_index <-> polygon face -> vertices index (O to n-1)
 *                          size: n_elements + 1
 *   connectivity       <-> element -> vertex connectivity
 *                          size: connectivity_index[n_elements]
 *
 *----------------------------------------------------------------------------*/

void PROCF(cwipi_define_mesh_cf,
           CWIPI_DEFINE_MESH_CF)
  (const char *coupling_name,
   const int  *l_coupling_name,
   const int *n_vertex,
   const int *n_element,
   double *coordinates,
   int *connectivity_index,
   int *connectivity
   ARGF_SUPP_CHAINE)
{
  char *coupling_nameC =
    _cwipi_fortran_to_c_string(coupling_name, *l_coupling_name);

  cwipi_define_mesh(coupling_nameC,
                        *n_vertex,
                        *n_element,
                        coordinates,
                        connectivity_index,
                        connectivity);
  free ( coupling_nameC);
}

void PROCF(cwipi_ho_define_mesh_cf,
           CWIPI_HO_DEFINE_MESH_CF)
  (const char *coupling_name,
   const int  *l_coupling_name,
   const int *n_vertex,
   const int *n_element,
   const int *order, 
   double *coordinates,
   int *connectivity_index,
   int *connectivity
   ARGF_SUPP_CHAINE)
{
  char *coupling_nameC =
    _cwipi_fortran_to_c_string(coupling_name, *l_coupling_name);

  cwipi_ho_define_mesh(coupling_nameC,
                       *n_vertex,
                       *n_element,
                       *order,
                       coordinates,
                       connectivity_index,
                       connectivity);
  free ( coupling_nameC);
}

void PROCF(cwipi_ho_options_set_cf,
           CWIPI_HO_OPTIONS_SET_CF)
  (const char *coupling_name,
   const int  *l_coupling_name,
   const char *option,
   const int  *l_option,
   const char *value,
   const int  *l_value
   ARGF_SUPP_CHAINE)
{
  char *coupling_nameC =
    _cwipi_fortran_to_c_string(coupling_name, *l_coupling_name);
  char *optionC =
    _cwipi_fortran_to_c_string(option, *l_option);
  char *valueC =
    _cwipi_fortran_to_c_string(value, *l_value);

  cwipi_ho_options_set(coupling_nameC,
                       optionC,
                       valueC);

  free ( coupling_nameC);
  free ( optionC);
  free ( valueC);
}

void PROCF(cwipi_add_polyhedra_cf,
           CWIPI_ADD_POLYHEDRA_CF)
  (const char *coupling_name,
   const int  *l_coupling_name,
   const int *n_element,
   int *face_index,
   int *cell_to_face_connectivity,
   const int *n_faces,
   int *face_connectivity_index,
   int *face_connectivity
   ARGF_SUPP_CHAINE)

{
  char *coupling_nameC =
    _cwipi_fortran_to_c_string(coupling_name, *l_coupling_name);

  cwipi_add_polyhedra(coupling_nameC,
                      *n_element,
                      face_index,
                      cell_to_face_connectivity,
                      *n_faces,
                      face_connectivity_index,
                      face_connectivity);
  free ( coupling_nameC);
}

/*----------------------------------------------------------------------------
 *
 * Define ho element ordering from the location in the (u, v, w) grid
 *
 * parameters:
 *   coupling_id     <-- coupling name
 *   t_elt           <-- element type
 *   n_nodes         <-- number of nodes
 *   uvw_grid        <-- user ordering to (u, v, w) grid (size = elt_dim * n_nodes)
 *
 *----------------------------------------------------------------------------*/

void PROCF (cwipi_ho_ordering_from_ijk_set_cf,
            CWIPI_HO_ORDERING_FROM_IJK_SET_CF)
(const char *coupling_id,
 const int *l_coupling_id,
 const int *t_elt,
 const int *n_nodes,
 const int *ijk
 ARGF_SUPP_CHAINE)
{
  char *coupling_idC =
    _cwipi_fortran_to_c_string(coupling_id, *l_coupling_id);

  cwipi_ho_ordering_from_IJK_set (coupling_idC,
                                 (cwipi_element_t) *t_elt,
                                 *n_nodes,
                                 ijk);

  free ( coupling_idC);
}

/*----------------------------------------------------------------------------
 *
 * Define ho element ordering from reference element (definition between 0 - 1)
 *
 *   coupling_id        <-- coupling name
 *   t_elt              <-- element type
 *   n_nodes            <-- number of nodes
 *   coords             <-- node coordinates of reference element
 *                                TODO: decrire ici les elements de reference
 *
 *----------------------------------------------------------------------------*/

void PROCF (cwipi_ho_ordering_from_ref_elt_set_cf,
            CWIPI_HO_ORDERING_FROM_REF_ELT_SET_CF)
(const char   *coupling_id,
 const int *l_coupling_id,
 const int *t_elt,
 const int *n_nodes,
 const double *coords
 ARGF_SUPP_CHAINE)
{
  char *coupling_idC =
    _cwipi_fortran_to_c_string(coupling_id, *l_coupling_id);

  cwipi_ho_ordering_from_ref_elt_set (coupling_idC,
                                     (cwipi_element_t) *t_elt,
                                     *n_nodes,
                                     coords);

  free ( coupling_idC);

}

/*----------------------------------------------------------------------------
 *
 * Location completion.
 * This is a synchronization point with the coupled application
 *
 * parameters
 *   coupling_name        <-- Coupling identifier
 *----------------------------------------------------------------------------*/

void PROCF(cwipi_locate_cf, CWIPI_LOCATE_CF) (const char *coupling_name,
      const int  *l_coupling_name
      ARGF_SUPP_CHAINE)
{

  char *coupling_nameC =
    _cwipi_fortran_to_c_string(coupling_name, *l_coupling_name);

  cwipi::CouplingDataBase & couplingDataBase =
    cwipi::CouplingDataBase::getInstance();

  const std::string &coupling_name_str = coupling_nameC;

  cwipi::oldCoupling& coupling = couplingDataBase.getCoupling(coupling_name_str);

  coupling.locate();

  free ( coupling_nameC);
}

/*----------------------------------------------------------------------------
 *
 * cwipi_set_location_index
 *
 * parameters
 *   coupling_id          <-- Coupling identifier
 *   index                <-- location index
 *----------------------------------------------------------------------------*/
void PROCF(cwipi_set_location_index_cf,
           CWIPI_SET_LOCATION_INDEX_CF) (const char *coupling_name,
                                         const int  *l_coupling_name,
                                         const int  *index
                                         ARGF_SUPP_CHAINE)
{
  char *coupling_nameC =
    _cwipi_fortran_to_c_string(coupling_name, *l_coupling_name);
  cwipi::CouplingDataBase & couplingDataBase =
    cwipi::CouplingDataBase::getInstance();

  const std::string &coupling_name_str = coupling_nameC;

  cwipi::oldCoupling& coupling = couplingDataBase.getCoupling(coupling_name_str);

  coupling.setLocationIndex(*index);

  free ( coupling_nameC);


}

/*----------------------------------------------------------------------------
 *
 * save/load  location 
 *
 * parameters:
 *   coupling_name           <-- Coupling identifier
 *----------------------------------------------------------------------------*/

void PROCF(cwipi_load_location_cf, 
           CWIPI_LOAD_LOCATION_CF)
(const char *coupling_name,
 const int  *l_coupling_name
 ARGF_SUPP_CHAINE)
{
  char *coupling_nameC =
    _cwipi_fortran_to_c_string(coupling_name, *l_coupling_name);
  cwipi::CouplingDataBase & couplingDataBase =
    cwipi::CouplingDataBase::getInstance();

  const std::string &coupling_name_str = coupling_nameC;

  cwipi::oldCoupling& coupling = couplingDataBase.getCoupling(coupling_name_str);

  coupling.loadLocation();

  free ( coupling_nameC);


}

void PROCF(cwipi_save_location_cf,
           CWIPI_SAVE_LOCATION_CF)
(const char *coupling_name, 
 const int  *l_coupling_name
 ARGF_SUPP_CHAINE)
{
  char *coupling_nameC =
    _cwipi_fortran_to_c_string(coupling_name, *l_coupling_name);
  cwipi::CouplingDataBase & couplingDataBase =
    cwipi::CouplingDataBase::getInstance();

  const std::string &coupling_name_str = coupling_nameC;

  cwipi::oldCoupling& coupling = couplingDataBase.getCoupling(coupling_name_str);

  coupling.saveLocation();

  free ( coupling_nameC);


}


/*----------------------------------------------------------------------------
 *
 * cwipi_open_location_file
 *
 * parameters
 *   coupling_id          <-- Coupling identifier
 *   filename             <-- file name 
 *   mode                 <-- "r" : read
 *                            "w" : write
 *----------------------------------------------------------------------------*/

void PROCF(cwipi_open_location_file_cf,
           CWIPI_OPEN_LOCATION_FILE_CF)
 (const char *coupling_name,
  const int  *l_coupling_name,
  char *filename,
  const int  *l_filename,
  const char *mode,
  const int  *l_mode
 ARGF_SUPP_CHAINE)
{
  char *coupling_nameC =
    _cwipi_fortran_to_c_string(coupling_name, *l_coupling_name);
  char *filenameC =
    _cwipi_fortran_to_c_string(filename, *l_filename);
  char *modeC =
    _cwipi_fortran_to_c_string(mode, *l_mode);

  cwipi::CouplingDataBase & couplingDataBase =
    cwipi::CouplingDataBase::getInstance();

  const std::string &coupling_name_str = coupling_nameC;

  cwipi::oldCoupling& coupling = couplingDataBase.getCoupling(coupling_name_str);

  coupling.openLocationFile(filenameC, modeC);

  free ( coupling_nameC);
  free ( filenameC);
  free ( modeC);
}

/*----------------------------------------------------------------------------
 *
 * cwipi_close_location_file
 *
 * parameters
 *   coupling_id          <-- Coupling identifier
 *----------------------------------------------------------------------------*/

 void PROCF(cwipi_close_location_file_cf,
            CWIPI_CLOSE_LOCATION_FILE_CF)
(const char *coupling_name,
  const int  *l_coupling_name
  ARGF_SUPP_CHAINE)
{
  char *coupling_nameC =
    _cwipi_fortran_to_c_string(coupling_name, *l_coupling_name);
  cwipi::CouplingDataBase & couplingDataBase =
    cwipi::CouplingDataBase::getInstance();

  const std::string &coupling_name_str = coupling_nameC;

  cwipi::oldCoupling& coupling = couplingDataBase.getCoupling(coupling_name_str);

  coupling.closeLocationFile();

  free ( coupling_nameC);


}


/*----------------------------------------------------------------------------
 *
 * Return located points distance to th interface
 *
 * parameters
 *   coupling_name        <-- Coupling identifier
 *   distance             <-- distance              
 *----------------------------------------------------------------------------*/

void PROCF(cwipi_dist_located_pts_get_cf, CWIPI_DIST_LOCATED_PTS_GET_CF) 
     (const char *coupling_name,
      const int  *l_coupling_name,
      float      *distance
      ARGF_SUPP_CHAINE)
{
  char *coupling_nameC =
    _cwipi_fortran_to_c_string(coupling_name, *l_coupling_name);

  cwipi::CouplingDataBase & couplingDataBase =
    cwipi::CouplingDataBase::getInstance();

  const std::string &coupling_name_str = coupling_nameC;

  cwipi::oldCoupling& coupling = couplingDataBase.getCoupling(coupling_name_str);
  const int npts = coupling.getNLocatedPoint();

  const float *dist_c = coupling.distance();

  for (int i = 0; i < npts; i++)
    distance[i] = dist_c[i];

  free ( coupling_nameC);
}

/*----------------------------------------------------------------------------
 *
 * Update Location
 *
 * parameters
 *   coupling_name        <-- Coupling identifier
 *----------------------------------------------------------------------------*/

void PROCF(cwipi_update_location_cf, CWIPI_UPDATE_LOCATION_CF) (const char *coupling_name,
                                                                const int  *l_coupling_name
                                                                ARGF_SUPP_CHAINE)
{
  char *coupling_nameC =
    _cwipi_fortran_to_c_string(coupling_name, *l_coupling_name);

  cwipi::CouplingDataBase & couplingDataBase =
    cwipi::CouplingDataBase::getInstance();

  const std::string &coupling_name_str = coupling_nameC;

  cwipi::oldCoupling& coupling = couplingDataBase.getCoupling(coupling_name_str);

  coupling.updateLocation();

  free ( coupling_nameC);
}

/*----------------------------------------------------------------------------
 *
 * Get distant points location
 *
 * parameters
 *   coupling_name        <-- Coupling identifier
 *   location             --> Distant points location
 *----------------------------------------------------------------------------*/

void PROCF(cwipi_get_distant_location_cf,
           CWIPI_GET_DISTANT_LOCATION_CF) (const char *coupling_name,
                                       const int  *l_coupling_name,
                                       int *location
                                       ARGF_SUPP_CHAINE)
{
  char *coupling_nameC =
    _cwipi_fortran_to_c_string(coupling_name, *l_coupling_name);

  cwipi::CouplingDataBase & couplingDataBase =
    cwipi::CouplingDataBase::getInstance();

  const std::string &coupling_name_str = coupling_nameC;

  cwipi::oldCoupling& coupling = couplingDataBase.getCoupling(coupling_name_str);

  const unsigned int *locationC = (const unsigned int *) coupling.getDistantLocation();
  const int nDistantPoint = coupling.getNDistantPoint();
  for (int i = 0; i < nDistantPoint; i++)
    location[i] = locationC[i];

  free ( coupling_nameC);
}

/*----------------------------------------------------------------------------
 *
 *  Get distant points distance to location element
 *
 * parameters
 *   coupling_name        <-- Coupling identifier
 *   distance             --> Distant points distance to location element
 *----------------------------------------------------------------------------*/

void PROCF(cwipi_get_distant_distance_cf,
           CWIPI_GET_DISTANT_DISTANCE_CF) (const char *coupling_name,
                                           const int  *l_coupling_name,
                                           float *distance
                                           ARGF_SUPP_CHAINE)
{
  char *coupling_nameC =
    _cwipi_fortran_to_c_string(coupling_name, *l_coupling_name);

  cwipi::CouplingDataBase & couplingDataBase =
    cwipi::CouplingDataBase::getInstance();

  const std::string &coupling_name_str = coupling_nameC;

  cwipi::oldCoupling& coupling = couplingDataBase.getCoupling(coupling_name_str);

  const float* distanceC = (const float *) coupling.getDistantDistance();
  const int nDistantPoint = coupling.getNDistantPoint();
  for (int i = 0; i < nDistantPoint; i++)
    distance[i] = distanceC[i];

  free ( coupling_nameC);
}

/*----------------------------------------------------------------------------
 *
 * Get number of located distant points
 *
 * parameters
 *   coupling_name           <-- Coupling identifier
 *   n_located_distant_Points --> Number of located distant points
 *
 *----------------------------------------------------------------------------*/

void PROCF(cwipi_get_n_located_dist_pts_cf,
           CWIPI_GET_N_LOCATED_DIST_PTS_CF) (const char *coupling_name,
                                                       const int  *l_coupling_name,
                                                       int *n_located_distant_Points
                                                       ARGF_SUPP_CHAINE)
{
  char *coupling_nameC =
    _cwipi_fortran_to_c_string(coupling_name, *l_coupling_name);

  cwipi::CouplingDataBase & couplingDataBase =
    cwipi::CouplingDataBase::getInstance();

  const std::string &coupling_name_str = coupling_nameC;

  cwipi::oldCoupling& coupling = couplingDataBase.getCoupling(coupling_name_str);

  *n_located_distant_Points = coupling.getNDistantPoint();

  free ( coupling_nameC);
}


/*----------------------------------------------------------------------------
 *
 * Get distant points coordinates
 *
 * parameters
 *   coupling_name        <-- Coupling identifier
 *   coordinates          --> Distant points coordinates 
 *
 *----------------------------------------------------------------------------*/

void PROCF(cwipi_get_dis_coord_cf,
           CWIPI_GET_DIS_COORD_CF)(const char *coupling_name,
                                   const int  *l_coupling_name,
                                   double *coordinates
                                   ARGF_SUPP_CHAINE)
{
  char *coupling_nameC =
    _cwipi_fortran_to_c_string(coupling_name, *l_coupling_name);

  cwipi::CouplingDataBase & couplingDataBase =
    cwipi::CouplingDataBase::getInstance();

  const std::string &coupling_name_str = coupling_nameC;

  cwipi::oldCoupling& coupling = couplingDataBase.getCoupling(coupling_name_str);

  const double* coordinatesC = coupling.getDistantPointCoordinates();
  const int nDistantPoint = coupling.getNDistantPoint();

  for (int i = 0; i < 3 * nDistantPoint; i++)
    coordinates[i] = coordinatesC[i];

  free ( coupling_nameC);
}


/*----------------------------------------------------------------------------
 *
 * Get located points barycentric coordinates index
 *
 * parameters
 *   coupling_name                <-- Coupling identifier
 *   barycentricCoordinatesIndex  --> located points barycentric
 *                                    coordinates index
 *----------------------------------------------------------------------------*/

void PROCF(cwipi_get_dis_bary_coord_idx_cf,
           CWIPI_GET_DIS_BARY_COORD_IDX_CF) (const char *coupling_name,
                                       const int  *l_coupling_name,
                                       int *barycentricCoordinatesIndex
                                       ARGF_SUPP_CHAINE)
{
  char *coupling_nameC =
    _cwipi_fortran_to_c_string(coupling_name, *l_coupling_name);

  cwipi::CouplingDataBase & couplingDataBase =
    cwipi::CouplingDataBase::getInstance();

  const std::string &coupling_name_str = coupling_nameC;

  cwipi::oldCoupling& coupling = couplingDataBase.getCoupling(coupling_name_str);

  const int* barycentricCoordinatesIndexC = coupling.getDistantBarycentricCoordinatesIndex();
  const int nDistantPoint = coupling.getNDistantPoint();

  for (int i = 0; i < nDistantPoint + 1; i++)
    barycentricCoordinatesIndex[i] = barycentricCoordinatesIndexC[i];

  free ( coupling_nameC);
}

/*----------------------------------------------------------------------------
 *
 * Get number of located points
 *
 * parameters
 *   coupling_name        <-- Coupling identifier
 *   n_located_Points     --> Number of located points
 *----------------------------------------------------------------------------*/

void PROCF(cwipi_get_n_located_pts_cf,
           CWIPI_GET_N_LOCATED_PTS_CF) (const char *coupling_name,
                                               const int  *l_coupling_name,
                                               int *n_located_points
                                               ARGF_SUPP_CHAINE)
{
  char *coupling_nameC =
    _cwipi_fortran_to_c_string(coupling_name, *l_coupling_name);

  cwipi::CouplingDataBase & couplingDataBase =
    cwipi::CouplingDataBase::getInstance();

  const std::string &coupling_name_str = coupling_nameC;

  cwipi::oldCoupling& coupling = couplingDataBase.getCoupling(coupling_name_str);

  *n_located_points = coupling.getNLocatedPoint();

  free ( coupling_nameC);
}

/*----------------------------------------------------------------------------
 *
 * Get number of not located points
 *
 * parameters
 *   coupling_name        <-- Coupling identifier
 *   n_not_located_Points --> Number of not located points
 *----------------------------------------------------------------------------*/

void PROCF(cwipi_get_n_not_located_pts_cf,
           CWIPI_GET_N_NOT_LOCATED_PTS_CF) (const char *coupling_name,
                                                   const int  *l_coupling_name,
                                                   int *n_not_located_points
                                                   ARGF_SUPP_CHAINE)
{
  char *coupling_nameC =
    _cwipi_fortran_to_c_string(coupling_name, *l_coupling_name);

  cwipi::CouplingDataBase & couplingDataBase =
    cwipi::CouplingDataBase::getInstance();

  const std::string &coupling_name_str = coupling_nameC;

  cwipi::oldCoupling& coupling = couplingDataBase.getCoupling(coupling_name_str);

  *n_not_located_points = coupling.getNNotlocatedPoint();

  free ( coupling_nameC);
}

/*----------------------------------------------------------------------------
 *
 * Get located points barycentric coordinates
 *
 * parameters
 *   coupling_name                <-- Coupling identifier
 *   barycentricCoordinatesIndex  --> located points barycentric
 *                                    coordinates
 *----------------------------------------------------------------------------*/

void PROCF(cwipi_get_dis_bary_coord_cf,
           CWIPI_GET_DIS_BARY_COORD_CF) (const char *coupling_name,
                                       const int  *l_coupling_name,
                                       double *barycentricCoordinates
                                       ARGF_SUPP_CHAINE)
{
  char *coupling_nameC =
    _cwipi_fortran_to_c_string(coupling_name, *l_coupling_name);

  cwipi::CouplingDataBase & couplingDataBase =
    cwipi::CouplingDataBase::getInstance();

  const std::string &coupling_name_str = coupling_nameC;

  cwipi::oldCoupling& coupling =
    couplingDataBase.getCoupling(coupling_name_str);

  const int* barycentricCoordinatesIndexC = coupling.getDistantBarycentricCoordinatesIndex();
  const double* barycentricCoordinatesC = coupling.getDistantBarycentricCoordinates();
  const int nDistantPoint = coupling.getNDistantPoint();

  for (int i = 0; i < barycentricCoordinatesIndexC[nDistantPoint]; i++)
    barycentricCoordinates[i] = barycentricCoordinatesC[i];

  free ( coupling_nameC);
}

/*----------------------------------------------------------------------------
 *
 * Get number of distant ranks 
 *
 * parameters
 *   coupling_id          <-- Coupling identifier
 *
 * return
 *                        --> Number of distant ranks
 *
 *----------------------------------------------------------------------------*/

void PROCF(cwipi_get_n_dis_ranks_cf, CWIPI_GET_N_DIS_RANKS_CF)
(const char *coupling_name,
 const int  *l_coupling_name,
 int *n_dis_ranks
 ARGF_SUPP_CHAINE)
{
  char *coupling_nameC =
    _cwipi_fortran_to_c_string(coupling_name, *l_coupling_name);

  cwipi::CouplingDataBase & couplingDataBase =
    cwipi::CouplingDataBase::getInstance();

  const std::string &coupling_name_str = coupling_nameC;

  cwipi::oldCoupling& coupling =
    couplingDataBase.getCoupling(coupling_name_str);

  *n_dis_ranks = coupling.getNDistantRank();

  free ( coupling_nameC);
}


/*----------------------------------------------------------------------------
 *
 * Get distant point distribution on distant ranks (size = n_distant_rank + 1)
 *
 * parameters
 *   coupling_id          <-- Coupling identifier
 *
 * return
 *                             Distant point distribution on distant ranks
 *
 *----------------------------------------------------------------------------*/

void PROCF(cwipi_get_dis_distrib_cf, CWIPI_GET_DIS_DISTRIB_CF) 
(const char *coupling_name,
 const int  *l_coupling_name,
 int *distrib
 ARGF_SUPP_CHAINE)
{
  char *coupling_nameC =
    _cwipi_fortran_to_c_string(coupling_name, *l_coupling_name);

  cwipi::CouplingDataBase & couplingDataBase =
    cwipi::CouplingDataBase::getInstance();

  const std::string &coupling_name_str = coupling_nameC;

  cwipi::oldCoupling& coupling =
    couplingDataBase.getCoupling(coupling_name_str);

  const int n_dis_ranks = coupling.getNDistantRank();
  const int *distrib_c = coupling.getDistantDistribution();

  for (int i = 0; i < n_dis_ranks + 1; i++)
    distrib[i] = distrib_c[i];

  free ( coupling_nameC);
}


/*----------------------------------------------------------------------------
 *
 * Get located points distribution on distant ranks (size = n_distant_rank + 1)
 *
 * parameters
 *   coupling_id          <-- Coupling identifier
 *
 * return
 *                            Located points distribution
 *
 *----------------------------------------------------------------------------*/

void PROCF(cwipi_get_loc_pts_distrib_cf, CWIPI_GET_LOC_PTS_DISTRIB_CF)
(const char *coupling_name,
 const int  *l_coupling_name,
 int *distrib
 ARGF_SUPP_CHAINE)
{
  char *coupling_nameC =
    _cwipi_fortran_to_c_string(coupling_name, *l_coupling_name);

  cwipi::CouplingDataBase & couplingDataBase =
    cwipi::CouplingDataBase::getInstance();

  const std::string &coupling_name_str = coupling_nameC;

  cwipi::oldCoupling& coupling =
    couplingDataBase.getCoupling(coupling_name_str);

  const int n_dis_ranks = coupling.getNDistantRank();
  const int *distrib_c = coupling.getLocatedPointsDistribution();

  for (int i = 0; i < n_dis_ranks + 1; i++)
    distrib[i] = distrib_c[i];

  free ( coupling_nameC);
}


/*----------------------------------------------------------------------------
 *
 * Exchange data with the coupled application. This is a synchronization point
 * with the coupled application
 *
 * parameters
 *   coupling_id          <-- Coupling identifier
 *   exchange_name        <-- Exchange name
 *   exchange_type        <-- Exchange type
 *   exchange_dimension   <-- Dimension of exchanged data :
 *                            - CWIPI_DIMENSION_SCALAR
 *                            - CWIPI_DIMENSION_INTERLACED_VECTOR
 *   time_step            <-- Time step  (only for visualization)
 *   time_value           <-- Time value (only for visualization)
 *   sending_field_name   <-- Sending field name
 *   sending_field        <-- Sending field (NULL -> no sending)
 *   receiving_field_name <-- Receiving field name
 *   receiving_field      --> Receiving field
 *
 * returns :
 *   1 if data were received
 *   0 else
 *
 *----------------------------------------------------------------------------*/

void PROCF(cwipi_exch_with_user_itp_cf,
           CWIPI_EXCH_WITH_USER_ITP_CF)
  (const char      *coupling_name,
   const int       *l_coupling_name,
   const char      *exchange_name,
   const int       *l_exchange_name,
   const int       *stride,
   const int       *n_step,
   const double    *time_value,
   const char      *sending_field_name,
   const int       *l_sending_field_name,
   const double    *sending_field,
   char            *receiving_field_name,
   const int       *l_receiving_field_name,
   double          *receiving_field,
   void            *ptFortranInterpolationFct,
   int             *n_not_located_points,
   int             *exchange_status
   ARGF_SUPP_CHAINE)
{
  char *coupling_nameC =
    _cwipi_fortran_to_c_string(coupling_name, *l_coupling_name);

  char *exchange_nameC =
    _cwipi_fortran_to_c_string(exchange_name, *l_exchange_name);

  char *sending_field_nameC =
    _cwipi_fortran_to_c_string(sending_field_name, *l_sending_field_name);

  char *receiving_field_nameC =
    _cwipi_fortran_to_c_string(receiving_field_name, *l_receiving_field_name);


  cwipi::CouplingDataBase & couplingDataBase =
    cwipi::CouplingDataBase::getInstance();

  const std::string &coupling_name_str = coupling_nameC;

  cwipi::oldCoupling& coupling = couplingDataBase.getCoupling(coupling_name_str);

  *exchange_status =  coupling.exchange(exchange_nameC,
                                        *stride,
                                        *n_step,
                                        *time_value,
                                        sending_field_nameC,
                                        sending_field,
                                        receiving_field_nameC,
                                        receiving_field,
                                        ptFortranInterpolationFct);

  *n_not_located_points = coupling.getNNotlocatedPoint();

  free ( coupling_nameC);
  free ( exchange_nameC);
  free ( sending_field_nameC);
  free ( receiving_field_nameC);
}


void PROCF(cwipi_exchange_cf,
           CWIPI_EXCHANGE_CF)
  (const char      *coupling_name,
   const int       *l_coupling_name,
   const char      *exchange_name,
   const int       *l_exchange_name,
   const int       *stride,
   const int       *n_step,
   const double    *time_value,
   const char      *sending_field_name,
   const int       *l_sending_field_name,
   const double    *sending_field,
   char            *receiving_field_name,
   const int       *l_receiving_field_name,
   double          *receiving_field,
   int             *n_not_located_points,
   int             *exchange_status
   ARGF_SUPP_CHAINE)
{
  char *coupling_nameC =
    _cwipi_fortran_to_c_string(coupling_name, *l_coupling_name);

  char *exchange_nameC =
    _cwipi_fortran_to_c_string(exchange_name, *l_exchange_name);

  char *sending_field_nameC =
    _cwipi_fortran_to_c_string(sending_field_name, *l_sending_field_name);

  char *receiving_field_nameC =
    _cwipi_fortran_to_c_string(receiving_field_name, *l_receiving_field_name);


  cwipi::CouplingDataBase & couplingDataBase =
    cwipi::CouplingDataBase::getInstance();

  const std::string &coupling_name_str = coupling_nameC;

  cwipi::oldCoupling& coupling = couplingDataBase.getCoupling(coupling_name_str);

  *exchange_status =  coupling.exchange(exchange_nameC,
                                        *stride,
                                        *n_step,
                                        *time_value,
                                        sending_field_nameC,
                                        sending_field,
                                        receiving_field_nameC,
                                        receiving_field,
                                        NULL);

  *n_not_located_points = coupling.getNNotlocatedPoint();

  free ( coupling_nameC);
  free ( exchange_nameC);
  free ( sending_field_nameC);
  free ( receiving_field_nameC);
}

void PROCF(cwipi_receive_cf,
           CWIPI_RECEIVE_CF)
  (const char      *coupling_name,
   const int       *l_coupling_name,
   const char      *exchange_name,
   const int       *l_exchange_name,
   const int       *stride,
   const int       *n_step,
   const double    *time_value,
   char            *receiving_field_name,
   const int       *l_receiving_field_name,
   double          *receiving_field,
   int             *n_not_located_points,
   int             *exchange_status
   ARGF_SUPP_CHAINE)
{
  char *coupling_nameC =
    _cwipi_fortran_to_c_string(coupling_name, *l_coupling_name);

  char *exchange_nameC =
    _cwipi_fortran_to_c_string(exchange_name, *l_exchange_name);

  char *receiving_field_nameC =
    _cwipi_fortran_to_c_string(receiving_field_name, *l_receiving_field_name);


  cwipi::CouplingDataBase & couplingDataBase =
    cwipi::CouplingDataBase::getInstance();

  const std::string &coupling_name_str = coupling_nameC;

  cwipi::oldCoupling& coupling = couplingDataBase.getCoupling(coupling_name_str);

  *exchange_status =  coupling.exchange(exchange_nameC,
                                        *stride,
                                        *n_step,
                                        *time_value,
                                        NULL,
                                        NULL,
                                        receiving_field_nameC,
                                        receiving_field,
                                        NULL);
  *n_not_located_points = coupling.getNNotlocatedPoint();

  free ( coupling_nameC);
  free ( exchange_nameC);
  free ( receiving_field_nameC);
}

void PROCF(cwipi_send_with_user_itp_cf,
           CWIPI_SEND_WITH_USER_ITP_CF)
  (const char      *coupling_name,
   const int       *l_coupling_name,
   const char      *exchange_name,
   const int       *l_exchange_name,
   const int       *stride,
   const int       *n_step,
   const double    *time_value,
   const char      *sending_field_name,
   const int       *l_sending_field_name,
   const double    *sending_field,
   void            *ptFortranInterpolationFct,
   int             *exchange_status
   ARGF_SUPP_CHAINE)
{
  char *coupling_nameC =
    _cwipi_fortran_to_c_string(coupling_name, *l_coupling_name);

  char *exchange_nameC =
    _cwipi_fortran_to_c_string(exchange_name, *l_exchange_name);

  char *sending_field_nameC =
    _cwipi_fortran_to_c_string(sending_field_name, *l_sending_field_name);

  cwipi::CouplingDataBase & couplingDataBase =
    cwipi::CouplingDataBase::getInstance();

  const std::string &coupling_name_str = coupling_nameC;

  cwipi::oldCoupling& coupling = couplingDataBase.getCoupling(coupling_name_str);

  *exchange_status =  coupling.exchange(exchange_nameC,
                                        *stride,
                                        *n_step,
                                        *time_value,
                                        sending_field_nameC,
                                        sending_field,
                                        NULL,
                                        NULL,
                                        ptFortranInterpolationFct);
  free ( coupling_nameC);
  free ( exchange_nameC);
  free ( sending_field_nameC);
}

void PROCF(cwipi_send_cf,
           CWIPI_SEND_CF)
  (const char      *coupling_name,
   const int       *l_coupling_name,
   const char      *exchange_name,
   const int       *l_exchange_name,
   const int       *stride,
   const int       *n_step,
   const double    *time_value,
   const char      *sending_field_name,
   const int       *l_sending_field_name,
   const double    *sending_field,
   int             *exchange_status
   ARGF_SUPP_CHAINE)
{
  char *coupling_nameC =
    _cwipi_fortran_to_c_string(coupling_name, *l_coupling_name);

  char *exchange_nameC =
    _cwipi_fortran_to_c_string(exchange_name, *l_exchange_name);

  char *sending_field_nameC =
    _cwipi_fortran_to_c_string(sending_field_name, *l_sending_field_name);

  cwipi::CouplingDataBase & couplingDataBase =
    cwipi::CouplingDataBase::getInstance();

  const std::string &coupling_name_str = coupling_nameC;

  cwipi::oldCoupling& coupling = couplingDataBase.getCoupling(coupling_name_str);

  *exchange_status =  coupling.exchange(exchange_nameC,
                                        *stride,
                                        *n_step,
                                        *time_value,
                                        sending_field_nameC,
                                        sending_field,
                                        NULL,
                                        NULL,
                                        NULL);
  free ( coupling_nameC);
  free ( exchange_nameC);
  free ( sending_field_nameC);
}


void PROCF(cwipi_issend_cf,
           CWIPI_ISSEND_CF)
  (const char      *coupling_name,
   const int       *l_coupling_name,
   const char      *exchange_name,
   const int       *l_exchange_name,
   const int       *tag,
   const int       *stride,
   const int       *n_step,
   const double    *time_value,
   const char      *sending_field_name,
   const int       *l_sending_field_name,
   const double    *sending_field,
   int             *request
   ARGF_SUPP_CHAINE)
{
  char *coupling_nameC =
    _cwipi_fortran_to_c_string(coupling_name, *l_coupling_name);

  char *exchange_nameC =
    _cwipi_fortran_to_c_string(exchange_name, *l_exchange_name);

  char *sending_field_nameC =
    _cwipi_fortran_to_c_string(sending_field_name, *l_sending_field_name);

  cwipi::CouplingDataBase & couplingDataBase =
    cwipi::CouplingDataBase::getInstance();

  const std::string &coupling_name_str = coupling_nameC;

  cwipi::oldCoupling& coupling = couplingDataBase.getCoupling(coupling_name_str);

  coupling.issend(exchange_nameC,
                  *tag,
                  *stride,
                  *n_step,
                  *time_value,
                  sending_field_nameC,
                  sending_field,
                  NULL,
                  request);

  free ( coupling_nameC);
  free ( exchange_nameC);
  free ( sending_field_nameC);
}

void PROCF(cwipi_issend_with_user_itp_cf,
           CWIPI_ISSEND_WITH_USER_ITP_CF)
  (const char      *coupling_name,
   const int       *l_coupling_name,
   const char      *exchange_name,
   const int       *l_exchange_name,
   const int       *tag,
   const int       *stride,
   const int       *n_step,
   const double    *time_value,
   const char      *sending_field_name,
   const int       *l_sending_field_name,
   const double    *sending_field,
   void            *ptFortranInterpolationFct,
   int             *request
   ARGF_SUPP_CHAINE)
{
  char *coupling_nameC =
    _cwipi_fortran_to_c_string(coupling_name, *l_coupling_name);

  char *exchange_nameC =
    _cwipi_fortran_to_c_string(exchange_name, *l_exchange_name);

  char *sending_field_nameC =
    _cwipi_fortran_to_c_string(sending_field_name, *l_sending_field_name);

  cwipi::CouplingDataBase & couplingDataBase =
    cwipi::CouplingDataBase::getInstance();

  const std::string &coupling_name_str = coupling_nameC;

  cwipi::oldCoupling& coupling = couplingDataBase.getCoupling(coupling_name_str);

  coupling.issend(exchange_nameC,
                  *tag,
                  *stride,
                  *n_step,
                  *time_value,
                  sending_field_nameC,
                  sending_field,
                  ptFortranInterpolationFct,
                  request);

  free ( coupling_nameC);
  free ( exchange_nameC);
  free ( sending_field_nameC);
}

void PROCF(cwipi_ireceive_cf, 
           CWIPI_IRECEIVE_CF)
  (const char      *coupling_name,
   const int       *l_coupling_name,
   const char      *exchange_name,
   const int       *l_exchange_name,
   const int       *tag,
   const int       *stride,
   const int       *n_step,
   const double    *time_value,
   char            *receiving_field_name,
   const int       *l_receiving_field_name,
   double          *receiving_field,
   int             *request
   ARGF_SUPP_CHAINE)
{
  char *coupling_nameC =
    _cwipi_fortran_to_c_string(coupling_name, *l_coupling_name);

  char *exchange_nameC =
    _cwipi_fortran_to_c_string(exchange_name, *l_exchange_name);

  char *receiving_field_nameC =
    _cwipi_fortran_to_c_string(receiving_field_name, *l_receiving_field_name);

  cwipi::CouplingDataBase & couplingDataBase =
    cwipi::CouplingDataBase::getInstance();

  const std::string &coupling_name_str = coupling_nameC;

  cwipi::oldCoupling& coupling = couplingDataBase.getCoupling(coupling_name_str);

  coupling.irecv(exchange_nameC,
                 *tag,
                 *stride,
                 *n_step,
                 *time_value,
                 receiving_field_nameC,
                 receiving_field,
                 request);

  free ( coupling_nameC);
  free ( exchange_nameC);
  free ( receiving_field_nameC);
}


void PROCF(cwipi_wait_irecv_cf, CWIPI_WAIT_IRECV_CF)
  (const char      *coupling_name,
   const int       *l_coupling_name,
   const int       *request
   ARGF_SUPP_CHAINE)
{
  char *coupling_nameC =
    _cwipi_fortran_to_c_string(coupling_name, *l_coupling_name);

  cwipi::CouplingDataBase & couplingDataBase =
    cwipi::CouplingDataBase::getInstance();

  const std::string &coupling_name_str = coupling_nameC;

  cwipi::oldCoupling& coupling = couplingDataBase.getCoupling(coupling_name_str);

  coupling.waitIrecv(*request);

  free ( coupling_nameC);
}

void PROCF(cwipi_wait_issend_cf, CWIPI_WAIT_ISSEND_CF)
  (const char      *coupling_name,
   const int       *l_coupling_name,
   const int       *request
   ARGF_SUPP_CHAINE)
{
  char *coupling_nameC =
    _cwipi_fortran_to_c_string(coupling_name, *l_coupling_name);

  cwipi::CouplingDataBase & couplingDataBase =
    cwipi::CouplingDataBase::getInstance();

  const std::string &coupling_name_str = coupling_nameC;

  cwipi::oldCoupling& coupling = couplingDataBase.getCoupling(coupling_name_str);

  coupling.waitIssend(*request);

  free ( coupling_nameC);
}



/*----------------------------------------------------------------------------
 *
 * Delete a coupling
 *
 * parameters
 *   coupling_id          <-- Coupling identifier
 *
 *----------------------------------------------------------------------------*/

void PROCF(cwipi_delete_coupling_cf,
           CWIPI_DELETE_COUPLING_CF)
  (const char *coupling_name,
   const int       *l_coupling_name
   ARGF_SUPP_CHAINE)
{
  char *coupling_nameC =
    _cwipi_fortran_to_c_string(coupling_name, *l_coupling_name);

  cwipi_delete_coupling(coupling_nameC);
  free ( coupling_nameC);
}

/*----------------------------------------------------------------------------
 *
 * Finalize cwipi. This is a synchronization point between all applications
 *
 *----------------------------------------------------------------------------*/


void PROCF(cwipi_finalize_cf,
           CWIPI_FINALIZE_CF) ()
{
  cwipi_finalize();
}

/*----------------------------------------------------------------------------
 *
 * Dump application properties
 *
 *----------------------------------------------------------------------------*/

void PROCF(cwipi_dump_appli_properties_cf,
           CWIPI_DUMP_APPLI_PROPERTIES_CF) ()
{
  cwipi_dump_application_properties();
}

/*----------------------------------------------------------------------------
 *
 * Get not located points
 *
 * parameters
 *   coupling_id          <-- Coupling identifier
 *   notLocatedPoints     --> Not located points
 *
 *----------------------------------------------------------------------------*/

void PROCF(cwipi_get_not_located_pts_cf,
           CWIPI_GET_NOT_LOCATED_PTS_CF)
  (const char *coupling_name,
   const int  *l_coupling_name,
   int *notLocatedPoints)
{
  cwipi::CouplingDataBase & couplingDataBase =
    cwipi::CouplingDataBase::getInstance();

  char *coupling_nameC =
    _cwipi_fortran_to_c_string(coupling_name, *l_coupling_name);

  const std::string &coupling_name_str = coupling_nameC;

  cwipi::oldCoupling& coupling = couplingDataBase.getCoupling(coupling_name_str);

  const int n_not_located_points = coupling.getNNotlocatedPoint();
  const int *notLocatedPointsC = coupling.getNotlocatedPoint();

  for( int i = 0; i <  n_not_located_points; i++)
    notLocatedPoints[i] = notLocatedPointsC[i];

  free ( coupling_nameC);

}

/*----------------------------------------------------------------------------
 *
 * Get located points
 *
 * parameters
 *   coupling_id          <-- Coupling identifier
 *   notLocatedPoints     --> Not located points
 *
 *----------------------------------------------------------------------------*/

void PROCF(cwipi_get_located_pts_cf,
           CWIPI_GET_LOCATED_PTS_CF)
  (const char *coupling_name,
   const int  *l_coupling_name,
   int *locatedPoints)
{
  cwipi::CouplingDataBase & couplingDataBase =
    cwipi::CouplingDataBase::getInstance();

  char *coupling_nameC =
    _cwipi_fortran_to_c_string(coupling_name, *l_coupling_name);

  const std::string &coupling_name_str = coupling_nameC;

  cwipi::oldCoupling& coupling = couplingDataBase.getCoupling(coupling_name_str);

  const int n_located_points = coupling.getNLocatedPoint();
  const int *locatedPointsC = coupling.getLocatedPoint();

  for( int i = 0; i < n_located_points; i++)
    locatedPoints[i] = locatedPointsC[i];

  free ( coupling_nameC);

}

/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
}
#endif /* __cplusplus */
