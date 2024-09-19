/*
  This file is part of the CWIPI library.

  Copyright (C) 2011  ONERA

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
#include <stdarg.h>

/*----------------------------------------------------------------------------
 * BFT library headers
 *----------------------------------------------------------------------------*/

#include "bftc_mem.h"
#include "bftc_printf.h"

/*----------------------------------------------------------------------------
 * FVM library headers
 *----------------------------------------------------------------------------*/

#include "fvmc_parall.h"
#include "fvmc_ho_basis.h"
#include "fvmc_ho_location.h"

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "cwipi.h"
#include "cwipi_config.h"
#include "applicationPropertiesDataBase.hxx"
#include "applicationPropertiesDataBase_i.hxx"
#include "couplingDataBase.hxx"
#include "couplingDataBase_i.hxx"
#include "oldCoupling.hxx"
#include "oldCoupling_i.hxx"
//#include "conservativeMesh.hxx"


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

/*============================================================================
 * Private function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 *
 * bftc_printf proxy setting for C interface
 *
 *----------------------------------------------------------------------------*/

static int _cwipi_print_with_c
(
 const char     *const format,
       va_list         arg_ptr
)
{
  return vfprintf(_cwipi_output_listing, format, arg_ptr);
}

static int _cwipi_flush_output_listing(void)
{
  return fflush(_cwipi_output_listing);
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
 * This is a synchronization point between all applications
 *----------------------------------------------------------------------------*/

void cwipi_init
(const MPI_Comm                           common_comm       ,
 const char                               *application_name ,
 MPI_Comm                                 *application_comm )

{

  cwipi::ApplicationPropertiesDataBase & properties =
    cwipi::ApplicationPropertiesDataBase::getInstance();

  bftc_printf("\ncwipi " CWIPI_VERSION " initializing\n");
  bftc_printf("------------------------\n\n");
  bftc_printf_flush();

  *application_comm = properties.init(application_name,
                                      common_comm);
}

/*----------------------------------------------------------------------------
 *
 * Set up the file used for the output listing
 *
 * parameters:
 *   output_listing      <-- Output listing file (C function)
 *----------------------------------------------------------------------------*/

void cwipi_set_output_listing
(FILE *output_listing)
{
  _cwipi_output_listing = output_listing;
  bftc_printf_proxy_set(_cwipi_print_with_c);
  bftc_printf_flush_proxy_set(_cwipi_flush_output_listing);
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

void cwipi_add_local_int_control_parameter(const char *name, int initial_value)
{
  cwipi::ApplicationPropertiesDataBase & properties =
    cwipi::ApplicationPropertiesDataBase::getInstance();

  properties.addLocalIntControlParameter(name, initial_value);
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

void cwipi_add_local_double_control_parameter
(const char *name,
 double initial_value)
{
  cwipi::ApplicationPropertiesDataBase & properties =
    cwipi::ApplicationPropertiesDataBase::getInstance();

  properties.addLocalDoubleControlParameter(name, initial_value);
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

void cwipi_add_local_string_control_parameter
(const char *name,
 const char *initial_value)
{
  cwipi::ApplicationPropertiesDataBase & properties =
    cwipi::ApplicationPropertiesDataBase::getInstance();

  properties.addLocalStringControlParameter(name, initial_value);
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

void cwipi_set_local_int_control_parameter(const char *name, int value)
{
  cwipi::ApplicationPropertiesDataBase & properties =
    cwipi::ApplicationPropertiesDataBase::getInstance();

  properties.setLocalIntControlParameter(name, value);
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

void cwipi_set_local_double_control_parameter(const char *name, double value)
{
  cwipi::ApplicationPropertiesDataBase & properties =
    cwipi::ApplicationPropertiesDataBase::getInstance();

  properties.setLocalDoubleControlParameter(name, value);
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

void cwipi_set_local_string_control_parameter(const char *name, const char *value)
{
  cwipi::ApplicationPropertiesDataBase & properties =
    cwipi::ApplicationPropertiesDataBase::getInstance();

  properties.setLocalStringControlParameter(name, value);
}

/*----------------------------------------------------------------------------
 *
 * Get a integer control parameter of the current application
 *
 * parameters
 *    name           <-- parameter name
 *
 *----------------------------------------------------------------------------*/

int cwipi_get_local_int_control_parameter(const char *name)
{
  cwipi::ApplicationPropertiesDataBase & properties =
    cwipi::ApplicationPropertiesDataBase::getInstance();

  return properties.getLocalIntControlParameter(name);
}

/*----------------------------------------------------------------------------
 *
 * Get a double control parameter of the current application
 *
 * parameters
 *    name           <-- parameter name
 *
 *----------------------------------------------------------------------------*/

double cwipi_get_local_double_control_parameter(const char *name)
{
  cwipi::ApplicationPropertiesDataBase & properties =
    cwipi::ApplicationPropertiesDataBase::getInstance();

  return properties.getLocalDoubleControlParameter(name);
}

/*----------------------------------------------------------------------------
 *
 * Get a string control parameter of the current application
 *
 * parameters
 *    name           <-- parameter name
 *
 *----------------------------------------------------------------------------*/

const char* cwipi_get_local_string_control_parameter(const char *name)
{
  cwipi::ApplicationPropertiesDataBase & properties =
    cwipi::ApplicationPropertiesDataBase::getInstance();

  return (properties.getLocalStringControlParameter(name).c_str());
}

/*----------------------------------------------------------------------------
 *
 * Delete a current application Int parameter
 *
 * parameters
 *    name           <-- parameter name
 *
 *----------------------------------------------------------------------------*/

void cwipi_delete_local_int_control_parameter(const char *name)
{
  cwipi::ApplicationPropertiesDataBase & properties =
    cwipi::ApplicationPropertiesDataBase::getInstance();

  properties.eraseLocalIntControlParameter(name);
}

/*----------------------------------------------------------------------------
 *
 * Delete a current application double parameter
 *
 * parameters
 *    name           <-- parameter name
 *
 *----------------------------------------------------------------------------*/

void cwipi_delete_local_double_control_parameter(const char *name)
{
  cwipi::ApplicationPropertiesDataBase & properties =
    cwipi::ApplicationPropertiesDataBase::getInstance();

  properties.eraseLocalDoubleControlParameter(name);
}


/*----------------------------------------------------------------------------
 *
 * Delete a current application string parameter
 *
 * parameters
 *    name           <-- parameter name
 *
 *----------------------------------------------------------------------------*/

void cwipi_delete_local_string_control_parameter(const char *name)
{
  cwipi::ApplicationPropertiesDataBase & properties =
    cwipi::ApplicationPropertiesDataBase::getInstance();

  properties.eraseLocalStringControlParameter(name);
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

int cwipi_get_distant_int_control_parameter
(const char *application_name,
 const char *name)
{
  cwipi::ApplicationPropertiesDataBase & properties =
    cwipi::ApplicationPropertiesDataBase::getInstance();

  return properties.getDistantIntControlParameter(application_name, name);
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

double cwipi_get_distant_double_control_parameter
(const char *application_name,
 const char *name)
{
  cwipi::ApplicationPropertiesDataBase & properties =
    cwipi::ApplicationPropertiesDataBase::getInstance();

  return properties.getDistantDoubleControlParameter(application_name, name);
}

/*----------------------------------------------------------------------------
 *
 * Get a string control parameter of a other application
 *
 * parameters
 *    application_name    <-- application name
 *    name                <-- parameter name
 *
 *----------------------------------------------------------------------------*/

const char *cwipi_get_distant_string_control_parameter
(const char *application_name,
 const char *name)
{
  cwipi::ApplicationPropertiesDataBase & properties =
    cwipi::ApplicationPropertiesDataBase::getInstance();

  return properties.getDistantStringControlParameter(application_name, name).c_str();
}

/*----------------------------------------------------------------------------
 *
 * Has int parameter
 *
 * parameters
 *    application_name       <-- application name
 *    name                   <-- parameter name
 *
 * return
 *    1 : true / 0 : false
 *----------------------------------------------------------------------------*/

int cwipi_has_int_parameter
(const char *application_name,
 const char *name)
{
  cwipi::ApplicationPropertiesDataBase & properties =
    cwipi::ApplicationPropertiesDataBase::getInstance();

  if (!strcmp(application_name,  properties.getLocalName().c_str()))
    return properties.hasLocalIntControlParameter(name);
  else
    return properties.hasDistantIntControlParameter(application_name,name);
}

/*----------------------------------------------------------------------------
 *
 * Has double parameter
 *
 * parameters
 *    application_name       <-- application name
 *    name                   <-- parameter name
 *
 * return
 *    1 : true / 0 : false
 *----------------------------------------------------------------------------*/

int cwipi_has_double_parameter
(const char *application_name,
 const char *name)
{
  cwipi::ApplicationPropertiesDataBase & properties =
    cwipi::ApplicationPropertiesDataBase::getInstance();

  if (!strcmp(application_name,  properties.getLocalName().c_str()))
    return properties.hasLocalDoubleControlParameter(name);
  else
    return properties.hasDistantDoubleControlParameter(application_name,name);
}

/*----------------------------------------------------------------------------
 *
 * Has string parameter
 *
 * parameters
 *    application_name       <-- application name
 *    name                   <-- parameter name
 *
 * return
 *    1 : true / 0 : false
 *----------------------------------------------------------------------------*/

int cwipi_has_string_parameter
(const char *application_name,
 const char *name)
{
  cwipi::ApplicationPropertiesDataBase & properties =
    cwipi::ApplicationPropertiesDataBase::getInstance();

  if (!strcmp(application_name,  properties.getLocalName().c_str()))
    return properties.hasLocalStringControlParameter(name);
  else
    return properties.hasDistantStringControlParameter(application_name,name);
}

/*----------------------------------------------------------------------------
 *
 * Get number of int parameters
 *
 * parameters
 *    application_name       <-- application name
 *
 * return
 *    Number of int parameters
 *----------------------------------------------------------------------------*/

int cwipi_get_n_int_parameters
(const char *application_name)
{
  cwipi::ApplicationPropertiesDataBase & properties =
    cwipi::ApplicationPropertiesDataBase::getInstance();

  if (!strcmp(application_name,  properties.getLocalName().c_str()))
    return properties.getLocalNIntControlParameter();
  else
    return properties.getDistantNIntControlParameter(application_name);
}

/*----------------------------------------------------------------------------
 *
 * Get number of double parameters
 *
 * parameters
 *    application_name       <-- application name
 *
 * return
 *    Number of double parameters
 *----------------------------------------------------------------------------*/

int cwipi_get_n_double_parameters
(const char *application_name)
{
  cwipi::ApplicationPropertiesDataBase & properties =
    cwipi::ApplicationPropertiesDataBase::getInstance();

  if (!strcmp(application_name,  properties.getLocalName().c_str()))
    return properties.getLocalNDoubleControlParameter();
  else
    return properties.getDistantNDoubleControlParameter(application_name);
}

/*----------------------------------------------------------------------------
 *
 * Get number of string parameters
 *
 * parameters
 *    application_name       <-- application name
 *
 * return
 *    Number of string parameters
 *----------------------------------------------------------------------------*/

int cwipi_get_n_string_parameters
(const char *application_name)
{
  cwipi::ApplicationPropertiesDataBase & properties =
    cwipi::ApplicationPropertiesDataBase::getInstance();

  if (!strcmp(application_name,  properties.getLocalName().c_str()))
    return properties.getLocalNStringControlParameter();
  else
    return properties.getDistantNStringControlParameter(application_name);
}

/*----------------------------------------------------------------------------
 *
 * Get list int parameters
 *
 * parameters
 *    application_name       <-- application name
 *
 * return
 *    parameters name
 *----------------------------------------------------------------------------*/

char** cwipi_get_list_int_parameters
(const char *application_name)
{
  cwipi::ApplicationPropertiesDataBase & properties =
    cwipi::ApplicationPropertiesDataBase::getInstance();

  if (!strcmp(application_name,  properties.getLocalName().c_str()))
    return properties.getLocalListIntControlParameter();
  else
    return properties.getDistantListIntControlParameter(application_name);
}

/*----------------------------------------------------------------------------
 *
 * Get list double parameters
 *
 * parameters
 *    application_name       <-- application name
 *
 * return
 *    parameters name
 *----------------------------------------------------------------------------*/

char** cwipi_get_list_double_parameters
(const char *application_name)
{
  cwipi::ApplicationPropertiesDataBase & properties =
    cwipi::ApplicationPropertiesDataBase::getInstance();

  if (!strcmp(application_name,  properties.getLocalName().c_str()))
    return properties.getLocalListDoubleControlParameter();
  else
    return properties.getDistantListDoubleControlParameter(application_name);
}

/*----------------------------------------------------------------------------
 *
 * Get list string parameters
 *
 * parameters
 *    application_name       <-- application name
 *
 * return
 *    parameters name
 *----------------------------------------------------------------------------*/

char** cwipi_get_list_string_parameters
(const char *application_name)
{
  cwipi::ApplicationPropertiesDataBase & properties =
    cwipi::ApplicationPropertiesDataBase::getInstance();

  if (!strcmp(application_name,  properties.getLocalName().c_str()))
    return properties.getLocalListStringControlParameter();
  else
    return properties.getDistantListStringControlParameter(application_name);
}

/*----------------------------------------------------------------------------
 *
 * Synchronize local control parameters with an other application.
 *  This is a synchronization point with this second application
 *
 * parameters
 *    application_name    <-- application name
 *    name                <-- parameter name
 *
 *----------------------------------------------------------------------------*/

void cwipi_synchronize_control_parameter(const char *application_name)
{
  cwipi::ApplicationPropertiesDataBase & properties =
    cwipi::ApplicationPropertiesDataBase::getInstance();

  const std::string &application_nameStr = application_name;
  return properties.mergeParameters(application_nameStr);
}

/*----------------------------------------------------------------------------
 *
 * Create a coupling object
 *
 * parameters:
 *   coupling_name           <-- Coupling identifier
 *   coupling_type           <-- Coupling type
 *   coupled_application     <-- Coupled application name
 *   entitiesDim             <-- Mesh entities dimension (1, 2 or 3)
 *   tolerance               <-- Geometric tolerance to locate
 *   mesh_type               <-- CWIPI_STATIC_MESH
 *                               CWIPI_MOBILE_MESH (not implemented yet)
 *                               CWIPI_CYCLIC_MESH
 *   solver_type             <-- CWIPI_SOLVER_CELL_CENTER
 *                               CWIPI_SOLVER_CELL_VERTEX
 *   output_frequency        <-- Output frequency
 *   output_format           <-- Output format to visualize exchanged fields
 *                               on the coupled mesh. Choice between :
 *                                 - "EnSight Gold"
 *                                 - "MED_fichier"
 *                                 - "CGNS"
 *   output_format_option    <-- Output options
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
 *
 *----------------------------------------------------------------------------*/

void cwipi_create_coupling
( const char  *coupling_name,
  const cwipi_coupling_type_t coupling_type,
  const char  *coupled_application,
  const int entities_dim,
  const double tolerance,
  const cwipi_mesh_type_t mesh_type,
  const cwipi_solver_type_t solver_type,
  const int    output_frequency,
  const char  *output_format,
  const char  *output_format_option
  ...)
{
  // traitement du parametre optionnel

  va_list args;
  int nb_locations = 1;
  va_start(args, output_format_option);
  if(mesh_type == CWIPI_CYCLIC_MESH)
    nb_locations = va_arg(args, int);
  va_end(args);

  cwipi::CouplingDataBase & couplingDataBase =
    cwipi::CouplingDataBase::getInstance();

  cwipi::ApplicationPropertiesDataBase & properties =
    cwipi::ApplicationPropertiesDataBase::getInstance();

  const std::string &coupled_application_str = coupled_application;

  couplingDataBase.createCoupling(coupling_name,
                                  coupling_type,
                                  properties.getLocalApplicationProperties(),
                                  properties.getDistantApplicationProperties(coupled_application_str),
                                  entities_dim,
                                  tolerance,
                                  solver_type,
                                  output_frequency,
                                  output_format,
                                  output_format_option,
				  nb_locations);
}

/*----------------------------------------------------------------------------
 *
 * Set data user (optional)
 *
 * parameters:
 *   coupling_name           <-- Coupling identifier
 *   data                    <-- data user
 *----------------------------------------------------------------------------*/

void
cwipi_set_data_user
(
 const char  *coupling_name,
       void  *data
)
{
  cwipi::CouplingDataBase & couplingDataBase =
    cwipi::CouplingDataBase::getInstance();

  const std::string &coupling_name_str = coupling_name;

  cwipi::oldCoupling& coupling = couplingDataBase.getCoupling(coupling_name_str);

 coupling.set_data_user(data);

}


/*----------------------------------------------------------------------------
 *
 * Get data user (optional)
 *
 * parameters:
 *   coupling_name           <-- Coupling identifier
 *
 * return :
 *   data
 *----------------------------------------------------------------------------*/

void *
cwipi_get_data_user
(
 const char  *coupling_name
)
{
  cwipi::CouplingDataBase & couplingDataBase =
    cwipi::CouplingDataBase::getInstance();

  const std::string &coupling_name_str = coupling_name;

  cwipi::oldCoupling& coupling = couplingDataBase.getCoupling(coupling_name_str);

  return coupling.get_data_user();
}


/*----------------------------------------------------------------------------
 *
 * cwipi_set_location_index
 *
 * parameters
 *   coupling_id          <-- Coupling identifier
 *   index                <-- location index
 *----------------------------------------------------------------------------*/

void cwipi_set_location_index (const char *coupling_name,
			       const int index)
{
  cwipi::CouplingDataBase & couplingDataBase =
    cwipi::CouplingDataBase::getInstance();

  const std::string &coupling_name_str = coupling_name;

  cwipi::oldCoupling& coupling = couplingDataBase.getCoupling(coupling_name_str);

  coupling.setLocationIndex(index);
}
/*----------------------------------------------------------------------------
 *
 * cwipi_save_location
 *
 * parameters
 *   coupling_id          <-- Coupling identifier
 *----------------------------------------------------------------------------*/

void cwipi_save_location(const char *coupling_name)
{
  cwipi::CouplingDataBase & couplingDataBase =
    cwipi::CouplingDataBase::getInstance();

  const std::string &coupling_name_str = coupling_name;

  cwipi::oldCoupling& coupling = couplingDataBase.getCoupling(coupling_name_str);

  coupling.saveLocation();
}


/*----------------------------------------------------------------------------
 *
 * cwipi_load_location
 *
 * parameters
 *   coupling_id          <-- Coupling identifier
 *----------------------------------------------------------------------------*/

void cwipi_load_location(const char *coupling_name)
{
  cwipi::CouplingDataBase & couplingDataBase =
    cwipi::CouplingDataBase::getInstance();

  const std::string &coupling_name_str = coupling_name;

  cwipi::oldCoupling& coupling = couplingDataBase.getCoupling(coupling_name_str);

  coupling.loadLocation();
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

void cwipi_open_location_file (const char *coupling_name,
			       const char *filename,
			       const char *mode)
{
  cwipi::CouplingDataBase & couplingDataBase =
    cwipi::CouplingDataBase::getInstance();

  const std::string &coupling_name_str = coupling_name;

  cwipi::oldCoupling& coupling = couplingDataBase.getCoupling(coupling_name_str);

  coupling.openLocationFile(filename,mode);
}
/*----------------------------------------------------------------------------
 *
 * cwipi_close_location_file
 *
 * parameters
 *   coupling_id          <-- Coupling identifier
 *----------------------------------------------------------------------------*/

void cwipi_close_location_file (const char *coupling_name)
{
  cwipi::CouplingDataBase & couplingDataBase =
    cwipi::CouplingDataBase::getInstance();

  const std::string &coupling_name_str = coupling_name;

  cwipi::oldCoupling& coupling = couplingDataBase.getCoupling(coupling_name_str);

  coupling.closeLocationFile();
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

void cwipi_set_points_to_locate
(const char  *coupling_name,
 const int    n_points,
 double coordinate[])
{
  cwipi::CouplingDataBase & couplingDataBase =
    cwipi::CouplingDataBase::getInstance();

  const std::string &coupling_name_str = coupling_name;

  cwipi::oldCoupling& coupling = couplingDataBase.getCoupling(coupling_name_str);

  coupling.setPointsToLocate(n_points, coordinate);
}

/*----------------------------------------------------------------------------
 *
 * Define the support mesh for a coupling. The connectivity is sorted if
 * necessary.
 *
 *
 * Order definition :
 *    1D : edges
 *    2D : triangles, quadrangles, polygons
 *    3D : tetrahedra, pyramids, prism, hexaedra
 *
 * Local connectivity for the following element type :
 *
 *  - edge :
 *
 *   1 x-------x 2
 *
 *  - triangle :
 *
 *   1 x-------x 3
 *      \     /
 *       \   /
 *        \ /
 *         x 2
 *
 *  - quadrangle :
 *
 *      4 x-------x 3
 *       /       /
 *      /       /
 *   1 x-------x2
 *
 *   - tetrahedra :
 *
 *         x 4
 *        /|\
 *       / | \
 *      /  |  \
 *   1 x- -|- -x 3
 *      \  |  /
 *       \ | /
 *        \|/
 *         x 2
 *
 *   - pyramid :
 *
 *          5 x
 *           /|\
 *          //| \
 *         // |  \
 *      4 x/--|---x 3
 *       //   |  /
 *      //    | /
 *   1 x-------x 2
 *
 *  - prism :
 *
 *   4 x-------x 6
 *     |\     /|
 *     | \   / |
 *   1 x- \-/ -x 3
 *      \ 5x  /
 *       \ | /
 *        \|/
 *         x 2
 *
 *  -  hexaedra :
 *
 *      8 x-------x 7
 *       /|      /|
 *      / |     / |
 *   5 x-------x6 |
 *     | 4x----|--x 3
 *     | /     | /
 *     |/      |/
 *   1 x-------x 2
 *
 * Order definition :
 *    1D : edges
 *    2D : triangles, quadrangles, polygons
 *    3D : tetrahedra, pyramids, prism, hexaedra
 *
 *
 * parameters:
 *   coupling_id        <-- coupling name
 *   n_vertex           <-- number of vertices
 *   n_elements         <-- number of elements
 *   coordinates        <-- vertex interlaced coordinates
 *   connectivity_index <-> element -> vertices index (O to n-1)
 *                          size: n_elements + 1
 *                          (out : ordered connectivity_index)
 *   connectivity       <-> element -> vertex connectivity
 *                          size: connectivity_index[n_elements]
 *                          (out : ordered connectivity)
 *
 *----------------------------------------------------------------------------*/

void cwipi_define_mesh(const char *coupling_name,
                       const int n_vertex,
                       const int n_element,
                       double coordinates[],
                       int connectivity_index[],
                       int connectivity[])
{
  cwipi::CouplingDataBase & couplingDataBase =
    cwipi::CouplingDataBase::getInstance();

  const std::string &coupling_name_str = coupling_name;

  cwipi::oldCoupling& coupling = couplingDataBase.getCoupling(coupling_name_str);

  coupling.defineMesh(n_vertex,
                      n_element,
                      coordinates,
                      connectivity_index,
                      connectivity);
}

void cwipi_ho_define_mesh(const char *coupling_name,
                          const int n_vertex,
                          const int n_element,
                          const int order,
                          double coordinates[],
                          int connectivity_index[],
                          int connectivity[])
{
  cwipi::CouplingDataBase & couplingDataBase =
    cwipi::CouplingDataBase::getInstance();

  const std::string &coupling_name_str = coupling_name;

  cwipi::oldCoupling& coupling = couplingDataBase.getCoupling(coupling_name_str);

  coupling.defineMesh(n_vertex,
                      n_element,
                      coordinates,
                      connectivity_index,
                      connectivity,
                      order);
}



void cwipi_shared_fvmc_nodal(const char *coupling_name,
                            void* fvmc_nodal)
{
  cwipi::CouplingDataBase & couplingDataBase =
    cwipi::CouplingDataBase::getInstance();

  const std::string &coupling_name_str = coupling_name;

  cwipi::oldCoupling& coupling = couplingDataBase.getCoupling(coupling_name_str);

  coupling.defineMesh((fvmc_nodal_t *) fvmc_nodal);

}

void cwipi_add_polyhedra(const char *coupling_name,
                         const int n_element,
                         int face_index[],
                         int cell_to_face_connectivity[],
                         const int n_faces,
                         int face_connectivity_index[],
                         int face_connectivity[])
{
  cwipi::CouplingDataBase & couplingDataBase =
    cwipi::CouplingDataBase::getInstance();

  const std::string &coupling_name_str = coupling_name;

  cwipi::oldCoupling& coupling = couplingDataBase.getCoupling(coupling_name_str);

  coupling.defineMeshAddPolyhedra(n_element,
                                  face_index,
                                  cell_to_face_connectivity,
                                  n_faces,
                                  face_connectivity_index,
                                  face_connectivity);
}

/*----------------------------------------------------------------------------
 *
 * Define specific options for ho elements
 *
 * parameters:
 *   coupling_id     <-- coupling name
 *   option          <-- option name, Choice between :
 *                          - "opt_bbox_step"
 *                              * Description : step of discretization used
 *                                              to compute the optimized element
 *                                              bounding boxes
 *                                              -1 to deactivate this computation
 *                              * Default     : 10
 *   value           <-- option value
 *
 *----------------------------------------------------------------------------*/

void cwipi_ho_options_set (const char *coupling_id,
                           const char *option,
                           const char *value)
{
  cwipi::CouplingDataBase & couplingDataBase =
    cwipi::CouplingDataBase::getInstance();

  const std::string &coupling_name_str = coupling_id;

  cwipi::oldCoupling& coupling = couplingDataBase.getCoupling(coupling_name_str);

  coupling.hoOptionsSet (option, value);
}

/*----------------------------------------------------------------------------
 *
 * Define ho element ordering from the location in the (u, v, w) grid
 *
 * parameters:
 *   coupling_id     <-- coupling name
 *   t_elt           <-- element type
 *   n_nodes         <-- number of nodes
 *   ijk_grid        <-- user ordering to (u, v, w) grid (size = elt_dim * n_nodes)
 *
 *----------------------------------------------------------------------------*/

void cwipi_ho_ordering_from_IJK_set (const char *coupling_id,
                                     const cwipi_element_t t_elt,
                                     const int n_nodes,
                                     const int *ijk_grid)
{


  cwipi::CouplingDataBase & couplingDataBase =
    cwipi::CouplingDataBase::getInstance();

  const std::string &coupling_name_str = coupling_id;

  cwipi::oldCoupling& coupling = couplingDataBase.getCoupling(coupling_name_str);

  coupling.hoOrderingSet (t_elt, n_nodes, ijk_grid);

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

void cwipi_ho_ordering_from_ref_elt_set (const char   *coupling_id,
                                         const cwipi_element_t t_elt,
                                         const int n_nodes,
                                         const double *coords)
{
  cwipi::CouplingDataBase & couplingDataBase =
    cwipi::CouplingDataBase::getInstance();

  const std::string &coupling_name_str = coupling_id;

  cwipi::oldCoupling& coupling = couplingDataBase.getCoupling(coupling_name_str);

  coupling.hoOrderingFromRefEltSet (t_elt, n_nodes, coords);

}

/*----------------------------------------------------------------------------
 *
 * Set elementary functions
 *
 * parameters:
 *   location_tetra    <-- Location in a tetrahedron
 *   location_prism    <-- Location in a prism
 *   location_pyramid  <-- Location in a pyramid
 *   location_hexa     <-- Location in a hexaedron
 *   location_tria     <-- Location on a triangle
 *   location_quad     <-- Location on a quandragle
 *   location_edge     <-- Location on a edge
 *   interp_tetra      <-- Interpolation in a tetrahedron
 *   interp_prism      <-- Interpolation in a prism
 *   interp_pyramid    <-- Interpolation in a pyramid
 *   interp_hexa       <-- Interpolation in a hexaedron
 *   interp_tria       <-- Interpolation on a triangle
 *   interp_quad       <-- Interpolation on a quandragle
 *   interp_edge       <-- Interpolation on a edge
 *
 *----------------------------------------------------------------------------*/

void
cwipi_ho_user_elt_set (cwipi_element_t elt_type,
                       cwipi_ho_basis_fct_t element_basis,
                       cwipi_ho_location_fct_t location_in_element)
{
  fvmc_element_t _elt_type = (fvmc_element_t) 0;

  switch (elt_type) {

  case CWIPI_EDGEHO:
    _elt_type = FVMC_EDGE;
    break;
  case CWIPI_FACE_TRIAHO:
    _elt_type = FVMC_FACE_TRIA;
    break;
  case CWIPI_FACE_QUADHO:
    _elt_type = FVMC_FACE_QUAD;
    break;
  case CWIPI_CELL_TETRAHO:
    _elt_type = FVMC_CELL_TETRA;
    break;
  case CWIPI_CELL_HEXAHO:
    _elt_type = FVMC_CELL_HEXA;
    break;
  case CWIPI_CELL_PRISMHO:
    _elt_type = FVMC_CELL_PRISM;
    break;
  case CWIPI_CELL_PYRAMHO:
    _elt_type = FVMC_CELL_PYRAM;
    break;
  default:
    bftc_error(__FILE__, __LINE__, 0,
               "_cwipi_ho_user_elt_set : unvailable element type\n");


  }

  FVMC_ho_basis_user_elt_set (_elt_type,
                              (fvmc_ho_basis_fct_t) element_basis);

  fvmc_ho_location_user_elt_set (_elt_type,
                                 (fvmc_ho_location_fct_t) location_in_element);
}


/*----------------------------------------------------------------------------
 *
 * Location completion.
 * This is a synchronization point with the coupled application
 *
 * parameters
 *   coupling_id          <-- Coupling identifier
 *----------------------------------------------------------------------------*/

void cwipi_locate (const char *coupling_name)
{
  cwipi::CouplingDataBase & couplingDataBase =
    cwipi::CouplingDataBase::getInstance();

  const std::string &coupling_name_str = coupling_name;

  cwipi::oldCoupling& coupling = couplingDataBase.getCoupling(coupling_name_str);

  coupling.locate();
}

/*----------------------------------------------------------------------------
 *
 * Update location
 *
 * parameters
 *   coupling_id          <-- Coupling identifier
 *----------------------------------------------------------------------------*/

void cwipi_update_location (const char *coupling_name)
{
  cwipi::CouplingDataBase & couplingDataBase =
    cwipi::CouplingDataBase::getInstance();

  const std::string &coupling_name_str = coupling_name;

  cwipi::oldCoupling& coupling = couplingDataBase.getCoupling(coupling_name_str);

  coupling.updateLocation();
}

/*----------------------------------------------------------------------------
 *
 * Get distant point Location
 *
 * parameters
 *   coupling_id          <-- Coupling identifier
 * return
 *   located points location
 *----------------------------------------------------------------------------*/

const int *cwipi_get_distant_location (const char *coupling_id)
{
  cwipi::CouplingDataBase & couplingDataBase =
    cwipi::CouplingDataBase::getInstance();

  const std::string &coupling_name_str = coupling_id;

  cwipi::oldCoupling& coupling = couplingDataBase.getCoupling(coupling_name_str);

  return coupling.getDistantLocation();
}

/*----------------------------------------------------------------------------
 *
 * Get barycentric coordinates index
 *
 * parameters
 *   coupling_id          <-- Coupling identifier
 * return
 *   barycentric coordinates index
 *----------------------------------------------------------------------------*/

const int *cwipi_get_distant_barycentric_coordinates_index (const char *coupling_id)
{
  cwipi::CouplingDataBase & couplingDataBase =
    cwipi::CouplingDataBase::getInstance();

  const std::string &coupling_name_str = coupling_id;

  cwipi::oldCoupling& coupling = couplingDataBase.getCoupling(coupling_name_str);

  return coupling.getDistantBarycentricCoordinatesIndex();
}

/*----------------------------------------------------------------------------
 *
 * Get barycentric coordinates
 *
 * parameters
 *   coupling_id          <-- Coupling identifier
 * return
 *   barycentric coordinates
 *----------------------------------------------------------------------------*/

const double *cwipi_get_distant_barycentric_coordinates (const char *coupling_id)
{
  cwipi::CouplingDataBase & couplingDataBase =
    cwipi::CouplingDataBase::getInstance();

  const std::string &coupling_name_str = coupling_id;

  cwipi::oldCoupling& coupling = couplingDataBase.getCoupling(coupling_name_str);

  return coupling.getDistantBarycentricCoordinates();
}

/*----------------------------------------------------------------------------
 *
 * Get number of located points
 *
 * parameters
 *   coupling_id          <-- Coupling identifier
 *
 * return
 *                        --> Number of located points
 *
 *----------------------------------------------------------------------------*/

int cwipi_get_n_located_points(const char *coupling_id)
{
  cwipi::CouplingDataBase & couplingDataBase =
    cwipi::CouplingDataBase::getInstance();

  const std::string &coupling_name_str = coupling_id;

  cwipi::oldCoupling& coupling = couplingDataBase.getCoupling(coupling_name_str);

  return coupling.getNLocatedPoint();
}

/*----------------------------------------------------------------------------
 *
 * Get number of not located points
 *
 * parameters
 *   coupling_id          <-- Coupling identifier
 *
 * return
 *                        --> Number of not located points
 *
 *----------------------------------------------------------------------------*/

int cwipi_get_n_not_located_points(const char *coupling_id)
{
  cwipi::CouplingDataBase & couplingDataBase =
    cwipi::CouplingDataBase::getInstance();

  const std::string &coupling_name_str = coupling_id;

  cwipi::oldCoupling& coupling = couplingDataBase.getCoupling(coupling_name_str);

  return coupling.getNNotlocatedPoint();
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
 *   stride               <-- Number of interlaced field
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

cwipi_exchange_status_t cwipi_exchange
(const char                *coupling_name,
 const char                *exchange_name,
 const int                  stride,
 const int                  time_step,
 const double               time_value,
 const char                *sending_field_name,
 const double              *sending_field,
 const char                *receiving_field_name,
 double                    *receiving_field,
 int                       *nNotLocatedPoints)

{
  cwipi::CouplingDataBase & couplingDataBase =
    cwipi::CouplingDataBase::getInstance();

  const std::string &coupling_name_str = coupling_name;

  cwipi::oldCoupling& coupling = couplingDataBase.getCoupling(coupling_name_str);

  cwipi_exchange_status_t status;

  status = coupling.exchange(exchange_name,
                             stride,
                             time_step,
                             time_value,
                             sending_field_name,
                             sending_field,
                             (char *) receiving_field_name,
                             receiving_field,
                             NULL);

  *nNotLocatedPoints = coupling.getNNotlocatedPoint();

  return status;
}

/*----------------------------------------------------------------------------
 *
 * Send interpolated data to the coupled application.
 * Non blocking comunication.
 *
 * parameters
 *   coupling_id          <-- Coupling identifier
 *   exchange_name        <-- Exchange name
 *   stride               <-- Number of interlaced field
 *   time_step            <-- Time step  (only for visualization)
 *   time_value           <-- Time value (only for visualization)
 *   sending_field_name   <-- Sending field name
 *   sending_field        <-- Sending field
 *   request              --> Request
 *
 *----------------------------------------------------------------------------*/

void cwipi_issend
(const char                *coupling_name,
 const char                *exchange_name,
 const int                 tag,
 const int                 stride,
 const int                 time_step,
 const double              time_value,
 const char                *sending_field_name,
 const double              *sending_field,
 int                       *request)
{
  cwipi::CouplingDataBase & couplingDataBase =
    cwipi::CouplingDataBase::getInstance();

  const std::string &coupling_name_str = coupling_name;

  cwipi::oldCoupling& coupling = couplingDataBase.getCoupling(coupling_name_str);

  coupling.issend(exchange_name,
                  tag,
                  stride,
                  time_step,
                  time_value,
                  sending_field_name,
                  sending_field,
                  NULL,
                  request);

}

/*----------------------------------------------------------------------------
 *
 * Receive interpolated data from the coupled application.
 * Non blocking comunication. receiving_field is fully updated after
 * cwipi_wait_irecv calling
 *
 * parameters
 *   coupling_id          <-- Coupling identifier
 *   exchange_name        <-- Exchange name
 *   stride               <-- Number of interlaced field
 *   time_step            <-- Time step  (only for visualization)
 *   time_value           <-- Time value (only for visualization)
 *   receiving_field_name <-- Receiving field name
 *   receiving_field      <-- Receiving field
 *   request              --> Request
 *
 *----------------------------------------------------------------------------*/

void cwipi_irecv
(const char                *coupling_name,
 const char                *exchange_name,
 const int                 tag,
 const int                 stride,
 const int                 time_step,
 const double              time_value,
 const char                *receiving_field_name,
 double              *receiving_field,
 int                       *request)
{
  cwipi::CouplingDataBase & couplingDataBase =
    cwipi::CouplingDataBase::getInstance();

  const std::string &coupling_name_str = coupling_name;

  cwipi::oldCoupling& coupling = couplingDataBase.getCoupling(coupling_name_str);

  coupling.irecv(exchange_name,
                 tag,
                 stride,
                 time_step,
                 time_value,
                 receiving_field_name,
                 receiving_field,
                 request);
}

/*----------------------------------------------------------------------------
 *
 * Get located point distance to exchange area
 *
 * parameters
 *   coupling_id          <-- Coupling identifier
 *
 *----------------------------------------------------------------------------*/

const float *cwipi_distance_located_pts_get(const char  *coupling_name)
{
  cwipi::CouplingDataBase & couplingDataBase =
    cwipi::CouplingDataBase::getInstance();

  const std::string &coupling_name_str = coupling_name;

  cwipi::oldCoupling& coupling = couplingDataBase.getCoupling(coupling_name_str);

  return coupling.distance();
}


/*----------------------------------------------------------------------------
 *
 * Wait for cwipi_issend.
 *
 * parameters
 *   coupling_id          <-- Coupling identifier
 *   request              <-- Request
 *
 *----------------------------------------------------------------------------*/

void cwipi_wait_issend(const char  *coupling_name,
                       int          request)
{
  cwipi::CouplingDataBase & couplingDataBase =
    cwipi::CouplingDataBase::getInstance();

  const std::string &coupling_name_str = coupling_name;

  cwipi::oldCoupling& coupling = couplingDataBase.getCoupling(coupling_name_str);

  coupling.waitIssend(request);
}

/*----------------------------------------------------------------------------
 *
 * Wait for cwipi_irecv.
 *
 * parameters
 *   coupling_id          <-- Coupling identifier
 *   request              <-- Request
 *
 *----------------------------------------------------------------------------*/

void cwipi_wait_irecv(const char  *coupling_name,
                     int          request)
{
  cwipi::CouplingDataBase & couplingDataBase =
    cwipi::CouplingDataBase::getInstance();

  const std::string &coupling_name_str = coupling_name;

  cwipi::oldCoupling& coupling = couplingDataBase.getCoupling(coupling_name_str);

  coupling.waitIrecv(request);
}

/*----------------------------------------------------------------------------
 *
 * Define the interpolation function
 *
 * parameters
 *   coupling_id          <-- Coupling identifier
 *   fct                  <-- Interpolation function
 *
 *----------------------------------------------------------------------------*/

void cwipi_set_interpolation_function
(const char *coupling_name,
 cwipi_interpolation_fct_t fct)
{
  cwipi::CouplingDataBase & couplingDataBase =
    cwipi::CouplingDataBase::getInstance();

  const std::string &coupling_name_str = coupling_name;

  cwipi::oldCoupling& coupling = couplingDataBase.getCoupling(coupling_name_str);

  coupling.set_interpolation_function(fct);
}

/*----------------------------------------------------------------------------
 *
 * Define the interpolation function for high order
 *
 * parameters
 *   coupling_id          <-- Coupling identifier
 *   fct                  <-- Interpolation function
 *
 *----------------------------------------------------------------------------*/

void cwipi_ho_set_interpolation_function
(const char *coupling_name,
 cwipi_user_interp_ho_fct_t fct)
{
  cwipi::CouplingDataBase & couplingDataBase =
    cwipi::CouplingDataBase::getInstance();

  const std::string &coupling_name_str = coupling_name;

  cwipi::oldCoupling& coupling = couplingDataBase.getCoupling(coupling_name_str);

  coupling.set_ho_interpolation_function(fct);
}


/*----------------------------------------------------------------------------
 *
 * Define the interpolation function
 *
 * parameters
 *   coupling_id          <-- Coupling identifier
 *   fct                  <-- Interpolation function
 *
 *----------------------------------------------------------------------------*/

void cwipi_set_interpolation_function_f
(const char *coupling_name,
 void * fct)
{
  cwipi::CouplingDataBase & couplingDataBase =
    cwipi::CouplingDataBase::getInstance();

  const std::string &coupling_name_str = coupling_name;

  cwipi::oldCoupling& coupling = couplingDataBase.getCoupling(coupling_name_str);

  coupling.set_interpolation_function_f(fct);
}


/*----------------------------------------------------------------------------
 *
 * Delete a coupling
 *
 * parameters
 *   coupling_id          <-- Coupling identifier
 *
 *----------------------------------------------------------------------------*/

void cwipi_delete_coupling(const char *coupling_name)
{
  cwipi::CouplingDataBase & couplingDataBase =
    cwipi::CouplingDataBase::getInstance();

  const std::string &coupling_name_str = coupling_name;

  couplingDataBase.deleteCoupling(coupling_name_str);
}

/*----------------------------------------------------------------------------
 *
 * Finalize cwipi. This is a synchronization point between all applications
 *
 *----------------------------------------------------------------------------*/

void cwipi_finalize(void)
{
  cwipi::CouplingDataBase & couplingDataBase =
    cwipi::CouplingDataBase::getInstance();

  cwipi::ApplicationPropertiesDataBase & properties =
    cwipi::ApplicationPropertiesDataBase::getInstance();

  const MPI_Comm globalComm = properties.getGlobalComm();

  fvmc_ho_location_user_elts_unset ();
  FVMC_ho_basis_user_elts_unset ();

  bftc_printf("Finalize cwipi\n");
  couplingDataBase.kill();
  properties.kill();

  int flag = 0;
  MPI_Initialized(&flag);

  if (flag != 0) {
    bftc_printf_flush();
    MPI_Barrier(globalComm);
  }

  FVMC_ho_basis_free ();

}

/*----------------------------------------------------------------------------
 *
 * Dump application properties
 *
 *----------------------------------------------------------------------------*/

void cwipi_dump_application_properties(void)
{
  cwipi::ApplicationPropertiesDataBase & properties =
    cwipi::ApplicationPropertiesDataBase::getInstance();
  properties.dump();
}

/*----------------------------------------------------------------------------
 *
 * Get not located points
 *
 * parameters
 *   coupling_id          <-- Coupling identifier
 *
 * return
 *   notLocatedPoints     <-- Not located points
 *
 *----------------------------------------------------------------------------*/

const int * cwipi_get_not_located_points(const char *coupling_id)
{
  cwipi::CouplingDataBase & couplingDataBase =
    cwipi::CouplingDataBase::getInstance();

  const std::string &coupling_name_str = coupling_id;

  cwipi::oldCoupling& coupling = couplingDataBase.getCoupling(coupling_name_str);

  return coupling.getNotlocatedPoint();
}

/*----------------------------------------------------------------------------
 *
 * Get not located points
 *
 * parameters
 *   coupling_id          <-- Coupling identifier
 *
 * return
 *   locatedPoints        <-- Located points
 *
 *----------------------------------------------------------------------------*/

const int * cwipi_get_located_points(const char *coupling_id)
{
  cwipi::CouplingDataBase & couplingDataBase =
    cwipi::CouplingDataBase::getInstance();

  const std::string &coupling_name_str = coupling_id;

  cwipi::oldCoupling& coupling = couplingDataBase.getCoupling(coupling_name_str);

  return coupling.getLocatedPoint();
}

/*----------------------------------------------------------------------------
 *
 * Get number of located distant point
 *
 * parameters
 *   coupling_id          <-- Coupling identifier
 *
 * return
 *                        --> Number of located distant points
 *
 *----------------------------------------------------------------------------*/

int cwipi_get_n_distant_points(const char *coupling_id)
{
  cwipi::CouplingDataBase & couplingDataBase =
    cwipi::CouplingDataBase::getInstance();

  const std::string &coupling_name_str = coupling_id;

  cwipi::oldCoupling& coupling = couplingDataBase.getCoupling(coupling_name_str);

  return coupling.getNDistantPoint();
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

int cwipi_get_n_distant_ranks(const char *coupling_id)
{
  cwipi::CouplingDataBase & couplingDataBase =
    cwipi::CouplingDataBase::getInstance();

  const std::string &coupling_name_str = coupling_id;

  cwipi::oldCoupling& coupling = couplingDataBase.getCoupling(coupling_name_str);

  return coupling.getNDistantRank();
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

const int *cwipi_get_distant_distribution(const char *coupling_id)
{
  cwipi::CouplingDataBase & couplingDataBase =
    cwipi::CouplingDataBase::getInstance();

  const std::string &coupling_name_str = coupling_id;

  cwipi::oldCoupling& coupling = couplingDataBase.getCoupling(coupling_name_str);

  return coupling.getDistantDistribution();
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

const int *cwipi_get_located_points_distribution(const char *coupling_id)
{
  cwipi::CouplingDataBase & couplingDataBase =
    cwipi::CouplingDataBase::getInstance();

  const std::string &coupling_name_str = coupling_id;

  cwipi::oldCoupling& coupling = couplingDataBase.getCoupling(coupling_name_str);

  return coupling.getLocatedPointsDistribution();
}


/*----------------------------------------------------------------------------
 *
 * Get distant point coordinates
 *
 * parameters
 *   coupling_id          <-- Coupling identifier
 * return
 *   coordinates
 *----------------------------------------------------------------------------*/

const double *cwipi_get_distant_coordinates (const char *coupling_id)
{
  cwipi::CouplingDataBase & couplingDataBase =
    cwipi::CouplingDataBase::getInstance();

  const std::string &coupling_name_str = coupling_id;

  cwipi::oldCoupling& coupling = couplingDataBase.getCoupling(coupling_name_str);

  return coupling.getDistantPointCoordinates();
}

/*----------------------------------------------------------------------------
 *
 * Get distance to distant location element
 *
 * parameters
 *   coupling_id          <-- Coupling identifier
 * return
 *   distance
 *----------------------------------------------------------------------------*/

const float *cwipi_get_distant_distance (const char *coupling_id)
{
  cwipi::CouplingDataBase & couplingDataBase =
    cwipi::CouplingDataBase::getInstance();

  const std::string &coupling_name_str = coupling_id;

  cwipi::oldCoupling& coupling = couplingDataBase.getCoupling(coupling_name_str);

  return coupling.distance();
}



/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
}
#endif /* __cplusplus */
