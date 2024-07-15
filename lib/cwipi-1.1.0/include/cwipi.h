/*
 * \file
 */

#ifndef __CWIPI_H__
#define __CWIPI_H__
/*
  This file is part of the CWIPI library. 

  Copyright (C) 2011-2017  ONERA

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

#include <mpi.h>

#include <stdio.h>

/*=============================================================================
 * Macro definitions
 *============================================================================*/

#if !defined (__hpux) && !defined (_AIX)
#define PROCF(x, y) x##_
#else
#define PROCF(x, y) x
#endif

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force back Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */

/*============================================================================
 * Type definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * MPI ranks used for the coupling
 *----------------------------------------------------------------------------*/

typedef enum {

  CWIPI_COUPLING_PARALLEL_WITH_PARTITIONING,
  CWIPI_COUPLING_PARALLEL_WITHOUT_PARTITIONING,
  CWIPI_COUPLING_SEQUENTIAL

} cwipi_coupling_type_t;

/*----------------------------------------------------------------------------
 * MPI ranks used for the coupling
 *----------------------------------------------------------------------------*/

typedef enum {

  CWIPI_BASIC_INFO,
  CWIPI_DISTANT_MESH_INFO

} cwipi_located_point_info_t;

/*----------------------------------------------------------------------------
 * Mesh type
 *----------------------------------------------------------------------------*/

typedef enum {

  CWIPI_STATIC_MESH,
  CWIPI_MOBILE_MESH,
  CWIPI_CYCLIC_MESH

} cwipi_mesh_type_t;

/*----------------------------------------------------------------------------
 * Solver type
 *----------------------------------------------------------------------------*/

typedef enum {

  CWIPI_SOLVER_CELL_CENTER,
  CWIPI_SOLVER_CELL_VERTEX

} cwipi_solver_type_t;

/*----------------------------------------------------------------------------
 * Coupling type
 *----------------------------------------------------------------------------*/

typedef enum {

  CWIPI_FIELD_TYPE_FLOAT,
  CWIPI_FIELD_TYPE_DOUBLE

} cwipi_field_type_t;

/*----------------------------------------------------------------------------
 * Coupling interpolation type
 *----------------------------------------------------------------------------*/

typedef enum {

  CWIPI_INTERPOLATION_DEFAULT,
  CWIPI_INTERPOLATION_USER

} cwipi_interpolation_t;

/*----------------------------------------------------------------------------
 * Coupling exchange status
 *----------------------------------------------------------------------------*/

typedef enum {

  CWIPI_EXCHANGE_OK,
  CWIPI_EXCHANGE_BAD_RECEIVING

} cwipi_exchange_status_t;

/*----------------------------------------------------------------------------
 * Element type
 *----------------------------------------------------------------------------*/

typedef enum {

  CWIPI_NODE,
  CWIPI_EDGE2,
  CWIPI_EDGEHO,
  CWIPI_FACE_TRIA3,
  CWIPI_FACE_TRIAHO,
  CWIPI_FACE_QUAD4,
  CWIPI_FACE_QUADHO,  
  CWIPI_FACE_POLY,
  CWIPI_CELL_TETRA4,
  CWIPI_CELL_TETRAHO,
  CWIPI_CELL_HEXA8,
  CWIPI_CELL_HEXAHO,
  CWIPI_CELL_PRISM6,
  CWIPI_CELL_PRISMHO,
  CWIPI_CELL_PYRAM5,
  CWIPI_CELL_PYRAMHO,
  CWIPI_CELL_POLY

} cwipi_element_t;

/*----------------------------------------------------------------------------
 * Function pointer to define an user interpolation method (callback)
 *
 * parameters:
 * ----------
 *
 * entities_dim                              <-- entities dimension of
 *                                               the local mesh (1, 2 or 3)
 * n_local_vertex                            <-- local mesh vertices number
 * n_local_element                           <-- local mesh elements number
 *                                               (without polyhedra)
 * n_local_polyhedra                         <-- local mesh elements number
 * n_distant_point                           <-- located distant point number
 * local_coordinates                         <-- local mesh vertex coordinates
 * local_parent_elt_num                      <-- pointer to parent element
 *                                               (or NULL if sorted elements)
 * local_connectivity_index                  <-- element -> vertices index
 *                                               (O to n-1)
 *                                               size:n_local_elements + 1
 * local_connectivity                        <-- element -> vertex connectivity
 *                                                       of the local mesh
 *                               size:local_connectivity_index[n_local_elements]
 * local_polyhedra_face_index                <-- polyhedra volume -> faces index
 * local_polyhedra_cell_to_face_connectivity <-- polyhedra -> face connectivity
 * local_polyhedra_face_connectivity_index   <-- polyhedra faces
 *                                               face -> vertices index
 * local_polyhedra_face_connectivity         <-- polyhedra
 *                                               face -> vertex connectivity
 * distant_points_coordinates                <-- distant point coordinates
 * distant_points_location                   <-- distant point location
 * distant_points_barycentric_coordinates_index
 *                                           <-- element -> barycentric coordinates
 *                                                (0 to n-1)
 *                                               size: n_distant_point + 1
 * distant_points_barycentric_coordinates    <-- distant point barycentric coordinates
 *                                             size: distant_points_barycentric_coordinates_index[n_distant_point]
 * stride                                    <-- interlaced field number
 * local_field                               <-- local field
 * distant_field                             --> distant field
 *
 *----------------------------------------------------------------------------*/

typedef void (*cwipi_interpolation_fct_t)
  (const int entities_dim,
   const int n_local_vertex,
   const int n_local_element,
   const int n_local_polhyedra,
   const int n_distant_point,
   const double local_coordinates[],
   const int local_connectivity_index[],
   const int local_connectivity[],
   const int local_polyhedra_face_index[],
   const int local_polyhedra_cell_to_face_connectivity[],
   const int local_polyhedra_face_connectivity_index[],
   const int local_polyhedra_face_connectivity[],
   const double distant_points_coordinates[],
   const int distant_points_location[],
   const float distant_points_distance[],
   const int distant_points_barycentric_coordinates_index[],
   const double distant_points_barycentric_coordinates[],
   const int stride,
   const cwipi_solver_type_t  solver_type,
   const void *local_field,
   void *distant_field
   );



/*----------------------------------------------------------------------------
 * Function pointer to define an user interpolation method (callback)
 *
 * parameters:
 * ----------
 *
 * entities_dim                              <-- entities dimension of
 *                                               the local mesh (1, 2 or 3)
 * order                                     <-- Mesh order
 * n_local_vertex                            <-- local mesh vertices number
 * n_local_element                           <-- local mesh elements number 
 *                                               (without polyhedra)
 * n_local_polyhedra                         <-- local mesh elements number
 * n_distant_point                           <-- located distant point number
 * local_coordinates                         <-- local mesh vertex coordinates
 * local_parent_elt_num                      <-- pointer to parent element
 *                                               (or NULL if sorted elements)
 * local_connectivity_index                  <-- element -> vertices index
 *                                               (O to n-1)
 *                                               size:n_local_elements + 1
 * local_connectivity                        <-- element -> vertex connectivity
 *                                                       of the local mesh
 *                                               size:local_connectivity_index[n_local_elements]
 * local_polyhedra_face_index                <-- polyhedra volume -> faces index
 * local_polyhedra_cell_to_face_connectivity <-- polyhedra -> face connectivity
 * local_polyhedra_face_connectivity_index   <-- polyhedra faces
 *                                               face -> vertices index
 * local_polyhedra_face_connectivity         <-- polyhedra
 *                                               face -> vertex connectivity
 * distant_points_coordinates                <-- distant point coordinates
 * distant_points_location                   <-- distant point location
 * distant_points_weights_index
 *                                           <-- element -> weights
 *                                                (0 to n-1)
 *                                               size: n_distant_point + 1
 * distant_points_weights                    <-- distant point weights
 *                                             size: distant_points_barycentric_coordinates_index[n_distant_point]
 * distant_points_uvw                        <-- parametric coordinates of distant points (size = uvw_size * n_distant_point)
 * stride                                    <-- interlaced field number
 * local_field                               <-- local field
 * distant_field                             --> distant field
 *
 *----------------------------------------------------------------------------*/

typedef void (*cwipi_user_interp_ho_fct_t)
  (const int entities_dim,
   const int order,
   const int n_local_vertex,
   const int n_local_element,
   const int n_local_polhyedra,
   const int n_distant_point,
   const double local_coordinates[],
   const int local_connectivity_index[],
   const int local_connectivity[],
   const int local_polyhedra_face_index[],
   const int local_polyhedra_cell_to_face_connectivity[],
   const int local_polyhedra_face_connectivity_index[],
   const int local_polyhedra_face_connectivity[],
   const double distant_points_coordinates[],
   const int distant_points_location[],
   const float distant_points_distance[],
   const int distant_points_weights_index[],
   const double distant_points_weights[],
   const double distant_points_uvw[],
   const int stride,
   const cwipi_solver_type_t solver_type,
   const void *local_field,
   void *distant_field
   );

/*----------------------------------------------------------------------------
 * 
 * Callback to define location in a high order element
 * 
 * parameters:
 *   order             <-- element order
 *   n_nodes           <-- number of nodes of the element
 *   nodes_coords      <-- nodes coordinates
 *   point_coords      <-- point to locate coordinates
 *   projected_coords  --> projected point coordinates (if point is outside) 
 *   projected_uvw     --> parametric coordinates of the projected point
 * 
 * return: 
 *   distance to the cell (distance <= 0 if point is inside)
 *
 *----------------------------------------------------------------------------*/

typedef double (*cwipi_ho_location_fct_t)
(const int entities_dim,
 const int order,
 const int n_nodes,
 const double *nodes_coords,
 const double *point_coords,
 double *projected_coords,
 double *projected_uvw);


/*----------------------------------------------------------------------------
 * 
 * Callback to define the basis functions of an high order 
 * element
 * 
 * parameters:
 *   order             <-- element order
 *   n_nodes           <-- number of nodes of the element
 *   n_pts             <-- number of points
 *   uvv               <-- Parametric coordinates of points
 *   weights           --> Interpolation weights (size = n_nodes * n_pts)
 * 
 *----------------------------------------------------------------------------*/

typedef void (*cwipi_ho_basis_fct_t)
(const int entities_dim,
 const int order,
 const int n_nodes,
 const int n_pts,
 const double *uvw,
 double *weights);


/*----------------------------------------------------------------------------
 * 
 * Callback to define parametric coordinates of the element nodes
 * 
 * parameters:
 *   order             <-- element order
 *   n_nodes           <-- number of nodes of the element
 *   xsi_uvv           --> Parametric coordinates of a the element nodes
 *                         (size = dim of element * n_nodes)
 * 
 *----------------------------------------------------------------------------*/

typedef void (*cwipi_ho_xsi_fct_t)
(const int entities_dim,
 const int order,
 const int n_nodes,
 double *xsi_uvw);


/*=============================================================================
 * Static global variables
 *============================================================================*/

/*=============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------
 *
 * Initialize the cwipi library.
 * Redirect outputs in a file (Standard output with output_listing = NULL or
 * output_logical_unit = -1)
 * Create the current communicator application from 'common_comm'.
 *
 * parameters:
 *   common_comm       <-- Common MPI communicator
 *   application_name  <-- Current application name
 *   application_comm  --> Internal MPI communicator for the current
 *                         application
 *
 * It is a synchronization point between all applications
 *----------------------------------------------------------------------------*/

void cwipi_init
(const MPI_Comm                           common_comm,
 const char                               *application_name,
 MPI_Comm                                 *application_comm);

/*----------------------------------------------------------------------------
 *
 * Set up the file used for the output listing
 *
 * parameters:
 *   output_listing      <-- Output listing file (C function)
 *----------------------------------------------------------------------------*/

void PROCF (cwipi_set_output_logical_unit, CWIPI_SET_OUTPUT_LOGICAL_UNIT) (int *iunit);

void cwipi_set_output_listing(FILE *output_listing);

/*----------------------------------------------------------------------------
 *
 * Add a integer control parameter
 *
 * parameters
 *    name           <-- parameter name
 *    initial_value  <-- initial value
 *
 *----------------------------------------------------------------------------*/

void cwipi_add_local_int_control_parameter(const char *name, int initial_value);

/*----------------------------------------------------------------------------
 *
 * Add a double control parameter
 *
 * parameters
 *    name           <-- parameter name
 *    initial_value  <-- initial value
 *
 *----------------------------------------------------------------------------*/

void cwipi_add_local_double_control_parameter(const char *name, double initial_value);

/*----------------------------------------------------------------------------
 *
 * Add a string control parameter
 *
 * parameters
 *    name           <-- parameter name
 *    initial_value  <-- initial value
 *
 *----------------------------------------------------------------------------*/

void cwipi_add_local_string_control_parameter(const char *name, const char *initial_value);

/*----------------------------------------------------------------------------
 *
 * Set a integer control parameter
 *
 * parameters
 *    name           <-- parameter name
 *    value          <-- value
 *
 *----------------------------------------------------------------------------*/

void cwipi_set_local_int_control_parameter(const char *name, int value);

/*----------------------------------------------------------------------------
 *
 * Set a double control parameter
 *
 * parameters
 *    name           <-- parameter name
 *    value          <-- value
 *
 *----------------------------------------------------------------------------*/

void cwipi_set_local_double_control_parameter(const char *name, double value);

/*----------------------------------------------------------------------------
 *
 * Set a string control parameter
 *
 * parameters
 *    name           <-- parameter name
 *    value          <-- value
 *
 *----------------------------------------------------------------------------*/

void cwipi_set_local_string_control_parameter(const char *name, const char *value);

/*----------------------------------------------------------------------------
 *
 * Get a integer control parameter of the current application
 *
 * parameters
 *    name           <-- parameter name
 *
 *----------------------------------------------------------------------------*/

int cwipi_get_local_int_control_parameter(const char *name);

/*----------------------------------------------------------------------------
 *
 * Get a double control parameter of the current application
 *
 * parameters
 *    name           <-- parameter name
 *
 *----------------------------------------------------------------------------*/

double cwipi_get_local_double_control_parameter(const char *name);

/*----------------------------------------------------------------------------
 *
 * Get a string control parameter of the current application
 *
 * parameters
 *    name           <-- parameter name
 *
 *----------------------------------------------------------------------------*/

const char* cwipi_get_local_string_control_parameter(const char *name);

/*----------------------------------------------------------------------------
 *
 * Delete a current application int parameter
 *
 * parameters
 *    name           <-- parameter name
 *
 *----------------------------------------------------------------------------*/

void cwipi_delete_local_int_control_parameter(const char *name);

/*----------------------------------------------------------------------------
 *
 * Delete a current application double parameter
 *
 * parameters
 *    name           <-- parameter name
 *
 *----------------------------------------------------------------------------*/

void cwipi_delete_local_double_control_parameter(const char *name);

/*----------------------------------------------------------------------------
 *
 * Delete a current application string parameter
 *
 * parameters
 *    name           <-- parameter name
 *
 *----------------------------------------------------------------------------*/

void cwipi_delete_local_string_control_parameter(const char *name);

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
 const char *name);

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
 const char *name);

/*----------------------------------------------------------------------------
 *
 * Get a string control parameter of a other application
 *
 * parameters
 *    application_name    <-- application name
 *    name                <-- parameter name
 *
 *----------------------------------------------------------------------------*/

const char* cwipi_get_distant_string_control_parameter
(const char *application_name,
 const char *name);

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
 const char *name);

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
 const char *name);

/*----------------------------------------------------------------------------
 *
 * Get a string control parameter of a other application
 *
 * parameters
 *    application_name    <-- application name
 *    name                <-- parameter name
 *
 *----------------------------------------------------------------------------*/

const char* cwipi_get_distant_string_control_parameter
(const char *application_name,
 const char *name);

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
 const char *name);

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
 const char *name);

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
 const char *name);

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
(const char *application_name);

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
(const char *application_name);

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
(const char *application_name);

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
(const char *application_name);

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
(const char *application_name);

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
(const char *application_name);

/*----------------------------------------------------------------------------
 *
 * Synchronize local control parameters with an other application.
 *  It is a synchronization point with this second application
 *
 * parameters
 *    application_name    <-- application name
 *
 *----------------------------------------------------------------------------*/

void cwipi_synchronize_control_parameter(const char *application_name);

/*----------------------------------------------------------------------------
 *
 * Dump application properties
 *
 *----------------------------------------------------------------------------*/

void cwipi_dump_application_properties(void);

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
 *   output_format_option    <-- Output options "opt1, opt2,
 *                             text             output text files
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
 *   nblocations             <-- Number of possible localisations with
 *                               CWIPI_CYCLIC_MESH, optional
 *
 *----------------------------------------------------------------------------*/

void cwipi_create_coupling
( const char  *coupling_name,
  const cwipi_coupling_type_t coupling_type,
  const char  *coupled_application,
  const int    entitiesDim,
  const double tolerance,
  const cwipi_mesh_type_t mesh_type,
  const cwipi_solver_type_t solver_type,
  const int    output_frequency,
  const char  *output_format,
  const char  *output_format_option,
  ...);

//TODO : ajouter un deuxieme coupling permettant d'eliminer certains procs
// du couplage

/* void cwipi_create_coupling_select_proc */
/* ( const char  *coupling_name, */
/*   const cwipi_coupling_type_t coupling_type, */
/*   const char  *coupled_application, */
/*   const int    entitiesDim, */
/*   const double tolerance, */
/*   const cwipi_mesh_type_t mesh_type, */
/*   const cwipi_solver_type_t solver_type, */
/*   const int    output_frequency, */
/*   const char  *output_format, */
/*   const char  *output_format_option, */
/*   int   isCoupled (TRUE/FALSE), */
/*   ...); */

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
);


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
);

/*----------------------------------------------------------------------------
 *
 * Set the index location for multiple location with CWIPI_CYCLIC_MESH
 *
 * parameters:
 *   coupling_name           <-- Coupling identifier
 *   index                   <-- location index
 *----------------------------------------------------------------------------*/

void cwipi_set_location_index
( const char  *coupling_name,
  const int index);

/*----------------------------------------------------------------------------
 *
 * save/load  location 
 *
 * parameters:
 *   coupling_name           <-- Coupling identifier
 *----------------------------------------------------------------------------*/

void cwipi_load_location(const char *coupling_name);
void cwipi_save_location(const char *coupling_name);


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
                               const char *mode);
/*----------------------------------------------------------------------------
 *
 * cwipi_close_location_file
 *
 * parameters
 *   coupling_id          <-- Coupling identifier
 *----------------------------------------------------------------------------*/

void cwipi_close_location_file (const char *coupling_name);

/*----------------------------------------------------------------------------
 *
 * Set points to locate. This function must be called if the points to locate
 * do not correspond to :
 *        - vertices for CELL_VERTEX nature
 *        - cell center for CELL_CENTER nature
 *
 * parameters:
 *   coupling_id        <-- coupling identifier
 *   n_points           <-- number of points to locate
 *   coordinates        <-- coordinates of points to locate (enterlaced)
 *
 *----------------------------------------------------------------------------*/

void cwipi_set_points_to_locate
(const char  *coupling_id,
 const int    n_points,
 double       coordinate[]);

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
 *   coordinates        <-- vertex interlaced coordinates (always 3 coordinates, even in 2D)
 *   connectivity_index <-> element -> vertices index (O to n-1)
 *                          size: n_elements + 1
 *                          (out : ordered connectivity_index)
 *   connectivity       <-> element -> vertex connectivity (1 to n)
 *                          size: connectivity_index[n_elements]
 *                          (out : ordered connectivity)
 *
 *----------------------------------------------------------------------------*/

void cwipi_define_mesh(const char *coupling_id,
                       const int n_vertex,
                       const int n_element,
                       double coordinates[],
                       int connectivity_index[],
                       int connectivity[]);

void cwipi_shared_fvmc_nodal(const char *coupling_name,
                            void * fvmc_nodal);


/*----------------------------------------------------------------------------
 *
 * Define a high order mesh interface for the current coupling. 
 *
 *    1D : edges (not implemented yet)
 *    2D : triangles, quadrangles
 *    3D : tetrahedra, pyramids, prism, hexaedra (not implemented yet)
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
 *   connectivity       <-> element -> vertex connectivity (1 to n)
 *                          size: connectivity_index[n_elements]
 *
 *----------------------------------------------------------------------------*/

void cwipi_ho_define_mesh(const char *coupling_id,
                          const int n_vertex,
                          const int n_element,
                          const int order,
                          double coordinates[],
                          int connectivity_index[],
                          int connectivity[]);



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
                           const char *value);

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
                                     const int *ijk_grid);

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
                                         const double *coords);

/*----------------------------------------------------------------------------
 * 
 * Set a user element
 * 
 * parameters:
 *   elt_type            <-- Element type
 *   element_basis       <-- Element basis function
 *   xsi_coordinates     <-- Parametric coordinates of nodes function
 *   location_in_element <-- Element location in element function
 *
 *----------------------------------------------------------------------------*/

void
cwipi_ho_user_elt_set (cwipi_element_t elt_type,
                       cwipi_ho_basis_fct_t element_basis,
                       cwipi_ho_location_fct_t location_in_element);

/*----------------------------------------------------------------------------
 *
 * Add polyhedra to the mesh
 *
 * parameters:
 *   coupling_id                  <-- Coupling identifier
 *   n_elements                   <-- Polyhedra number to add
 *   face_index                   <-- Face index (0 to n-1)
 *                                    size : n_elements + 1
 *   cell_to_face_connectivity    <-- Polyhedra -> face (1 to n)
 *                                    size : face_index[n_elements]
 *   n_faces                      <-- Faces number
 *   face_connectivity_index      <-- Face connectivity index (0 to n-1)
 *                                    size : n_faces + 1
 *   face_connectivity            <-- Face connectivity (1 to n)
 *                                    size : face_connectivity_index[n_faces]
 *
 *----------------------------------------------------------------------------*/

void cwipi_add_polyhedra(const char *coupling_id,
                         const int n_element,
                         int face_index[],
                         int cell_to_face_connectivity[],
                         const int n_faces,
                         int face_connectivity_index[],
                         int face_connectivity[]);

/*----------------------------------------------------------------------------
 *
 * Distant points location
 * (synchronization point with the coupled application)
 *
 * parameters
 *   coupling_id          <-- Coupling identifier
 *----------------------------------------------------------------------------*/

void cwipi_locate (const char *coupling_id);

/*----------------------------------------------------------------------------
 *
 * Update location
 *
 * parameters
 *   coupling_id          <-- Coupling identifier
 *----------------------------------------------------------------------------*/

void cwipi_update_location (const char *coupling_id);

/*----------------------------------------------------------------------------
 *
 * Get distant point Location
 *
 * parameters
 *   coupling_id          <-- Coupling identifier
 * return
 *   distant point location
 *----------------------------------------------------------------------------*/

const int *cwipi_get_distant_location (const char *coupling_id);

/*----------------------------------------------------------------------------
 *
 * Get distance to distant location element
 *
 * parameters
 *   coupling_id          <-- Coupling identifier
 * return
 *   distance
 *----------------------------------------------------------------------------*/

const float *cwipi_get_distant_distance (const char *coupling_id);

/*----------------------------------------------------------------------------
 *
 * Get barycentric coordinates index
 *
 * parameters
 *   coupling_id          <-- Coupling identifier
 * return
 *   barycentric coordinates index
 *----------------------------------------------------------------------------*/

const int *cwipi_get_distant_barycentric_coordinates_index (const char *coupling_id);

/*----------------------------------------------------------------------------
 *
 * Get barycentric coordinates index
 *
 * parameters
 *   coupling_id          <-- Coupling identifier
 * return
 *   barycentric coordinates
 *----------------------------------------------------------------------------*/

const double *cwipi_get_distant_barycentric_coordinates (const char *coupling_id);

/*----------------------------------------------------------------------------
 *
 * Get distant point coordinates
 *
 * parameters
 *   coupling_id          <-- Coupling identifier
 * return
 *   coordinates
 *----------------------------------------------------------------------------*/

const double *cwipi_get_distant_coordinates (const char *coupling_id);

/*----------------------------------------------------------------------------
 *
 * Exchange data with the coupled application.
 * It is a synchronization point with the coupled application
 *
 * parameters
 *   coupling_id          <-- Coupling identifier
 *   exchange_name        <-- Exchange name
 *   exchange_type        <-- Exchange type (not implemented yet)
 *   stride               <-- Number of interlaced fields
 *   time_step            <-- Time step  (only for visualization)
 *   time_value           <-- Time value (only for visualization)
 *   sending_field_name   <-- Sending field name
 *   sending_field        <-- Sending field (NULL -> no sending)
 *   receiving_field_name <-- Receiving field name
 *   receiving_field      --> Receiving field
 *   n_not_located_points --> Number of not located points
 *
 * returns :
 *   cwipi_exchange_status
 *
 *----------------------------------------------------------------------------*/

cwipi_exchange_status_t cwipi_exchange
(const char                          *coupling_id,
 const char                          *exchange_name,
 const int                            stride,
 const int                            time_step,
 const double                         time_value,
 const char                          *sending_field_name,
 const double                        *sending_field,
 const char                          *receiving_field_name,
 double                              *receiving_field,
 int                                 *n_not_located_points);

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
 int                       *request);

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
 double                    *receiving_field,
 int                       *request);

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
                       int          request);

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
                      int          request);

/*----------------------------------------------------------------------------
 *
 * Get located point distance to exchange area 
 *
 * parameters
 *   coupling_id          <-- Coupling identifier
 *
 *----------------------------------------------------------------------------*/

const float *cwipi_distance_located_pts_get(const char  *coupling_name);

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
(const char *coupling_id,
 cwipi_interpolation_fct_t fct);

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
(const char *coupling_id,
 cwipi_user_interp_ho_fct_t fct);

/*----------------------------------------------------------------------------
 *
 * Define a FORTRAN interpolation function
 *
 * parameters
 *   coupling_id          <-- Coupling identifier
 *   fct                  <-- Interpolation function
 *
 *----------------------------------------------------------------------------*/

void cwipi_set_interpolation_function_f
(const char *coupling_id,
 void* fct);

/*----------------------------------------------------------------------------
 *
 * Delete a coupling
 *
 * parameters
 *   coupling_id          <-- Coupling identifier
 *
 *----------------------------------------------------------------------------*/

void cwipi_delete_coupling(const char *coupling_id);

/*----------------------------------------------------------------------------
 *
 * Finalize cwipi. This is a synchronization point between all applications
 *
 *----------------------------------------------------------------------------*/

void cwipi_finalize(void);

/*----------------------------------------------------------------------------
 *
 * Get not located points
 *
 * parameters
 *   coupling_id          <-- Coupling identifier
 *
 * return
 *                        --> Not located points
 *----------------------------------------------------------------------------*/

const int * cwipi_get_not_located_points(const char *coupling_id);

/*----------------------------------------------------------------------------
 *
 * Get located points
 *
 * parameters
 *   coupling_id          <-- Coupling identifier
 *
 * return
 *   locatedPoints        <-- Located points
 *
 *----------------------------------------------------------------------------*/

const int * cwipi_get_located_points(const char *coupling_id);

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

int cwipi_get_n_located_points(const char *coupling_id);

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

int cwipi_get_n_not_located_points(const char *coupling_id);


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

int cwipi_get_n_distant_points(const char *coupling_id);


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

int cwipi_get_n_distant_ranks(const char *coupling_id);


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

const int *cwipi_get_distant_distribution(const char *coupling_id);


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

const int *cwipi_get_located_points_distribution(const char *coupling_id);

/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* __CWIPI_H__ */
