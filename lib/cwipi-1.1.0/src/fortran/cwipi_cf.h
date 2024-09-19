#ifndef __CWIPI_CF_H__
#define __CWIPI_CF_H__
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

#include <mpi.h>
#include <stdio.h>

/*=============================================================================
 * Macro definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Macro for handling of different symbol names (underscored or not,
 * lowercase or uppercase) between C and Fortran, for link resolution.
 *----------------------------------------------------------------------------*/

#if !defined (__hpux) && !defined (_AIX)
#define PROCF(x, y) x##_
#else
#define PROCF(x, y) x
#endif

/*----------------------------------------------------------------------------
 * Macro used to handle automatic "Fortran string length" arguments
 * (not used by CWIPI, but set by many compilers).
 * Some compilers, like the Fujitsu VPP 5000 compiler, may not
 * support the variable length lists in mixed C/Fortran calls.
 *----------------------------------------------------------------------------*/

#if defined (__uxpv__)  /* Fujitsu VPP 5000 case */
#define ARGF_SUPP_CHAINE
#else
#define ARGF_SUPP_CHAINE , ...
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
   ARGF_SUPP_CHAINE);

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
   ARGF_SUPP_CHAINE);

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
   ARGF_SUPP_CHAINE);


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
   ARGF_SUPP_CHAINE);

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
   int *initial_value
   ARGF_SUPP_CHAINE);

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
   double *initial_value
   ARGF_SUPP_CHAINE);

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
   ARGF_SUPP_CHAINE);

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
   ARGF_SUPP_CHAINE);

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
   ARGF_SUPP_CHAINE);


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
   char *value,
   int *l_value
   ARGF_SUPP_CHAINE);

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
   ARGF_SUPP_CHAINE);

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
   ARGF_SUPP_CHAINE);

/*----------------------------------------------------------------------------
 *
 * Get a integer control parameter of a other application
 *
 * parameters
 *    application_name       <-- application name
 *    name                   <-- parameter name
 *
 *----------------------------------------------------------------------------*/


void PROCF(cwipi_del_loc_str_ctrl_param_cf,
           CWIPI_DEL_LOC_STR_CTRL_PARAM_CF)
  (const char *name,
   const int  *l_name
   ARGF_SUPP_CHAINE);

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
   ARGF_SUPP_CHAINE);

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
   ARGF_SUPP_CHAINE);

/*----------------------------------------------------------------------------
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
   char *value,
   int *l_value
   ARGF_SUPP_CHAINE);

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
   ARGF_SUPP_CHAINE);

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
   ARGF_SUPP_CHAINE);


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
   ARGF_SUPP_CHAINE);


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
   ARGF_SUPP_CHAINE);


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
   ARGF_SUPP_CHAINE);

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
   ARGF_SUPP_CHAINE);


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
   ARGF_SUPP_CHAINE);

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
   ARGF_SUPP_CHAINE);

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
   ARGF_SUPP_CHAINE);

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

void PROCF(cwipi_synch_ctrl_param_cf,
           CWIPI_SYNCH_CTRL_PARAM_CF)
  (const char *application_name,
   const int  *l_application_name
   ARGF_SUPP_CHAINE);

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
 *  nbLocations             <-- maximun number of locations 
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
  ARGF_SUPP_CHAINE);

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
  (const char   *coupling_id,
   const int  *l_coupling_id,
   const int    *n_points,
   double *coordinate
   ARGF_SUPP_CHAINE);

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
   ARGF_SUPP_CHAINE);

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
   ARGF_SUPP_CHAINE);

void PROCF(cwipi_ho_options_set_cf,
           CWIPI_HO_OPTIONS_SET_CF)
  (const char *coupling_name,
   const int  *l_coupling_name,
   const char *option,
   const int  *l_option,
   const char *value,
   const int  *l_value
   ARGF_SUPP_CHAINE);

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
   ARGF_SUPP_CHAINE);

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
 const int *uvw_grid
 ARGF_SUPP_CHAINE);

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
 ARGF_SUPP_CHAINE);


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
 *   weight_tetra       <-- Weight computation in a tetrahedron
 *   weight_prism       <-- Weight computation in a prism
 *   weight_pyramid     <-- Weight computation in a pyramid
 *   weight_hexa        <-- Weight computation in a hexaedron
 *   weight_tria        <-- Weight computation on a triangle
 *   weight_quad        <-- Weight computation on a quandragle
 *   weight_edge        <-- Weight computation on a edge
 *   interp_tetra      <-- Interpolation in a tetrahedron
 *   interp_prism      <-- Interpolation in a prism
 *   interp_pyramid    <-- Interpolation in a pyramid
 *   interp_hexa       <-- Interpolation in a hexaedron
 *   interp_tria       <-- Interpolation on a triangle
 *   interp_quad       <-- Interpolation on a quandragle
 *   interp_edge       <-- Interpolation on a edge
 *
 *----------------------------------------------------------------------------*/

//void
//cwipi_ho_user_elementary_functions_set (cwipi_ho_location_fct_t location_tetra,
//                                        cwipi_ho_location_fct_t location_prism,
//                                        cwipi_ho_location_fct_t location_pyramid,
//                                        cwipi_ho_location_fct_t location_hexa,
//                                        cwipi_ho_location_fct_t location_tria,
//                                        cwipi_ho_location_fct_t location_quad,
//                                        cwipi_ho_location_fct_t location_edge,
//                                        cwipi_ho_weight_fct_t weight_tetra,
//                                        cwipi_ho_weight_fct_t weight_prism,
//                                        cwipi_ho_weight_fct_t weight_pyramid,
//                                        cwipi_ho_weight_fct_t weight_hexa,
//                                        cwipi_ho_weight_fct_t weight_tria,
//                                        cwipi_ho_weight_fct_t weight_quad,
//                                        cwipi_ho_weight_fct_t weight_edge,
//                                        cwipi_ho_interp_fct_t interp_tetra,
//                                        cwipi_ho_interp_fct_t interp_prism,
//                                        cwipi_ho_interp_fct_t interp_pyramid,
//                                        cwipi_ho_interp_fct_t interp_hexa,
//                                        cwipi_ho_interp_fct_t interp_tria,
//                                        cwipi_ho_interp_fct_t interp_quad,
//                                        cwipi_ho_interp_fct_t interp_edge);

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
                                                      ARGF_SUPP_CHAINE);


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
      ARGF_SUPP_CHAINE);

/*----------------------------------------------------------------------------
 *
 * Update Location
 *
 * parameters
 *   coupling_name        <-- Coupling identifier
 *----------------------------------------------------------------------------*/

void PROCF(cwipi_update_location_cf, CWIPI_UPDATE_LOCATION_CF) (const char *coupling_name,
                                                                const int  *l_coupling_name
                                                                ARGF_SUPP_CHAINE);


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
                                       ARGF_SUPP_CHAINE);

/*----------------------------------------------------------------------------
 *
 * Get distant points distance to location element
 *
 * parameters
 *   coupling_name        <-- Coupling identifier
 *   distance             --> Distant points distance to location element
 *----------------------------------------------------------------------------*/

void PROCF(cwipi_get_distant_distance_cf,
           CWIPI_GET_DISTANT_DISTANCE_CF) (const char *coupling_name,
                                           const int  *l_coupling_name,
                                           float *distance
                                           ARGF_SUPP_CHAINE);

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
                                               ARGF_SUPP_CHAINE);

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
                                                   ARGF_SUPP_CHAINE);

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
                                       ARGF_SUPP_CHAINE);

/*----------------------------------------------------------------------------
 *
 * Get number of located distant points
 *
 * parameters
 *   coupling_name           <-- Coupling identifier
 *   n_located_distant_Points --> Number of located distant points
 *
 *----------------------------------------------------------------------------*/

void PROCF(cwipi_get_n_located_dist_pts_cf, CWIPI_GET_N_LOCATED_DIST_PTS_CF) 
(const char *coupling_name,
 const int  *l_coupling_name,
 int *n_located_distant_Points
 ARGF_SUPP_CHAINE);

/*----------------------------------------------------------------------------
 *
 * Get distant points coordinates
 *
 * parameters
 *   coupling_name        <-- Coupling identifier
 *   coordinates          --> Distant points coordinates 
 *
 *----------------------------------------------------------------------------*/

void PROCF(cwipi_get_dis_coord_cf, CWIPI_GET_DIS_COORD_CF)
(const char *coupling_name,
 const int  *l_coupling_name,
 double *coordinates
 ARGF_SUPP_CHAINE
);

/*----------------------------------------------------------------------------
 *
 * Get located points barycentric coordinates index
 *
 * parameters
 *   coupling_name                <-- Coupling identifier
 *   barycentricCoordinatesIndex  --> located points barycentric
 *                                    coordinates index
 *----------------------------------------------------------------------------*/

void PROCF(cwipi_get_dis_bary_coord_idx_cf, CWIPI_GET_DIS_BARY_COORD_IDX_CF) 
(const char *coupling_name,
 const int  *l_coupling_name,
 int *barycentricCoordinatesIndex
 ARGF_SUPP_CHAINE);

/*----------------------------------------------------------------------------
 *
 * Get located points barycentric coordinates
 *
 * parameters
 *   coupling_name                <-- Coupling identifier
 *   barycentricCoordinatesIndex  --> located points barycentric
 *                                    coordinates
 *----------------------------------------------------------------------------*/

void PROCF(cwipi_get_distant_bary_coord_cf, CWIPI_GET_DISTANT_BARY_COORD_CF) 
(const char *coupling_name,
 const int  *l_coupling_name,
 double *barycentricCoordinates
 ARGF_SUPP_CHAINE);

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
 ARGF_SUPP_CHAINE);


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
 ARGF_SUPP_CHAINE);


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
 ARGF_SUPP_CHAINE);



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

void PROCF(cwipi_exchange_cf,
           CWIPI_EXCHANGE_CF)
  (const char      *coupling_name,
   const int       *l_coupling_name,
   const char      *exchange_name,
   const int       *l_exchange_name,
   const int       *exchange_dimension,
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
   ARGF_SUPP_CHAINE);

void PROCF(cwipi_exch_with_user_itp_cf,
           CWIPI_EXCH_WITH_USER_ITP_CF)
  (const char      *coupling_name,
   const int       *l_coupling_name,
   const char      *exchange_name,
   const int       *l_exchange_name,
   const int       *exchange_dimension,
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
   ARGF_SUPP_CHAINE);

void PROCF(cwipi_receive_cf,
           CWIPI_RECEIVE_CF)
  (const char      *coupling_name,
   const int       *l_coupling_name,
   const char      *exchange_name,
   const int       *l_exchange_name,
   const int       *exchange_dimension,
   const int       *n_step,
   const double    *time_value,
   char            *receiving_field_name,
   const int       *l_receiving_field_name,
   double          *receiving_field,
   int             *n_not_located_points,
   int             *exchange_status
   ARGF_SUPP_CHAINE);


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
   ARGF_SUPP_CHAINE);


void PROCF(cwipi_send_with_user_interpolation_cf,
           CWIPI_SEND_WITH_USER_INTERPOLATION_CF)
  (const char      *coupling_name,
   const int       *l_coupling_name,
   const char      *exchange_name,
   const int       *l_exchange_name,
   const int       *exchange_dimension,
   const int       *n_step,
   const double    *time_value,
   const char      *sending_field_name,
   const int       *l_sending_field_name,
   const double    *sending_field,
   void            *ptFortranInterpolationFct,
   int             *exchange_status
   ARGF_SUPP_CHAINE);


void PROCF(cwipi_send_cf,
           CWIPI_SEND_CF)
  (const char      *coupling_name,
   const int       *l_coupling_name,
   const char      *exchange_name,
   const int       *l_exchange_name,
   const int       *exchange_dimension,
   const int       *n_step,
   const double    *time_value,
   const char      *sending_field_name,
   const int       *l_sending_field_name,
   const double    *sending_field,
   int             *exchange_status
   ARGF_SUPP_CHAINE);

void PROCF(cwipi_wait_issend_cf, CWIPI_WAIT_ISSEND_CF)
  (const char      *coupling_name,
   const int       *l_coupling_name,
   const int       *request
   ARGF_SUPP_CHAINE);

void PROCF(cwipi_wait_irecv_cf, CWIPI_WAIT_IRECV_CF)
  (const char      *coupling_name,
   const int       *l_coupling_name,
   const int       *request
   ARGF_SUPP_CHAINE);

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
   ARGF_SUPP_CHAINE);

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
   ARGF_SUPP_CHAINE);

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
   ARGF_SUPP_CHAINE);

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
   const int  *l_coupling_name
   ARGF_SUPP_CHAINE);

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
   int *notLocatedPoints);

/*----------------------------------------------------------------------------
 *
 * Get located points
 *
 * parameters
 *   coupling_id          <-- Coupling identifier
 *   locatedPoints        --> Located points
 *
 *----------------------------------------------------------------------------*/

void PROCF(cwipi_get_located_pts_cf,
           CWIPI_GET_LOCATED_PTS_CF)
  (const char *coupling_name,
   const int  *l_coupling_name,
   int *locatedPoints);

/*----------------------------------------------------------------------------
 *
 * Finalize cwipi. This is a synchronization point between all applications
 *
 *----------------------------------------------------------------------------*/

void PROCF(cwipi_finalize_cf,
           CWIPI_FINALIZE_CF) ();

/*----------------------------------------------------------------------------
 *
 * Dump application properties
 *
 *----------------------------------------------------------------------------*/

void PROCF(cwipi_dump_appli_properties_cf,
           CWIPI_DUMP_APPLI_PROPERTIES_CF) ();


/*----------------------------------------------------------------------------
 *
 * Get distant elements that contain located points
 *
 * parameters
 *   coupling_id          <-- Coupling identifier
 *   elements             --> Element that contain located points
 *
 *----------------------------------------------------------------------------*/

void PROCF(cwipi_get_elt_cont_cf,
           CWIPI_GET_ELT_CONT_CF)
(const char *coupling_id,
 const int  *l_coupling_id,
 int *elements);

/*----------------------------------------------------------------------------
 *
 * Get number of vertices of distant elements that contain located points
 *
 * parameters
 *   coupling_id          <-- Coupling identifier
 *   n_vertices           --> Number of vertices of element that contain
 *                            located point
 *
 *----------------------------------------------------------------------------*/

void PROCF(cwipi_get_elt_cont_n_vtx_cf,
           CWIPI_GET_ELT_CONT_N_VTX_CF)
(const char *coupling_id,
 const int  *l_coupling_id,
 int *n_vertices);

/*----------------------------------------------------------------------------
 *
 * Get vertices id of distant elements that contain located points
 *
 * parameters
 *   coupling_id          <-- Coupling identifier
 *   vertices             --> Vertices id
 *
 *----------------------------------------------------------------------------*/

void PROCF(cwipi_get_elt_cont_vtx_cf,
           CWIPI_GET_ELT_CONT_VTX_CF)
(const char *coupling_id,
 const int  *l_coupling_id,
 int *vertices);

/*----------------------------------------------------------------------------
 *
 * Get vertices coordinates of distant elements that contain located points
 *
 * parameters
 *   coupling_id          <-- Coupling identifier
 *   coordinates          --> Vertices coordinates
 *
 *----------------------------------------------------------------------------*/

void PROCF(cwipi_get_elt_cont_vtx_coo_cf,
           CWIPI_GET_ELT_CONT_VTX_COO_CF)
(const char *coupling_id,
 const int  *l_coupling_id,
 double* coordinates);

/*----------------------------------------------------------------------------
 *
 * Get barycentric coords in distant elements for located points
 *
 * parameters
 *   coupling_id          <-- Coupling identifier
 *   coordinates          --> Barycentric coordinates
 *
 *----------------------------------------------------------------------------*/

void PROCF(cwipi_get_elt_cont_bar_coo_cf,
           CWIPI_GET_ELT_CONT_BAR_COO_CF)
(const char *coupling_id,
 const int  *l_coupling_id,
 double *coordinates);

/*----------------------------------------------------------------------------
 *
 * For each located point get the MPI rank of distant element
 *
 * parameters
 *   coupling_id          <-- Coupling identifier
 *   MPIranks             --> MPI ranks that contains located point
 *
 *----------------------------------------------------------------------------*/

void PROCF(cwipi_get_elt_cont_mpi_rank_cf,
           CWIPI_GET_ELT_CONT_MPI_RANK_CF)
(const char *coupling_id,
 const int  *l_coupling_id,
 int *MPIrank);

/*----------------------------------------------------------------------------
 *
 * Exchange Fields on vertices of element containing each located point
 *
 * parameters
 *   coupling_id          <-- Coupling identifier
 *   sendingField         <-- Field defined on local mesh vertices
 *   receivingField       --> Field defined on vertices of distant
 *                            elements that contain each located point
 *   stride               <-- Number of field component
 *
 *----------------------------------------------------------------------------*/

void PROCF(cwipi_exch_cellvtxfd_eltcont_cf,
           CWIPI_EXCH_CELLVTXFD_ELTCONT_CF)
(const char *coupling_id,
 const int  *l_coupling_id,
 double *sendingField,
 double *receivingField,
 const int *stride);

void PROCF(cwipi_send_cellvtxfd_eltcont_cf,
           CWIPI_SEND_CELLVTXFD_ELTCONT_CF)
(const char *coupling_id,
 const int  *l_coupling_id,
 double *sendingField,
 const int *stride);

void PROCF(cwipi_recv_cellvtxfd_eltcont_cf,
           CWIPI_RECV_CELLVTXFD_ELTCONT_CF)
(const char *coupling_id,
 const int  *l_coupling_id,
 double *receivingField,
 const int *stride);
/*----------------------------------------------------------------------------
 *
 * Exchange field on cells that contain each located points
 *
 * parameters
 *   coupling_id          <-- Coupling identifier
 *   sendingField         <-- Field defined on local mesh vertices
 *   receivingField       --> Field defined on vertices of distant
 *                            elements that contain each located point
 *   stride               <-- Number of field component
 *
 *----------------------------------------------------------------------------*/

void PROCF(cwipi_exch_cellcenfd_eltcont_cf,
           CWIPI_EXCH_CELLCENFD_ELTCONT_CF)
(const char *coupling_id,
 const int  *l_coupling_id,
 double *sendingField,
 double *receivingField,
 const int *stride);

void PROCF(cwipi_send_cellcenfd_eltcont_cf,
           CWIPI_SEND_CELLCENFD_ELTCONT_CF)
(const char *coupling_id,
 const int  *l_coupling_id,
 double *sendingField,
 const int *stride);

void PROCF(cwipi_recv_cellcenfd_eltcont_cf,
           CWIPI_RECV_CELLCENFD_ELTCONT_CF)
(const char *coupling_id,
 const int  *l_coupling_id,
 double *receivingField,
 const int *stride);

/*----------------------------------------------------------------------------
 *
 * Set coupling info
 *
 * parameters
 *   coupling_id          <-- Coupling identifier
 *   info                 <-- Coupling info
 *----------------------------------------------------------------------------*/

void PROCF(cwipi_set_info_cf,
           CWIPI_SET_INFO_CF)
(const char *coupling_id,
 const int  *l_coupling_id,
 const int *info);

/*----------------------------------------------------------------------------
 *
 * cwipi_set_location_index
 *
 * parameters
 *   coupling_id          <-- Coupling identifier
 *   index                <-- location index
 *----------------------------------------------------------------------------*/

void PROCF(cwipi_set_location_index_cf,
	   CWIPI_SET_LOCATION_INDEX_CF)
(const char *coupling_id,
 const int  *l_coupling_id,
 const int  *index
 ARGF_SUPP_CHAINE);


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
 ARGF_SUPP_CHAINE);

void PROCF(cwipi_save_location_cf,
           CWIPI_SAVE_LOCATION_CF)
(const char *coupling_name, 
 const int  *l_coupling_name
 ARGF_SUPP_CHAINE);


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
  const int  *l_filename_name,
  const char *mode,
  const int  *l_mode
 ARGF_SUPP_CHAINE);

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
  ARGF_SUPP_CHAINE);


/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* __CWIPI_H__ */
