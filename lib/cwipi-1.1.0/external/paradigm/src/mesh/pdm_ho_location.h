/*
 * \file
 */

#ifndef __PDM_HO_LOCATION_H__
#define __PDM_HO_LOCATION_H__

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "pdm.h"
#include "pdm_mesh_nodal.h"

/*=============================================================================
 * Macro definitions
 *============================================================================*/

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force back Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */

/*============================================================================
 * Type
 *============================================================================*/

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

typedef double (*PDM_ho_location_fct_t)
(const int     entities_dim,
 const int     order,
 const int     n_nodes,
 const double *nodes_coords,
 const double *point_coords,
 double       *projected_coords,
 double       *uvw);

/*=============================================================================
 * Static global variables
 *============================================================================*/

/*=============================================================================
 * Public function prototypes
 *============================================================================*/


/**
 * \brief Point location in a high order element
 *
 * \param [in]   type              Element type
 * \param [in]   order             Element order
 * \param [in]   n_nodes           Number of nodes
 * \param [in]   nodes_coords      Coordinates of the nodes (size = 3 * \ref n_nodes)
 * \param [in]   point_coords      Coordinates of the point to locate (size = 3)
 * \param [out]  projected_coords  Coordinates of the projection on the element (size = 3)
 * \param [out]  uvw               Parametric coordinates of the projection on the element
 *
 * \return       Squared distance from the point to the element
 *
 */

double
PDM_ho_location
(
 const PDM_Mesh_nodal_elt_t  type,
 const int                   order,
 const int                   n_nodes,
 const double               *nodes_coords,
 const double               *point_coords,
 double                     *projected_coords,
 double                     *uvw
 );


/**
 * \brief Get the parametric coordinates of the Lagrange nodes
 *
 * \param [in]   type              Element type
 * \param [in]   order             Element order
 * \param [in]   umin              Lower bound on u-coordinate
 * \param [in]   umax              Upper bound on u-coordinate
 * \param [in]   vmin              Lower bound on v-coordinate
 * \param [in]   vmax              Upper bound on v-coordinate
 * \param [in]   wmin              Lower bound on w-coordinate
 * \param [in]   wmax              Upper bound on w-coordinate
 * \param [out]  uvw_node          Parametric coordinates of the Lagrange nodes
 *
 */

void
PDM_ho_location_uvw_nodes
(
 const PDM_Mesh_nodal_elt_t  type,
 const int                   order,
 const double                umin,
 const double                umax,
 const double                vmin,
 const double                vmax,
 const double                wmin,
 const double                wmax,
       double               *uvw_node
 );


/**
 * \brief Point location in a high order element
 *
 * \param [in]   type              Element type
 * \param [in]   order             Element order
 * \param [in]   n_nodes           Number of nodes
 * \param [in]   nodes_coords      Coordinates of the nodes (size = 3 * \ref n_nodes)
 * \param [in]   point_coords      Coordinates of the point to locate (size = 3)
 * \param [in]   tolerance         Tolerance for Newton step in uvw-space
 * \param [out]  projected_coords  Coordinates of the projection on the element (size = 3)
 * \param [out]  uvw               Parametric coordinates of the projection on the element
 * \param [out]  converged         Convergence status (1 if successful, 0 else)
 * \param [in]   work_array        Optional work array (size = (elt_dim+1) * \ref n_nodes or NULL)
 *
 * \return       Squared distance from the point to the element
 *
 */

double
PDM_ho_location_newton
(
 const PDM_Mesh_nodal_elt_t  type,
 const int                   order,
 const int                   n_nodes,
 const double               *nodes_coords,
 const double               *point_coords,
 const double                tolerance,
       double               *projected_coords,
       double               *uvw,
       int                  *converged,
       double               *work_array
);

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* __PDM_HO_LOCATION_H__ */
