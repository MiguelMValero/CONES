/*
 * \file
 */

#ifndef __PDM_LAGRANGE_TO_BEZIER_H__
#define __PDM_LAGRANGE_TO_BEZIER_H__

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "pdm.h"
#include "pdm_mesh_nodal.h"
#include "pdm_mpi.h"
#include "pdm_io.h"

/*=============================================================================
 * Macro definition
 *============================================================================*/


/*============================================================================
 * Types definition
 *============================================================================*/

/*=============================================================================
 * Global variables
 *============================================================================*/

/*=============================================================================
 * Static global variables
 *============================================================================*/

/*=============================================================================
 * Public function interfaces
 *============================================================================*/


/**
 *
 * \brief Convert Lagrange nodes to Bézier nodes for a high-order bar
 *
 * \param [in]  order   Element order
 * \param [in]  lag     Coordinates of Lagrange nodes
 * \param [out] bez     Coordinates of Bézier nodes
 * \param [in]  matrix  Optional work array to store the transition matrix
 *
 */

void
PDM_lagrange_to_bezier_bar
(
 const int     order,
       double *lag,
       double *bez,
       double *matrix
);


/**
 *
 * \brief Convert Lagrange nodes to Bézier nodes for a high-order triangle
 *
 * \param [in]  order   Element order
 * \param [in]  lag     Coordinates of Lagrange nodes
 * \param [out] bez     Coordinates of Bézier nodes
 * \param [in]  matrix  Optional work array to store the transition matrix
 *
 */

void
PDM_lagrange_to_bezier_tria
(
 const int     order,
       double *lag,
       double *bez,
       double *matrix
);


/**
 *
 * \brief Convert Lagrange nodes to Bézier nodes for a high-order quadrangle
 *
 * \param [in]  order   Element order
 * \param [in]  lag     Coordinates of Lagrange nodes
 * \param [out] bez     Coordinates of Bézier nodes
 * \param [in]  matrix  Optional work array to store the transition matrix
 *
 */

void
PDM_lagrange_to_bezier_quad
(
 const int     order,
       double *lag,
       double *bez,
       double *matrix
);


/**
 *
 * \brief Convert Lagrange nodes to Bézier nodes for a high-order tetrahedron
 * /!\ Only the boundary nodes are actually converted (for bounding-box computation)
 *
 * \param [in]  order   Element order
 * \param [in]  lag     Coordinates of Lagrange nodes
 * \param [out] bez     Coordinates of Bézier nodes
 * \param [in]  matrix  Optional work array to store the transition matrix
 *
 */

void
PDM_lagrange_to_bezier_tetra
(
 const int     order,
       double *lag,
       double *bez,
       double *matrix
);


/**
 *
 * \brief Convert Lagrange nodes to Bézier nodes for a high-order pyramid
 * /!\ Only the boundary nodes are actually converted (for bounding-box computation)
 *
 * \param [in]  order   Element order
 * \param [in]  lag     Coordinates of Lagrange nodes
 * \param [out] bez     Coordinates of Bézier nodes
 * \param [in]  matrix  Optional work array to store the transition matrix
 *
 */

void
PDM_lagrange_to_bezier_pyramid
(
 const int     order,
       double *lag,
       double *bez,
       double *matrix
);


/**
 *
 * \brief Convert Lagrange nodes to Bézier nodes for a high-order prism
 * /!\ Only the boundary nodes are actually converted (for bounding-box computation)
 *
 * \param [in]  order   Element order
 * \param [in]  lag     Coordinates of Lagrange nodes
 * \param [out] bez     Coordinates of Bézier nodes
 * \param [in]  matrix  Optional work array to store the transition matrix
 *
 */

void
PDM_lagrange_to_bezier_prism
(
 const int     order,
       double *lag,
       double *bez,
       double *matrix
);


/**
 *
 * \brief Convert Lagrange nodes to Bézier nodes for a high-order hexahedron
 * /!\ Only the boundary nodes are actually converted (for bounding-box computation)
 *
 * \param [in]  order   Element order
 * \param [in]  lag     Coordinates of Lagrange nodes
 * \param [out] bez     Coordinates of Bézier nodes
 * \param [in]  matrix  Optional work array to store the transition matrix
 *
 */

void
PDM_lagrange_to_bezier_hexa
(
 const int     order,
       double *lag,
       double *bez,
       double *matrix
);

/* Get bezier bounding box (copied from pdm_t_dcube_nodal_gen.c) */
void
PDM_bezier_bounding_boxes
(
 const PDM_Mesh_nodal_elt_t  t_elt,
 const int                   order,
 const int                   n_nodes,
 int                         n_elt,
 double                     *lagrange_coord,
 double                    **extents
);

#endif /* __PDM_LAGRANGE_TO_BEZIER_H__ */
