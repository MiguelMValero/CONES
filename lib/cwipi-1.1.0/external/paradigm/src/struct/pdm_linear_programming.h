/*
 * \file
 */

#ifndef __PDM_LINEAR_PROGRAMMING_H__
#define __PDM_LINEAR_PROGRAMMING_H__

/*============================================================================
 * Search octrees and quadtrees of boxes.
 *============================================================================*/

/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <stdio.h>

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "pdm_box.h"

/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force back Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */


/*============================================================================
 * Type
 *============================================================================*/

typedef enum {

  PDM_LP_FEASIBLE,
  PDM_LP_UNFEASIBLE,
  PDM_LP_UNBOUNDED

} PDM_lp_status_t;

/*=============================================================================
 * Public function definitions
 *============================================================================*/

/**
 *
 * \brief Solve the d-dimensional linear optimization problem
 *        maximize c.x
 *        subject to constraints ai.x <= bi
 *
 * \param [in]     dim   Dimension
 * \param [in]     n     Number of inequality constraints
 * \param [in]     a     a in ax <= b
 * \param [in]     b     b in ax <= b
 * \param [in]     c     Constant in the objective function
 * \param [inout]  x     Initial point - Optimum
 *
 */

PDM_lp_status_t
PDM_lp_solve_nd
(
 const int  dim,
 const int  n,
 double    *a,
 double    *b,
 double    *c,
 double    *x
 );

/**
 *
 * \brief Determine if the current box intersects a given volume
 *
 * \param [in]   n_plane          Number of planes in the current volume
 * \param [in]   plane_origin     Coordinates of a point on each plane
 * \param [in]   plane_normal     Normal vector of each plane
 * \param [in]   box_extents      Extents of the box (x_min, y_min, z_min, x_max, y_max, z_max)
 *
 */

int
PDM_lp_intersect_volume_box
(
 const int  n_plane,
 double    *plane_origin,
 double    *plane_normal,
 double    *box_extents
);

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* __PDM_LINEAR_PROGRAMMING_H__ */
