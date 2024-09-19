/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "pdm.h"
#include "pdm_priv.h"
#include "pdm_logging.h"
#include "pdm_plane.h"
#include "pdm_line.h"
#include "pdm_polygon.h"
#include "pdm_geom_elem.h"
#include "pdm_predicate.h"
#include "pdm_mean_values.h"

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

/*=============================================================================
 * Private function definition
 *============================================================================*/

/**
 * \brief Return a random value in [min, max]
 *
 * \param [in]  min  Minimum
 * \param [in]  max  Maximum
 *
 * \return      Random value
 */

static double
_randomVal
(
 const double min,
 const double max
)
{
  double resultat = ((double)rand())/((double)RAND_MAX);
  resultat = min + resultat * (max - min);

  return resultat;
}


/*=============================================================================
 * Public function definition
 *============================================================================*/

/**
 * \brief Get bounds
 *
 * \param [in]  numPts   Number of polygon vertices
 * \param [in]  pts      Polygon vertices coordinates
 *
 * \return      Bounds
 *
 */


double *
PDM_polygon_bounds_get
(
 const int     numPts,
 const double *pts
)
{
  double *bounds = malloc (sizeof(double) * 6);

  bounds[0] = DBL_MAX;
  bounds[1] = -DBL_MAX;
  bounds[2] = DBL_MAX;
  bounds[3] = -DBL_MAX;
  bounds[4] = DBL_MAX;
  bounds[5] = -DBL_MAX;

  for (int isom = 0; isom < numPts; isom++) {
    for (int l = 0; l < 3; l++) {
      double coord = pts[3*isom + l];
      if (bounds[2*l] > coord) {
        bounds[2*l] = coord;
      }
      if (bounds[2*l+1] < coord) {
        bounds[2*l+1] = coord;
      }
    }
  }

  return bounds;
}


/**
 * \brief Evaluates the position in a polygon
 *
 * \param [in]  x        Point coordinates to evaluate position
 * \param [in]  numPts   Number of polygon vertices
 * \param [in]  pts      Polygon vertices coordinates
 * \param [out] closest  Closest Point in Polygon or NULL
 * \param [out] minDist2 Square of the distance
 *
 * \return      \ref PDM_POLYGON_INSIDE or \ref PDM_POLYGON_OUTSIDE
 *              if the projected is in the polygon or not
 */

/*  This function is derived from VTK                                      */
/*  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen                */
/*  All rights reserved.                                                   */
/*  See Copyright.txt or http://www.kitware.com/Copyright.htm for details. */

PDM_polygon_status_t
PDM_polygon_evaluate_position
(
 const double  x[3],
 const int     numPts,
 const double *pts,
 double        closestPoint[3],
 double        *minDist2
 )
{
  double p0[3];
  double p10[3];
  double l10;
  double p20[3];
  double l20;
  double n[3];
  double cp[3];
  double ray[3];
  double bary[3];

  double pts_p[3*numPts];

  /*
   * Average plane
   */

  PDM_plane_normal (numPts, pts, n);

  PDM_polygon_compute_barycenter (numPts, pts, bary);

  for (int k = 0; k < numPts; k++) {
    double *pt = (double *) pts + 3*k;
    double *pt_p = pts_p + 3*k;
    PDM_plane_projection2 (pt, bary, n, pt_p);
  }

  double *_pts_p = pts_p;

  PDM_polygon_parameterize (numPts, _pts_p, p0, p10, &l10, p20, &l20, n);

  PDM_plane_projection (x,p0,n,cp);

  for (int i = 0; i < 3; i++) {
    ray[i] = cp[i] - p0[i];
  }

  double pcoords[3];

  pcoords[0] = PDM_DOT_PRODUCT(ray,p10) / (l10*l10);
  pcoords[1] = PDM_DOT_PRODUCT(ray,p20) / (l20*l20);
  pcoords[2] = 0.0;

  // double bounds[6] = {DBL_MAX, -DBL_MAX,
  //                     DBL_MAX, -DBL_MAX,
  //                     DBL_MAX, -DBL_MAX};

  // for (int isom = 0; isom < numPts; isom++) {
  //   for (int l = 0; l < 3; l++) {
  //     double coord = _pts_p[3*isom + l];
  //     if (bounds[2*l] > coord) {
  //       bounds[2*l] = coord;
  //     }
  //     if (bounds[2*l+1] < coord) {
  //       bounds[2*l+1] = coord;
  //     }
  //   }
  // }
  double *bounds = NULL;

  if (pcoords[0] >= 0.0 && pcoords[0] <= 1.0 &&
      pcoords[1] >= 0.0 && pcoords[1] <= 1.0 &&
      (PDM_polygon_point_in_new (cp, numPts, _pts_p,
                                 bounds, n) == PDM_POLYGON_INSIDE) ) {
    if (closestPoint) {
      closestPoint[0] = cp[0];
      closestPoint[1] = cp[1];
      closestPoint[2] = cp[2];
      double v[3] = {x[0] - closestPoint[0],
                     x[1] - closestPoint[1],
                     x[2] - closestPoint[2]};

      *minDist2 = PDM_DOT_PRODUCT (v, v);
    }
    return PDM_POLYGON_INSIDE;
  }

  /*
   * If here, point is outside of polygon, so need to find distance to boundary
   */

  else {
    double t, dist2;
    double closest[3];
    double *pt1, *pt2;

    if (closestPoint) {
      *minDist2 = DBL_MAX;
      for (int i=0; i<numPts; i++) {
        pt1 = (double *) pts_p + 3 * i;
        pt2 = (double *) pts_p + 3 * ((i+1)%numPts);
        dist2 = PDM_line_distance (x, pt1, pt2, &t, closest);
        if ( dist2 < *minDist2 ) {
          closestPoint[0] = closest[0];
          closestPoint[1] = closest[1];
          closestPoint[2] = closest[2];
          *minDist2 = dist2;
        }
      }
    }
    return PDM_POLYGON_OUTSIDE;
  }
}


/**
 * \brief Computes polygon parametrization
 *
 * \param [in]  numPts  Number of polygon vertices
 * \param [in]  pts     Polygon vertices coordinates
 * \param [out] p0,     Origin vertex
 * \param [out] p10,    First edge vector
 * \param [out] l10,    First edge length
 * \param [out] p20,    First edge vector
 * \param [out] l20,    Second edge vector
 * \param [out] n       Normal
 *
 * \return      \ref PDM_TRUE except for a triangle
 *
 */

/*  This function is derived from VTK                                      */
/*  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen                */
/*  All rights reserved.                                                   */
/*  See Copyright.txt or http://www.kitware.com/Copyright.htm for details. */

PDM_bool_t
PDM_polygon_parameterize
(
 const int     numPts,
 const double *pts,
 double       *p0,
 double       *p10,
 double       *l10,
 double       *p20,
 double       *l20,
 double       *n
)
{
  double s, t, p[3], p1[3], p2[3], sbounds[2], tbounds[2];

  if (numPts < 3) {
    return PDM_FALSE;
  }

  /*
   *  This is a two pass process: first create a p' coordinate system
   *  that is then adjusted to insure that the polygon points are all in
   *  the range 0<=s,t<=1.  The p' system is defined by the polygon normal,
   *  first vertex and the first edge.
   */

  PDM_plane_normal (numPts, pts, n);

  double x1[3];
  x1[0] = pts[0];
  x1[1] = pts[1];
  x1[2] = pts[2];

  double x2[3];
  x2[0] = pts[3];
  x2[1] = pts[3+1];
  x2[2] = pts[3+2];

  for (int i = 0; i < 3; i++) {
    p0[i] = x1[i];
    p10[i] = x2[i] - x1[i];
  }

  PDM_CROSS_PRODUCT(p20,n,p10);

  /*
   * Determine lengths of edges
   */

PDM_GCC_SUPPRESS_WARNING_WITH_PUSH("-Wfloat-equal")

  if ( ((*l10)= PDM_DOT_PRODUCT(p10,p10)) == 0.0
    || ((*l20)= PDM_DOT_PRODUCT(p20,p20)) == 0.0 ) {
    return PDM_FALSE;
  }

PDM_GCC_SUPPRESS_WARNING_POP

  /*
   *  Now evalute all polygon points to determine min/max parametric
   *  coordinate values.
   *
   * first vertex has (s,t) = (0,0)
   */

  sbounds[0] = 0.0;
  sbounds[1] = 0.0;

  tbounds[0] = 0.0;
  tbounds[1] = 0.0;

  for (int i = 1; i < numPts; i++) {
    x1[0] = pts[3*i];
    x1[1] = pts[3*i+1];
    x1[2] = pts[3*i+2];
    for (int j = 0; j < 3; j++) {
      p[j] = x1[j] - p0[j];
    }

    s = PDM_DOT_PRODUCT (p,p10) / (*l10);
    t = PDM_DOT_PRODUCT (p,p20) / (*l20);

    sbounds[0] = (s<sbounds[0]?s:sbounds[0]);
    sbounds[1] = (s>sbounds[1]?s:sbounds[1]);
    tbounds[0] = (t<tbounds[0]?t:tbounds[0]);
    tbounds[1] = (t>tbounds[1]?t:tbounds[1]);
  }

  /*
   * Re-evaluate coordinate system
   */

  for (int i = 0; i < 3; i++) {
    p1[i] = p0[i] + sbounds[1]*p10[i] + tbounds[0]*p20[i];
    p2[i] = p0[i] + sbounds[0]*p10[i] + tbounds[1]*p20[i];
    p0[i] = p0[i] + sbounds[0]*p10[i] + tbounds[0]*p20[i];
    p10[i] = p1[i] - p0[i];
    p20[i] = p2[i] - p0[i];
  }

  (*l10) = PDM_MODULE(p10);
  (*l20) = PDM_MODULE(p20);

  return PDM_TRUE;
}


/**
 * \brief Computes polygon parametrization
 *
 * \param [in]  x        Point coordinates to evaluate position
 * \param [in]  numPts  Number of polygon vertices
 * \param [in]  pts     Polygon vertices coordinates
 * \param [out] p0,     Origin vertex
 * \param [out] p10,    First edge vector
 * \param [out] l10,    First edge length
 * \param [out] p20,    First edge vector
 * \param [out] l20,    Second edge vector
 * \param [out] n       Normal
 *
 * \return      \ref Status inside, outside or degenerated
 *
 */

/*  This function is derived from VTK                                      */
/*  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen                */
/*  All rights reserved.                                                   */
/*  See Copyright.txt or http://www.kitware.com/Copyright.htm for details. */


#define _TOL_POLY 1.e-05
#define _POLYGON_INTERSECTION 2
#define _POLYGON_ON_LINE 3

#define _POLYGON_CERTAIN 1
#define _POLYGON_UNCERTAIN 0
#define _POLYGON_RAY_TOL 1.e-03 //Tolerance for ray firing
#define _POLYGON_MAX_ITER 10    //Maximum iterations for ray-firing
#define _POLYGON_VOTE_THRESHOLD 2


PDM_polygon_status_t
PDM_polygon_point_in
(
 const double  x[3],
 const int     numPts,
 const double *pts,
 double       *bounds,
 double       *n
)
{
  double *x1, *x2, xray[3], u, v;
  double rayMag, mag=1, ray[3];
  int rayOK, status;
  int iterNumber;
  int maxComp, comps[2];
  int deltaVotes;

  /*
   * Do a quick bounds check
   */

  if ( x[0] < bounds[0] || x[0] > bounds[1] ||
       x[1] < bounds[2] || x[1] > bounds[3] ||
       x[2] < bounds[4] || x[2] > bounds[5]) {

    return PDM_POLYGON_OUTSIDE;
  }

  /*
   * Define a ray to fire. The ray is a random ray normal to the
   *  normal of the face. The length of the ray is a function of the
   *  size of the face bounding box.
   */

  for (int i = 0; i < 3; i++) {
    ray[i] = ( bounds[2*i+1] - bounds[2*i] )*1.1 +
          fabs((bounds[2*i+1] + bounds[2*i])/2.0 - x[i]);
  }

  if ((rayMag = PDM_MODULE(ray)) < 1.e-15 ) {
    return PDM_POLYGON_OUTSIDE;
  }

  /* Get the maximum component of the normal. */

  if (fabs(n[0]) > fabs(n[1])) {
    if (fabs(n[0]) > fabs(n[2])) {
      maxComp  = 0;
      comps[0] = 1;
      comps[1] = 2;
    }
    else {
      maxComp  = 2;
      comps[0] = 0;
      comps[1] = 1;
    }
  }
  else {
    if (fabs(n[1]) > fabs(n[2])) {
      maxComp  = 1;
      comps[0] = 0;
      comps[1] = 2;
    }
    else {
      maxComp = 2;
      comps[0] = 0;
      comps[1] = 1;
    }
  }

  /* Check that max component is non-zero */

  if (fabs(n[maxComp]) < 1.e-15) {
    return PDM_POLYGON_DEGENERATED;
  }

  /*
   * Enough information has been acquired to determine the random ray.
   * Random rays are generated until one is satisfactory (i.e.,
   * produces a ray of non-zero magnitude).  Also, since more than one
   * ray may need to be fired, the ray-firing occurs in a large loop.
   *
   * The variable iterNumber counts the number of iterations and is
   * limited by the defined variable _POLYGON_MAX_ITER.
   *
   * The variable deltaVotes keeps track of the number of votes for
   * "in" versus "out" of the face.  When delta_vote > 0, more votes
   * have counted for "in" than "out".  When delta_vote < 0, more votes
   * have counted for "out" than "in".  When the delta_vote exceeds or
   * equals the defined variable _POLYGON_VOTE_THRESHOLD, than the
   * appropriate "in" or "out" status is returned.
   */

  for (deltaVotes = 0, iterNumber = 1;
       (iterNumber < _POLYGON_MAX_ITER)
         && (PDM_ABS(deltaVotes) < _POLYGON_VOTE_THRESHOLD);
       iterNumber++) {

    /*
     * Generate ray
     */

    for (rayOK = PDM_FALSE; rayOK == PDM_FALSE; ) {
      ray[comps[0]] = _randomVal (-rayMag, rayMag);
      ray[comps[1]] = _randomVal (-rayMag, rayMag);
      ray[maxComp] = -(n[comps[0]]*ray[comps[0]] +
                       n[comps[1]]*ray[comps[1]]) / n[maxComp];
      if ( (mag = PDM_MODULE(ray)) > rayMag * _TOL_POLY ) {
        rayOK = PDM_TRUE;
      }
    }

    /*
     * The ray must be appropriately sized.
     */

    for (int i = 0; i < 3; i++) {
      xray[i] = x[i] + (rayMag/mag)*ray[i];
    }

    /*
     * The ray may now be fired against all the edges
     */


    int testResult = _POLYGON_CERTAIN;
    int numInts = 0;
    for (int i = 0; i < numPts; i++) {
      x1 = (double *) pts + 3*i;
      x2 = (double *) pts + 3*((i+1)%numPts);

      /*
       * Fire the ray and compute the number of intersections.  Be careful
       * of degenerate cases (e.g., ray intersects at vertex).
       */

      if ((status = PDM_line_intersection_mean_square(x,xray,x1,x2, &u,&v)) ==
          PDM_LINE_INTERSECT_YES) {

        /*
         * This test checks for vertex and edge intersections
         * For example
         *  Vertex intersection
         *    (u=0 v=0), (u=0 v=1), (u=1 v=0), (u=1 v=0)
         *  Edge intersection
         *    (u=0 v!=0 v!=1), (u=1 v!=0 v!=1)
         *    (u!=0 u!=1 v=0), (u!=0 u!=1 v=1)
         */

        if ( (0. < u) && (u < 1.0) &&
             (0. < v) && (v < 1.0) ) {
        /* if ( (_POLYGON_RAY_TOL < u) && (u < 1.0 - _POLYGON_RAY_TOL) && */
        /*      (_POLYGON_RAY_TOL < v) && (v < 1.0 - _POLYGON_RAY_TOL) ) { */
          numInts++;
        }
        else {
          testResult = _POLYGON_UNCERTAIN;
        }
      }

      else if ( status == _POLYGON_ON_LINE ) {
        testResult = _POLYGON_UNCERTAIN;
      }

    }

    if ( testResult == _POLYGON_CERTAIN ) {
      if ( numInts % 2 == 0) {
        --deltaVotes;
      }
      else {
        ++deltaVotes;
      }
    }
  } /* try another ray */

  /*
   * If the number of intersections is odd, the point is in the polygon.
   */

  if (deltaVotes <= 0) {
    return PDM_POLYGON_OUTSIDE;
  }
  else {
    return PDM_POLYGON_INSIDE;
  }
}

#undef _TOL_POLY
#undef _POLYGON_INTERSECTION
#undef _POLYGON_ON_LINE

#undef _POLYGON_CERTAIN
#undef _POLYGON_UNCERTAIN
#undef _POLYGON_RAY_TOL
#undef _POLYGON_MAX_ITER
#undef _POLYGON_VOTE_THRESHOLD


/**
 * \brief Computes polygon barycenter
 *
 * \param [in]   numPts  Number of polygon vertices
 * \param [in]   pts     Polygon vertices coordinates
 * \param [out]  bary    Barycenter
 *
 *
 */


void
PDM_polygon_compute_barycenter
(
const int numPts,
const double *pts,
double bary[3]
 )
{
  bary[0] = 0.;
  bary[1] = 0.;
  bary[2] = 0.;

  for (int i = 0; i < 3; i++) {
    for (int ipt = 0; ipt < numPts; ipt++) {
      bary[i] += pts[3*ipt+i];
    }
    bary[i] /= numPts;
  }
}

/**
 * Warning : unstable function (Rather use PDM_polygon_point_in )
 *
 * \brief Test if a point is inside a 2d polygon using the Winding Number method
 *        (see http://geomalgorithms.com/a03-_inclusion.html)
 *
 * \param [in]  xy            Point (x,y)-coordinates
 * \param [in]  n_vtx         Number of polygon vertices
 * \param [in]  vtx_xy        Polygon vertices (x,y)-coordinates
 * \param [in]  char_length   Characteristic length (used to scale tolerance)
 * \param [in]  bounds        Bounds (xmin, xmax, ymin, ymax)
 *
 * \return      \ref Status inside, outside or degenerated
 *
 */

PDM_polygon_status_t PDM_polygon_point_in_2d
(
 const double  xy[2],
 const int     n_vtx,
 const double *vtx_xy,
 // const double  char_length,
 double        bounds[4]
 )
{
  //const double eps_base = 1e-16;

  //  const double eps = eps_base * char_length;
  //const double eps2 = eps * eps;
  const double eps = 1e-16;
  const double eps2 = 1e-16;

  /*
   * Do a quick bounds check
   */
  if (bounds != NULL) {
    if ( xy[0] < bounds[0] || xy[0] > bounds[1] ||
         xy[1] < bounds[2] || xy[1] > bounds[3]) {
      return PDM_POLYGON_OUTSIDE;
    }
  }


  /*
   * Compute winding number
   */
  int wn = 0;

  /* Loop over edges */
  for (int i = 0; i < n_vtx; i++) {
    int ip = (i+1) % n_vtx;
    const double *p0 = vtx_xy + 2*i;
    const double *p1 = vtx_xy + 2*ip;

    double dx = xy[0] - p0[0];
    double dy = xy[1] - p0[1];

    /* Point coincident with vertex */
    if (dx*dx + dy*dy < eps2) {
      return PDM_POLYGON_INSIDE;
    }


    if (dy >= 0) {

      if (p1[1] > xy[1]) {
        /* Upward crossing */
        double ex = p1[0] - p0[0];
        double ey = p1[1] - p0[1];

        double s = ex*dy - ey*dx;
        if (s > eps) {
          wn++;
        } else if (s > -eps) {
          double denom = ex*ex + ey*ey;
          if (denom < eps2) {
            /* Vertices i and ip are coincident */
            return PDM_POLYGON_DEGENERATED;
          } else {
            /* Compute parameter along edge (i, ip) */
            double t = (ex*dx + ey*dy) / denom;
            if (t > -eps && t < 1 + eps) {
              return PDM_POLYGON_INSIDE;
            }
          }
        }
      }

    } else {// dy < 0

      if (p1[1] < xy[1]) {
        /* Downward crossing */
        double ex = p1[0] - p0[0];
        double ey = p1[1] - p0[1];

        double s = ex*dy - ey*dx;
        if (s < -eps) {
          wn--;
        } else if (s < eps) {
          double denom = ex*ex + ey*ey;
          if (denom < eps2) {
            /* Vertices i and ip are coincident */
            return PDM_POLYGON_DEGENERATED;
          } else {
            /* Compute parameter along edge (i, ip) */
            double t = (ex*dx + ey*dy) / denom;
            if (t > -eps && t < 1 + eps) {
              return PDM_POLYGON_INSIDE;
            }
          }
        }

      }

    }
  } // Loop over edges

  /* Result */
  if (wn == 0) {
    return PDM_POLYGON_OUTSIDE;
  } else {
    return PDM_POLYGON_INSIDE;
  }

}



/**
 * Warning : unstable function (Rather use PDM_polygon_point_in )
 *
 * \brief Test if a point is inside a 3d polygon using the Winding Number method
 *        (see http://geomalgorithms.com/a03-_inclusion.html)
 *
 * \param [in]  xyz           Point (x,y,z)-coordinates
 * \param [in]  n_vtx         Number of polygon vertices
 * \param [in]  vtx_xyz       Polygon vertices (x,y,z)-coordinates
 * \param [in]  char_length   Characteristic length (used to scale tolerance)
 * \param [in]  bounds        Bounds (xmin, xmax, ymin, ymax, zmin, zmax)
 *
 * \return      \ref Status inside, outside or degenerated
 *
 */

PDM_polygon_status_t PDM_polygon_point_in_3d
(
 const double  xyz[3],
 const int     n_vtx,
 const double *vtx_xyz,
 // const double  char_length,
 double        bounds[6],
 double        normal[3]
 )
{
  /*
   *  Do a quick bounds check
   */
  if (bounds != NULL) {
    if ( xyz[0] < bounds[0] || xyz[0] > bounds[1] ||
         xyz[1] < bounds[2] || xyz[1] > bounds[3] ||
         xyz[2] < bounds[4] || xyz[2] > bounds[5]) {
      return PDM_POLYGON_OUTSIDE;
    }
  }


  /*
   *  Projection onto median plane
   */
  double *xy = malloc (sizeof(double) * (n_vtx + 1) * 2);
  double *p = xy + 2*n_vtx;
  PDM_polygon_3d_to_2d (n_vtx,
                        vtx_xyz,
                        xy,
                        1,
                        xyz,
                        p,
                        normal);

  /*
   *  Perform 2D point-in-polygon test
   */

  //double char_length = 1e-6;

  PDM_polygon_status_t stat = PDM_polygon_point_in_2d (p,
                                                       n_vtx,
                                                       xy,
                                                       //char_length,
                                                       NULL);

  free (xy);
  return stat;
}


/**
 * \brief Compute parametric coordinates of a polygon's vertices
 * and a set of point inside the polygon's median plane.
 *
 * The uv-coordinate system is defined as follows:
 *   - origin at first polygon vertex ;
 *   - u-axis oriented from first to second polygon vertex.
 * Aspect-ratio and scale is preserved by the projection (no normalization).
 *
 * \param [in]  n_vtx    Number of polygon vertices
 * \param [in]  vtx_xyz  xyz-coordinates of polygon vertices (size = 3 * \ref n_vtx)
 * \param [out] vtx_uv   uv-coordinates of polygon vertices (size = 2 * \ref n_vtx)
 * \param [in]  n_pts    Number of points
 * \param [in]  pts_xyz  xyz-coordinates of points (size = 3 * \ref n_pts)
 * \param [out] pts_uv   uv-coordinates of points (size = 2 * \ref n_pts)
 * \param [in]  normal   Polygon's normal (size = 3 or NULL)
 *
 * \return      0 if the polygon is degenerate, 1 else.
 *
 */

int PDM_polygon_3d_to_2d
(
 const int    n_vtx,
 const double vtx_xyz[],
 double       vtx_uv[],
 const int    n_pts,
 const double pts_xyz[],
 double       pts_uv[],
 double       normal[3]
 )
{
  int idim, i;
  double *_normal = NULL;

  /* Compute polygon barycenter */
  double origin_xyz[3];
  PDM_polygon_compute_barycenter (n_vtx, vtx_xyz, origin_xyz);

  /* Build a suitable orthonormal frame */
  double tangent_u[3] = {0., 0., 0.}, tangent_v[3];
  if (normal == NULL) {
    _normal = malloc (sizeof(double) * 3);
    PDM_plane_normal (n_vtx,
                      vtx_xyz,
                      _normal);
  } else {
    _normal = normal;
  }

  /*   First tangent direction */
  int imin = -1;
  double nmin = HUGE_VAL;
  for (idim = 0; idim < 3; idim++) {
    if (_normal[idim] < 0. && nmin > -_normal[idim]) {
      imin = idim;
      nmin = -_normal[idim];
    } else if (_normal[idim] >= 0. && nmin > _normal[idim]) {
      imin = idim;
      nmin = _normal[idim];
    }
  }

  tangent_u[imin] = 1.;
  double mag = 0;
  for (idim = 0; idim < 3; idim++) {
    tangent_u[idim] -= _normal[imin] * _normal[idim];
    mag += tangent_u[idim] * tangent_u[idim];
  }

  if (mag < 1e-15) {
    return 0;
  }

  mag = 1. / sqrt(mag);
  for (idim = 0; idim < 3; idim++) {
    tangent_u[idim] *= mag;
  }

  /*   Second tangent direction */
  PDM_CROSS_PRODUCT (tangent_v, _normal, tangent_u);


  /* Compute coordinates in this frame */
  for (i = 0; i < n_vtx; i++) {
    double *uv = vtx_uv + 2*i;
    uv[0] = 0.;
    uv[1] = 0.;

    for (idim = 0; idim < 3; idim++) {
      double d = vtx_xyz[3*i + idim] - origin_xyz[idim];

      uv[0] += d * tangent_u[idim];
      uv[1] += d * tangent_v[idim];
    }
  }

  if (pts_xyz != NULL) {
    for (i = 0; i < n_pts; i++) {
      double *uv = pts_uv + 2*i;
      uv[0] = 0.;
      uv[1] = 0.;

      for (idim = 0; idim < 3; idim++) {
        double d = pts_xyz[3*i + idim] - origin_xyz[idim];

        uv[0] += d * tangent_u[idim];
        uv[1] += d * tangent_v[idim];
      }
    }
  }

  if (normal == NULL) {
    free (_normal);
  }

  return 1;
}








PDM_GCC_SUPPRESS_WARNING_WITH_PUSH("-Wfloat-equal")
inline static int _a_eq_b (const double a,
                           const double b) {
  return a == b;
  const double eps = 1e-15;
  return PDM_ABS (a - b) < eps;
}


inline static int _crossing (const double piy,
                             const double pipy,
                             const double ry) {
  return ((piy < ry) != (pipy < ry));
}

inline static int _right_crossing (const double det,
                                   const double piy,
                                   const double pipy)
{
  return ((det > 0) == (pipy > piy));
}

inline static void _modify_w (const double  piy,
                              const double  pipy,
                              int          *w) {
  *w += 2*(pipy > piy) - 1;
}
PDM_GCC_SUPPRESS_WARNING_POP

/**
 *  Algorithm 7 from
 *  "The point in polygon problem for arbitrary polygons", K.Hormann, A. Agathos (2001)
 **/

PDM_polygon_status_t
PDM_polygon_point_in_new
(
 const double  xyz[3],
 const int     n_vtx,
 const double *vtx_xyz,
 double        bounds[6],
 double        normal[3]
 )
{
  const double eps = 1e-15;

  /*
   *  Do a quick bounds check
   */
  if (bounds != NULL) {
    if (xyz[0] < bounds[0] || xyz[0] > bounds[1] ||
        xyz[1] < bounds[2] || xyz[1] > bounds[3] ||
        xyz[2] < bounds[4] || xyz[2] > bounds[5]) {
      return PDM_POLYGON_OUTSIDE;
    }
  }

  int imax = 0;
  double nmax = 0.;
  for (int i = 0; i < 3; i++) {
    double ni = PDM_ABS (normal[i]);
    if (ni > nmax) {
      imax = i;
      nmax = ni;
    }
  }
  if (nmax < eps) {
    return PDM_POLYGON_DEGENERATED;
  }

  int i1 = (imax + 1)%3;
  int i2 = (i1 + 1)%3;

  const double x = xyz[i1];
  const double y = xyz[i2];

  /*
   *  Perform 2d point-in-polygon test
   */
  const double *xi = vtx_xyz + i1;
  const double *yi = vtx_xyz + i2;

  if (_a_eq_b(*yi, y) && _a_eq_b(*xi, x)) {
    /* Point on vertex */
    return PDM_POLYGON_INSIDE;
  }

  const double *xj;
  const double *yj;
  int w = 0;

  for (int i = 0; i < n_vtx; i++) {
    int j = (i + 1)%n_vtx;

    xj = vtx_xyz + 3*j + i1;
    yj = vtx_xyz + 3*j + i2;

    if (_a_eq_b(*yj, y)) {
      if (_a_eq_b(*xj, x)) {
        /* Point on vertex */
        return PDM_POLYGON_INSIDE;
      } else {
        if (_a_eq_b(*yi, y) &&
            ((*xj > x) == (*xi < x))) {
          /* Point on edge */
          return PDM_POLYGON_INSIDE;
        }
      }
    }

    if (_crossing(*yi, *yj, y)) {
      //double det = (*xi - x)*(*yj - y) - (*yi - y)*(*xj - x);
      /*double a[2] = {x, y};
      double b[2] = {*xi, *yi};
      double c[2] = {*xj, *yj};*/
      double a[2], b[2], c[2];
      a[0] = x; a[1] = y; b[0] = *xi; b[1] = *yi; c[0] = *xj; c[1] = *yj;
      double det = PDM_predicate_orient2d (a, b, c);
      if (PDM_ABS(det) < 0) {//eps) {
        /* Point on edge */
        return PDM_POLYGON_INSIDE;
      }

      if (*xi >= x) {
        if (*xj > x) {
          _modify_w (*yi, *yj, &w);
        } else {
          if (_right_crossing(det, *yi, *yj)) {
            _modify_w (*yi, *yj, &w);
          }
        }
      } else {
        if (*xj > x) {
          if (_right_crossing(det, *yi, *yj)) {
            _modify_w (*yi, *yj, &w);
          }
        }
      }
    }


    xi += 3;
    yi += 3;
  }

  if (w == 0) {
    return PDM_POLYGON_OUTSIDE;
  } else {
    return PDM_POLYGON_INSIDE;
  }
}


/**
 * \brief Compute intersection point between a polygon and a semi-infinite ray
 *
 * \param[in]  ray_origin        Ray origin
 * \param[in]  ray_direction     Ray direction (need not be normalized)
 * \param[in]  n_vtx             Number of vertices
 * \param[in]  vtx_coord         Coordinates of the polygon's vertices
 * \param[in]  poly_center       Coordinates of the polygon's centroid (or NULL)
 * \param[in]  poly_normal       Normal to polygon's median plane (or NULL)
 * \param[in]  poly_bound        Polygon's bounds ([xmin, xmax, ymin, ymax, zmin, zmax] or NULL)
 * \param[out] intersection      Coordinates of the intersection point
 * \param[out] t                 Ray-parameter of the intersection point
 * \param[out] weight            Barycentric coordinates in polygon of intersection point (or NULL)
 *
 * \return Intersection status
 *
 */

PDM_polygon_status_t
PDM_polygon_ray_intersection
(
 const double  ray_origin[3],
 const double  ray_direction[3],
 const int     n_vtx,
 const double *vtx_coord,
       double *poly_center,
       double *poly_normal,
       double *poly_bound,
       double  intersection[3],
       double *t,
       double *weight
 )
{
  const double epsilon = 1e-12;

  double __poly_center[3], __poly_normal[3];
  double *_poly_center = poly_center;
  double *_poly_normal = poly_normal;

  if (_poly_center == NULL || _poly_normal == NULL) {
    int poly_vtx_idx[2] = {0, n_vtx};
    int poly_vtx[n_vtx];
    for (int i = 0; i < n_vtx; i++) {
      poly_vtx[i] = i+1;
    }
    int is_degenerate = 0;
    PDM_geom_elem_polygon_properties(1,
                                     poly_vtx_idx,
                                     poly_vtx,
                                     vtx_coord,
                                     __poly_normal,
                                     __poly_center,
                                     NULL,
                                     &is_degenerate);

    if (is_degenerate) {
      return PDM_POLYGON_DEGENERATED;
    }
  }

  if (_poly_center == NULL) {
    _poly_center = __poly_center;
  }

  if (_poly_normal == NULL) {
    _poly_normal = __poly_normal;
  }


  // intersect ray with polygon's median plane
  double vec[3] = {
    _poly_center[0] - ray_origin[0],
    _poly_center[1] - ray_origin[1],
    _poly_center[2] - ray_origin[2]
  };

  double denom = PDM_DOT_PRODUCT(ray_direction, _poly_normal);
  double numer = PDM_DOT_PRODUCT(vec,           _poly_normal);

  if (PDM_ABS(denom) < epsilon) {   // Epsilon is here to avoid division by 0
    if (PDM_ABS(numer) < epsilon) { // The ray is inside the polygon's median plane
      // 1) Check if ray origin is inside polygon
      PDM_polygon_status_t orig_in_poly = PDM_polygon_point_in_new(ray_origin,
                                                                   n_vtx,
                                                                   vtx_coord,
                                                                   poly_bound,
                                                                   _poly_normal);

      if (orig_in_poly != PDM_POLYGON_OUTSIDE) {

        if (orig_in_poly == PDM_POLYGON_INSIDE) {
          memcpy(intersection, ray_origin, sizeof(double) * 3);
          *t = 0.;

          if (weight != NULL) {
            PDM_mean_values_polygon_3d(n_vtx,
                                       vtx_coord,
                                       1,
                                       intersection,
                                       weight);
          }
        }

        return orig_in_poly;
      }

      // 2) Find first intersection between ray and the polygon's edges
      double ray_destination[3] = {
        ray_origin[0] + ray_direction[0],
        ray_origin[1] + ray_direction[1],
        ray_origin[2] + ray_direction[2]
      };

      PDM_polygon_status_t stat = PDM_POLYGON_OUTSIDE;
      *t = HUGE_VAL;
      double s, _s, _t;
      int iedge = -1;
      for (int i = 0; i < n_vtx; i++) {
        PDM_line_intersection_mean_square(vtx_coord + 3*i,
                                          vtx_coord + 3*((i+1)%n_vtx),
                                          ray_origin,
                                          ray_destination,
                                          &_s,
                                          &_t);
        if (_s >= 0. && _s <= 1. && _t >= 0.) {
          stat = PDM_POLYGON_INSIDE;
          if (_t < *t) {
            iedge = i;
            s  = _s;
            *t = _t;
          }
        }
      }

      if (stat == PDM_POLYGON_INSIDE) {
        intersection[0] = ray_origin[0] + (*t)*ray_direction[0];
        intersection[1] = ray_origin[1] + (*t)*ray_direction[1];
        intersection[2] = ray_origin[2] + (*t)*ray_direction[2];

        if (weight != NULL) {
          for (int i = 0; i < n_vtx; i++) {
            weight[i] = 0;
          }
          weight[iedge]           = 1. - s;
          weight[(iedge+1)%n_vtx] =      s;
        }
      }

      return stat;


    } else { // The ray is parallel but not coplanar
      return PDM_POLYGON_OUTSIDE;
    }
  }
  else { // General case
    *t = numer/denom;

    if (*t > 0) {
      intersection[0] = ray_origin[0] + (*t)*ray_direction[0];
      intersection[1] = ray_origin[1] + (*t)*ray_direction[1];
      intersection[2] = ray_origin[2] + (*t)*ray_direction[2];
    }
  }


  /* We found an intersection point, now check if it is inside the polygon */
  // double poly_bound[6] = {
  //   HUGE_VAL, -HUGE_VAL,
  //   HUGE_VAL, -HUGE_VAL,
  //   HUGE_VAL, -HUGE_VAL
  // };
  // for (int i = 0; i < poly_vtx_n; i++) {
  //   double *vc = vtx_coord + 3*i;
  //   for (int j = 0; j < 3; j++) {
  //     poly_bound[2*j  ] = PDM_MIN(poly_bound[2*j  ], vc[j]);
  //     poly_bound[2*j+1] = PDM_MAX(poly_bound[2*j+1], vc[j]);
  //   }
  // }

  // // Inflate the face's bounding box
  // double d = 0.;
  // for (int j = 0; j < 3; j++) {
  //   d += (face_bound[2*j+1] - face_bound[2*j])*(face_bound[2*j+1] - face_bound[2*j]);
  // }
  // d = 0.1*sqrt(d);
  // for (int j = 0; j < 3; j++) {
  //   face_bound[2*j  ] -= d;
  //   face_bound[2*j+1] += d;
  // }

  // PDM_log_trace_array_double(intersection, 3, "intersection : ");

  PDM_polygon_status_t in_poly = PDM_polygon_point_in_new(intersection,
                                                          n_vtx,
                                                          vtx_coord,
                                                          poly_bound,
                                                          _poly_normal);

  if (in_poly == PDM_POLYGON_INSIDE && weight != NULL) {
    PDM_mean_values_polygon_3d(n_vtx,
                               vtx_coord,
                               1,
                               intersection,
                               weight);
  }

  return in_poly;
}

#ifdef __cplusplus
}
#endif /* __cplusplus */
