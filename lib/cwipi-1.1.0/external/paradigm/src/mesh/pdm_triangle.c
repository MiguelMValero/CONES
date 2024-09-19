/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "pdm.h"
#include "pdm_priv.h"
#include "pdm_triangle.h"
#include "pdm_line.h"
#include "pdm_plane.h"
#include "pdm_logging.h"

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
 * Private function definitions
 *============================================================================*/

/* Computes the determinant of a 3x3 matrix defined by its columns */

static inline double
_determinant_3x3
(
 const double a[3],
 const double b[3],
 const double c[3]
 )
{
  return a[0] * (b[1]*c[2] - b[2]*c[1])
    +    a[1] * (b[2]*c[0] - b[0]*c[2])
    +    a[2] * (b[0]*c[1] - b[1]*c[0]);
}

/*============================================================================
 * Public function definitions
 *============================================================================*/


/**
 * \brief Evaluates the position in a triangle
 *
 * \param [in]  x        Point coordinates to evaluate position
 * \param [in]  pts      Triangle vertices coordinates
 * \param [out] closest  Closest Point in Triangle or NULL
 * \param [out] min_dist2 Square of the distance
 * \param [out] weights  Vertices weights or NULL
 *
 * \return      \ref PDM_TRIANGLE_INSIDE or \ref PDM_TRIANGLE_OUTSIDE
 *              if the projected is in the triangle or not
 *
 */

/*  This function is derived from VTK                                      */
/*  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen                */
/*  All rights reserved.                                                   */
/*  See Copyright.txt or http://www.kitware.com/Copyright.htm for details. */


PDM_triangle_status_t
PDM_triangle_evaluate_position_old
(
 const double  x[3],
 const double  pts[9],
       double *closest_point,
       double *min_dist2,
       double *weights
)
{
  double *pt1, *pt2, *pt3;
  double n[3], fabsn;
  double rhs[2], c1[2], c2[2];
  double det;
  double maxComponent;
  int idx=0, indices[2];
  double dist_to_point, dist_to_line1, dist_to_line2;
  double *closest, closest_point1[3], closest_point2[3], cp[3];
  double pcoords[3];

  pcoords[2] = 0.0;

  double weights_local[3];
  double *_weights = weights_local;
  if (weights != NULL) {
    _weights = weights;
  }

  /*
   * Get normal for triangle, only the normal direction is needed, i.e. the
   * normal need not be normalized (unit length)
   */

  PDM_plane_normal (3, pts, n);

  pt1 = (double *) pts;
  pt2 = (double *) pts + 3;
  pt3 = (double *) pts + 6;

  /*
   * Project point to plane
   */

  PDM_plane_projection2 (x, pt1, n, cp);

  /*
   * Construct matrices.  Since we have over determined system, need to find
   * which 2 out of 3 equations to use to develop equations. (Any 2 should
   * work since we've projected point to plane.)
   */

  maxComponent = 0.0;
  for (int i = 0; i < 3; i++) {

    /*
     * Trying to avoid an expensive call to fabs()
     */

    if (n[i] < 0) {
      fabsn = -n[i];
    }
    else {
      fabsn = n[i];
    }
    if (fabsn > maxComponent) {
      maxComponent = fabsn;
      idx = i;
    }
  }

  for (int j = 0, i = 0; i < 3; i++) {
    if (i != idx) {
      indices[j++] = i;
    }
  }

  for (int i = 0; i < 2; i++) {
    // rhs[i] = cp[indices[i]] - pt3[indices[i]];
    // c1[i] = pt1[indices[i]] - pt3[indices[i]];
    // c2[i] = pt2[indices[i]] - pt3[indices[i]];
    rhs[i] = cp[indices[i]] - pt1[indices[i]];
    c1[i] = pt2[indices[i]] - pt1[indices[i]];
    c2[i] = pt3[indices[i]] - pt1[indices[i]];
  }

PDM_GCC_SUPPRESS_WARNING_WITH_PUSH("-Wfloat-equal")

  if ((det = PDM_DETERMINANT2X2(c1,c2)) == 0.0) {
    pcoords[0] = 0.0;
    pcoords[1] = 0.0;
    return PDM_TRIANGLE_DEGENERATED;
  }

PDM_GCC_SUPPRESS_WARNING_POP

  pcoords[0] = PDM_DETERMINANT2X2(rhs,c2) / det;
  pcoords[1] = PDM_DETERMINANT2X2(c1,rhs) / det;

  /*
   * Okay, now find closest point to element
   */

  _weights[0] = 1 - (pcoords[0] + pcoords[1]);
  _weights[1] = pcoords[0];
  _weights[2] = pcoords[1];

  if ( _weights[0] >= 0.0 && _weights[0] <= 1.0 &&
       _weights[1] >= 0.0 && _weights[1] <= 1.0 &&
       _weights[2] >= 0.0 && _weights[2] <= 1.0 ) {

    /*
     * Projection distance
     */

    if (closest_point) {
      double v_cp_x[3];
      for (int i = 0; i < 3; i++) {
        v_cp_x[i] = cp[i] - x[i];
      }

      *min_dist2 = PDM_DOT_PRODUCT(v_cp_x, v_cp_x);
      closest_point[0] = cp[0];
      closest_point[1] = cp[1];
      closest_point[2] = cp[2];
    }
    return PDM_TRIANGLE_INSIDE;
  }

  else {
    double t1, t2;
    if (closest_point) {
      if (_weights[1] < 0.0 && _weights[2] < 0.0) {
        double v_pt3_x[3];
        for (int i = 0; i < 3; i++) {
          v_pt3_x[i] = pt3[i] - x[i];
        }
        dist_to_point = PDM_DOT_PRODUCT (v_pt3_x, v_pt3_x);
        dist_to_line1 = PDM_line_distance (x, pt1, pt3, &t1, closest_point1);
        dist_to_line2 = PDM_line_distance (x, pt3, pt2, &t2, closest_point2);
        if (dist_to_point < dist_to_line1) {
          *min_dist2 = dist_to_point;
          closest = pt3;
          _weights[0] = 0.;
          _weights[1] = 0.;
          _weights[2] = 1.;
        }
        else {
          *min_dist2 = dist_to_line1;
          closest = closest_point1;
          if (t1 < 0.) {
            t1 = 0.;
          } else if (t1 > 1.) {
            t1 = 1.;
          }
          _weights[0] = 1. - t1;
          _weights[1] = 0.;
          _weights[2] = t1;
        }
        if (dist_to_line2 < *min_dist2) {
          *min_dist2 = dist_to_line2;
          closest = closest_point2;
          if (t2 < 0.) {
            t2 = 0.;
          } else if (t2 > 1.) {
            t2 = 1.;
          }
          _weights[0] = 0.;
          _weights[1] = t2;
          _weights[2] = 1. - t2;
        }
        for (int i = 0; i < 3; i++) {
          closest_point[i] = closest[i];
        }

      }
      else if (_weights[2] < 0.0 && _weights[0] < 0.0) {
        double v_pt1_x[3];
        for (int i = 0; i < 3; i++) {
          v_pt1_x[i] = pt1[i] - x[i];
        }
        dist_to_point = PDM_DOT_PRODUCT(v_pt1_x, v_pt1_x);
        dist_to_line1 = PDM_line_distance (x, pt1, pt3, &t1, closest_point1);
        dist_to_line2 = PDM_line_distance (x, pt1, pt2, &t2, closest_point2);
        if (dist_to_point < dist_to_line1) {
          *min_dist2 = dist_to_point;
          closest = pt1;
          _weights[0] = 1.;
          _weights[1] = 0.;
          _weights[2] = 0.;
        }
        else {
          *min_dist2 = dist_to_line1;
          closest = closest_point1;
          if (t1 < 0.) {
            t1 = 0.;
          } else if (t1 > 1.) {
            t1 = 1.;
          }
          _weights[0] = 1. - t1;
          _weights[1] = 0.;
          _weights[2] = t1;
        }
        if (dist_to_line2 < *min_dist2) {
          *min_dist2 = dist_to_line2;
          closest = closest_point2;
           if (t2 < 0.) {
            t2 = 0.;
          } else if (t2 > 1.) {
            t2 = 1.;
          }
          _weights[0] = 1. - t2;
          _weights[1] = t2;
          _weights[2] = 0.;
        }
        for (int i = 0; i < 3; i++) {
          closest_point[i] = closest[i];
        }

      }
      else if ( _weights[1] < 0.0 && _weights[0] < 0.0 ) {
        double v_pt2_x[3];
        for (int i = 0; i < 3; i++) {
          v_pt2_x[i] = pt2[i] - x[i];
        }
        dist_to_point = PDM_DOT_PRODUCT (v_pt2_x, v_pt2_x);
        dist_to_line1 = PDM_line_distance (x, pt2, pt3, &t1, closest_point1);
        dist_to_line2 = PDM_line_distance (x, pt1, pt2, &t2, closest_point2);
        if (dist_to_point < dist_to_line1) {
          *min_dist2 = dist_to_point;
          closest = pt2;
          _weights[0] = 0.;
          _weights[1] = 1.;
          _weights[2] = 0.;
        }
        else {
          *min_dist2 = dist_to_line1;
          closest = closest_point1;
          if (t1 < 0.) {
            t1 = 0.;
          } else if (t1 > 1.) {
            t1 = 1.;
          }
          _weights[0] = 0.;
          _weights[1] = 1. - t1;
          _weights[2] = t1;
        }
        if (dist_to_line2 < *min_dist2) {
          *min_dist2 = dist_to_line2;
          closest = closest_point2;
           if (t2 < 0.) {
            t2 = 0.;
          } else if (t2 > 1.) {
            t2 = 1.;
          }
          _weights[0] = 1. - t2;
          _weights[1] = t2;
          _weights[2] = 0.;
        }
        for (int i = 0; i < 3; i++) {
          closest_point[i] = closest[i];
        }

      }
      else if (_weights[0] < 0.0) {
        *min_dist2 = PDM_line_distance (x, pt1, pt2, &t1, closest_point);
        if (t1 < 0.) {
            t1 = 0.;
          } else if (t1 > 1.) {
            t1 = 1.;
          }
          _weights[0] = 1. - t1;
          _weights[1] = t1;
          _weights[2] = 0.;
      }
      else if (_weights[1] < 0.0) {
          *min_dist2 = PDM_line_distance (x, pt2, pt3, &t1, closest_point);
          if (t1 < 0.) {
            t1 = 0.;
          } else if (t1 > 1.) {
            t1 = 1.;
          }
          _weights[0] = 0.;
          _weights[1] = 1. - t1;
          _weights[2] = t1;
        }
      else if (_weights[2] < 0.0) {
          *min_dist2 = PDM_line_distance (x, pt1, pt3, &t1, closest_point);
          if (t1 < 0.) {
            t1 = 0.;
          } else if (t1 > 1.) {
            t1 = 1.;
          }
          _weights[0] = 1. - t1;
          _weights[1] = 0.;
          _weights[2] = t1;
        }

        }
    return PDM_TRIANGLE_OUTSIDE;
  }
}

/**
 * \brief Computes the intersection between a line and a triangle
 *
 * \param [in]   line        Points of the line
 * \param [in]   tria_coord  Points of the triangle
 * \param [out]  ip          Intersection point
 *
 */

PDM_triangle_status_t
PDM_triangle_line_intersection
(
const double line[6],
const double tria_coord[9],
      double ip[3]
)
{

  PDM_plane_line_intersection(line, tria_coord, ip);

  PDM_triangle_status_t found = PDM_TRIANGLE_OUTSIDE;

  /* Check if intersection point is in the triangle */
  double c1[3];
  double c2[3];
  double c3[3];

  double det1, det2, det3;

  // ABD
  c1[0] = tria_coord[0]; c1[1] = tria_coord[3]; c1[2] = ip[0];
  c2[0] = tria_coord[1]; c2[1] = tria_coord[4]; c2[2] = ip[1];
  c3[0] = tria_coord[2]; c3[1] = tria_coord[5]; c3[2] = ip[2];

  det1 = _determinant_3x3(c1, c2, c3);

  // DBC
  c1[0] = ip[0]; c1[1] = tria_coord[3]; c1[2] = tria_coord[6];
  c2[0] = ip[1]; c2[1] = tria_coord[4]; c2[2] = tria_coord[7];
  c3[0] = ip[2]; c3[1] = tria_coord[5]; c3[2] = tria_coord[8];

  det2 = _determinant_3x3(c1, c2, c3);

  // ADC
  c1[0] = tria_coord[0]; c1[1] = ip[0]; c1[2] = tria_coord[6];
  c2[0] = tria_coord[1]; c2[1] = ip[1]; c2[2] = tria_coord[7];
  c3[0] = tria_coord[2]; c3[1] = ip[2]; c3[2] = tria_coord[8];

  det3 = _determinant_3x3(c1, c2, c3);

  if (det1 > 0 && det2 > 0 && det3 > 0) {

    found = PDM_TRIANGLE_INSIDE;

  } // end if point is inside triangle

  return found;

}


/**
 * \brief Computes triangle barycenter
 *
 * \param [in]   pts     Triangle vertices coordinates
 * \param [out]  bary    Barycenter
 *
 */

void
PDM_triangle_compute_barycenter
(
 const double pts[9],
       double bary[3]
)
{
  bary[0] = 0.;
  bary[1] = 0.;
  bary[2] = 0.;

  for (int i = 0; i < 3; i++) {
    for (int ipt = 0; ipt < 3; ipt++) {
      bary[i] += pts[3*ipt+i];
    }
    bary[i] /= 3;
  }
}



/**
 * \brief Computes closest point on a triangle
 *
 * \param [in]   p              Point coordinates
 * \param [in]   v              Triangle vertices coordinates
 * \param [out]  closest_point  Closest point coordinates
 * \param [out]  min_dist2      Squared distance from point to triangle
 * \param [out]  weights        Barycentric coordinates of closest point (or NULL)
 *
 */

PDM_triangle_status_t
PDM_triangle_closest_point
(
 const double  x[3],
 const double  v[9],
 double       *closest_point,
 double       *min_dist2,
 double       *weights
 )
{
  const double *v0 = v;
  const double *v1 = v + 3;
  const double *v2 = v + 6;

  // double v20[3] = {v0[0] - v2[0],
  //                  v0[1] - v2[1],
  //                  v0[2] - v2[2]};

  // double v21[3] = {v1[0] - v2[0],
  //                  v1[1] - v2[1],
  //                  v1[2] - v2[2]};

  // double v2x[3] = {x[0] - v2[0],
  //                  x[1] - v2[1],
  //                  x[2] - v2[2]};

  // double a = PDM_DOT_PRODUCT (v20, v20);
  // double b = PDM_DOT_PRODUCT (v20, v21);
  // double c = PDM_DOT_PRODUCT (v21, v21);
  double e1[3] = {v1[0] - v0[0],
                  v1[1] - v0[1],
                  v1[2] - v0[2]};
  double e2[3] = {v2[0] - v0[0],
                  v2[1] - v0[1],
                  v2[2] - v0[2]};
  double e[3]  = {x[0] - v0[0],
                  x[1] - v0[1],
                  x[2] - v0[2]};

  double a = PDM_DOT_PRODUCT(e1, e1);
  double b = PDM_DOT_PRODUCT(e1, e2);
  double c = PDM_DOT_PRODUCT(e2, e2);

  double det = a*c - b*b;

  if (det < 1e-14) {
    return PDM_TRIANGLE_DEGENERATED;
  }

  // double r = PDM_DOT_PRODUCT (v20, v2x);
  // double s = PDM_DOT_PRODUCT (v21, v2x);
  double r = PDM_DOT_PRODUCT(e1, e);
  double s = PDM_DOT_PRODUCT(e2, e);


  /* Solve for weights of orthogonal projection of point on triangle's plane */
  double weights_local[3];
  double *_weights = weights_local;
  if (weights != NULL) {
    _weights = weights;
  }
  _weights[1] = (r*c - s*b) / det;
  _weights[2] = (s*a - r*b) / det;
  // log_trace("  u = %f, v = %f\n", _weights[1], _weights[2]);
  _weights[0] = 1. - _weights[1] - _weights[2];

  /* Projection inside triangle (= closest point) */
  if ( _weights[0] >= 0.0 && _weights[0] <= 1.0 &&
       _weights[1] >= 0.0 && _weights[1] <= 1.0 &&
       _weights[2] >= 0.0 && _weights[2] <= 1.0 ) {

    *min_dist2 = 0.;
    for (int idim = 0; idim < 3; idim++) {
      // closest_point[idim] = v[6 + idim] + _weights[1]*v20[idim] + _weights[2]*v21[idim];
      closest_point[idim] = v[idim] + _weights[1]*e1[idim] + _weights[2]*e2[idim];
      double delta = x[idim] - closest_point[idim];
      *min_dist2 += delta * delta;
    }

    return PDM_TRIANGLE_INSIDE;
  }

  /* Projection inside triangle --> find closest point on triangle's boundary */
  else {

    double t01, t12, t20, d01, d12, d20, c01[3], c12[3], c20[3];

    int i, j, k;
    double t;
    d01 = PDM_line_distance (x, v0, v1, &t01, c01);
    d12 = PDM_line_distance (x, v1, v2, &t12, c12);
    d20 = PDM_line_distance (x, v2, v0, &t20, c20);
    // log_trace("d01/12/20 = %f / %f / %f\n", d01, d12, d20);

    if (d01 <= d12 && d01 <= d20) {
      // i = 0;
      // j = 1;
      // k = 2;
      i = 2; j = 0; k = 1;
      *min_dist2 = d01;
      t = t01;

      closest_point[0] = c01[0];
      closest_point[1] = c01[1];
      closest_point[2] = c01[2];
    }

    else if (d12 <= d01 && d12 <= d20) {
      // i = 1;
      // j = 2;
      // k = 0;
      i = 0; j = 1; k = 2;
      *min_dist2 = d12;
      t = t12;

      closest_point[0] = c12[0];
      closest_point[1] = c12[1];
      closest_point[2] = c12[2];
    }

    else {
      // i = 2;
      // j = 0;
      // k = 1;
      i = 1; j = 0; k = 2;
      *min_dist2 = d20;
      t = t20;

      closest_point[0] = c20[0];
      closest_point[1] = c20[1];
      closest_point[2] = c20[2];
    }

    if (t < 0.) {
      t = 0.;
    } else if (t > 1.) {
      t = 1.;
    }

    PDM_UNUSED(i);
    PDM_UNUSED(j);
    PDM_UNUSED(k);
    // _weights[i] = 0.;
    // _weights[j] = 1. - t;
    // _weights[k] = t;

    return PDM_TRIANGLE_OUTSIDE;
  }

}



/**
 * \brief Computes the center and radius of a triangle's circumcircle
 *
 * \param [in]   vtx_coord  Triangle vertices coordinates
 * \param [out]  center     Circumcircle center
 * \param [out]  radius     Circumcircle radius
 *
 */

void
PDM_triangle_circumcircle
(
 const double  vtx_coord[9],
 double        center[3],
 double       *radius
 )
{
  double a[3] = {vtx_coord[0] - vtx_coord[6],
                 vtx_coord[1] - vtx_coord[7],
                 vtx_coord[2] - vtx_coord[8]};

  double b[3] = {vtx_coord[3] - vtx_coord[6],
                 vtx_coord[4] - vtx_coord[7],
                 vtx_coord[5] - vtx_coord[8]};

  double ab[3] = {b[0] - a[0],
                  b[1] - a[1],
                  b[2] - a[2]};


  double axb[3];
  PDM_CROSS_PRODUCT (axb, a, b);

  double a2 = PDM_DOT_PRODUCT (a, a);
  double b2 = PDM_DOT_PRODUCT (b, b);
  double ab2 = PDM_DOT_PRODUCT (ab, ab);
  double axb2 = PDM_DOT_PRODUCT (axb, axb);

  PDM_GCC_SUPPRESS_WARNING_WITH_PUSH("-Wfloat-equal")

    if (axb2 == 0) {
      center[0] = vtx_coord[6];
      center[1] = vtx_coord[7];
      center[2] = vtx_coord[8];
      *radius = 0.;
    }

    else {
      axb2 = 0.25 / axb2;

      a[0] = a2*b[0] - b2*a[0];
      a[1] = a2*b[1] - b2*a[1];
      a[2] = a2*b[2] - b2*a[2];
      PDM_CROSS_PRODUCT (b, a, axb);

      center[0] = vtx_coord[6] + 2. * b[0] * axb2;
      center[1] = vtx_coord[7] + 2. * b[1] * axb2;
      center[2] = vtx_coord[8] + 2. * b[2] * axb2;
      *radius = sqrt(a2 * b2 * ab2 * axb2);
    }

  PDM_GCC_SUPPRESS_WARNING_POP
    }


/**
 * \brief Compute intersection point between a triangle and a semi-infinite ray
 *
 * \param[in]  origin        Ray origin
 * \param[in]  direction     Ray direction (need not be normalized)
 * \param[in]  tri_coord     Coordinates of the triangle's vertices
 * \param[out] intersection  Coordinates of the intersection point
 * \param[out] t             Ray-parameter of the intersection point
 * \param[out] weight        Barycentric coordinates in triangle of intersection point (or NULL)
 *
 * \return Intersection status
 *
 */

PDM_triangle_status_t
PDM_triangle_ray_intersection
(
 const double  origin[3],
 const double  direction[3],
 const double  tri_coord[9],
       double  intersection[3],
       double *t,
       double *weight
 )
{
  const double epsilon = 1e-12;

  double v01[3] = {
    tri_coord[3] - tri_coord[0],
    tri_coord[4] - tri_coord[1],
    tri_coord[5] - tri_coord[2]
  };

  double v02[3] = {
    tri_coord[6] - tri_coord[0],
    tri_coord[7] - tri_coord[1],
    tri_coord[8] - tri_coord[2]
  };

  double normal[3];
  PDM_CROSS_PRODUCT(normal, v01, v02);
  double det = PDM_DOT_PRODUCT(normal, normal);
  if (det < 1e-30) {
    return PDM_TRIANGLE_DEGENERATED;
  }
  // PDM_log_trace_array_double(normal, 3, "normal : ");

  double vec[3] = {
    tri_coord[0] - origin[0],
    tri_coord[1] - origin[1],
    tri_coord[2] - origin[2]
  };
  double denom = PDM_DOT_PRODUCT(direction, normal);
  double numer = PDM_DOT_PRODUCT(vec,       normal);

  // log_trace("denom = %e, numer = %e\n", denom, numer);

  if (PDM_ABS(denom) < epsilon) {
    // Ray parallel to plane

    if (PDM_ABS(numer) < epsilon) {
      // Ray inside plane

      // 1) Check if ray origin is inside triangle
      double a = PDM_DOT_PRODUCT(v01, v01);
      double b = PDM_DOT_PRODUCT(v01, v02);
      double c = PDM_DOT_PRODUCT(v02, v02);

      double e = -PDM_DOT_PRODUCT(vec, v01);
      double f = -PDM_DOT_PRODUCT(vec, v02);

      double u = e*c - f*b;
      double v = a*f - b*e;

      if (u >= 0. && v >= 0. && u + v <= det) {
        memcpy(intersection, origin, sizeof(double)*3);
        *t = 0.;
        if (weight != NULL) {
          double idet = 1./det;
          weight[1] = u*idet;
          weight[2] = v*idet;
          weight[0] = 1. - weight[1] - weight[2];
        }
        return PDM_TRIANGLE_INSIDE;
      }

      // 2) Find first intersection between ray and the triangle's edges
      double destination[3] = {
        origin[0] + direction[0],
        origin[1] + direction[1],
        origin[2] + direction[2]
      };

      PDM_triangle_status_t stat = PDM_TRIANGLE_OUTSIDE;
      *t = HUGE_VAL;
      double s, _s, _t;
      int iedge = -1;
      for (int i = 0; i < 3; i++) {
        PDM_line_intersection_mean_square(tri_coord + 3*i,
                                          tri_coord + 3*((i+1)%3),
                                          origin,
                                          destination,
                                          &_s,
                                          &_t);
        if (_s >= 0. && _s <= 1. && _t >= 0.) {
          stat = PDM_TRIANGLE_INSIDE;
          iedge = i;
          s = _s;
          *t = PDM_MIN(*t, _t);
        }
      }

      if (stat == PDM_TRIANGLE_INSIDE) {
        intersection[0] = origin[0] + (*t)*direction[0];
        intersection[1] = origin[1] + (*t)*direction[1];
        intersection[2] = origin[2] + (*t)*direction[2];

        if (weight != NULL) {
          if (iedge == 0) {
            weight[1] = s;
            weight[2] = 0;
          }
          else if (iedge == 1) {
            weight[1] = 1 - s;
            weight[2] = s;
          }
          else {
            weight[1] = 0;
            weight[2] = 1 - s;
          }
          weight[0] = 1. - weight[1] - weight[2];
        }
      }

      return stat;

    }
    else {
      return PDM_TRIANGLE_OUTSIDE;
    }

  }
  else {
    // General case
    *t = numer/denom;
    // log_trace("t = %f (%f %f %f)\n",
    //           t,
    //           origin[0] + t*direction[0],
    //           origin[1] + t*direction[1],
    //           origin[2] + t*direction[2]);

    if (*t < 0.) {
      return PDM_TRIANGLE_OUTSIDE;
    }
    else {
      // Check if ray-plane intersection is inside triangle
      intersection[0] = origin[0] + (*t)*direction[0];
      intersection[1] = origin[1] + (*t)*direction[1];
      intersection[2] = origin[2] + (*t)*direction[2];
      // log_trace("intersection = %f %f %f\n",
      //           intersection[0], intersection[1], intersection[2]);

      double a = PDM_DOT_PRODUCT(v01, v01);
      double b = PDM_DOT_PRODUCT(v01, v02);
      double c = PDM_DOT_PRODUCT(v02, v02);

      double vec2[3] = {
        intersection[0] - tri_coord[0],
        intersection[1] - tri_coord[1],
        intersection[2] - tri_coord[2]
      };

      double e = PDM_DOT_PRODUCT(vec2, v01);
      double f = PDM_DOT_PRODUCT(vec2, v02);

      double u = e*c - f*b;
      double v = a*f - b*e;

      if (weight != NULL) {
        double idet = 1./det;
        weight[1] = u*idet;
        weight[2] = v*idet;
        weight[0] = 1. - weight[1] - weight[2];
      }

      if (u >= 0. && v >= 0. && u + v <= det) {
        return PDM_TRIANGLE_INSIDE;
      }

    }

  }

  return PDM_TRIANGLE_OUTSIDE;
}


/* Compute squared euclidean distance between points a and b */
static inline double
_dist2_point_point
(
 const double a[3],
 const double b[3]
 )
 {
  double dist2 = 0;
  for (int i = 0; i < 3; i++) {
    dist2 += (a[i] - b[i])*(a[i] - b[i]);
  }
  return dist2;
 }

/* Compute squared euclidean distance between point p and segment ab */
static inline double
_dist2_point_segment
(
 const double  p[3],
 const double  a[3],
 const double  b[3],
       double *t
 )
{
  double ab[3] = {b[0] - a[0], b[1] - a[1], b[2] - a[2]};
  double abab = PDM_DOT_PRODUCT(ab, ab);

  if (abab <= 0.) {
    // degenerate segment
    double mid[3] = {0.5*(a[0] + b[0]), 0.5*(a[1] + b[1]), 0.5*(a[2] + b[2])};
    *t = 0.5;
    return _dist2_point_point(p, mid);
  }
  else {
    double ap[3] = {p[0] - a[0], p[1] - a[1], p[2] - a[2]};
    *t = PDM_DOT_PRODUCT(ap, ab) / abab;
    if (*t < 0) {
      *t = 0;
      return _dist2_point_point(p, a);
    }
    else if (*t > 1) {
      *t = 1;
      return _dist2_point_point(p, b);
    }
    else {
      double c[3] = {
        a[0] + (*t)*ab[0],
        a[1] + (*t)*ab[1],
        a[2] + (*t)*ab[2]
      };
      return _dist2_point_point(p, c);
    }
  }
}

PDM_triangle_status_t
PDM_triangle_evaluate_position
(
 const double  x[3],
 const double  pts[9],
       double *closest_point,
       double *min_dist2,
       double *weights
)
{
  PDM_triangle_status_t stat = PDM_TRIANGLE_OUTSIDE;

  double __weight[3];
  double *_weight = weights;
  if (weights == NULL) {
    _weight = __weight;
  }

  double v01[3];
  double v02[3];
  double v0x[3];
  for (int i = 0; i < 3; i++) {
    v01[i] = pts[i+3] - pts[i];
    v02[i] = pts[i+6] - pts[i];
    v0x[i] =   x[i  ] - pts[i];
  }

  double v01v01 = PDM_DOT_PRODUCT(v01, v01);
  double v01v02 = PDM_DOT_PRODUCT(v01, v02);
  double v02v02 = PDM_DOT_PRODUCT(v02, v02);

  double det = v01v01*v02v02 - v01v02*v01v02;

  if (det <= 0.) {
    return PDM_TRIANGLE_DEGENERATED;
  }

  double v0xv01 = PDM_DOT_PRODUCT(v0x, v01);
  double v0xv02 = PDM_DOT_PRODUCT(v0x, v02);

  _weight[1] = v0xv01*v02v02 - v0xv02*v01v02;
  _weight[2] = v0xv02*v01v01 - v0xv01*v01v02;
  _weight[0] = det - _weight[1] - _weight[2];

  if (_weight[1] >= 0 && _weight[1] <= det &&
      _weight[2] >= 0 && _weight[2] <= det &&
      _weight[0] >= 0) {
    // projection lies inside the triangle
    stat = PDM_TRIANGLE_INSIDE;

    double idet = 1./det;
    _weight[0] *= idet;
    _weight[1] *= idet;
    _weight[2] *= idet;


    *min_dist2 = 0;
    for (int i = 0; i < 3; i++) {
      double delta = x[i] - (_weight[0]*pts[i] + _weight[1]*pts[3+i] + _weight[2]*pts[6+i]);
      *min_dist2 += delta*delta;
    }
  }

  else {
    // projection lies outside the triangle
    // => find closest point on the triangle boundary
    *min_dist2 = HUGE_VAL;
    int    imin = -1;
    double tmin = 0.;
    for (int i = 0; i < 3; i++) {
      if (_weight[i] < 0) {
        double t;
        double dist2 = _dist2_point_segment(x, &pts[3*((i+1)%3)], &pts[3*((i+2)%3)], &t);

        if (dist2 < *min_dist2) {
          *min_dist2 = dist2;
          imin = (i+1)%3;
          tmin = t;
        }
      }
    }

    _weight[imin      ] = 1. - tmin;
    _weight[(imin+1)%3] = tmin;
    _weight[(imin+2)%3] = 0.;
  }

  if (closest_point != NULL) {
    for (int i = 0; i < 3; i++) {
      closest_point[i] = _weight[0]*pts[i] + _weight[1]*pts[3+i] + _weight[2]*pts[6+i];
    }
  }
  return stat;
}


/**
 * \brief Build triangle->vertex from triangle->edge and edge->vertex connectivities.
 *
 * \note In each triangle, edge #i is opposite to vertex #i:
 *         v2
 *         o
 *        / \
 *    e1 /   \ e0
 *      /     \
 *  v0 o-------o v1
 *         e2
 *
 * \param [in]  n_face     Number of faces
 * \param [in]  face_edge  Face -> edge (signed) connectivity (1-based, size : 3 * \p n_face)
 * \param [in]  face_edge  Edge -> vertex connectivity (1-based, size : 2 * *n_edge*)
 * \param [out] face_vtx   Face -> vertex (signed) connectivity (size : 3 * \p n_face)
 *
 */

void
PDM_triangle_ngon_to_nodal
(
 int   n_face,
 int  *face_edge,
 int  *edge_vtx,
 int **face_vtx
 )
{
  *face_vtx = malloc(sizeof(int) * n_face * 3);

  for (int iface = 0; iface < n_face; iface++) {

    int *fv = *face_vtx  + 3*iface;
    int *fe =  face_edge + 3*iface;

    int iedge = fe[0];

    if (iedge < 0) {
      iedge = -iedge - 1;
      fv[1] = edge_vtx[2*iedge+1];
      fv[2] = edge_vtx[2*iedge  ];
    } else {
      iedge =  iedge - 1;
      fv[1] = edge_vtx[2*iedge  ];
      fv[2] = edge_vtx[2*iedge+1];
    }

    iedge = PDM_ABS(fe[1]) - 1;
    int ivtx1 = edge_vtx[2*iedge  ];
    int ivtx2 = edge_vtx[2*iedge+1];
    if (ivtx1 == fv[1] || ivtx1 == fv[2]) {
      fv[0] = ivtx2;
    } else {
      assert(ivtx2 == fv[1] || ivtx2 == fv[2]);
      fv[0] = ivtx1;
    }

  }
}


#ifdef __cplusplus
}
#endif /* __cplusplus */
