/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/
#include <string.h>

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "pdm.h"
#include "pdm_priv.h"
#include "pdm_plane.h"

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
 * Static global variables
 *============================================================================*/


/*============================================================================
 * Private function definitions
 *============================================================================*/


/**
 * \brief Computes barycenter
 *
 * \param [in]   num_pts  Number of polygon vertices
 * \param [in]   pts      Polygon vertices coordinates
 * \param [out]  bary     Barycenter
 *
 */

static
void
_compute_barycenter
(
 const int     num_pts,
 const double *pts,
 double        bary[3]
)
{

  bary[0] = 0.;
  bary[1] = 0.;
  bary[2] = 0.;

  for (int i = 0; i < 3; i++) {
    for (int ipt = 0; ipt < num_pts; ipt++) {
      bary[i] += pts[3*ipt+i];
    }
    bary[i] /= num_pts;
  }

}


/*=============================================================================
 * Public function definitions
 *============================================================================*/


/**
 * \brief Computes barycenter
 *
 * \param [in]   num_pts  Number of polygon vertices
 * \param [in]   pts      Polygon vertices coordinates
 * \param [out]  bary     Barycenter
 *
 */

void
PDM_plane_barycenter
(
 const int     num_pts,
 const double *pts,
 double        n[3]
)
{
  _compute_barycenter (num_pts, pts, n);
}

/**
 * \brief Computes normal
 *
 * \param [in]   num_pts  Number of polygon vertices
 * \param [in]   pts      Polygon vertices coordinates
 * \param [out]  n        Normal
 *
 */

void
PDM_plane_normal
(
 const int     num_pts,
 const double *pts,
       double  n[3]
)
{
  double length = 0.;
  double bary[3]= {0., 0., 0.};

  n[0] = 0.;
  n[1] = 0.;
  n[2] = 0.;

  _compute_barycenter (num_pts, pts, bary);

  for (int ipt = 0; ipt < num_pts; ipt++) {

    const double *pt1 = pts + 3 * ipt;
    const double *pt2 = pts + 3 * ((ipt+1)%num_pts);
    double vect1[3];
    double vect2[3];

    for (int i = 0; i < 3; i++) {
      vect1[i] = pt1[i] - bary[i];
      vect2[i] = pt2[i] - bary[i];
    }

    n[0] += vect1[1] * vect2[2] - vect1[2] * vect2[1];
    n[1] += vect1[2] * vect2[0] - vect1[0] * vect2[2];
    n[2] += vect1[0] * vect2[1] - vect1[1] * vect2[0];

  } //over all points

  length = sqrt (n[0] * n[0] + n[1] * n[1] + n[2] * n[2]);
  if (length > 0.0) {
    n[0] /= length;
    n[1] /= length;
    n[2] /= length;
  }
}


/**
 * \brief Performs plane projection
 *
 * \param [in]   x       Point to project
 * \param [in]   origin  Plane origin
 * \param [in]   n       Plane normal
 * \param [out]  cp      Projected point
 *
 */

void
PDM_plane_projection
(
const double x[3],
const double origin[3],
const double n[3],
      double cp[3]
)
{
  double xo[3];

  xo[0] = x[0] - origin[0];
  xo[1] = x[1] - origin[1];
  xo[2] = x[2] - origin[2];

  double t = PDM_DOT_PRODUCT(n, xo);

  cp[0] = x[0] - t * n[0];
  cp[1] = x[1] - t * n[1];
  cp[2] = x[2] - t * n[2];
}


/**
 * \brief Performs plane projection
 *
 * \param [in]   x       Point to project
 * \param [in]   pt      Point inside the plane
 * \param [in]   n       Plane normal
 * \param [out]  cp      Projected point
 *
 */

void
PDM_plane_projection2
(
const double x[3],
const double pt[3],
const double n[3],
      double cp[3]
)
{
  double cst   = PDM_DOT_PRODUCT(n, pt);
  double cst1  = PDM_DOT_PRODUCT(n, x);
  double norm2 = PDM_DOT_PRODUCT(n, n);

  double t = - (cst1 - cst)/ norm2;

  cp[0] = x[0] + t * n[0];
  cp[1] = x[1] + t * n[1];
  cp[2] = x[2] + t * n[2];
}

/**
 * \brief Intersection between a plane and a line (taken from _intersect_faces_rays in pdm_inside_cloud_surf)
 *
 * \param [in]   line   Points of the line
 * \param [in]   plane  Points of the plane
 * \param [out]  ip     Intersection point
 *
 */

void
PDM_plane_line_intersection
(
const double line[6],
const double plane[9],
      double ip[3]
)
{
  // Plane normal
  double n[3];
  PDM_plane_normal(3, plane, n);

  // Line vector
  double line_vect[3] = {
    line[3] - line[0],
    line[4] - line[1],
    line[5] - line[2],
  };

  // Plane-Line vector
  double plane_line_vect[3]= {
    plane[0] - line[0],
    plane[1] - line[1],
    plane[2] - line[2],
  };

  // Intersection

  double denom = PDM_DOT_PRODUCT(line_vect, n);
  double numer = PDM_DOT_PRODUCT(plane_line_vect, n);

  const double epsilon = 1e-12;

  if (PDM_ABS(denom) < epsilon) { // epislon to avoid zero division
    if (PDM_ABS(numer) < epsilon) { // the line is in the plane
      memcpy(ip, line, sizeof(double) * 3);
    } else {
      // the line is parrallel to the plane
    }
  } else {  // the line intersects the plane
    double t = numer/denom;
    if ((t >= 0) && (t < 1)) {
      for (int j = 0; j < 3; j++) {
        ip[j] = line[j] + t*line_vect[j];
      }
    }
  }

}


#ifdef __cplusplus
}
#endif /* __cplusplus */
