/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <assert.h>
#include <math.h>
#include <float.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "pdm.h"
#include "pdm_priv.h"
#include "pdm_printf.h"
#include "pdm_logging.h"
#include "pdm_error.h"
#include "pdm_plane.h"
#include "pdm_line.h"
#include "pdm_triangle.h"
#include "pdm_polygon.h"
#include "pdm_triangulate.h"
#include "pdm_mesh_nodal.h"
#include "pdm_mean_values.h"
#include "pdm_geom_elem.h"
#include "pdm_binary_search.h"
#include "pdm_ho_location.h"
#include "pdm_ho_basis.h"
#include "pdm_ho_ordering.h"
#include "pdm_part_mesh_nodal_elmts.h"
#include "pdm_array.h"

#include "pdm_point_location.h"

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
 * Type definitions
 *============================================================================*/

/*============================================================================
 * Static global variables
 *============================================================================*/

static double _epsilon_denom = 1.e-30;       /* Minimum denominator */

/*=============================================================================
 * Private function definition
 *============================================================================*/

/**
 * \brief Compute 2x2 determinant
 *
 * \param [in]  a
 * \param [in]  b
 * \param [in]  c
 * \param [in]  d
 *
 * \return      \ref PDM_TRUE or \ref PDM_FALSE
 *
 */

static double
_det_2x2
(
 double a,
 double b,
 double c,
 double d
)
{
  return (a * d - b * c);
}


/**
 * \brief Solve 2x2 system
 *
 * \param [in]     A  Matrix
 * \param [inout]  x  right hand side in, solution out
 *
 * \return   \ref PDM_FALSE if matrix is singular, \ref PDM_TRUE otherwise
 *
 */

static int
_solve_2x2
(
 double A[2][2],
 double b[2],
 double x[2]
)
{
  const double _eps = 1e-15;

  double det = _det_2x2 (A[0][0], A[0][1], A[1][0], A[1][1]);

  if (PDM_ABS(det) < _eps) {
    return PDM_FALSE;
  }

  x[0] = (A[1][1]*b[0] - A[0][1]*b[1]) / det;
  x[1] = (-A[1][0]*b[0] + A[0][0]*b[1]) / det;

  return PDM_TRUE;
}

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

/*---------------------------------------------------------------------------
 * Solve the equation "A.x = b" with Cramer's rule.
 *
 * parameters:
 *   A[3][3] <-- equation matrix
 *   b[3]    <-- b equation right hand side
 *   x[3]    <-> equation solution (unchanged if matrix is singular)
 *
 * returns:
 *   PDM_FALSE if matrix is singular, PDM_TRUE otherwise
 *----------------------------------------------------------------------------*/

static int
_solve_3x3(double  A[3][3],
           double  b[3],
           double  x[3])
{
  double det, det_inv, x0, x1, x2;

  det = A[0][0]*(A[1][1]*A[2][2] - A[2][1]*A[1][2])
    -   A[1][0]*(A[0][1]*A[2][2] - A[2][1]*A[0][2])
    +   A[2][0]*(A[0][1]*A[1][2] - A[1][1]*A[0][2]);

  if (PDM_ABS(det) < _epsilon_denom) {
    printf("_solve_3x3: det = %e\n", det);
    return PDM_FALSE;
  }

  else {
    det_inv = 1./det;
  }

  /* Use local variables to ensure no aliasing */

  x0 = (  b[0]*(A[1][1]*A[2][2] - A[2][1]*A[1][2])
          - b[1]*(A[0][1]*A[2][2] - A[2][1]*A[0][2])
          + b[2]*(A[0][1]*A[1][2] - A[1][1]*A[0][2])) * det_inv;

  x1 = (  A[0][0]*(b[1]*A[2][2] - b[2]*A[1][2])
          - A[1][0]*(b[0]*A[2][2] - b[2]*A[0][2])
          + A[2][0]*(b[0]*A[1][2] - b[1]*A[0][2])) * det_inv;

  x2 = (  A[0][0]*(A[1][1]*b[2] - A[2][1]*b[1])
          - A[1][0]*(A[0][1]*b[2] - A[2][1]*b[0])
          + A[2][0]*(A[0][1]*b[1] - A[1][1]*b[0])) * det_inv;

  /* Copy local variables to output */
  x[0] = x0;
  x[1] = x1;
  x[2] = x2;

  return PDM_TRUE;
}


/*----------------------------------------------------------------------------
 * Compute 3d shape functions and their derivatives given element
 * parametric coordinates.
 *
 * This function is adapted from the CGNS interpolation tool.
 *
 * parameters:
 *   elt_type    <-- type of element
 *   uvw[]       <-- parametric coordinates
 *   shapef[]    --> barycenter's coordinates
 *   deriv [][]  --> derivative of shape function
 *----------------------------------------------------------------------------*/

static void
_compute_shapef_3d
(
 const PDM_Mesh_nodal_elt_t elt_type,
 const double               uvw[3],
 double                     shapef[8],
 double                     deriv[8][3]
 )
{
  double u = uvw[0];
  double v = uvw[1];
  double w = uvw[2];
  double u1 = 1. - u;
  double v1 = 1. - v;
  double w1 = 1. - w;

  switch (elt_type) {

  case PDM_MESH_NODAL_QUAD4: {
    shapef[0] = u1 * v1;
    shapef[1] = u  * v1;
    shapef[2] = u  * v ;
    shapef[3] = u1 * v ;

    if (deriv != NULL) {
      deriv[0][0] = -v1;
      deriv[0][1] = -u1;
      deriv[1][0] =  v1;
      deriv[1][1] = -u ;
      deriv[2][0] =  v ;
      deriv[2][1] =  u ;
      deriv[3][0] = -v ;
      deriv[3][1] =  u1;
    }

    break;
  }

  case PDM_MESH_NODAL_PYRAMID5: {

    if (fabs(w1) > 1.e-6) {
      w1 = 1. / w1;
    }

    shapef[0] = (1. - u - w) * (1. - v - w) * w1;
    shapef[1] =            u * (1. - v - w) * w1;
    shapef[2] =            u *            v * w1;
    shapef[3] = (1. - u - w) *            v * w1;
    shapef[4] = w;

    if (deriv != NULL) {
      deriv[0][0] = (v + w - 1.) * w1;
      deriv[0][1] = (u + w - 1.) * w1;
      deriv[0][2] = shapef[0]*w1 + deriv[0][0] + deriv[0][1];

      deriv[1][0] = (1. - v - w) * w1;
      deriv[1][1] = -u * w1;
      deriv[1][2] = shapef[1]*w1 + deriv[1][1];

      deriv[2][0] = v * w1;
      deriv[2][1] = u * w1;
      deriv[2][2] = shapef[2]*w1;

      deriv[3][0] = -v * w1;
      deriv[3][1] = (1. - u - w) * w1;
      deriv[3][2] = shapef[3]*w1 + deriv[3][0];

      deriv[4][0] = 0.;
      deriv[4][1] = 0.;
      deriv[4][2] = 1.;
    }

    break;
  }



  case PDM_MESH_NODAL_PRISM6: {

    shapef[0] = (1. - u - v) * w1;
    shapef[1] = u * w1;
    shapef[2] = v * w1;
    shapef[3] = (1. - u - v) * w;
    shapef[4] = u * w;
    shapef[5] = v * w;

    if (deriv != NULL) {
      deriv[0][0] = -w1;
      deriv[0][1] = -w1;
      deriv[0][2] = -(1. - u - v);
      deriv[1][0] =  w1;
      deriv[1][1] =  0.;
      deriv[1][2] = -u;
      deriv[2][0] =  0.;
      deriv[2][1] =  w1;
      deriv[2][2] = -v;
      deriv[3][0] = -w;
      deriv[3][1] = -w;
      deriv[3][2] =  (1. - u - v);
      deriv[4][0] =  w;
      deriv[4][1] =  0.;
      deriv[4][2] =  u;
      deriv[5][0] =  0.;
      deriv[5][1] =  w;
      deriv[5][2] =  v;
    }

    break;
  }



  case PDM_MESH_NODAL_HEXA8: {

    shapef[0] = u1 * v1 * w1;
    shapef[1] = u  * v1 * w1;
    shapef[2] = u  * v  * w1;
    shapef[3] = u1 * v  * w1;
    shapef[4] = u1 * v1 * w;
    shapef[5] = u  * v1 * w;
    shapef[6] = u  * v  * w;
    shapef[7] = u1 * v  * w;

    if (deriv != NULL) {
      deriv[0][0] = -v1 * w1;
      deriv[0][1] = -u1 * w1;
      deriv[0][2] = -u1 * v1;
      deriv[1][0] =  v1 * w1;
      deriv[1][1] = -u  * w1;
      deriv[1][2] = -u  * v1;
      deriv[2][0] =  v  * w1;
      deriv[2][1] =  u  * w1;
      deriv[2][2] = -u  * v;
      deriv[3][0] = -v  * w1;
      deriv[3][1] =  u1 * w1;
      deriv[3][2] = -u1 * v;
      deriv[4][0] = -v1 * w;
      deriv[4][1] = -u1 * w;
      deriv[4][2] =  u1 * v1;
      deriv[5][0] =  v1 * w;
      deriv[5][1] = -u  * w;
      deriv[5][2] =  u  * v1;
      deriv[6][0] =  v  * w;
      deriv[6][1] =  u  * w;
      deriv[6][2] =  u  * v;
      deriv[7][0] = -v  * w;
      deriv[7][1] =  u1 * w;
      deriv[7][2] =  u1 * v;
    }

    break;
  }



  default:
    PDM_error (__FILE__, __LINE__, 0, "Wrong element type\n");

  }

}


/*----------------------------------------------------------------------------
 * Compute hexahedron, pyramid, or prism parametric coordinates for a
 * given point.
 *
 * This function is adapted from the CGNS interpolation tool.
 *
 * parameters:
 *   elt_type            <-- type of element
 *   point_coords        <-- point coordinates
 *   vertex_coords[]     <-- pointer to element vertex coordinates
 *   tolerance           <-- location tolerance factor
 *   uvw[]               --> parametric coordinates of point in element
 *----------------------------------------------------------------------------*/

static PDM_bool_t
_compute_uvw
(
 const PDM_Mesh_nodal_elt_t elt_type,
 const double               point_coords[3],
 const double               vertex_coords[],
 const double               tolerance,
       double               uvw[3],
       double               init_uvw[3]
 )
{
  int dbg = 0;

  int i, j, n_elt_vertices, iter;
  const int max_iter = 30;
  const double tolerance2 = tolerance * tolerance;
  double dist;
  double a[3][3], b[3], x[3], shapef[8], dw[8][3];

  /* Get number of vertices */
  const int order = 1;
  n_elt_vertices = PDM_Mesh_nodal_n_vertices_element (elt_type, order);

  assert (elt_type == PDM_MESH_NODAL_QUAD4    ||
          elt_type == PDM_MESH_NODAL_PYRAMID5 ||
          elt_type == PDM_MESH_NODAL_PRISM6   ||
          elt_type == PDM_MESH_NODAL_HEXA8);

  int elt_dim = PDM_Mesh_nodal_elt_dim_get(elt_type);

  /* Use Newton-method to determine parametric coordinates and shape function */
  if (init_uvw == NULL) {
    for (i = 0; i < elt_dim; i++) {
      uvw[i] = 0.5;
    }
  }
  else {
    for (i = 0; i < elt_dim; i++) {
      uvw[i] = init_uvw[i];
    }
  }

  if (dbg) {
    log_trace(">> _compute_uvw, type %d, tolerance2 = %e\n", (int) elt_type, tolerance2);
  }
  for (iter = 0; iter < max_iter; iter++) {

    // if (dbg) {
    //   log_trace("  iter %d: uvw = ", iter);
    //   PDM_log_trace_array_double(uvw, elt_dim, "");
    // }
    _compute_shapef_3d (elt_type, uvw, shapef, dw);

    b[0] = - point_coords[0];
    b[1] = - point_coords[1];
    b[2] = - point_coords[2];

    for (i = 0; i < 3; i++) {
      for (j = 0; j < 3; j++) {
        a[i][j] = 0.0;
      }
    }

    for (i = 0; i < n_elt_vertices; i++) {

      b[0] += (shapef[i] * vertex_coords[3*i  ]);
      b[1] += (shapef[i] * vertex_coords[3*i+1]);
      b[2] += (shapef[i] * vertex_coords[3*i+2]);

      for (j = 0; j < elt_dim; j++) {
        a[0][j] -= (dw[i][j] * vertex_coords[3*i  ]);
        a[1][j] -= (dw[i][j] * vertex_coords[3*i+1]);
        a[2][j] -= (dw[i][j] * vertex_coords[3*i+2]);
      }

    }

    if (elt_dim == 3) {
      if (_solve_3x3(a, b, x) == PDM_FALSE) {
        printf("_compute_uvw: singular matrix\n");
        return PDM_FALSE;
      }
    }
    else {
      double mat[2][2];
      mat[0][0] = a[0][0]*a[0][0] + a[1][0]*a[1][0] + a[2][0]*a[2][0];
      mat[0][1] = a[0][0]*a[0][1] + a[1][0]*a[1][1] + a[2][0]*a[2][1];
      mat[1][0] = mat[0][1];
      mat[1][1] = a[0][1]*a[0][1] + a[1][1]*a[1][1] + a[2][1]*a[2][1];
      double rhs[2];
      rhs[0] = a[0][0]*b[0] + a[1][0]*b[1] + a[2][0]*b[2];
      rhs[1] = a[0][1]*b[0] + a[1][1]*b[1] + a[2][1]*b[2];
      if (_solve_2x2(mat, rhs, x) == PDM_FALSE) {
        printf("_compute_uvw: singular matrix\n");
        return PDM_FALSE;
      }
    }

    dist = 0.0;

    for (i = 0; i < elt_dim; i++) {
      dist   += x[i] * x[i];
      uvw[i] += x[i];
    }

    if (dbg) {
      log_trace("  iter %d: |b| = %e, dist = %e\n", iter, PDM_MODULE(b), dist);
    }

    if (dist <= tolerance2) {
      if (dbg) {
        _compute_shapef_3d (elt_type, uvw, shapef, dw);
        b[0] = - point_coords[0];
        b[1] = - point_coords[1];
        b[2] = - point_coords[2];
        for (i = 0; i < n_elt_vertices; i++) {
          b[0] += (shapef[i] * vertex_coords[3*i  ]);
          b[1] += (shapef[i] * vertex_coords[3*i+1]);
          b[2] += (shapef[i] * vertex_coords[3*i+2]);
        }
        log_trace("  converged : |b| = %e\n", PDM_MODULE(b));
      }
      return PDM_TRUE;
    }

  }

  return PDM_FALSE;
}


/*----------------------------------------------------------------------------
 * Locate points on an edge.
 *
 * parameters:
 *   edge_vtx           <-- ids of edge vertices
 *   parent_vertex_num  <-- pointer to parent vertex numbers (or NULL)
 *   vtx_coord          <-- pointer to vertex coordinates
 *   n_pts              <-- number of points to locate
 *   pts_coord          <-- coordinates of points to locate (size: dim * n_pts)
 *   distance           --> distance from point to edge (size: n_pts)
 *   proj_coord         --> coordinates of closest points (size: 3*n_pts)
 *   bar_coord          --> barycentric coordinates of closest points (size: 2*n_pts)
 *----------------------------------------------------------------------------*/
static void
_locate_on_edge
(
 const double       vtx_coord[],
 const int          n_pts,
 const double       pts_coord[],
 double             distance[],
 double             proj_coord[],
 double             bar_coord[]
 )
{
  int idim, ipt;

  /* Calculate edge vector and length */
  double u[3], uu = 0.;
  for (idim = 0; idim < 3; idim++) {
    u[idim] = vtx_coord[3 + idim] - vtx_coord[idim];
    uu += u[idim] * u[idim];
  }

  if (uu < _epsilon_denom){
    PDM_printf("warning _locate_on_edge : Reduce _epsilon_denom criteria : %12.5e < %12.5e\n", uu, _epsilon_denom);
    PDM_printf_flush();
    return;
  }

  double inv_uu = 1. / uu;

  /* Loop on points to locate on edge */
  double v[3], uv = 0., t;
  for (ipt = 0; ipt < n_pts; ipt++) {

    const double *_pt = pts_coord + 3 * ipt;
    double       *_bc = bar_coord + 2 * ipt;
    double       *_cp = proj_coord + 3 * ipt;

    /* Calculate linear coordinates of projection of point on edge axis */
    for (idim = 0; idim < 3; idim++) {
      v[idim] = _pt[idim] - vtx_coord[idim];
      uv += u[idim] * v[idim];
    }

    t = uv * inv_uu;


    /* Set v to be the vector from the point to the closest point on
       the segment (if t < 0, v is already that vector) */
    if (t >= 1.) {
      for (idim = 0; idim < 3; idim++) {
        v[idim] = _pt[idim] - vtx_coord[3 + idim];
        _cp[idim] = vtx_coord[3 + idim];
      }
    }

    else if (t > 0.) {
      for (idim = 0; idim < 3; idim++) {
        v[idim] -= t * u[idim];
        _cp[idim] = vtx_coord[idim] + v[idim];
      }
    }

    else {
      for (idim = 0; idim < 3; idim++) {
        _cp[idim] = vtx_coord[idim];
      }
    }

    /* Distance between point to locate and its projection */
    distance[ipt] = 0.;
    for (idim = 0; idim < 3; idim++) {
      distance[ipt] += v[idim] * v[idim];
    }

    /* Barycentric coordinates */
    _bc[0] = 1. - t;
    _bc[1] =      t;

  } // End of loop on points

}

/*----------------------------------------------------------------------------
 * Locate points in a given set of triangles.
 *
 * This function is called for sets of triangles belonging to the subdivision
 * of a given 3d face. Barycentric coordinates are used to locate the
 * projection of points.
 *
 * parameters:
 *   n_tri               <-- number of triangles
 *   tri_vtx             <-- triangles connectivity (size: n_tri * 3)
 *   parent_vertex_num   <-- pointer to parent vertex numbers (or NULL)
 *   vtx_coord           <-- pointer to vertex coordinates (size: 9)
 *   n_pts               <-- number of points to locate
 *   pts_coord           <-- point coordinates (size: 3 * n_pts)
 *   location            <-> lnum of element containing or closest to each point (size: n_pts)
 *   distance            <-> distance from point to element (size: n_pts)
 *   proj_coord          --> coordinates of closest points (size: 3*n_pts)
 *   bar_coord           <-> barcyentric coordinates of closest points (size: n_pts * 3)
 *----------------------------------------------------------------------------*/

static void
_locate_on_triangles
(
 const int          n_tri,
 const PDM_l_num_t  tri_vtx[],
 const double       vtx_coord[],
 const int          n_pts,
 const double       pts_coord[],
 int                location[],
 double             distance[],
 double             proj_coord[],
 double             bar_coord[]
 )
{
  int itri, ipt, ivtx, idim;

  /* Initialize distance of points to locate */
  for (ipt = 0; ipt < n_pts; ipt++) {
    distance[ipt] = HUGE_VAL;
  }

  if (location != NULL) {
    for (ipt = 0; ipt < n_pts; ipt++) {
      location[ipt] = -1;
    }
  }

  /* const int _order = 1;
     const int n_vtx_tri = (_order+1)*(_order+2)/2;*/
  const int n_vtx_tri = 3;

  double tri_coord[9];

  PDM_l_num_t id[3];

  double weights[3];

  /* Loop on triangles */
  for (itri = 0; itri < n_tri; itri++) {

    /* vertex index of current triangle */
    for (ivtx = 0; ivtx < 3; ivtx++) {
      id[ivtx] = tri_vtx[itri*n_vtx_tri + ivtx] - 1;
    }

    /* Calculate triangle-constant values for barycentric coordinates */
    for (ivtx = 0; ivtx < 3; ivtx++) {
      for (idim = 0; idim < 3; idim++) {
        tri_coord[3*ivtx + idim] = vtx_coord[3*id[ivtx] + idim];
      }
    }

    /* Loop on points to locate */
    for (ipt = 0; ipt < n_pts; ipt++) {

      const double *_pt = pts_coord + 3 * ipt;
      double       *_bc = bar_coord + 3 * ipt;
      double       *_cp = proj_coord + 3 * ipt;

      double closest_point[3];
      double dist2;

      PDM_triangle_status_t stat = PDM_triangle_evaluate_position(_pt,
                                                                  tri_coord,
                                                                  closest_point,
                                                                  &dist2,
                                                                  weights);

      if (stat == PDM_TRIANGLE_DEGENERATED) {
        continue;
      }

      if (dist2 < distance[ipt]) {
        if (bar_coord != NULL) {
          // PERMUTATION
          // _bc[0] = weights[1];
          // _bc[1] = weights[2];
          // _bc[2] = weights[0];
          memcpy(_bc, weights, sizeof(double)*3);
        }

        if (proj_coord != NULL) {
          _cp[0] = closest_point[0];
          _cp[1] = closest_point[1];
          _cp[2] = closest_point[2];
        }

        distance[ipt] = dist2;

        if (location != NULL) {
          location[ipt] = itri;
        }
      }

    } // End of loop on points

  } // End of loop on triangles
}


/*----------------------------------------------------------------------------
 * Locate points in a given quadrangle.
 *
 * Barycentric coordinates are used to locate the projection of points.
 *
 * parameters:
 *   quad_vtx            <-- quadrangle connectivity (size: 4)
 *   parent_vertex_num   <-- pointer to parent vertex numbers (or NULL)
 *   vtx_coord           <-- pointer to vertex coordinates (size: dim * 4)
 *   n_pts               <-- number of points to locate
 *   pts_coord           <-- point coordinates (size: dim * n_pts)
 *   distance            <-> distance from point to element (size: n_pts)
 *   proj_coord          --> coordinates of closest points (size: 3*n_pts)
 *   bar_coord           <-> barcyentric coordinates of closest points (size: n_pts * 4)
 *----------------------------------------------------------------------------*/



static void
_locate_on_quadrangle
(
 const double quad_coord[],
 const int    n_pts,
 const double pts_coord[],
 const double tolerance,
       double distance[],
       double proj_coord[],
       double bar_coord[],
       double uvw[]
 )
{
  PDM_mean_values_polygon_3d(4,
                             quad_coord,
                             n_pts,
                             pts_coord,
                             bar_coord);

  for (int ipt = 0; ipt < n_pts; ipt++) {

    const double *_pt = pts_coord + 3 * ipt;
    double *_bc  = bar_coord  + 4 * ipt;
    double *_cp  = proj_coord + 3 * ipt;
    double *_uvw = uvw        + 3 * ipt;

    /* Compute parametric coordinates with Newton method */
    double init_uvw[3 ] = {_bc[1] + _bc[2], _bc[2] + _bc[3], -1.};
    PDM_bool_t stat_uvw = _compute_uvw(PDM_MESH_NODAL_QUAD4,
                                       _pt,
                                       quad_coord,
                                       tolerance,
                                       _uvw,
                                       init_uvw);

    if (stat_uvw == PDM_TRUE) {
      if (_uvw[0] >= 0 && _uvw[0] <= 1 && _uvw[1] >= 0 && _uvw[1] <= 1) {
        _compute_shapef_3d(PDM_MESH_NODAL_QUAD4, _uvw, _bc, NULL);
      }
    }
    else {
      /* Failed to compute parametric coordinates */
      _uvw[0] = init_uvw[0];
      _uvw[1] = init_uvw[1];
    }
    _uvw[2] = -1;

    for (int idim = 0; idim < 3; idim++) {
      _cp[idim] = 0.;
    }

    for (int ivtx = 0; ivtx < 4; ivtx++) {
      for (int idim = 0; idim < 3; idim++) {
        _cp[idim] += _bc[ivtx] * quad_coord[3*ivtx + idim];
      }
    }

    double v_cp_p[3] = {_pt[0] - _cp[0], _pt[1] - _cp[1], _pt[2] - _cp[2]};
    distance[ipt] = PDM_DOT_PRODUCT (v_cp_p, v_cp_p);
  }

}


/*----------------------------------------------------------------------------
 * Locate points in a tetrahedron.
 *
 * parameters:
 *   tetra_coords        <-- tetrahedra vertex coordinates
 *   n_pts               <-- number of points to locate
 *   point_coords        <-- point coordinates
 *   distance            <-> distance from point to element (size: n_pts)
 *   proj_coord          --> coordinates of closest points (size: 3*n_pts)
 *   bar_coord           <-> barcyentric coordinates of closest points (size: n_pts * 4)
 *----------------------------------------------------------------------------*/
static void
_locate_in_tetrahedron
(
 const double        tetra_coord[12],
 const int           n_pts,
 const double        pts_coord[],
 double             *distance,
 double             *proj_coord,
 double             *bar_coord
 )
{
  int ivtx, idim, ipt, i, j, k;

  int n_pts_out = 0;
  int *pts_out = malloc (sizeof(int) * n_pts);

  double v[3][3];
  for (ivtx = 0; ivtx < 3; ivtx++) {
    for (idim = 0; idim < 3; idim++) {
      v[ivtx][idim] = tetra_coord[3*(ivtx+1) + idim] - tetra_coord[idim];
    }
  }

  double vol6 = v[0][0] * (v[1][1]*v[2][2] - v[1][2]*v[2][1]) +
    v[0][1] * (v[1][2]*v[2][0] - v[1][0]*v[2][2]) +
    v[0][2] * (v[1][0]*v[2][1] - v[1][1]*v[2][0]);


  if (fabs(vol6) < _epsilon_denom){
    PDM_printf("warning _locate_in_tetrahedron : Reduce _epsilon_denom criteria : %12.5e < %12.5e\n", vol6, _epsilon_denom);
    PDM_printf_flush();
    return;
  }

  int orientation = vol6 > 0.;

  double r[3][3];
  for (i = 0; i < 3; i++) {
    j = (i + 1) % 3;
    k = (i + 2) % 3;

    PDM_CROSS_PRODUCT (r[i], v[k], v[j]);

    for (idim = 0; idim < 3; idim++) {
      r[i][idim] /= vol6;
    }
  }

  /*
   *  Locate points inside and identify points outside
   */
  double vp0[3];
  for (ipt = 0; ipt < n_pts; ipt++) {
    const double *_pt = pts_coord + 3 * ipt;
    double       *_bc = bar_coord + 4 * ipt;
    double       *_cp = proj_coord + 3 * ipt;

    for (idim = 0; idim < 3; idim++) {
      vp0[idim] = tetra_coord[idim] - _pt[idim];
    }

    /* Compute barycentric coordinates of current point in tetrahedron */
    _bc[0] = 1.;
    for (i = 0; i < 3; i++) {
      _bc[i+1] = PDM_DOT_PRODUCT (vp0, r[i]);
      _bc[0] -= _bc[i+1];
    }

    /* Compute 'distance' from current point to tetrahedron */
    double min_bc = HUGE_VAL;
    for (ivtx = 0; ivtx < 4; ivtx++) {
      min_bc = PDM_MIN (min_bc, _bc[ivtx]);
    }

    distance[ipt] = min_bc * min_bc;
    if (min_bc > 0.) {
      distance[ipt] = -distance[ipt];
    }

    /* Point outside tetrahedron */
    if (distance[ipt] > 0.) {
      pts_out[n_pts_out++] = ipt;
    }

    /* Point inside tetrahedron */
    else {
      for (i = 0; i < 3; i++) {
        _cp[i] = _pt[i];
      }
    }

  } // End loop on points


  if (n_pts_out == 0) {
    free (pts_out);
    return;
  }


  /*
   *  Locate points outside (closest points on boundary)
   */
  double *pts_out_coord = malloc (sizeof(double) * n_pts_out * 3);
  for (ipt = 0; ipt < n_pts_out; ipt++) {
    int id_pt = pts_out[ipt];
    for (idim = 0; idim < 3; idim++) {
      pts_out_coord[3*ipt + idim] = pts_coord[3*id_pt + idim];
    }
  }

  /* Get tetrahedron's faces */
  const int cell_vtx[4] = {1, 2, 3, 4};
  int face_vtx_idx[5], face_vtx[12];
  PDM_geom_elem_tetra_faces (1,
                             orientation,
                             cell_vtx,
                             face_vtx_idx,
                             face_vtx);

  int *id_face = malloc (sizeof(int) * n_pts_out);
  double *bar_coord_face = malloc (sizeof(double) * n_pts_out * 3);
  double *closest_point_face = malloc (sizeof(double) * n_pts_out * 3);
  double *distance_face = malloc (sizeof(double) * n_pts_out);
  _locate_on_triangles (4,
                        face_vtx,
                        tetra_coord,
                        n_pts_out,
                        pts_out_coord,
                        id_face,
                        distance_face,
                        closest_point_face,
                        bar_coord_face);

  for (ipt = 0; ipt < n_pts_out; ipt++) {
    int id_pt = pts_out[ipt];
    double *_bc = bar_coord + 4 * id_pt;
    double *_cp = proj_coord + 3 * id_pt;

    for (ivtx = 0; ivtx < 4; ivtx++) {
      _bc[ivtx] = 0.;
    }

    for (ivtx = 0; ivtx < 3; ivtx++) {
      int id_vtx = face_vtx[face_vtx_idx[id_face[ipt]] + ivtx] - 1;
      _bc[id_vtx] = bar_coord_face[3*ipt + ivtx];
    }

    for (i = 0; i < 3; i++) {
      _cp[i] = closest_point_face[3*ipt + i];
    }

    distance[id_pt] = distance_face[ipt];
  }

  free (pts_out);
  free (pts_out_coord);
  free (id_face);
  free (bar_coord_face);
  free (distance_face);
  free (closest_point_face);
}


/*----------------------------------------------------------------------------
 * Locate points in a standard cell (i.e. tetrahedron, pyramid, prism or hexahedron).
 *
 * parameters:
 *   elt_type            <-- type of cell (tetrahedron, pyramid, prism or hexahedron)
 *   cell_vtx            <-- cell-vertex connectivity
 *   parent_vertex_num   <-- pointer to parent vertex numbers (or NULL)
 *   vtx_coord           <-- pointer to vertex coordinates
 *   n_pts               <-- number of points to locate
 *   pts_coord           <-- point coordinates (size: dim * n_pts)
 *   tolerance           <-- tolerance (used to check coincidence between points and vertices)
 *   distance            <-> distance from point to element (size: n_pts)
 *   proj_coord          --> coordinates of closest points (size: 3*n_pts)
 *   bar_coord           <-> barcyentric coordinates of closest points (size: n_pts * 4)
 *----------------------------------------------------------------------------*/
static void
_locate_in_cell_3d
(
 const PDM_Mesh_nodal_elt_t  elt_type,
 /*const PDM_l_num_t           cell_vtx[],
   const PDM_l_num_t          *parent_vertex_num,*/
 const double                cell_coord[],
 const int                   n_pts,
 const double                pts_coord[],
 const double                tolerance,
 double                     *distance,
 double                     *proj_coord,
 double                     *bar_coord,
 double                     *uvw
 )
{
  double _uvw[3];
  double eps_vtx2 = 1.e-6 * tolerance;
  eps_vtx2 *= eps_vtx2;

  const int order = 1;
  const int n_vtx = PDM_Mesh_nodal_n_vertices_element (elt_type, order);

  int *pts_out = malloc (sizeof(int) * n_pts);
  int n_pts_out = 0;



  /* Tetrahedron */
  if (elt_type == PDM_MESH_NODAL_TETRA4) {
    /* Shape functions may be computed directly */
    _locate_in_tetrahedron (cell_coord,
                            n_pts,
                            pts_coord,
                            distance,
                            proj_coord,
                            bar_coord);
    if (uvw != NULL) {
      for (int ipt = 0; ipt < n_pts; ipt++) {
        uvw[3*ipt  ] = bar_coord[4*ipt+1];
        uvw[3*ipt+1] = bar_coord[4*ipt+2];
        uvw[3*ipt+2] = bar_coord[4*ipt+3];
      }
    }
    return;
  }

  double *_cell_coord = NULL;
  if (elt_type == PDM_MESH_NODAL_PRISM6) {
    _cell_coord = (double *) cell_coord;
  }
  else {

    //FIXME: Allocation a chaque passage. Prevoir buffer en argument alloue a 8 * 3

    _cell_coord = malloc (sizeof(double) * n_vtx * 3);
    _cell_coord[ 0] = cell_coord[ 0];
    _cell_coord[ 1] = cell_coord[ 1];
    _cell_coord[ 2] = cell_coord[ 2];

    _cell_coord[ 3] = cell_coord[ 3];
    _cell_coord[ 4] = cell_coord[ 4];
    _cell_coord[ 5] = cell_coord[ 5];

    _cell_coord[ 6] = cell_coord[ 9];
    _cell_coord[ 7] = cell_coord[10];
    _cell_coord[ 8] = cell_coord[11];

    _cell_coord[ 9] = cell_coord[ 6];
    _cell_coord[10] = cell_coord[ 7];
    _cell_coord[11] = cell_coord[ 8];


    _cell_coord[12] = cell_coord[12];
    _cell_coord[13] = cell_coord[13];
    _cell_coord[14] = cell_coord[14];

    if (elt_type == PDM_MESH_NODAL_HEXA8) {
      _cell_coord[15] = cell_coord[15];
      _cell_coord[16] = cell_coord[16];
      _cell_coord[17] = cell_coord[17];

      _cell_coord[18] = cell_coord[21];
      _cell_coord[19] = cell_coord[22];
      _cell_coord[20] = cell_coord[23];

      _cell_coord[21] = cell_coord[18];
      _cell_coord[22] = cell_coord[19];
      _cell_coord[23] = cell_coord[20];
    }
  }


  /* Other cell types, shape functions must be computed iteratively */
  for (int ipt = 0; ipt < n_pts; ipt++) {

    const double *_pt = pts_coord + 3 * ipt;
    double *_bc = bar_coord + n_vtx * ipt;
    double *_cp = proj_coord + 3 * ipt;

    int dbg = (_pt[0] > 1.46 && _pt[0] < 1.47 &&
               _pt[1] > 14.7 && _pt[1] < 14.8 &&
               _pt[2] > 10.0 && _pt[2] < 10.1);
    dbg = 0;

    if (dbg) {
      log_trace("cell_coord = \n");
      for (int ivtx = 0; ivtx < n_vtx; ivtx++) {
        log_trace("  %f %f %f\n",
                  cell_coord[3*ivtx  ],
                  cell_coord[3*ivtx+1],
                  cell_coord[3*ivtx+2]);
      }
    }

    /* Check vertices (To avoid singularity with pyramids) */
    PDM_bool_t on_vtx = PDM_FALSE;
    for (int ivtx = 0; ivtx < n_vtx; ivtx++) {

      double dist2 = 0.;
      for (int idim = 0; idim < 3; idim++) {
        double delta = cell_coord[3*ivtx + idim] - _pt[idim];
        dist2 += delta * delta;
      }

      if (dist2 < eps_vtx2) {
        distance[ipt] = 0.;

        for (int i = 0; i < n_vtx; i++) {
          _bc[i] = 0.;
        }
        _bc[ivtx] = 1.;

        for (int i = 0; i < 3; i++) {
          _cp[i] = cell_coord[3*ivtx + i];
        }

        on_vtx = PDM_TRUE;
        break;
      }

    } // End loop on vertices

    if (on_vtx == PDM_TRUE) {
      continue;
    }

    /* Compute parametric coordinates : Try Newton */
    PDM_bool_t stat_uvw = _compute_uvw (elt_type,
                                        _pt,
                                        cell_coord,
                                        tolerance,
                                        _uvw,
                                        NULL);
    if (dbg) {
      log_trace("pt %f %f %f : stat_uvw = %d, _uvw = %f %f %f\n",
                _pt[0], _pt[1], _pt[2],
                stat_uvw,
                _uvw[0], _uvw[1], _uvw[2]);
    }

    if (stat_uvw == PDM_TRUE) {
      _compute_shapef_3d (elt_type, _uvw, _bc, NULL);

      double min_bc = HUGE_VAL;
      for (int ivtx = 0; ivtx < n_vtx; ivtx++) {
        min_bc = PDM_MIN (min_bc, _bc[ivtx]);
      }

      distance[ipt] = min_bc * min_bc;
      if (min_bc > 0.) {
        distance[ipt] = -distance[ipt];
      }
    }

    /* Failed to compute parametric coordinates */
    else {
      // Retrieve ho type
      PDM_Mesh_nodal_elt_t t_elt_ho;
      switch (elt_type) {
        case PDM_MESH_NODAL_BAR2:
          t_elt_ho = PDM_MESH_NODAL_BARHO;
          break;

        case PDM_MESH_NODAL_TRIA3:
          t_elt_ho = PDM_MESH_NODAL_TRIAHO;
          break;

        case PDM_MESH_NODAL_QUAD4:
          t_elt_ho = PDM_MESH_NODAL_QUADHO;
          break;

        case PDM_MESH_NODAL_TETRA4:
          t_elt_ho = PDM_MESH_NODAL_TETRAHO;
          break;

        case PDM_MESH_NODAL_PYRAMID5:
          t_elt_ho = PDM_MESH_NODAL_PYRAMIDHO;
          break;

        case PDM_MESH_NODAL_PRISM6:
          t_elt_ho = PDM_MESH_NODAL_PRISMHO;
          break;

        case PDM_MESH_NODAL_HEXA8:
          t_elt_ho = PDM_MESH_NODAL_HEXAHO;
          break;

        default:
          t_elt_ho = elt_type;
      }
      /*
       * Use hierarchical subdivision
       */
      double _proj_coord[3];
      double dist2 = PDM_ho_location (t_elt_ho,
                                      order,
                                      n_vtx,
                                      _cell_coord,
                                      _pt,
                                      _proj_coord,
                                      _uvw);
      /* Point inside */
      if (dist2 < 1.e-12) {
        _compute_shapef_3d (elt_type, _uvw, _bc, NULL);

        double min_bc = HUGE_VAL;
        for (int ivtx = 0; ivtx < n_vtx; ivtx++) {
          min_bc = PDM_MIN (min_bc, _bc[ivtx]);
        }

        distance[ipt] = -min_bc * min_bc;
      }

      /* Point outside */
      else {
        distance[ipt] = dist2;
      }
    }

    /* Point outside cell */
    if (distance[ipt] > 0.) {
      pts_out[n_pts_out++] = ipt;
      distance[ipt] = HUGE_VAL;
    }

     /* Point inside cell */
    else {
      for (int i = 0; i < 3; i++) {
        _cp[i] = _pt[i];
      }
    }

    if (uvw != NULL) {
      memcpy(uvw + 3*ipt, _uvw, sizeof(double) * 3);
    }

  } // End loop on points

  if (elt_type != PDM_MESH_NODAL_PRISM6) {
    free (_cell_coord);
  }


  /* Locate points outside cell (closest point on boundary) */
  if (n_pts_out > 0) {

    for (int ipt = 0; ipt < n_pts_out; ipt++) {
      int _ipt = pts_out[ipt];
      distance[_ipt] = HUGE_VAL;
    }

    // int    *closest_face  = malloc (sizeof(int)    * n_pts_out);
    // double *closest_point = malloc (sizeof(double) * n_pts_out * 3);
    // int    *inside_polygon = PDM_array_zeros_int(n_pts_out);

    int n_face = -1;
    int face_vtx_idx[7];
    int face_vtx[24];
    int _cell_vtx[8] = {1, 2, 3, 4, 5, 6, 7, 8};
    int orientation = 0;
    int n_vtx_face;

    switch (elt_type) {
    case PDM_MESH_NODAL_PYRAMID5:
      n_face = 5;

      PDM_geom_elem_pyramid_faces (1,
                                   orientation,
                                   _cell_vtx,
                                   face_vtx_idx,
                                   face_vtx);
      break;

    case PDM_MESH_NODAL_PRISM6:
      n_face = 5;

      PDM_geom_elem_prism_faces (1,
                                 orientation,
                                 _cell_vtx,
                                 face_vtx_idx,
                                 face_vtx);
      break;

    case PDM_MESH_NODAL_HEXA8:
      n_face = 6;

      PDM_geom_elem_hexa_faces (1,
                                orientation,
                                _cell_vtx,
                                face_vtx_idx,
                                face_vtx);
      break;

    default:
      PDM_error (__FILE__, __LINE__, 0, "Wrong standard element type\n");
      break;

    }

    PDM_l_num_t n_tri;
    double tri_coord[9];
    PDM_l_num_t _tri_vtx[3], tri_vtx[6];
    // PDM_triangulate_state_t *state = PDM_triangulate_state_create (4);

    /* Find closest face/point for each point outside the cell */
    for (int iface = 0; iface < n_face; iface++) {

      const PDM_l_num_t *_face_vtx = face_vtx + face_vtx_idx[iface];
      n_vtx_face = face_vtx_idx[iface+1] - face_vtx_idx[iface];

      /* Triangular face */
      if (n_vtx_face == 3) {
        n_tri = 1;
        for (int ivtx = 0; ivtx < 3; ivtx++) {
          tri_vtx[ivtx] = _face_vtx[ivtx];
        }
      }

      /* Quadrilateral face */
      else {
        n_tri = PDM_triangulate_quadrangle (3,
                                            (double *) cell_coord,
                                            NULL,
                                            _face_vtx,
                                            tri_vtx);
      }

      /* Loop on face triangles */
      for (int itri = 0; itri < n_tri; itri++) {

        for (int ivtx = 0; ivtx < 3; ivtx++) {
          _tri_vtx[ivtx] = tri_vtx[3*itri + ivtx] - 1;
        }

        for (int ivtx = 0; ivtx < 3; ivtx++) {
          for (int idim = 0; idim < 3; idim++) {
            tri_coord[3*ivtx + idim] = cell_coord[3*_tri_vtx[ivtx] + idim];
          }
        }

        /* Loop on points */
        for (int ipt = 0; ipt < n_pts_out; ipt++) {

          int _ipt = pts_out[ipt];
          const double *_pt = pts_coord     + 3 * _ipt;
          // double       *_cp = closest_point + 3 *  ipt;
          double       *_pp = proj_coord + 3 *  _ipt;
          double       *_bc = bar_coord  + n_vtx * _ipt;
          int dbg = (_pt[0] > 1.46 && _pt[0] < 1.47 &&
                     _pt[1] > 14.7 && _pt[1] < 14.8 &&
                     _pt[2] > 10.0 && _pt[2] < 10.1);
          dbg = 0;

          if (dbg) {
            setenv("DBG_TRIANGLE", "1", 1);
          }
          else {
            setenv("DBG_TRIANGLE", "0", 1);
          }

          double min_dist2, closest[3], weight[3];
          PDM_triangle_status_t error = PDM_triangle_evaluate_position(_pt,
                                                                       tri_coord,
                                                                       closest,
                                                                       &min_dist2,
                                                                       weight);

          if (dbg) {
            double err[3] = {closest[0], closest[1], closest[2]};
            for (int i = 0; i < 3; i++) {
              for (int j = 0; j < 3; j++) {
                err[j] -= weight[i]*tri_coord[3*i+j];
              }
            }
            log_trace("  err = %e\n", PDM_MODULE(err));

            log_trace("pt %20.16f %20.16f %20.16f at dist %20.16f from tri #%d of face %d\n",
                      _pt[0], _pt[1], _pt[2], sqrt(min_dist2), itri, iface);
            for (int i = 0; i < 3; i++) {
              log_trace("  %20.16f %20.16f %20.16f\n",
                        tri_coord[3*i], tri_coord[3*i+1], tri_coord[3*i+2]);
            }
          }


          if (error == PDM_TRIANGLE_DEGENERATED) {
            continue;
          }

          // if (error == PDM_TRIANGLE_INSIDE) {
          //   inside_polygon[ipt] = 1;
          // }

          if (distance[_ipt] > min_dist2) {
            distance[_ipt] = min_dist2;

            // closest_face[ipt] = iface;
            for (int idim = 0; idim < 3; idim++) {
              _pp[idim] = closest[idim];
            }

            for (int ivtx = 0; ivtx < n_vtx; ivtx++) {
              _bc[ivtx] = 0;
            }

            for (int ivtx = 0; ivtx < 3; ivtx++) {
              _bc[_tri_vtx[ivtx]] = weight[ivtx];
              if (dbg) {
                log_trace(" += %f * (%f %f %f)\n",
                          _bc[_tri_vtx[ivtx]],
                          cell_coord[3*_tri_vtx[ivtx] + 0],
                          cell_coord[3*_tri_vtx[ivtx] + 1],
                          cell_coord[3*_tri_vtx[ivtx] + 2]);
              }
            }
            if (dbg) {
              log_trace(" pp = %f %f %f\n", _pp[0], _pp[1], _pp[2]);
            }

          }

        } // End loop on points

      } // End loop on triangles

    } // End loop on faces





    // /* Find closest face/point for each point outside the cell */
    // for (int iface = 0; iface < n_face; iface++) {

    //   const PDM_l_num_t *_face_vtx = face_vtx + face_vtx_idx[iface];
    //   n_vtx_face = face_vtx_idx[iface+1] - face_vtx_idx[iface];

    //   /* Triangular face */
    //   if (n_vtx_face == 3) {
    //     n_tri = 1;
    //     for (int ivtx = 0; ivtx < 3; ivtx++) {
    //       tri_vtx[ivtx] = _face_vtx[ivtx];
    //     }
    //   }

    //   /* Quadrilateral face */
    //   else {
    //     n_tri = PDM_triangulate_quadrangle (3,
    //                                         (double *) cell_coord,
    //                                         NULL,
    //                                         _face_vtx,
    //                                         tri_vtx);
    //   }

    //   /* Loop on face triangles */
    //   for (int itri = 0; itri < n_tri; itri++) {

    //     for (int ivtx = 0; ivtx < 3; ivtx++) {
    //       _tri_vtx[ivtx] = tri_vtx[3*itri + ivtx] - 1;
    //     }

    //     for (int ivtx = 0; ivtx < 3; ivtx++) {
    //       for (int idim = 0; idim < 3; idim++) {
    //         tri_coord[3*ivtx + idim] = cell_coord[3*_tri_vtx[ivtx] + idim];
    //       }
    //     }

    //     /* Loop on points */
    //     for (int ipt = 0; ipt < n_pts_out; ipt++) {

    //       int _ipt = pts_out[ipt];
    //       const double *_pt = pts_coord     + 3 * _ipt;
    //       double       *_cp = closest_point + 3 *  ipt;
    //       int dbg = (_pt[0] > -8.323 && _pt[0] < -8.322 &&
    //                  _pt[1] > 13.33  && _pt[1] < 13.34  &&
    //                  _pt[2] > 3.89   && _pt[2] < 3.90);

    //       if (dbg) {
    //         setenv("DBG_TRIANGLE", "1", 1);
    //       }
    //       else {
    //         setenv("DBG_TRIANGLE", "0", 1);
    //       }

    //       double min_dist2, closest[3], weight[3];
    //       PDM_triangle_status_t error = PDM_triangle_evaluate_position(_pt,
    //                                                                    tri_coord,
    //                                                                    closest,
    //                                                                    &min_dist2,
    //                                                                    weight);

    //       if (dbg) {
    //         double err[3] = {closest[0], closest[1], closest[2]};
    //         for (int i = 0; i < 3; i++) {
    //           for (int j = 0; j < 3; j++) {
    //             err[j] -= weight[i]*tri_coord[3*i+j];
    //           }
    //         }
    //         log_trace("  err = %e\n", PDM_MODULE(err));

    //         log_trace("pt %20.16f %20.16f %20.16f at dist %20.16f from tri #%d of face %d\n",
    //                   _pt[0], _pt[1], _pt[2], sqrt(min_dist2), itri, iface);
    //         for (int i = 0; i < 3; i++) {
    //           log_trace("  %20.16f %20.16f %20.16f\n",
    //                     tri_coord[3*i], tri_coord[3*i+1], tri_coord[3*i+2]);
    //         }
    //       }


    //       if (error == PDM_TRIANGLE_DEGENERATED) {
    //         continue;
    //       }

    //       if (error == PDM_TRIANGLE_INSIDE) {
    //         inside_polygon[ipt] = 1;
    //       }

    //       if (distance[_ipt] > min_dist2) {
    //         distance[_ipt] = min_dist2;

    //         closest_face[ipt] = iface;
    //         for (int idim = 0; idim < 3; idim++) {
    //           _cp[idim] = closest[idim];
    //         }
    //       }

    //     } // End loop on points

    //   } // End loop on triangles

    // } // End loop on faces

    // state = PDM_triangulate_state_destroy(state);

    // /* Loctate closest points on closest faces */
    // double face_coord[12];
    // double bar_coord_face[4];
    // for (int ipt = 0; ipt < n_pts_out; ipt++) {

    //   int _ipt = pts_out[ipt];

    //   const double *_pt = pts_coord     + 3     * _ipt;
    //   double       *_cp = closest_point + 3     *  ipt;
    //   double       *_bc = bar_coord     + n_vtx * _ipt;

    //   int dbg = (_pt[0] > -8.323 && _pt[0] < -8.322 &&
    //              _pt[1] > 13.33  && _pt[1] < 13.34  &&
    //              _pt[2] > 3.89   && _pt[2] < 3.90);

    //   for (int ivtx = 0; ivtx < n_vtx; ivtx++) {
    //     _bc[ivtx] = 0.;
    //   }

    //   int iface = closest_face[ipt];
    //   n_vtx_face = face_vtx_idx[iface+1] - face_vtx_idx[iface];

    //   for (int ivtx = 0; ivtx < n_vtx_face; ivtx++) {
    //     int _ivtx = face_vtx[face_vtx_idx[iface] + ivtx] - 1;

    //     for (int idim = 0; idim < 3; idim++) {
    //       face_coord[3*ivtx + idim] = cell_coord[3*_ivtx + idim];
    //     }
    //   }

    //   if (dbg) {
    //     log_trace(">>> cp = %20.16e %20.16e %20.16e\n", _cp[0], _cp[1], _cp[2]);
    //   }
    //   if (inside_polygon[ipt]) {
    //     PDM_mean_values_polygon_3d (n_vtx_face,
    //                                 face_coord,
    //                                 1,
    //                                 _cp,
    //                                 bar_coord_face);
    //   }
    //   else {

    //   }

    //   if (dbg) {
    //     log_trace("pt %f %f %f : closest_face %d, dist = %f,  coord :\n",
    //               _pt[0], _pt[1], _pt[2], iface, distance[_ipt]);
    //     for (int i = 0; i < n_vtx_face; i++) {
    //       log_trace("  %20.16e %20.16e %20.16e\n",
    //                 face_coord[3*i], face_coord[3*i+1], face_coord[3*i+2]);
    //     }
    //     PDM_log_trace_array_double(bar_coord_face, n_vtx_face, "  mvc in face : ");
    //   }

    //   for (int ivtx = 0; ivtx < n_vtx_face; ivtx++) {
    //     int _ivtx = face_vtx[face_vtx_idx[iface] + ivtx] - 1;
    //     _bc[_ivtx] = bar_coord_face[ivtx];
    //   }

    //   // const double *_pt = pts_coord + 3 * _ipt;
    //   double *_pp = proj_coord + 3 * _ipt;
    //   // // for (int idim = 0; idim < 3; idim++) {
    //   // //   _pp[idim] = 0.;
    //   // // }

    //   for (int ivtx = 0; ivtx < n_vtx_face; ivtx++) {
    //     int _ivtx = face_vtx[face_vtx_idx[iface] + ivtx] - 1;
    //     if (dbg) {
    //       log_trace(" += %f * (%f %f %f)\n",
    //                 bar_coord_face[ivtx],
    //                 face_coord[3*ivtx + 0],
    //                 face_coord[3*ivtx + 1],
    //                 face_coord[3*ivtx + 2]);
    //                 // cell_coord[3*_ivtx + 0],
    //                 // cell_coord[3*_ivtx + 1],
    //                 // cell_coord[3*_ivtx + 2]);
    //     }
    //     for (int idim = 0; idim < 3; idim++) {
    //       // _pp[idim] += bar_coord_face[ivtx] * cell_coord[3*_ivtx + idim];
    //       _pp[idim] += bar_coord_face[ivtx] * face_coord[3*ivtx + idim];
    //     }
    //   }

    //   double err = 0;
    //   for (int idim = 0; idim < 3; idim++) {
    //     double delta = _cp[idim] - _pp[idim];
    //     err += delta*delta;
    //   }

    //   if (dbg) {
    //     log_trace(  "err = %e\n", sqrt(err));
    //   }

    //   double v_p_cp[3] = {_pt[0] - _pp[0], _pt[1] - _pp[1], _pt[2] - _pp[2]};
    //   distance[_ipt] = PDM_DOT_PRODUCT (v_p_cp, v_p_cp);
    // }

    // free (closest_face);
    // free (closest_point);
    // free (inside_polygon);
  }

  free (pts_out);

}



/*----------------------------------------------------------------------------
 * Locate points in a given 3d polygon.
 *
 * Barycentric coordinates are used to locate the projection of points.
 *
 * parameters:
 *   n_vtx               <-- number of vertices
 *   vtx_coord           <-- pointer to vertex coordinates (size: 3 * n_vtx)
 *   n_pts               <-- number of points to locate
 *   pts_coord           <-- point coordinates (size: 3 * n_pts)
 *   distance            <-> distance from point to element (size: n_pts)
 *   proj_coord          --> coordinates of closest points (size: 3*n_pts)
 *   bar_coord           <-> barcyentric coordinates of closest points (size: n_pts * n_vtx)
 *----------------------------------------------------------------------------*/
static void
_locate_in_polygon
(
 const PDM_l_num_t n_vtx,
 const double      vtx_coord[],
 const int         n_pts,
 const double      pts_coord[],
 double            distance[],
 double            proj_coord[],
 double            bar_coord[]
 )
{
  /* Compute mean value coordinates of closest points on polygon */
  PDM_mean_values_polygon_3d (n_vtx,
                              vtx_coord,
                              n_pts,
                              pts_coord,
                              bar_coord);

  /* Compute distances */
  for (int ipt = 0; ipt < n_pts; ipt++) {
    const double *_pt = pts_coord  + ipt * 3;
    double       *_bc = bar_coord  + ipt * n_vtx;
    double       *_cp = proj_coord + ipt * 3;

    for (int idim = 0; idim < 3; idim++) {
      _cp[idim] = 0.;
    }

    for (int ivtx = 0; ivtx < n_vtx; ivtx++) {
      for (int idim = 0; idim < 3; idim++) {
        _cp[idim] += _bc[ivtx] * vtx_coord[3*ivtx + idim];
      }
    }

    double v_cp_p[3] = {_pt[0] - _cp[0], _pt[1] - _cp[1], _pt[2] - _cp[2]};
    distance[ipt] = PDM_DOT_PRODUCT (v_cp_p, v_cp_p);
    int dbg = 0;
    // int dbg = (_pt[0] > -0.424127 && _pt[0] < -0.424125 &&
    //            _pt[1] > -0.029731 && _pt[1] < -0.029729);
    if (dbg) {
      log_trace("distance = %e\n", distance[ipt]);
    }

  }
}




static void
_compute_mean_value_coord_polyhedron_point_inside
(
 const int                      n_vtx,
 const double                   vtx_coord[],
 const int                      n_face,
 const int                      face_vtx_idx[],
 const int                      face_vtx[],
 const int                      face_orientation[],
 const double                   pt_coord[],
       double                   weight[],
       PDM_triangulate_state_t *tri_state,
       int                     *tri_vtx
 )
 {
  double eps_in_triangle = 1e-12;

  /* Initialize weights */
  for (int i = 0; i < n_vtx; i++) {
    weight[i] = 0.;
  }

  /*
   * Compute distances d between p and the vertices
   * and unit vectors u from point to vertex
   */
  double d[n_vtx];
  double u[n_vtx][3];

  for (int ivtx = 0; ivtx < n_vtx; ivtx++) {
    d[ivtx] = 0.;
    for (int j = 0; j < 3; j++) {
      u[ivtx][j] = vtx_coord[3*ivtx+j] - pt_coord[j];
      d[ivtx] += u[ivtx][j]*u[ivtx][j];
    }

    d[ivtx] = sqrt(d[ivtx]);
    double invd = 1/d[ivtx];
    for (int j = 0; j < 3; j++) {
      u[ivtx][j] *= invd;
    }
  } // End of loop on vertices

  /* Loop on faces */
  int on_triangle = 0;
  for (int iface = 0; iface < n_face; iface++) {

    const int *_face_vtx = face_vtx + face_vtx_idx[iface];
    int face_vtx_n = face_vtx_idx[iface+1] - face_vtx_idx[iface];

    /* Triangulate face if necessary */
    int n_tri;
    if (face_vtx_n == 3) {
      /* Triangular face */
      n_tri = 1;
      memcpy(tri_vtx, _face_vtx, sizeof(int) * 3);
    }
    else if (face_vtx_n == 4) {
      /* Quadrilateral face */
      n_tri = PDM_triangulate_quadrangle(3,
                                         vtx_coord,
                                         NULL,
                                         _face_vtx,
                                         tri_vtx);
    }
    else {
      /* Polygonal face */
      n_tri = PDM_triangulate_polygon(3,
                                      face_vtx_n,
                                      vtx_coord,
                                      NULL,
                                      _face_vtx,
                                      PDM_TRIANGULATE_MESH_DEF,
                                      tri_vtx,
                                      tri_state);
    }


    /* Loop on triangles */
    double *tri_u[3];
    // double tri_eps2[3];
    for (int itri = 0; itri < n_tri; itri++) {

      int __tri_vtx[3];
      if (face_orientation[iface] > 0) {
        for (int ivtx = 0; ivtx < 3; ivtx++) {
          __tri_vtx[ivtx] = tri_vtx[3*itri + ivtx] - 1;
        }
      }
      else {
        for (int ivtx = 0; ivtx < 3; ivtx++) {
          __tri_vtx[ivtx] = tri_vtx[3*itri + 2-ivtx] - 1;
        }
      }

      for (int ivtx = 0; ivtx < 3; ivtx++) {
        tri_u   [ivtx] = u[__tri_vtx[ivtx]];
        // tri_eps2[ivtx] = vtx_eps2[__tri_vtx[ivtx]];
      }

      double l[3];
      double theta[3];
      double sin_theta[3];
      double h = 0.;
      for (int ivtx = 0; ivtx < 3; ivtx++) {
        int ip = (ivtx+1)%3;
        int im = (ivtx+2)%3;
        l[ivtx] = 0.;
        for (int j = 0; j < 3; j++) {
          double delta = tri_u[ip][j] - tri_u[im][j];
          l[ivtx] += delta*delta;
        }

        // assert(l[ivtx] > tri_eps2[ip] && l[ivtx] > tri_eps2[im]);

        l[ivtx] = sqrt(l[ivtx]);

        theta[ivtx] = asin(PDM_MIN(1., 0.5*l[ivtx]));
        h += theta[ivtx];
        theta[ivtx] *= 2;

        sin_theta[ivtx] = sin(theta[ivtx]);
      }

      /* Check if p is coplanar with current triangle to avoid division by zero */
      on_triangle = (PDM_PI - h < eps_in_triangle);

      if (on_triangle) {
        // log_trace("!! on triangle, pt_coord = %f %f %f\n",
        //           pt_coord[0], pt_coord[1], pt_coord[2]);
        // log_trace("tri_coord :\n");
        // for (int ivtx = 0; ivtx < 3; ivtx++) {
        //   PDM_log_trace_array_double(vtx_coord + 3*__tri_vtx[ivtx], 3, "");
        // }

        double tri_coord[9];
        for (int ivtx = 0; ivtx < 3; ivtx++) {
          memcpy(tri_coord + 3*ivtx,
                 vtx_coord + 3*__tri_vtx[ivtx],
                 sizeof(double) * 3);
        }

        double tri_closest_point[3];
        double tri_dist2;
        double tri_weight[3];
        PDM_triangle_status_t stat = PDM_triangle_evaluate_position(pt_coord,
                                                                    tri_coord,
                                                                    tri_closest_point,
                                                                    &tri_dist2,
                                                                    tri_weight);
        assert(stat == PDM_TRIANGLE_INSIDE);
        for (int ivtx = 0; ivtx < n_vtx; ivtx++) {
          weight[ivtx] = 0;
        }
        for (int ivtx = 0; ivtx < 3; ivtx++) {
          weight[__tri_vtx[ivtx]] = tri_weight[ivtx];
        }

        break;
      }

      double sdet = PDM_SIGN(_determinant_3x3(tri_u[0],
                                              tri_u[1],
                                              tri_u[2]));

      double c[3];
      double s[3];
      double sin_h = sin(h);
      int coplanar_outside = 0;
      for (int ivtx = 0; ivtx < 3; ivtx++) {
        int ip = (ivtx+1)%3;
        int im = (ivtx+2)%3;

        // check denominator ?
        c[ivtx] = -1 + (2 * sin_h * sin(h - theta[ivtx])) /
        (sin_theta[ip] * sin_theta[im]);

        s[ivtx] = sdet * sqrt(PDM_MAX(0., 1 - c[ivtx]*c[ivtx]));

        if (PDM_ABS(s[ivtx]) < eps_in_triangle) {
          coplanar_outside = 1;
          break;
        }
      }

      if (coplanar_outside) {
        /* p lies outside current triangle, on the same plane */
        continue;
      }

      for (int i = 0; i < 3; i++) {
        int ip = (i+1)%3;
        int im = (i+2)%3;
        weight[__tri_vtx[i]] += (theta[i] - c[ip]*theta[im] - c[im]*theta[ip]) /
        (d[__tri_vtx[i]] * sin_theta[ip] * s[im]);
      }

    } // End of loop on triangles

    if (on_triangle) break;

  } // End of loop on faces


  /* Normalize weights */
  double sum_w = 0.;
  for (int ivtx = 0; ivtx < n_vtx; ivtx++) {
    sum_w += weight[ivtx];
  }
  double isum_w = 1./sum_w;
  for (int ivtx = 0; ivtx < n_vtx; ivtx++) {
    weight[ivtx] *= isum_w;
  }
 }

/*----------------------------------------------------------------------------
 * Locate points in a given polyhedron.
 *
 * parameters:
 *   n_vtx               <-- number of vertices
 *   vtx_coord           <-- pointer to vertex coordinates (size: 3 * n_vtx)
 *   n_face              <-- number of faces
 *   face_vtx_idx        <-- index of face-vertex connectivity (size: n_face + 1)
 *   face_vtx            <-- face-vertex connectivity (size: face_vtx_idx[n_face])
 *   face_orientation    <-- orientation of faces in current polyhedron (size: n_face)
 *   n_pts               <-- number of points to locate
 *   pts_coord           <-- point coordinates (size: 3 * n_pts)
 *   distance            <-> distance from point to element (size: n_pts)
 *   proj_coord          --> coordinates of closest points (size: 3*n_pts)
 *   bar_coord           <-> barcyentric coordinates of closest points (size: n_pts * n_vtx)
 *----------------------------------------------------------------------------*/

static int
_locate_in_polyhedron
(
 const int                      n_vtx,
 const double                   vtx_coord[],
 const int                      n_face,
 const int                      face_vtx_idx[],
 const int                      face_vtx[],
 const int                      face_orientation[],
 const int                      n_pts,
 const double                   pts_coord[],
       double                   distance[],
       double                   proj_coord[],
       double                   weight[],
       PDM_triangulate_state_t *tri_state,
       int                     *tri_vtx
 )
{
  int dbg = 0;

  double eps_on_face2 = 1.e-24;

  /* Prepare for triangulation */
  PDM_triangulate_state_t *_tri_state = tri_state;
  int                     *_tri_vtx   = tri_vtx;
  if (_tri_state == NULL || _tri_vtx == NULL) {
    int max_face_vtx_n = 0;
    for (int iface = 0; iface < n_face; iface++) {
      int face_vtx_n = face_vtx_idx[iface+1] - face_vtx_idx[iface];
      max_face_vtx_n = PDM_MAX(max_face_vtx_n, face_vtx_n);
    }

    if (_tri_state == NULL && max_face_vtx_n > 4) {
      _tri_state = PDM_triangulate_state_create(max_face_vtx_n);
    }

    if (_tri_vtx == NULL) {
      _tri_vtx = malloc(sizeof(int) * (max_face_vtx_n - 2) * 3);
    }
  }


  /*
   * First, identify points inside/outside polyhedron
   *  - compute sum of solid angles
   *  - compute distance from closest subtriangle
   */
  double solid_angle[n_pts];
  for (int ipt = 0; ipt < n_pts; ipt++) {
    solid_angle [ipt] = 0.;
    distance    [ipt] = HUGE_VAL;
  }

  /* Loop on faces */
  for (int iface = 0; iface < n_face; iface++) {

    const int *_face_vtx = face_vtx + face_vtx_idx[iface];
    int face_vtx_n = face_vtx_idx[iface+1] - face_vtx_idx[iface];

    /* Triangulate face if necessary */
    int n_tri;
    if (face_vtx_n == 3) {
      /* Triangular face */
      n_tri = 1;
      memcpy(_tri_vtx, _face_vtx, sizeof(int) * 3);
    }
    else if (face_vtx_n == 4) {
      /* Quadrilateral face */
      n_tri = PDM_triangulate_quadrangle(3,
                                         vtx_coord,
                                         NULL,
                                         _face_vtx,
                                         _tri_vtx);
    }
    else {
      /* Polygonal face */
      n_tri = PDM_triangulate_polygon(3,
                                      face_vtx_n,
                                      vtx_coord,
                                      NULL,
                                      _face_vtx,
                                      PDM_TRIANGULATE_MESH_DEF,
                                      _tri_vtx,
                                      _tri_state);
    }

    /* Loop on subtriangles */
    for (int itri = 0; itri < n_tri; itri++) {

      double tri_coord[9];
      // if (face_orientation[iface] > 0) {
      //   for (int idx_vtx = 0; idx_vtx < 3; idx_vtx++) {
      //     int vtx_id = _tri_vtx[3*itri + idx_vtx] - 1;
      //     memcpy(tri_coord + 3*idx_vtx, vtx_coord + 3*vtx_id, sizeof(double) * 3);
      //   }
      // } else {
      //   for (int idx_vtx = 0; idx_vtx < 3; idx_vtx++) {
      //     int vtx_id = _tri_vtx[3*itri + 2-idx_vtx] - 1;
      //     memcpy(tri_coord + 3*idx_vtx, vtx_coord + 3*vtx_id, sizeof(double) * 3);
      //   }
      // }
      if (face_orientation[iface] < 0) {
        int tmp = _tri_vtx[3*itri];
        _tri_vtx[3*itri  ] = _tri_vtx[3*itri+2];
        _tri_vtx[3*itri+2] = tmp;
      }


      for (int idx_vtx = 0; idx_vtx < 3; idx_vtx++) {
        int vtx_id = _tri_vtx[3*itri + idx_vtx] - 1;
        memcpy(tri_coord + 3*idx_vtx, vtx_coord + 3*vtx_id, sizeof(double) * 3);
      }

      /* Loop on points */
      for (int ipt = 0; ipt < n_pts; ipt++) {

        if (distance[ipt] < eps_on_face2) {
          continue;
        }

        const double *pt_coord = pts_coord + 3*ipt;

        /* Distance */
        double tri_dist2, tri_closest_point[3], tri_weight[3];
        PDM_triangle_status_t stat = PDM_triangle_evaluate_position(pt_coord,
                                                                    tri_coord,
                                                                    tri_closest_point,
                                                                    &tri_dist2,
                                                                    tri_weight);

        if (stat == PDM_TRIANGLE_DEGENERATED) {
          /* Raise error? */
          continue;
        }

        if (tri_dist2 < distance[ipt]) {
          distance    [ipt] = tri_dist2;

          memcpy(proj_coord + 3*ipt, tri_closest_point, sizeof(double) * 3);

          double *w = weight + ipt*n_vtx;
          for (int i = 0; i < n_vtx; i++) {
            w[i] = 0;
          }

          for (int i = 0; i < 3; i++) {
            int vtx_id = _tri_vtx[3*itri + i] - 1;
            w[vtx_id] = tri_weight[i];
          }

          if (distance[ipt] < eps_on_face2) {
            continue;
          }
        }

        /* Solid angle */
        double v[3][3], lv[3];
        double denom = 1.;
        for (int i = 0; i < 3; i++) {
          for (int j = 0; j < 3; j++) {
            v[i][j] = tri_coord[3*i + j] - pt_coord[j];
          }
          lv[i] = PDM_MODULE(v[i]);
          denom *= lv[i];
        }
        for (int i = 0; i < 3; i++) {
          denom += lv[i] * PDM_DOT_PRODUCT(v[(i+1)%3], v[(i+2)%3]);
        }

        double numer = _determinant_3x3(v[0], v[1], v[2]);

        double half_angle = atan2(numer, denom);
        if (half_angle < 0. && numer > 0) {
          half_angle = 2.*PDM_PI - half_angle;
        }
        else if (half_angle > 0. && numer < 0){
          half_angle = -2.*PDM_PI + half_angle;
        }

        solid_angle[ipt] += 2*half_angle;

      } // End of loop on points

    } // End of loop on subtriangles

  } // End of loop on faces


  /*
   * Locate points (compute mean value coordinates)
   */
  double eps_solid_angle = 1.e-6;
  double threshold_inside = 4*PDM_PI*(1 - eps_solid_angle);

  for (int ipt = 0; ipt < n_pts; ipt++) {

    const double *p  = pts_coord  + ipt*3;
    double       *w  = weight     + ipt*n_vtx;
    double       *pc = proj_coord + ipt*3;

    if (dbg) {
      log_trace("point %d (%f %f %f), solid_angle = %f*PI, dist2 = %e\n",
                ipt, p[0], p[1], p[2], solid_angle[ipt]/PDM_PI, distance[ipt]);
    }


    if (solid_angle[ipt] > threshold_inside) {
      /* Point strictly inside polyhedron */
      distance[ipt] = -distance[ipt];

      // Compute mean value coordinates
      _compute_mean_value_coord_polyhedron_point_inside(n_vtx,
                                                        vtx_coord,
                                                        n_face,
                                                        face_vtx_idx,
                                                        face_vtx,
                                                        face_orientation,
                                                        p,
                                                        w,
                                                        _tri_state,
                                                        _tri_vtx);

      for (int j = 0; j < 3; j++) {
        pc[j] = 0.;
      }

      for (int ivtx = 0; ivtx < n_vtx; ivtx++) {
        const double *vc = vtx_coord + 3*ivtx;
        for (int j = 0; j < 3; j++) {
          pc[j] += w[ivtx] * vc[j];
        }
      }

      if (dbg) {
        log_trace("  interior\n");
        double e[3] = {pc[0] - p[0], pc[1] - p[1], pc[2] - p[2]};
        double err = PDM_DOT_PRODUCT(e, e);
        // PDM_log_trace_array_double(w, n_vtx, "  mean_value_coord : ");
        log_trace("  proj_coord = %f %f %f, dist2 = %e\n", pc[0], pc[1], pc[2], err);
      }

    }

    else if (distance[ipt] > eps_on_face2) {
      if (dbg) {
        /* Point strictly outside polyhedron */
        log_trace("  exterior\n");
        // PDM_log_trace_array_double(w, n_vtx, "  mean_value_coord : ");
        log_trace("  proj_coord = %f %f %f, dist2 = %e\n", pc[0], pc[1], pc[2], distance[ipt]);
      }
    }

    else {
      if (dbg) {
        /* Point on a polyhedron face */
        log_trace("  on face\n");
        // PDM_log_trace_array_double(w, n_vtx, "  mean_value_coord : ");
        log_trace("  proj_coord = %f %f %f, dist2 = %e\n", pc[0], pc[1], pc[2], distance[ipt]);
      }
    }


  } // End of loop on points



  if (_tri_state != tri_state && _tri_state != NULL) {
    PDM_triangulate_state_destroy(_tri_state);
  }

  if (_tri_vtx != tri_vtx) {
    free(_tri_vtx);
  }


  return 0;
}

/*============================================================================
 * Public function definitions
 *============================================================================*/

/**
 * \brief Locate a set points inside a set of elements
 *
 * Elements are ordered by type (points, lines, triangles, quadrangles,
 * polygons, tetrahedra, pyramids, prisms, hexahedra, polyhedra).
 *
 * \param [in]   type_idx           Index for the element types (size = 11)
 * \param [in]   elt_vtx_idx        Index of the element-vertex connectivity
 * \param [in]   elt_vtx_coord      Coordinates of the elements' vertices
 * \param [in]   poly3d_face_idx    Index of the element-face connectivity (only for polyhedra)
 * \param [in]   face_vtx_idx       Index for the face-vertex connectivity
 * \param [in]   face_vtx           Face-vertex connectivity
 * \param [in]   face_orientation   Orientation of the faces
 * \param [in]   pts_idx            Index of points (size = n_elt + 1)
 * \param [in]   pts_coord          Coordinates of the points to locate
 * \param [in]   tolerance          Geometric tolerance
 * \param [out]  distance           Distance from points to elements (< 0 if inside, > 0 if outside)
 * \param [out]  projected_coord    Coordinates of the projection of the points on the elements
 * \param [out]  bar_coord_idx      Index for the mean-value coordinates of the projections
 * \param [out]  bar_coord          Mean-value coordinates of the projections
 *
 */


/* Awkward to pass n_part and pvtx_coord aside from pmne... */
void
PDM_point_location_nodal
(
 PDM_part_mesh_nodal_elmts_t   *pmne,
 const int                      n_part,
 const double                 **pvtx_coord,
 const int                    **pts_idx,
 const double                 **pts_coord,
 const double                   tolerance,
 double                      ***distance,
 double                      ***projected_coord,
 int                         ***bar_coord_idx,
 double                      ***bar_coord,
 double                      ***uvw
 )
{
  int n_section = PDM_part_mesh_nodal_elmts_n_section_get(pmne);

  int *sections_id = PDM_part_mesh_nodal_elmts_sections_id_get(pmne);

  *distance        = malloc(sizeof(double *) * n_part);
  *projected_coord = malloc(sizeof(double *) * n_part);
  *bar_coord_idx   = malloc(sizeof(int    *) * n_part);
  *bar_coord       = malloc(sizeof(double *) * n_part);
  *uvw             = malloc(sizeof(double *) * n_part);

  /* First loop to allocate */
  for (int ipart = 0; ipart < n_part; ipart++) {

    int n_cell = 0;
    for (int isection = 0; isection < n_section; isection++) {
      int n_elt = PDM_part_mesh_nodal_elmts_section_n_elt_get(pmne,
                                                            sections_id[isection],
                                                            ipart);

      n_cell += n_elt;
    }

    int n_pts = pts_idx[ipart][n_cell];

    (*distance       )[ipart] = malloc (sizeof(double) * n_pts);
    (*projected_coord)[ipart] = malloc (sizeof(double) * n_pts * 3);
    (*uvw            )[ipart] = malloc (sizeof(double) * n_pts * 3);
    (*bar_coord_idx  )[ipart] = PDM_array_zeros_int(n_pts+1);

    int *_bar_coord_idx = (*bar_coord_idx)[ipart];

    for (int isection = 0; isection < n_section; isection++) {

      int id_section = sections_id[isection];

      int n_elt = PDM_part_mesh_nodal_elmts_section_n_elt_get(pmne,
                                                            id_section,
                                                            ipart);

      PDM_Mesh_nodal_elt_t t_elt = PDM_part_mesh_nodal_elmts_section_type_get(pmne,
                                                                            id_section);

      int *parent_num = PDM_part_mesh_nodal_elmts_parent_num_get(pmne,
                                                                 id_section,
                                                                 ipart,
                                                                 PDM_OWNERSHIP_KEEP);

      if (t_elt == PDM_MESH_NODAL_POLY_2D) {
        /* Polygonal section */
        int *connec_idx;
        int *connec;
        PDM_part_mesh_nodal_elmts_section_poly2d_get(pmne,
                                                   id_section,
                                                   ipart,
                                                   &connec_idx,
                                                   &connec,
                                                   PDM_OWNERSHIP_KEEP);

        for (int ielt = 0; ielt < n_elt; ielt++) {
          int icell = ielt;
          if (parent_num != NULL) {
            icell = parent_num[ielt];
          }

          int n_vtx = connec_idx[ielt+1] - connec_idx[ielt];

          for (int idx_pt = pts_idx[ipart][icell]; idx_pt < pts_idx[ipart][icell+1]; idx_pt++) {
            _bar_coord_idx[idx_pt+1] = n_vtx;
          }
        }
      }

      else if (t_elt == PDM_MESH_NODAL_POLY_3D) {
        /* Polyhedral section */
        int *connec_idx;
        int *connec;
        PDM_part_mesh_nodal_elmts_section_poly3d_cell_vtx_connect_get(pmne,
                                                                    id_section,
                                                                    ipart,
                                                                    &connec_idx,
                                                                    &connec,
                                                                    PDM_OWNERSHIP_KEEP);

        for (int ielt = 0; ielt < n_elt; ielt++) {
          int icell = ielt;
          if (parent_num != NULL) {
            icell = parent_num[ielt];
          }

          int n_vtx = connec_idx[ielt+1] - connec_idx[ielt];
          for (int idx_pt = pts_idx[ipart][icell]; idx_pt < pts_idx[ipart][icell+1]; idx_pt++) {
            _bar_coord_idx[idx_pt+1] = n_vtx;
          }
        }
      }

      else {
        /* Standard section */
              int         *connec              = NULL;
              PDM_g_num_t *numabs              = NULL;
              int         *_parent_num         = NULL;
              PDM_g_num_t *parent_entity_g_num = NULL;
              int          order               = 0;
        const char        *ho_ordering         = NULL;
        PDM_part_mesh_nodal_elmts_section_std_ho_get(pmne,
                                                   id_section,
                                                   ipart,
                                                   &connec,
                                                   &numabs,
                                                   &_parent_num,
                                                   &parent_entity_g_num,
                                                   &order,
                                                   &ho_ordering,
                                                   PDM_OWNERSHIP_KEEP);

        int n_vtx = PDM_Mesh_nodal_n_vtx_elt_get(t_elt,
                                                 order);

        for (int ielt = 0; ielt < n_elt; ielt++) {
          int icell = ielt;
          if (parent_num != NULL) {
            icell = parent_num[ielt];
          }

          for (int idx_pt = pts_idx[ipart][icell]; idx_pt < pts_idx[ipart][icell+1]; idx_pt++) {
            _bar_coord_idx[idx_pt+1] = n_vtx;
          }
        }
      }

    } // End of loop on current part's sections

    for (int i = 0; i < n_pts; i++) {
      _bar_coord_idx[i+1] += _bar_coord_idx[i];
    }
    (*bar_coord)[ipart] = malloc(sizeof(double) * _bar_coord_idx[n_pts]);

  } // End of loop on parts



  /* Second loop to perform location */
  for (int ipart = 0; ipart < n_part; ipart++) {

    double *_distance        = (*distance)       [ipart];
    double *_projected_coord = (*projected_coord)[ipart];
    int    *_bar_coord_idx   = (*bar_coord_idx)  [ipart];
    double *_bar_coord       = (*bar_coord)      [ipart];
    double *_uvw             = (*uvw)            [ipart];

    const double *vtx_coord = pvtx_coord[ipart];

    for (int isection = 0; isection < n_section; isection++) {

      int id_section = sections_id[isection];

      int n_elt = PDM_part_mesh_nodal_elmts_section_n_elt_get(pmne,
                                                            id_section,
                                                            ipart);

      PDM_Mesh_nodal_elt_t t_elt = PDM_part_mesh_nodal_elmts_section_type_get(pmne,
                                                                            id_section);

      int *parent_num = PDM_part_mesh_nodal_elmts_parent_num_get(pmne,
                                                                 id_section,
                                                                 ipart,
                                                                 PDM_OWNERSHIP_KEEP);

      if (t_elt == PDM_MESH_NODAL_POLY_2D) {
        /* Polygonal section */
        int *connec_idx;
        int *connec;
        PDM_part_mesh_nodal_elmts_section_poly2d_get(pmne,
                                                   id_section,
                                                   ipart,
                                                   &connec_idx,
                                                   &connec,
                                                   PDM_OWNERSHIP_KEEP);

        int n_vtx_max = 0;
        for (int ielt = 0; ielt < n_elt; ielt++) {
          n_vtx_max = PDM_MAX(n_vtx_max, connec_idx[ielt+1] - connec_idx[ielt]);
        }

        double *poly_coord = malloc(sizeof(double) * n_vtx_max * 3);

        for (int ielt = 0; ielt < n_elt; ielt++) {
          int icell = ielt;
          if (parent_num != NULL) {
            icell = parent_num[ielt];
          }

          int n_vtx = connec_idx[ielt+1] - connec_idx[ielt];
          int *_connec = connec + connec_idx[ielt];
          for (int i = 0; i < n_vtx; i++) {
            int vtx_id = _connec[i] - 1;
            memcpy(poly_coord + 3*i, vtx_coord + 3*vtx_id, sizeof(double) * 3);
          }

          int idx_pt = pts_idx[ipart][icell];
          int _n_pts = pts_idx[ipart][icell+1] - idx_pt;

          _locate_in_polygon(n_vtx,
                             poly_coord,
                             _n_pts,
            (const double *) pts_coord[ipart] + idx_pt * 3,
                             _distance        + idx_pt,
                             _projected_coord + idx_pt * 3,
                             _bar_coord       + _bar_coord_idx[idx_pt]);

          for (int ipt = idx_pt; ipt < pts_idx[ipart][icell+1]; ipt++) {
            _uvw[3*ipt] = _uvw[3*ipt+1] = _uvw[3*ipt+2] = -1.;
          }
        } // End of loop on polygons
        free(poly_coord);

      } // end if PDM_MESH_NODAL_POLY_2D


      else if (t_elt == PDM_MESH_NODAL_POLY_3D) {
        /* Polyhedral section */
        int *cell_vtx_idx;
        int *cell_vtx;
        PDM_part_mesh_nodal_elmts_section_poly3d_cell_vtx_connect_get(pmne,
                                                                    id_section,
                                                                    ipart,
                                                                    &cell_vtx_idx,
                                                                    &cell_vtx,
                                                                    PDM_OWNERSHIP_KEEP);

        int *cell_face_idx;
        int *cell_face;
        int  n_face;
        int *face_vtx_idx;
        int *face_vtx;
        PDM_g_num_t *face_ln_to_gn            = NULL;
        int         *_parent_num              = NULL;
        PDM_g_num_t *numabs                   = NULL;
        PDM_g_num_t *parent_entitity_ln_to_gn = NULL;
        PDM_part_mesh_nodal_elmts_section_poly3d_get(pmne,
                                                   id_section,
                                                   ipart,
                                                   &n_face,
                                                   &face_ln_to_gn,
                                                   &face_vtx_idx,
                                                   &face_vtx,
                                                   &numabs,
                                                   &cell_face_idx,
                                                   &cell_face,
                                                   &_parent_num,
                                                   &parent_entitity_ln_to_gn,
                                                   PDM_OWNERSHIP_KEEP);

        int pn_vtx = 0;
        for (int i = 0; i < cell_vtx_idx[n_elt]; i++) {
          pn_vtx = PDM_MAX(pn_vtx, cell_vtx[i]);
        }
        int *_vtx_id = PDM_array_zeros_int(pn_vtx);

        int n_vtx_max      = 0;
        int n_face_max     = 0;
        int s_face_vtx_max = 0;
        for (int ielt = 0; ielt < n_elt; ielt++) {
          n_vtx_max  = PDM_MAX(n_vtx_max,  cell_vtx_idx [ielt+1] - cell_vtx_idx [ielt]);
          n_face_max = PDM_MAX(n_face_max, cell_face_idx[ielt+1] - cell_face_idx[ielt]);

          int s_face_vtx = 0;
          for (int idx = cell_face_idx[ielt]; idx < cell_face_idx[ielt+1]; idx++) {
            int face_id = PDM_ABS(cell_face[idx]) - 1;
            s_face_vtx += face_vtx_idx[face_id+1] - face_vtx_idx[face_id];
          }

          s_face_vtx_max = PDM_MAX(s_face_vtx_max, s_face_vtx);
        }

        int *_face_orientation = malloc(sizeof(int) * n_face_max);
        int *_face_vtx_idx = malloc(sizeof(int) * (n_face_max + 1));
        _face_vtx_idx[0] = 0;

        int *_face_vtx = malloc(sizeof(int) * s_face_vtx_max);
        double *_vtx_coord = malloc(sizeof(double) * n_vtx_max * 3);

        PDM_triangulate_state_t *_tri_state = NULL;
        _tri_state = PDM_triangulate_state_create(s_face_vtx_max);
        int n_tri_vtx = PDM_MAX((s_face_vtx_max - 2) * 3, 0);
        int *_tri_vtx = malloc(sizeof(int) * n_tri_vtx);


        for (int ielt = 0; ielt < n_elt; ielt++) {
          int icell = ielt;
          if (parent_num != NULL) {
            icell = parent_num[ielt];
          }

          int _n_vtx  = cell_vtx_idx [ielt+1] - cell_vtx_idx [ielt];
          int _n_face = cell_face_idx[ielt+1] - cell_face_idx[ielt];

          for (int ivtx = 0; ivtx < _n_vtx; ivtx++) {
            int vtx_id = cell_vtx[cell_vtx_idx[ielt] + ivtx] - 1;
            _vtx_id[vtx_id] = ivtx+1;
            memcpy(_vtx_coord + 3*ivtx, pvtx_coord[ipart] + 3*vtx_id, sizeof(double) * 3);
          }

          for (int iface = 0; iface < _n_face; iface++) {
            int idx_face = cell_face_idx[ielt] + iface;
            int face_id = PDM_ABS(cell_face[idx_face]) - 1;
            _face_orientation[iface] = PDM_SIGN(cell_face[idx_face]);
            _face_vtx_idx[iface+1] = _face_vtx_idx[iface];
            for (int idx_vtx = face_vtx_idx[face_id]; idx_vtx < face_vtx_idx[face_id+1]; idx_vtx++) {
              int vtx_id = face_vtx[idx_vtx] - 1;
              assert(_vtx_id[vtx_id] > 0);
              _face_vtx[_face_vtx_idx[iface+1]++] = _vtx_id[vtx_id];
            }
          }

          int idx_pt = pts_idx[ipart][icell];
          int _n_pts = pts_idx[ipart][icell+1] - idx_pt;

          _locate_in_polyhedron(_n_vtx,
                (const double *) _vtx_coord,
                                 _n_face,
                                 _face_vtx_idx,
                (const int    *) _face_vtx,
                (const int    *) _face_orientation,
                                 _n_pts,
                (const double *) pts_coord[ipart]  + idx_pt * 3,
                                 _distance         + idx_pt,
                                 _projected_coord  + idx_pt * 3,
                                 _bar_coord        + _bar_coord_idx[idx_pt],
                                 _tri_state,
                                 _tri_vtx);

          for (int ipt = idx_pt; ipt < pts_idx[ipart][icell+1]; ipt++) {
            _uvw[3*ipt] = _uvw[3*ipt+1] = _uvw[3*ipt+2] = -1.;
          }

          // Reset local vtx_id
          for (int idx = cell_vtx_idx[ielt]; idx < cell_vtx_idx[ielt+1]; idx++) {
            int vtx_id = cell_vtx[idx] - 1;
            _vtx_id[vtx_id] = 0;
          }
        } // End of loop on polyhedra

        free(_vtx_id);
        free(_face_orientation);
        free(_face_vtx_idx);
        free(_face_vtx);
        free(_vtx_coord);
        free(_tri_vtx);
        PDM_triangulate_state_destroy(_tri_state);

      } // end if PDM_MESH_NODAL_POLY_3D


      else {
        /* Standard section */
              int         *connec              = NULL;
              PDM_g_num_t *numabs              = NULL;
              int         *_parent_num         = NULL;
              PDM_g_num_t *parent_entity_g_num = NULL;
              int          order               = 0;
        const char        *ho_ordering         = NULL;
        PDM_part_mesh_nodal_elmts_section_std_ho_get(pmne,
                                                   id_section,
                                                   ipart,
                                                   &connec,
                                                   &numabs,
                                                   &_parent_num,
                                                   &parent_entity_g_num,
                                                   &order,
                                                   &ho_ordering,
                                                   PDM_OWNERSHIP_KEEP);


        switch (t_elt) {

          case PDM_MESH_NODAL_POINT: {
            for (int ielt = 0; ielt < n_elt; ielt++) {
              int icell = ielt;
              if (parent_num != NULL) {
                icell = parent_num[ielt];
              }

              const double *node_coord = pvtx_coord[ipart] + 3*(connec[ielt] - 1);
              for (int idx_pt = pts_idx[ipart][icell]; idx_pt < pts_idx[ipart][icell+1]; idx_pt++) {
                double dist2 = 0.;
                for (int j = 0; j < 3; j++) {
                  double dj = pts_coord[ipart][3*idx_pt+j] - node_coord[j];
                  dist2 += dj * dj;
                }
                _distance [idx_pt] = dist2;
                _bar_coord[_bar_coord_idx[idx_pt]] = 1.;
                _uvw      [3*idx_pt  ] = 0.;
                _uvw      [3*idx_pt+1] = -1.;
                _uvw      [3*idx_pt+2] = -1.;
              }
            }
            break;
          } // end case PDM_MESH_NODAL_POINT


          case PDM_MESH_NODAL_BAR2: {
            double edge_coord[6];
            for (int ielt = 0; ielt < n_elt; ielt++) {
              int icell = ielt;
              if (parent_num != NULL) {
                icell = parent_num[ielt];
              }

              int *_connec = connec + 2*ielt;
              for (int i = 0; i < 2; i++) {
                int vtx_id = _connec[i] - 1;
                memcpy(edge_coord + 3*i, vtx_coord + 3*vtx_id, sizeof(double) * 3);
              }

              int idx_pt = pts_idx[ipart][icell];

              _locate_on_edge(edge_coord,
                              pts_idx[ipart][icell+1] - idx_pt,
             (const double *) pts_coord[ipart] + idx_pt * 3,
                              _distance        + idx_pt,
                              _projected_coord + idx_pt * 3,
                              _bar_coord       + _bar_coord_idx[idx_pt]);

              for (int ipt = idx_pt; ipt < pts_idx[ipart][icell+1]; ipt++) {
                _uvw[3*ipt  ] = _bar_coord[_bar_coord_idx[ipt]+1];
                _uvw[3*ipt+1] = -1.;
                _uvw[3*ipt+2] = -1.;
              }

            } // End of loop on edges
            break;
          } // end case PDM_MESH_NODAL_BAR2


          case PDM_MESH_NODAL_TRIA3: {
            for (int ielt = 0; ielt < n_elt; ielt++) {
              int icell = ielt;
              if (parent_num != NULL) {
                icell = parent_num[ielt];
              }

              int *_connec = connec + 3*ielt;
              int idx_pt = pts_idx[ipart][icell];

              _locate_on_triangles(1,
                                   _connec,
                                   pvtx_coord[ipart],
                                   pts_idx[ipart][icell+1] - idx_pt,
                  (const double *) pts_coord[ipart] + idx_pt * 3,
                                   NULL,
                                   _distance        + idx_pt,
                                   _projected_coord + idx_pt * 3,
                                   _bar_coord       + _bar_coord_idx[idx_pt]);

              for (int ipt = idx_pt; ipt < pts_idx[ipart][icell+1]; ipt++) {
                _uvw[3*ipt  ] = _bar_coord[_bar_coord_idx[ipt]+1];
                _uvw[3*ipt+1] = _bar_coord[_bar_coord_idx[ipt]+2];
                _uvw[3*ipt+2] = -1.;
              }

            } // End of loop on triangles

            break;
          } // end case PDM_MESH_NODAL_TRIA3


          case PDM_MESH_NODAL_QUAD4: {
            double quad_coord[12];
            for (int ielt = 0; ielt < n_elt; ielt++) {
              int icell = ielt;
              if (parent_num != NULL) {
                icell = parent_num[ielt];
              }

              int *_connec = connec + 4*ielt;
              for (int i = 0; i < 4; i++) {
                int vtx_id = _connec[i] - 1;
                memcpy(quad_coord + 3*i, vtx_coord + 3*vtx_id, sizeof(double) * 3);
              }

              int idx_pt = pts_idx[ipart][icell];

              // _locate_on_quadrangle(quad_coord,
              //                       pts_idx[ipart][icell+1] - idx_pt,
              //      (const double *) pts_coord[ipart] + idx_pt * 3,
              //                       _distance        + idx_pt,
              //                       _projected_coord + idx_pt * 3,
              //                       _bar_coord       + _bar_coord_idx[idx_pt]);
              _locate_on_quadrangle(quad_coord,
                                     pts_idx[ipart][icell+1] - idx_pt,
                    (const double *) pts_coord[ipart] + idx_pt * 3,
                                     tolerance,
                                     _distance        + idx_pt,
                                     _projected_coord + idx_pt * 3,
                                     _bar_coord       + _bar_coord_idx[idx_pt],
                                     _uvw             + idx_pt * 3);

              // // TODO: compute uvw
              // for (int ipt = idx_pt; ipt < pts_idx[ipart][icell+1]; ipt++) {
              //   _uvw[3*ipt  ] = -1;
              //   _uvw[3*ipt+1] = -1;
              //   _uvw[3*ipt+2] = -1;
              // }

            } // End of loop on quadrangles
            break;
          } // end case PDM_MESH_NODAL_QUAD4

          case PDM_MESH_NODAL_TETRA4: {
            double tetra_coord[12];
            for (int ielt = 0; ielt < n_elt; ielt++) {
              int icell = ielt;
              if (parent_num != NULL) {
                icell = parent_num[ielt];
              }

              int *_connec = connec + 4*ielt;
              for (int i = 0; i < 4; i++) {
                int vtx_id = _connec[i] - 1;
                memcpy(tetra_coord + 3*i, vtx_coord + 3*vtx_id, sizeof(double) * 3);
              }

              int idx_pt = pts_idx[ipart][icell];

              _locate_in_tetrahedron(tetra_coord,
                                     pts_idx[ipart][icell+1] - idx_pt,
                    (const double *) pts_coord[ipart] + idx_pt * 3,
                                     _distance        + idx_pt,
                                     _projected_coord + idx_pt * 3,
                                     _bar_coord       + _bar_coord_idx[idx_pt]);

              for (int ipt = idx_pt; ipt < pts_idx[ipart][icell+1]; ipt++) {
                _uvw[3*ipt  ] = _bar_coord[_bar_coord_idx[ipt]+1];
                _uvw[3*ipt+1] = _bar_coord[_bar_coord_idx[ipt]+2];
                _uvw[3*ipt+2] = _bar_coord[_bar_coord_idx[ipt]+3];
              }

            } // End of loop on tetrahedra
            break;
          } // end case PDM_MESH_NODAL_TETRA4


          case PDM_MESH_NODAL_PYRAMID5:
          case PDM_MESH_NODAL_PRISM6  :
          case PDM_MESH_NODAL_HEXA8   : {
            double elt_coord[24];
            int n_vtx = PDM_Mesh_nodal_n_vtx_elt_get(t_elt,
                                                     order);
            for (int ielt = 0; ielt < n_elt; ielt++) {
              int icell = ielt;
              if (parent_num != NULL) {
                icell = parent_num[ielt];
              }

              int *_connec = connec + n_vtx*ielt;
              for (int i = 0; i < n_vtx; i++) {
                int vtx_id = _connec[i] - 1;
                memcpy(elt_coord + 3*i, vtx_coord + 3*vtx_id, sizeof(double) * 3);
              }

              int idx_pt = pts_idx[ipart][icell];

              _locate_in_cell_3d(t_elt,
                                 elt_coord,
                                 pts_idx[ipart][icell+1] - idx_pt,
                (const double *) pts_coord[ipart] + idx_pt * 3,
                                 tolerance,
                                 _distance        + idx_pt,
                                 _projected_coord + idx_pt * 3,
                                 _bar_coord       + _bar_coord_idx[idx_pt],
                                 _uvw             + idx_pt * 3);

            } // End of loop on pyramids/prisms/hexahedra
            break;
          } // end case PDM_MESH_NODAL_PYRAMID5/PRISM6/HEXA8

          case PDM_MESH_NODAL_BARHO    :
          case PDM_MESH_NODAL_TRIAHO   :
          case PDM_MESH_NODAL_QUADHO   :
          case PDM_MESH_NODAL_TETRAHO  :
          case PDM_MESH_NODAL_PYRAMIDHO:
          case PDM_MESH_NODAL_PRISMHO  :
          case PDM_MESH_NODAL_HEXAHO   : {

            int n_node = PDM_Mesh_nodal_n_vtx_elt_get(t_elt,
                                                      order);
            double *elt_coord = malloc(sizeof(double) * n_node * 3);

            int elt_dim = PDM_Mesh_nodal_elt_dim_get(t_elt);

            double *work_array = malloc(sizeof(double) * n_node * (elt_dim+1));

            int *ijk_to_user = NULL;
            if (ho_ordering != NULL) {
              ijk_to_user = PDM_ho_ordering_ijk_to_user_get(ho_ordering,
                                                            t_elt,
                                                            order);
              assert(ijk_to_user != NULL);
            }

            int count           = 0;
            int count_converged = 0;
            for (int ielt = 0; ielt < n_elt; ielt++) {
              int icell = ielt;
              if (parent_num != NULL) {
                icell = parent_num[ielt];
              }
              // log_trace("icell = %d\n", icell);

              int *_connec = connec + n_node*ielt;
              for (int idx_ijk = 0; idx_ijk < n_node; idx_ijk++) {
                int idx_user = idx_ijk;
                if (ijk_to_user != NULL) {
                  idx_user = ijk_to_user[idx_ijk];
                }
                int node_id = _connec[idx_user] - 1;
                memcpy(elt_coord + 3*idx_ijk, vtx_coord + 3*node_id, sizeof(double) * 3);
              }

              for (int idx_pt = pts_idx[ipart][icell]; idx_pt < pts_idx[ipart][icell+1]; idx_pt++) {
                count++;
                // First, try Newton method
                int converged = 0;
                if (1) {
                  _distance[idx_pt] = PDM_ho_location_newton(t_elt,
                                                             order,
                                                             n_node,
                                                             elt_coord,
                                            (const double *) pts_coord[ipart] + idx_pt * 3,
                                                             tolerance,
                                                             _projected_coord + idx_pt * 3,
                                                             _uvw             + idx_pt * 3,
                                                             &converged,
                                                             work_array);
                }
                // log_trace("converged? %d\n", converged);

                if (converged) {
                  count_converged++;
                  // log_trace("Newton converged :D (%f %f %f)\n",
                  //           pts_coord[ipart][3*idx_pt+0],
                  //           pts_coord[ipart][3*idx_pt+1],
                  //           pts_coord[ipart][3*idx_pt+2]);
                }
                else {
                  // log_trace("Newton HO failed to converge! (%f %f %f)\n",
                  //           pts_coord[ipart][3*idx_pt+0],
                  //           pts_coord[ipart][3*idx_pt+1],
                  //           pts_coord[ipart][3*idx_pt+2]);
                  // Newton failed to converge, try more robust method
                  _distance[idx_pt] = PDM_ho_location(t_elt,
                                                      order,
                                                      n_node,
                                                      elt_coord,
                                     (const double *) pts_coord[ipart] + idx_pt * 3,
                                                      _projected_coord + idx_pt * 3,
                                                      _uvw             + idx_pt * 3);
                }
                // log_trace("dist2 = %e\n", _distance[idx_pt]);

              } // End of loop on current elt's points

              /* /!\ uvw is stride-3, even for lower-dimensional elements */
              int idx_pt0 = pts_idx[ipart][icell];
              for (int i = idx_pt0; i < pts_idx[ipart][icell+1]; i++) {
                PDM_ho_basis(t_elt,
                             order,
                             n_node,
                             1,
                             (const double *) _uvw + i * 3,
                             _bar_coord + _bar_coord_idx[i]);
              }

            } // End of loop on elt
            // log_trace("%d converged / %d\n", count_converged, count);
            free(elt_coord);
            free(work_array);
            break;
          } // end case HO LAGRANGE


          // HO_BEZIER...

          default: {
            PDM_error(__FILE__, __LINE__, 0,
                      "PDM_point_location_nodal : element type %d not supported yet\n",
                      (int) t_elt);
          }

        } // End switch on std elt type


      } // End if standard section

    } // End of loop on current part's sections

  } // End of loop on parts

}


/**
 * \brief Compute quadrangle, hexahedron, pyramid, or prism parametric coordinates for a
 * given point.
 *
 * This function is adapted from the CGNS interpolation tool.
 *
 * \param [in]   elt_type       Type of element
 * \param [in]   point_coords   Point coordinates
 * \param [in]   vertex_coords  Pointer to element vertex coordinates
 * \param [in]   tolerance      Location tolerance factor
 * \param [out]  uvw            Parametric coordinates of point in element
 * \param [in]   init_uvw       Initial uvw guess for Newton method (or NULL)
 *
 *  \return Convergence status of Newton method
 */

PDM_bool_t
PDM_point_location_compute_uvw
(
 const PDM_Mesh_nodal_elt_t elt_type,
 const double               point_coords[3],
 const double               vertex_coords[],
 const double               tolerance,
       double               uvw[3],
       double               init_uvw[3]
 )
{
  return _compute_uvw (elt_type,
                       point_coords,
                       vertex_coords,
                       tolerance,
                       uvw,
                       init_uvw);
}

#ifdef __cplusplus
}
#endif /* __cplusplus */
