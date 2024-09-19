#include <math.h>
#include <sys/time.h>
#include <time.h>
#include <sys/resource.h>
#include <unistd.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>


#include "pdm.h"
#include "pdm_config.h"
#include "pdm_mpi.h"
#include "pdm_priv.h"
#include "pdm_error.h"
#include "pdm_printf.h"
#include "pdm_logging.h"
#include "pdm_vtk.h"
#include "pdm_sort.h"
#include "pdm_line.h"
#include "pdm_triangle.h"
#include "pdm_predicate.h"

/*============================================================================
 * Macro definitions
 *============================================================================*/

/*============================================================================
 * Type definitions
 *============================================================================*/

/* point classification wrt a triangle */
typedef enum {
  LOC_ON_VTX0           = 0,
  LOC_ON_VTX1           = 1,
  LOC_ON_VTX2           = 2,
  LOC_INTERIOR_EDGE0    = 3,
  LOC_INTERIOR_EDGE1    = 4,
  LOC_INTERIOR_EDGE2    = 5,
  LOC_INTERIOR_FACE     = 6,
  LOC_EXTERIOR_ON_PLANE = 7,
  LOC_EXTERIOR_ABOVE    = 8,
  LOC_EXTERIOR_BELOW    = 9
} _point_triangle_loc_t;


typedef struct {
  double coordA[3];
  double coordB[3];

  PDM_mesh_entities_t entity_typeA;
  PDM_mesh_entities_t entity_typeB;

  int entity_idA;
  int entity_idB;
} _intersection_point_t;

/*============================================================================
 * Private functions
 *============================================================================*/

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

/* Compute squared euclidean distance between point p and triangle abc */
static inline double
_dist2_point_triangle
(
 const double  p[3],
 const double  a[3],
 const double  b[3],
 const double  c[3],
       double *u,
       double *v
 )
{
  double ab[3] = {b[0] - a[0], b[1] - a[1], b[2] - a[2]};
  double ac[3] = {c[0] - a[0], c[1] - a[1], c[2] - a[2]};
  double ap[3] = {p[0] - a[0], p[1] - a[1], p[2] - a[2]};

  double abab = PDM_DOT_PRODUCT(ab, ab);
  double abac = PDM_DOT_PRODUCT(ab, ac);
  double acac = PDM_DOT_PRODUCT(ac, ac);
  double abap = PDM_DOT_PRODUCT(ab, ap);
  double acap = PDM_DOT_PRODUCT(ac, ap);

  double det = abab*acac - abac*abac;

  double _u = abap*acac - acap*abac;
  double _v = acap*abab - abap*abac;

  if (_u >= 0 && _u <= det &&
      _v >= 0 && _v <= det &&
      _u + _v <= det) {
    if (det >= 0.0) {
      double idet = 1./det;
      _u *= idet;
      _v *= idet;
    }
    *u = _u;
    *v = _v;
    double q[3] = {
      a[0] + _u*ab[0] + _v*ac[0],
      a[1] + _u*ab[1] + _v*ac[1],
      a[2] + _u*ab[2] + _v*ac[2]
    };
    return _dist2_point_point(p, q);
  }

  const double *x[3] = {a, b, c};
  double dist2 = HUGE_VAL;
  int    imin = -1;
  double tmin = 0;
  for (int i = 0; i < 3; i++) {
    double t;
    double d = _dist2_point_segment(p,
                                    x[i],
                                    x[(i+1)%3],
                                    &t);
    if (d < dist2) {
      dist2 = d;
      tmin  = t;
      imin  = i;
    }
  }

  if (imin == 0) {
    *u = tmin;
    *v = 0;
  }
  else if (imin == 1) {
    *u = 1 - tmin;
    *v = tmin;
  }
  else {
    *u = 0;
    *v = 1 - tmin;
  }

  return dist2;
}

static inline double
_rand01
(
 void
 )
{
  return (double) rand() / (double) RAND_MAX;
}

static const double isqrt3 = 0.5773502691; // slightly < 1/sqrt(3)

static void
_gen_config
(
 double                x[3][3],
 double                r[3],
 _point_triangle_loc_t loc[3],
 double                eps[3],
 double                noise_scale,
 double                coord[3][3]
 )
{
  int max_try = 20;

  for (int i = 0; i < 3; i++) {

    // generate a random vector with magnitude <= 1
    double noise[3];
    for (int j = 0; j < 3; j++) {
      noise[j] = 2*_rand01() - 1.;
      noise[j] *= isqrt3;
    }

    if (loc[i] < 3) {
      // on vtx
      int ivtx = loc[i];
      double e = PDM_MIN(r[ivtx], eps[i]);
      for (int j = 0; j < 3; j++) {
        coord[i][j] = x[ivtx][j] + noise_scale*e*noise[j];
      }
    }

    else if (loc[i] < 6) {
      // on edge
      int iedge = loc[i] - 3;

      for (int itry = 0; itry < max_try; itry++) {
        double t = _rand01();

        double e = PDM_MIN((1-t)*r[iedge] + t*r[(iedge+1)%3], eps[i]);
        for (int j = 0; j < 3; j++) {
          coord[i][j] = (1-t)*x[iedge][j] + t*x[(iedge+1)%3][j] + noise_scale*e*noise[j];
        }

        double d2 = noise_scale*e;
        d2 *= d2 * PDM_DOT_PRODUCT(noise, noise);

        int keep = 1;
        for (int j = 0; j < 2; j++) {
          double tj;
          double dj = _dist2_point_segment(coord[i], x[(iedge+j+1)%3], x[(iedge+j+2)%3], &tj);
          double ej = PDM_MIN((1-tj)*r[(iedge+j+1)%3] + tj*r[(iedge+j+2)%3], eps[i]);
          if (dj < ej*ej && dj < d2) {// apply Eric's rule for points close to adjacent edges?
            keep = 0;
            // log_trace("edge %d try %d failed (edge)\n", iedge, itry, j);
            break;
          }
        }

        if (!keep) {
          continue;
        }

        for (int j = 0; j < 2; j++) {
          double dj = _dist2_point_point(coord[i], x[(iedge+j)%3]);
          double ej = PDM_MIN(r[(iedge+j)%3], eps[i]);
          if (dj < ej*ej) {
            keep = 0;
            // log_trace("edge %d try %d failed (vtx %d)\n", iedge, itry, j);
            break;
          }
        }

        if (keep) {
          break;
        }
      }

    }

    else {

      for (int itry = 0; itry < max_try; itry++) {

        int keep = 1;
        double u = _rand01();
        double v = _rand01();
        double w = _rand01();
        double isum = 1./(u + v + w);
        u *= isum;
        v *= isum;
        w *= isum;

        for (int j = 0; j < 3; j++) {
          coord[i][j] = w*x[0][j] + u*x[1][j] + v*x[2][j];
        }

        double e = PDM_MIN(w*r[0] + u*r[1] + v*r[2], eps[i]);
        if (loc[i] == LOC_INTERIOR_FACE) {
          for (int j = 0; j < 3; j++) {
            coord[i][j] += noise_scale*e*noise[j];
          }

          for (int j = 0; j < 3; j++) {
              double tj;
              double dj = _dist2_point_segment(coord[i], x[j], x[(j+1)%3], &tj);
              double ej = PDM_MIN((1-tj)*r[j] + tj*r[(j+1)%3], eps[i]);
            if (dj < ej*ej) {
              keep = 0;
              // log_trace("face try %d failed (edge %d)\n", itry, j);
              break;
            }
          }

        }
        else {
          if (loc[i] == LOC_EXTERIOR_ON_PLANE) {
            double dmin = HUGE_VAL;
            double tmin;
            int    jmin = -1;
            for (int j = 0; j < 3; j++) {
              double t;
              double d = _dist2_point_segment(coord[i], x[j], x[(j+1)%3], &t);
              if (d < dmin) {
                dmin = d;
                tmin = t;
                jmin = j;
              }
            }
            if (dmin < 1e-32) {
              keep = 0;
              continue;
            }
            double vec[3] = {
              (1-tmin)*x[jmin][0] + tmin*x[(jmin+1)%3][0] - coord[i][0],
              (1-tmin)*x[jmin][1] + tmin*x[(jmin+1)%3][1] - coord[i][1],
              (1-tmin)*x[jmin][2] + tmin*x[(jmin+1)%3][2] - coord[i][2]
            };
            e = PDM_MAX(eps[i], (1-tmin)*r[jmin] + tmin*r[(jmin+1)%3]);
            double scale = 1 + PDM_MAX(1, 20*e/sqrt(dmin));
            for (int j = 0; j < 3; j++) {
              coord[i][j] += scale*vec[j];
            }
          }

          else {
            double xu[3] = {x[1][0] - x[0][0], x[1][1] - x[0][1], x[1][2] - x[0][2]};
            double xv[3] = {x[2][0] - x[0][0], x[2][1] - x[0][1], x[2][2] - x[0][2]};
            double normal[3];
            PDM_CROSS_PRODUCT(normal, xu, xv);

            double sign = PDM_SIGN(PDM_DOT_PRODUCT(normal, noise));
            double epsn = eps[i] / PDM_MODULE(normal);
            for (int j = 0; j < 3; j++) {
              noise[j] += sign*normal[j]*epsn;
            }

            double scale = sign / PDM_MODULE(noise);
            scale *= PDM_MAX(20*PDM_MAX(e, eps[i]), 1);
            if (loc[i] == LOC_EXTERIOR_BELOW) {
              scale = -scale;
            }
            for (int j = 0; j < 3; j++) {
              coord[i][j] += scale*noise[j];
            }
          }
        }

        if (keep) {
          break;
        }

      }

    }

  }
}



/**
 * \brief Compute classification of a set of (fat) points wrt a (fat) triangle
 *
 * \param [in ] x          Coordinates of the triangle vertices
 * \param [in ] r          Tolerance of the triangle vertices
 * \param [in ] n_pts      Number of points to classify
 * \param [in ] pts_coord  Coordinates of the points to classify (size = 3*n_pts)
 * \param [in ] pts_eps    Tolerance of the points to classify (size = n_pts)
 * \param [out] pts_loc    Classification of points (size = n_pts)
 * \param [out] pts_det    Signed heigt of points wrt triangle plane (size = n_pts or NULL)
 * \param [in ] verbose    Verbosity (1 or 0)
 *
 */

static void
_classify_pts_wrt_triangle
(
       double                 x[3][3],
       double                 r[3],
 const int                    n_pts,
       double                *pts_coord,
       double                *pts_eps,
       _point_triangle_loc_t *pts_loc,
       double                *pts_det,
 const int                    verbose
 )
{
  double x01[3];
  double x02[3];
  double max_r2 = 0;

  for (int i = 0; i < 3; i++) {
    x01[i] = x[1][i] - x[0][i];
    x02[i] = x[2][i] - x[0][i];
    max_r2 = PDM_MAX(max_r2, r[i]);
  }
  max_r2 *= max_r2;

  // double normal[3];
  // PDM_CROSS_PRODUCT(normal, x01, x02);


  double x01x01 = PDM_DOT_PRODUCT(x01, x01);
  double x01x02 = PDM_DOT_PRODUCT(x01, x02);
  double x02x02 = PDM_DOT_PRODUCT(x02, x02);

  double det = x01x01*x02x02 - x01x02*x01x02;

  // TODO : handle degenerate case more gracefully...
  assert(det >= 0.);

  double idet = 1./det;

  for (int ipt = 0; ipt < n_pts; ipt++) {
    double *p = pts_coord + 3*ipt;
    double epsp2 = pts_eps[ipt]*pts_eps[ipt];

    // compute orthogonal projection of p on the triangle's supporting plane
    double x0p[3] = {p[0] - x[0][0], p[1] - x[0][1], p[2] - x[0][2]};

    double x0px01 = PDM_DOT_PRODUCT(x0p, x01);
    double x0px02 = PDM_DOT_PRODUCT(x0p, x02);

    // double _pt_det = PDM_DOT_PRODUCT(x0p, normal)*idet;
    // double _pt_det = PDM_DOT_PRODUCT(x0p, normal);
    double _pt_det = PDM_predicate_orient3d(x[0],
                                            x[1],
                                            x[2],
                                            p);

    if (verbose) {
      // log_trace("_pt_det**2 = %e / %e\n", _pt_det*_pt_det, det*pts_eps[ipt]*pts_eps[ipt]);
      log_trace("_pt_det = %e * sqrt(det)\n", _pt_det/sqrt(det));
    }

    // /!\ symmetry of "on plane" is broken...
    // Find a good epsilon to ensure symmetry?
    // if (PDM_ABS(_pt_det) < eps_on_plane) {
    // if (_pt_det*_pt_det < det*pts_eps[ipt]*pts_eps[ipt]) {
    //   _pt_det = 0.;
    //   pts_loc[ipt] = LOC_EXTERIOR_ON_PLANE;
    // }
    // else if (_pt_det > 0) {
    //   pts_loc[ipt] = LOC_EXTERIOR_ABOVE;
    // }
    // else {
    //   pts_loc[ipt] = LOC_EXTERIOR_BELOW;
    // }
    if (_pt_det > 0) {
      pts_loc[ipt] = LOC_EXTERIOR_ABOVE;
    }
    else if (_pt_det < 0) {
      pts_loc[ipt] = LOC_EXTERIOR_BELOW;
    }
    else {
      pts_loc[ipt] = LOC_EXTERIOR_ON_PLANE;
    }


    if (pts_det != NULL) {
      pts_det[ipt] = _pt_det;
    }

    // solve for parametric coordinates of projection
    double u = (x0px01*x02x02 - x0px02*x01x02)*idet;
    double v = (x0px02*x01x01 - x0px01*x01x02)*idet;
    double w = 1 - u - v;

    if (verbose) {
      log_trace("point %f %f %f, eps = %e\n", p[0], p[1], p[2], pts_eps[ipt]);
      log_trace("  ort proj: u = %f, v = %f, w = %f\n", u, v, w);
    }


    double q[3];
    for (int i = 0; i < 3; i++) {
      q[i] = w*x[0][i] + u*x[1][i] + v*x[2][i];
    }

    double d = _dist2_point_point(p, q);

    double eps = w*r[0] + u*r[1] + v*r[2];

    if (d > epsp2 || d > max_r2) {
      // p cannot be inside the fat triangle
      if (verbose) {
        log_trace("  p clearly exterior\n");
      }
      continue; // move on to the next point
    }

    if (u >= 0 && u <= 1 &&
        v >= 0 && v <= 1 &&
        w >= 0) {
      // projection lies inside the (closed) triangle
      if (d <= eps*eps) {
        // p is inside the fat triangle
        pts_loc[ipt] = LOC_INTERIOR_FACE;
        if (verbose) {
          log_trace("  p interior face\n");
        }
      }
    }

    // find closest point on the triangle boundary
    int    on_edge[3] = {0, 0, 0};
    double d_edge [3];
    for (int iedge = 0; iedge < 3; iedge++) {
      double t;
      d_edge[iedge] = _dist2_point_segment(p,
                                           x[iedge],
                                           x[(iedge+1)%3],
                                           &t);
      eps = (1-t)*r[iedge] + t*r[(iedge+1)%3];
      if (d_edge[iedge] <= epsp2 && d_edge[iedge] <= eps*eps) {
        // p is inside the fat edge
        on_edge[iedge] = 1;
        if (verbose) {
          log_trace("  p on edge %d\n", iedge);
        }
      }
    }

    // determine if p lies inside the fat boundary and, if so, on which entity
    if (on_edge[0]) {
      if (on_edge[1]) {
        if (on_edge[2]) {
          // weird but why not?
          int    on_vtx[3] = {0, 0, 0};
          double d_vtx [3];
          for (int ivtx = 0; ivtx < 3; ivtx++) {
            d_vtx[ivtx] = _dist2_point_point(p, x[ivtx]);
            if (d_vtx[ivtx] <= epsp2 && d_vtx[ivtx] < r[ivtx]*r[ivtx]) {
              on_vtx[ivtx] = 1;
            }
          }

          if (on_vtx[0]) {
            if (on_vtx[1]) {
              if (on_vtx[2]) {
                // WTF
                double dmin = HUGE_VAL;
                int    imin = -1;
                for (int ivtx = 0; ivtx < 3; ivtx++) {
                  if (d_vtx[ivtx] < dmin) {
                    dmin = d_vtx[ivtx];
                    imin = ivtx;
                  }
                }

                pts_loc[ipt] = (_point_triangle_loc_t) imin;
              } // end if on 3 vtx
              else if (d_vtx[0] < d_vtx[1]) {
                pts_loc[ipt] = LOC_ON_VTX0;
              }
              else {
                pts_loc[ipt] = LOC_ON_VTX1;
              }
            } // end if on vtx 0 and 1
            else if (on_vtx[2]) {
              if (d_vtx[0] < d_vtx[2]) {
                pts_loc[ipt] = LOC_ON_VTX0;
              }
              else {
                pts_loc[ipt] = LOC_ON_VTX2;
              }
            } // end if on vtx 0 and 2
            else {
              pts_loc[ipt] = LOC_ON_VTX0;
            }
          } // end if on vtx 0
          else if (on_vtx[1]) {
            if (on_vtx[2]) {
              if (d_vtx[1] < d_vtx[2]) {
                pts_loc[ipt] = LOC_ON_VTX1;
              }
              else {
                pts_loc[ipt] = LOC_ON_VTX2;
              }
            } // end if on vtx 1 and 2
            else {
              pts_loc[ipt] = LOC_ON_VTX1;
            }
          } // end if on vtx 1
          else if (on_vtx[2]) {
            pts_loc[ipt] = LOC_ON_VTX2;
          } // end if on vtx 2
          else {
            double dmin = HUGE_VAL;
            int    imin = -1;
            for (int iedge = 0; iedge < 3; iedge++) {
              if (d_edge[iedge] < dmin) {
                dmin = d_edge[iedge];
                imin = iedge;
              }
            }

            pts_loc[ipt] = (_point_triangle_loc_t) (3 + imin);
          }

        } // end if on 3 edges
        else {
          double d1 = _dist2_point_point(p, x[1]);
          if (d1 <= epsp2 && d1 < r[1]*r[1]) {
            // p is inside the fat vertex 1
            pts_loc[ipt] = LOC_ON_VTX1;
          }
          else if (d_edge[0] < d_edge[1]) { // Eric's rule?
            pts_loc[ipt] = LOC_INTERIOR_EDGE0;
          }
          else {
            pts_loc[ipt] = LOC_INTERIOR_EDGE1;
          }
        }
      } // end if on edges 0 and 1
      else if (on_edge[2]) {
        double d0 = _dist2_point_point(p, x[0]);
        if (d0 <= epsp2 && d0 < r[0]*r[0]) {
          // p is inside the fat vertex 0
          pts_loc[ipt] = LOC_ON_VTX0;
        }
        else if (d_edge[0] < d_edge[2]) { // Eric's rule?
          pts_loc[ipt] = LOC_INTERIOR_EDGE0;
        }
        else {
          pts_loc[ipt] = LOC_INTERIOR_EDGE2;
        }
      } // end if on edges 0 and 2
      else {
        pts_loc[ipt] = LOC_INTERIOR_EDGE0;
      }
    } // end if on edge 0
    else if (on_edge[1]) {
      if (on_edge[2]) {
        double d2 = _dist2_point_point(p, x[2]);
        if (d2 <= epsp2 && d2 < r[2]*r[2]) {
          // p is inside the fat vertex 2
          pts_loc[ipt] = LOC_ON_VTX2;
        }
        else if (d_edge[1] < d_edge[2]) { // Eric's rule?
          pts_loc[ipt] = LOC_INTERIOR_EDGE1;
        }
        else {
          pts_loc[ipt] = LOC_INTERIOR_EDGE2;
        }
      } // end if on edges 1 and 2
      else {
        //  p is in the interior of the fat edge 1
        pts_loc[ipt] = LOC_INTERIOR_EDGE1;
      }
    } // end if on edge 1
    else if (on_edge[2]) {
      //  p is in the interior of the fat edge 2
      pts_loc[ipt] = LOC_INTERIOR_EDGE2;
    } // end if on edge 2

    if (pts_det != NULL && pts_loc[ipt] <= LOC_INTERIOR_FACE) {
      // pts_det[ipt] = 0; // ??
    }
  }

}



static inline int
_intersect_edges_plane
(
 _point_triangle_loc_t vtx_loc[3],
 double                vtx_det[3],
 double                tri_coord[3][3],
 double                inter_coord[6],
 _point_triangle_loc_t inter_loc_edge[2]
 )
{
  int n_inter = 0;

  for (int ivtx1 = 0; ivtx1 < 3; ivtx1++) {
    int ivtx2 = (ivtx1+1)%3;

    if ((vtx_loc[ivtx1] == LOC_EXTERIOR_ABOVE  &&
         vtx_loc[ivtx2] == LOC_EXTERIOR_BELOW) ||
        (vtx_loc[ivtx1] == LOC_EXTERIOR_BELOW  &&
         vtx_loc[ivtx2] == LOC_EXTERIOR_ABOVE)) {
      // intersection point in the interior of edgeA [ivtx1, ivtx2]
      double t = vtx_det[ivtx1] / (vtx_det[ivtx1] - vtx_det[ivtx2]);
      for (int i = 0; i < 3; i++) {
        inter_coord[3*n_inter+i] = (1-t)*tri_coord[ivtx1][i] + t*tri_coord[ivtx2][i];
      }
      inter_loc_edge[n_inter] = (_point_triangle_loc_t) (ivtx1 + 3);
      n_inter++;
    }

    else if (vtx_loc[ivtx1] <= LOC_EXTERIOR_ON_PLANE) {
      // vertex on plane
      memcpy(inter_coord + 3*n_inter,
             tri_coord[ivtx1],
             sizeof(double) * 3);
      inter_loc_edge[n_inter] = (_point_triangle_loc_t) ivtx1;
      n_inter++;
    }
  }


  return n_inter;
}

static inline double
difference_of_products
(
 double a,
 double b,
 double c,
 double d
 )
{
  double cd  = c * d;
  double err = fma(-c, d,  cd);
  double dop = fma( a, b, -cd);
  return dop + err;
}

static inline void
_triangle_normal
(
 const double *a,
 const double *b,
 const double *c,
 double       *normal
 )
{
  double u[3] = {b[0] - a[0], b[1] - a[1], b[2] - a[2]};
  double v[3] = {c[0] - a[0], c[1] - a[1], c[2] - a[2]};

  normal[0] = difference_of_products(u[1], v[2], u[2], v[1]);
  normal[1] = difference_of_products(u[2], v[0], u[0], v[2]);
  normal[2] = difference_of_products(u[0], v[1], u[1], v[0]);
}








static int
_intersect_triangle_triangle_coplanar
(
 // double  coordA[3][3],
 // double  coordB[3][3],
 double                 coordA[9],
 double                 coordB[9],
 double                 epsA[3],
 double                 epsB[3],
 _point_triangle_loc_t  vtxA_locB[3],
 _point_triangle_loc_t  vtxB_locA[3],
 double                *inter_coord, // keep inter_coordA and inter_coordB instead?
 _point_triangle_loc_t *inter_locA,
 _point_triangle_loc_t *inter_locB
 )
{
  int dbg = 0;
  PDM_UNUSED(epsA);
  PDM_UNUSED(epsB);

  if (dbg) {
    log_trace("coordA :\n");
    for (int i = 0; i < 3; i++) {
      log_trace("  %f %f %f\n", coordA[3*i], coordA[3*i+1], coordA[3*i+2]);
    }
    log_trace("coordB :\n");
    for (int i = 0; i < 3; i++) {
      log_trace("  %f %f %f\n", coordB[3*i], coordB[3*i+1], coordB[3*i+2]);
    }
  }


  double normalA[3];
  _triangle_normal(&coordA[0], &coordA[3], &coordA[6], normalA);
  double normalB[3];
  _triangle_normal(&coordB[0], &coordB[3], &coordB[6], normalB);

  double sign = PDM_SIGN(PDM_DOT_PRODUCT(normalA, normalB));

  int n_inter = 0;

  double inter_tA[6];
  double inter_tB[6];

  /* Intersect edges of A with edges of B */
  for (int iedgeA = 0; iedgeA < 3; iedgeA++) {

    int ivtxA1 = iedgeA;
    int ivtxA2 = (iedgeA+1)%3;

    for (int iedgeB = 0; iedgeB < 3; iedgeB++) {

      int ivtxB1, ivtxB2;
      if (sign > 0) {
        ivtxB1 = iedgeB;
        ivtxB2 = (iedgeB+1)%3;
      }
      else {
        ivtxB1 = (iedgeB+1)%3;
        ivtxB2 = iedgeB;
      }

      if (dbg) {
        log_trace("edgeA %d, edgeB %d (%d %d), vtxA_locB %d %d, vtxB_locA %d %d\n",
                  iedgeA, iedgeB,
                  ivtxB1, ivtxB2,
                  vtxA_locB[ivtxA1], vtxA_locB[ivtxA2],
                  vtxB_locA[ivtxB1], vtxB_locA[ivtxB2]);
      }

      /* Detect special cases (vtx on vtx, vtx on edge, collinear edges, ...) */
      if (vtxA_locB[ivtxA1] == (_point_triangle_loc_t) ivtxB1) {
        assert(vtxB_locA[ivtxB1] == (_point_triangle_loc_t) ivtxA1);
        if (dbg) {
          log_trace("  A1 on B1\n");
        }
        for (int i = 0; i < 3; i++) {
          inter_coord[3*n_inter+i] = 0.5*(coordA[3*ivtxA1+i] + coordB[3*ivtxB1+i]);
        }
        if (dbg) {
          log_trace("  inter_coord = %f %f %f\n",
                    inter_coord[3*n_inter], inter_coord[3*n_inter+1], inter_coord[3*n_inter+2]);
        }
        inter_locA[n_inter] = (_point_triangle_loc_t) (ivtxA1);
        inter_locB[n_inter] = (_point_triangle_loc_t) (ivtxB1);
        inter_tA  [n_inter] = 0;
        inter_tB  [n_inter] = 0;
        n_inter++;
      }
      else if (vtxA_locB[ivtxA1] == (_point_triangle_loc_t) iedgeB+3) {
        /* vtx A interior edge B */
        if (dbg) {
          log_trace("  A1 in edgeB\n");
        }
        double dist2;
        double weight[2];
        double _coordB[6];
        memcpy(_coordB    , coordB + 3*ivtxB1, sizeof(double)*3);
        memcpy(_coordB + 3, coordB + 3*ivtxB2, sizeof(double)*3);
        PDM_line_evaluate_position(&coordA[3*ivtxA1],
                                   _coordB,
                                   NULL,
                                   &dist2,
                                   weight);
        double tB = weight[1];

        for (int i = 0; i < 3; i++) {
          inter_coord[3*n_inter+i] = 0.5*(coordA[3*ivtxA1+i] +
                                          (1-tB)*coordB[3*ivtxB1+i] + tB*coordB[3*ivtxB2+i]);
        }
        if (dbg) {
          log_trace("  inter_coord = %f %f %f\n",
                    inter_coord[3*n_inter], inter_coord[3*n_inter+1], inter_coord[3*n_inter+2]);
        }
        inter_locA[n_inter] = (_point_triangle_loc_t) (ivtxA1);
        inter_locB[n_inter] = (_point_triangle_loc_t) (iedgeB + 3);
        inter_tA  [n_inter] = 0;
        inter_tB  [n_inter] = tB; // (1 - tB) if sign < 0 ??
        n_inter++;
      }
      else if (vtxB_locA[ivtxB1] == (_point_triangle_loc_t) iedgeA+3) {
        /* vtx B interior edge A */
        if (dbg) {
          log_trace("  B1 in edgeA\n");
        }
        double dist2;
        double weight[2];
        double _coordA[6];
        memcpy(_coordA    , coordA + 3*ivtxA1, sizeof(double)*3);
        memcpy(_coordA + 3, coordA + 3*ivtxA2, sizeof(double)*3);
        PDM_line_evaluate_position(&coordB[3*ivtxB1],
                                   _coordA,
                                   NULL,
                                   &dist2,
                                   weight);
        double tA = weight[1];

        for (int i = 0; i < 3; i++) {
          inter_coord[3*n_inter+i] = 0.5*((1-tA)*coordA[3*ivtxA1+i] + tA*coordA[3*ivtxA2+i] +
                                          coordB[3*ivtxB1+i]);
        }
        if (dbg) {
          log_trace("  inter_coord = %f %f %f\n",
                    inter_coord[3*n_inter], inter_coord[3*n_inter+1], inter_coord[3*n_inter+2]);
        }
        inter_locA[n_inter] = (_point_triangle_loc_t) (iedgeA + 3);
        inter_locB[n_inter] = (_point_triangle_loc_t) (ivtxB1);
        inter_tA  [n_inter] = tA;
        inter_tB  [n_inter] = 0;
        n_inter++;
      }
      else if (vtxA_locB[ivtxA1] == (_point_triangle_loc_t) ivtxB2   ||
               vtxA_locB[ivtxA2] == (_point_triangle_loc_t) ivtxB1   ||
               vtxA_locB[ivtxA2] == (_point_triangle_loc_t) ivtxB2   ||
               vtxA_locB[ivtxA2] == (_point_triangle_loc_t) iedgeB+3 ||
               vtxB_locA[ivtxB2] == (_point_triangle_loc_t) iedgeA+3) {
        /* Will be handled further down the loop to avoid duplicate points */
        continue;
      }
      else {
        /* General case (regular, tranversal (or empty) intersection) */
        double tA;
        double tB;
        PDM_line_intersect_t stat = PDM_line_intersection_mean_square(&coordA[3*ivtxA1],
                                                                      &coordA[3*ivtxA2],
                                                                      &coordB[3*ivtxB1],
                                                                      &coordB[3*ivtxB2],
                                                                      &tA,
                                                                      &tB);
        if (stat == PDM_LINE_INTERSECT_YES) {
          if (dbg) {
            log_trace("  regular intersection\n");
          }
          for (int i = 0; i < 3; i++) {
            inter_coord[3*n_inter+i] = 0.5*((1-tA)*coordA[3*ivtxA1+i] + tA*coordA[3*ivtxA2+i] +
                                            (1-tB)*coordB[3*ivtxB1+i] + tB*coordB[3*ivtxB2+i]);
          }
          if (dbg) {
            log_trace("  inter_coord = %f %f %f\n",
                      inter_coord[3*n_inter], inter_coord[3*n_inter+1], inter_coord[3*n_inter+2]);
          }
          inter_locA[n_inter] = (_point_triangle_loc_t) (iedgeA + 3);
          inter_locB[n_inter] = (_point_triangle_loc_t) (iedgeB + 3);
          inter_tA  [n_inter] = tA;
          inter_tB  [n_inter] = tB; // (1 - tB) if sign < 0 ??
          n_inter++;
        }
      }
    }
  }

  /* Add potential interior points */
  for (int ivtx = 0; ivtx < 3; ivtx++) {
    if (vtxA_locB[ivtx] == LOC_INTERIOR_FACE) {
      memcpy(inter_coord + 3*n_inter, coordA + 3*ivtx, sizeof(double) * 3);
      inter_locA[n_inter] = (_point_triangle_loc_t) ivtx;
      inter_locB[n_inter] = LOC_INTERIOR_FACE;
      n_inter++;
    }
  }

  for (int ivtx = 0; ivtx < 3; ivtx++) {
    if (vtxB_locA[ivtx] == LOC_INTERIOR_FACE) {
      memcpy(inter_coord + 3*n_inter, coordB + 3*ivtx, sizeof(double) * 3);
      inter_locA[n_inter] = LOC_INTERIOR_FACE;
      inter_locB[n_inter] = (_point_triangle_loc_t) ivtx;
      n_inter++;
    }
  }


  if (n_inter == 0) {
    /* Empty intersection */
    return 0;
  }

  /* Poly-clipping */
  PDM_UNUSED(inter_tA);
  PDM_UNUSED(inter_tB);
  //...

  return n_inter;
}



/**
 * \brief Compute intersection between two (fat) triangles
 *
 * Returns between 0 and 6 intersection points and their classification wrt both triangles.
 * The intersection points are ordered :
 *  - if 2   : the segment is oriented by the vector normalA x normalB
 *  - if > 2 : the polygon is oriented by normalA
 *
 * \param [in ] coordA       Coordinates of the vertices of triangle A
 * \param [in ] coordB       Coordinates of the vertices of triangle B
 * \param [in ] epsA         Tolerance of the vertices of triangle A
 * \param [in ] epsB         Tolerance of the vertices of triangle B
 * \param [out] inter_coord  Coordinates of the intersection points (size: up to 6*3)
 * \param [out] inter_locA   Classification of intersection points wrt A (size: up to 6)
 * \param [out] inter_locB   Classification of intersection points wrt B (size: up to 6)
 *
 * \return Number of intersection points (<= 6)
 */

static int
_intersect_triangle_triangle
(
 // double  coordA[3][3],
 // double  coordB[3][3],
 double                 coordA[9],
 double                 coordB[9],
 double                 epsA[3],
 double                 epsB[3],
 double                *inter_coord,
 _point_triangle_loc_t *inter_locA,
 _point_triangle_loc_t *inter_locB
 )
{
  int verbose = 0;
  int dbg     = 0;

  /* Classify vtx of A wrt B */
  _point_triangle_loc_t vtxA_locB[3];
  double vtxA_detB[3];

   double tri_coordB[3][3] = {
    {coordB[0], coordB[1], coordB[2]},
    {coordB[3], coordB[4], coordB[5]},
    {coordB[6], coordB[7], coordB[8]}
  };

  _classify_pts_wrt_triangle(tri_coordB,
                             epsB,
                             3,
                             coordA,
                             epsA,
                             vtxA_locB,
                             vtxA_detB,
                             verbose);
  if (dbg) {
    PDM_log_trace_array_int((int *) vtxA_locB, 3, "vtxA_locB : ");
  }

  if ((vtxA_locB[0] == LOC_EXTERIOR_ABOVE  &&
       vtxA_locB[1] == LOC_EXTERIOR_ABOVE  &&
       vtxA_locB[2] == LOC_EXTERIOR_ABOVE) ||
      (vtxA_locB[0] == LOC_EXTERIOR_BELOW  &&
       vtxA_locB[1] == LOC_EXTERIOR_BELOW  &&
       vtxA_locB[2] == LOC_EXTERIOR_BELOW)) {
    // Clearly empty intersection
    return 0;
  }

  else if (vtxA_locB[0] <= LOC_INTERIOR_FACE &&
           vtxA_locB[1] <= LOC_INTERIOR_FACE &&
           vtxA_locB[2] <= LOC_INTERIOR_FACE) {
    // A is inside B
    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++) {
        inter_coord[3*i+j] = coordA[3*i+j];
      }
      inter_locA[i] = (_point_triangle_loc_t) i;
      inter_locB[i] = vtxA_locB[i];
    }
    return 3;
  }


  /* Classify vtx of B wrt A */
  _point_triangle_loc_t vtxB_locA[3];
  double vtxB_detA[3];

  double tri_coordA[3][3] = {
    {coordA[0], coordA[1], coordA[2]},
    {coordA[3], coordA[4], coordA[5]},
    {coordA[6], coordA[7], coordA[8]}
  };

  _classify_pts_wrt_triangle(tri_coordA,
                             epsA,
                             3,
                             coordB,
                             epsB,
                             vtxB_locA,
                             vtxB_detA,
                             verbose);

  if (dbg) {
    PDM_log_trace_array_int((int *) vtxB_locA, 3, "vtxB_locA : ");
  }

  if ((vtxB_locA[0] == LOC_EXTERIOR_ABOVE  &&
       vtxB_locA[1] == LOC_EXTERIOR_ABOVE  &&
       vtxB_locA[2] == LOC_EXTERIOR_ABOVE) ||
      (vtxB_locA[0] == LOC_EXTERIOR_BELOW  &&
       vtxB_locA[1] == LOC_EXTERIOR_BELOW  &&
       vtxB_locA[2] == LOC_EXTERIOR_BELOW)) {
    // Clearly empty intersection
    return 0;
  }

  else if (vtxB_locA[0] <= LOC_INTERIOR_FACE &&
           vtxB_locA[1] <= LOC_INTERIOR_FACE &&
           vtxB_locA[2] <= LOC_INTERIOR_FACE) {
    // B is inside A
    double normalA[3];
    _triangle_normal(&coordA[0], &coordA[3], &coordA[6], normalA);
    double normalB[3];
    _triangle_normal(&coordB[0], &coordB[3], &coordB[6], normalB);

    double sign = PDM_SIGN(PDM_DOT_PRODUCT(normalA, normalB));

    for (int i = 0; i < 3; i++) {
      int vtx_id = i;
      if (sign < 0) {
        vtx_id = 2-i;
      }
      for (int j = 0; j < 3; j++) {
        inter_coord[3*i+j] = coordB[3*vtx_id+j];
      }
      inter_locA[i] = vtxB_locA[vtx_id];
      inter_locB[i] = (_point_triangle_loc_t) vtx_id;
    }
    return 3;
  }

  /* Handle coplanar case */
  // /!\ symmetry of "on plane" is broken...
  // if (vtxB_locA[0] <= LOC_EXTERIOR_ON_PLANE &&
  //     vtxB_locA[1] <= LOC_EXTERIOR_ON_PLANE &&
  //     vtxB_locA[2] <= LOC_EXTERIOR_ON_PLANE) {

  //   assert(vtxA_locB[0] <= LOC_EXTERIOR_ON_PLANE &&
  //          vtxA_locB[1] <= LOC_EXTERIOR_ON_PLANE &&
  //          vtxA_locB[2] <= LOC_EXTERIOR_ON_PLANE);
  // if (vtxB_locA[0] <= LOC_EXTERIOR_ON_PLANE &&
  //     vtxB_locA[1] <= LOC_EXTERIOR_ON_PLANE &&
  //     vtxB_locA[2] <= LOC_EXTERIOR_ON_PLANE &&
  //     vtxA_locB[0] <= LOC_EXTERIOR_ON_PLANE &&
  //     vtxA_locB[1] <= LOC_EXTERIOR_ON_PLANE &&
  //     vtxA_locB[2] <= LOC_EXTERIOR_ON_PLANE) {
  int coplanar = 0;
  if (vtxA_locB[0] <= LOC_EXTERIOR_ON_PLANE &&
      vtxA_locB[1] <= LOC_EXTERIOR_ON_PLANE &&
      vtxA_locB[2] <= LOC_EXTERIOR_ON_PLANE) {
    coplanar = 1;
    // is this legit??
    for (int i = 0; i < 3; i++) {
      vtxB_locA[i] = PDM_MIN(vtxB_locA[i], LOC_EXTERIOR_ON_PLANE);
    }
  }
  else if (vtxB_locA[0] <= LOC_EXTERIOR_ON_PLANE &&
           vtxB_locA[1] <= LOC_EXTERIOR_ON_PLANE &&
           vtxB_locA[2] <= LOC_EXTERIOR_ON_PLANE) {
    coplanar = 1;
    // is this legit??
    for (int i = 0; i < 3; i++) {
      vtxA_locB[i] = PDM_MIN(vtxA_locB[i], LOC_EXTERIOR_ON_PLANE);
    }
  }

  if (coplanar) {
    // coplanar triangles
    return _intersect_triangle_triangle_coplanar(coordA,
                                                 coordB,
                                                 epsA,
                                                 epsB,
                                                 vtxA_locB,
                                                 vtxB_locA,
                                                 inter_coord,
                                                 inter_locA,
                                                 inter_locB);
  }



  /* Intersect edges of A with plane B */
  double                edgeA_inter_planeB_coord[6];
  _point_triangle_loc_t edgeA_inter_planeB_locA [2];
  int n_edgeA_inter_planeB = _intersect_edges_plane(vtxA_locB,
                                                    vtxA_detB,
                                                    tri_coordA,
                                                    edgeA_inter_planeB_coord,
                                                    edgeA_inter_planeB_locA);


  if (dbg) {
    log_trace("n_edgeA_inter_planeB = %d\n", n_edgeA_inter_planeB);
    for (int i = 0; i < n_edgeA_inter_planeB; i++) {
      log_trace("  %d: %f %f %f (loc : %d)\n",
                i,
                edgeA_inter_planeB_coord[3*i  ],
                edgeA_inter_planeB_coord[3*i+1],
                edgeA_inter_planeB_coord[3*i+2],
                edgeA_inter_planeB_locA[i]);
    }

    if (n_edgeA_inter_planeB == 1) {
      PDM_vtk_write_point_cloud("edgeA_inter_planeB.vtk",
                                1,
                                edgeA_inter_planeB_coord,
                                NULL,
                                NULL);
    }
    else if (n_edgeA_inter_planeB == 2) {
      PDM_vtk_write_lines("edgeA_inter_planeB.vtk",
                          1,
                          edgeA_inter_planeB_coord,
                          NULL,
                          NULL);
    }
  }

  /* Intersect edges of B with plane A */
  double                edgeB_inter_planeA_coord[6];
  _point_triangle_loc_t edgeB_inter_planeA_locB [2];
  int n_edgeB_inter_planeA = _intersect_edges_plane(vtxB_locA,
                                                    vtxB_detA,
                                                    tri_coordB,
                                                    edgeB_inter_planeA_coord,
                                                    edgeB_inter_planeA_locB);

  if (dbg) {
    log_trace("n_edgeB_inter_planeA = %d\n", n_edgeB_inter_planeA);
    for (int i = 0; i < n_edgeB_inter_planeA; i++) {
      log_trace("  %d: %f %f %f (loc : %d)\n",
                i,
                edgeB_inter_planeA_coord[3*i  ],
                edgeB_inter_planeA_coord[3*i+1],
                edgeB_inter_planeA_coord[3*i+2],
                edgeB_inter_planeA_locB[i]);
    }

    if (n_edgeB_inter_planeA == 1) {
      PDM_vtk_write_point_cloud("edgeB_inter_planeA.vtk",
                                1,
                                edgeB_inter_planeA_coord,
                                NULL,
                                NULL);
    }
    else if (n_edgeB_inter_planeA == 2) {
      PDM_vtk_write_lines("edgeB_inter_planeA.vtk",
                          1,
                          edgeB_inter_planeA_coord,
                          NULL,
                          NULL);
    }
    // else if (n_edgeB_inter_planeA == 3) {
    //   PDM_vtk_write_polydata("edgeB_inter_planeA.vtk",
    //                          1,
    //                          edgeB_inter_planeA_coord,
    //                          NULL,
    //                          NULL);
    // }
  }

  /* Check for duplicate points between edgeA_inter_planeB_coord and edgeA_inter_planeB_coord */
  //....


  // quick and dirty
  if (1) {
    double normalA[3];
    _triangle_normal(tri_coordA[0], tri_coordA[1], tri_coordA[2], normalA);
    double normalB[3];
    _triangle_normal(tri_coordB[0], tri_coordB[1], tri_coordB[2], normalB);

    double direction[3];
    PDM_CROSS_PRODUCT(direction, normalA, normalB);

    // place line origin at iso-barycenter of intersection points
    // to center u-values on 0 and get maximal precision
    double origin[3] = {0., 0., 0.};
    for (int i = 0; i < n_edgeA_inter_planeB; i++) {
      for (int j = 0; j < 3; j++) {
        origin[j] += edgeA_inter_planeB_coord[3*i+j];
      }
    }
    for (int i = 0; i < n_edgeB_inter_planeA; i++) {
      for (int j = 0; j < 3; j++) {
        origin[j] += edgeB_inter_planeA_coord[3*i+j];
      }
    }
    double normalization = 1./(double) (n_edgeA_inter_planeB + n_edgeB_inter_planeA);
    for (int j = 0; j < 3; j++) {
      origin[j] *= normalization;
    }

    double u[4];
    double deviation[4];
    double idd = 1./PDM_DOT_PRODUCT(direction, direction);
    int    order[4];
    double u_minA =  HUGE_VAL;
    double u_maxA = -HUGE_VAL;
    for (int i = 0; i < n_edgeA_inter_planeB; i++) {
      // u    [i] = PDM_DOT_PRODUCT(direction, edgeA_inter_planeB_coord + 3*i);
      u[i] = 0.;
      deviation[i] = 0.;
      for (int j = 0; j < 3; j++) {
        double x = edgeA_inter_planeB_coord[3*i+j] - origin[j];
        u[i] += direction[j] * x;
        double delta = x - u[i]*direction[j]*idd;
        deviation[i] += delta*delta;
      }
      deviation[i] = sqrt(deviation[i]);
      order[i] = i;
      u_minA = PDM_MIN(u_minA, u[i]);
      u_maxA = PDM_MAX(u_maxA, u[i]);
    }
    double u_minB =  HUGE_VAL;
    double u_maxB = -HUGE_VAL;
    for (int i = 0; i < n_edgeB_inter_planeA; i++) {
      // u    [n_edgeA_inter_planeB+i] = PDM_DOT_PRODUCT(direction, edgeB_inter_planeA_coord + 3*i);
      u[n_edgeA_inter_planeB+i] = 0.;
      deviation[n_edgeA_inter_planeB+i] = 0.;
      for (int j = 0; j < 3; j++) {
        double x = edgeB_inter_planeA_coord[3*i+j] - origin[j];
        u[n_edgeA_inter_planeB+i] += direction[j] * x;
        double delta = x - u[n_edgeA_inter_planeB+i]*direction[j]*idd;
        deviation[n_edgeA_inter_planeB+i] += delta*delta;
      }
      deviation[n_edgeA_inter_planeB+i] = sqrt(deviation[n_edgeA_inter_planeB+i]);
      order[n_edgeA_inter_planeB+i] = i + n_edgeA_inter_planeB;
      u_minB = PDM_MIN(u_minB, u[n_edgeA_inter_planeB+i]);
      u_maxB = PDM_MAX(u_maxB, u[n_edgeA_inter_planeB+i]);
    }

    if (u_minA > u_maxB || u_minB > u_maxA) {
      return 0;
    }

    PDM_sort_double(u, order, n_edgeA_inter_planeB + n_edgeB_inter_planeA);

    for (int i = 0; i < 2; i++) {
      if (order[i+1] < n_edgeA_inter_planeB) {
        memcpy(inter_coord    + 3*i,
               edgeA_inter_planeB_coord + 3*order[i+1],
               sizeof(double) * 3);
        inter_locA[i] = edgeA_inter_planeB_locA[i];
        inter_locB[i] = -1; // ?
      }
      else {
        memcpy(inter_coord    + 3*i,
               edgeB_inter_planeA_coord + 3*(order[i+1] - n_edgeA_inter_planeB),
               sizeof(double) * 3);
        inter_locA[i] = -1; // ?
        inter_locB[i] = edgeB_inter_planeA_locB[i];
      }
    }

    if (dbg) {
      log_trace("deviation = %e %e %e %e\n",
                deviation[0],
                deviation[1],
                deviation[2],
                deviation[3]);
      PDM_log_trace_array_int((int *) inter_locA, 2, "inter_locA : ");
      PDM_log_trace_array_int((int *) inter_locB, 2, "inter_locB : ");
      PDM_vtk_write_lines("triA_inter_triB.vtk",
                          1,
                          inter_coord,
                          NULL,
                          NULL);

      double line[6];
      log_trace("%e %e\n", PDM_MIN(u_minA, u_minB), PDM_MAX(u_maxA, u_maxB));
      for (int i = 0; i < 3; i++) {
        line[  i] = origin[i] + PDM_MIN(u_minA, u_minB)*direction[i]*idd;
        line[3+i] = origin[i] + PDM_MAX(u_maxA, u_maxB)*direction[i]*idd;
      }
      PDM_vtk_write_lines("line.vtk",
                          1,
                          line,
                          NULL,
                          NULL);

    }
  }

  return 2; // !!!
}





static void
_dump_triangle
(
 const char   *filename,
       double  coord[3][3]
 )
{
  int connec[3] = {1, 2, 3};
  double vtx_coord[9] = {
    coord[0][0], coord[0][1], coord[0][2],
    coord[1][0], coord[1][1], coord[1][2],
    coord[2][0], coord[2][1], coord[2][2]
  };

  PDM_g_num_t vtx_id[3] = {1, 2, 3};

  PDM_vtk_write_std_elements(filename,
                             3,
                             vtx_coord,
                             vtx_id,
                             PDM_MESH_NODAL_TRIA3,
                             1,
                             connec,
                             NULL,
                             0, NULL, NULL);
}




/**
 *
 * \brief  Usage
 *
 */

static void
_usage
(
int exit_code
)
{
  PDM_printf
    ("\n"
     "  Usage: \n\n"
     "  -n      <level>  Number of points (default : 10).\n\n"
     "  -radius <level>  Radius of domain (default : 10).\n\n"
     "  -local           Number of points is local (default : global).\n\n"
     "  -post            Export vtk (default : false).\n\n"
     "  -h               This message.\n\n");

  exit (exit_code);
}



/**
 *
 * \brief  Read arguments from the command line
 *
 */

static void
_read_args
(
 int            argc,
 char         **argv,
 int           *seed,
 double        *eps,
 int           *n_repeat,
 int           *randA,
 double        *a0,
 double        *a1,
 double        *a2,
 double        *scale,
 double        *noise_scale,
 int           *verbose
 )
{
  int i = 1;

  /* Parse and check command line */

  while (i < argc) {

    if (strcmp(argv[i], "-h") == 0)
      _usage(EXIT_SUCCESS);

    else if (strcmp(argv[i], "-eps") == 0) {
      i++;
      if (i >= argc) {
        _usage(EXIT_FAILURE);
      }
      else {
        *eps = atof(argv[i]);
      }
    }
    else if (strcmp(argv[i], "-seed") == 0) {
      i++;
      if (i >= argc) {
        _usage(EXIT_FAILURE);
      }
      else {
        *seed = atoi(argv[i]);
      }
    }
    else if (strcmp(argv[i], "-rep") == 0) {
      i++;
      if (i >= argc) {
        _usage(EXIT_FAILURE);
      }
      else {
        *n_repeat = atoi(argv[i]);
      }
    }
    else if (strcmp(argv[i], "-randA") == 0) {
      *randA = 1;
    }
    else if (strcmp(argv[i], "-a0") == 0) {
      for (int j = 0; j < 3; j++) {
        i++;
        if (i >= argc) {
          _usage(EXIT_FAILURE);
        }
        else {
          a0[j] = atof(argv[i]);
        }
      }
    }
    else if (strcmp(argv[i], "-a1") == 0) {
      for (int j = 0; j < 3; j++) {
        i++;
        if (i >= argc) {
          _usage(EXIT_FAILURE);
        }
        else {
          a1[j] = atof(argv[i]);
        }
      }
    }
    else if (strcmp(argv[i], "-a2") == 0) {
      for (int j = 0; j < 3; j++) {
        i++;
        if (i >= argc) {
          _usage(EXIT_FAILURE);
        }
        else {
          a2[j] = atof(argv[i]);
        }
      }
    }
    else if (strcmp(argv[i], "-scale") == 0) {
      i++;
      if (i >= argc) {
        _usage(EXIT_FAILURE);
      }
      else {
        *scale = atof(argv[i]);
      }
    }
    else if (strcmp(argv[i], "-noise_scale") == 0) {
      i++;
      if (i >= argc) {
        _usage(EXIT_FAILURE);
      }
      else {
        *noise_scale = atof(argv[i]);
      }
    }
    else if (strcmp(argv[i], "-verbose") == 0) {
      *verbose = 1;
    }

    else {
      _usage(EXIT_FAILURE);
    }
    i++;
  }
}



/**
 *
 * \brief  Main
 *
 */

int main(int argc, char *argv[])
{
  /*
   *  Init
   */
  PDM_MPI_Init(&argc, &argv);

  PDM_predicate_exactinit();

  int    seed        = 0;
  double eps         = 1e-6;
  int    n_repeat    = 1;
  int    randA       = 0;
  double a0[3]       = {0, 0, 0};
  double a1[3]       = {1, 0, 0};
  double a2[3]       = {0, 1, 0};
  double scale       = 1.;
  double noise_scale = 0.;
  int    verbose     = 0;
  _read_args(argc,
             argv,
             &seed,
             &eps,
             &n_repeat,
             &randA,
             a0,
             a1,
             a2,
             &scale,
             &noise_scale,
             &verbose);

  if (seed < 0) {
    seed = time(NULL);
  }
  srand(seed);


  double vtx_coordA[3][3] = {
    {a0[0], a0[1], a0[2]},
    {a1[0], a1[1], a1[2]},
    {a2[0], a2[1], a2[2]}
  };

  if (randA) {
    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 2; j++) {
        vtx_coordA[i][j] = scale * (2*_rand01() - 1.);
      }
    }
  }

  double eps_abs = 1e-15;
  double vtx_epsA[3] = {HUGE_VAL, HUGE_VAL, HUGE_VAL};
  for (int i = 0; i < 3; i++) {
    int j = (i+1)%3;
    double d = _dist2_point_point(vtx_coordA[i], vtx_coordA[j]);
    vtx_epsA[i] = PDM_MIN(vtx_epsA[i], d);
    vtx_epsA[j] = PDM_MIN(vtx_epsA[j], d);
  }

  for (int i = 0; i < 3; i++) {
    vtx_epsA[i] = PDM_MAX(eps_abs, eps*sqrt(vtx_epsA[i]));
  }




  // double vtx_coordB[3][3] = {
  //   { 0.5,  0.5, 0.},
  //   {-0.5,  0.5, 0.},
  //   { 0.5, -0.5, 0.}
  // };

  // if (randA) {
  //   for (int i = 0; i < 3; i++) {
  //     for (int j = 0; j < 2; j++) {
  //       // vtx_coordB[i][j] = 1./3. + 0.3*(vtx_coordB[i][j] - 1./3.);
  //       vtx_coordB[i][j] += scale * (2*_rand01() - 1.);
  //     }
  //   }
  // }

  // double vtx_epsB[3] = {HUGE_VAL, HUGE_VAL, HUGE_VAL};
  // for (int i = 0; i < 3; i++) {
  //   int j = (i+1)%3;
  //   double d = _dist2_point_point(vtx_coordB[i], vtx_coordB[j]);
  //   vtx_epsB[i] = PDM_MIN(vtx_epsB[i], d);
  //   vtx_epsB[j] = PDM_MIN(vtx_epsB[j], d);
  // }

  // for (int i = 0; i < 3; i++) {
  //   vtx_epsB[i] = PDM_MAX(eps_abs, eps*sqrt(vtx_epsB[i]));
  // }

  double vtx_coordB[3][3];
  double vtx_epsB[3] = {scale*eps, scale*eps, scale*eps};


  _point_triangle_loc_t vtxB_locA[3];

  for (int i = 0; i < 3; i++) {
    while (1) {
      int keep = 1;
      vtxB_locA[i] = (_point_triangle_loc_t) (8*sqrt(_rand01()));

      for (int j = 0; j < i; j++) {
        if (vtxB_locA[j] == vtxB_locA[i] && vtxB_locA[i] < 3) {
          keep = 0;
          break;
        }
      }

      if (keep) {
        break;
      }
    }
  }

  // double min_epsB = HUGE_VAL;
  // for (int i = 0; i < 3; i++) {
  //   min_epsB = PDM_MIN(min_epsB, vtx_epsB[i]);
  // }

  // double noise_scale = min_epsB; // abs(noise_scale) <= 1
  // if (eps >= 0) {
  //   noise_scale /= eps;
  // }
  _gen_config(vtx_coordA,
              vtx_epsA,
              vtxB_locA,
              vtx_epsB,
              noise_scale,
              vtx_coordB);

  if (1) {
    for (int i = 0; i < 3; i++) {
      vtx_epsB[i] = HUGE_VAL;
    }

    for (int i = 0; i < 3; i++) {
      int j = (i+1)%3;
      double d = _dist2_point_point(vtx_coordB[i], vtx_coordB[j]);
      vtx_epsB[i] = PDM_MIN(vtx_epsB[i], d);
      vtx_epsB[j] = PDM_MIN(vtx_epsB[j], d);
    }

    for (int i = 0; i < 3; i++) {
      vtx_epsB[i] = PDM_MAX(eps_abs, eps*sqrt(vtx_epsB[i]));
    }
  }

  if (verbose) {
    log_trace("seed = %d, noise_scale = %f, scale = %f\n",
              seed, noise_scale, scale);
    log_trace("vtx_epsA : %e %e %e\n",
              vtx_epsA[0], vtx_epsA[1], vtx_epsA[2]);
    log_trace("vtx_epsB : %e %e %e\n",
              vtx_epsB[0], vtx_epsB[1], vtx_epsB[2]);
  }

  double coordA[9] = {
    vtx_coordA[0][0], vtx_coordA[0][1], vtx_coordA[0][2],
    vtx_coordA[1][0], vtx_coordA[1][1], vtx_coordA[1][2],
    vtx_coordA[2][0], vtx_coordA[2][1], vtx_coordA[2][2]
  };
  double coordB[9] = {
    vtx_coordB[0][0], vtx_coordB[0][1], vtx_coordB[0][2],
    vtx_coordB[1][0], vtx_coordB[1][1], vtx_coordB[1][2],
    vtx_coordB[2][0], vtx_coordB[2][1], vtx_coordB[2][2]
  };

  if (verbose) {
    _dump_triangle("triangleA.vtk", vtx_coordA);
    _dump_triangle("triangleB.vtk", vtx_coordB);
  }


  double                inter_coord[6*3];
  _point_triangle_loc_t inter_locA[6];
  _point_triangle_loc_t inter_locB[6];
  int n_inter = _intersect_triangle_triangle(coordA,
                                             coordB,
                                             vtx_epsA,
                                             vtx_epsB,
                                             inter_coord,
                                             inter_locA,
                                             inter_locB);

  if (verbose) {
    log_trace("n_inter = %d\n", n_inter);
    for (int i = 0; i < n_inter; i++) {
      log_trace("  %f %f %f A%d B%d\n",
                inter_coord[3*i], inter_coord[3*i+1], inter_coord[3*i+2],
                (int) inter_locA[i],
                (int) inter_locB[i]);
    }


    // PDM_vtk_write_point_cloud("interAB.vtk",
    //                           n_inter,
    //                           inter_coord,
    //                           NULL, NULL);
    const char *field_name [2] = {"locA", "locB"};
    const int  *field_value[2] = {(int *) inter_locA, (int *) inter_locB};

    PDM_vtk_write_std_elements("interAB.vtk",
                               n_inter,
                               inter_coord,
                               NULL,
                               PDM_MESH_NODAL_POINT,
                               n_inter,
                               NULL,
                               NULL,
                               2,
                               field_name,
                               field_value);
  }



  PDM_MPI_Finalize();

  return 0;
}


#if 0
/**
 *
 * \brief  Main
 *
 */

int main(int argc, char *argv[])
{
  /*
   *  Init
   */
  PDM_MPI_Init(&argc, &argv);

  int    seed     = 0;
  double eps      = 1e-6;
  int    n_repeat = 1;
  int    randA    = 0;
  double a0[3]    = {0, 0, 0};
  double a1[3]    = {1, 0, 0};
  double a2[3]    = {0, 1, 0};
  double scale    = 1.;
  _read_args(argc,
             argv,
             &seed,
             &eps,
             &n_repeat,
             &randA,
             a0,
             a1,
             a2,
             &scale);


  log_trace("seed = %d, eps = %e\n", seed, eps);


  double vtx_coordA[3][3] = {
    {a0[0], a0[1], a0[2]},
    {a1[0], a1[1], a1[2]},
    {a2[0], a2[1], a2[2]}
  };

  if (randA) {
    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++) {
        vtx_coordA[i][j] = scale * (2*_rand01() - 1.);
      }
    }
  }

  double eps_abs = 1e-15;
  double vtx_epsA[3] = {HUGE_VAL, HUGE_VAL, HUGE_VAL};
  for (int i = 0; i < 3; i++) {
    int j = (i+1)%3;
    double d = _dist2_point_point(vtx_coordA[i], vtx_coordA[j]);
    vtx_epsA[i] = PDM_MIN(vtx_epsA[i], d);
    vtx_epsA[j] = PDM_MIN(vtx_epsA[j], d);
  }

  for (int i = 0; i < 3; i++) {
    vtx_epsA[i] = PDM_MAX(eps_abs, eps*sqrt(vtx_epsA[i]));
  }


  double vtx_coordB[3][3];
  double vtx_epsB[3] = {scale*eps, scale*eps, scale*eps};

  log_trace("vtx_epsA : %e %e %e\n",
            vtx_epsA[0], vtx_epsA[1], vtx_epsA[2]);
  log_trace("vtx_epsB : %e %e %e\n",
            vtx_epsB[0], vtx_epsB[1], vtx_epsB[2]);


  _point_triangle_loc_t vtxB_locA[3];

  double noise_scale = 1.; // abs(noise_scale) <= 1

  int count    = 0;
  int count_ok = 0;

  for (int irepeat = 0; irepeat < n_repeat; irepeat++) {

    if (seed < 0) {
      srand(time(NULL));
    }
    else {
      srand(seed + irepeat);
    }

    for (int loc0 = 0; loc0 < 10; loc0++) {
      vtxB_locA[0] = (_point_triangle_loc_t) loc0;

      for (int loc1 = 0; loc1 < 10; loc1++) {
        vtxB_locA[1] = (_point_triangle_loc_t) loc1;

        for (int loc2 = 0; loc2 < 10; loc2++) {
          vtxB_locA[2] = (_point_triangle_loc_t) loc2;

          _gen_config(vtx_coordA,
                      vtx_epsA,
                      vtxB_locA,
                      vtx_epsB,
                      noise_scale,
                      vtx_coordB);



          double coordB[9] = {
            vtx_coordB[0][0], vtx_coordB[0][1], vtx_coordB[0][2],
            vtx_coordB[1][0], vtx_coordB[1][1], vtx_coordB[1][2],
            vtx_coordB[2][0], vtx_coordB[2][1], vtx_coordB[2][2]
          };
          _point_triangle_loc_t locB[3];
          double vtx_detB[3];
          _classify_pts_wrt_triangle(vtx_coordA,
                                     vtx_epsA,
                                     3,
                                     coordB,
                                     vtx_epsB,
                                     locB,
                                     vtx_detB,
                                     0);

          if (vtxB_locA[0] != locB[0] ||
              vtxB_locA[1] != locB[1] ||
              vtxB_locA[2] != locB[2]) {
            log_trace("error !!!\n");
            PDM_log_trace_array_int((int *) vtxB_locA, 3, "  vtxB_locA : ");
            PDM_log_trace_array_int((int *) locB,      3, "  locB      : ");
            log_trace("\n\n");

            char filename[999];
            sprintf(filename, "triangleB_%d%d%d_%d.vtk",
                    vtxB_locA[0], vtxB_locA[1], vtxB_locA[2],
                    irepeat);
            _dump_triangle(filename, vtx_coordB);

            sprintf(filename, "triangleA_%d.vtk",
                    irepeat);
            _dump_triangle(filename, vtx_coordA);
          }
          else {
            count_ok++;
          }
          count++;

          // if (loc0 == LOC_EXTERIOR_ON_PLANE &&
          //     loc1 == LOC_EXTERIOR_ON_PLANE &&
          //     loc2 == LOC_EXTERIOR_BELOW) {
          if (0) {
            _dump_triangle("triangleA.vtk", vtx_coordA);
            _dump_triangle("test.vtk", vtx_coordB);

            double coordA[9] = {
              vtx_coordA[0][0], vtx_coordA[0][1], vtx_coordA[0][2],
              vtx_coordA[1][0], vtx_coordA[1][1], vtx_coordA[1][2],
              vtx_coordA[2][0], vtx_coordA[2][1], vtx_coordA[2][2]
            };

            double                inter_coord[6*3];
            _point_triangle_loc_t inter_locA[6];
            _point_triangle_loc_t inter_locB[6];
            int n_inter = _intersect_triangle_triangle(coordA,
                                                       coordB,
                                                       vtx_epsA,
                                                       vtx_epsB,
                                                       inter_coord,
                                                       inter_locA,
                                                       inter_locB);
            PDM_UNUSED(n_inter);
          }

        }
      }
    }

  }



  printf("%d OK / %d\n", count_ok, count);

  if (count_ok < count) {
    _dump_triangle("triangleA.vtk", vtx_coordA);
  }


  PDM_MPI_Finalize();

  return 0;
}
#endif
