#include <math.h>
#include <sys/time.h>
#include <time.h>
#include <sys/resource.h>
#include <unistd.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>

#include "pdm_logging.h"
#include "pdm_vtk.h"
#include "pdm.h"
#include "pdm_error.h"
#include "pdm_priv.h"

#include "pdm_line.h"
#include "pdm_mesh_intersection_surf_surf_atomic.h"

/*============================================================================
 * Global variables
 *============================================================================*/

// static int dbg_enabled = 0;

/*============================================================================
 * Type definitions
 *============================================================================*/

/*=============================================================================
 * Static function definitions
 *============================================================================*/

PDM_GCC_SUPPRESS_WARNING_WITH_PUSH("-Wfloat-equal")
static int
_solve_quadratic
(
 const double a,
 const double b,
 const double c,
       double solutions[2]
 )
{
  if (a == 0) {
    // Linear equation
    if (b == 0) {
      if (c == 0) {
        // trivial equation 0 = 0
        return -1;
      }
      else {
        // no solution
        return 0;
      }
    }
    else {
      solutions[0] = -c/b;
      return 1;
    }
  }

  else {
    // True quadratic equation
    double ia = 1./a;
    double _b = b*ia;
    double _c = c*ia;

    if (_c == 0) {
      solutions[0] = PDM_MIN(0, -_b);
      solutions[1] = PDM_MAX(0, -_b);
      return 2;
    }
    else {
      double d = _b*_b - 4.*_c;
      if (d < 0) {
        // no real solution
        return 0;
      }
      else {
        // two real solutions (possibly one double solution)
        d = sqrt(d);
        double x1, x2;
        if (_b < 0) {
          x1 = 0.5*(-_b + d);
        }
        else {
          x1 = 0.5*(-_b - d);
        }
        x2 = _c/x1;
        solutions[0] = PDM_MIN(x1, x2);
        solutions[1] = PDM_MAX(x1, x2);
        return 2;
      }
    }
  }
}
PDM_GCC_SUPPRESS_WARNING_POP

static PDM_line_intersect_t
_line_intersection_projection
(
 const double  a0[3],
 const double  a1[3],
 const double  b0[3],
 const double  b1[3],
 const double  d0[3],
 const double  d1[3],
       double *ta,
       double *tb
)
{
  double a0a1[3] = {a1[0] - a0[0],
                    a1[1] - a0[1],
                    a1[2] - a0[2]};

  double b0b1[3] = {b1[0] - b0[0],
                    b1[1] - b0[1],
                    b1[2] - b0[2]};

  double d0d1[3] = {d1[0] - d0[0],
                    d1[1] - d0[1],
                    d1[2] - d0[2]};

  double a0b0[3] = {b0[0] - a0[0],
                    b0[1] - a0[1],
                    b0[2] - a0[2]};

  double a0a1_x_d0d1[3];
  PDM_CROSS_PRODUCT(a0a1_x_d0d1, a0a1, d0d1);

  double a0a1_x_d0[3];
  PDM_CROSS_PRODUCT(a0a1_x_d0, a0a1, d0);

  double c2 = PDM_DOT_PRODUCT(a0a1_x_d0d1, b0b1);
  double c1 = PDM_DOT_PRODUCT(a0a1_x_d0, b0b1) + PDM_DOT_PRODUCT(a0a1_x_d0d1, a0b0);
  double c0 = PDM_DOT_PRODUCT(a0a1_x_d0, a0b0);


  double solutions[2];
  int n_solutions = _solve_quadratic(c2, c1, c0, solutions);


  if (n_solutions < 0) {
    return PDM_LINE_INTERSECT_ON_LINE;
  }
  else if (n_solutions == 0) {
    return PDM_LINE_INTERSECT_NO;
  }
  else if (n_solutions == 1) {
    *tb = solutions[0];
    if (*tb < 0 || *tb > 1) {
      return PDM_LINE_INTERSECT_NO;
    }
  }
  else {
    if (solutions[0] >= 0 && solutions[0] <= 1) {
      if (solutions[1] >= 0 && solutions[0] <= 1) {
        // two valid solutions, which one do we choose??
        return PDM_LINE_INTERSECT_UNDEF;
      }
      else {
        *tb = solutions[0];
      }
    }
    else if (solutions[1] >= 0 && solutions[0] <= 1) {
      *tb = solutions[1];
    }
    else {
      return PDM_LINE_INTERSECT_NO;
    }
  }



  double d[3];
  for (int i = 0; i < 3; i++) {
    d[i] = (1-(*tb))*d0[i] + (*tb)*d1[i];
  }

  double l[3];
  PDM_CROSS_PRODUCT(l, b0b1, d);

  double denom = PDM_DOT_PRODUCT(l, a0a1);

  PDM_GCC_SUPPRESS_WARNING_WITH_PUSH("-Wfloat-equal")
  if (denom == 0) {
    return PDM_LINE_INTERSECT_UNDEF;
  }
  else {
    *ta = PDM_DOT_PRODUCT(l, a0b0) / denom;
    if (*ta >= 0 && *ta <= 1) {
      return PDM_LINE_INTERSECT_YES;
    }
    else {
      return PDM_LINE_INTERSECT_NO;
    }
  }
  PDM_GCC_SUPPRESS_WARNING_POP

}







static double
_line_signed_distance
(
 const double coord[2],
 const int    iline
 )
{
  switch (iline) {
  case 0: // x = 0
    return coord[0];
  case 1: // x = 1
    return 1. - coord[0];
  case 2: // y = 0
    return coord[1];
  case 3: // x + y = 1
    return 1. - coord[0] - coord[1];
  default:
    log_error("_line_signed_distance: wrong line number %d\n", iline);
  }
  return 0;
}


static void
_line_vtxA_id
(
 const int  iline,
       int *vtxA_id0,
       int *vtxA_id1
 )
{
  switch (iline) {
  case 0: // x = 0
    *vtxA_id0 = 2;
    *vtxA_id1 = 0;
    break;
  case 1: // x = 1
    *vtxA_id0 = -1;
    *vtxA_id1 = -1;
    break;
  case 2: // y = 0
    *vtxA_id0 = 0;
    *vtxA_id1 = 1;
    break;
  case 3: // x + y = 1
    *vtxA_id0 = 1;
    *vtxA_id1 = 2;
    break;
  default:
    log_error("_line_vtxA_id: wrong line number %d\n", iline);
  }
}



static inline void
_project_on_line3
(
       double *coord,
 const int     n
 )
{
  for (int i = 0; i < n; i++) {
    coord[2*i+1] = 1 - coord[2*i];
  }
}


PDM_GCC_SUPPRESS_WARNING_WITH_PUSH("-Wfloat-equal")
static int
_grandy2d
(
 double *coord,
 double  triaA_coord[9],
 double  edgeB_coord[6],
 double  edgeB_normal[6]
 )
{
  int dbg_enabled = 1;
  PDM_UNUSED(coord);
  PDM_UNUSED(triaA_coord);
  PDM_UNUSED(edgeB_coord);
  PDM_UNUSED(edgeB_normal);

  double uvA[6] = {
    0, 0,
    1, 0,
    0, 1
  };

  double tmin = 0;
  double tmax = 1;

  /* Check if the initial segment is inside the unit triangle */
  int inside = 1;
  for (int i = 0; i < 2; i++) {
    for (int j = 0; j < 2; j++) {
      if (coord[2*i+j] < 0 || coord[2*i+j] > 1) {
        inside = 0;
        break;
      }
    }
    if (!inside) {
      break;
    }

    if (1 - coord[2*i] - coord[2*i+1] < 0) {
        inside = 0;
        break;
      }
  }

  if (inside) {
    if (dbg_enabled) {
      log_trace("initial segment in triangle\n");
    }
    return 2;
  }


  /* Intersect/clip with all 4 lines */
  for (int iline = 0; iline < 4; iline++) {

    double f0 = _line_signed_distance(&coord[0], iline);
    double f1 = _line_signed_distance(&coord[2], iline);


    if (f0 == 0) {

      if (f1 > 0) {
        // 0 on, rest in => all in
        continue;
      }
      else {
        // 0 on, rest on or out
        if (iline == 3) {
          _project_on_line3(&coord[2], 1);
          return 2;
        }
        else {
          return 0;
        }
      }

    } // end if f0 == 0

    else if (f1 == 0) {

      if (f0 < 0) {
        // 1 on, rest out => all out
        if (iline == 3) {
          _project_on_line3(&coord[0], 1);
          return 2;
        }
        else {
          return 0;
        }
      }
      else {
        // 1 on, rest in => all in
        continue;
      }

    }

    else if (f0*f1 < 0) {
      // intersection
      int vtxA_id0 = -1;
      int vtxA_id1 = -1;
      _line_vtxA_id(iline, &vtxA_id0, &vtxA_id1);

      double tA, tB;
      double inter_coord[2];
      PDM_line_intersect_t stat = PDM_LINE_INTERSECT_UNDEF;
      if (edgeB_normal != NULL && vtxA_id0 >= 0 && vtxA_id1 >= 0) {
        stat = _line_intersection_projection(&triaA_coord[3*vtxA_id0],
                                             &triaA_coord[3*vtxA_id1],
                                             &edgeB_coord[0],
                                             &edgeB_coord[3],
                                             &edgeB_normal[0],
                                             &edgeB_normal[3],
                                             &tA,
                                             &tB);
        if (stat == PDM_LINE_INTERSECT_YES) {
          // intersection coord from tA or tB?????
          assert(tB >= tmin);
          assert(tB <= tmax);
          for (int i = 0; i < 2; i++) {
            inter_coord[i] = (1-tA)*uvA[2*vtxA_id0+i] + tA*uvA[2*vtxA_id1+i];
          }
        }
      }

      if (stat != PDM_LINE_INTERSECT_YES) {
        tB = f0 / (f0 - f1);
        for (int i = 0; i < 2; i++) {
          inter_coord[i] = (1-tB)*coord[i] + tB*coord[2+i];
        }

        tB = tmin + (tmax-tmin)*tB;
      }


      if (iline == 3) {
        memcpy(&coord[4], &coord[2],   sizeof(double)*2);
        memcpy(&coord[2], inter_coord, sizeof(double)*2);
        if (f0 < 0) {
          _project_on_line3(&coord[0], 1);
        }
        else {
          _project_on_line3(&coord[4], 1);
        }
        return 3;
      }
      else {
        if (f0 < 0) {
          tmin = tB;
          memcpy(&coord[0], inter_coord, sizeof(double)*2);
        }
        else {
          tmax = tB;
          memcpy(&coord[2], inter_coord, sizeof(double)*2);
        }
      }

    }

    else if (f0 < 0) {
      // no intersection, all out
      if (iline == 3) {
        _project_on_line3(coord, 2);
        return 2;
      }
      else {
        return 0;
      }
    }

    else {
      // no intersection, all in
      continue;
    }

  } // End of loop on lines

  return 0;
}
PDM_GCC_SUPPRESS_WARNING_POP





/*=============================================================================
 * Public function definitions
 *============================================================================*/

double PDM_mesh_intersection_surf_surf_atomic_compute
(
 double triaA_coord[9],
 double edgeB_coord[6],
 double edgeB_normal[6]
 )
{
  int dbg_enabled = 0;

  double mat[3][3];
  double rhs[3] = {0, 0, 0};

  for (int i = 0; i < 3; i++) {
    mat[i][0] = triaA_coord[3+i] - triaA_coord[i];
    mat[i][1] = triaA_coord[6+i] - triaA_coord[i];
  }

  double uv[6];
  for (int i = 0; i < 2; i++) {

PDM_GCC_SUPPRESS_WARNING_WITH_PUSH("-Wfloat-equal")
    if (edgeB_normal == NULL) {
      // projection = barycentric coordinates in A
      double m00 = 0;
      double m01 = 0;
      double m11 = 0;
      for (int j = 0; j < 3; j++) {
        m00 += mat[j][0] * mat[j][0];
        m01 += mat[j][0] * mat[j][1];
        m11 += mat[j][1] * mat[j][1];
        rhs[0] += mat[j][0] * (edgeB_coord[3*i+j] - triaA_coord[i]);
        rhs[1] += mat[j][1] * (edgeB_coord[3*i+j] - triaA_coord[i]);
      }

      double det = m00*m11 - m01*m01;
PDM_GCC_SUPPRESS_WARNING_WITH_PUSH("-Wfloat-equal")
      if (det == 0) {
        PDM_error(__FILE__, __LINE__, 0,
                  "cannot project (degenerate triangle A)\n");
      }
PDM_GCC_SUPPRESS_WARNING_POP

      double idet = 1./det;

      uv[2*i  ] = (rhs[0]*m11 - rhs[1]*m01) * idet;
      uv[2*i+1] = (rhs[1]*m00 - rhs[0]*m01) * idet;
    }

    else {

      for (int j = 0; j < 3; j++) {
        mat[j][2] = -edgeB_normal[3*i+j];
        rhs[j] = edgeB_coord[3*i+j] - triaA_coord[i];
      }

      double det = mat[0][0]*(mat[1][1]*mat[2][2] - mat[2][1]*mat[1][2])
      -            mat[1][0]*(mat[0][1]*mat[2][2] - mat[2][1]*mat[0][2])
      +            mat[2][0]*(mat[0][1]*mat[1][2] - mat[1][1]*mat[0][2]);

PDM_GCC_SUPPRESS_WARNING_WITH_PUSH("-Wfloat-equal")
      if (det == 0) {
        PDM_error(__FILE__, __LINE__, 0,
                  "cannot project : mat =\n"
                  "%f %f %f\n"
                  "%f %f %f\n"
                  "%f %f %f\n",
                  mat[0][0], mat[0][1], mat[0][2],
                  mat[1][0], mat[1][1], mat[1][2],
                  mat[2][0], mat[2][1], mat[2][2]);
      }
PDM_GCC_SUPPRESS_WARNING_POP

      double idet = 1./det;

      uv[2*i    ] = idet *
      (  rhs[0]*(mat[1][1]*mat[2][2] - mat[2][1]*mat[1][2])
       - rhs[1]*(mat[0][1]*mat[2][2] - mat[2][1]*mat[0][2])
       + rhs[2]*(mat[0][1]*mat[1][2] - mat[1][1]*mat[0][2]));

      uv[2*i + 1] = idet *
      (  mat[0][0]*(rhs[1]*mat[2][2] - rhs[2]*mat[1][2])
       - mat[1][0]*(rhs[0]*mat[2][2] - rhs[2]*mat[0][2])
       + mat[2][0]*(rhs[0]*mat[1][2] - rhs[1]*mat[0][2]));
    }

  }

  /* Clipping */
  if (dbg_enabled) {
    log_trace("initial segment :\n");
    for (int i = 0; i < 2; i++) {
      log_trace("%20.16f %20.16f\n", uv[2*i], uv[2*i+1]);
    }
  }
  int n_vtx = _grandy2d(uv,
                        triaA_coord,
                        edgeB_coord,
                        edgeB_normal);

  if (dbg_enabled) {
    log_trace("clipped :\n");
    for (int i = 0; i < n_vtx; i++) {
      log_trace("%20.16f %20.16f\n", uv[2*i], uv[2*i+1]);
    }
  }

  /* Compute column area */
  double area = 0;
  for (int i = 0; i < n_vtx-1; i++) {
    area += 0.5 * (uv[2*i+1] + uv[2*i+3]) * (uv[2*i] - uv[2*i+2]);
  }

  return area;
}
