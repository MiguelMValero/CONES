/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#include <float.h>

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "pdm.h"
#include "pdm_config.h"
#include "pdm_mpi.h"
#include "pdm_printf.h"
#include "pdm_error.h"
#include "pdm_priv.h"
#include "pdm_vtk.h"
#include "pdm_logging.h"
#include "pdm_array.h"
#include "pdm_linear_programming.h"

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force back Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */

/*=============================================================================
 * Static global variables
 *============================================================================*/

static const double eps_lp  = 1.e-15;
static const double big_lp  = 1.e15;
static const int    verbose = 0;

/*=============================================================================
 * Private function definitions
 *============================================================================*/

/* To avoid zero division */
static int
_is_zero
(
const double x
)
{
 return (PDM_ABS(x) < eps_lp);
}

/**
 * Solve the 1-dimensional linear programming problem:
 *  maximize c*x
 *  subject to constraints ai*x <= bi
 */
static PDM_lp_status_t
_lp_solve_1d
(
 const int  n,
 double    *a,
 double    *b,
 double     c,
 double    *x
 )
{
  double L = -big_lp;
  double R =  big_lp;

  if (verbose) PDM_log_trace_array_double(a, n, "  a : ");


  int n_L = 0;
  int n_R = 0;
  for (int i = 0; i < n; i++) {

    if (verbose) log_trace("i = %d, a[i] = %f\n", i, a[i]);
    if (_is_zero(a[i])) {
      if (b[i] < 0) {
        return PDM_LP_UNFEASIBLE;
      } else {
        continue;
      }
    }

    double si = b[i] / a[i];

    if (a[i] < 0) {
      L = PDM_MAX(L, si);
      n_L++;
    } else {
      R = PDM_MIN(R, si);
      n_R++;
    }
  }

  if (L > R) {
    // unfeasible problem
    if (verbose) log_trace("unfeasible (L = %f, R = %f)\n", L, R);
    return PDM_LP_UNFEASIBLE;
  }

  else {
    if (c > 0) {
      if (n_R > 0) {
        *x = R;
        if (verbose) log_trace("feasible (L = %f, R = %f), x = R\n", L, R);
        return PDM_LP_FEASIBLE;
      } else {
          *x = R; // ray [L; +inf) is a solution
          if (verbose) log_trace("unbounded (L = %f, R = %f), ray [L; +inf) is a solution\n", L, R);
          return PDM_LP_UNBOUNDED;
        }
      } else {
        if (n_L > 0) {
          *x = L;
          if (verbose) log_trace("feasible (L = %f, R = %f), x = L\n", L, R);
          return PDM_LP_FEASIBLE;
        } else {
          *x = L; // ray (-inf; R] is a solution
          if (verbose) log_trace("unbounded (L = %f, R = %f), ray (-inf; R] is a solution\n", L, R);
          return PDM_LP_UNBOUNDED;
        }
      }

    }

  return PDM_LP_FEASIBLE;


  // // Different cases
  // enum { BOTH,
  //        POSITIVE,
  //        NEGATIVE};

  // int cas; // cases of 1d Seidel's algorithm
  // int count; // count if ai coefficients have same sign

  // // Compute L and R
  // double L = - big_lp;
  // double R = big_lp;

  // for (int i = 0; i < n; i++) {

  //   if (!_is_zero(a[i])) {
  //     if (a[i] > 0) {
  //       double si = b[i] / a[i];
  //       R = PDM_MIN(si, R);
  //       count ++;
  //     }

  //     if (a[i] < 0) {
  //       double si = b[i] / a[i];
  //       L = PDM_MAX(si, L);
  //       count --;
  //     }
  //   } // a[i] is not null

  // }

  // // log_trace("L = %lf and R = %lf\n", L, R);

  // // Get case
  // if (PDM_ABS(count) != n) {
  //   cas = BOTH;
  // } else {
  //   if (count > 0) {
  //     cas = POSITIVE;
  //   } else {
  //     cas = NEGATIVE;
  //   }
  // }

  // // log_trace("case = %d\n", cas);

  // // Procede according to case
  // if (cas == BOTH) {
  //   if (L > R) {
  //     return PDM_LP_UNFEASIBLE;
  //   } else {
  //     if (c > 0) {
  //       *x = R;
  //       return PDM_LP_FEASIBLE;
  //     } else {
  //       *x = L;
  //       return PDM_LP_FEASIBLE;
  //     }
  //   }
  // }

  // if (cas == POSITIVE) {
  //   if (c > 0) {
  //     *x = R;
  //     return PDM_LP_FEASIBLE;
  //   } else {
  //     *x = R;
  //     return PDM_LP_UNBOUNDED;
  //   }
  // }

  // if (cas == NEGATIVE) {
  //   if (c < 0) {
  //     *x = L;
  //     return PDM_LP_FEASIBLE;
  //   } else {
  //     *x = L;
  //     return PDM_LP_UNBOUNDED;
  //   }
  // }

  // return PDM_LP_FEASIBLE;
}

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
 )
{
  if (dim == 1) {
    return _lp_solve_1d(n, a, b, c[0], x);
  }

  // Set for sub-problem
  double sub_a[n*(dim-1)];
  double sub_b[n];
  double sub_c[dim-1];
  double sub_x[dim-1];

    if (verbose) log_trace("\n\ndim = %d\n", dim);

  for (int i = 0; i < n; i++) {

    // i-th constraint
    double *ai = a + dim*i;
    double ax_b = -b[i];
    for (int j = 0; j < dim; j++) {
      ax_b += ai[j] * x[j];
    }

    if (verbose) {
      log_trace("constraint #%d\n", i);
      PDM_log_trace_array_double(x,  dim, "  x  : ");
      PDM_log_trace_array_double(ai, dim, "  ai : ");
      log_trace("  bi = %f\n", b[i]);
      log_trace("  ai.x - bi = %f\n", ax_b);
    }

    // case: i-th constraint is satisfied
    if (ax_b <= 0) {
      continue;
    }

    // case: i-th constraint not satisfied
    // --> find the largest component
    int    jmax   = -1;
    double aijmax = 0.;
    for (int j = 0; j < dim; j++) {
      double aij = PDM_ABS(ai[j]);
      if (aij > aijmax) {
        jmax = j;
        aijmax = aij;
      }
    }

    if (verbose) {
      if (jmax < 0) {
        log_trace("  jmax = %d\n", jmax);
      }
      else {
        log_trace("  jmax = %d, aijmax = %f\n", jmax, ai[jmax]);
      }
    }

    // --> null constraint
    if (jmax < 0) {
      if (b[i] < 0) {
        return PDM_LP_UNFEASIBLE;
      } else {
        continue;
      }
    }

    // --> project previous constraints to lower dimension
    double inv_aijmax = 1. / ai[jmax];

    for (int k = 0; k < i; k++) {
      double *ak     = a     + dim    *k;
      double *sub_ak = sub_a + (dim-1)*k;
      int jj = 0;
      for (int j = 0; j < dim; j++) {
        if (j == jmax) continue;
        sub_ak[jj] = ak[j] - ak[jmax]*ai[j]*inv_aijmax;
        jj++;
      }

      sub_b[k] = b[k] - ak[jmax]*b[i]*inv_aijmax;
    }

    // --> project c to lower dimension
    int jj = 0;
    double mag_sub_c = 0.;
    for (int j = 0; j < dim; j++) {
      if (j == jmax) continue;
      sub_c[jj] = c[j] - c[jmax]*ai[j]*inv_aijmax;
      mag_sub_c += sub_c[jj]*sub_c[jj];
      sub_x[jj] = x[j];
      jj++;
    }

    if (verbose) PDM_log_trace_array_double(sub_c, dim-1, "  sub_c : ");
    if (mag_sub_c < eps_lp) {
      if (verbose) log_trace("c is null !\n");
    }

    PDM_lp_status_t stat = PDM_lp_solve_nd(dim-1,
                                           i,
                                           sub_a,
                                           sub_b,
                                           sub_c,
                                           sub_x);

    if (stat == PDM_LP_UNFEASIBLE) {
      return PDM_LP_UNFEASIBLE;
    }

    // back substitution
    jj = 0;
    double aix = 0.;
    for (int j = 0; j < dim; j++) {
      if (j == jmax) continue;
      x[j] = sub_x[jj];
      aix += ai[j] * x[j];
      jj++;
    }

    x[jmax] = (b[i] - aix) * inv_aijmax;
    if (verbose) PDM_log_trace_array_double(x,  dim, "  new x : ");

  } // end loop on constraints

  return PDM_LP_FEASIBLE;

}

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
)
{
  // objective function
  double c[3] = {1., 1., 1.};

  // constraints
  double a[(n_plane+6)*3]; // 6 for the box
  double b[n_plane+6];
  a[ 0] = -1.; a[ 1] =  0.; a[ 2] =  0.;
  a[ 3] =  0.; a[ 4] = -1.; a[ 5] =  0.;
  a[ 6] =  0.; a[ 7] =  0.; a[ 8] = -1.;
  a[ 9] =  1.; a[10] =  0.; a[11] =  0.;
  a[12] =  0.; a[13] =  1.; a[14] =  0.;
  a[15] =  0.; a[16] =  0.; a[17] =  1.;
  for (int i = 0; i < 3; i++) {
    b[i]   = -box_extents[i];
    b[i+3] = box_extents[i+3];
  }
  for (int iplane = 0; iplane < n_plane; iplane++) {
    for (int i = 0; i < 3; i++) {
      a[18 + 3*iplane + i] = -plane_normal[3*iplane + i];
    }
    b[6 + iplane] = -PDM_DOT_PRODUCT(plane_normal + 3*iplane, plane_origin + 3*iplane);
  }

  // Solve LP problem
  double x[3];
  for (int i = 0; i < 3; i++) {
    if (c[i] > 0) {
      x[i] = box_extents[i+3];
    } else {
      x[i] = box_extents[i];
    }
  }

  PDM_lp_status_t stat =  PDM_lp_solve_nd(3, n_plane+6, a, b, c, x);

  return (stat != PDM_LP_UNFEASIBLE);

}
