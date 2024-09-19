/*============================================================================
 * Morton encoding for 2D or 3D coordinates.
 *============================================================================*/
/*
  This file is part of the "Finite Volume Mesh" library, intended to provide
  finite volume mesh and associated fields I/O and manipulation services.

  Copyright (C) 2011       ONERA
  Copyright (C) 2008-2010  EDF

  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License as published by the Free Software Foundation; either
  version 2.1 of the License, or (at your option) any later version.

  This library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public
  License along with this library; if not, write to the Free Software
  Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
*/

/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <assert.h>
#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stdbool.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "pdm_priv.h"
#include "pdm_morton.h"
#include "pdm_hilbert.h"
#include "pdm_binary_search.h"
#include "pdm_printf.h"
#include "pdm_error.h"
#include "pdm_logging.h"

/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force back Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */

/*=============================================================================
 * Local Macro definitions
 *============================================================================*/

// static const PDM_morton_int_t PDM_morton_max_level = 31u;

/*=============================================================================
 * Static global variables
 *============================================================================*/

static const double  PDM_morton_distrib_tol = 0.10;

/* Max. number of sub-iterations to get a well-balanced distribution */
static const int PDM_morton_distrib_n_iter_max = 5;

static const int _sampling_factors[4] = {1, /* OD */
                                         2, /* 1D */
                                         2, /* 2D */
                                         4, /* 3D */};

static const int _3d_children[8][3] = {{0, 0, 0},    /* child 1 */
                                       {0, 0, 1},    /* 2 */
                                       {0, 1, 0},    /* 3 */
                                       {0, 1, 1},    /* 4 */
                                       {1, 0, 0},    /* 5 */
                                       {1, 0, 1},    /* 6 */
                                       {1, 1, 0},    /* 7 */
                                       {1, 1, 1}};   /* 8 */

static const int _2d_children[4][2] = {{0, 0},   /* child 1 */
                                       {0, 1},   /* 2 */
                                       {1, 0},   /* 3 */
                                       {1, 1}};  /* 4 */

static const int _1d_children[2][1] = {{0},   /* child 1 */
                                       {1}};  /* 2 */

/*============================================================================
 * Private function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Transorm local extents to global extents.
 *
 * parameters:
 *   dim       <-- spatial dimension (1, 2, or 3)
 *   g_extents <-> global extents (size: dim*2)
 *   comm      <-- associated MPI communicator
 *---------------------------------------------------------------------------*/

static void
_local_to_global_extents(int         dim,
                         double  extents[],
                         PDM_MPI_Comm    comm)
{
  int i;
  double  l_min[3], l_max[3];

  for (i = 0; i < dim; i++) {
    l_min[i] = extents[i];
    l_max[i] = extents[i + dim];
  }

  PDM_MPI_Allreduce(l_min, extents, dim, PDM_MPI_DOUBLE, PDM_MPI_MIN, comm);
  PDM_MPI_Allreduce(l_max, extents + dim, dim, PDM_MPI_DOUBLE, PDM_MPI_MAX, comm);
}

/*----------------------------------------------------------------------------
 * Test if Morton code "a" is greater or equal to Morton code "b"
 *
 * parameters:
 *   code_a <-- first Morton code to compare
 *   code_b <-- second Morton code to compare
 *
 * returns:
 *  true or false
 *----------------------------------------------------------------------------*/

inline static _Bool
_a_ge_b(PDM_morton_code_t  code_a,
        PDM_morton_code_t  code_b)
{
  int i, a, b, a_diff, b_diff;
  int l = PDM_MAX(code_a.L, code_b.L);

  a_diff = l - code_a.L;
  b_diff = l - code_b.L;

  if (a_diff > 0) {
    code_a.L = l;
    code_a.X[0] = code_a.X[0] << a_diff;
    code_a.X[1] = code_a.X[1] << a_diff;
    code_a.X[2] = code_a.X[2] << a_diff;
  }

  if (b_diff > 0) {
    code_b.L = l;
    code_b.X[0] = ((code_b.X[0] + 1) << b_diff) - 1;
    code_b.X[1] = ((code_b.X[1] + 1) << b_diff) - 1;
    code_b.X[2] = ((code_b.X[2] + 1) << b_diff) - 1;
  }

  i = l - 1;
  while (i > 0) {
    if (   code_a.X[0] >> i != code_b.X[0] >> i
        || code_a.X[1] >> i != code_b.X[1] >> i
        || code_a.X[2] >> i != code_b.X[2] >> i)
      break;
    i--;
  }

  a =   ((code_a.X[0] >> i) % 2) * 4
      + ((code_a.X[1] >> i) % 2) * 2
      + ((code_a.X[2] >> i) % 2);
  b =   ((code_b.X[0] >> i) % 2) * 4
      + ((code_b.X[1] >> i) % 2) * 2
      + ((code_b.X[2] >> i) % 2);
 return (a >= b) ? true : false;
}

/*----------------------------------------------------------------------------
 * Test if Morton code "a" is equal to Morton code "b"
 *
 * parameters:
 *   code_a <-- first Morton code to compare
 *   code_b <-- second Morton code to compare
 *
 * returns:
 *  true or false
 *----------------------------------------------------------------------------*/

inline static _Bool
_a_eq_b(PDM_morton_code_t  code_a,
        PDM_morton_code_t  code_b)
{
  int i, a, b, a_diff, b_diff;
  int l = PDM_MAX(code_a.L, code_b.L);

  a_diff = l - code_a.L;
  b_diff = l - code_b.L;

  if (a_diff > 0) {
    code_a.L = l;
    code_a.X[0] = code_a.X[0] << a_diff;
    code_a.X[1] = code_a.X[1] << a_diff;
    code_a.X[2] = code_a.X[2] << a_diff;
  }

  if (b_diff > 0) {
    code_b.L = l;
    code_b.X[0] = ((code_b.X[0] + 1) << b_diff) - 1;
    code_b.X[1] = ((code_b.X[1] + 1) << b_diff) - 1;
    code_b.X[2] = ((code_b.X[2] + 1) << b_diff) - 1;
  }

  i = l - 1;
  while (i > 0) {
    if (   code_a.X[0] >> i != code_b.X[0] >> i
        || code_a.X[1] >> i != code_b.X[1] >> i
        || code_a.X[2] >> i != code_b.X[2] >> i)
      break;
    i--;
  }

  a =   ((code_a.X[0] >> i) % 2) * 4
      + ((code_a.X[1] >> i) % 2) * 2
      + ((code_a.X[2] >> i) % 2);
  b =   ((code_b.X[0] >> i) % 2) * 4
      + ((code_b.X[1] >> i) % 2) * 2
      + ((code_b.X[2] >> i) % 2);

  return (a == b) ? true : false;
}

/*----------------------------------------------------------------------------
 * Test if Morton code "a" is greater than Morton code "b"
 *
 * parameters:
 *   code_a <-- first Morton code to compare
 *   code_b <-- second Morton code to compare
 *
 * returns:
 *  true or false
 *----------------------------------------------------------------------------*/

inline static _Bool
_a_gt_b(PDM_morton_code_t  code_a,
        PDM_morton_code_t  code_b)
{
  int i, a, b, a_diff, b_diff;
  int l = PDM_MAX(code_a.L, code_b.L);

  a_diff = l - code_a.L;
  b_diff = l - code_b.L;

  if (a_diff > 0) {
    code_a.L = l;
    code_a.X[0] = code_a.X[0] << a_diff;
    code_a.X[1] = code_a.X[1] << a_diff;
    code_a.X[2] = code_a.X[2] << a_diff;
  }

  if (b_diff > 0) {
    code_b.L = l;
    code_b.X[0] = ((code_b.X[0] + 1) << b_diff) - 1;
    code_b.X[1] = ((code_b.X[1] + 1) << b_diff) - 1;
    code_b.X[2] = ((code_b.X[2] + 1) << b_diff) - 1;
  }

  i = l - 1;
  while (i > 0) {

    if (   code_a.X[0] >> i != code_b.X[0] >> i
        || code_a.X[1] >> i != code_b.X[1] >> i
        || code_a.X[2] >> i != code_b.X[2] >> i)
      break;
    i--;
  }

   a =   ((code_a.X[0] >> i) % 2) * 4
       + ((code_a.X[1] >> i) % 2) * 2
       + ((code_a.X[2] >> i) % 2);
   b =   ((code_b.X[0] >> i) % 2) * 4
       + ((code_b.X[1] >> i) % 2) * 2
       + ((code_b.X[2] >> i) % 2);

  return (a > b) ? true : false;
}



/*----------------------------------------------------------------------------
 * Test if Morton code "a" is greater than Morton code "b" (compare anchors)
 *
 * parameters:
 *   code_a <-- first Morton code to compare
 *   code_b <-- second Morton code to compare
 *
 * returns:
 *  true or false
 *----------------------------------------------------------------------------*/

inline static _Bool
_a_gtmin_b(PDM_morton_code_t  code_a,
           PDM_morton_code_t  code_b)
{
  int i, a, b, a_diff, b_diff;
  int l = PDM_MAX(code_a.L, code_b.L);

  a_diff = l - code_a.L;
  b_diff = l - code_b.L;

  if (a_diff > 0) {
    code_a.L = l;
    code_a.X[0] = code_a.X[0] << a_diff;
    code_a.X[1] = code_a.X[1] << a_diff;
    code_a.X[2] = code_a.X[2] << a_diff;
  }

  if (b_diff > 0) {
    code_b.L = l;
    code_b.X[0] = code_b.X[0] << b_diff;
    code_b.X[1] = code_b.X[1] << b_diff;
    code_b.X[2] = code_b.X[2] << b_diff;
  }

  i = l - 1;
  while (i > 0) {

    if (   code_a.X[0] >> i != code_b.X[0] >> i
        || code_a.X[1] >> i != code_b.X[1] >> i
        || code_a.X[2] >> i != code_b.X[2] >> i)
      break;
    i--;
  }

   a =   ((code_a.X[0] >> i) % 2) * 4
       + ((code_a.X[1] >> i) % 2) * 2
       + ((code_a.X[2] >> i) % 2);
   b =   ((code_b.X[0] >> i) % 2) * 4
       + ((code_b.X[1] >> i) % 2) * 2
       + ((code_b.X[2] >> i) % 2);

  return (a > b) ? true : false;
}



/*----------------------------------------------------------------------------
 * Build a heap structure or order a heap structure.
 *
 * parameters:
 *  parent       <--  parent id in the Morton code list
 *  n_codes      <--  number of codes to work with
 *  morton_codes <->  list of Morton codes to work with
 *----------------------------------------------------------------------------*/

static void
_descend_morton_heap(PDM_g_num_t        parent,
                     int                n_codes,
                     PDM_morton_code_t  morton_codes[])
{
  PDM_morton_code_t  tmp;
  PDM_g_num_t child = 2 * parent + 1;

  while (child < n_codes) {

    if (child + 1 < n_codes)
      if (_a_gt_b(morton_codes[child + 1], morton_codes[child]))
        child++;

    if (_a_ge_b(morton_codes[parent], morton_codes[child])) return;

    tmp = morton_codes[parent];
    morton_codes[parent] = morton_codes[child];
    morton_codes[child] = tmp;
    parent = child;
    child = 2 * parent + 1;

  } /* End of while */

}

/*----------------------------------------------------------------------------
 * Convert a double into a Morton code.
 *
 * parameters:
 *   dim       <-- 2D or 3D
 *   input     <-- double to convert
 *   level     <-- level of the grid on which the code has to be built
 *
 * returns:
 *  a Morton code associated to the input.
 *----------------------------------------------------------------------------*/

inline static PDM_morton_code_t
_double_to_code(int     dim,
                double  input,
                int     level)
{
  int  l, child_id;
  PDM_morton_code_t  code;
  double coords[3] = {0.0, 0.0, 0.0};
  double l_mult = 1.0;

  const int max_level = 15; /* no more than 52 bits in mantissa / 3 */

  /* Build associated Morton code */

  code.L = max_level;

  if (input <= 0.0) {
    coords[0] = 0.0;
    coords[1] = 0.0;
    coords[2] = 0.0;
  }

  else if (input >= 1.0) {
    coords[0] = 1.0;
    coords[1] = 1.0;
    coords[2] = 1.0;
  }

  else if (dim == 3) {
    for (l = 0; l < max_level; l++) {
      l_mult *= 0.5;
      child_id = (int)(input*8);
      if (child_id > 7) child_id = 7;
      input = input*8 - child_id;
      coords[0] += child_id/4 * l_mult;
      coords[1] += (child_id%4)/2 * l_mult;
      coords[2] += child_id%2 * l_mult;
    }
  }

  else if (dim == 2) {
    coords[2] = 0;
    for (l = 0; l < max_level; l++) {
      l_mult *= 0.5;
      child_id = (int)(input*4);
      if (child_id > 3) child_id = 3;
      input = input*4 - child_id;
      coords[0] += child_id/2 * l_mult;
      coords[1] += child_id%2 * l_mult;
    }
  }

  else if (dim == 1) {
    coords[1] = 0;
    coords[2] = 0;
    for (l = 0; l < max_level; l++) {
      l_mult *= 0.5;
      child_id = (int)(input*2);
      if (child_id > 1) child_id = 1;
      input = input*2 - child_id;
      coords[0] += child_id * l_mult;
    }
  }

  code = PDM_morton_encode(dim, level, coords);

  return code;
}

/*----------------------------------------------------------------------------
 * Build a heap structure or order a heap structure with a working array
 * to save the ordering.
 *
 * parameters:
 *  parent       <-- parent id in the Morton code list
 *  n_codes      <-- number of codes to work with
 *  morton_codes <-- list of Morton codes to work with
 *  order        <-> working array to save the ordering
 *----------------------------------------------------------------------------*/

static void
_descend_morton_heap_with_order(PDM_g_num_t                 parent,
                                int                 n_codes,
                                const PDM_morton_code_t   morton_codes[],
                                int                *order)
{
  int          tmp;
  PDM_g_num_t  child = 2 * parent + 1;

  while (child < n_codes) {

    if (child + 1 < n_codes) {
      if (_a_gt_b(morton_codes[order[child + 1]],
                  morton_codes[order[child]]))
        child++;
    }

    if (_a_ge_b(morton_codes[order[parent]],
                morton_codes[order[child]]))
      return;

    tmp = order[parent];
    order[parent] = order[child];
    order[child] = tmp;
    parent = child;
    child = 2 * parent + 1;

  } /* End while */
}

/*----------------------------------------------------------------------------
 * Evaluate a distribution array.
 *
 * parameters:
 *   n_ranges     <-- Number of ranges in the distribution
 *   distribution <-- Number of elements associated to each range of
 *                    the distribution
 *   optim        <-- Optimal count in each range
 *
 * returns:
 *   a fit associated to the distribution. If fit = 0,
 *   distribution is perfect.
 *----------------------------------------------------------------------------*/

static double
_evaluate_distribution(int          n_ranges,
                       double      *distribution,
                       double       optim)
{
  int  i;
  double  d_low = 0, d_up = 0, fit = 0;

  /*
     d_low is the max gap between the distribution count and the optimum when
     distribution is lower than optimum.
     d_up is the max gap between the distribution count and the optimum when
     distribution is greater than optimum.
  */

  for (i = 0; i < n_ranges; i++) {

    if (distribution[i] > optim)
      d_up = PDM_MAX(d_up, distribution[i] - optim);
    else
      d_low = PDM_MAX(d_low, optim - distribution[i]);

  }

  fit = (d_up + d_low) / optim;

#if 0 && defined(dbg_enabled) && !defined(Ndbg_enabled)
  /* if (cs_glob_rank_id <= 0) */
  /*   PDM_printf("<DISTRIBUTION EVALUATION> optim: %g, fit: %g\n", */
  /*              optim, fit); */
#endif

  return  fit;
}

static uint64_t
_interleave
(
 const int               dimension,
 const PDM_morton_code_t c
 )
{
  int k = 0;
  uint64_t i = 0;
  for (PDM_morton_int_t l = 0; l < c.L; l++) {
    for (int j = dimension-1; j >= 0; j--) {
      uint64_t a = (c.X[j] >> l) & 1l;
      i += a << k;
      k++;
    }
  }

  return i;
}

static double
_code_to_double
(
 const int               dimension,
 const PDM_morton_code_t code
 )
{
  uint64_t i = _interleave(dimension, code);
  return (double) i;
}


/*----------------------------------------------------------------------------
 * Define a global distribution associated to a sampling array i.e. count
 * the number of elements in each range.
 *
 * parameters:
 *   dim          <-- 2D or 3D
 *   n_ranks      <-- number of ranks (= number of ranges)
 *   gmax_level   <-- level on which Morton encoding is build
 *   gsum_weight  <-- global sum of all weightings
 *   n_codes      <-- local number of Morton codes
 *   morton_codes <-- local list of Morton codes to distribute
 *   weight       <-- weighting related to each code
 *   order        <-- ordering array
 *   sampling     <-- sampling array
 *   c_freq       <-> pointer to the cumulative frequency array
 *   g_distrib    <-> pointer to a distribution array
 *   comm         <-- mpi communicator
 *----------------------------------------------------------------------------*/

static void
_define_rank_distrib(int                      dim,
                     int                      n_ranks,
                     int                      gmax_level,
                     double                   gsum_weight,
                     int                      n_codes,
                     const PDM_morton_code_t  morton_codes[],
                     const double             weight[],
                     const int                order[],
                     const double             sampling[],
                     double                   cfreq[],
                     double                   g_distrib[],
                     PDM_MPI_Comm             comm)
{

  const int  sampling_factor = _sampling_factors[dim];
  const int  n_samples = sampling_factor * n_ranks;

  /* Initialization */
  double   *l_distrib = (double  *) malloc (n_samples * sizeof(double));

  for (int id = 0; id < n_samples; id++) {
    l_distrib[id] = 0;
    g_distrib[id] = 0;
  }

  if(order == NULL) {
    gmax_level = 21;
    uint64_t max     = 1l << (gmax_level * dim);
    double   inv_max = 1./((double) max);

    /* morton_codes codes if not ordered !!!! */
    for (int i = 0; i < n_codes; i++) {
      PDM_morton_code_t lcode;
      PDM_morton_copy (morton_codes[i], &lcode);
      PDM_morton_assign_level (&lcode, gmax_level);

      double val1d = _code_to_double(dim, lcode) * inv_max;
      size_t t_bucket = PDM_hilbert_quantile_search(n_samples, val1d, sampling);
      l_distrib[t_bucket] += weight[i];
    } /* End of loop on elements */
  } else {
    /* morton_codes are supposed to be ordered */
    PDM_morton_code_t  sample_code;
    int  bucket_id = 1;
    sample_code = _double_to_code(dim, sampling[bucket_id], gmax_level);
    for (int i = 0; i < n_codes; i++) {
      PDM_g_num_t   o_id = order[i];
      if (_a_ge_b(sample_code, morton_codes[o_id])) {
        l_distrib[bucket_id - 1] += weight[o_id];
      }
      else {
        while (_a_gt_b(morton_codes[o_id], sample_code)) {
          bucket_id++;
          assert(bucket_id < n_samples + 1);
          sample_code = _double_to_code(dim, sampling[bucket_id], gmax_level);
        }
        l_distrib[bucket_id - 1] += weight[o_id];
      }
    } /* End of loop on elements */
  }

  /* Define the global distribution */
  PDM_MPI_Allreduce(l_distrib, g_distrib, n_samples, PDM_MPI_DOUBLE, PDM_MPI_SUM, comm);

  free(l_distrib);

  /* Define the cumulative frequency related to g_distribution */

  cfreq[0] = 0.;
  for (int id = 0; id < n_samples; id++) {
    cfreq[id+1] = cfreq[id] + (double)g_distrib[id]/(double)gsum_weight;
  }
  cfreq[n_samples] = 1.0;

#if 0 && defined(dbg_enabled) && !defined(dbg_enabled) /* For dbg_enabledging purpose only */

  /* if (cs_glob_rank_id <= 0) { */

  /*   FILE  *dbg_file = NULL; */
  /*   char  *rfilename = NULL; */
  /*   int  len; */
  /*   static int  loop_id1 = 0; */

  /*   len = strlen("DistribOutput_l.dat")+1+2; */
  /*   rfilename = (char *) malloc(len * sizeof(char)); */
  /*   sprintf(rfilename, "DistribOutput_l%02d.dat", loop_id1); */

  /*   loop_id1++; */

  /*   dbg_file = fopen(rfilename, "w"); */

  /*   fprintf(dbg_file, */
  /*           "# Sample_id  |  OptCfreq  |  Cfreq  |  Sampling  |" */
  /*           "Global Distrib\n"); */
  /*   for (i = 0; i < n_samples; i++) */
  /*     fprintf(dbg_file, "%8d %15.5f %15.10f %15.10f %10u\n", */
  /*             i, (double)i/(double)n_samples, cfreq[i], */
  /*             sampling[i], g_distrib[i]); */
  /*   fprintf(dbg_file, "%8d %15.5f %15.10f %15.10f %10u\n", */
  /*           i, 1.0, 1.0, 1.0, 0); */

  /*   fclose(dbg_file); */
  /*   free(rfilename); */

  /* } */

#endif /* dbg_enabledging output */

  /* Convert global distribution from n_samples to n_ranks */

  for (int rank_id = 0; rank_id < n_ranks; rank_id++) {

    PDM_g_num_t   sum = 0;
    int   shift = rank_id * sampling_factor;

    for (int id = 0; id < sampling_factor; id++)
      sum += g_distrib[shift + id];
    g_distrib[rank_id] = sum;

  } /* End of loop on ranks */

#if 0 && defined(dbg_enabled) && !defined(Ndbg_enabled) /* Sanity check in dbg_enabled */
  {
    PDM_g_num_t   sum = 0;
    for (rank_id = 0; rank_id < n_ranks; rank_id++)
      sum += g_distrib[rank_id];

    if (sum != gsum_weight) {
      PDM_error(__FILE__, __LINE__, 0,
                "Error while computing global distribution.\n"
                "sum = %u and gsum_weight = %u\n",
                sum, gsum_weight);
      abort();
    }
  }
#endif /* sanity check */

}

/*----------------------------------------------------------------------------
 * Update a distribution associated to sampling to assume a well-balanced
 * distribution of the leaves of the tree.
 *
 * parameters:
 *   dim      <-- 1D, 2D or 3D
 *   n_ranks  <-- number of ranks (= number of ranges)
 *   c_freq   <-> cumulative frequency array
 *   sampling <-> pointer to pointer to a sampling array
 *   comm     <-- mpi communicator
 *----------------------------------------------------------------------------*/

static void
_update_sampling(int      dim,
                 int      n_ranks,
                 double   c_freq[],
                 double  *sampling[])
{
  int  i, j, next_id;
  double  target_freq, f_high, f_low, delta;
  double  s_low, s_high;

  double  *new_sampling = NULL, *_sampling = *sampling;

  const int  sampling_factor = _sampling_factors[dim];
  const int  n_samples = sampling_factor * n_ranks;
  const double  unit = 1/(double)n_samples;

  /* Compute new_sampling */

  new_sampling = (double *) malloc((n_samples + 1) * sizeof(double));

  new_sampling[0] = _sampling[0];
  next_id = 1;

  for (i = 0; i < n_samples; i++) {

    target_freq = (i+1)*unit;

    /* Find the next id such as c_freq[next_id] >= target_freq */

    for (j = next_id; j < n_samples + 1; j++) {
      if (c_freq[j] >= target_freq) {
        next_id = j;
        break;
      }
    }

    /* Find new s such as new_s is equal to target_freq by
       a linear interpolation */

    f_low = c_freq[next_id-1];
    f_high = c_freq[next_id];

    s_low = _sampling[next_id-1];
    s_high = _sampling[next_id];

    if (f_high - f_low > 0) {
      delta = (target_freq - f_low) * (s_high - s_low) / (f_high - f_low);
      new_sampling[i+1] = s_low + delta;
    }
    else /* f_high = f_low */
      new_sampling[i+1] = s_low + 0.5 * (s_low + s_high);

#if 0 && defined(dbg_enabled) && !defined(Ndbg_enabled)
    /* PDM_printf(" <_update_distrib> (rank: %d) delta: %g, target: %g," */
    /*            " next_id: %d, f_low: %g, f_high: %g, s_low: %g, s_high: %g\n" */
    /*            "\t => new_sampling: %g\n", */
    /*            cs_glob_rank_id, delta, target_freq, next_id, */
    /*            f_low, f_high, s_low, s_high, new_sampling[i+1]); */
#endif

  } /* End of loop on samples */

  new_sampling[n_samples] = 1.0;

  free(_sampling);

  /* Return pointers */

  *sampling = new_sampling;
}

/*----------------------------------------------------------------------------
 * Compute a sampling array which assumes a well-balanced distribution of
 * leaves of the tree among the ranks.
 *
 * parameters:
 *   dim          <--  2D or 3D
 *   n_ranks      <--  number of ranks
 *   gmax_level   <--  level on which Morton encoding is build
 *   n_codes      <--  local number of Morton ids
 *   morton_codes <--  local list of Morton ids to distribute
 *   weight       <--  weighting related to each code
 *   order        <--  ordering array
 *   sampling     <-> pointer to pointer to a sampling array
 *   comm         <--  mpi communicator
 *
 * returns:
 *   fit associated to the returned sampling array
 *----------------------------------------------------------------------------*/

static double
_bucket_sampling(int                      dim,
                 int                      n_ranks,
                 int                      gmax_level,
                 int                      n_codes,
                 const PDM_morton_code_t  morton_codes[],
                 const double             weight[],
                 const int                order[],
                 double                  *sampling[],
                 PDM_MPI_Comm             comm)
{
  int  i, n_iters;
  int   j;
  double  fit, best_fit, optim;

  double   lsum_weight = 0, gsum_weight = 0;
  double   *distrib = NULL;
  double  *cfreq = NULL, *best_sampling = NULL;
  double  *_sampling = *sampling;

  const int  sampling_factor = _sampling_factors[dim];
  const int  n_samples = PDM_MAX(1, sampling_factor * n_ranks);
  const double  unit = 1/(double)n_samples;

  /* Compute the global number of elements and the optimal number of elements
     on each rank */

  for (j = 0; j < n_codes; j++) {
    lsum_weight += weight[j];
  }

  PDM_MPI_Allreduce(&lsum_weight, &gsum_weight, 1,  PDM_MPI_DOUBLE, PDM_MPI_SUM, comm);

  optim = (double)gsum_weight / (double)n_ranks;

  /* Define a naive sampling (uniform distribution) */

  for (i = 0; i < n_samples + 1; i++) {
    _sampling[i] = i*unit;
  }

  /* Define the distribution associated to the current sampling array */

  distrib = (double      *) malloc(n_samples       * sizeof(double     ));
  cfreq   = (double      *) malloc((n_samples + 1) * sizeof(double     ));

  _define_rank_distrib(dim,
                       n_ranks,
                       gmax_level,
                       gsum_weight,
                       n_codes,
                       morton_codes,
                       weight,
                       order,
                       _sampling,
                       cfreq,
                       distrib,
                       comm);

  /* Initialize best choice */

  fit = _evaluate_distribution(n_ranks, distrib, optim);
  best_fit = fit;

  best_sampling = malloc((n_samples + 1) * sizeof(double));

  for (i = 0; i < n_samples + 1; i++) {
    best_sampling[i] = _sampling[i];
  }

  /* Loop to get a better sampling array */

  for (n_iters = 0;
       (   n_iters < PDM_morton_distrib_n_iter_max
        && fit > PDM_morton_distrib_tol);
       n_iters++)  {

    _update_sampling(dim, n_ranks, cfreq, &_sampling);

    /* Compute the new distribution associated to the new sampling */

    _define_rank_distrib(dim,
                         n_ranks,
                         gmax_level,
                         gsum_weight,
                         n_codes,
                         morton_codes,
                         weight,
                         order,
                         _sampling,
                         cfreq,
                         distrib,
                         comm);

    fit = _evaluate_distribution(n_ranks, distrib, optim);

    /* Save the best sampling array and its fit */

    if (fit < best_fit) {

      best_fit = fit;
      for (i = 0; i < n_samples + 1; i++)
        best_sampling[i] = _sampling[i];

    }

  } /* End of while */

#if 0 && defined(dbg_enabled) && !defined(Ndbg_enabled)
  /* if (cs_glob_rank_id <= 0) */
  /*   PDM_printf("\n  <_bucket_sampling> n_iter: %d, opt: %g, best_fit: %g\n", */
  /*              n_iters, optim, best_fit); */
#endif

  /* Free memory */

  free(cfreq);
  free(distrib);
  free(_sampling);

  *sampling = best_sampling;

  return best_fit;
}




static int
_intersect_node_box
(
 const int                dim,
 const PDM_morton_code_t  node,
 const PDM_morton_code_t  box_min,
 const PDM_morton_code_t  box_max,
 int                     *inside
 )
{
  const int dbg_enabled = 0;
  *inside = 1;

  assert (box_min.L >= node.L);
  if (dbg_enabled) {
    printf("node: L = %u, X = %u %u %u\n", node.L, node.X[0], node.X[1], node.X[2]);
  }

  const PDM_morton_int_t level_diff = box_min.L - node.L;

  const PDM_morton_int_t side = 1 << level_diff;

  for (int i = 0; i < dim; i++) {
    PDM_morton_int_t xmin = side * node.X[i];
    PDM_morton_int_t xmax = xmin + side;

    if (xmin > box_max.X[i]+1 || xmax < box_min.X[i]) {
      if (dbg_enabled) {
        //printf("\t not intersecting (dim %d, xmin = %u, box_max = %u, box_min = %u, xmax = %u\n", i, xmin, box_max.X[i]+1, box_min.X[i], xmax);
        double s = 1. / pow(2., box_min.L);
        printf("\t not intersecting (dim %d, xmin = %f, box_max = %f, box_min = %f, xmax = %f\n", i, xmin*s, (box_max.X[i]+1)*s, box_min.X[i]*s, xmax*s);
      }
      return 0;
    } else if (xmin < box_min.X[i] || xmax > box_max.X[i]+1) {
      *inside = 0;
    };
  }

  if (dbg_enabled) {
    printf("\t intersecting\n");
  }

  return 1;
}


inline static double
_min_dist2
(
 const int               dim,
 const PDM_morton_code_t node,
 const double            point[],
 const double            d[]
 )
{
  double min_dist2 = 0.;
  double delta, xmin, xmax;
  double side = 1. / ((double) (1 << node.L));

  for (int i = 0; i < dim; i++) {
    xmin = side * node.X[i];
    xmax = xmin + side;

    if (point[i] > xmax) {
      delta = d[i] * (point[i] - xmax);
      min_dist2 += delta * delta;
    } else if (point[i] < xmin) {
      delta = d[i] * (point[i] - xmin);
      min_dist2 += delta * delta;
    }
  }

  return min_dist2;
}

/*============================================================================
 * Public function definitions
 *============================================================================*/

/**
 * \brief Determine the global extents associated with a set of coordinates
 *
 * \param [in]   dim        Spatial dimension
 * \param [in]   n_coords   Local number of coordinates
 * \param [in]   coords     Coordinates, interleaved (size = \ref dim * \ref n_coords)
 * \param [out]  g_extents  Global extents (size = 2 * \ref dim)
 * \param [in]   comm       Associated MPI communicator
 *
 */

void
PDM_morton_get_coord_extents(int          dim,
                             size_t       n_coords,
                             const double coords[],
                             double       g_extents[],
                             PDM_MPI_Comm comm)
{
  size_t  i, j;

  /* Get global min/max coordinates */

  for (j = 0; j < (size_t)dim; j++) {
    g_extents[j]       = DBL_MAX;
    g_extents[j + dim] = -DBL_MAX;
  }

  for (i = 0; i < n_coords; i++) {
    for (j = 0; j < (size_t)dim; j++) {
      if (coords[i*dim + j] < g_extents[j])
        g_extents[j] = coords[i*dim + j];
      if (coords[i*dim + j] > g_extents[j + dim])
        g_extents[j + dim] = coords[i*dim + j];
    }
  }

  if (comm != PDM_MPI_COMM_NULL)
    _local_to_global_extents(dim, g_extents, comm);

}



/**
 * \brief Determine the global extents associated with a set of local extents
 *
 * \param [in]   dim        Spatial dimension
 * \param [in]   n_extents  Local number of extents
 * \param [in]   extents    Local extents (size = 2 * \ref dim * \ref n_extents)
 * \param [out]  g_extents  Global extents (size = 2 * \ref dim)
 * \param [in]   comm       Associated MPI communicator
 *
 */

void
PDM_morton_get_global_extents(int           dim,
                              size_t        n_extents,
                              const double  extents[],
                              double        g_extents[],
                              PDM_MPI_Comm  comm)
{
  size_t  i, j;

  /* Get global min/max coordinates */

  for (i = 0; i < (size_t)dim; i++) {
    g_extents[i]       = DBL_MAX;
    g_extents[i + dim] = -DBL_MAX;
  }

  for (i = 0; i < n_extents; i++) {
    for (j = 0; j < (size_t)dim; j++) {
      g_extents[j]     = PDM_MIN(g_extents[j],
                                extents[i*dim*2 + j]);
      g_extents[j+dim] = PDM_MAX(g_extents[j + dim],
                                extents[i*dim*2 + j + dim]);
    }
  }

  if (comm != PDM_MPI_COMM_NULL)
    _local_to_global_extents(dim, g_extents, comm);

}



/**
 * \brief Build a Morton code according to the level in an octree grid and its coordinates in the grid.
 *
 * \param [in]   dim        Spatial dimension (1, 2 or 3)
 * \param [in]   level      Level in the grid
 * \param [in]   coords     Coordinates in the grid (normalized) (size = \ref dim)
 *
 * \return                  a Morton code
 *
 */

PDM_morton_code_t
PDM_morton_encode(int               dim,
                  PDM_morton_int_t  level,
                  const double      coords[])
{
  int  i;
  PDM_morton_code_t  morton_code;

  PDM_morton_int_t  refinement = 1u << level;

  morton_code.L = level;

  /* Initialize last components for 1D or 2D case */

  morton_code.X[1] = 0;
  morton_code.X[2] = 0;

  for (i = 0; i < dim; i++)
    morton_code.X[i] = (PDM_morton_int_t) PDM_MIN(floor(coords[i]*refinement), refinement - 1);

  return morton_code;
}



/**
 * \brief Encode an array of coordinates.
 *
 * The caller is responsible for freeing the returned array once it is
 * no longer useful.
 *
 * \param [in]   dim        Spatial dimension (1, 2 or 3)
 * \param [in]   level      Level in the grid
 * \param [in]   extents    Coordinate extents for normalization (size = 2 * \ref dim)
 * \param [in]   n_coords   Local number of coordinates
 * \param [in]   coords     Coordinates, interleaved (size = \ref dim * \ref n_coords)
 * \param [out]  m_code     Array of corresponding Morton codes (size = \ref n_coords)
 * \param [out]  d          Normalization vector (dilatation component)
 * \param [out]  s          Normalization vector (translation component)
 *
 */

void
PDM_morton_encode_coords(int                dim,
                         PDM_morton_int_t   level,
                         const double       extents[],
                         size_t             n_coords,
                         const double       coords[],
                         PDM_morton_code_t  m_code[],
                         double             d[3],
                         double             s[3])
{
  size_t i, j;
  double n[3];
  double d_max = 0.0;

  PDM_morton_int_t  refinement = 1u << level;

  for (i = 0; i < (size_t)dim; i++) {
    s[i] = extents[i];
    d[i] = extents[i+dim] - extents[i];
    d_max = PDM_MAX(d_max, d[i]);
  }

  for (i = 0; i < (size_t)dim; i++) { /* Reduce effective dimension */
    if (d[i] < d_max * 1e-10)
      d[i] = d_max * 1e-10;
  }

  switch(dim) {

  case 3:
    for (i = 0; i < n_coords; i++) {
      m_code[i].L = level;
      for (j = 0; j < 3; j++) {
        n[j] = (coords[i*dim + j] - s[j]) / d[j];
        m_code[i].X[j] = (PDM_morton_int_t) PDM_MIN(floor(n[j]*refinement), refinement - 1);
      }
    }
    break;

  case 2:
    for (i = 0; i < n_coords; i++) {
      m_code[i].L = level;
      for (j = 0; j < 2; j++) {
        n[j] = (coords[i*dim + j] - s[j]) / d[j];
        m_code[i].X[j] = (PDM_morton_int_t) PDM_MIN(floor(n[j]*refinement), refinement - 1);
      }
      m_code[i].X[2] = 0;
    }
    break;

  case 1:
    for (i = 0; i < n_coords; i++) {
      m_code[i].L = level;
      n[0] = (coords[i] - s[0]) / d[0];
      m_code[i].X[0] = (PDM_morton_int_t) PDM_MIN(floor(n[0]*refinement), refinement - 1);
      m_code[i].X[1] = 0;
      m_code[i].X[2] = 0;
    }
    break;

  default:
    assert(dim > 0 && dim < 4);
    break;
  }
}



/**
 * \brief Compute the Morton codes of the children of a node
 *
 * Given a Morton code in the grid, compute the Morton codes of its
 * children when refining the grid by one level.
 *
 * \param [in]   dim        Spatial dimension (1, 2 or 3)
 * \param [in]   parent     Morton code associated with the parent
 * \param [out]  children   Array of children Morton codes (size = 2^\ref dim)
 *
 */

void
PDM_morton_get_children(int                dim,
                        PDM_morton_code_t  parent,
                        PDM_morton_code_t  children[])
{
  int  i;
  PDM_morton_code_t  anchor;

  if (dim == 3) {

    for (i = 0; i < 3; i++)
      anchor.X[i] = 2 * parent.X[i];

    for (i = 0; i < 8; i++) {
      children[i].L = parent.L + 1;
      children[i].X[0] = anchor.X[0] + _3d_children[i][0];
      children[i].X[1] = anchor.X[1] + _3d_children[i][1];
      children[i].X[2] = anchor.X[2] + _3d_children[i][2];
    }

  }
  else if (dim == 2) {

    for (i = 0; i < 2; i++)
      anchor.X[i] = 2 * parent.X[i];

    for (i = 0; i < 4; i++) {
      children[i].L = parent.L + 1;
      children[i].X[0] = anchor.X[0] + _2d_children[i][0];
      children[i].X[1] = anchor.X[1] + _2d_children[i][1];
      children[i].X[2] = 0;
    }

  }

  else if (dim == 1) {

    anchor.X[0] = 2 * parent.X[0];

    for (i = 0; i < 2; i++) {
      children[i].L = parent.L + 1;
      children[i].X[0] = anchor.X[0] + _1d_children[i][0];
      children[i].X[1] = 0;
      children[i].X[2] = 0;
    }

  }
}



/**
 * \brief Get local order in a list of Morton codes
 *
 * \param [in]   n_codes       Number of Morton codes to order
 * \param [in]   morton_codes  Array of Morton codes to order
 * \param [out]  order         Pointer to pre-allocated ordering table
 *
 */

void
PDM_morton_local_order(int                     n_codes,
                       const PDM_morton_code_t morton_codes[],
                       int                     order[])
{
  int   i;
  int   tmp;

  assert(n_codes == 0 || morton_codes != NULL);
  assert(n_codes == 0 || order != NULL);

  for (i = 0; i < n_codes; i++)
    order[i] = i;

  /* Build heap */

  for (i = n_codes/2 - 1; (int)i >= 0; i--)
    _descend_morton_heap_with_order(i,  n_codes, morton_codes, order);

  /* Sort array */

  for (i = n_codes - 1; (int)i >= 0; i--) {

    tmp = order[0];
    order[0] = order[i];
    order[i] = tmp;

    _descend_morton_heap_with_order(0, i, morton_codes, order);

  }

#if 0 && defined(dbg_enabled) && !defined(Ndbg_enabled)   /* Check ordering */
  for (i = 1; i < n_codes; i++) {
    if (_a_gt_b(morton_codes[order[i-1]], morton_codes[order[i]])) {
      PDM_error(__FILE__, __LINE__, 0,
              "Id: %u inconsistent: bad ordering of Morton codes.",
              (unsigned)i);
      abort();
    }
  }
#endif
}



/**
 * \brief Sort a local list of Morton codes
 *
 * \param [in]     n_codes       Number of Morton codes to sort
 * \param [in,out] morton_codes  Array of Morton codes to sort
 *
 */

void
PDM_morton_local_sort(int               n_codes,
                      PDM_morton_code_t morton_codes[])
{
  int   i;
  PDM_morton_code_t  tmp;

  /* Build heap */

  for (i = n_codes/2 - 1; (int)i >= 0; i--)
    _descend_morton_heap(i, n_codes, morton_codes);

  /* Sort array */

  for (i = n_codes - 1; (int)i >= 0; i--) {

    tmp = morton_codes[0];
    morton_codes[0] = morton_codes[i];
    morton_codes[i] = tmp;

    _descend_morton_heap(0, i, morton_codes);

  }

#if 0 && defined(dbg_enabled) && !defined(Ndbg_enabled)   /* Check good ordering */
  for (i = 1; i < n_codes; i++) {
    if (_a_gt_b(dim, morton_codes[i - 1], morton_codes[i])) {
      PDM_error(__FILE__, __LINE__, 0,
              "Id: %u inconsistent: bad ordering of Morton codes.",
              (unsigned)i);
      abort();
    }
  }
#endif

}



/**
 * \brief Compare two Morton codes
 *
 * Compare two Morton encoding and check if these two codes are equal,
 * different or shared the same anchor.
 *
 * \param [in]   dim        Spatial dimension (2 or 3)
 * \param [in]   code_a     First Morton code to compare
 * \param [in]   code_b     Second Morton code to compare
 *
 * \return                  A type on the kind of relation between the two Morton codes.
 *
 */

PDM_morton_compare_t
PDM_morton_compare(int                dim,
                   PDM_morton_code_t  code_a,
                   PDM_morton_code_t  code_b)
{
  int i;

  if (code_a.L == code_b.L) {

    for (i = 0; i < dim; i++)
      if (code_a.X[i] != code_b.X[i])
        return PDM_MORTON_DIFFERENT_ID;
    return PDM_MORTON_EQUAL_ID;

  }
  else {

    if (code_a.L < code_b.L) {

      PDM_morton_int_t  delta = code_b.L - code_a.L;

      for (i = 0; i < dim; i++)
        code_a.X[i] = code_a.X[i] << delta;

    }
    else {

      PDM_morton_int_t  delta = code_a.L - code_b.L;

      for (i = 0; i < dim; i++)
        code_b.X[i] = code_b.X[i] << delta;

    }

    for (i = 0; i < dim; i++)
      if (code_a.X[i] != code_b.X[i])
        return PDM_MORTON_DIFFERENT_ID;
    return PDM_MORTON_SAME_ANCHOR;

  }

}



/**
 * \brief Copy Morton code "a" into Morton code "b"
 *
 * \param [in]   a     Morton code to copy
 * \param [out]  b     Morton code receiving the copy
 *
 */

void
PDM_morton_copy (PDM_morton_code_t  a,
                 PDM_morton_code_t  *b)
{
  b->L = a.L;
  for (int i = 0; i < 3; i++) {
    b->X[i] = a.X[i];
  }
}



/**
 * \brief Get the nearest common ancestor of two Morton codes
 *
 * \param [in]   code_a     First Morton code
 * \param [in]   code_b     Second Morton code
 * \param [out]  c          Nearest common ancestor of the two codes
 *
 */

void
PDM_morton_nearest_common_ancestor (PDM_morton_code_t  code_a,
                                    PDM_morton_code_t  code_b,
                                    PDM_morton_code_t *code_c)
{
  int i, a_diff, b_diff;
  int l = PDM_MIN(code_a.L, code_b.L);

  a_diff = code_a.L - l;
  b_diff = code_b.L - l;

  if (a_diff > 0) {
    code_a.L = l;
    code_a.X[0] = code_a.X[0] >> a_diff;
    code_a.X[1] = code_a.X[1] >> a_diff;
    code_a.X[2] = code_a.X[2] >> a_diff;
  }

  if (b_diff > 0) {
    code_b.L = l;
    code_b.X[0] = code_b.X[0] >> b_diff;
    code_b.X[1] = code_b.X[1] >> b_diff;
    code_b.X[2] = code_b.X[2] >> b_diff;
  }

  i = 0;
  while (i <= l) {
    if (   code_a.X[0] >> i == code_b.X[0] >> i
        && code_a.X[1] >> i == code_b.X[1] >> i
        && code_a.X[2] >> i == code_b.X[2] >> i)
      break;
    i++;
  }

  code_c->L = l - i;
  code_c->X[0] = code_a.X[0] >> i;
  code_c->X[1] = code_a.X[1] >> i;
  code_c->X[2] = code_a.X[2] >> i;

}



/**
 * \brief Test if Morton code "a" is greater than Morton code "b"
 *
 * \param [in]   a     First Morton code to compare
 * \param [in]   b     Second Morton code to compare
 *
 * \return             True or false
 *
 */

_Bool
PDM_morton_a_gt_b(PDM_morton_code_t  a,
                  PDM_morton_code_t  b)
{
  return  _a_gt_b(a, b);
}




/**
 * \brief Test if Morton code "a" is greater than Morton code "b" (compare \em anchors)
 *
 * \param [in]   a     First Morton code to compare
 * \param [in]   b     Second Morton code to compare
 *
 * \return             True or false
 *
 */

_Bool
PDM_morton_a_gtmin_b(PDM_morton_code_t  a,
                     PDM_morton_code_t  b)
{
  return  _a_gtmin_b(a, b);
}



/**
 * \brief Test if Morton code "a" is greater than or equal to Morton code "b"
 *
 * \param [in]   a     First Morton code to compare
 * \param [in]   b     Second Morton code to compare
 *
 * \return             True or false
 *
 */

_Bool
PDM_morton_a_ge_b(PDM_morton_code_t  a,
                  PDM_morton_code_t  b)
{
  return  _a_ge_b(a, b);
}



/**
 * \brief Test if Morton code "a" is equal to Morton code "b" at the level = max (level a , level b)
 *
 * \param [in]   a     First Morton code to compare
 * \param [in]   b     Second Morton code to compare
 *
 * \return             True or false
 *
 */

_Bool
PDM_morton_a_eq_b(PDM_morton_code_t  a,
                  PDM_morton_code_t  b)
{
  return  _a_eq_b (a, b);
}



/**
 * \brief Assign a level to a Morton code
 *
 * \param [in,out]  a     Morton code
 * \param [in]      l     Level to assign
 *
 */

void
PDM_morton_assign_level (PDM_morton_code_t  *code,
                         int                l)
{
  int a_diff = l - code->L;

  if (a_diff > 0) {
    code->L = l;
    code->X[0] = code->X[0] << a_diff;
    code->X[1] = code->X[1] << a_diff;
    code->X[2] = code->X[2] << a_diff;
  }

  else if (a_diff < 0) {
    code->L = l;
    a_diff = -a_diff;
    code->X[0] = code->X[0] >> a_diff;
    code->X[1] = code->X[1] >> a_diff;
    code->X[2] = code->X[2] >> a_diff;
  }

}



/**
 * \brief Get the id associated to a Morton code in an array using a binary search.
 *
 * \param [in]  size   Size of the array
 * \param [in]  code   Morton code we are searching for
 * \param [in]  codes  Array of Morton codes
 *
 * \return             Id associated to the given code in the codes array.
 *
 */

int
PDM_morton_binary_search(int                size,
                         PDM_morton_code_t  code,
                         PDM_morton_code_t *codes)
{
  int start = 0;
  int end = size;

  while (end - start > 1) {

    int  middle = (end - start)/2 + start;

    /* PDM_morton_dump( 3, codes[middle]); */
    /* PDM_morton_dump( 3, code); */
    /* printf ("start end middle : %d %d %d\n", start, end, middle); */
    if (_a_gt_b(codes[middle], code))
      end = middle;
    else
      start = middle;

  }

  return start;
}



/**
 * \brief Get the quantile associated to a Morton code using a binary search.
 *
 * No check is done to ensure that the code is present in the quantiles.
 *
 * \param [in]  n_quantiles     Number of quantiles
 * \param [in]  code            Morton code we are searching for
 * \param [in]  quantile_start  First Morton code in each quantile (size = \ref n_quantiles)
 *
 * \return                      Quantile associated to the given code in the codes array.
 *
 */

size_t
PDM_morton_quantile_search(size_t             n_quantiles,
                           PDM_morton_code_t  code,
                           PDM_morton_code_t *quantile_start)
{
  size_t mid_id = 0;
  size_t start_id = 0;
  size_t end_id = n_quantiles;

  /* use binary search */

  while (start_id + 1 < end_id) {
    mid_id = start_id + ((end_id -start_id) / 2);
    if (_a_gt_b(quantile_start[mid_id], code))
      end_id = mid_id;
    else
      start_id = mid_id;
  }

  /* We may have stopped short of the required value,
     or have multiple occurences of a quantile start
     (in case of empty quantiles), of which we want to
     find the find highest one */

  while (   start_id < n_quantiles - 1
         && _a_ge_b(code, quantile_start[start_id+1]))
    start_id++;

  return start_id;
}



/**
 * \brief Get the quantiles intersected by a Morton code using a binary search.
 *
 * No check is done to ensure that the code is present in the quantiles.
 *
 * \param [in]  n_quantiles     Number of quantiles
 * \param [in]  code            Morton code we are searching for
 * \param [in]  quantile_start  First Morton code in each quantile (size = \ref n_quantiles)
 * \param [out] n_intersect     Number of intersected quantiles
 * \param [out] intersect       Intersected quantiles (size = \ref n_quantiles)
 *
 */

void
PDM_morton_quantile_intersect(size_t             n_quantiles,
                              PDM_morton_code_t  code,
                              PDM_morton_code_t *quantile_start,
                              size_t            *start,
                              size_t            *end )
{
  size_t mid_id = 0;
  size_t start_id = 0;
  size_t end_id = n_quantiles;

  while (start_id + 1 < end_id) {
    mid_id = start_id + ((end_id -start_id) / 2);
    if (_a_gt_b(code, quantile_start[mid_id]))
      start_id = mid_id;
    else
      end_id = mid_id;
  }

  size_t start_id_save = start_id;
  int assigned = 0, condition = 0;
  for (size_t i = start_id_save; i < n_quantiles; i++) {

    if (i == n_quantiles-1) {
      condition = 1;
    } else {
      condition = _a_gtmin_b(quantile_start[i+1], code);
    }

    if (condition) {
      if (_a_gt_b(quantile_start[i], code)) {
        break;
      } else {
        if (!assigned) {
          start_id = i;
          end_id = i;
          assigned = 1;
        }
        end_id++;
      }
    }
  }

  *start = start_id;
  *end   = end_id;
}



/**
 * \brief Get the Morton codes intersected by a given Morton code using a binary search.
 *
 * The array of Morton codes \emph must be sorted in ascending order.
 *
 * \param [in]  n_codes       Number of codes
 * \param [in]  code          Morton code we are searching for
 * \param [in]  codes         First Morton code in each code (size = \ref n_codes)
 * \param [out] start         Id of the first intersected Morton code
 * \param [out] end           Id of the last intersected Morton code + 1
 *
 */

void
PDM_morton_list_intersect(size_t             n_quantiles,
                          PDM_morton_code_t  code,
                          PDM_morton_code_t *codes,
                          size_t            *start,
                          size_t            *end)
{
  size_t mid_id = 0;
  size_t start_id = 0;
  size_t end_id = n_quantiles;

  /* use binary search */

  while (start_id + 1 < end_id) {
    mid_id = start_id + ((end_id -start_id) / 2);
    if (_a_gt_b(code, codes[mid_id]))
      start_id = mid_id;
    else
      end_id = mid_id;
  }

  int start_id_save = start_id;

  if (PDM_morton_ancestor_is (code, codes[start_id]) ||
      PDM_morton_ancestor_is (codes[start_id], code)) {
    // quantile #start_id_save does intersect
    // Place start_id to the leftmost duplicate of quantile #start_id_save
    while (start_id > 0
           && codes[start_id_save].L    == codes[start_id-1].L
           && codes[start_id_save].X[0] == codes[start_id-1].X[0]
           && codes[start_id_save].X[1] == codes[start_id-1].X[1]
           && codes[start_id_save].X[2] == codes[start_id-1].X[2]) {
      --start_id;
    }
  } else {
    // quantile #start_id_save does NOT intersect
    // Place start_id to the rightmost duplicate of quantile #start_id_save
    start_id++;
    while (start_id < n_quantiles-1
           && codes[start_id_save].L    == codes[start_id+1].L
           && codes[start_id_save].X[0] == codes[start_id+1].X[0]
           && codes[start_id_save].X[1] == codes[start_id+1].X[1]
           && codes[start_id_save].X[2] == codes[start_id+1].X[2]) {
      ++start_id;
    }
  }

  // Sweep right
  end_id = start_id;
  while (end_id < n_quantiles
         && (PDM_morton_ancestor_is (code, codes[end_id]) ||
             PDM_morton_ancestor_is (codes[end_id], code))) {
    end_id++;
  }

  *start = start_id;
  *end   = end_id;
}



/**
 * \brief Build a global Morton encoding rank index from sorted Morton codes
 *
 * The rank_index[i] contains the first Morton code assigned to rank i.
 *
 * \param [in]  dim           Spatial dimension (1, 2 or 3)
 * \param [in]  gmax_level    Level in octree used to build the Morton encoding
 * \param [in]  n_codes       Number of Morton codes to be indexed
 * \param [in]  ordered_code  Array of Morton codes to be indexed (size = \ref n_codes)
 * \param [in]  weight        Weight associated to each Morton code (size = \ref n_codes)
 * \param [out] rank_index    Pointer to the global Morton encoding rank index (size = n_rank + 1)
 * \param [in]  comm          MPI communicator on which we build the global index
 *
 */

void
PDM_morton_ordered_build_rank_index
(
 int                      dim,
 int                      gmax_level,
 PDM_l_num_t              n_codes,
 const PDM_morton_code_t  ordered_code[],
 const int                weight[],
 PDM_morton_code_t        rank_index[],
 PDM_MPI_Comm             comm
)
{
  PDM_UNUSED(dim);
  PDM_UNUSED(gmax_level);

  PDM_g_num_t *_weight = malloc(sizeof(PDM_g_num_t) * n_codes);

  int comm_size;
  PDM_MPI_Comm_size (comm, &comm_size);

  int comm_rank;
  PDM_MPI_Comm_rank (comm, &comm_rank);

  if (n_codes > 0) {
    _weight[0] = (PDM_g_num_t) weight[0];
    for (int i = 1; i < n_codes; i++) {
      _weight[i] = (PDM_g_num_t) weight[i] + _weight[i-1];
    }
  }

//  log_trace("n_codes = %d\n", n_codes);
//  PDM_log_trace_array_int ( weight, n_codes ,"weight  : ");
//  PDM_log_trace_array_long(_weight, n_codes ,"_weight : ");


  PDM_g_num_t scan;
  if  (n_codes > 0) {
    PDM_MPI_Scan (_weight + n_codes-1, &scan, 1, PDM__PDM_MPI_G_NUM, PDM_MPI_SUM, comm);
    scan += -_weight[n_codes-1];
    for (int i = 0; i < n_codes; i++) {
      _weight[i] += scan;
    }
  }
  else {
    PDM_g_num_t __weight_0 = 0;
    PDM_MPI_Scan (&__weight_0, &scan, 1, PDM__PDM_MPI_G_NUM, PDM_MPI_SUM, comm);
  }

  PDM_g_num_t total_weight = 0;

  if (comm_rank == (comm_size - 1)) {
    if (n_codes > 0) {
      total_weight = _weight[n_codes-1];
    }
    else {
      total_weight = scan;
    }
  }

  PDM_MPI_Bcast (&total_weight, 1, PDM__PDM_MPI_G_NUM,
                 comm_size - 1, comm);

  PDM_g_num_t mean = total_weight / comm_size;
  PDM_g_num_t k    = total_weight % comm_size;

  int idx = 0;

  PDM_g_num_t *quantiles = malloc (sizeof(PDM_g_num_t) * (comm_size + 1));

  for (int i = 0; i < comm_size; i++) {
    if (i < k) {
      quantiles[i] = i * (mean + 1);
    }
    else {
      quantiles[i] = i * mean + k;
    }
  }
  quantiles[comm_size] = total_weight + 1;

  int *send_count = malloc (sizeof(int) * comm_size);
  for (int i = 0; i < comm_size; i++) {
    send_count[i] = 0;
  }

  for (int i = 0; i < n_codes; i++) {
    int i_rank = PDM_binary_search_gap_long (_weight[i], quantiles, comm_size+1);
    send_count[i_rank]++;
  }

  int *recv_count = malloc (sizeof(int) * comm_size);
  PDM_MPI_Alltoall (send_count, 1, PDM_MPI_INT,
                    recv_count, 1, PDM_MPI_INT,
                    comm);

  int *send_idx = malloc (sizeof(int) * (comm_size + 1));
  int *recv_idx = malloc (sizeof(int) * (comm_size + 1));

  send_idx[0] = 0;
  recv_idx[0] = 0;
  for (int i = 0; i < comm_size; i++) {
    send_idx[i+1] = send_idx[i] + 4 * send_count[i];
    recv_count[i] *= 4;
    recv_idx[i+1] = recv_idx[i] + recv_count[i];
    send_count[i] = 0;
  }

  PDM_morton_int_t *send_data = malloc(sizeof(PDM_morton_int_t) * send_idx[comm_size]);
  PDM_morton_int_t *recv_data = malloc(sizeof(PDM_morton_int_t) * recv_idx[comm_size]);


  for (int i = 0; i < n_codes; i++) {
    int i_rank = PDM_binary_search_gap_long (_weight[i], quantiles, comm_size+1);
    send_data[send_idx[i_rank]+send_count[i_rank]++] = ordered_code[i].L;
    send_data[send_idx[i_rank]+send_count[i_rank]++] = ordered_code[i].X[0];
    send_data[send_idx[i_rank]+send_count[i_rank]++] = ordered_code[i].X[1];
    send_data[send_idx[i_rank]+send_count[i_rank]++] = ordered_code[i].X[2];
  }
  /*for (int i = 0; i < send_idx[comm_size]/4; i++) {
    log_trace("send_data[%d] : L = %zu, X = %zu %zu %zu\n",
              i,
              send_data[4*i],
              send_data[4*i+1],
              send_data[4*i+2],
              send_data[4*i+3]);
              }*/

  free (quantiles);
  free (_weight);

  PDM_MPI_Alltoallv(send_data, send_count, send_idx, PDM_MPI_UNSIGNED,
                    recv_data, recv_count, recv_idx, PDM_MPI_UNSIGNED,
                    comm);

  int n_recv_codes = recv_idx[comm_size] / 4;
  /*for (int i = 0; i < n_recv_codes; i++) {
    log_trace("recv_data[%d] : L = %zu, X = %zu %zu %zu\n",
              i,
              recv_data[4*i],
              recv_data[4*i+1],
              recv_data[4*i+2],
              recv_data[4*i+3]);
              }*/

  free (send_data);
  free (send_count);
  free (send_idx);
  free (recv_count);
  free (recv_idx);

  PDM_morton_code_t min_code;
  min_code.L = 31u;
  min_code.X[0] = (1u << min_code.L) - 1u;
  min_code.X[1] = (1u << min_code.L) - 1u;
  min_code.X[2] = (1u << min_code.L) - 1u;

  idx = 0;

  PDM_morton_int_t send_min_code[4] = {min_code.L, min_code.X[0], min_code.X[1], min_code.X[2]};//{31u, 0, 0, 0};

  for (int i = 0; i < n_recv_codes; i++) {
    PDM_morton_code_t tmp_code;
    tmp_code.L = recv_data[idx++];
    tmp_code.X[0] = recv_data[idx++];
    tmp_code.X[1] = recv_data[idx++];
    tmp_code.X[2] = recv_data[idx++];

    /*log_trace("i = %d, tmp_code : L = %zu, X = %zu %zu %zu, min_code : L = %zu, X = %zu %zu %zu\n",
      i,
      tmp_code.L,
      tmp_code.X[0],
      tmp_code.X[1],
      tmp_code.X[2],
      min_code.L,
      min_code.X[0],
      min_code.X[1],
      min_code.X[2]);
      log_trace("min gt tmp? %d\n", _a_gt_b( min_code, tmp_code));
      log_trace("min ge tmp? %d\n", _a_ge_b( min_code, tmp_code));*/

    //if (_a_gt_b( min_code, tmp_code)) {
    if (_a_ge_b( min_code, tmp_code)) {
      min_code.L = tmp_code.L;
      min_code.X[0] = tmp_code.X[0];
      min_code.X[1] = tmp_code.X[1];
      min_code.X[2] = tmp_code.X[2];
      send_min_code[0] = tmp_code.L;
      send_min_code[1] = tmp_code.X[0];
      send_min_code[2] = tmp_code.X[1];
      send_min_code[3] = tmp_code.X[2];
    }
  }

  /*log_trace("send_min_code : L = %zu, X = %zu %zu %zu\n",
            send_min_code[0],
            send_min_code[1],
            send_min_code[2],
            send_min_code[3]);*/

  free (recv_data);

  PDM_morton_int_t *buff_min_codes = malloc(sizeof(PDM_morton_int_t) * 4 * comm_size);

  int *n_nodes = malloc(sizeof(int) * comm_size);

  PDM_MPI_Allgather (&n_recv_codes, 1, PDM_MPI_INT,
                     n_nodes, 1, PDM_MPI_INT,
                     comm);
  PDM_MPI_Allgather (send_min_code, 4, PDM_MPI_UNSIGNED,
                     buff_min_codes, 4, PDM_MPI_UNSIGNED,
                     comm);
  /*for (int i = 0; i < comm_size; i++) {
    log_trace("buff_min_codes :: rank %d : L = %zu, X = %zu %zu %zu\n",
              i,
              buff_min_codes[4*i],
              buff_min_codes[4*i+1],
              buff_min_codes[4*i+2],
              buff_min_codes[4*i+3]);
              }*/

  idx = 0;
  for (int i = 0; i < comm_size; i++) {
    //log_trace("rank %d, n_nodes = %d, idx = %d\n", i, n_nodes[i], idx);
    if (n_nodes[i] > 0) {
      rank_index[i].L = buff_min_codes[idx++];
      rank_index[i].X[0] = buff_min_codes[idx++];
      rank_index[i].X[1] = buff_min_codes[idx++];
      rank_index[i].X[2] = buff_min_codes[idx++];
    }
    else {
      idx += 4;
    }
  }

  /*for (int i = 0; i < comm_size; i++) {
    log_trace("rank_index (1) [%d] : L = %zu, X = %zu %zu %zu\n",
              i,
              rank_index[i].L,
              rank_index[i].X[0],
              rank_index[i].X[1],
              rank_index[i].X[2]);
              }*/

  for (int i = comm_size - 1; i >= 0; i--) {
    if ((n_nodes[i] == 0) && (i != comm_size - 1)) {
      rank_index[i].L = rank_index[i+1].L;
      rank_index[i].X[0] = rank_index[i+1].X[0];
      rank_index[i].X[1] = rank_index[i+1].X[1];
      rank_index[i].X[2] = rank_index[i+1].X[2];
    }
  }

  rank_index[comm_size].L = 31u;
  rank_index[comm_size].X[0] = (1u << 31u) - 1u;
  rank_index[comm_size].X[1] = (1u << 31u) - 1u;
  rank_index[comm_size].X[2] = (1u << 31u) - 1u;

  for (int i = 0; i < comm_size; i++) {
    assert(_a_ge_b(rank_index[i+1], rank_index[i]));
  }

  free (buff_min_codes);
  free (n_nodes);

}



/**
 * \brief Build a global Morton encoding rank index
 *
 * The rank_index[i] contains the first Morton code assigned to rank i.
 *
 * \param [in]  dim           Spatial dimension (1, 2 or 3)
 * \param [in]  gmax_level    Level in octree used to build the Morton encoding
 * \param [in]  n_codes       Number of Morton codes to be indexed
 * \param [in]  ordered_code  Array of Morton codes to be indexed (size = \ref n_codes)
 * \param [in]  weight        Weight associated to each Morton code (size = \ref n_codes)
 * \param [in]  order         Ordering array (size = \ref n_codes)
 * \param [out] rank_index    Pointer to the global Morton encoding rank index (size = n_rank + 1)
 * \param [in]  comm          MPI communicator on which we build the global index
 *
 */

double
PDM_morton_build_rank_index(int                     dim,
                            int                     gmax_level,
                            PDM_l_num_t             n_codes,
                            const PDM_morton_code_t code[],
                            const double            weight[],
                            const int               order[],
                            PDM_morton_code_t       rank_index[],
                            PDM_MPI_Comm            comm)
{
  int  i, id, rank_id, n_ranks, n_samples;
  double  best_fit;

  double  *sampling = NULL;

  const int  sampling_factor = _sampling_factors[dim];

  /* Allocations and Initialization */

  PDM_MPI_Comm_size(comm, &n_ranks);

  n_samples = sampling_factor * n_ranks;

  sampling = (double *) malloc((n_samples + 1) * sizeof(double));

  for (i = 0; i < n_samples + 1; i++)
    sampling[i] = 0.0;

  best_fit = _bucket_sampling(dim,
                              n_ranks,
                              gmax_level,
                              n_codes,
                              code,
                              weight,
                              order,
                              &sampling,
                              comm);

  /* Define Morton index */

  for (rank_id = 0; rank_id < n_ranks + 1; rank_id++) {

    id = rank_id * sampling_factor;
    rank_index[rank_id] = _double_to_code(dim, sampling[id], gmax_level);

  }

#if 0 && defined(dbg_enabled) && !defined(Ndbg_enabled)
  { /* Dump Morton index and associated sampling on rank 0 */
    PDM_printf("\nMorton rank index:\n\n");
    for (rank_id = 0; rank_id < n_ranks + 1; rank_id++) {
      id = sampling_factor * rank_id;
      PDM_printf("rank: %5d (sampling: %7.4g)- ", rank_id, sampling[id]);
      PDM_morton_dump(dim, rank_index[rank_id]);

    }
    PDM_printf("\n");
    fflush(stdout);
  }
#endif

  /* Free memory */

  free(sampling);

  return best_fit;
}



/**
 * \brief Dump a Morton to standard output or to a file
 *
 * \param [in]  dim           Spatial dimension (2 or 3)
 * \param [in]  code          Morton code to dump
 *
 */

void
PDM_morton_dump(int                dim,
                PDM_morton_code_t  code)
{
  int  i;
  double  coord[3];

  const unsigned long   n = 1u << code.L;
  const double  stride = 1/(double)n;

  for (i = 0; i < dim; i++)
    coord[i] = stride * code.X[i];

  if (dim == 3)
    //PDM_printf("Morton Code:\n"
    log_trace("Morton Code:\n"
               "L =  %3u [X, Y, Z] - [%5u %5u %5u]"
               "[%6.5lf %6.5lf %6.5lf]\n",
               code.L, code.X[0], code.X[1], code.X[2],
               coord[0], coord[1], coord[2]);

  else if (dim == 2)
    PDM_printf("Morton Code\n"
               "L =  %3u [X, Y] - [%5u %5u] [%6.5lf %6.5lf]\n",
               code.L, code.X[0], code.X[1], coord[0], coord[1]);

  fflush(stdout);
}



/**
 * \brief Test if Morton code 'a' is an ancestor of Morton code 'b'
 *
 * \param [in]   a     First Morton code
 * \param [in]   b     Second Morton code
 *
 * \return             True or false
 *
 */

_Bool
PDM_morton_ancestor_is (PDM_morton_code_t  a,
                        PDM_morton_code_t  b)
{
  _Bool status = 0;

  if (a.L <= b.L) {
    PDM_morton_assign_level (&b,
                             a.L);

    assert (a.L == b.L);

    status = a.X[0] == b.X[0] && a.X[1] == b.X[1] && a.X[2] == b.X[2];
  }

  return status;
}







/**
 * \brief Intersect a box with an array of sorted Morton codes
 *
 * A recursive top-down approach is used, starting from the deepest common
 * ancestor of all the Morton codes in the array.
 *
 * \param [in]     dim           Spatial dimension
 * \param [in]     node          Morton code of current node
 * \param [in]     box_min       Morton code of the box's lower corner
 * \param [in]     box_max       Morton code of the box's upper corner
 * \param [in]     nodes         Array of Morton codes
 * \param [in]     n_points      Number of points stored inside each node
 * \param [in]     start         Id of the first descendant of the current node
 * \param [in]     end           Id of the last descendant of the current node + 1
 * \param [in,out] n_intersect   Number of intersected nodes
 * \param [in,out] intersect     Intersected nodes
 *
 */

static const size_t N_BRUTE_FORCE = 10;
void
PDM_morton_intersect_box
(
 const int                dim,
 const PDM_morton_code_t  node,
 const PDM_morton_code_t  box_min,
 const PDM_morton_code_t  box_max,
 const PDM_morton_code_t  nodes[],
 int                     *n_points,
 const size_t             start,
 const size_t             end,
 size_t                  *n_intersect,
 int                     *intersect
 )
{
  int dbg_enabled = 0;
  int inside;

  if (dbg_enabled) {
    printf("node: L = %u, X = %u %u %u, start = %zu, end = %zu\n",
           node.L, node.X[0], node.X[1], node.X[2], start, end);
  }

  /* If current range contains few octants, go brute force */
  if (end - start < N_BRUTE_FORCE) {

    if (n_points == NULL) {
      for (size_t i = start; i < end; i++) {
        if (_intersect_node_box (dim,
                                 nodes[i],
                                 box_min,
                                 box_max,
                                 &inside)) {
          intersect[(*n_intersect)++] = i;
        }
      }
    }

    else {
      for (size_t i = start; i < end; i++) {
        if (n_points[i] > 0) {
          if (_intersect_node_box (dim,
                                   nodes[i],
                                   box_min,
                                   box_max,
                                   &inside)) {
            intersect[(*n_intersect)++] = i;
          }
        }
      }
    }
    return;

  }

  else {

    if (_intersect_node_box (dim,
                             node,
                             box_min,
                             box_max,
                             &inside)) {

      if (inside) {
        /* Every descendant must intersect the box */
        if (n_points == NULL) {
          for (size_t i = start; i < end; i++) {
            intersect[(*n_intersect)++] = i;
          }
        }
        else {
          for (size_t i = start; i < end; i++) {
            if (n_points[i] > 0) {
              intersect[(*n_intersect)++] = i;
            }
          }
        }
      }

      else {
        /* Some descendants may intersect the box */
        const size_t n_children = 1 << dim;
        PDM_morton_code_t children[8];
        PDM_morton_get_children (dim,
                                 node,
                                 children);

        size_t new_start, new_end;
        size_t prev_end = start;
        for (size_t ichild = 0; ichild < n_children; ichild++) {
          if (dbg_enabled) {
            printf("  child: L = %u, X = %u %u %u\n",
                   children[ichild].L,
                   children[ichild].X[0], children[ichild].X[1], children[ichild].X[2]);
          }

          /* get start and end of range in list of nodes covered by current child */
          /* new_start <-- first descendant of child in list */
          new_start = prev_end; // end of previous child's range
          // linear search
          while (new_start < end) {
            if (PDM_morton_ancestor_is (children[ichild], nodes[new_start])) {
              break;
            } else if (_a_gt_b(nodes[new_start], children[ichild])) {
              /* all the following nodes are clearly not descendants of current child */
              new_start = end+1;
              break;
            }
            new_start++;
          }

          if (dbg_enabled) {
            printf("   new_start = %zu\n", new_start);
          }

          if (new_start > end) {
            /* no need to go further for that child
               because it has no descendants in the node list */
            continue;
          }

          /* new_end <-- next of last descendant of child in list */
          size_t l = new_start;
          new_end = end;
          while (new_end > l + 1) {
            size_t m = l + (new_end - l) / 2;
            if (PDM_morton_ancestor_is (children[ichild], nodes[m])) {
              l = m;
            } else {
              new_end = m;
            }
          }

          prev_end = new_end;
          if (dbg_enabled) {
            printf("   new_end   = %zu\n", new_end);
          }


          /* Carry on recursion */
          PDM_morton_intersect_box (dim,
                                    children[ichild],
                                    box_min,
                                    box_max,
                                    nodes,
                                    n_points,
                                    new_start,
                                    new_end,
                                    n_intersect,
                                    intersect);
        }
      }
    }
  }
}








/**
 * \brief Get the closest node (described by Morton codes) of a given point.
 *
 * A recursive top-down approach is used, starting from the deepest common
 * ancestor of all the Morton codes in the array.
 * (NOT USED ANYMORE)
 *
 * \param [in]     dim            Spatial dimension
 * \param [in]     node           Morton code of current node
 * \param [in]     nodes          Array of Morton codes
 * \param [in]     point          Coordinates of the point
 * \param [in]     d              Normalization vector (dilatation component)
 * \param [in]     start          Id of the first descendant of the current node
 * \param [in]     end            Id of the last descendant of the current node + 1
 * \param [in,out] closest_node   Id of the closest node
 * \param [in,out] closest_dist2  Squared distance between the point and the closest node
 *
 */

void
PDM_morton_closest_node
(
 const int                dim,
 const PDM_morton_code_t  node,
 const PDM_morton_code_t  nodes[],
 const double             point[],
 const double             d[],
 const size_t             start,
 const size_t             end,
 int                     *closest_node,
 double                  *closest_dist2
 )
{
  /* If current range contains few octants, go brute force */
  if (end - start < N_BRUTE_FORCE) {

    for (size_t i = start; i < end; i++) {
      double dist2 = _min_dist2 (dim,
                                 nodes[i],
                                 point,
                                 d);
      if (dist2 < *closest_dist2) {
        *closest_dist2 = dist2;
        *closest_node = i;
      }
    }
    return;
  }

  else {
    double dist2 = _min_dist2 (dim,
                               node,
                               point,
                               d);
    if (dist2 < *closest_dist2) {
      /* Inspect children of current node */
      const size_t n_children = 1 << dim;
      PDM_morton_code_t children[8];
      PDM_morton_get_children (dim,
                               node,
                               children);

      size_t new_start, new_end;
      size_t prev_end = start;
      for (size_t ichild = 0; ichild < n_children; ichild++) {
        /* get start and end of range in list of nodes covered by current child */
          /* new_start <-- first descendant of child in list */
          new_start = prev_end; // end of previous child's range
          // linear search
          while (new_start < end) {
            if (PDM_morton_ancestor_is (children[ichild], nodes[new_start])) {
              break;
            } else if (_a_gt_b(nodes[new_start], children[ichild])) {
              /* all the following nodes are clearly not descendants of current child */
              new_start = end+1;
              break;
            }
            new_start++;
          }

          if (new_start > end) {
            /* no need to go further for that child
               because it has no descendants in the node list */
            continue;
          }

          /* new_end <-- next of last descendant of child in list */
          size_t l = new_start;
          new_end = end;
          while (new_end > l + 1) {
            size_t m = l + (new_end - l) / 2;
            if (PDM_morton_ancestor_is (children[ichild], nodes[m])) {
              l = m;
            } else {
              new_end = m;
            }
          }

          prev_end = new_end;

          /* Carry on recursion */
          PDM_morton_closest_node (dim,
                                   children[ichild],
                                   nodes,
                                   point,
                                   d,
                                   new_start,
                                   new_end,
                                   closest_node,
                                   closest_dist2);

      } // End loop children
    }
  }
}


/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
}
#endif /* __cplusplus */
