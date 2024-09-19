/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <assert.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "pdm.h"
#include "pdm_priv.h"
#include "pdm_config.h"
#include "pdm_printf.h"
#include "pdm_logging.h"
#include "pdm_error.h"
#include "pdm_sort.h"
#include "pdm_binary_search.h"
#include "pdm_array.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "pdm_point_tree_seq.h"
#include "pdm_point_tree_seq_priv.h"

/*----------------------------------------------------------------------------*/


/*============================================================================
 * Local macro definitions
 *============================================================================*/

#if !defined(HUGE_VAL)
#define HUGE_VAL 1.0e+30
#endif

/*============================================================================
 * Local structure definitions
 *============================================================================*/

/*============================================================================
 * Global variable
 *============================================================================*/

static const double _eps_default = 1.e-12;
static const int    dbg_ptree    = 0;

/*=============================================================================
 * Private function definitions
 *============================================================================*/


/*=============================================================================
 * Private function definitions
 *============================================================================*/


/**
 * \brief   Evaluate a distribution array.
 *
 * \param [in]  n_ranges     Number of ranges in the distribution
 * \param [in]  distribution Number of elements associated to each range of the distribution
 * \param [in]  optim        Optimal count in each range
 *
 * \return  a fit associated to the distribution. If fit = 0, distribution is perfect.
 *
 */
static double
_evaluate_distribution(int          n_ranges,
                       int         *distribution,
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

#if 0 && defined(DEBUG) && !defined(NDEBUG)
  if (cs_glob_rank_id <= 0)
    PDM_printf( "<DISTRIBUTION EVALUATION> optim: %g, fit: %g\n",
               optim, fit);
#endif

  return  fit;
}




static double
_median_point
(
 int     split_direction,
 double  extents_min,
 double  extents_max,
 int     point_range[2],
 double *pts_coord
 )
{
  double mid = 0.5*(extents_min + extents_max);

  int n_pts = point_range[1] - point_range[0];
  if (n_pts > 2) {
    double *x = malloc(sizeof(double) * n_pts);

    for (int j = 0; j < n_pts; j++) {
      int i = point_range[0] + j;

      x[j] = pts_coord[3*i + split_direction];
    }

    PDM_sort_double(x, NULL, n_pts);

    int h = n_pts / 2;
    // log_trace("n_pts = %d, h = %d\n", n_pts, h);
    mid = 0.5*(x[h] + x[h+1]);
    free(x);
  }

  return mid;
}


static void
_define_rank_distrib
(
       int     split_direction,
       int     point_range[2],
 const double *pts_coord,
       int     n_sample,
       double *sampling,
       double *cfreq,
       int     distrib[2]
 )
{
  double i_npts = 1. / (double) (point_range[1] - point_range[0]);

  int l_distrib[n_sample];
  for (int i = 0; i < n_sample; i++) {
    l_distrib[i] = 0;
  }


  for (int i = point_range[0]; i < point_range[1]; i++) {
    double x = pts_coord[3*i + split_direction];

    int isample = PDM_binary_search_gap_double(x,
                                               sampling,
                                               n_sample+1);

    // if (isample >= n_sample || isample < 0) {
    //   PDM_log_trace_array_double(sampling, n_sample+1, "sampling : ");
    //   log_trace("!!! x = %f, isample = %d / %d\n", x, isample, n_sample);
    // }

    l_distrib[isample]++;
  }

  // PDM_log_trace_array_int(l_distrib, n_sample, "l_distrib : ");


  /* Define the cumulative frequency related to g_distribution */
  cfreq[0] = 0.;
  for (int id = 0; id < n_sample; id++) {
    cfreq[id+1] = cfreq[id] + l_distrib[id] * i_npts;
  }
  cfreq[n_sample] = 1.0;
  // PDM_log_trace_array_double(cfreq, n_sample+1, "cfreq : ");

  distrib[0] = 0;
  for (int id = 0; id < n_sample/2; id++) {
    distrib[0] += l_distrib[id];
  }
  distrib[1] = (point_range[1] - point_range[0]) - distrib[0];

  // PDM_log_trace_array_int(distrib, 2, "distrib : ");

}


static void
_update_sampling(int     n_sample,
                 double  c_freq[],
                 double *sampling[])
{
  int  i, j, next_id;
  double  target_freq, f_high, f_low, delta;
  double  s_low, s_high;

  // double new_sampling[n_sample+1];
  double *new_sampling = malloc(sizeof(double) * (n_sample+1));
  double *_sampling = *sampling;


  const double unit = 1/(double)n_sample;

  /* Compute new_sampling */
  new_sampling[0] = _sampling[0];
  next_id = 1;

  for (i = 0; i < n_sample; i++) {

    target_freq = (i+1)*unit;

    /* Find the next id such as c_freq[next_id] >= target_freq */

    for (j = next_id; j < n_sample + 1; j++) {
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

  } /* End of loop on samples */


  new_sampling[n_sample] = _sampling[n_sample];

  free(_sampling);

  /* Return pointers */
  *sampling = new_sampling;
}


static double
_approx_median_point
(
 const int     n_sample,
       int     split_direction,
       double  extents_min,
       double  extents_max,
       int     point_range[2],
 const double *pts_coord
 )
{
  double  fit, best_fit, optim;

  int n_pts = point_range[1] - point_range[0];

  // double sampling[n_sample+1];
  double *sampling = malloc(sizeof(double) * (n_sample+1));

   /* Define a naive sampling (uniform distribution) */
  double step = (extents_max - extents_min) / (double) n_sample;
  for (int i = 0; i <= n_sample; i++) {
    sampling[i] = extents_min + i*step;
  }
  sampling[n_sample] += 1e-3;

  // PDM_log_trace_array_double(sampling, n_sample+1, "sampling : ");


  int distrib[2];
  double cfreq[n_sample+1];
  _define_rank_distrib(split_direction,
                       point_range,
                       pts_coord,
                       n_sample,
                       sampling,
                       cfreq,
                       distrib);

  optim = 0.5*n_pts;

  /* Initialize best choice */

  fit = _evaluate_distribution(2, distrib, optim);
  best_fit = fit;

  double best_sampling[n_sample+1];
  for (int i = 0; i < (n_sample + 1); i++) {
    best_sampling[i] = sampling[i];
  }

  /* Loop to get a better sampling array */

  // log_trace(">> loop\n");
  for (int n_iters = 0; (n_iters < 5 && fit > 0.10); n_iters++) {

    _update_sampling(n_sample, cfreq, &sampling);

    /* Compute the new distribution associated to the new sampling */

    _define_rank_distrib(split_direction,
                         point_range,
                         pts_coord,
                         n_sample,
                         sampling,
                         cfreq,
                         distrib);

    fit = _evaluate_distribution(2, distrib, optim);
    // log_trace("n_iters = %d, fit = %f\n", n_iters, fit);
    // PDM_log_trace_array_double(sampling, n_sample+1, "sampling : ");

    /* Save the best sampling array and its fit */

    if (fit < best_fit) {

      best_fit = fit;
      for (int i = 0; i < (n_sample + 1); i++){
        best_sampling[i] = sampling[i];
      }
    }

  } /* End of while */
  free(sampling);

  // PDM_log_trace_array_double(best_sampling, n_sample+1, "best_sampling : ");


  double mid = best_sampling[n_sample/2];
  // log_trace("mid = %f\n", mid);

  return mid;
}



static void
_build_point_tree_seq_leaves
(
 const int                         ancestor_id,
 const PDM_point_tree_seq_child_t  location_in_ancestor,
 const int                         depth,
 const double                      extents[],
       PDM_point_tree_seq_t       *ptree,
       int                         point_range[2],
       int                        *new_to_old
 )
{
  _l_nodes_t *nodes = ptree->nodes;

  if (dbg_ptree) {
    log_trace("\nnode_id = %d\n",
              ptree->n_nodes);
    log_trace("ancestor_id = %d, location_in_ancestor = %d, depth = %d, point_range = %d/%d\n",
              ancestor_id, (int) location_in_ancestor, depth, point_range[0], point_range[1]);
    log_trace("extents = %f %f %f  %f %f %f\n",
              extents[0], extents[1], extents[2], extents[3], extents[4], extents[5]);
  }

  int n_children = PDM_point_tree_n_children_get(ptree);

  /* Resize point_tree if necessary */
  int _n_nodes = ptree->n_nodes;
  int tmp_size = ptree->n_nodes;

  if (ptree->n_nodes >= ptree->n_nodes_max) {
    if (ptree->n_nodes == 0) {
      ptree->n_nodes     = 1;
      ptree->n_nodes_max = 8;
    }
    ptree->n_nodes_max *= 2;

    nodes->ancestor_id = realloc(nodes->ancestor_id, sizeof(int   ) * ptree->n_nodes_max);
    nodes->is_leaf     = realloc(nodes->is_leaf,     sizeof(int   ) * ptree->n_nodes_max);
    nodes->depth       = realloc(nodes->depth,       sizeof(int   ) * ptree->n_nodes_max);
    nodes->children_id = realloc(nodes->children_id, sizeof(int   ) * ptree->n_nodes_max * n_children);
    nodes->range       = realloc(nodes->range,       sizeof(int   ) * ptree->n_nodes_max * 2);
    nodes->idx         = realloc(nodes->idx,         sizeof(int   ) * ptree->n_nodes_max * (n_children+1));
    nodes->n_points    = realloc(nodes->n_points,    sizeof(int   ) * ptree->n_nodes_max);
    nodes->extents     = realloc(nodes->extents,     sizeof(double) * ptree->n_nodes_max * 6);
    nodes->location_in_ancestor = realloc(nodes->location_in_ancestor, sizeof(PDM_point_tree_seq_child_t) * ptree->n_nodes_max);
  }


  /* Number of points */
  int _n_points = point_range[1] - point_range[0];

  int idx[9];
  int child_id[8] = {-1, -1, -1, -1, -1, -1, -1, -1};
  int octant_mask[3] = {4, 2, 1}; /* pow(2, 2), pow(2, 1), pow(2,0) */

  int is_leaf = 1;
  if (depth < ptree->depth_max && _n_points > ptree->points_in_leaf_max) {

    /* Choose split direction */
    double max_range = -1.;
    int split_direction = -1;

    if (ptree->tree_type == PDM_DOCTREE_LOCAL_TREE_KDTREE) {
      for (int direction = 0; direction < 3; direction++) {
        double range = extents[3+direction] - extents[direction];
        if (range > max_range) {
          max_range        = range;
          split_direction = direction;
        }
      }
      if (dbg_ptree) {
        log_trace("split_direction = %d\n", split_direction);
      }
    }

    /* Choose split point */
    double mid[3];
    if (ptree->tree_type == PDM_DOCTREE_LOCAL_TREE_KDTREE) {
      if (1) {
        mid[0] = _approx_median_point(6,
                                      split_direction,
                                      extents[split_direction],
                                      extents[3+split_direction],
                                      point_range,
                                      ptree->_pts_coord);
        // mid[0] = 0.5*(extents[split_direction] + extents[3+split_direction]);
      }
      else {
        mid[0] = _median_point(split_direction,
                               extents[split_direction],
                               extents[3+split_direction],
                               point_range,
                               ptree->_pts_coord);
      }
      if (dbg_ptree) {
        log_trace("mid = %f\n", mid[0]);
      }
    }
    else {
      for (int i = 0; i < 3; i++) {
        mid[i] = 0.5*(extents[i] + extents[3+i]);
      }
      if (dbg_ptree) {
        log_trace("mid = %f %f %f\n", mid[0], mid[1], mid[2]);
      }
    }


    /* Count and reorder points in each child node */
    int count[8] = {0, 0, 0, 0, 0, 0, 0, 0};
    for (int i = point_range[0]; i < point_range[1]; i++) {
      int ichild = 0;
      if (ptree->tree_type == PDM_DOCTREE_LOCAL_TREE_KDTREE) {
        ichild = (ptree->_pts_coord[3*i + split_direction] > mid[0]);
      }
      else {
        for (int j = 0; j < 3; j++) {
          if (ptree->_pts_coord[3*i + j] > mid[j]) {
            ichild += octant_mask[j];
          }
        }
      }
      count[ichild]++;
    }

    if (dbg_ptree) {
      PDM_log_trace_array_int(count, n_children, "count : ");
    }


    /* Build index */
    idx[0] = 0;
    for (int ichild = 0; ichild < n_children; ichild++) {
      idx[ichild+1] = idx[ichild] + count[ichild];
      count[ichild] = 0;
    }

    for (int i = point_range[0]; i < point_range[1]; i++) {
      int ichild = 0;
      if (ptree->tree_type == PDM_DOCTREE_LOCAL_TREE_KDTREE) {
        ichild = (ptree->_pts_coord[3*i + split_direction] > mid[0]);
      }
      else {
        for (int j = 0; j < 3; j++) {
          if (ptree->_pts_coord[3*i + j] > mid[j]) {
            ichild += octant_mask[j];
          }
        }
      }

      int old = ptree->new_to_old[i];
      int pos = point_range[0] + idx[ichild] + count[ichild];

      new_to_old[pos]        = old;
      ptree->old_to_new[old] = pos;

      count[ichild]++;
    }

    for (int i = point_range[0]; i < point_range[1]; i++) {
      ptree->new_to_old[i] = new_to_old[i];
    }
    if (dbg_ptree) {
      PDM_log_trace_array_int(ptree->new_to_old + point_range[0],
                              _n_points,
                              "new_to_old: ");
    }

    /* Reorder points */
    for (int i = point_range[0]; i < point_range[1]; i++) {
      memcpy(ptree->_pts_coord + 3*i,
             ptree->pts_coord  + 3*ptree->new_to_old[i],
             sizeof(double) * 3);
    }


    for (int i = 0; i <= n_children; i++) {
      idx[i] += point_range[0];
    }

    if (dbg_ptree) {
      PDM_log_trace_array_int(idx, n_children+1, "idx : ");
    }

    /* Build leaves recursively */
    double sub_extents[6];
    for (int ichild = 0; ichild < n_children; ichild++) {

      if (idx[ichild+1] <= idx[ichild]) {
        continue;
      }

      tmp_size++;

      child_id[ichild] = tmp_size;
      is_leaf = 0;

      memcpy(sub_extents, extents, sizeof(double) * 6);
      if (ptree->tree_type == PDM_DOCTREE_LOCAL_TREE_KDTREE) {
        if (ichild == 0) {
          sub_extents[3+split_direction] = mid[0];
        }
        else {
          sub_extents[split_direction] = mid[0];
        }
      }
      else if (ptree->tree_type == PDM_DOCTREE_LOCAL_TREE_OCTREE) {
        // log_trace("child #%d\n", ichild);
        if (ichild%2 == 0) {
          sub_extents[5] = mid[2];
        }
        else {
          sub_extents[2] = mid[2];
        }

        if (ichild%4 < 2) {
          sub_extents[4] = mid[1];
        }
        else {
          sub_extents[1] = mid[1];
        }

        if (ichild < 4) {
          sub_extents[3] = mid[0];
        }
        else {
          sub_extents[0] = mid[0];
        }
      }


      // int _is_leaf = (idx[ichild+1] - idx[ichild]) <=
      //                 ptree->points_in_leaf_max;
      if (dbg_ptree) {//_is_leaf) {
        log_trace("child %d, id %d, loose extents = %f %f %f  %f %f %f\n",
                  ichild, child_id[ichild],
                  sub_extents[0], sub_extents[1], sub_extents[2],
                  sub_extents[3], sub_extents[4], sub_extents[5]);
      }

      /* Tight extents (fit contained points) */
      if (1) {
        for (int j = 0; j < 3; j++) {
          sub_extents[j  ] =  HUGE_VAL;
          sub_extents[j+3] = -HUGE_VAL;
        }
        for (int ipt = idx[ichild]; ipt < idx[ichild+1]; ipt++) {
          for (int j = 0; j < 3; j++) {
            double x = ptree->_pts_coord[3*ipt+j];
            sub_extents[j  ] = PDM_MIN(sub_extents[j  ], x);
            sub_extents[j+3] = PDM_MAX(sub_extents[j+3], x);
          }
        }
      }

      /* Inflate slightly */
      for (int j = 0; j < 3; j++) {
        // if (sub_extents[j+3] < sub_extents[j] + _eps_default) {
        sub_extents[j  ] -= 0.5*_eps_default;
        sub_extents[j+3] += 0.5*_eps_default;
        // }
      }

      if (dbg_ptree) {//_is_leaf) {
        log_trace("child %d, id %d, tight extents = %f %f %f  %f %f %f\n",
                  ichild, child_id[ichild],
                  sub_extents[0], sub_extents[1], sub_extents[2],
                  sub_extents[3], sub_extents[4], sub_extents[5]);
      }

      ptree->n_nodes = tmp_size;

      _build_point_tree_seq_leaves(_n_nodes,
                                   (PDM_point_tree_seq_child_t) ichild,
                                   depth + 1,
                                   sub_extents,
                                   ptree,
                                   idx + ichild,
                                   new_to_old);

      tmp_size = ptree->n_nodes;
    }

  }

  /* Finalize node */
  for (int i = 0; i < 2; i++) {
    nodes->range[2*_n_nodes + i] = point_range[i];
  }

  for (int i = 0; i <= n_children; i++) {
    nodes->idx[(n_children+1)*_n_nodes + i] = idx[i];
  }

  for (int i = 0; i < 6; i++) {
    nodes->extents[6*_n_nodes + i] = extents[i];
  }

  for (int i = 0; i < n_children; i++) {
    nodes->children_id[n_children*_n_nodes + i] = child_id[i];
  }

  nodes->is_leaf[_n_nodes] = is_leaf;

  nodes->ancestor_id[_n_nodes] = ancestor_id;
  nodes->depth[_n_nodes]       = depth;

  nodes->n_points[_n_nodes] = _n_points;
  nodes->location_in_ancestor[_n_nodes] = location_in_ancestor;
}


inline static int
_intersect_box_box
(
 const int              dim,
 const double *restrict box_extents_a,
 const double *restrict box_extents_b
 )
{
  for (int i = 0; i < dim; i++) {
    if (box_extents_a[i] > box_extents_b[i+dim] || box_extents_b[i] > box_extents_a[i+dim]) {
      return 0;
    }
  }

  return 1;
}

static void
_build_point_tree_seq_leaves_from_boxes
(
 const int                         ancestor_id,
 const PDM_point_tree_seq_child_t  location_in_ancestor,
 const int                         depth,
 const double                      extents[],
 const double                     *box_extents,
 const int                         curr_n_box,
       int                        *curr_box_ids,
       PDM_point_tree_seq_t       *ptree,
       int                         point_range[2],
       int                        *new_to_old
 )
{
  _l_nodes_t *nodes = ptree->nodes;

  if (dbg_ptree) {
    log_trace("\nnode_id = %d\n",
              ptree->n_nodes);
    log_trace("ancestor_id = %d, location_in_ancestor = %d, depth = %d, point_range = %d/%d\n",
              ancestor_id, (int) location_in_ancestor, depth, point_range[0], point_range[1]);
    log_trace("extents = %f %f %f  %f %f %f\n",
              extents[0], extents[1], extents[2], extents[3], extents[4], extents[5]);
  }

  int n_children = PDM_point_tree_n_children_get(ptree);

  /* Resize point_tree if necessary */
  int _n_nodes = ptree->n_nodes;
  int tmp_size = ptree->n_nodes;

  if (ptree->n_nodes >= ptree->n_nodes_max) {
    if (ptree->n_nodes == 0) {
      ptree->n_nodes     = 1;
      ptree->n_nodes_max = 8;
    }
    ptree->n_nodes_max *= 2;

    nodes->ancestor_id = realloc(nodes->ancestor_id, sizeof(int   ) * ptree->n_nodes_max);
    nodes->is_leaf     = realloc(nodes->is_leaf,     sizeof(int   ) * ptree->n_nodes_max);
    nodes->depth       = realloc(nodes->depth,       sizeof(int   ) * ptree->n_nodes_max);
    nodes->children_id = realloc(nodes->children_id, sizeof(int   ) * ptree->n_nodes_max * n_children);
    nodes->range       = realloc(nodes->range,       sizeof(int   ) * ptree->n_nodes_max * 2);
    nodes->idx         = realloc(nodes->idx,         sizeof(int   ) * ptree->n_nodes_max * (n_children+1));
    nodes->n_points    = realloc(nodes->n_points,    sizeof(int   ) * ptree->n_nodes_max);
    nodes->extents     = realloc(nodes->extents,     sizeof(double) * ptree->n_nodes_max * 6);
    nodes->location_in_ancestor = realloc(nodes->location_in_ancestor, sizeof(PDM_point_tree_seq_child_t) * ptree->n_nodes_max);
  }

  int child_n_box[8] = {0, 0, 0, 0, 0, 0, 0, 0};
  int child_box_ids[8][curr_n_box];

  /* Number of points */
  int _n_points = point_range[1] - point_range[0];

  int idx[9];
  int child_id[8] = {-1, -1, -1, -1, -1, -1, -1, -1};
  int octant_mask[3] = {4, 2, 1}; /* pow(2, 2), pow(2, 1), pow(2,0) */

  int is_leaf = 1;
  if (depth < ptree->depth_max && _n_points > ptree->points_in_leaf_max && curr_n_box > ptree->points_in_leaf_max) {

    /* Choose split direction */
    double max_range = -1.;
    int split_direction = -1;

    if (ptree->tree_type == PDM_DOCTREE_LOCAL_TREE_KDTREE) {
      for (int direction = 0; direction < 3; direction++) {
        double range = extents[3+direction] - extents[direction];
        if (range > max_range) {
          max_range        = range;
          split_direction = direction;
        }
      }
      if (dbg_ptree) {
        log_trace("split_direction = %d\n", split_direction);
      }
    }

    /* Choose split point */
    double mid[3];
    if (ptree->tree_type == PDM_DOCTREE_LOCAL_TREE_KDTREE) {
      if (1) {
        mid[0] = _approx_median_point(6,
                                      split_direction,
                                      extents[split_direction],
                                      extents[3+split_direction],
                                      point_range,
                                      ptree->_pts_coord);
        // mid[0] = 0.5*(extents[split_direction] + extents[3+split_direction]);
      }
      else {
        mid[0] = _median_point(split_direction,
                               extents[split_direction],
                               extents[3+split_direction],
                               point_range,
                               ptree->_pts_coord);
      }
      if (dbg_ptree) {
        log_trace("mid = %f\n", mid[0]);
      }
    }
    else {
      for (int i = 0; i < 3; i++) {
        mid[i] = 0.5*(extents[i] + extents[3+i]);
      }
    }


    /* Count and reorder points in each child node */
    int count[8] = {0, 0, 0, 0, 0, 0, 0, 0};
    for (int i = point_range[0]; i < point_range[1]; i++) {
      int ichild = 0;
      if (ptree->tree_type == PDM_DOCTREE_LOCAL_TREE_KDTREE) {
        ichild = (ptree->_pts_coord[3*i + split_direction] > mid[0]);
      }
      else {
        for (int j = 0; j < 3; j++) {
          if (ptree->_pts_coord[3*i + j] > mid[j]) {
            ichild += octant_mask[j];
          }
        }
      }
      count[ichild]++;
    }

    if (dbg_ptree) {
      PDM_log_trace_array_int(count, n_children, "count : ");
    }


    /* Build index */
    idx[0] = 0;
    for (int ichild = 0; ichild < n_children; ichild++) {
      idx[ichild+1] = idx[ichild] + count[ichild];
      count[ichild] = 0;
    }

    for (int i = point_range[0]; i < point_range[1]; i++) {
      int ichild = 0;
      if (ptree->tree_type == PDM_DOCTREE_LOCAL_TREE_KDTREE) {
        ichild = (ptree->_pts_coord[3*i + split_direction] > mid[0]);
      }
      else {
        for (int j = 0; j < 3; j++) {
          if (ptree->_pts_coord[3*i + j] > mid[j]) {
            ichild += octant_mask[j];
          }
        }
      }

      int old = ptree->new_to_old[i];
      int pos = point_range[0] + idx[ichild] + count[ichild];

      new_to_old[pos]        = old;
      ptree->old_to_new[old] = pos;

      count[ichild]++;
    }

    for (int i = point_range[0]; i < point_range[1]; i++) {
      ptree->new_to_old[i] = new_to_old[i];
    }
    if (dbg_ptree) {
      PDM_log_trace_array_int(ptree->new_to_old + point_range[0],
                              _n_points,
                              "new_to_old: ");
    }

    /* Reorder points */
    for (int i = point_range[0]; i < point_range[1]; i++) {
      memcpy(ptree->_pts_coord + 3*i,
             ptree->pts_coord  + 3*ptree->new_to_old[i],
             sizeof(double) * 3);
    }


    for (int i = 0; i <= n_children; i++) {
      idx[i] += point_range[0];
    }

    if (dbg_ptree) {
      PDM_log_trace_array_int(idx, n_children+1, "idx : ");
    }

    /* Build leaves recursively */
    double sub_extents[6];
    for (int ichild = 0; ichild < n_children; ichild++) {

      if (idx[ichild+1] <= idx[ichild]) {
        continue;
      }

      // tmp_size++;
      // child_id[ichild] = tmp_size;

      memcpy(sub_extents, extents, sizeof(double) * 6);
      if (ichild == 0) {
        if (ptree->tree_type == PDM_DOCTREE_LOCAL_TREE_KDTREE) {
          sub_extents[3+split_direction] = mid[0];
        }
        else if (ptree->tree_type == PDM_DOCTREE_LOCAL_TREE_OCTREE) {
          sub_extents[3+0] = mid[0];
          sub_extents[3+1] = mid[1];
          sub_extents[3+2] = mid[2];
        }
      }
      else {
        if (ptree->tree_type == PDM_DOCTREE_LOCAL_TREE_KDTREE) {
          sub_extents[split_direction] = mid[0];
        }
        else if (ptree->tree_type == PDM_DOCTREE_LOCAL_TREE_OCTREE) {
          sub_extents[0] = mid[0];
          sub_extents[1] = mid[1];
          sub_extents[2] = mid[2];
        }
      }

      // int _is_leaf = (idx[ichild+1] - idx[ichild]) <=
      //                 ptree->points_in_leaf_max;
      if (dbg_ptree) {//_is_leaf) {
        log_trace("child %d, id %d, loose extents = %f %f %f  %f %f %f\n",
                  ichild, child_id[ichild],
                  sub_extents[0], sub_extents[1], sub_extents[2],
                  sub_extents[3], sub_extents[4], sub_extents[5]);
      }

      /* Tight extents (fit contained points) */
      if (1) {
        for (int j = 0; j < 3; j++) {
          sub_extents[j  ] =  HUGE_VAL;
          sub_extents[j+3] = -HUGE_VAL;
        }
        for (int ipt = idx[ichild]; ipt < idx[ichild+1]; ipt++) {
          for (int j = 0; j < 3; j++) {
            double x = ptree->_pts_coord[3*ipt+j];
            sub_extents[j  ] = PDM_MIN(sub_extents[j  ], x);
            sub_extents[j+3] = PDM_MAX(sub_extents[j+3], x);
          }
        }
      }

      /* Inflate slightly */
      for (int j = 0; j < 3; j++) {
        // if (sub_extents[j+3] < sub_extents[j] + _eps_default) {
        sub_extents[j  ] -= 0.5*_eps_default;
        sub_extents[j+3] += 0.5*_eps_default;
        // }
      }

      // if (dbg_ptree) {//_is_leaf) {
      //   log_trace("child %d, id %d, tight extents = %f %f %f  %f %f %f\n",
      //             ichild, child_id[ichild],
      //             sub_extents[0], sub_extents[1], sub_extents[2],
      //             sub_extents[3], sub_extents[4], sub_extents[5]);
      // }


      for (int idx_box = 0; idx_box < curr_n_box; idx_box++) {
        int box_id = curr_box_ids[idx_box];

        if (_intersect_box_box(3, &box_extents[6*box_id], sub_extents)) {
          child_box_ids[ichild][child_n_box[ichild]++] = box_id;
        }
      }

      if (dbg_ptree) {//_is_leaf) {
        log_trace("node_id %d, child %d, id %d, n_box = %d\n", _n_nodes, ichild, child_id[ichild], child_n_box[ichild]);
      }

      if (child_n_box[ichild] > 0) {

        tmp_size++;
        child_id[ichild] = tmp_size;

        ptree->n_nodes = tmp_size;
        is_leaf = 0;

        _build_point_tree_seq_leaves_from_boxes(_n_nodes,
                                                (PDM_point_tree_seq_child_t) ichild,
                                                depth + 1,
                                                sub_extents,
                                                box_extents,
                                                child_n_box[ichild],
                                                child_box_ids[ichild],
                                                ptree,
                                                idx + ichild,
                                                new_to_old);

        tmp_size = ptree->n_nodes;
      }
    }

  }

  /* Finalize node */
  for (int i = 0; i < 2; i++) {
    nodes->range[2*_n_nodes + i] = point_range[i];
  }

  for (int i = 0; i <= n_children; i++) {
    nodes->idx[(n_children+1)*_n_nodes + i] = idx[i];
  }

  for (int i = 0; i < 6; i++) {
    nodes->extents[6*_n_nodes + i] = extents[i];
  }

  for (int i = 0; i < n_children; i++) {
    nodes->children_id[n_children*_n_nodes + i] = child_id[i];
  }

  nodes->is_leaf[_n_nodes] = is_leaf;

  nodes->ancestor_id[_n_nodes] = ancestor_id;
  nodes->depth[_n_nodes]       = depth;

  nodes->n_points[_n_nodes] = _n_points;
  nodes->location_in_ancestor[_n_nodes] = location_in_ancestor;

  if (is_leaf) {
    if (ptree->n_leaf >= ptree->n_leaf_max) {
      if (ptree->n_leaf_max == 0) {
        ptree->n_leaf_max = 8;
      }
      else {
        ptree->n_leaf_max *= 2;
      }

      ptree->leaf_ids = realloc(ptree->leaf_ids, sizeof(int) * ptree->n_leaf_max);
      // ptree->box_ids  = realloc(ptree->leaf_ids, sizeof(int) * ptree->n_leaf_max);
      ptree->leaf_box_idx  = realloc(ptree->leaf_box_idx, sizeof(int) * (ptree->n_leaf_max+1));
    }

    if(ptree->n_leaf_box_max + curr_n_box >= ptree->n_leaf_box_max) {
      if (ptree->n_leaf_box_max == 0) {
        ptree->n_leaf_box_max = 8;
      }
      else {
        ptree->n_leaf_box_max *= 2;
      }
      ptree->n_leaf_box_max = PDM_MAX(ptree->n_leaf_box_max, ptree->leaf_box_idx[ptree->n_leaf] + curr_n_box);

      ptree->leaf_box_ids  = realloc(ptree->leaf_box_ids, sizeof(int) * ptree->n_leaf_box_max);
    }

    ptree->leaf_box_idx[ptree->n_leaf+1] = ptree->leaf_box_idx[ptree->n_leaf];
    for(int i = 0; i < curr_n_box; ++i) {
      ptree->leaf_box_ids[ptree->leaf_box_idx[ptree->n_leaf+1]++] = curr_box_ids[i];
    }

    // log_trace("leaf #%d: node id %d\n", ptree->n_leaf, _n_nodes);
    ptree->leaf_ids[ptree->n_leaf++] = _n_nodes;
  }
}








/**
 *
 * \brief Build a point_tree
 *
 * \param[in]  ptree    Current point_tree
 * .
 */

static void
_build_point_tree
(
 PDM_point_tree_seq_t *ptree
)
{
  int point_range[2];

  /* Initialization */

  ptree->n_nodes     = 0;
  ptree->n_nodes_max = 0;

  ptree->nodes = malloc(sizeof(_l_nodes_t));
  ptree->nodes->ancestor_id          = NULL;
  ptree->nodes->is_leaf              = NULL;
  ptree->nodes->location_in_ancestor = NULL;
  ptree->nodes->depth                = NULL;
  ptree->nodes->children_id          = NULL;
  ptree->nodes->range                = NULL;
  ptree->nodes->idx                  = NULL;
  ptree->nodes->n_points             = NULL;
  ptree->nodes->extents              = NULL;

  for (int i = 0; i < ptree->n_pts; i++) {
    for (int j = 0; j < 3; j++) {
      double x = ptree->pts_coord[3*i+j];
      ptree->extents[j  ] = PDM_MIN(ptree->extents[j  ], x);
      ptree->extents[j+3] = PDM_MAX(ptree->extents[j+3], x);
    }
  }

  ptree->new_to_old = malloc(sizeof(int) * ptree->n_pts);
  for (int i = 0; i < ptree->n_pts; i++) {
    ptree->new_to_old[i] = i;
  }

  ptree->old_to_new = malloc(sizeof(int) * ptree->n_pts);
  for (int i = 0; i < ptree->n_pts; i++) {
    ptree->old_to_new[i] = i;
  }

  ptree->_pts_coord = malloc(sizeof(double) * ptree->n_pts * 3);
  memcpy(ptree->_pts_coord, ptree->pts_coord, sizeof(double) * ptree->n_pts * 3);

  double delta = -1;
  for (int i = 0; i < 3; i++) {
    delta = PDM_MAX (ptree->tolerance * (ptree->extents[i + 3] - ptree->extents[i]),
                     delta);
  }
  delta = PDM_MAX (delta,_eps_default);

  for (int i = 0; i < 3; i++) {
    ptree->extents[i  ] += -1.01*delta;
    ptree->extents[i+3] += delta;;
  }

  point_range[0] = 0;
  point_range[1] = ptree->n_pts;


  /* Build point_tree recursively */

  if (dbg_ptree) {
    log_trace(">> _build_point_tree_seq_leaves\n");
  }
  int *tmp_new_to_old = malloc(sizeof(int) * ptree->n_pts);
  _build_point_tree_seq_leaves(-1,
                               (PDM_point_tree_seq_child_t) 0,
                               0,
                               ptree->extents,
                               ptree,
                               point_range,
                               tmp_new_to_old);
  free(tmp_new_to_old);


  if (ptree->n_nodes > 1) {
    ptree->n_nodes += 1;
  }


  /* Realloc */
  int n_children = PDM_point_tree_n_children_get(ptree);
  ptree->nodes->ancestor_id = realloc(ptree->nodes->ancestor_id, sizeof(int   ) * ptree->n_nodes);
  ptree->nodes->is_leaf     = realloc(ptree->nodes->is_leaf,     sizeof(int   ) * ptree->n_nodes);
  ptree->nodes->depth       = realloc(ptree->nodes->depth,       sizeof(int   ) * ptree->n_nodes);
  ptree->nodes->children_id = realloc(ptree->nodes->children_id, sizeof(int   ) * ptree->n_nodes * n_children);
  ptree->nodes->range       = realloc(ptree->nodes->range,       sizeof(int   ) * ptree->n_nodes * 2);
  ptree->nodes->idx         = realloc(ptree->nodes->idx,         sizeof(int   ) * ptree->n_nodes * (n_children+1));
  ptree->nodes->n_points    = realloc(ptree->nodes->n_points,    sizeof(int   ) * ptree->n_nodes);
  ptree->nodes->extents     = realloc(ptree->nodes->extents,     sizeof(double) * ptree->n_nodes * 6);
  ptree->nodes->location_in_ancestor = realloc(ptree->nodes->location_in_ancestor, sizeof(PDM_point_tree_seq_child_t) * ptree->n_nodes);

  _l_nodes_t *nodes = ptree->nodes;

  if(0) {
    int depth_max      = 0;
    int n_pts_leaf_max = -1;
    int n_pts_leaf_min = ptree->n_pts+1;
    int n_pts_mean     = 0;
    int n_tot_leaf     = 0;
    for (int i = 0; i < ptree->n_nodes; i++) {
      depth_max     = PDM_MAX(depth_max, nodes->depth[i]);
      if(nodes->is_leaf[i]) {
        n_pts_leaf_max = PDM_MAX(n_pts_leaf_max, nodes->range[2*i+1]-nodes->range[2*i]);
        n_pts_leaf_min = PDM_MIN(n_pts_leaf_min, nodes->range[2*i+1]-nodes->range[2*i]);
        n_pts_mean += nodes->range[2*i+1]-nodes->range[2*i];
        n_tot_leaf += 1;
      }
    }
    n_pts_mean = n_pts_mean/n_tot_leaf;
    // log_trace("point_tree stats (n_pts = %i) : depth_max = %i / n_pts_leaf_min = %i / n_pts_leaf_max = %i / n_pts_mean = %i \n",
    //           ptree->n_pts, depth_max, n_pts_leaf_min, n_pts_leaf_max, n_pts_mean );
  }

  if (dbg_ptree) {
    // PDM_log_trace_array_int(ptree->old_to_new,
    //                         ptree->n_pts,
    //                         "old_to_new : ");
    // PDM_log_trace_array_int(ptree->new_to_old,
    //                         ptree->n_pts,
    //                         "new_to_old : ");
    for (int i = 0; i < ptree->n_pts; i++) {
      if (ptree->new_to_old[ptree->old_to_new[i]] != i) {
        log_trace("!!! point %d error with old_to_new_to_old\n", i);
      }
    }

    // Dump kd-tree
    for (int i = 0; i < ptree->n_nodes; i++) {
      if (1) {//nodes->is_leaf[i]) {
        log_trace("\nNode %d :", i);
        log_trace("  depth = %d\n", nodes->depth[i]);
        log_trace("  is_leaf = %d\n", nodes->is_leaf[i]);
        // log_trace("  children_id = %d %d\n", nodes->children_id[2*i], nodes->children_id[2*i+1]);
        PDM_log_trace_array_int(nodes->children_id + n_children*i,
                                n_children,
                                "  children_id : ");
        log_trace("  extents = %f %f %f  %f %f %f\n",
                  nodes->extents[6*i+0],
                  nodes->extents[6*i+1],
                  nodes->extents[6*i+2],
                  nodes->extents[6*i+3],
                  nodes->extents[6*i+4],
                  nodes->extents[6*i+5]);
        log_trace("  point_range = %d / %d\n", nodes->range[2*i+0], nodes->range[2*i+1]);
      }
    }
  }

}





static void
_build_point_tree_from_boxes
(
       PDM_point_tree_seq_t *ptree,
 const int                   n_box,
       double               *box_extents
)
{
  int point_range[2];

  /* Initialization */

  ptree->n_nodes     = 0;
  ptree->n_nodes_max = 0;

  ptree->nodes = malloc(sizeof(_l_nodes_t));
  ptree->nodes->ancestor_id          = NULL;
  ptree->nodes->is_leaf              = NULL;
  ptree->nodes->location_in_ancestor = NULL;
  ptree->nodes->depth                = NULL;
  ptree->nodes->children_id          = NULL;
  ptree->nodes->range                = NULL;
  ptree->nodes->idx                  = NULL;
  ptree->nodes->n_points             = NULL;
  ptree->nodes->extents              = NULL;

  for (int i = 0; i < ptree->n_pts; i++) {
    for (int j = 0; j < 3; j++) {
      double x = ptree->pts_coord[3*i+j];
      ptree->extents[j  ] = PDM_MIN(ptree->extents[j  ], x);
      ptree->extents[j+3] = PDM_MAX(ptree->extents[j+3], x);
    }
  }

  ptree->new_to_old = malloc(sizeof(int) * ptree->n_pts);
  for (int i = 0; i < ptree->n_pts; i++) {
    ptree->new_to_old[i] = i;
  }

  ptree->old_to_new = malloc(sizeof(int) * ptree->n_pts);
  for (int i = 0; i < ptree->n_pts; i++) {
    ptree->old_to_new[i] = i;
  }

  ptree->_pts_coord = malloc(sizeof(double) * ptree->n_pts * 3);
  memcpy(ptree->_pts_coord, ptree->pts_coord, sizeof(double) * ptree->n_pts * 3);

  double delta = -1;
  for (int i = 0; i < 3; i++) {
    delta = PDM_MAX (ptree->tolerance * (ptree->extents[i + 3] - ptree->extents[i]),
                     delta);
  }
  delta = PDM_MAX (delta,_eps_default);

  for (int i = 0; i < 3; i++) {
    ptree->extents[i  ] += -1.01*delta;
    ptree->extents[i+3] += delta;;
  }

  point_range[0] = 0;
  point_range[1] = ptree->n_pts;


  /* Build point_tree recursively */

  int curr_n_box = 0;
  int curr_box_ids[n_box];
  for (int box_id = 0; box_id < n_box; box_id++) {
    if (_intersect_box_box(3, &box_extents[6*box_id], ptree->extents)) {
      curr_box_ids[curr_n_box++] = box_id;
    }
  }

  ptree->leaf_box_idx = malloc(2 * sizeof(int));
  ptree->leaf_box_idx[0] = 0;

  if (dbg_ptree) {
    log_trace(">> _build_point_tree_seq_leaves\n");
  }
  int *tmp_new_to_old = malloc(sizeof(int) * ptree->n_pts);
  _build_point_tree_seq_leaves_from_boxes(-1,
                                          (PDM_point_tree_seq_child_t) 0,
                                          0,
                                          ptree->extents,
                                          box_extents,
                                          curr_n_box,
                                          curr_box_ids,
                                          ptree,
                                          point_range,
                                          tmp_new_to_old);
  free(tmp_new_to_old);


  if (ptree->n_nodes > 1) {
    ptree->n_nodes += 1;
  }


  /* Realloc */
  int n_children = PDM_point_tree_n_children_get(ptree);
  ptree->nodes->ancestor_id = realloc(ptree->nodes->ancestor_id, sizeof(int   ) * ptree->n_nodes);
  ptree->nodes->is_leaf     = realloc(ptree->nodes->is_leaf,     sizeof(int   ) * ptree->n_nodes);
  ptree->nodes->depth       = realloc(ptree->nodes->depth,       sizeof(int   ) * ptree->n_nodes);
  ptree->nodes->children_id = realloc(ptree->nodes->children_id, sizeof(int   ) * ptree->n_nodes * n_children);
  ptree->nodes->range       = realloc(ptree->nodes->range,       sizeof(int   ) * ptree->n_nodes * 2);
  ptree->nodes->idx         = realloc(ptree->nodes->idx,         sizeof(int   ) * ptree->n_nodes * (n_children+1));
  ptree->nodes->n_points    = realloc(ptree->nodes->n_points,    sizeof(int   ) * ptree->n_nodes);
  ptree->nodes->extents     = realloc(ptree->nodes->extents,     sizeof(double) * ptree->n_nodes * 6);
  ptree->nodes->location_in_ancestor = realloc(ptree->nodes->location_in_ancestor, sizeof(PDM_point_tree_seq_child_t) * ptree->n_nodes);

  _l_nodes_t *nodes = ptree->nodes;

  if(0) {
    int depth_max      = 0;
    int n_pts_leaf_max = -1;
    int n_pts_leaf_min = ptree->n_pts+1;
    int n_pts_mean     = 0;
    int n_tot_leaf     = 0;
    for (int i = 0; i < ptree->n_nodes; i++) {
      depth_max     = PDM_MAX(depth_max, nodes->depth[i]);
      if(nodes->is_leaf[i]) {
        n_pts_leaf_max = PDM_MAX(n_pts_leaf_max, nodes->range[2*i+1]-nodes->range[2*i]);
        n_pts_leaf_min = PDM_MIN(n_pts_leaf_min, nodes->range[2*i+1]-nodes->range[2*i]);
        n_pts_mean += nodes->range[2*i+1]-nodes->range[2*i];
        n_tot_leaf += 1;
      }
    }
    n_pts_mean = n_pts_mean/n_tot_leaf;
    // log_trace("point_tree stats (n_pts = %i) : depth_max = %i / n_pts_leaf_min = %i / n_pts_leaf_max = %i / n_pts_mean = %i \n",
    //           ptree->n_pts, depth_max, n_pts_leaf_min, n_pts_leaf_max, n_pts_mean );
  }

  ptree->leaf_ids = realloc(ptree->leaf_ids, sizeof(int) * ptree->n_leaf);
  // PDM_log_trace_array_int(ptree->leaf_ids, ptree->n_leaf, "ptree->leaf_ids : ");

  if (dbg_ptree) {
    // PDM_log_trace_array_int(ptree->old_to_new,
    //                         ptree->n_pts,
    //                         "old_to_new : ");
    // PDM_log_trace_array_int(ptree->new_to_old,
    //                         ptree->n_pts,
    //                         "new_to_old : ");
    for (int i = 0; i < ptree->n_pts; i++) {
      if (ptree->new_to_old[ptree->old_to_new[i]] != i) {
        log_trace("!!! point %d error with old_to_new_to_old\n", i);
      }
    }

    // Dump kd-tree
    for (int i = 0; i < ptree->n_nodes; i++) {
      if (1) {//nodes->is_leaf[i]) {
        log_trace("\nNode %d :", i);
        log_trace("  depth = %d\n", nodes->depth[i]);
        log_trace("  is_leaf = %d\n", nodes->is_leaf[i]);
        // log_trace("  children_id = %d %d\n", nodes->children_id[2*i], nodes->children_id[2*i+1]);
        PDM_log_trace_array_int(nodes->children_id + n_children*i,
                                n_children,
                                "  children_id : ");
        log_trace("  extents = %f %f %f  %f %f %f\n",
                  nodes->extents[6*i+0],
                  nodes->extents[6*i+1],
                  nodes->extents[6*i+2],
                  nodes->extents[6*i+3],
                  nodes->extents[6*i+4],
                  nodes->extents[6*i+5]);
        log_trace("  point_range = %d / %d\n", nodes->range[2*i+0], nodes->range[2*i+1]);
      }
    }
  }

}



static void
_l_nodes_free
(
 PDM_point_tree_seq_t *ptree
 )
{
  if (ptree->nodes != NULL) {

    free(ptree->nodes->ancestor_id);
    free(ptree->nodes->is_leaf);
    free(ptree->nodes->depth);
    free(ptree->nodes->location_in_ancestor);
    free(ptree->nodes->children_id);
    free(ptree->nodes->range);
    free(ptree->nodes->idx);
    free(ptree->nodes->n_points);
    free(ptree->nodes->extents);

    free(ptree->nodes);

    ptree->nodes = NULL;
  }
}




inline static int
_box_dist2_min
(
 const int              dim,
 const double          *restrict extents,
 const double          *restrict coords,
 double                *restrict min_dist2
 )
{

  int inbox = 0;
  *min_dist2 = 0.;

  for (int i = 0; i < dim; i++) {
    if (coords[i] > extents[i+dim]) {
      double _min_dist2 = coords[i] - extents[dim+i];
      *min_dist2 += _min_dist2 * _min_dist2;
    }

    else if (coords[i] < extents[i]) {
      double _min_dist2 = coords[i] - extents[i];
      *min_dist2 += _min_dist2 * _min_dist2;
    }

    else {
      inbox += 1;
    }
  }

  return inbox == dim;
}


static int
_intersect_node_box_explicit
(
 int           dim,
 const double *node_extents,
 const double *box_extents,
 int          *inside
)
{
  *inside = 1;

  for (int i = 0; i < dim; i++) {
    if (node_extents[i]   > box_extents[i+3] ||
        node_extents[i+3] < box_extents[i]) {
      return 0;
    }
    else if (node_extents[i]  < box_extents[i] ||
             node_extents[i+3] > box_extents[i+3]) {
      *inside = 0;
    }
  }

  return 1;
}

/*============================================================================
 * Public function definitions
 *============================================================================*/

/**
 *
 * \brief Create a point_tree structure
 *
 * \param [in]   tree_type          Tree type
 * \param [in]   depth_max          Maximum depth
 * \param [in]   points_in_leaf_max Maximum points in a leaf
 * \param [in]   tolerance          Relative geometric tolerance
 *
 * \return     Pointer to \ref PDM_point_tree_seq object
 */

PDM_point_tree_seq_t *
PDM_point_tree_seq_create
(
 const PDM_doctree_local_tree_t tree_type,
 const int                      depth_max,
 const int                      points_in_leaf_max,
 const double                   tolerance
)
{
  if (tree_type != PDM_DOCTREE_LOCAL_TREE_OCTREE &&
      tree_type != PDM_DOCTREE_LOCAL_TREE_KDTREE) {
    PDM_error(__FILE__, __LINE__, 0,
              "Tree_type %d not implemented yet\n", (int) tree_type);
  }


  PDM_point_tree_seq_t *ptree = (PDM_point_tree_seq_t *) malloc(sizeof(PDM_point_tree_seq_t));

  ptree->tree_type = tree_type;

  ptree->depth_max          = depth_max;
  ptree->points_in_leaf_max = points_in_leaf_max;
  ptree->tolerance          = tolerance;

  ptree->nodes = NULL;

  ptree->n_nodes     = 0;
  ptree->n_nodes_max = 0;

  for (int i = 0; i < 3; i++) {
    ptree->extents[i]     =  HUGE_VAL;
    ptree->extents[i + 3] = -HUGE_VAL;
  }

  ptree->n_pts = 0;
  ptree->pts_coord  = NULL;
  ptree->_pts_coord = NULL;
  ptree->new_to_old = NULL;


  ptree->n_leaf = 0;
  ptree->n_leaf_max = 0;
  ptree->n_leaf_box_max = 0;
  ptree->leaf_ids = NULL;
  ptree->leaf_box_idx = NULL;
  ptree->leaf_box_ids = NULL;

  return ptree;
}


/**
 *
 * \brief Free a point_tree structure
 *
 * \param [in]   ptree             Pointer to \ref PDM_point_tree_seq object
 *
 */

void
PDM_point_tree_seq_free
(
 PDM_point_tree_seq_t *ptree
)
{
  if (ptree != NULL) {
    if (ptree->_pts_coord != NULL) {
      free(ptree->_pts_coord);
    }

    if (ptree->new_to_old != NULL) {
      free(ptree->new_to_old);
    }

    if (ptree->old_to_new != NULL) {
      free(ptree->old_to_new);
    }

    if(ptree->leaf_box_ids != NULL) {
      free(ptree->leaf_box_ids);
    }
    if(ptree->leaf_box_idx != NULL) {
      free(ptree->leaf_box_idx);
    }
    if(ptree->leaf_ids != NULL) {
      free(ptree->leaf_ids);
    }

    _l_nodes_free(ptree);

    free(ptree);
  }
}


/**
 *
 * \brief Set a point cloud
 *
 * \param [in]   ptree             Pointer to \ref PDM_point_tree_seq object
 * \param [in]   n_pts             Number of points in cloud
 * \param [in]   pts_coord         Point coordinates
 *
 */

void
PDM_point_tree_seq_point_cloud_set
(
       PDM_point_tree_seq_t *ptree,
 const int                   n_pts,
 const double               *pts_coord
)
{
  ptree->n_pts     = n_pts;
  ptree->pts_coord = pts_coord;
}


/**
 *
 * \brief Build point_tree
 *
 * \param [in]   ptree             Pointer to \ref PDM_point_tree_seq object
 *
 */

void
PDM_point_tree_seq_build
(
 PDM_point_tree_seq_t *ptree
)
{
  if (ptree->nodes == NULL) {
    if (dbg_ptree) {
      log_trace(">> _build_point_tree\n");
    }
    _build_point_tree(ptree);
  }

}


void
PDM_point_tree_seq_build_from_boxes
(
       PDM_point_tree_seq_t *ptree,
 const int                   n_box,
       double               *box_extents
)
{
  if (ptree->nodes == NULL) {
    if (dbg_ptree) {
      log_trace(">> _build_point_tree_from_boxes\n");
    }
    _build_point_tree_from_boxes(ptree,
                                 n_box,
                                 box_extents);
  }

}


/**
 *
 * \brief Write point_tree nodes in a VTK file
 *
 * \param [in]   ptree                 Pointer to \ref PDM_point_tree_seq object
 * \param [in]   filename              Output file name
 *
 */

void PDM_point_tree_seq_write_nodes
(
       PDM_point_tree_seq_t *ptree,
 const char                 *filename
 )
{
  _l_nodes_t *nodes = ptree->nodes;

  double tol_visu = 1e-3;
  double _ext[6];

  // write VTK
  FILE *f = fopen(filename, "w");

  fprintf(f, "# vtk DataFile Version 2.0\n");
  fprintf(f, "ptree_seq\n");
  fprintf(f, "ASCII\n");
  fprintf(f, "DATASET UNSTRUCTURED_GRID\n");

  fprintf(f, "POINTS %d double\n", 8*ptree->n_nodes);
  for (int inode = 0; inode < ptree->n_nodes; inode++) {
    double *ext = nodes->extents + 6*inode;

    if (0 && (ptree->nodes->range[2*inode+1] - ptree->nodes->range[2*inode] == 1)) {
      // Trick to visualize nodes with degenerate extents (single point)
      ext = _ext;
      for (int i = 0; i < 3; i++) {
        double x = nodes->extents[6*inode + i];
        double eps = tol_visu*(ptree->extents[i+3] - ptree->extents[i]);
        ext[i  ] = x - eps;
        ext[i+3] = x + eps;
      }
    }

    for (int k = 0; k < 2; k++) {
      for (int j = 0; j < 2; j++) {
        for (int i = 0; i < 2; i++) {
          int ii = (1-j)*i + j*(1-i);
          fprintf(f, "%f %f %f\n", ext[3*ii], ext[3*j+1], ext[3*k+2]);
        }
      }
    }
  }

  fprintf(f, "CELLS %d %d\n", ptree->n_nodes, 9*ptree->n_nodes);
  for (int inode = 0; inode < ptree->n_nodes; inode++) {
    fprintf(f, "8 ");
    for (int j = 0; j < 8; j++) {
      fprintf(f, "%d ", 8*inode+j);
    }
    fprintf(f, "\n");
  }

  fprintf(f, "CELL_TYPES %d\n", ptree->n_nodes);
  for (int i = 0; i < ptree->n_nodes; i++) {
    fprintf(f, "%d\n", 12);
  }

  fprintf(f, "CELL_DATA %d\n", ptree->n_nodes);

  fprintf(f, "FIELD node_field 3\n");
  fprintf(f, "depth 1 %d int\n", ptree->n_nodes);
  for (int i = 0; i < ptree->n_nodes; i++) {
    fprintf(f, "%d\n", ptree->nodes->depth[i]);
  }
  fprintf(f, "is_leaf 1 %d int\n", ptree->n_nodes);
  for (int i = 0; i < ptree->n_nodes; i++) {
    fprintf(f, "%d\n", ptree->nodes->is_leaf[i]);
  }
  fprintf(f, "n_pts 1 %d int\n", ptree->n_nodes);
  for (int i = 0; i < ptree->n_nodes; i++) {
    fprintf(f, "%d\n", ptree->nodes->range[2*i+1] - ptree->nodes->range[2*i]);
  }

  fclose(f);
}


/**
 *
 * \brief Get number of children per node in  point_tree
 *
 * \param [in]   ptree                 Pointer to \ref PDM_point_tree_seq object
 *
 * \return   Number of children per node in point_tree
 */

int
PDM_point_tree_n_children_get
(
 PDM_point_tree_seq_t *ptree
 )
{
  switch (ptree->tree_type) {
    case PDM_DOCTREE_LOCAL_TREE_OCTREE:
    case PDM_DOCTREE_LOCAL_TREE_LINEAR_OCTREE:
    return 8;

    case PDM_DOCTREE_LOCAL_TREE_KDTREE:
    return 2;

    default:
    PDM_error(__FILE__, __LINE__, 0,
              "Invalid tree_type %d\n", (int) ptree->tree_type);
    break;
  }

  return -1;
}


/**
 *
 * \brief Get point order in ptree
 *
 * \param [in]   ptree                 Pointer to \ref PDM_ptree_seq object
 * \param [out]  new_to_old             New to old order of points in ptree
 *
 */

void
PDM_point_tree_seq_point_new_to_old_get
(
 PDM_point_tree_seq_t  *ptree,
 int                  **new_to_old
)
{
  *new_to_old = ptree->new_to_old;
}


/**
 *
 * \brief Get point order in ptree
 *
 * \param [in]   ptree                 Pointer to \ref PDM_point_tree_seq object
 * \param [out]  old_to_new             Old to new order of points in ptree
 *
 */

void
PDM_point_tree_seq_point_old_to_new_get
(
 PDM_point_tree_seq_t  *ptree,
 int                  **old_to_new
)
{
  *old_to_new = ptree->old_to_new;
}


/**
 *
 * \brief Get point coords in point_tree's order
 *
 * \param [in]   ptree                 Pointer to \ref PDM_point_tree_seq object
 * \param [out]  pts_coord             Point coordinates
 *
 */

void
PDM_point_tree_seq_sorted_points_get
(
 PDM_point_tree_seq_t  *ptree,
 double               **pts_coord
)
{
  *pts_coord = ptree->_pts_coord;
}


/**
 *
 * \brief Get point range of a point_tree node
 *
 * \param [in]   ptree                 Pointer to \ref PDM_point_tree_seq object
 * \param [in]   node_id               Current node ID (zero-based)
 * \param [out]  point_range           Point range of current node
 *
 * \return   Number of points inside current node
 *
 */

int
PDM_point_tree_seq_point_range_get
(
       PDM_point_tree_seq_t *ptree,
 const int                   node_id,
       int                  *point_range
)
{
  assert(node_id < ptree->n_nodes);
  for (int i = 0; i < 2; i++) {
    point_range[i] = ptree->nodes->range[2*node_id + i];
  }

  return point_range[1] - point_range[0];
}



/**
 *
 * \brief Get node extents of subtree of given depth and starting from given root
 *
 * \param [in]  kdtree               Pointer to \ref PDM_kdtree_seq object
 * \param [in]  root_id              ID of subtree root
 * \param [in]  n_depth              Depth of subtree
 * \param [out] n_node               Number of subtree nodes
 * \param [out] node_ids             IDs of subtree nodes
 * \param [out] node_extents         Extents of subtree nodes
 * \param [out] node_weight          Weights of subtree nodes
 *
 */

void
PDM_point_tree_seq_extract_nodes
(
  PDM_point_tree_seq_t  *ptree,
  int                    root_id,
  int                    n_depth,
  int                   *n_node,
  int                  **node_ids,
  double               **node_extents,
  int                  **node_weight
)
{
  _l_nodes_t *nodes = ptree->nodes;

  int n_children   = PDM_point_tree_n_children_get(ptree);
  int s_pt_stack   = ((n_children - 1) * (ptree->depth_max - 1) + n_children);
  int *stack_id    = malloc (s_pt_stack * sizeof(int));
  int *stack_depth = malloc (s_pt_stack * sizeof(int));

  int *id_to_extract = malloc(ptree->n_nodes * sizeof(int));

  int n_extract = 0;
  int pos_stack = 0;
  stack_id   [pos_stack] = root_id;
  stack_depth[pos_stack] = 0;
  pos_stack++;
  while(pos_stack > 0) {

    /* Inspect node */
    --pos_stack;
    int node_id = stack_id   [pos_stack];
    int depth   = stack_depth[pos_stack];

    if(nodes->is_leaf[node_id] || depth == n_depth) {
      if(nodes->n_points[node_id] > 0) {
        id_to_extract[n_extract++] = node_id;
      }
    } else {
      for (int i = 0; i < n_children; i++) {
        int child_id = nodes->children_id[n_children*node_id+i];
        if (child_id < 0) {
          continue;
        }

        if(depth < n_depth) {
          stack_id   [pos_stack] = child_id;
          stack_depth[pos_stack] = depth + 1;
          pos_stack++;
        }
      }
    }
  }
  free(stack_id);
  free(stack_depth);

  double* _extents = malloc(n_extract * 6 * sizeof(double));
  int   * _n_pts   = malloc(n_extract *     sizeof(int   ));
  for(int i = 0; i < n_extract; ++i) {
    int node_id = id_to_extract[i];
    _n_pts[i] = nodes->n_points[node_id];
    for(int k = 0; k < 6; ++k) {
      _extents[6*i+k] = nodes->extents[6*node_id+k];
    }
  }

  *n_node       = n_extract;
  *node_ids     = id_to_extract;
  *node_extents = _extents;
  *node_weight  = _n_pts;
}


/**
 *
 * \brief Get node extents of subtree of given depth and starting from given root
 *
 * \param [in]  kdtree               Pointer to \ref PDM_kdtree_seq object
 * \param [in]  root_id              ID of subtree root
 * \param [in]  n_depth              Depth of subtree
 * \param [out] n_node               Number of subtree nodes
 * \param [out] node_ids             IDs of subtree nodes
 * \param [out] node_extents         Extents of subtree nodes
 * \param [out] node_weight          Weights of subtree nodes
 *
 */

void
PDM_point_tree_seq_extract_extents_by_child_ids
(
  PDM_point_tree_seq_t  *ptree,
  const int              n_node_to_extract,
  const int             *node_ids_to_extract,
        int             *n_extract_child,
        int            **node_to_child_idx,
        int            **extract_child_id,
        int            **extract_is_leaf,
        double         **extract_extents
)
{
  _l_nodes_t *nodes = ptree->nodes;

  int    *_node_to_child_idx = malloc((n_node_to_extract + 1 ) * sizeof(int   ));
  _node_to_child_idx[0] = 0;

  int n_children   = PDM_point_tree_n_children_get(ptree);
  double *_extract_extents   = malloc(n_node_to_extract * n_children * 6 * sizeof(double));
  int    *_extract_child_id  = malloc(n_node_to_extract * n_children     * sizeof(int   ));
  int    *_extract_is_leaf   = malloc(n_node_to_extract * n_children     * sizeof(int   ));
  int     _n_extract_child   = 0;

  for(int i_node_to_extract = 0; i_node_to_extract < n_node_to_extract; ++i_node_to_extract) {
    int node_id = node_ids_to_extract[i_node_to_extract];
    _node_to_child_idx[i_node_to_extract+1] = _node_to_child_idx[i_node_to_extract];

    if(nodes->is_leaf[node_id]) {
      _extract_child_id[_n_extract_child] = 1;
      for(int k = 0; k < 6; ++k) {
        _extract_extents[6*_n_extract_child+k] = nodes->extents[6*node_id+k];
      }
      _n_extract_child++;

      _node_to_child_idx[i_node_to_extract+1]++;

    } else {
      const int *_child_ids = nodes->children_id + n_children*node_id;
      for (int i = 0; i < n_children; i++) {
        int child_id = _child_ids[i];
        if (child_id < 0) {
          continue;
        }

        _extract_child_id[_n_extract_child] = child_id;
        _extract_is_leaf [_n_extract_child] = 1;

        for(int k = 0; k < 6; ++k) {
          _extract_extents[6*_n_extract_child+k] = nodes->extents[6*child_id+k];
        }

        _node_to_child_idx[i_node_to_extract+1]++;
        _n_extract_child++;
      }
    }
  }

  _extract_extents  = realloc(_extract_extents , _n_extract_child * 6 * sizeof(double));
  _extract_child_id = realloc(_extract_child_id, _n_extract_child     * sizeof(int   ));
  _extract_is_leaf  = realloc(_extract_is_leaf , _n_extract_child     * sizeof(int   ));

  *n_extract_child   = _n_extract_child;
  *node_to_child_idx = _node_to_child_idx;
  *extract_child_id  = _extract_child_id;
  *extract_is_leaf   = _extract_is_leaf;
  *extract_extents   = _extract_extents;
}


/**
 *
 * \brief Look for closest points stored inside a point_tree
 *
 * \param [in]   ptree                  Pointer to \ref PDM_point_tree_seq object
 * \param [in]   n_pts                  Number of points
 * \param [in]   pts                    Point Coordinates
 * \param [out]  closest_kdtree_pt_id   Closest point in kdtree ID (zero-based)
 * \param [out]  closest_kdtree_pt_dist Closest point in kdtree distance
 *
 */

void
PDM_point_tree_seq_closest_point
(
PDM_point_tree_seq_t *ptree,
const int             n_pts,
double               *pts,
int                  *closest_ptree_pt_id,
double               *closest_ptree_pt_dist2
)
{

  const int n_children = PDM_point_tree_n_children_get(ptree);


  int s_pt_stack = ((n_children - 1) * (ptree->depth_max - 1) + n_children);
  int sort_child[n_children];
  double dist_child[n_children];
  int inbox_child[n_children];

  int    *stack           = malloc (sizeof(int   ) * s_pt_stack);
  int    *inbox_stack     = malloc (sizeof(int   ) * s_pt_stack);
  double *min_dist2_stack = malloc (sizeof(double) * s_pt_stack);

  _l_nodes_t *nodes = ptree->nodes;

  int dim = 3;

  for (int i = 0; i < n_pts; i++) {

    int pos_stack = 0;
    const double *_pt = pts + dim * i;

    /* Init stack */

    closest_ptree_pt_id   [i] = -1;
    closest_ptree_pt_dist2[i] = HUGE_VAL;

    stack[pos_stack] = 0; /* push root in the stack */

    double _min_dist2;
    int inbox1 = _box_dist2_min(dim,
                                &nodes->extents[0],
                                _pt,
                                &_min_dist2);

    inbox_stack[pos_stack]     = inbox1;
    min_dist2_stack[pos_stack] = _min_dist2;
    pos_stack++;

    while (pos_stack > 0) {
      int node_id = stack[--pos_stack];

      double min_dist2 = min_dist2_stack[pos_stack];
      int    inbox     = inbox_stack[pos_stack];

      if ((min_dist2 <= closest_ptree_pt_dist2[i]) || (inbox == 1)) {

        if (!nodes->is_leaf[node_id]) {

          /* Sort children and store them into the stack */

          const int *_child_ids = nodes->children_id + n_children*node_id;

          for (int j = 0; j < n_children; j++) {
            dist_child[j] = HUGE_VAL;
          }

          int n_selec = 0;
          for (int j = 0; j < n_children; j++) {


            int child_id = _child_ids[j];


            int child_inbox = 0;

            if (child_id != -1) {

              double child_min_dist2;

              child_inbox = _box_dist2_min(dim,
                                           &nodes->extents[6*child_id],
                                           _pt,
                                           &child_min_dist2);

              int i1 = 0;
              for (i1 = n_selec;
                   (i1 > 0) && (dist_child[i1-1] > child_min_dist2) ; i1--) {
                dist_child[i1]  = dist_child[i1-1];
                sort_child[i1]  = sort_child[i1-1];
                inbox_child[i1] = inbox_child[i1-1];
              }

              sort_child[i1]  = child_id;
              dist_child[i1]  = child_min_dist2;
              inbox_child[i1] = child_inbox;

              n_selec += 1;

            }
          }

          for (int j = 0; j < n_selec; j++) {
            int j1 = n_selec- 1 - j;
            int child_id = sort_child[j1];
            if (child_id != -1) {
              if ((dist_child[j1] < closest_ptree_pt_dist2[i]) &&
                  (nodes->n_points[child_id] > 0)) {

                min_dist2_stack[pos_stack] = dist_child[j1];
                inbox_stack[pos_stack]     = inbox_child[j1];

                stack[pos_stack++] = child_id; /* push root in th stack */
              }
            }
          }
        }

        else {


          for (int j = 0; j < nodes->n_points[node_id]; j++) {

            double point_dist2 = 0;
            int ipt = nodes->range[2*node_id] + j;
            double *_coords = ptree->_pts_coord + 3*ipt;

            for (int k = 0; k < dim; k++) {
              point_dist2 += (_coords[k] - _pt[k]) *
                             (_coords[k] - _pt[k]);
            }

            if (point_dist2 < closest_ptree_pt_dist2[i]) {
              closest_ptree_pt_id[i]    = ptree->new_to_old[ipt];
              closest_ptree_pt_dist2[i] = point_dist2;
            }
          }
        }
      }
    }

  }

  free (inbox_stack);
  free (min_dist2_stack);
  free (stack);

}


/**
 *
 * \brief Get points located inside a set of boxes
 *
 * \param [in]   ptree                  Pointer to \ref PDM_point_tree_seq object
 * \param [in]   n_box                  Number of boxes
 * \param [in]   box_extents            Extents of boxes
 * \param [out]  box_pts_idx            Index of points located in boxes
 * \param [out]  box_pts                Local ids of points located in boxes (zero-based)
 *
 */

void
PDM_point_tree_seq_points_inside_boxes
(
       PDM_point_tree_seq_t  *ptree,
 const int                    n_box,
 const double                 box_extents[],
       int                  **box_pts_idx,
       int                  **box_pts
)
{
  *box_pts_idx = malloc (sizeof(int) * (n_box + 1));
  int *_box_pts_idx = *box_pts_idx;
  _box_pts_idx[0] = 0;

  if (n_box < 1) {
    *box_pts = malloc (sizeof(int) * _box_pts_idx[n_box]);
    return;
  }

  const int n_children = PDM_point_tree_n_children_get(ptree);


  _l_nodes_t *nodes = ptree->nodes;

  int s_pt_stack = ((n_children - 1) * (ptree->depth_max - 1) + n_children);
  int *stack_id  = malloc (s_pt_stack * sizeof(int));

  int node_inside_box;
  int intersect;

  int tmp_size = 4 * n_box;
  *box_pts = malloc (sizeof(int) * tmp_size);
  int *_box_pts = *box_pts;

  for (int ibox = 0; ibox < n_box; ibox++) {
    int dbg_enabled = 0;
    if (dbg_enabled) {
      log_trace("box %d\n", ibox);
    }

    _box_pts_idx[ibox+1] = _box_pts_idx[ibox];

    const double *_box_extents = box_extents + 6*ibox;
    const double *box_min      = box_extents + 6*ibox;
    const double *box_max      = box_min + 3;

    intersect = _intersect_node_box_explicit (3,
                                              &nodes->extents[0],
                                              _box_extents,
                                              &node_inside_box);

    if (!intersect) {
      continue;
    }


    if (node_inside_box) {
      /* The box must contain all points */
      if (dbg_enabled) {
        log_trace("    add pts with lnum %d through %d\n", nodes->range[0], nodes->range[1] + nodes->n_points[0]);
      }
      int new_size = _box_pts_idx[ibox+1] + nodes->n_points[0];

      if (tmp_size <= new_size) {
        tmp_size = PDM_MAX (2*tmp_size, new_size);
        *box_pts = realloc (*box_pts, sizeof(int) * tmp_size);
        _box_pts = *box_pts;

      }

      for (int j = 0; j < nodes->n_points[0]; j++) {
        _box_pts[_box_pts_idx[ibox+1]++] = ptree->new_to_old[nodes->range[0] + j];
      }
      continue;
    } /* End node_inside_box */


    /* Push root in stack */
    int pos_stack = 0;
    stack_id[pos_stack++] = 0;

    while (pos_stack > 0) {
      int node_id = stack_id[--pos_stack];

      if (dbg_enabled) {
        log_trace("  node %d, range=%d/%d, n_points=%d, leaf_id=%d\n",
                  node_id,
                  nodes->range[2*node_id], nodes->range[2*node_id+1],
                  nodes->n_points[node_id],
                  nodes->is_leaf[node_id]);
      }

      /* is leaf */
      if(nodes->is_leaf[node_id]) {


        for (int i = 0; i < nodes->n_points[node_id]; i++) {
          int ipt = nodes->range[2*node_id] + i;
          const double *_pt = ptree->_pts_coord + 3*ipt;

          int pt_inside_box = 1;
          for (int idim = 0; idim < 3; idim++) {
            if (_pt[idim] < box_min[idim] || _pt[idim] > box_max[idim]) {
              pt_inside_box = 0;
              break;
            }
          }

          if (pt_inside_box) {
            if (_box_pts_idx[ibox+1] >= tmp_size) {
              tmp_size = PDM_MAX (2*tmp_size, _box_pts_idx[ibox+1] + 1);
              *box_pts = realloc (*box_pts, sizeof(int) * tmp_size);
              _box_pts = *box_pts;
            }

            _box_pts[_box_pts_idx[ibox+1]++] = ptree->new_to_old[ipt];
          }
        }
      } else { /* Internal nodes */

        const int *_child_ids = nodes->children_id + n_children*node_id;
        for (int i = 0; i < n_children; i++) {
          int child_id = _child_ids[i];
          if (child_id < 0) {
            continue;
          }

          if (dbg_enabled) {
            log_trace("    child %d: id=%d, range=%d/%d, n_points=%d, leaf_id=%d\n",
                   i,
                   child_id,
                   nodes->range[2*child_id+0],
                   nodes->range[2*child_id+1],
                   nodes->n_points[child_id],
                   nodes->is_leaf[child_id]);
            log_trace("    pts_extents = %f %f %f %f %f %f\n",
                   nodes->extents[6*child_id+0],
                   nodes->extents[6*child_id+1],
                   nodes->extents[6*child_id+2],
                   nodes->extents[6*child_id+3],
                   nodes->extents[6*child_id+4],
                   nodes->extents[6*child_id+5]);
          }

          intersect = _intersect_node_box_explicit (3,
                                                    // child_node->extents,
                                                    nodes->extents + 6*child_id,
                                                    _box_extents,
                                                    &node_inside_box);

          if (dbg_enabled) {
            log_trace("    intersect = %d\n", intersect);
          }

          if (intersect) {
            if (node_inside_box) {
              /* The box must contain all points */
              if (dbg_enabled) {
                log_trace("    add pts with lnum %d through %d\n", nodes->range[2*child_id+0], nodes->range[2*child_id+1]);
              }

              int new_size = _box_pts_idx[ibox+1] + nodes->n_points[child_id];

              if (tmp_size <= new_size) {
                tmp_size = PDM_MAX (2*tmp_size, new_size);
                *box_pts = realloc (*box_pts, sizeof(int) * tmp_size);
                _box_pts = *box_pts;
              }

              for (int j = 0; j < nodes->n_points[child_id]; j++) {
                _box_pts[_box_pts_idx[ibox+1]++] = ptree->new_to_old[nodes->range[2*child_id] + j];
              }
            }

            else {
              /* Push child in stack */
              stack_id[pos_stack++] = child_id;
            }
          }
        } // End of loop on children
      }

    } /* End While */
  } /* End boxe loop */

  free (stack_id);
  *box_pts = realloc (*box_pts, sizeof(int) * _box_pts_idx[n_box]);
}



/**
 *
 * \brief Look for points inside at set of balls
 *
 * \param [in]  ptree                Pointer to \ref PDM_point_tree_seq object
 * \param [in]  n_ball               Number of balls
 * \param [in]  ball_center          Center of balls (size = \ref n_ball * 3)
 * \param [in]  ball_radius2         Squared radius of balls (size = \ref n_ball)
 * \param [out] ball_pts_idx         Index for ball->points graph (size \ref n_ball + 1)
 * \param [out] ball_pts             Ball->points graph (zero-based IDs)
 * \param [out] ball_pts_dist2       Distance from points to ball centers
 *
 */

void
PDM_point_tree_seq_points_inside_balls
(
       PDM_point_tree_seq_t  *ptree,
 const int                    n_ball,
       double                *ball_center,
       double                *ball_radius2,
       int                  **ball_pts_idx,
       int                  **ball_pts,
       double               **ball_pts_dist2
)
{
  int dbg = 0;
  const int n_children = PDM_point_tree_n_children_get(ptree);

  int s_pt_stack = ((n_children - 1) * (ptree->depth_max - 1) + n_children);


  *ball_pts_idx = malloc(sizeof(int) * (n_ball + 1));
  int *pib_idx = *ball_pts_idx;
  pib_idx[0] = 0;

  int s_pib = 4*n_ball;
  *ball_pts       = malloc(sizeof(int   ) * s_pib);
  *ball_pts_dist2 = malloc(sizeof(double) * s_pib);

  int    *pib_l_num = *ball_pts;
  double *pib_dist2 = *ball_pts_dist2;


  _l_nodes_t *nodes = ptree->nodes;


  int *stack = malloc(sizeof(int) * s_pt_stack);


  for (int iball = 0; iball < n_ball; iball++) {

    pib_idx[iball+1] = pib_idx[iball];

    double *_center  = ball_center + 3*iball;
    double  _radius2 = ball_radius2[iball];

    if (dbg) {
      log_trace("ball %d, center = %f %f %f, radius2 = %f\n",
                iball,
                _center[0], _center[1], _center[2], _radius2);
    }


    /* Start by root */
    int pos_stack = 0;
    double min_dist2;
    int inside_box = _box_dist2_min(3,
                                    &nodes->extents[0],
                                    _center,
                                    &min_dist2);

    if (inside_box || min_dist2 <= _radius2) {
      if (dbg) {
        log_trace("  push root\n");
      }
      stack[pos_stack++] = 0;
    }


    while (pos_stack > 0) {

      int node_id = stack[--pos_stack];
      if (dbg) {
        log_trace("  node_id = %d (is_leaf? %d)\n", node_id, nodes->is_leaf[node_id]);
      }

      if (nodes->is_leaf[node_id]) {
        /* Leaf node */
        for (int i = nodes->range[2*node_id]; i < nodes->range[2*node_id+1]; i++) {
          double *_pt = ptree->_pts_coord + 3*i;

          double dist2 = 0.;
          for (int j = 0; j < 3; j++) {
            double delta = _pt[j] - _center[j];
            dist2 += delta*delta;
          }

          if (dbg) {
            log_trace("    pt %d: dist2 = %f\n", i, dist2);
          }

          if (dist2 <= _radius2) {
            /* Check size and realloc if necessary */
            if (pib_idx[iball+1] >= s_pib) {
              s_pib *= 2;

              *ball_pts       = realloc(*ball_pts,       sizeof(int   ) * s_pib);
              *ball_pts_dist2 = realloc(*ball_pts_dist2, sizeof(double) * s_pib);

              pib_l_num = *ball_pts;
              pib_dist2 = *ball_pts_dist2;
            }

            /* Add point */
            pib_l_num[pib_idx[iball+1]] = ptree->new_to_old[i];
            pib_dist2[pib_idx[iball+1]] = dist2;

            pib_idx[iball+1]++;
          }
        } // End of loop on current leaf's points
      }

      else {
        /* Internal node */
        for (int ichild = 0; ichild < n_children; ichild++) {

          int child_id = nodes->children_id[n_children*node_id + ichild];

          if (child_id < 0) {
            continue;
          }

          if (nodes->n_points[child_id] == 0) {
            continue;
          }

          inside_box = _box_dist2_min(3,
                                      &nodes->extents[6*child_id],
                                      _center,
                                      &min_dist2);

          if (inside_box || min_dist2 <= _radius2) {
            stack[pos_stack++] = child_id;
            if (dbg) {
              log_trace("  push child %d (%d)", ichild, child_id);
            }
          }

        }

      }


    } // End of while loop


  } // End of loop on points
  free(stack);

  s_pib = pib_idx[n_ball];
  *ball_pts       = realloc(*ball_pts,       sizeof(int   ) * s_pib);
  *ball_pts_dist2 = realloc(*ball_pts_dist2, sizeof(double) * s_pib);

}



PDM_point_tree_seq_shm_t *
PDM_point_tree_make_shared
(
  PDM_point_tree_seq_t *local_ptree,
  PDM_MPI_Comm          comm_shared
)
{
  int n_rank_in_shm, i_rank_in_shm;
  PDM_MPI_Comm_rank (comm_shared, &i_rank_in_shm);
  PDM_MPI_Comm_size (comm_shared, &n_rank_in_shm);

  PDM_point_tree_seq_shm_t* shm_ptree = malloc(sizeof(PDM_point_tree_seq_shm_t));

  shm_ptree->tree_type = local_ptree->tree_type;

  shm_ptree->comm_shared = comm_shared;
  shm_ptree->ptrees      = malloc(n_rank_in_shm * sizeof(PDM_point_tree_seq_t));

  /*
   * Exchange size
   */
  int s_shm_data_in_rank[2] = {0};
  s_shm_data_in_rank[0] = local_ptree->n_nodes;
  s_shm_data_in_rank[1] = local_ptree->n_pts;
  int *s_shm_data_in_all_nodes = malloc(2 * n_rank_in_shm * sizeof(int));

  PDM_MPI_Allgather(s_shm_data_in_rank     , 2, PDM_MPI_INT,
                    s_shm_data_in_all_nodes, 2, PDM_MPI_INT, comm_shared);

  int *shared_nodes_idx = malloc((n_rank_in_shm+1) * sizeof(int));
  int *shared_pts_idx   = malloc((n_rank_in_shm+1) * sizeof(int));
  shared_nodes_idx[0] = 0;
  shared_pts_idx  [0] = 0;
  for(int i = 0; i < n_rank_in_shm; ++i) {
    shared_nodes_idx[i+1] = shared_nodes_idx[i] + s_shm_data_in_all_nodes[2*i  ];
    shared_pts_idx  [i+1] = shared_pts_idx  [i] + s_shm_data_in_all_nodes[2*i+1];
  }

  int n_nodes_shared_tot = shared_nodes_idx[n_rank_in_shm];
  int n_pts_shared_tot   = shared_pts_idx  [n_rank_in_shm];

  /* Nodes */
  int n_children = PDM_point_tree_n_children_get(local_ptree);
  shm_ptree->w_is_leaf     = PDM_mpi_win_shared_create(n_nodes_shared_tot,              sizeof(int   ), comm_shared);
  shm_ptree->w_children_id = PDM_mpi_win_shared_create(n_nodes_shared_tot * n_children, sizeof(int   ), comm_shared);
  shm_ptree->w_range       = PDM_mpi_win_shared_create(n_nodes_shared_tot * 2,          sizeof(int   ), comm_shared);
  shm_ptree->w_n_points    = PDM_mpi_win_shared_create(n_nodes_shared_tot,              sizeof(int   ), comm_shared);
  shm_ptree->w_extents     = PDM_mpi_win_shared_create(n_nodes_shared_tot * 6,          sizeof(double), comm_shared);
  int    *ptr_is_leaf      = PDM_mpi_win_shared_get(shm_ptree->w_is_leaf);
  int    *ptr_children_id  = PDM_mpi_win_shared_get(shm_ptree->w_children_id);
  int    *ptr_range        = PDM_mpi_win_shared_get(shm_ptree->w_range);
  int    *ptr_n_points     = PDM_mpi_win_shared_get(shm_ptree->w_n_points);
  double *ptr_extents      = PDM_mpi_win_shared_get(shm_ptree->w_extents);


  /* Points */
  shm_ptree->w_pts_coord  = PDM_mpi_win_shared_create(n_pts_shared_tot * 3, sizeof(double), comm_shared);
  shm_ptree->w_new_to_old = PDM_mpi_win_shared_create(n_pts_shared_tot,     sizeof(int   ), comm_shared);
  shm_ptree->w_old_to_new = PDM_mpi_win_shared_create(n_pts_shared_tot,     sizeof(int   ), comm_shared);

  double *ptr_pts_coord   = PDM_mpi_win_shared_get(shm_ptree->w_pts_coord);
  int    *ptr_new_to_old  = PDM_mpi_win_shared_get(shm_ptree->w_new_to_old);
  int    *ptr_old_to_new  = PDM_mpi_win_shared_get(shm_ptree->w_old_to_new);

  /*
   *  Set window pointers
   */
  shm_ptree->shm_n_nodes     = malloc(sizeof(int     ) * n_rank_in_shm);
  shm_ptree->shm_is_leaf     = malloc(sizeof(int    *) * n_rank_in_shm);
  shm_ptree->shm_children_id = malloc(sizeof(int    *) * n_rank_in_shm);
  shm_ptree->shm_range       = malloc(sizeof(int    *) * n_rank_in_shm);
  shm_ptree->shm_n_points    = malloc(sizeof(int    *) * n_rank_in_shm);
  shm_ptree->shm_extents     = malloc(sizeof(double *) * n_rank_in_shm);
  shm_ptree->shm_n_pts       = malloc(sizeof(int     ) * n_rank_in_shm);
  shm_ptree->shm_pts_coord   = malloc(sizeof(double *) * n_rank_in_shm);
  shm_ptree->shm_new_to_old  = malloc(sizeof(int    *) * n_rank_in_shm);
  shm_ptree->shm_old_to_new  = malloc(sizeof(int    *) * n_rank_in_shm);

  for(int i = 0; i < n_rank_in_shm; ++i) {

    /* Nodes */
    shm_ptree->shm_n_nodes    [i] = s_shm_data_in_all_nodes[2*i];
    shm_ptree->shm_is_leaf    [i] = &ptr_is_leaf    [shared_nodes_idx[i]];
    shm_ptree->shm_children_id[i] = &ptr_children_id[shared_nodes_idx[i] * n_children];
    shm_ptree->shm_range      [i] = &ptr_range      [shared_nodes_idx[i] * 2];
    shm_ptree->shm_n_points   [i] = &ptr_n_points   [shared_nodes_idx[i]];
    shm_ptree->shm_extents    [i] = &ptr_extents    [shared_nodes_idx[i] * 6];

    /* Points */
    shm_ptree->shm_n_pts      [i] = s_shm_data_in_all_nodes[2*i+1];
    shm_ptree->shm_pts_coord  [i] = &ptr_pts_coord         [shared_pts_idx[i] * 3];
    shm_ptree->shm_new_to_old [i] = &ptr_new_to_old        [shared_pts_idx[i]];
    shm_ptree->shm_old_to_new [i] = &ptr_old_to_new        [shared_pts_idx[i]];

  }


  free(s_shm_data_in_all_nodes);
  shm_ptree->shared_nodes_idx = shared_nodes_idx;
  shm_ptree->shared_pts_idx   = shared_pts_idx;


  /*
   * Copy in window
   */

  /* Nodes */
  _l_nodes_t *local_nodes = local_ptree->nodes;
  memcpy(shm_ptree->shm_is_leaf[i_rank_in_shm],
         local_nodes->is_leaf,
         sizeof(int) * local_ptree->n_nodes);
  // PDM_log_trace_array_int(local_nodes->is_leaf, local_ptree->n_nodes, "local_is_leaf : ");
  // PDM_log_trace_array_int(shm_ptree->shm_is_leaf[i_rank_in_shm], local_ptree->n_nodes, "  shm_is_leaf : ");
  memcpy(shm_ptree->shm_children_id[i_rank_in_shm],
         local_nodes->children_id,
         sizeof(int) * local_ptree->n_nodes * n_children);
  memcpy(shm_ptree->shm_range[i_rank_in_shm],
         local_nodes->range,
         sizeof(int) * local_ptree->n_nodes * 2);
  memcpy(shm_ptree->shm_n_points[i_rank_in_shm],
         local_nodes->n_points,
         sizeof(int) * local_ptree->n_nodes);
  memcpy(shm_ptree->shm_extents[i_rank_in_shm],
         local_nodes->extents,
         sizeof(double) * local_ptree->n_nodes * 6);

  /* Points */
  memcpy(shm_ptree->shm_pts_coord[i_rank_in_shm],
         local_ptree->_pts_coord,
         sizeof(double) * local_ptree->n_pts * 3);
  memcpy(shm_ptree->shm_new_to_old[i_rank_in_shm],
         local_ptree->new_to_old,
         sizeof(int) * local_ptree->n_pts);
  memcpy(shm_ptree->shm_old_to_new[i_rank_in_shm],
         local_ptree->old_to_new,
         sizeof(int) * local_ptree->n_pts);

  PDM_MPI_Barrier(comm_shared);

  return shm_ptree;
}

void
PDM_point_tree_seq_shm_free
(
 PDM_point_tree_seq_shm_t* shm_ptree
)
{

  free(shm_ptree->shm_n_nodes    );
  free(shm_ptree->shm_is_leaf    );
  free(shm_ptree->shm_children_id);
  free(shm_ptree->shm_range      );
  free(shm_ptree->shm_n_points   );
  free(shm_ptree->shm_extents    );
  free(shm_ptree->shm_n_pts      );
  free(shm_ptree->shm_pts_coord  );
  free(shm_ptree->shm_new_to_old );
  free(shm_ptree->shm_old_to_new );

  PDM_mpi_win_shared_free(shm_ptree->w_is_leaf    );
  PDM_mpi_win_shared_free(shm_ptree->w_children_id);
  PDM_mpi_win_shared_free(shm_ptree->w_range      );
  PDM_mpi_win_shared_free(shm_ptree->w_n_points   );
  PDM_mpi_win_shared_free(shm_ptree->w_extents    );

  PDM_mpi_win_shared_free(shm_ptree->w_pts_coord );
  PDM_mpi_win_shared_free(shm_ptree->w_new_to_old);
  PDM_mpi_win_shared_free(shm_ptree->w_old_to_new);

  free(shm_ptree->shared_nodes_idx);
  free(shm_ptree->shared_pts_idx  );

  free(shm_ptree->ptrees);
  free(shm_ptree);
}



/**
 *
 * \brief Get point coords in point_tree's order
 *
 * \param [in]   ptree                 Pointer to \ref PDM_point_tree_seq object
 * \param [out]  pts_coord             Point coordinates
 *
 */

void
PDM_point_tree_seq_shm_sorted_points_get
(
 PDM_point_tree_seq_shm_t  *shm_tree,
 int                        i_shm,
 double                   **pts_coord
)
{
  *pts_coord = shm_tree->shm_pts_coord[i_shm];
}


/**
 *
 * \brief Get point order in ptree
 *
 * \param [in]   ptree                 Pointer to \ref PDM_ptree_seq object
 * \param [out]  new_to_old             New to old order of points in ptree
 *
 */

void
PDM_point_tree_seq_shm_point_new_to_old_get
(
 PDM_point_tree_seq_shm_t  *shm_tree,
 int                        i_shm,
 int                      **new_to_old
)
{
  *new_to_old = shm_tree->shm_new_to_old[i_shm];
}


/**
 *
 * \brief Get point order in ptree
 *
 * \param [in]   ptree                 Pointer to \ref PDM_ptree_seq object
 * \param [out]  new_to_old             New to old order of points in ptree
 *
 */

void
PDM_point_tree_seq_shm_point_old_to_new_get
(
 PDM_point_tree_seq_shm_t  *shm_tree,
 int                        i_shm,
 int                      **old_to_new
)
{
  *old_to_new = shm_tree->shm_old_to_new[i_shm];
}




void
PDM_point_tree_seq_points_inside_boxes_shared
(
       PDM_point_tree_seq_shm_t  *shm_ptree,
 const int                        i_shm_rank,
 const int                        n_box,
 const double                     box_extents[],
 // const PDM_g_num_t                box_g_num[],
       int                      **box_pts_idx,
       int                      **box_pts
 )
{

  *box_pts_idx = malloc (sizeof(int) * (n_box + 1));
  int *_box_pts_idx = *box_pts_idx;
  _box_pts_idx[0] = 0;

  if (n_box < 1) {
    *box_pts = malloc (sizeof(int) * _box_pts_idx[n_box]);
    return;
  }

   // const int n_children = PDM_point_tree_n_children_get(ptree);
  int n_children = 2;
  if (shm_ptree->tree_type == PDM_DOCTREE_LOCAL_TREE_KDTREE) {
    n_children = 2;
  }
  else if (shm_ptree->tree_type == PDM_DOCTREE_LOCAL_TREE_OCTREE) {
    n_children = 8;
  }
  else {
    PDM_error(__FILE__, __LINE__, 0,
              "Invalid tree_type %d\n", (int) shm_ptree->tree_type);
  }

  const int depth_max = 31;
  // int     n_nodes     = shm_ptree->shm_n_nodes    [i_shm_rank];
  int    *is_leaf     = shm_ptree->shm_is_leaf    [i_shm_rank];
  int    *children_id = shm_ptree->shm_children_id[i_shm_rank];
  int    *range       = shm_ptree->shm_range      [i_shm_rank];
  int    *n_points    = shm_ptree->shm_n_points   [i_shm_rank];
  double *extents     = shm_ptree->shm_extents    [i_shm_rank];
  double *pts_coord   = shm_ptree->shm_pts_coord  [i_shm_rank];
  int    *new_to_old  = shm_ptree->shm_new_to_old [i_shm_rank];

  int s_pt_stack = ((n_children - 1) * (depth_max - 1) + n_children);
  int *stack_id  = malloc (s_pt_stack * sizeof(int));


  int node_inside_box;
  int intersect;

  int tmp_size = 4 * n_box;
  *box_pts = malloc (sizeof(int) * tmp_size);
  int *_box_pts = *box_pts;


  for (int ibox = 0; ibox < n_box; ibox++) {

    int dbg_enabled = 0;
    if (dbg_enabled) {
      log_trace("box %d\n", ibox);
    }

    _box_pts_idx[ibox+1] = _box_pts_idx[ibox];

    const double *_box_extents = box_extents + 6*ibox;
    const double *box_min      = box_extents + 6*ibox;
    const double *box_max      = box_min + 3;

    intersect = _intersect_node_box_explicit (3,
                                              &extents[0],
                                              _box_extents,
                                              &node_inside_box);

    if (!intersect) {
      continue;
    }

    if (node_inside_box) {
      /* The box must contain all points */
      if (dbg_enabled) {
        log_trace("    add pts with lnum %d through %d\n", range[0], range[1]);
      }
      int new_size = _box_pts_idx[ibox+1] + n_points[0];

      if (tmp_size <= new_size) {
        tmp_size = PDM_MAX (2*tmp_size, new_size);
        *box_pts = realloc (*box_pts, sizeof(int) * tmp_size);
        _box_pts = *box_pts;

      }

      for (int j = 0; j < n_points[0]; j++) {
        _box_pts[_box_pts_idx[ibox+1]++] = range[0] + j;
      }
      continue;
    } /* End node_inside_box */


    /* Push root in stack */
    int pos_stack = 0;
    stack_id[pos_stack++] = 0;

    while (pos_stack > 0) {

      int node_id = stack_id[--pos_stack];

      if (dbg_enabled) {
        log_trace("  node %d, range=%d/%d, n_points=%d, is_leaf=%d\n",
                  node_id,
                  range[2*node_id], range[2*node_id+1],
                  n_points[node_id],
                  is_leaf[node_id]);
      }

      if(is_leaf[node_id]) {
        /* Leaf node */

        for (int i = 0; i < n_points[node_id]; i++) {
          int ipt = range[2*node_id] + i;
          const double *_pt = pts_coord + 3*ipt;

          int pt_inside_box = 1;
          for (int idim = 0; idim < 3; idim++) {
            if (_pt[idim] < box_min[idim] || _pt[idim] > box_max[idim]) {
              pt_inside_box = 0;
              break;
            }
          }

          if (pt_inside_box) {
            if (_box_pts_idx[ibox+1] >= tmp_size) {
              tmp_size = PDM_MAX (2*tmp_size, _box_pts_idx[ibox+1] + 1);
              *box_pts = realloc (*box_pts, sizeof(int) * tmp_size);
              _box_pts = *box_pts;
            }

            _box_pts[_box_pts_idx[ibox+1]++] = new_to_old[ipt];
          }
        }
      }
      else {
        /* Internal nodes */

        const int *_child_ids = children_id + n_children*node_id;
        for (int i = 0; i < n_children; i++) {
          int child_id = _child_ids[i];
          if (child_id < 0) {
            continue;
          }

          if (dbg_enabled) {
            // log_trace("    child %d / %d\n", child_id, n_nodes);
            log_trace("    child %d: id=%d, range=%d/%d, n_points=%d, leaf_id=%d\n",
                   i,
                   child_id,
                   range[2*child_id+0],
                   range[2*child_id+1],
                   n_points[child_id],
                   is_leaf[child_id]);
            log_trace("    pts_extents = %f %f %f %f %f %f\n",
                   extents[6*child_id+0],
                   extents[6*child_id+1],
                   extents[6*child_id+2],
                   extents[6*child_id+3],
                   extents[6*child_id+4],
                   extents[6*child_id+5]);
          }

          intersect = _intersect_node_box_explicit (3,
                                                    extents + 6*child_id,
                                                    _box_extents,
                                                    &node_inside_box);

          if (dbg_enabled) {
            log_trace("    intersect = %d\n", intersect);
          }

          if (intersect) {
            if (node_inside_box) {
              /* The box must contain all points */
              if (dbg_enabled) {
                log_trace("    add pts with lnum %d through %d\n", range[2*child_id+0], range[2*child_id+1]);
              }

              int new_size = _box_pts_idx[ibox+1] + n_points[child_id];

              if (tmp_size <= new_size) {
                tmp_size = PDM_MAX (2*tmp_size, new_size);
                *box_pts = realloc (*box_pts, sizeof(int) * tmp_size);
                _box_pts = *box_pts;
              }

              for (int j = 0; j < n_points[child_id]; j++) {
                _box_pts[_box_pts_idx[ibox+1]++] = new_to_old[range[2*child_id] + j];
              }
            }

            else {
              /* Push child in stack */
              stack_id[pos_stack++] = child_id;
            }
          }
        } // End of loop on children
      }

    } // End of while loop

  } // End of loop on boxes

  free (stack_id);
  *box_pts = realloc (*box_pts, sizeof(int) * _box_pts_idx[n_box]);
}



void
PDM_point_tree_seq_write_nodes_shared
(
       PDM_point_tree_seq_shm_t *shm_ptree,
 const int                       i_shm_rank,
 const char                     *filename
 )
{
  // _l_nodes_t *nodes = shm_ptree->nodes;
  int     n_nodes     = shm_ptree->shm_n_nodes    [i_shm_rank];
  int    *is_leaf     = shm_ptree->shm_is_leaf    [i_shm_rank];
  // int    *children_id = shm_ptree->shm_children_id[i_shm_rank];
  int    *range       = shm_ptree->shm_range      [i_shm_rank];
  // int    *n_points    = shm_ptree->shm_n_points   [i_shm_rank];
  double *extents     = shm_ptree->shm_extents    [i_shm_rank];
  // double *pts_coord   = shm_ptree->shm_pts_coord  [i_shm_rank];
  // int    *new_to_old  = shm_ptree->shm_new_to_old [i_shm_rank];

  // PDM_log_trace_array_int(shm_ptree->shm_is_leaf[i_shm_rank], n_nodes, " *shm_is_leaf : ");

  // double tol_visu = 1e-3;
  // double _ext[6];

  // write VTK
  FILE *f = fopen(filename, "w");

  fprintf(f, "# vtk DataFile Version 2.0\n");
  fprintf(f, "ptree_seq_shared\n");
  fprintf(f, "ASCII\n");
  fprintf(f, "DATASET UNSTRUCTURED_GRID\n");

  fprintf(f, "POINTS %d double\n", 8*n_nodes);
  for (int inode = 0; inode < n_nodes; inode++) {
    double *ext = extents + 6*inode;

    // if (1 && (range[2*inode+1] - range[2*inode] == 1)) {
    //   // Trick to visualize nodes with degenerate extents (single point)
    //   ext = _ext;
    //   for (int i = 0; i < 3; i++) {
    //     double x = nodes->extents[6*inode + i];
    //     double eps = tol_visu*(ptree->extents[i+3] - ptree->extents[i]);
    //     ext[i  ] = x - eps;
    //     ext[i+3] = x + eps;
    //   }
    // }

    for (int k = 0; k < 2; k++) {
      for (int j = 0; j < 2; j++) {
        for (int i = 0; i < 2; i++) {
          int ii = (1-j)*i + j*(1-i);
          fprintf(f, "%f %f %f\n", ext[3*ii], ext[3*j+1], ext[3*k+2]);
        }
      }
    }
  }

  fprintf(f, "CELLS %d %d\n", n_nodes, 9*n_nodes);
  for (int inode = 0; inode < n_nodes; inode++) {
    fprintf(f, "8 ");
    for (int j = 0; j < 8; j++) {
      fprintf(f, "%d ", 8*inode+j);
    }
    fprintf(f, "\n");
  }

  fprintf(f, "CELL_TYPES %d\n", n_nodes);
  for (int i = 0; i < n_nodes; i++) {
    fprintf(f, "%d\n", 12);
  }

  fprintf(f, "CELL_DATA %d\n", n_nodes);

  fprintf(f, "FIELD node_field 2\n");
  // fprintf(f, "depth 1 %d int\n", n_nodes);
  // for (int i = 0; i < n_nodes; i++) {
  //   fprintf(f, "%d\n", depth[i]);
  // }
  fprintf(f, "is_leaf 1 %d int\n", n_nodes);
  for (int i = 0; i < n_nodes; i++) {
    fprintf(f, "%d\n", is_leaf[i]);
  }
  fprintf(f, "n_pts 1 %d int\n", n_nodes);
  for (int i = 0; i < n_nodes; i++) {
    fprintf(f, "%d\n", range[2*i+1] - range[2*i]);
  }

  fclose(f);
}


void
PDM_point_tree_seq_closest_point_shared
(
       PDM_point_tree_seq_shm_t *shm_ptree,
 const int                       i_shm_rank,
 const int                       n_pts,
       double                   *pts_coord,
       int                      *closest_point_id,
       double                   *closest_point_dist2
 )
{
  // const int n_children = PDM_point_tree_n_children_get(ptree);
  int n_children = 2;
  if (shm_ptree->tree_type == PDM_DOCTREE_LOCAL_TREE_KDTREE) {
    n_children = 2;
  }
  else if (shm_ptree->tree_type == PDM_DOCTREE_LOCAL_TREE_OCTREE) {
    n_children = 8;
  }
  else {
    PDM_error(__FILE__, __LINE__, 0,
              "Invalid tree_type %d\n", (int) shm_ptree->tree_type);
  }

  const int depth_max = 31;
  // int     n_nodes     = shm_ptree->shm_n_nodes    [i_shm_rank];
  int    *is_leaf     = shm_ptree->shm_is_leaf    [i_shm_rank];
  int    *children_id = shm_ptree->shm_children_id[i_shm_rank];
  int    *range       = shm_ptree->shm_range      [i_shm_rank];
  int    *n_points    = shm_ptree->shm_n_points   [i_shm_rank];
  double *extents     = shm_ptree->shm_extents    [i_shm_rank];
  double *_pts_coord  = shm_ptree->shm_pts_coord  [i_shm_rank];
  int    *new_to_old  = shm_ptree->shm_new_to_old [i_shm_rank];

  int s_pt_stack = ((n_children - 1) * (depth_max - 1) + n_children);
  int    *stack           = malloc(sizeof(int   ) * s_pt_stack);
  int    *inbox_stack     = malloc(sizeof(int   ) * s_pt_stack);
  double *min_dist2_stack = malloc(sizeof(double) * s_pt_stack);

  int    sort_child[8];
  double dist_child[8];
  int    inbox_child[8];

  const int dim = 3;

  for (int i = 0; i < n_pts; i++) {

    int pos_stack = 0;
    const double *_pt = pts_coord + dim*i;

    /* Init stack */
    closest_point_id   [i] = -1;
    closest_point_dist2[i] = HUGE_VAL;


    double _min_dist2;
    int inbox1 = _box_dist2_min(dim,
                                &extents[0],
                                _pt,
                                &_min_dist2);

    stack          [pos_stack] = 0;
    inbox_stack    [pos_stack] = inbox1;
    min_dist2_stack[pos_stack] = _min_dist2;
    pos_stack++;

    while (pos_stack > 0) {

      pos_stack--;
      int    node_id   = stack          [pos_stack];
      double min_dist2 = min_dist2_stack[pos_stack];
      int    inbox     = inbox_stack    [pos_stack];

      if ((min_dist2 > closest_point_dist2[i]) || (inbox == 0)) {
        continue;
      }

      if (is_leaf[node_id]) {
        /* Leaf node */
        for (int ipt = range[2*node_id]; ipt < range[2*node_id+1]; ipt++) {
          double *_coords = _pts_coord + 3*ipt;

          double point_dist2 = 0;
          for (int j = 0; j < dim; j++) {
            point_dist2 +=
            (_coords[j] - _pt[j]) *
            (_coords[j] - _pt[j]);
          }

          if (point_dist2 < closest_point_dist2[i]) {
            closest_point_id   [i] = new_to_old[ipt];
            closest_point_dist2[i] = point_dist2;
          }
        } // End of loop on current leaf's points
      }

      else {
        /* Internal node */
        /* Sort children */
        const int *_child_ids = children_id + n_children*node_id;

        for (int ichild = 0; ichild < n_children; ichild++) {
          dist_child[ichild] = HUGE_VAL;
        }

        int n_selec = 0;
        for (int ichild = 0; ichild < n_children; ichild++) {

          int child_id = _child_ids[ichild];

          if (child_id < 0) {
            continue;
          }

          double child_min_dist2;
          int child_inbox = _box_dist2_min(dim,
                                           &extents[6*child_id],
                                           _pt,
                                           &child_min_dist2);

          int i1 = 0;
          for (i1 = n_selec;
               (i1 > 0) && (dist_child[i1-1] > child_min_dist2) ; i1--) {
            dist_child [i1] = dist_child[i1-1];
            sort_child [i1] = sort_child[i1-1];
            inbox_child[i1] = inbox_child[i1-1];
          }

          sort_child [i1] = child_id;
          dist_child [i1] = child_min_dist2;
          inbox_child[i1] = child_inbox;

          n_selec++;

        } // End of loop on current node's children


        /* Push selected children into the stack */
        for (int j = 0; j < n_selec; j++) {
          int j1 = n_selec - 1 - j;
          int child_id = sort_child[j1];
          if (child_id != -1) {
            if ((dist_child[j1] < closest_point_dist2[i]) &&
                (n_points[child_id] > 0)) {

              min_dist2_stack[pos_stack] = dist_child[j1];
              inbox_stack    [pos_stack] = inbox_child[j1];

              stack[pos_stack++] = child_id; /* push root in th stack */
            }
          }
        } // End of loop on selected children

      }


    } // End of while loop

  } // End of loop on target points
}


inline static int
_point_inside_box
(
 const int              dim,
 const double *restrict extents,
 const double *restrict coords
 )
{
  for (int i = 0; i < dim; i++) {
    if (coords[i] > extents[i+dim] || coords[i] < extents[i]) {
      return 0;
    }
  }

  return 1;
}

static int
_binary_search
(
 const int  elem,
 const int *array,
 const int  n
 )
{
  int l = 0;
  int r = n;

  if (n < 1)
    return 0;

  while (l + 1 < r) {
    int m = l + (r - l)/2;

    if (elem < array[m])
      r = m;
    else
      l = m;
  }

  if (array[l] < elem)
    return l + 1;
  else
    return l;
}

static inline void
_insertion_sort
(
 const int   point_id,
 int        *box_pts_n,
 int        *box_pts_s,
 int       **box_pts
 )
{
  int dbg = 0;

  if ((*box_pts_n) == 0) {
    (*box_pts)[(*box_pts_n)++] = point_id;
  }
  else {
    int i = _binary_search(point_id,
                           *box_pts,
                           *box_pts_n);
    if (dbg) {
      log_trace("%d at pos %d in array ", point_id, i);
      PDM_log_trace_array_int((*box_pts), (*box_pts_n), "");
    }

    if ((*box_pts)[i] == point_id) {
      return;
    }

    if ((*box_pts_s) <= (*box_pts_n)) {
      (*box_pts_s) *= 2;
      (*box_pts)    = realloc((*box_pts),
                              sizeof(int) * (*box_pts_s));
    }

    for (int j = (*box_pts_n); j > i; j--) {
      (*box_pts)[j] = (*box_pts)[j-1];
    }

    (*box_pts)[i] = point_id;
    (*box_pts_n)++;

    if (dbg) {
      PDM_log_trace_array_int(*box_pts, *box_pts_n, "after insertion : ");
    }
  }
}


typedef enum {
  SUBDIVISION_CRITERION_VOLUME,
  SUBDIVISION_CRITERION_LENGTH,
  SUBDIVISION_CRITERION_DEPTH,
} _subdivision_criterion_t;

void
PDM_tree_intersection_point_box2
(
 PDM_point_tree_seq_t  *btree, // Really a box tree
 PDM_point_tree_seq_t  *ptree,
 double                *box_extents,
 int                  **box_pts_idx,
 int                  **box_pts
 )
{
  _subdivision_criterion_t subdiv_crit = SUBDIVISION_CRITERION_VOLUME;

  int dbg  = 0;
  // int visu = 0;

  // PDM_boxes_t *boxes;
  // PDM_box_tree_data_t *box_tree_data;

  // boxes         = btree->boxes->local_boxes;
  // box_tree_data = btree->local_data;

  // int n_boxes = boxes->n_boxes;

  // double *btree_s, *btree_d;
  // PDM_box_set_normalization_get((PDM_box_set_t *) btree->boxes,
  //                               &btree_s,
  //                               &btree_d);

  /* Get point_tree data (use gets!!!) */
  int n_boxes = btree->n_pts;
  int btree_n_children = PDM_point_tree_n_children_get(btree);
  int    *btree_depth       = btree->nodes->depth;
  int    *btree_is_leaf     = btree->nodes->is_leaf;
  int    *btree_range       = btree->nodes->range;
  int    *btree_children_id = btree->nodes->children_id;
  double *btree_extents     = btree->nodes->extents;

  int *btree_new_to_old = NULL;
  PDM_point_tree_seq_point_new_to_old_get(btree,
                                          &btree_new_to_old);



  /* Get point_tree data (use gets!!!) */
  int ptree_n_children = PDM_point_tree_n_children_get(ptree);
  int    *ptree_depth       = ptree->nodes->depth;
  int    *ptree_is_leaf     = ptree->nodes->is_leaf;
  int    *ptree_range       = ptree->nodes->range;
  int    *ptree_children_id = ptree->nodes->children_id;
  double *ptree_extents     = ptree->nodes->extents;

  // int n_pts = ptree->n_pts;
  double *ptree_pts_coord;
  PDM_point_tree_seq_sorted_points_get(ptree,
                                       &ptree_pts_coord);

  double *_pts_coord = ptree_pts_coord;

  int *ptree_new_to_old = NULL;
  PDM_point_tree_seq_point_new_to_old_get(ptree,
                                          &ptree_new_to_old);



  // double btree_extents[6];

  /* Start from both roots */
  int btree_node_id = 0;
  int ptree_node_id = 0;

  int intersect = _intersect_box_box(3,
                                     &btree_extents[6*btree_node_id],
                                     &ptree_extents[6*ptree_node_id]);

  if (!intersect) {
    *box_pts_idx = PDM_array_zeros_int(n_boxes + 1);
    if (dbg) {
      log_trace("roots do not intersect\n");
    }
    return;
  }

  int s_queue = 1000; // ?
  int *queue0 = malloc(sizeof(int) * s_queue * 2);
  int *queue1 = malloc(sizeof(int) * s_queue * 2);
  int *queues[2] = {queue0, queue1};

  int n_queue = 0;
  queues[0][2*n_queue  ] = btree_node_id;
  queues[0][2*n_queue+1] = ptree_node_id;
  n_queue++;



  int  *__box_pts_n = PDM_array_zeros_int(n_boxes);
  int  *__box_pts_s = PDM_array_const_int(n_boxes, 4);
  int **__box_pts   = malloc(sizeof(int *) * n_boxes);
  for (int i = 0; i < n_boxes; i++) {
    __box_pts[i] = malloc(sizeof(int) * __box_pts_s[i]);
  }


  int istep = -1;
  while (n_queue > 0) {

    int new_n_queue = 0;

    if (dbg) {
      PDM_log_trace_array_int(queues[0], 2*n_queue, "queue : ");
    }

    for (int ipair = 0; ipair < n_queue; ipair++) {

      istep++;

      btree_node_id = queues[0][2*ipair  ];
      ptree_node_id = queues[0][2*ipair+1];
      // if (subdiv_crit == SUBDIVISION_CRITERION_DEPTH) {
      //   btree_depth = queues_depth[0][ipair];
      // }

      if (dbg) {
        log_trace("  Step %d\n", istep);
        log_trace("  node ids %d %d\n", btree_node_id, ptree_node_id);
      }

      // if (visu) {
      //   char filename[999];
      //   sprintf(filename, "intersection_btree_ptree_step_%4.4d.vtk", istep);
      //   _visu_pair(filename,
      //              btree,
      //              ptree,
      //              btree_node_id,
      //              ptree_node_id);
      // }

      // _node_t *node = &(box_tree_data->nodes[btree_node_id]);

      int isubdiv;

      /* Which tree do we subdivide? */
      if (btree_is_leaf[btree_node_id]) {

        if (ptree_is_leaf[ptree_node_id]) {
          /* Both leaves */
          isubdiv = -1;

          if (dbg) {
            log_trace("  both leaves\n");
          }

          /* inspect boxes contained in current leaf node */
          //for (int ibox = 0; ibox < node->n_boxes; ibox++) {
          for (int ibox = btree_range[2*btree_node_id]; ibox < btree_range[2*btree_node_id+1]; ibox++) {

            // int box_id = box_tree_data->box_ids[node->start_id + ibox];
            int box_id = btree_new_to_old[ibox];

            if (dbg) {
              log_trace("    box_id = %d\n", box_id);
            }

            // const double *box_extents = boxes->extents + box_id*6;
            double *_box_extents = box_extents + box_id*6;

            for (int ipt = ptree_range[2*ptree_node_id]; ipt < ptree_range[2*ptree_node_id+1]; ipt++) {

              double *pt = _pts_coord + 3*ipt;
              int point_id = ptree_new_to_old[ipt];

              if (dbg) {
                log_trace("      point_id = %d (%f %f %f)\n",
                          point_id, pt[0], pt[1], pt[2]);
              }

              if (_point_inside_box(3, _box_extents, pt)) {
                if (dbg) {
                  log_trace("        inside box\n");
                }

                _insertion_sort(point_id,
                                &__box_pts_n[box_id],
                                &__box_pts_s[box_id],
                                &__box_pts  [box_id]);
              }

            } // End of loop on current ptree leaf's points
          } // End of loop on current btree leaf's boxes


        }
        else {
          /* Subdivide point tree */
          isubdiv = 1;
        }

      }

      else if (ptree_is_leaf[ptree_node_id]) {
        /* Subdivide box tree */
        isubdiv = 0;
      }

      else {

        // Decide which tree is subdivided
        // _extents_real(3,
        //               box_tree_data->nodes[btree_node_id].morton_code,
        //               btree_s,
        //               btree_d,
        //               btree_extents);

        double btree_crit = 1.;
        double ptree_crit = 1.;
        switch (subdiv_crit) {

          case SUBDIVISION_CRITERION_VOLUME: {
            btree_crit = 1.;
            ptree_crit = 1.;
            for (int i = 0; i < 3; i++) {
              btree_crit *= (btree_extents[6*btree_node_id + i+3] - btree_extents[6*btree_node_id + i]);
              ptree_crit *= (ptree_extents[6*ptree_node_id + i+3] - ptree_extents[6*ptree_node_id + i]);
            }
            break;
          }
          case SUBDIVISION_CRITERION_LENGTH: {
            btree_crit = 0.;
            ptree_crit = 0.;
            for (int i = 0; i < 3; i++) {
              btree_crit = PDM_MAX(btree_crit, (btree_extents[6*btree_node_id + i+3] - btree_extents[6*btree_node_id + i]));
              ptree_crit = PDM_MAX(ptree_crit, (ptree_extents[6*ptree_node_id + i+3] - ptree_extents[6*ptree_node_id + i]));
            }
            break;
          }
          case SUBDIVISION_CRITERION_DEPTH: {
            btree_crit = btree_depth[btree_node_id];
            ptree_crit = ptree_depth[ptree_node_id];
            break;
          }
          default: {
            PDM_error(__FILE__, __LINE__, 0,
                      "Subdivision criterion %d not implemented\n", (int) subdiv_crit);
            break;
          }
        }

        if (btree_crit > ptree_crit) {
          /* Subdivide box tree */
          isubdiv = 0;
        }
        else {
          /* Subdivide point tree */
          isubdiv = 1;
        }

      }


      /* Add children to new queue */
      if (isubdiv == 0) {
        /* Subdivide box tree */
        if (dbg) {
          log_trace("  subdivide box tree\n");
        }

        // int *children_id = box_tree_data->child_ids + btree_node_id*btree->n_children;
        int *children_id = btree_children_id + btree_node_id*btree_n_children;
        // for (int ichild = 0; ichild < btree->n_children; ichild++) {
        for (int ichild = 0; ichild < btree_n_children; ichild++) {
          int child_id = children_id[ichild];

          if (child_id < 0) continue;

          // _extents_real(3,
          //               box_tree_data->nodes[child_id].morton_code,
          //               btree_s,
          //               btree_d,
          //               btree_extents);
          // log_trace("btree_extents : %f %f %f  %f %f %f",
          //           btree_extents[0], btree_extents[1], btree_extents[2],
          //           btree_extents[3], btree_extents[4], btree_extents[5]);

          intersect = _intersect_box_box(3,
                                         &btree_extents[6*child_id],
                                         &ptree_extents[6*ptree_node_id]);

          if (dbg) {
            log_trace("    child %d, intersect? %d\n", child_id, intersect);
          }

          if (intersect) {
            // Check size!!!
            if (new_n_queue >= s_queue) {
              s_queue *= 2;
              queues[0] = realloc(queues[0], sizeof(int) * s_queue * 2);
              queues[1] = realloc(queues[1], sizeof(int) * s_queue * 2);
            }

            queues[1][2*new_n_queue  ] = child_id;
            queues[1][2*new_n_queue+1] = ptree_node_id;

            new_n_queue++;
          }
        }

      }
      else if (isubdiv == 1) {
        /* Subdivide point tree */
        if (dbg) {
          log_trace("  subdivide point tree\n");
        }


        int *children_id = ptree_children_id + ptree_node_id*ptree_n_children;
        for (int ichild = 0; ichild < ptree_n_children; ichild++) {
          int child_id = children_id[ichild];

          if (child_id < 0) continue;

          intersect = _intersect_box_box(3,
                                         &btree_extents[6*btree_node_id],
                                         &ptree_extents[6*child_id]);

          if (dbg) {
            log_trace("    child %d, intersect? %d\n", child_id, intersect);
          }

          if (intersect) {
            // Check size!!!
            if (new_n_queue >= s_queue) {
              s_queue *= 2;
              queues[0] = realloc(queues[0], sizeof(int) * s_queue * 2);
              queues[1] = realloc(queues[1], sizeof(int) * s_queue * 2);
            }

            queues[1][2*new_n_queue  ] = btree_node_id;
            queues[1][2*new_n_queue+1] = child_id;

            new_n_queue++;
          }
        }

      }


    } // End of loop in current queue

    if (dbg) {
      PDM_log_trace_array_int(queues[1], 2*new_n_queue, "new_queue : ");
    }

    /* Swap queues */
    int *tmp = queues[1];
    queues[1] = queues[0];
    queues[0] = tmp;

    n_queue = new_n_queue;

  } // End of while loop
  free(__box_pts_s);
  // free(_pts_coord);
  free(queues[0]);
  free(queues[1]);

  /* Re-arrange result */
  *box_pts_idx = PDM_array_new_idx_from_sizes_int(__box_pts_n, n_boxes);

  *box_pts = malloc(sizeof(int) * (*box_pts_idx)[n_boxes]);
  for (int i = 0; i < n_boxes; i++) {
    int *bp = *box_pts + (*box_pts_idx)[i];

    for (int j = 0; j < __box_pts_n[i]; j++) {
      bp[j] = __box_pts[i][j];
    }

    free(__box_pts[i]);
  }
  free(__box_pts_n);
  free(__box_pts);

  free(queue0);
  free(queue1);

}




void
PDM_point_tree_seq_intersect_box_leaf
(
       PDM_point_tree_seq_t  *ptree,
 const int                    n_box,
 const double                 box_extents[],
       int                  **box_leaf_idx,
       int                  **box_leaf
 )
{
  *box_leaf_idx = malloc (sizeof(int) * (n_box + 1));
  int *_box_leaf_idx = *box_leaf_idx;
  _box_leaf_idx[0] = 0;

  if (n_box < 1) {
    *box_leaf = malloc (sizeof(int) * _box_leaf_idx[n_box]);
    return;
  }

  const int n_children = PDM_point_tree_n_children_get(ptree);


  _l_nodes_t *nodes = ptree->nodes;

  int s_pt_stack = ((n_children - 1) * (ptree->depth_max - 1) + n_children);
  int *stack_id  = malloc (s_pt_stack * sizeof(int));

  int node_inside_box;
  int intersect;

  int tmp_size = 4 * n_box;
  *box_leaf = malloc (sizeof(int) * tmp_size);
  int *_box_leaf = *box_leaf;

  for (int ibox = 0; ibox < n_box; ibox++) {
    int dbg_enabled = 0;
    if (dbg_enabled) {
      log_trace("box %d\n", ibox);
    }

    _box_leaf_idx[ibox+1] = _box_leaf_idx[ibox];

    const double *_box_extents = box_extents + 6*ibox;

    intersect = _intersect_node_box_explicit (3,
                                              &nodes->extents[0],
                                              _box_extents,
                                              &node_inside_box);

    if (!intersect) {
      continue;
    }


    if (nodes->is_leaf[0]) {
      if (tmp_size <= _box_leaf_idx[ibox+1]+1) {
        tmp_size = PDM_MAX (2*tmp_size, _box_leaf_idx[ibox+1]+1);
        *box_leaf = realloc (*box_leaf, sizeof(int) * tmp_size);
        _box_leaf = *box_leaf;
      }
      _box_leaf[_box_leaf_idx[ibox+1]++] = 0;
      continue;
    }


    /* Push root in stack */
    int pos_stack = 0;
    stack_id[pos_stack++] = 0;

    while (pos_stack > 0) {
      int node_id = stack_id[--pos_stack];

      if (dbg_enabled) {
        log_trace("  node %d, range=%d/%d, n_points=%d, leaf_id=%d\n",
                  node_id,
                  nodes->range[2*node_id], nodes->range[2*node_id+1],
                  nodes->n_points[node_id],
                  nodes->is_leaf[node_id]);
      }

      const int *_child_ids = nodes->children_id + n_children*node_id;
      for (int i = 0; i < n_children; i++) {
        int child_id = _child_ids[i];
        if (child_id < 0) {
          continue;
        }

        if (dbg_enabled) {
          log_trace("    child %d: id=%d, range=%d/%d, n_points=%d, leaf_id=%d\n",
                    i,
                    child_id,
                    nodes->range[2*child_id+0],
                    nodes->range[2*child_id+1],
                    nodes->n_points[child_id],
                    nodes->is_leaf[child_id]);
          log_trace("    leaf_extents = %f %f %f %f %f %f\n",
                    nodes->extents[6*child_id+0],
                    nodes->extents[6*child_id+1],
                    nodes->extents[6*child_id+2],
                    nodes->extents[6*child_id+3],
                    nodes->extents[6*child_id+4],
                    nodes->extents[6*child_id+5]);
        }

        intersect = _intersect_node_box_explicit (3,
                                                  nodes->extents + 6*child_id,
                                                  _box_extents,
                                                  &node_inside_box);

        if (dbg_enabled) {
          log_trace("    intersect = %d\n", intersect);
        }

        if (intersect) {
          if (nodes->is_leaf[child_id]) {
            if (tmp_size <= _box_leaf_idx[ibox+1]+1) {
              tmp_size = PDM_MAX(2*tmp_size, _box_leaf_idx[ibox+1]+1);
              *box_leaf = realloc(*box_leaf, sizeof(int) * tmp_size);
              _box_leaf = *box_leaf;
            }
            _box_leaf[_box_leaf_idx[ibox+1]++] = child_id;
          }
          else {
            /* Push child in stack */
            stack_id[pos_stack++] = child_id;
          }
        }

      } // End of loop on children

    } /* End While */
  } /* End boxe loop */

  free (stack_id);
  *box_leaf = realloc(*box_leaf, sizeof(int) * _box_leaf_idx[n_box]);

}
