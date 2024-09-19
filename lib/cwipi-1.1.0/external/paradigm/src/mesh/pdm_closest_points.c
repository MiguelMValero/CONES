/*
  This file is part of the ParaDiGM library.

  Copyright (C) 2019       ONERA

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
/*----------------------------------------------------------------------------
 * System headers
 *----------------------------------------------------------------------------*/

#include <assert.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "pdm_printf.h"
#include "pdm_error.h"
#include "pdm.h"
#include "pdm_priv.h"
#include "pdm_mpi.h"
#include "pdm_timer.h"
#include "pdm_closest_points.h"
#include "pdm_closest_points_priv.h"
#include "pdm_para_octree.h"
#include "pdm_part_to_block.h"
#include "pdm_part_to_part.h"
#include "pdm_block_to_part.h"
#include "pdm_array.h"
#include "pdm_sort.h"
#include "pdm_binary_search.h"
#include "pdm_logging.h"

/*----------------------------------------------------------------------------*/

#ifdef	__cplusplus
extern "C" {
#if 0
} /* Fake brace to force Emacs auto-indentation back to column 0 */
#endif
#endif

/*============================================================================
 * Macro definitions
 *============================================================================*/


/*============================================================================
 * Type definitions
 *============================================================================*/

#define NTIMER 2

/**
 * \enum _timer_step_t
 *
 */

typedef enum {

  BEGIN    = 0,
  END      = 1,

} _timer_step_t;

/*============================================================================
 * Global variable
 *============================================================================*/

//static int idebug = 0;

/*=============================================================================
 * Private function definitions
 *============================================================================*/

/**
 *
 * \brief Reverse result
 *
 * \param [in]   cls                   Pointer to \ref PDM_closest_points object
 *
 */
static void
_closest_points_reverse_results
(
 PDM_closest_point_t  *cls
)
{

  assert (cls->tgt_cloud->closest_src_gnum != NULL);
  assert (cls->tgt_cloud->closest_src_dist != NULL);

  int* n_points = (int * ) malloc( cls->tgt_cloud->n_part * sizeof(int));
  PDM_g_num_t **tgt_g_num   = (PDM_g_num_t ** ) malloc( cls->tgt_cloud->n_part * sizeof(PDM_g_num_t *));
  int         **tgt_g_num_n = (int         ** ) malloc( cls->tgt_cloud->n_part * sizeof(int         *));
  for (int i_part = 0; i_part < cls->tgt_cloud->n_part; i_part++) {

    // if (1) {
    //   for (int i = 0; i < cls->tgt_cloud->n_points[i_part]; i++) {
    //     log_trace(PDM_FMT_G_NUM" : ", cls->tgt_cloud->gnum[i_part][i]);
    //     PDM_log_trace_array_long(cls->tgt_cloud->closest_src_gnum[i_part] + cls->n_closest*i,
    //                              cls->n_closest,
    //                              "");
    //   }
    // }

    n_points[i_part] = cls->tgt_cloud->n_points[i_part] * cls->n_closest;
    tgt_g_num  [i_part] = (PDM_g_num_t * ) malloc( n_points[i_part] * sizeof(PDM_g_num_t));
    tgt_g_num_n[i_part] = (int         * ) malloc( n_points[i_part] * sizeof(int        ));

    // PDM_log_trace_array_long(cls->tgt_cloud->closest_src_gnum[i_part], cls->tgt_cloud->n_points[i_part], "cls->tgt_cloud->closest_src_gnum:: " );

    for(int i = 0; i < cls->tgt_cloud->n_points[i_part]; ++i) {
      for(int ii = 0; ii < cls->n_closest; ++ii) {
        int idx = i * cls->n_closest + ii;
        tgt_g_num  [i_part][idx] = cls->tgt_cloud->gnum[i_part][i];
        tgt_g_num_n[i_part][idx] = 1;
      }
    }
  }

  /*
   * Compute the total number of target point to setup properly the part_to_block partial
   */
  PDM_g_num_t n_g_src = 0;
  for (int i_part = 0; i_part < cls->src_cloud->n_part; i_part++) {
    for(int i = 0; i < cls->src_cloud->n_points[i_part]; ++i) {
      n_g_src = PDM_MAX(n_g_src, cls->src_cloud->gnum[i_part][i]);
    }
  }
  PDM_g_num_t _n_g_src = 0;
  PDM_MPI_Allreduce (&n_g_src, &_n_g_src, 1, PDM__PDM_MPI_G_NUM, PDM_MPI_MAX, cls->comm);
  n_g_src = _n_g_src;

  if (0) {
    for (int i_part = 0; i_part < cls->tgt_cloud->n_part; i_part++) {
      PDM_log_trace_array_long(cls->tgt_cloud->closest_src_gnum[i_part],
                               n_points[i_part],
                               "closest_src_gnum : ");
    }

    for (int i_part = 0; i_part < cls->tgt_cloud->n_part; i_part++) {
      PDM_log_trace_array_double(cls->tgt_cloud->closest_src_dist[i_part],
                               n_points[i_part],
                               "closest_src_dist : ");
    }
  }

  /*
   *  First part to block to map in global numbering of src all target associate
   */
  PDM_part_to_block_t *ptb = PDM_part_to_block_create(PDM_PART_TO_BLOCK_DISTRIB_ALL_PROC,
                                                      PDM_PART_TO_BLOCK_POST_MERGE,
                                                      1.,
                                                      cls->tgt_cloud->closest_src_gnum,
                                                      NULL,
                                                      n_points,
                                                      cls->tgt_cloud->n_part,
                                                      cls->comm);

  /*
   * For each target (associated to a src gnum prepare an exhange of the current target g_num)
   */
  int *block_tgt_in_src_n = NULL;
  PDM_g_num_t *block_tgt_in_src_g_num = NULL;
  int blk_size = PDM_part_to_block_exch (ptb,
                                         sizeof(PDM_g_num_t),
                                         PDM_STRIDE_VAR_INTERLACED,
                                        1,
                                        tgt_g_num_n,
                              (void **) tgt_g_num,
                                        &block_tgt_in_src_n,
                              (void **) &block_tgt_in_src_g_num);

  if(0 == 1) {
    int block_n_elt = PDM_part_to_block_n_elt_block_get (ptb);
    PDM_log_trace_array_int(block_tgt_in_src_n     , block_n_elt, "block_tgt_in_src_n:: " );
    PDM_log_trace_array_long(block_tgt_in_src_g_num, blk_size   , "block_tgt_in_src_g_num:: " );
  }

  int *block_src_dist_n = NULL;
  double *block_tgt_in_src_dist = NULL;
  blk_size = PDM_part_to_block_exch(ptb,
                                    sizeof(double),
                                    PDM_STRIDE_VAR_INTERLACED,
                                    1,
                                    tgt_g_num_n,
                          (void **) cls->tgt_cloud->closest_src_dist,
                                    &block_src_dist_n,
                          (void **) &block_tgt_in_src_dist);
  free(block_src_dist_n); // Same than block_tgt_in_src_n

  for (int i_part = 0; i_part < cls->tgt_cloud->n_part; i_part++) {
    free(tgt_g_num  [i_part]);
    free(tgt_g_num_n[i_part]);
  }
  free(tgt_g_num  );
  free(tgt_g_num_n);
  free(n_points);

  int i_rank;
  int n_rank;
  PDM_MPI_Comm_rank(cls->comm, &i_rank);
  PDM_MPI_Comm_size(cls->comm, &n_rank);

  PDM_g_num_t *_block_distrib_idx = PDM_part_to_block_adapt_partial_block_to_block(ptb,
                                                                                   &block_tgt_in_src_n, /* Realloc inside */
                                                                                   n_g_src);

  PDM_part_to_block_free(ptb);

  PDM_block_to_part_t *btp = PDM_block_to_part_create(_block_distrib_idx,
                               (const PDM_g_num_t **) cls->src_cloud->gnum,
                                                      cls->src_cloud->n_points,
                                                      cls->src_cloud->n_part,
                                                      cls->comm);


  int** tgt_in_src_n;
  PDM_block_to_part_exch(btp,
                          sizeof(PDM_g_num_t),
                          PDM_STRIDE_VAR_INTERLACED,
                          block_tgt_in_src_n,
                          block_tgt_in_src_g_num,
                         &tgt_in_src_n,
              (void ***) &cls->src_cloud->tgt_in_src);

  int** useless_stride;
  PDM_block_to_part_exch(btp,
                          sizeof(double),
                          PDM_STRIDE_VAR_INTERLACED,
                          block_tgt_in_src_n,
                          block_tgt_in_src_dist,
                         &useless_stride, //Same than tgt_in_src_n
              (void ***) &cls->src_cloud->tgt_in_src_dist);

  /* Suppression des doublons de gnum */

  for (int i_part = 0; i_part < cls->src_cloud->n_part; i_part++) {
    int  max_tgt_in_src_n = 0;
    for (int i_point = 0; i_point < cls->src_cloud->n_points[i_part]; i_point++){
      max_tgt_in_src_n = PDM_MAX( max_tgt_in_src_n, tgt_in_src_n[i_part][i_point] );
    }
    int  idx_read  = 0;
    int  idx_write = 0;
    int    *tgt_order = (int    *) malloc(sizeof(int   ) * max_tgt_in_src_n);
    double *tgt_dist  = (double *) malloc(sizeof(double) * max_tgt_in_src_n);
    for (int i_point = 0; i_point < cls->src_cloud->n_points[i_part]; i_point++){
      if (tgt_in_src_n[i_part][i_point] > 0){

        /* Tri des gnum */
        int n_tgt = tgt_in_src_n[i_part][i_point];
        for (int i_tgt = 0; i_tgt < n_tgt; i_tgt++){
          tgt_order[i_tgt] = i_tgt;
          tgt_dist [i_tgt] = cls->src_cloud->tgt_in_src_dist[i_part][idx_read+i_tgt];
        }
        PDM_sort_long(&cls->src_cloud->tgt_in_src[i_part][idx_read], tgt_order, n_tgt);

        /* Ajout du premier point */
        tgt_in_src_n[i_part][i_point] = 1;
        cls->src_cloud->tgt_in_src     [i_part][idx_write] = cls->src_cloud->tgt_in_src[i_part][idx_read];
        cls->src_cloud->tgt_in_src_dist[i_part][idx_write] = tgt_dist[tgt_order[0]];
        idx_write++;

        /* Ajout des points suivants sans doublon */
        for (int i_tgt = 1; i_tgt < n_tgt; i_tgt++){
          if (cls->src_cloud->tgt_in_src[i_part][idx_read+i_tgt] > cls->src_cloud->tgt_in_src[i_part][idx_read+i_tgt-1]){
            cls->src_cloud->tgt_in_src     [i_part][idx_write] = cls->src_cloud->tgt_in_src[i_part][idx_read+i_tgt];
            cls->src_cloud->tgt_in_src_dist[i_part][idx_write] = tgt_dist[tgt_order[i_tgt]];
            tgt_in_src_n[i_part][i_point]++;
            idx_write++;
          }
        }

        idx_read = idx_read + n_tgt;
      }
    }
    free(tgt_order);
    free(tgt_dist);
  }

  cls->src_cloud->tgt_in_src_idx = (int **) malloc( cls->src_cloud->n_part * sizeof(int *));
  for (int i_part = 0; i_part < cls->src_cloud->n_part; i_part++) {
    // PDM_log_trace_array_int(tgt_in_src_n[i_part]     , cls->src_cloud->n_points[i_part], "cls->src_cloud->n_points[i_part]:: " );
    cls->src_cloud->tgt_in_src_idx[i_part] = PDM_array_new_idx_from_sizes_int(tgt_in_src_n[i_part], cls->src_cloud->n_points[i_part]);
    /* Réallocation à la bonne taille sans doublon */
    cls->src_cloud->tgt_in_src     [i_part] = realloc(cls->src_cloud->tgt_in_src     [i_part], sizeof(PDM_g_num_t) * cls->src_cloud->tgt_in_src_idx[i_part][cls->src_cloud->n_points[i_part]]);
    cls->src_cloud->tgt_in_src_dist[i_part] = realloc(cls->src_cloud->tgt_in_src_dist[i_part], sizeof(double)      * cls->src_cloud->tgt_in_src_idx[i_part][cls->src_cloud->n_points[i_part]]);
    free(tgt_in_src_n[i_part]);
    free(useless_stride[i_part]);
  }
  free(tgt_in_src_n);
  free(useless_stride);

  PDM_block_to_part_free(btp);
  free(_block_distrib_idx);
  free(block_tgt_in_src_n);
  free(block_tgt_in_src_g_num);
  free(block_tgt_in_src_dist);
}


/*============================================================================
 * Public function definitions
 *============================================================================*/

/**
 *
 * \brief Create a structure to look for the closest points of a point cloud
 * (target cloud) in an other point cloud (source cloud)
 *
 * \param [in]   comm           MPI communicator
 * \param [in]   n_closest      Number of closest source points to find for each
 *                              target point
 * \param [in ] owner           Ownership for \ref PDM_closest_point_t
 *
 * \return     Pointer to \ref PDM_closest_points object
 *
 */

PDM_closest_point_t*
PDM_closest_points_create
(
 const PDM_MPI_Comm    comm,
 const int             n_closest,
 const PDM_ownership_t owner
)
{
  PDM_closest_point_t *closest = (PDM_closest_point_t *) malloc(sizeof(PDM_closest_point_t));

  closest->comm                         = comm;
  closest->owner                        = owner;
  closest->results_is_getted            = PDM_FALSE;
  closest->tgt_in_src_results_is_getted = PDM_FALSE;
  closest->tgt_in_src_results_is_getted_d = PDM_FALSE;

  closest->n_closest = n_closest;
  closest->src_cloud = NULL;
  closest->tgt_cloud = NULL;

  closest->timer = PDM_timer_create ();

  for (int i = 0; i < NTIMER; i++) {
    closest->times_elapsed[i] = 0.;
    closest->times_cpu[i] = 0.;
    closest->times_cpu_u[i] = 0.;
    closest->times_cpu_s[i] = 0.;
  }

  closest->ptp = NULL;
  closest->ptp_ownership = PDM_OWNERSHIP_KEEP;

  return closest;
}



/**
 *
 * \brief Set the number of partitions of a point cloud
 *
 * \param [in]   cls               Pointer to \ref PDM_closest_points object
 * \param [in]   n_part_cloud_src  Number of partitions of the source cloud
 * \param [in]   n_part_cloud_tgt  Number of partitions of the target cloud
 *
 */

void
PDM_closest_points_n_part_cloud_set
(
       PDM_closest_point_t* cls,
 const int                  n_part_cloud_src,
 const int                  n_part_cloud_tgt
)
{
  assert(cls->src_cloud == NULL);
  assert(cls->tgt_cloud == NULL);

  cls->src_cloud = malloc (sizeof(_src_point_cloud_t));
  cls->tgt_cloud = malloc (sizeof(_tgt_point_cloud_t));

  cls->src_cloud->n_part   = n_part_cloud_src;
  cls->src_cloud->coords   = malloc (sizeof(double      *) * n_part_cloud_src);
  cls->src_cloud->gnum     = malloc (sizeof(PDM_g_num_t *) * n_part_cloud_src);
  cls->src_cloud->n_points = malloc (sizeof(int          ) * n_part_cloud_src);

  cls->src_cloud->tgt_in_src_idx = NULL;
  cls->src_cloud->tgt_in_src     = NULL;
  cls->src_cloud->tgt_in_src_dist= NULL;

  cls->tgt_cloud->n_part            = n_part_cloud_tgt;
  cls->tgt_cloud->coords            = malloc (sizeof(double      *) * n_part_cloud_tgt);
  cls->tgt_cloud->gnum              = malloc (sizeof(PDM_g_num_t *) * n_part_cloud_tgt);
  cls->tgt_cloud->n_points          = malloc (sizeof(int          ) * n_part_cloud_tgt);
  cls->tgt_cloud->closest_src_gnum  = NULL;
  cls->tgt_cloud->closest_src_dist  = NULL;
}


/**
 *
 * \brief Set the target point cloud
 *
 * \param [in]   cls             Pointer to \ref PDM_closest_points object
 * \param [in]   i_part          Index of partition
 * \param [in]   n_points        Number of points
 * \param [in]   coords          Point coordinates
 * \param [in]   gnum            Point global number
 *
 */

void
PDM_closest_points_tgt_cloud_set
(
       PDM_closest_point_t *cls,
 const int                  i_part,
 const int                  n_points,
       double              *coords,
       PDM_g_num_t         *gnum
)
{
  assert(cls->tgt_cloud != NULL);
  cls->tgt_cloud->n_points[i_part] = n_points;
  cls->tgt_cloud->coords[i_part] = coords;
  cls->tgt_cloud->gnum[i_part] = gnum;
}


/**
 *
 * \brief Set the source point cloud
 *
 * \param [in]   cls             Pointer to \ref PDM_closest_points object
 * \param [in]   i_part          Index of partition
 * \param [in]   n_points        Number of points
 * \param [in]   coords          Point coordinates
 * \param [in]   gnum            Point global number
 *
 */

void
PDM_closest_points_src_cloud_set
(
       PDM_closest_point_t *cls,
 const int                  i_part,
 const int                  n_points,
       double              *coords,
       PDM_g_num_t         *gnum
)
{
  assert(cls->src_cloud != NULL);
  cls->src_cloud->n_points[i_part] = n_points;
  cls->src_cloud->coords  [i_part] = coords;
  cls->src_cloud->gnum    [i_part] = gnum;
}

/**
 *
 * \brief Look for closest points
 *
 * \param [in]   cls Pointer to \ref PDM_closest_points object
 *
 */

void
PDM_closest_points_compute
(
PDM_closest_point_t *cls
)
{

  cls->times_elapsed[BEGIN] = PDM_timer_elapsed(cls->timer);
  cls->times_cpu[BEGIN]     = PDM_timer_cpu(cls->timer);
  cls->times_cpu_u[BEGIN]   = PDM_timer_cpu_user(cls->timer);
  cls->times_cpu_s[BEGIN]   = PDM_timer_cpu_sys(cls->timer);

  PDM_timer_resume(cls->timer);


  int i_rank;
  PDM_MPI_Comm_rank (cls->comm, &i_rank);
  int n_rank;
  PDM_MPI_Comm_rank (cls->comm, &n_rank);

  /*if (i_rank == 0) {
    printf(">> PDM_closest_points_compute\n");
    fflush(stdout);
    }*/

  /*
   *  Make sure we have at least as many source points as requested closest points
   */
  PDM_g_num_t ln_src_pts = 0;
  for (int i_part = 0; i_part < cls->src_cloud->n_part; i_part++) {
    ln_src_pts += cls->src_cloud->n_points[i_part];
  }

  PDM_g_num_t gn_src_pts;
  PDM_MPI_Allreduce (&ln_src_pts, &gn_src_pts, 1, PDM__PDM_MPI_G_NUM, PDM_MPI_SUM, cls->comm);

  assert (gn_src_pts >= cls->n_closest);



  const int depth_max = 31;
  const int points_in_leaf_max = cls->n_closest;
  const int build_leaf_neighbours = 0;


  /* Create empty parallel octree structure */
  PDM_para_octree_t *octree = PDM_para_octree_create (cls->src_cloud->n_part,
                                                      depth_max,
                                                      points_in_leaf_max,
                                                      build_leaf_neighbours,
                                                      cls->comm);


  /* Set source point clouds */
  for (int i_part = 0; i_part < cls->src_cloud->n_part; i_part++) {
    PDM_para_octree_point_cloud_set (octree,
                                     i_part,
                                     cls->src_cloud->n_points[i_part],
                                     cls->src_cloud->coords[i_part],
                                     cls->src_cloud->gnum[i_part]);
  }

  /* Build parallel octree */
  PDM_para_octree_build (octree, NULL);

  /*if (i_rank == 0) {
    printf("PDM_para_octree_build OK\n");
    fflush(stdout);
    }*/

  //PDM_para_octree_dump (octree);
  if (0) {
    PDM_para_octree_dump_times (octree);
  }
  //<--


  /* Concatenate partitions */
  int n_tgt = 0;
  for (int i_part = 0; i_part < cls->tgt_cloud->n_part; i_part++)
    n_tgt += cls->tgt_cloud->n_points[i_part];

  double      *tgt_coord = malloc (sizeof(double)      * n_tgt * 3);
  PDM_g_num_t *tgt_g_num = malloc (sizeof(PDM_g_num_t) * n_tgt);
  PDM_g_num_t *closest_src_gnum = malloc (sizeof(PDM_g_num_t) * n_tgt * cls->n_closest);
  double      *closest_src_dist = malloc (sizeof(double)      * n_tgt * cls->n_closest);

  n_tgt = 0;
  for (int i_part = 0; i_part < cls->tgt_cloud->n_part; i_part++) {
    for (int i = 0; i < cls->tgt_cloud->n_points[i_part]; i++) {
      for (int j = 0; j < 3; j++)
        tgt_coord[3*(n_tgt + i) + j] = cls->tgt_cloud->coords[i_part][3*i + j];
      tgt_g_num[n_tgt + i] = cls->tgt_cloud->gnum[i_part][i];
    }
    n_tgt += cls->tgt_cloud->n_points[i_part];
  }

  /* Search closest source points from target points */
  if (cls->n_closest == 1) {
    PDM_para_octree_single_closest_point (octree,
                                          n_tgt,
                                          tgt_coord,
                                          tgt_g_num,
                                          closest_src_gnum,
                                          closest_src_dist);
  } else {
    PDM_para_octree_closest_points (octree,
                                    cls->n_closest,
                                    n_tgt,
                                    tgt_coord,
                                    tgt_g_num,
                                    closest_src_gnum,
                                    closest_src_dist);
  }

  if (0) {
    for (int i = 0; i < n_tgt; i++) {
      log_trace(PDM_FMT_G_NUM" (%f %f %f) : ",
                tgt_g_num[i],
                tgt_coord[3*i], tgt_coord[3*i+1], tgt_coord[3*i+2]);
      PDM_log_trace_array_long(closest_src_gnum + cls->n_closest*i,
                               cls->n_closest,
                               "");
    }
  }


  // PDM_log_trace_array_long(tgt_g_num, n_tgt, "tgt_g_num:: " );
  // PDM_log_trace_array_double(tgt_coord, 3 * n_tgt, "tgt_coord:: " );
  // PDM_log_trace_array_long(closest_src_gnum, n_tgt * cls->n_closest, "closest_src_gnum:: " );
  // PDM_log_trace_array_double(closest_src_dist, n_tgt * cls->n_closest, "closest_src_dist:: " );

  /* Restore partitions */
  free (tgt_coord);
  free (tgt_g_num);
  n_tgt = 0;

  cls->tgt_cloud->closest_src_gnum = malloc (sizeof(PDM_g_num_t *) * cls->tgt_cloud->n_part);
  cls->tgt_cloud->closest_src_dist = malloc (sizeof(double *)      * cls->tgt_cloud->n_part);

  for (int i_part = 0; i_part < cls->tgt_cloud->n_part; i_part++) {
    int s_closest_src = cls->n_closest * cls->tgt_cloud->n_points[i_part];

    cls->tgt_cloud->closest_src_gnum[i_part] = malloc (sizeof(PDM_g_num_t) * s_closest_src);
    cls->tgt_cloud->closest_src_dist[i_part] = malloc (sizeof(double)      * s_closest_src);

    for (int i = 0; i < cls->tgt_cloud->n_points[i_part]; i++) {
      for (int j = 0; j < cls->n_closest; j++) {
        cls->tgt_cloud->closest_src_gnum[i_part][cls->n_closest*i+j] =
          closest_src_gnum[n_tgt + cls->n_closest*i + j];

        cls->tgt_cloud->closest_src_dist[i_part][cls->n_closest*i+j] =
          closest_src_dist[n_tgt + cls->n_closest*i + j];
      }
    }
    n_tgt += cls->n_closest * cls->tgt_cloud->n_points[i_part];
  }
  free (closest_src_gnum);
  free (closest_src_dist);


  /* Sort closest source points in ascending order of global id */
  int    *order = malloc(sizeof(int   ) * cls->n_closest);
  double *tmp   = malloc(sizeof(double) * cls->n_closest);
  for (int i_part = 0; i_part < cls->tgt_cloud->n_part; i_part++) {
    for (int i = 0; i < cls->tgt_cloud->n_points[i_part]; i++) {
      for (int j = 0; j < cls->n_closest; j++) {
        order[j] = j;
      }

      PDM_sort_long(cls->tgt_cloud->closest_src_gnum [i_part] + cls->n_closest*i,
                    order,
                    cls->n_closest);

      memcpy(tmp,
             cls->tgt_cloud->closest_src_dist[i_part] + cls->n_closest*i,
             sizeof(double) * cls->n_closest);
      for (int j = 0; j < cls->n_closest; j++) {
        cls->tgt_cloud->closest_src_dist[i_part][cls->n_closest*i+j] = tmp[order[j]];
      }
    }
  }
  free(tmp);
  free(order);


  //-->GPU
  /* Free parallel octree */
  PDM_para_octree_free (octree);
  //<--

  _closest_points_reverse_results(cls);


  /* Create ptp object */
  // TO DO: transport triplets to avoid costly gnum_location
  cls->ptp = PDM_part_to_part_create((const PDM_g_num_t **) cls->src_cloud->gnum,
                                     (const int          *) cls->src_cloud->n_points,
                                                            cls->src_cloud->n_part,
                                     (const PDM_g_num_t **) cls->tgt_cloud->gnum,
                                     (const int          *) cls->tgt_cloud->n_points,
                                                            cls->tgt_cloud->n_part,
                                     (const int         **) cls->src_cloud->tgt_in_src_idx,
                                     (const PDM_g_num_t **) cls->src_cloud->tgt_in_src,
                                                            cls->comm);


  PDM_timer_hang_on(cls->timer);

  cls->times_elapsed[END] = PDM_timer_elapsed(cls->timer);
  cls->times_cpu[END]     = PDM_timer_cpu(cls->timer);
  cls->times_cpu_u[END]   = PDM_timer_cpu_user(cls->timer);
  cls->times_cpu_s[END]   = PDM_timer_cpu_sys(cls->timer);

  PDM_timer_resume(cls->timer);
}


/**
 *
 * \brief Get closest source points global ids and (squared) distance
 *
 * \param [in]   cls                   Pointer to \ref PDM_closest_points object
 * \param [in]   i_part_tgt            Index of partition of the cloud
 * \param [out]  closest_src_gnum      Global number of the closest element (size = n_closest * n_tgt_points)
 * \param [out]  closest_src_distance  Distance (size = n_closest * n_tgt_points)
 *
 */

void
PDM_closest_points_get
(
       PDM_closest_point_t  *cls,
 const int                   i_part_tgt,
       PDM_g_num_t         **closest_src_gnum,
       double              **closest_src_distance
)
{

  assert (cls->tgt_cloud->closest_src_gnum != NULL);
  assert (cls->tgt_cloud->closest_src_dist != NULL);

  *closest_src_gnum     = cls->tgt_cloud->closest_src_gnum[i_part_tgt];
  *closest_src_distance = cls->tgt_cloud->closest_src_dist[i_part_tgt];

  cls->results_is_getted = PDM_TRUE;
}


/**
 *
 * \brief Get mesh distance
 *
 * \param [in]   cls                Pointer to \ref PDM_closest_points object
 * \param [in]   i_part_src         Index of partition of the cloud
 * \param [out]  tgt_in_src_idx     For each src point the number of target localised  (size = n_src_points )
 * \param [out]  tgt_in_src         For each src point the globla number of target point located (size = tgt_in_src_idx[n_src_points] )
 *
 */

void
PDM_closest_points_tgt_in_src_get
(
       PDM_closest_point_t  *cls,
 const int                   i_part_src,
       int                 **tgt_in_src_idx,
       PDM_g_num_t         **tgt_in_src
)
{

  assert (cls->src_cloud->tgt_in_src_idx != NULL);
  assert (cls->src_cloud->tgt_in_src != NULL);

  *tgt_in_src_idx = cls->src_cloud->tgt_in_src_idx[i_part_src];
  *tgt_in_src     = cls->src_cloud->tgt_in_src    [i_part_src];

  // int size = cls->src_cloud->n_points[i_part_src];
  // PDM_log_trace_array_long(cls->src_cloud->tgt_in_src_idx[i_part_src], size, "get -> tgt_in_src_idx :: " );
  // PDM_log_trace_array_long(cls->src_cloud->tgt_in_src    [i_part_src], (*tgt_in_src_idx)[size], "get -> tgt_in_src :: " );

  cls->tgt_in_src_results_is_getted = PDM_TRUE;
}

/**
 *
 * \brief Get mesh distance
 *
 * \param [in]   cls                Pointer to \ref PDM_closest_points object
 * \param [in]   i_part_src         Index of partition of the cloud
 * \param [out]  tgt_in_src_idx     For each src point the number of target localised  (size = n_src_points )
 * \param [out]  tgt_in_src_dist    For each src point the distance to the target point located (size = tgt_in_src_idx[n_src_points] )
 *
 */

void
PDM_closest_points_tgt_in_src_dist_get
(
       PDM_closest_point_t  *cls,
 const int                   i_part_src,
       int                 **tgt_in_src_idx,
       double              **tgt_in_src_dist
)
{

  assert (cls->src_cloud->tgt_in_src_idx != NULL);
  assert (cls->src_cloud->tgt_in_src_dist != NULL);

  *tgt_in_src_idx = cls->src_cloud->tgt_in_src_idx[i_part_src];
  *tgt_in_src_dist = cls->src_cloud->tgt_in_src_dist[i_part_src];

  cls->tgt_in_src_results_is_getted_d = PDM_TRUE;
}



/**
 *
 * \brief Free a closest points structure
 *
 * \param [in]  cls      Pointer to \ref PDM_closest_points object
 *
 */

void
PDM_closest_points_free
(
PDM_closest_point_t  *cls
)
{

  if(( cls->owner == PDM_OWNERSHIP_KEEP ) ||
     ( cls->owner == PDM_OWNERSHIP_UNGET_RESULT_IS_FREE && !cls->results_is_getted)){
    if (cls->tgt_cloud->closest_src_gnum != NULL) {
      for (int j = 0; j < cls->tgt_cloud->n_part ; j++) {
        if (cls->tgt_cloud->closest_src_gnum[j] != NULL) {
          free (cls->tgt_cloud->closest_src_gnum[j]);
        }
      }
    }

    if (cls->tgt_cloud->closest_src_dist != NULL) {
      for (int j = 0; j < cls->tgt_cloud->n_part ; j++) {
        if (cls->tgt_cloud->closest_src_dist[j] != NULL) {
          free (cls->tgt_cloud->closest_src_dist[j]);
        }
      }
    }
  }

  free (cls->tgt_cloud->closest_src_gnum);
  free (cls->tgt_cloud->closest_src_dist);

  int free_tgt_in_src_gnum = (cls->owner == PDM_OWNERSHIP_KEEP) ||
     ( cls->owner == PDM_OWNERSHIP_UNGET_RESULT_IS_FREE && !cls->tgt_in_src_results_is_getted)||
     ( cls->owner == PDM_OWNERSHIP_USER                 && !cls->tgt_in_src_results_is_getted); // Dernière condition pour le python essentiellement ou si un utilisateur n'a pas besoin de ce résultats
  int free_tgt_in_src_dist = (cls->owner == PDM_OWNERSHIP_KEEP) ||
     ( cls->owner == PDM_OWNERSHIP_UNGET_RESULT_IS_FREE && !cls->tgt_in_src_results_is_getted_d)||
     ( cls->owner == PDM_OWNERSHIP_USER                 && !cls->tgt_in_src_results_is_getted_d);

  if (free_tgt_in_src_gnum) {
    if (cls->src_cloud->tgt_in_src != NULL) {
      for (int j = 0; j < cls->src_cloud->n_part ; j++) {
        if (cls->src_cloud->tgt_in_src[j] != NULL) {
          free (cls->src_cloud->tgt_in_src[j]);
        }
      }
    }
  }
  if (free_tgt_in_src_dist) {
    if (cls->src_cloud->tgt_in_src_dist != NULL) {
      for (int j = 0; j < cls->src_cloud->n_part ; j++) {
        if (cls->src_cloud->tgt_in_src_dist[j] != NULL) {
          free (cls->src_cloud->tgt_in_src_dist[j]);
        }
      }
    }
  }
  if (free_tgt_in_src_gnum && free_tgt_in_src_dist) {
    if (cls->src_cloud->tgt_in_src_idx != NULL) {
      for (int j = 0; j < cls->src_cloud->n_part ; j++) {
        if (cls->src_cloud->tgt_in_src_idx[j] != NULL) {
          free (cls->src_cloud->tgt_in_src_idx[j]);
        }
      }
    }
  }

  free (cls->src_cloud->tgt_in_src_idx);
  free (cls->src_cloud->tgt_in_src);
  free (cls->src_cloud->tgt_in_src_dist);

  if (cls->tgt_cloud->gnum != NULL) {
    free (cls->tgt_cloud->gnum);
  }
  if (cls->tgt_cloud->coords != NULL) {
    free (cls->tgt_cloud->coords);
  }
  if (cls->tgt_cloud->n_points != NULL) {
    free (cls->tgt_cloud->n_points);
  }
  if (cls->tgt_cloud != NULL) {
    free (cls->tgt_cloud);
  }


  if (cls->src_cloud->gnum != NULL) {
    free (cls->src_cloud->gnum);
  }
  if (cls->src_cloud->coords != NULL) {
    free (cls->src_cloud->coords);
  }
  if (cls->src_cloud->n_points != NULL) {
    free (cls->src_cloud->n_points);
  }
  if (cls->src_cloud != NULL) {
    free (cls->src_cloud);
  }

  PDM_timer_free(cls->timer);

  if (cls->ptp_ownership == PDM_OWNERSHIP_KEEP) {
    PDM_part_to_part_free(cls->ptp);
    cls->ptp = NULL;
  }

  free (cls);

}


/**
 *
 * \brief  Dump elapsed and CPU time
 *
 * \param [in]  cls      Pointer to \ref PDM_closest_points object
 *
 */

void
PDM_closest_points_dump_times
(
PDM_closest_point_t  *cls
)
{
  double t1 = cls->times_elapsed[END] - cls->times_elapsed[BEGIN];
  double t2 = cls->times_cpu[END] - cls->times_cpu[BEGIN];

  double t1max;
  PDM_MPI_Allreduce (&t1, &t1max, 1, PDM_MPI_DOUBLE, PDM_MPI_MAX, cls->comm);

  double t2max;
  PDM_MPI_Allreduce (&t2, &t2max, 1, PDM_MPI_DOUBLE, PDM_MPI_MAX, cls->comm);

  int rank;
  PDM_MPI_Comm_rank (cls->comm, &rank);

  if (rank == 0) {

    PDM_printf( "closest_points timer : all (elapsed and cpu) : %12.5es %12.5es\n",
                t1max, t2max);
  }
}

/*
 * Reverse operation of child creation : from a parent global numbering
 * (parent_ln_to_gn) and a child global numbering (child_ln_to_gn) created
 * from it, take some child global numbers (gnum_to_transform) and retrieve
 * their original parent number.
 *
 * This function allow parent/child numbering to have a different partitionning than
 * gnum_to_transform. For both arrays, number of part and number of elt per part must
 * be provided.
 *
 * gnum_to_transform is modified inplace
*/

void
PDM_transform_to_parent_gnum
(
 const int           n_part_initial,
 const int          *n_elmt_initial,
 const PDM_g_num_t **child_ln_to_gn,
 const PDM_g_num_t **parent_ln_to_gn,
 const int           n_part_to_transform,
 const int          *n_elmt_to_transform,
       PDM_g_num_t **gnum_to_transform,
       PDM_MPI_Comm  comm
)
{
  PDM_part_to_block_t* ptb = PDM_part_to_block_create(PDM_PART_TO_BLOCK_DISTRIB_ALL_PROC,
                                                      PDM_PART_TO_BLOCK_POST_CLEANUP,
                                                      1.,
                                     (PDM_g_num_t **) child_ln_to_gn,
                                                      NULL,
                                              (int *) n_elmt_initial,
                                                      n_part_initial,
                                                      comm);

  int         *block_stride = NULL;
  PDM_g_num_t *block_parent = NULL;
  int s_block_data = PDM_part_to_block_exch (ptb,
                          sizeof(PDM_g_num_t),
                          PDM_STRIDE_CST_INTERLACED,
                          1,
                          NULL,
               (void **)  parent_ln_to_gn,
                         &block_stride,
               (void **) &block_parent);

  PDM_g_num_t *block_distrib_idx = PDM_part_to_block_distrib_index_get (ptb);

  if(0 == 1){
    PDM_log_trace_array_long(block_parent, s_block_data, "block_parent :: " );
  }

  PDM_block_to_part_t *btp = PDM_block_to_part_create(block_distrib_idx,
                               (const PDM_g_num_t **) gnum_to_transform,
                                                      n_elmt_to_transform,
                                                      n_part_to_transform,
                                                      comm);

  int stride_one = 1;
  PDM_block_to_part_exch_in_place(btp,
                         sizeof(PDM_g_num_t),
                         PDM_STRIDE_CST_INTERLACED,
                        &stride_one,
                         block_parent,
                         NULL,
               (void **) gnum_to_transform);


  PDM_part_to_block_free(ptb);
  PDM_block_to_part_free(btp);
  free(block_parent);
}

/**
 *
 * \brief  transfert _closest_pts var as it seems this static var is not readable
 *          when we switch to the nvcc compiler
 *
 */

PDM_closest_point_t *
PDM_closest_points_closest_transfert
(
  PDM_closest_point_t  *cls
)
{
  return cls;
}



/**
 *
 * \brief  Get the number of target points in a partition
 *
 * \param [in]  cls     Pointer to \ref PDM_closest_points object
 * \param [in]  i_part  Index of partition of the target cloud
 *
 * \return   Number of target point in the partition \ref i_part
 *
 */

int
PDM_closest_points_n_tgt_get
(
  PDM_closest_point_t  *cls,
  const int             i_part
)
{
  assert(cls->tgt_cloud != NULL);
  return cls->tgt_cloud->n_points[i_part];
}


/**
 *
 * \brief  Get the number of source points in a partition
 *
 * \param [in]  cls     Pointer to \ref PDM_closest_points object
 * \param [in]  i_part  Index of partition of the target cloud
 *
 * \return   Number of source point in the partition \ref i_part
 *
 */

int
PDM_closest_points_n_src_get
(
  PDM_closest_point_t  *cls,
  const int             i_part
)
{
  assert(cls->src_cloud != NULL);
  return cls->src_cloud->n_points[i_part];
}


/**
 *
 * \brief  Get the number of closest points
 *
 * \param [in]  cls     Pointer to \ref PDM_closest_points object
 *
 * \return   Number of closest points
 *
 */

int
PDM_closest_points_n_closest_get
(
  PDM_closest_point_t  *cls
)
{
  return cls->n_closest;
}


/**
 * \brief Get part_to_part object to exchange data between
 * the source and target point clouds (both in user frame)
 *
 * \param [in ] cls        Pointer to \ref PDM_closest_point_t object
 * \param [out] ptp        Pointer to \ref PDM_part_to_part_t object
 * \param [in ] ownership  Ownership for ptp
 *
 */

void
PDM_closest_points_part_to_part_get
(
 PDM_closest_point_t  *cls,
 PDM_part_to_part_t  **ptp,
 PDM_ownership_t       ownership
 )
{
  *ptp = cls->ptp;
  cls->ptp_ownership = ownership;
}

#ifdef	__cplusplus
}
#endif
#undef NTIMER
