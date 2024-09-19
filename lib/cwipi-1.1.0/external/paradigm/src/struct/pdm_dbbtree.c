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

#include "pdm_mpi.h"
#include "pdm.h"
#include "pdm_morton.h"
#include "pdm_mpi.h"
#include "pdm_box.h"
#include "pdm_sort.h"
#include "pdm_array.h"
#include "pdm_hash_tab.h"
#include "pdm_box_priv.h"
#include "pdm_box_tree.h"
#include "pdm_dbbtree.h"
#include "pdm_dbbtree_priv.h"
#include "pdm_printf.h"
#include "pdm_error.h"
#include "pdm_priv.h"
#include "pdm_distrib.h"
#include "pdm_vtk.h"
#include "pdm_part_to_block.h"
#include "pdm_block_to_part.h"
#include "pdm_box_tree_priv.h"
#include "pdm_timer.h"
#include "pdm_logging.h"
#include "pdm_binary_search.h"
#include "pdm_unique.h"

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force back Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */

/*=============================================================================
 * Macro definitions
 *============================================================================*/

#define _MIN(a,b)   ((a) < (b) ?  (a) : (b))  /* Minimum of a et b */

#define _MAX(a,b)   ((a) > (b) ?  (a) : (b))  /* Maximum of a et b */

/*============================================================================
 * Type
 *============================================================================*/

/*=============================================================================
 * Static global variables
 *============================================================================*/

/*=============================================================================
 * Static function definitions
 *============================================================================*/

/**
 * \brief  Normalize coorndinates
 */

static void
_normalize
(
 _PDM_dbbtree_t *dbbt,
 const double *pt_origin,
 double *pt_nomalized
 )
{
  for (int j = 0; j < dbbt->dim; j++) {
    pt_nomalized[j] = (pt_origin[j] - dbbt->s[j]) / dbbt->d[j];
  }
}

/**
 * \brief  Normalize normal vector coordinates
 */

static void
_normalize_normal_vector
(
 _PDM_dbbtree_t *dbbt,
 const double   *normal_vector,
 double         *normal_vector_nomalized
)
{
  for (int j = 0; j < dbbt->dim; j++) {
    normal_vector_nomalized[j] = normal_vector[j] * dbbt->d[j];
  }
}

/**
 * \brief  Initialize box_tree statistics
 *
 * \param [inout]  bts  pointer to box tree statistics structure
 *
 */

static void
_init_bt_statistics
(
 _box_tree_stats_t  *bts
 )
{
  size_t i;

  assert(bts != NULL);

  bts->dim = 0;

  for (i = 0; i < 3; i++) {
    bts->depth[i] = 0;
    bts->n_leaves[i] = 0;
    bts->n_boxes[i] = 0;
    bts->n_threshold_leaves[i] = 0;
    bts->n_leaf_boxes[i] = 0;
    bts->mem_used[i] = 0;
    bts->mem_required[i] = 0;
  }
}


/**
 * \brief Update box-tree statistics.
 *
 * For most fields, we replace previous values with the current ones.
 *
 * For memory required, we are interested in the maximum values over time
 * (i.e. algorthm steps); this is the case even for the minimal memory
 * required, we is thus the time maximum of the rank minimum.
 *
 * \param [inout]   bts   Pointer to box tree statistics structure
 * \param [inout]   bt    Pointer to box tree structure
 *
 */

static void
_update_bt_statistics
(
 _box_tree_stats_t     *bts,
 const PDM_box_tree_t  *bt
 )
{
  int dim;
  size_t i;
  size_t mem_required[3];

  assert(bts != NULL);

  dim = PDM_box_tree_get_stats (bt,
                                bts->depth,
                                bts->n_leaves,
                                bts->n_boxes,
                                bts->n_threshold_leaves,
                                bts->n_leaf_boxes,
                                bts->mem_used,
                                mem_required);

  bts->dim = dim;

  for (i = 0; i < 3; i++)
    bts->mem_required[i] = _MAX(bts->mem_required[i], mem_required[i]);
}


/**
 * \brief Distribute bounding boxes over the ranks according to a Morton encoding
 * index. Try to get a well-balanced distribution and spatially coherent.
 *
 * \param [inout]  n     <-> pointer to neighborhood management structure
 * \param [inout]  boxes <-> box set to redistribute
 *
 */

static void
_redistribute_boxes
(
 _PDM_dbbtree_t      *dbbt
 )
{

  /* Sanity checks */

  assert (dbbt != NULL);

  PDM_box_tree_t  *coarse_tree = PDM_box_tree_create (dbbt->maxTreeDepthCoarse,
                                                      dbbt->maxBoxesLeafCoarse,
                                                      dbbt->maxBoxRatioCoarse);

  /* Build a tree and associate boxes */

  PDM_box_tree_set_boxes (coarse_tree,
                          dbbt->boxes,
                          PDM_BOX_TREE_ASYNC_LEVEL);

  _update_bt_statistics(&(dbbt->btsCoarse), coarse_tree);

  if (1 == 0) {
    PDM_printf ("-- dump stats\n");

    PDM_box_tree_dump_statistics(coarse_tree);

    PDM_printf ("-- fin dump stats\n");

    PDM_printf ("-- dump \n");

    PDM_box_tree_dump(coarse_tree);

    PDM_printf ("-- fin dump\n");
  }

  if(0 == 1) {
    char filename[999];
    int i_rank;
    PDM_MPI_Comm_rank (dbbt->comm, &i_rank);
    sprintf(filename, "dbbt_coarse_tree_%3.3d.vtk",i_rank);
    PDM_vtk_write_boxes (filename,
                         dbbt->boxes->local_boxes->n_boxes,
                         dbbt->boxes->local_boxes->extents,
                         dbbt->boxes->local_boxes->g_num);
  }
  /*
   * Compute an index based on Morton encoding to ensure a good distribution
   * of bounding boxes among the ranks.
   */

  PDM_box_distrib_t  *distrib = PDM_box_tree_get_distrib (coarse_tree, dbbt->boxes);

  PDM_box_tree_destroy (&coarse_tree);

  if (1 == 0) {
    PDM_box_distrib_dump_statistics (distrib, dbbt->comm);
  }

  /* Define a new distribution of boxes according to the Morton
     encoding index */

  if (1 == 0) {
    PDM_printf("affichage 1\n");
    PDM_box_set_dump( dbbt->boxes,1);
    PDM_printf("fin affichage 1\n");
  }

  PDM_box_set_redistribute (distrib, dbbt->boxes);

  if (1 == 0) {
    PDM_printf("affichage 2\n");
    PDM_box_set_dump( dbbt->boxes,1);
    PDM_printf("fin affichage 2\n");
  }

  /* Delete intermediate structures */

  PDM_box_distrib_destroy (&distrib);

}



static
void
_adapt_tree_weight_for_intersect_line
(
 _PDM_dbbtree_t      *dbbt,
 PDM_box_tree_t      *coarse_tree,
 int                  n_line,
 double              *line_coord
)
{
  int debug = 0;
  int n_rank;
  int i_rank;
  PDM_MPI_Comm_size (dbbt->comm, &n_rank);
  PDM_MPI_Comm_rank (dbbt->comm, &i_rank);
  const int sExtents = dbbt->dim * 2;

  /*
   * Compute an index based on Morton encoding to ensure a good distribution
   * of bounding boxes among the ranks.
   */
  int     normalized       = 0;
  // int     n_depth_max      = 0;
  // int     n_extract_boxes  = 0;
  // double *extract_extents  = NULL;
  // int     n_extract_child  = 0;
  // int    *extract_child_id = NULL;
  // PDM_box_tree_extract_extents(coarse_tree,
  //                              normalized,
  //                              n_depth_max,
  //                              &n_extract_boxes,
  //                              &extract_extents,
  //                              &n_extract_child,
  //                              &extract_child_id);
  // PDM_log_trace_array_int(extract_child_id, n_extract_child, "extract_child_id ::");
  /*
   *  Normalize coordinates
   */
  if(0 == 1) {
    double *_line_coord = malloc (sizeof(double) * n_line * 6);
    for (int i = 0; i < 2*n_line; i++) {
      _normalize (dbbt, line_coord  + 3*i, _line_coord + 3*i);
    }

    if(debug) {
      char filename[999];
      sprintf(filename, "ray_normalize_%i.vtk", i_rank);
      PDM_vtk_write_lines(filename,
                          n_line,
                          _line_coord,
                          NULL,
                          NULL);
    }
    free(_line_coord);
  }

  PDM_g_num_t _n_g_line = n_line;
  PDM_g_num_t  n_g_line = 0;
  PDM_MPI_Allreduce(&_n_g_line, &n_g_line, 1,
                     PDM__PDM_MPI_G_NUM, PDM_MPI_SUM, dbbt->comm);

  double bucket_size_min = 0.2 * n_g_line / (n_rank);
  if (debug) {
    log_trace("bucket_size_min : %12.5e \n", bucket_size_min);
  }

  int n_max_it = 6;

  int  n_child_to_extract   = 1;
  int *child_ids_to_extract = malloc(n_child_to_extract * sizeof(int));

  int  n_all_child_ids = 0;
  int *all_child_ids                 = malloc(n_all_child_ids * sizeof(int));
  int *all_g_count_by_sampling_boxes = malloc(n_all_child_ids * sizeof(int));

  /* Init algo by root extract */
  child_ids_to_extract[0] = 0;

  for(int it = 0; it < n_max_it; ++it) {

    int     n_extract_boxes       = 0;
    double *extract_extents       = NULL;
    int     n_extract_child       = 0;
    int    *extract_child_id      = NULL;
    int    *extract_child_is_leaf = NULL;

    PDM_box_tree_extract_extents_by_child_ids(coarse_tree,
                                              normalized,
                                              n_child_to_extract,
                                              child_ids_to_extract,
                                              &n_extract_boxes,
                                              &extract_extents,
                                              &n_extract_child,
                                              &extract_child_id,
                                              &extract_child_is_leaf);

    if (debug) {
      PDM_log_trace_array_int(extract_child_id     , n_extract_child, "extract_child_id :");
      PDM_log_trace_array_int(extract_child_is_leaf, n_extract_child, "extract_child_is_leaf :");
    }

    if(debug) {
      char filename[999];
      sprintf(filename, "dbbt_extract_coarse_tree_it=%3.3d_%3.3d.vtk",it,i_rank);
      PDM_vtk_write_boxes (filename,
                           n_extract_boxes,
                           extract_extents,
                           NULL);
    }

    int *n_g_extract_boxes = (int *) malloc (sizeof(int) * n_rank);
    PDM_MPI_Allgather (&n_extract_boxes , 1, PDM_MPI_INT,
                       n_g_extract_boxes, 1, PDM_MPI_INT, dbbt->comm);

    /*
     *  Compute idx
     */
    int *g_extract_boxes_idx = (int *) malloc (sizeof(int) * (n_rank+1));
    g_extract_boxes_idx[0] = 0;
    for(int i = 0; i < n_rank; ++i) {
      g_extract_boxes_idx[i+1] = g_extract_boxes_idx[i] + n_g_extract_boxes[i];
    }

    if(0 == 1) {
      PDM_log_trace_array_int(g_extract_boxes_idx, n_rank+1, "g_extract_boxes_idx ::");
    }

    int n_g_extract_boxes_all = g_extract_boxes_idx[n_rank];
    PDM_g_num_t *g_num_sampling = (PDM_g_num_t * ) malloc( n_g_extract_boxes_all * sizeof(PDM_g_num_t));
    for(int i = 0; i < n_rank; ++i) {
      for(int j = g_extract_boxes_idx[i]; j < g_extract_boxes_idx[i+1]; ++j) {
        g_num_sampling[j] = i;
      }
    }

    double *g_sampling_extent = (double *) malloc (sizeof(double) * sExtents * n_g_extract_boxes_all);
    for(int i = 0; i < n_rank; ++i) {
      g_extract_boxes_idx[i] *= sExtents;
      n_g_extract_boxes  [i] *= sExtents;
    }

    PDM_MPI_Allgatherv(extract_extents  , n_extract_boxes   * sExtents , PDM__PDM_MPI_REAL,
                       g_sampling_extent, n_g_extract_boxes,
                       g_extract_boxes_idx,
                       PDM__PDM_MPI_REAL, dbbt->comm);
    for(int i = 0; i < n_rank; ++i) {
      g_extract_boxes_idx[i] /= sExtents;
      n_g_extract_boxes  [i] /= sExtents;
    }

    // PDM_log_trace_array_double(g_sampling_extent, n_g_extract_boxes_all * sExtents, "g_sampling_extent ::");

    int *init_location_proc = PDM_array_zeros_int(3 * n_g_extract_boxes_all);
    PDM_box_set_t* rank_boxes = PDM_box_set_create(3,
                                                   0,  // No normalization to preserve initial extents
                                                   0,  // No projection to preserve initial extents
                                                   n_g_extract_boxes_all,
                                                   g_num_sampling,
                                                   g_sampling_extent,
                                                   1,
                                                   &n_g_extract_boxes_all,
                                                   init_location_proc,
                                                   dbbt->comm);
    memcpy (rank_boxes->d, dbbt->d, sizeof(double) * 3);
    memcpy (rank_boxes->s, dbbt->s, sizeof(double) * 3);
    free   (init_location_proc);

    /*
     *  Build a shared box_tree to evaluate distribution
     *
     */
    PDM_box_tree_t* shared_box_tree = PDM_box_tree_create (dbbt->maxTreeDepthShared,
                                                           dbbt->maxBoxesLeafShared,
                                                           dbbt->maxBoxRatioShared);


    /* Build a tree and associate boxes */
    PDM_box_tree_set_boxes (shared_box_tree,
                            rank_boxes,
                            PDM_BOX_TREE_ASYNC_LEVEL);

    if(debug && i_rank == 0) {
      const char* filename = "dbbt_sampling_shared_tree.vtk";
      PDM_vtk_write_boxes (filename,
                           rank_boxes->local_boxes->n_boxes,
                           rank_boxes->local_boxes->extents,
                           rank_boxes->local_boxes->g_num);
    }

    /*
     * Compute for each line the number of intersection with bt shared
     *
     */
    int *line_rank_idx = NULL;
    int *line_rank     = NULL;
    PDM_box_tree_intersect_lines_boxes (shared_box_tree,
                                        -1,
                                        n_line,
                                        line_coord,
                                        &line_rank_idx,
                                        &line_rank);

    /* Count points to send to each rank */
    int* send_count_by_sampling_boxes = PDM_array_zeros_int (n_g_extract_boxes_all);
    for (int i = 0; i < line_rank_idx[n_line]; i++) {
      send_count_by_sampling_boxes[line_rank[i]]++;
    }

    if(0 == 1) {
      PDM_log_trace_array_int(send_count_by_sampling_boxes, n_g_extract_boxes_all, "send_count_by_sampling_boxes ::");
    }

    /*
     * For each extract_child_id associate a extra_weight
     */
    int* g_count_by_sampling_boxes = malloc(n_g_extract_boxes[i_rank] * sizeof(int));
    PDM_MPI_Reduce_scatter(send_count_by_sampling_boxes, g_count_by_sampling_boxes, n_g_extract_boxes, PDM_MPI_INT, PDM_MPI_SUM, dbbt->comm);

    if(debug) {
      PDM_log_trace_array_int(g_count_by_sampling_boxes, n_g_extract_boxes[i_rank], "g_count_by_sampling_boxes ::");
    }

    /*
     * Pour l'algo adaptatif il faut qu'on garde le lien extract_child_id + proc
     * A chaque passage on ne cherche une profondeur de plus uniquement sur le child id qui a trop de point
     *  Du coup on rafine l'arbre uniquement ou on a besoin
     */

    n_child_to_extract = 0;
    child_ids_to_extract          = realloc(child_ids_to_extract         ,                    n_g_extract_boxes[i_rank]  * sizeof(int));
    all_child_ids                 = realloc(all_child_ids                , (n_all_child_ids + n_g_extract_boxes[i_rank]) * sizeof(int));
    all_g_count_by_sampling_boxes = realloc(all_g_count_by_sampling_boxes, (n_all_child_ids + n_g_extract_boxes[i_rank]) * sizeof(int));

    for(int i = 0; i < n_g_extract_boxes[i_rank]; ++i) {
      // if(g_count_by_sampling_boxes[i] > 0) {
      if(g_count_by_sampling_boxes[i] > bucket_size_min) {
        child_ids_to_extract     [n_child_to_extract] = extract_child_id[i];
        g_count_by_sampling_boxes[n_child_to_extract] = g_count_by_sampling_boxes[i];
        n_child_to_extract++;

        // Keep all for global extraction
        if(extract_child_is_leaf[i] == 1 || it == n_max_it-1) {
          all_child_ids                [n_all_child_ids] = extract_child_id[i];
          all_g_count_by_sampling_boxes[n_all_child_ids] = g_count_by_sampling_boxes[i];
          n_all_child_ids++;
        }
      }
    }

    /*
     * On garde toutes les feuilles
     */
    // PDM_log_trace_array_int(g_count_by_sampling_boxes, n_child_to_extract, "slect_g_count_by_sampling_boxes ::");

    PDM_box_tree_assign_weight(coarse_tree,
                               n_child_to_extract,
                               extract_child_id,
                               g_count_by_sampling_boxes);

    PDM_box_tree_destroy(&shared_box_tree);
    PDM_box_set_destroy (&rank_boxes);
    free(n_g_extract_boxes);
    free(send_count_by_sampling_boxes);
    free(g_count_by_sampling_boxes);
    free(g_extract_boxes_idx);
    free(g_sampling_extent);
    free(g_num_sampling);

    free(extract_extents);
    free(extract_child_id);
    free(extract_child_is_leaf);
    free(line_rank_idx);
    free(line_rank);

  }
  free(child_ids_to_extract);


  /*
   * Setup implicit child_ln_to_gn
   */
  PDM_g_num_t *child_ln_to_gn = malloc(n_all_child_ids * sizeof(PDM_g_num_t));
  double      *child_weight   = malloc(n_all_child_ids * sizeof(double     ));

  PDM_g_num_t _n_all_child_ids = n_all_child_ids;
  PDM_g_num_t beg_num_abs;
  PDM_MPI_Scan(&_n_all_child_ids, &beg_num_abs, 1, PDM__PDM_MPI_G_NUM, PDM_MPI_SUM, dbbt->comm);
  beg_num_abs -= _n_all_child_ids;

  for(int i = 0; i < n_all_child_ids; ++i) {
    child_ln_to_gn[i] = beg_num_abs + i + 1;
    child_weight  [i] = (double) all_g_count_by_sampling_boxes[i];
  }

  // log_trace("beg_num_abs = "PDM_FMT_G_NUM" \n", beg_num_abs);

  int sampling_factor = 2;
  int n_iter_max      = 5;
  double tol = 0.10;
  PDM_g_num_t* equi_child_distrib = NULL;
  PDM_distrib_weight(sampling_factor,
                     n_rank,
                     1,
                     &n_all_child_ids,
  (const PDM_g_num_t **)  &child_ln_to_gn,
  (const double      **)  &child_weight,
                     n_iter_max,
                     tol,
                     dbbt->comm,
                     &equi_child_distrib);
  free(child_weight);

  PDM_g_num_t* raw_child_distrib = PDM_compute_entity_distribution(dbbt->comm, n_all_child_ids);
  if (debug) {
    PDM_log_trace_array_long(raw_child_distrib, n_rank+1, "raw_child_distrib ::");
    PDM_log_trace_array_long(equi_child_distrib, n_rank+1, "equi_child_distrib ::");
  }

  /*
   *  Extraction des boîtes --> Il faudra faire un  is_visited_box[box_id] + visited_box = box_id
   */
  if (debug) {
    // PDM_log_trace_array_int(all_child_ids, n_all_child_ids, "all_child_ids ::");
    PDM_log_trace_array_int(all_g_count_by_sampling_boxes, n_all_child_ids, "all_g_count_by_sampling_boxes ::");
  }

  /*
   * Equilibrate :
   *   - Get all boxes id for each required boxes
   *   - Send gnum / origin / extents
   */
  int* send_n = malloc((n_rank) * sizeof(int));
  int* recv_n = malloc((n_rank) * sizeof(int));
  for(int i = 0; i < n_rank; ++i) {
    send_n[i] = 0;
  }

  PDM_boxes_t *_local_boxes = dbbt->boxes->local_boxes;

  // int n_boxes = _local_boxes->n_boxes;

  // int* is_visited = malloc(n_boxes * sizeof(int));
  // for(int i = 0; i < n_boxes; ++i) {
  //   is_visited[i] = 0;
  // }

  int  *n_child_box_ids = (int  *) malloc(n_all_child_ids * sizeof(int  ));
  int **child_box_ids   = (int **) malloc(n_all_child_ids * sizeof(int *));

  for(int i = 0; i < n_all_child_ids; ++i) {
    PDM_g_num_t child_gnum = child_ln_to_gn[i];
    int t_rank = PDM_binary_search_gap_long(child_gnum-1, equi_child_distrib, n_rank+1);

    /* Get number of box in child box_tree */
    int i_child = all_child_ids[i];

    n_child_box_ids[i] = PDM_box_tree_get_box_ids(coarse_tree,
                                                  i_child,
                                                  &child_box_ids[i]);

    /* Reset is_visited */
    // for(int i = 0; i < n_boxes; ++i) {
    //   is_visited[i] = 0;
    // }

    // send_n[t_rank] += PDM_box_tree_get_node_n_boxes(coarse_tree, i_child);
    send_n[t_rank] += n_child_box_ids[i];

    // PDM_log_trace_array_int(child_box_ids[i],n_child_box_ids[i], "n_child_box_ids ::");

  }

  /*
   * TO MANAGE --> On peut envoyer plusieurs fois la même boite au même rang !!!
   */

  /*
   * Exchange size
   */
  PDM_MPI_Alltoall(send_n, 1, PDM_MPI_INT,
                   recv_n, 1, PDM_MPI_INT, dbbt->comm);

  if(debug) {
    PDM_log_trace_array_int(send_n, n_rank, "send_n ::");
    PDM_log_trace_array_int(recv_n, n_rank, "recv_n ::");
  }

  int* send_idx = malloc((n_rank+1) * sizeof(int));
  int* recv_idx = malloc((n_rank+1) * sizeof(int));
  send_idx[0] = 0;
  recv_idx[0] = 0;
  for(int i = 0; i < n_rank; ++i) {
    send_idx[i+1] = send_idx[i] + send_n[i];
    recv_idx[i+1] = recv_idx[i] + recv_n[i];
    send_n  [i  ] = 0;
  }

  /*
   * TODO :  Asynchrone
   */
  const int stride = dbbt->dim * 2;
  const int stride_origin= 3;

  /* Allocate */
  PDM_g_num_t *send_g_num   = malloc(send_idx[n_rank] *                 sizeof(PDM_g_num_t));
  double      *send_extents = malloc(send_idx[n_rank] * stride *        sizeof(double     ));
  int         *send_origin  = malloc(send_idx[n_rank] * stride_origin * sizeof(int        ));

  /* To be allocated before modification of recv_n / recv_idx */
  PDM_g_num_t *recv_g_num   = malloc(recv_idx[n_rank] *                 sizeof(PDM_g_num_t));
  double      *recv_extents = malloc(recv_idx[n_rank] * stride *        sizeof(double     ));
  int         *recv_origin  = malloc(recv_idx[n_rank] * stride_origin * sizeof(int        ));

  int n_recv_boxes = recv_idx[n_rank];

  for(int i = 0; i < n_all_child_ids; ++i) {
    PDM_g_num_t child_gnum = child_ln_to_gn[i];
    int t_rank = PDM_binary_search_gap_long(child_gnum-1, equi_child_distrib, n_rank+1);

    /*
     * Copy gnum
     */
    for(int j = 0; j < n_child_box_ids[i]; ++j) {

      int shift  = send_idx[t_rank] + send_n[t_rank];
      int box_id = child_box_ids[i][j];

      send_g_num[shift] = _local_boxes->g_num[box_id];

      for (int k = 0; k < stride; k++) {
        send_extents[shift*stride + k] = _local_boxes->extents[box_id*stride + k];
      }

      for (int k = 0; k < stride_origin; k++) {
        send_origin[shift*stride_origin + k] = _local_boxes->origin[box_id*stride_origin + k];
      }

      send_n[t_rank]++;
    }

  }

  /* Exchange boxes between processes */
  PDM_MPI_Alltoallv(send_g_num, send_n, send_idx, PDM__PDM_MPI_G_NUM,
                    recv_g_num, recv_n, recv_idx, PDM__PDM_MPI_G_NUM,
                    dbbt->comm);
  free(send_g_num);

  // PDM_log_trace_array_long(recv_g_num, recv_idx[n_rank], "recv_g_num :");

  /*  */
  for (int i = 0; i < n_rank; i++) {
    send_n  [i] *= stride;
    send_idx[i] *= stride;
    recv_n  [i] *= stride;
    recv_idx[i] *= stride;
  }

  /* Exchange extents */
  PDM_MPI_Alltoallv(send_extents, send_n, send_idx, PDM_MPI_DOUBLE,
                    recv_extents, recv_n, recv_idx, PDM_MPI_DOUBLE,
                    dbbt->comm);
  free(send_extents);

  for (int i = 0; i < n_rank; i++) {
    send_n  [i] = send_n  [i]/stride * stride_origin;
    send_idx[i] = send_idx[i]/stride * stride_origin;
    recv_n  [i] = recv_n  [i]/stride * stride_origin;
    recv_idx[i] = recv_idx[i]/stride * stride_origin;
  }

  PDM_MPI_Alltoallv(send_origin, send_n, send_idx, PDM_MPI_INT,
                    recv_origin, recv_n, recv_idx, PDM_MPI_INT,
                    dbbt->comm);
  free(send_origin);

  if(0 == 1) {
    char filename[999];
    sprintf(filename, "dbbt_experimental_tree_%3.3d.vtk",i_rank);
    PDM_vtk_write_boxes (filename,
                         n_recv_boxes,
                         recv_extents,
                         recv_g_num);
  }


  /*
   * Re-setup boxes
   */
  if (debug) {
    log_trace("Avant :  _local_boxes->n_boxes = %i --> %i \n",  _local_boxes->n_boxes, n_recv_boxes);
  }
  _local_boxes->n_boxes = n_recv_boxes;
  free(_local_boxes->g_num);
  free(_local_boxes->extents);
  free(_local_boxes->origin);

  _local_boxes->g_num   = recv_g_num  ;
  _local_boxes->extents = recv_extents;
  _local_boxes->origin  = recv_origin ;

  // free(recv_g_num  );
  // free(recv_extents);
  // free(recv_origin );

  free(send_idx);
  free(recv_idx);

  free(send_n);
  free(recv_n);

  for(int i = 0; i < n_all_child_ids; ++i) {
    free(child_box_ids[i]);
  }
  free(child_box_ids);
  free(n_child_box_ids);

  free(raw_child_distrib);
  free(equi_child_distrib);
  free(all_child_ids);
  free(child_ln_to_gn);
  free(all_g_count_by_sampling_boxes);


}


static void
_redistribute_boxes_for_intersect_line
(
 _PDM_dbbtree_t      *dbbt,
 int                  n_line,
 double              *line_coord
 )
{

  /* Sanity checks */

  assert (dbbt != NULL);

  PDM_box_tree_t  *coarse_tree = PDM_box_tree_create (dbbt->maxTreeDepthCoarse,
                                                      dbbt->maxBoxesLeafCoarse,
                                                      dbbt->maxBoxRatioCoarse);

  /* Build a tree and associate boxes */

  PDM_box_tree_set_boxes (coarse_tree,
                          dbbt->boxes,
                          PDM_BOX_TREE_ASYNC_LEVEL);

  _update_bt_statistics(&(dbbt->btsCoarse), coarse_tree);

  if (1 == 0) {
    PDM_printf ("-- dump stats\n");

    PDM_box_tree_dump_statistics(coarse_tree);

    PDM_printf ("-- fin dump stats\n");

    PDM_printf ("-- dump \n");

    PDM_box_tree_dump(coarse_tree);

    PDM_printf ("-- fin dump\n");
  }

  if(0 == 1) {
    char filename[999];
    int i_rank;
    PDM_MPI_Comm_rank (dbbt->comm, &i_rank);
    sprintf(filename, "dbbt_coarse_tree_%3.3d.vtk",i_rank);
    PDM_vtk_write_boxes (filename,
                         dbbt->boxes->local_boxes->n_boxes,
                         dbbt->boxes->local_boxes->extents,
                         dbbt->boxes->local_boxes->g_num);
  }

  // _adapt_tree_weight_for_intersect_line(dbbt, coarse_tree, n_line, line_coord);

  PDM_box_distrib_t  *distrib = PDM_box_tree_get_distrib (coarse_tree, dbbt->boxes);

  PDM_box_tree_destroy (&coarse_tree);

  if (0 == 1) {
    PDM_box_distrib_dump_statistics (distrib, dbbt->comm);
  }

  /* Define a new distribution of boxes according to the Morton
     encoding index */

  if (1 == 0) {
    PDM_printf("affichage 1\n");
    PDM_box_set_dump( dbbt->boxes,1);
    PDM_printf("fin affichage 1\n");
  }

  PDM_box_set_redistribute (distrib, dbbt->boxes);

  /*
   * Rebuild coarse tree but with new distrib
   */
  coarse_tree = PDM_box_tree_create (dbbt->maxTreeDepthCoarse,
                                     dbbt->maxBoxesLeafCoarse,
                                     dbbt->maxBoxRatioCoarse);

  PDM_box_tree_set_boxes (coarse_tree,
                          dbbt->boxes,
                          PDM_BOX_TREE_ASYNC_LEVEL);

  _update_bt_statistics(&(dbbt->btsCoarse), coarse_tree);
  _adapt_tree_weight_for_intersect_line(dbbt, coarse_tree, n_line, line_coord);

  PDM_box_tree_destroy (&coarse_tree);

  /*
   * Rebuild coarse tree but with new distrib
   */
  coarse_tree = PDM_box_tree_create (dbbt->maxTreeDepthCoarse,
                                     dbbt->maxBoxesLeafCoarse,
                                     dbbt->maxBoxRatioCoarse);

  PDM_box_tree_set_boxes (coarse_tree,
                          dbbt->boxes,
                          PDM_BOX_TREE_ASYNC_LEVEL);

  _update_bt_statistics(&(dbbt->btsCoarse), coarse_tree);
  _adapt_tree_weight_for_intersect_line(dbbt, coarse_tree, n_line, line_coord);

  PDM_box_tree_destroy (&coarse_tree);

  if(0 == 1) {
    char filename[999];
    int i_rank;
    PDM_MPI_Comm_rank (dbbt->comm, &i_rank);
    sprintf(filename, "dbbt_redistribute_coarse_tree_%3.3d.vtk",i_rank);
    PDM_vtk_write_boxes (filename,
                         dbbt->boxes->local_boxes->n_boxes,
                         dbbt->boxes->local_boxes->extents,
                         dbbt->boxes->local_boxes->g_num);
  }

  if (0 == 1) {
    PDM_printf("affichage 2\n");
    PDM_box_set_dump( dbbt->boxes,1);
    PDM_printf("fin affichage 2\n");
  }

  /* Delete intermediate structures */

  PDM_box_distrib_destroy (&distrib);

}


static void _export_point_cloud
(
 char         *filename,
 int           n_part,
 const int    *n_pts,
 double      **coord,
 PDM_g_num_t **g_num,
 PDM_g_num_t **parent_g_num
 )
{
  FILE *f = fopen(filename, "w");

  fprintf(f, "# vtk DataFile Version 2.0\noctree points\nASCII\nDATASET UNSTRUCTURED_GRID\n");

  int n_pts_t = 0;
  for (int ipart = 0; ipart < n_part; ipart++) {
    n_pts_t += n_pts[ipart];
  }

  fprintf(f, "POINTS %d double\n", n_pts_t);
  for (int ipart = 0; ipart < n_part; ipart++) {
    for (int i = 0; i < n_pts[ipart]; i++) {
      for (int j = 0; j < 3; j++) {
        fprintf(f, "%f ", coord[ipart][3*i + j]);
      }
      fprintf(f, "\n");
    }
  }

  fprintf(f, "CELLS %d %d\n", n_pts_t, 2*n_pts_t);
  for (int i = 0; i < n_pts_t; i++) {
    fprintf(f, "1 %d\n", i);
  }

  fprintf(f, "CELL_TYPES %d\n", n_pts_t);
  for (int i = 0; i < n_pts_t; i++) {
    fprintf(f, "1\n");
  }

  fprintf(f, "CELL_DATA %d\n", n_pts_t);
  fprintf(f, "SCALARS gnum int\n LOOKUP_TABLE default\n");
  for (int ipart = 0; ipart < n_part; ipart++) {
    for (int i = 0; i < n_pts[ipart]; i++) {
      fprintf(f, ""PDM_FMT_G_NUM"\n", g_num[ipart][i]);
    }
  }

  if (parent_g_num != NULL) {
    fprintf(f, "FIELD FieldData 1\n");
    fprintf(f, "parent_gnum 1 %d int\n", n_pts_t);
    for (int ipart = 0; ipart < n_part; ipart++) {
      for (int i = 0; i < n_pts[ipart]; i++) {
        fprintf(f, ""PDM_FMT_G_NUM"\n", parent_g_num[ipart][i]);
      }
    }
  }

  fclose(f);
}


/*=============================================================================
 * Public function definitions
 *============================================================================*/

/**
 * \brief Return an intialized \ref PDM_dbbtree_t structure
 *
 * This function returns an initialized \ref PDM_dbbtree_t structure
 *
 * \param [in]  comm             Associated communicator
 * \param [in]  dim              boxes dimension
 * \param [in]  global_extents   Globals of elements to storage into the tree
 *                               (automatic computation if NULL)
 *
 * \return      A new initialized \ref PDM_dbbtree_t structure
 *
 */

PDM_dbbtree_t *
PDM_dbbtree_create
(
 PDM_MPI_Comm  comm,
 int           dim,
 double       *global_extents
 )
{
  _PDM_dbbtree_t *_dbbt = (_PDM_dbbtree_t *) malloc(sizeof(_PDM_dbbtree_t));

  _dbbt->comm                 = comm;
  _dbbt->rankComm             = PDM_MPI_COMM_NULL;
  _dbbt->dim                  = dim;

  _dbbt->maxTreeDepth         = 30;
  _dbbt->maxBoxRatio          = 10.;
  _dbbt->maxBoxesLeaf         = 30;

  _dbbt->maxTreeDepthShared   = 10;
  _dbbt->maxBoxRatioShared    =  6;
  _dbbt->maxBoxesLeafShared   =  5;

  _dbbt->maxTreeDepthCoarse   = 10;
  _dbbt->maxBoxRatioCoarse    =  4;
  _dbbt->maxBoxesLeafCoarse   = 30;

  _dbbt->rankBoxes            = NULL;
  _dbbt->btShared             = NULL;

  _dbbt->nUsedRank            = 0;
  _dbbt->usedRank             = NULL;

  _dbbt->boxes                = NULL;
  _dbbt->btLoc                = NULL;

  _dbbt->global_extents       = NULL;

  if (global_extents != NULL) {
    _dbbt->global_extents = malloc (sizeof(double) * dim * 2);
    memcpy(_dbbt->global_extents, global_extents, sizeof(double) * dim * 2);
    for (int j = 0; j < dim; j++) {
      _dbbt->s[j] = _dbbt->global_extents[j];
      _dbbt->d[j] = _dbbt->global_extents[j+dim] - _dbbt->global_extents[j];
    }
  }

  else {
    for (int j = 0; j < dim; j++) {
      _dbbt->s[j] = 0.;
      _dbbt->d[j] = 1.;
    }
  }

  _init_bt_statistics (&(_dbbt->btsShared));
  _init_bt_statistics (&(_dbbt->btsLoc));
  _init_bt_statistics (&(_dbbt->btsCoarse));

  return (PDM_dbbtree_t *) _dbbt;
}



/**
 * \brief Free a \ref PDM_dbbtree_t structure
 *
 * \param [in]  dbbt   Pointer to a distributed bounding box tree
 *
 * \return      NULL
 *
 */

PDM_dbbtree_t *
PDM_dbbtree_free
(
 PDM_dbbtree_t     *dbbt
 )
{
  if (dbbt != NULL) {
    _PDM_dbbtree_t *_dbbt = (_PDM_dbbtree_t *) dbbt;

    //PDM_box_set_destroy (&_dbbt->boxes);

    if (_dbbt->global_extents != NULL) {
      free (_dbbt->global_extents);
    }

    free (_dbbt->usedRank);
    if (_dbbt->rankComm != PDM_MPI_COMM_NULL) {
      PDM_MPI_Comm_free (&(_dbbt->rankComm));
    }

    PDM_box_tree_destroy (&_dbbt->btShared);
    PDM_box_tree_destroy (&_dbbt->btLoc);
    PDM_box_set_destroy (&_dbbt->rankBoxes);

    free (_dbbt);
  }

  return NULL;
}


PDM_box_set_t *
PDM_dbbtree_boxes_set_with_init_location
(
 PDM_dbbtree_t      *dbbt,
 const int           n_part,
 const int          *nElts,
 const int         **init_location,
 const double      **extents,
 const PDM_g_num_t **gNum
)
{
  assert (dbbt != NULL);
  _PDM_dbbtree_t *_dbbt = (_PDM_dbbtree_t *) dbbt;

  int myRank;
  PDM_MPI_Comm_rank (_dbbt->comm, &myRank);
  int lComm;
  PDM_MPI_Comm_size (_dbbt->comm, &lComm);

  const int nInfoLocation = 3;
  const int sExtents = _dbbt->dim * 2;

  int nEltsProc = 0;
  for (int i = 0; i < n_part; i++) {
    nEltsProc += nElts[i];
  }

  PDM_g_num_t *_boxGnum = (PDM_g_num_t *) malloc (sizeof(PDM_g_num_t) * nEltsProc);
  double *_extents      = (double *) malloc (sizeof(double) * nEltsProc * sExtents);
  int *_initLocation   = (int *) malloc (sizeof(int) * nEltsProc * nInfoLocation);

  int idx = 0;
  int idx1 = 0;
  int idx2 = 0;

  if(init_location == NULL) {
    for (int i = 0; i < n_part; i++) {
      for (int j = 0; j < nElts[i]; j++) {
        _boxGnum[idx++] = gNum[i][j];

        for (int k = 0; k < sExtents; k++) {
          _extents[idx1++] = extents[i][sExtents*j+k];
        }

        _initLocation[idx2++] = myRank;
        _initLocation[idx2++] = i;
        _initLocation[idx2++] = j;
      }
    }
  } else {
    for (int i = 0; i < n_part; i++) {
      for (int j = 0; j < nElts[i]; j++) {
        _boxGnum[idx++] = gNum[i][j];

        for (int k = 0; k < sExtents; k++) {
          _extents[idx1++] = extents[i][sExtents*j+k];
        }

        _initLocation[idx2++] = init_location[i][3*j  ];
        _initLocation[idx2++] = init_location[i][3*j+1];
        _initLocation[idx2++] = init_location[i][3*j+2];
      }
    }
  }

  /*
   * Redistribute boxes of mesh A
   */

  if (0 == 1) {

    PDM_printf ("nEltsProc : %d\n", nEltsProc);

    PDM_printf ("_boxGnum :");
    for (int i = 0; i < nEltsProc; i++) {
      PDM_printf (" "PDM_FMT_G_NUM, _boxGnum[i]);
    }
    PDM_printf ("\n");

    PDM_printf ("_extents m2 :\n");
    idx1 = 0;
    for (int i = 0; i < nEltsProc; i++) {
      PDM_printf (" "PDM_FMT_G_NUM":", _boxGnum[i]);
      PDM_printf (" %12.5e", _extents[idx1++]);
      PDM_printf (" %12.5e", _extents[idx1++]);
      PDM_printf (" %12.5e", _extents[idx1++]);
      PDM_printf (" %12.5e", _extents[idx1++]);
      PDM_printf (" %12.5e", _extents[idx1++]);
      PDM_printf (" %12.5e", _extents[idx1++]);
      PDM_printf ("\n");
    }
    PDM_printf ("\n");

    PDM_printf ("_initLocation :");
    idx1 = 0;
    for (int i = 0; i < nEltsProc; i++) {
      PDM_printf (" %d", _initLocation[idx1++]);
      PDM_printf (" %d", _initLocation[idx1++]);
      PDM_printf (" %d", _initLocation[idx1++]);
      PDM_printf ("\n");
    }
    PDM_printf ("\n");
  }

  int ind_norm = 1;
  if (_dbbt->global_extents != NULL) {
    ind_norm = 0;
    for (int i = 0; i < nEltsProc; i++) {
      _normalize (_dbbt,
                  _extents+2*i*_dbbt->dim,
                  _extents+2*i*_dbbt->dim);
      _normalize (_dbbt,
                  _extents+(2*i+1)*_dbbt->dim,
                  _extents+(2*i+1)*_dbbt->dim);
    }
  }

  _dbbt->boxes = PDM_box_set_create(3,
                                    ind_norm,  // No normalization to preserve initial extents
                                    0,  // No projection to preserve initial extents
                                    nEltsProc,
                                    _boxGnum,
                                    _extents,
                                    n_part,
                                    nElts,
                                    _initLocation,
                                    _dbbt->comm);

  if (_dbbt->global_extents == NULL) {
    memcpy (_dbbt->d, _dbbt->boxes->d, sizeof(double) * 3);
    memcpy (_dbbt->s, _dbbt->boxes->s, sizeof(double) * 3);
  }

  free (_boxGnum);
  free (_extents);
  free (_initLocation);

  if (lComm > 1) {
    _redistribute_boxes(_dbbt);

    /*
     * Compute processus extents
     */

    int nBoxes = PDM_box_set_get_size (_dbbt->boxes);
    const double *extents2 = PDM_box_set_get_extents (_dbbt->boxes);

    double gExtents[sExtents];
    for (int i = 0; i < _dbbt->dim; i++) {
      gExtents[i]   =  DBL_MAX;
      gExtents[_dbbt->dim+i] = -DBL_MAX;
    }


    for (int i = 0; i < nBoxes; i++) {
      for (int k1 = 0; k1 < _dbbt->dim; k1++) {
        gExtents[k1]   = _MIN (gExtents[k1], extents2[sExtents * i + k1]);
        gExtents[_dbbt->dim+k1] = _MAX (gExtents[_dbbt->dim+k1], extents2[sExtents * i
                                                                          + _dbbt->dim + k1]);
      }
    }

    /*
     * Exchange extents to build shared boxes
     */

    int *allNBoxes = (int *) malloc (sizeof(int) * lComm);
    PDM_MPI_Allgather (&nBoxes, 1, PDM_MPI_INT,
                       allNBoxes, 1, PDM_MPI_INT,
                       _dbbt->comm);

    int nUsedRank = 0;
    for (int i = 0; i < lComm; i++) {
      if (allNBoxes[i] > 0) {
        nUsedRank += 1;
      }
    }


    double *allGExtents = (double *) malloc (sizeof(double) * sExtents * lComm);
    PDM_MPI_Allgather (gExtents, sExtents, PDM__PDM_MPI_REAL,
                       allGExtents, sExtents, PDM__PDM_MPI_REAL,
                       _dbbt->comm);

    int *numProc = (int *) malloc (sizeof(int *) * nUsedRank);

    _dbbt->usedRank = numProc;
    _dbbt->nUsedRank = nUsedRank;

    PDM_g_num_t *gNumProc = (PDM_g_num_t *) malloc (sizeof(PDM_g_num_t) * nUsedRank);

    idx = 0;
    for (int i = 0; i < lComm; i++) {
      if (allNBoxes[i] > 0) {
        gNumProc[idx] = idx;
        numProc[idx] = i;
        for (int j = 0; j < sExtents; j++) {
          allGExtents[idx*sExtents + j] = allGExtents[i*sExtents + j];
        }
        idx += 1;
      }
    }

    free (allNBoxes);

    allGExtents = (double *) realloc (allGExtents, sizeof(double) * sExtents * nUsedRank);

    int *initLocation_proc = PDM_array_zeros_int(3 * nUsedRank);

    //TODO: Faire un PDM_box_set et PDM_box_tree_create sequentiel ! Le comm split a u n 1 proc ici : pas terrible

    PDM_MPI_Comm_split(_dbbt->comm, myRank, 0, &(_dbbt->rankComm));

    _dbbt->rankBoxes = PDM_box_set_create(3,
                                          0,  // No normalization to preserve initial extents
                                          0,  // No projection to preserve initial extents
                                          nUsedRank,
                                          gNumProc,
                                          allGExtents,
                                          1,
                                          &nUsedRank,
                                          initLocation_proc,
                                          _dbbt->rankComm);

    memcpy (_dbbt->rankBoxes->d, _dbbt->d, sizeof(double) * 3);
    memcpy (_dbbt->rankBoxes->s, _dbbt->s, sizeof(double) * 3);

    _dbbt->btShared = PDM_box_tree_create (_dbbt->maxTreeDepthShared,
                                           _dbbt->maxBoxesLeafShared,
                                           _dbbt->maxBoxRatioShared);

    /* Build a tree and associate boxes */

    PDM_box_tree_set_boxes (_dbbt->btShared,
                            _dbbt->rankBoxes,
                            PDM_BOX_TREE_ASYNC_LEVEL);

    _update_bt_statistics(&(_dbbt->btsShared), _dbbt->btShared);

    if(0 == 1 && myRank == 0) {
      const char* filename = "dbbt_shared_tree.vtk";
      PDM_vtk_write_boxes (filename,
                           _dbbt->rankBoxes->local_boxes->n_boxes,
                           _dbbt->rankBoxes->local_boxes->extents,
                           _dbbt->rankBoxes->local_boxes->g_num);
    }

    free (allGExtents);
    free (gNumProc);
    free (initLocation_proc);

  }

  /*
   * Build local bt
   */

  _dbbt->btLoc = PDM_box_tree_create (_dbbt->maxTreeDepth,
                                      _dbbt->maxBoxesLeaf,
                                      _dbbt->maxBoxRatio);

  /* Build a tree and associate boxes */

  PDM_box_tree_set_boxes (_dbbt->btLoc,
                          _dbbt->boxes,
                          PDM_BOX_TREE_ASYNC_LEVEL);

  _update_bt_statistics(&(_dbbt->btsLoc), _dbbt->btLoc);

  char *env_var_oct = getenv ("OCTREE_SHARED");
  int use_shared_octree = 0;
  if (env_var_oct != NULL) {
    use_shared_octree = atoi(env_var_oct);
  }

  if(use_shared_octree == 1) {
    PDM_box_tree_copy_to_shm(_dbbt->btLoc);
  }

  // double dt = PDM_MPI_Wtime() - t1;
  // log_trace("PDM_dbbtree_boxes_set : %12.5e \n", dt);

  return _dbbt->boxes;
}

/**
 * \brief Assign a set of boxes to an empty \ref PDM_dbbtree_t structure.
 *
 * This function assigns a set of boxes to an empty \ref PDM_dbbtree_t structure.
 *
 * \param [in]  dbbt     Pointer to a distributed bounding box tree
 * \param [in]  n_part    Number of partitions
 * \param [in]  nElts    Number of elements of each partition
 * \param [in]  extents  Extents of each element of each partition
 * \param [in]  gNum     Global number of each element of each partition
 *
 * \return associated \ref PDM_box_set_t structure distributed according to
 * the tree location
 *
 */

PDM_box_set_t *
PDM_dbbtree_boxes_set
(
 PDM_dbbtree_t      *dbbt,
 const int           n_part,
 const int          *nElts,
 const double      **extents,
 const PDM_g_num_t **gNum
)
{
  return PDM_dbbtree_boxes_set_with_init_location(dbbt, n_part, nElts, NULL, extents, gNum);
}


/**
 * \brief Assign a set of boxes to an empty \ref PDM_dbbtree_t structure.
 *
 * This function assigns a set of boxes to an empty \ref PDM_dbbtree_t structure.
 *
 * \param [in]  dbbt     Pointer to a distributed bounding box tree
 * \param [in]  n_part    Number of partitions
 * \param [in]  nElts    Number of elements of each partition
 * \param [in]  extents  Extents of each element of each partition
 * \param [in]  gNum     Global number of each element of each partition
 *
 * \return associated \ref PDM_box_set_t structure distributed according to
 * the tree location
 *
 */

PDM_box_set_t *
PDM_dbbtree_boxes_set_for_intersect_line
(
 PDM_dbbtree_t      *dbbt,
 const int           n_part,
 const int          *nElts,
 const double      **extents,
 const PDM_g_num_t **gNum,
 const int           n_line,
 double             *line_coord
)
{
  assert (dbbt != NULL);
  _PDM_dbbtree_t *_dbbt = (_PDM_dbbtree_t *) dbbt;

  int myRank;
  PDM_MPI_Comm_rank (_dbbt->comm, &myRank);
  int lComm;
  PDM_MPI_Comm_size (_dbbt->comm, &lComm);

  const int nInfoLocation = 3;
  const int sExtents = _dbbt->dim * 2;

  int nEltsProc = 0;
  for (int i = 0; i < n_part; i++) {
    nEltsProc += nElts[i];
  }

  PDM_g_num_t *_boxGnum = (PDM_g_num_t *) malloc (sizeof(PDM_g_num_t) * nEltsProc);
  double *_extents      = (double *) malloc (sizeof(double) * nEltsProc * sExtents);
  int *_initLocation   = (int *) malloc (sizeof(int) * nEltsProc * nInfoLocation);

  int idx = 0;
  int idx1 = 0;
  int idx2 = 0;

  for (int i = 0; i < n_part; i++) {

    for (int j = 0; j < nElts[i]; j++) {
      _boxGnum[idx++] = gNum[i][j];

      for (int k = 0; k < sExtents; k++) {
        _extents[idx1++] = extents[i][sExtents*j+k];
      }

      _initLocation[idx2++] = myRank;
      _initLocation[idx2++] = i;
      _initLocation[idx2++] = j;
    }
  }

  /*
   * Redistribute boxes of mesh A
   */

  if (0 == 1) {

    PDM_printf ("nEltsProc : %d\n", nEltsProc);

    PDM_printf ("_boxGnum :");
    for (int i = 0; i < nEltsProc; i++) {
      PDM_printf (" "PDM_FMT_G_NUM, _boxGnum[i]);
    }
    PDM_printf ("\n");

    PDM_printf ("_extents m2 :\n");
    idx1 = 0;
    for (int i = 0; i < nEltsProc; i++) {
      PDM_printf (" "PDM_FMT_G_NUM":", _boxGnum[i]);
      PDM_printf (" %12.5e", _extents[idx1++]);
      PDM_printf (" %12.5e", _extents[idx1++]);
      PDM_printf (" %12.5e", _extents[idx1++]);
      PDM_printf (" %12.5e", _extents[idx1++]);
      PDM_printf (" %12.5e", _extents[idx1++]);
      PDM_printf (" %12.5e", _extents[idx1++]);
      PDM_printf ("\n");
    }
    PDM_printf ("\n");

    PDM_printf ("_initLocation :");
    idx1 = 0;
    for (int i = 0; i < nEltsProc; i++) {
      PDM_printf (" %d", _initLocation[idx1++]);
      PDM_printf (" %d", _initLocation[idx1++]);
      PDM_printf (" %d", _initLocation[idx1++]);
      PDM_printf ("\n");
    }
    PDM_printf ("\n");
  }

  int ind_norm = 1;
  if (_dbbt->global_extents != NULL) {
    ind_norm = 0;
    for (int i = 0; i < nEltsProc; i++) {
      _normalize (_dbbt,
                  _extents+2*i*_dbbt->dim,
                  _extents+2*i*_dbbt->dim);
      _normalize (_dbbt,
                  _extents+(2*i+1)*_dbbt->dim,
                  _extents+(2*i+1)*_dbbt->dim);
    }
  }

  _dbbt->boxes = PDM_box_set_create(3,
                                    ind_norm,  // No normalization to preserve initial extents
                                    0,  // No projection to preserve initial extents
                                    nEltsProc,
                                    _boxGnum,
                                    _extents,
                                    n_part,
                                    nElts,
                                    _initLocation,
                                    _dbbt->comm);

  if (_dbbt->global_extents == NULL) {
    memcpy (_dbbt->d, _dbbt->boxes->d, sizeof(double) * 3);
    memcpy (_dbbt->s, _dbbt->boxes->s, sizeof(double) * 3);
  }

  free (_boxGnum);
  free (_extents);
  free (_initLocation);

  if (lComm > 1) {
    _redistribute_boxes_for_intersect_line(_dbbt, n_line, line_coord);

    /*
     * Compute processus extents
     */

    int nBoxes = PDM_box_set_get_size (_dbbt->boxes);
    const double *extents2 = PDM_box_set_get_extents (_dbbt->boxes);

    double gExtents[sExtents];
    for (int i = 0; i < _dbbt->dim; i++) {
      gExtents[i]   =  DBL_MAX;
      gExtents[_dbbt->dim+i] = -DBL_MAX;
    }


    for (int i = 0; i < nBoxes; i++) {
      for (int k1 = 0; k1 < _dbbt->dim; k1++) {
        gExtents[k1]   = _MIN (gExtents[k1], extents2[sExtents * i + k1]);
        gExtents[_dbbt->dim+k1] = _MAX (gExtents[_dbbt->dim+k1], extents2[sExtents * i
                                                                          + _dbbt->dim + k1]);
      }
    }

    /*
     * Exchange extents to build shared boxes
     */

    int *allNBoxes = (int *) malloc (sizeof(int) * lComm);
    PDM_MPI_Allgather (&nBoxes, 1, PDM_MPI_INT,
                       allNBoxes, 1, PDM_MPI_INT,
                       _dbbt->comm);

    int nUsedRank = 0;
    for (int i = 0; i < lComm; i++) {
      if (allNBoxes[i] > 0) {
        nUsedRank += 1;
      }
    }


    double *allGExtents = (double *) malloc (sizeof(double) * sExtents * lComm);
    PDM_MPI_Allgather (gExtents, sExtents, PDM__PDM_MPI_REAL,
                       allGExtents, sExtents, PDM__PDM_MPI_REAL,
                       _dbbt->comm);

    int *numProc = (int *) malloc (sizeof(int *) * nUsedRank);

    _dbbt->usedRank = numProc;
    _dbbt->nUsedRank = nUsedRank;

    PDM_g_num_t *gNumProc = (PDM_g_num_t *) malloc (sizeof(PDM_g_num_t) * nUsedRank);

    idx = 0;
    for (int i = 0; i < lComm; i++) {
      if (allNBoxes[i] > 0) {
        gNumProc[idx] = idx;
        numProc[idx] = i;
        for (int j = 0; j < sExtents; j++) {
          allGExtents[idx*sExtents + j] = allGExtents[i*sExtents + j];
        }
        idx += 1;
      }
    }

    free (allNBoxes);

    allGExtents = (double *) realloc (allGExtents, sizeof(double) * sExtents * nUsedRank);

    int *initLocation_proc = PDM_array_zeros_int(3 * nUsedRank);

    //TODO: Faire un PDM_box_set et PDM_box_tree_create sequentiel ! Le comm split a u n 1 proc ici : pas terrible

    PDM_MPI_Comm_split(_dbbt->comm, myRank, 0, &(_dbbt->rankComm));

    _dbbt->rankBoxes = PDM_box_set_create(3,
                                          0,  // No normalization to preserve initial extents
                                          0,  // No projection to preserve initial extents
                                          nUsedRank,
                                          gNumProc,
                                          allGExtents,
                                          1,
                                          &nUsedRank,
                                          initLocation_proc,
                                          _dbbt->rankComm);

    memcpy (_dbbt->rankBoxes->d, _dbbt->d, sizeof(double) * 3);
    memcpy (_dbbt->rankBoxes->s, _dbbt->s, sizeof(double) * 3);

    _dbbt->btShared = PDM_box_tree_create (_dbbt->maxTreeDepthShared,
                                           _dbbt->maxBoxesLeafShared,
                                           _dbbt->maxBoxRatioShared);

    /* Build a tree and associate boxes */

    PDM_box_tree_set_boxes (_dbbt->btShared,
                            _dbbt->rankBoxes,
                            PDM_BOX_TREE_ASYNC_LEVEL);

    _update_bt_statistics(&(_dbbt->btsShared), _dbbt->btShared);

    if(0 == 1 && myRank == 0) {
      const char* filename = "dbbt_shared_tree.vtk";
      PDM_vtk_write_boxes (filename,
                           _dbbt->rankBoxes->local_boxes->n_boxes,
                           _dbbt->rankBoxes->local_boxes->extents,
                           _dbbt->rankBoxes->local_boxes->g_num);
    }

    free (allGExtents);
    free (gNumProc);
    free (initLocation_proc);

  }

  /*
   * Build local bt
   */

  _dbbt->btLoc = PDM_box_tree_create (_dbbt->maxTreeDepth,
                                      _dbbt->maxBoxesLeaf,
                                      _dbbt->maxBoxRatio);

  /* Build a tree and associate boxes */

  PDM_box_tree_set_boxes (_dbbt->btLoc,
                          _dbbt->boxes,
                          PDM_BOX_TREE_ASYNC_LEVEL);

  _update_bt_statistics(&(_dbbt->btsLoc), _dbbt->btLoc);
  // PDM_box_tree_dump(_dbbt->btLoc);

  return _dbbt->boxes;

}


/**
 * \brief Assign boxes to intersect to the tree.
 *
 * This function  assigns boxes to intersect to the tree.
 *
 * \param [in]  dbbt       Pointer to a distributed bounding box tree
 * \param [in]  n_part     Number of partitions
 * \param [in]  nElts      Number of elements of each partition
 * \param [in]  extents    Extents of each element of each partition
 * \param [in]  gNum       Global number of each element of each partition
 * \param [out] box_index  Pointer to the index array on associated tree bounding boxeq
 * \param [out] box_g_num  Pointer to the list of intersecting bounding boxes
 *
 * \return associated \ref PDM_box_set_t structure distributed according
 * to the tree intersection
 *
 */

PDM_box_set_t *
PDM_dbbtree_intersect_boxes_with_init_location_set
(
 PDM_dbbtree_t      *dbbt,
 const int           n_part,
 const int          *nElts,
 const int         **init_location,
 const double      **extents,
 const PDM_g_num_t **gNum,
 int                *box_index[],
 int                *box_l_num[]
)
{
  assert (dbbt != NULL);
  _PDM_dbbtree_t *_dbbt = (_PDM_dbbtree_t *) dbbt;

  int myRank;
  PDM_MPI_Comm_rank (_dbbt->comm, &myRank);
  int lComm;
  PDM_MPI_Comm_size (_dbbt->comm, &lComm);

  const int nInfoLocation = 3;
  const int sExtents = _dbbt->dim * 2;

  int nEltsProc = 0;
  for (int i = 0; i < n_part; i++) {
    nEltsProc += nElts[i];
  }

  PDM_g_num_t *_boxGnum = (PDM_g_num_t *) malloc (sizeof(PDM_g_num_t) * nEltsProc);
  double *_extents     = (double *) malloc (sizeof(double) * nEltsProc * sExtents);
  int *_initLocation   = (int *) malloc (sizeof(int) * nEltsProc * nInfoLocation);

  int idx = 0;
  int idx1 = 0;
  int idx2 = 0;

  if(init_location == NULL) {
    for (int i = 0; i < n_part; i++) {

      for (int j = 0; j < nElts[i]; j++) {
        _boxGnum[idx++] = gNum[i][j];

        for (int k = 0; k < sExtents; k++) {
          _extents[idx1++] = extents[i][sExtents*j+k];
        }

        _initLocation[idx2++] = myRank;
        _initLocation[idx2++] = i;
        _initLocation[idx2++] = j;
      }
    }
  } else {
    for (int i = 0; i < n_part; i++) {

      for (int j = 0; j < nElts[i]; j++) {
        _boxGnum[idx++] = gNum[i][j];

        for (int k = 0; k < sExtents; k++) {
          _extents[idx1++] = extents[i][sExtents*j+k];
        }

        _initLocation[idx2++] = init_location[i][3*j  ];
        _initLocation[idx2++] = init_location[i][3*j+1];
        _initLocation[idx2++] = init_location[i][3*j+2];
      }
    }
  }
  if (1 == 0) {
    PDM_printf ("_extents m1 :\n");
    idx1 = 0;
    for (int i = 0; i < nEltsProc; i++) {
      PDM_printf (" "PDM_FMT_G_NUM":", _boxGnum[i]);
      PDM_printf (" %12.5e", _extents[idx1++]);
      PDM_printf (" %12.5e", _extents[idx1++]);
      PDM_printf (" %12.5e", _extents[idx1++]);
      PDM_printf (" %12.5e", _extents[idx1++]);
      PDM_printf (" %12.5e", _extents[idx1++]);
      PDM_printf (" %12.5e", _extents[idx1++]);
      PDM_printf ("\n");
    }
    PDM_printf ("\n");

    PDM_printf ("_initLocation :");
    idx1 = 0;
    for (int i = 0; i < nEltsProc; i++) {
      PDM_printf (" %d", _initLocation[idx1++]);
      PDM_printf (" %d", _initLocation[idx1++]);
      PDM_printf (" %d", _initLocation[idx1++]);
      PDM_printf ("\n");
    }
    PDM_printf ("\n");
  }

  for (int i = 0; i < nEltsProc; i++) {
    _normalize (_dbbt,
                _extents+2*i*_dbbt->boxes->dim,
                _extents+2*i*_dbbt->boxes->dim);
    _normalize (_dbbt,
                _extents+(2*i+1)*_dbbt->boxes->dim,
                _extents+(2*i+1)*_dbbt->boxes->dim);
  }

  if (1 == 0) {

    PDM_printf ("nEltsProc : %d\n", nEltsProc);

    PDM_printf ("_boxGnum :");
    for (int i = 0; i < nEltsProc; i++) {
      PDM_printf (" "PDM_FMT_G_NUM, _boxGnum[i]);
    }
    PDM_printf ("\n");

    PDM_printf ("_extents m2 :\n");
    idx1 = 0;
    for (int i = 0; i < nEltsProc; i++) {
      PDM_printf (" "PDM_FMT_G_NUM":", _boxGnum[i]);
      PDM_printf (" %12.5e", _extents[idx1++]);
      PDM_printf (" %12.5e", _extents[idx1++]);
      PDM_printf (" %12.5e", _extents[idx1++]);
      PDM_printf (" %12.5e", _extents[idx1++]);
      PDM_printf (" %12.5e", _extents[idx1++]);
      PDM_printf (" %12.5e", _extents[idx1++]);
      PDM_printf ("\n");
    }
    PDM_printf ("\n");

    PDM_printf ("_initLocation :");
    idx1 = 0;
    for (int i = 0; i < nEltsProc; i++) {
      PDM_printf (" %d", _initLocation[idx1++]);
      PDM_printf (" %d", _initLocation[idx1++]);
      PDM_printf (" %d", _initLocation[idx1++]);
      PDM_printf ("\n");
    }
    PDM_printf ("\n");
  }

  PDM_box_set_t  *boxes = PDM_box_set_create (3,
                                              0,  // No normalization to preserve initial extents
                                              0,  // No projection to preserve initial extents
                                              nEltsProc,
                                              _boxGnum,
                                              _extents,
                                              n_part,
                                              nElts,
                                              _initLocation,
                                              _dbbt->comm);


  /*
   * Intersection boxes whith shared tree
   */

  if (_dbbt->btShared != NULL) {

    PDM_box_tree_get_boxes_intersects (_dbbt->btShared,
                                       boxes,
                                       box_index,
                                       box_l_num);

    int nUsedRank = PDM_box_set_get_size (_dbbt->rankBoxes);
    const int *usedRanks = _dbbt->usedRank;

    /*
     * Distribute boxes on intersection ranks
     */

    if (1 == 0){
      PDM_printf ("box_l_num_shared : \n");
      for (int i = 0; i < nUsedRank; i++) {
        printf("[%d] : %12.5e %12.5e %12.5e %12.5e %12.5e %12.5e\n",
               i, _dbbt->rankBoxes->local_boxes->extents[6*i]
               , _dbbt->rankBoxes->local_boxes->extents[6*i+1]
               , _dbbt->rankBoxes->local_boxes->extents[6*i+2]
               , _dbbt->rankBoxes->local_boxes->extents[6*i+3]
               , _dbbt->rankBoxes->local_boxes->extents[6*i+4]
               , _dbbt->rankBoxes->local_boxes->extents[6*i+5]);
        printf("[%d] : ",i);
        for (int j = (*box_index)[i]; j < (*box_index)[i+1]; j++) {
          PDM_printf (" "PDM_FMT_G_NUM"",  boxes->local_boxes->g_num[(*box_l_num)[j]]);
        }
        PDM_printf ("\n");
      }
    }

    PDM_box_distrib_t  *distrib = NULL;
    distrib = PDM_box_distrib_shared_create (boxes->local_boxes->n_boxes,//*
                                             boxes->n_g_boxes,
                                             1, // Don't use in this case
                                             boxes->comm);

    PDM_array_reset_int(distrib->index, lComm + 1, 0);

    for (int i = 0; i < nUsedRank; i++) {
      distrib->index[usedRanks[i]+1] = (*box_index)[i+1] - (*box_index)[i];
    }

    for (int i = 0; i < lComm; i++) {
      distrib->index[i+1] = distrib->index[i+1] + distrib->index[i];
    }

    distrib->list = (int *) malloc (sizeof(int) * distrib->index[lComm]);

    for (int i = 0; i < distrib->index[lComm]; i++) {
      distrib->list[i] = (*box_l_num)[i];
    }

    /*
     * Redistribute boxes on intersecting ranks
     */

    PDM_box_set_redistribute (distrib, boxes);

    if (1 == 0) {
      printf ("Boxes B apres redistribution : %d\n", boxes->local_boxes->n_boxes);
      for (int i = 0; i < boxes->local_boxes->n_boxes; i++) {
        printf (" "PDM_FMT_G_NUM"", boxes->local_boxes->g_num[i]);
      }
      printf("\n");
    }

    /*
     * Free
     */

    PDM_box_distrib_destroy (&distrib);

    free (*box_l_num);
    free (*box_index);

  }

  free (_boxGnum);
  free (_extents);
  free (_initLocation);

  /*
   * Intersection boxes whith local tree
   */

  //PDM_box_tree_dump_statistics(_dbbt->btLoc);

  PDM_box_tree_get_boxes_intersects (_dbbt->btLoc,
                                     boxes,
                                     box_index,
                                     box_l_num);

  /*
   * Sort boxes and remove double boxes
   */

  int nBoxesA = _dbbt->boxes->local_boxes->n_boxes;//*

  int *newIndex = PDM_array_zeros_int(nBoxesA + 1);

  int *_box_index = *box_index;
  int *_box_l_num = *box_l_num;

  idx = 0;
  for (int i = 0; i < nBoxesA; i++) {

    int *ideb = _box_l_num + _box_index[i];
    int length = _box_index[i+1] - _box_index[i];

    PDM_sort_int (ideb, NULL, length);

    int pre = -1;
    for (int k = 0; k < length; k++) {
      if (pre != ideb[k]) {
        (*box_l_num)[idx++] = ideb[k];
        newIndex[i+1] += 1;
        pre = ideb[k];
      }
    }
  }

  for (int i = 0; i < nBoxesA; i++) {
    newIndex[i+1] += newIndex[i];
  }

  free (*box_index);
  *box_index = newIndex;

  *box_l_num = (int *) realloc (*box_l_num, sizeof (int) * newIndex[nBoxesA]);

  if (1 == 0) {
    printf ("Intersections : %d\n", boxes->local_boxes->n_boxes);
    for (int i = 0; i < nBoxesA; i++) {
      printf ("A elt "PDM_FMT_G_NUM" :", _dbbt->boxes->local_boxes->g_num[i]);
      for (int j = newIndex[i]; j < newIndex[i+1]; j++) {
        printf (" "PDM_FMT_G_NUM"", boxes->local_boxes->g_num[(*box_l_num)[j]]);
      }
      printf("\n");
    }
  }
  return boxes;
}

PDM_box_set_t *
PDM_dbbtree_intersect_boxes_set
(
 PDM_dbbtree_t      *dbbt,
 const int           n_part,
 const int          *nElts,
 const double      **extents,
 const PDM_g_num_t **gNum,
 int                *box_index[],
 int                *box_l_num[]
)
{
  return PDM_dbbtree_intersect_boxes_with_init_location_set(dbbt, n_part, nElts, NULL, extents, gNum, box_index, box_l_num);
}



/**
 *
 * \brief Get the boxes closer than the upper bound distance
 *
 *   \param [in]  bt                 Pointer to box tree structure
 *   \param [in]  n_pts              Number of points
 *   \param [in]  pts                Point coordinates (size = 3 * \ref n_pts)
 *   \param [in]  pts_g_num          Point global ids
 *   \param [in]  upper_bound_dist2  Upper bound of the square of the distance (size = \ref n_pts)
 *   \param [out] box_index          Index of boxes (size = \ref n_pts + 1)
 *   \param [out] box_g_num          Global ids of boxes (size = \ref i_boxes[\ref n_pts])
 *
 */

void
PDM_dbbtree_closest_upper_bound_dist_boxes_get
(
 PDM_dbbtree_t   *dbbt,
 const int        n_pts,
 double           pts[],
 PDM_g_num_t      pts_g_num[],
 double           upper_bound_dist2[],
 int             *box_index[],
 PDM_g_num_t     *box_g_num[]
)
{
  PDM_UNUSED(pts_g_num);
  /*
   * RANK DATA COPY PARAMETERS
   */
  const double RANK_COPY_threshold  = 1.2;  // factor of the mean nb of requests
  const double RANK_COPY_max_copies = 0.15; // factor of the total nb of processes

  /*
   * Initialization
   */
  assert (dbbt != NULL);
  _PDM_dbbtree_t *_dbbt = (_PDM_dbbtree_t *) dbbt;

  double *_pts = pts;
  if (_dbbt->global_extents != NULL) {
    _pts = malloc (sizeof(double) * n_pts * 3);

    for (int i = 0; i < n_pts; i++) {
      _normalize (_dbbt,
                  pts + 3*i,
                  _pts + 3*i);
    }
  }

  int     n_pts_local            = n_pts;
  double *pts_local              = _pts;
  double *upper_bound_dist_local = upper_bound_dist2;


  int myRank;
  PDM_MPI_Comm_rank (_dbbt->comm, &myRank);
  int lComm;
  PDM_MPI_Comm_size (_dbbt->comm, &lComm);

  /*
   * Determine for each point the list of involved processes
   */
  int *n_send_pts = NULL;
  int *n_recv_pts = NULL;

  int *box_index_tmp = NULL;
  int *box_l_num_tmp = NULL;

  int *n_pts_rank = NULL;
  int *n_pts_send = NULL;
  int *n_pts_recv = NULL;
  int n_pts_recv_total = 0;

  int *i_pts_rank = NULL;
  int *i_pts_send = NULL;
  int *i_pts_recv = NULL;

  double *pts_rank              = NULL;
  double *upper_bound_dist_rank = NULL;
  double *pts_recv              = NULL;
  double *upper_bound_dist_recv = NULL;

  int n_copied_ranks = 0;
  int *copied_ranks  = NULL;
  int *rank_copy_num = NULL;

  const int *usedRanks = _dbbt->usedRank;

  const int idebug = 0;

  int i1 = 0, i2 = 0, i3 = 0;
  int i_rank = 0;
  if (_dbbt->btShared != NULL) {

    if (idebug) {
      printf ("  **** deb PDM_box_tree_closest_upper_bound_dist_boxes_get shared _pts : %d\n", n_pts);
    }

    PDM_box_tree_closest_upper_bound_dist_boxes_get (_dbbt->btShared,
                                                     n_pts,
                                                     pts,
                                                     upper_bound_dist2,
                                                     &box_index_tmp,
                                                     &box_l_num_tmp);

    if (idebug) {
      printf ("  **** fin PDM_box_tree_closest_upper_bound_dist_boxes_get shared n_pts : %d\n", n_pts);
      for (int i = 0; i < n_pts; i++) {
        printf ("%d : (%12.5e %12.5e %12.5e) %12.5e\n", i,
                pts[3*i], pts[3*i+1], pts[3*i+2],
                upper_bound_dist2[i]);
        printf ("  boxes %d :" , box_index_tmp[i+1] - box_index_tmp[i]);
        for (int j = box_index_tmp[i]; j < box_index_tmp[i+1]; j++) {
          printf (" %d", box_l_num_tmp[j]);
        }
        printf ("\n");
      }
    }

    /*
     * Count (provisional) nb of points to send to each process
     */
    n_send_pts = PDM_array_zeros_int(lComm);
    n_recv_pts = malloc (sizeof(int) * lComm);


    for (int i = 0; i < n_pts; i++) {
      for (int j = box_index_tmp[i]; j < box_index_tmp[i+1]; j++) {
        n_send_pts[usedRanks[box_l_num_tmp[j]]]++;
      }
    }

    PDM_MPI_Alltoall (n_send_pts, 1, PDM_MPI_INT,
                      n_recv_pts, 1, PDM_MPI_INT,
                      _dbbt->comm);
    free(n_send_pts);

    /*
     * Prepare copies
     */
    // total nb of requests received by current process
    int local_sum_nrecv = 0;
    for (int i = 0; i < lComm; i++) {
      local_sum_nrecv += n_recv_pts[i];
    }
    free(n_recv_pts);

    int *n_requests = malloc (lComm * sizeof(int));
    PDM_MPI_Allgather (&local_sum_nrecv, 1, PDM_MPI_INT,
                       n_requests,       1, PDM_MPI_INT,
                       _dbbt->comm);

    // mean nb of requests
    int mean_n_requests = 0;
    for (int i = 0; i < lComm; i++) {
      mean_n_requests += n_requests[i];
    }
    mean_n_requests /= lComm;

    /* sort the ranks in ascending order of
     * the total nb of points they are supposed to receive */
    int *order = malloc (lComm * sizeof(int));
    for (int i = 0; i < lComm; i ++) {
      order[i] = i;
    }

    PDM_sort_int (n_requests, order, lComm);

    // identify ranks to be copied
    double threshold_n_req = RANK_COPY_threshold*mean_n_requests;
    int max_copied_ranks   = (int) _MAX (1, RANK_COPY_max_copies*lComm);

    n_copied_ranks = 0;
    copied_ranks = malloc (max_copied_ranks * sizeof(int));

    for (int i = 0; i < max_copied_ranks; i++) {
      i_rank = lComm - 1 - i;

      if ( n_requests[i_rank] > threshold_n_req ) {
        copied_ranks[n_copied_ranks++] = order[i_rank];
      } else {
        break;
      }
    }
    free(order);
    free(n_requests);

    //------------->>>
    if ( idebug && myRank == 0 ) {
      if ( n_copied_ranks == 0 ) {
        printf("n_copied_ranks = 0\n");
      } else {
        printf("copied rank(s) = ");
        for (int i = 0; i < n_copied_ranks; i++) {
          printf("%d ", copied_ranks[i]);
        }
        printf("\n");
      }
    }
    //<<<-------------

    /*
     * Copy the data of selected ranks
     */
    rank_copy_num = (int *) malloc (sizeof(int) * lComm);
    PDM_box_tree_copy_to_ranks (_dbbt->btLoc,
                                &n_copied_ranks,
                                copied_ranks,
                                rank_copy_num);
    /* rank_copy_num[_dbbt->btLoc->copied_ranks[i]] (def)= i*/
    free(copied_ranks);


    /*
     * Distribution of points...
     *    ..._local --> search in local box tree
     *    ..._rank  --> search in copied box trees
     *    ..._send  --> search in distant box trees (send to other processes)
     */
    n_pts_local = 0;

    n_pts_rank = PDM_array_zeros_int(n_copied_ranks);

    n_pts_send = PDM_array_zeros_int(lComm);
    n_pts_recv = malloc (sizeof(int) * lComm);


    for (int i = 0; i < n_pts; i++) {
      for (int j = box_index_tmp[i]; j < box_index_tmp[i+1]; j++) {
        i_rank = usedRanks[box_l_num_tmp[j]];
        if ( i_rank == myRank ) {
          // ---> search in btLoc->local_data of current process
          n_pts_local++;
        } else if ( rank_copy_num[i_rank] >= 0 ) {
          // ---> search in btLoc->rank_data[rank_copy_num[i_rank]] of current process
          n_pts_rank[rank_copy_num[i_rank]]++;
        } else {
          // ---> search in btLoc->local_data of process with rank i_rank
          n_pts_send[i_rank]++;
        }
      }
    }

    PDM_MPI_Alltoall (n_pts_send, 1, PDM_MPI_INT,
                      n_pts_recv, 1, PDM_MPI_INT,
                      _dbbt->comm);

    i_pts_rank = PDM_array_new_idx_from_sizes_int(n_pts_rank, n_copied_ranks);

    i_pts_send = malloc (sizeof(int) * (lComm+1));
    i_pts_recv = malloc (sizeof(int) * (lComm+1));
    i_pts_send[0] = 0;
    i_pts_recv[0] = 0;
    for (int i = 0; i < lComm; i++) {
      i_pts_send[i+1] = i_pts_send[i] + 4 * n_pts_send[i];
      i_pts_recv[i+1] = i_pts_recv[i] + 4 * n_pts_recv[i];
      n_pts_recv[i] *= 4;
    }

    pts_local              = malloc (sizeof(double) * n_pts_local*3);
    upper_bound_dist_local = malloc (sizeof(double) * n_pts_local);

    pts_rank              = malloc (sizeof(double) * i_pts_rank[n_copied_ranks]*3);
    upper_bound_dist_rank = malloc (sizeof(double) * i_pts_rank[n_copied_ranks]);

    double *data_send = malloc (sizeof(double) * i_pts_send[lComm]);
    double *data_recv = malloc (sizeof(double) * i_pts_recv[lComm]);

    PDM_array_reset_int(n_pts_send, lComm, 0);
    PDM_array_reset_int(n_pts_rank, n_copied_ranks, 0);

    i1 = 0; i2 = 0; i3 = 0;
    for (int i = 0; i < n_pts; i++) {
      for (int j = box_index_tmp[i]; j < box_index_tmp[i+1]; j++) {
        i_rank = usedRanks[box_l_num_tmp[j]];
        if ( i_rank == myRank ) {
          // pts_local, upper_bound_dist_local (local points, local box tree)
          upper_bound_dist_local[i1] = upper_bound_dist2[i];
          pts_local[3*i1]            = _pts[3*i];
          pts_local[3*i1+1]          = _pts[3*i+1];
          pts_local[3*i1+2]          = _pts[3*i+2];
          i1++;
        } else if ( rank_copy_num[i_rank] >= 0 ) {
          // pts_rank, upper_bound_dist_rank (local points, distant (copied) box trees)
          int j_rank = rank_copy_num[i_rank];
          i2 = i_pts_rank[j_rank] + n_pts_rank[j_rank];
          upper_bound_dist_rank[i2] = upper_bound_dist2[i];
          pts_rank[3*i2]            = _pts[3*i];
          pts_rank[3*i2+1]          = _pts[3*i+1];
          pts_rank[3*i2+2]          = _pts[3*i+2];
          n_pts_rank[j_rank]++;
        } else {
          // data_send (local points, distant (not copied) box trees)
          i3 = i_pts_send[i_rank] + 4*n_pts_send[i_rank];
          data_send[i3++] = _pts[3*i];
          data_send[i3++] = _pts[3*i+1];
          data_send[i3++] = _pts[3*i+2];
          data_send[i3++] = upper_bound_dist2[i];

          n_pts_send[i_rank]++;
        }
      }
    }

    for (int i = 0; i < lComm; i++) {
      n_pts_send[i] *= 4;
    }


    // Send points to search in distant (not copied) box trees
    PDM_MPI_Alltoallv (data_send, n_pts_send, i_pts_send, PDM_MPI_DOUBLE,
                       data_recv, n_pts_recv, i_pts_recv, PDM_MPI_DOUBLE,
                       _dbbt->comm);
    free(data_send);

    n_pts_recv_total = i_pts_recv[lComm] / 4;

    pts_recv              = (double *) malloc (sizeof(double) * n_pts_recv_total * 3);
    upper_bound_dist_recv = (double *) malloc (sizeof(double) * n_pts_recv_total);

    for (int i = 0; i < n_pts_recv_total; i++) {
      for (int j = 0; j < 3; j++) {
        pts_recv[3*i+j] = data_recv[4*i+j];
      }
      upper_bound_dist_recv[i] = data_recv[4*i+3];
    }
    free(data_recv);
  }
  if (_pts != pts && _pts != pts_local) {
    free (_pts);
  }


  // Determine candidate boxes in local box tree (local points)
  int *box_index_local;
  int *box_l_num_local;
  PDM_box_tree_closest_upper_bound_dist_boxes_get_v2 (_dbbt->btLoc,
                                                      -1, // search in local box tree
                                                      n_pts_local,
                                                      pts_local,
                                                      upper_bound_dist_local,
                                                      &box_index_local,
                                                      &box_l_num_local,
                                                      _dbbt->d);
  if (pts_local != pts) {
    free(pts_local);
  }
  if (upper_bound_dist_local != upper_bound_dist2) {
    free(upper_bound_dist_local);
  }

  // conversion local --> global numbering
  const PDM_g_num_t *gnum_boxes_local = PDM_box_set_get_g_num (_dbbt->boxes);
  PDM_g_num_t *box_g_num_local = malloc(sizeof(PDM_g_num_t) * box_index_local[n_pts_local]);

  for (int i = 0; i < box_index_local[n_pts_local]; i++) {
    box_g_num_local[i] = gnum_boxes_local[box_l_num_local[i]];
  }
  free(box_l_num_local);


  if (_dbbt->btShared == NULL) {

    *box_index = box_index_local;
    *box_g_num = box_g_num_local;

  } else {
    // Determine candidate boxes in copied box trees (local points)
    int **box_index_rank;
    int **box_l_num_rank;

    box_index_rank = (int **) malloc (sizeof(int *) * n_copied_ranks);
    box_l_num_rank = (int **) malloc (sizeof(int *) * n_copied_ranks);

    for (int i_copied_rank = 0; i_copied_rank < n_copied_ranks; i_copied_rank++) {
      int n_pts_copied_rank = i_pts_rank[i_copied_rank+1] - i_pts_rank[i_copied_rank];
      double *pts_copied_rank = pts_rank + 3*i_pts_rank[i_copied_rank];
      double *upper_bound_dist_copied_rank = upper_bound_dist_rank + i_pts_rank[i_copied_rank];
      PDM_box_tree_closest_upper_bound_dist_boxes_get_v2 (_dbbt->btLoc,
                                                          i_copied_rank,
                                                          n_pts_copied_rank,
                                                          pts_copied_rank,
                                                          upper_bound_dist_copied_rank,
                                                          &(box_index_rank[i_copied_rank]),
                                                          &(box_l_num_rank[i_copied_rank]),
                                                          _dbbt->d);
    }

    free(pts_rank);
    free(upper_bound_dist_rank);

    // conversion local --> global numbering for each copied rank
    PDM_g_num_t **box_g_num_rank = malloc(sizeof(PDM_g_num_t *) * n_copied_ranks);

    PDM_g_num_t *gnum_boxes_rank = NULL;
    for (int i_copied_rank = 0; i_copied_rank < n_copied_ranks; i_copied_rank++) {
      gnum_boxes_rank = PDM_box_set_get_rank_boxes_g_num (_dbbt->boxes,
                                                          i_copied_rank);
      box_g_num_rank[i_copied_rank] = malloc(sizeof(PDM_g_num_t) * box_index_rank[i_copied_rank][n_pts_rank[i_copied_rank]]);
      for (int i = 0; i < box_index_rank[i_copied_rank][n_pts_rank[i_copied_rank]]; i++) {
        box_g_num_rank[i_copied_rank][i] = gnum_boxes_rank[box_l_num_rank[i_copied_rank][i]];
      }

      free(box_l_num_rank[i_copied_rank]);
    }
    free(box_l_num_rank);



    // Determine candidate boxes in local box tree (received points)
    int *n_pts_send2 = NULL;
    int *i_pts_send2 = NULL;
    int *n_box_l_num_per_pts = NULL;
    PDM_g_num_t *box_g_num_per_pts = NULL;


    int *box_index_recv = NULL;
    int *box_l_num_recv = NULL;

    PDM_box_tree_closest_upper_bound_dist_boxes_get_v2 (_dbbt->btLoc,
                                                        -1, // search in local box tree
                                                        n_pts_recv_total,
                                                        pts_recv,
                                                        upper_bound_dist_recv,
                                                        &box_index_recv,
                                                        &box_l_num_recv,
                                                        _dbbt->d);
    free(pts_recv);
    free(upper_bound_dist_recv);

    /*
     * Send back results for distant points to original processes:
     *     - nb of boxes for each point
     *     - global numbering of these boxes
     */

    int *n_box_l_num_recv = malloc (sizeof(int) * n_pts_recv_total);

    for (int i = 0; i < n_pts_recv_total; i++) {
      n_box_l_num_recv[i] = box_index_recv[i+1] - box_index_recv[i];
    }

    for (int i = 0; i < lComm; i++) {
      i_pts_send[i+1] = i_pts_send[i+1]/4;
      i_pts_recv[i+1] = i_pts_recv[i+1]/4;
      n_pts_send[i]   = n_pts_send[i]/4;
      n_pts_recv[i]   = n_pts_recv[i]/4;
    }

    n_box_l_num_per_pts = malloc (sizeof(int) * i_pts_send[lComm]);

    PDM_MPI_Alltoallv (n_box_l_num_recv,    n_pts_recv, i_pts_recv, PDM_MPI_INT,
                       n_box_l_num_per_pts, n_pts_send, i_pts_send, PDM_MPI_INT,
                       _dbbt->comm);

    n_pts_send2 = malloc (sizeof(int) * lComm);
    i_pts_send2 = malloc (sizeof(int) * (lComm+1));

    int *n_pts_recv2 = malloc (sizeof(int) * lComm);
    int *i_pts_recv2 = malloc (sizeof(int) * (lComm+1));

    for (int i = 0; i < lComm; i++) {
      n_pts_send2[i] = 0;
      n_pts_recv2[i] = 0;
    }

    for (int i = 0; i < lComm; i++) {
      for (int j = i_pts_recv[i]; j < i_pts_recv[i+1]; j++) {
        n_pts_recv2[i] += n_box_l_num_recv[j];
      }
      for (int j = i_pts_send[i]; j < i_pts_send[i+1]; j++) {
        n_pts_send2[i] += n_box_l_num_per_pts[j];
      }
    }

    free(n_box_l_num_recv);

    i_pts_send2[0] = 0;
    i_pts_recv2[0] = 0;
    for (int i = 0; i < lComm; i++) {
      i_pts_send2[i+1] = i_pts_send2[i] + n_pts_send2[i];
      i_pts_recv2[i+1] = i_pts_recv2[i] + n_pts_recv2[i];
    }


    // Conversion local --> global numbering
    PDM_g_num_t *box_g_num_recv = malloc(sizeof(PDM_g_num_t) * box_index_recv[n_pts_recv_total]);

    for (int i = 0; i < box_index_recv[n_pts_recv_total]; i++) {
      box_g_num_recv[i] = gnum_boxes_local[box_l_num_recv[i]];
    }

    box_g_num_per_pts = malloc(sizeof(PDM_g_num_t) * i_pts_send2[lComm]);
    PDM_MPI_Alltoallv (box_g_num_recv,    n_pts_recv2, i_pts_recv2, PDM__PDM_MPI_G_NUM,
                       box_g_num_per_pts, n_pts_send2, i_pts_send2, PDM__PDM_MPI_G_NUM,
                       _dbbt->comm);

    free(box_index_recv);
    free(box_l_num_recv);
    free(box_g_num_recv);

    free(n_pts_recv);
    free(i_pts_recv);
    free(n_pts_recv2);
    free(i_pts_recv2);




    /*
     * Merge all results and resolve duplicates
     */
    *box_index = malloc(sizeof(int) * (n_pts + 1));

    int max_n_box_g_num = box_index_local[n_pts_local];
    for (int i_copied_rank = 0; i_copied_rank < n_copied_ranks; i_copied_rank++) {
      max_n_box_g_num += box_index_rank[i_copied_rank][n_pts_rank[i_copied_rank]];
    }
    max_n_box_g_num += i_pts_send2[lComm];

    *box_g_num = malloc(sizeof(PDM_g_num_t) * max_n_box_g_num);


    int *rank_index = (int *) malloc(sizeof(int) * (n_pts+1));
    memcpy(rank_index, box_index_tmp, sizeof(int) * (n_pts+1));

    for (int i = 0; i < lComm; i++) {
      n_pts_send[i] = 0;
      n_pts_send2[i] = 0;
    }

    PDM_array_reset_int(n_pts_rank, n_copied_ranks, 0);



    int keyMax = 3 * n_pts;
    int key = 0;
    int found = 0;
    PDM_hash_tab_t *ht = PDM_hash_tab_create (PDM_HASH_TAB_KEY_INT,
                                              &keyMax);

    PDM_g_num_t i_box = 0;

    box_index_tmp[0] = 0;
    int idx = 0;
    i1 = 0; i2 = 0; i3 = 0;
    if ( 1 ) {
      for (int i = 0; i < n_pts; i++) { // loop over local points
        box_index_tmp[i+1] = box_index_tmp[i];
        for (int j = rank_index[i]; j < rank_index[i+1]; j++) { // loop over procs to which the current point was sent
          i_rank = usedRanks[box_l_num_tmp[j]]; // i_rank = rank j-th proc to which the current point was sent

          if ( i_rank == myRank ) { // boxes local to current process
            for (int k = box_index_local[i1]; k < box_index_local[i1+1]; k++) {
              i_box = box_g_num_local[k];

              found = PDM_hash_tab_check_collision (ht, i_box, keyMax, &key);

              if (!found) {
                PDM_hash_tab_data_add (ht, (void *) &key, *box_g_num + idx);
                (*box_g_num)[idx++] = i_box;
                box_index_tmp[i+1] += 1;
              }

            }
            i1++;

          } else if ( rank_copy_num[i_rank] >= 0 ) { // distant boxes copied in current process
            int j_rank = rank_copy_num[i_rank];
            i2 = n_pts_rank[j_rank];

            for (int k = box_index_rank[j_rank][i2]; k < box_index_rank[j_rank][i2+1]; k++) {
              i_box = box_g_num_rank[j_rank][k];

              found = PDM_hash_tab_check_collision (ht, i_box, keyMax, &key);

              if (!found) {
                PDM_hash_tab_data_add (ht, (void *) &key, *box_g_num + idx);
                (*box_g_num)[idx++] = i_box;
                box_index_tmp[i+1] += 1;
              }

            }
            n_pts_rank[j_rank]++;

          } else { // distant boxes (not copied)
            i3 = n_pts_send[i_rank];
            int i4 = i_pts_send2[i_rank] + n_pts_send2[i_rank];
            int i5 = i_pts_send[i_rank] + i3;
            for (int k = 0; k < n_box_l_num_per_pts[i5]; k++) {
              i_box = box_g_num_per_pts[i4++];

              found = PDM_hash_tab_check_collision (ht, i_box, keyMax, &key);

              if (!found) {
                PDM_hash_tab_data_add (ht, (void *) &key, *box_g_num + idx);
                (*box_g_num)[idx++] = i_box;
                box_index_tmp[i+1] += 1;
              }

            }
            n_pts_send2[i_rank] += n_box_l_num_per_pts[i5];
            n_pts_send[i_rank]++;

          }
        }
        PDM_hash_tab_purge (ht, PDM_FALSE);
      }
    }
    PDM_hash_tab_free (ht);

    free (n_box_l_num_per_pts);


    if ( *box_index != NULL ) {
      free (*box_index);
    }
    *box_index = box_index_tmp;


    *box_g_num = realloc (*box_g_num, sizeof(PDM_g_num_t) * box_index_tmp[n_pts]);

    /*
     * Deallocate stuff
     */
    free(box_l_num_tmp);
    free(box_g_num_per_pts);

    free(box_index_local);
    free(box_g_num_local);

    for (int i_copied_rank = 0; i_copied_rank < n_copied_ranks; i_copied_rank++) {
      free(box_index_rank[i_copied_rank]);
      free(box_g_num_rank[i_copied_rank]);
    }
    free(box_index_rank);
    free(box_g_num_rank);


    free(n_pts_send);
    free(i_pts_send);
    free(n_pts_send2);
    free(i_pts_send2);

    free(i_pts_rank);
    free(n_pts_rank);

    free(rank_index);
    free(rank_copy_num);
  }

  PDM_box_tree_free_copies(_dbbt->btLoc);
}

/**
 *
 * \brief Get the boxes closer than the upper bound distance
 *
 *   \param [in]  bt                 Pointer to box tree structure
 *   \param [in]  n_pts              Number of points
 *   \param [in]  pts                Point coordinates (size = 3 * \ref n_pts)
 *   \param [in]  pts_g_num          Point global ids
 *   \param [in]  upper_bound_dist2  Upper bound of the square of the distance (size = \ref n_pts)
 *   \param [out] n_extract_boxes    Number of extracted box
 *   \param [out] box_l_num          Index of boxes (size = \ref n_extract_boxes )
 *   \param [out] box_pts_idx        Index of boxes (size = \ref n_extract_boxes + 1)
 *   \param [out] box_g_num          Global ids of boxes (size = \ref box_pts_idx[\ref n_pts])
 *
 */
void
PDM_dbbtree_closest_upper_bound_dist_boxes_pts_shared_get
(
 PDM_dbbtree_t   *dbbt,
 const int        n_pts,
 double           pts[],
 PDM_g_num_t      pts_g_num[],
 double           upper_bound_dist2[],
 int             *out_n_extract_boxes,
 PDM_g_num_t     *out_box_gnum[],
 int             *out_box_init_location[],
 int             *out_dbox_pts_idx[],
 PDM_g_num_t     *out_dbox_pts_g_num[],
 double          *out_dbox_pts_coord[]
)
{

  assert (dbbt != NULL);
  _PDM_dbbtree_t *_dbbt = (_PDM_dbbtree_t *) dbbt;

  double *_pts = pts;
  if (_dbbt->global_extents != NULL) {
    _pts = malloc (sizeof(double) * n_pts * 3);

    for (int i = 0; i < n_pts; i++) {
      _normalize (_dbbt, pts + 3*i, _pts + 3*i);
    }
  }

  if(_dbbt->btLoc->shm_data == NULL) {
    PDM_box_tree_copy_to_shm(_dbbt->btLoc);
  }

  int i_rank;
  PDM_MPI_Comm_rank (_dbbt->comm, &i_rank);
  int n_rank;
  PDM_MPI_Comm_size (_dbbt->comm, &n_rank);

  PDM_MPI_Datatype mpi_pts_coords_type;
  PDM_MPI_Type_create_contiguous(3, PDM_MPI_DOUBLE, &mpi_pts_coords_type);
  PDM_MPI_Type_commit(&mpi_pts_coords_type);

  // Shared
  PDM_MPI_Comm comm_shared;
  PDM_MPI_Comm_split_type(_dbbt->comm, PDM_MPI_SPLIT_NUMA, &comm_shared);

  int n_rank_in_shm, i_rank_in_shm;
  PDM_MPI_Comm_rank (comm_shared, &i_rank_in_shm);
  PDM_MPI_Comm_size (comm_shared, &n_rank_in_shm);

  PDM_mpi_win_shared_t* wshared_recv_gnum              = NULL;
  PDM_mpi_win_shared_t* wshared_recv_pts_coord         = NULL;
  PDM_mpi_win_shared_t* wshared_recv_upper_bound_dist2 = NULL;

  int          n_pts1 = 0;
  PDM_g_num_t *pts_g_num1            = NULL;
  double      *pts_coord1            = NULL;
  double      *pts_upper_bound_dist2 = NULL;

  int *distrib_search_by_rank_idx = NULL;


  if (_dbbt->btShared != NULL) {
    int *pts_rank_idx = NULL;
    int *pts_rank     = NULL;

    int *send_count = NULL;
    int *send_shift = NULL;
    int *recv_count = NULL;
    int *recv_shift = NULL;

    PDM_g_num_t *send_g_num             = NULL;
    double      *send_coord             = NULL;
    double      *send_upper_bound_dist2 = NULL;

    PDM_box_tree_closest_upper_bound_dist_boxes_get (_dbbt->btShared,
                                                     n_pts,
                                                     pts,
                                                     upper_bound_dist2,
                                                     &pts_rank_idx,
                                                     &pts_rank);

    /*
     * Count (provisional) nb of points to send to each process
     */
    send_count = PDM_array_zeros_int(n_rank);
    recv_count = malloc (sizeof(int) * n_rank);

    for (int i = 0; i < pts_rank_idx[n_pts]; i++) {
      int t_rank = _dbbt->usedRank[pts_rank[i]];
      pts_rank[i] = t_rank;
      send_count[t_rank]++;
    }

    PDM_MPI_Alltoall (send_count, 1, PDM_MPI_INT,
                      recv_count, 1, PDM_MPI_INT,
                      _dbbt->comm);

    send_shift = malloc ( ( n_rank + 1) * sizeof(int));
    recv_shift = malloc ( ( n_rank + 1) * sizeof(int));

    // Deduce size of recv buffer shared inside the same node
    int* shared_recv_count  = malloc(n_rank_in_shm * sizeof(int));

    int n_tot_recv = 0;
    send_shift[0] = 0;
    recv_shift[0] = 0;
    for(int i = 0; i < n_rank; ++i) {
      n_tot_recv += recv_count[i];
      send_shift[i+1] = send_shift[i] + send_count[i];
      recv_shift[i+1] = recv_shift[i] + recv_count[i];
    }

    PDM_MPI_Allgather(&n_tot_recv      , 1, PDM_MPI_INT,
                      shared_recv_count, 1, PDM_MPI_INT,
                      comm_shared);


    int* shared_recv_idx = malloc((n_rank_in_shm+1) * sizeof(int));
    shared_recv_idx[0] = 0;
    for(int i = 0; i < n_rank_in_shm; ++i) {
      shared_recv_idx[i+1] = shared_recv_idx[i] + shared_recv_count[i];
    }

    if(0 == 1) {
      PDM_log_trace_array_int(shared_recv_count , n_rank_in_shm, "shared_recv_count  :: ");
      PDM_log_trace_array_int(shared_recv_idx , n_rank_in_shm+1, "shared_recv_idx  :: ");
    }

    int n_tot_recv_shared = shared_recv_idx[n_rank_in_shm];
    wshared_recv_gnum              = PDM_mpi_win_shared_create(    n_tot_recv_shared, sizeof(PDM_g_num_t), comm_shared);
    wshared_recv_pts_coord         = PDM_mpi_win_shared_create(3 * n_tot_recv_shared, sizeof(double     ), comm_shared);
    wshared_recv_upper_bound_dist2 = PDM_mpi_win_shared_create(    n_tot_recv_shared, sizeof(double     ), comm_shared);
    PDM_MPI_Barrier(comm_shared);


    PDM_g_num_t *shared_recv_gnum              = PDM_mpi_win_shared_get(wshared_recv_gnum);
    double      *shared_recv_pts_coord         = PDM_mpi_win_shared_get(wshared_recv_pts_coord);
    double      *shared_recv_upper_bound_dist2 = PDM_mpi_win_shared_get(wshared_recv_upper_bound_dist2);

    PDM_mpi_win_shared_lock_all (0, wshared_recv_gnum);
    PDM_mpi_win_shared_lock_all (0, wshared_recv_pts_coord);
    PDM_mpi_win_shared_lock_all (0, wshared_recv_upper_bound_dist2);

    PDM_g_num_t *lrecv_gnum              = &shared_recv_gnum             [    shared_recv_idx[i_rank_in_shm]];
    double      *lrecv_pts_coord         = &shared_recv_pts_coord        [3 * shared_recv_idx[i_rank_in_shm]];
    double      *lrecv_upper_bound_dist2 = &shared_recv_upper_bound_dist2[    shared_recv_idx[i_rank_in_shm]];


    /* Prepare send */
    send_g_num             = malloc(    send_shift[n_rank] * sizeof(PDM_g_num_t));
    send_coord             = malloc(3 * send_shift[n_rank] * sizeof(double     ));
    send_upper_bound_dist2 = malloc(    send_shift[n_rank] * sizeof(double     ));

    for(int i = 0; i < n_rank; ++i) {
      send_count[i] = 0;
    }


    for (int ipt = 0; ipt < n_pts; ipt++) {
      for (int i = pts_rank_idx[ipt]; i < pts_rank_idx[ipt+1]; i++) {
        int t_rank = pts_rank[i];

        int idx_write = send_shift[t_rank] + send_count[t_rank]++;
        send_g_num            [idx_write] = pts_g_num        [ipt];
        send_upper_bound_dist2[idx_write] = upper_bound_dist2[ipt];
        for (int j = 0; j < 3; j++) {
          send_coord[3*idx_write + j] = _pts[3*ipt + j];
        }
      }
    }
    free(pts_rank_idx);
    free(pts_rank);

    PDM_MPI_Alltoallv (send_g_num, send_count, send_shift, PDM__PDM_MPI_G_NUM,
                       lrecv_gnum, recv_count, recv_shift, PDM__PDM_MPI_G_NUM,
                       _dbbt->comm);
    free(send_g_num);

    PDM_MPI_Alltoallv (send_upper_bound_dist2 , send_count, send_shift, PDM_MPI_DOUBLE,
                       lrecv_upper_bound_dist2, recv_count, recv_shift, PDM_MPI_DOUBLE,
                       _dbbt->comm);
    free(send_upper_bound_dist2);


    PDM_MPI_Alltoallv (send_coord     , send_count, send_shift, mpi_pts_coords_type,
                       lrecv_pts_coord, recv_count, recv_shift, mpi_pts_coords_type,
                       _dbbt->comm);
    free(send_coord);


    free (send_count);
    free (send_shift);
    free (recv_count);
    free (recv_shift);

    PDM_MPI_Barrier (comm_shared);
    PDM_mpi_win_shared_sync (wshared_recv_gnum);
    PDM_mpi_win_shared_sync (wshared_recv_pts_coord);
    PDM_mpi_win_shared_sync (wshared_recv_upper_bound_dist2);

    /*
     * Redistribution by numa
     */

    PDM_g_num_t* distrib_search = PDM_compute_uniform_entity_distribution(comm_shared, n_tot_recv_shared);

    int  dn_search = distrib_search[i_rank_in_shm+1] - distrib_search[i_rank_in_shm];

    distrib_search_by_rank_idx = malloc((n_rank_in_shm+1) * sizeof(int));
    int* distrib_search_by_rank_n   = malloc((n_rank_in_shm  ) * sizeof(int));
    for(int i = 0; i < n_rank_in_shm; ++i) {
      distrib_search_by_rank_n[i] = 0;
    }

    // TODO : Faire un algo d'intersection de range pour ne pas faire la dicotomie x fois !
    for(int i = distrib_search[i_rank_in_shm]; i < distrib_search[i_rank_in_shm+1]; ++i) {
      int t_rank = PDM_binary_search_gap_int(i, shared_recv_idx, n_rank_in_shm+1);
      distrib_search_by_rank_n[t_rank]++;
    }

    distrib_search_by_rank_idx[0] = 0;
    for(int i = 0; i < n_rank_in_shm; ++i) {
      distrib_search_by_rank_idx[i+1] = distrib_search_by_rank_idx[i] + distrib_search_by_rank_n[i];
    }

    // PDM_log_trace_array_int(distrib_search_by_rank_idx, n_rank_in_shm+1, "distrib_search_by_rank_idx ::");

    n_pts1                = dn_search;
    pts_g_num1            = (PDM_g_num_t *) &shared_recv_gnum             [    distrib_search[i_rank_in_shm]];
    pts_coord1            = (double      *) &shared_recv_pts_coord        [3 * distrib_search[i_rank_in_shm]];
    pts_upper_bound_dist2 = (double      *) &shared_recv_upper_bound_dist2[    distrib_search[i_rank_in_shm]];

    free(distrib_search);
    free(distrib_search_by_rank_n);
    free(shared_recv_count );
    free(shared_recv_idx );


  } else {
    n_pts1 = n_pts;

    pts_g_num1            = (PDM_g_num_t *) pts_g_num;
    pts_coord1            = _pts;
    pts_upper_bound_dist2 = upper_bound_dist2;
  }

  /*
   * Solicitation local + shared
   */

  /* Management of size */
  int n_part = n_rank_in_shm;
  int          *pn_boxes           = malloc (sizeof(int         *) * n_part);
  double      **pbox_center        = malloc (sizeof(double      *) * n_part);
  double      **pbox_pts_coords    = malloc (sizeof(double      *) * n_part);
  double      **pbox_weight        = malloc (sizeof(double      *) * n_part);
  int         **pbox_init_location = malloc (sizeof(int         *) * n_part);
  int         **pstride_one        = malloc (sizeof(int         *) * n_part);
  int         **pbox_pts_n         = malloc (sizeof(int         *) * n_part);
  PDM_g_num_t **pbox_g_num         = malloc (sizeof(PDM_g_num_t *) * n_part);
  PDM_g_num_t **pbox_pts_g_num     = malloc (sizeof(PDM_g_num_t *) * n_part);


  // double t1 = PDM_MPI_Wtime();
  if(_dbbt->btShared != NULL) {
    for(int i_shm = 0; i_shm < n_rank_in_shm; ++i_shm) {

      int beg    = distrib_search_by_rank_idx[i_shm  ];
      int n_lpts = distrib_search_by_rank_idx[i_shm+1] - beg;

      PDM_g_num_t *lpts_gnum              = &pts_g_num1           [  beg];
      double      *lpts_coord             = &pts_coord1           [3*beg];
      double      *lpts_upper_bound_dist2 = &pts_upper_bound_dist2[  beg];

      // log_trace("Search in %i : n_lpts = %i \n", i_shm, n_lpts);
      if(0 == 1) {
        PDM_log_trace_array_long  (lpts_gnum,    n_lpts, "lpts_gnum ::");
        PDM_log_trace_array_double(lpts_coord, 3*n_lpts, "lpts_coord::");
      }


      // if (1) {
      //   char filename[999];
      //   sprintf(filename, "lpts_%d_%d.vtk", i_shm, i_rank);
      //   PDM_vtk_write_point_cloud(filename,
      //                             n_lpts,
      //                             lpts_coord,
      //                             lpts_gnum,
      //                             NULL);

      //   log_trace("i_shm = %d\n", i_shm);
      //   for (int i = 0; i < n_lpts; i++) {
      //     log_trace("point "PDM_FMT_G_NUM" : normalized: %f %f %f, real: %f %f %f\n",
      //               lpts_gnum[i],
      //               lpts_coord[3*i  ],
      //               lpts_coord[3*i+1],
      //               lpts_coord[3*i+2],
      //               _dbbt->s[0] + _dbbt->d[0] * lpts_coord[3*i  ],
      //               _dbbt->s[1] + _dbbt->d[1] * lpts_coord[3*i+1],
      //               _dbbt->s[2] + _dbbt->d[2] * lpts_coord[3*i+2]);
      //   }
      // }

      int *tmp_box_pts_idx = NULL;
      int *tmp_box_pts     = NULL;
      PDM_box_tree_closest_upper_bound_dist_boxes_get_shared_box_pov(_dbbt->btLoc,
                                                                     i_shm,
                                                                     n_lpts,
                                                                     lpts_coord,
                                                                     lpts_upper_bound_dist2,
                                                                     &tmp_box_pts_idx,
                                                                     &tmp_box_pts,
                                                                     _dbbt->d);

      int          n_boxes             = _dbbt->boxes->shm_boxes[i_shm].n_boxes;
      int         *boxes_init_location = _dbbt->boxes->shm_boxes[i_shm].origin;
      PDM_g_num_t *boxes_gnum          = _dbbt->boxes->shm_boxes[i_shm].g_num;
      double      *boxes_extents       = _dbbt->boxes->shm_boxes[i_shm].extents;

      // if (1) {
      //   log_trace("i_shm = %d\n", i_shm);
      //   for (int i = 0; i < n_boxes; i++) {
      //     log_trace("box "PDM_FMT_G_NUM": pts ", boxes_gnum[i]);
      //     for (int j = tmp_box_pts_idx[i]; j < tmp_box_pts_idx[i+1]; j++) {
      //       int pt_id = tmp_box_pts[j];
      //       log_trace(PDM_FMT_G_NUM" ", lpts_gnum[pt_id]);
      //     }
      //     log_trace("\n");
      //   }
      // }

      /*
       * Extract only boxes with some pts
       */
      pn_boxes[i_shm] = 0;
      for(int i_box = 0; i_box < n_boxes; ++i_box) {
        if(tmp_box_pts_idx[i_box+1] - tmp_box_pts_idx[i_box] > 0) {
          pn_boxes[i_shm]++;
        }
      }

      // pn_boxes[i_shm] = 0;
      pbox_center       [i_shm] = malloc (3 * pn_boxes[i_shm] * sizeof(double     ));
      pbox_weight       [i_shm] = malloc (    pn_boxes[i_shm] * sizeof(double     ));
      pbox_init_location[i_shm] = malloc (3 * pn_boxes[i_shm] * sizeof(int        ));
      pstride_one       [i_shm] = malloc (    pn_boxes[i_shm] * sizeof(int        ));
      pbox_pts_n        [i_shm] = malloc (    pn_boxes[i_shm] * sizeof(int        ));
      pbox_g_num        [i_shm] = malloc (    pn_boxes[i_shm] * sizeof(PDM_g_num_t));

      int n_box_pts_tot = tmp_box_pts_idx[n_boxes];
      pbox_pts_g_num    [i_shm] = malloc (    n_box_pts_tot * sizeof(PDM_g_num_t));
      pbox_pts_coords   [i_shm] = malloc (3 * n_box_pts_tot * sizeof(double     ));

      int idx_write = 0;
      pn_boxes[i_shm] = 0;
      for(int i_box = 0; i_box < n_boxes; ++i_box) {
        if(tmp_box_pts_idx[i_box+1] - tmp_box_pts_idx[i_box] == 0) {
          continue;
        }
        int i_box_e = pn_boxes[i_shm]++;
        pbox_center       [i_shm][3*i_box_e  ] = 0.5 * (boxes_extents[6*i_box  ] + boxes_extents[6*i_box+3]);
        pbox_center       [i_shm][3*i_box_e+1] = 0.5 * (boxes_extents[6*i_box+1] + boxes_extents[6*i_box+4]);
        pbox_center       [i_shm][3*i_box_e+2] = 0.5 * (boxes_extents[6*i_box+2] + boxes_extents[6*i_box+5]);

        pbox_weight       [i_shm][i_box_e] = tmp_box_pts_idx[i_box+1] - tmp_box_pts_idx[i_box];
        pbox_pts_n        [i_shm][i_box_e] = tmp_box_pts_idx[i_box+1] - tmp_box_pts_idx[i_box];
        pstride_one       [i_shm][i_box_e] = 1; // useles
        pbox_g_num        [i_shm][i_box_e] = boxes_gnum[i_box];
        pbox_init_location[i_shm][3*i_box_e  ] = boxes_init_location[3*i_box  ];
        pbox_init_location[i_shm][3*i_box_e+1] = boxes_init_location[3*i_box+1];
        pbox_init_location[i_shm][3*i_box_e+2] = boxes_init_location[3*i_box+2];

        for(int idx_pts = tmp_box_pts_idx[i_box]; idx_pts < tmp_box_pts_idx[i_box+1]; ++idx_pts) {
          int i_pts = tmp_box_pts[idx_pts];
          pbox_pts_g_num [i_shm][idx_write] = lpts_gnum[i_pts];
          pbox_pts_coords[i_shm][3*idx_write  ] = lpts_coord[3*i_pts  ];
          pbox_pts_coords[i_shm][3*idx_write+1] = lpts_coord[3*i_pts+1];
          pbox_pts_coords[i_shm][3*idx_write+2] = lpts_coord[3*i_pts+2];
          idx_write++;
        }

      }

      free(tmp_box_pts_idx);
      free(tmp_box_pts);
    }

    PDM_mpi_win_shared_unlock_all(wshared_recv_gnum);
    PDM_mpi_win_shared_unlock_all(wshared_recv_pts_coord);
    PDM_mpi_win_shared_unlock_all(wshared_recv_upper_bound_dist2);
    PDM_mpi_win_shared_free (wshared_recv_gnum);
    PDM_mpi_win_shared_free (wshared_recv_pts_coord);
    PDM_mpi_win_shared_free (wshared_recv_upper_bound_dist2);
    free(distrib_search_by_rank_idx);
  } else {

    int *tmp_box_pts_idx = NULL;
    int *tmp_box_pts     = NULL;
    PDM_box_tree_closest_upper_bound_dist_boxes_get_v2_box_pov(_dbbt->btLoc,
                                                               -1,
                                                               n_pts1,
                                                               pts_coord1,
                                                               pts_upper_bound_dist2,
                                                               &tmp_box_pts_idx,
                                                               &tmp_box_pts,
                                                               _dbbt->d);

    int          n_boxes             = _dbbt->boxes->local_boxes->n_boxes;
    int         *boxes_init_location = _dbbt->boxes->local_boxes->origin;
    PDM_g_num_t *boxes_gnum          = _dbbt->boxes->local_boxes->g_num;
    double      *boxes_extents       = _dbbt->boxes->local_boxes->extents;

    // if (1) {
    //   log_trace("local box_tree\n");
    //   for (int i = 0; i < n_boxes; i++) {
    //     log_trace("box "PDM_FMT_G_NUM": pts ", boxes_gnum[i]);
    //     for (int j = tmp_box_pts_idx[i]; j < tmp_box_pts_idx[i+1]; j++) {
    //       int pt_id = tmp_box_pts[j];
    //       log_trace(PDM_FMT_G_NUM" ", pts_g_num1[pt_id]);
    //     }
    //     log_trace("\n");
    //   }
    // }

    /*
     * Extract only boxes with some pts
     */
    pn_boxes[0] = 0;
    for(int i_box = 0; i_box < n_boxes; ++i_box) {
      if(tmp_box_pts_idx[i_box+1] - tmp_box_pts_idx[i_box] > 0) {
        pn_boxes[0]++;
      }
    }

    pbox_center       [0] = malloc (3 * pn_boxes[0] * sizeof(double     ));
    pbox_weight       [0] = malloc (    pn_boxes[0] * sizeof(double     ));
    pbox_init_location[0] = malloc (3 * pn_boxes[0] * sizeof(int        ));
    pstride_one       [0] = malloc (    pn_boxes[0] * sizeof(int        ));
    pbox_pts_n        [0] = malloc (    pn_boxes[0] * sizeof(int        ));
    pbox_g_num        [0] = malloc (    pn_boxes[0] * sizeof(PDM_g_num_t));

    int n_box_pts_tot = tmp_box_pts_idx[n_boxes];
    pbox_pts_g_num    [0] = malloc (    n_box_pts_tot * sizeof(PDM_g_num_t));
    pbox_pts_coords   [0] = malloc (3 * n_box_pts_tot * sizeof(double     ));

    int idx_write = 0;
    pn_boxes[0] = 0;
    for(int i_box = 0; i_box < n_boxes; ++i_box) {
      if(tmp_box_pts_idx[i_box+1] - tmp_box_pts_idx[i_box] == 0) {
        continue;
      }
      int i_box_e = pn_boxes[0]++;
      pbox_center       [0][3*i_box_e  ] = 0.5 * (boxes_extents[6*i_box  ] + boxes_extents[6*i_box+3]);
      pbox_center       [0][3*i_box_e+1] = 0.5 * (boxes_extents[6*i_box+1] + boxes_extents[6*i_box+4]);
      pbox_center       [0][3*i_box_e+2] = 0.5 * (boxes_extents[6*i_box+2] + boxes_extents[6*i_box+5]);

      pbox_weight       [0][i_box_e] = tmp_box_pts_idx[i_box+1] - tmp_box_pts_idx[i_box];
      pbox_pts_n        [0][i_box_e] = tmp_box_pts_idx[i_box+1] - tmp_box_pts_idx[i_box];
      pstride_one       [0][i_box_e] = 1; // useles
      pbox_g_num        [0][i_box_e] = boxes_gnum[i_box];
      pbox_init_location[0][3*i_box_e  ] = boxes_init_location[3*i_box  ];
      pbox_init_location[0][3*i_box_e+1] = boxes_init_location[3*i_box+1];
      pbox_init_location[0][3*i_box_e+2] = boxes_init_location[3*i_box+2];

      for(int idx_pts = tmp_box_pts_idx[i_box]; idx_pts < tmp_box_pts_idx[i_box+1]; ++idx_pts) {
        int i_pts = tmp_box_pts[idx_pts];
        pbox_pts_g_num[0][idx_write] = pts_g_num1[i_pts];
        pbox_pts_coords[0][3*idx_write  ] = pts_coord1[3*i_pts  ];
        pbox_pts_coords[0][3*idx_write+1] = pts_coord1[3*i_pts+1];
        pbox_pts_coords[0][3*idx_write+2] = pts_coord1[3*i_pts+2];
        idx_write++;
      }
    }

    free(tmp_box_pts_idx);
    free(tmp_box_pts);

  }
  free(_pts);

  // if (1) {
  //   for (int ipart = 0; ipart < n_part; ipart++) {
  //     int idx = 0;
  //     for (int i = 0; i < pn_boxes[ipart]; i++) {
  //       log_trace("box "PDM_FMT_G_NUM": pts ", pbox_g_num[ipart][i]);
  //       for (int j = 0; j < pbox_pts_n[ipart][i]; j++) {
  //         log_trace(PDM_FMT_G_NUM" ", pbox_pts_g_num[ipart][idx++]);
  //       }
  //       log_trace("\n");
  //     }
  //   }
  // }

  PDM_part_to_block_t *ptb = PDM_part_to_block_geom_create(PDM_PART_TO_BLOCK_DISTRIB_ALL_PROC,
                                                           PDM_PART_TO_BLOCK_POST_MERGE,
                                                           1.,
                                                           PDM_PART_GEOM_MORTON,
                                                           pbox_center,
                                                           pbox_g_num,
                                                           pbox_weight,
                                                           pn_boxes,
                                                           n_part,
                                                           _dbbt->comm);
  for(int i_part = 0; i_part < n_part; ++i_part) {
    free(pbox_center[i_part]);
    free(pbox_g_num [i_part]);
    free(pbox_weight[i_part]);
  }
  free(pbox_center);
  free(pbox_g_num );
  free(pbox_weight);

  /*
   * Exchange elmt_pts_g_num
   */
  int         *dbox_pts_n     = NULL;
  PDM_g_num_t *dbox_pts_g_num = NULL;
  int s_dbox_pts_g_num = PDM_part_to_block_exch(ptb,
                                                sizeof(PDM_g_num_t),
                                                PDM_STRIDE_VAR_INTERLACED,
                                                1,
                                                pbox_pts_n,
                                      (void **) pbox_pts_g_num,
                                                &dbox_pts_n,
                                      (void **) &dbox_pts_g_num);
  PDM_UNUSED(s_dbox_pts_g_num);

  for(int i_part = 0; i_part < n_part; ++i_part) {
    free(pbox_pts_g_num[i_part]);
  }
  free(pbox_pts_g_num);
  free(dbox_pts_n);
  dbox_pts_n = NULL;

  double *dbox_pts_coord = NULL;
  int n_recv_data = PDM_part_to_block_exch(ptb,
                                           3 * sizeof(double),
                                           PDM_STRIDE_VAR_INTERLACED,
                                           1,
                                           pbox_pts_n,
                                 (void **) pbox_pts_coords,
                                           &dbox_pts_n,
                                 (void **) &dbox_pts_coord);

  for(int i_part = 0; i_part < n_part; ++i_part) {
    free(pbox_pts_n     [i_part]);
    free(pbox_pts_coords[i_part]);
  }
  free(pbox_pts_n     );
  free(pbox_pts_coords);

  /*
   * Exchange elmt_init_location and compress
   */
  int *dbox_init_location_n = NULL;
  int *dbox_init_location   = NULL;
  PDM_part_to_block_exch(ptb,
                         3 * sizeof(int),
                         PDM_STRIDE_VAR_INTERLACED,
                         1,
                         pstride_one,
               (void **) pbox_init_location,
                         &dbox_init_location_n,
               (void **) &dbox_init_location);
  for(int i_part = 0; i_part < n_part; ++i_part) {
    free(pbox_init_location[i_part]);
    free(pstride_one[i_part]);
  }
  free(pbox_init_location);
  free(pstride_one);

  /* Merge */
  int n_equi_box = PDM_part_to_block_n_elt_block_get(ptb);
  PDM_g_num_t *ptb_equi_box_g_num = PDM_part_to_block_block_gnum_get(ptb);

  int max_box_pts_n = 0;
  for (int i = 0; i < n_equi_box; i++) {
    max_box_pts_n = PDM_MAX(max_box_pts_n, dbox_pts_n[i]);
  }

  int    *pts_unique_order = malloc(    max_box_pts_n * sizeof(int   ));
  double *tmp_coord        = malloc(3 * max_box_pts_n * sizeof(double));

  int idx_read  = 0;
  int idx_write = 0;
  for (int i = 0; i < n_equi_box; i++) {
    PDM_g_num_t *_g_num = dbox_pts_g_num + idx_read;
    double      *_coord = dbox_pts_coord + idx_read*3;

    // if (1) {
    //   for (int j = 0; j < dbox_pts_n[i]; j++) {
    //     log_trace("blk point "PDM_FMT_G_NUM" : %f %f %f\n",
    //               _g_num[j],
    //               _dbbt->s[0] + _dbbt->d[0] * _coord[3*j  ],
    //               _dbbt->s[1] + _dbbt->d[1] * _coord[3*j+1],
    //               _dbbt->s[2] + _dbbt->d[2] * _coord[3*j+2]);
    //   }
    // }

    for (int j = 0; j < dbox_pts_n[i]; j++) {
      pts_unique_order[j] = j;
    }

    // PDM_log_trace_array_long(_g_num, dbox_pts_n[i], "_g_num before unique : ");
    int n_unique_pts = PDM_inplace_unique_long2(_g_num,
                                                pts_unique_order,
                                                0,
                                                dbox_pts_n[i]-1);
    // PDM_log_trace_array_int(pts_unique_order, n_unique_pts, "pts_unique_order : ");
    // log_trace("box "PDM_FMT_G_NUM" n_unique_pts = %d / %d, idx_write = %d->%d / %d\n",
    //           ptb_equi_box_g_num[i],
    //           n_unique_pts, dbox_pts_n[i],
    //           idx_write, idx_write + n_unique_pts, s_dbox_pts_g_num);
    // memcpy(dbox_pts_g_num + idx_write, _g_num, sizeof(PDM_g_num_t) * n_unique_pts);
    for (int j = 0; j < n_unique_pts; j++) {
      // log_trace("  j = %d -> %d\n", j, idx_write+j);
      dbox_pts_g_num[idx_write+j] = _g_num[j];
    }

    memcpy(tmp_coord, _coord, sizeof(double) * dbox_pts_n[i]*3);
    for (int j = 0; j < dbox_pts_n[i]; j++) {
      int idx = idx_write + pts_unique_order[j];
      memcpy(dbox_pts_coord + idx*3, tmp_coord + j*3, sizeof(double)*3);
    }
    // for (int j = 0; j < n_unique_pts; j++) {
    //   int idx = pts_unique_order[j];
    //   log_trace("** "PDM_FMT_G_NUM" -> %f %f %f\n",
    //             _g_num[j],
    //             _dbbt->s[0] + _dbbt->d[0] * tmp_coord[3*idx],
    //             _dbbt->s[1] + _dbbt->d[1] * tmp_coord[3*idx+1],
    //             _dbbt->s[2] + _dbbt->d[2] * tmp_coord[3*idx+2]);
    //   memcpy(dbox_pts_coord + idx_write*3, tmp_coord + idx*3, sizeof(double)*3);
    //   // for (int k = 0; k < 3; k++) {
    //   //   dbox_pts_coord[3*idx_write+k] = _coord[3*idx+k];
    //   // }
    //   idx_write++;
    // }
    idx_write += n_unique_pts;
    idx_read  += dbox_pts_n[i];

    dbox_pts_n[i] = n_unique_pts;
  }
  free(pts_unique_order);
  free(tmp_coord);
  dbox_pts_g_num = realloc(dbox_pts_g_num, sizeof(PDM_g_num_t) * idx_write);
  dbox_pts_coord = realloc(dbox_pts_coord, sizeof(double     ) * idx_write*3);

  n_recv_data = idx_write;

  idx_read = 0;
  for (int i = 0; i < n_equi_box; i++) {
    // memcpy(dbox_init_location + 3*i, dbox_init_location + 3*idx_read, sizeof(int)*3);
    for (int j = 0; j < 3; j++) {
      dbox_init_location[3*i+j] = dbox_init_location[3*idx_read+j];
    }
    idx_read += dbox_init_location_n[i];
  }
  free(dbox_init_location_n);
  dbox_init_location = realloc(dbox_init_location, sizeof(int) * n_equi_box*3);



  PDM_g_num_t *equi_box_gnum = malloc(n_equi_box * sizeof(PDM_g_num_t));

  for(int i_box = 0; i_box < n_equi_box; ++i_box) {
    equi_box_gnum[i_box] = ptb_equi_box_g_num[i_box];
  }
  free(pn_boxes);

  PDM_part_to_block_free(ptb);


  //-->>
  for (int i = 0; i < n_recv_data; i++) {
    for (int j = 0; j < 3; j++) {
      dbox_pts_coord[3*i+j] = _dbbt->s[j] + _dbbt->d[j] * (dbox_pts_coord[3*i+j]);
    }
    // if (dbox_pts_g_num[i] == 9) {
    //   log_trace("dbox_pts_coord ("PDM_FMT_G_NUM") = %f %f %f\n",
    //             dbox_pts_g_num[i],
    //             dbox_pts_coord[3*i+0],
    //             dbox_pts_coord[3*i+1],
    //             dbox_pts_coord[3*i+2]);
    // }
  }
  //<<--

  PDM_MPI_Type_free(&mpi_pts_coords_type);

  PDM_MPI_Comm_free(&comm_shared);

  *out_n_extract_boxes   = n_equi_box;
  *out_box_gnum          = equi_box_gnum;
  *out_box_init_location = dbox_init_location;
  *out_dbox_pts_idx      = PDM_array_new_idx_from_sizes_int(dbox_pts_n, n_equi_box);
  *out_dbox_pts_g_num    = dbox_pts_g_num;
  *out_dbox_pts_coord    = dbox_pts_coord;
  free(dbox_pts_n);
}


/**
 *
 * \brief Get the boxes closer than the upper bound distance (Asynchronous)
 *
 *   \param [in]  bt                 Pointer to box tree structure
 *   \param [in]  n_pts              Number of points
 *   \param [in]  pts                Point coordinates (size = 3 * \ref n_pts)
 *   \param [in]  pts_g_num          Point global ids
 *   \param [in]  upper_bound_dist2  Upper bound of the square of the distance (size = \ref n_pts)
 *   \param [out] box_index          Index of boxes (size = \ref n_pts + 1)
 *   \param [out] box_g_num          Global ids of boxes (size = \ref i_boxes[\ref n_pts])
 *
 */

void
PDM_dbbtree_closest_upper_bound_dist_boxes_get_async
(
 PDM_dbbtree_t   *dbbt,
 const int        n_pts,
 double           pts[],
 PDM_g_num_t      pts_g_num[],
 double           upper_bound_dist2[],
 int             *box_index[],
 PDM_g_num_t     *box_g_num[]
 )
{
  PDM_UNUSED(pts_g_num);
  /*
   * RANK DATA COPY PARAMETERS
   */
  const double RANK_COPY_threshold  = 1.2;  // factor of the mean nb of requests
  const double RANK_COPY_max_copies = 0.15; // factor of the total nb of processes

  /*
   * Initialization
   */
  assert (dbbt != NULL);
  _PDM_dbbtree_t *_dbbt = (_PDM_dbbtree_t *) dbbt;

  double *_pts = pts;
  if (_dbbt->global_extents != NULL) {
    _pts = malloc (sizeof(double) * n_pts * 3);

    for (int i = 0; i < n_pts; i++) {
      _normalize (_dbbt,
                  pts + 3*i,
                  _pts + 3*i);
    }
  }

  int     n_pts_local            = n_pts;
  double *pts_local              = _pts;
  double *upper_bound_dist_local = upper_bound_dist2;


  int myRank;
  PDM_MPI_Comm_rank (_dbbt->comm, &myRank);
  int lComm;
  PDM_MPI_Comm_size (_dbbt->comm, &lComm);

  /*
   * Determine for each point the list of involved processes
   */
  int *n_send_pts = NULL;
  int *n_recv_pts = NULL;

  int *box_index_tmp = NULL;
  int *box_l_num_tmp = NULL;

  int *n_pts_rank = NULL;
  int *n_pts_send = NULL;
  int *n_pts_recv = NULL;
  int n_pts_recv_total = 0;

  int *i_pts_rank = NULL;
  int *i_pts_send = NULL;
  int *i_pts_recv = NULL;

  double *pts_rank              = NULL;
  double *upper_bound_dist_rank = NULL;
  double *pts_recv              = NULL;
  double *upper_bound_dist_recv = NULL;

  int n_copied_ranks = 0;
  int *copied_ranks  = NULL;
  int *rank_copy_num = NULL;

  const int *usedRanks = _dbbt->usedRank;

  const int idebug = 0;

  int i1 = 0, i2 = 0, i3 = 0;
  int i_rank = 0;
  PDM_MPI_Request request;
  double *data_send = NULL;
  double *data_recv = NULL;
  if (_dbbt->btShared != NULL) {

    if (idebug) {
      printf ("  **** deb PDM_box_tree_closest_upper_bound_dist_boxes_get shared _pts : %d\n", n_pts);
    }

    PDM_box_tree_closest_upper_bound_dist_boxes_get (_dbbt->btShared,
                                                     n_pts,
                                                     pts,
                                                     upper_bound_dist2,
                                                     &box_index_tmp,
                                                     &box_l_num_tmp);

    if (idebug) {
      printf ("  **** fin PDM_box_tree_closest_upper_bound_dist_boxes_get shared n_pts : %d\n", n_pts);
      for (int i = 0; i < n_pts; i++) {
        printf ("%d : (%12.5e %12.5e %12.5e) %12.5e\n", i,
                pts[3*i], pts[3*i+1], pts[3*i+2],
                upper_bound_dist2[i]);
        printf ("  boxes %d :" , box_index_tmp[i+1] - box_index_tmp[i]);
        for (int j = box_index_tmp[i]; j < box_index_tmp[i+1]; j++) {
          printf (" %d", box_l_num_tmp[j]);
        }
        printf ("\n");
      }
    }

    /*
     * Count (provisional) nb of points to send to each process
     */
    n_send_pts = PDM_array_zeros_int(lComm);
    n_recv_pts = malloc (sizeof(int) * lComm);


    for (int i = 0; i < n_pts; i++) {
      for (int j = box_index_tmp[i]; j < box_index_tmp[i+1]; j++) {
        n_send_pts[usedRanks[box_l_num_tmp[j]]]++;
      }
    }

    PDM_MPI_Alltoall (n_send_pts, 1, PDM_MPI_INT,
                      n_recv_pts, 1, PDM_MPI_INT,
                      _dbbt->comm);
    free(n_send_pts);

    /*
     * Prepare copies
     */
    // total nb of requests received by current process
    int local_sum_nrecv = 0;
    for (int i = 0; i < lComm; i++) {
      local_sum_nrecv += n_recv_pts[i];
    }
    free(n_recv_pts);

    int *n_requests = malloc (lComm * sizeof(int));
    PDM_MPI_Allgather (&local_sum_nrecv, 1, PDM_MPI_INT,
                       n_requests,       1, PDM_MPI_INT,
                       _dbbt->comm);

    // mean nb of requests
    int mean_n_requests = 0;
    for (int i = 0; i < lComm; i++) {
      mean_n_requests += n_requests[i];
    }
    mean_n_requests /= lComm;

    /* sort the ranks in ascending order of
     * the total nb of points they are supposed to receive */
    int *order = malloc (lComm * sizeof(int));
    for (int i = 0; i < lComm; i ++) {
      order[i] = i;
    }

    PDM_sort_int (n_requests, order, lComm);

    // identify ranks to be copied
    double threshold_n_req = RANK_COPY_threshold*mean_n_requests;
    int max_copied_ranks   = (int) _MAX (1, RANK_COPY_max_copies*lComm);

    n_copied_ranks = 0;
    copied_ranks = malloc (max_copied_ranks * sizeof(int));

    for (int i = 0; i < max_copied_ranks; i++) {
      i_rank = lComm - 1 - i;

      if ( n_requests[i_rank] > threshold_n_req ) {
        copied_ranks[n_copied_ranks++] = order[i_rank];
      } else {
        break;
      }
    }
    free(order);
    free(n_requests);

    //------------->>>
    if ( idebug && myRank == 0 ) {
      if ( n_copied_ranks == 0 ) {
        printf("n_copied_ranks = 0\n");
      } else {
        printf("copied rank(s) = ");
        for (int i = 0; i < n_copied_ranks; i++) {
          printf("%d ", copied_ranks[i]);
        }
        printf("\n");
      }
    }
    //<<<-------------

    /*
     * Copy the data of selected ranks
     */
    rank_copy_num = (int *) malloc (sizeof(int) * lComm);
    PDM_box_tree_copy_to_ranks (_dbbt->btLoc,
                                &n_copied_ranks,
                                copied_ranks,
                                rank_copy_num);
    /* rank_copy_num[_dbbt->btLoc->copied_ranks[i]] (def)= i*/
    free(copied_ranks);


    /*
     * Distribution of points...
     *    ..._local --> search in local box tree
     *    ..._rank  --> search in copied box trees
     *    ..._send  --> search in distant box trees (send to other processes)
     */
    n_pts_local = 0;

    n_pts_rank = PDM_array_zeros_int(n_copied_ranks);

    n_pts_send = PDM_array_zeros_int(lComm);
    n_pts_recv = malloc (sizeof(int) * lComm);


    for (int i = 0; i < n_pts; i++) {
      for (int j = box_index_tmp[i]; j < box_index_tmp[i+1]; j++) {
        i_rank = usedRanks[box_l_num_tmp[j]];
        if ( i_rank == myRank ) {
          // ---> search in btLoc->local_data of current process
          n_pts_local++;
        } else if ( rank_copy_num[i_rank] >= 0 ) {
          // ---> search in btLoc->rank_data[rank_copy_num[i_rank]] of current process
          n_pts_rank[rank_copy_num[i_rank]]++;
        } else {
          // ---> search in btLoc->local_data of process with rank i_rank
          n_pts_send[i_rank]++;
        }
      }
    }

    PDM_MPI_Alltoall (n_pts_send, 1, PDM_MPI_INT,
                      n_pts_recv, 1, PDM_MPI_INT,
                      _dbbt->comm);

    i_pts_rank = PDM_array_new_idx_from_sizes_int(n_pts_rank, n_copied_ranks);

    i_pts_send = malloc (sizeof(int) * (lComm+1));
    i_pts_recv = malloc (sizeof(int) * (lComm+1));
    i_pts_send[0] = 0;
    i_pts_recv[0] = 0;
    for (int i = 0; i < lComm; i++) {
      i_pts_send[i+1] = i_pts_send[i] + 4 * n_pts_send[i];
      i_pts_recv[i+1] = i_pts_recv[i] + 4 * n_pts_recv[i];
      n_pts_recv[i] *= 4;
    }

    pts_local              = malloc (sizeof(double) * n_pts_local*3);
    upper_bound_dist_local = malloc (sizeof(double) * n_pts_local);

    pts_rank              = malloc (sizeof(double) * i_pts_rank[n_copied_ranks]*3);
    upper_bound_dist_rank = malloc (sizeof(double) * i_pts_rank[n_copied_ranks]);

    data_send = malloc (sizeof(double) * i_pts_send[lComm]);
    data_recv = malloc (sizeof(double) * i_pts_recv[lComm]);

    PDM_array_reset_int(n_pts_send, lComm, 0);
    PDM_array_reset_int(n_pts_rank, n_copied_ranks, 0);

    i1 = 0; i2 = 0; i3 = 0;
    for (int i = 0; i < n_pts; i++) {
      for (int j = box_index_tmp[i]; j < box_index_tmp[i+1]; j++) {
        i_rank = usedRanks[box_l_num_tmp[j]];
        if ( i_rank == myRank ) {
          // pts_local, upper_bound_dist_local (local points, local box tree)
          upper_bound_dist_local[i1] = upper_bound_dist2[i];
          pts_local[3*i1]            = _pts[3*i];
          pts_local[3*i1+1]          = _pts[3*i+1];
          pts_local[3*i1+2]          = _pts[3*i+2];
          i1++;
        } else if ( rank_copy_num[i_rank] >= 0 ) {
          // pts_rank, upper_bound_dist_rank (local points, distant (copied) box trees)
          int j_rank = rank_copy_num[i_rank];
          i2 = i_pts_rank[j_rank] + n_pts_rank[j_rank];
          upper_bound_dist_rank[i2] = upper_bound_dist2[i];
          pts_rank[3*i2]            = _pts[3*i];
          pts_rank[3*i2+1]          = _pts[3*i+1];
          pts_rank[3*i2+2]          = _pts[3*i+2];
          n_pts_rank[j_rank]++;
        } else {
          // data_send (local points, distant (not copied) box trees)
          i3 = i_pts_send[i_rank] + 4*n_pts_send[i_rank];
          data_send[i3++] = _pts[3*i];
          data_send[i3++] = _pts[3*i+1];
          data_send[i3++] = _pts[3*i+2];
          data_send[i3++] = upper_bound_dist2[i];

          n_pts_send[i_rank]++;
        }
      }
    }

    for (int i = 0; i < lComm; i++) {
      n_pts_send[i] *= 4;
    }


    // Send points to search in distant (not copied) box trees
    // PDM_MPI_Alltoallv (data_send, n_pts_send, i_pts_send, PDM_MPI_DOUBLE,
    //                    data_recv, n_pts_recv, i_pts_recv, PDM_MPI_DOUBLE,
    //                    _dbbt->comm);
    // free(data_send);

    // n_pts_recv_total = i_pts_recv[lComm] / 4;

    // pts_recv              = (double *) malloc (sizeof(double) * n_pts_recv_total * 3);
    // upper_bound_dist_recv = (double *) malloc (sizeof(double) * n_pts_recv_total);

    // for (int i = 0; i < n_pts_recv_total; i++) {
    //   for (int j = 0; j < 3; j++) {
    //     pts_recv[3*i+j] = data_recv[4*i+j];
    //   }
    //   upper_bound_dist_recv[i] = data_recv[4*i+3];
    // }
    // free(data_recv);

    /* Asynchrone */
    PDM_MPI_Ialltoallv (data_send, n_pts_send, i_pts_send, PDM_MPI_DOUBLE,
                        data_recv, n_pts_recv, i_pts_recv, PDM_MPI_DOUBLE,
                        _dbbt->comm, &request);



  }
  if (_pts != pts && _pts != pts_local) {
    free (_pts);
  }


  // Determine candidate boxes in local box tree (local points)
  int *box_index_local;
  int *box_l_num_local;

  // Determine candidate boxes in local box tree (received points)
  int *n_pts_send2 = NULL;
  int *i_pts_send2 = NULL;
  int *n_box_l_num_per_pts = NULL;
  PDM_g_num_t *box_g_num_per_pts = NULL;
  int *box_index_recv = NULL;
  int *box_l_num_recv = NULL;
  int *n_box_l_num_recv = NULL;

  if (_dbbt->btShared != NULL) {
    PDM_MPI_Wait(&request);
    free(data_send);

    n_pts_recv_total = i_pts_recv[lComm] / 4;

    pts_recv              = (double *) malloc (sizeof(double) * n_pts_recv_total * 3);
    upper_bound_dist_recv = (double *) malloc (sizeof(double) * n_pts_recv_total);

    for (int i = 0; i < n_pts_recv_total; i++) {
      for (int j = 0; j < 3; j++) {
        pts_recv[3*i+j] = data_recv[4*i+j];
      }
      upper_bound_dist_recv[i] = data_recv[4*i+3];
    }
    free(data_recv);

    PDM_box_tree_closest_upper_bound_dist_boxes_get_v2 (_dbbt->btLoc,
                                                        -1, // search in local box tree
                                                        n_pts_recv_total,
                                                        pts_recv,
                                                        upper_bound_dist_recv,
                                                        &box_index_recv,
                                                        &box_l_num_recv,
                                                        _dbbt->d);
    free(pts_recv);
    free(upper_bound_dist_recv);

    /*
     * Send back results for distant points to original processes:
     *     - nb of boxes for each point
     *     - global numbering of these boxes
     */

    n_box_l_num_recv = malloc (sizeof(int) * n_pts_recv_total);

    // log_debug("n_pts_recv_total %i \n", n_pts_recv_total);
    for (int i = 0; i < n_pts_recv_total; i++) {
      n_box_l_num_recv[i] = box_index_recv[i+1] - box_index_recv[i];
    }

    for (int i = 0; i < lComm; i++) {
      i_pts_send[i+1] = i_pts_send[i+1]/4;
      i_pts_recv[i+1] = i_pts_recv[i+1]/4;
      n_pts_send[i]   = n_pts_send[i]/4;
      n_pts_recv[i]   = n_pts_recv[i]/4;
    }

    n_box_l_num_per_pts = malloc (sizeof(int) * i_pts_send[lComm]);

    // double t1i = PDM_MPI_Wtime();
    // PDM_MPI_Alltoallv (n_box_l_num_recv,    n_pts_recv, i_pts_recv, PDM_MPI_INT,
    //                    n_box_l_num_per_pts, n_pts_send, i_pts_send, PDM_MPI_INT,
    //                    _dbbt->comm);
    PDM_MPI_Ialltoallv (n_box_l_num_recv,    n_pts_recv, i_pts_recv, PDM_MPI_INT,
                        n_box_l_num_per_pts, n_pts_send, i_pts_send, PDM_MPI_INT,
                        _dbbt->comm, &request);
    // double dti = PDM_MPI_Wtime() - t1i;
    // log_debug("dt first all_to_all %12.5e \n", dti);
    // PDM_log_trace_array_int(n_pts_recv, lComm, "n_pts_recv :: ");
    // PDM_log_trace_array_int(n_pts_send, lComm, "n_pts_send :: ");
  }

  // double t1 = PDM_MPI_Wtime();
  // log_debug("n_pts_local %i \n", n_pts_local);
  PDM_box_tree_closest_upper_bound_dist_boxes_get_v2 (_dbbt->btLoc,
                                                      -1, // search in local box tree
                                                      n_pts_local,
                                                      pts_local,
                                                      upper_bound_dist_local,
                                                      &box_index_local,
                                                      &box_l_num_local,
                                                      _dbbt->d);
  if (pts_local != pts) {
    free(pts_local);
  }
  if (upper_bound_dist_local != upper_bound_dist2) {
    free(upper_bound_dist_local);
  }

  // conversion local --> global numbering
  const PDM_g_num_t *gnum_boxes_local = PDM_box_set_get_g_num (_dbbt->boxes);
  PDM_g_num_t *box_g_num_local = malloc(sizeof(PDM_g_num_t) * box_index_local[n_pts_local]);

  for (int i = 0; i < box_index_local[n_pts_local]; i++) {
    box_g_num_local[i] = gnum_boxes_local[box_l_num_local[i]];
  }
  free(box_l_num_local);

  if (_dbbt->btShared == NULL) {

    *box_index = box_index_local;
    *box_g_num = box_g_num_local;

  } else {
    // Determine candidate boxes in copied box trees (local points)
    int **box_index_rank;
    int **box_l_num_rank;

    box_index_rank = (int **) malloc (sizeof(int *) * n_copied_ranks);
    box_l_num_rank = (int **) malloc (sizeof(int *) * n_copied_ranks);

    for (int i_copied_rank = 0; i_copied_rank < n_copied_ranks; i_copied_rank++) {
      int n_pts_copied_rank = i_pts_rank[i_copied_rank+1] - i_pts_rank[i_copied_rank];
      double *pts_copied_rank = pts_rank + 3*i_pts_rank[i_copied_rank];
      double *upper_bound_dist_copied_rank = upper_bound_dist_rank + i_pts_rank[i_copied_rank];
      PDM_box_tree_closest_upper_bound_dist_boxes_get_v2 (_dbbt->btLoc,
                                                          i_copied_rank,
                                                          n_pts_copied_rank,
                                                          pts_copied_rank,
                                                          upper_bound_dist_copied_rank,
                                                          &(box_index_rank[i_copied_rank]),
                                                          &(box_l_num_rank[i_copied_rank]),
                                                          _dbbt->d);
    }

    free(pts_rank);
    free(upper_bound_dist_rank);

    // conversion local --> global numbering for each copied rank
    PDM_g_num_t **box_g_num_rank = malloc(sizeof(PDM_g_num_t *) * n_copied_ranks);

    PDM_g_num_t *gnum_boxes_rank = NULL;
    for (int i_copied_rank = 0; i_copied_rank < n_copied_ranks; i_copied_rank++) {
      gnum_boxes_rank = PDM_box_set_get_rank_boxes_g_num (_dbbt->boxes,
                                                          i_copied_rank);
      box_g_num_rank[i_copied_rank] = malloc(sizeof(PDM_g_num_t) * box_index_rank[i_copied_rank][n_pts_rank[i_copied_rank]]);
      for (int i = 0; i < box_index_rank[i_copied_rank][n_pts_rank[i_copied_rank]]; i++) {
        box_g_num_rank[i_copied_rank][i] = gnum_boxes_rank[box_l_num_rank[i_copied_rank][i]];
      }

      free(box_l_num_rank[i_copied_rank]);
    }
    free(box_l_num_rank);
    // double dt = PDM_MPI_Wtime() - t1;
    // log_debug("dt first management copie %12.5e \n", dt);

    if (_dbbt->btShared != NULL) {
      PDM_MPI_Wait(&request);
    }



    n_pts_send2 = malloc (sizeof(int) * lComm);
    i_pts_send2 = malloc (sizeof(int) * (lComm+1));

    int *n_pts_recv2 = malloc (sizeof(int) * lComm);
    int *i_pts_recv2 = malloc (sizeof(int) * (lComm+1));

    for (int i = 0; i < lComm; i++) {
      n_pts_send2[i] = 0;
      n_pts_recv2[i] = 0;
    }

    for (int i = 0; i < lComm; i++) {
      for (int j = i_pts_recv[i]; j < i_pts_recv[i+1]; j++) {
        n_pts_recv2[i] += n_box_l_num_recv[j];
      }
      for (int j = i_pts_send[i]; j < i_pts_send[i+1]; j++) {
        n_pts_send2[i] += n_box_l_num_per_pts[j];
      }
    }

    free(n_box_l_num_recv);

    i_pts_send2[0] = 0;
    i_pts_recv2[0] = 0;
    for (int i = 0; i < lComm; i++) {
      i_pts_send2[i+1] = i_pts_send2[i] + n_pts_send2[i];
      i_pts_recv2[i+1] = i_pts_recv2[i] + n_pts_recv2[i];
    }


    // Conversion local --> global numbering
    PDM_g_num_t *box_g_num_recv = malloc(sizeof(PDM_g_num_t) * box_index_recv[n_pts_recv_total]);

    for (int i = 0; i < box_index_recv[n_pts_recv_total]; i++) {
      box_g_num_recv[i] = gnum_boxes_local[box_l_num_recv[i]];
    }

    box_g_num_per_pts = malloc(sizeof(PDM_g_num_t) * i_pts_send2[lComm]);
    // t1 = PDM_MPI_Wtime();
    PDM_MPI_Alltoallv (box_g_num_recv,    n_pts_recv2, i_pts_recv2, PDM__PDM_MPI_G_NUM,
                       box_g_num_per_pts, n_pts_send2, i_pts_send2, PDM__PDM_MPI_G_NUM,
                       _dbbt->comm);
    // dt = PDM_MPI_Wtime() - t1;
    // log_debug("dt second all_to_all %12.5e \n", dt);
    // PDM_log_trace_array_int(n_pts_recv2, lComm, "n_pts_recv2 :: ");
    // PDM_log_trace_array_int(n_pts_send2, lComm, "n_pts_send2 :: ");

    free(box_index_recv);
    free(box_l_num_recv);
    free(box_g_num_recv);

    free(n_pts_recv);
    free(i_pts_recv);
    free(n_pts_recv2);
    free(i_pts_recv2);




    /*
     * Merge all results and resolve duplicates
     */
    *box_index = malloc(sizeof(int) * (n_pts + 1));

    int max_n_box_g_num = box_index_local[n_pts_local];
    for (int i_copied_rank = 0; i_copied_rank < n_copied_ranks; i_copied_rank++) {
      max_n_box_g_num += box_index_rank[i_copied_rank][n_pts_rank[i_copied_rank]];
    }
    max_n_box_g_num += i_pts_send2[lComm];

    *box_g_num = malloc(sizeof(PDM_g_num_t) * max_n_box_g_num);


    int *rank_index = (int *) malloc(sizeof(int) * (n_pts+1));
    memcpy(rank_index, box_index_tmp, sizeof(int) * (n_pts+1));

    for (int i = 0; i < lComm; i++) {
      n_pts_send[i] = 0;
      n_pts_send2[i] = 0;
    }

    PDM_array_reset_int(n_pts_rank, n_copied_ranks, 0);



    int keyMax = 3 * n_pts;
    int key = 0;
    int found = 0;
    PDM_hash_tab_t *ht = PDM_hash_tab_create (PDM_HASH_TAB_KEY_INT,
                                              &keyMax);

    PDM_g_num_t i_box = 0;

    box_index_tmp[0] = 0;
    int idx = 0;
    i1 = 0; i2 = 0; i3 = 0;
    if ( 1 ) {
      for (int i = 0; i < n_pts; i++) { // loop over local points
        box_index_tmp[i+1] = box_index_tmp[i];
        for (int j = rank_index[i]; j < rank_index[i+1]; j++) { // loop over procs to which the current point was sent
          i_rank = usedRanks[box_l_num_tmp[j]]; // i_rank = rank j-th proc to which the current point was sent

          if ( i_rank == myRank ) { // boxes local to current process
            for (int k = box_index_local[i1]; k < box_index_local[i1+1]; k++) {
              i_box = box_g_num_local[k];

              found = PDM_hash_tab_check_collision (ht, i_box, keyMax, &key);

              if (!found) {
                PDM_hash_tab_data_add (ht, (void *) &key, *box_g_num + idx);
                (*box_g_num)[idx++] = i_box;
                box_index_tmp[i+1] += 1;
              }

            }
            i1++;

          } else if ( rank_copy_num[i_rank] >= 0 ) { // distant boxes copied in current process
            int j_rank = rank_copy_num[i_rank];
            i2 = n_pts_rank[j_rank];

            for (int k = box_index_rank[j_rank][i2]; k < box_index_rank[j_rank][i2+1]; k++) {
              i_box = box_g_num_rank[j_rank][k];

              found = PDM_hash_tab_check_collision (ht, i_box, keyMax, &key);

              if (!found) {
                PDM_hash_tab_data_add (ht, (void *) &key, *box_g_num + idx);
                (*box_g_num)[idx++] = i_box;
                box_index_tmp[i+1] += 1;
              }

            }
            n_pts_rank[j_rank]++;

          } else { // distant boxes (not copied)
            i3 = n_pts_send[i_rank];
            int i4 = i_pts_send2[i_rank] + n_pts_send2[i_rank];
            int i5 = i_pts_send[i_rank] + i3;
            for (int k = 0; k < n_box_l_num_per_pts[i5]; k++) {
              i_box = box_g_num_per_pts[i4++];

              found = PDM_hash_tab_check_collision (ht, i_box, keyMax, &key);

              if (!found) {
                PDM_hash_tab_data_add (ht, (void *) &key, *box_g_num + idx);
                (*box_g_num)[idx++] = i_box;
                box_index_tmp[i+1] += 1;
              }

            }
            n_pts_send2[i_rank] += n_box_l_num_per_pts[i5];
            n_pts_send[i_rank]++;

          }
        }
        PDM_hash_tab_purge (ht, PDM_FALSE);
      }
    }
    PDM_hash_tab_free (ht);

    free (n_box_l_num_per_pts);


    if ( *box_index != NULL ) {
      free (*box_index);
    }
    *box_index = box_index_tmp;


    *box_g_num = realloc (*box_g_num, sizeof(PDM_g_num_t) * box_index_tmp[n_pts]);

    /*
     * Deallocate stuff
     */
    free(box_l_num_tmp);
    free(box_g_num_per_pts);

    free(box_index_local);
    free(box_g_num_local);

    for (int i_copied_rank = 0; i_copied_rank < n_copied_ranks; i_copied_rank++) {
      free(box_index_rank[i_copied_rank]);
      free(box_g_num_rank[i_copied_rank]);
    }
    free(box_index_rank);
    free(box_g_num_rank);


    free(n_pts_send);
    free(i_pts_send);
    free(n_pts_send2);
    free(i_pts_send2);

    free(i_pts_rank);
    free(n_pts_rank);

    free(rank_index);
    free(rank_copy_num);
  }

  PDM_box_tree_free_copies(_dbbt->btLoc);
}



/**
 *
 * \brief Export boxes and points to ASCII VTK format for debugging purposes
 *
 * A set of files (normalized/unnormalized points and boxes) is written
 * for each MPI rank.
 *
 *   \param [in]  filename_pattern  File name pattern (prefix)
 *   \param [in]  dbbt              Pointer to distributed box tree structure
 *   \param [in]  n_pts             Number of points
 *   \param [in]  pts_g_num         Point global ids (size = \ref n_pts)
 *   \param [in]  pts_coord         Point coordinates (size = 3 * \ref n_pts)
 *
 */

void
PDM_dbbtree_points_debug
(
 char*               filename_pattern,
 PDM_dbbtree_t      *dbbt,
 const int           n_pts,
 double             *pts_coord,
 PDM_g_num_t        *pts_g_num
)
{
  assert (dbbt != NULL);
  _PDM_dbbtree_t *_dbbt = (_PDM_dbbtree_t *) dbbt;

  int my_rank;
  PDM_MPI_Comm_rank (_dbbt->comm, &my_rank);
  char filename[999];
  //-->>
  double *_pts_coord = malloc (sizeof(double) * n_pts * 3);
  for (int i = 0; i < n_pts; i++) {
    _normalize (_dbbt,
                pts_coord + 3*i,
                _pts_coord + 3*i);
  }

  sprintf(filename, "%s_dbbt_pts_n_%3.3d.vtk", filename_pattern, my_rank);
  _export_point_cloud (filename,
                       1,
                       &n_pts,
                       &_pts_coord,
                       &pts_g_num,
                       NULL);

  sprintf(filename, "%s_dbbt_boxes_n_%3.3d.vtk", filename_pattern, my_rank);
  PDM_vtk_write_boxes (filename,
                       _dbbt->boxes->local_boxes->n_boxes,
                       _dbbt->boxes->local_boxes->extents,
                       _dbbt->boxes->local_boxes->g_num);


  double *_extents = malloc (sizeof(double) * _dbbt->boxes->local_boxes->n_boxes * 6);
  for (int i = 0; i < 2*_dbbt->boxes->local_boxes->n_boxes; i++) {
    for (int j = 0; j < 3; j++) {
      _extents[3*i+j] = _dbbt->s[j] + _dbbt->d[j] * _dbbt->boxes->local_boxes->extents[3*i+j];
    }
  }

  sprintf(filename, "%s_dbbt_pts_%3.3d.vtk", filename_pattern, my_rank);
  _export_point_cloud (filename,
                       1,
                       &n_pts,
                       &pts_coord,
                       &pts_g_num,
                       NULL);

  sprintf(filename, "%s_dbbt_boxes_%3.3d.vtk", filename_pattern, my_rank);
  PDM_vtk_write_boxes (filename,
                       _dbbt->boxes->local_boxes->n_boxes,
                       _extents,
                       _dbbt->boxes->local_boxes->g_num);

  free(_extents);
  free(_pts_coord);
}



/**
 *
 * \brief Get an indexed list of all points inside the boxes of a distributed box tree
 *
 *   \param [in]  dbbt               Pointer to distributed box tree structure
 *   \param [in]  n_pts              Number of points
 *   \param [in]  pts_g_num          Point global ids (size = \ref n_pts)
 *   \param [in]  pts_coord          Point coordinates (size = 3 * \ref n_pts)
 *   \param [in]  n_boxes            Number of boxes
 *   \param [in]  box_g_num          Global ids of boxes (size = \ref n_boxes)
 *   \param [out] pts_in_box_idx     Index of points in boxes (size = \ref n_boxes + 1, allocated inside function)
 *   \param [out] pts_in_box_g_num   Global ids of points in boxes (size = \ref pts_in_box_idx[\ref n_boxes], allocated inside function)
 *   \param [out] pts_in_box_coord   Coordinates of points in boxes (size = 3 *\ref pts_in_box_idx[\ref n_boxes], allocated inside function)
 *   \param [in]  ellipsoids         Consider boxes as axis-aligned ellipsoids (1 or 0)
 *
 */

void
PDM_dbbtree_points_inside_boxes
(
 PDM_dbbtree_t      *dbbt,
 const int           n_pts,
 PDM_g_num_t         pts_g_num[],
 double              pts_coord[],
 const int           n_boxes,
 const PDM_g_num_t   box_g_num[],
 int               **pts_in_box_idx,
 PDM_g_num_t       **pts_in_box_g_num,
 double            **pts_in_box_coord,
 const int           ellipsoids
 )
{
  _PDM_dbbtree_t *_dbbt = (_PDM_dbbtree_t *) dbbt;

  PDM_part_to_block_t *ptb            = NULL;
  int                 *dbox_pts_n     = NULL;
  PDM_g_num_t         *dbox_pts_g_num = NULL;
  double              *dbox_pts_coord = NULL;
  PDM_dbbtree_points_inside_boxes_block_frame(dbbt,
                                              n_pts,
                                              pts_g_num,
                                              pts_coord,
                                              &ptb,
                                              &dbox_pts_n,
                                              &dbox_pts_g_num,
                                              &dbox_pts_coord,
                                              ellipsoids);

  /*
   *  Block to part -> back to origin frame
   */
  int dn_box = PDM_part_to_block_n_elt_block_get(ptb);
  PDM_g_num_t *dbox_g_num = PDM_part_to_block_block_gnum_get(ptb);

  PDM_block_to_part_t *btp = PDM_block_to_part_create_from_sparse_block(dbox_g_num,
                                                                        dn_box,
                                                 (const PDM_g_num_t **) &box_g_num,
                                                                        &n_boxes,
                                                                        1,
                                                                        _dbbt->comm);

  int         **_tmp_pts_in_box_n     = NULL;
  PDM_g_num_t **_tmp_pts_in_box_g_num = NULL;
  PDM_block_to_part_exch(btp,
                         sizeof(PDM_g_num_t),
                         PDM_STRIDE_VAR_INTERLACED,
                         dbox_pts_n,
                         dbox_pts_g_num,
                         &_tmp_pts_in_box_n,
              (void ***) &_tmp_pts_in_box_g_num);
  free(dbox_pts_g_num);

  int *pts_in_box_n = _tmp_pts_in_box_n[0];
  free(_tmp_pts_in_box_n);

  *pts_in_box_idx = PDM_array_new_idx_from_sizes_int(pts_in_box_n, n_boxes);
  free(pts_in_box_n);

  *pts_in_box_g_num = _tmp_pts_in_box_g_num[0];
  free(_tmp_pts_in_box_g_num);


  double **_tmp_pts_in_box_coord = NULL;
  PDM_block_to_part_exch(btp,
                          3 * sizeof(double),
                          PDM_STRIDE_VAR_INTERLACED,
                          dbox_pts_n,
                 (void *) dbox_pts_coord,
                          &_tmp_pts_in_box_n,
               (void ***) &_tmp_pts_in_box_coord);
  free(dbox_pts_n);
  free(dbox_pts_coord);
  free(_tmp_pts_in_box_n[0]);
  free(_tmp_pts_in_box_n);

  *pts_in_box_coord = _tmp_pts_in_box_coord[0];
  free(_tmp_pts_in_box_coord);

  PDM_part_to_block_free(ptb);
  PDM_block_to_part_free(btp);
}



void
PDM_dbbtree_points_inside_boxes_block_frame
(
 PDM_dbbtree_t        *dbbt,
 const int             n_pts,
 PDM_g_num_t           pts_g_num[],
 double                pts_coord[],
 PDM_part_to_block_t **ptb_out,
 int                 **dbox_pts_n,
 PDM_g_num_t         **dbox_pts_g_num,
 double              **dbox_pts_coord,
 const int             ellipsoids
 )
{
  // double t1 = PDM_MPI_Wtime();
  int idebug = 0;

  const float f_threshold = 1.1;  // factor of the mean nb of requests
  const float f_max_copy  = 0.1;  // factor of the total nb of processes

  assert (dbbt != NULL);
  _PDM_dbbtree_t *_dbbt = (_PDM_dbbtree_t *) dbbt;

  int i_rank, n_rank;
  PDM_MPI_Comm_rank (_dbbt->comm, &i_rank);
  PDM_MPI_Comm_size (_dbbt->comm, &n_rank);

  /*
   *  Normalize target points coords
   */
  double *_pts_coord = malloc (sizeof(double) * n_pts * 3);
  for (int i = 0; i < n_pts; i++) {
    _normalize (_dbbt,
                pts_coord  + 3*i,
                _pts_coord + 3*i);
  }

  /*
   *  For each point, find all ranks that might have boxes containing that point
   */
  int *send_count = NULL;
  int *send_shift = NULL;
  int *recv_count = NULL;
  int *recv_shift = NULL;
  PDM_g_num_t *send_g_num = NULL;
  PDM_g_num_t *recv_g_num = NULL;
  double      *send_coord = NULL;
  double      *recv_coord = NULL;

  int n_copied_ranks = 0;
  int *copied_ranks = NULL;
  int n_pts_local  = 0;
  int n_pts_recv   = 0;
  int n_pts_copied = 0;
  int n_pts1;

  int *copied_shift = NULL;

  PDM_g_num_t *pts_g_num1 = NULL;
  double      *pts_coord1 = NULL;

  if (_dbbt->btShared != NULL) {
    int *pts_rank_idx = NULL;
    int *pts_rank = NULL;
    PDM_box_tree_boxes_containing_points(_dbbt->btShared,
                                         -1,
                                         n_pts,
                                         pts_coord,
                                         &pts_rank_idx,
                                         &pts_rank);

    /* Count points to send to each rank */
    send_count = PDM_array_zeros_int (n_rank);

    for (int i = 0; i < pts_rank_idx[n_pts]; i++) {
      int rank = _dbbt->usedRank[pts_rank[i]];
      pts_rank[i] = rank;
      send_count[rank]++;
    }

    recv_count = malloc (sizeof(int) * n_rank);
    PDM_MPI_Alltoall (send_count, 1, PDM_MPI_INT,
                      recv_count, 1, PDM_MPI_INT,
                      _dbbt->comm);

    // PDM_log_trace_array_int(send_count, n_rank, "send_count ::");
    // PDM_log_trace_array_int(recv_count, n_rank, "recv_count ::");

    n_pts_recv = 0;
    for (int i = 0; i < n_rank; i++) {
      n_pts_recv += recv_count[i];
    }
    int n_pts_recv_no_copies = n_pts_recv;

    // log_trace(" n_pts_recv_no_copies = %i \n", n_pts_recv_no_copies);
    /* Prepare copies */
    int n_max_copy = (int) (f_max_copy * n_rank);
    int *i_copied_rank = NULL;
    int *copied_count = NULL;
    int mean_n_pts_recv = 0;
    int *n_pts_recv_copied_ranks = NULL;

    if (n_max_copy > 0) {
      int *all_n_pts_recv = malloc (sizeof(int) * n_rank);
      PDM_MPI_Allgather (&n_pts_recv,    1, PDM_MPI_INT,
                         all_n_pts_recv, 1, PDM_MPI_INT,
                         _dbbt->comm);

      // Mean number of pts_recvs
      long l_mean_n_pts_recv = 0;
      for (int i = 0; i < n_rank; i++) {
        l_mean_n_pts_recv += all_n_pts_recv[i];
      }
      mean_n_pts_recv = (int) (l_mean_n_pts_recv / n_rank);

      float n_threshold = f_threshold * mean_n_pts_recv;

      // Sort ranks
      int *order = malloc (sizeof(int) * n_rank);
      for (int i = 0; i < n_rank; i++) {
        order[i] = i;
      }

      PDM_sort_int (all_n_pts_recv,
                    order,
                    n_rank);

      // Identify ranks to copy
      copied_ranks = malloc (sizeof(int) * n_max_copy);
      n_pts_recv_copied_ranks = malloc (sizeof(int) * n_max_copy);
      for (int i = 0; i < n_max_copy; i++) {
        int j = n_rank - i - 1;

        if (all_n_pts_recv[j] > n_threshold) {
          copied_ranks[n_copied_ranks] = order[j];
          n_pts_recv_copied_ranks[n_copied_ranks] = all_n_pts_recv[j];
          n_copied_ranks++;
        }
        else {
          break;
        }
      }
      free (all_n_pts_recv);
      free (order);

      if (n_copied_ranks > 0) {
        copied_ranks = realloc (copied_ranks, sizeof(int) * n_copied_ranks);
        n_pts_recv_copied_ranks = realloc (n_pts_recv_copied_ranks,
                                           sizeof(int) * n_copied_ranks);
      }
    }

    if (idebug && i_rank == 0) {
      if (n_copied_ranks > 0) {
        if (n_copied_ranks == 1) {
          printf("1 copied rank: %d\n", copied_ranks[0]);
        }
        else {
          printf("%d copied ranks:", n_copied_ranks);
          for (int i = 0; i < n_copied_ranks; i++) {
            printf(" %d", copied_ranks[i]);
          }
          printf("\n");
        }
      }
      else {
        printf("0 copied ranks\n");
      }
    }

    /*
     * Copy the data of selected ranks
     */
    i_copied_rank = malloc (sizeof(int) * n_rank);
    PDM_box_tree_copy_to_ranks (_dbbt->btLoc,
                                &n_copied_ranks,
                                copied_ranks,
                                i_copied_rank);

    /* Re-compute send/recv counts */
    copied_count = PDM_array_zeros_int (n_copied_ranks);
    n_pts_local = 0;
    for (int i = 0; i < n_copied_ranks; i++) {
      int rank = copied_ranks[i];
      if (rank != i_rank) {
        int si = send_count[rank];

        si = PDM_MIN (si, PDM_MAX (0, (n_pts_recv_copied_ranks[i] - n_pts_recv)/2));
        if (i_copied_rank[i_rank] < 0) {
          si = PDM_MIN (si, PDM_MAX (0, mean_n_pts_recv - n_pts_recv));
        }

        copied_count[i] = si;
        n_pts_recv += si;
      }
    }
    if (copied_ranks != NULL) {
      free (copied_ranks);
    }

    if (n_pts_recv_copied_ranks != NULL) {
      free (n_pts_recv_copied_ranks);
    }

    for (int i = 0; i < n_rank; i++) {
      if (i == i_rank) {
        n_pts_local += send_count[i];
        send_count[i] = 0;
      }
      else if (i_copied_rank[i] >= 0) {
        send_count[i] -= copied_count[i_copied_rank[i]];
      }
    }

    copied_shift = PDM_array_new_idx_from_sizes_int (copied_count, n_copied_ranks);
    int *copied_count_tmp = PDM_array_zeros_int (n_copied_ranks);
    n_pts_copied = copied_shift[n_copied_ranks];

    /* Exchange new send/recv counts */
    PDM_MPI_Alltoall (send_count, 1, PDM_MPI_INT,
                      recv_count, 1, PDM_MPI_INT,
                      _dbbt->comm);

    send_shift = PDM_array_new_idx_from_sizes_int (send_count, n_rank);
    recv_shift = PDM_array_new_idx_from_sizes_int (recv_count, n_rank);
    PDM_array_reset_int (send_count, n_rank, 0);

    n_pts_recv = recv_shift[n_rank];
    n_pts1 = n_pts_local + n_pts_recv + n_pts_copied;
    if (idebug) {
      printf("[%d] n_pts1 = %d (without copies : %d)\n", i_rank, n_pts1, n_pts_recv_no_copies);
    }

    pts_g_num1 = malloc (sizeof(PDM_g_num_t) * n_pts1);
    pts_coord1 = malloc (sizeof(double)      * n_pts1 * 3);


    /* Fill send buffers */
    send_g_num = malloc (sizeof(PDM_g_num_t) * send_shift[n_rank]);
    send_coord = malloc (sizeof(double)      * send_shift[n_rank] * 3);

    int idx_copied = n_pts_local + n_pts_recv;
    PDM_g_num_t *copied_g_num = pts_g_num1 + idx_copied;
    double      *copied_coord = pts_coord1 + idx_copied * 3;
    n_pts_local = 0;
    for (int ipt = 0; ipt < n_pts; ipt++) {
      for (int i = pts_rank_idx[ipt]; i < pts_rank_idx[ipt+1]; i++) {
        int rank = pts_rank[i];

        if (rank == i_rank) {
          pts_g_num1[n_pts_local] = pts_g_num[ipt];
          for (int j = 0; j < 3; j++) {
            pts_coord1[3*n_pts_local + j] = _pts_coord[3*ipt + j];
          }
          n_pts_local++;
        }

        else if (i_copied_rank[rank] >= 0) {
          int _rank = i_copied_rank[rank];

          if (copied_count_tmp[_rank] < copied_count[_rank]) {
            int k = copied_shift[_rank] + copied_count_tmp[_rank];
            copied_g_num[k] = pts_g_num[ipt];
            for (int j = 0; j < 3; j++) {
              copied_coord[3*k + j] = _pts_coord[3*ipt + j];
            }
            copied_count_tmp[_rank]++;
          }
          else {
            int k = send_shift[rank] + send_count[rank];
            send_g_num[k] = pts_g_num[ipt];
            for (int j = 0; j < 3; j++) {
              send_coord[3*k + j] = _pts_coord[3*ipt + j];
            }
            send_count[rank]++;
          }
        }

        else {
          int k = send_shift[rank] + send_count[rank];
          send_g_num[k] = pts_g_num[ipt];
          for (int j = 0; j < 3; j++) {
            send_coord[3*k + j] = _pts_coord[3*ipt + j];
          }
          send_count[rank]++;
        }
      }
    }
    free (_pts_coord);
    if (copied_count != NULL) {
      free (copied_count);
    }
    if (copied_count_tmp != NULL) {
      free (copied_count_tmp);
    }
    if (i_copied_rank != NULL) {
      free (i_copied_rank);
    }
    free (pts_rank);
    free (pts_rank_idx);


    /* Exchange points */
    recv_g_num = pts_g_num1 + n_pts_local;
    PDM_MPI_Alltoallv (send_g_num, send_count, send_shift, PDM__PDM_MPI_G_NUM,
                       recv_g_num, recv_count, recv_shift, PDM__PDM_MPI_G_NUM,
                       _dbbt->comm);
    free (send_g_num);

    for (int i = 0; i < n_rank; i++) {
      send_count[i] *= 3;
      recv_count[i] *= 3;
      send_shift[i+1] *= 3;
      recv_shift[i+1] *= 3;
    }

    recv_coord = pts_coord1 + 3*n_pts_local;
    PDM_MPI_Alltoallv (send_coord, send_count, send_shift, PDM_MPI_DOUBLE,
                       recv_coord, recv_count, recv_shift, PDM_MPI_DOUBLE,
                       _dbbt->comm);

    free (send_coord);
    free (send_count);
    free (send_shift);
    free (recv_count);
    free (recv_shift);
  }

  else {
    n_pts_local  = n_pts;
    n_pts_recv   = 0;
    n_pts_copied = 0;

    n_pts1 = n_pts;

    pts_g_num1 = (PDM_g_num_t *) pts_g_num;
    pts_coord1 = _pts_coord;
  }



  /*
   *  Get boxes containing points
   */
  int n_part = 1 + n_copied_ranks;
  int *part_n_pts = malloc (sizeof(int) * n_part);
  part_n_pts[0] = n_pts_local + n_pts_recv;

  int **pts_box_idx   = malloc (sizeof(int *) * n_part);
  int **pts_box_l_num = malloc (sizeof(int *) * n_part);

  int size_pts_box = 0;

  void (*_points_inside_volumes) (PDM_box_tree_t *bt,
                                  const int,
                                  const int,
                                  const double *,
                                  int **,
                                  int **);
  if (ellipsoids) {
    _points_inside_volumes = &PDM_box_tree_ellipsoids_containing_points;
  } else {
    _points_inside_volumes = &PDM_box_tree_boxes_containing_points;
  }

  /*
   *  Search in local tree
   */
  // double t1 = PDM_MPI_Wtime();
  _points_inside_volumes (_dbbt->btLoc,
                          -1,
                          part_n_pts[0],
                          pts_coord1,
                          &(pts_box_idx[0]),
                          &(pts_box_l_num[0]));
  size_pts_box += pts_box_idx[0][part_n_pts[0]];
  // double t2 = PDM_MPI_Wtime()-t1;
  // log_trace("Point in boxes  = %12.5e \n", t2);


  /*
   *  Search in copied trees
   */
  double      *pts_coord_copied = NULL;
  PDM_g_num_t *pts_g_num_copied = NULL;
  if (n_copied_ranks > 0) {
    pts_coord_copied = pts_coord1 + part_n_pts[0] * 3;
    pts_g_num_copied = pts_g_num1 + part_n_pts[0];
    for (int i = 0; i < n_copied_ranks; i++) {
      part_n_pts[i+1] = copied_shift[i+1] - copied_shift[i];

      _points_inside_volumes (_dbbt->btLoc,
                              i,
                              part_n_pts[i+1],
                              pts_coord_copied + copied_shift[i] * 3,
                              &(pts_box_idx[i+1]),
                              &(pts_box_l_num[i+1]));

      size_pts_box += pts_box_idx[i+1][part_n_pts[i+1]];
    }
  }
  if (copied_shift != NULL) free (copied_shift);


  PDM_g_num_t *part_box_g_num = malloc (sizeof(PDM_g_num_t) * size_pts_box);
  PDM_g_num_t *part_pts_g_num = malloc (sizeof(PDM_g_num_t) * size_pts_box);
  double      *part_pts_coord = malloc (sizeof(double)      * size_pts_box * 3);
  int idx = 0;
  for (int j = 0; j < part_n_pts[0]; j++) {
    for (int k = pts_box_idx[0][j]; k < pts_box_idx[0][j+1]; k++) {
      part_box_g_num[idx] = _dbbt->boxes->local_boxes->g_num[pts_box_l_num[0][k]];
      part_pts_g_num[idx] = pts_g_num1[j];
      for (int l = 0; l < 3; l++) {
        part_pts_coord[3*idx + l] = pts_coord1[3*j + l];
      }
      idx++;
    }
  }
  free (pts_box_idx[0]);
  free (pts_box_l_num[0]);

  for (int i = 0; i < n_copied_ranks; i++) {
    for (int j = 0; j < part_n_pts[i+1]; j++) {
      for (int k = pts_box_idx[i+1][j]; k < pts_box_idx[i+1][j+1]; k++) {
        part_box_g_num[idx] = _dbbt->boxes->rank_boxes[i].g_num[pts_box_l_num[i+1][k]];
        part_pts_g_num[idx] = pts_g_num_copied[j];
        for (int l = 0; l < 3; l++) {
          part_pts_coord[3*idx + l] = pts_coord_copied[3*j + l];
        }
        idx++;
      }
    }
    pts_coord_copied += part_n_pts[i+1] * 3;
    pts_g_num_copied += part_n_pts[i+1];
    free (pts_box_idx[i+1]);
    free (pts_box_l_num[i+1]);
  }
  // PDM_log_trace_array_int(part_n_pts, n_part, "part_n_pts ::");
  free (part_n_pts);
  free (pts_box_idx);
  free (pts_box_l_num);
  free (pts_coord1);
  if (pts_g_num1 != pts_g_num) free (pts_g_num1);

  /*
   *  Part-to-block
   */
  int    *part_stride = malloc (sizeof(int)    * size_pts_box);
  double *weight      = malloc (sizeof(double) * size_pts_box);
  for (int i = 0; i < size_pts_box; i++) {
    part_stride[i] = 1;
    weight[i]      = 1.;
  }

  // double dt = PDM_MPI_Wtime()-t1;
  // log_trace("point_in_boxes : %12.5e \n", dt);
  // PDM_MPI_Barrier(_dbbt->comm);

  PDM_part_to_block_t *ptb = PDM_part_to_block_create (PDM_PART_TO_BLOCK_DISTRIB_ALL_PROC,
                                                       PDM_PART_TO_BLOCK_POST_MERGE,
                                                       1.,
                                                       &part_box_g_num,
                                                       &weight,
                                                       &size_pts_box,
                                                       1,
                                                       _dbbt->comm);
  free(weight);
  free(part_box_g_num);

  int         *block_box_pts_n     = NULL;
  PDM_g_num_t *block_box_pts_g_num = NULL;
  PDM_part_to_block_exch(ptb,
                         sizeof(PDM_g_num_t),
                         PDM_STRIDE_VAR_INTERLACED,
                         1,
                         &part_stride,
               (void **) &part_pts_g_num,
                         &block_box_pts_n,
               (void **) &block_box_pts_g_num);
  free (part_pts_g_num);


  int    *block_stride        = NULL;
  double *block_box_pts_coord = NULL;
  PDM_part_to_block_exch(ptb,
                         3*sizeof(double),
                         PDM_STRIDE_VAR_INTERLACED,
                         1,
                         &part_stride,
               (void **) &part_pts_coord,
                         &block_stride,
               (void **) &block_box_pts_coord);
  free (block_stride);
  free (part_stride);
  free (part_pts_coord);
  PDM_box_tree_free_copies(_dbbt->btLoc);

  /* De-normalize pts coord */
  int dn_box = PDM_part_to_block_n_elt_block_get(ptb);

  idx = 0;
  for (int i = 0; i < dn_box; i++) {
    for (int j = 0; j < block_box_pts_n[i]; j++) {
      for (int k = 0; k < 3; k++) {
        block_box_pts_coord[idx] = _dbbt->s[k] + _dbbt->d[k] * (block_box_pts_coord[idx]);
        idx++;
      }
    }
  }

  *ptb_out = ptb;
  *dbox_pts_n     = block_box_pts_n;
  *dbox_pts_g_num = block_box_pts_g_num;
  *dbox_pts_coord = block_box_pts_coord;
}

/**
 *
 * \brief Get an indexed list of all points inside the boxes of a distributed box tree
 *
 *   \param [in]  dbbt               Pointer to distributed box tree structure
 *   \param [in]  n_pts              Number of points
 *   \param [in]  pts_g_num          Point global ids (size = \ref n_pts)
 *   \param [in]  pts_coord          Point coordinates (size = 3 * \ref n_pts)
 *   \param [in]  n_boxes            Number of boxes
 *   \param [in]  box_g_num          Global ids of boxes (size = \ref n_boxes)
 *   \param [out] pts_in_box_idx     Index of points in boxes (size = \ref n_boxes + 1, allocated inside function)
 *   \param [out] pts_in_box_g_num   Global ids of points in boxes (size = \ref pts_in_box_idx[\ref n_boxes], allocated inside function)
 *   \param [out] pts_in_box_coord   Coordinates of points in boxes (size = 3 *\ref pts_in_box_idx[\ref n_boxes], allocated inside function)
 *   \param [in]  ellipsoids         Consider boxes as axis-aligned ellipsoids (1 or 0)
 *
 */

void
PDM_dbbtree_points_inside_boxes_shared
(
 PDM_dbbtree_t      *dbbt,
 const int           n_pts,
 PDM_g_num_t         pts_g_num[],
 double              pts_coord[],
 const int           n_boxes,
 const PDM_g_num_t   box_g_num[],
 int               **pts_in_box_idx,
 PDM_g_num_t       **pts_in_box_g_num,
 double            **pts_in_box_coord,
 const int           ellipsoids
)
{
  // double t1 = PDM_MPI_Wtime();
  // int idebug = 0;

  assert (dbbt != NULL);
  _PDM_dbbtree_t *_dbbt = (_PDM_dbbtree_t *) dbbt;

  int i_rank, n_rank;
  PDM_MPI_Comm_rank (_dbbt->comm, &i_rank);
  PDM_MPI_Comm_size (_dbbt->comm, &n_rank);

  // Shared
  PDM_MPI_Comm comm_shared;
  PDM_MPI_Comm_split_type(_dbbt->comm, PDM_MPI_SPLIT_NUMA, &comm_shared);

  int n_rank_in_shm, i_rank_in_shm;
  PDM_MPI_Comm_rank (comm_shared, &i_rank_in_shm);
  PDM_MPI_Comm_size (comm_shared, &n_rank_in_shm);

  /*
   *  Normalize target points coords
   */
  double *_pts_coord = malloc (sizeof(double) * n_pts * 3);
  for (int i = 0; i < n_pts; i++) {
    _normalize (_dbbt, pts_coord  + 3*i, _pts_coord + 3*i);
  }

  // PDM_log_trace_array_double(pts_coord, 3 * n_pts, "pts_coord ::");
  // PDM_log_trace_array_long  (pts_g_num,     n_pts, "pts_g_num ::");

  int n_pts_local  = 0;
  int n_pts_recv   = 0;
  int n_pts_copied = 0;
  int n_pts1       = 0;

  PDM_UNUSED(n_pts_recv);
  PDM_UNUSED(n_pts_copied);
  PDM_UNUSED(n_pts1);
  PDM_g_num_t *pts_g_num1 = NULL;
  double      *pts_coord1 = NULL;

  PDM_mpi_win_shared_t* wshared_recv_gnum      = NULL;
  PDM_mpi_win_shared_t* wshared_recv_pts_coord = NULL;

  int *distrib_search_by_rank_idx = NULL;

  /*
   *  For each point, find all ranks that might have boxes containing that point
   */
  if (_dbbt->btShared != NULL) {
    int *send_count = NULL;
    int *send_shift = NULL;
    int *recv_count = NULL;
    int *recv_shift = NULL;
    PDM_g_num_t *send_g_num = NULL;
    // PDM_g_num_t *recv_g_num = NULL;
    double      *send_coord = NULL;
    // double      *recv_coord = NULL;

    int *pts_rank_idx = NULL;
    int *pts_rank     = NULL;
    PDM_box_tree_boxes_containing_points(_dbbt->btShared,
                                         -1,
                                         n_pts,
                                         pts_coord,
                                         &pts_rank_idx,
                                         &pts_rank);

    /* Count points to send to each rank */
    send_count = PDM_array_zeros_int (n_rank);

    for (int i = 0; i < pts_rank_idx[n_pts]; i++) {
      int t_rank = _dbbt->usedRank[pts_rank[i]];
      pts_rank[i] = t_rank;
      send_count[t_rank]++;
    }

    recv_count = malloc (sizeof(int) * n_rank);
    PDM_MPI_Alltoall (send_count, 1, PDM_MPI_INT,
                      recv_count, 1, PDM_MPI_INT,
                      _dbbt->comm);

    send_shift = malloc ( ( n_rank + 1) * sizeof(int));
    recv_shift = malloc ( ( n_rank + 1) * sizeof(int));

    // Deduce size of recv buffer shared inside the same node
    int* shared_recv_count  = malloc(n_rank_in_shm * sizeof(int));

    int n_tot_recv = 0;
    send_shift[0] = 0;
    recv_shift[0] = 0;
    for(int i = 0; i < n_rank; ++i) {
      n_tot_recv += recv_count[i];
      send_shift[i+1] = send_shift[i] + send_count[i];
      recv_shift[i+1] = recv_shift[i] + recv_count[i];
    }

    PDM_MPI_Allgather(&n_tot_recv      , 1, PDM_MPI_INT,
                      shared_recv_count, 1, PDM_MPI_INT,
                      comm_shared);


    int* shared_recv_idx = malloc((n_rank_in_shm+1) * sizeof(int));
    shared_recv_idx[0] = 0;
    for(int i = 0; i < n_rank_in_shm; ++i) {
      shared_recv_idx[i+1] = shared_recv_idx[i] + shared_recv_count[i];
    }

    if(0 == 1) {
      PDM_log_trace_array_int(shared_recv_count , n_rank_in_shm, "shared_recv_count  :: ");
      PDM_log_trace_array_int(shared_recv_idx , n_rank_in_shm+1, "shared_recv_idx  :: ");
    }

    int n_tot_recv_shared = shared_recv_idx[n_rank_in_shm];
    wshared_recv_gnum      = PDM_mpi_win_shared_create(    n_tot_recv_shared, sizeof(PDM_g_num_t), comm_shared);
    wshared_recv_pts_coord = PDM_mpi_win_shared_create(3 * n_tot_recv_shared, sizeof(double     ), comm_shared);

    PDM_g_num_t *shared_recv_gnum      = PDM_mpi_win_shared_get(wshared_recv_gnum);
    double      *shared_recv_pts_coord = PDM_mpi_win_shared_get(wshared_recv_pts_coord);

    PDM_mpi_win_shared_lock_all (0, wshared_recv_gnum);
    PDM_mpi_win_shared_lock_all (0, wshared_recv_pts_coord);

    PDM_g_num_t *lrecv_gnum      = &shared_recv_gnum     [    shared_recv_idx[i_rank_in_shm]];
    double      *lrecv_pts_coord = &shared_recv_pts_coord[3 * shared_recv_idx[i_rank_in_shm]];

    /* Prepare send */
    send_g_num = malloc (sizeof(PDM_g_num_t) * send_shift[n_rank]);
    send_coord = malloc (sizeof(double)      * send_shift[n_rank] * 3);

    for(int i = 0; i < n_rank; ++i) {
      send_count[i] = 0;
    }

    for (int ipt = 0; ipt < n_pts; ipt++) {
      for (int i = pts_rank_idx[ipt]; i < pts_rank_idx[ipt+1]; i++) {
        int t_rank = pts_rank[i];

        int idx_write = send_shift[t_rank] + send_count[t_rank]++;
        send_g_num[idx_write] = pts_g_num[ipt];
        for (int j = 0; j < 3; j++) {
          send_coord[3*idx_write + j] = _pts_coord[3*ipt + j];
        }
      }
    }


    free(pts_rank_idx);
    free(pts_rank);


    PDM_MPI_Alltoallv (send_g_num, send_count, send_shift, PDM__PDM_MPI_G_NUM,
                       lrecv_gnum, recv_count, recv_shift, PDM__PDM_MPI_G_NUM,
                       _dbbt->comm);
    free(send_g_num);

    for (int i = 0; i < n_rank; i++) {
      send_count[i  ] *= 3;
      recv_count[i  ] *= 3;
      send_shift[i+1] *= 3;
      recv_shift[i+1] *= 3;
    }

    PDM_MPI_Alltoallv (send_coord     , send_count, send_shift, PDM_MPI_DOUBLE,
                       lrecv_pts_coord, recv_count, recv_shift, PDM_MPI_DOUBLE,
                       _dbbt->comm);

    free (send_coord);
    free (send_count);
    free (send_shift);
    free (recv_count);
    free (recv_shift);

    PDM_MPI_Barrier (comm_shared);
    PDM_mpi_win_shared_sync (wshared_recv_gnum);
    PDM_mpi_win_shared_sync (wshared_recv_pts_coord);

    /*
     * Redistribution by numa
     */

    PDM_g_num_t* distrib_search = PDM_compute_uniform_entity_distribution(comm_shared, n_tot_recv_shared);

    int  dn_search = distrib_search[i_rank_in_shm+1] - distrib_search[i_rank_in_shm];

    distrib_search_by_rank_idx = malloc((n_rank_in_shm+1) * sizeof(int));
    int* distrib_search_by_rank_n   = malloc((n_rank_in_shm  ) * sizeof(int));
    for(int i = 0; i < n_rank_in_shm; ++i) {
      distrib_search_by_rank_n[i] = 0;
    }

    // TODO : Faire un algo d'intersection de range pour ne pas faire la dicotomie x fois !
    for(int i = distrib_search[i_rank_in_shm]; i < distrib_search[i_rank_in_shm+1]; ++i) {
      int t_rank = PDM_binary_search_gap_int(i, shared_recv_idx, n_rank_in_shm+1);
      distrib_search_by_rank_n[t_rank]++;
    }

    distrib_search_by_rank_idx[0] = 0;
    for(int i = 0; i < n_rank_in_shm; ++i) {
      distrib_search_by_rank_idx[i+1] = distrib_search_by_rank_idx[i] + distrib_search_by_rank_n[i];
    }

    // PDM_log_trace_array_long(check, dn_search, "check ::");
    // PDM_log_trace_array_int(distrib_search_by_rank_idx, n_rank_in_shm+1, "distrib_search_by_rank_idx ::");

    n_pts_local  = 0;
    n_pts_recv   = dn_search;
    n_pts_copied = 0;

    n_pts1       = dn_search;
    pts_g_num1   = (PDM_g_num_t *) &shared_recv_gnum     [    distrib_search[i_rank_in_shm]];
    pts_coord1   = (double      *) &shared_recv_pts_coord[3 * distrib_search[i_rank_in_shm]];

    free(distrib_search);
    free(distrib_search_by_rank_n);
    free(shared_recv_count );
    free(shared_recv_idx );


  } else {
    n_pts_local  = n_pts;
    n_pts_recv   = 0;
    n_pts_copied = 0;

    n_pts1 = n_pts;

    pts_g_num1 = (PDM_g_num_t *) pts_g_num;
    pts_coord1 = _pts_coord;
  }

  void (*_points_inside_volumes) (PDM_box_tree_t *bt,
                                  const int,
                                  const int,
                                  const double *,
                                  int **,
                                  int **);
  if (ellipsoids) {
    _points_inside_volumes = &PDM_box_tree_ellipsoids_containing_points;
  } else {
    _points_inside_volumes = &PDM_box_tree_boxes_containing_points;
  }

  /* Management of size */
  int n_part = n_rank_in_shm;
  int *part_n_pts_box = malloc (sizeof(int) * n_part);
  int *part_n_pts     = malloc (sizeof(int) * n_part);

  int **pts_box_idx   = malloc (sizeof(int *) * n_part);
  int **pts_box_l_num = malloc (sizeof(int *) * n_part);

  PDM_g_num_t **part_box_g_num  = malloc (sizeof(PDM_g_num_t *) * n_part);
  PDM_g_num_t **part_pts_g_num  = malloc (sizeof(PDM_g_num_t *) * n_part);
  double      **part_pts_coord  = malloc (sizeof(double      *) * n_part);
  int         **part_pts_strid  = malloc (sizeof(int         *) * n_part);
  double      **part_pts_weight = malloc (sizeof(double      *) * n_part);

  // double t1 = PDM_MPI_Wtime();
  if(_dbbt->btShared != NULL) {
    for(int i_shm = 0; i_shm < n_rank_in_shm; ++i_shm) {

      int beg    = distrib_search_by_rank_idx[i_shm  ];
      int n_lpts = distrib_search_by_rank_idx[i_shm+1] - beg;

      part_n_pts    [i_shm] = n_lpts;

      PDM_g_num_t *lpts_gnum  = &pts_g_num1[  beg];
      double      *lpts_coord = &pts_coord1[3*beg];

      if(0 == 1) {
        PDM_log_trace_array_long  (lpts_gnum,    n_lpts, "lpts_gnum ::");
        PDM_log_trace_array_double(lpts_coord, 3*n_lpts, "lpts_coord::");
      }

      PDM_box_tree_boxes_containing_points_shared(_dbbt->btLoc,
                                                  i_shm,
                                                  part_n_pts[i_shm],
                                                  lpts_coord,
                                                  &(pts_box_idx[i_shm]),
                                                  &(pts_box_l_num[i_shm]));
    // _points_inside_volumes (_dbbt->btLoc,
    //                         -1,
    //                         part_n_pts[i_shm],
    //                                               lpts_coord,
    //                                               &(pts_box_idx[i_shm]),
    //                                               &(pts_box_l_num[i_shm]));

      part_n_pts_box[i_shm] = pts_box_idx[i_shm][part_n_pts[i_shm]];

      // PDM_log_trace_connectivity_int(pts_box_idx[i_shm], pts_box_l_num[i_shm], part_n_pts[i_shm], "pts_box_l_num ::");

      part_box_g_num [i_shm] = malloc (sizeof(PDM_g_num_t) *     part_n_pts_box[i_shm]);
      part_pts_g_num [i_shm] = malloc (sizeof(PDM_g_num_t) *     part_n_pts_box[i_shm]);
      part_pts_coord [i_shm] = malloc (sizeof(double     ) * 3 * part_n_pts_box[i_shm]);
      part_pts_strid [i_shm] = malloc (sizeof(int        ) *     part_n_pts_box[i_shm]);
      part_pts_weight[i_shm] = malloc (sizeof(double     ) *     part_n_pts_box[i_shm]);

      for (int j = 0; j < part_n_pts[i_shm]; j++) {
        for (int k = pts_box_idx[i_shm][j]; k < pts_box_idx[i_shm][j+1]; k++) {
          part_box_g_num[i_shm][k] = _dbbt->boxes->shm_boxes[i_shm].g_num[pts_box_l_num[i_shm][k]];
          // part_box_g_num[i_shm][k] = _dbbt->boxes->local_boxes->g_num[pts_box_l_num[i_shm][k]];
          part_pts_g_num[i_shm][k] = lpts_gnum[j];
          for (int l = 0; l < 3; l++) {
            part_pts_coord[i_shm][3*k + l] = lpts_coord[3*j + l];
          }
          part_pts_strid [i_shm][k]= 1;
          part_pts_weight[i_shm][k]= 1.;
        }
      }
    }

    PDM_mpi_win_shared_unlock_all(wshared_recv_gnum);
    PDM_mpi_win_shared_unlock_all(wshared_recv_pts_coord);
    PDM_mpi_win_shared_free (wshared_recv_gnum);
    PDM_mpi_win_shared_free (wshared_recv_pts_coord);
    free(distrib_search_by_rank_idx);
  } else {
    part_n_pts[0] = n_pts_local + n_pts_recv;
    _points_inside_volumes (_dbbt->btLoc,
                            -1,
                            part_n_pts[0],
                            pts_coord1,
                            &(pts_box_idx[0]),
                            &(pts_box_l_num[0]));

    part_n_pts_box [0] = pts_box_idx[0][part_n_pts[0]];
    part_box_g_num [0] = malloc (sizeof(PDM_g_num_t) *     part_n_pts_box[0]);
    part_pts_g_num [0] = malloc (sizeof(PDM_g_num_t) *     part_n_pts_box[0]);
    part_pts_coord [0] = malloc (sizeof(double     ) * 3 * part_n_pts_box[0]);
    part_pts_strid [0] = malloc (sizeof(int        ) *     part_n_pts_box[0]);
    part_pts_weight[0] = malloc (sizeof(double     ) *     part_n_pts_box[0]);

    for (int j = 0; j < part_n_pts[0]; j++) {
      for (int k = pts_box_idx[0][j]; k < pts_box_idx[0][j+1]; k++) {
        part_box_g_num[0][k] = _dbbt->boxes->local_boxes->g_num[pts_box_l_num[0][k]];
        part_pts_g_num[0][k] = pts_g_num1[j];
        for (int l = 0; l < 3; l++) {
          part_pts_coord[0][3*k + l] = pts_coord1[3*j + l];
        }
        part_pts_strid [0][k]= 1;
        part_pts_weight[0][k]= 1.;
      }
    }
  }

  // double dt = PDM_MPI_Wtime() - t1;
  // log_trace("PDM_box_tree_boxes_containing_points_shared sol = %12.5e \n", dt);


  PDM_part_to_block_t *ptb = PDM_part_to_block_create (PDM_PART_TO_BLOCK_DISTRIB_ALL_PROC,
                                                       PDM_PART_TO_BLOCK_POST_MERGE,
                                                       1.,
                                                       part_box_g_num,
                                                       part_pts_weight,
                                                       part_n_pts_box,
                                                       n_part,
                                                       _dbbt->comm);
  for(int i = 0; i < n_part; ++i) {
    free(part_pts_weight[i]);
  }
  free(part_pts_weight);

  /*
   *  Exchange gnum
   */
  int         *block_box_pts_n = NULL;
  PDM_g_num_t *block_box_pts_g_num = NULL;
  PDM_part_to_block_exch (ptb,
                          sizeof(PDM_g_num_t),
                          PDM_STRIDE_VAR_INTERLACED,
                          1,
                          part_pts_strid,
                (void **) part_pts_g_num,
                          &block_box_pts_n,
                (void **) &block_box_pts_g_num);

  for(int i = 0; i < n_part; ++i) {
    free(part_box_g_num [i]);
  }
  free(part_box_g_num);

  int    *block_stride = NULL;
  double *block_box_pts_coord = NULL;
  PDM_part_to_block_exch (ptb,
                          3*sizeof(double),
                          PDM_STRIDE_VAR_INTERLACED,
                          1,
                          part_pts_strid,
                (void **) part_pts_coord,
                          &block_stride,
                (void **) &block_box_pts_coord);
  free(block_stride);

  int n_elt_block = PDM_part_to_block_n_elt_block_get (ptb);
  PDM_g_num_t* blk_gnum = PDM_part_to_block_block_gnum_get(ptb);

  if(0 == 1) {
    int size_tot = 0;
    for(int i = 0; i < n_elt_block; ++i) {
      size_tot += block_box_pts_n[i];
    }
    PDM_log_trace_array_long(blk_gnum, n_elt_block, "blk_gnum :: ");
    PDM_log_trace_array_long(block_box_pts_g_num, size_tot, "block_box_pts_g_num :: ");
    PDM_log_trace_array_int (block_box_pts_n    , n_elt_block, "block_box_pts_n :: ");
  }


  PDM_block_to_part_t *btp = PDM_block_to_part_create_from_sparse_block(blk_gnum,
                                                                        n_elt_block,
                                                (const PDM_g_num_t **) &box_g_num,
                                                                       &n_boxes,
                                                                       1,
                                                                       _dbbt->comm);

  int         **_tmp_pts_in_box_n     = NULL;
  PDM_g_num_t **_tmp_pts_in_box_g_num = NULL;
  PDM_block_to_part_exch(btp,
                         sizeof(PDM_g_num_t),
                         PDM_STRIDE_VAR_INTERLACED,
                         block_box_pts_n,
                         block_box_pts_g_num,
                         &_tmp_pts_in_box_n,
              (void ***) &_tmp_pts_in_box_g_num);
  free (block_box_pts_g_num);

  *pts_in_box_g_num = _tmp_pts_in_box_g_num[0];
  free(_tmp_pts_in_box_g_num);

  free(_tmp_pts_in_box_n[0]);
  free(_tmp_pts_in_box_n);

  /*
   * Coordinates
   */
  double **_tmp_pts_in_box_coord = NULL;
  PDM_block_to_part_exch(btp,
                         3*sizeof(double),
                         PDM_STRIDE_VAR_INTERLACED,
                         block_box_pts_n,
                         block_box_pts_coord,
                         &_tmp_pts_in_box_n,
              (void ***) &_tmp_pts_in_box_coord);
  *pts_in_box_coord = _tmp_pts_in_box_coord[0];
  free(_tmp_pts_in_box_coord);

  *pts_in_box_idx = PDM_array_new_idx_from_sizes_int(_tmp_pts_in_box_n[0], n_boxes);
  free(_tmp_pts_in_box_n[0]);
  free(_tmp_pts_in_box_n);

  free(block_box_pts_n);
  free(block_box_pts_coord);


  PDM_part_to_block_free(ptb);
  PDM_block_to_part_free(btp);

  for(int i = 0; i < n_part; ++i) {
    free(pts_box_idx    [i]);
    free(pts_box_l_num  [i]);
    free(part_pts_g_num [i]);
    free(part_pts_coord [i]);
    free(part_pts_strid [i]);
  }

  free(pts_box_idx    );
  free(pts_box_l_num  );
  free(part_pts_g_num );
  free(part_pts_coord );
  free(part_pts_strid );

  free(part_n_pts_box);
  free(part_n_pts    );

  //-->>
  for (int i = 0; i < (*pts_in_box_idx)[n_boxes]; i++) {
    for (int j = 0; j < 3; j++) {
      (*pts_in_box_coord)[3*i+j] = _dbbt->s[j] + _dbbt->d[j] * ((*pts_in_box_coord)[3*i+j]);
    }
  }
  //<<--

  free(_pts_coord);
  PDM_MPI_Comm_free(&comm_shared);
}

/**
 *
 * \brief Get an indexed list of all boxes containing points
 *
 *   \param [in]  dbbt        Pointer to distributed box tree structure
 *   \param [in]  n_pts       Number of points
 *   \param [in]  pts_g_num   Point global ids (size = \ref n_pts)
 *   \param [in]  pts_coord   Point coordinates (size = 3 * \ref n_pts)
 *   \param [out] box_idx     Index of boxes (size = \ref n_pts + 1, allocated inside function)
 *   \param [out] box_g_num   Global ids of boxes (size = \ref box_idx[\ref n_pts], allocated inside function)
 *   \param [in]  ellipsoids  Consider boxes as axis-aligned ellipsoids (1 or 0)
 *
 */

void
PDM_dbbtree_boxes_containing_points
(
 PDM_dbbtree_t      *dbbt,
 const int           n_pts,
 PDM_g_num_t         pts_g_num[],
 double              pts_coord[],
 int               **box_idx,
 PDM_g_num_t       **box_g_num,
 const int           ellipsoids
 )
{
  int idebug = 0;

  const float f_threshold = 1.1;  // factor of the mean nb of requests
  const float f_max_copy  = 0.1;  // factor of the total nb of processes

  assert (dbbt != NULL);
  _PDM_dbbtree_t *_dbbt = (_PDM_dbbtree_t *) dbbt;

  int i_rank, n_rank;
  PDM_MPI_Comm_rank (_dbbt->comm, &i_rank);
  PDM_MPI_Comm_size (_dbbt->comm, &n_rank);

  /*
   *  Normalize target points coords
   */
  double *_pts_coord = malloc (sizeof(double) * n_pts * 3);
  for (int i = 0; i < n_pts; i++) {
    _normalize (_dbbt,
                pts_coord  + 3*i,
                _pts_coord + 3*i);
  }

  /*
   *  For each point, find all ranks that might have boxes containing that point
   */
  int *send_count = NULL;
  int *send_shift = NULL;
  int *recv_count = NULL;
  int *recv_shift = NULL;
  PDM_g_num_t *send_g_num = NULL;
  PDM_g_num_t *recv_g_num = NULL;
  double      *send_coord = NULL;
  double      *recv_coord = NULL;

  int n_copied_ranks = 0;
  int *copied_ranks = NULL;
  int n_pts_local  = 0;
  int n_pts_recv   = 0;
  int n_pts_copied = 0;
  int n_pts1;

  int *copied_shift = NULL;

  PDM_g_num_t *pts_g_num1 = NULL;
  double      *pts_coord1 = NULL;

  if (_dbbt->btShared != NULL) {
    int *pts_rank_idx = NULL;
    int *pts_rank = NULL;
    PDM_box_tree_boxes_containing_points (_dbbt->btShared,
                                       -1,
                                       n_pts,
                                       pts_coord,
                                       &pts_rank_idx,
                                       &pts_rank);

    /* Count points to send to each rank */
    send_count = PDM_array_zeros_int (n_rank);

    for (int i = 0; i < pts_rank_idx[n_pts]; i++) {
      int rank = _dbbt->usedRank[pts_rank[i]];
      pts_rank[i] = rank;
      send_count[rank]++;
    }

    recv_count = malloc (sizeof(int) * n_rank);
    PDM_MPI_Alltoall (send_count, 1, PDM_MPI_INT,
                      recv_count, 1, PDM_MPI_INT,
                      _dbbt->comm);

    n_pts_recv = 0;
    for (int i = 0; i < n_rank; i++) {
      n_pts_recv += recv_count[i];
    }
    int n_pts_recv_no_copies = n_pts_recv;

    /* Prepare copies */
    int n_max_copy = (int) (f_max_copy * n_rank);
    int *i_copied_rank = NULL;
    int *copied_count = NULL;
    int mean_n_pts_recv = 0;
    int *n_pts_recv_copied_ranks = NULL;

    if (n_max_copy > 0) {
      int *all_n_pts_recv = malloc (sizeof(int) * n_rank);
      PDM_MPI_Allgather (&n_pts_recv,    1, PDM_MPI_INT,
                         all_n_pts_recv, 1, PDM_MPI_INT,
                         _dbbt->comm);

      // Mean number of pts_recvs
      long l_mean_n_pts_recv = 0;
      for (int i = 0; i < n_rank; i++) {
        l_mean_n_pts_recv += all_n_pts_recv[i];
      }
      mean_n_pts_recv = (int) (l_mean_n_pts_recv / n_rank);

      float n_threshold = f_threshold * mean_n_pts_recv;

      // Sort ranks
      int *order = malloc (sizeof(int) * n_rank);
      for (int i = 0; i < n_rank; i++) {
        order[i] = i;
      }

      PDM_sort_int (all_n_pts_recv,
                    order,
                    n_rank);

      // Identify ranks to copy
      copied_ranks = malloc (sizeof(int) * n_max_copy);
      n_pts_recv_copied_ranks = malloc (sizeof(int) * n_max_copy);
      for (int i = 0; i < n_max_copy; i++) {
        int j = n_rank - i - 1;

        if (all_n_pts_recv[j] > n_threshold) {
          copied_ranks[n_copied_ranks] = order[j];
          n_pts_recv_copied_ranks[n_copied_ranks] = all_n_pts_recv[j];
          n_copied_ranks++;
        }
        else {
          break;
        }
      }
      free (all_n_pts_recv);
      free (order);

      if (n_copied_ranks > 0) {
        copied_ranks = realloc (copied_ranks, sizeof(int) * n_copied_ranks);
        n_pts_recv_copied_ranks = realloc (n_pts_recv_copied_ranks,
                                           sizeof(int) * n_copied_ranks);
      }
    }

    if (idebug && i_rank == 0) {
      if (n_copied_ranks > 0) {
        if (n_copied_ranks == 1) {
          printf("1 copied rank: %d\n", copied_ranks[0]);
        }
        else {
          printf("%d copied ranks:", n_copied_ranks);
          for (int i = 0; i < n_copied_ranks; i++) {
            printf(" %d", copied_ranks[i]);
          }
          printf("\n");
        }
      }
      else {
        printf("0 copied ranks\n");
      }
    }

    /*
     * Copy the data of selected ranks
     */
    i_copied_rank = malloc (sizeof(int) * n_rank);
    PDM_box_tree_copy_to_ranks (_dbbt->btLoc,
                                &n_copied_ranks,
                                copied_ranks,
                                i_copied_rank);

    /* Re-compute send/recv counts */
    copied_count = PDM_array_zeros_int (n_copied_ranks);
    n_pts_local = 0;
    for (int i = 0; i < n_copied_ranks; i++) {
      int rank = copied_ranks[i];
      if (rank != i_rank) {
        int si = send_count[rank];

        si = PDM_MIN (si, PDM_MAX (0, (n_pts_recv_copied_ranks[i] - n_pts_recv)/2));
        if (i_copied_rank[i_rank] < 0) {
          si = PDM_MIN (si, PDM_MAX (0, mean_n_pts_recv - n_pts_recv));
        }

        copied_count[i] = si;
        n_pts_recv += si;
      }
    }
    if (copied_ranks != NULL) {
      free (copied_ranks);
    }

    if (n_pts_recv_copied_ranks != NULL) {
      free (n_pts_recv_copied_ranks);
    }

    for (int i = 0; i < n_rank; i++) {
      if (i == i_rank) {
        n_pts_local += send_count[i];
        send_count[i] = 0;
      }
      else if (i_copied_rank[i] >= 0) {
        send_count[i] -= copied_count[i_copied_rank[i]];
      }
    }

    copied_shift = PDM_array_new_idx_from_sizes_int (copied_count, n_copied_ranks);
    int *copied_count_tmp = PDM_array_zeros_int (n_copied_ranks);
    n_pts_copied = copied_shift[n_copied_ranks];

    /* Exchange new send/recv counts */
    PDM_MPI_Alltoall (send_count, 1, PDM_MPI_INT,
                      recv_count, 1, PDM_MPI_INT,
                      _dbbt->comm);

    send_shift = PDM_array_new_idx_from_sizes_int (send_count, n_rank);
    recv_shift = PDM_array_new_idx_from_sizes_int (recv_count, n_rank);
    PDM_array_reset_int (send_count, n_rank, 0);

    n_pts_recv = recv_shift[n_rank];
    n_pts1 = n_pts_local + n_pts_recv + n_pts_copied;
    if (idebug) {
      printf("[%d] n_pts1 = %d (without copies : %d)\n", i_rank, n_pts1, n_pts_recv_no_copies);
    }

    pts_g_num1 = malloc (sizeof(PDM_g_num_t) * n_pts1);
    pts_coord1 = malloc (sizeof(double)      * n_pts1 * 3);


    /* Fill send buffers */
    send_g_num = malloc (sizeof(PDM_g_num_t) * send_shift[n_rank]);
    send_coord = malloc (sizeof(double)      * send_shift[n_rank] * 3);

    int idx_copied = n_pts_local + n_pts_recv;
    PDM_g_num_t *copied_g_num = pts_g_num1 + idx_copied;
    double      *copied_coord = pts_coord1 + idx_copied * 3;
    n_pts_local = 0;
    for (int ipt = 0; ipt < n_pts; ipt++) {
      for (int i = pts_rank_idx[ipt]; i < pts_rank_idx[ipt+1]; i++) {
        int rank = pts_rank[i];

        if (rank == i_rank) {
          pts_g_num1[n_pts_local] = pts_g_num[ipt];
          for (int j = 0; j < 3; j++) {
            pts_coord1[3*n_pts_local + j] = _pts_coord[3*ipt + j];
          }
          n_pts_local++;
        }

        else if (i_copied_rank[rank] >= 0) {
          int _rank = i_copied_rank[rank];

          if (copied_count_tmp[_rank] < copied_count[_rank]) {
            int k = copied_shift[_rank] + copied_count_tmp[_rank];
            copied_g_num[k] = pts_g_num[ipt];
            for (int j = 0; j < 3; j++) {
              copied_coord[3*k + j] = _pts_coord[3*ipt + j];
            }
            copied_count_tmp[_rank]++;
          }
          else {
            int k = send_shift[rank] + send_count[rank];
            send_g_num[k] = pts_g_num[ipt];
            for (int j = 0; j < 3; j++) {
              send_coord[3*k + j] = _pts_coord[3*ipt + j];
            }
            send_count[rank]++;
          }
        }

        else {
          int k = send_shift[rank] + send_count[rank];
          send_g_num[k] = pts_g_num[ipt];
          for (int j = 0; j < 3; j++) {
            send_coord[3*k + j] = _pts_coord[3*ipt + j];
          }
          send_count[rank]++;
        }
      }
    }
    free (_pts_coord);
    if (copied_count != NULL) {
      free (copied_count);
    }
    if (copied_count_tmp != NULL) {
      free (copied_count_tmp);
    }
    if (i_copied_rank != NULL) {
      free (i_copied_rank);
    }
    free (pts_rank);
    free (pts_rank_idx);


    /* Exchange points */
    recv_g_num = pts_g_num1 + n_pts_local;
    PDM_MPI_Alltoallv (send_g_num, send_count, send_shift, PDM__PDM_MPI_G_NUM,
                       recv_g_num, recv_count, recv_shift, PDM__PDM_MPI_G_NUM,
                       _dbbt->comm);
    free (send_g_num);

    for (int i = 0; i < n_rank; i++) {
      send_count[i] *= 3;
      recv_count[i] *= 3;
      send_shift[i+1] *= 3;
      recv_shift[i+1] *= 3;
    }

    recv_coord = pts_coord1 + 3*n_pts_local;
    PDM_MPI_Alltoallv (send_coord, send_count, send_shift, PDM_MPI_DOUBLE,
                       recv_coord, recv_count, recv_shift, PDM_MPI_DOUBLE,
                       _dbbt->comm);

    free (send_coord);
    free (send_count);
    free (send_shift);
    free (recv_count);
    free (recv_shift);
  }

  else {
    n_pts_local  = n_pts;
    n_pts_recv   = 0;
    n_pts_copied = 0;

    n_pts1 = n_pts;

    pts_g_num1 = (PDM_g_num_t *) pts_g_num;
    pts_coord1 = _pts_coord;
  }



  /*
   *  Get boxes containing points
   */
  int n_part = 1 + n_copied_ranks;
  int *part_n_pts = malloc (sizeof(int) * n_part);
  part_n_pts[0] = n_pts_local + n_pts_recv;

  int **pts_box_idx   = malloc (sizeof(int *) * n_part);
  int **pts_box_l_num = malloc (sizeof(int *) * n_part);

  void (*_points_inside_volumes) (PDM_box_tree_t *bt,
                                  const int,
                                  const int,
                                  const double *,
                                  int **,
                                  int **);
  if (ellipsoids) {
    _points_inside_volumes = &PDM_box_tree_ellipsoids_containing_points;
  } else {
    _points_inside_volumes = &PDM_box_tree_boxes_containing_points;
  }

  /*
   *  Search in local tree
   */
  _points_inside_volumes (_dbbt->btLoc,
                          -1,
                          part_n_pts[0],
                          pts_coord1,
                          &(pts_box_idx[0]),
                          &(pts_box_l_num[0]));

  /*
   *  Search in copied trees
   */
  double      *pts_coord_copied = NULL;
  //PDM_g_num_t *pts_g_num_copied = NULL;
  if (n_copied_ranks > 0) {
    pts_coord_copied = pts_coord1 + part_n_pts[0] * 3;
    //pts_g_num_copied = pts_g_num1 + part_n_pts[0];
    for (int i = 0; i < n_copied_ranks; i++) {
      part_n_pts[i+1] = copied_shift[i+1] - copied_shift[i];

      _points_inside_volumes (_dbbt->btLoc,
                              i,
                              part_n_pts[i+1],
                              pts_coord_copied + copied_shift[i] * 3,
                              &(pts_box_idx[i+1]),
                              &(pts_box_l_num[i+1]));
    }
  }
  if (copied_shift != NULL) free (copied_shift);
  free (pts_coord1);

  PDM_g_num_t **pts_box_g_num = malloc (sizeof(PDM_g_num_t *) * n_part);
  pts_box_g_num[0] = malloc (sizeof(PDM_g_num_t) * pts_box_idx[0][part_n_pts[0]]);
  for (int j = 0; j < part_n_pts[0]; j++) {
    for (int k = pts_box_idx[0][j]; k < pts_box_idx[0][j+1]; k++) {
      pts_box_g_num[0][k] = _dbbt->boxes->local_boxes->g_num[pts_box_l_num[0][k]];
    }
  }

  for (int i = 0; i < n_copied_ranks; i++) {
    pts_box_g_num[i+1] = malloc (sizeof(PDM_g_num_t) * pts_box_idx[i+1][part_n_pts[i+1]]);
    for (int j = 0; j < part_n_pts[i+1]; j++) {
      for (int k = pts_box_idx[i+1][j]; k < pts_box_idx[i+1][j+1]; k++) {
        pts_box_g_num[i+1][k] = _dbbt->boxes->rank_boxes[i].g_num[pts_box_l_num[i+1][k]];
      }
    }
  }


  PDM_g_num_t **part_pts_g_num = malloc (sizeof(PDM_g_num_t *) * n_part);
  int    **part_stride = malloc (sizeof(int *)    * n_part);
  double **part_weight = malloc (sizeof(double *) * n_part);
  int idx = 0;
  for (int ipart = 0; ipart < n_part; ipart++) {
    part_pts_g_num[ipart] = pts_g_num1 + idx;
    idx += part_n_pts[ipart];

    part_stride[ipart] = malloc (sizeof(int)    * part_n_pts[ipart]);
    part_weight[ipart] = malloc (sizeof(double) * part_n_pts[ipart]);
    for (int i = 0; i < part_n_pts[ipart]; i++) {
      part_stride[ipart][i] = pts_box_idx[ipart][i+1] - pts_box_idx[ipart][i];
      part_weight[ipart][i] = (double) part_stride[ipart][i];
    }
    free (pts_box_idx[ipart]);
    free (pts_box_l_num[ipart]);
  }
  free (pts_box_idx);
  free (pts_box_l_num);


  if (0) {
    printf("[%d] n_pts1 = %d\n", i_rank, n_pts1);
    printf("[%d] before merge:\n", i_rank);
    for (int ipart = 0; ipart < n_part; ipart++) {
      for (int i = 0; i< part_n_pts[ipart]; i++) {
        printf("[%d] pt "PDM_FMT_G_NUM" inside boxes", i_rank, part_pts_g_num[ipart][i]);
        for (int j = pts_box_idx[ipart][i]; j < pts_box_idx[ipart][i+1]; j++) {
          printf(" "PDM_FMT_G_NUM, pts_box_g_num[ipart][j]);
        }
        printf("\n");
      }
    }
  }


  /*
   *  Merge results
   */
  /* 1) Part-to-Block */
  PDM_part_to_block_t *ptb = PDM_part_to_block_create (PDM_PART_TO_BLOCK_DISTRIB_ALL_PROC,
                                                       PDM_PART_TO_BLOCK_POST_MERGE,
                                                       1.,
                                                       part_pts_g_num,
                                                       part_weight,
                                                       part_n_pts,
                                                       n_part,
                                                       _dbbt->comm);

  int *block_box_n = NULL;
  PDM_g_num_t *block_box_g_num = NULL;
  PDM_part_to_block_exch (ptb,
                          sizeof(PDM_g_num_t),
                          PDM_STRIDE_VAR_INTERLACED,
                          1,
                          part_stride,
                          (void **) pts_box_g_num,
                          &block_box_n,
                          (void **) &block_box_g_num);

  for (int ipart = 0; ipart < n_part; ipart++) {
    free (part_stride[ipart]);
    free (part_weight[ipart]);
    free (pts_box_g_num[ipart]);
  }
  free (part_stride);
  free (part_weight);
  free (pts_box_g_num);


  /* Remove doubles */
  int idx1 = 0, idx2 = 0;
  int n_pts_block = PDM_part_to_block_n_elt_block_get (ptb);
  int max_n = 0;
  for (int i = 0; i < n_pts_block; i++) {
    max_n = PDM_MAX (max_n, block_box_n[i]);
  }

  int *order = malloc (sizeof(int) * max_n);
  idx1 = 0;
  idx2 = 0;
  for (int i = 0; i < n_pts_block; i++) {
    if (block_box_n[i] == 0) continue;

    PDM_g_num_t *_g_num1 = block_box_g_num + idx1;
    PDM_g_num_t *_g_num2 = block_box_g_num + idx2;

    for (int j = 0; j < block_box_n[i]; j++) {
      order[j] = j;
    }
    PDM_sort_long (_g_num1,
                   order,
                   block_box_n[i]);

    _g_num2[0] = _g_num1[0];
    int tmp_n = 1;
    for (int j = 1; j < block_box_n[i]; j++) {
      if (_g_num1[j] != _g_num2[tmp_n-1]) {
        _g_num2[tmp_n++] = _g_num1[j];
      }
    }

    idx1 += block_box_n[i];
    idx2 += tmp_n;
    block_box_n[i] = tmp_n;
  }
  free (order);

  /* Fix partial block stride */
  PDM_g_num_t l_max_g_num = 0;
  for (int i = 0; i < n_pts; i++) {
    l_max_g_num = PDM_MAX (l_max_g_num, pts_g_num[i]);
  }

  PDM_g_num_t g_max_g_num;
  PDM_MPI_Allreduce (&l_max_g_num, &g_max_g_num, 1,
                     PDM__PDM_MPI_G_NUM, PDM_MPI_MAX, _dbbt->comm);

  PDM_g_num_t *block_distrib_idx =
    PDM_part_to_block_adapt_partial_block_to_block (ptb,
                                                    &block_box_n,
                                                    g_max_g_num);

  /* 2) Block-to-Part */
  PDM_block_to_part_t *btp = PDM_block_to_part_create (block_distrib_idx,
                                                       (const PDM_g_num_t **) &pts_g_num,
                                                       &n_pts,
                                                       1,
                                                       _dbbt->comm);

  int *box_n = malloc (sizeof(int) * n_pts);
  int one = 1;
  PDM_block_to_part_exch_in_place (btp,
                          sizeof(int),
                          PDM_STRIDE_CST_INTERLACED,
                          &one,
                          block_box_n,
                          NULL,
                          (void **) &box_n);

  *box_idx = PDM_array_new_idx_from_sizes_int (box_n, n_pts);

  *box_g_num = malloc (sizeof(PDM_g_num_t) * (*box_idx)[n_pts]);
  PDM_block_to_part_exch_in_place (btp,
                          sizeof(PDM_g_num_t),
                          PDM_STRIDE_VAR_INTERLACED,
                          block_box_n,
                          block_box_g_num,
                          &box_n,
                          (void **) box_g_num);
  free (block_box_g_num);
  free (block_box_n);
  free (box_n);

  btp = PDM_block_to_part_free (btp);
  ptb = PDM_part_to_block_free (ptb);
  free (block_distrib_idx);
  free (part_n_pts);

  free (part_pts_g_num);
  if (pts_g_num1 != pts_g_num) free (pts_g_num1);

  PDM_box_tree_free_copies(_dbbt->btLoc);
}












static void
_lines_intersect_shared_box_tree
(
 PDM_dbbtree_t   *dbbt,
 const int        n_line,
 PDM_g_num_t     *line_g_num,
 double          *line_coord,
 float            f_threshold,
 float            f_max_copy,
 int            **redistrib_n_line,
 PDM_g_num_t   ***redistrib_line_g_num,
 double        ***redistrib_line_coord,
 int             *n_copied_ranks
 )
{
  int idebug = 0;

  assert (dbbt != NULL);
  _PDM_dbbtree_t *_dbbt = (_PDM_dbbtree_t *) dbbt;

  int i_rank, n_rank;
  PDM_MPI_Comm_rank (_dbbt->comm, &i_rank);
  PDM_MPI_Comm_size (_dbbt->comm, &n_rank);

  /*
   *  Normalize coordinates
   */
  double *_line_coord = malloc (sizeof(double) * n_line * 6);
  for (int i = 0; i < 2*n_line; i++) {
    _normalize (_dbbt,
                line_coord  + 3*i,
                _line_coord + 3*i);
  }


  /*
   *  For each line, find all ranks that might have boxes intersecting that line
   */
  int *send_count = NULL;
  int *send_shift = NULL;
  int *recv_count = NULL;
  int *recv_shift = NULL;
  PDM_g_num_t *send_g_num = NULL;
  PDM_g_num_t *recv_g_num = NULL;
  double      *send_coord = NULL;
  double      *recv_coord = NULL;

  *n_copied_ranks = 0;
  int *copied_ranks = NULL;
  int n_line_local  = 0;
  int n_line_recv   = 0;
  int n_line_copied = 0;
  int n_line1;

  int *copied_shift = NULL;

  // PDM_g_num_t *line_g_num1 = NULL;
  // double      *line_coord1 = NULL;

  if (_dbbt->btShared != NULL) {
    int *line_rank_idx = NULL;
    int *line_rank     = NULL;
    PDM_box_tree_intersect_lines_boxes (_dbbt->btShared,
                                        -1,
                                        n_line,
                                        line_coord,
                                        &line_rank_idx,
                                        &line_rank);

    /* Count points to send to each rank */
    send_count = PDM_array_zeros_int (n_rank);

    for (int i = 0; i < line_rank_idx[n_line]; i++) {
      int rank = _dbbt->usedRank[line_rank[i]];
      line_rank[i] = rank;
      send_count[rank]++;
    }

    // PDM_log_trace_array_int(send_count, _dbbt->nUsedRank, "send_count ::");

    recv_count = malloc (sizeof(int) * n_rank);
    PDM_MPI_Alltoall (send_count, 1, PDM_MPI_INT,
                      recv_count, 1, PDM_MPI_INT,
                      _dbbt->comm);

    n_line_recv = 0;
    for (int i = 0; i < n_rank; i++) {
      n_line_recv += recv_count[i];
    }
    int n_line_recv_no_copies = n_line_recv;

    /* Prepare copies */
    int n_max_copy = (int) (f_max_copy * n_rank);
    int *i_copied_rank = NULL;
    int *copied_count = NULL;
    int mean_n_line_recv = 0;
    int *n_line_recv_copied_ranks = NULL;

    if (n_max_copy > 0) {
      int *all_n_line_recv = malloc (sizeof(int) * n_rank);
      PDM_MPI_Allgather (&n_line_recv,    1, PDM_MPI_INT,
                         all_n_line_recv, 1, PDM_MPI_INT,
                         _dbbt->comm);

      // Mean number of line_recvs
      long l_mean_n_line_recv = 0;
      for (int i = 0; i < n_rank; i++) {
        l_mean_n_line_recv += all_n_line_recv[i];
      }
      mean_n_line_recv = (int) (l_mean_n_line_recv / n_rank);

      float n_threshold = f_threshold * mean_n_line_recv;

      // Sort ranks
      int *order = malloc (sizeof(int) * n_rank);
      for (int i = 0; i < n_rank; i++) {
        order[i] = i;
      }

      PDM_sort_int (all_n_line_recv,
                    order,
                    n_rank);

      // Identify ranks to copy
      copied_ranks = malloc (sizeof(int) * n_max_copy);
      n_line_recv_copied_ranks = malloc (sizeof(int) * n_max_copy);
      for (int i = 0; i < n_max_copy; i++) {
        int j = n_rank - i - 1;

        if (all_n_line_recv[j] > n_threshold) {
          copied_ranks[*n_copied_ranks] = order[j];
          n_line_recv_copied_ranks[*n_copied_ranks] = all_n_line_recv[j];
          (*n_copied_ranks)++;
        }
        else {
          break;
        }
      }
      free (all_n_line_recv);
      free (order);

      if (*n_copied_ranks > 0) {
        copied_ranks = realloc (copied_ranks, sizeof(int) * (*n_copied_ranks));
        n_line_recv_copied_ranks = realloc (n_line_recv_copied_ranks,
                                           sizeof(int) * (*n_copied_ranks));
      }
    }

    if (idebug && i_rank == 0) {
      if (*n_copied_ranks > 0) {
        if (*n_copied_ranks == 1) {
          printf("1 copied rank: %d\n", copied_ranks[0]);
        }
        else {
          printf("%d copied ranks:", *n_copied_ranks);
          for (int i = 0; i < *n_copied_ranks; i++) {
            printf(" %d", copied_ranks[i]);
          }
          printf("\n");
        }
      }
      else {
        printf("0 copied ranks\n");
      }
    }

    /*
     * Copy the data of selected ranks
     */
    i_copied_rank = malloc (sizeof(int) * n_rank);
    PDM_box_tree_copy_to_ranks (_dbbt->btLoc,
                                n_copied_ranks,
                                copied_ranks,
                                i_copied_rank);

    /* Re-compute send/recv counts */
    copied_count = PDM_array_zeros_int (*n_copied_ranks);
    n_line_local = 0;
    for (int i = 0; i < *n_copied_ranks; i++) {
      int rank = copied_ranks[i];
      if (rank != i_rank) {
        int si = send_count[rank];

        si = PDM_MIN (si, PDM_MAX (0, (n_line_recv_copied_ranks[i] - n_line_recv)/2));
        if (i_copied_rank[i_rank] < 0) {
          si = PDM_MIN (si, PDM_MAX (0, mean_n_line_recv - n_line_recv));
        }

        copied_count[i] = si;
        n_line_recv += si;
      }
    }
    if (copied_ranks != NULL) {
      free (copied_ranks);
    }

    if (n_line_recv_copied_ranks != NULL) {
      free (n_line_recv_copied_ranks);
    }

    for (int i = 0; i < n_rank; i++) {
      if (i == i_rank) {
        n_line_local += send_count[i];
        send_count[i] = 0;
      }
      else if (i_copied_rank[i] >= 0) {
        send_count[i] -= copied_count[i_copied_rank[i]];
      }
    }

    copied_shift = PDM_array_new_idx_from_sizes_int (copied_count, *n_copied_ranks);
    int *copied_count_tmp = PDM_array_zeros_int (*n_copied_ranks);
    n_line_copied = copied_shift[*n_copied_ranks];

    /* Exchange new send/recv counts */
    PDM_MPI_Alltoall (send_count, 1, PDM_MPI_INT,
                      recv_count, 1, PDM_MPI_INT,
                      _dbbt->comm);

    send_shift = PDM_array_new_idx_from_sizes_int (send_count, n_rank);
    recv_shift = PDM_array_new_idx_from_sizes_int (recv_count, n_rank);
    PDM_array_reset_int (send_count, n_rank, 0);

    n_line_recv = recv_shift[n_rank];
    n_line1 = n_line_local + n_line_recv + n_line_copied;
    if (idebug) {
      printf("[%d] n_line1 = %d (without copies : %d)\n", i_rank, n_line1, n_line_recv_no_copies);
    }

    // line_g_num1 = malloc (sizeof(PDM_g_num_t) * n_line1);
    // line_coord1 = malloc (sizeof(double)      * n_line1 * 6);

    int n_part = 1 + (*n_copied_ranks);
    *redistrib_n_line     = (int          *) malloc(sizeof(int          ) * n_part);
    *redistrib_line_g_num = (PDM_g_num_t **) malloc(sizeof(PDM_g_num_t *) * n_part);
    *redistrib_line_coord = (double      **) malloc(sizeof(double      *) * n_part);

    (*redistrib_n_line)    [0] = n_line_local + n_line_recv;
    (*redistrib_line_g_num)[0] = (PDM_g_num_t *) malloc(sizeof(PDM_g_num_t) * (n_line_local + n_line_recv));
    (*redistrib_line_coord)[0] = (double      *) malloc(sizeof(double     ) * (n_line_local + n_line_recv)*6);
    for (int i = 0; i < (*n_copied_ranks); i++) {
      int n = copied_shift[i+1] - copied_shift[i];
      (*redistrib_n_line)    [i+1] = n;
      (*redistrib_line_g_num)[i+1] = (PDM_g_num_t *) malloc(sizeof(PDM_g_num_t) * n);
      (*redistrib_line_coord)[i+1] = (double      *) malloc(sizeof(double     ) * n*6);
    }


    /* Fill send buffers */
    send_g_num = malloc (sizeof(PDM_g_num_t) * send_shift[n_rank]);
    send_coord = malloc (sizeof(double)      * send_shift[n_rank] * 6);

    // int idx_copied = n_line_local + n_line_recv;
    // PDM_g_num_t *copied_g_num = line_g_num1 + idx_copied;
    // double      *copied_coord = line_coord1 + idx_copied * 6;
    n_line_local = 0;
    for (int iline = 0; iline < n_line; iline++) {
      for (int i = line_rank_idx[iline]; i < line_rank_idx[iline+1]; i++) {
        int rank = line_rank[i];

        if (rank == i_rank) {
          // line_g_num1[n_line_local] = line_g_num[iline];
          (*redistrib_line_g_num)[0][n_line_local] = line_g_num[iline];
          for (int j = 0; j < 6; j++) {
            // line_coord1[6*n_line_local + j] = _line_coord[6*iline + j];
            (*redistrib_line_coord)[0][6*n_line_local + j] = _line_coord[6*iline + j];
          }
          n_line_local++;
        }

        else if (i_copied_rank[rank] >= 0) {
          int _rank = i_copied_rank[rank];

          if (copied_count_tmp[_rank] < copied_count[_rank]) {
            // int k = copied_shift[_rank] + copied_count_tmp[_rank];
            // copied_g_num[k] = line_g_num[iline];
            int k = copied_count_tmp[_rank];
            (*redistrib_line_g_num)[_rank+1][k] = line_g_num[iline];

            for (int j = 0; j < 6; j++) {
              // copied_coord[6*k + j] = _line_coord[6*iline + j];
              (*redistrib_line_coord)[_rank+1][6*k + j] = _line_coord[6*iline + j];
            }
            copied_count_tmp[_rank]++;
          }
          else {
            int k = send_shift[rank] + send_count[rank];
            send_g_num[k] = line_g_num[iline];
            for (int j = 0; j < 6; j++) {
              send_coord[6*k + j] = _line_coord[6*iline + j];
            }
            send_count[rank]++;
          }
        }

        else {
          int k = send_shift[rank] + send_count[rank];
          send_g_num[k] = line_g_num[iline];
          for (int j = 0; j < 6; j++) {
            send_coord[6*k + j] = _line_coord[6*iline + j];
          }
          send_count[rank]++;
        }
      }
    }
    free (_line_coord);
    if (copied_count != NULL) {
      free (copied_count);
    }
    if (copied_count_tmp != NULL) {
      free (copied_count_tmp);
    }
    if (i_copied_rank != NULL) {
      free (i_copied_rank);
    }
    if (copied_shift != NULL) {
      free(copied_shift);
    }
    free (line_rank);
    free (line_rank_idx);


    /* Exchange points */
    // recv_g_num = line_g_num1 + n_line_local;
    // recv_g_num = (PDM_g_num_t *) malloc(sizeof(PDM_g_num_t) * recv_shift[n_rank]);
    recv_g_num = (*redistrib_line_g_num)[0] + n_line_local;
    PDM_MPI_Alltoallv (send_g_num, send_count, send_shift, PDM__PDM_MPI_G_NUM,
                       recv_g_num, recv_count, recv_shift, PDM__PDM_MPI_G_NUM,
                       _dbbt->comm);

    free (send_g_num);

    for (int i = 0; i < n_rank; i++) {
      send_count[i]   *= 6;
      recv_count[i]   *= 6;
      send_shift[i+1] *= 6;
      recv_shift[i+1] *= 6;
    }

    // recv_coord = line_coord1 + 6*n_line_local;
    // recv_coord = (double *) malloc(sizeof(double) * recv_shift[n_rank]);
    recv_coord = (*redistrib_line_coord)[0] + 6*n_line_local;
    PDM_MPI_Alltoallv (send_coord, send_count, send_shift, PDM_MPI_DOUBLE,
                       recv_coord, recv_count, recv_shift, PDM_MPI_DOUBLE,
                       _dbbt->comm);

    free (send_coord);
    free (send_count);
    free (send_shift);
    free (recv_count);
    free (recv_shift);
  }

  else {
    // n_line_local  = n_line;
    // n_line_recv   = 0;
    // n_line_copied = 0;

    // n_line1 = n_line;

    // line_g_num1 = (PDM_g_num_t *) line_g_num;
    // line_coord1 = _line_coord;
    *redistrib_n_line = (int *) malloc(sizeof(int) * 1);
    (*redistrib_n_line)[0] = n_line;

    *redistrib_line_g_num = (PDM_g_num_t **) malloc(sizeof(PDM_g_num_t *) * 1);
    // (*redistrib_line_g_num)[0] = line_g_num;
    (*redistrib_line_g_num)[0] = malloc(sizeof(PDM_g_num_t) * n_line);
    memcpy((*redistrib_line_g_num)[0], line_g_num, sizeof(PDM_g_num_t) * n_line);

    *redistrib_line_coord = (double **) malloc(sizeof(double *) * 1);
    (*redistrib_line_coord)[0] = _line_coord;
    // (*redistrib_line_coord)[0] = malloc(sizeof(double) * n_line * 6);
    // memcpy((*redistrib_line_coord)[0], _line_coord, sizeof(double) * n_line * 6);

    // free (_line_coord);
  }

  // *redistrib_n_line_local = n_line_local + n_line_recv;
  // *redistrib_n_line_copy  = n_line_copied;
  // *redistrib_line_g_num   = line_g_num1;
  // *redistrib_line_coord   = line_coord1;
}






/**
 *
 * \brief Get an indexed list of all boxes intersecting lines
 *
 *   \param [in]  dbbt        Pointer to distributed box tree structure
 *   \param [in]  n_line      Number of points
 *   \param [in]  line_g_num  Line global ids (size = \ref n_line)
 *   \param [in]  line_coord  Line coordinates (size = 6 * \ref n_line)
 *   \param [out] box_idx     Index of boxes (size = \ref n_line + 1, allocated inside function)
 *   \param [out] box_g_num   Global ids of boxes (size = \ref box_line_idx[\ref n_line], allocated inside function)
 *
 */

void
PDM_dbbtree_lines_intersect_boxes
(
 PDM_dbbtree_t  *dbbt,
 const int       n_line,
 PDM_g_num_t    *line_g_num,
 double         *line_coord,
 int           **box_idx,
 PDM_g_num_t   **box_g_num
 )
{
  int idebug = 0;

  const float f_threshold = 1.1;  // factor of the mean nb of requests
  const float f_max_copy  = 0.1;  // factor of the total nb of processes

  assert (dbbt != NULL);
  _PDM_dbbtree_t *_dbbt = (_PDM_dbbtree_t *) dbbt;

  int i_rank, n_rank;
  PDM_MPI_Comm_rank (_dbbt->comm, &i_rank);
  PDM_MPI_Comm_size (_dbbt->comm, &n_rank);

  /*
   *  Normalize coordinates
   */
  double *_line_coord = malloc (sizeof(double) * n_line * 6);
  for (int i = 0; i < 2*n_line; i++) {
    _normalize (_dbbt,
                line_coord  + 3*i,
                _line_coord + 3*i);
  }



  /*
   *  For each line, find all ranks that might have boxes intersecting that line
   */
  int *send_count = NULL;
  int *send_shift = NULL;
  int *recv_count = NULL;
  int *recv_shift = NULL;
  PDM_g_num_t *send_g_num = NULL;
  PDM_g_num_t *recv_g_num = NULL;
  double      *send_coord = NULL;
  double      *recv_coord = NULL;

  int n_copied_ranks = 0;
  int *copied_ranks = NULL;
  int n_line_local  = 0;
  int n_line_recv   = 0;
  int n_line_copied = 0;
  int n_line1;

  int *copied_shift = NULL;

  PDM_g_num_t *line_g_num1 = NULL;
  double      *line_coord1 = NULL;

  if (_dbbt->btShared != NULL) {
    int *line_rank_idx = NULL;
    int *line_rank     = NULL;
    PDM_box_tree_intersect_lines_boxes (_dbbt->btShared,
                                        -1,
                                        n_line,
                                        line_coord,
                                        &line_rank_idx,
                                        &line_rank);

    /* Count points to send to each rank */
    send_count = PDM_array_zeros_int (n_rank);

    for (int i = 0; i < line_rank_idx[n_line]; i++) {
      int rank = _dbbt->usedRank[line_rank[i]];
      line_rank[i] = rank;
      send_count[rank]++;
    }

    recv_count = malloc (sizeof(int) * n_rank);
    PDM_MPI_Alltoall (send_count, 1, PDM_MPI_INT,
                      recv_count, 1, PDM_MPI_INT,
                      _dbbt->comm);

    n_line_recv = 0;
    for (int i = 0; i < n_rank; i++) {
      n_line_recv += recv_count[i];
    }
    int n_line_recv_no_copies = n_line_recv;

    /* Prepare copies */
    int n_max_copy = (int) (f_max_copy * n_rank);
    int *i_copied_rank = NULL;
    int *copied_count = NULL;
    int mean_n_line_recv = 0;
    int *n_line_recv_copied_ranks = NULL;

    if (n_max_copy > 0) {
      int *all_n_line_recv = malloc (sizeof(int) * n_rank);
      PDM_MPI_Allgather (&n_line_recv,    1, PDM_MPI_INT,
                         all_n_line_recv, 1, PDM_MPI_INT,
                         _dbbt->comm);

      // Mean number of line_recvs
      long l_mean_n_line_recv = 0;
      for (int i = 0; i < n_rank; i++) {
        l_mean_n_line_recv += all_n_line_recv[i];
      }
      mean_n_line_recv = (int) (l_mean_n_line_recv / n_rank);

      float n_threshold = f_threshold * mean_n_line_recv;

      // Sort ranks
      int *order = malloc (sizeof(int) * n_rank);
      for (int i = 0; i < n_rank; i++) {
        order[i] = i;
      }

      PDM_sort_int (all_n_line_recv,
                    order,
                    n_rank);

      // Identify ranks to copy
      copied_ranks = malloc (sizeof(int) * n_max_copy);
      n_line_recv_copied_ranks = malloc (sizeof(int) * n_max_copy);
      for (int i = 0; i < n_max_copy; i++) {
        int j = n_rank - i - 1;

        if (all_n_line_recv[j] > n_threshold) {
          copied_ranks[n_copied_ranks] = order[j];
          n_line_recv_copied_ranks[n_copied_ranks] = all_n_line_recv[j];
          n_copied_ranks++;
        }
        else {
          break;
        }
      }
      free (all_n_line_recv);
      free (order);

      if (n_copied_ranks > 0) {
        copied_ranks = realloc (copied_ranks, sizeof(int) * n_copied_ranks);
        n_line_recv_copied_ranks = realloc (n_line_recv_copied_ranks,
                                           sizeof(int) * n_copied_ranks);
      }
    }

    if (idebug && i_rank == 0) {
      if (n_copied_ranks > 0) {
        if (n_copied_ranks == 1) {
          printf("1 copied rank: %d\n", copied_ranks[0]);
        }
        else {
          printf("%d copied ranks:", n_copied_ranks);
          for (int i = 0; i < n_copied_ranks; i++) {
            printf(" %d", copied_ranks[i]);
          }
          printf("\n");
        }
      }
      else {
        printf("0 copied ranks\n");
      }
    }

    /*
     * Copy the data of selected ranks
     */
    i_copied_rank = malloc (sizeof(int) * n_rank);
    PDM_box_tree_copy_to_ranks (_dbbt->btLoc,
                                &n_copied_ranks,
                                copied_ranks,
                                i_copied_rank);

    /* Re-compute send/recv counts */
    copied_count = PDM_array_zeros_int (n_copied_ranks);
    n_line_local = 0;
    for (int i = 0; i < n_copied_ranks; i++) {
      int rank = copied_ranks[i];
      if (rank != i_rank) {
        int si = send_count[rank];

        si = PDM_MIN (si, PDM_MAX (0, (n_line_recv_copied_ranks[i] - n_line_recv)/2));
        if (i_copied_rank[i_rank] < 0) {
          si = PDM_MIN (si, PDM_MAX (0, mean_n_line_recv - n_line_recv));
        }

        copied_count[i] = si;
        n_line_recv += si;
      }
    }
    if (copied_ranks != NULL) {
      free (copied_ranks);
    }

    if (n_line_recv_copied_ranks != NULL) {
      free (n_line_recv_copied_ranks);
    }

    for (int i = 0; i < n_rank; i++) {
      if (i == i_rank) {
        n_line_local += send_count[i];
        send_count[i] = 0;
      }
      else if (i_copied_rank[i] >= 0) {
        send_count[i] -= copied_count[i_copied_rank[i]];
      }
    }

    copied_shift = PDM_array_new_idx_from_sizes_int (copied_count, n_copied_ranks);
    int *copied_count_tmp = PDM_array_zeros_int (n_copied_ranks);
    n_line_copied = copied_shift[n_copied_ranks];

    /* Exchange new send/recv counts */
    PDM_MPI_Alltoall (send_count, 1, PDM_MPI_INT,
                      recv_count, 1, PDM_MPI_INT,
                      _dbbt->comm);

    send_shift = PDM_array_new_idx_from_sizes_int (send_count, n_rank);
    recv_shift = PDM_array_new_idx_from_sizes_int (recv_count, n_rank);
    PDM_array_reset_int (send_count, n_rank, 0);

    n_line_recv = recv_shift[n_rank];
    n_line1 = n_line_local + n_line_recv + n_line_copied;
    if (idebug) {
      printf("[%d] n_line1 = %d (without copies : %d)\n", i_rank, n_line1, n_line_recv_no_copies);
    }

    line_g_num1 = malloc (sizeof(PDM_g_num_t) * n_line1);
    line_coord1 = malloc (sizeof(double)      * n_line1 * 6);


    /* Fill send buffers */
    send_g_num = malloc (sizeof(PDM_g_num_t) * send_shift[n_rank]);
    send_coord = malloc (sizeof(double)      * send_shift[n_rank] * 6);

    int idx_copied = n_line_local + n_line_recv;
    PDM_g_num_t *copied_g_num = line_g_num1 + idx_copied;
    double      *copied_coord = line_coord1 + idx_copied * 6;
    n_line_local = 0;
    for (int ipt = 0; ipt < n_line; ipt++) {
      for (int i = line_rank_idx[ipt]; i < line_rank_idx[ipt+1]; i++) {
        int rank = line_rank[i];

        if (rank == i_rank) {
          line_g_num1[n_line_local] = line_g_num[ipt];
          for (int j = 0; j < 6; j++) {
            line_coord1[6*n_line_local + j] = _line_coord[6*ipt + j];
          }
          n_line_local++;
        }

        else if (i_copied_rank[rank] >= 0) {
          int _rank = i_copied_rank[rank];

          if (copied_count_tmp[_rank] < copied_count[_rank]) {
            int k = copied_shift[_rank] + copied_count_tmp[_rank];
            copied_g_num[k] = line_g_num[ipt];
            for (int j = 0; j < 6; j++) {
              copied_coord[6*k + j] = _line_coord[6*ipt + j];
            }
            copied_count_tmp[_rank]++;
          }
          else {
            int k = send_shift[rank] + send_count[rank];
            send_g_num[k] = line_g_num[ipt];
            for (int j = 0; j < 6; j++) {
              send_coord[6*k + j] = _line_coord[6*ipt + j];
            }
            send_count[rank]++;
          }
        }

        else {
          int k = send_shift[rank] + send_count[rank];
          send_g_num[k] = line_g_num[ipt];
          for (int j = 0; j < 6; j++) {
            send_coord[6*k + j] = _line_coord[6*ipt + j];
          }
          send_count[rank]++;
        }
      }
    }
    free (_line_coord);
    if (copied_count != NULL) {
      free (copied_count);
    }
    if (copied_count_tmp != NULL) {
      free (copied_count_tmp);
    }
    if (i_copied_rank != NULL) {
      free (i_copied_rank);
    }
    free (line_rank);
    free (line_rank_idx);


    /* Exchange points */
    recv_g_num = line_g_num1 + n_line_local;
    PDM_MPI_Alltoallv (send_g_num, send_count, send_shift, PDM__PDM_MPI_G_NUM,
                       recv_g_num, recv_count, recv_shift, PDM__PDM_MPI_G_NUM,
                       _dbbt->comm);
    free (send_g_num);

    for (int i = 0; i < n_rank; i++) {
      send_count[i]   *= 6;
      recv_count[i]   *= 6;
      send_shift[i+1] *= 6;
      recv_shift[i+1] *= 6;
    }

    recv_coord = line_coord1 + 6*n_line_local;
    PDM_MPI_Alltoallv (send_coord, send_count, send_shift, PDM_MPI_DOUBLE,
                       recv_coord, recv_count, recv_shift, PDM_MPI_DOUBLE,
                       _dbbt->comm);

    free (send_coord);
    free (send_count);
    free (send_shift);
    free (recv_count);
    free (recv_shift);
  }

  else {
    n_line_local  = n_line;
    n_line_recv   = 0;
    n_line_copied = 0;

    n_line1 = n_line;

    line_g_num1 = (PDM_g_num_t *) line_g_num;
    line_coord1 = _line_coord;
  }



  /*
   *  Get boxes intersecting lines
   */
  int n_part = 1 + n_copied_ranks;
  int *part_n_line = malloc (sizeof(int) * n_part);
  part_n_line[0] = n_line_local + n_line_recv;

  int **line_box_idx   = malloc (sizeof(int *) * n_part);
  int **line_box_l_num = malloc (sizeof(int *) * n_part);

  /*
   *  Search in local tree
   */
  PDM_box_tree_intersect_lines_boxes (_dbbt->btLoc,
                                      -1,
                                      part_n_line[0],
                                      line_coord1,
                                      &(line_box_idx[0]),
                                      &(line_box_l_num[0]));

  /*
   *  Search in copied trees
   */
  double *line_coord_copied = NULL;
  if (n_copied_ranks > 0) {
    line_coord_copied = line_coord1 + part_n_line[0] * 6;
    for (int i = 0; i < n_copied_ranks; i++) {
      part_n_line[i+1] = copied_shift[i+1] - copied_shift[i];

      PDM_box_tree_intersect_lines_boxes (_dbbt->btLoc,
                                          i,
                                          part_n_line[i+1],
                                          line_coord_copied + copied_shift[i] * 6,
                                          &(line_box_idx[i+1]),
                                          &(line_box_l_num[i+1]));
    }
  }
  if (copied_shift != NULL) free (copied_shift);
  free (line_coord1);

  PDM_g_num_t **line_box_g_num = malloc (sizeof(PDM_g_num_t *) * n_part);
  line_box_g_num[0] = malloc (sizeof(PDM_g_num_t) * line_box_idx[0][part_n_line[0]]);
  for (int j = 0; j < part_n_line[0]; j++) {
    for (int k = line_box_idx[0][j]; k < line_box_idx[0][j+1]; k++) {
      line_box_g_num[0][k] = _dbbt->boxes->local_boxes->g_num[line_box_l_num[0][k]];
    }
  }

  for (int i = 0; i < n_copied_ranks; i++) {
    line_box_g_num[i+1] = malloc (sizeof(PDM_g_num_t) * line_box_idx[i+1][part_n_line[i+1]]);
    for (int j = 0; j < part_n_line[i+1]; j++) {
      for (int k = line_box_idx[i+1][j]; k < line_box_idx[i+1][j+1]; k++) {
        line_box_g_num[i+1][k] = _dbbt->boxes->rank_boxes[i].g_num[line_box_l_num[i+1][k]];
      }
    }
  }


  PDM_g_num_t **part_line_g_num = malloc (sizeof(PDM_g_num_t *) * n_part);
  int    **part_stride = malloc (sizeof(int *)    * n_part);
  double **part_weight = malloc (sizeof(double *) * n_part);
  int idx = 0;
  for (int ipart = 0; ipart < n_part; ipart++) {
    part_line_g_num[ipart] = line_g_num1 + idx;
    idx += part_n_line[ipart];

    part_stride[ipart] = malloc (sizeof(int)    * part_n_line[ipart]);
    part_weight[ipart] = malloc (sizeof(double) * part_n_line[ipart]);
    for (int i = 0; i < part_n_line[ipart]; i++) {
      part_stride[ipart][i] = line_box_idx[ipart][i+1] - line_box_idx[ipart][i];
      part_weight[ipart][i] = (double) part_stride[ipart][i];
    }
    free (line_box_idx[ipart]);
    free (line_box_l_num[ipart]);
  }
  free (line_box_idx);
  free (line_box_l_num);


  /*
   *  Merge results
   */
  /* 1) Part-to-Block */
  PDM_part_to_block_t *ptb = PDM_part_to_block_create (PDM_PART_TO_BLOCK_DISTRIB_ALL_PROC,
                                                       PDM_PART_TO_BLOCK_POST_MERGE,
                                                       1.,
                                                       part_line_g_num,
                                                       part_weight,
                                                       part_n_line,
                                                       n_part,
                                                       _dbbt->comm);

  int *block_box_n = NULL;
  PDM_g_num_t *block_box_g_num = NULL;
  PDM_part_to_block_exch (ptb,
                          sizeof(PDM_g_num_t),
                          PDM_STRIDE_VAR_INTERLACED,
                          1,
                          part_stride,
                          (void **) line_box_g_num,
                          &block_box_n,
                          (void **) &block_box_g_num);

  for (int ipart = 0; ipart < n_part; ipart++) {
    free (part_stride[ipart]);
    free (part_weight[ipart]);
    free (line_box_g_num[ipart]);
  }
  free (part_stride);
  free (part_weight);
  free (line_box_g_num);


  /* Remove doubles */
  int idx1 = 0, idx2 = 0;
  int n_line_block = PDM_part_to_block_n_elt_block_get (ptb);
  int max_n = 0;
  for (int i = 0; i < n_line_block; i++) {
    max_n = PDM_MAX (max_n, block_box_n[i]);
  }

  int *order = malloc (sizeof(int) * max_n);
  idx1 = 0;
  idx2 = 0;
  for (int i = 0; i < n_line_block; i++) {
    if (block_box_n[i] == 0) continue;

    PDM_g_num_t *_g_num1 = block_box_g_num + idx1;
    PDM_g_num_t *_g_num2 = block_box_g_num + idx2;

    for (int j = 0; j < block_box_n[i]; j++) {
      order[j] = j;
    }
    PDM_sort_long (_g_num1,
                   order,
                   block_box_n[i]);

    _g_num2[0] = _g_num1[0];
    int tmp_n = 1;
    for (int j = 1; j < block_box_n[i]; j++) {
      if (_g_num1[j] != _g_num2[tmp_n-1]) {
        _g_num2[tmp_n++] = _g_num1[j];
      }
    }

    idx1 += block_box_n[i];
    idx2 += tmp_n;
    block_box_n[i] = tmp_n;
  }
  free (order);

  /* Fix partial block stride */
  PDM_g_num_t l_max_g_num = 0;
  for (int i = 0; i < n_line; i++) {
    l_max_g_num = PDM_MAX (l_max_g_num, line_g_num[i]);
  }

  PDM_g_num_t g_max_g_num;
  PDM_MPI_Allreduce (&l_max_g_num, &g_max_g_num, 1,
                     PDM__PDM_MPI_G_NUM, PDM_MPI_MAX, _dbbt->comm);

  PDM_g_num_t *block_distrib_idx =
    PDM_part_to_block_adapt_partial_block_to_block (ptb,
                                                    &block_box_n,
                                                    g_max_g_num);

  /* 2) Block-to-Part */
  PDM_block_to_part_t *btp = PDM_block_to_part_create (block_distrib_idx,
                                                       (const PDM_g_num_t **) &line_g_num,
                                                       &n_line,
                                                       1,
                                                       _dbbt->comm);

  int *box_n = malloc (sizeof(int) * n_line);
  int one = 1;
  PDM_block_to_part_exch_in_place (btp,
                          sizeof(int),
                          PDM_STRIDE_CST_INTERLACED,
                          &one,
                          block_box_n,
                          NULL,
                          (void **) &box_n);

  *box_idx = PDM_array_new_idx_from_sizes_int (box_n, n_line);

  *box_g_num = malloc (sizeof(PDM_g_num_t) * (*box_idx)[n_line]);
  PDM_block_to_part_exch_in_place (btp,
                          sizeof(PDM_g_num_t),
                          PDM_STRIDE_VAR_INTERLACED,
                          block_box_n,
                          block_box_g_num,
                          &box_n,
                          (void **) box_g_num);
  free (block_box_g_num);
  free (block_box_n);
  free (box_n);

  btp = PDM_block_to_part_free (btp);
  ptb = PDM_part_to_block_free (ptb);
  free (block_distrib_idx);
  free (part_n_line);

  free (part_line_g_num);
  if (line_g_num1 != line_g_num) free (line_g_num1);

  PDM_box_tree_free_copies(_dbbt->btLoc);
}



/**
 *
 * \brief Get an indexed list of all lines intersecting boxes
 * /!\ Results conform to the dbbtree's box partitionning
 *
 *   \param [in]  dbbt             Pointer to distributed box tree structure
 *   \param [in]  n_line           Number of points
 *   \param [in]  line_g_num       Line global ids (size = \ref n_line)
 *   \param [in]  line_coord       Line coordinates (size = 6 * \ref n_line)
 *   \param [out] box_line_idx     Index of lines (size = \ref n_box + 1, allocated inside function)
 *   \param [out] box_line_g_num   Global ids of lines (size = \ref box_idx[\ref n_box], allocated inside function)
 *
 */

void
PDM_dbbtree_lines_intersect_boxes2
(
 PDM_dbbtree_t  *dbbt,
 const int       n_line,
 PDM_g_num_t    *line_g_num,
 double         *line_coord,
 int            *n_part,
 int           **redistrib_n_box,
 PDM_g_num_t  ***redistrib_box_ln_to_gn,
 int          ***redistrib_box_init_location,
 int          ***box_line_idx,
 PDM_g_num_t  ***box_line_g_num
)
{
  assert (dbbt != NULL);
  _PDM_dbbtree_t *_dbbt = (_PDM_dbbtree_t *) dbbt;

  const float f_threshold = 1.1;  // factor of the mean nb of requests
  const float f_max_copy  = 0.1;  // factor of the total nb of processes

  /*
   *  Redistribute lines (with potential box-tree copies for load balancing)
   */
  int          *redistrib_n_line     = NULL;
  PDM_g_num_t **redistrib_line_g_num = NULL;
  double      **redistrib_line_coord = NULL;
  int           n_copied_ranks       = 0;
  _lines_intersect_shared_box_tree(dbbt,
                                   n_line,
                                   line_g_num,
                                   line_coord,
                                   f_threshold,
                                   f_max_copy,
                                   &redistrib_n_line,
                                   &redistrib_line_g_num,
                                   &redistrib_line_coord,
                                   &n_copied_ranks);

  if (0 == 1) {
    char filename[999];
    int i_rank;
    PDM_MPI_Comm_rank(_dbbt->comm, &i_rank);
    sprintf(filename, "redistrib_ray%2.2d.vtk", i_rank);
    PDM_vtk_write_lines(filename,
                        redistrib_n_line[0],
                        redistrib_line_coord[0],
                        redistrib_line_g_num[0],
                        NULL);
  }

  *n_part = 1 + n_copied_ranks;
  int _n_part = *n_part;

  *box_line_idx   = malloc (sizeof(int *)         * _n_part);

  /* Search in local box tree */
  int *box_line_l_num = NULL;
  PDM_box_tree_intersect_boxes_lines(_dbbt->btLoc,
                                      -1,
                                      redistrib_n_line[0],
                                      redistrib_line_coord[0],
                                      &((*box_line_idx)[0]),
                                      &box_line_l_num);
  free(redistrib_line_coord[0]);

  /*
   * Allocate and setup shortcut
   */
  *redistrib_n_box             = malloc(sizeof(int          ) * _n_part);
  *redistrib_box_ln_to_gn      = malloc(sizeof(PDM_g_num_t *) * _n_part);
  *redistrib_box_init_location = malloc(sizeof(int         *) * _n_part);
  *box_line_g_num              = malloc(sizeof(PDM_g_num_t *) * _n_part);

  int          *_redistrib_n_box             = *redistrib_n_box;
  PDM_g_num_t **_redistrib_box_ln_to_gn      = *redistrib_box_ln_to_gn;
  int         **_redistrib_box_init_location = *redistrib_box_init_location;
  PDM_g_num_t **_box_line_g_num              = *box_line_g_num;
  int         **_box_line_idx                = *box_line_idx;

  int **extract_box_line_idx = malloc(sizeof(int *) * _n_part);

  /*
   * Compress local
   */
  int n_boxes_local = _dbbt->boxes->local_boxes->n_boxes;
  PDM_g_num_t *boxes_gnum          = _dbbt->boxes->local_boxes->g_num;
  int         *boxes_init_location = _dbbt->boxes->local_boxes->origin;

  /* Count */
  _redistrib_n_box[0] = 0;
  for(int i_box = 0; i_box < n_boxes_local; ++i_box) {
    if(_box_line_idx [0][i_box+1] - _box_line_idx [0][i_box] > 0){
      _redistrib_n_box[0]++;
    }
  }

  int _n_line = _box_line_idx[0][n_boxes_local];
  _redistrib_box_ln_to_gn     [0] = malloc(    _redistrib_n_box[0]  * sizeof(PDM_g_num_t));
  _redistrib_box_init_location[0] = malloc(3 * _redistrib_n_box[0]  * sizeof(int        ));
  _box_line_g_num             [0] = malloc(    _n_line              * sizeof(PDM_g_num_t));
  extract_box_line_idx        [0] = malloc( (_redistrib_n_box[0]+1) * sizeof(int        ));

  /* Fill */
  _redistrib_n_box[0] = 0;
  extract_box_line_idx[0][0] = 0;
  for(int i_box = 0; i_box < n_boxes_local; ++i_box) {
    if(_box_line_idx[0][i_box+1] - _box_line_idx[0][i_box] > 0){

      extract_box_line_idx[0][_redistrib_n_box[0]+1] = extract_box_line_idx[0][_redistrib_n_box[0]];

      _redistrib_box_ln_to_gn     [0][_redistrib_n_box[0]] = boxes_gnum[i_box];

      _redistrib_box_init_location[0][3*_redistrib_n_box[0]  ] = boxes_init_location[3*i_box  ];
      _redistrib_box_init_location[0][3*_redistrib_n_box[0]+1] = boxes_init_location[3*i_box+1];
      _redistrib_box_init_location[0][3*_redistrib_n_box[0]+2] = boxes_init_location[3*i_box+2];

      for(int j = _box_line_idx [0][i_box]; j < _box_line_idx[0][i_box+1]; ++j) {
        _box_line_g_num[0][j] = redistrib_line_g_num[0][box_line_l_num[j]];
        extract_box_line_idx[0][_redistrib_n_box[0]+1]++;
      }
      _redistrib_n_box[0]++;
    }
  }
  free(box_line_l_num);
  free(_box_line_idx[0]);
  _box_line_idx[0] = extract_box_line_idx[0];

  /* Search in copied box trees */
  for (int i = 0; i < n_copied_ranks; i++) {
    PDM_box_tree_intersect_boxes_lines(_dbbt->btLoc,
                                        i,
                                        redistrib_n_line    [i+1],
                                        redistrib_line_coord[i+1],
                                        &((*box_line_idx)   [i+1]),
                                        &box_line_l_num);
    free(redistrib_line_coord[i+1]);

    /*
     * Compress
     */
    _redistrib_n_box[i+1] = 0;
    int n_boxes_copied = _dbbt->boxes->rank_boxes[i].n_boxes;
    PDM_g_num_t *copied_boxes_gnum          = _dbbt->boxes->rank_boxes[i].g_num;
    int         *copied_boxes_init_location = _dbbt->boxes->rank_boxes[i].origin;

    /* Count */
    _redistrib_n_box[i+1] = 0;
    for(int i_box = 0; i_box < n_boxes_copied; ++i_box) {
      if(_box_line_idx [i+1][i_box+1] - _box_line_idx [i+1][i_box] > 0){
        _redistrib_n_box[i+1]++;
      }
    }

    int _n_line_copied = _box_line_idx[i+1][n_boxes_copied];
    _redistrib_box_ln_to_gn     [i+1] = malloc(    _redistrib_n_box[i+1]  * sizeof(PDM_g_num_t));
    _redistrib_box_init_location[i+1] = malloc(3 * _redistrib_n_box[i+1]  * sizeof(int        ));
    _box_line_g_num             [i+1] = malloc(      _n_line_copied       * sizeof(PDM_g_num_t));
    extract_box_line_idx        [i+1] = malloc( (_redistrib_n_box[i+1]+1) * sizeof(int        ));

    _redistrib_n_box[i+1] = 0;
    extract_box_line_idx[i+1][0] = 0;
    for(int i_box = 0; i_box < n_boxes_copied; ++i_box) {
      if(_box_line_idx[i+1][i_box+1] - _box_line_idx[i+1][i_box] > 0){

        extract_box_line_idx[i+1][_redistrib_n_box[i+1]+1] = extract_box_line_idx[i+1][_redistrib_n_box[i+1]];

        _redistrib_box_ln_to_gn     [i+1][_redistrib_n_box[i+1]] = copied_boxes_gnum[i_box];

        _redistrib_box_init_location[i+1][3*_redistrib_n_box[i+1]  ] = copied_boxes_init_location[3*i_box  ];
        _redistrib_box_init_location[i+1][3*_redistrib_n_box[i+1]+1] = copied_boxes_init_location[3*i_box+1];
        _redistrib_box_init_location[i+1][3*_redistrib_n_box[i+1]+2] = copied_boxes_init_location[3*i_box+2];

        for(int j = _box_line_idx [i+1][i_box]; j < _box_line_idx [i+1][i_box+1]; ++j) {
          _box_line_g_num[i+1][j] = redistrib_line_g_num[i+1][box_line_l_num[j]];
          extract_box_line_idx[i+1][_redistrib_n_box[i+1]+1]++;
        }
        _redistrib_n_box[i+1]++;
      }
    }

    free(_box_line_idx[i+1]);
    _box_line_idx[i+1] = extract_box_line_idx[i+1];

    free(box_line_l_num);
    free(redistrib_line_g_num[i+1]);
  }
  free(redistrib_n_line);
  free(redistrib_line_coord);
  free(redistrib_line_g_num[0]);
  free(redistrib_line_g_num);
  free(extract_box_line_idx);

  /*
   *  Free memory
   */
  PDM_box_tree_free_copies(_dbbt->btLoc);
}

/**
 *
 * \brief Get an indexed list of all lines intersecting boxes
 * /!\ Results conform to the dbbtree's box partitionning
 *
 *   \param [in]  dbbt             Pointer to distributed box tree structure
 *   \param [in]  n_line           Number of points
 *   \param [in]  line_g_num       Line global ids (size = \ref n_line)
 *   \param [in]  line_coord       Line coordinates (size = 6 * \ref n_line)
 *   \param [out] box_line_idx     Index of lines (size = \ref n_box + 1, allocated inside function)
 *   \param [out] box_line_g_num   Global ids of lines (size = \ref box_idx[\ref n_box], allocated inside function)
 *
 */

void
PDM_dbbtree_lines_intersect_boxes2_shared
(
 PDM_dbbtree_t  *dbbt,
 const int       n_line,
 PDM_g_num_t    *line_g_num,
 double         *line_coord,
 int            *n_part,
 int           **redistrib_n_box,
 PDM_g_num_t  ***redistrib_box_ln_to_gn,
 int          ***redistrib_box_init_location,
 int          ***box_line_idx,
 PDM_g_num_t  ***box_line_g_num
)
{
  assert (dbbt != NULL);
  _PDM_dbbtree_t *_dbbt = (_PDM_dbbtree_t *) dbbt;

  int i_rank, n_rank;
  PDM_MPI_Comm_rank (_dbbt->comm, &i_rank);
  PDM_MPI_Comm_size (_dbbt->comm, &n_rank);

  // Shared
  PDM_MPI_Comm comm_shared;
  PDM_MPI_Comm_split_type(_dbbt->comm, PDM_MPI_SPLIT_NUMA, &comm_shared);

  int n_rank_in_shm, i_rank_in_shm;
  PDM_MPI_Comm_rank (comm_shared, &i_rank_in_shm);
  PDM_MPI_Comm_size (comm_shared, &n_rank_in_shm);

  /*
   *  Normalize coordinates
   */
  double *_line_coord = malloc (sizeof(double) * n_line * 6);
  for (int i = 0; i < 2*n_line; i++) {
    _normalize (_dbbt,
                line_coord  + 3*i,
                _line_coord + 3*i);
  }

  int n_line_local  = 0;
  int n_line_recv   = 0;
  int n_line_copied = 0;
  int n_line1       = 0;

  PDM_UNUSED(n_line_local);
  PDM_UNUSED(n_line_recv);
  PDM_UNUSED(n_line_copied);
  PDM_UNUSED(n_line1);
  PDM_g_num_t *line_g_num1 = NULL;
  double      *line_coord1 = NULL;

  PDM_mpi_win_shared_t* wshared_recv_gnum       = NULL;
  PDM_mpi_win_shared_t* wshared_recv_line_coord = NULL;

  int *distrib_search_by_rank_idx = NULL;

  if (_dbbt->btShared != NULL) {
    int *send_count = NULL;
    int *send_shift = NULL;
    int *recv_count = NULL;
    int *recv_shift = NULL;
    PDM_g_num_t *send_g_num = NULL;
    double      *send_coord = NULL;

    int *line_rank_idx = NULL;
    int *line_rank     = NULL;

    double t1a = PDM_MPI_Wtime();
    PDM_box_tree_intersect_lines_boxes (_dbbt->btShared,
                                        -1,
                                        n_line,
                                        line_coord,
                                        &line_rank_idx,
                                        &line_rank);
    double t2a = PDM_MPI_Wtime();
    log_trace("PDM_box_tree_intersect_lines_boxes (Shared) = %12.5e \n", t2a - t1a);

    /* Count points to send to each rank */
    send_count = PDM_array_zeros_int (n_rank);

    for (int i = 0; i < line_rank_idx[n_line]; i++) {
      int rank = _dbbt->usedRank[line_rank[i]];
      line_rank[i] = rank;
      send_count[rank]++;
    }

    // PDM_log_trace_array_int(send_count, _dbbt->nUsedRank, "send_count ::");

    recv_count = malloc (sizeof(int) * n_rank);
    PDM_MPI_Alltoall (send_count, 1, PDM_MPI_INT,
                      recv_count, 1, PDM_MPI_INT,
                      _dbbt->comm);
    send_shift = malloc ( ( n_rank + 1) * sizeof(int));
    recv_shift = malloc ( ( n_rank + 1) * sizeof(int));

    // Deduce size of recv buffer shared inside the same node
    int* shared_recv_count  = malloc(n_rank_in_shm * sizeof(int));

    int n_tot_recv = 0;
    send_shift[0] = 0;
    recv_shift[0] = 0;
    for(int i = 0; i < n_rank; ++i) {
      n_tot_recv += recv_count[i];
      send_shift[i+1] = send_shift[i] + send_count[i];
      recv_shift[i+1] = recv_shift[i] + recv_count[i];
    }

    PDM_MPI_Allgather(&n_tot_recv      , 1, PDM_MPI_INT,
                      shared_recv_count, 1, PDM_MPI_INT,
                      comm_shared);

    int* shared_recv_idx = malloc((n_rank_in_shm+1) * sizeof(int));
    shared_recv_idx[0] = 0;
    for(int i = 0; i < n_rank_in_shm; ++i) {
      shared_recv_idx[i+1] = shared_recv_idx[i] + shared_recv_count[i];
    }

    int n_tot_recv_shared = shared_recv_idx[n_rank_in_shm];
    wshared_recv_gnum       = PDM_mpi_win_shared_create(    n_tot_recv_shared, sizeof(PDM_g_num_t), comm_shared);
    wshared_recv_line_coord = PDM_mpi_win_shared_create(6 * n_tot_recv_shared, sizeof(double     ), comm_shared);

    PDM_g_num_t *shared_recv_gnum       = PDM_mpi_win_shared_get(wshared_recv_gnum);
    double      *shared_recv_line_coord = PDM_mpi_win_shared_get(wshared_recv_line_coord);

    PDM_mpi_win_shared_lock_all (0, wshared_recv_gnum);
    PDM_mpi_win_shared_lock_all (0, wshared_recv_line_coord);

    PDM_g_num_t *lrecv_gnum       = &shared_recv_gnum      [    shared_recv_idx[i_rank_in_shm]];
    double      *lrecv_line_coord = &shared_recv_line_coord[6 * shared_recv_idx[i_rank_in_shm]];

    /* Prepare send */
    send_g_num = malloc (sizeof(PDM_g_num_t) * send_shift[n_rank]);
    send_coord = malloc (sizeof(double)      * send_shift[n_rank] * 6);

    for(int i = 0; i < n_rank; ++i) {
      send_count[i] = 0;
    }

    for (int ipt = 0; ipt < n_line; ipt++) {
      for (int i = line_rank_idx[ipt]; i < line_rank_idx[ipt+1]; i++) {
        int t_rank = line_rank[i];

        int idx_write = send_shift[t_rank] + send_count[t_rank]++;
        send_g_num[idx_write] = line_g_num[ipt];
        for (int j = 0; j < 6; j++) {
          send_coord[6*idx_write + j] = _line_coord[6*ipt + j];
        }
      }
    }

    free(line_rank_idx);
    free(line_rank);


    PDM_MPI_Alltoallv (send_g_num, send_count, send_shift, PDM__PDM_MPI_G_NUM,
                       lrecv_gnum, recv_count, recv_shift, PDM__PDM_MPI_G_NUM,
                       _dbbt->comm);
    free(send_g_num);

    for (int i = 0; i < n_rank; i++) {
      send_count[i  ] *= 6;
      recv_count[i  ] *= 6;
      send_shift[i+1] *= 6;
      recv_shift[i+1] *= 6;
    }

    PDM_MPI_Alltoallv (send_coord      , send_count, send_shift, PDM_MPI_DOUBLE,
                       lrecv_line_coord, recv_count, recv_shift, PDM_MPI_DOUBLE,
                       _dbbt->comm);

    free (send_coord);
    free (send_count);
    free (send_shift);
    free (recv_count);
    free (recv_shift);

    PDM_MPI_Barrier (comm_shared);
    PDM_mpi_win_shared_sync (wshared_recv_gnum);
    PDM_mpi_win_shared_sync (wshared_recv_line_coord);

    /*
     * Redistribution by numa
     */
    PDM_g_num_t* distrib_search = PDM_compute_uniform_entity_distribution(comm_shared, n_tot_recv_shared);

    int  dn_search = distrib_search[i_rank_in_shm+1] - distrib_search[i_rank_in_shm];

    distrib_search_by_rank_idx = malloc((n_rank_in_shm+1) * sizeof(int));
    int* distrib_search_by_rank_n   = malloc((n_rank_in_shm  ) * sizeof(int));
    for(int i = 0; i < n_rank_in_shm; ++i) {
      distrib_search_by_rank_n[i] = 0;
    }

    // TODO : Faire un algo d'intersection de range pour ne pas faire la dicotomie x fois !
    for(int i = distrib_search[i_rank_in_shm]; i < distrib_search[i_rank_in_shm+1]; ++i) {
      int t_rank = PDM_binary_search_gap_int(i, shared_recv_idx, n_rank_in_shm+1);
      distrib_search_by_rank_n[t_rank]++;
    }

    distrib_search_by_rank_idx[0] = 0;
    for(int i = 0; i < n_rank_in_shm; ++i) {
      distrib_search_by_rank_idx[i+1] = distrib_search_by_rank_idx[i] + distrib_search_by_rank_n[i];
    }

    // PDM_log_trace_array_long(check, dn_search, "check ::");
    PDM_log_trace_array_int(distrib_search_by_rank_idx, n_rank_in_shm+1, "distrib_search_by_rank_idx ::");

    n_line_local  = 0;
    n_line_recv   = dn_search;
    n_line_copied = 0;

    n_line1       = dn_search;
    line_g_num1   = (PDM_g_num_t *) &shared_recv_gnum      [    distrib_search[i_rank_in_shm]];
    line_coord1   = (double      *) &shared_recv_line_coord[6 * distrib_search[i_rank_in_shm]];

    free(distrib_search);
    free(distrib_search_by_rank_n);
    free(shared_recv_count );
    free(shared_recv_idx );

  } else {
    n_line_local  = n_line;
    n_line_recv   = 0;
    n_line_copied = 0;

    n_line1 = n_line;

    line_g_num1 = (PDM_g_num_t *) line_g_num;
    line_coord1 = _line_coord;
  }

  int n_part_box = n_rank_in_shm;
  *n_part = n_part_box;

  *box_line_idx        = malloc (sizeof(int *) * n_part_box);

  *redistrib_n_box             = (int          *) malloc(sizeof(int          ) * n_part_box);
  *redistrib_box_ln_to_gn      = (PDM_g_num_t **) malloc(sizeof(PDM_g_num_t *) * n_part_box);
  *redistrib_box_init_location = (int         **) malloc(sizeof(int         *) * n_part_box);
  *box_line_g_num              = (PDM_g_num_t **) malloc(sizeof(PDM_g_num_t *) * n_part_box);

  int         **_box_line_idx                = *box_line_idx;
  int          *_redistrib_n_box             = *redistrib_n_box;
  PDM_g_num_t **_redistrib_box_ln_to_gn      = *redistrib_box_ln_to_gn;
  int         **_redistrib_box_init_location = *redistrib_box_init_location;
  PDM_g_num_t **_box_line_g_num              = *box_line_g_num;

  if(_dbbt->btShared != NULL) {

    double t1a = PDM_MPI_Wtime();
    for(int i_shm = 0; i_shm < n_rank_in_shm; ++i_shm) {

      int beg     = distrib_search_by_rank_idx[i_shm  ];
      int n_lline = distrib_search_by_rank_idx[i_shm+1] - beg;

      // part_n_line    [i_shm] = n_lline;

      PDM_g_num_t *lline_gnum  = &line_g_num1[  beg];
      double      *lline_coord = &line_coord1[6*beg];

      int *box_line_l_num = NULL;

      PDM_box_tree_intersect_boxes_lines_shared(_dbbt->btLoc,
                                                i_shm,
                                                n_lline,
                                                lline_coord,
                                                &_box_line_idx [i_shm],
                                                &box_line_l_num);

      /*
       * Hook box_tree and compress information
       */
      int n_boxes = _dbbt->boxes->shm_boxes[i_shm].n_boxes;
      PDM_g_num_t *boxes_gnum   = _dbbt->boxes->shm_boxes[i_shm].g_num;
      int         *boxes_origin = _dbbt->boxes->shm_boxes[i_shm].origin;

      /* Count */
      _redistrib_n_box[i_shm] = 0;
      for(int i_box = 0; i_box < n_boxes; ++i_box) {
        if(_box_line_idx [i_shm][i_box+1] - _box_line_idx [i_shm][i_box] > 0){
          _redistrib_n_box[i_shm]++;
        }
      }

      _redistrib_box_ln_to_gn     [i_shm] = malloc(    _redistrib_n_box[i_shm]      * sizeof(PDM_g_num_t));
      _redistrib_box_init_location[i_shm] = malloc(3 * _redistrib_n_box[i_shm]      * sizeof(int        ));
      _box_line_g_num             [i_shm] = malloc(_box_line_idx   [i_shm][n_boxes] * sizeof(PDM_g_num_t));

      int *extract_box_line_idx = malloc((_redistrib_n_box[i_shm]+1) * sizeof(int));

      /*
       * Translate in g_num
       */
      _redistrib_n_box[i_shm] = 0;
      extract_box_line_idx[0] = 0;
      for(int i_box = 0; i_box < n_boxes; ++i_box) {

        if(_box_line_idx [i_shm][i_box+1] - _box_line_idx [i_shm][i_box] > 0){

          _redistrib_box_ln_to_gn[i_shm][_redistrib_n_box[i_shm]] = boxes_gnum[i_box];
          extract_box_line_idx[_redistrib_n_box[i_shm]+1] = extract_box_line_idx[_redistrib_n_box[i_shm]];

          _redistrib_box_init_location[i_shm][3 * _redistrib_n_box[i_shm]  ] = boxes_origin[3*i_box  ];
          _redistrib_box_init_location[i_shm][3 * _redistrib_n_box[i_shm]+1] = boxes_origin[3*i_box+1];
          _redistrib_box_init_location[i_shm][3 * _redistrib_n_box[i_shm]+2] = boxes_origin[3*i_box+2];

          for(int idx_line = _box_line_idx [i_shm][i_box]; idx_line < _box_line_idx [i_shm][i_box+1]; ++idx_line) {
            int idx_write = extract_box_line_idx[_redistrib_n_box[i_shm]+1]++;
            _box_line_g_num[i_shm][idx_write] = lline_gnum[box_line_l_num[idx_line]];
          }

          _redistrib_n_box[i_shm]++;
        }
      }

      /*
       *  Replace
       */
      free(_box_line_idx [i_shm]);
      _box_line_idx [i_shm] = extract_box_line_idx;
      free(box_line_l_num);
    }
    double t2a = PDM_MPI_Wtime();
    log_trace("PDM_box_tree_intersect_lines_boxes (shm) = %12.5e \n", t2a-t1a);

    PDM_mpi_win_shared_unlock_all(wshared_recv_gnum);
    PDM_mpi_win_shared_unlock_all(wshared_recv_line_coord);
    PDM_mpi_win_shared_free (wshared_recv_gnum);
    PDM_mpi_win_shared_free (wshared_recv_line_coord);
    free(distrib_search_by_rank_idx);


  } else {

    int *box_line_l_num = NULL;
    PDM_box_tree_intersect_boxes_lines(_dbbt->btLoc,
                                       -1,
                                       n_line_local,
                                       line_coord1,
                                       &_box_line_idx[0],
                                       &box_line_l_num);
    /*
     * Hook box_tree and compress information
     */
    int n_boxes = _dbbt->boxes->local_boxes->n_boxes;
    PDM_g_num_t *boxes_gnum   = _dbbt->boxes->local_boxes->g_num;
    int         *boxes_origin = _dbbt->boxes->local_boxes->origin;

    /* Count */
    _redistrib_n_box[0] = 0;
    for(int i_box = 0; i_box < n_boxes; ++i_box) {
      if(_box_line_idx [0][i_box+1] - _box_line_idx [0][i_box] > 0){
        _redistrib_n_box[0]++;
      }
    }

    _redistrib_box_ln_to_gn     [0] = malloc(    _redistrib_n_box[0] * sizeof(PDM_g_num_t));
    _redistrib_box_init_location[0] = malloc(3 * _redistrib_n_box[0] * sizeof(int        ));
    _box_line_g_num             [0] = malloc(_box_line_idx   [0][n_boxes] * sizeof(PDM_g_num_t));

    int *extract_box_line_idx = malloc((_redistrib_n_box[0]+1) * sizeof(int));

    /*
     * Translate in g_num
     */
    _redistrib_n_box[0] = 0;
    extract_box_line_idx[0] = 0;
    for(int i_box = 0; i_box < n_boxes; ++i_box) {

      if(_box_line_idx [0][i_box+1] - _box_line_idx [0][i_box] > 0){

        _redistrib_box_ln_to_gn[0][_redistrib_n_box[0]] = boxes_gnum[i_box];
        extract_box_line_idx[_redistrib_n_box[0]+1] = extract_box_line_idx[_redistrib_n_box[0]];

        _redistrib_box_init_location[0][3 * _redistrib_n_box[0]  ] = boxes_origin[3*i_box  ];
        _redistrib_box_init_location[0][3 * _redistrib_n_box[0]+1] = boxes_origin[3*i_box+1];
        _redistrib_box_init_location[0][3 * _redistrib_n_box[0]+2] = boxes_origin[3*i_box+2];

        for(int idx_line = _box_line_idx [0][i_box]; idx_line < _box_line_idx [0][i_box+1]; ++idx_line) {
          int idx_write = extract_box_line_idx[_redistrib_n_box[0]+1]++;
          _box_line_g_num[0][idx_write] = line_g_num1[box_line_l_num[idx_line]];
        }

        _redistrib_n_box[0]++;
      }
    }

    /*
     *  Replace
     */
    free(_box_line_idx [0]);
    _box_line_idx [0] = extract_box_line_idx;
    free(box_line_l_num);
  }


  free (_line_coord);

  PDM_MPI_Comm_free(&comm_shared);
}

/**
 *
 * \brief Get an indexed list of all boxes intersecting volumes
 *
 * \param [in]   dbbt                  Pointer to distributed box tree structure
 * \param [in]   n_volumes             Number of volumes
 * \param [in]   volume_g_num          Global number of volumes
 * \param [in]   volume_plane_idx      Index of the number of planes per volume
 * \param [in]   plane_normal          Oriented normal vector for a given plane (oriented toward the interior of the volume)
 * \param [in]   plane_pt_coord        Point on plane coordinates (xa0, ya0, za0, xb0, yb0, zb0, xa1, ...)
 * \param [out]  volume_box_idx        Index of boxes (size = \ref n_line + 1, allocated inside function)
 * \param [out]  volume_box_g_num      Global ids of boxes (size = \ref box_line_idx[\ref n_line], allocated inside function)
 *
 */

/**
 * TO DO Karmijn: - remplacer l'ensemble sur les volumes puis sur les plans
 *                  par une grosse boucle en mettant le volume_n_plan
 *                  au niveau de l.5840 irank_jsubtree_n_volume
 */

void
PDM_dbbtree_volumes_intersect_boxes
(
 PDM_dbbtree_t  *dbbt,
 const int       n_volumes,
 PDM_g_num_t    *volume_g_num,
 int            *volume_plane_idx,
 double         *plane_normal,
 double         *plane_pt_coord,
 int           **volume_box_idx,
 PDM_g_num_t   **out_volume_box_g_num
)
{
  const float f_threshold = 1.1;  // factor of the mean nb of requests
  const float f_max_copy  = 0.1;  // factor of the total nb of processes

  assert (dbbt != NULL);
  _PDM_dbbtree_t *_dbbt = (_PDM_dbbtree_t *) dbbt;

  int i_rank, n_rank;
  PDM_MPI_Comm_rank (_dbbt->comm, &i_rank);
  PDM_MPI_Comm_size (_dbbt->comm, &n_rank);

  int n_planes = volume_plane_idx[n_volumes];

  // Normalize plane points and normal vectors
  double *plane_pt_coord_normalized = malloc (sizeof(double) * n_planes * 3);
  double *plane_normal_normalized = malloc (sizeof(double) * n_planes * 3);
  for (int i = 0; i < n_planes; i++) {
    _normalize(_dbbt,
               plane_pt_coord            + 3*i,
               plane_pt_coord_normalized + 3*i);
    _normalize_normal_vector(_dbbt,
                             plane_normal            + 3*i,
                             plane_normal_normalized + 3*i);
  } // end loop on planes

  // Set up for counting the number of volumes that might intersect boxes on ranks
  int          isubtree_total_n_volume            = 0;    // Total number of volumes associated to sub-box_tree i
  int          n_copied_ranks                     = 0;    // Number of subtrees to copy
  int          n_copied_volumes                   = 0;    // number of volumes linked to copied subtrees on rank i
  int          n_receive_volumes                  = 0;    // number of received volumes linked to subtree i on rank i
  int         *irank_jsubtree_n_volume            = NULL; // Number of volumes per sub-box_tree on rank i
  int         *isubtree_jrank_n_volume            = NULL; // Number of volumes for sub-box_tree i on each rank j
  int         *all_jrank_jsubtree_total_n_volume  = NULL; // Total number of volumes associated to sub-box_tree j for each rank j
  int         *copied_ranks                       = NULL; // Index of copied subtrees
  PDM_g_num_t *send_volume_g_num                  = NULL;
  double      *send_plane_normal                  = NULL;
  double      *send_plane_pt_coord                = NULL;
  PDM_g_num_t *receive_volume_g_num               = NULL;
  double      *receive_plane_normal               = NULL;
  double      *receive_plane_pt_coord             = NULL;
  int         *i_copied_rank                      = NULL;
  int         *jsubtree_to_copy_total_n_volume    = NULL;
  int         *copied_volume_stride               = NULL; // portion of volumes linked to a copied subtree to keep on rank i
  int         *send_volume_stride                 = NULL; // stride of volumes that need to be sent to other ranks
  int         *receive_volume_stride              = NULL;
  int         *copied_volume_idx                  = NULL;
  int         *copied_plane_idx                   = NULL;
  int         *send_volume_idx                    = NULL;
  int         *receive_volume_idx                 = NULL;
  int         *send_plane_idx                     = NULL; // because there might not be the same amount of planes per volume
  int         *receive_plane_idx                  = NULL; // because there might not be the same amount of planes per volume
  int         *volume_subtree_idx                 = NULL;
  int         *volume_subtree                     = NULL;
  int         *copied_rank_n_plane                = NULL;
  PDM_g_num_t *copied_volume_g_num                = NULL;
  int         *copied_volume_plane_idx            = NULL;
  double      *copied_plane_normal                = NULL;
  double      *copied_plane_pt_coord              = NULL;
  int         *copied_rank_volume_head            = NULL;
  int         *send_rank_volume_head              = NULL;
  int         *send_rank_n_plane                  = NULL;
  int         *receive_rank_n_plane               = NULL;
  int         *send_rank_plane_head               = NULL;
  int         *copied_rank_plane_head             = NULL;
  int         *send_rank_volume_n_plane           = NULL;
  int         *receive_rank_volume_n_plane        = NULL;

  int   n_part = 0;
  int  *part_n_volumes          = NULL;
  int **n_part_volume_box_idx   = NULL;
  int **n_part_volume_box_l_num = NULL;

  if (_dbbt->btShared != NULL) {

    // For each volume, find all ranks that might have boxes intersecting that volume
    // The shared box tree shouldn't be normalized
    PDM_box_tree_intersect_volume_boxes(_dbbt->btShared,
                                        -1,
                                        n_volumes,
                                        volume_plane_idx,
                                        plane_normal,
                                        plane_pt_coord,
                                        &volume_subtree_idx,
                                        &volume_subtree);

    // Count on irank the number of volumes associated to each j sub-box_tree
    irank_jsubtree_n_volume = PDM_array_zeros_int(n_rank);

    for (int k = 0; k < volume_subtree_idx[n_volumes]; k++) {
      int jsubtree = _dbbt->usedRank[volume_subtree[k]];
      volume_subtree[k] = jsubtree;
      irank_jsubtree_n_volume[jsubtree]++;
    } // end loop on volumes associated to subtrees on rank i

    isubtree_jrank_n_volume = malloc(sizeof(int) * n_rank);
    PDM_MPI_Alltoall(irank_jsubtree_n_volume, 1, PDM_MPI_INT,
                     isubtree_jrank_n_volume, 1, PDM_MPI_INT,
                     _dbbt->comm);

    isubtree_total_n_volume = 0;
    for (int jrank = 0; jrank < n_rank; jrank++) {
      isubtree_total_n_volume += isubtree_jrank_n_volume[jrank];
    }

    // Determine mean value of received volumes
    int  n_max_copy                      = (int) (f_max_copy * n_rank);
    int  mean_jsubtree_total_n_volume    = 0;

    if (n_max_copy > 0) {
      all_jrank_jsubtree_total_n_volume = malloc(sizeof(int) * n_rank);
      PDM_MPI_Allgather(&isubtree_total_n_volume,    1, PDM_MPI_INT,
                        all_jrank_jsubtree_total_n_volume, 1, PDM_MPI_INT,
                        _dbbt->comm);

      long long_tmp_mean = 0;
      for (int j = 0; j < n_rank; j++) {
        long_tmp_mean += all_jrank_jsubtree_total_n_volume[j];
      }
      mean_jsubtree_total_n_volume = (int) (long_tmp_mean / n_rank);

      float n_threshold = f_threshold * mean_jsubtree_total_n_volume;

      // Sort subtrees (i on rank i) according to the total number of associated volumes
      int *order = malloc(sizeof(int) * n_rank);
      for (int i = 0; i < n_rank; i++) {
        order[i] = i;
      }

      PDM_sort_int(all_jrank_jsubtree_total_n_volume,
                   order,
                   n_rank);

      // Identify ranks to copy
      copied_ranks                    = malloc (sizeof(int) * n_max_copy);
      jsubtree_to_copy_total_n_volume = malloc (sizeof(int) * n_max_copy);
      for (int i = 0; i < n_max_copy; i++) {
        int j = n_rank - i - 1;

        if (all_jrank_jsubtree_total_n_volume[j] > n_threshold) {
          copied_ranks[n_copied_ranks] = order[j];
          jsubtree_to_copy_total_n_volume[n_copied_ranks] = all_jrank_jsubtree_total_n_volume[j];
          n_copied_ranks++;
        }
        else {
          break;
        }
      }
      free(all_jrank_jsubtree_total_n_volume);
      free(order);

      if (n_copied_ranks > 0) {
        copied_ranks = realloc(copied_ranks, sizeof(int) * n_copied_ranks);
        jsubtree_to_copy_total_n_volume = realloc(jsubtree_to_copy_total_n_volume,
                                           sizeof(int) * n_copied_ranks);
      }

    } // end if we want to copy subtrees

    // Copy subtree selected to be copied to all other ranks
    i_copied_rank = malloc(sizeof(int) * n_rank);
    PDM_box_tree_copy_to_ranks(_dbbt->btLoc,
                               &n_copied_ranks,
                               copied_ranks,
                               i_copied_rank);

    // Determine which fraction of copied rank we keep
    copied_volume_stride = PDM_array_zeros_int (n_copied_ranks);
    for (int i = 0; i < n_copied_ranks; i++) {
      int rank = copied_ranks[i];
      if (rank != i_rank) {
        int keep_portion = irank_jsubtree_n_volume[rank];

        keep_portion = PDM_MIN (keep_portion, PDM_MAX (0, (jsubtree_to_copy_total_n_volume[i] - isubtree_total_n_volume)/2));
        if (i_copied_rank[i_rank] < 0) {
          keep_portion = PDM_MIN (keep_portion, PDM_MAX (0, mean_jsubtree_total_n_volume - isubtree_total_n_volume));
        }

        copied_volume_stride[i] = keep_portion;
        isubtree_total_n_volume += keep_portion;
      } // edn if not i rank
    } // end loop on copied ranks

    if (copied_ranks != NULL) {
      free (copied_ranks);
    }

    if (jsubtree_to_copy_total_n_volume != NULL) {
      free (jsubtree_to_copy_total_n_volume);
    }

    // Create stride to send local data and copied_not_kept data
    for (int i = 0; i < n_rank; i++) {
      // i == i_rank: everything there is is sent
      if (i_copied_rank[i] >= 0) {
        irank_jsubtree_n_volume[i] -= copied_volume_stride[i_copied_rank[i]]; // remove copied_kept part
      }
    }

    // Create copy idx
    copied_volume_idx   = PDM_array_new_idx_from_sizes_int(copied_volume_stride, n_copied_ranks);
    n_copied_volumes    = copied_volume_idx[n_copied_ranks];
    copied_rank_n_plane = PDM_array_zeros_int(n_rank);

    // Get for subtree i the amount of volumes rank i will get from rank j
    // Note: I will receive some from myself (MPI doesn't put that on the network)
    PDM_MPI_Alltoall (irank_jsubtree_n_volume, 1, PDM_MPI_INT, // send volume stride
                      isubtree_jrank_n_volume, 1, PDM_MPI_INT, // receive volume stride
                      _dbbt->comm);

    send_volume_idx       = PDM_array_new_idx_from_sizes_int(irank_jsubtree_n_volume, n_rank); // irank_jsubtree_n_volume i as send_volume_stride
    receive_volume_idx    = PDM_array_new_idx_from_sizes_int(isubtree_jrank_n_volume, n_rank); // isubtree_jrank_n_volume is a receive_volume_stride
    send_volume_stride    = irank_jsubtree_n_volume; // updated with copied-kept
    receive_volume_stride = isubtree_jrank_n_volume;

    n_receive_volumes = receive_volume_idx[n_rank];

    // Allocate table for copied and local-rceived data on rank i
    int tmp_n_copied_planes                  = 3 * n_copied_volumes * 2;
    copied_volume_g_num                      = malloc(sizeof(PDM_g_num_t) * n_copied_volumes);
    copied_plane_normal                      = malloc(sizeof(double     ) * tmp_n_copied_planes);
    copied_plane_pt_coord                    = malloc(sizeof(double     ) * tmp_n_copied_planes);
    copied_rank_volume_head                  = PDM_array_zeros_int(n_copied_ranks);

    // Allocate table for data sent by rank i to ranks j
    send_volume_g_num                  = malloc(sizeof(PDM_g_num_t) * send_volume_idx[n_rank]);
    send_rank_volume_n_plane           = malloc(sizeof(int) * send_volume_idx[n_rank]);
    send_rank_volume_head              = PDM_array_zeros_int(n_rank); // for each rank points on current table head
    send_rank_n_plane                  = PDM_array_zeros_int(n_rank);

     // Initialize index
    copied_volume_plane_idx  = malloc(sizeof(int) * (n_copied_ranks+copied_volume_idx[n_copied_ranks]));
    for (int i = 0; i < n_rank; i++) {
      if (i_copied_rank[i] > -1) {
        int _rank = i_copied_rank[i];
        copied_volume_plane_idx[copied_volume_idx[_rank] + _rank] = 0;
      }
    }

    // Order data in table between what remains on rank i and what is sent to ranks j (on volumes)
    for (int ivol = 0; ivol < n_volumes; ivol++) {
      for (int i = volume_subtree_idx[ivol]; i < volume_subtree_idx[ivol+1]; i++) {
        int rank = volume_subtree[i];

        // copied rank
        if (i_copied_rank[rank] >= 0) {
          int _rank = i_copied_rank[rank]; // how manyth copied rank

          // part of copied that we keep
          if (copied_rank_volume_head[_rank] < copied_volume_stride[_rank]) {
            int _rank_idx = copied_volume_idx[_rank] + copied_rank_volume_head[_rank];
            copied_volume_plane_idx[_rank_idx+1+_rank] = copied_volume_plane_idx[_rank_idx+_rank] + (volume_plane_idx[ivol+1] - volume_plane_idx[ivol]);
            copied_volume_g_num[_rank_idx]  = volume_g_num[ivol];
            copied_rank_n_plane[_rank] += volume_plane_idx[ivol+1] - volume_plane_idx[ivol];
            copied_rank_volume_head[_rank]++;
          } // end if kept-copied part

          // part of copied that we send
          else {
            int rank_idx = send_volume_idx[rank] + send_rank_volume_head[rank];
            send_rank_volume_n_plane[rank_idx] = volume_plane_idx[ivol+1] - volume_plane_idx[ivol];
            send_volume_g_num[rank_idx] =  volume_g_num[ivol];
            send_rank_n_plane[rank]     += volume_plane_idx[ivol+1] - volume_plane_idx[ivol];
            send_rank_volume_head[rank]++;
          } // end if sent-copied part

        } // end if copied rank

        // not-copied rank
        else {
          int rank_idx = send_volume_idx[rank] + send_rank_volume_head[rank];
          send_rank_volume_n_plane[rank_idx] = volume_plane_idx[ivol+1] - volume_plane_idx[ivol];
          send_volume_g_num[rank_idx] =  volume_g_num[ivol];
          send_rank_n_plane[rank]     += volume_plane_idx[ivol+1] - volume_plane_idx[ivol];
          send_rank_volume_head[rank]++;
        } // end if not-copied rank

      } // end loop on subtrees associated to volume ivol
    } // end loop on volumes

    // Exchange volume g_num
    receive_volume_g_num = malloc(sizeof(PDM_g_num_t) * n_receive_volumes);
    PDM_MPI_Alltoallv(send_volume_g_num, send_volume_stride, send_volume_idx, PDM__PDM_MPI_G_NUM,
                      receive_volume_g_num, receive_volume_stride, receive_volume_idx, PDM__PDM_MPI_G_NUM,
                      _dbbt->comm);
    free (send_volume_g_num);

    // Exhange plane stride
    receive_rank_n_plane = malloc(sizeof(int) * n_rank);
    PDM_MPI_Alltoall(send_rank_n_plane, 1, PDM_MPI_INT,
                     receive_rank_n_plane, 1, PDM_MPI_INT,
                     _dbbt->comm);

    // Order data in table between what remains on rank i and what is sent to ranks j (on planes)
    send_rank_plane_head     = PDM_array_zeros_int(n_rank);
    copied_rank_plane_head   = PDM_array_zeros_int(n_copied_ranks);
    send_plane_idx           = PDM_array_new_idx_from_sizes_int(send_rank_n_plane, n_rank);
    receive_plane_idx        = PDM_array_new_idx_from_sizes_int(receive_rank_n_plane, n_rank);
    copied_plane_idx         = PDM_array_new_idx_from_sizes_int(copied_rank_n_plane, n_copied_ranks);

    int tmp_n_send_planes = 3 * send_plane_idx[n_rank];
    send_plane_normal     = malloc(sizeof(double     ) * tmp_n_send_planes);
    send_plane_pt_coord   = malloc(sizeof(double     ) * tmp_n_send_planes);

    PDM_array_reset_int(copied_rank_volume_head, n_copied_ranks, 0);

    for (int ivol = 0; ivol < n_volumes; ivol++) {
      for (int i = volume_subtree_idx[ivol]; i < volume_subtree_idx[ivol+1]; i++) {
        int rank = volume_subtree[i];

        // copied rank
        if (i_copied_rank[rank] >= 0) {
          int _rank = i_copied_rank[rank]; // how manyth copied rank

          // part of copied that we keep
          if (copied_rank_volume_head[_rank] < copied_volume_stride[_rank]) {
            for (int iplane = volume_plane_idx[ivol]; iplane < volume_plane_idx[ivol+1]; iplane++) {
              int _rank_idx_plane = copied_plane_idx[_rank] + copied_rank_plane_head[_rank];
              // Realloc if necessary
              if ( 3*(_rank_idx_plane+1) > tmp_n_copied_planes) {
                tmp_n_copied_planes *= 2;
                copied_plane_normal   = realloc(copied_plane_normal   , sizeof(double) * tmp_n_copied_planes);
                copied_plane_pt_coord = realloc(copied_plane_pt_coord , sizeof(double) * tmp_n_copied_planes);
              }
              for (int j = 0; j < 3; j++) {
                copied_plane_normal[3*_rank_idx_plane + j]   = plane_normal_normalized[3*iplane + j];
                copied_plane_pt_coord[3*_rank_idx_plane + j] = plane_pt_coord_normalized[3*iplane + j];
              }
              copied_rank_plane_head[_rank]++;
            } // end loop on ivol planes
            copied_rank_volume_head[_rank]++;
          } // end if kept-copied part

          // part of copied that we send
          else {
            for (int iplane = volume_plane_idx[ivol]; iplane < volume_plane_idx[ivol+1]; iplane++) {
              int rank_idx_plane = send_plane_idx[rank] + send_rank_plane_head[rank];
              // Realloc if necessary
              if ( 3*(rank_idx_plane+1) > tmp_n_send_planes) {
                tmp_n_send_planes *= 2;
                send_plane_normal   = realloc(send_plane_normal   , sizeof(double) * tmp_n_send_planes);
                send_plane_pt_coord = realloc(send_plane_pt_coord , sizeof(double) * tmp_n_send_planes);
              }
              for (int j = 0; j < 3; j++) {
                send_plane_normal[3* rank_idx_plane + j]   = plane_normal_normalized[3*iplane + j];
                send_plane_pt_coord[3* rank_idx_plane + j] = plane_pt_coord_normalized[3*iplane + j];
              }
              send_rank_plane_head[rank]++;
            } // end loop on ivol planes
          } // end if sent-copied part

        } // end if copied rank

        // not-copied rank
        else {
          for (int iplane = volume_plane_idx[ivol]; iplane < volume_plane_idx[ivol+1]; iplane++) {
            int rank_idx_plane = send_plane_idx[rank] + send_rank_plane_head[rank];
            // Realloc if necessary
            if ( 3*(rank_idx_plane+1) > tmp_n_send_planes) {
              tmp_n_send_planes *= 2;
              send_plane_normal   = realloc(send_plane_normal   , sizeof(double) * tmp_n_send_planes);
              send_plane_pt_coord = realloc(send_plane_pt_coord , sizeof(double) * tmp_n_send_planes);
            }
            for (int j = 0; j < 3; j++) {
              send_plane_normal[3*rank_idx_plane + j]   = plane_normal_normalized[3*iplane + j];
              send_plane_pt_coord[3*rank_idx_plane + j] = plane_pt_coord_normalized[3*iplane + j];
            }
            send_rank_plane_head[rank]++;
          } // end loop on ivol planes
        } // end if not-copied rank

      } // end loop on subtrees associated to volume ivol
    } // end loop on volumes

    // Exchange to recreate index
    receive_rank_volume_n_plane = malloc(sizeof(int) * receive_volume_idx[n_rank]);
    PDM_MPI_Alltoallv(send_rank_volume_n_plane, send_volume_stride, send_volume_idx, PDM_MPI_INT,
                      receive_rank_volume_n_plane, receive_volume_stride, receive_volume_idx, PDM_MPI_INT,
                     _dbbt->comm);

    for (int i = 0; i < n_rank; i++) {
      send_rank_n_plane[i]    *= 3;
      receive_rank_n_plane[i] *= 3;
      send_plane_idx[i+1]     *= 3;
      receive_plane_idx[i+1]  *= 3;
    }

    send_plane_normal   = realloc(send_plane_normal   , sizeof(double) * send_plane_idx[n_rank]);
    send_plane_pt_coord = realloc(send_plane_pt_coord , sizeof(double) * send_plane_idx[n_rank]);

    // Exchange plane data
    receive_plane_normal = malloc(sizeof(double) * 3 * receive_plane_idx[n_rank]);
    PDM_MPI_Alltoallv(send_plane_normal, send_rank_n_plane, send_plane_idx, PDM_MPI_DOUBLE,
                      receive_plane_normal, receive_rank_n_plane, receive_plane_idx, PDM_MPI_DOUBLE,
                      _dbbt->comm);

    receive_plane_pt_coord = malloc(sizeof(double) * 3 * receive_plane_idx[n_rank]);
    PDM_MPI_Alltoallv(send_plane_pt_coord, send_rank_n_plane, send_plane_idx, PDM_MPI_DOUBLE,
                      receive_plane_pt_coord, receive_rank_n_plane, receive_plane_idx, PDM_MPI_DOUBLE,
                      _dbbt->comm);



    // Get boxes intersecting volumes
    n_part = 1 + n_copied_ranks;
    part_n_volumes = malloc (sizeof(int) * n_part);
    n_part_volume_box_idx   = malloc (sizeof(int *) * n_part);
    n_part_volume_box_l_num = malloc (sizeof(int *) * n_part);

    // -> on copied subtrees
    if (n_copied_ranks > 0) {

      for (int i = 0; i < n_copied_ranks; i++) {
          part_n_volumes[i+1] = copied_volume_idx[i+1] - copied_volume_idx[i];

          PDM_box_tree_intersect_volume_boxes(_dbbt->btLoc,
                                              i,
                                              part_n_volumes[i+1],
                                              copied_volume_plane_idx + copied_volume_idx[i] + i,
                                              copied_plane_normal     + copied_plane_idx[i]  * 3,
                                              copied_plane_pt_coord   + copied_plane_idx[i]  * 3,
                                              &(n_part_volume_box_idx[i+1]),
                                              &(n_part_volume_box_l_num[i+1]));

      } // end loop on copied ranks
    } // end if there are copied ranks

    // --> on local subtree

    int *receive_volume_plane_idx = malloc( sizeof(int) * (n_receive_volumes +1));
    receive_volume_plane_idx[0] = 0;
    for (int ivol = 0; ivol < n_receive_volumes; ivol++) {
      receive_volume_plane_idx[ivol+1] = receive_volume_plane_idx[ivol] + receive_rank_volume_n_plane[ivol];
    }

    part_n_volumes[0] = n_receive_volumes;

    PDM_box_tree_intersect_volume_boxes(_dbbt->btLoc,
                                        -1,
                                        part_n_volumes[0],
                          (const int *) receive_volume_plane_idx,
                                        receive_plane_normal,
                                        receive_plane_pt_coord,
                                        &(n_part_volume_box_idx[0]),
                                        &(n_part_volume_box_l_num[0]));


    free(receive_rank_volume_n_plane);
    free(receive_volume_plane_idx);

  } // end if there is a shared bt
  else {

    n_part = 1;
    part_n_volumes = malloc (sizeof(int) * n_part);
    n_part_volume_box_idx   = malloc (sizeof(int *) * n_part);
    n_part_volume_box_l_num = malloc (sizeof(int *) * n_part);
    part_n_volumes[0] = n_volumes;

    PDM_box_tree_intersect_volume_boxes(_dbbt->btLoc,
                                        -1,
                                        part_n_volumes[0],
                                        volume_plane_idx,
                                        plane_normal_normalized,
                                        plane_pt_coord_normalized,
                                        &(n_part_volume_box_idx[0]),
                                        &(n_part_volume_box_l_num[0]));
  }

  // from l_num to g_num
  PDM_g_num_t **volume_box_g_num = malloc (sizeof(PDM_g_num_t *) * n_part);
  volume_box_g_num[0] = malloc (sizeof(PDM_g_num_t) * n_part_volume_box_idx[0][part_n_volumes[0]]);
  for (int j = 0; j < part_n_volumes[0]; j++) {
    for (int k = n_part_volume_box_idx[0][j]; k < n_part_volume_box_idx[0][j+1]; k++) {
      volume_box_g_num[0][k] = _dbbt->boxes->local_boxes->g_num[n_part_volume_box_l_num[0][k]];
    }
  }

  for (int i = 0; i < n_copied_ranks; i++) {
    volume_box_g_num[i+1] = malloc (sizeof(PDM_g_num_t) * n_part_volume_box_idx[i+1][part_n_volumes[i+1]]);
    for (int j = 0; j < part_n_volumes[i+1]; j++) {
      for (int k = n_part_volume_box_idx[i+1][j]; k < n_part_volume_box_idx[i+1][j+1]; k++) {
        volume_box_g_num[i+1][k] = _dbbt->boxes->rank_boxes[i].g_num[n_part_volume_box_l_num[i+1][k]];
      }
    }
  }

  // Get weights of elements
  PDM_g_num_t **part_volume_g_num = malloc (sizeof(PDM_g_num_t *) * n_part);
  int    **part_stride = malloc (sizeof(int *)    * n_part);
  double **part_weight = malloc (sizeof(double *) * n_part);
  int idx = 0;
  for (int ipart = 0; ipart < n_part; ipart++) {
    if (_dbbt->btShared != NULL) {
      if (ipart == 0) {
         part_volume_g_num[ipart] = receive_volume_g_num;
      } else {
          part_volume_g_num[ipart] = copied_volume_g_num + idx;
          idx += part_n_volumes[ipart];
      }
    } else {
      part_volume_g_num[ipart] = volume_g_num;
    }

    part_stride[ipart] = malloc (sizeof(int)    * part_n_volumes[ipart]);
    part_weight[ipart] = malloc (sizeof(double) * part_n_volumes[ipart]);
    for (int i = 0; i < part_n_volumes[ipart]; i++) {
      part_stride[ipart][i] = n_part_volume_box_idx[ipart][i+1] - n_part_volume_box_idx[ipart][i];
      part_weight[ipart][i] = (double) part_stride[ipart][i];
    }
  }

  // Merge results
  /* 1) Part-to-Block */
  PDM_part_to_block_t *ptb = PDM_part_to_block_create (PDM_PART_TO_BLOCK_DISTRIB_ALL_PROC,
                                                       PDM_PART_TO_BLOCK_POST_MERGE,
                                                       1.,
                                                       part_volume_g_num,
                                                       part_weight,
                                                       part_n_volumes,
                                                       n_part,
                                                       _dbbt->comm);

  int *block_box_n = NULL;
  PDM_g_num_t *block_box_g_num = NULL;
  PDM_part_to_block_exch (ptb,
                          sizeof(PDM_g_num_t),
                          PDM_STRIDE_VAR_INTERLACED,
                          1,
                          part_stride,
                          (void **) volume_box_g_num,
                          &block_box_n,
                          (void **) &block_box_g_num);

  for (int ipart = 0; ipart < n_part; ipart++) {
    free (part_stride[ipart]);
    free (part_weight[ipart]);
    free (volume_box_g_num[ipart]);
  }
  free (part_stride);
  free (part_weight);
  free (volume_box_g_num);


  /* Remove doubles */
  int idx1 = 0, idx2 = 0;
  int n_volume_block = PDM_part_to_block_n_elt_block_get (ptb);
  int max_n = 0;
  for (int i = 0; i < n_volume_block; i++) {
    max_n = PDM_MAX (max_n, block_box_n[i]);
  }

  int *order = malloc (sizeof(int) * max_n);
  idx1 = 0;
  idx2 = 0;
  for (int i = 0; i < n_volume_block; i++) {
    if (block_box_n[i] == 0) continue;

    PDM_g_num_t *_g_num1 = block_box_g_num + idx1;
    PDM_g_num_t *_g_num2 = block_box_g_num + idx2;

    for (int j = 0; j < block_box_n[i]; j++) {
      order[j] = j;
    }
    PDM_sort_long (_g_num1,
                   order,
                   block_box_n[i]);

    _g_num2[0] = _g_num1[0];
    int tmp_n = 1;
    for (int j = 1; j < block_box_n[i]; j++) {
      if (_g_num1[j] != _g_num2[tmp_n-1]) {
        _g_num2[tmp_n++] = _g_num1[j];
      }
    }

    idx1 += block_box_n[i];
    idx2 += tmp_n;
    block_box_n[i] = tmp_n;
  }
  free (order);

  /* Fix partial block stride */
  PDM_g_num_t l_max_g_num = 0;
  for (int i = 0; i < n_volumes; i++) {
    l_max_g_num = PDM_MAX (l_max_g_num, volume_g_num[i]);
  }

  PDM_g_num_t g_max_g_num;
  PDM_MPI_Allreduce (&l_max_g_num, &g_max_g_num, 1,
                     PDM__PDM_MPI_G_NUM, PDM_MPI_MAX, _dbbt->comm);

  PDM_g_num_t *block_distrib_idx =
    PDM_part_to_block_adapt_partial_block_to_block (ptb,
                                                    &block_box_n,
                                                    g_max_g_num);

  /* 2) Block-to-Part */
  PDM_block_to_part_t *btp = PDM_block_to_part_create (block_distrib_idx,
                                                       (const PDM_g_num_t **) &volume_g_num,
                                                       &n_volumes,
                                                       1,
                                                       _dbbt->comm);

  int *box_n = malloc (sizeof(int) * n_volumes);
  int one = 1;
  PDM_block_to_part_exch_in_place (btp,
                          sizeof(int),
                          PDM_STRIDE_CST_INTERLACED,
                          &one,
                          block_box_n,
                          NULL,
                          (void **) &box_n);

  *volume_box_idx = PDM_array_new_idx_from_sizes_int (box_n, n_volumes);

  *out_volume_box_g_num = malloc (sizeof(PDM_g_num_t) * (*volume_box_idx)[n_volumes]);
  PDM_block_to_part_exch_in_place (btp,
                          sizeof(PDM_g_num_t),
                          PDM_STRIDE_VAR_INTERLACED,
                          block_box_n,
                          block_box_g_num,
                          &box_n,
                          (void **) out_volume_box_g_num);
  free (block_box_g_num);
  free (block_box_n);
  free (box_n);

  PDM_block_to_part_free(btp);
  PDM_part_to_block_free(ptb);
  free (block_distrib_idx);
  free (part_volume_g_num);

  free (part_n_volumes);
  for (int i = 0; i < n_part; i++) {
    free(n_part_volume_box_idx[i]);
    free(n_part_volume_box_l_num[i]);
  }
  free(n_part_volume_box_idx);
  free(n_part_volume_box_l_num);

  if (copied_volume_g_num != NULL)     free (copied_volume_g_num);
  if (copied_volume_plane_idx != NULL) free (copied_volume_plane_idx);
  if (copied_plane_normal != NULL)     free (copied_plane_normal);
  if (copied_plane_pt_coord != NULL)   free (copied_plane_pt_coord);

  free(plane_pt_coord_normalized);
  free(plane_normal_normalized);
  if (_dbbt->btShared != NULL) {
    free(irank_jsubtree_n_volume);
    free(isubtree_jrank_n_volume);
    free(send_plane_normal);
    free(send_plane_pt_coord);
    free(receive_volume_g_num);
    free(receive_plane_normal);
    free(receive_plane_pt_coord);
    free(i_copied_rank);
    free(copied_volume_stride);
    free(copied_volume_idx);
    free(copied_plane_idx);
    free(send_volume_idx);
    free(receive_volume_idx);
    free(send_plane_idx);
    free(receive_plane_idx);
    free(volume_subtree_idx);
    free(volume_subtree);
    free(copied_rank_n_plane);
    free(copied_rank_volume_head);
    free(send_rank_volume_head);
    free(send_rank_n_plane);
    free(receive_rank_n_plane);
    free(send_rank_plane_head);
    free(copied_rank_plane_head);
    free(send_rank_volume_n_plane);
  }


  PDM_box_tree_free_copies(_dbbt->btLoc);
}

void
PDM_dbbtree_box_tree_write_vtk
(
 const char    *filename,
 PDM_dbbtree_t *dbbt,
 const int      i_copied_rank,
 const int      normalized
 )
{
  assert (dbbt != NULL);
  _PDM_dbbtree_t *_dbbt = (_PDM_dbbtree_t *) dbbt;

  PDM_box_tree_write_vtk(filename,
                         _dbbt->btLoc,
                         i_copied_rank,
                         normalized);
}


void
PDM_dbbtree_box_tree_write_vtk2
(
 const char    *filename,
 PDM_dbbtree_t *dbbt,
 const int      i_copied_rank
 )
{
  assert (dbbt != NULL);
  _PDM_dbbtree_t *_dbbt = (_PDM_dbbtree_t *) dbbt;

  PDM_box_tree_write_vtk2(filename,
                          _dbbt->btLoc,
                          i_copied_rank,
                          _dbbt->s,
                          _dbbt->d);
}

#undef _MIN
#undef _MAX

#ifdef __cplusplus
}
#endif /* __cplusplus */
