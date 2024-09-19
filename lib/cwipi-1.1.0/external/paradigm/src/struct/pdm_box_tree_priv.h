#ifndef __PDM_BOX_TREE_PRIV_H__
#define __PDM_BOX_TREE_PRIV_H__

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/
#include "pdm_priv.h"
#include "pdm_box_priv.h"
#include "pdm_morton.h"
/*=============================================================================
 * Macro definitions
 *============================================================================*/

/*=============================================================================
 * Static global variables
 *============================================================================*/

/*============================================================================
 * Type definitions
 *============================================================================*/

/* Structures for each octant or quadrant */
/*----------------------------------------*/

/* If the type is BOX_TREE_NODE, the ordering of children is defined as follows,
   using notation B: bottom, U: up, E: east, W: west, S: south,  N: north.

   octant:   0: BSW, 1: BSE, 2: BNW, 3: BNE, 4: USW, 5: USE, 6: UNW, 7: UNE
   quadrant: 0:  SW, 1:  SE, 2:  NW, 3:  NE
   segment:  0:   W, 1:   E
*/

typedef struct {

  _Bool              is_leaf;      /* True for leaf nodes */

  PDM_morton_code_t  morton_code;  /* Level and coordinates in the grid
                                      according to Morton encoding */

  int   n_boxes;             /* Number of associated bounding boxes */
  int   start_id;            /* Position of the first box_id */
  int   extra_weight;

} _node_t;

/* Structure used to manage statistics */

typedef struct {

  unsigned    max_level_reached;  /* Max level number reached */

  int   n_leaves;           /* Number of leaves in the tree */
  int   n_boxes;            /* Number of boxes to locate in the tree */
  int   n_linked_boxes;     /* Number of linked boxes in the tree */
  int   n_spill_leaves;     /* Number of leaves where n_boxes > threshold */

  int   min_linked_boxes;   /* Minimum number of boxes for a leaf */
  int   max_linked_boxes;   /* Maximum number of boxes for a leaf */

} PDM_box_tree_stats_t;



/* Box tree data */
struct _PDM_box_tree_data_t {

  int        n_max_nodes;     /* Current max. allocated nodes */
  int        n_nodes;         /* Number of nodes (including leaves) */

  _node_t   *nodes;           /* Array of nodes (root at index 0) */

  int       *child_ids;       /* Ids of associated children
         (size: 2^dim * n_max_nodes) */

  int       *box_ids;         /* List of associated box ids.
         size = stat.n_linked_boxes */

  int       *stack;           /* Stack for look for closest leaves */

  int       *pos_stack;       /* Current position in the stack */


  int        n_build_loops;   /* Number of loops required to build */

};

typedef struct {

  PDM_mpi_win_shared_t *w_nodes;
  PDM_mpi_win_shared_t *w_child_ids;
  PDM_mpi_win_shared_t *w_box_ids;

  int                n_max_nodes;
  int                n_nodes;
  int                n_linked_boxes;

} _w_box_tree_data_t;





/* Main box tree structure */
/*-------------------------*/

struct _PDM_box_tree_t {

  PDM_MPI_Comm         comm;         /* Associated MPI communicator */
  int                  n_children;    /* 8, 4, or 2 (2^dim) */

  int                  max_level;     /* Max. possible level */
  int                  threshold;     /* Max number of boxes linked to a
                                         node if max_level is not reached */
  float                max_box_ratio; /* Max n_linked_boxes / n_boxes value */

  PDM_box_tree_stats_t stats;         /* Statistics related to the structure */

  PDM_box_tree_data_t *local_data;    /* Local box tree data */

  int n_copied_ranks;                 /* Number of copies from other ranks */
  int *copied_ranks;                  /* Copied ranks */
  PDM_box_tree_data_t *rank_data;     /* Box tree data copied from other ranks */

  /* Shared memory */
  int                  n_rank_in_shm;
  PDM_box_tree_data_t *shm_data;
  _w_box_tree_data_t  *wbox_tree_data;


  PDM_box_set_t  *boxes;              /* Associated boxes */
};


#endif /* __PDM_BOX_TREE_PRIV_H__ */
