/*
 * \file
 */

#ifndef __PDM_CUTHILL_H__
#define __PDM_CUTHILL_H__

/*============================================================================
 * Hilbert space-filling curve construction for coordinates.
 *============================================================================*/

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "pdm.h"
#include "pdm_mpi.h"

#include "pdm_part.h"
#include "pdm_part_priv.h"

/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force back Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */

/*============================================================================
 * Macro and type definitions
 *============================================================================*/

/**
 * \enum PDM_rcm_graph_t
 * \brief Type of graph for building RCM numbering
 *
 */

typedef enum {

  PDM_CUTHILL_GR1 = 0,  /*!< Rank 1 Graph */
  PDM_CUTHILL_GR2 = 1,  /*!< Rank 2 Graph */

} PDM_cuthill_graph_t;


/*============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Determine the reverse Cut-Hill Mac-Kee numbering
 *
 * parameters:
 *   n_elm          --> Number of elements to reorder
 *   dual_graph_idx --> Element to element connectivity indexes (size=n_elm+1)
 *   dual_graph     --> Element to element connectivity (size=dual_graph_idx[n_elm])
 *   perm           <-- Array of permutations
 *---------------------------------------------------------------------------*/

void
PDM_cuthill_generate
(
 int                n_elm,
 int               *dual_graph_idx,
 int               *dual_graph,
 int               *perm
);

/*----------------------------------------------------------------------------
 * Compute the bandwidth of the current mesh
 *
 * parameters:
 *   n_elm          --> Number of elements to reorder
 *   dual_graph_idx --> Element to element connectivity indexes (size=n_elm+1)
 *   dual_graph     --> Element to element connectivity (size=dual_graph_idx[n_elm])
 *---------------------------------------------------------------------------*/
int
PDM_cuthill_checkbandwidth
(
 int                n_elm,
 int               *dual_graph_idx,
 int               *dual_graph
);

/*----------------------------------------------------------------------------
 * Compute the bandwidth of the current mesh
 *
 * parameters:
 *   part       --> Mesh Partition
 *---------------------------------------------------------------------------*/
void
PDM_genrcm
(
int node_num,
int adj_row[],
int adj[],
int perm[]
);


/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* __PDM_CUTHILL_H__ */
