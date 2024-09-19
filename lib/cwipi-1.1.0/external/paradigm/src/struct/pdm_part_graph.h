/*
 * \file
 */

#ifndef __PDM_PART_GRAPH_H__
#define __PDM_PART_GRAPH_H__

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


/*============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Perform a cells renumbering with cache blocking (Synchrone)
 *
 * parameters:
 *   part       --> Mesh Partition
 *---------------------------------------------------------------------------*/

void
PDM_part_graph_split
(
 int         method,
 int         n_part,
 _part_t    *part_ini,
 int        *cell_cell_idx,
 int        *cell_cell,
 int        *cell_weight,
 int        *face_weight,
 int       **cell_part
);

/*----------------------------------------------------------------------------
 * Perform a cells renumbering with cache blocking (Asynchrone)
 *
 * parameters:
 *   part       --> Mesh Partition
 *---------------------------------------------------------------------------*/

void
PDM_part_graph_compute_from_face_cell
(
  _part_t        *part_ini,
  int           **cell_cell_idxCompressed,
  int           **cell_cellCompressed
);


/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------
 * Perform a cells renumbering with cache blocking (Synchrone)
 *
 * parameters:
 *   part       --> Mesh Partition
 *---------------------------------------------------------------------------*/

void
PDM_part_graph_split_bis
(
 int         method,
 int         n_part,
 int         graphSize,
 int        *cell_cell_idx,
 int        *cell_cell,
 int        *cell_weight,
 int        *face_weight,
 int       **cell_part
);

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* __PDM_RENUM_CACHE_BLOCKING_H__ */
