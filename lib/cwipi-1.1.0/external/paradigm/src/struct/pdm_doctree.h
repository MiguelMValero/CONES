/*
 * \file
 */

#ifndef PDM_DOCTREE_H
#define PDM_DOCTREE_H

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "pdm.h"
#include "pdm_mpi.h"

/*----------------------------------------------------------------------------*/

#ifdef  __cplusplus
extern "C" {
#if 0
} /* Fake brace */
#endif
#endif

/*============================================================================
 * Macro definitions
 *============================================================================*/

/*============================================================================
 * Type definitions
 *============================================================================*/

typedef struct _pdm_doctree_t PDM_doctree_t;


/*============================================================================
 * Public function definitions
 *============================================================================*/



PDM_doctree_t*
PDM_doctree_create
(
 PDM_MPI_Comm              comm,
 int                       dim,
 int                       n_part_cloud,
 double                   *global_extents,
 PDM_doctree_local_tree_t  local_tree_kind
);

void
PDM_doctree_build
(
 PDM_doctree_t     *doct
);


void
PDM_doctree_point_set
(
 PDM_doctree_t     *doct,
 const int          i_part_cloud,
 const int          n_points,
 const int         *pts_init_location,
 const PDM_g_num_t *pts_g_num,
 const double      *pts_coords
);


void
PDM_doctree_free
(
  PDM_doctree_t   *doct
);

void
PDM_doctree_solicitation_set
(
 PDM_doctree_t             *doct,
 PDM_tree_solicitation_t    solicitation_kind,
 int                        n_part,
 int                       *n_entity,
 int                      **init_location_entity,
 PDM_g_num_t              **entity_gnum,
 double                   **entity_coords
);

void
PDM_doctree_results_in_orig_frame_get
(
 PDM_doctree_t       *doct,
 int                  n_boxes,
 PDM_g_num_t         *box_g_num,
 int                **box_pts_idx,
 PDM_g_num_t        **box_pts,
 double             **pts_coord
);

void
PDM_doctree_results_in_block_frame_get
(
 PDM_doctree_t       *doct,
 int                 *dn_box,
 PDM_g_num_t        **dbox_g_num,
 int                **dbox_pts_n,
 PDM_g_num_t        **dbox_pts,
 double             **pts_coord,
 PDM_ownership_t      ownership
);

void
PDM_doctree_dump_times
(
  PDM_doctree_t   *doct
);

void
PDM_doctree_init_pts_location_get
(
  PDM_doctree_t   *doct,
  int             *n_pts,
  int            **equi_pts_init_location_idx,
  int            **equi_pts_init_location
);


#ifdef  __cplusplus
}
#endif

#endif  /* PDM_DOCTREE_H */
