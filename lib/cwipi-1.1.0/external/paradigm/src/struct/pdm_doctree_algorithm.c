/*----------------------------------------------------------------------------
 *  System headers
 *----------------------------------------------------------------------------*/
#include <sys/resource.h>

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <string.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "pdm.h"
#include "pdm_priv.h"
#include "pdm_mpi.h"
#include "pdm_printf.h"
#include "pdm_error.h"
#include "pdm_doctree_priv.h"
#include "pdm_doctree.h"
#include "pdm_timer.h"
#include "pdm_part_to_block.h"
#include "pdm_block_to_part.h"
#include "pdm_sort.h"
#include "pdm_box.h"
#include "pdm_box_tree.h"
#include "pdm_box_priv.h"
#include "pdm_array.h"
#include "pdm_logging.h"
#include "pdm_distrib.h"
#include "pdm_binary_search.h"
#include "pdm_octree_seq.h"
#include "pdm_vtk.h"


#ifdef __cplusplus
extern "C"
#if 0
} /* Fake brace to force back Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */

/*============================================================================
 * Local macro definitions
 *============================================================================*/

/*============================================================================
 * Type definitions
 *============================================================================*/



/*=============================================================================
 * Public function definitions
 *============================================================================*/

// void
// PDM_doctree_points_inside_boxes
// (
//        PDM_doctree_t      *doct,
//  const int                 n_boxes,
//  const double             *box_extents,
//  const PDM_g_num_t        *box_g_num,
//        int               **pts_in_box_idx,
//        PDM_g_num_t       **pts_in_box_g_num,
//        double            **pts_in_box_coord
// )
// {

//   PDM_UNUSED(doct);
//   PDM_UNUSED(n_boxes);
//   PDM_UNUSED(box_extents);
//   PDM_UNUSED(box_g_num);
//   PDM_UNUSED(pts_in_box_idx);
//   PDM_UNUSED(pts_in_box_g_num);
//   PDM_UNUSED(pts_in_box_coord);

// }







#ifdef __cplusplus
}
#endif /* __cplusplus */
