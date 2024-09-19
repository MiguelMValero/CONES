/*----------------------------------------------------------------------------
 *  System headers
 *----------------------------------------------------------------------------*/

#include <stdlib.h>
#include <stdio.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "pdm.h"
#include "pdm_config.h"
#include "pdm_mpi.h"
#include "pdm_mpi_ext_dependencies.h"
#include "pdm_printf.h"
#include "pdm_error.h"

/*=============================================================================
 * Macro definitions
 *============================================================================*/

/*============================================================================
 * Type
 *============================================================================*/

/*=============================================================================
 * Static global variables
 *============================================================================*/

/*=============================================================================
 * Static function definitions
 *============================================================================*/

/*=============================================================================
 * Public function definitions
 *============================================================================*/

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force back Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */


#ifdef PDM_HAVE_PARMETIS

int
PDM_ParMETIS_V3_PartKway
(
const PDM_g_num_t *vtxdist,
const PDM_g_num_t *xadj,
const PDM_g_num_t *adjncy,
const int *vwgt,
const int *adjwgt,
const int *wgtflag,
const int *numflag,
const int *ncon,
const int *n_parts,
const double *tpwgts,
const double *ubvec,
const int *edgecut,
int *part,
const PDM_MPI_Comm comm
)
{
  PDM_UNUSED(vtxdist);
  PDM_UNUSED(xadj);
  PDM_UNUSED(adjncy);
  PDM_UNUSED(vwgt);
  PDM_UNUSED(adjwgt);
  PDM_UNUSED(wgtflag);
  PDM_UNUSED(numflag);
  PDM_UNUSED(ncon);
  PDM_UNUSED(n_parts);
  PDM_UNUSED(tpwgts);
  PDM_UNUSED(ubvec);
  PDM_UNUSED(edgecut);
  PDM_UNUSED(part);
  PDM_UNUSED(comm);

  PDM_error(__FILE__, __LINE__, 0,"PDM_ParMETIS_V3_PartKway : Unavailable function with pdm_no_mpi library\n" );
  abort();
  return 0;
}


#endif

#ifdef PDM_HAVE_PTSCOTCH

void
PDM_SCOTCH_dpart
(
const PDM_g_num_t dn_cell,
const PDM_g_num_t *ddual_graph_idx,
const PDM_g_num_t *ddual_graph,
const int *cell_weight,
const int *edgeWeight,
const int check,
const PDM_MPI_Comm comm,
const int n_part,
int *part
)
{
  PDM_UNUSED(dn_cell);
  PDM_UNUSED(ddual_graph_idx);
  PDM_UNUSED(ddual_graph);
  PDM_UNUSED(cell_weight);
  PDM_UNUSED(edgeWeight);
  PDM_UNUSED(check);
  PDM_UNUSED(comm);
  PDM_UNUSED(n_part);
  PDM_UNUSED(part);

  PDM_error(__FILE__, __LINE__, 0,"PDM_SCOTCH_dpart : Unavailable function with pdm_no_mpi library\n" );
  abort();
}


#endif



#ifdef __cplusplus
}
#endif /* __cplusplus */
