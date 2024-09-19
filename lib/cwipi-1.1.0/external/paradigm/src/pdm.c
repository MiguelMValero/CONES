
/*----------------------------------------------------------------------------
 *  Standar headers
 *----------------------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "pdm.h"
#include "pdm_part.h"
#include "pdm_part_renum.h"
#include "pdm_part_coarse_mesh.h"


#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force back Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */

/*============================================================================
 * Fortran function header
 *============================================================================*/

/*============================================================================
 * Local macro definitions
 *============================================================================*/

/*============================================================================
 * Type definitions
 *============================================================================*/

/*============================================================================
 * Global variable
 *============================================================================*/

/*============================================================================
 * Private function definitions
 *============================================================================*/


/*=============================================================================
 * Public function definitions
 *============================================================================*/

PDM_mesh_entities_t
PDM_connectivity_type_to_entity_type
(
 PDM_connectivity_type_t   connectivity_type
)
{
  PDM_mesh_entities_t mesh_entity = PDM_MESH_ENTITY_MAX;
  if( connectivity_type == PDM_CONNECTIVITY_TYPE_CELL_ELMT ||
      connectivity_type == PDM_CONNECTIVITY_TYPE_CELL_CELL ||
      connectivity_type == PDM_CONNECTIVITY_TYPE_CELL_FACE ||
      connectivity_type == PDM_CONNECTIVITY_TYPE_CELL_EDGE ||
      connectivity_type == PDM_CONNECTIVITY_TYPE_CELL_VTX)
  {
    mesh_entity = PDM_MESH_ENTITY_CELL;
  } else if( connectivity_type == PDM_CONNECTIVITY_TYPE_FACE_ELMT ||
             connectivity_type == PDM_CONNECTIVITY_TYPE_FACE_CELL ||
             connectivity_type == PDM_CONNECTIVITY_TYPE_FACE_FACE ||
             connectivity_type == PDM_CONNECTIVITY_TYPE_FACE_EDGE ||
             connectivity_type == PDM_CONNECTIVITY_TYPE_FACE_VTX )
  {
    mesh_entity = PDM_MESH_ENTITY_FACE;
  } else if( connectivity_type == PDM_CONNECTIVITY_TYPE_EDGE_ELMT ||
             connectivity_type == PDM_CONNECTIVITY_TYPE_EDGE_CELL ||
             connectivity_type == PDM_CONNECTIVITY_TYPE_EDGE_FACE ||
             connectivity_type == PDM_CONNECTIVITY_TYPE_EDGE_EDGE ||
             connectivity_type == PDM_CONNECTIVITY_TYPE_EDGE_VTX )
  {
    mesh_entity = PDM_MESH_ENTITY_EDGE;
  } else if( connectivity_type == PDM_CONNECTIVITY_TYPE_VTX_ELMT ||
             connectivity_type == PDM_CONNECTIVITY_TYPE_VTX_CELL ||
             connectivity_type == PDM_CONNECTIVITY_TYPE_VTX_FACE ||
             connectivity_type == PDM_CONNECTIVITY_TYPE_VTX_EDGE ||
             connectivity_type == PDM_CONNECTIVITY_TYPE_VTX_VTX )
  {
    mesh_entity = PDM_MESH_ENTITY_VTX;
  }
  return mesh_entity;
}


/**
 * \brief Finalize PDM
 *
 * This function frees all allocated global variables
 *
 */

void
PDM_Finalize
(
void
)
{

 /**
  *  Free global array inside part_renum
  */

 PDM_part_renum_method_purge();
 PDM_coarse_mesh_method_purge();

}


#ifdef __cplusplus
}
#endif /* __cplusplus */
