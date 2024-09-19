/*
  This file is part of the ParaDiGM library.

  Copyright (C) 2017       ONERA

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

/*============================================================================
 * Interface structure to represent a distributed mesh
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "pdm.h"
#include "pdm_mpi.h"
#include "pdm_printf.h"
#include "pdm_error.h"
#include "pdm_logging.h"
#include "pdm_part_mesh.h"
#include "pdm_part_mesh_priv.h"
#include "pdm_priv.h"
#include "pdm_writer.h"
#include "pdm_part_connectivity_transform.h"
#include "pdm_extract_part.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "pdm_dmesh.h"

/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */

/*============================================================================
 * Local structure definitions
 *============================================================================*/

/*============================================================================
 * Global variable
 *============================================================================*/

/*=============================================================================
 * Private function definitions
 *============================================================================*/

/* Local bound extractions for \ref PDM_part_mesh_dump_ensight */

static void
_build_extract_part_bound
(
 PDM_part_mesh_t     *pmesh,
 PDM_bound_type_t     bound_type,
 PDM_extract_part_t **extrp
 )
{
  int i_rank;
  PDM_MPI_Comm_rank(pmesh->comm, &i_rank);

  int dim;
  PDM_mesh_entities_t entity_type;
  switch (bound_type) {
    case PDM_BOUND_TYPE_VTX:
      dim = 0;
      entity_type = PDM_MESH_ENTITY_VTX;
      break;

    case PDM_BOUND_TYPE_EDGE:
      dim = 1;
      entity_type = PDM_MESH_ENTITY_EDGE;
      break;

    case PDM_BOUND_TYPE_FACE:
      dim = 2;
      entity_type = PDM_MESH_ENTITY_FACE;
      break;

    case PDM_BOUND_TYPE_CELL:
      dim = 3;
      entity_type = PDM_MESH_ENTITY_CELL;
      break;

    default:
      PDM_error(__FILE__, __LINE__, 0, "Invalid bound_type %d\n", bound_type);
  }

  PDM_UNUSED(entity_type);

  int   n_group          = pmesh->n_group_bnd[bound_type];
  int **group_entity_idx = pmesh->pconcat_bound_idx[bound_type];
  int **group_entity     = pmesh->pconcat_bound    [bound_type];

  int         **selected_l_num = malloc(sizeof(int         *) * pmesh->n_part);
  // PDM_g_num_t **selected_g_num = malloc(sizeof(PDM_g_num_t *) * pmesh->n_part);
  // int         **selected_loc   = malloc(sizeof(int         *) * pmesh->n_part);

  *extrp = PDM_extract_part_create(dim,
                                   pmesh->n_part,
                                   pmesh->n_part,
                                   // PDM_EXTRACT_PART_KIND_FROM_TARGET,
                                   PDM_EXTRACT_PART_KIND_LOCAL,
                                   PDM_SPLIT_DUAL_WITH_HILBERT, // unused
                                   PDM_TRUE,
                                   PDM_OWNERSHIP_KEEP,
                                   pmesh->comm);

  /* Set parts */
  for (int i_part = 0; i_part < pmesh->n_part; i_part++) {

    selected_l_num[i_part] = malloc(sizeof(int        ) * group_entity_idx[i_part][n_group]);
    // selected_g_num[i_part] = malloc(sizeof(PDM_g_num_t) * group_entity_idx[i_part][n_group]);
    // selected_loc  [i_part] = malloc(sizeof(int        ) * group_entity_idx[i_part][n_group] * 3);
    for (int i = 0; i < group_entity_idx[i_part][n_group]; i++) {
      selected_l_num[i_part][i] = group_entity[i_part][i] - 1;
      // selected_g_num[i_part][i] = pmesh->pentity_ln_to_gn[entity_type][i_part][group_entity[i_part][i] - 1];
      // selected_loc  [i_part][3*i  ] = i_rank;
      // selected_loc  [i_part][3*i+1] = i_part;
      // selected_loc  [i_part][3*i+2] = group_entity[i_part][i] - 1;
    }

    // PDM_log_trace_array_int(selected_l_num[i_part],
    //                         group_entity_idx[i_part][n_group],
    //                         "selected_l_num : ");


    PDM_extract_part_selected_lnum_set(*extrp,
                                       i_part,
                                       group_entity_idx[i_part][n_group],
                                       selected_l_num[i_part]);

    // PDM_log_trace_array_long(selected_g_num[i_part],
    //                          group_entity_idx[i_part][n_group],
    //                          "selected_g_num       : ");
    // PDM_extract_part_target_set(*extrp,
    //                             i,
    //                             group_entity_idx[i_part][n_group],
    //                             selected_g_num[i_part],
    //                             selected_loc  [i_part]);

    int n_cell                 = 0;
    int n_face                 = 0;
    int n_edge                 = 0;
    int *cell_face_idx         = NULL;
    int *cell_face             = NULL;
    int *face_edge_idx         = NULL;
    int *face_edge             = NULL;
    int *face_vtx_idx          = NULL;
    int *face_vtx              = NULL;
    int *edge_vtx              = NULL;
    PDM_g_num_t *cell_ln_to_gn = NULL;
    PDM_g_num_t *face_ln_to_gn = NULL;
    PDM_g_num_t *edge_ln_to_gn = NULL;

    if (dim >= 3 && pmesh->pn_entity[PDM_MESH_ENTITY_CELL] != NULL) {
      n_cell        = pmesh->pn_entity[PDM_MESH_ENTITY_CELL][i_part];
      cell_ln_to_gn = pmesh->pentity_ln_to_gn[PDM_MESH_ENTITY_CELL][i_part];

      if (pmesh->pconnectivity_idx[PDM_CONNECTIVITY_TYPE_CELL_FACE] != NULL) {
        cell_face_idx = pmesh->pconnectivity_idx[PDM_CONNECTIVITY_TYPE_CELL_FACE][i_part];
        cell_face     = pmesh->pconnectivity    [PDM_CONNECTIVITY_TYPE_CELL_FACE][i_part];
      }
    }

    if (pmesh->pn_entity[PDM_MESH_ENTITY_FACE] != NULL) {
      n_face        = pmesh->pn_entity[PDM_MESH_ENTITY_FACE][i_part];
      face_ln_to_gn = pmesh->pentity_ln_to_gn[PDM_MESH_ENTITY_FACE][i_part];

      if (pmesh->pconnectivity_idx[PDM_CONNECTIVITY_TYPE_FACE_EDGE] != NULL) {
        face_edge_idx = pmesh->pconnectivity_idx[PDM_CONNECTIVITY_TYPE_FACE_EDGE][i_part];
        face_edge     = pmesh->pconnectivity    [PDM_CONNECTIVITY_TYPE_FACE_EDGE][i_part];
      }

      if (pmesh->pconnectivity_idx[PDM_CONNECTIVITY_TYPE_FACE_VTX] != NULL) {
        face_vtx_idx = pmesh->pconnectivity_idx[PDM_CONNECTIVITY_TYPE_FACE_VTX][i_part];
        face_vtx     = pmesh->pconnectivity    [PDM_CONNECTIVITY_TYPE_FACE_VTX][i_part];
      }
    }

    if (dim >= 1 && pmesh->pn_entity[PDM_MESH_ENTITY_EDGE] != NULL) {
      n_edge        = pmesh->pn_entity[PDM_MESH_ENTITY_EDGE][i_part];
      edge_ln_to_gn = pmesh->pentity_ln_to_gn[PDM_MESH_ENTITY_EDGE][i_part];

      if (pmesh->pconnectivity_idx[PDM_CONNECTIVITY_TYPE_EDGE_VTX] != NULL) {
        edge_vtx = pmesh->pconnectivity[PDM_CONNECTIVITY_TYPE_EDGE_VTX][i_part];
      }
    }

    PDM_extract_part_part_set(*extrp,
                              i_part,
                              n_cell,
                              n_face,
                              n_edge,
                              pmesh->pn_entity[PDM_MESH_ENTITY_VTX][i_part],
                              cell_face_idx,
                              cell_face,
                              face_edge_idx,
                              face_edge,
                              edge_vtx,
                              face_vtx_idx,
                              face_vtx,
                              cell_ln_to_gn,
                              face_ln_to_gn,
                              edge_ln_to_gn,
                              pmesh->pentity_ln_to_gn[PDM_MESH_ENTITY_VTX][i_part],
                              pmesh->vtx_coords[i_part]);

  }

  PDM_extract_part_compute(*extrp);

  for (int i = 0; i < pmesh->n_part; i++) {
    free(selected_l_num[i]);
    // free(selected_g_num[i]);
    // free(selected_loc  [i]);
  }
  free(selected_l_num);
  // free(selected_g_num);
  // free(selected_loc  );
}

/*=============================================================================
 * Public function definitions
 *============================================================================*/

/**
 *
 * \brief Build a distributed mesh structure
 *
 * \param [in]   dn_cell             Number of distributed cells
 * \param [in]   dn_face             Number of distributed faces
 * \param [in]   dn_vtx              Number of distributed vertices
 * \param [in]   n_bnd               Number of boundaries
 * \param [in]   n_join              Number of interfaces with other domains
 *
 * \return     Identifier
 */

PDM_part_mesh_t*
PDM_part_mesh_create
(
 const int             n_part,
       PDM_MPI_Comm    comm
)
{
  PDM_part_mesh_t *pmesh = (PDM_part_mesh_t *) malloc(sizeof(PDM_part_mesh_t));

  pmesh->n_part = n_part;
  pmesh->comm   = comm;

  int tn_part;
  PDM_MPI_Allreduce(&pmesh->n_part, &tn_part, 1, PDM_MPI_INT, PDM_MPI_SUM, comm);
  pmesh->tn_part = tn_part;

  pmesh->pconnectivity          = (int         *** ) malloc( PDM_CONNECTIVITY_TYPE_MAX * sizeof(int         **) );
  pmesh->pconnectivity_idx      = (int         *** ) malloc( PDM_CONNECTIVITY_TYPE_MAX * sizeof(int         **) );
  pmesh->pentity_ln_to_gn       = (PDM_g_num_t *** ) malloc( PDM_MESH_ENTITY_MAX       * sizeof(PDM_g_num_t **) );
  pmesh->pentity_color          = (int         *** ) malloc( PDM_MESH_ENTITY_MAX       * sizeof(int         **) );

  pmesh->pn_entity              = (int          ** ) malloc( PDM_MESH_ENTITY_MAX       * sizeof(int          *) );

  pmesh->is_owner_connectivity  = malloc( PDM_CONNECTIVITY_TYPE_MAX * sizeof(PDM_bool_t   ) );
  pmesh->is_owner_ln_to_gn      = malloc( PDM_MESH_ENTITY_MAX       * sizeof(PDM_bool_t   ) );
  pmesh->is_owner_color         = malloc( PDM_MESH_ENTITY_MAX       * sizeof(PDM_bool_t   ) );

  for(int i = 0; i < PDM_CONNECTIVITY_TYPE_MAX; ++i) {
    pmesh->is_owner_connectivity[i] = PDM_FALSE;
    pmesh->pconnectivity        [i] = NULL;
    pmesh->pconnectivity_idx    [i] = NULL;
  }

  for(int i = 0; i < PDM_MESH_ENTITY_MAX; ++i) {
    pmesh->is_owner_ln_to_gn[i] = PDM_FALSE;
    pmesh->is_owner_color   [i] = PDM_FALSE;
    pmesh->pentity_ln_to_gn [i] = NULL;
    pmesh->pentity_color    [i] = NULL;
    pmesh->pn_entity        [i] = NULL;
  }

  pmesh->pn_bound                = malloc( PDM_BOUND_TYPE_MAX * sizeof(int          **) );
  pmesh->pbound                  = malloc( PDM_BOUND_TYPE_MAX * sizeof(int         ***) );
  pmesh->pbound_ln_to_gn         = malloc( PDM_BOUND_TYPE_MAX * sizeof(PDM_g_num_t ***) );
  pmesh->is_owner_bound          = malloc( PDM_BOUND_TYPE_MAX * sizeof(PDM_bool_t     ) );

  for(int i = 0; i < PDM_BOUND_TYPE_MAX; ++i ) {
    pmesh->n_group_bnd[i] = 0;
  }

  for(int i = 0; i < PDM_BOUND_TYPE_MAX; ++i) {
    pmesh->is_owner_bound         [i] = PDM_FALSE;

    pmesh->pn_bound       [i] = malloc( pmesh->n_part * sizeof(int          *) );
    pmesh->pbound         [i] = malloc( pmesh->n_part * sizeof(int         **) );
    pmesh->pbound_ln_to_gn[i] = malloc( pmesh->n_part * sizeof(PDM_g_num_t **) );
    for(int i_part = 0; i_part < n_part; ++i_part) {
      pmesh->pn_bound               [i][i_part] = NULL;
      pmesh->pbound                 [i][i_part] = NULL;
      pmesh->pbound_ln_to_gn        [i][i_part] = NULL;
    }
  }

  pmesh->pconcat_bound_idx       = malloc( PDM_BOUND_TYPE_MAX * sizeof(int         **) );
  pmesh->pconcat_bound           = malloc( PDM_BOUND_TYPE_MAX * sizeof(int         **) );
  pmesh->pconcat_bound_ln_to_gn  = malloc( PDM_BOUND_TYPE_MAX * sizeof(PDM_g_num_t **) );
  pmesh->is_owner_concat_bound   = malloc( PDM_BOUND_TYPE_MAX * sizeof(PDM_bool_t    ) );
  pmesh->is_compute_concat_bound = malloc( PDM_BOUND_TYPE_MAX * sizeof(PDM_bool_t    ) );
  for(int i = 0; i < PDM_BOUND_TYPE_MAX; ++i) {
    pmesh->is_owner_concat_bound         [i] = PDM_FALSE;
    pmesh->is_compute_concat_bound       [i] = PDM_FALSE;

    pmesh->pconcat_bound_idx     [i] = malloc( pmesh->n_part * sizeof(int         *) );
    pmesh->pconcat_bound         [i] = malloc( pmesh->n_part * sizeof(int         *) );
    pmesh->pconcat_bound_ln_to_gn[i] = malloc( pmesh->n_part * sizeof(PDM_g_num_t *) );
    for(int i_part = 0; i_part < n_part; ++i_part) {
      pmesh->pconcat_bound_idx     [i][i_part] = NULL;
      pmesh->pconcat_bound         [i][i_part] = NULL;
      pmesh->pconcat_bound_ln_to_gn[i][i_part] = NULL;
    }
  }


  pmesh->vtx_coords = malloc(pmesh->n_part * sizeof(double));
  for(int i_part = 0; i_part < n_part; ++i_part) {
    pmesh->vtx_coords[i_part] = NULL;
  }
  pmesh->is_owner_vtx_coord = PDM_FALSE;

  pmesh->ppart_bound_proc_idx = malloc( PDM_BOUND_TYPE_MAX * sizeof(int         **) );
  pmesh->ppart_bound_part_idx = malloc( PDM_BOUND_TYPE_MAX * sizeof(int         **) );
  pmesh->ppart_bound          = malloc( PDM_BOUND_TYPE_MAX * sizeof(int         **) );
  pmesh->is_owner_part_bound  = malloc( PDM_BOUND_TYPE_MAX * sizeof(PDM_bool_t    ) );

  for(int i = 0; i < PDM_BOUND_TYPE_MAX; ++i) {
    pmesh->is_owner_part_bound [i] = PDM_FALSE;
    pmesh->ppart_bound_proc_idx[i] = NULL;
    pmesh->ppart_bound_part_idx[i] = NULL;
    pmesh->ppart_bound         [i] = NULL;
  }

  return pmesh;
}


void
PDM_part_mesh_n_entity_set
(
 PDM_part_mesh_t          *pmesh,
 int                       i_part,
 PDM_mesh_entities_t       entity_type,
 int                       pn_entity
)
{
  if(pmesh->pn_entity[entity_type] == NULL) {
    pmesh->pn_entity[entity_type] = malloc(pmesh->n_part * sizeof(int));
    for(int i = 0; i < pmesh->n_part; ++i) {
      pmesh->pn_entity[entity_type][i] = 0;
    }
  }
  pmesh->pn_entity[entity_type][i_part] = pn_entity;
}


int
PDM_part_mesh_n_entity_get
(
 PDM_part_mesh_t          *pmesh,
 int                       i_part,
 PDM_mesh_entities_t       entity_type
)
{
  if(pmesh->pn_entity[entity_type] != NULL) {
    return pmesh->pn_entity[entity_type][i_part];
  } else {
    return 0;
  }

}


void
PDM_part_mesh_entity_ln_to_gn_set
(
 PDM_part_mesh_t          *pmesh,
 int                       i_part,
 PDM_mesh_entities_t       entity_type,
 PDM_g_num_t              *pentity_ln_to_gn,
 PDM_ownership_t           ownership
)
{
  if(pmesh->pentity_ln_to_gn [entity_type] == NULL) {
    pmesh->pentity_ln_to_gn [entity_type] = malloc(pmesh->n_part * sizeof(PDM_g_num_t *));
    for(int i = 0; i < pmesh->n_part; ++i) {
      pmesh->pentity_ln_to_gn [entity_type][i] = NULL;
    }
  }
  pmesh->pentity_ln_to_gn [entity_type][i_part] = pentity_ln_to_gn;
  if(ownership == PDM_OWNERSHIP_USER || ownership == PDM_OWNERSHIP_UNGET_RESULT_IS_FREE) {
    pmesh->is_owner_ln_to_gn[entity_type] = PDM_FALSE;
  } else {
    pmesh->is_owner_ln_to_gn[entity_type] = PDM_TRUE;
  }
}


void
PDM_part_mesh_entity_ln_to_gn_get
(
 PDM_part_mesh_t          *pmesh,
 int                       i_part,
 PDM_mesh_entities_t       entity_type,
 PDM_g_num_t             **pentity_ln_to_gn,
 PDM_ownership_t           ownership
)
{
  if(pmesh->pentity_ln_to_gn[entity_type] != NULL) {
    *pentity_ln_to_gn = pmesh->pentity_ln_to_gn[entity_type][i_part];
  } else {
    *pentity_ln_to_gn = NULL;
  }
  if(ownership == PDM_OWNERSHIP_USER || ownership == PDM_OWNERSHIP_UNGET_RESULT_IS_FREE) {
    pmesh->is_owner_ln_to_gn[entity_type] = PDM_FALSE;
  } else if (ownership == PDM_OWNERSHIP_KEEP) {
    pmesh->is_owner_ln_to_gn[entity_type] = PDM_TRUE;
  }
}


void
PDM_part_mesh_entity_color_set
(
 PDM_part_mesh_t          *pmesh,
 int                       i_part,
 PDM_mesh_entities_t       entity_type,
 int                      *pentity_color,
 PDM_ownership_t           ownership
)
{
  if(pmesh->pentity_color [entity_type] == NULL) {
    pmesh->pentity_color [entity_type] = malloc(pmesh->n_part * sizeof(int *));
    for(int i = 0; i < pmesh->n_part; ++i) {
      pmesh->pentity_color [entity_type][i] = NULL;
    }
  }
  pmesh->pentity_color [entity_type][i_part] = pentity_color;
  if(ownership == PDM_OWNERSHIP_USER || ownership == PDM_OWNERSHIP_UNGET_RESULT_IS_FREE) {
    pmesh->is_owner_color[entity_type] = PDM_FALSE;
  } else {
    pmesh->is_owner_color[entity_type] = PDM_TRUE;
  }
}

void
PDM_part_mesh_entity_color_get
(
 PDM_part_mesh_t          *pmesh,
 int                       i_part,
 PDM_mesh_entities_t       entity_type,
 int                     **pentity_color,
 PDM_ownership_t           ownership
)
{
  if(pmesh->pentity_color[entity_type] != NULL) {
    *pentity_color = pmesh->pentity_color[entity_type][i_part];
  } else {
    *pentity_color = NULL;
  }
  if(ownership == PDM_OWNERSHIP_USER || ownership == PDM_OWNERSHIP_UNGET_RESULT_IS_FREE) {
    pmesh->is_owner_color[entity_type] = PDM_FALSE;
  } else if (ownership == PDM_OWNERSHIP_KEEP) {
    pmesh->is_owner_color[entity_type] = PDM_TRUE;
  }
}



void
PDM_part_mesh_connectivity_set
(
 PDM_part_mesh_t          *pmesh,
 int                       i_part,
 PDM_connectivity_type_t   connectivity_type,
 int                      *connect,
 int                      *connect_idx,
 PDM_ownership_t           ownership
)
{
  if(pmesh->pconnectivity [connectivity_type] == NULL) {
    pmesh->pconnectivity    [connectivity_type] = malloc(pmesh->n_part * sizeof(int *));
    pmesh->pconnectivity_idx[connectivity_type] = malloc(pmesh->n_part * sizeof(int *));
    for(int i = 0; i < pmesh->n_part; ++i) {
      pmesh->pconnectivity    [connectivity_type][i] = NULL;
      pmesh->pconnectivity_idx[connectivity_type][i] = NULL;
    }
  }
  pmesh->pconnectivity        [connectivity_type][i_part] = connect;
  pmesh->pconnectivity_idx    [connectivity_type][i_part] = connect_idx;

  if(ownership == PDM_OWNERSHIP_USER || ownership == PDM_OWNERSHIP_UNGET_RESULT_IS_FREE) {
    pmesh->is_owner_connectivity[connectivity_type] = PDM_FALSE;
  } else {
    pmesh->is_owner_connectivity[connectivity_type] = PDM_TRUE;
  }
}

void
PDM_part_mesh_vtx_coord_set
(
 PDM_part_mesh_t   *pmesh,
 int                i_part,
 double            *vtx_coord,
 PDM_ownership_t    ownership
)
{
  pmesh->vtx_coords[i_part] = vtx_coord;
  if(ownership == PDM_OWNERSHIP_USER || ownership == PDM_OWNERSHIP_UNGET_RESULT_IS_FREE) {
    pmesh->is_owner_vtx_coord = PDM_FALSE;
  } else {
    pmesh->is_owner_vtx_coord = PDM_TRUE;
  }
}

void
PDM_part_mesh_vtx_coord_get
(
 PDM_part_mesh_t   *pmesh,
 int                i_part,
 double           **vtx_coord,
 PDM_ownership_t    ownership
)
{
  if(pmesh->vtx_coords[i_part] != NULL) {
    *vtx_coord = pmesh->vtx_coords[i_part];
  } else {
    *vtx_coord = NULL;
  }
  if(ownership == PDM_OWNERSHIP_USER || ownership == PDM_OWNERSHIP_UNGET_RESULT_IS_FREE) {
    pmesh->is_owner_vtx_coord = PDM_FALSE;
  } else if (ownership == PDM_OWNERSHIP_KEEP) {
    pmesh->is_owner_vtx_coord = PDM_TRUE;
  }
}

void
PDM_part_mesh_connectivity_get
(
 PDM_part_mesh_t           *pmesh,
 int                        i_part,
 PDM_connectivity_type_t    connectivity_type,
 int                      **connect,
 int                      **connect_idx,
 PDM_ownership_t           ownership
)
{
  if(pmesh->pconnectivity    [connectivity_type] != NULL && pmesh->pconnectivity    [connectivity_type][i_part] != NULL) {
    *connect     = pmesh->pconnectivity    [connectivity_type][i_part];
    *connect_idx = pmesh->pconnectivity_idx[connectivity_type][i_part];
  } else {
    *connect     = NULL;
    *connect_idx = NULL;
  }

  if(ownership == PDM_OWNERSHIP_USER || ownership == PDM_OWNERSHIP_UNGET_RESULT_IS_FREE) {
    pmesh->is_owner_connectivity[connectivity_type] = PDM_FALSE;
  } else if (ownership == PDM_OWNERSHIP_KEEP) {
    pmesh->is_owner_connectivity[connectivity_type] = PDM_TRUE;
  }
}


void
PDM_part_mesh_n_bound_set
(
 PDM_part_mesh_t          *pmesh,
 PDM_bound_type_t          bound_type,
 int                       n_bound
)
{
  assert(pmesh->n_group_bnd[bound_type] == 0);
  pmesh->n_group_bnd[bound_type] = n_bound;

  for(int i_part = 0; i_part < pmesh->n_part; ++i_part) {
    pmesh->pn_bound       [bound_type][i_part] = malloc( n_bound * sizeof(int          ) );
    pmesh->pbound         [bound_type][i_part] = malloc( n_bound * sizeof(int         *) );
    pmesh->pbound_ln_to_gn[bound_type][i_part] = malloc( n_bound * sizeof(PDM_g_num_t *) );
    for(int i_group = 0; i_group < n_bound; ++i_group) {
      pmesh->pn_bound       [bound_type][i_part][i_group] = 0;
      pmesh->pbound         [bound_type][i_part][i_group] = NULL;
      pmesh->pbound_ln_to_gn[bound_type][i_part][i_group] = NULL;
    }
  }
}

int
PDM_part_mesh_n_bound_get
(
 PDM_part_mesh_t          *pmesh,
 PDM_bound_type_t          bound_type
)
{
  return pmesh->n_group_bnd[bound_type];
}

int
PDM_part_mesh_tn_part_get
(
 PDM_part_mesh_t          *pmesh
)
{
  return pmesh->tn_part;
}

int
PDM_part_mesh_n_part_get
(
 PDM_part_mesh_t          *pmesh
)
{
  return pmesh->n_part;
}

void
PDM_part_mesh_bound_set
(
 PDM_part_mesh_t          *pmesh,
 int                       i_part,
 int                       i_group,
 PDM_bound_type_t          bound_type,
 int                       pn_bound,
 int                      *pbound,
 PDM_g_num_t              *pbound_ln_to_gn,
 PDM_ownership_t           ownership
)
{
  pmesh->pn_bound       [bound_type][i_part][i_group] = pn_bound;
  pmesh->pbound         [bound_type][i_part][i_group] = pbound;
  pmesh->pbound_ln_to_gn[bound_type][i_part][i_group] = pbound_ln_to_gn;

  if(ownership == PDM_OWNERSHIP_USER || ownership == PDM_OWNERSHIP_UNGET_RESULT_IS_FREE) {
    pmesh->is_owner_bound[bound_type] = PDM_FALSE;
  } else {
    pmesh->is_owner_bound[bound_type] = PDM_TRUE;
  }
}


void
PDM_part_mesh_bound_get
(
 PDM_part_mesh_t          *pmesh,
 int                       i_part,
 int                       i_group,
 PDM_bound_type_t          bound_type,
 int                      *pn_bound,
 int                     **pbound,
 PDM_g_num_t             **pbound_ln_to_gn,
 PDM_ownership_t           ownership
)
{
  if(pmesh->pbound[bound_type] != NULL) {
    *pn_bound        = pmesh->pn_bound       [bound_type][i_part][i_group];
    *pbound          = pmesh->pbound         [bound_type][i_part][i_group];
    *pbound_ln_to_gn = pmesh->pbound_ln_to_gn[bound_type][i_part][i_group];
  } else {
    *pn_bound        = 0;
    *pbound          = NULL;
    *pbound_ln_to_gn = NULL;
  }

  if(ownership == PDM_OWNERSHIP_USER || ownership == PDM_OWNERSHIP_UNGET_RESULT_IS_FREE) {
    pmesh->is_owner_bound[bound_type] = PDM_FALSE;
  } else if (ownership == PDM_OWNERSHIP_KEEP) {
    pmesh->is_owner_bound[bound_type] = PDM_TRUE;
  }
}


void
PDM_part_mesh_bound_concat_set
(
 PDM_part_mesh_t          *pmesh,
 int                       i_part,
 PDM_bound_type_t          bound_type,
 int                       n_bound,
 int                      *pbound_idx,
 int                      *pbound,
 PDM_g_num_t              *pbound_ln_to_gn,
 PDM_ownership_t           ownership
)
{
  pmesh->n_group_bnd           [bound_type]         = n_bound;
  pmesh->pconcat_bound_idx     [bound_type][i_part] = pbound_idx;
  pmesh->pconcat_bound         [bound_type][i_part] = pbound;
  pmesh->pconcat_bound_ln_to_gn[bound_type][i_part] = pbound_ln_to_gn;

  if(ownership == PDM_OWNERSHIP_USER || ownership == PDM_OWNERSHIP_UNGET_RESULT_IS_FREE) {
    pmesh->is_owner_concat_bound  [bound_type] = PDM_FALSE;
    pmesh->is_compute_concat_bound[bound_type] = PDM_FALSE;
  } else {
    pmesh->is_owner_concat_bound  [bound_type] = PDM_TRUE;
    pmesh->is_compute_concat_bound[bound_type] = PDM_TRUE;
  }

}


void
PDM_part_mesh_bound_concat_compute
(
 PDM_part_mesh_t          *pmesh,
 int                       i_part,
 PDM_bound_type_t          bound_type
)
{
  if(pmesh->is_compute_concat_bound[bound_type] == PDM_FALSE) {
    int n_group = pmesh->n_group_bnd[bound_type];
    assert(pmesh->pconcat_bound_idx[bound_type][i_part] == NULL);
    pmesh->pconcat_bound_idx[bound_type][i_part] = malloc( (n_group+1) * sizeof(int));
    int *_pconcat_bound_idx = pmesh->pconcat_bound_idx[bound_type][i_part];

    _pconcat_bound_idx[0] = 0;
    for(int i_group = 0; i_group < n_group; ++i_group) {
      _pconcat_bound_idx[i_group+1] = _pconcat_bound_idx[i_group] + pmesh->pn_bound[bound_type][i_part][i_group];
    }

    pmesh->pconcat_bound         [bound_type][i_part] = malloc(_pconcat_bound_idx[n_group] * sizeof(int        ));
    pmesh->pconcat_bound_ln_to_gn[bound_type][i_part] = malloc(_pconcat_bound_idx[n_group] * sizeof(PDM_g_num_t));
    int         *_pconcat_bound          = pmesh->pconcat_bound         [bound_type][i_part];
    PDM_g_num_t *_pconcat_bound_ln_to_gn = pmesh->pconcat_bound_ln_to_gn[bound_type][i_part];

    for(int i_group = 0; i_group < n_group; ++i_group) {
      int idx_write = _pconcat_bound_idx[i_group];
      for(int i = 0; i < pmesh->pn_bound[bound_type][i_part][i_group]; ++i) {
        _pconcat_bound         [idx_write+i] = pmesh->pbound         [bound_type][i_part][i_group][i];
        _pconcat_bound_ln_to_gn[idx_write+i] = pmesh->pbound_ln_to_gn[bound_type][i_part][i_group][i];
      }
    }
    pmesh->is_compute_concat_bound[bound_type] = PDM_TRUE;
  }
}

void
PDM_part_mesh_bound_concat_get
(
 PDM_part_mesh_t          *pmesh,
 int                       i_part,
 PDM_bound_type_t          bound_type,
 int                     **pbound_idx,
 int                     **pbound,
 PDM_g_num_t             **pbound_ln_to_gn,
 PDM_ownership_t           ownership
)
{

  *pbound_idx      = pmesh->pconcat_bound_idx     [bound_type][i_part];
  *pbound          = pmesh->pconcat_bound         [bound_type][i_part];
  *pbound_ln_to_gn = pmesh->pconcat_bound_ln_to_gn[bound_type][i_part];

  if(ownership == PDM_OWNERSHIP_USER || ownership == PDM_OWNERSHIP_UNGET_RESULT_IS_FREE) {
    pmesh->is_owner_concat_bound  [bound_type] = PDM_FALSE;
  } else if (ownership == PDM_OWNERSHIP_KEEP) {
    pmesh->is_owner_concat_bound  [bound_type] = PDM_TRUE;
  }
}


void
PDM_part_mesh_part_graph_comm_set
(
 PDM_part_mesh_t          *pmesh,
 int                       i_part,
 PDM_bound_type_t          bound_type,
 int                      *ppart_bound_proc_idx,
 int                      *ppart_bound_part_idx,
 int                      *ppart_bound,
 PDM_ownership_t           ownership
)
{
  if(pmesh->ppart_bound[bound_type] == NULL) {

    pmesh->ppart_bound_proc_idx[bound_type] = malloc( pmesh->n_part * sizeof(int         *) );
    pmesh->ppart_bound_part_idx[bound_type] = malloc( pmesh->n_part * sizeof(int         *) );
    pmesh->ppart_bound         [bound_type] = malloc( pmesh->n_part * sizeof(PDM_g_num_t *) );
    for(int j = 0; j < pmesh->n_part; ++j) {
      pmesh->ppart_bound_proc_idx[bound_type][j] = NULL;
      pmesh->ppart_bound_part_idx[bound_type][j] = NULL;
      pmesh->ppart_bound         [bound_type][j] = NULL;
    }
  }

  pmesh->ppart_bound_proc_idx[bound_type][i_part] = ppart_bound_proc_idx;
  pmesh->ppart_bound_part_idx[bound_type][i_part] = ppart_bound_part_idx;
  pmesh->ppart_bound         [bound_type][i_part] = ppart_bound;

  if(ownership == PDM_OWNERSHIP_USER || ownership == PDM_OWNERSHIP_UNGET_RESULT_IS_FREE) {
    pmesh->is_owner_part_bound[bound_type] = PDM_FALSE;
  } else {
    pmesh->is_owner_part_bound[bound_type] = PDM_TRUE;
  }
}


void
PDM_part_mesh_part_graph_comm_get
(
 PDM_part_mesh_t          *pmesh,
 int                       i_part,
 PDM_bound_type_t          bound_type,
 int                     **ppart_bound_proc_idx,
 int                     **ppart_bound_part_idx,
 int                     **ppart_bound,
 PDM_ownership_t           ownership
)
{
  if(pmesh->ppart_bound[bound_type] != NULL) {
    *ppart_bound_proc_idx = pmesh->ppart_bound_proc_idx[bound_type][i_part];
    *ppart_bound_part_idx = pmesh->ppart_bound_part_idx[bound_type][i_part];
    *ppart_bound          = pmesh->ppart_bound         [bound_type][i_part];
  } else {
    *ppart_bound_proc_idx = NULL;
    *ppart_bound_part_idx = NULL;
    *ppart_bound          = NULL;
  }

  if(ownership == PDM_OWNERSHIP_USER || ownership == PDM_OWNERSHIP_UNGET_RESULT_IS_FREE) {
    pmesh->is_owner_part_bound[bound_type] = PDM_FALSE;
  } else {
    pmesh->is_owner_part_bound[bound_type] = PDM_TRUE;
  }
}

void
PDM_part_mesh_free
(
 PDM_part_mesh_t        *pmesh
)
{
  if (pmesh != NULL) {
    /* Free connectivity */
    for(int i = 0; i < PDM_CONNECTIVITY_TYPE_MAX; ++i) {
      if(pmesh->is_owner_connectivity[i] == PDM_TRUE) {
        for(int i_part = 0; i_part < pmesh->n_part; ++i_part) {
          if(pmesh->pconnectivity[i] != NULL && pmesh->pconnectivity[i][i_part] != NULL) {
            free(pmesh->pconnectivity[i][i_part]);
          }
          if(pmesh->pconnectivity_idx[i] != NULL && pmesh->pconnectivity_idx[i][i_part] != NULL) {
            free(pmesh->pconnectivity_idx[i][i_part]);
          }
        }
      }

      if(pmesh->pconnectivity[i] != NULL) {
        free(pmesh->pconnectivity[i]);
        pmesh->pconnectivity[i] = NULL;
      }

      if(pmesh->pconnectivity_idx[i] != NULL) {
        free(pmesh->pconnectivity_idx[i]);
        pmesh->pconnectivity_idx[i] = NULL;
      }
    }

    /* Free ln_to_gn */
    for(int i = 0; i < PDM_MESH_ENTITY_MAX; ++i) {
      if(pmesh->is_owner_ln_to_gn[i] == PDM_TRUE) {
        for(int i_part = 0; i_part < pmesh->n_part; ++i_part) {
          if(pmesh->pentity_ln_to_gn[i] != NULL && pmesh->pentity_ln_to_gn[i][i_part] != NULL) {
            free(pmesh->pentity_ln_to_gn[i][i_part]);
          }
        }
      }

      if(pmesh->pentity_ln_to_gn[i] != NULL) {
        free(pmesh->pentity_ln_to_gn[i]);
        pmesh->pentity_ln_to_gn[i] = NULL;
      }
    }

    /* Free ln_to_gn */
    for(int i = 0; i < PDM_MESH_ENTITY_MAX; ++i) {
      if(pmesh->is_owner_color[i] == PDM_TRUE) {
        for(int i_part = 0; i_part < pmesh->n_part; ++i_part) {
          if(pmesh->pentity_color[i][i_part] != NULL) {
            free(pmesh->pentity_color[i][i_part]);
          }
        }
      }

      if(pmesh->pentity_color[i] != NULL) {
        free(pmesh->pentity_color[i]);
        pmesh->pentity_color[i] = NULL;
      }
    }

    /* Free vtx__coord */
    if(pmesh->is_owner_vtx_coord == PDM_TRUE) {
      for(int i_part = 0; i_part < pmesh->n_part; ++i_part) {
        free(pmesh->vtx_coords[i_part]);
        pmesh->vtx_coords[i_part] = NULL;
      }
    }
    free(pmesh->vtx_coords);

    /* Free group */
    for(int i = 0; i < PDM_BOUND_TYPE_MAX; ++i) {
      if(pmesh->is_owner_bound[i] == PDM_TRUE) {
        for(int i_part = 0; i_part < pmesh->n_part; ++i_part) {
          for(int i_group = 0; i_group < pmesh->n_group_bnd[i]; ++i_group) {
            if(pmesh->pbound[i][i_part][i_group] != NULL) {
              free(pmesh->pbound[i][i_part][i_group]);
            }
            if(pmesh->pbound_ln_to_gn[i][i_part][i_group] != NULL) {
              free(pmesh->pbound_ln_to_gn[i][i_part][i_group]);
            }
          }

          if(pmesh->pn_bound[i][i_part] != NULL) {
            free(pmesh->pn_bound[i][i_part]);
          }
          if(pmesh->pbound[i][i_part] != NULL) {
            free(pmesh->pbound[i][i_part]);
          }
          if(pmesh->pbound_ln_to_gn[i][i_part] != NULL) {
            free(pmesh->pbound_ln_to_gn[i][i_part]);
          }
        }
      }

      if(pmesh->pn_bound[i] != NULL) {
        free(pmesh->pn_bound[i]);
        pmesh->pn_bound[i] = NULL;
      }

      if(pmesh->pbound[i] != NULL) {
        free(pmesh->pbound[i]);
        pmesh->pbound[i] = NULL;
      }

      if(pmesh->pbound_ln_to_gn[i] != NULL) {
        free(pmesh->pbound_ln_to_gn[i]);
        pmesh->pbound_ln_to_gn[i] = NULL;
      }
    }


    /* Free group */
    for(int i = 0; i < PDM_BOUND_TYPE_MAX; ++i) {
      if(pmesh->is_owner_concat_bound[i] == PDM_TRUE) {
        for(int i_part = 0; i_part < pmesh->n_part; ++i_part) {
          if(pmesh->pconcat_bound_idx[i][i_part] != NULL) {
            free(pmesh->pconcat_bound_idx[i][i_part]);
          }
          if(pmesh->pconcat_bound[i][i_part] != NULL) {
            free(pmesh->pconcat_bound[i][i_part]);
          }
          if(pmesh->pconcat_bound_ln_to_gn[i][i_part] != NULL) {
            free(pmesh->pconcat_bound_ln_to_gn[i][i_part]);
          }
        }
      }

      if(pmesh->pconcat_bound_idx[i] != NULL) {
        free(pmesh->pconcat_bound_idx[i]);
        pmesh->pconcat_bound_idx[i] = NULL;
      }

      if(pmesh->pconcat_bound[i] != NULL) {
        free(pmesh->pconcat_bound[i]);
        pmesh->pconcat_bound[i] = NULL;
      }

      if(pmesh->pconcat_bound_ln_to_gn[i] != NULL) {
        free(pmesh->pconcat_bound_ln_to_gn[i]);
        pmesh->pconcat_bound_ln_to_gn[i] = NULL;
      }
    }

    /* Free group */
    for(int i = 0; i < PDM_BOUND_TYPE_MAX; ++i) {
      if(pmesh->is_owner_part_bound[i] == PDM_TRUE) {
        for(int i_part = 0; i_part < pmesh->n_part; ++i_part) {
          if(pmesh->ppart_bound_proc_idx[i] != NULL && pmesh->ppart_bound_proc_idx[i][i_part] != NULL) {
            free(pmesh->ppart_bound_proc_idx[i][i_part]);
          }
          if(pmesh->ppart_bound_part_idx[i] != NULL && pmesh->ppart_bound_part_idx[i][i_part] != NULL) {
            free(pmesh->ppart_bound_part_idx[i][i_part]);
          }
          if(pmesh->ppart_bound[i] != NULL && pmesh->ppart_bound[i][i_part] != NULL) {
            free(pmesh->ppart_bound[i][i_part]);
          }
        }
      }

      if(pmesh->ppart_bound_proc_idx[i] != NULL) {
        free(pmesh->ppart_bound_proc_idx[i]);
        pmesh->ppart_bound_proc_idx[i] = NULL;
      }

      if(pmesh->ppart_bound_part_idx[i] != NULL) {
        free(pmesh->ppart_bound_part_idx[i]);
        pmesh->ppart_bound_part_idx[i] = NULL;
      }

      if(pmesh->ppart_bound[i] != NULL) {
        free(pmesh->ppart_bound[i]);
        pmesh->ppart_bound[i] = NULL;
      }
    }


    for(int i = 0; i < PDM_MESH_ENTITY_MAX; ++i) {
      if(pmesh->pn_entity[i] !=NULL){
        free(pmesh->pn_entity[i]);
      }
    }

    free(pmesh->pn_entity);
    free(pmesh->pconnectivity);
    free(pmesh->pconnectivity_idx);
    free(pmesh->pentity_ln_to_gn);
    free(pmesh->pentity_color);
    free(pmesh->is_owner_connectivity);
    free(pmesh->is_owner_ln_to_gn    );
    free(pmesh->is_owner_color       );

    free(pmesh->pn_bound       );
    free(pmesh->pbound         );
    free(pmesh->pbound_ln_to_gn);
    free(pmesh->is_owner_bound );

    free(pmesh->pconcat_bound         );
    free(pmesh->pconcat_bound_ln_to_gn);
    free(pmesh->pconcat_bound_idx);
    free(pmesh->is_owner_concat_bound);
    free(pmesh->is_compute_concat_bound);

    free(pmesh->ppart_bound_proc_idx);
    free(pmesh->ppart_bound_part_idx);
    free(pmesh->ppart_bound         );
    free(pmesh->is_owner_part_bound );

    free(pmesh);
    pmesh = NULL;
  }
}



/**
 * \brief Export a partitioned mesh in Ensight format
 *
 * \param [in] pmesh          Pointer to \ref PDM_part_mesh_t object
 * \param [in] directory      Output directory
 * \param [in] name           Output name
 * \param [in] export_bounds  Option to export bounds
 *
 */

void
PDM_part_mesh_dump_ensight
(
 PDM_part_mesh_t *pmesh,
 const char      *directory,
 const char      *name,
 PDM_bool_t       export_bounds
)
{
  int i_rank;
  PDM_MPI_Comm_rank(pmesh->comm, &i_rank);

  /* Compute mesh highest dimension */
  int mesh_dimension = 2;
  int tn_cell = 0;
  if (pmesh->pn_entity[PDM_MESH_ENTITY_CELL] != NULL) {
    for (int i_part = 0; i_part < pmesh->n_part; i_part++) {
      tn_cell += pmesh->pn_entity[PDM_MESH_ENTITY_CELL][i_part];
    }
  }

  int max_tn_cell = 0;
  PDM_MPI_Allreduce(&tn_cell, &max_tn_cell, 1,
                    PDM_MPI_INT, PDM_MPI_MAX, pmesh->comm);

  if (max_tn_cell > 0) {
    mesh_dimension = 3;
  }



  PDM_writer_t *wrt = PDM_writer_create("Ensight",
                                        PDM_WRITER_FMT_BIN,
                                        PDM_WRITER_TOPO_CST,
                                        PDM_WRITER_OFF,
                                        directory,
                                        name,
                                        pmesh->comm,
                                        PDM_IO_KIND_MPI_SIMPLE,
                                        1.,
                                        NULL);

  /* Define main geometry */
  int id_geom = PDM_writer_geom_create(wrt,
                                       name,
                                       pmesh->n_part);

  /* Define bound geometries and perform extractions */
  int id_geom_bound[PDM_BOUND_TYPE_MAX];
  PDM_extract_part_t *extrp[PDM_BOUND_TYPE_MAX];

  // Idea: move bound_type_name to pdm.h and generalize to
  // PDM_mesh_entities_t, PDM_geometry_kind_t, PDM_Mesh_nodal_elt_t, etc ?
  const char *bound_type_name[PDM_BOUND_TYPE_MAX];
  bound_type_name[PDM_BOUND_TYPE_VTX ] = "corners";
  bound_type_name[PDM_BOUND_TYPE_EDGE] = "ridges";
  bound_type_name[PDM_BOUND_TYPE_FACE] = "surfaces";
  bound_type_name[PDM_BOUND_TYPE_CELL] = "volumes";
  bound_type_name[PDM_BOUND_TYPE_ELMT] = "???";

  for (int bound_type = 0; bound_type < PDM_BOUND_TYPE_MAX; bound_type++) {
    id_geom_bound[bound_type] = -1;
    extrp        [bound_type] = NULL;

    log_trace("%d %s\n", pmesh->n_group_bnd[bound_type], bound_type_name[bound_type]);

    if (export_bounds && pmesh->n_group_bnd[bound_type] > 0) {
      char geom_name[999];
      sprintf(geom_name, "%s_%s", name, bound_type_name[bound_type]);
      id_geom_bound[bound_type] = PDM_writer_geom_create(wrt,
                                                         geom_name,
                                                         pmesh->n_part);

      if (!pmesh->is_compute_concat_bound[bound_type]) {
        pmesh->is_owner_concat_bound[bound_type] = PDM_TRUE;
        for (int i_part = 0; i_part < pmesh->n_part; i_part++) {
          PDM_part_mesh_bound_concat_compute(pmesh, i_part, (PDM_bound_type_t) bound_type);
        }
      }

      _build_extract_part_bound(pmesh,
             (PDM_bound_type_t) bound_type,
                                &extrp[bound_type]);
    }
  }

  /* Define variables */
  int id_var_num_part = PDM_writer_var_create(wrt,
                                              PDM_WRITER_OFF,
                                              PDM_WRITER_VAR_SCALAR,
                                              PDM_WRITER_VAR_ELEMENTS,
                                              "num_part");

  int id_var_bound_id = PDM_writer_var_create(wrt,
                                              PDM_WRITER_OFF,
                                              PDM_WRITER_VAR_SCALAR,
                                              PDM_WRITER_VAR_ELEMENTS,
                                              "bound_id");

  int id_var_bound_type = PDM_writer_var_create(wrt,
                                                PDM_WRITER_OFF,
                                                PDM_WRITER_VAR_SCALAR,
                                                PDM_WRITER_VAR_ELEMENTS,
                                                "bound_type");

  PDM_writer_step_beg(wrt, 0.);

  /* Set coordinates */
  for (int i_part = 0; i_part < pmesh->n_part; i_part++) {
    PDM_writer_geom_coord_set(wrt,
                              id_geom,
                              i_part,
                              pmesh->pn_entity[PDM_MESH_ENTITY_VTX][i_part],
                              pmesh->vtx_coords[i_part],
                              pmesh->pentity_ln_to_gn[PDM_MESH_ENTITY_VTX][i_part],
                              PDM_OWNERSHIP_USER);

    /* Bounds */
    for (int bound_type = 0; bound_type < PDM_BOUND_TYPE_MAX; bound_type++) {
      if (id_geom_bound[bound_type] >= 0) {

        PDM_g_num_t *vtx_ln_to_gn = NULL;
        double      *vtx_coord    = NULL;
        int n_vtx = PDM_extract_part_ln_to_gn_get(extrp[bound_type],
                                                  i_part,
                                                  PDM_MESH_ENTITY_VTX,
                                                  &vtx_ln_to_gn,
                                                  PDM_OWNERSHIP_KEEP);
        PDM_extract_part_vtx_coord_get(extrp[bound_type],
                                       i_part,
                                       &vtx_coord,
                                       PDM_OWNERSHIP_KEEP);


        PDM_writer_geom_coord_set(wrt,
                                  id_geom_bound[bound_type],
                                  i_part,
                                  n_vtx,
                                  vtx_coord,
                                  vtx_ln_to_gn,
                                  PDM_OWNERSHIP_USER);

      }
    }
  }


  /* Set elements */
  int has_face_vtx = (pmesh->pconnectivity[PDM_CONNECTIVITY_TYPE_FACE_VTX] != NULL);
  if (!has_face_vtx) {
    assert(pmesh->pconnectivity[PDM_CONNECTIVITY_TYPE_FACE_EDGE] != NULL);
  }

  int **pface_vtx_idx = malloc(sizeof(int *) * pmesh->n_part);
  int **pface_vtx     = malloc(sizeof(int *) * pmesh->n_part);

  if (mesh_dimension == 3) {
    for (int i_part = 0; i_part < pmesh->n_part; i_part++) {
      if (has_face_vtx) {
        pface_vtx_idx[i_part] = pmesh->pconnectivity_idx[PDM_CONNECTIVITY_TYPE_FACE_VTX][i_part];
        pface_vtx    [i_part] = pmesh->pconnectivity    [PDM_CONNECTIVITY_TYPE_FACE_VTX][i_part];
      }
      else {
        // Compute face->vtx if missing
        pface_vtx_idx[i_part] = pmesh->pconnectivity_idx[PDM_CONNECTIVITY_TYPE_FACE_EDGE][i_part];
        PDM_compute_face_vtx_from_face_and_edge(pmesh->pn_entity[PDM_MESH_ENTITY_FACE][i_part],
                                                pmesh->pconnectivity_idx[PDM_CONNECTIVITY_TYPE_FACE_EDGE][i_part],
                                                pmesh->pconnectivity    [PDM_CONNECTIVITY_TYPE_FACE_EDGE][i_part],
                                                pmesh->pconnectivity    [PDM_CONNECTIVITY_TYPE_EDGE_VTX] [i_part],
                                                &pface_vtx[i_part]);
      }
    }

    for (int i_part = 0; i_part < pmesh->n_part; i_part++) {
      PDM_writer_geom_cell3d_cellface_add(wrt,
                                          id_geom,
                                          i_part,
                                          pmesh->pn_entity[PDM_MESH_ENTITY_CELL][i_part],
                                          pmesh->pn_entity[PDM_MESH_ENTITY_FACE][i_part],
                                          pface_vtx_idx[i_part],
                                          NULL,
                                          pface_vtx[i_part],
                                          pmesh->pconnectivity_idx[PDM_CONNECTIVITY_TYPE_CELL_FACE][i_part],
                                          NULL,
                                          pmesh->pconnectivity[PDM_CONNECTIVITY_TYPE_CELL_FACE][i_part],
                                          pmesh->pentity_ln_to_gn[PDM_MESH_ENTITY_CELL][i_part]);
    }
  }
  else {

    for (int i_part = 0; i_part < pmesh->n_part; i_part++) {
      if (has_face_vtx) {
        PDM_writer_geom_faces_facesom_add(wrt,
                                          id_geom,
                                          i_part,
                                          pmesh->pn_entity[PDM_MESH_ENTITY_FACE][i_part],
                                          pmesh->pconnectivity_idx[PDM_CONNECTIVITY_TYPE_FACE_VTX][i_part],
                                          NULL,
                                          pmesh->pconnectivity[PDM_CONNECTIVITY_TYPE_FACE_VTX][i_part],
                                          pmesh->pentity_ln_to_gn[PDM_MESH_ENTITY_FACE][i_part]);
      }
      else {
        PDM_writer_geom_cell2d_cellface_add(wrt,
                                            id_geom,
                                            i_part,
                                            pmesh->pn_entity[PDM_MESH_ENTITY_FACE][i_part],
                                            pmesh->pn_entity[PDM_MESH_ENTITY_EDGE][i_part],
                                            NULL,
                                            NULL,
                                            pmesh->pconnectivity[PDM_CONNECTIVITY_TYPE_EDGE_VTX][i_part],
                                            pmesh->pconnectivity_idx[PDM_CONNECTIVITY_TYPE_FACE_EDGE][i_part],
                                            NULL,
                                            pmesh->pconnectivity[PDM_CONNECTIVITY_TYPE_FACE_EDGE][i_part],
                                            pmesh->pentity_ln_to_gn[PDM_MESH_ENTITY_FACE][i_part]);
      }
    }
  }

  /* Set bound elements */
  int **vtx_vtx = NULL;
  for (int bound_type = 0; bound_type < PDM_BOUND_TYPE_MAX; bound_type++) {

    if (id_geom_bound[bound_type] < 0) continue;

    int id_block = -1;
    if (bound_type == PDM_BOUND_TYPE_EDGE) {
      id_block = PDM_writer_geom_bloc_add(wrt,
                                          id_geom_bound[bound_type],
                                          PDM_WRITER_BAR2,
                                          PDM_OWNERSHIP_USER);

      for (int i_part = 0; i_part < pmesh->n_part; i_part++) {
        PDM_g_num_t *edge_ln_to_gn = NULL;
        int n_edge = PDM_extract_part_ln_to_gn_get(extrp[bound_type],
                                                   i_part,
                                                   PDM_MESH_ENTITY_EDGE,
                                                   &edge_ln_to_gn,
                                                   PDM_OWNERSHIP_KEEP);

        int *edge_vtx     = NULL;
        int *edge_vtx_idx = NULL;
        PDM_extract_part_connectivity_get(extrp[bound_type],
                                          i_part,
                                          PDM_CONNECTIVITY_TYPE_EDGE_VTX,
                                          &edge_vtx,
                                          &edge_vtx_idx,
                                          PDM_OWNERSHIP_KEEP);

        PDM_writer_geom_bloc_std_set(wrt,
                                     id_geom_bound[bound_type],
                                     id_block,
                                     i_part,
                                     n_edge,
                                     edge_vtx,
                                     edge_ln_to_gn);
      }
    }

    else if (bound_type == PDM_BOUND_TYPE_FACE) {

      for (int i_part = 0; i_part < pmesh->n_part; i_part++) {
        PDM_g_num_t *face_ln_to_gn = NULL;
        int n_face = PDM_extract_part_ln_to_gn_get(extrp[bound_type],
                                                   i_part,
                                                   PDM_MESH_ENTITY_FACE,
                                                   &face_ln_to_gn,
                                                   PDM_OWNERSHIP_KEEP);

        int *edge_vtx     = NULL;
        int *edge_vtx_idx = NULL;
        int n_edge = PDM_extract_part_connectivity_get(extrp[bound_type],
                                                       i_part,
                                                       PDM_CONNECTIVITY_TYPE_EDGE_VTX,
                                                       &edge_vtx,
                                                       &edge_vtx_idx,
                                                       PDM_OWNERSHIP_KEEP);

        int *face_vtx     = NULL;
        int *face_vtx_idx = NULL;
        PDM_extract_part_connectivity_get(extrp[bound_type],
                                          i_part,
                                          PDM_CONNECTIVITY_TYPE_FACE_VTX,
                                          &face_vtx,
                                          &face_vtx_idx,
                                          PDM_OWNERSHIP_KEEP);

        int *face_edge     = NULL;
        int *face_edge_idx = NULL;
        PDM_extract_part_connectivity_get(extrp[bound_type],
                                          i_part,
                                          PDM_CONNECTIVITY_TYPE_FACE_EDGE,
                                          &face_edge,
                                          &face_edge_idx,
                                          PDM_OWNERSHIP_KEEP);

        if (face_vtx_idx != NULL) {
          PDM_writer_geom_faces_facesom_add(wrt,
                                            id_geom_bound[bound_type],
                                            i_part,
                                            n_face,
                                            face_vtx_idx,
                                            NULL,
                                            face_vtx,
                                            face_ln_to_gn);
        }
        else {
          assert(face_edge_idx != NULL);
          PDM_writer_geom_cell2d_cellface_add(wrt,
                                              id_geom_bound[bound_type],
                                              i_part,
                                              n_face,
                                              n_edge,
                                              NULL,
                                              NULL,
                                              edge_vtx,
                                              face_edge_idx,
                                              NULL,
                                              face_edge,
                                              face_ln_to_gn);
        }
      }

    }

    else if (bound_type == PDM_BOUND_TYPE_VTX) {

      id_block = PDM_writer_geom_bloc_add(wrt,
                                          id_geom_bound[bound_type],
                                          PDM_WRITER_POINT,
                                          PDM_OWNERSHIP_USER);

      vtx_vtx = malloc(sizeof(int *) * pmesh->n_part);

      for (int i_part = 0; i_part < pmesh->n_part; i_part++) {

        PDM_g_num_t *vtx_ln_to_gn = NULL;
        int n_vtx = PDM_extract_part_ln_to_gn_get(extrp[bound_type],
                                                  i_part,
                                                  PDM_MESH_ENTITY_VTX,
                                                  &vtx_ln_to_gn,
                                                  PDM_OWNERSHIP_KEEP);

        vtx_vtx[i_part] = malloc(sizeof(int) * n_vtx);
        for (int i = 0; i < n_vtx; i++) {
          vtx_vtx[i_part][i] = i + 1;
        }

        PDM_writer_geom_bloc_std_set(wrt,
                                     id_geom_bound[bound_type],
                                     id_block,
                                     i_part,
                                     n_vtx,
                                     vtx_vtx[i_part],
                                     vtx_ln_to_gn);
      }

    }
  }


  /* Write geometries */
  PDM_writer_geom_write(wrt,
                        id_geom);

  for (int i = 0; i < PDM_BOUND_TYPE_MAX; i++) {
    if (id_geom_bound[i] >= 0) {
      PDM_writer_geom_write(wrt,
                            id_geom_bound[i]);
    }
  }


  /* Set variables */
  PDM_real_t **val_num_part   = malloc(sizeof(PDM_real_t *) * pmesh->n_part);
  PDM_real_t **val_bound_id   = malloc(sizeof(PDM_real_t *) * pmesh->n_part);
  PDM_real_t **val_bound_type = malloc(sizeof(PDM_real_t *) * pmesh->n_part);
  for (int i_part = 0; i_part < pmesh->n_part; i_part++) {
    int n_entity = 0;
    if (mesh_dimension == 2) {
      n_entity = pmesh->pn_entity[PDM_MESH_ENTITY_FACE][i_part];
    }
    else {
      n_entity = pmesh->pn_entity[PDM_MESH_ENTITY_CELL][i_part];
    }

    val_num_part  [i_part] = malloc(sizeof(PDM_real_t) * n_entity);
    val_bound_id  [i_part] = malloc(sizeof(PDM_real_t) * n_entity);
    val_bound_type[i_part] = malloc(sizeof(PDM_real_t) * n_entity);
    for (int i = 0; i < n_entity; i++) {
      val_num_part  [i_part][i] = i_rank * pmesh->n_part + i_part; // !! works only if each rank has the same nb of partitions
      val_bound_id  [i_part][i] = 0;
      val_bound_type[i_part][i] = -1;
    }

    PDM_writer_var_set(wrt,
                       id_var_num_part,
                       id_geom,
                       i_part,
                       val_num_part[i_part]);
    free(val_num_part[i_part]);

    PDM_writer_var_set(wrt,
                       id_var_bound_id,
                       id_geom,
                       i_part,
                       val_bound_id[i_part]);
    free(val_bound_id[i_part]);

    PDM_writer_var_set(wrt,
                       id_var_bound_type,
                       id_geom,
                       i_part,
                       val_bound_type[i_part]);
    free(val_bound_type[i_part]);
  }

  for (int bound_type = 0; bound_type < PDM_BOUND_TYPE_MAX; bound_type++) {
    if (id_geom_bound[bound_type] >= 0) {

      for (int i_part = 0; i_part < pmesh->n_part; i_part++) {
        int n_entity = 0;
        if (bound_type == PDM_BOUND_TYPE_EDGE) {
          n_entity = PDM_extract_part_n_entity_get(extrp[bound_type],
                                                   i_part,
                                                   PDM_MESH_ENTITY_EDGE);
        }
        else if (bound_type == PDM_BOUND_TYPE_FACE) {
          n_entity = PDM_extract_part_n_entity_get(extrp[bound_type],
                                                   i_part,
                                                   PDM_MESH_ENTITY_FACE);
        }
        else if (bound_type == PDM_BOUND_TYPE_VTX) {
          n_entity = PDM_extract_part_n_entity_get(extrp[bound_type],
                                                   i_part,
                                                   PDM_MESH_ENTITY_VTX);
        }

        val_num_part  [i_part] = malloc(sizeof(PDM_real_t) * n_entity);
        val_bound_type[i_part] = malloc(sizeof(PDM_real_t) * n_entity);
        for (int k = 0; k < n_entity; k++) {
          val_num_part  [i_part][k] = i_rank * pmesh->n_part + i_part; // !! works only if each rank has the same nb of partitions
          val_bound_type[i_part][k] = bound_type;
        }

        val_bound_id[i_part] = malloc(sizeof(PDM_real_t) * n_entity);
        for (int bound_id = 0; bound_id < pmesh->n_group_bnd[bound_type]; bound_id++) {
          for (int k = pmesh->pconcat_bound_idx[bound_type][i_part][bound_id]; k < pmesh->pconcat_bound_idx[bound_type][i_part][bound_id+1]; k++) {
            val_bound_id[i_part][k] = bound_id+1; // !! assume reequilibrate 'local' yields same order as selected_l_num
          }
        }

        PDM_writer_var_set(wrt,
                           id_var_num_part,
                           id_geom_bound[bound_type],
                           i_part,
                           val_num_part[i_part]);
        free(val_num_part[i_part]);

        PDM_writer_var_set(wrt,
                           id_var_bound_id,
                           id_geom_bound[bound_type],
                           i_part,
                           val_bound_id[i_part]);
        free(val_bound_id[i_part]);

        PDM_writer_var_set(wrt,
                           id_var_bound_type,
                           id_geom_bound[bound_type],
                           i_part,
                           val_bound_type[i_part]);
        free(val_bound_type[i_part]);
      }
    }
  }

  /* Write variables */
  PDM_writer_var_write(wrt,
                       id_var_num_part);

  PDM_writer_var_write(wrt,
                       id_var_bound_id);

  PDM_writer_var_write(wrt,
                       id_var_bound_type);


  PDM_writer_step_end(wrt);

  PDM_writer_free(wrt);


  for (int i = 0; i < PDM_BOUND_TYPE_MAX; i++) {
    if (extrp[i] != NULL) {
      PDM_extract_part_free(extrp[i]);
    }
  }

  if (mesh_dimension == 3 && !has_face_vtx) {
    for (int i = 0; i < pmesh->n_part; i++) {
      free(pface_vtx[i]);
    }
  }
  free(pface_vtx_idx );
  free(pface_vtx     );
  free(val_num_part  );
  free(val_bound_id  );
  free(val_bound_type);

  if (vtx_vtx != NULL) {
    for (int i = 0; i < pmesh->n_part; i++) {
      free(vtx_vtx[i]);
    }
    free(vtx_vtx);
  }
}
