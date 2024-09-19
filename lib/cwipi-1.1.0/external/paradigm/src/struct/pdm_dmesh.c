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
#include "pdm_dmesh_priv.h"
#include "pdm_priv.h"

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

/*=============================================================================
 * Public function definitions
 *============================================================================*/

/**
 *
 * \brief Build a distributed mesh structure
 *
 * \param [in]   dn_cell   Number of distributed cells
 * \param [in]   dn_face   Number of distributed faces
 * \param [in]   dn_edge   Number of distributed edges
 * \param [in]   dn_vtx    Number of distributed vertices
 * \param [in]   comm      PDM_MPI communicator
 *
 * \return     Pointer to a new \ref PDM_dmesh_t object
 */
PDM_dmesh_t*
PDM_dmesh_create
(
       PDM_ownership_t owner,
 const int             dn_cell,
 const int             dn_face,
 const int             dn_edge,
 const int             dn_vtx,
       PDM_MPI_Comm    comm
)
{
  PDM_dmesh_t *dmesh = (PDM_dmesh_t *) malloc(sizeof(PDM_dmesh_t));

  dmesh->comm              = comm;
  dmesh->owner             = owner;
  dmesh->results_is_getted = malloc( PDM_CONNECTIVITY_TYPE_MAX * sizeof(PDM_bool_t) );

  dmesh->dn_cell           = dn_cell;
  dmesh->dn_face           = dn_face;
  dmesh->dn_edge           = dn_edge;
  dmesh->dn_vtx            = dn_vtx;

  dmesh->n_g_cell          = 0;
  dmesh->n_g_face          = 0;
  dmesh->n_g_edge          = 0;
  dmesh->n_g_vtx           = 0;

  PDM_g_num_t _dn_cell = dmesh->dn_cell;
  PDM_g_num_t _dn_face = dmesh->dn_face;
  PDM_g_num_t _dn_edge = dmesh->dn_edge;
  PDM_g_num_t _dn_vtx  = dmesh->dn_vtx;

  PDM_MPI_Allreduce(&_dn_cell, &dmesh->n_g_cell, 1, PDM__PDM_MPI_G_NUM, PDM_MPI_SUM, comm);
  PDM_MPI_Allreduce(&_dn_face, &dmesh->n_g_face, 1, PDM__PDM_MPI_G_NUM, PDM_MPI_SUM, comm);
  PDM_MPI_Allreduce(&_dn_edge, &dmesh->n_g_edge, 1, PDM__PDM_MPI_G_NUM, PDM_MPI_SUM, comm);
  PDM_MPI_Allreduce(&_dn_vtx , &dmesh->n_g_vtx , 1, PDM__PDM_MPI_G_NUM, PDM_MPI_SUM, comm);

  dmesh->cell_distrib      = NULL;
  dmesh->face_distrib      = NULL;
  dmesh->edge_distrib      = NULL;
  dmesh->vtx_distrib       = NULL;

  dmesh->_dvtx_coord       = NULL;
  dmesh->is_owner_vtx_coord  = PDM_TRUE;

  dmesh->dconnectivity         = malloc( PDM_CONNECTIVITY_TYPE_MAX * sizeof(PDM_g_num_t *) );
  dmesh->dconnectivity_idx     = malloc( PDM_CONNECTIVITY_TYPE_MAX * sizeof(int         *) );
  dmesh->is_owner_connectivity = malloc( PDM_CONNECTIVITY_TYPE_MAX * sizeof(PDM_bool_t   ) );

  for(int i = 0; i < PDM_CONNECTIVITY_TYPE_MAX; ++i) {
    dmesh->is_owner_connectivity[i] = PDM_FALSE;
    dmesh->dconnectivity        [i] = NULL;
    dmesh->dconnectivity_idx    [i] = NULL;
  }

  dmesh->dbound          = malloc( PDM_BOUND_TYPE_MAX * sizeof(PDM_g_num_t *) );
  dmesh->dbound_idx      = malloc( PDM_BOUND_TYPE_MAX * sizeof(int         *) );
  dmesh->is_owner_bound  = malloc( PDM_BOUND_TYPE_MAX * sizeof(PDM_bool_t   ) );

  for(int i = 0; i < PDM_BOUND_TYPE_MAX; ++i ) {
    dmesh->n_group_bnd[i] = 0;
  }

  for(int i = 0; i < PDM_BOUND_TYPE_MAX; ++i) {
    dmesh->is_owner_bound[i] = PDM_FALSE;
    dmesh->dbound        [i] = NULL;
    dmesh->dbound_idx    [i] = NULL;
  }

  dmesh->is_computed_g_extents = PDM_FALSE;

  return dmesh;
}

/**
 *
 * \brief Get the dimensions of the distributed mesh
 *
 * \param [in]    dmesh    Pointer to \ref PDM_dmesh_t object
 * \param [out]   dn_cell  Number of distributed cells
 * \param [out]   dn_face  Number of distributed faces
 * \param [out]   dn_edge  Number of distributed edges
 * \param [out]   dn_vtx   Number of distributed vertices
 */
void
PDM_dmesh_dims_get
(
 PDM_dmesh_t *dmesh,
 int         *dn_cell,
 int         *dn_face,
 int         *dn_edge,
 int         *dn_vtx
)
{
  *dn_cell = dmesh->dn_cell;
  *dn_face = dmesh->dn_face;
  *dn_edge = dmesh->dn_edge;
  *dn_vtx  = dmesh->dn_vtx;
}


/**
 *
 * \brief Get the dimensions of the distributed mesh
 *
 * \param [in]    dmesh             Pointer to \ref PDM_dmesh_t object
 * \param [in]    entity_type       Kind of entity
 * \param [in]    dn_cell           Number of distributed cells
 */
void
PDM_dmesh_dn_entity_set
(
 PDM_dmesh_t         *dmesh,
 PDM_mesh_entities_t  entity_type,
 int                  dn_entity
)
{
  if(entity_type == PDM_MESH_ENTITY_CELL) {
    dmesh->dn_cell = dn_entity;
  } else if(entity_type == PDM_MESH_ENTITY_FACE) {
    dmesh->dn_face = dn_entity;
  } else if(entity_type == PDM_MESH_ENTITY_EDGE) {
    dmesh->dn_edge = dn_entity;
  } else if(entity_type == PDM_MESH_ENTITY_VTX) {
    dmesh->dn_vtx = dn_entity;
  }
}

/**
 *
 * \brief Get the dimensions of the distributed mesh
 *
 * \param [in]    dmesh             Pointer to \ref PDM_dmesh_t object
 * \param [in]    entity_type       Kind of mesh entity \ref PDM_mesh_entities_t
 */
int
PDM_dmesh_dn_entity_get
(
 PDM_dmesh_t         *dmesh,
 PDM_mesh_entities_t  entity_type
)
{
  if(entity_type == PDM_MESH_ENTITY_CELL) {
    return dmesh->dn_cell;
  } else if(entity_type == PDM_MESH_ENTITY_FACE) {
    return dmesh->dn_face;
  } else if(entity_type == PDM_MESH_ENTITY_EDGE) {
    return dmesh->dn_edge;
  } else if(entity_type == PDM_MESH_ENTITY_VTX) {
    return dmesh->dn_vtx;
  } else {
    return -1;
  }
}

/**
 *
 * \brief Get the distributed coordinates array
 *
 * \param [in]   dmesh          Pointer to \ref PDM_dmesh_t object
 * \param [out]  dvtx_coord     Vertex coordinate (size = 3 * dn_vtx)
 * \param [in]   ownership      Ownership for color ( \ref PDM_ownership_t )
 */
void
PDM_dmesh_vtx_coord_get
(
 PDM_dmesh_t      *dmesh,
 double          **dvtx_coord,
 PDM_ownership_t   ownership
)
{
  *dvtx_coord      = dmesh->_dvtx_coord;
  if(ownership == PDM_OWNERSHIP_USER || ownership == PDM_OWNERSHIP_UNGET_RESULT_IS_FREE) {
    dmesh->is_owner_vtx_coord = PDM_FALSE;
  } else if(ownership == PDM_OWNERSHIP_KEEP) {
    dmesh->is_owner_vtx_coord = PDM_TRUE;
  }
}

/**
 *
 * \brief Set the distributed coordinates array
 *
 * \param [in]   dmesh          Pointer to \ref PDM_dmesh_t object
 * \param [out]  dvtx_coord     Vertex coordinate (size = 3 * dn_vtx)
 * \param [in]   ownership      Ownership for color ( \ref PDM_ownership_t )
 */
void
PDM_dmesh_vtx_coord_set
(
 PDM_dmesh_t      *dmesh,
 double           *dvtx_coord,
 PDM_ownership_t   ownership
)
{
  dmesh->_dvtx_coord = dvtx_coord;
  if(ownership == PDM_OWNERSHIP_USER || ownership == PDM_OWNERSHIP_UNGET_RESULT_IS_FREE) {
    dmesh->is_owner_vtx_coord = PDM_FALSE;
  } else if(ownership == PDM_OWNERSHIP_KEEP) {
    dmesh->is_owner_vtx_coord = PDM_TRUE;
  }
}


/**
 *
 * \brief Set the distributed connectivity array
 *
 * \param [in]  dmesh             Pointer to \ref PDM_dmesh_t object
 * \param [in]  connectivity_type Connectivity kind \ref PDM_connectivity_type_t
 * \param [in]  connect           Connectivity array (size = connect_idx[n_entity] )
 * \param [in]  connect_idx       Connectivity index (size = n_entity+1 )
 * \param [in]  ownership         Choice of ownership of the input arrays \ref PDM_ownership_t
 */
void
PDM_dmesh_connectivity_set
(
 PDM_dmesh_t              *dmesh,
 PDM_connectivity_type_t   connectivity_type,
 PDM_g_num_t              *connect,
 int                      *connect_idx,
 PDM_ownership_t           ownership
)
{
  assert(dmesh != NULL);
  dmesh->dconnectivity    [connectivity_type] = connect;
  dmesh->dconnectivity_idx[connectivity_type] = connect_idx;

  if(ownership == PDM_OWNERSHIP_USER || ownership == PDM_OWNERSHIP_UNGET_RESULT_IS_FREE) {
    dmesh->is_owner_connectivity[connectivity_type] = PDM_FALSE;
  } else if(ownership == PDM_OWNERSHIP_KEEP) {
    dmesh->is_owner_connectivity[connectivity_type] = PDM_TRUE;
  }
}

/**
 *
 * \brief Get the distributed connectivity array
 *
 * \param [in]  dmesh             Pointer to \ref PDM_dmesh_t object
 * \param [in]  connectivity_type Connectivity kind \ref PDM_connectivity_type_t
 * \param [out] connect           Connectivity array (size = connect_idx[n_entity] )
 * \param [out] connect_idx       Connectivity index (size = n_entity+1 )
 * \param [in]  ownership         Choice of ownership of the input arrays \ref PDM_ownership_t
 * \return Number of element of entity kind
 */
int
PDM_dmesh_connectivity_get
(
 PDM_dmesh_t              *dmesh,
 PDM_connectivity_type_t   connectivity_type,
 PDM_g_num_t             **connect,
 int                     **connect_idx,
 PDM_ownership_t           ownership
)
{
  assert(dmesh != NULL);

  *connect     = dmesh->dconnectivity    [connectivity_type];
  *connect_idx = dmesh->dconnectivity_idx[connectivity_type];

  if(ownership == PDM_OWNERSHIP_USER || ownership == PDM_OWNERSHIP_UNGET_RESULT_IS_FREE) {
    dmesh->is_owner_connectivity[connectivity_type] = PDM_FALSE;
  } else if(ownership == PDM_OWNERSHIP_KEEP) {
    dmesh->is_owner_connectivity[connectivity_type] = PDM_TRUE;
  }

  int dn_entity = -1;
  if( connectivity_type == PDM_CONNECTIVITY_TYPE_CELL_ELMT ||
      connectivity_type == PDM_CONNECTIVITY_TYPE_CELL_CELL ||
      connectivity_type == PDM_CONNECTIVITY_TYPE_CELL_FACE ||
      connectivity_type == PDM_CONNECTIVITY_TYPE_CELL_EDGE ||
      connectivity_type == PDM_CONNECTIVITY_TYPE_CELL_VTX)
  {
    dn_entity = dmesh->dn_cell;
  } else if( connectivity_type == PDM_CONNECTIVITY_TYPE_FACE_ELMT ||
             connectivity_type == PDM_CONNECTIVITY_TYPE_FACE_CELL ||
             connectivity_type == PDM_CONNECTIVITY_TYPE_FACE_FACE ||
             connectivity_type == PDM_CONNECTIVITY_TYPE_FACE_EDGE ||
             connectivity_type == PDM_CONNECTIVITY_TYPE_FACE_VTX )
  {
    dn_entity = dmesh->dn_face;
  } else if( connectivity_type == PDM_CONNECTIVITY_TYPE_EDGE_ELMT ||
             connectivity_type == PDM_CONNECTIVITY_TYPE_EDGE_CELL ||
             connectivity_type == PDM_CONNECTIVITY_TYPE_EDGE_FACE ||
             connectivity_type == PDM_CONNECTIVITY_TYPE_EDGE_EDGE ||
             connectivity_type == PDM_CONNECTIVITY_TYPE_EDGE_VTX )
  {
    dn_entity = dmesh->dn_edge;
  } else if( connectivity_type == PDM_CONNECTIVITY_TYPE_VTX_ELMT ||
             connectivity_type == PDM_CONNECTIVITY_TYPE_VTX_CELL ||
             connectivity_type == PDM_CONNECTIVITY_TYPE_VTX_FACE ||
             connectivity_type == PDM_CONNECTIVITY_TYPE_VTX_EDGE ||
             connectivity_type == PDM_CONNECTIVITY_TYPE_VTX_VTX )
  {
    dn_entity = dmesh->dn_vtx;
  } else if( connectivity_type == PDM_CONNECTIVITY_TYPE_ELMT_CELL ||
             connectivity_type == PDM_CONNECTIVITY_TYPE_ELMT_FACE ||
             connectivity_type == PDM_CONNECTIVITY_TYPE_ELMT_EDGE ||
             connectivity_type == PDM_CONNECTIVITY_TYPE_ELMT_VTX )
  {
    dn_entity = -1;
  }

  return dn_entity;
}



/**
 *
 * \brief Get the distributed connectivity bound array
 *
 * \param [in]  dmesh             Pointer to \ref PDM_dmesh_t object
 * \param [in]  bound_type        Connectivity kind \ref PDM_bound_type_t
 * \param [out] connect           Connectivity array (size = connect_idx[n_bound] )
 * \param [out] connect_idx       Connectivity index (size = n_bound+1 )
 * \param [in]  ownership         Choice of ownership of the input arrays \ref PDM_ownership_t
 * \return Number of group for the requested entity (n_bound)
 */
int
PDM_dmesh_bound_get
(
 PDM_dmesh_t       *dmesh,
 PDM_bound_type_t   bound_type,
 PDM_g_num_t      **connect,
 int              **connect_idx,
 PDM_ownership_t    ownership
)
{
  assert(dmesh != NULL);

  if(ownership == PDM_OWNERSHIP_USER || ownership == PDM_OWNERSHIP_UNGET_RESULT_IS_FREE) {
    dmesh->is_owner_bound[bound_type] = PDM_FALSE;
  } else if(ownership == PDM_OWNERSHIP_KEEP) {
    dmesh->is_owner_bound[bound_type] = PDM_TRUE;
  }

  // assert(dmesh->dbound[bound_type] != NULL);

  *connect     = dmesh->dbound    [bound_type];
  *connect_idx = dmesh->dbound_idx[bound_type];

  return dmesh->n_group_bnd[bound_type];
}


/**
 *
 * \brief Get the distribution of requested entity
 *
 * \param [in]  dmesh             Pointer to \ref PDM_dmesh_t object
 * \param [out] entity_type       entity_type       Kind of mesh entity \ref PDM_mesh_entities_t
 * \param [out] distrib           Distribution array (size = n_rank+1, numbering start at 0)
 * \return Number of process on this distribution ( n_rank )
 */
int
PDM_dmesh_distrib_get
(
 PDM_dmesh_t              *dmesh,
 PDM_mesh_entities_t       entity_type,
 PDM_g_num_t             **distrib
)
{
  switch (entity_type) {
   case PDM_MESH_ENTITY_CELL:
     *distrib = dmesh->cell_distrib;
     break;
   case PDM_MESH_ENTITY_FACE:
     *distrib = dmesh->face_distrib;
     break;
   case PDM_MESH_ENTITY_EDGE:
     *distrib = dmesh->edge_distrib;
     break;
   case PDM_MESH_ENTITY_VTX:
     *distrib = dmesh->vtx_distrib;
     break;
   default:
    PDM_error(__FILE__, __LINE__, 0, "PDM_dmesh_distrib_get invalid entity_type %d\n", entity_type);
    break;
   }
   int n_rank;
   PDM_MPI_Comm_size(dmesh->comm, &n_rank);
   return n_rank;
}


/**
 *
 * \brief Free
 *
 * \param [in]  dmesh             Pointer to \ref PDM_dmesh_t object
 *
 */
void
PDM_dmesh_free
(
 PDM_dmesh_t         *dmesh
)
{
  if (dmesh == NULL) {
    return;
  }
  dmesh->dn_cell           = 0;
  dmesh->dn_face           = 0;
  dmesh->dn_edge           = 0;
  dmesh->dn_vtx            = 0;

  if(dmesh->is_owner_vtx_coord ==  PDM_TRUE) {
    if(dmesh->_dvtx_coord != NULL) {
      free(dmesh->_dvtx_coord);
    }
  }
  dmesh->_dvtx_coord       = NULL;

  // On doit gérer les cas ou la structure est partagé en python et auquel cas
  // On est owner des resultats et il faut free le reste
  // Donc il faut un is_getted + is_owner pour s'en sortir

  if(( dmesh->owner == PDM_OWNERSHIP_KEEP ) ||
     ( dmesh->owner == PDM_OWNERSHIP_UNGET_RESULT_IS_FREE)){
    for(int i = 0; i < PDM_CONNECTIVITY_TYPE_MAX; ++i) {

      if(dmesh->is_owner_connectivity[i] == PDM_TRUE) {

        if(dmesh->dconnectivity[i] != NULL){
          free(dmesh->dconnectivity[i]);
        }
        if(dmesh->dconnectivity_idx[i] != NULL){
          free(dmesh->dconnectivity_idx[i]);
        }
        dmesh->dconnectivity    [i] = NULL;
        dmesh->dconnectivity_idx[i] = NULL;

      }
    }

    for(int i = 0; i < PDM_BOUND_TYPE_MAX; ++i) {

      if(dmesh->is_owner_bound[i] == PDM_TRUE) {

        //printf(" dmesh_free :: %i \n", i);
        if(dmesh->dbound[i] != NULL) {
          free(dmesh->dbound[i]);
        }
        if(dmesh->dbound_idx[i] != NULL){
          free(dmesh->dbound_idx[i]);
        }
        dmesh->dbound    [i] = NULL;
        dmesh->dbound_idx[i] = NULL;

      }
    }
  }

  free(dmesh->results_is_getted    );
  free(dmesh->dconnectivity        );
  free(dmesh->dconnectivity_idx    );
  free(dmesh->is_owner_connectivity);

  free(dmesh->dbound        );
  free(dmesh->dbound_idx    );
  free(dmesh->is_owner_bound);

  /* This result is never getted so we can free them */
  if(dmesh->cell_distrib != NULL) {
    free(dmesh->cell_distrib);
    dmesh->cell_distrib = NULL;
  }

  if(dmesh->face_distrib != NULL) {
    free(dmesh->face_distrib);
    dmesh->face_distrib = NULL;
  }

  if(dmesh->edge_distrib != NULL) {
    free(dmesh->edge_distrib);
    dmesh->edge_distrib = NULL;
  }

  if(dmesh->vtx_distrib != NULL) {
    free(dmesh->vtx_distrib);
    dmesh->vtx_distrib = NULL;
  }

  free (dmesh);
}



/**
 *
 * \brief Compute the bounding box extend of current distributed mesh
 *
 * \param [in]  dmesh             Pointer to \ref PDM_dmesh_t object
 * \return Extents of current mesh (6 components Xmin, Ymin, Zmin, Xmax, Ymax, Zmax )
 *
 */
const double *
PDM_dmesh_global_extents_get
(
 PDM_dmesh_t         *dmesh
 )
{
  if (dmesh->is_computed_g_extents == PDM_FALSE) {

    double l_min[3] = { HUGE_VAL};
    double l_max[3] = {-HUGE_VAL};

    for (int i = 0; i < dmesh->dn_vtx; i++) {
      for (int j = 0; j < 3; j++) {
        double x = dmesh->_dvtx_coord[3*i + j];
        l_min[j] = PDM_MIN(l_min[j], x);
        l_max[j] = PDM_MAX(l_max[j], x);
      }
    }

    PDM_MPI_Allreduce(l_min, dmesh->g_extents,   3, PDM_MPI_DOUBLE, PDM_MPI_MIN, dmesh->comm);
    PDM_MPI_Allreduce(l_max, dmesh->g_extents+3, 3, PDM_MPI_DOUBLE, PDM_MPI_MAX, dmesh->comm);

    dmesh->is_computed_g_extents = PDM_TRUE;
  }

  return dmesh->g_extents;
}


/**
 *
 * \brief Get the distributed connectivity bound array
 *
 * \param [in]  dmesh             Pointer to \ref PDM_dmesh_t object
 * \param [in]  bound_type        Connectivity kind \ref PDM_bound_type_t
 * \param [in]  n_bound           Number of bound for current entity
 * \param [in]  connect           Connectivity array (size = connect_idx[n_bound] )
 * \param [in]  connect_idx       Connectivity index (size = n_bound+1 )
 * \param [in]  ownership         Choice of ownership of the input arrays \ref PDM_ownership_t
 */
void
PDM_dmesh_bound_set
(
 PDM_dmesh_t      *dmesh,
 PDM_bound_type_t  bound_type,
 int               n_bound,
 PDM_g_num_t      *connect,
 int              *connect_idx,
 PDM_ownership_t   ownership
)
{
  assert(dmesh != NULL);

  dmesh->n_group_bnd[bound_type] = n_bound;
  dmesh->dbound     [bound_type] = connect;
  dmesh->dbound_idx [bound_type] = connect_idx;

  if(ownership == PDM_OWNERSHIP_USER || ownership == PDM_OWNERSHIP_UNGET_RESULT_IS_FREE) {
    dmesh->is_owner_bound[bound_type] = PDM_FALSE;
  } else if(ownership == PDM_OWNERSHIP_KEEP) {
    dmesh->is_owner_bound[bound_type] = PDM_TRUE;
  }
}

#ifdef __cplusplus
}
#endif /* __cplusplus */
