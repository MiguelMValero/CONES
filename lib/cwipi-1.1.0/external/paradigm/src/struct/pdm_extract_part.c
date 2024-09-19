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
#include "pdm_priv.h"
#include "pdm_sort.h"
#include "pdm_binary_search.h"
#include "pdm_sort.h"
#include "pdm_part_to_block.h"
#include "pdm_block_to_part.h"
#include "pdm_part_to_part.h"
#include "pdm_timer.h"
#include "pdm_printf.h"
#include "pdm_error.h"
#include "pdm_array.h"
#include "pdm_gnum.h"
#include "pdm_extract_part.h"
#include "pdm_extract_part_priv.h"
#include "pdm_partitioning_algorithm.h"
#include "pdm_dconnectivity_transform.h"
#include "pdm_para_graph_dual.h"
#include "pdm_vtk.h"
#include "pdm_logging.h"
#include "pdm_distrib.h"
#include "pdm_gnum_location.h"
#include "pdm_unique.h"
#include "pdm_part_mesh_nodal_elmts_priv.h"
#include "pdm_ho_ordering.h"

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force back Emacs auto-indentation back to column 0 */
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

static
void
_edge_center_2d
(
  int       n_part_in,
  int      *n_extract,
  int     **extract_lnum,
  int     **pedge_vtx,
  double  **pvtx_coord,
  double ***edge_center
)
{
  for(int i_part = 0; i_part < n_part_in; ++i_part) {
    assert(pedge_vtx     [i_part] != NULL);
    assert(pvtx_coord    [i_part] != NULL);
  }

  double** entity_center = malloc(n_part_in * sizeof(double * ));
  for(int i_part = 0; i_part < n_part_in; ++i_part) {
    entity_center[i_part] = (double *) malloc(3 * n_extract[i_part] * sizeof(double));

    double *_pvtx_coord = pvtx_coord    [i_part];
    int    *_pedge_vtx  = pedge_vtx     [i_part];

    for(int idx_edge = 0; idx_edge < n_extract[i_part]; ++idx_edge) {
      int i_edge = extract_lnum[i_part][idx_edge];
      int i_vtx1 = _pedge_vtx[2*i_edge  ]-1;
      int i_vtx2 = _pedge_vtx[2*i_edge+1]-1;
      entity_center[i_part][3*idx_edge  ] = 0.5 * (_pvtx_coord[3*i_vtx1  ] + _pvtx_coord[3*i_vtx2  ]);
      entity_center[i_part][3*idx_edge+1] = 0.5 * (_pvtx_coord[3*i_vtx1+1] + _pvtx_coord[3*i_vtx2+1]);
      entity_center[i_part][3*idx_edge+2] = 0.5 * (_pvtx_coord[3*i_vtx1+2] + _pvtx_coord[3*i_vtx2+2]);
    }
  }
  *edge_center = entity_center;
}


static
void
_face_center_2d
(
  int       n_part_in,
  int      *n_extract,
  int     **extract_lnum,
  int     **pface_edge_idx,
  int     **pface_edge,
  int     **pedge_vtx,
  double  **pvtx_coord,
  double ***face_center
)
{
  for(int i_part = 0; i_part < n_part_in; ++i_part) {
    assert(pface_edge    [i_part] != NULL);
    assert(pface_edge_idx[i_part] != NULL);
    assert(pedge_vtx     [i_part] != NULL);
    assert(pvtx_coord    [i_part] != NULL);
  }

  double** entity_center = malloc(n_part_in * sizeof(double * ));
  for(int i_part = 0; i_part < n_part_in; ++i_part) {
    entity_center[i_part] = (double *) malloc(3 * n_extract[i_part] * sizeof(double));

    int    *_pface_edge     = pface_edge    [i_part];
    int    *_pface_edge_idx = pface_edge_idx[i_part];
    int    *_pedge_vtx      = pedge_vtx     [i_part];
    double *_pvtx_coord     = pvtx_coord    [i_part];

    for(int idx_face = 0; idx_face < n_extract[i_part]; ++idx_face) {

      int i_face = extract_lnum[i_part][idx_face];
      entity_center[i_part][3*idx_face  ] = 0.;
      entity_center[i_part][3*idx_face+1] = 0.;
      entity_center[i_part][3*idx_face+2] = 0.;

      double inv = 1./((double) _pface_edge_idx[i_face+1] - _pface_edge_idx[i_face]);

      for(int idx_edge = _pface_edge_idx[i_face]; idx_edge < _pface_edge_idx[i_face+1]; ++idx_edge) {
        int i_edge = PDM_ABS(_pface_edge[idx_edge])-1;
        int i_vtx1 = _pedge_vtx[2*i_edge  ] - 1;
        int i_vtx2 = _pedge_vtx[2*i_edge+1] - 1;

        entity_center[i_part][3*idx_face  ] += 0.5 * (_pvtx_coord[3*i_vtx1  ] + _pvtx_coord[3*i_vtx2  ]);
        entity_center[i_part][3*idx_face+1] += 0.5 * (_pvtx_coord[3*i_vtx1+1] + _pvtx_coord[3*i_vtx2+1]);
        entity_center[i_part][3*idx_face+2] += 0.5 * (_pvtx_coord[3*i_vtx1+2] + _pvtx_coord[3*i_vtx2+2]);

      }
      entity_center[i_part][3*idx_face  ] = entity_center[i_part][3*idx_face  ] * inv;
      entity_center[i_part][3*idx_face+1] = entity_center[i_part][3*idx_face+1] * inv;
      entity_center[i_part][3*idx_face+2] = entity_center[i_part][3*idx_face+2] * inv;
    }
  }

  *face_center = entity_center;
}

static
void
_face_center_2d_from_vtx
(
  int       n_part_in,
  int      *n_extract,
  int     **extract_lnum,
  int     **pface_vtx_idx,
  int     **pface_vtx,
  double  **pvtx_coord,
  double ***face_center
)
{
  for(int i_part = 0; i_part < n_part_in; ++i_part) {
    assert(pface_vtx    [i_part] != NULL);
    assert(pface_vtx_idx[i_part] != NULL);
    assert(pvtx_coord   [i_part] != NULL);
  }

  double** entity_center = malloc(n_part_in * sizeof(double * ));
  for(int i_part = 0; i_part < n_part_in; ++i_part) {
    entity_center[i_part] = (double *) malloc(3 * n_extract[i_part] * sizeof(double));

    int    *_pface_vtx     = pface_vtx    [i_part];
    int    *_pface_vtx_idx = pface_vtx_idx[i_part];
    double *_pvtx_coord    = pvtx_coord   [i_part];

    for(int idx_face = 0; idx_face < n_extract[i_part]; ++idx_face) {

      int i_face = extract_lnum[i_part][idx_face];
      entity_center[i_part][3*idx_face  ] = 0.;
      entity_center[i_part][3*idx_face+1] = 0.;
      entity_center[i_part][3*idx_face+2] = 0.;

      double inv = 1./((double) _pface_vtx_idx[i_face+1] - _pface_vtx_idx[i_face]);

      for(int idx_vtx = _pface_vtx_idx[i_face]; idx_vtx < _pface_vtx_idx[i_face+1]; ++idx_vtx) {
        int i_vtx = _pface_vtx[idx_vtx] - 1;

        entity_center[i_part][3*idx_face  ] += _pvtx_coord[3*i_vtx  ];
        entity_center[i_part][3*idx_face+1] += _pvtx_coord[3*i_vtx+1];
        entity_center[i_part][3*idx_face+2] += _pvtx_coord[3*i_vtx+2];

      }
      entity_center[i_part][3*idx_face  ] = entity_center[i_part][3*idx_face  ] * inv;
      entity_center[i_part][3*idx_face+1] = entity_center[i_part][3*idx_face+1] * inv;
      entity_center[i_part][3*idx_face+2] = entity_center[i_part][3*idx_face+2] * inv;
    }
  }

  *face_center = entity_center;
}


static
void
_cell_center_3d
(
  int       n_part_in,
  int      *n_extract,
  int     **extract_lnum,
  int     **pcell_face_idx,
  int     **pcell_face,
  int     **pface_edge_idx,
  int     **pface_edge,
  int     **pface_vtx_idx,
  int     **pface_vtx,
  int     **pedge_vtx,
  double  **pvtx_coord,
  double ***cell_center
)
{
  int from_edge = 0;
  int from_face = 0;
  for(int i_part = 0; i_part < n_part_in; ++i_part) {
    if(pface_edge    [i_part] != NULL) {
      from_edge = 1;
    }
    if(pface_vtx    [i_part] != NULL) {
      from_face = 1;
    }
    assert(pvtx_coord    [i_part] != NULL);
  }

  double** entity_center = malloc(n_part_in * sizeof(double * ));

  if(from_face == 1) {
    for(int i_part = 0; i_part < n_part_in; ++i_part) {
      entity_center[i_part] = (double *) malloc(3 * n_extract[i_part] * sizeof(double));

      int    *_pcell_face     = pcell_face    [i_part];
      int    *_pcell_face_idx = pcell_face_idx[i_part];
      int    *_pface_vtx      = pface_vtx     [i_part];
      int    *_pface_vtx_idx  = pface_vtx_idx [i_part];
      double *_pvtx_coord     = pvtx_coord    [i_part];

      // PDM_log_trace_array_int(extract_lnum[i_part], n_extract[i_part], "extract_lnum ::");
      for(int idx_cell = 0; idx_cell < n_extract[i_part]; ++idx_cell) {
        int i_cell = extract_lnum[i_part][idx_cell];

        entity_center[i_part][3*idx_cell  ] = 0.;
        entity_center[i_part][3*idx_cell+1] = 0.;
        entity_center[i_part][3*idx_cell+2] = 0.;

        double inv = 1./((double)  _pcell_face_idx[idx_cell+1] - _pcell_face_idx[idx_cell]);

        for(int idx_face = _pcell_face_idx[i_cell]; idx_face < _pcell_face_idx[i_cell+1]; ++idx_face) {
          int i_face = PDM_ABS(_pcell_face[idx_face])-1;

          double fcx = 0;
          double fcy = 0;
          double fcz = 0;
          double inv2 = 1./((double)  _pface_vtx_idx[i_face+1] - _pface_vtx_idx[i_face]);

          for(int idx_vtx = _pface_vtx_idx[i_face]; idx_vtx < _pface_vtx_idx[i_face+1]; ++idx_vtx) {
            int i_vtx = _pface_vtx[idx_vtx]-1;
            fcx += _pvtx_coord[3*i_vtx  ];
            fcy += _pvtx_coord[3*i_vtx+1];
            fcz += _pvtx_coord[3*i_vtx+2];
          }
          fcx = fcx * inv2;
          fcy = fcy * inv2;
          fcz = fcz * inv2;

          entity_center[i_part][3*idx_cell  ] += fcx;
          entity_center[i_part][3*idx_cell+1] += fcy;
          entity_center[i_part][3*idx_cell+2] += fcz;
        }

        entity_center[i_part][3*idx_cell  ] = entity_center[i_part][3*idx_cell  ] * inv;
        entity_center[i_part][3*idx_cell+1] = entity_center[i_part][3*idx_cell+1] * inv;
        entity_center[i_part][3*idx_cell+2] = entity_center[i_part][3*idx_cell+2] * inv;
      } /* End cell */
    }
  } else if( from_edge == 1) {
    for(int i_part = 0; i_part < n_part_in; ++i_part) {
      entity_center[i_part] = (double *) malloc(3 * n_extract[i_part] * sizeof(double));

      int    *_pcell_face     = pcell_face    [i_part];
      int    *_pcell_face_idx = pcell_face_idx[i_part];
      int    *_pface_edge     = pface_edge    [i_part];
      int    *_pface_edge_idx = pface_edge_idx[i_part];
      int    *_pedge_vtx      = pedge_vtx     [i_part];
      double *_pvtx_coord     = pvtx_coord    [i_part];

      for(int idx_cell = 0; idx_cell < n_extract[i_part]; ++idx_cell) {
        int i_cell = extract_lnum[i_part][idx_cell];

        entity_center[i_part][3*idx_cell  ] = 0.;
        entity_center[i_part][3*idx_cell+1] = 0.;
        entity_center[i_part][3*idx_cell+2] = 0.;

        double inv = 1./((double)  _pcell_face_idx[idx_cell+1] - _pcell_face_idx[idx_cell]);

        double fcx = 0;
        double fcy = 0;
        double fcz = 0;
        for(int idx_face = _pcell_face_idx[i_cell]; idx_face < _pcell_face_idx[i_cell+1]; ++idx_face) {
          int i_face = PDM_ABS(_pcell_face[idx_face])-1;

          double inv2 = 1./((double)  _pface_edge_idx[i_face+1] - _pface_edge_idx[i_face]);

          for(int idx_edge = _pface_edge_idx[i_face]; idx_edge < _pface_edge_idx[i_face+1]; ++idx_edge) {
            int i_edge = PDM_ABS(_pface_edge[idx_edge])-1;
            int i_vtx1 = _pedge_vtx[2*i_edge  ] - 1;
            int i_vtx2 = _pedge_vtx[2*i_edge+1] - 1;
            fcx += 0.5 * (_pvtx_coord[3*i_vtx1  ] + _pvtx_coord[3*i_vtx2  ]);
            fcy += 0.5 * (_pvtx_coord[3*i_vtx1+1] + _pvtx_coord[3*i_vtx2+1]);
            fcz += 0.5 * (_pvtx_coord[3*i_vtx1+2] + _pvtx_coord[3*i_vtx2+2]);
          }
          fcx = fcx * inv2;
          fcy = fcy * inv2;
          fcz = fcz * inv2;

          entity_center[i_part][3*idx_cell  ] += fcx;
          entity_center[i_part][3*idx_cell+1] += fcy;
          entity_center[i_part][3*idx_cell+2] += fcz;
        }

        entity_center[i_part][3*idx_cell  ] = entity_center[i_part][3*idx_cell  ] * inv;
        entity_center[i_part][3*idx_cell+1] = entity_center[i_part][3*idx_cell+1] * inv;
        entity_center[i_part][3*idx_cell+2] = entity_center[i_part][3*idx_cell+2] * inv;
      } /* End cell */
    }
  }

  *cell_center = entity_center;
}


static
void
_extract_part_group
(
  PDM_extract_part_t        *extrp,
  PDM_bound_type_t           bound_type
)
{
  int i_rank;
  PDM_MPI_Comm_rank(extrp->comm, &i_rank);

  PDM_mesh_entities_t entity_type = PDM_MESH_ENTITY_MAX;
  int* n_entity = 0;
  if(bound_type == PDM_BOUND_TYPE_FACE) {
    entity_type = PDM_MESH_ENTITY_FACE;
    n_entity = extrp->n_face;
  } else if(bound_type == PDM_BOUND_TYPE_EDGE) {
    entity_type = PDM_MESH_ENTITY_EDGE;
    n_entity = extrp->n_edge;
  } else if(bound_type == PDM_BOUND_TYPE_VTX) {
    entity_type = PDM_MESH_ENTITY_VTX;
    n_entity = extrp->n_vtx;
  } else {
    return;
  }

  PDM_part_to_part_t* ptp = extrp->ptp_entity[entity_type];
  if(ptp == NULL) {
    return;
  }

  int n_group = extrp->n_group[bound_type];
  int          **n_group_entity        = extrp->n_group_entity       [bound_type];
  int         ***group_entity          = extrp->group_entity         [bound_type];
  PDM_g_num_t ***group_entity_ln_to_gn = extrp->group_entity_ln_to_gn[bound_type];

  int         **entity_tag           = malloc(extrp->n_part_in * sizeof(int         *));
  int         **entity_send_n        = malloc(extrp->n_part_in * sizeof(int         *));
  int         **entity_init_location = malloc(extrp->n_part_in * sizeof(int         *));
  PDM_g_num_t **entity_ln_to_gn      = malloc(extrp->n_part_in * sizeof(PDM_g_num_t *));

  for(int i_part = 0; i_part < extrp->n_part_in; ++i_part) {

    entity_send_n       [i_part] = malloc(    n_entity[i_part] * sizeof(int        ));

    for(int i_entity = 0; i_entity < n_entity[i_part]; ++i_entity) {
      entity_send_n[i_part][i_entity] = 0;
    }

    /* First loop to count */
    for(int i_group = 0; i_group < n_group; ++i_group) {
      for(int idx_entity = 0; idx_entity < n_group_entity[i_group][i_part]; ++idx_entity) {
        int i_entity = group_entity[i_group][i_part][idx_entity];
        entity_send_n[i_part][i_entity]++;
      }
    }

    int *entity_send_idx = malloc( (n_entity[i_part] +1) * sizeof(int        ));
    entity_send_idx[0] = 0;
    for(int i_entity = 0; i_entity < n_entity[i_part]; ++i_entity) {
      entity_send_idx[i_entity+1] = entity_send_idx[i_entity] + entity_send_n[i_part][i_entity];
      entity_send_n[i_part][i_entity] = 0;
    }

    entity_tag          [i_part] = malloc(    entity_send_idx[n_entity[i_part]] * sizeof(int        ));
    entity_init_location[i_part] = malloc(3 * entity_send_idx[n_entity[i_part]] * sizeof(int        ));
    entity_ln_to_gn     [i_part] = malloc(    entity_send_idx[n_entity[i_part]] * sizeof(PDM_g_num_t));

    for(int i_group = 0; i_group < n_group; ++i_group) {
      for(int idx_entity = 0; idx_entity < n_group_entity[i_group][i_part]; ++idx_entity) {
        int i_entity = group_entity[i_group][i_part][idx_entity];
        int idx_write = entity_send_idx[i_entity] + entity_send_n[i_part][i_entity]++;
        entity_tag          [i_part][  idx_write  ] = i_group;
        entity_init_location[i_part][3*idx_write  ] = i_rank;
        entity_init_location[i_part][3*idx_write+1] = i_part;
        entity_init_location[i_part][3*idx_write+2] = idx_entity;
        entity_ln_to_gn     [i_part][  idx_write  ] = group_entity_ln_to_gn[i_group][i_part][idx_entity];
      }
    }

    free(entity_send_idx);
  }

  int           exch_request = -1;
  int **pextract_entity_tag = NULL;
  int **pextract_entity_n   = NULL;
  PDM_part_to_part_reverse_iexch(ptp,
                                 PDM_MPI_COMM_KIND_P2P,
                                 PDM_STRIDE_VAR_INTERLACED,
                                 PDM_PART_TO_PART_DATA_DEF_ORDER_PART2,
                                 1,
                                 sizeof(int),
                (const int  **)  entity_send_n,
                (const void **)  entity_tag,
                                 &pextract_entity_n,
                    (void ***)   &pextract_entity_tag,
                                 &exch_request);
  PDM_part_to_part_reverse_iexch_wait(ptp, exch_request);

  int **pextract_entity_init_location   = NULL;
  int **pextract_entity_init_location_n = NULL;
  PDM_part_to_part_reverse_iexch(ptp,
                                 PDM_MPI_COMM_KIND_P2P,
                                 PDM_STRIDE_VAR_INTERLACED,
                                 PDM_PART_TO_PART_DATA_DEF_ORDER_PART2,
                                 1,
                                 3 * sizeof(int),
                (const int  **)  entity_send_n,
                (const void **)  entity_init_location,
                                 &pextract_entity_init_location_n,
                    (void ***)   &pextract_entity_init_location,
                                 &exch_request);
  PDM_part_to_part_reverse_iexch_wait(ptp, exch_request);

  PDM_g_num_t **pextract_entity_parent_ln_to_gn   = NULL;
  int         **pextract_entity_parent_ln_to_gn_n = NULL;
  PDM_part_to_part_reverse_iexch(ptp,
                                 PDM_MPI_COMM_KIND_P2P,
                                 PDM_STRIDE_VAR_INTERLACED,
                                 PDM_PART_TO_PART_DATA_DEF_ORDER_PART2,
                                 1,
                                 sizeof(PDM_g_num_t),
                (const int  **)  entity_send_n,
                (const void **)  entity_ln_to_gn,
                                 &pextract_entity_parent_ln_to_gn_n,
                    (void ***)   &pextract_entity_parent_ln_to_gn,
                                 &exch_request);
  PDM_part_to_part_reverse_iexch_wait(ptp, exch_request);

  int          **pextract_n_group_entity               = malloc(n_group * sizeof(int          *));
  int         ***pextract_group_entity                 = malloc(n_group * sizeof(int         **));
  int         ***pextract_group_entity_init_location   = malloc(n_group * sizeof(int         **));
  PDM_g_num_t ***pextract_group_entity_ln_to_gn        = malloc(n_group * sizeof(PDM_g_num_t **));
  PDM_g_num_t ***pextract_group_entity_parent_ln_to_gn = malloc(n_group * sizeof(PDM_g_num_t **));
  // int ***pextract_group_entity   = malloc(n_group * sizeof(int *));

  for(int i_group = 0; i_group < n_group; ++i_group) {
    pextract_n_group_entity              [i_group] = malloc(extrp->n_part_out * sizeof(int          ));
    pextract_group_entity                [i_group] = malloc(extrp->n_part_out * sizeof(int         *));
    pextract_group_entity_init_location  [i_group] = malloc(extrp->n_part_out * sizeof(int         *));
    pextract_group_entity_ln_to_gn       [i_group] = malloc(extrp->n_part_out * sizeof(PDM_g_num_t *));
    pextract_group_entity_parent_ln_to_gn[i_group] = malloc(extrp->n_part_out * sizeof(PDM_g_num_t *));
  }

  for(int i_part = 0; i_part < extrp->n_part_out; ++i_part) {
    int n_entity_out = extrp->pextract_n_entity[entity_type][i_part];

    // for(int i = 0; i < n_entity_out; ++i) {
    //   n_recv += pextract_entity_n[i_part][i];
    // }

    /* On doit trier par group */
    for(int i_group = 0; i_group < n_group; ++i_group) {
      pextract_n_group_entity[i_group][i_part] = 0;
    }

    int idx_read = 0;
    for(int i = 0; i < n_entity_out; ++i) {
      for(int j = 0; j < pextract_entity_n[i_part][i]; ++j) {
        int i_group = pextract_entity_tag[i_part][idx_read++];
        if(i_group > -1) {
          pextract_n_group_entity[i_group][i_part]++;
        }
      }
    }

    for(int i_group = 0; i_group < n_group; ++i_group) {
      pextract_group_entity                [i_group][i_part] = malloc(    pextract_n_group_entity[i_group][i_part] * sizeof(int         *));
      pextract_group_entity_init_location  [i_group][i_part] = malloc(3 * pextract_n_group_entity[i_group][i_part] * sizeof(int         *));
      pextract_group_entity_parent_ln_to_gn[i_group][i_part] = malloc(    pextract_n_group_entity[i_group][i_part] * sizeof(PDM_g_num_t *));
      pextract_n_group_entity              [i_group][i_part] = 0;
    }

    /* Fill */
    idx_read = 0;
    for(int i = 0; i < n_entity_out; ++i) {
      for(int j = 0; j < pextract_entity_n[i_part][i]; ++j) {
        int         i_group  = pextract_entity_tag            [i_part][  idx_read  ];
        int         t_rank   = pextract_entity_init_location  [i_part][3*idx_read  ];
        int         t_part   = pextract_entity_init_location  [i_part][3*idx_read+1];
        int         t_entity = pextract_entity_init_location  [i_part][3*idx_read+2];
        PDM_g_num_t g_num    = pextract_entity_parent_ln_to_gn[i_part][  idx_read  ];

        idx_read++;
        if(i_group > -1) {
          int idx_write = pextract_n_group_entity[i_group][i_part]++;

          pextract_group_entity                [i_group][i_part][  idx_write  ] = i;
          pextract_group_entity_init_location  [i_group][i_part][3*idx_write  ] = t_rank;
          pextract_group_entity_init_location  [i_group][i_part][3*idx_write+1] = t_part;
          pextract_group_entity_init_location  [i_group][i_part][3*idx_write+2] = t_entity;
          pextract_group_entity_parent_ln_to_gn[i_group][i_part][  idx_write  ] = g_num;

        }
      }
    }

    // PDM_log_trace_array_int(pextract_entity_tag[i_part], n_recv, "pextract_entity_tag ::");
    free(pextract_entity_init_location    [i_part]);
    free(pextract_entity_init_location_n  [i_part]);
    free(pextract_entity_n                [i_part]);
    free(pextract_entity_tag              [i_part]);
    free(pextract_entity_parent_ln_to_gn  [i_part]);
    free(pextract_entity_parent_ln_to_gn_n[i_part]);
  }
  free(pextract_entity_init_location  );
  free(pextract_entity_init_location_n);
  free(pextract_entity_parent_ln_to_gn  );
  free(pextract_entity_parent_ln_to_gn_n);
  free(pextract_entity_n  );
  free(pextract_entity_tag);

  for(int i_part = 0; i_part < extrp->n_part_in; ++i_part) {
    free(entity_tag          [i_part]);
    free(entity_send_n       [i_part]);
    free(entity_ln_to_gn     [i_part]);
    free(entity_init_location[i_part]);
  }
  free(entity_tag          );
  free(entity_send_n       );
  free(entity_ln_to_gn     );
  free(entity_init_location);

  /*
   * Prepare output
   */
  assert(extrp->ptp_group_entity                     [bound_type] == NULL);
  assert(extrp->ptp_group_ownership                  [bound_type] == NULL);
  assert(extrp->pn_extract_group_entity              [bound_type] == NULL);
  assert(extrp->pextract_group_entity                [bound_type] == NULL);
  assert(extrp->pextract_group_entity_ln_to_gn       [bound_type] == NULL);
  assert(extrp->pextract_group_entity_parent_ln_to_gn[bound_type] == NULL);
  assert(extrp->group_array_ownership                [bound_type] == NULL);
  assert(extrp->is_owner_extract_group               [bound_type] == NULL);

  extrp->ptp_group_entity      [bound_type] = malloc(n_group * sizeof(PDM_part_to_part_t  *));
  extrp->ptp_group_ownership   [bound_type] = malloc(n_group * sizeof(PDM_ownership_t      ));
  extrp->group_array_ownership [bound_type] = malloc(n_group * sizeof(PDM_ownership_t      ));
  extrp->is_owner_extract_group[bound_type] = malloc(n_group * sizeof(PDM_bool_t           ));

  /* Create all ptp */
  for(int i_group = 0; i_group < n_group; ++i_group) {

    /*
     * Compute child gnum
     */
    PDM_gen_gnum_t* gen_gnum = PDM_gnum_create(3,
                                               extrp->n_part_out,
                                               PDM_TRUE,
                                               1e-3,
                                               extrp->comm,
                                               PDM_OWNERSHIP_USER);

    for(int i_part = 0; i_part < extrp->n_part_out; ++i_part) {
      PDM_gnum_set_from_parents(gen_gnum,
                                i_part,
                                pextract_n_group_entity[i_group][i_part],
                                pextract_group_entity_parent_ln_to_gn[i_group][i_part]);
    }
    PDM_gnum_compute(gen_gnum);


    for(int i_part = 0; i_part < extrp->n_part_out; ++i_part) {
      pextract_group_entity_ln_to_gn[i_group][i_part] = PDM_gnum_get(gen_gnum, i_part);
    }

    PDM_gnum_free(gen_gnum);

    int **part2_entity1_to_part1_entity1_idx = malloc(extrp->n_part_out * sizeof(int *));
    for(int i_part = 0; i_part < extrp->n_part_out; ++i_part) {
      int n_entity_in_group = pextract_n_group_entity[i_group][i_part];
      part2_entity1_to_part1_entity1_idx[i_part] = malloc((n_entity_in_group+1) * sizeof(int));
      for(int i = 0; i < n_entity_in_group+1; ++i ) {
        part2_entity1_to_part1_entity1_idx[i_part][i] = 3 * i;
      }
    }
    PDM_part_to_part_t* ptp_group = PDM_part_to_part_create_from_num2_triplet((const PDM_g_num_t **) pextract_group_entity_ln_to_gn[i_group],
                                                                              (const int          *) pextract_n_group_entity[i_group],
                                                                                                     extrp->n_part_out,
                                                                              (const int          *) extrp->n_group_entity[bound_type][i_group],
                                                                                                     extrp->n_part_in,
                                                                              (const int         **) part2_entity1_to_part1_entity1_idx,
                                                                              (const int         **) NULL,
                                                                              (const int         **) pextract_group_entity_init_location[i_group],
                                                                                                     extrp->comm);
    if(0 == 1) {
      int          *n_ref_entity1     = NULL;
      int         **ref_l_num_entity1 = NULL;
      PDM_part_to_part_ref_lnum2_get(ptp_group, &n_ref_entity1, &ref_l_num_entity1);

      for(int i_part = 0; i_part < extrp->n_part_in; ++i_part) {
        PDM_log_trace_array_int(ref_l_num_entity1[i_part], n_ref_entity1[i_part], "ref_l_num_entity1 ::");
      }
    }

    for(int i_part = 0; i_part < extrp->n_part_out; ++i_part) {
      free(part2_entity1_to_part1_entity1_idx[i_part]);
    }
    free(part2_entity1_to_part1_entity1_idx);

    // PDM_part_to_part_free(ptp_group);

    extrp->ptp_group_ownership   [bound_type][i_group] = PDM_OWNERSHIP_KEEP;
    extrp->group_array_ownership [bound_type][i_group] = PDM_OWNERSHIP_KEEP;
    extrp->is_owner_extract_group[bound_type][i_group] = PDM_TRUE;
    extrp->ptp_group_entity      [bound_type][i_group] = ptp_group;
  }

  extrp->pn_extract_group_entity              [bound_type] = pextract_n_group_entity;
  extrp->pextract_group_entity                [bound_type] = pextract_group_entity;
  extrp->pextract_group_entity_ln_to_gn       [bound_type] = pextract_group_entity_ln_to_gn;
  extrp->pextract_group_entity_parent_ln_to_gn[bound_type] = pextract_group_entity_parent_ln_to_gn;


  for(int i_group = 0; i_group < n_group; ++i_group) {
    for(int i_part = 0; i_part < extrp->n_part_out; ++i_part) {
      // free(pextract_group_entity                [i_group][i_part]);
      free(pextract_group_entity_init_location  [i_group][i_part]);
      // free(pextract_group_entity_ln_to_gn       [i_group][i_part]);
      // free(pextract_group_entity_parent_ln_to_gn[i_group][i_part]);
    }
    // free(pextract_n_group_entity              [i_group]);
    // free(pextract_group_entity                [i_group]);
    free(pextract_group_entity_init_location  [i_group]);
    // free(pextract_group_entity_ln_to_gn       [i_group]);
    // free(pextract_group_entity_parent_ln_to_gn[i_group]);
  }
  // free(pextract_n_group_entity              );
  // free(pextract_group_entity                );
  free(pextract_group_entity_init_location  );
  // free(pextract_group_entity_ln_to_gn       );
  // free(pextract_group_entity_parent_ln_to_gn);


}


static
void
_compute_child
(
 PDM_MPI_Comm   comm,
 int            n_part,
 int           *n_child,
 PDM_g_num_t  **parent_g_num,
 PDM_g_num_t  **child_g_num
)
{
  PDM_gen_gnum_t* gnum_extract = PDM_gnum_create(3, n_part, PDM_FALSE,
                                                 1.e-6,
                                                 comm,
                                                 PDM_OWNERSHIP_USER);
  for(int i_part = 0; i_part < n_part; ++i_part) {
    PDM_gnum_set_from_parents(gnum_extract, i_part, n_child[i_part], parent_g_num[i_part]);
  }
  PDM_gnum_compute(gnum_extract);

  for (int i_part = 0; i_part < n_part; i_part++){
    child_g_num[i_part] = PDM_gnum_get(gnum_extract, i_part);
  }
  PDM_gnum_free(gnum_extract);
}

static
void
_extract_gnum_and_compute_child
(
 PDM_MPI_Comm    comm,
 int             compute_child,
 int             n_part,
 int            *n_extract,
 PDM_g_num_t   **entity_g_num,
 int           **extract_lnum,
 PDM_g_num_t  ***parent_selected_g_num_out,
 PDM_g_num_t  ***child_selected_g_num_out
)
{
  PDM_g_num_t** entity_extract_g_num = (PDM_g_num_t **) malloc( n_part * sizeof(PDM_g_num_t *));

  for(int i_part = 0; i_part < n_part; ++i_part) {
    entity_extract_g_num[i_part] = (PDM_g_num_t *) malloc( n_extract[i_part] * sizeof(PDM_g_num_t));
    for(int i_entity = 0; i_entity < n_extract[i_part]; ++i_entity) {
      // log_trace("extract_lnum[i_part][%d] = %d\n",
      //           i_entity,
      //           extract_lnum[i_part][i_entity]);
      entity_extract_g_num[i_part][i_entity] = entity_g_num[i_part][extract_lnum[i_part][i_entity]];
    }
  }

  PDM_g_num_t **child_selected_g_num = (PDM_g_num_t **) malloc( n_part * sizeof(PDM_g_num_t *));
  for(int i_part = 0; i_part < n_part; ++i_part)  {
    child_selected_g_num[i_part] = NULL;
  }

  if(compute_child == 1) {
    _compute_child(comm,
                   n_part,
                   n_extract,
                   entity_extract_g_num,
                   child_selected_g_num);
  }

  *parent_selected_g_num_out = entity_extract_g_num;
  *child_selected_g_num_out  = child_selected_g_num;
}

static
void
_compute_dual_graph
(
  PDM_extract_part_t    *extrp,
  PDM_part_to_block_t   *ptb_equi,
  PDM_g_num_t         **out_dual_graph_idx,
  PDM_g_num_t         **out_dual_graph
)
{
  int i_rank;
  PDM_MPI_Comm_rank(extrp->comm, &i_rank);

  // int dn_equi = PDM_part_to_block_n_elt_block_get(ptb_equi);
  // PDM_g_num_t *dequi_g_num = PDM_part_to_block_block_gnum_get(ptb_equi);

  PDM_g_num_t* cell_distri = PDM_part_to_block_distrib_index_get(ptb_equi);

  int dn_cell_equi = cell_distri[i_rank+1] - cell_distri[i_rank];
  int         **pelmt_to_arc_n     = malloc(extrp->n_part_in * sizeof(int         *));
  PDM_g_num_t **pelmt_to_arc_g_num = malloc(extrp->n_part_in * sizeof(PDM_g_num_t *));

  int from_face_edge = 0;
  int from_face_vtx  = 0;
  for(int i_part = 0; i_part < extrp->n_part_in; ++i_part) {
    if(extrp->pface_edge    [i_part] != NULL) {
      from_face_edge = 1;
    }
    if(extrp->pface_vtx    [i_part] != NULL) {
      from_face_vtx = 1;
    }
  }

  int         **pelmt_to_arc_idx = NULL;
  int         **pelmt_to_arc     = NULL;
  PDM_g_num_t **arc_ln_to_gn     = NULL;

  if(extrp->dim == 3) {
    pelmt_to_arc_idx = extrp->pcell_face_idx;
    pelmt_to_arc     = extrp->pcell_face;
    arc_ln_to_gn     = extrp->face_ln_to_gn;
  } else if (extrp->dim == 2) {
    if(from_face_edge == 1) {
      pelmt_to_arc_idx = extrp->pface_edge_idx;
      pelmt_to_arc     = extrp->pface_edge;
      arc_ln_to_gn     = extrp->edge_ln_to_gn;
    } else {
      assert(from_face_vtx == 1);
      pelmt_to_arc_idx = extrp->pface_vtx_idx;
      pelmt_to_arc     = extrp->pface_vtx;
      arc_ln_to_gn     = extrp->vtx_ln_to_gn;
    }
  } else if (extrp->dim == 1) {
    pelmt_to_arc_idx = malloc(extrp->n_part_in * sizeof(int *));
    for(int i_part = 0; i_part < extrp->n_part_in; ++i_part) {
      pelmt_to_arc_idx[i_part] = malloc((extrp->n_edge[i_part]+1) * sizeof(int));
      for(int i = 0; i < extrp->n_edge[i_part]+1; ++i) {
        pelmt_to_arc_idx[i_part][i] = 2*i;
      }
    }
    pelmt_to_arc     = extrp->pedge_vtx;
    arc_ln_to_gn     = extrp->vtx_ln_to_gn;
  } else {
    PDM_error(__FILE__, __LINE__, 0,"PDM_extract_part_compute : cannot not use split_method !=  PDM_SPLIT_DUAL_WITH_HILBERT with dim=0 (use PDM_SPLIT_DUAL_WITH_HILBERT instead)\n");
  }


  /*
   * Convert to gnum
   */
  for(int i_part = 0; i_part < extrp->n_part_in; ++i_part ) {
    pelmt_to_arc_n[i_part] = malloc(extrp->n_extract[i_part] * sizeof(int));

    int n_extract_tot_cell_face = 0;
    for(int idx_entity = 0; idx_entity < extrp->n_extract[i_part]; idx_entity++) {
      int lnum = extrp->extract_lnum[i_part][idx_entity];
      pelmt_to_arc_n[i_part][idx_entity] = pelmt_to_arc_idx[i_part][lnum+1]-pelmt_to_arc_idx[i_part][lnum];
      n_extract_tot_cell_face += pelmt_to_arc_n[i_part][idx_entity];
    }

    pelmt_to_arc_g_num[i_part] = malloc(n_extract_tot_cell_face * sizeof(PDM_g_num_t));
    int idx_write = 0;
    for(int idx_entity = 0; idx_entity < extrp->n_extract[i_part]; idx_entity++) {
      int lnum = extrp->extract_lnum[i_part][idx_entity];
      for(int idx_face = pelmt_to_arc_idx[i_part][lnum]; idx_face < pelmt_to_arc_idx[i_part][lnum+1]; ++idx_face) {
        int i_face = PDM_ABS(pelmt_to_arc[i_part][idx_face])-1;
        pelmt_to_arc_g_num[i_part][idx_write++] = arc_ln_to_gn[i_part][i_face];
      }
    }
  }

  if (extrp->dim == 1) {
    for(int i_part = 0; i_part < extrp->n_part_in; ++i_part) {
      free(pelmt_to_arc_idx[i_part]);
    }
    free(pelmt_to_arc_idx);
  }

  int request_elmt_to_arc = -1;
  int         *delmt_to_arc_n = NULL;
  PDM_g_num_t *delmt_to_arc = NULL;
  PDM_part_to_block_iexch(ptb_equi,
                          PDM_MPI_COMM_KIND_COLLECTIVE,
                          sizeof(PDM_g_num_t),
                          PDM_STRIDE_VAR_INTERLACED,
                          1,
                          pelmt_to_arc_n,
                (void **) pelmt_to_arc_g_num,
                          &delmt_to_arc_n,
                (void **) &delmt_to_arc,
                          &request_elmt_to_arc);

  for(int i_part = 0; i_part < extrp->n_part_in; ++i_part ) {
    free(pelmt_to_arc_n    [i_part]);
    free(pelmt_to_arc_g_num[i_part]);
  }
  free(pelmt_to_arc_n);
  free(pelmt_to_arc_g_num);


  int n_tot_arc = PDM_part_to_block_iexch_wait(ptb_equi, request_elmt_to_arc);

  if(0 == 1) {
    int *delmt_to_arc_idx = PDM_array_new_idx_from_sizes_int(delmt_to_arc_n, dn_cell_equi);
    PDM_log_trace_connectivity_long(delmt_to_arc_idx, delmt_to_arc, dn_cell_equi, "delmt_to_arc ::");
    free(delmt_to_arc_idx);
  }

  /*
   * Transpose connectivity
   */
  PDM_part_to_block_t* ptb_merge_arc = PDM_part_to_block_create(PDM_PART_TO_BLOCK_DISTRIB_ALL_PROC,
                                                                PDM_PART_TO_BLOCK_POST_MERGE,
                                                                1.,
                                                                &delmt_to_arc,
                                                                NULL,
                                                                &n_tot_arc,
                                                                1,
                                                                extrp->comm);

  int *stride_one = PDM_array_const_int(n_tot_arc, 1);
  PDM_g_num_t* send_elmt_ln_to_gn = malloc(n_tot_arc * sizeof(PDM_g_num_t));

  int n_elmt_to_arc_tot = 0;
  for(int i = 0; i < dn_cell_equi; ++i) {
    for(int j = 0; j < delmt_to_arc_n[i]; ++j) {
      send_elmt_ln_to_gn[n_elmt_to_arc_tot++] = cell_distri[i_rank] + i + 1;
    }
  }

  int request_arc_to_elmt = -1;
  int         *darc_to_elmt_n = NULL;
  PDM_g_num_t *darc_to_elmt   = NULL;
  PDM_part_to_block_iexch(ptb_merge_arc,
                          PDM_MPI_COMM_KIND_COLLECTIVE,
                          sizeof(PDM_g_num_t),
                          PDM_STRIDE_VAR_INTERLACED,
                          1,
                          &stride_one,
                (void **) &send_elmt_ln_to_gn,
                          &darc_to_elmt_n,
                (void **) &darc_to_elmt,
                          &request_arc_to_elmt);

  PDM_part_to_block_iexch_wait(ptb_merge_arc, request_arc_to_elmt);
  free(stride_one);
  free(send_elmt_ln_to_gn);

  int dn_arc_equi = PDM_part_to_block_n_elt_block_get(ptb_merge_arc);
  if(0 == 1) {
    int *darc_to_elmt_idx = PDM_array_new_idx_from_sizes_int(darc_to_elmt_n, dn_arc_equi);
    PDM_log_trace_connectivity_long(darc_to_elmt_idx, darc_to_elmt, dn_arc_equi, "darc_to_elmt ::");
    free(darc_to_elmt_idx);
  }

  /*
   * Renvoie des numero delmts associé à chaque arc
   */
  PDM_g_num_t *arc_gnum = PDM_part_to_block_block_gnum_get(ptb_merge_arc);

  PDM_block_to_part_t* btp = PDM_block_to_part_create_from_sparse_block(arc_gnum,
                                                                        dn_arc_equi,
                                                (const PDM_g_num_t **)  &delmt_to_arc,
                                                                        &n_elmt_to_arc_tot,
                                                                        1,
                                                                        extrp->comm);

  int         **tmp_parc_to_elmt_n = NULL;
  PDM_g_num_t **tmp_parc_to_elmt   = NULL;
  PDM_block_to_part_exch(btp,
                         sizeof(PDM_g_num_t),
                         PDM_STRIDE_VAR_INTERLACED,
                         darc_to_elmt_n,
                         darc_to_elmt,
                         &tmp_parc_to_elmt_n,
         (void ***)      &tmp_parc_to_elmt);
  int         *parc_to_elmt_n = tmp_parc_to_elmt_n[0];
  PDM_g_num_t *parc_to_elmt   = tmp_parc_to_elmt  [0];
  free(tmp_parc_to_elmt_n);
  free(tmp_parc_to_elmt  );

  /*
   * Merge results
   */
  int idx_read = 0;
  int n_data_tot = 0;
  int max_connect = 0;
  for(int i = 0; i < dn_cell_equi; ++i) {
    int ln_arc = delmt_to_arc_n[i];
    int n_data = 0;
    for(int i_arc = 0; i_arc < ln_arc; ++i_arc) {
      n_data += parc_to_elmt_n[idx_read++];
    }
    max_connect = PDM_MAX(max_connect, n_data);
    n_data_tot += n_data;
  }

  PDM_g_num_t *dual_graph     = malloc(n_data_tot       * sizeof(PDM_g_num_t));
  PDM_g_num_t *dual_graph_idx = malloc((dn_cell_equi+1) * sizeof(PDM_g_num_t));
  PDM_g_num_t *tmp_dual_graph = malloc(max_connect      * sizeof(PDM_g_num_t));

  idx_read = 0;
  int idx_read_data = 0;
  dual_graph_idx[0] = 0;
  for(int i = 0; i < dn_cell_equi; ++i) {
    PDM_g_num_t cell_g_num = cell_distri[i_rank] + i + 1;

    dual_graph_idx[i+1] = dual_graph_idx[i];

    int ln_arc = delmt_to_arc_n[i];
    int idx_cpy = 0;
    for(int j = 0; j < ln_arc; ++j) {
      for(int k = 0; k < parc_to_elmt_n[idx_read]; ++k) {
        PDM_g_num_t elmt_gnum = parc_to_elmt[idx_read_data++];
        if(elmt_gnum != cell_g_num) {
          tmp_dual_graph[idx_cpy++] = elmt_gnum;
        }

      }
      idx_read++;
    }

    int ln_elmt = PDM_inplace_unique_long(tmp_dual_graph, NULL, 0, idx_cpy-1);
    for(int j = 0; j < ln_elmt; ++j) {
      dual_graph[dual_graph_idx[i+1]++] = tmp_dual_graph[j]-1; // Graphe start at 0 for metis/scoth
    }
  }

  dual_graph = realloc(dual_graph, dual_graph_idx[dn_cell_equi] * sizeof(PDM_g_num_t));

  free(parc_to_elmt_n);
  free(parc_to_elmt  );
  PDM_block_to_part_free(btp);
  PDM_part_to_block_free(ptb_merge_arc);

  free(delmt_to_arc_n);
  free(delmt_to_arc);
  free(darc_to_elmt_n);
  free(darc_to_elmt);
  free(tmp_dual_graph);

  *out_dual_graph_idx = dual_graph_idx;
  *out_dual_graph     = dual_graph;

}



static
void
extract_entity1_entity2_new
(
int            n_part,
int           *n_entity2,
int           *n_extract_entity1,
int          **extract_entity1_lnum,
int          **entity1_entity2_idx,
int          **entity1_entity2,
PDM_g_num_t  **entity2_ln_to_gn,
int          **n_extract_entity2,
int         ***extract_entity2_lnum,
PDM_g_num_t ***extract_parent_entity2_ln_to_gn,
int         ***old_to_new_entity2_no
)
{
  int         **_extract_entity2_lnum            = malloc(n_part * sizeof(int         *));
  PDM_g_num_t **_extract_parent_entity2_ln_to_gn = malloc(n_part * sizeof(PDM_g_num_t *));
  int          *_n_extract_entity2               = malloc(n_part * sizeof(int          ));
  int         **_old_to_new_entity2_no           = malloc(n_part * sizeof(int         *));

  for(int i_part = 0; i_part < n_part; ++i_part) {
    int n_extract = n_extract_entity1[i_part];
    int         *_pentity1_entity2     = entity1_entity2    [i_part];
    int         *_pentity1_entity2_idx = entity1_entity2_idx[i_part];
    PDM_g_num_t *_pentity2_ln_to_gn    = entity2_ln_to_gn   [i_part];
    int          _pn_entity2           = n_entity2          [i_part];

    _extract_entity2_lnum           [i_part] = (int          *) malloc( (_pn_entity2 ) * sizeof(int        ));
    _extract_parent_entity2_ln_to_gn[i_part] = (PDM_g_num_t  *) malloc( (_pn_entity2 ) * sizeof(PDM_g_num_t));
    _old_to_new_entity2_no          [i_part] = (int          *) malloc( (_pn_entity2 ) * sizeof(int        ));

    int         *is_visited = (int *) malloc( _pn_entity2 * sizeof(int));
    for(int i = 0; i < _pn_entity2; ++i) {
      is_visited     [i] = 0;
      _old_to_new_entity2_no[i_part][i] = -1;
    }

    int idx_write = 0;
    _n_extract_entity2[i_part] = 0;
    for(int idx_entity = 0; idx_entity < n_extract; ++idx_entity) {
      int i_entity = extract_entity1_lnum[i_part][idx_entity];

      for(int idx_entity2 = _pentity1_entity2_idx[i_entity]; idx_entity2 < _pentity1_entity2_idx[i_entity+1]; ++idx_entity2) {
        int i_entity2 = PDM_ABS(_pentity1_entity2[idx_entity2])-1;
        if(is_visited[i_entity2] == 0) {
          int idx = _n_extract_entity2[i_part]++;
          _extract_entity2_lnum           [i_part][idx] = i_entity2;
          _extract_parent_entity2_ln_to_gn[i_part][idx] = _pentity2_ln_to_gn[i_entity2];
          is_visited                    [i_entity2] = 1;
          _old_to_new_entity2_no[i_part][i_entity2] = idx_write++;
        }
      }
    }
    _extract_entity2_lnum           [i_part] = realloc(_extract_entity2_lnum           [i_part], _n_extract_entity2[i_part] * sizeof(int        ));
    _extract_parent_entity2_ln_to_gn[i_part] = realloc(_extract_parent_entity2_ln_to_gn[i_part], _n_extract_entity2[i_part] * sizeof(PDM_g_num_t));

    free(is_visited);
  }

  *n_extract_entity2               = _n_extract_entity2;
  *extract_entity2_lnum            = _extract_entity2_lnum;
  *extract_parent_entity2_ln_to_gn = _extract_parent_entity2_ln_to_gn;
  *old_to_new_entity2_no           = _old_to_new_entity2_no;
}

static
void
extract_and_local_renum_entity1_entity2
(
PDM_MPI_Comm           comm,
int                    compute_child_gnum,
int                    n_part,
int                   *n_entity1,
int                   *n_entity2,
int                   *n_extract_entity1,
int                  **extract_entity1_lnum,
int                  **entity1_entity2_idx,
int                  **entity1_entity2,
PDM_g_num_t          **entity2_ln_to_gn,
int                  **n_extract_entity2,
int                 ***selected_entity1_entity2_idx,
int                 ***selected_entity1_entity2,
PDM_g_num_t         ***selected_entity2_ln_to_gn,
PDM_g_num_t         ***selected_parent_entity2_ln_to_gn,
int                 ***extract_entity2_lnum
)
{
  PDM_g_num_t **_extract_parent_entity2_ln_to_gn = NULL;
  int         **old_to_new_entity2_no            = NULL;

  /*
   *  Compute extract_entity2_lnum and extract_entity2_g_num
   */
  extract_entity1_entity2_new(n_part,
                              n_entity2,
                              n_extract_entity1,
                              extract_entity1_lnum,
                              entity1_entity2_idx,
                              entity1_entity2,
                              entity2_ln_to_gn,
                              n_extract_entity2,
                              extract_entity2_lnum,
                              &_extract_parent_entity2_ln_to_gn,
                              &old_to_new_entity2_no);


  int *_n_extract_entity2 = *n_extract_entity2;

  int         **_selected_entity1_entity2_idx = malloc(n_part * sizeof(int         *));
  int         **_selected_entity1_entity2     = malloc(n_part * sizeof(int         *));
  for(int i_part = 0; i_part < n_part; ++i_part) {

    int  _pn_entity1           = n_entity1          [i_part];
    int *_pentity1_entity2     = entity1_entity2    [i_part];
    int *_pentity1_entity2_idx = entity1_entity2_idx[i_part];

    _selected_entity1_entity2    [i_part] = (int  *) malloc( _pentity1_entity2_idx[_pn_entity1] * sizeof(int));
    _selected_entity1_entity2_idx[i_part] = (int  *) malloc( (n_extract_entity1[i_part]+1)      * sizeof(int));

    //
    int idx_write = 0;
    _selected_entity1_entity2_idx[i_part][0] = 0;
    for(int idx_entity = 0; idx_entity < n_extract_entity1[i_part]; ++idx_entity) {
      int i_entity = extract_entity1_lnum[i_part][idx_entity];

      int n_tmp = _pentity1_entity2_idx[i_entity+1] - _pentity1_entity2_idx[i_entity];
      _selected_entity1_entity2_idx[i_part][idx_entity+1] = _selected_entity1_entity2_idx[i_part][idx_entity] + n_tmp;

      for(int idx_entity2 = _pentity1_entity2_idx[i_entity]; idx_entity2 < _pentity1_entity2_idx[i_entity+1]; ++idx_entity2) {
        int i_entity2         = PDM_ABS (_pentity1_entity2[idx_entity2])-1;
        int sgn               = PDM_SIGN(_pentity1_entity2[idx_entity2]);
        int i_extract_entity2 = old_to_new_entity2_no[i_part][i_entity2];
        _selected_entity1_entity2[i_part][idx_write++] = sgn * (i_extract_entity2 + 1);
      }
    }

    assert(idx_write == _selected_entity1_entity2_idx[i_part][n_extract_entity1[i_part]]);

    _selected_entity1_entity2[i_part] = realloc(_selected_entity1_entity2[i_part], idx_write * sizeof(int));
    free(old_to_new_entity2_no[i_part]);
  }
  free(old_to_new_entity2_no);


  PDM_g_num_t **_child_entity2_ln_to_gn = malloc(n_part * sizeof(PDM_g_num_t *));

  for(int i_part = 0; i_part < n_part; ++i_part) {
    _child_entity2_ln_to_gn[i_part] = NULL;
  }
  if(compute_child_gnum  ==  1) {
    PDM_gen_gnum_t* gnum_extract_entity2 = PDM_gnum_create(3,
                                                           n_part,
                                                           PDM_FALSE,
                                                           1.e-6,
                                                           comm,
                                                           PDM_OWNERSHIP_USER);

    for(int i_part = 0; i_part < n_part; ++i_part) {
      PDM_gnum_set_from_parents(gnum_extract_entity2, i_part, _n_extract_entity2[i_part], _extract_parent_entity2_ln_to_gn[i_part]);
    }

    PDM_gnum_compute(gnum_extract_entity2);

    for(int i_part = 0; i_part < n_part; ++i_part) {
      _child_entity2_ln_to_gn[i_part] = PDM_gnum_get(gnum_extract_entity2, i_part);
      if(0 == 1) {
        PDM_log_trace_array_long(_child_entity2_ln_to_gn[i_part], _n_extract_entity2[i_part], "_child_entity2_ln_to_gn :: ");
      }
    }
    PDM_gnum_free(gnum_extract_entity2);
  }

  *selected_entity1_entity2_idx     = _selected_entity1_entity2_idx;
  *selected_entity1_entity2         = _selected_entity1_entity2;
  *selected_entity2_ln_to_gn        = _child_entity2_ln_to_gn;
  *selected_parent_entity2_ln_to_gn = _extract_parent_entity2_ln_to_gn;

}

static
void
_extract_part_nodal
(
  PDM_extract_part_t        *extrp
)
{
  int          *pn_entity    = 0;
  PDM_g_num_t **entity_g_num = NULL;
  if(extrp->dim == 3) {
    pn_entity    = extrp->n_cell;
    entity_g_num = extrp->cell_ln_to_gn;
  } else {
    pn_entity    = extrp->n_face;
    entity_g_num = extrp->face_ln_to_gn;
  }

  /*
   *  Create array selected in gnum
   */
  PDM_gen_gnum_t* gnum_extract = PDM_gnum_create(3, extrp->n_part_in, PDM_FALSE,
                                                 1.e-6,
                                                 extrp->comm,
                                                 PDM_OWNERSHIP_USER);
  PDM_g_num_t** entity_extract_g_num = (PDM_g_num_t **) malloc( extrp->n_part_in * sizeof(PDM_g_num_t *));

  for(int i_part = 0; i_part < extrp->n_part_in; ++i_part) {
    entity_extract_g_num[i_part] = (PDM_g_num_t *) malloc( extrp->n_extract[i_part] * sizeof(PDM_g_num_t));
    for(int i_entity = 0; i_entity < extrp->n_extract[i_part]; ++i_entity) {
      entity_extract_g_num[i_part][i_entity] = entity_g_num[i_part][extrp->extract_lnum[i_part][i_entity]];
    }

    PDM_gnum_set_from_parents(gnum_extract, i_part, extrp->n_extract[i_part], entity_extract_g_num[i_part]);
    if(0 == 1) {
      PDM_log_trace_array_long(entity_extract_g_num[i_part], extrp->n_extract[i_part], "entity_extract_g_num ::" );
    }
  }

  /*
   * Global numering computation
   */
  PDM_g_num_t **child_selected_g_num = (PDM_g_num_t **) malloc( extrp->n_part_in * sizeof(PDM_g_num_t *));
  PDM_gnum_compute(gnum_extract);

  for (int i_part = 0; i_part < extrp->n_part_in; i_part++){
    child_selected_g_num[i_part] = PDM_gnum_get(gnum_extract, i_part);
    // PDM_log_trace_array_long(child_selected_g_num[i_part], extrp->n_extract[i_part], "child_selected_g_num : ");
  }
  PDM_gnum_free(gnum_extract);

  if(extrp->dim == 3) {
    extrp->pextract_n_entity       [PDM_MESH_ENTITY_CELL] = malloc(extrp->n_part_out * sizeof(int));
    for(int i_part = 0; i_part < extrp->n_part_in; ++i_part) {
      extrp->pextract_n_entity[PDM_MESH_ENTITY_CELL][i_part] = extrp->n_extract[i_part];
    }
    extrp->pextract_entity_ln_to_gn       [PDM_MESH_ENTITY_CELL] = child_selected_g_num;
    extrp->pextract_entity_parent_ln_to_gn[PDM_MESH_ENTITY_CELL] = entity_extract_g_num;
  } else {
    extrp->pextract_n_entity       [PDM_MESH_ENTITY_FACE] = malloc(extrp->n_part_out * sizeof(int));
    for(int i_part = 0; i_part < extrp->n_part_in; ++i_part) {
      extrp->pextract_n_entity[PDM_MESH_ENTITY_FACE][i_part] = extrp->n_extract[i_part];
    }
    extrp->pextract_entity_ln_to_gn       [PDM_MESH_ENTITY_FACE] = child_selected_g_num;
    extrp->pextract_entity_parent_ln_to_gn[PDM_MESH_ENTITY_FACE] = entity_extract_g_num;
  }

  int n_section = PDM_part_mesh_nodal_elmts_n_section_get(extrp->pmne);

  int *sections_id = PDM_part_mesh_nodal_elmts_sections_id_get(extrp->pmne);

  assert(extrp->pmne->n_part == extrp->n_part_in);

  /*
   * For each section we extract the selected part
   *
   */
  assert(extrp->n_part_in == extrp->n_part_out);
  extrp->pextract_entity_parent_ln_to_gn[PDM_MESH_ENTITY_VTX] = malloc( extrp->n_part_in * sizeof(PDM_g_num_t *));
  extrp->pextract_entity_ln_to_gn       [PDM_MESH_ENTITY_VTX] = malloc( extrp->n_part_in * sizeof(PDM_g_num_t *));
  extrp->pextract_entity_parent_lnum    [PDM_MESH_ENTITY_VTX] = malloc( extrp->n_part_in * sizeof(int         *));

  int           *n_extract_vtx            = malloc( extrp->n_part_in * sizeof(int          ));
  int          **is_selected              = malloc( extrp->n_part_in * sizeof(int         *));
  int          **is_selected_vtx          = malloc( extrp->n_part_in * sizeof(int         *));
  int          **old_to_new_vtx           = malloc( extrp->n_part_in * sizeof(int         *));
  int          **extract_vtx_lnum         = malloc( extrp->n_part_in * sizeof(int         *));
  PDM_g_num_t  **extract_parent_vtx_g_num = malloc( extrp->n_part_in * sizeof(PDM_g_num_t *));

  int          **n_selected_section       = malloc( extrp->n_part_in * sizeof(int         *));
  int         ***idx_selected_section     = malloc( extrp->n_part_in * sizeof(int        **));
  int         ***extract_parent_num       = malloc( extrp->n_part_in * sizeof(int        **));

  for(int i_part = 0; i_part < extrp->pmne->n_part; ++i_part) {
    extrp->pextract_entity_parent_ln_to_gn[PDM_MESH_ENTITY_VTX][i_part] = NULL;
    extrp->pextract_entity_ln_to_gn       [PDM_MESH_ENTITY_VTX][i_part] = NULL;
    extrp->pextract_entity_parent_lnum    [PDM_MESH_ENTITY_VTX][i_part] = NULL;

    is_selected             [i_part] = malloc(    pn_entity[i_part] * sizeof(int        ));
    is_selected_vtx         [i_part] = malloc( extrp->n_vtx[i_part] * sizeof(int        ));
    old_to_new_vtx          [i_part] = malloc( extrp->n_vtx[i_part] * sizeof(int        ));
    extract_vtx_lnum        [i_part] = malloc( extrp->n_vtx[i_part] * sizeof(int        ));
    extract_parent_vtx_g_num[i_part] = malloc( extrp->n_vtx[i_part] * sizeof(PDM_g_num_t));

    /*
     * En polyhédrique il faut aussi les faces ou edges a extraire
     */
    for(int i = 0; i < pn_entity[i_part]; ++i) {
      is_selected[i_part][i] = -1;
    }
    for(int i = 0; i < extrp->n_vtx[i_part]; ++i) {
      is_selected_vtx[i_part][i] = 0;
      old_to_new_vtx [i_part][i] = -1;
    }

    for(int i = 0; i < extrp->n_extract[i_part]; ++i) {
      int s_num = extrp->extract_lnum[i_part][i];
      is_selected[i_part][s_num] = i;
    }

    n_extract_vtx[i_part] = 0;

    n_selected_section  [i_part] = malloc( n_section * sizeof(int  ));
    idx_selected_section[i_part] = malloc( n_section * sizeof(int *));
    extract_parent_num  [i_part] = malloc( n_section * sizeof(int *));

    /* First pass to hook all vtx and create gnum */
    for(int i_section = 0; i_section < n_section; ++i_section) {


      int n_elt = PDM_part_mesh_nodal_elmts_section_n_elt_get(extrp->pmne, sections_id[i_section], i_part);
      PDM_Mesh_nodal_elt_t t_elt = PDM_part_mesh_nodal_elmts_section_type_get(extrp->pmne, sections_id[i_section]);

      n_selected_section  [i_part][i_section] = 0;
      idx_selected_section[i_part][i_section] = malloc( n_elt * sizeof(int));
      extract_parent_num  [i_part][i_section] = malloc( n_elt * sizeof(int));

      int         *elt_vtx          = NULL;
      int         *_parent_num      = NULL;
      PDM_g_num_t *elt_ln_to_gn     = NULL;
      PDM_g_num_t *parent_elt_g_num = NULL;
      int          order            = 0;
      const char  *ho_ordering      = NULL;

      PDM_part_mesh_nodal_elmts_section_std_ho_get(extrp->pmne,
                                                 sections_id[i_section],
                                                 i_part,
                                                 &elt_vtx,
                                                 &elt_ln_to_gn,
                                                 &_parent_num,
                                                 &parent_elt_g_num,
                                                 &order,
                                                 &ho_ordering,
                                                 PDM_OWNERSHIP_KEEP);

      int n_vtx_per_elmt = PDM_Mesh_nodal_n_vtx_elt_get(t_elt, order);

      int *ijk_to_user = NULL;
      if (order > 1 && ho_ordering != NULL) {
        ijk_to_user = PDM_ho_ordering_ijk_to_user_get(ho_ordering,
                                                      t_elt,
                                                      order);
      }

      /* Selection */
      for(int i_elt = 0; i_elt < n_elt; ++i_elt) {
        int parent_elt = i_elt;
        if (_parent_num != NULL) {
          parent_elt = _parent_num[i_elt];
        }
        if(is_selected[i_part][parent_elt] != -1) {

          int idx_write = n_selected_section[i_part][i_section]++;
          idx_selected_section[i_part][i_section][idx_write] = i_elt;
          extract_parent_num  [i_part][i_section][idx_write] = is_selected[i_part][parent_elt];

          int beg = i_elt * n_vtx_per_elmt;
          for(int idx_vtx = 0; idx_vtx < n_vtx_per_elmt; ++idx_vtx) {
            int _idx_vtx = idx_vtx;
            if (ijk_to_user != NULL) {
              _idx_vtx = ijk_to_user[idx_vtx];
            }
            int i_vtx = elt_vtx[beg+_idx_vtx] - 1;
            if(is_selected_vtx[i_part][i_vtx] == 0) {
              is_selected_vtx         [i_part][i_vtx                ] = 1;
              old_to_new_vtx          [i_part][i_vtx                ] = n_extract_vtx[i_part];
              extract_parent_vtx_g_num[i_part][n_extract_vtx[i_part]] = extrp->vtx_ln_to_gn[i_part][i_vtx];
              extract_vtx_lnum        [i_part][n_extract_vtx[i_part]] = i_vtx;
              n_extract_vtx[i_part]++;
            }
          }
        }
      }


      idx_selected_section[i_part][i_section] = realloc(idx_selected_section[i_part][i_section], n_selected_section  [i_part][i_section] * sizeof(int));
      extract_parent_num  [i_part][i_section] = realloc(extract_parent_num  [i_part][i_section], n_selected_section  [i_part][i_section] * sizeof(int));

    } /* End section */

    extract_vtx_lnum        [i_part] = realloc(extract_vtx_lnum        [i_part], n_extract_vtx[i_part] * sizeof(int        ));
    extract_parent_vtx_g_num[i_part] = realloc(extract_parent_vtx_g_num[i_part], n_extract_vtx[i_part] * sizeof(PDM_g_num_t));

    extrp->pextract_entity_parent_ln_to_gn[PDM_MESH_ENTITY_VTX][i_part] = extract_parent_vtx_g_num[i_part];
    extrp->pextract_entity_parent_lnum    [PDM_MESH_ENTITY_VTX][i_part] = extract_vtx_lnum[i_part];

  }


  /*
   *  Create absolute numbering of vtx
   */
  PDM_gen_gnum_t* gnum_extract_vtx = PDM_gnum_create(3, extrp->n_part_in, PDM_FALSE,
                                                     1.e-6,
                                                     extrp->comm,
                                                     PDM_OWNERSHIP_USER);
  for(int i_part = 0; i_part < extrp->n_part_in; ++i_part) {
    PDM_gnum_set_from_parents(gnum_extract_vtx, i_part, n_extract_vtx[i_part], extract_parent_vtx_g_num[i_part]);
  }

  PDM_gnum_compute(gnum_extract_vtx);

  for(int i_part = 0; i_part < extrp->n_part_in; ++i_part) {
    extrp->pextract_entity_ln_to_gn[PDM_MESH_ENTITY_VTX][i_part] = PDM_gnum_get(gnum_extract_vtx, i_part);
  }
  PDM_gnum_free(gnum_extract_vtx);


  /*
   * Second pass to create the new part_mesh_nodal
   */
  PDM_part_mesh_nodal_elmts_t* extract_pmne = PDM_part_mesh_nodal_elmts_create(extrp->pmne->mesh_dimension,
                                                                               extrp->pmne->n_part,
                                                                               extrp->pmne->comm);



  for(int i_part = 0; i_part < extrp->pmne->n_part; ++i_part) {

    for(int i_section = 0; i_section < n_section; ++i_section) {

      PDM_Mesh_nodal_elt_t t_elt = PDM_part_mesh_nodal_elmts_section_type_get(extrp->pmne, sections_id[i_section]);

      int         *elt_vtx          = NULL;
      int         *_parent_num      = NULL;
      PDM_g_num_t *elt_ln_to_gn     = NULL;
      PDM_g_num_t *parent_elt_g_num = NULL;
      int          order            = 0;
      const char  *ho_ordering      = NULL;

      PDM_part_mesh_nodal_elmts_section_std_ho_get(extrp->pmne,
                                                 sections_id[i_section],
                                                 i_part,
                                                 &elt_vtx,
                                                 &elt_ln_to_gn,
                                                 &_parent_num,
                                                 &parent_elt_g_num,
                                                 &order,
                                                 &ho_ordering,
                                                 PDM_OWNERSHIP_KEEP);

      int n_vtx_per_elmt = PDM_Mesh_nodal_n_vtx_elt_get(t_elt, order);

      int *ijk_to_user = NULL;
      if (order > 1 && ho_ordering != NULL) {
        ijk_to_user = PDM_ho_ordering_ijk_to_user_get(ho_ordering,
                                                      t_elt,
                                                      order);
      }

      /* Allocate */
      int         *extract_elt_vtx      = malloc( n_selected_section  [i_part][i_section] * n_vtx_per_elmt * sizeof(int        ));
      PDM_g_num_t *extract_elt_ln_to_gn = malloc( n_selected_section  [i_part][i_section]                  * sizeof(PDM_g_num_t));

      // PDM_log_trace_array_int(idx_selected_section[i_part][i_section], n_selected_section  [i_part][i_section], "idx_selected_section :");
      // PDM_log_trace_array_int(old_to_new_vtx[i_part], extrp->n_vtx[i_part], "old_to_new_vtx :");

      int idx_write = 0;
      for(int i = 0; i < n_selected_section  [i_part][i_section]; ++i) {
        int ielt = idx_selected_section[i_part][i_section][i];
        int beg  = ielt * n_vtx_per_elmt;

        for(int k = 0; k < n_vtx_per_elmt; ++k) {
          int _k = k;
          if (ijk_to_user != NULL) {
            _k = ijk_to_user[k];
          }
          int old_vtx = elt_vtx[beg+_k]-1;
          extract_elt_vtx[idx_write++] = old_to_new_vtx[i_part][old_vtx]+1;
        }

        int idx_parent = extract_parent_num[i_part][i_section][i];
        extract_elt_ln_to_gn[i] = child_selected_g_num[i_part][idx_parent];
      }

      /* Fill up structure */
      // int extract_section_id = PDM_part_mesh_nodal_elmts_add(extract_pmne, t_elt);
      // PDM_part_mesh_nodal_elmts_std_set(extract_pmne,
      //                                   extract_section_id,
      //                                   i_part,
      //                                   n_selected_section[i_part][i_section],
      //                                   extract_elt_vtx,
      //                                   extract_elt_ln_to_gn,
      //                                   extract_parent_num[i_part][i_section],
      //                                   NULL,
      //                                   PDM_OWNERSHIP_KEEP);
      int extract_section_id = PDM_part_mesh_nodal_elmts_add(extract_pmne, t_elt);
      PDM_part_mesh_nodal_elmts_std_ho_set(extract_pmne,
                                           extract_section_id,
                                           i_part,
                                           n_selected_section[i_part][i_section],
                                           extract_elt_vtx,
                                           extract_elt_ln_to_gn,
                                           extract_parent_num[i_part][i_section],
                                           NULL,
                                           order,
                                           NULL,//ho_ordering,
                                           PDM_OWNERSHIP_KEEP);

      free(idx_selected_section[i_part][i_section]);

    }

    free(idx_selected_section[i_part]);
    free(n_selected_section  [i_part]);
    free(extract_parent_num  [i_part]);

  }


  free(idx_selected_section);
  free(n_selected_section  );
  free(extract_parent_num  );



  for(int i_part = 0; i_part < extrp->n_part_in; ++i_part) {
    free(is_selected    [i_part]);
    free(is_selected_vtx[i_part]);
    free(old_to_new_vtx [i_part]);
  }

  free(is_selected    );
  free(is_selected_vtx);
  free(old_to_new_vtx );

  /*
   * Extract coordinates
   */
  extrp->pextract_vtx_coord = (double **) malloc( extrp->n_part_in * sizeof(double *));
  for(int i_part = 0; i_part < extrp->n_part_in; ++i_part) {
    extrp->pextract_vtx_coord[i_part] = (double *) malloc( 3 * n_extract_vtx[i_part] * sizeof(double));

    for(int idx_vtx = 0; idx_vtx < n_extract_vtx[i_part]; ++idx_vtx) {
      int i_vtx = extract_vtx_lnum[i_part][idx_vtx];
      extrp->pextract_vtx_coord[i_part][3*idx_vtx  ] = extrp->pvtx_coord[i_part][3*i_vtx  ];
      extrp->pextract_vtx_coord[i_part][3*idx_vtx+1] = extrp->pvtx_coord[i_part][3*i_vtx+1];
      extrp->pextract_vtx_coord[i_part][3*idx_vtx+2] = extrp->pvtx_coord[i_part][3*i_vtx+2];
    }
  }

  extrp->pextract_n_entity[PDM_MESH_ENTITY_VTX] = n_extract_vtx;

  extrp->extract_pmne = extract_pmne;


  if(0 == 1){
    int i_rank;
    PDM_MPI_Comm_rank(extrp->comm, &i_rank);

    // int *extract_sections_id = PDM_part_mesh_nodal_elmts_sections_id_get(extrp->extract_pmne);
    for(int i_part = 0; i_part < extrp->n_part_in; ++i_part) {

      char filename[999];
      sprintf(filename, "out_extract_%i_%i.vtk", i_part, i_rank);

      int id_section = 0;
      PDM_Mesh_nodal_elt_t t_elt = PDM_part_mesh_nodal_elmts_section_type_get(extract_pmne, id_section);
      int         *elmt_vtx                 = NULL;
      int         *parent_num               = NULL;
      PDM_g_num_t *numabs                   = NULL;
      PDM_g_num_t *parent_entitity_ln_to_gn = NULL;
      PDM_part_mesh_nodal_elmts_section_std_get(extract_pmne,
                                                id_section,
                                                i_part,
                                                &elmt_vtx,
                                                &numabs,
                                                &parent_num,
                                                &parent_entitity_ln_to_gn,
                                                PDM_OWNERSHIP_KEEP);

      PDM_vtk_write_std_elements(filename,
                                 n_extract_vtx[i_part],
                                 extrp->pextract_vtx_coord[i_part],
                                 extrp->pextract_entity_ln_to_gn[PDM_MESH_ENTITY_VTX][i_part],
                                 t_elt,
                                 extrp->n_extract[i_part],
                                 elmt_vtx,
                                 child_selected_g_num[i_part],
                                 0,
                                 NULL,
                                 NULL);
    }
  }

  free(extract_vtx_lnum        );
  free(extract_parent_vtx_g_num);

}

static
void
_extract_part_and_reequilibrate_nodal_from_target
(
  PDM_extract_part_t        *extrp
)
{
  int          *pn_entity       = NULL;
  PDM_mesh_entities_t entity_type = PDM_MESH_ENTITY_MAX;
  PDM_g_num_t **entity_g_num = NULL;
  if(extrp->dim == 3) {
    pn_entity    = extrp->n_cell;
    entity_type = PDM_MESH_ENTITY_CELL;
    entity_g_num = extrp->cell_ln_to_gn;
  } else if (extrp->dim == 2){
    pn_entity    = extrp->n_face;
    entity_type = PDM_MESH_ENTITY_FACE;
    entity_g_num = extrp->face_ln_to_gn;
  } else if (extrp->dim == 1){
    pn_entity    = extrp->n_edge;
    entity_type = PDM_MESH_ENTITY_EDGE;
    entity_g_num = extrp->edge_ln_to_gn;
  } else {
    PDM_error(__FILE__, __LINE__, 0,"_extract_part_and_reequilibrate_nodal_from_target : wrong entity \n");
  }

  int i_rank;
  PDM_MPI_Comm_rank(extrp->comm, &i_rank);

  extrp->pextract_n_entity[entity_type] = (int *) malloc(extrp->n_part_out * sizeof(int));
  for(int i_part = 0; i_part < extrp->n_part_out; ++i_part) {
    extrp->pextract_n_entity[entity_type][i_part] = extrp->n_target[i_part];
  }

  // Target : reference des cellules
  assert(extrp->pmne != NULL);

  int **part2_cell_to_part1_cell_idx = (int **) malloc( extrp->n_part_out * sizeof(int * ));
  for(int i_part = 0; i_part < extrp->n_part_out; ++i_part) {
    part2_cell_to_part1_cell_idx[i_part] = (int * ) malloc( (extrp->n_target[i_part]+1) * sizeof(int));
    part2_cell_to_part1_cell_idx[i_part][0] = 0;
    for(int i = 0; i < extrp->n_target[i_part]; ++i) {
      part2_cell_to_part1_cell_idx[i_part][i+1] = part2_cell_to_part1_cell_idx[i_part][i] + 3;
    }
  }

  PDM_part_to_part_t* ptp = PDM_part_to_part_create_from_num2_triplet((const PDM_g_num_t **) extrp->target_gnum,
                                                                      (const int         * ) extrp->n_target,
                                                                                             extrp->n_part_out,
                                                                                             pn_entity,
                                                                                             extrp->n_part_in,
                                                                      (const int         **) part2_cell_to_part1_cell_idx,
                                                                                             NULL,
                                                                      (const int         **) extrp->target_location,
                                                                      extrp->comm);
  extrp->ptp_entity[entity_type] = ptp;

  for(int i_part = 0; i_part < extrp->n_part_out; ++i_part) {
    free(part2_cell_to_part1_cell_idx[i_part]);
  }
  free(part2_cell_to_part1_cell_idx);

  /*
   * Protocol are created then we can extract information in part1 to reverse send it to part2
   */
  int          *n_ref_entity1     = NULL;
  int         **ref_l_num_entity1 = NULL;
  PDM_part_to_part_ref_lnum2_get(ptp, &n_ref_entity1, &ref_l_num_entity1);

  int         **gnum1_come_from_idx = NULL;
  PDM_g_num_t **gnum1_come_from     = NULL;
  PDM_part_to_part_gnum1_come_from_get(ptp, &gnum1_come_from_idx, &gnum1_come_from);

  int n_section    = PDM_part_mesh_nodal_elmts_n_section_get(extrp->pmne);
  int *sections_id = PDM_part_mesh_nodal_elmts_sections_id_get(extrp->pmne);

  assert(extrp->pmne->n_part == extrp->n_part_in);

  int                   *n_extract_vtx     = malloc(extrp->n_part_in * sizeof(int                   ));
  int                  **is_selected       = malloc(extrp->n_part_in * sizeof(int                  *));

  int                  **elmt_vtx_n        = malloc(extrp->n_part_in * sizeof(int                  *));
  PDM_Mesh_nodal_elt_t **elmt_type         = malloc(extrp->n_part_in * sizeof(PDM_Mesh_nodal_elt_t *));
  PDM_g_num_t          **elmt_vtx          = malloc(extrp->n_part_in * sizeof(PDM_g_num_t          *));
  int                  **vtx_init_location = malloc(extrp->n_part_in * sizeof(int                  *));
  int                  **elmt_section_id   = malloc(extrp->n_part_in * sizeof(int                  *));

  int                  **elmt_face_n       = malloc(extrp->n_part_in * sizeof(int                  *));
  int                  **elmt_face_vtx_n   = malloc(extrp->n_part_in * sizeof(int                  *));
  PDM_g_num_t          **elmt_face         = malloc(extrp->n_part_in * sizeof(PDM_g_num_t          *));

  for(int i_part = 0; i_part < extrp->n_part_in; ++i_part) {

    if(0 == 1) {
      PDM_log_trace_array_int(ref_l_num_entity1[i_part], n_ref_entity1[i_part],"ref_l_num_entity1 :");
      PDM_log_trace_connectivity_long(gnum1_come_from_idx[i_part],
                                      gnum1_come_from    [i_part],
                                      n_ref_entity1  [i_part], "gnum1_come_from ::");
    }

    is_selected[i_part] = malloc(pn_entity[i_part] * sizeof(int));

    for(int i = 0; i < pn_entity[i_part]; ++i) {
      is_selected[i_part][i] = -1;
    }

    // Preparation des buffers d'envoi
    for(int j = 0; j < n_ref_entity1[i_part]; ++j) {
      int i_entity1 = ref_l_num_entity1[i_part][j]-1;
      is_selected[i_part][i_entity1] = j;
    }

    if(0 == 1) {
      PDM_log_trace_array_int(is_selected[i_part], pn_entity[i_part], "is_selected ::");
    }

    /* Compute buffer size */
    int n_elmt_to_send      = 0;
    int n_elmt_vtx_to_send  = 0;
    int n_elmt_face_to_send = 0;
    int parent_elt = -1;
    for(int i_section = 0; i_section < n_section; ++i_section) {

      int n_elt = PDM_part_mesh_nodal_elmts_section_n_elt_get(extrp->pmne, sections_id[i_section], i_part);
      PDM_Mesh_nodal_elt_t t_elt = PDM_part_mesh_nodal_elmts_section_type_get(extrp->pmne, sections_id[i_section]);

      int *parent_num = PDM_part_mesh_nodal_elmts_parent_num_get(extrp->pmne,
                                                                 sections_id[i_section],
                                                                 i_part,
                                                                 PDM_OWNERSHIP_KEEP);

      if (t_elt == PDM_MESH_NODAL_POLY_2D) {
        int *face_vtx_idx;
        int *face_vtx;
        PDM_part_mesh_nodal_elmts_section_poly2d_get(extrp->pmne,
                                                   sections_id[i_section],
                                                   i_part,
                                                   &face_vtx_idx,
                                                   &face_vtx,
                                                   PDM_OWNERSHIP_KEEP);

        // int* parent_num = NULL; // Il faut adpater tout le part_mesh_nodal_elmts

        /* Selection */
        for(int i_elt = 0; i_elt < n_elt; ++i_elt) {
          // int parent_elt = i_elt;
          if (parent_num != NULL) {
            parent_elt = parent_num[i_elt];
          }
          else {
            parent_elt++;
          }
          if(is_selected[i_part][parent_elt] != -1) {
            n_elmt_to_send     += 1;
            n_elmt_vtx_to_send += face_vtx_idx[i_elt+1] - face_vtx_idx[i_elt];
          }
        }

      }
      else if (t_elt == PDM_MESH_NODAL_POLY_3D) {

        /* Polyhedral section */
        // if (1) {
        //   // Check cell-vtx connectivity
        //   int *cell_vtx_idx;
        //   int *cell_vtx;
        //   PDM_part_mesh_nodal_elmts_section_poly3d_cell_vtx_connect_get(extrp->pmne,
        //                                                               sections_id[i_section],
        //                                                               i_part,
        //                                                               &cell_vtx_idx,
        //                                                               &cell_vtx);

        //   log_trace("---Avant---\n");
        //   for (int i = 0; i < n_elt; i++) {
        //     int parent_elt = i;
        //     if (parent_num != NULL) {
        //       parent_elt = parent_num[i];
        //     }
        //     if(is_selected[i_part][parent_elt] != -1) {
        //       // log_trace("cell "PDM_FMT_G_NUM"\n", ?);
        //       log_trace("cell %d: ", i);
        //       for (int j = cell_vtx_idx[i]; j < cell_vtx_idx[i+1]; j++) {
        //         int vtx_id = cell_vtx[j] - 1;
        //         log_trace(PDM_FMT_G_NUM" ", extrp->vtx_ln_to_gn[i_part][vtx_id]);
        //       }
        //       log_trace("\n");
        //     }
        //   }
        // }

        int *cell_face_idx;
        int *cell_face;
        int  n_face;
        int *face_vtx_idx;
        int *face_vtx;
        PDM_g_num_t *face_ln_to_gn    = NULL;
        int         *_parent_num      = NULL; // Il faut adpater tout le part_mesh_nodal_elmts
        PDM_g_num_t *elt_ln_to_gn     = NULL;
        PDM_g_num_t *parent_elt_g_num = NULL;
        PDM_part_mesh_nodal_elmts_section_poly3d_get(extrp->pmne,
                                                   sections_id[i_section],
                                                   i_part,
                                                   &n_face,
                                                   &face_ln_to_gn,
                                                   &face_vtx_idx,
                                                   &face_vtx,
                                                   &elt_ln_to_gn,
                                                   &cell_face_idx,
                                                   &cell_face,
                                                   &_parent_num,
                                                   &parent_elt_g_num,
                                                   PDM_OWNERSHIP_KEEP);

        /* Selection */
        for(int i_elt = 0; i_elt < n_elt; ++i_elt) {
          // int parent_elt = i_elt;
          if (parent_num != NULL) {
            parent_elt = parent_num[i_elt];
          }
          else {
            parent_elt++;
          }
          if(is_selected[i_part][parent_elt] != -1) {
            n_elmt_to_send      += 1;
            n_elmt_face_to_send += cell_face_idx[i_elt+1] - cell_face_idx[i_elt];
            for(int idx_face = cell_face_idx[i_elt]; idx_face < cell_face_idx[i_elt+1]; ++idx_face) {
              int i_face = PDM_ABS(cell_face[idx_face])-1;
              n_elmt_vtx_to_send += face_vtx_idx[i_face+1] - face_vtx_idx[i_face];
            }
          }
        }

      }
      else {
        int         *elt_vtx          = NULL;
        int         *_parent_num      = NULL;
        PDM_g_num_t *elt_ln_to_gn     = NULL;
        PDM_g_num_t *parent_elt_g_num = NULL;
        int          order;
        const char  *ho_ordering;

        PDM_part_mesh_nodal_elmts_section_std_ho_get(extrp->pmne,
                                                   sections_id[i_section],
                                                   i_part,
                                                   &elt_vtx,
                                                   &elt_ln_to_gn,
                                                   &_parent_num,
                                                   &parent_elt_g_num,
                                                   &order,
                                                   &ho_ordering,
                                                   PDM_OWNERSHIP_KEEP);

        int n_vtx_per_elmt = PDM_Mesh_nodal_n_vtx_elt_get(t_elt , order);
        if (0) {
          PDM_log_trace_array_int(parent_num, n_elt, "parent_num : ");
        }

        /* Selection */
        for(int i_elt = 0; i_elt < n_elt; ++i_elt) {
          // int parent_elt = i_elt;
          if (parent_num != NULL) {
            parent_elt = parent_num[i_elt];
          }
          else {
            parent_elt++;
          }
          if(is_selected[i_part][parent_elt] != -1) {
            n_elmt_to_send     += 1;
            n_elmt_vtx_to_send += n_vtx_per_elmt;
          }
        }
      }
    } /* End section */

    elmt_vtx_n       [i_part] = malloc(    n_elmt_to_send      * sizeof(int                 ));
    elmt_type        [i_part] = malloc(    n_elmt_to_send      * sizeof(PDM_Mesh_nodal_elt_t));
    elmt_vtx         [i_part] = malloc(    n_elmt_vtx_to_send  * sizeof(PDM_g_num_t         ));
    vtx_init_location[i_part] = malloc(3 * n_elmt_vtx_to_send  * sizeof(int                 ));
    elmt_section_id  [i_part] = malloc(    n_elmt_to_send      * sizeof(int                 ));
    elmt_face_n      [i_part] = malloc(    n_elmt_to_send      * sizeof(int                 ));
    elmt_face_vtx_n  [i_part] = malloc(    n_elmt_face_to_send * sizeof(int                 ));
    elmt_face        [i_part] = malloc(    n_elmt_vtx_to_send  * sizeof(PDM_g_num_t         ));

    PDM_g_num_t* _vtx_ln_to_gn = extrp->vtx_ln_to_gn[i_part];

    /* Remplissage */
    parent_elt = -1;
    for(int i_section = 0; i_section < n_section; ++i_section) {
      int n_elt = PDM_part_mesh_nodal_elmts_section_n_elt_get(extrp->pmne, sections_id[i_section], i_part);
      PDM_Mesh_nodal_elt_t t_elt = PDM_part_mesh_nodal_elmts_section_type_get(extrp->pmne, sections_id[i_section]);

      int *parent_num = PDM_part_mesh_nodal_elmts_parent_num_get(extrp->pmne,
                                                                 sections_id[i_section],
                                                                 i_part,
                                                                 PDM_OWNERSHIP_KEEP);

      if (t_elt == PDM_MESH_NODAL_POLY_2D) {
        int *face_vtx_idx;
        int *face_vtx;
        PDM_part_mesh_nodal_elmts_section_poly2d_get(extrp->pmne,
                                                   sections_id[i_section],
                                                   i_part,
                                                   &face_vtx_idx,
                                                   &face_vtx,
                                                   PDM_OWNERSHIP_KEEP);

        /* Selection */
        for(int i_elt = 0; i_elt < n_elt; ++i_elt) {
          // int parent_elt = i_elt;
          if (parent_num != NULL) {
            parent_elt = parent_num[i_elt];
          }
          else {
            parent_elt++;
          }
          if(is_selected[i_part][parent_elt] != -1) {
            int idx = is_selected[i_part][parent_elt];

            elmt_type      [i_part][idx] = t_elt;
            elmt_vtx_n     [i_part][idx] = face_vtx_idx[i_elt+1] - face_vtx_idx[i_elt];
            elmt_face_n    [i_part][idx] = 0;
            elmt_section_id[i_part][idx] = i_section;

            if (0) {
              log_trace("extract polygon "PDM_FMT_G_NUM" : vtx ", entity_g_num[i_part][parent_elt]);
              for (int idx_vtx = face_vtx_idx[i_elt]; idx_vtx < face_vtx_idx[i_elt+1]; idx_vtx++) {
                log_trace(" "PDM_FMT_G_NUM, _vtx_ln_to_gn[face_vtx[idx_vtx]-1]);
              }
              log_trace("\n");
            }
          }
        }

      }
      else if (t_elt == PDM_MESH_NODAL_POLY_3D) {
        int *cell_face_idx;
        int *cell_face;
        int  n_face;
        int *face_vtx_idx;
        int *face_vtx;
        PDM_g_num_t *face_ln_to_gn    = NULL;
        int         *_parent_num      = NULL; // Il faut adpater tout le part_mesh_nodal_elmts
        PDM_g_num_t *elt_ln_to_gn     = NULL;
        PDM_g_num_t *parent_elt_g_num = NULL;
        PDM_part_mesh_nodal_elmts_section_poly3d_get(extrp->pmne,
                                                   sections_id[i_section],
                                                   i_part,
                                                   &n_face,
                                                   &face_ln_to_gn,
                                                   &face_vtx_idx,
                                                   &face_vtx,
                                                   &elt_ln_to_gn,
                                                   &cell_face_idx,
                                                   &cell_face,
                                                   &_parent_num,
                                                   &parent_elt_g_num,
                                                   PDM_OWNERSHIP_KEEP);

        /* Selection */
        for(int i_elt = 0; i_elt < n_elt; ++i_elt) {
          // int parent_elt = i_elt;
          if (parent_num != NULL) {
            parent_elt = parent_num[i_elt];
          }
          else {
            parent_elt++;
          }
          if(is_selected[i_part][parent_elt] != -1) {
            int idx = is_selected[i_part][parent_elt];

            int n_vtx_per_elmt = 0;
            for(int idx_face = cell_face_idx[i_elt]; idx_face < cell_face_idx[i_elt+1]; ++idx_face) {
              int i_face = PDM_ABS(cell_face[idx_face])-1;
              n_vtx_per_elmt += face_vtx_idx[i_face+1] - face_vtx_idx[i_face];
            }

            elmt_type      [i_part][idx] = t_elt;
            elmt_vtx_n     [i_part][idx] = n_vtx_per_elmt;
            elmt_face_n    [i_part][idx] = cell_face_idx[i_elt+1] - cell_face_idx[i_elt];
            elmt_section_id[i_part][idx] = i_section;
          }
        }
      }
      else {
        int         *elt_vtx          = NULL;
        int         *_parent_num      = NULL;
        PDM_g_num_t *elt_ln_to_gn     = NULL;
        PDM_g_num_t *parent_elt_g_num = NULL;
        int          order            = 0;
        const char  *ho_ordering      = NULL;

        PDM_part_mesh_nodal_elmts_section_std_ho_get(extrp->pmne,
                                                   sections_id[i_section],
                                                   i_part,
                                                   &elt_vtx,
                                                   &elt_ln_to_gn,
                                                   &_parent_num,
                                                   &parent_elt_g_num,
                                                   &order,
                                                   &ho_ordering,
                                                   PDM_OWNERSHIP_KEEP);

        int n_vtx_per_elmt = PDM_Mesh_nodal_n_vtx_elt_get(t_elt, order);

        /* Selection */
        for(int i_elt = 0; i_elt < n_elt; ++i_elt) {
          // int parent_elt = i_elt;
          if (parent_num != NULL) {
            parent_elt = parent_num[i_elt];
          }
          else {
            parent_elt++;
          }
          if(is_selected[i_part][parent_elt] != -1) {

            int idx = is_selected[i_part][parent_elt];

            elmt_type      [i_part][idx] = t_elt;
            elmt_vtx_n     [i_part][idx] = n_vtx_per_elmt;
            elmt_face_n    [i_part][idx] = 0;
            elmt_section_id[i_part][idx] = i_section;
          }
        }
      } /* End section */
    }

    /* Remplissage vtx */
    parent_elt = -1;
    int *elt_vtx_idx  = PDM_array_new_idx_from_sizes_int(elmt_vtx_n [i_part], n_elmt_to_send);
    int *elt_face_idx = PDM_array_new_idx_from_sizes_int(elmt_face_n[i_part], n_elmt_to_send);
    for(int i_section = 0; i_section < n_section; ++i_section) {
      int n_elt = PDM_part_mesh_nodal_elmts_section_n_elt_get(extrp->pmne, sections_id[i_section], i_part);
      PDM_Mesh_nodal_elt_t t_elt = PDM_part_mesh_nodal_elmts_section_type_get(extrp->pmne, sections_id[i_section]);

      int *parent_num = PDM_part_mesh_nodal_elmts_parent_num_get(extrp->pmne,
                                                                 sections_id[i_section],
                                                                 i_part,
                                                                 PDM_OWNERSHIP_KEEP);

      if (t_elt == PDM_MESH_NODAL_POLY_2D) {
        int *face_vtx_idx;
        int *face_vtx;
        PDM_part_mesh_nodal_elmts_section_poly2d_get(extrp->pmne,
                                                   sections_id[i_section],
                                                   i_part,
                                                   &face_vtx_idx,
                                                   &face_vtx,
                                                   PDM_OWNERSHIP_KEEP);

        /* Selection */
        for(int i_elt = 0; i_elt < n_elt; ++i_elt) {
          // int parent_elt = i_elt;
          if (parent_num != NULL) {
            parent_elt = parent_num[i_elt];
          }
          else {
            parent_elt++;
          }
          if(is_selected[i_part][parent_elt] != -1) {
            int idx = is_selected[i_part][parent_elt];

            int n_vtx_per_elmt = face_vtx_idx[i_elt+1] - face_vtx_idx[i_elt];

            int idx_read  = face_vtx_idx[i_elt];;
            int idx_write = elt_vtx_idx[idx];
            for(int k = 0; k < n_vtx_per_elmt; ++k) {
              elmt_vtx         [i_part][idx_write+k] = _vtx_ln_to_gn[face_vtx[idx_read+k]-1];
              vtx_init_location[i_part][3*(idx_write+k)  ] = i_rank;
              vtx_init_location[i_part][3*(idx_write+k)+1] = i_part;
              vtx_init_location[i_part][3*(idx_write+k)+2] = face_vtx[idx_read+k]-1;
            }

            // n_elmt_to_send     += 1;
            // n_elmt_vtx_to_send += n_vtx_per_elmt;
          }
        }

      }
      else if (t_elt == PDM_MESH_NODAL_POLY_3D) {
        int *cell_face_idx;
        int *cell_face;
        int  n_face;
        int *face_vtx_idx;
        int *face_vtx;
        PDM_g_num_t *face_ln_to_gn    = NULL;
        int         *_parent_num      = NULL; // Il faut adpater tout le part_mesh_nodal_elmts
        PDM_g_num_t *elt_ln_to_gn     = NULL;
        PDM_g_num_t *parent_elt_g_num = NULL;
        PDM_part_mesh_nodal_elmts_section_poly3d_get(extrp->pmne,
                                                   sections_id[i_section],
                                                   i_part,
                                                   &n_face,
                                                   &face_ln_to_gn,
                                                   &face_vtx_idx,
                                                   &face_vtx,
                                                   &elt_ln_to_gn,
                                                   &cell_face_idx,
                                                   &cell_face,
                                                   &_parent_num,
                                                   &parent_elt_g_num,
                                                   PDM_OWNERSHIP_KEEP);
        // PDM_log_trace_connectivity_int(face_vtx_idx, face_vtx, n_face, "face_vtx : ");
        if (0) {
          char filename[999];
          sprintf(filename, "check_faces_%d_%d_%d.vtk", i_part, i_section, i_rank);
          PDM_vtk_write_polydata(filename,
                                 extrp->n_vtx[i_part],
                                 extrp->pvtx_coord[i_part],
                                 NULL,
                                 n_face,
                                 face_vtx_idx,
                                 face_vtx,
                                 face_ln_to_gn,
                                 NULL);
        }

        /* Selection */
        for(int i_elt = 0; i_elt < n_elt; ++i_elt) {
          // int parent_elt = i_elt;
          if (parent_num != NULL) {
            parent_elt = parent_num[i_elt];
          }
          else {
            parent_elt++;
          }
          if(is_selected[i_part][parent_elt] != -1) {
            int idx = is_selected[i_part][parent_elt];
            // log_trace("idx = %d (i_elt %d, parent_elt %d, cell "PDM_FMT_G_NUM")\n",
            //           idx,
            //           i_elt, parent_elt,
            //           entity_g_num[i_part][parent_elt]);

            int idx_write_vtx  = elt_vtx_idx [idx];
            int idx_write_face = elt_face_idx[idx];
            for(int idx_face = cell_face_idx[i_elt]; idx_face < cell_face_idx[i_elt+1]; ++idx_face) {
              int i_face = PDM_ABS(cell_face[idx_face])-1;

              elmt_face      [i_part][idx_write_face] = PDM_SIGN(cell_face[idx_face]) * face_ln_to_gn[i_face];
              elmt_face_vtx_n[i_part][idx_write_face] = face_vtx_idx[i_face+1] - face_vtx_idx[i_face];
              // log_trace("  face %d ("PDM_FMT_G_NUM"), has %d vtx\n",
              //           i_face,
              //           elmt_face      [i_part][idx_write_face],
              //           elmt_face_vtx_n[i_part][idx_write_face]);
              idx_write_face++;

              for (int idx_vtx = face_vtx_idx[i_face]; idx_vtx < face_vtx_idx[i_face+1]; idx_vtx++) {
                elmt_vtx         [i_part][  idx_write_vtx  ] = _vtx_ln_to_gn[face_vtx[idx_vtx]-1];
                vtx_init_location[i_part][3*idx_write_vtx  ] = i_rank;
                vtx_init_location[i_part][3*idx_write_vtx+1] = i_part;
                vtx_init_location[i_part][3*idx_write_vtx+2] = face_vtx[idx_vtx]-1;
                idx_write_vtx++;
              }
            }

            // n_elmt_to_send     += 1;
            // n_elmt_vtx_to_send += n_vtx_per_elmt;
          }
        }


      }
      else {
        int         *elt_vtx          = NULL;
        int         *_parent_num      = NULL;
        PDM_g_num_t *elt_ln_to_gn     = NULL;
        PDM_g_num_t *parent_elt_g_num = NULL;
        int          order            = 0;
        const char  *ho_ordering      = NULL;

        PDM_part_mesh_nodal_elmts_section_std_ho_get(extrp->pmne,
                                                   sections_id[i_section],
                                                   i_part,
                                                   &elt_vtx,
                                                   &elt_ln_to_gn,
                                                   &_parent_num,
                                                   &parent_elt_g_num,
                                                   &order,
                                                   &ho_ordering,
                                                   PDM_OWNERSHIP_KEEP);

        int n_vtx_per_elmt = PDM_Mesh_nodal_n_vtx_elt_get(t_elt, order);

        int *ijk_to_user = NULL;
        if (order > 1 && ho_ordering != NULL) {
          ijk_to_user = PDM_ho_ordering_ijk_to_user_get(ho_ordering,
                                                        t_elt,
                                                        order);
        }

        /* Selection */
        for(int i_elt = 0; i_elt < n_elt; ++i_elt) {
          // int parent_elt = i_elt;
          if (parent_num != NULL) {
            parent_elt = parent_num[i_elt];
          }
          else {
            parent_elt++;
          }
          if(is_selected[i_part][parent_elt] != -1) {

            int idx = is_selected[i_part][parent_elt];

            int idx_read  = i_elt * n_vtx_per_elmt;
            int idx_write = elt_vtx_idx[idx];
            for(int k = 0; k < n_vtx_per_elmt; ++k) {
              int _k = k;
              if (ijk_to_user != NULL) {
                _k = ijk_to_user[k];
              }
              elmt_vtx         [i_part][idx_write+k] = _vtx_ln_to_gn[elt_vtx[idx_read+_k]-1];
              vtx_init_location[i_part][3*(idx_write+k)  ] = i_rank;
              vtx_init_location[i_part][3*(idx_write+k)+1] = i_part;
              vtx_init_location[i_part][3*(idx_write+k)+2] = elt_vtx[idx_read+k]-1;
            }

            // n_elmt_to_send     += 1;
            // n_elmt_vtx_to_send += n_vtx_per_elmt;
          }
        }
      }
    } /* End section */
    // PDM_log_trace_connectivity_long(elt_face_idx, elmt_face      [i_part], n_elmt_to_send, "elmt_face       : ");
    // PDM_log_trace_connectivity_int (elt_face_idx, elmt_face_vtx_n[i_part], n_elmt_to_send, "elmt_face_vtx_n : ");

    free(elt_vtx_idx);
    free(elt_face_idx);
  } /* End i_part */


  PDM_Mesh_nodal_elt_t **recv_elmt_type    = NULL;
  int                    request_elmt_type = -1;
  PDM_part_to_part_reverse_iexch(ptp,
                                 PDM_MPI_COMM_KIND_P2P,
                                 PDM_STRIDE_CST_INTERLACED,
                                 PDM_PART_TO_PART_DATA_DEF_ORDER_GNUM1_COME_FROM,
                                 1,
                                 sizeof(PDM_Mesh_nodal_elt_t),
                                 NULL,
                (const void **)  elmt_type,
                                 NULL,
                    (void ***)   &recv_elmt_type,
                                 &request_elmt_type);
  PDM_part_to_part_reverse_iexch_wait(ptp, request_elmt_type);

  int **recv_elmt_section_id    = NULL;
  int   request_elmt_section_id = -1;
  PDM_part_to_part_reverse_iexch(ptp,
                                 PDM_MPI_COMM_KIND_P2P,
                                 PDM_STRIDE_CST_INTERLACED,
                                 PDM_PART_TO_PART_DATA_DEF_ORDER_GNUM1_COME_FROM,
                                 1,
                                 sizeof(int),
                                 NULL,
                (const void **)  elmt_section_id,
                                 NULL,
                    (void ***)   &recv_elmt_section_id,
                                 &request_elmt_section_id);
  PDM_part_to_part_reverse_iexch_wait(ptp, request_elmt_section_id);

  int         **recv_elmt_vtx_n  = NULL;
  PDM_g_num_t **recv_elmt_vtx    = NULL;
  int           request_elmt_vtx = -1;
  PDM_part_to_part_reverse_iexch(ptp,
                                 PDM_MPI_COMM_KIND_P2P,
                                 PDM_STRIDE_VAR_INTERLACED,
                                 PDM_PART_TO_PART_DATA_DEF_ORDER_GNUM1_COME_FROM,
                                 -1,
                                 sizeof(PDM_g_num_t),
                (const int  **)  elmt_vtx_n,
                (const void **)  elmt_vtx,
                                 &recv_elmt_vtx_n,
                    (void ***)   &recv_elmt_vtx,
                                 &request_elmt_vtx);
  PDM_part_to_part_reverse_iexch_wait(ptp, request_elmt_vtx);

  int         **recv_elmt_face_n  = NULL;
  PDM_g_num_t **recv_elmt_face    = NULL;
  int           request_elmt_face = -1;
  PDM_part_to_part_reverse_iexch(ptp,
                                 PDM_MPI_COMM_KIND_P2P,
                                 PDM_STRIDE_VAR_INTERLACED,
                                 PDM_PART_TO_PART_DATA_DEF_ORDER_GNUM1_COME_FROM,
                                 -1,
                                 sizeof(PDM_g_num_t),
                (const int  **)  elmt_face_n,
                (const void **)  elmt_face,
                                 &recv_elmt_face_n,
                    (void ***)   &recv_elmt_face,
                                 &request_elmt_face);
  PDM_part_to_part_reverse_iexch_wait(ptp, request_elmt_face);

  int **recv_elmt_face_vtx_n = NULL;
  PDM_part_to_part_reverse_iexch(ptp,
                                 PDM_MPI_COMM_KIND_P2P,
                                 PDM_STRIDE_VAR_INTERLACED,
                                 PDM_PART_TO_PART_DATA_DEF_ORDER_GNUM1_COME_FROM,
                                 -1,
                                 sizeof(int),
                (const int  **)  elmt_face_n,
                (const void **)  elmt_face_vtx_n,
                                 &recv_elmt_face_n,
                    (void ***)   &recv_elmt_face_vtx_n,
                                 &request_elmt_face);
  PDM_part_to_part_reverse_iexch_wait(ptp, request_elmt_face);

  int **recv_vtx_init_location_n  = NULL;
  int **recv_vtx_init_location    = NULL;
  int   request_vtx_init_location = -1;
  PDM_part_to_part_reverse_iexch(ptp,
                                 PDM_MPI_COMM_KIND_P2P,
                                 PDM_STRIDE_VAR_INTERLACED,
                                 PDM_PART_TO_PART_DATA_DEF_ORDER_GNUM1_COME_FROM,
                                 -1,
                                 3 * sizeof(int),
                (const int  **)  elmt_vtx_n,
                (const void **)  vtx_init_location,
                                 &recv_vtx_init_location_n,
                    (void ***)   &recv_vtx_init_location,
                                 &request_vtx_init_location);
  PDM_part_to_part_reverse_iexch_wait(ptp, request_vtx_init_location);

  for(int i_part = 0; i_part < extrp->n_part_out; ++i_part) {
    free(recv_vtx_init_location_n[i_part]);
  }
  free(recv_vtx_init_location_n);


  /*
   * Free
   */
  for(int i_part = 0; i_part < extrp->n_part_in; ++i_part) {
    free(is_selected      [i_part]);
    free(elmt_vtx_n       [i_part]);
    free(elmt_type        [i_part]);
    free(elmt_vtx         [i_part]);
    free(elmt_section_id  [i_part]);
    free(vtx_init_location[i_part]);
    free(elmt_face_n      [i_part]);
    free(elmt_face_vtx_n  [i_part]);
    free(elmt_face        [i_part]);
  }
  free(is_selected);
  free(n_extract_vtx);
  free(elmt_vtx_n     );
  free(elmt_type      );
  free(elmt_vtx       );
  free(elmt_section_id);
  free(vtx_init_location);
  free(elmt_face_n);
  free(elmt_face_vtx_n);
  free(elmt_face      );

  /*
   * Second pass to create the new part_mesh_nodal
   */
  PDM_part_mesh_nodal_elmts_t* extract_pmne = PDM_part_mesh_nodal_elmts_create(extrp->pmne->mesh_dimension,
                                                                               extrp->n_part_out,//extrp->pmne->n_part, // == n_part_out
                                                                               extrp->pmne->comm);

  extrp->extract_pmne = extract_pmne;
  /*
   * Post-traitement
   */
  assert(extrp->pextract_n_entity[PDM_MESH_ENTITY_VTX] == NULL);
  extrp->pextract_n_entity              [PDM_MESH_ENTITY_VTX] = malloc(extrp->n_part_out * sizeof(int         *));
  extrp->pextract_entity_parent_ln_to_gn[PDM_MESH_ENTITY_VTX] = malloc(extrp->n_part_out * sizeof(PDM_g_num_t *));
  int **target_vtx_to_part1_vtx = malloc(extrp->n_part_out * sizeof(int *));

  for(int i_part = 0; i_part < extrp->n_part_out; ++i_part) {

    int n_tot_size = 0;
    for(int i = 0; i < extrp->n_target[i_part]; ++i) {
      n_tot_size += recv_elmt_vtx_n[i_part][i];
    }

    // PDM_log_trace_array_long(recv_elmt_section_id, n_tot_size, " :");
    // PDM_log_trace_array_long(recv_elmt_section_id[i_part], extrp->n_target[i_part], "recv_elmt_section_id :");

    int *unique_order_entity2 = malloc( n_tot_size * sizeof(int));
    int n_lextract_vtx = PDM_inplace_unique_long2(recv_elmt_vtx[i_part], unique_order_entity2, 0, n_tot_size-1);
    recv_elmt_vtx[i_part] = realloc(recv_elmt_vtx[i_part], n_lextract_vtx * sizeof(PDM_g_num_t));

    extrp->pextract_n_entity              [PDM_MESH_ENTITY_VTX][i_part] = n_lextract_vtx;
    extrp->pextract_entity_parent_ln_to_gn[PDM_MESH_ENTITY_VTX][i_part] = recv_elmt_vtx[i_part];

    // Bucket sort by sections id
    int *n_elmt_by_section      = PDM_array_zeros_int(n_section);
    int *s_elmt_vtx_by_section  = PDM_array_zeros_int(n_section);
    int *s_elmt_face_by_section = PDM_array_zeros_int(n_section);

    for(int i = 0; i < extrp->n_target[i_part]; ++i) {
      n_elmt_by_section[recv_elmt_section_id[i_part][i]]++;
      s_elmt_vtx_by_section [recv_elmt_section_id[i_part][i]] += recv_elmt_vtx_n [i_part][i];
      s_elmt_face_by_section[recv_elmt_section_id[i_part][i]] += recv_elmt_face_n[i_part][i];
    }

    if(0 == 1) {
      PDM_log_trace_array_int(n_elmt_by_section, n_section, "n_elmt_by_section ::");
    }

    int         **elmt_face_idx_by_section   = malloc(n_section * sizeof(int         *));
    int         **elmt_vtx_idx_by_section    = malloc(n_section * sizeof(int         *));
    int         **elmt_vtx_by_section        = malloc(n_section * sizeof(int         *));
    PDM_g_num_t **elmt_face_by_section       = malloc(n_section * sizeof(PDM_g_num_t *));
    int         **elmt_face_sign_by_section  = malloc(n_section * sizeof(int         *));
    int         **elmt_face_vtx_n_by_section = malloc(n_section * sizeof(int         *));
    int         **extract_parent_num         = malloc(n_section * sizeof(int         *));
    PDM_g_num_t **extract_parent_g_num       = malloc(n_section * sizeof(PDM_g_num_t *));
    for(int i_section = 0; i_section < n_section; ++i_section) {

      PDM_Mesh_nodal_elt_t t_elt = PDM_part_mesh_nodal_elmts_section_type_get(extrp->pmne, sections_id[i_section]);
      // int n_vtx_per_elmt         = PDM_Mesh_nodal_n_vtx_elt_get            (t_elt    , 1);

      // elmt_vtx_by_section [i_section  ] = malloc( n_vtx_per_elmt * n_elmt_by_section[i_section] * sizeof(int        ));
      if (t_elt == PDM_MESH_NODAL_POLY_2D) {
        elmt_vtx_idx_by_section[i_section] = malloc(sizeof(int) * (n_elmt_by_section[i_section] + 1));
        elmt_vtx_idx_by_section[i_section][0] = 0;
      }
      else if (t_elt == PDM_MESH_NODAL_POLY_3D) {
        elmt_face_idx_by_section  [i_section] = malloc(sizeof(int) * (n_elmt_by_section[i_section] + 1));
        elmt_face_idx_by_section  [i_section][0] = 0;
        elmt_face_by_section      [i_section] = malloc(s_elmt_face_by_section[i_section] * sizeof(PDM_g_num_t));
        elmt_face_sign_by_section [i_section] = malloc(s_elmt_face_by_section[i_section] * sizeof(int        ));
        elmt_face_vtx_n_by_section[i_section] = malloc(s_elmt_face_by_section[i_section] * sizeof(int        ));
        elmt_vtx_idx_by_section   [i_section] = malloc(sizeof(int) * (n_elmt_by_section[i_section] + 1));
        elmt_vtx_idx_by_section   [i_section][0] = 0;
      }
      elmt_vtx_by_section [i_section] = malloc(s_elmt_vtx_by_section[i_section] * sizeof(int        ));
      extract_parent_num  [i_section] = malloc(    n_elmt_by_section[i_section] * sizeof(int        ));
      extract_parent_g_num[i_section] = malloc(    n_elmt_by_section[i_section] * sizeof(PDM_g_num_t));

      n_elmt_by_section   [i_section] = 0;
    }
    free(s_elmt_vtx_by_section);

    // On reclasse tout les éléments
    int idx_read = 0;
    int idx_read_face = 0;
    for(int i = 0; i < extrp->n_target[i_part]; ++i) {
      int lsection_id     = recv_elmt_section_id[i_part][i];
      int n_vtx_per_elmt  = recv_elmt_vtx_n     [i_part][i];
      int n_face_per_elmt = recv_elmt_face_n    [i_part][i];

      int idx_write = n_elmt_by_section[lsection_id]++;

      extract_parent_num  [lsection_id][idx_write] = i;
      extract_parent_g_num[lsection_id][idx_write] = extrp->target_gnum[i_part][i];

      PDM_Mesh_nodal_elt_t t_elt = PDM_part_mesh_nodal_elmts_section_type_get(extrp->pmne,
                                                                            sections_id[lsection_id]);
      int idx_vtx;
      if (t_elt == PDM_MESH_NODAL_POLY_2D) {
        elmt_vtx_idx_by_section[lsection_id][idx_write+1] = elmt_vtx_idx_by_section[lsection_id][idx_write] + n_vtx_per_elmt;
        idx_vtx = elmt_vtx_idx_by_section[lsection_id][idx_write];
      }
      else if (t_elt == PDM_MESH_NODAL_POLY_3D) {
        // log_trace("lsection_id = %d, type = %d\n", lsection_id, t_elt);
        elmt_face_idx_by_section[lsection_id][idx_write+1] = elmt_face_idx_by_section[lsection_id][idx_write] + n_face_per_elmt;
        idx_vtx = elmt_vtx_idx_by_section[lsection_id][idx_write];
        elmt_vtx_idx_by_section[lsection_id][idx_write+1] = idx_vtx + n_vtx_per_elmt;
        int idx0 = elmt_face_idx_by_section[lsection_id][idx_write];
        // log_trace("cell "PDM_FMT_G_NUM", idx_write = %d, idx0 = %d, idx_vtx = %d\n",
        //           extrp->target_gnum[i_part][i], idx_write, idx0, idx_vtx);
        for (int j = 0; j < n_face_per_elmt; j++) {
          PDM_g_num_t face_gnum = recv_elmt_face[i_part][idx_read_face];
          // log_trace("  face "PDM_FMT_G_NUM", %d vtx\n",
          //           face_gnum, recv_elmt_face_vtx_n[i_part][idx_read_face]);
          elmt_face_by_section      [lsection_id][idx0+j] = PDM_ABS (face_gnum);
          elmt_face_sign_by_section [lsection_id][idx0+j] = PDM_SIGN(face_gnum);
          elmt_face_vtx_n_by_section[lsection_id][idx0+j] = recv_elmt_face_vtx_n[i_part][idx_read_face];
          // log_trace("  idx_read_face = %d, gnum = "PDM_FMT_G_NUM", vtx_n = %d\n",
          //           idx_read_face,
          //           face_gnum,
          //           recv_elmt_face_vtx_n[i_part][idx_read_face]);
          idx_read_face++;
        }
      }
      else {
        idx_vtx = n_vtx_per_elmt*idx_write;
      }

      for(int j = 0; j < n_vtx_per_elmt; ++j) {
        int l_elmt     = unique_order_entity2[idx_read++];
        elmt_vtx_by_section[lsection_id][idx_vtx++] = (l_elmt+1);
      }
    }

    /*
     * Prepare vtx_init_location
     */
    target_vtx_to_part1_vtx[i_part] = malloc(3 * n_lextract_vtx * sizeof(int));
    for(int i = 0; i < n_tot_size; ++i) {
      int l_elmt = unique_order_entity2[i];
      // C'est maybe ecraser plusieurs fois
      target_vtx_to_part1_vtx[i_part][3*l_elmt  ] = recv_vtx_init_location[i_part][3*i  ];
      target_vtx_to_part1_vtx[i_part][3*l_elmt+1] = recv_vtx_init_location[i_part][3*i+1];
      target_vtx_to_part1_vtx[i_part][3*l_elmt+2] = recv_vtx_init_location[i_part][3*i+2];
    }

    /*
     * Fill up structure
     */
    for(int i_section = 0; i_section < n_section; ++i_section) {

      PDM_Mesh_nodal_elt_t t_elt = PDM_part_mesh_nodal_elmts_section_type_get(extrp->pmne, sections_id[i_section]);

      int extract_section_id = PDM_part_mesh_nodal_elmts_add(extract_pmne, t_elt);

      if (t_elt == PDM_MESH_NODAL_POLY_2D) {
        PDM_part_mesh_nodal_elmts_section_poly2d_set(extract_pmne,
                                                   extract_section_id,
                                                   i_part,
                                                   n_elmt_by_section[i_section],
                                                   elmt_vtx_idx_by_section[i_section],
                                                   elmt_vtx_by_section[i_section],
                                                   NULL,
                                                   extract_parent_num[i_section],
                                                   // extract_parent_g_num[i_section],
                                                   PDM_OWNERSHIP_KEEP);
        if (0) {
          for (int i_elt = 0; i_elt < n_elmt_by_section[i_section]; i_elt++) {
            log_trace("        polygon "PDM_FMT_G_NUM" : vtx ", extract_parent_g_num[i_section][i_elt]);
            for (int idx_vtx = elmt_vtx_idx_by_section[i_section][i_elt]; idx_vtx < elmt_vtx_idx_by_section[i_section][i_elt+1]; idx_vtx++) {
              log_trace(" "PDM_FMT_G_NUM,
                        extrp->pextract_entity_parent_ln_to_gn[PDM_MESH_ENTITY_VTX][i_part][elmt_vtx_by_section[i_section][idx_vtx]-1]);
            }
            log_trace("\n");
          }
        }
        free(extract_parent_g_num[i_section]);// pass to extract_pmne?
      }
      else if (t_elt == PDM_MESH_NODAL_POLY_3D) {

        int *unique_order_face = malloc(sizeof(int) * s_elmt_face_by_section[i_section]);

        PDM_g_num_t *tmp = malloc(sizeof(PDM_g_num_t) * s_elmt_face_by_section[i_section]);
        memcpy(tmp, elmt_face_by_section[i_section],
               sizeof(PDM_g_num_t) * s_elmt_face_by_section[i_section]);

        int n_lextract_face = PDM_inplace_unique_long2(tmp,//elmt_face_by_section[i_section],
                                                       unique_order_face,
                                                       0,
                                                       s_elmt_face_by_section[i_section]-1);
        free(tmp);
        // PDM_log_trace_connectivity_long(elmt_face_idx_by_section[i_section],
        //                                 elmt_face_by_section[i_section],
        //                                 n_elmt_by_section[i_section],
        //                                 "elmt_face_by_section : ");
        // PDM_log_trace_array_long(elmt_face_by_section[i_section], s_elmt_face_by_section[i_section], "elmt_face_by_section : ");
        // PDM_log_trace_array_int (unique_order_face, s_elmt_face_by_section[i_section], "unique_order_face : ");

        int *cell_face = malloc(sizeof(int) * elmt_face_idx_by_section[i_section][n_elmt_by_section[i_section]]);
        PDM_g_num_t *face_ln_to_gn = malloc(sizeof(PDM_g_num_t) * n_lextract_face);

        int *face_vtx_idx = malloc(sizeof(int) * (n_lextract_face + 1));
        face_vtx_idx[0] = 0;
        for (int i = 0; i < n_elmt_by_section[i_section]; i++) {
          for (int j = elmt_face_idx_by_section[i_section][i]; j < elmt_face_idx_by_section[i_section][i+1]; j++) {
            face_vtx_idx[unique_order_face[j]+1] = elmt_face_vtx_n_by_section[i_section][j];
          }
        }

        for (int i = 0; i < n_lextract_face; i++) {
          face_vtx_idx[i+1] += face_vtx_idx[i];
        }

        int *face_vtx = malloc(sizeof(int) * face_vtx_idx[n_lextract_face]);

        int idx = 0;
        for (int i = 0; i < n_elmt_by_section[i_section]; i++) {
          // log_trace("cell %d ("PDM_FMT_G_NUM")\n", i,
          //           extrp->target_gnum[i_part][extract_parent_num[i_section][i]]);
          idx = elmt_vtx_idx_by_section[i_section][i];
          for (int j = elmt_face_idx_by_section[i_section][i]; j < elmt_face_idx_by_section[i_section][i+1]; j++) {
            int face_id = unique_order_face[j];
            cell_face[j] = elmt_face_sign_by_section[i_section][j] * (face_id+1);
            face_ln_to_gn[face_id] = elmt_face_by_section[i_section][j];
            int face_vtx_n = face_vtx_idx[face_id+1] - face_vtx_idx[face_id];
            // log_trace("  j = %d, face %d ("PDM_FMT_G_NUM"), %d vtx, idx = %d\n",
            //           j, cell_face[j], face_ln_to_gn[face_id], face_vtx_n, idx);

            for (int k = 0; k < face_vtx_n; k++) {
              // log_trace("    idx = %d ("PDM_FMT_G_NUM" -> %d)\n",
              //           idx, elmt_vtx_by_section[i_section][idx], unique_order_entity2[idx]);
              face_vtx[face_vtx_idx[face_id]+k] = elmt_vtx_by_section[i_section][idx++];
            }
          }
        }
        free(unique_order_face);
        free(elmt_face_by_section      [i_section]);
        free(elmt_face_sign_by_section [i_section]);
        free(elmt_face_vtx_n_by_section[i_section]);
        free(elmt_vtx_by_section[i_section]);
        free(elmt_vtx_idx_by_section[i_section]);
        // PDM_log_trace_connectivity_int(face_vtx_idx,
        //                                face_vtx,
        //                                n_lextract_face,
        //                                "face_vtx : ");

        PDM_part_mesh_nodal_elmts_section_poly3d_set(extract_pmne,
                                                     extract_section_id,
                                                     i_part,
                                                     n_elmt_by_section[i_section],
                                                     n_lextract_face,
                                                     face_vtx_idx,
                                                     face_vtx,
                                                     face_ln_to_gn,
                                                     elmt_face_idx_by_section[i_section],
                                                     cell_face,
                                                     NULL,
                                                     extract_parent_num  [i_section],
                                                     extract_parent_g_num[i_section],
                                                     PDM_OWNERSHIP_KEEP);

        // if (1) {
        //   // Check cell-vtx connectivity
        //   int *cell_vtx_idx;
        //   int *cell_vtx;
        //   PDM_part_mesh_nodal_elmts_section_poly3d_cell_vtx_connect_get(extract_pmne,
        //                                                               extract_section_id,
        //                                                               i_part,
        //                                                               &cell_vtx_idx,
        //                                                               &cell_vtx);

        //   log_trace("---Après---\n");
        //   for (int i = 0; i < n_elmt_by_section[i_section]; i++) {
        //     // log_trace("cell "PDM_FMT_G_NUM"\n", ?);
        //     log_trace("cell %d: ", i);
        //     for (int j = cell_vtx_idx[i]; j < cell_vtx_idx[i+1]; j++) {
        //       int vtx_id = cell_vtx[j] - 1;
        //       log_trace(PDM_FMT_G_NUM" ", recv_elmt_vtx[i_part][vtx_id]);
        //     }
        //     log_trace("\n");
        //   }
        // }

        // free(extract_parent_g_num[i_section]);// pass to extract_pmne?
      }
      else {
        if (PDM_Mesh_nodal_elmt_is_ho(t_elt)) {
          int order = extrp->pmne->sections_std[sections_id[i_section]]->order;
          const char *ho_ordering = extrp->pmne->sections_std[sections_id[i_section]]->ho_ordering;
          PDM_part_mesh_nodal_elmts_std_ho_set(extract_pmne,
                                               extract_section_id,
                                               i_part,
                                               n_elmt_by_section[i_section],
                                               elmt_vtx_by_section[i_section],
                                               NULL,
                                               extract_parent_num  [i_section],
                                               extract_parent_g_num[i_section],
                                               order,
                                               ho_ordering,
                                               PDM_OWNERSHIP_KEEP);

        }
        else {
          PDM_part_mesh_nodal_elmts_std_set(extract_pmne,
                                            extract_section_id,
                                            i_part,
                                            n_elmt_by_section[i_section],
                                            elmt_vtx_by_section[i_section],
                                            NULL,
                                            extract_parent_num  [i_section],
                                            extract_parent_g_num[i_section],
                                            PDM_OWNERSHIP_KEEP);
        }
      }

    }

    free(n_elmt_by_section);
    free(elmt_vtx_idx_by_section);
    free(elmt_vtx_by_section);
    free(elmt_face_idx_by_section);
    free(elmt_face_by_section);
    free(elmt_face_sign_by_section);
    free(elmt_face_vtx_n_by_section);
    free(extract_parent_num);
    free(extract_parent_g_num);
    free(unique_order_entity2);
    free(s_elmt_face_by_section);

    free(recv_elmt_section_id  [i_part]);
    free(recv_elmt_type        [i_part]);
    free(recv_elmt_vtx_n       [i_part]);
    free(recv_vtx_init_location[i_part]);
    free(recv_elmt_face_n      [i_part]);
    free(recv_elmt_face        [i_part]);
    free(recv_elmt_face_vtx_n  [i_part]);
  }

  free(recv_elmt_section_id  );
  free(recv_elmt_type        );
  free(recv_elmt_vtx         );
  free(recv_elmt_vtx_n       );
  free(recv_vtx_init_location);
  free(recv_elmt_face_n);
  free(recv_elmt_face);
  free(recv_elmt_face_vtx_n);

  /*
   * Vtx only
   */
  int **part2_vtx_to_part1_vtx_idx = malloc(extrp->n_part_out * sizeof(int *));
  for(int i_part = 0; i_part < extrp->n_part_out; ++i_part) {
    int n_vtx = extrp->pextract_n_entity[PDM_MESH_ENTITY_VTX][i_part];
    part2_vtx_to_part1_vtx_idx[i_part] =  PDM_array_new_idx_from_const_stride_int(3, n_vtx);;
  }

  PDM_part_to_part_t* ptp_vtx    = NULL;
  ptp_vtx = PDM_part_to_part_create_from_num2_triplet((const PDM_g_num_t **) extrp->pextract_entity_parent_ln_to_gn[PDM_MESH_ENTITY_VTX],
                                                      (const int          *) extrp->pextract_n_entity              [PDM_MESH_ENTITY_VTX],
                                                      extrp->n_part_out,
                                                      extrp->n_vtx,
                                                      extrp->n_part_in,
                                                      (const int **) part2_vtx_to_part1_vtx_idx,
                                                      NULL,
                                                      (const int **) target_vtx_to_part1_vtx,
                                                      extrp->comm);
  extrp->ptp_entity[PDM_MESH_ENTITY_VTX] = ptp_vtx;

  int           exch_request = -1;
  PDM_part_to_part_reverse_iexch(ptp_vtx,
                                 PDM_MPI_COMM_KIND_P2P,
                                 PDM_STRIDE_CST_INTERLACED,
                                 PDM_PART_TO_PART_DATA_DEF_ORDER_PART2,
                                 3,
                                 sizeof(double),
                                 NULL,
                (const void **)  extrp->pvtx_coord,
                                 NULL,
                    (void ***)   &extrp->pextract_vtx_coord,
                                 &exch_request);
  PDM_part_to_part_reverse_iexch_wait(ptp_vtx, exch_request);

  /*
   * Free
   */
  for(int i_part = 0; i_part < extrp->n_part_out; ++i_part) {

    if(part2_vtx_to_part1_vtx_idx[i_part] != NULL) {
      free(part2_vtx_to_part1_vtx_idx[i_part]);
    }

    if(target_vtx_to_part1_vtx[i_part] != NULL) {
      free(target_vtx_to_part1_vtx[i_part]);
    }
  }
  free(part2_vtx_to_part1_vtx_idx);
  free(target_vtx_to_part1_vtx);


  // if(ptp_vtx != NULL) {
  //   PDM_part_to_part_free(ptp_vtx);
  // }

}

static
void
_extract_part_and_reequilibrate_nodal
(
  PDM_extract_part_t        *extrp
)
{
  PDM_UNUSED(extrp);
  abort();
}

static
void
_extract_part
(
  PDM_extract_part_t        *extrp
)
{
  int                  *pn_entity    = 0;
  PDM_g_num_t         **entity_g_num = NULL;
  PDM_mesh_entities_t   entity_type;
  if (extrp->dim == 3) {
    entity_type  = PDM_MESH_ENTITY_CELL;
    pn_entity    = extrp->n_cell;
    entity_g_num = extrp->cell_ln_to_gn;
  }
  else if (extrp->dim == 2) {
    entity_type  = PDM_MESH_ENTITY_FACE;
    pn_entity    = extrp->n_face;
    entity_g_num = extrp->face_ln_to_gn;
  }
  else if (extrp->dim == 1) {
    entity_type  = PDM_MESH_ENTITY_EDGE;
    pn_entity    = extrp->n_edge;
    entity_g_num = extrp->edge_ln_to_gn;
  }
  else {
    entity_type  = PDM_MESH_ENTITY_VTX;
    pn_entity    = extrp->n_vtx;
    entity_g_num = extrp->vtx_ln_to_gn;
  }
  PDM_UNUSED(pn_entity);

  /*
   *  Create array selected in gnum
   */
  PDM_g_num_t** entity_extract_g_num = NULL;
  PDM_g_num_t** child_selected_g_num = NULL;

  // PDM_log_trace_array_long(entity_g_num[0], pn_entity[0], "entity_g_num : ");

  _extract_gnum_and_compute_child(extrp->comm,
                                  extrp->compute_child_gnum,
                                  extrp->n_part_in,
                                  extrp->n_extract,
                                  entity_g_num,
                                  extrp->extract_lnum,
                                  &entity_extract_g_num,
                                  &child_selected_g_num);

  extrp->pextract_n_entity[entity_type] = malloc(extrp->n_part_out * sizeof(int));
  for(int i_part = 0; i_part < extrp->n_part_in; ++i_part) {
    extrp->pextract_n_entity[entity_type][i_part] = extrp->n_extract[i_part];
  }
  extrp->pextract_entity_ln_to_gn       [entity_type] = child_selected_g_num;
  extrp->pextract_entity_parent_ln_to_gn[entity_type] = entity_extract_g_num;


  /*
   * Extraction des connectivités
   */
  int from_face_edge = 0;
  int from_face_vtx  = 0;
  for(int i_part = 0; i_part < extrp->n_part_in; ++i_part) {
    if(extrp->pface_edge    [i_part] != NULL) {
      from_face_edge = 1;
    }
    if(extrp->pface_vtx    [i_part] != NULL) {
      from_face_vtx = 1;
    }
  }

  if (extrp->dim == 3) {
    extract_and_local_renum_entity1_entity2(extrp->comm,
                                            extrp->compute_child_gnum,
                                            extrp->n_part_in,
                                            extrp->n_cell,
                                            extrp->n_face,
                                            extrp->n_extract,
                                            extrp->extract_lnum,
                                            extrp->pcell_face_idx,
                                            extrp->pcell_face,
                                            extrp->face_ln_to_gn,
                                            &extrp->pextract_n_entity              [PDM_MESH_ENTITY_FACE],
                                            &extrp->pextract_connectivity_idx      [PDM_CONNECTIVITY_TYPE_CELL_FACE],
                                            &extrp->pextract_connectivity          [PDM_CONNECTIVITY_TYPE_CELL_FACE],
                                            &extrp->pextract_entity_ln_to_gn       [PDM_MESH_ENTITY_FACE],
                                            &extrp->pextract_entity_parent_ln_to_gn[PDM_MESH_ENTITY_FACE],
                                            &extrp->pextract_entity_parent_lnum    [PDM_MESH_ENTITY_FACE]);

    if(from_face_edge == 1) {

      extract_and_local_renum_entity1_entity2(extrp->comm,
                                              extrp->compute_child_gnum,
                                              extrp->n_part_in,
                                              extrp->n_face,
                                              extrp->n_edge,
                                              extrp->pextract_n_entity               [PDM_MESH_ENTITY_FACE],
                                              extrp->pextract_entity_parent_lnum     [PDM_MESH_ENTITY_FACE],
                                              extrp->pface_edge_idx,
                                              extrp->pface_edge,
                                              extrp->edge_ln_to_gn,
                                              &extrp->pextract_n_entity              [PDM_MESH_ENTITY_EDGE],
                                              &extrp->pextract_connectivity_idx      [PDM_CONNECTIVITY_TYPE_FACE_EDGE],
                                              &extrp->pextract_connectivity          [PDM_CONNECTIVITY_TYPE_FACE_EDGE],
                                              &extrp->pextract_entity_ln_to_gn       [PDM_MESH_ENTITY_EDGE],
                                              &extrp->pextract_entity_parent_ln_to_gn[PDM_MESH_ENTITY_EDGE],
                                              &extrp->pextract_entity_parent_lnum    [PDM_MESH_ENTITY_EDGE]);

      int **pedge_vtx_idx = (int **) malloc(extrp->n_part_in * sizeof(int *));
      for(int i_part = 0; i_part < extrp->n_part_in; ++i_part) {
        pedge_vtx_idx[i_part] = (int *) malloc((extrp->n_edge[i_part]+1) * sizeof(int));
        pedge_vtx_idx[i_part][0] = 0;
        for(int i_edge = 0; i_edge < extrp->n_edge[i_part]; ++i_edge) {
          pedge_vtx_idx[i_part][i_edge+1] = pedge_vtx_idx[i_part][i_edge] + 2;
        }
      }

      extract_and_local_renum_entity1_entity2(extrp->comm,
                                              extrp->compute_child_gnum,
                                              extrp->n_part_in,
                                              extrp->n_edge,
                                              extrp->n_vtx,
                                              extrp->pextract_n_entity               [PDM_MESH_ENTITY_EDGE],
                                              extrp->pextract_entity_parent_lnum     [PDM_MESH_ENTITY_EDGE],
                                              pedge_vtx_idx,
                                              extrp->pedge_vtx,
                                              extrp->vtx_ln_to_gn,
                                              &extrp->pextract_n_entity              [PDM_MESH_ENTITY_VTX],
                                              &extrp->pextract_connectivity_idx      [PDM_CONNECTIVITY_TYPE_EDGE_VTX],
                                              &extrp->pextract_connectivity          [PDM_CONNECTIVITY_TYPE_EDGE_VTX],
                                              &extrp->pextract_entity_ln_to_gn       [PDM_MESH_ENTITY_VTX],
                                              &extrp->pextract_entity_parent_ln_to_gn[PDM_MESH_ENTITY_VTX],
                                              &extrp->pextract_entity_parent_lnum    [PDM_MESH_ENTITY_VTX]);

      for(int i_part = 0; i_part < extrp->n_part_in; ++i_part) {
        free(pedge_vtx_idx   [i_part]);
      }
      free(pedge_vtx_idx      );

    } else if(from_face_vtx == 1){

      extract_and_local_renum_entity1_entity2(extrp->comm,
                                              extrp->compute_child_gnum,
                                              extrp->n_part_in,
                                              extrp->n_face,
                                              extrp->n_vtx,
                                              extrp->pextract_n_entity               [PDM_MESH_ENTITY_FACE],
                                              extrp->pextract_entity_parent_lnum     [PDM_MESH_ENTITY_FACE],
                                              extrp->pface_vtx_idx,
                                              extrp->pface_vtx,
                                              extrp->vtx_ln_to_gn,
                                              &extrp->pextract_n_entity              [PDM_MESH_ENTITY_VTX],
                                              &extrp->pextract_connectivity_idx      [PDM_CONNECTIVITY_TYPE_FACE_VTX],
                                              &extrp->pextract_connectivity          [PDM_CONNECTIVITY_TYPE_FACE_VTX],
                                              &extrp->pextract_entity_ln_to_gn       [PDM_MESH_ENTITY_VTX],
                                              &extrp->pextract_entity_parent_ln_to_gn[PDM_MESH_ENTITY_VTX],
                                              &extrp->pextract_entity_parent_lnum    [PDM_MESH_ENTITY_VTX]);
    }
  }
  else if (extrp->dim == 2) {

    if(from_face_edge == 1) {
      extract_and_local_renum_entity1_entity2(extrp->comm,
                                              extrp->compute_child_gnum,
                                              extrp->n_part_in,
                                              extrp->n_face,
                                              extrp->n_edge,
                                              extrp->n_extract,
                                              extrp->extract_lnum,
                                              extrp->pface_edge_idx,
                                              extrp->pface_edge,
                                              extrp->edge_ln_to_gn,
                                              &extrp->pextract_n_entity              [PDM_MESH_ENTITY_EDGE],
                                              &extrp->pextract_connectivity_idx      [PDM_CONNECTIVITY_TYPE_FACE_EDGE],
                                              &extrp->pextract_connectivity          [PDM_CONNECTIVITY_TYPE_FACE_EDGE],
                                              &extrp->pextract_entity_ln_to_gn       [PDM_MESH_ENTITY_EDGE],
                                              &extrp->pextract_entity_parent_ln_to_gn[PDM_MESH_ENTITY_EDGE],
                                              &extrp->pextract_entity_parent_lnum    [PDM_MESH_ENTITY_EDGE]);
      int **pedge_vtx_idx = (int **) malloc(extrp->n_part_in * sizeof(int *));
      for(int i_part = 0; i_part < extrp->n_part_in; ++i_part) {
        pedge_vtx_idx[i_part] = (int *) malloc((extrp->n_edge[i_part]+1) * sizeof(int));
        pedge_vtx_idx[i_part][0] = 0;
        for(int i_edge = 0; i_edge < extrp->n_edge[i_part]; ++i_edge) {
          pedge_vtx_idx[i_part][i_edge+1] = pedge_vtx_idx[i_part][i_edge] + 2;
        }
      }

      extract_and_local_renum_entity1_entity2(extrp->comm,
                                              extrp->compute_child_gnum,
                                              extrp->n_part_in,
                                              extrp->n_edge,
                                              extrp->n_vtx,
                                              extrp->pextract_n_entity               [PDM_MESH_ENTITY_EDGE],
                                              extrp->pextract_entity_parent_lnum     [PDM_MESH_ENTITY_EDGE],
                                              pedge_vtx_idx,
                                              extrp->pedge_vtx,
                                              extrp->vtx_ln_to_gn,
                                              &extrp->pextract_n_entity              [PDM_MESH_ENTITY_VTX],
                                              &extrp->pextract_connectivity_idx      [PDM_CONNECTIVITY_TYPE_EDGE_VTX],
                                              &extrp->pextract_connectivity          [PDM_CONNECTIVITY_TYPE_EDGE_VTX],
                                              &extrp->pextract_entity_ln_to_gn       [PDM_MESH_ENTITY_VTX],
                                              &extrp->pextract_entity_parent_ln_to_gn[PDM_MESH_ENTITY_VTX],
                                              &extrp->pextract_entity_parent_lnum    [PDM_MESH_ENTITY_VTX]);

      for(int i_part = 0; i_part < extrp->n_part_in; ++i_part) {
        free(pedge_vtx_idx   [i_part]);
      }
      free(pedge_vtx_idx      );

    } else if(from_face_vtx == 1) {
      extract_and_local_renum_entity1_entity2(extrp->comm,
                                              extrp->compute_child_gnum,
                                              extrp->n_part_in,
                                              extrp->n_face,
                                              extrp->n_vtx,
                                              extrp->n_extract,
                                              extrp->extract_lnum,
                                              extrp->pface_vtx_idx,
                                              extrp->pface_vtx,
                                              extrp->vtx_ln_to_gn,
                                              &extrp->pextract_n_entity              [PDM_MESH_ENTITY_VTX],
                                              &extrp->pextract_connectivity_idx      [PDM_CONNECTIVITY_TYPE_FACE_VTX],
                                              &extrp->pextract_connectivity          [PDM_CONNECTIVITY_TYPE_FACE_VTX],
                                              &extrp->pextract_entity_ln_to_gn       [PDM_MESH_ENTITY_VTX],
                                              &extrp->pextract_entity_parent_ln_to_gn[PDM_MESH_ENTITY_VTX],
                                              &extrp->pextract_entity_parent_lnum    [PDM_MESH_ENTITY_VTX]);
    }
  }
  else if (extrp->dim == 1) {
    int **pedge_vtx_idx = malloc(sizeof(int *) * extrp->n_part_in);
    for (int i = 0; i < extrp->n_part_in; i++) {
      pedge_vtx_idx[i] = PDM_array_new_idx_from_const_stride_int(2, extrp->n_edge[i]);
    }
    extract_and_local_renum_entity1_entity2(extrp->comm,
                                            extrp->compute_child_gnum,
                                            extrp->n_part_in,
                                            extrp->n_edge,
                                            extrp->n_vtx,
                                            extrp->n_extract,
                                            extrp->extract_lnum,
                                            pedge_vtx_idx,
                                            extrp->pedge_vtx,
                                            extrp->vtx_ln_to_gn,
                                            &extrp->pextract_n_entity              [PDM_MESH_ENTITY_VTX],
                                            &extrp->pextract_connectivity_idx      [PDM_CONNECTIVITY_TYPE_EDGE_VTX],
                                            &extrp->pextract_connectivity          [PDM_CONNECTIVITY_TYPE_EDGE_VTX],
                                            &extrp->pextract_entity_ln_to_gn       [PDM_MESH_ENTITY_VTX],
                                            &extrp->pextract_entity_parent_ln_to_gn[PDM_MESH_ENTITY_VTX],
                                            &extrp->pextract_entity_parent_lnum    [PDM_MESH_ENTITY_VTX]);
    for (int i = 0; i < extrp->n_part_in; i++) {
      free(pedge_vtx_idx[i]);
    }
    free(pedge_vtx_idx);
  }
  else {
    int **pvtx_vtx_idx = malloc(sizeof(int *) * extrp->n_part_in);
    int **pvtx_vtx     = malloc(sizeof(int *) * extrp->n_part_in);
    for (int i = 0; i < extrp->n_part_in; i++) {
      pvtx_vtx_idx[i] = PDM_array_new_idx_from_const_stride_int(1, extrp->n_vtx[i]);
      pvtx_vtx    [i] = pvtx_vtx_idx[i] + 1;
    }

    int **pextract_connectivity_idx = NULL;
    int **pextract_connectivity     = NULL;
    extract_and_local_renum_entity1_entity2(extrp->comm,
                                            extrp->compute_child_gnum,
                                            extrp->n_part_in,
                                            extrp->n_vtx,
                                            extrp->n_vtx,
                                            extrp->n_extract,
                                            extrp->extract_lnum,
                                            pvtx_vtx_idx,
                                            pvtx_vtx,
                                            extrp->vtx_ln_to_gn,
                                            &extrp->pextract_n_entity              [PDM_MESH_ENTITY_VTX],
                                            &pextract_connectivity_idx,// &extrp->pextract_connectivity_idx      [PDM_CONNECTIVITY_TYPE_EDGE_VTX],
                                            &pextract_connectivity,// &extrp->pextract_connectivity          [PDM_CONNECTIVITY_TYPE_EDGE_VTX],
                                            &extrp->pextract_entity_ln_to_gn       [PDM_MESH_ENTITY_VTX],
                                            &extrp->pextract_entity_parent_ln_to_gn[PDM_MESH_ENTITY_VTX],
                                            &extrp->pextract_entity_parent_lnum    [PDM_MESH_ENTITY_VTX]);
    for (int i = 0; i < extrp->n_part_in; i++) {
      free(pvtx_vtx_idx[i]);
      free(pextract_connectivity_idx[i]);
    }
    free(pvtx_vtx_idx);
    free(pvtx_vtx);
    free(pextract_connectivity_idx);

    // PDM_error(__FILE__, __LINE__, 0, "Not yet implemented for dimension %d\n", extrp->dim);
  }

  int  *n_extract_vtx    = extrp->pextract_n_entity          [PDM_MESH_ENTITY_VTX];
  int **extract_vtx_lnum = extrp->pextract_entity_parent_lnum[PDM_MESH_ENTITY_VTX];

  /*
   * Exchange coordinates
   */
  extrp->pextract_vtx_coord = (double **) malloc( extrp->n_part_in * sizeof(double *));
  for(int i_part = 0; i_part < extrp->n_part_in; ++i_part) {
    // log_trace("n_extract_vtx[%d] = %d\n", i_part, n_extract_vtx[i_part]);
    extrp->pextract_vtx_coord[i_part] = (double *) malloc( 3 * n_extract_vtx[i_part] * sizeof(double));

    for(int idx_vtx = 0; idx_vtx < n_extract_vtx[i_part]; ++idx_vtx) {
      int i_vtx = extract_vtx_lnum[i_part][idx_vtx];
      extrp->pextract_vtx_coord[i_part][3*idx_vtx  ] = extrp->pvtx_coord[i_part][3*i_vtx  ];
      extrp->pextract_vtx_coord[i_part][3*idx_vtx+1] = extrp->pvtx_coord[i_part][3*i_vtx+1];
      extrp->pextract_vtx_coord[i_part][3*idx_vtx+2] = extrp->pvtx_coord[i_part][3*i_vtx+2];
    }
  }

  /*
  * Extract face groups for test
  */
  /*if (extrp->n_group[PDM_BOUND_TYPE_FACE] > 0){
    int **part2_face_to_part1_face_idx = (int **) malloc( extrp->n_part_out * sizeof(int * ));

    for(int i_part = 0; i_part < extrp->n_part_out; ++i_part) {
      int n_face = extrp->pextract_n_entity[PDM_MESH_ENTITY_FACE][i_part];
      part2_face_to_part1_face_idx[i_part] = PDM_array_new_idx_from_const_stride_int(1, n_face);;
    }

    PDM_part_to_part_t *ptp_fac = PDM_part_to_part_create((const PDM_g_num_t **) extrp->pextract_entity_parent_ln_to_gn[PDM_MESH_ENTITY_FACE],
                                                          (const int          *) extrp->pextract_n_entity              [PDM_MESH_ENTITY_FACE],
                                                          extrp->n_part_out,
                                                          (const PDM_g_num_t **) extrp->face_ln_to_gn,
                                                          extrp->n_face,
                                                          extrp->n_part_in,
                                                          (const int **) part2_face_to_part1_face_idx,
                                                          (const PDM_g_num_t **) extrp->pextract_entity_parent_ln_to_gn[PDM_MESH_ENTITY_FACE],
                                                          extrp->comm);
    extrp->ptp_entity[PDM_MESH_ENTITY_FACE] = ptp_fac;
    free(part2_face_to_part1_face_idx);

    for(int i_kind = 0; i_kind < PDM_BOUND_TYPE_MAX; ++i_kind) {
      _extract_part_group(extrp,
      (PDM_bound_type_t) i_kind);
    }
  }*/

}


static
void
_extract_part_from_target_rebuild_connectivities_3d
(
  PDM_extract_part_t        *extrp,
  int                       *pn_entity,
  int                      **entity_target_location,
  int                     ***pextract_face_to_face_location_out,
  int                     ***pextract_edge_to_edge_location_out,
  int                     ***pextract_vtx_to_vtx_location_out
)
{

  int have_cell_face_l = 0;
  int have_edge_vtx_l  = 0;
  int from_face_edge_l = 0;
  int from_face_vtx_l  = 0;
  for(int i_part = 0; i_part < extrp->n_part_in; ++i_part) {
    if(extrp->pface_edge_idx    [i_part] != NULL) {
      from_face_edge_l = 1;
    }
    if(extrp->pface_vtx_idx    [i_part] != NULL) {
      from_face_vtx_l = 1;
    }
    if(extrp->pcell_face_idx    [i_part] != NULL) {
      have_cell_face_l = 1;
    }
    if(extrp->pedge_vtx[i_part] != NULL) {
      have_edge_vtx_l = 1;
    }
  }
  // Procs can have n_part == 0 so we have to recover information
  int have_cell_face, have_edge_vtx, from_face_edge, from_face_vtx;
  PDM_MPI_Allreduce(&have_cell_face_l, &have_cell_face, 1, PDM_MPI_INT, PDM_MPI_MAX, extrp->comm);
  PDM_MPI_Allreduce(&have_edge_vtx_l,  &have_edge_vtx,  1, PDM_MPI_INT, PDM_MPI_MAX, extrp->comm);
  PDM_MPI_Allreduce(&from_face_edge_l, &from_face_edge, 1, PDM_MPI_INT, PDM_MPI_MAX, extrp->comm);
  PDM_MPI_Allreduce(&from_face_vtx_l,  &from_face_vtx,  1, PDM_MPI_INT, PDM_MPI_MAX, extrp->comm);


  int **part2_cell_to_part1_cell_idx = (int **) malloc( extrp->n_part_out * sizeof(int * ));
  int **part2_face_to_part1_face_idx = (int **) malloc( extrp->n_part_out * sizeof(int * ));
  int **part2_edge_to_part1_edge_idx = (int **) malloc( extrp->n_part_out * sizeof(int * ));
  int **part2_vtx_to_part1_vtx_idx   = (int **) malloc( extrp->n_part_out * sizeof(int * ));

  int **pextract_face_to_face_location = NULL;
  int **pextract_edge_to_edge_location = NULL;
  int **pextract_vtx_to_vtx_location   = NULL;

  for(int i_part = 0; i_part < extrp->n_part_out; ++i_part) {
    part2_cell_to_part1_cell_idx[i_part] = NULL;
    part2_face_to_part1_face_idx[i_part] = NULL;
    part2_edge_to_part1_edge_idx[i_part] = NULL;
    part2_vtx_to_part1_vtx_idx  [i_part] = NULL;
  }

  /* Create the link */
  for(int i_part = 0; i_part < extrp->n_part_out; ++i_part) {
    part2_cell_to_part1_cell_idx[i_part] = (int * ) malloc( (extrp->n_target[i_part]+1) * sizeof(int));
    part2_cell_to_part1_cell_idx[i_part][0] = 0;
    for(int i = 0; i < extrp->n_target[i_part]; ++i) {
      part2_cell_to_part1_cell_idx[i_part][i+1] = part2_cell_to_part1_cell_idx[i_part][i] + 3;
    }
  }

  if (extrp->pextract_n_entity[PDM_MESH_ENTITY_CELL] == NULL) {
    extrp->pextract_n_entity[PDM_MESH_ENTITY_CELL] = malloc(extrp->n_part_out * sizeof(int));
    for(int i_part = 0; i_part < extrp->n_part_out; ++i_part) {
      extrp->pextract_n_entity[PDM_MESH_ENTITY_CELL][i_part] = extrp->n_target[i_part];
    }
  }


  /*
   *  cell->face
   */
  if(have_cell_face == 1) {
    PDM_pconnectivity_to_pconnectivity_from_location_keep(extrp->comm,
                                                          extrp->n_part_in,
                                 (const int            *) pn_entity,
                                 (const int           **) extrp->pcell_face_idx,
                                 (const int           **) extrp->pcell_face,
                                 (const PDM_g_num_t   **) extrp->face_ln_to_gn, // TODO 2D
                                                          extrp->n_part_out,
                                 (const int            *) extrp->n_target,
                                 (const PDM_g_num_t   **) extrp->target_gnum,
                                 (const int           **) part2_cell_to_part1_cell_idx,
                                 (const int           **) entity_target_location,
                                                          &extrp->pextract_n_entity              [PDM_MESH_ENTITY_FACE],
                                                          &extrp->pextract_connectivity_idx      [PDM_CONNECTIVITY_TYPE_CELL_FACE],
                                                          &extrp->pextract_connectivity          [PDM_CONNECTIVITY_TYPE_CELL_FACE],
                                                          &extrp->pextract_entity_parent_ln_to_gn[PDM_MESH_ENTITY_FACE],
                                                          &pextract_face_to_face_location,
                                                          &extrp->ptp_entity[PDM_MESH_ENTITY_CELL]);

    if (extrp->compute_child_gnum) {
      extrp->pextract_entity_ln_to_gn[PDM_MESH_ENTITY_FACE] = malloc(sizeof(PDM_g_num_t *) * extrp->n_part_out);
      _compute_child(extrp->comm,
                     extrp->n_part_out,
                     extrp->pextract_n_entity              [PDM_MESH_ENTITY_FACE],
                     extrp->pextract_entity_parent_ln_to_gn[PDM_MESH_ENTITY_FACE],
                     extrp->pextract_entity_ln_to_gn       [PDM_MESH_ENTITY_FACE]);
    }

    /*
     * face -> edge
     */
    for(int i_part = 0; i_part < extrp->n_part_out; ++i_part) {
      int n_face = extrp->pextract_n_entity[PDM_MESH_ENTITY_FACE][i_part];
      part2_face_to_part1_face_idx[i_part] = PDM_array_new_idx_from_const_stride_int(3, n_face);
    }
  }

  if(from_face_edge == 1 && have_cell_face == 1) {
    PDM_pconnectivity_to_pconnectivity_from_location_keep(extrp->comm,
                                                          extrp->n_part_in,
                                 (const int            *) extrp->n_face,
                                 (const int           **) extrp->pface_edge_idx,
                                 (const int           **) extrp->pface_edge,
                                 (const PDM_g_num_t   **) extrp->edge_ln_to_gn, // TODO 2D
                                                          extrp->n_part_out,
                                 (const int            *) extrp->pextract_n_entity              [PDM_MESH_ENTITY_FACE],
                                 (const PDM_g_num_t   **) extrp->pextract_entity_parent_ln_to_gn[PDM_MESH_ENTITY_FACE],
                                 (const int           **) part2_face_to_part1_face_idx,
                                 (const int           **) pextract_face_to_face_location,
                                                          &extrp->pextract_n_entity              [PDM_MESH_ENTITY_EDGE],
                                                          &extrp->pextract_connectivity_idx      [PDM_CONNECTIVITY_TYPE_FACE_EDGE],
                                                          &extrp->pextract_connectivity          [PDM_CONNECTIVITY_TYPE_FACE_EDGE],
                                                          &extrp->pextract_entity_parent_ln_to_gn[PDM_MESH_ENTITY_EDGE],
                                                          &pextract_edge_to_edge_location,
                                                          &extrp->ptp_entity[PDM_MESH_ENTITY_FACE]);

    if (extrp->compute_child_gnum) {
      extrp->pextract_entity_ln_to_gn[PDM_MESH_ENTITY_EDGE] = malloc(sizeof(PDM_g_num_t *) * extrp->n_part_out);
      _compute_child(extrp->comm,
                     extrp->n_part_out,
                     extrp->pextract_n_entity              [PDM_MESH_ENTITY_EDGE],
                     extrp->pextract_entity_parent_ln_to_gn[PDM_MESH_ENTITY_EDGE],
                     extrp->pextract_entity_ln_to_gn       [PDM_MESH_ENTITY_EDGE]);
    }

    if(have_edge_vtx == 1) {
      int **pedge_vtx_idx = (int **) malloc(extrp->n_part_in * sizeof(int *));
      for(int i_part = 0; i_part < extrp->n_part_in; ++i_part) {
        pedge_vtx_idx[i_part] = (int *) malloc((extrp->n_edge[i_part]+1) * sizeof(int));
        pedge_vtx_idx[i_part][0] = 0;
        for(int i_edge = 0; i_edge < extrp->n_edge[i_part]; ++i_edge) {
          pedge_vtx_idx[i_part][i_edge+1] = pedge_vtx_idx[i_part][i_edge] + 2;
        }
      }

      for(int i_part = 0; i_part < extrp->n_part_out; ++i_part) {
        free(part2_face_to_part1_face_idx[i_part]);
        part2_face_to_part1_face_idx[i_part] = NULL;
      }

      for(int i_part = 0; i_part < extrp->n_part_out; ++i_part) {
        int n_edge = extrp->pextract_n_entity[PDM_MESH_ENTITY_EDGE][i_part];
        part2_edge_to_part1_edge_idx[i_part] =  PDM_array_new_idx_from_const_stride_int(3, n_edge);
      }

      PDM_pconnectivity_to_pconnectivity_from_location_keep(extrp->comm,
                                                            extrp->n_part_in,
                                   (const int            *) extrp->n_edge,
                                   (const int           **) pedge_vtx_idx,
                                   (const int           **) extrp->pedge_vtx,
                                   (const PDM_g_num_t   **) extrp->vtx_ln_to_gn, // TODO 2D
                                                            extrp->n_part_out,
                                   (const int            *) extrp->pextract_n_entity              [PDM_MESH_ENTITY_EDGE],
                                   (const PDM_g_num_t   **) extrp->pextract_entity_parent_ln_to_gn[PDM_MESH_ENTITY_EDGE],
                                   (const int           **) part2_edge_to_part1_edge_idx,
                                   (const int           **) pextract_edge_to_edge_location,
                                                            &extrp->pextract_n_entity              [PDM_MESH_ENTITY_VTX],
                                                            &extrp->pextract_connectivity_idx      [PDM_CONNECTIVITY_TYPE_EDGE_VTX],
                                                            &extrp->pextract_connectivity          [PDM_CONNECTIVITY_TYPE_EDGE_VTX],
                                                            &extrp->pextract_entity_parent_ln_to_gn[PDM_MESH_ENTITY_VTX],
                                                            &pextract_vtx_to_vtx_location,
                                                            &extrp->ptp_entity[PDM_MESH_ENTITY_EDGE]);

      if (extrp->compute_child_gnum) {
        extrp->pextract_entity_ln_to_gn[PDM_MESH_ENTITY_VTX] = malloc(sizeof(PDM_g_num_t *) * extrp->n_part_out);
        _compute_child(extrp->comm,
                       extrp->n_part_out,
                       extrp->pextract_n_entity              [PDM_MESH_ENTITY_VTX],
                       extrp->pextract_entity_parent_ln_to_gn[PDM_MESH_ENTITY_VTX],
                       extrp->pextract_entity_ln_to_gn       [PDM_MESH_ENTITY_VTX]);
      }

      for(int i_part = 0; i_part < extrp->n_part_in; ++i_part) {
        free(pedge_vtx_idx   [i_part]);
      }
      free(pedge_vtx_idx);

      for(int i_part = 0; i_part < extrp->n_part_out; ++i_part) {
        free(part2_edge_to_part1_edge_idx[i_part]);
        part2_edge_to_part1_edge_idx[i_part] = NULL;
      }
    }
  } else if(from_face_vtx == 1 && have_cell_face == 1) {
    PDM_pconnectivity_to_pconnectivity_from_location_keep(extrp->comm,
                                                          extrp->n_part_in,
                                 (const int            *) extrp->n_face,
                                 (const int           **) extrp->pface_vtx_idx,
                                 (const int           **) extrp->pface_vtx,
                                 (const PDM_g_num_t   **) extrp->vtx_ln_to_gn, // TODO 2D
                                                          extrp->n_part_out,
                                 (const int            *) extrp->pextract_n_entity              [PDM_MESH_ENTITY_FACE],
                                 (const PDM_g_num_t   **) extrp->pextract_entity_parent_ln_to_gn[PDM_MESH_ENTITY_FACE],
                                 (const int           **) part2_face_to_part1_face_idx,
                                 (const int           **) pextract_face_to_face_location,
                                                          &extrp->pextract_n_entity              [PDM_MESH_ENTITY_VTX],
                                                          &extrp->pextract_connectivity_idx      [PDM_CONNECTIVITY_TYPE_FACE_VTX],
                                                          &extrp->pextract_connectivity          [PDM_CONNECTIVITY_TYPE_FACE_VTX],
                                                          &extrp->pextract_entity_parent_ln_to_gn[PDM_MESH_ENTITY_VTX],
                                                          &pextract_vtx_to_vtx_location,
                                                          &extrp->ptp_entity[PDM_MESH_ENTITY_FACE]);

    if (extrp->compute_child_gnum) {
      extrp->pextract_entity_ln_to_gn[PDM_MESH_ENTITY_VTX] = malloc(sizeof(PDM_g_num_t *) * extrp->n_part_out);
      _compute_child(extrp->comm,
                     extrp->n_part_out,
                     extrp->pextract_n_entity              [PDM_MESH_ENTITY_VTX],
                     extrp->pextract_entity_parent_ln_to_gn[PDM_MESH_ENTITY_VTX],
                     extrp->pextract_entity_ln_to_gn       [PDM_MESH_ENTITY_VTX]);
    }

    for(int i_part = 0; i_part < extrp->n_part_out; ++i_part) {
      free(part2_face_to_part1_face_idx[i_part]);
      part2_face_to_part1_face_idx[i_part] = NULL;
    }
  }
  *pextract_face_to_face_location_out = pextract_face_to_face_location;
  *pextract_edge_to_edge_location_out = pextract_edge_to_edge_location;
  *pextract_vtx_to_vtx_location_out   = pextract_vtx_to_vtx_location;

  for(int i_part = 0; i_part < extrp->n_part_out; ++i_part) {
    free(part2_cell_to_part1_cell_idx[i_part]);
    part2_cell_to_part1_cell_idx[i_part] = NULL;
  }

  free(part2_cell_to_part1_cell_idx);
  free(part2_face_to_part1_face_idx);
  free(part2_edge_to_part1_edge_idx);
  free(part2_vtx_to_part1_vtx_idx  );
}

static
void
_extract_part_from_target_rebuild_connectivities_2d
(
  PDM_extract_part_t        *extrp,
  int                       *pn_entity,
  int                      **entity_target_location,
  int                     ***pextract_face_to_face_location_out,
  int                     ***pextract_edge_to_edge_location_out,
  int                     ***pextract_vtx_to_vtx_location_out
)
{
  int have_edge_vtx_l  = 0;
  int from_face_edge_l = 0;
  int from_face_vtx_l  = 0;
  for(int i_part = 0; i_part < extrp->n_part_in; ++i_part) {
    if(extrp->pface_edge_idx    [i_part] != NULL) {
      from_face_edge_l = 1;
    }
    if(extrp->pface_vtx_idx    [i_part] != NULL) {
      from_face_vtx_l = 1;
    }
    if(extrp->pedge_vtx[i_part] != NULL) {
      have_edge_vtx_l = 1;
    }
  }

  // Procs can have n_part == 0 so we have to recover information
  int have_edge_vtx, from_face_edge, from_face_vtx;
  PDM_MPI_Allreduce(&have_edge_vtx_l,  &have_edge_vtx,  1, PDM_MPI_INT, PDM_MPI_MAX, extrp->comm);
  PDM_MPI_Allreduce(&from_face_edge_l, &from_face_edge, 1, PDM_MPI_INT, PDM_MPI_MAX, extrp->comm);
  PDM_MPI_Allreduce(&from_face_vtx_l,  &from_face_vtx,  1, PDM_MPI_INT, PDM_MPI_MAX, extrp->comm);

  int **part2_face_to_part1_face_idx = (int **) malloc( extrp->n_part_out * sizeof(int * ));
  int **part2_edge_to_part1_edge_idx = (int **) malloc( extrp->n_part_out * sizeof(int * ));
  int **part2_vtx_to_part1_vtx_idx   = (int **) malloc( extrp->n_part_out * sizeof(int * ));

  int **pextract_face_to_face_location = NULL;
  int **pextract_edge_to_edge_location = NULL;
  int **pextract_vtx_to_vtx_location   = NULL;

  for(int i_part = 0; i_part < extrp->n_part_out; ++i_part) {
    part2_face_to_part1_face_idx[i_part] = NULL;
    part2_edge_to_part1_edge_idx[i_part] = NULL;
    part2_vtx_to_part1_vtx_idx  [i_part] = NULL;
  }

  /* Create the link */
  for(int i_part = 0; i_part < extrp->n_part_out; ++i_part) {
    part2_face_to_part1_face_idx[i_part] = (int * ) malloc( (extrp->n_target[i_part]+1) * sizeof(int));
    part2_face_to_part1_face_idx[i_part][0] = 0;
    for(int i = 0; i < extrp->n_target[i_part]; ++i) {
      part2_face_to_part1_face_idx[i_part][i+1] = part2_face_to_part1_face_idx[i_part][i] + 3;
    }
  }

  if (extrp->pextract_n_entity[PDM_MESH_ENTITY_FACE] == NULL) {
    extrp->pextract_n_entity[PDM_MESH_ENTITY_FACE] = malloc(extrp->n_part_out * sizeof(int));
    for(int i_part = 0; i_part < extrp->n_part_out; ++i_part) {
      extrp->pextract_n_entity[PDM_MESH_ENTITY_FACE][i_part] = extrp->n_target[i_part];
    }
  }

  if(from_face_edge == 1) {
    PDM_pconnectivity_to_pconnectivity_from_location_keep(extrp->comm,
                                                          extrp->n_part_in,
                                 (const int            *) pn_entity,
                                 (const int           **) extrp->pface_edge_idx,
                                 (const int           **) extrp->pface_edge,
                                 (const PDM_g_num_t   **) extrp->edge_ln_to_gn, // TODO 2D
                                                          extrp->n_part_out,
                                 (const int            *) extrp->n_target,
                                 (const PDM_g_num_t   **) extrp->target_gnum,
                                 (const int           **) part2_face_to_part1_face_idx,
                                 (const int           **) entity_target_location,
                                                          &extrp->pextract_n_entity              [PDM_MESH_ENTITY_EDGE],
                                                          &extrp->pextract_connectivity_idx      [PDM_CONNECTIVITY_TYPE_FACE_EDGE],
                                                          &extrp->pextract_connectivity          [PDM_CONNECTIVITY_TYPE_FACE_EDGE],
                                                          &extrp->pextract_entity_parent_ln_to_gn[PDM_MESH_ENTITY_EDGE],
                                                          &pextract_edge_to_edge_location,
                                                          &extrp->ptp_entity[PDM_MESH_ENTITY_FACE]);

    if (extrp->compute_child_gnum) {
      extrp->pextract_entity_ln_to_gn[PDM_MESH_ENTITY_EDGE] = malloc(sizeof(PDM_g_num_t *) * extrp->n_part_out);
      _compute_child(extrp->comm,
                     extrp->n_part_out,
                     extrp->pextract_n_entity              [PDM_MESH_ENTITY_EDGE],
                     extrp->pextract_entity_parent_ln_to_gn[PDM_MESH_ENTITY_EDGE],
                     extrp->pextract_entity_ln_to_gn       [PDM_MESH_ENTITY_EDGE]);
    }


    if(have_edge_vtx == 1) {

      int **pedge_vtx_idx = (int **) malloc(extrp->n_part_in * sizeof(int *));
      for(int i_part = 0; i_part < extrp->n_part_in; ++i_part) {
        pedge_vtx_idx[i_part] = (int *) malloc((extrp->n_edge[i_part]+1) * sizeof(int));
        pedge_vtx_idx[i_part][0] = 0;
        for(int i_edge = 0; i_edge < extrp->n_edge[i_part]; ++i_edge) {
          pedge_vtx_idx[i_part][i_edge+1] = pedge_vtx_idx[i_part][i_edge] + 2;
        }
      }

      for(int i_part = 0; i_part < extrp->n_part_out; ++i_part) {
        int n_edge = extrp->pextract_n_entity[PDM_MESH_ENTITY_EDGE][i_part];
        part2_edge_to_part1_edge_idx[i_part] =  PDM_array_new_idx_from_const_stride_int(3, n_edge);;
      }

      PDM_pconnectivity_to_pconnectivity_from_location_keep(extrp->comm,
                                                            extrp->n_part_in,
                                   (const int            *) extrp->n_edge,
                                   (const int           **) pedge_vtx_idx,
                                   (const int           **) extrp->pedge_vtx,
                                   (const PDM_g_num_t   **) extrp->vtx_ln_to_gn, // TODO 2D
                                                            extrp->n_part_out,
                                   (const int            *) extrp->pextract_n_entity              [PDM_MESH_ENTITY_EDGE],
                                   (const PDM_g_num_t   **) extrp->pextract_entity_parent_ln_to_gn[PDM_MESH_ENTITY_EDGE],
                                   (const int           **) part2_edge_to_part1_edge_idx,
                                   (const int           **) pextract_edge_to_edge_location,
                                                            &extrp->pextract_n_entity              [PDM_MESH_ENTITY_VTX],
                                                            &extrp->pextract_connectivity_idx      [PDM_CONNECTIVITY_TYPE_EDGE_VTX],
                                                            &extrp->pextract_connectivity          [PDM_CONNECTIVITY_TYPE_EDGE_VTX],
                                                            &extrp->pextract_entity_parent_ln_to_gn[PDM_MESH_ENTITY_VTX],
                                                            &pextract_vtx_to_vtx_location,
                                                            &extrp->ptp_entity[PDM_MESH_ENTITY_EDGE]);

      if (extrp->compute_child_gnum) {
        extrp->pextract_entity_ln_to_gn[PDM_MESH_ENTITY_VTX] = malloc(sizeof(PDM_g_num_t *) * extrp->n_part_out);
        _compute_child(extrp->comm,
                       extrp->n_part_out,
                       extrp->pextract_n_entity              [PDM_MESH_ENTITY_VTX],
                       extrp->pextract_entity_parent_ln_to_gn[PDM_MESH_ENTITY_VTX],
                       extrp->pextract_entity_ln_to_gn       [PDM_MESH_ENTITY_VTX]);
      }

      for(int i_part = 0; i_part < extrp->n_part_out; ++i_part) {
        free(part2_edge_to_part1_edge_idx[i_part]);
        part2_edge_to_part1_edge_idx[i_part] = NULL;
      }

      for(int i_part = 0; i_part < extrp->n_part_in; ++i_part) {
        free(pedge_vtx_idx   [i_part]);
      }
      free(pedge_vtx_idx);
    }

  } else if(from_face_vtx == 1) { // from_face_vtx
    PDM_pconnectivity_to_pconnectivity_from_location_keep(extrp->comm,
                                                          extrp->n_part_in,
                                 (const int            *) pn_entity,
                                 (const int           **) extrp->pface_vtx_idx,
                                 (const int           **) extrp->pface_vtx,
                                 (const PDM_g_num_t   **) extrp->vtx_ln_to_gn, // TODO 2D
                                                          extrp->n_part_out,
                                 (const int            *) extrp->n_target,
                                 (const PDM_g_num_t   **) extrp->target_gnum,
                                 (const int           **) part2_face_to_part1_face_idx,
                                 (const int           **) entity_target_location,
                                                          &extrp->pextract_n_entity              [PDM_MESH_ENTITY_VTX],
                                                          &extrp->pextract_connectivity_idx      [PDM_CONNECTIVITY_TYPE_FACE_VTX],
                                                          &extrp->pextract_connectivity          [PDM_CONNECTIVITY_TYPE_FACE_VTX],
                                                          &extrp->pextract_entity_parent_ln_to_gn[PDM_MESH_ENTITY_VTX],
                                                          &pextract_vtx_to_vtx_location,
                                                          &extrp->ptp_entity[PDM_MESH_ENTITY_FACE]);

    if (extrp->compute_child_gnum) {
      extrp->pextract_entity_ln_to_gn[PDM_MESH_ENTITY_VTX] = malloc(sizeof(PDM_g_num_t *) * extrp->n_part_out);
      _compute_child(extrp->comm,
                     extrp->n_part_out,
                     extrp->pextract_n_entity              [PDM_MESH_ENTITY_VTX],
                     extrp->pextract_entity_parent_ln_to_gn[PDM_MESH_ENTITY_VTX],
                     extrp->pextract_entity_ln_to_gn       [PDM_MESH_ENTITY_VTX]);
    }

    for(int i_part = 0; i_part < extrp->n_part_out; ++i_part) {
      free(part2_face_to_part1_face_idx[i_part]);
      part2_face_to_part1_face_idx[i_part] = NULL;
    }
  }

  for(int i_part = 0; i_part < extrp->n_part_out; ++i_part) {
    free(part2_face_to_part1_face_idx[i_part]);
    part2_face_to_part1_face_idx[i_part] = NULL;
  }

  *pextract_face_to_face_location_out = pextract_face_to_face_location;
  *pextract_edge_to_edge_location_out = pextract_edge_to_edge_location;
  *pextract_vtx_to_vtx_location_out   = pextract_vtx_to_vtx_location;

  free(part2_face_to_part1_face_idx);
  free(part2_edge_to_part1_edge_idx);
  free(part2_vtx_to_part1_vtx_idx  );

}


static
void
_extract_part_and_reequilibrate_from_target
(
  PDM_extract_part_t        *extrp
)
{
  int                *pn_entity              = NULL;
  PDM_g_num_t       **entity_g_num           = NULL;
  int               **entity_target_location = NULL;
  PDM_mesh_entities_t entity_type;
  if(extrp->dim == 3) {
    pn_entity    = extrp->n_cell;
    entity_g_num = extrp->cell_ln_to_gn;
    entity_type = PDM_MESH_ENTITY_CELL;
  } else if (extrp->dim == 2){
    pn_entity    = extrp->n_face;
    entity_g_num = extrp->face_ln_to_gn;
    entity_type = PDM_MESH_ENTITY_FACE;
  } else if(extrp->dim == 1){
    pn_entity    = extrp->n_edge;
    entity_g_num = extrp->edge_ln_to_gn;
    entity_type = PDM_MESH_ENTITY_EDGE;
  } else {
    entity_g_num = extrp->vtx_ln_to_gn;
    entity_type = PDM_MESH_ENTITY_VTX;
  }

  /* Not very bright... let the user be in charge of keeping track of the target g_num instead? */
  // -->>
  // (Copy to avoid double free)
//  if (1) {
  if (extrp->target_ownership == PDM_OWNERSHIP_KEEP) {
    extrp->pextract_entity_parent_ln_to_gn[entity_type] = malloc(sizeof(PDM_g_num_t * ) * extrp->n_part_out);
    for (int ipart = 0; ipart < extrp->n_part_out; ipart++) {
      // log_trace("extrp->n_target[%d] = %d\n", ipart, extrp->n_target[ipart]);
      extrp->pextract_entity_parent_ln_to_gn[entity_type][ipart] = extrp->target_gnum[ipart];
      // extrp->pextract_entity_parent_ln_to_gn[entity_type][ipart] = malloc(sizeof(PDM_g_num_t) * extrp->n_target[ipart]);
      // memcpy(extrp->pextract_entity_parent_ln_to_gn[entity_type][ipart],
      //        extrp->target_gnum[ipart],
      //        sizeof(PDM_g_num_t) * extrp->n_target[ipart]);
    }
  }
  // <<--

  int have_init_location_l = 1;
  int have_init_location   = 1;
  for(int i_part = 0; i_part < extrp->n_part_out; ++i_part) {
    if(extrp->target_location[i_part] == NULL) {
      have_init_location = 0;
    }
  }
  PDM_MPI_Allreduce(&have_init_location_l, &have_init_location, 1, PDM_MPI_INT, PDM_MPI_MAX, extrp->comm);
  // log_trace("have_init_location = %i \n", have_init_location);

  if(have_init_location == 0) {
    PDM_gnum_location_t* gnum_loc = PDM_gnum_location_create(extrp->n_part_in,
                                                             extrp->n_part_out,
                                                             extrp->comm,
                                                             PDM_OWNERSHIP_USER);
    for(int i_part = 0; i_part < extrp->n_part_in; ++i_part) {
      PDM_gnum_location_elements_set(gnum_loc,
                                     i_part,
                                     pn_entity[i_part],
                                     entity_g_num[i_part]);
    }
    for(int i_part = 0; i_part < extrp->n_part_out; ++i_part) {
      PDM_gnum_location_requested_elements_set(gnum_loc,
                                               i_part,
                                               extrp->n_target   [i_part],
                                               extrp->target_gnum[i_part]);
    }
    PDM_gnum_location_compute(gnum_loc);

    entity_target_location = malloc(extrp->n_part_out * sizeof(int *));

    for(int i_part = 0; i_part < extrp->n_part_out; ++i_part) {
      int *location_idx = NULL;
      PDM_gnum_location_get(gnum_loc, i_part, &location_idx, &entity_target_location[i_part]);
      free(location_idx);
    }

    PDM_gnum_location_free(gnum_loc);
  } else {
    entity_target_location = extrp->target_location;
  }


  int **part2_cell_to_part1_cell_idx = (int **) malloc( extrp->n_part_out * sizeof(int * ));
  int **part2_face_to_part1_face_idx = (int **) malloc( extrp->n_part_out * sizeof(int * ));
  int **part2_edge_to_part1_edge_idx = (int **) malloc( extrp->n_part_out * sizeof(int * ));
  int **part2_vtx_to_part1_vtx_idx   = (int **) malloc( extrp->n_part_out * sizeof(int * ));

  for(int i_part = 0; i_part < extrp->n_part_out; ++i_part) {
    part2_cell_to_part1_cell_idx[i_part] = NULL;
    part2_face_to_part1_face_idx[i_part] = NULL;
    part2_edge_to_part1_edge_idx[i_part] = NULL;
    part2_vtx_to_part1_vtx_idx  [i_part] = NULL;
  }

  /*
   * Extraction des connectivités
   */
  PDM_part_to_part_t* ptp_vtx          = NULL;
  int **pextract_face_to_face_location = NULL;
  int **pextract_edge_to_edge_location = NULL;
  int **pextract_vtx_to_vtx_location   = NULL;

  if(extrp->dim == 3) {
    _extract_part_from_target_rebuild_connectivities_3d(extrp,
                                                        pn_entity,
                                                        entity_target_location,
                                                        &pextract_face_to_face_location,
                                                        &pextract_edge_to_edge_location,
                                                        &pextract_vtx_to_vtx_location);
  } else if(extrp->dim == 2){
    _extract_part_from_target_rebuild_connectivities_2d(extrp,
                                                        pn_entity,
                                                        entity_target_location,
                                                        &pextract_face_to_face_location,
                                                        &pextract_edge_to_edge_location,
                                                        &pextract_vtx_to_vtx_location);

  } else if(extrp->dim == 1){

    /* Create the link */
    for(int i_part = 0; i_part < extrp->n_part_out; ++i_part) {
      part2_edge_to_part1_edge_idx[i_part] = (int * ) malloc( (extrp->n_target[i_part]+1) * sizeof(int));
      part2_edge_to_part1_edge_idx[i_part][0] = 0;
      for(int i = 0; i < extrp->n_target[i_part]; ++i) {
        part2_edge_to_part1_edge_idx[i_part][i+1] = part2_edge_to_part1_edge_idx[i_part][i] + 3;
      }
    }

    int **pedge_vtx_idx = (int **) malloc(extrp->n_part_in * sizeof(int *));
    for(int i_part = 0; i_part < extrp->n_part_in; ++i_part) {
      pedge_vtx_idx[i_part] = (int *) malloc((extrp->n_edge[i_part]+1) * sizeof(int));
      pedge_vtx_idx[i_part][0] = 0;
      for(int i_edge = 0; i_edge < extrp->n_edge[i_part]; ++i_edge) {
        pedge_vtx_idx[i_part][i_edge+1] = pedge_vtx_idx[i_part][i_edge] + 2;
      }
    }

    if (extrp->pextract_n_entity[PDM_MESH_ENTITY_EDGE] == NULL) {
      extrp->pextract_n_entity[PDM_MESH_ENTITY_EDGE] = malloc(extrp->n_part_out * sizeof(int));
      for(int i_part = 0; i_part < extrp->n_part_out; ++i_part) {
        extrp->pextract_n_entity[PDM_MESH_ENTITY_EDGE][i_part] = extrp->n_target[i_part];
      }
    }

    PDM_pconnectivity_to_pconnectivity_from_location_keep(extrp->comm,
                                                          extrp->n_part_in,
                                 (const int            *) extrp->n_edge,
                                 (const int           **) pedge_vtx_idx,
                                 (const int           **) extrp->pedge_vtx,
                                 (const PDM_g_num_t   **) extrp->vtx_ln_to_gn,
                                                          extrp->n_part_out,
                                 (const int            *) extrp->n_target,
                                 (const PDM_g_num_t   **) extrp->target_gnum,
                                 (const int           **) part2_edge_to_part1_edge_idx,
                                 (const int           **) entity_target_location,
                                                          &extrp->pextract_n_entity              [PDM_MESH_ENTITY_VTX],
                                                          &extrp->pextract_connectivity_idx      [PDM_CONNECTIVITY_TYPE_EDGE_VTX],
                                                          &extrp->pextract_connectivity          [PDM_CONNECTIVITY_TYPE_EDGE_VTX],
                                                          &extrp->pextract_entity_parent_ln_to_gn[PDM_MESH_ENTITY_VTX],
                                                          &pextract_vtx_to_vtx_location,
                                                          &extrp->ptp_entity[PDM_MESH_ENTITY_EDGE]);

    if (extrp->compute_child_gnum) {
      extrp->pextract_entity_ln_to_gn[PDM_MESH_ENTITY_VTX] = malloc(sizeof(PDM_g_num_t *) * extrp->n_part_out);
      _compute_child(extrp->comm,
                     extrp->n_part_out,
                     extrp->pextract_n_entity              [PDM_MESH_ENTITY_VTX],
                     extrp->pextract_entity_parent_ln_to_gn[PDM_MESH_ENTITY_VTX],
                     extrp->pextract_entity_ln_to_gn       [PDM_MESH_ENTITY_VTX]);
    }

    for(int i_part = 0; i_part < extrp->n_part_out; ++i_part) {
      free(part2_edge_to_part1_edge_idx[i_part]);
      part2_edge_to_part1_edge_idx[i_part] = NULL;
    }

    for(int i_part = 0; i_part < extrp->n_part_in; ++i_part) {
      free(pedge_vtx_idx   [i_part]);
    }
    free(pedge_vtx_idx);

  } else { // dim == 0
    assert(ptp_vtx == NULL);

    /* Create the link */
    for(int i_part = 0; i_part < extrp->n_part_out; ++i_part) {
      part2_vtx_to_part1_vtx_idx[i_part] = (int * ) malloc( (extrp->n_target[i_part]+1) * sizeof(int));
      part2_vtx_to_part1_vtx_idx[i_part][0] = 0;
      for(int i = 0; i < extrp->n_target[i_part]; ++i) {
        part2_vtx_to_part1_vtx_idx[i_part][i+1] = part2_vtx_to_part1_vtx_idx[i_part][i] + 3;
      }
    }

    if (extrp->pextract_n_entity[PDM_MESH_ENTITY_VTX] == NULL) {
      extrp->pextract_n_entity[PDM_MESH_ENTITY_VTX] = malloc(extrp->n_part_out * sizeof(int));
      for(int i_part = 0; i_part < extrp->n_part_out; ++i_part) {
        extrp->pextract_n_entity[PDM_MESH_ENTITY_VTX][i_part] = extrp->n_target[i_part];
      }
    }

    ptp_vtx = PDM_part_to_part_create_from_num2_triplet((const PDM_g_num_t **) extrp->target_gnum,
                                                        extrp->n_target,
                                                        extrp->n_part_out,
                                                        extrp->n_vtx,
                                                        extrp->n_part_in,
                                                        (const int **) part2_vtx_to_part1_vtx_idx,
                                                        NULL,
                                                        (const int **) entity_target_location,
                                                        extrp->comm);
  }


  if(have_init_location == 0) {
    for(int i_part = 0; i_part < extrp->n_part_out; ++i_part) {
      free(entity_target_location[i_part]);
    }
    free(entity_target_location);
  }

  /*
   * Vtx only
   */
  if(extrp->dim != 0 && pextract_vtx_to_vtx_location != NULL) {
    for(int i_part = 0; i_part < extrp->n_part_out; ++i_part) {
      int n_vtx = extrp->pextract_n_entity[PDM_MESH_ENTITY_VTX][i_part];
      part2_vtx_to_part1_vtx_idx[i_part] =  PDM_array_new_idx_from_const_stride_int(3, n_vtx);;
    }

    ptp_vtx = PDM_part_to_part_create_from_num2_triplet((const PDM_g_num_t **) extrp->pextract_entity_parent_ln_to_gn[PDM_MESH_ENTITY_VTX],
                                                        (const int          *) extrp->pextract_n_entity              [PDM_MESH_ENTITY_VTX],
                                                        extrp->n_part_out,
                                                        extrp->n_vtx,
                                                        extrp->n_part_in,
                                                        (const int **) part2_vtx_to_part1_vtx_idx,
                                                        NULL,
                                                        (const int **) pextract_vtx_to_vtx_location,
                                                        extrp->comm);
  }
  extrp->ptp_entity[PDM_MESH_ENTITY_VTX] = ptp_vtx;

  if(ptp_vtx != NULL) {
    int           exch_request = -1;
    PDM_part_to_part_reverse_iexch(ptp_vtx,
                                   PDM_MPI_COMM_KIND_P2P,
                                   PDM_STRIDE_CST_INTERLACED,
                                   PDM_PART_TO_PART_DATA_DEF_ORDER_PART2,
                                   3,
                                   sizeof(double),
                                   NULL,
                  (const void **)  extrp->pvtx_coord,
                                   NULL,
                      (void ***)   &extrp->pextract_vtx_coord,
                                   &exch_request);
    PDM_part_to_part_reverse_iexch_wait(ptp_vtx, exch_request);
  }

  for(int i_part = 0; i_part < extrp->n_part_out; ++i_part) {
    free(part2_vtx_to_part1_vtx_idx[i_part]);
  }

  /*
   *
   */
  if(extrp->dim != 0 && pextract_vtx_to_vtx_location != NULL) {
    for(int i_part = 0; i_part < extrp->n_part_out; ++i_part) {
      free(pextract_vtx_to_vtx_location[i_part]);
    }
    free(pextract_vtx_to_vtx_location);
  }

  if(pextract_face_to_face_location != NULL) {
    for(int i_part = 0; i_part < extrp->n_part_out; ++i_part) {
      if(pextract_face_to_face_location[i_part] != NULL) {
        free(pextract_face_to_face_location[i_part]);
      }
    }
    free(pextract_face_to_face_location);
  }

  if(pextract_edge_to_edge_location != NULL) {
    for(int i_part = 0; i_part < extrp->n_part_out; ++i_part) {
      if(pextract_edge_to_edge_location[i_part] != NULL) {
        free(pextract_edge_to_edge_location[i_part]);
      }
    }
    free(pextract_edge_to_edge_location);
  }

  free(part2_cell_to_part1_cell_idx);
  free(part2_face_to_part1_face_idx);
  free(part2_edge_to_part1_edge_idx);
  free(part2_vtx_to_part1_vtx_idx  );

  /*
   * Manage group if any
   */
  for(int i_kind = 0; i_kind < PDM_BOUND_TYPE_MAX; ++i_kind) {
    _extract_part_group(extrp,
     (PDM_bound_type_t) i_kind);
  }
}

static
void
_extract_part_and_reequilibrate
(
  PDM_extract_part_t        *extrp
)
{
  int i_rank;
  PDM_MPI_Comm_rank(extrp->comm, &i_rank);

  PDM_g_num_t **entity_g_num = NULL;
  PDM_mesh_entities_t entity_type = PDM_MESH_ENTITY_CELL;
  if(extrp->dim == 3) {
    entity_g_num = extrp->cell_ln_to_gn;
    entity_type  = PDM_MESH_ENTITY_CELL;
  } else if(extrp->dim == 2) {
    entity_g_num = extrp->face_ln_to_gn;
    entity_type  = PDM_MESH_ENTITY_FACE;
  } else if(extrp->dim == 1) {
    entity_g_num = extrp->edge_ln_to_gn;
    entity_type  = PDM_MESH_ENTITY_EDGE;
  } else if(extrp->dim == 0) {
    entity_g_num = extrp->vtx_ln_to_gn;
    entity_type  = PDM_MESH_ENTITY_VTX;
  } else {
    abort();
  }

  int from_face_edge = 0;
  for(int i_part = 0; i_part < extrp->n_part_in; ++i_part) {
    if(extrp->pface_edge[i_part] != NULL) {
      from_face_edge = 1;
    }
  }

  /*
   * Calcul des coordonnées to setup hilbert ordering (independant of parallelism )
   */
  double **entity_center = NULL;
  if(extrp->have_user_entity_center == 1) {
    entity_center = malloc(extrp->n_part_in * sizeof(double * ));
    for (int i_part = 0; i_part < extrp->n_part_in; i_part++) {
      entity_center[i_part] = malloc(3 * extrp->n_extract[i_part] * sizeof(double));
      for(int i = 0; i < extrp->n_extract[i_part]; ++i) {
        int lnum = extrp->extract_lnum[i_part][i];
        entity_center[i_part][3*i  ] = extrp->entity_center[i_part][3*lnum  ];
        entity_center[i_part][3*i+1] = extrp->entity_center[i_part][3*lnum+1];
        entity_center[i_part][3*i+2] = extrp->entity_center[i_part][3*lnum+2];
      }
    }
  } else {
    if(extrp->dim == 3) {
      _cell_center_3d(extrp->n_part_in,
                      extrp->n_extract,
                      extrp->extract_lnum,
                      extrp->pcell_face_idx,
                      extrp->pcell_face,
                      extrp->pface_edge_idx,
                      extrp->pface_edge,
                      extrp->pface_vtx_idx,
                      extrp->pface_vtx,
                      extrp->pedge_vtx,
                      extrp->pvtx_coord,
                      &entity_center);
    } else if(extrp->dim == 2) {
      // Compute entity_center with face_edge + edge_vtx
      if (from_face_edge) {
        _face_center_2d(extrp->n_part_in,
                        extrp->n_extract,
                        extrp->extract_lnum,
                        extrp->pface_edge_idx,
                        extrp->pface_edge,
                        extrp->pedge_vtx,
                        extrp->pvtx_coord,
                        &entity_center);
      }
      else {
       _face_center_2d_from_vtx(extrp->n_part_in,
                                extrp->n_extract,
                                extrp->extract_lnum,
                                extrp->pface_vtx_idx,
                                extrp->pface_vtx,
                                extrp->pvtx_coord,
                                &entity_center);
      }
    } else if (extrp->dim == 1) {
      _edge_center_2d(extrp->n_part_in,
                      extrp->n_extract,
                      extrp->extract_lnum,
                      extrp->pedge_vtx,
                      extrp->pvtx_coord,
                      &entity_center);

    } else {
      entity_center = malloc(extrp->n_part_in * sizeof(double * ));
      for (int i_part = 0; i_part < extrp->n_part_in; i_part++) {
        entity_center[i_part] = malloc(3 * extrp->n_extract[i_part] * sizeof(double));

        for(int i = 0; i < extrp->n_extract[i_part]; ++i) {
          int lnum = extrp->extract_lnum[i_part][i];
          entity_center[i_part][3*i  ] = extrp->pvtx_coord[i_part][3*lnum  ];
          entity_center[i_part][3*i+1] = extrp->pvtx_coord[i_part][3*lnum+1];
          entity_center[i_part][3*i+2] = extrp->pvtx_coord[i_part][3*lnum+2];
        }
      }
    }
  }

  double      **weight              = malloc(sizeof(double      *) * extrp->n_part_in);
  PDM_g_num_t **extract_entity_gnum = malloc(sizeof(PDM_g_num_t *) * extrp->n_part_in);
  for (int i_part = 0; i_part < extrp->n_part_in; i_part++) {
    weight             [i_part] = PDM_array_const_double(extrp->n_extract[i_part], 1.);
    extract_entity_gnum[i_part] = malloc(extrp->n_extract[i_part] * sizeof(PDM_g_num_t));

    for(int i = 0; i < extrp->n_extract[i_part]; ++i) {
      int lnum = extrp->extract_lnum[i_part][i];
      extract_entity_gnum[i_part][i] = entity_g_num[i_part][lnum];
    }
  }

  /*
   *  Remake equilibrate block -> Block is not partial because we use child_gnum
   */
  PDM_part_to_block_t *ptb_equi = PDM_part_to_block_geom_create(PDM_PART_TO_BLOCK_DISTRIB_ALL_PROC,
                                                                PDM_PART_TO_BLOCK_POST_CLEANUP,
                                                                1.,
                                                                PDM_PART_GEOM_HILBERT,
                                                                entity_center,
                                                                extract_entity_gnum,
                                                                weight,
                                                                extrp->n_extract,
                                                                extrp->n_part_in,
                                                                extrp->comm);


  for (int i_part = 0; i_part < extrp->n_part_in; i_part++) {
    free(weight       [i_part]);
    free(entity_center[i_part]);
    free(extract_entity_gnum[i_part]);
  }
  free(weight       );
  free(entity_center);
  free(extract_entity_gnum);

  int dn_equi = PDM_part_to_block_n_elt_block_get(ptb_equi);
  PDM_g_num_t *dequi_g_num = PDM_part_to_block_block_gnum_get(ptb_equi);
  if (0 == 1) {
    PDM_log_trace_array_long(dequi_g_num, dn_equi, "dequi_g_num        : ");
  }

  /*
   * Deduce a target gnum from graph librarie or directly with hilbert
   */
  if(extrp->split_dual_method == PDM_SPLIT_DUAL_WITH_HILBERT) {
    assert(extrp->n_part_out == 1);

    PDM_g_num_t* _parent_entity_ln_to_gn = malloc(dn_equi * sizeof(PDM_g_num_t));
    for(int i = 0; i < dn_equi; ++i) {
      _parent_entity_ln_to_gn[i] = dequi_g_num[i];
    }

    int i_part0 = 0;
    extrp->pextract_entity_parent_ln_to_gn[entity_type] = malloc(sizeof(PDM_g_num_t *) * extrp->n_part_out);
    extrp->pextract_entity_parent_ln_to_gn[entity_type][i_part0] = _parent_entity_ln_to_gn;

    int **extract_init_location = malloc(sizeof(int *) * extrp->n_part_in);
    for (int i_part = 0; i_part < extrp->n_part_in; i_part++) {
      extract_init_location[i_part] = malloc(3 * extrp->n_extract[i_part] * sizeof(int));
      for(int i = 0; i < extrp->n_extract[i_part]; ++i) {
        extract_init_location[i_part][3*i  ] = i_rank;
        extract_init_location[i_part][3*i+1] = i_part;
        extract_init_location[i_part][3*i+2] = extrp->extract_lnum[i_part][i];
      }
    }

    int request_entity_init_location = -1;
    int *dequi_init_location = NULL;
    PDM_part_to_block_iexch(ptb_equi,
                            PDM_MPI_COMM_KIND_COLLECTIVE,
                            3 * sizeof(int),
                            PDM_STRIDE_CST_INTERLACED,
                            1,
                            NULL,
                  (void **) extract_init_location,
                            NULL,
                  (void **) &dequi_init_location,
                            &request_entity_init_location);

    PDM_g_num_t* cell_distri = PDM_part_to_block_distrib_index_get(ptb_equi);
    int dn_cell_equi = cell_distri[i_rank+1] - cell_distri[i_rank];

    extrp->pextract_n_entity       [entity_type] = (int          *) malloc(extrp->n_part_out * sizeof(int          ));
    extrp->pextract_entity_ln_to_gn[entity_type] = (PDM_g_num_t **) malloc(extrp->n_part_out * sizeof(PDM_g_num_t *));

    extrp->pextract_n_entity       [entity_type][i_part0] = dn_cell_equi;

    extrp->pextract_entity_ln_to_gn[entity_type][i_part0] = NULL;
    if(extrp->compute_child_gnum == 1) {
      extrp->pextract_entity_ln_to_gn[entity_type][i_part0] = (PDM_g_num_t *) malloc( extrp->pextract_n_entity[entity_type][i_part0] * sizeof(PDM_g_num_t));
      for(int i = 0; i < dn_cell_equi; ++i) {
        extrp->pextract_entity_ln_to_gn[entity_type][i_part0][i] = cell_distri[i_rank] + i + 1;
      }
    }

    // Early return
    PDM_part_to_block_iexch_wait(ptb_equi, request_entity_init_location);
    for (int i_part = 0; i_part < extrp->n_part_in; i_part++) {
      free(extract_init_location[i_part]);
    }
    free(extract_init_location);

    // Call extraction from target !
    extrp->n_target       [i_part0] = dn_cell_equi;
    // extrp->target_gnum    [i_part0] = extrp->pextract_entity_ln_to_gn[entity_type][i_part0];
    extrp->target_gnum    [i_part0] = _parent_entity_ln_to_gn;
    extrp->target_location[i_part0] = dequi_init_location;
    _extract_part_and_reequilibrate_from_target(extrp);

    extrp->n_target       [i_part0] = 0;
    extrp->target_gnum    [i_part0] = NULL;
    extrp->target_location[i_part0] = NULL;
    free(dequi_init_location);

    PDM_part_to_block_free(ptb_equi);

    return;
  }


  PDM_g_num_t *distrib_elmt     = PDM_part_to_block_distrib_index_get(ptb_equi);

  PDM_g_num_t *dual_graph_idx = NULL;
  PDM_g_num_t *dual_graph     = NULL;


  int **extract_init_location = malloc(sizeof(int *) * extrp->n_part_in);
  for (int i_part = 0; i_part < extrp->n_part_in; i_part++) {
    extract_init_location[i_part] = malloc(3 * extrp->n_extract[i_part] * sizeof(int));
    for(int i = 0; i < extrp->n_extract[i_part]; ++i) {
      extract_init_location[i_part][3*i  ] = i_rank;
      extract_init_location[i_part][3*i+1] = i_part;
      extract_init_location[i_part][3*i+2] = extrp->extract_lnum[i_part][i];
    }
  }

  int request_entity_init_location = -1;
  int *dequi_init_location = NULL;
  PDM_part_to_block_iexch(ptb_equi,
                          PDM_MPI_COMM_KIND_COLLECTIVE,
                          3 * sizeof(int),
                          PDM_STRIDE_CST_INTERLACED,
                          1,
                          NULL,
                (void **) extract_init_location,
                          NULL,
                (void **) &dequi_init_location,
                          &request_entity_init_location);
  PDM_part_to_block_iexch_wait(ptb_equi, request_entity_init_location);
  for (int i_part = 0; i_part < extrp->n_part_in; i_part++) {
    free(extract_init_location[i_part]);
  }
  free(extract_init_location);


  _compute_dual_graph(extrp, ptb_equi, &dual_graph_idx, &dual_graph);

  int tn_part;
  PDM_MPI_Allreduce(&extrp->n_part_out, &tn_part, 1, PDM_MPI_INT, PDM_MPI_SUM, extrp->comm);

  int dn_elmt = distrib_elmt[i_rank+1] - distrib_elmt[i_rank];
  int *_elmt_part = malloc(dn_elmt * sizeof(int));
  PDM_para_graph_split(extrp->split_dual_method,
                       distrib_elmt,
                       dual_graph_idx,
                       dual_graph,
                       NULL,
                       NULL,
                       tn_part,
                       NULL,
                       _elmt_part,
                       extrp->comm);

  PDM_g_num_t *distrib_partition = PDM_compute_entity_distribution(extrp->comm, extrp->n_part_out);

  int **pinit_location = NULL;
  PDM_g_num_t **_target_parent_ln_to_gn = NULL;
  PDM_part_assemble_partitions(extrp->comm,
                               distrib_partition,
                               distrib_elmt,
                               _elmt_part,
                               dequi_g_num,
                               dequi_init_location,
                               &extrp->pextract_n_entity[entity_type],
                               &_target_parent_ln_to_gn,
                               &pinit_location);
  free(distrib_partition);
  free(dequi_init_location);
  free(_elmt_part);
  free(dual_graph);
  free(dual_graph_idx);

  /*
   * Compute new global numbering for extract_part
   */
  extrp->pextract_entity_ln_to_gn[entity_type] = (PDM_g_num_t **) malloc(extrp->n_part_out * sizeof(PDM_g_num_t *));
  _compute_child(extrp->comm,
                 extrp->n_part_out,
                 extrp->pextract_n_entity[entity_type],
                 _target_parent_ln_to_gn,
                 extrp->pextract_entity_ln_to_gn[entity_type]);


  /*
   * Fake from target
   */
  for(int i_part = 0; i_part < extrp->n_part_out; ++i_part){
    extrp->n_target       [i_part] = extrp->pextract_n_entity       [entity_type][i_part];
    // extrp->target_gnum    [i_part] = extrp->pextract_entity_ln_to_gn[entity_type][i_part];
    extrp->target_gnum    [i_part] = _target_parent_ln_to_gn[i_part];
    extrp->target_location[i_part] = pinit_location[i_part];
  }
  _extract_part_and_reequilibrate_from_target(extrp);

  for(int i_part = 0; i_part < extrp->n_part_out; ++i_part){
    extrp->n_target       [i_part] = 0;
    extrp->target_gnum    [i_part] = NULL;
    extrp->target_location[i_part] = NULL;
    free(_target_parent_ln_to_gn[i_part]);
  }
  free(_target_parent_ln_to_gn);

  PDM_part_to_block_free(ptb_equi);

  // Il faut echanger les parent a travers le ptp OU a travers les blocks
  int           exch_request = -1;
  PDM_part_to_part_reverse_iexch(extrp->ptp_entity[entity_type],
                                 PDM_MPI_COMM_KIND_P2P,
                                 PDM_STRIDE_CST_INTERLACED,
                                 PDM_PART_TO_PART_DATA_DEF_ORDER_PART2,
                                 1,
                                 sizeof(PDM_g_num_t),
                                 NULL,
                (const void **)  entity_g_num,
                                 NULL,
                   (void ***)   &extrp->pextract_entity_parent_ln_to_gn[entity_type],
                                &exch_request);
  PDM_part_to_part_reverse_iexch_wait(extrp->ptp_entity[entity_type], exch_request);

  if(0 == 1) {
    for(int i_part = 0; i_part < extrp->n_part_out; ++i_part) {
      PDM_log_trace_array_long(extrp->pextract_entity_parent_ln_to_gn[entity_type][i_part], extrp->pextract_n_entity[entity_type][i_part], "pextract_entity_parent_ln_to_gn ::");
    }
  }

  for(int i_part = 0; i_part < extrp->n_part_out; ++i_part) {
    free(pinit_location[i_part]);
  }
  free(pinit_location);

}

/*=============================================================================
 * Public function definitions
 *============================================================================*/

/**
 *
 * \brief Build an extract_part struct
 *
 * \param [in]   dim                 Extraction dimension
 * \param [in]   n_part_in           Number of initial partition
 * \param [in]   n_part_out          Number of final partition
 * \param [in]   extract_kind        Extraction kind : (local/requilibrate/from target)
 * \param [in]   split_dual_method   Split method if requilibrate extract_kind
 * \param [in]   compute_child_gnum  Yes/No computation of a newest global numbering
 * \param [in]   ownership           Tell if you want ownership of resulting
 * \param [in]   comm                MPI communicator
 *
 * \return   Initialized \ref PDM_extract_part_t instance
 */
PDM_extract_part_t*
PDM_extract_part_create
(
 const int                     dim,
 const int                     n_part_in,
 const int                     n_part_out,
       PDM_extract_part_kind_t extract_kind,
       PDM_split_dual_t        split_dual_method,
       PDM_bool_t              compute_child_gnum,
       PDM_ownership_t         ownership,
       PDM_MPI_Comm            comm
)
{

  PDM_extract_part_t *extrp = (PDM_extract_part_t *) malloc(sizeof(PDM_extract_part_t));

  if(extract_kind == PDM_EXTRACT_PART_KIND_LOCAL) {
    if(n_part_in != n_part_out) {
      PDM_error(__FILE__, __LINE__, 0,"PDM_extract_part_create : cannot not equilibrate with not same number of n_part_in / n_part_out \n");
    }
  }

  extrp->extract_kind       = extract_kind;
  extrp->compute_child_gnum = compute_child_gnum;

  extrp->dim                   = dim;
  extrp->n_part_in             = n_part_in;
  extrp->n_part_out            = n_part_out;
  extrp->split_dual_method     = split_dual_method;
  extrp->ownership             = ownership;
  extrp->comm                  = comm;
  extrp->pmne                  = NULL;
  extrp->is_owner_extract_pmne = PDM_TRUE;
  extrp->extract_pmne          = NULL;

  extrp->n_cell         = (int          *) malloc(n_part_in * sizeof(int          ));
  extrp->n_face         = (int          *) malloc(n_part_in * sizeof(int          ));
  extrp->n_edge         = (int          *) malloc(n_part_in * sizeof(int          ));
  extrp->n_vtx          = (int          *) malloc(n_part_in * sizeof(int          ));
  extrp->n_extract      = (int          *) malloc(n_part_in * sizeof(int          ));

  extrp->pcell_face     = (int         **) malloc(n_part_in * sizeof(int         *));
  extrp->pcell_face_idx = (int         **) malloc(n_part_in * sizeof(int         *));
  extrp->pface_edge     = (int         **) malloc(n_part_in * sizeof(int         *));
  extrp->pface_edge_idx = (int         **) malloc(n_part_in * sizeof(int         *));
  extrp->pedge_vtx      = (int         **) malloc(n_part_in * sizeof(int         *));

  extrp->extract_lnum   = (int         **) malloc(n_part_in * sizeof(int         *));

  extrp->cell_ln_to_gn  = (PDM_g_num_t **) malloc(n_part_in * sizeof(PDM_g_num_t *));
  extrp->face_ln_to_gn  = (PDM_g_num_t **) malloc(n_part_in * sizeof(PDM_g_num_t *));
  extrp->edge_ln_to_gn  = (PDM_g_num_t **) malloc(n_part_in * sizeof(PDM_g_num_t *));
  extrp->vtx_ln_to_gn   = (PDM_g_num_t **) malloc(n_part_in * sizeof(PDM_g_num_t *));

  extrp->pface_vtx_idx  = (int         **) malloc(n_part_in * sizeof(int         *));
  extrp->pface_vtx      = (int         **) malloc(n_part_in * sizeof(int         *));

  extrp->pvtx_coord     = (double      **) malloc(n_part_in * sizeof(double      *));
  extrp->entity_center  = (double      **) malloc(n_part_in * sizeof(double      *));
  extrp->have_user_entity_center  = 0;

  for(int i_part = 0; i_part < n_part_in; ++i_part) {
    extrp->n_cell        [i_part] = 0;
    extrp->n_face        [i_part] = 0;
    extrp->n_edge        [i_part] = 0;
    extrp->n_vtx         [i_part] = 0;
    extrp->n_extract     [i_part] = 0;
    extrp->pcell_face    [i_part] = NULL;
    extrp->pcell_face_idx[i_part] = NULL;
    extrp->pface_edge    [i_part] = NULL;
    extrp->pface_edge_idx[i_part] = NULL;
    extrp->pedge_vtx     [i_part] = NULL;
    extrp->extract_lnum  [i_part] = NULL;
    extrp->cell_ln_to_gn [i_part] = NULL;
    extrp->face_ln_to_gn [i_part] = NULL;
    extrp->edge_ln_to_gn [i_part] = NULL;
    extrp->vtx_ln_to_gn  [i_part] = NULL;
    extrp->pface_vtx_idx [i_part] = NULL;
    extrp->pface_vtx     [i_part] = NULL;
    extrp->pvtx_coord    [i_part] = NULL;
    extrp->entity_center [i_part] = NULL;
  }

  for(int i = 0; i < PDM_BOUND_TYPE_MAX; ++i) {
    extrp->n_group[i] = 0;
  }

  for(int i = 0; i < PDM_BOUND_TYPE_MAX; ++i) {
    extrp->n_group_entity       [i] = NULL;
    extrp->group_entity         [i] = NULL;
    extrp->group_entity_ln_to_gn[i] = NULL;
  }

  extrp->from_target     = 0;
  extrp->n_target        = (int          *) malloc(n_part_out * sizeof(int          ));
  extrp->target_gnum     = (PDM_g_num_t **) malloc(n_part_out * sizeof(PDM_g_num_t *));
  extrp->target_ownership = PDM_OWNERSHIP_USER;
  extrp->target_location = (int         **) malloc(n_part_out * sizeof(int         *));

  for(int i_part = 0; i_part < n_part_out; ++i_part) {
    extrp->n_target       [i_part] = 0;
    extrp->target_gnum    [i_part] = NULL;
    extrp->target_location[i_part] = NULL;
  }

  if (dim == 3) {
    extrp->master_entity = PDM_MESH_ENTITY_CELL;
  }
  else if (dim == 2) {
    extrp->master_entity = PDM_MESH_ENTITY_FACE;
  }
  else if (dim == 1) {
    extrp->master_entity = PDM_MESH_ENTITY_EDGE;
  }
  else {
    extrp->master_entity = PDM_MESH_ENTITY_VTX;
  }

  extrp->is_owner_connectivity    = malloc( PDM_CONNECTIVITY_TYPE_MAX * sizeof(PDM_bool_t   ) );
  extrp->is_owner_ln_to_gn        = malloc( PDM_MESH_ENTITY_MAX       * sizeof(PDM_bool_t   ) );
  extrp->is_owner_parent_ln_to_gn = malloc( PDM_MESH_ENTITY_MAX       * sizeof(PDM_bool_t   ) );
  extrp->is_owner_parent_lnum     = malloc( PDM_MESH_ENTITY_MAX       * sizeof(PDM_bool_t   ) );
  extrp->is_owner_vtx_coord       = PDM_TRUE;

  for(int i = 0; i < PDM_CONNECTIVITY_TYPE_MAX; ++i) {
    extrp->is_owner_connectivity    [i] = PDM_TRUE;
    extrp->pextract_connectivity    [i] = NULL;
    extrp->pextract_connectivity_idx[i] = NULL;
  }
  for(int i = 0; i < PDM_MESH_ENTITY_MAX; ++i) {
    extrp->is_owner_ln_to_gn              [i] = PDM_TRUE;
    extrp->is_owner_parent_ln_to_gn       [i] = PDM_TRUE;
    extrp->is_owner_parent_lnum           [i] = PDM_TRUE;
    extrp->pextract_entity_ln_to_gn       [i] = NULL;
    extrp->pextract_entity_parent_ln_to_gn[i] = NULL;
    extrp->pextract_entity_parent_lnum    [i] = NULL;
    extrp->ptp_entity                     [i] = NULL;
    extrp->ptp_ownership                  [i] = PDM_OWNERSHIP_KEEP;
  }
  extrp->pextract_vtx_coord = NULL;

  for(int i = 0; i < PDM_MESH_ENTITY_MAX; ++i) {
    extrp->pextract_n_entity[i] = NULL;
  }

  for(int i = 0; i < PDM_BOUND_TYPE_MAX; ++i) {
    extrp->ptp_group_entity                     [i] = NULL;
    extrp->ptp_group_ownership                  [i] = NULL;
    extrp->group_array_ownership                [i] = NULL;
    extrp->pn_extract_group_entity              [i] = NULL;
    extrp->pextract_group_entity                [i] = NULL;
    extrp->pextract_group_entity_ln_to_gn       [i] = NULL;
    extrp->pextract_group_entity_parent_ln_to_gn[i] = NULL;
    extrp->is_owner_extract_group               [i] = NULL;
  }

  return extrp;
}


/**
 *
 * \brief Compute extraction
 *
 * \param [in]   extrp      PDM_extract_part_t
 *
 */
void
PDM_extract_part_compute
(
  PDM_extract_part_t        *extrp
)
{
  if(extrp->extract_kind == PDM_EXTRACT_PART_KIND_LOCAL) {
    if(extrp->pmne == NULL) {
      _extract_part(extrp);
    } else {
      _extract_part_nodal(extrp);
    }
  } else if(extrp->extract_kind == PDM_EXTRACT_PART_KIND_REEQUILIBRATE) {
    if(extrp->pmne == NULL) {
      _extract_part_and_reequilibrate(extrp);
    } else {
      PDM_error(__FILE__, __LINE__, 0,"PDM_extract_part_compute : PDM_EXTRACT_PART_KIND_REEQUILIBRATE for nodal mesh not yet implemented \n");
      _extract_part_and_reequilibrate_nodal(extrp);
    }
  } else if (extrp->extract_kind == PDM_EXTRACT_PART_KIND_FROM_TARGET) {
    if(extrp->pmne == NULL) {
      _extract_part_and_reequilibrate_from_target(extrp);
    } else {
      _extract_part_and_reequilibrate_nodal_from_target(extrp);
    }
  }
}


void
PDM_extract_part_part_set
(
  PDM_extract_part_t       *extrp,
  int                       i_part,
  int                       n_cell,
  int                       n_face,
  int                       n_edge,
  int                       n_vtx,
  int                      *cell_face_idx,
  int                      *cell_face,
  int                      *face_edge_idx,
  int                      *face_edge,
  int                      *edge_vtx,
  int                      *face_vtx_idx,
  int                      *face_vtx,
  PDM_g_num_t              *cell_ln_to_gn,
  PDM_g_num_t              *face_ln_to_gn,
  PDM_g_num_t              *edge_ln_to_gn,
  PDM_g_num_t              *vtx_ln_to_gn,
  double                   *vtx_coord
)
{
  extrp->n_cell        [i_part] = n_cell;
  extrp->n_face        [i_part] = n_face;
  extrp->n_edge        [i_part] = n_edge;
  extrp->n_vtx         [i_part] = n_vtx;
  extrp->pcell_face    [i_part] = cell_face;
  extrp->pcell_face_idx[i_part] = cell_face_idx;
  extrp->pface_edge    [i_part] = face_edge;
  extrp->pface_edge_idx[i_part] = face_edge_idx;
  extrp->pedge_vtx     [i_part] = edge_vtx;
  extrp->cell_ln_to_gn [i_part] = cell_ln_to_gn;
  extrp->face_ln_to_gn [i_part] = face_ln_to_gn;
  extrp->edge_ln_to_gn [i_part] = edge_ln_to_gn;
  extrp->vtx_ln_to_gn  [i_part] = vtx_ln_to_gn;
  extrp->pface_vtx_idx [i_part] = face_vtx_idx;
  extrp->pface_vtx     [i_part] = face_vtx;
  extrp->pvtx_coord    [i_part] = vtx_coord;
}

/**
 *
 * \brief Set partition group (optional)
 * \param [in]   extrp             PDM_extract_part_t
 * \param [in]   bound_type        Kind of group
 * \param [in]   n_group           Number of group of kind bound_type
 */
void
PDM_extract_part_n_group_set
(
  PDM_extract_part_t        *extrp,
  PDM_bound_type_t           bound_type,
  int                        n_group
)
{
  extrp->n_group[bound_type] = n_group;

  extrp->n_group_entity       [bound_type] = malloc(n_group * sizeof(int          * ));
  extrp->group_entity         [bound_type] = malloc(n_group * sizeof(int         ** ));
  extrp->group_entity_ln_to_gn[bound_type] = malloc(n_group * sizeof(PDM_g_num_t ** ));

  for(int i_group = 0; i_group < n_group; ++i_group) {
    extrp->n_group_entity       [bound_type][i_group] = malloc(extrp->n_part_in * sizeof(int           ));
    extrp->group_entity         [bound_type][i_group] = malloc(extrp->n_part_in * sizeof(int         * ));
    extrp->group_entity_ln_to_gn[bound_type][i_group] = malloc(extrp->n_part_in * sizeof(PDM_g_num_t * ));
    for(int i_part = 0; i_part < extrp->n_part_in; ++i_part) {
      extrp->n_group_entity       [bound_type][i_group][i_part] = 0;
      extrp->group_entity         [bound_type][i_group][i_part] = NULL;
      extrp->group_entity_ln_to_gn[bound_type][i_group][i_part] = NULL;
    }
  }
}

/**
 *
 * \brief Set partition group (optional)
 *
 * \param [in]   extrp             PDM_extract_part_t
 * \param [in]   i_part            part identifier
 * \param [in]   i_group           group identifier
 * \param [in]   bound_type        Kind of group
 * \param [in]   n_group_entity    Number of entity in current group
 * \param [in]   group_entity      List of entity in group (size = n_group_entity)
 *
 */
void
PDM_extract_part_part_group_set
(
  PDM_extract_part_t        *extrp,
  int                       i_part,
  int                       i_group,
  PDM_bound_type_t          bound_type,
  int                       n_group_entity,
  int                      *group_entity,
  PDM_g_num_t              *group_entity_ln_to_gn
)
{
  assert(extrp->n_group_entity       [bound_type][i_group][i_part] == 0   );
  assert(extrp->group_entity         [bound_type][i_group][i_part] == NULL);
  assert(extrp->group_entity_ln_to_gn[bound_type][i_group][i_part] == NULL);

  extrp->n_group_entity       [bound_type][i_group][i_part] = n_group_entity;
  extrp->group_entity         [bound_type][i_group][i_part] = group_entity;
  extrp->group_entity_ln_to_gn[bound_type][i_group][i_part] = group_entity_ln_to_gn;

}


/**
 *
 * \brief Set PDM_part_mesh_nodal_elmts_t
 *
 * \param [in]   extrp            PDM_extract_part_t structure
 * \param [in]   pmne             PDM_part_mesh_nodal_elmts_t corresponding of dimenstion
 *
 */
void
PDM_extract_part_part_nodal_set
(
  PDM_extract_part_t          *extrp,
  PDM_part_mesh_nodal_elmts_t *pmne
)
{
  assert(extrp->dim == pmne->mesh_dimension);

  extrp->pmne = pmne;
}

/**
 *
 * \brief Set the extract number
 *
 * \param [in]   extrp         PDM_extract_part_t
 * \param [in]   i_part        part identifier
 * \param [in]   n_extract     Number of entity to select
 * \param [in]   extract_lnum  List of id to extract (starting at 0)
 *
 */
void
PDM_extract_part_selected_lnum_set
(
  PDM_extract_part_t       *extrp,
  int                       i_part,
  int                       n_extract,
  int                      *extract_lnum
)
{
  extrp->n_extract   [i_part] = n_extract;
  extrp->extract_lnum[i_part] = extract_lnum;
}

/**
 *
 * \brief Set the extract target number
 *
 * \param [in]   extrp             PDM_extract_part_t
 * \param [in]   i_part            part identifier
 * \param [in]   n_target          Number of target to select
 * \param [in]   target_gnum       List of global id to extract
 * \param [in]   target_location   Init location (optional NULL pointer accepted and computed internaly)
 *
 */
void
PDM_extract_part_target_set
(
  PDM_extract_part_t       *extrp,
  int                       i_part,
  int                       n_target,
  PDM_g_num_t              *target_gnum,
  int                      *target_location
)
{
  extrp->from_target = 1;
  assert(extrp->extract_kind == PDM_EXTRACT_PART_KIND_FROM_TARGET);
  extrp->n_target       [i_part] = n_target;
  extrp->target_gnum    [i_part] = target_gnum;
  extrp->target_location[i_part] = target_location;
}


/**
 *
 * \brief Keep target_gnum data ownership inside extrp
 *
 * \param [in]   extrp             PDM_extract_part_t
 *
 */
void
PDM_extract_part_target_gnum_keep_ownnership
(
  PDM_extract_part_t       *extrp
)
{
  extrp->target_ownership = PDM_OWNERSHIP_KEEP;
}


/**
 *
 * \brief Set entity center (ussefull for equilibrate / hilbert ordering )
 *
 * \param [in]   extrp            PDM_extract_part_t structure
 * \param [in]   i_part           part identifier
 * \param [in]   entity_center    Center of entity (relative of dim throught create)
 *
 */
void
PDM_extract_part_entity_center_set
(
  PDM_extract_part_t          *extrp,
  int                          i_part,
  double                      *entity_center
)
{
  extrp->have_user_entity_center = 1;
  extrp->entity_center[i_part] = entity_center;
}




/**
 *
 * \brief Set PDM_part_mesh_nodal_elmts_t
 *
 * \param [in]   extrp               PDM_extract_part_t structure
 * \param [in]   i_part_out          part identifier
 * \param [in]   PDM_mesh_entities_t Kind of entity required \ref PDM_mesh_entities_t
 *
 * \return Number of entity
 */
int
PDM_extract_part_n_entity_get
(
 PDM_extract_part_t       *extrp,
 int                       i_part_out,
 PDM_mesh_entities_t       entity_type
)
{
  if(extrp->pextract_n_entity[entity_type] != NULL) {
    return extrp->pextract_n_entity[entity_type][i_part_out];
  } else {
    return 0;
  }
}


/**
 * \brief Return size of leading connectivity on current partition ( n_entity )
 * \param [in]  extrp               Pointer to \ref PDM_extract_part_t object
 * \param [in]  i_part              Id of part
 * \param [in]  connectivity_type   Connectivity kind \ref PDM_connectivity_type_t
 * \param [in]  connect             Connectivity array (size = connect_idx[n_entity] )
 * \param [in]  connect_idx         Connectivity index (size = n_entity+1 )
 * \param [in]  ownership           Choice of ownership of the resulting arrays \ref PDM_ownership_t
 */
int
PDM_extract_part_connectivity_get
(
 PDM_extract_part_t        *extrp,
 int                        i_part_out,
 PDM_connectivity_type_t    connectivity_type,
 int                      **connect,
 int                      **connect_idx,
 PDM_ownership_t           ownership
)
{
  assert(i_part_out < extrp->n_part_out);

  PDM_mesh_entities_t entity_type = PDM_connectivity_type_to_entity_type(connectivity_type);

  if(extrp->pextract_connectivity[connectivity_type] != NULL) {
    *connect     = extrp->pextract_connectivity    [connectivity_type][i_part_out];
  } else {
    *connect     = NULL;
  }

  if(extrp->pextract_connectivity_idx[connectivity_type] != NULL) {
    *connect_idx = extrp->pextract_connectivity_idx[connectivity_type][i_part_out];
  } else {
    *connect_idx = NULL; // edge_vtx / face_cell
  }

  if(ownership == PDM_OWNERSHIP_USER || ownership == PDM_OWNERSHIP_UNGET_RESULT_IS_FREE) {
    extrp->is_owner_connectivity[connectivity_type] = PDM_FALSE;
  } else {
    extrp->is_owner_connectivity[connectivity_type] = PDM_TRUE;
  }

  if(extrp->pextract_n_entity[entity_type] != NULL) {
    return extrp->pextract_n_entity[entity_type][i_part_out];
  } else {
    return 0;
  }
}



/**
 *
 * \brief Return size of entity_type on current partition ( n_entity )
 * \param [in]  extrp             Pointer to \ref PDM_extract_part_t object
 * \param [in]  i_part            Id of part
 * \param [in]  entity_type       Entity kind \ref PDM_mesh_entities_t)
 * \param [out] entity_ln_to_gn   Entity local numbering to global numbering (size = n_entity, numbering : 1 to n)
 * \param [in]  ownership         Ownership for entity_ln_to_gn ( \ref PDM_ownership_t )
 */
int
PDM_extract_part_ln_to_gn_get
(
 PDM_extract_part_t        *extrp,
 int                        i_part_out,
 PDM_mesh_entities_t        entity_type,
 PDM_g_num_t              **pentity_ln_to_gn,
 PDM_ownership_t            ownership
)
{
  if(extrp->pextract_n_entity[entity_type] != NULL) {
    *pentity_ln_to_gn = NULL;
    if(extrp->pextract_entity_ln_to_gn[entity_type] != NULL) {
      *pentity_ln_to_gn = extrp->pextract_entity_ln_to_gn[entity_type][i_part_out];
    }
    if(ownership == PDM_OWNERSHIP_USER || ownership == PDM_OWNERSHIP_UNGET_RESULT_IS_FREE) {
      extrp->is_owner_ln_to_gn[entity_type] = PDM_FALSE;
    } else {
      extrp->is_owner_ln_to_gn[entity_type] = PDM_TRUE;
    }

    return extrp->pextract_n_entity[entity_type][i_part_out];
  }
  else {
    return 0;
  }
}


/**
 *
 * \brief Return size of entity_type on current partition ( n_entity )
 * \param [in]  extrp                  Pointer to \ref PDM_extract_part_t object
 * \param [in]  i_part                 Id of part
 * \param [in]  entity_type            Entity kind \ref PDM_mesh_entities_t)
 * \param [out] parent_entity_ln_to_gn Entity local numbering to global numbering of parent entity, correspond to the input mesh (size = n_entity, numbering : 1 to n)
 * \param [in]  ownership              Ownership for entity_ln_to_gn ( \ref PDM_ownership_t )
 */
int
PDM_extract_part_parent_ln_to_gn_get
(
 PDM_extract_part_t        *extrp,
 int                        i_part_out,
 PDM_mesh_entities_t        entity_type,
 PDM_g_num_t              **parent_entity_ln_to_gn,
 PDM_ownership_t            ownership
)
{

  int n_entity = 0;

  if (entity_type == extrp->master_entity && extrp->target_ownership != PDM_OWNERSHIP_KEEP && extrp->extract_kind == PDM_EXTRACT_PART_KIND_FROM_TARGET) {
    PDM_error(__FILE__, __LINE__, 0, "Error PDM_extract_part_parent_ln_to_gn_get : parent_ln_to_gn is not available for the master entity,"
                                     " call PDM_extract_part_target_gnum_keep_ownnership to get it\n");
  }

  if(extrp->pextract_n_entity[entity_type] != NULL) {

    if(extrp->pextract_entity_parent_ln_to_gn[entity_type] != NULL) {
      *parent_entity_ln_to_gn = extrp->pextract_entity_parent_ln_to_gn[entity_type][i_part_out];
    }

    if(ownership == PDM_OWNERSHIP_USER || ownership == PDM_OWNERSHIP_UNGET_RESULT_IS_FREE) {
      extrp->is_owner_parent_ln_to_gn[entity_type] = PDM_FALSE;
    }
    else {
      extrp->is_owner_parent_ln_to_gn[entity_type] = PDM_TRUE;
    }

    n_entity = extrp->pextract_n_entity[entity_type][i_part_out];
  }
  else {
   *parent_entity_ln_to_gn = NULL;
  }


  return n_entity;
}

/**
 *
 * \brief Return size of entity_type on current partition ( n_entity )
 * \param [in]  extrp                  Pointer to \ref PDM_extract_part_t object
 * \param [in]  i_part                 Id of part
 * \param [in]  entity_type            Entity kind \ref PDM_mesh_entities_t)
 * \param [out] parent_entity_lnum     Local indexes of parent entity, correspond to the input mesh (size = n_entity)
 * \param [in]  ownership              Ownership for entity_ln_to_gn ( \ref PDM_ownership_t )
 */
int
PDM_extract_part_parent_lnum_get
(
 PDM_extract_part_t        *extrp,
 int                        i_part_out,
 PDM_mesh_entities_t        entity_type,
 int                      **parent_entity_lnum,
 PDM_ownership_t            ownership
)
{
  *parent_entity_lnum = extrp->pextract_entity_parent_lnum[entity_type][i_part_out];
  if(ownership == PDM_OWNERSHIP_USER || ownership == PDM_OWNERSHIP_UNGET_RESULT_IS_FREE) {
    extrp->is_owner_parent_lnum[entity_type] = PDM_FALSE;
  } else {
    extrp->is_owner_parent_lnum[entity_type] = PDM_TRUE;
  }

  return extrp->pextract_n_entity[entity_type][i_part_out];
}


/**
 *
 * \brief Get the vertex coordinates on current i_part partition and return number of vertices
 *
 * \param [in]   extrp      Pointer to \ref PDM_extract_part_t object
 * \param [in]   i_part     Id of part
 * \param [out]  vtx_coord  Vertex coordinate (size = 3 * n_vtx)
 * \param [in]   ownership  Ownership for color ( \ref PDM_ownership_t )
 */
int
PDM_extract_part_vtx_coord_get
(
 PDM_extract_part_t         *extrp,
 int                        i_part_out,
 double                   **pvtx_coord,
 PDM_ownership_t            ownership
)
{
  if(extrp->pextract_vtx_coord != NULL){
    *pvtx_coord = extrp->pextract_vtx_coord[i_part_out];
    if(ownership == PDM_OWNERSHIP_USER || ownership == PDM_OWNERSHIP_UNGET_RESULT_IS_FREE) {
      extrp->is_owner_vtx_coord = PDM_FALSE;
    } else {
      extrp->is_owner_vtx_coord = PDM_TRUE;
    }
    return extrp->pextract_n_entity[PDM_MESH_ENTITY_VTX][i_part_out];
  } else {
    *pvtx_coord = NULL;
    return 0;
  }
}


void
PDM_extract_part_part_mesh_nodal_get
(
  PDM_extract_part_t           *extrp,
  PDM_part_mesh_nodal_elmts_t **extract_pmne,
  PDM_ownership_t               ownership
)
{
  *extract_pmne = extrp->extract_pmne;

  if(ownership == PDM_OWNERSHIP_USER || ownership == PDM_OWNERSHIP_UNGET_RESULT_IS_FREE) {
    extrp->is_owner_extract_pmne = PDM_FALSE;
  } else {
    extrp->is_owner_extract_pmne = PDM_TRUE;
  }
}

void
PDM_extract_part_free
(
  PDM_extract_part_t  *extrp
)
{
  free(extrp->n_cell        );
  free(extrp->n_face        );
  free(extrp->n_edge        );
  free(extrp->n_vtx         );
  free(extrp->n_extract     );
  free(extrp->n_target      );

  free(extrp->pcell_face    );
  free(extrp->pcell_face_idx);

  free(extrp->pface_edge    );
  free(extrp->pface_edge_idx);

  free(extrp->pface_vtx     );
  free(extrp->pface_vtx_idx );
  free(extrp->pedge_vtx     );
  free(extrp->entity_center );

  if(extrp->from_target == 1) {
    for(int i_part = 0; i_part < extrp->n_part_in; ++i_part) {
      free(extrp->extract_lnum[i_part]);
    }
  }

  for(int i = 0; i < PDM_BOUND_TYPE_MAX; ++i) {
    for(int i_group = 0; i_group < extrp->n_group[i]; ++i_group) {
      if(extrp->n_group_entity[i][i_group] != NULL) {
        free(extrp->n_group_entity[i][i_group]);
      }
      if(extrp->group_entity[i][i_group] != NULL) {
        free(extrp->group_entity[i][i_group]);
      }
      if(extrp->group_entity_ln_to_gn[i][i_group]) {
        free(extrp->group_entity_ln_to_gn[i][i_group]);
      }
    }
    if(extrp->n_group_entity[i] != NULL) {
      free(extrp->n_group_entity   [i]);
    }
    if(extrp->group_entity[i] != NULL) {
      free(extrp->group_entity[i]);
    }
    if(extrp->group_entity_ln_to_gn[i]) {
      free(extrp->group_entity_ln_to_gn[i]);
    }
  }

  free(extrp->extract_lnum   );

  free(extrp->target_gnum    );
  free(extrp->target_location);

  free(extrp->cell_ln_to_gn );
  free(extrp->face_ln_to_gn );
  free(extrp->edge_ln_to_gn );
  free(extrp->vtx_ln_to_gn  );

  free(extrp->pvtx_coord);

  /*
   * Free extracted partition
   */
  PDM_extract_part_partial_free(extrp);

  free(extrp->is_owner_connectivity   );
  free(extrp->is_owner_ln_to_gn       );
  free(extrp->is_owner_parent_ln_to_gn);
  free(extrp->is_owner_parent_lnum);

  free(extrp);
}

/**
 *
 * \brief Free all resulting array if not owner
 *
 * \param [in]   extrp      Pointer to \ref PDM_extract_part_t object
 */
void
PDM_extract_part_partial_free
(
  PDM_extract_part_t  *extrp
)
{
  /*
   * Free extracted partition if owner
   */
  /* Free connectivity */
  for(int i = 0; i < PDM_CONNECTIVITY_TYPE_MAX; ++i) {
    if(extrp->is_owner_connectivity[i] == PDM_TRUE) {
      if(extrp->pextract_connectivity[i] != NULL) {
        for(int i_part = 0; i_part < extrp->n_part_out; ++i_part) {
          if(extrp->pextract_connectivity[i][i_part] != NULL) {
            free(extrp->pextract_connectivity[i][i_part]);
          }
          if(extrp->pextract_connectivity_idx[i] != NULL) {
            if(extrp->pextract_connectivity_idx[i][i_part] != NULL) {
              free(extrp->pextract_connectivity_idx[i][i_part]);
            }
          }
        }
      }
    }

    if(extrp->pextract_connectivity[i] != NULL) {
      free(extrp->pextract_connectivity[i]);
      extrp->pextract_connectivity[i] = NULL;
    }

    if(extrp->pextract_connectivity_idx[i] != NULL) {
      free(extrp->pextract_connectivity_idx[i]);
      extrp->pextract_connectivity_idx[i] = NULL;
    }
  }

  /* Free ln_to_gn */
  for(int i = 0; i < PDM_MESH_ENTITY_MAX; ++i) {
    if(extrp->pextract_entity_ln_to_gn[i] != NULL) {
      if(extrp->is_owner_ln_to_gn[i] == PDM_TRUE) {
        for(int i_part = 0; i_part < extrp->n_part_out; ++i_part) {
          if(extrp->pextract_entity_ln_to_gn[i][i_part] != NULL) {
            free(extrp->pextract_entity_ln_to_gn[i][i_part]);
          }
        }
      }
      free(extrp->pextract_entity_ln_to_gn[i]);
      extrp->pextract_entity_ln_to_gn[i] = NULL;
    }
  }


  /* Free parent_ln_to_gn */
  for(int i = 0; i < PDM_MESH_ENTITY_MAX; ++i) {

    if ((extrp->from_target == 1) &&
        (extrp->target_ownership == PDM_OWNERSHIP_KEEP) &&
        ((int) extrp->master_entity == i)) {
      if(extrp->pextract_entity_parent_ln_to_gn[i] != NULL) {
        if(extrp->is_owner_parent_ln_to_gn[i] == PDM_TRUE) {
          for(int i_part = 0; i_part < extrp->n_part_out; ++i_part) {
            if(extrp->pextract_entity_parent_ln_to_gn[i][i_part] != NULL) {
              free(extrp->pextract_entity_parent_ln_to_gn[i][i_part]);
            }
          }
        }
        free(extrp->pextract_entity_parent_ln_to_gn[i]);
        extrp->pextract_entity_parent_ln_to_gn[i] = NULL;
      }
    }
    else {
      if(extrp->pextract_entity_parent_ln_to_gn[i] != NULL) {
        if(extrp->is_owner_parent_ln_to_gn[i] == PDM_TRUE) {
          for(int i_part = 0; i_part < extrp->n_part_out; ++i_part) {
            if(extrp->pextract_entity_parent_ln_to_gn[i][i_part] != NULL) {
              free(extrp->pextract_entity_parent_ln_to_gn[i][i_part]);
            }
          }
        }
        free(extrp->pextract_entity_parent_ln_to_gn[i]);
        extrp->pextract_entity_parent_ln_to_gn[i] = NULL;
      }
    }
  }

  /* Free parent_ln_to_gn */
  for(int i = 0; i < PDM_MESH_ENTITY_MAX; ++i) {
    if(extrp->pextract_entity_parent_lnum[i] != NULL) {
      if(extrp->is_owner_parent_lnum[i] == PDM_TRUE) {
        for(int i_part = 0; i_part < extrp->n_part_out; ++i_part) {
          if(extrp->pextract_entity_parent_lnum[i][i_part] != NULL) {
            free(extrp->pextract_entity_parent_lnum[i][i_part]);
          }
        }
      }
      free(extrp->pextract_entity_parent_lnum[i]);
      extrp->pextract_entity_parent_lnum[i] = NULL;
    }
  }

  if(extrp->is_owner_vtx_coord == PDM_TRUE) {
    if(extrp->pextract_vtx_coord != NULL){
      for(int i_part = 0; i_part < extrp->n_part_out; ++i_part) {
        if(extrp->pextract_vtx_coord[i_part] != NULL){
          free(extrp->pextract_vtx_coord[i_part]);
        }
      }
    }
  }
  if(extrp->pextract_vtx_coord != NULL){
    free(extrp->pextract_vtx_coord);
  }

  for(int i = 0; i < PDM_MESH_ENTITY_MAX; ++i) {
    free(extrp->pextract_n_entity[i]);
  }

  /*
   * Free bound
   */
  for(int i = 0; i < PDM_BOUND_TYPE_MAX; ++i) {
    if(extrp->ptp_group_entity[i] != NULL) {
      for(int i_group = 0; i_group < extrp->n_group[i]; ++i_group) {
        if(extrp->ptp_group_entity   [i][i_group] != NULL &&
           extrp->ptp_group_ownership[i][i_group] == PDM_OWNERSHIP_KEEP) {
          PDM_part_to_part_free(extrp->ptp_group_entity   [i][i_group]);
        }

        /* Free array */
        // if(extrp->is_owner_extract_group[i][i_group] == PDM_TRUE) {
        if(extrp->group_array_ownership[i][i_group] == PDM_OWNERSHIP_KEEP) {
          for(int i_part = 0; i_part < extrp->n_part_out; ++i_part) {
            free(extrp->pextract_group_entity                [i][i_group][i_part]);
            free(extrp->pextract_group_entity_ln_to_gn       [i][i_group][i_part]);
            free(extrp->pextract_group_entity_parent_ln_to_gn[i][i_group][i_part]);
          }
        }
        free(extrp->pn_extract_group_entity              [i][i_group]);
        free(extrp->pextract_group_entity                [i][i_group]);
        free(extrp->pextract_group_entity_ln_to_gn       [i][i_group]);
        free(extrp->pextract_group_entity_parent_ln_to_gn[i][i_group]);
      }
      free(extrp->pn_extract_group_entity              [i]);
      free(extrp->pextract_group_entity                [i]);
      free(extrp->pextract_group_entity_ln_to_gn       [i]);
      free(extrp->pextract_group_entity_parent_ln_to_gn[i]);

      free(extrp->is_owner_extract_group  [i]);
      free(extrp->ptp_group_entity        [i]);
      free(extrp->ptp_group_ownership     [i]);
      free(extrp->group_array_ownership   [i]);

    }
  }

  if(extrp->is_owner_extract_pmne == PDM_TRUE && extrp->extract_pmne != NULL) {
    PDM_part_mesh_nodal_elmts_free(extrp->extract_pmne);
  }

  for (int i = 0; i < PDM_MESH_ENTITY_MAX; ++i) {
    if (extrp->ptp_ownership[i] == PDM_OWNERSHIP_KEEP) {
      PDM_part_to_part_free(extrp->ptp_entity[i]);
    }
  }

}

/**
 *
 * \brief Get for entity_type (cells/faces/edge/vertices) the associated part_to_part (\ref PDM_part_to_part_t ). The part to part exchange protocol allow user to
 * exchange easily data from input mesh to the extract one.
 *
 * \param [in]   extrp        Pointer to \ref PDM_extract_part_t object
 * \param [in]   entity_type  Bound type \ref PDM_mesh_entities_t
 * \param [out]  ptp          Part to part protocol exchange, to exchange betwenn the input mesh and the output one (\ref PDM_part_to_part_t)
 * \param [in]   ownership    Ownership for color ( \ref PDM_ownership_t )
 *
 */
void
PDM_extract_part_part_to_part_get
(
       PDM_extract_part_t   *extrp,
 const PDM_mesh_entities_t   entity_type,
       PDM_part_to_part_t  **ptp,
       PDM_ownership_t       ownership

 )
{
  *ptp = extrp->ptp_entity[entity_type];

  extrp->ptp_ownership[entity_type] = ownership;
}


/**
 *
 * \brief Get for bound_type the associated part_to_part (\ref PDM_part_to_part_t ). The part to part exchange protocol allow user to
 * exchange easily data from input mesh to the extract one.
 *
 * \param [in]   extrp        Pointer to \ref PDM_extract_part_t object
 * \param [in]   bound_type   Bound type \ref PDM_bound_type_t
 * \param [in]   i_group      Id of group
 * \param [out]  ptp          Part to part protocol exchange, to exchange betwenn the input mesh and the output one (\ref PDM_part_to_part_t)
 * \param [in]   ownership    Ownership for color ( \ref PDM_ownership_t )
 *
 */
void
PDM_extract_part_part_to_part_group_get
(
       PDM_extract_part_t   *extrp,
 const PDM_bound_type_t      bound_type,
       int                   i_group,
       PDM_part_to_part_t  **ptp,
       PDM_ownership_t       ownership

)
{
  *ptp = extrp->ptp_group_entity[bound_type][i_group];

  extrp->ptp_group_ownership[bound_type][i_group] = ownership;
}


/**
 *
 * \brief Get the bound description for the entity (cell/face/edge/vertices)
 *
 * \param [in]   extrp                              Pointer to \ref PDM_extract_part_t object
 * \param [in]   bound_type                             Bound type \ref PDM_bound_type_t
 * \param [in]   i_part                                 Id of part
 * \param [in]   i_group                                Id of group
 * \param [out]  pn_extract_group_entity                Number of entity in current group
 * \param [out]  pextract_group_entity_ln_to_gn         Entity global id in current partition (size = pn_extract_group_entity)
 * \param [out]  pextract_group_entity_parent_ln_to_gn  Entity global id in parent partition (size = pn_extract_group_entity)
 * \param [in]   ownership      Ownership for color ( \ref PDM_ownership_t )
 *
 */
void
PDM_extract_part_group_get
(
       PDM_extract_part_t   *extrp,
 const PDM_bound_type_t      bound_type,
       int                   i_part,
       int                   i_group,
       int                  *pn_extract_group_entity,
       int                 **pextract_group_entity,
       PDM_g_num_t         **pextract_group_entity_ln_to_gn,
       PDM_g_num_t         **pextract_group_entity_parent_ln_to_gn,
       PDM_ownership_t       ownership
)
{

  *pn_extract_group_entity               = extrp->pn_extract_group_entity              [bound_type][i_group][i_part];
  *pextract_group_entity                 = extrp->pextract_group_entity                [bound_type][i_group][i_part];
  *pextract_group_entity_ln_to_gn        = extrp->pextract_group_entity_ln_to_gn       [bound_type][i_group][i_part];
  *pextract_group_entity_parent_ln_to_gn = extrp->pextract_group_entity_parent_ln_to_gn[bound_type][i_group][i_part];
  extrp->group_array_ownership   [bound_type][i_group] = ownership;

}



#ifdef __cplusplus
}
#endif /* __cplusplus */
